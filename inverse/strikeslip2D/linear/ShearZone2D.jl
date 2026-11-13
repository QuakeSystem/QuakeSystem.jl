# make sure we use the local environment
using Pkg; Pkg.activate(@__DIR__)

const isCUDA = false

@static if isCUDA
    using CUDA
end

using JustRelax, JustRelax.JustRelax2D

const backend = @static if isCUDA
    CUDABackend # Options: CPUBackend, CUDABackend, AMDGPUBackend
else
    JustRelax.CPUBackend
end

using ParallelStencil, ParallelStencil.FiniteDifferences2D

@static if isCUDA
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end

using JustPIC
const backend_JP = @static if isCUDA
    CUDA.CUDABackend # Options: JustPIC.CPU, CUDA.CUDABackend, AMDGPU.ROCBackend
else
    JustPIC.CPU
end

using GeoParams, CairoMakie, LaTeXStrings, Printf

const SEC_PER_YR = 3.15576e7
const MS_TO_CMYR = 100 * SEC_PER_YR

# Console output: 0 silent, 1 one line per step, 2 diagnostics and solver iteration log
const VERBOSE = Ref(1)
log_at(level) = VERBOSE[] >= level

# LinearTicks: the default WilkinsonTicks prints "No strict ticks found" for some ranges
set_theme!(
    theme_latexfonts();
    Axis = (xticks = LinearTicks(5), yticks = LinearTicks(5)),
    Colorbar = (ticks = LinearTicks(5),),
)

include(joinpath(@__DIR__, "ShearZone2D_rheology.jl"))
include(joinpath(@__DIR__, "ShearZone2D_setup.jl"))
include(joinpath(@__DIR__, "localModules", "geometry.jl"))
include(joinpath(@__DIR__, "localModules", "sensitivities.jl"))
include(joinpath(@__DIR__, "localModules", "plotting.jl"))
include(joinpath(@__DIR__, "localModules", "sensitivityAnalysis.jl"))

is_step(it, every, nt) = it == nt || (!isnothing(every) && every > 0 && it % every == 0)

## BEGIN OF MAIN SCRIPT --------------------------------------------------------------
function main(setup, igg)
    (; domain, zone, material, rheology, loading, time, adjoint, solver, output) = setup
    (; dt, nt) = time
    VERBOSE[] = output.verbose
    FIG_SETTINGS.show[] = output.show
    dir_fwd = joinpath(output.figdir, "forward")
    dir_adj = joinpath(output.figdir, "adjoint")
    zone_kw = (x_min = zone.x[1], x_max = zone.x[2], y_min = zone.y[1], y_max = zone.y[2])

    # Grid -----------------------------------------------
    ni = domain.nx, domain.ny
    li = domain.lx, domain.ly
    grid = Geometry(ni, li; origin = domain.origin)
    (; xvi) = grid
    # ----------------------------------------------------

    # Phases -------------------------------------------
    phase_ratios = sharp_phase_ratios(grid, zone, length(rheology))
    # ----------------------------------------------------

    # Stokes ---------------------------------------------
    # Periodic in x
    flow_bcs = VelocityBoundaryConditions(;
        no_slip = (left = false, right = false, top = false, bot = false),
        free_slip = (left = false, right = false, top = false, bot = false),
        periodic = (left = true, right = true, top = false, bot = false),
    )
    stokes = StokesArrays(backend, ni, flow_bcs)
    stokes_ad = AdjointStokesArrays(backend, ni, flow_bcs)
    gradients = material_controls(backend, ni, adjoint.target_parameters; nphases = length(rheology))
    observation = box_observation(
        stokes_ad.V.Vx, grid, adjoint.center, adjoint.half_width;
        antisymmetric = adjoint.antisymmetric,
    )
    # ----------------------------------------------------

    # Body force (uniform density), pressure, viscosity ---
    ρg = (@zeros(ni...), @zeros(ni...))
    T = @ones(ni .+ 2...) .* material.T
    stokes.P .= material.P_bg
    stokes.P0 .= material.P_bg
    compute_viscosity!(stokes, phase_ratios, (T = T, P = stokes.P, dt = dt), rheology, material.viscosity_cutoff)
    # ----------------------------------------------------

    # Boundary conditions ---------------------------------
    # Simple shear Vx = 2εbg·y on every Vx row, the two ghost rows half a cell outside the walls
    # included: the wall velocity is then ±εbg·ly. The top and bottom carry no condition, so
    # the solver never overwrites the ghost rows.
    yVx = grid.xi_vel[1][2]
    stokes.V.Vx .= PTArray(backend)([2 * loading.εbg * y for _ in xvi[1], y in yVx])
    fill!(stokes.V.Vy, 0.0)
    flow_bcs!(stokes, flow_bcs)
    update_halo!(@velocity(stokes)...)
    # ----------------------------------------------------

    # DYREL ----------------------------------------------
    # The penalty γ_num = γfact·mean(η) uses the creep viscosity, ≫ η_ve in the matrix:
    # by default γfact gives γ_num = 20·mean(η_ve). It only affects convergence.
    η = filter(isfinite, Array(stokes.viscosity.η))
    η_mean = sum(η) / length(η)
    η_ve_mean = sum(@. 1 / (1 / η + 1 / (material.G * dt))) / length(η)
    γfact = something(solver.γfact, 20 * η_ve_mean / η_mean)
    dyrel = DYREL(backend, stokes, rheology, phase_ratios, grid.di, dt; ϵ = solver.ϵ, γfact = γfact)
    # ----------------------------------------------------

    if log_at(1)
        print_setup(setup)
        @printf("DYREL penalty: γfact = %.3e (γ_num = %.3e Pa s)\n", γfact, γfact * η_mean)
        @printf("Objective: %d Vx nodes in the box, adjoint = %s\n", count(!iszero, Array(observation.weights)), adjoint.run)
    end
    output.plot_every === nothing || mkpath.((output.figdir, dir_fwd, dir_adj))

    # Time loop ------------------------------------------
    history = NamedTuple[]
    args = (; T = T, P = stokes.P, dt = dt, ΔT = @zeros(ni .+ 2...))
    for it in 1:nt
        t = it * dt
        run_adjoint = adjoint.run && is_step(it, adjoint.every, nt)

        t_solve = @elapsed out = solve_DYREL!(
            stokes, stokes_ad, ρg, dyrel, flow_bcs, phase_ratios, rheology, args, grid, dt, igg;
            solver.kwargs...,
            adjoint = run_adjoint,
            observation = run_adjoint ? observation : nothing,
            gradients = gradients,
        )

        J = objective_value(stokes, observation)
        if run_adjoint
            totals = sensitivity_totals(gradients, stokes_ad, phase_ratios, rheology, adjoint.target_parameters)
            push!(history, (; step = it, t, J, totals))
            log_at(2) && print_sensitivity_totals(totals, J; phase_names = material.phase_names)
            plot_results = output.plot_every !== nothing
            types = output.sensitivity_types
            run_sensitivity_analysis(
                stokes, stokes_ad, gradients, grid, observation, dir_adj;
                ρ_matrix = material.ρ, ρ_weak = material.ρ, G = material.G, K = material.Kb,
                zone_kw..., step = it, t = t, τM = time.τM_weak,
                plot = plot_results ? types : (),
                combined = plot_results && :combined in types,
            )
            if plot_results && :absolute in types
                for mode in output.sensitivity_modes
                    plot_sensitivities(
                        stokes_ad, gradients, grid, dir_adj;
                        mode = mode, zone_kw...,
                        observation_center = adjoint.center,
                        observation_half_width = adjoint.half_width, step = it,
                    )
                end
            end
        end
        if output.plot_every !== nothing && is_step(it, output.plot_every, nt)
            plot_forward_fields(stokes, grid, it, t, dir_fwd; zone_kw...)
        end

        if log_at(1)
            err = isempty(out.err_evo_tot) ? NaN : last(out.err_evo_tot)
            @printf(
                "Step %3d/%d%s | %7d iterations | err = %.2e | J = %.4f cm/yr | %.1f s\n",
                it, nt, run_adjoint ? " + adjoint" : "          ", out.iter, err,
                J * MS_TO_CMYR, t_solve
            )
        end
    end
    # ----------------------------------------------------

    if !isempty(history)
        h = last(history)
        println("\nSensitivities at step $(h.step):")
        print_sensitivity_totals(h.totals, h.J; phase_names = material.phase_names)
        dG, dKb = (sum(r.dJdp for r in h.totals if r.name === n) for n in (:G, :Kb))
        nus = nu_sensitivities(dG, dKb, material.G, material.ν)
        @printf("  elasticity from (G, ν): dJ/dν = % .6e, dJ/dG at fixed ν = % .6e\n", nus.ν, nus.G_ν)
    end

    return (; stokes, stokes_ad, gradients, grid, phase_ratios, observation, history)
end
## END OF MAIN SCRIPT ----------------------------------------------------------------

# Run ---------------------------------------------------------------------------------
setup = setup_linear()

if ImplicitGlobalGrid.grid_is_initialized()
    finalize_global_grid(; finalize_MPI = false)
end
igg = IGG(init_global_grid(setup.domain.nx, setup.domain.ny, 1; init_MPI = !JustRelax.MPI.Initialized())...)

result = main(setup, igg)
