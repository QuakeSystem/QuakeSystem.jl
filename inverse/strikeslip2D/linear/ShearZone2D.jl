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

## OUTPUT SETTINGS --------------------------------------------------------------------

"""
    configure_output!(output)

Set the console verbosity and whether saved figures are also displayed.
"""
function configure_output!(output)
    VERBOSE[] = output.verbose
    FIG_SETTINGS.show[] = output.show
    return nothing
end

figure_dirs(output) = (; fwd = joinpath(output.figdir, "forward"), adj = joinpath(output.figdir, "adjoint"))

zone_keywords(zone) = (x_min = zone.x[1], x_max = zone.x[2], y_min = zone.y[1], y_max = zone.y[2])

"""
    prepare_figure_dirs(output)

Create the figure folders, unless figures are off (`plot_every = nothing`).
"""
function prepare_figure_dirs(output)
    dirs = figure_dirs(output)
    output.plot_every === nothing || mkpath.((output.figdir, dirs.fwd, dirs.adj))
    return nothing
end

## MODEL CONSTRUCTION -----------------------------------------------------------------

"""
    allocate_stokes(grid, ni, nphases, adjoint)

Forward and adjoint Stokes arrays (periodic in x), gradient buffers and the observation box.
"""
function allocate_stokes(grid, ni, nphases, adjoint)
    flow_bcs = VelocityBoundaryConditions(;
        no_slip = (left = false, right = false, top = false, bot = false),
        free_slip = (left = false, right = false, top = false, bot = false),
        periodic = (left = true, right = true, top = false, bot = false),
    )
    stokes = StokesArrays(backend, ni, flow_bcs)
    stokes_ad = AdjointStokesArrays(backend, ni, flow_bcs)
    gradients = material_controls(backend, ni, adjoint.target_parameters; nphases = nphases)
    observation = box_observation(
        stokes_ad.V.Vx, grid, adjoint.center, adjoint.half_width;
        antisymmetric = adjoint.antisymmetric,
    )
    return (; flow_bcs, stokes, stokes_ad, gradients, observation)
end

"""
    initialise_material_state!(stokes, phase_ratios, rheology, material, ni, dt)

Zero body force (uniform density), background pressure, temperature and viscosity.
"""
function initialise_material_state!(stokes, phase_ratios, rheology, material, ni, dt)
    ρg = (@zeros(ni...), @zeros(ni...))
    T = @ones(ni .+ 2...) .* material.T
    stokes.P .= material.P_bg
    stokes.P0 .= material.P_bg
    compute_viscosity!(stokes, phase_ratios, (T = T, P = stokes.P, dt = dt), rheology, material.viscosity_cutoff)
    return (; ρg, T)
end

"""
    apply_simple_shear!(stokes, flow_bcs, grid, εbg)

Simple shear Vx = 2εbg·y on every Vx row, the two ghost rows half a cell outside the walls
included: the wall velocity is then ±εbg·ly. The top and bottom carry no condition, so
the solver never overwrites the ghost rows.
"""
function apply_simple_shear!(stokes, flow_bcs, grid, εbg)
    yVx = grid.xi_vel[1][2]
    stokes.V.Vx .= PTArray(backend)([2 * εbg * y for _ in grid.xvi[1], y in yVx])
    fill!(stokes.V.Vy, 0.0)
    flow_bcs!(stokes, flow_bcs)
    update_halo!(@velocity(stokes)...)
    return nothing
end

"""
    build_dyrel(stokes, rheology, phase_ratios, grid, material, solver, dt)

DYREL solver state. The penalty γ_num = γfact·mean(η) uses the creep viscosity, ≫ η_ve in
the matrix: by default γfact gives γ_num = 20·mean(η_ve). It only affects convergence.
"""
function build_dyrel(stokes, rheology, phase_ratios, grid, material, solver, dt)
    η = filter(isfinite, Array(stokes.viscosity.η))
    η_mean = sum(η) / length(η)
    η_ve_mean = sum(@. 1 / (1 / η + 1 / (material.G * dt))) / length(η)
    γfact = something(solver.γfact, 20 * η_ve_mean / η_mean)
    dyrel = DYREL(backend, stokes, rheology, phase_ratios, grid.di, dt; ϵ = solver.ϵ, γfact = γfact)
    return (; dyrel, γfact, η_mean)
end

"""
    build_model(setup)

Grid, phases, Stokes arrays, initial state and solver, ready for the first time step.
The returned named tuple also carries `setup`, so the step functions below need only `model`.
"""
function build_model(setup)
    (; domain, zone, material, rheology, loading, time, adjoint, solver) = setup
    (; dt) = time

    ni = domain.nx, domain.ny
    li = domain.lx, domain.ly
    grid = Geometry(ni, li; origin = domain.origin)
    phase_ratios = sharp_phase_ratios(grid, zone, length(rheology))

    (; flow_bcs, stokes, stokes_ad, gradients, observation) = allocate_stokes(grid, ni, length(rheology), adjoint)
    (; ρg, T) = initialise_material_state!(stokes, phase_ratios, rheology, material, ni, dt)
    apply_simple_shear!(stokes, flow_bcs, grid, loading.εbg)
    (; dyrel, γfact, η_mean) = build_dyrel(stokes, rheology, phase_ratios, grid, material, solver, dt)

    args = (; T = T, P = stokes.P, dt = dt, ΔT = @zeros(ni .+ 2...))

    return (;
        setup, grid, phase_ratios, flow_bcs, stokes, stokes_ad, gradients, observation,
        ρg, T, args, dyrel, γfact, η_mean,
    )
end

function print_model_summary(model)
    log_at(1) || return nothing
    (; setup, γfact, η_mean, observation) = model
    print_setup(setup)
    @printf("DYREL penalty: γfact = %.3e (γ_num = %.3e Pa s)\n", γfact, γfact * η_mean)
    @printf("Objective: %d Vx nodes in the box, adjoint = %s\n", count(!iszero, Array(observation.weights)), setup.adjoint.run)
    return nothing
end

## FORWARD + ADJOINT SOLVE ------------------------------------------------------------

"""
    solve_step!(model, igg; run_adjoint)

Forward DYREL solve of one time step and, with `run_adjoint`, the adjoint solve. Both are
one `solve_DYREL!` call: the adjoint runs on the converged forward state before the solver
advances the stress history (τ_o ← τ), so it cannot be split off as a later call.
Returns the solver output `out`, the objective `J` [m/s] and the solve wall time `t_solve`.
"""
function solve_step!(model, igg; run_adjoint)
    (; setup, stokes, stokes_ad, ρg, dyrel, flow_bcs, phase_ratios, args, grid, observation, gradients) = model
    (; rheology, solver, time) = setup
    t_solve = @elapsed out = solve_DYREL!(
        stokes, stokes_ad, ρg, dyrel, flow_bcs, phase_ratios, rheology, args, grid, time.dt, igg;
        solver.kwargs...,
        adjoint = run_adjoint,
        observation = run_adjoint ? observation : nothing,
        gradients = gradients,
    )
    J = objective_value(stokes, observation)
    return (; out, J, t_solve)
end

## GRADIENTS --------------------------------------------------------------------------

"""
    compute_sensitivity_totals(model, J)

Domain-integrated dJ/dp per phase and target parameter, from the last adjoint solve.
"""
function compute_sensitivity_totals(model, J)
    (; setup, gradients, stokes_ad, phase_ratios) = model
    (; rheology, adjoint, material) = setup
    totals = sensitivity_totals(gradients, stokes_ad, phase_ratios, rheology, adjoint.target_parameters)
    log_at(2) && print_sensitivity_totals(totals, J; phase_names = material.phase_names)
    return totals
end

"""
    analyse_sensitivities(model, it, t; plot)

Sensitivity fields of step `it` (`raw`, `logarithmic` p·∂J/∂p, `relative` (p/J)·∂J/∂p).
With `plot`, also saves the figures selected by `output.sensitivity_types`.
"""
function analyse_sensitivities(model, it, t; plot)
    (; setup, stokes, stokes_ad, gradients, grid, observation) = model
    (; material, time, output, zone) = setup
    types = output.sensitivity_types
    return run_sensitivity_analysis(
        stokes, stokes_ad, gradients, grid, observation, figure_dirs(output).adj;
        ρ_matrix = material.ρ, ρ_weak = material.ρ, G = material.G, K = material.Kb,
        zone_keywords(zone)..., step = it, t = t, τM = time.τM_weak,
        plot = plot ? types : (),
        combined = plot && :combined in types,
    )
end

"""
    print_final_sensitivities(history, material)

Sensitivities of the last adjoint step, plus dJ/dν and dJ/dG at fixed ν from dJ/dG and dJ/dKb.
"""
function print_final_sensitivities(history, material)
    isempty(history) && return nothing
    h = last(history)
    println("\nSensitivities at step $(h.step):")
    print_sensitivity_totals(h.totals, h.J; phase_names = material.phase_names)
    dG, dKb = (sum(r.dJdp for r in h.totals if r.name === n) for n in (:G, :Kb))
    nus = nu_sensitivities(dG, dKb, material.G, material.ν)
    @printf("  elasticity from (G, ν): dJ/dν = % .6e, dJ/dG at fixed ν = % .6e\n", nus.ν, nus.G_ν)
    return nothing
end

## PLOTTING ---------------------------------------------------------------------------

"""
    plot_absolute_sensitivities(model, it)

∂J/∂p maps of step `it`, one figure per style in `output.sensitivity_modes`.
"""
function plot_absolute_sensitivities(model, it)
    (; setup, stokes_ad, gradients, grid) = model
    (; output, adjoint, zone) = setup
    for mode in output.sensitivity_modes
        plot_sensitivities(
            stokes_ad, gradients, grid, figure_dirs(output).adj;
            mode = mode, zone_keywords(zone)...,
            observation_center = adjoint.center,
            observation_half_width = adjoint.half_width, step = it,
        )
    end
    return nothing
end

plot_forward_step(model, it, t) = plot_forward_fields(
    model.stokes, model.grid, it, t, figure_dirs(model.setup.output).fwd; zone_keywords(model.setup.zone)...
)

## TIME STEPPING ----------------------------------------------------------------------

function print_step(sol, it, nt, run_adjoint)
    log_at(1) || return nothing
    (; out, J, t_solve) = sol
    err = isempty(out.err_evo_tot) ? NaN : last(out.err_evo_tot)
    @printf(
        "Step %3d/%d%s | %7d iterations | err = %.2e | J = %.4f cm/yr | %.1f s\n",
        it, nt, run_adjoint ? " + adjoint" : "          ", out.iter, err,
        J * MS_TO_CMYR, t_solve
    )
    return nothing
end

"""
    time_step!(model, igg, it, history)

Time step `it`: solve, then on adjoint steps the gradients (appended to `history`) and
sensitivity figures, then the forward figures and the console line.
"""
function time_step!(model, igg, it, history)
    (; time, adjoint, output) = model.setup
    t = it * time.dt
    run_adjoint = adjoint.run && is_step(it, adjoint.every, time.nt)
    plot_results = output.plot_every !== nothing

    sol = solve_step!(model, igg; run_adjoint)

    if run_adjoint
        totals = compute_sensitivity_totals(model, sol.J)
        push!(history, (; step = it, t, J = sol.J, totals))
        analyse_sensitivities(model, it, t; plot = plot_results)
        plot_results && :absolute in output.sensitivity_types && plot_absolute_sensitivities(model, it)
    end
    plot_results && is_step(it, output.plot_every, time.nt) && plot_forward_step(model, it, t)

    print_step(sol, it, time.nt, run_adjoint)
    return sol
end

function run_time_loop!(model, igg)
    history = NamedTuple[]
    for it in 1:model.setup.time.nt
        time_step!(model, igg, it, history)
    end
    return history
end

## MAIN -------------------------------------------------------------------------------
function main(setup, igg)
    configure_output!(setup.output)
    model = build_model(setup)
    print_model_summary(model)
    prepare_figure_dirs(setup.output)

    history = run_time_loop!(model, igg)
    print_final_sensitivities(history, setup.material)

    (; stokes, stokes_ad, gradients, grid, phase_ratios, observation) = model
    return (; stokes, stokes_ad, gradients, grid, phase_ratios, observation, history)
end
## END OF MAIN ------------------------------------------------------------------------

# Run ---------------------------------------------------------------------------------
setup = setup_linear()

if ImplicitGlobalGrid.grid_is_initialized()
    finalize_global_grid(; finalize_MPI = false)
end
igg = IGG(init_global_grid(setup.domain.nx, setup.domain.ny, 1; init_MPI = !JustRelax.MPI.Initialized())...)

result = main(setup, igg)
nothing   # to keep `result` quiet / to not print result
