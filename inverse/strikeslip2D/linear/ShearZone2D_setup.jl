# Linear viscoelastic shear zone in periodic simple shear: every model setting.
# The setup only holds data; ShearZone2D.jl builds and runs the model from it.
# Derived quantities (εbg, Kb, t_r, dt, nt) are computed here from the primary ones.

"""
    setup_linear(; figdir, kwargs...)

Settings of the linear shear-zone model, as one named tuple of groups:
`domain`, `zone`, `material`, `rheology`, `loading`, `time`, `adjoint`, `solver`, `output`.
Keyword arguments override the defaults below, so a parameter sweep is a loop over
`setup_linear(; η_weak = ...)` without editing this file.
"""
function setup_linear(;
        # Domain and grid -----------------------------------------------------------
        lx = 1000.0, ly = 1000.0,       # [m], origin at (-lx/2, -ly/2)
        nx = 40, ny = 80,               # dx = 25 m, dy = 12.5 m: zone = 8 cells thick
        # Weak zone -----------------------------------------------------------------
        zone_x = (-500.0, 500.0),       # [m]: spans the periodic domain in x
        zone_y = (-50.0, 50.0),
        # Material ------------------------------------------------------------------
        η_matrix = 1.124e25,            # [Pa s]
        η_weak = 2.171e18,
        G = 50.0e9,                     # [Pa]
        ν = 0.45,                       # Kb follows from G and ν
        ρ = 3300.0,                     # [kg/m³], same in both phases
        g = 9.81,                       # [m/s²], gravity along -y
        depth = 40.0e3,                 # [m]: background pressure P_bg = ρ g depth
        T = 1.0,                        # [K]: no effect on η (E = 0), but ∂J/∂V ∝ P/(nRT)
        viscosity_cutoff = (1.0e17, 1.0e30),
        # Loading -------------------------------------------------------------------
        VEL = 6.0e-2 / SEC_PER_YR,      # [m/s]: top - bottom velocity difference
        # Time stepping -------------------------------------------------------------
        N_tr = 1.0,                     # run time in elastic loading times t_r
        steps_per_tr = 10,
        # Adjoint -------------------------------------------------------------------
        run_adjoint = true,
        adjoint_every = nothing,        # adjoint every n steps; nothing → last step only
        observation_center = (0.0, 0.0),
        observation_half_width = (100.0, 100.0),
        antisymmetric = true,           # J = Σ Vx(upper half) - Σ Vx(lower half) of the box
        target_parameters = (:G, :Kb, :A, :V),   # parameters to differentiate
        # Solver --------------------------------------------------------------------
        ϵ = 1.0e-7,
        γfact = nothing,                # nothing → γ_num = 20 · mean(η_ve)
        # Output --------------------------------------------------------------------
        figdir = joinpath(@__DIR__, "Figures_linear"),
        plot_every = 1,                 # forward figures every n steps; nothing → none
        # sensitivity figures at every adjoint step, any subset of
        #   :absolute    ∂J/∂p              (adjoint/dJdp_heatmap, dJdp_contour)
        #   :logarithmic p·∂J/∂p            (adjoint/logarithmic, logarithmic_shared)
        #   :relative    (p/J)·∂J/∂p        (adjoint/relative, relative_shared, max_relative)
        #   :combined    ∂J/∂p and p·∂J/∂p in one figure (adjoint/combined)
        sensitivity_types = (:absolute, :logarithmic, :relative, :combined),
        sensitivity_modes = (:heatmap, :contour),   # style of the :absolute figures
        verbose = 1,                    # 0 silent, 1 one line per step, 2 solver log
        show = false,                   # also display every saved figure
    )

    domain = (; lx, ly, nx, ny, origin = (-lx / 2, -ly / 2))
    zone = (; x = zone_x, y = zone_y)

    # GeoParams laws in ShearZone2D_rheology.jl
    rheology = init_rheology_linear(; η_matrix, η_weak, G, ν, ρ, g)
    Kb = bulk_modulus(rheology)
    material = (;
        η_matrix, η_weak, G, ν, Kb, ρ, g, T, P_bg = ρ * g * depth, viscosity_cutoff,
        phase_names = ("matrix", "weak zone"),
    )

    loading = (; VEL, εbg = VEL / (2ly))

    # Elastic loading time t_r = ly/(G S), S = W/η_weak + (ly - W)/η_matrix
    W = zone_y[2] - zone_y[1]
    t_r = ly / (G * (W / η_weak + (ly - W) / η_matrix))
    dt = t_r / steps_per_tr
    time = (; N_tr, steps_per_tr, t_r, dt, nt = round(Int, N_tr * steps_per_tr), τM_weak = η_weak / G)

    adjoint = (;
        run = run_adjoint, every = adjoint_every,
        center = observation_center, half_width = observation_half_width,
        antisymmetric, target_parameters,
    )

    solver = (;
        ϵ, γfact,
        kwargs = (;
            verbose_PH = verbose >= 2, verbose_DR = false,
            iterMax = 5.0e5, iterMax_PH = 1.0e3, total_iterMax = Inf, nout = 100,
            rel_drop = 1.0e-1,
            λ_relaxation_PH = 1, λ_relaxation_DR = 1,
            linear_viscosity = true, viscosity_relaxation = 1, η_multiplier = nothing,
            viscosity_cutoff,
        ),
    )

    output = (; figdir, plot_every, sensitivity_types, sensitivity_modes, verbose, show)

    return (; domain, zone, material, rheology, loading, time, adjoint, solver, output)
end

function print_setup(s)
    (; material, time, loading) = s
    η_ve(η) = 1 / (1 / η + 1 / (material.G * time.dt))
    @printf(
        "Grid %d × %d, domain %g × %g m; weak zone x ∈ %s, y ∈ %s m\n",
        s.domain.nx, s.domain.ny, s.domain.lx, s.domain.ly, s.zone.x, s.zone.y
    )
    @printf(
        "η_matrix = %.3e, η_weak = %.3e Pa s; G = %.3e Pa, ν = %.2f (Kb = %.3e Pa)\n",
        material.η_matrix, material.η_weak, material.G, material.ν, material.Kb
    )
    @printf(
        "τ_M = %.3g yr, t_r = %.3g yr; nt = %d, Δt = %.3g yr = %.3g τ_M; ΔV = %.2f cm/yr\n",
        time.τM_weak / SEC_PER_YR, time.t_r / SEC_PER_YR, time.nt,
        time.dt / SEC_PER_YR, time.dt / time.τM_weak, loading.VEL * MS_TO_CMYR
    )
    @printf(
        "η_ve matrix = %.3e, η_ve weak = %.3e Pa s (contrast %.1f)\n",
        η_ve(material.η_matrix), η_ve(material.η_weak),
        η_ve(material.η_matrix) / η_ve(material.η_weak)
    )
    return nothing
end
