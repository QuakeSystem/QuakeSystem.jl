# Prescribed viscosity field for the ηv inversion
#
# The parameter is ηv on the vertices (optimised as ln ηv). η on the centres follows from it as
# the harmonic mean of the four corner vertices. Both are written into `stokes` before the solve,
# and the solver runs with `update_material = false`: otherwise `solve_DYREL!` overwrites them
# from the rheology (compute_viscosity!) at the start of every time step. The rheology then
# only supplies the elasticity (G, Kb); its creep laws no longer set the viscosity.
#
# Periodic in x: vertex columns 1 and nx+1 are the same physical vertex (the seam), so they
# must hold the same value, and their gradients are summed.

"""
    setup_ηv(steps_per_tr; run_adjoint = false, η_matrix_ref = 1.0e25, η_weak_ref = 1.0e18)

Setup for runs with a prescribed ηv. The time step comes from the reference viscosities, so
it stays fixed while ηv changes (as in `setup_G`).
"""
function setup_ηv(steps_per_tr; run_adjoint = false, η_matrix_ref = 1.0e25, η_weak_ref = 1.0e18)
    s = setup_linear(; η_matrix = η_matrix_ref, η_weak = η_weak_ref, steps_per_tr,
                     plot_every = nothing, verbose = 1, run_adjoint)
    solver = merge(s.solver, (; kwargs = merge(s.solver.kwargs, (; update_material = false))))
    return merge(s, (; solver))
end

"""
    rheology_viscosity(; η_matrix, η_weak)

The (ηv, η) the two-phase rheology gives: the starting point of the inversion, as CPU arrays.
"""
function rheology_viscosity(; η_matrix, η_weak)
    s = setup_linear(; η_matrix, η_weak, plot_every = nothing, verbose = 0, run_adjoint = false)
    model = build_model(s)
    return (; ηv = Array(model.stokes.viscosity.ηv), ηc = Array(model.stokes.viscosity.η))
end

# ηv[1, :] and ηv[nx+1, :] are the same vertex on the periodic seam
enforce_seam!(ηv) = (@views ηv[end, :] .= ηv[1, :]; ηv)

# The seam vertex is one parameter: sum both gradient entries and give both columns the total,
# so a step on the full array keeps the seam consistent
fold_seam!(g) = (@views g[1, :] .+= g[end, :]; @views g[end, :] .= g[1, :]; g)

"""
    vertex2center_harm(ηv)

η on the centres: harmonic mean of the four corner vertices of each cell.
"""
function vertex2center_harm(ηv)
    nx, ny = size(ηv) .- 1
    return [4 / (1 / ηv[i, j] + 1 / ηv[i + 1, j] + 1 / ηv[i, j + 1] + 1 / ηv[i + 1, j + 1])
            for i in 1:nx, j in 1:ny]
end

"""
    vertex2center_harm_adjoint(ηv, ηc, gc)

Pull a centre sensitivity gc = ∂J/∂η back to the vertices through `vertex2center_harm`:
∂ηc/∂ηv_k = ηc² / (4 ηv_k²) for each of the four corners k of a cell.
"""
function vertex2center_harm_adjoint(ηv, ηc, gc)
    gv = zero(ηv)
    nx, ny = size(ηc)
    for j in 1:ny, i in 1:nx
        w = gc[i, j] * ηc[i, j]^2 / 4
        for (a, b) in ((i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1))
            gv[a, b] += w / ηv[a, b]^2
        end
    end
    return gv
end

"""
    set_viscosity!(model, ηv; ηc = nothing)

Write ηv (seam enforced, clamped to the viscosity cutoff) and η on the centres into
`model.stokes`, then rebuild the DYREL state from them. `ηc` defaults to
`vertex2center_harm(ηv)`; pass it only to prescribe both fields (e.g. for the sanity check).
Returns the updated model: use the returned one, it carries the new `dyrel`.
"""
function set_viscosity!(model, ηv; ηc = nothing)
    (; setup, stokes, phase_ratios, grid) = model
    lo, hi = setup.material.viscosity_cutoff
    ηv = enforce_seam!(clamp.(ηv, lo, hi))
    ηc = clamp.(something(ηc, vertex2center_harm(ηv)), lo, hi)
    copyto!(stokes.viscosity.ηv, ηv)
    copyto!(stokes.viscosity.η, ηc)
    # γ_num = γfact·mean(η) and the preconditioner depend on η
    dyrel_state = build_dyrel(stokes, setup.rheology, phase_ratios, grid, setup.material, setup.solver, setup.time.dt)
    return merge(model, dyrel_state)
end

"""
    forward_Vx_ηv(setup, igg, ηv; ηc = nothing)

Forward run (no adjoint) with the prescribed viscosity; `setup` from `setup_ηv`.
Returns the final Vx and the box mask, like `forward_Vx`.
"""
function forward_Vx_ηv(setup, igg, ηv; ηc = nothing)
    setup.solver.kwargs.update_material && error("setup must have update_material = false (use setup_ηv)")
    configure_output!(setup.output)
    model = set_viscosity!(build_model(setup), ηv; ηc)
    run_time_loop!(model, igg)
    (; center, half_width) = setup.adjoint
    box = box_observation(model.stokes.V.Vx, model.grid, center, half_width; antisymmetric = false)
    return Array(model.stokes.V.Vx), Array(box.weights) .!= 0
end

"""
    misfit_gradient_lnηv(setup, igg, ηv, Vx, Vx_obs, mask)

∂m/∂ln ηv on the vertices, with η on the centres derived from ηv:
∂m/∂ηv = (∂m/∂ηv)_vertex + vertex2center_harm_adjoint(∂m/∂η)_centre, seam folded,
then times ηv. Like `misfit_gradient_G`, it is the gradient of the last time step with the
stress history of the earlier steps frozen.
"""
function misfit_gradient_lnηv(setup, igg, ηv, Vx, Vx_obs, mask)
    s = merge(setup, (; adjoint = merge(setup.adjoint, (; run = true))))
    configure_output!(s.output)
    model = set_viscosity!(build_model(s), ηv)
    model.observation.weights .= mask .* (Vx .- Vx_obs)
    run_time_loop!(model, igg)
    ηv_used = Array(model.stokes.viscosity.ηv)
    ηc = Array(model.stokes.viscosity.η)
    gv = Array(model.stokes_ad.viscosity.ηv)
    gc = Array(model.stokes_ad.viscosity.η)
    dmdηv = fold_seam!(gv .+ vertex2center_harm_adjoint(ηv_used, ηc, gc))
    return ηv_used .* dmdηv
end
