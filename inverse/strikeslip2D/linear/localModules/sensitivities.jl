# Per-phase totals of the adjoint gradients

# Value of the rheology parameter `name` (:G, :Kb, :A, :V, ...) of one phase, or nothing
function parameter_value(phase, name)
    for law in phase.CompositeRheology[1].elements
        hasproperty(law, name) && return Float64(getproperty(law, name).val)
    end
    return nothing
end

"""
    sensitivity_totals(gradients, stokes_ad, phase_ratios, rheology, target_parameters)

∂J/∂p for a uniform change of phase k's parameter p: the sum of its cell gradients (the
vertex contributions are already folded into the cells). Also p·∂J/∂p = ∂J/∂ln p. ∂J/∂ρ
is split over the phases by their cell fractions. Returns a vector of named tuples
(name, phase, p, dJdp, pdJdp).
"""
function sensitivity_totals(gradients, stokes_ad, phase_ratios, rheology, target_parameters)
    rows = NamedTuple[]
    for name in target_parameters, (k, phase) in enumerate(rheology)
        p = parameter_value(phase, name)
        isnothing(p) && continue
        dJdp = sum(Array(@view gradients[name].center[k, :, :]))
        push!(rows, (; name, phase = k, p, dJdp, pdJdp = p * dJdp))
    end
    dJdρ = Array(stokes_ad.ρ)
    ratios = phase_ratios.center          # CPU CellArray: ratios[i, j][k]
    for (k, phase) in enumerate(rheology)
        ρ = Float64(phase.Density[1].ρ.val)
        dJdp = sum(dJdρ[I] * ratios[I][k] for I in CartesianIndices(dJdρ))
        push!(rows, (; name = :ρ, phase = k, p = ρ, dJdp, pdJdp = ρ * dJdp))
    end
    return rows
end

"""
    nu_sensitivities(dJdG, dJdKb, G, ν)

Elasticity given by (G, ν) instead of (G, Kb), Kb = 2G(1+ν)/(3(1-2ν)). The solver only uses
Kb, so the adjoint gives ∂J/∂Kb and ∂J/∂G at fixed Kb; by the chain rule
  ∂J/∂ν        = ∂J/∂Kb · 2G/(1-2ν)²
  ∂J/∂G at ν   = ∂J/∂G at Kb + ∂J/∂Kb · Kb/G
"""
function nu_sensitivities(dJdG, dJdKb, G, ν)
    Kb = 2G * (1 + ν) / (3 * (1 - 2ν))
    return (; ν = dJdKb * 2G / (1 - 2ν)^2, G_ν = dJdG + dJdKb * Kb / G)
end

function print_sensitivity_totals(rows, J; phase_names)
    @printf("  J = %.6e m/s = %.4f cm/yr\n", J, J * MS_TO_CMYR)
    for r in rows
        @printf(
            "  %-3s %-16s p = % .4e   dJ/dp = % .6e   p·dJ/dp = % .6e cm/yr\n",
            r.name, phase_names[r.phase], r.p, r.dJdp, r.pdJdp * MS_TO_CMYR
        )
    end
    return nothing
end
