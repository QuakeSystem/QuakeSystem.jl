# Rheology of the shear-zone models: phase 1 = matrix, phase 2 = weak zone

"""
    init_rheology_linear(; η_matrix, η_weak, G, ν, ρ, g)

Linear viscoelastic rheology. The creep law is DislocationCreep with n = 1, A = 1/(3η) and
E = V = 0, which is a linear viscosity η, so ∂J/∂A and ∂J/∂V are available as in the
power-law model. Elasticity from (G, ν); the solver uses the bulk modulus Kb derived from them.
"""
function init_rheology_linear(; η_matrix, η_weak, G, ν, ρ, g)
    elasticity = ConstantElasticity(; G = G, ν = ν)
    creep_matrix = DislocationCreep(; n = 1.0, r = 0.0, A = 1 / (3η_matrix), E = 0.0, V = 0.0)
    creep_weak = DislocationCreep(; n = 1.0, r = 0.0, A = 1 / (3η_weak), E = 0.0, V = 0.0)
    return (
        material_phase("Matrix", 1, creep_matrix, elasticity, ρ, g),
        material_phase("Weak zone", 2, creep_weak, elasticity, ρ, g),
    )
end

function material_phase(name, id, creep, elasticity, ρ, g)
    return SetMaterialParams(;
        Name = name, Phase = id,
        Density = ConstantDensity(; ρ = ρ),
        CompositeRheology = CompositeRheology((creep, elasticity)),
        Elasticity = elasticity,
        Gravity = ConstantGravity(; g = g),
    )
end

# Bulk modulus the solver uses (the same in every phase)
bulk_modulus(rheology) = Float64(rheology[1].Elasticity[1].Kb.val)
