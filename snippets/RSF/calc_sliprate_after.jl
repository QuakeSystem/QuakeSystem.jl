# recalculate plastic slip rate after solving N-S
begin
    function calc_sliprate_RSF(a, b, mu0, V0, state, P, tauii, lambda=0, cohesion=0)
        if tauii > 0
            Vp_after = 2 * V0 * sinh(max((tauii - cohesion), 0) / (a * P * (1 - lambda))) * exp(-(mu0 + b * state) / a)
        else
            Vp_after = 0
        end
        return Vp_after
    end
    function calc_sliprate_RSF_qd(a, b, mu0, V0, state, P, tauii, G, rho, lambda=0, cohesion=0)
        cs = sqrt(G / rho) # shear wave velocity
        etaqd = G / (2 * cs) # radiation damping term
        if tauii > 0
            Vp_after = 2 * V0 * sinh(max((tauii - cohesion - etaqd), 0) / (a * P * (1 - lambda))) * exp(-(mu0 + b * state) / a)
        else
            Vp_after = 0
        end
        return Vp_after
    end
end