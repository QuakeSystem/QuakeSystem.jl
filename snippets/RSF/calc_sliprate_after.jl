# recalculate plastic slip rate after solving N-S
function calc_sliprate_RSF(a, b, mu0, V0, state, P, tauii, lambda=0, cohesion=0)
    if tauii > 0
        Vp_after = 2 * V0 * sinh(max((tauii - cohesion), 0) / (a * P * (1 - lambda))) * exp(-(mu0 + b * state) / a)
    else
        Vp_after = 0
    end
    return Vp_after
end