# calculate new (dimensionless) state, dynamic friction coefficient, and yield strength
function calcyield(a, b, L, mu0, V0, Vp, dt, oldstate, P, lambda=0, cohesion=0)
    if (Vp * dt / L <= 1e-6)
        newstate = log(exp(oldstate) * (1 - Vp * dt / L) + V0 * dt / L)
    else
        newstate = log(V0 / Vp + (exp(oldstate) - V0 / Vp) * exp(-Vp * dt / L))
    end
    mud = a * asinh(Vp / (2 * V0) * exp((mu0 + b * newstate) / a))
    syield = P * (1 - lambda) * mud + cohesion
    return newstate, mud, syield
end