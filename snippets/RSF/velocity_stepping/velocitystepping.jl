# Synthetic velocity stepping experiment as initial test of RSF state evolution and yield strength calculation
#
using Plots

# simpler version of function using only one analytical solution for state update
function calcyield(a, b, L, mu0, V0, Vp, dt, oldstate, P, lambda=0, cohesion=0)
    newstate = log(V0 / Vp + (exp(oldstate) - V0 / Vp) * exp(-Vp * dt / L))
    mud = a * asinh(Vp / (2 * V0) * exp((mu0 + b * newstate) / a))
    syield = P * (1 - lambda) * mud + cohesion
    return newstate, mud, syield
end
# actual version with case seperation for simpler anyltical solution for low Vp
function calcyield2(a, b, L, mu0, V0, Vp, dt, oldstate, P, lambda=0, cohesion=0)
    if (Vp * dt / L <= 1e-6)
        newstate = log(exp(oldstate) * (1 - Vp * dt / L) + V0 * dt / L)
    else
        newstate = log(V0 / Vp + (exp(oldstate) - V0 / Vp) * exp(-Vp * dt / L))
    end
    mud = a * asinh(Vp / (2 * V0) * exp((mu0 + b * newstate) / a))
    syield = P * (1 - lambda) * mud + cohesion
    return newstate, mud, syield
end
#phi = 5
#phi, syield = calcyield(0.011, 0.015, 0.09, 0.5, 2e-11, 2e-11, 1e5, phi, 50e6);
#print(phi)
# initiate arrays
phi = zeros(10000)
syield = zeros(10000)
mud = zeros(10000)
slip = zeros(10000)
phi[1] = 5
slip[1] = 0
# steady state creep
for i = 1:100
    phi[i+1], mud[i+1], syield[i+1] = calcyield2(0.011, 0.015, 0.09, 0.5, 2e-11, 2e-11, 5e8, phi[i], 50e6)
end
slip[2:100] = slip[1].+2e-11*5e8*1:99
# velocity step
for i = 101:2000
    phi[i+1], mud[i+1], syield[i+1] = calcyield2(0.011, 0.015, 0.09, 0.5, 2e-11, 1, 1e-2, phi[i], 50e6)
end
slip[101:2000] = slip[100] .+ 1 * 1e-2 * (1:1900)
# arrest
for i = 2001:3000
    phi[i+1], mud[i+1], syield[i+1] = calcyield2(0.011, 0.015, 0.09, 0.5, 2e-11, 1e-10, 1e7, phi[i], 50e6)
end
slip[2001:3000] = slip[2000] .+ 2e-11 * 1e7 * (1:1000)

p = Plots.plot(slip[2:3000], mud[2:3000])
xlabel!("Slip in m")
ylabel!("μ")
display(p)
p = Plots.plot(slip[2:3000], phi[2:3000])
xlabel!("Slip in m")
ylabel!("Ω")
display(p)