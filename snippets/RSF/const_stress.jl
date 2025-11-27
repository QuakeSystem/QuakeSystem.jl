#force slip by applying constant stress (fun expontential growth!)
using Plots
include("calcyield_rsf.jl")
include("calcvpvisc.jl")
include("calc_sliprate_after.jl")
#cd(@__DIR__)

# numerical paramters
dx = 500
dy = 500

# material parameters
a = 0.011
b = 0.015
L = 0.09
mu0 = 0.5
V0 = 2e-11
lambda = 0
cohesion = 0

# variables
#Vp = 1e-8
P = 50e6
etav = 1e24

# constant stress
tauii = 26e6

# initial velocity and state
Vp = V0
state = 1

# time step duration
dt = 100000

function timestep(oldstate, Vp)
    newstate, mud, syield = calcyield(a, b, L, mu0, V0, Vp, dt, oldstate, P, lambda, cohesion)
    etavp = calcvpvisc(syield, Vp, etav, dx, dy)
    # NAVIER STOKES WOULD BE SOLVED HERE
    #println(etavp)

    Vp = calc_sliprate_RSF(a, b, mu0, V0, newstate, P, tauii, lambda, cohesion)
    return newstate, mud, Vp
end

# time step count
niter = 25500


# initialize arrays
Vps = zeros(niter)
states = zeros(niter)
muds = zeros(niter)

for i = 1:niter
    global state, Vp
    state, mud, Vp = timestep(state, Vp)
    Vps[i] = Vp
    states[i] = state
    muds[i] = mud
end

plot(states)
plot(Vps)
plot(muds[2:end])