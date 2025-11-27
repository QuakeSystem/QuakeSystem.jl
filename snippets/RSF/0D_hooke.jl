# assume constant stressing rate
# THIS IMPLEMENTATION IS WRONG BUT IT SOMEHOW WORKS REALLY WELL
using Plots
include("calcyield_rsf.jl")
include("calcvpvisc.jl")
include("calc_sliprate_after.jl")
include("calc_dt.jl")

# numerical paramters
dx = 500
dy = 500

# material parameters
a = 0.011 # rate parameter
b = 0.015 # state paramter
L = 0.05 # characteristic slip distance [m]
mu0 = 0.5 # reference friction
V0 = 4e-9 # reference slip rate [m/s]
lambda = 0 # fluid pressure ratio
cohesion = 0 # cohesion [Pa]
poisson = 0.25 # poisson ratio
G = 20e9 # shear modulus [Pa]

# constant variables
#Vp = 1e-8
P = 50e6 # pressure [Pa]
etav = 1e24 # ductile viscosity [Pa]
vy = 0

# stressing rate in 0D
# simplified hookes law: stress rate = 2*G*strainrate
# strainrate = sliprate/(2*dx)
platerate = 4e-9
stressrate = 2 * G * platerate / (2 * dx)

# initial velocity and state
Vp = V0
state = 0
# initial viscoplastic viscosity
etavp = etav
# initial stress
tau = 0e6
function timestep(oldstate, Vp, etavp, tau)
    vx = Vp / 2 # update vx for dt calculation (assume all slip is on x-axis)
    dt = calc_dt(P, lambda, G, etavp, L, V0, poisson, dx, dy, vx, vy, Vp, oldstate)
    # add stress
    tau = tau + stressrate * dt
    newstate, mud, syield = calcyield(a, b, L, mu0, V0, Vp, dt, oldstate, P, lambda, cohesion)
    etavp = calcvpvisc(syield, Vp, etav, dx, dy)
    # NAVIER STOKES WOULD BE SOLVED HERE
    Vp = calc_sliprate_RSF(a, b, mu0, V0, newstate, P, tau, lambda, cohesion)
    # reduce stress according to slip
    # this would need to be multiplied by dt to make sense, but that doesn't work while this does
    taudiff = -2 * G * Vp / (2 * dx)
    #println(taudiff)
    tau = tau + taudiff
    return dt, newstate, mud, etavp, Vp, tau
end

# time step count
niter = 1500
#niter = 80000

# initialize arrays
Vps = zeros(niter)
states = zeros(niter)
muds = zeros(niter)
dts = zeros(niter)
taus = zeros(niter)

for i = 1:niter
    global state, Vp, etavp, tau
    dt, state, mud, etavp, Vp, tau = timestep(state, Vp, etavp, tau)
    Vps[i] = Vp
    states[i] = state
    muds[i] = mud
    dts[i] = dt
    taus[i] = tau
end

times = cumsum(dts)

# interface strength
istrength = P .* (mu0 .+ b * states)

#plot(times, states)
p = plot(times[1:end] / 3600 / 24 / 365.25, Vps[1:end], yscale=:log10)
xlabel!("Time in years")
ylabel!("Vp in m/s")
display(p)
p = plot(times[1:end] / 3600 / 24 / 365.25, taus[1:end] / 1e6, label="stress")
p = plot!(times[1:end] / 3600 / 24 / 365.25, istrength[1:end] / 1e6, label="interface strength")
xlabel!("Time in years")
ylabel!("stress in MPa")
display(p)
p = plot(dts, yscale=:log10)
xlabel!("Time step no.")
ylabel!("dt in s")
display(p)