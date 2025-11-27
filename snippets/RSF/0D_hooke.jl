# assume constant stressing rate
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
L = 0.09 # characteristic slip distance [m]
mu0 = 0.5 # reference friction
V0 = 2e-11 # reference slip rate [m/s]
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
platerate = 1e-9
stressrate = 2 * G * platerate / (2 * dx)

# initial velocity and state
Vp = V0
state = 1
# initial viscoplastic viscosity
etavp = etav
# initial stress
tau = 0

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
    taudiff = -2 * G * Vp / (2 * dx)
    tau = tau + taudiff
    return dt, newstate, mud, etavp, Vp, tau
end

# time step count
niter = 1000

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

#plot(times, states)
p = plot(times / 3600 / 24 / 365.25, Vps, yscale=:log10)
xlabel!("Time in years")
ylabel!("Vp in m/s")
display(p)
p = plot(times / 3600 / 24 / 365.25, taus / 1e6)
xlabel!("Time in years")
ylabel!("stress in MPa")
display(p)