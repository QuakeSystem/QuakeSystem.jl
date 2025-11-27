#force slip by applying constant stress (fun expontential growth!)
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
b = 0.012 # state paramter
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


# constant stress
tauii = 26e6

# initial velocity and state
Vp = V0
state = 1
# initial viscoplastic viscosity
etavp = etav
# initial time step duration
dt = 100000

function timestep(oldstate, Vp, etavp)
    vx = Vp / 2 # update vx for dt calculation (assume all slip is on x-axis)
    dt = calc_dt(P, lambda, G, etavp, L, V0, poisson, dx, dy, vx, vy, Vp, oldstate)
    newstate, mud, syield = calcyield(a, b, L, mu0, V0, Vp, dt, oldstate, P, lambda, cohesion)
    etavp = calcvpvisc(syield, Vp, etav, dx, dy)
    # NAVIER STOKES WOULD BE SOLVED HERE

    Vp = calc_sliprate_RSF(a, b, mu0, V0, newstate, P, tauii, lambda, cohesion)
    return dt, newstate, mud, etavp, Vp
end

# time step count
niter = 100


# initialize arrays
Vps = zeros(niter)
states = zeros(niter)
muds = zeros(niter)
dts = zeros(niter)

for i = 1:niter
    global state, Vp, etavp
    dt, state, mud, etavp, Vp = timestep(state, Vp, etavp)
    Vps[i] = Vp
    states[i] = state
    muds[i] = mud
    dts[i] = dt
end

times = cumsum(dts)

#plot(times, states)
plot(times, Vps, yscale=:log10)
xlabel!("Time in s")
ylabel!("Vp in m/s")

plot(dts, yscale=:log10)
#plot(muds[2:end])