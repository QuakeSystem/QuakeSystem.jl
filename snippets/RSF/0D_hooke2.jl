# assume constant stressing rate
using Plots
include("calcyield_rsf.jl")
include("calcvpvisc.jl")
include("calc_sliprate_after.jl")
include("calc_dt.jl")

# numerical paramters
dx = 500
dy = 500

# material parameters (has to be very specific or it decays/explodes, no adequate hinderance to extreme slip rates (add qd or inertia?))
a = 0.011 # rate parameter
b = 0.015 # state paramter
L = 0.0047 # characteristic slip distance [m]
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
# stress = 2*G*strain
# strain = displacement/(2*dx)
platerate = 4e-9

# initial velocity and state
Vp = V0 / 10
state = 0
# initial viscoplastic viscosity
etavp = etav
# initial stress
tau = 0e6

# initial displacement
x = 0.5

function timestep(oldstate, Vp, etavp, tau, x)
    vx = Vp / 2 # update vx for dt calculation (assume all slip is on x-axis)
    dt = calc_dt(P, lambda, G, etavp, L, V0, poisson, dx, dy, vx, vy, Vp, oldstate)
    if dt < 1e-4
        dt = 1e-4
    elseif dt > 1e7
        dt = 1e7
    end
    # predefine variables
    newstate = nothing
    mud = nothing
    syield = nothing
    negative = false

    while true
        # move block
        x = x + platerate * dt
        # compute stess
        tau = 2 * G * x / (2 * dx)
        if tau < 0
            println("negative")
            tau = abs(tau)
            negative = true
        end
        #println(tau)
        newstate, mud, syield = calcyield(a, b, L, mu0, V0, Vp, dt, oldstate, P, lambda, cohesion)
        etavp = calcvpvisc(syield, Vp, etav, dx, dy)
        # NAVIER STOKES WOULD BE SOLVED HERE
        Vp = calc_sliprate_RSF(a, b, mu0, V0, newstate, P, tau, lambda, cohesion)
        dtnew = calc_dt(P, lambda, G, etavp, L, V0, poisson, dx, dy, vx, vy, Vp, oldstate)
        # println(dt - dtnew)
        #if dtnew * 1.1 < dt
        #    println("reset")
        #    # reset displacement
        #    x = x - platerate * dt
        #    dt = dt * 0.5
        #if abs(Vp*dt)>0.4

        #else
        break
        #end
    end

    #println(x)
    # update position
    if negative == true
        x = x + Vp * dt
    else
        x = x - Vp * dt
    end
    # recompute stress
    #tau = 2 * G * x / (2 * dx)
    #println(x)
    #println(Vp)
    return dt, newstate, mud, etavp, Vp, tau, x
end

# time step count
niter = 5000

# initialize arrays
Vps = zeros(niter)
states = zeros(niter)
muds = zeros(niter)
dts = zeros(niter)
taus = zeros(niter)
xs = zeros(niter)

for i = 1:niter
    global state, Vp, etavp, tau, x
    dt, state, mud, etavp, Vp, tau, x = timestep(state, Vp, etavp, tau, x)
    #println(tau)
    Vps[i] = Vp
    states[i] = state
    muds[i] = mud
    dts[i] = dt
    taus[i] = tau
    xs[i] = x
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
#p = plot!(times[1:end] / 3600 / 24 / 365.25, istrength[1:end] / 1e6, label="interface strength")
xlabel!("Time in years")
ylabel!("stress in MPa")
display(p)
p = plot(dts, yscale=:log10)
xlabel!("Time step no.")
ylabel!("dt in s")
display(p)
p = plot(times / 3600 / 24 / 365.25, xs)
xlabel!("Time in years")
ylabel!("displacement in m")
display(p)