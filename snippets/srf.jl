# 0D code for Strongly Slip-Rate dependent Friction (SSRF)
# ? No need to set all to zero?

# Define parameters needed in functions ( - indicates name in I2ELVIS)
bbrit = 0.6;        # friction coefficient (-) - bbrit 
lamb  = 0.4;        # pore fluid pressure ratio    1-?
coh = 10e9;         # cohesion (MPA) - abrit 
relvw = 1;          # relative amount of Velocity-Weakening (VW) vs Velocity-Strengthening (VS); 1 = VW full, 0 = VS full
                    # could calculate as a function of temperature of depth
v_c = 1e-7;         # characteristic slip rate, at which half of friction change occurs
                    # in test set equal to v_slip, such that can see if friction is halfway 

sigin = 1e6;        # second invariant of deviatoric stress (Pa) 
epsin = 1e-10;      # second invariant of strain rate (/s)
pres = 1e9;         # pressure (Pa) - 
dx = 500;           # grid size (m) - res_high

v_slip = calcfriction(epsin,dx);

# Function to calculate SSRF
function calcfriction(epsin,dx)

    # Calculate local slip rate following van Dinther et al., 2013a
    v_slip = 2.0*epsin*dx;

    # 

    return vslip#, mu
end

# Function to calculate plastic strength 
#function calcstrength(mu,coh,lamb,pressure)


    #return strength
#end
