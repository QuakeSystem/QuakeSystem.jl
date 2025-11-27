# Individual time step calculation functions (seperated for ease of testing and reading)
begin
    # displacement time step (2D)
    function calc_dtd(vx, dx, vy, dy, ddmax=1e-3)
        # inverse of formulation from paper because vx/dx is more stable for very small vx
        maxdr = max(abs(vx / dx), abs(vy / dy))
        dtd = ddmax * 1.0 / maxdr
        return dtd
    end

    function max_state_change(P, lambda, G, poisson, dx, dy)
        # following Lapusta et al. (2000)
        Peff = P * (1 - lambda)
        faultwidth = 0.5 * (dx + dy)
        k = 2.0 / pi * (G / (1 - poisson)) / faultwidth
        xi = 1.0 / 4.0 * (k * L / (a * Peff) - (b - a) / a)^2 - k * L / (a * Peff)
        if xi > 0
            thetamax = a * Peff / (k * L - (b - a) * Peff)
        elseif xi < 0
            thetamax = 1 - (b - a) * Peff / (k * L)
        else
            error("xi=0 in time step calculation")
        end
        thetamax = min(thetamax, 0.2) # limit maximum according to Lapusta and Liu (2009)
        thetamax = max(thetamax, 0.1) # this is not mentioned in the paper but it's in the code
        return thetamax
    end
    # weakening time step
    function calc_dtw(thetamax, L, Vp)
        # STM-RSF limits Vp to 1e-4 for this calculation
        dtw = thetamax * L / Vp
        return dtw
    end
    # healing time step
    function calc_dth(thetamax)
        dth = 0.2 * thetamax
        return dth
    end

    # viscoelastic relaxation time
    function calc_dvep(etavp, G, fmax=0.2)
        dvep = fmax * etavp / G
        return dvep
    end
end