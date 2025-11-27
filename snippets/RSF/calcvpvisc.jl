function calcvpvisc(syield, Vp, etav, dx, dy)
    faultwidth = 0.5 * (dx + dy)
    etavp = etav * syield / (2 * etav * Vp / faultwidth + syield)
    return etavp
end

