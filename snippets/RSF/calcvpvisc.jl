# calculate viscoplastic viscosity (2/2 could be simplified but it makes it less readable)
function calcvpvisc(syield, Vp, etav, dx, dy)
    faultwidth = 0.5 * (dx + dy)
    etavp = etav * syield / (2 * etav * Vp / (2 * faultwidth) + syield)
    return etavp
end

