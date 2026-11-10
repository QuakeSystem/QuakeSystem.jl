# Weak-zone phases and the adjoint observation box

in_zone(x, y, zone) = zone.x[1] <= x <= zone.x[2] && zone.y[1] <= y <= zone.y[2]

# Sharp phase field: a cell is weak (phase 2) iff its centre lies in the zone; a vertex gets
# the mean weak fraction of its four cells (periodic in x, clamped in y). With harmonic phase
# averaging this gives ηv = harmonic mean of the cell viscosities.
@parallel_indices (i, j) function _sharp_center_phases!(ratios, wc)
    @index ratios[1, i, j] = 1.0 - wc[i, j]
    @index ratios[2, i, j] = wc[i, j]
    return nothing
end

@parallel_indices (i, j) function _sharp_vertex_phases!(ratios, wc)
    nx, ny = size(wc)
    iW, iE = i == 1 ? nx : i - 1, i == nx + 1 ? 1 : i
    jS, jN = max(j - 1, 1), min(j, ny)
    w = 0.25 * (wc[iW, jS] + wc[iE, jS] + wc[iW, jN] + wc[iE, jN])
    @index ratios[1, i, j] = 1.0 - w
    @index ratios[2, i, j] = w
    return nothing
end

function sharp_phase_ratios(grid, zone, nphases)
    ni = size(grid.xci[1], 1), size(grid.xci[2], 1)
    phase_ratios = PhaseRatios(backend_JP, nphases, ni)
    xc, yc = grid.xci
    wc = PTArray(backend)([in_zone(x, y, zone) ? 1.0 : 0.0 for x in xc, y in yc])
    @parallel (@idx ni) _sharp_center_phases!(phase_ratios.center, wc)
    @parallel (@idx ni .+ 1) _sharp_vertex_phases!(phase_ratios.vertex, wc)
    return phase_ratios
end

"""
    box_observation(Vx, grid, center, half_width; antisymmetric = true)

J = Σ w·Vx over the Vx nodes with |x - xc| ≤ hx and |y - yc| ≤ hy. With `antisymmetric`,
w = +1 above and -1 below the box centre (a symmetric box centred on the zone would cancel,
since Vx is antisymmetric across it). Vx column 1, the periodic image of the last column,
is left out. W[i, j + 1] ↔ Vx node (grid.xvi[1][i], grid.xci[2][j]) (one ghost row in y).
"""
function box_observation(Vx, grid, center, half_width; antisymmetric = true)
    W = zeros(Float64, size(Vx))
    xg, yg = collect(grid.xvi[1]), collect(grid.xci[2])
    (xc, yc), (hx, hy) = center, half_width
    for j in eachindex(yg), i in 2:length(xg)
        (abs(xg[i] - xc) <= hx && abs(yg[j] - yc) <= hy) || continue
        W[i, j + 1] = antisymmetric ? sign(yg[j] - yc) : 1.0
    end
    count(!iszero, W) == 0 && throw(ArgumentError("the observation box contains no Vx nodes"))
    return (; field = :Vx, weights = PTArray(backend)(W))
end

objective_value(stokes, observation) = sum(Array(observation.weights) .* Array(stokes.V.Vx))
