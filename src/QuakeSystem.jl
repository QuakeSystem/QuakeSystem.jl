"""
QuakeSystem module

see https://github.com/QuakeSystem/QuakeSystem.jl
"""

module QuakeSystem

include("sens/Sensitivities.jl")
include("optimisation/Optimization.jl")
include("da/DataAssimilation.jl")
include("vis/Plotting.jl")

end