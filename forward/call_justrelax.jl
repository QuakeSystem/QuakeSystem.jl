using Pkg
Pkg.activate(@__DIR__)
Pkg.add(url="https://github.com/QuakeSystem/JustRelax.jl", rev="v0.8.1-QS") # add the JustRelax package from our fork
#Pkg.add(url="https://github.com/QuakeSystem/JustRelax.jl", rev="dev_iska") # or certain branch

using JustRelax
