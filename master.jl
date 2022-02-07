using Arpack, Oscar, HomotopyContinuation, LinearAlgebra
import Pkg; Pkg.add("MultiFloats")
using MultiFloats
using GenericSVD
using Plots, PolynomialAmoebas

Pkg.develop(path=raw"/Users/kemal/github/PolynomialAmoebas.jl/")
MultiFloats.use_bigfloat_transcendentals()
#Pkg.add("https://github.com/JuliaLinearAlgebra/GenericSVD.jl.git")
https://github.com/saschatimme/PolynomialAmoebas.jl

include("get_vtx.jl")
include("get_newton.jl")
