using SymPyPythonCall
using Test

include("tests.jl")
include("test-math.jl")
include("test-matrix.jl")
include("test-ode.jl")
#=
=#
include("test-logical.jl")
include("test-specialfuncs.jl")
#include("test-permutations.jl")
#include("test-physics.jl")
#include("test-external-module.jl")
include("test-latexify.jl")

if VERSION >= v"1.9.0-"
    @testset "Symbolics integration" begin include("symbolics-integration.jl") end
end
