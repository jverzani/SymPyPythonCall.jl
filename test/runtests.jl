using SymPyPythonCall

using SymPyCore
using LinearAlgebra
using SpecialFunctions
using Test


path = joinpath(pathof(SymPyPythonCall.SymPyCore), "../../test")

include(joinpath(path, "test-legacy.jl")) # need to clean these up!
include(joinpath(path, "test-core.jl"))
include(joinpath(path, "test-math.jl"))
include(joinpath(path, "test-matrix.jl"))
include(joinpath(path, "test-specialfuncs.jl"))
include(joinpath(path, "test-ode.jl"))
#include(joinpath(path, "test-lambdify.jl"))
#include(joinpath(path, "test-logical.jl"))
#include(joinpath(path, "test-permuations.jl"))
#include(joinpath(path, "test-physics.jl"))
#include(joinpath(path, "test-external-module.jl"))
