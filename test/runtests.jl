using SymPyPythonCall

path = joinpath(pathof(SymPyPythonCall.SymPyCore), "../../test")
include(joinpath(path, "runtests-sympycore.jl"))

## VERSION >= v"1.9.0" && include("test-extensions.jl")
