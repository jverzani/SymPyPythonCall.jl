using SymPyPythonCall

path = joinpath(pathof(SymPyPythonCall.SymPyCore), "../../test")
include(joinpath(path, "runtests-sympycore.jl"))
