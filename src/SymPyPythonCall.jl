"""
    SymPyPythonCall

Module to call Python's SymPy package from Julia using PythonCall
"""
module SymPyPythonCall

using SymPyCore
using PythonCall

const _PyType = PythonCall.Py
_pynull() = PythonCall.pynew()
_copy!(a, b) = PythonCall.pycopy!(a,b)
_pyimport(a) = PythonCall.pyimport(a)
_pyimport_conda(a,b) = PythonCall.pyimport(a)  # XXX lose things
_pyobject(x) = PythonCall.pyconvert(Py, x)
_pytype_mapping(typ,a)  = nothing

core_src_path = joinpath(pathof(SymPyCore), "../../src/SymPy")
include(joinpath(core_src_path, "sympy.jl"))

include("python_connection.jl")


end
