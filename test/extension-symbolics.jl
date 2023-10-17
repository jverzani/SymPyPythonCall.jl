using SymPyPythonCall
import Symbolics
using Test

@test isa(SymPyPythonCall.PythonCall.pyconvert(Symbolics.Num, sympy.sympify("x")), Symbolics.Num)
