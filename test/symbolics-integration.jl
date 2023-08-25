using SymPyPythonCall
using Symbolics

@test isa(SymPyPythonCall.PythonCall.pyconvert(Symbolics.Num, sympy.sympify("x")), Symbolics.Num)
