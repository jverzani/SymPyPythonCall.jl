"""
    SymPyPythonCall

Module to call Python's SymPy package from Julia using PythonCall
"""
module SymPyPythonCall

# usings/imports
using PythonCall
import PythonCall: Py
using SpecialFunctions
using LinearAlgebra
using Markdown
import CommonSolve
import CommonSolve: solve
using CommonEq
using RecipesBase
using Latexify

# includes
include("types.jl")
include("utils.jl")
include("decl.jl")
include("constructors.jl")
include("generic.jl")
include("arithmetic.jl")
include("equations.jl")
include("gen_methods.jl")
include("conveniences.jl")
include("logical.jl")
include("matrix.jl")
include("assumptions.jl")
include("lambdify.jl")
include("introspection.jl")
include("patternmatch.jl")
include("plot_recipes.jl")
include("latexify_recipe.jl")
# export
export sympy,
    Sym,
    @syms, symbols,
    N,
    solve, Eq, Lt, Le, Ne, Ge, Gt,
    PI, E, oo, zoo, IM, TRUE, FALSE,
    𝑄, ask, refine,
    rewrite,
    Differential


# define sympy in int
const _sympy_ = PythonCall.pynew()
const sympy = Sym(_sympy_)

# core.sympy.numbers
const _PI_ = PythonCall.pynew()
const PI = Sym(_PI_)

const _E_ = PythonCall.pynew()
const E = Sym(_E_)


const _IM_ = PythonCall.pynew()
const IM = Sym(_IM_)

const _oo_ = PythonCall.pynew()
const oo = Sym(_oo_)

const _zoo_ = PythonCall.pynew()
const zoo = Sym(_zoo_)

const _TRUE_ = PythonCall.pynew()
const TRUE = Sym(_TRUE_)

const _FALSE_ = PythonCall.pynew()
const FALSE = Sym(_FALSE_)

const _𝑄_ = PythonCall.pynew()
const 𝑄 = Sym(_𝑄_)

function __init__()

    PythonCall.pycopy!(_sympy_, PythonCall.pyimport("sympy"))

    PythonCall.pycopy!(_PI_, _sympy_.pi)
    PythonCall.pycopy!(_E_, _sympy_.E)
    PythonCall.pycopy!(_IM_, _sympy_.I)
    PythonCall.pycopy!(_oo_, _sympy_.oo)
    PythonCall.pycopy!(_zoo_, _sympy_.zoo)
    PythonCall.pycopy!(_TRUE_, pyconvert(Py, true))
    PythonCall.pycopy!(_FALSE_, pyconvert(Py, false))
    PythonCall.pycopy!(_𝑄_, _sympy_.Q)


end

end
