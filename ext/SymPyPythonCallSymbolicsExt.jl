module SymPyPythonCallSymbolicsExt

# from https://github.com/JuliaSymbolics/Symbolics.jl/pull/957/
# by @jClugstor
import SymPyPythonCall
sp = SymPyPythonCall._sympy_
const PythonCall = SymPyPythonCall.PythonCall
import PythonCall: pyconvert, pyimport, pyisinstance

import Symbolics
import Symbolics: @variables

PythonCall.pyconvert(::Type{T}, x::SymPyPythonCall.Sym) where {T} = pyconvert(T, x.o)

# rule functions
function pyconvert_rule_sympy_symbol(::Type{Symbolics.Num}, x)
    if !pyisinstance(x,sp.Symbol)
        return PythonCall.pyconvert_unconverted()
    end
    name = PythonCall.pyconvert(Symbol,x.name)
    return PythonCall.pyconvert_return(Symbolics.variable(name))
end

function pyconvert_rule_sympy_pow(::Type{Symbolics.Num}, x)
    if !pyisinstance(x,sp.Pow)
        return PythonCall.pyconvert_unconverted()
    end
    expbase = pyconvert(Symbolics.Num,x.base)
    exp = pyconvert(Symbolics.Num,x.exp)
    return PythonCall.pyconvert_return(expbase^exp)
end

function pyconvert_rule_sympy_mul(::Type{Symbolics.Num}, x)
    if !pyisinstance(x,sp.Mul)
        return PythonCall.pyconvert_unconverted()
    end
    mult = reduce(*,PythonCall.pyconvert.(Symbolics.Num,x.args))
    return PythonCall.pyconvert_return(mult)
end

function pyconvert_rule_sympy_add(::Type{Symbolics.Num}, x)
    if !pyisinstance(x,sp.Add)
        return PythonCall.pyconvert_unconverted()
    end
    sum = reduce(+, PythonCall.pyconvert.(Symbolics.Num,x.args))
    return PythonCall.pyconvert_return(sum)
end

function pyconvert_rule_sympy_derivative(::Type{Symbolics.Num}, x)
    if !pyisinstance(x,sp.Derivative)
        return PythonCall.pyconvert_unconverted()
    end
    variables = pyconvert.(Symbolics.Num,x.variables)
    derivatives = prod(var -> Differential(var), variables)
    expr = pyconvert(Symbolics.Num, x.expr)
    return PythonCall.pyconvert_return(derivatives(expr))
end

function pyconvert_rule_sympy_function(::Type{Symbolics.Num}, x)
    if !pyisinstance(x,sp.Function)
        return PythonCall.pyconvert_unconverted()
    end
    nm = PythonCall.pygetattr(x, "func", nothing)
    isnothing(nm) && return PythonCall.pyconvert_unconverted() # XXX
    name = pyconvert(Symbol, nm)
    args = pyconvert.(Symbolics.Num, x.args)
    func = @variables $name(..)
    return PythonCall.pyconvert_return(first(func)(args...))
end

function pyconvert_rule_sympy_equality(::Type{Symbolics.Equation}, x)
    if !pyisinstance(x,sp.Equality)
         return PythonCall.pyconvert_unconverted()
    end
    rhs = pyconvert(Symbolics.Num,x.rhs)
    lhs = pyconvert(Symbolics.Num,x.lhs)
    return PythonCall.pyconvert_return(rhs ~ lhs)
end


function __init__()
    # added rules
    # T = Symbolics.Num
    PythonCall.pyconvert_add_rule("sympy.core.symbol:Symbol", Symbolics.Num, pyconvert_rule_sympy_symbol)

    PythonCall.pyconvert_add_rule("sympy.core.power:Pow",     Symbolics.Num, pyconvert_rule_sympy_pow)

    PythonCall.pyconvert_add_rule("sympy.core.mul:Mul",       Symbolics.Num, pyconvert_rule_sympy_mul)

    PythonCall.pyconvert_add_rule("sympy.core.add:Add",       Symbolics.Num, pyconvert_rule_sympy_add)

    PythonCall.pyconvert_add_rule("sympy.core.function:Derivative", Symbolics.Num, pyconvert_rule_sympy_derivative)

    PythonCall.pyconvert_add_rule("sympy.core.function:Function",   Symbolics.Num, pyconvert_rule_sympy_function)

    # T = Symbolics.Equation
    PythonCall.pyconvert_add_rule("sympy.core.relational:Equality", Symbolics.Equation, pyconvert_rule_sympy_equality)

    # core numbers
    add_pyconvert_rule(f, cls, T=Symbolics.Num) = PythonCall.pyconvert_add_rule(cls, T, f)

    add_pyconvert_rule("sympy.core.numbers:Pi") do T::Type{Symbolics.Num}, x
        PythonCall.pyconvert_return(Symbolics.Num(pi))
    end
    add_pyconvert_rule("sympy.core.numbers:Exp1") do T::Type{Symbolics.Num}, x
        PythonCall.pyconvert_return(Symbolics.Num(ℯ))
    end
    add_pyconvert_rule("sympy.core.numbers:Infinity") do T::Type{Symbolics.Num}, x
        PythonCall.pyconvert_return(Symbolics.Num(Inf))
    end
    # Complex{Num}
    add_pyconvert_rule("sympy.core.numbers:ImaginaryUnit", Complex{Symbolics.Num}) do T::Type{Complex{Symbolics.Num}}, x
        PythonCall.pyconvert_return(Complex(Symbolics.Num(0), Symbolics.Num{1}))
    end
    add_pyconvert_rule("sympy.core.numbers:ComplexInfinity", Complex{Symbolics.Num}) do T::Type{Complex{Symbolics.Num}}, x
        PythonCall.pyconvert_return(Complex(Symbolics.Num(0), Symbolics.Num(Inf)) )
    end
end

end
