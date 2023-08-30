# we have this picture for calling SymPy methods
# Julia ⋯⋯⋯⋯ >  Julia
#  |                ^
#  |                |
#  | (unSym, Py, ↓) | (asSymbolic, ↑)
#  |                |
#  v      (SymPy)   |
# Python  ------> Python


"""
    asSymbolic(x)

Convert `Python` object into symbolic `Julia` object. Aliased to `↑`.
"""
asSymbolic(x) = Sym(pyconvert(Py, x))
asSymbolic(x::Sym) = x
asSymbolic(x::String) = Sym(x)

PyTypeName(x) = Val(Symbol(PythonCall.pyconvert_typename(pytype(x))))
asSymbolic(x::Py) = asSymbolic(PyTypeName(x), x) # special case based on typename

asSymbolic(::Val{T}, x::Py) where {T} = Sym(x) # fallback

asSymbolic(::Val{Symbol("builtins:list")}, x::Py) = isempty(x) ? Sym[] : [asSymbolic(xᵢ) for xᵢ ∈ x] # XXX are lists Tuples or vectors???
asSymbolic(::Val{Symbol("builtins:tuple")}, x::Py) = Tuple(asSymbolic(xᵢ) for xᵢ ∈ x)
asSymbolic(::Val{Symbol("builtins:dict")}, x::Py) = Dict(asSymbolic(k) => asSymbolic(v) for (k,v) ∈ pyconvert(PyDict, x))
asSymbolic(::Val{Symbol("builtins:set")}, x::Py) = Set(asSymbolic(xᵢ) for xᵢ ∈ x)
asSymbolic(::Val{Symbol("sympy.sets.sets:FiniteSet")}, x::Py) = Set(asSymbolic(xᵢ) for xᵢ ∈ x)
asSymbolic(::Val{Symbol("sympy.core.containers:Tuple")}, x::Py) = Tuple(asSymbolic(xᵢ) for xᵢ ∈ x)

function _as_symbolic_dense_matrix(x)
    sz = pyconvert(Tuple, x.shape)
    return [asSymbolic(x.__getitem__((i-1, j-1))) for i ∈ 1:sz[1], j ∈ 1:sz[2]]
end
asSymbolic(::Val{Symbol("sympy.matrices.dense:MutableDenseMatrix")}, x::Py) = _as_symbolic_dense_matrix(x)
asSymbolic(::Val{Symbol("sympy.matrices.dense:ImmutableDenseMatrix")}, x::Py) = _as_symbolic_dense_matrix(x)

function _as_symbolic_symbolic_matrix(x)
    sz = pyconvert(Tuple, x.shape)
    !isa(sz[1], Real) && return Sym(x) # MatrixSymbol special case
    return [asSymbolic(x.__getitem__((i-1, j-1))) for i ∈ 1:sz[1], j ∈ 1:sz[2]]
end
asSymbolic(::Val{Symbol("sympy.matrices.expressions.matexpr:MatrixSymbol")}, x::Py) = _as_symbolic_symbolic_matrix(x)
asSymbolic(::Val{Symbol("sympy.matrices.expressions.inverse:Inverse")}, x::Py) = _as_symbolic_symbolic_matrix(x)

# used to pass arguments down into calls
unSym(x) = x
unSym(x::Sym) = getfield(x, :py)
unSym(x::Tuple) = unSym.(x)
unSym(x::Vector) = Tuple(unSym.(x)) # <-- is this an issue?
unSym(x::Matrix) = sympy.py.Matrix(Tuple(↓(mᵢ) for mᵢ ∈ eachrow(x))) # mutable dense matrix
unSym(x::Dict) = convert(PyDict,Dict(unSym(k) => unSym(v) for (k,v) ∈ x)) # AbstractDict?
unSym(x::Set) = sympy.py.Set(unSym(collect(x)))
unSym(x::Irrational{:π}) = unSym(sympy.pi)
unSymkwargs(kw) = (k=>unSym(v) for (k,v) ∈ kw)


#PythonCall.Py(x::Sym) = ↓(x)
PythonCall.Py(x::NTuple{N,T}) where {N, T <: SymbolicObject} = ↓(x)
PythonCall.Py(x::Vector{T}) where {T <: SymbolicObject} = unSym.(x)
PythonCall.Py(x::Matrix{T}) where {T <: SymbolicObject} = unSym.(x)
PythonCall.Py(x::Set{T}) where {T <: SymbolicObject} = unSym(x)

# use ↑, ↓ shortcuts
↑(args...; kwargs...) = asSymbolic(args...; kwargs...)
↓(args...; kwargs...) = unSym(args...; kwargs...)


## --------------------------------------------------

# generic method pattern: (pymodule, pymethod, juliamodule, juliamethod)
# these don't need to be exported
# in SymPy we use getmembers to generate this; not sure what is better
#const
generic_methods = (
    (:sympy, :cos, :Base, :cos),
    (:sympy, :sin, :Base, :sin),
    (:sympy, :tan, :Base, :tan),
    (:sympy, :sec, :Base, :sec),
    (:sympy, :csc, :Base, :csc),
    (:sympy, :cot, :Base, :cot),
    (:sympy, :acos, :Base, :acos),
    (:sympy, :asin, :Base, :asin),
#    (:sympy, :atan, :Base, :atan),
    (:sympy, :asec, :Base, :asec),
    (:sympy, :acsc, :Base, :acsc),
    (:sympy, :acot, :Base, :acot),
    #
    (:sympy, :cosh, :Base, :cosh),
    (:sympy, :sinh, :Base, :sinh),
    (:sympy, :tanh, :Base, :tanh),
    (:sympy, :sech, :Base, :sech),
    (:sympy, :csch, :Base, :csch),
    (:sympy, :coth, :Base, :coth),
    (:sympy, :acosh, :Base, :acosh),
    (:sympy, :asinh, :Base, :asinh),
    (:sympy, :atanh, :Base, :atanh),
#    (:sympy, :asech, :Base, :asech),
#    (:sympy, :acsch, :Base, :acsch),
    (:sympy, :acoth, :Base, :acoth),
    #
    (:sympy, :sqrt, :Base, :sqrt),
    (:sympy, :exp, :Base, :exp),
    (:sympy, :log, :Base, :log),
    #
    (:sympy, :Mod, :Base, :mod),
    (:sympy, :floor, :Base, :floor),
    (:sympy, :ceiling, :Base, :ceil),
    #
    (:sympy, :numer, :Base, :numerator),
    (:sympy, :denom, :Base, :denominator),
    (:sympy, :Max, :Base, :max),
    (:sympy, :Abs, :Base, :abs),
    (:sympy, :Min, :Base, :min),
    #
    (:sympy, :re, :Base, :real),
    (:sympy, :im, :Base, :imag),

    # solve
    (:sympy, :solve, :CommonSolve, :solve),

    # Eq
    (:sympy, :Eq, :CommonEq, :Eq),
    (:sympy, :Lt, :CommonEq, :Lt),
    (:sympy, :Le, :CommonEq, :Le),
    (:sympy, :Ne, :CommonEq, :Ne),
    (:sympy, :Ge, :CommonEq, :Ge),
    (:sympy, :Gt, :CommonEq, :Gt),

    # diff
    (:sympy, :diff, :CommonEq, :diff),

    # collect
    (:sympy, :collect, :Base, :collect),

    # SpecialFunctions
    (:sympy, :airyai ,      :SpecialFunctions, :airyai),
    (:sympy, :airyaiprime , :SpecialFunctions, :airyaiprime),
    (:sympy, :airybi ,      :SpecialFunctions, :airybi),
    (:sympy, :besseli ,     :SpecialFunctions, :besseli),
    (:sympy, :besselj ,     :SpecialFunctions, :besselj),
    (:sympy, :besselk ,     :SpecialFunctions, :besselk),
    (:sympy, :bessely ,     :SpecialFunctions, :bessely),
    (:sympy, :beta ,        :SpecialFunctions, :beta),
    (:sympy, :erf ,         :SpecialFunctions, :erf),
    (:sympy, :erfc ,        :SpecialFunctions, :erfc),
    (:sympy, :erfi ,        :SpecialFunctions, :erfi),
    (:sympy, :erfinv ,      :SpecialFunctions, :erfinv),
    (:sympy, :erfcinv ,     :SpecialFunctions, :erfcinv),
    (:sympy, :gamma ,       :SpecialFunctions, :gamma),
    (:sympy, :digamma ,     :SpecialFunctions, :digamma),
    (:sympy, :polygamma ,   :SpecialFunctions, :polygamma),
    (:sympy, :hankel1,      :SpecialFunctions, :hankelh1),
    (:sympy, :hankel2,      :SpecialFunctions, :hankelh2),
    (:sympy, :zeta ,        :SpecialFunctions, :zeta),
)

# comment, Py(...) speeds things up a bit.
for (pmod, pmeth, jmod, jmeth) ∈ generic_methods
    @eval begin
        ($(jmod).$(jmeth))(x::Sym, args...; kwargs...) =
            ↑(Py($(pmod)).$(pmeth)(↓(x), ↓(args)...; unSymkwargs(kwargs)...))
    end
end

## --------------------------------------------------
# pmod, pmeth, meth
#const
new_exported_methods = (
    (:sympy, :simplify,    :simplify),
    (:sympy, :expand_trig, :expand_trig),
    (:sympy, :expand,      :expand),
    (:sympy, :together,    :together),
    (:sympy, :apart,       :apart),
    (:sympy, :factor,      :factor),
    (:sympy, :cancel,      :cancel),
    #
    (:sympy, :degree,      :degree),
    #
    (:sympy, :integrate,   :integrate),
    #
    (:sympy, :real_roots,  :real_roots),
    (:sympy, :roots,       :roots),
    (:sympy, :nroots,      :nroots),
    (:sympy, :dsolve,      :dsolve),
    (:sympy, :nsolve,      :nsolve),
    (:sympy, :linsolve,    :linsolve),
    (:sympy, :nonlinsolve, :nonlinsolve),
    (:sympy, :solveset,    :solveset),
    #
    (:sympy, :series,      :series),
    (:sympy, :summation,   :summation),
    (:sympy, :hessian,     :hessian),
)

for (pmod, pmeth, jmeth) ∈ new_exported_methods
    @eval begin
        ($(jmeth))(x::Sym, args...; kwargs...) =
            ↑(Py($(pmod)).$(pmeth)(↓(x), ↓(args)...; unSymkwargs(kwargs)...))
        export $(jmeth)
    end
end

## --------------------------------------------------

object_methods = (
    (:conjugate, :Base, :conj),
)

for (ometh, jmod, jmeth) ∈ object_methods
    @eval begin
        ($(jmod).$(jmeth))(x::Sym, args...; kwargs...) =
            ↑(Py(x).$(ometh)(↓(args)...; unSymkwargs(kwargs)...))
    end
end
