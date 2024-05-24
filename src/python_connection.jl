## PythonCall specific usage
Base.convert(::Type{S}, x::Sym{T}) where {T <: PythonCall.Py, S<:Sym} = x
Base.convert(::Type{S}, x::Sym{T}) where {T<:PythonCall.Py S<:Sym{PythonCall.Py}} = x
Base.convert(::Type{S}, x::T) where {T<:PythonCall.Py, S <: SymbolicObject} = Sym(x)

SymPyCore._convert(::Type{T}, x) where {T} = pyconvert(T, x)
function SymPyCore._convert(::Type{Bool}, x::Py)
    pyisinstance(x, _sympy_.logic.boolalg.BooleanTrue) && return true
    pyisinstance(x, _sympy_.logic.boolalg.BooleanFalse) && return false
    pyconvert(Bool, pybool(x))
end


SymPyCore.Bool3(::Sym{Nothing}) = nothing

function SymPyCore.Bool3(x::Sym{T}) where {T <: PythonCall.Py}
    y = ↓(x)
    isnothing(y) && return nothing
    if hasproperty(y, "is_Boolean")
        if pyconvert(Bool, y.is_Boolean)
            return SymPyCore._convert(Bool, y)
        end
    elseif hasproperty(y, "__bool__")
        if pyconvert(Bool, y.__bool__ != ↓(Sym(nothing)))
            return pyconvert(Bool, y.__bool__())
        end
    end
    return nothing
end

## Modifications for ↓, ↑
Sym(x::Nothing) = Sym(pybuiltins.None)

SymPyCore.:↓(x::PythonCall.Py) = x
SymPyCore.:↓(d::Dict) = pydict((↓(k) => ↓(v) for (k,v) ∈ pairs(d)))
SymPyCore.:↓(x::Set) = _sympy_.sympify(pyset(↓(sᵢ) for sᵢ ∈ x))

SymPyCore.:↑(::Type{<:AbstractString}, x) = Sym(Py(x))
function SymPyCore.:↑(::Type{PythonCall.Py}, x)
    # this lower level approach shouldn't allocate
    pyisinstance(x, pybuiltins.set) && return Set(Sym.(collect(x))) #Set(↑(xᵢ) for xᵢ ∈ x)
    pyisinstance(x, pybuiltins.tuple) && return Tuple(↑(xᵢ) for xᵢ ∈ x)
    pyisinstance(x, pybuiltins.list) && return [↑(xᵢ) for xᵢ ∈ x]
    pyisinstance(x, pybuiltins.dict) && return Dict(↑(k) => ↑(x[k]) for k ∈ x)

    # add more sympy containers in sympy.jl and here
    pyisinstance(x, _FiniteSet_) && return Set(Sym.(collect(x)))
    pyisinstance(x, _MutableDenseMatrix_) && return _up_matrix(x) #map(↑, x.tolist())

    # fallback
    Sym(x)
end


function _up_matrix(m) # ↑ for matrices
    sh = m.shape
    r, c = SymPyCore._convert(Int, sh.__getitem__(0)), SymPyCore._convert(Int, sh.__getitem__(1))
    out = [↑(m.__getitem__(i*c + j)) for i ∈ 0:(r-1), j ∈ 0:(c-1)]
    out
end


function Base.hash(x::SymbolicObject{T}, h::UInt) where {T <: PythonCall.Py}
    o = ↓(x)
    hash(o, h)
end


# should we also have different code path for a::String like  PyCall?
function Base.getproperty(x::SymbolicObject{T}, a::Symbol) where {T <: PythonCall.Py}

    a == :o && return getfield(x, a)

    if a == :py
        Base.depwarn("The field `.py` has been renamed `.o`", :getproperty)
        return getfield(x,:o)
    end

    val = ↓(x)
    𝑎 = string(a)

    hasproperty(val, 𝑎) || return nothing # not a property
    meth = getproperty(val, 𝑎)

    pyis(meth, pybuiltins.None) && return nothing

    ## __call__
    if pycallable(meth) # "__call__") #hasproperty(meth, "__call__")
        return SymPyCore.SymbolicCallable(meth)
    end

    # __class__ dispatch
    if pyisinstance(meth, _bool_)
        return pyconvert(Bool, meth)
    end

    if pyisinstance(meth, _ModuleType_)
        return Sym(meth)
    end


    # just convert
    return ↑(convert(PythonCall.Py, meth))

end


# do we need this conversion?
#Base.convert(::Type{T}, o::Py) where {T <: Sym} = T(o)
