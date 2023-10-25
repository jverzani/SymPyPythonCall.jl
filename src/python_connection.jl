## PythonCall specific usage
Base.convert(::Type{Complex{S}}, x::Sym{T}) where {S<:Number, T<:PythonCall.Py} =
    Complex(convert(S, real(x)), convert(S, imag(x)))

Base.convert(::Type{S}, x::Sym{T}) where {T <: PythonCall.Py, S<:Sym} = x
Base.convert(::Type{S}, x::T) where {T<:PythonCall.Py, S <: SymbolicObject} = Sym(x)

SymPyCore._convert(::Type{T}, x) where {T} = pyconvert(T, x)
function SymPyCore._convert(::Type{Bool}, x::Py)
    pyconvert(Bool, x == _sympy_.logic.boolalg.BooleanTrue) && return true
    pyconvert(Bool, x == _sympy_.logic.boolalg.BooleanFalse) && return false
    pyconvert(Bool, pybool(x))
end


## Modifications for ↓, ↑
Sym(x::Nothing) = Sym(pybuiltins.None)
SymPyCore.:↓(x::PythonCall.Py) = x
SymPyCore.:↓(d::Dict) = pydict((↓(k) => ↓(v) for (k,v) ∈ pairs(d)))
SymPyCore.:↓(x::Set) = _sympy_.sympify(pyset(↓(sᵢ) for sᵢ ∈ x))

SymPyCore.:↑(::Type{<:AbstractString}, x) = Py(x)
function SymPyCore.:↑(::Type{PythonCall.Py}, x)
    class_nm = SymPyCore.classname(x)
    class_nm == "set"       && return Set(Sym.(collect(x)))
    class_nm == "tuple" && return Tuple(↑(xᵢ) for xᵢ ∈ x)
    class_nm == "list" && return [↑(xᵢ) for xᵢ ∈ x]
    class_nm == "dict" && return Dict(↑(k) => ↑(x[k]) for k ∈ x)

    class_nm == "FiniteSet" && return Set(Sym.(collect(x)))
    class_nm == "MutableDenseMatrix" && return _up_matrix(x) #map(↑, x.tolist())

    # others ... more hands on than pytype_mapping

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
    !hasproperty(val, a) && return nothing # not a property

    meth = getproperty(val, a)

    pyconvert(Bool, meth == pybuiltins.None) && return nothing

    if hasproperty(meth, "is_Boolean")
        if pyconvert(Bool, meth.is_Boolean == true)
            return meth == _sympy_.logic.boolalg.BooleanFalse
        end
    end

    # __class__ dispatch
    if hasproperty(meth, :__class__)
        cnm = string(meth.__class__.__name__)
        if cnm == "bool"
            return pyconvert(Bool, meth)
        end
        if cnm == "module"
            # treat modules, callsm others differently
            return Sym(meth)
        end
    end

    ## __call__
    if hasproperty(meth, "__call__")
        return SymPyCore.SymbolicCallable(meth)
    end

    return ↑(convert(PythonCall.Py, meth))

end


# do we need this conversion?
#Base.convert(::Type{T}, o::Py) where {T <: Sym} = T(o)
