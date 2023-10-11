## PythonCall specific usage
SymPyCore._convert(::Type{T}, x) where {T} = pyconvert(T, x)
SymPyCore._convert(::Type{Bool}, x)  = pyconvert(Bool, pybool(x))


## Modifications for ↓, ↑
SymPyCore.:↓(x::PythonCall.Py) = x
SymPyCore.:↓(d::Dict) = pydict((↓(k) => ↓(v) for (k,v) ∈ pairs(d)))
SymPyCore.:↓(x::Set) = _sympy_.sympify(pyset(↓(sᵢ) for sᵢ ∈ x))
SymPyCore.:↑(::Type{<:AbstractString}, x) = Py(x)

Sym(x::Nothing) = Sym(pybuiltins.None)

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
    val = ↓(x)
    if hasproperty(val, a)
        meth = getproperty(val, a)
        pyconvert(Bool, meth == pybuiltins.None) && return nothing

        if hasproperty(meth, "is_Boolean")
            o = Sym(getproperty(meth, "is_Boolean"))
            o == Sym(true) && return true
            a == :is_Boolean && return o == Sym(False) ? false : nothing
        end

        # __class__ dispath
        if hasproperty(meth, :__class__)
            cnm = string(meth.__class__.__name__)
            if cnm == "bool"
                a = Sym(meth)
                return a == Sym(true) ? true :
                    a == Sym(false) ? false : nothing
            end
            if cnm == "module"
                # treat modules, callsm others differently
                return Sym(meth)
            end
        end
        ## __function__
        if hasproperty(meth, "__call__")
            meth = getproperty(meth, "__call__")
            return SymPyCore.SymbolicCallable(meth)
        end

        return Sym(convert(PythonCall.Py, meth))

    end
    return nothing
    @show :huh_getproperty_not_found, string(a)
end


# do we need this conversion?
#Base.convert(::Type{T}, o::Py) where {T <: Sym} = T(o)
