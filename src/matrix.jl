# matrix
# SymMatrix has methods to expose.
# XXX has issue with @inferred
# XXX link into pyconvert...

Base.eachrow(M::Matrix{T}) where {T <: SymbolicObject} = (M[i,:] for i ∈ 1:size(M,1))
function Base.view(A::AbstractArray{T,N}, I::Vararg{Any,M}) where {T <: SymbolicObject, N,M}
    A[I...] # can't take view!!
end



Base.promote_op(::T, ::Type{S}, ::Type{Sym}) where {T, S <: Number} = Sym
Base.promote_op(::T, ::Type{Sym}, ::Type{S}) where {T, S <: Number} = Sym
Base.promote_op(::T, ::Type{Sym}, ::Type{Sym}) where {T} = Sym # This helps out linear alg


function unSym(M::AbstractMatrix{T}) where {T <: Sym}
    sympy.py.Matrix(Tuple(Mᵢ for Mᵢ ∈ eachrow(M)))
end

call_matrix_meth(M::Matrix{T}, meth::Symbol, args...; kwargs...) where {T <: Sym} =
    ↑(getproperty(↓(M), meth)(↓(args)...; unSymkwargs(kwargs)...))

# include those exported from base or LinearAlgebra with corresponding sympy
matrix_meths = (
    (:adjoint, :Base, :adjoint),
    (:exp, :Base, :exp),
    (:inv, :Base, :inv),
    (:transpose, :Base, :transpose),
    (:det, :LinearAlgebra, :det),
    (:diag, :LinearAlgebra, :diag),
    (:eigenvects, :LinearAlgebra, :eigvecs),
    (:eigenvals, :LinearAlgebra, :eigvals),
    (:norm, :LinearAlgebra, :norm),
    (:pinv, :LinearAlgebra, :pinv),
    (:rank, :LinearAlgebra, :rank),
)

for (smeth, jmod, jmeth) ∈ matrix_meths
    @eval begin
        ($(jmod).$(jmeth))(M::Matrix{T}, args...; kwargs...) where {T <: Sym} =
            ↑(↓(M).$(smeth)(↓(args)...; unSymkwargs(kwargs)...))
    end
end


function Base.getproperty(M::AbstractArray{<:Sym}, prop::Symbol)
    MM = sympy.Matrix.val(↓(M))
    val = pygetattr(MM, string(prop))
    #val = pygetattr(↓(M), string(prop))
    meth = pygetattr(val, "__call__", nothing)
    meth == nothing && return Sym(val)
    return SymbolicCallable(meth)
end

# function Base.getproperty(X::Vector{Sym}, prop::Symbol)
#     val = pygetattr(sympy.Matrix.val(Tuple(xᵢ for xᵢ ∈ X)), string(prop))
#     meth = pygetattr(val, "__call__", nothing)
#     meth == nothing && return Sym(val)
#     return SymbolicCallable(meth)
# end

LinearAlgebra.qr(A::AbstractArray{Sym,2}) = ↑(↓(A).QRdecomposition())

# solve Ax=b for x, avoiding generic `lu`, which can be very slow for bigger n values
# fix suggested by @olof3 in issue 355
function LinearAlgebra.:\(A::AbstractArray{Sym,2}, b::AbstractArray{S,1}) where {S}

    m,n  = size(A)
    x =  [Sym("x$i") for  i in 1:n]
    out = solve(A*x-b, x)
    isempty(out) && throw(SingularException(0)) # Could also return out here?
    ret = Vector{Sym}(undef, n)
    for (i,xᵢ)  in enumerate(x)
        ret[i] =  get(out,  xᵢ, xᵢ)
    end

    return ret

end

# function LinearAlgebra.:\(A::AbstractArray{T,2}, B::AbstractArray{S,2}) where {T <: Sym, S}
#     hcat([A \ bⱼ for bⱼ in eachcol(B)]...)
# end



LinearAlgebra.norm(x::AbstractVector{T}, args...; kwargs...) where {T <: SymbolicObject} =
    ↑(getproperty(sympy.Matrix(Tuple(xᵢ for xᵢ ∈ x)), :norm)(↓(args)...; unSymkwargs(kwargs)...))



## Issue #359 so that A  +  λI is of type Sym
Base.:+(A::AbstractMatrix{T}, J::UniformScaling) where {T <: SymbolicObject}    = _sym_plus_I(A,J)
Base.:+(A::AbstractMatrix, J::UniformScaling{T}) where {T <: SymbolicObject}    = _sym_plus_I(A,J)
Base.:+(A::AbstractMatrix{T}, J::UniformScaling{T}) where {T <: SymbolicObject} = _sym_plus_I(A,J)

Base.:-(J::UniformScaling, A::AbstractMatrix{T}) where {T <: SymbolicObject}    = (-A) + J
Base.:-(J::UniformScaling{T}, A::AbstractMatrix) where {T <: SymbolicObject}    = (-A) + J
Base.:-(J::UniformScaling{T}, A::AbstractMatrix{T}) where {T <: SymbolicObject} = (-A) + J

function _sym_plus_I(A::AbstractArray{T,N}, J::UniformScaling{S}) where {T, N, S}
    n = LinearAlgebra.checksquare(A)
    B = convert(AbstractArray{promote_type(T,S),N}, copy(A))
    for i ∈ 1:n
        B[i,i] += J.λ
    end
    B
end
