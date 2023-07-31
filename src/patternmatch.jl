## Pattern matching modifications


"""
    Wild(x)

create a "wild card" for pattern matching
"""
Wild(x::AbstractString) = sympy.Wild(x)
Wild(x::Symbol) = Wild(string(x))
export Wild


"""
    match(pattern, expression, ...)

Match a pattern against an expression; returns a dictionary of matches.

If a match is unsuccesful, returns an *empty* dictionary. (SymPy returns "nothing")

The order of the arguments follows `Julia`'s `match` function, not `sympy.match`, which can be used directly, otherwise.
"""
function Base.match(pat::Sym, ex::Sym, args...; kwargs...)
    out = ex.match(pat, args...; kwargs...)
    out == nothing && return Dict()
    out
end

"""
    replace(expression, pattern, value, ...)

From: [SymPy Docs](http://docs.sympy.org/dev/modules/core.html)

Traverses an expression tree and performs replacement of matching
subexpressions from the bottom to the top of the tree. The default
approach is to do the replacement in a simultaneous fashion so changes
made are targeted only once. If this is not desired or causes
problems, `simultaneous` can be set to `false`. In addition, if an
expression containing more than one `Wild` symbol is being used to match
subexpressions and the `exact` flag is `true`, then the match will only
succeed if non-zero values are received for each `Wild` that appears in
the match pattern.


Differences from SymPy:

* "types" are specified via calling `func` on the head of an expression: `func(sin(x))` -> `sin`, or directly through `sympy.sin`

* functions are not supported


Examples (from the SymPy docs)

```jldoctest replace
julia> using SymPyCall

julia> @syms x, y, z
(x, y, z)

julia> f = log(sin(x)) + tan(sin(x^2)); print(f) # `print(f)` only so doctest can run
log(sin(x)) + tan(sin(x^2))

```

## "type" -> "type"

Types are specified through `func`:

```jldoctest replace
julia> func = SymPyCall.Introspection.func
func (generic function with 1 method)

julia> replace(f, func(sin(x)), func(cos(x))) |> print  # type -> type
log(cos(x)) + tan(cos(x^2))

julia> #replace(f, sympy.sin, sympy.cos)  |>  print  # "log(cos(x)) + tan(cos(x^2))" **but fails**

julia> sin(x).replace(func(sin(x)), func(cos(x)), map=true)
(cos(x), Dict{Sym, Sym}(sin(x) => cos(x)))

```


## "type" -> "function"

XXX

## "pattern" -> "expression"

Using "`Wild`" variables allows a pattern to be replaced by an expression:

```jldoctest replace
julia> a, b = Wild("a"), Wild("b")
(a_, b_)

julia> replace(f, sin(a), tan(2a)) |> print
log(tan(2*x)) + tan(tan(2*x^2))

julia> replace(f, sin(a), tan(a/2)) |> print
log(tan(x/2)) + tan(tan(x^2/2))

julia> f.replace(sin(a), a) |> print
log(x) + tan(x^2)

julia> (x*y).replace(a*x, a)
y

```

In the SymPy docs we have:

Matching is exact by default when more than one Wild symbol is used: matching fails unless the match gives non-zero values for all Wild symbols."

```jldoctest replace
julia> replace(2x + y, a*x+b, b-a)  # y - 2
y - 2

julia> replace(2x + y, a*x+b, b-a, exact=false)  # y + 2/x
y + 2/x
```

## "pattern" -> "func"

XXX

## "func" -> "func"

XXX

"""
function Base.replace(ex::Sym, query::Sym, fn::Function; exact=true, kwargs...)
    ## XXX this is failing!
    ex.replace(query, PyCall.PyObject((args...) ->fn(args...)); exact=exact, kwargs...)
end


function Base.replace(ex::Sym, query::Any, value; exact=true, kwargs...)
    ex.replace(query, value; exact=exact, kwargs...)
end
