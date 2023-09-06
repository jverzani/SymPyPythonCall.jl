# Basic overview

An interface between `Julia` and SymPy requires a connection between the two languages. The `PythonCall` package provides a means to  call into an underlying sympy library in `Python`. For example

```julia
julia> import PythonCall

julia> const _sympy_ = PythonCall.pynew()
Python: NULL

julia> PythonCall.pycopy!(_sympy_, PythonCall.pyimport("sympy"));
```

The `_sympy_` object holds references to the underlying sympy library. As an example, the following creates a symbolic variable, `x`, and calls the `sin` function on it:

```julia
julia> x = _sympy_.symbols("x")
Python: x

julia> _sympy_.sin(x)
Python: sin(x)

julia> x.is_commutative
Python: True

julia> typeof(x)
PythonCall.Py
```

The `PythonCall` package provides some basic operations for `Py` objects, such as basic math operations:

```julia
julia> x + x
Python: 2*x

julia> x * x
Python: x**2

julia> x ^ x
Python: x**x
```

`SymPyPythonCall` wraps the `Py` objects in its `Sym` class to provide a means to dispatch familiar `Julia` generics:

```jldoctest overview
julia> using SymPyPythonCall

```

```jldoctest overview
julia> x = symbols("x") # or @syms x
x

julia> simplify(2sin(x)*cos(x))
sin(2*x)
```

The package also provides methods for some sympy methods, such as `simplify` above. To make this work, there needs to be a means to take `Sym` objects to their `Py` counterparts and a means to take `Py` objects to a symbolic type. As these conversions may be type dependent two operators (with the unicode equivalents `↓` and  `↑`) are used internally to allow the definition along these lines:

```
simplify(x::Sym, args...; kwargs...) = ↑(sympy.simplify(↓(x), ↓(args)...; ↓(kwargs)...))
```

(Though there is some overhead introduced, it does not seem to be significant compared to computational cost of most symbolic computations.)

The `expand_log` function is not wrapped as such, but can still be called from the `sympy` object exported by `SymPyPythonCall`:

```jldoctest overview
julia> @syms x::real
(x,)

julia> simplify(log(2x))
log(2*x)

julia> sympy.expand_log(log(2x))
log(x) + log(2)
```

Methods of `sympy` are also called using the conversion operators above.

## Using other SymPy modules

We follow part of the `SymPy` docs to see how to access one of the numerous external modules of `sympy` beyond those exposed immediately by `SymPy`. In this case, the `stats` module.



```jldoctest overview
julia> import SymPyPythonCall.PythonCall: pyimport

julia> stats = pyimport("sympy.stats");
```

The `stats` module holds several probability functions, similar to the `Distributions` package of `Julia`. This set of commands creates a normally distributed random variable, `X`, with symbolic parameters:

```jldoctest overview
julia> @syms μ, σ::positive;

julia> X = stats.Normal("X", μ, σ)
Python: X

julia> stats.E(X)
Python: μ

julia> stats.E(X^2)
Python: μ**2 + σ**2

julia> stats.variance(X)
Python: σ**2
```

The one thing to note is the method calls return `Py` objects, as there is no intercepting of the method calls done the way there is for the `sympy` module. For `Sym` objects, like `sympy`, the `getproperty` method is overridden to wrap any underlying Python method in an internal `SymbolicCallable` type, obviating the need to call `Sym` as above. This can be done here, as follows:

```jldoctest overview
julia> stats = Sym(stats);

julia> stats.variance(X)
σ^2
```

Next statements like $P(X > \mu)$ can be answered by specifying the inequality using `Gt` in the following:

```jldoctest overview
julia> stats.P(Gt(X, μ))
1/2
```

The unicode `≧` operator (`\geqq[tab]`) is an infix alternative to `Gt`.

A typical calculation for the normal distribution is the area one or more standard deviations larger than the mean:

```jldoctest overview
julia> stats.P(X ≧ μ + 1 * σ)
sqrt(2)*(-sqrt(2)*pi*exp(1/2)*erf(sqrt(2)/2)/2 + sqrt(2)*pi*exp(1/2)/2)*exp(-1/2)/(2*pi)
```

The familiar  answer could be found by calling `N` or `evalf`.

One more distribution is illusrated, the uniform distribution over a symbolic interval $[a,b]$:

```jldoctest overview
julia> @syms a::real b::real
(a, b)

julia> U = stats.Uniform("U", a, b)
U

julia> stats.E(U)
a/2 + b/2

julia> stats.variance(U) |> factor
(a - b)^2/12
```

## Different output types

`SymPyPythonCall` provides a few conversions into containers of symbolic objects, like for lists, tuples, finite sets, and matrices.
Not all outputs are so simple to incorporate and are simply wrapped in the `Sym` type.
Conversion to a workable `Julia` structure can require some massaging. This example shows how to get at the pieces of the `Piecewise` type.

The output of many integration problems is a piecewise function:

```jldoctest overview
julia> @syms n::integer x::real
(n, x)

julia> u = integrate(x^n, x)
Piecewise((x^(n + 1)/(n + 1), Ne(n, -1)), (log(x), True))
```


The `u` object is of type `Sym`, but there are no methods for working with it. The `.args` call will break this into its argument, which again will by symbolic. The `SymPyPythonCall.Introspection.args` function will perform the same thing.

```jldoctest overview
julia> as = SymPyPythonCall.Introspection.args(u)
((x^(n + 1)/(n + 1), Ne(n, -1)), (log(x), True))
```

The `as` object is a tuple of `Sym` objects. Consider the first one:

```jldoctest overview
julia> c1 = first(as)
(x^(n + 1)/(n + 1), Ne(n, -1))
```

The value of `c1` prints as a tuple, but is of type `Sym` and sympy type `CondExprPair`. Though the underlying python type is iterable or indexable, the wrapped type is not. It can  be made iterable in a manner of ways: by finding the underlying Python object through the attribute `.py`, as in `c1.py`, or by calling the `Py` method of `PythonCall`, as in `PythonCall.Py(c1)`. More generically, `PythonCall` provides a `PyIterable` wrapper. As there is no defined `lastindex`, the latter is a bit more cumbersome to use. This pattern below, expands the conditions into a dictionary:

```jldoctest overview
julia> [Sym(a.py[1]) => Sym(a.py[0]) for a ∈ as]
2-element Vector{Pair{Py, Py}}:
 Ne(n, -1) => x**(n + 1)/(n + 1)
      True => log(x)
```

The python objects are iterated over, so 0-based indexing is used above. These pieces are then converted to `Sym` objects for familiarity.
