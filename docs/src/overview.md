# Basic overview

An interface between `Julia` and SymPy requires a few things. `PythonCall` can call into an underlying sympy library, for example

```julia
julia> import PythonCall

julia> const _sympy_ = PythonCall.pynew()
Python: NULL

julia> PythonCall.pycopy!(_sympy_, PythonCall.pyimport("sympy"));
```

The `_sympy_` object holds references to the underlying sympy library. For example, we can create a symbolic variable and call the `sin` function on it as follows:

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

`SymPyCall` wraps the `Py` objects in its `Sym` class to provide a means to dispatch familiar `Julia` generics:

```jldoctest overview
julia> using SymPyCall

```

```jldoctest overview
julia> x = symbols("x")
x

julia> simplify(2sin(x)*cos(x))
sin(2*x)
```

The package also provides methods for some sympy methods, such as `simplify` above. To make this work, there needs to be a means to take `Sym` objects to the `Py` counterparts and to take `Py` objects to `Sym`. As these conversions may be type dependent two operators (with the unicode equivalents `↓` and  `↑`) are used internally to allow the definition along these lines:

```
simplify(x::Sym, args...; kwargs...) = ↑(sympy.simplify(↓(x), ↓(args)...; ↓(kwargs)...))
```

The `expand_log` function is not wrapped as such, but can be called from the `sympy` object exported by `SymPyCall`:

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

We follow part of the `SymPy` docs to see how to access one of the numerous external modules of `sympy` beyond those exposed immediately by `SymPy`.



```jldoctest overview
julia> import SymPyCall.PythonCall: pyimport

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

The one thing to note is the method calls return `Py` objects, as there is no intercepting of the method calls done the way there is for the `sympy` module. The `Sym` command can be used to wrap the output:

```jldoctest overview
julia> stats.variance(X) |> Sym
σ^2
```

The methods in the `stats` module are qualified with the module name above.

Next we see that statements like $P(X > \mu)$ can be answered specifying the inequality using `sympy.Gt` in the following:

```jldoctest overview
julia> stats.P(sympy.Gt(X, μ))
Python: 1/2
```

A typical calculation for the normal distribution is the area one or more standard deviations larger than the mean:

```jldoctest overview
julia> stats.P(sympy.Gt(X, μ + 1σ)) |> Sym
sqrt(2)*(-sqrt(2)*pi*exp(1/2)*erf(sqrt(2)/2)/2 + sqrt(2)*pi*exp(1/2)/2)*exp(-1/2)/(2*pi)
```

The familiar  answer could be found by calling `N` or `evalf`.

We show one more distribution, the uniform over $[a,b]$:

```jldoctest overview
julia> @syms a::real b::real
(a, b)

julia> U = stats.Uniform("U", a, b)
Python: U

julia> stats.E(U) |> simplify |> println
a/2 + b/2

julia> stats.variance(U) |> Sym |> simplify |> factor |> println
(a - b)^2/12
```

----

Not all outputs are so simple to incorporate. `PythonCall` does a good job of converting the arguments from `Julia` to Python, but the conversion from a Python (SymPy) structure back to a workable `Julia` structure can require some massaging.


The output of many integration problems is a piecewise function:

```jldoctest overview
julia> @syms n::integer x::real
(n, x)

julia> u = integrate(x^n, x)
Piecewise((x^(n + 1)/(n + 1), Ne(n, -1)), (log(x), True))
```


The `u` object is of type `Sym`, but there are no methods for working with it. The `.args` call will break this into its argument, which again will by symbolic. We call the `.py` attribute of the `Sym` object to get an iterable:




```jldoctest overview
julia> [Sym(t[1]) => Sym(t[0]) for t ∈ u.args.py]
2-element Vector{Pair{Sym, Sym}}:
 Ne(n, -1) => x^(n + 1)/(n + 1)
      True => log(x)
```

The python objects are iterated over, so 0-based indexing is used above. These pieces are then converted to `Sym` objects for familiarity.
