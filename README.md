# SymPyPythonCall

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://jverzani.github.io/SymPyPythonCall.jl/dev)
[![Build Status](https://github.com/jverzani/SymPyPythonCall.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jverzani/SymPyPythonCall.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jverzani/SymPyPythonCall.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jverzani/SymPyPythonCall.jl)


This package allows access to the [SymPy](https://www.sympy.org/en/index.html) Python library to `Julia` users through [PythonCall](https://github.com/cjdoris/PythonCall.jl).

(The more established [SymPy.jl](https://github.com/JuliaPy/SymPy.jl) uses [PyCall.jl](https://github.com/JuliaPy/PyCall.jl).)

At the moment, the expectation is that *if* that change proves desirable, this would become `SymPy`, but for now this is a standalone package. This may be of interest for those having difficulty installing the underlying `sympy` library using `PyCall`.

----

Though nearly the same as `SymPy.jl`, for now, there are some small design decisions differing from `SymPy`:

* `@vars` would be deprecated; use `@syms` only

* `elements` for sets is deprecated (conversion to a `Set` is the newdefault)

* `sympy.poly` *not* `sympy.Poly`

* `limit(ex, x, c)` deprecated; use `limit(ex, x=>c)` or `sympy.limit`

* `Base.show` isn't *currently* using pretty printing

* Would `Q` be ported? (Use `\itQ` for now)

* What to do with matrices? Using `Matrix{Sym}` with no `SymMatrix` type expected. Views seem off, so for now a copy is made.

* no new special functions exported, just the ones in `SpecialFunctions.jl`


## Installing `sympy`,

The Python `sympy` package should install with the package through `PythonCall` and `CondaPkg`. If not,
installing `sympy` in `PythonCall` isn't hard; shown below how it may be done.

```
using PythonCall
]
conda add python sympy
conda resolve
[delete key]
sympy = pyimport("sympy")
```
