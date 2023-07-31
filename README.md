# SymPyCall

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jverzani.github.io/SymPyCall.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jverzani.github.io/SymPyCall.jl/dev)
[![Build Status](https://github.com/jverzani/SymPyCall.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jverzani/SymPyCall.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jverzani/SymPyCall.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jverzani/SymPyCall.jl)


This is a start on what is needed to use `PythonCall` instead of `PyCall` for `SymPy.jl`.
At the moment, the expectation is that *if* that change proves desirable, this would become `SymPy`.

For now, there are some small design decisions from `SymPy` reflected here:

There would be a few deprecations:

* `@vars` would be deprecated; use `@syms` only

* `elements` for sets would be removed (convert to a `Set` by default)

* `sympy.poly` *not* `sympy.Poly`

* `limit(ex, x, c)` deprecated; use `limit(ex, x=>c)` or `sympy.limit`

* `Base.show` isn't *currently* using pretty printing

* Would `Q` be ported? (Use `\itQ` for now)

* What to do with matrices? Using `Matrix{Sym}` with no `SymMatrix` type expected. Views seem off, so for now a copy is made.

* no new special functions exported, just the ones in SpecialFunctions.jl


## Installing `sympy`,

The `sympy` package for Python should install with the package through `PythonCall` and `CondaPkg`. If not,
to install `sympy` in `PythonCall` isn't hard; below shows how it may be done.

```
using PythonCall
]
conda add python sympy
conda resolve
[delete key]
sympy = pyimport("sympy")
```
