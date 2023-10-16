# SymPyPythonCall

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://jverzani.github.io/SymPyCore.jl/dev)
[![Build Status](https://github.com/jverzani/SymPyPythonCall.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jverzani/SymPyPythonCall.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jverzani/SymPyPythonCall.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jverzani/SymPyPythonCall.jl)

[SymPyCore](https://github.com/jverzani/SymPyCore.jl) provides a `Julia`n interface to the [SymPy](https://www.sympy.org/) library of Python.


[SymPyPythonCall](https://github.com/jverzani/SymPyPythonCall.jl) utilizes `SymPyCore` and the `PythonCall` package (to provide the interop between `Julia` and `Python`) to enable access to Python's SymPy library using the practices and idioms of `Julia`.

The package [SymPyPyCall](https://github.com/jverzani/SymPyPyCall.jl) does a similar thing with the `PyCall` package providing the interop.



## Installation

Installing this package should install the `SymPyCore` package and the `PythonCall` package. The `PythonCall` package installs a `CondaPkg` package that installs the underlying SymPy library for Python.
