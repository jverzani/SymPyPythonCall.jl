
## Plotting

The `Plots` package allows many 2-dimensional plots of `SymPy` objects
to be agnostic as to a backend plotting package.  `SymPy` provides
recipes that allow symbolic expressions to be used where functions are
part of the `Plots` interface.
[See the help page for `sympy_plotting`.]

In particular, the following methods of `plot` are defined:

* `plot(ex::Sym, a, b)` will plot the expression of single variable over the interval `[a,b]`
* `plot!(ex::Sym, a, b)` will add to the current plot a plot of  the expression of single variable over the interval `[a,b]`
* `plot(exs::Vector{Sym}, a, b)` will plot each expression over `[a,b]`
* `plot(ex1, ex2, a, b)` will plot a parametric plot of the two expressions over the interval `[a,b]`.
* `contour(xs, ys, ex::Sym)` will make a contour plot of the expression of two variables over the grid specifed by the `xs` and `ys`.
* `surface(xs, ys, ex::Sym)` will make a surface plot of the expression of two variables over the grid specifed by the `xs` and `ys`.


For example:

```@example plots
using SymPy, Plots
@syms x
plot(x^2 - 2, -2,2)
savefig("plot-1.svg"); nothing  # hide
```

![](plot-1.svg)

Or a parametric plot:

```@example plots
plot(sin(2x), cos(3x), 0, 4pi);
savefig("plot-2.svg"); nothing # hide
```

![](plot-2.svg)

For plotting with other plotting packages, it is generally faster to
first call `lambdify` on the expression and then generate `y` values
with the resulting `Julia` function. An example might follow this pattern:

```@example plots
ex = cos(x)^2  +  cos(x^2)
fn = lambdify(ex)
xs = range(0, stop=10, length=256)
plot(xs, fn.(xs))
savefig("plot-3.svg"); nothing #hide
```

![](plot-3.svg)

----

In addition, with `PyPlot` a few other plotting functions from `SymPy` are available from its interface to `MatplotLib`:

* `plot3d_parametric_surface(ex1::Sym, ex2::Sym, ex3::Sym), (uvar, a0,
  b0), (vvar, a1, b1))` -- make a surface plot of the expressions
  parameterized by the region `[a0,b0] x [a1,b1]`. The default region
  is `[-5,5]x[-5,5]` where the ordering of the variables is given by
  `free_symbols(ex)`.

* `plot_implicit(predictate, (xvar, a0, b0), (yvar, a1, b1))` -- make
an implicit equation plot of the expressions over the region `[a0,b0]
x [a1,b1]`. The default region is `[-5,5]x[-5,5]` where the ordering
of the variables is given by `free_symbols(ex)`.  To create predicates
from the variable, the functions `Lt`, `Le`, `Eq`, `Ge`, and `Gt` can
be used, as with `Lt(x*y, 1)`. For infix notation, unicode operators
can be used: `\ll<tab>`, `\leqq<tab>`, `\Equal<tab>`, `\geqq<tab>`, and
`\gg<tab>`. For example, `x*y ≪ 1`.  To combine terms, the unicode
`\vee<tab>` (for "or"), `\wedge<tab>` (for "and") can be used.
