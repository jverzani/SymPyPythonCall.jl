# A SymPyCall introduction

!!! note "SymPyCall"

    This is a test of what works in `SymPyCall`.
	`SymPyCall` is just a temporary package for now, with the expectation that it may become `SymPy`.
	There would be a few deprecations:
	* `@vars` would be deprecated; use `@syms` only
	* `Base.show` isn't *currently* using pretty printing
	* `elements` for sets would be removed (convert to a `Set` by default)
	* Would `Q` be exported? (only the slant italic Q is currently)
	* What to do with matrices? E.g. `A\b` is failing now
	* `sympy.poly` *not* `sympy.Poly`
    * would we export `CommonEq`?
	* `limit(ex, x, c)` deprecated; use `limit(ex, x=>c)` or `sympy.limit`
    * no new special functions exported, just the ones in SpecialFunctions.jl

This document provides an introduction to using `SymPy` within `Julia`.
It owes an enormous debt to the tutorial for using SymPy within Python which may be found
[here](http://docs.sympy.org/dev/tutorial/index.html). The overall structure and many examples are taken from that work, with adjustments and additions to illustrate the differences due to using `SymPy` within `Julia`.


After installing `SymPy`, which is discussed in the package's `README`
file, we must first load it into `Julia` with the standard command
`using`:


```@setup introduction
using SymPyCall
sympy.init_printing(use_latex=true)
```

```jldoctest introduction
julia> using SymPyCall
```


## Symbols

At the core of `SymPy` is the introduction of symbolic variables that
differ quite a bit from `Julia`'s variables. Symbolic variables do not
immediately evaluate to a value, rather the "symbolicness" propagates
when interacted with. To keep things manageable, SymPy does some
simplifications along the way.

Symbolic expressions are primarily of the `Sym` type and can be constructed in the standard way:

```jldoctest introduction
julia> x = Sym("x")
x
```

This creates a symbolic object `x`, which can be manipulated through further function calls.


There is the `@syms` macro that makes creating multiple variables a
bit less typing, as it creates variables in the local scope -- no
assignment is necessary. Compare these similar ways to create symbolic
variables:

```jldoctest introduction
julia> @syms a b c
(a, b, c)

julia> a,b,c = Sym("a,b,c")
3-element Vector{Sym}:
 a
 b
 c
```

Here are two ways to make related variables:

```jldoctest introduction
julia> @syms xs[1:5]
(Sym[xs₁, xs₂, xs₃, xs₄, xs₅],)

julia> ys = [Sym("y$i") for i in 1:5]
5-element Vector{Sym}:
 y1
 y2
 y3
 y4
 y5
```

The former much more succinct, but the latter pattern of use when the number of terms is a variable.


The `@syms` macro is recommended, and will be modeled in the following, as it makes the specification of assumptions and symbolic functions more natural.

### Assumptions

Finally, there is the `symbols` constructor for producing symbolic objects. With `symbols` it is
possible to pass assumptions onto the variables. A list of possible
assumptions is
[here](http://docs.sympy.org/dev/modules/core.html#module-sympy.core.assumptions). Some
examples are:

```jldoctest introduction
julia> u = symbols("u")
u

julia> x = symbols("x", real=true)
x

julia> y1, y2 = symbols("y1, y2", positive=true)
2-element Vector{Sym}:
 y1
 y2

julia> alpha = symbols("alpha", integer=true, positive=true)
alpha

```

As seen, the `symbols` function can be used to make one or more variables with zero, one or more assumptions.

We jump ahead for a second to illustrate, but here we see that `solve` will respect these assumptions, by failing to find solutions to these equations:

```jldoctest introduction
julia> solve(x^2 + 1)   # ±i are not real
Any[]

```

```jldoctest introduction
julia> solve(y1 + 1)    # -1 is not positive
Any[]

```

The `@syms` macro allows annotations, akin to type annotations, to specify assumptions on new variables:

```jldoctest introduction
julia> @syms u1::positive u2::positive
(u1, u2)

julia> solve(u1 + u2)  # empty, though solving u1 - u2 is not.
Any[]
```

Additionally you can rename arguments using pair notation:
```
julia> @syms a1=>"α₁" a2=>"α₂"
(α₁, α₂)
```

In this example, the Julia variables `a1` and `a2` are defined to store SymPy
symbols with the "pretty" names `α₁` and `α₂` respectively.

As can be seen, there are several ways to create symbolic values, but
the recommended way is to use `@syms`. One caveat is that one can't
use `Sym` to create a variable from a function name in Base.

### Special constants

`Julia` has its math constants, like `pi` and `e`, `SymPy` as well. A few of these have `Julia` counterparts provided by `SymPy`. For example, these two constants are defined (where `oo` is for infinity):

```jldoctest introduction
julia> PI,  oo
(pi, oo)

```

(The pretty printing of SymPy objects does not work for tuples.)

Numeric values themselves can be symbolic. This example shows the
difference. The first `asin` call dispatches to `Julia`'s `asin`
function, the second to `SymPy`'s:

```jldoctest introduction
julia> [asin(1), asin(Sym(1))]
2-element Vector{Sym}:
 1.57079632679490
             pi/2

```


## Substitution

SymPy provides a means to substitute values in for the symbolic expressions. The specification requires an expression, a variable in the expression to substitute in for, and a new value. For example, this is one way to make a polynomial in a new variable:

```jldoctest introduction
julia> @syms x y
(x, y)

julia> ex = x^2 + 2x + 1
x^2 + 2*x + 1

julia> ex.subs(x, y)
y^2 + 2*y + 1

```


Substitution can also be numeric:

```jldoctest introduction
julia> ex.subs(x, 0)
1


```

The output has no free variables, but is still symbolic.

Expressions with more than one variable can have multiple substitutions, where each is expressed as a tuple:

```jldoctest introduction
julia> @syms x,y,z
(x, y, z)

julia> ex = x + y + z
x + y + z

julia> ex.subs([(x,1), (y, pi)])
z + 1 + pi
```

!!! note

    The calling pattern for `subs` is different from a typical `Julia` function call. The `subs` call is `object.method(arguments)` whereas a more "`Julia`n" function call is `method(objects, other objects....)`, as `Julia` offers multiple dispatch of methods. `SymPy` uses the Python calling method, adding in `Julia`n style when appropriate for generic usage within `Julia`. `SymPy` imports most all generic functions from the underlying `sympy` module and specializes them on a symbolic first argument.

    For `subs`, the simple substitution `ex.object(x,a)` is similar to simple function evaluation, so `Julia`'s call notation will work. To specify the pairing off of `x` and `a`, the `=>`  pairs notation is used.

This calling style will be equivalent to the last:

```jldoctest introduction
julia> ex(x=>1, y=>pi)
z + 1 + pi
```

A straight call is also possble, where the order of the variables is determined by `free_symbols`.
This is useful for expressions of a single variable, but being more explicit through the use of paired values is recommended.


## Conversion from symbolic to numeric

SymPy provides two identical means to convert a symbolic math
expression to a number. One is `evalf`, the other `N`. Within `Julia`
we decouple this, using `N` to also convert to a `Julian` value and
`evalf` to leave the conversion as a symbolic object.  The `N`
function converts symbolic integers, rationals, irrationals, and
complex values, while attempting to find an appropriate `Julia` type
for the value.

To see the difference, we use both on `PI`:

```jldoctest introduction
julia> N(PI)  # converts to underlying pi irrational
π = 3.1415926535897...

```

Whereas, `evalf` will produce a symbolic numeric value:

```jldoctest introduction
julia> (PI).evalf()
3.14159265358979

```


The `evalf` call allows for a precision argument to be passed through the second argument. This is how 30 digits of $\pi$ can be extracted:

```jldoctest introduction
julia> PI.evalf(30)
3.14159265358979323846264338328

```

This is a SymPy, symbolic number, not a `Julia` object. Composing with `N`

```jldoctest introduction
julia> N(PI.evalf(30))
3.14159265358979311599796346854419

```

will produce a `Julia` number,


Explicit conversion via `convert(T, ex)` can also be done, and is
necessary at times if `N` does not give the desired type.

## Algebraic expressions

`SymPy` overloads many of `Julia`'s functions to work with symbolic objects, such as seen above with `asin`. The usual mathematical operations such as `+`, `*`, `-`, `/` etc. work through `Julia`'s promotion mechanism, where numbers are promoted to symbolic objects, others dispatch internally to related `SymPy` functions.

In most all  cases, thinking about this distinction between numbers and symbolic numbers is unnecessary, as numeric values passed to `SymPy` functions are typically promoted to symbolic expressions. This conversion will take math constants to their corresponding `SymPy` counterpart, rational expressions to rational expressions, and floating point values to floating point values. However there are edge cases. An expression like `1//2 * pi * x` will differ from the seemingly identical  `1//2 * (pi * x)`. The former will produce a floating point value from `1//2 * pi` before being promoted to a symbolic instance. Using the symbolic value `PI` makes this expression work either way.

Most of `Julia`'s
[mathematical](http://julia.readthedocs.org/en/latest/manual/mathematical-operations/#elementary-functions)
functions are overloaded to work with symbolic expressions. `Julia`'s
generic definitions are used, as possible. This also introduces some
edge cases. For example, `x^(-2)` will balk due to the negative,
integer exponent, but either `x^(-2//1)` or `x^Sym(-2)` will work as
expected, as the former call first dispatches to a generic defintion,
but the latter two expressions do not.


`SymPy` makes it very easy to work with polynomial and rational expressions. First we create some variables:

```jldoctest introduction
julia> @syms x y z
(x, y, z)

```

### The expand, factor, collect, and simplify functions

A typical polynomial expression in a single variable can be written in two common ways, expanded or factored form. Using `factor` and `expand` can move between the two.

For example,

```jldoctest introduction
julia> p = x^2 + 3x + 2
x^2 + 3*x + 2

julia> factor(p)
(x + 1)*(x + 2)
```

Or

```jldoctest introduction
julia> expand(prod((x-i) for i in 1:5))
x^5 - 15*x^4 + 85*x^3 - 225*x^2 + 274*x - 120

```

The `factor` function factors over the rational numbers, so something like this with obvious factors is not finished:

```jldoctest introduction
julia> factor(x^2 - 2)
x^2 - 2

```

When expressions involve one or more variables, it can be convenient to be able to manipulate them. For example,
if we define `q` by:

```jldoctest introduction
julia> q = x*y + x*y^2 + x^2*y + x
x^2*y + x*y^2 + x*y + x

```

Then we can collect the terms by the variable `x`:

```jldoctest introduction
julia> collect(q, x)
x^2*y + x*(y^2 + y + 1)

```

or the variable `y`:

```jldoctest introduction
julia> collect(q, y)
x*y^2 + x + y*(x^2 + x)

```

These are identical expressions, though viewed differently.

A more broad-brush approach is to let `SymPy` simplify the values. In this case, the common value of `x` is factored out:

```jldoctest introduction
julia> simplify(q)
x*(x*y + y^2 + y + 1)

```

The `simplify` function attempts to apply the dozens of functions related to simplification that are part of SymPy. It is also possible to apply these functions one at a time, for example `sympy.trigsimp` does trigonometric simplifications.

The SymPy tutorial illustrates that `expand` can also result in simplifications through this example:

```jldoctest introduction
julia> expand((x + 1)*(x - 2) - (x - 1)*x)
-2

```


These methods are not restricted to polynomial expressions and will
work with other expressions. For example, `factor` identifies the
following as a factorable object in terms of the variable `exp(x)`:

```jldoctest introduction
julia> factor(exp(2x) + 3exp(x) + 2)
(exp(x) + 1)*(exp(x) + 2)

```

## Rational expressions: apart, together, cancel

When working with rational expressions, SymPy does not do much
simplification unless asked. For example this expression is not
simplified:

```jldoctest introduction
julia> r = 1/x + 1/x^2
1/x + x^(-2)

```

To put the terms of `r` over a common denominator, the `together` function is available:

```jldoctest introduction
julia> together(r)
(x + 1)/x^2

```

The `apart` function does the reverse, creating a partial fraction decomposition from a ratio of polynomials:

```jldoctest introduction
julia> apart( (4x^3 + 21x^2 + 10x + 12) /  (x^4 + 5x^3 + 5x^2 + 4x))
(2*x - 1)/(x^2 + x + 1) - 1/(x + 4) + 3/x

```

Some times SymPy will cancel factors, as here:

```jldoctest introduction
julia> top = (x-1)*(x-2)*(x-3)
(x - 3)*(x - 2)*(x - 1)

julia> bottom = (x-1)*(x-4)
(x - 4)*(x - 1)

julia> top/bottom
(x - 3)*(x - 2)/(x - 4)

```

(This might make math faculty a bit upset, but it is in line with student thinking.)

However, with expanded terms, the common factor of `(x-1)` is not cancelled:

```jldoctest introduction
julia> r = expand(top) / expand(bottom)
(x^3 - 6*x^2 + 11*x - 6)/(x^2 - 5*x + 4)

```

The `cancel` function instructs SymPy to perform cancellations. It
takes rational functions and puts them in a canonical $p/q$ form with
no common (rational) factors and leading terms which are integers:

```jldoctest introduction
julia> cancel(r)
(x^2 - 5*x + 6)/(x - 4)

```
## Powers

The SymPy [tutorial](http://docs.sympy.org/dev/tutorial/simplification.html#powers) offers a thorough explanation on powers and which get simplified and under what conditions. Basically

* $x^a x^b = x^{a+b}$ is always true. However

* $x^a y^a=(xy)^a$ is only true with assumptions, such as $x,y \geq 0$ and $a$ is real, but not in general. For example, $x=y=-1$ and $a=1/2$ has $x^a \cdot y^a = i \cdot i =  -1$, where as $(xy)^a = 1$.

* $(x^a)^b = x^{ab}$ is only true with assumptions. For example $x=-1, a=2$, and $b=1/2$ gives $(x^a)^b = 1^{1/2} = 1$, whereas $x^{ab} = -1^1 = -1$.


We see that with assumptions, the following expression does simplify to $0$:

```jldoctest introduction
julia> @syms x::nonnegatve y::nonnegative  a::real
(x, y, a)

julia> simplify(x^a * y^a - (x*y)^a)
0

```

However, without assumptions this is not the case

```jldoctest introduction
julia> @syms x,y,a
(x, y, a)

julia> simplify(x^a * y^a - (x*y)^a)
x^a*y^a - (x*y)^a

```

The `simplify` function calls `powsimp` to simplify powers, as above. The `powsimp` function has the keyword argument `force=true` to force simplification even if assumptions are not specified:

```jldoctest introduction
julia> sympy.powsimp(x^a * y^a - (x*y)^a, force=true)
0

```

## Trigonometric simplification

For trigonometric expressions, `simplify` will use `trigsimp` to simplify:

```jldoctest introduction
julia> @syms theta::real
(theta,)

julia> p = cos(theta)^2 + sin(theta)^2
sin(theta)^2 + cos(theta)^2

```

Calling either `simplify` or `trigsimp` will apply the Pythagorean identity:

```jldoctest introduction
julia> simplify(p)
1

```

While often forgotten,  the `trigsimp` function is, of course,  aware of the double angle formulas:

```jldoctest introduction
julia> simplify(sin(2theta) - 2sin(theta)*cos(theta))
0

```

The `expand_trig` function will expand such expressions:

```jldoctest introduction
julia> expand_trig(sin(2theta))
2*sin(theta)*cos(theta)

```


## Coefficients

Returning to polynomials, there are a few functions to find various pieces of the polynomials. First we make a general quadratic polynomial:

```jldoctest introduction
julia> @syms a,b,c,x
(a, b, c, x)

julia> p = a*x^2 + b*x + c
a*x^2 + b*x + c

```

If given a polynomial, like `p`, there are different means to extract the coefficients:

* SymPy provides a `coeffs` method for `Poly` objects, but `p` must first be converted to one.

* SymPy provides the `coeff` method for expressions, which allows extration of a coeffiecient for a given monomial




The `ex.coeff(monom)` call will return the corresponding coefficient of the monomial:

```jldoctest introduction
julia> p.coeff(x^2) # a
a

julia> p.coeff(x)   # b
b

```

The constant can be found through substitution:

```jldoctest introduction
julia> p(x=>0)
c

```

Though one could use some trick like this to find all the coefficients, that is cumbersome, at best.

```jldoctest introduction
julia> vcat([p.coeff(x^i) for i in N(degree(p,gen=x)):-1:1], [p(x=>0)])
3-element Vector{Sym}:
 a
 b
 c

```



SymPy has a function `coeffs`, but it is defined for polynomial types, so will fail on `p`:


```jldoctest introduction
julia> try p.coeffs() catch err "ERROR: KeyError: key `coeffs` not found" end # wrap p.coeffs() for doctest of error
"ERROR: KeyError: key `coeffs` not found"
```

Polynomials are a special class in SymPy and must be constructed. The `poly` constructor can be used. As there is more than one free variable in `p`, we specify the variable `x` below:

```jldoctest introduction
julia> q = sympy.poly(p, x)
Poly(a*x^2 + b*x + c, x, domain='ZZ[a,b,c]')

julia> q.coeffs()
3-element Vector{Sym}:
 a
 b
 c

```

!!! note

    XXX This is SYmPy not SymPyCall Either `sympy.poly` or `sympy.Poly` may be used. The `Poly` constructor from SymPy is *not* a function, so is not exported when `SymPy` is loaded. To access it, the object must be qualified by its containing module, in this case `Poly`. Were it to be used frequently, an alias could be used, as in `const Poly=sympy.Poly` *or* the `import_from` function, as in `import_from(sympy, :Poly)`. The latter has some attempt to avoid naming collisions.


## Polynomial roots: solve, real_roots, polyroots, nroots

SymPy provides functions to find the roots of a polynomial. In
general, a polynomial with real coefficients of degree $n$ will have
$n$ roots when multiplicities and complex roots are accounted for. The
number of real roots is consequently between $0$ and $n$.

For a *univariate* polynomial expression (a single variable), the real
roots, when available, are returned by `real_roots`. For example,

```jldoctest introduction
julia> real_roots(x^2 - 2)
2-element Vector{Sym}:
 -sqrt(2)
  sqrt(2)

```

Unlike `factor` -- which only factors over rational factors --
`real_roots` finds the two irrational roots here. It is well known
(the
[Abel-Ruffini theorem](http://en.wikipedia.org/wiki/Abel%E2%80%93Ruffini_theorem))
that for degree 5 polynomials, or higher, it is not always possible to
express the roots in terms of radicals. However, when the roots are
rational `SymPy` can have success:


```jldoctest introduction
julia> p = (x-3)^2*(x-2)*(x-1)*x*(x+1)*(x^2 + x + 1)
x*(x - 3)^2*(x - 2)*(x - 1)*(x + 1)*(x^2 + x + 1)

julia> real_roots(p)
6-element Vector{Sym}:
 -1
  0
  1
  2
  3
  3

```

!!! note "Why `string`?"

    XXX THIS ISN"T CURRENTLY NEEDED, so we have commented them out. XXX The uses of `string(p)` above and elsewhere throughout the introduction is only for technical reasons related to doctesting and how `Documenter.jl` parses  the expected output. This usage is not idiomatic, or suggested; it  only allows the cell  to  be tested programatically for  regressions. Similarly, expected errors  are  wrapped in `try`-`catch` blocks just  for testing purposes.


In this example, the degree of `p` is 8, but only the 6 real roots
returned, the double root of $3$ is accounted for. The two complex
roots of `x^2 + x+ 1` are not considered by this function. The complete set
of distinct roots can be found with `solve`:

```jldoctest introduction
julia> solve(p)
7-element Vector{Sym}:
                 -1
                  0
                  1
                  2
                  3
 -1/2 - sqrt(3)*I/2
 -1/2 + sqrt(3)*I/2

```

This finds the complex roots, but does not account for the double
root. The `roots` function of SymPy does.



The output of calling `roots` will be a dictionary whose keys are the roots and values the multiplicity.

```julia
julia> roots(p)
Dict{Any, Any} with 7 entries:
  -1                 => 1
  3                  => 2
  1                  => 1
  0                  => 1
  -1/2 - sqrt(3)*I/2 => 1
  2                  => 1
  -1/2 + sqrt(3)*I/2 => 1

```

When exact answers are not provided, the `roots` call is contentless:

```jldoctest introduction
julia> p = x^5 - x + 1
x^5 - x + 1

julia> sympy.roots(p)
Dict{Any, Any}()

```

Calling `solve` seems to produce very little as well:

```jldoctest introduction
julia> rts = solve(p)
5-element Vector{Sym}:
 CRootOf(x^5 - x + 1, 0)
 CRootOf(x^5 - x + 1, 1)
 CRootOf(x^5 - x + 1, 2)
 CRootOf(x^5 - x + 1, 3)
 CRootOf(x^5 - x + 1, 4)

```

But in fact, `rts` contains lots of information. We can extract numeric values quite easily with `N`:

```jldoctest introduction
julia> N.(rts)
5-element Vector{Number}:
                     -1.1673039782614187
 -0.18123244446987538 - 1.0839541013177107im
 -0.18123244446987538 + 1.0839541013177107im
   0.7648844336005847 - 0.35247154603172626im
   0.7648844336005847 + 0.35247154603172626im

```

These are numeric approximations to irrational values. For numeric
approximations to polynomial roots, the `nroots` function is also
provided. The answers are still symbolic:

```jldoctest introduction
julia> nroots(p)
5-element Vector{Sym}:
                       -1.16730397826142
 -0.181232444469875 - 1.08395410131771*I
 -0.181232444469875 + 1.08395410131771*I
 0.764884433600585 - 0.352471546031726*I
 0.764884433600585 + 0.352471546031726*I

```


## The solve function

The `solve` function is more general purpose than just finding roots of univariate polynomials. The function tries to solve for when an expression is 0, or a set of expressions are all 0.

For example, it can be used to solve when $\cos(x) = \sin(x)$:

```jldoctest introduction
julia> solve(cos(x) - sin(x))
1-element Vector{Sym}:
 pi/4

```

Though there are infinitely many correct solutions, these are within a certain range.

The
[solveset](http://docs.sympy.org/latest/modules/solvers/solveset.html)
function appears in version 1.0 of SymPy and is an intended
replacement for `solve`. Here we see it describes all solutions:

```jldoctest introduction
julia> u = solveset(cos(x) - sin(x))
Union(ImageSet(Lambda(_n, 2*_n*pi + 5*pi/4), Integers), ImageSet(Lambda(_n, 2*_n*pi + pi/4), Integers))

```

The output of `solveset` is a set, rather than a vector or
dictionary. To get the values requires some work. For *finite sets* we collect the elements
with `collect`, but first we must convert to a `Julia` `Set`:

```jldoctest introduction
julia> v = solveset(x^2 - 4)
Set{Sym} with 2 elements:
  2
  -2
```

```
XXX This is not current
julia> collect(Set(v...))
```

This composition is done in the `elements` function:

```
XXX jldoctest introduction not current
julia> elements(v)

```

!!! note "no elements"
    In `SymPyCall` the `elements` function is not used, just use `collect`, as in `collect(v)

```jldoctest introduction
julia> collect(v)
2-element Vector{Sym}:
  2
 -2
```


The `elements` function does not work for more complicated (non-finite) sets, such as `u`. For these, the `contains` method may be useful to query the underlying elements



Solving within Sympy has limits. For example, there is no symbolic solution here:

```jldoctest introduction
julia> try  solve(cos(x) - x)  catch err "error" end # wrap command for doctest of error
"error"
```

(And hence the error message generated.)

For such an equation, a numeric method would be needed, similar to the `Roots` package. For example:

```jldoctest introduction
julia> nsolve(cos(x) - x, 1)
0.739085133215161

```

Though it can't solve everything, the `solve` function can also solve
equations of a more general type. For example, here it is used to
derive the quadratic equation:

```jldoctest introduction
julia> @syms a::real, b::real, c::real
(a, b, c)

julia> p = a*x^2 + b*x + c
a*x^2 + b*x + c

julia> xs = solve(p, x);

julia> [simplify(p(x => xᵢ)) for xᵢ ∈ xs]
2-element Vector{Sym}:
 0
 0

```

The extra argument `x` is passed to `solve` so that `solve` knows
which variable to solve for.

The `solveset` function is similar:

```jldoctest introduction
julia> solveset(p, x); # Set with two elements
```


If the `x` value is not given, `solveset` will error and  `solve` will try to find a
solution over all the free variables:

```jldoctest introduction
julia> solve(p)
1-element Vector{Dict{Sym, Sym}}:
 Dict(a => (-b*x - c)/x^2)
```

Systems of equations can be solved as well. We specify them within a
vector of expressions, `[ex1, ex2, ..., exn]` where a found solution
is one where all the expressions are 0. For example, to solve this
linear system: $2x + 3y = 6, 3x - 4y=12$, we have:

```jldoctest introduction
julia> @syms x::real, y::real
(x, y)

julia> exs = [2x+3y-6, 3x-4y-12]
2-element Vector{Sym}:
  2*x + 3*y - 6
 3*x - 4*y - 12
```

```jldoctest introduction
julia> d = solve(exs); # Dict(x=>60/17, y=>-6/17)
```



We can "check our work" by plugging into each equation. We take advantage of how the `subs` function allows us to pass in a dictionary:

```jldoctest introduction
julia> map(ex -> ex.subs(d), exs)
2-element Vector{Sym}:
 0
 0

```



The more `Julia`n way to solve a linear  equation, like this   would be as follows:

```jldoctest introduction
julia> A = Sym[2 3; 3  -4]; b = Sym[6, 12]
2-element Vector{Sym}:
  6
 12
julia> A \ b
2-element Vector{Sym}:
 Fraction(60, 17)
 Fraction(-6, 17)

```

(Rather than use a generic  `lu` solver through `Julia` (which  proved slow for larger  systems),  the `\` operator utilizes  `solve` to perform this  computation.)



In the previous example, the system had two equations and two
unknowns. When that is not the case, one can specify the variables to
solve for as a vector. In this example, we find a quadratic polynomial
that approximates $\cos(x)$ near $0$:

```julia
julia> a,b,c,h = symbols("a,b,c,h", real=true)
(a, b, c, h)

julia> p = a*x^2 + b*x + c
   2
a⋅x  + b⋅x + c

julia> fn = cos
cos (generic function with 14 methods)

julia> exs = [fn(0*h)-p(x=>0), fn(h)-p(x => h), fn(2h)-p(x => 2h)]
3-element Vector{Sym}:
                           1 - c
       -a*h^2 - b*h - c + cos(h)
 -4*a*h^2 - 2*b*h - c + cos(2*h)

julia> d = solve(exs, (a,b,c))
Dict{Any, Any} with 3 entries:
  a => -cos(h)/h^2 + cos(2*h)/(2*h^2) + 1/(2*h^2)
  c => 1
  b => 2*cos(h)/h - cos(2*h)/(2*h) - 3/(2*h)

```

Again, a dictionary is returned. The polynomial itself can be found by
substituting back in for `a`, `b`, and `c`:

```julia
julia> quad_approx = p.subs(d)#; string(quad_approx)
"x^2*(-cos(h)/h^2 + cos(2*h)/(2*h^2) + 1/(2*h^2)) + x*(2*cos(h)/h - cos(2*h)/(2*h) - 3/(2*h)) + 1"

```

Taking the "limit" as $h$ goes to 0 produces the answer $1 - x^2/2$, as  will be shown.

Finally for `solve`, we show one way to re-express the polynomial $a_2x^2 + a_1x + a_0$
as $b_2(x-c)^2 + b_1(x-c) + b_0$ using `solve` (and not, say, an
expansion theorem.)

```julia
julia> n = 3
3

julia> @syms x, c
(x, c)

julia> @syms as[1:3]
(Sym[as₁, as₂, as₃],)

julia> @syms bs[1:3]
(Sym[bs₁, bs₂, bs₃],)

julia> p = sum([as[i+1]*x^i for i in 0:(n-1)]);

julia> q = sum([bs[i+1]*(x-c)^i for i in 0:(n-1)]);

julia> solve(p-q, bs)
Dict{Any, Any} with 3 entries:
  bs₁ => as₁ + as₂*c + as₃*c^2
  bs₂ => as₂ + 2*as₃*c
  bs₃ => as₃

```


### Solving using logical operators

The `solve` function does not need to just solve `ex = 0`. There are other means to specify an equation. Ideally, it would be nice to say `ex1 == ex2`, but the interpretation of `==` is not for this. Rather, `SymPy` introduces `Eq` for equality. So this expression

```jldoctest introduction
julia> solve(Eq(x, 1))
1-element Vector{Sym}:
 1

```

gives 1, as expected from solving `x == 1`.

In addition to `Eq`, there are `Lt`, `Le`, `Ge`, `Gt`. The Unicode
operators (e.g., `\leq`  and not  `\leq`)  are not aliased to these, but there are alternatives
`\ll[tab]`, `\leqq[tab]`, `\Equal[tab]`, `\geqq[tab]`, `\gg[tab]` and
`\neg[tab]` to negate.

So, the above could have been written with the following nearly identical expression, though it is entered with `\Equal[tab]`:

```jldoctest introduction
julia> solve(x ⩵ 1)
1-element Vector{Sym}:
 1

```

Here is an alternative way of asking a previous question on a pair of linear equations:

```julia
julia> x, y = symbols("x,y", real=true)
(x, y)

julia> exs = [2x+3y ⩵ 6, 3x-4y ⩵ 12]    ## Using \Equal[tab]
2-element Vector{Sym}:
  2⋅x + 3⋅y = 6
 3⋅x - 4⋅y = 12

julia> d = solve(exs)
Dict{Any, Any} with 2 entries:
  x => 60/17
  y => -6/17

```

Here  is  one other way  to  express  the same

```jldoctest introduction
julia> Eq.( [2x+3y,3x-4y], [6,12]) |>  solve == d
true
```
