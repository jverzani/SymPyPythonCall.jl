using Test
using SymPyPythonCall
import SymPyPythonCall.PythonCall: PyException

using SpecialFunctions
using LinearAlgebra
using Base.MathConstants

import PythonCall: Py

SymPy = SymPyPythonCall
func = SymPyPythonCall.Introspection.func
args = SymPyPythonCall.Introspection.args
import SymPyPythonCall: SymFunction
sympify = sympy.sympify
And = sympy.And
Heaviside = sympy.Heaviside
SymMatrix = Matrix{Sym}

## notes -- @vars => @syms; fix syntax

@testset "Core" begin
    ## Symbol creation
    x = Sym("x")
    #x = sym"x" # deprecated
    x = Sym(:x)
    #XXX x,y = Sym(:x, :y)
    x,y = symbols("x,y")

    @syms u1 u2 u3
    @syms u::positive
    @test length(solve(u+1)) == 0
    # make sure @syms defines in a local scope
    let
        @syms w
    end
    @test_throws UndefVarError isdefined(w)
    @syms a b c

    # Renaming with @syms
    @syms a=>"α₁"
    @test string((a.name).py) == "α₁"  # XXX <<--- modified

    ## extract symbols
    @syms z
    ex = x*y*z
    @test isa(free_symbols(ex), Vector{Sym})
    @test free_symbols(ex) == [x, y, z]

    ## number conversions
    @test Sym(2) == 2
    @test Sym(2.0) == 2.0
    @test Sym(2//1) == 2
    @test Sym(im) == 1im
    @test Sym(2im) == 2im
    @test Sym(1 + 2im) == 1 + 2im

    pi, e, catalan = Base.MathConstants.pi, Base.MathConstants.e, Base.MathConstants.catalan
    @test N(Sym(pi)) == pi
    @test N(Sym(ℯ)) ==  ℯ
    # @test N(Sym(catalan)) == catalan XXX <<--- expose catalan?


    ## function conversion
    f1 = convert(Function, x^2)
    @test f1(2) == 4

    ### Subs
    ## subs, |> (x == number)
    f = x -> x^2 - 2
    y = f(x)
    @test float(y.subs( x, 1)) == f(1)
    ### @test float( y |> subs(x,1) ) == f(1) no subs method

    ## interfaces
    ex = (x-1)*(y-2)
    @test ex.subs(x, 1) == 0
    @test ex.subs(((x,1),)) == 0
    @test ex.subs(((x,2),(y,2))) == 0

    @test subs(ex, x=>1) == 0       # removed
    @test subs(ex, x=>2, y=>2) == 0
    @test subs(ex, Dict(x=>1)) == 0
    @test ex(x=>1) == 0
    @test ex(x=>2, y=>2) == 0
    @test ex.subs(Dict(x=>1)) == 0

    ### doit
    @syms x g()
    @syms f()

    D = Differential(x)
    df = D(f(x))
    dfx = subs(df, f(x) => x^2)
    @test dfx.doit() == 2*x
    @test doit(dfx) == 2*x
    @test dfx |> doit == 2*x
    # use deep=true to force nested evaluations
    dgfx = g(dfx)
    @test dgfx.doit(deep=true) == g(2*x)
    @test doit(dgfx, deep=true) == g(2*x)
    @test dgfx |> doit(deep=true) == g(2*x)

    ## match, replace, xreplace, rewrite
    x,y,z = symbols("x, y, z")
    #a,b,c = map(sympy.Wild, (:a,:b,:c)) ## XXX Wild needs strings -- not symbols
    a,b,c = map(sympy.Wild, ("a", "b", "c")) ## XXX Wild needs strings -- not symbols
    # ## match: we have pattern, expression to follow `match`
    d = match(a^a, (x+y)^(x+y))
    @test d[a] == x+y
    d = match(a^b, (x+y)^(x+y))
    @test d[b] == x + y
    ex = (2x)^2
    pat = a*b^c
    d = match(pat, ex)
    @test d[a] == 4 && d[b] == x && d[c] == 2
    @test pat.xreplace(d) == 4x^2

    ## replace
    a = Wild("a")
    ex = log(sin(x)) + tan(sin(x^2))
    ##XXX    @test replace(ex, func(sin(x)), u ->  sin(2u)) == log(sin(2x)) + tan(sin(2x^2))
    @test replace(ex, func(sin(x)), func(tan(x))) == log(tan(x)) + tan(tan(x^2))
    @test replace(ex, sin(a), tan(a)) ==  log(tan(x)) + tan(tan(x^2))
    @test replace(ex, sin(a), a) == log(x) + tan(x^2)
    @test replace(x*y, a*x, a) == y

    ## xreplace
    @test (1 + x*y).xreplace(Dict(x => PI)) == 1 + PI*y
    @test (x*y + z).xreplace(Dict(x*y => PI)) == z + PI
    @test (x*y * z).xreplace(Dict(x*y => PI)) == x* y * z
    @test (x + 2 + exp(x + 2)).xreplace(Dict(x+2=>y)) == x + exp(y) + 2


    # Test subs on simple numbers
    @syms x y
    @test Sym(2)(x=>300, y=>1.2) == 2

    #Test subs for pars and dicts
    ex = 1
    dict1 = Dict{String,Any}()
    dict2 = Dict{Any,Any}()
    #test subs
    for i=1:4
        x = Sym("x$i")
        ex=ex*x
        dict2[x] = i
        dict1[string(x)] = i
    end
    for d in [dict1, dict2]
        @test ex |> subs(d) == factorial(4)
        @test subs(ex, d) == factorial(4)
        @test subs(ex, d...) == factorial(4)
        @test ex |> subs(d...) == factorial(4)
        @test ex(d) == factorial(4)
        @test ex(d...) == factorial(4)
    end

    a = Sym("a")
    b = Sym("b")
    line = x -> a + b * x
    sol = solve([line(0)-1, line(1)-2],[a,b])
    ex = line(10)
    @test ex(sol) == ex(sol...) == 11

    ## Simplify (issue 343)
    @syms x
    @test simplify(sin(x)^2 + cos(x)^2) == 1
    @test simplify(sympy.gamma(x)/sympy.gamma(x-2)) == (x-1)*(x-2)
    _ones = (1, 1.0, big"1", "1", one)
    @test simplify.(_ones) == _ones

    ## Conversion
    x = Sym("x")
    p = x.subs(x,pi)
    q = x.subs(x,Sym(1)/2)
    r = x.subs(x,1.2)
    z = x.subs(x,1)
    @test isa(N(p), Irrational) # isa(N(p), Float64)
##XXX    @test isa(N(p, 60), BigFloat)
    @test isa(p.evalf(), Sym)
    @test isa(N(x), Sym)
    @test isa(N(q), Rational)
    @test isa(N(r), Float64)
    @test isa(N(z), Integer)

    ## method calls via getproperty
    p = (x-1)*(x-2)
    @test sympy.roots(p) == Dict{Any,Any}(Sym(1) => 1, Sym(2)=> 1) # sympy.roots
    p = sympy.poly(p, x) # XXX sympy.Poly(p, x)
    @test p.coeffs() == Any[1,-3,2] # p.coeffs

    ## algebra
    @test expand((x + 1)*(x + 2)) == x^2 + 3x + 2  # v0.7 deprecates expand, in v1.0 this is fine w/o qualifacation
    x1 = (x + 1)*(x + 2)
    @test expand(x1) == x^2 + 3x + 2
    @test sympy.expand_trig(sin(2x)) == 2sin(x)*cos(x)

    ## math functions
    u = abs(x^2 - 2)
    @test u(x=>0) == 2
    u = min(x, x^2, x^3, x^4)
    @test u(x=>2) == 2
    @test u(x=>1//2) == 1//2^4

    ## solve
    x,y,a = symbols("x,y,a", real=true)
    solve(x^2 - 2x)
    solve(x^2 - 2a, x)
    solve(x^2 - 2a, a)
    solve(Lt(x-2, 0))
    solve( x-2 ≪ 0)
    exs = [x-y-1, x+y-2]
    di = solve(exs)
    @test di[x] == 3//2
    @test map(ex -> subs(ex, di), exs) == [0,0]
    solve([x-y-a, x+y], [x,y])

    ## linsolve
    @syms x y
    M=Sym[1 2 3; 2 3 4]
    as = linsolve(M, x, y)

    # XXX elements -> collect
    @test length(collect(as)) == 1 # XXXlength(elements(as)) == 1
    @syms a b; eqs = (a*x+2y-3, 2b*x + 3y - 4)
    as = linsolve(eqs, x, y)
    @test length(collect(as)) == 1 # XXX length(elements(as)) == 1

    ## nsolve -- not method for arrays, issue 268
    @syms z1::positive z2z1::positive
    # XXX No error? @test_throws MethodError nsolve([z1^2-1, z1+z2z1-2], [z1,z2z1], [1,1])  # non symbolic first argument
    @test all(N.(sympy.nsolve([z1^2-1, z1+z2z1-2], [z1,z2z1], (1,1))) .≈ [1.0, 1.0])

    ## issue 355: direct definition of LinearAlgebra.:\
    @test_throws SingularException Sym[1 1; 1 1] \ [1, 2]
    @syms  a b c d e f
    A, b= [a b; c d], [e, f]
    x = A \ b
    @test  simplify.(A*x-b) == [0,0]
    #out =  Sym[1 1; 1 1] \ [1,1] # XXX non-singular
    #u =   free_symbols(out)[1]
    #@test out ==  [1-u,u]


    # Just a made-up example to test  if manageable
    @syms a1 a2
    n = 7
    A = diagm(0 => ones(Int, n), 1 => fill(a1, n-1), 2 => fill(1, n-2), -1 => fill(a2, n-1))
    b = Vector{Sym}(1:n)
    x = A \ b
    @test length(free_symbols(x)) ==  2


    ## limits
    @syms x
    # @test limit(x -> sin(x)/x, 0) == 1 # XXX deprected; use x=>c
    #@test limit(sin(x)/x, x, 0) |> float == 1 #XXX deprecated
    @test limit(sin(x)/x, x => 0) == 1
    @syms x h
    out = limit((sin(x+h) - sin(x))/h, h => 0) ## XXX limit((sin(x+h) - sin(x))/h, h, 0)
    @test (out.replace(x, pi) |> float) == -1.0


    ## diff
    diff(sin(x), x)
    out = diff(sin(x), x, 2)
    @test abs((out.replace(x, pi/4) |> float) - - sin(pi/4)) < sqrt(eps())

    # partial derivatives
    @syms x y
    @test diff(x^2 + x*y^2, x, 1) == 2x + y^2

    t = symbols("t", real=true)     # vector-valued functions
    r1(t) = [sin(t), cos(t), t]
    u = r1(t)
    kappa = norm(diff.(u) × diff.(u,t,2)) / norm(diff.(u))^3 |> simplify
    @test kappa == 1//2 # XXX convert(Rational,kappa) == 1//2

    u = SymFunction("u")
    eqn = Eq(x^2 + u(x)^2, x^3 - u(x))
    # diff(eqn) doesn't evaluate over Eq:
    @test Bool(func(eqn)(diff.(args(eqn))...) == Eq(2x + 2u(x) * diff(u(x),x), 3x^2 - diff(u(x),x))) ## XXXfunc(eqn)(diff.(args(eqn))...) == Eq(2x + 2u(x) * diff(u(x),x), 3x^2 - diff(u(x),x))

    ## integrate
    @test integrate(sin(x)) == -cos(x)
    @test integrate(sin(x), (x, 0, pi)) == 2.0
    a, b, t = symbols("a, b, t")
    @test integrate(sin(x), (x, a, b)) == cos(a) - cos(b)
    @test integrate(sin(x), (x, a, b)).replace(a, 0).replace(b, pi) == 2.0
    @test integrate(sin(x) * sympy.DiracDelta(x)) == sin(Sym(0)) # XXX integrate(sin(x) * DiracDelta(x)) == sin(Sym(0))
    @test integrate(sympy.Heaviside(x), (x, -1, 1)) == 1 ## XXX integrate(Heaviside(x), (x, -1, 1)) == 1
    #= not working XXX
    curv = sympy.Curve([exp(t)-1, exp(t)+1], (t, 0, log(Sym(2))))
    @test sympy.line_integrate(x + y, curv, [x,y]) == 3 * sqrt(Sym(2)) ## XXX line_integrate(x + y, curv, [x,y]) == 3 * sqrt(Sym(2))
    =#

    ## summation
    summation(1/x^2, (x, 1, 10))
    out = summation(1/x^2, (x, 1, 10))
    out1 = sum([1//x^2 for  x in 1:10])
    @test out.p == out1.num #XXX round(Integer, out.p) == out1.num
    @test out.q == out1.den ## XXX round(Integer, out.q) == out1.den


    ## Ops
    s = 3
    x = Sym("x")
    v = [x, 1]
    rv = [x 1]
    a = [x 1; 1 x]
    b = [x 1 2; 1 2 x]

    DIMERROR =  DimensionMismatch
    DimensionOrMethodError =  Union{MethodError, DimensionMismatch}
    ## scalar, [vector, matrix]
    @test s .+ v == [x+3, 4]
    @test v .+ s == [x+3, 4]
    @test s .+ rv == [x+3 4]
    @test rv .+ s == [x+3 4]
    @test s .+ a == [x+3 4; 4 x+3]
    @test a .+ s == [x+3 4; 4 x+3]

    @test s .- v == [3-x, 2]
    @test v .- s == [x-3, -2]
    @test s .- rv == [3-x 2]
    @test rv .- s == [x-3 -2]
    @test s .- a == [3-x 2; 2 3-x]
    @test a .- s == [x-3 -2; -2 x-3]

    2v
    2rv
    2a
    s .* v
    v .* s
    s .* rv
    rv .* s
    s .* a
    a .* s

    ## s / v  ## broadcasts s Depreacated
    @test s ./ v == [3/x, 3]
    @test v / s == [x/3, Sym(1)/3]
    @test v .\ s == s ./ v
    @test s \ v == v / s
    ## s / rv ## broadcasts s Deprecated
    @test s ./ rv == [3/x 3]
    @test rv / s == [x/3 Sym(1)/3]
    @test rv .\ s == s ./ rv
    @test s \ rv == rv / s
    ## s / a ## broadcasts s Deprecated
    @test s ./ a == [3/x 3; 3 3/x]
    @test a / s == [x/3 Sym(1)/3; Sym(1)/3 x/3]
    @test a .\ s == s ./ a
    @test s \ a == a / s

    @test_throws MethodError  s ^ v ## error
    @test s .^ v == [3^x, 3]
    @test_throws DimensionOrMethodError v ^ s ## error
    v .^ s
    @test_throws MethodError  s ^ rv ## error
    s .^ rv
    @test_throws DimensionMismatch  rv ^ s ## error
    rv .^ s
    @test_throws MethodError  s ^ a ## error
    s .^ a
    a ^ s
    a .^ s


    ## vector vector
    @test v .+ v == 2v
    ##@test_throws MethodError  v .+ rv ##  no longer an error, broadcase
    @test v .- v == [0, 0]
    @test_throws DIMERROR  v - rv
    @test_throws DimensionOrMethodError   v * v ## error
    @test v .* v == [x^2,1]
    @test dot(v, v) == 1 + conj(x)*x
    v * rv
    rv * v ## 1x2 2 x 1 == 1x1
    v .* rv ## XXX ?? should be what? -- not 2 x 2
    rv .* v ## XXX ditto
    ## @test_throws MethodError  v / v ## error
    v ./ v ## ones()
    v .\ v
    ## @test_throws MethodError  v / rv ## error
    v ./ rv  ## ??
    rv .\ v
    @test_throws MethodError  v ^ v ## error
    v .^ v
    @test_throws MethodError  v ^ rv ## error
    v .^ rv ## ??


    ## vector matrix
    @test_throws DIMERROR v + a ## error (Broadcast?)
    @test_throws DIMERROR a + v ## error
    v .+ a ## broadcasts
    a .+ v
    @test_throws DIMERROR  v - a ## error
    v .- a
    @test_throws DimensionMismatch  v * a ## error
    v .* a
    #@test_throws MethodError  v / a ## error
    v ./ a
    a .\ v
    @test_throws MethodError  v ^ a ## error
    v .^ a
    v
    ## matrix matrix
    a + a
    @test_throws DIMERROR  a + b ## error
    a + 2a
    a - a
    @test_throws DIMERROR  a - b ## error
    a * a
    a .* a
    a * b ## 2x2 * 2*3 -- 2x3
    @test_throws DIMERROR  a .* b ## error -- wrong size
    #@test_throws MethodError  a / a
    a ./ a ## ones
    a .\ a
    ##@test_throws MethodError  a / b ## error
    @test_throws DIMERROR  a ./ b ## error
    @test_throws DIMERROR  b .\ a ## error
    @test_throws MethodError  a ^ a ## error
    a .^ a
    @test_throws MethodError  a ^ b ## error
    @test_throws DIMERROR  a .^ b ## error


    ## Number theory
    #@test isprime(100) == isprime(Sym(100))
    #@test factorint(Sym(100)) == factor(100)
    @test sympy.prime(Sym(100)) == 541 # XXX prime(Sym(100)) == 541
    @test sympy.multiplicity(Sym(10), 100) == 2 # XXX multiplicity(Sym(10), 100) == 2


    ## polynomials
    @syms x y
    f1 = 5x^2  + 10x + 3
    g1 = 2x + 2
    q,r = sympy.div(f1,g1, domain="QQ") # use sympy.div to dispatch; o/w we can't disambiguate div(Sym(7), 5)) to do integer division
    @test r == Sym(-2)
    @test simplify(q*g1 + r - f1) == Sym(0)
    ## sympy.interpolate as first arg is not symbolic
    @test sympy.interpolate([1,2,4], x) == sympy.interpolate(collect(zip([1,2,3], [1,2,4])), x)
    @test sympy.interpolate(collect(zip([-1,0,1], [0,1,0])), x) == 1 - x^2

    ## piecewise
    x = Sym("x")
    # sympy.Piecewise is a FunctionClass,  we qualify, as args not Symbolic
    p = sympy.Piecewise((x, Ge(x,0)), (0, Lt(x,0)), (1, Eq(x,0)))
    ## using infix \ll<tab>, \gt<tab>, \Equal<tab>
    p = sympy.Piecewise((x, (x ≫ 0)), (0, x ≪ 0), (1, x ⩵ 0))
    @test p.subs(x,2) == 2
    @test p.subs(x,-1) == 0
    @test p.subs(x,0) == 1

    ## if VERSION < v"0.7.0-" # ifelse changed
    ##     u = ifelse(Lt(x, 0), "neg", ifelse(Gt(x, 0), "pos", "zero"))
    ##     @test u.subs(u,x,-1) == Sym("neg")
    ##     @test subs(u,x, 0) == Sym("zero")
    ##     @test subs(u,x, 1) == Sym("pos")
    ## end
    p = sympy.Piecewise((-x, x ≪ 0), (x, x ≧ 0))


    ## relations
    x, y = symbols("x, y")
    ex  = Eq(x^2, x)
    @test ex.lhs == x^2
    @test ex.rhs == x
    @test args(ex) == (x^2, x)

    # alternative operators
    for (ex1,ex2) ∈ ((Eq(x^2, x), x^2 ⩵ x),
                   (Eq(x^2, x), x^2 ~ x),
                   (Lt(x^2, x), x^2 ≪ x),
                   (Le(x^2, x), x^2 ≦ x),
                   (Ge(x^2, x), x^2 ≧ x),
                   (Gt(x^2, x), x^2 ≫ x))
        @test lhs(ex1) == lhs(ex2)
        @test rhs(ex1) == rhs(ex2)
    end

    ## mpmath functions
#    if @isdefined mpmath
#    if isdefined(SymPy, :mpmath)
        x = Sym("x")
        Sym(big(2))
        Sym(big(2.0))                   # may need mpmath (e.g., conda install mpmath)

    # XXXsympy.__version__ != "1.9" && (@test limit(besselj(Sym(1),1/x), x => 0) == Sym(0))
    string(sympy.py.__version__) != "1.9" && (@test limit(besselj(Sym(1),1/x), x => 0) == Sym(0))
    complex(N(sympy.hankel2(2, pi)))    # XXX complex(N(SymPy.mpmath.hankel2(2, pi)))

    #XXX        SymPy.mpmath.bei(2, 3.5)
    # XXX       SymPy.mpmath.bei(1+im, 2+3im)
#    end

    ## Assumptions
    @test ask(𝑄.even(Sym(2))) == true
    @test ask(𝑄.even(Sym(3))) == false
    @test ask(𝑄.nonzero(Sym(3))) == true
    @syms x_real::real
    @syms x_real_positive::(real, positive)
    @test ask(𝑄.positive(x_real)) == nothing
    @test ask(𝑄.positive(x_real_positive)) == true
    @test ask(𝑄.nonnegative(x_real^2)) == true
    ## XXX @test ask(𝑄.upper_triangular([x_real 1; 0 x_real])) == true
    @test ask(𝑄.positive_definite([x_real 1; 1 x_real])) == nothing


    ## sets
    ## XXX -- Sets need work
    s = sympy.FiniteSet("H","T")
    s1 = Py(s).powerset() # XXX s1 = s.powerset()
    # XXX VERSION >= v"0.4.0" && @test length(collect(convert(Set, s1))) == length(collect(s1.__pyobject__))
    a, b = sympy.Interval(0,1), sympy.Interval(2,3)
    @test a.is_disjoint(b) == true
    @test a.union(b).measure == 2



    ## test cse output
    @test sympy.cse(x) == (Any[], Sym[x])
    @test sympy.cse([x]) == (Any[], [x]) # XXXsympy.cse([x]) == (Any[], [reshape([x],1,1)])
    @test sympy.cse([x, x]) == sympy.cse([x, x]) == (Any[],  [x,x] ) # XXX (Any[],  [reshape([x,x], 2, 1)] )
    @test sympy.cse([x x; x x]) == (Any[], [[x x; x x]])

    ## sympy"..."(...)
    ## removed
    #@syms x
    #@test sympy"sin"(1) == sin(Sym(1))
end

@testset "Syms macro" begin
    @syms u
    @test isa(u, Sym)

    ret = @syms a, b, c
    @test isa(ret, Tuple{Sym, Sym, Sym})

    @syms x::(real,positive)=>"x₀", y, z::complex, n::integer
    @test isa(x, Sym)
    @test ask(sympy.And(𝑄.real(x), 𝑄.positive(x))) # XXX ask(And(𝑄.real(x), 𝑄.positive(x)))
    @test string(x) == "x₀"

    @test isa(y, Sym)

    @test isa(z, Sym)
    @test ask(𝑄.complex(z))

    @test isa(n, Sym)
    @test ask(𝑄.integer(n))

    @syms f()::(real, positive), g(), h()::complex=>"h̄"
    @test isa(f, SymFunction)
    @test ask(And(𝑄.real(f(x)), 𝑄.positive(f(x))))

    @test isa(g, SymFunction)

    @test isa(h, SymFunction)
    @test ask(𝑄.complex(h(x)))
    @test string(h) == "h̄"

    @syms X[1:20]
    @test isa(X, Vector{Sym})
    @test size(X) == (20,)
    @test string(X[11]) == "X₁₁"

    @syms bigy[1:5]=>"Y"
    @test string(bigy[3]) == "Y₃"

    @syms Z[1:5, 1:6]
    @test isa(Z, Matrix{Sym})
    @test size(Z) == (5, 6)
    @test string(Z[2,4]) == "Z₂_₄"

    ## XXX
    #@syms F[1:2](), G()[1:2]
    #@test isa(F, Vector{SymFunction})
    #@test isa(G, Vector{SymFunction})

    #@syms WOW[1:3, 1:2:4]()::(real, positive)=>"f"
    #@test isa(WOW, Matrix{SymFunction})
    #@test size(WOW) == (3, 2)
    #@test ask(And(𝑄.real(WOW[1,2](x)), 𝑄.positive(WOW[1,2](x))))
    #@test string(WOW[1,2]) == "f₁_₃"
end

@testset "SymFunctions" begin
    @syms x::real y::real
    #@symfuns f g real=true
    @syms f()::real g()::real

    @test isreal(f(x))
    @test isreal(g(y))

    #@symfuns h real=true positive=true
    @syms h()::(real, positive)
    @test ask(𝑄.positive(h(x))) # XXX h(x)>0
end

@testset "Fix past issues" begin
    @syms x y z
    ## Issue # 56
    @test Sym(1+2im) == 1+2IM
    @test convert(Sym, 1 + 2im) == 1 + 2IM


    ## Issue #59
    sympy.cse(sin(x)+sin(x)*cos(x))
    sympy.cse([sin(x), sin(x)*cos(x)])
    sympy.cse( [sin(x), sin(x)*cos(x), cos(x), sin(x)*cos(x)])

    ## Issue #60, lambidfy
    x, y = symbols("x, y")
    lambdify(sin(x)*cos(2x) * exp(x^2/2))
    fn = lambdify(sin(x)*asin(x)*sinh(x)); fn(0.25)
    lambdify(real(x)*imag(x))
    @test lambdify(min(x,y))(3,2) == 2 # XXX lambdify(Min(x,y))(3,2) == 2

    ex = 2 * x^2/(3-x)*exp(x)*sin(x)*sind(x)
    fn = lambdify(ex); map(fn, rand(10))
    ex = x - y
    @test lambdify(ex, (x,y))(3,2) == 1

    Indicator(x, a, b) = sympy.Piecewise((1, Lt(x, b) & Gt(x,a)), (0, Le(x,a)), (0, Ge(x,b)))
    i = Indicator(x, 0, 1)
    u = lambdify(i)
    @test u(.5) == 1
    @test u(1.5) == 0

    # SymPy issue 567; constants
    u = lambdify(Sym(1//2))
    @test u() == u(1,2,3) == 1/2
    @syms x
    ex = integrate(sqrt(1 + (1/x)^2), (x, 1/sympy.E, sympy.E))
    @test lambdify(ex)() ≈ 3.1961985135995072

#    i2 = SymPy.lambdify_expr(x^2,name=:square)
#    @test i2.head == :function
#    @test i2.args[1].args[1] == :square
    ## @test i2.args[2] == :(x.^2) # too fussy


    ## issue #67
    @test N(Sym(4)/3) == 4//3
    @test N(convert(Sym, 4//3)) == 4//3

    ## issue #71
    @test log(Sym(3), Sym(4)) == log(Sym(4)) / log(Sym(3))

    ## issue #103 # this does not work for `x` (which has `classname(x) == "Symbol"`), but should work for other expressions
    for ex in (sin(x), x*y^2*x, sqrt(x^2 - 2y))
        @test Bool(func(ex)(SymPy.Introspection.args(ex)...) == ex) # XXX func(ex)(SymPy.Introspection.args(ex)...) == ex
    end

    ## properties (Issue #119)
    @test (sympify(3).is_odd) == true
    @test sympy.poly(x^2 -2, x).is_monic == true # not Poly

    ## test round (Issue #153)
    y = Sym(eps())
    @test round(N(y), digits=5) == 0
    @test round(N(y), digits=16) != 0

    ## lambdify over a matrix #218
    @syms x y
    s = [1 0; 0 1]
    @test lambdify(x*s)(2) == 2 * s
    U = [x-y x+y; x+y x-y]
    @test lambdify(U, [x,y])(2,4) == [-2 6;6 -2]
    @test lambdify(U, [y,x])(2,4) == [ 2 6;6  2]

    @test eltype(lambdify([x 0; 1 x])(0)) <: Integer
    @test eltype(lambdify([x 0; 1 x], T=Float64)(0)) == Float64

    # issue 222 type of eigvals
    A = [Sym("a") 1; 0 1]
    @test typeof(eigvals(A)) <: Dict ## XXX typeof(eigvals(A)) <: Vector{Sym}

    # issue 231 Q.complex
    @syms x_maybecomplex
    @syms x_imag::imaginary
    @test ask(𝑄.complex(x_maybecomplex)) == nothing
    @test ask(𝑄.complex(x_imag)) == true

    # issue 242 lambdify and conjugate
    @syms x
    #expr = conjugate(x)
    expr = conj(x)
    fn = lambdify(expr)
    @test fn(1.0im) == 0.0 - 1.0im
    fn = lambdify(expr, use_julia_code=true)
    @test fn(1.0im) == 0.0 - 1.0im

    # issue 245 missing sincos
    @test applicable(sincos, x)
    @test sincos(x)[1] == sin(x)

    # issue 256 det
    @syms rho::real phi::real theta::real
    xs = [rho*cos(theta)*sin(phi), rho*sin(theta)*sin(phi), rho*cos(phi)]
    J = [diff(x, u) for x in xs, u in (rho, phi, theta)]
    J.det()

    # issue #273 x[i]
    x = sympy.IndexedBase("x")
    i,j = sympy.symbols("i j", integer=true)
    ## XXX not defined
    #@test x[i] == PyCall.py"$x[$i]"


    # issue 298 lambdify for results of dsolve
    @syms t
    F = SymFunction("F")
    diffeq = diff(F(t),t) - 3*F(t)
    res = dsolve(diffeq, F(t), ics=Dict(F(0) => 2))  # 2exp(3t)
    @test lambdify(res)(1) ≈ 2*exp(3*1)

    # issue 304 wrong values for sind, ...
    a = Sym(45)
    @test sind(a) == sin(PI/4)

    # issue #319  with   use   of  Dummy, but
    #  really  a  lambdify issue
    dummy = sympy.Dummy
    # Symbolic differentiation of functions
    function D(f)
        x = dummy("x")
        lambdify(diff.(f(x), x), (x,))
    end
    @test D(t -> t^2)(1) == 2

    # issue   #320  with integrate(f) when
    # f  is consant
    # XXX @test integrate(x -> 1, 0, 1)  == 1 <-- want to deprecate
    # XXX @test limit(x->1,  x, 0) == 1
    # xxx @test diff(x->1)  ==   0

    ## Issue 324 with inference of matrix operations
    A = fill(Sym("a"), 2, 2)
    @test eltype(A*A) == Sym
    @test eltype(A*ones(2,2)) == Sym
    @test eltype(A*Diagonal([1,1])) == Sym
    VERSION >= v"1.2.0"  && @test eltype(A * I(2)) == Sym

    ## Issue 328 with E -> e
    @syms x
    ex = 3 * sympy.E * x
    fn = lambdify(ex)
    @test fn(1) ≈ 3*exp(1) * 1

    ## Issue 332 with `abs2`
    @syms x::real
    @test abs2(x) == x*x
    @syms x
    @test abs2(x) == x*conj(x)

    ## Issue 376 promote to Sym Before pycall
    f = x -> x^2 + 1 +log(abs( 11*x-15 ))/99
    @test limit(f(x), x=>15//11) == -oo

    ## Issue #390 on div (__div__ was depracated, use __truediv__)
    #XXX@test Sym(2):-Sym(2):-Sym(2) |> collect == [2, 0, -2]

    ## Lambda function to create a lambda
    # XXX not working!
    #@syms x
    #ex = x^2 - 2
    #fn1 = Lambda(x, ex)
    #fn2 = lambdify(ex)
    #@test fn1(3) == fn2(3)

    ## issue 402 with lamdify and Order
    @syms x
    t = series(exp(x), x, 0, 2)
    @test lambdify(t)(1/2) == 1 + 1/2

    ## issue #411 with Heaviside
    @syms t
    u = Heaviside(t)
    λ = lambdify(u)
    @test all((iszero(λ(-1)), isone(λ(1))))
    VersionNumber(string(sympy.py.__version__)) >= v"1.9" && @test λ(0) == 1//2
    u = Heaviside(t, 1)
    λ = lambdify(u)
    @test all((iszero(λ(-1)), isone(λ(0)), isone(λ(1))))

    ## issue #295 with piecewise function
    @syms x
    p = sympy.Piecewise((x,Gt(x,0)), (x^2, Le(x,0)))
    @test lambdify(p)(2) == 2
    @test lambdify(p)(-2) == (-2)^2

    ## Issue catch all for N; floats only
    @syms x
    ex = integrate(sqrt(1 + (1/x)^2), (x, 1/sympy.E, sympy.E))
    @test N(ex) ≈ 3.196198513599507
end

@testset "generic programming, issue 223" begin
    # arose in issue 223
    @syms xreal::real
    @syms xcomplex
    zreal = sympify(1)
    zcomplex = sympify(1) + sympify(2)*IM

    @test isreal(xreal)     # is_real(xreal) is also true, but xreal is Sym, not a Julia object
    @test !isreal(xcomplex) # is_real(xcomplex) is nothing
    @test isreal(zreal)
    @test !isreal(zcomplex)

    # conversions
    @test complex(xreal) == xreal
    @test complex(xreal, xreal) == xreal + IM*xreal
    @test complex(xcomplex) != xcomplex
    @test complex(zreal) == zreal
    @test complex(zreal) !== zreal      # Complex{Int} !== Sym
    @test complex(zcomplex) == zcomplex
    @test complex(zcomplex) !== zcomplex

    ## issue 284 N(PI,50)
    ## XXX @test N(PI, 50) ≈ pi <--- issue with pyconvert????
    @test length(string(N(PI,50))) == 50 # XXX 2 + 50
    ## issue 284 lambdify of Pi
    ## XXX mpi = SymPy.PyCall.pyimport("sympy.parsing.mathematica")."mathematica"("Pi")
    #= XXX
    mpi = PythonCall.pyimport("sympy.parsing.mathematica").mathematica("Pi")
    @test SymPy.walk_expression(mpi) == :pi
    =#
    @test lambdify(PI^4*xreal)(256) == 256 * pi^4


    ## Issue 351 booleans  and arithmetic operations
    @test Sym(1) + true == Sym(2) == true +  Sym(1)
    @test Sym(1) - true == Sym(0) == true -  Sym(1)
    @test Sym(1) * true == Sym(1) == true * Sym(1)
    @test Sym(1) / true == Sym(1) == true / Sym(1)
    @test true^Sym(1)   == Sym(1) == Sym(1)^true

    ## issue with `pycall_hasproperty` and nothing values.
    @test !SymPyPythonCall.is_(:rational, Sym(2.5)) ## !SymPy.is_rational(Sym(2.5))

    ## Issue #405 with ambigous methods
    @syms α
    M = Sym[1 2; 3 4] ## XXX M = SymMatrix([1 2; 3 4])
    @test α * M == M * α
    @test 2 * M == M * 2
    @test isa(M/α, SymMatrix)
    @test isa(α * inv(M), SymMatrix)

    ## issue #408 with inv
    @syms n::(integer,positive)
    A = sympy.MatrixSymbol("A", n, n)
    @test A.inv() == A.I ## XXX inv(A) == A.I

    # ceil broken
    @syms x
    @test limit(ceil(x), x=>0, dir="+") != limit(ceil(x), x=>0, dir="-")
    @test limit(floor(x), x=>0, dir="+") != limit(floor(x), x=>0, dir="-")

    ## Issue #433 add sympy docstrings, clean up docstring
    # XX sprint(io -> show(io, SymPy.Doc(:sin)))
    Base.Docs.getdoc(sin(x))

    ## ???
    @test SymPy.convert_expr(sympy.Indexed(sympy.IndexedBase("x"), 1, -2)) == :(x[1, -2]) ## XXX SymPy.convert_expr(sympy.Indexed(sympy.IndexedBase(:x), 1, -2)) == :(x[1, -2])

end
