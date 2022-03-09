using SymPyCall
using LinearAlgebra
using Test

## SymPyCall only is adding methods for those found in SpecialFunctions
## so not the special polynomials, etc.
using SpecialFunctions


@testset "SpecialFuns" begin
    @syms a b x y
    @syms n m θ ϕ
    @syms ν::integer


    @test sympy.fresnels(Sym(0)) == 0
    @test sympy.fresnels(Sym(oo)) == Sym(1)/2
    @test diff(sympy.fresnels(x), x) == sin(PI*x^2/2)
    #@test evalf(fresnels(Sym(2)), 30) == Sym(parse(BigFloat, "0.343415678363698242195300815958"))


    @test  diff(sympy.fresnelc(x), x) == cos(PI*x^2/2)


    @test sympy.Ei(Sym(-1)) == sympy.Ei(exp(1im*PI))
    @test diff(sympy.Ei(x), x) == exp(x)/x
    @test diff(sympy.Si(x), x) == sin(x)/x
    @test diff(sympy.Ci(x), x) == cos(x)/x


    @test airyai(Sym(0)) == 3^(Sym(1)/3)/(3*gamma(Sym(2)/3))
    @test airybi(Sym(0)) == 3^(Sym(5)/6)/(3*gamma(Sym(2)/3))


    @test sympy.jacobi(0, a, b, x) == 1
    @test sympy.jacobi(1, a, b, x) == a/2 - b/2 + x*(a/2 + b/2 + 1)
    @test sympy.jacobi(n, 0, 0, x) == sympy.legendre(n, x)
    @test sympy.jacobi(n, a, b, -x) == (-1)^n*sympy.jacobi(n, b, a, x)
    @test conj(sympy.jacobi(n, a, b, x)) == sympy.jacobi(n, conj(a), conj(b), conj(x))


    @test sympy.gegenbauer(0, a, x) == 1
    @test sympy.gegenbauer(1, a, x) == 2a*x
    @test sympy.gegenbauer(2, a, x) == -a + x^2*(2a^2 + 2a)
    @test conj(sympy.gegenbauer(n, a, x)) == sympy.gegenbauer(n, conj(a), conj(x))
    @test diff(sympy.gegenbauer(n, a, x), x) == 2a*sympy.gegenbauer(n-1, a+1, x)


    @test sympy.chebyshevt(0, x) == 1
    @test sympy.chebyshevt(1, x) == x
    @test sympy.chebyshevt(2, x) == 2x^2-1


    @test sympy.chebyshevu(0, x) == 1
    @test sympy.chebyshevu(1, x) == 2x
    @test sympy.chebyshevu(2, x) == 4x^2 - 1
    @test sympy.chebyshevu(n, 0) == cos(π*n/2)
    @test sympy.chebyshevu(n, 1) == n + 1


    @test sympy.legendre(0, x) == 1
    @test sympy.legendre(1, x) == x
    @test sympy.legendre(2, x) == 3x^2/2 - 1//2
    @test sympy.legendre(n, x) == sympy.legendre(n, x)
    @test diff(sympy.legendre(n,x), x) == n*(x*sympy.legendre(n, x) - sympy.legendre(n - 1, x))/(x^2 - 1)


    @test sympy.assoc_legendre(0,0, x) == 1
    @test sympy.assoc_legendre(1,0, x) == x
    @test sympy.assoc_legendre(1,1, x) == -sqrt(-x^2 + 1)
    #@test assoc_legendre(n,m,x) == assoc_legendre(n, m, x)


    @test sympy.hermite(0, x) == 1
    @test sympy.hermite(1, x) == 2x
    @test sympy.hermite(2, x) == 4x^2 - 2
    @test sympy.hermite(n, x) == sympy.hermite(n, x)
    @test diff(sympy.hermite(n,x), x) == 2n*sympy.hermite(n - 1, x)
    @test sympy.hermite(n, -x) == (-1)^n*sympy.hermite(n, x)


    @test sympy.laguerre(0, x) == 1
    @test sympy.laguerre(1, x) == -x + 1
    @test sympy.laguerre(2, x) == x^2/2 - 2x + 1
    @test sympy.laguerre(3, x) == -x^3/6 + 3x^2/2 - 3x + 1
    @test diff(sympy.laguerre(n, x), x) == -sympy.assoc_laguerre(n - 1, 1, x)


    @test sympy.assoc_laguerre(0, a, x) == 1
    @test sympy.assoc_laguerre(1, a, x) == a - x + 1
    @test sympy.assoc_laguerre(2, a, x) == a^2/2 + 3a/2 + x^2/2 + x*(-a - 2) + 1
    @test sympy.assoc_laguerre(n, 0, x) == sympy.laguerre(n, x)
    @test diff(sympy.assoc_laguerre(n, a, x), x) == -sympy.assoc_laguerre(n - 1, a + 1, x)


    @test sympy.Ynm(n, m, θ, ϕ) == sympy.Ynm(n, m, θ, ϕ)
    @test sympy.Ynm(n, -m, θ, ϕ) == (-1)^m*exp(-2im*m*ϕ)*sympy.Ynm(n, m, θ, ϕ)
    @test sympy.Ynm(n, m, -θ, ϕ) == sympy.Ynm(n, m, θ, ϕ)
    @test sympy.Ynm(n, m, θ, -ϕ) == exp(-2im*m*ϕ)*sympy.Ynm(n, m, θ, ϕ)
    @test expand(simplify(sympy.Ynm(0, 0, θ, ϕ)), func=true) == 1/(2*sqrt(PI))

    @test sympy.Ynm_c(n, m, θ, ϕ) == conj(sympy.Ynm(n, m, θ, ϕ))

    @test diff(hankelh1(n, x), x) == hankelh1(n - 1, x)/2 - hankelh1(n + 1, x)/2
    @test diff(hankelh2(n, x), x) == hankelh2(n - 1, x)/2 - hankelh2(n + 1, x)/2


    @test sympy.expand_func(sympy.jn(0, x)) == sin(x)/x
    @test sympy.expand_func(sympy.jn(1, x)) == sin(x)/x^2 - cos(x)/x
    @test rewrite(sympy.jn(ν, x), "besselj") == sqrt(2PI/x)*besselj(ν + Sym(1)/2, x)/2
    VERSION < v"0.6.0-dev" && @test rewrite(sympy.jn(ν, x), "bessely") == (-1)^ν*sqrt(2PI/x)*bessely(-ν - Sym(1)/2, x)/2
    u = N(sympy.jn(2, 5.2+0.3im).evalf(20))
    @test norm(real(u) - 0.099419756723640344491) <= 1e-15 && norm(imag(u) + 0.054525080242173562897) <= 1e-15


    @test sympy.expand_func(sympy.yn(0, x)) == -cos(x)/x
    @test sympy.expand_func(sympy.yn(1, x)) == -cos(x)/x^2-sin(x)/x
    VERSION < v"0.6.0-dev" && @test rewrite(sympy.yn(ν, x), "besselj") == (-1)^(ν + 1)*sqrt(2PI/x)*besselj(-ν - Sym(1)/2, x)/2
    @test rewrite(sympy.yn(ν, x), "bessely") == sqrt(2PI/x)*bessely(ν + Sym(1)/2, x)/2
    @test N(sympy.yn(2, 5.2+0.3im)) ≈ 0.18525034196069722536 + 0.014895573969924817587im


    # gamma, beta and related functions
    @test gamma(Sym(1)) == 1
    @test gamma(Sym(3)/2) == sqrt(PI)/2


    @test diff(polygamma(Sym(0), x), x) == polygamma(Sym(1), x)
    @test diff(polygamma(Sym(0), x), x, 2) == polygamma(Sym(2), x)


    @test diff(beta(x, y), x) == (polygamma(Sym(0), x) - polygamma(Sym(0), x + y)) * beta(x, y)
    @test diff(beta(x, y), y) == (polygamma(Sym(0), y) - polygamma(Sym(0), x + y)) * beta(x, y)


    # test numerical consistency with Julia functions
    @test N(gamma(Sym(4.1))) ≈ gamma(4.1)
    @test N(polygamma(Sym(2), Sym(3.2))) ≈ polygamma(2, 3.2)
    VERSION >= v"0.5.0" && @test N(beta(Sym(1)+1im, Sym(1)+1im)) ≈ beta(1.0+1im, 1.0+1im)


    # Elliptic-type functions

    @test sympy.elliptic_k(Sym(0)) == PI/2
    @test N(sympy.elliptic_k(Sym(1.0 + im))) ≈ 1.50923695405127 + 0.625146415202697im

    @test N(sympy.elliptic_f(Sym(3.0 + im/2), Sym(1.0 + im))) ≈ 2.909449841483 + 1.74720545502474im

    @test sympy.elliptic_e(Sym(0)) == PI/2
    @test N(sympy.elliptic_e(Sym(2.0 - im))) ≈ 0.991052601328069 + 0.81879421395609im

    @test sympy.elliptic_pi(Sym(0), Sym(0)) == PI/2
    @test N(sympy.elliptic_pi(Sym(1.0 - im/3), Sym(2.0 + im))) ≈ 3.29136443417283 + 0.32555634906645im


    # Bessel-type functions
    @test diff(besselj(n, x), x) == (besselj(n - 1, x) - besselj(n + 1, x))/2
    @test rewrite(besselj(n, x), "jn") == sqrt(2x/PI)*sympy.jn(n - 1//2, x)


    @test diff(bessely(n, x), x) == (bessely(n - 1, x) - bessely(n + 1, x))/2
    @test rewrite(bessely(n, x), "yn") == sqrt(2x/PI)*sympy.yn(n - 1//2, x)


    @test diff(besseli(n, x), x) == (besseli(n - 1, x) + besseli(n + 1, x))/2


    @test diff(besselk(n, x), x) == -(besselk(n - 1, x) + besselk(n + 1, x))/2


    # test numerical consistency with Julia functions
    @test N(besselj(Sym(3.2), Sym(1.5))) ≈ besselj(3.2, 1.5)
    @test N(bessely(Sym(3.2), Sym(1.5))) ≈ bessely(3.2, 1.5)
    @test N(besseli(Sym(3.2), Sym(1.5))) ≈ besseli(3.2, 1.5)
    @test N(besselk(Sym(3.2), Sym(1.5))) ≈ besselk(3.2, 1.5)


    @test sympy.expand_func(x*sympy.hyper([1, 1], [2], -x)) == log(x + 1)
end
