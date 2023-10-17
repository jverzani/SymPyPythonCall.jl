# deprecations
# exported in 1.0.3
# [Symbol("@syms"), :Differential, :E, :Eq, :FALSE, :Ge, :Gt, :IM, :Le, :Lt, :N, :Ne, :PI, :Sym, :SymPyPythonCall, :TRUE, :VectorField, :Wild, :apart, :ask, :cancel, :degree, :doit, :dsolve, :elements, :expand, :expand_trig, :factor, :free_symbols, :hessian, :integrate, :lambdify, :lhs, :limit, :linsolve, :nonlinsolve, :nroots, :nsolve, :oo, :plot_implicit, :plot_parametric_surface, :real_roots, :refine, :rewrite, :rhs, :roots, :series, :simplify, :solve, :solveset, :subs, :summation, :symbols, :sympy, :sympy_plotting, :together, :zoo, :¬, :∧, :∨, :≦, :≧, :≪, :≫, :≶, :≷, :⩵, :𝑄]

Base.@deprecate expand_trig(x) x.expand_trig() true

Base.@deprecate_binding _sympy_ sympy.o
Base.@deprecate_binding Q 𝑄

function unSym(x)
    Base.depwarn("`unSym` is deprecated. Use `↓(x)`.", :unSym)
    ↓(x)
end

function asSymbolic(x)
    Base.depwarn("`asSymbolic` is deprecated. Use `↑(x)`.", :asSymbolic)
    ↑(x)
end

function VectorField(args...; kwargs...)
    Base.depwarn("`VectorField` has been removed", :VectorField)
    nothing
end

function plot_implicit(args...; kwargs...)
    Base.depwarn("`plot_implicit` has been removed. See `sympy.plotting.plot_implicit`", :VectorField)
    nothing
end

function plot_parametric_surface(args...; kwargs...)
    Base.depwarn("`VectorField` has been removed. See `sympy.plotting.plot3d_parametric_surface`", :plot_parametric_surface)
    nothing
end
