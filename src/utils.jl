## --------------------------------------------------
## no Doc trick
Base.Docs.getdoc(u::Sym) = Base.Docs.getdoc(getpy(u))
## --------------------------------------------------

is_(k::Symbol, x::Sym) = is_(k, getpy(x))
function is_(k::Symbol, x::Py)::Bool
    key = "is_$(string(k))"
    val = pygetattr(x, key, false)
    Bool(val)
end
