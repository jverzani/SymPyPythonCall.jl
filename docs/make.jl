ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"


using SymPyPythonCall
using Documenter

makedocs(
    sitename = "SymPyPythonCall",
    format = Documenter.HTML(),
    modules = [SymPyPythonCall],
    warnonly = [:missing_docs],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/jverzani/SymPyPythonCall.jl.git"
)
