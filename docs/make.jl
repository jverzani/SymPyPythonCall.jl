ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"


using SymPyPythonCall
using Documenter

makedocs(
    sitename = "SymPyPythonCall",
    format = Documenter.HTML(),
    modules = [SymPyPythonCall]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/jverzani/SymPyPythonCall.jl.git"
)


#DocMeta.setdocmeta!(SymPyPythonCall, :DocTestSetup, :(using SymPyPythonCall); recursive=true)

# makedocs(;
#     modules=[SymPyPythonCall],
#     authors="jverzani <jverzani@gmail.com> and contributors",
#     repo="https://github.com/jverzani/SymPyPythonCall.jl/blob/{commit}{path}#{line}",
#     sitename="SymPyPythonCall.jl",
#     format=Documenter.HTML(;
#         prettyurls=get(ENV, "CI", "false") == "true",
#         canonical="https://jverzani.github.io/SymPyPythonCall.jl",
#         assets=String[],
#     ),
#     pages=[
#         "Home" => "index.md",
#         "Examples" => "introduction.md"
#     ],
# )

# # deploydocs(;
#     repo="github.com/jverzani/SymPyPythonCall.jl",
#     devbranch="main",
# )
