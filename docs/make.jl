ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"


using SymPyCall
using Documenter

makedocs(sitename="My Documentation")


#DocMeta.setdocmeta!(SymPyCall, :DocTestSetup, :(using SymPyCall); recursive=true)

# makedocs(;
#     modules=[SymPyCall],
#     authors="jverzani <jverzani@gmail.com> and contributors",
#     repo="https://github.com/jverzani/SymPyCall.jl/blob/{commit}{path}#{line}",
#     sitename="SymPyCall.jl",
#     format=Documenter.HTML(;
#         prettyurls=get(ENV, "CI", "false") == "true",
#         canonical="https://jverzani.github.io/SymPyCall.jl",
#         assets=String[],
#     ),
#     pages=[
#         "Home" => "index.md",
#         "Examples" => "introduction.md"
#     ],
# )

# # deploydocs(;
#     repo="github.com/jverzani/SymPyCall.jl",
#     devbranch="main",
# )
