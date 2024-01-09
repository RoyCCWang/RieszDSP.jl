

# if Base.HOME_PROJECT[] !== nothing
#     # JuliaLang/julia/pull/28625
#     Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
# end

using Documenter, RieszDSP

# Generate examples
#include("generate.jl")

makedocs(
    sitename = "RieszDSP",
    modules = [RieszDSP],
    #format = Documenter.HTML(),
    pages = [
        "Overview" => "index.md",
        "Public API" => "api.md",
        "Demo: image analysis" => "demo.md",
    ],
)

# makedocs(
#     sitename="RieszDSP.jl",
#     modules=[RieszDSP],
#     format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing)=="true"),
#     pages=[
#         "Overview" => "index.md",
#         "Public API" => "api.md",
#         "Demo: code walk-through" => "demo_code.md",
#         "Demo: return variables" => "return_vars.md",
#     ],
#     #strict=true,
# )
# deploydocs(
#     repo = "github.com/RoyCCWang/RieszDSP.jl.git",
#     target = "build",
#     branch = "gh-pages",
#     versions = ["stable" => "v^", "v#.#" ],
# )