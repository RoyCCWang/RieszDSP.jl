
using Documenter, RieszDSP

# # local.
# makedocs(
#     sitename = "RieszDSP",
#     modules = [RieszDSP],
#     #format = Documenter.HTML(),
#     pages = [
#         "Overview" => "index.md",
#         "Public API" => "api.md",
#         "Demo: image analysis" => "demo.md",
#     ],
# )

makedocs(
    sitename="RieszDSP.jl",
    modules=[RieszDSP],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing)=="true"),
    pages = [
        "Overview" => "index.md",
        "Public API" => "api.md",
        "Demo: image analysis" => "demo.md",
    ],
)
deploydocs(
    repo = "github.com/RoyCCWang/RieszDSP.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ],
)