using Documenter, DynaWAVE

makedocs(
    format = :html,
    sitename = "DynaWAVE",
    pages = [
        "index.md", "Functions" => "funs.md"
    ]
)

deploydocs(
    repo   = "github.com/vvjn/DynaWAVE.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
    julia = "0.6"
)