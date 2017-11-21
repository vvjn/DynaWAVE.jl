using Documenter, DynaWAVE

makedocs(
    format = :html,
    sitename = "DynaWAVE",
    pages = [
             "Tutorial" => "index.md", "Static NA" => "static.md", "Functions" => "funs.md",
             "Internals" => "internals.md"
    ]
)

deploydocs(
    repo   = "github.com/vvjn/DynaWAVE.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    julia = "0.6"
)
