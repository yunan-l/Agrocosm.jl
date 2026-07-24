using Agrocosm
using Documenter

DocMeta.setdocmeta!(Agrocosm, :DocTestSetup, :(using Agrocosm); recursive = true)

makedocs(
    modules = [Agrocosm],
    sitename = "Agrocosm.jl",
    authors = "Yunan Lin and contributors",
    format = Documenter.HTML(
        canonical = "https://yunan-l.github.io/Agrocosm.jl",
        edit_link = "main",
        prettyurls = get(ENV, "CI", "false") == "true",
        size_threshold = 500_000,
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "API reference" => "api.md",
    ],
    warnonly = [:missing_docs],
)

deploydocs(
    repo = "github.com/yunan-l/Agrocosm.jl.git",
    devbranch = "main",
    push_preview = true,
)
