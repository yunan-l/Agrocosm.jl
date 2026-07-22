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
        "Model concepts" => [
            "Overview" => "concepts/overview.md",
            "State lifecycle" => "concepts/state_lifecycle.md",
            "Daily process order" => "concepts/daily_processes.md",
        ],
        "Using Agrocosm" => [
            "Inputs and outputs" => "guide/inputs_outputs.md",
            "CPU, GPU, and precision" => "guide/backends.md",
            "Checkpoints" => "guide/checkpoints.md",
        ],
        "Science and development" => [
            "Validation and limitations" => "science/validation.md",
            "Roadmap" => "development/roadmap.md",
        ],
        "API reference" => "api.md",
    ],
    warnonly = [:missing_docs],
)

deploydocs(
    repo = "github.com/yunan-l/Agrocosm.jl.git",
    devbranch = "main",
    push_preview = true,
)
