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
            "State variables" => "concepts/state_lifecycle.md",
            "Process contracts" => "concepts/process_contracts.md",
            "Daily process order" => "concepts/daily_processes.md",
            "Model processes" => "science/model_equations.md",
            "Crop processes" => "science/crop.md",
            "Soil processes" => "science/soil.md",
            "Climate and surface processes" => "science/climate_surface.md",
            "Numerics and conservation" => "science/numerics.md",
            "Initialization and warm-up" => "science/initialization_warmup.md",
        ],
        "Using Agrocosm" => [
            "Inputs and outputs" => "guide/inputs_outputs.md",
            "Multi-CFT simulations" => "guide/multi_cft.md",
            "HWSD soil initialization" => "guide/hwsd_initialization.md",
            "Global wheat test data" => "guide/global_wheat_subset.md",
            "CPU, GPU, and precision" => "guide/backends.md",
            "Checkpoints" => "guide/checkpoints.md",
        ],
        "Science and development" => [
            "Validation and limitations" => "science/validation.md",
            "Citations" => "science/citations.md",
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
