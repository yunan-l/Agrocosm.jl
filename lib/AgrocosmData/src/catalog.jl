function _resolve_catalog_path(base::AbstractString, path::AbstractString)
    expanded = expanduser(path)
    return isabspath(expanded) ? normpath(expanded) : normpath(joinpath(base, expanded))
end

"""Load dataset paths and PFT metadata from a TOML catalog."""
function load_catalog(path::AbstractString)
    raw = TOML.parsefile(path)
    haskey(raw, "datasets") || throw(ArgumentError("catalog must contain a [datasets] section"))
    haskey(raw, "pfts") || throw(ArgumentError("catalog must contain a [pfts] section"))

    base = dirname(abspath(path))
    specs = Dict{Symbol, DatasetSpec}()
    for (name, entry) in raw["datasets"]
        entry isa AbstractDict || throw(ArgumentError("dataset '$name' must be a TOML table"))
        haskey(entry, "path") || throw(ArgumentError("dataset '$name' is missing path"))
        haskey(entry, "variable") || throw(ArgumentError("dataset '$name' is missing variable"))
        specs[Symbol(name)] = DatasetSpec(
            _resolve_catalog_path(base, entry["path"]),
            entry["variable"];
            units = get(entry, "units", ""),
            pft_ids = get(entry, "pft_ids", Int[]),
            rainfed_bands = get(entry, "rainfed_bands", Int[]),
            irrigated_bands = get(entry, "irrigated_bands", Int[]),
        )
    end

    pfts = raw["pfts"]
    haskey(pfts, "ids") || throw(ArgumentError("[pfts] is missing ids"))
    haskey(pfts, "names") || throw(ArgumentError("[pfts] is missing names"))
    registry = PFTRegistry(pfts["ids"], pfts["names"])
    registry_ids = Set(registry.ids)
    for (name, spec) in specs
        all(id -> id in registry_ids, spec.pft_ids) ||
            throw(ArgumentError("dataset '$name' contains PFT ids absent from the global registry"))
        if spec.management_bands !== nothing
            length(spec.management_bands.rainfed) == length(registry.ids) ||
                throw(ArgumentError("dataset '$name' must map every crop PFT"))
        end
    end
    return DatasetCatalog(specs, registry)
end

function dataset(catalog::DatasetCatalog, name::Symbol)
    haskey(catalog.datasets, name) || throw(KeyError(name))
    return catalog.datasets[name]
end

function pft_index(registry::PFTRegistry, id::Integer)
    index = findfirst(==(Int32(id)), registry.ids)
    isnothing(index) && throw(ArgumentError("unknown PFT id $id"))
    return index
end

function pft_index(registry::PFTRegistry, name::AbstractString)
    index = findfirst(==(String(name)), registry.names)
    isnothing(index) && throw(ArgumentError("unknown PFT name '$name'"))
    return index
end

pft_name(registry::PFTRegistry, id::Integer) = registry.names[pft_index(registry, id)]
