# Run in fresh processes against the unchanged checkout and isolated checkout.
# The same synthetic inputs exercise C3/C4 and both nitrogen branches, without AD.
using Agrocosm, SHA, TOML
include(joinpath(@__DIR__, "weather_attribution_fixture.jl"))

function array_fingerprints!(result, value, path)
    if value isa AbstractArray
        result[path] = bytes2hex(sha256(reinterpret(UInt8, vec(Array(value)))))
    elseif !(value isa Number || value isa Symbol || value isa Function)
        for name in fieldnames(typeof(value))
            array_fingerprints!(result, getfield(value, name), "$path.$name")
        end
    end
    return result
end

length(ARGS) == 1 || error("usage: weather_production_fingerprint.jl OUTPUT.toml")
ispath(ARGS[1]) && error("refusing to overwrite fingerprint")
results = Dict{String, Any}()
for T in (Float32, Float64), cft_id in (1, 3), nitrogen in (false, true)
    case = weather_attribution_fixture(cft_id; T)
    driver = cft_id == 1 ? Agrocosm.daily_crop_C3! : Agrocosm.daily_crop_C4!
    driver(first(case.days), case.harvest_day + 2,
        Agrocosm.ProcessModules(case.cft, case.parameters), case.climate, case.state;
        nitrogen_limit_vcmax = nitrogen, crop_resp_fix = true,
        fertilizer = :yes, manure = true, with_tillage = true,
        update_vernalization_requirement = false, reuse_output = true)
    label = "$(T)_cft$(cft_id)_nitrogen$(nitrogen)"
    results[label] = array_fingerprints!(Dict{String, Any}(), case.state, "state")
end
open(ARGS[1], "w") do io
    TOML.print(io, results; sorted = true)
end
println("Production fingerprints written for 8 precision/CFT/nitrogen cases")
