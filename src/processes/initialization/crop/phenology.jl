"""Integrated phenology coordinates and discrete prognostic modes."""
mutable struct CropPhenology{A, B, I}
    vdsum::A               # Accumulated effective vernalization days (day equivalent).
    husum::A               # Accumulated heat units since cultivation (°C day).
    senescence::B          # Current phenological senescence-mode flag.
    senescence_previous::B # Previous-day senescence-mode flag used at transitions.
    harvesting::B          # Current phenological harvest-readiness flag.
    harvesting_previous::B # Previous-day harvest-readiness flag used to detect harvest.
    growing_days::I        # Number of simulated days since cultivation (day).
    is_growing::I          # Active crop-presence/growth mode (0/1).
    # Surviving fraction of potential grain set (0-1), reduced irreversibly by
    # heat exposure during flowering and reset to one at sowing. Prognostic
    # rather than diagnostic precisely so the damage cannot be undone by a cool
    # spell after the event.
    grain_set_fraction::A
    # Surviving fraction of potential grain FILLING (0-1), reduced irreversibly
    # by heat exposure after anthesis and reset to one at sowing. Separate from
    # `grain_set_fraction` because the two damage different things at different
    # times: grain set fixes the NUMBER of grains around flowering, filling
    # fixes their WEIGHT afterwards, and a model carrying only the first cannot
    # express terminal heat at all. Prognostic for the same reason as its
    # sibling - starch not deposited is not deposited later.
    grain_fill_fraction::A
    # Accumulated rainfall EXCESS above this crop's heavy-rain threshold over the
    # season (mm), reset to zero at sowing, and converted into a recovery loss at
    # harvest once it passes the crop's tolerance.
    #
    # An accumulator with a tolerance rather than a daily fraction, because the
    # first version charged for every heavy day and so took 54% of rice yield in
    # an AVERAGE season - the term was pricing a normal climate, not a wet year.
    # A field tolerates ordinary heavy rain; damage begins when the season's
    # accumulation exceeds what the drainage, the standing crop and the harvest
    # window can absorb. That is the same structure as the FAO-56 supply plateau,
    # which is the one mechanism in this project that lowered variance and raised
    # correlation at once. A third pathway, distinct from both siblings: grain set fixes the
    # number of grains and grain filling their weight, and neither can express a
    # crop that grew normally and was then lodged, sprouted, diseased or left
    # unharvestable in a wet field.
    #
    # It reads the FORCING rather than a soil state deliberately. The documented
    # wet-year pathways - lodging, sprouting, grain disease, harvest and field
    # access loss - act on the canopy and the grain, not through root-zone
    # anoxia, so they need no saturated soil. That matters here because this
    # model cannot produce a saturated root zone at all: under 120 mm/day for
    # five days a 75%-clay soil reaches 23% of its gravitational pore space,
    # since rejected rainfall leaves instantly as surface runoff rather than
    # ponding. See `docs/20`.
    heavy_rain_excess::A
end

"""Static and current-day algebraically derived phenology variables."""
mutable struct CropPhenologyAuxiliary{A, B}
    phu::A         # Potential heat units required for maturity (°C day).
    winter_type::B # Winter-crop/vernalization requirement flag.
    fphu::A        # Current heat-unit fraction derived from `husum / phu` (0–1).
end

init_crop_phenology(cell_size::Int, device) = init_crop_phenology(Float32, cell_size, device)
function init_crop_phenology(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))
    bool_state(value = false) = device(fill(value, cell_size))

    return CropPhenology(
        float_state(),
        float_state(),
        bool_state(),
        bool_state(),
        bool_state(true),
        bool_state(true),
        device(zeros(Int32, cell_size)),
        device(zeros(Int32, cell_size)),
        device(ones(T, cell_size)),
        device(ones(T, cell_size)),
        device(zeros(T, cell_size)),
    )
end

function init_crop_phenology_auxiliary(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    return CropPhenologyAuxiliary(
        device(zeros(T, cell_size)),
        device(fill(false, cell_size)),
        device(zeros(T, cell_size)),
    )
end
