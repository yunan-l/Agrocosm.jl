"""Static or slowly varying soil physical and chemical properties."""
mutable struct SoilProperties{A, M}
    sand_fraction::M # Soil sand mass fraction by layer (0–1).
    clay_fraction::M # Soil clay mass fraction by layer (0–1).
    ph::A            # Soil pH used by nitrogen transformations.
    anion_exclusion::A # Soil-class-specific anion-excluded porosity fraction.
    nitrification_a::A # Soil-class-specific nitrification moisture parameter a.
    nitrification_b::A # Soil-class-specific nitrification moisture parameter b.
    nitrification_c::A # Soil-class-specific nitrification moisture parameter c.
    nitrification_d::A # Soil-class-specific nitrification moisture parameter d.
    layer_depth::A   # Thickness of each model soil layer (mm).
end

const LPJML_NITRIFICATION_FINE = (a = 0.45, b = 1.27, c = 0.0012, d = 2.84)
const LPJML_NITRIFICATION_SANDY = (a = 0.55, b = 1.7, c = -0.007, d = 3.22)

@inline function _uses_sandy_nitrification_parameters(soilcode::Integer)
    return soilcode == 3 || soilcode == 6 || soilcode == 9 ||
        soilcode == 11 || soilcode == 12
end

@inline function _uses_high_anion_exclusion(soilcode::Integer)
    return _uses_sandy_nitrification_parameters(soilcode) ||
        soilcode == 13 || soilcode == 14
end

@kernel inbounds = true function initialize_soil_class_parameters_kernel!(
    anion_exclusion::AbstractVector{T},
    nitrification_a::AbstractVector{T},
    nitrification_b::AbstractVector{T},
    nitrification_c::AbstractVector{T},
    nitrification_d::AbstractVector{T},
    soilcode::AbstractVector{I},
) where {T <: AbstractFloat, I <: Integer}
    cell = @index(Global)
    parameters = _uses_sandy_nitrification_parameters(soilcode[cell]) ?
        LPJML_NITRIFICATION_SANDY : LPJML_NITRIFICATION_FINE
    anion_exclusion[cell] = _uses_high_anion_exclusion(soilcode[cell]) ?
        T(0.4) : T(0.3)
    nitrification_a[cell] = T(parameters.a)
    nitrification_b[cell] = T(parameters.b)
    nitrification_c[cell] = T(parameters.c)
    nitrification_d[cell] = T(parameters.d)
end

init_soil_properties(cell_size::Int, soildepth, device) =
    init_soil_properties(Float32, cell_size, soildepth, device)
function init_soil_properties(::Type{T}, cell_size::Int, soildepth, device) where {T <: AbstractFloat}
    return SoilProperties(
        device(zeros(T, 1, cell_size)),
        device(zeros(T, 1, cell_size)),
        device(zeros(T, cell_size)),
        device(fill(T(0.3), cell_size)),
        device(fill(T(LPJML_NITRIFICATION_FINE.a), cell_size)),
        device(fill(T(LPJML_NITRIFICATION_FINE.b), cell_size)),
        device(fill(T(LPJML_NITRIFICATION_FINE.c), cell_size)),
        device(fill(T(LPJML_NITRIFICATION_FINE.d), cell_size)),
        device(T.(soildepth)),
    )
end

function initialize_soil_class_parameters!(properties::SoilProperties,
                                           soilcode::AbstractVector)
    length(soilcode) == length(properties.ph) || throw(DimensionMismatch(
        "soilcode must contain one LPJ soil class per cell",
    ))
    launch_1D!(
        initialize_soil_class_parameters_kernel!,
        properties.anion_exclusion,
        properties.nitrification_a,
        properties.nitrification_b,
        properties.nitrification_c,
        properties.nitrification_d,
        soilcode,
    )
    return nothing
end
