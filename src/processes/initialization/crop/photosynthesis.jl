"""Current-day crop photosynthetic capacities and limiting factors."""
mutable struct CropPhotosynthesisAuxiliary{A}
    potential_vcmax::A    # Potential Rubisco carboxylation capacity before N limitation (gC m⁻² day⁻¹).
    vcmax::A              # Active Rubisco carboxylation capacity (gC m⁻² day⁻¹).
    nitrogen_limitation::A # Retained fraction of potential Vcmax after N limitation (0–1).
    lambda::A             # Ratio of intercellular to ambient CO₂ partial pressure (0–1).
    temperature_stress::A # Photosynthetic temperature-response multiplier (0–1).
    # The P-model's least-cost optimal ci/ca, used as the UPPER BOUND on the
    # lambda solve rather than as lambda itself. Both matter: LPJmL's solved
    # lambda encodes SOIL WATER limitation - it sits at its 0.85 cap most days and
    # crashes to 0.02 under stress - while the P-model's optimum encodes
    # ATMOSPHERIC demand and moves smoothly with vapour deficit and temperature.
    # Measured on real cells the two correlate only 0.03 to 0.05, so they are not
    # substitutes; bounding the solve by the optimum keeps the water response and
    # adds the demand response it never had. Zero when the P-model is off, and the
    # solver then falls back to its own 0.85.
    pmodel_chi::A         # P-model optimal ci/ca; 0 when the P-model is off.
end

function init_crop_photosynthesis_auxiliary(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_auxiliary() = device(zeros(T, cell_size))
    return CropPhotosynthesisAuxiliary(ntuple(_ -> float_auxiliary(), 6)...)
end
