"""
PftParameters{T,S}

Plant functional type parameter bundle for one crop type.
Contains phenology, photosynthesis, allocation, and nutrient traits.
"""
struct Temp{T} # lower and upper coldest monthly mean temperature(deg C)
    low::T  # Lower coldest-month temperature limit (°C).
    high::T # Upper coldest-month temperature limit (°C).
end

struct TempCO2{T}  # lower and upper temperature limit for co2 (deg C)
    low::T  # Lower temperature limit of CO₂ response (°C).
    high::T # Upper temperature limit of CO₂ response (°C).
end

struct TempPhotos{T} # lower and upper limit of temperature optimum for photosynthesis(deg C)
    low::T  # Lower optimum temperature for photosynthesis (°C).
    high::T # Upper optimum temperature for photosynthesis (°C).
end

struct TvEff{T} # min & max tv: lower and upper temperature threshold under which vernalization is possible (deg C)
    low::T  # Lower temperature allowing effective vernalization (°C).
    high::T # Upper temperature allowing effective vernalization (°C).
end

struct TvOpt{T}  # min & max tv: lower and upper temperature threshold under which vernalization is optimal (deg C)
    low::T  # Lower optimum vernalization temperature (°C).
    high::T # Upper optimum vernalization temperature (°C).
end

struct BaseTemp{T} # min & max basetemp: base temperature
    low::T  # Base temperature before flowering (°C).
    high::T # Base temperature after flowering (°C).
end

struct nc_ratio{T} # N:C mass ratio
    root::T # Root nitrogen-to-carbon mass ratio (gN gC⁻¹).
    sto::T  # Storage-organ nitrogen-to-carbon mass ratio (gN gC⁻¹).
    pool::T # Mobile-pool nitrogen-to-carbon mass ratio (gN gC⁻¹).
end

struct ratio{T} # relative C:N ratios
    root::T # Root-to-leaf relative C:N ratio.
    sto::T  # Storage-organ-to-leaf relative C:N ratio.
    pool::T # Mobile-pool-to-leaf relative C:N ratio.
end

struct ncleaf{T} # relative C:N ratios
    low::T    # Minimum leaf nitrogen-to-carbon mass ratio (gN gC⁻¹).
    median::T # Reference leaf nitrogen-to-carbon mass ratio (gN gC⁻¹).
    high::T   # Maximum leaf nitrogen-to-carbon mass ratio (gN gC⁻¹).
end

struct K_Litter10{T}
    leaf::T # annual turnover rate at 10 °C (yr⁻¹)
    root::T # annual below-ground turnover rate at 10 °C (yr⁻¹)
end

struct NuptakeKinetics{T}
    vmax::T # Maximum uptake per unit fine-root carbon (gN kgC-1 day-1)
    kmin::T # Uptake term independent of Michaelis-Menten saturation
    Km::T   # Half-saturation concentration (gN m-3)
end

_convert_precision(::Type{T}, value::AbstractFloat) where {T <: AbstractFloat} = T(value)
_convert_precision(::Type{T}, value::Integer) where {T <: AbstractFloat} = value
_convert_precision(::Type{T}, value::Temp) where {T <: AbstractFloat} = Temp{T}(T(value.low), T(value.high))
_convert_precision(::Type{T}, value::TempCO2) where {T <: AbstractFloat} = TempCO2{T}(T(value.low), T(value.high))
_convert_precision(::Type{T}, value::TempPhotos) where {T <: AbstractFloat} = TempPhotos{T}(T(value.low), T(value.high))
_convert_precision(::Type{T}, value::TvEff) where {T <: AbstractFloat} = TvEff{T}(T(value.low), T(value.high))
_convert_precision(::Type{T}, value::TvOpt) where {T <: AbstractFloat} = TvOpt{T}(T(value.low), T(value.high))
_convert_precision(::Type{T}, value::BaseTemp) where {T <: AbstractFloat} = BaseTemp{T}(T(value.low), T(value.high))
_convert_precision(::Type{T}, value::nc_ratio) where {T <: AbstractFloat} = nc_ratio{T}(T(value.root), T(value.sto), T(value.pool))
_convert_precision(::Type{T}, value::ratio) where {T <: AbstractFloat} = ratio{T}(T(value.root), T(value.sto), T(value.pool))
_convert_precision(::Type{T}, value::ncleaf) where {T <: AbstractFloat} = ncleaf{T}(T(value.low), T(value.median), T(value.high))
_convert_precision(::Type{T}, value::K_Litter10) where {T <: AbstractFloat} = K_Litter10{T}(T(value.leaf), T(value.root))
_convert_precision(::Type{T}, value::NuptakeKinetics) where {T <: AbstractFloat} =
    NuptakeKinetics{T}(T(value.vmax), T(value.kmin), T(value.Km))

@kwdef struct PftParameters{T <: AbstractFloat, S <: Integer}
    name::S                 # Numeric crop/PFT identifier.
    plant_type::S           # LPJmL plant-type category identifier.
    path::S                 # Photosynthetic pathway: 1 = C3, 2 = C4.
    temp::Temp{T}           # Bioclimatic cold-temperature limits (°C).
    temp_co2::TempCO2{T}    # Temperature limits of CO₂ response (°C).
    temp_photos::TempPhotos{T} # Optimum photosynthesis temperature interval (°C).
    tv_eff::TvEff{T}        # Effective vernalization temperature interval (°C).
    tv_opt::TvOpt{T}        # Optimum vernalization temperature interval (°C).
    psens::T                # Photoperiod sensitivity coefficient.
    pb::T                   # Lower/short-day photoperiod threshold (h).
    ps::T                   # Upper/long-day photoperiod threshold (h).
    basetemp::BaseTemp{T}   # Development base temperatures (°C).
    fphuc::T                # Heat-unit fraction ending initial LAI phase.
    flaimaxc::T             # LAI fraction at end of initial LAI phase.
    fphuk::T                # Heat-unit fraction at LAI sigmoid inflection.
    flaimaxk::T             # LAI fraction at sigmoid inflection.
    fphusen::T              # Heat-unit fraction at onset of senescence.
    flaimaxharvest::T       # Fraction of maximum LAI retained at harvest.
    laimax::T               # Maximum potential leaf-area index (m² m⁻²).
    laimin::T               # Minimum/reference LAI trait used by allocation (m² m⁻²).
    hlimit::S               # Maximum crop-cycle duration (days).
    pvd_max::S              # Maximum vernalizing-day requirement (days).
    b::T                    # Leaf maintenance respiration as fraction of `vcmax`.
    albedo_leaf::T          # Leaf shortwave albedo (0–1).
    albedo_litter::T        # Surface-litter shortwave albedo (0–1).
    alphaa::T               # Canopy allometry/allocation coefficient.
    lightextcoeff::T        # Lambert–Beer canopy light-extinction coefficient.
    longevity::T            # Characteristic leaf longevity (years).
    sla::T                  # Specific leaf area (m² leaf gC⁻¹).
    respcoeff::T            # Root/storage maintenance-respiration coefficient.
    shapesenescencenorm::T   # Normalized senescence curve shape parameter.
    fpc::T                   # Foliage projective cover scaling (0–1).
    nc_ratio::nc_ratio{T}   # Organ nitrogen-to-carbon ratios.
    ratio::ratio{T}         # Relative organ C:N ratios used in N partitioning.
    ncleaf::ncleaf{T}       # Minimum/reference/maximum leaf N:C ratios.
    k_litter10::K_Litter10{T} # Litter turnover rates at 10 °C.
    beta_root::T            # Exponential root-depth distribution parameter.
    intc::T                 # Canopy interception storage parameter.
    emax::T                 # Maximum transpiration/conductance scaling parameter.
    gmin::T                 # Minimum canopy conductance (mm s⁻¹).
    knstore::T              # Fraction of N demand allocated to storage reserve.
    no3_uptake::NuptakeKinetics{T} # Root NO₃ uptake kinetic parameters.
    nh4_uptake::NuptakeKinetics{T} # Root NH₄ uptake kinetic parameters.
    hiopt::T                # Optimal harvest index.
    himin::T                # Minimum harvest index under stress.
end

"""Return a PFT parameter set whose floating fields consistently use `T`."""
function convert_precision(::Type{T}, pft::PftParameters{<:AbstractFloat, S}) where {T <: AbstractFloat, S <: Integer}
    names = fieldnames(typeof(pft))
    values = map(name -> _convert_precision(T, getfield(pft, name)), names)
    kwargs = NamedTuple{names}(values)
    return PftParameters{T, S}(; kwargs...)
end

"""
cft1

Default PFT preset for temperate cereals (wheat-like C3 crop).
"""
cft1 = PftParameters{Float32, Int32}(
    name = 1,
    plant_type = 1,
    path = 1,
    temp = Temp{Float32}(-1000.0, 1000.0),
    temp_co2 = TempCO2{Float32}(0.0, 40.0),
    temp_photos = TempPhotos{Float32}(12.0, 17.0),
    tv_eff = TvEff{Float32}(-4.0, 17.0),
    tv_opt = TvOpt{Float32}(3.0, 10.0),
    psens = 1.0,
    pb = 8.0,
    ps = 20,
    basetemp = BaseTemp{Float32}(0.0, 0.0),
    fphuc = 0.05,
    flaimaxc = 0.05,
    fphuk = 0.45,
    flaimaxk = 0.95,
    fphusen = 0.7,
    flaimaxharvest = 0.0,
    laimax = 5.0,
    laimin = 2.0,
    hlimit = 330,
    pvd_max = 70,
    b = 0.015,
    albedo_leaf = 0.18,
    albedo_litter = 0.06,
    alphaa = 1.0,
    lightextcoeff = 0.5,
    longevity = 0.5,
    sla = 0.0364661,
    respcoeff = 1.0,
    shapesenescencenorm = 2.0,
    fpc = 1.0,
    nc_ratio = nc_ratio{Float32}(1/25.0, 1/25.0, 1/25.0),
    ratio = ratio{Float32}(1.16, 0.99, 3),
    ncleaf = ncleaf{Float32}(1/58.8, 1/25.0, 1/14.3),
    k_litter10 = K_Litter10{Float32}(0.97, 0.3),
    beta_root = 0.969,
    intc = 0.01,
    emax = 8.0,
    gmin = 1.0,
    knstore = 0.1,
    no3_uptake = NuptakeKinetics{Float32}(1.5, 0.05, 0.70),
    nh4_uptake = NuptakeKinetics{Float32}(4.7, 0.05, 0.45),
    hiopt = 0.60,
    himin = 0.20
)

"""
cft2

Default PFT preset for rice (C3).
"""
cft2 = PftParameters{Float32, Int32}(
    name = 2,
    plant_type = 1,
    path = 1,
    temp = Temp{Float32}(-1000.0, 1000.0),
    temp_co2 = TempCO2{Float32}(6.0, 55.0),
    temp_photos = TempPhotos{Float32}(20.0, 45.0),
    tv_eff = TvEff{Float32}(1000.0, 1000.0),
    tv_opt = TvOpt{Float32}(1000.0, 1000.0),
    psens = 1.0,
    pb = 24,
    ps = 0,
    basetemp = BaseTemp{Float32}(10.0, 10.0),
    fphuc = 0.1,
    flaimaxc = 0.05,
    fphuk = 0.5,
    flaimaxk = 0.95,
    fphusen = 0.8,
    flaimaxharvest = 0.0,
    laimax = 5.0,
    laimin = 5.0,
    hlimit = 180,
    pvd_max = 0,
    b = 0.015,
    albedo_leaf = 0.18,
    albedo_litter = 0.06,
    alphaa = 1.0,
    lightextcoeff = 0.5,
    longevity = 0.33,
    sla = 0.0430598,
    respcoeff = 1.0,
    shapesenescencenorm = 2.0,
    fpc = 1.0,
    nc_ratio = nc_ratio{Float32}(1/25.0, 1/25.0, 1/25.0),
    ratio = ratio{Float32}(1.16, 1.3, 3),
    ncleaf = ncleaf{Float32}(1/58.8, 1/25.0, 1/14.3),
    k_litter10 = K_Litter10{Float32}(0.97, 0.3),
    beta_root = 0.969,
    intc = 0.01,
    emax = 8.0,
    gmin = 1.0,
    knstore = 0.1,
    no3_uptake = NuptakeKinetics{Float32}(1.5, 0.05, 0.70),
    nh4_uptake = NuptakeKinetics{Float32}(4.7, 0.05, 0.45),
    hiopt = 0.60,
    himin = 0.25
)

"""
cft3

Default PFT preset for maize (C4).
"""
cft3 = PftParameters{Float32, Int32}(
    name = 3,
    plant_type = 1,
    path = 2, # C4
    temp = Temp{Float32}(-1000.0, 1000.0),
    temp_co2 = TempCO2{Float32}(8.0, 42.0),
    temp_photos = TempPhotos{Float32}(21.0, 26.0),
    tv_eff = TvEff{Float32}(1000.0, 1000.0),
    tv_opt = TvOpt{Float32}(1000.0, 1000.0),
    psens = 1.0,
    pb = 0,
    ps = 24,
    basetemp = BaseTemp{Float32}(5.0, 15.0),
    fphuc = 0.1,
    flaimaxc = 0.05,
    fphuk = 0.5,
    flaimaxk = 0.95,
    fphusen = 0.75,
    flaimaxharvest = 0.0,
    laimax = 5.0,
    laimin = 4.0,
    hlimit = 240,
    pvd_max = 0,
    b = 0.035,
    albedo_leaf = 0.18,
    albedo_litter = 0.06,
    alphaa = 1.0,
    lightextcoeff = 0.5,
    longevity = 0.33,
    sla = 0.0430598,
    respcoeff = 1.0,
    shapesenescencenorm = 2.0,
    fpc = 1.0,
    nc_ratio = nc_ratio{Float32}(1/25.0, 1/25.0, 1/25.0),
    ratio = ratio{Float32}(1.16, 0.83, 3),
    ncleaf = ncleaf{Float32}(1/58.8, 1/25.0, 1/14.3),
    k_litter10 = K_Litter10{Float32}(0.97, 0.3),
    beta_root = 0.969,
    intc = 0.01,
    emax = 8.0,
    gmin = 1.2,
    knstore = 0.1,
    no3_uptake = NuptakeKinetics{Float32}(1.5, 0.05, 0.70),
    nh4_uptake = NuptakeKinetics{Float32}(4.7, 0.05, 0.45),
    hiopt = 0.60,
    himin = 0.30
)

"""
cft4

Default PFT preset for soybean (C3).
"""
cft4 = PftParameters{Float32, Int32}(
    name = 4,
    plant_type = 1,
    path = 1, # C3
    temp = Temp{Float32}(-1000.0, 1000.0),
    temp_co2 = TempCO2{Float32}(5.0, 45.0),
    temp_photos = TempPhotos{Float32}(28.0, 32.0),
    tv_eff = TvEff{Float32}(1000.0, 1000.0),
    tv_opt = TvOpt{Float32}(1000.0, 1000.0),
    psens = 1.0,
    pb = 0,
    ps = 24,
    basetemp = BaseTemp{Float32}(10.0, 10.0),
    fphuc = 0.15,
    flaimaxc = 0.05,
    fphuk = 0.5,
    flaimaxk = 0.95,
    fphusen = 0.7,
    flaimaxharvest = 0.0,
    laimax = 5.0,
    laimin = 5.0,
    hlimit = 240,
    pvd_max = 70,
    b = 0.015,
    albedo_leaf = 0.18,
    albedo_litter = 0.06,
    alphaa = 1.0,
    lightextcoeff = 0.5,
    longevity = 0.66,
    sla = 0.0326332,
    respcoeff = 1.0,
    shapesenescencenorm = 0.5,
    fpc = 1.0,
    nc_ratio = nc_ratio{Float32}(1/25.0, 1/25.0, 1/25.0),
    ratio = ratio{Float32}(1.16, 0.42, 3),
    ncleaf = ncleaf{Float32}(1/58.8, 1/25.0, 1/14.3),
    k_litter10 = K_Litter10{Float32}(0.97, 0.3),
    beta_root = 0.969,
    intc = 0.01,
    emax = 8.0,
    gmin = 1.2,
    knstore = 0.1,
    no3_uptake = NuptakeKinetics{Float32}(1.5, 0.05, 0.70),
    nh4_uptake = NuptakeKinetics{Float32}(4.7, 0.05, 0.45),
    hiopt = 0.40,
    himin = 0.10
)
