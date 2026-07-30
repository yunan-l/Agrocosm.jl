# Compatibility only for pre-lifecycle unit tests. Production kernels and the
# package module define selectors exclusively for `ModelState`.
Agrocosm.crop_prognostic(x::Agrocosm.Crop) = x.state
Agrocosm.crop_fluxes(x::Agrocosm.Crop) = x.fluxes
Agrocosm.crop_events(x::Agrocosm.Crop) = x.events
Agrocosm.crop_canopy_auxiliary(x::Agrocosm.Crop) = x.auxiliary.canopy
Agrocosm.crop_photosynthesis_auxiliary(x::Agrocosm.Crop) = x.auxiliary.photosynthesis
Agrocosm.crop_stress_auxiliary(x::Agrocosm.Crop) = x.auxiliary.stress
Agrocosm.crop_phenology_auxiliary(x::Agrocosm.Crop) = x.auxiliary.phenology
Agrocosm.crop_phenology_input(x::Agrocosm.Crop) = x.auxiliary.phenology
Agrocosm.crop_calendar_input(x::Agrocosm.Crop) = x.auxiliary.calendar
Agrocosm.crop_root_auxiliary(x::Agrocosm.Crop) = x.auxiliary.root
Agrocosm.crop_root_input(x::Agrocosm.Crop) = x.auxiliary.root
Agrocosm.soil_properties(x::Agrocosm.Soil) = x.properties

for (selector, field) in (
    (:soil_water_prognostic, :water),
    (:soil_water_fluxes, :water),
    (:soil_water_auxiliary, :water),
    (:soil_thermal_prognostic, :thermal),
    (:soil_thermal_fluxes, :thermal),
    (:soil_thermal_input, :thermal),
    (:soil_carbon_prognostic, :carbon),
    (:soil_carbon_fluxes, :carbon),
    (:soil_carbon_auxiliary, :carbon),
    (:soil_nitrogen_prognostic, :nitrogen),
    (:soil_nitrogen_fluxes, :nitrogen),
    (:soil_nitrogen_auxiliary, :nitrogen),
    (:soil_decomposition_auxiliary, :decomposition),
    (:soil_decomposition_input, :decomposition),
    (:soil_decomposition_workspace, :decomposition),
    (:soil_management_prognostic, :management),
    (:soil_management_fluxes, :management),
    (:soil_management_input, :management),
    (:soil_surface_litter_prognostic, :surface_litter),
    (:soil_surface_litter_fluxes, :surface_litter),
    (:soil_surface_litter_auxiliary, :surface_litter),
    (:soil_snow_prognostic, :snow),
    (:soil_snow_fluxes, :snow),
)
    @eval Agrocosm $selector(x::Soil) = getfield(x, $(QuoteNode(field)))
end
