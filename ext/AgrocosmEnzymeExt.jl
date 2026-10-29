module AgrocosmEnzymeExt

import Agrocosm
import Enzyme

function __init__()
    # Float64 CFT tail fields extend past Enzyme's default 512-byte type offset.
    Enzyme.API.maxtypeoffset!(1024)
    return nothing
end

include("AgrocosmEnzymeExt/core.jl")
include("AgrocosmEnzymeExt/physiology.jl")
include("AgrocosmEnzymeExt/soil_water.jl")
include("AgrocosmEnzymeExt/continuous_transition.jl")
include("AgrocosmEnzymeExt/parameter_rules.jl")
include("AgrocosmEnzymeExt/seasonal_loss.jl")
include("AgrocosmEnzymeExt/seasonal_soil_loss.jl")
include("AgrocosmEnzymeExt/seasonal_gradients.jl")
include("AgrocosmEnzymeExt/management_adaptation.jl")
include("AgrocosmEnzymeExt/yield_sensitivity.jl")

end
