using Agrocosm
using Enzyme
using Test

# Regression guard for a SILENT reverse-mode zero.
#
# Enzyme's type analysis walks a struct only up to a byte-offset limit, and a
# field past that limit differentiates to an exact 0.0 in reverse mode with no
# error and a correct primal. `CFTParameters{Float64, Int32}` is 640 bytes, and
# the truncated tail is exactly this project's new physics - `leaf_dimension`,
# `leaf_emissivity`, `flowering_start`, `flowering_end`,
# `sterility_temperature`, `sterility_rate` - while `b`, `gmin` and `knstore`,
# the only parameters the rest of the AD suite differentiates, sit below the cut
# and stay correct in the same sweep. That is why nothing here caught it.
#
# `AgrocosmEnzymeExt.__init__` raises the limit. This file exists so that if the
# call is ever removed, the failure is a red test rather than a sink whose
# calibration gradient is structurally zero and reads as a converged
# insensitivity.
#
# The assertions below are deliberately about the FIELD LAYOUT rather than about
# a fixed list of names: a field added to the end of `CFTParameters` inherits
# the same hazard, so the test has to find it on its own.

const _EXT = Base.get_extension(Agrocosm, :AgrocosmEnzymeExt)

"""Reverse-mode d/dθ of `θ -> field(θ)^2` routed through the CFT replacement."""
function _replacement_gradient(cft, name::Symbol, value::Float64)
    objective(theta, base, names) = begin
        replaced = _EXT._replace_cft_parameters(base, theta, names)
        v = getfield(replaced, first(names))
        return v * v
    end
    theta = [value]
    dtheta = [0.0]
    Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal), objective,
        Enzyme.Active, Enzyme.Duplicated(theta, dtheta),
        Enzyme.Const(cft), Enzyme.Const((name,)),
    )
    return dtheta[1]
end

"""Scalar `T` fields of `CFTParameters{Float64, Int32}` with their byte offsets."""
function _scalar_fields()
    type = Agrocosm.CFTParameters{Float64, Int32}
    return [
        (name, Int(fieldoffset(type, index)))
        for (index, name) in enumerate(fieldnames(type))
        if fieldtype(type, index) === Float64
    ]
end

@testset "Every scalar CFT parameter is reverse-differentiable" begin
    cft = convert_precision(Float64, Agrocosm.cft1)
    fields = _scalar_fields()
    @test !isempty(fields)

    # The test is only meaningful if the struct actually extends past the
    # default offset limit. If a future layout change brings every field under
    # it, this fails and says so rather than passing vacuously.
    beyond = [(name, offset) for (name, offset) in fields if offset >= 512]
    @test !isempty(beyond)

    for (name, offset) in fields
        value = Float64(getfield(cft, name))
        # A field whose value is zero has an analytic derivative of zero too,
        # so it cannot distinguish a working gradient from a truncated one.
        # `sterility_rate` is exactly this case: its CFT default is 0.0.
        probe = iszero(value) ? 0.25 : value
        gradient = _replacement_gradient(cft, name, probe)
        @test isapprox(gradient, 2 * probe; rtol = 1e-8) ||
              error("reverse gradient for $name at byte offset $offset is " *
                    "$gradient, expected $(2 * probe); if it is exactly zero, " *
                    "AgrocosmEnzymeExt.__init__ is no longer raising " *
                    "Enzyme.API.maxtypeoffset!")
    end
end

@testset "The sink's own parameters carry a gradient" begin
    # Named explicitly, on top of the layout sweep above, because these are the
    # ones the reproductive sink is calibrated through. A zero here is the
    # failure that would look like a converged fit.
    cft = convert_precision(Float64, Agrocosm.cft1)
    for (name, probe) in ((:sterility_rate, 0.01), (:sterility_temperature, 35.0),
                          (:flowering_start, 0.45), (:flowering_end, 0.70),
                          (:leaf_dimension, 0.08), (:leaf_emissivity, 0.97))
        gradient = _replacement_gradient(cft, name, probe)
        @test gradient ≈ 2 * probe rtol = 1e-8
        @test gradient != 0
    end
end
