using Agrocosm
using Enzyme
using Test

const _CFT_OFFSET_EXTENSION = Base.get_extension(Agrocosm, :AgrocosmEnzymeExt)

function _cft_offset_square(theta, base_cft, parameter_names)
    cft = _CFT_OFFSET_EXTENSION._replace_cft_parameters(
        base_cft, theta, parameter_names,
    )
    value = getfield(cft, first(parameter_names))
    return value * value
end

@testset "Enzyme CFT tail-field gradients" begin
    cft = convert_precision(Float64, cft1)
    for name in (:hiopt, :himin)
        index = findfirst(==(name), fieldnames(typeof(cft)))
        @test fieldoffset(typeof(cft), index) >= 512

        value = getfield(cft, name)
        theta = [value]
        gradient = [0.0]
        Enzyme.autodiff(
            Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
            _cft_offset_square,
            Enzyme.Active,
            Enzyme.Duplicated(theta, gradient),
            Enzyme.Const(cft),
            Enzyme.Const((name,)),
        )
        @test gradient[1] ≈ 2value rtol = 1e-8
    end
end
