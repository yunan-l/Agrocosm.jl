# Experimental CPU-MPI building block. Every member must share its start and
# harvest calendar. Global heterogeneous-event scheduling is a separate gate.
function weather_cohort_replay(forcing, initial_state, cft, parameters, climate, days, harvest_day;
    irrigation::Bool=false)
    T=eltype(forcing)
    cells=size(forcing,2)
    ndims(forcing)==3 && size(forcing,3)==5 && cells>0 || error("invalid cohort forcing")
    length(Agrocosm.crop_prognostic(initial_state).canopy.lai)==cells || error("cohort cell mismatch")
    all(!iszero,Array(Agrocosm.crop_prognostic(initial_state).phenology.is_growing)) || error("cohort must already be established")
    first(days)>=1 && last(days)+1==harvest_day<=size(forcing,1) || error("invalid cohort calendar")
    size(climate.temp)==size(forcing)[1:2] || error("cohort climate mismatch")
    all(isfinite,forcing) || error("cohort forcing must be finite")
    all(variable->all(>=(zero(T)),view(forcing,:,:,variable)),(2,3,5)) || error("negative cohort forcing")
    state=deepcopy(initial_state)
    controlled=_controlled_climate(climate,forcing)
    driver=cft.path==1 ? Agrocosm.daily_crop_C3! : cft.path==2 ? Agrocosm.daily_crop_C4! : error("unsupported crop pathway")
    processes=Agrocosm.ProcessModules(cft,parameters)
    counts=zeros(Int,cells); dates=zeros(Int,cells)
    yields=fill(T(NaN),cells); failed=falses(cells)
    for day in first(days):harvest_day
        driver(day,day,processes,controlled,state;fertilizer=:yes,manure=true,with_tillage=true,
            irrigation,nitrogen_limit_vcmax=true,crop_resp_fix=true,
            update_vernalization_requirement=false,reuse_output=true)
        # Host reads are an output-validation boundary, outside the AD path.
        events=vec(Array(state.output.calendar.harvest_event))
        carbon=Array(Agrocosm.crop_fluxes(state).carbon.yield)
        failure=Array(Agrocosm.crop_events(state).harvest)
        for cell in 1:cells
            if events[cell]!=0
                counts[cell]+=1; dates[cell]=day
                yields[cell]=carbon[cell]/T(0.45)*T(0.01)
                failed[cell]=failure[cell]!=0
            end
        end
    end
    matches=all(==(1),counts) && all(==(harvest_day),dates) && !any(failed)
    return (;yields,harvest_dates=dates,harvest_counts=counts,failed,schedule_matches=matches)
end

function _weather_cohort_block(forcing,state,cft,parameters,climate,days,depth,terminal,pathway,irrigation)
    fill!(Agrocosm.crop_events(state).sowing,0)
    for day in days
        _enzyme_continuous_transition!(state,cft,parameters,climate,day,:gpp,depth,
            irrigation,true,true,true,forcing,pathway,Val(:production))
    end
    T=eltype(forcing)
    return terminal ? sum(Agrocosm.crop_prognostic(state).carbon.storage)/T(0.45)*T(0.01) : zero(T)
end

function weather_cohort_gradient(forcing::Array{T,3},initial_state,cft,parameters,climate,days,harvest_day;
    block_days::Int=7,irrigation::Bool=false) where {T<:AbstractFloat}
    block_days>0 || error("block days must be positive")
    reference=weather_cohort_replay(forcing,initial_state,cft,parameters,climate,days,harvest_day;irrigation)
    reference.schedule_matches || error("cohort members do not share one valid harvest calendar")
    state=deepcopy(initial_state)
    Agrocosm.enzyme_prepare_daily_state!(state)
    depth=Tuple(state.inputs.soil.properties.layer_depth)
    pathway=cft.path==1 ? Val(:C3) : Val(:C4)
    ranges=[day:min(day+block_days-1,last(days)) for day in first(days):block_days:last(days)]
    snapshots=Vector{typeof(state)}(undef,length(ranges))
    primal=zero(T)
    for index in eachindex(ranges)
        snapshots[index]=deepcopy(state)
        primal+=_weather_cohort_block(forcing,state,cft,parameters,climate,ranges[index],depth,index==length(ranges),pathway,irrigation)
    end
    isapprox(primal,sum(reference.yields);rtol=1e-10,atol=1e-10) || error("cohort primal mismatch")
    gradient=zeros(T,size(forcing))
    shadow=Agrocosm.enzyme_zero_tangent(state)
    reverse_primal=zero(T)
    for index in reverse(eachindex(ranges))
        block_state=deepcopy(snapshots[index])
        result=Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),_weather_cohort_block,
            Enzyme.Duplicated(forcing,gradient),Enzyme.Duplicated(block_state,shadow),
            Enzyme.Const(cft),Enzyme.Const(parameters),Enzyme.Const(climate),Enzyme.Const(ranges[index]),
            Enzyme.Const(depth),Enzyme.Const(index==length(ranges)),Enzyme.Const(pathway),Enzyme.Const(irrigation))
        reverse_primal+=result[2]
    end
    all(isfinite,gradient) || error("nonfinite cohort gradient")
    return (;primal,reverse_primal,gradient,yields=reference.yields,
        harvest_day,attribution_scope=:common_calendar_cohort)
end
