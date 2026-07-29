# Climate buffer updates: monthly aggregation, rolling means, and vernalization metrics.
"""
annual_climbuf!(daily_temp, climbuf, CFT; n=5, kk=0.05)

Update annual climate-buffer diagnostics used by phenology, including rolling
monthly means and vernalization requirements.
"""
function annual_climbuf!(daily_temp::AbstractArray{T},
                         climbuf::ClimBuf,
                         CFT::CFTParameters;
                         n::Int = 5,
                         kk = T(0.05)
) where {T <: AbstractFloat}
    kk = T(kk)
    # Calculate the average temperature for each month.
    # update_monthly
    # length(daily_temp) = 365
    # n = 5, the first n coldest months
    # kk is to rescale the 20-year average monthly temprerature
    
    monthlytemp!(daily_temp, climbuf.mtemp)
    
    # 20-year moving monthly climatology (month, cell).
    launch_2D!(
        climbuf_mtemp20_kernel!,
        climbuf.mtemp20,
        climbuf.mtemp,
        kk,
    )
    n == size(climbuf.min_temp, 1) || throw(ArgumentError(
        "n must match the first dimension of climbuf.min_temp",
    ))
    launch_custom!(
        climbuf_annual_diagnostics_kernel!,
        climbuf.V_req_a,
        length(climbuf.V_req_a),
        climbuf.min_temp,
        climbuf.V_req,
        climbuf.atemp_mean,
        climbuf.mtemp20,
        daily_temp,
        CFT,
        n,
        kk,
    )

end


@kernel inbounds = true function climbuf_mtemp20_kernel!(
                                         climbuf_mtemp20::AbstractArray{T},
                                         climbuf_mtemp::AbstractArray{T},
                                         kk
) where {T <: AbstractFloat}
    
    month, cell = @index(Global, NTuple)
    
    if climbuf_mtemp20[month, cell] < -9998
        climbuf_mtemp20[month, cell] = climbuf_mtemp[month, cell]
    else
        climbuf_mtemp20[month, cell] = (1 - kk) * climbuf_mtemp20[month, cell] + kk * climbuf_mtemp[month, cell]
    end
    
end


@kernel inbounds = true function climbuf_annual_diagnostics_kernel!(
    climbuf_V_req_a::AbstractVector{T},
    climbuf_min_temp::AbstractMatrix{T},
    climbuf_V_req::AbstractVector{T},
    climbuf_atemp_mean::AbstractVector{T},
    climbuf_mtemp20::AbstractMatrix{T},
    daily_temp::AbstractMatrix{T},
    CFT::CFTParameters,
    n::Integer,
    kk::T,
) where {T <: AbstractFloat}
    cell = @index(Global)
    @unpack tv_opt, pvd_max = CFT

    # One thread owns one cell. Maintain the n smallest monthly values in
    # sorted order without allocating or sorting a 12-by-cells temporary.
    for rank in 1:n
        climbuf_min_temp[rank, cell] = typemax(T)
    end
    for month in axes(climbuf_mtemp20, 1)
        candidate = climbuf_mtemp20[month, cell]
        for rank in 1:n
            if candidate < climbuf_min_temp[rank, cell]
                candidate, climbuf_min_temp[rank, cell] =
                    climbuf_min_temp[rank, cell], candidate
            end
        end
    end

    sum_v_req = zero(T)
    for rank in 1:n
        temperature = climbuf_min_temp[rank, cell]
        if temperature <= tv_opt.low && temperature > T(-9999)
            sum_v_req += pvd_max / T(n)
        elseif temperature > tv_opt.low && temperature < tv_opt.high
            sum_v_req += pvd_max / T(n) *
                (one(T) - (temperature - tv_opt.low) / (tv_opt.high - tv_opt.low))
        end
    end
    climbuf_V_req_a[cell] = sum_v_req
    if climbuf_V_req[cell] < -9998
        climbuf_V_req[cell] = sum_v_req
    else
        climbuf_V_req[cell] = (one(T) - kk) * climbuf_V_req[cell] + kk * sum_v_req
    end

    annual_temperature = zero(T)
    for day in axes(daily_temp, 1)
        annual_temperature += daily_temp[day, cell]
    end
    climbuf_atemp_mean[cell] = annual_temperature / T(size(daily_temp, 1))
end

function monthlytemp!(daily_temp::AbstractArray{T},
                      climbuf_mtemp::AbstractArray{T}
) where {T <: AbstractFloat}
    """
    Calculate the average temperature for each month.

    Args:
        daily_temps::Vector{Float64}: Daily temperature data (length is 365)

    Return:
        A vector of length 12 representing the average temperature for each month.
    """
    # Month metadata is copied to the active device to avoid host reads inside kernels.
    # ndaymonth = device([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    # start_indices = device([1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335])
    
    # start_indices = cumsum(vcat(1, ndaymonth[1:end-1]))
    # mtemp = similar(daily_temp, (12, cell_size))  # Store mean temperatures for each month
#     start_idx = 1  # Start index for each month in daily_temps

#     for month = 1:12
#         end_idx = start_idx + ndaymonth[month] - 1  # End index for the current month
#         mtemp[month] = mean(daily_temp[start_idx : end_idx])  # Calculate the monthly average
#         start_idx = end_idx + 1  # Update start index for the next month
#     end
    
    launch_2D!(
        monthlytemp_kernel!,
        climbuf_mtemp,
        daily_temp
    )
    
end


@kernel inbounds = true function monthlytemp_kernel!(
                                     climbuf_mtemp::AbstractArray{T}, 
                                     daily_temp::AbstractArray{T}
) where {T <: AbstractFloat}
    
    # launch layout is (month, cell).
    month, cell = @index(Global, NTuple)

    # compile-time constants
    ndaymonth = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    start_indices  = (1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)

    start_idx = start_indices[month]
    days = ndaymonth[month]

    sum_temp = zero(T)
    
    for i in 0:(days - 1)
        sum_temp += daily_temp[start_idx + i, cell]
    end

    climbuf_mtemp[month, cell] = sum_temp / days 

end


"""
daily_climbuf!(temp, climbuf_temp)

Advance the rolling daily temperature buffer by one day.
"""
function daily_climbuf!(temp::AbstractArray{T},
                        climbuf_temp::AbstractArray{T};
                        annual_temperature::AbstractArray{T} = climbuf_temp,
                        annual_day::Integer = 0,
) where {T <: AbstractFloat}

    kernel_params = (NDAYS = 31,)

    launch_1D!(
        daily_climbuf_kernel!,
        temp,
        climbuf_temp,
        annual_temperature,
        annual_day,
        kernel_params
    )

end


@kernel inbounds = true function daily_climbuf_kernel!(
                                       temp::AbstractArray{T},
                                       climbuf_temp::AbstractArray{T},
                                       annual_temperature::AbstractArray{T},
                                       annual_day::Integer,
                                       kernel_params
) where {T <: AbstractFloat}

    cell = @index(Global)

    @unpack NDAYS = kernel_params

    # Shift the rolling daily climate buffer left and append today's temperature.
    for day in 2:NDAYS
        climbuf_temp[day-1, cell] = climbuf_temp[day, cell]
    end
    climbuf_temp[NDAYS, cell] = temp[cell]
    if annual_day != 0
        annual_temperature[annual_day, cell] = temp[cell]
    end

end
