"""
spin_up_climbuf!(PFT, climate, climbuf, day, lat, temp, lwnet, swdown, output)

Run climate-buffer spin-up for one step before full crop process integration.
"""
function spin_up_climbuf!(PFT::PftParameters, 
                          temp_spinup::AbstractArray{T}, 
                          climbuf::ClimBuf;
                          year_spinup = 1
) where {T <: AbstractFloat}
    for i = 1 : year_spinup
        year_temp = temp_spinup[365*(i-1)+1 : 365*i, :]
        for day in axes(year_temp, 1)
            daily_climbuf!(year_temp[day, :], climbuf.temp)
        end
        annual_climbuf!(year_temp, climbuf, PFT)
    end

end


"""
update_climbuf!(PFT, climbuf, day, lat, temp, lwnet, swdown)

Update climate-buffer and PET diagnostics during daily simulation.
"""
function update_climbuf!(PFT::PftParameters, 
                         temp::AbstractArray{T},
                         climbuf::ClimBuf,
                         day::Integer
)where {T <: AbstractFloat}

    if day > 1 && day % 365 == 1
        annual_climbuf!(climbuf.atemp, climbuf, PFT)
    end
    
    if day % 365 == 0
        day = 365
    else
        day = day % 365
    end
    
    daily_climbuf!(
        temp, climbuf.temp;
        annual_temperature = climbuf.atemp,
        annual_day = day,
    )

end
