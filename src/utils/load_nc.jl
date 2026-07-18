"""
load_nc_file_one_dimension(file_path, variable)

Load a 1D variable from a NetCDF file.
"""
function load_nc_file_one_dimension(file_path::String, 
                                    variable::String,
                                    timerange::UnitRange
)

    ds = NCDataset(file_path, "r")

    dataset = ds[variable][:, :, timerange]

    close(ds)

    return dataset
end

"""
load_nc_file_dimensions(file_path, variable)

Load a multi-dimensional variable from a NetCDF file.
"""
function load_nc_file_dimensions(file_path::String, 
                                 variable::String,
                                 timerange::UnitRange
)

    ds = NCDataset(file_path, "r")

    dataset = ds[variable][:, :, :, timerange]

    close(ds)

    return dataset
end
