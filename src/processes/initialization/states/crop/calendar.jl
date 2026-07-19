"""Crop sowing and harvest calendar state."""
mutable struct CropCalendar{I}
    sowing_date::I
    harvest_date::I
    sowing_callback::I
    harvest_callback::I
    harvesting_year::I
end

function init_crop_calendar(cell_size::Int, device)
    int_state() = device(zeros(Int32, cell_size))
    return CropCalendar(ntuple(_ -> int_state(), 5)...)
end
