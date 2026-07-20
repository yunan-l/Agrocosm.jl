"""Persistent crop-calendar state."""
mutable struct CropCalendarState{I}
    sowing_date::I
    harvest_date::I
    harvesting_year::I
end

"""Discrete crop events that occur during the current day."""
mutable struct CropEvents{I}
    sowing::I
    harvest::I
end

function init_crop_calendar_state(cell_size::Int, device)
    int_state() = device(zeros(Int32, cell_size))
    return CropCalendarState(ntuple(_ -> int_state(), 3)...)
end

function init_crop_events(cell_size::Int, device)
    int_event() = device(zeros(Int32, cell_size))
    return CropEvents(ntuple(_ -> int_event(), 2)...)
end
