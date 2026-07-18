"""
Standardize KernelAbstractions launches.
These wrappers keep backend/ndrange/synchronize logic in one place.
"""

function launch_1D!(kernelfun, ref_array, args...)
    # Convention: ref_array is both the launch reference and kernel arg #1.
    # This keeps call sites concise while preserving explicit launch shape.
    backend = KernelAbstractions.get_backend(ref_array)
    kernel = kernelfun(backend)
    kernel(ref_array, args..., ndrange = length(ref_array))
    KernelAbstractions.synchronize(backend)
    return nothing
end

function launch_2D!(kernelfun, ref_array, args...)
    # 2D launch assumes ref_array uses (dim1, dim2) layout consistent with @index(..., NTuple).
    backend = KernelAbstractions.get_backend(ref_array)
    kernel = kernelfun(backend)
    kernel(ref_array, args..., ndrange = (size(ref_array, 1), size(ref_array, 2)))
    KernelAbstractions.synchronize(backend)
    return nothing
end

function launch_custom!(kernelfun, ref_array, ndrange, args...)
    # Use custom launch only when ndrange does not match length/size(ref_array).
    backend = KernelAbstractions.get_backend(ref_array)
    kernel = kernelfun(backend)
    kernel(ref_array, args..., ndrange = ndrange)
    KernelAbstractions.synchronize(backend)
    return nothing
end
