# Checkpoints

Save only at a completed daily boundary:

```julia
save_checkpoint("checkpoint.jld2", simulation)
```

Restore into a newly initialized, compatible simulation:

```julia
resumed = initialize_simulation(
    cft1, initial_data;
    indices = [1], T = Float32, days = total_days,
    auto_fertilizer = false,
)
restore_checkpoint!(resumed, "checkpoint.jld2")
run_simulation!(resumed, remaining_climate; spinup = false)
```

The v2 checkpoint records prognostic state, restart-relevant inputs, partial
outputs, diagnostics, process parameters, and run metadata. Arrays are stored
on the host and copied to the target backend during restoration.

Precision, cell count, configured duration, photosynthetic pathway, and major
run options must match. Current v1 development checkpoints are not supported.
Checkpoint compatibility is not yet guaranteed across arbitrary package
versions; archive the exact Agrocosm commit with long-lived experiments.
