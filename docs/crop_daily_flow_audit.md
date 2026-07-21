# Prescribed crop lifecycle and daily-flow audit

This audit covers the current single-crop, prescribed-sowing pathway against
LPJmL commit `572e2b906ac2c55b2ee6661a93e4633b126254e4`.

## Completed stages

- **S0 — event timeline:** two consecutive years produce one sowing and one
  harvest transition per year. Growth is active only between those events and
  `fphu` is monotonic within each crop season.
- **S1 — cultivation and management:** cultivation, tillage, seed input and
  the first prescribed fertilizer fraction occur only on the sowing event.
  The pending fraction is applied once after the previous day's `fphu` exceeds
  0.25, matching `daily_agriculture.c` before its current-day phenology call.
- **S2 — prescribed phenology:** heat units, senescence and the `hlimit`
  fallback follow `phenology_crop.c`. Reaching PHU sets `fphu = 1`; harvest is
  triggered at the start of the following daily phenology step, as in LPJmL.
- **S3 — harvest routing:** storage is harvested product, the removed fraction
  of leaf and pool biomass is exported, retained shoots enter surface litter,
  and roots enter root litter. Carbon and nitrogen close independently.
- **S4 — cross-year continuity:** sowing reconstructs seasonal crop state but
  does not reset soil water, heat, carbon or nitrogen. One 730-day climate block
  and two continuous 365-day blocks produce the same crop events and soil state.
- **S5 — API regression:** the high-level API passes the 20-year rainfed-wheat
  notebook and the complete CPU test suite. A separate CUDA lifecycle test
  covers the two-year event sequence and CPU/GPU agreement.

## Persistent struct versus LPJmL allocation

LPJmL deletes the harvested crop PFT and allocates a new `Pftcrop` at the next
cultivation. Agrocosm keeps one structure-of-arrays allocation for efficient GPU
execution. The sowing kernel therefore explicitly reconstructs the seasonal
fields initialized by `new_crop.c`: phenology sums and flags, canopy seed state,
seed C/N, organ N, pending management inputs, and seasonal water/N stress sums.
Annual output buffers and all soil states deliberately remain persistent.

After harvest, `harvesting = true` is Agrocosm's inactive/killed sentinel.
From the following day until sowing, active crop C/N stocks, phenology sums,
LAI, photosynthetic fluxes, and seasonal water/N accumulators are zero. Static
configuration (`phu`, root distribution and seed-organ templates) remains
allocated. Bare-soil albedo, root-zone soil water, temperature response and the
environmental response factors are diagnostics/workspace rather
than living-crop stocks and therefore need not be zero.

## Input safety

A one-dimensional CO₂ forcing is an annual global series; a matrix is daily
cell-resolved forcing. The high-level API now checks the required number of
annual values or daily rows before launching an `@inbounds` CPU/GPU kernel, so a
multi-year run cannot silently read beyond the forcing array.

## Deferred scope

Dynamic sowing dates, multiple simultaneous stands/crops, double cropping and
stand-fraction competition remain second-stage work. They are not approximated
inside the current prescribed single-crop pathway.
