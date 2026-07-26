# Crop processes

This page collects the crop-side equations and management logic used by the
current default pathway.

## Radiation and albedo

Day length, absorbed PAR, and equilibrium evaporation are computed from solar
geometry and daily radiation forcing. For most crops, the absorbed green
fraction is

```math
f_g=1-\exp(-k_LL_{act}),
```

with actual LAI

```math
L_{act}=\max(0,L-L_{deficit}).
```

Maize uses the LPJmL empirical APAR parameterization instead. Surface albedo is
assembled from green canopy, litter cover, bare soil, and snow when present.

## Phenology and LAI

Heat units accumulate above the base temperature, winter crops may require
vernalization, and the photoperiod factor is applied until the senescence
threshold is reached. Fractional development is

```math
f_{PHU}=\min\left(1,\frac{HU_\Sigma}{PHU}\right).
```

Potential LAI follows the LPJmL-style growth curve before senescence and is
then scaled from the LAI reached at senescence onset. Water and nitrogen stress
reduce the daily LAI increment.

## C3 and C4 photosynthesis

C3 photosynthesis uses the Haxeltine–Prentice style co-limited formulation.
Temperature-adjusted Michaelis constants, Rubisco specificity, and internal CO₂
ratio `p_i=\lambda p_a` are used to compute light-limited and Rubisco-limited
rates. Gross assimilation is the non-rectangular hyperbola combination of those
two limits.

C4 photosynthesis uses the same co-limitation equation, but the internal-CO₂
response is capped by the C4 threshold and the Rubisco branch is saturated in
the current formulation.

## Respiration and conductance

Daily leaf dark respiration is

```math
R_d=bV_{c\max}.
```

Daytime net assimilation subtracts the daylight-scaled leaf respiration from
gross assimilation. That daytime net assimilation is then converted into the
canopy conductance target used by the water-limitation solver.

The internal-CO₂ ratio `\lambda` is solved by bounded bisection from the water
supply and assimilation closure equation.

## Transpiration and water stress

Root-weighted relative soil water provides potential supply:

```math
w_r=\sum_l r_lw_l,\qquad
S=e_{max}w_r(1-e^{-0.04C_{root}}).
```

Atmospheric demand depends on equilibrium evaporation, canopy wetness, and
canopy conductance. Actual transpiration is the minimum of supply and demand,
with each soil layer capped by its available storage. The seasonal water-deficit
diagnostic is reported on the LPJmL 0–100 scale.

## Nitrogen demand and uptake

Leaf nitrogen demand combines Rubisco demand and minimum structural leaf N.
Total crop demand then adds root, pool, and storage demand through the PFT C:N
ratios.

Mineral uptake is split across nitrate and ammonium with layer-wise
Michaelis–Menten kinetics. Uptake is limited by root distribution, soil
temperature, soil water, and available mineral pools. When automatic fertilizer
is enabled, any remaining deficit is supplied as a boundary input.

## Fertilizer and manure

`fertilizer` is configured separately from `manure`.

- `fertilizer = :no` means no prescribed or automatic mineral fertilizer.
- `fertilizer = :yes` means use prescribed mineral-fertilizer input.
- `fertilizer = :auto` means satisfy crop mineral-N demand automatically.

Manure is controlled by a separate Boolean option and is applied independently
of the mineral-fertilizer mode.

Prescribed fertilizer and manure are split between sowing and a later split
application. Mineral fertilizer is partitioned between nitrate and ammonium.
Manure contributes mineral ammonium plus incorporated organic litter C/N.

## Respiration, NPP, and allocation

Daily NPP is

```math
NPP=A_g-R_d-(R_{root}+R_{storage}+R_{pool}+R_g),
```

with growth respiration a fixed fraction of the remaining carbon after the
maintenance terms. Negative NPP is retained. If biomass plus NPP becomes
non-positive, Agrocosm terminates the stand using the LPJmL-style failed-crop
pathway.

Living biomass is partitioned into root, leaf, storage, and mobile-pool carbon.
Root fraction declines with development and stress, leaf carbon is limited by
SLA and LAI, storage carbon follows the harvest-index curve, and the mobile
pool closes the biomass balance exactly.

## Harvest and failure

Normal harvest routes storage carbon to yield and sends the configured residue
fraction to litter. Failed-crop termination removes the stand immediately,
exports the appropriate above-ground material, routes roots and residues to
litter, and resets crop state.

## Code map

- `src/processes/crop/radiation.jl`
- `src/processes/crop/albedo.jl`
- `src/processes/climate/snow.jl`
- `src/processes/crop/phenology.jl`
- `src/processes/crop/photosynthesis.jl`
- `src/processes/crop/lambda_solver.jl`
- `src/processes/crop/transpiration.jl`
- `src/processes/crop/respiration.jl`
- `src/processes/crop/carbon_allocation.jl`
- `src/processes/crop/nitrogen_demand.jl`
- `src/processes/crop/nitrogen_uptake.jl`
- `src/processes/crop/fertilizer.jl`
- `src/processes/crop/harvesting.jl`
