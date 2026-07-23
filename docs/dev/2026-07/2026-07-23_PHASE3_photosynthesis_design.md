# Phase 3 — Crop C3/C4 photosynthesis: port design

> Companion to the migration plan. Records the exact legacy LPJmL photosynthesis physics and the
> continuous-time reformulation used for the Terrarium `CropPhotosynthesis <: AbstractPhotosynthesis`
> port. Legacy reference: `src/processes/crop/{photosynthesis,lambda_solver,radiation}.jl`,
> `src/processes/climate/temp_stress.jl` (retained on disk).

## Key finding: LPJmL C3 ≈ Terrarium `LUEPhotosynthesis`

Terrarium's `LUEPhotosynthesis` (`Terrarium/src/processes/vegetation/photosynthesis/lue_photosynthesis.jl`)
already implements the BIOME3 / Haxeltine & Prentice (1996) mechanistic C3 photosynthesis that LPJmL
derives from — the **instantaneous (per-second)** form of the same co-limitation biochemistry:

- `compute_kinetic_parameters` (τ, Kc, Ko via Q10) ↔ LPJmL `ko/kc/tau` (photosynthesis.jl:236-239).
- `compute_Γ_star = pO₂/(2τ)` ↔ LPJmL `gammastar` (:240).
- `compute_temperature_stress` (double sigmoid with `k1,k2,k3`) is **identical** to LPJmL `temp_stress`
  (temp_stress.jl:38-55) — same `k1 = 2 log(1/0.99−1)/(T_lo−T_photlo)`, `k2 = (T_lo+T_photlo)/2`,
  `k3 = log(0.99/0.01)/(T_hi−T_phothi)`. LPJmL adds a daylength gate and a C3/C4 hard cutoff (`tmc3=45`,
  `tmc4=55`).
- `compute_assimilation_factors` (c₁, c₂), `compute_JE_JC`, `compute_Ag` (θ co-limitation quadratic) ↔
  LPJmL `c1,c2,je,jc,agd` (:261-270).
- `Rd = α·Vcmax·β` ↔ LPJmL `leaf = b·vcmax` (:273).

So the continuous crop C3 photosynthesis **is** `LUEPhotosynthesis` with crop parameters. The genuinely
crop-specific additions are (1) the **C4 pathway** and (2) LPJmL's **λ water-coupling solver**.

## Continuous-time reformulation (the daily→per-second decision)

LPJmL is daily-integrated; the daily artifacts must be dropped, not reproduced:

- **Instantaneous PAR.** Use `PAR = ½·swdown·(1−α_leaf)·cq` (mol/m²/s) as in `compute_PAR`, not the
  daily `par = 86400·swdown/2` (radiation.jl:103). `APAR = α_a·PAR·(1−exp(−k_ext·LAI))`.
- **Drop the `/daylength … ·daylength` bracketing.** In LPJmL `je = c₁·apar/daylength` then
  `agd = (…)·daylength`; the daylength factors cancel in the instantaneous limit. The continuous rate
  is `Ag = (JE+JC−√((JE+JC)²−4θ·JE·JC))/(2θ)·β` on instantaneous JE, JC (exactly `compute_Ag`).
- **`vcmax`/`Rd` are instantaneous** (gC/m²/s), no `/24` (`hour2day`) or `daylength/24` night scaling.
- Outputs match the interface: `net_assimilation` (An, gC/m²/s), `leaf_respiration` (Rd), and
  `gross_primary_production` (GPP = An·1e-3, kgC/m²/s).

Consequence: the port does **not** reproduce the legacy *daily totals* bit-for-bit (different quantity —
instantaneous rate vs day-integral). Validation is at the **primitive** level (τ, Kc, Ko, Γ\*, c₁, c₂,
temperature stress, the θ co-limitation) against the legacy scalar equations, plus physical invariants
(An ≥ 0; An = 0 below −3 °C, without light, or when `T_stress < 1e-2`; monotonic in APAR; C4 CO₂
saturation via `φ(pᵢ) = min(1, λ/λ_mc4)`).

## C4 pathway (new vs `LUEPhotosynthesis`)

C4 drops the Γ\*/Kc/Ko chain: `c₂ ≡ 1`, and the light term is scaled by `φ = min(1, λ/λ_mc4)`
(`lambdamc4 = 0.4`): `c₁ = T_stress·φ·α_C4`, `JE = c₁·APAR`, `JC = Vcmax` (instantaneous),
`Ag = θ-co-limitation(JE, JC)`. `α_C4 = 0.053`. High-T cutoff `tmc4 = 55 °C`.

## λ water-coupling → a crop stomatal-conductance process

Terrarium separates stomatal conductance (produces `leaf_to_air_co2_ratio` λc) from photosynthesis;
LPJmL couples them by solving λ so diffusive **supply** equals biochemical **demand**:

```
find λ ∈ [0.02, 0.85] :  fac·(1−λ) − A_mm(λ) = 0
fac   = gpd/1.6·co2·1e-5,   gpd = daylength·3600·(g_canopy − gmin·fpar)   (Pa→bar via 1e-5)
A_mm(λ) = net assimilation at trial λ, converted gC→mm by ideal gas ·8.314·(T+273.15)/p·1000
```
solved by the LPJmL 30-step "return best sample" bisection (`lpj_bisect`, already retained). This maps
to a `CropStomatalConductance <: AbstractStomatalConductance` producing `leaf_to_air_co2_ratio` (λ) and
`canopy_water_conductance`. **Staging:** first port `CropPhotosynthesis` taking λc as input (assembled
with Terrarium's `MedlynStomatalConductance` for λc); then replace Medlyn with the LPJmL λ solver.

## Interface contract (from `LUEPhotosynthesis`)

```
variables(::CropPhotosynthesis) = (
    auxiliary(:net_assimilation, XY(), units=u"g/m^2/s"),
    auxiliary(:leaf_respiration, XY(), units=u"g/m^2/s"),
    auxiliary(:gross_primary_production, XY(), units=u"kg/m^2/s"),
    input(:soil_moisture_limiting_factor, XY(), default=1),
    input(:leaf_area_index, XY()),
)
compute_photosynthesis(i, j, grid, fields, photo, constants, atmos) -> (Rd, An, GPP)
compute_auxiliary!(state, grid, photo, stomcond, constants, atmos, args...)   # launches XY kernel
```
Reads: `air_temperature`, `air_pressure`, `shortwave_down` (atmos accessors), `fields.CO2`,
`fields.soil_moisture_limiting_factor` (β), `fields.leaf_area_index`, `fields.leaf_to_air_co2_ratio` (λc).
C3/C4 selected by a type parameter or `PFT.path`; crop temperature thresholds and `α_C3/α_C4`, `tmc3/tmc4`
come from the CFT registry (Phase 5 wires per-CFT presets).

## Legacy parameter defaults (for the crop presets)

`PhotoParams`: po2=20.9e3 Pa, q10ko=1.2, q10kc=2.1, q10tau=0.57, tau25=2600, cmass=12, cq=4.6e-6,
α_C3(alphac3)=0.08, α_C4(alphac4)=0.053, λ_mc4(lambdamc4)=0.4, tmc3=45, tmc4=55, θ(theta)=0.9.
`PftParameters`: b=0.031 (leaf resp fraction), gmin (0.99–1.6 per CFT), lightextcoeff=0.5,
albedo_leaf=0.18, temp_co2/temp_photos per CFT.

## Status and integration finding (2026-07-23)

`CropPhotosynthesis` is implemented (`src/crop/photosynthesis.jl`) and unit-validated:
`test/crop/test_photosynthesis.jl` proves the C3 scalar `(Rd, An)` matches Terrarium's
`LUEPhotosynthesis` across 72 input combinations at `rtol = 1e-10`, and checks C4 gating, light
response, φ saturation, and the 55 °C/45 °C cutoffs.

**Integration blocker (upstream, small):** assembling `VegetationCarbon(NF; photosynthesis =
CropPhotosynthesis(NF))` fails because `MedlynStomatalConductance.compute_auxiliary!`
(`Terrarium/src/processes/vegetation/stomatal_conductance/medlyn_stomatal_conductance.jl:87-90`) and
`compute_stomatal_conductance` (:109, :131) dispatch on `photo::LUEPhotosynthesis` **specifically**,
not `AbstractPhotosynthesis` — even though they only read the `net_assimilation` field that every
`AbstractPhotosynthesis` produces. Two resolutions, both clean:

1. **Upstream (preferred, one-line):** widen those three Medlyn signatures from `LUEPhotosynthesis`
   to `AbstractPhotosynthesis`. Then any crop photosynthesis injects into `VegetationCarbon`
   unchanged. Candidate for a small Terrarium PR.
2. **Downstream:** implement the paired crop stomatal conductance (the LPJmL λ water-coupling solver,
   §"λ … → a crop stomatal-conductance process") as `CropStomatalConductance <:
   AbstractStomatalConductance`, and inject both together. This is the planned next Phase-3 step and
   removes the dependence on Medlyn entirely.

Until one of these lands, the crop photosynthesis is validated at the process/unit level (which is
the physics that matters); the end-to-end `LandModel` assembly with the crop photosynthesis is
deferred to the crop stomatal-conductance port.
