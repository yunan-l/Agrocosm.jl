# Initialization and warm-up

This page summarizes how source soil observations become Agrocosm state and
how the finite agricultural warm-up changes that state before production.
Implementation and operational details are in
[HWSD soil initialization](../guide/hwsd_initialization.md).

## Canonical spatial alignment

Climate, management, soil type, and HWSD-derived C/N are aligned through the
canonical `grid.nc` mapping

```math
i_{cell}\longleftrightarrow(i_{lon},i_{lat})
\longleftrightarrow(\lambda_{geo},\varphi_{geo}).
```

Source dimension order is not assumed. Every compact array retains canonical
`cellid`, and outputs use that key when reconstructed to the 720 × 280 grid.

## HWSD concentrations to stocks

For organic-carbon concentration ``OC`` in wt%, bulk density ``\rho_b`` in
g cm⁻³, layer thickness ``\Delta z`` in cm, and coarse-fragment percentage
``f_{cf}``, layer SOC stock is

```math
C_{SOC}=OC\,\rho_b\,\Delta z
\left(1-\frac{f_{cf}}{100}\right)100
\quad[\mathrm{gC\ m^{-2}}].
```

For total nitrogen ``N_{conc}`` in g kg⁻¹,

```math
N_{tot}=N_{conc}\,\rho_b\,\Delta z
\left(1-\frac{f_{cf}}{100}\right)10
\quad[\mathrm{gN\ m^{-2}}].
```

Multiple HWSD soil components are combined with their original component
shares. Missing layers below a component's declared root depth are structural
zeros; valid deep components are not renormalized upward.

## Horizontal and vertical aggregation

The 30-arc-second HWSD pixels inside an Agrocosm cell are aggregated with
spherical pixel-area weights:

```math
\bar X_g=\frac{\sum_{p\in g}A_pX_p}{\sum_{p\in g}A_p}.
```

Water, glacier, and NODATA pixels are excluded from the soil-area denominator.
HWSD's seven layers are remapped to Agrocosm's five layers by depth overlap:

```math
X_j=\sum_i X_i
\frac{|[z_{i,0},z_{i,1}]\cap[z_{j,0},z_{j,1}]|}
{z_{i,1}-z_{i,0}}.
```

HWSD ends at 2 m. The default 2--3 m rule extends the 1.5--2 m stock density
and marks the result uncertain. True unresolved gaps may use a bounded nearest
complete profile, with donor coordinates and distance retained as provenance.

## Constructing model pools

HWSD supplies total SOC and total N, not Agrocosm litter, fast, and slow pools.
The current interim initialization is

```math
C_{fast,l}=0.4C_{SOC,l},\qquad
C_{slow,l}=0.6C_{SOC,l},\qquad
C_{litter,l}=0,
```

with the same 40:60 split for organic nitrogen after reserving the configured
initial mineral pools:

```math
N_{org,l}=\max(0,N_{tot,l}-N_{NO_3,l}-N_{NH_4,l}),
```

```math
N_{fast,l}=0.4N_{org,l},\qquad
N_{slow,l}=0.6N_{org,l},\qquad
N_{litter,l}=0.
```

Liquid soil water begins at field capacity:

```math
W_l=W_{fc,l}.
```

The 40:60 partition approximates the mean legacy ten-cell partition but is an
explicit initialization assumption, not an HWSD measurement.

## Adaptive target-constrained agricultural warm-up

The 40:60 partition is a starting condition. Before reported production,
Agrocosm runs at least ten agricultural years using the same crop, fertilizer,
manure, irrigation, residue, and tillage configuration. If ``n_f`` complete
forcing years are supplied, warm-up year ``y`` selects

```math
y_f=1+\operatorname{mod}(y-1,n_f).
```

All crop, water, heat, C, and N processes run normally. At each year end, the
mineral-soil pools are returned to the initial HWSD layer targets while their
internal fractions remain process-determined:

```math
C_{fast,l}+C_{slow,l}=C_{SOC,l},
```

```math
N_{fast,l}+N_{slow,l}+N_{NO_3,l}+N_{NH_4,l}=N_{tot,l}.
```

Litter is not constrained. Only lifecycle state is retained:

```math
x_{warm}^{(y+1)}=
\mathcal{T}_{365}\left(x_{warm}^{(y)},u^{(y_f)},m^{(y_f)};\vartheta\right).
```

Production time and output remain unchanged,

```math
d_{production}=0,\qquad Y_{production}=\varnothing,
```

until the formal run begins. After the ten-year minimum, annual total-C/N
changes, fast-pool fraction changes, and target corrections must remain below
their thresholds for the configured consecutive years. Otherwise warm-up
continues to its maximum duration. Annual reports retain convergence status,
target corrections, litter, fast, slow, total C/N, mineral N, and soil water.
The final state is saved as a native Agrocosm checkpoint.

## Interpretation

Ten years is the minimum, not an equilibrium claim. A run that reaches its
maximum duration without satisfying the per-cell criteria is reported as
`target_constrained_maximum_years` and remains a baseline rather than an
accepted equilibrium initialization. The annual correction itself is an
initialization diagnostic and is never introduced as a production flux.

Warm-up, HWSD preprocessing, input loading, checkpoint I/O, and reporting stay
outside the future Enzyme-differentiable one-day transition.

## Code map

- `lib/AgrocosmData/src/hwsd.jl`
- `lib/AgrocosmData/src/soil.jl`
- `src/simulations/agricultural_warmup.jl`
- `docs/src/guide/hwsd_initialization.md`
