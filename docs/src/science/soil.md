# Soil processes

This page documents the five-layer soil water, thermal, carbon, and nitrogen
pathway used by the current production configuration.

## Soil water state

Each layer separates water below wilting point, plant-available water, free
water above field capacity, and the corresponding ice pools. Relative
plant-available water is diagnosed as

```math
w_l=\frac{W_l-W_{wilt,l}}{W_{fc,l}-W_{wilt,l}},
```

with storage-aware clipping in the implementation. Rain and snowmelt enter the
top boundary after interception. Infiltration is processed in bounded slugs for
numerical stability. Every layer enforces

```math
0\le W_l\le W_{sat,l},
```

and partitions excess water into vertical percolation and runoff. Water leaving
the bottom layer is drainage. Nitrate moves with percolating water using the
source-layer concentration and cannot exceed the source nitrate pool.

Full irrigation restores each rooted layer to field capacity every day:

```math
I_l=\max(0,W_{fc,l}-W_l),\qquad W_l\leftarrow W_l+I_l.
```

The added water is an unconstrained boundary input; this option is not a water
allocation, river, reservoir, or irrigation-infrastructure model. Irrigated
patches therefore do not report a closed daily water-balance residual; rainfed
patches continue to use the water-balance diagnostic.

## Soil thermal state

Temperature is recovered from layer enthalpy rather than advanced as an
independent diagnostic. Conceptually,

```math
H_l=C_lT_l-L_fW_{ice,l},
```

where effective heat capacity ``C_l`` depends on mineral soil, organic matter,
liquid water, and ice. Interface conductivity uses the harmonic mean

```math
k_{l+1/2}=\frac{2k_lk_{l+1}}{k_l+k_{l+1}}.
```

Liquid water carries enthalpy between layers. Rain enters at air temperature,
snowmelt at 0 °C, evaporation and transpiration remove source-layer water
enthalpy, and bottom drainage exports source-layer enthalpy. After every water
transfer, the implementation reconstructs temperature, liquid/ice partition,
heat capacity, latent heat, and conductivity.

Consequently saturation fraction and previous-step thermal properties belong
to the numerical state lifecycle. They are not pure end-of-day diagnostics.

## Decomposition response

The liquid pore-space moisture factor is

```math
\omega_l=\operatorname{clamp}\left(
\frac{w_lW_{hc,l}+W_{wilt,l}-W_{wilt,ice,l}+W_{free,l}}
{W_{sat,l}-W_{wilt,ice,l}-W_{avail,ice,l}-W_{free,ice,l}},
\epsilon,1\right).
```

The temperature response is

```math
f_T(T)=\begin{cases}
0,&T<-15,\\
\exp\left[e_0\left(\dfrac{1}{T_r+10}-
\dfrac{1}{\min(T,40)+T_r}\right)\right],&T\ge-15.
\end{cases}
```

Temperature and moisture combine as

```math
R_l=\operatorname{clamp}\left[
f_T(T_l)(a_0+a_1\omega_l+a_2\omega_l^2+a_3\omega_l^3),0,1\right].
```

Surface litter uses surface-litter temperature and wetness **without** the
upper bound of one, following LPJmL 5.10/6.1's methane-disabled pathway. Its
response is zero when the topsoil temperature response is zero. Incorporated
and root litter retain the bounded topsoil response above.

## Soil carbon

For stock ``C`` with rate ``k`` and environmental response ``R``, exact daily
first-order loss is

```math
\Delta C=-\operatorname{expm1}(-kR)C
=\left(1-e^{-kR}\right)C.
```

The `expm1` form avoids cancellation for small Float32 rates. Litter loss is
divided between atmospheric respiration and retained fast/slow material:

```math
C_{fast,in,l}=s_{f,l}f_{fast}(1-f_{atm})\Delta C_{lit},
```

```math
C_{slow,in,l}=s_{s,l}(1-f_{fast})(1-f_{atm})\Delta C_{lit},
```

where ``s_f`` and ``s_s`` are normalized vertical routing fractions.
Heterotrophic respiration is

```math
R_h=f_{atm}\Delta C_{lit}+
\sum_l(\Delta C_{fast,l}+\Delta C_{slow,l}).
```

Tillage changes topsoil hydraulic properties and routes incorporated litter;
bioturbation provides the no-tillage vertical pathway. Harvest inputs are
routed after decomposition, so new residue first decomposes the following day.
Surface-litter cover and water capacity are refreshed immediately after
decomposition, before the day's interception and evaporation. Water exceeding
the reduced litter capacity returns to the topsoil, conserving total water.

## Organic nitrogen and immobilization

Organic litter, fast, and slow nitrogen use the same response, decay fractions,
and vertical routing as carbon. Gross mineralization is

```math
M_l=M_{lit,l}+\Delta N_{fast,l}+\Delta N_{slow,l},
```

and enters ammonium. If retained litter carbon requires more nitrogen than
decomposed litter supplies,

```math
N_{def}=\frac{\Delta C_{lit}}{(C:N)_{soil}}-\Delta N_{lit}>0,
```

the model immobilizes mineral nitrogen into fast and slow organic pools. The
response to mineral-N concentration ``c_N`` is

```math
f_{immob}=\frac{c_N}{K_{immob}+c_N}.
```

Removal is divided between ammonium and nitrate according to availability and
cannot make either pool negative.

## Mineral nitrogen transformations

Nitrification converts ammonium to nitrate using pH, temperature, and
water-filled pore-space responses. A configured fraction becomes N₂O and the
remainder enters nitrate. Denitrification acts after crop uptake and is limited
by nitrate, soil carbon, temperature, pH, and anaerobic moisture response. Its
products are partitioned between N₂O and N₂.

Surface NH₃ volatilization acts on the **entire top-layer ammonium pool each
day**, not only on fertilizer applied that day. Its response uses top-layer pH,
temperature, and wind speed. Cumulative loss can therefore remain large when
mineralization and fertilization repeatedly replenish ammonium.

The daily mineral-N sequence is:

1. organic N mineralization and immobilization;
2. nitrification;
3. fertilizer/manure application and crop uptake at their positions in the
   daily process order;
4. nitrate transport with water;
5. denitrification and NH₃ volatilization.

## Tillage and litter routing

Tillage is a configuration switch rather than a gridded input. When enabled,
it modifies the topsoil hydraulic state and incorporates eligible surface
litter. Without tillage, bioturbation controls downward redistribution. Both
pathways preserve litter C/N accounting and feed the same decomposition pools.

## Code map

- `src/processes/soil/soil_water.jl`
- `src/processes/soil/infil_perc.jl`
- `src/processes/soil/evaporation.jl`
- `src/processes/soil/soil_temp.jl`
- `src/processes/soil/soil_response.jl`
- `src/processes/soil/soil_carbon.jl`
- `src/processes/soil/nitrogen_transform.jl`
- `src/processes/soil/litter_routing.jl`
- `src/processes/soil/tillage.jl`
