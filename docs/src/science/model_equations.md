# Model equations

This page documents the equations and numerical choices in the **current
default Agrocosm process pathway**. The formulation is primarily informed by
LPJmL crop and soil processes, but the implementation uses Agrocosm's own
lifecycle state and cell-local kernels. This is an implementation description:
when a general textbook formulation differs from the code, the equations below
describe the code.

## Scope and conventions

The present production configuration assumes:

- one crop PFT per simulated stand and grid cell;
- a five-layer soil column;
- daily forcing on a strict 365-day, no-leap calendar;
- rainfed conditions or unconstrained full irrigation;
- C3 or C4 photosynthesis;
- coupled crop carbon, crop nitrogen, soil water, soil heat, and soil C/N;
- no natural-vegetation competition, groundwater abstraction, or lateral
  exchange between grid cells;
- no dedicated frozen-soil infiltration impedance or advective heat transport
  through frozen soil.

Unless stated otherwise, carbon and nitrogen stocks are expressed per unit
ground area. The main units are:

| Quantity | Symbol | Unit |
|:--|:--|:--|
| Carbon stock or daily flux | ``C``, ``F_C`` | g C m⁻²; g C m⁻² day⁻¹ |
| Nitrogen stock or daily flux | ``N``, ``F_N`` | g N m⁻²; g N m⁻² day⁻¹ |
| Water storage or daily flux | ``W``, ``F_W`` | mm; mm day⁻¹ |
| Temperature | ``T`` | °C |
| Radiation forcing | ``R`` | W m⁻² |
| Leaf-area index | ``L`` | m² leaf m⁻² ground |
| CO₂ partial pressure | ``p_a``, ``p_i`` | Pa |

The model distinguishes prognostic state (carried to the next day), fluxes
(integrated over the current day), auxiliary quantities (recomputed from state
and forcing), inputs, and discrete events. This distinction is important both
for restart correctness and for the future differentiable transition API.

## Daily coupled transition

For state ``x_d``, daily forcing ``u_d``, management ``m_d``, and parameters
``\vartheta``, one model day can be written abstractly as

```math
x_{d+1}=\mathcal{T}(x_d,u_d,m_d;\vartheta),
```

with diagnostics and boundary fluxes

```math
y_d=\mathcal{H}(x_d,u_d,m_d;\vartheta).
```

The operator is split into ordered physical and biological processes. It
updates climate memory and management, surface properties and snow, soil
hydraulic and thermal properties, existing soil C/N pools, crop phenology and
harvest, infiltration, canopy radiation, photosynthesis and water supply,
nitrogen uptake, respiration and allocation, evapotranspiration, and post-crop
nitrogen losses.

Order matters. For example, mineral nitrogen released from existing organic
matter is available for crop uptake on the same day. New harvest residue is
routed after decomposition of the existing litter, so it first decomposes on
the following day. The exact executable sequence is listed in
[Daily process order](../concepts/daily_processes.md).

## Solar geometry, radiation, and equilibrium evaporation

For day of year ``d`` and latitude ``\varphi``, solar declination is

```math
\delta=-23.4^\circ\cos\left(\frac{2\pi(d+10)}{365}\right).
```

Define

```math
u=\sin\varphi\sin\delta,\qquad
v=\cos\varphi\cos\delta.
```

Day length ``h_d`` in hours is

```math
h_d=\begin{cases}
24, & u\ge v,\\
0, & u\le -v,\\
\dfrac{24}{\pi}\cos^{-1}(-u/v), & \text{otherwise}.
\end{cases}
```

With downward shortwave radiation ``R_{sw}``, Agrocosm approximates daily PAR
energy as

```math
I_{PAR}=\frac{t_{day}R_{sw}}{2},\qquad t_{day}=86400\ \mathrm{s}.
```

Net shortwave radiation uses the effective surface albedo ``\alpha_s``:

```math
R_{sw,net}=(1-\alpha_s)R_{sw}.
```

The equilibrium evaporation diagnostic is

```math
E_{eq}=t_{day}\frac{s(T)}{s(T)+\gamma(T)}
\frac{R_{sw,net}+R_{lw,net}h_d/24}{\Lambda(T)},
```

where

```math
s(T)=\frac{2.503\times10^6\exp[17.269T/(237.3+T)]}{(237.3+T)^2},
```

```math
\gamma(T)=65.05+0.064T,\qquad
\Lambda(T)=2.495\times10^6-2380T.
```

For numerical protection, ``E_{eq}`` is clipped to ``[0,15]`` mm day⁻¹.
Actual canopy interception, transpiration, litter evaporation, and soil
evaporation are subsequently constrained separately.

## Surface albedo and absorbed PAR

Actual LAI accounts for carbon shortage:

```math
L_{act}=\max(0,L-L_{deficit}).
```

For most crops, the green-canopy fraction is Beer–Lambert absorption,

```math
f_g=1-\exp(-k_L L_{act}),
```

where ``k_L`` is the PFT light-extinction coefficient. Maize uses the LPJmL
empirical relationship

```math
f_g=\operatorname{clamp}(0.2558\max(0.01,L_{act})-0.0024,0,1).
```

Surface litter dry matter and cover are

```math
M_{lit}=\frac{C_{lit,surf}}{f_{C,lit}},\qquad
f_{lit}=1-\exp(-0.006M_{lit}),
```

with the current default ``f_{C,lit}=0.42``. The crop-covered surface albedo
combines mutually exclusive green, litter, and soil background fractions:

```math
\alpha_{crop}=f_g\alpha_g
+f_{lit}(1-f_g)\alpha_{lit}
+(1-f_{lit})(1-f_g)\alpha_{soil}^{*}.
```

Snow replaces leaf and litter albedo when snow is present and blends soil and
snow using snow-cover fraction. Outside the crop fractional cover ``f_{pc}``,
bare-surface albedo is added explicitly. With no active crop, only the bare
surface contributes.

Absorbed PAR is

```math
I_a=I_{PAR}(1-\alpha_{leaf})\alpha_a f_g.
```

When snow depth is positive, crop ``f_g`` and ``I_a`` are set to zero.

## Snow

Precipitation falls as snow when ``T<T_{snow}`` and is removed from liquid
precipitation. Snow water equivalent is capped at ``S_{max}``; excess is snow
runoff. The current scheme removes 0.1 mm day⁻¹ as sublimation whenever the
pack exceeds 0.1 mm.

For positive air temperature, conductive energy through a bounded snow skin
depth supplies melt energy. Melt water is

```math
M=\frac{\min(Q_{melt},\,D_s L_f)}{L_f},
```

where ``D_s`` is the active snow-water depth and ``L_f`` is latent heat of
fusion in the units used by the kernel. Melt is added to the day's liquid
precipitation before interception and infiltration. Snow height ``H_s`` and
fractional cover are

```math
H_s=c_{ws}\frac{S}{1000},\qquad
f_s=\frac{H_s}{H_s+0.5c_r}.
```

## Crop establishment and phenology

A stand is established on its prescribed sowing date if its cultivation
conditions are satisfied. During growth, daily heat units are

```math
HU_d=\max(0,T_d-T_b).
```

For winter crops, the vernalization increment is a piecewise-linear response:

```math
\Delta V_d=\begin{cases}
\dfrac{T-T_{v,e}^{lo}}{T_{v,o}^{lo}-T_{v,e}^{lo}},
&T_{v,e}^{lo}\le T<T_{v,o}^{lo},\\
1,&T_{v,o}^{lo}\le T<T_{v,o}^{hi},\\
\dfrac{T_{v,e}^{hi}-T}{T_{v,e}^{hi}-T_{v,o}^{hi}},
&T_{v,o}^{hi}\le T\le T_{v,e}^{hi},\\
0,&\text{otherwise}.
\end{cases}
```

Let ``V_{req}`` be the vernalization requirement and ``V_b=V_{req}/5``. The
vernalization reduction factor is

```math
f_V=\begin{cases}
0,&V<V_b,\\
\operatorname{clamp}\left(\dfrac{V-V_b}{V_{req}-V_b},0,1\right),
&V_b\le V<V_{req},\\
1,&V\ge V_{req}.
\end{cases}
```

Before senescence, the photoperiod factor is

```math
f_P=(1-p_{sens})\operatorname{clamp}
\left(\frac{h_d-p_b}{p_s-p_b},0,1\right)+p_{sens};
```

after senescence begins, ``f_P=1``. Accumulated heat units and fractional
phenological development are

```math
HU_\Sigma\leftarrow HU_\Sigma+HU_df_Vf_P,
\qquad
f_{PHU}=\min\left(1,\frac{HU_\Sigma}{PHU}\right).
```

Before the senescence threshold ``f_s``, the potential LAI fraction follows
the PFT logistic-like curve

```math
f_L=\frac{f_{PHU}}{f_{PHU}+c(c/k)^z},
```

with

```math
c=\frac{f_{PHU,c}}{f_{L,c}}-f_{PHU,c},\qquad
k=\frac{f_{PHU,k}}{f_{L,k}}-f_{PHU,k},\qquad
z=\frac{f_{PHU,c}-f_{PHU}}{f_{PHU,k}-f_{PHU,c}}.
```

During senescence,

```math
x_s=\left(\frac{1-f_{PHU}}{1-f_s}\right)^{q_s},\qquad
f_L=x_s(1-f_{L,harvest})+f_{L,harvest},
```

which is algebraically equivalent to the implemented interpolation between 1
and ``f_{L,harvest}``. Harvest is triggered at completed PHU or at the PFT
maximum growing duration.

Before senescence, potential LAI change is reduced by water and nitrogen:

```math
\Delta L=\Delta L_{pot}\min\left(\frac{f_W}{1.5},f_N\right).
```

During senescence, the model applies ``f_L`` to the LAI reached at the onset of
senescence rather than to the original PFT maximum.

## C3 photosynthesis

Temperature-adjusted Michaelis constants and Rubisco specificity are

```math
K_o(T)=K_{o,25}q_{10,o}^{(T-25)/10},\qquad
K_c(T)=K_{c,25}q_{10,c}^{(T-25)/10},
```

```math
\tau(T)=\tau_{25}q_{10,\tau}^{(T-25)/10},\qquad
\Gamma^*(T)=\frac{p_{O_2}}{2\tau(T)}.
```

Internal CO₂ partial pressure is represented by

```math
p_i=\lambda p_a.
```

The light- and Rubisco-response factors are

```math
c_1=f_T\alpha_{C3}\frac{p_i-\Gamma^*}{p_i+2\Gamma^*},\qquad
c_2=\frac{p_i-\Gamma^*}{p_i+K_c(1+p_{O_2}/K_o)}.
```

For absorbed PAR ``I_a``, daylight duration ``h_d``, carbon molar conversion
``M_C``, quantum conversion ``q_C``, and daily maximum carboxylation capacity
``V_{c\max}``, the two limiting rates are

```math
J_e=\frac{c_1I_aM_Cq_C}{h_d+\epsilon},\qquad
J_c=\frac{c_2V_{c\max}}{24}.
```

Agrocosm uses a non-rectangular hyperbola rather than a hard minimum:

```math
A_g=h_d\frac{J_e+J_c-
\sqrt{\max[0,(J_e+J_c)^2-4\theta J_eJ_c]}}{2\theta}.
```

Temperature-inactive cells and negative round-off are set to zero. Potential
``V_{c\max}`` is diagnosed at the configured reference internal-CO₂ ratio. The
optional nitrogen limitation of ``V_{c\max}`` is currently disabled by default.

## C4 photosynthesis

The C4 pathway uses the same co-limitation equation, with

```math
\phi_{p_i}=\min\left(1,\frac{\lambda}{\lambda_{m,C4}}\right),\qquad
c_1=f_T\alpha_{C4}\phi_{p_i},\qquad c_2=1.
```

Thus C4 assimilation saturates once the internal-CO₂ response reaches its
configured threshold. C3 and C4 use separate biochemical kernels but share the
same daily coupling, water limitation, respiration, allocation, and harvest
interfaces.

## Leaf dark respiration and conductance

Daily leaf dark respiration is

```math
R_d=bV_{c\max}.
```

The leaf respiration subtracted from daytime gross assimilation is scaled by
the daylight fraction:

```math
A_{dt}=\max\left(0,A_g-\frac{h_d}{24}R_d\right).
```

This distinction is important: the photosynthesis water-demand pathway uses
daytime dark respiration, while daily NPP subtracts the full daily ``R_d``.

At the potential internal-CO₂ ratio ``\lambda_{opt}``, canopy conductance is

```math
g_p=\frac{1.6A_{dt,mm}}
{p_a10^{-5}(1-\lambda_{opt})t_h}+g_{min}f_g,
```

where ``t_h`` is daylight duration in seconds and ``A_{dt,mm}`` is net
assimilation converted to the water-equivalent gas-volume quantity used by the
LPJmL conductance formulation.

## Transpiration and water limitation

Root-weighted relative soil water is

```math
w_r=\sum_l r_lw_l,
```

where ``r_l`` is the PFT root fraction and ``w_l`` is relative plant-available
water. Potential root water supply and atmospheric demand are

```math
S=e_{max}w_r\left(1-e^{-0.04C_{root}}\right),
```

```math
D=(1-f_{wet})E_{eq}\frac{\alpha_M}
{1+G_M\alpha_M/g_p}.
```

Potential transpiration is ``\min(S,D)`` multiplied by crop fractional cover.
Each layer's extraction is then capped by its liquid plant-available storage.
After this cap, actual conductance is recomputed when actual supply is below
demand:

```math
g_{actual}=G_M\alpha_M
\frac{E_{actual}}{(1-f_{wet})E_{eq}\alpha_M-E_{actual}}.
```

The seasonal diagnostic used in allocation is

```math
W_{deficit}=100\frac{\sum_d\min(S_d,D_d)}{\sum_dD_d},
```

bounded to 0–100% while the crop is active. Despite its historical name, a
larger value represents better water supply.

After layer-limited extraction, define

```math
g_{pd}=t_h(g_{actual}-g_{min}f_g),\qquad
f_g^*=\frac{g_{pd}}{1.6}p_a10^{-5}.
```

The realized internal-CO₂ ratio is the root of

```math
F(\lambda)=f_g^*(1-\lambda)-A_{dt,mm}(\lambda)=0.
```

Agrocosm solves this scalar equation with a bounded, fixed-maximum-iteration
bisection independently in every cell. The resulting ``\lambda`` is used for
the final photosynthesis calculation.

## Canopy interception and evaporation

Canopy interception storage fraction is

```math
f_{wet}=\operatorname{clamp}\left(
\frac{\min(i_cL_{act},0.9999)P}
{E_{eq}\alpha_{PT}},0,0.9999\right),
```

and interception evaporation is

```math
E_i=E_{eq}\alpha_{PT}f_{wet}f_{pc}.
```

The remaining atmospheric evaporation capacity after interception and
transpiration constrains litter and soil evaporation. Wet-litter evaporation is

```math
E_{lit}=\min\left[
E_{pot}w_{lit}^2f_{lit},\ W_{lit},\ E_{remaining}
\right],
```

where ``w_{lit}=W_{lit}/W_{lit,max}``. Potential soil evaporation is reduced by
canopy absorption and surface litter. Its moisture limitation is logistic:

```math
E_{soil}=\frac{E_{pot}}
{1+\exp(5-10W_{evap}/W_{hc,evap})}\max[0.05,(1-f_{lit})].
```

This flux is distributed over the configured evaporation depth and cannot
remove more liquid water than is available above wilting point.

## Nitrogen demand, uptake, and limitation

Leaf nitrogen demand combines Rubisco demand and minimum structural nitrogen:

```math
N_{leaf}^{dem}=p_N10^{-3}\frac{V_{c\max}}
{86400\times12\times10^{-6}}
\exp[-k_T(T-25)]+\rho_{N:C,leaf}^{min}C_{leaf}.
```

The resulting leaf N:C ratio is clipped to the PFT range. Total crop demand is

```math
N_{tot}^{dem}=N_{leaf}^{dem}+\rho_{N:C}
\left(\frac{C_{root}}{r_{root}}+
\frac{C_{pool}}{r_{pool}}+
\frac{C_{storage}}{r_{storage}}\right).
```

For mineral form ``j\in\{NO_3,NH_4\}`` in layer ``l``, potential uptake uses a
Michaelis–Menten response:

```math
U_{j,l}^{pot}=V_{max,j}\left(k_{min,j}+
\frac{N_{j,l}f_w}{N_{j,l}f_w+K_{m,j}W_{sat,l}D_l/1000}
\right)f_Tf_{NC}C_{root}\frac{r_l}{1000}.
```

Uptake is zero without liquid water, is temperature- and plant-N-status
limited, cannot exceed the layer mineral pool, and is globally scaled so that

```math
U_N=\min\left(\sum_{j,l}U_{j,l}^{pot},
\max[0,N_{tot}^{dem}-N_{plant}]\right).
```

Accepted uptake is removed proportionally from the contributing nitrate and
ammonium pools. With automatic fertilization, any remaining demand is supplied
as a boundary input and N sufficiency is set to one. Otherwise the crop N
sufficiency factor is derived from achieved versus optimal leaf demand.

## Fertilizer and manure

Prescribed fertilizer and manure are split between sowing and a second
application after ``f_{PHU}>0.25``. At application, mineral fertilizer is
partitioned as

```math
N_{fert}=f_{NO_3}N_{fert}\rightarrow NO_3
+(1-f_{NO_3})N_{fert}\rightarrow NH_4.
```

Manure sends fraction ``f_{NH_4}`` directly to top-layer ammonium. The organic
fraction enters incorporated litter nitrogen, accompanied by carbon determined
from the manure C:N ratio. Mineral-fertilizer configuration modes are:

- `no`: no prescribed or automatic mineral fertilizer;
- `yes`: use prescribed mineral-fertilizer input;
- `auto`: satisfy crop mineral-N demand automatically.

Manure is controlled by a separate Boolean option. When enabled, prescribed
manure is applied independently of the mineral-fertilizer mode.

## Autotrophic respiration and NPP

The acclimated temperature multiplier for air or topsoil temperature is

```math
g(T)=\exp\left[e_0\left(
\frac{1}{T_r+10}-\frac{1}{\min(T,40)+T_r}
\right)\right],
```

and is zero below −15 °C. Root maintenance respiration uses topsoil
temperature; storage and mobile-pool maintenance respiration use air
temperature:

```math
R_x=C_xr_{coeff}k\rho_{N:C,x}g(T_x).
```

Growth respiration is a fixed fraction of carbon remaining after full daily
leaf dark respiration and maintenance respiration:

```math
R_g=\max\left[0,
(A_g-R_d-R_{root}-R_{storage}-R_{pool})r_g\right].
```

Daily NPP is therefore

```math
NPP=A_g-R_d-(R_{root}+R_{storage}+R_{pool}+R_g).
```

Negative NPP is retained. This allows maintenance costs to reduce biomass on
days with little photosynthesis. If biomass plus NPP is non-positive, Agrocosm
follows the LPJmL failed-crop pathway and terminates the stand instead of
allowing negative living pools to persist.

## Carbon allocation

Living crop biomass closes over four pools:

```math
C_{bio}=C_{root}+C_{leaf}+C_{storage}+C_{pool}.
```

The stress metric used for root allocation is

```math
d_f=\min(W_{deficit},N_{sufficiency}),
```

and the root fraction is

```math
f_{root}=0.4-
0.3f_{PHU}\frac{d_f}{d_f+\exp(6.13-0.0883d_f)}.
```

Thus ``C_{root}=f_{root}C_{bio}``. Before senescence, leaf carbon is the lesser
of the LAI requirement and carbon remaining after roots:

```math
C_{leaf}=\min\left(\frac{L}{SLA},C_{bio}-C_{root}\right).
```

If carbon is insufficient to support potential LAI, the difference is stored
as ``L_{deficit}`` rather than changing the phenological LAI trajectory.

The potential harvest-index development factor is

```math
f_{HI}=\frac{100f_{PHU}}
{100f_{PHU}+\exp(11.1-10f_{PHU})}.
```

PFT optimal and minimum harvest indices are transformed according to whether
they are defined below or above one, then interpolated using the same
stress-response function. Storage carbon is computed from the resulting
harvest index and cannot exceed biomass remaining after leaf and root carbon.
The mobile pool is the exact residual:

```math
C_{pool}=C_{bio}-C_{leaf}-C_{root}-C_{storage}.
```

During senescence, a negative residual is repaired conservatively in the order
storage, root, then leaf, keeping living pools non-negative while preserving
total biomass.

## Harvest and failed crops

At normal harvest, storage carbon becomes yield. A configured fraction of leaf
and mobile-pool material is retained as surface or incorporated residue; root
material becomes root litter. The remainder crosses the model boundary as
harvest export. Nitrogen follows the corresponding organ-specific routing.

The failed-crop condition is

```math
C_{bio}+NPP\le10^{-4}
\quad\text{or}\quad
(L_{act}\le0\ \text{before senescence}).
```

A failed crop is terminated immediately. Existing storage is recorded as
yield, removable above-ground material is exported according to the residue
fraction, roots and retained residues enter litter, and all living crop state
and pending fertilizer/manure state are reset.

## Soil water

Each layer separates water below wilting point, plant-available water, free
water above field capacity, and corresponding ice pools. A commonly used
diagnostic is relative plant-available water

```math
w_l=\frac{W_l-W_{wilt,l}}{W_{fc,l}-W_{wilt,l}},
```

with storage-aware clipping in the implementation. Rain and snowmelt enter the
top boundary after canopy interception. Infiltration is processed in bounded
slugs for numerical stability. The layer update enforces

```math
0\le W_l\le W_{sat,l}
```

and partitions excess into vertical percolation and runoff. Water exceeding
the bottom-layer capacity exits as drainage. Nitrate is transported with
percolating water using the source-layer concentration and cannot exceed the
source nitrate pool.

Full irrigation, when enabled, restores each rooted soil layer to field
capacity every day. It is intentionally unconstrained by water availability;
the added water is a model boundary input. This option does not represent an
irrigation infrastructure or allocation model.

## Soil thermal state

Soil temperature is recovered from layer enthalpy rather than updated as an
independent diagnostic. A conceptual layer enthalpy is

```math
H_l=C_lT_l-L_fW_{ice,l},
```

where effective heat capacity ``C_l`` depends on mineral soil, organic matter,
liquid water, and ice. Thermal conductivity is reconstructed from the current
phase composition, and interface conductivity uses a harmonic mean:

```math
k_{l+1/2}=\frac{2k_lk_{l+1}}{k_l+k_{l+1}}.
```

Liquid water carries enthalpy between layers. Rain enters at air temperature,
snowmelt at 0 °C, evaporation and transpiration remove source-layer water
enthalpy, and bottom drainage exports source-layer enthalpy. After transfer,
temperature, liquid/ice partition, heat capacity, latent heat, and conductivity
are reconstructed. Consequently, saturation fraction and the previous-step
thermal properties belong to the numerical state lifecycle rather than being
pure daily output diagnostics.

## Soil decomposition response

The liquid pore-space moisture factor first diagnoses

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

The combined decomposition response is

```math
R_l=\operatorname{clamp}\left[
f_T(T_l)(a_0+a_1\omega_l+a_2\omega_l^2+a_3\omega_l^3),0,1\right].
```

Surface litter uses surface-litter temperature and wetness; incorporated and
root litter use the topsoil response.

## Soil carbon

For a pool stock ``C`` with rate ``k`` and environmental response ``R``, the
exact daily first-order loss is

```math
\Delta C=-\operatorname{expm1}(-kR)C
=\left(1-e^{-kR}\right)C.
```

The `expm1` form avoids cancellation for small Float32 decay rates. Litter loss
is divided between atmospheric respiration and retained material:

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
bioturbation provides the non-tillage vertical litter pathway. Harvest inputs
are routed only after the current day's decomposition.

## Soil nitrogen transformations

Organic litter, fast, and slow nitrogen use the same environmental response,
decay fractions, and vertical routing as their carbon counterparts. Gross
mineralization comprises atmospheric-fraction litter N release plus decomposed
fast and slow SOM nitrogen:

```math
M_l=M_{lit,l}+\Delta N_{fast,l}+\Delta N_{slow,l},
```

and enters ammonium. If retained litter carbon requires more N than decomposed
litter supplies,

```math
N_{def}=\frac{\Delta C_{lit}}{(C:N)_{soil}}-\Delta N_{lit}>0,
```

the model immobilizes mineral nitrogen into fast and slow organic pools. The
immobilization response to mineral-N concentration ``c_N`` is

```math
f_{immob}=\frac{c_N}{K_{immob}+c_N}.
```

Removal is shared between ammonium and nitrate in proportion to their
availability and cannot make either pool negative.

Nitrification converts ammonium to nitrate using pH, temperature, and
water-filled pore-space responses. A configured fraction becomes N₂O; the
remainder enters nitrate. Denitrification acts after crop uptake and is limited
by nitrate, available soil carbon, temperature, pH, and anaerobic moisture
response. Its products are partitioned between N₂O and N₂.

Surface NH₃ volatilization acts on the **entire top-layer ammonium pool each
day**, not only on fertilizer applied that day. Its response uses top-layer pH,
temperature, and wind speed. This is why cumulative volatilization can remain
large when fertilization and mineralization continuously replenish ammonium.

The daily mineral-N sequence is:

1. organic N mineralization and immobilization;
2. nitrification;
3. fertilizer/manure application and crop uptake according to daily ordering;
4. nitrate transport with water;
5. denitrification and NH₃ volatilization.

## Conservation diagnostics

For a conserved quantity ``Q`` over one day, Agrocosm evaluates the residual

```math
\varepsilon_Q=
Q_{d+1}-Q_d-\sum F_{in}+\sum F_{out}.
```

Carbon inputs and outputs include assimilation, harvest export, autotrophic
respiration, and heterotrophic respiration. Nitrogen includes prescribed or
automatic fertilizer/manure, harvest export, leaching, N₂, N₂O, and NH₃.
Water includes precipitation, irrigation, runoff, drainage, interception,
transpiration, and evaporation. Energy closure includes surface forcing and
advective enthalpy carried by water.

Absolute residuals can look nonzero for large energy stocks in Float32, so
both absolute and scale-aware relative residuals should be inspected. CPU and
accelerator equivalence is a separate property from physical closure: two
backends can agree while sharing the same missing boundary term.

## Numerical execution

The reference and backend implementations use the same equations. Cell-local
scalar operations are launched through KernelAbstractions. One work item owns
one horizontal grid cell and loops over the fixed soil column where vertical
dependence exists. This avoids races between layers of the same column and
keeps the same process implementation available across supported array
backends.

The main numerical safeguards are:

- non-negative square-root discriminants in photosynthesis;
- bounded bisection iterations for ``\lambda``;
- exact exponential decay with `expm1`;
- storage-limited water, carbon, and nitrogen withdrawals;
- explicit clipping of fractions and environmental response functions;
- bounded infiltration iterations;
- failed-crop termination before negative living biomass persists.

I/O, agricultural warm-up orchestration, checkpoints, and streamed output stay
outside the daily numerical kernels. Agricultural warm-up is a preprocessing
stage and is not required to enter the future Enzyme-differentiable production
transition.

## Implementation map

| Process | Primary implementation |
|:--|:--|
| Radiation, day length, APAR | `src/processes/crop/radiation.jl` |
| Surface albedo | `src/processes/crop/albedo.jl` |
| Snow | `src/processes/climate/snow.jl` |
| Phenology and LAI | `src/processes/crop/phenology.jl`, `lai_crop.jl` |
| C3/C4 photosynthesis | `src/processes/crop/photosynthesis.jl` |
| Water supply and transpiration | `src/processes/crop/transpiration.jl` |
| Internal CO₂ solver | `src/processes/crop/lambda_solver.jl` |
| Respiration and allocation | `src/processes/crop/respiration.jl`, `carbon_allocation.jl` |
| Crop nitrogen | `src/processes/crop/nitrogen_*.jl` |
| Fertilizer and manure | `src/processes/crop/fertilizer.jl` |
| Harvest and failure | `src/processes/crop/harvesting.jl` |
| Infiltration and percolation | `src/processes/soil/infil_perc.jl` |
| Soil evaporation | `src/processes/soil/evaporation.jl` |
| Soil temperature and enthalpy | `src/processes/soil/soil_temp.jl` |
| Soil C decomposition | `src/processes/soil/soil_carbon.jl` |
| Soil N transformations | `src/processes/soil/nitrogen_transform.jl` |

Parameter values are intentionally not duplicated on this page because they
depend on crop PFT and model configuration. The symbols above map directly to
the corresponding `PftParameters`, `LPJmLParams`, `PhotoParams`, snow, and soil
decomposition parameter fields.
