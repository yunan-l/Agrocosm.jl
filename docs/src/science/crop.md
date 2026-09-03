# Crop processes

This page documents the equations used by the current default crop pathway.
Parameter values depend on the selected CFT; symbols map to fields in the
model parameter structures rather than defining a second parameter source.

## Establishment, phenology, and LAI

By default (`sowing_mode=:prescribed_sdate`), a stand is established on its
prescribed sowing date when its cultivation conditions are satisfied. The
optional `sowing_mode=:dynamic_sdate` resolves an annual sowing month from
the retained climate history and then emits a one-day establishment event:
temperature-threshold crossing for wheat and maize, or the first wet day of
the resolved rainy-season month for rice and soybean (with a month-end
fallback). The prescribed date remains the deterministic bootstrap month until
the precipitation climatology is available. Daily heat units are

```math
HU_d=\max(0,T_d-T_b).
```

For winter crops, the vernalization increment is

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

With vernalization requirement ``V_{req}`` and ``V_b=V_{req}/5``, the
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

```math
f_L=\frac{f_{PHU}}{f_{PHU}+c(c/k)^z},
```

where

```math
c=\frac{f_{PHU,c}}{f_{L,c}}-f_{PHU,c},\qquad
k=\frac{f_{PHU,k}}{f_{L,k}}-f_{PHU,k},\qquad
z=\frac{f_{PHU,c}-f_{PHU}}{f_{PHU,k}-f_{PHU,c}}.
```

During senescence,

```math
x_s=\left(\frac{1-f_{PHU}}{1-f_s}\right)^{q_s},\qquad
f_L=x_s(1-f_{L,harvest})+f_{L,harvest}.
```

Before senescence, water and nitrogen reduce potential LAI growth:

```math
\Delta L=\Delta L_{pot}\min\left(\frac{f_W}{1.5},f_N\right).
```

During senescence the model scales from LAI at senescence onset. Harvest is
triggered by completed PHU or the CFT maximum growing duration.

## Canopy cover and absorbed PAR

Carbon shortage reduces actual LAI:

```math
L_{act}=\max(0,L-L_{deficit}).
```

For most crops, the green fraction follows Beer--Lambert absorption,

```math
f_g=1-\exp(-k_LL_{act}),
```

where ``k_L`` is the CFT extinction coefficient. Maize uses

```math
f_g=\operatorname{clamp}(0.2558\max(0.01,L_{act})-0.0024,0,1).
```

Absorbed PAR is

```math
I_a=I_{PAR}(1-\alpha_{leaf})\alpha_af_g.
```

Snow sets crop ``f_g`` and ``I_a`` to zero. The complete surface-albedo
mixture is documented in [Climate and surface processes](climate_surface.md).

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

For absorbed PAR ``I_a``, day length ``h_d``, carbon molar conversion ``M_C``,
quantum conversion ``q_C``, and daily maximum carboxylation capacity
``V_{c\max}``, the limiting rates are

```math
J_e=\frac{c_1I_aM_Cq_C}{h_d+\epsilon},\qquad
J_c=\frac{c_2V_{c\max}}{24}.
```

Gross assimilation is their non-rectangular-hyperbola co-limitation:

```math
A_g=h_d\frac{J_e+J_c-
\sqrt{\max[0,(J_e+J_c)^2-4\theta J_eJ_c]}}{2\theta}.
```

Temperature-inactive cells and negative round-off are set to zero. Potential
``V_{c\max}`` is diagnosed at the configured reference internal-CO₂ ratio.
Nitrogen limitation of ``V_{c\max}`` is disabled by default.

## C4 photosynthesis

The C4 pathway uses the same co-limitation equation, with

```math
\phi_{p_i}=\min\left(1,\frac{\lambda}{\lambda_{m,C4}}\right),\qquad
c_1=f_T\alpha_{C4}\phi_{p_i},\qquad c_2=1.
```

Assimilation therefore saturates when the internal-CO₂ response reaches its
configured threshold. C3 and C4 share the daily water limitation, respiration,
allocation, and harvest interfaces.

## Leaf dark respiration and conductance

Daily leaf dark respiration is

```math
R_d=bV_{c\max}.
```

The photosynthesis water-demand pathway subtracts daylight-scaled dark
respiration,

```math
A_{dt}=\max\left(0,A_g-\frac{h_d}{24}R_d\right),
```

whereas daily NPP subtracts the full daily ``R_d``. At the potential
internal-CO₂ ratio ``\lambda_{opt}``, canopy conductance is

```math
g_p=\frac{1.6A_{dt,mm}}
{p_a10^{-5}(1-\lambda_{opt})t_h}+g_{min}f_g,
```

where ``t_h`` is daylight duration in seconds and ``A_{dt,mm}`` is net
assimilation converted to the gas-volume quantity used by the conductance
formulation.

## Transpiration and water limitation

Root-weighted relative soil water is

```math
w_r=\sum_l r_lw_l,
```

where ``r_l`` is root fraction and ``w_l`` relative plant-available water.
Potential root supply and atmospheric demand are

```math
S=e_{max}w_r\left(1-e^{-0.04C_{root}}\right),
```

```math
D=(1-f_{wet})E_{eq}\frac{\alpha_M}
{1+G_M\alpha_M/g_p}.
```

Potential transpiration is ``\min(S,D)`` times crop fractional cover. Each
layer's extraction is capped by its liquid plant-available storage. If actual
supply is below demand, conductance is recomputed as

```math
g_{actual}=G_M\alpha_M
\frac{E_{actual}}{(1-f_{wet})E_{eq}\alpha_M-E_{actual}}.
```

The seasonal diagnostic used in allocation is

```math
W_{deficit}=100\frac{\sum_d\min(S_d,D_d)}{\sum_dD_d},
```

bounded to 0--100% while the crop is active. Despite its historical name, a
larger value means better water supply.

After layer-limited extraction,

```math
g_{pd}=t_h(g_{actual}-g_{min}f_g),\qquad
f_g^*=\frac{g_{pd}}{1.6}p_a10^{-5}.
```

The realized internal-CO₂ ratio is the root of

```math
F(\lambda)=f_g^*(1-\lambda)-A_{dt,mm}(\lambda)=0.
```

Agrocosm solves this equation by bounded, fixed-maximum-iteration bisection in
each cell. The resulting ``\lambda`` is used for final photosynthesis.

## Interception and evaporation

Canopy wetness is

```math
f_{wet}=\operatorname{clamp}\left(
\frac{\min(i_cL_{act},0.9999)P}
{E_{eq}\alpha_{PT}},0,0.9999\right),
```

and interception evaporation is

```math
E_i=E_{eq}\alpha_{PT}f_{wet}f_{pc}.
```

The remaining atmospheric capacity constrains litter and soil evaporation.
Wet-litter evaporation is

```math
E_{lit}=\min\left[E_{pot}w_{lit}^2f_{lit},W_{lit},E_{remaining}\right].
```

Potential soil evaporation is moisture-limited by

```math
E_{soil}=\frac{E_{pot}}
{1+\exp(5-10W_{evap}/W_{hc,evap})}\max[0.05,(1-f_{lit})],
```

distributed over the configured evaporation depth and capped by liquid water
above wilting point.

## Nitrogen demand and uptake

Leaf nitrogen demand combines Rubisco demand and structural nitrogen:

```math
N_{leaf}^{dem}=p_N10^{-3}\frac{V_{c\max}}
{86400\times12\times10^{-6}}
\exp[-k_T(T-25)]+\rho_{N:C,leaf}^{min}C_{leaf}.
```

After clipping leaf N:C to the CFT range, total demand is

```math
N_{tot}^{dem}=N_{leaf}^{dem}+\rho_{N:C}
\left(\frac{C_{root}}{r_{root}}+
\frac{C_{pool}}{r_{pool}}+
\frac{C_{storage}}{r_{storage}}\right).
```

For mineral form ``j\in\{NO_3,NH_4\}`` and layer ``l``, potential uptake is

```math
U_{j,l}^{pot}=V_{max,j}\left(k_{min,j}+
\frac{N_{j,l}f_w}{N_{j,l}f_w+K_{m,j}W_{sat,l}D_l/1000}
\right)f_Tf_{NC}C_{root}\frac{r_l}{1000}.
```

Uptake is water-, temperature-, root-, plant-status-, and pool-limited, then
scaled so that

```math
U_N=\min\left(\sum_{j,l}U_{j,l}^{pot},
\max[0,N_{tot}^{dem}-N_{plant}]\right).
```

Accepted uptake is removed proportionally from nitrate and ammonium. Automatic
fertilization supplies remaining demand as a boundary input; otherwise the
achieved leaf demand determines nitrogen sufficiency.

## Fertilizer and manure

Prescribed fertilizer and manure are split between sowing and a second
application after ``f_{PHU}>0.25``. Mineral fertilizer is partitioned as

```math
N_{fert}=f_{NO_3}N_{fert}\rightarrow NO_3
+(1-f_{NO_3})N_{fert}\rightarrow NH_4.
```

Manure sends fraction ``f_{NH_4}`` to top-layer ammonium. Its organic fraction
enters incorporated litter N with carbon determined by the manure C:N ratio.
Mineral-fertilizer modes are `no`, `yes`, and `auto`; manure is controlled
independently.

## Autotrophic respiration and NPP

The acclimated temperature multiplier is

```math
g(T)=\exp\left[e_0\left(
\frac{1}{T_r+10}-\frac{1}{\min(T,40)+T_r}
\right)\right],
```

and is zero below -15 °C. Root maintenance respiration uses topsoil
temperature; storage and mobile-pool respiration use air temperature:

```math
R_x=C_xr_{coeff}k\rho_{N:C,x}g(T_x).
```

Growth respiration is a fixed fraction of carbon remaining after leaf dark
respiration and maintenance respiration:

```math
R_g=\max\left[0,
(A_g-R_d-R_{root}-R_{storage}-R_{pool})r_g\right].
```

Daily NPP is

```math
NPP=A_g-R_d-(R_{root}+R_{storage}+R_{pool}+R_g).
```

Negative NPP is retained so maintenance can reduce biomass on low-assimilation
days. A stand is terminated before negative living pools can persist.

## Carbon allocation

Living biomass closes over four pools:

```math
C_{bio}=C_{root}+C_{leaf}+C_{storage}+C_{pool}.
```

With stress metric ``d_f=\min(W_{deficit},N_{sufficiency})``, root fraction is

```math
f_{root}=0.4-
0.3f_{PHU}\frac{d_f}{d_f+\exp(6.13-0.0883d_f)}.
```

Thus ``C_{root}=f_{root}C_{bio}``. Before senescence,

```math
C_{leaf}=\min\left(\frac{L}{SLA},C_{bio}-C_{root}\right).
```

Carbon insufficient for potential LAI is recorded in ``L_{deficit}``. Harvest
index develops according to

```math
f_{HI}=\frac{100f_{PHU}}
{100f_{PHU}+\exp(11.1-10f_{PHU})}.
```

Storage is limited by biomass remaining after leaf and root carbon. The mobile
pool closes the balance exactly:

```math
C_{pool}=C_{bio}-C_{leaf}-C_{root}-C_{storage}.
```

During senescence a negative residual is repaired conservatively in storage,
root, then leaf order while preserving total biomass.

## Harvest and failed crops

At normal harvest, storage carbon becomes yield. The configured residue
fraction sends leaf and mobile-pool material to surface or incorporated litter;
root material becomes root litter. Remaining material crosses the boundary as
harvest export, with nitrogen following organ-specific routing.

The model output `crop_yield` is harvested storage-organ carbon in
``\mathrm{gC\,m^{-2}}``, not dry matter. For LPJmL/GGCMI-compatible reporting,
convert every crop with the common LPJmL carbon fraction of 0.45:

```math
Y_{\mathrm{tDM\,ha^{-1}}}=\frac{Y_{\mathrm{gC\,m^{-2}}}}{0.45}\times0.01.
```

The factor applies uniformly across CFTs for this reporting convention. It is
distinct from `CCpDM = 0.4763`, which describes leaf carbon per unit leaf dry
matter and must not be used to convert harvested storage-organ yield.

The failed-crop condition is

```math
C_{bio}+NPP\le10^{-4}
\quad\text{or}\quad
(L_{act}\le0\ \text{before senescence}).
```

Failure immediately records existing storage as yield, routes retained roots
and residues to litter, exports removable above-ground material, and resets
living crop and pending management state.

## Code map

- `src/processes/crop/radiation.jl`
- `src/processes/crop/albedo.jl`
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
