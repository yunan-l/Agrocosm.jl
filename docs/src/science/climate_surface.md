# Climate and surface processes

This page documents the radiation, equilibrium evaporation, surface albedo,
snow, and rolling climate-memory terms that feed the crop and soil processes.

## Solar geometry and day length

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
24,&u\ge v,\\
0,&u\le-v,\\
\dfrac{24}{\pi}\cos^{-1}(-u/v),&\text{otherwise}.
\end{cases}
```

The explicit polar-day and polar-night branches prevent invalid inverse-cosine
arguments.

## Radiation and equilibrium evaporation

With downward shortwave radiation ``R_{sw}``, daily PAR energy is approximated
as

```math
I_{PAR}=\frac{t_{day}R_{sw}}{2},\qquad t_{day}=86400\ \mathrm{s}.
```

Net shortwave radiation uses effective surface albedo ``\alpha_s``:

```math
R_{sw,net}=(1-\alpha_s)R_{sw}.
```

The equilibrium evaporation diagnostic is

```math
E_{eq}=t_{day}\frac{s(T)}{s(T)+\gamma(T)}
\frac{R_{sw,net}+R_{lw,net}h_d/24}{\Lambda(T)},
```

with

```math
s(T)=\frac{2.503\times10^6\exp[17.269T/(237.3+T)]}{(237.3+T)^2},
```

```math
\gamma(T)=65.05+0.064T,\qquad
\Lambda(T)=2.495\times10^6-2380T.
```

For numerical protection, ``E_{eq}`` is clipped to 0--15 mm day⁻¹. Canopy
interception, transpiration, litter evaporation, and soil evaporation are then
limited separately.

## Surface albedo

Surface litter dry mass and cover are

```math
M_{lit}=\frac{C_{lit,surf}}{f_{C,lit}},\qquad
f_{lit}=1-\exp(-0.006M_{lit}),
```

with current default ``f_{C,lit}=0.42``. For green fraction ``f_g``, the
crop-covered surface albedo is

```math
\alpha_{crop}=f_g\alpha_g
+f_{lit}(1-f_g)\alpha_{lit}
+(1-f_{lit})(1-f_g)\alpha_{soil}^{*}.
```

Snow replaces leaf and litter albedo when present and blends the soil and snow
background through snow-cover fraction. Outside crop fractional cover
``f_{pc}``, bare-surface albedo is added explicitly. With no active crop, only
the bare soil/snow mixture contributes.

## Snow accumulation and melt

Precipitation falls as snow when ``T<T_{snow}`` and is removed from liquid
precipitation. Snow water equivalent is capped at ``S_{max}``; excess becomes
snow runoff. The current scheme removes 0.1 mm day⁻¹ as sublimation whenever
the pack exceeds 0.1 mm.

For positive air temperature, conduction through a bounded snow-skin depth
supplies melt energy. Melt water is

```math
M=\frac{\min(Q_{melt},D_sL_f)}{L_f},
```

where ``D_s`` is active snow-water depth and ``L_f`` latent heat of fusion in
the kernel's units. Melt joins liquid precipitation before interception and
infiltration. Snow height and fractional cover are

```math
H_s=c_{ws}\frac{S}{1000},\qquad
f_s=\frac{H_s}{H_s+0.5c_r}.
```

## Climate memory

The rolling climate buffers retain recent temperature history used by
phenology and temperature diagnostics. Their updates are part of the daily
state transition. Annual extrema and summaries are computed by cell-local
kernels at the year boundary, so each work item owns all mutations for its
grid cell and no host-side array sort participates in backend execution.

These histories influence later days and are therefore lifecycle state, not
disposable daily output.

## Forcing contract

The current default process pathway requires daily mean temperature,
precipitation, net longwave radiation, downward shortwave radiation, and annual
CO₂. Production input is normalized to a strict 365-day calendar before the
process kernels run. Wind enters the NH₃ volatilization input path. Humidity,
VPD, and pressure will only become mandatory if the optional
Farquhar--Medlyn--Penman--Monteith pathway is implemented.

## Code map

- `src/processes/crop/radiation.jl`
- `src/processes/crop/albedo.jl`
- `src/processes/climate/snow.jl`
- `src/processes/climate/climbuf.jl`
- `src/processes/climate/temp_stress.jl`
- `src/processes/climate/readclimate.jl`
