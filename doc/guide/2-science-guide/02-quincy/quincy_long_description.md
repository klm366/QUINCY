---
title: QUINCY model documentation
---

# Detailed model description

This appendix includes a detailed model description with equations.
Section <a href="#sec:structure" data-reference-type="ref"
data-reference="sec:structure">1.1</a> with parameter
Tab <a href="#tab:processes:mavg_period" data-reference-type="ref"
data-reference="tab:processes:mavg_period">[tab:processes:mavg_period]</a>
summarises the general structure and vertical discretisation of
vegetation and soil. Section
<a href="#ch:assimi" data-reference-type="ref"
data-reference="ch:assimi">1.2</a> describes the canopy processes, such
as photosynthesis and stomatal coupling, with parameters in
Tab. <a href="#tab:assimi:param" data-reference-type="ref"
data-reference="tab:assimi:param">[tab:assimi:param]</a>. Section
(<a href="#sec:veg" data-reference-type="ref"
data-reference="sec:veg">1.3</a>) introduces vegetation growth, turnover
and dynamics and the corresponding parameters are in Table
<a href="#tab:veg:param" data-reference-type="ref"
data-reference="tab:veg:param">[tab:veg:param]</a>. The soil
biochemistry is described in
Section <a href="#sec:soil_biogeochemistry" data-reference-type="ref"
data-reference="sec:soil_biogeochemistry">1.4</a>), and its parameters
are in Tab. <a href="#tab:sb:param" data-reference-type="ref"
data-reference="tab:sb:param">[tab:sb:param]</a>. Section
 <a href="#sec:rad" data-reference-type="ref"
data-reference="sec:rad">1.6</a> describes the radiation scheme, surface
energy balance and soil hydrology, with parameters described in
Tab. <a href="#tab:radsse:param" data-reference-type="ref"
data-reference="tab:radsse:param">[tab:radsse:param]</a>. The
PFT-specific parameters are in the Table
<a href="#tab:lctlib:description" data-reference-type="ref"
data-reference="tab:lctlib:description">[tab:lctlib:description]</a>.

## General model structure and discretisation

Each gridcell of the models is subdivided into nested tiles, each of
which is occupied by one specific vegetation-type, representing a plant
functional type (PFT). The number of tiles per gridcell is flexible,
making it is easy to implement more/different PFTs in the future. In the
model, vegetation is represented by an average individual composed of a
range of structural pools (leaves, sapwood, heartwood, coarse roots,
fine roots, and fruit), a fast overturning, respiring non-structural
pool (labile), as well as a seasonal, non-respiring, and non-structural
storage pool (reserve). Tree vegetation types are furthermore
characterised by their height (m), diameter (m), and stand density
(m<sup>−2</sup>). Soil biogeochemistry is represented using five organic
pools: metabolic (met), structural (str) and and woody (wl) litter, as
well as fast (f) and slow (s) overturning soil organic matter. Each of
these pools contains carbon (C), nitrogen (N) and phosphorus (P), as
well as , , and . The unit of the pools is mol   X   m<sup>−2</sup> for
vegetation and mol   X   m<sup>−3</sup> for soil biogeochemical pools,
where $X$ represents any of these elements. In addition, the model
represents the following soil biogeochemical pools (, , , , , and ),
with equivalent units.

The model operates on a half-hourly time-scale (denoted as $dt$).
Vegetation processes are assumed to respond to these instantaneous
conditions and associated fluxes with a process-specific lag time
($\tau_{mavg}^{process}$, see Tab.
<a href="#tab:processes:mavg_period" data-reference-type="ref"
data-reference="tab:processes:mavg_period">[tab:processes:mavg_period]</a>),
representing a form of memory for instance in the calculation of
allocation or vegetation dynamic responses. Where appropriate, the
fluxes or pool sizes are calculated as running means with a
time-averaging filter as

$$\begin{align}
X_{mavg,new}^{process} &= X_{mavg,old}^{process} \times (1 - \omega) + X_{current} \times \omega, where\\
&\omega = \frac{dt}{\tau_{mavg}^{process}}
\end{align}$$

where $X_{current}$ is the instantaneous state or flux of interest, and
$X_{mavg,old}^{process}$, as well as $X_{mavg,new}^{process}$ are the
averaged values of the previous and current time step, respectively.

### Vertical discretisation

The canopy is discretised into 10 layers (denoted by subscript $cl$),
with exponentially increasing layer depth ($LAI_{cl}$) to allow for a
better resolution of canopy processes in areas of high light and
nitrogen gradients. Total foliage N ($N_{leaf}$) is distributed to each
canopy layer $cl$, following $$\label{eq:veg:canopy_n_distribution}
\begin{align}
N_{leaf,cl} &= N_{leaf,cl=1} \times e^{-k_n \times LAI_c}, where \\
&N_{leaf,cl=1} = \frac{k_n}{1-e^{-k_n \times LAI}} N_{leaf}  
\end{align}$$

where $LAI_c$ is the cumulative leaf area above the centre-point of the
canopy layer, and $LAI$ the total leaf area, such that less N is
allocated to the lower, darker canopy layers, in accordance with
observations of canopy N distribution (Niinemets, Kull, and Tenhunen
1998).

All soil state-variables (temperature, moisture, texture, soil
biogeochemical pools) and fluxes are discretised into 15 soil layers
(denoted by subscript ($sl$). Layer thickness increases exponentially
with increasing layer depth up to a total depth of 9.5 m, and with a
minimum layer thickness for the top layer of 0.065 m. Fine roots and
coarse roots are assumed to be distributed in exponentially decreasing
density along the soil profile according to

$$\label{eq:veg:root_distribution}
\begin{align}
X_{root,sl} &= X_{root,sl=1} \times e^{-k_{rd} \times depth_{sl}}, where \\
&X_{root,sl=1} = \frac{k_{rd}}{1-e^{-k_{rd} \times depth_r}} X_{root} 
\end{align}$$

where $k_{rd}$ is a PFT-specific parameter (see Tab.
<a href="#tab:lctlib:description" data-reference-type="ref"
data-reference="tab:lctlib:description">[tab:lctlib:description]</a>,
Jackson et al. (1996)), $depth_{sl}$ the depth of the soil layer’s mid
point, $depth_r$ the rooting depth, derived from site characteristics
and $X_{root}$ the respective fine or coarse root mass. For clarity in
the following, the subscript for canopy and soil layer is ignored if
processes are treated similarly across layers.

## Canopy processes

### Canopy nitrogen allocation

The leaf area index ($LAI$) and canopy nitrogen content ($N_{leaf}$) are
dynamic properties of the model, as described in Section
<a href="#sec:veg" data-reference-type="ref"
data-reference="sec:veg">1.3</a>, and are discretised to canopy layers
given by Eq.
<a href="#eq:veg:canopy_n_distribution" data-reference-type="ref"
data-reference="eq:veg:canopy_n_distribution">[eq:veg:canopy_n_distribution]</a>.
$N_{leaf,cl}$ is partitioned into photosynthetic and non-photosynthetic,
or structural, N. The fraction of structural N ($fN_{struc,cl}$) is
calculated as a function of the total leaf N in the respective canopy
layer:

$$\label{eq:veg:struc_n}
fN_{struc,cl}=k^{struc}_0-k^{struc}_1 N_{leaf,cl}$$

where $k^{struc}_0$ and $k^{struc}_1$ are empirical parameters.

The photosynthetic N is further separated into the fraction that is
associated with Rubisco ($fN_{rub}$), electron transport ($fN_{et}$),
chlorophyll ($fN_{chl}$) and in the case of C4 plants, a fourth fraction
for PEP carboxylase $fN_{pepc}$.

The fraction of N in chlorophyll for each layer is calculated as
decreasing with canopy depth: $$\label{eq:canopy_chl}
fN_{chl}=\frac{k_0^{chl} - k_1^{chl} e^{-k_{fn}^{chl} LAI_c}} {a_{chl}^n},$$

where $k_0^{chl}$, $k_1^{chl}$ and $k^{chl}$ are empirical parameters,
$a_{chl}^n$ is the molecular N content of chlorophyll.

The values of $fN_{rub}$ and $fN_{et}$ are calculated assuming a fixed
ratio of the $V_{cmax}$ and $J_{max}$ photosynthetic parameters at 25°C,
$r_{J2V}$, given the calculated values of the structural and
photosynthetic fractions. The PEP carboxylase fraction, $fN_{pep}$, is
considered to be a constant.

### Leaf-level net photosynthesis

The calculation of leaf-level photosynthesis is based on (Kull and
Kruijt 1998), extended for C4 photosynthesis according to (Friend et al.
2009), and with the temperature response curves as described in
(Bernacchi et al. 2001). Unless stated otherwise, temperature
sensitivities follow the form:

$$\label{eq:ps_t_sensitivity}
f_x(T_{air}) = e^{E_0^x - E_1^x/(R\times T_{air})}$$

where $T_{air}$ is the air temperature (K), $R$ is the universal gas
constant (J mol<sup>−1</sup> K<sup>−1</sup>), and the process-wise
$E_0^x$ and $E_1^x$ are given in Table
<a href="#tab:assimi:param" data-reference-type="ref"
data-reference="tab:assimi:param">[tab:assimi:param]</a>.

Photosynthesis and stomatal conductance are calculated for the
mid-points of each canopy layer and light-quality class (sunlit and
shaded; as defined in Sect.
<a href="#sec:rad:vegetation" data-reference-type="ref"
data-reference="sec:rad:vegetation">1.6.1</a>). For clarity, the
subscript $_{cl}$ is omitted in this section.

In light-saturated conditions, gross photosynthesis ($A_g$) in C3 plants
is calculated as the minimum of two potential rates, The
electron-transport capacity limited carboxylation ($A_j$) and the
Rubisco-limited rate of photosynthesis ($A_v$). $A_j$ is given by
$$\label{eq:jmax}
\begin{align}
A_j=& m_1 \times J_{max}, where \\
&m_1=\frac{c_i}{c_i+2\times \Gamma^{*}}, \\
&J_{max} = n_1 \times N_{leaf},  \\
&n_1 = g_{jmax}(T) \times \beta_{soil}^{ps} \times \beta_{sinklim}^{ps} \times  \beta_{soa}^{ps} \times j_{max}^n \times fN_{et}
\end{align}$$

where $c_i$ is the intercellular partial pressure of (Pa,
Eq. <a href="#eq:stomatal_coupling" data-reference-type="ref"
data-reference="eq:stomatal_coupling">[eq:stomatal_coupling]</a>), and
$\Gamma^{*}$ is the compensation point in the absence of dark
respiration. $\beta_{sinklim}^{ps}$ is a signal to reduce photosynthesis
in the case of C sink limitation
(Eq. <a href="#eq:sinklim_c" data-reference-type="ref"
data-reference="eq:sinklim_c">[eq:sinklim_c]</a>) and $\beta_{soa}^{ps}$
is accounting for the effect of low-temperature acclimation in the
evergreen species (Eq. <a href="#eq:veg:soa" data-reference-type="ref"
data-reference="eq:veg:soa">[eq:veg:soa]</a>). Excessive soil moisture
constraints (as discussed in Rogers et al. (2017)) are assumed to reduce
light-saturated photosynthetic activity by:
$$\beta_{soil}^{ps} = 1 - \frac{\Psi_{soil}}{\Psi_{leaf}^{min}},$$

where $\Psi_{soil}$ is the soil water potential in the root zone
(Eq. <a href="#eq:sse:soil_matrix_potential" data-reference-type="ref"
data-reference="eq:sse:soil_matrix_potential">[eq:sse:soil_matrix_potential]</a>)
and $\Psi_{leaf}^{min}$ is the PFT-specific minimum leaf water
potential. The temperature sensitivity of electron transport is assumed
to follow the bell-shaped form described by (June, Evans, and Farquhar
2004), where $T_{jmax}^{opt}$ is the optimum temperature for $J_{max}$
according to (Friend 2010), as follows:

$$\label{eq:tjmax}
\begin{align}
g_{jmax}(T_{air}) =& e^{-\left({\frac{T_{air}-\bar{T}_{jmax}^{opt}}{T_{\Omega}}}\right)^2},\\
&T_{jmax}^{opt} = k0_{jmax}^{topt} + k1_{jmax}^{topt} \times T_{air}; T_{jmax,min}^{opt} < T_{jmax}^{opt} \le T_{jmax,max}^{opt},
\end{align}$$

where $T_{\Omega}$ is a PFT-specific parameter, $k0_{jmax}^{topt}$ and
$k1_{jmax}^{topt}$ are parameters, $T_{air}$ is the air temperature, and
$\bar{T}_{jmax}^{opt}$ in Eq.
<a href="#eq:tjmax" data-reference-type="ref"
data-reference="eq:tjmax">[eq:tjmax]</a> is the mean of the daytime
$T_{jmax}^{opt}$ over the past few days ($\tau_{mavg}^{jmax}$), thereby
accounting for temperature acclimation of photosynthesis as in (Friend
2010).

$A_v$ is given by $$\label{eq:vcmax}
\begin{align}
A_v=&  m_2 \times V_{cmax}, where \\
&m_2 = \frac{c_i}{c_i + k_c (1+O_i/k_o)},\\
&V_{cmax} = n_2 \times N_{leaf} \\
&n_2 =  f_{vcmax}(T) \times \beta_{soil}^{ps} \times \beta_{sinklim}^{ps}\times \beta_{soa}^{ps} \times v_{cmax}^n \times fN_{rub}
\end{align}$$

where $O_i$ is the intercellular partial pressure of , and $k_c$ and
$k_o$ are the Michaelis-Menten constants for and respectively, derived
$E^{kc}_0$, $E^{kc}_1$, or $E^{ko}_0$, $E^{ko}_0$ using
Eq. <a href="#eq:ps_t_sensitivity" data-reference-type="ref"
data-reference="eq:ps_t_sensitivity">[eq:ps_t_sensitivity]</a>.

The N-specific light-saturated rate of C3 photosynthesis can then be
calculated as:

$$m_{sat} = min(n_1 \times m_1, n_2 \times m_2)$$

(Friend et al. 2009) adjusted the scheme by (Kull and Kruijt 1998) using
the Collatz formulation of C4 photosynthesis (Collatz, Ribas-Carbo, and
Berry 1992). The simplified assumption is that $A_j$ and $A_v$ can be
calculated as above, but at saturating $c_i$ ($c_{i,max}$).
Bundle-sheath transport limitation ($A_p$) is then further limiting C4
photosynthesis, given by $$\label{eq:pepc}
\begin{align}
A_p=& V_{pmax} \times m_3, where \\
&V_{pmax} = n_3 \times N_{leaf} \\
&n_3 =  g_{pepc}(T) \times v_{pepc}^n \times fN_{pepc}, \\
&m_3 = \frac{c_i}{p},
\end{align}$$

where the temperature response is $$\label{eq:pepc_t_sensitivity}
g_{pepc}(T) = 2^{(T - T_{ref}^{pepc})/ T_{base}^{pepc}}$$

The N-specific light-saturated rate of C4 photosynthesis can be
calculated as
$$m_{sat} = min(n_1 \times m_1, n_2 \times m_2, n_3 \times m_3)$$

The light-harvesting limited rate of photosynthesis ($A_h$) can be
written as $$\label{eq:ps_lightharvesting}
\begin{align}
A_h =&  m_1 \times \alpha_i \times PPFD_a,\\
&PPFD_a = PPFD_0 (1-e^{-k_a \times C_{chl}}), where\\
&C_{chl} = a_{chl}^n \times fN_{chl} \times N_{leaf}
\end{align}$$

where $\alpha_i$ is the intrinsic quantum efficiency for uptake,
$PPFD_0$ is the photosynthetic photon flux density (PPFD) penetrating
sunlit or shaded foliage, converted from the adsorbed radiation of the
canopy layer (Eq.
<a href="#eq:rad:a_shadesunlit" data-reference-type="ref"
data-reference="eq:rad:a_shadesunlit">[eq:rad:a_shadesunlit]</a>), and
$k_a$ as well as $a_{chl}^n$ are parameters specified in Table
<a href="#tab:assimi:param" data-reference-type="ref"
data-reference="tab:assimi:param">[tab:assimi:param]</a>.

As (Kull and Kruijt 1998) show, this system of equations (Eq.
<a href="#eq:jmax" data-reference-type="ref"
data-reference="eq:jmax">[eq:jmax]</a>-<a href="#eq:ps_lightharvesting" data-reference-type="ref"
data-reference="eq:ps_lightharvesting">[eq:ps_lightharvesting]</a>) can
be solved to yield gross photosynthesis ($A _g$) for one canopy layer
and light-class as:

$$\begin{align}
A _g = (1-\frac{\Gamma^{*}}{c_i}) &[m_{sat}N_{sat}+\alpha_i PPFD_a (e^{-k_a n_3 N_{sat}}-e^{-k_a n_3 N_{leaf}}), where\\
N_{lim} &=  - \frac{ln(m_{sat}/[\alpha_i PPFD k_a n_3 m_1])}{k_a n_3}, and \\
&N_{lim} < 0: N_{sat}=0\\
&N_{lim} < N_{leaf}: N_{sat}=N_{lim}\\
&N_{lim} \ge N_{leaf}: N_{sat}=N_{leaf}
\end{align}$$

### Stomatal coupling

The combination of leaf-level net photosynthesis ($A_{n,cl}$), stomatal
conductance ($g_{s,cl}$), and leaf internal concentration ($c_{i,cl}$)
satisfying Eq. <a href="#eq:stomatal_coupling" data-reference-type="ref"
data-reference="eq:stomatal_coupling">[eq:stomatal_coupling]</a> is
sought iteratively for each canopy layer:

$$\label{eq:stomatal_coupling}
\begin{align}
A_n &= (A_{g,sunlit} \times f_{sunlit} + A_{g,shaded} \times (1- f_{sunlit})) - R_l\\
g_s &= [g_0 + g_1 \frac{A_n\beta_{air}\beta_{soil}^{gs}}{c_a}] \times R  \times T_{air} / p\\
c_i &= c_1 \times c_a - c_2 \times A_n \times ( \frac{D_{air}^{wv2co2}}{g_s} + \frac{D_{turb}^{wv2co2}}{g_a} ) \times R \times T_{air} 
\end{align}$$

where $f_{sunlit}$ is the fraction of sunlit leaves in a canopy layer
(see Eq. <a href="#eq:rad:fsunlit" data-reference-type="ref"
data-reference="eq:rad:fsunlit">[eq:rad:fsunlit]</a>), $R_l$ is the
maintenance respiration of leaves (see Sect.
<a href="#sec:veg:MaintResp" data-reference-type="ref"
data-reference="sec:veg:MaintResp">1.3.2</a>), $g_0$ and $g_1$ are
PFT-specific parameters, $R$ is the molar gas constant, $T_{air}$ is air
temperature (K), $p$ air pressure (Pa), $c_1$ converts concentration
from ppm to Pa, $c_2$ converts  mol m<sup>−2</sup> s<sup>−1</sup> to
mol m<sup>−2</sup> s<sup>−1</sup>, $g_a$ is the aerodynamic conductance,
and the $D$’s are the diffusion coefficient corrections for and water
(Bonan 2015).

$\beta_{air}$ is taken as relative humidity (Ball, Wooddrow, and Berry
1987; Knauer, Werner, and Zaehle 2015) and $\beta_{soil}^{gs}$ is the
stomatal response to soil moisture, described by:

$$\label{eq:veg:beta_soil_gs}
\beta_{soil}^{gs} = 1 - \frac{\Psi_{soil}}{\Psi_{leaf}^{min}}$$

where $\Psi_{soil}$ is the soil water potential in the root zone
(Eq. <a href="#eq:sse:soil_matrix_potential" data-reference-type="ref"
data-reference="eq:sse:soil_matrix_potential">[eq:sse:soil_matrix_potential]</a>)
and $\Psi_{leaf}^{min}$ is the PFT-specific minimum leaf water
potential.

### Canopy integration

Canopy-level fluxes are derived by summing the product of layer-level
fluxes or state-variables and the depth of the layer
$$F = \sum_{cl=1}^{ncanopy} F_{cl} \times LAI_{cl},$$

where $F$ is the canopy-level equivalent of a leaf-level variable
$F_{cl}$ (per unit leaf-area), such as gross photosynthesis ($A_g$), net
photosynthesis ($A_n$), and leaf-level stomatal conductance ($g_s$),
(see Sect. <a href="#sec:assimi:leaf_ps" data-reference-type="ref"
data-reference="sec:assimi:leaf_ps">1.2.2</a> and
<a href="#sec:assimi:stomatal_coupling" data-reference-type="ref"
data-reference="sec:assimi:stomatal_coupling">1.2.3</a>) and $LAI_{cl}$
is the leaf area index of the canopy layer. The resulting canopy net
assimilation $A_g$ is used as input to the vegetation model
(Eq. <a href="#eq:veg:labile_pool_dynamics" data-reference-type="ref"
data-reference="eq:veg:labile_pool_dynamics">[eq:veg:labile_pool_dynamics]</a>),
the canopy conductance ($G_s$) is used for the calculation of ecosystem
transpiration
(Eq. <a href="#eq:sse_surface_energy_balance" data-reference-type="ref"
data-reference="eq:sse_surface_energy_balance">[eq:sse_surface_energy_balance]</a>).

## Vegetation growth, turnover and dynamics

### Labile pool dynamics

The growth of a plant is dependent on the partitioning of its labile
resource to new tissue growth, storage production, as well as – in the
case of C – respiration for maintenance and resource uptake:

$$\label{eq:veg:labile_pool_dynamics}
\begin{align}
\frac{dC_{labile}}{dt} &= A_g + \Delta S_C - R_m - R_r - (1 + f_{resp,growth}) \times G_C \\
\frac{dN_{labile}}{dt} &= U_{root,N} + \Delta S_N  - G_N  \\
\frac{dP_{labile}}{dt} &= U_{root,P} + \Delta S_P - G_P 
\end{align}$$

where $R_m$ is maintenance respiration, $R_r$ is resource uptake
respiration (for both see Sect.
<a href="#sec:veg:MaintResp" data-reference-type="ref"
data-reference="sec:veg:MaintResp">1.3.2</a>), $\Delta S_X$ is the net
exchange between the labile and reserve pool (Sect.
<a href="#sec:veg:reserve" data-reference-type="ref"
data-reference="sec:veg:reserve">1.3.6</a>), $f_{resp,growth}$ is the
constant fraction of respiration associated with growth, $G_X$ are the
growth rates to build new tissues, $U_{root,X}$ are nutrient rates of
root uptake.

### Maintenance respiration

Maintenance respiration for every vegetation pool ($i$) is estimated
from its N content ($N_i$) as $$\label{eq:veg:respiration}
R_{m,i}=f_{temp} \times f_{resp,maint}^{i} \times N_i$$

where $f_{resp, maint}^{i}$ is the maintenance respiration per unit N,
which differs between woody and non-woody pools, and $f_{temp}$ is the
instantaneous temperature response of respiration ((Lloyd and Taylor
1994)) $$\label{eq:veg:t_resp}
f_{temp}=e^{t_{k1} (\frac{1}{t_{k2}} - \frac{1} {T - t_{k3}}})$$

where $t_{k1}, t_{k2}$ and $t_{k3}$ are temperature sensitivity
parameters and $T$ is the instantaneous air or soil temperature for
above- and belowground tissues, respectively. Following (Atkin, Meir,
and Turnbull 2014), the basal respiration rate acclimates to temperature
$$\label{eq:veg:resp_acclim}
f_{maint\_rate}=f_{maint\_rate\_ref}10^{f_{resp\_acclim}(T_{acclim}-T_{acclim,ref})}$$

where $f_{maint\_rate,ref}$ is the N-specific maintenance respiration
rate at the temperature $T_{acclim,ref}$, $f_{resp\_acclim}$ is the
slope of the temperature acclimation and $T_{acclim}$ is the running
average of air or soil temperature ($\tau_{mavg}^{resp}$), respectively.

Resource uptake respiration for nutrients is given by specific costs
($cost_i$, (Zerihun, McKenzie, and Morton 1998)) to transform nutrients
from mineral sources ($i$) into organic material and the actual plants
uptake $U_{plant,i}$
(<a href="#sec:sb:simple_uptake" data-reference-type="ref"
data-reference="sec:sb:simple_uptake">1.4.5</a>) as
$$\label{eq:veg:nutrient_costs}
R_{r,i}=cost_{i} \times U_{plant,i}$$

where $i$ is either or .

### Growth

The potential, source-limited growth rate ($G^{*}_X$) is given by the
product of the maximum turnover rate of the labile pool
($1/\tau_{labile}$) and the labile pool size ($X_{labile}$). Three sink
limitation processes operate, which control the down-regulation of this
potential growth rate to the actual growth rates ($G_X$):

- the temperature and moisture sensitivity of the meristem (Eq.
  <a href="#eq:veg:k0_lab" data-reference-type="ref"
  data-reference="eq:veg:k0_lab">[eq:veg:k0_lab]</a>);

- the requirement for maintenance respiration (Eq.
  <a href="#eq:veg:respiration" data-reference-type="ref"
  data-reference="eq:veg:respiration">[eq:veg:respiration]</a>), which
  takes priority over new growth; and

- the co-limitation of growth by the nutrients required to grow specific
  pools (Eq. <a href="#eq:veg:growth_rate" data-reference-type="ref"
  data-reference="eq:veg:growth_rate">[eq:veg:growth_rate]</a>).

The meteorological/environmental control on the meristem mechanism is
represented by a reduction of the maximum turnover rate at low
temperatures and low soil moisture conditions.

$$\label{eq:veg:k0_lab}
k^{*}_{labile} = \frac{1}{\tau_{labile}} \times e^{-(\lambda^{labile}_{temp}\times T_{air})^{k^{labile}_{temp}}} \times e^{-(\lambda^{labile}_{\Theta}\times \Theta)^{k^{labile}_{\Theta}}}$$

where $T_{air}$ is air temperature in degrees Celsius, $\Theta$ is the
fractional soil moisture content (Eq.
<a href="#eq:sse:theta" data-reference-type="ref"
data-reference="eq:sse:theta">[eq:sse:theta]</a>), and the $\lambda$ and
$k$ are parameters. $k^{*}_{labile}$ is set to zero outside the growing
season (see Sect. <a href="#sec:veg:phen" data-reference-type="ref"
data-reference="sec:veg:phen">1.3.8</a>).

Respiration is assumed to have priority over growth. However, under
severe C deficit, the meristem activity also down-regulates maintenance
and resource uptake respiration. To ensure that carbon growth ($G_C$)
does not exceed the size of the labile carbon pool, the turnover rate of
the labile pool to growth is corrected by the current respiration rate,
and constrained to positive solutions:

$$\label{eq:growth_c}
G^{*}_C = k^{*}_{labile} \times C_{labile} \times dt - R_m - R_r; G^{*}_C \ge 0$$

Given $G^{*}_C$, and the stoichiometric requirements for biomass growth
($req^{growth}_{NC}$, and $req^{growth}_{PN}$, respectively):

$$\label{eq:veg:req_growth}
req^{growth}_{XY} = \sum^{pools}_{i}{f^{i}_{alloc} \times \frac{X_i}{Y_i}}$$

where $f^{i}_{alloc}$ are the allocation fractions (Sect.
<a href="#sec:veg:allocation" data-reference-type="ref"
data-reference="sec:veg:allocation">1.3.4</a>) to each pool $i$, and
$\frac{X_i}{Y_i}$ are the target stoichiometries of C:N:P (Sect.
<a href="#sec:veg:stoichem" data-reference-type="ref"
data-reference="sec:veg:stoichem">1.3.5</a>) of the leaf, fine root,
coarse root, sapwood, and fruits pools. The actual growth rates can be
calculated as

$$\label{eq:veg:growth_rate}
\begin{align}
G_P &= req^{growth}_{N:P} \times G_N = req^{growth}_{P:N} \times req^{growth}_{N:C} \times G_C \label{eq:veg:growth_rate_x} \\
&G_C \le G^{*}_C \\
&G_N \le \frac{k^{nut}_{labile}}{\tau_{labile}} \times N_{labile} \times dt \\
&G_P \le \frac{(k^{nut}_{labile})^2}{\tau_{labile}} \times P_{labile} \times dt,
\end{align}$$

Note that only the minimum of the three rates in eq.
<a href="#eq:veg:growth_rate_x" data-reference-type="ref"
data-reference="eq:veg:growth_rate_x">[eq:veg:growth_rate_x]</a> can
actually be realised. The other two growth rates are adjusted, implying
a relative accumulation of these elements in the labile pool. The use of
$k^{nut}_{labile}$ implies a stronger mobilisation capacity for
nutrients than for carbon growth. The assumption behind this is that the
temperature and moisture control of the meristem is already accounted
for by $G_C$, and that the plant is able to mobilise the required
nutrients from the labile pool to support this growth.

Outside the growing season, all growth fluxes are set to zero (see Sect.
<a href="#sec:veg:phen" data-reference-type="ref"
data-reference="sec:veg:phen">1.3.8</a>).

### Growth partitioning

The labile pool partitioned to growth is first split into reproductive
(fruit pool) and structural (leaves, fine and coarse roots and sap wood)
growth. The latter is then further subdivided following a functional
balance relationship between the structural pools (Zaehle and Friend
2010).

The fraction of carbon growth allocated to fruit depends on the ability
of the plant to build reserves ($\Delta S_C$, see Sect.
<a href="#sec:veg:reserve" data-reference-type="ref"
data-reference="sec:veg:reserve">1.3.6</a>), implying that fruit growth
is suppressed during phases of rapid leaf growth (beginning of the
growing season), as well as periods of C starvation (e.g. severe
drought).

$$f^{fruit}_{alloc,C} = k1^{fruit}_{alloc} + (k2^{fruit}_{alloc} - k1^{fruit}_{alloc}) \times exp^{-(\lambda^{fruit}_{alloc}\times (\Delta S_C+k3^{fruit}_{alloc}))^{k4^{fruit}_{alloc}}} \label{eq:veg:falloc_fruit}$$

where the maximum fraction of allocation to fruits
($k2^{fruit}_{alloc}$) is a PFT-specific parameter. The allocation of
the remaining growth ($(1-f^{fruit}_{alloc,C}) \times G_C$) to the
structural pools leaves, fine root, coarse roots, and sapwood follows a
set of following allometric relationships (Zaehle and Friend 2010).

For both trees and grasses, fine root and leaves are assumed to be in
homeostatic balance described as $$\label{eq:veg:ltor}
C_{leaf} =  f^{ltor}_{N|P|H_20} \times k_{rtos} \times \frac{C_{leaf}}{C_{sap\_wood}} \times C_{fine\_root}$$

where $k_{rtos}$ is a PFT-specific parameter and $f^{ltor}$ is a
function describing the response of the leaf to root ratio to nutrient
and water limitation. It is calculated as the minimum of three functions
describing N, P and water limitation respectively, calculated as:

$$\label{eq:veg:dynamic_alloc}
f^{ltor} = min(\frac{N_{labile}/C_{labile}}{req^{growth}_{NC}},\frac{P_{labile}/N_{labile}}{req^{growth}_{PN}},\frac{W_{soil,root}}{W^{alloc}_{soil,crit}})$$

$X_{labile}$ refers to the content of the respective element in the
labile pool and $req^{growth}_{NC}$ and $req^{growth}_{PN}$ are the N:C
and P:N ratios required for growth, respectively (see eqn.
<a href="#eq:veg:req_growth" data-reference-type="ref"
data-reference="eq:veg:req_growth">[eq:veg:req_growth]</a>). For the
water limited allocation, $W^{alloc}_{soil,crit}$ is the critical level
of soil moisture below which root allocation increases and
$W_{soil,root}$ is calculated as the ratio between the current water
content in the root zone and the water content in the root zone at field
capacity.

In grasses, halms are assumed to be a proportion of leaf mass, and no
height restrictions apply. $$\label{eq:l2s_grass}
C_{leaf} = k_{htol} \times C_{sap\_wood}$$

In trees and shrubs, leaf and woody biomass are linked through the
pipe-model hypothesis (requiring a constant ratio of leaf area, $LA$, to
sapwood area, $SA$) $$\label{eq:veg:latosa}
C_{leaf} = \frac{k_{latosa} \times C_{sap\_wood}}{sla \times \rho_{wood} \times H} <=> LA = k_{latosa} \times SA,$$

where $sla$, $k_{latosa}$, and $\rho_{wood}$ are the PFT-specific
specific leaf area, leaf to sapwood area ratio, and wood density,
respectively. $H$ is the mean forest canopy height calculated as:
$$\label{eq:veg:allometry}
H = k1_{allom} \times D^{k2_{allom}},$$

where the diameter at breast height ($D$) is determined from woody
biomass, assuming that the entire trunk is a cylinder. As an extension
to the pipe-model theory below-ground, coarse root biomass is assumed to
be proportional to sap wood mass:

$$\label{eq:coarse_root}
C_{coarse\_root} = k_{ctos} \times C_{sap\_wood}$$

where $k_{ctos}$ is a PFT-specific parameters.

### Tissue stoichiometry

Following (Meyerholt and Zaehle 2015), C:N:P stoichiometry for
slow-overturning structural tissues (sap wood, coarse roots) as well as
fruits, is assumed to be time-invariant and modelled as dependent on the
PFT-specific mean foliar stoichiometry ($\chi^{C:N}_{leaf}$, and
$\chi^{N:P}_{leaf}$) and set ratios (see Table
<a href="#tab:veg:param" data-reference-type="ref"
data-reference="tab:veg:param">[tab:veg:param]</a>). Heartwood
stoichiometry differs from sapwood stoichiometry because a fraction of
the nutrients are retranslocated ($k^{wood}_{resorb}$) to the labile
pool upon heartwood formation.

Following (Zaehle and Friend 2010), the C:N and N:P ratios of leaves are
varied in response to the nutrient demand and supply so that:

$$\chi^{X:Y}_{leaf|t+1} = \chi^{X:Y}_{leaf} \times (1 + \delta^{\chi}_{leaf} \times \Gamma_{XY} )$$

where $\chi^{X:Y}_{leaf}$ denotes either the C:N or N:P ratio of the
leaves, $\delta^{\chi}_{leaf}$ is a parameter denoting the maximum
amount that leaf nutrients can change per timestep and $\Gamma_{X:Y}$
heuristically accounts for limits to the plasticity of foliar
stoichiometry as: $$\label{eq:veg:stoichem_dynamic}
    \Gamma_{X:Y} = 
    \left\{
    \begin{array}{ll}
        e^{-(\lambda^{\chi}_{leaf}\frac{ \chi^{X:Y}_{leaf}}{\chi^{X:Y}_{leaf,min}+\chi^{X:Y}_{leaf,max}})^{k^{\chi}_{leaf}}}  & 
        \mbox{if } \chi^{Y:X}_{labile} \leq req^{growth}_{Y:X} \\
        -(1 - e^{-(\lambda^{\chi}_{leaf}\frac{ \chi^{X:Y}_{leaf}}{\chi^{X:Y}_{leaf,min}+\chi^{X:Y}_{leaf,max}})^{k^{\chi}_{leaf}}}) & \mbox{if } \chi^{Y:X}_{labile} > req^{growth}_{Y:X}
    \end{array}
\right.$$ In the above, $\chi^{X:Y}_{leaf,min}$ and
$\chi^{X:Y}_{leaf,max}$ are PFT-specific parameters. The condition
refers to the relationship between the nutrients available for growth in
the labile pool and the nutrients required for growth (Eq.
<a href="#eq:veg:req_growth" data-reference-type="ref"
data-reference="eq:veg:req_growth">[eq:veg:req_growth]</a>), averaged at
the time-scale of $\tau_{mavg}^{\chi}$. The stoichiometric ratios of the
fine roots vary proportionally to those of the leaves, whereas the
stoichiometry of wood is assumed time-invariant (Meyerholt and Zaehle
2015). The stoichiometry of the labile and reserve pools are prognostic
properties, as described in
Sect. <a href="#sec:veg:reserve" data-reference-type="ref"
data-reference="sec:veg:reserve">1.3.6</a>.

### Long-term reserve dynamics

The target leaf carbon pool ($C^{target}_{leaf}$) is determined by the
current allometry, and corresponds to the leaf area index implied by
current sapwood area.

The target labile carbon pool size buffers short term fluctuations in
GPP, and is assumed to correspond to the maximum of the cumulated $GPP$
or $R_m$ over the turnover time of the labile pool, while the target of
the labile nitrogen and phosphorus pool corresponds to the average
stoichiometric growth requirement over the turnover time of the labile
pool: $$\label{eq:target_labile}
\begin{align}
C^{target}_{labile} &= max(\int_{t=-{\tau^{0}_{labile}}}^{t=0} GPP \times dt, \int_{t=-{\tau^{0}_{labile}}}^{t=0} R_m\times dt) \\
N^{target}_{labile} &= req^{growth}_{NC} \times C^{target}_{labile}  \\
P^{target}_{labile} &= req^{growth}_{PN} \times N^{target}_{labile} 
\end{align}$$

The target size of the reserve pool depends on the C required to replace
the annual growth of leaves and fine roots. It is scaled by a
PFT-specific constant as a measure of risk avoidance
($k^{target}_{reserve}$), with larger values indicating a preference for
storage over growth.

$$\label{eq:target_storage}
\begin{align}
C^{target}_{reserve} &= min(k_{reserve} \times (1+f_{resp,growth}) \times \frac{LAI^{target}}{sla}, \sum_{i=l,f,s} f_{reserve,max,i} \times C_{i}), where \\
&k_{store} = k^{target}_{reserve} \times ( min(1,\frac{1}{\tau_{leaf}}) + \frac{1}{\tau_{fine\_root} \times R_{leaf:fine\_root}} ) 
\end{align}$$

where $LAI^{target}$ is the target leaf area index, which is constrained
to values below $LAI^{target}_{max}$ (see Sect.
<a href="#sec:veg:allocation" data-reference-type="ref"
data-reference="sec:veg:allocation">1.3.4</a>) and $sla$ the
PFT-specific specific leaf area, $\tau_{leaf}$ and $\tau_{fine\_root}$
are the PFT-specific turnover times of foliage and fine roots, and
$R_{leaf:fine\_root}$ is the leaf to root ratio (Eq.
<a href="#eq:veg:dynamic_alloc" data-reference-type="ref"
data-reference="eq:veg:dynamic_alloc">[eq:veg:dynamic_alloc]</a>)
averaged over the lifetime of the fine roots ($\tau_{mavg}^{root}$). The
N and P target pools are defined in an equivalent manner, respecting the
current target stoichiometry of leaves and fine roots (Sect.
<a href="#sec:veg:stoichem" data-reference-type="ref"
data-reference="sec:veg:stoichem">1.3.5</a>).

The net exchange between the labile and reserve pool is calculated as

$$\label{eq:veg:phimaint}
\begin{align}
\Delta S_X &= \frac{1}{\tau_{labile}} \times (\Phi^{X}_{maint} \times X_{reserve} - \Phi^{X}_{store} \times X_{labile}) \times dt; with \\
&\Phi^{X}_{maint} = e^{-(\lambda^{\Phi}_{maint} \times \frac{X_{labile}}{X^{target}_{labile}})^{k^{\Phi}_{maint}}}, and \\
&\Phi^{X}_{store} = 1 - e^{-(\lambda^{\Phi}_{store} \times \frac{X_{reserve}}{X^{target}_{reserve}})^{k^{\Phi}_{store}}}\label{eq:veg:phistore}
\end{align}$$

where $\lambda^{\Phi}_{maint}$, $k^{\Phi}_{maint}$,
$\lambda^{\Phi}_{store}$, and $k^{\Phi}_{store}$ are parameters of a
Weibull-type function. Under conditions of severe resource stress (i.e.
low labile pool size corresponding to its target size), the build up of
reserves, is reduced according to

$$\label{eq:veg:phi_interact}
  \Phi^{X}_{store} = \frac{ 1 - \Phi^{X}_{maint}}{1 - k_{\Phi,inter}} \Phi^{X}_{store},\ if\ \Phi^{X}_{maint} > k_{\Phi,inter}$$

where $k_{\Phi,inter}$ is a parameter. To support leaf and fine-root
growth at the beginning of the growing season, $\Delta S_X$ is further
modified by the phenological pull ($\Phi^{X}_{phen}$) during the growing
season as follows:

$$\label{eq:veg:phiphen}
\begin{align}
\Delta S_X  &= \Delta S_X + k^{*}_{labile} \times \Phi^{X}_{phen} \times X_{reserve} \times dt, with \\
&\Phi^{X}_{phen} = e^{-(\lambda^{\Phi}_{phen} \times \frac{X_{leaf}}{X^{target}_{leaf}})^{k^{\Phi}_{phen}}} 
\end{align}$$

where $\lambda^{\Phi}_{phen}$, and $k^{\Phi}_{phen}$ are parameters of a
Weibull-type function.

### Photosynthetic sink limitation

In the case that the labile carbon pool exceeds its target size
substantially because growth is limited by temperature, moisture, or
because sufficient nutrients are lacking to allow growth (Eq.
<a href="#eq:veg:growth_rate_x" data-reference-type="ref"
data-reference="eq:veg:growth_rate_x">[eq:veg:growth_rate_x]</a>),
sink-limitation down-regulates photosynthetic activity so that:

$$\label{eq:sinklim_c}
\begin{align}
    \beta_{sinklim}^{ps} &= \beta_{sinklim,min}^{ps} + (1 - \beta_{sinklim,min}^{ps}) \times e^{-(\lambda_{sinklim}^{ps} \times X)^{k_{sinklim}^{ps}}}, where \\
    &X= \frac{C_{labile} - C_{labile}^{target}}{C_{labile}^{target}}
\end{align}$$

Here, $\beta_{sinklim, min}$, $\lambda_{sinklim}^{ps}$ and
$k_{sinklim}^{ps}$ are parameters (see Table
<a href="#tab:veg:param" data-reference-type="ref"
data-reference="tab:veg:param">[tab:veg:param]</a>) and
$C_{labile}^{target}$ is the target value for the labile pool (Eq.
<a href="#eq:target_labile" data-reference-type="ref"
data-reference="eq:target_labile">[eq:target_labile]</a>).

In addition, if the C:N or N:P ratios of the labile pool exceed those of
the target labile pool, the sink limitation factor is further modified
as a function of the stoichiometric ratio of the labile pool and that of
the labile target as:

$$\label{eq:sinklim_cnp}
\beta_{sinklim}^{ps} = \beta_{sinklim}^{ps} \times min(1,\frac{\chi^{N:C}_{labile}}{k^{CNP}_{sinklim} \times \chi^{N:C}_{labile,target}}, \frac{\chi^{P:N}_{labile}}{k^{CNP}_{sinklim} \times \chi^{P:N}_{labile,target}})$$

where $k^{CNP}_{sinklim}$ is a parameter.

### Phenology

The phenology of vegetation, describing the seasonal development of
foliage biomass, is simulated prognostically given the ability of the
plant to grow new tissues, which depends on the size and turnover of the
meristems (Eq. <a href="#eq:veg:growth_rate" data-reference-type="ref"
data-reference="eq:veg:growth_rate">[eq:veg:growth_rate]</a>), as well
as the fractional allocation of growth to plant organs (see Sect.
<a href="#sec:growth" data-reference-type="ref"
data-reference="sec:growth">1.3.3</a>). The start and end of the the
growing season are determined by meteorological triggers and soil
moisture, with plant growth set to zero outside the growing season (Eq.
<a href="#eq:veg:growth_rate" data-reference-type="ref"
data-reference="eq:veg:growth_rate">[eq:veg:growth_rate]</a>). The
meteorological variables determining these phenological triggers are
averaged over $\tau_{mavg}^{phen}$, to smooth out the effect of
day-to-day climate variability. While the beginning and ending mark the
start and end of tissue production, only the turnover of the leaves is
directly affected by phenological triggers. The turnover of all other
tissues is assumed to be constant (see Sect.
<a href="#veg:sec:turnover" data-reference-type="ref"
data-reference="veg:sec:turnover">1.3.9</a>).

The model differentiates evergreen, cold deciduous, rain deciduous tree
and shrub phenological strategies, as well as herbaceous perennial
phenological strategies.

The growing season start for cold deciduous and herbaceous PFTs is
described as a function of the accumulated growing degree days
($GDD_{acc}$) as:

$$\label{eq:phen:gdd}
\begin{align}
 GDD_{acc} &> GDD_{req}^{max} \times exp^{- k^{GDD}_{dormance}\times NDD}, where \\
 &\frac{GDD_{acc}}{dt} = GDD_{acc} + MAX(t_{air} - t_{tair}^{GDD},0.0)
\end{align}$$

where $GDD_{acc}$ denotes the current growing degree days above the
temperature threshold ($t_{air}^{GDD}$), $NDD$ is the number of dormancy
days, taken as days since the last growing season, and
$k^{GDD}_{dormance}$ is a PFT-specific parameter relating dormancy to
the PFT-specific maximum growing degree days requirement
($GDD_{req}^{max}$) to account for the chilling requirements of the buds
(Krinner et al. 2005).

For rain deciduous phenology, the start of the growing season is
triggered when the soil moisture stress factor ($\beta_{soil}^{gs}$, see
Eq. <a href="#eq:veg:beta_soil_gs" data-reference-type="ref"
data-reference="eq:veg:beta_soil_gs">[eq:veg:beta_soil_gs]</a>) is
larger than a PFT-specific threshold ($\beta_{soil}^{flush}$). This
criterion is also applied for herbaceous PFTs in addition to the
$GDD$-criterion.

For the evergreen phenology, recovery of photosynthesis in spring is
delayed according to the state of acclimation ($S$) to air temperature,
which reduces photosynthesis in spring until acclimation is reached
(Mäkelä et al. 2004). $S$ is calculated as

$$\label{eq:veg:soa}
\frac{dS}{dt} = \frac{1}{\tau_{soa}}(T_{air}-S)$$,

where $\tau_{soa}$ is a time constant. The reduction factor for Rubisco-
and electron transport limited photosynthesis ($\beta_{soa}$) is
calculated as

$$\label{eq:veg_beta_soa}
\beta_{soa} = \frac{(S - T_{min}^{soa})}{(T_{max}^{soa}-T_{min}^{soa})}$$

where $T_{min}^{soa}$ and $T_{max}^{soa}$ are parameters and
$\beta_{soa}$ is constrained to the range 0.1 and 1.

The end of the growing season for cold deciduous and herbaceous PFTs is
triggered by decreasing average air temperatures below a PFT-specific
temperature threshold ($t_{air}^{sen}$). For raingreen and herbaceous
PFTs the end of the growing season is triggered when the soil moisture
stress factor ($\beta_{soil}^{gs}$) becomes lower than a PFT-specific
threshold ($\beta_{soil}^{sen}$). In addition, herbaceous PFTs end their
growing season, once the weekly carbon balance ($GPP - R_m$) becomes
negative. Senescence is generally only introduced once the leaf age has
become larger than a PFT-specific threshold ($age^{leaf}_{min}$).

### Turnover

The breakdown of leaf and fine-root nutrients occurs at the time-scale
of $\tau_{{nut}\_{recycle}}$. The freed nutrients enter the labile pool
and are replaced by new nutrients of the labile pool according to the
current target C:N:P of the respective pool.

$$\label{eq:veg:nut_recycling}
\frac{X_{pool}}{dt} = (Y_{pool} \frac{X_{pool}^{target}}{Y_{pool}^{target}} - X_{pool}) \times 1/\tau_{{nut}\_{recycle}},$$

where $X$ is either N or P, and $Y$ C or N, respectively. The flux from
the labile pool is limited by the turnover rate and size of the labile
pool to ascertain that the latter cannot be exhausted.

The turnover time of most tissue types (fine and coarse roots, sapwood,
and fruits) is assumed constant for each PFT ($\tau_{fine\_root}$,
$\tau_{coarse\_root}$, $\tau_{sap\_wood}$, and $\tau_{fruit}$,
respectively). While roots and fruit turn directly into litter, only a
small fraction of sapwood ($f_{sap\_wood}^{branch}$) is turned to
litter, assuming it is lost as branches ($\tau_{branches}$), whereas the
predominant fraction of sapwood turns into non-respiring hardwood at the
timescale of $\tau_{sap\_wood}$. In evergreen trees, foliar turnover to
litter is assumed to be constant ($\tau_{leaves}$). For deciduous and
herbaceous PFTs only minor turnover happens at $\tau_{leaves}$ during
the growing season. At the end of the growing season (see Sect.
<a href="#sec:veg:phen" data-reference-type="ref"
data-reference="sec:veg:phen">1.3.8</a>), foliar turnover is set to a
constant rate $$\label{eq:veg:dec_leaf_turnover}
    f^{leaf}_{turn} = min(f_{shed,max} \times \frac{LAI^{target}}{LAI}, 1)$$

Resorption of nutrients to the labile pool during litterfall is assumed
to only occur during foliage turnover and heartwood formation at a
constant fraction ($k_{resorb}^{X}$), so that:

$$\label{eq:veg:resorb}
\begin{align}
&flux^{X}_{pool \rightarrow litter} = (1-k_{resorb}^{X}) \times \frac{  X_{pool}}{\tau_{pool}} \times dt, and \\
&flux^{X}_{pool \rightarrow labile} = k_{resorb}^{X} \times \frac{X_{pool}}{\tau_{pool}} \times dt,
\end{align}$$

where $flux^{X}_{pool \rightarrow litter}$ is the litterfall from any
one pool, and $flux^{X}_{pool \rightarrow labile}$ the retranslocated
flux into the labile pool.

### Vegetation dynamics

To assess stand density, we define for tree functional types the crown
area as:

$$\label{eq:veg:crownarea}
CA = k_{CA} \times D^{k_{rp}}$$

where $k_{CA}$ and $k_{rp}$ are parameters, $D$ is the tree diameter at
breast hight and $CA$ is constrained to be less then a maximum crown
area ($CA_{max}$). Using stand-scale LAI, individuum density and crown
area, the LAI of an individual tree is defined ($LAI_{ind}$), which is
used to calculated the foliage projective cover ($FPC$) as:

$$\label{eq:veg:fpc}
FPC = CA \times dens_{ind} ( 1 - e^{-k_{fpc} \times LAI_{ind}})$$

where $k_{fpc}$ is a paramter. To avoid strong seasonal cycles in
foliage projective cover for the calculation of vegetation dynamics
(Krinner et al. 2005; Zaehle and Friend 2010), $LAI_{ind}$ is diagnosed
from the sapwood area implied by the pipe-model (Eq.
<a href="#eq:veg:latosa" data-reference-type="ref"
data-reference="eq:veg:latosa">[eq:veg:latosa]</a>), implying that $FPC$
is essentially representing last year’s maximum LAI. For grasses, the
calculation of $FPC$ is not required and it is set to zero.

The establishment flux for a PFT is dependent on the size of the
seed-bed pool and the average, PFT-specific seed-bed turnover time
($\tau_{seed,est}$):

$$\label{eq:veg:flux_est}
flux_{est,X} = f_{dens} \times f_{temp} \times f_{moist} \times \frac{X_{seed\_bed}}{\tau_{seed,est}}$$

where $f_{temp}$ and $f_{moist}$ represent limitations for establishment
at low temperature and low moisture availability in the form of
Weibull-functions with parameters $\lambda_{est}^{env}$ and
$k_{est}^{env}$, where $env$ refers to either weekly air temperature
($T_{air}$) or weekly top-soil moisture ($\Theta_1$). Density dependency
of establishment ($f_{dens}$) is modelled as in (Sitch et al. 2003):

$$\label{eq:veg:est_dens}
f_{dens} = MAX(FPC_{max} - FPC, 0)$$

Three types of mortality are considered as additive processes,
growth-efficiency related mortality ($mort_{greff}$), density dependent
mortality ($mort_{dens}$), and a PFT-specific background mortality,
representing currently unaccounted for processes such as disturbance or
grazing.

$$\label{eq:mort_tot}
f_{mort} = MIN(mort_{greff} + mort_{dens} + mort_{bg,PFT}, 1)$$

Growth-efficiency mortality, represents any kind of mortality associated
with trees lacking the ability to defend themselves against stress (e.g.
pathogens) and is calculated as:

$$\label{eq:veg:mort_greff}
\begin{align}
mort_{greff} &= \frac{k1_{mort\_greff}}{1 + k2_{mort\_greff} \times eff_{growth}}, where \\
&eff_{growth} = \frac{NPP - \sum{Turnover_i}}{LAI} 
\end{align}$$

where, as in Sitch et al. (2003), $k1_{mort\_greff}$ and
$k2_{mort\_greff}$ are parameters and growth efficiency depends on net
primary production minus tissue turnover (of all tissues $i$) per unit
leaf area, calculated as running means over $\tau_{mavg}^{dynamics}$.

Space constraints in tree populations (as for grasses $FPC$ is zero) are
considered by constraining the foliage projected cover to a prescribed
maximum ($FPC_{max}$):

$$\label{eq:veg:mort_selfthining}
mort_{dens} = MAX(FPC - FPC_{max},0)$$

Litterfall from vegetation dynamics is then the product of the current
pool size and $f_{mort}$, scaled to the timestep of the model.

For trees, the appropriate number of individuals is also removed
following mortality. This does not affect the size of trees, as woody
biomass and stand density are modified proportionally. In total, the
change in vegetation individual density following establishment and
mortality is written as:

$$\label{eq:veg:dens_ind}
\frac{dens_{ind}}{dt} = f_{mort} \times dens_{ind} + \frac{flux_{est,C}}{k_{seed}}$$

where $flux_{est,C}$ is the carbon flux defined by Eq.
<a href="#eq:veg:flux_est" data-reference-type="ref"
data-reference="eq:veg:flux_est">[eq:veg:flux_est]</a>, and $k_{seed}$
is the PFT-specific seed size.

## Soil biogeochemistry

The dynamics of the soil organic pools ($X$) are described in general
as: $$\label{eq:sb:simple_om_dynamics}
\begin{align}
\frac{\partial}{\partial t} X_{met} &= \sum(f_{vp \rightarrow met} F_{L_{vp}}) + f_{wl \rightarrow met} \eta_{wl \rightarrow met} \frac{X_{wl}}{\tau_{wl}} - \frac{X_{met}}{\tau_{met}}\\
\frac{\partial}{\partial t} X_{str} &= \sum(f_{vp \rightarrow str} F_{L_{vp}}) + f_{wl \rightarrow str} \eta_{wl \rightarrow str} \frac{X_{wl}}{\tau_{wl}}- \frac{X_{str}}{\tau_{str}}\\
\frac{\partial}{\partial t} X_{wl} &= \sum(f_{vp \rightarrow wl} F_{L_{vp}}) - \frac{X_{wl}}{\tau_{wl}}\\
\frac{\partial}{\partial t} X_{fast} &= \eta_{l \rightarrow {fast}} (\frac{X_{met}}{\tau_{met}} + \frac{X_{str}}{\tau_{str}}) + \eta_{{slow} \rightarrow {fast}} \frac{X_{slow}}{\tau_{slow}} - \frac{X_{fast}}{\tau_{fast}} + \Phi_{l \rightarrow {fast}} + \Phi_{{slow} \rightarrow {fast}} + \frac{\partial}{\partial z} (D_b \frac{\partial X_{fast}}{\partial z}) \\
\frac{\partial}{\partial t} X_{slow} &= \eta_{{fast} \rightarrow {slow}} \frac{X_{fast}}{\tau_{fast}} - \frac{X_{slow}}{\tau_{slow}}  + \Phi_{{fast} \rightarrow {slow}} + \frac{\partial}{\partial z} (D_b \frac{\partial X_{slow}}{\partial z})
\end{align}$$

where $F_{L_{vp}}$ is the litterfall of the various plant tissue types,
$f_{vp \rightarrow i}$ are the coefficients determining the partitioning
of this litterfall to the litter pools (see Section
<a href="#sec:sb:litter_partitioning" data-reference-type="ref"
data-reference="sec:sb:litter_partitioning">1.4.1</a>), $\tau_i$ are
temperature and moisture adjusted, nitrogen-limited turnover times of
the respective pools (X; i = met, str, wl, fast, slow; see Section
<a href="#sec:sb:formation" data-reference-type="ref"
data-reference="sec:sb:formation">1.4.3</a>), $\eta_{i \rightarrow j }$
are the mass transfer from pool $i$ to $j$ (see Section
<a href="#sec:sb:formation" data-reference-type="ref"
data-reference="sec:sb:formation">1.4.3</a>), the $\Phi_x$ are the net
mineralisation terms for N and P, respectively, required to balance the
carbon inflow to the fast and slow SOM pools and their respective C:N:P
stoichiometry (see Sect.
<a href="#sec:sb:formation" data-reference-type="ref"
data-reference="sec:sb:formation">1.4.3</a>). The transfer of soil
organic matter through bioturbation is represented with a prescribed
diffusion constant (see Sect.
<a href="#sec:sb:bioturbation" data-reference-type="ref"
data-reference="sec:sb:bioturbation">1.4.4</a>).

The dynamics of the inorganic nitrogen pools are described as:
$$\label{eq:sb:simple_im_dynamics}
\begin{align}
\frac{\partial}{\partial t} NH_4 &= F_{dep,NH_4} - U_{plant,NH_4} - \sum(\Phi_{i,NH_4}) - U_{nit} - \frac{\partial v_{NH_4} NH_4}{\partial z}\\
\frac{\partial}{\partial t} NO_3 &= F_{dep,NO_3} + F_{nit,NO3} - U_{plant,NO_3} - \sum(\Phi_{i,NO_3}) - U_{denit} - \frac{\partial v_{NO_3} NO_3}{\partial z}\\
\frac{\partial}{\partial t} NO_y &= F_{nit,NO_y} + F_{denit,NO_y} - E_{NO_y} \\
\frac{\partial}{\partial t} N_2O &= F_{nit,N_2O} + F_{denit,N_2O} - E_{N_2O} \\
\frac{\partial}{\partial t} N_2 &=  F_{nit,N_2} + F_{denit,N_2} - E_{N_2}
\end{align}$$

where $U$ are the uptake rates of plants, or (de-)nitrifying bacteria,
respectively (see Section
<a href="#sec:sb:simple_uptake" data-reference-type="ref"
data-reference="sec:sb:simple_uptake">1.4.5</a> and
<a href="#sec:sb:nitrification_denitrification"
data-reference-type="ref"
data-reference="sec:sb:nitrification_denitrification">1.4.7</a>,
respectively); the $F_{dep}$ are the atmospheric deposition fluxes; the
$F_{nit,i}$, and $F_{denit,i}$ are the production of , and by
nitrification and denitrification, respectively; , and
$\frac{\partial v_{x} X}{\partial z}$ the vertical transport loss term
given by the product of ion concentration and water mass flow between
soil layers (see Sect.
<a href="#sec:sse:surface_hydrology" data-reference-type="ref"
data-reference="sec:sse:surface_hydrology">1.6.3</a>). Sorption of is
not explicitly modelled, and is accounted for by a reduced mobility in
water ($f_{leach,{NH_4}}$).

The dynamics of the inorganic phosphorus pools are described as:
$$\label{eq:sb:pi_dynamics}
\begin{align}
\frac{\partial}{\partial t} PO_4 &= F_{dep,PO_4} + F_{weath,PO_4} + F_{biomin,PO_4} - U_{plant,PO_4} - F_{adsorp,PO_4} - \sum(\Phi_{i,PO_4}) -  \frac{\partial v_{PO_4} PO_4}{\partial z}\\
\frac{\partial}{\partial t} P_{lab} &= F_{adsorp,PO_4} - F_{desorp,PO_4} + \frac{\partial}{\partial z} (D_b \frac{\partial P_{lab}}{\partial z})\\
\frac{\partial}{\partial t} P_{sorb} &= F_{desorp,PO_4} - F_{occlusion,PO_4} + \frac{\partial}{\partial z} (D_b \frac{\partial P_{sorb}}{\partial z})\\
\frac{\partial}{\partial t} P_{ocl} &= k_{ocl} P_{sorb} + \frac{\partial}{\partial z} (D_b \frac{\partial P_{ocl}}{\partial z})\\
\frac{\partial}{\partial t} P_{primary} &= - F_{weath,PO_4}
\end{align}$$

where $P_{lab}$, $P_{sorb}$, $P_{ocl}$, and $P_{primary}$ are labile,
absorbed, occluded, and primary P, respectively; the $F_{dep,PO_4}$,
$F_{weath,PO_4}$, $F_{biomin,PO_4}$, $F_{adsorp,PO_4}$, and
$F_{adsorp,PO_4}$ are the atmospheric deposition, weathering, fast
adsorption, and phosphorus fluxes, respectively (see Section
<a href="#sec:sb:p_flux" data-reference-type="ref"
data-reference="sec:sb:p_flux">1.4.8</a>). All pools except the primary
phosphorus pool are assumed to be affected by bioturbation (see
Sect. <a href="#sec:sb:bioturbation" data-reference-type="ref"
data-reference="sec:sb:bioturbation">1.4.4</a>).

### Partitioning of litterfall to litter pools

Non-woody litterfall is partitioned to the metabolic and structural
litter according to the CENTURY approach (Parton et al. 1993). Litter
from labile and reserve pools is assumed to enter the metabolic pools,
litter from sap- and heartwood enters the woody pool. The metabolic
fraction of litterfall from each vegetation pool ($_{vp}$, i.e. leaves,
fine and coarse roots, fruits and seed-bed) is determined as:
$$\label{eq:sb:simple_litterpartitioning_c}
f_{vp \rightarrow met,C} = f_{met,max,C} - k_{met,C} \times LC_{vp} \frac{C_{vp}}{N_{vp}}$$

where $f_{vp \rightarrow met,C}$ is constrained to positive solutions,
$f_{met,max,C}$ is the maximum fraction allocated to the metabolic pool,
$k_{met,C}$ a factor relating the metabolic litter fraction to the
lignin to nitrogen ratio, $LC_{vp}$ the tissue-specific fraction of the
lignin content of that tissue type, and $\frac{C_{vp}}{N_{vp}}$ the C:N
ratio of litterfall from that tissue. The lignin content is assumed
constant for all but the leaf tissues. For the latter, an empirical
dependency between lignin content and specific leaf-area ($sla$) is used
(White et al. 2000). $$\label{eq:sb:simple_leaf_lc}
LC_{leaf} = LC_{leaf,max} + k_{leaf2sla} \times sla$$

The remainder of litterfall is allocated to the structural pool. For N
and P, the partitioning assumes that the relative proportions of C:N and
N:P are preserved in the partitioning according to:
$$\label{eq:sb:simple_litterpartitioning_np}
f_{vp \rightarrow met,X} = \frac{1}{1+\frac{1-f_{vp \rightarrow met,C}}{k_{met,vp,X} \times f_{vp \rightarrow met,C}}}$$

Woody decomposition is assumed to be a two-stage process. The first step
implies physical destabilisation and a first level of biochemical
processing, which releases a constant fraction of carbon (1 -
$\eta_{C,wl \rightarrow met,str}$) to heterotrophic respiration. During
this step, a fraction of the nutrients (1 - $\eta_{N}|\eta_{P}$) is
leached to the mineral phase to account for inefficiencies of the
microbiota in mineral processing decomposing wood. The remainder is
assumed to enter the metabolic and structural litter
(Eq. <a href="#eq:sb:simple_litterpartitioning_c" data-reference-type="ref"
data-reference="eq:sb:simple_litterpartitioning_c">[eq:sb:simple_litterpartitioning_c]</a>
and
<a href="#eq:sb:simple_litterpartitioning_np" data-reference-type="ref"
data-reference="eq:sb:simple_litterpartitioning_np">[eq:sb:simple_litterpartitioning_np]</a>).

### SOM and litter turnover rates

The turnover times ($\tau^{base}_i$) of the litter and SOM pools respond
to soil temperature ($T_{soil}$) and soil matrix potential
($\Psi_{soil}$) as follows: $$\label{eq:sb:simple_som_turnover}
\begin{align}
\tau^{*}_i &= \tau^{base}_i \times f(T_{soil}) \times g(\Psi_{soil}),\ where\\
f(T_{soil}) &= \frac{E_{d,decomp} \times e^{E_{a,decomp} \times T}}{E_{d,decomp} - E_{a,decomp} \times ( 1 - e^{E_{d,decomp} \times T})},\ with\\
T &= \frac{T_{soil} - T_{opt,decomp}}{T_{soil} \times T_{opt,decomp} \times R},\ and \\ 
g(\Psi_{soil}) &= 1 - \Psi_{soil} / \Psi_{dec,min}
\end{align}$$

### SOM formation

Matter entering the fast and slow SOM pool (Eq.
<a href="#eq:sb:simple_om_dynamics" data-reference-type="ref"
data-reference="eq:sb:simple_om_dynamics">[eq:sb:simple_om_dynamics]</a>)
is required to fulfill the prescribed stoichiometry of the SOM pools
($\chi_{SOM}$). These are assumed constant with the exception of the
fast SOM C:N ratio, which varies with available following CENTURY
(Parton et al. 1993): $$\label{eq:sb:SOMCNfast}
\chi_{SOM_{fast}^{C:N}} = MAX(\chi_{SOM_{fast,max}^{C:N}} - f_{\chi} \times NH_4,\chi_{SOM_{fast,min}^{C:N}})$$

where $\chi_{SOM_{fast,min}^{C:N}}$, $\chi_{SOM_{fast,max}^{C:N}}$ and
$f_{\chi}$ are parameters.

The difference in stoichiometry of the matter entering the pool and the
required stoichiometry of the pool leads to the estimate of the
potential immobilisation flux:
$$\label{eq:sb:simple_pot_n_immobilisation}
\Phi^{*}_{l \rightarrow {fast}, NH_4} = \frac{\eta_{C,litter \rightarrow fast}}{\chi_{SOM_{fast}^{C:N}}} (\frac{C_{met}}{\tau^{*}_{met}} + \frac{C_{str}}{\tau^{*}_{str}}) - \eta_N (\frac{N_{met}}{\tau^{*}_{met}} + \frac{N_{str}}{\tau^{*}_{str}})$$

where $\tau^{*}_i$ are the temperature and moisture constrained turnover
times (Eq.
<a href="#eq:sb:simple_som_turnover" data-reference-type="ref"
data-reference="eq:sb:simple_som_turnover">[eq:sb:simple_som_turnover]</a>).

The actual immobilisation rate ($\Phi_{l \rightarrow {fast}, NH_4}$) is
limited to the amount of ammonium ( available, subject to co-occurring
potential N uptake from plants ($U^{*}_{NH_4,plant}$, Eq.
<a href="#eq:sb:simple_potential_uptake" data-reference-type="ref"
data-reference="eq:sb:simple_potential_uptake">[eq:sb:simple_potential_uptake]</a>),
and nitrifiers ($U^{*}_{nit}$; Eq.
<a href="#eq:sb:simple_dynamic_nitrification" data-reference-type="ref"
data-reference="eq:sb:simple_dynamic_nitrification">[eq:sb:simple_dynamic_nitrification]</a>).
Note that, similar as for the plant uptake, the uptake of ammonium is
limited with a Michaelis-Menten kinetic to account for reduced
accessibility of N at very low values.

$$\label{eq:sb:simple_act_n_immobilisation}
\Phi_{l \rightarrow {fast}, NH_4} = \frac{NH_4}{max(NH_4, \ U^{*}_{NH_4,plant} + U^{*}_{nit} + \Phi^{*}_{l \rightarrow {fast}, NH_4})} \times {\Phi^{*}_{l \rightarrow {fast}, NH_4}}$$

In the case that $\Phi_{l \rightarrow {fast}, NH_4}$ is insufficient to
ensure that the newly formed fast SOM has a C:N ratio of
$\chi_{SOM_{fast}^{C:N}}$, the turnover times of the metabolic and
structural litter pool is increased to match the N available for
immobilisation ((Parton et al. 1993)):

$$\tau_{met|str} = \tau^{*}_{met|str} \frac{\frac{\eta_{l \rightarrow {fast},C}}{\chi_{SOM_{fast}^{C:N}}} (\frac{C_{met}}{\tau^{*}_{met}} + \frac{C_{str}}{\tau^{*}_{str}}) - \eta_N (\frac{N_{met}}{\tau^{*}_{met}} + \frac{N_{str}}{\tau^{*}_{str}})}{\Phi_{l \rightarrow {fast},NH_4}}$$

Should the available be insufficient to maintain the uptake rates of
plants
(Eq. <a href="#eq:sb:simple_potential_uptake" data-reference-type="ref"
data-reference="eq:sb:simple_potential_uptake">[eq:sb:simple_potential_uptake]</a>)
and nitrifiers
(Eq. <a href="#eq:sb:simple_dynamic_nitrification" data-reference-type="ref"
data-reference="eq:sb:simple_dynamic_nitrification">[eq:sb:simple_dynamic_nitrification]</a>),
these fluxes are downregulated in proportion.

The potential immobilisation flux of phosphorus
($\Phi^{*}_{P,l \rightarrow {fast}}$) is defined in a similar manner as
potential $NH_4$ immobilisation, but now considering the actual turnover
time of the litter pools: $$\label{eq:sb:simple_pot_p_immobilisation}
\begin{align}
\Phi^{*}_{P,l \rightarrow {fast}} &= \frac{\eta_{C,litter \rightarrow {fast}}}{\chi_{{SOM}_{fast}^{C:N}}\chi_{{SOM}_{fast}^{N:P}}} (\frac{C_{met}}{\tau_{met}} + \frac{C_{str}}{\tau_{str}}) -
\eta_P (\frac{P_{met}}{\tau_{met}} + \frac{P_{str}}{\tau_{str}})\\
\Phi_{P,l \rightarrow {fast}} &= \frac{PO_4}{max(PO_4, \ U^{*}_{PO_4,plant} + \Phi^{*}_{P,l \rightarrow {fast}})} \times {\Phi^{*}_{P,l \rightarrow {fast}}}
\end{align}$$

Because the C:N:P stoichiometry and uptake use-efficiencies are
organised such that decomposition of these pools is always leading to
net mineralisation of nutrients, the $\Phi_{{slow} \rightarrow {fast}}$
values are negative and do not require special treatment to affect the
carbon-use efficiency or turnover rates (i.e.
$\tau_{fast}=\tau_{fast}^{*}$, and $\tau_{slow}=\tau_{slow}^{*}$). The
processing of fast and slow SOM is assumed to also include higher-order
trophic levels of heterotrophic respiration (Parton et al. 1993),
therefore only a fraction of the respired material
($\eta_{C, {fast} \rightarrow {slow}}$ and
$\eta_{C, {slow} \rightarrow {fast}}$) is assumed to enter the
subsequent pool:

$$\label{eq:sb:mineralization}
\begin{align}
\Phi_{{fast} \rightarrow {slow},NH_4} &= \frac{\frac{\eta_{C, {fast} \rightarrow {slow}}}{\chi_{{SOM}_{slow}^{C:N}}} C_{fast} - N_{fast}}{\tau_{fast}} \\
\Phi_{{slow} \rightarrow {fast},NH_4} &= \frac{\frac{\eta_{C, {slow} \rightarrow {fast}}}{\chi_{{SOM}_{fast}^{C:N}}} C_{slow} - N_{slow}}{\tau_{slow}} \\
\Phi_{{fast} \rightarrow {slow},PO_4} &= \frac{\frac{\eta_{C, {fast} \rightarrow {slow}}}{\chi_{{SOM}_{slow}^{C:N}}\chi_{{SOM}_{slow}^{N:P}}} C_{fast} - P_{fast}}{\tau_{fast}} \\
\Phi_{{slow}\rightarrow {fast},PO_4} &= \frac{\frac{\eta_{C, slow \rightarrow fast}}{\chi_{{SOM}_{fast}^{C:N}}\chi_{{SOM}_{fast}^{N:P}}} C_{slow} - P_{slow}}{\tau_{slow}} - F_{biomin,PO_4}
\end{align}$$

### Bioturbation

Bioturbation is treated as simple diffusive flux with a rate constant
$D_b$, as in (Ahrens et al. 2015), but declining with soil depth in
proportion to the fraction of roots in the layer to account for reduced
biological activity with increasing soil depth:

$$\label{eq:sb:bioturbation}
\begin{align}
D_b &= \frac{root_{frac}}{dz} \times \frac{k_{org}^{diff}}{\rho_{soil}^{cor}}, \ and\\
\rho_{soil}^{cor} &= MAX(\rho_{org}^{bulk}, \rho_{OM} + \rho_{soil} - \rho_{OM}\frac{\rho_{soil}}{\rho_{org}^{bulk}})
\end{align}$$

where $root_{frac}$ and $dz$ are the root fraction and depth of the soil
layer, $\rho_{soil}^{cor}$ is the soil bulk density corrected with soil
organic matter, $\rho_{org}^{bulk}$ is the bulk density of organic
material, $\rho_{OM}$ is the organic matter density of the soil layer
which depends on the organic matter content in the soil layer,
$\rho_{soil}$ is the bulk density of fine mineral soil, and
$k_{org}^{diff}$ is the diffusion coefficient for organic material due
to bioturbation.

### Plant uptake rates

The potential uptake rates of plants for X= , , and follow an extended
Michaelis-Menten kinetics:

$$\label{eq:sb:simple_potential_uptake}
\begin{align}
U^{*}_X &= v_{max,X}(T_{soil},\Psi) \times X \times (K_{m1,X}(T_{soil},\Theta)+ \frac{1}{K_{m2,X}(T_{soil},\Theta) + X}) \times f^X_{demand} \times C_{fine\_root},\ where\\
&v_{max,X}(T_{soil},\Psi) = v_{max,X} \frac{E_{d,uptake} \times e^{E_{a,uptake} \times T}}{E_{d,uptake} - E_{a,uptake} \times ( 1 - e^{E_{d,uptake} \times T})} \times \frac{\Psi_{fine\_root}}{\Psi_{leaf,min}},\ with\\
&T = \frac{T_{soil} - T_{opt,uptake}}{T_{soil} \times T_{opt,uptake} \times R},\ and \\ 
&K_{m1,X}(T_{soil},\Theta) = K_{m1,X} / (e^{-\frac{Ea,hsc}{R} \times (\frac{1}{T_{soil}} - \frac{1}{T^{\tau}_{ref}})} \times (\frac{\Theta}{\Theta_{fc}})^{k_{hsc}}),\ and\\
&K_{m2,X}(T_{soil},\Theta) = K_{m2,X} \times e^{-\frac{Ea,hsc}{R} \times (\frac{1}{T_{soil}} - \frac{1}{T^{\tau}_{ref}})} \times (\frac{\Theta}{\Theta_{fc}})^{k_{hsc}},\ and\\
&f^X_{demand} = 1 - e^{-(\frac{\chi^{X:Y}_{max}-\chi^{X:Y}_{labile}}{\chi^{X:Y}_{max} \times (1 - K_{demand}^{half,X})})^{k_{demand}}}  
\end{align}$$

where $v_{max,X}$ is the PFT-specific temperature-sensitive maximum
uptake rate per unit biomass, adjusted by the current root zone moisture
potential ($\Psi_{fine_root}$) to account for limited transport of
nutrients towards the roots in dry soils, $C_{fine\_root}$ is the
biomass density of fine roots (mol   C   m<sup>−3</sup>, see Eq.
<a href="#eq:veg:root_distribution" data-reference-type="ref"
data-reference="eq:veg:root_distribution">[eq:veg:root_distribution]</a>),
$T_{soil}$ is the soil temperature and the $K_m$ parameters are nutrient
sensitivities of the low and high affinity transporters. These
affinities are assumed to be temperature sensitive and are adjusted to
soil moisture to account for the difference between mass-based and soil
solution concentrations (Ahrens et al. 2015). The potential uptake of
nutrients can be down-regulated by plants given their internal demand
$f^X_{demand}$, where $X$ refers to either N or P and $X:Y$ refers to
either N:C or P:N ratios. $\chi^{X:Y}_{max}$ corresponds to the X:Y
ratio of growing a unit of leaves and fine roots at the current
leaf-to-root ratio (see
<a href="#sec:veg:allocation" data-reference-type="ref"
data-reference="sec:veg:allocation">1.3.4</a>, $K_{demand}^{half, X}$ is
a parameter denoting the fraction of $\chi^{X:Y}_{max}$ at which uptake
is reduced to 50% and $k_{demand}$ is shape parameter.

### Asymbiotic biological nitrogen fixation

The asymbiotic biological nitrogen fixation (BNF) is represented as:

$$\label{eq:sb:asymb_bnf}
F_{BNF}^{NH_4} = v_{max,BNF} \times f(T_{soil})$$

where $v_{max,BNF}$ is a parameter representing the base rate of
fixation and the temperature response is calculated as above (Eq.
<a href="#eq:sb:simple_som_turnover" data-reference-type="ref"
data-reference="eq:sb:simple_som_turnover">[eq:sb:simple_som_turnover]</a>).
$BNF$ is suppressed if the sum of and exceeds a critical threshold
$N_{limit}^{BNF}$ (Zaehle, Friedlingstein, and Friend 2010). The
distribution of $F_{BNF}^{NH_4}$ across soil layers follows the
distribution of fine roots, as indicator for C inputs into the soil. All
N fixed through this mechanism is added to the mineral soil pool.

### Nitrification and denitrification

Calculation of nitrification and denitrification follows (Xu-Ri and
Prentice 2008; Zaehle et al. 2011), which relies on the separation of
the soil into aerobic and anaerobic volume fractions ($anvf$):
$$\label{eq:sb:simple_anaerobic_vf}
\begin{align}
anvf &= e^{-(\lambda_{anvf} \times (1-afps))^{k_{anvf}}}, where\\
&afps = \frac{W_{fc}-W_{soil}}{W_{fc}} 
\end{align}$$

where $\lambda_{anvf}$ and $k_{anvf}$ are parameters, $afps$ is the air
filled pore space, and the $W_x$ are the soil moisture contents as
defined in Sect.
<a href="#sec:sse:surface_hydrology" data-reference-type="ref"
data-reference="sec:sse:surface_hydrology">1.6.3</a>.

The potential rate of nitrification ($U^{*}_{nit}$) in the aerobic
fraction of the soil is modified by temperature and soil moisture
according to:

$$\label{eq:sb:simple_dynamic_nitrification}
\begin{align}
U^{*}_{nit} &= v_{max,nit} \times f(T_{soil}) \times g(\Theta) \times NH_4, where\\
&f(T_{soil}) = E_{d,nit} \frac{e^{\frac{E_{a,nit}\times kt}{R_{gas}}}}{E_{d,nit}-E_{a,nit}\times (1 - e^{\frac{E_{d,nit} \times kt}{R_{gas}}})} \\
&kt = \frac{T_{soil}-T_{opt,nit}}{T_{soil} \times T_{opt,nit}}\\
&g(\Theta) = 1-afps
\end{align}$$

The actual rate of nitrification ($U_{nit}$), given the potential rate
and competing demands from plant and microbial uptake
(Sect. <a href="#sec:sb:formation" data-reference-type="ref"
data-reference="sec:sb:formation">1.4.3</a>, is partitioned into its
products (, , and ) according to
$$\label{eq:sb:simple_dynamic_nitrification_prod}
\begin{align}
F_{nit,NO_3} &= (1 - f_{nit}^{NO_y} - f_{nit}^{N_2O}) \times U_{nit}\\
F_{nit,NO_y} &= f_{nit}^{NO_y}  \times U_{nit}\\
F_{nit,N_2O} &= f_{nit}^{N_2O}  \times U_{nit}
\end{align}$$

where the $F_{nit,X}$ are the nitrification fluxes for , , and ,
respectively
(Eq. <a href="#eq:sb:simple_im_dynamics" data-reference-type="ref"
data-reference="eq:sb:simple_im_dynamics">[eq:sb:simple_im_dynamics]</a>).

The potential rate of denitrification ($U^{*}_{denit}$) in the anaerobic
fraction of the soil is modified by temperature:

$$\label{eq:sb:simple_dynamic_denitrification}
\begin{align}
U^{*}_{denit} &= anvf \times v_{max,denit}(T_{soil}) \times \frac{C_{fast}}{K_{m,denit}^{fast} + C_{fast}} \frac{NO_3}{K_{m,denit}^{NO_3} + NO_3}, where\\
&v_{max,denit}(T_{soil}) = v_{max,denit} \times e^{-\frac{E_{a,denit}}{R} \times (\frac{1}{T_{soil}} - \frac{1}{T^{\tau}_{ref}})}
\end{align}$$

The actual rate of denitrification ($U_{denit}$), given the potential
rate and competing demands from plant uptake
(Sect. <a href="#sec:sb:formation" data-reference-type="ref"
data-reference="sec:sb:formation">1.4.3</a>, is partitioned into its
products (, , and ) according to
$$\label{eq:sb:simple_dynamic_denitrification_prod}
\begin{align}
F_{denit,NO_y} &= f_{denit}^{NO_y}  \times U_{denit}\\
F_{denit,N_2O} &= f_{denit}^{N_2O}  \times U_{denit}\\
F_{denit,N_2} &= (1 - f_{denit}^{NO_y} - f_{denit}^{NO_y}) \times U_{denit}
\end{align}$$

where the $F_{denit,X}$ are the denitrification fluxes for , , and ,
respectively
(Eq. <a href="#eq:sb:simple_im_dynamics" data-reference-type="ref"
data-reference="eq:sb:simple_im_dynamics">[eq:sb:simple_im_dynamics]</a>).
The model currently ignores the effect of ammonia volatilisation, which
is of low relevance for natural, unfertilised ecosystems.

### Phosphorus weathering and biomineralisation

Weathering is modelled following (Wang, Law, and Pak 2010) as:

$$\label{eq:sb:weathering}
\begin{align}
F_{weath,PO_4} &=f(T_{soil}) \times g(\Theta) \times f(C_{fine\_root}) \times k_{weath} \times \rho_{soil}^{cor},\ where\\
&f(T_{soil}) = e^{-\frac{Ea,hsc}{R} \times (\frac{1}{T_{soil}} - \frac{1}{T^{\tau}_{ref}})},\\
&g(\Theta) =(\frac{\Theta}{\Theta_{fc}})^3, and \\
&f(C_{fine\_root}) =\frac{C_{fine\_root}}{K_{m,weath}^{root}+C_{fine\_root}}
\end{align}$$ where $k_{weath}$ is the rate constant for weathering, and
$\rho_{soil}^{cor}$ is the soil bulk density corrected by SOM content.
The weathering rate decreases with soil depth as the fine root C
decreases and is modified by soil temperature and moisture.

The potential biomineralisation rate of is determined as an additional
turnover of the P contained in the slow SOM pool, modified by
temperature and moisture modifiers, and affected by the concentration of
$PO_4$ and the root biomass:

$$\label{eq:sb:p_biomin}
\begin{align}
F^{*}_{biomin,PO_4} &= \frac{C_{slow}}{\chi_{{SOM}_{slow}^{C:N}} \times \chi_{{SOM}_{slow}^{N:P}} \times \tau_{biomin}} \times f(C_{fine\_root})  \times f(PO_4) \times  f(T_{soil}) \times  g(\Theta),\ where\\
&f(C_{fine\_root}) = \frac{C_{fine\_root}}{K_{m,biomin}^{root} + C_{fine\_root}},\ and\\
&f(PO_4) = \frac{K_{m,biomin}^{PO_4}}{K_{m,biomin}^{PO_4}+PO_4}
\end{align}$$ where $K_{m,biomin}^{root}$ and $K_{m,biomin}^{PO_4}$ are
constants constraining the biomineralisation rate under low root biomass
and high $PO_4$ concentration, respectively; the temperature and
moisture responses are calculated as those in Eq.
<a href="#eq:sb:simple_som_turnover" data-reference-type="ref"
data-reference="eq:sb:simple_som_turnover">[eq:sb:simple_som_turnover]</a>.
The biomineralisation rate is further constrained so that it does not
alter the stoichiometry of the fast pool.
$$F_{biomin,PO_4} = MIN(F^{*}_{biomin,PO_4}, \frac{\frac{\eta_{C, {slow} \rightarrow {fast}}}{\chi_{{SOM}_{slow}^{C:N}}\chi_{{SOM}_{slow}^{N:P}}} C_{slow} - P_{slow}}{\tau_{slow}})$$

### Phosphorus adsorption and (ab)sorption

desorption is modelled as: $$\label{eq:sb:pi_desorp}
\begin{align}
F_{desorp,PO_4} &=f(T_{soil},E_{a,abs}) \times  k_{abs}   \times P_{lab} - f(T_{soil},E_{a,des})   \times k_{des}   \times P_{sorb},\ where\\
&f(T_{soil,E_a}) = e^{-\frac{E_{a}}{R} \times (\frac{1}{T_{soil}} - \frac{1}{T^{\tau}_{ref}})}
\end{align}$$ where $k_{abs}$ and $k_{des}$ are the rate constants of
(ab)sorption and desorption, and $E_{a,abs}$ and $E_{a,des}$ the
respective activation energies.

The adsorption ($F_{adsorp,PO_4}$) flux from soil solution to the soil
adsorption sites is calculated assuming constant Langmuir equilibrium
(Barrow 1978) between soluble and adsorbed P: $$\label{eq:sb:langmuir}
\begin{align}
PO_4 &= \frac{S_{max} \times PO_4}{K_S + PO_4},\ thus\\
&\frac{\partial P_{lab}}{\partial t} = \frac{S_{max}\times PO_4}{(K_S + PO_4)^2} \frac{\partial PO_4}{\partial t},
\end{align}
\text{by rearranging Eq. \ref{eq:sb:langmuir}b}
\begin{align}
&\frac{\partial P_{lab}}{\partial t} = k_{p}\frac{\partial {(P_{lab}+PO_4)}}{\partial t}\\
&\frac{\partial PO_4}{\partial t} = (1-k_{p})\frac{\partial {(P_{lab}+PO_4)}}{\partial t},\ where\\
&k_p = \frac{S_{max}\times PO_4}{(K_S + PO_4)^2+S_{max}PO_4},
\end{align}$$

where $S_{max}$ and $K_S$ are the maximum sorption capacity, and the
half-saturation concentration coefficient of the soil, and are modified
by soil moisture and SOM content as follows:

$$\label{eq:sb:Psorbpar}
\begin{align}
S_{max} &= \Theta_{soil} \times (S^{max}_{om} V^{frac}_{om} \rho_{org}^{bulk} + S^{max}_{mineral} V^{frac}_{mineral} \rho_{soil}), \ and\\
K_S &=  K^{sorb}_{m,om} V^{frac}_{om} \rho_{org}^{bulk} + K^{sorb}_{m,mineral} V^{frac}_{mineral} \rho_{soil}
\end{align}$$

where $V^{frac}_{om}$ and $V^{frac}_{mineral}$ are volumetric fractiosn
of organic matter and fine soil minerals, respectively. $S^{max}_{om}$
and $S^{max}_{mineral}$ are the maximum sorption capacity of pure
organic matter and pure fine soil, respectively. $K^{sorb}_{m,om}$ and
$K^{sorb}_{m,mineral}$ are the half-saturation concentration coefficient
of pure organic matter and pure fine soil, respectively.

Based on
Eq.<a href="#eq:sb:simple_im_dynamics" data-reference-type="ref"
data-reference="eq:sb:simple_im_dynamics">[eq:sb:simple_im_dynamics]</a>f
and Eq.<a href="#eq:sb:pi_dynamics" data-reference-type="ref"
data-reference="eq:sb:pi_dynamics">[eq:sb:pi_dynamics]</a>a, the
equilibrium in Eq.<a href="#eq:sb:langmuir" data-reference-type="ref"
data-reference="eq:sb:langmuir">[eq:sb:langmuir]</a> could be solved .
$$\begin{split}
\frac{\partial {(P_{lab}+PO_4)}}{\partial t} &= F_{dep,PO_4} + F_{weath,PO_4} + F_{biomin,PO_4} - U_{plant,PO_4} - U_{myc,PO_4} - F_{desorp,PO_4}\\
&- \sum(\Phi_{i,PO_4}) -  \frac{\partial v_{PO_4} PO_4}{\partial z} + \frac{\partial}{\partial z} (D_b \frac{\partial P_{lab}}{\partial z})\\
\end{split}$$

### Soil in- and outfluxes

Currently, gas diffusion is not modelled explicitly. Instead, is assumed
to be directly released to the atmosphere. The carbon efflux per soil
layer is described as:

$$F^{\uparrow}_{CO_2} = ( (1-\eta_{C,litter \rightarrow fast}) (\frac{C_{met}}{\tau_{met}} + \frac{C_{str}}{\tau_{str}}) + (1 - \eta_{C,wl \rightarrow met|str}) \frac{C_{wl}}{\tau_{wl}} + (1 - \eta_{C,fast \rightarrow slow}) \frac{C_{fast}}{\tau_{fast}} + (1 - \eta_{C,slow \rightarrow fast}) \frac{C_{slow}}{\tau_{slow}} ) \Delta t,$$

and similar for and fluxes. Slightly differently from that, the emission
of gaseous N species is assumed to follow (Xu-Ri and Prentice 2008),
which considers the effect of temperature and moisture on gas loss.
However, transfer between soil layers is equally not treated explicitly.

$$\label{eq:sb:simple_ngasloss}
\begin{align}
F^{\uparrow}_{X} &= f(T_{soil}) \times afps \times X, where \\
f(T_{soil}) &= e^{-\frac{E_{a,diff}}{R} \times (\frac{1}{T_{soil}} - \frac{1}{T^{\tau}_{ref}})}
\end{align}$$

and $afps$ is the air-filled pore fraction of the soil (see Eq.
<a href="#eq:sb:simple_anaerobic_vf" data-reference-type="ref"
data-reference="eq:sb:simple_anaerobic_vf">[eq:sb:simple_anaerobic_vf]</a>.)

## Isotopic composition and fractionation

The carbon (C) and nitrogen (N) flows and pool tracked in the model are
comprising all major isotopes (i.e. C = + + , and N = + ). The model
explicitly tracks the mass flow of , and as separate entities for all
biogeochemical pools and fluxes. The molar mixing ratio ($R_X$) of the
isotope (e.g. to the main element (e.g. ) of each biogeochemical pool
can be calculated as

$$R_{13C} = \frac{^{13}C}{^{12}C} = \frac{^{13}C}{C - ^{13}C}$$

and by convention

$$\label{eq:isotope:delta_notation}
\delta_{13C} = ( \frac{R_{13C}}{ R_{ref,13C}} - 1 ) \times 1000$$

where $R_{ref,13C}$ is the reference isotopic molar mixing ratio and
$\delta$ is in . Similar calculations are done for the ratio of to and
the reference value $R_{ref,15N}$. By convention, the delta notation of
is dependent on the content, see (Levin et al. 2010), and the molar mass
of is ignored in the calculations of because of the extremely low
concentrations.

Biogeochemical processes discriminate against the heavier isotope, and
this fractionation process is treated by calculating the mixing ratio of
the isotope of the resulting flux as

$$\label{eq:isotope:fractionation}
R_{sink} = \frac{R_{source}}{\frac{\epsilon_{process}}{1000}+1}$$

where $R_{source}$ is the molar mixing ratio of the source pool of the
reaction, $R_{sink}$ is the molar mixing ratio of the resulting matter
flux, and $\epsilon_{process}$ is a process and isotope specific
discrimination rate.

Isotopic discrimination of and by photosynthesis is modelled according
to the general equation derived by (Farquhar, OLeary, and Berry 1982;
Drake 2014), so that

$$\label{eq:ps_discrimination}
D_x = a_x + ( c_x + \phi_{C4} \times b_x - a_x )  \frac{c_i}{c_a}$$

where $a_x$, and $b_x$ are isotope-specific constants ( and ,
respectively, see Table
<a href="#tab:assimi:param" data-reference-type="ref"
data-reference="tab:assimi:param">[tab:assimi:param]</a>). $c_x$ and
$\phi_{C4}$ account for the additional bundle-sheath processes in C4
plants. For C3 plants, these processes do not play a role and these
parameters are 0 and 1, respectively. In the model, currently only
photosynthesis is assumed to result in C-isotope discrimination,
ignoring the effect of the smaller and uncertain discrimination by
tissue construction, storage formation and respiration (Brüggemann et
al. 2011).

Isotopic discrimination for various nitrogen cycle processes (biological
nitrogen fixation, ammonification, plant and microbial N uptake, and
processes associated with nitrification and denitrification) are taken
from (Robinson 2001). According to (Robinson 2001), in case of
near-complete consumption of the source pool, the discrimination is
reduced as

$$\epsilon_{process}^{corr} = \epsilon_{process} (f_{source}-1)log(\frac{1-f_{source}}{f_{source}})$$

where $f_{source}$ is the ratio of the source consumption to the source
pool size.

## Radiation, surface energy balance and soil hydrology

### Net surface shortwave radiation budget

Canopy radiation interception is calculated with a multi-layer scheme
following (Spitters 1986), with radiation levels calculated at the
mid-point of each canopy layer. The scheme uses up to 20 (default 10)
canopy layers, with exponentially increasing layer thickness as the
canopy depth increases. The original scheme has been extended to
diagnose canopy albedo, to account for clumping (see eq.
<a href="#eq:rad:omega" data-reference-type="ref"
data-reference="eq:rad:omega">[eq:rad:omega]</a>), and to approximate
the attenuation of the shortwave radiation back-scatter from the soil to
allow for a smooth transition of surface albedo from soil to vegetation
values with incrasing leaf coverage.

The scheme is applied separately to the visible ($vis$) and
near-infrared ($nir$) radiation band, where the parameterisation of the
visible radiation is based on the assumption that the radiation
interception and reflection is similar to that of the photosynthetically
active range (i.e. 400-700 nm). In the following, the subscripts for
visible and near-infrared are omitted for readability.

The reflection coefficient ($\rho_{rad}$) of the green canopy is given
by: $$\label{eq:rad:leaf_albedo}
\rho_{rad} = \frac{1-\sqrt{1-\sigma}}{1+\sqrt{1-\sigma}}\times \frac{2}{1+ \rho^{sbeta} \times cos(\gamma^{*})},$$

where the first term on the right hand side is the reflection of a
horizontally oriented canopy, and the second term empirically adjusts
the reflection to a spherical distribution. $\sigma$ is the PFT-specific
single leaf scattering coefficient and $\rho^{sbeta}$ is a conversion
constant. Because all equations for leaf reflection and absorption
coefficients are only valid for high solar elevation, the true zenith
angle ($\gamma$) is constrained to values larger than 10$^\circ$
($\gamma^{*}$). Note that (Spitters 1986) use the sine of the solar
elevation angle.

For a spherical leaf angle distribution with leaves distributed randomly
within the canopy volume, the extinction coefficients of the diffuse
flux ($k_{df}$) and that of the direct component of the direct flux
($k_{bl}$) are approximated, respectively, by:

$$\label{eq:rad:extinction_coef}
\begin{align}
k_{df} &=  k_{df,0} \sqrt{1-\sigma} \times \Omega\\
k_{bl} &= \frac{k_{bl,0} \times \Omega}{cos(\gamma^{*})}
\end{align}$$

where $\Omega$ is the clumping index according to (Campbell and Norman
1998), which is calculated as:

$$\label{eq:rad:omega}
\Omega = \Omega_{0} / (\Omega_{0} + (1 - \Omega_{0}) \times e^{-k_{csf} \times acos(\gamma^{*})^{\phi_{crown}}}),$$

where $\Omega_{0}$ and $\phi_{crown}$ are the PFT-specific clumping
factor at nadir and crown shape factor, respectively, and $k_{csf}$ is a
correction factor.

Light levels decrease exponentially in the canopy, such that the
attenuation of direct ($dr$) and diffuse ($df$) top-of-the-canopy
irradiance ($I_{dr\downarrow,0}$ and $I_{df\downarrow,0}$, respectively)
at any cumulative leaf area index ($LAI_c$; from the top) is given by:

$$\label{eq:rad:light_profile}
\begin{align}
&I_{df \downarrow} = (1-\rho_{rad}) \times I_{df \downarrow,0} \times e^{-k_{df}\times LAI_{c}}\\
&I_{dr \downarrow} = (1-\rho_{rad}) \times I_{dr \downarrow,0} \times e^{-\sqrt{1-\sigma} \times k_{bl}\times LAI_{c}}\\
&I_{dr,dr \downarrow} = I_{dr \downarrow,0} \times e^{- k_{bl}\times LAI_{c}}\label{eq:rad:idrdr}\\
&I_{dr,df \downarrow} = I_{dr \downarrow} - I_{dr,dr \downarrow}\\
\end{align}$$

where $I_{dr,df \downarrow}$ is the diffuse part of the direct beam
resulting from scattering of the direct beam and $I_{dr,dr \downarrow}$
is the direct beam remaining direct at the canopy depth $LAI_{c}$.

Below the canopy ($bc$), i.e. at the soil surface, the downwelling
energy flux ($I_{bc \downarrow}$) is divided into a part that is
absorbed by the soil ($I_{a,soil}$) and a part that is backscattered as
diffuse radiation ($I_{bc \uparrow}$), depending on the soil’s albedo
($alb_{soil}$): $$\label{eq:rad:soil_netbalance}
\begin{align}
&I_{bc \downarrow} = (1-\rho_{rad}) \times ( I_{df \downarrow,0} \times e^{-k_{df}\times LAI} + I_{dr \downarrow,0} \times e^{-\sqrt{1-\sigma} \times k_{bl}\times LAI})\\
&I_{a,soil} = (1-alb_{soil}) \times I_{bc \downarrow} \\
&I_{bc \uparrow} = alb_{soil} \times I_{bc \downarrow}
\end{align}$$

To first-order, the diffuse light profile of the canopy accounting for
the backscatter of diffuse radiation from the soil can be approximated
as $$\label{eq:rad:omega2}
I_{df \downarrow \uparrow} = (1-\rho_{rad}) \times (I_{df \downarrow,0} \times e^{-k_{df}\times LAI_{c}} + \times I_{bc \uparrow} \times e^{-k_{df} \times (LAI-LAI_{c}))}$$

Absorption ($A^{rad}$) is taken to be complementary to transmission,
therefore the absorbed diffuse and direct energy flux at a cumulative
leaf area $LAI_c$ can be written as:
$$\label{eq:rad:canopy_absorbed_radiation_components}
\begin{align}
&A^{rad}_{df} = k_{df} \times I_{df \downarrow \uparrow}\\
&A^{rad}_{dr} = (1-\sigma) \times k_{bl} \times I_{dr}\\
&A^{rad}_{dr,dr} = (1-\sigma) \times k_{bl} \times I_{dr \downarrow,0} \times e^{-k_{bl}\times LAI_c}\label{eq:rad:idrdra}\\
&A^{rad}_{dr,df} = A^{rad}_{dr} - A^{rad}_{dr,dr}
\end{align}$$

The canopy is then split into a sunlit and a shaded part, with the
sunlit fraction defined as: $$\label{eq:rad:fsunlit}
f_{sunlit} = e^{-k_{bl}\times LAI_c}$$

following eq. <a href="#eq:rad:idrdr" data-reference-type="ref"
data-reference="eq:rad:idrdr">[eq:rad:idrdr]</a>. The sunlit part
receives both diffuse and direct radiation, whereas the shaded part only
received diffuse radiation. Thus, $$\label{eq:rad:a_shadesunlit}
\begin{align}
&A^{rad}_{shaded} = A^{rad}_{df} + A^{rad}_{dr,df}\\
&A^{rad}_{sunlit} = A^{rad}_{shaded} + (1-\sigma) \times k_{bl} \times I_{dr \downarrow,0} 
\end{align}$$

The canopy albedo is diagnosed (rather than simply taken as
$\rho_{rad}$) as:

$$\label{eq:rad:canopy_albedo}
alb_{can} = 1 - \frac{I_{bc \downarrow} +  (1-f_{sunlit}) \times A^{rad}_{shaded} + f_{sunlit} \times A^{rad}_{sunlit} }{I_{df\downarrow,0}+I_{dr \downarrow, 0}}$$

The total shortwave upward flux is diffuse and calculated as the
backscattered flux of the canopy plus the backscattered flux from the
soil, which is transmitted diffusely through the canopy:

$$\label{eq:rad:idfup}
I_{df \uparrow} = alb_{can} \times (I_{df\downarrow,0}+I_{dr \downarrow, 0} - I_{bc \downarrow}) + (1 - \rho_{rad}) \times I_{bc \uparrow} \times e^{-k_{df}\times LAI}$$

Based on this, the total surface albedo $alb_{surf}$ (i.e. the albedo
derived from vegetation and soil radiation transfer, absorption, and
reflection) and net shortwave flux $I_{net}$ can be calculated as:

$$\label{eq:rad:canopy_absorbed_radiation_components2}
\begin{align}
&alb_{surf} =  \frac{I_{df \uparrow}}{I_{df\downarrow,0}+I_{dr \downarrow, 0}} \\
&I_{net} = (1-alb_{surf}) \times (I_{df\downarrow,0}+I_{dr \downarrow, 0})
\end{align}$$

### Surface energy balance

The representation of the surface energy balance including the turbulent
momentum and heat exchanges, surface, and soil temperature calculations
follows largely the scheme of JSBACH 3, as described by (Roeckner et al.
2003). The net radiation ($R_{net}$) at the surface consists of the
following components:
$$R_{net} = (1 - \alpha_s ) R_{sd} + R_{ld} - \epsilon \sigma {T_s}^4$$

where the term $\alpha_s \times R\_{sd}$ is the net surface shortwave
balance, denoted as $I_{net}$ in Eq.
<a href="#eq:rad:canopy_absorbed_radiation_components2"
data-reference-type="ref"
data-reference="eq:rad:canopy_absorbed_radiation_components2">[eq:rad:canopy_absorbed_radiation_components2]</a>,
$R_{ld}$ the downwelling longwave radiation, $\epsilon$ the surface
emissivity, $\sigma$ the Stefan-Boltzmann constant, and $T_s$ is the
surface temperature.

The surface energy balance can then be written as
$$\label{eq:sse_surface_energy_balance}
C_{s,l} \frac{\partial T_s}{\partial t} = R_{net} + LE + H + G$$

where $C_{s,l}$ is the heat capacity of the surface layer, $H$ is the
sensible heat flux, $LE$ the latent heat flux, composed of interception
loss $E_i$, soil evaporation $E_s$ and transpiration $E_v$, as described
in Section
<a href="#sec:sse:surface_hydrology" data-reference-type="ref"
data-reference="sec:sse:surface_hydrology">1.6.3</a>. $G$ is the ground
heat flux, which is obtained from the solution of the thermal diffusion
equation, which is used to diagnose the temperature profile within the
soil $$\label{eq:sse:ground_heat_flux}
C_s\frac{\partial T}{\partial t}=-\frac{\partial G}{\partial z}=-\frac{\partial}{\partial z}\left(-\lambda_s\frac{\partial T}{\partial z}\right)$$

where $C_s$ is the volumetric heat capacity of the soil, $G$ is the
thermal heat flux (positive downward), $\lambda_s=C_s \kappa_s$ is the
thermal conductivity, $\kappa_s$ the thermal diffusivity, and $z$ the
depth. This equation is solved down to a prescribed depth assuming zero
flux conditions at the bottom and surface temperature $T_s$ at the top
as obtained from Eq.
<a href="#eq:sse_surface_energy_balance" data-reference-type="ref"
data-reference="eq:sse_surface_energy_balance">[eq:sse_surface_energy_balance]</a>.

The sensible heat flux $H$ is obtained from the temperature gradient
according to

$$H=\rho_{air} C_p\frac{T_{air}-T_s}{r_a}\label{eq:sensibleHeatFlux}$$

where $\rho_{air}$ and $C_p$ are density and heat capacity of the
surface atmospheric layer, and $T_{air}$ and $T_s$ are air and surface
temperature. $r_a$ is the aerodynamic resistance.

$$r_{a}=(C_{h}\left|v_{h}\right|)^{-1}\label{eq:aerodynResistance}$$ and
depends on the transfer coefficient for heat $C_{h}$ (Roeckner et al.
2003). and the absolute value of horizontal wind velocity $v_{h}$.

### Surface and soil hydrology

The soil model keeps track of water in terms of amount of water
($W_{skin|soil}$; m), which in the soil can be converted to the layers
fractional water content ($\Theta$) as: $$\label{eq:sse:theta}
\Theta_{soil,sl} = \frac{W_{soil,sl}}{dz_{sl}}$$

with a soil water matrix potential ($\Psi_{soil,sl}$) derived from a
pedotransfer function (Saxton and Rawls 2006) as follows:
$$\label{eq:sse:soil_matrix_potential}
\Psi_{soil,sl} = k^{\Psi}_A \Theta_{soil,sl}^{k^{\Psi}_B}$$

where $k^{\Psi}_A$, and $k^{\Psi}_B$ are soil texture dependent
parameters (Saxton and Rawls 2006).

Precipitation ($Precip$) becomes partially intercepted by the canopy
before reaching the soil: $$\label{eq:sse_rain_interception}
F_{inter} = MIN( k_{eff,inter} \times Precip , \frac{w_{skin,max} \times LAI - W_{skin}}{dt})$$

Different from (Roeckner et al. 2003), the remaining throughfall is
partitioned in the first soil layer into infiltration ($F_{inf}$),
within the limits of its field capacity ($\Theta_{fc}$), surface runoff
($F_{runoff}$), and a fraction that is leaked to the next layer along
preferential flow pathways ($F_{pref}$):

$$\label{eq:sse:surface_water_balance}
\begin{align}
F_{input} &= MIN(Precip - F_{inter} ,\frac{W_{fc,sl=1} - W_{act,sl=1}}{dt}),  \\
F_{inf} &= (1 - \frac{k_{pref}}{dz_{sl}}) F_{input}\\
F_{runoff} &= (1 - k_{pref,runoff}) (F_{through} - F_{input}), and \\
F_{pref,sl=1} &= \frac{k_{pref}}{dz_{sl}} F_{inf} + k_{pref,runoff} (F_{through} - F_{input})
\end{align}$$

Different from (Roeckner et al. 2003), the diffusive flux between two
layers of depth $dz$ is given by the Richards-equation:

$$\label{eq:sse:water_diffusion}
\begin{align}
F_{diff,sl-1 \rightarrow sl} &= \frac{\bar{K}_{diff,sl-1 \rightarrow sl}}{0.5\times(dz_{sl-1}+dz_{sl})} (\Psi_{soil,sl} - \Psi_{soil,sl-1}), where \\
\bar{K}_{diff,sl-1 \rightarrow sl} &= \frac{dz_{sl-1} K_{diff,sat,sl-1} (\frac{\Theta_{act,sl-1}}{\Theta_{sat,sl-1}})^{k^{diff}_{C,sl}} + dz_{sl} K_{diff,sat,sl} (\frac{\Theta_{act,sl}}{\Theta_{sat,sl}})^{k^{diff}_{C,sl}} }{dz_{sl-1} + dz_{sl}}
\end{align}$$

where $K_{diff,sat,sl}$ and ${k^{diff}_{C,sl}}$ are derived from soil
texture according to (Saxton and Rawls 2006).

$$E_{i}=\rho_{air}\frac{q_{a}-q_{s}(T_{s},p_{s})}{r_{a}} \frac{W_{skin}}{w_{skin,max} \times LAI},$$

where $\rho_{air}$ is air density, $q_{a}$ specific humidity in lowest
atmospheric level, $q_{s}$ saturation specific humidity at surface
temperature $T_{s}$ and pressure $p_{s}$, and ratio of the size of the
canopy water storage ($W_{skin}$) and the maximum size of this
reservoir, where $w_{skin,max}$ is a parameter, and $LAI$ is the current
leaf area index of the vegetation.

Evaporation from the soil surface is calculated as:

$$E_{s}=\rho_{air}\frac{q_{a}-q_{s}(T_{s},p_{s})}{r^*_{a}} \Theta_{soil,1} (1-exp(k_{fpc} \times LAI))$$

where $\Theta_{soil,1}$ the fractional soil water content of the first
soil layer. The term $(1-exp(k_{fpc} \times LAI))$ has been added to the
model of (Roeckner et al. 2003) accounts for the reduced energy
available for evaporation underneath a closed canopy.

Transpiration from the dry vegetation surfaces is
$$E_{v}=\rho_{air}\frac{q_{a}-q_{s}(T_{s},p_{s})}{r_{a}+r}$$ where
stomatal resistance $r$ of the canopy is the inverse of the stomatal
conductance of the canopy (see Eq.
<a href="#eq:stomatal_coupling" data-reference-type="ref"
data-reference="eq:stomatal_coupling">[eq:stomatal_coupling]</a>).

The partitioning of the transpiration flux across soil layers
$f_{trans,sl}$ is calculated based on the layered soil water potential
($\Psi_{soil,sl}$), the fractional root distribution ($f_{root,sl}$), as
well as a PFT-specific minimum soil water potential,

$$\begin{align}
f_{trans,sl} &= \frac{f_{root,sl} \times \beta_{gs,soil,sl}}{ \sum\limits_{i=1}^n f_{root,i} \times \beta_{gs,soil,i}} , where \\
 \beta_{gs,soil,sl} &= \Psi_{soil,sl} / \Psi_{soil,min}
\end{align}$$

Based on these fluxes, the water budget can be described as follows

$$\begin{align}
\frac{W_{skin}}{dt} &= F_{inter} - E_{i} \\
\frac{W_{soil,sl=1}}{dt} &= F_{inf} - E_{a} - f_{frans,sl=1} E_{v} - F_{diff,sl=1 \rightarrow sl=2} \\
\frac{W_{soil,sl=2,n}}{dt} &= (1 - \frac{k_{pref}}{dz_{sl}}) F_{pref,sl-1} - f_{frans,sl} F_{trans} + F_{diff,sl-1 \rightarrow sl} - F_{diff,sl \rightarrow sl+1} \\
\end{align}$$

Saturation of a soil layer leads to increased percolation to the next
lower layer. The lower boundary layer in the soil is modelled as a zero
gradient boundary, i.e. only percolation flow leads to drainage from the
soil column. The net water transport between layers is used as input to
the vertical flow of soluable biogeochemical pools.

<div class="figure*">

<embed src="figures/Figure1.eps" style="width:12cm" />

</div>

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Ahrens:2015ie" class="csl-entry">

Ahrens, Bernhard, Maarten C Braakhekke, Georg Guggenberger, Marion
Schrumpf, and Markus Reichstein. 2015. “Contribution of Sorption, DOC
Transport and Microbial Interactions to the 14C Age of a Soil Organic
Carbon Profile: Insights from a Calibrated Process Model.” *Soil Biology
and Biochemistry* 88 (C): 390–402.

</div>

<div id="ref-Atkin:2014" class="csl-entry">

Atkin, Owen K., Patrick Meir, and Matthew H. Turnbull. 2014. “Improving
Representation of Leaf Respiration in Large-Scale Predictive
Climate–Vegetation Models.” *New Phytologist* 202 (3): 743–48.
<https://doi.org/10.1111/nph.12686>.

</div>

<div id="ref-Ball:1987uv" class="csl-entry">

Ball, JT, IE Wooddrow, and Joseph A Berry. 1987. “<span class="nocase">A
model predicting stomatal conductance and its contribution to the
control of photosynthesis under different environmental
conditions</span>.” In *Progress in Photosynthesis Research*, edited by
J Biggens, 1–4. Dordrecht, The Netherlands.

</div>

<div id="ref-Barrow:1978" class="csl-entry">

Barrow, N. J. 1978. “The Description of Phosphate Adsorption Curves.”
Journal Article. *Journal of Soil Science* 29: 447–62.

</div>

<div id="ref-Bernacchi:2001wv" class="csl-entry">

Bernacchi, C J, E L Singsaas, C Pimentel, A R Protis JR, and S P Long.
2001. “<span class="nocase">Improved temperature response functions for
models of Rubisco-limited photosynthesis</span>.” *Plant, Cell &
Environment* 24 (January): 253–59.

</div>

<div id="ref-Bonan:2015" class="csl-entry">

Bonan, Gordon. 2015. *Ecological Climatology: Concepts and
Applications*. 3rd ed. Cambridge University Press.
<https://doi.org/10.1017/CBO9781107339200>.

</div>

<div id="ref-Bruggemann:2011co" class="csl-entry">

Brüggemann, N, A Gessler, Z Kayler, S G Keel, F Badeck, M Barthel, P
Boeckx, et al. 2011. “<span class="nocase">Carbon allocation and carbon
isotope fluxes in the plant-soil-atmosphere continuum: a review</span>.”
*Biogeosciences* 8 (11): 3457–89.

</div>

<div id="ref-Campbell:1998vb" class="csl-entry">

Campbell, G S, and J M Norman. 1998. *<span class="nocase">An
Introduction to Environmental Biophysics</span>*. Springer.

</div>

<div id="ref-Collatz:1992us" class="csl-entry">

Collatz, G James, M Ribas-Carbo, and Joseph A Berry. 1992.
“<span class="nocase">Coupled photosynthesis-stomatal conductance model
for leaves of C4 plants</span>.” *Australian Journal of Plant
Physiology* 19 (May): 519–38.

</div>

<div id="ref-Drake:2014dl" class="csl-entry">

Drake, Brandon Lee. 2014. “<span class="nocase">Using Models of Carbon
Isotope Fractionation during Photosynthesis to Understand the Natural
Fractionation Ratio</span>.” *Radiocarbon* 56 (1): 29–38.

</div>

<div id="ref-Farquhar:1982ks" class="csl-entry">

Farquhar, G D, M H OLeary, and J A Berry. 1982. “<span class="nocase">On
the Relationship Between Carbon Isotope Discrimination and the
Intercellular Carbon Dioxide Concentration in Leaves</span>.”
*Functional Plant Biology* 9 (2): 121–37.

</div>

<div id="ref-Friend:2010je" class="csl-entry">

Friend, Andrew D. 2010. “<span class="nocase">Terrestrial plant
production and climate change</span>.” *Journal of Experimental Botany*
61 (5): 1293–1309.

</div>

<div id="ref-Friend:2009tf" class="csl-entry">

Friend, Andrew D, R J Geider, M J Behrenfeld, and C J Still. 2009.
“<span class="nocase">Photosynthesis in Global-Scale Models</span>.” In
*Photosynthesis in Silico: Understanding Complexity from Molecules to
Ecosystems*, edited by A Laisk, L Nedbal, and Govindjee, 465–97.
Springer.

</div>

<div id="ref-Jackson:1996tr" class="csl-entry">

Jackson, R B, J Canadell, J R Ehleringer, H A Mooney, O E Sala, and E D
Schulze. 1996. “<span class="nocase">A global analysis of root
distributions for terrestrial biomes</span>.” *Oecologia* 108 (3):
389–411.

</div>

<div id="ref-June:2004kj" class="csl-entry">

June, Tania, John R Evans, and Graham D Farquhar. 2004.
“<span class="nocase">A simple new equation for the reversible
temperature dependence of photosynthetic electron transport: a study on
soybean leaf</span>.” *Functional Plant Biology* 31 (3): 275.

</div>

<div id="ref-Knauer:2015dh" class="csl-entry">

Knauer, Jürgen, Christiane Werner, and Sönke Zaehle. 2015.
“<span class="nocase">Evaluating stomatal models and their atmospheric
drought response in a land surface scheme: A multibiome
analysis</span>.” *Journal of Geophysical Research: Biogeosciences* 120
(10): 1894–1911.

</div>

<div id="ref-Krinner:2005wg" class="csl-entry">

Krinner, G, Nicolas Viovy, Nathalie de Noblet-Ducoudré, Jérôme Ogée, Jan
Polcher, Pierre Friedlingstein, Philippe Ciais, Stephen A Sitch, and I C
Prentice. 2005. “<span class="nocase">A dynamic global vegetation model
for studies of the coupled atmosphere-biosphere system</span>.” *Global
Biogeochemical Cycles* 19 (1): GB1015.

</div>

<div id="ref-Kull:1998cr" class="csl-entry">

Kull, O, and B Kruijt. 1998. “<span class="nocase">Leaf photosynthetic
light response: a mechanistic model for scaling photosynthesis to leaves
and canopies</span>.” *Functional Ecology* 12 (5): 767–77.

</div>

<div id="ref-Levin:2010fu" class="csl-entry">

Levin, Ingeborg, Tobias Naegler, Bernd Kromer, Moritz Diehl, Roger
Francey, Angel Gomez-Pelaez, Paul Steele, Dietmar Wagenbach, Rolf
Weller, and Douglas Worthy. 2010. “<span class="nocase">Observations and
modelling of the global distribution and long-term trend of atmospheric
</span>.” *Tellus Series B-Chemical and Physical Meteorology* 62 (1):
26–46.

</div>

<div id="ref-Lloyd1994-jn" class="csl-entry">

Lloyd, John, and J A Taylor. 1994. “On the Temperature Dependence of
Soil Respiration.” *Functional Ecology*, 315–23.
<https://www.jstor.org/stable/2389824>.

</div>

<div id="ref-Makela:2004" class="csl-entry">

Mäkelä, Annikki, Pertti Hari, Frank Berninger, Heikki Hänninen, and Eero
Nikinmaa. 2004. “Acclimation of Photosynthetic Capacity in Scots Pine to
the Annual Cycle of Temperature.” *Tree Physiology* 24 (4): 369–76.
<https://doi.org/10.1093/treephys/24.4.369>.

</div>

<div id="ref-Meyerholt:2015hh" class="csl-entry">

Meyerholt, Johannes, and Sönke Zaehle. 2015. “<span class="nocase">The
role of stoichiometric flexibility in modelling forest ecosystem
responses to nitrogen fertilization</span>.” *New Phytologist*, July,
n/a–.

</div>

<div id="ref-Niinemets:1998" class="csl-entry">

Niinemets, Ülo, Olevi Kull, and John D. Tenhunen. 1998. “An Analysis of
Light Effects on Foliar Morphology, Physiology, and Light Interception
in Temperate Deciduous Woody Species of Contrasting Shade Tolerance.”
*Tree Physiology* 18 (10): 681–96.
<https://doi.org/10.1093/treephys/18.10.681>.

</div>

<div id="ref-Parton:1993ud" class="csl-entry">

Parton, W J, J M O Scurlock, D S Ojima, T G Gilmanov, R J Scholes, D S
Schimmel, T Kirchner, et al. 1993. “<span class="nocase">Observations
and modelling of biomass and soil organic matter dynamics for the
grassland biome worldwide</span>.” *Global Biogeochemical Cycles* 7 (4):
785–809.

</div>

<div id="ref-Robinson2001-hy" class="csl-entry">

Robinson, David. 2001. “$\delta$<!-- -->15N as an Integrator of the
Nitrogen Cycle.” *Trends in Ecology & Evolution* 16 (3): 153–62.
<https://doi.org/10.1016/S0169-5347(00)02098-X>.

</div>

<div id="ref-Roeckner:2003wc" class="csl-entry">

Roeckner, E, A Tompkins, G Bäuml, L Bonaventura, R Brokopf, M Esch, M
Giorgetta, et al. 2003. “<span class="nocase">The atmospheric general
circulation model ECHAM5 - Part I: Model description</span>.” 349.
Hamburg: Max Planck Institute for Meteorology.

</div>

<div id="ref-Rogers:2017np" class="csl-entry">

Rogers, Alistair, Belinda E Medlyn, Jeffrey S Dukes, Gordon B Bonan,
Susanne von Caemmerer, Michael C Dietze, J Kattge, et al. 2017.
“<span class="nocase">A roadmap for improving the representation of
photosynthesis in Earth system models</span>.” *New Phytologist* 213
(1): 22–42.

</div>

<div id="ref-Saxton:2006jb" class="csl-entry">

Saxton, K E, and W J Rawls. 2006. “<span class="nocase">Soil Water
Characteristic Estimates by Texture and Organic Matter for Hydrologic
Solutions</span>.” *Soil Science Society of America Journal* 70 (5):
1569–10.

</div>

<div id="ref-Sitch:2003wg" class="csl-entry">

Sitch, Stephen A, Benjamin Smith, I C Prentice, Almut Arneth, Alberte
Bondeau, W Cramer, J O Kaplan, et al. 2003.
“<span class="nocase">Evaluation of ecosystem dynamics, plant geography
and terrestrial carbon cycling in the LPJ Dynamic Global Vegetation
Model</span>.” *Global Change Biology* 9: 161–85.

</div>

<div id="ref-Spitters:1986tv" class="csl-entry">

Spitters, C J T. 1986. “<span class="nocase">Separating the Diffuse and
Direct Component of Global Radiation and Its Implications for Modeling
Canopy Photosynthesis .2. Calculation of Canopy Photosynthesis</span>.”
*Agricultural and Forest Meteorology* 38 (1-3): 231–42.

</div>

<div id="ref-Wang:2010ba" class="csl-entry">

Wang, Y P, R M Law, and B Pak. 2010. “<span class="nocase">A global
model of carbon, nitrogen and phosphorus cycles for the terrestrial
biosphere</span>.” *Biogeosciences* 7 (7): 2261–82.

</div>

<div id="ref-White:2000vl" class="csl-entry">

White, M A, P E Thornton, S Running, and R Nemani. 2000.
“<span class="nocase">Parameterization and Sensitivity Analysis of the
BIOME-BGC Terrestrial Ecosystem Model: Net Primary Production
Controls</span>.” *Earth Interactions* 4: 1–55.

</div>

<div id="ref-XURI:2008jq" class="csl-entry">

Xu-Ri, and I C Prentice. 2008. “<span class="nocase">Terrestrial
nitrogen cycle simulation with a dynamic global vegetation
model</span>.” *Global Change Biology* 14 (8): 1745–64.

</div>

<div id="ref-Zaehle:2011gy" class="csl-entry">

Zaehle, S, Philippe Ciais, Andrew D Friend, and Vincent Prieur. 2011.
“<span class="nocase">Carbon benefits of anthropogenic reactive nitrogen
offset by nitrous oxide emissions</span>.” *Nature Geoscience* 4 (8):
601–5.

</div>

<div id="ref-Zaehle:2010a" class="csl-entry">

Zaehle, S, Pierre Friedlingstein, and Andrew D. Friend. 2010.
“Terrestrial Nitrogen Feedbacks May Accelerate Future Climate Change.”
*Geophysical Research Letters* 37 (1).
<https://doi.org/10.1029/2009GL041345>.

</div>

<div id="ref-Zaehle:2010cz" class="csl-entry">

Zaehle, S, and Andrew D Friend. 2010. “<span class="nocase">Carbon and
nitrogen cycle dynamics in the O-CN land surface model: 1. Model
description, site-scale evaluation, and sensitivity to parameter
estimates</span>.” *Global Biogeochemical Cycles* 24 (1): GB 1005.

</div>

<div id="ref-Zerihun:1998" class="csl-entry">

Zerihun, Ayalsew, Bruce A McKenzie, and James D Morton. 1998.
“Photosynthate Costs Associated with the Utilization of Different
Nitrogen–Forms: Influence on the Carbon Balance of Plants and Shoot–Root
Biomass Partitioning.” *The New Phytologist* 138 (1): 1–11.

</div>

</div>
