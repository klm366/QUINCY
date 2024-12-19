---
title: Scientific documentation for QUINCY: Overview 
---

### The QUINCY model 

<figure>
  <img src="../../../media/QUINCY-high-res.png" 
  style="float: left; border:thin" width="50%" name="fig_quincy_schematic" 
  alt="QUINCY schematic">
  <figcaption> **Fig 1:** Schematic representation of biogeochemical, water and energy flows in QUINCY</figcaption>
</figure>

**QUINCY (QUantifying Interactions between terrestrial Nutrient
  CYcles and the climate system)** has been developed to include a 
  consistent representation of element cycling in terrestrial
  ecosystems suitable for the modelling of land-atmosphere interactions. 
  QUINCY includes i) a representation of plant
  growth that separates source (e.g. photosynthesis) and sink (growth
  rate of individual tissues, constrained by nutrients, temperature, and
  water availability) processes; ii) the acclimation of many
  ecophysiological processes to meteorological conditions and/or
  nutrient availabilities; iii) an explicit representation of vertical
  soil processes to separate litter and soil organic matter dynamics. QUINCY 
  was originally described by Thum et al. (2019) as a stand-alone ecosystem model. 
  It is now integrated into the ICON-Land infrastructure and shares the physical 
  process representation as well as the software infrastructure. 
  It follows a modular structure to facilitate the consistent testing of multiple 
  hypothesis for one particular process. It also includes a range of new diagnostics 
  (leaf chlorophyll; isotopic tracers) that allow for a more in-depth model evaluation. 
  
  
This page gives a broad overview of 
model features, whereas the detailed description can be found [here](./quincy_long_description.html).

### Overview

The QUINCY model v1.0 represents the coupled biogeochemical cycles of
carbon (C), nitrogen (N), and phosphorus (P) in terrestrial ecosystems,
as well as their interaction with the terrestrial water and energy
balance (<a href="fig_quincy_schematic">[Fig. 1]</a>). The model traces
the flow of these elements as well as C and N isotopes through
vegetation and soil at a half-hourly time-step. The model considers
eight broadly defined plant functional types (PFTs, see Tab.
<a href="./quincy_long_description.html#tab_pft_list">[Table: A1]</a>), representing
different leaf types (herbaceous, broadleaved, coniferous),
phenology (evergreen, rain green, and summer green), as well as growth
form (grasses and trees), with a set of traits describing time-invariant
functional differences across types (see Tab.
<a href="#tab:lctlib:description" data-reference-type="ref"
data-reference="tab:lctlib:description">[tab:lctlib:description]</a>.
Different to other TBMs, certain "soft" plant functional type
characteristics, such as the photosynthetic capacity of leaves or the
temperature sensitivity of photosynthesis and respiration acclimate to
local conditions, and therefore lead to a smoother transition of
ecosystem functioning across eco-climatic gradients. A brief overview of
the model structure and key processes is given in this section, a
detailed mathematical description of the model processes can be found in
Appendix <a href="#sec:appendix" data-reference-type="ref"
data-reference="sec:appendix">[sec:appendix]</a>.

Although conceived to be coupled to the land-surface scheme of an Earth
system model, the model is currently applied stand-alone at the
ecosystem scale. The stand-alone version requires half-hourly short- and
longwave radiation, air temperature, precipitation and snowfall, wind
velocity, atmospheric , , mole fractions, as well as rates of , , and
deposition as time-dependent driving data. In addition, the model
requires input describing the geographical coordinates, vegetation type,
and soil physical and chemical parameters (texture, bulk density,
rooting and soil depth, as well as inorganic soil P content).

### Vegetation processes

Vegetation is characterised by an average individual representative of a
PFT, composed of three fast-lived structural tissue types (pools:
leaves, fine roots, and fruits), a respiring non-structural pool
(labile), as well as a seasonal, non-respiring, and non-structural
storage pool (reserve). In the case of tree vegetation types, three
longer-lived structural tissue types (sapwood, heartwood, and coarse
roots) are represented in addition, as are stand characteristics
(height, diameter, density).

Importantly, plant growth is modelled as the balance of source processes
(photosynthesis, nutrient uptake) and the capacity of the plant to
create sinks (production of new biomass tissue, respiration, and
storage; Sect. <a href="#sec:growth" data-reference-type="ref"
data-reference="sec:growth">[sec:growth]</a>). The capacity of the plant
to generate sinks is dependent on (i) the activity of the plant’s
meristem controlling the usage of the labile pool for growth, which
becomes reduced at low air temperature and/or soil moisture; (ii) the
partitioning of its labile resource to new tissue growth and the
availability of nutrients to create these tissues; (iii) the need for
longer-term storage production; and (iv) respiration for maintenance,
which is given priority over growth of new tissue.

Photosynthesis is calculated for sunlit and shaded leaves separately,
explicitly taking account of the vertical canopy gradient of light,
foliar chlorophyll and photosynthetic N according to (Kull and Kruijt
1998), with extensions for C4 photosynthesis (Friend et al. 2009), and
the temperature response as in (Bernacchi et al. 2001) (Sect.
<a href="#ch:assimi" data-reference-type="ref"
data-reference="ch:assimi">[ch:assimi]</a>). The stomatal conductance is
coupled to the calculation of photosynthesis according to the Ball-Berry
formulation (Ball, Wooddrow, and Berry 1987; Knauer, Werner, and Zaehle
2015). Maintenance respiration is calculated for each tissue type as a
function of its N content (Sect.
<a href="#sec:veg:MaintResp" data-reference-type="ref"
data-reference="sec:veg:MaintResp">[sec:veg:MaintResp]</a>) and tissue
temperature (Lloyd and Taylor 1994). Photosynthesis and respiration
acclimate to growth temperature according to (Friend 2010) and (Atkin,
Meir, and Turnbull 2014), respectively.

Plant nutrient uptake is dependent on fine root biomass density, soil
soluble , , and concentrations as well as plant internal nutrient demand
(Sect. <a href="#sec:sb:simple_uptake" data-reference-type="ref"
data-reference="sec:sb:simple_uptake">[sec:sb:simple_uptake]</a>,
(Zaehle and Friend 2010)). Tissue growth follows a set of allometric
equations describing size-dependent relationships of foliar, fine root,
coarse root, sapwood and fruit biomass (Sect.
<a href="#sec:veg:allocation" data-reference-type="ref"
data-reference="sec:veg:allocation">[sec:veg:allocation]</a>). The
allometric relationship between leaves and fine root responds to N, P
and water limitation and therefore increases uptake capacity under
nutrient limitation.

The allometric equations and the N and P content of each tissue give the
stoichiometric constraint of plant growth. The C:N and N:P ratios of the
leaves and fine roots are flexible within empirical bounds and respond
dynamically to an imbalance between C and N, or N and P availability,
respectively (Sect.
<a href="#sec:veg:stoichem" data-reference-type="ref"
data-reference="sec:veg:stoichem">[sec:veg:stoichem]</a>).
Photosynthesis can be down-regulated through sink-limitation (Sect.
<a href="#sec:veg:sinklim" data-reference-type="ref"
data-reference="sec:veg:sinklim">[sec:veg:sinklim]</a>), if nutrient
shortage leads to the accumulation of carbon in the labile pool.

Short-term fluctuations in the balance between carbon acquisition
through photosynthesis and nutrient acquisition by roots are balanced by
the time-scale of the labile pool turnover. Seasonal and longer-term
fluctuations are balanced to some extend by the build-up and use of the
reserve pool. The flow of compounds between the labile and reserve pools
(Sect. <a href="#sec:veg:reserve" data-reference-type="ref"
data-reference="sec:veg:reserve">[sec:veg:reserve]</a>) is described by
a set of pull (from the reserve pool) and push (to the reserve pool)
mechanisms attempting to simulate the signaling related to the beginning
of the growing season (phenological pull), the need for maintenance and
growth (maintenance pull), as well as to provide resources for the next
growing season or resilience against interannual variability (reserve
push).

Turnover is considered at two time-scales: (i) the rapid breakdown of
enzymes (N) and associated RNA (P) from tissue, and their ensuing
replacement from the labile N and P pools (recycling; at the time scale
of days; (Zaehle and Friend 2010)); and (ii) the senescence of entire
tissues and their conversion to litter, including the resorption of
nutrients from foliage (Sect.
<a href="#veg:sec:turnover" data-reference-type="ref"
data-reference="veg:sec:turnover">[veg:sec:turnover]</a>). Tissue
senescence is calculated given a fixed turnover time for each tissue,
with the exception of the leaves, where turnover is determined by the
PFT-specific phenological timing
<a href="#sec:veg:phen" data-reference-type="ref"
data-reference="sec:veg:phen">[sec:veg:phen]</a>.

Stand-level vegetation dynamics (Sect.
<a href="#veg:sec:dynamics" data-reference-type="ref"
data-reference="veg:sec:dynamics">[veg:sec:dynamics]</a>) are
represented through density-dependent mortality, as well as a dynamic
background mortality related to the overall growth efficiency of the
plant (Sitch et al. 2003; Zaehle and Friend 2010). Re-establishment of
vegetation depends directly on the past vegetation ability to produce
fruits, as the re-establishment is explicitly taking account of the
available seed pool at the forest / grassland floor.

### Soil biogeochemical processes

The vertical soil profile of biogeochemical pools and processes is
explicitly represented, with exponentially increasing layer thickness
with increasing soil depth. For each of these layers, the model
represents different organic pools (metabolic, structural, and woody
litter, as well as fast and slow overturning soil organic matter (SOM)),
as well as inorganic pools of N and P (soluble inorganic , , as well as
, , and , soluble inorganic , adsorbed , absorbed , occluded and primary
). For each soil layer, temperature and moisture are calculated based on
soil physical characteristics, and the transport and atmospheric
exchange of energy and water (Sect.
<a href="#sec:sse:surface_hydrology" data-reference-type="ref"
data-reference="sec:sse:surface_hydrology">[sec:sse:surface_hydrology]</a>).

The turnover and formation of SOM follows, with some modifications, the
CENTURY approach (Parton et al. 1993) (Sect.
<a href="#sec:sb:turnover" data-reference-type="ref"
data-reference="sec:sb:turnover">[sec:sb:turnover]</a> and Sect.
<a href="#sec:sb:formation" data-reference-type="ref"
data-reference="sec:sb:formation">[sec:sb:formation]</a>): The turnover
of litter and SOM pools are all calculated using first-order kinetics
with temperature and moisture dependence. The litter stoichiometry is
determined by the stoichiometry of senescent plant tissue and
tissue-specific allocation fractions. The stoichiometry of the fast SOM
pool adjusts to available inorganic nutrients, whereas the slow SOM pool
is assumed to have fixed C:N:P stoichiometry.

The soluble inorganic , and are assumed to be available for plant and
microbial uptake, as concluded by White (2006) and many others. Plants
and soil microbes compete for these nutrients based on their respective
demand and uptake capacity, which varies in time and with soil depth. In
the aerobic part of soil, is oxidised to through nitrification and in
anaerobic part of soil is reduced to through denitrification (Zaehle et
al. 2011). Both processes produce and as by-products (Sect.
<a href="#sec:sb:nitrification_denitrification"
data-reference-type="ref"
data-reference="sec:sb:nitrification_denitrification">[sec:sb:nitrification_denitrification]</a>).
The representation of soil inorganic P dynamics (Sect.
<a href="#sec:sb:p_flux" data-reference-type="ref"
data-reference="sec:sb:p_flux">[sec:sb:p_flux]</a>) follow the concept
from (Yang et al. 2014) and (Wang, Law, and Pak 2010). The available
soluable is exchanged via ad-/desorption with the mineral surface,
leached with soil water, and replenished by atmospheric deposition,
weathering, and biomineralisation of P that is bound in the slow SOM
pool. The adsorbed P is further absorbed into the soil matrix and
gradually becomes an inactive form of P, occluded P, as proposed by
(Walker and Syers 1976).

### Literature

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="1">

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

<div id="ref-Bernacchi:2001wv" class="csl-entry">

Bernacchi, C J, E L Singsaas, C Pimentel, A R Protis JR, and S P Long.
2001. “<span class="nocase">Improved temperature response functions for
models of Rubisco-limited photosynthesis</span>.” *Plant, Cell &
Environment* 24 (January): 253–59.

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

<div id="ref-Knauer:2015dh" class="csl-entry">

Knauer, Jürgen, Christiane Werner, and Sönke Zaehle. 2015.
“<span class="nocase">Evaluating stomatal models and their atmospheric
drought response in a land surface scheme: A multibiome
analysis</span>.” *Journal of Geophysical Research: Biogeosciences* 120
(10): 1894–1911.

</div>

<div id="ref-Kull:1998cr" class="csl-entry">

Kull, O, and B Kruijt. 1998. “<span class="nocase">Leaf photosynthetic
light response: a mechanistic model for scaling photosynthesis to leaves
and canopies</span>.” *Functional Ecology* 12 (5): 767–77.

</div>

<div id="ref-Lloyd1994-jn" class="csl-entry">

Lloyd, John, and J A Taylor. 1994. “On the Temperature Dependence of
Soil Respiration.” *Functional Ecology*, 315–23.
<https://www.jstor.org/stable/2389824>.

</div>

<div id="ref-Parton:1993ud" class="csl-entry">

Parton, W J, J M O Scurlock, D S Ojima, T G Gilmanov, R J Scholes, D S
Schimmel, T Kirchner, et al. 1993. “<span class="nocase">Observations
and modelling of biomass and soil organic matter dynamics for the
grassland biome worldwide</span>.” *Global Biogeochemical Cycles* 7 (4):
785–809.

</div>

<div id="ref-Sitch:2003wg" class="csl-entry">

Sitch, Stephen A, Benjamin Smith, I C Prentice, Almut Arneth, Alberte
Bondeau, W Cramer, J O Kaplan, et al. 2003.
“<span class="nocase">Evaluation of ecosystem dynamics, plant geography
and terrestrial carbon cycling in the LPJ Dynamic Global Vegetation
Model</span>.” *Global Change Biology* 9: 161–85.

</div>

<div id="ref-Walker:1976" class="csl-entry">

Walker, T. W., and J. K. Syers. 1976. “The Fate of Phosphorus During
Pedogenesis.” Journal Article. *Geoderma* 15 (1): 1–19.
<https://doi.org/10.1016/0016-7061(76)90066-5>.

</div>

<div id="ref-Wang:2010ba" class="csl-entry">

Wang, Y P, R M Law, and B Pak. 2010. “<span class="nocase">A global
model of carbon, nitrogen and phosphorus cycles for the terrestrial
biosphere</span>.” *Biogeosciences* 7 (7): 2261–82.

</div>

<div id="ref-White:2006" class="csl-entry">

White, Robert E. 2006. *Principles and Practice of Soil Science the Soil
as a Natural Resource*. Book. 4th ed. Vol. 32. Malden: Blackwell
Publishing. <https://doi.org/10.1177/0340035206070163>.

</div>

<div id="ref-Yang:2014" class="csl-entry">

Yang, X., P. E. Thornton, D. M. Ricciuto, and W. M. Post. 2014. “The
Role of Phosphorus Dynamics in Tropical Forests – a Modeling Study Using
CLM-CNP.” *Biogeosciences* 11 (6): 1667–81.
<https://doi.org/10.5194/bg-11-1667-2014>.

</div>

<div id="ref-Zaehle:2011gy" class="csl-entry">

Zaehle, S, Philippe Ciais, Andrew D Friend, and Vincent Prieur. 2011.
“<span class="nocase">Carbon benefits of anthropogenic reactive nitrogen
offset by nitrous oxide emissions</span>.” *Nature Geoscience* 4 (8):
601–5.

</div>

<div id="ref-Zaehle:2010cz" class="csl-entry">

Zaehle, S, and Andrew D Friend. 2010. “<span class="nocase">Carbon and
nitrogen cycle dynamics in the O-CN land surface model: 1. Model
description, site-scale evaluation, and sensitivity to parameter
estimates</span>.” *Global Biogeochemical Cycles* 24 (1): GB 1005.

</div>

</div>
