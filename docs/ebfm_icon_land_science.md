# The EBFM–ICON-Land coupling: a science description

This document describes the physics of the coupling between the energy balance and firn model
(EBFM) and ICON-Land (JSBACH) in the TerraDT coupled system — what each model computes, what is
exchanged and why, and how the coupled system differs from EBFM run standalone. The technical
counterpart (YAC configuration, namelists, file formats, run scripts) is
[ebfm_icon_land_coupling.md](https://gitlab.dkrz.de/jsbach/jsbach/-/blob/feature/ebfm-coupling/docs/ebfm_icon_land_coupling.md).

## The models and their roles

The TerraDT system couples the ICON atmosphere (ICON-A, with the land model ICON-Land/JSBACH as
part of its physics), the ICON ocean (ICON-O), the ice-sheet model Elmer/Ice, and EBFM. EBFM — the
energy balance and firn model of Ward van Pelt and co-workers (van Pelt et al., 2012, 2019) —
simulates the surface energy balance, the snow/firn column and the climatic mass balance of
glaciated terrain, here on the (static) Elmer mesh covering Greenland.

In the coupled system the division of labour is:

* **ICON-Land computes the surface energy balance (SEB)** of the EBFM-covered ice — every
  atmospheric time step (40 s at R2B4), on a dedicated land-ice tile (`lice`) — because the
  surface fluxes and radiative properties seen by the atmosphere must be consistent with the
  atmospheric state at every step.
* **EBFM keeps the firn model**: heat conduction, percolation, refreezing, water storage,
  densification, runoff, and the climatic mass balance.
* Conceptually, EBFM acts as the "firn engine" of ICON-Land's ice tile: it provides the surface
  and subsurface *state*, ICON-Land computes the *fluxes* from it. This is the same pattern as the
  ICON atmosphere–ocean coupling, where the ocean provides the SST and the atmosphere computes the
  turbulent and radiative fluxes.

## How EBFM computes the surface energy balance (standalone)

Run standalone (or before this coupling), EBFM solves at every time step (3 h) and mesh cell for
the surface temperature $T_s$ that closes the surface energy balance

$$ E(T_s) \;=\; SW_{in} - SW_{out} + LW_{in} - LW_{out}(T_s) + SHF(T_s) + LHF(T_s) + GHF(T_s) \;=\; 0 ,$$

found by bisection over a ±40 K interval. The individual terms are parameterised from
near-surface meteorology (air temperature, humidity, pressure, wind, cloud cover, precipitation):

* **Shortwave**: top-of-atmosphere insolation with solar geometry and terrain shading, attenuated
  by parameterised atmospheric transmissivities (Rayleigh scattering, water vapour, aerosols;
  after Iqbal, 1983; Meyers & Dale, 1983; Atwater & Brown, 1974) and a cloud-dependent split into
  direct and diffuse radiation (Oerlemans, 1992).
* **Albedo**: a prognostic snow albedo that decays from a fresh-snow value towards a firn value
  with an e-folding time depending on surface temperature (wet vs. cold snow), reset by snowfall;
  bare ice has a constant albedo (Oerlemans & Knap, 1998; Bougamont et al., 2005).
* **Longwave**: incoming longwave from an atmospheric emissivity combining a clear-sky part
  (function of vapour pressure and temperature; Konzelmann et al., 1994) and a cloud part;
  outgoing longwave as black-body emission at $T_s$.
* **Turbulent fluxes**: bulk formulations for sensible and latent heat with an exchange
  coefficient composed of a constant background part and a katabatic contribution that grows with
  the surface-air temperature deficit — a formulation suited to the stably stratified katabatic
  boundary layers of ice sheets. The latent heat flux uses the saturation humidity at $T_s$ over
  ice or water and yields sublimation/deposition (below melting) or evaporation/condensation.
* **Ground (subsurface) heat flux**: heat conduction between the surface and the uppermost firn
  layer, $GHF = k\,(T_1 - T_s)$, with the conductance $k$ derived from the firn density–dependent
  thermal conductivity over the depth of the surface layer.

If the solved $T_s$ exceeds the melting point, $T_s$ is capped at 273.15 K and the residual energy
becomes **melt energy**: $M = E(T_0)\,\Delta t / (\rho_w L_m)$ [m w.e.]. The latent heat flux is
converted into the corresponding **moisture fluxes** (sublimation, deposition, evaporation,
condensation) that enter the mass balance.

## The firn model and runoff

Below the surface, EBFM maintains a multi-layer firn column at every mesh cell (also on ice-free
cells, where a standard column is assumed):

* **Heat conduction** through the layered column with density-dependent conductivity and heat
  capacity.
* **Melt water and rain percolate** downward (selectable schemes: instantaneous bucket
  percolation, or deep percolation with a normal/linear distribution of infiltration depths,
  representing preferential flow).
* Water **refreezes** where the column is below melting point (limited by cold content and pore
  space), releasing latent heat; remaining water is held as **irreducible water** against gravity
  up to a density-dependent capacity, and beyond that forms **slush** on impermeable horizons.
* **Runoff** is generated from water that can neither refreeze nor be stored: a surface component
  (melt/rain on bare ice), a slush component released with a linear-reservoir time scale, and the
  slow release of irreducible and deep water. Runoff at a given moment therefore reflects the
  firn's storage history, not the instantaneous melt.
* **Densification** by destructive metamorphism, overburden pressure and wind packing gradually
  turns snow into firn and ice; layers are merged/split accordingly.
* The **climatic mass balance** (smb) is the sum of snowfall, deposition, rain that does not run
  off, minus sublimation and runoff — the quantity handed to the ice-dynamics side.

The essential scientific point for the coupling: between the atmosphere and the runoff stands a
**reservoir with memory** (the firn), which retains melt water and precipitation, refreezes part
of it, and releases the rest with delays from hours (slush) to seasons (irreducible/deep water).

## What changes in the coupled system

### ICON-Land computes the SEB with EBFM's firn underneath

On the `lice` tile, ICON-Land's implicit surface energy balance (TMX turbulence, semi-explicit
coupling to the atmosphere) runs every 40 s as on any other tile, but three ingredients come from
EBFM's state, received at every 3 h coupling instant:

* the **ground heat flux** is computed against the first firn layer:
  $GHF = k\,(T_1 - T_s^{old})$, with conductance $k$ (`ghf_cond`) and firn temperature $T_1$
  (`t_sub`) from EBFM — replacing JSBACH's own soil/snow column, whose energy scheme (SSE) is
  switched off on this tile;
* the **surface-layer heat capacity** in the implicit solver is EBFM's firn surface-layer heat
  capacity (`hcap_sub`), so the thermal inertia of the surface matches the firn it sits on;
* the **albedo** (both radiation bands) is EBFM's prognostic snow/ice albedo, so the radiation
  balance carries EBFM's snow-aging memory.

The surface temperature of the lice tile is also initialised from the firn temperature at
experiment start, so the coupled system starts on the firn's thermal state rather than a generic
initial value.

### EBFM's firn is driven by ICON-Land's surface balance

In return, EBFM receives at each of its 3 h steps the time-averaged results of ICON-Land's SEB:
surface temperature, snow/ice melt, and total evapotranspiration (which EBFM partitions into
sublimation/deposition/evaporation/condensation by surface temperature, as it does with its own
latent heat flux). These replace EBFM's own SEB results as input to the firn model — EBFM's SEB is
still evaluated every step, but purely as a **diagnostic** (`ebm_*` in its output), which makes
the two energy balances directly comparable cell by cell.

Precipitation and the meteorological forcing continue to come from ICON-A directly (they are
inputs to both models' surface schemes and to EBFM's mass balance).

### Runoff takes the physical route

Uncoupled, JSBACH has no firn: on its glacier tiles all water (P − E) is passed to the river
routing (HD) immediately, which conserves water trivially but ignores storage. Coupled, the runoff
handed to HD on the lice tile is **EBFM's firn runoff** — melt and rain filtered through
percolation, refreezing and storage. The difference between the instantaneous ICON-Land water flux
(P + E − runoff) and EBFM's runoff is booked as an explicit firn water storage (`weq_firn`) in
ICON-Land's water balance, which therefore stays closed while representing the physically delayed
freshwater release. Accumulated snow no longer reaches the ocean immediately: it stays in the firn
until EBFM releases it — or until it leaves as ice through the dynamics (next section).

## The EBFM–Elmer/Ice side

EBFM also closes the loop with the ice dynamics. To Elmer/Ice it sends the **climatic mass
balance** (the surface forcing of the ice sheet), the **runoff**, and the **near-surface ice
temperature** (a boundary condition for the thermomechanics). From Elmer it receives the evolving
**surface elevation** of the ice sheet, which feeds back on EBFM's surface climate. Elmer/Ice, in
turn, exchanges with the ICON ocean, where dynamic mass loss enters the ocean's freshwater and
heat budgets.

Together this closes the ice-sheet mass loop physically: accumulation and melt are computed
against an atmosphere-consistent energy balance (ICON-Land), buffered and partitioned by the firn
(EBFM), the retained part feeds the ice dynamics as smb (Elmer), and mass leaves either as liquid
runoff through the land hydrology (HD) or as ice discharge through the ocean.

## Known issues and open science questions

* **Turbulent flux differences**: ICON-Land's (TMX) sublimation over the ice sheet is a factor
  3–6 larger than EBFM's bulk/katabatic formulation. Attribution (Jan 1-day tests, 2026-08-24):
  the dominant factor is the roughness treatment — TMX used z0h = z0m = 5 mm (bare-soil value;
  the properly computed scalar roughness of 1e-4 m was diagnostic only — a single roughness is a
  property of ICON's whole Smagorinsky surface-transfer lineage, LES → vdiff-Smagorinsky → TMX,
  inherited from prescribed-flux LES benchmarks, while ICON's other land schemes — TTE vdiff via
  z0h_lnd, turbtran via the laminar resistance rlam_heat — keep the momentum/scalar distinction,
  and the stability functions' own source, Zeng et al. 1998, prescribes a distinct z0t); at 10 m
  the neutral C_H is overestimated by 1.5× from z0h alone, ~1.8× incl. an ice-sheet z0m of
  1e-3 m (measured: 1.6× in the mean exchange coefficient). Fixing z0m/z0h
  (`roughness_lice`, `l_rough_h_lice`) halves the excess, and substituting EBFM's katabatic
  closure for the scalar exchange (`l_katabatic_lice`) brings the ratio to ~1.7× (the remainder
  is the 40 s instantaneous vs. 3 h mean evaluation and EBFM's dynamically fitted lapse rate).
  All three are namelist switches, default off — which combination becomes the default for
  production is an open decision, as is the general question which formulation is more
  appropriate over Greenland's katabatic boundary layer.
* **Static ice geometry in ICON-Land**: the lice tile fraction is fixed at initialisation (from
  Elmer's ice mask); EBFM's own ice mask is likewise frozen at its initialisation. Retreating or
  advancing ice margins are not yet reflected in the tile fractions (this requires ICON-Land's
  land-cover-change machinery with water/energy conservation).
* **Firn spin-up**: EBFM's initial firn column drains for a while (elevated runoff and negative
  smb in the first days) and its deep-water reservoir carries a climatological memory; multi-year
  spin-up is needed before mass-balance numbers are meaningful.
* **Forcing representativeness**: EBFM receives the ICON forcing by nearest-neighbour
  interpolation at the ICON grid's resolution; how the elevation difference between the coarse
  ICON orography and the (much finer) Elmer mesh topography is treated in the coupled forcing —
  standalone EBFM applies lapse-rate corrections — deserves a closer look, especially for the
  steep margins.
* **One energy balance, two melt paths**: melt is computed by ICON-Land, but refreezing happens in
  EBFM; the latent heat of refreezing warms the firn (EBFM) but is not returned to the
  atmosphere's energy budget within the same step. The 3 h exchange makes this a small,
  time-lagged inconsistency; it should be kept in mind when analysing the coupled energy budget.
* **Runoff is delivered to two models**: EBFM sends its runoff both to ICON-Land (where it drives
  the river routing, HD — the active freshwater path to the ocean) and to Elmer/Ice, where it is
  currently only received into an exported variable and used by no solver, so there is no double
  counting today (EBFM's `smb` also already has runoff subtracted, so the ice-mass budget is
  clean). This changes the moment Elmer starts consuming the field — e.g. as meltwater input to a
  GlaDS subglacial hydrology (fine as a sliding forcing, but its marginal discharge must then not
  also reach the ocean) or as direct subglacial/frontal discharge into the ocean coupling
  (immediate double counting with HD). At that point a decision is needed on **which path carries
  the mass to the ocean** — with a scientific case for the Elmer path, since HD delivers at the
  river mouths of a coarse routing network while discharge through Elmer would enter at the
  actual ice margin/fjord locations that matter for fjord circulation.

## References

* van Pelt, W. J. J., Oerlemans, J., Reijmer, C. H., Pohjola, V. A., Pettersson, R., and van
  Angelen, J. H. (2012): Simulating melt, runoff and refreezing on Nordenskiöldbreen, Svalbard,
  using a coupled snow and energy balance model, The Cryosphere, 6, 641–659.
* van Pelt, W. J. J., et al. (2019): A long-term dataset of climatic mass balance, snow
  conditions, and runoff in Svalbard (1957–2018), The Cryosphere, 13, 2259–2280.
* Oerlemans, J. and Knap, W. H. (1998): A 1 year record of global radiation and albedo in the
  ablation zone of Morteratschgletscher, Switzerland, J. Glaciol., 44, 231–238.
* Bougamont, M., Bamber, J. L., and Greuell, W. (2005): A surface mass balance model for the
  Greenland Ice Sheet, J. Geophys. Res., 110, F04018.
* Konzelmann, T., et al. (1994): Parameterization of global and longwave incoming radiation for
  the Greenland Ice Sheet, Global Planet. Change, 9, 143–164.
