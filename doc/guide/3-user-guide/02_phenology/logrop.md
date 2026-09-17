---
title: Logrop phenology scheme
order: 1
---

## Overview

`logrop` is the default prognostic phenology scheme. It calculates LAI dynamically at
every model timestep from environmental drivers (temperature, NPP, soil moisture) using
plant-functional-type-specific growth and senescence parameterisations.

Activate it with:
```fortran
&jsb_pheno_nml
  scheme = 'logrop'
/
```

---

## Phenology types

The scheme selects a sub-model for each PFT based on the `PhenologyType` entry in the
land-cover type library (`lctlib`):

| `PhenologyType` | Vegetation type | Key drivers |
|---|---|---|
| 0 | None (bare/glacier) | — |
| 1 | Evergreen | Day of year, pseudo soil temperature, heat sums |
| 2 | Summergreen (deciduous) | Day of year, temperature, chilling and heat sums, latitude |
| 3 | Raingreen | Root-zone soil moisture, NPP |
| 4 | Grasses | Air temperature, soil moisture, NPP |
| 5 | Crops | Temperature, NPP, soil moisture, day length |

Each sub-model computes the growth and senescence phases of its respective vegetation type
by tracking internal state variables (heat sums, chill days, growth phase indicators) that
are propagated across timesteps.

---

## Dependencies

`logrop` requires the assimilation process (`ASSIMI_`) to be active, because it uses
potential net primary productivity (`NPP_pot_rate_ca`) as a driver. The model will abort at
startup if `ASSIMI_` is not enabled alongside `PHENO_`.

Process dependencies:

| Process | Role |
|---|---|
| `ASSIMI_` | Provides `NPP_pot_rate_ca` (daily accumulated and averaged) |
| `SEB_` | Provides daily min/max/mean temperature and pseudo soil temperature |
| `HYDRO_` | Provides relative root-zone soil moisture |

---

## Forest regrowth option

Setting `l_forestRegrowth = .TRUE.` in `jsb_pheno_nml` enables an optional coupling to the
carbon cycle: for forest PFTs the maximum LAI (`LaiMax`) is derived dynamically from
allometric relationships based on standing wood carbon biomass. This requires the `CARBON_`
process to be active.

```fortran
&jsb_pheno_nml
  scheme          = 'logrop'
  l_forestRegrowth = .TRUE.
/
```

When enabled:
- At the start of each year the maximum LAI is recomputed from allometric relationships.
- Carbon pools are rescaled to account for changes in the vegetation fraction correction
  factor.
- Carbon conservation is checked and any residual imbalance is recorded in
  `cconservation_allom`.

---

## Namelist settings

```fortran
&jsb_pheno_nml
  scheme           = 'logrop'                ! Required
  l_forestRegrowth = .FALSE.                 ! Enable biomass-dependent max LAI
  bc_filename      = 'bc_land_phenology.nc'  ! Monthly LAI climatology (used for IC)
  ic_filename      = 'ic_land_phenology.nc'  ! Initial values for heat sums and state vars
/
```

| Key | Default | Description |
|---|---|---|
| `scheme` | `'logrop'` | Must be `'logrop'` for this scheme |
| `l_forestRegrowth` | `.FALSE.` | Enable allometric max-LAI for forest PFTs |
| `bc_filename` | `'bc_land_phenology.nc'` | Source of climatological LAI for IC |
| `ic_filename` | `'ic_land_phenology.nc'` | Initial state (heat sums, growth phases, …) |
