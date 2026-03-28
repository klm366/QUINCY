---
title: Phenology
---

## Overview

The phenology process (`PHENO_`) controls the leaf area index (LAI) of each plant
functional type (PFT) tile and derives the foliage projected cover (`fract_fpc`) from it.
Three mutually exclusive schemes are available, selected via the `scheme` key in the
`jsb_pheno_nml` namelist:

| Scheme | Description | Default? |
|---|---|---|
| [`logrop`](logrop.md) | Prognostic phenology driven by temperature and NPP | Yes |
| [`climatology`](climatology.md) | Prescribed monthly-mean LAI from a static file | No |
| [`external_lai`](external_lai.md) | Prescribed time-varying per-PFT LAI from yearly files | No |

The scheme can be changed without recompilation — only a namelist edit is required.

---

## Shared namelist settings

All schemes share two namelist keys for boundary-condition and initial-condition files.
The `external_lai`-specific keys are documented on the
[external_lai page](external_lai.md).

```fortran
&jsb_pheno_nml
  scheme       = 'logrop'                ! 'logrop' | 'climatology' | 'external_lai'
  bc_filename  = 'bc_land_phenology.nc' ! Boundary conditions (climatology and logrop)
  ic_filename  = 'ic_land_phenology.nc' ! Initial conditions for prognostic variables
/
```

| Key | Default | Description |
|---|---|---|
| `scheme` | `'logrop'` | Phenology scheme to use |
| `bc_filename` | `'bc_land_phenology.nc'` | Monthly LAI/fpc climatology (used by `climatology`; also read for `logrop` initial LAI) |
| `ic_filename` | `'ic_land_phenology.nc'` | Initial values for heat sums and other prognostic PHENO state variables |

---

## Foliage projected cover

All three schemes feed their LAI output into `get_foliage_projected_cover`, which derives
`fract_fpc` from LAI, the maximum foliage projected cover `fract_fpc_max`, and the
PFT-specific clumpiness factor from the land-cover type library. The `climatology` scheme
additionally uses a prescribed `fract_fpc` monthly climatology from `bc_filename` for
general vegetation tiles that have no associated PFTs.
