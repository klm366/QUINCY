---
title: Climatology phenology scheme
order: 2
---

## Overview

The `climatology` scheme prescribes LAI (and foliage projected cover) from a fixed monthly
climatology read from a boundary-condition file. The climatology is applied identically
every year: it does not change from year to year, and no prognostic state is evolved. LAI
is linearly interpolated between the two nearest monthly-mean values at every timestep.

Activate it with:
```fortran
&jsb_pheno_nml
  scheme = 'climatology'
/
```

---

## Input file

The monthly-mean LAI and foliage projected cover are read from `bc_filename`
(default `bc_land_phenology.nc`) during model initialisation. The file must contain:

| Variable | Dimensions | Units | Description |
|---|---|---|---|
| `lai` | `(nproma, nblks, 12)` | m² m⁻² | Monthly-mean LAI |
| `fract_fpc` | `(nproma, nblks, 12)` | — | Monthly-mean foliage projected cover |

The 12 monthly slices are stored internally in slots 1–12 of a 0:13 buffer. Slots 0
(= December of previous year) and 13 (= January of next year) are copies of the adjacent
months and enable seamless linear interpolation across the December–January boundary.

---

## Timestep interpolation

At every model timestep the LAI is computed as:

$$\text{LAI}(t) = w_1 \cdot \text{LAI}_{m_1} + w_2 \cdot \text{LAI}_{m_2}$$

where $m_1$ and $m_2$ are the two bracketing monthly indices and $w_1$, $w_2$ are the
linear interpolation weights returned by `get_time_interpolation_weights`. The same
interpolation is applied to `fract_fpc` for general vegetation tiles that have no
associated PFTs.

---

## Dependencies and fallback behaviour

The `climatology` scheme has no runtime dependencies on other processes; it only needs the
boundary-condition file to have been read during initialisation.

When a vegetation tile has no associated PFTs (i.e. `tile%lcts(1)%lib_id == 0`) and any
other scheme is selected, the model falls back silently to `climatology` for that tile.
A message is written to the log in this case.

---

## Namelist settings

```fortran
&jsb_pheno_nml
  scheme      = 'climatology'
  bc_filename = 'bc_land_phenology.nc'  ! Required: monthly LAI and fract_fpc
  ic_filename = 'ic_land_phenology.nc'  ! Read but only affects prognostic state vars
/
```

| Key | Default | Description |
|---|---|---|
| `scheme` | `'logrop'` | Must be set to `'climatology'` |
| `bc_filename` | `'bc_land_phenology.nc'` | File containing monthly LAI and fract_fpc |
| `ic_filename` | `'ic_land_phenology.nc'` | Initial-condition file (heat sums unused but still read) |
