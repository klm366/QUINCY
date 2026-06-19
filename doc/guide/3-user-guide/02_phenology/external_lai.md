---
title: External LAI forcing
order: 3
---

## Overview

The `external_lai` scheme prescribes the leaf area index (LAI) of each plant functional
type (PFT) tile directly from yearly NetCDF input files, bypassing the prognostic `logrop`
and climatological `climatology` schemes. This is useful for offline LAI reconstructions,
sensitivity studies, and atmospheric reanalysis driving.

Activate it with:
```fortran
&jsb_pheno_nml
  scheme = 'external_lai'
/
```

---

## Temporal modes

### Monthly mode (`extlai_mode = 'monthly'`)

Each yearly file contains 12 monthly-mean LAI fields per PFT. At every model timestep the
LAI is linearly interpolated between the two bracketing monthly-mean values, using the same
time-interpolation weights as the built-in `climatology` scheme. This produces a smooth,
continuously varying LAI within each year.

The current-year file is read at experiment start and at the first timestep of each new
year.

### Annual mode (`extlai_mode = 'annual'`)

Each yearly file contains a single annual target LAI per PFT. Following the ALCC
convention, the file for year *Y*+1 is read at the beginning of year *Y* and used as the
end-of-year target. A constant per-second increment is computed from
`(target_lai - current_lai) / year_length_in_seconds` and applied every model timestep, so
that LAI reaches the prescribed target value by 31 December. The increment is stored in
the restart file so that mid-year restarts resume the same linear ramp.

---

## Input file conventions

### File naming

Files must follow the pattern `{extlai_filename_prefix}{YYYY}.nc`, e.g.
`fc_land_lai_2010.nc` with the default prefix `fc_land_lai_`.

| Mode | Year read | Purpose |
|---|---|---|
| `monthly` | *Y* | Current-year monthly LAI values |
| `annual` | *Y* + 1 | End-of-year LAI target for year *Y* |

### Variable names

Each PFT tile must have its own variable in the file, named `pheno_lai_ta_{tile_name}`,
where `{tile_name}` is the name of the corresponding child tile of the `veg` tile in the
model tile tree. For example, if your setup has tiles `trbee`, `trdc`, and `grass`, the
file must contain variables `pheno_lai_ta_trbee`, `pheno_lai_ta_trdc`, and
`pheno_lai_ta_grass`.

### Dimensions and units

| Mode | Expected variable dimensions | Units |
|---|---|---|
| `monthly` | `(nproma, nblks, 12)` | m² m⁻² |
| `annual` | `(nproma, nblks)` | m² m⁻² |

All LAI values must be **non-negative** and must not exceed the configured upper bound
`extlai_lai_max` (default 20 m² m⁻²). The model aborts with a clear error message if
either constraint is violated.

---

## Namelist configuration

```fortran
&jsb_pheno_nml
  scheme                  = 'external_lai'   ! Select this scheme
  extlai_mode             = 'monthly'        ! 'monthly' (default) or 'annual'
  extlai_filename_prefix  = 'fc_land_lai_'   ! Prefix for yearly input files
  extlai_lai_max          = 20.0             ! Upper bound guard [m2 m-2]
/
```

| Key | Type | Default | Description |
|---|---|---|---|
| `scheme` | string | `'logrop'` | Must be set to `'external_lai'` |
| `extlai_mode` | string | `'monthly'` | `'monthly'` or `'annual'` |
| `extlai_filename_prefix` | string | `'fc_land_lai_'` | File path prefix; year and `.nc` are appended |
| `extlai_lai_max` | real | `20.0` | Maximum allowed LAI value; model aborts if exceeded |

---

## Constraints and known limitations

- **QUINCY not supported.** The `external_lai` scheme cannot be used when
  `model_scheme = 'quincy'`. The model aborts during configuration with a clear error
  message.

- **Forest age classes not supported.** The `external_lai` scheme cannot be combined with
  the `jsbach_forest_age_classes` usecase, because the age-class tile structure requires
  prognostic phenology. The model aborts during configuration.

- **Foliage projected cover.** The `fract_fpc` variable is derived from the prescribed LAI
  using `get_foliage_projected_cover` (same as in `logrop`), driven by `fract_fpc_max` and
  the clumpiness factor of each PFT.

- **No post-logrop override.** This scheme is an alternative to `logrop`, not a
  post-processing override of it. To prescribe LAI while keeping prognostic carbon
  cycling, a separate mechanism would be needed.
