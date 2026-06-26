# CLAUDE.md — branch `feature/support-external-lai`

## What this branch does

Adds a third phenology scheme `external_lai` that prescribes time-varying, per-PFT leaf area
index (LAI) from yearly NetCDF files, bypassing the prognostic `logrop` and static
`climatology` schemes.

Activate with `scheme = 'external_lai'` in `&jsb_pheno_nml`.

## Key files changed

| File | Role |
|---|---|
| `src/phenology/mo_pheno_config_class.f90` | New namelist parameters: `extlai_mode`, `extlai_filename_prefix`, `extlai_lai_max`; validation against QUINCY / forest-age-class usecases |
| `src/phenology/mo_pheno_memory_class.f90` | No new fields; `external_lai` reuses `lai_mon` from the climatology scheme |
| `src/phenology/mo_pheno_init.f90` | `read_extlai_data()`: reads yearly NetCDF files via `mo_jsb_extvar_reader` context API |
| `src/phenology/mo_pheno_interface.f90` | `update_phenology_external_lai()`: per-timestep LAI interpolation (monthly mode) |
| `src/shared/mo_jsb_extvar_reader.f90` | **NEW** — general-purpose reader and interpolator for external yearly BC variables; see API below |
| `src/base/mo_jsb_time.f90` | Added `t_time_interpolation_weights` and `calculate_time_interpolation_weights` to PUBLIC (re-exported from `mo_jsb_time_iface`) |
| `doc/guide/3-user-guide/02_phenology/external_lai.md` | User-facing documentation |

## Temporal modes

### Monthly mode (`extlai_mode = 'monthly'`)

LAI is linearly interpolated between bracketing monthly means using
`get_time_interpolation_weights` (same mechanism as `climatology`).

The interpolation buffer `lai_mon(0:13)` stores:
- Slot 0  : December of year Y−1
- Slots 1–12 : January–December of year Y
- Slot 13 : January of year Y+1

**File reading strategy (2 files/year after init, no preload buffer):**

At experiment start, three files are read (Y−1, Y, Y+1):

```
lai_mon(:,:,0)    = file(Y-1), month 12     ! Dec of Y-1
lai_mon(:,:,1:12) = file(Y),   months 1-12  ! full year Y
lai_mon(:,:,13)   = file(Y+1), month 1      ! Jan of Y+1
```

At each year transition Y → Y+1, two files are read (Y+1 and Y+2):

```
lai_mon(:,:,0)    = lai_mon(:,:,12)          ! Dec of Y, free
lai_mon(:,:,1:12) = file(Y+1), months 1-12  ! full year Y+1
lai_mon(:,:,13)   = file(Y+2), month 1      ! Jan of Y+2
```

`lai_mon` is read-only on the GPU; the CPU copy is always authoritative.
No `ACC UPDATE HOST` is needed before the year-transition copy.

## Input file conventions

- Filename pattern: `{extlai_filename_prefix}{YYYY}.nc`
- Variable names: `pheno_lai_{tile_name}` for each PFT child tile of `veg`
- Monthly mode dimensions: `(nproma, nblks, 12)`

## General-purpose extvar reader (`src/shared/mo_jsb_extvar_reader.f90`)

Abstracts file I/O and buffer management for any variable read from yearly NetCDF BC files
(`{prefix}{YYYY}.nc`).  Callers supply the tile loop and memory access; all temporal logic
(experiment start vs. year transition, slot shuffling, preloading) is hidden inside.

### File-reading primitives (caller opens/closes file)

- `extvar_filename(prefix, year)` → filename string; aborts on non-4-digit year
- `t_input_file%Read_2d_time(var, [fill_array], [start/end_time_step], [val_min/max])` — monthly slice range; `fill_array` writes directly into an existing array (no dealloc needed); without it the caller must `DEALLOCATE` the returned pointer
- `t_input_file%Read_2d_one_time(var, time, [fill_array], [val_min/max])` — single time slice, returns `POINTER(:,:)`; if `fill_array` is given the pointer aliases it (no `DEALLOCATE` needed), otherwise caller must `DEALLOCATE`
- `read_extvar_2d(file, var, out(:,:), [val_min/max])` — timeless 2D field

### Monthly-buffer context API (preferred for monthly mode)

Opaque type `t_extvar_monthly_ctx` holds the open files and an `is_init` flag.

```fortran
CALL extvar_monthly_open(ctx, current_datetime, filename_prefix, grid_id)
DO WHILE (ASSOCIATED(current_tile))
  CALL extvar_monthly_update(ctx, varname, tile_lai_mon, val_min=, val_max=)
  ! !$ACC UPDATE DEVICE(monthly_buf) ASYNC(1) is issued inside extvar_monthly_update
  ...
END DO
!$ACC WAIT(1)
CALL extvar_monthly_close(ctx)
```

`extvar_monthly_open` calls `is_time_experiment_start` internally: at experiment start it opens
three files (Y-1, Y, Y+1); at every subsequent year transition it opens two files (Y+1, Y+2).
`extvar_monthly_update` dispatches to `extvar_monthly_buffer_init` or `extvar_monthly_buffer_advance`
accordingly.  The lower-level `buffer_init` / `buffer_advance` routines are also public for
callers that need finer control.

### Time-interpolation API

- `extvar_interpolate_monthly(buf, out [,acc_stream])` — generic; fetches cached weights via
  `get_time_interpolation_weights` internally, no weight arguments at the call site:
  - full-field: `buf(:,:,0:)`, `out(:,:)` (`INTENT(OUT)`)
  - on-chunk:   `buf(:,0:)`,   `out(:)`   (`INTENT(INOUT)`, writes in-place)

## Constraints

- Not supported with `model_scheme = 'quincy'`
- Not supported with `jsbach_forest_age_classes` usecase
- `fract_fpc` is derived from LAI via `get_foliage_projected_cover` (same as `logrop`)
