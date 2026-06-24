---
title: Using `mo_jsb_extvar_reader`
order: 1
---

## Purpose

`mo_jsb_extvar_reader` (`src/shared/mo_jsb_extvar_reader.f90`) abstracts the generic
machinery for reading 2D fields from yearly NetCDF boundary-condition files.  It covers:

- building filenames from a prefix and a calendar year,
- reading monthly or timeless 2D fields with optional range validation,
- managing the 0:13 monthly interpolation buffer across year boundaries (experiment start
  and year transitions), and
- linearly interpolating the buffer to a single 2D field at any model timestep.

---

## File and variable naming conventions

### File naming

All yearly files must follow the pattern:

```
{filename_prefix}{YYYY}.nc
```

where `filename_prefix` is a caller-supplied string (e.g. `'/data/bc_lai_'`) and `YYYY`
is a zero-padded 4-digit year.  Years outside 1000–9999 are rejected with a hard abort.

Use `extvar_filename(prefix, year)` to obtain the filename string:

```fortran
USE mo_jsb_extvar_reader, ONLY: extvar_filename
CHARACTER(len=filename_max) :: fname
fname = extvar_filename('/data/bc_lai_', 2010)  ! → '/data/bc_lai_2010.nc'
```

### Variable naming

Variable names within the file are left entirely to the caller.  A recommended convention
(used by `external_lai`) is `{process_prefix}_{tile_name}`, e.g. `pheno_lai_ta_trbee`.

---

## Monthly-mode context API

This is the preferred API when a process reads monthly-mean fields and needs them
throughout the year with seamless cross-year interpolation.

### Buffer layout

The interpolation buffer has shape `(nproma, nblks, 0:13)` with the following slot assignment:

| Slot | Content |
|---|---|
| `0` | December of year *Y*−1 |
| `1`–`12` | January–December of year *Y* |
| `13` | January of year *Y*+1 |

A separate preload buffer `(nproma, nblks, 1:12)` holds all 12 months of year *Y*+1 so
that the year transition requires only **one file read** (the file for *Y*+2).

### Usage pattern

The context API reduces all year-management logic to three calls that bracket the tile loop:

```fortran
USE mo_jsb_extvar_reader, ONLY: t_extvar_monthly_ctx,     &
  &                             extvar_monthly_open,       &
  &                             extvar_monthly_update,     &
  &                             extvar_monthly_close

TYPE(t_extvar_monthly_ctx) :: ctx

! --- called once per year (or at experiment start) ---
CALL extvar_monthly_open(ctx, current_datetime, filename_prefix, model%grid_id)

current_tile => tile%Get_first_child_tile()
DO WHILE (ASSOCIATED(current_tile))
  ! ... obtain memory handle for current_tile ...
  CALL extvar_monthly_update(ctx, varname,          &
    &   mem%my_monthly_buf,                         &  ! (nproma, nblks, 0:13)
    &   mem%my_preload_buf,                         &  ! (nproma, nblks, 1:12)
    &   val_min=0._wp, val_max=upper_bound)
  ! !$ACC UPDATE DEVICE is issued inside extvar_monthly_update ASYNC(1)
  current_tile => current_tile%Get_next_sibling_tile()
END DO
!$ACC WAIT(1)  ! wait for all per-tile device updates

CALL extvar_monthly_close(ctx)
```

`extvar_monthly_open` internally calls `is_time_experiment_start`: at experiment start it
opens three files (years *Y*−1, *Y*, *Y*+1); at each subsequent year transition it opens
only one file (year *Y*+1).  `extvar_monthly_update` dispatches accordingly — the caller
has no branching on experiment start.

### What `extvar_monthly_open` / `extvar_monthly_update` do

**At experiment start** (`extvar_monthly_buffer_init` is called internally):

| Action | Buffer slot(s) |
|---|---|
| Read month 12 from file *Y*−1 | `monthly_buf(:,:,0)` |
| Read all 12 months from file *Y* | `monthly_buf(:,:,1:12)` |
| Read all 12 months from file *Y*+1 | `preload_buf(:,:,:)` |
| Copy Jan of *Y*+1 from preload | `monthly_buf(:,:,13)` |

**At year transition** (`extvar_monthly_buffer_advance` is called internally):

| Action | Buffer slot(s) |
|---|---|
| Copy old December from `monthly_buf(:,:,12)` | `monthly_buf(:,:,0)` |
| Copy all 12 months from `preload_buf` | `monthly_buf(:,:,1:12)` |
| Read all 12 months from file *Y*+2 | `preload_buf(:,:,:)` |
| Copy Jan of *Y*+2 from preload | `monthly_buf(:,:,13)` |

Only one file is read at each year transition; all other data moves are free in-memory copies.

---

## Time interpolation

At every model timestep, call `calculate_time_interpolation_weights` once and then
`extvar_interpolate_monthly` once per tile (or domain):

```fortran
USE mo_jsb_extvar_reader, ONLY: t_time_interpolation_weights,          &
  &                             calculate_time_interpolation_weights,   &
  &                             extvar_interpolate_monthly

TYPE(t_time_interpolation_weights) :: tiw

tiw = calculate_time_interpolation_weights(options%current_datetime)

CALL extvar_interpolate_monthly(                  &
  &   mem%my_monthly_buf,                         &  ! (nproma, nblks, 0:13)
  &   tiw%weight1, tiw%month1_index,              &
  &   tiw%weight2, tiw%month2_index,              &
  &   my_output_field,                            &  ! (nproma, nblks)
  &   lacc=lzacc, acc_stream=acc_stream)
```

`t_time_interpolation_weights` fields:

| Field | Description |
|---|---|
| `weight1`, `weight2` | Linear interpolation weights; `weight1 + weight2 = 1` |
| `month1_index`, `month2_index` | Buffer slot indices in MPIM 0–13 range |

`extvar_interpolate_monthly` is GPU-capable: when `lacc=.TRUE.`, `monthly_buffer` and
`data_out` must already be present on the device (e.g. via an enclosing `!$ACC DATA`
region or module-level `ENTER DATA`).  Scalar arguments are copied to the device automatically.

Both `t_time_interpolation_weights` and `calculate_time_interpolation_weights` are
re-exported from `mo_jsb_extvar_reader` and from `mo_jsb_time`; prefer `mo_jsb_extvar_reader`
in modules that already USE it.

---

## File-reading primitives

For use cases where the context API does not fit (e.g. annual mode, or reading from
externally managed file handles), the low-level primitives are also public.

All read routines operate on a `t_input_file` that the caller opens with
`jsb_netcdf_open_input` and closes with `%Close()`:

```fortran
USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
USE mo_jsb_extvar_reader, ONLY: extvar_filename, read_extvar_2d

TYPE(t_input_file) :: f
f = jsb_netcdf_open_input(extvar_filename(prefix, year), grid_id)
CALL read_extvar_2d(f, varname, my_array, val_min=0._wp)
CALL f%Close()
```

### `Read_2d_time` (method on `t_input_file`)

```fortran
ptr => input_file%Read_2d_time(variable_name, [fill_array(:,:,:)], &
  &                            [start_time_step], [end_time_step], &
  &                            [val_min], [val_max])
```

Reads a contiguous range of monthly time slices; returns a `(nproma, nblks, n)` pointer.
Pass `fill_array` to write directly into an existing array without an extra allocation — the
returned pointer then aliases `fill_array` and must not be deallocated.
Without `fill_array` the caller owns the returned allocation and must `DEALLOCATE(ptr)`.

### `Read_2d_one_time` (method on `t_input_file`)

```fortran
CALL input_file%Read_2d_one_time(variable_name, data_out(:,:), month, &
  &                               [val_min], [val_max])
```

Reads a single monthly time slice directly into a 2D `(nproma, nblks)` array.
Equivalent to calling `Read_2d_time` with `start=end=month` and extracting `ptr(:,:,1)`.

### `Read_2d` (method on `t_input_file`)

```fortran
ptr => input_file%Read_2d(variable_name, [fill_array(:,:)], [val_min], [val_max])
```

Reads a variable with no time dimension; returns a `(nproma, nblks)` pointer.
Pass `fill_array` to write directly into an existing array — the returned pointer then
aliases `fill_array` and must not be deallocated (`NULLIFY(ptr)` after use).
`fill_array` must have the `TARGET` attribute at the call site.

### Range validation

All read routines accept optional `val_min` and `val_max` arguments.  If any read
value falls outside the specified range the model aborts with an error message naming the
variable and file.

---

## Lower-level buffer management

For callers that need finer control than the context API provides,
`extvar_monthly_buffer_init` and `extvar_monthly_buffer_advance` are also public.  Both
require the caller to open and close files explicitly:

```fortran
USE mo_jsb_extvar_reader, ONLY: extvar_monthly_buffer_init,    &
  &                             extvar_monthly_buffer_advance,  &
  &                             extvar_filename
USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input

! --- experiment start ---
TYPE(t_input_file) :: file_prev, file_curr, file_next
file_prev = jsb_netcdf_open_input(extvar_filename(prefix, year - 1), grid_id)
file_curr = jsb_netcdf_open_input(extvar_filename(prefix, year    ), grid_id)
file_next = jsb_netcdf_open_input(extvar_filename(prefix, year + 1), grid_id)
CALL extvar_monthly_buffer_init(file_prev, file_curr, file_next, varname, &
  &                             monthly_buf, preload_buf, val_min=0._wp)
CALL file_prev%Close();  CALL file_curr%Close();  CALL file_next%Close()

! --- year transition ---
TYPE(t_input_file) :: file_next2
file_next2 = jsb_netcdf_open_input(extvar_filename(prefix, year + 2), grid_id)
CALL extvar_monthly_buffer_advance(file_next2, varname, &
  &                                monthly_buf, preload_buf, val_min=0._wp)
CALL file_next2%Close()
```

The context API (`extvar_monthly_open` / `extvar_monthly_update` / `extvar_monthly_close`)
handles all of the above automatically and is preferred for new implementations.

---

## Annual-mode recipe

For annual-mode external BC (single target value per year), file management is simpler —
one file per year is read without a preload buffer.  The standard pattern:

```fortran
USE mo_jsb_extvar_reader, ONLY: extvar_filename
USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
USE mo_jsb_time,          ONLY: get_year

TYPE(t_input_file)          :: input_file
REAL(wp), POINTER           :: ptr2d(:,:)
CHARACTER(len=filename_max) :: filename
INTEGER                     :: current_year

current_year = get_year(current_datetime)
! ALCC convention: read year Y+1 at the start of year Y
filename = extvar_filename(filename_prefix, current_year + 1)
input_file = jsb_netcdf_open_input(TRIM(filename), grid_id)
ptr2d => input_file%Read_2d(variable_name=varname,     &
  &   fill_array=my_target_array,                      &
  &   val_min=0._wp, val_max=upper_bound)
NULLIFY(ptr2d)
CALL input_file%Close()
```

See `src/phenology/mo_pheno_init.f90` (`read_extlai_data`, annual branch) for the
complete reference example including the per-second increment computation.

---

## Fortran note: lower-bound preservation

The monthly buffer argument `monthly_buf(:,:,0:)` uses a declared lower bound of `0` on
the third dimension.  For the lower bound to be preserved by the compiler, the actual
argument must be passed **without subscripts** — not as a section like `lai_mon(:,:,:)`:

```fortran
! correct — lower bound 0 is preserved
CALL extvar_monthly_update(ctx, varname, mem%lai_mon, mem%lai_mon_next, ...)

! wrong — section resets lower bound to 1
CALL extvar_monthly_update(ctx, varname, mem%lai_mon(:,:,:), mem%lai_mon_next, ...)
```

The memory field `lai_mon` must be declared with `DIMENSION(nproma, nblks, 0:13)` (or
equivalent) so that the compiler knows the lower bound when passing it.
