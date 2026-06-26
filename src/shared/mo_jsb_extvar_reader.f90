!> General-purpose reader and time-interpolator for external time-varying boundary-condition
!> variables stored in yearly NetCDF files.
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Reading and interpolating 2D boundary-condition fields from yearly NetCDF files
!>
!> Yearly files follow the naming pattern `{filename_prefix}{YYYY}.nc`.
!>
!> **File-reading API** -- all read routines operate on a `t_input_file` that the caller has already
!> opened with `jsb_netcdf_open_input`.  Optional `val_min` / `val_max` arguments trigger a hard
!> abort when any read value falls outside the specified range.
!>
!> - `extvar_filename`       -- build the filename for a given year (aborts on non-4-digit year)
!>
!> **Time-interpolation API** -- fetches the cached interpolation weights via
!> `get_time_interpolation_weights` and applies mid-month linear interpolation.
!>
!> - `extvar_interpolate_monthly(buf, out [,acc_stream])` -- generic; dispatch is by rank:
!>   - full-field: `buf(:,:,0:)`, `out(:,:)` (INTENT OUT)
!>   - on-chunk:   `buf(:,0:)`,   `out(:)`   (INTENT INOUT, writes in-place)
!>   Both overloads call `get_time_interpolation_weights` internally; no weight arguments needed.
!>
MODULE mo_jsb_extvar_reader
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message, finish
  USE mo_jsb_control,            ONLY: debug_on
  USE mo_jsb_io_netcdf, ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_io_units,      ONLY: filename_max
  USE mo_jsb_time,      ONLY: t_datetime,                        &
    &                         get_time_interpolation_weights,    &
    &                         get_year, is_time_experiment_start

  IMPLICIT NONE
  PRIVATE

  !> Opaque context for one monthly-buffer update cycle (experiment start or year transition).
  !> Create with `extvar_monthly_open`, use in the tile loop with `extvar_monthly_update`,
  !> then release with `extvar_monthly_close`.
  !>
  !> At experiment start: `file_prev` (Y-1), `file_curr` (Y), `file_next` (Y+1) are all open.
  !> At a year transition: only `file_curr` (Y+1) and `file_next` (Y+2) are open.
  TYPE, PUBLIC :: t_extvar_monthly_ctx
    PRIVATE
    TYPE(t_input_file) :: file_prev  !< Open at experiment start only (year Y-1)
    TYPE(t_input_file) :: file_curr  !< Year Y at init; year Y+1 at advance
    TYPE(t_input_file) :: file_next  !< Year Y+1 at init; year Y+2 at advance (Jan only)
    LOGICAL            :: is_init = .FALSE.
  END TYPE t_extvar_monthly_ctx

  ! File-reading utilities
  PUBLIC :: extvar_filename

  ! Monthly-buffer management -- high-level context API
  PUBLIC :: extvar_monthly_open
  PUBLIC :: extvar_monthly_update
  PUBLIC :: extvar_monthly_close
  ! Monthly-buffer management -- lower-level (single-phase) routines
  PUBLIC :: extvar_monthly_buffer_init
  PUBLIC :: extvar_monthly_buffer_advance

  ! Time-interpolation utility - generic: 3D/2D (full-field) or 2D/1D (on-chunk)
  INTERFACE extvar_interpolate_monthly
    MODULE PROCEDURE extvar_interpolate_monthly_2d
    MODULE PROCEDURE extvar_interpolate_monthly_1d
  END INTERFACE extvar_interpolate_monthly
  PUBLIC :: extvar_interpolate_monthly

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_extvar_reader'

CONTAINS

  ! ======================================================================================================= !
  !> Build the path for a yearly boundary-condition file.
  !>
  !> Returns `'{filename_prefix}{YYYY}.nc'`.  Aborts if `year` is outside the 4-digit range 1000-9999.
  !>
  FUNCTION extvar_filename(filename_prefix, year) RESULT(filename)

    CHARACTER(len=*), INTENT(in) :: filename_prefix  !< Common prefix, e.g. '/data/lai_'
    INTEGER,          INTENT(in) :: year              !< Calendar year (must be 4-digit)
    CHARACTER(len=filename_max)  :: filename

    CHARACTER(len=*), PARAMETER :: routine = modname//':extvar_filename'

    IF ((year < 1000) .OR. (year > 9999)) &
      CALL finish(TRIM(routine), 'Yearly BC filenames require a 4-digit year (1000-9999).')

    WRITE(filename, '(a,I4.4,a)') TRIM(filename_prefix), year, '.nc'

  END FUNCTION extvar_filename

  ! ======================================================================================================= !
  !> Open the yearly BC file(s) needed for one monthly-buffer update cycle.
  !>
  !> At experiment start (`is_time_experiment_start`) three files are opened: year Y-1, Y, and Y+1.
  !> At every subsequent year transition two files are opened: year Y+1 (full) and Y+2 (Jan only).
  !> The resulting context is passed unchanged to `extvar_monthly_update` (once per tile) and finally
  !> to `extvar_monthly_close`.
  !>
  SUBROUTINE extvar_monthly_open(ctx, current_datetime, filename_prefix, grid_id)

    TYPE(t_extvar_monthly_ctx), INTENT(out) :: ctx
    TYPE(t_datetime), POINTER,  INTENT(in)  :: current_datetime
    CHARACTER(len=*),           INTENT(in)  :: filename_prefix  !< Prefix for `{prefix}{YYYY}.nc`
    INTEGER,                    INTENT(in)  :: grid_id

    INTEGER                    :: year
    CHARACTER(len=filename_max) :: fname_prev, fname_curr, fname_next
    CHARACTER(len=*), PARAMETER :: routine = modname//':extvar_monthly_open'

    year          = get_year(current_datetime)
    ctx%is_init   = is_time_experiment_start(current_datetime)

    IF (ctx%is_init) THEN
      fname_prev = extvar_filename(filename_prefix, year - 1)
      fname_curr = extvar_filename(filename_prefix, year)
      fname_next = extvar_filename(filename_prefix, year + 1)
      IF (debug_on()) CALL message(TRIM(routine), 'Experiment start: opening ' &
        & //TRIM(fname_prev)//', '//TRIM(fname_curr)//', '//TRIM(fname_next))
      ctx%file_prev = jsb_netcdf_open_input(TRIM(fname_prev), grid_id)
      ctx%file_curr = jsb_netcdf_open_input(TRIM(fname_curr), grid_id)
      ctx%file_next = jsb_netcdf_open_input(TRIM(fname_next), grid_id)
    ELSE
      fname_curr = extvar_filename(filename_prefix, year)
      fname_next = extvar_filename(filename_prefix, year + 1)
      IF (debug_on()) CALL message(TRIM(routine), 'Year transition: opening ' &
        & //TRIM(fname_curr)//', '//TRIM(fname_next))
      ctx%file_curr = jsb_netcdf_open_input(TRIM(fname_curr), grid_id)
      ctx%file_next = jsb_netcdf_open_input(TRIM(fname_next), grid_id)
    END IF

  END SUBROUTINE extvar_monthly_open

  ! ======================================================================================================= !
  !> Update a per-tile monthly buffer using the context created by `extvar_monthly_open`.
  !>
  !> Dispatches to `extvar_monthly_buffer_init` (experiment start) or
  !> `extvar_monthly_buffer_advance` (year transition) depending on the context flag.
  !> Call once per tile inside the tile loop between `extvar_monthly_open` and
  !> `extvar_monthly_close`.
  !>
  SUBROUTINE extvar_monthly_update(ctx, varname, monthly_buf, val_min, val_max)

    TYPE(t_extvar_monthly_ctx), INTENT(inout)        :: ctx
    CHARACTER(len=*),           INTENT(in)           :: varname
    REAL(wp),                   INTENT(inout)        :: monthly_buf(:,:,0:)  !< (nproma, nblks, 0:13)
    REAL(wp),                   INTENT(in), OPTIONAL :: val_min
    REAL(wp),                   INTENT(in), OPTIONAL :: val_max

    IF (ctx%is_init) THEN
      CALL extvar_monthly_buffer_init(ctx%file_prev, ctx%file_curr, ctx%file_next, varname, &
        &                             monthly_buf, val_min=val_min, val_max=val_max)
    ELSE
      CALL extvar_monthly_buffer_advance(ctx%file_curr, ctx%file_next, varname, &
        &                                monthly_buf, val_min=val_min, val_max=val_max)
    END IF

    !$ACC UPDATE DEVICE(monthly_buf) ASYNC(1)

  END SUBROUTINE extvar_monthly_update

  ! ======================================================================================================= !
  !> Close the file(s) opened by `extvar_monthly_open`.
  !>
  SUBROUTINE extvar_monthly_close(ctx)

    TYPE(t_extvar_monthly_ctx), INTENT(inout) :: ctx

    IF (ctx%is_init) CALL ctx%file_prev%Close()
    CALL ctx%file_curr%Close()
    CALL ctx%file_next%Close()

  END SUBROUTINE extvar_monthly_close

  ! ======================================================================================================= !
  !> Fill a `(nproma, nblks, 0:13)` monthly interpolation buffer at experiment start.
  !>
  !> Three yearly files must be open.  Slots are populated as:
  !>   - `monthly_buf(:,:,0)`    = month 12 of `file_prev`  (December of year Y-1)
  !>   - `monthly_buf(:,:,1:12)` = all months of `file_curr` (year Y)
  !>   - `monthly_buf(:,:,13)`   = month 1 of `file_next`   (January of year Y+1)
  !>
  SUBROUTINE extvar_monthly_buffer_init(file_prev, file_curr, file_next, varname, &
    &                                   monthly_buf, val_min, val_max)

    TYPE(t_input_file), INTENT(inout)        :: file_prev            !< Open file for year Y-1
    TYPE(t_input_file), INTENT(inout)        :: file_curr            !< Open file for year Y
    TYPE(t_input_file), INTENT(inout)        :: file_next            !< Open file for year Y+1
    CHARACTER(len=*),   INTENT(in)           :: varname
    REAL(wp),           INTENT(inout)        :: monthly_buf(:,:,0:)  !< (nproma, nblks, 0:13)
    REAL(wp),           INTENT(in), OPTIONAL :: val_min
    REAL(wp),           INTENT(in), OPTIONAL :: val_max

    REAL(wp), POINTER :: ptr_2d(:,:), ptr_3d(:,:,:)

    ptr_2d => file_prev%Read_2d_one_time(varname, time=12, fill_array=monthly_buf(:,:, 0), &
      & val_min=val_min, val_max=val_max)
    ptr_3d => file_curr%Read_2d_time(varname, fill_array=monthly_buf(:,:, 1:12), start_time_step=1, end_time_step=12, &
      & val_min=val_min, val_max=val_max)
    ptr_2d => file_next%Read_2d_one_time(varname, time=1, fill_array=monthly_buf(:,:, 13), &
      & val_min=val_min, val_max=val_max)

    NULLIFY(ptr_2d)
    NULLIFY(ptr_3d)

  END SUBROUTINE extvar_monthly_buffer_init

  ! ======================================================================================================= !
  !> Advance a `(nproma, nblks, 0:13)` monthly buffer at a year transition Y -> Y+1.
  !>
  !> Copies December in-place, then reads the new year and the next January from file:
  !>   - `monthly_buf(:,:,0)`    <- old `monthly_buf(:,:,12)` (December of Y, free)
  !>   - `monthly_buf(:,:,1:12)` = all months of `file_curr` (year Y+1; one file read)
  !>   - `monthly_buf(:,:,13)`   = month 1 of `file_next`    (January of Y+2; one time-step read)
  !>
  SUBROUTINE extvar_monthly_buffer_advance(file_curr, file_next, varname, monthly_buf, val_min, val_max)

    TYPE(t_input_file), INTENT(inout)        :: file_curr            !< Open file for year Y+1
    TYPE(t_input_file), INTENT(inout)        :: file_next            !< Open file for year Y+2
    CHARACTER(len=*),   INTENT(in)           :: varname
    REAL(wp),           INTENT(inout)        :: monthly_buf(:,:,0:)  !< (nproma, nblks, 0:13)
    REAL(wp),           INTENT(in), OPTIONAL :: val_min
    REAL(wp),           INTENT(in), OPTIONAL :: val_max

    REAL(wp), POINTER :: ptr_2d(:,:), ptr_3d(:,:,:)

    monthly_buf(:,:, 0) = monthly_buf(:,:, 12)
    ptr_3d => file_curr%Read_2d_time(TRIM(varname), fill_array=monthly_buf(:,:, 1:12), &
      & start_time_step=1, end_time_step=12, val_min=val_min, val_max=val_max)
    NULLIFY(ptr_3d)
    ptr_2d => file_next%Read_2d_one_time(TRIM(varname), time=1, fill_array=monthly_buf(:,:, 13), &
      & val_min=val_min, val_max=val_max)
    NULLIFY(ptr_2d)

  END SUBROUTINE extvar_monthly_buffer_advance

  ! ======================================================================================================= !
  !> Mid-month linear interpolation of a monthly buffer using the current cached time weights.
  !>
  !> Generic interface; dispatch is by rank of the first two arguments:
  !>   - full-field: `monthly_buffer(:,:,0:)`, `data_out(:,:)` (INTENT OUT)
  !>   - on-chunk:   `monthly_buffer(:,0:)`,   `data_out(:)`   (INTENT INOUT, writes in-place)
  !>
  !> Interpolation weights and month indices are obtained internally from
  !> `get_time_interpolation_weights`; no weight arguments are needed at the call site.
  !>
  !> OpenACC: `monthly_buffer` and `data_out` must already be present on the device
  !> (e.g. via prior `!$ACC ENTER DATA` or `DEFAULT(PRESENT)` context).
  !>
  SUBROUTINE extvar_interpolate_monthly_2d(monthly_buffer, data_out, acc_stream)

    REAL(wp), INTENT(in)           :: monthly_buffer(:,:,0:)  !< Monthly buffer, 3rd dim bounds 0:13
    REAL(wp), INTENT(out)          :: data_out(:,:)           !< Interpolated field
    INTEGER,  INTENT(in), OPTIONAL :: acc_stream              !< OpenACC async stream (default: 1)

    REAL(wp) :: w1, w2
    INTEGER  :: n1, n2, iacc_stream, ic, ib

    iacc_stream = 1;  IF (PRESENT(acc_stream)) iacc_stream = acc_stream

    !$ACC DATA CREATE(w1, w2, n1, n2) ASYNC(iacc_stream)
    CALL get_time_interpolation_weights(w1, w2, n1, n2, iacc_stream)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) PRESENT(w1, w2, n1, n2) GANG VECTOR ASYNC(iacc_stream) COLLAPSE(2)
    DO ib = 1, SIZE(data_out, 2)
      DO ic = 1, SIZE(data_out, 1)
        data_out(ic,ib) = w1 * monthly_buffer(ic,ib,n1) + w2 * monthly_buffer(ic,ib,n2)
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC END DATA

  END SUBROUTINE extvar_interpolate_monthly_2d

  SUBROUTINE extvar_interpolate_monthly_1d(monthly_buffer, data_out, acc_stream)

    REAL(wp), INTENT(in)           :: monthly_buffer(:,0:)   !< Monthly buffer, 2nd dim bounds 0:13
    REAL(wp), INTENT(inout)        :: data_out(:)            !< Interpolated field (on-chunk)
    INTEGER,  INTENT(in), OPTIONAL :: acc_stream             !< OpenACC async stream (default: 1)

    REAL(wp) :: w1, w2
    INTEGER  :: n1, n2, iacc_stream, ic

    iacc_stream = 1;  IF (PRESENT(acc_stream)) iacc_stream = acc_stream

    !$ACC DATA CREATE(w1, w2, n1, n2) ASYNC(iacc_stream)
    CALL get_time_interpolation_weights(w1, w2, n1, n2, iacc_stream)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) PRESENT(w1, w2, n1, n2) GANG VECTOR ASYNC(iacc_stream)
    DO ic = 1, SIZE(data_out)
      data_out(ic) = w1 * monthly_buffer(ic, n1) + w2 * monthly_buffer(ic, n2)
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC END DATA

  END SUBROUTINE extvar_interpolate_monthly_1d

#endif
END MODULE mo_jsb_extvar_reader
