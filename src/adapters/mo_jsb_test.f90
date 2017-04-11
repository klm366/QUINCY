!>
!! @brief Module with test functions and parameters
!!
!! --> ATTENTION: currently only implemented for ECHAM: therefore it is located in adapters.
!!     TODO: implement for use with ICON and then move to e.g. atmland folder?!
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! Most parts copied from jsbach3 mo_test.f90
!!
MODULE mo_jsb_test
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_jsb_time,           ONLY: is_time_experiment_start

!#ifdef __ICON__
!! TODO: fillvalue_double is private in icon
!  USE mo_write_netcdf,       ONLY: fillvalue_double
!#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: write_interface_variables, write_interface_vars, read_interface_vars, read_interface_variables, fillValue

  ! For test with direct echam output
  LOGICAL :: read_interface_vars, write_interface_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_test'

!#ifndef __ICON__
  REAL(wp), PARAMETER :: fillValue = -9.e+33_wp
!#else
!  REAL(wp), PARAMETER :: fillValue = fillvalue_double
!#endif

  CHARACTER(len=*), PARAMETER :: fillValueAttrName = '_FillValue'

CONTAINS

  !>
  !! @brief write all variables from the call to the full jsbach interface
  !!
  !! Currently only working for ECHAM
  !!
  !! @par History
  !! Initial version copied from jsbach3 mo_test.f90: write_interface_variables
  !!
  SUBROUTINE write_interface_variables( ics, ice, iblk, t_air, q_air, rain, snow, wind_air, &
    & wind_10m, lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, press_srf, drag_srf, &
    & t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle )

    USE mo_exception,     ONLY: finish, message_text

    USE mo_jsb_parallel_iface, ONLY: my_process_is_mpi_parallel

#ifndef __ICON__
    USE mo_time_control, ONLY: lresume, get_date_components, next_date
    USE mo_filename,     ONLY: compose_filenames, standard_output_file
    USE mo_netcdf,       ONLY: nf_check

    USE mo_decomposition, ONLY: local_decomposition
    USE mo_gaussgrid,     ONLY: philon,philat
#endif

    INTEGER,  INTENT(in) :: ics, ice, iblk
    REAL(wp), INTENT(in) :: t_air(:)
    REAL(wp), INTENT(in) :: q_air(:)
    REAL(wp), INTENT(in) :: rain(:)
    REAL(wp), INTENT(in) :: snow(:)
    REAL(wp), INTENT(in) :: wind_air(:)
    REAL(wp), INTENT(in) :: wind_10m(:)
    REAL(wp), INTENT(in) :: lw_srf_down(:)
    REAL(wp), INTENT(in) :: swvis_srf_down(:)
    REAL(wp), INTENT(in) :: swnir_srf_down(:)
    REAL(wp), INTENT(in) :: swpar_srf_down(:)
    REAL(wp), INTENT(in) :: press_srf(:)
    REAL(wp), INTENT(in) :: drag_srf(:)
    REAL(wp), INTENT(in) :: t_acoef(:)
    REAL(wp), INTENT(in) :: t_bcoef(:)
    REAL(wp), INTENT(in) :: q_acoef(:)
    REAL(wp), INTENT(in) :: q_bcoef(:)
    REAL(wp), INTENT(in) :: pch(:)
    REAL(wp), INTENT(in) :: cos_zenith_angle(:)
  !----------------------------------------------------------------------------

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':write_interface_variables'

#ifndef __ICON__
    INCLUDE 'netcdf.inc'

    INTEGER, SAVE :: ncid, time_id          ! netcdf file id and id of dim time
    INTEGER :: lon_dim, lat_dim, time_dim   ! dimension IDs
    INTEGER :: lon_id, lat_id               ! other dimension variable IDs
    INTEGER :: dim3(3)
    INTEGER :: tstart(1), tcount(1), fstart(3), fcount(3)

    INTEGER  :: year, month, day, hour, minute, second
    INTEGER  :: nlon, nlat, lonCount, nproma
    REAL(wp) :: date

    CHARACTER(len=nf_max_name), SAVE :: filename

    INTEGER, SAVE :: icount = 0             ! time step counter

    INTEGER, SAVE :: t_air_id
    INTEGER, SAVE :: q_air_id
    INTEGER, SAVE :: rain_id
    INTEGER, SAVE :: snow_id
    INTEGER, SAVE :: wind_air_id
    INTEGER, SAVE :: wind_10m_id
    INTEGER, SAVE :: lw_srf_down_id
    INTEGER, SAVE :: swvis_srf_down_id
    INTEGER, SAVE :: swnir_srf_down_id
    INTEGER, SAVE :: swpar_srf_down_id
    INTEGER, SAVE :: press_srf_id
    INTEGER, SAVE :: drag_srf_id
    INTEGER, SAVE :: t_acoef_id
    INTEGER, SAVE :: t_bcoef_id
    INTEGER, SAVE :: q_acoef_id
    INTEGER, SAVE :: q_bcoef_id
    INTEGER, SAVE :: pch_id
    INTEGER, SAVE :: cos_zenith_angle_id

    nlon = local_decomposition%nlon
    nlat = local_decomposition%nlat
    nproma = ice

    !-- as in jsbach3 standalone, this routine is not implemented to work in parallel
    IF (my_process_is_mpi_parallel()) THEN
      CALL finish(TRIM(routine), 'writing interface variables not possible in a parallel run')
    END IF
    IF (.NOT. nproma == nlon) THEN
      WRITE (message_text,*) 'Writing interface variables only possible if nproma == nlon; here: ', nproma
      CALL finish(TRIM(routine), message_text)
    END IF

    !-- Create netcdf file
    IF ((is_time_experiment_start() .OR. lresume) .AND. iblk == 1) THEN
      CALL compose_filenames
      filename = TRIM(standard_output_file)//'_interface.nc'

      CALL nf_check(nf_create(TRIM(filename), NF_CLOBBER, ncid))

      ! define dimensions
      CALL nf_check(nf_def_dim(ncid, 'lon', nlon, lon_dim))
      CALL nf_check(nf_def_dim(ncid, 'lat', nlat, lat_dim))
      CALL nf_check(nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim))

      ! define variables
      CALL nf_check(nf_def_var(ncid, 'lon',  NF_DOUBLE, 1, (/lon_dim/),  lon_id))
      CALL nf_check(nf_def_var(ncid, 'lat',  NF_DOUBLE, 1, (/lat_dim/),  lat_id))
      CALL nf_check(nf_def_var(ncid, 'time', NF_DOUBLE, 1, (/time_dim/), time_id))

      dim3(1:3) = (/ lon_dim, lat_dim, time_dim /)

      CALL nf_check(nf_def_var(ncid, 't_air',            NF_DOUBLE, 3, dim3, t_air_id))
      CALL nf_check(nf_def_var(ncid, 'q_air',            NF_DOUBLE, 3, dim3, q_air_id))
      CALL nf_check(nf_def_var(ncid, 'rain',             NF_DOUBLE, 3, dim3, rain_id))
      CALL nf_check(nf_def_var(ncid, 'snow',             NF_DOUBLE, 3, dim3, snow_id))
      CALL nf_check(nf_def_var(ncid, 'wind_air',         NF_DOUBLE, 3, dim3, wind_air_id))
      CALL nf_check(nf_def_var(ncid, 'wind_10m',         NF_DOUBLE, 3, dim3, wind_10m_id))
      CALL nf_check(nf_def_var(ncid, 'lw_srf_down',      NF_DOUBLE, 3, dim3, lw_srf_down_id))
      CALL nf_check(nf_def_var(ncid, 'swvis_srf_down',   NF_DOUBLE, 3, dim3, swvis_srf_down_id))
      CALL nf_check(nf_def_var(ncid, 'swnir_srf_down',   NF_DOUBLE, 3, dim3, swnir_srf_down_id))
      CALL nf_check(nf_def_var(ncid, 'swpar_srf_down',   NF_DOUBLE, 3, dim3, swpar_srf_down_id))
      CALL nf_check(nf_def_var(ncid, 'press_srf',        NF_DOUBLE, 3, dim3, press_srf_id))
      CALL nf_check(nf_def_var(ncid, 'drag_srf',         NF_DOUBLE, 3, dim3, drag_srf_id))
      CALL nf_check(nf_def_var(ncid, 't_acoef',          NF_DOUBLE, 3, dim3, t_acoef_id))
      CALL nf_check(nf_def_var(ncid, 't_bcoef',          NF_DOUBLE, 3, dim3, t_bcoef_id))
      CALL nf_check(nf_def_var(ncid, 'q_acoef',          NF_DOUBLE, 3, dim3, q_acoef_id))
      CALL nf_check(nf_def_var(ncid, 'q_bcoef',          NF_DOUBLE, 3, dim3, q_bcoef_id))
      CALL nf_check(nf_def_var(ncid, 'pch',              NF_DOUBLE, 3, dim3, pch_id))
      CALL nf_check(nf_def_var(ncid, 'cos_zenith_angle', NF_DOUBLE, 3, dim3, cos_zenith_angle_id))

      CALL nf_check(nf_put_att_double(ncid, t_air_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, q_air_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, rain_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, snow_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, wind_air_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, wind_10m_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, lw_srf_down_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, swvis_srf_down_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, swnir_srf_down_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, swpar_srf_down_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, press_srf_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, drag_srf_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, t_acoef_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, t_bcoef_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, q_acoef_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, q_bcoef_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, pch_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))
      CALL nf_check(nf_put_att_double(ncid, cos_zenith_angle_id, fillValueAttrName, NF_DOUBLE, 1, fillValue))

      ! assign global attributes
      CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'conventions', 6, 'CF-1.0'))
      CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'source', 15, 'ECHAM6'))
      CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'institution', 36, 'Max-Planck-Institute for Meteorology'))

      ! leave define mode
      CALL nf_check(nf_enddef(ncid))

      CALL nf_check(nf_put_var_double(ncid, lat_id, philat(:)))
      CALL nf_check(nf_put_var_double(ncid, lon_id, philon(:)))


    ELSE
      CALL nf_check(nf_open(TRIM(filename), NF_WRITE, ncid))
    END IF

    IF (iblk == 1) THEN
      !-- write data of the current time step
      icount = icount+1

      CALL get_date_components(next_date, year, month, day, hour, minute, second)
      date = ISIGN(1,year)*(IABS(year))*10000 + month*100 + day  &
              + (hour*3600 + minute*60 + second)/86400._wp

      tstart(1) = icount
      tcount(1) = 1
      CALL nf_check(nf_put_vara_double(ncid, time_id, tstart, tcount, date))
    ENDIF

    lonCount = ice-ics + 1

    fstart(1:3) = (/ ics, iblk, icount /)
    fcount(1:3) = (/ lonCount, 1, 1 /)

    CALL nf_check(nf_put_vara_double(ncid, t_air_id,            fstart, fcount, t_air))
    CALL nf_check(nf_put_vara_double(ncid, q_air_id,            fstart, fcount, q_air))
    CALL nf_check(nf_put_vara_double(ncid, rain_id,             fstart, fcount, rain))
    CALL nf_check(nf_put_vara_double(ncid, snow_id,             fstart, fcount, snow))
    CALL nf_check(nf_put_vara_double(ncid, wind_air_id,         fstart, fcount, wind_air))
    CALL nf_check(nf_put_vara_double(ncid, wind_10m_id,         fstart, fcount, wind_10m))
    CALL nf_check(nf_put_vara_double(ncid, lw_srf_down_id,      fstart, fcount, lw_srf_down))
    CALL nf_check(nf_put_vara_double(ncid, swvis_srf_down_id,   fstart, fcount, swvis_srf_down))
    CALL nf_check(nf_put_vara_double(ncid, swnir_srf_down_id,   fstart, fcount, swnir_srf_down))
    CALL nf_check(nf_put_vara_double(ncid, swpar_srf_down_id,   fstart, fcount, swpar_srf_down))
    CALL nf_check(nf_put_vara_double(ncid, press_srf_id,        fstart, fcount, press_srf))
    CALL nf_check(nf_put_vara_double(ncid, drag_srf_id,         fstart, fcount, drag_srf))
    CALL nf_check(nf_put_vara_double(ncid, t_acoef_id,          fstart, fcount, t_acoef))
    CALL nf_check(nf_put_vara_double(ncid, t_bcoef_id,          fstart, fcount, t_bcoef))
    CALL nf_check(nf_put_vara_double(ncid, q_acoef_id,          fstart, fcount, q_acoef))
    CALL nf_check(nf_put_vara_double(ncid, q_bcoef_id,          fstart, fcount, q_bcoef))
    CALL nf_check(nf_put_vara_double(ncid, pch_id,              fstart, fcount, pch))
    CALL nf_check(nf_put_vara_double(ncid, cos_zenith_angle_id, fstart, fcount, cos_zenith_angle))

    !-- close the file
    CALL nf_check(nf_close(ncid))

#else
    CALL finish(TRIM(routine), 'Not implemented for use with ICON')
#endif

  END SUBROUTINE write_interface_variables

  !>
  !! @brief read all variables required for a call to the full jsbach interface
  !!
  !! Currently only working for ECHAM
  !!
  !! @par History
  !! Initial version copied from jsbach3 mo_test.f90: read_interface_variables
  !!
  SUBROUTINE read_interface_variables( ics, ice, iblk, ncid, t_air, q_air, rain, snow, wind_air, wind_10m, &
             & lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, press_srf, drag_srf, &
             & t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle )

    USE mo_exception,     ONLY: finish, message_text

    USE mo_jsb_parallel_iface, ONLY: my_process_is_mpi_parallel

#ifndef __ICON__
    USE mo_time_control, ONLY: lresume, get_date_components, next_date
    USE mo_netcdf,       ONLY: nf_check

    USE mo_decomposition, ONLY: local_decomposition
#endif

    INTEGER,  INTENT(in)  :: ics, ice, iblk
    INTEGER,  INTENT(inout) :: ncid
    REAL(wp), INTENT(out) :: t_air(:)
    REAL(wp), INTENT(out) :: q_air(:)
    REAL(wp), INTENT(out) :: rain(:)
    REAL(wp), INTENT(out) :: snow(:)
    REAL(wp), INTENT(out) :: wind_air(:)
    REAL(wp), INTENT(out) :: wind_10m(:)
    REAL(wp), INTENT(out) :: lw_srf_down(:)
    REAL(wp), INTENT(out) :: swvis_srf_down(:)
    REAL(wp), INTENT(out) :: swnir_srf_down(:)
    REAL(wp), INTENT(out) :: swpar_srf_down(:)
    REAL(wp), INTENT(out) :: press_srf(:)
    REAL(wp), INTENT(out) :: drag_srf(:)
    REAL(wp), INTENT(out) :: t_acoef(:)
    REAL(wp), INTENT(out) :: t_bcoef(:)
    REAL(wp), INTENT(out) :: q_acoef(:)
    REAL(wp), INTENT(out) :: q_bcoef(:)
    REAL(wp), INTENT(out) :: pch(:)
    REAL(wp), INTENT(out) :: cos_zenith_angle(:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':read_interface_variables'
    CHARACTER(len=*), PARAMETER :: filename = 'interface_variables.nc' 

#ifndef __ICON__
    INCLUDE 'netcdf.inc'

    INTEGER, SAVE :: lon_dim, lat_dim, time_dim   ! dimensions
    INTEGER, SAVE :: lon_id, lat_id, time_id      ! dimension IDs

    INTEGER  :: fstart(3), fcount(3)
    INTEGER  :: nlon, nlat, nproma, currentTimestep 
    INTEGER  :: year, month, day, hour, minute, second
    REAL(wp) :: date
    REAL(wp), SAVE, ALLOCATABLE :: timevals(:)

    INTEGER, SAVE :: t_air_id
    INTEGER, SAVE :: q_air_id
    INTEGER, SAVE :: rain_id
    INTEGER, SAVE :: snow_id
    INTEGER, SAVE :: wind_air_id
    INTEGER, SAVE :: wind_10m_id
    INTEGER, SAVE :: lw_srf_down_id
    INTEGER, SAVE :: swvis_srf_down_id
    INTEGER, SAVE :: swnir_srf_down_id
    INTEGER, SAVE :: swpar_srf_down_id
    INTEGER, SAVE :: press_srf_id
    INTEGER, SAVE :: drag_srf_id
    INTEGER, SAVE :: t_acoef_id
    INTEGER, SAVE :: t_bcoef_id
    INTEGER, SAVE :: q_acoef_id
    INTEGER, SAVE :: q_bcoef_id
    INTEGER, SAVE :: pch_id
    INTEGER, SAVE :: cos_zenith_angle_id

    nlon = local_decomposition%nlon
    nlat = local_decomposition%nlat
    nproma = local_decomposition%nproma

    !-- as in jsbach3 standalone, this routine is not implemented to work in parallel
    IF (my_process_is_mpi_parallel()) THEN
      CALL finish(TRIM(routine), 'reading interface variables not implemented for parallel runs')
    END IF
    IF (.NOT. nproma == nlon) THEN
      WRITE (message_text,*) 'reading interface variables only implemented when nlon=', nlon, ', and nproma=', nproma, ' match'
      CALL finish (TRIM(routine), message_text)
    END IF

    IF ((is_time_experiment_start() .OR. lresume) .AND. iblk == 1) THEN
      CALL nf_check(nf_open(TRIM(filename), NF_NOWRITE, ncid))

      CALL nf_check(nf_inq_dimid(ncid, 'lon', lon_id))
      CALL nf_check(nf_inq_dimlen(ncid, lon_id, lon_dim))
      IF (lon_dim /= nlon) THEN
          WRITE (message_text,*) 'longitude dimensions do not match: nlon=', nlon, ', nlon in file: ', lon_dim
          CALL finish (TRIM(routine), message_text)
      END IF

      CALL nf_check(nf_inq_dimid(ncid, 'lat', lat_id))
      CALL nf_check(nf_inq_dimlen(ncid, lat_id, lat_dim))
      IF (lat_dim /= nlat) THEN
          WRITE (message_text,*) 'longitude dimensions do not match: nlon=', nlon, ', nlon in file: ', lat_dim
          CALL finish (TRIM(routine), message_text)
      END IF
      CALL nf_check(nf_inq_dimid(ncid, 'time', time_id))
      CALL nf_check(nf_inq_dimlen(ncid, time_id, time_dim))
      CALL nf_check(nf_inq_varid(ncid, 'time', time_id))
      ALLOCATE(timevals(time_dim))
      CALL nf_check(nf_get_var(ncid, time_id, timevals))

      CALL nf_check(nf_inq_varid(ncid, 't_air',            t_air_id))
      CALL nf_check(nf_inq_varid(ncid, 'q_air',            q_air_id))
      CALL nf_check(nf_inq_varid(ncid, 'rain',             rain_id))
      CALL nf_check(nf_inq_varid(ncid, 'snow',             snow_id))
      CALL nf_check(nf_inq_varid(ncid, 'wind_air',         wind_air_id))
      CALL nf_check(nf_inq_varid(ncid, 'wind_10m',         wind_10m_id))
      CALL nf_check(nf_inq_varid(ncid, 'lw_srf_down',      lw_srf_down_id))
      CALL nf_check(nf_inq_varid(ncid, 'swvis_srf_down',   swvis_srf_down_id))
      CALL nf_check(nf_inq_varid(ncid, 'swnir_srf_down',   swnir_srf_down_id))
      CALL nf_check(nf_inq_varid(ncid, 'swpar_srf_down',   swpar_srf_down_id))
      CALL nf_check(nf_inq_varid(ncid, 'press_srf',        press_srf_id))
      CALL nf_check(nf_inq_varid(ncid, 'drag_srf',         drag_srf_id))
      CALL nf_check(nf_inq_varid(ncid, 't_acoef',          t_acoef_id))
      CALL nf_check(nf_inq_varid(ncid, 't_bcoef',          t_bcoef_id))
      CALL nf_check(nf_inq_varid(ncid, 'q_acoef',          q_acoef_id))
      CALL nf_check(nf_inq_varid(ncid, 'q_bcoef',          q_bcoef_id))
      CALL nf_check(nf_inq_varid(ncid, 'pch',              pch_id))
      CALL nf_check(nf_inq_varid(ncid, 'cos_zenith_angle', cos_zenith_angle_id))
    ENDIF

    CALL get_date_components(next_date, year, month, day, hour, minute, second)
    date = REAL(ISIGN(1,year)*(IABS(year)*10000+month*100+day), wp)  &
        + REAL(hour*3600 + minute*60 + second, wp)/86400._wp

    currentTimestep = 1
    DO WHILE (currentTimestep .LE. time_dim)
      IF ( date == timevals(currentTimestep) ) EXIT
      currentTimestep = currentTimestep + 1
    END DO

    IF (currentTimestep > time_dim) THEN 
      WRITE (message_text,*) 'date ', date, ' not available in ', filename 
      CALL finish (TRIM(routine), message_text)
    END IF
 
    fstart(1:3) = (/ ics, iblk, currentTimestep /)
    fcount(1:3) = (/ ice-ics + 1, 1, 1 /)

    CALL nf_check(nf_get_vara_double(ncid, t_air_id,            fstart, fcount, t_air))
    CALL nf_check(nf_get_vara_double(ncid, q_air_id,            fstart, fcount, q_air))
    CALL nf_check(nf_get_vara_double(ncid, rain_id,             fstart, fcount, rain))
    CALL nf_check(nf_get_vara_double(ncid, snow_id,             fstart, fcount, snow))
    CALL nf_check(nf_get_vara_double(ncid, wind_air_id,         fstart, fcount, wind_air))
    CALL nf_check(nf_get_vara_double(ncid, wind_10m_id,         fstart, fcount, wind_10m))
    CALL nf_check(nf_get_vara_double(ncid, lw_srf_down_id,      fstart, fcount, lw_srf_down))
    CALL nf_check(nf_get_vara_double(ncid, swvis_srf_down_id,   fstart, fcount, swvis_srf_down))
    CALL nf_check(nf_get_vara_double(ncid, swnir_srf_down_id,   fstart, fcount, swnir_srf_down))
    CALL nf_check(nf_get_vara_double(ncid, swpar_srf_down_id,   fstart, fcount, swpar_srf_down))
    CALL nf_check(nf_get_vara_double(ncid, press_srf_id,        fstart, fcount, press_srf))
    CALL nf_check(nf_get_vara_double(ncid, drag_srf_id,         fstart, fcount, drag_srf))
    CALL nf_check(nf_get_vara_double(ncid, t_acoef_id,          fstart, fcount, t_acoef))
    CALL nf_check(nf_get_vara_double(ncid, t_bcoef_id,          fstart, fcount, t_bcoef))
    CALL nf_check(nf_get_vara_double(ncid, q_acoef_id,          fstart, fcount, q_acoef))
    CALL nf_check(nf_get_vara_double(ncid, q_bcoef_id,          fstart, fcount, q_bcoef))
    CALL nf_check(nf_get_vara_double(ncid, pch_id,              fstart, fcount, pch)) 
    CALL nf_check(nf_get_vara_double(ncid, cos_zenith_angle_id, fstart, fcount, cos_zenith_angle))
#else
    CALL finish(TRIM(routine), 'Not implemented for use with ICON')
#endif

  END SUBROUTINE read_interface_variables

#endif
END MODULE mo_jsb_test
