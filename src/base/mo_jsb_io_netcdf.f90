MODULE mo_jsb_io_netcdf
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish, message

  USE mo_jsb_io_netcdf_iface, ONLY: nf_open, nf_inq_dimid, nf_inq_varid, nf_close, &
    NF_NOWRITE, NF_NOERR, nf, t_input_file, netcdf_open_input, netcdf_close_file, &
    netcdf_read_2d, netcdf_read_2d_time, netcdf_read_2d_extdim
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_grid,           ONLY: Get_grid

  IMPLICIT NONE

  PRIVATE

!!$  PUBLIC :: t_input_file, jsb_netcdf_open_input, &
!!$            jsb_netcdf_read_2d_time, jsb_netcdf_read_2d, jsb_netcdf_read_2d_extdim, &
!!$            ncfile_has_dim, ncfile_has_var
  PUBLIC :: t_input_file, jsb_netcdf_open_input, ncfile_has_dim, ncfile_has_var


  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io_netcdf'

CONTAINS

  TYPE(t_input_file) FUNCTION jsb_netcdf_open_input(filename, grid_id)

    CHARACTER(LEN=*), INTENT(in) :: filename
    INTEGER,          INTENT(in) :: grid_id

    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':jsb_netcdf_open_input'

    grid => Get_grid(grid_id)

    jsb_netcdf_open_input = netcdf_open_input(filename, patch=grid%patch, decomp=grid%decomposition)

  END FUNCTION jsb_netcdf_open_input

!!$  SUBROUTINE jsb_netcdf_close_file(input_file)
!!$
!!$    TYPE(t_input_file), INTENT(inout) :: input_file
!!$
!!$    CALL netcdf_close_file(input_file)
!!$
!!$  END SUBROUTINE jsb_netcdf_close_file
!!$
!!$  FUNCTION jsb_netcdf_read_2d(input_file, variable_name, fill_array)
!!$
!!$    TYPE(t_input_file), INTENT(inout) :: input_file
!!$    CHARACTER(LEN=*),   INTENT(in) :: variable_name
!!$    REAL(wp), TARGET, OPTIONAL     :: fill_array(:,:)
!!$    REAL(wp), POINTER              :: jsb_netcdf_read_2d(:,:)
!!$
!!$    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':jsb_netcdf_read_2d'
!!$
!!$    jsb_netcdf_read_2d => netcdf_read_2d( &
!!$      input_file, variable_name,          &
!!$      fill_array=fill_array)
!!$
!!$  END FUNCTION jsb_netcdf_read_2d
!!$
!!$  FUNCTION jsb_netcdf_read_2d_time(input_file, variable_name, fill_array, &
!!$    start_time_step, end_time_step)
!!$
!!$    TYPE(t_input_file), INTENT(inout) :: input_file
!!$    CHARACTER(LEN=*),   INTENT(in) :: variable_name
!!$    REAL(wp), TARGET, OPTIONAL     :: fill_array(:,:,:)
!!$    INTEGER, OPTIONAL,  INTENT(in) :: start_time_step, end_time_step
!!$    REAL(wp), POINTER              :: jsb_netcdf_read_2d_time(:,:,:)
!!$
!!$    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':jsb_netcdf_read_2d_time'
!!$
!!$    jsb_netcdf_read_2d_time => netcdf_read_2d_time( &
!!$      input_file, variable_name,                    &
!!$      fill_array=fill_array,                        &
!!$      start_time_step=start_time_step,              &
!!$      end_time_step=end_time_step)
!!$
!!$  END FUNCTION jsb_netcdf_read_2d_time
!!$
!!$  FUNCTION jsb_netcdf_read_2d_extdim(input_file, variable_name, fill_array, &
!!$    start_extdim, end_extdim, extdim_name)
!!$
!!$    TYPE(t_input_file),         INTENT(inout) :: input_file
!!$    CHARACTER(LEN=*),           INTENT(in) :: variable_name
!!$    REAL(wp), TARGET, OPTIONAL             :: fill_array(:,:,:)
!!$    INTEGER, OPTIONAL,          INTENT(in) :: start_extdim, end_extdim
!!$    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: extdim_name
!!$    REAL(wp), POINTER                      :: jsb_netcdf_read_2d_extdim(:,:,:)
!!$
!!$    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':jsb_netcdf_read_2d_extdim'
!!$
!!$    jsb_netcdf_read_2d_extdim => netcdf_read_2d_extdim( &
!!$      input_file, variable_name,                        &
!!$      fill_array=fill_array,                            &
!!$      start_extdim=start_extdim,                        &
!!$      end_extdim=end_extdim,                            &
!!$      extdim_name=extdim_name)
!!$
!!$  END FUNCTION jsb_netcdf_read_2d_extdim

  LOGICAL FUNCTION ncfile_has_dim(filename, dimname)

    CHARACTER(LEN=*), INTENT(in) :: filename
    CHARACTER(LEN=*), INTENT(in) :: dimname

    INTEGER :: ncid, dimid, status

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':ncfile_has_dim'

    CALL nf(nf_open(TRIM(filename), NF_NOWRITE, ncid), routine)

    status = nf_inq_dimid(ncid, TRIM(dimname), dimid)
    ncfile_has_dim = (status == NF_NOERR)

    CALL nf(nf_close(ncid), routine)

  END FUNCTION ncfile_has_dim

  LOGICAL FUNCTION ncfile_has_var(filename, varname)

    CHARACTER(LEN=*), INTENT(in) :: filename
    CHARACTER(LEN=*), INTENT(in) :: varname

    INTEGER :: ncid, varid, status

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':ncfile_has_var'

    CALL nf(nf_open(TRIM(filename), NF_NOWRITE, ncid), routine)

    status = nf_inq_varid(ncid, TRIM(varname), varid)
    ncfile_has_var = (status == NF_NOERR)

    CALL nf(nf_close(ncid), routine)

  END FUNCTION ncfile_has_var

#endif
END MODULE mo_jsb_io_netcdf
