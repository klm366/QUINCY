MODULE mo_jsb_io_netcdf
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish, message

  USE mo_jsb_io_netcdf_iface, ONLY: nf_open, nf_inq_dimid, nf_inq_varid, nf_close, &
    & NF_NOWRITE, NF_NOERR, nf, t_input_file, netcdf_open_input, netcdf_close_file
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_grid,           ONLY: Get_grid

  IMPLICIT NONE

  PRIVATE

!!$  PUBLIC :: t_input_file, jsb_netcdf_open_input, &
!!$            jsb_netcdf_read_2d_time, jsb_netcdf_read_2d, jsb_netcdf_read_2d_extdim, &
!!$            ncfile_has_dim, ncfile_has_var
  PUBLIC :: t_input_file, jsb_netcdf_open_input

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io_netcdf'

CONTAINS

  TYPE(t_input_file) FUNCTION jsb_netcdf_open_input(filename, grid_id)

    CHARACTER(LEN=*), INTENT(in)           :: filename
    INTEGER,          INTENT(in), OPTIONAL :: grid_id

    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsb_netcdf_open_input'

    IF (PRESENT(grid_id)) THEN
      grid => Get_grid(grid_id)
      jsb_netcdf_open_input = netcdf_open_input(filename, patch=grid%patch, decomp=grid%decomposition)
    ELSE
      jsb_netcdf_open_input = netcdf_open_input(filename)
    END IF

  END FUNCTION jsb_netcdf_open_input

#endif
END MODULE mo_jsb_io_netcdf
