MODULE mo_jsb_io_netcdf
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish, message

  USE mo_jsb_io_netcdf_iface, ONLY: nf_open, nf_inq_dimid, nf_inq_varid, nf_close, &
    NF_NOWRITE, NF_NOERR, nf, t_input_file, netcdf_open_input, netcdf_close_file
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

    CHARACTER(LEN=*), INTENT(in)           :: filename
    INTEGER,          INTENT(in), OPTIONAL :: grid_id

    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':jsb_netcdf_open_input'

    IF (PRESENT(grid_id)) THEN
      grid => Get_grid(grid_id)
      jsb_netcdf_open_input = netcdf_open_input(filename, patch=grid%patch, decomp=grid%decomposition)
    ELSE
      jsb_netcdf_open_input = netcdf_open_input(filename)
    END IF

  END FUNCTION jsb_netcdf_open_input

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
