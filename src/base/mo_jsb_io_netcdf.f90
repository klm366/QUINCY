MODULE mo_jsb_io_netcdf
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish, message

  USE mo_jsb_io_netcdf_iface, ONLY: nf_open, nf_inq_dimid, nf_close, &
    NF_NOWRITE, NF_NOERR, nf,                                        &
    netcdf_read_2d, netcdf_read_2d_time, netcdf_read_2d_extdim
  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_base,          ONLY: get_model

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: jsb_netcdf_read_2d_time, jsb_netcdf_read_2d, jsb_netcdf_read_2d_extdim, &
    ncfile_has_dim


  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io_netcdf'

CONTAINS

  FUNCTION jsb_netcdf_read_2d(filename, variable_name, fill_array, model_id)

    CHARACTER(LEN=*), INTENT(in) :: filename
    CHARACTER(LEN=*), INTENT(in) :: variable_name
    REAL(wp), TARGET, OPTIONAL   :: fill_array(:,:)
    INTEGER,          INTENT(in) :: model_id
    REAL(wp), POINTER            :: jsb_netcdf_read_2d(:,:)

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':jsb_netcdf_read_2d'

    model => get_model(model_id)

    jsb_netcdf_read_2d => netcdf_read_2d( &
      filename, variable_name,            &
      fill_array=fill_array,              &
      patch=model%grid%patch,             &
      decomp=model%grid%decomposition)

  END FUNCTION jsb_netcdf_read_2d

  FUNCTION jsb_netcdf_read_2d_time(filename, variable_name, fill_array, model_id, &
    start_time_step, end_time_step)

    CHARACTER(LEN=*), INTENT(in) :: filename
    CHARACTER(LEN=*), INTENT(in) :: variable_name
    REAL(wp), TARGET, OPTIONAL   :: fill_array(:,:,:)
    INTEGER,          INTENT(in) :: model_id
    INTEGER, OPTIONAL, INTENT(in) :: start_time_step, end_time_step
    REAL(wp), POINTER            :: jsb_netcdf_read_2d_time(:,:,:)

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':jsb_netcdf_read_2d_time'

    model => get_model(model_id)

    jsb_netcdf_read_2d_time => netcdf_read_2d_time( &
      filename, variable_name,                      &
      fill_array=fill_array,                        &
      patch=model%grid%patch,                       &
      decomp=model%grid%decomposition,              &
      start_time_step=start_time_step,              &
      end_time_step=end_time_step)

  END FUNCTION jsb_netcdf_read_2d_time

  FUNCTION jsb_netcdf_read_2d_extdim(filename, variable_name, fill_array, model_id, &
    start_extdim, end_extdim, extdim_name)

    CHARACTER(LEN=*), INTENT(in) :: filename
    CHARACTER(LEN=*), INTENT(in) :: variable_name
    REAL(wp), TARGET, OPTIONAL   :: fill_array(:,:,:)
    INTEGER,          INTENT(in) :: model_id
    INTEGER, OPTIONAL, INTENT(in) :: start_extdim, end_extdim
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: extdim_name
    REAL(wp), POINTER            :: jsb_netcdf_read_2d_extdim(:,:,:)

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':jsb_netcdf_read_2d_extdim'

    model => get_model(model_id)

    jsb_netcdf_read_2d_extdim => netcdf_read_2d_extdim( &
      filename, variable_name,                          &
      fill_array=fill_array,                            &
      patch=model%grid%patch,                           &
      decomp=model%grid%decomposition,                  &
      start_extdim=start_extdim,                        &
      end_extdim=end_extdim,                            &
      extdim_name=extdim_name)

  END FUNCTION jsb_netcdf_read_2d_extdim

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

#endif
END MODULE mo_jsb_io_netcdf
