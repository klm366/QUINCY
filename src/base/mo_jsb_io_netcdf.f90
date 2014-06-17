MODULE mo_jsb_io_netcdf
#ifndef __NO_JSBACH__

  USE mo_exception, ONLY: finish, message, em_warn
  USE mo_io_units,         ONLY: filename_max

#if !defined(NOMPI) && defined(__ICON__)
  USE mo_read_netcdf_parallel, ONLY:               &
    nf_nowrite, nf_global, nf_noerr, nf_strerror,  &
    nf_open            => p_nf_open,               &
    nf_close           => p_nf_close,              &
    nf_inq_dimid       => p_nf_inq_dimid,          &
    nf_get_att_int     => p_nf_get_att_int ,       &
    nf_get_att_double  => p_nf_get_att_double
#endif

  IMPLICIT NONE

  PUBLIC

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io_netcdf'

#if defined(NOMPI) || !defined(__ICON__)
  INCLUDE 'netcdf.inc'
#endif

CONTAINS

  LOGICAL FUNCTION ncfile_has_dim(filename, dimname)

    CHARACTER(LEN=*), INTENT(in) :: filename
    CHARACTER(LEN=*), INTENT(in) :: dimname

    INTEGER :: ncid, dimid, status

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':ncfile_has_dim'

    CALL nf(nf_open(TRIM(filename), NF_NOWRITE, ncid), routine)

    status = nf_inq_dimid(ncid, TRIM(dimname), dimid)
    ncfile_has_dim = (status == nf_noerr)

    CALL nf(nf_close(ncid), routine)

  END FUNCTION ncfile_has_dim

  SUBROUTINE nf(status, routine, warnonly, silent)
    
    INTEGER, INTENT(in)           :: status
    CHARACTER(len=*), INTENT(in)  :: routine
    LOGICAL, INTENT(in), OPTIONAL :: warnonly
    LOGICAL, INTENT(in), OPTIONAL :: silent
    
    LOGICAL :: lwarnonly, lsilent
    
    lwarnonly = .FALSE.
    lsilent   = .FALSE.
    IF(PRESENT(warnonly)) lwarnonly = .TRUE.
    IF(PRESENT(silent))   lsilent   = silent
    
    IF (lsilent) RETURN

    IF (status /= nf_noerr) THEN
      IF (lwarnonly) THEN
        CALL message( TRIM(routine)//' netCDF error', nf_strerror(status), &
          & level=em_warn)
      ELSE
        CALL finish( TRIM(routine)//' netCDF error', nf_strerror(status))
      ENDIF
    ENDIF
    
  END SUBROUTINE nf

#endif
END MODULE mo_jsb_io_netcdf
