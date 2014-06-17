!>
!! @brief Contains the interface to io 
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_io_netcdf_iface
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message
  USE mo_jsb_domain_iface,   ONLY: t_patch, pe_decomposed
  USE mo_jsb_parallel_iface, ONLY: my_process_is_mpi_parallel, my_process_is_stdio, p_bcast, p_io
#ifdef __ICON__
  USE mo_netcdf_read, ONLY: netcdf_read_oncells_2d,       &
                            netcdf_read_oncells_2d_time,  &
                            netcdf_read_oncells_2d_extdim
#else
  USE mo_tr_scatter, ONLY: scatter_field
#endif

  IMPLICIT NONE

  PUBLIC

  INTERFACE netcdf_read_2d
    MODULE PROCEDURE netcdf_read_real_2d
  END INTERFACE netcdf_read_2d

  INTERFACE netcdf_read_2d_time
    MODULE PROCEDURE netcdf_read_real_2d_time
  END INTERFACE netcdf_read_2d_time

  INTERFACE netcdf_read_2d_extdim
    MODULE PROCEDURE netcdf_read_real_2d_extdim
  END INTERFACE netcdf_read_2d_extdim

!#ifndef __ICON__
  INCLUDE 'netcdf.inc'
!#endif

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io_netcdf_iface'

CONTAINS

  FUNCTION netcdf_read_real_2d(filename, variable_name, fill_array, patch, decomp)

    CHARACTER(LEN=*), INTENT(in) :: filename
    CHARACTER(LEN=*), INTENT(in) :: variable_name
    REAL(wp), TARGET, OPTIONAL   :: fill_array(:,:)
    TYPE(t_patch),       INTENT(in) :: patch
    TYPE(pe_decomposed), INTENT(in) :: decomp
    REAL(wp), POINTER            :: netcdf_read_real_2d(:,:)

    REAL(wp), POINTER :: tmp_array(:,:)
    INTEGER :: file_id, var_id

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':netcdf_read_real_2d'

    NULLIFY(netcdf_read_real_2d)

    IF (patch%id > 0) THEN
#ifdef __ICON__
      netcdf_read_real_2d => &
        netcdf_read_oncells_2d(TRIM(filename), TRIM(variable_name), fill_array, patch)
#else
      CALL finish(TRIM(routine), 'Incompatible combination patch/no-ICON')
#endif
    ELSE IF (decomp%nlon > 0) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'Incompatible combination decomp/ICON')
#else
      ALLOCATE(tmp_array(decomp%nlon, decomp%nlat))

      IF (my_process_is_stdio()) THEN
        CALL nf(nf_open(TRIM(filename), NF_NOWRITE, file_id), TRIM(routine)//':'//TRIM(filename))

        CALL nf(nf_inq_varid(file_id, TRIM(variable_name), var_id), TRIM(routine))

        CALL nf(nf_get_var_double(file_id, var_id, tmp_array), TRIM(routine))
      END IF

      IF (PRESENT(fill_array)) THEN
        netcdf_read_real_2d => fill_array
      ELSE
        ALLOCATE(netcdf_read_real_2d(decomp%nproma, decomp%ngpblks))
        netcdf_read_real_2d(:,:) = 0._wp
      END IF

      CALL scatter_field(tmp_array, netcdf_read_real_2d)

      DEALLOCATE(tmp_array)

      IF (my_process_is_stdio()) CALL nf(nf_close(file_id), TRIM(routine))
#endif
    ELSE
#ifdef __ICON__
      CALL finish(TRIM(routine), 'No patch')
#else
      CALL finish(TRIM(routine), 'No decomposition')
#endif
    END IF

  END FUNCTION netcdf_read_real_2d

  FUNCTION netcdf_read_real_2d_time(filename, variable_name, fill_array, patch, decomp, &
    start_time_step, end_time_step)

    CHARACTER(LEN=*), INTENT(in) :: filename
    CHARACTER(LEN=*), INTENT(in) :: variable_name
    REAL(wp), TARGET, OPTIONAL   :: fill_array(:,:,:)
    TYPE(t_patch),       INTENT(in) :: patch
    TYPE(pe_decomposed), INTENT(in) :: decomp
    INTEGER, OPTIONAL, INTENT(in) :: start_time_step, end_time_step
    REAL(wp), POINTER            :: netcdf_read_real_2d_time(:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':netcdf_read_real_2d_time'

    netcdf_read_real_2d_time => netcdf_read_real_2d_extdim( &
      filename, variable_name, fill_array, patch, decomp,   &
      start_extdim=start_time_step, end_extdim=end_time_step, extdim_name="time")
                                 
  END FUNCTION netcdf_read_real_2d_time

  FUNCTION netcdf_read_real_2d_extdim(filename, variable_name, fill_array, patch, decomp, &
    start_extdim, end_extdim, extdim_name)

    CHARACTER(LEN=*), INTENT(in) :: filename
    CHARACTER(LEN=*), INTENT(in) :: variable_name
    REAL(wp), TARGET, OPTIONAL   :: fill_array(:,:,:)
    TYPE(t_patch),       INTENT(in) :: patch
    TYPE(pe_decomposed), INTENT(in) :: decomp
    INTEGER, OPTIONAL, INTENT(in) :: start_extdim, end_extdim
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: extdim_name
    REAL(wp), POINTER            :: netcdf_read_real_2d_extdim(:,:,:)

    INTEGER :: i, start_allocated_step, end_allocated_step, steps, start_step, end_step
    REAL(wp), POINTER :: tmp_array(:,:,:), ptr(:,:)
    INTEGER :: file_id, var_id, var_ndims, var_dim_ids(1:3), var_extdim_size, start_read_index(1:3), count_read_index(1:3)
    CHARACTER(LEN=NF_MAX_NAME) :: var_extdim_name

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':netcdf_read_real_2d_extdim'

    NULLIFY(netcdf_read_real_2d_extdim)

    IF (patch%id > 0) THEN
#ifdef __ICON__
      netcdf_read_real_2d_extdim => &
        netcdf_read_oncells_2d_extdim(TRIM(filename), TRIM(variable_name), fill_array, patch, &
                                      start_extdim=start_extdim, end_extdim=end_extdim, extdim_name=extdim_name)
#else
      CALL finish(TRIM(routine), 'Incompatible combination patch/no-ICON')
#endif
    ELSE IF (decomp%nlon > 0) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'Incompatible combination decomp/ICON')
#else
      IF (my_process_is_stdio()) THEN
        CALL nf(nf_open(TRIM(filename), NF_NOWRITE, file_id), TRIM(routine)//':'//TRIM(filename))

        CALL nf(nf_inq_varid(file_id, TRIM(variable_name), var_id), TRIM(routine))

        CALL nf(nf_inq_varndims(file_id, var_id, var_ndims), TRIM(routine)//':'//TRIM(variable_name))
        IF (var_ndims /= 3) CALL finish(TRIM(routine), 'Dimension mismatch')

        CALL nf(nf_inq_vardimid(file_id, var_id, var_dim_ids(:)), TRIM(routine)//':'//TRIM(variable_name))
        CALL nf(nf_inq_dimlen(file_id, var_dim_ids(3), var_extdim_size), TRIM(routine)//':'//TRIM(variable_name))

        IF (PRESENT(extdim_name)) THEN
          CALL nf(nf_inq_dimname(file_id, var_dim_ids(3), var_extdim_name), TRIM(routine)//':'//TRIM(variable_name))
          IF (TRIM(extdim_name) /= TRIM(var_extdim_name)) CALL finish(TRIM(routine), 'extdim_name mismatch')
        END IF
      END IF

      CALL p_bcast(var_extdim_size, p_io)

      IF (PRESENT(start_extdim)) THEN
        start_step = start_extdim
      ELSE
        start_step = 1
      END IF
      IF (PRESENT(end_extdim)) THEN
        end_step = end_extdim
      ELSE
        end_step = var_extdim_size
      END IF
      steps = end_step - start_step + 1
      IF (steps < 1) CALL finish(TRIM(routine), 'extdim_size < 1')

      IF (PRESENT(fill_array)) THEN
        netcdf_read_real_2d_extdim => fill_array
      ELSE
        ALLOCATE(netcdf_read_real_2d_extdim(decomp%nproma, decomp%ngpblks, steps))
        netcdf_read_real_2d_extdim(:,:,:) = 0._wp
      END IF

      start_allocated_step = LBOUND(netcdf_read_real_2d_extdim, 3)
      end_allocated_step = UBOUND(netcdf_read_real_2d_extdim, 3)
      ALLOCATE(tmp_array(decomp%nlon, decomp%nlat, start_allocated_step:end_allocated_step))
      IF (SIZE(netcdf_read_real_2d_extdim,3) < steps) CALL finish(TRIM(routine), 'allocated size < steps')

      IF (my_process_is_stdio()) THEN
        start_read_index = (/ 1, 1, start_step /)
        count_read_index = (/ decomp%nlon, decomp%nlat, steps /)
        CALL nf(nf_get_vara_double(file_id, var_id, start_read_index, count_read_index, tmp_array), &
          TRIM(routine)//':'//TRIM(variable_name))
      END IF

      DO i=start_allocated_step,end_allocated_step
        IF (my_process_is_stdio()) ptr => tmp_array(:,:,i)
        CALL scatter_field(ptr, netcdf_read_real_2d_extdim(:,:,i))
      END DO

      DEALLOCATE(tmp_array)

      IF (my_process_is_stdio()) CALL nf(nf_close(file_id), TRIM(routine))
#endif
    ELSE
#ifdef __ICON__
      CALL finish(TRIM(routine), 'No patch')
#else
      CALL finish(TRIM(routine), 'No decomposition')
#endif
    END IF

  END FUNCTION netcdf_read_real_2d_extdim

  SUBROUTINE nf(status, routine, warnonly, silent)
    
    USE mo_exception, ONLY: em_warn

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
        CALL message( TRIM(routine)//' netCDF error', NF_STRERROR(status), &
          & level=em_warn)
      ELSE
        CALL finish( TRIM(routine)//' netCDF error', NF_STRERROR(status))
      ENDIF
    ENDIF
    
  END SUBROUTINE nf

#endif
END MODULE mo_jsb_io_netcdf_iface
