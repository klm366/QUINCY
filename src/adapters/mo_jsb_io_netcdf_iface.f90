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
  USE mo_io_units,           ONLY: filename_max
  USE mo_jsb_domain_iface,   ONLY: t_patch, pe_decomposed
  USE mo_jsb_parallel_iface, ONLY: my_process_is_mpi_parallel, my_process_is_stdio, p_bcast, p_io
#ifdef __ICON__
  USE mo_read_interface, ONLY: read_1D, read_2D, read_2D_time, read_2D_extdim, openInputFile, closeFile, on_cells, &
                               t_stream_id, read_netcdf_broadcast_method
#else
  USE mo_tr_scatter, ONLY: scatter_field
#endif

  IMPLICIT NONE

  PUBLIC

  TYPE t_input_file
    INTEGER           :: TYPE                ! 1: ICON infrastructure, 2: ECHAM infrastructure
    CHARACTER(len=filename_max) :: filename  ! Name of input file
#ifdef __ICON__
    TYPE(t_stream_id) :: stream_id  ! ICON stream_id structure (for > 1D)
    INTEGER :: file_id              ! netcdf file id (for 1D)
#else
    INTEGER :: file_id          ! netcdf file id
    INTEGER :: return_status    ! the latest operation return status
    INTEGER :: nlon    = -1
    INTEGER :: nlat    = -1
    INTEGER :: nproma  = -1
    INTEGER :: ngpblks = -1
#endif
  CONTAINS
    PROCEDURE :: Close          => netcdf_close_file
    PROCEDURE :: Read_1d        => netcdf_read_real_1d
    PROCEDURE :: Read_2d        => netcdf_read_real_2d
    PROCEDURE :: Read_2d_time   => netcdf_read_real_2d_time
    PROCEDURE :: Read_2d_extdim => netcdf_read_real_2d_extdim
  END TYPE t_input_file

  INTERFACE netcdf_open_input
    MODULE PROCEDURE netcdf_open_input_generic
    MODULE PROCEDURE netcdf_open_input_patch
  END INTERFACE netcdf_open_input

  ! INTERFACE netcdf_read_2d
  !   MODULE PROCEDURE netcdf_read_real_2d
  ! END INTERFACE netcdf_read_2d

  ! INTERFACE netcdf_read_2d_time
  !   MODULE PROCEDURE netcdf_read_real_2d_time
  ! END INTERFACE netcdf_read_2d_time

  ! INTERFACE netcdf_read_2d_extdim
  !   MODULE PROCEDURE netcdf_read_real_2d_extdim
  ! END INTERFACE netcdf_read_2d_extdim

!#ifndef __ICON__
  INCLUDE 'netcdf.inc'

  INTEGER, PARAMETER :: MAX_VAR_DIMS = 16 ! NF_MAX_VAR_DIMS
!#endif

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io_netcdf_iface'

CONTAINS

  TYPE(t_input_file) FUNCTION netcdf_open_input_generic(filename)

    CHARACTER(LEN=*),    INTENT(in) :: filename

    INTEGER :: file_id, return_status

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':netcdf_open_input_generic'

    netcdf_open_input_generic%filename = TRIM(filename)

#ifdef __ICON__
      netcdf_open_input_generic%type      = 1
      ! netcdf_open_input%stream_id = openInputFile(TRIM(filename), patch, read_netcdf_broadcast_method)
      netcdf_open_input_generic%file_id = openInputFile(TRIM(filename))
#else
      netcdf_open_input_generic%type = 2

      IF (my_process_is_stdio()) THEN
        return_status = nf_open(TRIM(filename), NF_NOWRITE, file_id)
        CALL nf(return_status, TRIM(routine)//':'//TRIM(filename))
      ELSE
        return_status = 1
        file_id = -1
      END IF
      netcdf_open_input_generic%file_id       = file_id
      netcdf_open_input_generic%return_status = return_status
#endif

  END FUNCTION netcdf_open_input_generic

  TYPE(t_input_file) FUNCTION netcdf_open_input_patch(filename, patch, decomp)

    CHARACTER(LEN=*),    INTENT(in) :: filename
    TYPE(t_patch),       INTENT(in) :: patch
    TYPE(pe_decomposed), INTENT(in) :: decomp

    INTEGER :: file_id, return_status

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':netcdf_open_input_patch'

    netcdf_open_input_patch%filename = TRIM(filename)

    IF (patch%id > 0) THEN
#ifdef __ICON__
      netcdf_open_input_patch%type      = 1
      ! netcdf_open_input_patch%stream_id = openInputFile(TRIM(filename), patch, read_netcdf_broadcast_method)
      netcdf_open_input_patch%stream_id = openInputFile(TRIM(filename), patch)
#else
      CALL finish(TRIM(routine), 'Incompatible combination patch/no-ICON')
#endif
    ELSE IF (decomp%nlon > 0) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'Incompatible combination decomp/ICON')
#else
      netcdf_open_input_patch%type = 2

      IF (my_process_is_stdio()) THEN
        return_status = nf_open(TRIM(filename), NF_NOWRITE, file_id)
        CALL nf(return_status, TRIM(routine)//':'//TRIM(filename))
      ELSE
        return_status = 1
        file_id = -1
      END IF
      netcdf_open_input_patch%file_id       = file_id
      netcdf_open_input_patch%nlon          = decomp%nlon
      netcdf_open_input_patch%nlat          = decomp%nlat
      netcdf_open_input_patch%nproma        = decomp%nproma
      netcdf_open_input_patch%ngpblks       = decomp%ngpblks
      netcdf_open_input_patch%return_status = return_status
#endif
    ELSE
#ifdef __ICON__
      CALL finish(TRIM(routine), 'No patch')
#else
      CALL finish(TRIM(routine), 'No decomposition')
#endif
    END IF

  END FUNCTION netcdf_open_input_patch

  SUBROUTINE netcdf_close_file(input_file)

    CLASS(t_input_file), INTENT(inout) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':netcdf_close_file'

#ifdef __ICON__
    IF (ALLOCATED(input_file%stream_id%read_info)) THEN
      CALL closeFile(input_file%stream_id)
    ELSE
      CALL closeFile(input_file%file_id)
    END IF
#else
    IF (my_process_is_stdio()) &
      CALL nf(nf_close(input_file%file_id), TRIM(routine)//':'//TRIM(input_file%filename))
#endif

  END SUBROUTINE netcdf_close_file

  SUBROUTINE netcdf_close_stream(input_file)

    CLASS(t_input_file), INTENT(inout) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':netcdf_close_stream'

#ifdef __ICON__
    CALL closeFile(input_file%stream_id)
#else
    IF (my_process_is_stdio()) &
      CALL nf(nf_close(input_file%file_id), TRIM(routine)//':'//TRIM(input_file%filename))
#endif

  END SUBROUTINE netcdf_close_stream

  FUNCTION netcdf_read_real_1d(input_file, variable_name, fill_array)

    CLASS(t_input_file), INTENT(inout) :: input_file
    CHARACTER(LEN=*),   INTENT(in)     :: variable_name
    REAL(wp), TARGET, OPTIONAL         :: fill_array(:)
    REAL(wp), POINTER                  :: netcdf_read_real_1d(:)

#ifndef __ICON__
    REAL(wp), POINTER :: tmp_array(:)
    INTEGER :: file_id, var_id, dim_id(MAX_VAR_DIMS), ndims
    INTEGER :: varid, var_type
    INTEGER :: var_size(MAX_VAR_DIMS)
#endif

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':netcdf_read_real_1d'

    NULLIFY(netcdf_read_real_1d)

    IF (input_file%type == 1) THEN
#ifdef __ICON__
      CALL read_1D(input_file%file_id, TRIM(variable_name), fill_array, netcdf_read_real_1d)
#else
      CALL finish(TRIM(routine), 'Incompatible input file type')
#endif
    ELSE IF (input_file%type == 2) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'Incompatible input file type')
#else
      file_id = input_file%file_id

      IF (my_process_is_stdio()) THEN
        CALL nf(nf_inq_varid(file_id, TRIM(variable_name), varid), TRIM(routine))
        CALL nf(nf_inq_varndims(file_id, varid, ndims), TRIM(routine))
        IF (ndims /= 1) CALL finish(TRIM(routine), 'only one variable dimension expected')
        CALL nf(nf_inq_vardimid(file_id, varid, dim_id(1:ndims)), TRIM(routine))
        CALL nf(nf_inq_dimlen(file_id, dim_id(1), var_size(1)), TRIM(routine))
      END IF
      
      CALL p_bcast(var_size(1), p_io)

      ALLOCATE(tmp_array(var_size(1)))

      IF (my_process_is_stdio()) THEN
        CALL nf(nf_inq_varid(file_id, TRIM(variable_name), var_id), TRIM(routine))

        CALL nf(nf_get_var_double(file_id, var_id, tmp_array), TRIM(routine))
      END IF
      CALL p_bcast(tmp_array, p_io)

      IF (PRESENT(fill_array)) THEN
        fill_array(:) = tmp_array(:)
        netcdf_read_real_1d => fill_array
      ELSE
        ALLOCATE(netcdf_read_real_1d(var_size(1)))
        netcdf_read_real_1d(:) = tmp_array(:)
      END IF

      DEALLOCATE(tmp_array)
#endif
    ELSE
      CALL finish(TRIM(routine), 'Input file type not recognized.')
    END IF

  END FUNCTION netcdf_read_real_1d

  FUNCTION netcdf_read_real_2d(input_file, variable_name, fill_array)

    CLASS(t_input_file), INTENT(inout) :: input_file
    CHARACTER(LEN=*),   INTENT(in)    :: variable_name
    REAL(wp), TARGET, OPTIONAL        :: fill_array(:,:)
    REAL(wp), POINTER                 :: netcdf_read_real_2d(:,:)

    REAL(wp), POINTER :: tmp_array(:,:)
    INTEGER :: file_id, var_id

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':netcdf_read_real_2d'

    NULLIFY(netcdf_read_real_2d)

    IF (input_file%type == 1) THEN
#ifdef __ICON__
      CALL read_2D(input_file%stream_id, on_cells, TRIM(variable_name), fill_array, netcdf_read_real_2d)
#else
      CALL finish(TRIM(routine), 'Incompatible input file type')
#endif
    ELSE IF (input_file%type == 2) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'Incompatible input file type')
#else
      ALLOCATE(tmp_array(input_file%nlon, input_file%nlat))

      file_id = input_file%file_id

      IF (my_process_is_stdio()) THEN
        CALL nf(nf_inq_varid(file_id, TRIM(variable_name), var_id), TRIM(routine))

        CALL nf(nf_get_var_double(file_id, var_id, tmp_array), TRIM(routine))
      END IF

      IF (PRESENT(fill_array)) THEN
        netcdf_read_real_2d => fill_array
      ELSE
        ALLOCATE(netcdf_read_real_2d(input_file%nproma, input_file%ngpblks))
        netcdf_read_real_2d(:,:) = 0._wp
      END IF

      CALL scatter_field(tmp_array, netcdf_read_real_2d)

      DEALLOCATE(tmp_array)
#endif
    ELSE
      CALL finish(TRIM(routine), 'Input file type not recognized.')
    END IF

  END FUNCTION netcdf_read_real_2d

  FUNCTION netcdf_read_real_2d_time(input_file, variable_name, fill_array, &
    start_time_step, end_time_step)

    CLASS(t_input_file), INTENT(inout) :: input_file
    CHARACTER(LEN=*),   INTENT(in)    :: variable_name
    REAL(wp), TARGET, OPTIONAL        :: fill_array(:,:,:)
    INTEGER, OPTIONAL,  INTENT(in)    :: start_time_step, end_time_step
    REAL(wp), POINTER                 :: netcdf_read_real_2d_time(:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':netcdf_read_real_2d_time'

    netcdf_read_real_2d_time => netcdf_read_real_2d_extdim( &
      input_file, variable_name, fill_array,    &
      start_extdim=start_time_step, end_extdim=end_time_step, extdim_name="time")

  END FUNCTION netcdf_read_real_2d_time

  FUNCTION netcdf_read_real_2d_extdim(input_file, variable_name, fill_array, &
    start_extdim, end_extdim, extdim_name)

    CLASS(t_input_file),         INTENT(inout) :: input_file
    CHARACTER(LEN=*),           INTENT(in)    :: variable_name
    REAL(wp), TARGET, OPTIONAL                :: fill_array(:,:,:)
    INTEGER, OPTIONAL,          INTENT(in)    :: start_extdim, end_extdim
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)    :: extdim_name
    REAL(wp), POINTER                         :: netcdf_read_real_2d_extdim(:,:,:)

    INTEGER :: i, start_allocated_step, end_allocated_step, steps, start_step, end_step
    REAL(wp), POINTER :: tmp_array(:,:,:), ptr(:,:)
    INTEGER :: file_id, var_id, var_ndims, var_dim_ids(1:3), var_extdim_size, start_read_index(1:3), count_read_index(1:3)
    CHARACTER(LEN=NF_MAX_NAME) :: var_extdim_name

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':netcdf_read_real_2d_extdim'

    NULLIFY(netcdf_read_real_2d_extdim)

    IF (input_file%type == 1) THEN
#ifdef __ICON__
      CALL read_2D_extdim(input_file%stream_id, on_cells, TRIM(variable_name), fill_array, netcdf_read_real_2d_extdim, &
                          start_extdim=start_extdim, end_extdim=end_extdim, extdim_name=extdim_name)
#else
      CALL finish(TRIM(routine), 'Incompatible input file type')
#endif
    ELSE IF (input_file%type == 2) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'Incompatible input file type')
#else
      file_id = input_file%file_id

      IF (my_process_is_stdio()) THEN
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
        ALLOCATE(netcdf_read_real_2d_extdim(input_file%nproma, input_file%ngpblks, steps))
        netcdf_read_real_2d_extdim(:,:,:) = 0._wp
      END IF

      start_allocated_step = LBOUND(netcdf_read_real_2d_extdim, 3)
      end_allocated_step = UBOUND(netcdf_read_real_2d_extdim, 3)
      ALLOCATE(tmp_array(input_file%nlon, input_file%nlat, start_allocated_step:end_allocated_step))
      IF (SIZE(netcdf_read_real_2d_extdim,3) < steps) CALL finish(TRIM(routine), 'allocated size < steps')

      IF (my_process_is_stdio()) THEN
        start_read_index = (/ 1, 1, start_step /)
        count_read_index = (/ input_file%nlon, input_file%nlat, steps /)
        CALL nf(nf_get_vara_double(file_id, var_id, start_read_index, count_read_index, tmp_array), &
          TRIM(routine)//':'//TRIM(variable_name))
      END IF

      DO i=start_allocated_step,end_allocated_step
        IF (my_process_is_stdio()) ptr => tmp_array(:,:,i)
        CALL scatter_field(ptr, netcdf_read_real_2d_extdim(:,:,i))
      END DO

      DEALLOCATE(tmp_array)
#endif
    ELSE
      CALL finish(TRIM(routine), 'Input file type not recognized.')
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
