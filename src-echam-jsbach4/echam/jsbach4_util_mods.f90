!>
!! Contains modules that provide interfaces to ECHAM infrastructure for JSBACH4
!!
!! @par Revision History
!! Moved from adapters directory of JSBACH4 to ECHAM     by Reiner Schnur (2019-01-30)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
#ifdef JSBACH4

!! ==============================================================================================================================
!! @brief Contains interface to parallel infrastructure
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_parallel_iface

  USE mo_mpi,             ONLY: p_all_comm, p_global_comm
  USE mo_mpi, ONLY: p_io,    & !< Processor performing I/O
                    p_bcast    !< Broadcasting routine

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_parallel_iface, init_mpi_communicators, my_process_is_mpi_parallel, my_process_is_stdio, get_my_global_mpi_id
  PUBLIC :: p_io, mpi_comm, p_global_comm, p_bcast

  INTEGER :: mpi_comm

CONTAINS

  SUBROUTINE init_parallel_iface

    mpi_comm = p_all_comm

  END SUBROUTINE init_parallel_iface

!>
!! Only called for jsbach offline
  SUBROUTINE init_mpi_communicators

    USE mo_mpi, ONLY: p_init_communicators
    USE mo_control, ONLY: nproca, nprocb, nprocio

    CALL p_init_communicators(nproca, nprocb, nprocio)

  END SUBROUTINE init_mpi_communicators

  LOGICAL FUNCTION my_process_is_stdio()

  USE mo_mpi, ONLY: p_parallel_io

#ifdef NOMPI
    my_process_is_stdio = .TRUE.
#else
    my_process_is_stdio = p_parallel_io
#endif

  END FUNCTION my_process_is_stdio

  LOGICAL FUNCTION my_process_is_mpi_parallel()

    USE mo_mpi, ONLY: p_parallel

    my_process_is_mpi_parallel = p_parallel

  END FUNCTION my_process_is_mpi_parallel

  INTEGER FUNCTION get_my_global_mpi_id()

    USE mo_mpi, ONLY: p_pe

    get_my_global_mpi_id = p_pe

  END FUNCTION get_my_global_mpi_id

END MODULE mo_jsb_parallel_iface

!! ==============================================================================================================================
!>
!! @brief Contains interface to ECHAM decomposition for JSBACH4
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                              by Reiner Schnur (2013-04-04)
!! New subroutines and some adaptations for jsb4 standalone   by Julia Nabel (2016-02-20)
!!
MODULE mo_jsb_domain_iface

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
  USE mo_decomposition,      ONLY: t_patch => pe_decomposed, global_decomposition

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_patch, get_host_patch_id
  PUBLIC :: get_grid_filename, get_ntotal, get_ntotal_g, get_dims_g, get_nlat_g, get_nproma, get_nblks
  PUBLIC :: get_blk_start, get_blk_end, get_col_start, get_col_end
  PUBLIC :: jsbalone_init_decomp

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_domain_iface'

CONTAINS

  INTEGER FUNCTION get_host_patch_id(decomp)

    TYPE(t_patch), INTENT(in)   :: decomp

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_host_patch_id'

    get_host_patch_id = 1

  END FUNCTION get_host_patch_id

  FUNCTION get_grid_filename(decomp) RESULT(filename)

    USE mo_io_units, ONLY: filename_max

    TYPE(t_patch), INTENT(in)   :: decomp
    CHARACTER(len=filename_max) :: filename

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_grid_filename'

    filename = ''

  END FUNCTION get_grid_filename

  INTEGER FUNCTION get_ntotal(decomp)

    TYPE(t_patch), INTENT(in) :: decomp

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_ntotal'

    get_ntotal = decomp%nglon * decomp%nglat

  END FUNCTION get_ntotal

  INTEGER FUNCTION get_ntotal_g(decomp)

    TYPE(t_patch), INTENT(in) :: decomp

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_ntotal_g'

    get_ntotal_g = decomp%nlon * decomp%nlat

  END FUNCTION get_ntotal_g

  FUNCTION get_dims_g(decomp)

    TYPE(t_patch), INTENT(in) :: decomp
    INTEGER                                   :: get_dims_g(2)

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_dims_g'

    get_dims_g(1) = decomp%nlon
    get_dims_g(2) = decomp%nlat

  END FUNCTION get_dims_g

  FUNCTION get_nlat_g(decomp)

    TYPE(t_patch), INTENT(in) :: decomp
    INTEGER                   :: get_nlat_g

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_nlat_g'

    get_nlat_g = decomp%nlat

  END FUNCTION get_nlat_g

  INTEGER FUNCTION get_nblks(decomp)

    TYPE(t_patch), INTENT(in) :: decomp

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_nblks'

    get_nblks = decomp%ngpblks

  END FUNCTION get_nblks

  INTEGER FUNCTION get_nproma(decomp)

    TYPE(t_patch), INTENT(in) :: decomp

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_nproma'

    get_nproma = decomp%nproma

  END FUNCTION get_nproma

  INTEGER FUNCTION get_blk_start(decomp)

    TYPE(t_patch), INTENT(in) :: decomp

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_blk_start'

    get_blk_start = 1

  END FUNCTION get_blk_start

  INTEGER FUNCTION get_blk_end(decomp)

    TYPE(t_patch), INTENT(in) :: decomp

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_blk_end'

    get_blk_end = decomp%ngpblks

  END FUNCTION get_blk_end

  INTEGER FUNCTION get_col_start(jb, decomp)

    INTEGER,       INTENT(in) :: jb
    TYPE(t_patch), INTENT(in) :: decomp

    INTEGER :: i_startblk, i_endblk, jcs, jce, rl_start, rl_end

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_col_start'

    get_col_start = 1

  END FUNCTION get_col_start

  INTEGER FUNCTION get_col_end(jb, decomp)

    INTEGER,       INTENT(in) :: jb
    TYPE(t_patch), INTENT(in) :: decomp

    INTEGER :: i_startblk, i_endblk, jcs, jce, rl_start, rl_end

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_col_end'

    IF (jb == decomp%ngpblks) THEN
      get_col_end = decomp%npromz
    ELSE
      get_col_end = decomp%nproma
   END IF

  END FUNCTION get_col_end

  !>
  !! @brief Get a decomposition for jsbach standalone
  !!
  !! Currently only working for ECHAM
  !!
  !! @par History
  !! JN: Initial version copied from jsbach3 jsbalone_init_decomp.f90
  !!
  SUBROUTINE jsbalone_init_decomp (nlon, nlat, nproca, nprocb, npedim, decomposition)

    USE mo_decomposition, ONLY : decompose, local_decomposition
    USE mo_mpi,           ONLY : p_nprocs, p_parallel_io, p_pe
    USE mo_transpose,     ONLY : indx
    USE mo_io_units,      ONLY : nout

    INTEGER, INTENT(in) :: nlon, nlat, nproca, nprocb, npedim
    TYPE(t_patch), INTENT(OUT) :: decomposition
    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbalone_init_decomp'

    INTEGER :: debug_parallel, i, p
    LOGICAL :: lrot

    lrot          = .TRUE.  ! true: no rotation of longitudes

    ! nprocb has to be 1, otherwise the consistency checks in
    ! mo_decomposition will fail (nprocb - must be < nlev).
    IF (nprocb > 1) THEN
      CALL finish(TRIM(routine), 'nprocb has to be 1')
    ENDIF

    IF (p_nprocs == nproca*nprocb) THEN
      debug_parallel = -1
    ELSE IF (p_nprocs == nproca*nprocb+1) THEN
      debug_parallel = 0
    ELSE
      CALL finish(TRIM(routine), 'Number of runtime PEs doesn''t fit nproca*nprocb(+1)')
    ENDIF

    ! Derive decomposition
    ALLOCATE(global_decomposition(1:p_nprocs))

    ! The truncation is set to nproca, otherwise the consistency checks in
    ! mo_decomposition will fail
    CALL decompose(global_decomposition, npedim, nproca, nprocb, nlat, nlon, norot=lrot, nlev=1, &
                   nm=nproca, nn=nproca, nk=nproca, nadvec=0, debug=debug_parallel)

    ! Keep index values, not id values
    DO p=1,p_nprocs
      global_decomposition(p)%mapmesh(:,:) = &
            indx(global_decomposition(p)%mapmesh(:,:),global_decomposition)
    ENDDO

    ! Copy global decomposition table entry to local decomposition
    DO i=1,p_nprocs
       IF (global_decomposition(i)% pe == p_pe) THEN
          local_decomposition = global_decomposition(i)
       END IF
    END DO

    IF (p_parallel_io) THEN
       WRITE(nout,*) '---------------------------------------------'
       WRITE(nout,*) ' Blocking information:'
       WRITE(nout,*) '   lreg  nproma ngpblks  npromz'
       WRITE(nout,'(L5,3I8)') local_decomposition% lreg, &
            local_decomposition%nproma, &
            local_decomposition%ngpblks, &
            local_decomposition%npromz
       WRITE(nout,*) '---------------------------------------------'
    ENDIF

    ! Copy content of local_decomposition to output
    decomposition = local_decomposition
  END SUBROUTINE jsbalone_init_decomp

END MODULE mo_jsb_domain_iface

!! ==============================================================================================================================
!>
!! @brief Contains interface to namelist handling
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-02-06)
!!
MODULE mo_jsb_namelist_iface

  USE mo_exception, ONLY: finish

  USE mo_namelist, ONLY: open_nml, close_nml, position_nml_echam => position_nml
  USE mo_namelist, ONLY: POSITIONED

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: open_nml, POSITIONED, position_nml, close_nml

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_namelist_iface'

CONTAINS

  !>
  !! Wrapper function to make ICON subroutine position_nml a function, as in ECHAM.
  !! Should be replaced when open_nml in ICON gets updated.
  INTEGER FUNCTION position_nml(name, nml_handler, lrewind, status)

    CHARACTER(len=*), INTENT(in)            :: name         ! namelist group name
    INTEGER,          INTENT(in)            :: nml_handler  ! file handler
    LOGICAL,          INTENT(in)  ,OPTIONAL :: lrewind      ! default: true
    INTEGER,          INTENT(out) ,OPTIONAL :: status       ! error return value

    CHARACTER(len=*), PARAMETER :: routine = modname//':position_nml'

    position_nml = position_nml_echam(name, nml_handler, rewind=lrewind, status=status)

  END FUNCTION position_nml

END MODULE mo_jsb_namelist_iface

!! ==============================================================================================================================
!! @brief Contains the interface to io
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_io_netcdf_iface

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message
  USE mo_io_units,           ONLY: filename_max
  USE mo_jsb_domain_iface,   ONLY: t_patch
  USE mo_jsb_parallel_iface, ONLY: my_process_is_mpi_parallel, my_process_is_stdio, p_bcast, p_io, mpi_comm
  USE mo_tr_scatter, ONLY: scatter_field

  IMPLICIT NONE

  PUBLIC

  TYPE t_input_file
    INTEGER           :: TYPE                ! 1: ICON infrastructure, 2: ECHAM infrastructure
    CHARACTER(len=filename_max) :: filename  ! Name of input file
    LOGICAL           :: is_open = .FALSE.
    INTEGER :: file_id          ! netcdf file id
    INTEGER :: return_status    ! the latest operation return status
    INTEGER :: nlon    = -1
    INTEGER :: nlat    = -1
    INTEGER :: nproma  = -1
    INTEGER :: ngpblks = -1
  CONTAINS
    PROCEDURE :: Close          => netcdf_close_file
    PROCEDURE :: Read_1d        => netcdf_read_real_1d
    PROCEDURE :: Read_2d        => netcdf_read_real_2d
    PROCEDURE :: Read_2d_time   => netcdf_read_real_2d_time
    PROCEDURE :: Read_2d_1lev_1time => netcdf_read_real_2d_1lev_1time
    PROCEDURE :: Read_2d_extdim => netcdf_read_real_2d_extdim
    PROCEDURE :: Has_dim        => netcdf_file_has_dim
    PROCEDURE :: Has_var        => netcdf_file_has_var
    PROCEDURE :: Get_dimlen     => netcdf_file_get_dimlen
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

  INCLUDE 'netcdf.inc'

  INTEGER, PARAMETER :: MAX_VAR_DIMS = 16 ! NF_MAX_VAR_DIMS

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io_netcdf_iface'

CONTAINS

  TYPE(t_input_file) FUNCTION netcdf_open_input_generic(filename)

    CHARACTER(LEN=*),    INTENT(in) :: filename

    INTEGER :: file_id, return_status

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_open_input_generic'

    netcdf_open_input_generic%filename = TRIM(filename)

    netcdf_open_input_generic%type = 2

    IF (my_process_is_stdio()) THEN
      return_status = nf_open(TRIM(filename), NF_NOWRITE, file_id)
      CALL nf(return_status, TRIM(routine)//':'//TRIM(filename))
      netcdf_open_input_generic%is_open = .TRUE.
    ELSE
      return_status = 1
      file_id = -1
    END IF
    netcdf_open_input_generic%file_id       = file_id
    netcdf_open_input_generic%return_status = return_status

  END FUNCTION netcdf_open_input_generic

  TYPE(t_input_file) FUNCTION netcdf_open_input_patch(filename, decomp)

    CHARACTER(LEN=*),  INTENT(in) :: filename
    TYPE(t_patch),     INTENT(in) :: decomp

    INTEGER :: file_id, return_status

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_open_input_patch'

    netcdf_open_input_patch%filename = TRIM(filename)

    IF (decomp%nlon > 0) THEN
      netcdf_open_input_patch%type = 2

      IF (my_process_is_stdio()) THEN
        return_status = nf_open(TRIM(filename), NF_NOWRITE, file_id)
        CALL nf(return_status, TRIM(routine)//':'//TRIM(filename))
        netcdf_open_input_patch%is_open = .TRUE.
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
    ELSE
      CALL finish(TRIM(routine), 'No decomposition')
    END IF

  END FUNCTION netcdf_open_input_patch

  SUBROUTINE netcdf_close_file(input_file)

    CLASS(t_input_file), INTENT(inout) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_close_file'

    IF (my_process_is_stdio()) &
      CALL nf(nf_close(input_file%file_id), TRIM(routine)//':'//TRIM(input_file%filename))

      input_file%is_open = .FALSE.

  END SUBROUTINE netcdf_close_file

  FUNCTION netcdf_read_real_1d(input_file, variable_name, fill_array)

    CLASS(t_input_file), INTENT(inout) :: input_file
    CHARACTER(LEN=*),   INTENT(in)     :: variable_name
    REAL(wp), TARGET, OPTIONAL         :: fill_array(:)
    REAL(wp), POINTER                  :: netcdf_read_real_1d(:)

    REAL(wp), POINTER :: tmp_array(:)
    INTEGER :: file_id, var_id, dim_id(MAX_VAR_DIMS), ndims
    INTEGER :: varid, var_type
    INTEGER :: var_size(MAX_VAR_DIMS)

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_read_real_1d'

    NULLIFY(netcdf_read_real_1d)

    IF (input_file%type == 1) THEN
      CALL finish(TRIM(routine), 'Incompatible input file type')
    ELSE IF (input_file%type == 2) THEN
      file_id = input_file%file_id

      IF (my_process_is_stdio()) THEN
        CALL nf(nf_inq_varid(file_id, TRIM(variable_name), varid), TRIM(routine))
        CALL nf(nf_inq_varndims(file_id, varid, ndims), TRIM(routine))
        IF (ndims /= 1) CALL finish(TRIM(routine), 'only one variable dimension expected')
        CALL nf(nf_inq_vardimid(file_id, varid, dim_id(1:ndims)), TRIM(routine))
        CALL nf(nf_inq_dimlen(file_id, dim_id(1), var_size(1)), TRIM(routine))
      END IF

      CALL p_bcast(var_size(1), p_io, mpi_comm)

      ALLOCATE(tmp_array(var_size(1)))

      IF (my_process_is_stdio()) THEN
        CALL nf(nf_inq_varid(file_id, TRIM(variable_name), var_id), TRIM(routine))

        CALL nf(nf_get_var_double(file_id, var_id, tmp_array), TRIM(routine))
      END IF
      CALL p_bcast(tmp_array, p_io, mpi_comm)

      IF (PRESENT(fill_array)) THEN
        fill_array(:) = tmp_array(:)
        netcdf_read_real_1d => fill_array
      ELSE
        ALLOCATE(netcdf_read_real_1d(var_size(1)))
        netcdf_read_real_1d(:) = tmp_array(:)
      END IF

      DEALLOCATE(tmp_array)
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
    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_read_real_2d'

    NULLIFY(netcdf_read_real_2d)

    IF (input_file%type == 1) THEN
      CALL finish(TRIM(routine), 'Incompatible input file type')
    ELSE IF (input_file%type == 2) THEN
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

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_read_real_2d_time'

    netcdf_read_real_2d_time => netcdf_read_real_2d_extdim( &
      input_file, variable_name, fill_array,    &
      start_extdim=start_time_step, end_extdim=end_time_step, extdim_name="time")

  END FUNCTION netcdf_read_real_2d_time

  ! Dummy routine, not needed with ECHAM
  FUNCTION netcdf_read_real_2d_1lev_1time(input_file, variable_name, fill_array)

    CLASS(t_input_file), INTENT(inout) :: input_file
    CHARACTER(LEN=*),   INTENT(in)    :: variable_name
    REAL(wp), TARGET, OPTIONAL        :: fill_array(:,:)
    REAL(wp), POINTER                 :: netcdf_read_real_2d_1lev_1time(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_read_real_2d_1lev_1time'

    netcdf_read_real_2d_1lev_1time(:,:) = 0._wp

    CALL finish(TRIM(routine), 'This function should not be called with ECHAM (only with ICON and init_from_ifs=.TRUE.')

  END FUNCTION netcdf_read_real_2d_1lev_1time

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
    INTEGER :: file_id, var_id
    INTEGER :: var_ndims, var_dim_ids(1:3), var_extdim_size, start_read_index(1:3), count_read_index(1:3)
    CHARACTER(LEN=NF_MAX_NAME) :: var_extdim_name

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_read_real_2d_extdim'

    NULLIFY(netcdf_read_real_2d_extdim)

    IF (input_file%type == 1) THEN
      CALL finish(TRIM(routine), 'Incompatible input file type')
    ELSE IF (input_file%type == 2) THEN
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

      CALL p_bcast(var_extdim_size, p_io, mpi_comm)

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
    ELSE
      CALL finish(TRIM(routine), 'Input file type not recognized.')
    END IF

  END FUNCTION netcdf_read_real_2d_extdim

  LOGICAL FUNCTION netcdf_file_has_dim(input_file, dimname)

    CLASS(t_input_file), INTENT(in) :: input_file
    CHARACTER(LEN=*),    INTENT(in) :: dimname

    INTEGER :: dimid, status

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_file_has_dim'

    IF (my_process_is_stdio()) THEN
      IF (.NOT. input_file%is_open) CALL finish(TRIM(routine), 'NetCDF file not open')

      status = nf_inq_dimid(input_file%file_id, TRIM(dimname), dimid)
      netcdf_file_has_dim = (status == NF_NOERR)
    END IF

    CALL p_bcast(netcdf_file_has_dim, p_io, mpi_comm)

  END FUNCTION netcdf_file_has_dim

  LOGICAL FUNCTION netcdf_file_has_var(input_file, varname)

    CLASS(t_input_file), INTENT(in) :: input_file
    CHARACTER(LEN=*),    INTENT(in) :: varname

    INTEGER :: varid, status

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_file_has_var'

    IF (my_process_is_stdio()) THEN
      IF(.NOT. input_file%is_open) CALL finish(TRIM(routine), 'NetCDF file not open')

      status = nf_inq_varid(input_file%file_id, TRIM(varname), varid)
      netcdf_file_has_var = (status == NF_NOERR)
    END IF

    CALL p_bcast(netcdf_file_has_var, p_io, mpi_comm)

  END FUNCTION netcdf_file_has_var

  INTEGER FUNCTION netcdf_file_get_dimlen(input_file, dimname)

    CLASS(t_input_file), INTENT(in) :: input_file
    CHARACTER(LEN=*),    INTENT(in) :: dimname

    INTEGER :: IO_dim_id, status

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_file_get_dimlen'

    IF (my_process_is_stdio()) THEN
      IF(.NOT. input_file%is_open) CALL finish(TRIM(routine), 'NetCDF file not open')

      CALL nf(nf_inq_dimid  (input_file%file_id, TRIM(dimname), IO_dim_id), TRIM(routine))
      CALL nf(nf_inq_dimlen (input_file%file_id, IO_dim_id, netcdf_file_get_dimlen), TRIM(routine))
    ELSE
      CALL finish(TRIM(routine), 'should only be called on io process!')
    END IF

  END FUNCTION netcdf_file_get_dimlen

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

END MODULE mo_jsb_io_netcdf_iface

!! ==============================================================================================================================
!>
!>
!! @brief Contains interface to time control
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                              by Reiner Schnur (2013-02-06)
!! New subroutines and some adaptations for jsb4 standalone   by Julia Nabel (2016-02-20)
!! Added ltimer although this doesn't really belong here      by Reiner Schnur (2018-06-04)
!!
MODULE mo_jsb_time_iface

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message, message_text

  USE mo_time_conversion, ONLY: t_datetime => time_days, &
                                get_day_length => day_len
  USE mo_time_control,    ONLY: get_year_length => get_year_len, get_month_length => get_month_len, &
                                get_year_day_echam => get_year_day, get_date_components
  USE mo_time_control,    ONLY: start_date, resume_date, stop_date, current_date, previous_date, next_date, &
    &                           lstart, lresume, l_trigradm1
  USE mo_interpo,         ONLY: wgt1, wgt2, nmw1, nmw2
  USE mo_semi_impl,       ONLY: eps
  USE mo_control,         ONLY: l_timer_host => ltimer

#ifdef STANDALONE
  USE mo_jsb_parallel_iface, ONLY: my_process_is_stdio, p_io, p_bcast, my_process_is_mpi_parallel, mpi_comm
  USE mo_time_event,      ONLY: io_time_event, TIME_INC_DAYS, TIME_INC_MONTHS, TRIG_FIRST, TRIG_LAST
  USE mo_time_control,    ONLY: lresume, dt_start, dt_resume, dt_stop, delta_time, no_cycles, &
                                no_days, no_steps, putdata, putrerun, trigfiles, p_bcast_event

  USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, close_nml, position_nml
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_datetime, deallocateDatetime
  PUBLIC :: get_year_length, get_month_length, get_day_length, get_year_day
  PUBLIC :: get_time_dt, get_time_nsteps, &
            get_time_experiment_start, get_time_start, get_time_stop, &
            get_time_current, get_time_previous, get_time_next, is_time_experiment_start, is_time_restart, &
            is_time_ltrig_rad_m1
  PUBLIC :: read_time_namelist !, configure_time
  PUBLIC :: start_timestep, finish_timestep
  PUBLIC :: configure_time_and_events
  PUBLIC :: get_date_components
  PUBLIC :: get_time_interpolation_weights
  PUBLIC :: get_asselin_coef
  PUBLIC :: l_timer_host

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_time_iface'

CONTAINS

  SUBROUTINE read_time_namelist(filename, restart_jsbach)

    CHARACTER(LEN=*), INTENT(in) :: filename
    LOGICAL, INTENT(IN) :: restart_jsbach

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_time_namelist'

#ifdef STANDALONE
    INTEGER :: nml_handler, nml_unit, istat

    NAMELIST /jsb_time_nml/ &
      no_cycles,            &
      dt_start,             &
      dt_resume,            &
      dt_stop,              &
      no_days,              &
      no_steps,             &
      delta_time,           &
      putdata,              &
      putrerun,             &
      trigfiles

    !-----------------------------------------------------------------------------------------------

    !NOTE: lresume is not read here, but was already read as 'restart_jsbach' from the jsb_control_nml
    ! in read_jsb_control_namelist. However, it needs to be set here to properly work with echam
    IF(restart_jsbach) THEN
      lresume = .TRUE.
    ELSE
      lresume = .FALSE.
    ENDIF

    IF(my_process_is_stdio()) THEN
      ! default values
      no_cycles = 1
      dt_start(:) = 0
      dt_resume(:) = 0
      dt_stop(:) = 0
      no_days = -1
      no_steps = -1
      delta_time = 0
      putdata  = io_time_event(1, TIME_INC_DAYS, TRIG_FIRST, 0)
      putrerun = io_time_event(1, TIME_INC_MONTHS, TRIG_LAST, 0)
      trigfiles = io_time_event(1, TIME_INC_MONTHS, TRIG_FIRST, 0)

      nml_handler = open_nml(TRIM(filename))

      nml_unit =  position_nml('jsb_time_nml', nml_handler, STATUS=istat)
      IF (istat == POSITIONED) READ(nml_unit, jsb_time_nml)

      CALL close_nml(nml_handler)

      !Assertion: in case of lresume = .TRUE. dt_start should not be specified, because this
      !   can lead to problemes (diverging results) if it does not match what echam calculates from
      !   fdate/vdate i.e. do not specifiy dt_start if running with lresume but let echam calculate it!
      IF ((lresume) .AND. ANY(dt_start(:) .NE. 0)) CALL finish(routine, 'If lresume=true dt_start should not be specified' &
        &//' (it seems to lead to inconsistencies in the echam event management).')

    ENDIF

    IF (my_process_is_mpi_parallel()) THEN
    ! jsb3 text: The following are used from ECHAM5 modules and ... need
    ! to be broadcast if running in standalone mode so that the modules can use them
      CALL p_bcast(no_cycles, p_io, comm=mpi_comm)
      CALL p_bcast(dt_start, p_io, comm=mpi_comm)
      CALL p_bcast(dt_resume, p_io, comm=mpi_comm)
      CALL p_bcast(dt_stop, p_io, comm=mpi_comm)
      CALL p_bcast(no_days, p_io, comm=mpi_comm)
      CALL p_bcast(no_steps, p_io, comm=mpi_comm)
      CALL p_bcast(delta_time, p_io, comm=mpi_comm)
      CALL p_bcast_event(putdata, p_io, comm=mpi_comm)
      CALL p_bcast_event(trigfiles, p_io, comm=mpi_comm)
      CALL p_bcast_event(putrerun, p_io, comm=mpi_comm)
    ENDIF
#endif

  END SUBROUTINE read_time_namelist

  SUBROUTINE configure_time_and_events(model_shortname)

#ifdef STANDALONE
    USE mo_time_control, ONLY: lfirst_cycle, ec_manager_init, &
                               init_manager, init_times, init_events, time_set
#endif
    USE mo_time_control, ONLY: time_weights

    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: model_shortname

    CHARACTER(len=*), PARAMETER :: routine = modname//':configure_time_and_events'
    INTEGER :: IO_timestep, istep

#ifdef STANDALONE
    !Assertion:
    IF (.NOT. PRESENT(model_shortname)) CALL finish(routine, &
      & 'Short name missing! In standalone model, a model_short name needs to be passed, please check!')

    CALL get_dates(model_shortname, IO_timestep, istep)           ! get timestep
    CALL ec_manager_init(IO_timestep, istep)     ! time manager initialization
    IF (lfirst_cycle) CALL init_manager
    CALL init_times                              ! initialize all dates
    CALL init_events
#endif
    CALL time_weights(l_compute_rad=.FALSE.) ! radiation_time is not initialized yet at the time of JSBACH initialization


  END SUBROUTINE configure_time_and_events

  !>
  !! Get time step length
  !!
  !! Function returns the length of the integration time step
  !!
  REAL(wp) FUNCTION get_time_dt(model_id)

    USE mo_time_control, ONLY: delta_time

    INTEGER, INTENT(in) :: model_id

    REAL(wp) :: ztime

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_time_dt'

    IF (model_id > 1) CALL finish(routine, 'Only one JSBACH model instance allowed')

    ztime = delta_time

    IF (ztime <= 0.0_wp) &
      CALL finish(routine, 'time step not configured yet.')

    get_time_dt = ztime

  END FUNCTION get_time_dt

  !>
  !! Get number of time steps
  !!
  !! Function returns the number of time steps for current model integration
  !!
  INTEGER FUNCTION get_time_nsteps()

    USE mo_time_control, ONLY: start_date, stop_date
    USE mo_time_base,    ONLY: IDAYLEN
    USE mo_time_conversion, ONLY: TC_get

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_time_nsteps'

    INTEGER :: start_day, start_second, stop_day, stop_second, idt
    INTEGER :: isteps

    CALL TC_get(start_date, start_day, start_second)
    CALL TC_get(stop_date, stop_day, stop_second)

    idt = INT(get_time_dt(1))
    isteps =                                   &
      & (( IDAYLEN/idt * (stop_day - start_day - 1))   & ! Number of seconds in full days
      &    + (stop_second                              & ! Number of seconds in last day
      &       + (IDAYLEN - start_second - 1)           & ! Number of seconds in first day
      &      ) / idt)

!@todo: error message is misleading! should be: error in istep calculation or comparable!?
    IF (isteps <= 0) &
      CALL finish(routine, 'number of time steps not configured yet.')

    get_time_nsteps = isteps

  END FUNCTION get_time_nsteps

  FUNCTION get_time_experiment_start()

    TYPE(t_datetime), POINTER :: get_time_experiment_start

    ALLOCATE(get_time_experiment_start, source=start_date)

  END FUNCTION get_time_experiment_start

  FUNCTION get_time_start()

    TYPE(t_datetime), POINTER :: get_time_start

    ALLOCATE(get_time_start, source=resume_date)

  END FUNCTION get_time_start

  FUNCTION get_time_stop()

    TYPE(t_datetime), POINTER :: get_time_stop

    ALLOCATE(get_time_stop, source=stop_date)

  END FUNCTION get_time_stop

  !>
  !! @brief Get current time
  !!
  FUNCTION get_time_current() RESULT(current)

    TYPE(t_datetime), POINTER :: current

    ALLOCATE(current, source=current_date)

  END FUNCTION get_time_current

  !>
  !! @brief Get time at previous time step
  !!
  FUNCTION get_time_previous(current, dt) RESULT(previous)

    TYPE(t_datetime), POINTER, INTENT(in) :: current
    REAL(wp),                  INTENT(in) :: dt
    TYPE(t_datetime), POINTER             :: previous

    ALLOCATE(previous, source=previous_date)

  END FUNCTION get_time_previous

  !>
  !! @brief Get time at next time step
  !!
  FUNCTION get_time_next(current, dt) RESULT(next)

    TYPE(t_datetime), POINTER, INTENT(in) :: current
    REAL(wp),                  INTENT(in) :: dt
    TYPE(t_datetime), POINTER             :: next

    ALLOCATE(next, source=next_date)

  END FUNCTION get_time_next

  !
  ! Check whether we are one time step before radiation is calculated in atmosphere
  !
  FUNCTION is_time_ltrig_rad_m1(current, dt, model_id) RESULT(ltrig_rad_m1)

    TYPE(t_datetime), POINTER, INTENT(in) :: current
    REAL(wp),                  INTENT(in) :: dt
    INTEGER,                   INTENT(in) :: model_id
    LOGICAL                               :: ltrig_rad_m1

    ltrig_rad_m1 = l_trigradm1

  END FUNCTION is_time_ltrig_rad_m1

  LOGICAL FUNCTION is_time_experiment_start(current)

    TYPE(t_datetime), POINTER, INTENT(in) :: current

    ! In ECHAM, the first computed time step is the start date
    is_time_experiment_start = lstart

  END FUNCTION is_time_experiment_start

  LOGICAL FUNCTION is_time_restart(current)

  TYPE(t_datetime), POINTER, INTENT(in) :: current

  is_time_restart = lresume

  END FUNCTION is_time_restart

  SUBROUTINE deallocateDatetime(datetime_ptr)

    TYPE(t_datetime), POINTER, INTENT(inout) :: datetime_ptr

    DEALLOCATE(datetime_ptr)

  END SUBROUTINE deallocateDatetime

  SUBROUTINE start_timestep(istep)

    USE mo_time_control, ONLY: time_set

    INTEGER, INTENT(in) :: istep

    CHARACTER(len=*), PARAMETER :: routine = modname//':start_timestep'

    CALL time_set

  END SUBROUTINE start_timestep

  SUBROUTINE finish_timestep(istep)

    USE mo_time_control, ONLY: time_reset

    INTEGER, INTENT(in) :: istep

    CHARACTER(len=*), PARAMETER :: routine = modname//':finish_timestep'

    CALL time_reset

  END SUBROUTINE finish_timestep
  !>
  !! Get day in year
  !!
  !! Function returns the day of the year for a given date.
  !! The seconds in the day are given as fraction of the day.
  !!
  REAL(wp) FUNCTION get_year_day(date)

    ! t_datetime points to mtimes datetime
    TYPE(t_datetime), POINTER, INTENT(in)  :: date

    TYPE(t_datetime) :: date_tmp

    date_tmp = date
    get_year_day = get_year_day_echam(date_tmp)

  END FUNCTION get_year_day

  SUBROUTINE  get_time_interpolation_weights(w1, w2, n1, n2)

    REAL(wp), INTENT(out) :: w1, w2
    INTEGER, INTENT(out)  :: n1, n2

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_time_interpolation_weights'

    w1 = wgt1
    w2 = wgt2
    n1 = nmw1
    n2 = nmw2

  END SUBROUTINE get_time_interpolation_weights

  REAL(wp) FUNCTION get_asselin_coef()

    get_asselin_coef = eps

  END FUNCTION get_asselin_coef

#ifdef STANDALONE
  !>
  !! @brief Get dates for jsbach standalone
  !!
  !! Currently only working for ECHAM
  !!
  !! @par History
  !! JN: Initial version copied from jsbach3 mo_jsbach_interface.f90
  !!
  SUBROUTINE get_dates(model_shortname, IO_timestep, istep)

    ! For standalone mode only.
    ! Corresponds to ECHAM5 subroutine IO_init in mo_io.f90
    ! This routine doesn't need to be called from ECHAM in coupled mode.

    USE mo_filename, ONLY: out_expname
    USE mo_jsb_io_netcdf_iface, ONLY: t_input_file, netcdf_open_input
    USE mo_netCDF, ONLY: NF_GLOBAL, io_get_att_int
    USE mo_time_control, ONLY: INIT_STEP, delta_time, dt_start, &
                               lresume, start_date, inp_convert_date, write_date
    USE mo_time_conversion, ONLY : time_native, TC_set, TC_convert

    USE mo_mpi,              ONLY : p_start, p_bcast, p_parallel_io, p_parallel, p_io

    CHARACTER(len=*), INTENT(IN) :: model_shortname
    INTEGER, INTENT(out) :: IO_timestep, istep

    TYPE(t_input_file) :: input_file
    INTEGER :: IO_file_id
    INTEGER :: forecast_date, verification_date   ! YYYYMMDD
    INTEGER :: forecast_time, verification_time   ! HHMMSS
    TYPE(time_native) :: date_nat

    ! Initialize time
    IF (p_parallel_io) THEN
       IF (lresume) THEN
         ! copied from jsb3 get_dates:
         ! Get time information from restart file
          input_file = netcdf_open_input('restart_'//TRIM(out_expname)//'_'//TRIM(model_shortname)//'.nc')
          IO_file_id = input_file%file_id
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'nstep', istep)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'timestep', IO_timestep)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'fdate', forecast_date)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'ftime', forecast_time)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'vdate', verification_date)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'vtime', verification_time)
          CALL input_file%Close()
       ELSE
          istep = INIT_STEP
          IO_timestep = delta_time
       ENDIF
    ENDIF

    IF (p_parallel) THEN
       CALL p_bcast(istep, p_io, comm=mpi_comm)
       CALL p_bcast(IO_timestep, p_io, comm=mpi_comm)
       IF (lresume) THEN
          CALL p_bcast(forecast_date, p_io, comm=mpi_comm)
          CALL p_bcast(forecast_time, p_io, comm=mpi_comm)
          CALL p_bcast(verification_date, p_io, comm=mpi_comm)
          CALL p_bcast(verification_time, p_io, comm=mpi_comm)
       ENDIF
    ENDIF

    IF (lresume) THEN
       CALL inp_convert_date(forecast_date, forecast_time, start_date)
       CALL inp_convert_date(verification_date, verification_time, resume_date)
    ELSE
       IF (SUM(dt_start(:)) /= 0) THEN
          CALL TC_set(dt_start(1),dt_start(2),dt_start(3), &
               dt_start(4),dt_start(5),dt_start(6), date_nat)
          CALL TC_convert(date_nat, start_date)
          resume_date = start_date
       ELSE
          CALL finish('mo_jsb_time_iface', 'Start date not set')
       ENDIF
    ENDIF

    IF (p_parallel_io) THEN
       CALL write_date(start_date, 'Start date (initial/restart): ')
       CALL write_date(resume_date, 'Resume date (initial/restart): ')
    ENDIF

  END SUBROUTINE get_dates

#endif

END MODULE mo_jsb_time_iface

!! ==============================================================================================================================
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

  USE mo_kind,               ONLY: wp
  USE mo_jsb_time_iface,     ONLY: is_time_experiment_start, t_datetime

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
  SUBROUTINE write_interface_variables( ics, ice, iblk, current_datetime, t_air, q_air, rain, snow, wind_air, &
    & wind_10m, lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, press_srf, drag_srf, &
    & t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle )

    USE mo_exception,     ONLY: finish, message_text

    USE mo_jsb_parallel_iface, ONLY: my_process_is_mpi_parallel

    USE mo_time_control, ONLY: lresume, get_date_components, next_date
    USE mo_filename,     ONLY: compose_filenames, standard_output_file
    USE mo_netcdf,       ONLY: nf_check

    USE mo_decomposition, ONLY: local_decomposition
    USE mo_gaussgrid,     ONLY: philon,philat

    INTEGER,  INTENT(in) :: ics, ice, iblk
    TYPE(t_datetime), POINTER, INTENT(in) :: current_datetime
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

    CHARACTER(len=*), PARAMETER :: routine = modname//':write_interface_variables'

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
    IF ((is_time_experiment_start(current_datetime) .OR. lresume) .AND. iblk == 1) THEN
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

  END SUBROUTINE write_interface_variables

  !>
  !! @brief read all variables required for a call to the full jsbach interface
  !!
  !! Currently only working for ECHAM
  !!
  !! @par History
  !! Initial version copied from jsbach3 mo_test.f90: read_interface_variables
  !!
  SUBROUTINE read_interface_variables( ics, ice, iblk, ncid, current_datetime, t_air, q_air, rain, snow, wind_air, wind_10m, &
             & lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, press_srf, drag_srf, &
             & t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle )

    USE mo_exception,     ONLY: finish, message_text

    USE mo_jsb_parallel_iface, ONLY: my_process_is_mpi_parallel

    USE mo_time_control, ONLY: lresume, get_date_components, next_date
    USE mo_netcdf,       ONLY: nf_check

    USE mo_decomposition, ONLY: local_decomposition

    INTEGER,  INTENT(in)  :: ics, ice, iblk
    INTEGER,  INTENT(inout) :: ncid
    TYPE(t_datetime), POINTER, INTENT(in) :: current_datetime
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

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_interface_variables'
    CHARACTER(len=*), PARAMETER :: filename = 'interface_variables.nc'

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

    IF ((is_time_experiment_start(current_datetime) .OR. lresume) .AND. iblk == 1) THEN
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

  END SUBROUTINE read_interface_variables

END MODULE mo_jsb_test

!! ==============================================================================================================================
!>
!! @brief Contains the interface to io
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                 by Reiner Schnur (2013-04-04)
!! IO namelist subroutine for jsb4 standalone    by Julia Nabel (2016-02-20)
!!
MODULE mo_jsb_io_iface

  USE mo_kind,               ONLY: wp
  USE mo_linked_list,   ONLY: NETCDF2, NETCDF4, GRIB
  USE mo_netcdf,        ONLY: add_dim
  USE mo_control,       ONLY: ldebugio

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_cf_var, t_grib2_var, grib2_var
  PUBLIC :: Get_netcdf_precision
  PUBLIC :: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, DATATYPE_PACK24,  &
            ZAXIS_SURFACE, ZAXIS_GENERIC, ZAXIS_DEPTH_BELOW_LAND,               &
            FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB, FILETYPE_GRB2, &
            GRID_CELL, GRID_UNSTRUCTURED, GRID_UNSTRUCTURED_CELL, GRID_LONLAT, &
            TSTEP_CONSTANT, TSTEP_INSTANT, &
            Create_zaxis, cdiDefMissval, read_jsb_io_namelist
            ! Create_zaxis, Destroy_zaxis, cdiDefMissval, read_jsb_io_namelist
  PUBLIC :: ldebugio

  TYPE t_cf_var
    CHARACTER(len=128) :: standard_name = ''
    CHARACTER(len=128) :: units         = ''
    CHARACTER(len=128) :: long_name     = ''
    INTEGER            :: datatype      = -1
  END TYPE t_cf_var

  TYPE t_grib2_var
    INTEGER :: discipline
    INTEGER :: category
    INTEGER :: number
    INTEGER :: bits
    INTEGER :: gridtype
    INTEGER :: subgridtype
  END TYPE t_grib2_var

  INCLUDE 'cdi.inc'

  ! From ICON mo_cdi_constants.f90
  INTEGER, PARAMETER :: GRID_CELL              =  1
  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_CELL =  1

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io_iface'

CONTAINS

  FUNCTION Get_netcdf_precision() RESULT(return_value)

    INTEGER :: return_value

    ! This might be overwritten by ECHAM in mo_output.f90 based on the value of grib precision.
    return_value = DATATYPE_FLT32

  END FUNCTION Get_netcdf_precision

  ! constructor for GRIB2 derived data type
  FUNCTION grib2_var(discipline, category, number, bits, gridtype, subgridtype)
    INTEGER, INTENT(IN) :: discipline
    INTEGER, INTENT(IN) :: category
    INTEGER, INTENT(IN) :: number
    INTEGER, INTENT(IN) :: bits
    INTEGER, INTENT(IN) :: gridtype
    INTEGER, INTENT(IN) :: subgridtype
    TYPE(t_grib2_var) :: grib2_var

    grib2_var%discipline  = discipline
    grib2_var%category    = category
    grib2_var%number      = number
    grib2_var%bits        = bits
    grib2_var%gridtype    = gridtype
    grib2_var%subgridtype = subgridtype

!!$    grib2_var%additional_keys%nint_keys = 0
  END FUNCTION grib2_var

  SUBROUTINE Create_zaxis(ZaxisID, echamZaxisIdx, cdi_axis_type, length, name, &
    & levels, units, longname, lbounds, ubounds)

    INTEGER,           INTENT(out)   :: ZaxisID
    INTEGER,           INTENT(out)   :: echamZaxisIdx
    INTEGER,           INTENT(in)    :: cdi_axis_type
    INTEGER,           INTENT(in)    :: length
    CHARACTER(len=*),  INTENT(in)    :: name
    REAL(wp),          INTENT(in)    :: levels(:)
    CHARACTER(len=*),  INTENT(in)    :: units
    CHARACTER(len=*),  INTENT(in)    :: longname
    REAL(wp),          INTENT(in)    :: lbounds(:)
    REAL(wp),          INTENT(in)    :: ubounds(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Create_zaxis'

    ZaxisID = zaxisCreate(cdi_axis_type, length)
    CALL zaxisDefName    (ZaxisID, TRIM(name))
    CALL zaxisDefLevels  (ZaxisID, levels(1:length))
    CALL zaxisDefLbounds (ZaxisID, lbounds(1:length))
    CALL zaxisDefUbounds (ZaxisID, ubounds(1:length))
    CALL zaxisDefUnits   (ZaxisID, TRIM(units))
    CALL zaxisDefLongname(ZaxisID, TRIM(longname))

    CALL add_dim(TRIM(name), length, &
      & longname=TRIM(longname), units=TRIM(units), VALUE=levels(1:length), &
      & indx=echamZaxisIdx)

  END SUBROUTINE Create_zaxis

!   SUBROUTINE Destroy_zaxis(ZaxisID)

!     INTEGER, INTENT(in) :: ZaxisID

!     ! In ICON, mo_vertical_axis_types:t_verticalAxisList_finalize takes care of destroying the vertical
!     ! axes and corresponding CDI zaxes. It is called in mo_name_list_output:close_name_list_output
! #ifndef __ICON__
!     CALL zaxisDestroy(ZaxisID)
! #endif

!   END SUBROUTINE Destroy_zaxis

  SUBROUTINE read_jsb_io_namelist(filename)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, close_nml, position_nml
    USE mo_jsb_parallel_iface, ONLY: my_process_is_stdio, p_io, p_bcast, my_process_is_mpi_parallel, mpi_comm
    USE mo_jsb_test,           ONLY: read_interface_vars, write_interface_vars
#ifdef STANDALONE
    USE mo_filename, ONLY: out_expname
#endif

    CHARACTER(LEN=*), INTENT(in) :: filename

    INTEGER :: nml_handler, nml_unit, istat

    NAMELIST /jsb_io_nml/ &
#ifdef STANDALONE
      out_expname,          &
      out_expname,          &
      read_interface_vars,  &
#endif
      write_interface_vars

    IF(my_process_is_stdio()) THEN
      ! default values
#ifdef STANDALONE
      out_expname = 'undefined_expname'
#endif
      read_interface_vars = .FALSE.
      write_interface_vars = .FALSE.

      nml_handler = open_nml(TRIM(filename))

      nml_unit =  position_nml('jsb_io_nml', nml_handler, STATUS=istat)
      IF (istat == POSITIONED) READ(nml_unit, jsb_io_nml)

      CALL close_nml(nml_handler)

    ENDIF

    IF (my_process_is_mpi_parallel()) THEN
#ifdef STANDALONE
      CALL p_bcast(out_expname, p_io, comm=mpi_comm)
#endif
      CALL p_bcast(read_interface_vars, p_io, comm=mpi_comm)
      CALL p_bcast(write_interface_vars, p_io, comm=mpi_comm)
    ENDIF

  END SUBROUTINE read_jsb_io_namelist

END MODULE mo_jsb_io_iface

!! ==============================================================================================================================
!>
!! @brief Contains grid utilities for JSBACH4
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                              by Reiner Schnur (2013-04-04)
!! New subroutines and some adaptations for jsb4 standalone   by Julia Nabel (2016-02-20)
!!
MODULE mo_jsb_grid_iface

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
  USE mo_jsb_domain_iface,   ONLY: t_patch, get_nproma, get_nblks
  USE mo_decomposition,      ONLY: global_decomposition

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_lon, get_lat, get_area, get_nlon_nlat_from_file, get_lons_and_lats_from_file

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_grid_iface'

CONTAINS

FUNCTION get_lon(decomp) RESULT(lon)

#ifdef STANDALONE
    USE mo_control,   ONLY: nlon
    USE mo_gaussgrid,  ONLY: philon
    USE mo_transpose,     ONLY: reorder
    USE mo_decomposition, ONLY: ldc => local_decomposition
#endif

    USE mo_geoloc,     ONLY: philon_2d

    TYPE(t_patch), INTENT(in)  :: decomp
    REAL(wp),      POINTER     :: lon(:,:)

    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_lon'

#ifdef STANDALONE
    INTEGER :: jlat, jlon
    REAL(wp) :: tmp (ldc%nglon, ldc%nglat)
#endif

    nproma = get_nproma(decomp)
    nblks  = get_nblks (decomp)
    ALLOCATE(lon(nproma,nblks))

    IF (.NOT. ASSOCIATED(philon_2d)) &
      CALL finish(TRIM(routine), 'Longitudes of decomposition not computed yet')
    lon(:,:) = philon_2d(:,:)

  END FUNCTION get_lon

  FUNCTION get_lat(decomp) RESULT(lat)

    USE mo_geoloc, ONLY: philat_2d

    TYPE(t_patch),       INTENT(in)  :: decomp
    REAL(wp), POINTER                :: lat(:,:)

    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_lat'

    nproma = get_nproma(decomp)
    nblks  = get_nblks (decomp)
    ALLOCATE(lat(nproma,nblks))

    IF (.NOT. ASSOCIATED(philat_2d)) &
      CALL finish(TRIM(routine), 'Latitudes of decomposition not computed yet')
    lat(:,:) = philat_2d(:,:)

  END FUNCTION get_lat

  FUNCTION get_area(decomp) RESULT(area)

#ifdef STANDALONE
    USE mo_transpose,          ONLY: reorder
    USE mo_decomposition,      ONLY: ldc => local_decomposition
    USE mo_gaussgrid,          ONLY: gridarea
#else
    USE mo_geoloc, ONLY: gboxarea_2d
#endif

    TYPE(t_patch), INTENT(in)  :: decomp
    REAL(wp),      POINTER     :: area(:,:)

    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_area'

    nproma = get_nproma(decomp)
    nblks  = get_nblks (decomp)
    ALLOCATE(area(nproma,nblks))

#ifdef STANDALONE
    ! in jsb standalone init_geoloc is never called, therefore gboxarea_2d has not been calculated.
    CALL reorder (area ,SPREAD(gridarea(ldc%glat(:)), 1, ldc%nglon))
#else
    IF (.NOT. ASSOCIATED(gboxarea_2d)) &
      CALL finish(TRIM(routine), 'grid box surface area of decomposition not computed yet')
    area(:,:) = gboxarea_2d(:,:)
#endif

  END FUNCTION get_area

  !>
  !! @brief Read number of lons and lats from a given netcdf file
  !!
  !! Currently only working for ECHAM
  !!
  !! @par History
  !! JN: Initial parts copied from jsbach3 mo_jsbach_grid.f90:init_grid
  !!
  SUBROUTINE get_nlon_nlat_from_file(grid_filename, return_nlon, return_nlat) !filename

    USE mo_jsb_parallel_iface,  ONLY: my_process_is_stdio, p_io, p_bcast, my_process_is_mpi_parallel, mpi_comm
    USE mo_jsb_io_netcdf_iface, ONLY: t_input_file, netcdf_open_input
    USE mo_netCDF,              ONLY: io_inq_dimid, io_inq_dimlen
    USE mo_control,             ONLY: nlon, ngl, nhgl, nvclev, vct

    CHARACTER(len=*), INTENT(in) :: grid_filename
    INTEGER, INTENT(out) :: return_nlon, return_nlat

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_nlon_nlat_from_file'

    TYPE(t_input_file) :: input_file
    INTEGER :: IO_dim_id

    INTEGER :: read_nlon, read_nlat

    IF (my_process_is_stdio()) THEN
      input_file =  netcdf_open_input(TRIM(grid_filename))
      CALL io_inq_dimid  (input_file%file_id, 'lon', IO_dim_id)
      CALL io_inq_dimlen (input_file%file_id, IO_dim_id, read_nlon)
      CALL io_inq_dimid  (input_file%file_id, 'lat', IO_dim_id)
      CALL io_inq_dimlen (input_file%file_id, IO_dim_id, read_nlat)
      CALL input_file%Close()
    END IF

    IF (my_process_is_mpi_parallel()) THEN
      CALL p_bcast (read_nlon, p_io, comm=mpi_comm)
      CALL p_bcast (read_nlat, p_io, comm=mpi_comm)
    ENDIF

    !Assertion: number of lats needs to be even
    IF (.NOT. MOD(read_nlat,2) .EQ. 0) CALL finish(TRIM(routine), 'number of lats needs to be even.')


    return_nlon = read_nlon
    return_nlat = read_nlat

    !Required for standalone model for mo_gaussgrid::inigau()
    nlon = read_nlon
    ngl = read_nlat
    nhgl = ngl / 2
    nvclev = 0

  END SUBROUTINE get_nlon_nlat_from_file

  !>
  !! @brief Read actual lons and lats from the given netcdf file
  !!
  !! Currently only working for ECHAM
  !!    -> read lons and lats and write to philon_2d and philat_2d
  !!       (in the standalone case init_geoloc is not called, but
  !!        philon_2d and philat_2d are required for get_lon and get_lat
  !!        and to get the output right in the site-level/part grid case)
  !!
  !! @par History
  !! First version: JN (2016-02-15)
  !!
  SUBROUTINE get_lons_and_lats_from_file(grid_filename, nproma, nblks) !filename

    USE mo_jsb_parallel_iface,  ONLY: my_process_is_stdio !, p_io, my_process_is_mpi_parallel
    USE mo_jsb_io_netcdf_iface, ONLY: t_input_file, netcdf_open_input
    USE mo_netCDF,              ONLY: nf_get_vara_double, nf_inq_varid, nf_check
    USE mo_tr_scatter,          ONLY: scatter_field
    USE mo_control,             ONLY: nlon, ngl

    USE mo_geoloc,              ONLY: philon_2d, philat_2d

    CHARACTER(len=*), INTENT(IN) :: grid_filename
    INTEGER,          INTENT(IN) :: nproma, nblks
    CHARACTER(len=*), PARAMETER :: routine = modname//':get_lons_and_lats_from_file'

    REAL(wp), ALLOCATABLE :: read_lons(:), read_lats(:)
    REAL(wp), POINTER :: lons_2D(:,:), lats_2D(:,:)
    TYPE(t_input_file) :: input_file
    INTEGER :: IO_var_id, jlat, jlon

    ALLOCATE(read_lons(nlon))
    ALLOCATE(lons_2D(nlon,ngl))
    ALLOCATE(philon_2d(nproma, nblks))

    ALLOCATE(read_lats(ngl))
    ALLOCATE(lats_2D(nlon,ngl))
    ALLOCATE(philat_2d(nproma, nblks))

    IF (my_process_is_stdio()) THEN
      input_file = netcdf_open_input(TRIM(grid_filename))
      CALL nf_check(nf_inq_varid(input_file%file_id, 'lon', IO_var_id))
      CALL nf_check(nf_get_vara_double(input_file%file_id, IO_var_id, (/ 1 /), &
        & (/ nlon /), read_lons))
      CALL nf_check(nf_inq_varid(input_file%file_id, 'lat', IO_var_id))
      CALL nf_check(nf_get_vara_double(input_file%file_id, IO_var_id, (/ 1 /), &
        & (/ ngl /), read_lats))
      CALL input_file%Close()

      lons_2D = SPREAD(read_lons, DIM=2, NCOPIES=ngl)
      lats_2D = SPREAD(read_lats, DIM=1, NCOPIES=nlon)
    ENDIF

    CALL scatter_field(lons_2D, philon_2d)
    CALL scatter_field(lats_2D, philat_2d)
    DEALLOCATE(read_lons, read_lats, lons_2D, lats_2D)

  END SUBROUTINE get_lons_and_lats_from_file

END MODULE mo_jsb_grid_iface

!! ==============================================================================================================================
!>
!! @brief Contains adapter to different utilities
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2016-02-29)
!!
MODULE mo_jsb_utils_iface

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: assign_if_present, assign_if_present_allocatable

  INTERFACE assign_if_present
    MODULE PROCEDURE assign_if_present_character
    MODULE PROCEDURE assign_if_present_logical
    MODULE PROCEDURE assign_if_present_integer
    MODULE PROCEDURE assign_if_present_integers
    MODULE PROCEDURE assign_if_present_real
  END INTERFACE

  INTERFACE assign_if_present_allocatable
    MODULE PROCEDURE assign_if_present_logical_allocatable_1d
    MODULE PROCEDURE assign_if_present_integer_allocatable
    MODULE PROCEDURE assign_if_present_integer_allocatable_1d
    MODULE PROCEDURE assign_if_present_real_allocatable
    MODULE PROCEDURE assign_if_present_real_allocatable_1d
    MODULE PROCEDURE assign_if_present_character_allocatable
  END INTERFACE assign_if_present_allocatable

  INTEGER, PARAMETER :: SUCCESS = 0

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_grid_iface'

CONTAINS

  SUBROUTINE assign_if_present_character (y,x)
    CHARACTER(len=*) ,INTENT(inout)        :: y
    CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT. PRESENT(x)) RETURN
    IF ( x == ' ' )       RETURN
    y = x
  END SUBROUTINE assign_if_present_character
!------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_logical (y,x)
    LOGICAL ,INTENT(inout)        :: y
    LOGICAL ,INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_logical
!------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_integer (y,x)
    INTEGER ,INTENT(inout)        :: y
    INTEGER ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT. PRESENT(x)) RETURN
    IF ( x == -HUGE(x)  ) RETURN
    y = x
  END SUBROUTINE assign_if_present_integer
!------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_integers (y,x)
    INTEGER ,INTENT(inout)        :: y (:)
    INTEGER ,INTENT(in) ,OPTIONAL :: x (:)
    INTEGER :: n
    IF (PRESENT(x)) THEN
      n = MIN(SIZE(x), SIZE(y))
      y(1:n) = x(1:n)
    ENDIF
  END SUBROUTINE assign_if_present_integers
!------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_real (y,x)
    REAL(wp), INTENT(inout)        :: y
    REAL(wp) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    IF ( x == -HUGE(x) ) RETURN
    y = x
  END SUBROUTINE assign_if_present_real

  SUBROUTINE assign_if_present_logical_allocatable_1d(y, x)
    LOGICAL, ALLOCATABLE, INTENT(INOUT) :: y(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: x(:)

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_logical_allocatable_1d"

    IF(.NOT.PRESENT(x)) RETURN
    IF(ALLOCATED(y)) THEN
        IF(SIZE(y) /= SIZE(x)) DEALLOCATE(y)
    END IF
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y(SIZE(x)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y(:) = x(:)
  END SUBROUTINE assign_if_present_logical_allocatable_1d

  SUBROUTINE assign_if_present_integer_allocatable(y, x)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: y
    INTEGER, OPTIONAL, INTENT(IN) :: x

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_integer_allocatable"

    IF(.NOT.PRESENT(x)) RETURN
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y = x
  END SUBROUTINE assign_if_present_integer_allocatable

  SUBROUTINE assign_if_present_integer_allocatable_1d(y, x)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: y(:)
    INTEGER, OPTIONAL, INTENT(IN) :: x(:)

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_integer_allocatable_1d"

    IF(.NOT.PRESENT(x)) RETURN
    IF(ALLOCATED(y)) THEN
        IF(SIZE(y) /= SIZE(x)) DEALLOCATE(y)
    END IF
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y(SIZE(x)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y(:) = x(:)
  END SUBROUTINE assign_if_present_integer_allocatable_1d

  SUBROUTINE assign_if_present_real_allocatable(y, x)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: y
    REAL(wp), OPTIONAL, INTENT(IN) :: x

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_real_allocatable"

    IF(.NOT.PRESENT(x)) RETURN
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y = x
  END SUBROUTINE assign_if_present_real_allocatable

  SUBROUTINE assign_if_present_real_allocatable_1d(y, x)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: y(:)
    REAL(wp), OPTIONAL, INTENT(IN) :: x(:)

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_real_allocatable_1d"

    IF(.NOT.PRESENT(x)) RETURN
    IF(ALLOCATED(y)) THEN
        IF(SIZE(y) /= SIZE(x)) DEALLOCATE(y)
    END IF
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y(SIZE(x)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y(:) = x(:)
  END SUBROUTINE assign_if_present_real_allocatable_1d

  SUBROUTINE assign_if_present_character_allocatable(y,x)
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(INOUT) :: y
    CHARACTER(LEN=*), OPTIONAL,    INTENT(IN)    :: x

    IF (.NOT. PRESENT(x)) RETURN
    IF ( TRIM(x) == '' )  RETURN
    y = x
  END SUBROUTINE assign_if_present_character_allocatable

END MODULE mo_jsb_utils_iface

!! ==============================================================================================================================
!>
!! @brief Contains the interface to varlists/streams
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_varlist_iface

  USE mo_kind,               ONLY: wp, dp
  USE mo_exception,          ONLY: finish
  USE mo_memory_base,        ONLY: new_stream,                       &
                                   get_var_list => get_stream,       &
                                   add_stream_element,               &
                                   default_stream_setting
  USE mo_linked_list,        ONLY: t_var_list => t_stream,           &
                                   t_var_metadata => memory_info,    &
                                   t_list_element => list_element,   &
                                   find_list_element

  USE mo_jsb_io_iface, ONLY: t_cf_var, t_grib2_var

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: VARNAME_LEN
  PUBLIC :: t_var_list, t_var_metadata, t_list_element, get_var_list
  PUBLIC :: new_var_list
  PUBLIC :: add_var_list_element_r2d, add_var_list_element_r3d
  PUBLIC :: is_variable_in_output   ! a helper function that returns always TRUE [originally implemented in ICON with more code]

  INTEGER, PARAMETER :: VARNAME_LEN = 64

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_varlist_iface'

CONTAINS

  SUBROUTINE new_var_list (this_list, name, patch_id, output_type, restart_type, &
       &                       post_suf, rest_suf, init_suf, loutput, lrestart,  &
       &                       linitial, table)
    !
    TYPE(t_var_list), POINTER, INTENT(inout) :: this_list    ! anchor
    CHARACTER(len=*), INTENT(in)             :: name         ! name of output var_list
    INTEGER,          INTENT(in)             :: patch_id     ! patch ID
    INTEGER,          INTENT(in), OPTIONAL   :: output_type  ! 'GRIB1', 'GRIB2' or 'NetCDF[12]'
    INTEGER,          INTENT(in), OPTIONAL   :: restart_type ! 'NetCDF[12]'
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: post_suf     ! suffix of output file
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: rest_suf     ! suffix of restart file
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: init_suf     ! suffix of initial file
    LOGICAL,          INTENT(in), OPTIONAL   :: loutput      ! write to  output file
    LOGICAL,          INTENT(in), OPTIONAL   :: lrestart     ! write to restart file
    LOGICAL,          INTENT(in), OPTIONAL   :: linitial     ! read from initial file
    INTEGER,          INTENT(in), OPTIONAL   :: table        ! used only for ECHAM

    INTEGER :: table_loc

    IF(PRESENT(table)) THEN
      table_loc = table
    ELSE
      table_loc = 180
    END IF
    CALL new_stream(this_list, name,                                         &
                    filetype=output_type,                                    &
                    post_suf=post_suf, rest_suf=rest_suf, init_suf=init_suf, &
                    lpost=loutput, lpout=loutput, lrerun=lrestart            &
                   )

    !@todo: hardcoded FALSE for now to reduce output size when running with ECHAM
    !       until the "in_groups" behavior is also implemented for ECHAM!
    ! - now: for required vars loutput needs to be set to true! (+ see namelist selector in jsb3 sitelevel)
    !@todo: implement the "in_groups" behavior for ECHAM!
    CALL default_stream_setting(this_list, table=table_loc, lpost=.FALSE., lrerun=lrestart)

  END SUBROUTINE new_var_list

  SUBROUTINE add_var_list_element_r2d(this_list, name, ptr,                             &
    hgrid, vgrid, cf, grib2, code, table, ldims, gdims, levelindx, loutput, lcontainer, &
    lrestart, lrestart_cont, initval_r, isteptype,                                      &
    resetval_r, lmiss, missval_r, tlev_source, info, p5,                                &
    in_groups, verbose, new_element                                                     &
    )

    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:,:)            ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type
                                                                      ! ICON: zaxis type (ZA_*)
                                                                      ! ECHAM: cdi zaxis ID
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in)           :: code                ! GRIB1 code number
    INTEGER,              INTENT(in)           :: table               ! GRIB1 table number
    INTEGER,              INTENT(in), OPTIONAL :: ldims(2)            ! local dimensions
    INTEGER,              INTENT(in), OPTIONAL :: gdims(2)            ! global dimensions
    INTEGER,              INTENT(in), OPTIONAL :: levelindx           ! info%levelindx for ECHAM
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(wp),             TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    CHARACTER(len=VARNAME_LEN), INTENT(in), OPTIONAL :: in_groups(:)  ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element           ! pointer to new var list element

    REAL(wp), POINTER :: p4(:,:,:,:)

    TYPE (t_list_element), POINTER :: element

    CHARACTER(len=*), PARAMETER :: routine = modname//':add_var_list_element_r2d'

    IF (PRESENT(p5)) THEN
      IF (SIZE(p5, DIM=5) > 1) &
        CALL finish(TRIM(routine), 'p5: only four dimensions allowed currently because of ECHAM compatibility')
      p4 => p5(:,:,:,:,1)
    END IF

    IF (PRESENT(p5)) THEN
      CALL add_stream_element(this_list, TRIM(name), ptr, &
        ldims=ldims, gdims=gdims, leveltype=levelindx, &
        longname=cf%long_name, units=cf%units, &
        code=code, table=table, &
        lpost=loutput, lrerun=lrestart, contnorest=lrestart_cont, &
        lmiss=lmiss, missval=missval_r, p4=p4, verbose=verbose &
        )
    ELSE
      CALL add_stream_element(this_list, TRIM(name), ptr, &
        ldims=ldims, gdims=gdims, leveltype=levelindx, &
        longname=cf%long_name, units=cf%units, &
        code=code, table=table, &
        lpost=loutput, lrerun=lrestart, contnorest=lrestart_cont, &
        lmiss=lmiss, missval=missval_r, verbose=verbose &
        )
    END IF
    element => find_list_element(this_list, TRIM(name))
    element%field%info%zaxisID = vgrid

    IF (PRESENT(initval_r)) THEN
      ptr(:,:) = initval_r
    END IF

  END SUBROUTINE add_var_list_element_r2d

  SUBROUTINE add_var_list_element_r3d(this_list, name, ptr,                             &
    hgrid, vgrid, cf, grib2, code, table, ldims, gdims, levelindx, loutput, lcontainer, &
    lrestart, lrestart_cont, initval_r, isteptype,                                      &
    resetval_r, lmiss, missval_r, tlev_source, info, p5,                                &
    in_groups, verbose, new_element                                                     &
    )

    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:,:,:)          ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type
                                                                      ! ICON: zaxis type (ZA_*)
                                                                      ! ECHAM: cdi zaxis ID
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in)           :: code                ! GRIB1 code number
    INTEGER,              INTENT(in)           :: table               ! GRIB1 table number
    INTEGER,              INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions
    INTEGER,              INTENT(in), OPTIONAL :: gdims(3)            ! global dimensions
    INTEGER,              INTENT(in), OPTIONAL :: levelindx           ! info%levelindx for ECHAM
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(wp),             TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     ,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    CHARACTER(len=VARNAME_LEN), INTENT(in), OPTIONAL :: in_groups(:)  ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element           ! pointer to new var list element

    REAL(wp), POINTER :: p4(:,:,:,:)

    TYPE (t_list_element), POINTER :: element

    CHARACTER(len=*), PARAMETER :: routine = modname//':add_var_list_element_r3d'

    IF (PRESENT(p5)) THEN
      IF (SIZE(p5, DIM=5) > 1) &
        CALL finish(TRIM(routine), 'p5: only four dimensions allowed currently because of ECHAM compatibility')
      p4 => p5(:,:,:,:,1)
    END IF

    IF (PRESENT(p5)) THEN
      CALL add_stream_element(this_list, TRIM(name), ptr, &
        ldims=ldims, gdims=gdims, leveltype=levelindx, &
        longname=cf%long_name, units=cf%units, &
        code=code, table=table, &
        lpost=loutput, lmiss=lmiss, missval=missval_r, lrerun=lrestart, contnorest=lrestart_cont, &
        p4=p4 &
      !mem_info=info &
        )
    ELSE
      CALL add_stream_element(this_list, TRIM(name), ptr, &
        ldims=ldims, gdims=gdims, leveltype=levelindx, &
        longname=cf%long_name, units=cf%units, &
        code=code, table=table, &
        lpost=loutput, lmiss=lmiss, missval=missval_r, lrerun=lrestart, contnorest=lrestart_cont &
        )
    END IF

    element => find_list_element(this_list, TRIM(name))
    element%field%info%zaxisID = vgrid

    IF (PRESENT(initval_r)) THEN
      ptr(:,:,:) = initval_r
    END IF

  END SUBROUTINE add_var_list_element_r3d


  !! helper function that always returns TRUE
  !! needed in jsbach4:mo_jsb_varlist after implementing GPU code
  !! in ICON this function checks whether a variable creates output
  !! in ECHAM this functionality is not implemented, and, thus, this function returns always TRUE, stating the var creates output
  FUNCTION is_variable_in_output(full_name, in_groups) RESULT(l_var_output)
    
    IMPLICIT NONE

    ! ---------------------------
    ! 0.1 InOut
    CHARACTER(len=*),                     INTENT(in) :: full_name
    CHARACTER(len=VARNAME_LEN), OPTIONAL, INTENT(in) :: in_groups(:) ! groups to which this variable belongs to
    LOGICAL                                          :: l_var_output
    
    
    l_var_output = .TRUE.
  
  END FUNCTION is_variable_in_output

END MODULE mo_jsb_varlist_iface

!! ==============================================================================================================================
!>
!! @brief Contains the adapter for mo_jsb_vertical_axes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2017-03-14)
!!
MODULE mo_jsb_vertical_axes_iface

  USE mo_kind, ONLY: wp !, dp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_verticalAxisList, Setup_jsb_vertical_axis, t_Vgrid, Set_jsb_restart_vgrid

  ! Dummy
  TYPE t_verticalAxis
  END TYPE t_verticalAxis
  TYPE t_verticalAxisList
  END TYPE t_verticalAxisList
  TYPE t_Vgrid
  END TYPE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_vertical_axes_iface'

CONTAINS

  SUBROUTINE Setup_jsb_vertical_axis(verticalAxisList, zaxis_id, name, length, longname, units, levels, lbounds, ubounds)
    ! Dummy
    TYPE(t_verticalAxisList), INTENT(inout) :: verticalAxisList
    INTEGER,                  INTENT(in)    :: zaxis_id
    CHARACTER(len=*),         INTENT(in)    :: name
    INTEGER,                  INTENT(in)    :: length
    CHARACTER(len=*),         INTENT(in)    :: longname
    CHARACTER(len=*),         INTENT(in)    :: units
    REAL(wp),                 INTENT(in)    :: levels(:), lbounds(:), ubounds(:)

  END SUBROUTINE Setup_jsb_vertical_axis

  SUBROUTINE Set_jsb_restart_vgrid(vgrid_defs, count, zaxis_id, levels)
    ! Dummy
    TYPE(t_Vgrid), INTENT(inout) :: vgrid_defs(:)
    INTEGER,       INTENT(inout) :: count
    INTEGER, VALUE               :: zaxis_id
    REAL(wp),      INTENT(in)    :: levels(:)
  END SUBROUTINE Set_jsb_restart_vgrid

END MODULE mo_jsb_vertical_axes_iface

!! ==============================================================================================================================
!>
!! @brief Contains adapters to physcial constants
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2016-03-18)
!!
MODULE mo_physical_constants_iface

  USE mo_kind, ONLY: wp

  USE mo_physical_constants, ONLY: &
    grav,            & !< Avg. gravitational acceleration        [m/s2]
    stbo,            & !< Stephan-Boltzmann constant             [W/(m2 K4)]
    argas,           & !< Molar/universal/ideal gas constant     [J/(K mol)]
    zemiss_def,      & !< Long-wave surface emissivity factor    []
    rd,              & !< Gas constant for dry air               [J/(K kg)]
    rv,              & !< Gas constant for water vapor           [J/(K*kg)]
    cpd,             & !< Specific heat at constant pressure     [J/(K kg)]
    cvd,             & !< Specific heat at constant volume       [J/(K kg)]
    rd_o_cpd,        & !< rd/cpd
    rhoh2o,          & !< Density of liquid water                [kg/m3]
    alv,             & !< Latent heat for vaporisation           [J/kg]
    als,             & !< Latent heat for sublimation            [J/kg]
    alf,             & !< Latent heat for fusion                 [J/kg]
    tmelt,           & !< Melting temperature of ice/snow        [K]
    amd,             & !< Molar weight of dry air                [g/mol]
    amco2,           & !< Molar weight of CO2                    [g/mol]
    clw => cpliq,    & !< Specific heat capacity of liquid water [J/(K kg)]
    ci => cpice,     & !< Specific heat capacity of ice          [J/(K kg)]
    cs => cpsno,     & !< Specific heat capacity of snow         [J/(K kg)]
    rhoi => rhoice,  & !< Density of (sea) ice                   [kg/m3]
    rhos => rhosno,  & !< Density of snow                        [kg/m3]
    ki => alice,     & !< Heat conductivity of ice               [J/(m s K)]
    ks => alsno,     & !< Heat conductivity of snow              [J/(m s K)]
    ! Auxiliary constants
    rvd1  => vtmpc1, & !< = rv/rd-1                              []
    cpvd1 => vtmpc2    !< = cpv/cpd-1                            []

! cvdifts in echam is not a PARAMETER, therefore it can't be USEd from modules here
! TBD: Use cvdifts=1 for standalone model?
!#ifndef STANDALONE
!#ifdef __ICON__
!  USE mo_echam_vdiff_params, ONLY: cvdifts
!#else
!  USE mo_physc2, ONLY: cvdifts
!#endif
!#endif

  IMPLICIT NONE
  PUBLIC

!#ifdef STANDALONE
!  REAL(wp), PARAMETER :: cvdifts = 1._wp
!#else
  REAL(wp), PARAMETER :: cvdifts = 1.5_wp
!#endif
  REAL(wp), PARAMETER :: tpfac1 = cvdifts
  REAL(wp), PARAMETER :: tpfac2 = 1._wp / tpfac1
  REAL(wp), PARAMETER :: tpfac3 = 1._wp - tpfac2

END MODULE mo_physical_constants_iface

!! ==============================================================================================================================
!>
!! @brief Contains interfaces to ICON convect_tables
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                 by Reiner Schnur (2022-03-07)
!!
MODULE mo_jsb_convect_tables_iface

  USE mo_kind,                 ONLY: wp
  USE mo_echam_convect_tables, ONLY: init_convect_tables, tlucua, jptlucu1, jptlucu2

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_convect_tables, tlucua, jptlucu1, jptlucu2

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_convect_tables_iface'

! CONTAINS

END MODULE mo_jsb_convect_tables_iface

!! ==============================================================================================================================
!>
!! @brief Contains interfaces to ICON orbital and solar functions
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                 by Reiner Schnur (2022-08-26)
!!
MODULE mo_jsb_orbit_solar_iface

  USE mo_orbit,                ONLY: inquire_declination, orbit_kepler, orbit_vsop87

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: inquire_declination, compute_cos_zenith_angle

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_orbit_solar_iface'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Compute cosine of zenith angle (only used for standalone land model)
  !
  SUBROUTINE compute_cos_zenith_angle(sinlon, sinlat, coslon, coslat, cos_zenith_angle)

    USE mo_kind,                 ONLY: wp
    USE mo_radiation_parameters, ONLY: lyr_perp, nmonth, yr_perp, solar_parameters
    USE mo_time_control,         ONLY: l_orbvsop87, get_orbit_times

    REAL(wp), INTENT(in)  :: sinlon(:,:), sinlat(:,:), coslon(:,:), coslat(:,:)
    REAL(wp), INTENT(out) :: cos_zenith_angle(:,:)

    REAL(wp) :: daylight_frc(SIZE(cos_zenith_angle,1),SIZE(cos_zenith_angle,2))
    LOGICAL  :: l_rad_call
    REAL(wp) :: time_of_day, orbit_date, rasc_sun, decl_sun, dist_sun, flx_ratio

    l_rad_call = .FALSE.
    CALL get_orbit_times(l_rad_call, lyr_perp, nmonth, yr_perp, time_of_day, orbit_date)

    IF (l_orbvsop87) THEN
      CALL orbit_vsop87(orbit_date, rasc_sun, decl_sun, dist_sun) ! in, out, out, out
    ELSE
      CALL orbit_kepler(orbit_date, rasc_sun, decl_sun, dist_sun) ! in, out, out, out
    END IF

    CALL solar_parameters(decl_sun, dist_sun, time_of_day,   &
      & sinlon(:,:), sinlat(:,:), coslon(:,:), coslat(:,:),  &
      & flx_ratio, cos_zenith_angle(:,:), daylight_frc(:,:))

  END SUBROUTINE compute_cos_zenith_angle

END MODULE mo_jsb_orbit_solar_iface

!! ==============================================================================================================================

#else

MODULE mo_dummy_module

END MODULE mo_dummy_module

#endif
