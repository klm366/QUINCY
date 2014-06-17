!>
!! @brief Contains JSBACH grid structure and methods
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_grid_class
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max
  USE mo_exception,          ONLY: message_text, message, print_value, finish
  USE mo_jsb_domain_iface,   ONLY: t_patch, p_patch

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_grid
  PUBLIC :: new_grid, delete_grid, print_grid
  PUBLIC :: get_nblks, get_nproma

  TYPE t_jsb_grid
    INTEGER                     :: id = 0
    CHARACTER(LEN=10)           :: name = ''
    CHARACTER(LEN=filename_max) :: filename = ''
    INTEGER                     :: ntotal = 0          !< Number of cells on this PE
    INTEGER                     :: nland = 0
    REAL(wp), POINTER           :: lon(:,:) => NULL()
    REAL(wp), POINTER           :: lat(:,:) => NULL()
    LOGICAL,  POINTER           :: lsm(:,:) => NULL()
    REAL(wp), POINTER           :: lsf(:,:) => NULL()
    INTEGER                     :: nblks = 1
    INTEGER                     :: nproma = 0
    INTEGER                     :: npromz = 0
    INTEGER                     :: ntotal_g            !< Total number of cells in grid
    LOGICAL                     :: configured = .FALSE.
    TYPE(t_patch), POINTER      :: patch
  END TYPE t_jsb_grid

  INTERFACE new_grid
    MODULE PROCEDURE new_grid_from_file
    MODULE PROCEDURE new_grid_from_patch
    MODULE PROCEDURE new_grid_from_lonlat
  END INTERFACE new_grid

  INTEGER, SAVE :: ngrid_count = 0

CONTAINS

  !>
  !! Grid constructor

  FUNCTION new_grid_from_file(filename, name) RESULT(grid)

    USE mo_jsb_io_netcdf,      ONLY: ncfile_has_dim
    USE mo_jsb_parallel_iface, ONLY: my_process_is_stdio, p_io, p_bcast, mpi_comm

    CHARACTER(LEN=*), INTENT(in)    :: filename
    CHARACTER(LEN=*), INTENT(in)    :: name
    TYPE(t_jsb_grid), POINTER       :: grid

    TYPE(t_jsb_grid), POINTER       :: new

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid_class:new_grid_from_file'

    CALL message(TRIM(routine), 'starting construction of new JSBACH grid instance from file:')
    CALL message('     ', TRIM(filename))

    ALLOCATE(new)

    ngrid_count = ngrid_count + 1
    new%id = ngrid_count
    new%name = name
    new%filename = filename

    IF (my_process_is_stdio()) THEN

      IF (ncfile_has_dim(TRIM(filename), 'lon_cell_centre')) THEN
        ! ICON grid
      ELSE
        ! ECHAM grid
      END IF

    END IF

    !CALL p_bcast(, p_io, mpi_comm)

    new%configured = .TRUE.

    grid => new

    CALL message(TRIM(routine), 'construction of JSBACH grid completed.')

  END FUNCTION new_grid_from_file

  FUNCTION new_grid_from_patch(patch) RESULT(grid)

!!$    TYPE(t_patch),    TARGET, INTENT(in) :: patch
    TYPE(t_patch),             INTENT(in) :: patch
    TYPE(t_jsb_grid), POINTER             :: grid

    TYPE(t_jsb_grid), POINTER       :: new
    INTEGER :: ncells, nblks, nproma, npromz

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid_class:new_grid_from_patch'

    CALL message(TRIM(routine), 'starting construction of new JSBACH grid instance from patch')

    ALLOCATE(new)

    ngrid_count = ngrid_count + 1
    new%id = ngrid_count
!!$    new%name = name
    new%filename = patch%grid_filename

    ncells = patch%n_patch_cells
    nblks = patch%nblks_c
    npromz = patch%npromz_c
    nproma = INT( (ncells-npromz)/(nblks-1) )
    new%ntotal = ncells
    new%nblks = nblks
    new%nproma = nproma
    new%npromz = npromz
    ALLOCATE(new%lon(nproma,nblks), new%lat(nproma,nblks), &
             new%lsm(nproma,nblks), new%lsf(nproma,nblks))
    new%lon(:,:) = patch%cells%center(:,:)%lon
    new%lat(:,:) = patch%cells%center(:,:)%lat
    new%lsm(:,:) = .TRUE.
    new%lsf(:,:) = 1.0_wp
    new%ntotal_g = patch%n_patch_cells_g

!!$    ALLOCATE(new%patch)
!!$    new%patch%id = patch%id
!!$    new%patch%n_patch_cells = patch%n_patch_cells
!!$    new%patch%n_patch_cells_g = patch%n_patch_cells_g
!!$    new%patch%nblks_c = patch%nblks_c
!!$    new%patch%npromz_c = patch%npromz_c
!!$    new%patch%grid_filename = patch%grid_filename
!!$    new%patch%comm_pat_gather_c = patch%comm_pat_gather_c
!!$    new%patch%cells%glb_index = patch%cells%glb_index
    new%patch => p_patch(patch%id)

    new%configured = .TRUE.

    grid => new

    CALL message(TRIM(routine), 'construction of JSBACH grid completed.')

  END FUNCTION new_grid_from_patch

  FUNCTION new_grid_from_lonlat() RESULT(grid)

    TYPE(t_jsb_grid), POINTER       :: grid

    TYPE(t_jsb_grid), POINTER       :: new

    ALLOCATE(new)

    grid => new

  END FUNCTION new_grid_from_lonlat

  SUBROUTINE delete_grid(self)

#ifdef HAVE_F2003
    TYPE(t_jsb_grid), POINTER, INTENT(inout) :: self
#else
    TYPE(t_jsb_grid), POINTER :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid_class:delete_grid'

    CALL message(TRIM(routine), 'starting destruction of JSBACH grid instance.')

    IF (ASSOCIATED(self%lsm)) DEALLOCATE(self%lsm)
    IF (ASSOCIATED(self%lsf)) DEALLOCATE(self%lsf)
    IF (ASSOCIATED(self%lon)) DEALLOCATE(self%lon)
    IF (ASSOCIATED(self%lat)) DEALLOCATE(self%lat)

    DEALLOCATE(self)

    CALL message(TRIM(routine), 'destruction of JSBACH grid completed.')

  END SUBROUTINE delete_grid

  !>
  !! Print method for grid instance.
  !!
  !! Subroutine to print description of grid instance.
  !!
  SUBROUTINE print_grid(self)

#ifdef HAVE_F2003
    TYPE(t_jsb_grid), POINTER, INTENT(in) :: self
#else
    TYPE(t_jsb_grid), POINTER             :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid_class:print_grid'

    CALL message(TRIM(routine), 'print JSBACH grid description')

    WRITE(message_text,'(A,I2,A,A)') 'Grid id: ', self%id, ', grid name: ', TRIM(self%name)
    CALL message('  ', message_text)

    IF (.NOT. self%configured) THEN
      CALL message(TRIM(routine), 'Grid not configured yet!')
      RETURN
    END IF

    IF (self%filename /= '') THEN
      WRITE(message_text,'(A)') 'Grid filename: '//TRIM(self%filename)
      CALL message('  ', message_text)
    END IF

    IF (self%patch%id >= 0) THEN
      WRITE(message_text,'(A,I2)') 'Patch ID: ', self%patch%id
      CALL message('  ', message_text)
    END IF

  END SUBROUTINE print_grid

  INTEGER FUNCTION get_nblks(self)

#ifdef HAVE_F2003
    TYPE(t_jsb_grid), POINTER, INTENT(in) :: self
#else
    TYPE(t_jsb_grid), POINTER             :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid_class:get_nblks'

    get_nblks = self%nblks

  END FUNCTION get_nblks

  INTEGER FUNCTION get_nproma(self)

#ifdef HAVE_F2003
    TYPE(t_jsb_grid), POINTER, INTENT(in) :: self
#else
    TYPE(t_jsb_grid), POINTER             :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid_class:get_nproma'

    get_nproma = self%nproma

  END FUNCTION get_nproma

  INTEGER FUNCTION get_blk_start(grid)

    USE mo_jsb_domain_iface, ONLY: get_blk_start_domain => get_blk_start

    TYPE(t_jsb_grid), INTENT(in) :: grid

    get_blk_start = get_blk_start_domain(grid%patch)

  END FUNCTION get_blk_start

  INTEGER FUNCTION get_blk_end(grid)

    USE mo_jsb_domain_iface, ONLY: get_blk_end_domain => get_blk_end

    TYPE(t_jsb_grid), INTENT(in) :: grid

    get_blk_end = get_blk_end_domain(grid%patch)

  END FUNCTION get_blk_end

  INTEGER FUNCTION get_col_start(iblk, grid)

    USE mo_jsb_domain_iface, ONLY: get_col_start_domain => get_col_start

    INTEGER,          INTENT(in) :: iblk
    TYPE(t_jsb_grid), INTENT(in) :: grid

    get_col_start = get_col_start_domain(iblk, grid%patch)

  END FUNCTION get_col_start

  INTEGER FUNCTION get_col_end(iblk, grid)

    USE mo_jsb_domain_iface, ONLY: get_col_end_domain => get_col_end

    INTEGER,          INTENT(in) :: iblk
    TYPE(t_jsb_grid), INTENT(in) :: grid

    get_col_end = get_col_end_domain(iblk, grid%patch)

  END FUNCTION get_col_end

#endif
END MODULE mo_jsb_grid_class
