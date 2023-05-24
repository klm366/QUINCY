!>
!! @brief Contains JSBACH grid structure and methods
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_grid_class
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_jsb_domain_iface,   ONLY: t_patch, get_host_patch_id,                    &
                                   get_grid_filename_domain => get_grid_filename, &
                                   get_ntotal_domain        => get_ntotal,        &
                                   get_ntotal_g_domain      => get_ntotal_g,      &
                                   get_dims_g_domain        => get_dims_g,        &
                                   get_nlat_g_domain        => get_nlat_g,        &
                                   get_nproma_domain        => get_nproma,        &
                                   get_nblks_domain         => get_nblks
  USE mo_jsb_grid_iface,     ONLY: get_lon_domain           => get_lon,           &
                                   get_lat_domain           => get_lat,           &
                                   get_area_domain          => get_area
  USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_grid, t_jsb_grid_ptr
  PUBLIC :: new_grid, delete_grid

  TYPE t_jsb_grid
    INTEGER                       :: id = 0
    CHARACTER(LEN=SHORT_NAME_LEN) :: name = ''
    CHARACTER(LEN=filename_max)   :: filename = ''
    INTEGER                       :: ntotal = 0          !< Number of cells on this PE
    INTEGER                       :: nland = 0
    REAL(wp), POINTER             :: lon(:,:) => NULL()  !< Longitudes of grid in degrees
    REAL(wp), POINTER             :: lat(:,:) => NULL()  !< Latitudes of grid in degrees
    REAL(wp), POINTER             :: coslat(:,:) => NULL()
    REAL(wp), POINTER             :: sinlat(:,:) => NULL()
    REAL(wp), POINTER             :: coslon(:,:) => NULL()
    REAL(wp), POINTER             :: sinlon(:,:) => NULL()
    REAL(wp), POINTER             :: area(:,:) => NULL()
    LOGICAL,  POINTER             :: lsm(:,:) => NULL()
    REAL(wp), POINTER             :: lsf(:,:) => NULL()
    INTEGER                       :: nblks = 1
    INTEGER                       :: nproma = 0
    INTEGER                       :: npromz = 0
    INTEGER                       :: ntotal_g              !< Total number of cells in grid
    INTEGER                       :: dims_g(2)             !< Dimensions of global grid
    INTEGER                       :: nlat_g = 0            !< Number of latitudes in global grid (effective for ICON)
    LOGICAL                       :: configured = .FALSE.
    CHARACTER(LEN=:), ALLOCATABLE :: type
    TYPE(t_patch), POINTER        :: patch => NULL()
    INTEGER                       :: host_patch_id
  CONTAINS
    PROCEDURE :: Get_nblks
    PROCEDURE :: Get_nproma
    PROCEDURE :: Get_npromz
    PROCEDURE :: Get_blk_start
    PROCEDURE :: Get_blk_end
    PROCEDURE :: Get_col_start
    PROCEDURE :: Get_col_end
    PROCEDURE :: Print             => Print_grid
  END TYPE t_jsb_grid

  TYPE t_jsb_grid_ptr
    TYPE(t_jsb_grid), POINTER :: p => NULL()
  END TYPE t_jsb_grid_ptr

  INTERFACE new_grid
    MODULE PROCEDURE new_grid_from_file
    MODULE PROCEDURE new_grid_from_host
    ! MODULE PROCEDURE new_grid_from_lonlat ! Not implemented, yet!
  END INTERFACE new_grid

  INTEGER, SAVE :: ngrid_count = 0

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_grid_class'

CONTAINS

  !>
  !! Grid constructor

  FUNCTION new_grid_from_host(patch, type) RESULT(grid)

    USE mo_jsb_math_constants, ONLY: pi

    TYPE(t_patch),    TARGET, INTENT(in) :: patch
    CHARACTER(LEN=*),         INTENT(in) :: type
    TYPE(t_jsb_grid),         POINTER    :: grid

    TYPE(t_jsb_grid), POINTER       :: new

    CHARACTER(len=*), PARAMETER :: routine = modname//':new_grid_from_host'

    ! CALL message(routine, 'starting construction of new JSBACH grid instance from host model')

    ALLOCATE(new)

    ngrid_count = ngrid_count + 1
    new%id = ngrid_count
!!$    new%name = name

    SELECT CASE (TRIM(type))
    CASE ('icon')
      new%patch => patch
    CASE ('echam')
      ALLOCATE(new%patch)
      new%patch = patch
    END SELECT
    new%type = TRIM(type)
    new%host_patch_id = get_host_patch_id(patch)

    new%filename = get_grid_filename_domain(patch)
    new%ntotal   = get_ntotal_domain       (patch)
    new%ntotal_g = get_ntotal_g_domain     (patch)
    new%dims_g(1:2) = get_dims_g_domain    (patch)
    new%nlat_g   = get_nlat_g_domain       (patch)
    new%nblks    = get_nblks_domain        (patch)
    new%nproma   = get_nproma_domain       (patch)
    new%npromz   = new%ntotal - new%nproma * (new%nblks - 1)
    new%lon     => get_lon_domain          (patch)
    new%lat     => get_lat_domain          (patch)
    ALLOCATE(new%coslat(new%nproma,new%nblks), new%sinlat(new%nproma,new%nblks), &
      &      new%coslon(new%nproma,new%nblks), new%sinlon(new%nproma,new%nblks))
    new%coslat(:,:) = COS(pi * new%lat(:,:) / 180._wp)
    new%sinlat(:,:) = SIN(pi * new%lat(:,:) / 180._wp)
    new%coslon(:,:) = COS(pi * new%lon(:,:) / 180._wp)
    new%sinlon(:,:) = SIN(pi * new%lon(:,:) / 180._wp)
    new%area    => get_area_domain         (patch)
    ALLOCATE(new%lsm(new%nproma,new%nblks), new%lsf(new%nproma,new%nblks))
    new%lsm(:,:) = .TRUE.
    new%lsf(:,:) = 1.0_wp

    new%configured = .TRUE.
    !$ACC ENTER DATA COPYIN(new%lat)

    grid => new

    ! CALL message(routine, 'construction of JSBACH grid completed.')

  END FUNCTION new_grid_from_host

  FUNCTION new_grid_from_file(filename, name) RESULT(grid)

    ! USE mo_jsb_io_netcdf,      ONLY: ncfile_has_dim
    ! USE mo_jsb_parallel, ONLY: my_process_is_stdio

    CHARACTER(LEN=*), INTENT(in)    :: filename
    CHARACTER(LEN=*), INTENT(in)    :: name
    TYPE(t_jsb_grid), POINTER       :: grid

    TYPE(t_jsb_grid), POINTER       :: new

    CHARACTER(len=*), PARAMETER :: routine = modname//':new_grid_from_file'

    CALL message(routine, 'starting construction of new JSBACH grid instance from file:')
    CALL message('     ', TRIM(filename))

    ALLOCATE(new)

    ngrid_count = ngrid_count + 1
    new%id = ngrid_count
    new%name = name
    new%filename = filename
    new%host_patch_id = 0

    ! IF (my_process_is_stdio()) THEN

    !   IF (ncfile_has_dim(TRIM(filename), 'lon_cell_centre')) THEN
    !     ! ICON grid
    !   ELSE
    !     ! ECHAM grid
    !   END IF

    ! END IF

    new%configured = .TRUE.

    grid => new

    CALL message(routine, 'construction of JSBACH grid completed.')

  END FUNCTION new_grid_from_file

  ! Not implemented, yet!
  ! FUNCTION new_grid_from_lonlat(lon, lat) RESULT(grid)

  !   REAL(wp)                        :: lon
  !   REAL(wp)                        :: lat
  !   TYPE(t_jsb_grid), POINTER       :: grid

  !   TYPE(t_jsb_grid), POINTER       :: new

  !   CHARACTER(len=*), PARAMETER :: routine = modname//':new_grid_from_lonlat'

  !   ALLOCATE(new)

  !   grid => new

  ! END FUNCTION new_grid_from_lonlat

  SUBROUTINE delete_grid(self)

#ifdef HAVE_F2003
    TYPE(t_jsb_grid), POINTER, INTENT(inout) :: self
#else
    TYPE(t_jsb_grid), POINTER :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':delete_grid'

    CALL message(routine, 'starting destruction of JSBACH grid instance.')

    IF (ASSOCIATED(self%lsm)) DEALLOCATE(self%lsm)
    IF (ASSOCIATED(self%lsf)) DEALLOCATE(self%lsf)
    IF (ASSOCIATED(self%lon)) DEALLOCATE(self%lon)
    IF (ASSOCIATED(self%lat)) DEALLOCATE(self%lat)
    IF (ASSOCIATED(self%coslat)) DEALLOCATE(self%coslat)
    IF (ASSOCIATED(self%sinlat)) DEALLOCATE(self%sinlat)
    IF (ASSOCIATED(self%coslon)) DEALLOCATE(self%coslon)
    IF (ASSOCIATED(self%sinlon)) DEALLOCATE(self%sinlon)

    DEALLOCATE(self)

    CALL message(routine, 'destruction of JSBACH grid completed.')

  END SUBROUTINE delete_grid

  !>
  !! Print method for grid instance.
  !!
  !! Subroutine to print description of grid instance.
  !!
  SUBROUTINE Print_grid(this)

    CLASS(t_jsb_grid), INTENT(in) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':Print_grid'

    CALL message(routine, 'print JSBACH grid description')

    WRITE(message_text,'(A,I2,A,A)') 'Grid id: ', this%id, ', grid name: ', TRIM(this%name)
    CALL message('  ', message_text)

    IF (.NOT. this%configured) THEN
      CALL message(routine, 'Grid not configured yet!')
      RETURN
    END IF

    IF (this%filename /= '') THEN
      WRITE(message_text,'(A)') 'Grid filename: '//TRIM(this%filename)
      CALL message('  ', message_text)
    END IF

    IF (this%host_patch_id > 0) THEN
      WRITE(message_text,'(A,I2)') 'Patch ID in host model: ', this%host_patch_id
      CALL message('  ', message_text)
    END IF

  END SUBROUTINE Print_grid

  INTEGER FUNCTION get_nblks(self)

#ifdef HAVE_F2003
    CLASS(t_jsb_grid), INTENT(in) :: self
#else
    CLASS(t_jsb_grid)             :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_nblks'

    get_nblks = self%nblks

  END FUNCTION get_nblks

  INTEGER FUNCTION get_nproma(self)

#ifdef HAVE_F2003
    CLASS(t_jsb_grid), INTENT(in) :: self
#else
    CLASS(t_jsb_grid)             :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_nproma'

    get_nproma = self%nproma

  END FUNCTION get_nproma

  INTEGER FUNCTION get_npromz(self)

#ifdef HAVE_F2003
    CLASS(t_jsb_grid), INTENT(in) :: self
#else
    CLASS(t_jsb_grid)             :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_npromz'

    get_npromz = self%npromz

  END FUNCTION get_npromz

  INTEGER FUNCTION get_blk_start(grid)

    USE mo_jsb_domain_iface, ONLY: get_blk_start_domain => get_blk_start

    CLASS(t_jsb_grid), INTENT(in) :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_blk_start'

    get_blk_start = get_blk_start_domain(grid%patch)

  END FUNCTION get_blk_start

  INTEGER FUNCTION get_blk_end(grid)

    USE mo_jsb_domain_iface, ONLY: get_blk_end_domain => get_blk_end

    CLASS(t_jsb_grid), INTENT(in) :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_blk_end'

    get_blk_end = get_blk_end_domain(grid%patch)

  END FUNCTION get_blk_end

  INTEGER FUNCTION get_col_start(grid, iblk)

    USE mo_jsb_domain_iface, ONLY: get_col_start_domain => get_col_start

    CLASS(t_jsb_grid), INTENT(in) :: grid
    INTEGER,          INTENT(in) :: iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_col_start'

    get_col_start = get_col_start_domain(iblk, grid%patch)

  END FUNCTION get_col_start

  INTEGER FUNCTION get_col_end(grid, iblk)

    USE mo_jsb_domain_iface, ONLY: get_col_end_domain => get_col_end

    CLASS(t_jsb_grid), INTENT(in) :: grid
    INTEGER,          INTENT(in) :: iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_col_end'

    get_col_end = get_col_end_domain(iblk, grid%patch)

  END FUNCTION get_col_end

#endif
END MODULE mo_jsb_grid_class
