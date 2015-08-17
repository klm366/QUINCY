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
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_jsb_domain_iface,   ONLY: t_patch, p_patch, pe_decomposed,               &
                                   get_grid_filename_domain => get_grid_filename, &
                                   get_ntotal_domain        => get_ntotal,        &
                                   get_ntotal_g_domain      => get_ntotal_g,      &
                                   get_dims_g_domain        => get_dims_g,        &
                                   get_nproma_domain        => get_nproma,        &
                                   get_nblks_domain         => get_nblks,         &
                                   get_lon_domain           => get_lon,           &
                                   get_lat_domain           => get_lat,           &
                                   get_area_domain          => get_area

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_grid, t_jsb_vgrid, t_jsb_vgrid_ptr
  PUBLIC :: new_grid, delete_grid, print_grid, new_vgrid, delete_vgrid
  PUBLIC :: get_nblks, get_nproma

  TYPE t_jsb_grid
    INTEGER                      :: id = 0
    CHARACTER(LEN=10)            :: name = ''
    CHARACTER(LEN=filename_max)  :: filename = ''
    INTEGER                      :: ntotal = 0          !< Number of cells on this PE
    INTEGER                      :: nland = 0
    REAL(wp), POINTER            :: lon(:,:) => NULL()
    REAL(wp), POINTER            :: lat(:,:) => NULL()
    REAL(wp), POINTER            :: area(:,:) => NULL()
    LOGICAL,  POINTER            :: lsm(:,:) => NULL()
    REAL(wp), POINTER            :: lsf(:,:) => NULL()
    INTEGER                      :: nblks = 1
    INTEGER                      :: nproma = 0
    INTEGER                      :: npromz = 0
    INTEGER                      :: ntotal_g            !< Total number of cells in grid
    INTEGER                      :: dims_g(2)            !< Dimensions of global grid
    LOGICAL                      :: configured = .FALSE.
    TYPE(t_patch), POINTER       :: patch
    TYPE(pe_decomposed), POINTER :: decomposition
  END TYPE t_jsb_grid

  TYPE t_jsb_vgrid
    INTEGER ::           id           = 0
    CHARACTER(len=20) :: name         = ''
    CHARACTER(len=50) :: longname     = ''
    CHARACTER(len=10) :: units        = ''
    INTEGER :: n_levels               = 0
    REAL(wp), POINTER :: levels(:)    => NULL()
    REAL(wp), POINTER :: lbounds(:)   => NULL()
    REAL(wp), POINTER :: ubounds(:)   => NULL()
    INTEGER :: cdiZaxisID             = 0
  END TYPE t_jsb_vgrid

  TYPE t_jsb_vgrid_ptr
    TYPE(t_jsb_vgrid), POINTER :: p
  END TYPE t_jsb_vgrid_ptr

  INTERFACE new_grid
    MODULE PROCEDURE new_grid_from_file
    MODULE PROCEDURE new_grid_from_host
    MODULE PROCEDURE new_grid_from_lonlat
  END INTERFACE new_grid

  INTEGER, SAVE :: ngrid_count = 0

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_grid_class'

CONTAINS

  !>
  !! Grid constructor

  FUNCTION new_grid_from_host(patch, decomp) RESULT(grid)

    TYPE(t_patch),       OPTIONAL, INTENT(in) :: patch          !< used with icon
    TYPE(pe_decomposed), OPTIONAL, INTENT(in) :: decomp         !< used with echam
    TYPE(t_jsb_grid), POINTER                 :: grid

    TYPE(t_jsb_grid), POINTER       :: new

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_grid_from_host'

    CALL message(TRIM(routine), 'starting construction of new JSBACH grid instance from host model')

    IF (PRESENT(patch) .AND. PRESENT(decomp)) &
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')

    IF (.NOT. PRESENT(patch) .AND. .NOT. PRESENT(decomp)) &
      CALL finish(TRIM(routine), 'Either patch or decomposition must be specified')

    ALLOCATE(new)

    ngrid_count = ngrid_count + 1
    new%id = ngrid_count
!!$    new%name = name

    IF (PRESENT(patch)) THEN
      new%patch => p_patch(patch%id)
      ALLOCATE(new%decomposition)
      new%decomposition%nlon = -1
    ELSE
      ALLOCATE(new%decomposition)
      new%decomposition = decomp
      ALLOCATE(new%patch)
      new%patch%id = -1
    END IF

    new%filename = get_grid_filename_domain(patch, decomp)
    new%ntotal   = get_ntotal_domain       (patch, decomp)
    new%ntotal_g = get_ntotal_g_domain     (patch, decomp)
    new%dims_g(1:2) = get_dims_g_domain    (patch, decomp)
    new%nblks    = get_nblks_domain        (patch, decomp)
    new%nproma   = get_nproma_domain       (patch, decomp)
    new%npromz   = new%ntotal - new%nproma * (new%nblks - 1)
    new%lon     => get_lon_domain          (patch, decomp)
    new%lat     => get_lat_domain          (patch, decomp)
    new%area    => get_area_domain         (patch, decomp)
    ALLOCATE(new%lsm(new%nproma,new%nblks), new%lsf(new%nproma,new%nblks))
    new%lsm(:,:) = .TRUE.
    new%lsf(:,:) = 1.0_wp

!!$    ALLOCATE(new%patch)
!!$    new%patch%id = patch%id
!!$    new%patch%n_patch_cells = patch%n_patch_cells
!!$    new%patch%n_patch_cells_g = patch%n_patch_cells_g
!!$    new%patch%nblks_c = patch%nblks_c
!!$    new%patch%npromz_c = patch%npromz_c
!!$    new%patch%grid_filename = patch%grid_filename
!!$    new%patch%comm_pat_gather_c = patch%comm_pat_gather_c
!!$    new%patch%cells%glb_index = patch%cells%glb_index


    new%configured = .TRUE.

    grid => new

    CALL message(TRIM(routine), 'construction of JSBACH grid completed.')

  END FUNCTION new_grid_from_host

  FUNCTION new_grid_from_file(filename, name) RESULT(grid)

!!$    USE mo_jsb_io_netcdf,      ONLY: ncfile_has_dim
    USE mo_jsb_parallel_iface, ONLY: my_process_is_stdio, p_io, p_bcast, mpi_comm

    CHARACTER(LEN=*), INTENT(in)    :: filename
    CHARACTER(LEN=*), INTENT(in)    :: name
    TYPE(t_jsb_grid), POINTER       :: grid

    TYPE(t_jsb_grid), POINTER       :: new

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_grid_from_file'

    CALL message(TRIM(routine), 'starting construction of new JSBACH grid instance from file:')
    CALL message('     ', TRIM(filename))

    ALLOCATE(new)

    ngrid_count = ngrid_count + 1
    new%id = ngrid_count
    new%name = name
    new%filename = filename

!!$    IF (my_process_is_stdio()) THEN
!!$
!!$      IF (ncfile_has_dim(TRIM(filename), 'lon_cell_centre')) THEN
!!$        ! ICON grid
!!$      ELSE
!!$        ! ECHAM grid
!!$      END IF
!!$
!!$    END IF

    !CALL p_bcast(, p_io, mpi_comm)

    new%configured = .TRUE.

    grid => new

    CALL message(TRIM(routine), 'construction of JSBACH grid completed.')

  END FUNCTION new_grid_from_file

  FUNCTION new_grid_from_lonlat(lon, lat) RESULT(grid)

    REAL(wp)                        :: lon
    REAL(wp)                        :: lat
    TYPE(t_jsb_grid), POINTER       :: grid

    TYPE(t_jsb_grid), POINTER       :: new

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_grid_from_lonlat'

    ALLOCATE(new)

    grid => new

  END FUNCTION new_grid_from_lonlat

  SUBROUTINE delete_grid(self)

#ifdef HAVE_F2003
    TYPE(t_jsb_grid), POINTER, INTENT(inout) :: self
#else
    TYPE(t_jsb_grid), POINTER :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':delete_grid'

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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':print_grid'

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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_nblks'

    get_nblks = self%nblks

  END FUNCTION get_nblks

  INTEGER FUNCTION get_nproma(self)

#ifdef HAVE_F2003
    TYPE(t_jsb_grid), POINTER, INTENT(in) :: self
#else
    TYPE(t_jsb_grid), POINTER             :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_nproma'

    get_nproma = self%nproma

  END FUNCTION get_nproma

  INTEGER FUNCTION get_blk_start(grid)

    USE mo_jsb_domain_iface, ONLY: get_blk_start_domain => get_blk_start

    TYPE(t_jsb_grid), INTENT(in) :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_blk_start'

    IF (grid%patch%id >= 0) THEN
      get_blk_start = get_blk_start_domain(patch=grid%patch)
    ELSE IF (grid%decomposition%nlon >= 0) THEN
      get_blk_start = get_blk_start_domain(decomp=grid%decomposition)
    ELSE
      CALL finish(TRIM(routine), 'Something went wrong.')
    END IF

  END FUNCTION get_blk_start

  INTEGER FUNCTION get_blk_end(grid)

    USE mo_jsb_domain_iface, ONLY: get_blk_end_domain => get_blk_end

    TYPE(t_jsb_grid), INTENT(in) :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_blk_end'

    IF (grid%patch%id >= 0) THEN
      get_blk_end = get_blk_end_domain(patch=grid%patch)
    ELSE IF (grid%decomposition%nlon >= 0) THEN
      get_blk_end = get_blk_end_domain(decomp=grid%decomposition)
    ELSE
      CALL finish(TRIM(routine), 'Something went wrong.')
    END IF

  END FUNCTION get_blk_end

  INTEGER FUNCTION get_col_start(iblk, grid)

    USE mo_jsb_domain_iface, ONLY: get_col_start_domain => get_col_start

    INTEGER,          INTENT(in) :: iblk
    TYPE(t_jsb_grid), INTENT(in) :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_col_start'

    IF (grid%patch%id >= 0) THEN
      get_col_start = get_col_start_domain(iblk, patch=grid%patch)
    ELSE IF (grid%decomposition%nlon >= 0) THEN
      get_col_start = get_col_start_domain(iblk, decomp=grid%decomposition)
    ELSE
      CALL finish(TRIM(routine), 'Something went wrong.')
    END IF

  END FUNCTION get_col_start

  INTEGER FUNCTION get_col_end(iblk, grid)

    USE mo_jsb_domain_iface, ONLY: get_col_end_domain => get_col_end

    INTEGER,          INTENT(in) :: iblk
    TYPE(t_jsb_grid), INTENT(in) :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_col_end'

    IF (grid%patch%id >= 0) THEN
      get_col_end = get_col_end_domain(iblk, patch=grid%patch)
    ELSE IF (grid%decomposition%nlon >= 0) THEN
      get_col_end = get_col_end_domain(iblk, decomp=grid%decomposition)
    ELSE
      CALL finish(TRIM(routine), 'Something went wrong.')
    END IF

  END FUNCTION get_col_end

  FUNCTION new_vgrid(name, cdi_axis_type, length, longname, units, levels, lbounds, ubounds) RESULT(new)

    USE mo_jsb_io,       ONLY: cdi_zaxis_ids, zaxisCreate, zaxisInqName, zaxisInqSize
    USE mo_jsb_io_iface, ONLY: zaxisDefName, zaxisDefLevels, zaxisDefLbounds,          &
                               zaxisDefUbounds, zaxisDefUnits, zaxisDefLongname

    CHARACTER(len=*),           INTENT(in) :: name
    INTEGER,                    INTENT(in) :: cdi_axis_type
    INTEGER,                    INTENT(in) :: length
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: longname
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: units
    REAL(wp),         OPTIONAL, INTENT(in) :: levels(:), lbounds(:), ubounds(:)

    TYPE(t_jsb_vgrid), POINTER :: new

    INTEGER :: cdi_size
    CHARACTER(len=100) :: cdi_name
    INTEGER :: i, idx

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_vgrid'

    IF(length < 1) CALL finish(TRIM(routine), TRIM(name)//' - length < 1')

    ALLOCATE(new)

    new%name = TRIM(name)
    new%n_levels = length
    ALLOCATE(new%levels(length), new%lbounds(length), new%ubounds(length))
    IF (PRESENT(levels)) THEN
      new%levels(:) = levels(1:length)
    ELSE
      DO i=1,length
        new%levels(:) = REAL(i,wp)
      END DO
    END IF
    IF (PRESENT(lbounds)) THEN
      new%lbounds(:) = lbounds(1:length)
    ELSE
      new%lbounds(1) = 0._wp
      DO i=2,length
        new%lbounds(i) = new%levels(i-1) + (new%levels(i-1) - new%lbounds(i-1))
      END DO
    END IF
    IF (PRESENT(ubounds)) THEN
      new%ubounds(:) = ubounds(1:length)
    ELSE
      DO i=1,length
        new%ubounds(i) = new%levels(i) + (new%levels(i) - new%lbounds(i))
      END DO
    END IF
    IF (PRESENT(longname)) new%longname = TRIM(longname)
    IF (PRESENT(units)) new%units  = TRIM(units)

    ! Find free slot in cdi_zaxis_ids
    idx = 0
    DO i=1,SIZE(cdi_zaxis_ids)
      IF (cdi_zaxis_ids(i) == 0) THEN
        idx = i
        EXIT
      END IF
    END DO
    IF (idx == 0) CALL finish(TRIM(routine), 'Increase size of cdi_zaxis_ids')
    new%id = idx

    ! Find CDI zaxis with same name and same size
    idx = 0
    DO i=1,SIZE(cdi_zaxis_ids)
      IF (cdi_zaxis_ids(i) == 0) CYCLE
      CALL zaxisInqName(cdi_zaxis_ids(i), cdi_name)
      IF (TRIM(cdi_name) /= name) THEN
        CYCLE
      ELSE
        cdi_size = zaxisInqSize(cdi_zaxis_ids(i))
        IF (cdi_size /= length) THEN
          CYCLE
        ELSE
          idx = i
          EXIT
        END IF
      END IF
    END DO
    ! If none found create a new CDI zaxis, else use existing one
    IF (idx == 0) THEN
      new%cdiZaxisID = zaxisCreate(cdi_axis_type, length)
      CALL zaxisDefName    (new%cdiZaxisID, TRIM(new%name))
      CALL zaxisDefLevels  (new%cdiZaxisID, new%levels(1:length))
      CALL zaxisDefLbounds (new%cdiZaxisID, new%lbounds(1:length))
      CALL zaxisDefUbounds (new%cdiZaxisID, new%ubounds(1:length))
      CALL zaxisDefUnits   (new%cdiZaxisID, TRIM(new%units))
      CALL zaxisDefLongname(new%cdiZaxisID, TRIM(new%longname))
      WRITE(message_text, '(A,A,2I3)') 'Created new vgrid: name, ID, length= ',TRIM(new%name), new%cdiZaxisID, length
      CALL message(TRIM(routine), message_text)
    ELSE
      new%cdiZaxisID = cdi_zaxis_ids(idx)
    END IF
    cdi_zaxis_ids(new%id) = new%cdiZaxisID

  END FUNCTION new_vgrid

  SUBROUTINE delete_vgrid(self)

    USE mo_jsb_io, ONLY: cdi_zaxis_ids, zaxisDestroy

#ifdef HAVE_F2003
    TYPE(t_jsb_vgrid), POINTER, INTENT(inout) :: self
#else
    TYPE(t_jsb_vgrid), POINTER :: self
#endif

    INTEGER :: i, idx

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':delete_vgrid'

    IF (self%id == 0) RETURN

    CALL message(TRIM(routine), 'starting destruction of JSBACH vgrid instance.')

    IF (ASSOCIATED(self%levels)) DEALLOCATE(self%levels)
    IF (ASSOCIATED(self%lbounds)) DEALLOCATE(self%lbounds)
    IF (ASSOCIATED(self%ubounds)) DEALLOCATE(self%ubounds)

    cdi_zaxis_ids(self%id) = 0
    ! If this was the last instance of a JSBACH vgrid that uses a particular cdiZaxisID, also destroy CDI zaxis
    idx = 0
    DO i=1,SIZE(cdi_zaxis_ids)
      IF (cdi_zaxis_ids(i) == self%cdiZaxisID) THEN
        idx = i
        EXIT
      END IF
    END DO
    IF (idx == 0) CALL zaxisDestroy(self%cdiZaxisID)

    DEALLOCATE(self)
!!$    self%cdiZaxisID = 0
!!$    self%id = 0
!!$    self%n_levels = 0
!!$    self%name = ''

    CALL message(TRIM(routine), 'destruction of JSBACH vgrid completed.')

  END SUBROUTINE delete_vgrid

#endif
END MODULE mo_jsb_grid_class
