!>
!! @brief Contains interface to spatial domains from ICON/ECHAM6
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_domain_iface
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
#ifdef __ICON__
  USE mo_model_domain,       ONLY: t_patch, p_patch
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_impl_constants,     ONLY: min_rlcell_int
  USE mo_loopindices,        ONLY: get_indices_c
#else
  USE mo_decomposition,      ONLY: pe_decomposed, global_decomposition
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_patch, p_patch, pe_decomposed
  PUBLIC :: get_grid_filename, get_ntotal, get_ntotal_g, get_dims_g, get_nproma, get_nblks, get_lon, get_lat, get_area
  PUBLIC :: get_blk_start, get_blk_end, get_col_start, get_col_end

#ifdef __ICON__
  ! Dummy type for ECHAM decomposition
  TYPE pe_decomposed
    INTEGER :: nlon = -1
  END TYPE pe_decomposed

  TYPE(pe_decomposed), SAVE :: global_decomposition(1)
#else
  ! Dummy type for ICON patch
  TYPE t_patch
     INTEGER :: id = -1
  END type t_patch

  TYPE(t_patch), TARGET, SAVE :: p_patch(1)
#endif

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_domain_iface'

CONTAINS

  FUNCTION get_grid_filename(patch, decomp) RESULT(filename)

    USE mo_io_units, ONLY: filename_max

    TYPE(t_patch),       OPTIONAL, INTENT(in) :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in) :: decomp
    CHARACTER(len=filename_max)               :: filename

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_grid_filename'

    IF (PRESENT(patch) .AND. PRESENT(decomp)) THEN
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')
    ELSE IF(PRESENT(patch)) THEN
#ifdef __ICON__
      filename = patch%grid_filename
#else
      CALL finish(TRIM(routine), 'patch not allowed')
#endif
    ELSE IF(PRESENT(decomp)) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'decomposition not allowed')
#else
      filename = ''
#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_grid_filename

  INTEGER FUNCTION get_ntotal(patch, decomp)

    TYPE(t_patch),       OPTIONAL, INTENT(in) :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in) :: decomp

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_ntotal'

    IF (PRESENT(patch) .AND. PRESENT(decomp)) THEN
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')
    ELSE IF(PRESENT(patch)) THEN
#ifdef __ICON__
      get_ntotal = patch%n_patch_cells
#else
      CALL finish(TRIM(routine), 'patch not allowed')
#endif
    ELSE IF(PRESENT(decomp)) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'decomposition not allowed')
#else
      get_ntotal = decomp%nglon * decomp%nglat
#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_ntotal

  INTEGER FUNCTION get_ntotal_g(patch, decomp)

    TYPE(t_patch),       OPTIONAL, INTENT(in) :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in) :: decomp

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_ntotal_g'

    IF (PRESENT(patch) .AND. PRESENT(decomp)) THEN
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')
    ELSE IF(PRESENT(patch)) THEN
#ifdef __ICON__
      get_ntotal_g = patch%n_patch_cells_g
#else
      CALL finish(TRIM(routine), 'patch not allowed')
#endif
    ELSE IF(PRESENT(decomp)) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'decomposition not allowed')
#else
      get_ntotal_g = decomp%nlon * decomp%nlat
#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_ntotal_g

  FUNCTION get_dims_g(patch, decomp)

    TYPE(t_patch),       OPTIONAL, INTENT(in) :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in) :: decomp
    INTEGER                                   :: get_dims_g(2)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_dims_g'

    IF (PRESENT(patch) .AND. PRESENT(decomp)) THEN
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')
    ELSE IF(PRESENT(patch)) THEN
#ifdef __ICON__
      get_dims_g(1) = patch%n_patch_cells_g
      get_dims_g(2) = 1
#else
      CALL finish(TRIM(routine), 'patch not allowed')
#endif
    ELSE IF(PRESENT(decomp)) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'decomposition not allowed')
#else
      get_dims_g(1) = decomp%nlon
      get_dims_g(2) = decomp%nlat
#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_dims_g

  INTEGER FUNCTION get_nblks(patch, decomp)

    TYPE(t_patch),       OPTIONAL, INTENT(in) :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in) :: decomp

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_nblks'

    IF (PRESENT(patch) .AND. PRESENT(decomp)) THEN
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')
    ELSE IF(PRESENT(patch)) THEN
#ifdef __ICON__
      get_nblks = patch%nblks_c
#else
      CALL finish(TRIM(routine), 'patch not allowed')
#endif
    ELSE IF(PRESENT(decomp)) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'decomposition not allowed')
#else
      get_nblks = decomp%ngpblks
#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_nblks

  INTEGER FUNCTION get_nproma(patch, decomp)

    TYPE(t_patch),       OPTIONAL, INTENT(in) :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in) :: decomp

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_nproma'

    IF (PRESENT(patch) .AND. PRESENT(decomp)) THEN
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')
    ELSE IF(PRESENT(patch)) THEN
#ifdef __ICON__
      get_nproma = INT( (patch%n_patch_cells - patch%npromz_c)/(patch%nblks_c-1))
#else
      CALL finish(TRIM(routine), 'patch not allowed')
#endif
    ELSE IF(PRESENT(decomp)) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'decomposition not allowed')
#else
      get_nproma = decomp%nproma
#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_nproma

  INTEGER FUNCTION get_blk_start(patch, decomp)

    TYPE(t_patch),       OPTIONAL, INTENT(in) :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in) :: decomp

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_blk_start'

    IF (PRESENT(patch) .AND. PRESENT(decomp)) THEN
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')
    ELSE IF(PRESENT(patch)) THEN
#ifdef __ICON__
      get_blk_start = patch%cells%start_blk(grf_bdywidth_c + 1, 1)
#else
      CALL finish(TRIM(routine), 'patch not allowed')
#endif
    ELSE IF(PRESENT(decomp)) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'decomposition not allowed')
#else
      get_blk_start = 1
#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_blk_start

  INTEGER FUNCTION get_blk_end(patch, decomp)

    TYPE(t_patch),       OPTIONAL, INTENT(in) :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in) :: decomp

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_blk_end'

    IF (PRESENT(patch) .AND. PRESENT(decomp)) THEN
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')
    ELSE IF(PRESENT(patch)) THEN
#ifdef __ICON__
      get_blk_end = patch%cells%end_blk(min_rlcell_int, MAX(1, patch%n_childdom))
#else
      CALL finish(TRIM(routine), 'patch not allowed')
#endif
    ELSE IF(PRESENT(decomp)) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'decomposition not allowed')
#else
      get_blk_end = decomp%ngpblks
#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_blk_end

  INTEGER FUNCTION get_col_start(jb, patch, decomp)

    INTEGER,                       INTENT(in) :: jb
    TYPE(t_patch),       OPTIONAL, INTENT(in) :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in) :: decomp

    INTEGER :: i_startblk, i_endblk, jcs, jce, rl_start, rl_end

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_col_start'

    IF (PRESENT(patch) .AND. PRESENT(decomp)) THEN
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')
    ELSE IF(PRESENT(patch)) THEN
#ifdef __ICON__
      i_startblk = get_blk_start(patch)
      i_endblk   = get_blk_end(patch)
      rl_start   = grf_bdywidth_c + 1
      rl_end     = min_rlcell_int

      CALL get_indices_c(patch, jb, i_startblk, i_endblk, jcs, jce, rl_start, rl_end)

      get_col_start = jcs
#else
      CALL finish(TRIM(routine), 'patch not allowed')
#endif
    ELSE IF(PRESENT(decomp)) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'decomp not allowed')
#else
      get_col_start = 1
#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_col_start

  INTEGER FUNCTION get_col_end(jb, patch, decomp)

    INTEGER,                       INTENT(in) :: jb
    TYPE(t_patch),       OPTIONAL, INTENT(in) :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in) :: decomp

    INTEGER :: i_startblk, i_endblk, jcs, jce, rl_start, rl_end

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_col_end'

    IF (PRESENT(patch) .AND. PRESENT(decomp)) THEN
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')
    ELSE IF(PRESENT(patch)) THEN
#ifdef __ICON__
      i_startblk = get_blk_start(patch)
      i_endblk   = get_blk_end(patch)
      rl_start   = grf_bdywidth_c + 1
      rl_end     = min_rlcell_int

      CALL get_indices_c(patch, jb, i_startblk, i_endblk, jcs, jce, rl_start, rl_end)

      get_col_end = jce
#else
      CALL finish(TRIM(routine), 'patch not allowed')
#endif
    ELSE IF(PRESENT(decomp)) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'decomp not allowed')
#else
      IF (jb == decomp%ngpblks) THEN
        get_col_end = decomp%npromz
      ELSE
        get_col_end = decomp%nproma
     END IF
#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_col_end

  FUNCTION get_lon(patch, decomp) RESULT(lon)

#ifndef __ICON__
    USE mo_geoloc, ONLY: philon_2d
#endif

    TYPE(t_patch),       OPTIONAL, INTENT(in)  :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in)  :: decomp
    REAL(wp), POINTER                          :: lon(:,:)

    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_lon'

    nproma = get_nproma(patch, decomp)
    nblks  = get_nblks (patch, decomp)
    ALLOCATE(lon(nproma,nblks))

    IF (PRESENT(patch) .AND. PRESENT(decomp)) THEN
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')
    ELSE IF(PRESENT(patch)) THEN
#ifdef __ICON__
      lon(:,:) = patch%cells%center(:,:)%lon
#else
      CALL finish(TRIM(routine), 'patch not allowed')
#endif
    ELSE IF(PRESENT(decomp)) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'decomposition not allowed')
#else
      IF (.NOT. ASSOCIATED(philon_2d)) &
        CALL finish(TRIM(routine), 'Longitudes of decomposition not computed yet')
      lon(:,:) = philon_2d(:,:)
#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_lon

  FUNCTION get_lat(patch, decomp) RESULT(lat)

#ifndef __ICON__
    USE mo_geoloc, ONLY: philat_2d
#endif

    TYPE(t_patch),       OPTIONAL, INTENT(in)  :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in)  :: decomp
    REAL(wp), POINTER                          :: lat(:,:)

    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_lat'

    nproma = get_nproma(patch, decomp)
    nblks  = get_nblks (patch, decomp)
    ALLOCATE(lat(nproma,nblks))

    IF (PRESENT(patch) .AND. PRESENT(decomp)) THEN
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')
    ELSE IF(PRESENT(patch)) THEN
#ifdef __ICON__
      lat(:,:) = patch%cells%center(:,:)%lat
#else
      CALL finish(TRIM(routine), 'patch not allowed')
#endif
    ELSE IF(PRESENT(decomp)) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'decomposition not allowed')
#else
      IF (.NOT. ASSOCIATED(philat_2d)) &
        CALL finish(TRIM(routine), 'Latitudes of decomposition not computed yet')
      lat(:,:) = philat_2d(:,:)
#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_lat

  FUNCTION get_area(patch, decomp) RESULT(area)

#ifndef __ICON__
    USE mo_geoloc, ONLY: gboxarea_2d
#endif

    TYPE(t_patch),       OPTIONAL, INTENT(in)  :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in)  :: decomp
    REAL(wp), POINTER                          :: area(:,:)

    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_area'

    nproma = get_nproma(patch, decomp)
    nblks  = get_nblks (patch, decomp)
    ALLOCATE(area(nproma,nblks))

    IF (PRESENT(patch) .AND. PRESENT(decomp)) THEN
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')
    ELSE IF(PRESENT(patch)) THEN
#ifdef __ICON__
      area(:,:) = patch%cells%area(:,:)
#else
      CALL finish(TRIM(routine), 'patch not allowed')
#endif
    ELSE IF(PRESENT(decomp)) THEN
#ifdef __ICON__
      CALL finish(TRIM(routine), 'decomposition not allowed')
#else
      IF (.NOT. ASSOCIATED(gboxarea_2d)) &
        CALL finish(TRIM(routine), 'grid box surface area of decomposition not computed yet')
      area(:,:) = gboxarea_2d(:,:)
#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_area

#endif
END MODULE mo_jsb_domain_iface
