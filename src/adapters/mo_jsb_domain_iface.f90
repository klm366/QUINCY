!>
!! @brief Contains interface to spatial domains from ICON/ECHAM6
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! JN (2016-02-20) new subroutines and some adaptations for jsb4 standalone
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
  PUBLIC :: get_nlon_nlat_from_file, jsbalone_init_decomp

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

#ifdef __ICON__
    USE mo_parallel_config,    ONLY: nproma
#endif

    TYPE(t_patch),       OPTIONAL, INTENT(in) :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in) :: decomp

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_nproma'

    IF (PRESENT(patch) .AND. PRESENT(decomp)) THEN
      CALL finish(TRIM(routine), 'Only one of patch or decomposition can be specified')
    ELSE IF(PRESENT(patch)) THEN
#ifdef __ICON__
      get_nproma = nproma
!!$      IF (patch%nblks_c > 1) THEN
!!$        get_nproma = INT( (patch%n_patch_cells - patch%npromz_c)/(patch%nblks_c-1))
!!$      ELSE
!!$        get_nproma = patch%n_patch_cells
!!$      END IF
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

#ifdef STANDALONE
    USE mo_control,   ONLY: nlon
    USE mo_gaussgrid,  ONLY: philon
    USE mo_transpose,     ONLY: reorder
    USE mo_decomposition, ONLY: ldc => local_decomposition
#else
    USE mo_geoloc,     ONLY: philon_2d
#endif

#endif

    TYPE(t_patch),       OPTIONAL, INTENT(in)  :: patch
    TYPE(pe_decomposed), OPTIONAL, INTENT(in)  :: decomp
    REAL(wp), POINTER                          :: lon(:,:)

    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_lon'

#if !defined(__ICON__) && defined(STANDALONE)
    INTEGER :: jlat, jlon
    REAL(wp) :: tmp (ldc%nglon, ldc%nglat)
#endif


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

#ifdef STANDALONE
      ! in jsb standalone init_geoloc is never called, therefore philon_2d is never calculated.
      ! -> alternatively lons could be read from a file (grid_filename) 
      !     - this solution would probably be more suitable for a driver for ECHAM and ICON?
      DO jlat = 1, ldc%nglat     ! local  latitude index N -> S
        DO jlon = 1, ldc%nglon   ! local  longitude index
          tmp(jlon,jlat) = philon (jlon + ldc%glon(jlat))
        END DO
      END DO

      CALL reorder (lon, tmp)
#else
      IF (.NOT. ASSOCIATED(philon_2d)) &
        CALL finish(TRIM(routine), 'Longitudes of decomposition not computed yet')
      lon(:,:) = philon_2d(:,:)
#endif

#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_lon

  FUNCTION get_lat(patch, decomp) RESULT(lat)

#ifndef __ICON__

#ifdef STANDALONE
    USE mo_gaussgrid,  ONLY: philat
    USE mo_transpose,     ONLY: reorder
    USE mo_decomposition, ONLY: ldc => local_decomposition
#else
    USE mo_geoloc, ONLY: philat_2d
#endif

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

#ifdef STANDALONE
      ! in jsb standalone init_geoloc is never called, therefore philat_2d is never calculated.
      ! -> alternatively lats could be read from a file (grid_filename) 
      !     - this solution would probably be more suitable for a driver for ECHAM and ICON?
      CALL reorder (lat, SPREAD(philat(ldc%glat(:)),1 ,ldc%nglon))
#else
      IF (.NOT. ASSOCIATED(philat_2d)) &
        CALL finish(TRIM(routine), 'Latitudes of decomposition not computed yet')
      lat(:,:) = philat_2d(:,:)
#endif

#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_lat

  FUNCTION get_area(patch, decomp) RESULT(area)

#ifndef __ICON__

#ifdef STANDALONE
    USE mo_transpose,     ONLY: reorder
    USE mo_decomposition, ONLY: ldc => local_decomposition
    USE mo_gaussgrid,     ONLY: gridarea
#else
    USE mo_geoloc, ONLY: gboxarea_2d
#endif

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

#ifdef STANDALONE
      ! in jsb standalone init_geoloc is never called, therefore gboxarea_2d has not been calculated.
      CALL reorder (area ,SPREAD(gridarea(ldc%glat(:)), 1, ldc%nglon))
#else
      IF (.NOT. ASSOCIATED(gboxarea_2d)) &
        CALL finish(TRIM(routine), 'grid box surface area of decomposition not computed yet')
      area(:,:) = gboxarea_2d(:,:)
#endif

#endif
    ELSE
      CALL finish(TRIM(routine), 'Need to specify either patch or decomposition')
    END IF

  END FUNCTION get_area

  !>
  !! @brief Read number of lons and lats from a given netcdf file
  !!
  !! Currently only working for ECHAM
  !!
  !! @par History
  !! JN: Initial parts copied from jsbach3 mo_jsbach_grid.f90:init_grid
  !!
  SUBROUTINE get_nlon_nlat_from_file(grid_filename, read_nlon, read_nlat) !filename
  
    USE mo_jsb_parallel_iface, ONLY: my_process_is_stdio, p_io, p_bcast, my_process_is_mpi_parallel

#ifndef __ICON__  
    USE mo_control,   ONLY: nlon, ngl, nhgl, nvclev, vct
    USE mo_io,     ONLY: IO_READ, IO_open, IO_close
    USE mo_netCDF, ONLY: file_info, io_inq_dimid, io_inq_dimlen
#endif

    CHARACTER(len=*), INTENT(in) :: grid_filename
    INTEGER, INTENT(out) :: read_nlon, read_nlat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_nlon_nlat_from_file'

#ifndef __ICON__  
    TYPE(file_info) :: grid_file
    INTEGER :: IO_dim_id
#endif


#ifdef __ICON__  
    CALL finish(TRIM(routine), 'not implemented for use with ICON')
#else

    IF (my_process_is_stdio()) THEN
      grid_file%opened = .FALSE.

      CALL io_open(grid_filename, grid_file, IO_READ)
      CALL io_inq_dimid  (grid_file%file_id, 'lon', IO_dim_id)
      CALL io_inq_dimlen (grid_file%file_id, IO_dim_id, read_nlon)
      CALL io_inq_dimid  (grid_file%file_id, 'lat', IO_dim_id)
      CALL io_inq_dimlen (grid_file%file_id, IO_dim_id, read_nlat)

      CALL IO_close(grid_file)
    ENDIF

    IF (my_process_is_mpi_parallel()) THEN
      CALL p_bcast (read_nlon, p_io)
      CALL p_bcast (read_nlat, p_io)
    ENDIF

    !Required for standalone model for mo_gaussgrid::inigau()
    nlon = read_nlon
    ngl = read_nlat
    nhgl = ngl / 2
    nvclev = 0

#endif
  END SUBROUTINE get_nlon_nlat_from_file


  !>
  !! @brief Get a decomposition for jsbach standalone
  !!
  !! Currently only working for ECHAM
  !!
  !! @par History
  !! JN: Initial version copied from jsbach3 jsbalone_init_decomp.f90
  !!
  SUBROUTINE jsbalone_init_decomp (nlon, nlat, nproca, nprocb, npedim, decomposition)

#ifndef __ICON__  
    USE mo_decomposition, ONLY : decompose, local_decomposition, global_decomposition
    USE mo_mpi,           ONLY : p_nprocs, p_parallel_io, p_pe
    USE mo_transpose,     ONLY : indx
    USE mo_io_units,      ONLY : nout
#endif

    INTEGER, INTENT(in) :: nlon, nlat, nproca, nprocb, npedim
    TYPE(pe_decomposed), INTENT(OUT) :: decomposition
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':jsbalone_init_decomp'


#ifndef __ICON__  
    INTEGER :: debug_parallel, i, p
    LOGICAL :: lrot
#endif


#ifdef __ICON__  
    CALL finish(TRIM(routine), 'not implemented for use with ICON')
#else
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
#endif
  END SUBROUTINE jsbalone_init_decomp

#endif
END MODULE mo_jsb_domain_iface
