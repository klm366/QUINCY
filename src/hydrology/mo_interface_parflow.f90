!>
!! @brief Interface between JSBACH hydrology and ParFlow using YAC
!!
!! @author Jan Weinkaemmerer (FZ Juelich)
!!
!! @par Revision History
!! 

!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_interface_parflow
#if !defined(__NO_JSBACH__)

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
#ifdef YAC_coupling
  USE mo_model_domain,        ONLY: t_patch
  USE mo_mpi,                 ONLY: p_pe_work, p_comm_work, p_sum
  USE mo_parallel_config,     ONLY: nproma
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, iaes
  USE mo_master_control,      ONLY: get_my_process_name
  USE mo_yac_finterface,      ONLY: yac_fget_version, yac_fdef_comp, &
       &                            yac_fdef_datetime, yac_fdef_grid, &
       &                            yac_fdef_points, yac_fset_global_index, &
       &                            yac_fset_core_mask, yac_fdef_mask,      &
       &                            yac_fdef_field_mask, yac_fenddef,       &
       &                            yac_fput, yac_fget, &
       &                            YAC_LOCATION_CELL,        &
       &                            YAC_TIME_UNIT_ISO_FORMAT, &
       &                            YAC_ACTION_PUT_FOR_RESTART
  USE mo_run_config,          ONLY: iforcing
  USE mo_time_config,         ONLY: time_config
  USE mtime,                  ONLY: datetimeToString, MAX_DATETIME_STR_LEN, &
       &                            timedeltaToString, MAX_TIMEDELTA_STR_LEN
#endif

  USE mo_jsb_task_class,      ONLY: t_jsb_task_options

  dsl4jsb_Use_processes HYDRO_
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

#ifdef YAC_coupling
  PUBLIC :: construct_parflow_coupling
#endif
  PUBLIC :: interface_parflow
  PUBLIC :: lyac_very_1st_get, nbr_inner_cells, mask_checksum, field_id, &
       & is_valid_parflow, parflow_buffer

  CHARACTER(len=*), PARAMETER :: modname = 'mo_interface_parflow'

  INTEGER :: no_of_fields = 2

  INTEGER, ALLOCATABLE :: field_id(:)
  LOGICAL, ALLOCATABLE, TARGET :: is_valid_parflow(:,:)
  REAL(wp), ALLOCATABLE, TARGET :: parflow_buffer(:,:,:)

  INTEGER, SAVE :: nbr_inner_cells
  INTEGER, SAVE :: mask_checksum = -1
  LOGICAL, SAVE :: lyac_very_1st_get = .TRUE.

CONTAINS

#ifdef YAC_coupling
  SUBROUTINE construct_parflow_coupling(model_id)

    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: Get_model
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid

    INTEGER, INTENT(in) :: model_id
    
    ! Local variables

    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_grid),  POINTER    :: grid
    TYPE(t_jsb_vgrid), POINTER    :: soil_w

    CLASS(t_jsb_tile_abstract), POINTER :: tile

    TYPE(t_patch), POINTER :: patch_horz

    CHARACTER(LEN=max_char_length) :: grid_name
    CHARACTER(LEN=max_char_length) :: comp_name

    INTEGER :: comp_id           ! yac component id
    INTEGER :: comp_ids(1)
    INTEGER :: cell_point_ids(1) ! yac point ids
    INTEGER :: grid_id           ! yac grid id

    INTEGER :: cell_mask_ids(1)
    INTEGER :: jb, jc, nn
    INTEGER :: nblks
    INTEGER :: nsoil
    INTEGER :: nbr_vertices_per_cell

    REAL(wp), ALLOCATABLE :: buffer_lon(:)
    REAL(wp), ALLOCATABLE :: buffer_lat(:)
    INTEGER,  ALLOCATABLE :: buffer_c(:,:)

    LOGICAL,  ALLOCATABLE :: is_valid(:)

    REAL(wp), ALLOCATABLE :: lsmnolake(:,:)

    REAL(wp), PARAMETER :: eps = 1.E-10_wp

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: startdatestring
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: stopdatestring
    CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: modelTimeStep

    CHARACTER(len=40), ALLOCATABLE :: field_name(:)
    INTEGER :: collection_size(no_of_fields)

    dsl4jsb_Def_config(HYDRO_)

    model => Get_model(model_id)

    dsl4jsb_Get_config(HYDRO_)
    IF (.NOT. dsl4jsb_Config(HYDRO_)%active) RETURN

    CALL model%Get_top_tile(tile)

    grid => Get_grid(model%grid_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    ! YAC initialization steps

    patch_horz => grid%patch

    ! Print the YAC version
    CALL message('Running ICON atmosphere in coupled mode with YAC version ', TRIM(yac_fget_version()) )

    ! Overwrite job start and end date with component data
    CALL datetimeToString(time_config%tc_startdate, startdatestring)
    CALL datetimeToString(time_config%tc_stopdate, stopdatestring)
    CALL timedeltaToString(time_config%tc_dt_model, modelTimeStep)

    CALL yac_fdef_datetime ( start_datetime = TRIM(startdatestring), &
         &                   end_datetime   = TRIM(stopdatestring)   )

    ! Inform the coupler about what we are
    comp_name = TRIM(get_my_process_name())
    CALL yac_fdef_comp ( TRIM(comp_name), comp_id )
    comp_ids(1) = comp_id

    ! Announce one grid (patch) to the coupler
    grid_name = "icon_atmos_grid"

    ! Extract cell information

    nblks = max(patch_horz%nblks_c,patch_horz%nblks_v)

    ALLOCATE(buffer_lon(nproma*nblks))
    ALLOCATE(buffer_lat(nproma*nblks))
    ALLOCATE(buffer_c(3,nproma*nblks))

    ALLOCATE(lsmnolake(nproma,nblks))

    nbr_vertices_per_cell = 3

    !ICON_OMP_PARALLEL
    !ICON_OMP_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, patch_horz%nblks_v
      DO jc = 1, nproma
        nn = (jb-1)*nproma+jc
        buffer_lon(nn) = patch_horz%verts%vertex(jc,jb)%lon
        buffer_lat(nn) = patch_horz%verts%vertex(jc,jb)%lat
      ENDDO
    ENDDO
    !ICON_OMP_END_DO NOWAIT

    !ICON_OMP_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, patch_horz%nblks_c
      DO jc = 1, nproma
        nn = (jb-1)*nproma+jc
        buffer_c(1,nn) = (patch_horz%cells%vertex_blk(jc,jb,1)-1)*nproma + &
          &               patch_horz%cells%vertex_idx(jc,jb,1)
        buffer_c(2,nn) = (patch_horz%cells%vertex_blk(jc,jb,2)-1)*nproma + &
          &               patch_horz%cells%vertex_idx(jc,jb,2)
        buffer_c(3,nn) = (patch_horz%cells%vertex_blk(jc,jb,3)-1)*nproma + &
                          patch_horz%cells%vertex_idx(jc,jb,3)
      ENDDO
    ENDDO
    !ICON_OMP_END_DO
    !ICON_OMP_END_PARALLEL

    ! Description of elements, here as unstructured grid
    CALL yac_fdef_grid(           &
      & TRIM(grid_name),          &
      & patch_horz%n_patch_verts, &
      & patch_horz%n_patch_cells, &
      & nbr_vertices_per_cell,    &
      & buffer_lon,               &
      & buffer_lat,               &
      & buffer_c,                 &
      & grid_id)

    ! Define cell center points (location = 0)

    !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, patch_horz%nblks_c
      DO jc = 1, nproma
        nn = (jb-1)*nproma+jc
        buffer_lon(nn) = patch_horz%cells%center(jc,jb)%lon
        buffer_lat(nn) = patch_horz%cells%center(jc,jb)%lat
      ENDDO
    ENDDO
    !ICON_OMP_END_PARALLEL_DO

    ! center points in cells (needed e.g. for patch recovery and nearest neighbour interpolation)
    CALL yac_fdef_points (        &
      & grid_id,                  &
      & patch_horz%n_patch_cells, &
      & YAC_LOCATION_CELL,        &
      & buffer_lon,               &
      & buffer_lat,               &
      & cell_point_ids(1) )

    DEALLOCATE (buffer_lon, buffer_lat, buffer_c)

    CALL yac_fset_global_index (                &
      & patch_horz%cells%decomp_info%glb_index, &
      & YAC_LOCATION_CELL,                      &
      & grid_id )

    ALLOCATE(is_valid(nproma*patch_horz%nblks_c))
    ALLOCATE(is_valid_parflow(nproma,patch_horz%nblks_c))

    nbr_inner_cells = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jc) REDUCTION(+:nbr_inner_cells) ICON_OMP_RUNTIME_SCHEDULE
    DO jc = 1, patch_horz%n_patch_cells
       IF ( p_pe_work == patch_horz%cells%decomp_info%owner_local(jc) ) THEN
         is_valid(jc) = .TRUE.
         nbr_inner_cells = nbr_inner_cells + 1
       ELSE
         is_valid(jc) = .FALSE.
       ENDIF
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    CALL yac_fset_core_mask ( &
      & is_valid,             &
      & YAC_LOCATION_CELL,    &
      & grid_id )

    ALLOCATE(field_name(no_of_fields))
    ALLOCATE(field_id(no_of_fields))
    ALLOCATE(parflow_buffer(nproma,nsoil,patch_horz%nblks_c))

    field_name(1)             = "water_to_parflow"
    collection_size(1)        = nsoil

    field_name(2)             = "soil_water"
    collection_size(2)        = nsoil

    ! Define the mask for YAC

    SELECT CASE ( iforcing )

      CASE ( iaes )
#ifdef __NO_AES__
        CALL finish ('mo_interface_parflow:construct_parflow_coupling', &
            & 'coupled model needs aes; remove --disable-aes and reconfigure')
#else
        !ICON_OMP_PARALLEL_DO PRIVATE(jb,jc) ICON_OMP_RUNTIME_SCHEDULE
        DO jb = 1, patch_horz%nblks_c
            DO jc = 1, nproma
              lsmnolake(jc, jb) = grid%lsf(jc,jb)
            ENDDO
        ENDDO
        !ICON_OMP_END_PARALLEL_DO

        mask_checksum = 0
        !ICON_OMP_PARALLEL_DO PRIVATE(jb,jc) REDUCTION(+:mask_checksum) ICON_OMP_RUNTIME_SCHEDULE
        DO jb = 1, patch_horz%nblks_c
          DO jc = 1, nproma
            mask_checksum = mask_checksum + ABS( lsmnolake(jc,jb) )
          ENDDO
        ENDDO
        !ICON_OMP_END_PARALLEL_DO

        mask_checksum = p_sum(mask_checksum, comm=p_comm_work)

        !
        ! Define cell_mask_ids(1): all ocean and coastal points are valid
        ! This is the standard for the coupling of atmospheric fields listed below
        !
        IF ( mask_checksum > 0 ) THEN
          !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
          DO jb = 1, patch_horz%nblks_c
            DO jc = 1, nproma

              !IF ( lsmnolake(jc, jb) .LT. (1.0_wp - eps) &
              !  & .OR. ABS( grid%lat(jc,jb) ) .GT. model%config%parflow_lat ) THEN
              IF ( lsmnolake(jc, jb) .LT. (1.0_wp - eps) ) THEN
                ! ocean point (fraction of ocean is >0., lsmnolake .lt. 1.) is valid
                ! or grid point outside the Parflow latitude belt
                is_valid((jb-1)*nproma+jc) = .FALSE.
                is_valid_parflow(jc,jb) = .FALSE.
              ELSE
                ! land point (fraction of land is one, no sea water, lsmnolake=1.) is undef
                is_valid((jb-1)*nproma+jc) = .TRUE.
                is_valid_parflow(jc,jb) = .TRUE.
              ENDIF

            ENDDO
          ENDDO
          !ICON_OMP_END_PARALLEL_DO
        ELSE
          !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
          DO jb = 1, patch_horz%nblks_c
            DO jc = 1, nproma
              is_valid((jb-1)*nproma+jc) = .FALSE.
              is_valid_parflow(jc,jb) = .FALSE.
            ENDDO
          ENDDO
          !ICON_OMP_END_PARALLEL_DO
        ENDIF
#endif
      CASE DEFAULT

        CALL finish ('Please mask handling for new forcing.')

    END SELECT

    ! DO jb = 1, patch_horz%nblks_c
    !   DO jc = 1, nproma
    !     is_valid((jb-1)*nproma+jc) = .TRUE.
    !     is_valid_parflow(jc,jb) = .TRUE.
    !   ENDDO
    ! ENDDO

    CALL yac_fdef_mask (          &
      & grid_id,                  &
      & patch_horz%n_patch_cells, &
      & YAC_LOCATION_CELL,        &
      & is_valid,                 &
      & cell_mask_ids(1) )

    CALL yac_fsync_def ()

    DO jc = 1, no_of_fields
      CALL yac_fdef_field_mask (       &
        & TRIM(field_name(jc)),        &
        & comp_id,                     &
        & cell_point_ids,              &
        & cell_mask_ids(1),            & 
        & 1,                           &
        & collection_size(jc),         &
        & modelTimeStep,               &
        & YAC_TIME_UNIT_ISO_FORMAT,    &
        & field_id(jc) )
    ENDDO

    ! End of YAC initialization steps

    DEALLOCATE (is_valid)

    DEALLOCATE (lsmnolake)

    ! End definition of coupling fields and search

    CALL yac_fenddef ()

  END SUBROUTINE construct_parflow_coupling
#endif

  SUBROUTINE interface_parflow(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile
    TYPE(t_jsb_task_options),   INTENT(in) :: options

    ! Local variables

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_w

    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)

    CHARACTER(len=*),  PARAMETER :: routine = modname//':interface_parflow'

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: water_to_parflow
    dsl4jsb_Real3D_onDomain :: evapotrans_to_parflow

    ! Local variables

    LOGICAL               :: write_coupler_restart
    INTEGER               :: nbr_hor_cells
    REAL(wp), ALLOCATABLE :: buffer(:,:)
    INTEGER               :: info, ierror !< return values form cpl_put/get calls
    INTEGER               :: iblk, n, nn !, nlen
    INTEGER               :: nsoil, is

    n = options%nc ! only to avoid compiler warnings

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(HYDRO_)
    IF (.NOT. dsl4jsb_Config(HYDRO_)%active) RETURN

    ! What is that line good for?
    !IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'HD model works on root tile only!')

    grid => Get_grid(model%grid_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    ! Get reference to variables in memory
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onDomain(HYDRO_, water_to_parflow)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, evapotrans_to_parflow)

    ! Send water_to_parflow
    ! --------------------

    write_coupler_restart = .FALSE.

    ALLOCATE(buffer(grid%nproma * grid%nblks, nsoil))
    buffer(:,:) = 0.0_wp
    DO iblk = 1, grid%nblks
      nn = (iblk - 1) * grid%nproma
      !IF (iblk /= grid%nblks) THEN
      !  nlen = grid%nproma
      !ELSE
      !  nlen = grid%npromz
      !END IF
      DO n = 1, grid%nproma
        buffer(nn+n,1) = water_to_parflow(n,iblk) + evapotrans_to_parflow(n,1,iblk)
      END DO
    END DO
    DO is = 2, nsoil
      DO iblk = 1, grid%nblks
        nn = (iblk - 1) * grid%nproma
        DO n = 1, grid%nproma
          buffer(nn+n,is) = evapotrans_to_parflow(n,is,iblk)
        END DO
      END DO
    END DO

    nbr_hor_cells = grid%ntotal
#ifdef YAC_coupling
    CALL yac_fput(field_id(1), nbr_hor_cells, nsoil, buffer(1:nbr_hor_cells,1:nsoil), info, ierror)
    IF ( info == YAC_ACTION_PUT_FOR_RESTART ) write_coupler_restart = .TRUE.
#endif

    IF ( write_coupler_restart ) THEN
       CALL message(TRIM(routine), 'YAC says it is put for restart')
    END IF

    ! Receive soil_water
    ! --------------------

    buffer(:,:) = 0.0_wp

#ifdef YAC_coupling
    CALL yac_fget(field_id(2), nbr_hor_cells, nsoil, buffer(1:nbr_hor_cells,1:nsoil), info, ierror)

    IF ( info > 0 ) THEN
      DO is = 1, nsoil
        DO iblk = 1, grid%nblks
          nn = (iblk - 1) * grid%nproma
          DO n = 1, grid%nproma
            parflow_buffer(n,is,iblk) = buffer(nn+n,is)
          END DO
        END DO
      END DO

      lyac_very_1st_get = .FALSE.
    END IF
#endif

    DEALLOCATE(buffer)

  END SUBROUTINE interface_parflow

#endif
END MODULE mo_interface_parflow
