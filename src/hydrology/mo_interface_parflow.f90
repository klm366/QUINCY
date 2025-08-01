!> Interface between JSBACH hydrology and ParFlow using YAC
!>
!> ICON-Land
!>
!> Jan Weinkaemmerer (FZ Juelich)
!>
MODULE mo_interface_parflow
#if !defined(__NO_JSBACH__)

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_coupling_config,     ONLY: is_coupled_run
#ifdef YAC_coupling
  USE mo_model_domain,        ONLY: t_patch
  USE yac,                    ONLY: yac_fdef_field,           &
       &                            yac_fput, yac_fget,       &
       &                            YAC_TIME_UNIT_ISO_FORMAT, &
       &                            YAC_ACTION_PUT_FOR_RESTART
  USE mo_time_config,         ONLY: time_config
  USE mtime,                  ONLY: timedeltaToString, MAX_TIMEDELTA_STR_LEN
#endif
  USE mo_jsb_task_class,      ONLY: t_jsb_task_options
  USE mo_jsb_control,         ONLY: acc_stream

  dsl4jsb_Use_processes HYDRO_
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

#ifdef YAC_coupling
  PUBLIC :: construct_parflow_coupling
#endif
  PUBLIC :: interface_parflow
  PUBLIC :: field_id
  PUBLIC :: received_soil_water, is_valid_parflow, parflow_buffer

  CHARACTER(len=*), PARAMETER :: modname = 'mo_interface_parflow'

  INTEGER :: no_of_fields = 3

  LOGICAL :: received_fract_parflow = .FALSE.

  INTEGER, ALLOCATABLE :: field_id(:)
  LOGICAL, ALLOCATABLE, TARGET :: received_soil_water(:,:,:)
  REAL(wp), ALLOCATABLE, TARGET :: is_valid_parflow(:,:)
  REAL(wp), ALLOCATABLE, TARGET :: parflow_buffer(:,:,:)

CONTAINS

#ifdef YAC_coupling
  SUBROUTINE construct_parflow_coupling(p_patch, comp_id, cell_point_id, grid_id)

    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: Get_model
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid

    TYPE(t_patch), TARGET, INTENT(in) :: p_patch(:)
    INTEGER, INTENT(in)               :: comp_id           ! yac component id
    INTEGER, INTENT(in)               :: cell_point_id     ! yac point is
    INTEGER, INTENT(in)               :: grid_id           ! yac grid id
    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_grid),  POINTER    :: grid
    TYPE(t_jsb_vgrid), POINTER    :: soil_w

    CLASS(t_jsb_tile_abstract), POINTER :: tile

    TYPE(t_patch), POINTER :: patch_horz

    INTEGER :: model_id
    INTEGER :: jg, jb, jc
    INTEGER :: nproma!, nblks
    INTEGER :: nsoil

    CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: modelTimeStep

    CHARACTER(len=40), ALLOCATABLE :: field_name(:)
    INTEGER :: collection_size(no_of_fields)

    dsl4jsb_Def_config(HYDRO_)

    ! It is assumed that only one domain is active (jg=1) i.e. no nesting.
    ! Currently, this requirement is fulfilled since ocean coupling is only set up for one domain.

    jg = 1
    patch_horz => p_patch(jg)

    model_id = 1
    model => Get_model(model_id)

    dsl4jsb_Get_config(HYDRO_)
    IF (.NOT. dsl4jsb_Config(HYDRO_)%active) RETURN

    CALL model%Get_top_tile(tile)

    grid => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    !nblks  = grid%Get_nblks()

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    CALL timedeltaToString(time_config%tc_dt_model, modelTimeStep)

    ALLOCATE(field_name(no_of_fields))
    ALLOCATE(field_id(no_of_fields))

    field_name(1)             = "water_to_parflow"
    collection_size(1)        = nsoil

    field_name(2)             = "fract_parflow"
    collection_size(2)        = 1

    field_name(3)             = "soil_water"
    collection_size(3)        = nsoil

    ALLOCATE(received_soil_water(nproma,nsoil,patch_horz%nblks_c))
    ALLOCATE(is_valid_parflow(nproma,patch_horz%nblks_c))
    ALLOCATE(parflow_buffer(nproma,nsoil,patch_horz%nblks_c))

    received_soil_water(:,:,:) = .FALSE.
    is_valid_parflow(:,:) = 0.0_wp

    DO jc = 1, no_of_fields
      CALL yac_fdef_field (            &
        & TRIM(field_name(jc)),        &
        & comp_id,                     &
        & cell_point_id,               &
        & 1,                           &
        & collection_size(jc),         &
        & modelTimeStep,               &
        & YAC_TIME_UNIT_ISO_FORMAT,    &
        & field_id(jc) )
    ENDDO

  END SUBROUTINE construct_parflow_coupling
#endif

  SUBROUTINE interface_parflow(tile, options)

    USE mo_fortran_tools, ONLY: init

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
    INTEGER               :: iblk, n, nn, nlen
    INTEGER               :: nsoil, is
    REAL(wp), POINTER     :: dz(:)

    n = options%nc ! only to avoid compiler warnings

    IF ( .NOT. is_coupled_run() ) RETURN

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(HYDRO_)
    IF (.NOT. dsl4jsb_Config(HYDRO_)%active) RETURN

    !IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'HD model works on root tile only!')

    grid => Get_grid(model%grid_id)
    nbr_hor_cells = grid%ntotal

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels
    dz => soil_w%dz(:)

    ! Get reference to variables in memory
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onDomain(HYDRO_, water_to_parflow)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, evapotrans_to_parflow)

    ALLOCATE(buffer(grid%nproma * grid%nblks, nsoil))
    !$ACC DATA CREATE(buffer)

    ! Send water_to_parflow
    ! --------------------

!$OMP PARALLEL
    CALL init(buffer(:,:), lacc=.TRUE.)
!$OMP END PARALLEL

    write_coupler_restart = .FALSE.

!$OMP PARALLEL DO PRIVATE(iblk, nn, nlen, n)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO iblk = 1, grid%nblks
      nn = (iblk - 1) * grid%nproma
      IF (iblk /= grid%nblks) THEN
        nlen = grid%nproma
      ELSE
        nlen = grid%npromz
      END IF
      !$ACC LOOP GANG VECTOR
      DO n = 1, nlen
        buffer(nn+n,1) = ( water_to_parflow(n,iblk) + evapotrans_to_parflow(n,1,iblk) ) * 3600.0_wp / dz(1)
      END DO
    END DO
    !$ACC END PARALLEL
    !$ACC UPDATE HOST(buffer(:,1)) ASYNC(acc_stream)
    !$ACC WAIT(acc_stream)
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(iblk, nn, nlen, is, n)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP SEQ
    DO iblk = 1, grid%nblks
      nn = (iblk - 1) * grid%nproma
      IF (iblk /= grid%nblks) THEN
        nlen = grid%nproma
      ELSE
        nlen = grid%npromz
      END IF
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO is = 2, nsoil
        DO n = 1, nlen
          buffer(nn+n,is) = evapotrans_to_parflow(n,is,iblk) * 3600.0_wp / dz(is)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL
    !$ACC UPDATE HOST(buffer(:,2:nsoil)) ASYNC(acc_stream)
    !$ACC WAIT(acc_stream)
!$OMP END PARALLEL DO

#ifdef YAC_coupling
    CALL yac_fput(field_id(1), nbr_hor_cells, nsoil, buffer(1:nbr_hor_cells,1:nsoil), info, ierror)
    IF ( info == YAC_ACTION_PUT_FOR_RESTART ) write_coupler_restart = .TRUE.
#endif

    IF ( write_coupler_restart ) THEN
       CALL message(TRIM(routine), 'YAC says it is put for restart')
    END IF

    ! Receive fract_parflow
    ! --------------------

!$OMP PARALLEL
    CALL init(buffer(:,:), lacc=.TRUE.)
!$OMP END PARALLEL

#ifdef YAC_coupling
    IF ( .NOT. received_fract_parflow ) THEN
      CALL yac_fget(field_id(2), nbr_hor_cells, 1, buffer(1:nbr_hor_cells,1:1), info, ierror)
      IF ( info > 0 ) THEN
        CALL message(TRIM(routine), 'Received fract_parflow')
      END IF

      IF ( info > 0 ) THEN
!$OMP PARALLEL DO PRIVATE(iblk, nn, nlen, n)
        !$ACC UPDATE DEVICE(buffer(1:nbr_hor_cells,1:1)) ASYNC(acc_stream)
        !$ACC DATA COPYOUT(is_valid_parflow)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream) NO_CREATE(is_valid_parflow)
        !$ACC LOOP SEQ
        DO iblk = 1, grid%nblks
          nn = (iblk - 1) * grid%nproma
          IF (iblk /= grid%nblks) THEN
            nlen = grid%nproma
          ELSE
            nlen = grid%npromz
          END IF
          !$ACC LOOP GANG VECTOR
          DO n = 1, nlen
            is_valid_parflow(n,iblk) = buffer(nn+n,1)
          END DO
        END DO
        !$ACC END PARALLEL
        !$ACC WAIT(acc_stream)
        !$ACC END DATA
!$OMP END PARALLEL DO

        received_fract_parflow = .TRUE.
      END IF
    END IF
#endif

    ! Receive soil_water
    ! --------------------

!$OMP PARALLEL
    CALL init(buffer(:,:), lacc=.TRUE.)
!$OMP END PARALLEL

#ifdef YAC_coupling
    CALL yac_fget(field_id(3), nbr_hor_cells, nsoil, buffer(1:nbr_hor_cells,1:nsoil), info, ierror)

    IF ( info > 0 ) THEN
!$OMP PARALLEL DO PRIVATE(iblk, nn, nlen, is, n)
      !$ACC UPDATE DEVICE(buffer(1:nbr_hor_cells,1:nsoil)) ASYNC(acc_stream)
      !$ACC DATA COPYOUT(parflow_buffer)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream) NO_CREATE(parflow_buffer)
      !$ACC LOOP SEQ
      DO iblk = 1, grid%nblks
        nn = (iblk - 1) * grid%nproma
        IF (iblk /= grid%nblks) THEN
          nlen = grid%nproma
        ELSE
          nlen = grid%npromz
        END IF
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO is = 1, nsoil
          DO n = 1, nlen
            parflow_buffer(n,is,iblk) = buffer(nn+n,is)
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
      !$ACC WAIT(acc_stream)
      !$ACC END DATA
!$OMP END PARALLEL DO

      received_soil_water(:,:,:) = .TRUE.
    END IF
#endif

    !$ACC END DATA
    DEALLOCATE(buffer)

  END SUBROUTINE interface_parflow

#endif
END MODULE mo_interface_parflow
