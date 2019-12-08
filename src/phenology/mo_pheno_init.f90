!>
!! @brief Initialization of the the phenology memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-09-29)
!! Julia Nabel (2018-01-10): l_forestRegrowth flag
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_pheno_init
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_process_class,   ONLY: ON_TILE_, AGGREGATE_
  USE mo_jsb_task_class,      ONLY: t_jsb_task_options

  dsl4jsb_Use_processes PHENO_
  dsl4jsb_Use_config(PHENO_)
  dsl4jsb_Use_memory(PHENO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: pheno_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_init'

CONTAINS

  SUBROUTINE pheno_init(tile)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_init'

    IF (tile%process_action(PHENO_) == ON_TILE_) THEN
      CALL pheno_init_bc(tile)
      CALL pheno_init_ic(tile)
    ELSE IF (tile%process_action(PHENO_) == AGGREGATE_) THEN
      CALL pheno_aggregate_bc(tile)
      CALL pheno_aggregate_ic(tile)
    END IF

  END SUBROUTINE pheno_init

  SUBROUTINE pheno_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model

    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)

    REAL(wp), POINTER ::  &
      & ptr_mon(:,:,:),   &
      & ptr_2D(:,:)

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_init_bc'

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)

    IF (debug_on()) CALL message(TRIM(routine), 'Setting boundary conditions of pheno memory for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(PHENO_)%bc_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(PHENO_)%bc_filename), model%grid_id)

    ! Leaf area index climatology
    IF (debug_on()) CALL message(TRIM(routine), 'reading lai_clim ...')
    ptr_mon => input_file%Read_2d_time(    &
      & variable_name='lai_clim')
    ptr_mon(:,:,:) = MERGE(ptr_mon(:,:,:), 0._wp, ptr_mon(:,:,:) > 0._wp)
    dsl4jsb_memory(PHENO_)%lai_mon(:,:,1:12) = ptr_mon(:,:,1:12)
    dsl4jsb_memory(PHENO_)%lai_mon(:,:,0   ) = ptr_mon(:,:,12  )
    dsl4jsb_memory(PHENO_)%lai_mon(:,:,13  ) = ptr_mon(:,:,1   )
    DEALLOCATE(ptr_mon)
    !$acc enter data attach(dsl4jsb_memory(PHENO_)%lai_mon)
    !$acc update device(dsl4jsb_memory(PHENO_)%lai_mon)

    ! Foliage projected cover climatology
    IF (debug_on()) CALL message(TRIM(routine), 'reading veg_fract ...')
    ptr_mon => input_file%Read_2d_time(    &
      & variable_name='veg_fract')
    ptr_mon(:,:,:) = MERGE(ptr_mon(:,:,:), 0._wp, ptr_mon(:,:,:) > 0._wp)
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,1:12) = ptr_mon(:,:,1:12)
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,0   ) = ptr_mon(:,:,12  )
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,13  ) = ptr_mon(:,:,1   )
    DEALLOCATE(ptr_mon)
    !$acc update device(dsl4jsb_memory(PHENO_)%fract_fpc_mon)

    ! Forest fraction
    ! @todo: only used for use_alb_veg_simple=.TRUE.
    IF (debug_on()) CALL message(TRIM(routine), 'reading forest_fract ...')
    ptr_2D => input_file%Read_2d(      &
      & variable_name='forest_fract',  &
      & fill_array = dsl4jsb_var_ptr(PHENO_, fract_forest))
    dsl4jsb_var2D_onDomain(PHENO_, fract_forest) = &
      MERGE(dsl4jsb_var2D_onDomain(PHENO_, fract_forest), 0._wp, dsl4jsb_var2D_onDomain(PHENO_, fract_forest) >= 0._wp)
    !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, fract_forest))

    CALL input_file%Close()

    input_file = jsb_netcdf_open_input(TRIM(model%config%fract_filename), model%grid_id)     ! Temporary

    ! Maximum foliage projected cover
    IF (debug_on()) CALL message(TRIM(routine), 'reading veg_ratio_max ...')
    ptr_2D => input_file%Read_2d(      &
      & variable_name='veg_ratio_max', &
      & fill_array = dsl4jsb_var_ptr(PHENO_, fract_fpc_max))

    dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max) = &
      MERGE(dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max), 0._wp, dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max) >= 0._wp)
    !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max))

    CALL input_file%Close()

    ! avoid compiler warnings about pointers not being dereferenced
    IF (ASSOCIATED(ptr_2D)) CONTINUE
    NULLIFY(ptr_2D)

    IF (tile%lcts(1)%lib_id /= 0)  THEN ! only if the present tile is a pft
      !JN 10.01.18: for now: set the maxLai for non forest types to the one from the Lctlib
      !             later one might think about allometric relation ships in other multi annual plants, too?
      IF ( (dsl4jsb_Config(PHENO_)%l_forestRegrowth) .AND. (.NOT. dsl4jsb_Lctlib_param(ForestFlag)) ) THEN
        dsl4jsb_var2D_onDomain(PHENO_, maxLAI_allom) = dsl4jsb_Lctlib_param(MaxLAI)
        !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, maxLAI_allom))
      END IF

      !JN 10.01.18: with l_forestRegrowth the maxLAI can change over time for forest pfts
      IF ( (dsl4jsb_Config(PHENO_)%l_forestRegrowth) .AND. (dsl4jsb_Lctlib_param(ForestFlag)) ) THEN
        dsl4jsb_var2D_onDomain(PHENO_, veg_fract_correction) =   &
          & 1.0_wp - exp(-dsl4jsb_var2D_onDomain(PHENO_, maxLAI_allom) / dsl4jsb_Lctlib_param(clumpinessFactor))
      !@todo: do not rely on the lib_id but use a flag instead?
      ELSE
        dsl4jsb_var2D_onDomain(PHENO_, veg_fract_correction) =   &
          & 1.0_wp - exp(-dsl4jsb_Lctlib_param(MaxLAI) / dsl4jsb_Lctlib_param(clumpinessFactor))
      END IF
      !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, veg_fract_correction))
    END IF

  END SUBROUTINE pheno_init_bc

  SUBROUTINE pheno_aggregate_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_memory(PHENO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract
    TYPE(t_jsb_model),       POINTER :: model
    TYPE(t_jsb_grid),        POINTER :: grid
    TYPE(t_jsb_task_options)         :: options

    INTEGER :: iblk, ics, ice, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_aggregate_bc'

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    options%ics = 1
    options%ice = grid%nproma

    ics = options%ics
    ice = options%ice
    nc = ice - ics + 1

    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(PHENO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    DO iblk = 1, grid%nblks

      CALL model%Set_options(iblk=iblk)
      options%iblk = iblk
      IF (iblk == grid%nblks) THEN
        options%ice = grid%npromz
        ice = grid%npromz
        ! TODO: ics and ice should be retrieved for each iblk (are not necessarily constant for each block in ICON)
      END IF

      dsl4jsb_Aggregate_onChunk(PHENO_, fract_fpc_max,      weighted_by_fract)

    END DO

    IF (debug_on()) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE pheno_aggregate_bc

  SUBROUTINE pheno_init_ic(tile)

    USE mo_pheno_process, ONLY: get_foliage_projected_cover
    USE mo_jsb_time,      ONLY: get_time_interpolation_weights

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model

    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)

    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_init_ic'

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)

    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of pheno memory (lai) for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(PHENO_)%bc_filename))

    ! Set initial leaf area index from climatology (read in pheno_init_bc)
    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
    dsl4jsb_var2D_onDomain(PHENO_, lai) = &
      & wgt1 * dsl4jsb_memory(PHENO_)%lai_mon(:,:,nmw1) + wgt2 * dsl4jsb_memory(PHENO_)%lai_mon(:,:,nmw2)

    IF (model%config%l_compat401) THEN
      !$claw expand
      dsl4jsb_var2D_onDomain(PHENO_, fract_fpc) = get_foliage_projected_cover( &
        & dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max), dsl4jsb_var2D_onDomain(PHENO_, lai), 2._wp)
      !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, fract_fpc))
    ELSE
      IF (       TRIM(dsl4jsb_Config(PHENO_)%scheme) == 'climatology' &
          & .OR. tile%lcts(1)%lib_id == 0) THEN         ! lib_id=0: eneral vegetation tile, no specific PFT
        dsl4jsb_var2D_onDomain(PHENO_, fract_fpc) = &
          & wgt1 * dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,nmw1) + wgt2 * dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,nmw2)
      ELSE
        IF ( .NOT. dsl4jsb_Config(PHENO_)%l_forestRegrowth )  THEN
          dsl4jsb_var2D_onDomain(PHENO_, lai) = MIN(dsl4jsb_var2D_onDomain(PHENO_, lai), dsl4jsb_Lctlib_param(MaxLai))
        ELSE
          dsl4jsb_var2D_onDomain(PHENO_, lai) = &
            & MIN(dsl4jsb_var2D_onDomain(PHENO_, lai), dsl4jsb_var2D_onDomain(PHENO_, maxLAI_allom))
        END IF
        !$claw expand
        dsl4jsb_var2D_onDomain(PHENO_, fract_fpc) = get_foliage_projected_cover(                &
          & dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max), dsl4jsb_var2D_onDomain(PHENO_, lai), &
          & dsl4jsb_Lctlib_param(ClumpinessFactor))
      END IF
      !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, fract_fpc))
    END IF

    !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, lai))

  END SUBROUTINE pheno_init_ic

  SUBROUTINE pheno_aggregate_ic(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_memory(PHENO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract
    TYPE(t_jsb_model),       POINTER :: model
    TYPE(t_jsb_grid),        POINTER :: grid
    TYPE(t_jsb_task_options)         :: options

    INTEGER :: iblk, ics, ice, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_aggregate_ic'

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    options%ics = 1
    options%ice = grid%nproma

    ics = options%ics
    ice = options%ice
    nc = ice - ics + 1

    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(PHENO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    DO iblk = 1, grid%nblks

      CALL model%Set_options(iblk=iblk)
      options%iblk = iblk
      IF (iblk == grid%nblks) THEN
        options%ice = grid%npromz
        ice = grid%npromz
        ! TODO: ics and ice should be retrieved for each iblk (are not necessarily constant for each block in ICON)
      END IF

      dsl4jsb_Aggregate_onChunk(PHENO_, lai,          weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(PHENO_, fract_fpc,    weighted_by_fract)

    END DO

    IF (debug_on()) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE pheno_aggregate_ic

#endif
END MODULE mo_pheno_init
