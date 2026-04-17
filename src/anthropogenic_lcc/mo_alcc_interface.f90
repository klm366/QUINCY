!> Interface to the alcc process (anthropogenic land cover change)
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Contains the interfaces to the alcc process, determining area changes and moving to-be-conserved matter
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_alcc_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish, message_text

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
!  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes ALCC_, PPLCC_, FAGE_, NLCC_

  ! Use of process configurations
  dsl4jsb_Use_config(ALCC_)

  ! Use of process memories
  dsl4jsb_Use_memory(ALCC_)
  dsl4jsb_Use_memory(NLCC_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_alcc_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_alcc_interface'

  ! -------------------------------------------------------------------------------------------------------
  !> Type definition: task to determine ALCC cover fraction changes
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_alcc_delta_cover_fractions
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_alcc_delta_cover_fractions     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_alcc_delta_cover_fractions  !< Aggregates computed task variables
  END TYPE tsk_alcc_delta_cover_fractions

  !> Constructor interface for ALCC task deriving cover fraction changes
  INTERFACE tsk_alcc_delta_cover_fractions
    PROCEDURE Create_task_alcc_delta_cover_fractions                        !< Constructor function for task
  END INTERFACE tsk_alcc_delta_cover_fractions

  ! -------------------------------------------------------------------------------------------------------
  !> Type definition: alcc translocation task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_alcc_translocation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_alcc_translocation     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_alcc_translocation  !< Aggregates computed task variables
  END TYPE tsk_alcc_translocation

  !> Constructor interface for ALCC translocation task
  INTERFACE tsk_alcc_translocation
    PROCEDURE Create_task_alcc_translocation                        !< Constructor function for task
  END INTERFACE tsk_alcc_translocation

CONTAINS

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for ALCC task task deriving cover fraction changes
  !>
  FUNCTION Create_task_alcc_delta_cover_fractions(model_id) RESULT(return_ptr)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id    !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "alcc_delta_cover_fractions"
    ! -------------------------------------------------------------------------------------------------- !
    ALLOCATE(tsk_alcc_delta_cover_fractions::return_ptr)
    CALL return_ptr%Construct(name='alcc_delta_cover_fractions', process_id=ALCC_, owner_model_id=model_id)
  END FUNCTION Create_task_alcc_delta_cover_fractions

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for ALCC translocation task
  !>
  FUNCTION Create_task_alcc_translocation(model_id) RESULT(return_ptr)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id    !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "alcc_translocation"
    ! -------------------------------------------------------------------------------------------------- !
    ALLOCATE(tsk_alcc_translocation::return_ptr)
    CALL return_ptr%Construct(name='alcc_translocation', process_id=ALCC_, owner_model_id=model_id)
  END FUNCTION Create_task_alcc_translocation

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for ALCC process
  !>
  SUBROUTINE Register_alcc_tasks(this, model_id)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_process), INTENT(inout) :: this        !< Instance of alcc process class
    INTEGER,               INTENT(in)   :: model_id    !< Model id
    ! -------------------------------------------------------------------------------------------------- !
    CALL this%Register_task(tsk_alcc_delta_cover_fractions(model_id))
    CALL this%Register_task(tsk_alcc_translocation(model_id))

  END SUBROUTINE Register_alcc_tasks


  ! ====================================================================================================== !
  !>
  !> #### Update ALCC cover fraction changes
  !>
  !> The area fraction changes for each PFT tile are derived from prescribed anthropogenic land cover change (ALCC) forcing, read
  !> in [[mo_alcc_init:read_land_use_data]].
  !> Updating of cover fraction changes is executed at the beginning of each year.
  !> The calculation of land cover fractions and the used forcing depends on the configuration (annual or daily changes).
  !> Note: while this lcc process changes the PFT/age class tiles, it runs and is controlled on their parent tile (the VEG tile).
  !>
  SUBROUTINE update_alcc_delta_cover_fractions(tile, options)

    USE mo_util,              ONLY: one_of
    USE mo_jsb_time,          ONLY: is_newyear, get_year, get_year_length
    USE mo_jsb_lcc_class,     ONLY: t_jsb_lcc_proc

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(ALCC_)                 !< Memory of the anthropogenic land use change process
    dsl4jsb_Def_config(ALCC_)                 !< Configuration of the anthropogenic land use change process
    dsl4jsb_Real3D_onChunk :: cf_target_year  !< Target cover fractions derived from annually read land use data []
    dsl4jsb_Real3D_onChunk :: cf_delta        !< Change in cover fractions derived from cf_target_year (annual or daily change) []
    dsl4jsb_Real3D_onChunk :: transition_rate !< transition rate as read in mo_alcc_init:read_land_use_data (annual or daily change) []

    CLASS(t_jsb_tile_abstract), POINTER :: current_tile !< Pointer to the tile for which fraction changes are determined
    TYPE(t_jsb_model), POINTER          :: model        !< Current instance of the model

    REAL(wp), DIMENSION(options%nc) :: current_fract    !< Current cover fraction
    CHARACTER(len=15), DIMENSION(:), ALLOCATABLE :: trans_names

    INTEGER  :: &
      & iblk, &         !< Current block index
      & ics, &          !< Index of first cell of block
      & ice, &          !< Index of last cell of block
      & nc, &           !< Number of cells in current block
      & ic, &           !< Cell index
      & i_cf_tile, &    !< Index used for forest age classes
      & current_year, & !< Index of a PFT tile in the array with cover fractions changes
      & number_of_days, & !< Number of tiles participating in the ALCC process
      & idx             !< Index for transition rate loop
    REAL(wp) :: dtime   !< Time step length

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_alcc_delta_cover_fractions'
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    current_year  = get_year(options%current_datetime)
    number_of_days = get_year_length(current_year)

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(ALCC_)
    dsl4jsb_Get_memory(ALCC_)
    dsl4jsb_Get_var3D_onChunk(ALCC_, cf_target_year)
    dsl4jsb_Get_var3D_onChunk(ALCC_, cf_delta)
    dsl4jsb_Get_var3D_onChunk(ALCC_, transition_rate)

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(alcc_)) RETURN

    ! Assert that PPLCC process is active
    IF (.NOT. model%Is_process_enabled(PPLCC_)) THEN
      CALL finish(TRIM(routine), 'Violation of precondition: lcc processes need pplcc to be active')
    END IF

    ! Cover fraction changes only need to be calculated once at the beginning of a new year
    IF( .NOT. is_newyear(options%current_datetime,dtime)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')


    ! Assert, that this tile is the VEG tile
    IF (.NOT. tile%name .EQ. 'veg') THEN
      CALL finish(TRIM(routine), 'Violation of precondition: alcc processes is expected to run on the veg tile, instead' &
        & //' tried to run on '// trim(tile%name))
    END IF

    ! Calculate cover fraction changes
    current_tile => tile%Get_first_child_tile()
    i_cf_tile = 0

    !start case here
    SELECT CASE (TRIM(dsl4jsb_Config(ALCC_)%scheme))
    CASE ('maps')

      DO WHILE (ASSOCIATED(current_tile))
        i_cf_tile = i_cf_tile + 1

        CALL current_tile%Get_fraction(ics, ice, iblk, fract=current_fract)

        ! Derive area changes
        DO ic = 1,nc
          IF (cf_target_year(ic, i_cf_tile) < 0.0_wp) THEN
            ! This should only happen upon init with -0.1
            cf_delta(ic, i_cf_tile) = 0.0_wp
          ELSE
            cf_delta(ic, i_cf_tile) = cf_target_year(ic, i_cf_tile) - current_fract(ic)
          END IF

          IF (dsl4jsb_Config(ALCC_)%l_daily_alcc) THEN
            ! In case that area changes are daily the difference has to be distributed over the year
            ! such that the target is reached at the end of the year
            cf_delta(ic, i_cf_tile) = cf_delta(ic, i_cf_tile) / REAL(number_of_days, wp)
          END IF
        END DO

        current_tile => current_tile%Get_next_sibling_tile()
      END DO ! current_tile

    CASE('transition')
      
      IF (dsl4jsb_Config(ALCC_)%l_daily_alcc) THEN
        ! In case that area changes are daily the transition rates have to be divided by the number of days in the year
        ALLOCATE(trans_names(SIZE(dsl4jsb_Config(ALCC_)%trans_names)))
        trans_names = dsl4jsb_Config(ALCC_)%trans_names
        
        DO idx = 1, SIZE(trans_names)
          DO ic = 1,nc

            transition_rate(ic, idx) = transition_rate(ic, idx) / REAL(number_of_days, wp)
            !Todo-AKN: This only works if transition rates are given as change relative to whole area, not relative to pft area
            ! "Zinseszins Problem"

          END DO
        END DO

      END IF

    END SELECT

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_alcc_delta_cover_fractions


  ! ====================================================================================================== !
  !>
  !> #### Aggregation of variables from calculation of ALCC cover fraction change -- currently not used.
  !>
  SUBROUTINE aggregate_alcc_delta_cover_fractions(tile, options)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_alcc_delta_cover_fractions'
    ! -------------------------------------------------------------------------------------------------- !

  END SUBROUTINE aggregate_alcc_delta_cover_fractions


  ! ====================================================================================================== !
  !
  !> #### Calculate daily/annual relocation of fractions and matter
  !>
  !> This task applies the actual fraction changes as derived in [[update_alcc_delta_cover_fractions]] for each PFT/ or age class,
  !> as well as the corresponding translocation of matter from tiles loosing area to tiles gaining area. If running with natural
  !> land cover changes (NLCC) dynamic pfts are handled separately [[handle_dynamic_pfts]].
  !>
  !> Note:
  !>
  !>   - While this lcc process changes the PFT/age class tiles, it runs and is controlled on their parent tile (the VEG tile).
  !>   - The ALCC process requires that the PPLCC process is active, too (translocated matter is assumed to be relative to the box).
  !>
  SUBROUTINE update_alcc_translocation(tile, options)

    USE mo_util,              ONLY: one_of
    USE mo_jsb_time,          ONLY: is_newday, is_newyear
    USE mo_jsb_lcc_class,     ONLY: t_jsb_lcc_proc, max_tolerated_cf_mismatch
    USE mo_jsb_lcc,           ONLY: init_lcc_reloc
    USE mo_nlcc_process,      ONLY: fract_small

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(ALCC_)          !< Memory of the anthropogenic land use change process
    dsl4jsb_Def_memory(NLCC_)          !< Memory of the natural land cover change process
    dsl4jsb_Def_config(ALCC_)          !< Configuration of the anthropogenic land use change process
    dsl4jsb_Real3D_onChunk :: cf_delta        !< Derived change in cover fractions (annual or daily change) []
    dsl4jsb_Real3D_onChunk :: transition_rate !< Transition rate (annual or daily change) []
    dsl4jsb_Real3D_onChunk :: cover_fract_pot !< Potential cover fractions (without land use) from NLCC []

    CLASS(t_jsb_tile_abstract), POINTER :: current_tile    !< Pointer to the tile for which relocations are determined
    CLASS(t_jsb_tile_abstract), POINTER :: age_class_tile  !< Pointer to for age class tiles
    TYPE(t_jsb_model), POINTER          :: model           !< Current instance of the model
    TYPE(t_jsb_lcc_proc), POINTER       :: lcc_relocations !< lcc relocations instance of the ALCC process

    REAL(wp), POINTER ::   &
      & initial_area(:,:), & !< Pointer to initial area in lcc relocations instance of the ALCC process
      & lost_area(:,:),    & !< Pointer to lost area in lcc relocations instance of the ALCC process
      & gained_area(:,:)     !< Pointer to gained area in lcc relocations instance of the ALCC process

    REAL(wp), DIMENSION(options%nc) :: &
      & cf_diff,      & !< Change in cover fractions
      & current_fract,& !< Current cover fraction
      & from_cf_sum,  & !< Sum of cover fractions of all tiles losing area for a specific transition
      & to_cf_sum       !< Sum of cover fractions of all tiles gaining area for a specific transition

    INTEGER  :: &
      & iblk, &       !< Current block index
      & ics, &        !< Index of first cell of block
      & ice, &        !< Index of last cell of block
      & nc, &         !< Number of cells in current block
      & ic, &         !< Cell index
      & i_tile, &     !< Index of a tile in the lcc_relocations
      & i_ac_index, & !< Index used for forest age classes
      & i_cf_tile, &  !< Index of a PFT tile in the array with cover fractions changes
      & nr_of_tiles, & !< Number of tiles participating in the ALCC process
      & idx, &        !< Index for transition loop
      & idxe          !< End index for transition loop

    REAL(wp) :: dtime !< Time step length

    LOGICAL :: dynamic_tile !< Boolean indicating if a tile is a dynamic tile (participating in NLCC; only if NLCC is active)
    
    CHARACTER(len=15), DIMENSION(:), ALLOCATABLE :: trans_names
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_alcc_translocation'
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(ALCC_)
    dsl4jsb_Get_memory(ALCC_)
    dsl4jsb_Get_var3D_onChunk(ALCC_, cf_delta)
    dsl4jsb_Get_var3D_onChunk(ALCC_, transition_rate)
    IF (model%Is_process_enabled(NLCC_)) THEN
      dsl4jsb_Get_memory(NLCC_)
      dsl4jsb_Get_var3D_onChunk(NLCC_, cover_fract_pot)  ! in: potential cover fractions (no land use)
    END IF
    ! Todo-AKN: How to get potential cover fractions if NLCC is not active? Maybe read in potential cover 
    ! fractions as part of the ALCC process' land use data in this case?

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(alcc_)) RETURN

    ! Assert that PPLCC process is active
    IF (.NOT. model%Is_process_enabled(PPLCC_)) THEN
      CALL finish(TRIM(routine), 'Violation of precondition: lcc processes need pplcc to be active')
    END IF

    ! Redistributions are either only conducted at the start of each year or each day
    IF(.NOT. dsl4jsb_Config(ALCC_)%l_daily_alcc) THEN
      IF( .NOT. is_newyear(options%current_datetime,dtime)) RETURN
    ELSE
      IF( .NOT. is_newday(options%current_datetime,dtime)) RETURN
    END IF

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Assert, that this tile is the VEG tile
    IF (.NOT. tile%name .EQ. 'veg') THEN
      CALL finish(TRIM(routine), 'Violation of precondition: alcc processes is expected to run on the veg tile, instead' &
        & //' tried to run on '// trim(tile%name))
    END IF

    ! 1. Get lcc structure
    dsl4jsb_Get_lcc_relocations(ALCC_, lcc_relocations)
    nr_of_tiles = lcc_relocations%nr_of_tiles

    ! Init area vectors
    dsl4jsb_Get_lcc_area_matrix(initial_area)
    dsl4jsb_Get_lcc_area_matrix(lost_area)
    dsl4jsb_Get_lcc_area_matrix(gained_area)

    SELECT CASE (TRIM(dsl4jsb_Config(ALCC_)%scheme))
    CASE('transition')
      ALLOCATE(trans_names(SIZE(dsl4jsb_Config(ALCC_)%trans_names)))
      trans_names = dsl4jsb_Config(ALCC_)%trans_names
      idxe = SIZE(trans_names)
    CASE('maps')
      idxe = 1
    END SELECT
  
    DO idx = 1, idxe
      ! loop over all transition types for transitions scheme; only one iteration for 'maps' scheme
      initial_area(:,:) = 0.0_wp
      lost_area(:,:) = 0.0_wp
      gained_area(:,:) = 0.0_wp

      ! Collect current area of all involved tiles
      CALL init_lcc_reloc(lcc_relocations, options, tile, initial_area)

      ! For transition scheme, loop over tiles and calculate from_cf_sum and to_cf_sum
      SELECT CASE (TRIM(dsl4jsb_Config(ALCC_)%scheme))
      CASE('transition')
        current_tile => tile%Get_first_child_tile()
        i_tile = 0

        from_cf_sum(:) = 0.0_wp
        to_cf_sum(:) = 0.0_wp

        DO WHILE (ASSOCIATED(current_tile))
          i_tile = i_tile + 1

           ! Assert: pft without children
          IF (current_tile%Has_children()) THEN
            CALL finish(TRIM(routine), 'Violation of assertion: alcc only expects children of a lcc tile '  &
              & //'in case of the forest age class usecase, but usecase was '// TRIM(model%config%usecase)  &
              & // '. Tile with children: ' // trim(current_tile%name) // '. Please check.')
          END IF

          ! Assert: In mo_alcc_init read_land_use_data assumes that all child tiles of the veg tile
          ! are part of the lcc structure
          IF (one_of(current_tile%name, lcc_relocations%tile_names) <= 0) THEN
            CALL finish(TRIM(routine), 'Violation of assertion: child tile of veg tile(' //TRIM(current_tile%name) &
            & //') is not part of the alcc lcc-structure. Current implementation assumes that all pfts are part of it.')
          END IF

          DO ic = 1,nc
            ! Assert: all cover fractions are >= fract_small - tolerance
            !IF (initial_area(ic,i_tile) < fract_small - 1.0E-12_wp) THEN
            !  CALL finish(TRIM(routine), 'Violation of assertion: cover fraction is smaller than allowed minimum. '  &
            !    & // 'Please check initial areas and transition rates. Tile: ' // TRIM(current_tile%name))
            !END IF
       
            ! Now actually calculate from_cf_sum and to_cf_sum
            IF (dsl4jsb_Config(ALCC_)%from_class(idx) == dsl4jsb_Config(ALCC_)%pft_class(i_tile)) THEN
              from_cf_sum(ic) = from_cf_sum(ic) + max(0.0_wp, initial_area(ic, i_tile) - fract_small)
            END IF
            IF (dsl4jsb_Config(ALCC_)%to_class(idx) == dsl4jsb_Config(ALCC_)%pft_class(i_tile)) THEN
              IF (model%Is_process_enabled(NLCC_) .AND. dsl4jsb_Config(ALCC_)%to_class(idx) == 1) THEN
                ! Back-conversion rule for agricultural lands to natural scales with the difference between 
                ! potential cover fraction and current cover fraction
                to_cf_sum(ic) = to_cf_sum(ic) + max(0.0_wp, cover_fract_pot(ic, i_tile) - initial_area(ic, i_tile))
              ELSE
                ! For 'equal relative gain and loss' scale with the sum of the current cover fraction of the PFTs 
                to_cf_sum(ic) = to_cf_sum(ic) + initial_area(ic, i_tile)
              END IF
            END IF
          END DO

          current_tile => current_tile%Get_next_sibling_tile()
        END DO ! current_tile

        ! Assure some conditions:
        ! (1) transition rates <= from_cf_sum, i.e. that the transition does not remove more area than 
        !     currently available in all tiles beloning to one transition class, taking into account that
        !     a minimum fraction needs to be kept to avoid numerical issues in NLCC
        !     - Because the lost area for each tile is calculated as 
        !       transition_rate * initial_area / from_cf_sum, 
        !       this also ensures that the lost area does not exceed the initial area in each tile.
        ! (2) for back-conversion of agricultural lands to natural: transition rates <= to_cf_sum, i.e. 
        !     that the transition does not add more area than the potentially available in all natural tiles
        ! (3) if not back-conversion of agricultural lands to natural: transition_rate + to_cf_sum <= 1, 
        !     i.e. that the transition does not add more area than available in the vegetated grid box
        DO ic = 1, nc
          !transition_rate(ic, idx) = MAX(0.0_wp, MIN(transition_rate(ic, idx), from_cf_sum(ic) - n_tile_fc * fract_small))
          transition_rate(ic, idx) = MIN(transition_rate(ic, idx), from_cf_sum(ic))
          IF ((model%Is_process_enabled(NLCC_)) .AND. (dsl4jsb_Config(ALCC_)%to_class(idx) == 1)) THEN
            transition_rate(ic, idx) = MIN(transition_rate(ic, idx), to_cf_sum(ic))
          ELSE
            transition_rate(ic, idx) = MIN(transition_rate(ic, idx), 1.0_wp - to_cf_sum(ic))
          END IF
        END DO

      END SELECT

      ! 2. Calculate area changes
      current_tile => tile%Get_first_child_tile()
      i_tile = 0
      i_cf_tile = 0

      DO WHILE (ASSOCIATED(current_tile))
        i_tile = i_tile + 1
        i_cf_tile = i_cf_tile + 1

        ! Find out if there are dynamic PFT tiles, whose fractions should not be replaced by the map's fractions.
        dynamic_tile = .FALSE.
        IF (model%Is_process_enabled(NLCC_) .AND. dsl4jsb_Lctlib_param_tile(dynamic_PFT, current_tile)) THEN
          dynamic_tile = .TRUE.
        END IF

        ! Area changes need to be determined differently for forest age classes than for the pft usecase
        IF ((TRIM(model%config%usecase) == 'jsbach_forest_age_classes') .AND. (current_tile%Has_children())) THEN

          SELECT CASE (TRIM(dsl4jsb_Config(ALCC_)%scheme))
          CASE ('transition')
            CALL finish(TRIM(routine), 'Transition scheme not implemented in combination with forest age classes yet.' )
          END SELECT

          ! Assertion: age-classes are subsequent parts of the lcc structure
          age_class_tile => current_tile%Get_first_child_tile()
          i_ac_index = i_tile
          DO WHILE (ASSOCIATED(age_class_tile))
            IF (.NOT. TRIM(age_class_tile%name) .EQ. TRIM(lcc_relocations%tile_names(i_ac_index))) THEN
              CALL finish(TRIM(routine), 'Violation of assertion: age class (' //TRIM(age_class_tile%name) &
                & //') is not listed in the expected place in the alcc lcc-structure. Found ' &
                & // TRIM(lcc_relocations%tile_names(i_ac_index)) // ' instead. Please check!')
            END IF
            age_class_tile => age_class_tile%Get_next_sibling_tile()
            i_ac_index = i_ac_index + 1
          END DO

          ! Forest pfts are not part of the lcc structure (i.e. not involved tiles)
          ! and thus the initial area has not yet been automatically collected above
          CALL current_tile%Get_fraction(ics, ice, iblk, fract=current_fract)

          ! The area loss or gain of a forest pft needs to be redistributed to its children (the age classes)
          CALL redistribute_cover_fraction_changes_of_forest_pft( &
              & current_tile, options, cf_delta(:, i_cf_tile), i_tile, gained_area, lost_area)
          i_tile = i_tile + current_tile%Get_no_of_children() - 1

          ! Finally, update the fraction on the forest pft tile itself
          current_fract(:) = current_fract(:) + cf_delta(:, i_cf_tile)
          CALL current_tile%Set_fraction(ics, ice, iblk, fract=current_fract(:))
        ELSE
          ! Pft without children (either pft usecase or non forest pfts in fage usecase!)
          ! in the pft usecase i_cf_tile equals i_tile, however, this is not the case for non forest pfts in fage usecase
          ! therefore cf_delta is are used with pft index (i_cf_tile), while the lcc arrays (initial_area, lost_area and gained_area)
          ! are used with the current tile index (i_tile)

          ! Assert: pft without children
          IF (current_tile%Has_children()) THEN
            CALL finish(TRIM(routine), 'Violation of assertion: alcc only expects children of a lcc tile '  &
              & //'in case of the forest age class usecase, but usecase was '// TRIM(model%config%usecase)  &
              & // '. Tile with children: ' // trim(current_tile%name) // '. Please check.')
          END IF

          ! Assert: In mo_alcc_init read_land_use_data assumes that all child tiles of the veg tile
          ! are part of the lcc structure
          IF (one_of(current_tile%name, lcc_relocations%tile_names) <= 0) THEN
            CALL finish(TRIM(routine), 'Violation of assertion: child tile of veg tile(' //TRIM(current_tile%name) &
            & //') is not part of the alcc lcc-structure. Current implementation assumes that all pfts are part of it.')
          END IF

          SELECT CASE (TRIM(dsl4jsb_Config(ALCC_)%scheme))
          CASE ('maps')
            DO ic = 1,nc
              cf_diff(ic) = cf_delta(ic, i_cf_tile)

              ! Determine new cover fractions, while keeping within the boundaries valid for cover fractions [0,1]
              IF ((cf_diff(ic) < 0.0_wp) .AND. (ABS(cf_diff(ic)) >= initial_area(ic, i_tile))) THEN
                ! If area is lost it should not exceed the inital area
                current_fract(ic) = 0.0_wp
                cf_diff(ic) = -1.0_wp * initial_area(ic, i_tile)
              ELSE IF ((cf_diff(ic) > 0.0_wp) .AND. ((cf_diff(ic) + initial_area(ic, i_tile)) >= 1.0_wp)) THEN
                ! If area is gained the sum of initial area and gained area should not exeed 1 (total available area)
                current_fract(ic) = 1.0_wp
                cf_diff(ic) = 1.0_wp - initial_area(ic, i_tile)
              ELSE
                current_fract(ic) = initial_area(ic, i_tile) + cf_diff(ic)
              END IF

              ! Derive gained and lost area
              IF(cf_diff(ic) >= 0.0_wp) THEN
                gained_area(ic, i_tile) = cf_diff(ic)
              ELSE
                ! The lcc framework expects also lost area to be a positive value, therefore the sign is swaped here
                lost_area(ic, i_tile) = -1.0_wp * cf_diff(ic)
              END IF
            END DO

            IF (.NOT. dynamic_tile) THEN
              CALL current_tile%Set_fraction(ics, ice, iblk, fract=current_fract(:))
            END IF
          CASE ('transition')
            ! This implementation assumes that transition rates are given as change relative to the vegetated grid box, 
            ! not relative to each pft's area.
            ! It applies (1) the rule of equal relative gain and loss, and (2) the rule for back-conversion of agricultural lands.
            ! Todo-AKN: Add correction to make sure that lost area equals gained area

            DO ic = 1,nc
              
              IF (dsl4jsb_Config(ALCC_)%from_class(idx) == dsl4jsb_Config(ALCC_)%pft_class(i_tile)) THEN
                IF (from_cf_sum(ic) < 0.0_wp) THEN
                  CALL finish(TRIM(routine), 'Violation of assertion: from_cf_sum is negative. Please check transition rates' &
                  & // ' and initial areas. Tile: ' // TRIM(current_tile%name))
                END IF
                IF (from_cf_sum(ic) > 0.0_wp) THEN
                  ! Rule of equal relative loss
                  lost_area(ic, i_tile) = transition_rate(ic, idx) *  &
                      max(0.0_wp, (initial_area(ic, i_tile) - fract_small)) / from_cf_sum(ic)
                END IF
                current_fract(ic) = initial_area(ic, i_tile) - lost_area(ic, i_tile)
              ELSE IF (dsl4jsb_Config(ALCC_)%to_class(idx) == dsl4jsb_Config(ALCC_)%pft_class(i_tile)) THEN
                IF (to_cf_sum(ic) < 0.0_wp) THEN
                  CALL finish(TRIM(routine), 'Violation of assertion: to_cf_sum is negative. Please check transition rates' &
                   & // ' and initial areas. Tile: ' // TRIM(current_tile%name))
                END IF
                IF (to_cf_sum(ic) > 0.0_wp) THEN
                  IF (model%Is_process_enabled(NLCC_) .AND.dsl4jsb_Config(ALCC_)%to_class(idx) == 1) THEN
                    ! Back-conversion rule for agricultural lands to natural scales with the difference 
                    ! between potential cover fraction and current cover fraction
                    gained_area(ic, i_tile) = transition_rate(ic, idx) * &
                        max(0.0_wp, cover_fract_pot(ic, i_tile) - initial_area(ic, i_tile)) / to_cf_sum(ic)
                  ELSE
                    ! Rule of equal relative gain 
                    gained_area(ic, i_tile) = transition_rate(ic, idx) * initial_area(ic, i_tile) / to_cf_sum(ic)
                  END IF
                END IF
                ! Ensuring a minimum cover fraction should avoid this situation, but let's check.
                ! If it occured, area would not be conserved.
                IF (to_cf_sum(ic) == 0.0_wp .AND. transition_rate(ic, idx) > 0.0_wp) THEN
                  CALL finish(TRIM(routine), 'Violation of assertion: to_cf_sum is zero while transition_rate is positive.' &
                   & // ' Please check transition rates and initial areas. Tile: ' // TRIM(current_tile%name))
                END IF
                current_fract(ic) = initial_area(ic, i_tile) + gained_area(ic, i_tile)
              ELSE
                current_fract(ic) = initial_area(ic, i_tile)
              END IF

            END DO

            !IF (.NOT. dynamic_tile) THEN
            IF (ANY(current_fract(:) < 0.0_wp) .OR. ANY(current_fract(:) > 1.0_wp)) THEN
              CALL finish(TRIM(routine), 'Violation of assertion: cover fractions need to be between 0 and 1.' &
               & // ' Please check transition rates and initial areas. Tile: ' // TRIM(current_tile%name))
            END IF 
            CALL current_tile%Set_fraction(ics, ice, iblk, fract=current_fract(:))
            !END IF

          END SELECT
        END IF ! If running with jsbach_forest_age_classes usecase and a forest pft or ELSE a pft without leaves as current tile

        current_tile => current_tile%Get_next_sibling_tile()
      END DO ! Current_tile

      ! If demanded: handle the dynamic PFT tiles in maps scheme
      IF (model%Is_process_enabled(NLCC_) .AND. TRIM(dsl4jsb_Config(ALCC_)%scheme) == 'maps') THEN
        CALL handle_dynamic_pfts(tile, options, initial_area, gained_area, lost_area)
      END IF

      ! Assertion: the total gained area needs to equal the total lost area to enable matter conservation
      cf_diff(:) = SUM(gained_area(:,:),DIM=2) - SUM(lost_area(:,:),DIM=2)
      ! Check for violations outside of the tolerated range (which are probably rather due to errors than precision)
      ! todo-AKN: need to reapply the check once mismatches have been treated correctly
      IF (ANY(ABS(cf_diff(:)) > max_tolerated_cf_mismatch)) THEN
        WRITE (message_text,*) 'Violation of assertion: gained area does not equal lost area! Please check input data. ' &
          & // 'Max difference of: ', MAXVAL(cf_diff(:))
        CALL finish(TRIM(routine), TRIM(message_text))
      END IF
      IF (model%config%check_child_cf) THEN
        ! Check and correct remaining mismatches
        CALL check_and_correct_cf_mismatches(tile, options, nr_of_tiles, gained_area, lost_area)
      END IF

      ! Now do the matter translocation
      CALL alcc_matter_translocation_jsb(tile, options)
    END DO ! idx for transition names

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_alcc_translocation


  ! ====================================================================================================== !
  !>
  !> #### Aggregation of variables from calculation of ALCC translocation
  !>
  SUBROUTINE aggregate_alcc_translocation(tile, options)

    dsl4jsb_Use_memory(FAGE_)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER  :: &
      & iblk, &       !< Current block index
      & ics, &        !< Index of first cell of block
      & ice           !< Index of last cell of block

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract  !< Aggregation method: area weighted fractions
    dsl4jsb_Def_memory(FAGE_)                              !< Memory of the forest age classes process
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_alcc_translocation'
    ! -------------------------------------------------------------------------------------------------- !
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    IF(tile%Is_process_active(FAGE_)) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
      dsl4jsb_Get_memory(FAGE_)
      dsl4jsb_Aggregate_onChunk(FAGE_, mean_age, weighted_by_fract)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_alcc_translocation


  ! ====================================================================================================== !
  !
  !> Do the matter translocation for jsbach.
  !
  SUBROUTINE alcc_matter_translocation_jsb(tile, options)

    USE mo_jsb_time,          ONLY: is_newday, is_newyear
    USE mo_jsb_lcc_class,     ONLY: t_jsb_lcc_proc
    USE mo_jsb_lcc,           ONLY: start_lcc_reloc, end_lcc_reloc, transfer_active_to_passive_onChunk

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(ALCC_) !< Memory of the anthropogenic land use change process
    dsl4jsb_Def_config(ALCC_) !< Configuration of the anthropogenic land use change process

    CLASS(t_jsb_tile_abstract), POINTER :: current_tile    !< Pointer to the tile for which relocations are applied
    CLASS(t_jsb_tile_abstract), POINTER :: forest_pft_tile !< Pointer for forest pft tiles
    TYPE(t_jsb_model), POINTER          :: model           !< Current instance of the model
    TYPE(t_jsb_lcc_proc), POINTER       :: lcc_relocations !< lcc relocations instance of the ALCC process

    REAL(wp), POINTER ::   &
      & initial_area(:,:), & !< Pointer to initial area in lcc relocations instance of the ALCC process
      & lost_area(:,:),    & !< Pointer to lost area in lcc relocations instance of the ALCC process
      & gained_area(:,:)     !< Pointer to gained area in lcc relocations instance of the ALCC process

    INTEGER  :: &
      & iblk, &       !< Current block index
      & ics, &        !< Index of first cell of block
      & ice, &        !< Index of last cell of block
      & nc, &         !< Number of cells in current block
      & ic, &         !< Cell index
      & i_tile        !< Index of a tile in the lcc_relocations
    REAL(wp) :: dtime !< Time step length

    LOGICAL  :: is_age_class !< Boolean indicating if a tile is a forest age class

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_alcc_translocation'
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(ALCC_)
    dsl4jsb_Get_memory(ALCC_)

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(alcc_)) RETURN

    !>
    !> Assert that PPLCC process is active
    !>
    IF (.NOT. model%Is_process_enabled(PPLCC_)) THEN
      CALL finish(TRIM(routine), 'Violation of precondition: lcc processes need pplcc to be active')
    END IF

    ! Redistributions are either only conducted at the start of each year or each day
    IF(.NOT. dsl4jsb_Config(ALCC_)%l_daily_alcc) THEN
      IF( .NOT. is_newyear(options%current_datetime,dtime)) RETURN
    ELSE
      IF( .NOT. is_newday(options%current_datetime,dtime)) RETURN
    END IF

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Assert, that this tile is the VEG tile
    IF (.NOT. tile%name .EQ. 'veg') THEN
      CALL finish(TRIM(routine), 'Violation of precondition: alcc processes is expected to run on the veg tile, instead' &
        & //' tried to run on '// trim(tile%name))
    END IF

    !>
    !> 1. Get lcc structure
    !>
    dsl4jsb_Get_lcc_relocations(ALCC_, lcc_relocations)

    ! get area vectors
    dsl4jsb_Get_lcc_area_matrix(initial_area)
    dsl4jsb_Get_lcc_area_matrix(lost_area)
    dsl4jsb_Get_lcc_area_matrix(gained_area)

    !>
    !> 2. Collect to be transferred matter
    !>
    CALL start_lcc_reloc(lcc_relocations, options, lost_area, initial_area)

    !>
    !> 3. Transfer matter from active to passive vars
    !>
    current_tile => tile%Get_first_child_tile()
    is_age_class = .FALSE.
    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile))

      ! In case of forest age classes we have to decent one level further
      IF ((TRIM(model%config%usecase) == 'jsbach_forest_age_classes') .AND. (current_tile%Has_children())) THEN
        is_age_class = .TRUE.
        forest_pft_tile => current_tile
        current_tile => current_tile%Get_first_child_tile()
      END IF

      IF (ANY(current_tile%name .EQ. lcc_relocations%tile_names)) THEN
        i_tile = i_tile + 1
        IF (ANY(lost_area(:, i_tile) > 0.0_wp)) THEN
          CALL transfer_active_to_passive_onChunk(lcc_relocations, current_tile, i_tile, options)
        END IF
      END IF

      current_tile => current_tile%Get_next_sibling_tile()

      ! In case of forest age classes we have to go up one level and go to the next tile there
      IF (is_age_class .AND. .NOT. ASSOCIATED(current_tile)) THEN
        current_tile => forest_pft_tile
        current_tile => current_tile%Get_next_sibling_tile()
        is_age_class = .FALSE.
      END IF
    END DO

    !>
    !> 4. make the passive matter transfer
    !>
    CALL end_lcc_reloc(lcc_relocations, options, gained_area)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE alcc_matter_translocation_jsb


  ! ====================================================================================================== !
  !
  !> Redistribution of losses and gains of a forest pft to its age-classes if running with age classes.
  !
  SUBROUTINE redistribute_cover_fraction_changes_of_forest_pft(tile, options, cf_diff, i_first_child, gained_area, lost_area)
    USE mo_fage_interface,         ONLY: apply_fract_per_age_change_for_alcc_forest_gain, &
      &                                  distribute_forest_loss_from_alcc
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile          !< Forest pft for which the routine is called
    TYPE(t_jsb_task_options),   INTENT(in)    :: options       !< Additional run-time parameters
    REAL(wp),                   INTENT(in)    :: cf_diff(:)    !< To be distributed changes of cover fractions
    INTEGER,                    INTENT(in)    :: i_first_child
       !< Index of first child tile in gained and lost area matrices
    REAL(wp),                   INTENT(inout) :: gained_area(:,:) !< Gained area matrix for lcc calculations
    REAL(wp),                   INTENT(inout) :: lost_area(:,:)   !< Lost area matrix for lcc calculations
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':redistribute_cover_fraction_changes_of_forest_pft'

    REAL(wp), DIMENSION(options%nc) :: this_loss !< Area loss of forest PFT to be distributed to its children (age classes)
    REAL(wp), DIMENSION(options%nc) :: this_gain !< Area gains of forest PFT to be distributed to its children (age classes)
    ! -------------------------------------------------------------------------------------------------- !
    this_loss(:) = -1_wp * MIN(cf_diff(:), 0._wp)
    this_gain(:) = MAX(cf_diff(:), 0._wp)

    ! per definition, all gains are distributed to the first age class
    gained_area(:, i_first_child) = gained_area(:, i_first_child) + this_gain(:)
    CALL apply_fract_per_age_change_for_alcc_forest_gain(tile, options, this_gain(:))

    ! how losses are distributed depends on configurations of the FAGE process
    CALL distribute_forest_loss_from_alcc(tile, options, this_loss(:), i_first_child, lost_area(:, :))

  END SUBROUTINE redistribute_cover_fraction_changes_of_forest_pft


  ! ====================================================================================================== !
  !> #### Handle antropogenic LCC if NLCC is active
  !
  !>  In case natural landcover change (nlcc) is active, the natural PFT fractions should not be
  !>  set to the map's fractions. The fractions of the dynamic PFT tiles need to be rescaled all
  !>  together, corresponding to the changes of the non-dynamic PFTs.
  !
  SUBROUTINE handle_dynamic_pfts(tile, options, initial_area, gained_area, lost_area)
    USE mo_nlcc_process,      ONLY: fract_small
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile              !< Parent tile for which the routine is called
    TYPE(t_jsb_task_options),   INTENT(in)    :: options           !< Additional run-time parameters
    REAL(wp),                   INTENT(inout) :: initial_area(:,:) !< Initial_area matrix for lcc calculations
    REAL(wp),                   INTENT(inout) :: gained_area(:,:)  !< Gained_area matrix for lcc calculations
    REAL(wp),                   INTENT(inout) :: lost_area(:,:)    !< Lost_area matrix for lcc calculations
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER  :: ic, i_tile                                !< Looping indices
    TYPE(t_jsb_model), POINTER          :: model          !< Current instance of the model
    CLASS(t_jsb_tile_abstract), POINTER :: current_tile   !< Pointer on currently investigated tile7
    REAL(wp), DIMENSION(options%nc)     :: current_fract  !< Cover fraction if current tile
    REAL(wp), DIMENSION(options%nc)     :: area_dynamic   !< Cover fraction correction
    REAL(wp), DIMENSION(options%nc)     :: gained_dynamic !< Cover fraction correction
    REAL(wp), DIMENSION(options%nc)     :: lost_dynamic   !< Cover fraction correction
    CHARACTER(len=*), PARAMETER :: routine = modname//':handle_dynamic_pfts'
    ! -------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)

    ! Init area vectors
    lost_dynamic(:)   = 0.0_wp
    gained_dynamic(:) = 0.0_wp
    area_dynamic(:)   = 0.0_wp

    current_tile => tile%Get_first_child_tile()
    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile))
      i_tile = i_tile + 1

      IF (dsl4jsb_Lctlib_param_tile(dynamic_PFT, current_tile)) THEN
        ! Sum up the total area fraction of the dynamic PFT tiles
        area_dynamic(:) = area_dynamic(:) + initial_area(:,i_tile)
      ELSE
        ! The area losses/gains of the non-dynamic PFTs define the joint gains/losses
        ! of the dynamic PFTs
        gained_dynamic(:) = gained_dynamic(:) + lost_area(:,i_tile)
        lost_dynamic(:)   = lost_dynamic(:)   + gained_area(:,i_tile)
      END IF
      current_tile => current_tile%Get_next_sibling_tile()
    END DO

    current_tile => tile%Get_first_child_tile()
    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile))
      i_tile = i_tile + 1

      IF (dsl4jsb_Lctlib_param_tile(dynamic_PFT, current_tile)) THEN
        DO ic = 1, options%nc
          ! The dynamic PFTs are gaining or losing area proportional to their current area fractions
          IF (area_dynamic(ic) > fract_small) THEN
            gained_area(ic,i_tile) = MIN(gained_dynamic(ic) * initial_area(ic,i_tile) / area_dynamic(ic), &
              &                          1._wp - initial_area(ic,i_tile))
            lost_area(ic,i_tile)   = MIN(lost_dynamic(ic)   * initial_area(ic,i_tile) / area_dynamic(ic), &
              &                          initial_area(ic,i_tile))
          ELSE
            gained_area(ic,i_tile) = 0._wp
            lost_area(ic,i_tile)   = 0._wp
          END IF

          current_fract(ic) = initial_area(ic, i_tile) + gained_area(ic,i_tile) - lost_area(ic,i_tile)
        END DO
        ! Update the fractions of dynamic PFTs
        CALL current_tile%Set_fraction(options%ics, options%ice, options%iblk, fract=current_fract(:))

      END IF
      current_tile => current_tile%Get_next_sibling_tile()
    END DO

  END SUBROUTINE handle_dynamic_pfts


  ! ====================================================================================================== !
  !> #### Check and if required correct gained and lost areas
  !>
  !> Verify that the total gained area fraction equals the lost area fractions and apply corrections in
  !> case of a mismatch.
  !
  SUBROUTINE check_and_correct_cf_mismatches(tile, options, nr_of_tiles, gained_area, lost_area)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile             !< Parent tile for which the routine is called
    TYPE(t_jsb_task_options),   INTENT(in)    :: options          !< Additional run-time parameters
    INTEGER,                    INTENT(in)    :: nr_of_tiles      !< Number of tiles involved in the lcc calculation
    REAL(wp),                   INTENT(inout) :: gained_area(:,:) !< Gained_area matrix for lcc calculations
    REAL(wp),                   INTENT(inout) :: lost_area(:,:)   !< Lost_area matrix for lcc calculations
    ! -------------------------------------------------------------------------------------------------- !
    LOGICAL  :: gained_exceeds_loss !< Logical indicating if an area gain exceeds an area loss; used to assert carbon conservation
    INTEGER  :: ic, i_tile, i_sort  !< Looping indices
    REAL(wp) :: remaining_diff      !< Remaining difference between area losses and gains
    INTEGER  :: i_max_changing_tile !< Tile with the largest area change left
    CLASS(t_jsb_tile_abstract), POINTER            :: current_tile  !< Pointer on currently investigated tile
    REAL(wp), DIMENSION(options%nc)                :: current_fract !< Cover fraction if current tile
    REAL(wp), DIMENSION(options%nc)                :: cf_diff       !< Difference between gained and lost areas
    LOGICAL,  DIMENSION(nr_of_tiles)               :: masked_tiles  !< Mask of already considered tiles
    REAL(wp), DIMENSION(nr_of_tiles)               :: area_change   !< Area change to sort for (either losses or gains)
    REAL(wp), DIMENSION(options%nc, nr_of_tiles) :: cf_correction   !< Cover fraction correction
    CHARACTER(len=*), PARAMETER :: routine = modname//':check_and_correct_cf_mismatches'
    ! -------------------------------------------------------------------------------------------------- !
    cf_diff(:) = SUM(gained_area(:,:),DIM=2) - SUM(lost_area(:,:),DIM=2)
    cf_correction(:,:) = 0.0_wp
    area_change(:) = 0.0_wp

    IF (ANY(ABS(cf_diff(:)) > EPSILON(1._wp))) THEN  !EPSILON(1._wp) ! 1.E-15_wp
      DO ic = 1, options%nc
        IF (ABS(cf_diff(ic)) > EPSILON(1._wp)) THEN
          masked_tiles(:) = .FALSE.
          remaining_diff = cf_diff(ic)
          DO i_sort = 1, nr_of_tiles
            IF (remaining_diff > EPSILON(1._wp)) THEN
              ! I.e. gain > loss
              area_change(:) = gained_area(ic,:)
              i_tile = MAXLOC(area_change(:), DIM=1, MASK=masked_tiles(:))
              cf_correction(ic,i_tile) = MIN(remaining_diff, gained_area(ic,i_tile))
              remaining_diff = MAX(0.0_wp, remaining_diff - gained_area(ic,i_tile))
              gained_area(ic,i_tile) = gained_area(ic,i_tile) - cf_correction(ic,i_tile)
              masked_tiles(i_tile) = .TRUE.
            ELSE IF (remaining_diff < -EPSILON(1._wp)) THEN
              ! I.e. loss > gain
              area_change(:) = lost_area(ic,:)
              i_tile = MAXLOC(area_change(:), DIM=1, MASK=masked_tiles(:))
              cf_correction(ic,i_tile) = - MIN(ABS(remaining_diff), lost_area(ic,i_tile))
              remaining_diff = - MAX(0.0_wp, ABS(remaining_diff) - lost_area(ic,i_tile))
              lost_area(ic,i_tile) = lost_area(ic,i_tile) - cf_correction(ic,i_tile)
              masked_tiles(i_tile) = .TRUE.
            ELSE
              EXIT
            END IF
          END DO
        END IF
      END DO ! ic = 1, nc

      ! If requried modify tiles and area changes accordingly
      current_tile => tile%Get_first_child_tile()
      i_tile = 0
      DO WHILE (ASSOCIATED(current_tile))
        i_tile = i_tile + 1
        CALL current_tile%Get_fraction(options%ics, options%ice, options%iblk, fract=current_fract(:))
        current_fract(:) = current_fract(:) - cf_correction(:,i_tile)
        CALL current_tile%Set_fraction(options%ics, options%ice, options%iblk, fract=current_fract(:))
        current_tile => current_tile%Get_next_sibling_tile()
      END DO ! WHILE (ASSOCIATED(current_tile))
    END IF ! (ANY(ABS(cf_diff) > EPSILON(1._wp))) THEN

  END SUBROUTINE check_and_correct_cf_mismatches

#endif
END MODULE mo_alcc_interface
