! Meri's to-be-deleted notes:
! Copied flcc process: flcc->hlcc, burned -> harvested
! to do: make a dummy "harvested fraction" in the beginning to start testing
! to do: now this is only daily, but could also have annual like in alcc (but could start with just daily cos not a priority)
! Need to add some check (not here though) for harvesting age classes and model age classes so that they match
! to do: Now there is no "separate" option for harvesting (see e.g. wlcc config) but it is attempted to read here
! (Why alcc is 3D...? ("Get_var3D_onChunk"))
! does proper harvesting need pplc? 
!
!---------------------------------------
!> interface to the hlcc process (lcc due to (wood) harvesting))
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
MODULE mo_hlcc_interface
    #ifndef __NO_JSBACH__
    
      USE mo_jsb_control,     ONLY: debug_on
      USE mo_kind,            ONLY: wp
      USE mo_exception,       ONLY: message, finish
    
      USE mo_jsb_model_class,    ONLY: t_jsb_model
      USE mo_jsb_class,          ONLY: Get_model
      USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
      USE mo_jsb_process_class,  ONLY: t_jsb_process
      USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
    
      ! Use of processes in this module
      dsl4jsb_Use_processes HLCC_, CARBON_, FAGE_ ! MR-TODO : PPLCC_? (in alcc <- also without C) - CARBON I think is needed, PLCC probs later when input data
    
      ! Use of process configurations
      !dsl4jsb_Use_config(HLCC_) !MR-TODO : Not sure if needed here, will be done later too
      !dsl4jsb_Use_config(CARBON_)

      ! Use of process memories !MR-TODO : These come later again, around line 110 - unnecessary for the whole thing and enought there in the update routine?
      !dsl4jsb_Use_memory(HLCC_)
      !dsl4jsb_Use_memory(CARBON_)

      IMPLICIT NONE
      PRIVATE
      PUBLIC :: Register_hlcc_tasks
    
      CHARACTER(len=*), PARAMETER :: modname = 'mo_hlcc_interface'
    
      !> Type definition for hlcc pre task
      TYPE, EXTENDS(t_jsb_process_task) :: tsk_hlcc
      CONTAINS
        PROCEDURE, NOPASS :: Integrate => update_hlcc     !< Advances task computation for one timestep
        PROCEDURE, NOPASS :: Aggregate => aggregate_hlcc  !< Aggregates computed task variables
      END TYPE tsk_hlcc
    
      !> Constructor interface for hlcc pre task
      INTERFACE tsk_hlcc
        PROCEDURE Create_task_hlcc                        !< Constructor function for task
      END INTERFACE tsk_hlcc
    
    CONTAINS
    
      ! ====================================================================================================== !
      !
      !> Constructor for hlcc task
      !
      FUNCTION Create_task_hlcc(model_id) RESULT(return_ptr)
    
        ! -------------------------------------------------------------------------------------------------- !
        INTEGER,                   INTENT(in) :: model_id    !< Model id
        CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "hlcc"
        ! -------------------------------------------------------------------------------------------------- !
        ALLOCATE(tsk_hlcc::return_ptr)
        CALL return_ptr%Construct(name='hlcc', process_id=HLCC_, owner_model_id=model_id) ! MR-TODO hlcc and HLCC_ to relevant places
    
      END FUNCTION Create_task_hlcc
    
      ! ====================================================================================================== !
      !
      !> Register tasks for hlcc process
      !
      SUBROUTINE Register_hlcc_tasks(this, model_id)
    
        ! -------------------------------------------------------------------------------------------------- !
        CLASS(t_jsb_process), INTENT(inout) :: this        !< Instance of hlcc process class
        INTEGER,               INTENT(in)   :: model_id    !< Model id
        ! -------------------------------------------------------------------------------------------------- !
        CALL this%Register_task(tsk_hlcc(model_id))
    
      END SUBROUTINE Register_hlcc_tasks
    
      ! ====================================================================================================== !
      !
      !> Implementation of "update" for task "hlcc"
      !
      ! Harvesting of a fraction implemented here ! MR-TODO
      !>
      !> Task "hlcc" shifts harvested area (from youngest AC or "from tile back to the same tile") and cares for quantities conserved upon hlcc
      !
      SUBROUTINE update_hlcc(tile, options)
        USE mo_hlcc_process,       ONLY: harvested_fract_jsbach
        !USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract ! MR-TODO : not repeated here in flcc either, since already done
        USE mo_jsb_time,           ONLY: is_newday
        USE mo_jsb_lcc_class,      ONLY: min_daily_cf_change, t_jsb_lcc_proc, &
          &                              transfer_from_active_to_passive_var_onChunk, copy_from_active_to_passive_var_onChunk
        USE mo_jsb_lcc,            ONLY: init_lcc_reloc, start_lcc_reloc, end_lcc_reloc, transfer_active_to_passive_onChunk
        USE mo_carbon_constants,   ONLY: fract_green_aboveGround, fract_wood_aboveGround, sec_per_day, molarMassCO2_kg
        USE mo_fage_interface,     ONLY: recalc_fract_per_age_upon_disturbance
    
        dsl4jsb_Use_memory(HLCC_)   !MR-TODO : could also do in beginning
        dsl4jsb_Use_memory(CARBON_)
        dsl4jsb_Use_config(HLCC_)
        dsl4jsb_Use_config(CARBON_)
        
    
        ! -------------------------------------------------------------------------------------------------- !
        CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
        TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
        ! -------------------------------------------------------------------------------------------------- !
        CHARACTER(len=*), PARAMETER :: routine = modname//':update_hlcc'
        LOGICAL, PARAMETER :: allAtOnce = .FALSE. ! MR-TODO : what's this? can't see where it would be used, but it exists in also fage and some other interface
    
        CLASS(t_jsb_tile_abstract), POINTER :: current_tile, d_tile
        TYPE(t_jsb_model), POINTER          :: model
        TYPE(t_jsb_lcc_proc), POINTER       :: lcc_relocations
    
        INTEGER  :: i, j, i_tile, iblk, ics, ice, nc
        INTEGER  :: nr_of_tiles, nr_of_lcc_events, nr_of_combined_tiles
        REAL(wp) :: dtime
        REAL(wp), ALLOCATABLE :: lost_area(:,:), gained_area(:,:), initial_area(:,:)
    
        REAL(wp), DIMENSION(options%nc) :: current_fract
        REAL(wp) :: fire_fract_wood_2_atmos, conversion_fact, transfer_fract ! MR-TODO : need to replace fice_frac... with something else at this stage (and possibly other fracs too)
    
        LOGICAL  :: d_tile_is_ac       ! if the disturbed tile is a ac tile (only if running with age classes)
        INTEGER  :: i_first_ac         ! index of first ac (only if running with age classes)
        CLASS(t_jsb_tile_abstract), POINTER :: first_ac ! pointer to first ac (only if running with age classes)
    
        ! ---------
            ! Declare process memories
             dsl4jsb_Def_memory(HLCC_)
             dsl4jsb_Def_memory(CARBON_)
             dsl4jsb_Def_config(HLCC_)
             dsl4jsb_Def_config(CARBON_)

             
             dsl4jsb_Real2D_onChunk :: dtime !MR-TODO : same name as internal variable? problem?
             dsl4jsb_Real2D_onChunk :: current_datetime
             dsl4jsb_Real2D_onChunk :: harvested_fract

             dsl4jsb_Real2D_onChunk :: co2flux_fire_all_2_atm !MR-TODO : carbon process?
             dsl4jsb_Real2D_onChunk :: cflux_dist_green_2_soil
             dsl4jsb_Real2D_onChunk :: cflux_dist_woods_2_soil

             ! dsl4jsb_Real2D_onChunk :: cover_fract_pot ! MR-TODO : what is this?
        ! ---------

             dsl4jsb_Def_memory_tile(HLCC_, d_tile) ! MR-TODO : these retrieved later with Get_memory_tile... BUT I don't know why/if this is needed
             dsl4jsb_Real2D_onChunk :: harvested_fract_d_tile !MR-TODO : just declaring variable (why?)
        ! -------------------------------------------------------------------------------------------------- !
    
        ! Get local variables from options argument
        iblk    = options%iblk
        ics     = options%ics
        ice     = options%ice
        nc      = options%nc
        dtime   = options%dtime
    
        model => Get_model(tile%owner_model_id)
    
        ! Get process configurations ! MR-TODO : already done?
        dsl4jsb_Get_config(HLCC_)
        dsl4jsb_Get_config(CARBON_)
    
        ! If process is not to be calculated on this tile, do nothing
        IF (.NOT. tile%Is_process_calculated(HLCC_)) RETURN
    
        ! If not newday, do nothing
        IF( .NOT. is_newday(options%current_datetime,dtime)) RETURN
    
        IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    
        ! Assertion: hlcc is currently implemented to run on   !MR-TODO : Not sure about this - anyway supposed to be done on woody pft's anyway (later) so does this make a difference?
        ! - the veg tile in the pft use case
        ! - the pfts or pftXX_fcXX tiles (if running with age classes)
        IF(.NOT. INDEX(tile%name, 'pft') > 0) THEN
          IF(.NOT. (TRIM(model%config%usecase) == 'jsbach_pfts' .AND. (INDEX(tile%name, 'veg') > 0))) THEN
            CALL finish(TRIM(routine), 'Violation of assertion: hlcc (currently) is implemented to run on ' &
              & //' the pft tiles or the veg tile in the pft use case or - on the pft tiles or the '        &
              & //' acs (if with age classes), but tried to run for: ' // trim(tile%name) // '. Please check.')
          END IF
        END IF

        !>
        !> Get harvested fraction 
        !>

        ! Get process memories
        dsl4jsb_Get_memory(HLCC_)
        dsl4jsb_Get_memory(CARBON_)

        ! Get process variables - MR-TODO : not sure if these are correct (Get?). Also same variables as elsewhere

        dsl4jsb_Get_var2D_onChunk(HLCC_,   dtime)         ! in
        dsl4jsb_Get_var2D_onChunk(HLCC_,   current_datetime)    ! in

        dsl4jsb_Get_var2D_onChunk(HLCC_,   harvested_fract)             ! out

        dsl4jsb_Get_var2D_onChunk(CARBON_,    co2flux_fire_all_2_atm)   ! out
        dsl4jsb_Get_var2D_onChunk(CARBON_,    cflux_dist_green_2_soil)  ! out
        dsl4jsb_Get_var2D_onChunk(CARBON_,    cflux_dist_woods_2_soil)  ! out

        IF (is_newday(options%current_datetime, options%dtime)) THEN ! only once per day - MR-TODO: but there is already a check fot this earlier??
          harvested_fract(:) = 0._wp
          IF (dsl4jsb_Lctlib_param(dynamic_PFT)) THEN ! MR-TODO: what's this
    
    
            ! Calculate harvested area. ("Dummy harvesting")
            IF (dsl4jsb_Lctlib_param(woody_PFT)) THEN  ! MR-TODO: Woody Type only, harvested should otherwise be 0 (defined earlier)

    
              CALL harvested_fract_jsbach(                            &
                & options%dtime,                                      & ! in
                & options%current_datetime,                           & ! in
                & harvested_fract(:)                                  & ! inout
                & )


            ENDIF

          END IF
    
    
        END IF ! is_newday  



        !--------------------------------------------------------
        !>
        !> 1. Get lcc structure
        !>
        dsl4jsb_Get_lcc_relocations(HLCC_, lcc_relocations)
        nr_of_tiles = lcc_relocations%nr_of_tiles
    
        ! Allocate area vectors
        ALLOCATE(initial_area(nc, nr_of_tiles))
        ALLOCATE(lost_area(nc, nr_of_tiles))
        ALLOCATE(gained_area(nc, nr_of_tiles))
        initial_area(:,:) = 0.0_wp
    
        !>
        !> 2. Two possibilities:
        !>
        !>     2.a make all area movements at once (in that case there is no directed movement of matter
        !>         from one source to one sink but only a movement of matter from all shrinking to all extending tiles)
        !>
        !>     2.b make an outer loop around init, start, end, dealloc and only move area from one source to one sink tile
        !>         [TR: suggested to permute the order of tiles in the loop for alcc with transitions]
        !>
        IF (dsl4jsb_Config(HLCC_)%separate) THEN ! => 2b
          nr_of_lcc_events = nr_of_tiles
          nr_of_combined_tiles = 1
        ELSE ! not in a do loop but all at once    => 2a
          nr_of_lcc_events = 1
          nr_of_combined_tiles = nr_of_tiles
        ENDIF
    
        ! initiate pointers, logicals and indices to enable lcc, particularly if running with forest age classes
        d_tile_is_ac = .FALSE.
        IF (tile%Has_children()) THEN
          d_tile => tile%Get_first_child_tile()
          IF (d_tile%Has_children()) THEN
            CALL finish(TRIM(routine), 'Violation of assertion: hlcc (currently) is not implemented to run on ' &
              & //'a hierarchy level with more descendants than direct children. Tried to run on ' &
              & // trim(tile%name) // '. Please check.')
          END IF
          IF (TRIM(model%config%usecase) == 'jsbach_forest_age_classes') THEN
            ! Assert: for hlcc (currently) first age class is expected first in the list of involved tiles
            IF (.NOT. (d_tile%name .EQ. lcc_relocations%tile_names(1))) THEN
              CALL finish(TRIM(routine), 'Violation of assertion: for hlcc (currently) the first age class is ' &
                & //'expected to be listed as first involved tile, instead found: '// trim(d_tile%name))
            END IF
            first_ac => d_tile
            i_first_ac = 1
            d_tile_is_ac = .TRUE.
          END IF
        ELSE
          d_tile => tile%self_tile
        END IF
    
        DO i = 1,nr_of_lcc_events ! < outer loop enclosing one full lcc-cycle
          ! Collect initial area of all involved tiles
          CALL init_lcc_reloc(lcc_relocations, options, tile, initial_area)
    
          ! init gained and lost area
          lost_area(:,:) = 0.0_wp
          gained_area(:,:) = 0.0_wp
    
          ! Collect lost and gained area of involved tiles
          DO j =1,nr_of_combined_tiles ! < inner loop collecting lcc movements within this lcc-cycle
    
            IF (dsl4jsb_Config(HLCC_)%separate) THEN
              ! in this case the area redistribution follows the counter of the outer loop, i.e. one pft a time
              i_tile = i
            ELSE
              ! in this case the area redistribution follows the counter of the inner loop,
              ! i.e. all area is aggregated in this one loop (i.e. mixed) and then distributed at once
              i_tile = j
            ENDIF
    
            IF (TRIM(model%config%usecase) == 'jsbach_forest_age_classes') THEN
              ! Assertion: if running with age classes all classes are expected to be involved with hlcc
              IF (.NOT. (d_tile%name .EQ. lcc_relocations%tile_names(i_tile))) THEN
                CALL finish(TRIM(routine), 'Violation of assertion: for hlcc (currently) all age classes are expected ' &
                  & //'to be listed as involved, but found not involved age-class: '// trim(d_tile%name))
              END IF
            END IF
    
            ! - add harvested area of current child to lost area --------------------------------------------------------------------MR-TODO!!!
            CALL d_tile%Get_fraction(ics, ice, iblk, fract=current_fract(:))
            dsl4jsb_Get_memory_tile(HLCC_, d_tile) !MR-TODO : this d_tile structure is from flcc (with DISTURB_) not sure if needs to be done this way...
            dsl4jsb_Get_var2d_onChunk_tile_name(HLCC_, harvested_fract, d_tile) !MR-TODO : Not sure if works (basically allocating the just calculated harvested_frac to this d_tile)... (burned_frac is not "called" in flcc tho)
            !----MR-TODO
            ! where does harvested_fract come from? Apparently "defined" in process, and then comes through here through memory...there: harvested_fract is "TYPE(t_jsb_var_real2d)"      
            ! TODO : make the harvesting "process" in the process module, and make memory to match... then probs no other major things to do here? 

            !----MR-TODO
            lost_area(:,i_tile) = harvested_fract_d_tile(:) * current_fract(:)
            WHERE(lost_area(:,i_tile) < min_daily_cf_change)
              lost_area(:,i_tile) = 0.0_wp
            END WHERE
    
            IF (ANY(lost_area(:,i_tile) > 0._wp)) THEN
              ! - care for the gained area  -> @todo change when bare land tile available
              ! ... if forest age classes are used the area of any disturbed age class is added to the first class
              IF (d_tile_is_ac) THEN
                gained_area(:,i_first_ac) = gained_area(:,i_first_ac) + lost_area(:,i_tile)
    
                ! for all but the first age class thus the fractions need to be adapted
                IF (.NOT. d_tile%Is_first_child()) THEN
                  current_fract(:) = current_fract(:) -  lost_area(:,i_tile)
                  CALL d_tile%Set_fraction(ics, ice, iblk, fract=current_fract(:))
                  CALL first_ac%Get_fraction(ics, ice, iblk, fract=current_fract(:))
                  current_fract(:) = current_fract(:) + lost_area(:,i_tile)
                  CALL first_ac%Set_fraction(ics, ice, iblk, fract=current_fract(:))
    
                  ! Updated tracked age
                  CALL recalc_fract_per_age_upon_disturbance(tile, options, i_tile, lost_area(:,i_tile))
                END IF
    
              ELSE
                ! else (currently) back to the tile itself
                gained_area(:,i_tile) = harvested_fract_d_tile(:) * current_fract(:)
              END IF
            END IF
    
            ! - go to next tile (if associated)
            d_tile => d_tile%Get_next_sibling_tile()
    
          END DO ! inner loop collecting lcc movements within this lcc-cycle
    
          ! Only needs to be executed if there is area to be moved
          IF (ANY(lost_area(:, :) > 0.0_wp)) THEN
    
            !>
            !> 3. Each lcc process needs to call start lcc to collect to be transferred matter
            !>
            CALL start_lcc_reloc(lcc_relocations, options, lost_area, initial_area)
    
            !>
            !> 4. Afterwards active content needs to be relocated to passive vars by manipulating the arrays in lcc_relocations
            !>
            !>     4.a This can be done by explicitly specifying how active variables should be treated
            !>         e.g. by calling transfer_from_active_to_passive_var_onChunk for named variables
            !>
            !>     4.b and finally by calling the more general transfer_active_to_passive_onChunk lcc routine where calls to
            !>         process specific relocation routines can be specified
            !>
            IF (tile%Has_children()) THEN
              current_tile => tile%Get_first_child_tile()
            ELSE
              current_tile => tile%self_tile
            END IF
    
            i_tile = 0
            DO WHILE (i_tile < nr_of_tiles)
    
              i_tile = i_tile + 1
              IF (ANY(lost_area(:, i_tile) > 0.0_wp)) THEN
                ! 4a -- 'special transfer'
                IF (current_tile%Has_process_memory(CARBON_)) THEN
                  conversion_fact = molarMassCO2_kg / sec_per_day
                  ! all litter pools (ag yasso pools) are harvested
                  CALL transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
                    &  'c_acid_ag1', 'co2flux_fire_all_2_atm', conversion_fact = conversion_fact)
                  CALL transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
                    &  'c_water_ag1', 'co2flux_fire_all_2_atm', conversion_fact = conversion_fact)
                  CALL transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
                    &  'c_ethanol_ag1', 'co2flux_fire_all_2_atm', conversion_fact = conversion_fact)
                  CALL transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
                    &  'c_nonsoluble_ag1', 'co2flux_fire_all_2_atm', conversion_fact = conversion_fact)
                  CALL transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
                    &  'c_acid_ag2', 'co2flux_fire_all_2_atm', conversion_fact = conversion_fact)
                  CALL transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
                    &  'c_water_ag2', 'co2flux_fire_all_2_atm', conversion_fact = conversion_fact)
                  CALL transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
                    &  'c_ethanol_ag2', 'co2flux_fire_all_2_atm', conversion_fact = conversion_fact)
                  CALL transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
                    &  'c_nonsoluble_ag2', 'co2flux_fire_all_2_atm', conversion_fact = conversion_fact)
                  ! treat vegetation pools: these need to be treated explicitly here...
                  ! - first: copy matter for diagnostics
                  !   IMPORTANT: of course needs to happen before the ag fraction is transferred!
                  ! - second: actively transfer harvested parts
                  ! - finally: below ground parts will be actively transferred to soils by the generic routine
                  !            in the carbon interface
                  fire_fract_wood_2_atmos = dsl4jsb_Config(CARBON_)%fire_fract_wood_2_atmos
                  transfer_fract = 1._wp - (fire_fract_wood_2_atmos * fract_wood_aboveGround)
                  conversion_fact = 1.0_wp / sec_per_day
                  CALL copy_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
                    & 'c_woods', 'cflux_dist_woods_2_soil', copy_fract = transfer_fract,                &
                    & conversion_fact = conversion_fact)
                  transfer_fract = 1._wp - fract_green_aboveGround
                  CALL copy_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
                    & 'c_green', 'cflux_dist_green_2_soil', copy_fract = transfer_fract,                &
                    & conversion_fact = conversion_fact)
                  CALL copy_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
                    & 'c_reserve', 'cflux_dist_green_2_soil', copy_fract =  transfer_fract,             &
                    & conversion_fact = conversion_fact)
                  ! - parts of woody ag part and all green and reserve ag parts are harvested
                  transfer_fract = fire_fract_wood_2_atmos * fract_wood_aboveGround
                  conversion_fact = molarMassCO2_kg / sec_per_day
                  CALL transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
                    &  'c_woods', 'co2flux_fire_all_2_atm', transfer_fract, conversion_fact)
                  transfer_fract = fract_green_aboveGround
                  CALL transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
                    &  'c_green', 'co2flux_fire_all_2_atm', transfer_fract, conversion_fact)
                  CALL transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
                    &  'c_reserve', 'co2flux_fire_all_2_atm', transfer_fract, conversion_fact)
                END IF
    
                ! 4b -- general transfer of active vars (transfer the 'rest' of active vars)
                CALL transfer_active_to_passive_onChunk(lcc_relocations, current_tile, i_tile, options)
              END IF
    
              ! - go to next tile (if associated)
              current_tile => current_tile%Get_next_sibling_tile()
            END DO
    
            !>
            !> 5. Each lcc process needs to call end lcc to make the passive matter transfer
            !>
            CALL end_lcc_reloc(lcc_relocations, options, gained_area)
    
          END IF ! ... only to be conducted if there is any area movement
        END DO ! outer loop enclosing one full lcc-cycle: depending on if are movement is to be "separate" or mixed
    
        !>
        !> 6. Deallocate area vectors
        !>
        DEALLOCATE(lost_area, gained_area, initial_area)
    
        IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
    
      END SUBROUTINE update_hlcc
    
      ! ====================================================================================================== !
      !
      !> Implementation of "aggregate" for task "hlcc"
      !
      SUBROUTINE aggregate_hlcc(tile, options)
    
        dsl4jsb_Use_memory(FAGE_)
        ! -------------------------------------------------------------------------------------------------- !
        CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
        TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
        ! -------------------------------------------------------------------------------------------------- !
        INTEGER  :: iblk, ics, ice
        !X if necessary: REAL(wp) :: dtime, steplen
        CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract
        dsl4jsb_Def_memory(FAGE_)
        CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_hlcc'
        ! -------------------------------------------------------------------------------------------------- !
    
        ! Get local variables from options argument
        iblk = options%iblk
        ics  = options%ics
        ice  = options%ice
    
        !X if necessary: dtime   = options%dtime
        !X if necessary: steplen = options%steplen
    
        IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    
        IF(tile%Is_process_active(FAGE_)) THEN
          weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
          dsl4jsb_Get_memory(FAGE_)
          dsl4jsb_Aggregate_onChunk(FAGE_, mean_age, weighted_by_fract)
        ENDIF
    
        IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
    
      END SUBROUTINE aggregate_hlcc
    
    #endif
    END MODULE mo_hlcc_interface