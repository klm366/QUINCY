!> Summary: interface to the wlcc process (lcc due to wind)
!>
!>### Contains the interface to the wlcc process (lcc due to wind)
!>
!>
!> Copyright (c) 2022 Julia Nabel
!>
!>#### Revision History
!>
!> First version and ongoing changes   by Julia Nabel (2022-01-06)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
MODULE mo_wlcc_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_varlist,     ONLY: VARNAME_LEN

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model    
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes WLCC_, PPLCC_, DISTURB_, CARBON_, PHENO_

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_wlcc_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_wlcc_interface'
  CHARACTER(len=*), PARAMETER :: procname = 'wlcc'
  
  !> Type definition for wlcc pre task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_wlcc
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_wlcc     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_wlcc  !< Aggregates computed task variables
  END TYPE tsk_wlcc
  
  !> Constructor interface for wlcc pre task
  INTERFACE tsk_wlcc
    PROCEDURE Create_task_wlcc                        !< Constructor function for task
  END INTERFACE tsk_wlcc
  
CONTAINS

  ! ====================================================================================================== !
  !
  !> Constructor for wlcc task
  !
  FUNCTION Create_task_wlcc(model_id) RESULT(return_ptr)

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id    !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "wlcc"
    ! -------------------------------------------------------------------------------------------------- !
    ALLOCATE(tsk_wlcc::return_ptr)
    CALL return_ptr%Construct(name=procname, process_id=WLCC_, owner_model_id=model_id)

  END FUNCTION Create_task_wlcc
    
  ! ====================================================================================================== !
  !
  !> Register tasks for wlcc process
  !
  SUBROUTINE Register_wlcc_tasks(this, model_id)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_process), INTENT(inout) :: this        !< Instance of wlcc process class
    INTEGER,               INTENT(in)   :: model_id    !< Model id
    ! -------------------------------------------------------------------------------------------------- !
    CALL this%Register_task(tsk_wlcc(model_id))
    
  END SUBROUTINE Register_wlcc_tasks

  ! ====================================================================================================== !
  !
  !> Implementation of "update" for task "wlcc"
  !>
  !> Task "wlcc" shifts burned area and cares for quantities conserved upon wlcc
  !
  SUBROUTINE update_wlcc(tile, options)

    USE mo_jsb_time,          ONLY: is_newday
    USE mo_jsb_lcc_class,     ONLY: t_jsb_lcc_proc, &
      &                             transfer_from_active_to_passive_var_onChunk, copy_from_active_to_passive_var_onChunk
    USE mo_jsb_lcc,           ONLY: init_lcc_reloc, start_lcc_reloc, end_lcc_reloc, transfer_active_to_passive_onChunk
    USE mo_carbon_constants,  ONLY: fract_green_aboveGround, fract_wood_aboveGround, sec_per_day

    dsl4jsb_Use_memory(PHENO_)
    dsl4jsb_Use_memory(DISTURB_)
    dsl4jsb_Use_config(WLCC_)
    dsl4jsb_Use_config(CARBON_)

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_wlcc'
    LOGICAL, PARAMETER :: allAtOnce = .FALSE.

    CLASS(t_jsb_tile_abstract), POINTER :: current_tile, disturbed_tile
    TYPE(t_jsb_model), POINTER          :: model
    TYPE(t_jsb_lcc_proc), POINTER       :: lcc_relocations

    INTEGER  :: i, j, i_tile, iblk, ics, ice, nc
    INTEGER  :: nr_of_tiles, nr_of_lcc_events, nr_of_combined_tiles
    REAL(wp) :: dtime, conversion_fact
    REAL(wp), ALLOCATABLE :: lost_area(:,:), gained_area(:,:), initial_area(:,:)

    REAL(wp), DIMENSION(options%nc) :: current_fract

    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_config(WLCC_)
    dsl4jsb_Def_memory_tile(DISTURB_, disturbed_tile)
    dsl4jsb_Def_memory_tile(PHENO_, disturbed_tile)
    dsl4jsb_Real2D_onChunk :: damaged_fract_disturbed_tile
    dsl4jsb_Real2D_onChunk :: veg_fract_correction_disturbed_tile
    dsl4jsb_Real2D_onChunk :: fract_fpc_max_disturbed_tile
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    
    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(WLCC_)
    dsl4jsb_Get_config(CARBON_)

    !>
    !> 1. Get lcc structure
    !>
    dsl4jsb_Get_lcc_relocations(WLCC_, lcc_relocations)
    nr_of_tiles = lcc_relocations%nr_of_tiles
    
    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_active(WLCC_)) RETURN

    ! If not newday, do nothing
    IF( .NOT. is_newday(options%current_datetime,dtime)) RETURN
    
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Assert, that this tile is the VEG tile 
    ! - this lcc process works on the PFT tiles but is controlled on their parent tile
    IF (.NOT. tile%name .EQ. 'veg') THEN
      CALL finish(TRIM(routine), 'Violation of precondition: wlcc processes is expected to run on the veg tile, instead' &
        & //' tried to run on '// trim(tile%name))
    ENDIF

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
    IF (dsl4jsb_Config(WLCC_)%separate) THEN ! => 2b
      nr_of_lcc_events = tile%Get_no_of_children()
      nr_of_combined_tiles = 1
    ELSE ! not in a do loop but all at once    => 2a
      nr_of_lcc_events = 1
      nr_of_combined_tiles = tile%Get_no_of_children()
    ENDIF

    ! In both cases: start with the first child
    disturbed_tile => tile%Get_first_child_tile()

    DO i =1,nr_of_lcc_events ! < outer loop enclosing one full lcc-cycle

      ! Collect initial area of all involved tiles
      CALL init_lcc_reloc(lcc_relocations, options, tile, initial_area)

      ! init gained and lost area
      lost_area(:,:) = 0.0_wp
      gained_area(:,:) = 0.0_wp

      ! Collect lost and gained area of involved tiles
      DO j =1,nr_of_combined_tiles ! < inner loop collecting lcc movements within this lcc-cycle

        IF (dsl4jsb_Config(WLCC_)%separate) THEN
          ! in this case the area redistribution follows the counter of the outer loop, i.e. one pft a time
          i_tile = i
        ELSE
          ! in this case the area redistribution follows the counter of the inner loop,
          ! i.e. all area is aggregated in this one loop (i.e. mixed) and then distributed at once
          i_tile = j
        ENDIF

        ! - add damaged area of current child to lost area
        CALL disturbed_tile%Get_fraction(ics, ice, iblk, fract=current_fract(:))
        dsl4jsb_Get_memory_tile(DISTURB_, disturbed_tile)
        dsl4jsb_Get_var2d_onChunk_tile_name(DISTURB_, damaged_fract, disturbed_tile)
        dsl4jsb_Get_memory_tile(PHENO_, disturbed_tile)
        dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, veg_fract_correction, disturbed_tile)
        dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, fract_fpc_max, disturbed_tile)

        lost_area(:,i_tile) = damaged_fract_disturbed_tile(:) * current_fract(:)
        gained_area(:,i_tile) = damaged_fract_disturbed_tile(:) * current_fract(:)

        disturbed_tile => disturbed_tile%Get_next_sibling_tile()
      END DO

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
        current_tile => tile%Get_first_child_tile()
        i_tile = 0
        DO WHILE (ASSOCIATED(current_tile))
          IF (ANY(current_tile%name .EQ. lcc_relocations%tile_names)) THEN
            i_tile = i_tile + 1
            IF (ANY(lost_area(:, i_tile) > 0.0_wp)) THEN
              ! 4a -- 'special transfer'
              IF (current_tile%Has_process_memory(CARBON_)) THEN
                ! treat vegetation pools: copy matter for diagnostics
                ! -> matter itself will be actively transferred to soils by the generic routine in the carbon interface
                conversion_fact = 1.0_wp / sec_per_day
                CALL copy_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, options, &
                  & 'c_woods', 'cflux_dist_woods_2_soil', conversion_fact = conversion_fact )
                CALL copy_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, options, &
                  & 'c_green', 'cflux_dist_green_2_soil', conversion_fact = conversion_fact )
                CALL copy_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, options, &
                  & 'c_reserve', 'cflux_dist_green_2_soil', conversion_fact = conversion_fact )
              ENDIF
              ! 4b -- general transfer of active vars
              CALL transfer_active_to_passive_onChunk(lcc_relocations, current_tile, i_tile, options)
            ENDIF
            !last_involved_tile => current_tile
          ENDIF

          current_tile => current_tile%Get_next_sibling_tile()
        ENDDO

        !>
        !> 5. Each lcc process needs to call end lcc to make the passive matter transfer
        !>
        CALL end_lcc_reloc(lcc_relocations, options, gained_area)
      ENDIF ! ... only to be conducted if there is any area movement
    ENDDO ! Number of lcc-cycles: depending on if are movement is to be "separate" or mixed

    !>
    !> 6. Deallocate area vectors
    !>
    DEALLOCATE(lost_area, gained_area, initial_area)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
    
  END SUBROUTINE update_wlcc
  
  ! ====================================================================================================== !
  !
  !> Implementation of "aggregate" for task "wlcc"
  !
  SUBROUTINE aggregate_wlcc(tile, options)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER  :: iblk
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_wlcc'
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk = options%iblk
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    
    !> Currently nothing to do
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_wlcc

#endif
END MODULE mo_wlcc_interface
