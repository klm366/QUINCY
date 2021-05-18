!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Contains the interface to the soil_biogeochemistry processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20170104)
!! 2nd version by Jan Engel (20200618), add to jsbach4
!!
MODULE mo_sb_interface
#ifndef __NO_JSBACH__

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process

  ! Use of processes in this module
  dsl4jsb_Use_processes SB_
  ! Use of process configurations
  dsl4jsb_Use_config(SB_)
  ! Use of process memories
  dsl4jsb_Use_memory(SB_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_sb_tasks

  ! ======================================================================================================= !
  !>Type definition: reset_sb_fluxes task
  !> 
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_sb_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_reset_sb_fluxes  !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_reset_sb_fluxes    !< Aggregates computed task variables
  END TYPE tsk_reset_sb_fluxes

  ! !-----------------------------------------------------------------------------------------------------
  ! !> Type definition: update_sb_rate_modifier task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_sb_rate_modifier
  ! CONTAINS
  !   PROCEDURE, NOPASS :: Integrate => update_sb_rate_modifier      !< Advances task computation for one timestep
  !   PROCEDURE, NOPASS :: Aggregate => aggr_update_sb_rate_modifier !< Aggregates computed task variables
  ! END TYPE tsk_update_sb_rate_modifier
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Type definition: update_sb_simple_model task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_sb_simple_model
  ! CONTAINS
  !   PROCEDURE, NOPASS :: Integrate => update_sb_simple_model      !< Advances task computation for one timestep
  !   PROCEDURE, NOPASS :: Aggregate => aggr_update_sb_simple_model !< Aggregates computed task variables
  ! END TYPE tsk_update_sb_simple_model
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Type definition: update_sb_jsm task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_sb_jsm
  ! CONTAINS
  !   PROCEDURE, NOPASS :: Integrate => update_sb_jsm      !< Advances task computation for one timestep
  !   PROCEDURE, NOPASS :: Aggregate => aggr_update_sb_jsm !< Aggregates computed task variables
  ! END TYPE tsk_update_sb_jsm

  ! ======================================================================================================= !
  !>Type definition: update_sb_pools task
  !> 
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_sb_pools
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sb_pools      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_update_sb_pools !< Aggregates computed task variables
  END TYPE tsk_update_sb_pools  

  ! ======================================================================================================= !
  !>Type definition: time_average_soilbiogeochemistry task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_time_average_soilbiogeochemistry
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_time_average_soilbiogeochemistry    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_time_average_soilbiogeochemistry !< Aggregates computed task variables
  END TYPE tsk_time_average_soilbiogeochemistry


  ! ======================================================================================================= !
  !>Constructor interface: reset_sb_fluxes task
  !> 
  INTERFACE tsk_reset_sb_fluxes
    PROCEDURE Create_task_reset_sb_fluxes         !< Constructor function for task
  END INTERFACE tsk_reset_sb_fluxes

  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor interface: update_sb_rate_modifier task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! INTERFACE tsk_update_sb_rate_modifier
  !   PROCEDURE Create_task_update_sb_rate_modifier         !< Constructor function for task
  ! END INTERFACE tsk_update_sb_rate_modifier
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor interface: update_sb_simple_model task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! INTERFACE tsk_update_sb_simple_model
  !   PROCEDURE Create_task_update_sb_simple_model         !< Constructor function for task
  ! END INTERFACE tsk_update_sb_simple_model
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor interface: update_sb_jsm task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! INTERFACE tsk_update_sb_jsm
  !   PROCEDURE Create_task_update_sb_jsm         !< Constructor function for task
  ! END INTERFACE tsk_update_sb_jsm

  ! ======================================================================================================= !
  !>Constructor interface: update_sb_pools task
  !> 
  INTERFACE tsk_update_sb_pools
    PROCEDURE Create_task_update_sb_pools         !< Constructor function for task
  END INTERFACE tsk_update_sb_pools

  ! ======================================================================================================= !
  !>Constructor interface: time_average_soilbiogeochemistry task
  !>
  INTERFACE tsk_time_average_soilbiogeochemistry
    PROCEDURE Create_task_update_time_average_soilbiogeochemistry         !< Constructor function for task
  END INTERFACE tsk_time_average_soilbiogeochemistry

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_interface'


CONTAINS


  ! ======================================================================================================= !
  !>Register tasks: SB_ process
  !> 
  SUBROUTINE Register_sb_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_reset_sb_fluxes(model_id))
    ! CALL this%Register_task(tsk_update_sb_rate_modifier(model_id))
    ! CALL this%Register_task(tsk_update_sb_simple_model(model_id))
    ! CALL this%Register_task(tsk_update_sb_jsm(model_id))
    CALL this%Register_task(tsk_update_sb_pools(model_id))
    CALL this%Register_task(tsk_time_average_soilbiogeochemistry(model_id))

  END SUBROUTINE Register_sb_tasks


  ! ======================================================================================================= !
  !>Constructor: reset_sb_fluxes task
  !>
  FUNCTION Create_task_reset_sb_fluxes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_reset_sb_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_sb_fluxes', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_sb_fluxes

  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor: update_sb_rate_modifier task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! FUNCTION Create_task_update_sb_rate_modifier(model_id) RESULT(return_ptr)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   INTEGER,                   INTENT(in) :: model_id
  !   CLASS(t_jsb_process_task), POINTER    :: return_ptr
  ! 
  ! 
  !   ALLOCATE(tsk_update_sb_rate_modifier::return_ptr)
  !   CALL return_ptr%Construct(name='update_sb_rate_modifier', process_id=SB_, owner_model_id=model_id)
  ! 
  ! END FUNCTION Create_task_update_sb_rate_modifier
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor: update_sb_simple_model task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! FUNCTION Create_task_update_sb_simple_model(model_id) RESULT(return_ptr)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   INTEGER,                   INTENT(in) :: model_id
  !   CLASS(t_jsb_process_task), POINTER    :: return_ptr
  ! 
  ! 
  !   ALLOCATE(tsk_update_sb_simple_model::return_ptr)
  !   CALL return_ptr%Construct(name='update_sb_simple_model', process_id=SB_, owner_model_id=model_id)
  ! 
  ! END FUNCTION Create_task_update_sb_simple_model
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor: update_sb_jsm task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! FUNCTION Create_task_update_sb_jsm(model_id) RESULT(return_ptr)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   INTEGER,                   INTENT(in) :: model_id
  !   CLASS(t_jsb_process_task), POINTER    :: return_ptr
  ! 
  ! 
  !   ALLOCATE(tsk_update_sb_jsm::return_ptr)
  !   CALL return_ptr%Construct(name='update_sb_jsm', process_id=SB_, owner_model_id=model_id)
  ! 
  ! END FUNCTION Create_task_update_sb_jsm

  ! ======================================================================================================= !
  !>Constructor: update_sb_pools task
  !> 
  FUNCTION Create_task_update_sb_pools(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_update_sb_pools::return_ptr)
    CALL return_ptr%Construct(name='update_sb_pools', process_id=SB_, owner_model_id=model_id)
    
  END FUNCTION Create_task_update_sb_pools

  ! ======================================================================================================= !
  !>Constructor: tsk_time_average_soilbiogeochemistry task
  !> 
  FUNCTION Create_task_update_time_average_soilbiogeochemistry(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_time_average_soilbiogeochemistry::return_ptr)
    CALL return_ptr%Construct(name='time_average_soilbiogeochemistry', process_id=SB_, owner_model_id=model_id)
    
  END FUNCTION Create_task_update_time_average_soilbiogeochemistry


  ! ------------------------------------------------------------------------------------------------------- !
  ! update_* routines
  ! ------------------------------------------------------------------------------------------------------- !
  

  ! ======================================================================================================= !
  !>Implementation of "integrate": reset_sb_fluxes task
  !>
  SUBROUTINE update_reset_sb_fluxes(tile, options)

    USE mo_sb_util, ONLY: real_reset_sb_fluxes => reset_sb_fluxes
    dsl4jsb_Use_processes SB_

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_reset_sb_fluxes'

    iblk = options%iblk
    IF (.NOT. tile%Is_process_active(SB_)) RETURN
    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL real_reset_sb_fluxes(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_reset_sb_fluxes

  ! ======================================================================================================= !
  !>Implementation of "integrate": update_sb_pools task
  !>
  SUBROUTINE update_sb_pools(tile, options)

    USE mo_sb_update_pools, ONLY: real_update_sb_pools => update_sb_pools
    dsl4jsb_Use_processes SB_

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk
    CHARACTER(len=*), PARAMETER :: routine = modname//':integrate_update_sb_pools'

    iblk = options%iblk
    IF (.NOT. tile%Is_process_active(SB_)) RETURN
    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL real_update_sb_pools(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_sb_pools

  ! ======================================================================================================= !
  !>calculate time moving averages and daytime averages for SB_
  !>
  SUBROUTINE update_time_average_soilbiogeochemistry(tile, options)

    USE mo_jsb_control,         ONLY: debug_on
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,      ONLY: t_jsb_task_options
    USE mo_sb_util,             ONLY: real_update_time_average_soilbiogeochemistry => update_time_average_soilbiogeochemistry
    dsl4jsb_Use_processes SB_
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_time_average_soilbiogeochemistry'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_active(SB_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)          RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. options%iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')


    CALL real_update_time_average_soilbiogeochemistry(tile, options)
    
    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_time_average_soilbiogeochemistry


  ! ------------------------------------------------------------------------------------------------------- !
  ! aggregate_* routines
  ! ------------------------------------------------------------------------------------------------------- !
  

  ! ======================================================================================================= !
  !>Implementation of "aggregate": reset_sb_fluxes task
  !>
  SUBROUTINE aggr_reset_sb_fluxes(tile, options)

    dsl4jsb_Use_processes SB_

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_reset_sb_fluxes'

    ! dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    ! Get process config
    !dsl4jsb_Get_config(SB_)
    ! Get process memories
    dsl4jsb_Get_memory(SB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    ! nothing to aggregate, yet

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE aggr_reset_sb_fluxes

  ! !-----------------------------------------------------------------------------------------------------
  ! !> Implementation of "aggregate": update_sb_rate_modifier task
  ! !!  actually NOT using this routine, i.e., no var is aggregated here
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE aggr_update_sb_rate_modifier(tile, options)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
  !   TYPE(t_jsb_task_options),   INTENT(in)    :: options
  !   ! ---------------------------
  !   ! 0.2 Local
  !   TYPE(t_jsb_model),        POINTER         :: model
  !   CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
  !   INTEGER                                   :: iblk, ics, ice, nc
  !   !X if necessary: REAL(wp) :: dtime, steplen
  !   CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_sb_rate_modifier'
  !   ! ---------------------------
  !   ! 0.3 Declare Memory
  !   ! Declare process configuration and memory Pointers
  !   dsl4jsb_Def_config(SB_)
  !   dsl4jsb_Def_memory(SB_)
  !   ! Get local variables from options argument
  !   iblk    = options%iblk
  !   ics     = options%ics
  !   ice     = options%ice
  !   nc      = options%nc
  !   !X if necessary: dtime   = options%dtime
  !   !X if necessary: steplen = options%steplen
  !   ! ---------------------------
  !   ! 0.4 Process Activity, Debug Option
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
  !   ! ---------------------------
  !   ! 0.5 Get Memory
  !   ! Get process variables (Set pointers to variables in memory)
  !   model => Get_model(tile%owner_model_id)
  !   ! Get process config
  !   dsl4jsb_Get_config(SB_)
  !   ! Get process memories
  !   dsl4jsb_Get_memory(SB_)
  ! 
  ! 
  !   ! ------------------------------------------------------------------------------------------------------------
  !   ! Go aggregate
  !   ! ------------------------------------------------------------------------------------------------------------
  ! 
  ! 
  !   IF (tile%Has_children()) THEN
  !     weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
  !     ! do not aggregate after this routine, only after update_sb_pools()
  !   END IF
  ! 
  !   !X Implementation: Start your process scheme here...
  !   !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.
  ! 
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  ! 
  ! END SUBROUTINE aggr_update_sb_rate_modifier
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Implementation of "aggregate": update_sb_simple_model task
  ! !!  actually NOT using this routine, i.e., no var is aggregated here
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE aggr_update_sb_simple_model(tile, options)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
  !   TYPE(t_jsb_task_options),   INTENT(in)    :: options
  !   ! ---------------------------
  !   ! 0.2 Local
  !   TYPE(t_jsb_model),        POINTER         :: model
  !   CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
  !   INTEGER                                   :: iblk, ics, ice, nc
  !   !X if necessary: REAL(wp) :: dtime, steplen
  !   CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_sb_simple_model'
  !   ! ---------------------------
  !   ! 0.3 Declare Memory
  !   ! Declare process configuration and memory Pointers
  !   dsl4jsb_Def_config(SB_)
  !   dsl4jsb_Def_memory(SB_)
  !   ! Get local variables from options argument
  !   iblk    = options%iblk
  !   ics     = options%ics
  !   ice     = options%ice
  !   nc      = options%nc
  !   !X if necessary: dtime   = options%dtime
  !   !X if necessary: steplen = options%steplen
  !   ! ---------------------------
  !   ! 0.4 Process Activity, Debug Option
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
  !   ! ---------------------------
  !   ! 0.5 Get Memory
  !   ! Get process variables (Set pointers to variables in memory)
  !   model => Get_model(tile%owner_model_id)
  !   ! Get process config
  !   dsl4jsb_Get_config(SB_)
  !   ! Get process memories
  !   dsl4jsb_Get_memory(SB_)
  ! 
  ! 
  !   ! ------------------------------------------------------------------------------------------------------------
  !   ! Go aggregate
  !   ! ------------------------------------------------------------------------------------------------------------
  ! 
  ! 
  !   IF (tile%Has_children()) THEN
  !     weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
  !     ! do not aggregate after this routine, only after update_sb_pools()
  !   END IF
  ! 
  !   !X Implementation: Start your process scheme here...
  !   !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.
  ! 
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  ! 
  ! END SUBROUTINE aggr_update_sb_simple_model
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Implementation of "aggregate": update_sb_jsm task
  ! !!  actually NOT using this routine, i.e., no var is aggregated here
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE aggr_update_sb_jsm(tile, options)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
  !   TYPE(t_jsb_task_options),   INTENT(in)    :: options
  !   ! ---------------------------
  !   ! 0.2 Local
  !   TYPE(t_jsb_model),        POINTER         :: model
  !   CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
  !   INTEGER                                   :: iblk, ics, ice, nc
  !   !X if necessary: REAL(wp) :: dtime, steplen
  !   CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_sb_jsm'
  !   ! ---------------------------
  !   ! 0.3 Declare Memory
  !   ! Declare process configuration and memory Pointers
  !   dsl4jsb_Def_config(SB_)
  !   dsl4jsb_Def_memory(SB_)
  !   ! Get local variables from options argument
  !   iblk    = options%iblk
  !   ics     = options%ics
  !   ice     = options%ice
  !   nc      = options%nc
  !   !X if necessary: dtime   = options%dtime
  !   !X if necessary: steplen = options%steplen
  !   ! ---------------------------
  !   ! 0.4 Process Activity, Debug Option
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
  !   ! ---------------------------
  !   ! 0.5 Get Memory
  !   ! Get process variables (Set pointers to variables in memory)
  !   model => Get_model(tile%owner_model_id)
  !   ! Get process config
  !   dsl4jsb_Get_config(SB_)
  !   ! Get process memories
  !   dsl4jsb_Get_memory(SB_)
  ! 
  ! 
  !   ! ------------------------------------------------------------------------------------------------------------
  !   ! Go aggregate
  !   ! ------------------------------------------------------------------------------------------------------------
  ! 
  ! 
  !   IF (tile%Has_children()) THEN
  !     weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
  !     ! do not aggregate after this routine, only after update_sb_pools()
  !   END IF
  ! 
  !   !X Implementation: Start your process scheme here...
  !   !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.
  ! 
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  ! 
  ! END SUBROUTINE aggr_update_sb_jsm


  ! ======================================================================================================= !
  !>Implementation of "aggregate": update_sb_pools task
  !> 
  SUBROUTINE aggr_update_sb_pools(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_sb_pools'

    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)

    ! ------------------------------------------------------------------------------------------------------------
    ! Go aggregate
    ! ------------------------------------------------------------------------------------------------------------
    
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    ! dsl4jsb_Aggregate_onChunk(SB_,       total_flux_carbon          , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_,       total_flux_nitrogen        , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_,       total_flux_phosphorus      , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_,       total_flux_carbon13        , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_,       total_flux_carbon14        , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_,       total_flux_nitrogen15      , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_update_sb_pools


  ! ======================================================================================================= !
  !> Implementation of "aggregate": time_average_assimilation task
  !>
  SUBROUTINE aggregate_time_average_soilbiogeochemistry(tile, options)

    USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
    USE mo_jsb_task_class,     ONLY: t_jsb_task_options
    dsl4jsb_Use_processes SB_
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggregate_time_average_soilbiogeochemistry'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk      :: microbial_cue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk      :: microbial_nue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk      :: microbial_pue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk      :: myc_export_c_tmyc_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_n_tmyc_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_n_tlabile_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_p_tlabile_mavg_sl
    dsl4jsb_Real3D_onChunk      :: dom_cn_mavg
    dsl4jsb_Real3D_onChunk      :: dom_cp_mavg
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_c_mavg
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_n_mavg
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_p_mavg
    ! ----------------------------------------------------------------------------------------------------- !
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! SB_ 3D
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_cue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_nue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_pue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_c_tmyc_mavg_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_n_tmyc_mavg_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_n_tlabile_mavg_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_p_tlabile_mavg_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cn_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cp_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_c_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_n_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_p_mavg)


    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(SB_, microbial_cue_eff_tmic_mavg    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, microbial_nue_eff_tmic_mavg    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, microbial_pue_eff_tmic_mavg    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_c_tmyc_mavg_sl      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_n_tmyc_mavg_sl      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_n_tlabile_mavg_sl   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_p_tlabile_mavg_sl   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, dom_cn_mavg                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, dom_cp_mavg                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, enzyme_frac_poly_c_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, enzyme_frac_poly_n_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, enzyme_frac_poly_p_mavg        , weighted_by_fract)


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_time_average_soilbiogeochemistry

#endif
END MODULE mo_sb_interface
