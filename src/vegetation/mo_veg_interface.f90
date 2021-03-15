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
!! @brief Contains the interface for plant growth and turnover processes 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!! 2nd version by Jan Engel (20160719)
!!
MODULE mo_veg_interface
#ifndef __NO_JSBACH__


  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process
  USE mo_veg_turnover,                ONLY: update_veg_turnover
  USE mo_veg_dynamics,                ONLY: update_veg_dynamics
  USE mo_veg_growth,                  ONLY: update_veg_growth
  USE mo_veg_update_pools,            ONLY: update_veg_pools

  ! Use of processes in this module (USE mo_jsb_process_class, ONLY: )
  dsl4jsb_Use_processes VEG_

  ! Use of process configurations (t_PROC_config)
  dsl4jsb_Use_config(VEG_)

  ! Use of process memories (t_PROC_memory)
  dsl4jsb_Use_memory(VEG_)
  

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_veg_tasks


  !-----------------------------------------------------------------------------------------------------
  !> Type definition: reset_veg_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_veg_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_reset_veg_fluxes    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_reset_veg_fluxes !< Aggregates computed task variables
  END TYPE tsk_reset_veg_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_veg_turnover task
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_veg_turnover
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_veg_turnover      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_update_veg_turnover !< Aggregates computed task variables
  END TYPE tsk_update_veg_turnover

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_veg_dynamics task
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_veg_dynamics
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_veg_dynamics      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_update_veg_dynamics !< Aggregates computed task variables
  END TYPE tsk_update_veg_dynamics

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_veg_growth task
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_veg_growth
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_veg_growth      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_update_veg_growth !< Aggregates computed task variables
  END TYPE tsk_update_veg_growth

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_veg_pools task
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_veg_pools
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_veg_pools      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_update_veg_pools !< Aggregates computed task variables
  END TYPE tsk_update_veg_pools  

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: time_average_vegetation task
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_time_average_vegetation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_time_average_vegetation    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_time_average_vegetation !< Aggregates computed task variables
  END TYPE tsk_time_average_vegetation


  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: reset_veg_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_reset_veg_fluxes
    PROCEDURE Create_task_reset_veg_fluxes         !< Constructor function for task
  END INTERFACE tsk_reset_veg_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_veg_turnover task
  !! 
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_update_veg_turnover
    PROCEDURE Create_task_update_veg_turnover         !< Constructor function for task
  END INTERFACE tsk_update_veg_turnover

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_veg_dynamics task
  !! 
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_update_veg_dynamics
    PROCEDURE Create_task_update_veg_dynamics         !< Constructor function for task
  END INTERFACE tsk_update_veg_dynamics

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_veg_growth task
  !! 
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_update_veg_growth
    PROCEDURE Create_task_update_veg_growth         !< Constructor function for task
  END INTERFACE tsk_update_veg_growth

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_veg_pools task
  !! 
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_update_veg_pools
    PROCEDURE Create_task_update_veg_pools         !< Constructor function for task
  END INTERFACE tsk_update_veg_pools

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_time_average_vegetation task
  !! 
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_time_average_vegetation
    PROCEDURE Create_task_update_time_average_vegetation         !< Constructor function for task
  END INTERFACE tsk_time_average_vegetation

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_veg_interface'


CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: VEG_ process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Register_veg_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_reset_veg_fluxes(model_id))
    CALL this%Register_task(tsk_update_veg_turnover(model_id))
    CALL this%Register_task(tsk_update_veg_dynamics(model_id))
    CALL this%Register_task(tsk_update_veg_growth(model_id))
    CALL this%Register_task(tsk_update_veg_pools(model_id))
    CALL this%Register_task(tsk_time_average_vegetation(model_id))

  END SUBROUTINE Register_veg_tasks


  !-----------------------------------------------------------------------------------------------------
  !> Constructor: reset_veg_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_reset_veg_fluxes(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr


    ALLOCATE(tsk_reset_veg_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_veg_fluxes', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_veg_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_veg_turnover task
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_veg_turnover(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr


    ALLOCATE(tsk_update_veg_turnover::return_ptr)
    CALL return_ptr%Construct(name='update_veg_turnover', process_id=VEG_, owner_model_id=model_id)
    
  END FUNCTION Create_task_update_veg_turnover

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_veg_dynamics task
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_veg_dynamics(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr


    ALLOCATE(tsk_update_veg_dynamics::return_ptr)
    CALL return_ptr%Construct(name='update_veg_dynamics', process_id=VEG_, owner_model_id=model_id)
    
  END FUNCTION Create_task_update_veg_dynamics

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_veg_growth task
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_veg_growth(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr


    ALLOCATE(tsk_update_veg_growth::return_ptr)
    CALL return_ptr%Construct(name='update_veg_growth', process_id=VEG_, owner_model_id=model_id)
    
  END FUNCTION Create_task_update_veg_growth

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_veg_pools task
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_veg_pools(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr


    ALLOCATE(tsk_update_veg_pools::return_ptr)
    CALL return_ptr%Construct(name='update_veg_pools', process_id=VEG_, owner_model_id=model_id)
    
  END FUNCTION Create_task_update_veg_pools

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: tsk_time_average_vegetation task
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_time_average_vegetation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_time_average_vegetation::return_ptr)
    CALL return_ptr%Construct(name='tsk_time_average_vegetation', process_id=VEG_, owner_model_id=model_id)
    
  END FUNCTION Create_task_update_time_average_vegetation
  

  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "update": reset_veg_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_reset_veg_fluxes(tile, options)
  
    USE mo_veg_util, ONLY: reset_veg_fluxes
  
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
  
    INTEGER :: iblk
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_reset_veg_fluxes'
  
    iblk    = options%iblk
  
    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')
  
    CALL reset_veg_fluxes(tile, options)
  
    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  
  END SUBROUTINE update_reset_veg_fluxes
  
  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": reset_veg_fluxes task
  !!  actually NOT using this routine, i.e., no var is aggregated here
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggregate_reset_veg_fluxes(tile, options)
  
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
  
    TYPE(t_jsb_model),        POINTER :: model
    CLASS(t_jsb_aggregator),  POINTER :: weighted_by_fract
    INTEGER                           :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_reset_veg_fluxes'
  
    ! Declare process configuration and memory Pointers
    !dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
      IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')

      model => Get_model(tile%owner_model_id)

      ! Get process config
      !dsl4jsb_Get_config(VEG_)
      ! Get process memories
      dsl4jsb_Get_memory(VEG_)

      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

      ! nothing to aggregate, yet

      IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE aggregate_reset_veg_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": update_veg_turnover task
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggr_update_veg_turnover(tile, options)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_veg_turnover'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process variables (Set pointers to variables in memory)
    model => Get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(VEG_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go aggregate
    ! ------------------------------------------------------------------------------------------------------------
    

    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
      ! do not aggregate after this routine, only after update_veg_pools()
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_update_veg_turnover

  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": update_veg_dynamics task
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggr_update_veg_dynamics(tile, options)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_veg_dynamics'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process variables (Set pointers to variables in memory)
    model => Get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(VEG_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go aggregate
    ! ------------------------------------------------------------------------------------------------------------
    

    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
      ! do not aggregate after this routine, only after update_veg_pools()
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_update_veg_dynamics

  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": update_veg_growth task
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggr_update_veg_growth(tile, options)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_veg_growth'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process variables (Set pointers to variables in memory)
    model => Get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(VEG_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go aggregate
    ! ------------------------------------------------------------------------------------------------------------
    

    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
      ! do not aggregate after this routine, only after update_veg_pools()
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_update_veg_growth

  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": update_veg_pools task
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggr_update_veg_pools(tile, options)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_veg_pools'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process variables (Set pointers to variables in memory)
    model => Get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(VEG_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go aggregate
    ! ------------------------------------------------------------------------------------------------------------
    

    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

      ! *slmdev 
      ! the below list is not complete, it is just an initial placeholder for all var that need to be aggregated !
      dsl4jsb_Aggregate_onChunk(VEG_, diameter            , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, height              , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, mean_leaf_age       , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, dphi                , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, npp                 , weighted_by_fract)
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_update_veg_pools


  ! ======================================================================================================= !
  !> calculate time moving averages and daytime averages for VEG_
  !>
  SUBROUTINE update_time_average_vegetation(tile, options)

    USE mo_jsb_control,         ONLY: debug_on
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,      ONLY: t_jsb_task_options
    USE mo_veg_util,            ONLY: real_update_time_average_vegetation => update_time_average_vegetation
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_time_average_vegetation'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_active(VEG_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)           RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. options%iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')


    CALL real_update_time_average_vegetation(tile, options)
    
    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_time_average_vegetation


  ! ======================================================================================================= !
  !> Implementation of "aggregate": time_average_vegetation task
  !>
  SUBROUTINE aggregate_time_average_vegetation(tile, options)

    USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
    USE mo_jsb_task_class,     ONLY: t_jsb_task_options
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggregate_time_average_vegetation'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! VEG 2D
    dsl4jsb_Real2D_onChunk      :: daytime_counter
    dsl4jsb_Real2D_onChunk      :: t_air_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: press_srf_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: ga_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: t_air_daytime
    dsl4jsb_Real2D_onChunk      :: press_srf_daytime
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_daytime
    dsl4jsb_Real2D_onChunk      :: ga_daytime
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_daytime
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_daytime
    dsl4jsb_Real2D_onChunk      :: gpp_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: maint_respiration_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_p_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: npp_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: demand_uptake_n_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: demand_uptake_p_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_p_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: press_srf_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: ga_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: press_srf_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: ga_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: uptake_n_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: fmaint_rate_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cn_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: growth_np_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: npp_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: labile_carbon_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: labile_nitrogen_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: labile_phosphorus_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: cost_n_fixation_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: npp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: n_fixation_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_npp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: transpiration_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_transpiration_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: dphi_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_p_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_carbon_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_carbon_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: w_root_lim_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cn_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_np_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: exudation_c_tmyc_mavg
    dsl4jsb_Real2D_onChunk      :: leaf2root_troot_mavg
    dsl4jsb_Real2D_onChunk      :: unit_npp_troot_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_pot_troot_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_p_pot_troot_mavg
    dsl4jsb_Real2D_onChunk      :: fmaint_rate_troot_mavg
    dsl4jsb_Real2D_onChunk      :: an_boc_tvegdyn_mavg
    dsl4jsb_Real2D_onChunk      :: net_growth_tvegdyn_mavg
    dsl4jsb_Real2D_onChunk      :: lai_tvegdyn_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tacclim_mavg
    dsl4jsb_Real2D_onChunk      :: t_soil_root_tacclim_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_week_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_month_mavg
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    ! VEG 3D
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_daytime_dacc_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_daytime_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tcnl_mavg_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tfrac_mavg_cl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! VEG 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, daytime_counter)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, gpp_tlabile_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_tlabile_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tlabile_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tlabile_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tphen_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, demand_uptake_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, demand_uptake_p_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_n_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, fmaint_rate_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cn_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_np_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_carbon_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_nitrogen_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_phosphorus_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, cost_n_fixation_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_npp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, dphi_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_carbon_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_nitrogen_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_phosphorus_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_carbon_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_nitrogen_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_phosphorus_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, w_root_lim_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cn_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_np_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, exudation_c_tmyc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, leaf2root_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_npp_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_pot_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_pot_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, fmaint_rate_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, an_boc_tvegdyn_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, net_growth_tvegdyn_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, lai_tvegdyn_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tacclim_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_soil_root_tacclim_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_week_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_month_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)
    ! VEG 3D
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_daytime_dacc_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_daytime_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tcnl_mavg_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl)


    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(VEG_, daytime_counter        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_daytime_dacc         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, press_srf_daytime_dacc         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, co2_mixing_ratio_daytime_dacc        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, ga_daytime_dacc        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, beta_sinklim_ps_daytime_dacc         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_jmax_opt_daytime_dacc        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_daytime        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, press_srf_daytime        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, co2_mixing_ratio_daytime         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, ga_daytime         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, beta_sinklim_ps_daytime        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_jmax_opt_daytime         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, gpp_tlabile_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, maint_respiration_tlabile_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_req_n_tlabile_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_req_p_tlabile_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_tphen_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, npp_tuptake_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, demand_uptake_n_tuptake_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, demand_uptake_p_tuptake_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_req_n_tuptake_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_req_p_tuptake_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_tfrac_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, press_srf_tfrac_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, co2_mixing_ratio_tfrac_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, ga_tfrac_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, beta_sinklim_ps_tfrac_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_tcnl_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, press_srf_tcnl_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, co2_mixing_ratio_tcnl_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, ga_tcnl_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, beta_sinklim_ps_tcnl_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_n_tcnl_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fmaint_rate_tcnl_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_cn_tcnl_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_np_tcnl_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, npp_tcnl_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, labile_carbon_tcnl_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, labile_nitrogen_tcnl_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, labile_phosphorus_tcnl_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, cost_n_fixation_talloc_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, npp_talloc_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_fixation_talloc_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_npp_talloc_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, transpiration_talloc_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_transpiration_talloc_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, dphi_talloc_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_uptake_n_talloc_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_uptake_p_talloc_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, labile_carbon_talloc_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, labile_nitrogen_talloc_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, labile_phosphorus_talloc_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, reserve_carbon_talloc_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, reserve_nitrogen_talloc_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, reserve_phosphorus_talloc_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, w_root_lim_talloc_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_cn_talloc_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_cp_talloc_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_np_talloc_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, exudation_c_tmyc_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, leaf2root_troot_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_npp_troot_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_uptake_n_pot_troot_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_uptake_p_pot_troot_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fmaint_rate_troot_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, an_boc_tvegdyn_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, net_growth_tvegdyn_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, lai_tvegdyn_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_tacclim_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_soil_root_tacclim_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_week_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_month_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_jmax_opt_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fleaf_sunlit_daytime_dacc_cl         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fleaf_sunlit_daytime_cl        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fleaf_sunlit_tcnl_mavg_cl        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl         , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_time_average_vegetation

#endif
END MODULE mo_veg_interface
