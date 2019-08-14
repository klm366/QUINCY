!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., and Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Contains the interface to the sse processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-19)
!!
MODULE mo_sse_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process, &
                                            A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
  USE mo_sse_util,                    ONLY: reset_sse_fluxes
  USE mo_sse_update_turbulence,       ONLY: update_soil_turbulence
  USE mo_sse_update_physics,          ONLY: update_sse_physics

  ! Use of process configurations (t_PROC_config)
  dsl4jsb_Use_config(SSE_)

  ! Use of process memories (t_PROC_memory)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(L2A_)
  

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_sse_tasks, reset_sse_fluxes, update_soil_turbulence, update_sse_physics

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: reset_sse_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_sse_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => reset_sse_fluxes      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_reset_sse_fluxes !< Aggregates computed task variables
  END TYPE tsk_reset_sse_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_soil_turbulence task
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_soil_turbulence
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_turbulence      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_update_soil_turbulence !< Aggregates computed task variables
  END TYPE tsk_update_soil_turbulence

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_sse_physics task
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_sse_physics
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sse_physics      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_update_sse_physics !< Aggregates computed task variables
  END TYPE tsk_update_sse_physics


  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: reset_sse_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_reset_sse_fluxes
    PROCEDURE Create_task_reset_sse_fluxes         !< Constructor function for task
  END INTERFACE tsk_reset_sse_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_soil_turbulence task
  !! 
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_update_soil_turbulence
    PROCEDURE Create_task_update_soil_turbulence         !< Constructor function for task
  END INTERFACE tsk_update_soil_turbulence

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_sse_physics task
  !! 
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_update_sse_physics
    PROCEDURE Create_task_update_sse_physics         !< Constructor function for task
  END INTERFACE tsk_update_sse_physics


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_interface'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: SSE_ process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Register_sse_tasks(this, model_id)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id


    CALL this%Register_task(tsk_reset_sse_fluxes(model_id))
    CALL this%Register_task(tsk_update_soil_turbulence(model_id))
    CALL this%Register_task(tsk_update_sse_physics(model_id))

  END SUBROUTINE Register_sse_tasks


  !-----------------------------------------------------------------------------------------------------
  !> Constructor: reset_sse_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_reset_sse_fluxes(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr


    ALLOCATE(tsk_reset_sse_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_sse_fluxes', process_id=SSE_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_sse_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_soil_turbulence task
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_soil_turbulence(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr


    ALLOCATE(tsk_update_soil_turbulence::return_ptr)
    CALL return_ptr%Construct(name='update_soil_turbulence', process_id=SSE_, owner_model_id=model_id)
    
  END FUNCTION Create_task_update_soil_turbulence

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_sse_physics task
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_sse_physics(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr


    ALLOCATE(tsk_update_sse_physics::return_ptr)
    CALL return_ptr%Construct(name='update_sse_physics', process_id=SSE_, owner_model_id=model_id)
    
  END FUNCTION Create_task_update_sse_physics


  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": reset_sse_fluxes task
  !!  actually NOT using this routine, i.e., no var is aggregated here
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggr_reset_sse_fluxes(tile, options)

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
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_reset_sse_fluxes'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
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
    dsl4jsb_Get_config(SSE_)
    ! Get process memories
    dsl4jsb_Get_memory(SSE_)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go aggregate
    ! ------------------------------------------------------------------------------------------------------------
    

    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
      ! do not aggregate after this routine
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen pheno Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_reset_sse_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": update_soil_turbulence task
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggr_update_soil_turbulence(tile, options)

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
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_soil_turbulence'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(L2A_)
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
    dsl4jsb_Get_config(SSE_)
    ! Get process memories
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(L2A_)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go aggregate
    ! ------------------------------------------------------------------------------------------------------------
    

    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

      dsl4jsb_Aggregate_onChunk(A2L_, drag_srf             , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(A2L_, ga                   , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(A2L_, pch                  , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(A2L_, t_acoef              , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(A2L_, t_bcoef              , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(A2L_, q_acoef              , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(A2L_, q_bcoef              , weighted_by_fract)
      ! ---------------------------
      dsl4jsb_Aggregate_onChunk(L2A_, net_radiation        , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(L2A_, ground_heat_flx      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(L2A_, sensible_heat_flx    , weighted_by_fract)
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_update_soil_turbulence

  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": update_sse_physics task
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggr_update_sse_physics(tile, options)

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
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_sse_physics'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(L2A_)
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
    dsl4jsb_Get_config(SSE_)
    ! Get process memories
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(L2A_)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go aggregate
    ! ------------------------------------------------------------------------------------------------------------
    

    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

      dsl4jsb_Aggregate_onChunk(A2L_, ga                    , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(A2L_, pch                   , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(A2L_, drag_srf              , weighted_by_fract)
      ! ---------------------------
      dsl4jsb_Aggregate_onChunk(L2A_, ground_heat_flx_old   , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(L2A_, net_radiation         , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(L2A_, ground_heat_flx       , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(L2A_, latent_heat_flx       , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(L2A_, sensible_heat_flx     , weighted_by_fract)
      ! ---------------------------
      dsl4jsb_Aggregate_onChunk(SSE_, t_srf_old             , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, zril_old              , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, transpiration         , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, w_soil_root_theta     , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, w_soil_root_pot       , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, w_skin                , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, w_soil_root           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, t_soil_root           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, pet                   , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, interception          , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, evaporation           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, runoff                , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, drainage              , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, drainage_fraction     , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, w_soil_sl             , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, w_soil_pot_sl         , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, percolation_sl        , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, drainage_sl           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_, t_soil_sl             , weighted_by_fract)
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_update_sse_physics


#endif
END MODULE mo_sse_interface
