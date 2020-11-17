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
!! @brief Contains the interface to the assimilation processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-19)
!!
MODULE mo_assimi_interface
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
  USE mo_assimi_util,                 ONLY: reset_assimi_fluxes
  USE mo_assimi_update_canopy_fluxes, ONLY: update_canopy_fluxes

  ! Use of process configurations (t_PROC_config)
  dsl4jsb_Use_config(ASSIMI_)

  ! Use of process memories (t_PROC_memory)
  dsl4jsb_Use_memory(ASSIMI_)
  

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_assimi_tasks, reset_assimi_fluxes, update_canopy_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: reset_assimi_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_assimi_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => reset_assimi_fluxes      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_reset_assimi_fluxes !< Aggregates computed task variables
  END TYPE tsk_reset_assimi_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_canopy_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_canopy_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_fluxes      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_update_canopy_fluxes !< Aggregates computed task variables
  END TYPE tsk_update_canopy_fluxes


  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: reset_assimi_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_reset_assimi_fluxes
    PROCEDURE Create_task_reset_assimi_fluxes         !< Constructor function for task
  END INTERFACE tsk_reset_assimi_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_canopy_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_update_canopy_fluxes
    PROCEDURE Create_task_update_canopy_fluxes         !< Constructor function for task
  END INTERFACE tsk_update_canopy_fluxes


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_assimi_interface'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: ASSIMI_ process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Register_assimi_tasks(this, model_id)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id


    CALL this%Register_task(tsk_reset_assimi_fluxes(model_id))
    CALL this%Register_task(tsk_update_canopy_fluxes(model_id))

  END SUBROUTINE Register_assimi_tasks


  !-----------------------------------------------------------------------------------------------------
  !> Constructor: reset_assimi_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_reset_assimi_fluxes(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr


    ALLOCATE(tsk_reset_assimi_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_assimi_fluxes', process_id=ASSIMI_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_assimi_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_canopy_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_canopy_fluxes(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr


    ALLOCATE(tsk_update_canopy_fluxes::return_ptr)
    CALL return_ptr%Construct(name='update_canopy_fluxes', process_id=ASSIMI_, owner_model_id=model_id)
    
  END FUNCTION Create_task_update_canopy_fluxes


  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": reset_assimi_fluxes task
  !!  actually NOT using this routine, i.e., no var is aggregated here
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggr_reset_assimi_fluxes(tile, options)

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
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_reset_assimi_fluxes'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_memory(ASSIMI_)
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
    dsl4jsb_Get_config(ASSIMI_)
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)


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

  END SUBROUTINE aggr_reset_assimi_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": update_canopy_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggr_update_canopy_fluxes(tile, options)

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
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_canopy_fluxes'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_memory(ASSIMI_)
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
    dsl4jsb_Get_config(ASSIMI_)
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go aggregate
    ! ------------------------------------------------------------------------------------------------------------
    

    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

      dsl4jsb_Aggregate_onChunk(ASSIMI_, gross_assimilation          , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, gross_assimilation_C13      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, gross_assimilation_C14      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, net_assimilation            , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, net_assimilation_boc        , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, maint_respiration_leaf      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, canopy_cond                 , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, co2_conc_leaf               , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, beta_air                    , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, beta_soa                    , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, soa_tsoa_mavg               , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, beta_soil_gs                , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, beta_soil_ps                , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, t_jmax_opt                  , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, net_assimilation_cl         , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, gross_assimilation_cl       , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, maint_respiration_leaf_cl   , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, canopy_cond_cl              , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, co2_conc_leaf_cl            , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, jmax_cl                     , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(ASSIMI_, vcmax_cl                    , weighted_by_fract)
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_update_canopy_fluxes

#endif
END MODULE mo_assimi_interface
