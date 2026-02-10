!> QUINCY plant hydraulics process interface
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### definition and init of tasks for the plant hydraulic process incl. update and aggregate routines
!>
MODULE mo_phyd_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_control,        ONLY: debug_on
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_phyd_process,     ONLY: update_plant_hydraulics

  dsl4jsb_Use_processes PHYD_
  dsl4jsb_Use_config(PHYD_)
  dsl4jsb_Use_memory(PHYD_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_phyd_tasks

  ! ======================================================================================================= !
  !> Type definition for the plant_hydraulics task
  !>
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_plant_hydraulics
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_plant_hydraulics    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_plant_hydraulics !< Aggregates computed task variables
  END TYPE   tsk_plant_hydraulics
  !> Constructor interface for update_plant_hydraulics task
  !>
  INTERFACE tsk_plant_hydraulics
    PROCEDURE Create_task_plant_hydraulics         !< Constructor function for task
  END INTERFACE tsk_plant_hydraulics

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pyhdro_interface'

CONTAINS

  ! ======================================================================================================= !
  !> Constructor for phenology task
  !>
  FUNCTION Create_task_plant_hydraulics(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_plant_hydraulics::return_ptr)
    CALL return_ptr%Construct(name='plant_hydraulics', process_id=PHYD_, owner_model_id=model_id)
  END FUNCTION Create_task_plant_hydraulics

  ! ======================================================================================================= !
  !> Register tasks for pheno process
  !>
  SUBROUTINE Register_phyd_tasks(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_plant_hydraulics(model_id))
  END SUBROUTINE Register_phyd_tasks

  ! ======================================================================================================= !
  !> Implementation of "aggregate" of the "plant_hydraulics" task
  !>
  SUBROUTINE aggregate_plant_hydraulics(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(PHYD_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_plant_hydraulics'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(PHYD_)

    ! PHYD_
    dsl4jsb_Aggregate_onChunk(PHYD_,   psi_leaf   ,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(PHYD_,   psi_stem   ,    weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_plant_hydraulics

#endif
END MODULE mo_phyd_interface
