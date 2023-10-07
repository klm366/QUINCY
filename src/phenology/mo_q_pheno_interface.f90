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
!! @brief Contains the interface to the phenology processes of the quincy model
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20161214)
!! Second version by Steffen Richter (June 2018)
!! Jan Engel MPI-BGC (2023-09-26): re-write for quincy model based on the jsbach module
!!

MODULE mo_q_pheno_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_control,        ONLY: debug_on
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes PHENO_, VEG_
  dsl4jsb_Use_config(PHENO_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(VEG_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_pheno_tasks_quincy
  PUBLIC :: update_phenology

  ! ======================================================================================================= !
  !> Type definition for the phenology task
  !>
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_phenology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_phenology    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_phenology !< Aggregates computed task variables
  END TYPE   tsk_phenology
  !> Constructor interface for update_phenology task
  !>
  INTERFACE   tsk_phenology
    PROCEDURE Create_task_phenology         !< Constructor function for task
  END INTERFACE   tsk_phenology

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_pheno_interface'

CONTAINS

  ! ======================================================================================================= !
  !> Constructor for phenology task
  !>
  FUNCTION Create_task_phenology(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_phenology::return_ptr)
    CALL return_ptr%Construct(name='phenology', process_id=PHENO_, owner_model_id=model_id)
  END FUNCTION Create_task_phenology

  ! ======================================================================================================= !
  !> Register tasks for pheno process
  !>
  SUBROUTINE Register_pheno_tasks_quincy(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_phenology(model_id))
  END SUBROUTINE Register_pheno_tasks_quincy

  ! ======================================================================================================= !
  !> Implementation of "update" of the "phenology" task
  !>
  !>
  SUBROUTINE update_phenology(tile, options)
    USE mo_q_pheno_process,         ONLY: update_phenology_quincy
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_phenology'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(PHENO_)) RETURN

    CALL update_phenology_quincy(tile, options)
  END SUBROUTINE update_phenology

  ! ======================================================================================================= !
  !> Implementation of "aggregate" of the "phenology" task
  !>
  SUBROUTINE aggregate_phenology(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_phenology'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(VEG_)

    ! PHENO_
    dsl4jsb_Aggregate_onChunk(PHENO_, growing_season              ,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(PHENO_, gdd                         ,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(PHENO_, nd_dormance                 ,    weighted_by_fract)
    ! VEG_
    dsl4jsb_Aggregate_onChunk(VEG_,   growth_req_n_tlabile_mavg   ,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_,   growth_req_p_tlabile_mavg   ,    weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_phenology

#endif
END MODULE mo_q_pheno_interface
