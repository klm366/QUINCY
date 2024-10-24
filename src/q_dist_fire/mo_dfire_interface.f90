!>
!>  @par Copyright
!>    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!>    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!>    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!>
!>    This file is part of the QUINCY model.
!>    The QUINCY model is free software: you can redistribute it and/or modify
!>    it under the terms of the GNU General Public License version 3,
!>    as published by the Free Software Foundation.
!>
!>    The QUINCY model is distributed in the hope that it will be useful,
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!>    GNU General Public License for more details.
!>
!>    You should have received a copy of the GNU General Public License
!>    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!> 
!> @brief Contains the interface to the dfire processes
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_dfire_interface
#ifndef __NO_JSBACH__

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process, DFIRE_
  USE mo_dfire_spitfire,              ONLY: update_dfire_spitfire

  dsl4jsb_Use_memory(DFIRE_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_dfire_tasks
  PUBLIC :: update_time_average_dist_fire, update_dfire_spitfire


  ! ======================================================================================================= !
  !>
  !> Type definition: dfire_sf_fuel task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_dfire_spitfire
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_dfire_spitfire      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_dfire_spitfire        !< Aggregates computed task variables
  END TYPE tsk_dfire_spitfire

  ! ======================================================================================================= !
  !>
  !> Constructor interface: dfire_sf_fuel task
  !>
  INTERFACE tsk_dfire_spitfire
    PROCEDURE Create_task_dfire_spitfire         !< Constructor function for task
  END INTERFACE tsk_dfire_spitfire

  ! ======================================================================================================= !
  !>
  !> Constructor interface: time_average_dfire task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_time_average_dist_fire
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_time_average_dist_fire    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_time_average_dist_fire !< Aggregates computed task variables
  END TYPE tsk_time_average_dist_fire

  ! ======================================================================================================= !
  !>
  !> Constructor interface: update_time__average_dfire task
  !>
  INTERFACE tsk_time_average_dist_fire
    PROCEDURE Create_task_update_time_average_dist_fire         !< Constructor function for task
  END INTERFACE tsk_time_average_dist_fire

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_interface'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Register tasks: DFIRE_ process
  !>
  SUBROUTINE Register_dfire_tasks(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_dfire_spitfire(model_id))
  END SUBROUTINE Register_dfire_tasks

  ! ======================================================================================================= !
  !>
  !> Constructor: dfire_sf_fuel task
  !>
  FUNCTION Create_task_dfire_spitfire(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_dfire_spitfire::return_ptr)
    CALL return_ptr%Construct(name='dfire_spitfire', process_id=DFIRE_, owner_model_id=model_id)
  END FUNCTION Create_task_dfire_spitfire

  ! ======================================================================================================= !
  !>
  !> Implementation of "aggregate": dfire_spitfire task
  !>
  SUBROUTINE aggr_dfire_spitfire(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER         :: model
    INTEGER                                   :: iblk, ics, ice, nc
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_dfire_sf_fuel'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(DFIRE_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go aggregate
    ! ----------------------------------------------------------------------------------------------------- !

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggr_dfire_spitfire


  ! ======================================================================================================= !
  !>
  !> calculate time moving averages and daytime averages for DFIRE_
  !>
  FUNCTION Create_task_update_time_average_dist_fire(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_time_average_dist_fire::return_ptr)
    CALL return_ptr%Construct(name='tavrg_dist_fire', process_id=DFIRE_, owner_model_id=model_id)

  END FUNCTION Create_task_update_time_average_dist_fire


  ! ======================================================================================================= !
  !>
  !> calculate time moving averages and daytime averages for DFIRE_
  !>
  SUBROUTINE update_time_average_dist_fire(tile, options)

    USE mo_jsb_control,         ONLY: debug_on
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,      ONLY: t_jsb_task_options
!    USE mo_dfire_util,          ONLY: calculate_time_average_dist_fire,  &
!                                      calculate_daily_weather, &
!                                      calculate_pure_moving_averages

    dsl4jsb_Use_processes DFIRE_
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_time_average_dist_fire'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(DFIRE_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)           RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. options%iblk == 1) &
            & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

!    CALL calculate_time_average_dist_fire(tile, options)
!    CALL calculate_daily_weather(tile, options)
!    CALL calculate_pure_moving_averages(tile, options)

    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_time_average_dist_fire


  ! ======================================================================================================= !
  !>
  !> Implementation of "aggregate": time_average_dfire task
  !>
  SUBROUTINE aggregate_time_average_dist_fire(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_time_average_dist_fire'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(DFIRE_)

    ! DFIRE_ 2D
!    dsl4jsb_Aggregate_onChunk(DFIRE_, precip_qc_mavg                    , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_time_average_dist_fire

#endif
END MODULE mo_dfire_interface
