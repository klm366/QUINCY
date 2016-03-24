  !>
  !! @brief Contains the interfaces for the surface energy balance process.
  !!
  !! @author
  !!  Reiner Schnur, MPI-M Hamburg
  !!
  !! @par Revision History
  !! First version in mo_srf_interface.f90 by Reiner Schnur (2013-04-04)
  !! Separated out the surface energy balance part by Reiner Schnur (2016-03-22)
  !!
MODULE mo_seb_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_lct_class,      ONLY: LAND_TYPE, LAKE_TYPE
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes SEB
  dsl4jsb_Use_config(SEB)

  dsl4jsb_Use_memory(SEB)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_seb_process
  !
  ! ================================================================================================================================
  !
  TYPE, EXTENDS(t_jsb_process) :: t_seb_process
  CONTAINS
    PROCEDURE :: Register_tasks => Register_seb_tasks
  END TYPE t_seb_process

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_energy
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_energy
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_energy
  END TYPE tsk_surface_energy

  INTERFACE tsk_surface_energy
    PROCEDURE Create_task_surface_energy
  END INTERFACE tsk_surface_energy

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_interface'

CONTAINS

  FUNCTION Create_task_surface_energy(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_energy::return_ptr)

    return_ptr%name = 'surface_energy'
    return_ptr%process_id = SEB
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_surface_energy

  SUBROUTINE Register_seb_tasks(this, model_id)

    CLASS(t_seb_process), INTENT(inout) :: this
    INTEGER, INTENT(in) :: model_id

    CALL this%Register_task(tsk_surface_energy(model_id))

  END SUBROUTINE Register_seb_tasks
  !
  ! ================================================================================================================================
  !
  SUBROUTINE update_surface_energy(tile, config, options)

    USE mo_seb_land, ONLY: update_surface_energy_land
    USE mo_seb_lake, ONLY: update_surface_energy_lake

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_energy'

    IF (.NOT. tile%Is_process_active(SEB)) RETURN

    SELECT CASE (tile%lct_type)
    CASE (LAKE_TYPE)
      CALL update_surface_energy_lake(tile, config, options)
    CASE (LAND_TYPE)
      CALL update_surface_energy_land(tile, config, options)
    CASE DEFAULT
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END SELECT

  END SUBROUTINE update_surface_energy

  SUBROUTINE aggregate_surface_energy(tile, config, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_surface_energy'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_energy
  !
  ! ================================================================================================================================
  !
#endif
END MODULE mo_seb_interface
