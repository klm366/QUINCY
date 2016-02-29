!>
!! @brief Contains the interfaces to the lake processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_lake_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task

  dsl4jsb_Use_processes SRF, SOIL, ATM2LAND
  dsl4jsb_Use_config(SRF)
  dsl4jsb_Use_config(SOIL)

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SRF)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lake_process

  TYPE, EXTENDS(t_jsb_process) :: t_lake_process
  CONTAINS
    PROCEDURE :: Register_tasks => Register_lake_tasks
  END TYPE t_lake_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lake_interface'

CONTAINS

  SUBROUTINE Register_lake_tasks(this, model_id)

    CLASS(t_lake_process), INTENT(inout) :: this
    INTEGER, INTENT(in) :: model_id

    !CALL this%Register_task(tsk_radiation_surface(name="radiation_surface", process_id=LAKE))

  END SUBROUTINE Register_lake_tasks

#endif
END MODULE mo_lake_interface
