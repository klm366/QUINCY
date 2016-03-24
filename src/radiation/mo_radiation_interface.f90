!>
!! @brief Contains the interfaces to the radiation processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-14)
!!

MODULE mo_radiation_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module
  
  ! Use of basic structures
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish

  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_lct_class,      ONLY: LAND_TYPE, LAKE_TYPE, GLACIER_TYPE, BARE_TYPE, PFT_TYPE
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes SRF, SOIL, VEG, ATM2LAND, RAD, LAKE

  ! Use process configurations
  dsl4jsb_Use_config(SRF)
  dsl4jsb_Use_config(SOIL)
  dsl4jsb_Use_config(RAD)

  ! Use of process memories
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SRF)
  dsl4jsb_Use_memory(VEG)
  dsl4jsb_Use_memory(RAD)
  dsl4jsb_Use_memory(LAKE)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_radiation_process

  !> Type definition for radiation process
  TYPE, EXTENDS(t_jsb_process) :: t_radiation_process
  CONTAINS
    PROCEDURE :: Register_tasks => Register_radiation_tasks
  END TYPE t_radiation_process

  !> Type definition for radiation task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_radiation_surface
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_radiation_surface
    PROCEDURE, NOPASS :: Aggregate => aggregate_radiation_surface
  END TYPE tsk_radiation_surface

  INTERFACE tsk_radiation_surface
    PROCEDURE Create_task_radiation_surface
  END INTERFACE tsk_radiation_surface

  CHARACTER(len=*), PARAMETER :: modname = 'mo_radiation_interface'

CONTAINS
  ! ================================================================================================================================
  !! Constructors for tasks
  
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for radiation task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "radiation"
  !!
  FUNCTION Create_task_radiation_surface(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_radiation_surface::return_ptr)

    return_ptr%name = 'radiation_surface'
    return_ptr%process_id = RAD
    return_ptr%owner_model_id = model_id
! !!$    return_ptr%Integrate => update_radiation_surface
! !!$    return_ptr%Aggregate => aggregate_radiation_surface

  END FUNCTION Create_task_radiation_surface

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for radiation process
  !!
  !! @param[in,out] this      Instance of <PROCESS_NAME_LOWER_CASE> process class  
  !! @param[in]     model_id  Model id  
  !!  
  SUBROUTINE Register_radiation_tasks(this, model_id)

    CLASS(t_Radiation_process), INTENT(inout) :: this
    INTEGER, INTENT(in) :: model_id

!!$    CALL this%Register_task(tsk_Radiation_surface(name="radiation_surface", process_id=RAD))
    CALL this%Register_task(tsk_Radiation_surface(model_id))

  END SUBROUTINE Register_radiation_tasks


#endif
END MODULE mo_radiation_interface
