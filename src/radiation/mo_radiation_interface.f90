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

  USE mo_jsb_control,        ONLY: debug_on
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish

  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes SRF, SOIL, VEG, ATM2LAND, RAD
  dsl4jsb_Use_config(SRF)
  dsl4jsb_Use_config(SOIL)
  dsl4jsb_Use_config(RAD)

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SRF)
  dsl4jsb_Use_memory(VEG)
  dsl4jsb_Use_memory(RAD)

  IMPLICIT NONE
  dsl4jsb_Def_config(RAD)
  dsl4jsb_Def_memory(RAD)
  !CLASS(t_rad_memory), POINTER :: rad__mem
  dsl4jsb_Def_memory(VEG)
  dsl4jsb_Def_memory(ATM2LAND)

  PRIVATE

  PUBLIC :: t_Radiation_process

  TYPE, EXTENDS(t_jsb_process) :: t_Radiation_process
  CONTAINS
    PROCEDURE :: Register_tasks => Register_radiation_tasks
  END TYPE t_Radiation_process

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_Radiation_surface
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_radiation_surface
    PROCEDURE, NOPASS :: Aggregate => aggregate_radiation_surface
  END TYPE tsk_Radiation_surface

  INTERFACE tsk_Radiation_surface
    PROCEDURE Create_task_radiation_surface
  END INTERFACE tsk_Radiation_surface

  CHARACTER(len=*), PARAMETER :: modname = 'mo_radiation_interface'

CONTAINS

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

  SUBROUTINE Register_radiation_tasks(this, model_id)

    CLASS(t_Radiation_process), INTENT(inout) :: this
    INTEGER, INTENT(in) :: model_id

!!$    CALL this%Register_task(tsk_Radiation_surface(name="radiation_surface", process_id=RAD))
    CALL this%Register_task(tsk_Radiation_surface(model_id))

  END SUBROUTINE Register_radiation_tasks

!!$  SUBROUTINE update_radiation_surface(model_id, tile_id, iblk, ics, ice, nc)
  SUBROUTINE update_radiation_surface(tile, config, options)

    !-----------------------------------------------------------------------
    !  DECLARATIONS
    USE mo_radiation_process, ONLY: calc_radiation_surface_net

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p), INTENT(in)    :: config(:)
    TYPE(t_jsb_task_options), INTENT(in) :: options

    !-----------------------------------------------------------------------
    ! POINTERS TO VARIABLES IN MEMORY

    dsl4jsb_Def_memory(ATM2LAND)
    dsl4jsb_Def_memory(SRF)
    !dsl4jsb_Def_memory(RAD) ! siehe oben im Modul

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: swvis_srf_net
    dsl4jsb_Real2D_onChunk :: swnir_srf_net
    dsl4jsb_Real2D_onChunk :: sw_srf_net
    dsl4jsb_Real2D_onChunk :: lw_srf_net
    dsl4jsb_Real2D_onChunk :: rad_srf_net
    dsl4jsb_Real2D_onChunk :: lw_srf_down
    dsl4jsb_Real2D_onChunk :: swvis_srf_down
    dsl4jsb_Real2D_onChunk :: swnir_srf_down
    dsl4jsb_Real2D_onChunk :: alb_vis_srf
    dsl4jsb_Real2D_onChunk :: alb_nir_srf
    dsl4jsb_Real2D_onChunk :: t_srf

    INTEGER :: iblk, ics, ice, nc

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_radiation_surface'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    !
    !-----------------------------------------------------------------------
    !
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of surface radiation on tile '//TRIM(tile%name))

    dsl4jsb_Get_memory(ATM2LAND)
    dsl4jsb_Get_memory(RAD)
    dsl4jsb_Get_memory(SRF)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swvis_srf_down, IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swnir_srf_down, IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, lw_srf_down,    IN)

    dsl4jsb_Get_var2D_onChunk(RAD, swvis_srf_net, OUT)
    dsl4jsb_Get_var2D_onChunk(RAD, swnir_srf_net, OUT)
    dsl4jsb_Get_var2D_onChunk(RAD, sw_srf_net,    OUT)
    dsl4jsb_Get_var2D_onChunk(RAD, lw_srf_net,    OUT)
    dsl4jsb_Get_var2D_onChunk(RAD, rad_srf_net,   OUT)
    dsl4jsb_Get_var2D_onChunk(RAD, alb_vis_srf,   IN)
    dsl4jsb_Get_var2D_onChunk(RAD, alb_nir_srf,   IN)
    dsl4jsb_Get_var2D_onChunk(SRF,       t_srf,         IN)

    CALL calc_radiation_surface_net(swvis_srf_down(:),    & ! input
                                    swnir_srf_down(:),    & ! input
                                    alb_vis_srf(:),       & ! input
                                    alb_nir_srf(:),       & ! input
                                    lw_srf_down(:),       & ! input
                                    t_srf(:),             & ! input
                                    swvis_srf_net(:),     & ! output
                                    swnir_srf_net(:),     & ! output
                                    sw_srf_net(:),        & ! output
                                    lw_srf_net(:),        & ! output
                                    rad_srf_net(:) )        ! output

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished update of surface radiation')

  END SUBROUTINE update_radiation_surface

  SUBROUTINE aggregate_radiation_surface(tile, config, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

!!$    TYPE(t_jsb_model), POINTER :: model
!!$    CLASS(t_jsb_tile), POINTER :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_radiation_surface'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

  END SUBROUTINE aggregate_radiation_surface


#endif
END MODULE mo_radiation_interface
