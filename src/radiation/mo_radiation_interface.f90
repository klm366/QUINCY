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

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_hsm_class,          ONLY: t_hsm
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile
!!$  USE mo_jsb_memory_class,   ONLY: t_jsb_memory_p, t_jsb_memory
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_process_class,  ONLY: t_jsb_process_task
!  USE mo_jsb_parallel_iface, ONLY: get_my_global_mpi_id
!  USE mo_util,               ONLY: one_of, int2string

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
  dsl4jsb_Def_memory(VEG)
  dsl4jsb_Def_memory(ATM2LAND)

  PRIVATE

  PUBLIC :: radiation_surface

  TYPE, EXTENDS(t_jsb_process_task) :: radiation_surface
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_radiation_surface
    PROCEDURE, NOPASS :: Aggregate => aggregate_radiation_surface
  END TYPE radiation_surface

  CHARACTER(len=*), PARAMETER :: modname = 'mo_radiation_interface'

CONTAINS

  SUBROUTINE update_radiation_surface(hsm, iblk, ics, ice, nc)

    !-----------------------------------------------------------------------
    !  DECLARATIONS
    USE mo_radiation_process, ONLY: calc_radiation_surface_net

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_hsm),              INTENT(inout) :: hsm
    INTEGER,           INTENT(in) :: iblk, ics, ice
    INTEGER, OPTIONAL, INTENT(in) :: nc

    !-----------------------------------------------------------------------
    ! POINTERS TO VARIABLES IN MEMORY
    TYPE(t_jsb_model), POINTER    :: model
    CLASS(t_jsb_tile), POINTER :: tile

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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_radiation_surface'

    !-----------------------------------------------------------------------
    !
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of surface radiation')

    SELECT TYPE(t=>hsm%current)
    CLASS IS(t_jsb_tile)
      tile => t
      model => Get_model(tile%owner_model_id)
    END SELECT

    ! Get reference to variables for current block
    !
!!$    CALL Get_process_memory(mem_p, ATM2LAND, atm2land__mem)
!!$    CALL Get_process_memory(mem_p, RAD, rad__mem)
!!$    CALL Get_process_memory(mem_p, SRF, srf__mem)
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


  END SUBROUTINE update_radiation_surface

  SUBROUTINE aggregate_radiation_surface(hsm, iblk, ics, ice, nc)

    CLASS(t_hsm),              INTENT(inout) :: hsm
    INTEGER,             INTENT(in)    :: iblk, ics, ice
    INTEGER, OPTIONAL,   INTENT(in)    :: nc

    TYPE(t_jsb_model), POINTER    :: model
    CLASS(t_jsb_tile), POINTER :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_radiation_surface'

    SELECT TYPE(t=>hsm%current)
      CLASS IS(t_jsb_tile)
      tile => t
      model => Get_model(tile%owner_model_id)
    END SELECT

  END SUBROUTINE aggregate_radiation_surface


#endif
END MODULE mo_radiation_interface
