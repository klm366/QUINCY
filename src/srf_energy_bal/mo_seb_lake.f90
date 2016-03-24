!>
!! @brief Contains the routines for the surface energy balance on LAKE lct_type.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2016-03-22)
!!
MODULE mo_seb_lake
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options

  dsl4jsb_Use_processes SEB, RAD, ATM2LAND
  dsl4jsb_Use_config(SEB)

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(RAD)

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_lake'

CONTAINS

  SUBROUTINE update_surface_energy_lake(tile, config, options)

    USE mo_phy_schemes,            ONLY: qsat, q_effective, surface_dry_static_energy, heat_transfer_coef
    USE mo_srf_process,            ONLY: surface_temp_implicit
    USE mo_echam_vdiff_params,     ONLY: tpfac2, tpfac3
    USE mo_jsb_physical_constants, ONLY: tmelt

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(LAKE)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t

    ! Local variables
    INTEGER :: iblk, ics, ice, nc
    LOGICAL                       :: lstart

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_energy_lake'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_active(SRF)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(LAKE)

    dsl4jsb_Get_var2D_onChunk(SEB, t) ! IN

    t(:) = dsl4jsb_var2D_onChunk(LAKE,t_srf)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_energy_lake

#endif
END MODULE mo_seb_lake
