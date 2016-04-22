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

#include "dsl_jsbach.inc"

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
  USE mo_jsb_control,        ONLY: debug_on

  dsl4jsb_Use_processes SEB, RAD, HYDRO, ATM2LAND
  dsl4jsb_Use_config(SEB)

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(RAD)
  dsl4jsb_Use_memory(HYDRO)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_surface_energy_lake

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_lake'

CONTAINS

  SUBROUTINE update_surface_energy_lake(tile, options)

    USE mo_phy_schemes,            ONLY: qsat, q_effective, surface_dry_static_energy, heat_transfer_coef
    USE mo_jsb_lct_class,          ONLY: LAKE_TYPE, Contains_lct
    USE mo_jsb_physical_constants, ONLY: tmelt

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(SEB)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(RAD)
    dsl4jsb_Def_memory(HYDRO)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & t, &
      & t_wtr, &
      & t_ice, &
      & heat_cap, &
      & latent_hflx, &
      & sensible_hflx, &
      & forc_hflx, &
      & rad_srf_net, &
      & fract_lice, &
      & depth_lice

    ! Local variables
    REAL(wp) :: &
      hcap_wtr(options%nc), &
      hcap_ice(options%nc), &
      t_freeze(options%nc), &
      hflx_ice_melting(options%nc), &
      zero(options%nc)
    LOGICAL :: &
      l_ice_formation(options%nc), &
      l_ice_melting(options%nc)

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: dtime
    REAL(wp) :: &
      & ml_depth, &      !< Depth of lake mixed layer
      & min_ice_depth    !< Minimum ice thickness to start ice formation
    LOGICAL  :: lstart

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_energy_lake'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (.NOT. tile%Is_process_active(SEB) .OR. .NOT. Contains_lct(tile%lcts, LAKE_TYPE)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB)

    ml_depth      = dsl4jsb_Config(SEB)%lake_mixed_layer_depth
    min_ice_depth = dsl4jsb_Config(SEB)%lake_min_ice_depth

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SEB)
    dsl4jsb_Get_memory(RAD)
    dsl4jsb_Get_memory(HYDRO)

    dsl4jsb_Get_var2D_onChunk(SEB, t)         ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB, t_wtr)     ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB, t_ice)     ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB, heat_cap)     ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB, latent_hflx)         ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB, sensible_hflx)         ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB, forc_hflx)         ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO, fract_lice) ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO, depth_lice) ! OUT

    dsl4jsb_Get_var2D_onChunk(RAD, rad_srf_net)         ! OUT

    zero            (:) = 0._wp
    hflx_ice_melting(:) = 0._wp
    hcap_wtr        (:) = 0._wp
    hcap_ice        (:) = 0._wp
    l_ice_formation (:) = .FALSE.
    l_ice_melting   (:) = .FALSE.

    ! Calculate new water temperature and residual heat flux for ice formation over open water!
    CALL calc_lake_wtr_temp( &
      & t_wtr(:),                                       & !< inout; Surface temperature over open water
      & rad_srf_net(:)+latent_hflx(:)+sensible_hflx(:), & !< in;    Sum of all radiative and turbulent surface heat fluxes
      & forc_hflx(:),                                   & !< inout; Residual heat flux from ice formation (realized or unrealized)
      & hcap_wtr(:),                                    & !< inout; Heat capacity of lake over open water
      & l_ice_formation(:),                             & !< inout; Indicator for ice formation
      & l_ice_melting(:),                               & !< in;    Indicator for ice melting = .FALSE.
      & fract_lice(:) > 0.5_wp,                         & !< in;    Indicator for ice cover
      & dtime,                                          & !< in;    Time step
      & ml_depth,                                       & !< in;    Mixed-layer depth
      & min_ice_depth)                                    !< in;    Minimum ice thickness

    IF (ANY(hcap_wtr(:) < EPSILON(1._wp))) THEN
      CALL finish(TRIM(routine),'hcap_wtr can"t be zero')
    END IF

    t_ice(:) = tmelt

    ! Calculate new ice depth and residual heat flux from melting over ice

    
    ! Calculate new water temperature for the case of complete ice melting
    CALL calc_lake_wtr_temp( &
      & t_wtr(:),                                       & !< inout; Surface temperature over open water
      & zero(:),                                        & !< in;    Sum of radiative and turbulent surface heat fluxes (not used)
      & hflx_ice_melting(:),                            & !< inout; Residual heat flux from unrealized ice melting (not changed)
      & hcap_wtr(:),                                    & !< inout; Heat capacity of lake over open water (not changed)
      & l_ice_formation(:),                             & !< inout; Indicator for ice formation (not changed)
      & l_ice_melting(:),                               & !< in;    Indicator for ice melting = .TRUE.
      & fract_lice(:) > 0.5_wp,                         & !< in;    Indicator for ice cover
      & dtime,                                          & !< in;    Time step
      & ml_depth,                                       & !< in;    Mixed-layer depth
      & min_ice_depth)                                    !< in;    Minimum ice thickness

    t       (:) = (1._wp - fract_lice(:)) * t_wtr   (:) + fract_lice(:) * t_ice   (:)
    heat_cap(:) = (1._wp - fract_lice(:)) * hcap_wtr(:) + fract_lice(:) * hcap_ice(:)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_energy_lake

  ! Calculate new surface water temperature
  ELEMENTAL PURE SUBROUTINE calc_lake_wtr_temp(t_wtr, hflx, hflx_resid, hcap, l_ice_formation, l_ice_melting, &
    &                                          l_ice, dtime, ml_depth, min_ice_depth)

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, clw, rhoi, alf

    REAL(wp), INTENT(inout) :: t_wtr
    REAL(wp), INTENT(in)    :: hflx
    REAL(wp), INTENT(inout) :: hflx_resid
    REAL(wp), INTENT(inout) :: hcap
    LOGICAL,  INTENT(inout) :: l_ice_formation
    LOGICAL,  INTENT(in)    :: l_ice_melting
    LOGICAL,  INTENT(in)    :: l_ice
    REAL(wp), INTENT(in)    :: dtime
    REAL(wp), INTENT(in)    :: ml_depth
    REAL(wp), INTENT(in)    :: min_ice_depth

    REAL(wp) :: &
      & t_star, &
      & t_freeze

    IF (l_ice) RETURN

    IF (l_ice_melting) THEN
      t_wtr = t_wtr + dtime * hflx_resid / hcap
    ELSE
      ! Heat capacity of lake over open water
      ! Note: This is currently constant but in the future, different mixed-layer depths might be used
      hcap = rhoh2o * clw * ml_depth

      ! Water temperature below which a new slab of ice is formed
      t_freeze = tmelt - min_ice_depth * rhoi * alf / hcap

      ! Preliminary water temperature
      t_star = t_wtr + dtime * (hflx + hflx_resid) / hcap

      IF (t_star >= tmelt) THEN              ! No ice formation
        t_wtr = t_star
        hflx_resid = 0._wp
        l_ice_formation = .FALSE.
      ELSE                                   ! Realized or unrealized ice formation
        t_wtr = tmelt
        hflx_resid = (t_star - tmelt) * hcap / dtime      ! < 0
        l_ice_formation = t_star <= t_freeze
        ! l_ice_formation = .TRUE.  : hflx_resid will be used to form a new slab of ice
        !                 = .FALSE. : cooling is not large enough to form a new slab of ice and
        !                             hflx_resid will be used as forcing flux in the next time
        !                             step to compute water temperature t_wtr
      END IF

    END IF

  END SUBROUTINE calc_lake_wtr_temp

  ! Calculation of ice depth and fraction
  ELEMENTAL PURE SUBROUTINE calc_lake_ice(t_ice, hflx_resid, l_ice_formation, ice_depth, ice_fract, snow_depth, t_wtr, &
    & dtime, min_ice_depth)

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoi, rhos, alf, ki, ks

    REAL(wp), INTENT(in)    :: t_ice
    REAL(wp), INTENT(inout) :: hflx_resid       !< heat flux for ice formation (<0; in) or from ice melting (>0; out)
    LOGICAL,  INTENT(inout) :: l_ice_formation
    REAL(wp), INTENT(inout) :: ice_depth
    REAL(wp), INTENT(inout) :: ice_fract
    REAL(wp), INTENT(in)    :: snow_depth
    REAL(wp), INTENT(inout) :: t_wtr
    REAL(wp), INTENT(in)    :: dtime
    REAL(wp), INTENT(in)    :: min_ice_depth

    REAL(wp) :: zhflx_cond  !< Conductive heat flux
    REAL(wp) :: zice_depth_eff !< Effective ice thickness which takes into account snow layer on ice

    IF (l_ice_formation) THEN                              ! Formation of a new slab of ice
      ice_depth = - hflx_resid * dtime / (rhoi * alf)
      ice_fract = 1._wp
      hflx_resid = 0._wp
    ELSE IF (ice_depth >= min_ice_depth) THEN              ! Ice already present
      zice_depth_eff = ice_depth + (ki / ks) * (rhoi / rhos) * snow_depth
      zhflx_cond = ki * (t_ice - tmelt) / zice_depth_eff
    END IF

  END SUBROUTINE calc_lake_ice
  !
  ! ================================================================================================================================
  !
!!$  SUBROUTINE update_asselin_lake(tile, config, options)
!!$
!!$    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
!!$    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)
!!$    TYPE(t_jsb_task_options),   INTENT(in)    :: options
!!$
!!$    ! Local variables
!!$    !
!!$    dsl4jsb_Def_memory(SEB)
!!$
!!$    ! Pointers to variables in memory
!!$    dsl4jsb_Real2D_onChunk :: t
!!$    dsl4jsb_Real2D_onChunk :: t_old
!!$    dsl4jsb_Real2D_onChunk :: t_unfilt
!!$    dsl4jsb_Real2D_onChunk :: t_unfilt_old
!!$
!!$    ! Locally allocated variables
!!$    !
!!$    INTEGER :: iblk, ics, ice
!!$    REAL(wp) :: eps
!!$
!!$    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_asselin_lake'
!!$
!!$    iblk = options%iblk
!!$    ics  = options%ics
!!$    ice  = options%ice
!!$
!!$    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
!!$
!!$    ! Get reference to variables for current block
!!$    !
!!$    dsl4jsb_Get_memory(SEB)
!!$
!!$    dsl4jsb_Get_var2D_onChunk(SEB,      t)              ! OUT
!!$    dsl4jsb_Get_var2D_onChunk(SEB,      t_unfilt)       ! IN
!!$
!!$    ! Do not apply Asselin time filter
!!$    t(:) = t_unfilt(:)
!!$
!!$    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')
!!$
!!$  END SUBROUTINE update_asselin_lake

#endif
END MODULE mo_seb_lake
