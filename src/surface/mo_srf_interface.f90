!>
!! @brief Contains the interfaces to the srf processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_srf_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_class,       ONLY: get_model
  USE mo_jsb_tile_class,  ONLY: t_jsb_tile

  dsl4jsb_Use_processes SRF, SOIL, VEG, ATM2LAND
  dsl4jsb_Use_config(SRF)
  dsl4jsb_Use_config(SOIL)

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SRF)
  dsl4jsb_Use_memory(VEG)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_surface_radiation, update_snow_and_skin_fraction, update_surface_energy
  PUBLIC :: update_surface_fluxes, update_t_srf_asselin, update_humidity_scaling, update_roughness

  CHARACTER(len=*), PARAMETER :: modname = 'mo_srf_interface'

CONTAINS

  SUBROUTINE update_surface_radiation(tile, iblk, ics, ice, nc)

    USE mo_phy_schemes, ONLY: lwnet_from_lwdown

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    dsl4jsb_Def_memory(ATM2LAND)
    dsl4jsb_Def_memory(SRF)

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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_radiation'

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of surface radiation')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(ATM2LAND)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swvis_srf_down, IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swnir_srf_down, IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, lw_srf_down,    IN)

    dsl4jsb_Get_var2D_onChunk(SRF, swvis_srf_net, OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, swnir_srf_net, OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, sw_srf_net,    OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, lw_srf_net,    OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, rad_srf_net,   OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, alb_vis_srf,   IN)
    dsl4jsb_Get_var2D_onChunk(SRF, alb_nir_srf,   IN)
    dsl4jsb_Get_var2D_onChunk(SRF, t_srf,         IN)

    ! Compute net SW radiation from downward SW and albedo
    swvis_srf_net(:) = swvis_srf_down(:) * (1._wp - alb_vis_srf(:))
    swnir_srf_net(:) = swnir_srf_down(:) * (1._wp - alb_nir_srf(:))
    sw_srf_net   (:) = swvis_srf_net(:) + swnir_srf_net(:)

    ! Compute LW net radiation from incoming and the thermal radiation from the surface
    lw_srf_net(:) = lwnet_from_lwdown(lw_srf_down(:), t_srf(:))

    ! Compute net radiation
    rad_srf_net(:) = sw_srf_net(:) + lw_srf_net(:)

  END SUBROUTINE update_surface_radiation

  SUBROUTINE update_snow_and_skin_fraction( tile, iblk, ics, ice, nc, dtime, steplen)

    USE mo_phy_schemes,    ONLY: heat_transfer_coef, qsat
    USE mo_soil_process,   ONLY: snow_fract
    USE mo_srf_constants,  ONLY: rhoh2o

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc
    REAL(wp),          INTENT(in) :: dtime, steplen

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model

    dsl4jsb_Def_config(SOIL)
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(ATM2LAND)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: glacier
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: sfract_can
    dsl4jsb_Real2D_onChunk :: wsn_can
    dsl4jsb_Real2D_onChunk :: sfract_soil
    dsl4jsb_Real2D_onChunk :: wsn_soil
    dsl4jsb_Real2D_onChunk :: wfract_can
    dsl4jsb_Real2D_onChunk :: wsr_can
    dsl4jsb_Real2D_onChunk :: wfract_soil
    dsl4jsb_Real2D_onChunk :: wsr_soil
    dsl4jsb_Real2D_onChunk :: sfract_srf
    dsl4jsb_Real2D_onChunk :: wfract_srf
    dsl4jsb_Real2D_onChunk :: oro_stddev
    dsl4jsb_Real2D_onChunk :: drag_srf
    dsl4jsb_Real2D_onChunk :: t_srf
    dsl4jsb_Real2D_onChunk :: press_srf
    dsl4jsb_Real2D_onChunk :: q_air

    ! Locally allocated vectors
    REAL(wp) :: &
      wsr_max(nc), &    !< Water holding capacity of skin reservoir (soil and canopy)
      evap_pot(nc), &
      qsat_srf_old(nc), &
      t_srf_old(nc)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_snow_and_skin_fraction'

    IF (.NOT. tile%Is_process_active(SOIL)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of snow and skin reservoir fraction')

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SOIL)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(ATM2LAND)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, drag_srf,      IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_air,         IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, press_srf,     IN)

    dsl4jsb_Get_var2D_onChunk(SOIL, sfract_soil,   OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, wsn_soil,      OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, wfract_soil,   OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, wsr_soil,      OUT)

    dsl4jsb_Get_var2D_onChunk(SRF, glacier,       IN)
    dsl4jsb_Get_var2D_onChunk(SRF, sfract_srf,    OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, wfract_srf,    OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, oro_stddev,    IN)
    dsl4jsb_Get_var2D_onChunk(SRF, t_srf,         IN)

    dsl4jsb_Get_var2D_onChunk(VEG, lai,           IN)
    dsl4jsb_Get_var2D_onChunk(VEG, fract_fpc_max, IN)
    dsl4jsb_Get_var2D_onChunk(VEG, sfract_can,    OUT)
    dsl4jsb_Get_var2D_onChunk(VEG, wsn_can,       OUT)
    dsl4jsb_Get_var2D_onChunk(VEG, wfract_can,    OUT)
    dsl4jsb_Get_var2D_onChunk(VEG, wsr_can,       OUT)

    ! Maximum capacity of skin reservoir (soil + canopy)
    wsr_max(:) = dsl4jsb_Config(SOIL)%wsr_capacity * (1._wp + lai(:) * fract_fpc_max(:))

    ! Wet skin fraction on soil
    wfract_soil(:) = MERGE(MIN(1._wp, (wsr_soil(:) + wsr_can(:)) / wsr_max(:)), 0._wp, wsr_max(:) > EPSILON(1._wp))

    ! Wet skin fraction on canopy
!!$    wfract_can(:) = MERGE(MIN(1._wp, wsr_can(:) / wsr_max(:)), 0._wp, wsr_max(:) > EPSILON(1._wp))
    wfract_can(:) = 0._wp

    ! Snow fraction on soil
    sfract_soil(:) = snow_fract(wsn_soil(:), oro_stddev(:))
    sfract_soil(:) = MERGE(sfract_can(:), sfract_soil(:), sfract_soil(:) < EPSILON(1._wp) .AND. sfract_can(:) > EPSILON(1._wp))

    ! Snow fraction on canopy
    sfract_can(:) = MERGE(MIN(1._wp, wsn_can(:) / wsr_max(:)), 0._wp, wsr_max(:) > EPSILON(1._wp))

    ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than
    ! equivalent snow water content from soil and canopy; same for skin reservoir
    !
    ! Potential evaporation using old values of air and surface humidity
    !
    t_srf_old(:) = t_srf(:)
    qsat_srf_old(:) = qsat(t_srf_old(:), press_srf(:))
    evap_pot(:) =  heat_transfer_coef(drag_srf(:), steplen) * (qsat_srf_old(:) - q_air(:)) ! Positive upwards
    !
    WHERE (sfract_soil(:) > 0._wp .AND. wsn_soil(:)+wsn_can(:) > EPSILON(1._wp))
      sfract_soil(:) = sfract_soil(:) / MAX(1._wp, sfract_soil(:) * evap_pot(:) * dtime / (rhoh2o * &
        & MAX(EPSILON(1._wp),(wsn_soil(:) + wsn_can(:)))))
    END WHERE
    WHERE (wfract_soil(:) > 0._wp .AND. wsr_soil(:)+wsr_can(:) > EPSILON(1._wp))
      wfract_soil(:) = wfract_soil(:) / MAX(1._wp, (1._wp - sfract_soil(:)) * evap_pot(:) * dtime / (rhoh2o * &
        & MAX(EPSILON(1._wp), (wsr_soil(:) + wsr_can(:)))))
    END WHERE

    wfract_srf(:) = wfract_soil(:) + wfract_can(:)
    sfract_srf(:) = sfract_soil(:) !+ sfract_can(:)

    WHERE (glacier(:) > 0.5_wp)
      wfract_srf(:)  = 0._wp
      sfract_srf(:)  = 1._wp
      wfract_soil(:) = 0._wp
      sfract_soil(:) = 1._wp
      wfract_can(:)  = 0._wp
      sfract_can(:)  = 0._wp
    END WHERE

  END SUBROUTINE update_snow_and_skin_fraction

  SUBROUTINE update_surface_energy( tile, iblk, ics, ice, nc, dtime, steplen)

    USE mo_phy_schemes,        ONLY: qsat, q_effective, surface_dry_static_energy, heat_transfer_coef
    USE mo_srf_process,        ONLY: surface_temp_implicit
!!$    USE mo_srf_constants,      ONLY: eps
    USE mo_echam_vdiff_params, ONLY: tpfac2, tpfac3
    USE mo_physical_constants, ONLY: tmelt

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc
    REAL(wp),          INTENT(in) :: dtime, steplen

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    dsl4jsb_Def_config(SRF)
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_memory(ATM2LAND)

    LOGICAL                       :: lstart

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_srf
    dsl4jsb_Real2D_onChunk :: t_srf_old
    dsl4jsb_Real2D_onChunk :: t_srf_unfilt
    dsl4jsb_Real2D_onChunk :: t_srf_unfilt_old
    dsl4jsb_Real2D_onChunk :: t_eff_srf
    dsl4jsb_Real2D_onChunk :: qsat_srf
!!$   dsl4jsb_Real2D_onChunk :: qsat_srf_old
    dsl4jsb_Real2D_onChunk :: s_srf
!!$   dsl4jsb_Real2D_onChunk :: s_srf_old
    dsl4jsb_Real2D_onChunk :: press_srf
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: fact_qsat_srf
    dsl4jsb_Real2D_onChunk :: fact_q_air
    dsl4jsb_Real2D_onChunk :: drag_srf
    dsl4jsb_Real2D_onChunk :: rad_srf_net
    dsl4jsb_Real2D_onChunk :: grnd_hflx
    dsl4jsb_Real2D_onChunk :: hcap_grnd
    dsl4jsb_Real2D_onChunk :: wsn_srf
    dsl4jsb_Real2D_onChunk :: sfract_srf
    dsl4jsb_Real2D_onChunk :: t_acoef
    dsl4jsb_Real2D_onChunk :: t_bcoef
    dsl4jsb_Real2D_onChunk :: q_acoef
    dsl4jsb_Real2D_onChunk :: q_bcoef
    dsl4jsb_Real2D_onChunk :: glacier

    ! Locally allocated vectors
    !
    REAL(wp) ::             &
      qsat_srf_old(nc),     &
      s_srf_old(nc),        &
      dQdT(nc),             & !< Sensitivity of saturated surface specific humidity to temperature
      t2s_conv(nc),         & !< Conversion factor from temperature to dry static energy (C_pd * (1+(delta-1)*q_v))
      heat_tcoef(nc)          !< Heat transfer coefficient (rho*C_h*|v|)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_energy'

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of surface temperature')

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SRF)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(ATM2LAND)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, t_acoef,   IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, t_bcoef,   IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_acoef,   IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_bcoef,   IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, press_srf, IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_air,     IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, drag_srf,  IN)

    dsl4jsb_Get_var2D_onChunk(SRF, t_srf,            OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, t_srf_old,        OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, t_srf_unfilt,     OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, t_srf_unfilt_old, OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, t_eff_srf,        OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, qsat_srf,         OUT)
!!$    dsl4jsb_Get_var2D_onChunk( SRF,      qsat_srf_old,     OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, s_srf,            OUT)
!!$    dsl4jsb_Get_var2D_onChunk( SRF,      s_srf_old,        OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, fact_q_air,       OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, fact_qsat_srf,    OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, rad_srf_net,      OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, wsn_srf,          OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, sfract_srf,       OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, glacier,          OUT)

    dsl4jsb_Get_var2D_onChunk( SOIL, grnd_hflx, OUT)
    dsl4jsb_Get_var2D_onChunk( SOIL, hcap_grnd, OUT)

    ! Save old surface temperature and saturation spec. humidity
!!$    t_srf_unfilt_old(:) = t_srf_unfilt(:)
!!$    t_srf_old   (:) = t_srf(:)
    ! Temporary fix for temperatures over ocean getting too warm or too cold, leading to crash in convect_tables
    t_srf_unfilt_old(:) = MAX(MIN(t_srf_unfilt(:), 380._wp), 100._wp)
    t_srf_old(:)        = MAX(MIN(t_srf(:),        380._wp), 100._wp)
    qsat_srf_old(:) = qsat(t_srf_old(:), press_srf(:))

    ! Calculate surface saturation specific humidity

    ! Update surface temperature
    !
    ! 1. Compute sensitivity of surface saturation specific humidity to temperature
    dQdT(:) = ( qsat(t_srf(:) + 0.001_wp, press_srf(:)) - qsat_srf_old(:) ) * 1000._wp

    ! 2. Old dry static energy
    s_srf_old(:) = surface_dry_static_energy( t_srf_old(:), &
                                            & q_effective(qsat_srf_old(:), q_air(:), fact_qsat_srf(:), fact_q_air(:)))

    t2s_conv(:) = s_srf_old(:) / t_srf_old(:)  !

    heat_tcoef(:) = heat_transfer_coef(drag_srf(:), steplen)

    ! Compute new surface temperature and saturation humidity
    IF (ALL(qsat_srf(:) < 0._wp)) THEN          ! Start of experiment
      lstart = .TRUE.
      s_srf(:) = s_srf_old(:)
    ELSE
      lstart = .FALSE.
      CALL surface_temp_implicit( nc, steplen, &
        t2s_conv(:), &
        t_acoef(:), t_bcoef(:), q_acoef(:), q_bcoef(:), &
        s_srf_old(:), qsat_srf_old(:), dQdT(:), &
        rad_srf_net(:), grnd_hflx(:), heat_tcoef(:), &
        fact_q_air(:), fact_qsat_srf(:), &
        sfract_srf(:), hcap_grnd(:), &
        & s_srf(:))
    END IF

    ! New unfiltered surface temperature (\tilde(X)^(t+1)) (see below Eq. 3.3.1.5 in ECHAM3 manual, alpha = 1/tpfac2)
    t_srf_unfilt(:) = ( tpfac2 * s_srf(:) + tpfac3 * s_srf_old(:) ) / t2s_conv(:)

    ! Correct new dry static energy for snowmelt
    WHERE (wsn_srf(:) > dsl4jsb_Config(SRF)%wsn_critical .OR. glacier(:) > 0.5_wp)
      s_srf(:) = MIN(s_srf(:) / t2s_conv(:), tmelt) * t2s_conv(:)
    END WHERE

    ! Compute new surface qsat (after correction of dry static energy for snow melt!)
    qsat_srf(:) = 0._wp
    IF (lstart) THEN
      qsat_srf(:) = qsat_srf_old(:)
    ELSE
      qsat_srf(:) = qsat_srf_old(:) + dQdT(:) * (s_srf(:) - s_srf_old(:)) / t2s_conv(:)
    END IF

    ! 'Effective temperature' used in radheat of the atmosphere model (see Eq. 6.3 in ECHAM5 manual)
    t_eff_srf(:) = (t_srf_old(:)**3 * (4._wp * s_srf(:)/t2s_conv(:) - 3._wp * t_srf_old(:))) ** 0.25_wp

  END SUBROUTINE update_surface_energy

  ! Asselin time filter for surface temperature
  SUBROUTINE update_t_srf_asselin( tile, iblk, ics, ice, nc)

    USE mo_jsb_time, ONLY: get_asselin_coef

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    dsl4jsb_Def_memory(SRF)

    LOGICAL                       :: lstart
    REAL(wp)                      :: eps

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_srf
    dsl4jsb_Real2D_onChunk :: t_srf_old
    dsl4jsb_Real2D_onChunk :: t_srf_unfilt
    dsl4jsb_Real2D_onChunk :: t_srf_unfilt_old

    dsl4jsb_Get_memory(SRF)

    dsl4jsb_Get_var2D_onChunk( SRF, t_srf,            OUT)
    dsl4jsb_Get_var2D_onChunk( SRF, t_srf_old,        OUT)
    dsl4jsb_Get_var2D_onChunk( SRF, t_srf_unfilt,     OUT)
    dsl4jsb_Get_var2D_onChunk( SRF, t_srf_unfilt_old, OUT)

    eps = get_asselin_coef()
    IF (eps > 0._wp) THEN
      t_srf(:) = t_srf_unfilt_old(:) + eps * (t_srf_old(:) - 2._wp * t_srf_unfilt_old(:) + t_srf_unfilt(:))
    ELSE
      t_srf(:) = t_srf_unfilt(:)
    END IF

  END SUBROUTINE update_t_srf_asselin

  SUBROUTINE update_surface_fluxes( tile, iblk, ics, ice, nc, steplen)

    USE mo_phy_schemes,        ONLY: q_effective, surface_dry_static_energy, heat_transfer_coef
    USE mo_physical_constants, ONLY: alv, als

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc
    REAL(wp),          INTENT(in) :: steplen

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(ATM2LAND)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: qsat_srf
    dsl4jsb_Real2D_onChunk :: s_srf
    dsl4jsb_Real2D_onChunk :: t_acoef
    dsl4jsb_Real2D_onChunk :: t_bcoef
    dsl4jsb_Real2D_onChunk :: q_acoef
    dsl4jsb_Real2D_onChunk :: q_bcoef
    dsl4jsb_Real2D_onChunk :: fact_qsat_srf
    dsl4jsb_Real2D_onChunk :: fact_qsat_trans_srf
    dsl4jsb_Real2D_onChunk :: fact_q_air
    dsl4jsb_Real2D_onChunk :: drag_srf
    dsl4jsb_Real2D_onChunk :: latent_hflx
    dsl4jsb_Real2D_onChunk :: sensible_hflx
    dsl4jsb_Real2D_onChunk :: evapopot
    dsl4jsb_Real2D_onChunk :: evapotrans
    dsl4jsb_Real2D_onChunk :: trans
    dsl4jsb_Real2D_onChunk :: sfract_srf

    ! Locally allocated vectors
    !
    REAL(wp) ::                &
      & s_air          (nc),   &  !< Dry static energy at lowest atmospheric level
      & q_air          (nc),   &  !< Humidity at lowest atmospheric level
      & q_air_eff      (nc),   &  !< Effective humidity at lowest atmospheric level
      & qsat_srf_eff   (nc),   &  !< Effective surface saturation specific humidity
      heat_tcoef       (nc)       !< Heat transfer coefficient (rho*C_h*|v|)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_fluxes'

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of surface fluxes')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(ATM2LAND)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, t_acoef,  IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, t_bcoef,  IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_acoef,  IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_bcoef,  IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, drag_srf, IN)

    dsl4jsb_Get_var2D_onChunk(SRF, qsat_srf,            OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, s_srf,               OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, fact_q_air,          OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, fact_qsat_srf,       OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, fact_qsat_trans_srf, OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, sensible_hflx,       OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, latent_hflx,         OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, evapotrans,          OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, trans,               OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, evapopot,            OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, sfract_srf,          OUT)

    ! Compute new dry static energy and moisture at lowest atmospheric level by back-substitution
    !
    s_air(:) = t_acoef(:) * s_srf   (:) + t_bcoef(:)
    q_air(:) = q_acoef(:) * qsat_srf(:) + q_bcoef(:)

    ! Transfer coefficient
    heat_tcoef(:) = heat_transfer_coef(drag_srf(:), steplen)

    ! Compute effective air moisture and surface saturation humidity
    q_air_eff(:) = q_effective( SPREAD(0._wp, NCOPIES=nc, DIM=1), q_air(:), &
                              & SPREAD(1._wp, NCOPIES=nc, DIM=1), SPREAD(0._wp, NCOPIES=nc, DIM=1))
    qsat_srf_eff(:) = q_effective(qsat_srf(:), q_air(:), fact_qsat_srf(:), fact_q_air(:))

    ! Compute moisture fluxes
    !
    evapotrans(:) = heat_tcoef(:) * (q_air_eff(:) - qsat_srf_eff(:))        ! Evapotranspiratin
    evapopot  (:) = heat_tcoef(:) * (q_air    (:) - qsat_srf    (:))        ! Potential evaporation
    trans     (:) = fact_qsat_trans_srf(:) * evapopot(:)                    ! Transpiration

    ! Compute sensible heat flux
    !
    sensible_hflx(:) = heat_tcoef(:) * (s_air(:) - s_srf(:))

    ! Compute latent heat flux
    !
    latent_hflx(:) = alv * evapotrans(:) + (als - alv) * sfract_srf(:) * evapopot(:)

  END SUBROUTINE update_surface_fluxes

  SUBROUTINE update_humidity_scaling( tile, iblk, ics, ice, nc)

    USE mo_srf_constants,  ONLY: rhoh2o
    USE mo_soil_constants, ONLY: ws_wilt_fract
    USE mo_soil_process,   ONLY: relative_humidity

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(ATM2LAND)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: fract_fpc
    dsl4jsb_Real3D_onChunk :: ws_l
    dsl4jsb_Real3D_onChunk :: ws_fc_l
    dsl4jsb_Real3D_onChunk :: ws_pwp_l
!!$   dsl4jsb_Real2D_onChunk :: sfract_can
!!$   dsl4jsb_Real2D_onChunk :: sfract_soil
!!$   dsl4jsb_Real2D_onChunk :: wfract_can
!!$   dsl4jsb_Real2D_onChunk :: wfract_soil
    dsl4jsb_Real2D_onChunk :: sfract_srf
    dsl4jsb_Real2D_onChunk :: wfract_srf
    dsl4jsb_Real2D_onChunk :: fact_q_air
    dsl4jsb_Real2D_onChunk :: fact_qsat_srf
    dsl4jsb_Real2D_onChunk :: fact_qsat_trans_srf
    dsl4jsb_Real2D_onChunk :: qsat_srf
    dsl4jsb_Real2D_onChunk :: canopy_cond
    dsl4jsb_Real2D_onChunk :: pch
    dsl4jsb_Real2D_onChunk :: wind_air
    dsl4jsb_Real2D_onChunk :: q_air

    ! Locally allocated vectors
    REAL(wp) :: &
      fact_qsat_veg(nc),   &
      fact_qsat_trans(nc), &
      fact_qsat_soil(nc),  &
      fact_qair_veg(nc),   &
      fact_qair_soil(nc),  &
      rel_hum(nc)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_humidity_scaling'

    IF (.NOT. tile%Is_process_active(SOIL)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of humidity weights')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(ATM2LAND)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, pch,      IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, wind_air, IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_air,    IN)

    dsl4jsb_Get_var2D_onChunk(SRF, qsat_srf,            OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, sfract_srf,          OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, wfract_srf,          OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, fact_q_air,          OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, fact_qsat_srf,       OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, fact_qsat_trans_srf, OUT)

    dsl4jsb_Get_var3D_onChunk(SOIL, ws_l,     IN)
    dsl4jsb_Get_var3D_onChunk(SOIL, ws_fc_l,  IN)
    dsl4jsb_Get_var3D_onChunk(SOIL, ws_pwp_l, IN)

    dsl4jsb_Get_var2D_onChunk(VEG, fract_fpc,   IN)
    dsl4jsb_Get_var2D_onChunk(VEG, canopy_cond, IN)
!!$    sfract_can       => tile%mem%veg  %sfract_can%ptr          (ics:ice,    iblk)
!!$    sfract_soil      => tile%mem%soil %sfract_soil%ptr         (ics:ice,    iblk)
!!$    wfract_can       => tile%mem%veg  %wfract_can%ptr          (ics:ice,    iblk)
!!$    wfract_soil      => tile%mem%soil %wfract_soil%ptr         (ics:ice,    iblk)

    WHERE (ws_l(:,1) > ws_pwp_l(:,1))
      fact_qsat_veg(:) = sfract_srf(:) + (1._wp - sfract_srf(:)) *   &
        & (wfract_srf(:) + ( 1._wp - wfract_srf(:)) /                &
        & (1._wp + pch(:) * MAX(1._wp, wind_air(:)) / MAX(1.E-20_wp, canopy_cond(:))))
      fact_qsat_trans(:) = (1._wp - sfract_srf(:)) * &
        & (1._wp - wfract_srf(:)) /                  &
        & (1._wp + pch(:) * MAX(1._wp, wind_air(:)) / MAX(1.E-20_wp, canopy_cond(:)))
    ELSEWHERE
      fact_qsat_veg(:)   = sfract_srf(:) + (1._wp - sfract_srf(:)) * wfract_srf(:)
      fact_qsat_trans(:) = 0._wp
    END WHERE
    fact_qair_veg(:) = fact_qsat_veg(:)

    rel_hum(:) = relative_humidity(ws_l(:,1), ws_fc_l(:,1))

    WHERE (qsat_srf(:) * rel_hum(:) > q_air(:) .AND. rel_hum(:) > 1.e-10_wp)
      fact_qsat_soil(:) = sfract_srf(:) + (1._wp - sfract_srf(:)) * &
        & (wfract_srf(:) + (1._wp - wfract_srf(:)) * rel_hum(:))
      fact_qair_soil(:) = 1._wp
    ELSEWHERE
      fact_qsat_soil(:) = sfract_srf(:) + (1._wp - sfract_srf(:)) * wfract_srf(:)
      fact_qair_soil(:) = fact_qsat_soil(:)
    END WHERE
    WHERE (q_air(:) > qsat_srf(:))
      fact_qsat_soil(:) = 1._wp
      fact_qair_soil(:) = 1._wp
    END WHERE

    fact_qsat_srf(:) = fract_fpc(:) * fact_qsat_veg(:) + (1._wp - fract_fpc(:)) * fact_qsat_soil(:)
    fact_q_air(:)    = fract_fpc(:) * fact_qair_veg(:) + (1._wp - fract_fpc(:)) * fact_qair_soil(:)
    fact_qsat_trans_srf(:) = fract_fpc(:) * fact_qsat_trans(:)

  END SUBROUTINE update_humidity_scaling

  SUBROUTINE update_roughness( tile, iblk, ics, ice, nc)

    USE mo_srf_constants, ONLY: blending_height, roughness_snow

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc

    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: rough_m_srf
    dsl4jsb_Real2D_onChunk :: rough_h_srf
    dsl4jsb_Real2D_onChunk :: sfract_soil

    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)

    dsl4jsb_Get_var2D_onChunk( SRF,  rough_m_srf, OUT)
    dsl4jsb_Get_var2D_onChunk( SRF,  rough_h_srf, OUT)
    dsl4jsb_Get_var2D_onChunk( SOIL, sfract_soil, OUT)

    rough_h_srf(:) =            sfract_soil(:)  / (LOG(blending_height /             roughness_snow) ** 2)  &
                   & + (1._wp - sfract_soil(:)) / (LOG(blending_height / MIN(1._wp, rough_m_srf(:))) ** 2)

    ! Limit roughness length for heat by 1m (see below Eq.5.6 in ECHAM5 manual)
    rough_h_srf(:) = MIN(1._wp, blending_height * EXP(-1._wp / SQRT(rough_h_srf(:))))

  END SUBROUTINE update_roughness

!!$  SUBROUTINE update_srf( tile, iblk, ics, ice)
!!$
!!$    TYPE(t_jsb_tile),  POINTER    :: tile
!!$    INTEGER,           INTENT(in) :: iblk, ics, ice
!!$
!!$    INTEGER :: nc
!!$
!!$    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_srf'
!!$
!!$    IF (iblk == 1) CALL message(TRIM(routine), 'Starting update')
!!$
!!$    nc = ice - ics + 1
!!$
!!$    model => get_model(tile%owner_model_id)
!!$
!!$    ! Get reference to variables for current block
!!$    !
!!$    fract                => tile%mem%fract                               (ics:ice,    iblk)
!!$    IF (ALL(fract(:) < EPSILON(1._wp))) RETURN
!!$
!!$
!!$
!!$
!!$  END SUBROUTINE update_srf

#endif
END MODULE mo_srf_interface
