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

  dsl$jsb Use_processes SRF, SOIL, VEG, ATM2LAND
  dsl$jsb Use_config(srf)
  dsl$jsb Use_config(soil)

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl$jsb Use_memory(atm2land)
  dsl$jsb Use_memory(soil)
  dsl$jsb Use_memory(srf)
  dsl$jsb Use_memory(veg)

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
    dsl$jsb Def_memory(atm2land_mem, atm2land)
    dsl$jsb Def_memory(srf_mem, srf)

    ! Pointers to variables in memory
    dsl$jsb onChunk_Real2D :: swvis_srf_net
    dsl$jsb onChunk_Real2D :: swnir_srf_net
    dsl$jsb onChunk_Real2D :: sw_srf_net
    dsl$jsb onChunk_Real2D :: lw_srf_net
    dsl$jsb onChunk_Real2D :: rad_srf_net
    dsl$jsb onChunk_Real2D :: lw_srf_down
    dsl$jsb onChunk_Real2D :: swvis_srf_down
    dsl$jsb onChunk_Real2D :: swnir_srf_down
    dsl$jsb onChunk_Real2D :: alb_vis_srf
    dsl$jsb onChunk_Real2D :: alb_nir_srf
    dsl$jsb onChunk_Real2D :: t_srf

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_radiation'

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of surface radiation')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl$jsb Get_memory(srf_mem, SRF)
    dsl$jsb Get_memory(atm2land_mem, ATM2LAND)

    dsl$jsb Get_var2D_onChunk( atm2land_mem, swvis_srf_down, IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, swnir_srf_down, IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, lw_srf_down,    OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      swvis_srf_net,  OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      swnir_srf_net,  OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      sw_srf_net,     OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      lw_srf_net,     OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      rad_srf_net,    OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      alb_vis_srf,    OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      alb_nir_srf,    OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      t_srf,          OUT)

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

    dsl$jsb Def_config(soil_conf, soil)
    dsl$jsb Def_memory(srf_mem, srf)
    dsl$jsb Def_memory(soil_mem, soil)
    dsl$jsb Def_memory(veg_mem, veg)
    dsl$jsb Def_memory(atm2land_mem, atm2land)

    ! Pointers to variables in memory
    dsl$jsb onChunk_Real2D :: glacier
    dsl$jsb onChunk_Real2D :: lai
    dsl$jsb onChunk_Real2D :: fract_fpc_max
    dsl$jsb onChunk_Real2D :: sfract_can
    dsl$jsb onChunk_Real2D :: wsn_can
    dsl$jsb onChunk_Real2D :: sfract_soil
    dsl$jsb onChunk_Real2D :: wsn_soil
    dsl$jsb onChunk_Real2D :: wfract_can
    dsl$jsb onChunk_Real2D :: wsr_can
    dsl$jsb onChunk_Real2D :: wfract_soil
    dsl$jsb onChunk_Real2D :: wsr_soil
    dsl$jsb onChunk_Real2D :: sfract_srf
    dsl$jsb onChunk_Real2D :: wfract_srf
    dsl$jsb onChunk_Real2D :: oro_stddev
    dsl$jsb onChunk_Real2D :: drag_srf
    dsl$jsb onChunk_Real2D :: t_srf
    dsl$jsb onChunk_Real2D :: press_srf
    dsl$jsb onChunk_Real2D :: q_air

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

    dsl$jsb Get_config(soil_conf, SOIL)

    ! Get reference to variables for current block
    !
    dsl$jsb Get_memory(srf_mem, SRF)
    dsl$jsb Get_memory(soil_mem, SOIL)
    dsl$jsb Get_memory(veg_mem, VEG)
    dsl$jsb Get_memory(atm2land_mem, ATM2LAND)

    dsl$jsb Get_var2D_onChunk( atm2land_mem, drag_srf,      IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, q_air,         IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, press_srf,     IN)
    dsl$jsb Get_var2D_onChunk( soil_mem,     sfract_soil,   OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     wsn_soil,      OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     wfract_soil,   OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     wsr_soil,      OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      glacier,       OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      sfract_srf,    OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      wfract_srf,    OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      oro_stddev,    IN)
    dsl$jsb Get_var2D_onChunk( srf_mem,      t_srf,         OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,      lai,           OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,      fract_fpc_max, OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,      sfract_can,    OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,      wsn_can,       OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,      wfract_can,    OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,      wsr_can,       OUT)

    ! Maximum capacity of skin reservoir (soil + canopy)
    wsr_max(:) = soil_conf%wsr_capacity * (1._wp + lai(:) * fract_fpc_max(:))

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
    dsl$jsb Def_config(srf_conf, srf)
    dsl$jsb Def_memory(srf_mem, srf)
    dsl$jsb Def_memory(soil_mem, soil)
    dsl$jsb Def_memory(atm2land_mem, atm2land)

    LOGICAL                       :: lstart

    ! Pointers to variables in memory
    dsl$jsb onChunk_Real2D :: t_srf
    dsl$jsb onChunk_Real2D :: t_srf_old
    dsl$jsb onChunk_Real2D :: t_srf_unfilt
    dsl$jsb onChunk_Real2D :: t_srf_unfilt_old
    dsl$jsb onChunk_Real2D :: t_eff_srf
    dsl$jsb onChunk_Real2D :: qsat_srf
!!$   dsl$jsb onChunk_Real2D :: qsat_srf_old
    dsl$jsb onChunk_Real2D :: s_srf
!!$   dsl$jsb onChunk_Real2D :: s_srf_old
    dsl$jsb onChunk_Real2D :: press_srf
    dsl$jsb onChunk_Real2D :: q_air
    dsl$jsb onChunk_Real2D :: fact_qsat_srf
    dsl$jsb onChunk_Real2D :: fact_q_air
    dsl$jsb onChunk_Real2D :: drag_srf
    dsl$jsb onChunk_Real2D :: rad_srf_net
    dsl$jsb onChunk_Real2D :: grnd_hflx
    dsl$jsb onChunk_Real2D :: hcap_grnd
    dsl$jsb onChunk_Real2D :: wsn_srf
    dsl$jsb onChunk_Real2D :: sfract_srf
    dsl$jsb onChunk_Real2D :: t_acoef
    dsl$jsb onChunk_Real2D :: t_bcoef
    dsl$jsb onChunk_Real2D :: q_acoef
    dsl$jsb onChunk_Real2D :: q_bcoef
    dsl$jsb onChunk_Real2D :: glacier

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

    dsl$jsb Get_config(srf_conf, SRF)

    ! Get reference to variables for current block
    !
    dsl$jsb Get_memory(srf_mem, SRF)
    dsl$jsb Get_memory(soil_mem, SOIL)
    dsl$jsb Get_memory(atm2land_mem, ATM2LAND)

    dsl$jsb Get_var2D_onChunk( atm2land_mem, t_acoef,          IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, t_bcoef,          IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, q_acoef,          IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, q_bcoef,          IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, press_srf,        IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, q_air,            IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, drag_srf,         IN)
    dsl$jsb Get_var2D_onChunk( srf_mem,      t_srf,            OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      t_srf_old,        OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      t_srf_unfilt,     OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      t_srf_unfilt_old, OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      t_eff_srf,        OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      qsat_srf,         OUT)
!!$    dsl$jsb Get_var2D_onChunk( srf_mem,      qsat_srf_old,     OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      s_srf,            OUT)
!!$    dsl$jsb Get_var2D_onChunk( srf_mem,      s_srf_old,        OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      fact_q_air,       OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      fact_qsat_srf,    OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      rad_srf_net,      OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      wsn_srf,          OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      sfract_srf,       OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      glacier,          OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     grnd_hflx,        OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     hcap_grnd,        OUT)

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
    WHERE (wsn_srf(:) > srf_conf%wsn_critical .OR. glacier(:) > 0.5_wp)
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
    dsl$jsb Def_memory(srf_mem, srf)

    LOGICAL                       :: lstart
    REAL(wp)                      :: eps

    ! Pointers to variables in memory
    dsl$jsb onChunk_Real2D :: t_srf
    dsl$jsb onChunk_Real2D :: t_srf_old
    dsl$jsb onChunk_Real2D :: t_srf_unfilt
    dsl$jsb onChunk_Real2D :: t_srf_unfilt_old

    dsl$jsb Get_memory(srf_mem, SRF)

    dsl$jsb Get_var2D_onChunk( srf_mem,  t_srf,            OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,  t_srf_old,        OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,  t_srf_unfilt,     OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,  t_srf_unfilt_old, OUT)

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
    dsl$jsb Def_memory(srf_mem, srf)
    dsl$jsb Def_memory(atm2land_mem, ATM2LAND)

    ! Pointers to variables in memory
    dsl$jsb onChunk_Real2D :: qsat_srf
    dsl$jsb onChunk_Real2D :: s_srf
    dsl$jsb onChunk_Real2D :: t_acoef
    dsl$jsb onChunk_Real2D :: t_bcoef
    dsl$jsb onChunk_Real2D :: q_acoef
    dsl$jsb onChunk_Real2D :: q_bcoef
    dsl$jsb onChunk_Real2D :: fact_qsat_srf
    dsl$jsb onChunk_Real2D :: fact_qsat_trans_srf
    dsl$jsb onChunk_Real2D :: fact_q_air
    dsl$jsb onChunk_Real2D :: drag_srf
    dsl$jsb onChunk_Real2D :: latent_hflx
    dsl$jsb onChunk_Real2D :: sensible_hflx
    dsl$jsb onChunk_Real2D :: evapopot
    dsl$jsb onChunk_Real2D :: evapotrans
    dsl$jsb onChunk_Real2D :: trans
    dsl$jsb onChunk_Real2D :: sfract_srf

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
    dsl$jsb Get_memory(srf_mem, SRF)
    dsl$jsb Get_memory(atm2land_mem, ATM2LAND)

    dsl$jsb Get_var2D_onChunk( atm2land_mem, t_acoef,             IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, t_bcoef,             IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, q_acoef,             IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, q_bcoef,             IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, drag_srf,            IN)
    dsl$jsb Get_var2D_onChunk( srf_mem,      qsat_srf,            OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      s_srf,               OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      fact_q_air,          OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      fact_qsat_srf,       OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      fact_qsat_trans_srf, OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      sensible_hflx,       OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      latent_hflx,         OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      evapotrans,          OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      trans,               OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      evapopot,            OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      sfract_srf,          OUT)

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
    dsl$jsb Def_memory(srf_mem, srf)
    dsl$jsb Def_memory(soil_mem, soil)
    dsl$jsb Def_memory(veg_mem, veg)
    dsl$jsb Def_memory(atm2land_mem, atm2land)

    ! Pointers to variables in memory
    dsl$jsb onChunk_Real2D :: fract_fpc
    dsl$jsb onChunk_Real3D :: ws_l
    dsl$jsb onChunk_Real3D :: ws_fc_l
    dsl$jsb onChunk_Real3D :: ws_pwp_l
!!$   dsl$jsb onChunk_Real2D :: sfract_can
!!$   dsl$jsb onChunk_Real2D :: sfract_soil
!!$   dsl$jsb onChunk_Real2D :: wfract_can
!!$   dsl$jsb onChunk_Real2D :: wfract_soil
    dsl$jsb onChunk_Real2D :: sfract_srf
    dsl$jsb onChunk_Real2D :: wfract_srf
    dsl$jsb onChunk_Real2D :: fact_q_air
    dsl$jsb onChunk_Real2D :: fact_qsat_srf
    dsl$jsb onChunk_Real2D :: fact_qsat_trans_srf
    dsl$jsb onChunk_Real2D :: qsat_srf
    dsl$jsb onChunk_Real2D :: canopy_cond
    dsl$jsb onChunk_Real2D :: pch
    dsl$jsb onChunk_Real2D :: wind_air
    dsl$jsb onChunk_Real2D :: q_air

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
    dsl$jsb Get_memory(srf_mem, SRF)
    dsl$jsb Get_memory(soil_mem, SOIL)
    dsl$jsb Get_memory(veg_mem, VEG)
    dsl$jsb Get_memory(atm2land_mem, ATM2LAND)

    dsl$jsb Get_var2D_onChunk( atm2land_mem, pch,                 IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, wind_air,            IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, q_air,               IN)
    dsl$jsb Get_var2D_onChunk( srf_mem,      qsat_srf,            OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      sfract_srf,          OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      wfract_srf,          OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      fact_q_air,          OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      fact_qsat_srf,       OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      fact_qsat_trans_srf, OUT)
    dsl$jsb Get_var3D_onChunk( soil_mem,     ws_l,                OUT)
    dsl$jsb Get_var3D_onChunk( soil_mem,     ws_fc_l,             OUT)
    dsl$jsb Get_var3D_onChunk( soil_mem,     ws_pwp_l,            OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,      fract_fpc,           OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,      canopy_cond,         OUT)
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

    dsl$jsb Def_memory(srf_mem, srf)
    dsl$jsb Def_memory(soil_mem, soil)

    ! Pointers to variables in memory
    dsl$jsb onChunk_Real2D :: rough_m_srf
    dsl$jsb onChunk_Real2D :: rough_h_srf
    dsl$jsb onChunk_Real2D :: sfract_soil

    dsl$jsb Get_memory(srf_mem, SRF)
    dsl$jsb Get_memory(soil_mem, SOIL)

    dsl$jsb Get_var2D_onChunk( srf_mem,  rough_m_srf, OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,  rough_h_srf, OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem, sfract_soil, OUT)

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
