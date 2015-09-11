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

  Use_processes SRF, SOIL, VEG, ATM2LAND
  Use_config(srf)
  Use_config(soil)

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  Use_memory(atm2land)
  Use_memory(soil)
  Use_memory(srf)
  Use_memory(veg)

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
    Def_memory(atm2land_mem, atm2land)
    Def_memory(srf_mem, srf)

    ! Pointers to variables in memory
    REAL(wp), POINTER ::      &
      & swvis_srf_net    (:), &
      & swnir_srf_net    (:), &
      & sw_srf_net       (:), &
      & lw_srf_net       (:), &
      & rad_srf_net      (:), &
      & lw_srf_down      (:), &
      & swvis_srf_down   (:), &
      & swnir_srf_down   (:), &
      & alb_vis_srf      (:), &
      & alb_nir_srf      (:), &
      & t_srf            (:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_radiation'

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of surface radiation')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    Get_memory(srf_mem, SRF)
    Get_memory(atm2land_mem, ATM2LAND)

    swvis_srf_down => atm2land_mem     %swvis_srf_down%ptr (ics:ice, iblk)
    swnir_srf_down => atm2land_mem     %swnir_srf_down%ptr (ics:ice, iblk)
    lw_srf_down    => atm2land_mem     %lw_srf_down%ptr    (ics:ice, iblk)
    swvis_srf_net  => srf_mem %swvis_srf_net%ptr  (ics:ice, iblk)
    swnir_srf_net  => srf_mem %swnir_srf_net%ptr  (ics:ice, iblk)
    sw_srf_net     => srf_mem %sw_srf_net%ptr     (ics:ice, iblk)
    lw_srf_net     => srf_mem %lw_srf_net%ptr     (ics:ice, iblk)
    rad_srf_net    => srf_mem %rad_srf_net%ptr    (ics:ice, iblk)
    alb_vis_srf    => srf_mem %alb_vis_srf%ptr    (ics:ice, iblk)
    alb_nir_srf    => srf_mem %alb_nir_srf%ptr    (ics:ice, iblk)
    t_srf          => srf_mem %t_srf%ptr          (ics:ice, iblk)

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

    Def_config(soil_conf, soil)
    Def_memory(srf_mem, srf)
    Def_memory(soil_mem, soil)
    Def_memory(veg_mem, veg)
    Def_memory(atm2land_mem, atm2land)

    ! Pointers to variables in memory
    REAL(wp), POINTER :: &
                                !wsn           (:), &
      glacier       (:), &
      lai           (:), &
      fract_fpc_max (:), &
      sfract_can    (:), &
      wsn_can       (:), &
      sfract_soil   (:), &
      wsn_soil      (:), &
      wfract_can    (:), &
      wsr_can       (:), &
      wfract_soil   (:), &
      wsr_soil      (:), &
      sfract_srf    (:), &
      wfract_srf    (:), &
      oro_stddev    (:), &
      drag_srf      (:), &
      t_srf         (:), &
      press_srf     (:), &
      q_air         (:)

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

    Get_config(soil_conf, SOIL)

    ! Get reference to variables for current block
    !
    Get_memory(srf_mem, SRF)
    Get_memory(soil_mem, SOIL)
    Get_memory(veg_mem, VEG)
    Get_memory(atm2land_mem, ATM2LAND)

    glacier       => srf_mem  %glacier%ptr          (ics:ice, iblk)
    lai           => veg_mem  %lai%ptr              (ics:ice, iblk)
    fract_fpc_max => veg_mem  %fract_fpc_max%ptr    (ics:ice, iblk)
    sfract_can    => veg_mem  %sfract_can%ptr       (ics:ice, iblk)
    wsn_can       => veg_mem  %wsn_can%ptr          (ics:ice, iblk)
    sfract_soil   => soil_mem %sfract_soil%ptr      (ics:ice, iblk)
    wsn_soil      => soil_mem %wsn_soil%ptr         (ics:ice, iblk)
    wfract_can    => veg_mem  %wfract_can%ptr       (ics:ice, iblk)
    wsr_can       => veg_mem  %wsr_can%ptr          (ics:ice, iblk)
    wfract_soil   => soil_mem %wfract_soil%ptr      (ics:ice, iblk)
    wsr_soil      => soil_mem %wsr_soil%ptr         (ics:ice, iblk)
    sfract_srf    => srf_mem  %sfract_srf%ptr       (ics:ice, iblk)
    wfract_srf    => srf_mem  %wfract_srf%ptr       (ics:ice, iblk)
    oro_stddev    => srf_mem  %oro_stddev%ptr       (ics:ice, iblk)
    drag_srf      => atm2land_mem     %drag_srf%ptr         (ics:ice, iblk)
    q_air         => atm2land_mem     %q_air%ptr            (ics:ice, iblk)
    press_srf     => atm2land_mem     %press_srf%ptr        (ics:ice, iblk)
    t_srf         => srf_mem  %t_srf%ptr            (ics:ice, iblk)

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
    Def_config(srf_conf, srf)
    Def_memory(srf_mem, srf)
    Def_memory(soil_mem, soil)
    Def_memory(atm2land_mem, atm2land)

    LOGICAL                       :: lstart

    ! Pointers to variables in memory
    REAL(wp), POINTER ::       &
      & t_srf          (:),    &
      & t_srf_old      (:),    &
      & t_srf_unfilt   (:),    &
      & t_srf_unfilt_old(:),    &
      & t_eff_srf      (:),    &
      & qsat_srf       (:),    &
!!$      & qsat_srf_old   (:),    &
      & s_srf          (:),    &
!!$      & s_srf_old      (:),    &
      & press_srf      (:),    &
      & q_air          (:),    &
      & fact_qsat_srf  (:),    &
      & fact_q_air     (:),    &
      & drag_srf       (:),    &
      & rad_srf_net    (:),    &
      & grnd_hflx      (:),    &
      & hcap_grnd      (:),    &
      & wsn_srf        (:),    &
      & sfract_srf     (:),    &
      & t_acoef        (:),    &
      & t_bcoef        (:),    &
      & q_acoef        (:),    &
      & q_bcoef        (:),    &
      & glacier        (:)

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

    Get_config(srf_conf, SRF)

    ! Get reference to variables for current block
    !
    Get_memory(srf_mem, SRF)
    Get_memory(soil_mem, SOIL)
    Get_memory(atm2land_mem, ATM2LAND)

    t_srf                => srf_mem  %t_srf%ptr            (ics:ice, iblk)
    t_srf_old            => srf_mem  %t_srf_old%ptr        (ics:ice, iblk)
    t_srf_unfilt         => srf_mem  %t_srf_unfilt%ptr     (ics:ice, iblk)
    t_srf_unfilt_old     => srf_mem  %t_srf_unfilt_old%ptr (ics:ice, iblk)
    t_eff_srf            => srf_mem  %t_eff_srf%ptr        (ics:ice, iblk)
    qsat_srf             => srf_mem  %qsat_srf%ptr         (ics:ice, iblk)
!!$    qsat_srf_old         => srf_mem  %qsat_srf_old%ptr     (ics:ice, iblk)
    s_srf                => srf_mem  %s_srf%ptr            (ics:ice, iblk)
!!$    s_srf_old            => srf_mem  %s_srf_old%ptr        (ics:ice, iblk)
    press_srf            => atm2land_mem     %press_srf%ptr        (ics:ice, iblk)
    q_air                => atm2land_mem     %q_air%ptr            (ics:ice, iblk)
    fact_q_air           => srf_mem  %fact_q_air%ptr       (ics:ice, iblk)
    fact_qsat_srf        => srf_mem  %fact_qsat_srf%ptr    (ics:ice, iblk)
    drag_srf             => atm2land_mem     %drag_srf%ptr         (ics:ice, iblk)
    rad_srf_net          => srf_mem  %rad_srf_net%ptr      (ics:ice, iblk)
    grnd_hflx            => soil_mem %grnd_hflx%ptr        (ics:ice, iblk)
    hcap_grnd            => soil_mem %hcap_grnd%ptr        (ics:ice, iblk)
    wsn_srf              => srf_mem  %wsn_srf%ptr          (ics:ice, iblk)
    sfract_srf           => srf_mem  %sfract_srf%ptr       (ics:ice, iblk)
    t_acoef              => atm2land_mem     %t_acoef%ptr          (ics:ice, iblk)
    t_bcoef              => atm2land_mem     %t_bcoef%ptr          (ics:ice, iblk)
    q_acoef              => atm2land_mem     %q_acoef%ptr          (ics:ice, iblk)
    q_bcoef              => atm2land_mem     %q_bcoef%ptr          (ics:ice, iblk)
    glacier              => srf_mem  %glacier%ptr          (ics:ice, iblk)

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
    Def_memory(srf_mem, srf)
    LOGICAL                       :: lstart
    REAL(wp)                      :: eps

    ! Pointers to variables in memory
    REAL(wp), POINTER ::       &
      & t_srf          (:),    &
      & t_srf_old      (:),    &
      & t_srf_unfilt   (:),    &
      & t_srf_unfilt_old(:)

    Get_memory(srf_mem, SRF)

    t_srf                => srf_mem  %t_srf%ptr            (ics:ice, iblk)
    t_srf_old            => srf_mem  %t_srf_old%ptr        (ics:ice, iblk)
    t_srf_unfilt         => srf_mem  %t_srf_unfilt%ptr     (ics:ice, iblk)
    t_srf_unfilt_old     => srf_mem  %t_srf_unfilt_old%ptr (ics:ice, iblk)

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
    Def_memory(srf_mem, srf)
    Def_memory(atm2land_mem, ATM2LAND)

    ! Pointers to variables in memory
    REAL(wp), POINTER ::       &
      & qsat_srf       (:),    &
      & s_srf          (:),    &
      & t_acoef        (:),    &
      & t_bcoef        (:),    &
      & q_acoef        (:),    &
      & q_bcoef        (:),    &
      & fact_qsat_srf  (:),    &
      & fact_qsat_trans_srf  (:),    &
      & fact_q_air     (:),    &
      & drag_srf       (:),    &
      & latent_hflx    (:),    &
      & sensible_hflx  (:),    &
      & evapopot       (:),    &
      & evapotrans     (:),    &
      & trans          (:),    &
      & sfract_srf     (:)

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
    Get_memory(srf_mem, SRF)
    Get_memory(atm2land_mem, ATM2LAND)

    qsat_srf             => srf_mem  %qsat_srf%ptr         (ics:ice, iblk)
    s_srf                => srf_mem  %s_srf%ptr            (ics:ice, iblk)
    t_acoef              => atm2land_mem     %t_acoef%ptr          (ics:ice, iblk)
    t_bcoef              => atm2land_mem     %t_bcoef%ptr          (ics:ice, iblk)
    q_acoef              => atm2land_mem     %q_acoef%ptr          (ics:ice, iblk)
    q_bcoef              => atm2land_mem     %q_bcoef%ptr          (ics:ice, iblk)
    fact_q_air           => srf_mem  %fact_q_air%ptr       (ics:ice, iblk)
    fact_qsat_srf        => srf_mem  %fact_qsat_srf%ptr    (ics:ice, iblk)
    fact_qsat_trans_srf  => srf_mem  %fact_qsat_trans_srf%ptr (ics:ice, iblk)
    drag_srf             => atm2land_mem     %drag_srf%ptr         (ics:ice, iblk)
    sensible_hflx        => srf_mem  %sensible_hflx%ptr    (ics:ice, iblk)
    latent_hflx          => srf_mem  %latent_hflx%ptr      (ics:ice, iblk)
    evapotrans           => srf_mem  %evapotrans%ptr       (ics:ice, iblk)
    trans                => srf_mem  %trans%ptr            (ics:ice, iblk)
    evapopot             => srf_mem  %evapopot%ptr         (ics:ice, iblk)
    sfract_srf           => srf_mem  %sfract_srf%ptr       (ics:ice, iblk)

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
    Def_memory(srf_mem, srf)
    Def_memory(soil_mem, soil)
    Def_memory(veg_mem, veg)
    Def_memory(atm2land_mem, atm2land)

    ! Pointers to variables in memory
    REAL(wp), POINTER :: &
      fract_fpc     (:), &
      ws_l          (:,:), &
      ws_fc_l       (:,:), &
      ws_pwp_l      (:,:), &
!!$      sfract_can    (:), &
!!$      sfract_soil   (:), &
!!$      wfract_can    (:), &
!!$      wfract_soil   (:), &
      sfract_srf    (:), &
      wfract_srf    (:), &
      fact_q_air    (:), &
      fact_qsat_srf (:), &
      fact_qsat_trans_srf (:), &
      qsat_srf      (:), &
      canopy_cond   (:), &
      pch           (:), &
      wind_air      (:), &
      q_air         (:)

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
    Get_memory(srf_mem, SRF)
    Get_memory(soil_mem, SOIL)
    Get_memory(veg_mem, VEG)
    Get_memory(atm2land_mem, ATM2LAND)

    fract_fpc        => veg_mem  %fract_fpc%ptr           (ics:ice,    iblk)
    canopy_cond      => veg_mem  %canopy_cond%ptr         (ics:ice,    iblk)
    ws_l             => soil_mem %ws_l%ptr                (ics:ice, :, iblk)
    ws_fc_l          => soil_mem %ws_fc_l%ptr             (ics:ice, :, iblk)
    ws_pwp_l         => soil_mem %ws_pwp_l%ptr            (ics:ice, :, iblk)
!!$    sfract_can       => tile%mem%veg  %sfract_can%ptr          (ics:ice,    iblk)
!!$    sfract_soil      => tile%mem%soil %sfract_soil%ptr         (ics:ice,    iblk)
!!$    wfract_can       => tile%mem%veg  %wfract_can%ptr          (ics:ice,    iblk)
!!$    wfract_soil      => tile%mem%soil %wfract_soil%ptr         (ics:ice,    iblk)
    qsat_srf         => srf_mem  %qsat_srf%ptr            (ics:ice,    iblk)
    sfract_srf       => srf_mem  %sfract_srf%ptr          (ics:ice,    iblk)
    wfract_srf       => srf_mem  %wfract_srf%ptr          (ics:ice,    iblk)
    fact_q_air       => srf_mem  %fact_q_air%ptr          (ics:ice,    iblk)
    fact_qsat_srf    => srf_mem  %fact_qsat_srf%ptr       (ics:ice,    iblk)
    fact_qsat_trans_srf  => srf_mem  %fact_qsat_trans_srf%ptr (ics:ice,    iblk)
    pch              => atm2land_mem     %pch%ptr                 (ics:ice,    iblk)
    wind_air         => atm2land_mem     %wind_air%ptr            (ics:ice,    iblk)
    q_air            => atm2land_mem     %q_air%ptr               (ics:ice,    iblk)

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

    Def_memory(srf_mem, srf)
    Def_memory(soil_mem, soil)

    ! Pointers to variables in memory
    REAL(wp), POINTER ::       &
      & rough_m_srf    (:),    &
      & rough_h_srf    (:),    &
      & sfract_soil    (:)

    Get_memory(srf_mem, SRF)
    Get_memory(soil_mem, SOIL)

    rough_m_srf          => srf_mem  %rough_m_srf%ptr      (ics:ice, iblk)
    rough_h_srf          => srf_mem  %rough_h_srf%ptr      (ics:ice, iblk)
    sfract_soil          => soil_mem %sfract_soil%ptr      (ics:ice, iblk)

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
