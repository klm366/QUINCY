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

  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_tile_class,  ONLY: t_jsb_tile
  USE mo_jsb_base,        ONLY: get_model
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message
  USE mo_util_string,     ONLY: one_of

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_surface_radiation, update_snow_and_skin_fraction, update_surface_energy !, update_srf

  CHARACTER(len=*), PARAMETER :: modname = 'mo_srf_interface'

CONTAINS

  SUBROUTINE update_surface_radiation(tile, iblk, ics, ice, nc)

    USE mo_phy_schemes, ONLY: thermal_radiation

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model

    ! Pointers to variables in memory
    REAL(wp), POINTER ::      &
      & swrad_srf_net    (:), &
      & lwrad_srf_net    (:), &
      & rad_srf_net      (:), &
      & swrad_srf_down   (:), &
      & lwrad_srf_down   (:), &
      & alb_vis_srf      (:), &
      & alb_nir_srf      (:), &
      & t_srf            (:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_radiation'

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of surface radiation')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    swrad_srf_net  => tile%mem%srf %energy_flux  %swrad_srf_net  (ics:ice, iblk)
    swrad_srf_down => tile%mem     %atm2land     %swrad_srf_down (ics:ice, iblk)
    lwrad_srf_net  => tile%mem%srf %energy_flux  %lwrad_srf_net  (ics:ice, iblk)
    lwrad_srf_down => tile%mem     %atm2land     %lwrad_srf_down (ics:ice, iblk)
    rad_srf_net    => tile%mem%srf %energy_flux  %rad_srf_net    (ics:ice, iblk)
    alb_vis_srf    => tile%mem%srf %phys_state   %alb_vis_srf    (ics:ice, iblk)
    alb_nir_srf    => tile%mem%srf %phys_state   %alb_nir_srf    (ics:ice, iblk)
    lwrad_srf_net  => tile%mem%srf %energy_flux  %lwrad_srf_net  (ics:ice, iblk)
    t_srf          => tile%mem%srf %energy_state %t_srf          (ics:ice, iblk)

    ! Compute SW net radiation from incoming and albedo
    swrad_srf_net(:) = swrad_srf_down(:) * (1._wp - (alb_vis_srf(:) + alb_nir_srf(:)) / 2._wp)

    ! Compute LW net radiation from incoming and the thermal radiation from the surface
    lwrad_srf_net(:) = lwrad_srf_down(:) - thermal_radiation(t_srf(:))

    ! Compute net radiation
    rad_srf_net(:) = swrad_srf_net(:) + lwrad_srf_net(:)

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

    ! Pointers to variables in memory
    REAL(wp), POINTER :: &
                                !wsn           (:), &
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

    IF (one_of('soil', tile%components) < 1) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of snow and skin reservoir fraction')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    lai           => tile%mem%veg  %phys_state   %lai              (ics:ice, iblk)
    fract_fpc_max => tile%mem%veg  %phys_state   %fract_fpc_max    (ics:ice, iblk)
    sfract_can    => tile%mem%veg  %water_state  %sfract_can       (ics:ice, iblk)
    wsn_can       => tile%mem%veg  %water_state  %wsn_can          (ics:ice, iblk)
    sfract_soil   => tile%mem%soil %water_state  %sfract_soil      (ics:ice, iblk)
    wsn_soil      => tile%mem%soil %water_state  %wsn_soil         (ics:ice, iblk)
    wfract_can    => tile%mem%veg  %water_state  %wfract_can       (ics:ice, iblk)
    wsr_can       => tile%mem%veg  %water_state  %wsr_can          (ics:ice, iblk)
    wfract_soil   => tile%mem%soil %water_state  %wfract_soil      (ics:ice, iblk)
    wsr_soil      => tile%mem%soil %water_state  %wsr_soil         (ics:ice, iblk)
    sfract_srf    => tile%mem%srf  %water_state  %sfract_srf       (ics:ice, iblk)
    wfract_srf    => tile%mem%srf  %water_state  %wfract_srf       (ics:ice, iblk)
    oro_stddev    => tile%mem%srf  %phys_state   %oro_stddev       (ics:ice, iblk)
    drag_srf      => tile%mem      %atm2land     %drag_srf         (ics:ice, iblk)
    q_air         => tile%mem      %atm2land     %q_air            (ics:ice, iblk)
    press_srf     => tile%mem      %atm2land     %press_srf        (ics:ice, iblk)
    t_srf         => tile%mem%srf  %energy_state %t_srf            (ics:ice, iblk)

    ! Maximum capacity of skin reservoir (soil + canopy)
    wsr_max(:) = model%config%soil%wsr_capacity * (1._wp + lai(:) * fract_fpc_max(:))

    ! Wet skin fraction on soil
    wfract_soil(:) = MERGE(MIN(1._wp, wsr_soil(:) / wsr_max(:)), 0._wp, wsr_max(:) > EPSILON(1._wp))

    ! Wet skin fraction on canopy
    wfract_can(:) = MERGE(MIN(1._wp, wsr_can(:) / wsr_max(:)), 0._wp, wsr_max(:) > EPSILON(1._wp))

    ! Snow fraction on soil
    sfract_soil(:) = snow_fract(wsr_soil(:), oro_stddev(:))
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
    evap_pot(:) =  heat_transfer_coef(drag_srf(:), steplen) * (qsat_srf_old(:) - q_air(:))
    !
    WHERE (sfract_soil(:) > 0._wp)
      sfract_soil(:) = sfract_soil(:) / MAX(1._wp, sfract_soil(:) * evap_pot(:) * dtime / (rhoh2o * (wsn_soil(:) + wsn_can(:))))
    END WHERE
    WHERE (wfract_soil(:) > 0._wp)
      wfract_soil(:) = wfract_soil(:) / MAX(1._wp, wfract_soil(:) * evap_pot(:) * dtime / (rhoh2o * (wsr_soil(:) + wsr_can(:))))
    END WHERE

    wfract_srf(:) = wfract_soil(:) + wfract_can(:)
    sfract_srf(:) = sfract_soil(:) + sfract_can(:)

  END SUBROUTINE update_snow_and_skin_fraction

  SUBROUTINE update_surface_energy( tile, iblk, ics, ice, nc, dtime, steplen)

    USE mo_phy_schemes,        ONLY: qsat, surface_dry_static_energy, heat_transfer_coef
    USE mo_srf_process,        ONLY: surface_temp_and_qsat_implicit
    USE mo_srf_constants,      ONLY: eps
    USE mo_echam_vdiff_params, ONLY: tpfac2, tpfac3
    USE mo_physical_constants, ONLY: tmelt

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc
    REAL(wp),          INTENT(in) :: dtime, steplen

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    LOGICAL                       :: lstart

    ! Pointers to variables in memory
    REAL(wp), POINTER ::       &
      & fract          (:),    &
      & t_srf          (:),    &
      & t_srf_unfilt   (:),    &
      & t_rad_srf      (:),    &
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
      & q_bcoef        (:)

    ! Locally allocated vectors
    !
    REAL(wp) ::             &
      t_srf_old(nc),        &
      t_srf_unfilt_old(nc), &
      qsat_srf_old(nc),     &
      s_srf_old(nc),        &
      dQdT(nc),             & !< Sensitivity of saturated surface specific humidity to temperature
      t2s_conv(nc),         & !< Conversion factor from temperature to dry static energy (C_pd * (1+(delta-1)*q_v))
      heat_tcoef(nc)          !< Heat transfer coefficient (rho*C_h*|v|)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_energy'

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of surface temperature')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    t_srf                => tile%mem%srf  %energy_state %t_srf            (ics:ice, iblk)
    t_srf_unfilt         => tile%mem%srf  %energy_state %t_srf_unfilt     (ics:ice, iblk)
    t_rad_srf            => tile%mem%srf  %energy_state %t_rad_srf        (ics:ice, iblk)
    qsat_srf             => tile%mem%srf  %energy_state %qsat_srf         (ics:ice, iblk)
!!$    qsat_srf_old         => tile%mem%srf  %energy_state %qsat_srf_old     (ics:ice, iblk)
    s_srf                => tile%mem%srf  %energy_state %s_srf            (ics:ice, iblk)
!!$    s_srf_old            => tile%mem%srf  %energy_state %s_srf_old        (ics:ice, iblk)
    press_srf            => tile%mem      %atm2land     %press_srf        (ics:ice, iblk)
    q_air                => tile%mem      %atm2land     %q_air            (ics:ice, iblk)
    fact_q_air           => tile%mem%srf  %energy_state %fact_q_air       (ics:ice, iblk)
    fact_qsat_srf        => tile%mem%srf  %energy_state %fact_qsat_srf    (ics:ice, iblk)
    drag_srf             => tile%mem      %atm2land     %drag_srf         (ics:ice, iblk)
    rad_srf_net          => tile%mem%srf  %energy_flux  %rad_srf_net      (ics:ice, iblk)
    grnd_hflx            => tile%mem%soil %energy_flux  %grnd_hflx        (ics:ice, iblk)
    hcap_grnd            => tile%mem%soil %energy_state %hcap_grnd        (ics:ice, iblk)
    wsn_srf              => tile%mem%srf  %water_state  %wsn_srf          (ics:ice, iblk)
    sfract_srf           => tile%mem%srf  %water_state  %sfract_srf       (ics:ice, iblk)
    t_acoef              => tile%mem      %atm2land     %t_acoef          (ics:ice, iblk)
    t_bcoef              => tile%mem      %atm2land     %t_bcoef          (ics:ice, iblk)
    q_acoef              => tile%mem      %atm2land     %q_acoef          (ics:ice, iblk)
    q_bcoef              => tile%mem      %atm2land     %q_bcoef          (ics:ice, iblk)

    ! Save old surface temperature and saturation spec. humidity
    t_srf_unfilt_old(:) = t_srf_unfilt(:)
    t_srf_old   (:) = t_srf(:)
    !qsat_srf_old(:) = qsat_srf(:)   ! Is this the same as next line?
    qsat_srf_old(:) = qsat(t_srf_old(:), press_srf(:))

    ! Calculate surface saturation specific humidity

    ! Update surface temperature
    !
    ! 1. Compute sensitivity of surface saturation specific humidity to temperature
    dQdT(:) = ( qsat(t_srf(:) + 0.001_wp, press_srf(:)) - qsat_srf_old(:) ) * 1000._wp

    ! 2. Old dry static energy
    s_srf_old(:) = surface_dry_static_energy(t_srf_old(:), qsat_srf_old(:), q_air(:), fact_qsat_srf(:), fact_q_air(:))

    t2s_conv(:) = s_srf_old(:) / t_srf_old(:)  !

    heat_tcoef(:) = heat_transfer_coef(drag_srf(:), steplen)

    ! Compute new surface temperature and saturation humidity
    IF (ALL(qsat_srf(:) < 0._wp)) THEN          ! Start of experiment
      lstart = .TRUE.
      s_srf(:) = s_srf_old(:)
      qsat_srf(:) = qsat_srf_old(:)
    ELSE
      lstart = .FALSE.
      CALL surface_temp_and_qsat_implicit( nc, steplen, &
        t2s_conv(:), &
        t_acoef(:), t_bcoef(:), q_acoef(:), q_bcoef(:), &
        s_srf_old(:), qsat_srf_old(:), dQdT(:), &
        rad_srf_net(:), grnd_hflx(:), heat_tcoef(:), &
        fact_q_air(:), fact_qsat_srf(:), &
        sfract_srf(:), hcap_grnd(:), &
        & s_srf(:), qsat_srf(:))
    END IF

    ! Back out unfiltered surface temperature (see below Eq. 3.3.1.5 in ECHAM3 manual, alpha = 1/tpfac2)
    t_srf_unfilt(:) = ( tpfac2 * s_srf(:) + tpfac3 * s_srf_old(:) ) / t2s_conv(:)

    ! Correct new dry static energy for snowmelt
    WHERE (wsn_srf(:) > model%config%srf%wsn_critical)
      s_srf(:) = ( MIN(t_srf_unfilt(:), tmelt) * t2s_conv(:) - tpfac3 * s_srf_old(:) ) / tpfac2
    END WHERE

    ! 'Radiative temperature' used in the radiation scheme of the atmosphere model (see Eq. 6.3 in ECHAM5 manual)
    t_rad_srf(:) = s_srf(:) / t2s_conv(:)

    ! Time filter for surface temperature
    IF (.NOT. lstart) THEN
      t_srf(:) = t_srf_unfilt_old(:) + eps * (t_srf_old(:) - 2._wp * t_srf_unfilt_old(:) + t_srf_unfilt(:))
    ELSE
      t_srf(:) = t_srf_unfilt_old(:)
    END IF

  END SUBROUTINE update_surface_energy

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

END MODULE mo_srf_interface
