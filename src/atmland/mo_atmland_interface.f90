!>
!! @brief Interface to run JSBACH for one time step
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_atmland_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_tile_class, ONLY: t_jsb_tile
  USE mo_exception,      ONLY: finish
  USE mo_kind,           ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_atm2land, update_land2atm

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_interface'

CONTAINS

  SUBROUTINE update_atm2land(tile, iblk, ics, ice, &
    & t_air, q_air, rain, snow, wind_air, wind_10m, lwrad_srf_down, swrad_srf_net, &
    & press_srf, drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle)

    TYPE(t_jsb_tile),   INTENT(inout) :: tile
    INTEGER,            INTENT(in)    :: iblk, ics, ice
    REAL(wp), OPTIONAL, INTENT(in)    :: &
      & t_air(:), q_air(:), rain(:), snow(:), wind_air(:), wind_10m(:), lwrad_srf_down(:), swrad_srf_net(:), &
      & press_srf(:), drag_srf(:), t_acoef(:), t_bcoef(:), q_acoef(:), q_bcoef(:), pch(:), cos_zenith_angle(:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_atm2land'

    IF (ice-ics+1 /= SIZE(t_air,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    IF (PRESENT(t_air))            tile%mem%atm2land%t_air            (ics:ice, iblk) = t_air            (:)
    IF (PRESENT(q_air))            tile%mem%atm2land%q_air            (ics:ice, iblk) = q_air            (:)
    IF (PRESENT(rain))             tile%mem%atm2land%rain             (ics:ice, iblk) = rain             (:)
    IF (PRESENT(snow))             tile%mem%atm2land%snow             (ics:ice, iblk) = snow             (:)
    IF (PRESENT(wind_air))         tile%mem%atm2land%wind_air         (ics:ice, iblk) = wind_air         (:)
    IF (PRESENT(wind_10m))         tile%mem%atm2land%wind_10m         (ics:ice, iblk) = wind_10m         (:)
    IF (PRESENT(lwrad_srf_down))   tile%mem%atm2land%lwrad_srf_down   (ics:ice, iblk) = lwrad_srf_down   (:)
    IF (PRESENT(swrad_srf_net))    tile%mem%atm2land%swrad_srf_net    (ics:ice, iblk) = swrad_srf_net    (:)
    IF (PRESENT(press_srf))        tile%mem%atm2land%press_srf        (ics:ice, iblk) = press_srf        (:)
    IF (PRESENT(drag_srf))         tile%mem%atm2land%drag_srf         (ics:ice, iblk) = drag_srf         (:)
    IF (PRESENT(t_acoef))          tile%mem%atm2land%t_acoef          (ics:ice, iblk) = t_acoef          (:)
    IF (PRESENT(t_bcoef))          tile%mem%atm2land%t_bcoef          (ics:ice, iblk) = t_bcoef          (:)
    IF (PRESENT(q_acoef))          tile%mem%atm2land%q_acoef          (ics:ice, iblk) = q_acoef          (:)
    IF (PRESENT(q_bcoef))          tile%mem%atm2land%q_bcoef          (ics:ice, iblk) = q_bcoef          (:)
    IF (PRESENT(pch))              tile%mem%atm2land%pch              (ics:ice, iblk) = pch              (:)
    IF (PRESENT(cos_zenith_angle)) tile%mem%atm2land%cos_zenith_angle (ics:ice, iblk) = cos_zenith_angle (:)

  END SUBROUTINE update_atm2land

  SUBROUTINE update_land2atm(tile, iblk, ics, ice,           &
    & t_srf, t_eff_srf, qsat_srf, s_srf,                     &
    & fact_q_air, fact_qsat_srf,                             &
    & evapotrans, latent_hflx, sensible_hflx,                &
    & rough_h_srf, rough_m_srf, tte_corr,                    &
    & alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif)

    TYPE(t_jsb_tile),   INTENT(in) :: tile
    INTEGER,            INTENT(in) :: iblk, ics, ice
    REAL(wp), OPTIONAL, INTENT(out) ::                                        &
      & t_srf(:), t_eff_srf(:), qsat_srf(:), s_srf(:),                        &
      & fact_q_air(:), fact_qsat_srf(:),                                      &
      & evapotrans(:), latent_hflx(:), sensible_hflx(:),                      &
      & rough_h_srf(:), rough_m_srf(:), tte_corr(:),                          &
      & alb_vis_dir(:), alb_nir_dir(:), alb_vis_dif(:), alb_nir_dif(:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_land2atm'

    IF (ice-ics+1 /= SIZE(t_srf,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    IF (PRESENT(t_srf))           t_srf(:)            = tile%mem%srf  %energy_state %t_srf          (ics:ice, iblk)
    IF (PRESENT(t_eff_srf))       t_eff_srf(:)        = tile%mem%srf  %energy_state %t_eff_srf      (ics:ice, iblk)
    IF (PRESENT(qsat_srf))        qsat_srf(:)         = tile%mem%srf  %energy_state %qsat_srf       (ics:ice, iblk)
    IF (PRESENT(s_srf))           s_srf(:)            = tile%mem%srf  %energy_state %s_srf          (ics:ice, iblk)
    IF (PRESENT(fact_q_air))      fact_q_air(:)       = tile%mem%srf  %energy_state %fact_q_air     (ics:ice, iblk)
    IF (PRESENT(fact_qsat_srf))   fact_qsat_srf(:)    = tile%mem%srf  %energy_state %fact_qsat_srf  (ics:ice, iblk)
    IF (PRESENT(evapotrans))      evapotrans(:)       = tile%mem%srf  %water_flux   %evapotrans     (ics:ice, iblk)
    IF (PRESENT(latent_hflx))     latent_hflx(:)      = tile%mem%srf  %energy_flux  %latent_hflx    (ics:ice, iblk)
    IF (PRESENT(sensible_hflx))   sensible_hflx(:)    = tile%mem%srf  %energy_flux  %sensible_hflx  (ics:ice, iblk)
    IF (PRESENT(rough_h_srf))     rough_h_srf(:)      = tile%mem%srf  %phys_state   %rough_h_srf    (ics:ice, iblk)
    IF (PRESENT(rough_m_srf))     rough_m_srf(:)      = tile%mem%srf  %phys_state   %rough_m_srf    (ics:ice, iblk)
    IF (PRESENT(tte_corr))        tte_corr(:)         = tile%mem%srf  %energy_state %tte_corr       (ics:ice, iblk)
    IF (PRESENT(alb_vis_dir))     alb_vis_dir(:)      = tile%mem%srf  %phys_state   %alb_vis_srf    (ics:ice, iblk)  ! Temporary
    IF (PRESENT(alb_vis_dif))     alb_vis_dif(:)      = tile%mem%srf  %phys_state   %alb_vis_srf    (ics:ice, iblk)  ! Temporary
    IF (PRESENT(alb_nir_dir))     alb_nir_dir(:)      = tile%mem%srf  %phys_state   %alb_nir_srf    (ics:ice, iblk)  ! Temporary
    IF (PRESENT(alb_nir_dif))     alb_nir_dif(:)      = tile%mem%srf  %phys_state   %alb_nir_srf    (ics:ice, iblk)  ! Temporary

  END SUBROUTINE update_land2atm

#endif
END MODULE mo_atmland_interface
