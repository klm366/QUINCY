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

#include "dsl_jsbach.inc"

  USE mo_exception,      ONLY: finish
  USE mo_kind,           ONLY: wp
  USE mo_jsb_tile_class, ONLY: t_jsb_tile

  dsl$jsb Use_processes ATM2LAND, SRF, SOIL
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl$jsb Use_memory(atm2land)
  dsl$jsb Use_memory(srf)
  dsl$jsb Use_memory(soil)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_atm2land, update_land2atm

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_interface'

CONTAINS

  SUBROUTINE update_atm2land(tile, iblk, ics, ice, &
    & t_air, q_air, rain, snow, wind_air, wind_10m, lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, &
    & press_srf, drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle)

    TYPE(t_jsb_tile),   INTENT(inout) :: tile
    INTEGER,            INTENT(in)    :: iblk, ics, ice
    REAL(wp), OPTIONAL, INTENT(in)    :: &
      & t_air(:), q_air(:), rain(:), snow(:), wind_air(:), wind_10m(:), lw_srf_down(:), swvis_srf_down(:), swnir_srf_down(:), &
      & swpar_srf_down(:), press_srf(:), drag_srf(:), t_acoef(:), t_bcoef(:), q_acoef(:), q_bcoef(:), pch(:), cos_zenith_angle(:)

    dsl$jsb Def_memory(atm2land_mem, atm2land)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_atm2land'

    IF (ice-ics+1 /= SIZE(t_air,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl$jsb Get_memory(atm2land_mem, ATM2LAND)

    IF (PRESENT(t_air))            atm2land_mem%t_air%ptr            (ics:ice, iblk) = t_air            (:)
    IF (PRESENT(q_air))            atm2land_mem%q_air%ptr            (ics:ice, iblk) = q_air            (:)
    IF (PRESENT(rain))             atm2land_mem%rain%ptr             (ics:ice, iblk) = rain             (:)
    IF (PRESENT(snow))             atm2land_mem%snow%ptr             (ics:ice, iblk) = snow             (:)
    IF (PRESENT(wind_air))         atm2land_mem%wind_air%ptr         (ics:ice, iblk) = wind_air         (:)
    IF (PRESENT(wind_10m))         atm2land_mem%wind_10m%ptr         (ics:ice, iblk) = wind_10m         (:)
    IF (PRESENT(lw_srf_down))      atm2land_mem%lw_srf_down%ptr      (ics:ice, iblk) = lw_srf_down      (:)
    IF (PRESENT(swvis_srf_down))   atm2land_mem%swvis_srf_down%ptr   (ics:ice, iblk) = swvis_srf_down   (:)
    IF (PRESENT(swnir_srf_down))   atm2land_mem%swnir_srf_down%ptr   (ics:ice, iblk) = swnir_srf_down   (:)
    IF (PRESENT(swpar_srf_down))   atm2land_mem%swpar_srf_down%ptr   (ics:ice, iblk) = swpar_srf_down   (:)
    IF (PRESENT(press_srf))        atm2land_mem%press_srf%ptr        (ics:ice, iblk) = press_srf        (:)
    IF (PRESENT(drag_srf))         atm2land_mem%drag_srf%ptr         (ics:ice, iblk) = drag_srf         (:)
    IF (PRESENT(t_acoef))          atm2land_mem%t_acoef%ptr          (ics:ice, iblk) = t_acoef          (:)
    IF (PRESENT(t_bcoef))          atm2land_mem%t_bcoef%ptr          (ics:ice, iblk) = t_bcoef          (:)
    IF (PRESENT(q_acoef))          atm2land_mem%q_acoef%ptr          (ics:ice, iblk) = q_acoef          (:)
    IF (PRESENT(q_bcoef))          atm2land_mem%q_bcoef%ptr          (ics:ice, iblk) = q_bcoef          (:)
    IF (PRESENT(pch))              atm2land_mem%pch%ptr              (ics:ice, iblk) = pch              (:)
    IF (PRESENT(cos_zenith_angle)) atm2land_mem%cos_zenith_angle%ptr (ics:ice, iblk) = cos_zenith_angle (:)

  END SUBROUTINE update_atm2land

  SUBROUTINE update_land2atm(tile, iblk, ics, ice,                      &
    & t_srf, t_eff_srf, qsat_srf, s_srf,                                &
    & fact_q_air, fact_qsat_srf,                                        &
    & evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap,     &
    & rough_h_srf, rough_m_srf, tte_corr,                               &
    & alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif)

    TYPE(t_jsb_tile),   INTENT(in) :: tile
    INTEGER,            INTENT(in) :: iblk, ics, ice
    REAL(wp), OPTIONAL, INTENT(out) ::                                               &
      & t_srf(:), t_eff_srf(:), qsat_srf(:), s_srf(:),                               &
      & fact_q_air(:), fact_qsat_srf(:),                                             &
      & evapotrans(:), latent_hflx(:), sensible_hflx(:), grnd_hflx(:), grnd_hcap(:), &
      & rough_h_srf(:), rough_m_srf(:), tte_corr(:),                                 &
      & alb_vis_dir(:), alb_nir_dir(:), alb_vis_dif(:), alb_nir_dif(:)

    dsl$jsb Def_memory(srf_mem, srf)
    dsl$jsb Def_memory(soil_mem, soil)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_land2atm'

    IF (ice-ics+1 /= SIZE(t_srf,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl$jsb Get_memory(srf_mem, SRF)
    dsl$jsb Get_memory(soil_mem, SOIL)

    IF (PRESENT(t_srf))           t_srf(:)            = srf_mem  %t_srf%ptr          (ics:ice, iblk)
    IF (PRESENT(t_eff_srf))       t_eff_srf(:)        = srf_mem  %t_eff_srf%ptr      (ics:ice, iblk)
    IF (PRESENT(qsat_srf))        qsat_srf(:)         = srf_mem  %qsat_srf%ptr       (ics:ice, iblk)
    IF (PRESENT(s_srf))           s_srf(:)            = srf_mem  %s_srf%ptr          (ics:ice, iblk)
    IF (PRESENT(fact_q_air))      fact_q_air(:)       = srf_mem  %fact_q_air%ptr     (ics:ice, iblk)
    IF (PRESENT(fact_qsat_srf))   fact_qsat_srf(:)    = srf_mem  %fact_qsat_srf%ptr  (ics:ice, iblk)
    IF (PRESENT(evapotrans))      evapotrans(:)       = srf_mem  %evapotrans%ptr     (ics:ice, iblk)
    IF (PRESENT(latent_hflx))     latent_hflx(:)      = srf_mem  %latent_hflx%ptr    (ics:ice, iblk)
    IF (PRESENT(sensible_hflx))   sensible_hflx(:)    = srf_mem  %sensible_hflx%ptr  (ics:ice, iblk)
    IF (PRESENT(grnd_hflx))       grnd_hflx(:)        = soil_mem %grnd_hflx_old%ptr  (ics:ice, iblk)
    IF (PRESENT(grnd_hcap))       grnd_hcap(:)        = soil_mem %hcap_grnd_old%ptr  (ics:ice, iblk)
    IF (PRESENT(rough_h_srf))     rough_h_srf(:)      = srf_mem  %rough_h_srf%ptr    (ics:ice, iblk)
    IF (PRESENT(rough_m_srf))     rough_m_srf(:)      = srf_mem  %rough_m_srf%ptr    (ics:ice, iblk)
    IF (PRESENT(tte_corr))        tte_corr(:)         = srf_mem  %tte_corr%ptr       (ics:ice, iblk)
    IF (PRESENT(alb_vis_dir))     alb_vis_dir(:)      = srf_mem  %alb_vis_srf%ptr    (ics:ice, iblk)  ! Temporary
    IF (PRESENT(alb_vis_dif))     alb_vis_dif(:)      = srf_mem  %alb_vis_srf%ptr    (ics:ice, iblk)  ! Temporary
    IF (PRESENT(alb_nir_dir))     alb_nir_dir(:)      = srf_mem  %alb_nir_srf%ptr    (ics:ice, iblk)  ! Temporary
    IF (PRESENT(alb_nir_dif))     alb_nir_dif(:)      = srf_mem  %alb_nir_srf%ptr    (ics:ice, iblk)  ! Temporary

  END SUBROUTINE update_land2atm

#endif
END MODULE mo_atmland_interface
