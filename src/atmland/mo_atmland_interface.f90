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

  dsl4jsb_Use_processes ATM2LAND, SRF, SOIL
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SRF)
  dsl4jsb_Use_memory(SOIL)

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

    dsl4jsb_Def_memory(ATM2LAND)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_atm2land'

    IF (ice-ics+1 /= SIZE(t_air,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl4jsb_Get_memory(ATM2LAND)

    IF (PRESENT(t_air))            dsl4jsb_var2D_onChunk(ATM2LAND,t_air)            = t_air            (:)
    IF (PRESENT(q_air))            dsl4jsb_var2D_onChunk(ATM2LAND,q_air)            = q_air            (:)
    IF (PRESENT(rain))             dsl4jsb_var2D_onChunk(ATM2LAND,rain)             = rain             (:)
    IF (PRESENT(snow))             dsl4jsb_var2D_onChunk(ATM2LAND,snow)             = snow             (:)
    IF (PRESENT(wind_air))         dsl4jsb_var2D_onChunk(ATM2LAND,wind_air)         = wind_air         (:)
    IF (PRESENT(wind_10m))         dsl4jsb_var2D_onChunk(ATM2LAND,wind_10m)         = wind_10m         (:)
    IF (PRESENT(lw_srf_down))      dsl4jsb_var2D_onChunk(ATM2LAND,lw_srf_down)      = lw_srf_down      (:)
    IF (PRESENT(swvis_srf_down))   dsl4jsb_var2D_onChunk(ATM2LAND,swvis_srf_down)   = swvis_srf_down   (:)
    IF (PRESENT(swnir_srf_down))   dsl4jsb_var2D_onChunk(ATM2LAND,swnir_srf_down)   = swnir_srf_down   (:)
    IF (PRESENT(swpar_srf_down))   dsl4jsb_var2D_onChunk(ATM2LAND,swpar_srf_down)   = swpar_srf_down   (:)
    IF (PRESENT(press_srf))        dsl4jsb_var2D_onChunk(ATM2LAND,press_srf)        = press_srf        (:)
    IF (PRESENT(drag_srf))         dsl4jsb_var2D_onChunk(ATM2LAND,drag_srf)         = drag_srf         (:)
    IF (PRESENT(t_acoef))          dsl4jsb_var2D_onChunk(ATM2LAND,t_acoef)          = t_acoef          (:)
    IF (PRESENT(t_bcoef))          dsl4jsb_var2D_onChunk(ATM2LAND,t_bcoef)          = t_bcoef          (:)
    IF (PRESENT(q_acoef))          dsl4jsb_var2D_onChunk(ATM2LAND,q_acoef)          = q_acoef          (:)
    IF (PRESENT(q_bcoef))          dsl4jsb_var2D_onChunk(ATM2LAND,q_bcoef)          = q_bcoef          (:)
    IF (PRESENT(pch))              dsl4jsb_var2D_onChunk(ATM2LAND,pch)              = pch              (:)
    IF (PRESENT(cos_zenith_angle)) dsl4jsb_var2D_onChunk(ATM2LAND,cos_zenith_angle) = cos_zenith_angle (:)

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

    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_land2atm'

    IF (ice-ics+1 /= SIZE(t_srf,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)

    IF (PRESENT(t_srf))           t_srf(:)            = dsl4jsb_var2D_onChunk(SRF,t_srf)
    IF (PRESENT(t_eff_srf))       t_eff_srf(:)        = dsl4jsb_var2D_onChunk(SRF,t_eff_srf)
    IF (PRESENT(qsat_srf))        qsat_srf(:)         = dsl4jsb_var2D_onChunk(SRF,qsat_srf)
    IF (PRESENT(s_srf))           s_srf(:)            = dsl4jsb_var2D_onChunk(SRF,s_srf)
    IF (PRESENT(fact_q_air))      fact_q_air(:)       = dsl4jsb_var2D_onChunk(SRF,fact_q_air)
    IF (PRESENT(fact_qsat_srf))   fact_qsat_srf(:)    = dsl4jsb_var2D_onChunk(SRF,fact_qsat_srf)
    IF (PRESENT(evapotrans))      evapotrans(:)       = dsl4jsb_var2D_onChunk(SRF,evapotrans)
    IF (PRESENT(latent_hflx))     latent_hflx(:)      = dsl4jsb_var2D_onChunk(SRF,latent_hflx)
    IF (PRESENT(sensible_hflx))   sensible_hflx(:)    = dsl4jsb_var2D_onChunk(SRF,sensible_hflx)
    IF (PRESENT(grnd_hflx))       grnd_hflx(:)        = dsl4jsb_var2D_onChunk(SOIL,grnd_hflx_old)
    IF (PRESENT(grnd_hcap))       grnd_hcap(:)        = dsl4jsb_var2D_onChunk(SOIL,hcap_grnd_old)
    IF (PRESENT(rough_h_srf))     rough_h_srf(:)      = dsl4jsb_var2D_onChunk(SRF,rough_h_srf)
    IF (PRESENT(rough_m_srf))     rough_m_srf(:)      = dsl4jsb_var2D_onChunk(SRF,rough_m_srf)
    IF (PRESENT(tte_corr))        tte_corr(:)         = dsl4jsb_var2D_onChunk(SRF,tte_corr)
    IF (PRESENT(alb_vis_dir))     alb_vis_dir(:)      = dsl4jsb_var2D_onChunk(SRF,alb_vis_srf)  ! Temporary
    IF (PRESENT(alb_vis_dif))     alb_vis_dif(:)      = dsl4jsb_var2D_onChunk(SRF,alb_vis_srf)  ! Temporary
    IF (PRESENT(alb_nir_dir))     alb_nir_dir(:)      = dsl4jsb_var2D_onChunk(SRF,alb_nir_srf)  ! Temporary
    IF (PRESENT(alb_nir_dif))     alb_nir_dif(:)      = dsl4jsb_var2D_onChunk(SRF,alb_nir_srf)  ! Temporary

  END SUBROUTINE update_land2atm

#endif
END MODULE mo_atmland_interface
