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

  USE mo_exception,       ONLY: finish
  USE mo_kind,            ONLY: wp

  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_memory_class,   ONLY: t_jsb_memory, t_jsb_memory_p
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
  USE mo_jsb_lct_class,      ONLY: LAKE_TYPE

  dsl4jsb_Use_processes ATM2LAND, SEB, TURB, SOIL, HYDRO, RAD
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(TURB)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(HYDRO)
  dsl4jsb_Use_memory(RAD)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_atm2land, update_land2atm

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_interface'

CONTAINS

  SUBROUTINE update_atm2land(tile, config, options,                &
    & t_air, q_air, rain, snow, wind_air, wind_10m,                &
    & lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, &
    & press_srf, drag_srf,                                         &
    & t_acoef, t_bcoef, q_acoef, q_bcoef,                          &
    & pch, cos_zenith_angle)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p), INTENT(in)         :: config(:)
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    REAL(wp), OPTIONAL, INTENT(in) ::                                            &
      & t_air(:), q_air(:), rain(:), snow(:), wind_air(:), wind_10m(:),          &
      & lw_srf_down(:), swvis_srf_down(:), swnir_srf_down(:), swpar_srf_down(:), &
      & press_srf(:), drag_srf(:),                                               &
      & t_acoef(:), t_bcoef(:), q_acoef(:), q_bcoef(:),                          &
      & pch(:), cos_zenith_angle(:)

    dsl4jsb_Def_memory(ATM2LAND)

    INTEGER ::  iblk, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_atm2land'

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

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

  SUBROUTINE update_land2atm(tile, config, options,                     &
    & t_srf, t_eff_srf, qsat_srf, s_srf,                                &
    & fact_q_air, fact_qsat_srf,                                        &
    & evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap,     &
    & rough_h_srf, rough_m_srf, tte_corr,                               &
    & alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif,               &
    & t_wtr, t_ice, ice_fract_lake)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    REAL(wp), OPTIONAL, INTENT(out) ::                                               &
      & t_srf(:), t_eff_srf(:), qsat_srf(:), s_srf(:),                               &
      & fact_q_air(:), fact_qsat_srf(:),                                             &
      & evapotrans(:), latent_hflx(:), sensible_hflx(:), grnd_hflx(:), grnd_hcap(:), &
      & rough_h_srf(:), rough_m_srf(:), tte_corr(:),                                 &
      & alb_vis_dir(:), alb_nir_dir(:), alb_vis_dif(:), alb_nir_dif(:),              &
      & t_wtr(:), t_ice(:), ice_fract_lake(:)

    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(TURB)
    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_memory(HYDRO)
    dsl4jsb_Def_memory(RAD)

    INTEGER ::  iblk, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_land2atm'

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (ice-ics+1 /= SIZE(t_srf,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl4jsb_Get_memory(SEB)
    dsl4jsb_Get_memory(TURB)
    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(HYDRO)
    dsl4jsb_Get_memory(RAD)

    IF (PRESENT(t_srf))           t_srf(:)            = dsl4jsb_var2D_onChunk(SEB,t)
    IF (PRESENT(t_eff_srf))       t_eff_srf(:)        = dsl4jsb_var2D_onChunk(SEB,t_eff)
    IF (PRESENT(qsat_srf))        qsat_srf(:)         = dsl4jsb_var2D_onChunk(SEB,qsat_star)
    IF (PRESENT(s_srf))           s_srf(:)            = dsl4jsb_var2D_onChunk(SEB,s_star)
    IF (PRESENT(fact_q_air))      fact_q_air(:)       = dsl4jsb_var2D_onChunk(SEB,fact_q_air)
    IF (PRESENT(fact_qsat_srf))   fact_qsat_srf(:)    = dsl4jsb_var2D_onChunk(SEB,fact_qsat_srf)
    IF (PRESENT(evapotrans))      evapotrans(:)       = dsl4jsb_var2D_onChunk(HYDRO,evapotrans)
    IF (PRESENT(latent_hflx))     latent_hflx(:)      = dsl4jsb_var2D_onChunk(SEB,latent_hflx)
    IF (PRESENT(sensible_hflx))   sensible_hflx(:)    = dsl4jsb_var2D_onChunk(SEB,sensible_hflx)
    IF (PRESENT(grnd_hflx))       grnd_hflx(:)        = dsl4jsb_var2D_onChunk(SEB,forc_hflx)
    IF (PRESENT(grnd_hcap))       grnd_hcap(:)        = dsl4jsb_var2D_onChunk(SEB,heat_cap_old)
    IF (PRESENT(rough_h_srf))     rough_h_srf(:)      = dsl4jsb_var2D_onChunk(TURB,rough_h)
    IF (PRESENT(rough_m_srf))     rough_m_srf(:)      = dsl4jsb_var2D_onChunk(TURB,rough_m)
    IF (PRESENT(tte_corr))        tte_corr(:)         = dsl4jsb_var2D_onChunk(HYDRO,tte_corr)
    IF (PRESENT(alb_vis_dir))     alb_vis_dir(:)      = dsl4jsb_var2D_onChunk(RAD,alb_vis)  ! Temporary
    IF (PRESENT(alb_vis_dif))     alb_vis_dif(:)      = dsl4jsb_var2D_onChunk(RAD,alb_vis)  ! Temporary
    IF (PRESENT(alb_nir_dir))     alb_nir_dir(:)      = dsl4jsb_var2D_onChunk(RAD,alb_nir)  ! Temporary
    IF (PRESENT(alb_nir_dif))     alb_nir_dif(:)      = dsl4jsb_var2D_onChunk(RAD,alb_nir)  ! Temporary

    IF (tile%lct%Contains(LAKE_TYPE)) THEN
      IF (PRESENT(t_wtr))           t_wtr(:)            = dsl4jsb_var2D_onChunk(SEB,t_wtr)        ! Temporary
      IF (PRESENT(t_ice))           t_ice(:)            = dsl4jsb_var2D_onChunk(SEB,t_ice)        ! Temporary
      IF (PRESENT(ice_fract_lake))  ice_fract_lake(:)   = 0._wp                                   ! Temporary
    ELSE
      IF (PRESENT(t_wtr))           t_wtr(:)            = 0._wp                                   ! Temporary
      IF (PRESENT(t_ice))           t_ice(:)            = 0._wp                                   ! Temporary
      IF (PRESENT(ice_fract_lake))  ice_fract_lake(:)   = 0._wp                                   ! Temporary
    END IF

  END SUBROUTINE update_land2atm

#endif
END MODULE mo_atmland_interface
