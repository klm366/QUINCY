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
  USE mo_jsb_config_class,   ONLY: t_jsb_config
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
  USE mo_jsb_lct_class,      ONLY: LAKE_TYPE, Contains_lct

  dsl4jsb_Use_processes ATM2LAND, SEB, TURB, SOIL, HYDRO, RAD
  dsl4jsb_Use_memory(ATM2LAND)
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

  SUBROUTINE update_atm2land( &
    & tile, options,          &
    & t_air,                  &
    & q_air,                  &
    & rain,                   &
    & snow,                   &
    & wind_air,               & 
    & wind_10m,               &
    & lw_srf_down,            &
    & swvis_srf_down,         &
    & swnir_srf_down,         &
    & swpar_srf_down,         &
    & frac_par_diffuse,       &
    & press_srf,              &
    & drag_srf,               &
    & t_acoef,                &
    & t_bcoef,                &
    & q_acoef,                &
    & q_bcoef,                &
    & pch,                    &
    & cos_zenith_angle,       &
    & CO2_air,                &
    & drag_wtr,               &
    & drag_ice,               &
    & t_acoef_wtr,            &
    & t_bcoef_wtr,            &
    & q_acoef_wtr,            &
    & q_bcoef_wtr,            &
    & t_acoef_ice,            &
    & t_bcoef_ice,            &
    & q_acoef_ice,            &
    & q_bcoef_ice,            &
    & albedo_lwtr             &
    & )

    USE mo_jsb_physical_constants, ONLY: molarMassDryAir, molarMassCO2

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    REAL(wp), OPTIONAL, INTENT(in) ::                                            &
      & t_air(:), q_air(:), rain(:), snow(:), wind_air(:), wind_10m(:),          &
      & lw_srf_down(:), swvis_srf_down(:), swnir_srf_down(:), swpar_srf_down(:), &
      & frac_par_diffuse(:), press_srf(:), drag_srf(:),                                               &
      & t_acoef(:), t_bcoef(:), q_acoef(:), q_bcoef(:),                          &
      & pch(:), cos_zenith_angle(:),                                             &
      & CO2_air(:),                                                              &
      & drag_wtr(:), drag_ice(:),                                                &
      & t_acoef_wtr(:), t_bcoef_wtr(:), q_acoef_wtr(:), q_bcoef_wtr(:),          &
      & t_acoef_ice(:), t_bcoef_ice(:), q_acoef_ice(:), q_bcoef_ice(:),          &
      & albedo_lwtr(:)

    dsl4jsb_Def_memory(ATM2LAND)

    INTEGER ::  iblk, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_atm2land'

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (ice-ics+1 /= SIZE(t_air,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl4jsb_Get_memory(ATM2LAND)

    IF (PRESENT(t_air))             dsl4jsb_var2D_onChunk(ATM2LAND,t_air)             = t_air             (:)
    IF (PRESENT(q_air))             dsl4jsb_var2D_onChunk(ATM2LAND,q_air)             = q_air             (:)
    IF (PRESENT(rain))              dsl4jsb_var2D_onChunk(ATM2LAND,rain)              = rain              (:)
    IF (PRESENT(snow))              dsl4jsb_var2D_onChunk(ATM2LAND,snow)              = snow              (:)
    IF (PRESENT(wind_air))          dsl4jsb_var2D_onChunk(ATM2LAND,wind_air)          = wind_air          (:)
    IF (PRESENT(wind_10m))          dsl4jsb_var2D_onChunk(ATM2LAND,wind_10m)          = wind_10m          (:)
    IF (PRESENT(lw_srf_down))       dsl4jsb_var2D_onChunk(ATM2LAND,lw_srf_down)       = lw_srf_down       (:)
    IF (PRESENT(swvis_srf_down))    dsl4jsb_var2D_onChunk(ATM2LAND,swvis_srf_down)    = swvis_srf_down    (:)
    IF (PRESENT(swnir_srf_down))    dsl4jsb_var2D_onChunk(ATM2LAND,swnir_srf_down)    = swnir_srf_down    (:)
    IF (PRESENT(swpar_srf_down))    dsl4jsb_var2D_onChunk(ATM2LAND,swpar_srf_down)    = swpar_srf_down    (:)
    IF (PRESENT(frac_par_diffuse))  dsl4jsb_var2D_onChunk(ATM2LAND,frac_par_diffuse)  = frac_par_diffuse  (:)
    IF (PRESENT(press_srf))         dsl4jsb_var2D_onChunk(ATM2LAND,press_srf)         = press_srf         (:)
    IF (PRESENT(drag_srf))          dsl4jsb_var2D_onChunk(ATM2LAND,drag_srf)          = drag_srf          (:)
    IF (PRESENT(t_acoef))           dsl4jsb_var2D_onChunk(ATM2LAND,t_acoef)           = t_acoef           (:)
    IF (PRESENT(t_bcoef))           dsl4jsb_var2D_onChunk(ATM2LAND,t_bcoef)           = t_bcoef           (:)
    IF (PRESENT(q_acoef))           dsl4jsb_var2D_onChunk(ATM2LAND,q_acoef)           = q_acoef           (:)
    IF (PRESENT(q_bcoef))           dsl4jsb_var2D_onChunk(ATM2LAND,q_bcoef)           = q_bcoef           (:)
    IF (PRESENT(pch))               dsl4jsb_var2D_onChunk(ATM2LAND,pch)               = pch               (:)
    IF (PRESENT(cos_zenith_angle))  dsl4jsb_var2D_onChunk(ATM2LAND,cos_zenith_angle)  = cos_zenith_angle  (:)
    IF (PRESENT(CO2_air)) THEN
      dsl4jsb_var2D_onChunk(ATM2LAND, CO2_air)               = CO2_air(:)
      ! Convert CO2 mass mixing ratio [kg/kg] to particle mixing ratio [mol/mol]
      dsl4jsb_var2D_onChunk(ATM2LAND, CO2_air_mol) = CO2_air(:) * molarMassDryAir / molarMassCO2
    END IF

    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN

      IF (PRESENT(drag_wtr))        dsl4jsb_var2D_onChunk(ATM2LAND,drag_wtr)          = drag_wtr          (:)
      IF (PRESENT(drag_ice))        dsl4jsb_var2D_onChunk(ATM2LAND,drag_ice)          = drag_ice          (:)
      IF (PRESENT(t_acoef_wtr))     dsl4jsb_var2D_onChunk(ATM2LAND,t_acoef_wtr)       = t_acoef_wtr       (:)
      IF (PRESENT(t_bcoef_wtr))     dsl4jsb_var2D_onChunk(ATM2LAND,t_bcoef_wtr)       = t_bcoef_wtr       (:)
      IF (PRESENT(q_acoef_wtr))     dsl4jsb_var2D_onChunk(ATM2LAND,q_acoef_wtr)       = q_acoef_wtr       (:)
      IF (PRESENT(q_bcoef_wtr))     dsl4jsb_var2D_onChunk(ATM2LAND,q_bcoef_wtr)       = q_bcoef_wtr       (:)
      IF (PRESENT(t_acoef_ice))     dsl4jsb_var2D_onChunk(ATM2LAND,t_acoef_ice)       = t_acoef_ice       (:)
      IF (PRESENT(t_bcoef_ice))     dsl4jsb_var2D_onChunk(ATM2LAND,t_bcoef_ice)       = t_bcoef_ice       (:)
      IF (PRESENT(q_acoef_ice))     dsl4jsb_var2D_onChunk(ATM2LAND,q_acoef_ice)       = q_acoef_ice       (:)
      IF (PRESENT(q_bcoef_ice))     dsl4jsb_var2D_onChunk(ATM2LAND,q_bcoef_ice)       = q_bcoef_ice       (:)
      IF (PRESENT(albedo_lwtr))     dsl4jsb_var2D_onChunk(ATM2LAND,albedo_lwtr)       = albedo_lwtr       (:)

    END IF

  END SUBROUTINE update_atm2land

  SUBROUTINE update_land2atm(tile, options,                                          &
    & t_srf, t_eff_srf, qsat_srf, s_srf,                                             &
    & fact_q_air, fact_qsat_srf,                                                     &
    & evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap,                  &
    & rough_h_srf, rough_m_srf, tte_corr,                                            &
    & alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif,                            &
    & t_lwtr, t_lice,                                                                &
    & evapo_wtr, latent_hflx_wtr, sensible_hflx_wtr,                                 &
    & evapo_ice, latent_hflx_ice, sensible_hflx_ice,                                 &
    & ice_fract_lake,                                                                &
    & alb_vis_dir_wtr, alb_vis_dif_wtr, alb_nir_dir_wtr, alb_nir_dif_wtr, albedo_lice)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    REAL(wp), OPTIONAL, INTENT(out) ::                                               &
      & t_srf(:), t_eff_srf(:), qsat_srf(:), s_srf(:),                               &
      & fact_q_air(:), fact_qsat_srf(:),                                             &
      & evapotrans(:), latent_hflx(:), sensible_hflx(:), grnd_hflx(:), grnd_hcap(:), &
      & rough_h_srf(:), rough_m_srf(:), tte_corr(:),                                 &
      & t_lwtr(:), t_lice(:),                                                        &
      & evapo_wtr(:), latent_hflx_wtr(:), sensible_hflx_wtr(:),                      &
      & evapo_ice(:), latent_hflx_ice(:), sensible_hflx_ice(:),                      &
      & ice_fract_lake(:),                                                           &
      & alb_vis_dir_wtr(:), alb_vis_dif_wtr(:), alb_nir_dir_wtr(:),                  &
      & alb_nir_dif_wtr(:), albedo_lice(:),                                             &
      & alb_vis_dir(:), alb_nir_dir(:), alb_vis_dif(:), alb_nir_dif(:)

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

    IF (PRESENT(t_srf))           t_srf(:)            = dsl4jsb_var2D_onChunk(SEB,   t_filt)
    IF (PRESENT(t_eff_srf))       t_eff_srf(:)        = dsl4jsb_var2D_onChunk(SEB,   t_eff)
    IF (PRESENT(qsat_srf))        qsat_srf(:)         = dsl4jsb_var2D_onChunk(SEB,   qsat_star)
    IF (PRESENT(s_srf))           s_srf(:)            = dsl4jsb_var2D_onChunk(SEB,   s_star)
    IF (PRESENT(fact_q_air))      fact_q_air(:)       = dsl4jsb_var2D_onChunk(TURB,  fact_q_air)
    IF (PRESENT(fact_qsat_srf))   fact_qsat_srf(:)    = dsl4jsb_var2D_onChunk(TURB,  fact_qsat_srf)
    IF (PRESENT(evapotrans))      evapotrans(:)       = dsl4jsb_var2D_onChunk(HYDRO, evapo_lnd)
    IF (PRESENT(latent_hflx))     latent_hflx(:)      = dsl4jsb_var2D_onChunk(SEB,   latent_hflx)
    IF (PRESENT(sensible_hflx))   sensible_hflx(:)    = dsl4jsb_var2D_onChunk(SEB,   sensible_hflx)
    IF (PRESENT(grnd_hflx))       grnd_hflx(:)        = dsl4jsb_var2D_onChunk(SEB,   forc_hflx)
    IF (PRESENT(grnd_hcap))       grnd_hcap(:)        = dsl4jsb_var2D_onChunk(SEB,   heat_cap_old)
    IF (PRESENT(rough_h_srf))     rough_h_srf(:)      = dsl4jsb_var2D_onChunk(TURB,  rough_h)
    IF (PRESENT(rough_m_srf))     rough_m_srf(:)      = dsl4jsb_var2D_onChunk(TURB,  rough_m)
    IF (PRESENT(tte_corr))        tte_corr(:)         = dsl4jsb_var2D_onChunk(HYDRO, tte_corr)
    IF (PRESENT(alb_vis_dir))     alb_vis_dir(:)      = dsl4jsb_var2D_onChunk(RAD,   alb_vis)  ! Temporary
    IF (PRESENT(alb_vis_dif))     alb_vis_dif(:)      = dsl4jsb_var2D_onChunk(RAD,   alb_vis)  ! Temporary
    IF (PRESENT(alb_nir_dir))     alb_nir_dir(:)      = dsl4jsb_var2D_onChunk(RAD,   alb_nir)  ! Temporary
    IF (PRESENT(alb_nir_dif))     alb_nir_dif(:)      = dsl4jsb_var2D_onChunk(RAD,   alb_nir)  ! Temporary

    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN
      IF (PRESENT(t_lwtr))            t_lwtr(:)            = dsl4jsb_var2D_onChunk(SEB,   t_lwtr)
      IF (PRESENT(t_lice))            t_lice(:)            = dsl4jsb_var2D_onChunk(SEB,   t_lice)
      IF (PRESENT(evapo_wtr))         evapo_wtr(:)         = dsl4jsb_var2D_onChunk(HYDRO, evapo_wtr)
      IF (PRESENT(latent_hflx_wtr))   latent_hflx_wtr(:)   = dsl4jsb_var2D_onChunk(SEB,   latent_hflx_wtr)
      IF (PRESENT(sensible_hflx_wtr)) sensible_hflx_wtr(:) = dsl4jsb_var2D_onChunk(SEB,   sensible_hflx_wtr)
      IF (PRESENT(evapo_ice))         evapo_ice(:)         = dsl4jsb_var2D_onChunk(HYDRO, evapo_ice)
      IF (PRESENT(latent_hflx_ice))   latent_hflx_ice(:)   = dsl4jsb_var2D_onChunk(SEB,   latent_hflx_ice)
      IF (PRESENT(sensible_hflx_ice)) sensible_hflx_ice(:) = dsl4jsb_var2D_onChunk(SEB,   sensible_hflx_ice)
      IF (PRESENT(ice_fract_lake))    ice_fract_lake(:)    = dsl4jsb_var2D_onChunk(SEB,   fract_lice)
      IF (PRESENT(albedo_lice))       albedo_lice(:)       = dsl4jsb_var2D_onChunk(RAD,   albedo_lice)
    ELSE
      IF (PRESENT(t_lwtr))            t_lwtr(:)             = 0._wp                                   ! Temporary
      IF (PRESENT(t_lice))            t_lice(:)             = 0._wp                                   ! Temporary
      IF (PRESENT(evapo_wtr))         evapo_wtr(:)         = 0._wp
      IF (PRESENT(latent_hflx_wtr))   latent_hflx_wtr(:)   = 0._wp
      IF (PRESENT(sensible_hflx_wtr)) sensible_hflx_wtr(:) = 0._wp
      IF (PRESENT(evapo_ice))         evapo_ice(:)         = 0._wp
      IF (PRESENT(latent_hflx_ice))   latent_hflx_ice(:)   = 0._wp
      IF (PRESENT(sensible_hflx_ice)) sensible_hflx_ice(:) = 0._wp
      IF (PRESENT(ice_fract_lake))    ice_fract_lake(:)    = 0._wp                                   ! Temporary
      IF (PRESENT(albedo_lice))       albedo_lice(:)       = 0.55_wp
    END IF

  END SUBROUTINE update_land2atm

#endif
END MODULE mo_atmland_interface
