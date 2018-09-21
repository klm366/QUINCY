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

  USE mo_exception,       ONLY: finish
  USE mo_kind,            ONLY: wp

  USE mo_jsb_control,        ONLY: jsbach_runs_standalone
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_lct_class,      ONLY: LAKE_TYPE
  USE mo_jsb_memory_class,   ONLY: t_jsb_memory, t_jsb_memory_p
  USE mo_jsb_config_class,   ONLY: t_jsb_config
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options

  dsl4jsb_Use_processes A2L_, SEB_, TURB_, HYDRO_, RAD_, CARBON_
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(CARBON_)

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
    & fract_par_diffuse,       &
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
    & q_bcoef_ice             &
    & )

    USE mo_jsb_physical_constants, ONLY: molarMassDryAir, molarMassCO2

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    REAL(wp), OPTIONAL, INTENT(in) ::                                            &
      & t_air(:), q_air(:), rain(:), snow(:), wind_air(:), wind_10m(:),          &
      & lw_srf_down(:), swvis_srf_down(:), swnir_srf_down(:), swpar_srf_down(:), &
      & fract_par_diffuse(:), press_srf(:), drag_srf(:),                          &
      & t_acoef(:), t_bcoef(:), q_acoef(:), q_bcoef(:),                          &
      & pch(:), cos_zenith_angle(:),                                             &
      & CO2_air(:),                                                              &
      & drag_wtr(:), drag_ice(:),                                                &
      & t_acoef_wtr(:), t_bcoef_wtr(:), q_acoef_wtr(:), q_bcoef_wtr(:),          &
      & t_acoef_ice(:), t_bcoef_ice(:), q_acoef_ice(:), q_bcoef_ice(:)

    dsl4jsb_Def_memory(A2L_)

    INTEGER ::  iblk, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_atm2land'

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (ice-ics+1 /= SIZE(t_air,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl4jsb_Get_memory(A2L_)

    IF (PRESENT(t_air))             dsl4jsb_var2D_onChunk(A2L_,t_air)             = t_air             (:)
    IF (PRESENT(q_air))             dsl4jsb_var2D_onChunk(A2L_,q_air)             = q_air             (:)
    IF (PRESENT(rain))              dsl4jsb_var2D_onChunk(A2L_,rain)              = rain              (:)
    IF (PRESENT(snow))              dsl4jsb_var2D_onChunk(A2L_,snow)              = snow              (:)
    IF (PRESENT(wind_air))          dsl4jsb_var2D_onChunk(A2L_,wind_air)          = wind_air          (:)
    IF (PRESENT(wind_10m))          dsl4jsb_var2D_onChunk(A2L_,wind_10m)          = wind_10m          (:)
    IF (PRESENT(lw_srf_down))       dsl4jsb_var2D_onChunk(A2L_,lw_srf_down)       = lw_srf_down       (:)
    IF (PRESENT(swvis_srf_down))    dsl4jsb_var2D_onChunk(A2L_,swvis_srf_down)    = swvis_srf_down    (:)
    IF (PRESENT(swnir_srf_down))    dsl4jsb_var2D_onChunk(A2L_,swnir_srf_down)    = swnir_srf_down    (:)
    IF (PRESENT(swpar_srf_down))    dsl4jsb_var2D_onChunk(A2L_,swpar_srf_down)    = swpar_srf_down    (:)
    IF (PRESENT(fract_par_diffuse))  dsl4jsb_var2D_onChunk(A2L_,fract_par_diffuse)  = fract_par_diffuse  (:)
    IF (PRESENT(press_srf))         dsl4jsb_var2D_onChunk(A2L_,press_srf)         = press_srf         (:)
    IF (PRESENT(drag_srf))          dsl4jsb_var2D_onChunk(A2L_,drag_srf)          = drag_srf          (:)
    IF (PRESENT(t_acoef))           dsl4jsb_var2D_onChunk(A2L_,t_acoef)           = t_acoef           (:)
    IF (PRESENT(t_bcoef))           dsl4jsb_var2D_onChunk(A2L_,t_bcoef)           = t_bcoef           (:)
    IF (PRESENT(q_acoef))           dsl4jsb_var2D_onChunk(A2L_,q_acoef)           = q_acoef           (:)
    IF (PRESENT(q_bcoef))           dsl4jsb_var2D_onChunk(A2L_,q_bcoef)           = q_bcoef           (:)
    IF (PRESENT(pch))               dsl4jsb_var2D_onChunk(A2L_,pch)               = pch               (:)
    IF (PRESENT(cos_zenith_angle))  dsl4jsb_var2D_onChunk(A2L_,cos_zenith_angle)  = cos_zenith_angle  (:)
    IF (PRESENT(CO2_air)) THEN
      dsl4jsb_var2D_onChunk(A2L_, CO2_air) = CO2_air(:)
      ! Convert CO2 mass mixing ratio [kg/kg] to particle mixing ratio [mol/mol]
      dsl4jsb_var2D_onChunk(A2L_, CO2_air_mol) = CO2_air(:) * molarMassDryAir / molarMassCO2
    END IF

    IF (tile%contains_lake) THEN

      IF (PRESENT(drag_wtr))        dsl4jsb_var2D_onChunk(A2L_,drag_wtr)          = drag_wtr          (:)
      IF (PRESENT(drag_ice))        dsl4jsb_var2D_onChunk(A2L_,drag_ice)          = drag_ice          (:)
      IF (PRESENT(t_acoef_wtr))     dsl4jsb_var2D_onChunk(A2L_,t_acoef_wtr)       = t_acoef_wtr       (:)
      IF (PRESENT(t_bcoef_wtr))     dsl4jsb_var2D_onChunk(A2L_,t_bcoef_wtr)       = t_bcoef_wtr       (:)
      IF (PRESENT(q_acoef_wtr))     dsl4jsb_var2D_onChunk(A2L_,q_acoef_wtr)       = q_acoef_wtr       (:)
      IF (PRESENT(q_bcoef_wtr))     dsl4jsb_var2D_onChunk(A2L_,q_bcoef_wtr)       = q_bcoef_wtr       (:)
      IF (PRESENT(t_acoef_ice))     dsl4jsb_var2D_onChunk(A2L_,t_acoef_ice)       = t_acoef_ice       (:)
      IF (PRESENT(t_bcoef_ice))     dsl4jsb_var2D_onChunk(A2L_,t_bcoef_ice)       = t_bcoef_ice       (:)
      IF (PRESENT(q_acoef_ice))     dsl4jsb_var2D_onChunk(A2L_,q_acoef_ice)       = q_acoef_ice       (:)
      IF (PRESENT(q_bcoef_ice))     dsl4jsb_var2D_onChunk(A2L_,q_bcoef_ice)       = q_bcoef_ice       (:)

    END IF

  END SUBROUTINE update_atm2land

  SUBROUTINE update_land2atm(tile, options,                                          &
    & t_srf, t_eff_srf, qsat_srf, s_srf,                                             &
    & fact_q_air, fact_qsat_srf, evapopot,                                           &
    & evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap,                  &
    & rough_h_srf, rough_m_srf, q_snocpymlt,                                         &
    & alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif,                            &
    & CO2_flux,                                                                      &
    & t_lwtr, t_lice, qsat_lwtr, qsat_lice,                                          &
    & evapo_wtr, latent_hflx_wtr, sensible_hflx_wtr,                                 &
    & evapo_ice, latent_hflx_ice, sensible_hflx_ice,                                 &
    & ice_fract_lake,                                                                &
    & alb_vis_dir_wtr, alb_vis_dif_wtr, alb_nir_dir_wtr, alb_nir_dif_wtr,            &
    & albedo_lwtr, albedo_lice)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    REAL(wp), OPTIONAL, INTENT(out) ::                                               &
      & t_srf(:), t_eff_srf(:), qsat_srf(:), s_srf(:),                               &
      & fact_q_air(:), fact_qsat_srf(:), evapopot(:),                                &
      & evapotrans(:), latent_hflx(:), sensible_hflx(:), grnd_hflx(:), grnd_hcap(:), &
      & rough_h_srf(:), rough_m_srf(:), q_snocpymlt(:),                              &
      & t_lwtr(:), t_lice(:), qsat_lwtr(:), qsat_lice(:),                            &
      & evapo_wtr(:), latent_hflx_wtr(:), sensible_hflx_wtr(:),                      &
      & evapo_ice(:), latent_hflx_ice(:), sensible_hflx_ice(:),                      &
      & ice_fract_lake(:),                                                           &
      & alb_vis_dir_wtr(:), alb_vis_dif_wtr(:), alb_nir_dir_wtr(:),                  &
      & alb_nir_dif_wtr(:), albedo_lwtr(:), albedo_lice(:),                          &
      & alb_vis_dir(:), alb_nir_dir(:), alb_vis_dif(:), alb_nir_dif(:),              &
      & CO2_flux(:)

    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(CARBON_)

    INTEGER ::  iblk, ics, ice, i
    REAL(wp) :: fact_lake(options%nc), fact_land(options%nc)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_land2atm'

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (ice-ics+1 /= SIZE(t_srf,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(RAD_)
    IF (tile%Is_process_active(CARBON_)) THEN
      dsl4jsb_Get_memory(CARBON_)
    END IF

    IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone()) THEN
      ! With lakes, re-scale _lnd, _lice and _lwtr fluxes given back to the atmosphere so that they are relative to the whole grid
      ! box. This is necessary because, in the atmosphere, they are weighted with the fractions of lake/nolake in the whole
      ! grid box, while on the root tile of JSBACH they are relative to the land fraction of the grid box including lakes.
      ! Note that, if lakes are considered, the lake tile must be a direct child of the root tile.
      DO i=1,SIZE(tile%lcts)
        IF (tile%lcts(i)%id == LAKE_TYPE) EXIT
      END DO
      fact_lake(:) = tile%lcts(i)%fract(ics:ice,iblk)
      fact_land(:) = 1._wp - fact_lake(:)
      WHERE (fact_lake(:) > 0._wp)
        fact_lake(:) = 1._wp / fact_lake(:)
      ELSE WHERE
        fact_lake(:) = 1._wp
      END WHERE
      WHERE (fact_land(:) > 0._wp)
        fact_land(:) = 1._wp / fact_land(:)
      ELSE WHERE
        fact_land(:) = 1._wp
      END WHERE

      IF (PRESENT(t_srf))         t_srf(:)            =                dsl4jsb_var2D_onChunk(SEB_,   t_filt)
      IF (PRESENT(t_eff_srf))     t_eff_srf(:)        =                dsl4jsb_var2D_onChunk(SEB_,   t_eff)
      IF (PRESENT(qsat_srf))      qsat_srf(:)         =                dsl4jsb_var2D_onChunk(SEB_,   qsat_star)
      IF (PRESENT(s_srf))         s_srf(:)            =                dsl4jsb_var2D_onChunk(SEB_,   s_star)
      IF (PRESENT(fact_q_air))    fact_q_air(:)       =                dsl4jsb_var2D_onChunk(TURB_,  fact_q_air)
      IF (PRESENT(fact_qsat_srf)) fact_qsat_srf(:)    =                dsl4jsb_var2D_onChunk(TURB_,  fact_qsat_srf)
      IF (PRESENT(evapopot))      evapopot(:)         =                dsl4jsb_var2D_onChunk(HYDRO_, evapopot) ! For offline only
      IF (PRESENT(evapotrans))    evapotrans(:)       = fact_land(:) * dsl4jsb_var2D_onChunk(HYDRO_, evapotrans_lnd)
      IF (PRESENT(latent_hflx))   latent_hflx(:)      = fact_land(:) * dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_lnd)
      IF (PRESENT(sensible_hflx)) sensible_hflx(:)    = fact_land(:) * dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_lnd)
      IF (PRESENT(grnd_hflx))     grnd_hflx(:)        =                dsl4jsb_var2D_onChunk(SEB_,   forc_hflx)    ! TODO Not used in atm yet, and
      IF (PRESENT(grnd_hcap))     grnd_hcap(:)        =                dsl4jsb_var2D_onChunk(SEB_,   heat_cap_old) ! not computed for lake
      IF (PRESENT(rough_h_srf))   rough_h_srf(:)      =                dsl4jsb_var2D_onChunk(TURB_,  rough_h)
      IF (PRESENT(rough_m_srf))   rough_m_srf(:)      =                dsl4jsb_var2D_onChunk(TURB_,  rough_m)
      IF (PRESENT(q_snocpymlt))   q_snocpymlt(:)      = fact_land(:) * dsl4jsb_var2D_onChunk(HYDRO_, q_snocpymlt)
      IF (PRESENT(alb_vis_dir))   alb_vis_dir(:)      =                dsl4jsb_var2D_onChunk(RAD_,   alb_vis_lnd)  ! TODO
      IF (PRESENT(alb_vis_dif))   alb_vis_dif(:)      =                dsl4jsb_var2D_onChunk(RAD_,   alb_vis_lnd)  ! TODO
      IF (PRESENT(alb_nir_dir))   alb_nir_dir(:)      =                dsl4jsb_var2D_onChunk(RAD_,   alb_nir_lnd)  ! TODO
      IF (PRESENT(alb_nir_dif))   alb_nir_dif(:)      =                dsl4jsb_var2D_onChunk(RAD_,   alb_nir_lnd)  ! TODO
      IF (PRESENT(CO2_flux)) THEN
        IF (tile%Is_process_active(CARBON_)) THEN
                                  CO2_flux(:)         = fact_land(:) *                              &
          &                             (   dsl4jsb_var2D_onChunk(CARBON_,cflux_npp_2_atm_ta)       &
          &                               + dsl4jsb_var2D_onChunk(CARBON_,cflux_soilresp_2_atm_ta)  &
          &                               + dsl4jsb_var2D_onChunk(CARBON_,cflux_herb_2_atm_ta) &
          &                             )
        ELSE
                                  CO2_flux(:)         = 0._wp
        END IF
      END IF

      IF (PRESENT(t_lwtr))            t_lwtr(:)            =                dsl4jsb_var2D_onChunk(SEB_,   t_lwtr)
      IF (PRESENT(t_lice))            t_lice(:)            =                dsl4jsb_var2D_onChunk(SEB_,   t_lice)
      IF (PRESENT(qsat_lwtr))         qsat_lwtr(:)         =                dsl4jsb_var2D_onChunk(SEB_,   qsat_lwtr)
      IF (PRESENT(qsat_lice))         qsat_lice(:)         =                dsl4jsb_var2D_onChunk(SEB_,   qsat_lice)
      IF (PRESENT(evapo_wtr))         evapo_wtr(:)         = fact_lake(:) * dsl4jsb_var2D_onChunk(HYDRO_, evapo_wtr)
      IF (PRESENT(latent_hflx_wtr))   latent_hflx_wtr(:)   = fact_lake(:) * dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_wtr)
      IF (PRESENT(sensible_hflx_wtr)) sensible_hflx_wtr(:) = fact_lake(:) * dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_wtr)
      IF (PRESENT(evapo_ice))         evapo_ice(:)         = fact_lake(:) * dsl4jsb_var2D_onChunk(HYDRO_, evapo_ice)
      IF (PRESENT(latent_hflx_ice))   latent_hflx_ice(:)   = fact_lake(:) * dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_ice)
      IF (PRESENT(sensible_hflx_ice)) sensible_hflx_ice(:) = fact_lake(:) * dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_ice)
      IF (PRESENT(ice_fract_lake))    ice_fract_lake(:)    = fact_lake(:) * dsl4jsb_var2D_onChunk(SEB_,   fract_lice)
      IF (PRESENT(albedo_lwtr))       albedo_lwtr(:)       =                dsl4jsb_var2D_onChunk(RAD_,   albedo_lwtr)
      IF (PRESENT(albedo_lice))       albedo_lice(:)       =                dsl4jsb_var2D_onChunk(RAD_,   albedo_lice)

    ELSE

      IF (PRESENT(t_srf))           t_srf(:)            = dsl4jsb_var2D_onChunk(SEB_,   t_filt)
      IF (PRESENT(t_eff_srf))       t_eff_srf(:)        = dsl4jsb_var2D_onChunk(SEB_,   t_eff)
      IF (PRESENT(qsat_srf))        qsat_srf(:)         = dsl4jsb_var2D_onChunk(SEB_,   qsat_star)
      IF (PRESENT(s_srf))           s_srf(:)            = dsl4jsb_var2D_onChunk(SEB_,   s_star)
      IF (PRESENT(fact_q_air))      fact_q_air(:)       = dsl4jsb_var2D_onChunk(TURB_,  fact_q_air)
      IF (PRESENT(fact_qsat_srf))   fact_qsat_srf(:)    = dsl4jsb_var2D_onChunk(TURB_,  fact_qsat_srf)
      IF (PRESENT(evapopot))        evapopot(:)         = dsl4jsb_var2D_onChunk(HYDRO_, evapopot) ! For offline only
      IF (PRESENT(evapotrans))      evapotrans(:)       = dsl4jsb_var2D_onChunk(HYDRO_, evapotrans_lnd)
      IF (PRESENT(latent_hflx))     latent_hflx(:)      = dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_lnd)
      IF (PRESENT(sensible_hflx))   sensible_hflx(:)    = dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_lnd)
      IF (PRESENT(grnd_hflx))       grnd_hflx(:)        = dsl4jsb_var2D_onChunk(SEB_,   forc_hflx)    ! TODO Not used in atm yet, and
      IF (PRESENT(grnd_hcap))       grnd_hcap(:)        = dsl4jsb_var2D_onChunk(SEB_,   heat_cap_old) ! not computed for lake
      IF (PRESENT(rough_h_srf))     rough_h_srf(:)      = dsl4jsb_var2D_onChunk(TURB_,  rough_h)
      IF (PRESENT(rough_m_srf))     rough_m_srf(:)      = dsl4jsb_var2D_onChunk(TURB_,  rough_m)
      IF (PRESENT(q_snocpymlt))     q_snocpymlt(:)      = dsl4jsb_var2D_onChunk(HYDRO_, q_snocpymlt)
      IF (PRESENT(alb_vis_dir))     alb_vis_dir(:)      = dsl4jsb_var2D_onChunk(RAD_,   alb_vis_lnd)  ! TODO
      IF (PRESENT(alb_vis_dif))     alb_vis_dif(:)      = dsl4jsb_var2D_onChunk(RAD_,   alb_vis_lnd)  ! TODO
      IF (PRESENT(alb_nir_dir))     alb_nir_dir(:)      = dsl4jsb_var2D_onChunk(RAD_,   alb_nir_lnd)  ! TODO
      IF (PRESENT(alb_nir_dif))     alb_nir_dif(:)      = dsl4jsb_var2D_onChunk(RAD_,   alb_nir_lnd)  ! TODO
      IF (PRESENT(CO2_flux)) THEN
        IF (tile%Is_process_active(CARBON_)) THEN
                                  CO2_flux(:)           =                                             &
          &                               (   dsl4jsb_var2D_onChunk(CARBON_,cflux_npp_2_atm_ta)       &
          &                                 + dsl4jsb_var2D_onChunk(CARBON_,cflux_soilresp_2_atm_ta)  &
          &                                 + dsl4jsb_var2D_onChunk(CARBON_,cflux_herb_2_atm_ta) &
          &                               )
        ELSE
                                    CO2_flux(:)         = 0._wp
        END IF
      END IF

      IF (PRESENT(t_lwtr))            t_lwtr(:)            = 0._wp                                   ! Temporary
      IF (PRESENT(t_lice))            t_lice(:)            = 0._wp                                   ! Temporary
      IF (PRESENT(qsat_lwtr))         qsat_lwtr(:)         = 0._wp
      IF (PRESENT(qsat_lice))         qsat_lice(:)         = 0._wp
      IF (PRESENT(evapo_wtr))         evapo_wtr(:)         = 0._wp
      IF (PRESENT(latent_hflx_wtr))   latent_hflx_wtr(:)   = 0._wp
      IF (PRESENT(sensible_hflx_wtr)) sensible_hflx_wtr(:) = 0._wp
      IF (PRESENT(evapo_ice))         evapo_ice(:)         = 0._wp
      IF (PRESENT(latent_hflx_ice))   latent_hflx_ice(:)   = 0._wp
      IF (PRESENT(sensible_hflx_ice)) sensible_hflx_ice(:) = 0._wp
      IF (PRESENT(ice_fract_lake))    ice_fract_lake(:)    = 0._wp                                   ! Temporary
      IF (PRESENT(albedo_lwtr))       albedo_lwtr(:)       = 0.07_wp
      IF (PRESENT(albedo_lice))       albedo_lice(:)       = 0.55_wp
    END IF

  END SUBROUTINE update_land2atm

#endif
END MODULE mo_atmland_interface
