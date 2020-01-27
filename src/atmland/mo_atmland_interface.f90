!>
!! @brief Interface to run JSBACH for one time step
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_atmland_interface
#ifndef __NO_JSBACH__

  USE mo_exception,       ONLY: finish
  USE mo_kind,            ONLY: wp

  USE mo_jsb_control,        ONLY: jsbach_runs_standalone
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_lct_class,      ONLY: LAKE_TYPE
  !USE mo_jsb_config_class,   ONLY: t_jsb_config
  !USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options

  dsl4jsb_Use_processes A2L_, SEB_, TURB_, HYDRO_, RAD_, CARBON_
  dsl4jsb_Use_config(SEB_)
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
    REAL(wp), OPTIONAL, INTENT(in) ::                                                                                &
      & t_air(options%nc), q_air(options%nc), rain(options%nc), snow(options%nc),                                    &
      & wind_air(options%nc), wind_10m(options%nc),                                                                  &
      & lw_srf_down(options%nc), swvis_srf_down(options%nc), swnir_srf_down(options%nc), swpar_srf_down(options%nc), &
      & fract_par_diffuse(options%nc), press_srf(options%nc), drag_srf(options%nc),                                  &
      & t_acoef(options%nc), t_bcoef(options%nc), q_acoef(options%nc), q_bcoef(options%nc),                          &
      & pch(options%nc), cos_zenith_angle(options%nc),                                                               &
      & CO2_air(options%nc),                                                                                         &
      & drag_wtr(options%nc), drag_ice(options%nc),                                                                  &
      & t_acoef_wtr(options%nc), t_bcoef_wtr(options%nc), q_acoef_wtr(options%nc), q_bcoef_wtr(options%nc),          &
      & t_acoef_ice(options%nc), t_bcoef_ice(options%nc), q_acoef_ice(options%nc), q_bcoef_ice(options%nc)

    dsl4jsb_Def_memory(A2L_)

    INTEGER ::  iblk, ics, ice, nc, i
    LOGICAL ::  tile_contains_lake

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_atm2land'

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    ! make this logical accessible in the OpenACC code directly
    tile_contains_lake = tile%contains_lake

    IF (nc /= SIZE(t_air,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl4jsb_Get_memory(A2L_)

    !$acc parallel default(present)

    IF (PRESENT(t_air)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, t_air, ics+i-1) = t_air(i)
      END DO
    END IF
    IF (PRESENT(q_air)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, q_air, ics+i-1) = q_air(i)
      END DO
    END IF
    IF (PRESENT(rain)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, rain, ics+i-1)= rain(i)
      END DO
    END IF
    IF (PRESENT(snow)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, snow, ics+i-1) = snow(i)
      END DO
    END IF
    IF (PRESENT(wind_air)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, wind_air, ics+i-1) = wind_air(i)
      END DO
    END IF
    IF (PRESENT(wind_10m)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, wind_10m, ics+i-1) = wind_10m(i)
      END DO
    END IF
    IF (PRESENT(lw_srf_down)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, lw_srf_down, ics+i-1) = lw_srf_down(i)
      END DO
    END IF
    IF (PRESENT(swvis_srf_down)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, swvis_srf_down, ics+i-1) = swvis_srf_down(i)
      END DO
    END IF
    IF (PRESENT(swnir_srf_down)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, swnir_srf_down, ics+i-1) = swnir_srf_down(i)
      END DO
    END IF
    IF (PRESENT(swpar_srf_down)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, swpar_srf_down, ics+i-1) = swpar_srf_down(i)
      END DO
    END IF
    IF (PRESENT(fract_par_diffuse)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, fract_par_diffuse, ics+i-1)  = fract_par_diffuse(i)
      END DO
    END IF
    IF (PRESENT(press_srf)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, press_srf, ics+i-1) = press_srf(i)
      END DO
    END IF
    IF (PRESENT(drag_srf)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, drag_srf, ics+i-1) = drag_srf(i)
      END DO
    END IF
    IF (PRESENT(t_acoef)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, t_acoef, ics+i-1) = t_acoef(i)
      END DO
    END IF
    IF (PRESENT(t_bcoef)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, t_bcoef, ics+i-1) = t_bcoef(i)
      END DO
    END IF
    IF (PRESENT(q_acoef)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, q_acoef, ics+i-1) = q_acoef(i)
      END DO
    END IF
    IF (PRESENT(q_bcoef)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, q_bcoef, ics+i-1) = q_bcoef(i)
      END DO
    END IF
    IF (PRESENT(pch)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, pch, ics+i-1) = pch(i)
      END DO
    END IF
    IF (PRESENT(cos_zenith_angle)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, cos_zenith_angle, ics+i-1) = cos_zenith_angle(i)
      END DO
    END IF
    IF (PRESENT(CO2_air)) THEN
      !$acc loop gang vector
      DO i=1,nc
        dsl4jsb_var2D_indexed(A2L_, CO2_air, ics+i-1) = CO2_air(i)
        ! Convert CO2 mass mixing ratio [kg/kg] to particle mixing ratio [mol/mol]
        dsl4jsb_var2D_indexed(A2L_, CO2_air_mol, ics+i-1) = CO2_air(i) * molarMassDryAir / molarMassCO2
      END DO
    END IF

    IF (tile_contains_lake) THEN

      IF (PRESENT(drag_wtr)) THEN
      !$acc loop gang vector
        DO i=1,nc
          dsl4jsb_var2D_indexed(A2L_, drag_wtr, ics+i-1) = drag_wtr(i)
        END DO
      END IF
      IF (PRESENT(drag_ice)) THEN
      !$acc loop gang vector
        DO i=1,nc
          dsl4jsb_var2D_indexed(A2L_, drag_ice, ics+i-1) = drag_ice(i)
        END DO
      END IF
      IF (PRESENT(t_acoef_wtr)) THEN
      !$acc loop gang vector
        DO i=1,nc
          dsl4jsb_var2D_indexed(A2L_, t_acoef_wtr, ics+i-1) = t_acoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(t_bcoef_wtr)) THEN
      !$acc loop gang vector
        DO i=1,nc
          dsl4jsb_var2D_indexed(A2L_, t_bcoef_wtr, ics+i-1) = t_bcoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(q_acoef_wtr)) THEN
      !$acc loop gang vector
        DO i=1,nc
          dsl4jsb_var2D_indexed(A2L_, q_acoef_wtr, ics+i-1) = q_acoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(q_bcoef_wtr)) THEN
      !$acc loop gang vector
        DO i=1,nc
          dsl4jsb_var2D_indexed(A2L_, q_bcoef_wtr, ics+i-1) = q_bcoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(t_acoef_ice)) THEN
      !$acc loop gang vector
        DO i=1,nc
          dsl4jsb_var2D_indexed(A2L_, t_acoef_ice, ics+i-1) = t_acoef_ice(i)
        END DO
      END IF
      IF (PRESENT(t_bcoef_ice)) THEN
      !$acc loop gang vector
        DO i=1,nc
          dsl4jsb_var2D_indexed(A2L_, t_bcoef_ice, ics+i-1) = t_bcoef_ice(i)
        END DO
      END IF
      IF (PRESENT(q_acoef_ice)) THEN
      !$acc loop gang vector
        DO i=1,nc
          dsl4jsb_var2D_indexed(A2L_, q_acoef_ice, ics+i-1) = q_acoef_ice(i)
        END DO
      END IF
      IF (PRESENT(q_bcoef_ice)) THEN
      !$acc loop gang vector
        DO i=1,nc
          dsl4jsb_var2D_indexed(A2L_, q_bcoef_ice, ics+i-1) = q_bcoef_ice(i)
        END DO
      END IF

    END IF
    !$acc end parallel

  END SUBROUTINE update_atm2land

  SUBROUTINE update_land2atm(tile, options,                                          &
    & t_srf, t_eff_srf, qsat_srf, s_srf,                                             &
    & fact_q_air, fact_qsat_srf, evapopot,                                           &
    & evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap,                  &
    & rough_h_srf, rough_m_srf, q_snocpymlt,                                         &
    & alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif,                            &
    & CO2_flux,                                                                      &
    & t_lwtr, t_lice, qsat_lwtr, qsat_lice, s_lwtr, s_lice,                          &
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
      & t_lwtr(:), t_lice(:), qsat_lwtr(:), qsat_lice(:), s_lwtr(:), s_lice(:),      &
      & evapo_wtr(:), latent_hflx_wtr(:), sensible_hflx_wtr(:),                      &
      & evapo_ice(:), latent_hflx_ice(:), sensible_hflx_ice(:),                      &
      & ice_fract_lake(:),                                                           &
      & alb_vis_dir_wtr(:), alb_vis_dif_wtr(:), alb_nir_dir_wtr(:),                  &
      & alb_nir_dif_wtr(:), albedo_lwtr(:), albedo_lice(:),                          &
      & alb_vis_dir(:), alb_nir_dir(:), alb_vis_dif(:), alb_nir_dif(:),              &
      & CO2_flux(:)

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(CARBON_)

    TYPE(t_jsb_model), POINTER :: model

    INTEGER ::  iblk, ics, ice, nc, i, j
    REAL(wp) :: fact_lake(options%nc), fact_land(options%nc)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_land2atm'

    dsl4jsb_Real2D_onChunk :: &
      & t_filt_ptr, &
      & t_eff_ptr, &
      & qsat_star_ptr, &
      & s_star_ptr, &
      & fact_q_air_ptr, &
      & fact_qsat_srf_ptr, &
      & evapopot_ptr, &
      & evapotrans_lnd_ptr, &
      & latent_hflx_lnd_ptr, &
      & sensible_hflx_lnd_ptr, &
      & forc_hflx_ptr, &
      & heat_cap_old_ptr, &
      & rough_h_ptr, &
      & rough_m_ptr, &
      & q_snocpymlt_ptr, &
      & alb_vis_lnd_ptr, &
      & alb_nir_lnd_ptr, &
      & cflux_npp_2_atm_ta_ptr, &
      & cflux_soilresp_2_atm_ta_ptr, &
      & cflux_herb_2_atm_ta_ptr, &
      & t_lwtr_ptr, &
      & t_lice_ptr, &
      & qsat_lwtr_ptr, &
      & qsat_lice_ptr, &
      & s_lwtr_ptr, &
      & s_lice_ptr, &
      & evapo_wtr_ptr, &
      & latent_hflx_wtr_ptr, &
      & sensible_hflx_wtr_ptr, &
      & evapo_ice_ptr, &
      & latent_hflx_ice_ptr, &
      & fract_lice_ptr, &
      & albedo_lwtr_ptr, &
      & albedo_lice_ptr, &
      & sensible_hflx_ice_ptr, &
      & fract_ptr

    ! avoid compiler warnings about dummy arguments not being used
    IF (PRESENT(alb_nir_dif_wtr)) CONTINUE
    IF (PRESENT(alb_nir_dir_wtr)) CONTINUE
    IF (PRESENT(alb_vis_dif_wtr)) CONTINUE
    IF (PRESENT(alb_vis_dir_wtr)) CONTINUE

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (nc /= SIZE(t_srf,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
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
      DO j=1,SIZE(tile%lcts)
        IF (tile%lcts(j)%id == LAKE_TYPE) EXIT
      END DO

      !$acc data pcreate(fact_lake, fact_land, fract_ptr) present(tile%lcts)
      fract_ptr => tile%lcts(j)%fract(ics:ice,iblk)
      
      !$acc data present(fract_ptr)
      !$acc parallel loop
      DO i=1,nc
        fact_lake(i) = fract_ptr(i)
        fact_land(i) = 1._wp - fact_lake(i)
        IF (fact_lake(i) > 0._wp) THEN
          fact_lake(i) = 1._wp / fact_lake(i)
        ELSE
          fact_lake(i) = 1._wp
        END IF
        IF (fact_land(i) > 0._wp) THEN
          fact_land(i) = 1._wp / fact_land(i)
        ELSE
          fact_land(i) = 1._wp
        END IF
      END DO
      !$acc end data
     

      IF (PRESENT(t_srf)) THEN
        t_filt_ptr => dsl4jsb_var2D_onChunk(SEB_, t_filt)
        !$acc data present(t_srf, t_filt_ptr)
        !$acc parallel loop
        DO i=1,nc
          t_srf(i) = t_filt_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(t_eff_srf)) THEN
        t_eff_ptr => dsl4jsb_var2D_onChunk(SEB_,   t_eff)
        !$acc data present(t_eff_srf, t_eff_ptr)
        !$acc parallel loop
        DO i=1,nc
          t_eff_srf(i) = t_eff_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(qsat_srf)) THEN
        qsat_star_ptr => dsl4jsb_var2D_onChunk(SEB_,   qsat_star)
        !$acc data present(qsat_srf, qsat_star_ptr)
        !$acc parallel loop
        DO i=1,nc
          qsat_srf(i) = qsat_star_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(s_srf)) THEN
        s_star_ptr => dsl4jsb_var2D_onChunk(SEB_,   s_star)
        !$acc data present(s_srf, s_star_ptr)
        !$acc parallel loop
        DO i=1,nc
          s_srf(i) = s_star_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(fact_q_air)) THEN
        fact_q_air_ptr => dsl4jsb_var2D_onChunk(TURB_,  fact_q_air)
        !$acc data present(fact_q_air, fact_q_air_ptr)
        !$acc parallel loop
        DO i=1,nc
          fact_q_air(i) = fact_q_air_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(fact_qsat_srf)) THEN
        fact_qsat_srf_ptr => dsl4jsb_var2D_onChunk(TURB_,  fact_qsat_srf)
        !$acc data present(fact_qsat_srf, fact_qsat_srf_ptr)
        !$acc parallel loop
        DO i=1,nc
          fact_qsat_srf(i) = fact_qsat_srf_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(evapopot)) THEN
        evapopot_ptr => dsl4jsb_var2D_onChunk(HYDRO_, evapopot) ! For offline only
        !$acc data present(evapopot, evapopot_ptr)
        !$acc parallel loop
        DO i=1,nc
          evapopot(i) = evapopot_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(evapotrans)) THEN
        evapotrans_lnd_ptr => dsl4jsb_var2D_onChunk(HYDRO_, evapotrans_lnd)
        !$acc data present(fact_land, evapotrans)
        !$acc parallel loop
        DO i=1,nc
          evapotrans(i) = fact_land(i) * evapotrans_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(latent_hflx)) THEN
        latent_hflx_lnd_ptr => dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_lnd)
        !$acc data present(latent_hflx, latent_hflx_lnd_ptr, fact_land)
        !$acc parallel loop
        DO i=1,nc
          latent_hflx(i) = fact_land(i) * latent_hflx_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(sensible_hflx)) THEN
        sensible_hflx_lnd_ptr => dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_lnd)
        !$acc data present(sensible_hflx, fact_land, sensible_hflx_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          sensible_hflx(i) = fact_land(i) * sensible_hflx_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(grnd_hflx)) THEN
        forc_hflx_ptr => dsl4jsb_var2D_onChunk(SEB_,   forc_hflx)    ! TODO Not used
        !$acc data present(grnd_hflx, forc_hflx_ptr)
        !$acc parallel loop
        DO i=1,nc
          grnd_hflx(i) = forc_hflx_ptr(i)            
        END DO
        !$acc end data
      END IF
      IF (PRESENT(grnd_hcap)) THEN
        heat_cap_old_ptr => dsl4jsb_var2D_onChunk(SEB_,   heat_cap_old)
        !$acc data present(grnd_hcap, heat_cap_old_ptr)
        !$acc parallel loop
        DO i=1,nc
          grnd_hcap(i) = heat_cap_old_ptr(i)
        END DO
        !$acc end data
                                                                                                                   ! not computed
                                                                                                                   ! for lake
      END IF
      IF (PRESENT(rough_h_srf)) THEN
        rough_h_ptr => dsl4jsb_var2D_onChunk(TURB_,  rough_h)
        !$acc data present(rough_h_srf, rough_h_ptr)
        !$acc parallel loop
        DO i=1,nc
        rough_h_srf(i) = rough_h_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(rough_m_srf)) THEN
        rough_m_ptr => dsl4jsb_var2D_onChunk(TURB_,  rough_m)
        !$acc data present(rough_m_ptr, rough_m_srf)
        !$acc parallel loop
        DO i=1,nc
          rough_m_srf(i) = rough_m_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(q_snocpymlt)) THEN
        q_snocpymlt_ptr => dsl4jsb_var2D_onChunk(HYDRO_, q_snocpymlt)
        !$acc data present(q_snocpymlt, q_snocpymlt_ptr, fact_land)
        !$acc parallel loop
        DO i=1,nc
          q_snocpymlt(i) = fact_land(i) * q_snocpymlt_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_vis_dir)) THEN
        alb_vis_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_, alb_vis_lnd)  ! TODO
        !$acc data present(alb_vis_lnd_ptr, alb_vis_dir)
        !$acc parallel loop
        DO i=1,nc
          alb_vis_dir(i) = alb_vis_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_vis_dif)) THEN
        alb_vis_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_,   alb_vis_lnd)  ! TODO
        !$acc data present(alb_vis_dif, alb_vis_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          alb_vis_dif(i) = alb_vis_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_nir_dir)) THEN
        alb_nir_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_,   alb_nir_lnd)  ! TODO
        !$acc data present(alb_nir_dir, alb_nir_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          alb_nir_dir(i) = alb_nir_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_nir_dif)) THEN
        alb_nir_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_,   alb_nir_lnd)
        !$acc data present(alb_nir_dif, alb_nir_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          alb_nir_dif(i)      = alb_nir_lnd_ptr(i)
        END DO
        !$acc end data
      END IF


      IF (PRESENT(CO2_flux)) THEN
        IF (tile%Is_process_active(CARBON_)) THEN
          cflux_npp_2_atm_ta_ptr => dsl4jsb_var2D_onChunk(CARBON_,cflux_npp_2_atm_ta)
          cflux_soilresp_2_atm_ta_ptr => dsl4jsb_var2D_onChunk(CARBON_,cflux_soilresp_2_atm_ta)
          cflux_herb_2_atm_ta_ptr => dsl4jsb_var2D_onChunk(CARBON_,cflux_herb_2_atm_ta)
          !$acc data present(CO2_flux, fact_land, cflux_npp_2_atm_ta_ptr, cflux_soilresp_2_atm_ta_ptr, cflux_herb_2_atm_ta_ptr)
          !$acc parallel loop
          DO i=1,nc
            CO2_flux(i) = fact_land(i) *            &
              & (   cflux_npp_2_atm_ta_ptr(i)       &
              &   + cflux_soilresp_2_atm_ta_ptr(i)  &
              &   + cflux_herb_2_atm_ta_ptr(i)      &
              & )
          END DO
          !$acc end data
        ELSE
          !$acc kernels default(present)
          CO2_flux(:)         = 0._wp
          !$acc end kernels
        END IF
      END IF

      IF (PRESENT(t_lwtr)) THEN
        t_lwtr_ptr => dsl4jsb_var2D_onChunk(SEB_,   t_lwtr)
        !$acc data present(t_lwtr_ptr, t_lwtr)
        !$acc parallel loop
        DO i=1,nc
          t_lwtr(i) = t_lwtr_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(qsat_lwtr)) THEN
        qsat_lwtr_ptr => dsl4jsb_var2D_onChunk(SEB_,   qsat_lwtr)
        !$acc data present(qsat_lwtr, qsat_lwtr_ptr)
        !$acc parallel loop
        DO i=1,nc
          qsat_lwtr(i) = qsat_lwtr_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(s_lwtr)) THEN
        s_lwtr_ptr => dsl4jsb_var2D_onChunk(SEB_,   s_lwtr)
        !$acc data present(s_lwtr, s_lwtr_ptr)
        !$acc parallel loop
        DO i=1,nc
          s_lwtr(i) = s_lwtr_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(evapo_wtr)) THEN
        evapo_wtr_ptr => dsl4jsb_var2D_onChunk(HYDRO_, evapo_wtr)
        !$acc data present(evapo_wtr, evapo_wtr_ptr, fact_lake)
        !$acc parallel loop
        DO i=1,nc
          evapo_wtr(i) = fact_lake(i) * evapo_wtr_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(latent_hflx_wtr)) THEN
        latent_hflx_wtr_ptr => dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_wtr)
        !$acc data present(latent_hflx_wtr, latent_hflx_wtr_ptr, fact_lake)
        !$acc parallel loop
        DO i=1,nc
          latent_hflx_wtr(i)   = fact_lake(i) * latent_hflx_wtr_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(sensible_hflx_wtr)) THEN
        sensible_hflx_wtr_ptr => dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_wtr)
        !$acc data present(sensible_hflx_wtr, sensible_hflx_wtr_ptr, fact_lake)
        !$acc parallel loop
        DO i=1,nc
          sensible_hflx_wtr(i) = fact_lake(i) * sensible_hflx_wtr_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(albedo_lwtr)) THEN
        albedo_lwtr_ptr => dsl4jsb_var2D_onChunk(RAD_,   albedo_lwtr)
        !$acc data present(albedo_lwtr, albedo_lwtr_ptr)
        !$acc parallel loop
        DO i=1,nc
          albedo_lwtr(i) = albedo_lwtr_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(ice_fract_lake)) THEN
        fract_lice_ptr => dsl4jsb_var2D_onChunk(SEB_,   fract_lice)
        !$acc data present(ice_fract_lake, fract_lice_ptr)
        !$acc parallel loop
        DO i=1,nc
          ice_fract_lake(i) = fract_lice_ptr(i)
        END DO
        !$acc end data
      END IF

      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        IF (PRESENT(t_lice)) THEN
          t_lice_ptr => dsl4jsb_var2D_onChunk(SEB_,   t_lice)
          !$acc data present(t_lice, t_lice_ptr)
          !$acc parallel loop
          DO i=1,nc
            t_lice(i) = t_lice_ptr(i)
          END DO
          !$acc end data
        END IF
        IF (PRESENT(qsat_lice)) THEN
          qsat_lice_ptr => dsl4jsb_var2D_onChunk(SEB_,   qsat_lice)
          !$acc data present(qsat_lice, qsat_lice_ptr)
          !$acc parallel loop
          DO i=1,nc
            qsat_lice(i) = qsat_lice_ptr(i)
          END DO
          !$acc end data
        END IF
        IF (PRESENT(s_lice)) THEN
          s_lice_ptr => dsl4jsb_var2D_onChunk(SEB_,   s_lice)
          !$acc data present(s_lice, s_lice_ptr)
          !$acc parallel loop
          DO i=1,nc
            s_lice(i) = s_lice_ptr(i)
          END DO
          !$acc end data
        END IF
        IF (PRESENT(evapo_ice)) THEN
          evapo_ice_ptr => dsl4jsb_var2D_onChunk(HYDRO_, evapo_ice)
          !$acc data present(evapo_ice_ptr, evapo_ice, fact_lake)
          !$acc parallel loop
          DO i=1,nc
            evapo_ice(i) = fact_lake(i) * evapo_ice_ptr(i)
          END DO
          !$acc end data
        END IF
        IF (PRESENT(latent_hflx_ice)) THEN
          latent_hflx_ice_ptr => dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_ice)
          !$acc data present(latent_hflx_ice, latent_hflx_ice_ptr, fact_lake)
          !$acc parallel loop
          DO i=1,nc
            latent_hflx_ice(i) = fact_lake(i) * latent_hflx_ice_ptr(i)
          END DO
          !$acc end data
        END IF
        IF (PRESENT(sensible_hflx_ice)) THEN
          sensible_hflx_ice_ptr  => dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_ice)
          !$acc data present(sensible_hflx_ice_ptr, sensible_hflx_ice, fact_lake)
          !$acc parallel loop
          DO i=1,nc
            sensible_hflx_ice(i) = fact_lake(i) * sensible_hflx_ice_ptr(i)
          END DO
          !$acc end data
        END IF
        IF (PRESENT(albedo_lice)) THEN
          albedo_lice_ptr => dsl4jsb_var2D_onChunk(RAD_,   albedo_lice)
          !$acc data present(albedo_lice, albedo_lice_ptr)
          !$acc parallel loop
          DO i=1,nc
            albedo_lice(i) = albedo_lice_ptr(i)
          END DO
          !$acc end data
        END IF
      ELSE
        IF (PRESENT(t_lice)) THEN
          !$acc kernels default(present)
          t_lice(:)             = 273.15_wp
          !$acc end kernels
        END IF
        IF (PRESENT(qsat_lice)) THEN
          !$acc kernels default(present)
          qsat_lice(:)          = 0.0075_wp
          !$acc end kernels
        END IF
        IF (PRESENT(s_lice)) THEN
          !$acc kernels default(present)
          s_lice(:)             = 2.9E5_wp
          !$acc end kernels
        END IF
        IF (PRESENT(evapo_ice)) THEN
          !$acc kernels default(present)
          evapo_ice(:)          = 0._wp
          !$acc end kernels
        END IF
        IF (PRESENT(latent_hflx_ice)) THEN
          !$acc kernels default(present)
          latent_hflx_ice(:)    = 0._wp
          !$acc end kernels
        END IF
        IF (PRESENT(sensible_hflx_ice)) THEN
          !$acc kernels default(present)
          sensible_hflx_ice(:)  = 0._wp
          !$acc end kernels
        END IF
        IF (PRESENT(albedo_lice)) THEN
          !$acc kernels default(present)
          albedo_lice(:)        = 0._wp
          !$acc end kernels
        END IF
      END IF

      !$acc end data

    ELSE

      IF (PRESENT(t_srf)) THEN
        t_filt_ptr => dsl4jsb_var2D_onChunk(SEB_,   t_filt)
        !$acc data present(t_filt_ptr, t_srf)
        !$acc parallel loop
        DO i=1,nc
          t_srf(i) = t_filt_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(t_eff_srf)) THEN
        t_eff_ptr => dsl4jsb_var2D_onChunk(SEB_,   t_eff)
        !$acc data present(t_eff_srf, t_eff_ptr)
        !$acc parallel loop
        DO i=1,nc
          t_eff_srf(i) = t_eff_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(qsat_srf)) THEN
        qsat_star_ptr => dsl4jsb_var2D_onChunk(SEB_,   qsat_star)
        !$acc data present(qsat_srf, qsat_star_ptr)
        !$acc parallel loop
        DO i=1,nc
          qsat_srf(i)         = qsat_star_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(s_srf)) THEN
        s_star_ptr => dsl4jsb_var2D_onChunk(SEB_,   s_star)
        !$acc data present(s_star_ptr,s_srf)
        !$acc parallel loop
        DO i=1,nc
          s_srf(i)            = s_star_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(fact_q_air)) THEN
        fact_q_air_ptr => dsl4jsb_var2D_onChunk(TURB_,  fact_q_air)
        !$acc data present(fact_q_air, fact_q_air_ptr)
        !$acc parallel loop
        DO i=1,nc
          fact_q_air(i)       = fact_q_air_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(fact_qsat_srf)) THEN
        fact_qsat_srf_ptr => dsl4jsb_var2D_onChunk(TURB_,  fact_qsat_srf)
        !$acc data present(fact_qsat_srf, fact_qsat_srf_ptr)
        !$acc parallel loop
        DO i=1,nc
          fact_qsat_srf(i) = fact_qsat_srf_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(evapopot)) THEN
        evapopot_ptr => dsl4jsb_var2D_onChunk(HYDRO_, evapopot) ! For offline only
        !$acc data present(evapopot, evapopot_ptr)
        !$acc parallel loop
        DO i=1,nc
          evapopot(i) = evapopot_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(evapotrans)) THEN
        evapotrans_lnd_ptr => dsl4jsb_var2D_onChunk(HYDRO_, evapotrans_lnd)
        !$acc data present(evapotrans, evapotrans_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          evapotrans(i)= evapotrans_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(latent_hflx)) THEN
        latent_hflx_lnd_ptr => dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_lnd)
        !$acc data present(latent_hflx, latent_hflx_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          latent_hflx(i) = latent_hflx_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(sensible_hflx)) THEN
        sensible_hflx_lnd_ptr => dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_lnd)
        !$acc data present(sensible_hflx, sensible_hflx_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          sensible_hflx(i)    = sensible_hflx_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(grnd_hflx)) THEN
        forc_hflx_ptr => dsl4jsb_var2D_onChunk(SEB_,   forc_hflx)    ! TODO Not used in atm yet,
        !$acc data present(grnd_hflx, forc_hflx_ptr)
        !$acc parallel loop
        DO i=1,nc
          grnd_hflx(i) = forc_hflx_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(grnd_hcap)) THEN
        heat_cap_old_ptr => dsl4jsb_var2D_onChunk(SEB_,   heat_cap_old) ! and not computed for lake
        !$acc data present(grnd_hcap, heat_cap_old_ptr)
        !$acc parallel loop
        DO i=1,nc
          grnd_hcap(i)        = heat_cap_old_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(rough_h_srf)) THEN
        rough_h_ptr => dsl4jsb_var2D_onChunk(TURB_,  rough_h)
        !$acc data present(rough_h_srf, rough_h_ptr)
        !$acc parallel loop
        DO i=1,nc
          rough_h_srf(i) = rough_h_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(rough_m_srf)) THEN
        rough_m_ptr => dsl4jsb_var2D_onChunk(TURB_,  rough_m)
        !$acc data present(rough_m_srf, rough_m_ptr)
        !$acc parallel loop
        DO i=1,nc
          rough_m_srf(i)      = rough_m_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(q_snocpymlt)) THEN
        q_snocpymlt_ptr => dsl4jsb_var2D_onChunk(HYDRO_, q_snocpymlt)
        !$acc data present(q_snocpymlt, q_snocpymlt_ptr)
        !$acc parallel loop
        DO i=1,nc
          q_snocpymlt(i)      = q_snocpymlt(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_vis_dir)) THEN
        alb_vis_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_,   alb_vis_lnd)  ! TODO
       !$acc data present(alb_vis_lnd_ptr, alb_vis_dir)
       !$acc parallel loop
       DO i=1,nc
         alb_vis_dir(i) = alb_vis_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_vis_dif)) THEN
        alb_vis_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_,   alb_vis_lnd)  ! TODO
        !$acc data present(alb_vis_dif, alb_vis_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          alb_vis_dif(i) = alb_vis_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_nir_dir)) THEN
        alb_nir_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_,   alb_nir_lnd)  ! TODO
        !$acc data present(alb_nir_dir, alb_nir_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          alb_nir_dir(i) = alb_nir_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_nir_dif)) THEN
        alb_nir_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_,   alb_nir_lnd)  ! TODO
        !$acc data present(alb_nir_dif)
        !$acc parallel loop
        DO i=1,nc
          alb_nir_dif(i) = alb_nir_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(CO2_flux)) THEN
        IF (tile%Is_process_active(CARBON_)) THEN
          cflux_npp_2_atm_ta_ptr => dsl4jsb_var2D_onChunk(CARBON_,cflux_npp_2_atm_ta)
          cflux_soilresp_2_atm_ta_ptr => dsl4jsb_var2D_onChunk(CARBON_,cflux_soilresp_2_atm_ta)
          cflux_herb_2_atm_ta_ptr => dsl4jsb_var2D_onChunk(CARBON_,cflux_herb_2_atm_ta)
          !$acc data present(CO2_flux, cflux_npp_2_atm_ta_ptr, cflux_soilresp_2_atm_ta_ptr, cflux_herb_2_atm_ta_ptr)
          !$acc parallel loop
          DO i=1,nc
            CO2_flux(i) =                          &
              & (   cflux_npp_2_atm_ta_ptr(i)      &
              &   + cflux_soilresp_2_atm_ta_ptr(i) &
              &   + cflux_herb_2_atm_ta_ptr(i)     &
              & )
          END DO
          !$acc end data
        ELSE
          !$acc kernels default(present)
                                    CO2_flux(:)         = 0._wp
          !$acc end kernels
        END IF
      END IF

      IF (PRESENT(t_lwtr)) THEN 
        !$acc kernels default(present)
        t_lwtr(:)            = 0._wp                                   ! Temporary
        !$acc end kernels
      END IF
      IF (PRESENT(t_lice)) THEN 
        !$acc kernels default(present)
        t_lice(:)            = 0._wp                                   ! Temporary
        !$acc end kernels
      END IF
      IF (PRESENT(qsat_lwtr)) THEN 
        !$acc kernels default(present)
        qsat_lwtr(:)         = 0._wp
        !$acc end kernels
      END IF
      IF (PRESENT(qsat_lice)) THEN 
        !$acc kernels default(present)
        qsat_lice(:)         = 0._wp
        !$acc end kernels
      END IF
      IF (PRESENT(s_lwtr)) THEN 
        !$acc kernels default(present)
        s_lwtr(:)            = 0._wp                                   ! Temporary
        !$acc end kernels
      END IF
      IF (PRESENT(s_lice)) THEN 
        !$acc kernels default(present)
        s_lice(:)            = 0._wp                                   ! Temporary
        !$acc end kernels
      END IF
      IF (PRESENT(evapo_wtr)) THEN 
        !$acc kernels default(present)
        evapo_wtr(:)         = 0._wp
        !$acc end kernels
      END IF
      IF (PRESENT(latent_hflx_wtr)) THEN
        !$acc kernels default(present) 
        latent_hflx_wtr(:)   = 0._wp
        !$acc end kernels
      END IF
      IF (PRESENT(sensible_hflx_wtr)) THEN
        !$acc kernels default(present)
        sensible_hflx_wtr(:) = 0._wp
        !$acc end kernels
      END IF
      IF (PRESENT(evapo_ice)) THEN 
        !$acc kernels default(present)
        evapo_ice(:)         = 0._wp
        !$acc end kernels
      END IF
      IF (PRESENT(latent_hflx_ice)) THEN
        !$acc kernels default(present)
        latent_hflx_ice(:)   = 0._wp
        !$acc end kernels
      END IF
      IF (PRESENT(sensible_hflx_ice)) THEN
        !$acc kernels default(present)
        sensible_hflx_ice(:) = 0._wp
        !$acc end kernels
      END IF
      IF (PRESENT(ice_fract_lake)) THEN 
        !$acc kernels default(present)
        ice_fract_lake(:)    = 0._wp                                   ! Temporary
        !$acc end kernels
      END IF
      IF (PRESENT(albedo_lwtr)) THEN
        !$acc kernels default(present)
        albedo_lwtr(:)       = 0.07_wp
        !$acc end kernels
      END IF
      IF (PRESENT(albedo_lice)) THEN 
        !$acc kernels default(present)
        albedo_lice(:)       = 0.55_wp
        !$acc end kernels
      END IF
    END IF

  END SUBROUTINE update_land2atm

#endif
END MODULE mo_atmland_interface
