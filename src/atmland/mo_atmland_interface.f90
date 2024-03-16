!> Interface to run JSBACH for one time step
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_atmland_interface
#ifndef __NO_JSBACH__

  USE mo_exception,       ONLY: message, finish
  USE mo_kind,            ONLY: wp

  USE mo_jsb_control,        ONLY: jsbach_runs_standalone, debug_on
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
    & t_air_in,                  &
    & q_air_in,                  &
    & press_air_in,              &
    & rain_in,                   &
    & snow_in,                   &
    & wind_air_in,               &
    & wind_10m_in,               &
    & lw_srf_down_in,            &
    & swvis_srf_down_in,         &
    & swnir_srf_down_in,         &
    & swpar_srf_down_in,         &
    & fract_par_diffuse_in,      &
    & dz_srf_in,                 &
    & press_srf_in,              &
    & rho_srf_in,                &
    & drag_srf_in,               &
    & t_acoef_in,                &
    & t_bcoef_in,                &
    & q_acoef_in,                &
    & q_bcoef_in,                &
    & pch_in,                    &
    & cos_zenith_angle_in,       &
    & CO2_air_in,                &
    & DEBUG_VAR_in,              &
    ! For lakes:
    & drag_wtr_in,               &
    & drag_ice_in,               &
    & t_acoef_wtr_in,            &
    & t_bcoef_wtr_in,            &
    & q_acoef_wtr_in,            &
    & q_bcoef_wtr_in,            &
    & t_acoef_ice_in,            &
    & t_bcoef_ice_in,            &
    & q_acoef_ice_in,            &
    & q_bcoef_ice_in             &
    & )

    USE mo_jsb_physical_constants, ONLY: molarMassDryAir, molarMassCO2

    CLASS(t_jsb_tile_abstract), INTENT(inout)             :: tile
    TYPE(t_jsb_task_options),   INTENT(in)                :: options
    ! REAL(wp), OPTIONAL, DIMENSION(options%nc), INTENT(in) ::                                                       &
    REAL(wp), DIMENSION(options%nc), INTENT(in) ::                                                       &
      & t_air_in, q_air_in, rain_in, snow_in, wind_air_in, wind_10m_in, lw_srf_down_in, swvis_srf_down_in, swnir_srf_down_in, swpar_srf_down_in, &
      & fract_par_diffuse_in, press_srf_in, drag_srf_in, t_acoef_in, t_bcoef_in, q_acoef_in, q_bcoef_in, pch_in, cos_zenith_angle_in, &
      & CO2_air_in, drag_wtr_in, drag_ice_in,                                                                    &
      & t_acoef_wtr_in, t_bcoef_wtr_in, q_acoef_wtr_in, q_bcoef_wtr_in,                                                        &
      & t_acoef_ice_in, t_bcoef_ice_in, q_acoef_ice_in, q_bcoef_ice_in
    REAL(wp), OPTIONAL, DIMENSION(options%nc), INTENT(in) ::                                                       &
      & debug_var_in, dz_srf_in, press_air_in, rho_srf_in

    dsl4jsb_Def_memory(A2L_)

    dsl4jsb_real2D_onChunk :: t_air
    dsl4jsb_real2D_onChunk :: q_air
    dsl4jsb_real2D_onChunk :: rain
    dsl4jsb_real2D_onChunk :: snow
    dsl4jsb_real2D_onChunk :: wind_air
    dsl4jsb_real2D_onChunk :: wind_10m
    dsl4jsb_real2D_onChunk :: lw_srf_down
    dsl4jsb_real2D_onChunk :: swvis_srf_down
    dsl4jsb_real2D_onChunk :: swnir_srf_down
    dsl4jsb_real2D_onChunk :: swpar_srf_down
    dsl4jsb_real2D_onChunk :: fract_par_diffuse
    dsl4jsb_real2D_onChunk :: press_srf
    dsl4jsb_real2D_onChunk :: drag_srf
    dsl4jsb_real2D_onChunk :: t_acoef
    dsl4jsb_real2D_onChunk :: t_bcoef
    dsl4jsb_real2D_onChunk :: q_acoef
    dsl4jsb_real2D_onChunk :: q_bcoef
    dsl4jsb_real2D_onChunk :: pch
    dsl4jsb_real2D_onChunk :: cos_zenith_angle
    dsl4jsb_real2D_onChunk :: CO2_air
    dsl4jsb_real2D_onChunk :: CO2_air_mol
    dsl4jsb_real2D_onChunk :: drag_wtr
    dsl4jsb_real2D_onChunk :: drag_ice
    dsl4jsb_real2D_onChunk :: t_acoef_wtr
    dsl4jsb_real2D_onChunk :: t_acoef_ice
    dsl4jsb_real2D_onChunk :: t_bcoef_wtr
    dsl4jsb_real2D_onChunk :: t_bcoef_ice
    dsl4jsb_real2D_onChunk :: q_acoef_wtr
    dsl4jsb_real2D_onChunk :: q_acoef_ice
    dsl4jsb_real2D_onChunk :: q_bcoef_wtr
    dsl4jsb_real2D_onChunk :: q_bcoef_ice

    INTEGER ::  iblk, ics, ice, nc, i, ic
    TYPE(t_jsb_model), POINTER :: model
    LOGICAL :: use_quincy
    LOGICAL :: tile_contains_lake

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_atm2land'

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    use_quincy = model%config%use_quincy

    ! make this logical accessible in the OpenACC code directly
    tile_contains_lake = tile%contains_lake

    IF (nc /= SIZE(t_air_in,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl4jsb_Get_memory(A2L_)

    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, rain)
    dsl4jsb_Get_var2D_onChunk(A2L_, snow)
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_10m)
    dsl4jsb_Get_var2D_onChunk(A2L_, lw_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_, swvis_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_, swnir_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_, fract_par_diffuse)
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)
    dsl4jsb_Get_var2D_onChunk(A2L_, drag_srf)
    dsl4jsb_Get_var2D_onChunk(A2L_, t_acoef)
    dsl4jsb_Get_var2D_onChunk(A2L_, t_bcoef)
    dsl4jsb_Get_var2D_onChunk(A2L_, q_acoef)
    dsl4jsb_Get_var2D_onChunk(A2L_, q_bcoef)
    dsl4jsb_Get_var2D_onChunk(A2L_, pch)
    dsl4jsb_Get_var2D_onChunk(A2L_, cos_zenith_angle)
    dsl4jsb_Get_var2D_onChunk(A2L_, CO2_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, CO2_air_mol)
    dsl4jsb_Get_var2D_onChunk(A2L_, drag_wtr)
    dsl4jsb_Get_var2D_onChunk(A2L_, drag_ice)
    dsl4jsb_Get_var2D_onChunk(A2L_, t_acoef_wtr)
    dsl4jsb_Get_var2D_onChunk(A2L_, t_acoef_ice)
    dsl4jsb_Get_var2D_onChunk(A2L_, t_bcoef_wtr)
    dsl4jsb_Get_var2D_onChunk(A2L_, t_bcoef_ice)
    dsl4jsb_Get_var2D_onChunk(A2L_, q_acoef_wtr)
    dsl4jsb_Get_var2D_onChunk(A2L_, q_acoef_ice)
    dsl4jsb_Get_var2D_onChunk(A2L_, q_bcoef_wtr)
    dsl4jsb_Get_var2D_onChunk(A2L_, q_bcoef_ice)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

    ! IF (PRESENT(DEBUG_VAR)) THEN
    !   !$ACC LOOP GANG(STATIC: 1) VECTOR
    !   DO i=1,nc
    !     DEBUG_VAR_ptr(i) = DEBUG_VAR(i)
    !   END DO
    ! END IF

      !$ACC LOOP GANG VECTOR PRIVATE(ic)
      DO i=1,nc
        t_air(i) = t_air_in(i)
        q_air(i) = q_air_in(i)
    ! IF (PRESENT(press_air)) THEN
    !   !$ACC LOOP GANG(STATIC: 1) VECTOR
    !   DO i=1,nc
    !     press_air(i) = press_air_in(i)
    !   END DO
    ! END IF
        rain(i) = rain_in(i)
        snow(i) = snow_in(i)
        wind_air(i) = wind_air_in(i)
        wind_10m(i) = wind_10m_in(i)
        lw_srf_down(i) = lw_srf_down_in(i)
        swvis_srf_down(i) = swvis_srf_down_in(i)
        swnir_srf_down(i) = swnir_srf_down_in(i)
        swpar_srf_down(i) = swpar_srf_down_in(i)
        fract_par_diffuse(i) = fract_par_diffuse_in(i)
    ! IF (PRESENT(dz_srf)) THEN
    !   !$ACC LOOP GANG(STATIC: 1) VECTOR
    !   DO i=1,nc
    !     dz_srf(i) = dz_srf_in(i)
    !   END DO
    ! END IF
        press_srf(i) = press_srf_in(i)
        drag_srf(i) = drag_srf_in(i)
    ! IF (PRESENT(rho_srf)) THEN
    !   !$ACC LOOP GANG(STATIC: 1) VECTOR
    !   DO i=1,nc
    !     rho_srf(i) = rho_srf_in(i)
    !   END DO
    ! END IF
        t_acoef(i) = t_acoef_in(i)
        t_bcoef(i) = t_bcoef_in(i)
        q_acoef(i) = q_acoef_in(i)
        q_bcoef(i) = q_bcoef_in(i)
        pch(i) = pch_in(i)
        cos_zenith_angle(i) = cos_zenith_angle_in(i)
        CO2_air(i) = CO2_air_in(i)
        ! Convert CO2 mass mixing ratio [kg/kg] to particle mixing ratio [mol/mol]
        CO2_air_mol(i) = CO2_air_in(i) * molarMassDryAir / molarMassCO2
        ! convert CO2 from "molar ratio (volume)" to "co2 mixing ratio ppmv" (quincy | used in e.g. update_canopy_fluxes)
        ! IF (use_quincy) THEN
        !   CO2_mixing_ratio(i) = CO2_air_mol, ic) * 1000000._wp
        ! END IF

          drag_wtr(i) = drag_wtr_in(i)
          drag_ice(i) = drag_ice_in(i)
          t_acoef_wtr(i) = t_acoef_wtr_in(i)
          t_bcoef_wtr(i) = t_bcoef_wtr_in(i)
          q_acoef_wtr(i) = q_acoef_wtr_in(i)
          q_bcoef_wtr(i) = q_bcoef_wtr_in(i)
          t_acoef_ice(i) = t_acoef_ice_in(i)
          t_bcoef_ice(i) = t_bcoef_ice_in(i)
          q_acoef_ice(i) = q_acoef_ice_in(i)
          q_bcoef_ice(i) = q_bcoef_ice_in(i)

  END DO

    !$ACC END PARALLEL

    !$ACC WAIT(1)

  END SUBROUTINE update_atm2land

  SUBROUTINE update_land2atm(tile, options,                                          &
    & t_srf, t_srf_rad, t_eff_srf, qsat_srf, s_srf,                                  &
    & fact_q_air, fact_qsat_srf, evapopot,                                           &
    & evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap,                  &
    & rough_h_srf, rough_m_srf, q_snocpymlt,                                         &
    & alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif,                            &
    & kh, km, kh_neutral, km_neutral, CO2_flux,                                      &
    & t_lwtr, t_lice, qsat_lwtr, qsat_lice, s_lwtr, s_lice,                          &
    & evapo_wtr, latent_hflx_wtr, sensible_hflx_wtr,                                 &
    & evapo_ice, latent_hflx_ice, sensible_hflx_ice,                                 &
    & ice_fract_lake,                                                                &
    & alb_vis_dir_wtr, alb_vis_dif_wtr, alb_nir_dir_wtr, alb_nir_dif_wtr,            &
    & albedo_lwtr, albedo_lice)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! REAL(wp), OPTIONAL, DIMENSION(options%nc), INTENT(out) ::         &
    REAL(wp), DIMENSION(options%nc), INTENT(out) ::         &
      & t_srf, t_eff_srf, qsat_srf, s_srf,                 &
      & fact_q_air, fact_qsat_srf,                          &
      & evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap, &
      & rough_h_srf, rough_m_srf, q_snocpymlt,                        &
      & t_lwtr, t_lice, qsat_lwtr, qsat_lice, s_lwtr, s_lice,         &
      & evapo_wtr, latent_hflx_wtr, sensible_hflx_wtr,                &
      & evapo_ice, latent_hflx_ice, sensible_hflx_ice,                &
      & ice_fract_lake,                                               &
      & albedo_lwtr, albedo_lice,                    &
      & alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif,           &
      & CO2_flux
    REAL(wp), OPTIONAL, DIMENSION(options%nc), INTENT(out) ::         &
      & alb_vis_dir_wtr, alb_vis_dif_wtr, alb_nir_dir_wtr,            &
      & alb_nir_dif_wtr,                     &
      & t_srf_rad, evapopot,                          &
      & kh, km, kh_neutral, km_neutral

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(CARBON_)

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), POINTER :: land_tile
    dsl4jsb_Def_memory_tile(TURB_, land_tile)

    INTEGER  :: iblk, ics, ice, nc, i, j, it, nt
    LOGICAL  :: use_tmx, tile_is_carbon_active
    REAL(wp) :: fact_lake(options%nc), fact_land(options%nc)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_land2atm'

    dsl4jsb_Real2D_onChunk :: &
      & t_ptr, &                  ! tmx
      & t_filt_ptr, &
      & t_rad_ptr, &
      & t_eff_ptr, &
      & qsat_star_ptr, &
      & s_star_ptr, &
      & fact_q_air_ptr, &
      & fact_q_air_land_tile, &
      & fact_qsat_srf_ptr, &
      & fact_qsat_srf_land_tile, &
      & evapopot_ptr, &
      & evapotrans_ptr, &         ! tmx
      & latent_hflx_ptr, &        ! tmx
      & sensible_hflx_ptr, &      ! tmx
      & evapotrans_lnd_ptr, &
      & latent_hflx_lnd_ptr, &
      & sensible_hflx_lnd_ptr, &
      & forc_hflx_ptr, &
      & heat_cap_ptr, &
      & rough_h_ptr, &
      & rough_m_ptr, &
      & q_snocpymlt_ptr, &
      & alb_vis_lnd_ptr, &
      & alb_nir_lnd_ptr, &
      & kh_ptr, &
      & km_ptr, &
      & kh_neutral_ptr, &
      & km_neutral_ptr, &
      & co2flux_npp_2_atm_ta_ptr, &
      & co2flux_soilresp_2_atm_ta_ptr, &
      & co2flux_herb_2_atm_ta_ptr, &
      & co2flux_fire_all_2_atm_ta_ptr, &
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
      & sensible_hflx_ice_ptr

    ! avoid compiler warnings about dummy arguments not being used
    ! IF (PRESENT(alb_nir_dif_wtr)) CONTINUE
    ! IF (PRESENT(alb_nir_dir_wtr)) CONTINUE
    ! IF (PRESENT(alb_vis_dif_wtr)) CONTINUE
    ! IF (PRESENT(alb_vis_dir_wtr)) CONTINUE

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (nc /= SIZE(t_srf,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    model => Get_model(tile%owner_model_id)

    use_tmx = model%config%use_tmx
    tile_is_carbon_active = tile%Is_process_active(CARBON_)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(RAD_)
    IF (tile_is_carbon_active) THEN
      dsl4jsb_Get_memory(CARBON_)
    END IF

    !$ACC DATA CREATE(fact_land, fact_lake)

    IF (tile%contains_lake .AND. .NOT. use_tmx .AND. .NOT. jsbach_runs_standalone()) THEN

      ! With lakes, re-scale _lnd, _lice and _lwtr fluxes given back to the atmosphere so that they are relative to the
      ! jsbach grid box, not/only including lakes (JSBACH grid boxes do not include ocean fractions.). Re-scaling is
      ! necessary, as in the atmosphere fluxes are regarded as land-only (or lake-only) fluxes, and lake, land and ocean
      ! fractions are considered separately.
      ! Note that, if lakes are considered, the lake tile must be a direct child of the root tile.

      DO j=1,SIZE(tile%lcts)
        IF (tile%lcts(j)%id == LAKE_TYPE) EXIT         ! We get index j of the lake tile
      END DO

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) PRESENT(tile%lcts)
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        fact_lake(i) = tile%lcts(j)%fract(ics+i-1,iblk)      ! Lake fraction
        fact_land(i) = 1._wp - fact_lake(i)            ! Land fraction
        IF (fact_lake(i) > 0._wp) THEN
          fact_lake(i) = 1._wp / fact_lake(i)          ! Lake factor: 1/lake_fraction
        ELSE
          fact_lake(i) = 1._wp
        END IF
        IF (fact_land(i) > 0._wp) THEN
          fact_land(i) = 1._wp / fact_land(i)          ! Land factor: 1/land_fraction
        ELSE
          fact_land(i) = 1._wp
        END IF
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL

    ELSE    ! without lakes

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO i=1,nc
        fact_land(i) = 1._wp
        fact_lake(i) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP

    END IF

    ! Exchange fields on the land tile
    ! --------------------------------

    ! For some variables the atmosphere expects values for the land tile (excluding lakes).
    ! Thus we define a pointer to the land tile.
    nt = tile%Get_no_of_children()
    DO it=1,nt
      IF (it == 1) THEN
        land_tile => tile%Get_first_child_tile()
      ELSE
        land_tile => land_tile%Get_next_sibling_tile()
      END IF
      ! Exit if the land tile was found
      IF (TRIM(land_tile%name) == 'land') EXIT
    END DO
    dsl4jsb_Get_memory_tile(TURB_, land_tile)

    t_eff_ptr        => dsl4jsb_var2D_onChunk(SEB_,   t_eff4)
    qsat_star_ptr    => dsl4jsb_var2D_onChunk(SEB_,   qsat_star)
    s_star_ptr       => dsl4jsb_var2D_onChunk(SEB_,   s_star)
    dsl4jsb_Get_var2d_onChunk_tile_name(TURB_, fact_q_air, land_tile)
    dsl4jsb_Get_var2d_onChunk_tile_name(TURB_, fact_qsat_srf, land_tile)
    evapopot_ptr     => dsl4jsb_var2D_onChunk(HYDRO_, evapopot) ! For offline only
    forc_hflx_ptr    => dsl4jsb_var2D_onChunk(SEB_,   forc_hflx)    ! TODO Not used
    heat_cap_ptr     => dsl4jsb_var2D_onChunk(SEB_,   heat_cap)
    rough_h_ptr      => dsl4jsb_var2D_onChunk(TURB_,  rough_h)
    rough_m_ptr      => dsl4jsb_var2D_onChunk(TURB_,  rough_m)
    q_snocpymlt_ptr  => dsl4jsb_var2D_onChunk(HYDRO_, q_snocpymlt)
    IF (use_tmx) THEN
      t_ptr                 => dsl4jsb_var2D_onChunk(SEB_,   t)
      t_rad_ptr             => dsl4jsb_var2D_onChunk(SEB_,   t_rad4)
      evapotrans_ptr        => dsl4jsb_var2D_onChunk(HYDRO_, evapotrans)
      latent_hflx_ptr       => dsl4jsb_var2D_onChunk(SEB_,   latent_hflx)
      sensible_hflx_ptr     => dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx)
      alb_vis_lnd_ptr       => dsl4jsb_var2D_onChunk(RAD_,   alb_vis)      ! TODO
      alb_vis_lnd_ptr       => dsl4jsb_var2D_onChunk(RAD_,   alb_vis)      ! TODO
      alb_nir_lnd_ptr       => dsl4jsb_var2D_onChunk(RAD_,   alb_nir)      ! TODO
      alb_nir_lnd_ptr       => dsl4jsb_var2D_onChunk(RAD_,   alb_nir)
      kh_ptr                => dsl4jsb_var2D_onChunk(TURB_,  kh)
      km_ptr                => dsl4jsb_var2D_onChunk(TURB_,  km)
      kh_neutral_ptr        => dsl4jsb_var2D_onChunk(TURB_,  kh_neutral)
      km_neutral_ptr        => dsl4jsb_var2D_onChunk(TURB_,  km_neutral)
    ELSE
      t_filt_ptr            => dsl4jsb_var2D_onChunk(SEB_,   t_filt)
      evapotrans_lnd_ptr    => dsl4jsb_var2D_onChunk(HYDRO_, evapotrans_lnd)
      latent_hflx_lnd_ptr   => dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_lnd)
      sensible_hflx_lnd_ptr => dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_lnd)
      alb_vis_lnd_ptr       => dsl4jsb_var2D_onChunk(RAD_,   alb_vis_lnd)  ! TODO
      alb_vis_lnd_ptr       => dsl4jsb_var2D_onChunk(RAD_,   alb_vis_lnd)  ! TODO
      alb_nir_lnd_ptr       => dsl4jsb_var2D_onChunk(RAD_,   alb_nir_lnd)  ! TODO
      alb_nir_lnd_ptr       => dsl4jsb_var2D_onChunk(RAD_,   alb_nir_lnd)
    END IF

    IF (tile_is_carbon_active) THEN
      co2flux_npp_2_atm_ta_ptr      => dsl4jsb_var2D_onChunk(CARBON_,co2flux_npp_2_atm_ta)
      co2flux_soilresp_2_atm_ta_ptr => dsl4jsb_var2D_onChunk(CARBON_,co2flux_soilresp_2_atm_ta)
      co2flux_herb_2_atm_ta_ptr     => dsl4jsb_var2D_onChunk(CARBON_,co2flux_herb_2_atm_ta)
      co2flux_fire_all_2_atm_ta_ptr => dsl4jsb_var2D_onChunk(CARBON_,co2flux_fire_all_2_atm_ta)
    END IF

    IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone()) THEN
      t_lwtr_ptr            => dsl4jsb_var2D_onChunk(SEB_,   t_lwtr)
      qsat_lwtr_ptr         => dsl4jsb_var2D_onChunk(SEB_,   qsat_lwtr)
      s_lwtr_ptr            => dsl4jsb_var2D_onChunk(SEB_,   s_lwtr)
      evapo_wtr_ptr         => dsl4jsb_var2D_onChunk(HYDRO_, evapo_wtr)
      latent_hflx_wtr_ptr   => dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_wtr)
      sensible_hflx_wtr_ptr => dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_wtr)
      albedo_lwtr_ptr       => dsl4jsb_var2D_onChunk(RAD_,   albedo_lwtr)
      fract_lice_ptr        => dsl4jsb_var2D_onChunk(SEB_,   fract_lice)
    END IF

    IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone() .AND. dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      t_lice_ptr            => dsl4jsb_var2D_onChunk(SEB_,   t_lice)
      qsat_lice_ptr         => dsl4jsb_var2D_onChunk(SEB_,   qsat_lice)
      s_lice_ptr            => dsl4jsb_var2D_onChunk(SEB_,   s_lice)
      evapo_ice_ptr         => dsl4jsb_var2D_onChunk(HYDRO_, evapo_ice)
      latent_hflx_ice_ptr   => dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_ice)
      sensible_hflx_ice_ptr => dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_ice)
      albedo_lice_ptr       => dsl4jsb_var2D_onChunk(RAD_,   albedo_lice)
    END IF

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

    ! IF (use_tmx) THEN
    !   IF (PRESENT(t_srf)) THEN
    !     !$ACC LOOP GANG(STATIC: 1) VECTOR
    !     DO i=1,nc
    !       t_srf(i) = t_ptr(i)
    !     END DO
    !   END IF
    !   IF (PRESENT(t_srf_rad)) THEN
    !     !$ACC LOOP GANG(STATIC: 1) VECTOR
    !     DO i=1,nc
    !       t_srf_rad(i) = t_rad_ptr(i)**0.25_wp
    !     END DO
    !   END IF
    !   IF (PRESENT(evapotrans)) THEN
    !     !$ACC LOOP GANG(STATIC: 1) VECTOR
    !     DO i=1,nc
    !       evapotrans(i) = evapotrans_ptr(i)
    !     END DO
    !   END IF
    !   IF (PRESENT(latent_hflx)) THEN
    !     !$ACC LOOP GANG(STATIC: 1) VECTOR
    !     DO i=1,nc
    !       latent_hflx(i) = latent_hflx_ptr(i)
    !     END DO
    !   END IF
    !   IF (PRESENT(sensible_hflx)) THEN
    !     !$ACC LOOP GANG(STATIC: 1) VECTOR
    !     DO i=1,nc
    !       sensible_hflx(i) = sensible_hflx_ptr(i)
    !     END DO
    !   END IF
    ! ELSE
      ! IF (PRESENT(t_srf)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          t_srf(i) = t_filt_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(evapotrans)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          evapotrans(i) = fact_land(i) * evapotrans_lnd_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(latent_hflx)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          latent_hflx(i) = fact_land(i) * latent_hflx_lnd_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(sensible_hflx)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          sensible_hflx(i) = fact_land(i) * sensible_hflx_lnd_ptr(i)
        END DO
      ! END IF
    ! END IF
    ! IF (PRESENT(t_eff_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        t_eff_srf(i) = t_eff_ptr(i)**0.25_wp
      END DO
    ! END IF
    ! IF (PRESENT(qsat_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        qsat_srf(i) = qsat_star_ptr(i)
      END DO
    ! END IF
    ! IF (PRESENT(s_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        s_srf(i) = s_star_ptr(i)
      END DO
    ! END IF
    ! IF (PRESENT(fact_q_air)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        fact_q_air(i) = fact_q_air_land_tile(i)
      END DO
    ! END IF
    ! IF (PRESENT(fact_qsat_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        fact_qsat_srf(i) = fact_qsat_srf_land_tile(i)
      END DO
    ! END IF
    IF (PRESENT(evapopot)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        evapopot(i) = evapopot_ptr(i)
      END DO
    END IF
    ! IF (PRESENT(grnd_hflx)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        grnd_hflx(i) = forc_hflx_ptr(i)
      END DO
    ! END IF
    ! IF (PRESENT(grnd_hcap)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        grnd_hcap(i) = heat_cap_ptr(i)
      END DO
    ! END IF
    ! IF (PRESENT(rough_h_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        rough_h_srf(i) = rough_h_ptr(i)
      END DO
    ! END IF
    ! IF (PRESENT(rough_m_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        rough_m_srf(i) = rough_m_ptr(i)
      END DO
    ! END IF
    ! IF (PRESENT(q_snocpymlt)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        q_snocpymlt(i) = fact_land(i) * q_snocpymlt_ptr(i)
      END DO
    ! END IF
    ! IF (PRESENT(alb_vis_dir)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        alb_vis_dir(i) = alb_vis_lnd_ptr(i)
      END DO
    ! END IF
    ! IF (PRESENT(alb_vis_dif)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        alb_vis_dif(i) = alb_vis_lnd_ptr(i)
      END DO
    ! END IF
    ! IF (PRESENT(alb_nir_dir)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        alb_nir_dir(i) = alb_nir_lnd_ptr(i)
      END DO
    ! END IF
    ! IF (PRESENT(alb_nir_dif)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        alb_nir_dif(i) = alb_nir_lnd_ptr(i)
      END DO
    ! END IF
    ! IF (PRESENT(CO2_flux)) THEN
      IF (tile_is_carbon_active) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          CO2_flux(i) = fact_land(i) *              &
            & (   co2flux_npp_2_atm_ta_ptr(i)       &
            &   + co2flux_soilresp_2_atm_ta_ptr(i)  &
            &   + co2flux_herb_2_atm_ta_ptr(i)      &
            &   + co2flux_fire_all_2_atm_ta_ptr(i)  &
            & )
        END DO
      ELSE
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          CO2_flux(i) = 0._wp
        END DO
      END IF
    ! END IF
    ! IF (PRESENT(kh)) THEN
    !   !$ACC LOOP GANG(STATIC: 1) VECTOR
    !   DO i=1,nc
    !     kh(i) = kh_ptr(i)
    !   END DO
    ! END IF
    ! IF (PRESENT(km)) THEN
    !   !$ACC LOOP GANG(STATIC: 1) VECTOR
    !   DO i=1,nc
    !     km(i) = km_ptr(i)
    !   END DO
    ! END IF
    ! IF (PRESENT(kh_neutral)) THEN
    !   !$ACC LOOP GANG(STATIC: 1) VECTOR
    !   DO i=1,nc
    !     kh_neutral(i) = kh_neutral_ptr(i)
    !   END DO
    ! END IF
    ! IF (PRESENT(km_neutral)) THEN
    !   !$ACC LOOP GANG(STATIC: 1) VECTOR
    !   DO i=1,nc
    !     km_neutral(i) = km_neutral_ptr(i)
    !   END DO
    ! END IF

    !$ACC END PARALLEL

    ! Exchange fields on lake water fractions
    ! ---------------------------------------

    IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone()) THEN

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

      ! IF (PRESENT(t_lwtr)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          t_lwtr(i) = t_lwtr_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(qsat_lwtr)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          qsat_lwtr(i) = qsat_lwtr_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(s_lwtr)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          s_lwtr(i) = s_lwtr_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(evapo_wtr)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          evapo_wtr(i) = fact_lake(i) * evapo_wtr_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(latent_hflx_wtr)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          latent_hflx_wtr(i) = fact_lake(i) * latent_hflx_wtr_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(sensible_hflx_wtr)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          sensible_hflx_wtr(i) = fact_lake(i) * sensible_hflx_wtr_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(albedo_lwtr)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          albedo_lwtr(i) = albedo_lwtr_ptr(i)
        END DO
      ! END IF

    ! Exchange fields on the ice fraction of lakes
    ! --------------------------------------------

      ! IF (PRESENT(ice_fract_lake)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          ice_fract_lake(i) = fract_lice_ptr(i)
        END DO
      ! END IF

      !$ACC END PARALLEL

    ELSE

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  
      ! IF (PRESENT(t_lwtr)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          t_lwtr(i) = 0._wp
        END DO
      ! END IF
      ! IF (PRESENT(qsat_lwtr)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          qsat_lwtr(i) = 0._wp
        END DO
      ! END IF
      ! IF (PRESENT(s_lwtr)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          s_lwtr(i) = 0._wp
        END DO
      ! END IF
      ! IF (PRESENT(evapo_wtr)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          evapo_wtr(i) = 0._wp
        END DO
      ! END IF
      ! IF (PRESENT(latent_hflx_wtr)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          latent_hflx_wtr(i) = 0._wp
        END DO
      ! END IF
      ! IF (PRESENT(sensible_hflx_wtr)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          sensible_hflx_wtr(i) = 0._wp
        END DO
      ! END IF
      ! IF (PRESENT(albedo_lwtr)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          albedo_lwtr(i) = 0.07_wp
        END DO
      ! END IF

    ! Exchange fields on the ice fraction of lakes
    ! --------------------------------------------

      ! IF (PRESENT(ice_fract_lake)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          ice_fract_lake(i) = 0._wp
        END DO
      ! END IF

      !$ACC END PARALLEL

    END IF

    IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone() .AND. dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

      ! IF (PRESENT(t_lice)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          t_lice(i) = t_lice_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(qsat_lice)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          qsat_lice(i) = qsat_lice_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(s_lice)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          s_lice(i) = s_lice_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(evapo_ice)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          evapo_ice(i) = fact_lake(i) * evapo_ice_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(latent_hflx_ice)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          latent_hflx_ice(i) = fact_lake(i) * latent_hflx_ice_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(sensible_hflx_ice)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          sensible_hflx_ice(i) = fact_lake(i) * sensible_hflx_ice_ptr(i)
        END DO
      ! END IF
      ! IF (PRESENT(albedo_lice)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          albedo_lice(i) = albedo_lice_ptr(i)
        END DO
      ! END IF

      !$ACC END PARALLEL

    ELSE

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

      ! IF (PRESENT(t_lice)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          t_lice(i) = 273.15_wp
        END DO
      ! END IF
      ! IF (PRESENT(qsat_lice)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          qsat_lice(i) = 0.0075_wp
        END DO
      ! END IF
      ! IF (PRESENT(s_lice)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          s_lice(i) = 2.9E5_wp
        END DO
      ! END IF
      ! IF (PRESENT(evapo_ice)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          evapo_ice(i) = 0._wp
        END DO
      ! END IF
      ! IF (PRESENT(latent_hflx_ice)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          latent_hflx_ice(i) = 0._wp
        END DO
      ! END IF
      ! IF (PRESENT(sensible_hflx_ice)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          sensible_hflx_ice(i) = 0._wp
        END DO
      ! END IF
      ! IF (PRESENT(albedo_lice)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          albedo_lice(i) = 0.55_wp
        END DO
      ! END IF

      !$ACC END PARALLEL

    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE update_land2atm

#endif
END MODULE mo_atmland_interface
