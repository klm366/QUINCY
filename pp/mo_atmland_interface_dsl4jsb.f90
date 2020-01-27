!
! This file was automatically generated from JSBACH source file:
!    src/atmland/mo_atmland_interface.f90
! with md5 hash: 518091d1d59afcd2933a812b1adf1167
!
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

  USE mo_jsb_process_class, ONLY: a2l_, seb_, turb_, hydro_, rad_, carbon_
  USE mo_seb_config_class, ONLY: t_seb_config
  USE mo_a2l_memory_class, ONLY: t_a2l_memory
  USE mo_seb_memory_class, ONLY: t_seb_memory
  USE mo_turb_memory_class, ONLY: t_turb_memory
  USE mo_hydro_memory_class, ONLY: t_hydro_memory
  USE mo_rad_memory_class, ONLY: t_rad_memory
  USE mo_carbon_memory_class, ONLY: t_carbon_memory


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

    TYPE(t_a2l_memory), POINTER :: a2l__mem

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

    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT

    !$acc parallel default(present)

    IF (PRESENT(t_air)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%t_air%ptr(ics+i-1,iblk) = t_air(i)
      END DO
    END IF
    IF (PRESENT(q_air)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%q_air%ptr(ics+i-1,iblk) = q_air(i)
      END DO
    END IF
    IF (PRESENT(rain)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%rain%ptr(ics+i-1,iblk)= rain(i)
      END DO
    END IF
    IF (PRESENT(snow)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%snow%ptr(ics+i-1,iblk) = snow(i)
      END DO
    END IF
    IF (PRESENT(wind_air)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%wind_air%ptr(ics+i-1,iblk) = wind_air(i)
      END DO
    END IF
    IF (PRESENT(wind_10m)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%wind_10m%ptr(ics+i-1,iblk) = wind_10m(i)
      END DO
    END IF
    IF (PRESENT(lw_srf_down)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%lw_srf_down%ptr(ics+i-1,iblk) = lw_srf_down(i)
      END DO
    END IF
    IF (PRESENT(swvis_srf_down)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%swvis_srf_down%ptr(ics+i-1,iblk) = swvis_srf_down(i)
      END DO
    END IF
    IF (PRESENT(swnir_srf_down)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%swnir_srf_down%ptr(ics+i-1,iblk) = swnir_srf_down(i)
      END DO
    END IF
    IF (PRESENT(swpar_srf_down)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%swpar_srf_down%ptr(ics+i-1,iblk) = swpar_srf_down(i)
      END DO
    END IF
    IF (PRESENT(fract_par_diffuse)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%fract_par_diffuse%ptr(ics+i-1,iblk)  = fract_par_diffuse(i)
      END DO
    END IF
    IF (PRESENT(press_srf)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%press_srf%ptr(ics+i-1,iblk) = press_srf(i)
      END DO
    END IF
    IF (PRESENT(drag_srf)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%drag_srf%ptr(ics+i-1,iblk) = drag_srf(i)
      END DO
    END IF
    IF (PRESENT(t_acoef)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%t_acoef%ptr(ics+i-1,iblk) = t_acoef(i)
      END DO
    END IF
    IF (PRESENT(t_bcoef)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%t_bcoef%ptr(ics+i-1,iblk) = t_bcoef(i)
      END DO
    END IF
    IF (PRESENT(q_acoef)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%q_acoef%ptr(ics+i-1,iblk) = q_acoef(i)
      END DO
    END IF
    IF (PRESENT(q_bcoef)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%q_bcoef%ptr(ics+i-1,iblk) = q_bcoef(i)
      END DO
    END IF
    IF (PRESENT(pch)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%pch%ptr(ics+i-1,iblk) = pch(i)
      END DO
    END IF
    IF (PRESENT(cos_zenith_angle)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%cos_zenith_angle%ptr(ics+i-1,iblk) = cos_zenith_angle(i)
      END DO
    END IF
    IF (PRESENT(CO2_air)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a2l__mem%co2_air%ptr(ics+i-1,iblk) = co2_air(i)
        ! Convert CO2 mass mixing ratio [kg/kg] to particle mixing ratio [mol/mol]
        a2l__mem%co2_air_mol%ptr(ics+i-1,iblk) = co2_air(i) * molarMassDryAir / molarMassCO2
      END DO
    END IF

    IF (tile_contains_lake) THEN

      IF (PRESENT(drag_wtr)) THEN
      !$acc loop gang vector
        DO i=1,nc
          a2l__mem%drag_wtr%ptr(ics+i-1,iblk) = drag_wtr(i)
        END DO
      END IF
      IF (PRESENT(drag_ice)) THEN
      !$acc loop gang vector
        DO i=1,nc
          a2l__mem%drag_ice%ptr(ics+i-1,iblk) = drag_ice(i)
        END DO
      END IF
      IF (PRESENT(t_acoef_wtr)) THEN
      !$acc loop gang vector
        DO i=1,nc
          a2l__mem%t_acoef_wtr%ptr(ics+i-1,iblk) = t_acoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(t_bcoef_wtr)) THEN
      !$acc loop gang vector
        DO i=1,nc
          a2l__mem%t_bcoef_wtr%ptr(ics+i-1,iblk) = t_bcoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(q_acoef_wtr)) THEN
      !$acc loop gang vector
        DO i=1,nc
          a2l__mem%q_acoef_wtr%ptr(ics+i-1,iblk) = q_acoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(q_bcoef_wtr)) THEN
      !$acc loop gang vector
        DO i=1,nc
          a2l__mem%q_bcoef_wtr%ptr(ics+i-1,iblk) = q_bcoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(t_acoef_ice)) THEN
      !$acc loop gang vector
        DO i=1,nc
          a2l__mem%t_acoef_ice%ptr(ics+i-1,iblk) = t_acoef_ice(i)
        END DO
      END IF
      IF (PRESENT(t_bcoef_ice)) THEN
      !$acc loop gang vector
        DO i=1,nc
          a2l__mem%t_bcoef_ice%ptr(ics+i-1,iblk) = t_bcoef_ice(i)
        END DO
      END IF
      IF (PRESENT(q_acoef_ice)) THEN
      !$acc loop gang vector
        DO i=1,nc
          a2l__mem%q_acoef_ice%ptr(ics+i-1,iblk) = q_acoef_ice(i)
        END DO
      END IF
      IF (PRESENT(q_bcoef_ice)) THEN
      !$acc loop gang vector
        DO i=1,nc
          a2l__mem%q_bcoef_ice%ptr(ics+i-1,iblk) = q_bcoef_ice(i)
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

    TYPE(t_seb_config), POINTER :: seb__conf
    TYPE(t_seb_memory), POINTER :: seb__mem
    TYPE(t_turb_memory), POINTER :: turb__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem
    TYPE(t_rad_memory), POINTER :: rad__mem
    TYPE(t_carbon_memory), POINTER :: carbon__mem

    TYPE(t_jsb_model), POINTER :: model

    INTEGER ::  iblk, ics, ice, nc, i, j
    REAL(wp) :: fact_lake(options%nc), fact_land(options%nc)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_land2atm'

    REAL(wp), POINTER, DIMENSION(:) :: &
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

    SELECT TYPE (cc => model%processes(seb_)%p%config); TYPE IS (t_seb_config) ; seb__conf => cc ; END SELECT
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(turb_)%p); TYPE IS (t_turb_memory); turb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(rad_)%p); TYPE IS (t_rad_memory); rad__mem => mm ; END SELECT
    IF (tile%Is_process_active(CARBON_)) THEN
      SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT
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
        t_filt_ptr => seb__mem%t_filt%ptr(ics:ice,iblk)
        !$acc data present(t_srf, t_filt_ptr)
        !$acc parallel loop
        DO i=1,nc
          t_srf(i) = t_filt_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(t_eff_srf)) THEN
        t_eff_ptr => seb__mem%t_eff%ptr(ics:ice,iblk)
        !$acc data present(t_eff_srf, t_eff_ptr)
        !$acc parallel loop
        DO i=1,nc
          t_eff_srf(i) = t_eff_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(qsat_srf)) THEN
        qsat_star_ptr => seb__mem%qsat_star%ptr(ics:ice,iblk)
        !$acc data present(qsat_srf, qsat_star_ptr)
        !$acc parallel loop
        DO i=1,nc
          qsat_srf(i) = qsat_star_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(s_srf)) THEN
        s_star_ptr => seb__mem%s_star%ptr(ics:ice,iblk)
        !$acc data present(s_srf, s_star_ptr)
        !$acc parallel loop
        DO i=1,nc
          s_srf(i) = s_star_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(fact_q_air)) THEN
        fact_q_air_ptr => turb__mem%fact_q_air%ptr(ics:ice,iblk)
        !$acc data present(fact_q_air, fact_q_air_ptr)
        !$acc parallel loop
        DO i=1,nc
          fact_q_air(i) = fact_q_air_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(fact_qsat_srf)) THEN
        fact_qsat_srf_ptr => turb__mem%fact_qsat_srf%ptr(ics:ice,iblk)
        !$acc data present(fact_qsat_srf, fact_qsat_srf_ptr)
        !$acc parallel loop
        DO i=1,nc
          fact_qsat_srf(i) = fact_qsat_srf_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(evapopot)) THEN
        evapopot_ptr => hydro__mem%evapopot%ptr(ics:ice,iblk) ! For offline only
        !$acc data present(evapopot, evapopot_ptr)
        !$acc parallel loop
        DO i=1,nc
          evapopot(i) = evapopot_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(evapotrans)) THEN
        evapotrans_lnd_ptr => hydro__mem%evapotrans_lnd%ptr(ics:ice,iblk)
        !$acc data present(fact_land, evapotrans)
        !$acc parallel loop
        DO i=1,nc
          evapotrans(i) = fact_land(i) * evapotrans_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(latent_hflx)) THEN
        latent_hflx_lnd_ptr => seb__mem%latent_hflx_lnd%ptr(ics:ice,iblk)
        !$acc data present(latent_hflx, latent_hflx_lnd_ptr, fact_land)
        !$acc parallel loop
        DO i=1,nc
          latent_hflx(i) = fact_land(i) * latent_hflx_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(sensible_hflx)) THEN
        sensible_hflx_lnd_ptr => seb__mem%sensible_hflx_lnd%ptr(ics:ice,iblk)
        !$acc data present(sensible_hflx, fact_land, sensible_hflx_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          sensible_hflx(i) = fact_land(i) * sensible_hflx_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(grnd_hflx)) THEN
        forc_hflx_ptr => seb__mem%forc_hflx%ptr(ics:ice,iblk)    ! TODO Not used
        !$acc data present(grnd_hflx, forc_hflx_ptr)
        !$acc parallel loop
        DO i=1,nc
          grnd_hflx(i) = forc_hflx_ptr(i)            
        END DO
        !$acc end data
      END IF
      IF (PRESENT(grnd_hcap)) THEN
        heat_cap_old_ptr => seb__mem%heat_cap_old%ptr(ics:ice,iblk)
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
        rough_h_ptr => turb__mem%rough_h%ptr(ics:ice,iblk)
        !$acc data present(rough_h_srf, rough_h_ptr)
        !$acc parallel loop
        DO i=1,nc
        rough_h_srf(i) = rough_h_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(rough_m_srf)) THEN
        rough_m_ptr => turb__mem%rough_m%ptr(ics:ice,iblk)
        !$acc data present(rough_m_ptr, rough_m_srf)
        !$acc parallel loop
        DO i=1,nc
          rough_m_srf(i) = rough_m_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(q_snocpymlt)) THEN
        q_snocpymlt_ptr => hydro__mem%q_snocpymlt%ptr(ics:ice,iblk)
        !$acc data present(q_snocpymlt, q_snocpymlt_ptr, fact_land)
        !$acc parallel loop
        DO i=1,nc
          q_snocpymlt(i) = fact_land(i) * q_snocpymlt_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_vis_dir)) THEN
        alb_vis_lnd_ptr => rad__mem%alb_vis_lnd%ptr(ics:ice,iblk)  ! TODO
        !$acc data present(alb_vis_lnd_ptr, alb_vis_dir)
        !$acc parallel loop
        DO i=1,nc
          alb_vis_dir(i) = alb_vis_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_vis_dif)) THEN
        alb_vis_lnd_ptr => rad__mem%alb_vis_lnd%ptr(ics:ice,iblk)  ! TODO
        !$acc data present(alb_vis_dif, alb_vis_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          alb_vis_dif(i) = alb_vis_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_nir_dir)) THEN
        alb_nir_lnd_ptr => rad__mem%alb_nir_lnd%ptr(ics:ice,iblk)  ! TODO
        !$acc data present(alb_nir_dir, alb_nir_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          alb_nir_dir(i) = alb_nir_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_nir_dif)) THEN
        alb_nir_lnd_ptr => rad__mem%alb_nir_lnd%ptr(ics:ice,iblk)
        !$acc data present(alb_nir_dif, alb_nir_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          alb_nir_dif(i)      = alb_nir_lnd_ptr(i)
        END DO
        !$acc end data
      END IF


      IF (PRESENT(CO2_flux)) THEN
        IF (tile%Is_process_active(CARBON_)) THEN
          cflux_npp_2_atm_ta_ptr => carbon__mem%cflux_npp_2_atm_ta%ptr(ics:ice,iblk)
          cflux_soilresp_2_atm_ta_ptr => carbon__mem%cflux_soilresp_2_atm_ta%ptr(ics:ice,iblk)
          cflux_herb_2_atm_ta_ptr => carbon__mem%cflux_herb_2_atm_ta%ptr(ics:ice,iblk)
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
        t_lwtr_ptr => seb__mem%t_lwtr%ptr(ics:ice,iblk)
        !$acc data present(t_lwtr_ptr, t_lwtr)
        !$acc parallel loop
        DO i=1,nc
          t_lwtr(i) = t_lwtr_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(qsat_lwtr)) THEN
        qsat_lwtr_ptr => seb__mem%qsat_lwtr%ptr(ics:ice,iblk)
        !$acc data present(qsat_lwtr, qsat_lwtr_ptr)
        !$acc parallel loop
        DO i=1,nc
          qsat_lwtr(i) = qsat_lwtr_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(s_lwtr)) THEN
        s_lwtr_ptr => seb__mem%s_lwtr%ptr(ics:ice,iblk)
        !$acc data present(s_lwtr, s_lwtr_ptr)
        !$acc parallel loop
        DO i=1,nc
          s_lwtr(i) = s_lwtr_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(evapo_wtr)) THEN
        evapo_wtr_ptr => hydro__mem%evapo_wtr%ptr(ics:ice,iblk)
        !$acc data present(evapo_wtr, evapo_wtr_ptr, fact_lake)
        !$acc parallel loop
        DO i=1,nc
          evapo_wtr(i) = fact_lake(i) * evapo_wtr_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(latent_hflx_wtr)) THEN
        latent_hflx_wtr_ptr => seb__mem%latent_hflx_wtr%ptr(ics:ice,iblk)
        !$acc data present(latent_hflx_wtr, latent_hflx_wtr_ptr, fact_lake)
        !$acc parallel loop
        DO i=1,nc
          latent_hflx_wtr(i)   = fact_lake(i) * latent_hflx_wtr_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(sensible_hflx_wtr)) THEN
        sensible_hflx_wtr_ptr => seb__mem%sensible_hflx_wtr%ptr(ics:ice,iblk)
        !$acc data present(sensible_hflx_wtr, sensible_hflx_wtr_ptr, fact_lake)
        !$acc parallel loop
        DO i=1,nc
          sensible_hflx_wtr(i) = fact_lake(i) * sensible_hflx_wtr_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(albedo_lwtr)) THEN
        albedo_lwtr_ptr => rad__mem%albedo_lwtr%ptr(ics:ice,iblk)
        !$acc data present(albedo_lwtr, albedo_lwtr_ptr)
        !$acc parallel loop
        DO i=1,nc
          albedo_lwtr(i) = albedo_lwtr_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(ice_fract_lake)) THEN
        fract_lice_ptr => seb__mem%fract_lice%ptr(ics:ice,iblk)
        !$acc data present(ice_fract_lake, fract_lice_ptr)
        !$acc parallel loop
        DO i=1,nc
          ice_fract_lake(i) = fract_lice_ptr(i)
        END DO
        !$acc end data
      END IF

      IF (seb__conf%l_ice_on_lakes) THEN
        IF (PRESENT(t_lice)) THEN
          t_lice_ptr => seb__mem%t_lice%ptr(ics:ice,iblk)
          !$acc data present(t_lice, t_lice_ptr)
          !$acc parallel loop
          DO i=1,nc
            t_lice(i) = t_lice_ptr(i)
          END DO
          !$acc end data
        END IF
        IF (PRESENT(qsat_lice)) THEN
          qsat_lice_ptr => seb__mem%qsat_lice%ptr(ics:ice,iblk)
          !$acc data present(qsat_lice, qsat_lice_ptr)
          !$acc parallel loop
          DO i=1,nc
            qsat_lice(i) = qsat_lice_ptr(i)
          END DO
          !$acc end data
        END IF
        IF (PRESENT(s_lice)) THEN
          s_lice_ptr => seb__mem%s_lice%ptr(ics:ice,iblk)
          !$acc data present(s_lice, s_lice_ptr)
          !$acc parallel loop
          DO i=1,nc
            s_lice(i) = s_lice_ptr(i)
          END DO
          !$acc end data
        END IF
        IF (PRESENT(evapo_ice)) THEN
          evapo_ice_ptr => hydro__mem%evapo_ice%ptr(ics:ice,iblk)
          !$acc data present(evapo_ice_ptr, evapo_ice, fact_lake)
          !$acc parallel loop
          DO i=1,nc
            evapo_ice(i) = fact_lake(i) * evapo_ice_ptr(i)
          END DO
          !$acc end data
        END IF
        IF (PRESENT(latent_hflx_ice)) THEN
          latent_hflx_ice_ptr => seb__mem%latent_hflx_ice%ptr(ics:ice,iblk)
          !$acc data present(latent_hflx_ice, latent_hflx_ice_ptr, fact_lake)
          !$acc parallel loop
          DO i=1,nc
            latent_hflx_ice(i) = fact_lake(i) * latent_hflx_ice_ptr(i)
          END DO
          !$acc end data
        END IF
        IF (PRESENT(sensible_hflx_ice)) THEN
          sensible_hflx_ice_ptr  => seb__mem%sensible_hflx_ice%ptr(ics:ice,iblk)
          !$acc data present(sensible_hflx_ice_ptr, sensible_hflx_ice, fact_lake)
          !$acc parallel loop
          DO i=1,nc
            sensible_hflx_ice(i) = fact_lake(i) * sensible_hflx_ice_ptr(i)
          END DO
          !$acc end data
        END IF
        IF (PRESENT(albedo_lice)) THEN
          albedo_lice_ptr => rad__mem%albedo_lice%ptr(ics:ice,iblk)
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
        t_filt_ptr => seb__mem%t_filt%ptr(ics:ice,iblk)
        !$acc data present(t_filt_ptr, t_srf)
        !$acc parallel loop
        DO i=1,nc
          t_srf(i) = t_filt_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(t_eff_srf)) THEN
        t_eff_ptr => seb__mem%t_eff%ptr(ics:ice,iblk)
        !$acc data present(t_eff_srf, t_eff_ptr)
        !$acc parallel loop
        DO i=1,nc
          t_eff_srf(i) = t_eff_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(qsat_srf)) THEN
        qsat_star_ptr => seb__mem%qsat_star%ptr(ics:ice,iblk)
        !$acc data present(qsat_srf, qsat_star_ptr)
        !$acc parallel loop
        DO i=1,nc
          qsat_srf(i)         = qsat_star_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(s_srf)) THEN
        s_star_ptr => seb__mem%s_star%ptr(ics:ice,iblk)
        !$acc data present(s_star_ptr,s_srf)
        !$acc parallel loop
        DO i=1,nc
          s_srf(i)            = s_star_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(fact_q_air)) THEN
        fact_q_air_ptr => turb__mem%fact_q_air%ptr(ics:ice,iblk)
        !$acc data present(fact_q_air, fact_q_air_ptr)
        !$acc parallel loop
        DO i=1,nc
          fact_q_air(i)       = fact_q_air_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(fact_qsat_srf)) THEN
        fact_qsat_srf_ptr => turb__mem%fact_qsat_srf%ptr(ics:ice,iblk)
        !$acc data present(fact_qsat_srf, fact_qsat_srf_ptr)
        !$acc parallel loop
        DO i=1,nc
          fact_qsat_srf(i) = fact_qsat_srf_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(evapopot)) THEN
        evapopot_ptr => hydro__mem%evapopot%ptr(ics:ice,iblk) ! For offline only
        !$acc data present(evapopot, evapopot_ptr)
        !$acc parallel loop
        DO i=1,nc
          evapopot(i) = evapopot_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(evapotrans)) THEN
        evapotrans_lnd_ptr => hydro__mem%evapotrans_lnd%ptr(ics:ice,iblk)
        !$acc data present(evapotrans, evapotrans_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          evapotrans(i)= evapotrans_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(latent_hflx)) THEN
        latent_hflx_lnd_ptr => seb__mem%latent_hflx_lnd%ptr(ics:ice,iblk)
        !$acc data present(latent_hflx, latent_hflx_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          latent_hflx(i) = latent_hflx_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(sensible_hflx)) THEN
        sensible_hflx_lnd_ptr => seb__mem%sensible_hflx_lnd%ptr(ics:ice,iblk)
        !$acc data present(sensible_hflx, sensible_hflx_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          sensible_hflx(i)    = sensible_hflx_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(grnd_hflx)) THEN
        forc_hflx_ptr => seb__mem%forc_hflx%ptr(ics:ice,iblk)    ! TODO Not used in atm yet,
        !$acc data present(grnd_hflx, forc_hflx_ptr)
        !$acc parallel loop
        DO i=1,nc
          grnd_hflx(i) = forc_hflx_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(grnd_hcap)) THEN
        heat_cap_old_ptr => seb__mem%heat_cap_old%ptr(ics:ice,iblk) ! and not computed for lake
        !$acc data present(grnd_hcap, heat_cap_old_ptr)
        !$acc parallel loop
        DO i=1,nc
          grnd_hcap(i)        = heat_cap_old_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(rough_h_srf)) THEN
        rough_h_ptr => turb__mem%rough_h%ptr(ics:ice,iblk)
        !$acc data present(rough_h_srf, rough_h_ptr)
        !$acc parallel loop
        DO i=1,nc
          rough_h_srf(i) = rough_h_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(rough_m_srf)) THEN
        rough_m_ptr => turb__mem%rough_m%ptr(ics:ice,iblk)
        !$acc data present(rough_m_srf, rough_m_ptr)
        !$acc parallel loop
        DO i=1,nc
          rough_m_srf(i)      = rough_m_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(q_snocpymlt)) THEN
        q_snocpymlt_ptr => hydro__mem%q_snocpymlt%ptr(ics:ice,iblk)
        !$acc data present(q_snocpymlt, q_snocpymlt_ptr)
        !$acc parallel loop
        DO i=1,nc
          q_snocpymlt(i)      = q_snocpymlt(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_vis_dir)) THEN
        alb_vis_lnd_ptr => rad__mem%alb_vis_lnd%ptr(ics:ice,iblk)  ! TODO
       !$acc data present(alb_vis_lnd_ptr, alb_vis_dir)
       !$acc parallel loop
       DO i=1,nc
         alb_vis_dir(i) = alb_vis_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_vis_dif)) THEN
        alb_vis_lnd_ptr => rad__mem%alb_vis_lnd%ptr(ics:ice,iblk)  ! TODO
        !$acc data present(alb_vis_dif, alb_vis_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          alb_vis_dif(i) = alb_vis_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_nir_dir)) THEN
        alb_nir_lnd_ptr => rad__mem%alb_nir_lnd%ptr(ics:ice,iblk)  ! TODO
        !$acc data present(alb_nir_dir, alb_nir_lnd_ptr)
        !$acc parallel loop
        DO i=1,nc
          alb_nir_dir(i) = alb_nir_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(alb_nir_dif)) THEN
        alb_nir_lnd_ptr => rad__mem%alb_nir_lnd%ptr(ics:ice,iblk)  ! TODO
        !$acc data present(alb_nir_dif)
        !$acc parallel loop
        DO i=1,nc
          alb_nir_dif(i) = alb_nir_lnd_ptr(i)
        END DO
        !$acc end data
      END IF
      IF (PRESENT(CO2_flux)) THEN
        IF (tile%Is_process_active(CARBON_)) THEN
          cflux_npp_2_atm_ta_ptr => carbon__mem%cflux_npp_2_atm_ta%ptr(ics:ice,iblk)
          cflux_soilresp_2_atm_ta_ptr => carbon__mem%cflux_soilresp_2_atm_ta%ptr(ics:ice,iblk)
          cflux_herb_2_atm_ta_ptr => carbon__mem%cflux_herb_2_atm_ta%ptr(ics:ice,iblk)
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
