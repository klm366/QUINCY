!>
!! @brief Contains methods for the srf memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_srf_memory
#ifndef __NO_JSBACH__

  USE mo_srf_memory_class, ONLY: t_srf_memory, t_srf_phys_state, &
                                  t_srf_energy_state, t_srf_energy_flux,          &
                                  t_srf_water_state,  t_srf_water_flux
  USE mo_jsb_grid_class,    ONLY: get_nproma, get_nblks
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: new_srf_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_srf_memory'

CONTAINS

  FUNCTION new_srf_memory(varlist_name, prefix, suffix, tile) RESULT(mem)

    USE mo_jsb_varlist,       ONLY: t_var_list, jsb_get_var_list, jsb_new_var_list, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, nc_bits, t_cf, t_grib1, t_grib2, &
                                    FILETYPE_NC2, TSTEP_CONSTANT, ZAXIS_SURFACE, ZAXIS_DEPTH_BELOW_LAND, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid, new_vgrid
    USE mo_jsb_tile_class,    ONLY: t_jsb_tile
    USE mo_jsb_model_class,   ONLY: t_jsb_model, get_vgrid, register_vgrid
    USE mo_jsb_base,          ONLY: get_model

    CHARACTER(len=*),           INTENT(in) :: varlist_name
    CHARACTER(len=*),           INTENT(in) :: prefix
    CHARACTER(len=*),           INTENT(in) :: suffix
#ifdef HAVE_F2003
    TYPE(t_jsb_tile), POINTER , INTENT(in) :: tile
#else
    TYPE(t_jsb_tile), POINTER              :: tile
#endif
    TYPE(t_srf_memory), POINTER           :: mem

    TYPE(t_srf_memory), POINTER :: new
    TYPE(t_var_list),   POINTER :: this_list
    TYPE(t_jsb_model), POINTER  :: model

    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_srf_memory'

    model => get_model(tile%owner_model_id)

    ALLOCATE(new)

    table = tables(1)

    !this_list => get_var_list(TRIM(varlist_name))
    CALL jsb_get_var_list(this_list, TRIM(varlist_name))
    IF (.NOT. ASSOCIATED(this_list)) THEN
      !this_list => new_var_list(TRIM(varlist_name), patch_id=patch_id)
      ALLOCATE(this_list)
      CALL jsb_new_var_list(this_list, TRIM(varlist_name),              &
                            loutput=.TRUE., output_type=FILETYPE_NC2,   &
                            lrestart=.TRUE., restart_type=FILETYPE_NC2, &
                            patch_id=model%grid%patch%id, table=table)
    END IF

    hgrid        => model%grid

    surface => get_vgrid(model, 'surface')
    IF (.NOT. ASSOCIATED(surface)) THEN
      surface => new_vgrid('surface', ZAXIS_SURFACE, 1)
      CALL register_vgrid(model, surface)
    END IF

    IF (ASSOCIATED(tile%parent_tile)) THEN
       new%phys_state => tile%parent_tile%mem%srf%phys_state
    ELSE
       new%phys_state => new_phys_state()
    END IF

    new%energy_state => new_energy_state()
    new%energy_flux  => new_energy_flux()
    new%water_state  => new_water_state()
    new%water_flux   => new_water_flux()

    mem => new

  CONTAINS

    FUNCTION new_phys_state() RESULT(pstate)

      TYPE(t_srf_phys_state), POINTER :: pstate

      ALLOCATE(pstate)

      CALL jsb_add_var( this_list, prefix//'glacier'//suffix, pstate%glacier,                               &
        & hgrid, surface,                                                           &
        & t_cf('glacier_mask', '', '', nc_bits),                                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'oro_stddev'//suffix, pstate%oro_stddev,                         &
        & hgrid, surface,                                                           &
        & t_cf('orography_standard_deviation', '', '', nc_bits),                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'rough_m_srf'//suffix, pstate%rough_m_srf,                       &
        & hgrid, surface,                                                           &
        & t_cf('roughness_length_momentum', 'm', '', nc_bits),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'rough_h_srf'//suffix, pstate%rough_h_srf,                       &
        & hgrid, surface,                                                           &
        & t_cf('roughness_length_heat', 'm', '', nc_bits),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'albedo_srf'//suffix, pstate%albedo_srf,                         &
        & hgrid, surface,                                                           &
        & t_cf('surface_albedo_background', '', '', nc_bits),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'alb_vis_srf'//suffix, pstate%alb_vis_srf,                       &
        & hgrid, surface,                                                           &
        & t_cf('surface_albedo_visible', '', '', nc_bits),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'alb_nir_srf'//suffix, pstate%alb_nir_srf,                       &
        & hgrid, surface,                                                           &
        & t_cf('surface_albedo_nir', '', '', nc_bits),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

    END FUNCTION new_phys_state

    FUNCTION new_energy_state() RESULT(estate)

      TYPE(t_srf_energy_state), POINTER :: estate

      ALLOCATE(estate)

      CALL jsb_add_var( this_list, prefix//'t_srf'//suffix, estate%t_srf,                                   &
        & hgrid, surface,                                                           &
        & t_cf('surface_temperature', 'K', '', nc_bits),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'t_srf_old'//suffix, estate%t_srf_old,                           &
        & hgrid, surface,                                                           &
        & t_cf('surface_temperature_old', 'K', '', nc_bits),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'t_srf_unfilt'//suffix, estate%t_srf_unfilt,                     &
        & hgrid, surface,                                                           &
        & t_cf('surface_temperature_unfiltered', 'K', '', nc_bits),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'t_srf_unfilt_old'//suffix, estate%t_srf_unfilt_old,             &
        & hgrid, surface,                                                           &
        & t_cf('surface_temperature_unfiltered_old', 'K', '', nc_bits),                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'t_eff_srf'//suffix, estate%t_eff_srf,                           &
        & hgrid, surface,                                                           &
        & t_cf('surface_temperature_eff', 'K', '', nc_bits),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

!!$      CALL jsb_add_var( this_list, prefix//'t_eff_srf'//suffix, estate%t_eff_srf,                           &
!!$        & hgrid, surface,                                                           &
!!$        & t_cf('surface_temperature_eff', 'K', '', nc_bits),                                 &
!!$        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
!!$        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'s_srf'//suffix, estate%s_srf,                                   &
        & hgrid, surface,                                                           &
        & t_cf('dry_static_energy', '', '', nc_bits),                                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

!!$      CALL jsb_add_var( this_list, prefix//'s_srf_old'//suffix, estate%s_srf_old,                           &
!!$        & hgrid, surface,                                                           &
!!$        & t_cf('dry_static_energy_old', '', '', nc_bits),                                    &
!!$        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
!!$        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'qsat_srf'//suffix, estate%qsat_srf,                             &
        & hgrid, surface,                                                           &
        & t_cf('sat_specific_humidity', '', '', nc_bits),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

!!$      CALL jsb_add_var( this_list, prefix//'qsat_srf_old'//suffix, estate%qsat_srf_old,                     &
!!$        & hgrid, surface,                                                           &
!!$        & t_cf('sat_specific_humidity_old', '', '', nc_bits),                                &
!!$        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
!!$        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'tte_corr'//suffix, estate%tte_corr,                             &
        & hgrid, surface,                                                           &
        & t_cf('tte_correction_snowmelt', '', '', nc_bits),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'fact_q_air'//suffix, estate%fact_q_air,                         &
        & hgrid, surface,                                                           &
        & t_cf('factor_air_moisture', '', '', nc_bits),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'fact_qsat_srf'//suffix, estate%fact_qsat_srf,                   &
        & hgrid, surface,                                                           &
        & t_cf('factor_qsat_srf', '', '', nc_bits),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'fact_qsat_trans_srf'//suffix, estate%fact_qsat_trans_srf,       &
        & hgrid, surface,                                                           &
        & t_cf('factor_qsat_trans_srf', '', '', nc_bits),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

    END FUNCTION new_energy_state

    FUNCTION new_energy_flux() RESULT(eflux)

      TYPE(t_srf_energy_flux), POINTER :: eflux

      ALLOCATE(eflux)

      CALL jsb_add_var( this_list, prefix//'swrad_srf_down'//suffix, eflux%swrad_srf_down,                  &
        & hgrid, surface,                                                           &
        & t_cf('swrad_srf_down', '', 'Downward surface shortwave radiation', nc_bits),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'lwrad_srf_net'//suffix, eflux%lwrad_srf_net,                    &
        & hgrid, surface,                                                           &
        & t_cf('lwrad_srf_net', '', 'Net surface longwave radiation', nc_bits),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'rad_srf_net'//suffix, eflux%rad_srf_net,                        &
        & hgrid, surface,                                                           &
        & t_cf('rad_srf_net', '', 'Net surface radiation', nc_bits),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'latent_hflx'//suffix, eflux%latent_hflx,                        &
        & hgrid, surface,                                                           &
        & t_cf('latent_hflx', '', 'Latent heat flux', nc_bits),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'sensible_hflx'//suffix, eflux%sensible_hflx,                    &
        & hgrid, surface,                                                           &
        & t_cf('sensible_hflx', '', 'Sensible heat flux', nc_bits),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

    END FUNCTION new_energy_flux

    FUNCTION new_water_state() RESULT(wstate)

      TYPE(t_srf_water_state), POINTER :: wstate

      ALLOCATE(wstate)

      CALL jsb_add_var( this_list, prefix//'wfract_srf'//suffix, wstate%wfract_srf,                         &
        & hgrid, surface,                                                           &
        & t_cf('surface_wet_fraction', '', '', nc_bits),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'wsr_srf'//suffix, wstate%wsr_srf,                               &
        & hgrid, surface,                                                           &
        & t_cf('surface_skin_reservoir', 'm', '', nc_bits),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'sfract_srf'//suffix, wstate%sfract_srf,                         &
        & hgrid, surface,                                                           &
        & t_cf('surface_snow_fraction', '', '', nc_bits),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'wsn_srf'//suffix, wstate%wsn_srf,                               &
        & hgrid, surface,                                                           &
        & t_cf('surface_snow_water_equivalent', 'm', '', nc_bits),                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'wsn_glacier'//suffix, wstate%wsn_glacier,                       &
        & hgrid, surface,                                                           &
        & t_cf('glacier_snow_water_equivalent', 'm', '', nc_bits),                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

    END FUNCTION new_water_state

    FUNCTION new_water_flux() RESULT(wflux)

      TYPE(t_srf_water_flux), POINTER :: wflux

      ALLOCATE(wflux)

      CALL jsb_add_var( this_list, prefix//'evapopot'//suffix, wflux%evapopot,                                     &
        & hgrid, surface,                                                                  &
        & t_cf('surface_evaporation_potential', '', 'Potential evaporation from surface', nc_bits), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'evapotrans'//suffix, wflux%evapotrans,                                 &
        & hgrid, surface,                                                                  &
        & t_cf('surface_evapotranspiration', '', 'Evapotranspiration from surface', nc_bits),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'trans'//suffix, wflux%trans,                                           &
        & hgrid, surface,                                                                  &
        & t_cf('surface_transpiration', '', 'Transpiration from surface', nc_bits),                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & initval_r=0.0_wp )

    END FUNCTION new_water_flux

  END FUNCTION new_srf_memory

#endif
END MODULE mo_srf_memory
