!>
!! @brief Contains methods for the atmosphere-land interface memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_atmland_memory
#ifndef __NO_JSBACH__

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory, t_land2atm_memory
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: new_atm2land_memory, new_land2atm_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_memory'

CONTAINS

  FUNCTION new_atm2land_memory(varlist_name, prefix, suffix, tile) RESULT(atm2land)

    USE mo_jsb_varlist,       ONLY: t_var_list, jsb_get_var_list, jsb_new_var_list, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, nc_bits, t_cf, t_grib1, t_grib2, &
                                    FILETYPE_NC2, ZAXIS_SURFACE, tables
    USE mo_jsb_tile_class,    ONLY: t_jsb_tile
    USE mo_jsb_model_class,   ONLY: t_jsb_model, register_vgrid, get_vgrid
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid, new_vgrid
    USE mo_jsb_base,          ONLY: get_model

    CHARACTER(len=*),           INTENT(in) :: varlist_name
    CHARACTER(len=*),           INTENT(in) :: prefix
    CHARACTER(len=*),           INTENT(in) :: suffix
#ifdef HAVE_F2003
    TYPE(t_jsb_tile), POINTER , INTENT(in) :: tile
#else
    TYPE(t_jsb_tile), POINTER              :: tile
#endif
    TYPE(t_atm2land_memory), POINTER       :: atm2land

    TYPE(t_var_list),        POINTER :: this_list
    TYPE(t_jsb_grid),        POINTER :: hgrid
    TYPE(t_jsb_vgrid),       POINTER :: surface
    TYPE(t_jsb_model),       POINTER :: model

    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_atm2land_memory'

    model => get_model(tile%owner_model_id)

    ALLOCATE(atm2land)

    table = tables(1)

    !this_list => get_var_list(TRIM(varlist_name))
    CALL jsb_get_var_list(this_list, TRIM(varlist_name))
    IF (.NOT. ASSOCIATED(this_list)) THEN
      !this_list => new_var_list(TRIM(varlist_name), patch_id=patch_id)
      ALLOCATE(this_list)
      CALL jsb_new_var_list(this_list, TRIM(varlist_name),              &
                            lrestart=.TRUE., restart_type=FILETYPE_NC2, &
                            patch_id=model%grid%patch%id, table=table)
    END IF

    hgrid => model%grid
    surface => get_vgrid(model, 'surface')
    IF (.NOT. ASSOCIATED(surface)) THEN
      surface => new_vgrid('surface', ZAXIS_SURFACE, 1)
      CALL register_vgrid(model, surface)
    END IF

    CALL jsb_add_var( this_list, prefix//'t_air'//suffix, atm2land%t_air,       &
         & hgrid, surface,                                                      &
         & t_cf('air_temperature', 'K', '', nc_bits),                           &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'q_air'//suffix, atm2land%q_air,                      &
         & hgrid, surface,                                                      &
         & t_cf('air_specific_humidity', '', '', nc_bits),                      &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'rain'//suffix, atm2land%rain,                        &
         & hgrid, surface,                                                      &
         & t_cf('rain_fall', 'kg/m2/s', '', nc_bits),                           &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'snow'//suffix, atm2land%snow,                        &
         & hgrid, surface,                                                      &
         & t_cf('snow_fall', 'kg/m2/s', '', nc_bits),                           &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'wind_air'//suffix, atm2land%wind_air,                &
         & hgrid, surface,                                                      &
         & t_cf('air_wind_speed', 'm/s', '', nc_bits),                          &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'wind_10m'//suffix, atm2land%wind_10m,                &
         & hgrid, surface,                                                      &
         & t_cf('10m_wind_speed', 'm/s', '', nc_bits),                          &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'lwrad_srf_down'//suffix, atm2land%lwrad_srf_down,    &
         & hgrid, surface,                                                      &
         & t_cf('downward_longwave_radiation', '', '', nc_bits),                &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'swrad_srf_net'//suffix, atm2land%swrad_srf_net,      &
         & hgrid, surface,                                                      &
         & t_cf('net', 'K', '', nc_bits),                                       &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'press_srf'//suffix, atm2land%press_srf,              &
         & hgrid, surface,                                                      &
         & t_cf('surface_pressure', '', '', nc_bits),                           &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'drag_srf'//suffix, atm2land%drag_srf,                &
         & hgrid, surface,                                                      &
         & t_cf('surface_drag', '', '', nc_bits),                               &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'t_acoef'//suffix, atm2land%t_acoef,                  &
         & hgrid, surface,                                                      &
         & t_cf('temperature_acoef', '', '', nc_bits),                          &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'t_bcoef'//suffix, atm2land%t_bcoef,                  &
         & hgrid, surface,                                                      &
         & t_cf('temperature_acoef', '', '', nc_bits),                          &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'q_acoef'//suffix, atm2land%q_acoef,                  &
         & hgrid, surface,                                                      &
         & t_cf('specific_humidity_acoef', '', '', nc_bits),                    &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'q_bcoef'//suffix, atm2land%q_bcoef,                  &
         & hgrid, surface,                                                      &
         & t_cf('specific_humidity_bcoef', '', '', nc_bits),                    &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'pch'//suffix, atm2land%pch,                          &
         & hgrid, surface,                                                      &
         & t_cf('pch', '', '', nc_bits),                                        &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'cos_zenith_angle'//suffix, atm2land%cos_zenith_angle,&
         & hgrid, surface,                                                      &
         & t_cf('cosine_zenith_angle', '', '', nc_bits),                        &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
         & initval_r=0.0_wp )

  END FUNCTION new_atm2land_memory

  FUNCTION new_land2atm_memory(varlist_name, prefix, suffix, tile) RESULT(land2atm)

    USE mo_jsb_varlist,       ONLY: t_var_list, jsb_get_var_list, jsb_new_var_list, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, nc_bits, t_cf, t_grib1, t_grib2, &
                                    FILETYPE_NC2, ZAXIS_SURFACE, tables
    USE mo_jsb_tile_class,    ONLY: t_jsb_tile
    USE mo_jsb_model_class,   ONLY: t_jsb_model, register_vgrid, get_vgrid
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid, new_vgrid
    USE mo_jsb_base,          ONLY: get_model

    CHARACTER(len=*),           INTENT(in) :: varlist_name
    CHARACTER(len=*),           INTENT(in) :: prefix
    CHARACTER(len=*),           INTENT(in) :: suffix
#ifdef HAVE_F2003
    TYPE(t_jsb_tile), POINTER , INTENT(in) :: tile
#else
    TYPE(t_jsb_tile), POINTER              :: tile
#endif
    TYPE(t_land2atm_memory), POINTER       :: land2atm

    TYPE(t_var_list),        POINTER :: this_list
    TYPE(t_jsb_grid),        POINTER :: hgrid
    TYPE(t_jsb_vgrid),       POINTER :: surface
    TYPE(t_jsb_model),       POINTER :: model

    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_land2atm_memory'

    model => get_model(tile%owner_model_id)

    ALLOCATE(land2atm)

    table = tables(1)

    !this_list => get_var_list(TRIM(varlist_name))
    CALL jsb_get_var_list(this_list, TRIM(varlist_name))
    IF (.NOT. ASSOCIATED(this_list)) THEN
      !this_list => new_var_list(TRIM(varlist_name), patch_id=patch_id)
      ALLOCATE(this_list)
      CALL jsb_new_var_list(this_list, TRIM(varlist_name),              &
                            lrestart=.TRUE., restart_type=FILETYPE_NC2, &
                            patch_id=model%grid%patch%id)
    END IF

    hgrid => model%grid
    surface => get_vgrid(model, 'surface')
    IF (.NOT. ASSOCIATED(surface)) THEN
      surface => new_vgrid('surface', ZAXIS_SURFACE, 1)
      CALL register_vgrid(model, surface)
    END IF

    CALL jsb_add_var( this_list, prefix//'t_srf'//suffix, land2atm%t_srf,                      &
      & hgrid, surface,                                                      &
      & t_cf('surface_temperature', 'K', '', nc_bits),                          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'t_eff_srf'//suffix, land2atm%t_eff_srf,              &
      & hgrid, surface,                                                      &
      & t_cf('effective_surface_temperature', 'K', '', nc_bits),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'qsat_srf'//suffix, land2atm%qsat_srf,                &
      & hgrid, surface,                                                      &
      & t_cf('saturated_surface_specific_humidity', '', '', nc_bits),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'fact_q_air'//suffix, land2atm%fact_q_air,            &
      & hgrid, surface,                                                      &
      & t_cf('factor_air_moisture', '', '', nc_bits),                           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'fact_qsat_srf'//suffix, land2atm%fact_qsat_srf,      &
      & hgrid, surface,                                                      &
      & t_cf('factor_qsat_srf', '', '', nc_bits),                               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'zh_srf'//suffix, land2atm%zh_srf,                    &
      & hgrid, surface,                                                      &
      & t_cf('surface_roughness_humidity', '', '', nc_bits),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'zm_srf'//suffix, land2atm%zm_srf,                    &
      & hgrid, surface,                                                      &
      & t_cf('surface_roughness_momentum', '', '', nc_bits),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'alb_vis_dif'//suffix, land2atm%alb_vis_dif,          &
      & hgrid, surface,                                                      &
      & t_cf('albedo_visible_diffuse', '', '', nc_bits),                        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'alb_nir_dif'//suffix, land2atm%alb_nir_dif,          &
      & hgrid, surface,                                                      &
      & t_cf('albedo_near_infrared_diffuse', '', '', nc_bits),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'alb_vis_dir'//suffix, land2atm%alb_vis_dir,          &
      & hgrid, surface,                                                      &
      & t_cf('albedo_visible_direct', '', '', nc_bits),                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & initval_r=0.0_wp )

    CALL jsb_add_var( this_list, prefix//'alb_nir_dir'//suffix, land2atm%alb_nir_dir,          &
      & hgrid, surface,                                                      &
      & t_cf('albedo_near_infrared_dir', '', '', nc_bits),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & initval_r=0.0_wp )

  END FUNCTION new_land2atm_memory

#endif
END MODULE mo_atmland_memory
