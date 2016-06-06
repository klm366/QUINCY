!>
!! @brief Contains types for the atmosphere-land interface memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_atm2land_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_memory_class,  ONLY: t_jsb_memory
  USE mo_jsb_var_class,     ONLY: t_jsb_var_real2d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_atm2land_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 30

  TYPE, EXTENDS(t_jsb_memory) :: t_atm2land_memory
    TYPE(t_jsb_var_real2d) :: &
      t_air              , &
      q_air              , &
      rain               , &
      snow               , &
      wind_air           , &
      wind_10m           , &
      lw_srf_down        , &
      swvis_srf_down     , &
      swnir_srf_down     , &
      swpar_srf_down     , &
      press_srf          , &
      drag_srf           , &
      t_acoef            , &
      t_bcoef            , &
      q_acoef            , &
      q_bcoef            , &
      pch                , &
      cos_zenith_angle   , &
      ! For lakes:
      albedo_wtr         , &
      latent_hflx_wtr    , &
      sensible_hflx_wtr  , &
      latent_hflx_ice    , &
      sensible_hflx_ice
  CONTAINS
    PROCEDURE :: Init => Init_atm2land_memory
  END TYPE t_atm2land_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atm2land_memory_class'

CONTAINS

  SUBROUTINE Init_atm2land_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_atm2land_memory),  INTENT(inout), TARGET :: mem
    CHARACTER(len=*),          INTENT(in)    :: prefix
    CHARACTER(len=*),          INTENT(in)    :: suffix
    INTEGER,                   INTENT(in)    :: lct_ids(:)
    INTEGER,                   INTENT(in)    :: model_id

    TYPE(t_jsb_model),       POINTER :: model
    TYPE(t_jsb_grid),        POINTER :: hgrid
    TYPE(t_jsb_vgrid),       POINTER :: surface

    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_atm2land_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    CALL jsb_add_var( mem, 't_air', mem%t_air,       &
      & hgrid, surface,                                                      &
      & t_cf('air_temperature', 'K', ''),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'q_air', mem%q_air,                      &
      & hgrid, surface,                                                      &
      & t_cf('air_specific_humidity', '', ''),                               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'rain', mem%rain,                        &
      & hgrid, surface,                                                      &
      & t_cf('rain_fall', 'kg/m2/s', ''),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'snow', mem%snow,                        &
      & hgrid, surface,                                                      &
      & t_cf('snow_fall', 'kg/m2/s', ''),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'wind_air', mem%wind_air,                &
      & hgrid, surface,                                                      &
      & t_cf('air_wind_speed', 'm/s', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'wind_10m', mem%wind_10m,                &
      & hgrid, surface,                                                      &
      & t_cf('10m_wind_speed', 'm/s', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'lw_srf_down', mem%lw_srf_down,          &
      & hgrid, surface,                                                      &
      & t_cf('downward_longwave_radiation', 'W/m2', ''),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'swvis_srf_down', mem%swvis_srf_down,    &
      & hgrid, surface,                                                      &
      & t_cf('downward_visible_shortwave_radiation', 'W/m2', ''),            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'swnir_srf_down', mem%swnir_srf_down,    &
      & hgrid, surface,                                                      &
      & t_cf('downward_nir_shortwave_radiation', 'W/m2', ''),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'swpar_srf_down', mem%swpar_srf_down,    &
      & hgrid, surface,                                                      &
      & t_cf('downward_par_shortwave_radiation', 'W/m2', ''),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'press_srf', mem%press_srf,              &
      & hgrid, surface,                                                      &
      & t_cf('surface_pressure', '', ''),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'drag_srf', mem%drag_srf,                &
      & hgrid, surface,                                                      &
      & t_cf('surface_drag', '', ''),                                        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 't_acoef', mem%t_acoef,                  &
      & hgrid, surface,                                                      &
      & t_cf('temperature_acoef', '', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 't_bcoef', mem%t_bcoef,                  &
      & hgrid, surface,                                                      &
      & t_cf('temperature_acoef', '', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'q_acoef', mem%q_acoef,                  &
      & hgrid, surface,                                                      &
      & t_cf('specific_humidity_acoef', '', ''),                             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'q_bcoef', mem%q_bcoef,                  &
      & hgrid, surface,                                                      &
      & t_cf('specific_humidity_bcoef', '', ''),                             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'pch', mem%pch,                          &
      & hgrid, surface,                                                      &
      & t_cf('pch', '', ''),                                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'cos_zenith_angle', mem%cos_zenith_angle,&
      & hgrid, surface,                                                      &
      & t_cf('cosine_zenith_angle', '', ''),                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'albedo_wtr', mem%albedo_wtr,                     &
      & hgrid, surface,                                                      &
      & t_cf('albedo_wtr', '', 'Albedo of open lake water'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & loutput=.TRUE., lrestart=.TRUE., initval_r=0.07_wp )

    CALL jsb_add_var( mem, 'latent_hflx_wtr', mem%latent_hflx_wtr,           &
      & hgrid, surface,                                                      &
      & t_cf('latent_hflx_wtr', '', 'Latent heat flux over lake water'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'sensible_hflx_wtr', mem%sensible_hflx_wtr,       &
      & hgrid, surface,                                                      &
      & t_cf('sensible_hflx_wtr', '', 'Sensible heat flux over lake water'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'latent_hflx_ice', mem%latent_hflx_ice,           &
      & hgrid, surface,                                                      &
      & t_cf('latent_hflx_ice', '', 'Latent heat flux over lake ice'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'sensible_hflx_ice', mem%sensible_hflx_ice,       &
      & hgrid, surface,                                                      &
      & t_cf('sensible_hflx_ice', '', 'Sensible heat flux over lake ice'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

  END SUBROUTINE Init_atm2land_memory

#endif
END MODULE mo_atm2land_memory_class
