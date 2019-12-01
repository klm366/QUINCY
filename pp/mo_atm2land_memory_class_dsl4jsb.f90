!
! This file was automatically generated from:
!    /scratch/snx3000/rschnur/icon-rc-testing-plus-mc8.test/externals/jsbach/src/atmland/mo_atm2land_memory_class.f90
!    with md5 hash: 756cd68e0992cc26db2d47eea3560593
!
!>
!! @brief Contains types for the atmosphere-land interface memory.
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

MODULE mo_a2l_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  !USE mo_exception,         ONLY: message
  USE mo_util,              ONLY: One_of

  !USE mo_jsb_model_class,   ONLY: t_jsb_model
  !USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_memory_class,  ONLY: t_jsb_memory
  USE mo_jsb_var_class,     ONLY: t_jsb_var_real2d
  USE mo_jsb_lct_class,     ONLY: LAKE_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_a2l_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 35

  TYPE, EXTENDS(t_jsb_memory) :: t_a2l_memory
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
      fract_par_diffuse   , &
      press_srf          , &
      drag_srf           , &
      t_acoef            , &
      t_bcoef            , &
      q_acoef            , &
      q_bcoef            , &
      pch                , &
      cos_zenith_angle   , &
      CO2_air            , &
      CO2_air_mol        , &
      ! For lakes:
      drag_wtr           , &
      drag_ice           , &
      t_acoef_wtr        , &
      t_bcoef_wtr        , &
      q_acoef_wtr        , &
      q_bcoef_wtr        , &
      t_acoef_ice        , &
      t_bcoef_ice        , &
      q_acoef_ice        , &
      q_bcoef_ice
  CONTAINS
    PROCEDURE :: Init => Init_a2l_memory
  END TYPE t_a2l_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_a2l_memory_class'

CONTAINS

  SUBROUTINE Init_a2l_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_varlist,       ONLY: BASIC !, MEDIUM, FULL
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_a2l_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_grid),    POINTER :: hgrid
    TYPE(t_jsb_vgrid),   POINTER :: surface

    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':new_a2l_memory'
    IF (model_id > 0) CONTINUE ! avoid compiler warning about dummy argument not being used
    table = tables(1)

    hgrid => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    CALL jsb_add_var( mem, 't_air', mem%t_air,       &
      & hgrid, surface,                                                      &
      & t_cf('air_temperature', 'K', ''),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=283.15_wp )  ! R: I just invented a value to start with...

    CALL jsb_add_var( mem, 'q_air', mem%q_air,                                            &
      & hgrid, surface,                                                                   &
      & t_cf('air_specific_humidity', 'kg kg-1', 'Specific humidity of air at surface.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.FALSE.,                                                                 &
      & output_level=BASIC, &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'rain', mem%rain,                        &
      & hgrid, surface,                                                      &
      & t_cf('rain_fall', 'kg/m2/s', ''),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'snow', mem%snow,                        &
      & hgrid, surface,                                                      &
      & t_cf('snow_fall', 'kg/m2/s', ''),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'wind_air', mem%wind_air,                &
      & hgrid, surface,                                                      &
      & t_cf('air_wind_speed', 'm/s', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'wind_10m', mem%wind_10m,                &
      & hgrid, surface,                                                      &
      & t_cf('10m_wind_speed', 'm/s', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'lw_srf_down', mem%lw_srf_down,          &
      & hgrid, surface,                                                      &
      & t_cf('downward_longwave_radiation', 'W/m2', ''),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'swvis_srf_down', mem%swvis_srf_down,             &
      & hgrid, surface,                                                      &
      & t_cf('downward_visible_shortwave_radiation', 'W/m2', ''),            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'swnir_srf_down', mem%swnir_srf_down,             &
      & hgrid, surface,                                                      &
      & t_cf('downward_nir_shortwave_radiation', 'W/m2', ''),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'swpar_srf_down', mem%swpar_srf_down,    &
      & hgrid, surface,                                                      &
      & t_cf('downward_par_shortwave_radiation', 'W/m2', ''),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'fract_par_diffuse', mem%fract_par_diffuse,&
      & hgrid, surface,                                                      &
      & t_cf('fract_par_diffuse', '', ''),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'press_srf', mem%press_srf,                       &
      & hgrid, surface,                                                      &
      & t_cf('surface_pressure', 'N/m2', ''),                                & ! N/m2 =Pa
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=101325.0_wp )  ! R: I just took average atmospheric pressure
                                 !    at sea level as value to start with...

    CALL jsb_add_var( mem, 'drag_srf', mem%drag_srf,                         &
      & hgrid, surface,                                                      &
      & t_cf('surface_drag', '', ''),                                        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 't_acoef', mem%t_acoef,                           &
      & hgrid, surface,                                                      &
      & t_cf('temperature_acoef', '', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 't_bcoef', mem%t_bcoef,                           &
      & hgrid, surface,                                                      &
      & t_cf('temperature_acoef', '', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'q_acoef', mem%q_acoef,                           &
      & hgrid, surface,                                                      &
      & t_cf('specific_humidity_acoef', '', ''),                             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'q_bcoef', mem%q_bcoef,                           &
      & hgrid, surface,                                                      &
      & t_cf('specific_humidity_bcoef', '', ''),                             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'pch', mem%pch,                                   &
      & hgrid, surface,                                                      &
      & t_cf('pch', '', ''),                                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'cos_zenith_angle', mem%cos_zenith_angle,         &
      & hgrid, surface,                                                      &
      & t_cf('cosine_zenith_angle', '', ''),                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'co2_air', mem%CO2_air, &
      & hgrid, surface,                                                                            &
      & t_cf('CO2_air', '[kg(CO2)/kg(air)]', 'CO2 mass mixing ratio of lowest atmosphere level'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                         &
      & prefix, suffix,                                                                            &
      & output_level=BASIC, &
      & lrestart=.FALSE., initval_r=6.077E-4_wp ) ! R: I just freely choosed this value. Corresponds to 400 ppm:
                                                  !    CO2:44.011 g/mol; dry air: 28,97 g/mol
                                                  !    400 ppm => 6.077*10^-4

    CALL jsb_add_var( mem, 'co2_air_mol', mem%CO2_air_mol, &
      & hgrid, surface,                                                                                &
      & t_cf('CO2_air_mol', '[mol(CO2)/mol(air)]', 'CO2 mol mixing ratio of lowest atmosphere level'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=BASIC, &
      & lrestart=.FALSE., initval_r=4.0E-4_wp ) ! R: I just choosed this value for 400 ppm

    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0 .AND. .NOT. ASSOCIATED(mem%t_acoef_wtr%ptr)) THEN

      CALL jsb_add_var( mem, 'drag_wtr', mem%drag_wtr,                &
        & hgrid, surface,                                                      &
        & t_cf('surface_drag_wtr', '', ''),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL  jsb_add_var( mem, 'drag_ice', mem%drag_ice,                &
        & hgrid, surface,                                                      &
        & t_cf('surface_drag_ice', '', ''),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 't_acoef_wtr', mem%t_acoef_wtr,                   &
        & hgrid, surface,                                                      &
        & t_cf('temperature_acoef_wtr', '', ''),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 't_bcoef_wtr', mem%t_bcoef_wtr,                   &
        & hgrid, surface,                                                      &
        & t_cf('temperature_acoef_wtr', '', ''),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'q_acoef_wtr', mem%q_acoef_wtr,                   &
        & hgrid, surface,                                                      &
        & t_cf('specific_humidity_acoef_wtr', '', ''),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'q_bcoef_wtr', mem%q_bcoef_wtr,                   &
        & hgrid, surface,                                                      &
        & t_cf('specific_humidity_bcoef_wtr', '', ''),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 't_acoef_ice', mem%t_acoef_ice,                   &
        & hgrid, surface,                                                      &
        & t_cf('temperature_acoef_ice', '', ''),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 't_bcoef_ice', mem%t_bcoef_ice,                   &
        & hgrid, surface,                                                      &
        & t_cf('temperature_acoef_ice', '', ''),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'q_acoef_ice', mem%q_acoef_ice,                   &
        & hgrid, surface,                                                      &
        & t_cf('specific_humidity_acoef_ice', '', ''),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'q_bcoef_ice', mem%q_bcoef_ice,                   &
        & hgrid, surface,                                                      &
        & t_cf('specific_humidity_bcoef_ice', '', ''),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

    END IF

  END SUBROUTINE Init_a2l_memory

#endif
END MODULE mo_a2l_memory_class
