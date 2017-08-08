!>
!! @brief Contains types for the atmosphere-land interface memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! Integrated into SLM by Jan Engel (20160823) using JSBACH4:git
!!
MODULE mo_a2l_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later

! jsbach4
!  USE mo_util,              ONLY: One_of

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: Get_model
!  USE mo_jsb_memory_class,  ONLY: t_jsb_memory
!  USE mo_jsb_var_class,     ONLY: t_jsb_var_real2d
!  USE mo_jsb_lct_class,     ONLY: LAKE_TYPE

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_a2l_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

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
                              swpar_fdiffuse     , &  ! jsbach4: frac_par_diffuse
                              press_srf          , &
                              drag_srf           , &
                              t_acoef            , &
                              t_bcoef            , &
                              q_acoef            , &
                              q_bcoef            , &
                              pch                , &
                              cos_zenith_angle   , &
!                              CO2_air            , &   ! jsbach4 variables, replaced by co2_mixing_ratio  _c13  _c14
!                              CO2_air_mol        , &   !
                              co2_mixing_ratio     , &  ! from SLM, replacing CO2_air and CO2_air_mol from jsbach4
                              co2_mixing_ratio_C13 , &  ! from SLM, replacing CO2_air and CO2_air_mol from jsbach4
                              co2_mixing_ratio_C14 , &  ! from SLM, replacing CO2_air and CO2_air_mol from jsbach4
                              ga                   , &  !> @todo ga in t_a2l_memory: is a variable from SLM, which may be moved elsewhere
                              nhx_deposition     , &    ! mumol/m2/sec  - new with r507 ticket #175
                              noy_deposition     , &    ! mumol/m2/sec  - new with r507 ticket #175
                              p_deposition       , &    ! mumol/m2/sec  - new with r... ticket #261
                              ! For lakes:
                              albedo_lwtr        , &
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

    TYPE(t_jsb_var_real2d) :: &
                              prev_one_cos_zenith_angle   , &   ! cos_zenith_angle of the previous timestep
                              prev_two_cos_zenith_angle         ! cos_zenith_angle of the before previous timestep


  CONTAINS
    PROCEDURE :: Init => Init_a2l_memory
  END TYPE t_a2l_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_a2l_memory_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize memory (variables) for: a2l
  !!
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_a2l_memory(mem, prefix, suffix, lct_ids, model_id) ! jsbach4

! jsbach4
!    USE mo_jsb_varlist,       ONLY: jsb_add_var
!    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
!    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    ! needed to calc SLM init values:
    USE mo_atmland_constants, ONLY: def_wind, def_co2_deltaC13, def_co2_deltaC14, def_t_air, standard_press_srf, def_vpd, &
                                    eps_vpd, def_co2_mixing_ratio, def_swdown, frac_vis_swtot_srf, frac_par_swvis_srf, &
                                    def_fdiffuse, def_cos_angle
    USE mo_shared_isotope,    ONLY: calc_mixing_ratio_C13C12, calc_mixing_ratio_C14C
    USE mo_atmland_util,      ONLY: calc_spec_humidity_sat

    ! 0.1 InOut
    CLASS(t_a2l_memory),   INTENT(inout), TARGET :: mem
    CHARACTER(len=*),           INTENT(in)    :: prefix
    CHARACTER(len=*),           INTENT(in)    :: suffix
    INTEGER,                    INTENT(in)    :: lct_ids
!    INTEGER,                    INTENT(in)    :: lct_ids(:)  ! jsbach4
    INTEGER,                   INTENT(in)    :: model_id
    ! 0.2 Local
    REAL(wp) :: C13C12_mixing_ratio, &
                C14C_mixing_ratio
    REAL(wp) :: wind_air_init, &
                wind_10m_init, &
                drag_srf_init, &
                ga_init, &
                q_air_init, &
                t_air_init, &
                press_srf_init, &
                swvis_srf_down_init, &
                swnir_srf_down_init, &
                swpar_srf_down_init, &
                swpar_fdiffuse_init, &
                cos_zenith_angle_init, &
                co2_mixing_ratio_init, &
                co2_mixing_ratio_c13_init, &
                co2_mixing_ratio_c14_init


! jsbach4
    TYPE(t_jsb_model),       POINTER :: model 
!    TYPE(t_jsb_grid),        POINTER :: hgrid 
!    TYPE(t_jsb_vgrid),       POINTER :: surface  
!    INTEGER :: table 

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_a2l_memory'

! jsbach4
    model => Get_model(model_id) 
!    table = tables(1)  
!    hgrid => Get_grid(mem%grid_id) 
!    surface => Get_vgrid('surface') 


    !! SLM specific: initially calculate some init values from constants (also using routines from mo_shared_isotope)
    !! 
    !! 
    C13C12_mixing_ratio   = calc_mixing_ratio_C13C12(def_co2_deltaC13) 
    C14C_mixing_ratio     = calc_mixing_ratio_C14C(def_co2_deltaC13,def_co2_deltaC14)

    t_air_init            = def_t_air
    press_srf_init        = standard_press_srf
    q_air_init            = calc_spec_humidity_sat(t_air_init,press_srf_init) - def_vpd / press_srf_init * eps_vpd

    wind_air_init         = def_wind
    wind_10m_init         = wind_air_init
    drag_srf_init         = 1000.0_wp
    ga_init               = 40._wp / drag_srf_init * wind_10m_init

    swvis_srf_down_init   = def_swdown * frac_vis_swtot_srf
    swnir_srf_down_init   = def_swdown * (1._wp - frac_vis_swtot_srf)
    swpar_srf_down_init   = swvis_srf_down_init * frac_par_swvis_srf
    swpar_fdiffuse_init   = def_fdiffuse
    cos_zenith_angle_init = def_cos_angle

    co2_mixing_ratio_init       = def_co2_mixing_ratio
    co2_mixing_ratio_c13_init   = co2_mixing_ratio_init / (1._wp + 1._wp / C13C12_mixing_ratio)
    co2_mixing_ratio_c14_init   = C14C_mixing_ratio * co2_mixing_ratio_init 


    !! allocate and init variables
    !! 
    !! 
    CALL slm_add_var( mem, 't_air', mem%t_air,       &
! changed throughout this module
!    CALL jsb_add_var( mem, 't_air', mem%t_air,       &
!      & hgrid, surface,                                                      &
!      & t_cf('air_temperature', 'K', ''),                                    &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r= t_air_init )

    CALL slm_add_var( mem, 'q_air', mem%q_air,                      &
!      & hgrid, surface,                                                      &
!      & t_cf('air_specific_humidity', '', ''),                               &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r= q_air_init )

    CALL slm_add_var( mem, 'rain', mem%rain,                        &
!      & hgrid, surface,                                                      &
!      & t_cf('rain_fall', 'kg/m2/s', ''),                                    &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL slm_add_var( mem, 'snow', mem%snow,                        &
!      & hgrid, surface,                                                      &
!      & t_cf('snow_fall', 'kg/m2/s', ''),                                    &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL slm_add_var( mem, 'wind_air', mem%wind_air,                &
!      & hgrid, surface,                                                      &
!      & t_cf('air_wind_speed', 'm/s', ''),                                   &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r= wind_air_init )

    CALL slm_add_var( mem, 'wind_10m', mem%wind_10m,                &
!      & hgrid, surface,                                                      &
!      & t_cf('10m_wind_speed', 'm/s', ''),                                   &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r= wind_10m_init )

    CALL slm_add_var( mem, 'lw_srf_down', mem%lw_srf_down,          &
!      & hgrid, surface,                                                      &
!      & t_cf('downward_longwave_radiation', 'W/m2', ''),                     &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL slm_add_var( mem, 'swvis_srf_down', mem%swvis_srf_down,    &
!      & hgrid, surface,                                                      &
!      & t_cf('downward_visible_shortwave_radiation', 'W/m2', ''),            &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r= swvis_srf_down_init )

    CALL slm_add_var( mem, 'swnir_srf_down', mem%swnir_srf_down,    &
!      & hgrid, surface,                                                      &
!      & t_cf('downward_nir_shortwave_radiation', 'W/m2', ''),                &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r= swnir_srf_down_init )

    CALL slm_add_var( mem, 'swpar_srf_down', mem%swpar_srf_down,    &
!      & hgrid, surface,                                                      &
!      & t_cf('downward_par_shortwave_radiation', 'W/m2', ''),                &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r= swpar_srf_down_init )

    CALL slm_add_var( mem, 'swpar_fdiffuse', mem%swpar_fdiffuse,&
!    CALL slm_add_var( mem, 'frac_par_diffuse', mem%frac_par_diffuse,&  ! jsbach4
!      & hgrid, surface,                                                      &
!      & t_cf('frac_par_diffuse', '', ''),                                    &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r= swpar_fdiffuse_init )

    CALL slm_add_var( mem, 'press_srf', mem%press_srf,              &
!      & hgrid, surface,                                                      &
!      & t_cf('surface_pressure', '', ''),                                    &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r= press_srf_init )

    CALL slm_add_var( mem, 'drag_srf', mem%drag_srf,                &
!      & hgrid, surface,                                                      &
!      & t_cf('surface_drag', '', ''),                                        &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r= drag_srf_init )

    CALL slm_add_var( mem, 't_acoef', mem%t_acoef,                  &
!      & hgrid, surface,                                                      &
!      & t_cf('temperature_acoef', '', ''),                                   &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL slm_add_var( mem, 't_bcoef', mem%t_bcoef,                  &
!      & hgrid, surface,                                                      &
!      & t_cf('temperature_acoef', '', ''),                                   &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL slm_add_var( mem, 'q_acoef', mem%q_acoef,                  &
!      & hgrid, surface,                                                      &
!      & t_cf('specific_humidity_acoef', '', ''),                             &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL slm_add_var( mem, 'q_bcoef', mem%q_bcoef,                  &
!      & hgrid, surface,                                                      &
!      & t_cf('specific_humidity_bcoef', '', ''),                             &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL slm_add_var( mem, 'pch', mem%pch,                          &
!      & hgrid, surface,                                                      &
!      & t_cf('pch', '', ''),                                                 &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL slm_add_var( mem, 'cos_zenith_angle', mem%cos_zenith_angle,&
!      & hgrid, surface,                                                      &
!      & t_cf('cosine_zenith_angle', '', ''),                                 &
!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!      & prefix, suffix,                                                      &
      & initval_r= cos_zenith_angle_init )

!!    CALL slm_add_var( mem, 'CO2_air', mem%CO2_air, &
!!!      & hgrid, surface,                                                                                &
!!!      & t_cf('CO2_air', '[kg(CO2)/kg(air)]', 'CO2 mass mixing ratio of lowest atmosphere level'),      &
!!!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
!!!      & prefix, suffix,                                                                                &
!!!      & lrestart=.FALSE., &
!!      & initval_r=0.0_wp )

!!    CALL slm_add_var( mem, 'CO2_air_mol', mem%CO2_air_mol, &
!!!      & hgrid, surface,                                                                                &
!!!      & t_cf('CO2_air_mol', '[mol(CO2)/mol(air)]', 'CO2 mol mixing ratio of lowest atmosphere level'), &
!!!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
!!!      & prefix, suffix,                                                                                &
!!!      & lrestart=.FALSE., &
!!      & initval_r=0.0_wp )

    CALL slm_add_var( mem, 'co2_mixing_ratio', mem%co2_mixing_ratio,&
            & initval_r= co2_mixing_ratio_init )

    CALL slm_add_var( mem, 'co2_mixing_ratio_c13', mem%co2_mixing_ratio_c13,&
            & initval_r= co2_mixing_ratio_c13_init )

    CALL slm_add_var( mem, 'co2_mixing_ratio_c14', mem%co2_mixing_ratio_c14,&
            & initval_r= co2_mixing_ratio_c14_init )

    CALL slm_add_var( mem, 'ga', mem%ga,&
            & initval_r= ga_init )

    CALL slm_add_var( mem, 'nhx_deposition', mem%nhx_deposition,&
            & initval_r=0.0_wp )

    CALL slm_add_var( mem, 'noy_deposition', mem%noy_deposition,&
            & initval_r=0.0_wp )

    CALL slm_add_var( mem, 'p_deposition', mem%p_deposition,&
            & initval_r=0.0_wp )

!    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0 .AND. .NOT. ASSOCIATED(mem%t_acoef_wtr%ptr)) THEN

      CALL slm_add_var( mem, 'drag_wtr', mem%drag_wtr,                &
!        & hgrid, surface,                                                      &
!        & t_cf('surface_drag_wtr', '', ''),                                    &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )

      CALL  slm_add_var( mem, 'drag_ice', mem%drag_ice,                &
!        & hgrid, surface,                                                      &
!        & t_cf('surface_drag_ice', '', ''),                                    &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )

      CALL slm_add_var( mem, 't_acoef_wtr', mem%t_acoef_wtr,                   &
!        & hgrid, surface,                                                      &
!        & t_cf('temperature_acoef_wtr', '', ''),                               &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )

      CALL slm_add_var( mem, 't_bcoef_wtr', mem%t_bcoef_wtr,                   &
!        & hgrid, surface,                                                      &
!        & t_cf('temperature_acoef_wtr', '', ''),                               &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )

      CALL slm_add_var( mem, 'q_acoef_wtr', mem%q_acoef_wtr,                   &
!        & hgrid, surface,                                                      &
!        & t_cf('specific_humidity_acoef_wtr', '', ''),                         &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )

      CALL slm_add_var( mem, 'q_bcoef_wtr', mem%q_bcoef_wtr,                   &
!        & hgrid, surface,                                                      &
!        & t_cf('specific_humidity_bcoef_wtr', '', ''),                         &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )

      CALL slm_add_var( mem, 't_acoef_ice', mem%t_acoef_ice,                   &
!        & hgrid, surface,                                                      &
!        & t_cf('temperature_acoef_ice', '', ''),                               &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )

      CALL slm_add_var( mem, 't_bcoef_ice', mem%t_bcoef_ice,                   &
!        & hgrid, surface,                                                      &
!        & t_cf('temperature_acoef_ice', '', ''),                               &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )

      CALL slm_add_var( mem, 'q_acoef_ice', mem%q_acoef_ice,                   &
!        & hgrid, surface,                                                      &
!        & t_cf('specific_humidity_acoef_ice', '', ''),                         &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )

      CALL slm_add_var( mem, 'q_bcoef_ice', mem%q_bcoef_ice,                   &
!        & hgrid, surface,                                                      &
!        & t_cf('specific_humidity_bcoef_ice', '', ''),                         &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )

      CALL slm_add_var( mem, 'albedo_lwtr', mem%albedo_lwtr,                   &
!        & hgrid, surface,                                                      &
!        & t_cf('albedo_lwtr', '', 'Albedo of open lake water'),                &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!        & prefix, suffix,                                                      &
!        & loutput=.TRUE., lrestart=.TRUE., &
        & initval_r=0.07_wp )

!    END IF

      CALL slm_add_var( mem, 'prev_one_cos_zenith_angle', mem%prev_one_cos_zenith_angle,                   &
!        & hgrid, surface,                                                      &
!        & t_cf('prev_one_cos_zenith_angle', '', 'Albedo of open lake water'),                &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!        & prefix, suffix,                                                      &
!        & loutput=.TRUE., lrestart=.TRUE., &
        & initval_r=0.0_wp )

      CALL slm_add_var( mem, 'prev_two_cos_zenith_angle', mem%prev_two_cos_zenith_angle,                   &
!        & hgrid, surface,                                                      &
!        & t_cf('prev_two_cos_zenith_angle', '', 'Albedo of open lake water'),                &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!        & prefix, suffix,                                                      &
!        & loutput=.TRUE., lrestart=.TRUE., &
        & initval_r=0.0_wp )


  END SUBROUTINE Init_a2l_memory

#endif
END MODULE mo_a2l_memory_class
