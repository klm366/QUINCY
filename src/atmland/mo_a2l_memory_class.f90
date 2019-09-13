!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
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
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d              ! *slmjsb4

! jsbach4 needed for lakes
!  USE mo_util,              ONLY: One_of
!  USE mo_jsb_lct_class,     ONLY: LAKE_TYPE


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_a2l_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 60


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for a2l memory
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_a2l_memory

    TYPE(t_jsb_var_real2d) :: &
                              t_air                 , &  !> 2m air temperature [K] 
                              q_air                 , &  !> 2m air specific humidity [g g-1] 
                              rain                  , &  !> liquid water precipitation [kg m-2 s-1]
                              snow                  , &  !> solid water precipitation [kg m-2 s-1]
                              wind_air              , &  !> 2m air wind speed [m s-1] 
                              wind_10m              , &  !> 10m air wind speed [m s-1] 
                              lw_srf_down           , &  !> longwave surface downward flux [W m-2] 
                              swvis_srf_down        , &  !> shortwave surface downward flux in the visible band [W m-2]
                              swnir_srf_down        , &  !> shortwave surface downward flux in the near infrared band [W m-2] 
                              swpar_srf_down        , &  !> shortwave surface downward flux in the PAR band [W m-2] 
                              swpar_fdiffuse        , &  !> fraction of shortwave surface downward flux that is diffuse [fraction]
                              press_srf             , &  !> surface air pressure [Pa]
                              drag_srf              , &  !> total surface drag [unitless] 
                              t_acoef               , &  !> Richtmeyer-Morten coefficient A for temperature [unitless] 
                              t_bcoef               , &  !> Richtmeyer-Morten coefficient B for temperature [unitless]
                              q_acoef               , &  !> Richtmeyer-Morten coefficient A for specific humidity [unitless] 
                              q_bcoef               , &  !> Richtmeyer-Morten coefficient B for specific humidity [unitless]
                              pch                   , &  !> total surface drag [unitless] 
                              cos_zenith_angle           !> Cosine of solar zenith angle [unitless]

    ! CO2
    TYPE(t_jsb_var_real2d) :: &
                              ! jsbach4 variables, replaced by co2_mixing_ratio  _c13  _c14
                              ! CO2_air               , &  !> 
                              ! CO2_air_mol           , &  !> 
                              
                              ! from SLM, replacing CO2_air and CO2_air_mol from jsbach4
                              co2_mixing_ratio      , &  !> surface air co2 mixing ratio [ppm] 
                              co2_mixing_ratio_C13  , &  !> surface air 13co2 mixing ratio [ppm]
                              co2_mixing_ratio_C14       !> surface air 14co2 mixing ratio [scaled]
                              
    ! ga    @todo  ga in t_a2l_memory is a variable from SLM, which may be moved elsewhere
    TYPE(t_jsb_var_real2d) :: &
                              ga                         !> aerodynamic conductance [m s-1]

    ! Deposition
    TYPE(t_jsb_var_real2d) :: &                                                            
                              nhx_deposition        , &  !> surface downward reduced nitrogen deposition velocity [mumol m-2 sec-1]
                              noy_deposition        , &  !> surface downward oxidised nitrogen deposition velocity [mumol m-2 sec-1]
                              nhx_n15_deposition    , &  !> surface downward reduced nitrogen-15 deposition velocity [mumol m-2 sec-1]
                              noy_n15_deposition    , &  !> surface downward oxidised nitrogen-15 deposition velocity [mumol m-2 sec-1]
                              p_deposition               !> surface downward reduced phosphoruns deposition velocity [mumol m-2 sec-1]

    ! For lakes
    TYPE(t_jsb_var_real2d) :: &                              
                              albedo_lwtr           , &  !> lake albedo [unitless]
                              drag_wtr              , &  !> drag coefficient over water [unitless]
                              drag_ice              , &  !> drag coefficient over ice [unitless]
                              t_acoef_wtr           , &  !> Richtmeyer-Morten coefficient A for temperature over water [unitless]
                              t_bcoef_wtr           , &  !> Richtmeyer-Morten coefficient A for temperature over water [unitless]
                              q_acoef_wtr           , &  !> Richtmeyer-Morten coefficient A for specific humidity over water [unitless]
                              q_bcoef_wtr           , &  !> Richtmeyer-Morten coefficient A for specific humidity over water [unitless]
                              t_acoef_ice           , &  !> Richtmeyer-Morten coefficient A for temperature over ice [unitless]
                              t_bcoef_ice           , &  !> Richtmeyer-Morten coefficient A for temperature over ice [unitless]
                              q_acoef_ice           , &  !> Richtmeyer-Morten coefficient A for specific humidity over ice [unitless]
                              q_bcoef_ice                !> Richtmeyer-Morten coefficient A for specific humidity over ice [unitless]

    ! For calculation of the timestep after local midnight (1st timestep of the new day)
    TYPE(t_jsb_var_real2d) :: &
                              prev_one_cos_zenith_angle   , &   !> cos_zenith_angle of the previous timestep [unitless]
                              prev_two_cos_zenith_angle         !> cos_zenith_angle of the before previous timestep [unitless]

  CONTAINS
    PROCEDURE :: Init => Init_a2l_memory
  END TYPE t_a2l_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_a2l_memory_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize a2l process memory variables
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_a2l_memory(mem, prefix, suffix, lct_ids, model_id) ! jsbach4

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_slm_varlist,       ONLY: jsb_add_var   ! replace by mo_jsb_varlist later    *quincytemp
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2   ! , TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    ! needed to calc SLM init values:
    USE mo_atmland_constants, ONLY: def_wind, def_co2_deltaC13, def_co2_deltaC14, def_t_air, standard_press_srf, def_vpd, &
                                    eps_vpd, def_co2_mixing_ratio, def_swdown, frac_vis_swtot_srf, frac_par_swvis_srf, &
                                    def_fdiffuse, def_cos_angle
    USE mo_isotope_util,      ONLY: calc_mixing_ratio_C13C12, calc_mixing_ratio_C14C
    USE mo_atmland_util,      ONLY: calc_spec_humidity_sat

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_a2l_memory),  INTENT(inout), TARGET :: mem             !> a2l memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ---------------------------
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
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_surface                ! Vertical grid
    INTEGER                    :: table
    CHARACTER(len=*),PARAMETER :: routine = TRIM(modname)//':Init_a2l_memory'
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(model_id)
    table         = 1     ! in jsbach4 this line is: table = tables(1)     *quincytemp
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface => Get_vgrid('surface')


    !> 0.9 SLM specific: initially calculate some init values from constants (also using routines from mo_isotope_util)
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


    !> 1.0  init memory
    !! 
    CALL jsb_add_var( mem, 't_air', mem%t_air, &
            hgrid, vgrid_surface, &
            t_cf('t_air', 'K', '2m air temperature'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= t_air_init )

    CALL jsb_add_var( mem, 'q_air', mem%q_air, &
            hgrid, vgrid_surface, &
            t_cf('q_air', 'g g-1', '2m air specific humidity'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= q_air_init )

    CALL jsb_add_var( mem, 'rain', mem%rain, &
            hgrid, vgrid_surface, &
            t_cf('rain', 'kg m-2 s-1', 'liquid water precipitation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'snow', mem%snow, &
            hgrid, vgrid_surface, &
            t_cf('snow', 'kg m-2 s-1', 'solid water precipitation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'wind_air', mem%wind_air, &
            hgrid, vgrid_surface, &
            t_cf('wind_air', 'm s-1', '2m air wind speed'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= wind_air_init )

    CALL jsb_add_var( mem, 'wind_10m', mem%wind_10m, &
            hgrid, vgrid_surface, &
            t_cf('wind_10m', 'm s-1', '10m air wind speed'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= wind_10m_init )

    CALL jsb_add_var( mem, 'lw_srf_down', mem%lw_srf_down, &
            hgrid, vgrid_surface, &
            t_cf('lw_srf_down', 'W m-2', 'longwave surface downward flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'swvis_srf_down', mem%swvis_srf_down, &
            hgrid, vgrid_surface, &
            t_cf('swvis_srf_down', 'W m-2', 'shortwave surface downward flux in the visible band'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= swvis_srf_down_init )

    CALL jsb_add_var( mem, 'swnir_srf_down', mem%swnir_srf_down, &
            hgrid, vgrid_surface, &
            t_cf('swnir_srf_down', 'W m-2', 'shortwave surface downward flux in the near infrared band'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= swnir_srf_down_init )

    CALL jsb_add_var( mem, 'swpar_srf_down', mem%swpar_srf_down, &
            hgrid, vgrid_surface, &
            t_cf('swpar_srf_down', 'W m-2', 'shortwave surface downward flux in the PAR band'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= swpar_srf_down_init )

    CALL jsb_add_var( mem, 'swpar_fdiffuse', mem%swpar_fdiffuse, &
            hgrid, vgrid_surface, &
            t_cf('swpar_fdiffuse', 'fraction', 'fraction of shortwave surface downward flux that is diffuse'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= swpar_fdiffuse_init )

    CALL jsb_add_var( mem, 'press_srf', mem%press_srf, &
            hgrid, vgrid_surface, &
            t_cf('press_srf', 'Pa', 'surface air pressure'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= press_srf_init )

    CALL jsb_add_var( mem, 'drag_srf', mem%drag_srf, &
            hgrid, vgrid_surface, &
            t_cf('drag_srf', 'unitless', 'total surface drag'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= drag_srf_init )

    CALL jsb_add_var( mem, 't_acoef', mem%t_acoef, &
            hgrid, vgrid_surface, &
            t_cf('t_acoef', 'unitless', 'Richtmeyer-Morten coefficient A for temperature'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 't_bcoef', mem%t_bcoef, &
            hgrid, vgrid_surface, &
            t_cf('t_bcoef', 'unitless', 'Richtmeyer-Morten coefficient B for temperature'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'q_acoef', mem%q_acoef, &
            hgrid, vgrid_surface, &
            t_cf('q_acoef', 'unitless', 'Richtmeyer-Morten coefficient A for specific humidity'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'q_bcoef', mem%q_bcoef, &
            hgrid, vgrid_surface, &
            t_cf('q_bcoef', 'unitless', 'Richtmeyer-Morten coefficient B for specific humidity'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'pch', mem%pch, &
            hgrid, vgrid_surface, &
            t_cf('pch', 'unitless', 'total surface drag'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'cos_zenith_angle', mem%cos_zenith_angle, &
            hgrid, vgrid_surface, &
            t_cf('cos_zenith_angle', 'unitless', 'Cosine of solar zenith angle'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= cos_zenith_angle_init )

!!    CALL jsb_add_var( mem, 'CO2_air', mem%CO2_air, &
!!!      & hgrid, vgrid_surface,                                                                                &
!!!      & t_cf('CO2_air', '[kg(CO2)/kg(air)]', 'CO2 mass mixing ratio of lowest atmosphere level'),      &
!!!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
!!!      & prefix, suffix,                                                                                &
!!!      & lrestart=.FALSE., &
!!      & initval_r=0.0_wp )

!!    CALL jsb_add_var( mem, 'CO2_air_mol', mem%CO2_air_mol, &
!!!      & hgrid, vgrid_surface,                                                                                &
!!!      & t_cf('CO2_air_mol', '[mol(CO2)/mol(air)]', 'CO2 mol mixing ratio of lowest atmosphere level'), &
!!!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
!!!      & prefix, suffix,                                                                                &
!!!      & lrestart=.FALSE., &
!!      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'co2_mixing_ratio', mem%co2_mixing_ratio, &
            hgrid, vgrid_surface, &
            t_cf('co2_mixing_ratio', 'ppm', 'surface air co2 mixing ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= co2_mixing_ratio_init )

    CALL jsb_add_var( mem, 'co2_mixing_ratio_c13', mem%co2_mixing_ratio_c13, &
            hgrid, vgrid_surface, &
            t_cf('co2_mixing_ratio_C13', 'ppm', 'surface air 13co2 mixing ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= co2_mixing_ratio_c13_init )

    CALL jsb_add_var( mem, 'co2_mixing_ratio_C14', mem%co2_mixing_ratio_c14, &
            hgrid, vgrid_surface, &
            t_cf('co2_mixing_ratio_c14', 'scaled', 'surface air 14co2 mixing ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= co2_mixing_ratio_c14_init )

    CALL jsb_add_var( mem, 'ga', mem%ga, &
            hgrid, vgrid_surface, &
            t_cf('ga', 'm s-1', 'aerodynamic conductance'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= ga_init )

    CALL jsb_add_var( mem, 'nhx_deposition', mem%nhx_deposition, &
            hgrid, vgrid_surface, &
            t_cf('nhx_deposition', 'mumol m-2 sec-1', 'surface downward reduced nitrogen deposition velocity'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'noy_deposition', mem%noy_deposition, &
            hgrid, vgrid_surface, &
            t_cf('noy_deposition', 'mumol m-2 sec-1', 'surface downward oxidised nitrogen deposition velocity'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'nhx_n15_deposition', mem%nhx_n15_deposition, &
            hgrid, vgrid_surface, &
            t_cf('nhx_n15_deposition', 'mumol m-2 sec-1', 'surface downward reduced nitrogen-15 deposition velocity'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'noy_n15_deposition', mem%noy_n15_deposition, &
            hgrid, vgrid_surface, &
            t_cf('noy_n15_deposition', 'mumol m-2 sec-1', 'surface downward oxidised nitrogen-15 deposition velocity'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'p_deposition', mem%p_deposition, &
            hgrid, vgrid_surface, &
            t_cf('p_deposition', 'mumol m-2 sec-1', 'surface downward reduced phosphoruns deposition velocity'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

!    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0 .AND. .NOT. ASSOCIATED(mem%t_acoef_wtr%ptr)) THEN

      CALL jsb_add_var( mem, 'albedo_lwtr', mem%albedo_lwtr, &
            hgrid, vgrid_surface, &
            t_cf('albedo_lwtr', 'unitless', 'lake albedo'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.07_wp )

      CALL jsb_add_var( mem, 'drag_wtr', mem%drag_wtr, &
            hgrid, vgrid_surface, &
            t_cf('drag_wtr', 'unitless', 'drag coefficient over water'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL  jsb_add_var( mem, 'drag_ice', mem%drag_ice, &
            hgrid, vgrid_surface, &
            t_cf('drag_ice', 'unitless', 'drag coefficient over ice'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 't_acoef_wtr', mem%t_acoef_wtr, &
            hgrid, vgrid_surface, &
            t_cf('t_acoef_wtr', 'unitless', 'Richtmeyer-Morten coefficient A for temperature over water'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 't_bcoef_wtr', mem%t_bcoef_wtr, &
            hgrid, vgrid_surface, &
            t_cf('t_bcoef_wtr', 'unitless', 'Richtmeyer-Morten coefficient A for temperature over water'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'q_acoef_wtr', mem%q_acoef_wtr, &
            hgrid, vgrid_surface, &
            t_cf('q_acoef_wtr', 'unitless', 'Richtmeyer-Morten coefficient A for specific humidity over water'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'q_bcoef_wtr', mem%q_bcoef_wtr, &
            hgrid, vgrid_surface, &
            t_cf('q_bcoef_wtr', 'unitless', 'Richtmeyer-Morten coefficient A for specific humidity over water'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 't_acoef_ice', mem%t_acoef_ice, &
            hgrid, vgrid_surface, &
            t_cf('t_acoef_ice', 'unitless', 'Richtmeyer-Morten coefficient A for temperature over ice'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 't_bcoef_ice', mem%t_bcoef_ice, &
            hgrid, vgrid_surface, &
            t_cf('t_bcoef_ice', 'unitless', 'Richtmeyer-Morten coefficient A for temperature over ice'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'q_acoef_ice', mem%q_acoef_ice, &
            hgrid, vgrid_surface, &
            t_cf('q_acoef_ice', 'unitless', 'Richtmeyer-Morten coefficient A for specific humidity over ice'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'q_bcoef_ice', mem%q_bcoef_ice, &
            hgrid, vgrid_surface, &
            t_cf('q_bcoef_ice', 'unitless', 'Richtmeyer-Morten coefficient A for specific humidity over ice'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

!    END IF   ! IF lakes

    CALL jsb_add_var( mem, 'prev_one_cos_zenith_angle', mem%prev_one_cos_zenith_angle, &
            hgrid, vgrid_surface, &
            t_cf('prev_one_cos_zenith_angle', 'unitless', 'cos_zenith_angle of the previous timestep'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'prev_two_cos_zenith_angle', mem%prev_two_cos_zenith_angle, &
            hgrid, vgrid_surface, &
            t_cf('prev_two_cos_zenith_angle', 'unitless', 'cos_zenith_angle of the before previous timestep'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )


  END SUBROUTINE Init_a2l_memory

#endif
END MODULE mo_a2l_memory_class
