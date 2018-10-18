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
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d              ! *slmjsb4
  USE mo_slm_varlist,            ONLY: jsb_add_var   ! replace by mo_jsb_varlist later  *slmtemp

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
                              t_air                 , &  !> 2m air temperature (K) 
                              q_air                 , &  !> 2m air specific humidity (g g-1) 
                              rain                  , &  !> liquid water precipitation (kg m-2 s-1)
                              snow                  , &  !> solid water precipitation (kg m-2 s-1)
                              wind_air              , &  !> 2m air wind speed (m s-1) 
                              wind_10m              , &  !> 10m air wind speed (m s-1) 
                              lw_srf_down           , &  !> longwave surface downward flux (W m-2) 
                              swvis_srf_down        , &  !> shortwave surface downward flux in the visible band (W m-2)
                              swnir_srf_down        , &  !> shortwave surface downward flux in the near infrared band (W m-2) 
                              swpar_srf_down        , &  !> shortwave surface downward flux in the PAR band (W m-2) 
                              swpar_fdiffuse        , &  !> fraction of shortwave surface downward flux that is diffuse 
                              press_srf             , &  !> surface air pressure (Pa)
                              drag_srf              , &  !> total surface drag () 
                              t_acoef               , &  !> Richtmeyer-Morten coefficient A for temperature () 
                              t_bcoef               , &  !> Richtmeyer-Morten coefficient B for temperature ()
                              q_acoef               , &  !> Richtmeyer-Morten coefficient A for specific humidity () 
                              q_bcoef               , &  !> Richtmeyer-Morten coefficient B for specific humidity ()
                              pch                   , &  !> total surface drag () 
                              cos_zenith_angle           !> Cosine of solar zenith angle 

    ! CO2
    TYPE(t_jsb_var_real2d) :: &
                              ! jsbach4 variables, replaced by co2_mixing_ratio  _c13  _c14
                              ! CO2_air               , &  !> 
                              ! CO2_air_mol           , &  !> 
                              
                              ! from SLM, replacing CO2_air and CO2_air_mol from jsbach4
                              co2_mixing_ratio      , &  !> surface air co2 mixing ratio (ppm) 
                              co2_mixing_ratio_C13  , &  !> surface air 13co2 mixing ratio (ppm)
                              co2_mixing_ratio_C14       !> surface air 14co2 mixing ratio (scaled)
                              
    ! ga    @todo  ga in t_a2l_memory is a variable from SLM, which may be moved elsewhere
    TYPE(t_jsb_var_real2d) :: &
                              ga                         !> aerodynamic conductance (m s-1)

    ! Deposition
    TYPE(t_jsb_var_real2d) :: &                                                            
                              nhx_deposition        , &  !> surface downward reduced nitrogen deposition velocity [mumol m-2 sec-1]
                              noy_deposition        , &  !> surface downward oxidised nitrogen deposition velocity [mumol m-2 sec-1]
                              nhx_n15_deposition    , &  !> surface downward reduced nitrogen-15 deposition velocity [mumol m-2 sec-1]
                              noy_n15_deposition    , &  !> surface downward oxidised nitrogen-15 deposition velocity [mumol m-2 sec-1]
                              p_deposition               !> surface downward reduced phosphoruns deposition velocity [mumol m-2 sec-1]

    ! For lakes
    TYPE(t_jsb_var_real2d) :: &                              
                              albedo_lwtr           , &  !> lake albedo ()
                              drag_wtr              , &  !> drag coefficient over water ()
                              drag_ice              , &  !> drag coefficient over ice ()
                              t_acoef_wtr           , &  !> Richtmeyer-Morten coefficient A for temperature over water ()
                              t_bcoef_wtr           , &  !> Richtmeyer-Morten coefficient A for temperature over water ()
                              q_acoef_wtr           , &  !> Richtmeyer-Morten coefficient A for specific humidity over water ()
                              q_bcoef_wtr           , &  !> Richtmeyer-Morten coefficient A for specific humidity over water ()
                              t_acoef_ice           , &  !> Richtmeyer-Morten coefficient A for temperature over ice ()
                              t_bcoef_ice           , &  !> Richtmeyer-Morten coefficient A for temperature over ice ()
                              q_acoef_ice           , &  !> Richtmeyer-Morten coefficient A for specific humidity over ice ()
                              q_bcoef_ice                !> Richtmeyer-Morten coefficient A for specific humidity over ice ()

    ! For calculation of the timestep after local midnight (1st timestep of the new day)
    TYPE(t_jsb_var_real2d) :: &
                              prev_one_cos_zenith_angle   , &   ! cos_zenith_angle of the previous timestep
                              prev_two_cos_zenith_angle         ! cos_zenith_angle of the before previous timestep

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

!    USE mo_jsb_varlist,       ONLY: jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2   ! , tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model


    ! needed to calc SLM init values:
    USE mo_atmland_constants, ONLY: def_wind, def_co2_deltaC13, def_co2_deltaC14, def_t_air, standard_press_srf, def_vpd, &
                                    eps_vpd, def_co2_mixing_ratio, def_swdown, frac_vis_swtot_srf, frac_par_swvis_srf, &
                                    def_fdiffuse, def_cos_angle
    USE mo_shared_isotope,    ONLY: calc_mixing_ratio_C13C12, calc_mixing_ratio_C14C
    USE mo_atmland_util,      ONLY: calc_spec_humidity_sat

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_a2l_memory),   INTENT(inout), TARGET :: mem
    CHARACTER(len=*),      INTENT(in)            :: prefix
    CHARACTER(len=*),      INTENT(in)            :: suffix
    INTEGER,               INTENT(in)            :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,               INTENT(in)            :: model_id
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
    TYPE(t_jsb_model),      POINTER   :: model 
    TYPE(t_jsb_grid),       POINTER   :: hgrid 
    TYPE(t_jsb_vgrid),      POINTER   :: surface  
    INTEGER                           :: table
    CHARACTER(len=*),       PARAMETER :: routine = TRIM(modname)//':Init_a2l_memory'

    model => Get_model(model_id) 

!    table = tables(1)  
    table = 1

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
    CALL jsb_add_var( mem, 't_air', mem%t_air, &
            hgrid, surface, &
            t_cf('t_air', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= t_air_init )

    CALL jsb_add_var( mem, 'q_air', mem%q_air, &
            hgrid, surface, &
            t_cf('q_air', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= q_air_init )

    CALL jsb_add_var( mem, 'rain', mem%rain, &
            hgrid, surface, &
            t_cf('rain', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'snow', mem%snow, &
            hgrid, surface, &
            t_cf('snow', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'wind_air', mem%wind_air, &
            hgrid, surface, &
            t_cf('wind_air', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= wind_air_init )

    CALL jsb_add_var( mem, 'wind_10m', mem%wind_10m, &
            hgrid, surface, &
            t_cf('wind_10m', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= wind_10m_init )

    CALL jsb_add_var( mem, 'lw_srf_down', mem%lw_srf_down, &
            hgrid, surface, &
            t_cf('lw_srf_down', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'swvis_srf_down', mem%swvis_srf_down, &
            hgrid, surface, &
            t_cf('swvis_srf_down', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= swvis_srf_down_init )

    CALL jsb_add_var( mem, 'swnir_srf_down', mem%swnir_srf_down, &
            hgrid, surface, &
            t_cf('swnir_srf_down', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= swnir_srf_down_init )

    CALL jsb_add_var( mem, 'swpar_srf_down', mem%swpar_srf_down, &
            hgrid, surface, &
            t_cf('swpar_srf_down', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= swpar_srf_down_init )

    CALL jsb_add_var( mem, 'swpar_fdiffuse', mem%swpar_fdiffuse, &
            hgrid, surface, &
            t_cf('swpar_fdiffuse', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= swpar_fdiffuse_init )

    CALL jsb_add_var( mem, 'press_srf', mem%press_srf, &
            hgrid, surface, &
            t_cf('press_srf', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= press_srf_init )

    CALL jsb_add_var( mem, 'drag_srf', mem%drag_srf, &
            hgrid, surface, &
            t_cf('drag_srf', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= drag_srf_init )

    CALL jsb_add_var( mem, 't_acoef', mem%t_acoef, &
            hgrid, surface, &
            t_cf('t_acoef', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 't_bcoef', mem%t_bcoef, &
            hgrid, surface, &
            t_cf('t_bcoef', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'q_acoef', mem%q_acoef, &
            hgrid, surface, &
            t_cf('q_acoef', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'q_bcoef', mem%q_bcoef, &
            hgrid, surface, &
            t_cf('q_bcoef', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'pch', mem%pch, &
            hgrid, surface, &
            t_cf('pch', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'cos_zenith_angle', mem%cos_zenith_angle, &
            hgrid, surface, &
            t_cf('cos_zenith_angle', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= cos_zenith_angle_init )

!!    CALL jsb_add_var( mem, 'CO2_air', mem%CO2_air, &
!!!      & hgrid, surface,                                                                                &
!!!      & t_cf('CO2_air', '[kg(CO2)/kg(air)]', 'CO2 mass mixing ratio of lowest atmosphere level'),      &
!!!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
!!!      & prefix, suffix,                                                                                &
!!!      & lrestart=.FALSE., &
!!      & initval_r=0.0_wp )

!!    CALL jsb_add_var( mem, 'CO2_air_mol', mem%CO2_air_mol, &
!!!      & hgrid, surface,                                                                                &
!!!      & t_cf('CO2_air_mol', '[mol(CO2)/mol(air)]', 'CO2 mol mixing ratio of lowest atmosphere level'), &
!!!      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
!!!      & prefix, suffix,                                                                                &
!!!      & lrestart=.FALSE., &
!!      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'co2_mixing_ratio', mem%co2_mixing_ratio, &
            hgrid, surface, &
            t_cf('co2_mixing_ratio', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= co2_mixing_ratio_init )

    CALL jsb_add_var( mem, 'co2_mixing_ratio_c13', mem%co2_mixing_ratio_c13, &
            hgrid, surface, &
            t_cf('co2_mixing_ratio_c13', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= co2_mixing_ratio_c13_init )

    CALL jsb_add_var( mem, 'co2_mixing_ratio_c14', mem%co2_mixing_ratio_c14, &
            hgrid, surface, &
            t_cf('co2_mixing_ratio_c14', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= co2_mixing_ratio_c14_init )

    CALL jsb_add_var( mem, 'ga', mem%ga, &
            hgrid, surface, &
            t_cf('ga', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r= ga_init )

    CALL jsb_add_var( mem, 'nhx_deposition', mem%nhx_deposition, &
            hgrid, surface, &
            t_cf('nhx_deposition', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'noy_deposition', mem%noy_deposition, &
            hgrid, surface, &
            t_cf('noy_deposition', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'nhx_n15_deposition', mem%nhx_n15_deposition, &
            hgrid, surface, &
            t_cf('nhx_n15_deposition', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'noy_n15_deposition', mem%noy_n15_deposition, &
            hgrid, surface, &
            t_cf('noy_n15_deposition', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'p_deposition', mem%p_deposition, &
            hgrid, surface, &
            t_cf('p_deposition', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

!    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0 .AND. .NOT. ASSOCIATED(mem%t_acoef_wtr%ptr)) THEN

      CALL jsb_add_var( mem, 'drag_wtr', mem%drag_wtr, &
            hgrid, surface, &
            t_cf('drag_wtr', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL  jsb_add_var( mem, 'drag_ice', mem%drag_ice, &
            hgrid, surface, &
            t_cf('drag_ice', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 't_acoef_wtr', mem%t_acoef_wtr, &
            hgrid, surface, &
            t_cf('t_acoef_wtr', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 't_bcoef_wtr', mem%t_bcoef_wtr, &
            hgrid, surface, &
            t_cf('t_bcoef_wtr', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'q_acoef_wtr', mem%q_acoef_wtr, &
            hgrid, surface, &
            t_cf('q_acoef_wtr', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'q_bcoef_wtr', mem%q_bcoef_wtr, &
            hgrid, surface, &
            t_cf('q_bcoef_wtr', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 't_acoef_ice', mem%t_acoef_ice, &
            hgrid, surface, &
            t_cf('t_acoef_ice', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 't_bcoef_ice', mem%t_bcoef_ice, &
            hgrid, surface, &
            t_cf('t_bcoef_ice', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'q_acoef_ice', mem%q_acoef_ice, &
            hgrid, surface, &
            t_cf('q_acoef_ice', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'q_bcoef_ice', mem%q_bcoef_ice, &
            hgrid, surface, &
            t_cf('q_bcoef_ice', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'albedo_lwtr', mem%albedo_lwtr, &
            hgrid, surface, &
            t_cf('albedo_lwtr', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.07_wp )

!    END IF   ! IF lakes

    CALL jsb_add_var( mem, 'prev_one_cos_zenith_angle', mem%prev_one_cos_zenith_angle, &
            hgrid, surface, &
            t_cf('prev_one_cos_zenith_angle', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'prev_two_cos_zenith_angle', mem%prev_two_cos_zenith_angle, &
            hgrid, surface, &
            t_cf('prev_two_cos_zenith_angle', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )


  END SUBROUTINE Init_a2l_memory

#endif
END MODULE mo_a2l_memory_class
