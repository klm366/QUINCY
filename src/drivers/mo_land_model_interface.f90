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
!! @brief
!!  interface for land model \n
!!  run QUINCY for one time step
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!!  First version by Sönke Zaehle (2015-12-04) \n
!!  Reworked by Jan Engel (2016 July)
!!
MODULE mo_land_model_interface

#include "dsl_jsbach.inc"

  USE mo_kind
  USE mo_exception

  IMPLICIT NONE


  PRIVATE
  PUBLIC :: land_model_interface

  CHARACTER(len=*), PARAMETER :: modname = 'mo_land_model_interface'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> land model interface
  !!
  !! selection of the model via the namlist.slm - BASE_CTL:slm_model_name
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE land_model_interface(tile                            , &
                                  l_forcing                       , &
                                  l_output                        , &
                                  l_output_interval_match         , &
                                  fluxnet_diagnostics_sim_state   , &
                                  output_interval_flux            , &
                                  time_day, time_year             , &
                                  timestep_prop_day               , &
                                  !! a2l memory - INTENT(in) - start
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
                                  swpar_fdiffuse     , &
                                  press_srf          , &
                                  drag_srf           , &
                                  t_acoef            , &
                                  t_bcoef            , &
                                  q_acoef            , &
                                  q_bcoef            , &
                                  pch                , &
                                  cos_zenith_angle   , &
                                  co2_mixing_ratio     , &
                                  co2_mixing_ratio_C13 , &
                                  co2_mixing_ratio_C14 , &
                                  ga                   , &
                                  nhx_deposition     , &
                                  noy_deposition     , &
                                  nhx_n15_deposition , &
                                  noy_n15_deposition , &
                                  p_deposition       , &
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
                                  q_bcoef_ice        , &
                                  !! a2l memory - end
                                  !! start - variables for passing values from land to atmosphere
                                  !! INTENT(out)
                                  t_srf             ,   &
                                  t_eff_srf         ,   &
                                  qsat_srf          ,   &
                                  fact_q_air        ,   &
                                  fact_qsat_srf     ,   &
                                  zh_srf            ,   &
                                  zm_srf            ,   &
                                  alb_vis_dif       ,   &
                                  alb_nir_dif       ,   &
                                  alb_vis_dir       ,   &
                                  alb_nir_dir       ,   &
                                  albedo            ,   &
                                  sensible_heat_flx ,   &
                                  latent_heat_flx   ,   &
                                  net_radiation     ,   &
                                  vpd               ,   &
                                  ! For lakes:
                                  latent_hflx_wtr   ,   &
                                  sensible_hflx_wtr ,   &
                                  latent_hflx_ice   ,   &
                                  sensible_hflx_ice     &
                                  !! end - variables for passing values from land to atmosphere
                                  )


    USE mo_jsb_class,                     ONLY: Get_model
    USE mo_jsb_tile_class,                ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,                ONLY: t_jsb_process_task
    USE mo_jsb_model_class,               ONLY: t_jsb_model
    USE mo_jsb_process_class,             ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_, &
                                                Get_process_name
    USE mo_hsm_class,                     ONLY: t_Message
    USE mo_jsb4_forcing,                  ONLY: forcing_options
    USE mo_quincy_update_time_averages,   ONLY: update_time_averages  ! includes daytime and moving average
    USE mo_atmland_interface,             ONLY: update_a2l, update_l2a
    USE mo_quincy_output,                 ONLY: write_output_pool, write_output_flux, write_soil_forcing
    USE mo_test_canopy_model,             ONLY: run_test_canopy_model
    USE mo_test_radiation_model,          ONLY: run_test_radiation_model

    ! Process interfaces:
    USE mo_rad_interface
    USE mo_veg_interface
    USE mo_assimi_interface
    USE mo_sse_interface
    USE mo_sb_interface
    USE mo_pheno_interface

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile                           !< one tile with data structure for one lct
    LOGICAL,                    INTENT(in)        :: l_forcing, &                   !< IF .TRUE. THEN update_forcing(), update_a2l(), update_l2a()
                                                     l_output, &                    !< TRUE of output is enabled at all (i.e., fluxes are accumulated)
                                                     l_output_interval_match        !< TRUE if output may be created at the current timestep
    INTEGER,                    INTENT(in)        :: fluxnet_diagnostics_sim_state  !< possible values: 0, 1, 2, 3; see land_model_driver
    REAL(wp),                   INTENT(in)        :: output_interval_flux           !< flux output interval
    REAL(wp),                   INTENT(in)        :: time_day                       !< displays the 'real simulation time' in days; time_day == 1.0 after one day has passed
    INTEGER,                    INTENT(in)        :: time_year                      !< displays the 'real simulation time' in years; the 1st year of a simulation == 1 (not zero)
    REAL(wp),                   INTENT(in)        :: timestep_prop_day              !< proportion of one day that is simulated with one timestep (dtime/one_day)
    REAL(wp),                   INTENT(in)        :: &                 !! a2l memory - INTENT(in) - start
                                  t_air(:)              , &
                                  q_air(:)              , &
                                  rain(:)               , &
                                  snow(:)               , &
                                  wind_air(:)           , &
                                  wind_10m(:)           , &
                                  lw_srf_down(:)        , &
                                  swvis_srf_down(:)     , &
                                  swnir_srf_down(:)     , &
                                  swpar_srf_down(:)     , &
                                  swpar_fdiffuse(:)     , &
                                  press_srf(:)          , &
                                  drag_srf(:)           , &
                                  t_acoef(:)            , &
                                  t_bcoef(:)            , &
                                  q_acoef(:)            , &
                                  q_bcoef(:)            , &
                                  pch(:)                , &
                                  cos_zenith_angle(:)   , &
                                  co2_mixing_ratio(:)     , &
                                  co2_mixing_ratio_C13(:) , &
                                  co2_mixing_ratio_C14(:) , &
                                  ga(:)                   , &
                                  nhx_deposition(:)     , &
                                  noy_deposition(:)     , &
                                  nhx_n15_deposition(:) , &
                                  noy_n15_deposition(:) , &
                                  p_deposition(:)       , &
                                  ! For lakes:
                                  albedo_lwtr(:)        , &
                                  drag_wtr(:)           , &
                                  drag_ice(:)           , &
                                  t_acoef_wtr(:)        , &
                                  t_bcoef_wtr(:)        , &
                                  q_acoef_wtr(:)        , &
                                  q_bcoef_wtr(:)        , &
                                  t_acoef_ice(:)        , &
                                  t_bcoef_ice(:)        , &
                                  q_acoef_ice(:)        , &
                                  q_bcoef_ice(:)
                                                     !! a2l memory - end
    REAL(wp),INTENT(out)          :: &                      !! start - variables for passing values from land to atmosphere
                                  t_srf(:)             ,   &
                                  t_eff_srf(:)         ,   &
                                  qsat_srf(:)          ,   &
                                  fact_q_air(:)        ,   &
                                  fact_qsat_srf(:)     ,   &
                                  zh_srf(:)            ,   &
                                  zm_srf(:)            ,   &
                                  alb_vis_dif(:)       ,   &
                                  alb_nir_dif(:)       ,   &
                                  alb_vis_dir(:)       ,   &
                                  alb_nir_dir(:)       ,   &
                                  albedo(:)            ,   &
                                  sensible_heat_flx(:) ,   &
                                  latent_heat_flx(:)   ,   &
                                  net_radiation(:)     ,   &
                                  vpd(:)               ,   &
                                  ! For lakes:
                                  latent_hflx_wtr(:)   ,   &
                                  sensible_hflx_wtr(:) ,   &
                                  latent_hflx_ice(:)   ,   &
                                  sensible_hflx_ice(:)      !! end - variables for passing values from land to atmosphere
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),          POINTER   :: model
    CLASS(t_jsb_process_task),  POINTER   :: task
    CLASS(t_Message),           POINTER   :: msg
    CHARACTER(len=*),           PARAMETER :: routine = TRIM(modname)//':land_model_interface'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)    
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(SB_)



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 pass atmland variables from forcing input (driver) to a2l
    !!
    !!
    IF (l_forcing) THEN       ! l_forcing is true
      CALL update_a2l( &
      ! the sequence of input variables is consistent with the sequence of var in t_a2l_memory & atmland_interface
      tile                                         , &
      model%options                                , &
      t_air                   = t_air(:)              , &
      q_air                   = q_air(:)              , &
      rain                    = rain(:)               , &
      snow                    = snow(:)               , &
      wind_air                = wind_air(:)           , &
      wind_10m                = wind_10m(:)           , &
      lw_srf_down             = lw_srf_down(:)        , &
      swvis_srf_down          = swvis_srf_down(:)     , &
      swnir_srf_down          = swnir_srf_down(:)     , &
      swpar_srf_down          = swpar_srf_down(:)     , &
      swpar_fdiffuse          = swpar_fdiffuse(:)     , &
      press_srf               = press_srf(:)          , &
      drag_srf                = drag_srf(:)           , &
      t_acoef                 = t_acoef(:)            , &
      t_bcoef                 = t_bcoef(:)            , &
      q_acoef                 = q_acoef(:)            , &
      q_bcoef                 = q_bcoef(:)            , &
      pch                     = pch(:)                , &
      cos_zenith_angle        = cos_zenith_angle(:)   , &
      co2_mixing_ratio        = co2_mixing_ratio(:)     , &
      co2_mixing_ratio_C13    = co2_mixing_ratio_C13(:) , &
      co2_mixing_ratio_C14    = co2_mixing_ratio_C14(:) , &
      ga                      = ga(:)                   , &
      nhx_deposition          = nhx_deposition(:)     , &
      noy_deposition          = noy_deposition(:)     , &
      nhx_n15_deposition      = nhx_n15_deposition(:) , &
      noy_n15_deposition      = noy_n15_deposition(:) , &
      p_deposition            = p_deposition(:)       , &
      ! For lakes:
      albedo_lwtr             = albedo_lwtr(:)        , &
      drag_wtr                = drag_wtr(:)           , &
      drag_ice                = drag_ice(:)           , &
      t_acoef_wtr             = t_acoef_wtr(:)        , &
      t_bcoef_wtr             = t_bcoef_wtr(:)        , &
      q_acoef_wtr             = q_acoef_wtr(:)        , &
      q_bcoef_wtr             = q_bcoef_wtr(:)        , &
      t_acoef_ice             = t_acoef_ice(:)        , &
      t_bcoef_ice             = t_bcoef_ice(:)        , &
      q_acoef_ice             = q_acoef_ice(:)        , &
      q_bcoef_ice             = q_bcoef_ice(:)          &
      & )
    END IF


    !> 2.0 Run the tasks defined for the specific SLM model
    !! 
    !! land plant soil canopy \n 
    !! test_canopy is a special case and runs for only one timestep
    !! 
    SELECT CASE(TRIM(model%config%slm_model_name))
    !! land ---------------------------------------------------------------------------------------
    CASE('land')

      !! init/reset fluxes (with/to zero)
      !!
      CALL reset_assimi_fluxes    (tile, model%options)   ! mo_assimi_util
      CALL reset_veg_fluxes       (tile, model%options)   ! mo_veg_util
      CALL reset_sb_fluxes        (tile, model%options)   ! mo_sb_util
      CALL reset_sse_fluxes       (tile, model%options)   ! mo_sse_util

      !! run Tasks
      !!
      CALL update_sb_rate_modifier(tile, model%options)
      CALL update_soil_turbulence (tile, model%options)
      CALL update_radiation       (tile, model%options)
      CALL update_canopy_fluxes   (tile, model%options)
      CALL update_sse_physics     (tile, model%options)
      CALL update_phenology       (tile, model%options)
      CALL update_veg_turnover    (tile, model%options)
      CALL update_veg_dynamics    (tile, model%options)
      SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
        CASE("simple_0d", "simple_1d")
          CALL update_sb_simple_model (tile, model%options)
        CASE('jsm')
          CALL update_sb_jsm          (tile, model%options)
      ENDSELECT
      CALL update_veg_growth      (tile, model%options)
      CALL update_veg_pools       (tile, model%options)
      CALL update_sb_pools        (tile, model%options)

      !! variable accumulation over time - daytime & moving average
      !!
      CALL update_time_averages(tile, model%options, time_day, timestep_prop_day)    ! see shared

      
      !! create soil forcing for quincy-jsm standalone
      !! the output-routine is called each timestep, hence, output interval is timestep
      !! output litterfall and soilphysics from a static simulation
      IF (model%config%flag_generate_soil_forcing_stat) THEN
        IF(.NOT. model%config%l_transient_spinup .AND. &
           & time_year >= model%config%output_start_first_day_year .AND. time_year <= model%config%output_end_last_day_year) THEN
          CALL write_soil_forcing(tile, model%options, time_day)
        ENDIF
      ENDIF
      !! output litterfall and soilphysics from a transient fluxnet/3 diagnostics simulation
      IF (model%config%flag_generate_soil_forcing_trans .AND. l_output) THEN
        ! write soil forcing only for the last full forcing period of the particular fluxnet site (where static forcing is available)
        IF (fluxnet_diagnostics_sim_state == 2) THEN
          CALL write_soil_forcing(tile, model%options, time_day)
        ENDIF
      ENDIF

    CASE('plant')
    !! plant --------------------------------------------------------------------------------------

      !! init/reset fluxes (with/to zero)
      !!
      CALL reset_assimi_fluxes    (tile, model%options)  ! mo_assimi_util
      CALL reset_veg_fluxes       (tile, model%options)  ! mo_veg_util
      CALL reset_sse_fluxes       (tile, model%options)  ! mo_sse_util

      !! run Tasks
      !!
      CALL update_soil_turbulence (tile, model%options)
      CALL update_radiation       (tile, model%options)
      CALL update_canopy_fluxes   (tile, model%options)
      CALL update_sse_physics     (tile, model%options)
      CALL update_phenology       (tile, model%options)
      CALL update_veg_turnover    (tile, model%options)
      CALL update_veg_dynamics    (tile, model%options)
      CALL update_veg_growth      (tile, model%options)
      CALL update_veg_pools       (tile, model%options)

      !! variable accumulation over time - daytime & moving average
      !!
      CALL update_time_averages(tile, model%options, time_day, timestep_prop_day)    ! see shared

    CASE('soil')
    !! soil ---------------------------------------------------------------------------------------

      !! init/reset fluxes (with/to zero)
      !!
      CALL reset_sb_fluxes            (tile, model%options)   ! mo_sb_util
      CALL reset_sse_fluxes           (tile, model%options)   ! mo_sse_util

      CALL update_sb_rate_modifier    (tile, model%options)
      SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
        CASE("simple_0d", "simple_1d")
          CALL update_sb_simple_model (tile, model%options)
        CASE('jsm')
          CALL update_sb_jsm          (tile, model%options)
      ENDSELECT
      CALL update_sb_pools            (tile, model%options)

      !! accumulation over time - daytime & moving average
      !!
      CALL update_time_averages(tile, model%options, time_day, timestep_prop_day)    ! see shared

    CASE('canopy')
    !! canopy -------------------------------------------------------------------------------------

      !! init/reset fluxes (with/to zero)
      !!
      CALL reset_assimi_fluxes    (tile, model%options)  ! mo_assimi_util
      CALL reset_veg_fluxes       (tile, model%options)  ! mo_veg_util
      CALL reset_sse_fluxes       (tile, model%options)  ! mo_sse_util

      !! run Tasks
      !!
      CALL update_soil_turbulence (tile, model%options)
      CALL update_radiation       (tile, model%options)
      CALL update_canopy_fluxes   (tile, model%options)
      CALL update_sse_physics     (tile, model%options)
      CALL update_phenology       (tile, model%options)
      CALL update_veg_pools       (tile, model%options)

      !! accumulation over time - daytime & moving average
      !!
      CALL update_time_averages(tile, model%options, time_day, timestep_prop_day)    ! see shared
      
    CASE('test_canopy')
    !! test_canopy --------------------------------------------------------------------------------
      CALL run_test_canopy_model  (tile, model%options)

    CASE('test_radiation')
    !! test_radiation -----------------------------------------------------------------------------
      CALL run_test_radiation_model (tile, model%options)

    CASE DEFAULT
      CALL finish('land_model_interface','Invalid SLM model name. Try either land, plant, soil, canopy or test_canopy.')
    END SELECT  ! select SLM model !


    !> 2.1 Run output routines for land, plant, soil, canopy
    !! 
    !!
    IF (l_output) THEN 
      SELECT CASE(TRIM(model%config%slm_model_name))
      CASE('land','plant','soil','canopy')

        !! output pool: write output in certain intervals, see mo_quincy_output
        !!  note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
        !!  note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
        !! 
        IF (l_output_interval_match) THEN
          CALL write_output_pool(tile, model%options, time_day)
        ENDIF
        !! output flux: write output in certain intervals, see mo_quincy_output
        !!  if output is enabled, flux variables are accumulated
        !!  output is written when the interval is a match with current timestep
        !! 
        CALL write_output_flux(tile, model%options, time_day, output_interval_flux, &
                              l_flux_output_interval=l_output_interval_match, l_zero_acc_variables=.FALSE.)

      END SELECT  ! select SLM model !
    ENDIF                            

    !> 3.0 variables for passing values from land to atmosphere (via the driver)
    !! 
    !! 
    IF (l_forcing) THEN       ! l_forcing is true if update_a2l() is called (depends on the SLM model selected)
      CALL update_l2a( &
      tile                                              ,   &
      model%options                                     ,   &
      t_srf                   = t_srf(:)             ,   &
      t_eff_srf               = t_eff_srf(:)         ,   &
      qsat_srf                = qsat_srf(:)          ,   &
      fact_q_air              = fact_q_air(:)        ,   &
      fact_qsat_srf           = fact_qsat_srf(:)     ,   &
      zh_srf                  = zh_srf(:)            ,   &
      zm_srf                  = zm_srf(:)            ,   &
      alb_vis_dif             = alb_vis_dif(:)       ,   &
      alb_nir_dif             = alb_nir_dif(:)       ,   &
      alb_vis_dir             = alb_vis_dir(:)       ,   &
      alb_nir_dir             = alb_nir_dir(:)       ,   &
      albedo                  = albedo(:)            ,   &
      sensible_heat_flx       = sensible_heat_flx(:) ,   &
      latent_heat_flx         = latent_heat_flx(:)   ,   &
      net_radiation           = net_radiation(:)     ,   &
      vpd                     = vpd(:)               ,   &
      ! For lakes:
      latent_hflx_wtr         = latent_hflx_wtr(:)   ,   &
      sensible_hflx_wtr       = sensible_hflx_wtr(:) ,   &
      latent_hflx_ice         = latent_hflx_ice(:)   ,   &
      sensible_hflx_ice       = sensible_hflx_ice(:)     &
      & )
    END IF

  END SUBROUTINE land_model_interface

END MODULE mo_land_model_interface
