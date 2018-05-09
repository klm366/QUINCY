!>
!! @brief
!!  interface for land model \n
!!  run SLM for one time step
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
  SUBROUTINE land_model_interface(tile, &                ! one tile
                                  l_forcing, &           ! forcing on/off
                                  time_day, time_year, &
                                  timestep_prop_day  , &
                                  l_transient_spinup , &
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
                                 !! l2a memory - INTENT(out) - start
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
                                  !! l2a memory - end
                                  )


    USE mo_jsb_class,                     ONLY: Get_model
    USE mo_base_memory,                   ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,               ONLY: t_jsb_model
    USE mo_jsb_process_class,             ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_base_config_class,             ONLY: base_config
    USE mo_jsb4_forcing,                  ONLY: forcing_options
    USE mo_shared_update_time_averages,   ONLY: update_time_averages  ! includes daytime and moving average
    USE mo_atmland_interface,             ONLY: update_a2l, update_l2a
    USE mo_shared_output,                 ONLY: write_output_pool, write_output_flux, write_soil_forcing
    USE mo_veg_util,                      ONLY: reset_veg_fluxes
    USE mo_sb_util,                       ONLY: reset_sb_fluxes
    USE mo_soil_util,                     ONLY: reset_soil_fluxes
    USE mo_assimi_util,                   ONLY: reset_assimi_fluxes
    USE mo_test_canopy_model,             ONLY: run_test_canopy_model

    ! Process interfaces:
    USE mo_rad_interface
    USE mo_veg_interface
    USE mo_assimi_interface
    USE mo_soil_interface
    USE mo_sb_interface
    USE mo_pheno_interface

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)      :: tile              !< data structure for dynamic vegetation and soil variables
    LOGICAL,INTENT(in)                  :: l_forcing         !< IF .TRUE. THEN update_forcing(), update_a2l(), update_l2a()
    REAL(wp),INTENT(in)                 :: time_day          !< displays the 'real simulation time' in days; time_day == 1.0 after one day has passed
    INTEGER,INTENT(in)                  :: time_year         !< displays the 'real simulation time' in years; the 1st year of a simulation == 1 (not zero)
    REAL(wp),INTENT(in)                 :: timestep_prop_day !< proportion of one day that is simulated with one timestep (dtime/one_day)
    LOGICAL,INTENT(in)                  :: l_transient_spinup!< spinup phase of the transient forcing_mode
    REAL(wp),INTENT(in)           :: &               !! a2l memory - INTENT(in) - start
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
    REAL(wp),INTENT(out)          :: &              !! l2a memory - INTENT(out) - start
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
                                  sensible_hflx_ice(:) !! l2a memory - end
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER     :: model
    CHARACTER(len=*),       PARAMETER   :: routine = TRIM(modname)//':land_model_interface'
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


    !> 2.0 check stand harvest year (if this functionality is enabled: base_config\%flag_stand_harvest)
    !! 
    !! but enable harvesting only after transient spinup, i.e., when l_transient_spinup=FALSE
    IF (base_config%flag_stand_harvest .AND. .NOT. l_transient_spinup) THEN
      IF (base_config%l_do_stand_replacing_harvest) THEN
        base_config%l_do_stand_replacing_harvest = .FALSE.    ! per default this flag is FALSE
      END IF
      IF (forcing_options%forcing_year .EQ. base_config%stand_replacing_year .AND. &
          forcing_options%forcing_doy  .EQ. 1                                .AND. &
          forcing_options%forcing_hour .LT. eps4                                    ) THEN
        base_config%l_do_stand_replacing_harvest = .TRUE.     ! TRUE only for this specific timestep, per default this flag is FALSE
      END IF
    END IF


    !> 3.0 Run the tasks defined for the specific SLM model
    !! 
    !! land plant soil canopy \n 
    !! test_canopy is a species case and runs for only one timestep
    !! 
    SELECT CASE(TRIM(base_config%slm_model_name))
    !! land ---------------------------------------------------------------------------------------
    CASE('land')

      !! init/reset fluxes (with/to zero)
      !!
      CALL reset_assimi_fluxes(tile,model%options)! mo_assimi_util
      CALL reset_veg_fluxes(tile, model%options)  ! mo_veg_util
      CALL reset_sb_fluxes(tile, model%options)   ! mo_sb_util
      CALL reset_soil_fluxes(tile, model%options) ! mo_soil_util

      !! run Tasks
      !!
      CALL update_sb_rate_modifier(tile, model%options)
      CALL update_soil_turbulence (tile, model%options)
      CALL update_radiation       (tile, model%options)
      CALL update_canopy_fluxes   (tile, model%options)
      CALL update_soil_physics    (tile, model%options)
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

      !! accumulation over time - daytime & moving average
      !!
      CALL update_time_averages(tile, model%options, time_day, timestep_prop_day)    ! see shared

      !! output pool: write output in certain intervals, see mo_shared_output
      !! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
      !! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
      !!
      CALL write_output_pool(tile, model%options, time_day, time_year, l_transient_spinup)

      !! output flux: write output in certain intervals, see mo_shared_output
      !!
      !!
      CALL write_output_flux(tile, model%options, time_day, time_year, l_transient_spinup)

      !! output litterfall and soilphysics: usually this would be output at each timstep
      !!
      !!
      IF (base_config%flag_generate_soil_forcing) THEN
        CALL write_soil_forcing(tile, model%options, time_day, time_year, l_transient_spinup)
      END IF

    CASE('plant')
    !! plant --------------------------------------------------------------------------------------

      !! init/reset fluxes (with/to zero)
      !!
      CALL reset_assimi_fluxes(tile,model%options)! mo_assimi_util
      CALL reset_veg_fluxes(tile, model%options)  ! mo_veg_util
      CALL reset_soil_fluxes(tile, model%options) ! mo_soil_util

      !! run Tasks
      !!
      CALL update_soil_turbulence (tile, model%options)
      CALL update_radiation       (tile, model%options)
      CALL update_canopy_fluxes   (tile, model%options)
      CALL update_soil_physics    (tile, model%options)
      CALL update_phenology       (tile, model%options)
      CALL update_veg_turnover    (tile, model%options)
      CALL update_veg_dynamics    (tile, model%options)
      CALL update_veg_growth      (tile, model%options)
      CALL update_veg_pools       (tile, model%options)

      !! accumulation over time - daytime & moving average
      !!
      CALL update_time_averages(tile, model%options, time_day, timestep_prop_day)    ! see shared

      !! output pool: write output in certain intervals, see mo_shared_output
      !! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
      !! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
      !!
      CALL write_output_pool(tile, model%options, time_day, time_year, l_transient_spinup)

      !! output flux: write output in certain intervals, see mo_shared_output
      !!
      !!
      CALL write_output_flux(tile, model%options, time_day, time_year, l_transient_spinup)

      !! output litterfall and soilphysics: usually this would be output at each timstep
      !!
      !!
      IF (base_config%flag_generate_soil_forcing) THEN
        CALL write_soil_forcing(tile, model%options, time_day, time_year, l_transient_spinup)
      END IF

    CASE('soil')
    !! soil ---------------------------------------------------------------------------------------

      !! init/reset fluxes (with/to zero)
      !!
      CALL reset_sb_fluxes(tile, model%options)   ! mo_sb_util
      CALL reset_soil_fluxes(tile, model%options) ! mo_soil_util

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

      !! output pool: write output in certain intervals, see mo_shared_output
      !! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
      !! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
      !!
      CALL write_output_pool(tile, model%options, time_day, time_year, l_transient_spinup)

      !! output flux: write output in certain intervals, see mo_shared_output
      !!
      !!
      CALL write_output_flux(tile, model%options, time_day, time_year, l_transient_spinup)

    CASE('canopy')
    !! canopy -------------------------------------------------------------------------------------

      !! init/reset fluxes (with/to zero)
      !!
       CALL reset_veg_fluxes(tile, model%options)  ! mo_veg_util
      CALL reset_soil_fluxes(tile, model%options) ! mo_soil_util



      !! run Tasks
      !!
      CALL update_soil_turbulence (tile, model%options)
      CALL update_radiation       (tile, model%options)
      CALL update_canopy_fluxes   (tile, model%options)
      CALL update_soil_physics    (tile, model%options)
      CALL update_phenology       (tile, model%options)
      CALL update_veg_pools       (tile, model%options)

      !! accumulation over time - daytime & moving average
      !!
      CALL update_time_averages(tile, model%options, time_day, timestep_prop_day)    ! see shared
      
      !! output pool: write output in certain intervals, see mo_shared_output
      !! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
      !! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
      !!
      CALL write_output_pool(tile, model%options, time_day, time_year, l_transient_spinup)
      
      !! output flux: write output in certain intervals, see mo_shared_output
      !! ! only some flux files are written, see CASE statements in shared_output
      !!
      CALL write_output_flux(tile, model%options, time_day, time_year, l_transient_spinup)

    CASE('test_canopy')
    !! test_canopy --------------------------------------------------------------------------------
      CALL run_test_canopy_model  (tile, model%options)

    CASE DEFAULT
      CALL finish('land_model_interface','Invalid SLM model name. Try either land, plant, soil, canopy or test_canopy.')
    END SELECT  ! select SLM model !


    !> 4.0 pass atmland variables from l2a memory back to the driver
    !! 
    !! 
    IF (l_forcing) THEN       ! l_forcing is true if update_a2l() is called (depends on the SLM model selected)
      CALL update_l2a( &
      ! the sequence of input variables is consistent with the sequence of var in t_l2a_memory !
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
