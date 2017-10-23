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


    !> pass atmland variables from forcing input (driver) to a2l
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


    !> check stand harvest year (if this functionality is enabled: base_config\%flag_stand_harvest)
    !! 
    IF (base_config%flag_stand_harvest) THEN
      IF (base_config%l_do_stand_replacing_harvest) THEN
        base_config%l_do_stand_replacing_harvest = .FALSE.    ! per default this flag is FALSE
      END IF
      IF (forcing_options%forcing_year .EQ. base_config%stand_replacing_year .AND. &
          forcing_options%forcing_doy  .EQ. 1                                .AND. &
          forcing_options%forcing_hour .LT. eps4                                    ) THEN
        base_config%l_do_stand_replacing_harvest = .TRUE.     ! TRUE only for this specific timestep, per default this flag is FALSE
      END IF
    END IF


    !> @todo mo_land_model_interface:land_model_interface put all Tasks into a Queue accordingly to the model selection
    !! (land plant soil canopy test_canopy)
    SELECT CASE(TRIM(base_config%slm_model_name))
    !! land ---------------------------------------------------------------------------------------
    CASE('land')

      !! init/reset fluxes (with/to zero)
      !!
      CALL reset_veg_fluxes(tile, model%options)  ! mo_veg_util
      CALL reset_sb_fluxes(tile, model%options)   ! mo_sb_util

      !! run Tasks
      !!
      CALL update_sb_rate_modifier(tile, model%options)
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
      CALL reset_veg_fluxes(tile, model%options)  ! mo_veg_util

      !! run Tasks
      !!
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

      CALL update_sb_rate_modifier    (tile, model%options)
      SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
        CASE("simple_0d", "simple_1d")
          CALL update_sb_simple_model (tile, model%options)
        CASE('jsm')
          CALL update_sb_jsm          (tile, model%options)
      ENDSELECT
      CALL update_sb_pools            (tile, model%options)


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

      !! run Tasks
      !!
      CALL update_radiation       (tile, model%options)
      CALL update_canopy_fluxes   (tile, model%options)
      CALL update_soil_physics    (tile, model%options)
      CALL update_veg_pools       (tile, model%options)

      !! accumulation over time - daytime & moving average
      !!
      CALL update_time_averages(tile, model%options, time_day, timestep_prop_day)    ! see shared

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

    !> pass atmland variables from l2a memory back to the driver
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


  !-----------------------------------------------------------------------------------------------------
  !> Run the test_canopy model from the land model interface
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE run_test_canopy_model(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_rad_constants,         ONLY: rad2ppfd
    USE mo_atmland_constants,     ONLY: frac_vis_swtot_srf, frac_par_swvis_srf,eps_vpd, standard_press_srf, def_co2_mixing_ratio
    USE mo_atmland_util,          ONLY: calc_spec_humidity_sat
    USE mo_veg_util,              ONLY: reset_veg_fluxes
    USE mo_shared_constants,      ONLY: molar_mass_N, Tzero

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(ASSIMI_)    
    
    ! Use of process memories
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory

    ! Process interfaces:
    USE mo_atmland_interface                ! update_atmland
    USE mo_rad_interface                    ! update_radiation
    USE mo_assimi_interface                 ! update_canopy_fluxes
    USE mo_veg_canopy
    ! get parameters defined in init_slm_test_canopy()

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),           INTENT(inout)     :: tile         !< one tile
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    INTEGER, PARAMETER                        :: npft  = 8         !< number of PFTs simulated by the test_canopy model
    TYPE(t_jsb_model),      POINTER           :: model
    TYPE(t_lctlib_element), POINTER           :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp), DIMENSION(options%nc)           :: &
                  t_air           , &
                  press_srf       , &
                  q_air           , &
                  co2_mixing_ratio, &
                  rain            , &
                  snow           , &
                  swvis_srf_down , &
                  swnir_srf_down , &
                  swpar_srf_down , &
                  swpar_fdiffuse , &
                  wind_air       , &
                  wind_10m       , &
                  lw_srf_down    , &
                  drag_srf       , &
                  t_acoef        , &
                  t_bcoef        , &
                  q_acoef        , &
                  q_bcoef        , &
                  pch            , &
                  cos_zenith_angle, &
                  ga             , &
                  !latitude      , &
                  !longitude     , &
                  vpd
    REAL(wp),DIMENSION(options%nc)      :: hlp1
    INTEGER       :: i,lun
    INTEGER       :: ref_pft                          ! PFT identifier
    ! parameters (specific to the test_canopy model)
    REAL(wp),SAVE :: ref_t_air = 25.0_wp              ! temperature in deg C
    REAL(wp),SAVE :: ref_vpd = 1000._wp               ! vpd in Pa
    REAL(wp),SAVE :: ref_swdown = 652.173913_wp       ! PAR in W m-2
    REAL(wp),SAVE :: ref_fdiffuse = 0.2_wp            ! diffuse fraction
    REAL(wp),SAVE :: ref_cos_angle = 0.8_wp           ! cosine of zenith angle
    REAL(wp),SAVE :: ref_co2_mixing_ratio = 380.0_wp  ! co2 in ppm
    REAL(wp),SAVE :: ref_wind = 1.0_wp                ! wind in m/s
    REAL(wp),SAVE :: ref_theta = 1._wp                ! soil water availability in fraction of AWC
    REAL(wp),SAVE :: ref_lai = 3._wp                  ! LAI
    REAL(wp),SAVE :: ref_nleaf = 2.0_wp               ! N leaf per unit leaf area in gN/m2
    INTEGER       :: ipoi, ichunk                     ! number of the chunk: ipoi 1:nc / ichunk ics:ice
    INTEGER       :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':run_test_canopy_model'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_memory(ASSIMI_)                       !   TYPE(t_assimi_memory), POINTER      :: assimi__mem
    dsl4jsb_Def_memory(RAD_)                          !   TYPE(t_rad_memory),    POINTER      :: rad__mem
    dsl4jsb_Def_memory(VEG_)                          !   TYPE(t_veg_memory),    POINTER      :: veg__mem
    dsl4jsb_Def_memory(SOIL_)                         !   TYPE(t_soil_memory),   POINTER      :: soil__mem
    ! Declare pointers to variables in memory
    ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: net_assimilation 
    dsl4jsb_Real2D_onChunk      :: canopy_cond
    dsl4jsb_Real2D_onChunk      :: co2_conc_leaf
    ! RAD_
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tfrac_mavg_cl
    ! SOIL_
    dsl4jsb_Real2D_onChunk      :: w_soil_root_theta 
    ! VEG_
    dsl4jsb_Real2D_onChunk      :: lai 
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    dsl4jsb_Real3D_onChunk      :: lai_cl
    dsl4jsb_Real3D_onChunk      :: cumm_lai_cl
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tfrac_mavg_cl
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    ! if pftid is set in the namelist use this pftid
    IF(tile%lcts%lib_id>0)THEN
      ref_pft = tile%lcts%lib_id
    ! else use pftid 1
    ELSE
      ref_pft = 0
      tile%lcts%lib_id = 1
    ENDIF
    lctlib  => model%lctlib(tile%lcts%lib_id)
    ! Get process config
    dsl4jsb_Get_config(ASSIMI_)
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SOIL_)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, net_assimilation)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, canopy_cond)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, co2_conc_leaf)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tfrac_mavg_cl)      ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tfrac_mavg_cl)      ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_theta)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, cumm_lai_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_rub_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_et_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_pepc_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_oth_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl)     ! in

    ! define chunk for output: ipoi 1:nc / ichunk ics:ice
    ipoi   = 1
    ichunk = ics+ipoi-1

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 Init
    !!

    !> 1.1 Initalisation of met-forcing which is assumed constant
    !!
    !!
    t_air(:)            = ref_t_air + TZero
    press_srf(:)        = standard_press_srf
    q_air(:)            = calc_spec_humidity_sat(t_air(:),press_srf(:)) - ref_vpd / press_srf(:) * eps_vpd
    co2_mixing_ratio(:) = def_co2_mixing_ratio   ! def_co2_mixing_ratio is used instead of ref_co2_mixing_ratio
    rain(:)             = 0.0_wp
    snow(:)             = 0.0_wp
    swpar_srf_down(:)   = ref_swdown
    swvis_srf_down(:)   = swpar_srf_down(:) / frac_par_swvis_srf
    swnir_srf_down(:)   = swvis_srf_down(:) * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf
    swpar_fdiffuse(:)   = ref_fdiffuse
    wind_air(:)         = ref_wind
    wind_10m(:)         = wind_air(:)
    lw_srf_down(:)      = 0.0_wp
    drag_srf(:)         = 1.0_wp
    t_acoef(:)          = 0.0_wp
    t_bcoef(:)          = 0.0_wp
    q_acoef(:)          = 0.0_wp
    q_bcoef(:)          = 0.0_wp
    pch(:)              = 0.0_wp
    cos_zenith_angle(:) = ref_cos_angle
    ga(:)               = 0.4_wp / drag_srf(:) * wind_10m(:)
    !latitude                  = 50.0_wp
    !longitude                 = 0.0_wp


    CALL update_a2l( &
           tile                                            , &
           options                                         , &
           t_air                   = t_air(:)              , &
           q_air                   = q_air(:)              , &
           wind_air                = wind_air(:)           , &
           wind_10m                = wind_10m(:)           , &
           swvis_srf_down          = swvis_srf_down(:)     , &
           swnir_srf_down          = swnir_srf_down(:)     , &
           swpar_srf_down          = swpar_srf_down(:)     , &
           swpar_fdiffuse          = swpar_fdiffuse(:)     , &
           press_srf               = press_srf(:)          , &
           drag_srf                = drag_srf(:)           , &
           cos_zenith_angle        = cos_zenith_angle(:)   , &
           co2_mixing_ratio        = co2_mixing_ratio(:)   , &
           ga                      = ga(:)                 )

    !> 1.2 Initialisation of soil and vegetation (so far relevant)
    !!
    !!
    w_soil_root_theta(:) = ref_theta
    lai(:)   = ref_lai
    IF(ref_pft==0) THEN
      veg__mem%pool%leaf%nitrogen%ptr(:,iblk)     = ref_nleaf / molar_mass_N * lai(:)
    ELSE
      veg__mem%pool%leaf%nitrogen%ptr(:,iblk)     = lai(:) / lctlib%sla / lctlib%cn_leaf
      veg__mem%pool%leaf%phosphorus%ptr(:,iblk)   = veg__mem%pool%leaf%nitrogen%ptr(:,iblk) / lctlib%np_leaf
    ENDIF
    CALL calc_canopy_layers(&
                            nc                                             , & ! in
                            lctlib                                         , &
                            dsl4jsb_Config(ASSIMI_)%ncanopy                , &
                            dsl4jsb_Config(ASSIMI_)%flag_optimal_Nfraction , &
                            dsl4jsb_Config(ASSIMI_)%canopy_layer_depth(:)  , &
                            leaf_nitrogen       = veg__mem%pool%leaf%nitrogen%ptr(:,iblk)      , &
                            lai                 = lai(:)                          , &
                            ppfd_sunlit_cl      = ppfd_sunlit_tfrac_mavg_cl(:,:)  , &
                            ppfd_shaded_cl      = ppfd_shaded_tfrac_mavg_cl(:,:)  , & 
                            fleaf_sunlit_cl     = fleaf_sunlit_tfrac_mavg_cl(:,:) , &
                            lai_cl              = lai_cl(:,:), &
                            cumm_lai_cl         = cumm_lai_cl(:,:), &
                            leaf_nitrogen_cl    = leaf_nitrogen_cl(:,:), &
                            fn_rub_cl           = fn_rub_cl(:,:), &
                            fn_et_cl            = fn_et_cl(:,:), &
                            fn_pepc_cl          = fn_pepc_cl(:,:), &
                            fn_chl_cl           = fn_chl_cl(:,:), &
                            fn_oth_cl           = fn_oth_cl(:,:), &
                            inquire_n_fractions = .TRUE.)
  t_jmax_opt_mavg(:) = 20._wp+0.35_wp*(t_air(:) - Tzero)

  !> 2.0 Sensitivity plots
  !!
  !!

  !> 2.1 PAR (0-1500 µmol m-2 s-1)
  !!
  !!
  lun = 10
  WRITE(lun,*) 'PAR A GC CI'

  swpar_srf_down(:) = 0.0_wp
  DO i = 1,16

     swpar_srf_down(:) = swpar_srf_down(:) + 100.0_wp / rad2ppfd
     swvis_srf_down(:) = swpar_srf_down(:) / frac_par_swvis_srf
     swnir_srf_down(:) = swvis_srf_down(:) * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf
     CALL update_a2l(tile, &
           options            , &
           swvis_srf_down = swvis_srf_down(:), &
           swnir_srf_down = swnir_srf_down(:), &
           swpar_srf_down = swpar_srf_down(:))

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)      !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') swpar_srf_down(ipoi) * rad2ppfd, &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi) !,&
!           net_assimilation(ipoi),
!           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi),veg__mem%pool%leaf%nitrogen%ptr(ipoi,iblk)*molar_mass_N,&
!           net_assimilation(ipoi)/(veg__mem%pool%leaf%nitrogen%ptr(ipoi,iblk)*molar_mass_N)

  END DO
  swpar_srf_down(:) = ref_swdown
  swvis_srf_down(:) = swpar_srf_down(:) / frac_par_swvis_srf
  swnir_srf_down(:) = swvis_srf_down(:) * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf
  CALL update_a2l(tile, &
           options            , &
           swvis_srf_down = swvis_srf_down, &
           swnir_srf_down = swnir_srf_down, &
           swpar_srf_down = swpar_srf_down)

  !> 2.2 Temperature (-5,60 °C)
  !!
  !!
  lun = lun + 1
  WRITE(lun,*) 'T A GC CI'
  t_air(:) = - 5.0_wp + Tzero
  hlp1(:)  = t_jmax_opt_mavg(:)

  DO i = 1,70

     t_air(:) = t_air(:) + 1.0_wp
     q_air(:) = 0.8_wp * calc_spec_humidity_sat(t_air(:),press_srf(:))
     t_jmax_opt_mavg(:) = 20._wp+0.35_wp*(t_air(:) - Tzero)
     CALL update_a2l(tile, &
           options            , &
           t_air = t_air, &
           q_air = q_air)

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)      !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(5F16.7)') t_air(ipoi), &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)

  END DO
  t_air(:)    = ref_t_air + Tzero
  q_air(:)    = calc_spec_humidity_sat(t_air(:),press_srf(:)) - ref_vpd / press_srf(:) * eps_vpd
  CALL update_a2l(tile, &
        options            , &
        t_air = t_air, &
        q_air = q_air)
  t_jmax_opt_mavg(:) = hlp1(:)

  !> 2.3 atmospheric CO2 (175 - 1000 ppm)
  !!
  !!
  lun = lun + 1
  WRITE(lun,*) 'CO2 A GC CI'
  co2_mixing_ratio(:) = 150._wp
  DO i = 1,34

     co2_mixing_ratio(:) = co2_mixing_ratio(:) + 25._wp
     CALL update_a2l(tile, &
           options            , &
           co2_mixing_ratio = co2_mixing_ratio)

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)      !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') co2_mixing_ratio(ipoi), &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)

  END DO
  co2_mixing_ratio(:) = ref_co2_mixing_ratio
  CALL update_a2l(tile    , &
                       options , &
                       co2_mixing_ratio = co2_mixing_ratio)

  !> 2.4 VPD (0-3.3 kPa)
  !!
  !!
  lun = lun + 1
  WRITE(lun,*) 'VPD A GC CI'
  DO i = 1,30
     q_air(:) = calc_spec_humidity_sat(t_air(:),press_srf(:)) - i * 100.0_wp / press_srf(:) * eps_vpd
     CALL update_a2l(tile    , &
                          options , &
                          q_air = q_air)

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)      !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     vpd(:) = ( calc_spec_humidity_sat(t_air(:),press_srf(:)) - q_air(:) ) / eps_vpd * press_srf(:) / 1000.0_wp

     WRITE(lun,'(4F16.7)') vpd(ipoi), &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)

  END DO
  q_air(:) = calc_spec_humidity_sat(t_air(:),press_srf(:)) - ref_vpd / press_srf(:) * eps_vpd
  vpd(:) = ( calc_spec_humidity_sat(t_air(:),press_srf(:)) - q_air(:) ) / eps_vpd * press_srf(:) / 1000.0_wp
     CALL update_a2l(tile    , &
                          options , &
                          q_air = q_air)

  !> 2.5 Theta (0-100 %)
  !!
  !!
  lun = lun + 1
  WRITE(lun,*) 'THETA A GC CI'
  w_soil_root_theta(:) = 0.0_wp
  DO i = 1,100

     w_soil_root_theta(:) = w_soil_root_theta(:) + 0.01_wp

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)      !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') w_soil_root_theta(ipoi) *100._wp, &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)

  END DO
  w_soil_root_theta(:) = ref_theta

  !> 2.6 LAI (0-11)
  !!
  !!
  lun = lun + 1
  WRITE(lun,*) 'LAI A GC CI'
  lai(:) = 0.0_wp
! @todo test_canopy_model -> over what are you looping here ? "i" does not change any as far as I can see (Jan 20160927)
  DO i = 1,22

     lai(:) = lai(:) + 0.5_wp
     veg__mem%pool%leaf%nitrogen%ptr(:,iblk)=ref_nleaf/molar_mass_N*lai(:)
     CALL calc_canopy_layers(&
                              nc                                             , & ! in
                              lctlib                                         , &
                              dsl4jsb_Config(ASSIMI_)%ncanopy                , &
                              dsl4jsb_Config(ASSIMI_)%flag_optimal_Nfraction , &
                              dsl4jsb_Config(ASSIMI_)%canopy_layer_depth(:)  , &
                              leaf_nitrogen       = veg__mem%pool%leaf%nitrogen%ptr(:,iblk)      , &
                              lai                 = lai(:)                          , &
                              ppfd_sunlit_cl      = ppfd_sunlit_tfrac_mavg_cl(:,:)  , &
                              ppfd_shaded_cl      = ppfd_shaded_tfrac_mavg_cl(:,:)  , & 
                              fleaf_sunlit_cl     = fleaf_sunlit_tfrac_mavg_cl(:,:) , &
                              lai_cl              = lai_cl(:,:), &
                              cumm_lai_cl         = cumm_lai_cl(:,:), &
                              leaf_nitrogen_cl    = leaf_nitrogen_cl(:,:), &
                              fn_rub_cl           = fn_rub_cl(:,:), &
                              fn_et_cl            = fn_et_cl(:,:), &
                              fn_pepc_cl          = fn_pepc_cl(:,:), &
                              fn_chl_cl           = fn_chl_cl(:,:), &
                              fn_oth_cl           = fn_oth_cl(:,:), &
                              inquire_n_fractions = .TRUE.)

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)      !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') lai(ipoi), &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)

  END DO
  lai(:) = ref_lai

  !> 2.7 leaf N (0.5-3 g N / m2(LAI))
  !!
  !!
  lun = lun + 1
  WRITE(lun,*) 'Nleaf A GC CI'
  DO i = 1,24

     veg__mem%pool%leaf%nitrogen%ptr(:,iblk)=0.125_wp*(i)/molar_mass_N*lai(:)
     CALL calc_canopy_layers(&
                              nc                                             , & ! in
                              lctlib                                         , &
                              dsl4jsb_Config(ASSIMI_)%ncanopy                , &
                              dsl4jsb_Config(ASSIMI_)%flag_optimal_Nfraction , &
                              dsl4jsb_Config(ASSIMI_)%canopy_layer_depth(:)  , &
                              leaf_nitrogen       = veg__mem%pool%leaf%nitrogen%ptr(:,iblk)      , &
                              lai                 = lai(:)                          , &
                              ppfd_sunlit_cl      = ppfd_sunlit_tfrac_mavg_cl(:,:)  , &
                              ppfd_shaded_cl      = ppfd_shaded_tfrac_mavg_cl(:,:)  , & 
                              fleaf_sunlit_cl     = fleaf_sunlit_tfrac_mavg_cl(:,:) , &
                              lai_cl              = lai_cl(:,:), &
                              cumm_lai_cl         = cumm_lai_cl(:,:), &
                              leaf_nitrogen_cl    = leaf_nitrogen_cl(:,:), &
                              fn_rub_cl           = fn_rub_cl(:,:), &
                              fn_et_cl            = fn_et_cl(:,:), &
                              fn_pepc_cl          = fn_pepc_cl(:,:), &
                              fn_chl_cl           = fn_chl_cl(:,:), &
                              fn_oth_cl           = fn_oth_cl(:,:), &
                              inquire_n_fractions = .TRUE.)

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)      !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') veg__mem%pool%leaf%nitrogen%ptr(ichunk,iblk)/lai(ipoi)*molar_mass_N, &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)

  END DO
  veg__mem%pool%leaf%nitrogen%ptr(:,iblk)=ref_nleaf/molar_mass_N*lai(:)

  !> 2.8 PFT
  !!
  !!
  lun = lun + 1
  WRITE(lun,*) 'PFT A GC CI'
  DO i = 1,npft
     veg__mem%pool%leaf%nitrogen%ptr(:,iblk)  =lai(:) /model%lctlib(i)%sla /model%lctlib(i)%cn_leaf
     veg__mem%pool%leaf%phosphorus%ptr(:,iblk)=veg__mem%pool%leaf%nitrogen%ptr(:,iblk) /model%lctlib(i)%np_leaf
     ! change lctlib id
     tile%lcts%lib_id = i
     ! set pointer to the new PFT
     lctlib  => model%lctlib(tile%lcts%lib_id)
     CALL calc_canopy_layers(&
                              nc                                             , & ! in
                              lctlib                                         , &
                              dsl4jsb_Config(ASSIMI_)%ncanopy                , &
                              dsl4jsb_Config(ASSIMI_)%flag_optimal_Nfraction , &
                              dsl4jsb_Config(ASSIMI_)%canopy_layer_depth(:)  , &
                              leaf_nitrogen       = veg__mem%pool%leaf%nitrogen%ptr(:,iblk)      , &
                              lai                 = lai(:)                          , &
                              ppfd_sunlit_cl      = ppfd_sunlit_tfrac_mavg_cl(:,:)  , &
                              ppfd_shaded_cl      = ppfd_shaded_tfrac_mavg_cl(:,:)  , & 
                              fleaf_sunlit_cl     = fleaf_sunlit_tfrac_mavg_cl(:,:) , &
                              lai_cl              = lai_cl(:,:), &
                              cumm_lai_cl         = cumm_lai_cl(:,:), &
                              leaf_nitrogen_cl    = leaf_nitrogen_cl(:,:), &
                              fn_rub_cl           = fn_rub_cl(:,:), &
                              fn_et_cl            = fn_et_cl(:,:), &
                              fn_pepc_cl          = fn_pepc_cl(:,:), &
                              fn_chl_cl           = fn_chl_cl(:,:), &
                              fn_oth_cl           = fn_oth_cl(:,:), &
                              inquire_n_fractions = .TRUE.)

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)      !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(I1,3F16.7)') i, &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)

  END DO
  !STOP
  END SUBROUTINE run_test_canopy_model

END MODULE mo_land_model_interface
