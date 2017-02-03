!>
!! @brief 
!!  driver for the stand-alone land model \n
!!  kind of replacement for JSBACH4
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!!  First version by Sönke Zaehle (2015-12-04) \n
!!  Reworked by Jan Engel (2016 July)
!!
PROGRAM land_model_driver

  ! echam
  USE mo_kind
  ! land model
  USE mo_land_model_interface,    ONLY: land_model_interface 
  ! base
  USE mo_base_constants,          ONLY: ntile, one_year
  USE mo_base_memory,             ONLY: t_slm_tile !, t_model
  USE mo_base_config_class,       ONLY: base_config  !< TYPE, EXTENDS(t_jsb_config) :: t_base_config
  USE mo_jsb_model_init,          ONLY: init_model, &          ! 'init_model' is an interface
                                        time_day, time_year, ntimestep, itimestep, timestep_prop_day, &
                                        fluxnet_output_start, fluxnet_output_stop, fluxnet_output_status
  ! atmland
  USE mo_atmland_forcing,         ONLY: update_forcing, open_forcing_input_file, close_forcing_input_file ! valid for land plant canopy and soil
                                        
  ! vegetation
  USE mo_shared_output,           ONLY: init_veg_output_variables, &
                                        open_pool_output_files, &
                                        open_flux_output_files, &
                                        open_output_files_soil_forcing, &
                                        close_output_files, &
                                        close_output_files_soil_forcing



  IMPLICIT NONE

  !! 0.2 -----------------------------------------------------------------------------------------------
  !! 
  !! Local variables
  !! 
  LOGICAL                       :: l_forcing  ! IF .TRUE. THEN run forcing routines (e.g. update_forcing(), update_atm2land(), update_land2atm())
  REAL(wp), ALLOCATABLE         :: &     !! atm2land memory - INTENT(in) - start
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
                                         !! atm2land memory - end
  REAL(wp), ALLOCATABLE           :: &   !! land2atm memory - start
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
                                  sensible_hflx_ice
                                         !! land2atm memory - end
  CHARACTER(len=*), PARAMETER :: modname = 'land_model_driver'


  !! 1.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Create Tile & Model init
  !! 
  TYPE(t_slm_tile),DIMENSION(ntile), SAVE  :: tile          !< data structure for dynamic vegetation and soil variables
                                                            !! in jsbach4 tile is created in mo_jsb_model_usecases
  CALL init_model(tile(1), l_forcing)                       !< model init : mo_jsb_model_init 


  !! 2.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Forcing and Output files init
  !! 
  ! allocate forcing variables (atmland)
  !! atm2land memory
  ALLOCATE(t_air)                                            
  ALLOCATE(q_air)                                            
  ALLOCATE(rain)                                             
  ALLOCATE(snow)                                             
  ALLOCATE(wind_air)                                         
  ALLOCATE(wind_10m)                                         
  ALLOCATE(lw_srf_down)                                      
  ALLOCATE(swvis_srf_down)                                   
  ALLOCATE(swnir_srf_down)                                   
  ALLOCATE(swpar_srf_down)                                   
  ALLOCATE(swpar_fdiffuse)                                   
  ALLOCATE(press_srf)                                        
  ALLOCATE(drag_srf)                                         
  ALLOCATE(t_acoef)                                          
  ALLOCATE(t_bcoef)                                          
  ALLOCATE(q_acoef)                                          
  ALLOCATE(q_bcoef)                                          
  ALLOCATE(pch)                                              
  ALLOCATE(cos_zenith_angle)                                 
  ALLOCATE(co2_mixing_ratio)                                   
  ALLOCATE(co2_mixing_ratio_C13)                              
  ALLOCATE(co2_mixing_ratio_C14)                               
  ALLOCATE(ga)                                                 
  ! For lakes:
  ALLOCATE(albedo_lwtr)                                      
  ALLOCATE(drag_wtr)                                         
  ALLOCATE(drag_ice)                                         
  ALLOCATE(t_acoef_wtr)                                      
  ALLOCATE(t_bcoef_wtr)                                      
  ALLOCATE(q_acoef_wtr)                                      
  ALLOCATE(q_bcoef_wtr)                                      
  ALLOCATE(t_acoef_ice)                                      
  ALLOCATE(t_bcoef_ice)                                      
  ALLOCATE(q_acoef_ice)                                      
  ALLOCATE(q_bcoef_ice)
  !! land2atm memory
  ALLOCATE(t_srf)                                           
  ALLOCATE(t_eff_srf)                                       
  ALLOCATE(qsat_srf)                                        
  ALLOCATE(fact_q_air)                                      
  ALLOCATE(fact_qsat_srf)                                   
  ALLOCATE(zh_srf)                                          
  ALLOCATE(zm_srf)                                          
  ALLOCATE(alb_vis_dif)                                     
  ALLOCATE(alb_nir_dif)                                     
  ALLOCATE(alb_vis_dir)                                     
  ALLOCATE(alb_nir_dir)                                     
  ALLOCATE(albedo)                                          
  ALLOCATE(sensible_heat_flx)
  ALLOCATE(latent_heat_flx)
  ALLOCATE(net_radiation)                                   
  ALLOCATE(vpd)                                             
  ! For lakes:
  ALLOCATE(latent_hflx_wtr)                                 
  ALLOCATE(sensible_hflx_wtr)                               
  ALLOCATE(latent_hflx_ice)                                 
  ALLOCATE(sensible_hflx_ice)
  ! open forcing input files (and read the header into dummy variables); see mo_atmland_forcing
  ! l_forcing is set in init_slm() depending on the specific SLM model
  ! for SLM soil different input is used than for land plant canopy
  IF (l_forcing) THEN
    CALL open_forcing_input_file
  END IF
  ! Initialise output and flux accumulation, and open output files; see mo_shared_output
  SELECT CASE(TRIM(base_config%slm_model_name))
    CASE('land','plant')
      CALL init_veg_output_variables     ! mo_shared_output : interval variables & flux accumulation per tile
      CALL open_pool_output_files(1)     ! 1 == itile 
      CALL open_flux_output_files(1)     ! 1 == itile 
      IF (base_config%flag_generate_soil_forcing) THEN
        CALL open_output_files_soil_forcing(1)
      END IF
    CASE('soil')
      CALL open_pool_output_files(1)     ! 1 == itile 
      CALL open_flux_output_files(1)     ! 1 == itile 
!    CASE('canopy')
! 
    CASE('test_canopy')
      ntimestep = 1
   END SELECT


  !! 3.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Time loop calling the land model interface
  !! 
  DO itimestep = 1, ntimestep

    !! 3.1 advance timer
    !! 
    !! Note \n 
    !!  for each timestep the model simulates what happens during the period:
    !!  time_day (end of prevous timestep)  -->  time_day  + timestep_prop_day
    !! 
    time_day  = time_day  + timestep_prop_day         ! time in days

    !! 3.2 set orbital parameters and update time-variant forcing
    !! 
    !! l_forcing is set in init_slm() depending on the specific SLM model
    !! 
    IF (l_forcing) THEN
      CALL update_forcing(  tile(1) , & ! necessary for soil: veg%litterfall vegsoil%moisture
                            time_day, & 
                                  !! forcing (read from file and passed to atm2land memory via the land_model_interface)
                                  !! in update_forcing() defined as: REAL(wp), OPTIONAL, INTENT(out) :: 
                                  t_air                   = t_air              , &
                                  q_air                   = q_air              , &
                                  rain                    = rain               , &
                                  snow                    = snow               , &
                                  wind_air                = wind_air           , &
                                  wind_10m                = wind_10m           , &
                                  lw_srf_down             = lw_srf_down        , &
                                  swvis_srf_down          = swvis_srf_down     , &
                                  swnir_srf_down          = swnir_srf_down     , &
                                  swpar_srf_down          = swpar_srf_down     , &
                                  swpar_fdiffuse          = swpar_fdiffuse     , &
                                  press_srf               = press_srf          , &
                                  drag_srf                = drag_srf           , &
                                  t_acoef                 = t_acoef            , &
                                  t_bcoef                 = t_bcoef            , &
                                  q_acoef                 = q_acoef            , &
                                  q_bcoef                 = q_bcoef            , &
                                  pch                     = pch                , &
                                  cos_zenith_angle        = cos_zenith_angle   , &
                                  co2_mixing_ratio        = co2_mixing_ratio     , &
                                  co2_mixing_ratio_C13    = co2_mixing_ratio_C13 , &
                                  co2_mixing_ratio_C14    = co2_mixing_ratio_C14 , &
                                  ga                      = ga                   , &
                                  ! For lakes:
                                  albedo_lwtr             = albedo_lwtr        , &
                                  drag_wtr                = drag_wtr           , &
                                  drag_ice                = drag_ice           , &
                                  t_acoef_wtr             = t_acoef_wtr        , &
                                  t_bcoef_wtr             = t_bcoef_wtr        , &
                                  q_acoef_wtr             = q_acoef_wtr        , &
                                  q_bcoef_wtr             = q_bcoef_wtr        , &
                                  t_acoef_ice             = t_acoef_ice        , &
                                  t_bcoef_ice             = t_bcoef_ice        , &
                                  q_acoef_ice             = q_acoef_ice        , &
                                  q_bcoef_ice             = q_bcoef_ice          &
                                  &  )
    END IF

    
    !! 3.3 land model interface
    !! 
    CALL land_model_interface(tile(1), &
                              l_forcing, &
                              time_day, time_year, timestep_prop_day, &
                                  !! atm2land memory - INTENT(in) - start
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
                                 !! atm2land memory - end
                                 !! land2atm memory - INTENT(out) - start
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
                                  !! land2atm memory - end
                                &  ) ! land_model_interface()

    !! 3.4 check timer & adjust year counter
    !! 
    IF(time_day > one_year) THEN
      time_day  = time_day - one_year   ! largest value of time ~ 365; see mo_base_constants
      time_year = time_year + 1
    END IF

    !! 3.5 special case: fluxnet output
    !! 
    IF(base_config%fluxnet_output) THEN
      IF (itimestep == fluxnet_output_start) THEN         ! start of fluxnet output (at each timestep)
          fluxnet_output_status = 2                       ! status update: currently writing output data
          CALL close_output_files(1)                      ! 1 == itile
          base_config%output_interval_pool   = 'timestep' 
          base_config%output_interval_flux   = 'timestep' 
          CALL open_pool_output_files(1)                  ! 1 == itile 
          CALL open_flux_output_files(1)                  ! 1 == itile 

      ELSE IF (itimestep == fluxnet_output_stop) THEN     ! stop of fluxnet output; NO output after this timestep anyway
          fluxnet_output_status = 3                       ! status update: writing output data finished
          CALL close_output_files(1)                      ! 1 == itile

      END IF
    END IF

  END DO  ! end time loop 


  !! 4.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Close output and forcing files
  !! 
  IF (fluxnet_output_status == 0) CALL close_output_files(1)     ! 1 == itile 
  IF (l_forcing) CALL close_forcing_input_file ! only if forcing is enabled at all (land plant canopy soil)
  IF (base_config%flag_generate_soil_forcing) THEN
    CALL close_output_files_soil_forcing(1)
  END IF


END PROGRAM land_model_driver
