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
  USE mo_base_constants,          ONLY: one_year
  USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
  USE mo_jsb_model_init,          ONLY: init_model, &          ! 'init_model' is an interface
                                        time_day, time_year, ntimestep, itimestep, timestep_prop_day, &
                                        fluxnet_output_start, fluxnet_output_stop, fluxnet_output_status
  USE mo_jsb_class,               ONLY: get_model
  USE mo_jsb_model_class,         ONLY: t_jsb_model
  USE mo_jsb_base,                ONLY: init_base
  ! forcing
  USE mo_jsb4_forcing,            ONLY: forcing_options
  ! atmland
  USE mo_atmland_forcing,         ONLY: update_forcing, open_forcing_input_file, close_forcing_input_file ! valid for land plant canopy and soil
  ! vegetation
  USE mo_shared_output,           ONLY: open_pool_output_files, &
                                        open_flux_output_files, &
                                        open_output_files_soil_forcing, &
                                        close_output_files, &
                                        close_output_files_soil_forcing
  USE mo_veg_memory_class,        ONLY: Deallocate_var_veg_pools_elements_a
  ! sb
  USE mo_sb_memory_class,         ONLY: Deallocate_var_sb_pools_elements_a


  IMPLICIT NONE

  !! 0.2 -----------------------------------------------------------------------------------------------
  !! 
  !! Local variables
  !! 
  LOGICAL       :: l_forcing                          ! IF .TRUE. THEN run forcing routines (e.g. update_forcing(), update_a2l(), update_l2a())
  INTEGER       :: transient_spinup_itimestep_period  ! count the number of timesteps covering one repeated spinup period
  INTEGER       :: iblk, nproma, ics, ice, nc         !< *slmtmp may be an input from JSBACH4 later in SLM dev
  INTEGER       :: nsoil_sb                           !< necessary for open_pool_output_files() & init_output_writing_pool()
                                                        !!  number of soil layers used in SB_ - equal to nsoil_water
  INTEGER       :: nproca                             ! number of processors for jsbach
  INTEGER       :: nprocb                             ! has to be one
  INTEGER       :: nprocio                            ! number of processors for I/O server
  INTEGER       :: npedim                             ! Working dimension for blocks in each domain.
                                                      ! Default (-1): each domain is processed in one call
  INTEGER       :: model_id
  REAL(wp), ALLOCATABLE, DIMENSION(:)  :: &     !! a2l memory - INTENT(in) - start
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
                                  q_bcoef_ice
                                         !! a2l memory - end
  REAL(wp), ALLOCATABLE, DIMENSION(:)  :: &   !! l2a memory - start
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
                                         !! l2a memory - end
  CHARACTER(len=*), PARAMETER :: modname = 'land_model_driver'



  !! 1.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Create Tile and Model & run Init & Define iblk, ics, ice, ...
  !! 
  TYPE(t_jsb_model)         , POINTER :: model  !< local pointer to model
  CLASS(t_jsb_tile_abstract), POINTER :: tile   !< local data structure for dynamic vegetation and soil variables
                                                !! in jsbach4 tile is created in mo_jsb_model_usecases
  !! 1.1  define default values
  nproca   = 1
  nprocb   = 1
  nprocio  = 0
  npedim   = -1
  model_id = 1
  iblk     = 1
  nproma   = 1      ! note: nproma may not be smaller than ice; this is differently in jsbach4, and only a restriction to SLM !
  ics      = 1      ! 
  ice      = 1      ! note: ice must be smaller or equal to nproma; this is differently in jsbach4, and only a restriction to SLM !
  nc       = ice - ics + 1                 ! (end-number - start-number + 1)
  nsoil_sb = -1                            ! namelist parameter

  !! 1.2  create new model: jsbach%models(model_id)%m => current    ('current' is of TYPE(t_jsb_model))
  CALL init_base() 
  
  !! 1.3  model init
  CALL init_model(nproca, nprocb, npedim, model_id, l_forcing, nsoil_sb, iblk, nproma, ics, ice, nc)

  !! 1.4  get model and tile
  model => get_model(model_id)
  tile  => model%Get_top_tile()


  !! 2.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Forcing and Output files init
  !! 
  ! allocate forcing variables (atmland)
  !! a2l memory
  ALLOCATE(t_air(nc))                                           
  ALLOCATE(q_air(nc))                                           
  ALLOCATE(rain(nc))                                            
  ALLOCATE(snow(nc))                                            
  ALLOCATE(wind_air(nc))                                        
  ALLOCATE(wind_10m(nc))                                        
  ALLOCATE(lw_srf_down(nc))                                     
  ALLOCATE(swvis_srf_down(nc))                                  
  ALLOCATE(swnir_srf_down(nc))                                  
  ALLOCATE(swpar_srf_down(nc))                                  
  ALLOCATE(swpar_fdiffuse(nc))                                  
  ALLOCATE(press_srf(nc))                                       
  ALLOCATE(drag_srf(nc))                                        
  ALLOCATE(t_acoef(nc))                                         
  ALLOCATE(t_bcoef(nc))                                         
  ALLOCATE(q_acoef(nc))                                         
  ALLOCATE(q_bcoef(nc))                                         
  ALLOCATE(pch(nc))                                             
  ALLOCATE(cos_zenith_angle(nc))                                
  ALLOCATE(co2_mixing_ratio(nc))                                  
  ALLOCATE(co2_mixing_ratio_C13(nc))                             
  ALLOCATE(co2_mixing_ratio_C14(nc))                              
  ALLOCATE(ga(nc))                                                
  ALLOCATE(nhx_deposition(nc))                                                
  ALLOCATE(noy_deposition(nc))      
  ALLOCATE(p_deposition(nc))                                          
  ALLOCATE(nhx_n15_deposition(nc))                                                
  ALLOCATE(noy_n15_deposition(nc))                                                
  ! For lakes:
  ALLOCATE(albedo_lwtr(nc))                                     
  ALLOCATE(drag_wtr(nc))                                        
  ALLOCATE(drag_ice(nc))                                        
  ALLOCATE(t_acoef_wtr(nc))                                     
  ALLOCATE(t_bcoef_wtr(nc))                                     
  ALLOCATE(q_acoef_wtr(nc))                                     
  ALLOCATE(q_bcoef_wtr(nc))                                     
  ALLOCATE(t_acoef_ice(nc))                                     
  ALLOCATE(t_bcoef_ice(nc))                                     
  ALLOCATE(q_acoef_ice(nc))                                     
  ALLOCATE(q_bcoef_ice(nc))
  !! l2a memory
  ALLOCATE(t_srf(nc))                                          
  ALLOCATE(t_eff_srf(nc))                                      
  ALLOCATE(qsat_srf(nc))                                       
  ALLOCATE(fact_q_air(nc))                                     
  ALLOCATE(fact_qsat_srf(nc))                                  
  ALLOCATE(zh_srf(nc))                                         
  ALLOCATE(zm_srf(nc))                                         
  ALLOCATE(alb_vis_dif(nc))                                    
  ALLOCATE(alb_nir_dif(nc))                                    
  ALLOCATE(alb_vis_dir(nc))                                    
  ALLOCATE(alb_nir_dir(nc))                                    
  ALLOCATE(albedo(nc))                                         
  ALLOCATE(sensible_heat_flx(nc))
  ALLOCATE(latent_heat_flx(nc))
  ALLOCATE(net_radiation(nc))                                  
  ALLOCATE(vpd(nc))                                            
  ! For lakes:
  ALLOCATE(latent_hflx_wtr(nc))                                
  ALLOCATE(sensible_hflx_wtr(nc))                              
  ALLOCATE(latent_hflx_ice(nc))                                
  ALLOCATE(sensible_hflx_ice(nc))
  ! open forcing input files (and read the header into dummy variables); see mo_atmland_forcing
  ! l_forcing is set in init_slm() depending on the specific SLM model
  ! for SLM soil different input is used than for land plant canopy
  IF (l_forcing) THEN
    CALL open_forcing_input_file(TRIM(model%config%slm_model_name))
  END IF
  ! Initialise output and flux accumulation, and open output files; see mo_shared_output
  ! the output routine for SLM canopy is located in mo_land_model_interface -> canopy_model_write_data()
  SELECT CASE(TRIM(model%config%slm_model_name))
    CASE('land','plant','soil','canopy')
      CALL open_pool_output_files(tile, iblk, ics, ice, nc, nsoil_sb)
      CALL open_flux_output_files(tile, iblk, ics, ice, nc, nsoil_sb)
      IF (model%config%flag_generate_soil_forcing) THEN
        CALL open_output_files_soil_forcing(tile, iblk, ics, ice, nc, nsoil_sb)
      END IF
    CASE('test_canopy')
      ntimestep = 0     ! simulation starts with zero, and loop asks for 'IF (itimestep .GT. ntimestep) EXIT'
  END SELECT


  !! 3.0 -----------------------------------------------------------------------------------------------
  !! 
  !! settings for the transient forcing_mode
  !! 
  ! init with default values
  model%config%l_transient_forcing_mode          = .FALSE.
  model%config%l_transient_spinup                = .FALSE.
  model%config%l_transient_spinup_rewind_forcing = .FALSE.
  itimestep                         = 0         ! default start of time loop
  transient_spinup_itimestep_period = 0  
  ! set simulation values
  IF ( TRIM(forcing_options%forcing_mode) == "transient" ) THEN
    ! set logicals to TRUE
    model%config%l_transient_forcing_mode  = .TRUE.
    model%config%l_transient_spinup        = .TRUE.  ! transient forcing_mode always starts with a spinup
    ! number of spinup timesteps
    itimestep = itimestep - forcing_options%transient_spinup_timesteps_all
    ! init spinup period counter
    transient_spinup_itimestep_period = 0
  END IF


  !! 4.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Time loop calling the land model interface
  !! 
  DO
    IF (itimestep .GT. ntimestep) EXIT
    
    !! 4.1 advance simulation-timer and check timer for new year (year counter)
    !! 
    !! for timing data of the forcing input see forcing_options\%forcing_*
    !!
    !! Note: \n 
    !!  for each timestep the model simulates what happens during the period:
    !!  time_day (end of prevous timestep)  -->  time_day  + timestep_prop_day
    !!  time_day is initialized with 0.0
    !! 
    time_day  = time_day  + timestep_prop_day         ! time in days
    IF(time_day > one_year) THEN
      time_day  = time_day - one_year   ! largest value of time ~ 365; see mo_base_constants
      time_year = time_year + 1
    END IF


    !! 4.2 advance timestep of this loop
    !! init itimestep with 0, the timestep 1 means running from 0 to 1, equally to time_day (sec 4.1)
    itimestep = itimestep + 1   ! increment DO loop counter


    !! 4.3 control variables for transient forcing_mode
    !! 
    !! 
    IF (model%config%l_transient_forcing_mode) THEN
      model%config%l_transient_spinup_rewind_forcing = .FALSE.
      ! control spinup
      IF (model%config%l_transient_spinup) THEN
        transient_spinup_itimestep_period = transient_spinup_itimestep_period + 1
        ! if one spinup period has finished
        IF (transient_spinup_itimestep_period .GT. forcing_options%transient_spinup_timesteps_repeat) THEN
          transient_spinup_itimestep_period             = 1
          model%config%l_transient_spinup_rewind_forcing = .TRUE.
        END IF
        ! if complete spinup has ended and the simulation should start
        IF (itimestep .EQ. 1) THEN
          model%config%l_transient_spinup_rewind_forcing = .TRUE.
          model%config%l_transient_spinup                = .FALSE.           ! do not enter the "control spinup" IF statement anymore
          time_day                                      = timestep_prop_day ! equals 1.0 at the last simulated timestep of the day
          time_year                                     = 1                 ! the 1st year of a simulation == 1 (not zero)
        END IF
      END IF
    END IF
    
    
    !! 4.4 set orbital parameters and update time-variant forcing
    !! 
    !! l_forcing is set in init_slm() depending on the specific SLM model
    !! 
    IF (l_forcing) THEN
      CALL update_forcing(  tile , & ! necessary for soil: veg%litterfall vegsoil%moisture
                            ics     , &
                            ice     , &
                            nc      , &
                            iblk    , &
                            nsoil_sb, &
                            time_day, & 
                                  !! forcing (read from file and passed to a2l memory via the land_model_interface)
                                  !! in update_forcing() defined as: REAL(wp), OPTIONAL, INTENT(out) :: VAR(:)  ! assumed shape
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
                                  &  )
    END IF

    
    !! 4.5 land model interface
    !! 
    CALL land_model_interface(tile, &
                              l_forcing, &
                              time_day, &
                              time_year, &
                              timestep_prop_day, &
                                  !! a2l memory - INTENT(in) - start
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
                                  q_bcoef_ice(:)        , &
                                 !! a2l memory - end
                                 !! l2a memory - INTENT(out) - start
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
                                  sensible_hflx_ice(:)     &
                                  !! l2a memory - end
                                &  ) ! land_model_interface()

    !! 4.6 special case: fluxnet output
    !! 
    IF(model%config%fluxnet_output) THEN
      IF (itimestep == fluxnet_output_start) THEN         ! start of fluxnet output (at each timestep)
          fluxnet_output_status = 2                       ! status update: currently writing output data
          CALL close_output_files 
          model%config%output_interval_pool   = 'timestep' 
          model%config%output_interval_flux   = 'timestep' 
          CALL open_pool_output_files(tile, iblk, ics, ice, nc, nsoil_sb)
          CALL open_flux_output_files(tile, iblk, ics, ice, nc, nsoil_sb)

      ELSE IF (itimestep == fluxnet_output_stop) THEN     ! stop of fluxnet output; NO output after this timestep anyway
          fluxnet_output_status = 3                       ! status update: writing output data finished
          CALL close_output_files

      END IF
    END IF

  END DO  ! end time loop 


  !! 5.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Close output and forcing files
  !! 
  IF (fluxnet_output_status == 0) CALL close_output_files
  IF (l_forcing) CALL close_forcing_input_file(TRIM(model%config%slm_model_name)) ! only if forcing is enabled at all (land plant canopy soil)
  IF (model%config%flag_generate_soil_forcing) THEN
    CALL close_output_files_soil_forcing
  END IF


  !! 6.0 -----------------------------------------------------------------------------------------------
  !! 
  !! DEALLOCATE variables
  CALL Deallocate_var_sb_pools_elements_a     ! mo_sb_memory_class
  CALL Deallocate_var_veg_pools_elements_a    ! mo_veg_memory_class

END PROGRAM land_model_driver
