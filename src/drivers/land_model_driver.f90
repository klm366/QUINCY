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
!! @brief driver for the stand-alone QUINCY model
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!!  First version by Sönke Zaehle (2015-12-04) 
!!  Reworked by Jan Engel (multiple times)
!!
PROGRAM land_model_driver

  ! echam
  USE mo_kind
  ! land model
  USE mo_land_model_interface,    ONLY: land_model_interface 
  ! base
  USE mo_jsb_math_constants,      ONLY: one_year, eps4
  USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
  USE mo_jsb_model_init,          ONLY: jsbach_setup_grid, &      ! interface to e.g. setup_grid_standalone
                                        quincy_init, &            ! JSBACH4: jsbach_init
                                        time_day, time_year, ntimestep, itimestep, timestep_prop_day, &
                                        fluxnet_output_start, fluxnet_output_stop, fluxnet_output_status
  USE mo_jsb_class,               ONLY: get_model
  USE mo_jsb_model_class,         ONLY: t_jsb_model
  USE mo_jsb_base,                ONLY: setup_model_vgrid_procconf, quincy_proc_parameter_init, jsbach_setup_tiles
  ! forcing
  USE mo_jsb4_forcing,            ONLY: forcing_options, calc_output_time_interval_pool
  ! atmland
  USE mo_atmland_forcing,         ONLY: update_forcing, open_forcing_input_file, close_forcing_input_file ! valid for land plant canopy and soil
  ! vegetation
  USE mo_veg_memory_class,        ONLY: Deallocate_var_veg_pools_elements_a
  ! sb
  USE mo_sb_memory_class,         ONLY: Deallocate_var_sb_pools_elements_a
  ! output
  USE mo_quincy_output,           ONLY: open_pool_output_files, &
                                        open_flux_output_files, &
                                        open_output_files_soil_forcing, &
                                        close_output_files, &
                                        close_output_files_soil_forcing


  IMPLICIT NONE

  !! 0.2 -----------------------------------------------------------------------------------------------
  !! 
  !! Local variables
  !! 
  LOGICAL       :: l_forcing                          ! IF .TRUE. THEN run forcing routines (e.g. update_forcing(), update_a2l(), update_l2a())
  INTEGER       :: transient_spinup_itimestep_period  ! count the number of timesteps covering one repeated spinup period
  INTEGER       :: iblk, nproma, ics, ice, nc         !< *slmtmp may be an input from JSBACH4 later in SLM dev
  INTEGER       :: nproca                             ! number of processors for jsbach
  INTEGER       :: nprocb                             ! has to be one
  INTEGER       :: nprocio                            ! number of processors for I/O server
  INTEGER       :: npedim                             ! Working dimension for blocks in each domain.
                                                      ! Default (-1): each domain is processed in one call
  INTEGER       :: model_id
  ! output control
  LOGICAL       :: l_output, &                        ! output is enabled (flux accumulation is done and l_output_interval_match is checked for being TRUE)
                   l_output_interval_match            ! the specific output interval is a match with current time, hence output is written
  REAL(wp)      :: output_interval_pool_act, &        ! the actual interval number used at the current timestep
                   output_interval_pool, &            ! output interval of pools as a number (the same is used for flux output), converted from character
                   output_interval_pool_spinup        ! output interval of pools as a number valid for spinup period
  REAL(wp)      :: transient_spinup_timesteps         ! timesteps of spinup (transient forcing) as positive number, needed for output control
  
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
  !! start - variables for passing values from land to atmosphere
  REAL(wp), ALLOCATABLE, DIMENSION(:)  :: & 
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
                                  albedo            ,   &   ! not present in jsb4 | RAD_
                                  sensible_heat_flx ,   &   ! not present in jsb4 | 
                                  latent_heat_flx   ,   &   ! not present in jsb4 | 
                                  net_radiation     ,   &   ! not present in jsb4 | RAD_
                                  vpd               ,   &   ! not present in jsb4 | 
                                  ! For lakes:
                                  latent_hflx_wtr   ,   &
                                  sensible_hflx_wtr ,   &
                                  latent_hflx_ice   ,   &
                                  sensible_hflx_ice
  !! end - variables for passing values from land to atmosphere
  CHARACTER(len=*), PARAMETER :: modname = 'land_model_driver'



  !! 1.0 -----------------------------------------------------------------------------------------------
  !! 
  !! setup model, configure and init processes and grids
  !! 
  TYPE(t_jsb_model)         , POINTER :: model  !< local pointer to model
  CLASS(t_jsb_tile_abstract), POINTER :: tile   !< local data structure for dynamic vegetation and soil variables
                                                !! in jsbach4 tile is created in mo_jsb_model_usecases
  !! 1.1  define default values
  nproca   = 1
  nprocb   = 1
  nprocio  = 0
  npedim   = -1
  model_id = 1      ! must be consistent with model_idx set in setup_model_vgrid_procconf()
  iblk     = 1
  nproma   = 1      ! note: nproma may not be smaller than ice; this is differently in jsbach4, and only a restriction to SLM !
  ics      = 1      ! 
  ice      = 1      ! note: ice must be smaller or equal to nproma; this is differently in jsbach4, and only a restriction to SLM !
  nc       = ice - ics + 1                 ! (end-number - start-number + 1)

  !! 1.2  set up the model with vgrid and process configs (jsbach%models(model_id)%m => model  with 'model' of TYPE(t_jsb_model))
  CALL setup_model_vgrid_procconf("namelist.slm")
  model => get_model(model_id)

  !! 1.3 init of horizontal grid
  CALL jsbach_setup_grid(model%id, "hgrid0")    ! model%id, horizontal grid name

  !! 1.4 init process parameters (and modify parameter values for param-sensi)
  CALL quincy_proc_parameter_init(model%id)

  !! 1.5 set up tiles, read lctlib (and modify parameter values for param-sensi) and init usecase
  CALL jsbach_setup_tiles(model%id)

  !! 1.6 quincy-model init
  !!     and init of: forcing, SB_ & VEG_ pools, tiles, time, fluxnet diagnostics
  CALL quincy_init(model%id, iblk, nproma, ics, ice, nc, &    ! in
                   l_forcing)                                 ! out
  tile => model%Get_top_tile()



  !! 2.0 -----------------------------------------------------------------------------------------------
  !! 
  !! allocate forcing variables (atmland)
  !! 
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

  !! 2.1 -----------------------------------------------------------------------------------------------
  !! 
  !! open forcing files
  !! 
  ! open forcing input files (and read the header into dummy variables); see mo_atmland_forcing
  ! l_forcing is set in init_slm() depending on the specific SLM model
  ! for SLM soil different input is used than for land plant canopy
  IF (l_forcing) THEN
    CALL open_forcing_input_file(TRIM(model%config%slm_model_name))
  END IF


  !! 3.0 -----------------------------------------------------------------------------------------------
  !! 
  !! settings for the static/transient forcing_mode
  !! 
  ! get output interval of pools as a number
  output_interval_pool_act    = 0.0_wp
  output_interval_pool        = calc_output_time_interval_pool(model%config%output_interval_pool)
  output_interval_pool_spinup = calc_output_time_interval_pool(model%config%output_interval_pool_spinup)
  ! init with default values
  model%config%l_transient_forcing_mode          = .FALSE.
  model%config%l_transient_spinup                = .FALSE.
  model%config%l_transient_spinup_rewind_forcing = .FALSE.
  model%config%l_static_forcing_mode             = .FALSE.
  itimestep                         = 0         ! default start of time loop
  transient_spinup_itimestep_period = 0  
  transient_spinup_timesteps        = 0         ! is used only during l_transient_spinup
  ! set simulation values
  IF ( TRIM(forcing_options%forcing_mode) == "static" ) THEN
    model%config%l_static_forcing_mode     = .TRUE.
    output_interval_pool_act               = output_interval_pool
  ELSEIF ( TRIM(forcing_options%forcing_mode) == "transient" ) THEN
    ! set logicals to TRUE
    model%config%l_transient_forcing_mode  = .TRUE.
    model%config%l_transient_spinup        = .TRUE.  ! transient forcing_mode starts with spinup; is set to FALSE during runtime when spinup is done
    ! output interval for spinup
    output_interval_pool_act          = output_interval_pool_spinup
    ! number of spinup timesteps
    itimestep                         = itimestep - forcing_options%transient_spinup_timesteps_all
    transient_spinup_timesteps        = itimestep * (-1.0_wp)
    ! init spinup period counter
    transient_spinup_itimestep_period = 0
  ENDIF


  !! 4.0 -----------------------------------------------------------------------------------------------
  !! 
  !! open output files
  !! need to be done after 'output_interval_pool_act' is defined above
  !! test_canopy & test_radiation: special case, runs for just one timestep and pool/flux output routines are not called
  !! 
  ! Initialise output and flux accumulation, and open output files; see mo_quincy_output
  ! the output routine for SLM canopy is located in mo_land_model_interface -> canopy_model_write_data()
  SELECT CASE(TRIM(model%config%slm_model_name))
    CASE('land','plant','soil','canopy')
      ! open output files for static more or for transient spinup if output may be created (i.e. output interval is not 'none')
      IF (model%config%l_static_forcing_mode) THEN
        CALL open_pool_output_files(tile, iblk, ics, ice, nc, TRIM(model%config%output_interval_pool), '')
        CALL open_flux_output_files(tile, iblk, ics, ice, nc, TRIM(model%config%output_interval_flux), '')
      ELSEIF (model%config%l_transient_forcing_mode .AND. output_interval_pool_act > 0.0_wp) THEN
        CALL open_pool_output_files(tile, iblk, ics, ice, nc, TRIM(model%config%output_interval_pool_spinup), 'spinup_')
        CALL open_flux_output_files(tile, iblk, ics, ice, nc, TRIM(model%config%output_interval_flux_spinup), 'spinup_')
      ENDIF
      ! create soil forcing (output is created for a year after spinup; control is in mo_shared:output:write_soil_forcing)
      IF (model%config%flag_generate_soil_forcing) THEN
        CALL open_output_files_soil_forcing(tile, iblk, ics, ice, nc)
      ENDIF
    CASE('test_canopy','test_radiation')
      ntimestep = 0     ! simulation starts with zero, and loop asks for 'IF (itimestep .GT. ntimestep) EXIT'
  END SELECT


  !! 5.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Time loop calling the land model interface
  !! 
  DO
    IF (itimestep .GT. ntimestep) EXIT
    
    !! 5.1 advance simulation-timer and check timer for new year (year counter)
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
      time_day  = time_day - one_year   ! largest value of time ~ 365; see mo_jsb_math_constants
      time_year = time_year + 1
    END IF


    !! 5.2 advance timestep of this loop
    !! init itimestep with 0, the timestep 1 means running from 0 to 1, equally to time_day (sec 4.1)
    itimestep = itimestep + 1   ! increment DO loop counter


    !! 5.3 control for transient forcing_mode
    !! 
    !! 
    IF (model%config%l_transient_forcing_mode) THEN
      model%config%l_transient_spinup_rewind_forcing = .FALSE.
      ! control spinup
      IF (model%config%l_transient_spinup) THEN
        transient_spinup_itimestep_period = transient_spinup_itimestep_period + 1
        ! if one spinup period (i.e. repetitive part of the spinup) has finished
        IF (transient_spinup_itimestep_period .GT. forcing_options%transient_spinup_timesteps_repeat) THEN
          transient_spinup_itimestep_period              = 1
          model%config%l_transient_spinup_rewind_forcing = .TRUE.
        END IF
        ! if spinup is completed and the simulation should start
        IF (itimestep .EQ. 1) THEN
          model%config%l_transient_spinup_rewind_forcing = .TRUE.
          model%config%l_transient_spinup                = .FALSE.           ! do not enter the "control spinup" IF statement anymore
          time_day                                       = timestep_prop_day ! equals 1.0 at the last simulated timestep of the day
          time_year                                      = 1                 ! the 1st year of a simulation == 1 (not zero)
          ! close spinup output files if output was created (i.e. spinup output interval is not 'none')
          IF (output_interval_pool_act > 0.0_wp) THEN
            CALL close_output_files
          ENDIF
          ! change output interval
          output_interval_pool_act                       = output_interval_pool
          ! open new output files
          CALL open_pool_output_files(tile, iblk, ics, ice, nc, TRIM(model%config%output_interval_pool), '')
          CALL open_flux_output_files(tile, iblk, ics, ice, nc, TRIM(model%config%output_interval_flux), '')
        END IF
      END IF
    END IF


    !! 5.4 output control
    !!  does consider:  static & transient
    !!                  spinup & actual simulation
    !!  
    l_output                = .FALSE.      ! is generally set to FALSE tested below if TRUE at this timestep
    !! static
    IF (model%config%l_static_forcing_mode) THEN
      ! output period
      IF (time_year >= model%config%output_start_first_day_year .AND. time_year <= model%config%output_end_last_day_year) THEN
        ! might be set back to FALSE if fluxnet_output_status == 3 (i.e. fluxnet diagnostics output finished)
        l_output = .TRUE.
      ENDIF        
      ! fluxnet diagnostics output
      IF (model%config%fluxnet_output) THEN
        ! activate output during the last full forcing period (within 'ntimestep'), independent of the general output 'output_start_..' & 'output_end_...' options
        IF (fluxnet_output_status < 3.0_wp) THEN
          l_output = .TRUE.
        ! set l_output (back) to FALSE if the last full forcing period is finished (after the last full forcing period no further output is created and fluxnet_output_status = 3)
        ELSE
          l_output = .FALSE.
        ENDIF
      ENDIF
    !! transient
    ELSEIF (model%config%l_transient_forcing_mode) THEN
      ! spinup
      IF (model%config%l_transient_spinup) THEN
        ! if output may be generated during spinup (i.e. output interval is not 'none')
        IF (output_interval_pool_act > 0.0_wp) THEN
          l_output = .TRUE.
        ENDIF
      ! actual simulation after spinup
      ELSE
      ! differ between fluxnet diagnostics output and other transient runs
      ! fluxnet is special because the last full forcing period (within 'ntimestep') should be written out at timestep interval, and this period differs across sites
        ! fluxnet diagnostics - after spinup output is always created but at different intervals 
        IF (model%config%fluxnet_output) THEN
          ! if active and not finished yet (finished == after the last full forcing period no further output is created)
          IF (fluxnet_output_status < 3.0_wp) THEN
            l_output = .TRUE.
          ENDIF
        ! any other than fluxnet diagnostics, create output for a given start/end output period
        ELSE
          ! output period
          IF (time_year >= model%config%output_start_first_day_year .AND. time_year <= model%config%output_end_last_day_year) THEN
            l_output = .TRUE.
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    !!  match output interval
    !! 
    l_output_interval_match = .FALSE.      ! is generally set to FALSE tested below if TRUE at this timestep
    IF (l_output) THEN
      ! if interval is smaller than 'daily' then set to TRUE always, because then it can only be timestep
      IF (output_interval_pool_act < 1.0_wp) THEN
        l_output_interval_match = .TRUE.
      ! transient spinup
      ! during spinup the timesteps are a negative number, which is compensated by adding the spinup length to time_day
      ELSEIF (model%config%l_transient_spinup) THEN
        IF (MODULO(time_day + transient_spinup_timesteps, output_interval_pool_act) < eps4 .OR. &
          ABS(output_interval_pool_act - MODULO(time_day + transient_spinup_timesteps, output_interval_pool_act)) < eps4) THEN
          l_output_interval_match = .TRUE.
        ENDIF
      ! static simulation or transient simulation after spinup
      ELSE
        IF (MODULO(time_day, output_interval_pool_act) < eps4 .OR. &
          ABS(output_interval_pool_act - MODULO(time_day, output_interval_pool_act)) < eps4) THEN
          l_output_interval_match = .TRUE.
        ENDIF
      ENDIF
    ENDIF

    
    !! 5.5 set orbital parameters and update time-variant forcing
    !! 
    !! l_forcing is set in init_slm() depending on the specific SLM model
    !! 
    IF (l_forcing) THEN
      CALL update_forcing(  tile      , & ! necessary for soil: veg%litterfall vegsoil%moisture
                            ics       , &
                            ice       , &
                            nc        , &
                            iblk      , &
                            itimestep , &
                            time_day  , & 
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

    
    !! 5.6 land model interface
    !! 
    CALL land_model_interface(tile, &
                              l_forcing, &
                              l_output, &
                              l_output_interval_match, &
                              output_interval_pool_act, &  ! used for flux output: output_interval_flux
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
                                  !! start - variables for passing values from land to atmosphere
                                  !! INTENT(out)
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
                                  !! end - variables for passing values from land to atmosphere
                                &  ) ! land_model_interface()

    !! 5.7 special case: fluxnet output
    !! 
    IF(model%config%fluxnet_output) THEN
      IF (itimestep == fluxnet_output_start) THEN         ! start of fluxnet output (at each timestep)
          fluxnet_output_status = 2                       ! status update: currently writing output data
          CALL close_output_files 
          model%config%output_interval_pool   = 'timestep' 
          model%config%output_interval_flux   = 'timestep' 
          output_interval_pool                = calc_output_time_interval_pool(model%config%output_interval_pool)
          output_interval_pool_act            = output_interval_pool
          CALL open_pool_output_files(tile, iblk, ics, ice, nc, TRIM(model%config%output_interval_pool), '')
          CALL open_flux_output_files(tile, iblk, ics, ice, nc, TRIM(model%config%output_interval_flux), '')

      ELSE IF (itimestep == fluxnet_output_stop) THEN     ! stop of fluxnet output; NO output after this timestep anyway
          fluxnet_output_status = 3                       ! status update: writing output data finished
          CALL close_output_files

      END IF
    END IF

  END DO  ! end time loop 


  !! 6.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Close output and forcing files
  !! 
  IF (fluxnet_output_status == 0) CALL close_output_files
  IF (l_forcing) CALL close_forcing_input_file(TRIM(model%config%slm_model_name)) ! only if forcing is enabled at all (land plant canopy soil)
  IF (model%config%flag_generate_soil_forcing) THEN
    CALL close_output_files_soil_forcing
  END IF


  !! 7.0 -----------------------------------------------------------------------------------------------
  !! 
  !! DEALLOCATE variables
  CALL Deallocate_var_sb_pools_elements_a     ! mo_sb_memory_class
  CALL Deallocate_var_veg_pools_elements_a    ! mo_veg_memory_class

END PROGRAM land_model_driver
