!>
!! @brief driver for the stand-alone plant growth 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
  !! @par Initialization
  !! 1.0 Read config paramters from namelist.plant section \n
  !! 1.1 Initialisation of vegetation parameters, constants, forcing section, and daily & moving average \n
  !! 1.2 Initialise time step \n
  !! 1.3 Initialise tile \n
  !! 1.4 open input files (and read the header into dummy variables) \n
  !! 1.5 Initialise output and flux accumulation, and open output files \n
  !! @par Time loop
  !! 2.0 time loop \n
    !! 2.1 advance timer \n
    !! 2.2 set orbital parameters and update time-variant forcing
    !! @par Tile loop [nested within Time loop]
    !! 3.0 tile loop \n
    !! 3.1 set vegetation fluxes to zero \n
    !! 3.2 canopy calculation fluxes \n
    !! 3.3 update soil (biogeochemical part - here only N/P uptake) \n
    !! 3.4 calculate vegetation growth and turnover \n
    !! 3.5 update soil water bucket \n
    !! 3.6 Accumulation: daily & moving average \n
    !! 3.7 pool: write output in certain intervals \n
    !! 3.8 flux: write output in certain intervals \n
  !! 2.3 [Time loop] check timer; adjust year counter
  !! 
PROGRAM plant_model_driver
  ! shared
  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish ! message_text
  ! base
  USE mo_base_memory,           ONLY: t_slm_tile
  USE mo_base_constants
!  USE mo_base_util 
  USE mo_base_config_class,     ONLY: base_config  !< TYPE, EXTENDS(t_jsb_config) :: t_base_config
  ! atmland
  USE mo_atmland_constants
  USE mo_atmland_memory
  USE mo_atmland_util
  USE mo_atmland_forcing
  USE mo_atmland_config_class,  ONLY: atmland_config
  ! land
!  USE mo_land_constants  
!  USE mo_land_memory   
  USE mo_land_util,             ONLY: init_element_zero
  ! vegetation
  USE mo_veg_parameters
  USE mo_veg_radiation
  USE mo_veg_photosynthesis
  USE mo_veg_interface,         ONLY: update_vegetation
  USE mo_veg_output
  USE mo_veg_memory
  USE mo_veg_config_class,      ONLY: veg_config
  USE mo_veg_soil
  USE mo_veg_util

  IMPLICIT NONE

  TYPE(t_slm_veg_param),SAVE               :: veg_param        ! data structure for vegetation parameters
  TYPE(t_slm_tile),DIMENSION(ntile), SAVE  :: tile             ! data structure for dynamic vegetation and soil variables

  ! Local variables
  REAL(wp) :: time_day                     !< displays the 'real time' in days; time_day == 1.0 after one day has passed
  INTEGER  :: time_year                    !< displays the 'real time' in years; the 1st year of a simulation == 1 (not zero)
  INTEGER  :: ntimestep                    !< number of timesteps the model runs
  INTEGER  :: itimestep                    !< current timestep i, looping over timesteps ntimestep
  REAL(wp) :: timestep_prop_day            !< proportion of one day that is simulated with one timestep (dt/one_day)
  INTEGER  :: itile                        !< the current tile(i), for looping over ntile

  ! check number of tiles, exit if more than output routines can handle (for now)
  IF(ntile > 49) CALL finish('plant_model_driver','Number of tiles is too large for veg_output to work reliable; see lun')

  !! 1.0 Read config paramters from namelist.plant section
  !!
  !!  
  base_config%namelist_filename       = 'namelist.plant'
  veg_config%namelist_filename        = 'namelist.plant'
  atmland_config%namelist_filename    = 'namelist.plant'
  CALL base_config%Init
  CALL veg_config%Init
  CALL atmland_config%Init

  !! 1.1 Initialisation of vegetation parameters, constants and forcing section
  !! 
  !! 
  CALL init_base_constants
  CALL init_atmland_constants
  CALL init_element_zero                               ! mo_land_util; run before init_plant() & set_zero_fluxes_of_plant_pools()
  ! vegetation
  DO itile = 1,ntile
    CALL init_plant(tile(itile))                       ! mo_veg_util 
    CALL set_zero_fluxes_of_plant_pools(tile(itile))   ! mo_veg_util; called each timestep (see Tile loop below)
  END DO
  CALL init_veg_parameters(veg_param)                  ! mo_veg_parameters
  ! forcing
  CALL init_forcing                                !< needs atmland_config%Init to run before
  ! accumulation
  DO itile = 1,ntile
    CALL tile(itile)%plant%daily_average%init      !< initially sets all variables to zero or .FALSE.
    CALL tile(itile)%plant%moving_average%init     !< init variables with reasonable default values and sets accumulation periods
  END DO

  !! 1.2 Initialise time step
  !! 
  !! 
  time_day          = 0.0_wp                      ! equals 1.0 after one day has passed (i.e. start with zero)
  time_year         = 1                           ! the 1st year of a simulation == 1 (not zero)
  dt                = base_config%dt              ! time step length in seconds, default: 1800.0 == 30 minutes
  ntimestep         = base_config%Timesteps()     ! number of time steps
  timestep_prop_day = dt/one_day                  ! proportion of one day that is simulated with one timestep

  !! 1.3 Initialise tile
  !! 
  !! 
  veg_param%tjmax_opt(:)  = 20._wp+0.35_wp*(atm2land%t_air - Tzero)   ! conversion from K to °C
  DO itile = 1,ntile
    tile(itile)%soil%moisture                 = veg_config%soil_awc_prescribe
    tile(itile)%soil%theta                    = veg_config%soil_theta_prescribe 
    tile(itile)%soil%ammonium                 = 0.2_wp
    tile(itile)%soil%phosphate                = 0.01_wp
    tile(itile)%plant%pftid                   = veg_config%pft_id
    tile(itile)%plant%pool%reserve%carbon     = 50.0_wp
    tile(itile)%plant%pool%reserve%nitrogen   = tile(itile)%plant%pool%reserve%carbon / &
                                                (1._wp + fresp_growth) / veg_param%cn_leaf(tile(itile)%plant%pftid)
    tile(itile)%plant%pool%reserve%phosphorus = tile(itile)%plant%pool%reserve%nitrogen / &
                                                veg_param%np_leaf(tile(itile)%plant%pftid)
  END DO 

  !! 1.4 open input files (and read the header into dummy variables)
  !! 
  !! 
  CALL open_forcing_input_file  ! mo_atmland_forcing; only if: option_forcing == 'readfile'

  !! 1.5 Initialise output and flux accumulation, and open output files; see mo_veg_output
  !! 
  !! 
  CALL init_veg_output_variables          ! mo_veg_output : interval variables & flux accumulation per tile
  DO itile = 1,ntile
    CALL open_pool_output_files(itile)
    CALL open_flux_output_files(itile)
  END DO

  !! 2.0 time loop
  !! 
  !! 
  DO itimestep = 1,ntimestep

    !! 2.1 advance timer
    !! 
    !! 
    time_day  = time_day  + timestep_prop_day         ! time in days

    !! 2.2 set orbital parameters and update time-variant forcing
    !! 
    !! 
    CALL update_forcing(time_day)               !< mo_atmland_forcing

    !! 3.0 tile loop
    !! 
    !! 
    DO itile = 1,ntile

      !! 3.1 set vegetation fluxes to zero
      !! 
      !! 
      CALL set_zero_fluxes_of_plant_pools(tile(itile))  !< set to zero all fluxes of the plant pools; mo_veg_util

      !! 3.2 canopy calculation fluxes
      !! 
      !! 
      CALL update_radiation(tile(itile))  !< mo_veg_radiation

      !! 3.3 update soil (biogeochemical part - here only N/P uptake)
      !! 
      !! 
      CALL update_vegsoil_bgc(tile(itile), veg_param)  !< mo_veg_soil

      !! 3.4 calculate vegetation growth and turnover
      !! 
      !! 
      CALL update_vegetation(tile(itile), veg_param)   !< mo_veg_interface

      !! 3.5 update soil water bucket 
      !! 
      !! 
      CALL update_vegsoil_physics(tile(itile))  !< mo_veg_soil

      !! 3.6 Accumulation
      !! 
      !! daily & moving average
      CALL tile(itile)%plant%daily_average%calc_avrg(tile(itile)%plant, time_day, timestep_prop_day) !< daily avrg of the prev day
      CALL tile(itile)%plant%moving_average%update_veg(tile(itile)%plant)          !< moving averages - vegetation
      CALL tile(itile)%plant%moving_average%update_soil(tile(itile)%plant)         !< moving averages - soil -> plant N & P uptake
      CALL tile(itile)%plant%moving_average%update_tile(atm2land%t_air, atm2land%press_srf)  !< moving averages - atmland

      !! 3.7 output pool: write output in certain intervals, see mo_veg_output
      !! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
      !! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
      !! 
      CALL write_output_pool(tile(itile), veg_param, time_day, time_year, itile)

      !! 3.8 output flux: write output in certain intervals, see mo_veg_output
      !! 
      !! 
      CALL write_output_flux(tile(itile), time_day, time_year, itile)

    END DO ! end tile loop

    !! 2.3 check timer; adjust year counter
    !! 
    !! 
    IF(time_day > one_year) THEN
      time_day  = time_day - one_year   ! largest value of time ~ 365; see mo_base_constants
      time_year = time_year + 1
    END IF

  END DO ! end time loop

  !! 4.0 Close output files
  DO itile = 1,ntile
    CALL close_output_files(itile)
  END DO
    CALL close_forcing_input_file  ! only if: option_forcing == 'readfile'

END PROGRAM plant_model_driver
