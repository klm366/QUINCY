!>
!! @brief driver for the stand-alone plant growth 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
PROGRAM plant_model_driver
  USE mo_kind,                  ONLY: wp
  USE mo_base_memory
  USE mo_base_constants
  USE mo_base_config_class,     ONLY: base_config  !< TYPE, EXTENDS(t_jsb_config) :: t_base_config
  ! atmland
  USE mo_atmland_constants
  USE mo_atmland_memory
  USE mo_atmland_util
  USE mo_atmland_forcing
  USE mo_atmland_config_class,  ONLY: atmland_config
  ! land
  USE mo_land_memory,           ONLY: init_tile_variable_accumulation
  USE mo_land_constants         ! still necessary ?
  USE mo_land_util
  ! vegetation
  USE mo_veg_parameters
  USE mo_veg_radiation
  USE mo_veg_photosynthesis
  USE mo_veg_interface,         ONLY: update_vegetation
  USE mo_veg_output
  USE mo_veg_memory,            ONLY: init_plant_variable_accumulation
  USE mo_veg_util
  USE mo_veg_config_class,      ONLY: veg_config
  USE mo_veg_soil


  IMPLICIT NONE

  TYPE(VegetationParameterClass),SAVE     :: veg_param        ! data structure for vegetation parameters
  TYPE(TileClass),DIMENSION(ntile), SAVE  :: tile             ! data structure for dynamic vegetation and soil variables

  INTEGER  :: i
  REAL(wp) :: time_day                     ! displays the 'real time' in days; time_day == 1.0 after one day has passed
  INTEGER  :: time_year                    ! displays the 'real time' in years; the 1st year of a simulation == 1 (not zero)
  INTEGER  :: ntimestep
  REAL(wp) :: ref_theta  = 100._wp         ! soil water availability in %
  REAL(wp) :: ref_lai    = 3._wp           ! LAI
  REAL(wp) :: ref_nleaf  = 2.0_wp          ! N leaf per unit leaf area in gN/m2

  ! 
  ! 1.0 Read config paramters from namelist.plant 
  ! 
  base_config%namelist_filename       = 'namelist.plant'
  veg_config%namelist_filename        = 'namelist.plant'
  atmland_config%namelist_filename    = 'namelist.plant'
  CALL base_config%Init
  CALL veg_config%Init
  CALL atmland_config%Init
  ! 
  ! 1.1 Initialisation of vegetation parameters, constants and forcing
  ! 
  CALL init_base_constants
  CALL init_atmland_constants
  CALL init_veg_parameters(veg_param)
  CALL init_plant_variable_accumulation(tile)   ! set initial average values, defined in mo_veg_memory
  CALL init_tile_variable_accumulation(tile)    ! set initial average values, defined in mo_land_memory
  !CALL init_soil_variable_accumulation(tile) follows later
  CALL init_forcing   ! needs atmland_config%Init to run before
  ! 
  ! 1.2 Initialise time step
  ! 
  time_day   = 0.0_wp                      ! equals 1.0 after one day has passed (i.e. start with zero)
  time_year  = 1                           ! the 1st year of a simulation == 1 (not zero)
  dt         = base_config%dt              ! time step length in seconds, default: 1800.0 == 30 minutes
  ntimestep  = base_config%Timesteps()     ! number of time steps
  ! 
  ! 1.3 Initialise tile
  ! 
  veg_param%tjmax_opt(:)  = 20._wp+0.35_wp*(atm2land%t_air - Tzero)   ! conversion from K to °C
  tile(:)%soil%theta      = ref_theta
  tile(:)%soil%ammonium = 1.0_wp
  tile(:)%soil%phosphate  = 1.0_wp
  tile(:)%plant%lai       = 1.0_wp
  tile(1)%plant%pftid     = veg_config%pft_id
  DO i = 1,ntile
     tile(i)%plant%pool%leaf%carbon = tile(i)%plant%lai / veg_param%sla(tile(i)%plant%pftid)
     tile(i)%plant%pool%leaf%nitrogen = tile(i)%plant%pool%leaf%carbon / veg_param%cn_leaf(tile(i)%plant%pftid)
     tile(i)%plant%pool%leaf%phosphorus = tile(i)%plant%pool%leaf%nitrogen / veg_param%np_leaf(tile(i)%plant%pftid)
  ENDDO 
  tile(:)%plant%pool%labile%carbon      = 100.0_wp
  tile(:)%plant%pool%labile%nitrogen    = 10.0_wp
  tile(:)%plant%pool%labile%phosphorus  = 1.0_wp
  tile(:)%soil%ammonium = 1.0_wp
  tile(:)%soil%phosphate  = 1.0_wp
  ! 
  ! 1.4 open input files (and read the header into dummy variables); see mo_atmland_forcing
  ! 
  CALL open_forcing_input_file  ! only if: option_forcing == 'readfile'
  ! 
  ! 1.5 Initialise output and accumulation, and open output files; see mo_veg_output
  ! 
  CALL open_pool_output_files
  CALL open_flux_output_files
  ! 
  ! 2.0 time loop
  ! 
  DO i = 1,ntimestep
     ! 
     ! 2.1 advance timer
     ! 
     time_day  = time_day  + dt/one_day         ! time in days
     ! 
     ! 2.2 set orbital parameters and update time-variant forcing; 
     ! 
     CALL update_forcing(time_day)  
     ! 
     ! 3.0 canopy calculation fluxes
     ! 
     CALL update_radiation(tile)    
     ! 
     ! 4.0 update soil (biogeochemical part - here only N/P uptake)
     ! 
     CALL calc_nutrient_uptake(tile, veg_param)  
     !  
     ! 5.0 calculate vegetation growth and turnover
     ! 
     CALL update_vegetation(tile,veg_param)
     ! 
     ! 6.0 Accumulation
     ! 
     CALL update_vegetation_moving_averages(tile)
     CALL update_tile_moving_averages(tile)
     ! 
     ! 7.0 Output
     ! 
     ! 7.1 pool: write output in certain intervals, see mo_veg_output
     ! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
     ! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
     CALL write_output_pool(time_day, time_year, tile)
     ! 7.2 flux: write output in certain intervals, see mo_veg_output
     CALL write_output_flux(time_day, time_year, tile)
     ! 
     ! 8.0 check timer; adjust year counter
     ! 
     IF(time_day > one_year) THEN
       time_day  = time_day - one_year   ! largest value of time ~ 365; see mo_base_constants
       time_year = time_year + 1
     END IF

  ENDDO

  CALL close_output_files
  CALL close_forcing_input_file  ! only if: option_forcing == 'readfile'


END PROGRAM plant_model_driver
