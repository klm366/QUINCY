!!
!! @brief driver for the stand-alone plant growth 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
PROGRAM plant_model_driver
  USE mo_kind,               ONLY: wp
  USE mo_base_memory
  USE mo_base_constants
  USE mo_atmland_constants
  USE mo_land_constants
  USE mo_land_memory,        ONLY: init_tile_variable_accumulation
  USE mo_atmland_util
  USE mo_atmland_forcing
  USE mo_veg_parameters
  USE mo_veg_radiation
  USE mo_veg_photosynthesis
  USE mo_veg_interface,      ONLY: update_vegetation
  USE mo_veg_output
  USE mo_veg_memory,         ONLY: init_plant_variable_accumulation
  USE mo_veg_util
  USE mo_namelist,           ONLY: POSITIONED, open_nml, position_nml  ! read model options from file


  IMPLICIT NONE

  TYPE(VegetationParameterClass),SAVE :: veg_param           ! data structure for vegetation parameters
  TYPE(TileClass),DIMENSION(ntile), SAVE :: tile             ! data structure for dynamic vegetation and soil variables

  REAL(wp) :: &                                        ! model forcing
      t_air            , &
      q_air            , &
      rain             , &
      snow             , &
      wind_air         , &
      wind_10m         , &
      lw_srf_down      , &
      swvis_srf_down   , &
      swnir_srf_down   , &
      swpar_srf_down   , &
      swpar_fdiffuse   , &
      press_srf        , &
      drag_srf         , &
      t_acoef          , &
      t_bcoef          , &
      q_acoef          , &
      q_bcoef          , &
      pch              , &
      cos_zenith_angle , &
      co2_mixing_ratio
  REAL(wp) :: &                                       ! tentative model output for atmosphere from radiation routine
      alb_vis_dif ,   &
      alb_nir_dif ,   &
      alb_vis_dir ,   &
      alb_nir_dir ,   &
      ga, &
      vpd
  REAL(wp) :: time                         ! displays the 'real time' in days; time == 1.0 after one day has passed
  INTEGER :: ntimestep
  REAL(wp) :: sw_toa_down                  ! potential radiation at the top of the atmosphere
  REAL(wp) :: longitude, latitude          ! coordinates of point to be simulated
  REAL(wp) :: annual_t_air_min, annual_t_air_max         ! annual minimum and maximum daily t_air
  REAL(wp) :: diurnal_temp_range           ! diurnal temperature range

  REAL(wp) :: ref_t_air = 25.0_wp          ! temperature in deg C
  REAL(wp) :: ref_vpd = 1000._wp           ! vpd in Pa
  REAL(wp) :: ref_rh = 0.8_wp              ! relative humidity in fractions
  REAL(wp) :: ref_swdown = 652.173913_wp   ! PAR in W m-2
  REAL(wp) :: ref_fdiffuse = 0.2_wp        ! diffuse fraction
  REAL(wp) :: ref_cos_angle = 1.0_wp       ! cosine of zenith angle
  REAL(wp) :: ref_co2 = 380.0_wp           ! co2 in ppm
  REAL(wp) :: ref_wind = 1.0_wp            ! wind in m/s
  REAL(wp) :: ref_theta = 100._wp          ! soil water availability in %
  REAL(wp) :: ref_lai = 3._wp              ! LAI
  REAL(wp) :: ref_nleaf = 2.0_wp           ! N leaf per unit leaf area in gN/m2

  INTEGER :: i

  ! variables for reading model-options from file (namelist.plant)
  INTEGER           :: f_unit      = 0
  INTEGER           :: read_status = 0
  INTEGER, SAVE     :: nml_unit    = 0     ! copied from JSBACH3 mo_jsbach.f90 | Unit of the jsbach namelist
  ! define variables, and set their default values, used for read model-options from file, and passed to options% ... then
  CHARACTER(8) :: read_output_interval_pool  = 'monthly'   ! used for read from file
  CHARACTER(8) :: read_output_interval_flux  = 'monthly'   ! 
  INTEGER      :: output_start_date          = 1
  INTEGER      :: output_end_date            = 100
  LOGICAL      :: flag_optimal_Nfraction     = .FALSE.  ! .TRUE. .FALSE.
  LOGICAL      :: flag_optimal_Nleaf         = .FALSE.  ! .TRUE. .FALSE.

  ! 
  ! 1.0 Initialisation of vegetation parameters and constants
  ! 
  CALL init_base_constants
  CALL init_atmland_constants
  CALL init_pft_param_table(veg_param)
  CALL init_veg_parameters
  CALL init_plant_variable_accumulation(tile)   ! set initial average values, defined in mo_veg_memory
  CALL init_tile_variable_accumulation(tile)    ! set initial average values, defined in mo_land_memory
  !CALL init_soil_variable_accumulation(tile) follows later
  ! 
  ! 1.1 Initialise time step
  ! 
  time = 0.0_wp
  dt = 1800.0_wp     ! time step length in seconds
  ntimestep = NINT(10.0_wp * one_year * one_day / dt)  ! number of time steps
  latitude = 50.0_wp                       ! latitude in degrees
  longitude = 0.0_wp                       ! longitude in degrees
  annual_t_air_min = -5.0_wp + Tzero       ! annual minumum average day temperature
  annual_t_air_max = 25.0_wp + Tzero       ! annual maximum average day temperature
  diurnal_temp_range = .1_wp             ! diurnal temperature range
  ! 
  ! 1.2 Initalisation of met-forcing to standard values
  ! 
  t_air = ref_t_air + Tzero
  press_srf = standard_press_srf
  q_air = calc_spec_humidity_sat(t_air,press_srf) - ref_vpd / press_srf * eps_vpd
  co2_mixing_ratio = ref_co2
  rain = 0.0_wp
  snow = 0.0_wp
  swvis_srf_down = ref_swdown
  swnir_srf_down = swvis_srf_down * 0.1_wp
  swpar_srf_down = swvis_srf_down
  swpar_fdiffuse = ref_fdiffuse
  wind_air = ref_wind
  wind_10m = wind_air
  lw_srf_down = 0.0_wp
  drag_srf = 1.0_wp
  t_acoef = 0.0_wp
  t_bcoef = 0.0_wp
  q_acoef = 0.0_wp
  q_bcoef = 0.0_wp
  pch = 0.0_wp
  cos_zenith_angle = ref_cos_angle
  ga = .1_wp / ( drag_srf * wind_10m )
  ! 
  ! 1.3 Initialise tile
  ! 
  veg_param%tjmax_opt(:) = 20._wp+0.35_wp*(t_air - Tzero) 
  tile(:)%soil%theta = ref_theta
  tile(:)%plant%lai = 1.0_wp
  DO i = 1,ntile
     tile(i)%plant%pool%leaf%carbon = tile(i)%plant%lai / veg_param%sla(tile(i)%plant%pftid)
  ENDDO 
  tile(:)%plant%pool%leaf%nitrogen = tile(:)%plant%pool%leaf%carbon / & 
        ( 25.0_wp * molar_mass_N / molar_mass_C )
  tile(:)%plant%pool%leaf%phosphorus = tile(:)%plant%pool%leaf%nitrogen / &
        ( 14.0_wp * molar_mass_P / molar_mass_N )
  tile(:)%plant%pool%labile%carbon = 100.0_wp
  tile(:)%plant%pool%labile%nitrogen = 10.0_wp
  tile(:)%plant%pool%labile%phosphorus = 1.0_wp
  ! 
  ! 1.4 Read Model options from file (see src/shared & test_bed/plant/namelist.plant)
  ! 
  NAMELIST /plant_model_options/  &    ! see file namelist & READ() function below
    read_output_interval_pool, &
    read_output_interval_flux, &
    output_start_date, &
    output_end_date, &
    flag_optimal_Nfraction, &
    flag_optimal_Nleaf
  ! read the namelist
  nml_unit = open_nml  ('namelist.plant')
  f_unit = position_nml ('PLANT_MODEL_OPTIONS', nml_unit, status=read_status)
    SELECT CASE (read_status)
      CASE (POSITIONED)
        READ (f_unit, plant_model_options)
    END SELECT
  ! pass values as read from file
  options%model%output_interval_pool   = read_output_interval_pool
  options%model%output_interval_flux   = read_output_interval_flux
  options%model%output_start_date      = output_start_date
  options%model%output_end_date        = output_end_date
  options%model%flag_optimal_Nfraction = flag_optimal_Nfraction
  options%model%flag_optimal_Nleaf     = flag_optimal_Nleaf
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
     ! 2.1 advance timer and set orbital parameters
     ! 
     time = time + dt/one_day 
     CALL update_orbital_parameters(time,declination,sw_toa_down)
     ! 
     ! 2.2 update time-variant forcing
     ! 
     SELECT CASE(options%forcing%option_forcing)  ! see mo_base_configuration
       CASE('constant')
         ! default input values are not changed ?
       CASE('readfile')
         ! implement later
       CASE('generate')
         t_air = generate_inst_temperature(latitude,time,declination,annual_t_air_min, &
         annual_t_air_max,diurnal_temp_range)
         cos_zenith_angle = calc_solar_angle(latitude,longitude,time)
         swvis_srf_down = MAX(0.0,ref_swdown * cos_zenith_angle)
         swnir_srf_down = swvis_srf_down * 0.1_wp
         swpar_srf_down = swvis_srf_down
         swpar_fdiffuse = calc_sw_fdiffuse(swvis_srf_down,sw_toa_down,cos_zenith_angle)
         q_air = ref_rh * calc_spec_humidity_sat(t_air,press_srf)
     END SELECT
     ! 
     ! 3.0 canopy calculation fluxes
     ! 
     CALL update_radiation(tile,swvis_srf_down,swnir_srf_down,swpar_srf_down, &
                              swpar_fdiffuse,cos_zenith_angle,alb_vis_dif,alb_vis_dir, &
                              alb_nir_dif,alb_nir_dir)

     CALL update_canopy_fluxes(tile,t_air,q_air,press_srf,co2_mixing_ratio,ga,veg_param)
     ! 
     ! 4.0 calculate vegetation growth and turnover
     ! 
     CALL update_vegetation(tile,t_air,veg_param)
     ! 
     ! 5.0 Accumulation
     ! 
     CALL update_vegetation_moving_averages(tile)
     ! 
     ! 6.0 Output
     ! 
     ! 6.1 pool: write output in certain intervals, see mo_veg_output
     ! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
     ! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
     CALL write_output_pool(time, tile)
     ! 6.2 flux: write output in certain intervals, see mo_veg_output
     CALL write_output_flux(time, tile)
     ! 
     ! 7.0 check timer
     ! 
     IF(time > one_year) time = time - one_year   ! largest value of time ~ 365; see mo_base_constants

  ENDDO

  CALL close_output_files


END PROGRAM plant_model_driver
