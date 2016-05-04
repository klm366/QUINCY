!>
!! @brief driver for the stand-alone canopy model (fapar and photosynthesis) 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
PROGRAM canopy_model_driver
  USE mo_kind
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
  USE mo_land_constants,        ONLY: molar_mass_N
  ! vegetation
  USE mo_veg_parameters
  USE mo_veg_radiation
  USE mo_veg_photosynthesis
  USE mo_veg_soil
  USE mo_veg_config_class,      ONLY: veg_config
  USE mo_veg_memory,            ONLY: set_zero_fluxes_of_plant_pools
  ! shared
  !USE mo_namelist,              ONLY: POSITIONED, open_nml, position_nml  ! read model options from file

  IMPLICIT NONE

  TYPE(t_slm_veg_param),SAVE               :: veg_param   ! data structure for vegetation parameters
  TYPE(t_slm_tile),DIMENSION(ntile), SAVE  :: tile        ! data structure for dynamic vegetation and soil variables

  INTEGER  :: i
  REAL(wp) :: time_day                     ! displays the 'real time' in days; time_day == 1.0 after one day has passed
  INTEGER  :: time_year                    ! displays the 'real time' in years; the 1st year of a simulation == 1 (not zero)
  INTEGER  :: ntimestep 
  REAL(wp) :: timestep_prop_day            ! proportion of one day that is simulated with one timestep (dt/one_day)
  REAL(wp) :: ref_lai = 3._wp              ! LAI
  REAL(wp) :: ref_nleaf = 2.0_wp           ! N leaf per unit leaf area in gN/m2

  ! 
  ! 1.0 Read config paramters from namelist.canopy
  ! 
  base_config%namelist_filename       = 'namelist.canopy'
  veg_config%namelist_filename        = 'namelist.canopy'
  atmland_config%namelist_filename    = 'namelist.canopy'
  CALL base_config%Init
  CALL veg_config%Init
  CALL atmland_config%Init
  ! 
  ! 1.1 Initialisation of vegetation parameters, constants and forcing
  ! 
  CALL init_base_constants
  CALL init_atmland_constants
  ! vegetation (2do: may combine into one function init_vegetation() )
  CALL init_element_zero                                    ! mo_base_memory
  CALL set_zero_fluxes_of_plant_pools(tile(1)%plant%flux)   ! mo_veg_memory ; called each timestep in update_vegetation()
  CALL init_plant(tile(1)%plant)                            ! mo_base_memory; 
  CALL init_veg_parameters(veg_param)
  ! forcing
  CALL init_forcing    ! needs atmland_config%Init to run before
  ! 
  ! 1.2 Initialise time step
  ! 
  time_day   = 0.0_wp                      ! equals 1.0 after one day has passed (i.e. start with zero)
  time_year  = 1                           ! the 1st year of a simulation == 1 (not zero)
  dt        = base_config%dt               ! time step length in seconds, default: 1800.0 == 30 minutes
  ntimestep = base_config%Timesteps()      ! number of time steps
  timestep_prop_day = dt/one_day                  ! proportion of one day that is simulated with one timestep
  ! 
  ! 1.3 Initialise tile (Initialisation of soil and vegetation, so far relevant) 
  ! 
  veg_param%tjmax_opt(:) = 20._wp+0.35_wp*(atm2land%t_air - Tzero)  ! conversion from K to °C
  tile(:)%soil%moisture   = veg_config%soil_awc_prescribe
  tile(:)%soil%theta      = veg_config%soil_theta_prescribe
  tile(1)%plant%lai = ref_lai 
  IF(veg_config%pft_id==0) THEN
    tile(1)%plant%pool%leaf%nitrogen=ref_nleaf/molar_mass_N*tile(1)%plant%lai
  ELSE
    tile(1)%plant%pool%leaf%nitrogen=tile(1)%plant%lai & 
             / veg_param%sla(veg_config%pft_id)/veg_param%cn_leaf(veg_config%pft_id)
    tile(1)%plant%pool%leaf%phosphorus=tile(1)%plant%pool%leaf%nitrogen & 
             / veg_param%np_leaf(veg_config%pft_id)
    tile(1)%plant%pftid=veg_config%pft_id
  ENDIF
  ! 
  ! 1.4 open input files (and read the header into dummy variables); see mo_atmland_forcing
  ! 
  CALL open_forcing_input_file  ! only if: options%forcing%option_forcing == 'readfile'
  ! 
  ! 2.0 Time loop 
  ! 
  DO i = 1,ntimestep
     ! 
     ! 2.1 advance timer
     ! 
     time_day  = time_day  + timestep_prop_day         ! time in days
     ! 
     ! 2.2 set vegetation fluxes to zero
     ! 
     CALL set_zero_fluxes_of_plant_pools(tile(1)%plant%flux)   ! mo_veg_memory
     ! 
     ! 2.3 set orbital parameters and update time-variant forcing; 
     ! 
     CALL update_forcing(time_day)  ! see mo_atmland_forcing
     ! 
     ! 3.0 canopy calculation fluxes
     ! 
     CALL update_radiation(tile)                 ! atm2land land2atm are global TYPES, see mo_atmland_forcing
     CALL update_canopy_fluxes(tile(1), veg_param)  ! atm2land is a global TYPE, see mo_atmland_forcing
     CALL update_vegsoil_physics(tile)
     ! 
     ! 4.0 write output  -> filename: fort.20
     ! 
     WRITE(20,'(6F16.7)') time_day+time_year*one_year,atm2land%t_air, atm2land%swpar_srf_down * sw2par, & 
           tile(1)%plant%flux%gpp,tile(1)%plant%gc,tile(1)%plant%flux%ci
     ! 
     ! 5.0 check timer; adjust year counter
     ! 
     IF(time_day > one_year) THEN
       time_day  = time_day - one_year   ! largest value of time ~ 365; see mo_base_constants
       time_year = time_year + 1
     END IF
  ENDDO

END PROGRAM canopy_model_driver
