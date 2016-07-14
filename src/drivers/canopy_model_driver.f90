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
  ! shared
  USE mo_kind
  USE mo_exception,             ONLY: finish ! message_text
  ! base
  USE mo_base_memory,           ONLY: t_slm_tile
  USE mo_base_constants
  USE mo_base_config_class,     ONLY: base_config  !< TYPE, EXTENDS(t_jsb_config) :: t_base_config
  ! atmland
  USE mo_atmland_constants
  USE mo_atmland_memory
  USE mo_atmland_util
  USE mo_atmland_forcing
  USE mo_atmland_config_class,  ONLY: atmland_config
  ! shared
  USE mo_shared_constants,      ONLY: molar_mass_N
  USE mo_shared_util,           ONLY: init_element_zero
  ! vegetation
  USE mo_veg_parameters
  USE mo_veg_radiation
  USE mo_veg_photosynthesis
  USE mo_veg_soil
  USE mo_veg_config_class,      ONLY: veg_config
  USE mo_veg_memory
  USE mo_veg_util

  IMPLICIT NONE

  TYPE(t_slm_veg_param),SAVE               :: veg_param   ! data structure for vegetation parameters
  TYPE(t_slm_tile),DIMENSION(ntile), SAVE  :: tile        ! data structure for dynamic vegetation and soil variables

  REAL(wp) :: direction    ! this is a temporary solution to make the canopy model compile

  REAL(wp) :: time_day                     !< displays the 'real time' in days; time_day == 1.0 after one day has passed
  INTEGER  :: time_year                    !< displays the 'real time' in years; the 1st year of a simulation == 1 (not zero)
  INTEGER  :: ntimestep                    !< number of timesteps the model runs
  INTEGER  :: itimestep                    !< current timestep i, looping over timesteps ntimestep
  REAL(wp) :: timestep_prop_day            !< proportion of one day that is simulated with one timestep (dt/one_day)
  INTEGER  :: itile     = 1                !< the current tile(i), for looping over ntile
  REAL(wp) :: ref_lai   = 3._wp            !< LAI
  REAL(wp) :: ref_nleaf = 2.0_wp           !< N leaf per unit leaf area in gN/m2

  direction = 1.0_wp         ! this is a temporary solution to make the canopy model compile 

  ! check number of tiles, exit if more than one tile was defined (canopy output routines can handle just one; for now)
  IF(ntile > 1) CALL finish('canopy_model_driver','Number of tiles is too large for output routines to work reliable')

  !! 1.0 Read config paramters from namelist.canopy section
  !!
  !!  
  base_config%namelist_filename       = 'namelist.canopy'
  veg_config%namelist_filename        = 'namelist.canopy'
  atmland_config%namelist_filename    = 'namelist.canopy'
  CALL base_config%Init
  CALL veg_config%Init
  CALL atmland_config%Init

  !! 1.1 Initialisation of vegetation parameters, constants and forcing section
  !! 
  !! 
  CALL init_base_constants
  CALL init_atmland_constants
  CALL init_element_zero                               ! mo_shared_util; run before init_plant() & set_zero_fluxes_of_plant_pools()
  ! vegetation
  CALL tile(itile)%plant%init                       ! mo_veg_memory  ==> may move to veg init (or model init ?)
  CALL set_zero_fluxes_of_plant_pools(tile(itile))     ! mo_veg_util; called each timestep in update_vegetation()
  CALL init_veg_parameters(veg_param)                  ! mo_veg_parameters
  ! forcing
  CALL init_forcing                                             ! needs atmland_config%Init to run before

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
  !! Initialisation of soil and vegetation, so far relevant
  veg_param%tjmax_opt(:)      = 20._wp+0.35_wp*(atm2land%t_air - Tzero)  ! conversion from K to °C
  tile(itile)%soil%moisture   = veg_config%soil_awc_prescribe
  tile(itile)%soil%theta      = veg_config%soil_theta_prescribe
  tile(itile)%plant%lai       = ref_lai 
  IF(veg_config%pft_id==0) THEN
    tile(itile)%plant%pool%leaf%nitrogen   = ref_nleaf / molar_mass_N * tile(itile)%plant%lai
  ELSE
    tile(itile)%plant%pool%leaf%nitrogen   = tile(itile)%plant%lai / & 
                                              veg_param%sla(veg_config%pft_id) / veg_param%cn_leaf(veg_config%pft_id)
    tile(itile)%plant%pool%leaf%phosphorus = tile(itile)%plant%pool%leaf%nitrogen / & 
                                              veg_param%np_leaf(veg_config%pft_id)
    tile(itile)%plant%pftid                = veg_config%pft_id
  ENDIF

  !! 1.4 open input files (and read the header into dummy variables); see mo_atmland_forcing
  !! 
  !! 
  CALL open_forcing_input_file  ! only if: options%forcing%option_forcing == 'readfile'

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
    CALL update_forcing(time_day)  ! see mo_atmland_forcing

    !! NO tile loop - compared to the plant model driver !!
    !! the below section numbers equal the plant_model_driver section numbers !!

      !! 3.1 set vegetation fluxes to zero
      !! 
      !! 
      CALL set_zero_fluxes_of_plant_pools(tile(itile))   !< set to zero all fluxes of the plant pools; mo_veg_util

      !! 3.2 canopy calculation fluxes
      !! 
      !! 
      CALL update_radiation(tile(itile))                                 !< mo_veg_radiation
    !> update_canopy_fluxes needs the direction variable to work
    !! @todo: ticket #52 , proper implementation of "direction" consistent with plant model
      CALL update_canopy_fluxes(tile(itile), veg_param)       !< mo_veg_photosynthesis

      !! 3.5 update soil water bucket 
      !! 
      !! 
      CALL update_vegsoil_physics(tile(itile))                !< mo_veg_soil

      !! 3.8 output fluxes -> filename: fort.20
      !! 
      !! output also some atm2land variables
      WRITE(20,'(6F16.7)') time_day+time_year*one_year, atm2land%t_air, atm2land%swpar_srf_down * sw2par, & 
                           tile(itile)%plant%flux%gpp, tile(itile)%plant%gc, tile(itile)%plant%flux%ci

    !! END tile loop in the plant model driver !!

    !! 2.3 check timer; adjust year counter
    !! 
    !! 
    IF(time_day > one_year) THEN
      time_day  = time_day - one_year   ! largest value of time ~ 365; see mo_base_constants
      time_year = time_year + 1
    END IF
  END DO ! end time loop

END PROGRAM canopy_model_driver
