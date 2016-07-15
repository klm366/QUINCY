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
  USE mo_exception,               ONLY: finish, message, message_text
  ! land model
  USE mo_land_model_interface,    ONLY: land_model_interface 
  ! base
  USE mo_base_constants
  USE mo_base_memory,             ONLY: t_slm_tile, t_model
  USE mo_base_config_class,       ONLY: base_config  !< TYPE, EXTENDS(t_jsb_config) :: t_base_config
  ! shared
  USE mo_shared_isotope,          ONLY: calc_fractionation 
  ! atmland
  USE mo_atmland_constants
  USE mo_atmland_memory
  USE mo_atmland_forcing,         ONLY: init_forcing, update_forcing, open_forcing_input_file, close_forcing_input_file
  USE mo_atmland_config_class,    ONLY: atmland_config
  ! vegetation
  USE mo_veg_parameters,          ONLY: t_slm_veg_param, init_veg_parameters, &
                                        fresp_growth, ic3phot, CiCa_default_C3, CiCa_default_C4
  USE mo_veg_config_class,        ONLY: veg_config
  USE mo_veg_photosynthesis,      ONLY: discrimination_ps 
  USE mo_veg_output,              ONLY: init_veg_output_variables, open_pool_output_files, open_flux_output_files, &
                                        close_output_files

  IMPLICIT NONE

  TYPE(t_model), SAVE                      :: model            !< very basic version of one JSBACH4 model
  TYPE(t_slm_veg_param),SAVE               :: veg_param        !< data structure for vegetation parameters
  TYPE(t_slm_tile),DIMENSION(ntile), SAVE  :: tile             !< data structure for dynamic vegetation and soil variables

  ! Local variables
  REAL(wp) :: time_day                     !< displays the 'real time' in days; time_day == 1.0 after one day has passed
  INTEGER  :: time_year                    !< displays the 'real time' in years; the 1st year of a simulation == 1 (not zero)
  INTEGER  :: ntimestep                    !< number of timesteps the model runs
  INTEGER  :: itimestep                    !< current timestep i, looping over timesteps ntimestep
  REAL(wp) :: timestep_prop_day            !< proportion of one day that is simulated with one timestep (dt/one_day)
  INTEGER  :: itile                        !< the current tile(i), for looping over ntile
  INTEGER  :: icanop                       !< the current canopy layer, for looping over ncanop
  REAL(wp) :: hlp1
  REAl(wp), DIMENSION(ncanop) :: canopy_layer_depth_act !< actual canopy layer depth ?


  ! check number of tiles, exit if more than output routines can handle
  IF(ntile > 1) CALL finish('land_model_driver','Number of tiles is too large, i.e. larger than 1')
  ! 

  !! 1.0 Read config paramters from namelist.plant section
  !!
  !!  
  base_config%namelist_filename       = 'namelist.plant'
  veg_config%namelist_filename        = 'namelist.plant'
  atmland_config%namelist_filename    = 'namelist.plant'
  CALL base_config%Init
  CALL veg_config%Init
  CALL atmland_config%Init

  !! 1.1 Model init   ===> move to extra file
  !!
  !!  
  !! see email of Reiner from 20.6.16 \n
  !! in jsbach4 this is part of the model_init routine
  !! 
  model%id = 0                                    ! this is just an arbitrary value consistent with tile%owner_model_id
  model%shortname = "slm"                         ! arbitrary name

  !! 1.2 Initialise tile & Parameter Init & Time Init
  !! 
  !! 
  !! @todo land model driver: update Intel compiler to enable: init lcts of the tile via PROCEDURE Construct_jsb_lct
  ! lcts
  tile(1)%lcts%id        = 2                  ! VEG_TYPE     !< Vegetated
  tile(1)%lcts%name      = "veg"              ! Short name for lct
  tile(1)%lcts%long_name = "slm vegetation"   ! Long name for lct
  tile(1)%lcts%lib_id    = veg_config%pft_id  ! pftID == column number in model%slm_parameter (SLM) == Id of lct in lctlib (JSBACH4)

  ! vegetation parameters, constants and forcing section
  CALL init_base_constants
  CALL init_atmland_constants

  ! init time
  time_day          = 0.0_wp                      ! equals 1.0 after one day has passed (i.e. start with zero)
  time_year         = 1                           ! the 1st year of a simulation == 1 (not zero)
  dt                = base_config%dt              ! time step length in seconds, default: 1800.0 == 30 minutes
  ntimestep         = base_config%Timesteps()     ! number of time steps
  timestep_prop_day = dt/one_day                  ! proportion of one day that is simulated with one timestep

  ! vegetation
  DO itile = 1,ntile
    CALL tile(itile)%plant%init                       ! mo_veg_memory  ==> may move to veg init (or model init ?)
    CALL tile(itile)%plant%flux%set_zero              ! mo_veg_memory  ==> may move to veg init (or model init ?)
  END DO
  CALL init_veg_parameters(veg_param)                  ! mo_veg_parameters
  ! forcing
  CALL init_forcing                                !< needs atmland_config%Init to run before
  ! accumulation
  DO itile = 1,ntile
    CALL tile(itile)%plant%daily_average%init      !< initially sets all variables to zero or .FALSE.
    CALL tile(itile)%plant%moving_average%init     !< init variables with reasonable default values and sets accumulation periods
  END DO

  ! lctlib - replacing init_veg_parameters(veg_param)
!  CALL model%slm_parameter%init   ! should be deleted in a later version prior to slm-jsbach4 linking
!  model%lctlib => slm_lctlib_wrapper(model%slm_parameter)     ! mo_base_slm_parameter_class

  veg_param%tjmax_opt(:)  = 20._wp+0.35_wp*(atm2land%t_air - Tzero)   ! conversion from K to °C
  DO itile = 1,ntile
    tile(itile)%soil%moisture                 = veg_config%soil_awc_prescribe
    tile(itile)%soil%theta                    = veg_config%soil_theta_prescribe 
    tile(itile)%soil%ammonium                 = 0.1_wp
    tile(itile)%soil%phosphate                = 0.01_wp
    tile(itile)%soil%pool%litter%nitrogen     = 5.0_wp
    tile(itile)%plant%pftid                   = veg_config%pft_id
    tile(itile)%plant%target_cn_leaf          = veg_param%cn_leaf(tile(itile)%plant%pftid)
    tile(itile)%plant%target_np_leaf          = veg_param%np_leaf(tile(itile)%plant%pftid)
    tile(itile)%plant%pool%reserve%carbon     = 50.0_wp
    tile(itile)%plant%pool%reserve%nitrogen   = tile(itile)%plant%pool%reserve%carbon / &
                                                (1._wp + fresp_growth) / veg_param%cn_leaf(tile(itile)%plant%pftid)
    tile(itile)%plant%pool%reserve%phosphorus = tile(itile)%plant%pool%reserve%nitrogen / &
                                                veg_param%np_leaf(tile(itile)%plant%pftid)
    IF(veg_param%ps_pathway(tile(itile)%plant%pftid) == ic3phot)THEN
       hlp1 = discrimination_ps(CiCa_default_C3*atm2land%co2_mixing_ratio,atm2land%co2_mixing_ratio, & 
                                'C13',veg_param%ps_pathway(tile(itile)%plant%pftid))
    ELSE
       hlp1 = discrimination_ps(CiCa_default_C4*atm2land%co2_mixing_ratio,atm2land%co2_mixing_ratio, & 
                                'C13',veg_param%ps_pathway(tile(itile)%plant%pftid))
    ENDIF
    tile(itile)%plant%pool%reserve%carbon13   = calc_fractionation(atm2land%co2_mixing_ratio, & 
                                                                   atm2land%co2_mixing_ratio_C13,hlp1) & 
                                                * tile(itile)%plant%pool%reserve%carbon 
    IF(veg_param%ps_pathway(tile(itile)%plant%pftid) == ic3phot)THEN
       hlp1 = discrimination_ps(CiCa_default_C3*atm2land%co2_mixing_ratio,atm2land%co2_mixing_ratio, & 
                                'C14',veg_param%ps_pathway(tile(itile)%plant%pftid))
    ELSE
       hlp1 = discrimination_ps(CiCa_default_C4*atm2land%co2_mixing_ratio,atm2land%co2_mixing_ratio, & 
                                'C14',veg_param%ps_pathway(tile(itile)%plant%pftid))
    ENDIF
    tile(itile)%plant%pool%reserve%carbon14   = calc_fractionation(atm2land%co2_mixing_ratio, & 
                                                                   atm2land%co2_mixing_ratio_C14,hlp1) & 
                                                * tile(itile)%plant%pool%reserve%carbon 

  END DO ! end loop over tiles

  !! 1.3 open input files (and read the header into dummy variables); see mo_atmland_forcing
  !! 
  !! 
  CALL open_forcing_input_file  ! mo_atmland_forcing; only if: option_forcing == 'readfile'

  !! 1.4 Initialise output and flux accumulation, and open output files; see mo_veg_output
  !! 
  !! 
  CALL init_veg_output_variables     ! mo_veg_output : interval variables & flux accumulation per tile
  CALL open_pool_output_files(1)     ! 1 == itile 
  CALL open_flux_output_files(1)     ! 1 == itile 


  !! 2.0 time loop
  !! 
  DO itimestep = 1, ntimestep

    !! 2.1 advance timer
    !! 
    !! 
    time_day  = time_day  + timestep_prop_day         ! time in days

    !! 2.2 set orbital parameters and update time-variant forcing
    !! 
    !! 
    CALL update_forcing(time_day)               !< mo_atmland_forcing

    !! 2.3 land model interface
    !! 
    !! 
    !! ==> implement later: diff between models by: land_model_interface(tile, model, ....)
    !! model would be one of the strings land plant canopy soil canopy_test -> use enum to attach numbers to strings !!
    !! ==> 
    CALL land_model_interface(tile(1), veg_param, &
                              atm2land%t_air, atm2land%press_srf, atm2land%q_air, &
                              time_day, time_year, timestep_prop_day &
                             )

    !! 2.4 check timer & adjust year counter
    !! 
    !! 
    IF(time_day > one_year) THEN
      time_day  = time_day - one_year   ! largest value of time ~ 365; see mo_base_constants
      time_year = time_year + 1
    END IF

  END DO  ! end time loop 

  !! 3.0 Close output files
  !! 
  CALL close_output_files(1)     ! 1 == itile 
  CALL close_forcing_input_file  ! only if: option_forcing == 'readfile'

END PROGRAM land_model_driver
