!>
!! @brief Contains structures and methods for basic model config, temporary solution until the JSBACH4 config style is
!!        adopted in full
!!
!! @author
!!  Jan Engel, MPI-BGC
!!
!! @par Revision History
!! First version by Jan Engel (2016-04-11)
!! based on the mo_veg_config_class, which was initially developed by Reiner Schnur MPI-M Hamburg
!!
MODULE mo_base_config_class

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: base_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the model base, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_base_config
    !! BASE_CTL
    REAL(wp)          :: dtime                             !< length of one time step in seconds; base_ctl namelist: dtime_step_length_sec
    CHARACTER(len=20) :: slm_model_name
    CHARACTER(8)      :: output_interval_pool, &
                         output_interval_flux
    INTEGER           :: output_start_first_day_year, &  !< output start at the 1st  day of this year
                         output_end_last_day_year        !< output stop  at the last day of this year
    LOGICAL           :: fluxnet_output
    INTEGER           :: fluxnet_output_length
    LOGICAL           :: flag_generate_soil_forcing       !< all litterfall pools + some soilphysics 
                                                          !! at period: base_config%output_start_first_day_year .. output_end_last_day_year
    LOGICAL           :: set_parameter_values_from_file
    LOGICAL           :: flag_climate_forcing_one_year
    LOGICAL           :: flag_stand_harvest
    INTEGER           :: stand_replacing_year
    !! GRID_CTL
    REAL(wp)          :: lat, lon                           !< latitude longitude currently stored in mo_base_memory:t_slm_grid
    
    !! other parameters than namelist parameters
    LOGICAL           :: l_do_stand_replacing_harvest
    LOGICAL           :: l_transient_forcing_mode, &        !< set to TRUE if forcing_options%forcing_mode == transient
                         l_transient_spinup, &              !< TRUE for forcing_mode=='transient' and model spinup, otherwise FALSE
                         l_transient_spinup_rewind_forcing  !< TRUE at the timestep after one transient spinup period has finished (e.g. after 1930 for 1901-1930)
  CONTAINS
    PROCEDURE :: Init => Init_base_config
  END Type t_base_config

  TYPE(t_base_config), SAVE :: base_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_base_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_base_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_base_config(config)  ! this is a temporary solution until the JSBACH4 config style is adopted in full

    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml
    ! 0.1 InOut
    CLASS(t_base_config), INTENT(inout) :: config    !< config type for base
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    ! BASE_CTL
    REAL(wp)          :: dtime_step_length_sec   ! model-variable: dtime 
    CHARACTER(len=20) :: slm_model_name
    CHARACTER(8)      :: output_interval_pool          , &
                         output_interval_flux
    INTEGER           :: output_start_first_day_year   , &
                         output_end_last_day_year
    LOGICAL           :: fluxnet_output
    INTEGER           :: fluxnet_output_length
    LOGICAL           :: flag_generate_soil_forcing                     !< all litterfall pools + some soilphysics
    LOGICAL           :: set_parameter_values_from_file
    LOGICAL           :: flag_climate_forcing_one_year
    LOGICAL           :: flag_stand_harvest
    INTEGER           :: stand_replacing_year
    ! GRID_CTL
    REAL(wp)          :: latitude
    REAL(wp)          :: longitude
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_base_config'

    NAMELIST /base_ctl/                         &  
        dtime_step_length_sec                 , &
        slm_model_name                        , &
        output_interval_pool                  , &
        output_interval_flux                  , &
        output_start_first_day_year           , &
        output_end_last_day_year              , &
        fluxnet_output                        , &
        fluxnet_output_length                 , &
        flag_generate_soil_forcing            , &
        set_parameter_values_from_file        , &
        flag_climate_forcing_one_year         , &
        flag_stand_harvest                    , &
        stand_replacing_year

    NAMELIST /grid_ctl/                         &  
        latitude                              , &
        longitude


    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting base configuration')

    ! Set defaults
    ! 
    ! BASE_CTL
    dtime_step_length_sec                   = 1800.0_wp
    slm_model_name                          = 'land'
    output_interval_pool                    = 'weekly' ! timestep daily weekly monthly yearly
    output_interval_flux                    = 'weekly' ! timestep daily weekly monthly yearly
    output_start_first_day_year             = 1
    output_end_last_day_year                = 10
    fluxnet_output                          = .FALSE.     ! .TRUE. .FALSE.
    fluxnet_output_length                   = 2
    flag_generate_soil_forcing              = .FALSE.     ! .TRUE. .FALSE.
    set_parameter_values_from_file          = .FALSE.     ! .TRUE. .FALSE.
    flag_climate_forcing_one_year           = .FALSE.     ! .TRUE. .FALSE.
    flag_stand_harvest                      = .FALSE.     ! .TRUE. .FALSE.
    stand_replacing_year                    = 1500        ! the year 1500
    ! GRID_CTL
    latitude                                = 50.91       ! about the MPI-BGC Jena
    longitude                               = 11.56       ! about the MPI-BGC Jena


    ! read the namelist
    ! 
    ! &BASE_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('BASE_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, base_ctl)
    CALL close_nml(nml_handler)
    ! &GRID_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('GRID_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, grid_ctl)
    CALL close_nml(nml_handler)

    ! pass values as read from file
    ! 
    ! BASE_CTL
    config%dtime                                  = dtime_step_length_sec
    config%slm_model_name                         = slm_model_name
    config%output_interval_pool                   = output_interval_pool
    config%output_interval_flux                   = output_interval_flux
    config%output_start_first_day_year            = output_start_first_day_year
    config%output_end_last_day_year               = output_end_last_day_year
    config%fluxnet_output                         = fluxnet_output
    config%fluxnet_output_length                  = fluxnet_output_length
    config%flag_generate_soil_forcing             = flag_generate_soil_forcing
    config%set_parameter_values_from_file         = set_parameter_values_from_file
    config%flag_climate_forcing_one_year          = flag_climate_forcing_one_year
    config%flag_stand_harvest                     = flag_stand_harvest              ! might be modified below IF (config%stand_replacing_year .LT. 1901)
    config%stand_replacing_year                   = stand_replacing_year
    
    ! GRID_CTL
    config%lat                                    = latitude
    config%lon                                    = longitude
    
    
    ! limit the stand replacement functionality to years after 1900
    ! i.e., if there was stand replacement before 1901 this may not be considered by SLM (see site_list.dat column 21: PlantYear)
    IF (config%stand_replacing_year .LT. 1901) THEN
      config%flag_stand_harvest = .FALSE.
    END IF
    
    ! set default for other config% parameter that are changed during runtime not via namelist
    config%l_do_stand_replacing_harvest           = .FALSE.
    
  END SUBROUTINE Init_base_config


!#endif
END MODULE mo_base_config_class
