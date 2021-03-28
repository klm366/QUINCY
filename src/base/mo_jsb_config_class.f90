!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
!! @brief Contains structures and methods for JSBACH model config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! integration into QUINCY by Jan Engel MPI-BGC (2016-04-01) & update 2017-05-16
!!
MODULE mo_jsb_config_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp                              ! *slmjsb4
  USE mo_io_units,          ONLY: filename_max
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_config, t_jsb_config_p, t_jsb_model_config, new_model_config


  !-----------------------------------------------------------------------------------------------------
  !> model configuration TYPE
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_jsb_model_config
    CHARACTER(len=filename_max)    :: grid_filename
    CHARACTER(len=filename_max)    :: fract_filename
    CHARACTER(len=30)              :: grid_name
    CHARACTER(len=40)              :: usecase
    LOGICAL                        :: use_lakes
    LOGICAL                        :: use_glacier
    CHARACTER(len=10)              :: tpe_scheme       ! For terraplanet setup: open/closed
    CHARACTER(len=10)              :: hsm_mode
    LOGICAL                        :: l_compat401      ! Use configuration compatible with jsbach 4.01
    CHARACTER(len=10), ALLOCATABLE :: output_tiles(:)  ! List of tile names for which output should be generated

    !! QUINCY added !!
    !! BASE_CTL
    REAL(wp)          :: dtime                                    !< length of one time step in seconds; base_ctl namelist: dtime_step_length_sec
    CHARACTER(len=20) :: slm_model_name
    CHARACTER(8)      :: output_interval_pool, &                  !< interval at which output is generated during model runs not including spinup
                         output_interval_flux                     !< interval at which output is generated during model runs not including spinup - the pool value is used
    CHARACTER(8)      :: output_interval_pool_spinup, &           !< interval at which output is generated during spinup
                         output_interval_flux_spinup              !< interval at which output is generated during spinup - the pool value is used
    INTEGER           :: output_start_first_day_year, &           !< output start at the 1st  day of this year
                         output_end_last_day_year                 !< output stop  at the last day of this year
    LOGICAL           :: fluxnet_type_transient_timestep_output   !< transient simulation with fluxnet-type site-set and timestep-output for the last full forcing period (static forcing)
    INTEGER           :: fluxnet_static_forc_start_yr             !< fluxnet sites: first year of (static) forcing available
    INTEGER           :: fluxnet_static_forc_last_yr              !< fluxnet sites: last year with (static) forcing available
    CHARACTER(len=20) :: sim_siteset_specific                     !< specific simulation setup for particular site-set, e.g., EucFACE transient simulations
    LOGICAL           :: flag_generate_soil_forcing_stat          !< usecase 01: all litterfall pools + some soilphysics | static simulation
                                                                  !! at period: model%config%output_start_first_day_year .. output_end_last_day_year
    LOGICAL           :: flag_generate_soil_forcing_trans         !< usecase 08: all litterfall pools + some soilphysics | transient simaulation
    LOGICAL           :: set_parameter_values_from_file           !< set values of (some) parameters to values given by namelist (parameter_slm_run.list)
    LOGICAL           :: set_parameter_values_proportional        !< if TRUE: change the default (hardcoded) parameter values by proportion as defined in namelist
    LOGICAL           :: flag_climate_forcing_one_year
    LOGICAL           :: flag_stand_harvest
    INTEGER           :: stand_replacing_year
    LOGICAL           :: flag_spp1685
    ! GRID_CTL
    REAL(wp)          :: lat
    REAL(wp)          :: lon
    !! other parameters than namelist parameters => not set via namelist but at model runtime
    LOGICAL           :: l_do_stand_replacing_harvest
    LOGICAL           :: l_transient_forcing_mode, &        !< set to TRUE if forcing_options%forcing_mode == transient
                         l_transient_spinup, &              !< TRUE for forcing_mode=='transient' and model spinup, otherwise FALSE
                         l_transient_spinup_rewind_forcing  !< TRUE at the timestep after one transient spinup period has finished (e.g. after 1930 for 1901-1930)
    LOGICAL           :: l_static_forcing_mode              !< set to TRUE if forcing_options%forcing_mode == static
  END TYPE t_jsb_model_config


  !-----------------------------------------------------------------------------------------------------
  !> Abstract TYPE for process and model configuration
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, ABSTRACT :: t_jsb_config
    LOGICAL :: active
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename
    CHARACTER(LEN=filename_max) :: namelist_filename
    CHARACTER(len=10), ALLOCATABLE :: output_tiles(:)  ! List of tile names for which output should be generated
  CONTAINS
    PROCEDURE (Init_config), DEFERRED, PASS(config) :: Init
  END TYPE t_jsb_config

  TYPE t_jsb_config_p
    CLASS(t_jsb_config), POINTER :: p
  END TYPE t_jsb_config_p

  ! note from JanEngel *slmtemp
  ! https://de.wikibooks.org/wiki/Fortran:_Fortran_2003:_Zeiger#Abstraktes_Interface
  ! defines that the TYPE itself is passed to the procedure (i.e. t_*_config -> init_config subroutine) 
  ABSTRACT INTERFACE
    SUBROUTINE Init_config(config)
      IMPORT :: t_jsb_config    ! makes t_jsb_config accessible within this (abstract) interface
      CLASS(t_jsb_config), INTENT(inout) :: config
    END SUBROUTINE Init_config
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_config_class'


CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> create and init model config
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION new_model_config(namelist_filename) RESULT(model_config)

    !USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_namelist, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CHARACTER(len=*),         INTENT(in) :: namelist_filename
    TYPE(t_jsb_model_config), POINTER    :: model_config

    ! GRID_CTL
    REAL(wp)                    :: latitude
    REAL(wp)                    :: longitude
    CHARACTER(len=filename_max) :: grid_filename
    CHARACTER(len=30)           :: grid_name

    NAMELIST /grid_ctl/                         &  
        latitude                              , &
        longitude
    ! NAMELIST /jsb_grid_nml/   &
    !   grid_filename,          &
    !   grid_name

    ! BASE_CTL
    REAL(wp)                    :: dtime_step_length_sec                      !< model-variable: dtime 
    CHARACTER(len=20)           :: slm_model_name
    CHARACTER(8)                :: output_interval_pool          , &
                                   output_interval_flux
    CHARACTER(8)                :: output_interval_pool_spinup   , & 
                                   output_interval_flux_spinup       
    INTEGER                     :: output_start_first_day_year   , &
                                   output_end_last_day_year
    LOGICAL                     :: fluxnet_type_transient_timestep_output
    INTEGER                     :: fluxnet_static_forc_start_yr
    INTEGER                     :: fluxnet_static_forc_last_yr
    CHARACTER(len=20)           :: sim_siteset_specific
    LOGICAL                     :: flag_generate_soil_forcing_stat
    LOGICAL                     :: flag_generate_soil_forcing_trans
    LOGICAL                     :: set_parameter_values_from_file
    LOGICAL                     :: set_parameter_values_proportional
    LOGICAL                     :: flag_climate_forcing_one_year
    LOGICAL                     :: flag_stand_harvest
    INTEGER                     :: stand_replacing_year
    LOGICAL                     :: flag_spp1685

    LOGICAL                     :: use_lakes
    LOGICAL                     :: use_glacier
    LOGICAL                     :: l_compat401
    CHARACTER(len=filename_max) :: fract_filename
    CHARACTER(len=40)           :: usecase
    CHARACTER(len=10)           :: tpe_scheme
    CHARACTER(len=10)           :: hsm_mode
    CHARACTER(len=10)           :: output_tiles(99)

    NAMELIST /base_ctl/                           &  
        dtime_step_length_sec                   , &
        slm_model_name                          , &
        output_interval_pool                    , &
        output_interval_flux                    , &
        output_interval_pool_spinup             , &
        output_interval_flux_spinup             , &
        output_start_first_day_year             , &
        output_end_last_day_year                , &
        fluxnet_type_transient_timestep_output  , &
        fluxnet_static_forc_start_yr            , &
        fluxnet_static_forc_last_yr             , &
        sim_siteset_specific                    , &
        flag_generate_soil_forcing_stat         , &
        flag_generate_soil_forcing_trans        , &
        set_parameter_values_from_file          , &
        set_parameter_values_proportional       , &
        flag_climate_forcing_one_year           , &
        flag_stand_harvest                      , &
        stand_replacing_year                    , &
        flag_spp1685                          
                          
    ! NAMELIST /jsb_model_nml/  &
    !   use_lakes, &
    !   use_glacier, &
    !   l_compat401, &
    !   usecase, &
    !   tpe_scheme, &
    !   hsm_mode, &
    !   fract_filename, &
    !   output_tiles

    INTEGER :: nml_handler, nml_unit, istat, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':new_model_config'

    CALL message(TRIM(routine), 'Starting model configuration from '//TRIM(namelist_filename))

    ALLOCATE(model_config)

    ! Set defaults
    ! BASE_CTL
    dtime_step_length_sec                   = 1800.0_wp   ! length of one time step in seconds
    slm_model_name                          = 'land'      ! default is the coupled veg-soil model, a.k.a. land
    output_interval_pool                    = 'weekly'    ! timestep daily weekly monthly yearly
    output_interval_flux                    = 'weekly'    ! timestep daily weekly monthly yearly
    output_interval_pool_spinup             = 'yearly'    ! usually ~500 yrs are covered, so 'yearly' would be reasonable
    output_interval_flux_spinup             = 'yearly'    ! usually ~500 yrs are covered, so 'yearly' would be reasonable
    output_start_first_day_year             = 1
    output_end_last_day_year                = 10
    fluxnet_type_transient_timestep_output  = .FALSE.     ! .TRUE. .FALSE.
    fluxnet_static_forc_start_yr            = 1500        ! usually a year after 1990
    fluxnet_static_forc_last_yr             = 1500        ! usually a year after 1990
    sim_siteset_specific                    = 'none'      ! default is 'none'
    flag_generate_soil_forcing_stat         = .FALSE.     ! .TRUE. .FALSE.  << is set to T for usecase 01 only
    flag_generate_soil_forcing_trans        = .FALSE.     ! .TRUE. .FALSE.  << is set to T for usecase 08 only
    set_parameter_values_from_file          = .FALSE.     ! .TRUE. .FALSE.
    set_parameter_values_proportional       = .FALSE.     ! .TRUE. .FALSE.
    flag_climate_forcing_one_year           = .FALSE.     ! .TRUE. .FALSE.
    flag_stand_harvest                      = .FALSE.     ! .TRUE. .FALSE.
    stand_replacing_year                    = 1500        ! the year 1500
    ! GRID_CTL
    latitude                                = 50.91       ! about the MPI-BGC Jena
    longitude                               = 11.56       ! about the MPI-BGC Jena
    ! JSBACH4
    use_lakes                               = .FALSE.
    use_glacier                             = .FALSE.
    l_compat401                             = .FALSE.
    grid_filename                           = ''
    grid_name                               = ''
    usecase                                 = ''
    tpe_scheme                              = ''
    hsm_mode                                = 'simple'
    fract_filename                          = 'bc_land_frac.nc'
    output_tiles(:)                         = ''
    flag_spp1685                            = .FALSE.       ! if TRUE: soil layer specific texture information for SPP_1685 sites: simple_1d and jsm

    ! Read namelist
    nml_handler = open_nml(TRIM(namelist_filename))

    nml_unit = position_nml('GRID_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, grid_ctl)

    ! Write resulting values in config memory
    model_config%lat           = latitude
    model_config%lon           = longitude
    ! model_config%grid_filename = grid_filename
    ! model_config%grid_name     = grid_name

    nml_unit = position_nml('BASE_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, base_ctl)
    CALL close_nml(nml_handler)

    model_config%dtime                                  = dtime_step_length_sec
    model_config%slm_model_name                         = slm_model_name
    model_config%output_interval_pool                   = output_interval_pool
    model_config%output_interval_flux                   = output_interval_flux
    model_config%output_interval_pool_spinup            = output_interval_pool_spinup
    model_config%output_interval_flux_spinup            = output_interval_flux_spinup
    model_config%output_start_first_day_year            = output_start_first_day_year
    model_config%output_end_last_day_year               = output_end_last_day_year
    model_config%fluxnet_type_transient_timestep_output = fluxnet_type_transient_timestep_output
    model_config%fluxnet_static_forc_start_yr           = fluxnet_static_forc_start_yr
    model_config%fluxnet_static_forc_last_yr            = fluxnet_static_forc_last_yr
    model_config%sim_siteset_specific                   = sim_siteset_specific
    model_config%flag_generate_soil_forcing_stat        = flag_generate_soil_forcing_stat
    model_config%flag_generate_soil_forcing_trans       = flag_generate_soil_forcing_trans
    model_config%set_parameter_values_from_file         = set_parameter_values_from_file
    model_config%set_parameter_values_proportional      = set_parameter_values_proportional
    model_config%flag_climate_forcing_one_year          = flag_climate_forcing_one_year
    model_config%flag_stand_harvest                     = flag_stand_harvest              ! might be modified below IF (config%stand_replacing_year .LT. 1901)
    model_config%stand_replacing_year                   = stand_replacing_year
    model_config%flag_spp1685                           = flag_spp1685
    ! model_config%use_lakes      = use_lakes
    ! model_config%use_glacier    = use_glacier
    ! model_config%l_compat401    = l_compat401
    ! model_config%usecase        = usecase
    ! model_config%tpe_scheme     = tpe_scheme
    ! model_config%hsm_mode       = hsm_mode
    ! model_config%fract_filename = fract_filename

    ! limit the stand replacement functionality to years after 1900
    ! i.e., if there was stand replacement before 1901 this may not be considered by SLM (see site_list.dat column 21: PlantYear)
    IF (model_config%stand_replacing_year .LT. 1901) THEN
      model_config%flag_stand_harvest = .FALSE.
    END IF
    
    ! set default for other config% parameter that are changed during runtime not via namelist
    model_config%l_do_stand_replacing_harvest           = .FALSE.

    DO i=1,99
      IF (output_tiles(i) == '') EXIT
    END DO
    i = i - 1  ! Number of tile names in namelist variable output_tiles
    IF (i == 0) THEN             ! Default is 'box' if no tiles specified in namelist
      i = 1
      output_tiles(1) = 'box'
    END IF
    ALLOCATE(model_config%output_tiles(i))
    model_config%output_tiles(1:i) = output_tiles(1:i)

  END FUNCTION new_model_config

#endif
END MODULE mo_jsb_config_class
