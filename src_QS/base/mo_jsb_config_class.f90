!> Contains structures and methods for JSBACH model config
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_config_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_io_units,          ONLY: filename_max
  USE mo_exception,         ONLY: message, finish, message_text
  USE mo_jsb_io_netcdf,     ONLY: t_input_file

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
    CHARACTER(len=20)              :: model_scheme_char   !< model_scheme from namelist (used in mo_jsb_lctlib_class)
    INTEGER                        :: model_scheme        !< model_scheme as integer
    LOGICAL                        :: use_tmx
    LOGICAL                        :: use_lakes
    LOGICAL                        :: use_glacier
    INTEGER                        :: enforce_water_budget!< descriptor for water balance 'ignore', 'logging', or 'error'
    CHARACTER(len=10)              :: tpe_scheme       ! For terraplanet setup: open/closed
    CHARACTER(len=10)              :: hsm_mode
    CHARACTER(len=10), ALLOCATABLE :: output_tiles(:)  ! List of tile names for which output should be generated
    
    ! Characteristics of bgc material in this model: flags specifying used elements and bookkeeping vector for elements used
    ! Note: no logical for carbon because we presuppose that it is always used!
    LOGICAL                        :: force_c_isotopes      !< If c isotopes forcing should be used (else use atmland constant)
    LOGICAL                        :: include_nitrogen      !< If bgc materials contain nitrogen as element
    LOGICAL                        :: include_phosphorus    !< If bgc materials contain phosphorus as element
    LOGICAL                        :: include_carbon13      !< If bgc materials contain C13 as element
    LOGICAL                        :: include_carbon14      !< If bgc materials contain C14 as element
    LOGICAL                        :: include_nitrogen15    !< If bgc materials contain C15 as element
    INTEGER                        :: nr_of_elements        !< for bookkeeping: number of elements
    INTEGER                        :: nr_of_used_elements   !< for bookkeeping: number of actually used elements
    INTEGER, ALLOCATABLE           :: elements_index_map(:) !< vector for ID -> IDX mapping of the elements
    LOGICAL, ALLOCATABLE           :: is_element_used(:)    !< vector(ID) == .TRUE. if "element with ID" is used

    ! slow sb pools spin-up accelerator
    LOGICAL                     :: flag_slow_sb_pool_spinup_accelerator !< accelerating slow pool turnover during spin-up
    INTEGER                     :: slow_sb_pool_spinup_accelerator_frequency !< The freqency [years] to execute spin-up accelerator
    INTEGER                     :: slow_sb_pool_spinup_accelerator_length !< The length (loop times) of spin-up accelerator
    INTEGER                     :: slow_sb_pool_spinup_accelerator_start_year !< The year in which the spin-up accelerator should first be executed
    INTEGER                     :: slow_sb_pool_spinup_accelerator_max_executions !< maximum number of executions of spin-up accelerator
    
    !! QUINCY added !!
    !! BASE_CTL
    REAL(wp)          :: dtime                                    !< length of one time step in seconds; base_ctl namelist: dtime_step_length_sec
    CHARACTER(len=20) :: quincy_model_name                        !< quincy model-mode [land, plant, soil, canopy, test_canopy, test_radiation]
    LOGICAL           :: use_soil_phys_jsbach                     !< use either quincy SPQ_ or jsbach soil-physic processes
    CHARACTER(len=50) :: git_branch                               !< branch for which the binary was compiled
    CHARACTER(len=8)  :: git_commit_SHA                           !< sha of commit used for this binary
    LOGICAL           :: code_has_changed_since_SHA               !< logical indicating if there has been any changes since last commit
    CHARACTER(len=20) :: simulation_run_by                        !< name of the user making the simulation
    CHARACTER(len=filename_max) :: file_sel_output_variables      !< file (and path) with list of selected variables for netcdf output
    CHARACTER(len=20) :: output_txt_file_set                      !< set of ascii txt files for model output (options: mo_qs_output_txt_util)
    CHARACTER(8)      :: output_interval_pool, &                  !< interval at which output is generated during model runs not including spinup
                         output_interval_flux                     !< interval at which output is generated during model runs not including spinup - the pool value is used
    CHARACTER(8)      :: output_interval_pool_spinup, &           !< interval at which output is generated during spinup
                         output_interval_flux_spinup              !< interval at which output is generated during spinup - the pool value is used
    INTEGER           :: output_start_first_day_year, &           !< output start at the 1st  day of this year
                         output_end_last_day_year                 !< output stop  at the last day of this year
    INTEGER           :: forcing_file_start_yr                    !< first year of forcing data available with the forcing file
    INTEGER           :: forcing_file_last_yr                     !< last  year of forcing data available with the forcing file
    LOGICAL           :: fluxnet_type_transient_timestep_output   !< transient simulation with fluxnet-type site-set and timestep-output for the last full forcing period (static forcing)
    INTEGER           :: fluxnet_static_forc_start_yr             !< fluxnet sites: first year of (static) forcing available  - used for transient simulations
    INTEGER           :: fluxnet_static_forc_last_yr              !< fluxnet sites: last year with (static) forcing available - used for transient simulations
    CHARACTER(len=20) :: sim_siteset_specific                     !< specific simulation setup for particular site-set, e.g., EucFACE transient simulations
    LOGICAL           :: flag_generate_soil_forcing_stat          !< usecase 01: all litterfall pools + some soilphysics | static simulation
                                                                  !! at period: model%config%output_start_first_day_year .. output_end_last_day_year
    LOGICAL           :: flag_generate_soil_forcing_trans         !< usecase 08: all litterfall pools + some soilphysics | transient simaulation
    LOGICAL           :: flag_write_sb_output                     !< with netcdf output this can be used to switch off SB output
    LOGICAL           :: flag_write_veg_output                    !< with netcdf output this can be used to switch off VEG output
    LOGICAL           :: set_parameter_values_from_file           !< set values of (some) parameters to values given by namelist (parameter_slm_run.list)
    LOGICAL           :: set_parameter_values_proportional        !< if TRUE: change the default (hardcoded) parameter values by proportion as defined in namelist
    LOGICAL           :: flag_climate_forcing_one_year
    LOGICAL           :: flag_spp1685
    ! GRID_CTL
    REAL(wp)          :: lat
    REAL(wp)          :: lon
    !! other parameters than namelist parameters => not set via namelist but at model runtime
    INTEGER           :: qmodel_id                          !< defines the applied model configuration (quincy_model_name)
    LOGICAL           :: l_transient_forcing_mode, &        !< set to TRUE if forcing_options%forcing_mode == transient
                         l_transient_spinup, &              !< TRUE for forcing_mode=='transient' and model spinup, otherwise FALSE
                         l_transient_spinup_rewind_forcing  !< TRUE at the timestep after one transient spinup period has finished (e.g. after 1930 for 1901-1930)
    LOGICAL           :: l_static_forcing_mode              !< set to TRUE if forcing_options%forcing_mode == static

    LOGICAL                        :: init_from_ifs    ! Initialize from IFS analysis; quincy standalone: hardcoded FALSE
    CHARACTER(len=:), ALLOCATABLE  :: ifs_filename
    TYPE(t_input_file)             :: ifs_input_file
  END TYPE t_jsb_model_config


  !-----------------------------------------------------------------------------------------------------
  !> Abstract TYPE for process and model configuration
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, ABSTRACT :: t_jsb_config
    LOGICAL :: active
    LOGICAL :: lrestart_cont
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename
    CHARACTER(LEN=filename_max) :: namelist_filename
    TYPE(t_jsb_model_config), POINTER :: model_config => NULL()
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

    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_impl_constants,    ONLY: WB_IGNORE, WB_LOGGING, WB_ERROR
    USE mo_util,                  ONLY: tolower
    USE mo_quincy_model_config,   ONLY: Get_quincy_model_config_id

    CHARACTER(len=*),         INTENT(in) :: namelist_filename
    TYPE(t_jsb_model_config), POINTER    :: model_config

    ! GRID_CTL
    REAL(wp)                    :: latitude
    REAL(wp)                    :: longitude
    CHARACTER(len=filename_max) :: grid_filename
    CHARACTER(len=30)           :: grid_name
    INTEGER                     :: MODEL_JSBACH, MODEL_QUINCY

    NAMELIST /grid_ctl/                         &
        latitude                              , &
        longitude
    ! NAMELIST /jsb_grid_nml/   &
    !   grid_filename,          &
    !   grid_name

    ! BASE_CTL
    REAL(wp)                    :: dtime_step_length_sec                      !< model-variable: dtime 
    CHARACTER(len=20)           :: quincy_model_name
    LOGICAL                     :: use_soil_phys_jsbach
    CHARACTER(len=80)           :: git_branch                               !< branch for which the binary was compiled
    CHARACTER(len=8)            :: git_commit_SHA                           !< sha of commit used for this binary
    LOGICAL                     :: code_has_changed_since_SHA               !< logical indicating if there has been any changes since last commit
    CHARACTER(len=20)           :: simulation_run_by
    CHARACTER(len=filename_max) :: file_sel_output_variables
    CHARACTER(len=20)           :: output_txt_file_set
    CHARACTER(8)                :: output_interval_pool          , &
                                   output_interval_flux
    CHARACTER(8)                :: output_interval_pool_spinup   , & 
                                   output_interval_flux_spinup
    INTEGER                     :: output_start_first_day_year   , &
                                   output_end_last_day_year
    INTEGER                     :: forcing_file_start_yr
    INTEGER                     :: forcing_file_last_yr
    LOGICAL                     :: fluxnet_type_transient_timestep_output
    INTEGER                     :: fluxnet_static_forc_start_yr
    INTEGER                     :: fluxnet_static_forc_last_yr
    CHARACTER(len=20)           :: sim_siteset_specific
    LOGICAL                     :: flag_generate_soil_forcing_stat
    LOGICAL                     :: flag_generate_soil_forcing_trans
    LOGICAL                     :: flag_write_sb_output
    LOGICAL                     :: flag_write_veg_output
    LOGICAL                     :: set_parameter_values_from_file
    LOGICAL                     :: set_parameter_values_proportional
    LOGICAL                     :: flag_climate_forcing_one_year
    LOGICAL                     :: flag_spp1685

    LOGICAL                     :: include_nitrogen           !< element variables in bgc_material (infrastructure)
    LOGICAL                     :: include_phosphorus         !< element variables in bgc_material (infrastructure)
    LOGICAL                     :: include_carbon13           !< element variables in bgc_material (infrastructure)
    LOGICAL                     :: include_carbon14           !< element variables in bgc_material (infrastructure)
    LOGICAL                     :: include_nitrogen15         !< element variables in bgc_material (infrastructure)

    LOGICAL                     :: flag_slow_sb_pool_spinup_accelerator !< accelerating slow pool turnover during spin-up
    INTEGER                     :: slow_sb_pool_spinup_accelerator_frequency !< The freqency [years] to execute spin-up accelerator
    INTEGER                     :: slow_sb_pool_spinup_accelerator_length !< The length (loop times) of spin-up accelerator
    INTEGER                     :: slow_sb_pool_spinup_accelerator_start_year !< The year in which the spin-up accelerator should first be executed
    INTEGER                     :: slow_sb_pool_spinup_accelerator_max_executions !< maximum number of executions of spin-up accelerator

    CHARACTER(len=20)           :: model_scheme     !< model_scheme: jsbach / jsbach_quincy / quincy
    LOGICAL                     :: use_tmx
    LOGICAL                     :: use_lakes
    LOGICAL                     :: use_glacier
    CHARACTER(len=10)           :: enforce_water_budget
    CHARACTER(len=filename_max) :: fract_filename
    CHARACTER(len=40)           :: usecase
    CHARACTER(len=10)           :: tpe_scheme
    CHARACTER(len=10)           :: hsm_mode
    CHARACTER(len=10)           :: output_tiles(99)
    LOGICAL                     :: init_from_ifs
    CHARACTER(len=filename_max) :: ifs_filename

    NAMELIST /base_ctl/                           &  
        dtime_step_length_sec                   , &
        quincy_model_name                       , &
        use_soil_phys_jsbach                    , &
        git_branch                              , &
        git_commit_SHA                          , &
        code_has_changed_since_SHA              , &
        simulation_run_by                       , &
        file_sel_output_variables               , &
        output_txt_file_set                     , &
        output_interval_pool                    , &
        output_interval_flux                    , &
        output_interval_pool_spinup             , &
        output_interval_flux_spinup             , &
        output_start_first_day_year             , &
        output_end_last_day_year                , &
        forcing_file_start_yr                   , &
        forcing_file_last_yr                    , &
        fluxnet_type_transient_timestep_output  , &
        fluxnet_static_forc_start_yr            , &
        fluxnet_static_forc_last_yr             , &
        sim_siteset_specific                    , &
        flag_generate_soil_forcing_stat         , &
        flag_generate_soil_forcing_trans        , &
        flag_write_sb_output                    , &
        flag_write_veg_output                   , &
        set_parameter_values_from_file          , &
        set_parameter_values_proportional       , &
        flag_climate_forcing_one_year           , &
        flag_spp1685, &
        include_nitrogen, &
        include_phosphorus, &
        include_carbon13, &
        include_carbon14, &
        include_nitrogen15, &
        flag_slow_sb_pool_spinup_accelerator, &
        slow_sb_pool_spinup_accelerator_frequency, &
        slow_sb_pool_spinup_accelerator_length, &
        slow_sb_pool_spinup_accelerator_start_year, &
        slow_sb_pool_spinup_accelerator_max_executions, &
        enforce_water_budget
                          
    ! NAMELIST /jsb_model_nml/  &
    !   use_tmx, &
    !   use_lakes, &
    !   use_glacier, &
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
    quincy_model_name                       = 'land'      ! default is the coupled veg-soil model, a.k.a. land
#ifdef __QUINCY_STANDALONE__
    use_soil_phys_jsbach                    = .FALSE.     ! with quincy standalone the default is to not use the "jsbach soil physic" but SPQ_
#else
    use_soil_phys_jsbach                    = .TRUE.
#endif
    git_branch                              = 'not_set'   ! branch for which the binary was compiled
    git_commit_SHA                          = 'not_set'   ! commit SHA of the executed binary
    simulation_run_by                       = 'not_set'   ! person who executed the simulation
    file_sel_output_variables               = 'not_set'   ! file (with path) containing a list of selected output variables
    output_txt_file_set                     = 'basic'     ! set of ascii txt files for model output (options: mo_qs_output_txt_util)
    code_has_changed_since_SHA              = .TRUE.      ! whether the code of the binary changed since SHA
    output_interval_pool                    = 'weekly'    ! timestep daily weekly monthly yearly
    output_interval_flux                    = 'weekly'    ! timestep daily weekly monthly yearly
    output_interval_pool_spinup             = 'yearly'    ! usually ~500 yrs are covered, so 'yearly' would be reasonable
    output_interval_flux_spinup             = 'yearly'    ! usually ~500 yrs are covered, so 'yearly' would be reasonable
    output_start_first_day_year             = 1
    output_end_last_day_year                = 10
    forcing_file_start_yr                   = 1500        ! arbitrary default
    forcing_file_last_yr                    = 1500        ! arbitrary default
    fluxnet_type_transient_timestep_output  = .FALSE.     ! .TRUE. .FALSE.
    fluxnet_static_forc_start_yr            = 1500        ! usually a year after 1990
    fluxnet_static_forc_last_yr             = 1500        ! usually a year after 1990
    sim_siteset_specific                    = 'none'      ! default is 'none'
    flag_generate_soil_forcing_stat         = .FALSE.     ! .TRUE. .FALSE.  << is set to T for usecase 01 only
    flag_generate_soil_forcing_trans        = .FALSE.     ! .TRUE. .FALSE.  << is set to T for usecase 08 only
    flag_write_sb_output                    = .TRUE.      ! with netcdf output this can be used to switch off SB output
    flag_write_veg_output                   = .TRUE.      ! with netcdf output this can be used to switch off VEG output
    set_parameter_values_from_file          = .FALSE.     ! .TRUE. .FALSE.
    set_parameter_values_proportional       = .FALSE.     ! .TRUE. .FALSE.
    flag_climate_forcing_one_year           = .FALSE.     ! .TRUE. .FALSE.
    ! GRID_CTL
    latitude                                = 50.91       ! about the MPI-BGC Jena
    longitude                               = 11.56       ! about the MPI-BGC Jena
    ! JSBACH4
    ! model_scheme values
    ! define here as INTEGER replacing the ENUM "USE mo_jsb_model_class,    ONLY: MODEL_JSBACH, MODEL_QUINCY"
    ! this is needed to avoid a dependency cycle, creating an compile error (tested with NAG):
    !  stating that in mo_jsb_memory_class the mo_jsb_config_class is not available for "USE mo_jsb_config_class,  ONLY: t_jsb_config"
    MODEL_JSBACH               = 1
    MODEL_QUINCY               = 2
    model_scheme               = "quincy"  ! "jsbach"
    use_tmx                                 = .FALSE.
    use_lakes                               = .FALSE.
    use_glacier                             = .FALSE.
    grid_filename                           = ''
    grid_name                               = ''
    usecase                                 = ''
    tpe_scheme                              = ''
    hsm_mode                                = 'simple'
    fract_filename                          = 'bc_land_frac.nc'
    output_tiles(:)                         = ''
    flag_spp1685                            = .FALSE.       ! if TRUE: soil layer specific texture information for SPP_1685 sites: simple_1d and jsm
    include_nitrogen                        = .TRUE.
    include_phosphorus                      = .TRUE.
    include_carbon13                        = .TRUE.
    include_carbon14                        = .TRUE.
    include_nitrogen15                      = .TRUE.

    flag_slow_sb_pool_spinup_accelerator    = .FALSE.
    slow_sb_pool_spinup_accelerator_frequency = 100
    slow_sb_pool_spinup_accelerator_length = 1000
    slow_sb_pool_spinup_accelerator_start_year = 300
    slow_sb_pool_spinup_accelerator_max_executions = 4

    enforce_water_budget                    = 'ignore'
    init_from_ifs                           = .FALSE.
    ifs_filename                            = 'ifs2icon.nc'

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
    model_config%quincy_model_name                      = quincy_model_name
    model_config%use_soil_phys_jsbach                   = use_soil_phys_jsbach
    model_config%git_branch                             = git_branch
    model_config%git_commit_SHA                         = git_commit_SHA
    model_config%code_has_changed_since_SHA             = code_has_changed_since_SHA
    model_config%simulation_run_by                      = simulation_run_by
    model_config%file_sel_output_variables              = file_sel_output_variables
    model_config%output_txt_file_set                    = TRIM(output_txt_file_set)
    model_config%output_interval_pool                   = output_interval_pool
    model_config%output_interval_flux                   = output_interval_flux
    model_config%output_interval_pool_spinup            = output_interval_pool_spinup
    model_config%output_interval_flux_spinup            = output_interval_flux_spinup
    model_config%output_start_first_day_year            = output_start_first_day_year
    model_config%output_end_last_day_year               = output_end_last_day_year
    model_config%forcing_file_start_yr                  = forcing_file_start_yr 
    model_config%forcing_file_last_yr                   = forcing_file_last_yr  
    model_config%fluxnet_type_transient_timestep_output = fluxnet_type_transient_timestep_output
    model_config%fluxnet_static_forc_start_yr           = fluxnet_static_forc_start_yr
    model_config%fluxnet_static_forc_last_yr            = fluxnet_static_forc_last_yr
    model_config%sim_siteset_specific                   = sim_siteset_specific
    model_config%flag_generate_soil_forcing_stat        = flag_generate_soil_forcing_stat
    model_config%flag_generate_soil_forcing_trans       = flag_generate_soil_forcing_trans
    model_config%flag_write_sb_output                   = flag_write_sb_output
    model_config%flag_write_veg_output                  = flag_write_veg_output
    model_config%set_parameter_values_from_file         = set_parameter_values_from_file
    model_config%set_parameter_values_proportional      = set_parameter_values_proportional
    model_config%flag_climate_forcing_one_year          = flag_climate_forcing_one_year
    model_config%flag_spp1685                           = flag_spp1685
    ! organise elements in bgc_material
    model_config%include_nitrogen                       = include_nitrogen
    model_config%include_phosphorus                     = include_phosphorus
    model_config%include_carbon13                       = include_carbon13
    model_config%include_carbon14                       = include_carbon14
    model_config%include_nitrogen15                     = include_nitrogen15
    CALL setup_element_index_map_in_config(model_config, routine)

    model_config%flag_slow_sb_pool_spinup_accelerator   = flag_slow_sb_pool_spinup_accelerator
    model_config%slow_sb_pool_spinup_accelerator_frequency = slow_sb_pool_spinup_accelerator_frequency
    model_config%slow_sb_pool_spinup_accelerator_length = slow_sb_pool_spinup_accelerator_length
    model_config%slow_sb_pool_spinup_accelerator_start_year = slow_sb_pool_spinup_accelerator_start_year
    model_config%slow_sb_pool_spinup_accelerator_max_executions = slow_sb_pool_spinup_accelerator_max_executions
    
    model_config%model_scheme_char = tolower(model_scheme)
    SELECT CASE (TRIM(model_config%model_scheme_char))
    CASE ("jsbach")
      model_config%model_scheme = MODEL_JSBACH
      CALL message(TRIM(routine), 'Use model_scheme JSBACH')
    CASE ("quincy")
      model_config%model_scheme = MODEL_QUINCY
      CALL message(TRIM(routine), 'Use model_scheme QUINCY')
    CASE DEFAULT
      CALL finish(TRIM(routine), 'Selected model_scheme '//TRIM(model_scheme)//' not available.')
    END SELECT
    model_config%use_tmx        = use_tmx
    model_config%init_from_ifs  = init_from_ifs
    ! model_config%use_lakes      = use_lakes
    ! model_config%use_glacier    = use_glacier
    ! model_config%usecase        = usecase
    ! model_config%tpe_scheme     = tpe_scheme
    ! model_config%hsm_mode       = hsm_mode
    ! model_config%fract_filename = fract_filename

    ! set the model configuration id
    model_config%qmodel_id = Get_quincy_model_config_id(model_config%quincy_model_name)

    SELECT CASE (tolower(TRIM(enforce_water_budget)))
    CASE ("ignore")
      model_config%enforce_water_budget = WB_IGNORE
      CALL message(TRIM(routine), 'WARNING: Land surface water balance will not be checked during simulation '// &
        & '(unless set differently in hydro or hd config).')
    CASE ("logging")
      model_config%enforce_water_budget = WB_LOGGING
      CALL message(TRIM(routine), 'WARNING: Simulation will not stop due to any land water balance violation '// &
        & '(unless set differently in hydro or hd config) but information will be added to the log file')
    CASE ("error")
      model_config%enforce_water_budget = WB_ERROR
      CALL message(TRIM(routine), 'WARNING: Simulation will stop due to any land water balance violation '// &
        & '(unless set differently in hydro or hd config).')
    CASE DEFAULT
      CALL finish(TRIM(routine), 'enforce_water_budget == '//tolower(TRIM(enforce_water_budget))//' not available.')
    END SELECT

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

  ! ======================================================================================================= !
  !>
  !> Uses the element flags to setup the element map and to identify the number of used elements
  !>
  SUBROUTINE setup_element_index_map_in_config(config, caller)
    USE mo_lnd_bgcm_class,        ONLY: ELEM_C_ID, ELEM_N_ID, ELEM_P_ID,  &
      &                                 ELEM_C13_ID, ELEM_C14_ID, ELEM_N15_ID, LAST_ELEM_ID
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model_config), INTENT(INOUT) :: config  !< model config to be setup for used elements
    CHARACTER(len=*),            INTENT(IN) :: caller  !< calling routine
    ! -------------------------------------------------------------------------------------------------- !
    ! allocate map with size "number of element IDs defined in mo_lnd_bgcm_class"
    config%nr_of_elements = LAST_ELEM_ID
    ALLOCATE(config%elements_index_map(config%nr_of_elements))
    config%elements_index_map(:) = -1  ! set all IND to '-1' by default, -1 for "element not used"
    config%nr_of_used_elements = 0     ! init for use in setup_element
    ! allocate vector with info about whether the element is used
    ALLOCATE(config%is_element_used(config%nr_of_elements))
    config%is_element_used(:) = .FALSE. ! set all elements to "not used" by default
    ! information for logfile
    CALL message(TRIM(caller), 'Mapping of element ID -> IND:')

    ! carbon 
    ! is presumed to always be used, with ID = IND = 1 (thus call first for carbon!)
    CALL setup_element(config, .TRUE., ELEM_C_ID, 'C', caller)
    ! nitrogen
    CALL setup_element(config, config%include_nitrogen, ELEM_N_ID, 'N', caller)
    ! phosphorus
    CALL setup_element(config, config%include_phosphorus, ELEM_P_ID, 'P', caller)
    ! carbon13
    CALL setup_element(config, config%include_carbon13, ELEM_C13_ID, 'C13', caller)
    ! carbon14
    CALL setup_element(config, config%include_carbon14, ELEM_C14_ID, 'C14', caller)
    ! nitrogen15
    CALL setup_element(config, config%include_nitrogen15, ELEM_N15_ID, 'N15', caller)
    
  END SUBROUTINE setup_element_index_map_in_config
  
  ! ======================================================================================================= !
  !>
  !> Small helper to setup a given element according to its given flag
  !>
  SUBROUTINE setup_element(config, elem_flag, elem_id, short_name, caller)
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model_config), INTENT(INOUT) :: config       !< model config to be setup for used elements
    LOGICAL,                     INTENT(IN) :: elem_flag    !< flag for element to be setup
    INTEGER,                     INTENT(IN) :: elem_id      !< id of element to be setup
    CHARACTER(len=*),            INTENT(IN) :: short_name   !< a short name for the element
    CHARACTER(len=*),            INTENT(IN) :: caller       !< calling routine
    ! -------------------------------------------------------------------------------------------------- !
    ! information for logfile
    CALL message(TRIM(caller), 'Setup element:')

    IF (elem_flag) THEN
      config%nr_of_used_elements         = config%nr_of_used_elements + 1
      config%elements_index_map(elem_id) = config%nr_of_used_elements
      config%is_element_used(elem_id)    = .TRUE.
    END IF

    WRITE (message_text,*) elem_id, " -> ", config%elements_index_map(elem_id)
    CALL message(TRIM(short_name), TRIM(message_text))

  END SUBROUTINE setup_element

#endif
END MODULE mo_jsb_config_class
