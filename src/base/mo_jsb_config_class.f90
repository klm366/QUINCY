!> Contains structures and methods for the ICON-Land model configuration
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
!> Contains structures and methods for the ICON-Land model configuration
!>
!> Two configuration types are handled in this module: the basic model configuration ([[t_jsb_model_config]])
!> and [[t_jsb_config]], the part of the individual process configurations common to all ICON-Land processes.
!>
MODULE mo_jsb_config_class
#ifndef __NO_JSBACH__
  !>
  !> @TODO `filename_max` and explicit string length definitions should not be needed with
  !>       `CHARACTER(len=:), ALLOCATABLE :: ` syntax.
  !>
  USE mo_kind,               ONLY: lp
  USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN
  USE mo_io_units,           ONLY: filename_max
  USE mo_exception,          ONLY: message, finish, message_text
  USE mo_jsb_io_netcdf,      ONLY: t_input_file

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_config, t_jsb_config_p, t_jsb_model_config, new_model_config

  !>
  !> Type for basic model configuration
  !>
  !> @TODO Are `nr_of_elements`, `nr_of_used_elements`, `elements_index_map(:)`, and `is_element_used(:)` needed in the model config? They are not configurable!
  !>
  !> @TODO Check if `grid_filename` and `grid_name` variables are needed at all.
  !>
  !> @TODO remove `relative_fractions_in_file` function. We agreed to always use relative fractions in the initial file, although absolute (grid box)
  !>      fractions are used elsewhere in the code.
  !>
  !> @TODO: Are the `nr_of_elements`, `nr_of_used_elements`,
  !> `elements_index_map(:)`, and `is_element_used(:)` variables needed
  !> in the model config? They are not configurable!
  !>
  !>
  TYPE t_jsb_model_config

    ! TODO Check if `grid_filename` and `grid_name` variables are needed at all.
    CHARACTER(len=filename_max)    :: grid_filename       !< Model grid file name
    CHARACTER(len=filename_max)    :: fract_filename      !< Name of file containing tile fractions
    CHARACTER(len=30)              :: grid_name           !< Model grid name
    CHARACTER(len=40)              :: usecase             !< Model use case
    CHARACTER(len=20)              :: model_scheme_char   !< Model variant: 'jsbach' or 'quincy'
    INTEGER                        :: model_scheme        !< Model variant identifier
    INTEGER                        :: npft                !< Number of PFTs
    LOGICAL                        :: use_tmx             !< Use turbulent mixing scheme
    LOGICAL                        :: use_lakes           !< Lakes are represented
    LOGICAL                        :: use_glacier         !< Glaciers are represented
    LOGICAL                        :: check_child_cf      !< Assert that cover fractions of children tiles sum up to the tile cf
    INTEGER                        :: enforce_water_budget!< How to deal with water balance errors?</br>
                                                          !< 'ignore' (0):  short, one-line warning,</br>
                                                          !< 'logging' (1): detailed debugging output,</br>
                                                          !< 'error' (2):   debugging output and model stop.
                                                          !< Note: HYDRO parameter [[t_hydro_config:enforce_water_budget]]
                                                          !< can override this setting for the HYDRO process.
    CHARACTER(len=10)              :: quincy_model        !< Quincy model configuration: 'canopy', 'plant'
                                                          !< or 'land'
    INTEGER                        :: qmodel_id           !< Identifier for the quincy model configuration
    CHARACTER(len=SHORT_NAME_LEN)  :: tpe_scheme          !< Switch for Terraplanet setup: Handling of water
                                                          !< cycle: 'open' or 'closed'
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: output_tiles(:)  !< List of the tiles for which output will
                                                                   !< be generated (default: 'box')
    LOGICAL                        :: relative_fractions_in_file   !< Tile fractions in [[t_jsb_model_config:fract_filename]]
                                                                   !< are relative to the parent tiles (true)
                                                                   !< or relative to the grid box (false).
    ! TODO remove `relative_fractions_in_file` function. We agreed to always use relative fractions in the initial file, although absolute (grid box)
    !      fractions are used elsewhere in the code.

    ! Characteristics of biogeochemical (bgc) material in quincy: flags specifying the elements used and
    ! corresponding bookkeeping vectors.
    ! Note: No logical for carbon because we presuppose that it is always used.
    LOGICAL                        :: force_c_isotopes      !< `.TRUE.`: Read carbon isotopes forcing from file. `.FALSE.`: use constant values defined in [[mo_atmland_constants]]
    LOGICAL                        :: include_nitrogen      !< bgc materials contain nitrogen
    LOGICAL                        :: include_phosphorus    !< bgc materials contain phosphorus
    LOGICAL                        :: include_carbon13      !< bgc materials contain C13
    LOGICAL                        :: include_carbon14      !< bgc materials contain C14
    LOGICAL                        :: include_nitrogen15    !< bgc materials contain N15
    LOGICAL                        :: l_spinup              !< Running a spin-up simulation
    ! TODO: Are the `nr_of_elements`, `nr_of_used_elements`,
    ! `elements_index_map(:)`, and `is_element_used(:)` variables needed
    ! in the model config? They are not configurable!
    INTEGER                        :: nr_of_elements        !< For bookkeeping: number of elements in bgc materials
    INTEGER                        :: nr_of_used_elements   !< For bookkeeping: number of actually used elements
    INTEGER, ALLOCATABLE           :: elements_index_map(:) !< Vector for ID -> IDX mapping of the elements
    LOGICAL, ALLOCATABLE           :: is_element_used(:)    !< Vector(ID) == .TRUE. if "element with ID" is used

    ! Spinup accelerator for the slow soil biogeorchemical (sb) pools
    LOGICAL ::  flag_slow_sb_pool_spinup_accelerator           !< Accelerate slow pool turnover (for spin-up)
    INTEGER ::  slow_sb_pool_spinup_accelerator_frequency      !< Frequency of spin-up accelerator runs [years]
    INTEGER ::  slow_sb_pool_spinup_accelerator_length         !< Number of time loops in spin-up accelerator
    INTEGER ::  slow_sb_pool_spinup_accelerator_start_year     !< Year of the first spin-up accelerator runs
    INTEGER ::  slow_sb_pool_spinup_accelerator_max_executions !< Maximum number of spin-up accelerator runs

    ! Initialization from ECMWF Integrated Forecast System (IFS) analysis data
    LOGICAL                        :: init_from_ifs            !< Initialize from IFS analysis
    CHARACTER(len=:), ALLOCATABLE  :: ifs_filename             !< IFS input data file name
    TYPE(t_input_file)             :: ifs_input_file           !< Structure for the IFS initial file

    ! Single site setup
    LOGICAL                        :: is_single_site        !< Single site model
    CHARACTER(len=:), ALLOCATABLE  :: single_site_name      !< Single site name
    CHARACTER(len=:), ALLOCATABLE  :: single_site_shortname !< Single site short name
    REAL(lp)                       :: single_site_lonlat(2) !< Lon/lat of site in degrees

  END TYPE t_jsb_model_config

  !> Abstract type for process configurations
  TYPE, ABSTRACT :: t_jsb_config
    LOGICAL :: active                                          !< The process is active in the current setup.
    INTEGER :: run_case
    LOGICAL :: lrestart_cont                                   !< Continue although process variables are
                                                               !< missing in the restart file.
    CHARACTER(len=filename_max) :: ic_filename                 !< Name of the initial file
    CHARACTER(len=filename_max) :: bc_filename                 !< Name of the boundary conditions file
    CHARACTER(LEN=filename_max) :: namelist_filename           !< File containing the process namelist
                                                               !< **Note**: generally all ICON-Land namelists are
                                                               !< listed in one common namelist file
    TYPE(t_jsb_model_config), POINTER :: model_config => NULL() !< Initialize model configuration pointer
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: output_tiles(:) !< List of tiles for which this process
                                                                !< will write output
  CONTAINS
    !> Each process needs to define an init function as extension of this abstract type.
    PROCEDURE (Init_config), DEFERRED, PASS(config) :: Init
  END TYPE t_jsb_config

  !> Pointer type fot a `t_jsb_config` type
  TYPE t_jsb_config_p
    CLASS(t_jsb_config), POINTER :: p
  END TYPE t_jsb_config_p

  ABSTRACT INTERFACE
    !> Definition of the interface of the init function (needed with deferred procedures).
    !> The users are enforced to use the identical calling structure for all processes.
    SUBROUTINE Init_config(config)
      IMPORT :: t_jsb_config
      CLASS(t_jsb_config), INTENT(inout) :: config
    END SUBROUTINE Init_config
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_config_class' !< Name of the module

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Creates the model configuration
  !>
  !> The function reads the model configuration from namelist "jsb_model_nml" defined in the given namelist
  !> file and returns a pointer to the model configuration structure.
  !>
  !> @TODO Are `grid_filename`, `grid_name` and the `jsb_grid_nml` needed at all?
  !>
  !> @TODO Dissolve the dependency<br>
  !> `model_scheme` values are
  !> defined here as `INTEGER` replacing the ENUM `USE mo_jsb_model_class,    ONLY: MODEL_JSBACH, MODEL_QUINCY`
  !> This is needed to avoid a dependency cycle, creating a compile error (tested with NAG):
  !>  stating that in `mo_jsb_memory_class` the `mo_jsb_config_class` is not available for `USE mo_jsb_config_class,  ONLY: t_jsb_config`
  FUNCTION new_model_config(namelist_filename) RESULT(model_config)

    USE mo_jsb_namelist_iface,  ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_impl_constants,  ONLY: WB_IGNORE, WB_LOGGING, WB_ERROR    ! handling of water balance errors
    USE mo_jsb_control,         ONLY: jsbach_runs_standalone
    USE mo_util,                ONLY: tolower
    USE mo_quincy_model_config, ONLY: QLAND, Get_quincy_model_config_id

    CHARACTER(len=*),         INTENT(in) :: namelist_filename     !< File containing the model namelist
    TYPE(t_jsb_model_config), POINTER    :: model_config          !< Model configuration: structure
                                                                  !< containing the configurable parameters
    ! TODO Are `grid_filename`, `grid_name` and the `jsb_grid_nml` needed at all?
    CHARACTER(len=filename_max) :: grid_filename                  !< File with the model grid file
    CHARACTER(len=30)           :: grid_name                      !< Model grid name

    INTEGER                     :: MODEL_JSBACH, MODEL_QUINCY     !< Identifiers for the model scheme

    NAMELIST /jsb_grid_nml/   &
      grid_filename,          &
      grid_name

    !> Variables defined in [[t_jsb_model_config]]
    CHARACTER(len=20)             :: model_scheme                 !< Model scheme: 'jsbach' or 'quincy'
    LOGICAL                       :: use_tmx                      !< Use turbulent mixing scheme
    LOGICAL                       :: use_lakes                    !< Represent lakes
    LOGICAL                       :: use_glacier                  !< Represent glaciers
    LOGICAL                       :: check_child_cf               !< Assert that cover fractions of children tiles sum up to the tile cf
    CHARACTER(len=10)             :: enforce_water_budget         !< Consequence of water balance errors
    CHARACTER(len=10)             :: quincy_model                 !< Quincy model version: 'canopy', 'plant'
                                                                  !< or 'land'
    INTEGER                       :: qmodel_id                    !< Identifier for the quincy model version
    CHARACTER(len=filename_max)   :: fract_filename               !< File with the tile fractions
    CHARACTER(len=40)             :: usecase                      !< ICON-Land use case
    INTEGER                       :: npft                         !< Number of PFTs
    CHARACTER(len=SHORT_NAME_LEN) :: tpe_scheme                   !< Terraplanet setup: Handling of water
                                                                  !< cycle: 'open' or 'closed'
    CHARACTER(len=SHORT_NAME_LEN) :: output_tiles(99)             !< List of the output tiles
    LOGICAL                       :: relative_fractions_in_file   !< Tile fractions in [[t_jsb_model_config:fract_filename]]
                                                                  !< are relative to the parent tiles.
    LOGICAL :: force_c_isotopes                                   !< Read carbon isotopes forcing
    LOGICAL :: include_nitrogen                                   !< Nitrogen is element in bgc_material
    LOGICAL :: include_phosphorus                                 !< Phosphorus is element in bgc_material
    LOGICAL :: include_carbon13                                   !< C13 is element in bgc_material
    LOGICAL :: include_carbon14                                   !< C14 is element in bgc_material
    LOGICAL :: include_nitrogen15                                 !< N15 is element in bgc_material
    LOGICAL :: l_spinup                                           !< Running a spin-up simulation
    LOGICAL :: flag_slow_sb_pool_spinup_accelerator               !< Accelerate slow pool turnover during spin-up
    INTEGER :: slow_sb_pool_spinup_accelerator_frequency          !< Frequency of spin-up accelerator runs [years]
    INTEGER :: slow_sb_pool_spinup_accelerator_length             !< The length (loop times) of spin-up accelerator
    INTEGER :: slow_sb_pool_spinup_accelerator_start_year         !< Year in which spin-up accelerator is
                                                                  !< executed for the first time
    INTEGER :: slow_sb_pool_spinup_accelerator_max_executions     !< Maximum number of spin-up accelerator runs
    LOGICAL                       :: init_from_ifs                !< Initialize ICON-Land from IFS data
    CHARACTER(len=filename_max)   :: ifs_filename                 !< File containing IFS data
    LOGICAL                       :: is_single_site
    CHARACTER(len=filename_max)   :: single_site_name
    CHARACTER(len=SHORT_NAME_LEN) :: single_site_shortname
    REAL(lp)                      :: single_site_lonlat(2)

    ! Definition of the namelist
    NAMELIST /jsb_model_nml/  &
      model_scheme, &
      use_tmx, &
      use_lakes, &
      use_glacier, &
      check_child_cf, &
      enforce_water_budget, &
      quincy_model, &
      usecase, &
      npft, &
      tpe_scheme, &
      fract_filename, &
      output_tiles, &
      relative_fractions_in_file, &
      force_c_isotopes, &
      include_nitrogen, &
      include_phosphorus, &
      include_carbon13, &
      include_carbon14, &
      include_nitrogen15, &
      l_spinup, &
      flag_slow_sb_pool_spinup_accelerator, &
      slow_sb_pool_spinup_accelerator_frequency, &
      slow_sb_pool_spinup_accelerator_length, &
      slow_sb_pool_spinup_accelerator_start_year, &
      slow_sb_pool_spinup_accelerator_max_executions, &
      init_from_ifs, &
      ifs_filename, &
      is_single_site, &
      single_site_name, &
      single_site_shortname, &
      single_site_lonlat


    INTEGER ::       &
      & nml_handler, &     !< File handler for the namelist
      & nml_unit,    &     !< Unit number of the namelist file
      & istat,       &     !< Status reading the namelist
      & i                  !< Looping index

    CHARACTER(len=*), PARAMETER :: routine = modname//':new_model_config' !< Name of the module

    CALL message(TRIM(routine), 'Starting model configuration from '//TRIM(namelist_filename))

    ALLOCATE(model_config)

    ! TODO Dissolve the dependency<br>
    ! `model_scheme` values are
    ! defined here as `INTEGER` replacing the ENUM `USE mo_jsb_model_class,    ONLY: MODEL_JSBACH, MODEL_QUINCY`
    ! This is needed to avoid a dependency cycle, creating a compile error (tested with NAG):
    !  stating that in `mo_jsb_memory_class` the `mo_jsb_config_class` is not available for `USE mo_jsb_config_class,  ONLY: t_jsb_config`
    MODEL_JSBACH               = 1
    MODEL_QUINCY               = 2

    ! Set defaults
    model_scheme               = "jsbach"
    use_tmx                    = .FALSE.
    use_lakes                  = .TRUE.
    use_glacier                = .TRUE.
    check_child_cf             = .FALSE.
    enforce_water_budget       = 'ignore'
    quincy_model               = 'land'
    qmodel_id                  = QLAND
    grid_filename              = ''
    grid_name                  = ''
    usecase                    = ''
    npft                       = 0
    tpe_scheme                 = ''
    fract_filename             = 'bc_land_frac.nc'
    output_tiles(:)            = ''
    relative_fractions_in_file = .TRUE.
    force_c_isotopes           = .TRUE.
    include_nitrogen           = .TRUE.
    include_phosphorus         = .TRUE.
    include_carbon13           = .TRUE.
    include_carbon14           = .TRUE.
    include_nitrogen15         = .TRUE.
    l_spinup                   = .FALSE.
    flag_slow_sb_pool_spinup_accelerator    = .FALSE.
    slow_sb_pool_spinup_accelerator_frequency = 100
    slow_sb_pool_spinup_accelerator_length = 1000
    slow_sb_pool_spinup_accelerator_start_year = 300
    slow_sb_pool_spinup_accelerator_max_executions = 4
    init_from_ifs              = .FALSE.
    ifs_filename               = 'ifs2icon.nc'
    is_single_site             = .FALSE.
    single_site_name           = ''
    single_site_shortname      = ''
    single_site_lonlat         = [0._lp, 0._lp]

    ! Read namelist
    nml_handler = open_nml(TRIM(namelist_filename))

    ! jsb_grid_nml
    nml_unit = position_nml('jsb_grid_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_grid_nml)

    ! Write resulting values in config memory
    model_config%grid_filename = grid_filename
    model_config%grid_name     = grid_name

    ! jsb_model_nml
    nml_unit =  position_nml('jsb_model_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_model_nml)
    CALL close_nml(nml_handler)

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
    !>
    !> @TODO To be deleted: `relative_fractions_in_file` will be removed (clean up task)
    !>
    ! No absolute fraction files are available, therefore, assert accidental use of flag:
    IF(.NOT. relative_fractions_in_file) THEN
      CALL finish(TRIM(routine), 'relative_fractions_in_file was set to true, however, currently' &
      & // ' no absolute fraction files are available, please check.')
    ENDIF

    IF (model_config%model_scheme /= MODEL_JSBACH .AND. use_tmx) THEN
      CALL finish(routine, 'tmx currently only implemented for JSBACH.')
    END IF

    model_config%relative_fractions_in_file = relative_fractions_in_file
    ! set other namelist options
    model_config%use_tmx                    = use_tmx
    model_config%use_lakes                  = use_lakes
    model_config%use_glacier                = use_glacier
    model_config%check_child_cf             = check_child_cf
    model_config%quincy_model               = TRIM(quincy_model)
    model_config%qmodel_id                  = qmodel_id
    model_config%usecase                    = usecase
    model_config%npft                       = npft
    model_config%tpe_scheme                 = tpe_scheme
    model_config%fract_filename             = fract_filename
    model_config%force_c_isotopes           = force_c_isotopes
    ! organise bgcm elements
    model_config%include_nitrogen           = include_nitrogen
    model_config%include_phosphorus         = include_phosphorus
    model_config%include_carbon13           = include_carbon13
    model_config%include_carbon14           = include_carbon14
    model_config%include_nitrogen15         = include_nitrogen15
    model_config%l_spinup                   = l_spinup
    ! The logicals are used to set up a bookkeeping list for the elements
    CALL setup_element_index_map_in_config(model_config, routine)
    model_config%flag_slow_sb_pool_spinup_accelerator   = flag_slow_sb_pool_spinup_accelerator
    model_config%slow_sb_pool_spinup_accelerator_frequency = slow_sb_pool_spinup_accelerator_frequency
    model_config%slow_sb_pool_spinup_accelerator_length = slow_sb_pool_spinup_accelerator_length
    model_config%slow_sb_pool_spinup_accelerator_start_year = slow_sb_pool_spinup_accelerator_start_year
    model_config%slow_sb_pool_spinup_accelerator_max_executions = slow_sb_pool_spinup_accelerator_max_executions

    !>
    !> @TODO To be deleted `relative_fractions_in_file` (clean up task).
    !>
    ! Currently, no absolute fraction files are available. Assert, the flag was not set accidentally.
    IF(.NOT. relative_fractions_in_file) THEN
      CALL finish(TRIM(routine), 'relative_fractions_in_file was set to false, however, currently' &
      & // ' no absolute fraction files are available, please check your settings.')
    ENDIF

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

    model_config%qmodel_id = Get_quincy_model_config_id(model_config%quincy_model)  ! QLAND, QPLANT, ...

    !>
    !> `model_config%init_from_ifs`, and `model_config%ifs_filename`.
    !>
    model_config%init_from_ifs              = init_from_ifs
    model_config%ifs_filename               = TRIM(ifs_filename)
    model_config%is_single_site             = is_single_site
    model_config%single_site_name           = TRIM(single_site_name)
    model_config%single_site_shortname      = TRIM(single_site_shortname)
    model_config%single_site_lonlat(1:2)    = single_site_lonlat(1:2)

    IF (is_single_site .AND. .NOT. jsbach_runs_standalone()) &
        & CALL finish(routine, 'Single site only allowed in standalone experiment')

    ! Get the number of valid entries in the vector
    !>
    !> @TODO use `SIZE(output_tiles)` instead of `99` to explore the
    !> `IF (output_tiles(i) == '') EXIT` cases.
    !>
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
  !> The logicals that indicate which bgcm elements are to be included are used to set up
  !> a bookkeeping list for the elements and to determine the number of elements actually in use.
  !> The bookkeeping list contains the IDs of the elements in use.
  !>
  !> **Note**: It is assumed that carbon is always used (ID = index = 1).
  !>
  SUBROUTINE setup_element_index_map_in_config(config, caller)
    USE mo_lnd_bgcm_class,        ONLY: ELEM_C_ID, ELEM_N_ID, ELEM_P_ID,  &
      &                                 ELEM_C13_ID, ELEM_C14_ID, ELEM_N15_ID, LAST_ELEM_ID
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model_config), INTENT(INOUT) :: config  !< Model config to be setup for used elements
    CHARACTER(len=*),            INTENT(IN) :: caller  !< Calling routine
    ! -------------------------------------------------------------------------------------------------- !
    ! allocate map with size "number of element IDs defined in mo_lnd_bgcm_class"
    config%nr_of_elements = LAST_ELEM_ID
    ALLOCATE(config%elements_index_map(config%nr_of_elements))
    config%elements_index_map(:) = -1  ! set all indices (IND) to '-1' by default, -1 for "element not used"
    config%nr_of_used_elements = 0     ! init for use in setup_element
    ! allocate vector with info about whether the element is used
    ALLOCATE(config%is_element_used(config%nr_of_elements))
    config%is_element_used(:) = .FALSE. ! set all elements to "not used" by default
    CALL message(TRIM(caller), 'Mapping of element ID -> IND:')

    ! carbon is presumed to always be used, with ID = IND = 1 (thus call first for carbon!)
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
  !> Small helper to add a given element ID to the bookkeeping list - according to the given (include_) flag
  !>
  SUBROUTINE setup_element(config, elem_flag, elem_id, short_name, caller)
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model_config), INTENT(INOUT) :: config       !< Model config to be setup for used elements
    LOGICAL,                     INTENT(IN) :: elem_flag    !< Flag for element to be setup
    INTEGER,                     INTENT(IN) :: elem_id      !< ID of element to be setup
    CHARACTER(len=*),            INTENT(IN) :: short_name   !< A short name for the element
    CHARACTER(len=*),            INTENT(IN) :: caller       !< Calling routine
    ! -------------------------------------------------------------------------------------------------- !
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
