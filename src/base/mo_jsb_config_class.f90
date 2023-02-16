!>
!! @brief Contains structures and methods for JSBACH model config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! code additions for bgcm handling Julia Nabel, MPI BGC (20230208)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_config_class
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN
  USE mo_io_units,           ONLY: filename_max
  USE mo_exception,          ONLY: message, finish, message_text
  USE mo_jsb_io_netcdf,      ONLY: t_input_file
  ! USE mo_jsb_model_class,    ONLY: MODEL_JSBACH, MODEL_JSBACH_QUINCY, MODEL_QUINCY   ! comment out to avoid dependency cycle

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_config, t_jsb_config_p, t_jsb_model_config, new_model_config

  ! Type for model configuration
  TYPE t_jsb_model_config
    CHARACTER(len=filename_max)    :: grid_filename
    CHARACTER(len=filename_max)    :: fract_filename
    CHARACTER(len=30)              :: grid_name
    CHARACTER(len=40)              :: usecase
    CHARACTER(len=20)              :: model_scheme_char   !< model_scheme from namelist (used in mo_jsb_lctlib_class)
    INTEGER                        :: model_scheme        !< model_scheme as integer
    LOGICAL                        :: use_lakes
    LOGICAL                        :: use_glacier
    CHARACTER(len=10)              :: quincy_model        !< canopy/plant/land - selection of Tasks and Code within subroutines
    CHARACTER(len=SHORT_NAME_LEN)  :: tpe_scheme          !< For terraplanet setup: open/closed
    CHARACTER(len=SHORT_NAME_LEN)  :: hsm_mode
    LOGICAL                        :: l_compat401         !< Use configuration compatible with jsbach 4.01
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: output_tiles(:)  ! List of tile names for which output should be generated
    LOGICAL                        :: relative_fractions_in_file  ! true: tile fractions in fract_filename are relative to parent tile
    ! Characteristics of bgc material in this model: flags specifying used elements and bookkeeping vector for elements used
    ! Note: no logical for carbon because we presuppose that it is always used!
    LOGICAL                        :: include_nitrogen      !< If bgc materials contain nitrogen as element
    LOGICAL                        :: include_phosphorus    !< If bgc materials contain phosphorus as element
    LOGICAL                        :: include_carbon13      !< If bgc materials contain C13 as element
    LOGICAL                        :: include_carbon14      !< If bgc materials contain C14 as element
    LOGICAL                        :: include_nitrogen15    !< If bgc materials contain C15 as element
    INTEGER                        :: nr_of_used_elements   !< for bookkeeping: number of actually used elements
    INTEGER, ALLOCATABLE           :: elements_index_map(:) !< vector for ID -> IDX mapping of the elements

    ! quincy @TODO temporary solution for model%config%options
    REAL(wp)          :: dtime                   !< length of one time step in seconds; defaults to 1800  @TODO
    LOGICAL           :: flag_stand_harvest
    INTEGER           :: stand_replacing_year
    LOGICAL           :: l_do_stand_replacing_harvest
    LOGICAL           :: l_transient_spinup
    ! quincy - end
    LOGICAL                        :: init_from_ifs    ! Initialize from IFS analysis
    CHARACTER(len=:), ALLOCATABLE  :: ifs_filename
    TYPE(t_input_file)             :: ifs_input_file

  END TYPE t_jsb_model_config

  ! Abstract type for process configurations
  TYPE, ABSTRACT :: t_jsb_config
    LOGICAL :: active
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

  ABSTRACT INTERFACE
    SUBROUTINE Init_config(config)
      IMPORT :: t_jsb_config
      CLASS(t_jsb_config), INTENT(inout) :: config
    END SUBROUTINE Init_config
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_config_class'

CONTAINS

  FUNCTION new_model_config(namelist_filename) RESULT(model_config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_util_string,        ONLY: tolower

    CHARACTER(len=*),         INTENT(in) :: namelist_filename
    TYPE(t_jsb_model_config), POINTER    :: model_config

    CHARACTER(len=filename_max) :: grid_filename
    CHARACTER(len=30)           :: grid_name
    INTEGER                     :: MODEL_JSBACH, MODEL_JSBACH_QUINCY, MODEL_QUINCY

    NAMELIST /jsb_grid_nml/   &
      grid_filename,          &
      grid_name

    CHARACTER(len=20)             :: model_scheme     !< model_scheme: jsbach / jsbach_quincy / quincy
    LOGICAL                       :: use_lakes
    LOGICAL                       :: use_glacier
    CHARACTER(len=10)             :: quincy_model
    LOGICAL                       :: l_compat401
    CHARACTER(len=filename_max)   :: fract_filename
    CHARACTER(len=40)             :: usecase
    CHARACTER(len=SHORT_NAME_LEN) :: tpe_scheme
    CHARACTER(len=SHORT_NAME_LEN) :: hsm_mode
    CHARACTER(len=SHORT_NAME_LEN) :: output_tiles(99)
    LOGICAL                       :: relative_fractions_in_file
    LOGICAL                       :: include_nitrogen
    LOGICAL                       :: include_phosphorus
    LOGICAL                       :: include_carbon13
    LOGICAL                       :: include_carbon14
    LOGICAL                       :: include_nitrogen15
    LOGICAL                       :: init_from_ifs
    CHARACTER(len=filename_max)   :: ifs_filename



    NAMELIST /jsb_model_nml/  &
      model_scheme, &
      use_lakes, &
      use_glacier, &
      quincy_model, &
      l_compat401, &
      usecase, &
      tpe_scheme, &
      hsm_mode, &
      fract_filename, &
      output_tiles, &
      relative_fractions_in_file, &
      include_nitrogen, &
      include_phosphorus, &
      include_carbon13, &
      include_carbon14, &
      include_nitrogen15, &
      init_from_ifs, &
      ifs_filename


    INTEGER :: nml_handler, nml_unit, istat, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':new_model_config'

    CALL message(TRIM(routine), 'Starting model configuration from '//TRIM(namelist_filename))

    ALLOCATE(model_config)

    ! model_scheme values
    ! define here as INTEGER replacing the ENUM "USE mo_jsb_model_class,    ONLY: MODEL_JSBACH, MODEL_JSBACH_QUINCY, MODEL_QUINCY"
    ! this is needed to avoid a dependency cycle, creating an compile error (tested with NAG):
    !  stating that in mo_jsb_memory_class the mo_jsb_config_class is not available for "USE mo_jsb_config_class,  ONLY: t_jsb_config"
    MODEL_JSBACH               = 1
    MODEL_JSBACH_QUINCY        = 2
    MODEL_QUINCY               = 3

    ! Set defaults
    model_scheme               = "jsbach"
    use_lakes                  = .TRUE.
    use_glacier                = .TRUE.
    quincy_model               = 'land'
    l_compat401                = .FALSE.
    grid_filename              = ''
    grid_name                  = ''
    usecase                    = ''
    tpe_scheme                 = ''
    hsm_mode                   = 'simple'
    fract_filename             = 'bc_land_frac.nc'
    output_tiles(:)            = ''
    relative_fractions_in_file = .TRUE.
    include_nitrogen           = .TRUE.
    include_phosphorus         = .TRUE.
    include_carbon13           = .TRUE.
    include_carbon14           = .TRUE.
    include_nitrogen15         = .TRUE.
    init_from_ifs              = .FALSE.
    ifs_filename               = 'ifs2icon.nc'

    ! Read namelist
    nml_handler = open_nml(TRIM(namelist_filename))

    nml_unit = position_nml('jsb_grid_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_grid_nml)

    ! Write resulting values in config memory
    model_config%grid_filename = grid_filename
    model_config%grid_name     = grid_name

    nml_unit =  position_nml('jsb_model_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_model_nml)
    CALL close_nml(nml_handler)

    model_config%model_scheme_char = tolower(model_scheme)
    SELECT CASE (TRIM(model_config%model_scheme_char))
    CASE ("jsbach")
      model_config%model_scheme = MODEL_JSBACH
      CALL message(TRIM(routine), 'Use model_scheme JSBACH')
    CASE ("jsbach_quincy")
      model_config%model_scheme = MODEL_JSBACH_QUINCY
      CALL message(TRIM(routine), 'Use model_scheme JSBACH with QUINCY')
      CALL finish(TRIM(routine),  'The model_scheme JSBACH with QUINCY (MODEL_JSBACH_QUINCY) is not yet implemented.')
    CASE ("quincy")
      model_config%model_scheme = MODEL_QUINCY
      CALL message(TRIM(routine), 'Use model_scheme QUINCY')
    CASE DEFAULT
      CALL finish(TRIM(routine), 'Selected model_scheme '//TRIM(model_scheme)//' not available.')
    END SELECT
    ! currently, no absolute fraction files are available, therefore, assert accidental use of flag:
    IF(.NOT. relative_fractions_in_file) THEN
      CALL finish(TRIM(routine), 'relative_fractions_in_file was set to true, however, currently' &
      & // ' no absolute fraction files are available, please check.')
    ENDIF
    model_config%relative_fractions_in_file = relative_fractions_in_file
    ! set other namelist options
    model_config%use_lakes                  = use_lakes
    model_config%use_glacier                = use_glacier
    model_config%quincy_model               = TRIM(quincy_model)
    model_config%l_compat401                = l_compat401
    model_config%usecase                    = usecase
    model_config%tpe_scheme                 = tpe_scheme
    model_config%hsm_mode                   = hsm_mode
    model_config%fract_filename             = fract_filename
    ! organise elements
    model_config%include_nitrogen           = include_nitrogen
    model_config%include_phosphorus         = include_phosphorus
    model_config%include_carbon13           = include_carbon13
    model_config%include_carbon14           = include_carbon14
    model_config%include_nitrogen15         = include_nitrogen15
    CALL setup_element_index_map_in_config(model_config, routine)

    ! quincy - temporary solution @TODO
    model_config%dtime                        = 1800.0_wp   ! length of one time step in seconds
    model_config%flag_stand_harvest           = .FALSE.
    model_config%stand_replacing_year         = 1500
    model_config%l_do_stand_replacing_harvest = .FALSE.
    model_config%l_transient_spinup           = .FALSE.
    ! quincy - end
    model_config%init_from_ifs              = init_from_ifs
    model_config%ifs_filename               = TRIM(ifs_filename)

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
    ALLOCATE(config%elements_index_map(LAST_ELEM_ID - 1))
    config%elements_index_map(:) = -1  ! set all IND to '-1' by default, -1 for "element not used"
    config%nr_of_used_elements = 0     ! init for use in setup_element
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
    IF (elem_flag) THEN
      config%nr_of_used_elements = config%nr_of_used_elements + 1
      config%elements_index_map(elem_id) = config%nr_of_used_elements
    END IF

    WRITE (message_text,*) elem_id, " -> ", config%elements_index_map(elem_id)
    CALL message(TRIM(short_name), TRIM(message_text))

  END SUBROUTINE setup_element


#endif
END MODULE mo_jsb_config_class
