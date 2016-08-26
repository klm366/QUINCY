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
    REAL(wp)          :: dt        !< length of one time step in seconds; base_ctl namelist: read_dt_time_step_length_sec
    CHARACTER(1)      :: simulation_length_unit
    INTEGER           :: simulation_length_number
    CHARACTER(len=20) :: slm_model_name
    LOGICAL           :: fluxnet_output
    INTEGER           :: fluxnet_output_length
    REAL(wp)          :: lat, lon  !< latitude longitude currently stored in mo_base_memory:t_slm_grid
  CONTAINS
    PROCEDURE :: Init => Init_base_config
    PROCEDURE :: Timesteps => Get_timstep_number
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
    !LOGICAL          :: active
    REAL(wp)          :: dt_time_step_length_sec   ! model-variable: dt 
    CHARACTER(1)      :: simulation_length_unit
    INTEGER           :: simulation_length_number
    CHARACTER(len=20) :: slm_model_name
    LOGICAL           :: fluxnet_output
    INTEGER           :: fluxnet_output_length
    REAL(wp)          :: latitude, longitude


    NAMELIST /base_ctl/            &  
        dt_time_step_length_sec,   &
        simulation_length_unit,    &
        simulation_length_number,  &
        slm_model_name,            &
        fluxnet_output,            &
        fluxnet_output_length

    NAMELIST /grid_ctl/            &  
        latitude,                  &
        longitude


    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_base_config'

    CALL message(TRIM(routine), 'Starting base configuration')

    ! Set defaults
    dt_time_step_length_sec       = 1800.0_wp
    simulation_length_unit        = 'y'         ! y (years), w (weeks), d (days)
    simulation_length_number      = 20 
    slm_model_name                = 'land'
    fluxnet_output                = .FALSE.     ! .TRUE. .FALSE.
    fluxnet_output_length         = 2
    latitude                      = 50.91       ! about the MPI-BGC Jena
    longitude                     = 11.56       ! about the MPI-BGC Jena


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
    config%dt                         = dt_time_step_length_sec
    config%simulation_length_unit     = simulation_length_unit
    config%simulation_length_number   = simulation_length_number
    config%slm_model_name             = slm_model_name
    config%fluxnet_output             = fluxnet_output
    config%fluxnet_output_length      = fluxnet_output_length
    config%lat                        = latitude
    config%lon                        = longitude
  END SUBROUTINE Init_base_config

  !-----------------------------------------------------------------------------------------------------
  !> calculates the number of timesteps the model may run based on config-options set in the namelist
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Get_timstep_number(this) RESULT(timesteps)

    USE mo_base_constants, ONLY: one_day, one_year  ! length of one year in days, length of one day in seconds

    IMPLICIT NONE 
    ! 0.1 InOut
    CLASS(t_base_config), INTENT(in) :: this        !< config type for base
    ! 0.2 Local
    INTEGER                          :: timesteps    

    SELECT CASE(this%simulation_length_unit)
      CASE('y') ! years
        timesteps = NINT(this%simulation_length_number * one_year   * one_day / this%dt)
      CASE('w') ! weeks
        timesteps = NINT(this%simulation_length_number * 7.0_wp     * one_day / this%dt)
      CASE('d') ! days
        timesteps = NINT(this%simulation_length_number              * one_day / this%dt)
      CASE DEFAULT
        timesteps = 1
    END SELECT 
  END FUNCTION Get_timstep_number


!#endif
END MODULE mo_base_config_class
