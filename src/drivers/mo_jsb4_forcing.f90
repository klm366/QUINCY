!>
!! @brief Module to init, read and convert external forcing data
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2017-06-30): based on the git:JSBACH4 mo_jsb4_forcing module integrating only forcing_options and a basic version of Init_forcing()
!! 

MODULE mo_jsb4_forcing

#include "dsl_jsbach.inc"

  USE mo_kind

  IMPLICIT NONE


  PRIVATE
  PUBLIC :: Init_forcing, calc_output_time_interval_pool
  PUBLIC :: forcing_options


  !-----------------------------------------------------------------------------------------------------
  !> Structure for handling forcing/input data
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_forcing_options
    ! namelist
    CHARACTER(15)     :: n_deposition_scheme
    CHARACTER(15)     :: p_deposition_scheme
    CHARACTER(15)     :: forcing_mode
    LOGICAL           :: flag_read_dC13
    LOGICAL           :: flag_read_DC14
    CHARACTER(1)      :: simulation_length_unit
    INTEGER           :: simulation_length_number
    INTEGER           :: transient_spinup_start_year
    INTEGER           :: transient_spinup_end_year
    INTEGER           :: transient_spinup_years
    INTEGER           :: transient_simulation_start_year
    LOGICAL           :: flag_transient_sim_co2_const

    ! calculated
    INTEGER           :: transient_spinup_timesteps_all     !< number of timesteps for spinup with forcing_mode=='transient': calculated from transient_spinup_years
    INTEGER           :: transient_spinup_timesteps_repeat  !< number of timesteps of the forcing period repeatedly simulated for spinup

    ! saved at runtime
    REAL(wp)          :: co2_mixing_ratio_1st_timestep
    REAL(wp)          :: co2_dC13_1st_timestep
    REAL(wp)          :: co2_DC14_1st_timestep

    ! read from forcing files
    INTEGER           :: forcing_year
    INTEGER           :: forcing_doy
    REAL(wp)          :: forcing_hour       ! e.g. 0, 0.5, 1.0, 1.5, ... 23.5  for 30 min timesteps
  CONTAINS
    PROCEDURE :: Timesteps => Get_timstep_number
  END TYPE t_forcing_options
  TYPE(t_forcing_options), SAVE :: forcing_options




  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb4_forcing'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_sb_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_forcing(model_id)

    USE mo_namelist,                        ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file
    USE mo_exception,                       ONLY: message_text, message, finish
    USE mo_jsb_class,                       ONLY: Get_model
    USE mo_jsb_model_class,                 ONLY: t_jsb_model

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER, INTENT(in)                 :: model_id
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER          :: model
    CHARACTER(len=*), PARAMETER         :: routine = TRIM(modname)//':Init_forcing'
    CHARACTER(15)                       :: n_deposition_scheme
    CHARACTER(15)                       :: p_deposition_scheme
    CHARACTER(15)                       :: forcing_mode                 !< static / transient
    LOGICAL                             :: flag_read_dC13
    LOGICAL                             :: flag_read_DC14
    CHARACTER(1)                        :: simulation_length_unit       !< d days, w weeks, y years
    INTEGER                             :: simulation_length_number
    INTEGER                             :: transient_spinup_start_year
    INTEGER                             :: transient_spinup_end_year
    INTEGER                             :: transient_spinup_years
    INTEGER                             :: transient_simulation_start_year
    LOGICAL                             :: flag_transient_sim_co2_const

    NAMELIST /jsb_forcing_ctl/                          & 
                      n_deposition_scheme             , &
                      p_deposition_scheme             , &
                      forcing_mode                    , &
                      flag_read_dC13                  , &
                      flag_read_DC14                  , &
                      simulation_length_unit          , &
                      simulation_length_number        , &
                      transient_spinup_start_year     , &
                      transient_spinup_end_year       , &
                      transient_spinup_years          , &
                      transient_simulation_start_year , &
                      flag_transient_sim_co2_const


    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting Init_forcing')

    ! get model and namelist_filename
    model => Get_model(model_id)

    !> 1.0 Read namelist
    !! 

    ! Set defaults
    n_deposition_scheme                 = "none"
    p_deposition_scheme                 = "none"
    forcing_mode                        = "static"    ! static / transient
    flag_read_dC13                      = .FALSE.
    flag_read_DC14                      = .FALSE.
    simulation_length_unit              = 'y'         ! y (years), w (weeks), d (days)
    simulation_length_number            = 30 
    transient_spinup_start_year         = 1901
    transient_spinup_end_year           = 1930
    transient_spinup_years              = 500
    transient_simulation_start_year     = 1901
    flag_transient_sim_co2_const        = .FALSE.

    ! read the namelist
    ! 
    ! &JSB_FORCING_CTL
    nml_handler = open_nml(TRIM(model%namelist_filename))
    nml_unit = position_nml('JSB_FORCING_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_forcing_ctl)
    CALL close_nml(nml_handler)

    ! test some input values
    SELECT CASE (TRIM(n_deposition_scheme))
      CASE ("none","constant","dynamic")
        ! do nothing 
      CASE DEFAULT
        WRITE(message_text,'(2a)') 'invalid n_deposition_scheme: ',n_deposition_scheme
        CALL finish(routine,message_text)
    ENDSELECT
    SELECT CASE (TRIM(p_deposition_scheme))
      CASE ("none","constant","dynamic")
        ! do nothing 
      CASE DEFAULT
        WRITE(message_text,'(2a)') 'invalid p_deposition_scheme: ',p_deposition_scheme
        CALL finish(routine,message_text)
    ENDSELECT


    ! pass values as read from file
    forcing_options%n_deposition_scheme               = n_deposition_scheme
    forcing_options%p_deposition_scheme               = p_deposition_scheme
    forcing_options%forcing_mode                      = forcing_mode
    forcing_options%flag_read_dC13                    = flag_read_dC13
    forcing_options%flag_read_DC14                    = flag_read_DC14
    forcing_options%simulation_length_unit            = simulation_length_unit
    forcing_options%simulation_length_number          = simulation_length_number
    forcing_options%transient_spinup_start_year       = transient_spinup_start_year
    forcing_options%transient_spinup_end_year         = transient_spinup_end_year
    forcing_options%transient_spinup_years            = transient_spinup_years
    forcing_options%transient_simulation_start_year   = transient_simulation_start_year
    forcing_options%flag_transient_sim_co2_const      = flag_transient_sim_co2_const

    !> 2.0 Modify .. Configure ..
    !! 




    !> 3.0 Init forcing input 
    !! 
    forcing_options%forcing_year    = 0
    forcing_options%forcing_doy     = 0
    forcing_options%forcing_hour    = 0.0_wp


  END SUBROUTINE Init_forcing


  !-----------------------------------------------------------------------------------------------------
  !> calculates the number of timesteps the model may run based on config-options set in the namelist
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Get_timstep_number(this, dtime) RESULT(timesteps)

    USE mo_jsb_math_constants,              ONLY: one_day, one_year  ! length of one year in days, length of one day in seconds

    IMPLICIT NONE 
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_forcing_options), INTENT(inout)  :: this 
    REAL(wp)                , INTENT(in)     :: dtime
    ! ---------------------------
    ! 0.2 Local
    INTEGER                               :: timesteps   
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_timstep_number' 


    !< 1.0
    !! timesteps of the simulation length
    SELECT CASE(this%simulation_length_unit)
      CASE('y') ! years
        timesteps = NINT(this%simulation_length_number * one_year   * one_day / dtime)
      CASE('w') ! weeks
        timesteps = NINT(this%simulation_length_number * 7.0_wp     * one_day / dtime)
      CASE('d') ! days
        timesteps = NINT(this%simulation_length_number              * one_day / dtime)
      CASE DEFAULT
        timesteps = 1
    END SELECT 

    !< 2.0
    !! timing of transient spinup simulations, 
    !! namelist values are provided in unit: years
    this%transient_spinup_timesteps_all     = 1
    this%transient_spinup_timesteps_repeat  = 1
    IF ( TRIM(this%forcing_mode) == "transient" ) THEN
      this%transient_spinup_timesteps_all     = NINT(this%transient_spinup_years * one_year * one_day / dtime)
      this%transient_spinup_timesteps_repeat  = NINT( (this%transient_spinup_end_year - this%transient_spinup_start_year + 1 ) * &
                                                     one_year * one_day / dtime)
      
    END IF


  END FUNCTION Get_timstep_number


  !-----------------------------------------------------------------------------------------------------
  !> the input interval string is transformed into a number, which is used for timely model output
  !! 
  !! valid for pool output, but applied for pools and flux output
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_output_time_interval_pool(output_interval_pool) RESULT(interval_number)

    USE mo_jsb_math_constants,      ONLY: one_day, dtime

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CHARACTER(len=*)  ,INTENT(in) :: output_interval_pool
    REAL(wp)                      :: interval_number        !< number output



    !> 1.0 create interval_number depending on pool/flux
    !!
    SELECT CASE(output_interval_pool)
      CASE('none')
        interval_number = 0.0_wp
      CASE('timestep')
        interval_number = dtime/one_day  ! length of one time step as advanced in 'time' 
      CASE('daily')
        interval_number = 1.0_wp
      CASE('weekly')
        interval_number = 7.0_wp
      CASE('monthly')
        interval_number = 30.0_wp
      CASE('yearly')
        interval_number = 365.0_wp
    END SELECT
  END FUNCTION calc_output_time_interval_pool

  
END MODULE mo_jsb4_forcing
