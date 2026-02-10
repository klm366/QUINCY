!> QUINCY module to read and convert external forcing data
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!> - expects annual files with year in filename
!> - current implementation only for daily, timestep or constant data
!>
MODULE mo_qs_forcing
#ifndef __NO_JSBACH__

  USE mo_kind

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Init_forcing_qs, calc_output_time_interval_pool
  PUBLIC :: forcing_options

  !-----------------------------------------------------------------------------------------------------
  !> Structures for handling forcing/input data
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_forcing_variables
    REAL(wp)          :: tair   ! timestep forcing
    REAL(wp)          :: tmin   ! daily forcing
    REAL(wp)          :: tmax   ! daily forcing
    REAL(wp)          :: qair
    REAL(wp)          :: sw
    REAL(wp)          :: lw
    REAL(wp)          :: rain
    REAL(wp)          :: snow
    REAL(wp)          :: precipitation
    REAL(wp)          :: pressure
    REAL(wp)          :: wind
  END TYPE t_forcing_variables

  TYPE t_forcing_options
    ! namelist
    CHARACTER(15)     :: n_deposition_scheme
    CHARACTER(15)     :: p_deposition_scheme
    CHARACTER(15)     :: forcing_mode
    LOGICAL           :: flag_read_dC13
    LOGICAL           :: flag_read_DC14
    LOGICAL           :: flag_forcing_with_snow_data
    CHARACTER(1)      :: simulation_length_unit
    INTEGER           :: simulation_length_number
    INTEGER           :: transient_spinup_start_year
    INTEGER           :: transient_spinup_end_year
    INTEGER           :: transient_spinup_years
    INTEGER           :: transient_simulation_start_year
    LOGICAL           :: flag_forcing_co2_const
    LOGICAL           :: is_daily_forcing   ! default = .FALSE.; i.e., default is timestep forcing
    LOGICAL           :: read_precipitation ! default = .FALSE.; i.e., default is to read snow and rain separately | can be .TRUE. only if daily forcing is enabled

    ! calculated
    INTEGER           :: transient_spinup_timesteps_all     !< number of timesteps for spinup with forcing_mode=='transient': calculated from transient_spinup_years
    INTEGER           :: transient_spinup_timesteps_repeat  !< number of timesteps of the forcing period repeatedly simulated for spinup

    ! read from forcing files
    INTEGER           :: forcing_year
    INTEGER           :: forcing_doy
    REAL(wp)          :: forcing_hour       ! e.g. 0, 0.5, 1.0, 1.5, ... 23.5  for 30 min timesteps
    REAL(wp)          :: last_timestep_tair
    REAL(wp)          :: co2_mixing_ratio
    REAL(wp)          :: co2_dC13
    REAL(wp)          :: co2_dC14
    REAL(wp)          :: nhx_srf_down
    REAL(wp)          :: noy_srf_down
    REAL(wp)          :: p_srf_down

    TYPE(t_forcing_variables) :: current
    TYPE(t_forcing_variables) :: yDay

  CONTAINS
    PROCEDURE :: Timesteps => Get_timstep_number
  END TYPE t_forcing_options
  TYPE(t_forcing_options), SAVE :: forcing_options




  CHARACTER(len=*), PARAMETER :: modname = 'mo_qs_forcing'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_sb_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_forcing_qs(model_id)

    USE mo_jsb_namelist_iface,              ONLY: open_nml, POSITIONED, position_nml, close_nml
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
    CHARACTER(len=*), PARAMETER         :: routine = TRIM(modname)//':Init_forcing_qs'
    CHARACTER(15)                       :: n_deposition_scheme
    CHARACTER(15)                       :: p_deposition_scheme
    CHARACTER(15)                       :: forcing_mode                 !< static / transient
    LOGICAL                             :: flag_read_dC13
    LOGICAL                             :: flag_read_DC14
    LOGICAL                             :: flag_forcing_with_snow_data  !< forcing data includes data about snowfall
    CHARACTER(1)                        :: simulation_length_unit       !< d days, w weeks, y years
    INTEGER                             :: simulation_length_number
    INTEGER                             :: transient_spinup_start_year
    INTEGER                             :: transient_spinup_end_year
    INTEGER                             :: transient_spinup_years
    INTEGER                             :: transient_simulation_start_year
    LOGICAL                             :: flag_forcing_co2_const
    LOGICAL                             :: is_daily_forcing
    LOGICAL                             :: read_precipitation

    NAMELIST /jsb_forcing_ctl/                          & 
                      n_deposition_scheme             , &
                      p_deposition_scheme             , &
                      forcing_mode                    , &
                      flag_read_dC13                  , &
                      flag_read_DC14                  , &
                      flag_forcing_with_snow_data     , &
                      simulation_length_unit          , &
                      simulation_length_number        , &
                      transient_spinup_start_year     , &
                      transient_spinup_end_year       , &
                      transient_spinup_years          , &
                      transient_simulation_start_year , &
                      flag_forcing_co2_const          , &
                      is_daily_forcing                , &
                      read_precipitation


    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting Init_forcing_qs')

    ! get model and namelist_filename
    model => Get_model(model_id)

    !> 1.0 Read namelist
    !! 

    ! Set defaults
    n_deposition_scheme                 = "dynamic"
    p_deposition_scheme                 = "dynamic"
    forcing_mode                        = "static"    ! static / transient
    flag_read_dC13                      = .TRUE.
    flag_read_DC14                      = .TRUE.
    flag_forcing_with_snow_data         = .TRUE.
    simulation_length_unit              = 'y'         ! y (years), w (weeks), d (days)
    simulation_length_number            = 30 
    transient_spinup_start_year         = 1901
    transient_spinup_end_year           = 1930
    transient_spinup_years              = 500
    transient_simulation_start_year     = 1901
    flag_forcing_co2_const              = .FALSE.
    is_daily_forcing                    = .FALSE.
    read_precipitation                  = .FALSE.

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
    forcing_options%flag_forcing_with_snow_data       = flag_forcing_with_snow_data
    forcing_options%simulation_length_unit            = simulation_length_unit
    forcing_options%simulation_length_number          = simulation_length_number
    forcing_options%transient_spinup_start_year       = transient_spinup_start_year
    forcing_options%transient_spinup_end_year         = transient_spinup_end_year
    forcing_options%transient_spinup_years            = transient_spinup_years
    forcing_options%transient_simulation_start_year   = transient_simulation_start_year
    forcing_options%flag_forcing_co2_const            = flag_forcing_co2_const
    forcing_options%is_daily_forcing                  = is_daily_forcing
    forcing_options%read_precipitation                = read_precipitation

    !> 2.0 Modify .. Configure ..
    !! 




    !> 3.0 Init forcing input 
    !! 
    forcing_options%forcing_year          = 0
    forcing_options%forcing_doy           = 0
    forcing_options%forcing_hour          = 0.0_wp

    forcing_options%yDay%tair             = -999.0_wp
    forcing_options%last_timestep_tair    = -999.0_wp

    forcing_options%current%sw            = -999.0_wp
    forcing_options%current%lw            = -999.0_wp
    forcing_options%current%tmin          = -999.0_wp
    forcing_options%current%tmax          = -999.0_wp
    forcing_options%current%qair          = -999.0_wp
    forcing_options%current%wind          = -999.0_wp
    forcing_options%current%rain          = -999.0_wp
    forcing_options%current%snow          = -999.0_wp
    forcing_options%current%precipitation = -999.0_wp
    forcing_options%current%pressure      = -999.0_wp

    !> 4.0 output info about settings to logfile
    !>
    ! convert rain to snow is enabled
    IF (.NOT. flag_forcing_with_snow_data) THEN
      CALL message(TRIM(routine), 'Conversion of rain to snow is enabled (flag_forcing_with_snow_data == TRUE).')
    END IF

  END SUBROUTINE Init_forcing_qs


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
  FUNCTION calc_output_time_interval_pool(output_interval_pool, dtime) RESULT(interval_number)

    USE mo_jsb_math_constants,      ONLY: one_day

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CHARACTER(len=*)  ,INTENT(in) :: output_interval_pool
    REAL(wp)          ,INTENT(in) :: dtime                  !< timestep length
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

#endif
END MODULE mo_qs_forcing
