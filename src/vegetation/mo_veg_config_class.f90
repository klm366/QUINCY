!>
!! @brief Contains structures and methods for veg config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04);
!! replaces the previous SLM mo_veg_config_class created 2016-03-17 that was based on the JSBACH4 mo_veg_config_class anyway;
!! integration into SLM by Jan Engel, MPI BGC (20160401)
!!
MODULE mo_veg_config_class
!#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
!  PUBLIC :: t_veg_config
  PUBLIC :: veg_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the vegetation process, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_veg_config
    CHARACTER(8)     :: output_interval_pool, &
                        output_interval_flux
    INTEGER          :: output_start_first_day_year, &  !< output start at the 1st  day of this year
                        output_end_last_day_year        !< output stop  at the last day of this year
    INTEGER          :: pft_id
    LOGICAL          :: flag_calc_soil_moisture    !< calculate simple water bucket model
    LOGICAL          :: flag_calc_soil_energy      !< calculate simple soil energy model 
    LOGICAL          :: flag_tile_veg_dynamics     !< calculate within-tile vegetation dynamics + start from bareground
    LOGICAL          :: flag_optimal_Nleaf         !< optimise leaf C:N 
    LOGICAL          :: flag_optimal_allocation    !< optimise plant allocation
    CHARACTER(15)    :: n_uptake_scheme            !< select scheme to simulate N uptake
                                                   ! "standard|resistance|optimal"           
    REAL(wp)         :: soil_awc_prescribe, &
                        soil_theta_prescribe
    CONTAINS
      PROCEDURE :: Init => Init_veg_config
  END Type t_veg_config

  TYPE(t_veg_config), SAVE :: veg_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_veg_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_veg_config(config)   ! 'config' == t_veg_config, see ABSTRACT INTERFACE -> SUBROUTINE Init_config(config)
!    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
!    using the old JSBACH3 functions instead of wrappers implemented in JSBACH4 -> mo_jsb_namelist_iface
    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file
    ! 0.1 InOut
    CLASS(t_veg_config), INTENT(inout) :: config    !< config type for veg
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    !LOGICAL      :: active
    CHARACTER(8)  :: output_interval_pool          , &
                     output_interval_flux
    INTEGER       :: output_start_first_day_year   , &
                     output_end_last_day_year      , &
                     plant_functional_type_id
    LOGICAL       :: flag_calc_soil_moisture
    LOGICAL       :: flag_calc_soil_energy
    LOGICAL       :: flag_tile_veg_dynamics
    LOGICAL       :: flag_optimal_Nleaf
    LOGICAL       :: flag_optimal_allocation
    CHARACTER(15) :: n_uptake_scheme
    REAL(wp)      :: soil_awc_prescribe             , &
                     soil_theta_prescribe
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_veg_config'

    NAMELIST /vegetation_ctl/       & 
      output_interval_pool,         &
      output_interval_flux,         &
      output_start_first_day_year,  &
      output_end_last_day_year,     &
      flag_calc_soil_moisture,      &
      flag_calc_soil_energy,        &
      flag_tile_veg_dynamics,       &
      flag_optimal_Nleaf,           &
      flag_optimal_allocation,      &
      n_uptake_scheme,              & 
      soil_awc_prescribe,           &
      soil_theta_prescribe,         &
      plant_functional_type_id

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting veg configuration')

    ! Set defaults
    output_interval_pool           = 'monthly' ! timestep daily weekly monthly yearly
    output_interval_flux           = 'monthly' ! timestep daily weekly monthly yearly
    output_start_first_day_year    = 1
    output_end_last_day_year       = 10
    flag_calc_soil_moisture        = .FALSE.  ! .TRUE. .FALSE.
    flag_calc_soil_energy          = .FALSE.  ! .TRUE. .FALSE.
    flag_tile_veg_dynamics         = .FALSE.  ! .TRUE. .FALSE.
    flag_optimal_Nleaf             = .FALSE.  ! .TRUE. .FALSE.
    flag_optimal_allocation        = .FALSE.  ! .TRUE. .FALSE.
    n_uptake_scheme                = "standard"  ! prescribed standard resistance optimal
    soil_awc_prescribe             = 300.0_wp
    soil_theta_prescribe           = 100.0_wp 
    plant_functional_type_id       = 1

    ! read the namelist
    ! 
    ! &VEGETATION_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('VEGETATION_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, vegetation_ctl)
    CALL close_nml(nml_handler)

    ! check input for n_uptake scheme
    SELECT CASE (TRIM(n_uptake_scheme))
     CASE ("prescribed","standard","resistance","optimal")
     CASE DEFAULT
       WRITE(message_text,'(2a)') 'invalid n_uptake_scheme',n_uptake_scheme
       CALL finish("mo_veg_config_class:init_veg_config",message_text) 
    ENDSELECT

    ! pass values as read from file
    config%output_interval_pool          = output_interval_pool
    config%output_interval_flux          = output_interval_flux
    config%output_start_first_day_year   = output_start_first_day_year
    config%output_end_last_day_year      = output_end_last_day_year
    config%flag_calc_soil_moisture       = flag_calc_soil_moisture
    config%flag_calc_soil_energy         = flag_calc_soil_energy
    config%flag_tile_veg_dynamics        = flag_tile_veg_dynamics  
    config%flag_optimal_Nleaf            = flag_optimal_Nleaf
    config%flag_optimal_allocation       = flag_optimal_allocation
    config%n_uptake_scheme               = TRIM(n_uptake_scheme)
    config%soil_awc_prescribe            = soil_awc_prescribe 
    config%soil_theta_prescribe          = soil_theta_prescribe 
    config%pft_id                        = plant_functional_type_id

  END SUBROUTINE Init_veg_config

!#endif
END MODULE mo_veg_config_class
