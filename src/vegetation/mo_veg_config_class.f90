!>
!! @brief Contains structures and methods for veg config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04);
!! replaces the previous mo_veg_config_class created 2016-03-17 that was based on the JSBACH4 mo_veg_config_class anyway;
!! integration into SLM by Jan Engel, MPI BGC (20160401)
!!
MODULE mo_veg_config_class
!#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish
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
    CHARACTER(8)     :: output_interval_pool, output_interval_flux
    INTEGER          :: output_start_first_day_year, output_end_last_day_year !< output: 1st to last day of the year(s) provided
    INTEGER          :: pft_id
    LOGICAL          :: flag_optimal_Nfraction, flag_optimal_Nleaf, flag_calc_soil_uptake, flag_calc_soil_moisture
    REAL(wp)         :: soil_awc_prescribe, soil_theta_prescribe
    CONTAINS
      PROCEDURE :: Init => Init_veg_config
  END Type t_veg_config

  TYPE(t_veg_config), SAVE :: veg_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_config_class'

CONTAINS

  SUBROUTINE Init_veg_config(config)   ! 'config' == t_veg_config, see ABSTRACT INTERFACE -> SUBROUTINE Init_config(config)
!    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
!    using the old JSBACH3 functions instead of wrappers implemented in JSBACH4 -> mo_jsb_namelist_iface
    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file

    CLASS(t_veg_config), INTENT(inout) :: config

    ! variables for reading from namlist, identical to variable-name in namelist
    !LOGICAL      :: active
    CHARACTER(8) :: output_interval_pool, output_interval_flux
    INTEGER      :: output_start_first_day_year, output_end_last_day_year, plant_functional_type_id
    LOGICAL      :: flag_optimal_Nfraction, flag_optimal_Nleaf, flag_calc_soil_uptake, flag_calc_soil_moisture
    REAL(wp)     :: soil_awc_prescribe, soil_theta_prescribe

    NAMELIST /vegetation_ctl/       & 
      output_interval_pool,         &
      output_interval_flux,         &
      output_start_first_day_year,  &
      output_end_last_day_year,     &
      flag_optimal_Nfraction,       &
      flag_optimal_Nleaf,           &
      flag_calc_soil_uptake,        &
      flag_calc_soil_moisture,      &
      soil_awc_prescribe,           &
      soil_theta_prescribe,         &
      plant_functional_type_id

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_veg_config'

    CALL message(TRIM(routine), 'Starting veg configuration')

    ! Set defaults
    output_interval_pool           = 'monthly' 
    output_interval_flux           = 'monthly' 
    output_start_first_day_year    = 1
    output_end_last_day_year       = 10
    flag_optimal_Nfraction         = .FALSE.  ! .TRUE. .FALSE.
    flag_optimal_Nleaf             = .FALSE.  ! .TRUE. .FALSE.
    flag_calc_soil_uptake          = .FALSE.  ! .TRUE. .FALSE.
    flag_calc_soil_moisture        = .FALSE.  ! .TRUE. .FALSE.
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

    ! pass values as read from file
    config%output_interval_pool          = output_interval_pool
    config%output_interval_flux          = output_interval_flux
    config%output_start_first_day_year   = output_start_first_day_year
    config%output_end_last_day_year      = output_end_last_day_year
    config%flag_optimal_Nfraction        = flag_optimal_Nfraction
    config%flag_optimal_Nleaf            = flag_optimal_Nleaf
    config%flag_calc_soil_uptake         = flag_calc_soil_uptake
    config%flag_calc_soil_moisture       = flag_calc_soil_moisture
    config%soil_awc_prescribe            = soil_awc_prescribe 
    config%soil_theta_prescribe          = soil_theta_prescribe 
    config%pft_id                        = plant_functional_type_id

  END SUBROUTINE Init_veg_config

!#endif
END MODULE mo_veg_config_class
