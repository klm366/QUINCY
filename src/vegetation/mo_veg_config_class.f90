!>
!! @brief Contains structures and methods for veg config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! replacing the previous mo_veg_config_class created 2016-03-17, which was based on the JSBACH4 mo_veg_config_class anyway
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

  TYPE, EXTENDS(t_jsb_config) :: t_veg_config
    CHARACTER(8)     :: output_interval_pool, output_interval_flux
    INTEGER          :: output_start_date, output_end_date
    LOGICAL          :: flag_optimal_Nfraction, flag_optimal_Nleaf, flag_calc_soil_uptake
    CHARACTER(8)     :: option_forcing
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

    LOGICAL      :: active
    CHARACTER(8) :: read_output_interval_pool, read_output_interval_flux
    INTEGER      :: output_start_date, output_end_date
    LOGICAL      :: flag_optimal_Nfraction, flag_optimal_Nleaf, flag_calc_soil_uptake
    CHARACTER(8) :: create_forcing

    NAMELIST /vegetation_ctl/  & 
      read_output_interval_pool, &
      read_output_interval_flux, &
      output_start_date, &
      output_end_date, &
      flag_optimal_Nfraction, &
      flag_optimal_Nleaf, &
      flag_calc_soil_uptake
    NAMELIST /forcing_ctl/  &     ! --> 2do: move to atm2land
        create_forcing

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_veg_config'

    CALL message(TRIM(routine), 'Starting veg configuration')

    ! Set defaults
    read_output_interval_pool  = 'monthly'   ! used for read from file
    read_output_interval_flux  = 'monthly'   ! 
    output_start_date          = 1
    output_end_date            = 100
    flag_optimal_Nfraction     = .FALSE.  ! .TRUE. .FALSE.
    flag_optimal_Nleaf         = .FALSE.  ! .TRUE. .FALSE.
    flag_calc_soil_uptake      = .FALSE.  ! .TRUE. .FALSE.
    ! Set defaults
    create_forcing             = 'generate'

    ! &VEGETATION_CTL
    ! 
    ! read the namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('VEGETATION_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, vegetation_ctl)
    CALL close_nml(nml_handler)
    ! pass values as read from file
    config%output_interval_pool   = read_output_interval_pool
    config%output_interval_flux   = read_output_interval_flux
    config%output_start_date      = output_start_date
    config%output_end_date        = output_end_date
    config%flag_optimal_Nfraction = flag_optimal_Nfraction
    config%flag_optimal_Nleaf     = flag_optimal_Nleaf
    config%flag_calc_soil_uptake  = flag_calc_soil_uptake

    ! &FORCING_CTL  --> 2do: move to atm2land config file
    ! 
    ! read the namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('FORCING_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, forcing_ctl)
    CALL close_nml(nml_handler)
    ! pass values as read from file
    config%option_forcing = create_forcing

  END SUBROUTINE Init_veg_config

!#endif
END MODULE mo_veg_config_class
