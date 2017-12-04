!>
!! @brief Contains structures and methods for carbon config
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (20.05.2017)
!!
MODULE mo_carbon_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_carbon_config

  TYPE, EXTENDS(t_jsb_config) :: t_carbon_config
    LOGICAL                     :: read_cpools     ! Read Cpools from external File when model was 
                                                   ! newly startet (Cpools were initialized with fixed values) or
                                                   ! restarted (Cpools were initialized and afterwards overwritten with values from restartfile).
    INTEGER                     :: fuel_algorithm  ! If the fire_algorithm in the process disturbance is set to 1 (jsbach default),
                                                   ! this switch here selects the corresponding fuel algorithm: 0 : none, 1 : jsbach (default), 
                                                   ! 2 : Arora and Boer (2005), 3 : Thonicke et al.(2001), 4 : Read/GFED

    REAL(wp)                    :: fire_frac_wood_2_atmos 

   CONTAINS
     PROCEDURE :: Init => Init_carbon_config
  END type t_carbon_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize carbon process
  !!
  !! @param[inout]     config     Configuration type of process (t_carbon_config)
  !!
  SUBROUTINE Init_carbon_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_carbon_constants,   ONLY: FireFracWood2Atmos

    CLASS(t_carbon_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: read_cpools
    INTEGER                     :: fuel_algorithm  ! Select fuel algorithm: 0 : none, 1 : jsbach (default), 
                                                   ! 2 : Arora and Boer (2005), 3 : Thonicke et al.(2001), 4 : Read/GFED
    REAL(wp)                    :: fire_frac_wood_2_atmos
    ! LOGICAL           :: lemissions 

    NAMELIST /jsb_carbon_nml/      &
         active,                   &
         ic_filename,              &
         bc_filename,              &
         read_cpools,              &
         fuel_algorithm,           &
         fire_frac_wood_2_atmos

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_carbon_config'

    CALL message(TRIM(routine), 'Starting carbon configuration')

    ! Set defaults
    active                 = .FALSE.
    bc_filename            = 'bc_land_carbon.nc'
    ic_filename            = 'ic_land_carbon.nc'
    read_cpools            = .FALSE.
    fuel_algorithm         = 0
    fire_frac_wood_2_atmos = FireFracWood2Atmos

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_carbon_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_carbon_nml)

    CALL close_nml(nml_handler)

    ! Write resulting values in config memory
    config%active                 = active
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename
    config%read_cpools            = read_cpools
    config%fuel_algorithm         = fuel_algorithm
    config%fire_frac_wood_2_atmos = fire_frac_wood_2_atmos

  END SUBROUTINE Init_carbon_config

#endif
END MODULE mo_carbon_config_class
