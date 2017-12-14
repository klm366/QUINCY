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
  PUBLIC :: t_veg_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the vegetation process, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_veg_config
    INTEGER          :: pft_id
    CHARACTER(15)    :: bnf_scheme                 !< select scheme to simulate BNF (biological nitrogen fixation)
                                                   ! "none|fixed|standard|resistance|optimal"   
    LOGICAL          :: flag_tile_veg_dynamics     !< calculate within-tile vegetation dynamics + start from bareground
    LOGICAL          :: flag_optimal_Nleaf         !< optimise leaf C:N 
    LOGICAL          :: flag_optimal_allocation    !< optimise plant allocation
   CONTAINS
     PROCEDURE :: Init => Init_veg_config
  END type t_veg_config

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
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    INTEGER       :: plant_functional_type_id
    CHARACTER(15) :: bnf_scheme       
    LOGICAL       :: flag_tile_veg_dynamics
    LOGICAL       :: flag_optimal_Nleaf
    LOGICAL       :: flag_optimal_allocation
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_veg_config'

    NAMELIST /vegetation_ctl/       & 
         active,                 &
         ic_filename,            &
         bc_filename,            &
      bnf_scheme,                   & 
      plant_functional_type_id,     &
      flag_tile_veg_dynamics,       &
      flag_optimal_Nleaf,           &
      flag_optimal_allocation

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting veg configuration')

    ! Set defaults
    active              = .TRUE.
    bc_filename         = 'bc_land_phys.nc'
    ic_filename         = 'ic_land_phys.nc'
    bnf_scheme                     = "standard"  ! standard std_plus_costs resistance optimal
    plant_functional_type_id       = 1
    flag_tile_veg_dynamics         = .FALSE.  ! .TRUE. .FALSE.
    flag_optimal_Nleaf             = .FALSE.  ! .TRUE. .FALSE.
    flag_optimal_allocation        = .FALSE.  ! .TRUE. .FALSE.

    ! read the namelist
    ! &VEGETATION_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('VEGETATION_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, vegetation_ctl)

    CALL close_nml(nml_handler)

    ! check input for BNF scheme
    SELECT CASE (TRIM(bnf_scheme))
     CASE ("none","fixed","standard","resistance","optimal")
     CASE DEFAULT
       WRITE(message_text,'(2a)') 'invalid bnf_scheme',bnf_scheme
       CALL finish("mo_veg_config_class:init_veg_config",message_text) 
    ENDSELECT

    ! pass values as read from file
    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%bnf_scheme                    = TRIM(bnf_scheme)
    config%pft_id                        = plant_functional_type_id
    config%flag_tile_veg_dynamics        = flag_tile_veg_dynamics  
    config%flag_optimal_Nleaf            = flag_optimal_Nleaf
    config%flag_optimal_allocation       = flag_optimal_allocation

  END SUBROUTINE Init_veg_config

!#endif
END MODULE mo_veg_config_class
