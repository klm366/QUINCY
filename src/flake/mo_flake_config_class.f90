!>
!! @brief Contains structures and methods for flake config
!!
!! @author
!!  Zoé Rehder, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Zoé Rehder (2023-06-01)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_flake_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module
  
  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USe mo_util_string,       ONLY: real2string
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables 
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_flake_config

  TYPE, EXTENDS(t_jsb_config) :: t_flake_config
    REAL(wp) :: depth  !< Lake depth [m]
    ! Parameters to initialize the flake model with constant values. They are currently the only
    ! option for initialization and can be used with the site-level configuration.
    ! TODO: For the full configuration, these parameters will be read from the input file.
    REAL(wp) :: h_ice  !< Ice thickness [m]
    REAL(wp) :: h_ML   !< Mixed layer thickness [m]
    REAL(wp) :: h_snow !< Snow thickness [m]
    REAL(wp) :: T_bot  !< Bottom temperature [K]
    REAL(wp) :: T_lice !< Ice temperature [K]
    REAL(wp) :: T_mnw  !< Mixed layer temperature [K]
    ! REAL(wp) :: T_sfc  !< Surface temperature [K]
    ! REAL(wp) :: T_snow !< Snow temperature [K]
    REAL(wp) :: T_wML  !< Water temperature in mixed layer [K]
    REAL(wp) :: c_t    !< Shape factor with respect to the temperature profile in the thermocline
    REAL(wp) :: fetch  !< Fetch length [m]
  CONTAINS
      PROCEDURE :: Init => Init_flake_config
  END type t_flake_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_flake_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize flake process
  !!
  !! @param[inout]     config     Configuration type of process (t_flake_config)
  !!  
  SUBROUTINE Init_flake_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_flake_constants, ONLY : &
      & tpl_T_r, & ! Reference temperature [K]
      & tpl_T_f    ! Freezing temperature [K]

    CLASS(t_flake_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    REAL(wp)                    :: depth, &
      & h_ice, h_ML, h_snow, T_bot, T_lice, T_mnw, T_wML, fetch, c_t

    NAMELIST /jsb_flake_nml/      &
         active,                      &
         ic_filename,                 &
         bc_filename, &
         depth,       &
         h_ice,      &
         h_ML,       &
         h_snow,     &
         T_bot,      &
         T_lice,    &
         T_mnw,     &
         ! T_sfc,    &
         ! T_snow,   &
         T_wML,     &
         fetch,     &
         c_t

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_flake_config'

    CALL message(TRIM(routine), 'Starting flake configuration')

    ! Set defaults
    active      = .FALSE.
    bc_filename = 'bc_land_flake.nc'
    ic_filename = 'ic_land_flake.nc'
    depth       = 10.0_wp     !< Lake depth [m]
    h_ice       = 0.0_wp      !< Ice thickness [m]
    h_ML        = 1.0_wp      !< Mixed layer thickness [m]
    h_snow      = 0.0_wp      !< Snow thickness [m]
    T_bot       = tpl_T_r     !< Bottom temperature [K]
    T_lice      = tpl_T_f     !< Ice temperature [K]
    T_mnw       = 280._wp     !< Mean temperature of the water column [K]
    ! T_sfc       = tpl_T_f
    ! T_snow      = tpl_T_f
    T_wML       = 280._wp     !< Water temperature in mixed layer [K]
    fetch       = 1.0E+04_wp  !< Fetch length - 10km [m]
    c_t         = 0.5_wp      !< Shape factor with respect to the temperature profile in the thermocline

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_flake_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_flake_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%depth               = depth
    config%h_ice               = h_ice
    config%h_ML                = h_ML
    config%h_snow              = h_snow
    config%T_bot               = T_bot
    config%T_lice              = T_lice
    config%T_mnw               = T_mnw
    config%T_wML               = T_wML
    config%fetch               = fetch
    config%c_t                 = c_t

    CALL message(routine, 'Using lake depth = '//TRIM(real2string(depth))//' m')
    ! Snow on lake ice is not supported in the current implementation of flake
    IF (h_snow > 0.0_wp) CALL finish(routine, 'Snow thickness is greater than 0, which is not supported')
    
  END SUBROUTINE Init_flake_config

#endif
END MODULE mo_flake_config_class
