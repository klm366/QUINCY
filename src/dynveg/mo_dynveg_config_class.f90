!>
!! @brief Contains structures and methods for dynveg config
!!
!! @author
!!  Thomas Raddatz, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Thomas Raddatz (24.03.2020)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_dynveg_config_class
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
  PUBLIC :: t_dynveg_config

  TYPE, EXTENDS(t_jsb_config) :: t_dynveg_config
    LOGICAL             :: &
      & init_running_means
     !X e.g: LOGICAL                     :: use_alb_veg_simple

   CONTAINS
     PROCEDURE :: Init => Init_dynveg_config
  END type t_dynveg_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dynveg_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize dynveg process
  !!
  !! @param[inout]     config     Configuration type of process (t_dynveg_config)
  !!  
  SUBROUTINE Init_dynveg_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC

    CLASS(t_dynveg_config), INTENT(inout) :: config
    TYPE(t_jsb_vgrid), POINTER  :: pft_layer
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: init_running_means

    !X e.g: LOGICAL                     :: use_alb_veg_simple

    NAMELIST /jsb_dynveg_nml/      &
         active,                   &
         ic_filename,              &
         bc_filename,              &
         init_running_means
         !X e.g: use_<PROCESS_NAME_LOWER_CASE>_veg_simple

    INTEGER :: nml_handler, nml_unit, istat, npft, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_dynveg_config'

    CALL message(TRIM(routine), 'Starting dynveg configuration')

    ! Set defaults
    active              = .TRUE.
    bc_filename         = 'bc_land_dynveg.nc'
    ic_filename         = 'ic_land_dynveg.nc'
    init_running_means  = .FALSE.
    !X e.g: use_alb_veg_simple        = .FALSE.

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_dynveg_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_dynveg_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%init_running_means  = init_running_means
    IF (init_running_means) CALL message(TRIM(routine), 'Initializing climate buffer')
    !X e.g: config%use_alb_veg_simple        = use_alb_veg_simple

    ! Create vertical axis
    npft=11 ! TODO here we can not grab it from tile%no_of_children or mem%no_of_children
    pft_layer  => new_vgrid('pfts', ZAXIS_GENERIC, npft, &
      & levels=(/ (REAL(i,kind=wp),i=1,npft) /), units='')
    CALL Register_vgrid(pft_layer)

  END SUBROUTINE Init_dynveg_config

#endif
END MODULE mo_dynveg_config_class
