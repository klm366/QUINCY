!>
!! @brief Contains structures and methods for dynveg config
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2018-11-28)
!!
MODULE mo_dynveg_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
!  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_dynveg_config, Init_dynveg_config

  TYPE, EXTENDS(t_jsb_config) :: t_dynveg_config
    CHARACTER(len=15) :: scheme   !< Phenology scheme to use ('logrop', 'climatology')
     !LOGICAL :: first_call_dynveg_this_day  ! TRUE marks if a new day has just started and update_dynveg is called for the top tile
    CONTAINS
      PROCEDURE :: Init => Init_dynveg_config
  END type t_dynveg_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dynveg_config_class'

CONTAINS

  SUBROUTINE Init_dynveg_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid !, new_vgrid
!    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC
       
    CLASS(t_dynveg_config), INTENT(inout) :: config

    LOGICAL :: active
    CHARACTER(len=15) :: scheme
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    TYPE(t_jsb_vgrid), POINTER  :: pft_layer

    NAMELIST /jsb_dynveg_nml/ &
      & active,               &
      & scheme,               &
      & ic_filename,          &
      & bc_filename             ! this file is NOT used as a (maybe very big) maps file if DYNVEG_%scheme is set to 'maps'
                                ! because we already have yearly files in the pool and can skip merging them.
                                ! However, we have to rethink this because this is structurally wrong.

    INTEGER :: nml_handler, nml_unit, istat !, npft, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_dynveg_config'

    CALL message(TRIM(routine), 'Starting dynveg configuration')

    ! Set defaults
    active              = .FALSE.
    scheme              = 'maps'
    bc_filename         = 'bc_land_dynveg.nc'
    ic_filename         = 'ic_land_dynveg.nc'
    !first_call_dynveg_this_day = .FALSE.

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_dynveg_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_dynveg_nml)

    ! Write resulting values in config memory
    config%active              = active
    config%scheme              = TRIM(scheme)
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename

    CALL close_nml(nml_handler)

!     ! Create vertical axis
!     npft=11 ! R: here we can not grab it from tile%no_of_children or mem%no_of_children
!     pft_layer  => new_vgrid('pfts', ZAXIS_GENERIC, npft, &
!       & levels=(/ (REAL(i,kind=wp),i=1,npft) /), units='')
!     CALL register_vgrid(pft_layer)


  END SUBROUTINE Init_dynveg_config

#endif
END MODULE mo_dynveg_config_class
