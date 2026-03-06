!> Contains structures and methods for land_ice config
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Contains namelist information for the land-ice coupling process (e.g. Elmer/Ice)
!>
MODULE mo_land_ice_config_class
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_LAND_ICE__)
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
  PUBLIC :: t_land_ice_config

  TYPE, EXTENDS(t_jsb_config) :: t_land_ice_config
    INTEGER :: coupling_interval     !< Coupling interval with ice-sheet model [timesteps]
  CONTAINS
    PROCEDURE :: Init => Init_land_ice_config
  END TYPE t_land_ice_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_land_ice_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize land_ice process configuration
  !!
  SUBROUTINE Init_land_ice_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_land_ice_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    INTEGER                     :: coupling_interval

    NAMELIST /jsb_land_ice_nml/       &
      active,                         &
      ic_filename,                    &
      bc_filename,                    &
      coupling_interval

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_land_ice_config'

    CALL message(TRIM(routine), 'Starting land_ice configuration')

    ! Set defaults
    active              = .FALSE.
    bc_filename         = 'bc_land_ice.nc'
    ic_filename         = 'ic_land_ice.nc'
    coupling_interval   = 1

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_land_ice_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_land_ice_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%coupling_interval   = coupling_interval

  END SUBROUTINE Init_land_ice_config

#endif
END MODULE mo_land_ice_config_class
