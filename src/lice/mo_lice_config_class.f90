!> Contains structures and methods for the configuration of the land ice (lice) process
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
!>#### Configuration of the land ice (lice) process
!>
!> The lice process handles the coupling of the ICON-Land glacier tile to external
!> land ice models (currently the EBFM firn model) via YAC. In the first version,
!> only the ice fraction is received from EBFM.
!>
MODULE mo_lice_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lice_config

  TYPE, EXTENDS(t_jsb_config) :: t_lice_config
    LOGICAL :: l_get_icefract   !< receive ice fraction from the external land ice model (EBFM)
    LOGICAL :: l_get_albedo     !< receive surface albedo from the external land ice model (EBFM)
    LOGICAL :: l_put_fields     !< send fields to the external land ice model (not yet implemented)
  CONTAINS
    PROCEDURE :: Init => Init_lice_config
  END TYPE t_lice_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lice_config_class'

CONTAINS

  ! ======================================================================================================= !
  !> Initialize lice process configuration (read namelist jsb_lice_nml)
  !>
  SUBROUTINE Init_lice_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_lice_config), INTENT(inout) :: config

    LOGICAL                     :: active          !< lice process is active
    LOGICAL                     :: lrestart_cont   !< continue although lice variables are missing in restart file
    CHARACTER(len=filename_max) :: ic_filename     !< file with initial conditions (not used)
    CHARACTER(len=filename_max) :: bc_filename     !< file with boundary conditions (not used)
    LOGICAL                     :: l_get_icefract  !< receive ice fraction from external land ice model
    LOGICAL                     :: l_get_albedo    !< receive surface albedo from external land ice model
    LOGICAL                     :: l_put_fields    !< send fields to external land ice model

    NAMELIST /jsb_lice_nml/    &
      active,                  &
      lrestart_cont,           &
      ic_filename,             &
      bc_filename,             &
      l_get_icefract,          &
      l_get_albedo,            &
      l_put_fields

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_lice_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting lice configuration')

    ! Set defaults
    active         = .FALSE.
    lrestart_cont  = .TRUE.     ! Allow restarts from files written without the lice process
    ic_filename    = 'ic_land_lice.nc'
    bc_filename    = 'bc_land_lice.nc'
    l_get_icefract = .TRUE.
    l_get_albedo   = .TRUE.
    l_put_fields   = .FALSE.

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_lice_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_lice_nml)

    CALL close_nml(nml_handler)

    config%active         = active
    config%lrestart_cont  = lrestart_cont
    config%ic_filename    = ic_filename
    config%bc_filename    = bc_filename
    config%l_get_icefract = l_get_icefract
    config%l_get_albedo   = l_get_albedo
    config%l_put_fields   = l_put_fields

    IF (.NOT. active) RETURN

    IF (l_put_fields) CALL message(TRIM(routine), &
      & 'WARNING: jsb_lice_nml: l_put_fields=.TRUE. is not implemented yet and has no effect.')

  END SUBROUTINE Init_lice_config

#endif
END MODULE mo_lice_config_class
