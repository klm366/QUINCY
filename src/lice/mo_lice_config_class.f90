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
!> The lice process handles the coupling of the ICON-Land land ice tile 'lice' (a
!> second glacier tile for the cells covered by the external land ice model, currently
!> the EBFM firn model) to that model via YAC.
!>
MODULE mo_lice_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lice_config

  TYPE, EXTENDS(t_jsb_config) :: t_lice_config
    ! The set of exchanged fields is defined by the coupling configuration (coupling.yaml),
    ! not by namelist switches.
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

    NAMELIST /jsb_lice_nml/    &
      active,                  &
      lrestart_cont,           &
      ic_filename,             &
      bc_filename

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_lice_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting lice configuration')

    ! Set defaults
    active         = .FALSE.
    lrestart_cont  = .TRUE.     ! Allow restarts from files written without the lice process
    ic_filename    = 'ic_land_lice.nc'
    bc_filename    = 'bc_land_lice.nc'

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_lice_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_lice_nml)

    CALL close_nml(nml_handler)

    config%active         = active
    config%lrestart_cont  = lrestart_cont
    config%ic_filename    = ic_filename
    config%bc_filename    = bc_filename

    IF (.NOT. active) RETURN

    ! The lice tile and the EBFM-specific branches of the processes exist only in the TMX usecases
    IF (.NOT. config%model_config%use_tmx) THEN
      CALL finish(TRIM(routine), 'The lice process requires use_tmx=.TRUE. (jsb_model_nml)')
    END IF
    IF (.NOT. config%model_config%use_glacier) THEN
      CALL finish(TRIM(routine), 'The lice process requires use_glacier=.TRUE. (jsb_model_nml)')
    END IF

  END SUBROUTINE Init_lice_config

#endif
END MODULE mo_lice_config_class
