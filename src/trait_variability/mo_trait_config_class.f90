! Author: Dr. Sri Nandini-Weiß
! Date: 2026-07-15
!> Trait-variability process configuration
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
!>#### Contains namelist information for the trait-variability process
!>
MODULE mo_trait_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,        ONLY: message
  USE mo_io_units,         ONLY: filename_max
  USE mo_jsb_config_class, ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_trait_config

  TYPE, EXTENDS(t_jsb_config) :: t_trait_config
    ! Step 6 active-trait switch.
    ! Author: Dr. Sri Nandini-Weiß
    ! Date: 2026-07-27
    ! .FALSE.: retain the validated diagnostic-only behavior.
    ! .TRUE.: use bounded annual traits in the JSBACH process calculations, for those
    ! plant functional types that have a fitted model.
    ! Active annual trait use is opt-in: FALSE preserves the
    ! diagnostic-only configuration and prevents unintended physics changes.
    ! Active annual PFT13 trait use is opt-in: FALSE preserves the
    ! diagnostic-only configuration and prevents unintended physics changes.
    LOGICAL :: l_apply_dynamic_traits = .FALSE.
  CONTAINS
    PROCEDURE :: Init => Init_trait_config
  END TYPE t_trait_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_trait_config_class'

CONTAINS

  SUBROUTINE Init_trait_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_trait_config), INTENT(inout) :: config

    LOGICAL :: active, lrestart_cont, l_apply_dynamic_traits
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    NAMELIST /jsb_trait_nml/ &
      & active,                  &
      & lrestart_cont,           &
      & l_apply_dynamic_traits,    &
      & ic_filename,         &
      & bc_filename

    INTEGER :: nml_handler, nml_unit, istat
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_trait_config'

    CALL message(TRIM(routine), 'Starting trait configuration')

    active               = .FALSE.
    lrestart_cont        = .FALSE.
    l_apply_dynamic_traits = .FALSE.
    bc_filename   = 'bc_land_trait.nc'
    ic_filename   = 'ic_land_trait.nc'

    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('jsb_trait_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_trait_nml)
    CALL close_nml(nml_handler)

    config%active               = active
    config%lrestart_cont        = lrestart_cont
    config%l_apply_dynamic_traits = l_apply_dynamic_traits
    config%ic_filename   = ic_filename
    config%bc_filename   = bc_filename

  END SUBROUTINE Init_trait_config

#endif
END MODULE mo_trait_config_class
