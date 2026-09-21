!> Contains the YAC component of ICON-Land and its registration with the coupler
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
!>#### YAC component of ICON-Land
!>
!> ICON-Land takes part in a coupled run as its own YAC component (`land_comp_name`),
!> separate from the component of the host model, so that its coupling partners
!> (currently the EBFM firn model, see mo_interface_lice_ebfm) address the land model and
!> not the atmosphere in the coupling configuration.
!>
!> Inside ICON, ICON-Land runs on the processes of the atmosphere: the component is
!> then *predefined* with the coupler by jsb_construct_coupling, which the host calls
!> before it defines its own components (cpl_def_main); YAC completes the definition of
!> all components of the processes in that later call. The fields of the land component
!> are defined afterwards on the grid and points of the host (e.g. jsb_fdef_lice_fields),
!> so nothing about the land model has to be known by the coupling utilities of the host.
!>
!> In a standalone ICON-Land run the land component would be the main component of its
!> own processes and defined by cpl_def_main instead (not implemented yet).
!>
MODULE mo_jsb_coupling
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: finish, message
  USE mo_coupling_config,     ONLY: is_coupled_run, is_coupled_to_ebfm_land
  USE mo_coupling_utils,      ONLY: cpl_predef_comp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: land_comp_name, land_comp_id
  PUBLIC :: jsb_construct_coupling

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_coupling'

  !> Name of the YAC component of ICON-Land.
  !> Must match the component name used in the coupling configuration (coupling.yaml;
  !> see [[coupling.yaml]] [[[land]]] in ICON's run/mkexp/defaults/DEFAULT.config).
  CHARACTER(len=*), PARAMETER :: land_comp_name = 'icon-land'

  !> YAC component id of the land component (-1 until jsb_construct_coupling has run)
  INTEGER, PROTECTED :: land_comp_id = -1

CONTAINS

  ! ======================================================================================================= !
  !> Register the YAC component of ICON-Land on the processes of the host model
  !>
  !> To be called by the host before it defines its own components with cpl_def_main.
  !> Does nothing unless a coupling of ICON-Land is active.
  !>
  SUBROUTINE jsb_construct_coupling(caller)

    CHARACTER(len=*), INTENT(in) :: caller  !< name of the calling routine (for messages)

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsb_construct_coupling'

    IF (.NOT. is_coupled_run()) RETURN
    ! Couplings of ICON-Land that need the land component (so far only ICON-Land <-> EBFM)
    IF (.NOT. is_coupled_to_ebfm_land()) RETURN

    IF (land_comp_id >= 0) &
      & CALL finish(routine, 'The YAC component ' // land_comp_name // ' has already been defined')

    CALL cpl_predef_comp(TRIM(caller) // ':' // routine, land_comp_name, land_comp_id)

    CALL message(routine, 'Predefined the YAC component ' // land_comp_name // ' of ICON-Land')

  END SUBROUTINE jsb_construct_coupling

#endif
END MODULE mo_jsb_coupling
