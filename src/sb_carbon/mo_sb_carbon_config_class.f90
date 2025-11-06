!> Contains structures and methods for carbon config
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
MODULE mo_sb_carbon_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_sb_carbon_config

  TYPE, EXTENDS(t_jsb_config) :: t_sb_carbon_config
    LOGICAL  :: read_sb_cpools ! Read sb Cpools from external file when experiment is newly started (sb cpools are
                            ! initialized with values from file) or restarted (sb cpools from restart file are
                            ! overwritten by values from file).
                            ! ATTENTION: The run script has to take care that read_sbcpools is set to .FALSE.
                            !            in the experiment jobs following the first one! Or re-set it manually!
    LOGICAL  :: diag_maoc_fluxes  ! Activate diagnostic carbon fluxes to and from the maoc pool (jsm only)

   CONTAINS
     PROCEDURE :: Init => Init_sb_carbon_config
  END type t_sb_carbon_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_carbon_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize carbon process
  !!
  !! @param[inout]     config     Configuration type of process (t_sb_carbon_config)
  !!
  SUBROUTINE Init_sb_carbon_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_sb_carbon_config), INTENT(inout) :: config

    LOGICAL                     :: active, lrestart_cont
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: read_cpools, diag_humus_fluxes
    ! LOGICAL           :: lemissions

    NAMELIST /jsb_sb_carbon_nml/      &
         active,                   &
         lrestart_cont,            &
         ic_filename,              &
         bc_filename,              &
         read_sbcpools,            &
         diag_maoc_fluxes

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_sb_carbon_config'

    CALL message(TRIM(routine), 'Starting sb carbon configuration')

    ! Set defaults
    active                 = .FALSE.
    lrestart_cont          = .FALSE.       ! TRUE: Continue although CARBON variables are missing in restart file
    bc_filename            = 'bc_land_carbon.nc'
    ic_filename            = 'ic_land_carbon.nc'
    read_sbcpools          = .FALSE.
    diag_maoc_fluxes       = .FALSE.

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_sb_carbon_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_sb_carbon_nml)

    CALL close_nml(nml_handler)

    ! Write resulting values in config memory
    config%active                 = active
    config%lrestart_cont          = lrestart_cont
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename
    config%read_sbcpools          = read_sbcpools
    config%diag_maoc_fluxes       = diag_maoc_fluxes

  END SUBROUTINE Init_sb_carbon_config

#endif
END MODULE mo_sb_carbon_config_class
