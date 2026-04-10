!> Climate buffer process configuration
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
!>### Contains the namelist configuration for the climate buffer process
!>
!> 
!>
MODULE mo_climbuf_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_climbuf_config

  TYPE, EXTENDS(t_jsb_config) :: t_climbuf_config
    LOGICAL :: &
      & init_running_means

   CONTAINS
     PROCEDURE :: Init => Init_climbuf_config
  END type t_climbuf_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_climbuf_config_class'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialize climbuf process
  !
  SUBROUTINE Init_climbuf_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_climbuf_config), INTENT(inout) :: config !<Configuration type of process (t_climbuf_config)
    ! -------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    LOGICAL                     :: lrestart_cont
    LOGICAL                     :: init_running_means
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    ! -------------------------------------------------------------------------------------------------- !

    NAMELIST /jsb_climbuf_nml/      &
     &   active,                    &
     &   lrestart_cont,             &
     &   ic_filename,               &
     &   bc_filename,               &
     &   init_running_means

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_climbuf_config'

    CALL message(TRIM(routine), 'Starting climbuf configuration')

    ! Set defaults
    active              = .FALSE.
    lrestart_cont       = .FALSE.      ! TRUE: Continue although climbuf variables are missing in restart file
    bc_filename         = 'bc_land_climbuf.nc'
    ic_filename         = 'ic_land_climbuf.nc'
    init_running_means  = .FALSE.

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_climbuf_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_climbuf_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%lrestart_cont       = lrestart_cont
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%init_running_means  = init_running_means
    IF (init_running_means) CALL message(TRIM(routine), 'Initializing climate buffer')

  END SUBROUTINE Init_climbuf_config

#endif
END MODULE mo_climbuf_config_class
