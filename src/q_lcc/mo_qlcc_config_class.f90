!> Contains structures and methods for qlcc config
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
MODULE mo_qlcc_config_class
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
  PUBLIC :: t_qlcc_config

  TYPE, EXTENDS(t_jsb_config) :: t_qlcc_config
    LOGICAL             :: &
      & init_running_means, &
      & init_qlcc
    REAL(wp)            :: &
      & accelerate_qlcc

   CONTAINS
    PROCEDURE :: Init => Init_qlcc_config
  END type t_qlcc_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_qlcc_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize qlcc process
  !!
  !! @param[inout]     config     Configuration type of process (t_qlcc_config)
  !!
  SUBROUTINE Init_qlcc_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC

    CLASS(t_qlcc_config), INTENT(inout) :: config
    TYPE(t_jsb_vgrid), POINTER  :: pft_layer
    LOGICAL                     :: active, lrestart_cont
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: init_running_means
    LOGICAL                     :: init_qlcc
    REAL(wp)                    :: accelerate_qlcc

    NAMELIST /jsb_qlcc_nml/  &
      &  active,             &
      &  lrestart_cont,      &
      &  ic_filename,        &
      &  bc_filename,        &
      &  init_running_means, &
      &  init_qlcc,          &
      &  accelerate_qlcc

    INTEGER :: nml_handler, nml_unit, istat, npft, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_qlcc_config'

    CALL message(TRIM(routine), 'Starting qlcc configuration')

    ! Set defaults
    active              = .FALSE.
    lrestart_cont       = .FALSE.      ! TRUE: Continue although QLCC variables are missing in restart file
    bc_filename         = 'bc_land_qlcc.nc'
    ic_filename         = 'ic_land_qlcc.nc'
    init_running_means  = .FALSE.
    init_qlcc           = .FALSE.
    accelerate_qlcc     = 1._wp

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_qlcc_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_qlcc_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%lrestart_cont       = lrestart_cont
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%init_running_means  = init_running_means
    IF (init_running_means) CALL message(TRIM(routine), 'Initializing climate buffer')
    config%init_qlcc         = init_qlcc
    IF (init_qlcc) CALL message(TRIM(routine), 'Initializing QLCC')
    config%accelerate_qlcc   = accelerate_qlcc
    IF (accelerate_qlcc /= 1._wp) THEN
      WRITE(message_text, *) 'QLCC calculations accelerated by factor ', accelerate_qlcc
      CALL message(TRIM(routine), message_text)
    END IF

    ! Create vertical axis
    npft=11 ! TODO here we can not grab it from tile%no_of_children or mem%no_of_children
    pft_layer  => new_vgrid('pfts', ZAXIS_GENERIC, npft, &
      & levels=(/ (REAL(i,kind=wp),i=1,npft) /), units='')
    CALL Register_vgrid(pft_layer)

  END SUBROUTINE Init_qlcc_config

#endif
END MODULE mo_qlcc_config_class
