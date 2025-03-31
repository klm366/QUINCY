!> QUINCY agriculture process config
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### define simple harvest config structure, read simple harvest namelist and init configuration parameters
!>
MODULE mo_q_agr_config_class
#ifndef __NO_QUINCY__

  USE mo_exception,         ONLY: message, message_text, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_math_constants,ONLY: eps4
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_q_agr_config

  ! ======================================================================================================= !
  !> Q_AGR_ configuration parameters
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_q_agr_config
    LOGICAL   :: l_read_xyz               !< determines if x is to be read from a file (default FALSE) (ignored in QS)
    REAL(wp)  :: xyz                      !< default xyz used if not read from a map

  CONTAINS
    PROCEDURE :: Init => Init_q_agr_config
  END TYPE t_q_agr_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_agr_config_class'

CONTAINS

  ! ======================================================================================================= !
  !> read Q_AGR_ namelist and init configuration parameters
  !>
  SUBROUTINE Init_q_agr_config(config)
    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,              ONLY: Register_vgrid
    USE mo_jsb_io,                ONLY: ZAXIS_GENERIC
    USE mo_jsb_math_constants,    ONLY: eps4
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_agr_config), INTENT(inout) :: config
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename

    LOGICAL   :: l_read_xyz
    REAL(wp)  :: xyz

#ifdef __QUINCY_STANDALONE__
    NAMELIST /q_agr_ctl/  &
#else
    NAMELIST /q_agr_nml/  &
#endif
      & active,                   &
      & ic_filename,              &
      & bc_filename,              &
      & l_read_xyz,               &
      & xyz

    INTEGER                     :: nml_handler, nml_unit, istat
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_q_agr_config'
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting agriculture configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    active                        = .FALSE.                 ! what should be the default?
    bc_filename                   = 'bc_land_agr.nc'
    ic_filename                   = 'ic_land_agr.nc'
    !
    l_read_xyz            = .FALSE.
    xyz                   = 0.05_wp

    ! ----------------------------------------------------------------------------------------------------- !
    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
#ifdef __QUINCY_STANDALONE__
    nml_unit = position_nml('q_agr_ctl', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, q_agr_ctl)
#else
    nml_unit = position_nml('q_agr_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, q_agr_nml)
#endif

    ! ----------------------------------------------------------------------------------------------------- !
    ! Write namelist values into config
    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%l_read_xyz                = l_read_xyz
    config%xyz                       = xyz

    CALL close_nml(nml_handler)

  END SUBROUTINE Init_q_agr_config

#endif
END MODULE mo_q_agr_config_class
