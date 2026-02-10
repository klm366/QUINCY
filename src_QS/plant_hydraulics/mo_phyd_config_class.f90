!> QUINCY plant hydraulics process config
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### define phenology config structure, read phenology namelist and init configuration parameters
!>
MODULE mo_phyd_config_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_phyd_config

  ! ======================================================================================================= !
  !> PHYD_ configuration parameters
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_phyd_config
    LOGICAL            :: use_plant_hydraulics           !< on/off whether using plant hydraulics
  CONTAINS
      PROCEDURE :: Init => Init_phyd_config
  END type t_phyd_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_phyd_config_class'

CONTAINS

  ! ======================================================================================================= !
  !> read P_HYDRO_ namelist and init configuration parameters
  !>
  SUBROUTINE Init_phyd_config(config)
    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_phyd_config), INTENT(inout) :: config
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename
    LOGICAL                     :: use_plant_hydraulics

    NAMELIST /lnd_phyd_nml/           &
      & active,                     &
      & ic_filename,                &
      & bc_filename,                &
      & use_plant_hydraulics
    INTEGER :: nml_handler, nml_unit, istat
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_phyd_config'
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting plant hydraulic architecture configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    active                    = .FALSE.
    bc_filename               = 'bc_land_plant_hydraulics.nc'
    ic_filename               = 'ic_land_plant_hydraulics.nc'
    use_plant_hydraulics      = .FALSE.
    ! ----------------------------------------------------------------------------------------------------- !
    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
    ! ----------------------------------------------------------------------------------------------------- !
    nml_unit = position_nml('lnd_phyd_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_phyd_nml)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Write namelist values into config
    config%active                     = active
    config%ic_filename                = ic_filename
    config%bc_filename                = bc_filename
    config%use_plant_hydraulics       = use_plant_hydraulics
    CALL close_nml(nml_handler)
  END SUBROUTINE Init_phyd_config

#endif
END MODULE mo_phyd_config_class
