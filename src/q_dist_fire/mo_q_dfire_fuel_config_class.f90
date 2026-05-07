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
MODULE mo_q_dfire_fuel_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message, finish, message_text
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: tolower
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_q_dfire_fuel_config

  ! ======================================================================================================= !
  !>
  !> config TYPE of the dist_fire process, derived from t_jsb_config
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_q_dfire_fuel_config
    LOGICAL           :: flag_dfire_fuel              !< on/off disturbance fire process [boolean]

  CONTAINS
    PROCEDURE :: Init => Init_q_dfire_fuel_config
  END type t_q_dfire_fuel_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_fuel_config_class'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> configuration routine of t_q_dfire_fuel_config
  !>
  SUBROUTINE Init_q_dfire_fuel_config(config)
    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_io_iface,       ONLY: zaxis_generic
    USE mo_jsb_grid,           ONLY: register_vgrid
    USE mo_q_dfire_constants,  ONLY: df_sf_vgrid_dead_fuel_pools_name
    USE mo_q_dfire_constants,  ONLY: df_sf_vgrid_live_fuel_pools_name, &
                                      df_sf_vgrid_elements_name, &
                                      df_sf_vgrid_species_name
    USE mo_q_dfire_constants,  ONLY: sf_fuel_dead_fuel_nr, sf_fuel_live_fuel_nr
    USE mo_q_dfire_constants,  ONLY: df_sf_elements_nr, df_sf_species_nr


    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_dfire_fuel_config), INTENT(inout) :: config    !< config type for dfire_fuel
    ! ----------------------------------------------------------------------------------------------------- !
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                       :: active
    CHARACTER(len=filename_max)   :: ic_filename, &
                                     bc_filename
    LOGICAL                       :: flag_dfire_fuel            !< on/off disturbance fire process [boolean]

    TYPE(t_jsb_vgrid), POINTER    :: vgrid_dead_fuel_pools
    TYPE(t_jsb_vgrid), POINTER    :: vgrid_live_fuel_pools
    TYPE(t_jsb_vgrid), POINTER    :: vgrid_species
    TYPE(t_jsb_vgrid), POINTER    :: vgrid_elements

    CHARACTER(len=*), PARAMETER   :: routine = TRIM(modname)//':Init_q_dfire_fuel_config'
    ! ----------------------------------------------------------------------------------------------------- !
    NAMELIST /lnd_q_dist_fire_fuel_nml/        &
      & active,                     &
      & ic_filename,                &
      & bc_filename,                &
      & flag_dfire_fuel

    ! ----------------------------------------------------------------------------------------------------- !
    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting lnd_q_dist_fire_fuel configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults, see namelist.slm for more details
    active                   = .FALSE.
    bc_filename              = 'bc_land_dfire_fuel.nc'
    ic_filename              = 'ic_land_dfire_fuel.nc'
    flag_dfire_fuel          = .FALSE.

    ! ----------------------------------------------------------------------------------------------------- !
    ! read the namelist
    ! &Q_DIST_FIRE_FUEL_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))
    ! ----------------------------------------------------------------------------------------------------- !
    nml_unit = position_nml('lnd_q_dist_fire_fuel_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_q_dist_fire_fuel_nml)
    ! ----------------------------------------------------------------------------------------------------- !
    ! pass values as read from file
    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%flag_dfire_fuel           = flag_dfire_fuel
    ! ----------------------------------------------------------------------------------------------------- !
    CALL close_nml(nml_handler)
    ! ----------------------------------------------------------------------------------------------------- !
    vgrid_dead_fuel_pools  => new_vgrid(df_sf_vgrid_dead_fuel_pools_name, zaxis_generic, sf_fuel_dead_fuel_nr, &
        longname='Total number of dead fuel pools', &
        & units='g DM m-2')
    vgrid_live_fuel_pools  => new_vgrid(df_sf_vgrid_live_fuel_pools_name, zaxis_generic, sf_fuel_live_fuel_nr, &
        longname='Total number of live fuel pools', &
        & units='g DM m-2')
    vgrid_species  => new_vgrid(df_sf_vgrid_species_name, zaxis_generic, df_sf_species_nr, &
        longname='Emission Species', &
        & units='g X Kg-1 DM-1 m-2')
    vgrid_elements  => new_vgrid(df_sf_vgrid_elements_name, zaxis_generic, df_sf_elements_nr, &
        longname='Elements', &
        & units='mol X m-2')
    CALL register_vgrid(vgrid_dead_fuel_pools)
    CALL register_vgrid(vgrid_live_fuel_pools)
    CALL register_vgrid(vgrid_species)
    CALL register_vgrid(vgrid_elements)
  END SUBROUTINE Init_q_dfire_fuel_config

#endif
END MODULE mo_q_dfire_fuel_config_class
