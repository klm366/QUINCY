!>
!>  @par Copyright
!>    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!>    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!>    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!>
!>    This file is part of the QUINCY model.
!>    The QUINCY model is free software: you can redistribute it and/or modify
!>    it under the terms of the GNU General Public License version 3,
!>    as published by the Free Software Foundation.
!>
!>    The QUINCY model is distributed in the hope that it will be useful,
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!>    GNU General Public License for more details.
!>
!>    You should have received a copy of the GNU General Public License
!>    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!> 
!> @brief Contains structures and methods for config of the dist_fire_fuel process
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_dfire_fuel_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message, finish, message_text
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: tolower
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_dfire_fuel_config

  ! ======================================================================================================= !
  !>
  !> config TYPE of the dist_fire process, derived from t_jsb_config
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_dfire_fuel_config
    LOGICAL           :: flag_dfire_fuel              !< on/off disturbance fire process [boolean]

  CONTAINS
    PROCEDURE :: Init => Init_dfire_fuel_config
  END type t_dfire_fuel_config
  
  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_fuel_config_class'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> configuration routine of t_dfire_fuel_config
  !>
  SUBROUTINE Init_dfire_fuel_config(config)
#ifdef __QUINCY_STANDALONE__
    USE mo_namelist,            ONLY: open_nml, POSITIONED, position_nml, close_nml
#else
    USE mo_jsb_namelist_iface,  ONLY: open_nml, POSITIONED, position_nml, close_nml
#endif
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_io_iface,       ONLY: zaxis_generic
    USE mo_jsb_grid,           ONLY: register_vgrid
    USE mo_dfire_constants,    ONLY: df_sf_vgrid_dead_fuel_pools_name

    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_dfire_fuel_config), INTENT(inout) :: config    !< config type for dfire_fuel
    ! ----------------------------------------------------------------------------------------------------- !
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                       :: active
    CHARACTER(len=filename_max)   :: ic_filename, &
                                     bc_filename
    LOGICAL                       :: flag_dfire_fuel            !< on/off disturbance fire process [boolean]

    CHARACTER(len=*), PARAMETER   :: routine = TRIM(modname)//':Init_dfire_fuel_config'
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_STANDALONE__
    NAMELIST /dist_fire_fuel_ctl/        &
#else
    NAMELIST /dist_fire_fuel_nml/        &
#endif
      & active,                     &
      & ic_filename,                &
      & bc_filename,                &
      & flag_dfire_fuel

    ! ----------------------------------------------------------------------------------------------------- !
    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting dist_fire_fuel configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults, see namelist.slm for more details
    active                   = .FALSE.
    bc_filename              = 'bc_land_dfire_fuel.nc'
    ic_filename              = 'ic_land_dfire_fuel.nc'
    flag_dfire_fuel          = .FALSE.

    ! ----------------------------------------------------------------------------------------------------- !
    ! read the namelist
    ! &DIST_FIRE_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_STANDALONE__
    nml_unit = position_nml('dist_fire_fuel_ctl', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, dist_fire_fuel_ctl)
#else
    nml_unit = position_nml('dist_fire_fuel_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, dist_fire_fuel_nml)
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    ! pass values as read from file
    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%flag_dfire_fuel           = flag_dfire_fuel
    ! ----------------------------------------------------------------------------------------------------- !
    CALL close_nml(nml_handler)
    ! ----------------------------------------------------------------------------------------------------- !
    vgrid_dead_fuel_pools  => new_vgrid(df_sf_vgrid_dead_fuel_pools_name, zaxis_generic, t_sf_fuel_dead_fuel_nr, &
        longname='Total number of dead fuel pools', &
        & units='g DM m-2')
    CALL register_vgrid(vgrid_fuel_class)
  END SUBROUTINE Init_dfire_fuel_config

#endif
END MODULE mo_dfire_fuel_config_class
