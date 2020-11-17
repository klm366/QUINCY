!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Contains structures and methods for radiation config
!! 
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!! 
!! @par Revision History
!! First version by Rainer Schneck (2015-04-27)
!! implementation to SLM Steffen Richter (2017-12-14)
!! 
MODULE mo_rad_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message, finish, message_text
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_rad_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the radiation process, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_rad_config
    ! ..
  CONTAINS
    PROCEDURE :: Init => Init_rad_config
  END type t_rad_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_rad_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_rad_config(config)

    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_grid_class,     ONLY: ZAXIS_GENERIC   ! in JSBACH4 used from mo_jsb_io    *slmtemp

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_rad_config), INTENT(inout) :: config
    ! ---------------------------
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    INTEGER                     :: nspec
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_spectral
    

    NAMELIST /radiation_ctl/      &
         active,                      &
         ic_filename,                 &
         bc_filename


    ! variables for reading model-options from namelist
    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_rad_config'

    CALL message(TRIM(routine), 'Starting radiation configuration')

    ! Set defaults
    active              = .TRUE.
    bc_filename         = 'bc_land_rad.nc'  !R: gibts im Moment noch nicht
    ic_filename         = 'ic_land_rad.nc'  !R: gibts im Moment noch nicht
    nspec               = 2162

    ! read the namelist
    ! &RADIATION_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('RADIATION_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, radiation_ctl )

    CALL close_nml(nml_handler)

    ! pass values as read from file
    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename


    !< Create the 3rd dimension of spectral bands of the light for RAD_ variables using a vertical axis
    !! 
    vgrid_spectral  => new_vgrid('spectral_bands', ZAXIS_GENERIC, nspec, &
                      longname='Number of spectral bands the light is divided into, from RAD_ QUINCY') !, &
                      ! units=''  , &
                      ! levels=   , &
                      ! lbounds=  , & 
                      ! ubounds=  , &
                      ! dz=       )
    CALL register_vgrid(vgrid_spectral)

  END SUBROUTINE Init_rad_config

#endif
END MODULE mo_rad_config_class
