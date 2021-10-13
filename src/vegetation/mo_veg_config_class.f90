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
!! @brief Contains structures and methods for veg config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04);
!! replaces the previous SLM mo_veg_config_class created 2016-03-17 that was based on the JSBACH4 mo_veg_config_class anyway;
!! integration into SLM by Jan Engel, MPI BGC (20160401)
!!
MODULE mo_veg_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_veg_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the vegetation process, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_veg_config
    ! INTEGER          :: pft_id
    ! REAL(wp)         :: cohort_harvest_interval    !< interval of harvesting plants [yr], used with veg_dynamics_scheme="cohort"
    ! CHARACTER(15)    :: bnf_scheme                 !< select scheme to simulate BNF (biological nitrogen fixation)
    !                                                !! "none fixed standard optimal resistance unlimited"   
    ! CHARACTER(15)    :: veg_dynamics_scheme        !< select scheme to calculate within-tile vegetation dynamics + start from bareground
    !                                                !! none: constant mortality
    !                                                !! "none population cohort"   
    ! CHARACTER(15)    :: biomass_alloc_scheme       !< select scheme for veg biomass allocation: fixed dynamic optimal
    ! CHARACTER(15)    :: leaf_stoichom_scheme       !< select scheme for C/N leaf stoichometry:  fixed dynamic optimal
   CONTAINS
     PROCEDURE :: Init => Init_veg_config
  END type t_veg_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_config_class'

CONTAINS
  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_veg_config
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_veg_config(config)

    USE mo_jsb_namelist_iface, ONLY: process_nml, use_namelist_broadcast
    USE mo_jsb_parallel_iface, ONLY: p_bcast, process_mpi_root_id, mpi_comm

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_veg_config), INTENT(inout) :: config    !< config type for veg
    ! ---------------------------
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename,   &
                                   bc_filename
    ! INTEGER                     :: plant_functional_type_id
    ! REAL(wp)                    :: cohort_harvest_interval
    ! CHARACTER(15)               :: bnf_scheme       
    ! CHARACTER(15)               :: veg_dynamics_scheme
    ! CHARACTER(15)               :: biomass_alloc_scheme
    ! CHARACTER(15)               :: leaf_stoichom_scheme
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_veg_config'

    NAMELIST /jsb_veg_nml/       & 
      active,                       &
      bc_filename,                  &
      ic_filename!,                  &
      ! plant_functional_type_id,     &
      ! cohort_harvest_interval,      &
      ! bnf_scheme,                   & 
      ! veg_dynamics_scheme,          &
      ! biomass_alloc_scheme,         &
      ! leaf_stoichom_scheme

    CALL message(TRIM(routine), 'Starting veg configuration')

    ! Set defaults
    active                         = .TRUE.
    bc_filename                    = 'bc_land_vegetation.nc'
    ic_filename                    = 'ic_land_vegetation.nc'
    ! plant_functional_type_id       = 1
    ! cohort_harvest_interval        = 80.0_wp
    ! bnf_scheme                     = "standard"
    ! veg_dynamics_scheme            = "population"
    ! biomass_alloc_scheme           = "fixed"
    ! leaf_stoichom_scheme           = "fixed"

    ! read the namelist
    ! &jsb_veg_nml
    CALL process_nml(TRIM(config%namelist_filename), 'jsb_veg_nml', read_nml)

    IF (use_namelist_broadcast()) THEN
      CALL p_bcast(active,                   process_mpi_root_id, mpi_comm)
      CALL p_bcast(bc_filename,              process_mpi_root_id, mpi_comm)
      CALL p_bcast(ic_filename,              process_mpi_root_id, mpi_comm)
      ! CALL p_bcast(plant_functional_type_id, process_mpi_root_id, mpi_comm)
      ! CALL p_bcast(cohort_harvest_interval,  process_mpi_root_id, mpi_comm)
      ! CALL p_bcast(bnf_scheme,               process_mpi_root_id, mpi_comm)
      ! CALL p_bcast(veg_dynamics_scheme,      process_mpi_root_id, mpi_comm)
      ! CALL p_bcast(biomass_alloc_scheme,     process_mpi_root_id, mpi_comm)
      ! CALL p_bcast(leaf_stoichom_scheme,     process_mpi_root_id, mpi_comm)
    END IF

    ! pass values as read from file
    config%active                        = active
    config%bc_filename                   = bc_filename
    config%ic_filename                   = ic_filename
    ! config%pft_id                        = plant_functional_type_id
    ! config%cohort_harvest_interval       = cohort_harvest_interval
    ! config%bnf_scheme                    = TRIM(bnf_scheme)
    ! config%veg_dynamics_scheme           = TRIM(veg_dynamics_scheme)
    ! config%biomass_alloc_scheme          = biomass_alloc_scheme
    ! config%leaf_stoichom_scheme          = leaf_stoichom_scheme

  CONTAINS

    SUBROUTINE read_nml(nml_unit)
      INTEGER, INTENT(in) :: nml_unit
      READ(nml_unit, NML=jsb_veg_nml)
    END SUBROUTINE read_nml

  END SUBROUTINE Init_veg_config

#endif
END MODULE mo_veg_config_class
