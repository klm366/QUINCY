! Meri's to-be-deleted notes:
! Done so far: 
! - copied the flcc process and swapped all flcc -> hlcc
!
!
!> hlcc (lcc due to harvest) config
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
MODULE mo_hlcc_config_class
  #ifndef __NO_JSBACH__
  
    USE mo_exception,         ONLY: message
    USE mo_jsb_control,       ONLY: debug_on
    USE mo_io_units,          ONLY: filename_max
    USE mo_jsb_config_class,  ONLY: t_jsb_config
  
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: t_hlcc_config
  
    TYPE, EXTENDS(t_jsb_config) :: t_hlcc_config
        LOGICAL :: separate
      CONTAINS
        PROCEDURE :: Init => Init_hlcc_config
    END type t_hlcc_config
  
    CHARACTER(len=*), PARAMETER :: modname = 'mo_hlcc_config_class'
  
  CONTAINS
  
    ! ====================================================================================================== !
    !
    !> Initialize hlcc process
    !
    ! -------------------------------------------------------------------------------------------------------
    SUBROUTINE Init_hlcc_config(config)
  
      USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
  
      ! -------------------------------------------------------------------------------------------------- !
      CLASS(t_hlcc_config), INTENT(inout) :: config !<Configuration type of process (t_hlcc_config)
      ! -------------------------------------------------------------------------------------------------- !
      LOGICAL                       :: active
      CHARACTER(len=filename_max)   :: ic_filename, bc_filename
      LOGICAL                       :: separate !< whether burned tiles should be treated separately or at once
  
      NAMELIST /jsb_hlcc_nml/      &
           active,                 &
           ic_filename,            &
           bc_filename,            &
           separate
  
      INTEGER :: nml_handler, nml_unit, istat
  
      CHARACTER(len=*), PARAMETER :: routine = modname//':Init_hlcc_config'
      ! -------------------------------------------------------------------------------------------------- !
  
      IF (debug_on()) CALL message(TRIM(routine), 'Starting hlcc configuration')
  
      ! Set defaults
      separate            = .TRUE.
      active              = .FALSE.
      bc_filename         = 'bc_land_hlcc.nc'
      ic_filename         = 'ic_land_hlcc.nc'
  
      nml_handler = open_nml(TRIM(config%namelist_filename))
  
      nml_unit = position_nml('jsb_hlcc_nml', nml_handler, STATUS=istat)
      IF (istat == POSITIONED) READ(nml_unit, jsb_hlcc_nml)
  
      CALL close_nml(nml_handler)
  
      config%separate            = separate
      config%active              = active
      config%ic_filename         = ic_filename
      config%bc_filename         = bc_filename
  
    END SUBROUTINE Init_hlcc_config
  
  #endif
  END MODULE mo_hlcc_config_class