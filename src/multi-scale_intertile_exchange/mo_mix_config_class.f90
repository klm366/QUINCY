!> Contains structures and methods for the multi-scale intertile exchange config
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
MODULE mo_mix_config_class
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
  PUBLIC :: t_mix_config

  TYPE, EXTENDS(t_jsb_config) :: t_mix_config
     LOGICAL                  :: &
       & l_hydro_mix               !< enable lateral transport of water between grid tiles

   CONTAINS
     PROCEDURE :: Init => Init_mix_config
  END type t_mix_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_mix_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize multi-scale intertile exchange process
  !!  
  SUBROUTINE Init_mix_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_mix_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: l_hydro_mix

    NAMELIST /jsb_mix_nml/       &
         active,                 &
         l_hydro_mix

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_mix_config'

    CALL message(TRIM(routine), 'Starting multi-scale intertile exchange configuration')

    ! Set defaults
    active              = .TRUE.
    l_hydro_mix         = .FALSE.

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_mix_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_mix_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%l_hydro_mix         = l_hydro_mix
    IF (l_hydro_mix) THEN
      CALL message(TRIM(routine), 'Lateral water flow between tiles is enabled')
    ELSE
      CALL message(TRIM(routine), 'Lateral water flow between tiles is disabled')
    END IF

  END SUBROUTINE Init_mix_config

#endif
END MODULE mo_mix_config_class
