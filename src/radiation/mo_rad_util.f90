!>
!! @brief Contains the routines for the radiation processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-27)
!! move functions here from mo_jsb_rad_process by Jan Engel MPI-BGC (20230404)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_rad_util
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Has_minimal_radiation

  ! ======================================================================================================= !
  !> is day/night based on radiation
  !>
  INTERFACE Has_minimal_radiation
    MODULE PROCEDURE Has_minimal_radiation_sw
    MODULE PROCEDURE Has_minimal_radiation_swvis_swnir
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_util'

CONTAINS

  ! ======================================================================================================= !
  !> is day/night based on shortwave radiation
  !>
  ELEMENTAL FUNCTION Has_minimal_radiation_sw(sw) RESULT(l_day)
    USE mo_rad_parameters,    ONLY: MinRadiation ! Minimum amount of radiation for albedo calculations

    REAL(wp), INTENT(in) :: sw
    LOGICAL              :: l_day

    !$claw sca routine
    l_day = sw > MinRadiation
  END FUNCTION Has_minimal_radiation_sw

  ! ======================================================================================================= !
  !> is day/night based on shortwave "vis + nir" radiation
  !>
  ELEMENTAL FUNCTION Has_minimal_radiation_swvis_swnir(swvis, swnir) RESULT(l_day)
    USE mo_rad_parameters,    ONLY: MinRadiation ! Minimum amount of radiation for albedo calculations

    REAL(wp), INTENT(in) :: swvis
    REAL(wp), INTENT(in) :: swnir
    LOGICAL              :: l_day

    !$claw sca routine
    l_day = swvis + swnir > MinRadiation
  END FUNCTION Has_minimal_radiation_swvis_swnir

#endif
END MODULE mo_rad_util
