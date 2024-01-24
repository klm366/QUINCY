!> Summary: JSBACH & QUINCY helper routines for radiation
!>
!>#### ICON-Land
!>
!> ------------------------------------------------------------------------
!> Copyright (C) 2004-2024, MPI-M, MPI-BGC, DWD, DKRZ, KIT, ETH, MeteoSwiss  
!> Contact information: <https://icon-model.org>
!>
!> See AUTHORS.TXT for a list of authors  
!> See LICENSES/ for license information  
!> SPDX-License-Identifier: BSD-3-Clause
!> ------------------------------------------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### various helper routines for the radiation process
!>
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
