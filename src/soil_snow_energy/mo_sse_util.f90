!>
!! @brief Contains the utilities for the soil and snow energy process
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_sse_util
#ifndef __NO_JSBACH__

  USE mo_kind,          ONLY: wp
  USE mo_exception,     ONLY: finish, message, message_text
  USE mo_jsb_control,   ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calc_qmax_texture

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_util'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION calc_qmax_texture
  !!
  !! Output: sorption capacity for organic material
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_qmax_texture(qmax_fine_particle, silt_sl, clay_sl) RESULT(qmax_texture)


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: qmax_fine_particle    !< maximum sorption capacity of fine soil particle [mol / kg fine particle (silt or clay)]
    REAL(wp), INTENT(in) :: silt_sl               !< silt content of the mineral soil [kg silt / kg mineral soil]
    REAL(wp), INTENT(in) :: clay_sl               !< clay content of the mineral soil [kg clay / kg mineral soil]
    REAL(wp)             :: qmax_texture          !< maximum sorption capacity of the soil [mol / kg]
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_qmax_texture'

    ! maximum sorption capacity per kg mineral soil (mol/kg)
    qmax_texture = qmax_fine_particle * (silt_sl + clay_sl)

  END FUNCTION calc_qmax_texture

#endif
END MODULE mo_sse_util
