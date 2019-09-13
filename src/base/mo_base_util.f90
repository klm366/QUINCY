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
!! @brief generic functions and specific TYPE's useful across the whole model
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-02)
!!
MODULE mo_base_util

  USE mo_kind,            ONLY: wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: test_for_NaN_Inf

  !-----------------------------------------------------------------------------------------------------
  !> returns TRUE if the passed variable is NaN or Infinity
  !!
  !! accepts REAL and INT
  !-----------------------------------------------------------------------------------------------------
  INTERFACE test_for_NaN_Inf
    MODULE PROCEDURE test_int_for_NaN_and_Infinity
    MODULE PROCEDURE test_real_for_NaN_and_Infinity
  END INTERFACE test_for_NaN_Inf

  CHARACTER(len=*), PARAMETER :: modname = 'mo_base_util'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> test if the passed INTEGER variable is NaN or Infinity
  !!
  !! returns TRUE if the variable is NaN; IF (X /= X) THEN X is an NaN beacuse NaN is equal to nothing including itself \n
  !! returns TRUE if the variable is Infinity; in this case the subtraction of 1 does not change the value
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION test_int_for_NaN_and_Infinity(test_var) RESULT(l_result)
    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in)       :: test_var    !< INT variable to be tested
    LOGICAL                   :: l_result
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':test_int_for_NaN_and_Infinity'

    l_result = .FALSE.
    IF(test_var /= test_var) l_result = .TRUE.
    IF(.NOT. l_result) THEN
      IF(test_var - 1 .eq. test_var) l_result = .TRUE.
    END IF
  END FUNCTION test_int_for_NaN_and_Infinity

  !-----------------------------------------------------------------------------------------------------
  !> test if the passed REAL variable is NaN or Infinity
  !! 
  !! returns TRUE if the variable is NaN; IF (X /= X) THEN X is an NaN beacuse NaN is equal to nothing including itself \n
  !! returns TRUE if the variable is Infinity; in this case the subtraction of 1 does not change the value
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION test_real_for_NaN_and_Infinity(test_var) RESULT(l_result)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in)      :: test_var    !< REAL variable to be tested
    LOGICAL                   :: l_result
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':test_real_for_NaN_and_Infinity'

    l_result = .FALSE.
    IF(test_var /= test_var) l_result = .TRUE.
    IF(.NOT. l_result) THEN
      IF(test_var - 1.0_wp .eq. test_var) l_result = .TRUE.
    END IF
  END FUNCTION test_real_for_NaN_and_Infinity



END MODULE mo_base_util
