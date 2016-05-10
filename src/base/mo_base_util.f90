!>
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
 
!  PRIVATE
 
!  PUBLIC :: test_real_for_NaN, test_int_for_NaN



CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> test if the passed REAL variable is NaN
  !! 
  !! returns TRUE if the variable is NaN; IF (X /= X) THEN X is an NaN beacuse NaN is equal to nothing including itself 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION test_real_for_NaN(test_var) RESULT(l_result)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in)      :: test_var    !< variable to be tested
    LOGICAL                   :: l_result

    l_result = .FALSE.
    IF(test_var /= test_var) l_result = .TRUE.
  END FUNCTION test_real_for_NaN

  !-----------------------------------------------------------------------------------------------------
  !> test if the passed INTEGER variable is NaN
  !!
  !! returns TRUE if the variable is NaN; IF (X /= X) THEN X is an NaN beacuse NaN is equal to nothing including itself 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION test_int_for_NaN(test_var) RESULT(l_result)
    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in)       :: test_var    !< variable to be tested
    LOGICAL                   :: l_result

    l_result = .FALSE.
    IF(test_var /= test_var) l_result = .TRUE.
  END FUNCTION test_int_for_NaN



END MODULE mo_base_util
