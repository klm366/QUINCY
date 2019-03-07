!>
!! @brief Contains utilities for the forest-age process
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2017-11-21)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
MODULE mo_fage_util
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_mean_age

  INTERFACE get_mean_age
    MODULE PROCEDURE get_mean_age_2d
    MODULE PROCEDURE get_mean_age_3d
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_util'

CONTAINS

  FUNCTION get_mean_age_2d(fract_per_age) RESULT(mean_age)

    REAL(wp), INTENT(in) :: &
      & fract_per_age (:,:)     !< fraction for each age (nidx,z)
    REAL(wp) :: mean_age(SIZE(fract_per_age,1))  !< mean age

    INTEGER :: i,n

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_mean_age_2d'

    n  = SIZE(fract_per_age,2)       ! max tracked age

    !TODO-JN: in case of absolute fractions this would need to be normalised
    mean_age (:) = 0._wp
    DO i=1,n
      mean_age (:) = mean_age (:) + fract_per_age (:,i) * i
    END DO

  END FUNCTION get_mean_age_2d

  FUNCTION get_mean_age_3d(fract_per_age) RESULT(mean_age)

    REAL(wp), INTENT(in) :: &
      & fract_per_age (:,:,:)     !< fraction for each age (x,y,z)
    REAL(wp) :: mean_age(SIZE(fract_per_age,1),SIZE(fract_per_age,2))  !< mean age

    INTEGER :: i,n

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_mean_age_2d'

    n  = SIZE(fract_per_age,3)       ! max tracked age

    !TODO-JN: in case of absolute fractions this would need to be normalised
    mean_age (:,:) = 0._wp
    DO i=1,n
      mean_age (:,:) = mean_age (:,:) + fract_per_age (:,:,i) * i
    END DO

  END FUNCTION get_mean_age_3d


#endif
END MODULE mo_fage_util
