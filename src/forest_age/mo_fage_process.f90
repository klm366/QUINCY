!>
!! @brief Contains the routines for the forest-age process.
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2017-11-21)
!!
MODULE mo_fage_process
#ifndef __NO_JSBACH__
 
   USE mo_kind,      ONLY: wp
   USE mo_exception, ONLY: message
 
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: update_frac_per_age_tracking !TODO-JN: merge_cohorts update_age?

   CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_process'

 CONTAINS   


  SUBROUTINE update_frac_per_age_tracking(&
                    nidx, nage,                   &   !in
                    fract_per_age,                &   !inout
                    mean_age                      &   !inout
                    )

    USE mo_fage_util,             ONLY: get_mean_age

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    INTEGER,  INTENT(in)    ::  nidx, nage
    REAL(wp), INTENT(inout), DIMENSION(nidx, nage) :: fract_per_age
    REAL(wp), INTENT(inout), DIMENSION(nidx) :: mean_age

    !-----------------------------------------------------------------------
    !  LOCAL
    INTEGER :: i

    !-----------------------------------------------------------------------
    ! CONTENT

    ! --- update age
    IF (nage > 1) THEN
      ! oldest age accumulates
      fract_per_age (:,nage) = fract_per_age (:,nage) + fract_per_age (:,nage-1)

      DO i=nage-1,2,-1
        fract_per_age (:,i) = fract_per_age (:,i-1)
      END DO

      ! youngest is emptied
      fract_per_age (:,1) = 0._wp

      !recalculate mean age
      mean_age (:) = get_mean_age(fract_per_age)
    END IF

  END SUBROUTINE update_frac_per_age_tracking


#endif
END MODULE mo_fage_process
