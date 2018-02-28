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
   USE mo_exception, ONLY: message, message_text, finish
 
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: update_frac_per_age_upon_ageing, update_frac_per_age_equally_upon_disturbance 

   !TODO-JN: merge_cohorts update_age?

   CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_process'

 CONTAINS   


  SUBROUTINE update_frac_per_age_upon_ageing( &
                    nidx, nage,               &   !in
                    fract_per_age,            &   !inout
                    mean_age                  &   !inout
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

  END SUBROUTINE update_frac_per_age_upon_ageing


  !TODO-JN: make faster
  SUBROUTINE update_frac_per_age_equally_upon_disturbance( &
                    nidx, nage, ncohorts, iThisCohort,     &   !in
                    disturbed_area,                        &   !in
                    cohort_max_age,                        &   !in
                    fract_per_age,                         &   !inout
                    mean_age                               &   !inout
                    )

    USE mo_fage_util,             ONLY: get_mean_age

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    INTEGER,  INTENT(in)    ::  nidx, nage, ncohorts, iThisCohort
    REAL(wp), INTENT(inout), DIMENSION(nidx) :: disturbed_area
    REAL(wp), INTENT(inout), DIMENSION(nidx, ncohorts) :: cohort_max_age
    REAL(wp), INTENT(inout), DIMENSION(nidx, nage) :: fract_per_age
    REAL(wp), INTENT(inout), DIMENSION(nidx) :: mean_age

    !-----------------------------------------------------------------------
    !  LOCAL
    INTEGER :: j, k, thisMinAge, thisMaxAge
    REAL(wp) :: thisAreaSum

    !-----------------------------------------------------------------------
    ! CONTENT

    ! damaged area is assigned to zero age and taken evenly from all ages in this cohort
    ! TODO-JN-COHORTS: could be much simpler, if cohort boundaries will not be dynamic
    DO j=1,nidx

      thisMinAge = int(cohort_max_age(j,iThisCohort-1)+1)
      thisMaxAge = int(cohort_max_age(j,iThisCohort))

      thisAreaSum = SUM(fract_per_age(j,thisMinAge:thisMaxAge))

      IF (thisAreaSum .GT. 0._wp) THEN
        !Assertion: sum of available area needs to be larger than disturbed area
        IF(disturbed_area(j) .GT. thisAreaSum) THEN
          WRITE (message_text,*) 'more area disturbed (', disturbed_area(j), &
            & ') than available (', thisAreaSum , ') at cohort ', iThisCohort
          CALL finish ('update_frac_per_age_equally_upon_disturbance', message_text)
        ENDIF

        fract_per_age(j,1) = fract_per_age(j,1) + disturbed_area(j)

        DO k=thisMinAge,thisMaxAge
          fract_per_age(j,k) = fract_per_age(j,k) - (disturbed_area(j) * (fract_per_age(j,k)/thisAreaSum))
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE update_frac_per_age_equally_upon_disturbance


#endif
END MODULE mo_fage_process
