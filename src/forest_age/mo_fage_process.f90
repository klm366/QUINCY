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

   PUBLIC :: update_fract_per_age_upon_ageing, update_fract_per_age_equally_upon_disturbance, &
     & update_fract_per_age_and_determine_area_to_harvest_from_cohort, correct_and_synchronise_fracts, &
     & weighted_avg_var_upon_area_change

   !TODO-JN: merge_cohorts update_age?

   CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_process'

 CONTAINS   


  SUBROUTINE update_fract_per_age_upon_ageing( &
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

  END SUBROUTINE update_fract_per_age_upon_ageing


  !TODO-JN: make faster?
  SUBROUTINE correct_and_synchronise_fracts( &
                    nidx, nage, ncohorts, iThisAgeCohort,  &   !in
                    pftid, cohort_max_age,                 &   !in
                    fract_per_age,                         &   !inout
                    fract                                  &   !inout
                    )

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    INTEGER,  INTENT(in)    ::  nidx, nage, ncohorts, iThisAgeCohort, pftid
    REAL(wp), INTENT(inout), DIMENSION(nidx, ncohorts) :: cohort_max_age
    REAL(wp), INTENT(inout), DIMENSION(nidx, nage) :: fract_per_age
    REAL(wp), INTENT(inout), DIMENSION(nidx) :: fract

    !-----------------------------------------------------------------------
    !  LOCAL
    INTEGER :: j, k, thisMinAge, thisMaxAge
    REAL(wp) :: thisAreaSum, stillMissing, distribute, difference

    !-----------------------------------------------------------------------
    ! CONTENT

    ! TODO-JN-COHORTS: could be much simpler, if age cohort boundaries will not be dynamic
    DO j=1,nidx

      IF (iThisAgeCohort .GT. 1) THEN
        thisMinAge = int(cohort_max_age(j,iThisAgeCohort-1)+1)
        thisMaxAge = int(cohort_max_age(j,iThisAgeCohort))
        thisAreaSum = SUM(fract_per_age(j,thisMinAge:thisMaxAge))

        IF(thisAreaSum - fract(j) .GT. 1.e-14_wp) THEN !TODO-JN make 1.e-14_wp a constant somewhere
          WRITE (message_text,*) ' fract (', fract(j), &
             & ') and fage sum (', thisAreaSum , ') missmatch for age cohort ', iThisAgeCohort, &
             & '-- pft ', pftid
          CALL finish ('correct_and_synchronise_fracts', message_text)
        ENDIF

        fract(j) = thisAreaSum
      ELSE
        fract_per_age(j,1) = fract(j)
      ENDIF
    ENDDO
  END SUBROUTINE correct_and_synchronise_fracts


  !TODO-JN: make faster?
  SUBROUTINE update_fract_per_age_equally_upon_disturbance( &
                    nidx, nage, ncohorts, iThisAgeCohort,  &   !in
                    disturbed_area,                        &   !in
                    cohort_max_age,                        &   !in
                    fract_per_age                          &   !inout
                    )

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    INTEGER,  INTENT(in)    ::  nidx, nage, ncohorts, iThisAgeCohort
    REAL(wp), INTENT(inout), DIMENSION(nidx) :: disturbed_area
    REAL(wp), INTENT(inout), DIMENSION(nidx, ncohorts) :: cohort_max_age
    REAL(wp), INTENT(inout), DIMENSION(nidx, nage) :: fract_per_age

    !-----------------------------------------------------------------------
    !  LOCAL
    INTEGER :: j, k, thisMinAge, thisMaxAge
    REAL(wp) :: thisAreaSum

    !-----------------------------------------------------------------------
    ! CONTENT
    IF(iThisAgeCohort .EQ. 1) RETURN !area is only distributed to the first child - not from

    ! damaged area is assigned to zero age and taken evenly from all ages in this age cohort
    ! TODO-JN-COHORTS: could be much simpler, if age cohort boundaries will not be dynamic
    DO j=1,nidx
      thisMinAge = int(cohort_max_age(j,iThisAgeCohort-1)+1)
      thisMaxAge = int(cohort_max_age(j,iThisAgeCohort))

      thisAreaSum = SUM(fract_per_age(j,thisMinAge:thisMaxAge))

      IF (thisAreaSum .GT. 0._wp) THEN
        !Assertion: sum of available area needs to be larger than disturbed area
        IF(disturbed_area(j) .GT. thisAreaSum) THEN
          WRITE (message_text,*) 'more area disturbed (', disturbed_area(j), &
            & ') than available (', thisAreaSum , ') at age cohort ', iThisAgeCohort
          CALL finish ('update_fract_per_age_equally_upon_disturbance', message_text)
        ENDIF

        fract_per_age(j,1) = fract_per_age(j,1) + disturbed_area(j)

        DO k=thisMinAge,thisMaxAge
          fract_per_age(j,k) = fract_per_age(j,k) - (disturbed_area(j) * (fract_per_age(j,k)/thisAreaSum))
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE update_fract_per_age_equally_upon_disturbance


  !TODO-JN: make faster?
  SUBROUTINE update_fract_per_age_and_determine_area_to_harvest_from_cohort( &
                    harvest_kind, nidx, nage, &   !in
                    this_cohort_max_age,                             &   !in
                    prev_cohort_max_age,                             &   !in
                    fract_per_age,                                   &   !inout
                    harvested_area,                                  &   !inout
                    to_be_harvested_area                             &   !inout
                    )

    USE mo_harvest_constants,    ONLY: HARVEST_OLDEST_TREES

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    INTEGER,  INTENT(in)    ::  harvest_kind, nidx, nage
    REAL(wp), INTENT(inout), DIMENSION(nidx) :: harvested_area, to_be_harvested_area
    REAL(wp), INTENT(inout), DIMENSION(nidx) :: this_cohort_max_age, prev_cohort_max_age
    REAL(wp), INTENT(inout), DIMENSION(nidx, nage) :: fract_per_age

    !-----------------------------------------------------------------------
    !  LOCAL
    INTEGER :: j, k, thisMinAge, thisMaxAge
    REAL(wp) :: thisHarvest

    !-----------------------------------------------------------------------
    ! CONTENT

    ! TODO-JN-COHORTS: could be much simpler, if age cohort boundaries will not be dynamic
    DO j=1,nidx
      thisMinAge = int(prev_cohort_max_age(j)+1)
      thisMaxAge = int(this_cohort_max_age(j))

      IF (harvest_kind .EQ. HARVEST_OLDEST_TREES) THEN
        DO k=thisMaxAge,thisMinAge,-1
          IF(to_be_harvested_area(j) .GT. 0._wp) THEN
            thisHarvest = MIN(fract_per_age(j,k), to_be_harvested_area(j))
            fract_per_age(j,k) = fract_per_age(j,k) - thisHarvest
            fract_per_age(j,1) = fract_per_age(j,1) + thisHarvest
            to_be_harvested_area(j) = to_be_harvested_area(j) - thisHarvest
            harvested_area(j) = harvested_area(j) + thisHarvest
          ELSE
            EXIT
          ENDIF
        ENDDO
      ELSE
        CALL finish('update_fract_per_age_and_determine_area_to_harvest_from_cohort','Not implemented harvest kind')
      ENDIF
    ENDDO

  END SUBROUTINE update_fract_per_age_and_determine_area_to_harvest_from_cohort


  ELEMENTAL SUBROUTINE weighted_avg_var_upon_area_change( & 
                                          & old_fract,     & ! <- in
                                          & delta_fract,   & ! <- in
                                          & new_var,       & ! <- in
                                          & old_var)         ! <- inout

     !-----------------------------------------------------------------------
     !  ARGUMENTS
     REAL(wp), intent(in)    :: old_fract, delta_fract, new_var
     REAL(wp), intent(inout) :: old_var

     IF (delta_fract .GT. 0._wp) THEN 
       old_var = ((old_var * old_fract) + (new_var * delta_fract)) / (old_fract + delta_fract)
     ENDIF

     RETURN

  END SUBROUTINE weighted_avg_var_upon_area_change


#endif
END MODULE mo_fage_process
