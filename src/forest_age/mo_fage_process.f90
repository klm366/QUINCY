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

   PUBLIC :: ageing !TODO-JN: merge_cohorts update_age?

   CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_process'

 CONTAINS   


  SUBROUTINE ageing(nidx, nage, is_newyear,       &   !in
                    fract_per_age,                &   !inout
                    mean_age                      &   !inout
                    )

    USE mo_fage_util,             ONLY: get_mean_age

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    INTEGER,  INTENT(in) :: &
      & nidx, nage
    LOGICAL,   INTENT(in)    :: is_newyear
    REAL(wp), INTENT(inout), DIMENSION(nidx, nage) :: &
      & fract_per_age
    REAL(wp), INTENT(inout), DIMENSION(nidx) :: &
      & mean_age

    !-----------------------------------------------------------------------
    !  LOCAL
    INTEGER :: i

    !-----------------------------------------------------------------------
    ! CONTENT

    ! --- update age

    !TODO-JN: check comment
    ! Note that updating is done globally at the same time step, i.e. for different longitudes the updating happens at different
    ! local times.

    IF(is_newyear) THEN  ! year has changed - i.e. age - i.e. shift fractions!

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

    END IF

  END SUBROUTINE ageing


#endif
END MODULE mo_fage_process
