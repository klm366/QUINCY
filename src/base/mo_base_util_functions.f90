!!
!! @brief generic functions useful across the whole model
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-02)
!!
MODULE mo_base_util_functions

  USE mo_kind, ONLY : wp

  IMPLICIT NONE


CONTAINS

  FUNCTION moving_avrg(cav,nv,w,l) RESULT(nav)
  ! calculates the average across a certain period preceding to the current value
  ! actually, the average gets rather estimated than exactly calculated using the weighted arithmetic mean 
  ! cav: current average value
  ! nv:  new value
  ! w:   weight of the new value
  ! l:   bool, whether the new average is calculated from the nv, or av is used as nav
  ! nav: output of the function, the new average
    IMPLICIT NONE
    REAL(wp), INTENT(in)    :: cav,nv,w  ! input
    LOGICAL, INTENT(in) 		:: l         ! input
    REAL(wp)                :: nav       ! output

    IF(l) THEN
      nav = cav * (1-w) + nv * w         ! weighted mean
    ELSE
      nav = cav
    END IF
  END FUNCTION moving_avrg


END MODULE mo_base_util_functions
