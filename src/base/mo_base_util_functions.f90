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
 
  PRIVATE
 
  PUBLIC :: moving_avrg

CONTAINS

PURE ELEMENTAL FUNCTION moving_avrg(cav,nv,w,mask) RESULT(nav)
  ! calculates the average across a certain period preceding to the current value
  ! actually, the average gets rather estimated than exactly calculated using the weighted arithmetic mean 
  ! cav: current average value
  ! nv:  new value
  ! w:   weight of the new value
  ! mask:   bool, whether the new average is calculated from the nv, or av is used as nav !! this is an OPTIONAL argument
  ! nav: output of the function, the new average
    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp), INTENT(in)            :: cav,nv,w  ! input
    LOGICAL, OPTIONAL, INTENT(in)   :: mask    ! input
    ! 0.3 Outout
    REAL(wp)                        :: nav       ! output
    ! 0.4 Local
    LOGICAL                         :: l_mask

    l_mask = .FALSE.
    IF(PRESENT(mask))  l_mask = mask

    IF(l_mask) THEN
      nav = cav
    ELSE
      nav = cav * (1.0_wp - w) + nv * w  ! weighted mean
    END IF

END FUNCTION moving_avrg


END MODULE mo_base_util_functions
