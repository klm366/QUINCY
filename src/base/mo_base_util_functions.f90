!>
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

  !> calculates the average across a certain period preceding to the current value;
  !! actually, the average gets rather estimated than exactly calculated using the weighted arithmetic mean 
  PURE ELEMENTAL FUNCTION moving_avrg(cav,nv,w,mask) RESULT(nav)
    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp), INTENT(in)            :: cav     !< current average value
    REAL(wp), INTENT(in)            :: nv      !< new value
    REAL(wp), INTENT(in)            :: w       !< weight of the new value
    LOGICAL, OPTIONAL, INTENT(in)   :: mask    !< bool, whether the new average is calculated from the nv, or av is used as nav
    ! 0.2 Output
    REAL(wp)                        :: nav     !< function output: the new average
    ! 0.3 Local
    LOGICAL                         :: l_mask  ! True/ False of the Optional variable mask present

    l_mask = .FALSE.
    IF(PRESENT(mask))  l_mask = mask

    IF(l_mask) THEN
      nav = cav
    ELSE
      nav = cav * (1.0_wp - w) + nv * w  ! weighted mean
    END IF

  END FUNCTION moving_avrg


END MODULE mo_base_util_functions
