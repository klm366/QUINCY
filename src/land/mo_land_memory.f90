!>
!! @brief Contains the routines for vegetation and soil memory operations 
!!
!! @author
!!  Jan Engel, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-10)
!!
MODULE mo_land_memory

  USE mo_kind,    ONLY: wp

  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------
  !> elements used across veg & soil of the model (C N P 13C 14C 15N)
  !! 
  !! this type is used for all tissues, pools, fluxes in vegetation and soil;
  !! changes made here affect the whole model
  !-----------------------------------------------------------------------------------------------------
  TYPE t_element 
     REAL(wp) :: carbon      !< total carbon concentration of pool in mol / m2
     REAL(wp) :: nitrogen    !< totol nitrogen concentration of pool in mol / m2
     REAL(wp) :: phosphorus  !< phosphorus concentration of pool in mol / m2
     REAL(wp) :: carbon13    !< 13carbon concentration of pool in mol / m2
     REAL(wp) :: carbon14    !< 14carbon concentration of pool in mol / m2
     REAL(wp) :: nitrogen15  !< 15nitrogen concentration of pool in mol / m2
  END TYPE t_element


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> addition of all elements of two tissues is returned
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION add_elements(e1,e2) RESULT(e3)
     IMPLICIT NONE
     TYPE(t_element), INTENT(in) :: e1,e2
     TYPE(t_element)             :: e3
     e3%carbon     = e1%carbon     + e2%carbon 
     e3%nitrogen   = e1%nitrogen   + e2%nitrogen 
     e3%phosphorus = e1%phosphorus + e2%phosphorus 
     e3%carbon13   = e1%carbon13   + e2%carbon13 
     e3%carbon14   = e1%carbon14   + e2%carbon14
     e3%nitrogen15 = e1%nitrogen15 + e2%nitrogen15
  END FUNCTION add_elements  

  !-----------------------------------------------------------------------------------------------------
  !> subtraction of all elements of one tissue e1 from e2 is returned
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION subtract_elements(e1,e2) RESULT(e3)
     IMPLICIT NONE
     TYPE(t_element), INTENT(in) :: e1,e2
     TYPE(t_element)             :: e3
     e3%carbon     = e1%carbon     - e2%carbon 
     e3%nitrogen   = e1%nitrogen   - e2%nitrogen 
     e3%phosphorus = e1%phosphorus - e2%phosphorus 
     e3%carbon13   = e1%carbon13   - e2%carbon13 
     e3%carbon14   = e1%carbon14   - e2%carbon14
     e3%nitrogen15 = e1%nitrogen15 - e2%nitrogen15
  END FUNCTION subtract_elements  

  !-----------------------------------------------------------------------------------------------------
  !> multiplication of all elements of two tissues is returned
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION multiply_elements(e1,e2) RESULT(e3)
     IMPLICIT NONE
     TYPE(t_element), INTENT(in) :: e1,e2
     TYPE(t_element)             :: e3
     e3%carbon     = e1%carbon     * e2%carbon 
     e3%nitrogen   = e1%nitrogen   * e2%nitrogen 
     e3%phosphorus = e1%phosphorus * e2%phosphorus 
     e3%carbon13   = e1%carbon13   * e2%carbon13 
     e3%carbon14   = e1%carbon14   * e2%carbon14
     e3%nitrogen15 = e1%nitrogen15 * e2%nitrogen15
  END FUNCTION multiply_elements

  !-----------------------------------------------------------------------------------------------------
  !> division of all elements of tissues e1 by e2 is returned
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION divide_elements(e1,e2) RESULT(e3)
     IMPLICIT NONE
     TYPE(t_element), INTENT(in) :: e1,e2
     TYPE(t_element)             :: e3
     e3%carbon     = e1%carbon     / e2%carbon 
     e3%nitrogen   = e1%nitrogen   / e2%nitrogen 
     e3%phosphorus = e1%phosphorus / e2%phosphorus 
     e3%carbon13   = e1%carbon13   / e2%carbon13 
     e3%carbon14   = e1%carbon14   / e2%carbon14
     e3%nitrogen15 = e1%nitrogen15 / e2%nitrogen15
  END FUNCTION divide_elements  


END MODULE mo_land_memory
