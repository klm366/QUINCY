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
  CONTAINS
    PRIVATE
    PROCEDURE :: add_elements_type_type
    PROCEDURE :: add_elements_type_scalar
    PROCEDURE, PASS(this2) :: add_elements_scalar_type
    GENERIC, PUBLIC :: OPERATOR(+) => add_elements_type_type, add_elements_type_scalar, add_elements_scalar_type
    PROCEDURE :: subtract_elements_type_type
    PROCEDURE :: subtract_elements_type_scalar
    PROCEDURE, PASS(this2) :: subtract_elements_scalar_type
    GENERIC, PUBLIC :: OPERATOR(-) => subtract_elements_type_type, subtract_elements_type_scalar, subtract_elements_scalar_type
    PROCEDURE :: multiply_elements_type_type
    PROCEDURE :: multiply_elements_type_scalar
    PROCEDURE, PASS(this2) :: multiply_elements_scalar_type
    GENERIC, PUBLIC :: OPERATOR(*) => multiply_elements_type_type, multiply_elements_type_scalar, multiply_elements_scalar_type
    PROCEDURE :: divide_elements_type_type
    PROCEDURE :: divide_elements_type_scalar
    PROCEDURE, PASS(this2) :: divide_elements_scalar_type
    GENERIC, PUBLIC :: OPERATOR(/) => divide_elements_type_type, divide_elements_type_scalar, divide_elements_scalar_type
  END TYPE t_element

  TYPE(t_element) :: element_zero   !< all variables are zero, is used for setting other elements zero frequently


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> 1.0 addition of all elements of two tissues: type + type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION add_elements_type_type(this,t2) RESULT(t3)
    IMPLICIT NONE
    CLASS(t_element), INTENT(in) :: this, t2
    TYPE(t_element)              :: t3
    t3%carbon     = this%carbon     + t2%carbon 
    t3%nitrogen   = this%nitrogen   + t2%nitrogen 
    t3%phosphorus = this%phosphorus + t2%phosphorus 
    t3%carbon13   = this%carbon13   + t2%carbon13 
    t3%carbon14   = this%carbon14   + t2%carbon14
    t3%nitrogen15 = this%nitrogen15 + t2%nitrogen15
  END FUNCTION add_elements_type_type  

  !-----------------------------------------------------------------------------------------------------
  !> 1.1 addition of all elements of a tissues: type + scalar
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION add_elements_type_scalar(this,s2) RESULT(t3)
    IMPLICIT NONE
    CLASS(t_element), INTENT(in) :: this
    REAL, INTENT(in)             :: s2
    TYPE(t_element)              :: t3
    t3%carbon     = this%carbon     + s2
    t3%nitrogen   = this%nitrogen   + s2
    t3%phosphorus = this%phosphorus + s2
    t3%carbon13   = this%carbon13   + s2
    t3%carbon14   = this%carbon14   + s2
    t3%nitrogen15 = this%nitrogen15 + s2
  END FUNCTION add_elements_type_scalar  

  !-----------------------------------------------------------------------------------------------------
  !> 1.2 addition of all elements of a tissues: scalar + type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION add_elements_scalar_type(s1,this2) RESULT(t3)
    IMPLICIT NONE
    REAL, INTENT(in)             :: s1
    CLASS(t_element), INTENT(in) :: this2
    TYPE(t_element)              :: t3
    t3%carbon     = s1 + this2%carbon
    t3%nitrogen   = s1 + this2%nitrogen
    t3%phosphorus = s1 + this2%phosphorus
    t3%carbon13   = s1 + this2%carbon13
    t3%carbon14   = s1 + this2%carbon14
    t3%nitrogen15 = s1 + this2%nitrogen15
  END FUNCTION add_elements_scalar_type  

  !-----------------------------------------------------------------------------------------------------
  !> 2.0 subtraction of all elements of two tissues: type - type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION subtract_elements_type_type(this,t2) RESULT(t3)
    IMPLICIT NONE
    CLASS(t_element), INTENT(in) :: this, t2
    TYPE(t_element)              :: t3
    t3%carbon     = this%carbon     - t2%carbon 
    t3%nitrogen   = this%nitrogen   - t2%nitrogen 
    t3%phosphorus = this%phosphorus - t2%phosphorus 
    t3%carbon13   = this%carbon13   - t2%carbon13 
    t3%carbon14   = this%carbon14   - t2%carbon14
    t3%nitrogen15 = this%nitrogen15 - t2%nitrogen15
  END FUNCTION subtract_elements_type_type  

  !-----------------------------------------------------------------------------------------------------
  !> 2.1 subtraction of all elements of a tissues: type - scalar
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION subtract_elements_type_scalar(this,s2) RESULT(t3)
    IMPLICIT NONE
    CLASS(t_element), INTENT(in) :: this
    REAL, INTENT(in)             :: s2
    TYPE(t_element)              :: t3
    t3%carbon     = this%carbon     - s2
    t3%nitrogen   = this%nitrogen   - s2
    t3%phosphorus = this%phosphorus - s2
    t3%carbon13   = this%carbon13   - s2
    t3%carbon14   = this%carbon14   - s2
    t3%nitrogen15 = this%nitrogen15 - s2
  END FUNCTION subtract_elements_type_scalar  

  !-----------------------------------------------------------------------------------------------------
  !> 2.2 subtraction of all elements of a tissues: scalar - type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION subtract_elements_scalar_type(s1,this2) RESULT(t3)
    IMPLICIT NONE
    REAL, INTENT(in)             :: s1
    CLASS(t_element), INTENT(in) :: this2
    TYPE(t_element)              :: t3
    t3%carbon     = s1 - this2%carbon
    t3%nitrogen   = s1 - this2%nitrogen
    t3%phosphorus = s1 - this2%phosphorus
    t3%carbon13   = s1 - this2%carbon13
    t3%carbon14   = s1 - this2%carbon14
    t3%nitrogen15 = s1 - this2%nitrogen15
  END FUNCTION subtract_elements_scalar_type  

  !-----------------------------------------------------------------------------------------------------
  !> 3.0 multiplication of all elements of two tissues: type * type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION multiply_elements_type_type(this,t2) RESULT(t3)
    IMPLICIT NONE
    CLASS(t_element), INTENT(in) :: this, t2
    TYPE(t_element)              :: t3
    t3%carbon     = this%carbon     * t2%carbon 
    t3%nitrogen   = this%nitrogen   * t2%nitrogen 
    t3%phosphorus = this%phosphorus * t2%phosphorus 
    t3%carbon13   = this%carbon13   * t2%carbon13 
    t3%carbon14   = this%carbon14   * t2%carbon14
    t3%nitrogen15 = this%nitrogen15 * t2%nitrogen15
  END FUNCTION multiply_elements_type_type  

  !-----------------------------------------------------------------------------------------------------
  !> 3.1 multiplication of all elements of a tissues: type * scalar
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION multiply_elements_type_scalar(this,s2) RESULT(t3)
    IMPLICIT NONE
    CLASS(t_element), INTENT(in) :: this
    REAL, INTENT(in)             :: s2
    TYPE(t_element)              :: t3
    t3%carbon     = this%carbon     * s2
    t3%nitrogen   = this%nitrogen   * s2
    t3%phosphorus = this%phosphorus * s2
    t3%carbon13   = this%carbon13   * s2
    t3%carbon14   = this%carbon14   * s2
    t3%nitrogen15 = this%nitrogen15 * s2
  END FUNCTION multiply_elements_type_scalar  

  !-----------------------------------------------------------------------------------------------------
  !> 3.2 multiplication of all elements of a tissues: scalar * type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION multiply_elements_scalar_type(s1,this2) RESULT(t3)
    IMPLICIT NONE
    REAL, INTENT(in)             :: s1
    CLASS(t_element), INTENT(in) :: this2
    TYPE(t_element)              :: t3
    t3%carbon     = s1 * this2%carbon
    t3%nitrogen   = s1 * this2%nitrogen
    t3%phosphorus = s1 * this2%phosphorus
    t3%carbon13   = s1 * this2%carbon13
    t3%carbon14   = s1 * this2%carbon14
    t3%nitrogen15 = s1 * this2%nitrogen15
  END FUNCTION multiply_elements_scalar_type  

  !-----------------------------------------------------------------------------------------------------
  !> 4.0 division of all elements of two tissues: type / type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION divide_elements_type_type(this,t2) RESULT(t3)
    IMPLICIT NONE
    CLASS(t_element), INTENT(in) :: this, t2
    TYPE(t_element)              :: t3
    t3%carbon     = this%carbon     / t2%carbon 
    t3%nitrogen   = this%nitrogen   / t2%nitrogen 
    t3%phosphorus = this%phosphorus / t2%phosphorus 
    t3%carbon13   = this%carbon13   / t2%carbon13 
    t3%carbon14   = this%carbon14   / t2%carbon14
    t3%nitrogen15 = this%nitrogen15 / t2%nitrogen15
  END FUNCTION divide_elements_type_type  

  !-----------------------------------------------------------------------------------------------------
  !> 4.1 division of all elements of a tissues: type / scalar
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION divide_elements_type_scalar(this,s2) RESULT(t3)
    IMPLICIT NONE
    CLASS(t_element), INTENT(in) :: this
    REAL, INTENT(in)             :: s2
    TYPE(t_element)              :: t3
    t3%carbon     = this%carbon     / s2
    t3%nitrogen   = this%nitrogen   / s2
    t3%phosphorus = this%phosphorus / s2
    t3%carbon13   = this%carbon13   / s2
    t3%carbon14   = this%carbon14   / s2
    t3%nitrogen15 = this%nitrogen15 / s2
  END FUNCTION divide_elements_type_scalar  

  !-----------------------------------------------------------------------------------------------------
  !> 4.2 division of all elements of a tissues: scalar / type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION divide_elements_scalar_type(s1,this2) RESULT(t3)
    IMPLICIT NONE
    REAL, INTENT(in)             :: s1
    CLASS(t_element), INTENT(in) :: this2
    TYPE(t_element)              :: t3
    t3%carbon     = s1 / this2%carbon
    t3%nitrogen   = s1 / this2%nitrogen
    t3%phosphorus = s1 / this2%phosphorus
    t3%carbon13   = s1 / this2%carbon13
    t3%carbon14   = s1 / this2%carbon14
    t3%nitrogen15 = s1 / this2%nitrogen15
  END FUNCTION divide_elements_scalar_type  

END MODULE mo_land_memory
