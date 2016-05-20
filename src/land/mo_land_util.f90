!>
!! @brief diverse subroutines and functions for tile & forcing variables
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-23)
!!
MODULE mo_land_util

  USE mo_kind,         ONLY: wp

  IMPLICIT NONE




CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize the derived TYPE t_element element_zero with zero
  !! 
  !! this is done only one time, this element_zero is used to set other derived TYPE's of t_element to zero
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_element_zero

    USE mo_land_memory,  ONLY: element_zero   ! TYPE(t_element) :: element_zero | defined & declared in mo_land_memory

    element_zero%carbon     = 0.0_wp
    element_zero%nitrogen   = 0.0_wp
    element_zero%phosphorus = 0.0_wp
    element_zero%carbon13   = 0.0_wp
    element_zero%carbon14   = 0.0_wp
    element_zero%nitrogen15 = 0.0_wp
  END SUBROUTINE init_element_zero









END MODULE mo_land_util
