!>
!! @brief Contains the routines for land memory operations 
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



END MODULE mo_land_memory
