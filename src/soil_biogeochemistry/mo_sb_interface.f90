!>
!! @brief Contains the interface to the soil_biogeochemistry processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20170104)
!!
MODULE mo_sb_interface

#include "dsl_jsbach.inc"

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish

  USE mo_sb_update_pools,     ONLY: update_sb_pools

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_sb_pools 

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_interface'

!CONTAINS



END MODULE mo_sb_interface
