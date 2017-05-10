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
  USE mo_sb_rate_modifier,    ONLY: update_sb_rate_modifier
  USE mo_sb_simple_model,     ONLY: update_sb_simple_model
  USE mo_sb_jsm_main,         ONLY: update_sb_jsm


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_sb_pools, update_sb_rate_modifier
  PUBLIC :: update_sb_simple_model, update_sb_jsm

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_interface'

!CONTAINS



END MODULE mo_sb_interface
