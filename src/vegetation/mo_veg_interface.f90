!>
!! @brief Contains the interface for plant growth and turnover processes 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!! 2nd version by Jan Engel (20160719)
!!
MODULE mo_veg_interface

#include "dsl_jsbach.inc"

  USE mo_kind
  USE mo_exception,             ONLY: message, finish

  USE mo_veg_growth,            ONLY: update_growth
  USE mo_veg_turnover,          ONLY: update_veg_turnover
  USE mo_veg_dynamics,          ONLY: update_veg_dynamics
  USE mo_veg_update_pools,      ONLY: update_veg_pools
  USE mo_veg_soil,              ONLY: update_vegsoil_bgc


  IMPLICIT NONE
!  PRIVATE
!  PUBLIC :: 

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_veg_interface'

!CONTAINS


END MODULE mo_veg_interface
