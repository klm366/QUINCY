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

  USE mo_kind


! Tasks 
  USE mo_veg_phenology,         ONLY: update_phenology
  USE mo_veg_growth,            ONLY: update_growth
  USE mo_veg_turnover,          ONLY: update_veg_turnover
  USE mo_veg_dynamics,          ONLY: update_veg_dynamics
  USE mo_veg_update_veg_pools,  ONLY: update_veg_pools
  USE mo_veg_soil,              ONLY: update_vegsoil_bgc, update_vegsoil_physics, update_vegsoil_pool


  IMPLICIT NONE

!  PRIVATE
!  PUBLIC :: 

!CONTAINS


END MODULE mo_veg_interface
