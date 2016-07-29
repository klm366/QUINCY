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

! necessary ?
!  USE mo_base_memory,           ONLY: t_slm_soil, t_slm_tile
!  USE mo_veg_parameters,        ONLY: t_slm_veg_param

! Tasks 
  USE mo_veg_phenology,         ONLY: update_phenology
  USE mo_veg_growth,            ONLY: calc_growth
  USE mo_veg_turnover,          ONLY: calc_turnover
  USE mo_veg_dynamics,          ONLY: update_veg_dynamics
  USE mo_veg_update_veg_pools,  ONLY: update_veg_pools
  USE mo_veg_soil,              ONLY: update_vegsoil_bgc, update_vegsoil_physics, update_soil_pool


  IMPLICIT NONE

!  PRIVATE
!  PUBLIC :: 

CONTAINS

  !! Tasks as called in the land_model_interface !!

    !> Update soil biogeochemistry
    !! 
    !CALL update_vegsoil_bgc(tile, model%lctlib(tile%lcts%lib_id))

    !> update phenological state
    !! 
    !CALL update_phenology(tile, model%lctlib(tile%lcts%lib_id)%phenology_type)

    !> Calculate growth and maintenance
    !! 
    !CALL calc_growth(tile%plant,model%lctlib(tile%lcts%lib_id),tile%soil%theta)

    !> Calculate turnover from tissue aging and senescence
    !! 
    !CALL calc_turnover(tile%plant,model%lctlib(tile%lcts%lib_id))

    !> Calculate vegetation dynamics
    !! 
    !CALL update_veg_dynamics(tile,model%lctlib(tile%lcts%lib_id))

    !> Update soil pools
    !! 
    !CALL update_soil_pool(tile, model%lctlib(tile%lcts%lib_id))

    !> Update soil properties influenced by vegetation ?
    !! 
    !CALL update_vegsoil_physics(tile)

    !> Update Plant pools
    !! 
    !CALL update_veg_pools(tile,model%lctlib(tile%lcts%lib_id))



END MODULE mo_veg_interface
