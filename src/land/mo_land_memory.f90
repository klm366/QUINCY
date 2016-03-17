!!
!! @brief Contains the routines for land memory operations 
!!
!! @author
!!  Jan Engel, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-10)
!!
MODULE mo_land_memory

  USE mo_kind, ONLY : wp

  IMPLICIT NONE


CONTAINS

  SUBROUTINE init_tile_variable_accumulation(tile)
    ! set intial values for parameters that accumulate values (average) across a certain period
    ! called in plant_model_driver section 1.0
    ! the interval is defined in mo_land_parameters
    USE mo_base_memory
    IMPLICIT NONE
    ! 0.2 Output
    TYPE(TileClass), DIMENSION(ntile), INTENT(out) :: tile 
    ! 
    ! 1.0 initialize default values
    ! 
    tile(:)%accumulated%t_air             = 1._wp
    tile(:)%accumulated%q_air             = 1._wp
    tile(:)%accumulated%press_surf        = 1._wp
  END SUBROUTINE init_tile_variable_accumulation

END MODULE mo_land_memory
