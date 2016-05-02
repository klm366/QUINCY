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

  USE mo_kind, ONLY: wp

  IMPLICIT NONE

  PUBLIC :: update_tile_moving_averages 

CONTAINS

  !> calculate/ accumulate the moving average of certain vegetation variables across a certain interval of timesteps
  !! for initial average values of vegetation variables see mo_veg_memory; init subroutine called in plant_model_driver
  SUBROUTINE update_tile_moving_averages(tile)

    USE mo_base_memory,     ONLY: t_slm_tile
    USE mo_base_util,       ONLY: moving_avrg
    USE mo_base_constants,  ONLY: dt, one_day, ntile
    USE mo_atmland_memory , ONLY: atm2land
    USE mo_land_parameters, ONLY: AccumulationIntervalTileClass
    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile), DIMENSION(ntile), INTENT(inout) :: tile 
    ! 0.2 Local
    TYPE(AccumulationIntervalTileClass) :: tile_acc_interval ! intervals over which variables are averaged mo_land_parameters

    ! moving_avrg( current average, new value, weight of new value, use function TRUE FALSE [OPTIONAL argument] )
    tile(:)%accumulated%t_air_week          = moving_avrg( &
                                              tile(:)%accumulated%t_air_week, atm2land%t_air, &
                                              1.0_wp / (tile_acc_interval%weekly * one_day/ dt))
    tile(:)%accumulated%t_air_month         = moving_avrg( &
                                              tile(:)%accumulated%t_air_month, atm2land%t_air, &
                                              1.0_wp / (tile_acc_interval%monthly * one_day/ dt))
    tile(:)%accumulated%t_air_longterm      = moving_avrg( &
                                              tile(:)%accumulated%t_air_longterm, atm2land%t_air, &
                                              1.0_wp / (tile_acc_interval%longterm * one_day/ dt))
  END SUBROUTINE update_tile_moving_averages

END MODULE mo_land_util
