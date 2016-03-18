!!
!! @brief diverse subroutines and functions for the plant-model
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-17)
!!
MODULE mo_veg_util

  USE mo_kind, ONLY: wp
  USE mo_base_memory, ONLY: TileClass
  USE mo_base_util_functions, ONLY: moving_avrg
  USE mo_base_constants, ONLY: dt, one_day
  USE mo_land_parameters, ONLY: AccumulationIntervalTileClass
  USE mo_veg_parameters, ONLY: AccumulationIntervalPlantClass
  IMPLICIT NONE

  TYPE(AccumulationIntervalPlantClass) :: plant_acc_interval ! intervals over which variables are averaged mo_veg_parameters
  TYPE(AccumulationIntervalTileClass)  :: tile_acc_interval  ! intervals over which variables are averaged mo_land_parameters


CONTAINS

  SUBROUTINE accumulate_moving_averages_plant(tile, ntile, dt)
    ! calculate/ accumulate the moving average of certain vegetation variables across a certain interval of timesteps
    ! for initial average values of vegetation variables see mo_veg_memory; init subroutine called in plant_model_driver
    IMPLICIT NONE
    ! 0.1 Input
    INTEGER, INTENT(in)  :: ntile
    REAL(wp), INTENT(in) :: dt
    ! 0.2. In- & Output
    TYPE(TileClass), DIMENSION(ntile), INTENT(inout) :: tile 

    ! moving_avrg( current average, new value, weight of new value, use function TRUE FALSE [OPTIONAL argument] )
    tile(1)%plant%accumulated%par_sunlit    = moving_avrg( &
                                              tile(1)%plant%accumulated%par_sunlit, tile(1)%plant%canopy_layer(1)%par_sunlit, &
                                              1.0_wp / (plant_acc_interval%par_sunlit * one_day/ dt), .TRUE.) 
    tile(1)%plant%accumulated%par_shaded    = moving_avrg( &
                                              tile(1)%plant%accumulated%par_shaded, tile(1)%plant%canopy_layer(1)%par_shaded, &
                                              1.0_wp / (plant_acc_interval%par_shaded * one_day/ dt), .TRUE.) 
    tile(1)%plant%accumulated%fleaf_sunlit  = moving_avrg( &
                                              tile(1)%plant%accumulated%fleaf_sunlit, tile(1)%plant%canopy_layer(1)%fleaf_sunlit, &
                                              1.0_wp / (plant_acc_interval%fleaf_sunlit * one_day/ dt) )!, .TRUE.) 
    tile(1)%plant%accumulated%gpp_week      = moving_avrg( & 
                                              tile(1)%plant%accumulated%gpp_week, tile(1)%plant%flux%gpp, & 
                                              1.0_wp / (plant_acc_interval%weekly * one_day/dt))
    tile(1)%plant%accumulated%maint_respiration_week = moving_avrg(tile(1)%plant%accumulated%maint_respiration_week,& 
                                                       tile(1)%plant%flux%maint_respiration, & 
                                                       1._wp/(plant_acc_interval%weekly*one_day/dt))

  END SUBROUTINE accumulate_moving_averages_plant

END MODULE mo_veg_util
