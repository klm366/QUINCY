!>
!! @brief 
!!  interface for land model \n  
!!  run SLM for one time step
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!!  First version by Sönke Zaehle (2015-12-04) \n
!!  Reworked by Jan Engel (2016 July)
!!
MODULE mo_land_model_interface

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: land_model_interface


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> land model interface
  !! 
  !! @todo land_model: \n
  !!  current work: making it equivalent to plant model functionalilties, but 
  !!  putting all JSBACH4 stuff into the land model driver
  !! 
  !! ==> veg_param is only temp. used here but will be replaced by lctlib soon !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE land_model_interface(tile, veg_param, &
                                  t_air, press_srf, q_air, &
                                  time_day, time_year, timestep_prop_day &
                                 )

    USE mo_kind
    USE mo_base_memory,             ONLY: t_slm_tile
    USE mo_veg_util,                ONLY: set_zero_fluxes_of_plant_pools
    USE mo_veg_parameters,          ONLY: t_slm_veg_param
    USE mo_veg_interface,           ONLY: update_vegetation
    USE mo_veg_radiation,           ONLY: update_radiation
    USE mo_veg_soil,                ONLY: update_vegsoil_bgc, update_vegsoil_physics, update_soil_pool
    USE mo_veg_output,              ONLY: write_output_pool, write_output_flux

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)      :: tile              !< data structure for dynamic vegetation and soil variables
    TYPE(t_slm_veg_param),INTENT(in)    :: veg_param         !< vegetation parameters - replaced by the use of lctlib soon !!
    REAL(wp)                            :: time_day          !< displays the 'real time' in days; time_day == 1.0 after one day has passed
    INTEGER                             :: time_year         !< displays the 'real time' in years; the 1st year of a simulation == 1 (not zero)
    REAL(wp)                            :: timestep_prop_day !< proportion of one day that is simulated with one timestep (dt/one_day)
    REAL(wp)                            :: t_air, &          !< air temperature - atm2land
                                           press_srf, &      !< air pressure at surface - atm2land
                                           q_air             !< air ... - atm2land
    ! 0.2 Local



    !! no tile loop !!

      !! 1 set vegetation fluxes to zero
      !! 
      !! ==> integrate into vegetation process
      CALL set_zero_fluxes_of_plant_pools(tile)  !< set to zero all fluxes of the plant pools; mo_veg_util

      !! 2 canopy calculation fluxes
      !! 
      !! ==> process radiation
      CALL update_radiation(tile)  !< mo_veg_radiation

      !! 3 update soil (biogeochemical part - here only N/P uptake)
      !! 
      !! ==> 
      CALL update_vegsoil_bgc(tile, veg_param)  !< mo_veg_soil

      !! 4 calculate vegetation growth and turnover
      !! 
      !! 
      CALL update_vegetation(tile, veg_param)   !< mo_veg_interface

      !! 5 update soil litter and mineral N pools if using basic soil model
      !! 
      !! 
      !IF (veg_config%flag_base_soil_model) CALL update_soil_pool(tile, veg_param)  !< mo_veg_soil
      CALL update_soil_pool(tile, veg_param)  !< mo_veg_soil

      !! 6 update soil water bucket 
      !! 
      !! 
      CALL update_vegsoil_physics(tile)  !< mo_veg_soil

      !! 7 Accumulation
      !! 
      !! daily & moving average
      CALL tile%plant%daily_average%calc_avrg(tile%plant, time_day, timestep_prop_day) !< daily avrg of the prev day
      CALL tile%plant%moving_average%update_veg(tile%plant)          !< moving averages - vegetation
      !CALL tile%plant%moving_average%update_soil(tile%plant)         !< moving averages - soil -> plant N & P uptake
      CALL tile%plant%moving_average%update_tile(t_air, press_srf, q_air) !< moving averages - atmland

      !! 8 output pool: write output in certain intervals, see mo_veg_output
      !! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
      !! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
      !! 
      CALL write_output_pool(tile, veg_param, time_day, time_year, 1) ! 1 == itile

      !! 9 output flux: write output in certain intervals, see mo_veg_output
      !! 
      !! 
      CALL write_output_flux(tile, time_day, time_year, 1) ! 1 == itile 


  END SUBROUTINE land_model_interface
END MODULE mo_land_model_interface
