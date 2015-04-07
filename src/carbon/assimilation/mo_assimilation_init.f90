!>
!! @brief Initialization of the the assimilation memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-29)
!!



MODULE mo_assimilation_init
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile
  USE mo_carbon_memory_class, ONLY: t_carbon_memory, t_soil_pools, t_veg_pools, t_photosyn, t_canopy_layer
  USE mo_jsb_base,            ONLY: get_model
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_close_file, jsb_netcdf_open_input, &
                                    jsb_netcdf_read_2d_time, jsb_netcdf_read_2d
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: assimilation_init_canopylayers

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimilation_init'

CONTAINS



   SUBROUTINE assimilation_init_canopylayers
   Sollte der SUBROUTINE bethy_init_io entsprechen...

 END MODULE mo_assimilation_init
