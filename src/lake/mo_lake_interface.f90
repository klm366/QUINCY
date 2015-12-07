!>
!! @brief Contains the interfaces to the lake processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_lake_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_class,       ONLY: get_model
  USE mo_jsb_tile_class,  ONLY: t_jsb_tile

  dsl4jsb_Use_processes SRF, SOIL, ATM2LAND
  dsl4jsb_Use_config(SRF)
  dsl4jsb_Use_config(SOIL)

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SRF)

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lake_interface'

!CONTAINS

#endif
END MODULE mo_lake_interface
