!>
!! @brief Contains the interface to the soil processes (the Commission model)
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-19)
!!
MODULE mo_soil_interface

#include "dsl_jsbach.inc"

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message, finish

  USE mo_soil_update_turbulence,  ONLY: update_soil_turbulence
  USE mo_soil_update_physics,     ONLY: update_soil_physics

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_soil_physics, update_soil_turbulence


  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_interface'

!CONTAINS



END MODULE mo_soil_interface
