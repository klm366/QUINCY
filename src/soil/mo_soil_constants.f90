!>
!! @brief Contains constants for the soil process
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-21)
!!
MODULE mo_soil_constants

  USE mo_kind, ONLY: wp

  IMPLICIT NONE

  PUBLIC


!  REAL(wp), PARAMETER :: &      ! In Fortran 2003 and Fortran 2008, a maximum of 255 continuation lines are allowed

  REAL(wp), SAVE :: &
    soil_heat_cap,  &           !< soil heat capacity (temporarily for SLM offline) J/kg/K 
    rho_bulk_soil,  &           !< soil bulk denisty (kg / m3) 
    soil_therm_cond             !< soil thermal conductivity (W m-2 K-1)


CONTAINS
  !-----------------------------------------------------------------------------------------------------
  !> initialize constants from soil 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_soil_constants

    soil_heat_cap   = 1400.0_wp   ! textbook for wet sandy soil
    rho_bulk_soil   = 1695.86_wp  ! Bernhard Ahrens for Hainich
    soil_therm_cond = 0.3_wp      ! Text book for wet sandy soil

  END SUBROUTINE init_soil_constants

END MODULE mo_soil_constants
