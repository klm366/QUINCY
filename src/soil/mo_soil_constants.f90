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
    soil_therm_cond, &          !< soil thermal conductivity (W m-2 K-1)
    theta_sat, &                !< soil volumetric water content at saturation (fraction)
    theta_wilt, &               !< soil volumetric water content at wilting (fraction)
    k_soil_phi1,  &             !< coefficient for calculating soil water potential from soil moisture
    k_soil_phi2                 !< coefficient for calculating soil water potential from soil moisture 

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_soil_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants from soil 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_soil_constants

    IMPLICIT NONE
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_soil_constants'

    soil_heat_cap   = 1400.0_wp   ! textbook for wet sandy soil
    rho_bulk_soil   = 1695.86_wp  ! Bernhard Ahrens for Hainich
    soil_therm_cond = 0.3_wp      ! Text book for wet sandy soil
    theta_sat   = 0.49_wp         ! based on Rawls et al 1982 for loam
    theta_wilt  = 0.1_wp          ! based on Rawls et al 1982 for loam
    k_soil_phi1 = -0.0011_wp      ! based on Rawls et al 1982 for loam
    k_soil_phi2 = -4.55_wp        ! based on Rawls et al 1982 for loam

  END SUBROUTINE init_soil_constants

END MODULE mo_soil_constants
