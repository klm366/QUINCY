!>
!! @brief Contains constants for the soil process
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-21)
!!
MODULE mo_soil_constants

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  REAL(wp), SAVE :: &
    soil_heat_cap,  &           !< soil heat capacity (temporarily for SLM offline) J/kg/K 
    rho_bulk_soil,  &           !< soil bulk denisty (kg / m3) 
    soil_therm_cond, &          !< soil thermal conductivity (W m-2 K-1)
    theta_sat, &                !< soil volumetric water content at saturation (fraction)
    theta_wilt, &               !< soil volumetric water content at wilting (fraction)
    k_soil_phi1, &              !< coefficient for calculating soil water potential from soil moisture
    k_soil_phi2, &              !< coefficient for calculating soil water potential from soil moisture 
! constants in simple multi-layer soil moisture scheme to describe the 
    InterceptionEfficiency, &   !< Efficiency of interception of precipitation as rain
    w_skin_max, &               !< maximum water storage per unit LAI (m)
    drain_frac, &               !< preferential flow fraction of infiltrating water per meter soil 
    fdrain_runoff, &            !< preferential flow fraction of runoff
    max_infiltration_rate, &    !< maximum infiltration rate for bucket model 
    kdiff_theta, &              !< soil diffusivity for water [m-1]
    min_theta_wpot              !< minimum fractional water capacity for the calculation of soil water potential

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_soil_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants from soil 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_soil_constants

    IMPLICIT NONE
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_soil_constants'

    soil_heat_cap   = 1400.0_wp    ! textbook for wet sandy soil
    rho_bulk_soil   = 1695.86_wp   ! Bernhard Ahrens for Hainich
    soil_therm_cond = 0.3_wp       ! Text book for wet sandy soil
    theta_sat       = 0.49_wp      ! based on Rawls et al 1982 for loam
    theta_wilt      = 0.1_wp       ! based on Rawls et al 1982 for loam
    k_soil_phi1     = -0.0011_wp   ! based on Rawls et al 1982 for loam
    k_soil_phi2     = -4.55_wp     ! based on Rawls et al 1982 for loam

    InterceptionEfficiency = 0.25_wp ! JSBACH4 
    w_skin_max = 0.0002_wp           ! JSBACH4 
    drain_frac      = 0.01_wp      ! tuned
    fdrain_runoff   = 0.95_wp      ! set from ORCHIDEE
    max_infiltration_rate = 0.0025 ! mm / s (from 10 mm/hour according to FAO loam)
    kdiff_theta     = 0.00000005_wp ! set 
    min_theta_wpot  = 0.1_wp       ! for numerical stability, implies minimim water potential of -39 MPa

  END SUBROUTINE init_soil_constants

END MODULE mo_soil_constants
