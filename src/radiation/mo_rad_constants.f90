!>
!! @brief Contains constants for the radiation process
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-21)
!!
MODULE mo_rad_constants

  USE mo_kind, ONLY: wp


  IMPLICIT NONE

  PUBLIC


  REAL(wp), SAVE :: &      ! In Fortran 2003 and Fortran 2008, a maximum of 255 continuation lines are allowed

  ! Parameters required for albedo and fapar calculation
      ! Extinction coefficient for PAR on chlorophyll
      ka     ,     &   !< Extinction coefficient for PAR on chlorophyll
      ! single leaf scattering albedo
      sigma  ,     &   !< single leaf scattering albedo 
      ! extinction coefficient over black leaves
      kbl0   ,     &   !< extinction coefficient over black leaves 
      ! extinction coefficient for diffuse PAR (Spitters eq. 7)
      kdf0   ,     &   !< extinction coefficient for diffuse PAR (Spitters eq. 7) 
      ! scaling factor of solar angle in reflection calculation (mo_veg_fapar)
      rho2sbeta        !< scaling factor of solar angle in reflection calculation (mo_veg_fapar ?) 

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants used in the radiation process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_rad_constants

      ka         = 0.005_wp 
      sigma      = 0.2_wp 
      kbl0       = 0.5_wp 
      kdf0       = 0.8_wp 
      rho2sbeta  = 1.6_wp    

  END SUBROUTINE init_rad_constants
END MODULE mo_rad_constants
