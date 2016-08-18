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
      ka        ,     &   !< Extinction coefficient for PAR on chlorophyll
      sigma_par ,     &   !< single leaf scattering albedo in VIS/PAR range 
      kbl0_par  ,     &   !< extinction coefficient over black leaves 
      kdf0_par  ,     &   !< extinction coefficient for diffuse PAR (Spitters eq. 7) 
      sigma_nir ,     &   !< single leaf scattering albedo in NIR range 
      kbl0_nir  ,     &   !< extinction coefficient over black leaves in NIR range 
      kdf0_nir  ,     &   !< extinction coefficient for diffuse NIR (Spitters eq. 7) 
      rho2sbeta           !< scaling factor of solar angle in reflection calculation 

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants used in the radiation process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_rad_constants

      ka         = 0.005_wp 
      sigma_par  = 0.2_wp    ! Spitters 1986 
      kbl0_par   = 0.5_wp    ! Spitters 1986 
      kdf0_par   = 0.8_wp    ! Spitters 1986
      sigma_nir  = 0.8_wp    ! Spitters 1986
      kbl0_nir   = 0.05_wp   ! following Montheith 1995 in that NIR absorbance is close to none 
      kdf0_nir   = 0.05_wp   ! following Montheith 1995 in that NIR absorbance is close to none 
      rho2sbeta  = 1.6_wp    ! Spitters 1986 

  END SUBROUTINE init_rad_constants
END MODULE mo_rad_constants
