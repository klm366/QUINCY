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
      kbl0_vis  ,     &   !< extinction coefficient over black leaves 
      kdf0_vis  ,     &   !< extinction coefficient for diffuse PAR (Spitters eq. 7) 
      kbl0_nir  ,     &   !< extinction coefficient over black leaves in NIR range 
      kdf0_nir  ,     &   !< extinction coefficient for diffuse NIR (Spitters eq. 7) 
      rho2sbeta ,     &   !< scaling factor of solar angle in reflection calculation 
      min_cos_zenith_angle_rad  !< minimum coszine zenith angle for radiation calculation

  REAL(wp), SAVE :: & 
      def_alb_vis_soil, & !< default VIS albedo for soil
      def_alb_nir_soil    !< default NIR albedo for soil

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants used in the radiation process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_rad_constants

      kbl0_vis   = 0.5_wp    ! Spitters 1986 
      kdf0_vis   = 0.8_wp    ! Spitters 1986
      kbl0_nir   = 0.5_wp    ! Spitters 1986 
      kdf0_nir   = 0.8_wp    ! Spitters 1986
      rho2sbeta  = 1.6_wp    ! Spitters 1986 
      min_cos_zenith_angle_rad = 0.1736482_wp ! corresponding to 10° solar elevation, below which the 
                                              ! radiation equations are no longer valid
      def_alb_vis_soil = 0.15_wp ! detault from Bonan 2008
      def_alb_nir_soil = 0.30_wp ! default from Bonan 2008

  END SUBROUTINE init_rad_constants
END MODULE mo_rad_constants
