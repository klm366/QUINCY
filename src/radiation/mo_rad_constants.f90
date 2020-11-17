!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Contains constants for the radiation process
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-21)
!!
MODULE mo_rad_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  ! Parameters required for albedo and fapar calculation
  REAL(wp), SAVE ::   & 
      kbl0_vis  ,     &     !< extinction coefficient over black leaves 
      kdf0_vis  ,     &     !< extinction coefficient for diffuse PAR (Spitters eq. 7) 
      kbl0_nir  ,     &     !< extinction coefficient over black leaves in NIR range 
      kdf0_nir  ,     &     !< extinction coefficient for diffuse NIR (Spitters eq. 7) 
      rho2sbeta ,     &     !< scaling factor of solar angle in reflection calculation 
      k_csf     ,     &     !< scaling factor in the crown-shape effect on clumping (Campbell & Norman, 1998, eq 15.35)
      min_cos_zenith_angle_rad  !< minimum coszine zenith angle for radiation calculation

  REAL(wp), SAVE ::     & 
      def_alb_vis_soil, &   !< default VIS albedo for soil
      def_alb_nir_soil      !< default NIR albedo for soil

  ! conversion factors
  REAL(wp), SAVE ::   &
      rad2ppfd              !< µmol /J

  ! red to far red ratio calculations
  REAL(wp), SAVE ::  &
      rfr_ratio_toc, &      !< top of the canopy red:farred ratio  (Kull & Kruij 1998) 
      k_r2fr_chl            !< extinction of red:farred ratio on leaf chlorophyll (1/mmol N); 
                            !! (Kull & Krujit 1998 for the far red value, tuned to get a R:FR ratio
                            !! of 0.15 at full canopy closure (a common value) 

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_rad_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants used in the radiation process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_rad_constants

    IMPLICIT NONE
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_rad_constants'

    kbl0_vis                  = 0.5_wp        !< Spitters 1986 
    kdf0_vis                  = 0.8_wp        !< Spitters 1986 
    kbl0_nir                  = 0.5_wp        !< Spitters 1986 
    kdf0_nir                  = 0.8_wp        !< Spitters 1986 
    rho2sbeta                 = 1.6_wp        !< Spitters 1986 
    k_csf                     = 2.2_wp        !< Campbell & Norman 1998, eq 15.35
    min_cos_zenith_angle_rad  = 0.1736482_wp  !< corresponding to 10° solar elevation, below which the 
                                              !! radiation equations are no longer valid

    def_alb_vis_soil          = 0.15_wp       !< default from Bonan 2008
    def_alb_nir_soil          = 0.30_wp       !< default from Bonan 2008

    rad2ppfd                  = 4.6_wp        !< Monteith & Unsworth 1995 

    rfr_ratio_toc             = 1.2_wp        !< Kull & Kruij 1998
    k_r2fr_chl                = 0.0012_wp - 0.012_wp !<  Kull & Kruijt 1998

  END SUBROUTINE init_rad_constants

#endif
END MODULE mo_rad_constants
