!!
!! @brief fundamental constants of the land memory 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_base_constants

  USE mo_kind, ONLY : wp

  INTEGER, PARAMETER :: npft = 2           ! number of PFTs simulated by the model
  INTEGER, PARAMETER :: ntile = 1          ! number of tiles simulated by the model 
  INTEGER, PARAMETER :: ncanop = 20        ! number of canopy layer simulated by the model
  INTEGER, PARAMETER :: nsoil = 15         ! number of soil layer simulated by the model
  ! convenience 
  REAL(wp), PARAMETER :: zero = 0.0_wp 
  ! time control
  REAL(wp), SAVE :: dt
  REAL(wp), SAVE :: one_year
  ! Canopy layering
  REAL(wp), DIMENSION(ncanop), SAVE :: canopy_layer_depth 

CONTAINS

  SUBROUTINE init_canopy_layer_depth

  canopy_layer_depth(:) = 0.5_wp

  END SUBROUTINE init_canopy_layer_depth

END MODULE mo_base_constants
