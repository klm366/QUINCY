!>
!! @brief fundamental constants of the land memory 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_base_constants

  USE mo_kind, ONLY : wp,eps4

  INTEGER, PARAMETER    :: npft       = 8         !< number of PFTs simulated by the model
  INTEGER, PARAMETER    :: ntile      = 1         !< number of tiles simulated by the model; see also itile in the driver
  INTEGER, PARAMETER    :: ncanop     = 20        !< number of canopy layer simulated by the model (default for FAPAR based canopy layers would be 10)
  INTEGER, PARAMETER    :: nsoil      = 15        !< number of soil layer simulated by the model
  ! convenience 
  REAL(wp), PARAMETER   :: zero       = 0.0_wp    !< just zero
  ! time control
  REAL(wp), SAVE        :: dtime                  !< length of one time step in seconds
  REAL(wp), SAVE        :: one_year               !< length of one year in days 
  REAL(wp), SAVE        :: one_day                !< length of one day in seconds
  ! conversion of degree to radians
  REAL(wp), SAVE        :: pi                     !< the number PI
  REAL(wp), SAVE        :: deg2rad                !< conversion of degree to radians
  REAL(wp), PARAMETER   :: degcircle  = 360._wp   !< degrees per circle
  REAL(wp), PARAMETER   :: hourday    = 24.0_wp   !< hours in a day
  ! Canopy layering
  REAL(wp), DIMENSION(ncanop), SAVE :: canopy_layer_depth            !< ...
  REAL(wp)                          :: k_canopy_layer = 0.5_wp       !< extinction coefficient to calculate canopy layer depth

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialise base constants 
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_base_constants

  IMPLICIT NONE

    REAL(wp) :: dAPAR
    INTEGER  :: icanop

    pi      = 4.0_wp * ATAN(1.0_wp)
    deg2rad = pi / 180.0_wp

    ! alternative canopy layer depth calculation, approximately proportional to FAPAR to  
    ! reduce impact of actual canopy layering on PS calculation
    dAPAR = 1._wp/REAL(ncanop,KIND=wp)
    DO icanop=1,ncanop
      canopy_layer_depth(icanop) = -1.0_wp/k_canopy_layer * LOG(1._wp-(dAPAR*icanop)+eps4)
    ENDDO
    ! standard canopy layer calculation
    canopy_layer_depth(:) = 0.5_wp

    one_day   = 86400._wp
    one_year  = 365._wp
  END SUBROUTINE init_base_constants

END MODULE mo_base_constants
