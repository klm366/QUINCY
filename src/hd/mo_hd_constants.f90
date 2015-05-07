!>
!! @brief Contains constants for the hd processes
!!
!! @author
!!  Veronika Gayler, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Veronika Gayler (2013-04-20)
!!
MODULE mo_hd_constants
#ifndef __NO_JSBACH__

  USE mo_physical_constants, ONLY: &
    rhoh2o              !< Density of liquid water [kg/m3]
  USE mo_math_constants,     ONLY: &
    pi,               & !< PI
    deg2rad             !< conversion factor from degree to radians: pi/180
  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  !
  ! HD grid parameters
  !
  INTEGER,  PARAMETER :: &
    nl = 720,                          & !< number of longitudes
    nb = 360                             !< number of latitudes

  REAL(wp), PARAMETER :: &
    fullcirc    = 360.0_wp,            & !< full circle: 360 degree
    hd_scal_lon = fullcirc/nl,         & !<
    hd_scal_lat = 0.5_wp*fullcirc/nb,  & !<
    florg       = -180.0_wp,           & !< north-west corner longitude of gridbox(1,1)
    fborg       =   90.0_wp,           & !< north-west corner latitude of gridbox(1,1)
    fscal       =    0.5_wp              !< grid resolution [degree]

  !
  ! riverflow parameters
  !
  INTEGER, PARAMETER  :: &
    nmemrf = 5                           !< number of riverflow reservoir cascades


#endif
END MODULE mo_hd_constants
