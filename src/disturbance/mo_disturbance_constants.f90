!>
!! @brief Contains constants for the natural disturbances processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-15)
!! 2018-04-25, StW: Added used constants and names - removed unused
!!
MODULE mo_disturbance_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  ! Constants for controlling windbreak algorithms
  CHARACTER(len=*), PARAMETER :: windbreak_names = 'none jsb3'
  INTEGER,          PARAMETER :: WINDBREAK_NONE  = 0
  INTEGER,          PARAMETER :: WINDBREAK_JSB3  = 1

  ! Constants for controlling fire algorithms
  CHARACTER(len=*), PARAMETER :: fire_names      = 'none jsb3 spitfire'
  INTEGER,          PARAMETER :: FIRE_NONE       = 0
  INTEGER,          PARAMETER :: FIRE_JSB3       = 1
  INTEGER,          PARAMETER :: FIRE_SPITFIRE   = 2

  ! Constants regarding temporal smooting of disturbance relevant forcing field
  ! StW: would it make sense to move these to jsb_disturb_fire_jsb3/jsb_disturb_windbreak_jsb3 namelists respectively?
  REAL(wp), PARAMETER         :: persist_rel_hum  = 0.95_wp   ! factor concerning the smoothing of relative air humidity (~14 days)
  REAL(wp), PARAMETER         :: persist_wind_10m = 0.9995_wp ! factor concerning the smoothing of daily maximum wind speed

  ! Physical/chemical constants for Spitfire
  REAL(wp), PARAMETER :: MINER_TOT             =     0.055_wp
  REAL(wp), PARAMETER :: MINER_DAMP            =     0.41739_wp
  REAL(wp), PARAMETER :: fbd10hr_2_1hr         =     1.2_wp
  REAL(wp), PARAMETER :: fbd100hr_2_1hr        =     1.4_wp     ! from fuel bulk densities to 1hr fuel type equivs
  REAL(wp), PARAMETER :: fbd_C3_livegrass      =     4.0_wp     ! fuel density of livegrass [kg/m3]
  REAL(wp), PARAMETER :: fbd_C4_livegrass      =     4.0_wp     ! fuel density of livegrass [kg/m3]
  REAL(wp), PARAMETER :: moistfactor_livegrass =     0.2_wp
  REAL(wp), PARAMETER :: fract_1hr_wood_new    =     0.0_wp     ! frac of    1hr fuel in new carbon in woody litter pools
  REAL(wp), PARAMETER :: fract_10hr_wood_new   =     0.12_wp    ! frac of   10hr fuel in new carbon in woody litter pools
  REAL(wp), PARAMETER :: fract_100hr_wood_new  =     0.21_wp    ! frac of  100hr fuel in new carbon in woody litter pools
  REAL(wp), PARAMETER :: fract_1000hr_wood_new =     0.67_wp    ! frac of 1000hr fuel in new carbon in woody litter pools
  REAL(wp), PARAMETER :: HeatContentFuel       = 18000.0_wp     ! heat content of fuel in [kJ kg-1]
  REAL(wp), PARAMETER :: hmax                  =    24.19_wp    ! maximum vegetation height
  REAL(wp), PARAMETER :: hl                    =     0.19_wp    ! parameter to compute vegetation height
  REAL(wp), PARAMETER :: ft_per_m              =     3.281_wp   ! Conversion from m -> feet
  REAL(wp), PARAMETER :: part_dens             =   513.0_wp
  REAL(wp), PARAMETER :: alpha_1hr             =     0.001_wp
  REAL(wp), PARAMETER :: alpha_10hr            =     0.00005424_wp
  REAL(wp), PARAMETER :: alpha_100hr           =     0.00001485_wp
  REAL(wp), PARAMETER :: alpha_1000hr          =     0.000001_wp
  REAL(wp), PARAMETER :: allom2                =    40._wp
  REAL(wp), PARAMETER :: allom3                =     0.67_wp

#endif
END MODULE mo_disturbance_constants
