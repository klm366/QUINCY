!!
!! @brief functions needed to calculate variables affecting the atm-land interface 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-19)
!!
MODULE mo_atmland_util

  USE mo_kind

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calc_solar_angle,calc_declination

CONTAINS

  FUNCTION calc_solar_angle(latitude,longitude,time) RESULT(cos_zenith_angle)
    ! function to calculate solar angle (expressed as cosine of the solar zenith angle) 
    ! for a given latitude and longitude and time of the day
    USE mo_base_constants, ONLY: deg2rad,degcircle,hourday,pi 
    USE mo_atmland_constants, ONLY: declination,angle_at_GMT,reference_longitude_GMT

    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp) :: latitude                    ! latitude in degree
    REAL(wp) :: longitude                   ! longitude in degree
    REAL(wp) :: time                        ! days since 1. January
    ! 0.3 Output
    REAL(wp) :: cos_zenith_angle            ! cosine of solar zenith angle 
    ! 0.4 Local variables
    REAL(wp) :: sin_delta                   ! sine of declination 
    REAL(wp) :: cos_delta                   ! cosine of declination 
    REAL(wp) :: hour                        ! current hour 
    REAL(wp) :: hour_angle                  ! solar hour angle
    REAL(wp) :: angle_offset                ! offset for solar angle
    !
    ! 1.0 calculate the solar hour angle in radians
    !  
    hour = ( time - NINT(time-0.5_wp) ) * hourday
    IF(angle_at_GMT) THEN
       angle_offset = degcircle/hourday*(longitude - reference_longitude_GMT)
    ELSE
       angle_offset = 0.0_wp 
    ENDIF
    hour_angle = 2.0_wp*pi*(hour-hourday/2.0_wp+angle_offset)/hourday
    ! 
    ! 2.0 actual solar angle accounting for current solar declination 
    ! 
    cos_delta = COS(declination)
    sin_delta = SIN(declination)
    cos_zenith_angle=SIN(latitude*deg2rad)*sin_delta & 
                     + COS(latitude*deg2rad)*cos_delta*COS(hour_angle)

  END FUNCTION calc_solar_angle

  FUNCTION calc_declination(time) RESULT(declination)
    ! function to calculate declination for a given day as a function of orbital parameters
    ! the implementation here is taken from ORCHIDEE, as described in 
    ! the PhD thesis of Marie-France Loutre, ASTR-UCL, Belgium, 1993

    USE mo_base_constants, ONLY : deg2rad 
    USE mo_atmland_constants, ONLY : obliquity,perihel,eccentricity

    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp) :: time                         ! days since January 1st
    ! 0.3 Output
    REAL(wp) :: declination                  ! solar declination 
    ! 0.4 Local
    REAL(wp) :: so                           ! sine of obliquity
    REAL(wp) :: sd,cd                        ! sine and cosine of declination
    REAL(wp) :: xl,xllp                      ! perihel in radians
    REAL(wp) :: xee,xse                      ! helper variables for eccentricity calculation
    REAL(wp) :: ranm,rlam,xlam,dlamm,ranv    ! helper variables for anomaly calculation of eath orbit
    ! 
    ! 1.0 calculate orbital parameters
    ! 
    ! sinus of obliquity
    so = SIN(obliquity*deg2rad)
    !
    ! perihel in radians
    xl = perihel + 180.0_wp
    xllp = xl * deg2rad
    !
    ! helpers of eccentricity
    xee = eccentricity * eccentricity
    xse = SQRT(1.0_wp - xee)
    !
    ! xlam : true long. sun for mean long. = 0
    xlam = (eccentricity/2.0_wp + eccentricity * xee/8.0_wp) * (1.0_wp + xse) * SIN(xllp) & 
            - xee/4.0_wp * (0.5_wp + xse) * SIN(2.0_wp*xllp) & 
            + eccentricity * xee/8.0_wp*(1.0_wp/3.0_wp + xse) * SIN(3.0_wp*xllp)
    xlam  = 2.0_wp*xlam/deg2rad
    ! dlamm : mean long. sun for ma-ja
    dlamm = xlam+time-79.0_wp

    ranm  = (dlamm-xl) * deg2rad
    xee = xee*eccentricity
    ! ranv : true anomaly (radians)
    ranv = ranm+(2.0_wp*eccentricity-xee/4.0_wp)*SIN(ranm) & 
           + 5.0_wp/4.0_wp*eccentricity*eccentricity * SIN(2.0_wp*ranm) & 
           + 13.0_wp/12.0_wp*xee*SIN(3.0_wp*ranm)
    ! rlam : true longitude (radians)
    rlam =  ranv + xllp
    ! 
    ! 2.0 calculate current solar declination 
    ! 
    sd = so*SIN(rlam)
    cd = SQRT(1.0_wp-sd*sd)
    ! delta : Solar Declination (radians; angle sun at equator)
    declination = ATAN(sd/cd)

  END FUNCTION calc_declination

END MODULE mo_atmland_util
