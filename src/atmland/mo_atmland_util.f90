!>
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

  PUBLIC :: calc_solar_angle, update_orbital_parameters, calc_sw_fdiffuse, calc_spec_humidity_sat

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

  FUNCTION calc_sw_fdiffuse(sw_srf_down,sw_toa_down,cos_zenith_angle) RESULT(sw_fdiffuse)
    ! function to calculate the fraction of diffuse light in the shortwave downward flux
    ! at the surface
    USE mo_atmland_constants, ONLY : k1_fdiff,k2_fdiff,k3_fdiff,k4_fdiff,max_transmission_fdiff, & 
                                     min_sw_fdiffuse, min_cos_zenith_angle
    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp) :: sw_srf_down           ! surface downward shortwave radiation flux
    REAL(wp) :: sw_toa_down           ! top-of-atmoshere downward shortwave radiation flux
    REAL(wp) :: cos_zenith_angle      ! cosine of solar zenith angle
    ! 0.3 Output
    REAL(wp) :: sw_fdiffuse           ! fraction of surface downward radiation flux that is diffuse
    ! 0.4 Local
    REAL(wp) :: transmission          ! fraction of radiation reaching surface

    IF(cos_zenith_angle >= min_cos_zenith_angle) THEN  ! only if zenith angle is larger then min_cos_zenith_angle

      transmission = sw_srf_down / ( sw_toa_down * cos_zenith_angle )     

      IF(transmission >= max_transmission_fdiff) THEN
         sw_fdiffuse = min_sw_fdiffuse
      ELSE
         sw_fdiffuse = MIN(1.0_wp, k1_fdiff + transmission * & 
            (k2_fdiff + transmission * ( k3_fdiff + k4_fdiff * transmission )))
      ENDIF

    ELSE ! assume that only diffuse light is available
  
      sw_fdiffuse = 1.0_wp
  
    ENDIF

  END FUNCTION calc_sw_fdiffuse

  PURE ELEMENTAL FUNCTION calc_spec_humidity_sat(temp,press_srf) RESULT(q_sat)
        ! constants for the calculation of atmospheric water content at saturation
        ! from Campell and Norman's Introduction to Environmental Biophysics
     USE mo_kind, ONLY: wp
     USE mo_atmland_constants, ONLY: eps_vpd, Tzero
     IMPLICIT NONE 
     ! 0.1 Input
     REAL(wp), INTENT(in) :: temp            ! temperature (K)
     REAL(wp), INTENT(in) :: press_srf       ! atmospheric pressure (Pa) 
     ! 0.3 Output 
     REAL(wp) :: q_sat           ! saturating specific humidity (g/g) 
     ! 0.4 Local
     REAL(wp) :: esat            ! saturating vapour pressure (Pa)
     REAL(wp) :: a               ! Pa
     REAL(wp) :: b               ! einheitenlos
     REAL(wp) :: c               ! °C

     a = 611.0_wp    
     b = 17.502_wp  
     c = 240.97_wp    
     ! calculating saturating vapour pressure [Pa] (Campbell & Norman)
     esat = a * exp(b * (temp - Tzero)/(temp - Tzero + c))
     ! conversion to saturating specific humidity [kg kg-1]
     q_sat = eps_vpd * esat/(press_srf - (1.0_wp - eps_vpd) * esat)

  END FUNCTION calc_spec_humidity_sat

  SUBROUTINE update_orbital_parameters(time,declination,sw_toa_down) 
    ! subroutine to calculate declination for a given day as a function of orbital parameters
    ! given this, the potential shortwave radiation at the top of the atmosphere is also calculated
    ! the implementation here is taken from ORCHIDEE, as described in 
    ! the PhD thesis of Marie-France Loutre, ASTR-UCL, Belgium, 1993

    USE mo_base_constants, ONLY : deg2rad 
    USE mo_atmland_constants, ONLY : obliquity,perihel,eccentricity,solar_constant

    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp),INTENT(in) :: time              ! days since January 1st
    ! 0.3 Output
    REAL(wp),INTENT(out) :: declination      ! solar declination 
    REAL(wp),INTENT(out) :: sw_toa_down      ! shortwave downward flux at the top of the atmosphere (W m-2) 
    ! 0.4 Local
    REAL(wp) :: so                           ! sine of obliquity
    REAL(wp) :: sd,cd                        ! sine and cosine of declination
    REAL(wp) :: xl,xllp                      ! perihel in radians
    REAL(wp) :: xee,xse                      ! helper variables for eccentricity calculation
    REAL(wp) :: ranm,rlam,xlam,dlamm,ranv    ! helper variables for anomaly calculation of eath orbit
    REAL(wp) :: dist                         ! distance sun-Earth for potential radiation calculation 
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
    !
    ! 3.0 calculate potential radiation normal to the atmsphere
    ! 
    ! eccentricity effect: normalised Sun-Earth distance
    dist = 1.0_wp / (( 1.0_wp - eccentricity * eccentricity ) / & 
           ( 1.0_wp + eccentricity * cos(ranv)))
    ! potential insolation normal to the atmosphere (W m-2)
    sw_toa_down = dist * dist * solar_constant

  END SUBROUTINE update_orbital_parameters

END MODULE mo_atmland_util
