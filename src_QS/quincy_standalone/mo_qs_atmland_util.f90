!> Functions needed to calculate variables affecting the atm-land interface
!>
!> QUINCY standalone
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: https://doi.org/10.17871/quincy-model-2019
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_qs_atmland_util

  USE mo_kind

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calc_solar_angle, update_orbital_parameters, calc_fract_par_diffuse

  CHARACTER(len=*), PARAMETER :: modname = 'mo_qs_atmland_util'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> function to calculate solar angle (expressed as cosine of the solar zenith angle)
  !! 
  !! the solar angle is calculated for a given latitude and longitude and time of the day
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_solar_angle(latitude,longitude,time) RESULT(cos_zenith_angle)

    USE mo_jsb_math_constants,    ONLY: deg2rad,degcircle,hourday,pi 
    USE mo_qs_atmland_constants,  ONLY: declination, angle_at_GMT, reference_longitude_GMT

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 Input
    REAL(wp), INTENT(in) :: latitude                    !< latitude in degree
    REAL(wp), INTENT(in) :: longitude                   !< longitude in degree
    REAL(wp), INTENT(in) :: time                        !< days since 1. January
    REAL(wp)             :: cos_zenith_angle            !< cosine of solar zenith angle 
    ! ---------------------------
    ! 0.2 Local
    REAL(wp) :: sin_delta                   ! sine of declination 
    REAL(wp) :: cos_delta                   ! cosine of declination 
    REAL(wp) :: hour                        ! current hour 
    REAL(wp) :: hour_angle                  ! solar hour angle
    REAL(wp) :: angle_offset                ! offset for solar angle
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_solar_angle'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 calculate the solar hour angle in radians
    !! 
    !! 
    hour = ( time - NINT(time-0.5_wp) ) * hourday
    IF(angle_at_GMT) THEN
       angle_offset = degcircle/hourday*(longitude - reference_longitude_GMT)
    ELSE
       angle_offset = 0.0_wp 
    ENDIF
    hour_angle = 2.0_wp*pi*(hour-hourday/2.0_wp+angle_offset)/hourday

    !> 2.0 actual solar angle accounting for current solar declination 
    !! 
    !! 
    cos_delta         = COS(declination)
    sin_delta         = SIN(declination)
    cos_zenith_angle  = SIN(latitude*deg2rad)*sin_delta + & 
                        COS(latitude*deg2rad)*cos_delta*COS(hour_angle)

  END FUNCTION calc_solar_angle


  ! ====================================================================================================== !
  !>
  !> calculate the fraction of diffuse light in the shortwave downward flux at the surface
  !> Note: shortwave diffuse fraction is assumed to be identical to PAR diffuse radiation
  !>
  ELEMENTAL FUNCTION calc_fract_par_diffuse(sw_srf_down, sw_toa_down, cos_zenith_angle) RESULT(fract_par_diffuse)

    USE mo_qs_atmland_constants,    ONLY: k1_fdiff,k2_fdiff,k3_fdiff,k4_fdiff,max_transmission_fdiff, & 
                                          min_fract_par_diffuse, min_cos_zenith_angle

    ! -------------------------------------------------------------------------------------------------- !
    REAL(wp), INTENT(in) :: sw_srf_down       !< surface downward shortwave radiation flux
    REAL(wp), INTENT(in) :: sw_toa_down       !< top-of-atmoshere downward shortwave radiation flux
    REAL(wp), INTENT(in) :: cos_zenith_angle  !< cosine of solar zenith angle
    REAL(wp)             :: fract_par_diffuse 
      !< fraction of surface downward radiation flux that is diffuse (assumed to be identical to PAR diffuse radiation)
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)             :: transmission          ! fraction of radiation reaching surface
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fract_par_diffuse'
    ! -------------------------------------------------------------------------------------------------- !

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    IF(cos_zenith_angle >= min_cos_zenith_angle) THEN  ! only if zenith angle is larger then min_cos_zenith_angle

      transmission = sw_srf_down / ( sw_toa_down * cos_zenith_angle )

      IF(transmission >= max_transmission_fdiff) THEN
         fract_par_diffuse = min_fract_par_diffuse
      ELSE
         fract_par_diffuse = MIN(1.0_wp, k1_fdiff + transmission * & 
                       (k2_fdiff + transmission * ( k3_fdiff + k4_fdiff * transmission )))
      ENDIF

    ELSE ! assume that only diffuse light is available
  
      fract_par_diffuse = 1.0_wp
  
    ENDIF

  END FUNCTION calc_fract_par_diffuse

  !-----------------------------------------------------------------------------------------------------
  !> calculate declination for a given day as a function of orbital parameters
  !!
  !! given this, the potential shortwave radiation at the top of the atmosphere is also calculated \n
  !! the implementation here is taken from ORCHIDEE, as described in 
  !! the PhD thesis of Marie-France Loutre, ASTR-UCL, Belgium, 1993
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_orbital_parameters(time,declination,sw_toa_down) 

    USE mo_jsb_math_constants,    ONLY: deg2rad 
    USE mo_qs_atmland_constants,  ONLY: obliquity, perihel, eccentricity, solar_constant

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp),INTENT(in)  :: time             !< days since January 1st
    REAL(wp),INTENT(out) :: declination      !< solar declination 
    REAL(wp),INTENT(out) :: sw_toa_down(:)   !< shortwave downward flux at the top of the atmosphere (W m-2) with DIMENSION(nc)
    ! ---------------------------
    ! 0.2 Local
    REAL(wp) :: so                           ! sine of obliquity
    REAL(wp) :: sd,cd                        ! sine and cosine of declination
    REAL(wp) :: xl,xllp                      ! perihel in radians
    REAL(wp) :: xee,xse                      ! helper variables for eccentricity calculation
    REAL(wp) :: ranm,rlam,xlam,dlamm,ranv    ! helper variables for anomaly calculation of eath orbit
    REAL(wp) :: dist                         ! distance sun-Earth for potential radiation calculation 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_orbital_parameters'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 calculate orbital parameters
    !! 
    ! sinus of obliquity
    so    = SIN(obliquity*deg2rad)
    !
    ! perihel in radians
    xl    = perihel + 180.0_wp
    xllp  = xl * deg2rad
    !
    ! helpers of eccentricity
    xee   = eccentricity * eccentricity
    xse   = SQRT(1.0_wp - xee)
    !
    ! xlam : true long. sun for mean long. = 0
    xlam  = (eccentricity/2.0_wp + eccentricity * xee/8.0_wp) * (1.0_wp + xse) * SIN(xllp) - & 
             xee/4.0_wp * (0.5_wp + xse) * SIN(2.0_wp*xllp) + & 
             eccentricity * xee/8.0_wp*(1.0_wp/3.0_wp + xse) * SIN(3.0_wp*xllp)
    xlam  = 2.0_wp*xlam/deg2rad
    ! dlamm : mean long. sun for ma-ja
    dlamm = xlam+time-79.0_wp

    ranm  = (dlamm-xl) * deg2rad
    xee   = xee*eccentricity
    ! ranv : true anomaly (radians)
    ranv  = ranm+(2.0_wp*eccentricity-xee/4.0_wp)*SIN(ranm) + & 
            5.0_wp/4.0_wp*eccentricity*eccentricity * SIN(2.0_wp*ranm) + & 
            13.0_wp/12.0_wp*xee*SIN(3.0_wp*ranm)
    ! rlam : true longitude (radians)
    rlam  = ranv + xllp

    !> 2.0 calculate current solar declination 
    !! 
    sd    = so*SIN(rlam)
    cd    = SQRT(1.0_wp-sd*sd)
    ! delta : Solar Declination (radians; angle sun at equator)
    declination = ATAN(sd/cd)

    !> 3.0 calculate potential radiation normal to the atmsphere
    !! 
    ! eccentricity effect: normalised Sun-Earth distance
    dist  = 1.0_wp / (( 1.0_wp - eccentricity * eccentricity ) / & 
            (1.0_wp + eccentricity * cos(ranv)))
    ! potential insolation normal to the atmosphere (W m-2)
    sw_toa_down(:) = dist * dist * solar_constant

  END SUBROUTINE update_orbital_parameters

END MODULE mo_qs_atmland_util
