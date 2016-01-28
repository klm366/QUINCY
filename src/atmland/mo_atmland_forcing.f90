!!
!! @brief dummy forcing routine for generic forcing 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2016-01-27)
!!
MODULE mo_atmland_forcing

  USE mo_kind

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: generate_inst_temperature 

CONTAINS

  FUNCTION generate_inst_temperature(latitude,time,declination,annual_t_air_min, &
               annual_t_air_max,diurnal_temp_range) RESULT(temperature)
    ! function to generate instantaneous temperature from annual characteristics 
    ! note: meant to be used for generic testing only. Only works for present-day
    USE mo_atmland_constants, ONLY: solstice
    USE mo_atmland_util, ONLY: update_orbital_parameters
    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp),INTENT(in) :: latitude            ! latitude (degree)
    REAL(wp),INTENT(in) :: time                ! time since 1st Januar (days)
    REAL(wp),INTENT(in) :: declination         ! solar declination (radians?)  
    REAL(wp),INTENT(in) :: annual_t_air_min    ! annual minimum daily temperature (K)
    REAL(wp),INTENT(in) :: annual_t_air_max    ! annual maximum daily temperature (K)
    REAL(wp),INTENT(in) :: diurnal_temp_range  ! diurnal temperature range (K)
    ! 0.3 Output
    REAL(wp) :: temperature                    ! current air temperature (K)
    ! 0.4 Local
    REAL(wp) :: declination_max                ! declination at solstice 
    REAL(wp) :: dtmp                           ! today's temperature (K)
    REAL(wp) :: hlp1                           ! dummy variable

    CALL update_orbital_parameters(solstice,declination_max,hlp1)

    dtmp = annual_t_air_min + (annual_t_air_max - annual_t_air_min) * & 
                  ( 1.0_wp - ( declination_max - declination ) / ( 2.0_wp * declination_max ))

    temperature = diurnal_from_daily_temperature(latitude,time,declination,dtmp,diurnal_temp_range)

  END FUNCTION 

  FUNCTION diurnal_from_daily_temperature(latitude,time,declination,dtmp,diurnal_temp_range) & 
           RESULT(temperature)
    ! Function to calculate the temperature at a given time of the day given 
    ! daily minimum and maximum temperatures, solar declination and time of the day
    USE mo_base_constants, ONLY: deg2rad,pi
    USE mo_atmland_constants, ONLY: min_daylength,tmp_peak_time,dusk_time

    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp),INTENT(in) :: latitude             ! latitude (degree)
    REAL(wp),INTENT(in) :: time                 ! time since 1st Januar (days)
    REAL(wp),INTENT(in) :: declination          ! solar declination (radians?)  
    REAL(wp),INTENT(in) :: dtmp                 ! daily mean temperature (K)
    REAL(wp),INTENT(in) :: diurnal_temp_range   ! daily temperature range (K)
    ! 0.3 Output
    REAL(wp) :: temperature                     ! current air temperature 
    ! 0.4 Local
    REAL(wp) :: frac_time                       ! time of day as fraction 
    REAL(wp) :: spds,cpds,arg,sd                ! factors to calculate daylight hourse
    REAL(wp) :: frac_daylight                   ! fraction of the day that is light
    REAL(wp) :: sunrise,sunset                  ! fractional time of the day at which sunrise and sunset occurs
    REAL(wp) :: time_since_sunset               ! helper variable (temperature at sunset)
    REAL(wp) :: dtmp_min                        ! daily minimum temperature (K)
    REAL(wp) :: dtmp_max                        ! daily maximum temperature (K)

    ! work out today's daylength
    spds = SIN(latitude*deg2rad)*SIN(declination)
    cpds = COS(latitude*deg2rad)*COS(declination)
    arg = spds/(cpds+eps8)
    IF (arg > 1._wp) THEN
      frac_daylight = 0.0_wp
    ELSE IF(arg < -1.0_wp ) THEN 
      frac_daylight = 1.0_wp
    ELSE
      frac_daylight = ACOS(arg) / pi
    ENDIF
    IF( frac_daylight + 2.0_wp * dusk_time > 1.0_wp) frac_daylight = 1.0_wp - 2.0_wp * dusk_time

    ! normal day
    IF(frac_daylight >= min_daylength) THEN
       ! calculate current time as fraction of day
       frac_time = ( time - NINT(time-0.5_wp) ) 

       ! mean daily temperature as average of minimum and maximum temperature
       dtmp_max = dtmp + 0.56_wp * diurnal_temp_range 
       dtmp_min = dtmp_max - diurnal_temp_range
  
       ! relative times of sunrise and sunset
       sunrise = 0.5_wp - frac_daylight / 2.0_wp
       sunset = sunrise + frac_daylight 

       ! daytime
       IF(frac_time > sunrise .AND. frac_time < sunset ) THEN 

          sd = COS(pi*(frac_time-tmp_peak_time)/(frac_daylight/2.0_wp+dusk_time))

          IF(sd >= 0.0_wp) THEN
             temperature = dtmp + sd * (dtmp_max - dtmp)
          ELSE 
             temperature = dtmp + sd * (dtmp - dtmp_min)
          ENDIF

       ELSE ! night-time: linear decline of temperatures since sunset 

          sd = COS(pi*(sunset-tmp_peak_time)/(frac_daylight/2.0_wp+dusk_time))

          IF(sd >= 0) THEN
             arg = dtmp + sd * (dtmp_max - dtmp)
          ELSE 
             arg = dtmp + sd * (dtmp - dtmp_min)
           ENDIF

           time_since_sunset = MOD(frac_time-sunset+1._wp,1._wp)
           temperature = dtmp_min + (arg - dtmp_min) * & 
               (1._wp - time_since_sunset/(1._wp - frac_daylight))

       ENDIF 

    ELSE 

       temperature = dtmp

    ENDIF

  END FUNCTION diurnal_from_daily_temperature

END MODULE mo_atmland_forcing
