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
  USE mo_base_configuration
  USE mo_land_constants
  USE mo_atmland_memory
  USE mo_atmland_util
  USE mo_atmland_constants


  IMPLICIT NONE

  TYPE(t_atm2land_memory), SAVE           :: atm2land         ! atmostphere to land 
  TYPE(t_land2atm_memory), SAVE           :: land2atm         ! land to atmosphere

  REAL(wp) :: time                         ! displays the 'real time' in days; time == 1.0 after one day has passed
  INTEGER :: ntimestep
  REAL(wp) :: sw_toa_down                            ! potential radiation at the top of the atmosphere
  REAL(wp) :: longitude, latitude                    ! coordinates of point to be simulated
  REAL(wp) :: annual_t_air_min, annual_t_air_max     ! annual minimum and maximum daily t_air
  REAL(wp) :: diurnal_temp_range                     ! diurnal temperature range

  REAL(wp) :: ref_t_air = 25.0_wp          ! temperature in deg C
  REAL(wp) :: ref_vpd = 1000._wp           ! vpd in Pa
  REAL(wp) :: ref_rh = 0.8_wp              ! relative humidity in fractions
  REAL(wp) :: ref_swdown = 652.173913_wp   ! PAR in W m-2
  REAL(wp) :: ref_fdiffuse = 0.2_wp        ! diffuse fraction
  REAL(wp) :: ref_cos_angle = 1.0_wp       ! cosine of zenith angle
  REAL(wp) :: ref_co2 = 380.0_wp           ! co2 in ppm
  REAL(wp) :: ref_wind = 1.0_wp            ! wind in m/s
  REAL(wp) :: ref_theta = 100._wp          ! soil water availability in %
  REAL(wp) :: ref_lai = 3._wp              ! LAI
  REAL(wp) :: ref_nleaf = 2.0_wp           ! N leaf per unit leaf area in gN/m2


CONTAINS

  SUBROUTINE init_atmland_forcing
    ! docu
    IMPLICIT NONE

    ! Initialization
    annual_t_air_min = -5.0_wp + Tzero       ! annual minumum average day temperature
    annual_t_air_max = 25.0_wp + Tzero       ! annual maximum average day temperature
    diurnal_temp_range = .1_wp               ! diurnal temperature range
    ! 
    ! Initalisation of met-forcing to standard values
    ! 
    atm2land%t_air = ref_t_air + Tzero
    atm2land%press_srf = standard_press_srf
    atm2land%q_air = calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf) - ref_vpd / atm2land%press_srf * eps_vpd
    atm2land%co2_mixing_ratio = ref_co2
    atm2land%rain = 0.0_wp
    atm2land%snow = 0.0_wp
    atm2land%swvis_srf_down = ref_swdown
    atm2land%swnir_srf_down = atm2land%swvis_srf_down * 0.1_wp
    atm2land%swpar_srf_down = atm2land%swvis_srf_down
    atm2land%swpar_fdiffuse = ref_fdiffuse
    atm2land%wind_air = ref_wind
    atm2land%wind_10m = atm2land%wind_air
    atm2land%lw_srf_down = 0.0_wp
    atm2land%drag_srf = 1.0_wp
    atm2land%t_acoef = 0.0_wp
    atm2land%t_bcoef = 0.0_wp
    atm2land%q_acoef = 0.0_wp
    atm2land%q_bcoef = 0.0_wp
    atm2land%pch = 0.0_wp
    atm2land%cos_zenith_angle = ref_cos_angle
    atm2land%ga = .1_wp / ( atm2land%drag_srf * atm2land%wind_10m )
  END SUBROUTINE init_atmland_forcing

  SUBROUTINE update_forcing
    ! called in the time loop in plant_model_driver subsequently to advancing the timer
    IMPLICIT NONE

     CALL update_orbital_parameters(time,declination,sw_toa_down)
     ! 
     ! 2.2 update time-variant forcing
     ! 
     SELECT CASE(options%forcing%option_forcing)  ! see mo_base_configuration
       CASE('constant')
         ! default input values are not changed ?
       CASE('readfile')
         ! implement later
       CASE('generate')
         atm2land%t_air = generate_inst_temperature(latitude,time,declination,annual_t_air_min, &
         annual_t_air_max,diurnal_temp_range)
         atm2land%cos_zenith_angle = calc_solar_angle(latitude,longitude,time)
         atm2land%swvis_srf_down = MAX(0.0, ref_swdown * atm2land%cos_zenith_angle)
         atm2land%swnir_srf_down = atm2land%swvis_srf_down * 0.1_wp
         atm2land%swpar_srf_down = atm2land%swvis_srf_down
         atm2land%swpar_fdiffuse = calc_sw_fdiffuse(atm2land%swvis_srf_down, sw_toa_down, atm2land%cos_zenith_angle)
         atm2land%q_air = ref_rh * calc_spec_humidity_sat(atm2land%t_air, atm2land%press_srf)
     END SELECT

  END SUBROUTINE update_forcing

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

  END FUNCTION generate_inst_temperature

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
