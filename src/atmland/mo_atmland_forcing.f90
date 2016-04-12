!>
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
  USE mo_atmland_memory, ONLY: atm2land,land2atm
  USE mo_atmland_util
  USE mo_atmland_constants
  USE mo_atmland_config_class,  ONLY: atmland_config


  IMPLICIT NONE

  REAL(wp) :: sw_toa_down                            ! potential radiation at the top of the atmosphere
  REAL(wp) :: annual_t_air_min, annual_t_air_max     ! annual minimum and maximum daily t_air
  REAL(wp) :: diurnal_temp_range                     ! diurnal temperature range


  PRIVATE

  PUBLIC :: init_forcing,update_forcing,open_forcing_input_file,close_forcing_input_file


CONTAINS

  SUBROUTINE init_forcing
    ! docu
    USE mo_base_memory
    IMPLICIT NONE

    ! Initialization
    annual_t_air_min = -5.0_wp + Tzero       ! annual minumum average day temperature
    annual_t_air_max = 25.0_wp + Tzero       ! annual maximum average day temperature
    diurnal_temp_range = .1_wp               ! diurnal temperature range
    ! 
    ! Initalisation of met-forcing to standard values
    ! 
    atm2land%t_air = def_t_air
    atm2land%press_srf = standard_press_srf
    atm2land%q_air = calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf) - def_vpd / atm2land%press_srf * eps_vpd
    atm2land%co2_mixing_ratio = def_co2
    atm2land%rain = 0.0_wp
    atm2land%snow = 0.0_wp
    atm2land%swvis_srf_down = def_swdown
    atm2land%swnir_srf_down = atm2land%swvis_srf_down * 0.1_wp
    atm2land%swpar_srf_down = atm2land%swvis_srf_down
    atm2land%swpar_fdiffuse = def_fdiffuse
    atm2land%wind_air = def_wind
    atm2land%wind_10m = atm2land%wind_air
    atm2land%lw_srf_down = 0.0_wp
    atm2land%drag_srf = 1.0_wp
    atm2land%t_acoef = 0.0_wp
    atm2land%t_bcoef = 0.0_wp
    atm2land%q_acoef = 0.0_wp
    atm2land%q_bcoef = 0.0_wp
    atm2land%pch = 0.0_wp
    atm2land%cos_zenith_angle = def_cos_angle
    atm2land%ga = .1_wp / ( atm2land%drag_srf * atm2land%wind_10m )

    grid%lat = atmland_config%lat
    grid%lon = atmland_config%lon

  END SUBROUTINE init_forcing

  SUBROUTINE update_forcing(time)
    ! called in the time loop in plant_model_driver subsequently to advancing the timer
    USE mo_base_memory
    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp),INTENT(in) :: time 

     ! 
     ! 2.2 update time-variant forcing
     ! 
     SELECT CASE(atmland_config%option_forcing)  ! see mo_base_config
       CASE('constant')
         ! default input values are not changed ?
       CASE('readfile')
         ! read climate data from file climate.dat in test_bed/plant
         ! climate.dat is a copy from /Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET/file...
         CALL update_orbital_parameters(time,declination,sw_toa_down)
         CALL read_forcing
         atm2land%cos_zenith_angle = calc_solar_angle(grid%lat,grid%lon,time)
         atm2land%swnir_srf_down = atm2land%swvis_srf_down * 0.1_wp
         atm2land%swpar_srf_down = atm2land%swvis_srf_down
         atm2land%swpar_fdiffuse = calc_sw_fdiffuse(atm2land%swvis_srf_down, sw_toa_down, atm2land%cos_zenith_angle)
       CASE('generate')
         CALL update_orbital_parameters(time,declination,sw_toa_down)
         atm2land%t_air = generate_inst_temperature(grid%lat,time,declination,annual_t_air_min, &
                          annual_t_air_max,diurnal_temp_range)
         atm2land%cos_zenith_angle = calc_solar_angle(grid%lat,grid%lon,time)
         atm2land%swvis_srf_down = MAX(0.0, def_swdown * atm2land%cos_zenith_angle)
         atm2land%swnir_srf_down = atm2land%swvis_srf_down * 0.1_wp
         atm2land%swpar_srf_down = atm2land%swvis_srf_down
         atm2land%swpar_fdiffuse = calc_sw_fdiffuse(atm2land%swvis_srf_down, sw_toa_down, atm2land%cos_zenith_angle)
         atm2land%q_air = def_rh * calc_spec_humidity_sat(atm2land%t_air, atm2land%press_srf)
     END SELECT

  END SUBROUTINE update_forcing

  SUBROUTINE read_forcing
    ! read one line from file per timestep
    ! if file ends, start again at top of file with reading lines
    ! parameter units are re-scaled accordingly to the needs of SLM
    ! the two header-lines may be neglected
    ! write column headers of output files for fluxes
    IMPLICIT NONE
    ! 0.3 Local
    INTEGER  :: stat  = 0
    INTEGER  :: lun   = 100
    REAL(wp) :: year  = 0  ! dummy variable
    REAL(wp) :: day   = 0  ! dummy variable
    REAL(wp) :: hour  = 0  ! dummy variable
    REAL(wp) :: nhx   = 0  ! dummy variable for nhx_srf_down
    REAL(wp) :: nox   = 0  ! dummy variable for noy_srf_down
    ! 1.0 Set write options
    lun=100
    ! the read statement, by default, reads a list of values in a line and then advances to the beginning of the next line
    READ(lun,*, IOSTAT=stat) year, day, hour, & 
                     atm2land%swvis_srf_down,   atm2land%lw_srf_down, atm2land%t_air,    atm2land%q_air,    &
                     atm2land%press_srf,        atm2land%rain,        atm2land%snow,     atm2land%wind_air, &
                     atm2land%co2_mixing_ratio,                                                             &
                     nhx, nox
    IF(IS_IOSTAT_END(stat)) THEN  ! if end of file (IS_IOSTAT_END(stat) == Test for end-of-file value)
      ! start reading from line one & read header into dummy
      REWIND(lun)
      CALL read_forcing_header
      READ(lun,*) year, day, hour, &
                     atm2land%swvis_srf_down,   atm2land%lw_srf_down, atm2land%t_air,    atm2land%q_air,    &
                     atm2land%press_srf,        atm2land%rain,        atm2land%snow,     atm2land%wind_air, &
                     atm2land%co2_mixing_ratio,                                                             &
                     nhx, nox
    END IF

    ! parameter units are re-scaled
    atm2land%q_air     = atm2land%q_air / 1000._wp
    atm2land%press_srf = atm2land%press_srf * 100._wp
    atm2land%rain      = atm2land%rain / 86400._wp
    atm2land%snow      = atm2land%snow / 86400._wp

  END SUBROUTINE read_forcing

  SUBROUTINE open_forcing_input_file
    ! open the file test_bed/plant/climate.dat
    ! only IF model option is set to READ !
    USE mo_exception, ONLY: finish
    IMPLICIT NONE
    ! 0.3 local
    LOGICAL :: file_exists
    INTEGER :: lun = 100
    CHARACTER(100) :: fname = "filename"
    CHARACTER(len=*), PARAMETER :: routine = 'mo_atmland_forcing:open_forcing_input_file'
    ! 1.0 open file for READ
    IF(atmland_config%option_forcing == 'readfile') THEN
      lun=100
      fname = 'climate.dat'
      ! check if file exists
      INQUIRE(FILE=fname, EXIST=file_exists)
      ! open if fname == True
      IF(file_exists) THEN
        OPEN (UNIT=lun, FILE=fname, ACTION="read", STATUS="old") ! status=old: file exists already
        ! 2.0 read the two header lines into dummy variables
        CALL read_forcing_header
      ELSE
        CALL finish(TRIM(routine), 'climate forcing: climate.dat could not be opened, file does not exist')
      END IF
    END IF
  END SUBROUTINE open_forcing_input_file

  SUBROUTINE close_forcing_input_file
    ! simply close the file climate.dat
    IMPLICIT NONE
    IF(atmland_config%option_forcing == 'readfile') THEN
      CLOSE(100)
    END IF
  END SUBROUTINE close_forcing_input_file

  SUBROUTINE read_forcing_header
    ! read only the header into a dummy variable to advanve the READ to the 3rd line in the file
    IMPLICIT NONE
    ! 0.3 Local
    CHARACTER(200) :: header_dummy1 = "header"
    CHARACTER(200) :: header_dummy2 = "header"
    INTEGER :: lun = 100
    ! 1.0 read header
    lun=100
    READ(lun,'(a)') header_dummy1   ! header line 1
    READ(lun,'(a)') header_dummy2   ! header line 2

  END SUBROUTINE read_forcing_header

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
