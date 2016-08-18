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
  USE mo_shared_isotope
  USE mo_atmland_memory,        ONLY: atm2land, land2atm
  USE mo_atmland_util
  USE mo_atmland_constants
  USE mo_atmland_config_class,  ONLY: atmland_config


  IMPLICIT NONE

  REAL(wp) :: sw_toa_down                            ! potential radiation at the top of the atmosphere
  REAL(wp) :: annual_t_air_min, annual_t_air_max     ! annual minimum and maximum daily t_air
  REAL(wp) :: diurnal_temp_range                     ! diurnal temperature range


  PRIVATE

  PUBLIC :: init_forcing, update_forcing, open_forcing_input_file, close_forcing_input_file


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> intialize forcing variables
  !! 
  !! set standard values for atm2land, plus some more variables
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_forcing

    USE mo_base_memory,             ONLY: grid
    USE mo_shared_constants,        ONLY: Tzero

    IMPLICIT NONE

    REAL(wp) :: C13C12_mixing_ratio, &
                C14C_mixing_ratio

    !> Initialization of temperatures
    !! 
    !! 
    annual_t_air_min          = -5.0_wp + Tzero       ! annual minumum average day temperature
    annual_t_air_max          = 25.0_wp + Tzero       ! annual maximum average day temperature
    diurnal_temp_range        = .1_wp                 ! diurnal temperature range

    !> Initalisation of met-forcing to standard values
    !! 
    !! 
    atm2land%t_air                = def_t_air
    atm2land%press_srf            = standard_press_srf
    atm2land%q_air                = calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf) & 
                                    - def_vpd / atm2land%press_srf * eps_vpd
    atm2land%co2_mixing_ratio     = def_co2
    C13C12_mixing_ratio           = calc_mixing_ratio_C13C12(def_co2_deltaC13) 
    atm2land%co2_mixing_ratio_c13 = atm2land%co2_mixing_ratio / (1._wp + 1._wp/C13C12_mixing_ratio)
    C14C_mixing_ratio             = calc_mixing_ratio_C14C(def_co2_deltaC13,def_co2_deltaC14)
    atm2land%co2_mixing_ratio_c14 = C14C_mixing_ratio * atm2land%co2_mixing_ratio 
    atm2land%rain                 = 0.0_wp
    atm2land%snow                 = 0.0_wp
    atm2land%swvis_srf_down       = def_swdown * frac_vis_swtot_srf
    atm2land%swnir_srf_down       = def_swdown * (1._wp - frac_vis_swtot_srf)
    atm2land%swpar_srf_down       = atm2land%swvis_srf_down * frac_par_swvis_srf
    atm2land%swpar_fdiffuse       = def_fdiffuse
    atm2land%wind_air             = def_wind
    atm2land%wind_10m             = atm2land%wind_air
    atm2land%lw_srf_down          = 0.0_wp
    atm2land%drag_srf             = 1.0_wp
    atm2land%t_acoef              = 0.0_wp
    atm2land%t_bcoef              = 0.0_wp
    atm2land%q_acoef              = 0.0_wp
    atm2land%q_bcoef              = 0.0_wp
    atm2land%pch                  = 0.0_wp
    atm2land%cos_zenith_angle     = def_cos_angle
    atm2land%ga                   = 0.1_wp / ( atm2land%drag_srf * atm2land%wind_10m )

    !> Initialization of of grid\% lon and lat
    !! 
    !! 
    grid%lat                  = atmland_config%lat
    grid%lon                  = atmland_config%lon
  END SUBROUTINE init_forcing

  !-----------------------------------------------------------------------------------------------------
  !> update forcing variables per timestep
  !! 
  !! called in the time loop in plant_model_driver subsequently to advancing the timer
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_forcing(time)

    USE mo_base_memory,   ONLY: grid
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),INTENT(in) :: time     !< current time (not timestep !)

   !> 1.0 update time-variant forcing
   !!
   !! as one of three cases specified in the namelist: constant, readfile, generate
   SELECT CASE(atmland_config%option_forcing)
     CASE('constant')
       ! default input values are not changed ?
     CASE('readfile')
       ! read climate data from file climate.dat in test_bed/plant
       ! climate.dat is a copy from /Net/Groups/BGI/data/Model/SLM1/input/point/FLUXNET/file...
       CALL update_orbital_parameters(time,declination,sw_toa_down)
       CALL read_forcing
       atm2land%cos_zenith_angle = calc_solar_angle(grid%lat,grid%lon,time)
       atm2land%swnir_srf_down   = (1._wp - frac_vis_swtot_srf ) * atm2land%swvis_srf_down 
       atm2land%swvis_srf_down   = frac_vis_swtot_srf * atm2land%swvis_srf_down 
       atm2land%swpar_srf_down   = frac_par_swvis_srf * atm2land%swvis_srf_down
       atm2land%swpar_fdiffuse   = calc_sw_fdiffuse(atm2land%swvis_srf_down+atm2land%swnir_srf_down, & 
                                      sw_toa_down, atm2land%cos_zenith_angle)
     CASE('generate')
       CALL update_orbital_parameters(time,declination,sw_toa_down)
       atm2land%t_air            = generate_inst_temperature(grid%lat,time,declination,annual_t_air_min, &
                                   annual_t_air_max,diurnal_temp_range)
       atm2land%cos_zenith_angle = calc_solar_angle(grid%lat,grid%lon,time)
       atm2land%swvis_srf_down   = frac_vis_swtot_srf * MAX(0.0_wp, def_swdown * atm2land%cos_zenith_angle)
       atm2land%swnir_srf_down   = (1._wp - frac_vis_swtot_srf ) * MAX(0.0_wp, def_swdown * atm2land%cos_zenith_angle) 
       atm2land%swpar_srf_down   = frac_par_swvis_srf * atm2land%swvis_srf_down
       atm2land%swpar_fdiffuse   = calc_sw_fdiffuse(atm2land%swvis_srf_down+atm2land%swnir_srf_down, & 
                                      sw_toa_down, atm2land%cos_zenith_angle)
       atm2land%q_air            = def_rh * calc_spec_humidity_sat(atm2land%t_air, atm2land%press_srf)
   END SELECT
  END SUBROUTINE update_forcing

  !-----------------------------------------------------------------------------------------------------
  !> read forcing data per timestep from file climate.dat
  !! 
  !! read one line from file per timestep \n
  !! if file ends, start again at top of file with reading lines \n
  !! parameter units are re-scaled accordingly to the needs of SLM \n
  !! the two header-lines may be neglected \n
  !! write column headers of output files for fluxes \n
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE read_forcing

    IMPLICIT NONE
    ! 0.2 Local
    INTEGER  :: stat 
    INTEGER  :: lun  
    REAL(wp) :: year      ! dummy variable
    REAL(wp) :: day       ! dummy variable
    REAL(wp) :: hour      ! dummy variable
    REAL(wp) :: nhx       ! dummy variable for nhx_srf_down
    REAL(wp) :: nox       ! dummy variable for noy_srf_down
    REAL(wp) :: hlp1      ! dummy variable 

    ! 1.0 Set read options
    lun=9
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
    atm2land%wind_10m  = atm2land%wind_air
    atm2land%ga        = 0.1_wp / ( atm2land%drag_srf * atm2land%wind_10m )
    hlp1 = calc_mixing_ratio_C13C12(def_co2_deltaC13)
    atm2land%co2_mixing_ratio_c13 = atm2land%co2_mixing_ratio / (1._wp + 1._wp/hlp1)
    hlp1 = calc_mixing_ratio_C14C(def_co2_deltaC13,def_co2_deltaC14)
    atm2land%co2_mixing_ratio_c14 = hlp1 * atm2land%co2_mixing_ratio 
  END SUBROUTINE read_forcing

  !-----------------------------------------------------------------------------------------------------
  !> open the file test_bed/plant/climate.dat, only IF model option is set to READ
  !! 
  !! called in the driver
  !! 
  !! note: lun is set to 9 (see ticket 31)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE open_forcing_input_file

    USE mo_exception, ONLY: finish

    IMPLICIT NONE
    ! 0.3 local
    LOGICAL :: file_exists
    INTEGER :: lun 
    CHARACTER(100) :: fname 
    CHARACTER(len=*), PARAMETER :: routine = 'mo_atmland_forcing:open_forcing_input_file'

    !> 1.0 open file for READ
    IF(atmland_config%option_forcing == 'readfile') THEN
      lun=9
      fname = 'climate.dat'
      ! check if file exists
      INQUIRE(FILE=fname, EXIST=file_exists)
      ! open if fname == True
      IF(file_exists) THEN
        OPEN (UNIT=lun, FILE=fname, ACTION="read", STATUS="old") ! status=old: file exists already
    !> 2.0 read the two header lines into dummy variables
        CALL read_forcing_header
      ELSE
        CALL finish(TRIM(routine), 'climate forcing: climate.dat could not be opened, file does not exist')
      END IF
    END IF
  END SUBROUTINE open_forcing_input_file

  !-----------------------------------------------------------------------------------------------------
  !> closes the file test_bed/plant/climate.dat, only IF model option is set to READ
  !! 
  !! called in the driver
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE close_forcing_input_file

    IMPLICIT NONE
    IF(atmland_config%option_forcing == 'readfile') THEN
      CLOSE(9)
    END IF
  END SUBROUTINE close_forcing_input_file

  !-----------------------------------------------------------------------------------------------------
  !> read only the header into a dummy variable to advanve the READ to the 3rd line in the file
  !! 
  !! consider this a helper subroutine to the read_forcing subroutine
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE read_forcing_header

    IMPLICIT NONE
    ! 0.2 Local
    CHARACTER(200) :: header_dummy1 = "header"
    CHARACTER(200) :: header_dummy2 = "header"
    INTEGER :: lun 

    !> 1.0 read header
    lun=9
    READ(lun,'(a)') header_dummy1   ! header line 1
    READ(lun,'(a)') header_dummy2   ! header line 2
  END SUBROUTINE read_forcing_header

  !-----------------------------------------------------------------------------------------------------
  !> generate instantaneous temperature from annual characteristics 
  !! 
  !! note: meant to be used for generic testing only. Only works for present-day
  !-----------------------------------------------------------------------------------------------------
  FUNCTION generate_inst_temperature &
              (latitude ,time ,declination ,annual_t_air_min ,annual_t_air_max ,diurnal_temp_range) &
              RESULT(temperature)

    USE mo_atmland_constants, ONLY: solstice
    USE mo_atmland_util, ONLY: update_orbital_parameters

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),INTENT(in) :: latitude            !< latitude (degree)
    REAL(wp),INTENT(in) :: time                !< time since 1st Januar (days)
    REAL(wp),INTENT(in) :: declination         !< solar declination (radians?)  
    REAL(wp),INTENT(in) :: annual_t_air_min    !< annual minimum daily temperature (K)
    REAL(wp),INTENT(in) :: annual_t_air_max    !< annual maximum daily temperature (K)
    REAL(wp),INTENT(in) :: diurnal_temp_range  !< diurnal temperature range (K)
    REAL(wp) :: temperature                    !< current air temperature (K)
    ! 0.2 Local
    REAL(wp) :: declination_max                ! declination at solstice 
    REAL(wp) :: dtmp                           ! today's temperature (K)
    REAL(wp) :: hlp1                           ! dummy variable

    CALL update_orbital_parameters(solstice,declination_max,hlp1)

    dtmp = annual_t_air_min + (annual_t_air_max - annual_t_air_min) * & 
                  ( 1.0_wp - ( declination_max - declination ) / ( 2.0_wp * declination_max ))

    temperature = diurnal_from_daily_temperature(latitude,time,declination,dtmp,diurnal_temp_range)
  END FUNCTION generate_inst_temperature

  !-----------------------------------------------------------------------------------------------------
  !> calculate the temperature at a given time of the day 
  !! 
  !! calculations consider: daily minimum and maximum temperatures, solar declination and time of the day
  !-----------------------------------------------------------------------------------------------------
  FUNCTION diurnal_from_daily_temperature &
               (latitude,time,declination,dtmp,diurnal_temp_range) &
               RESULT(temperature)

    USE mo_base_constants, ONLY: deg2rad,pi
    USE mo_atmland_constants, ONLY: min_daylength,tmp_peak_time,dusk_time

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),INTENT(in) :: latitude             !< latitude (degree)
    REAL(wp),INTENT(in) :: time                 !< time since 1st Januar (days)
    REAL(wp),INTENT(in) :: declination          !< solar declination (radians?)  
    REAL(wp),INTENT(in) :: dtmp                 !< daily mean temperature (K)
    REAL(wp),INTENT(in) :: diurnal_temp_range   !< daily temperature range (K)
    REAL(wp) :: temperature                     !< current air temperature 
    ! 0.2 Local
    REAL(wp) :: frac_time                       ! time of day as fraction 
    REAL(wp) :: spds,cpds,arg,sd                ! factors to calculate daylight hourse
    REAL(wp) :: frac_daylight                   ! fraction of the day that is light
    REAL(wp) :: sunrise,sunset                  ! fractional time of the day at which sunrise and sunset occurs
    REAL(wp) :: time_since_sunset               ! helper variable (temperature at sunset)
    REAL(wp) :: dtmp_min                        ! daily minimum temperature (K)
    REAL(wp) :: dtmp_max                        ! daily maximum temperature (K)

    !> work out today's daylength
    !!
    !!
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

    !> IF: normal day
    !!
    !!
    IF(frac_daylight >= min_daylength) THEN
       ! calculate current time as fraction of day
       frac_time = ( time - NINT(time-0.5_wp) ) 

       ! mean daily temperature as average of minimum and maximum temperature
       dtmp_max = dtmp + 0.56_wp * diurnal_temp_range 
       dtmp_min = dtmp_max - diurnal_temp_range
  
       ! relative times of sunrise and sunset
       sunrise = 0.5_wp - frac_daylight / 2.0_wp
       sunset = sunrise + frac_daylight 

       !> IF: daytime
       IF(frac_time > sunrise .AND. frac_time < sunset ) THEN 

          sd = COS(pi*(frac_time-tmp_peak_time)/(frac_daylight/2.0_wp+dusk_time))

          IF(sd >= 0.0_wp) THEN
             temperature = dtmp + sd * (dtmp_max - dtmp)
          ELSE 
             temperature = dtmp + sd * (dtmp - dtmp_min)
          ENDIF

       ELSE !> ELSE: night-time: linear decline of temperatures since sunset 

          sd = COS(pi*(sunset-tmp_peak_time)/(frac_daylight/2.0_wp+dusk_time))

          IF(sd >= 0) THEN
             arg = dtmp + sd * (dtmp_max - dtmp)
          ELSE 
             arg = dtmp + sd * (dtmp - dtmp_min)
           ENDIF

           time_since_sunset = MOD(frac_time-sunset+1._wp,1._wp)
           temperature = dtmp_min + (arg - dtmp_min) * & 
               (1._wp - time_since_sunset/(1._wp - frac_daylight))

       END IF !> END IF

    !> ELSE: temperature = dtmp
    !!
    !!
    ELSE 
       temperature = dtmp

    ENDIF
  END FUNCTION diurnal_from_daily_temperature

END MODULE mo_atmland_forcing
