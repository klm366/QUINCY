!>
!! @brief dummy forcing routine for generic forcing 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2016-01-27)
!! 2nd version Jan Engel (20160825) moved memory from SLM to JSBACH4 style
!!
MODULE mo_atmland_forcing

  USE mo_kind
  USE mo_shared_isotope,        ONLY: calc_mixing_ratio_C13C12, calc_mixing_ratio_C14C
  USE mo_atmland_util,          ONLY: update_orbital_parameters, calc_solar_angle, calc_sw_fdiffuse, calc_spec_humidity_sat
  USE mo_atmland_constants


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_forcing, open_forcing_input_file, close_forcing_input_file


CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> update forcing variables per timestep & read forcing data per timestep from file climate.dat
  !! 
  !! called in the time loop in plant_model_driver subsequently to advancing the timer
  !! 
  !! read one line from file per timestep \n
  !! if file ends, start again at top of file with reading lines \n
  !! parameter units are re-scaled accordingly to the needs of SLM \n
  !! the two header-lines may be neglected \n
  !! write column headers of output files for fluxes \n
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_forcing(      time               , &  ! in
                                  t_air              , &  ! out
                                  q_air              , &  ! out
                                  rain               , &  ! out
                                  snow               , &  ! out
                                  wind_air           , &  ! out
                                  wind_10m           , &  ! out
                                  lw_srf_down        , &  ! out
                                  swvis_srf_down     , &  ! out
                                  swnir_srf_down     , &  ! out
                                  swpar_srf_down     , &  ! out
                                  swpar_fdiffuse     , &  ! out
                                  press_srf          , &  ! out
                                  drag_srf           , &  ! out
                                  t_acoef            , &  ! out
                                  t_bcoef            , &  ! out
                                  q_acoef            , &  ! out
                                  q_bcoef            , &  ! out
                                  pch                , &  ! out
                                  cos_zenith_angle   , &  ! out
                                  co2_mixing_ratio     , &! out
                                  co2_mixing_ratio_C13 , &! out
                                  co2_mixing_ratio_C14 , &! out
                                  ga                   , &! out
                                  ! For lakes:
                                  albedo_lwtr        , &  ! out
                                  drag_wtr           , &  ! out
                                  drag_ice           , &  ! out
                                  t_acoef_wtr        , &  ! out
                                  t_bcoef_wtr        , &  ! out
                                  q_acoef_wtr        , &  ! out
                                  q_bcoef_wtr        , &  ! out
                                  t_acoef_ice        , &  ! out
                                  t_bcoef_ice        , &  ! out
                                  q_acoef_ice        , &  ! out
                                  q_bcoef_ice          &  ! out
                                  &  )

    USE mo_base_memory,   ONLY: grid  ! used for calc_solar_angle()
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp)          , INTENT(in)  :: time     !< current time (not timestep !)
    REAL(wp), OPTIONAL, INTENT(out) :: &  !! these are all variables from atm2land memory
                                      t_air              , &    !< ..
                                      q_air              , &    !< ..
                                      rain               , &    !< ..
                                      snow               , &    !< ..
                                      wind_air           , &    !< ..
                                      wind_10m           , &    !< ..
                                      lw_srf_down        , &    !< ..
                                      swvis_srf_down     , &    !< ..
                                      swnir_srf_down     , &    !< ..
                                      swpar_srf_down     , &    !< ..
                                      swpar_fdiffuse     , &    !< ..
                                      press_srf          , &    !< ..
                                      drag_srf           , &    !< ..
                                      t_acoef            , &    !< ..
                                      t_bcoef            , &    !< ..
                                      q_acoef            , &    !< ..
                                      q_bcoef            , &    !< ..
                                      pch                , &    !< ..
                                      cos_zenith_angle   , &    !< ..
                                      co2_mixing_ratio     , &  !< ..
                                      co2_mixing_ratio_C13 , &  !< ..
                                      co2_mixing_ratio_C14 , &  !< ..
                                      ga                   , &  !< ..
                                      ! For lakes:
                                      albedo_lwtr        , &    !< ..
                                      drag_wtr           , &    !< ..
                                      drag_ice           , &    !< ..
                                      t_acoef_wtr        , &    !< ..
                                      t_bcoef_wtr        , &    !< ..
                                      q_acoef_wtr        , &    !< ..
                                      q_bcoef_wtr        , &    !< ..
                                      t_acoef_ice        , &    !< ..
                                      t_bcoef_ice        , &    !< ..
                                      q_acoef_ice        , &    !< ..
                                      q_bcoef_ice               !< ..
    ! 0.2 Local
    INTEGER  :: stat          ! ..
    INTEGER  :: lun           ! logical unit number for file input stream, see '2.0 read data'
    REAL(wp) :: sw_toa_down   ! potential radiation at the top of the atmosphere
    REAL(wp) :: year          ! dummy variable
    REAL(wp) :: day           ! dummy variable
    REAL(wp) :: hour          ! dummy variable
    REAL(wp) :: nhx           ! dummy variable for nhx_srf_down
    REAL(wp) :: nox           ! dummy variable for noy_srf_down
    REAL(wp) :: hlp1          ! dummy variable 

    ! give atm2land variables a value to avoid compiler warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value
    t_air           = 0.0_wp
    q_air           = 0.0_wp
    rain            = 0.0_wp
    snow            = 0.0_wp
    wind_air        = 0.0_wp
    wind_10m        = 0.0_wp
    lw_srf_down     = 0.0_wp
    swvis_srf_down  = 0.0_wp
    swnir_srf_down  = 0.0_wp
    swpar_srf_down  = 0.0_wp
    swpar_fdiffuse  = 0.0_wp
    press_srf       = 0.0_wp
    drag_srf        = 1.0_wp    ! 1.0_wp
    t_acoef         = 0.0_wp
    t_bcoef         = 0.0_wp
    q_acoef         = 0.0_wp
    q_bcoef         = 0.0_wp
    pch             = 0.0_wp
    cos_zenith_angle      = 0.0_wp
    co2_mixing_ratio      = 0.0_wp
    co2_mixing_ratio_C13  = 0.0_wp
    co2_mixing_ratio_C14  = 0.0_wp
    ga                    = 0.0_wp
    ! For lakes:
    albedo_lwtr     = 0.0_wp
    drag_wtr        = 0.0_wp
    drag_ice        = 0.0_wp
    t_acoef_wtr     = 0.0_wp
    t_bcoef_wtr     = 0.0_wp
    q_acoef_wtr     = 0.0_wp
    q_bcoef_wtr     = 0.0_wp
    t_acoef_ice     = 0.0_wp
    t_bcoef_ice     = 0.0_wp
    q_acoef_ice     = 0.0_wp
    q_bcoef_ice     = 0.0_wp


    !> 1.0 orbital parameters
    !!
    !! 'declination' is calculated as output paramater, but not needed here
    CALL update_orbital_parameters(time,declination,sw_toa_down)  ! mo_atmland_util

    !> 2.0 read data
    !!
    !! the read statement, by default, reads a list of values in a line and then advances to the beginning of the next line
    lun=9
    READ(lun,*, IOSTAT=stat) year, day, hour, & 
                     swvis_srf_down,   lw_srf_down, t_air,    q_air,    &
                     press_srf,        rain,        snow,     wind_air, &
                     co2_mixing_ratio,                                                             &
                     nhx, nox
    IF(IS_IOSTAT_END(stat)) THEN  ! if end of file (IS_IOSTAT_END(stat) == Test for end-of-file value)
      ! start reading from line one & read header into dummy
      REWIND(lun)
      CALL read_forcing_header
      READ(lun,*) year, day, hour, &
                     swvis_srf_down,   lw_srf_down, t_air,    q_air,    &
                     press_srf,        rain,        snow,     wind_air, &
                     co2_mixing_ratio,                                                             &
                     nhx, nox
    END IF

    !> 3.0 parameter units are re-scaled
    !!
    q_air     = q_air / 1000._wp
    press_srf = press_srf * 100._wp
    rain      = rain / 86400._wp
    snow      = snow / 86400._wp

    !> 4.0 update forcing variables based on 'read data'
    !!
    wind_10m             = wind_air
    ga                   = 0.04_wp / drag_srf * wind_10m
    hlp1                 = calc_mixing_ratio_C13C12(def_co2_deltaC13)                 ! function: calc_mixing_ratio_C13C12
    co2_mixing_ratio_C13 = co2_mixing_ratio / (1._wp + 1._wp/hlp1)
    hlp1                 = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14) ! function: calc_mixing_ratio_C14C
    co2_mixing_ratio_C14 = hlp1 * co2_mixing_ratio 

    cos_zenith_angle = calc_solar_angle(grid%lat,grid%lon,time)                       ! function: calc_solar_angle
    swnir_srf_down   = (1._wp - frac_vis_swtot_srf ) * swvis_srf_down 
    swvis_srf_down   = frac_vis_swtot_srf * swvis_srf_down 
    swpar_srf_down   = frac_par_swvis_srf * swvis_srf_down
    swpar_fdiffuse   = calc_sw_fdiffuse(swvis_srf_down + swnir_srf_down, sw_toa_down, cos_zenith_angle) ! function: calc_sw_fdiffuse

  END SUBROUTINE update_forcing


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
  END SUBROUTINE open_forcing_input_file


  !-----------------------------------------------------------------------------------------------------
  !> closes the file test_bed/plant/climate.dat, only IF model option is set to READ
  !! 
  !! called in the driver
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE close_forcing_input_file

    IMPLICIT NONE

      CLOSE(9)

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

END MODULE mo_atmland_forcing
