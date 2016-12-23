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

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp, eps8
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
  SUBROUTINE update_forcing(      tile               , &  ! inout
                                  time               , &  ! in
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

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
    USE mo_base_memory,           ONLY: t_slm_tile, t_model, grid  ! grid is used for calc_solar_angle()
    USE mo_base_config_class,     ONLY: base_config
    USE mo_soil_config_class,     ONLY: soil_config

    ! Use of process memories
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),   INTENT(inout) :: tile     !< one tile
    REAL(wp)          , INTENT(in)    :: time     !< current time (not timestep !)
    REAL(wp), OPTIONAL, INTENT(out)   :: &        !! these are all variables from atm2land memory
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
    INTEGER  :: stat          ! status of I/O (https://software.intel.com/en-us/node/691976)
    INTEGER  :: lun           ! logical unit number for file input stream, see '2.0 read data'
    REAL(wp) :: sw_toa_down   ! potential radiation at the top of the atmosphere
    REAL(wp) :: year          ! dummy variable
    REAL(wp) :: hour          ! dummy variable
    REAL(wp) :: nhx           ! dummy variable for nhx_srf_down
    REAL(wp) :: nox           ! dummy variable for noy_srf_down
    REAL(wp) :: hlp1          ! dummy variable 
    REAL(wp) :: time_hlp      ! dummy variable for soil forcing
    REAL(wp),         SAVE        :: day                    ! read input dummy variable - save value
    REAL(wp)                      :: doy_old, &             ! check for new year (column 2 in input file would change from 365 -> 1)
                                     doy_new
    TYPE(t_model),    POINTER     :: model                  ! the model structure
    INTEGER                       :: iblk         = -1      ! Number of current block (chunk)
    INTEGER                       :: ics          = -1      ! Index of chunk start
    INTEGER                       :: ice          = -1      ! Index of chunk end
    INTEGER                       :: nlevel_soil  = -1      ! number of levels at the Z-Axis (soil layer)
    INTEGER                       :: nlevel_veg   = -1      ! number of levels at the Z-Axis (canopy layer)
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(SOIL_)                       !     TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    TYPE(t_soil_memory),     POINTER :: soil        ! local pointer to soil__mem for SLM usability
    TYPE(t_veg_memory),      POINTER :: veg         ! local pointer to veg__mem for SLM usability
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(VEG_)

    soil    => soil__mem 
    veg     => veg__mem 
    model   => Get_model(tile%owner_model_id)


    doy_old = 0.0_wp  ! doy in the previous timestep, used for detecting a new year

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
    drag_srf        = 1000.0_wp    ! 1.0_wp
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


    !! LAND PLANT CANOPY
    SELECT CASE(TRIM(base_config%slm_model_name))
    CASE('land','plant','canopy')

    !> 1.0 orbital parameters
    !!
    !! 'declination' is calculated as output paramater, but not needed here
    CALL update_orbital_parameters(time,declination,sw_toa_down)  ! mo_atmland_util


    !> 2.0 read data
    !!
    !! read line by line from file until file end and start again at beginning of file \n
    !! if base_config%flag_climate_forcing_one_year = TRUE then 
    !!        read only the data of the 1st year of climate forcing and start again at beginning of file

    ! day of the year: previous time step
    doy_old = day 

    ! the read statement, by default a list of values in a line is read and then advanced to the beginning of the next line
    lun=9
    READ(lun,*, IOSTAT=stat) year, day, hour, & 
                     swvis_srf_down,   lw_srf_down, t_air,    q_air,    &
                     press_srf,        rain,        snow,     wind_air, &
                     co2_mixing_ratio,                                                             &
                     nhx, nox

    ! day of the year: current time step
    doy_new = day 
    
    ! if this flag is set TRUE, check whether the current time step is a new year, if so 
    ! set stat = -1 to indicate "end of file" which invokes the re-reading of the climate forcing input file at line 1
    IF(base_config%flag_climate_forcing_one_year) THEN
      IF( ABS(doy_old - doy_new - 364.0_wp) < eps8 ) THEN
        stat = -1
      END IF
    END IF

    ! if end of file (IS_IOSTAT_END(stat) == Test for end-of-file value)
    IF(IS_IOSTAT_END(stat)) THEN 
      ! start reading from line one & read header into dummy
      REWIND(lun)
      CALL read_forcing_header(lun)
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
    ga                   = 40.0_wp / drag_srf * wind_10m
    hlp1                 = calc_mixing_ratio_C13C12(def_co2_deltaC13)                 ! function: calc_mixing_ratio_C13C12
    co2_mixing_ratio_C13 = co2_mixing_ratio / (1._wp + 1._wp/hlp1)
    hlp1                 = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14) ! function: calc_mixing_ratio_C14C
    co2_mixing_ratio_C14 = hlp1 * co2_mixing_ratio 

    cos_zenith_angle = calc_solar_angle(grid%lat,grid%lon,time)                       ! function: calc_solar_angle
    swnir_srf_down   = (1._wp - frac_vis_swtot_srf ) * swvis_srf_down 
    swvis_srf_down   = frac_vis_swtot_srf * swvis_srf_down 
    swpar_srf_down   = frac_par_swvis_srf * swvis_srf_down
    swpar_fdiffuse   = calc_sw_fdiffuse(swvis_srf_down + swnir_srf_down, sw_toa_down, cos_zenith_angle) ! function: calc_sw_fdiffuse


    !! SOIL
    CASE('soil')

    iblk        = model%options%iblk
    ics         = model%options%ics
    ice         = model%options%ice
    nlevel_soil = soil_config%nsoil_water

    ! the read statement, by default a list of values in a line is read and then advanced to the beginning of the next line
    lun=9
    READ(lun,*, IOSTAT=stat) time_hlp, & 
        ! leaf
        veg%litterfall%leaf%carbon%ptr, veg%litterfall%leaf%nitrogen%ptr, veg%litterfall%leaf%phosphorus%ptr    , &
        veg%litterfall%leaf%carbon13%ptr, veg%litterfall%leaf%carbon14%ptr, veg%litterfall%leaf%nitrogen15%ptr  , &
        ! fine_root
        veg%litterfall%fine_root%carbon%ptr, veg%litterfall%fine_root%nitrogen%ptr, veg%litterfall%fine_root%phosphorus%ptr   , &
        veg%litterfall%fine_root%carbon13%ptr, veg%litterfall%fine_root%carbon14%ptr, veg%litterfall%fine_root%nitrogen15%ptr , &
        ! coarse_root
        veg%litterfall%coarse_root%carbon%ptr, veg%litterfall%coarse_root%nitrogen%ptr, veg%litterfall%coarse_root%phosphorus%ptr   , &
        veg%litterfall%coarse_root%carbon13%ptr, veg%litterfall%coarse_root%carbon14%ptr, veg%litterfall%coarse_root%nitrogen15%ptr , &
        ! sap_wood
        veg%litterfall%sap_wood%carbon%ptr, veg%litterfall%sap_wood%nitrogen%ptr, veg%litterfall%sap_wood%phosphorus%ptr    , &
        veg%litterfall%sap_wood%carbon13%ptr, veg%litterfall%sap_wood%carbon14%ptr, veg%litterfall%sap_wood%nitrogen15%ptr  , &
        ! heart_wood
        veg%litterfall%heart_wood%carbon%ptr, veg%litterfall%heart_wood%nitrogen%ptr, veg%litterfall%heart_wood%phosphorus%ptr    , &
        veg%litterfall%heart_wood%carbon13%ptr, veg%litterfall%heart_wood%carbon14%ptr, veg%litterfall%heart_wood%nitrogen15%ptr  , &
        ! labile
        veg%litterfall%labile%carbon%ptr, veg%litterfall%labile%nitrogen%ptr, veg%litterfall%labile%phosphorus%ptr    , &
        veg%litterfall%labile%carbon13%ptr, veg%litterfall%labile%carbon14%ptr, veg%litterfall%labile%nitrogen15%ptr  , &
        ! reserve
        veg%litterfall%reserve%carbon%ptr, veg%litterfall%reserve%nitrogen%ptr, veg%litterfall%reserve%phosphorus%ptr   , &
        veg%litterfall%reserve%carbon13%ptr, veg%litterfall%reserve%carbon14%ptr, veg%litterfall%reserve%nitrogen15%ptr , &
        ! fruit
        veg%litterfall%fruit%carbon%ptr, veg%litterfall%fruit%nitrogen%ptr, veg%litterfall%fruit%phosphorus%ptr   , &
        veg%litterfall%fruit%carbon13%ptr, veg%litterfall%fruit%carbon14%ptr, veg%litterfall%fruit%nitrogen15%ptr , &
        ! seed_bed
        veg%litterfall%seed_bed%carbon%ptr, veg%litterfall%seed_bed%nitrogen%ptr, veg%litterfall%seed_bed%phosphorus%ptr  , &
        veg%litterfall%seed_bed%carbon13%ptr, veg%litterfall%seed_bed%carbon14%ptr, veg%litterfall%seed_bed%nitrogen15%ptr


    ! if this flag is set TRUE, check whether the current time step is a new year, if so 
    ! set stat = -1 to indicate "end of file" which invokes the re-reading of the climate forcing input file at line 1
    IF(IS_IOSTAT_END(stat)) THEN 
      ! start reading from line one & read header into dummy
      REWIND(lun)
      CALL read_forcing_header(lun)
      READ(lun,*) time_hlp, &
        ! leaf
        veg%litterfall%leaf%carbon%ptr, veg%litterfall%leaf%nitrogen%ptr, veg%litterfall%leaf%phosphorus%ptr    , &
        veg%litterfall%leaf%carbon13%ptr, veg%litterfall%leaf%carbon14%ptr, veg%litterfall%leaf%nitrogen15%ptr  , &
        ! fine_root
        veg%litterfall%fine_root%carbon%ptr, veg%litterfall%fine_root%nitrogen%ptr, veg%litterfall%fine_root%phosphorus%ptr   , &
        veg%litterfall%fine_root%carbon13%ptr, veg%litterfall%fine_root%carbon14%ptr, veg%litterfall%fine_root%nitrogen15%ptr , &
        ! coarse_root
        veg%litterfall%coarse_root%carbon%ptr, veg%litterfall%coarse_root%nitrogen%ptr, veg%litterfall%coarse_root%phosphorus%ptr   , &
        veg%litterfall%coarse_root%carbon13%ptr, veg%litterfall%coarse_root%carbon14%ptr, veg%litterfall%coarse_root%nitrogen15%ptr , &
        ! sap_wood
        veg%litterfall%sap_wood%carbon%ptr, veg%litterfall%sap_wood%nitrogen%ptr, veg%litterfall%sap_wood%phosphorus%ptr    , &
        veg%litterfall%sap_wood%carbon13%ptr, veg%litterfall%sap_wood%carbon14%ptr, veg%litterfall%sap_wood%nitrogen15%ptr  , &
        ! heart_wood
        veg%litterfall%heart_wood%carbon%ptr, veg%litterfall%heart_wood%nitrogen%ptr, veg%litterfall%heart_wood%phosphorus%ptr    , &
        veg%litterfall%heart_wood%carbon13%ptr, veg%litterfall%heart_wood%carbon14%ptr, veg%litterfall%heart_wood%nitrogen15%ptr  , &
        ! labile
        veg%litterfall%labile%carbon%ptr, veg%litterfall%labile%nitrogen%ptr, veg%litterfall%labile%phosphorus%ptr    , &
        veg%litterfall%labile%carbon13%ptr, veg%litterfall%labile%carbon14%ptr, veg%litterfall%labile%nitrogen15%ptr  , &
        ! reserve
        veg%litterfall%reserve%carbon%ptr, veg%litterfall%reserve%nitrogen%ptr, veg%litterfall%reserve%phosphorus%ptr   , &
        veg%litterfall%reserve%carbon13%ptr, veg%litterfall%reserve%carbon14%ptr, veg%litterfall%reserve%nitrogen15%ptr , &
        ! fruit
        veg%litterfall%fruit%carbon%ptr, veg%litterfall%fruit%nitrogen%ptr, veg%litterfall%fruit%phosphorus%ptr   , &
        veg%litterfall%fruit%carbon13%ptr, veg%litterfall%fruit%carbon14%ptr, veg%litterfall%fruit%nitrogen15%ptr , &
        ! seed_bed
        veg%litterfall%seed_bed%carbon%ptr, veg%litterfall%seed_bed%nitrogen%ptr, veg%litterfall%seed_bed%phosphorus%ptr  , &
        veg%litterfall%seed_bed%carbon13%ptr, veg%litterfall%seed_bed%carbon14%ptr, veg%litterfall%seed_bed%nitrogen15%ptr
    END IF

    lun=10
    READ(lun,*, IOSTAT=stat) time_hlp, & 
        soil%t_soil%ptr(ics:ice,1:nlevel_soil,iblk), tile%vegsoil%moisture   ! *slmdev t_soil uses t_air

    IF(IS_IOSTAT_END(stat)) THEN 
      ! start reading from line one & read header into dummy
      REWIND(lun)
      CALL read_forcing_header(lun)
      READ(lun,*) time_hlp, &
        soil%t_soil%ptr(ics:ice,1:nlevel_soil,iblk), tile%vegsoil%moisture   ! *slmdev t_soil uses t_air
    END IF

    END SELECT  ! land plant canopy vs soil


  END SUBROUTINE update_forcing


  !-----------------------------------------------------------------------------------------------------
  !> open the file test_bed/plant/climate.dat, only IF model option is set to READ
  !! 
  !! called in the driver
  !! 
  !! note: lun is set to 9 (see ticket 31)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE open_forcing_input_file

    USE mo_exception,                       ONLY: finish
    USE mo_base_config_class,               ONLY: base_config

    IMPLICIT NONE
    ! 0.2 Local
    LOGICAL :: file_exists
    INTEGER :: lun 
    CHARACTER(100) :: fname
    CHARACTER(len=*), PARAMETER :: routine = 'mo_atmland_forcing:open_forcing_input_file'


    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant','canopy')
        !> 1.1 open file for READ
        !! 
        lun=9
        fname = 'climate.dat'
        ! check if file exists
        INQUIRE(FILE=fname, EXIST=file_exists)
        ! open if fname == True
        IF(file_exists) THEN
          OPEN (UNIT=lun, FILE=fname, ACTION="read", STATUS="old") ! status=old: file exists already
        !> 1.2 read the two header lines into dummy variables
        !! to advance to line 3
          CALL read_forcing_header(lun)
        ELSE
          CALL finish(TRIM(routine), 'climate forcing: climate.dat could not be opened, file does not exist')
        END IF

      CASE('soil')
        !> 2.1 open file for READ
        !! 
        lun=9
        fname = 'soil_forcing_litterfall.dat'
        ! check if file exists
        INQUIRE(FILE=fname, EXIST=file_exists)
        ! open if fname == True
        IF(file_exists) THEN
          OPEN (UNIT=lun, FILE=fname, ACTION="read", STATUS="old") ! status=old: file exists already
        !> 2.2 read the header line into dummy variable
        !! to advance to line 2
          CALL read_forcing_header(lun)
        ELSE
          CALL finish(TRIM(routine), 'soil forcing: soil_forcing_litterfall.dat could not be opened, file does not exist')
        END IF

        !> 3.1 open file for READ
        !! 
        lun=10
        fname = 'soil_forcing_soilphysics.dat'
        ! check if file exists
        INQUIRE(FILE=fname, EXIST=file_exists)
        ! open if fname == True
        IF(file_exists) THEN
          OPEN (UNIT=lun, FILE=fname, ACTION="read", STATUS="old") ! status=old: file exists already
        !> 3.2 read the header line into dummy variable
        !! to advance to line 2
          CALL read_forcing_header(lun)
        ELSE
          CALL finish(TRIM(routine), 'soil forcing: soil_forcing_soilphysics.dat could not be opened, file does not exist')
        END IF

      CASE DEFAULT
        CALL finish(TRIM(routine), 'climate forcing: either land, plant, canopy or soil should lead you here, but others not')
    END SELECT

  END SUBROUTINE open_forcing_input_file


  !-----------------------------------------------------------------------------------------------------
  !> closes the file test_bed/plant/climate.dat, only IF model option is set to READ
  !! 
  !! called in the driver
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE close_forcing_input_file

    USE mo_base_config_class,               ONLY: base_config

    IMPLICIT NONE

    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant','canopy')
        CLOSE(9)    ! climate

      CASE('soil')
        CLOSE(9)    ! litterfall
        CLOSE(10)   ! soilphysics
    END SELECT

  END SUBROUTINE close_forcing_input_file


  !-----------------------------------------------------------------------------------------------------
  !> read only the header into a dummy variable to advanve the READ to the 3rd line in the file
  !! 
  !! consider this a helper subroutine to the read_forcing subroutine
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE read_forcing_header(lun)

    USE mo_base_config_class,               ONLY: base_config

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER,INTENT(in)  :: lun 
    ! 0.2 Local
    CHARACTER(200)      :: header_dummy1 = "header"
    CHARACTER(500)      :: header_dummy2 = "header_long"


    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant','canopy')
      !> 1.0 read 2 header lines of one file
      !! lun = 9
        READ(lun,'(a)') header_dummy1   ! header line 1
        READ(lun,'(a)') header_dummy1   ! header line 2

      CASE('soil')
      !> 2.0 read 1 header line of two files
      !! lun = 9 or 10
        READ(lun,'(a)') header_dummy2   ! header line 1
    END SELECT

  END SUBROUTINE read_forcing_header

END MODULE mo_atmland_forcing
