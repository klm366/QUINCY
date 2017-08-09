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
  USE mo_shared_isotope,        ONLY: calc_mixing_ratio_C13C12, calc_mixing_ratio_C14C, calc_delta_C13, calc_delta_C14, calc_mixing_ratio_N15N14
  USE mo_atmland_util,          ONLY: update_orbital_parameters, calc_solar_angle, calc_sw_fdiffuse, calc_spec_humidity_sat
  USE mo_atmland_constants


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_forcing, open_forcing_input_file, close_forcing_input_file

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_forcing'

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
                                  iblk               , &  ! in
                                  ics                , &  ! in
                                  ice                , &  ! in
                                  nc                 , &  ! in
                                  nsoil_sb           , &  ! in
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
                                  nhx_deposition     , &  ! out
                                  noy_deposition     , &  ! out
                                  nhx_n15_deposition , &  ! out
                                  noy_n15_deposition , &  ! out
                                  p_deposition       , &  ! out
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
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_base_memory,           ONLY: t_slm_tile, grid  ! grid is used for calc_solar_angle()
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_base_config_class,     ONLY: base_config
    USE mo_shared_constants,      ONLY: molar_mass_N, molar_mass_P, grav
    USE mo_base_constants,        ONLY: one_day
    USE mo_jsb4_forcing,          ONLY: forcing_options

    ! Use of process memories
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile)  , INTENT(inout) :: tile     !< one tile
    INTEGER           , INTENT(in)    :: iblk, ics, ice, nc, nsoil_sb
    REAL(wp)          , INTENT(in)    :: time     !< current time (not timestep !)
    REAL(wp), OPTIONAL, INTENT(out)   :: &        !! these are all variables from a2l memory
                                      t_air(:)              , &    !< ..
                                      q_air(:)              , &    !< ..
                                      rain(:)               , &    !< ..
                                      snow(:)               , &    !< ..
                                      wind_air(:)           , &    !< ..
                                      wind_10m(:)           , &    !< ..
                                      lw_srf_down(:)        , &    !< ..
                                      swvis_srf_down(:)     , &    !< ..
                                      swnir_srf_down(:)     , &    !< ..
                                      swpar_srf_down(:)     , &    !< ..
                                      swpar_fdiffuse(:)     , &    !< ..
                                      press_srf(:)          , &    !< ..
                                      drag_srf(:)           , &    !< ..
                                      t_acoef(:)            , &    !< ..
                                      t_bcoef(:)            , &    !< ..
                                      q_acoef(:)            , &    !< ..
                                      q_bcoef(:)            , &    !< ..
                                      pch(:)                , &    !< ..
                                      cos_zenith_angle(:)   , &    !< ..
                                      co2_mixing_ratio(:)     , &  !< ..
                                      co2_mixing_ratio_C13(:) , &  !< ..
                                      co2_mixing_ratio_C14(:) , &  !< ..
                                      ga(:)                   , &  !< ..
                                      nhx_deposition(:)     , &    !< ..
                                      noy_deposition(:)     , &    !< ..
                                      nhx_n15_deposition(:) , &    !< ..
                                      noy_n15_deposition(:) , &    !< ..
                                      p_deposition(:)       , &    !< ..
                                      ! For lakes:
                                      albedo_lwtr(:)        , &    !< ..
                                      drag_wtr(:)           , &    !< ..
                                      drag_ice(:)           , &    !< ..
                                      t_acoef_wtr(:)        , &    !< ..
                                      t_bcoef_wtr(:)        , &    !< ..
                                      q_acoef_wtr(:)        , &    !< ..
                                      q_bcoef_wtr(:)        , &    !< ..
                                      t_acoef_ice(:)        , &    !< ..
                                      t_bcoef_ice(:)        , &    !< ..
                                      q_acoef_ice(:)        , &    !< ..
                                      q_bcoef_ice(:)               !< ..
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER     :: model                  ! the model structure
    ! READ() input variables for: land plant canopy (this variable names can differ from column names in forcing files)
    REAL(wp) ::      swvis_srf_down_r,   lw_srf_down_r, t_air_r,    q_air_r,    &
                     press_srf_r,        rain_r,        snow_r,     wind_air_r, &
                     co2_mixing_ratio_r, co2_dc13_r,    co2_dc14_r,             &
                     nhx_deposition_r, noy_deposition_r, &   ! variables reading nhx_srf_down & noy_srf_down from the climate forcing file (maybe rename there to)
                     nhx_dn15_deposition_r, noy_dn15_deposition_r, & 
                     p_deposition_r 
    ! READ() input variables for: soil
    REAL(wp) ::      leaf_carbon_r, leaf_nitrogen_r, leaf_phosphorus_r    , &
                     leaf_dC13_r, leaf_DC14_r, leaf_dN15_r  , &
                     fine_root_carbon_r, fine_root_nitrogen_r, fine_root_phosphorus_r   , &
                     fine_root_dC13_r, fine_root_DC14_r, fine_root_dN15_r , &
                     coarse_root_carbon_r, coarse_root_nitrogen_r, coarse_root_phosphorus_r   , &
                     coarse_root_dC13_r, coarse_root_DC14_r, coarse_root_dN15_r , &
                     sap_wood_carbon_r, sap_wood_nitrogen_r, sap_wood_phosphorus_r    , &
                     sap_wood_dC13_r, sap_wood_DC14_r, sap_wood_dN15_r  , &
                     heart_wood_carbon_r, heart_wood_nitrogen_r, heart_wood_phosphorus_r    , &
                     heart_wood_dC13_r, heart_wood_DC14_r, heart_wood_dN15_r  , &
                     labile_carbon_r, labile_nitrogen_r, labile_phosphorus_r    , &
                     labile_dC13_r, labile_DC14_r, labile_dN15_r  , &
                     reserve_carbon_r, reserve_nitrogen_r, reserve_phosphorus_r   , &
                     reserve_dC13_r, reserve_DC14_r, reserve_dN15_r , &
                     fruit_carbon_r, fruit_nitrogen_r, fruit_phosphorus_r   , &
                     fruit_dC13_r, fruit_DC14_r, fruit_dN15_r , &
                     seed_bed_carbon_r, seed_bed_nitrogen_r, seed_bed_phosphorus_r  , &
                     seed_bed_dC13_r, seed_bed_DC14_r, seed_bed_dN15_r
    REAL(wp), DIMENSION(nsoil_sb) :: & 
                     t_soil_sl_r, w_soil_sl_r, &
                     percolation_sl_r, root_fraction_sl_r
    REAL(wp) ::      veg_pool_fine_root_carbon_r


    INTEGER                   :: stat          ! status of I/O (https://software.intel.com/en-us/node/691976)
    INTEGER                   :: lun           ! logical unit number for file input stream, see '2.0 read data'
    REAL(wp), DIMENSION(nc)   :: sw_toa_down   ! potential radiation at the top of the atmosphere
    REAL(wp), DIMENSION(nc)   :: hlp1          ! dummy variable 
    REAL(wp)                  :: year          ! dummy variable
    REAL(wp)                  :: hour          ! dummy variable
    REAL(wp)                  :: time_hlp      ! dummy variable for soil forcing
    REAL(wp), SAVE            :: day                    ! read input dummy variable - save value
    REAL(wp)                  :: doy_old, &             ! check for new year (column 2 in input file would change from 365 -> 1)
                                 doy_new
    INTEGER                   :: nlevel_soil  = -1      ! number of levels at the Z-Axis (soil layer)
    INTEGER                   :: nlevel_veg   = -1      ! number of levels at the Z-Axis (canopy layer)
    INTEGER                   :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_forcing'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(SOIL_)                       !     TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    ! Declare pointers to variables in memory
    ! SOIL_
    dsl4jsb_Real3D_onChunk      :: t_soil_sl
    dsl4jsb_Real3D_onChunk      :: w_soil_sl
    dsl4jsb_Real3D_onChunk      :: percolation_sl
    dsl4jsb_Real3D_onChunk      :: root_fraction_sl
! 2 be removed: pointers to veg assimi *slmtemp
    TYPE(t_veg_memory),      POINTER :: veg         ! local pointer to veg__mem for SLM usability
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
!    lctlib => model%lctlib(tile%lcts%lib_id) 
    ! Get process config
!    dsl4jsb_Get_config(xx_)
    ! Get process memories
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(VEG_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var3D_onChunk(SOIL_, t_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, w_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, percolation_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, root_fraction_sl)


! 2 be removed  *slmtemp
    veg     => veg__mem 


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    doy_old = 0.0_wp  ! doy in the previous timestep, used for detecting a new year

    ! give a2l variables a value to avoid compiler warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value
    t_air(nc)           = 0.0_wp
    q_air(nc)           = 0.0_wp
    rain(nc)            = 0.0_wp
    snow(nc)            = 0.0_wp
    wind_air(nc)        = 0.0_wp
    wind_10m(nc)        = 0.0_wp
    lw_srf_down(nc)     = 0.0_wp
    swvis_srf_down(nc)  = 0.0_wp
    swnir_srf_down(nc)  = 0.0_wp
    swpar_srf_down(nc)  = 0.0_wp
    swpar_fdiffuse(nc)  = 0.0_wp
    press_srf(nc)       = 0.0_wp
    drag_srf(nc)        = 1000.0_wp    ! 1.0_wp
    t_acoef(nc)         = 0.0_wp
    t_bcoef(nc)         = 0.0_wp
    q_acoef(nc)         = 0.0_wp
    q_bcoef(nc)         = 0.0_wp
    pch(nc)             = 0.0_wp
    cos_zenith_angle(nc)      = 0.0_wp
    co2_mixing_ratio(nc)      = 0.0_wp
    co2_mixing_ratio_C13(nc)  = 0.0_wp
    co2_mixing_ratio_C14(nc)  = 0.0_wp
    ga(nc)                    = 0.0_wp
    nhx_deposition(nc)        = 0.0_wp
    noy_deposition(nc)        = 0.0_wp
    nhx_n15_deposition(nc)        = 0.0_wp
    noy_n15_deposition(nc)        = 0.0_wp
    p_deposition(nc)    = 0.0_wp
    ! For lakes:
    albedo_lwtr(nc)     = 0.0_wp
    drag_wtr(nc)        = 0.0_wp
    drag_ice(nc)        = 0.0_wp
    t_acoef_wtr(nc)     = 0.0_wp
    t_bcoef_wtr(nc)     = 0.0_wp
    q_acoef_wtr(nc)     = 0.0_wp
    q_bcoef_wtr(nc)     = 0.0_wp
    t_acoef_ice(nc)     = 0.0_wp
    t_bcoef_ice(nc)     = 0.0_wp
    q_acoef_ice(nc)     = 0.0_wp
    q_bcoef_ice(nc)     = 0.0_wp


    !!#######################
    !! LAND PLANT CANOPY
    !!#######################
    SELECT CASE(TRIM(base_config%slm_model_name))
    CASE('land','plant','canopy')

    !> 1.0 orbital parameters
    !!
    !! 'declination' is calculated as output paramater, but not needed here
    CALL update_orbital_parameters(time,declination,sw_toa_down(1))  ! mo_atmland_util


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
                     swvis_srf_down_r,   lw_srf_down_r, t_air_r,    q_air_r,    &
                     press_srf_r,        rain_r,        snow_r,     wind_air_r, &
                     co2_mixing_ratio_r, co2_dc13_r,    co2_dc14_r,             &
                     nhx_deposition_r, noy_deposition_r, p_deposition_r 

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
                     swvis_srf_down_r,   lw_srf_down_r, t_air_r,    q_air_r,    &
                     press_srf_r,        rain_r,        snow_r,     wind_air_r, &
                     co2_mixing_ratio_r, co2_dc13_r,    co2_dc14_r,             &
                     nhx_deposition_r, noy_deposition_r, p_deposition_r 
    END IF

    !> 3.0 parameter units are re-scaled
    !!
    q_air_r     = q_air_r / 1000._wp
    press_srf_r = press_srf_r * 100._wp
    rain_r      = rain_r / 86400._wp
    snow_r      = snow_r / 86400._wp

    !> 4.0 update forcing variables based on 'read data'
    !!
    wind_10m(1)             = wind_air_r
    ! neutral drag coefficient, assuminging a 12 metre tall forest
    ! In JSBACH 4, this is divided by stepen and grav, so these are multiplied here... 
    drag_srf(1)             = (( 0.41 / log( (22._wp - 2.0_wp/3._wp * 12._wp)/1.0_wp) )**2)
    ga(1)                   = drag_srf(1) * wind_10m(1) 
    drag_srf(1)             = drag_srf(1) * 1800._wp * grav 
    IF (forcing_options%flag_read_dC13) THEN
      hlp1(1)               = calc_mixing_ratio_C13C12(co2_dc13_r)
    ELSE
      hlp1(1)               = calc_mixing_ratio_C13C12(def_co2_deltaC13)
    END IF
    co2_mixing_ratio_C13(1) = co2_mixing_ratio_r / (1._wp + 1._wp/hlp1(1))
    IF (forcing_options%flag_read_DC14) THEN
      hlp1(1)               = calc_mixing_ratio_C14C(co2_dc13_r, co2_dc14_r)
    ELSE
      hlp1(1)               = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14)
    END IF
    co2_mixing_ratio_C14(1) = hlp1(1) * co2_mixing_ratio_r
    cos_zenith_angle(1) = calc_solar_angle(grid%lat,grid%lon,time)                      
    swnir_srf_down(1)   = (1._wp - frac_vis_swtot_srf ) * swvis_srf_down_r 
    swvis_srf_down_r    = frac_vis_swtot_srf * swvis_srf_down_r 
    swpar_srf_down(1)   = frac_par_swvis_srf * swvis_srf_down_r
    swpar_fdiffuse(1)   = calc_sw_fdiffuse(swvis_srf_down_r + swnir_srf_down(1), sw_toa_down(1), cos_zenith_angle(1)) ! function: calc_sw_fdiffuse

    !> 5.0 pass values read from file (i.e. scalar values) to ics-ice
    !!
    t_air(nc)           =         t_air_r
    q_air(nc)           =         q_air_r
    rain(nc)            =         rain_r
    snow(nc)            =         snow_r
    wind_air(nc)        =         wind_air_r
    wind_10m(nc)        =         wind_10m(1)       
    lw_srf_down(nc)     =         lw_srf_down_r
    swvis_srf_down(nc)  =         swvis_srf_down_r
    swnir_srf_down(nc)  =         swnir_srf_down(1) 
    swpar_srf_down(nc)  =         swpar_srf_down(1) 
    swpar_fdiffuse(nc)  =         swpar_fdiffuse(1) 
    press_srf(nc)       =         press_srf_r
    drag_srf(nc)        =         drag_srf(1)       
    t_acoef(nc)         =         t_acoef(1)        
    t_bcoef(nc)         =         t_bcoef(1)        
    q_acoef(nc)         =         q_acoef(1)        
    q_bcoef(nc)         =         q_bcoef(1)        
    pch(nc)             =         pch(1)            
    cos_zenith_angle(nc)      =   cos_zenith_angle(1)     
    co2_mixing_ratio(nc)      =   co2_mixing_ratio_r
    co2_mixing_ratio_C13(nc)  =   co2_mixing_ratio_C13(1) 
    co2_mixing_ratio_C14(nc)  =   co2_mixing_ratio_C14(1) 
    ga(nc)                    =   ga(1)         
    SELECT CASE (TRIM(forcing_options%n_deposition_scheme))
      CASE ("none")
        nhx_deposition(nc)    =   0.0_wp
        noy_deposition(nc)    =   0.0_wp
      CASE ("constant")
        nhx_deposition(nc)    =   def_n_deposition / molar_mass_N * 1000.0_wp / one_day    ! mg/m2/day -> mumol/m2/sec
        noy_deposition(nc)    =   def_n_deposition / molar_mass_N * 1000.0_wp / one_day    ! mg/m2/day -> mumol/m2/sec
      CASE ("dynamic")
        nhx_deposition(nc)    =   nhx_deposition_r / molar_mass_N * 1000.0_wp / one_day    ! mg/m2/day -> mumol/m2/sec
        noy_deposition(nc)    =   noy_deposition_r / molar_mass_N * 1000.0_wp / one_day    ! mg/m2/day -> mumol/m2/sec
    ENDSELECT
    nhx_n15_deposition(nc)    =   nhx_deposition(nc) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
    noy_n15_deposition(nc)    =   noy_deposition(nc) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 

    SELECT CASE (TRIM(forcing_options%p_deposition_scheme))
      CASE ("none")
        p_deposition(nc)      =   0.0_wp
      CASE ("constant")
        p_deposition(nc)      =   def_p_deposition / molar_mass_P * 1000.0_wp / one_day    ! mg/m2/day -> mumol/m2/sec
      CASE ("dynamic")
        p_deposition(nc)      =   p_deposition_r / molar_mass_P * 1000.0_wp / one_day      ! mg/m2/day -> mumol/m2/sec
    ENDSELECT
    ! For lakes:    
    albedo_lwtr(nc)     =         albedo_lwtr(1)    
    drag_wtr(nc)        =         drag_wtr(1)       
    drag_ice(nc)        =         drag_ice(1)       
    t_acoef_wtr(nc)     =         t_acoef_wtr(1)    
    t_bcoef_wtr(nc)     =         t_bcoef_wtr(1)    
    q_acoef_wtr(nc)     =         q_acoef_wtr(1)    
    q_bcoef_wtr(nc)     =         q_bcoef_wtr(1)    
    t_acoef_ice(nc)     =         t_acoef_ice(1)    
    t_bcoef_ice(nc)     =         t_bcoef_ice(1)    
    q_acoef_ice(nc)     =         q_acoef_ice(1)    
    q_bcoef_ice(nc)     =         q_bcoef_ice(1)    



    !!#######################
    !! SOIL
    !!#######################
    CASE('soil')

    nlevel_soil = nsoil_sb

    !> 1.0 read data
    !!
    !! the read statement, by default a list of values in a line is read and then advanced to the beginning of the next line
    lun=9
    READ(lun,*, IOSTAT=stat) time_hlp, & 
                     leaf_carbon_r, leaf_nitrogen_r, leaf_phosphorus_r    , &
                     leaf_dC13_r, leaf_DC14_r, leaf_dN15_r  , &
                     fine_root_carbon_r, fine_root_nitrogen_r, fine_root_phosphorus_r   , &
                     fine_root_dC13_r, fine_root_DC14_r, fine_root_dN15_r , &
                     coarse_root_carbon_r, coarse_root_nitrogen_r, coarse_root_phosphorus_r   , &
                     coarse_root_dC13_r, coarse_root_DC14_r, coarse_root_dN15_r , &
                     sap_wood_carbon_r, sap_wood_nitrogen_r, sap_wood_phosphorus_r    , &
                     sap_wood_dC13_r, sap_wood_DC14_r, sap_wood_dN15_r  , &
                     heart_wood_carbon_r, heart_wood_nitrogen_r, heart_wood_phosphorus_r    , &
                     heart_wood_dC13_r, heart_wood_DC14_r, heart_wood_dN15_r  , &
                     labile_carbon_r, labile_nitrogen_r, labile_phosphorus_r    , &
                     labile_dC13_r, labile_DC14_r, labile_dN15_r  , &
                     reserve_carbon_r, reserve_nitrogen_r, reserve_phosphorus_r   , &
                     reserve_dC13_r, reserve_DC14_r, reserve_dN15_r , &
                     fruit_carbon_r, fruit_nitrogen_r, fruit_phosphorus_r   , &
                     fruit_dC13_r, fruit_DC14_r, fruit_dN15_r , &
                     seed_bed_carbon_r, seed_bed_nitrogen_r, seed_bed_phosphorus_r  , &
                     seed_bed_dC13_r, seed_bed_DC14_r, seed_bed_dN15_r


    ! if this flag is set TRUE, check whether the current time step is a new year, if so 
    ! set stat = -1 to indicate "end of file" which invokes the re-reading of the climate forcing input file at line 1
    IF(IS_IOSTAT_END(stat)) THEN 
      ! start reading from line one & read header into dummy
      REWIND(lun)
      CALL read_forcing_header(lun)
      READ(lun,*) time_hlp, &
                     leaf_carbon_r, leaf_nitrogen_r, leaf_phosphorus_r    , &
                     leaf_dC13_r, leaf_DC14_r, leaf_dN15_r  , &
                     fine_root_carbon_r, fine_root_nitrogen_r, fine_root_phosphorus_r   , &
                     fine_root_dC13_r, fine_root_DC14_r, fine_root_dN15_r , &
                     coarse_root_carbon_r, coarse_root_nitrogen_r, coarse_root_phosphorus_r   , &
                     coarse_root_dC13_r, coarse_root_DC14_r, coarse_root_dN15_r , &
                     sap_wood_carbon_r, sap_wood_nitrogen_r, sap_wood_phosphorus_r    , &
                     sap_wood_dC13_r, sap_wood_DC14_r, sap_wood_dN15_r  , &
                     heart_wood_carbon_r, heart_wood_nitrogen_r, heart_wood_phosphorus_r    , &
                     heart_wood_dC13_r, heart_wood_DC14_r, heart_wood_dN15_r  , &
                     labile_carbon_r, labile_nitrogen_r, labile_phosphorus_r    , &
                     labile_dC13_r, labile_DC14_r, labile_dN15_r  , &
                     reserve_carbon_r, reserve_nitrogen_r, reserve_phosphorus_r   , &
                     reserve_dC13_r, reserve_DC14_r, reserve_dN15_r , &
                     fruit_carbon_r, fruit_nitrogen_r, fruit_phosphorus_r   , &
                     fruit_dC13_r, fruit_DC14_r, fruit_dN15_r , &
                     seed_bed_carbon_r, seed_bed_nitrogen_r, seed_bed_phosphorus_r  , &
                     seed_bed_dC13_r, seed_bed_DC14_r, seed_bed_dN15_r
    END IF

    lun=10
    READ(lun,*, IOSTAT=stat) time_hlp, & 
                     t_soil_sl_r(1:nlevel_soil), w_soil_sl_r(1:nlevel_soil), &
                     percolation_sl_r(1:nlevel_soil), root_fraction_sl_r(1:nlevel_soil), &
                     veg_pool_fine_root_carbon_r, &
                     nhx_deposition_r, noy_deposition_r, &
                     nhx_dn15_deposition_r, noy_dn15_deposition_r

    IF(IS_IOSTAT_END(stat)) THEN 
      ! start reading from line one & read header into dummy
      REWIND(lun)
      CALL read_forcing_header(lun)
      READ(lun,*) time_hlp, &
                     t_soil_sl_r(1:nlevel_soil), w_soil_sl_r(1:nlevel_soil), &
                     percolation_sl_r(1:nlevel_soil), root_fraction_sl_r(1:nlevel_soil), &
                     veg_pool_fine_root_carbon_r, &
                     nhx_deposition_r, noy_deposition_r, &
                     nhx_dn15_deposition_r, noy_dn15_deposition_r
    END IF

    !> 2.0 pass values read from file (i.e. scalar values) to ics-ice
    !!
    veg%litterfall%leaf%carbon%ptr     = leaf_carbon_r / 1000._wp
    veg%litterfall%leaf%nitrogen%ptr   = leaf_nitrogen_r / 1000000._wp
    veg%litterfall%leaf%phosphorus%ptr = leaf_phosphorus_r / 1000000._wp   
    veg%litterfall%leaf%carbon13%ptr   = veg%litterfall%leaf%carbon%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_C13C12(leaf_dC13_r))  
    veg%litterfall%leaf%carbon14%ptr   = veg%litterfall%leaf%carbon%ptr * calc_mixing_ratio_C14C(leaf_dC13_r,leaf_DC14_r) 
    veg%litterfall%leaf%nitrogen15%ptr = veg%litterfall%leaf%nitrogen%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(leaf_dN15_r))

    veg%litterfall%fine_root%carbon%ptr     = fine_root_carbon_r / 1000._wp
    veg%litterfall%fine_root%nitrogen%ptr   = fine_root_nitrogen_r / 1000000._wp
    veg%litterfall%fine_root%phosphorus%ptr = fine_root_phosphorus_r / 1000000._wp
    veg%litterfall%fine_root%carbon13%ptr   = veg%litterfall%fine_root%carbon%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_C13C12(fine_root_dC13_r)) 
    veg%litterfall%fine_root%carbon14%ptr   = veg%litterfall%fine_root%carbon%ptr * calc_mixing_ratio_C14C(fine_root_dC13_r,fine_root_DC14_r) 
    veg%litterfall%fine_root%nitrogen15%ptr = veg%litterfall%fine_root%nitrogen%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(fine_root_dN15_r)) 

    veg%litterfall%coarse_root%carbon%ptr     = coarse_root_carbon_r / 1000._wp
    veg%litterfall%coarse_root%nitrogen%ptr   = coarse_root_nitrogen_r / 1000000._wp
    veg%litterfall%coarse_root%phosphorus%ptr = coarse_root_phosphorus_r / 1000000._wp
    veg%litterfall%coarse_root%carbon13%ptr   = veg%litterfall%coarse_root%carbon%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_C13C12(coarse_root_dC13_r)) 
    veg%litterfall%coarse_root%carbon14%ptr   = veg%litterfall%coarse_root%carbon%ptr * calc_mixing_ratio_C14C(coarse_root_dC13_r,coarse_root_DC14_r) 
    veg%litterfall%coarse_root%nitrogen15%ptr = veg%litterfall%coarse_root%nitrogen%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(coarse_root_dN15_r))

    veg%litterfall%sap_wood%carbon%ptr     = sap_wood_carbon_r / 1000._wp
    veg%litterfall%sap_wood%nitrogen%ptr   = sap_wood_nitrogen_r / 1000000._wp
    veg%litterfall%sap_wood%phosphorus%ptr = sap_wood_phosphorus_r / 1000000._wp
    veg%litterfall%sap_wood%carbon13%ptr   = veg%litterfall%sap_wood%carbon%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_C13C12(sap_wood_dC13_r))
    veg%litterfall%sap_wood%carbon14%ptr   = veg%litterfall%sap_wood%carbon%ptr * calc_mixing_ratio_C14C(sap_wood_dC13_r,sap_wood_DC14_r) 
    veg%litterfall%sap_wood%nitrogen15%ptr = veg%litterfall%sap_wood%nitrogen%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(sap_wood_dN15_r))

    veg%litterfall%heart_wood%carbon%ptr     = heart_wood_carbon_r / 1000._wp
    veg%litterfall%heart_wood%nitrogen%ptr   = heart_wood_nitrogen_r / 1000000._wp
    veg%litterfall%heart_wood%phosphorus%ptr = heart_wood_phosphorus_r / 1000000._wp 
    veg%litterfall%heart_wood%carbon13%ptr   = veg%litterfall%heart_wood%carbon%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_C13C12(heart_wood_dC13_r)) 
    veg%litterfall%heart_wood%carbon14%ptr   = veg%litterfall%heart_wood%carbon%ptr * calc_mixing_ratio_C14C(heart_wood_dC13_r,heart_wood_DC14_r) 
    veg%litterfall%heart_wood%nitrogen15%ptr = veg%litterfall%heart_wood%nitrogen%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(heart_wood_dN15_r))

    veg%litterfall%labile%carbon%ptr     = labile_carbon_r / 1000._wp
    veg%litterfall%labile%nitrogen%ptr   = labile_nitrogen_r / 1000000._wp
    veg%litterfall%labile%phosphorus%ptr = labile_phosphorus_r / 1000000._wp 
    veg%litterfall%labile%carbon13%ptr   = veg%litterfall%labile%carbon%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_C13C12(labile_dC13_r))
    veg%litterfall%labile%carbon14%ptr   = veg%litterfall%labile%carbon%ptr * calc_mixing_ratio_C14C(labile_dC13_r,labile_DC14_r) 
    veg%litterfall%labile%nitrogen15%ptr = veg%litterfall%labile%nitrogen%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(labile_dN15_r))

    veg%litterfall%reserve%carbon%ptr     = reserve_carbon_r / 1000._wp
    veg%litterfall%reserve%nitrogen%ptr   = reserve_nitrogen_r / 1000000._wp
    veg%litterfall%reserve%phosphorus%ptr = reserve_phosphorus_r / 1000000._wp
    veg%litterfall%reserve%carbon13%ptr   = veg%litterfall%reserve%carbon%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_C13C12(reserve_dC13_r)) 
    veg%litterfall%reserve%carbon14%ptr   = veg%litterfall%reserve%carbon%ptr * calc_mixing_ratio_C14C(reserve_dC13_r,reserve_DC14_r) 
    veg%litterfall%reserve%nitrogen15%ptr = veg%litterfall%reserve%nitrogen%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(reserve_dN15_r))

    veg%litterfall%fruit%carbon%ptr     = fruit_carbon_r / 1000._wp
    veg%litterfall%fruit%nitrogen%ptr   = fruit_nitrogen_r / 1000000._wp
    veg%litterfall%fruit%phosphorus%ptr = fruit_phosphorus_r / 1000000._wp
    veg%litterfall%fruit%carbon13%ptr   = veg%litterfall%fruit%carbon%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_C13C12(fruit_dC13_r))
    veg%litterfall%fruit%carbon14%ptr   = veg%litterfall%fruit%carbon%ptr * calc_mixing_ratio_C14C(fruit_dC13_r,fruit_DC14_r) 
    veg%litterfall%fruit%nitrogen15%ptr = veg%litterfall%fruit%nitrogen%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(fruit_dN15_r))

    veg%litterfall%seed_bed%carbon%ptr     = seed_bed_carbon_r / 1000._wp
    veg%litterfall%seed_bed%nitrogen%ptr   = seed_bed_nitrogen_r / 1000000._wp
    veg%litterfall%seed_bed%phosphorus%ptr = seed_bed_phosphorus_r / 1000000._wp
    veg%litterfall%seed_bed%carbon13%ptr   = veg%litterfall%seed_bed%carbon%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_C13C12(seed_bed_dC13_r)) 
    veg%litterfall%seed_bed%carbon14%ptr   = veg%litterfall%seed_bed%carbon%ptr * calc_mixing_ratio_C14C(seed_bed_dC13_r,seed_bed_DC14_r) 
    veg%litterfall%seed_bed%nitrogen15%ptr = veg%litterfall%seed_bed%nitrogen%ptr / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(seed_bed_dN15_r))

    DO isoil = 1,nlevel_soil

       t_soil_sl(:,isoil)                              = t_soil_sl_r(isoil)
       w_soil_sl(:,isoil)                              = w_soil_sl_r(isoil)
       percolation_sl(:,isoil)                         = percolation_sl_r(isoil)
       root_fraction_sl(:,isoil)                       = root_fraction_sl_r(isoil)
    ENDDO
    veg%pool%fine_root%carbon%ptr                      = veg_pool_fine_root_carbon_r

    nhx_deposition(:)    = nhx_deposition_r
    noy_deposition(:)    = noy_deposition_r
    nhx_n15_deposition(:) = nhx_deposition(:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(nhx_dn15_deposition_r))
    nhx_n15_deposition(:) = noy_deposition(:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(noy_dn15_deposition_r))

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
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':open_forcing_input_file'


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
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':close_forcing_input_file'

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
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':read_forcing_header'

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
