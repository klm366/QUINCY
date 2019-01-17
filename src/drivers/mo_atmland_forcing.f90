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

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message

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
  SUBROUTINE update_forcing(  tile               , &  ! inout
                              ics                , &  ! in
                              ice                , &  ! in
                              nc                 , &  ! in
                              iblk               , &  ! in
                              nsoil_sb           , &  ! in
                              timestep           , &  ! in
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
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_grid_class,        ONLY: t_jsb_grid
    USE mo_jsb_grid,              ONLY: Get_Grid
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_phys_constants,    ONLY: molar_mass_N, molar_mass_P, grav
    USE mo_jsb_math_constants,    ONLY: dtime, one_day, eps8, eps4
    USE mo_jsb4_forcing,          ONLY: forcing_options
    USE mo_isotope_util,          ONLY: calc_mixing_ratio_C13C12, calc_mixing_ratio_C14C, calc_delta_C13, calc_delta_C14, &
                                        calc_mixing_ratio_N15N14
    USE mo_atmland_util,          ONLY: update_orbital_parameters, calc_solar_angle, calc_sw_fdiffuse, calc_spec_humidity_sat
    USE mo_atmland_constants,     ONLY: def_co2_deltaC13, def_co2_deltaC14, def_n_deposition, def_p_deposition, &
                                        declination, frac_vis_swtot_srf, frac_par_swvis_srf

    ! Use of process memories
    dsl4jsb_Use_memory(SSE_)       !  USE mo_sse_memory_class,     ONLY: t_sse_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4slm_Use_memory_a_veg
    USE mo_veg_memory_class,      ONLY: veg_pool, veg_litterfall

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile                           !< one tile with data structure for one lct
    INTEGER                   , INTENT(in)    :: ics, ice, nc, iblk, nsoil_sb   !< info on block, chunk and soil layer number
    INTEGER                   , INTENT(in)    :: timestep                       !< timestep; first timestep of static simulation and simulation after transient spinup is 1
    REAL(wp)                  , INTENT(in)    :: time                           !< current time in days (not timestep !)
    REAL(wp), OPTIONAL        , INTENT(out)   :: &                 !! these are all variables from a2l memory; with assumed shape VAR(1:nc)
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
    TYPE(t_jsb_grid),  POINTER     :: hgrid                  ! horizontal grid
    ! READ() input variables for: land plant canopy and soil (these variable names can differ from column names in forcing files)
    REAL(wp) ::      swvis_srf_down_r,   lw_srf_down_r, t_air_r,    q_air_r,    &
                     press_srf_r,        rain_r,        snow_r,     wind_air_r, &
                     co2_mixing_ratio_r, co2_dC13_r,    co2_DC14_r,             &
                     nhx_deposition_r,   noy_deposition_r, &   ! variables reading nhx_srf_down & noy_srf_down from the climate forcing file (maybe rename there)
                     nhx_dn15_deposition_r, noy_dn15_deposition_r, & 
                     p_deposition_r, f_n_demand_r, f_p_demand_r 
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
    REAL(wp), DIMENSION(nsoil_sb) ::  & 
                     t_soil_sl_r,     &
                     w_soil_sl_r,     &
                     percolation_sl_r
    REAL(wp), DIMENSION(nsoil_sb) ::  & 
                     root_fraction_sl_r
    REAL(wp) ::      veg_pool_fine_root_carbon_r
    ! other local variables
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
    INTEGER                   :: nlevel_veg   = -1      ! number of levels at the Z-Axis (canopy layer)
    INTEGER                   :: isoil, ichunk
    LOGICAL                   :: l_pulse_exp_add_15N    ! 15N pulse experiment
    REAL(wp)                  :: pulse_exp_add          ! 15N pulse experiment
    REAL(wp)                  :: pulse_exp_isofrac      ! 15N pulse experiment
    REAL(wp)                  :: pulse_exp_delta        ! 15N pulse experiment
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_forcing'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(SSE_)                       !     TYPE(t_sse_memory),     POINTER :: sse__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    ! Declare pointers to variables in memory
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: f_n_demand
    dsl4jsb_Real2D_onChunk      :: f_p_demand
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk      :: root_fraction_sl
    ! SSE_
    dsl4jsb_Real3D_onChunk      :: t_soil_sl
    dsl4jsb_Real3D_onChunk      :: w_soil_sl
    dsl4jsb_Real3D_onChunk      :: percolation_sl
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
!    lctlib => model%lctlib(tile%lcts%lib_id) 
    hgrid  => Get_grid(model%grid_id)
    ! Get process config
!    dsl4jsb_Get_config(xx_)
    ! Get process memories
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(VEG_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(VEG_,  f_n_demand)
    dsl4jsb_Get_var2D_onChunk(VEG_,  f_p_demand)
    dsl4jsb_Get_var3D_onChunk(VEG_,  root_fraction_sl)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, w_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, percolation_sl)
    ! Get process memory pool structures
    dsl4slm_get_vegpool_var(pool,fine_root,carbon)
#define Get_veg_litterfall
#include "slm_pools_a_incl.f90"
#undef Get_veg_litterfall


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    doy_old = 0.0_wp  ! doy in the previous timestep, used for detecting a new year

    ! give a2l variables a value to avoid compiler warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value
    ! these variables are declared in land_model_driver as: REAL(wp), ALLOCATABLE, DIMENSION(:)
    ! and allocated using 'nc'   ALLOCATE(t_air(nc)) -> i.e. no ics:ice,iblk here !
    t_air(:)           = 0.0_wp
    q_air(:)           = 0.0_wp
    rain(:)            = 0.0_wp
    snow(:)            = 0.0_wp
    wind_air(:)        = 0.0_wp
    wind_10m(:)        = 0.0_wp
    lw_srf_down(:)     = 0.0_wp
    swvis_srf_down(:)  = 0.0_wp
    swnir_srf_down(:)  = 0.0_wp
    swpar_srf_down(:)  = 0.0_wp
    swpar_fdiffuse(:)  = 0.0_wp
    press_srf(:)       = 0.0_wp
    drag_srf(:)        = 1000.0_wp
    t_acoef(:)         = 0.0_wp
    t_bcoef(:)         = 0.0_wp
    q_acoef(:)         = 0.0_wp
    q_bcoef(:)         = 0.0_wp
    pch(:)             = 0.0_wp
    cos_zenith_angle(:)      = 0.0_wp
    co2_mixing_ratio(:)      = 0.0_wp
    co2_mixing_ratio_C13(:)  = 0.0_wp
    co2_mixing_ratio_C14(:)  = 0.0_wp
    ga(:)                    = 0.0_wp
    nhx_deposition(:)        = 0.0_wp
    noy_deposition(:)        = 0.0_wp
    nhx_n15_deposition(:)    = 0.0_wp
    noy_n15_deposition(:)    = 0.0_wp
    p_deposition(:)          = 0.0_wp
    ! For lakes:
    albedo_lwtr(:)     = 0.0_wp
    drag_wtr(:)        = 0.0_wp
    drag_ice(:)        = 0.0_wp
    t_acoef_wtr(:)     = 0.0_wp
    t_bcoef_wtr(:)     = 0.0_wp
    q_acoef_wtr(:)     = 0.0_wp
    q_bcoef_wtr(:)     = 0.0_wp
    t_acoef_ice(:)     = 0.0_wp
    t_bcoef_ice(:)     = 0.0_wp
    q_acoef_ice(:)     = 0.0_wp
    q_bcoef_ice(:)     = 0.0_wp


    !! #############################################################################################
    !! LAND PLANT CANOPY
    !! #############################################################################################
    SELECT CASE(TRIM(model%config%slm_model_name))
    CASE('land','plant','canopy')

      !> 1.0 orbital parameters
      !!
      !! "declination" is a parameter from mo_atmland_constants and is calculated here as output paramater, though, not used here in this module
      CALL update_orbital_parameters(time, declination, sw_toa_down(:))  ! mo_atmland_util


      !> 2.0 read data
      !!
      !! read line by line from file until file end and start again at beginning of file \n
      !! if model%config%flag_climate_forcing_one_year = TRUE then 
      !!        read only the data of the 1st year of climate forcing and start again at beginning of file

      ! day of the year: previous time step
      doy_old = day 

      ! the read statement, by default a list of values in a line is read and then advanced to the beginning of the next line
      lun=20
      READ(lun,*, IOSTAT=stat) year, day, hour, & 
                       swvis_srf_down_r,   lw_srf_down_r, t_air_r,    q_air_r,    &
                       press_srf_r,        rain_r,        snow_r,     wind_air_r, &
                       co2_mixing_ratio_r, co2_dC13_r,    co2_DC14_r,             &
                       nhx_deposition_r,   noy_deposition_r, p_deposition_r 

      ! day of the year: current time step
      doy_new = day 
      
      ! if this flag is set TRUE, check whether the current time step is a new year, if so 
      ! set stat = -1 to indicate "end of file" which invokes the re-reading of the climate forcing input file at line 1
      IF(model%config%flag_climate_forcing_one_year) THEN
        IF( ABS(doy_old - doy_new - 364.0_wp) < eps8 ) THEN
          stat = -1
        END IF
      END IF

      ! if transient forcing_mode is enabled, the climate forcing file may be re-read at line 1 when:
      ! a) the spinup forcing period ends and restarts at start-year, and b) the spinup has finished
      IF(model%config%l_transient_spinup_rewind_forcing) THEN
        stat = -1
      END IF

      ! if end of file (IS_IOSTAT_END(stat) == Test for end-of-file value)
      IF(IS_IOSTAT_END(stat)) THEN 
        ! start reading from line one & read header into dummy
        REWIND(lun)
        CALL read_forcing_header(lun, TRIM(model%config%slm_model_name))
        READ(lun,*) year, day, hour, &
                       swvis_srf_down_r,   lw_srf_down_r, t_air_r,    q_air_r,    &
                       press_srf_r,        rain_r,        snow_r,     wind_air_r, &
                       co2_mixing_ratio_r, co2_dC13_r,    co2_DC14_r,             &
                       nhx_deposition_r, noy_deposition_r, p_deposition_r 
      END IF

      !> 3.0 parameter units are re-scaled
      !!
      q_air_r     = q_air_r     / 1000._wp
      press_srf_r = press_srf_r * 100._wp
      rain_r      = rain_r      / one_day
      snow_r      = snow_r      / one_day

      !> 4.0 forcing timing data and actual forcing values
      !!
      forcing_options%forcing_year  = NINT(year)
      forcing_options%forcing_doy   = NINT(day)
      forcing_options%forcing_hour  = hour

      !> 5.0 update forcing variables based on 'read data'
      !!
      t_air(:)                = t_air_r
      q_air(:)                = q_air_r
      rain(:)                 = rain_r
      snow(:)                 = snow_r
      wind_air(:)             = wind_air_r
      wind_10m(:)             = wind_air_r
      press_srf(:)            = press_srf_r

      lw_srf_down(:)          = lw_srf_down_r
      swvis_srf_down(:)       = swvis_srf_down_r
      cos_zenith_angle(:)     = calc_solar_angle(hgrid%lat(1,1), hgrid%lon(1,1), time) 
      swnir_srf_down(:)       = (1._wp - frac_vis_swtot_srf ) * swvis_srf_down_r 
      swvis_srf_down_r        = frac_vis_swtot_srf * swvis_srf_down_r 
      swpar_srf_down(:)       = frac_par_swvis_srf * swvis_srf_down_r
      swpar_fdiffuse(:)       = calc_sw_fdiffuse(swnir_srf_down(:) + swvis_srf_down_r, sw_toa_down(:), cos_zenith_angle(:)) ! elemental function: calc_sw_fdiffuse()
      
      !! CO2
      ! save the values of the first timestep 
      IF (timestep == 1) THEN
        forcing_options%co2_mixing_ratio_1st_timestep = co2_mixing_ratio_r
        forcing_options%co2_dC13_1st_timestep         = co2_dC13_r
        forcing_options%co2_DC14_1st_timestep         = co2_DC14_r
      ENDIF
      ! if CO2 is constant during simulation: overwrite the forcing values with the value read at the first timestep
      ! (note: '.AND. .NOT. model%config%l_transient_spinup' mean 'and l_transient_spinup is false', i.e., after spinup)
      IF (forcing_options%flag_transient_sim_co2_const .AND. .NOT. model%config%l_transient_spinup) THEN
        co2_mixing_ratio_r  = forcing_options%co2_mixing_ratio_1st_timestep    
        co2_dC13_r          = forcing_options%co2_dC13_1st_timestep
        co2_DC14_r          = forcing_options%co2_DC14_1st_timestep
      ENDIF
      ! pass CO2 value, either the one from file or the one from 1st timestep
      co2_mixing_ratio(:)     = co2_mixing_ratio_r      
      ! calc dC13 and DC14, using either a) the CO2 from file or the CO2 from 1st timestep, or b) the default dC13 & DC14 values
      IF (forcing_options%flag_read_dC13 .OR. &
          (forcing_options%flag_transient_sim_co2_const .AND. .NOT. model%config%l_transient_spinup)) THEN
        hlp1(:)               = calc_mixing_ratio_C13C12(co2_dC13_r)
      ELSE
        hlp1(:)               = calc_mixing_ratio_C13C12(def_co2_deltaC13)
      END IF
      co2_mixing_ratio_C13(:) = co2_mixing_ratio(:) / (1._wp + 1._wp/hlp1(:))
      IF (forcing_options%flag_read_DC14 .OR. &
          (forcing_options%flag_transient_sim_co2_const .AND. .NOT. model%config%l_transient_spinup)) THEN
        hlp1(:)               = calc_mixing_ratio_C14C(co2_dC13_r, co2_DC14_r)
      ELSE
        hlp1(:)               = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14)
      END IF
      co2_mixing_ratio_C14(:) = hlp1(:) * co2_mixing_ratio(:) 

      !! N deposition experiment: add 15N pulse (see ticket #511; aimed at running cruncep sites with u85)
      !!  works only with n_deposition_scheme = dynamic
      !!  no options or such elsewhere, this code is working but is commented out by default
      !!   see also the local variable defined above (all named pulse_exp_* or l_pulse_exp_*)
      !!  no selection of sites just dates, hence, this code runs for all sites, if enabled
      !! -> to run the pulse experiment uncomment the below lines spanning "IF(.NOT. model%config%l_transient_spinup)" & "IF (l_pulse_exp_add_15N)"
      l_pulse_exp_add_15N = .FALSE.
      ! IF(.NOT. model%config%l_transient_spinup) THEN
      !   ! date: 30 Apr 2007, 12PM midday, 120th doy
      !   IF (ABS(forcing_options%forcing_year - 2007.0_wp) < eps4 .AND. &
      !       ABS(forcing_options%forcing_doy  - 120.0_wp ) < eps4 .AND. &
      !       ABS(forcing_options%forcing_hour - 12.0_wp )  < eps4) THEN
      !     l_pulse_exp_add_15N = .TRUE.
      !   ENDIF
      !   ! date: 31 July 2007, 12PM midday, 212th doy
      !   IF (ABS(forcing_options%forcing_year - 2007.0_wp) < eps4 .AND. &
      !       ABS(forcing_options%forcing_doy  - 212.0_wp ) < eps4 .AND. &
      !       ABS(forcing_options%forcing_hour - 12.0_wp )  < eps4) THEN
      !     l_pulse_exp_add_15N = .TRUE.
      !   ENDIF
      !   ! date: 30 Oct 2007, 12PM midday, 303th doy
      !   IF (ABS(forcing_options%forcing_year - 2007.0_wp) < eps4 .AND. &
      !       ABS(forcing_options%forcing_doy  - 303.0_wp ) < eps4 .AND. &
      !       ABS(forcing_options%forcing_hour - 12.0_wp )  < eps4) THEN
      !     l_pulse_exp_add_15N = .TRUE.
      !   ENDIF
      ! ENDIF
      ! ! is always calculated (during the experiement) but not always needed/used
      ! pulse_exp_add     = 7.0_wp * one_day / dtime                                    ! mg N / m2 / timestep !! from 0.07 kg N / ha in 3 doses
      ! pulse_exp_isofrac = pulse_exp_add / (pulse_exp_add + noy_deposition_r)
      ! pulse_exp_delta   = 0._wp * (1 - pulse_exp_isofrac) + pulse_exp_isofrac * 990.0 ! add 99%Natom
      !! N deposition
      SELECT CASE (TRIM(forcing_options%n_deposition_scheme))
        CASE ("none")
          nhx_deposition(:)   = 0.0_wp
          noy_deposition(:)   = 0.0_wp
        CASE ("constant")
          nhx_deposition(:)   = def_n_deposition / molar_mass_N * 1000.0_wp / one_day    ! mg/m2/day -> mumol/m2/sec
          noy_deposition(:)   = def_n_deposition / molar_mass_N * 1000.0_wp / one_day    ! mg/m2/day -> mumol/m2/sec
        CASE ("dynamic")
          nhx_deposition(:)   = nhx_deposition_r / molar_mass_N * 1000.0_wp / one_day    ! mg/m2/day -> mumol/m2/sec
          noy_deposition(:)   = noy_deposition_r / molar_mass_N * 1000.0_wp / one_day    ! mg/m2/day -> mumol/m2/sec
          ! 15N pulse experiment (re-calc noy_deposition)
          IF (l_pulse_exp_add_15N) THEN
            noy_deposition(:) = (noy_deposition_r + pulse_exp_add) / molar_mass_N * 1000.0_wp / one_day ! mg/m2/day -> mumol/m2/sec
          ENDIF
      ENDSELECT
      nhx_n15_deposition(:)   = nhx_deposition(:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
      noy_n15_deposition(:)   = noy_deposition(:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
      ! 15N pulse experiment (re-calc noy_n15_deposition)
      IF (l_pulse_exp_add_15N) THEN
        noy_n15_deposition(:) = noy_deposition(:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(pulse_exp_delta))
      ENDIF

      !! P deposition
      SELECT CASE (TRIM(forcing_options%p_deposition_scheme))
        CASE ("none")
          p_deposition(:)     = 0.0_wp
        CASE ("constant")
          p_deposition(:)     = def_p_deposition / molar_mass_P * 1000.0_wp / one_day    ! mg/m2/day -> mumol/m2/sec
        CASE ("dynamic")
          p_deposition(:)     = p_deposition_r / molar_mass_P * 1000.0_wp / one_day      ! mg/m2/day -> mumol/m2/sec
      ENDSELECT

      !! For lakes
      ! currently no values read from file or calculated, just zet to zero
      ! albedo_lwtr(:)     =         albedo_lwtr(:)    
      ! drag_wtr(:)        =         drag_wtr(:)       
      ! drag_ice(:)        =         drag_ice(:)       
      ! t_acoef_wtr(:)     =         t_acoef_wtr(:)    
      ! t_bcoef_wtr(:)     =         t_bcoef_wtr(:)    
      ! q_acoef_wtr(:)     =         q_acoef_wtr(:)    
      ! q_bcoef_wtr(:)     =         q_bcoef_wtr(:)    
      ! t_acoef_ice(:)     =         t_acoef_ice(:)    
      ! t_bcoef_ice(:)     =         t_bcoef_ice(:)    
      ! q_acoef_ice(:)     =         q_acoef_ice(:)    
      ! q_bcoef_ice(:)     =         q_bcoef_ice(:)    



    !! #############################################################################################
    !! SOIL
    !! #############################################################################################
    CASE('soil')

      !> 1.0 read data
      !!
      !! the read statement, by default a list of values in a line is read and then advanced to the beginning of the next line
      lun=20
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

      ! if transient forcing_mode is enabled, the climate forcing file may be re-read at line 1 when:
      ! a) the spinup forcing period ends and restarts at start-year, and b) the spinup has finished
      IF(model%config%l_transient_spinup_rewind_forcing) THEN
        stat = -1
      END IF

      ! set stat = -1 to indicate "end of file" which invokes the re-reading of the climate forcing input file at line 1
      IF(IS_IOSTAT_END(stat)) THEN 
        ! start reading from line one & read header into dummy
        REWIND(lun)
        CALL read_forcing_header(lun, TRIM(model%config%slm_model_name))
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

      lun=21
      READ(lun,*, IOSTAT=stat) time_hlp, & 
                       t_soil_sl_r(1:nsoil_sb), w_soil_sl_r(1:nsoil_sb), &
                       percolation_sl_r(1:nsoil_sb), root_fraction_sl_r(1:nsoil_sb), &
                       veg_pool_fine_root_carbon_r, &
                       nhx_deposition_r, noy_deposition_r, &
                       nhx_dn15_deposition_r, noy_dn15_deposition_r, &
                       p_deposition_r, &
                       f_n_demand_r, f_p_demand_r

      ! if transient forcing_mode is enabled, the climate forcing file may be re-read at line 1 when:
      ! a) the spinup forcing period ends and restarts at start-year, and b) the spinup has finished
      IF(model%config%l_transient_spinup_rewind_forcing) THEN
        stat = -1
      END IF

      IF(IS_IOSTAT_END(stat)) THEN 
        ! start reading from line one & read header into dummy
        REWIND(lun)
        CALL read_forcing_header(lun, TRIM(model%config%slm_model_name))
        READ(lun,*) time_hlp, &
                       t_soil_sl_r(1:nsoil_sb), w_soil_sl_r(1:nsoil_sb), &
                       percolation_sl_r(1:nsoil_sb), root_fraction_sl_r(1:nsoil_sb), &
                       veg_pool_fine_root_carbon_r, &
                       nhx_deposition_r, noy_deposition_r, &
                       nhx_dn15_deposition_r, noy_dn15_deposition_r, &
                       p_deposition_r, &
                       f_n_demand_r, f_p_demand_r
      END IF

      !> 2.0 pass values read from file (i.e. scalar values) to ics-ice
      !!
      veg_litterfall%leaf%carbon(:)     = leaf_carbon_r / 1000._wp
      veg_litterfall%leaf%nitrogen(:)   = leaf_nitrogen_r / 1000000._wp
      veg_litterfall%leaf%phosphorus(:) = leaf_phosphorus_r / 1000000._wp   
      veg_litterfall%leaf%carbon13(:)   = veg_litterfall%leaf%carbon(:) / ( 1._wp + 1._wp / & 
                                            calc_mixing_ratio_C13C12(leaf_dC13_r))  
      veg_litterfall%leaf%carbon14(:)   = veg_litterfall%leaf%carbon(:) * & 
                                            calc_mixing_ratio_C14C(leaf_dC13_r,leaf_DC14_r) 
      veg_litterfall%leaf%nitrogen15(:) = veg_litterfall%leaf%nitrogen(:) / ( 1._wp + 1._wp / & 
                                            calc_mixing_ratio_N15N14(leaf_dN15_r))

      veg_litterfall%fine_root%carbon(:)     = fine_root_carbon_r / 1000._wp
      veg_litterfall%fine_root%nitrogen(:)   = fine_root_nitrogen_r / 1000000._wp
      veg_litterfall%fine_root%phosphorus(:) = fine_root_phosphorus_r / 1000000._wp
      veg_litterfall%fine_root%carbon13(:)   = veg_litterfall%fine_root%carbon(:) / ( 1._wp + 1._wp / & 
                                                 calc_mixing_ratio_C13C12(fine_root_dC13_r)) 
      veg_litterfall%fine_root%carbon14(:)   = veg_litterfall%fine_root%carbon(:) * & 
                                                 calc_mixing_ratio_C14C(fine_root_dC13_r,fine_root_DC14_r) 
      veg_litterfall%fine_root%nitrogen15(:) = veg_litterfall%fine_root%nitrogen(:) / ( 1._wp + 1._wp / & 
                                                 calc_mixing_ratio_N15N14(fine_root_dN15_r)) 

      veg_litterfall%coarse_root%carbon(:)     = coarse_root_carbon_r / 1000._wp
      veg_litterfall%coarse_root%nitrogen(:)   = coarse_root_nitrogen_r / 1000000._wp
      veg_litterfall%coarse_root%phosphorus(:) = coarse_root_phosphorus_r / 1000000._wp
      veg_litterfall%coarse_root%carbon13(:)   = veg_litterfall%coarse_root%carbon(:) / ( 1._wp + 1._wp / & 
                                                   calc_mixing_ratio_C13C12(coarse_root_dC13_r)) 
      veg_litterfall%coarse_root%carbon14(:)   = veg_litterfall%coarse_root%carbon(:) * & 
                                                   calc_mixing_ratio_C14C(coarse_root_dC13_r,coarse_root_DC14_r) 
      veg_litterfall%coarse_root%nitrogen15(:) = veg_litterfall%coarse_root%nitrogen(:) / ( 1._wp + 1._wp / & 
                                                   calc_mixing_ratio_N15N14(coarse_root_dN15_r))

      veg_litterfall%sap_wood%carbon(:)        = sap_wood_carbon_r / 1000._wp
      veg_litterfall%sap_wood%nitrogen(:)      = sap_wood_nitrogen_r / 1000000._wp
      veg_litterfall%sap_wood%phosphorus(:)    = sap_wood_phosphorus_r / 1000000._wp
      veg_litterfall%sap_wood%carbon13(:)      = veg_litterfall%sap_wood%carbon(:) / ( 1._wp + 1._wp / & 
                                                   calc_mixing_ratio_C13C12(sap_wood_dC13_r))
      veg_litterfall%sap_wood%carbon14(:)      = veg_litterfall%sap_wood%carbon(:) * & 
                                                   calc_mixing_ratio_C14C(sap_wood_dC13_r,sap_wood_DC14_r) 
      veg_litterfall%sap_wood%nitrogen15(:)    = veg_litterfall%sap_wood%nitrogen(:) / ( 1._wp + 1._wp / & 
                                                   calc_mixing_ratio_N15N14(sap_wood_dN15_r))

      veg_litterfall%heart_wood%carbon(:)      = heart_wood_carbon_r / 1000._wp
      veg_litterfall%heart_wood%nitrogen(:)    = heart_wood_nitrogen_r / 1000000._wp
      veg_litterfall%heart_wood%phosphorus(:)  = heart_wood_phosphorus_r / 1000000._wp 
      veg_litterfall%heart_wood%carbon13(:)    = veg_litterfall%heart_wood%carbon(:) / ( 1._wp + 1._wp / & 
                                                   calc_mixing_ratio_C13C12(heart_wood_dC13_r)) 
      veg_litterfall%heart_wood%carbon14(:)    = veg_litterfall%heart_wood%carbon(:) * & 
                                                   calc_mixing_ratio_C14C(heart_wood_dC13_r,heart_wood_DC14_r) 
      veg_litterfall%heart_wood%nitrogen15(:)  = veg_litterfall%heart_wood%nitrogen(:) / ( 1._wp + 1._wp / & 
                                                   calc_mixing_ratio_N15N14(heart_wood_dN15_r))

      veg_litterfall%labile%carbon(:)       = labile_carbon_r / 1000._wp
      veg_litterfall%labile%nitrogen(:)     = labile_nitrogen_r / 1000000._wp
      veg_litterfall%labile%phosphorus(:)   = labile_phosphorus_r / 1000000._wp 
      veg_litterfall%labile%carbon13(:)     = veg_litterfall%labile%carbon(:) / ( 1._wp + 1._wp / & 
                                                calc_mixing_ratio_C13C12(labile_dC13_r))
      veg_litterfall%labile%carbon14(:)     = veg_litterfall%labile%carbon(:) * & 
                                                calc_mixing_ratio_C14C(labile_dC13_r,labile_DC14_r) 
      veg_litterfall%labile%nitrogen15(:)   = veg_litterfall%labile%nitrogen(:) / ( 1._wp + 1._wp / & 
                                                calc_mixing_ratio_N15N14(labile_dN15_r))

      veg_litterfall%reserve%carbon(:)      = reserve_carbon_r / 1000._wp
      veg_litterfall%reserve%nitrogen(:)    = reserve_nitrogen_r / 1000000._wp
      veg_litterfall%reserve%phosphorus(:)  = reserve_phosphorus_r / 1000000._wp
      veg_litterfall%reserve%carbon13(:)    = veg_litterfall%reserve%carbon(:) / ( 1._wp + 1._wp / & 
                                                calc_mixing_ratio_C13C12(reserve_dC13_r)) 
      veg_litterfall%reserve%carbon14(:)    = veg_litterfall%reserve%carbon(:) * & 
                                                calc_mixing_ratio_C14C(reserve_dC13_r,reserve_DC14_r) 
      veg_litterfall%reserve%nitrogen15(:)  = veg_litterfall%reserve%nitrogen(:) / ( 1._wp + 1._wp / & 
                                                calc_mixing_ratio_N15N14(reserve_dN15_r))

      veg_litterfall%fruit%carbon(:)        = fruit_carbon_r / 1000._wp
      veg_litterfall%fruit%nitrogen(:)      = fruit_nitrogen_r / 1000000._wp
      veg_litterfall%fruit%phosphorus(:)    = fruit_phosphorus_r / 1000000._wp
      veg_litterfall%fruit%carbon13(:)      = veg_litterfall%fruit%carbon(:) / ( 1._wp + 1._wp / & 
                                                calc_mixing_ratio_C13C12(fruit_dC13_r))
      veg_litterfall%fruit%carbon14(:)      = veg_litterfall%fruit%carbon(:) * & 
                                                calc_mixing_ratio_C14C(fruit_dC13_r,fruit_DC14_r) 
      veg_litterfall%fruit%nitrogen15(:)    = veg_litterfall%fruit%nitrogen(:) / ( 1._wp + 1._wp / & 
                                                calc_mixing_ratio_N15N14(fruit_dN15_r))

      veg_litterfall%seed_bed%carbon(:)     = seed_bed_carbon_r / 1000._wp
      veg_litterfall%seed_bed%nitrogen(:)   = seed_bed_nitrogen_r / 1000000._wp
      veg_litterfall%seed_bed%phosphorus(:) = seed_bed_phosphorus_r / 1000000._wp
      veg_litterfall%seed_bed%carbon13(:)   = veg_litterfall%seed_bed%carbon(:) / ( 1._wp + 1._wp / & 
                                                calc_mixing_ratio_C13C12(seed_bed_dC13_r)) 
      veg_litterfall%seed_bed%carbon14(:)   = veg_litterfall%seed_bed%carbon(:) * & 
                                                calc_mixing_ratio_C14C(seed_bed_dC13_r,seed_bed_DC14_r) 
      veg_litterfall%seed_bed%nitrogen15(:) = veg_litterfall%seed_bed%nitrogen(:) / ( 1._wp + 1._wp / & 
                                                calc_mixing_ratio_N15N14(seed_bed_dN15_r))

      DO isoil = 1,nsoil_sb
        t_soil_sl(:,isoil)                  = t_soil_sl_r(isoil)
        w_soil_sl(:,isoil)                  = w_soil_sl_r(isoil)
        percolation_sl(:,isoil)             = percolation_sl_r(isoil)
        root_fraction_sl(:,isoil)           = root_fraction_sl_r(isoil)
      ENDDO
      veg_pool%fine_root%carbon(:)          = veg_pool_fine_root_carbon_r

      nhx_deposition(:)     = nhx_deposition_r
      noy_deposition(:)     = noy_deposition_r
      nhx_n15_deposition(:) = nhx_deposition(:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(nhx_dn15_deposition_r))
      nhx_n15_deposition(:) = noy_deposition(:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(noy_dn15_deposition_r))
      p_deposition(:)       = p_deposition_r
      f_n_demand(:)         = f_n_demand_r
      f_p_demand(:)         = f_p_demand_r

    END SELECT  ! land plant canopy vs soil

! Pass back values to process memory
dsl4slm_passback_vegpool_var(pool,fine_root,carbon)
#define Passback_veg_litterfall
#include "slm_pools_a_incl.f90"
#undef Passback_veg_litterfall

  END SUBROUTINE update_forcing


  !-----------------------------------------------------------------------------------------------------
  !> open the file test_bed/plant/climate.dat, only IF model option is set to READ
  !! 
  !! called in the driver
  !! 
  !! note: lun is set to 20 for land plant canopy input; and 20/21 for soil input
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE open_forcing_input_file(slm_model_name)

    USE mo_exception,                       ONLY: finish

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CHARACTER(len=*)    ,INTENT(in)     :: slm_model_name
    ! ---------------------------
    ! 0.2 Local
    LOGICAL         :: file_exists
    INTEGER         :: lun 
    CHARACTER(100)  :: fname
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':open_forcing_input_file'


    SELECT CASE(slm_model_name)
      CASE('land','plant','canopy')
        !> 1.1 open file for READ
        !! 
        lun   = 20
        fname = 'climate.dat'
        ! check if file exists
        INQUIRE(FILE=fname, EXIST=file_exists)
        ! open if fname == True
        IF(file_exists) THEN
          OPEN (UNIT=lun, FILE=fname, ACTION="read", STATUS="old") ! status=old: file exists already
        !> 1.2 read the two header lines into dummy variables
        !! to advance to line 3
          CALL read_forcing_header(lun, slm_model_name)
        ELSE
          CALL finish(TRIM(routine), 'climate forcing: climate.dat could not be opened, file does not exist')
        END IF

      CASE('soil')
        !> 2.1 open file for READ
        !! 
        lun   = 20
        fname = 'soil_forcing_litterfall.dat'
        ! check if file exists
        INQUIRE(FILE=fname, EXIST=file_exists)
        ! open if fname == True
        IF(file_exists) THEN
          OPEN (UNIT=lun, FILE=fname, ACTION="read", STATUS="old") ! status=old: file exists already
        !> 2.2 read the header line into dummy variable
        !! to advance to line 2
          CALL read_forcing_header(lun, slm_model_name)
        ELSE
          CALL finish(TRIM(routine), 'soil forcing: soil_forcing_litterfall.dat could not be opened, file does not exist')
        END IF

        !> 3.1 open file for READ
        !! 
        lun   = 21
        fname = 'soil_forcing_soilphysics.dat'
        ! check if file exists
        INQUIRE(FILE=fname, EXIST=file_exists)
        ! open if fname == True
        IF(file_exists) THEN
          OPEN (UNIT=lun, FILE=fname, ACTION="read", STATUS="old") ! status=old: file exists already
        !> 3.2 read the header line into dummy variable
        !! to advance to line 2
          CALL read_forcing_header(lun, slm_model_name)
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
  SUBROUTINE close_forcing_input_file(slm_model_name)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CHARACTER(len=*)    ,INTENT(in)     :: slm_model_name
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':close_forcing_input_file'

    SELECT CASE(slm_model_name)
      CASE('land','plant','canopy')
        CLOSE(20)    ! climate

      CASE('soil')
        CLOSE(20)   ! litterfall
        CLOSE(21)   ! soilphysics
    END SELECT

  END SUBROUTINE close_forcing_input_file


  !-----------------------------------------------------------------------------------------------------
  !> read only the header into a dummy variable to advanve the READ to the 3rd line in the file
  !! 
  !! consider this a helper subroutine to the read_forcing subroutine
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE read_forcing_header(lun, slm_model_name)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER             ,INTENT(in)     :: lun
    CHARACTER(len=*)    ,INTENT(in)     :: slm_model_name
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(200)      :: header_dummy1 = "header"
    CHARACTER(500)      :: header_dummy2 = "header_long"
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':read_forcing_header'

    SELECT CASE(slm_model_name)
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
