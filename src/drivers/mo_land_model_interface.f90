!>
!! @brief 
!!  interface for land model \n  
!!  run SLM for one time step
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!!  First version by Sönke Zaehle (2015-12-04) \n
!!  Reworked by Jan Engel (2016 July)
!!
MODULE mo_land_model_interface

#include "dsl_jsbach.inc"

  USE mo_kind

  IMPLICIT NONE


  PRIVATE
  PUBLIC :: land_model_interface


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> land model interface
  !! 
  !! selection of the model via the namlist.slm - BASE_CTL:slm_model_name
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE land_model_interface(tile, &                ! one specific tile
                                  l_forcing, &           ! forcing on/off
                                  time_day, time_year, timestep_prop_day, &
                                  !! atm2land memory - INTENT(in) - start
                                  t_air              , &
                                  q_air              , &
                                  rain               , &
                                  snow               , &
                                  wind_air           , &
                                  wind_10m           , &
                                  lw_srf_down        , &
                                  swvis_srf_down     , &
                                  swnir_srf_down     , &
                                  swpar_srf_down     , &
                                  swpar_fdiffuse     , &
                                  press_srf          , &
                                  drag_srf           , &
                                  t_acoef            , &
                                  t_bcoef            , &
                                  q_acoef            , &
                                  q_bcoef            , &
                                  pch                , &
                                  cos_zenith_angle   , &
                                  co2_mixing_ratio     , &
                                  co2_mixing_ratio_C13 , &
                                  co2_mixing_ratio_C14 , &
                                  ga                   , &
                                  ! For lakes:
                                  albedo_lwtr        , &
                                  drag_wtr           , &
                                  drag_ice           , &
                                  t_acoef_wtr        , &
                                  t_bcoef_wtr        , &
                                  q_acoef_wtr        , &
                                  q_bcoef_wtr        , &
                                  t_acoef_ice        , &
                                  t_bcoef_ice        , &
                                  q_acoef_ice        , &
                                  q_bcoef_ice        , &
                                 !! atm2land memory - end
                                 !! land2atm memory - INTENT(out) - start
                                  t_srf             ,   &
                                  t_eff_srf         ,   &
                                  qsat_srf          ,   &
                                  fact_q_air        ,   &
                                  fact_qsat_srf     ,   &
                                  zh_srf            ,   &
                                  zm_srf            ,   &
                                  alb_vis_dif       ,   &
                                  alb_nir_dif       ,   &
                                  alb_vis_dir       ,   &
                                  alb_nir_dir       ,   &
                                  albedo            ,   &
                                  sensible_heat_flx ,   &
                                  latent_heat_flx   ,   &
                                  net_radiation     ,   &
                                  vpd               ,   &
                                  ! For lakes:
                                  latent_hflx_wtr   ,   &
                                  sensible_hflx_wtr ,   &
                                  latent_hflx_ice   ,   &
                                  sensible_hflx_ice     &
                                  !! land2atm memory - end
                                  )


    USE mo_base_memory,                   ONLY: t_slm_tile
    USE mo_base_config_class,             ONLY: base_config
    USE mo_shared_update_time_averages,   ONLY: update_time_averages  ! includes daytime and moving average
    USE mo_atmland_interface,             ONLY: update_atm2land, update_land2atm
    USE mo_veg_output,                    ONLY: write_output_pool, write_output_flux

    ! Process interfaces:
    USE mo_rad_interface
    USE mo_veg_interface
    USE mo_assimi_interface
    !USE mo_soil_interface


    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)      :: tile              !< data structure for dynamic vegetation and soil variables
    LOGICAL,INTENT(in)                  :: l_forcing         !< IF .TRUE. THEN update_forcing(), update_atm2land(), update_land2atm()
    REAL(wp),INTENT(in)                 :: time_day          !< displays the 'real time' in days; time_day == 1.0 after one day has passed
    INTEGER,INTENT(in)                  :: time_year         !< displays the 'real time' in years; the 1st year of a simulation == 1 (not zero)
    REAL(wp),INTENT(in)                 :: timestep_prop_day !< proportion of one day that is simulated with one timestep (dt/one_day)
    REAL(wp),INTENT(in)           :: &               !! atm2land memory - INTENT(in) - start
                                  t_air              , &
                                  q_air              , &
                                  rain               , &
                                  snow               , &
                                  wind_air           , &
                                  wind_10m           , &
                                  lw_srf_down        , &
                                  swvis_srf_down     , &
                                  swnir_srf_down     , &
                                  swpar_srf_down     , &
                                  swpar_fdiffuse     , &
                                  press_srf          , &
                                  drag_srf           , &
                                  t_acoef            , &
                                  t_bcoef            , &
                                  q_acoef            , &
                                  q_bcoef            , &
                                  pch                , &
                                  cos_zenith_angle   , &
                                  co2_mixing_ratio     , &
                                  co2_mixing_ratio_C13 , &
                                  co2_mixing_ratio_C14 , &
                                  ga                   , &
                                  ! For lakes:
                                  albedo_lwtr        , &
                                  drag_wtr           , &
                                  drag_ice           , &
                                  t_acoef_wtr        , &
                                  t_bcoef_wtr        , &
                                  q_acoef_wtr        , &
                                  q_bcoef_wtr        , &
                                  t_acoef_ice        , &
                                  t_bcoef_ice        , &
                                  q_acoef_ice        , &
                                  q_bcoef_ice        
                                                     !! atm2land memory - end
    REAL(wp),INTENT(out)          :: &              !! land2atm memory - INTENT(out) - start
                                  t_srf             ,   &
                                  t_eff_srf         ,   &
                                  qsat_srf          ,   &
                                  fact_q_air        ,   &
                                  fact_qsat_srf     ,   &
                                  zh_srf            ,   &
                                  zm_srf            ,   &
                                  alb_vis_dif       ,   &
                                  alb_nir_dif       ,   &
                                  alb_vis_dir       ,   &
                                  alb_nir_dir       ,   &
                                  albedo            ,   &
                                  sensible_heat_flx ,   &
                                  latent_heat_flx   ,   &
                                  net_radiation     ,   &
                                  vpd               ,   &
                                  ! For lakes:
                                  latent_hflx_wtr   ,   &
                                  sensible_hflx_wtr ,   &
                                  latent_hflx_ice   ,   &
                                  sensible_hflx_ice !! land2atm memory - end
    ! 0.2 Local
    ! ...


    !> pass atmland variables from forcing input (driver) to atm2land
    !! 
    IF (l_forcing) THEN       ! l_forcing is true
      CALL update_atm2land( &
      ! the sequence of input variables is consistent with the sequence of var in t_atm2land_memory & atmland_interface
      tile                                         , &
      t_air                   = t_air              , &
      q_air                   = q_air              , &
      rain                    = rain               , &
      snow                    = snow               , &
      wind_air                = wind_air           , &
      wind_10m                = wind_10m           , &
      lw_srf_down             = lw_srf_down        , &
      swvis_srf_down          = swvis_srf_down     , &
      swnir_srf_down          = swnir_srf_down     , &
      swpar_srf_down          = swpar_srf_down     , &
      swpar_fdiffuse          = swpar_fdiffuse     , &
      press_srf               = press_srf          , &
      drag_srf                = drag_srf           , &
      t_acoef                 = t_acoef            , &
      t_bcoef                 = t_bcoef            , &
      q_acoef                 = q_acoef            , &
      q_bcoef                 = q_bcoef            , &
      pch                     = pch                , &
      cos_zenith_angle        = cos_zenith_angle   , &
      co2_mixing_ratio        = co2_mixing_ratio     , &
      co2_mixing_ratio_C13    = co2_mixing_ratio_C13 , &
      co2_mixing_ratio_C14    = co2_mixing_ratio_C14 , &
      ga                      = ga                   , &
      ! For lakes:
      albedo_lwtr             = albedo_lwtr        , &
      drag_wtr                = drag_wtr           , &
      drag_ice                = drag_ice           , &
      t_acoef_wtr             = t_acoef_wtr        , &
      t_bcoef_wtr             = t_bcoef_wtr        , &
      q_acoef_wtr             = q_acoef_wtr        , &
      q_bcoef_wtr             = q_bcoef_wtr        , &
      t_acoef_ice             = t_acoef_ice        , &
      t_bcoef_ice             = t_bcoef_ice        , &
      q_acoef_ice             = q_acoef_ice        , &
      q_bcoef_ice             = q_bcoef_ice          &
      & )
    END IF


    !> @todo mo_land_model_interface:land_model_interface put all Tasks into a Queue accordingly to the model selection
    !! (land plant soil canopy test_canopy)
    SELECT CASE(TRIM(base_config%slm_model_name))
    !! land ---------------------------------------------------------------------------------------
    CASE('land')
    !! set vegetation fluxes to zero
    !! 
      CALL set_zero_vegetation_fluxes(tile)   ! mo_assimi_update_canopy_fluxes

    !! run Tasks
    !!
      CALL update_radiation       (tile)
      CALL update_vegsoil_bgc     (tile)
      CALL update_phenology       (tile)
      CALL update_canopy_fluxes   (tile)
      CALL calc_growth            (tile)
      CALL calc_turnover          (tile)
      CALL update_veg_dynamics    (tile)
      CALL update_vegsoil_pool    (tile)
      CALL update_vegsoil_physics (tile)
      CALL update_veg_pools       (tile)

    !! accumulation over time - daytime & moving average
    !! 
      CALL update_time_averages(tile, time_day, timestep_prop_day)    ! see shared

    !! output pool: write output in certain intervals, see mo_veg_output
    !! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
    !! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
    !! 
      CALL write_output_pool(tile, time_day, time_year, 1) ! 1 == itile

    !! output flux: write output in certain intervals, see mo_veg_output
    !! 
    !! 
      CALL write_output_flux(tile, time_day, time_year, 1) ! 1 == itile 

    CASE('plant')
    !! plant --------------------------------------------------------------------------------------
    ! currently same as land !
    !! set vegetation fluxes to zero
    !! 
      CALL set_zero_vegetation_fluxes(tile)   ! mo_assimi_update_canopy_fluxes

    !! run Tasks
    !!
      CALL update_radiation       (tile)
      CALL update_vegsoil_bgc     (tile)
      CALL update_phenology       (tile)
      CALL update_canopy_fluxes   (tile)
      CALL calc_growth            (tile)
      CALL calc_turnover          (tile)
      CALL update_veg_dynamics    (tile)
      CALL update_vegsoil_pool    (tile)
      CALL update_vegsoil_physics (tile)
      CALL update_veg_pools       (tile)

    !! accumulation over time - daytime & moving average
    !! 
      CALL update_time_averages(tile, time_day, timestep_prop_day)    ! see shared

    !! output pool: write output in certain intervals, see mo_veg_output
    !! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
    !! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
    !! 
      CALL write_output_pool(tile, time_day, time_year, 1) ! 1 == itile

    !! output flux: write output in certain intervals, see mo_veg_output
    !! 
    !! 
      CALL write_output_flux(tile, time_day, time_year, 1) ! 1 == itile 

    CASE('soil')
    !! soil ---------------------------------------------------------------------------------------

    CASE('canopy')
    !! canopy -------------------------------------------------------------------------------------
    !! set vegetation fluxes to zero
    !! 
      CALL set_zero_vegetation_fluxes(tile)   ! mo_assimi_update_canopy_fluxes

    !! run Tasks
    !!
      CALL update_veg_pools       (tile)
      CALL update_radiation       (tile)
      CALL update_canopy_fluxes   (tile)
      CALL update_vegsoil_physics (tile)

    !! accumulation over time - daytime & moving average
    !! 
      CALL update_time_averages(tile, time_day, timestep_prop_day)    ! see shared

    !! output fluxes -> filename: fort.20
    !! output also some atm2land variables
    !! 
      CALL canopy_model_write_data(tile, time_day, time_year, t_air, swpar_srf_down)                           

    CASE('test_canopy')
    !! test_canopy --------------------------------------------------------------------------------
      CALL run_test_canopy_model  (tile)

    CASE DEFAULT
      CALL finish('land_model_interface','Invalid SLM model name. Try either land, plant, soil, canopy or test_canopy.')
    END SELECT  ! select SLM model !

    !> pass atmland variables from land2atm memory back to the driver
    !! 
    IF (l_forcing) THEN       ! l_forcing is true if update_atm2land() is called (depends on the SLM model selected)
      CALL update_land2atm( &
      ! the sequence of input variables is consistent with the sequence of var in t_land2atm_memory !
      tile                                        ,   &
      t_srf                   = t_srf             ,   &
      t_eff_srf               = t_eff_srf         ,   &
      qsat_srf                = qsat_srf          ,   &
      fact_q_air              = fact_q_air        ,   &
      fact_qsat_srf           = fact_qsat_srf     ,   &
      zh_srf                  = zh_srf            ,   &
      zm_srf                  = zm_srf            ,   &
      alb_vis_dif             = alb_vis_dif       ,   &
      alb_nir_dif             = alb_nir_dif       ,   &
      alb_vis_dir             = alb_vis_dir       ,   &
      alb_nir_dir             = alb_nir_dir       ,   &
      albedo                  = albedo            ,   &
      sensible_heat_flx       = sensible_heat_flx ,   &
      latent_heat_flx         = latent_heat_flx   ,   &
      net_radiation           = net_radiation     ,   &
      vpd                     = vpd               ,   &
      ! For lakes:
      latent_hflx_wtr         = latent_hflx_wtr   ,   &
      sensible_hflx_wtr       = sensible_hflx_wtr ,   &
      latent_hflx_ice         = latent_hflx_ice   ,   &
      sensible_hflx_ice       = sensible_hflx_ice     &
      & )
    END IF
  END SUBROUTINE land_model_interface


  !-----------------------------------------------------------------------------------------------------
  !> Write results and data from canopy model
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE canopy_model_write_data(tile, time_day, time_year, t_air, swpar_srf_down)

    USE mo_atmland_constants,     ONLY: rad2par
    USE mo_base_constants,        ONLY: one_year
    USE mo_base_memory,           ONLY: t_slm_tile
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)      :: tile              !< data structure for dynamic vegetation and soil variables
    REAL(wp),        INTENT(in)         :: time_day          !< displays the 'real time' in days; time_day == 1.0 after one day has passed
    INTEGER,         INTENT(in)         :: time_year         !< displays the 'real time' in years; the 1st year of a simulation == 1 (not zero)
    REAL(wp),        INTENT(in)         :: t_air             !< air temp
    REAL(wp),        INTENT(in)         :: swpar_srf_down    !< short wave light at surface level ?
    ! 0.2 Local
    ! ..
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          !   TYPE(t_veg_memory), POINTER      :: veg__mem  
    TYPE(t_veg_memory), POINTER         :: veg        ! local pointer to veg__mem for SLM usability
    dsl4jsb_Get_memory(VEG_)
!    SELECT TYPE (mm => tile%mem(VEG_)%p) 
!      TYPE IS(t_veg_memory)
!        veg__mem => mm
!    END SELECT
    veg => veg__mem 


      WRITE(20,'(7F16.7)')  time_day+time_year*one_year, t_air, swpar_srf_down * rad2par, &
                            veg%gpp%ptr, veg%apar%ptr, veg%gc%ptr, veg%ci%ptr
      WRITE(21,'(20F16.7)') veg%par_sunlit_cl%ptr
      WRITE(22,'(20F16.7)') veg%par_shaded_cl%ptr
      WRITE(23,'(20F16.7)') veg%fleaf_sunlit_cl%ptr
      WRITE(24,'(20F16.7)') veg%par_sunlit_cl%ptr(:,1) * veg%fleaf_sunlit_cl%ptr(:,1) + &
                            & veg%par_shaded_cl%ptr(:,1) + (1._wp-veg%fleaf_sunlit_cl%ptr(:,1))

  END SUBROUTINE canopy_model_write_data


  !-----------------------------------------------------------------------------------------------------
  !> Run the test_canopy model from the land model interface
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE run_test_canopy_model(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model
    USE mo_base_constants,        ONLY: npft,ncanop
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
    USE mo_atmland_constants,     ONLY: rad2par, frac_vis_swtot_srf, frac_par_swvis_srf,eps_vpd, standard_press_srf, &
                                        def_co2_mixing_ratio
    USE mo_atmland_util,          ONLY: calc_spec_humidity_sat
    USE mo_shared_constants,      ONLY: molar_mass_N, Tzero 

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    ! Process interfaces:
    USE mo_atmland_interface                ! update_atmland
    USE mo_rad_interface                    ! update_radiation
    USE mo_assimi_interface                 ! update_canopy_fluxes
    USE mo_veg_canopy
    ! get parameters defined in init_slm_test_canopy()

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_slm_tile),       INTENT(inout)     :: tile       !< data structure for dynamic vegetation and soil variables
    ! 0.2 Local
    TYPE(t_model),          POINTER           :: model
    TYPE(t_lctlib_element), POINTER           :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp) ::   t_air , &
                  press_srf , &
                  q_air      , &
                  co2_mixing_ratio , &
                  rain          , &
                  snow           , &
                  swvis_srf_down , &
                  swnir_srf_down , &
                  swpar_srf_down , &
                  swpar_fdiffuse , &
                  wind_air       , &
                  wind_10m       , &
                  lw_srf_down    , &
                  drag_srf       , &
                  t_acoef        , &
                  t_bcoef        , &
                  q_acoef        , &
                  q_bcoef        , &
                  pch            , &
                  cos_zenith_angle , &
                  ga             !, &
                  !latitude      , &
                  !longitude     , &
    REAL(wp)      :: vpd
    REAL(wp)      :: hlp1
    INTEGER       :: i,lun
    INTEGER       :: ref_pft                          ! PFT identifier

    ! parameters (specific to the test_canopy model)
    REAL(wp),SAVE :: ref_t_air = 25.0_wp              ! temperature in deg C 
    REAL(wp),SAVE :: ref_vpd = 1000._wp               ! vpd in Pa
    REAL(wp),SAVE :: ref_swdown = 652.173913_wp       ! PAR in W m-2
    REAL(wp),SAVE :: ref_fdiffuse = 0.2_wp            ! diffuse fraction
    REAL(wp),SAVE :: ref_cos_angle = 0.8_wp           ! cosine of zenith angle
    REAL(wp),SAVE :: ref_co2_mixing_ratio = 380.0_wp  ! co2 in ppm
    REAL(wp),SAVE :: ref_wind = 1.0_wp                ! wind in m/s
    REAL(wp),SAVE :: ref_theta = 100._wp              ! soil water availability in %
    REAL(wp),SAVE :: ref_lai = 3._wp                  ! LAI
    REAL(wp),SAVE :: ref_nleaf = 2.0_wp               ! N leaf per unit leaf area in gN/m2
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          !   TYPE(t_veg_memory), POINTER      :: veg__mem  
    TYPE(t_veg_memory), POINTER         :: veg        ! local pointer to veg__mem for SLM usability
    dsl4jsb_Get_memory(VEG_)
!    SELECT TYPE (mm => tile%mem(VEG_)%p) 
!      TYPE IS(t_veg_memory)
!        veg__mem => mm
!    END SELECT
    veg    => veg__mem 
    model  => Get_model(tile%owner_model_id)


    ! if pftid is set in the namelist use this pftid
    IF(tile%lcts%lib_id>0)THEN
      ref_pft = tile%lcts%lib_id
    ! else use pftid 1
    ELSE
      ref_pft = 0
      tile%lcts%lib_id = 1
    ENDIF
    lctlib  => model%lctlib(tile%lcts%lib_id)


    !> 1.0 Init
    !! 

    !> 1.1 Initalisation of met-forcing which is assumed constant 
    !! 
    !!  
    t_air            = ref_t_air + TZero 
    press_srf        = standard_press_srf
    q_air            = calc_spec_humidity_sat(t_air,press_srf) - ref_vpd / press_srf * eps_vpd
    co2_mixing_ratio = def_co2_mixing_ratio   ! def_co2_mixing_ratio is used instead of ref_co2_mixing_ratio 
    rain             = 0.0_wp
    snow             = 0.0_wp
    swpar_srf_down   = ref_swdown 
    swvis_srf_down   = swpar_srf_down / frac_par_swvis_srf 
    swnir_srf_down   = swvis_srf_down * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf
    swpar_fdiffuse   = ref_fdiffuse
    wind_air         = ref_wind
    wind_10m         = wind_air
    lw_srf_down      = 0.0_wp
    drag_srf         = 1.0_wp
    t_acoef          = 0.0_wp
    t_bcoef          = 0.0_wp
    q_acoef          = 0.0_wp
    q_bcoef          = 0.0_wp
    pch              = 0.0_wp
    cos_zenith_angle = ref_cos_angle
    ga               = 0.4_wp / ( drag_srf * wind_10m )
    !latitude                  = 50.0_wp
    !longitude                 = 0.0_wp
    CALL update_atm2land( &
           tile                                         , &
           t_air                   = t_air              , &
           q_air                   = q_air              , &
           wind_air                = wind_air           , &
           wind_10m                = wind_10m           , &
           swvis_srf_down          = swvis_srf_down     , &
           swnir_srf_down          = swnir_srf_down     , &
           swpar_srf_down          = swpar_srf_down     , &
           swpar_fdiffuse          = swpar_fdiffuse     , &
           press_srf               = press_srf          , &
           drag_srf                = drag_srf           , &
           cos_zenith_angle        = cos_zenith_angle   , &
           co2_mixing_ratio        = co2_mixing_ratio     , &
           ga                      = ga                   )

    !> 1.2 Initialisation of soil and vegetation (so far relevant) 
    !!
    !!  
    tile%soil%theta  = ref_theta
    veg%lai%ptr   = ref_lai 
    IF(ref_pft==0) THEN
      veg%pool%leaf%nitrogen%ptr     = ref_nleaf / molar_mass_N * veg%lai%ptr
    ELSE
      veg%pool%leaf%nitrogen%ptr     = veg%lai%ptr / lctlib%sla / lctlib%cn_leaf
      veg%pool%leaf%phosphorus%ptr   = veg%pool%leaf%nitrogen%ptr / lctlib%np_leaf
    ENDIF
    CALL calc_canopy_layers(lctlib = lctlib, & 
                          leaf_nitrogen = veg%pool%leaf%nitrogen%ptr, &
                          lai = veg%lai%ptr,&
                          lai_cl = veg%lai_cl%ptr(:,1), &
                          cumm_lai_cl = veg%cumm_lai_cl%ptr(:,1), &
                          n_tot_cl = veg%nl_tot_cl%ptr(:,1), &
                          fn_rub_cl = veg%fnl_rub_cl%ptr(:,1), &
                          fn_et_cl = veg%fnl_et_cl%ptr(:,1), &
                          fn_pepc_cl = veg%fnl_pepc_cl%ptr(:,1), &
                          fn_chl_cl = veg%fnl_chl_cl%ptr(:,1), &
                          fn_oth_cl = veg%fnl_oth_cl%ptr(:,1), &
                          inquire_n_fractions=.TRUE.)
    veg%t_jmax_opt_mavg%ptr = 20._wp+0.35_wp*(t_air - Tzero)
  !> 2.0 Sensitivity plots
  !!
  !!  

  !> 2.1 PAR (0-1500 µmol m-2 s-1)
  !! 
  !! 
  lun = 10
  WRITE(lun,*) 'PAR A GC CI' 

  swpar_srf_down = 0.0_wp 
  DO i = 1,16

     swpar_srf_down = swpar_srf_down + 100.0_wp / rad2par 
     swvis_srf_down = swpar_srf_down / frac_par_swvis_srf 
     swnir_srf_down = swvis_srf_down * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf
     CALL update_atm2land(tile, &
           swvis_srf_down = swvis_srf_down, &
           swnir_srf_down = swnir_srf_down, &
           swpar_srf_down = swpar_srf_down)

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    ! mo_assimi_update_canopy_fluxes
     CALL update_radiation(tile)              !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile)          !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') swpar_srf_down * rad2par, & 
           veg%an%ptr,veg%gc%ptr*1000._wp, & 
           veg%ci%ptr/co2_mixing_ratio !,&
!           veg%an%ptr,
!           veg%ci%ptr/co2_mixing_ratio,veg%pool%leaf%nitrogen%ptr*molar_mass_N,&
!           veg%an%ptr/(veg%pool%leaf%nitrogen%ptr*molar_mass_N)

  END DO
  swpar_srf_down = ref_swdown 
  swvis_srf_down = swpar_srf_down / frac_par_swvis_srf 
  swnir_srf_down = swvis_srf_down * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf
  CALL update_atm2land(tile, &
        swvis_srf_down = swvis_srf_down, &
        swnir_srf_down = swnir_srf_down, &
        swpar_srf_down = swpar_srf_down)

  !> 2.2 Temperature (-5,60 °C)
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'T A GC CI' 
  t_air = - 5.0_wp + Tzero
  hlp1 = veg%t_jmax_opt_mavg%ptr
  DO i = 1,70

     t_air = t_air + 1.0_wp 
     q_air = 0.8_wp * calc_spec_humidity_sat(t_air,press_srf)
     veg%t_jmax_opt_mavg%ptr = 20._wp+0.35_wp*(t_air - Tzero)
     CALL update_atm2land(tile, &
           t_air = t_air, &
           q_air = q_air)

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    ! mo_assimi_update_canopy_fluxes
     CALL update_radiation(tile)              !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile)          !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') t_air, &  
           veg%an%ptr,veg%gc%ptr*1000._wp, &
           veg%ci%ptr/co2_mixing_ratio 

  END DO
  t_air    = ref_t_air + Tzero
  q_air    = calc_spec_humidity_sat(t_air,press_srf) - ref_vpd / press_srf * eps_vpd
  CALL update_atm2land(tile, &
        t_air = t_air, &
        q_air = q_air)
  veg%t_jmax_opt_mavg%ptr = hlp1 

  !> 2.3 atmospheric CO2 (175 - 1000 ppm) 
  !! 
  !! 
  lun = lun + 1 
  WRITE(lun,*) 'CO2 A GC CI' 
  co2_mixing_ratio = 150._wp
  DO i = 1,34

     co2_mixing_ratio = co2_mixing_ratio + 25._wp
     CALL update_atm2land(tile,co2_mixing_ratio = co2_mixing_ratio)

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    ! mo_assimi_update_canopy_fluxes
     CALL update_radiation(tile)              !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile)          !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') co2_mixing_ratio, &
           veg%an%ptr,veg%gc%ptr*1000._wp, &
           veg%ci%ptr/co2_mixing_ratio 

  END DO
  co2_mixing_ratio = ref_co2_mixing_ratio
  CALL update_atm2land(tile,co2_mixing_ratio = co2_mixing_ratio)

  !> 2.4 VPD (0-3.3 kPa) 
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'VPD A GC CI' 
  DO i = 1,30
     q_air = calc_spec_humidity_sat(t_air,press_srf) - i * 100.0_wp / press_srf * eps_vpd
     CALL update_atm2land(tile, q_air = q_air) 

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    ! mo_assimi_update_canopy_fluxes
     CALL update_radiation(tile)              !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile)          !< mo_assimi_update_canopy_fluxes

     vpd = ( calc_spec_humidity_sat(t_air,press_srf) - q_air ) / eps_vpd * press_srf / 1000.0_wp

     WRITE(lun,'(4F16.7)') vpd, &
           veg%an%ptr,veg%gc%ptr*1000._wp, &
           veg%ci%ptr/co2_mixing_ratio 

  END DO
  q_air = calc_spec_humidity_sat(t_air,press_srf) - ref_vpd / press_srf * eps_vpd
  vpd = ( calc_spec_humidity_sat(t_air,press_srf) - q_air ) / eps_vpd * press_srf / 1000.0_wp
  CALL update_atm2land(tile, q_air = q_air) 

  !> 2.5 Theta (0-100 %)
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'THETA A GC CI' 
  tile%soil%theta = 0.0_wp
  DO i = 1,100

     tile%soil%theta = tile%soil%theta + 1.0_wp 

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    ! mo_assimi_update_canopy_fluxes
     CALL update_radiation(tile)              !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile)          !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') tile%soil%theta, &
           veg%an%ptr,veg%gc%ptr*1000._wp, &
           veg%ci%ptr/co2_mixing_ratio 

  END DO
  tile%soil%theta = ref_theta

  !> 2.6 LAI (0-11)
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'LAI A GC CI'
  veg%lai%ptr = 0.0_wp
! @todo test_canopy_model -> over what are you looping here ? "i" does not change any as far as I can see (Jan 20160927)
  DO i = 1,22 

     veg%lai%ptr = veg%lai%ptr + 0.5_wp 
     veg%pool%leaf%nitrogen%ptr=ref_nleaf/molar_mass_N*veg%lai%ptr
     CALL calc_canopy_layers(lctlib = lctlib, & 
                          leaf_nitrogen = veg%pool%leaf%nitrogen%ptr, &
                          lai = veg%lai%ptr,&
                          lai_cl = veg%lai_cl%ptr(:,1), &
                          cumm_lai_cl = veg%cumm_lai_cl%ptr(:,1), &
                          n_tot_cl = veg%nl_tot_cl%ptr(:,1), &
                          fn_rub_cl = veg%fnl_rub_cl%ptr(:,1), &
                          fn_et_cl = veg%fnl_et_cl%ptr(:,1), &
                          fn_pepc_cl = veg%fnl_pepc_cl%ptr(:,1), &
                          fn_chl_cl = veg%fnl_chl_cl%ptr(:,1), &
                          fn_oth_cl = veg%fnl_oth_cl%ptr(:,1), &
                          inquire_n_fractions=.TRUE.)

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    ! mo_assimi_update_canopy_fluxes
     CALL update_radiation(tile)              !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile)          !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') veg%lai%ptr, &
           veg%an%ptr,veg%gc%ptr*1000._wp, &
           veg%ci%ptr/co2_mixing_ratio 

  END DO
  veg%lai%ptr = ref_lai 

  !> 2.7 leaf N (0.5-3 g N / m2(LAI))
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'Nleaf A GC CI'
  DO i = 1,24

     veg%pool%leaf%nitrogen%ptr=0.125_wp*(i)/molar_mass_N*veg%lai%ptr
     CALL calc_canopy_layers(lctlib = lctlib, & 
                          leaf_nitrogen = veg%pool%leaf%nitrogen%ptr, &
                          lai = veg%lai%ptr,&
                          lai_cl = veg%lai_cl%ptr(:,1), &
                          cumm_lai_cl = veg%cumm_lai_cl%ptr(:,1), &
                          n_tot_cl = veg%nl_tot_cl%ptr(:,1), &
                          fn_rub_cl = veg%fnl_rub_cl%ptr(:,1), &
                          fn_et_cl = veg%fnl_et_cl%ptr(:,1), &
                          fn_pepc_cl = veg%fnl_pepc_cl%ptr(:,1), &
                          fn_chl_cl = veg%fnl_chl_cl%ptr(:,1), &
                          fn_oth_cl = veg%fnl_oth_cl%ptr(:,1), &
                          inquire_n_fractions=.TRUE.)

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    ! mo_assimi_update_canopy_fluxes
     CALL update_radiation(tile)              !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile)          !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') veg%pool%leaf%nitrogen%ptr/veg%lai%ptr*molar_mass_N, & 
           veg%an%ptr,veg%gc%ptr*1000._wp, &
           veg%ci%ptr/co2_mixing_ratio 

  END DO
  veg%pool%leaf%nitrogen%ptr=ref_nleaf/molar_mass_N*veg%lai%ptr

  !> 2.8 PFT
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'PFT A GC CI'
  DO i = 1,npft
     veg%pool%leaf%nitrogen%ptr  =veg%lai%ptr /model%lctlib(i)%sla /model%lctlib(i)%cn_leaf
     veg%pool%leaf%phosphorus%ptr=veg%pool%leaf%nitrogen%ptr /model%lctlib(i)%np_leaf
     tile%lcts%lib_id = i
     CALL calc_canopy_layers(lctlib = model%lctlib(i), & 
                          leaf_nitrogen = veg%pool%leaf%nitrogen%ptr, &
                          lai = veg%lai%ptr,&
                          lai_cl = veg%lai_cl%ptr(:,1), &
                          cumm_lai_cl = veg%cumm_lai_cl%ptr(:,1), &
                          n_tot_cl = veg%nl_tot_cl%ptr(:,1), &
                          fn_rub_cl = veg%fnl_rub_cl%ptr(:,1), &
                          fn_et_cl = veg%fnl_et_cl%ptr(:,1), &
                          fn_pepc_cl = veg%fnl_pepc_cl%ptr(:,1), &
                          fn_chl_cl = veg%fnl_chl_cl%ptr(:,1), &
                          fn_oth_cl = veg%fnl_oth_cl%ptr(:,1), &
                          inquire_n_fractions=.TRUE.)

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    ! mo_assimi_update_canopy_fluxes
     CALL update_radiation(tile)              !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile)          !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(I1,3F16.7)') i, &
           veg%an%ptr,veg%gc%ptr*1000._wp, &
           veg%ci%ptr/co2_mixing_ratio 

  END DO 
  !STOP
  END SUBROUTINE run_test_canopy_model

END MODULE mo_land_model_interface
