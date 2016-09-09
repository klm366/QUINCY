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
                                  model, &               ! mode & model%lctlib
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

    USE mo_kind
    USE mo_base_memory,                   ONLY: t_slm_tile, t_model
    USE mo_base_config_class,             ONLY: base_config
    USE mo_shared_update_time_averages,   ONLY: update_time_averages  ! includes daytime and moving average
    USE mo_atmland_interface,             ONLY: update_atm2land, update_land2atm
    USE mo_veg_output,                    ONLY: write_output_pool, write_output_flux

    ! used only for the canopy model:
    USE mo_atmland_constants,             ONLY: rad2par
    USE mo_base_constants,                ONLY: one_year

    ! Process interfaces:
    USE mo_rad_interface
    USE mo_veg_interface
    USE mo_assimi_interface
    !USE mo_soil_interface


    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)      :: tile              !< data structure for dynamic vegetation and soil variables
    TYPE(t_model),INTENT(inout)         :: model             !< model structure
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
    !! ==> integrate into vegetation process
      CALL tile%plant%set_zero_fluxes              ! mo_veg_memory

    !! run Tasks
    !!
      CALL update_radiation       (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_vegsoil_bgc     (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_phenology       (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_canopy_fluxes   (tile, model%lctlib(tile%lcts%lib_id))
      CALL calc_growth            (tile, model%lctlib(tile%lcts%lib_id))
      CALL calc_turnover          (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_veg_dynamics    (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_vegsoil_pool    (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_vegsoil_physics (tile)
      CALL update_veg_pools       (tile, model%lctlib(tile%lcts%lib_id))

    !! accumulation over time - daytime & moving average
    !! 
      CALL update_time_averages(tile, model%lctlib(tile%lcts%lib_id), time_day, timestep_prop_day)    ! see shared

    !! output pool: write output in certain intervals, see mo_veg_output
    !! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
    !! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
    !! 
      CALL write_output_pool(tile, model%lctlib(tile%lcts%lib_id), time_day, time_year, 1) ! 1 == itile

    !! output flux: write output in certain intervals, see mo_veg_output
    !! 
    !! 
      CALL write_output_flux(tile, time_day, time_year, 1) ! 1 == itile 

    CASE('plant')
    !! plant --------------------------------------------------------------------------------------
    ! currently same as land !
    !! set vegetation fluxes to zero
    !! 
    !! ==> integrate into vegetation process
      CALL tile%plant%set_zero_fluxes              ! mo_veg_memory
    !! run Tasks
    !!
      CALL update_radiation       (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_vegsoil_bgc     (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_phenology       (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_canopy_fluxes   (tile, model%lctlib(tile%lcts%lib_id))
      CALL calc_growth            (tile, model%lctlib(tile%lcts%lib_id))
      CALL calc_turnover          (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_veg_dynamics    (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_vegsoil_pool    (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_vegsoil_physics (tile)
      CALL update_veg_pools       (tile, model%lctlib(tile%lcts%lib_id))

    !! accumulation over time - daytime & moving average
    !! 
      CALL update_time_averages(tile,  model%lctlib(tile%lcts%lib_id), time_day, timestep_prop_day)    ! see shared

    !! output pool: write output in certain intervals, see mo_veg_output
    !! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
    !! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
    !! 
      CALL write_output_pool(tile, model%lctlib(tile%lcts%lib_id), time_day, time_year, 1) ! 1 == itile

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
    !! ==> integrate into vegetation process
      CALL tile%plant%set_zero_fluxes              ! mo_veg_memory

    !! run Tasks
    !!
      CALL update_veg_pools       (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_radiation       (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_canopy_fluxes   (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_vegsoil_physics (tile)

    !! accumulation over time - daytime & moving average
    !! 
      CALL update_time_averages(tile, model%lctlib(tile%lcts%lib_id), time_day, timestep_prop_day)    ! see shared

    !! output fluxes -> filename: fort.20
    !! output also some atm2land variables
    !! 
      WRITE(20,'(7F16.7)') time_day+time_year*one_year, t_air, swpar_srf_down * rad2par, & 
                           tile%plant%gpp, tile%plant%apar, tile%plant%gc, tile%plant%ci
      WRITE(21,'(20F16.7)') tile%plant%par_sunlit_cl
      WRITE(22,'(20F16.7)') tile%plant%par_shaded_cl
      WRITE(23,'(20F16.7)') tile%plant%fleaf_sunlit_cl
      WRITE(24,'(20F16.7)') tile%plant%par_sunlit_cl*tile%plant%fleaf_sunlit_cl+&
                                tile%plant%par_shaded_cl+(1._wp-tile%plant%fleaf_sunlit_cl)
                           

    CASE('test_canopy')
    !! test_canopy --------------------------------------------------------------------------------
      CALL run_test_canopy_model  (tile, model)

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
  !> Run the test_canopy model from the land model interface
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE run_test_canopy_model(tile, model)

    USE mo_base_memory,                     ONLY: t_slm_tile, t_model
    USE mo_base_constants,                  ONLY: npft,ncanop
    USE mo_jsb_lctlib_class,                ONLY: t_lctlib_element
    USE mo_atmland_constants,               ONLY: rad2par, frac_vis_swtot_srf, frac_par_swvis_srf,eps_vpd, standard_press_srf, &
                                                  def_co2_mixing_ratio
    USE mo_atmland_util,                    ONLY: calc_spec_humidity_sat
    USE mo_shared_constants,                ONLY: molar_mass_N, Tzero 
    ! Process interfaces:
    USE mo_atmland_interface                ! update_atmland
    USE mo_rad_interface                    ! update_radiation
    USE mo_assimi_interface                 ! update_canopy_fluxes
    USE mo_veg_canopy
    ! get parameters defined in init_slm_test_canopy()

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_model),INTENT(inout),TARGET           :: model      !< data structure for dynamic vegetation and soil variables
    TYPE(t_slm_tile),INTENT(inout)               :: tile       !< data structure for dynamic vegetation and soil variables
    ! 0.2 Local
    TYPE(t_lctlib_element),POINTER               :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp) :: vpd
    REAL(wp) :: hlp1
    INTEGER :: i,lun
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
    INTEGER ,SAVE :: ref_pft = 1                      ! PFT identifier

    !  bi-directional pointer
    lctlib => model%lctlib(ref_pft)

    !> 1.0 Init
    !! 
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
    tile%plant%lai   = ref_lai 
    IF(ref_pft==0) THEN
      tile%plant%pool%leaf%nitrogen     = ref_nleaf / molar_mass_N * tile%plant%lai
    ELSE
      tile%plant%pool%leaf%nitrogen     = tile%plant%lai / lctlib%sla / lctlib%cn_leaf
      tile%plant%pool%leaf%phosphorus   = tile%plant%pool%leaf%nitrogen / lctlib%np_leaf
      tile%plant%pftid                  = ref_pft
    ENDIF
    CALL calc_canopy_layers(lctlib = lctlib, & 
                          leaf_nitrogen = tile%plant%pool%leaf%nitrogen, &
                          lai = tile%plant%lai,&
                          lai_cl = tile%plant%lai_cl, &
                          cumm_lai_cl = tile%plant%cumm_lai_cl, &
                          n_tot_cl = tile%plant%nl_tot_cl, &
                          fn_rub_cl = tile%plant%fnl_rub_cl, &
                          fn_et_cl = tile%plant%fnl_et_cl, &
                          fn_pepc_cl = tile%plant%fnl_pepc_cl, &
                          fn_chl_cl = tile%plant%fnl_chl_cl, &
                          fn_oth_cl = tile%plant%fnl_oth_cl, &
                          inquire_n_fractions=.TRUE.)
    tile%plant%t_jmax_opt_mavg = 20._wp+0.35_wp*(t_air - Tzero)
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
     CALL tile%plant%set_zero_fluxes
     CALL update_radiation(tile, lctlib)            !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, lctlib)        !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') swpar_srf_down * rad2par, & 
           tile%plant%an,tile%plant%gc*1000._wp, & 
           tile%plant%ci/co2_mixing_ratio !,&
!           tile%plant%an,
!           tile%plant%ci/co2_mixing_ratio,tile%plant%pool%leaf%nitrogen*molar_mass_N,&
!           tile%plant%an/(tile%plant%pool%leaf%nitrogen*molar_mass_N)

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
  hlp1 = tile%plant%t_jmax_opt_mavg
  DO i = 1,70

     t_air = t_air + 1.0_wp 
     q_air = 0.8_wp * calc_spec_humidity_sat(t_air,press_srf)
     tile%plant%t_jmax_opt_mavg = 20._wp+0.35_wp*(t_air - Tzero)
     CALL update_atm2land(tile, &
           t_air = t_air, &
           q_air = q_air)

     !! canopy calculation fluxes
     CALL tile%plant%set_zero_fluxes
     CALL update_radiation(tile, lctlib)            !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, lctlib)        !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') t_air, &  
           tile%plant%an,tile%plant%gc*1000._wp, &
           tile%plant%ci/co2_mixing_ratio 

  END DO
  t_air    = ref_t_air + Tzero
  q_air    = calc_spec_humidity_sat(t_air,press_srf) - ref_vpd / press_srf * eps_vpd
  CALL update_atm2land(tile, &
        t_air = t_air, &
        q_air = q_air)
  tile%plant%t_jmax_opt_mavg = hlp1 

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
     CALL tile%plant%set_zero_fluxes
     CALL update_radiation(tile, lctlib)            !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, lctlib)        !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') co2_mixing_ratio, &
           tile%plant%an,tile%plant%gc*1000._wp, &
           tile%plant%ci/co2_mixing_ratio 

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
     CALL tile%plant%set_zero_fluxes
     CALL update_radiation(tile, lctlib)            !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, lctlib)        !< mo_assimi_update_canopy_fluxes

     vpd = ( calc_spec_humidity_sat(t_air,press_srf) - q_air ) / eps_vpd * press_srf / 1000.0_wp

     WRITE(lun,'(4F16.7)') vpd, &
           tile%plant%an,tile%plant%gc*1000._wp, &
           tile%plant%ci/co2_mixing_ratio 

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
     CALL tile%plant%set_zero_fluxes
     CALL update_radiation(tile, lctlib)            !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, lctlib)        !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') tile%soil%theta, &
           tile%plant%an,tile%plant%gc*1000._wp, &
           tile%plant%ci/co2_mixing_ratio 

  END DO
  tile%soil%theta = ref_theta

  !> 2.6 LAI (0-11)
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'LAI A GC CI'
  tile%plant%lai = 0.0_wp
  DO i = 1,22

     tile%plant%lai = tile%plant%lai + 0.5_wp 
     tile%plant%pool%leaf%nitrogen=ref_nleaf/molar_mass_N*tile%plant%lai
     CALL calc_canopy_layers(lctlib = lctlib, & 
                          leaf_nitrogen = tile%plant%pool%leaf%nitrogen, &
                          lai = tile%plant%lai,&
                          lai_cl = tile%plant%lai_cl, &
                          cumm_lai_cl = tile%plant%cumm_lai_cl, &
                          n_tot_cl = tile%plant%nl_tot_cl, &
                          fn_rub_cl = tile%plant%fnl_rub_cl, &
                          fn_et_cl = tile%plant%fnl_et_cl, &
                          fn_pepc_cl = tile%plant%fnl_pepc_cl, &
                          fn_chl_cl = tile%plant%fnl_chl_cl, &
                          fn_oth_cl = tile%plant%fnl_oth_cl, &
                          inquire_n_fractions=.TRUE.)

     !! canopy calculation fluxes
     CALL tile%plant%set_zero_fluxes
     CALL update_radiation(tile, lctlib)            !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, lctlib)        !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') tile%plant%lai, &
           tile%plant%an,tile%plant%gc*1000._wp, &
           tile%plant%ci/co2_mixing_ratio 

  END DO
  tile%plant%lai = ref_lai 

  !> 2.7 leaf N (0.5-3 g N / m2(LAI))
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'Nleaf A GC CI'
  DO i = 1,24

     tile%plant%pool%leaf%nitrogen=0.125_wp*(i)/molar_mass_N*tile%plant%lai
     CALL calc_canopy_layers(lctlib = lctlib, & 
                          leaf_nitrogen = tile%plant%pool%leaf%nitrogen, &
                          lai = tile%plant%lai,&
                          lai_cl = tile%plant%lai_cl, &
                          cumm_lai_cl = tile%plant%cumm_lai_cl, &
                          n_tot_cl = tile%plant%nl_tot_cl, &
                          fn_rub_cl = tile%plant%fnl_rub_cl, &
                          fn_et_cl = tile%plant%fnl_et_cl, &
                          fn_pepc_cl = tile%plant%fnl_pepc_cl, &
                          fn_chl_cl = tile%plant%fnl_chl_cl, &
                          fn_oth_cl = tile%plant%fnl_oth_cl, &
                          inquire_n_fractions=.TRUE.)

     !! canopy calculation fluxes
     CALL tile%plant%set_zero_fluxes
     CALL update_radiation(tile, lctlib)            !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, lctlib)        !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') tile%plant%pool%leaf%nitrogen/tile%plant%lai*molar_mass_N, & 
           tile%plant%an,tile%plant%gc*1000._wp, &
           tile%plant%ci/co2_mixing_ratio 

  END DO
  tile%plant%pool%leaf%nitrogen=ref_nleaf/molar_mass_N*tile%plant%lai

  !> 2.8 PFT
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'PFT A GC CI'
  DO i = 1,npft
     tile%plant%pool%leaf%nitrogen  =tile%plant%lai /model%lctlib(i)%sla /model%lctlib(i)%cn_leaf
     tile%plant%pool%leaf%phosphorus=tile%plant%pool%leaf%nitrogen /model%lctlib(i)%np_leaf
     tile%plant%pftid=i
     CALL calc_canopy_layers(lctlib = model%lctlib(i), & 
                          leaf_nitrogen = tile%plant%pool%leaf%nitrogen, &
                          lai = tile%plant%lai,&
                          lai_cl = tile%plant%lai_cl, &
                          cumm_lai_cl = tile%plant%cumm_lai_cl, &
                          n_tot_cl = tile%plant%nl_tot_cl, &
                          fn_rub_cl = tile%plant%fnl_rub_cl, &
                          fn_et_cl = tile%plant%fnl_et_cl, &
                          fn_pepc_cl = tile%plant%fnl_pepc_cl, &
                          fn_chl_cl = tile%plant%fnl_chl_cl, &
                          fn_oth_cl = tile%plant%fnl_oth_cl, &
                          inquire_n_fractions=.TRUE.)

     !! canopy calculation fluxes
     CALL tile%plant%set_zero_fluxes
     CALL update_radiation(tile, model%lctlib(i))            !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, model%lctlib(i))        !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(I1,3F16.7)') tile%plant%pftid, &
           tile%plant%an,tile%plant%gc*1000._wp, &
           tile%plant%ci/co2_mixing_ratio 

  END DO 
  !STOP
  END SUBROUTINE run_test_canopy_model

END MODULE mo_land_model_interface
