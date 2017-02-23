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

  CHARACTER(len=*), PARAMETER :: modname = 'mo_land_model_interface'

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


    USE mo_jsb_class,                     ONLY: Get_model
    USE mo_base_memory,                   ONLY: t_slm_tile, t_model
    USE mo_base_config_class,             ONLY: base_config
    USE mo_shared_update_time_averages,   ONLY: update_time_averages  ! includes daytime and moving average
    USE mo_atmland_interface,             ONLY: update_atm2land, update_land2atm
    USE mo_shared_output,                 ONLY: write_output_pool, write_output_flux, write_soil_forcing
    USE mo_veg_util,                      ONLY: set_zero_vegetation_fluxes

    ! Process interfaces:
    USE mo_rad_interface
    USE mo_veg_interface
    USE mo_assimi_interface
    USE mo_soil_interface
    USE mo_sb_interface
    USE mo_pheno_interface


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)      :: tile              !< data structure for dynamic vegetation and soil variables
    LOGICAL,INTENT(in)                  :: l_forcing         !< IF .TRUE. THEN update_forcing(), update_atm2land(), update_land2atm()
    REAL(wp),INTENT(in)                 :: time_day          !< displays the 'real time' in days; time_day == 1.0 after one day has passed
    INTEGER,INTENT(in)                  :: time_year         !< displays the 'real time' in years; the 1st year of a simulation == 1 (not zero)
    REAL(wp),INTENT(in)                 :: timestep_prop_day !< proportion of one day that is simulated with one timestep (dtime/one_day)
    REAL(wp),INTENT(in)           :: &               !! atm2land memory - INTENT(in) - start
                                  t_air(:)              , &
                                  q_air(:)              , &
                                  rain(:)               , &
                                  snow(:)               , &
                                  wind_air(:)           , &
                                  wind_10m(:)           , &
                                  lw_srf_down(:)        , &
                                  swvis_srf_down(:)     , &
                                  swnir_srf_down(:)     , &
                                  swpar_srf_down(:)     , &
                                  swpar_fdiffuse(:)     , &
                                  press_srf(:)          , &
                                  drag_srf(:)           , &
                                  t_acoef(:)            , &
                                  t_bcoef(:)            , &
                                  q_acoef(:)            , &
                                  q_bcoef(:)            , &
                                  pch(:)                , &
                                  cos_zenith_angle(:)   , &
                                  co2_mixing_ratio(:)     , &
                                  co2_mixing_ratio_C13(:) , &
                                  co2_mixing_ratio_C14(:) , &
                                  ga(:)                   , &
                                  ! For lakes:
                                  albedo_lwtr(:)        , &
                                  drag_wtr(:)           , &
                                  drag_ice(:)           , &
                                  t_acoef_wtr(:)        , &
                                  t_bcoef_wtr(:)        , &
                                  q_acoef_wtr(:)        , &
                                  q_bcoef_wtr(:)        , &
                                  t_acoef_ice(:)        , &
                                  t_bcoef_ice(:)        , &
                                  q_acoef_ice(:)        , &
                                  q_bcoef_ice(:)        
                                                     !! atm2land memory - end
    REAL(wp),INTENT(out)          :: &              !! land2atm memory - INTENT(out) - start
                                  t_srf(:)             ,   &
                                  t_eff_srf(:)         ,   &
                                  qsat_srf(:)          ,   &
                                  fact_q_air(:)        ,   &
                                  fact_qsat_srf(:)     ,   &
                                  zh_srf(:)            ,   &
                                  zm_srf(:)            ,   &
                                  alb_vis_dif(:)       ,   &
                                  alb_nir_dif(:)       ,   &
                                  alb_vis_dir(:)       ,   &
                                  alb_nir_dir(:)       ,   &
                                  albedo(:)            ,   &
                                  sensible_heat_flx(:) ,   &
                                  latent_heat_flx(:)   ,   &
                                  net_radiation(:)     ,   &
                                  vpd(:)               ,   &
                                  ! For lakes:
                                  latent_hflx_wtr(:)   ,   &
                                  sensible_hflx_wtr(:) ,   &
                                  latent_hflx_ice(:)   ,   &
                                  sensible_hflx_ice(:) !! land2atm memory - end
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_model),          POINTER     :: model
    CHARACTER(len=*),       PARAMETER   :: routine = TRIM(modname)//':land_model_interface'
    ! ---------------------------
    ! 0.5 Get Memory
    ! 
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> pass atmland variables from forcing input (driver) to atm2land
    !! 
    IF (l_forcing) THEN       ! l_forcing is true
      CALL update_atm2land( &
      ! the sequence of input variables is consistent with the sequence of var in t_atm2land_memory & atmland_interface
      tile                                         , &
      model%options                                , &
      t_air                   = t_air(:)              , &
      q_air                   = q_air(:)              , &
      rain                    = rain(:)               , &
      snow                    = snow(:)               , &
      wind_air                = wind_air(:)           , &
      wind_10m                = wind_10m(:)           , &
      lw_srf_down             = lw_srf_down(:)        , &
      swvis_srf_down          = swvis_srf_down(:)     , &
      swnir_srf_down          = swnir_srf_down(:)     , &
      swpar_srf_down          = swpar_srf_down(:)     , &
      swpar_fdiffuse          = swpar_fdiffuse(:)     , &
      press_srf               = press_srf(:)          , &
      drag_srf                = drag_srf(:)           , &
      t_acoef                 = t_acoef(:)            , &
      t_bcoef                 = t_bcoef(:)            , &
      q_acoef                 = q_acoef(:)            , &
      q_bcoef                 = q_bcoef(:)            , &
      pch                     = pch(:)                , &
      cos_zenith_angle        = cos_zenith_angle(:)   , &
      co2_mixing_ratio        = co2_mixing_ratio(:)     , &
      co2_mixing_ratio_C13    = co2_mixing_ratio_C13(:) , &
      co2_mixing_ratio_C14    = co2_mixing_ratio_C14(:) , &
      ga                      = ga(:)                   , &
      ! For lakes:
      albedo_lwtr             = albedo_lwtr(:)        , &
      drag_wtr                = drag_wtr(:)           , &
      drag_ice                = drag_ice(:)           , &
      t_acoef_wtr             = t_acoef_wtr(:)        , &
      t_bcoef_wtr             = t_bcoef_wtr(:)        , &
      q_acoef_wtr             = q_acoef_wtr(:)        , &
      q_bcoef_wtr             = q_bcoef_wtr(:)        , &
      t_acoef_ice             = t_acoef_ice(:)        , &
      t_bcoef_ice             = t_bcoef_ice(:)        , &
      q_acoef_ice             = q_acoef_ice(:)        , &
      q_bcoef_ice             = q_bcoef_ice(:)          &
      & )
    END IF


    !> @todo mo_land_model_interface:land_model_interface put all Tasks into a Queue accordingly to the model selection
    !! (land plant soil canopy test_canopy)
    SELECT CASE(TRIM(base_config%slm_model_name))
    !! land ---------------------------------------------------------------------------------------
    CASE('land')
    !! set vegetation fluxes to zero
    !! 
      CALL set_zero_vegetation_fluxes(tile)   ! mo_veg_util

    !! run Tasks
    !!
      CALL update_radiation       (tile, model%options)
      CALL update_vegsoil_bgc     (tile, model%options)
      CALL update_phenology       (tile, model%options)
      CALL update_canopy_fluxes   (tile, model%options)
      CALL update_growth          (tile, model%options)
      CALL update_veg_turnover    (tile, model%options)
      CALL update_veg_dynamics    (tile, model%options)
      CALL update_sb_pools        (tile, model%options)
      CALL update_soil_physics    (tile, model%options)
      CALL update_veg_pools       (tile, model%options)

    !! accumulation over time - daytime & moving average
    !! 
      CALL update_time_averages(tile, model%options, time_day, timestep_prop_day)    ! see shared

    !! output pool: write output in certain intervals, see mo_shared_output
    !! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
    !! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
    !! 
      CALL write_output_pool(tile, model%options, time_day, time_year, 1) ! 1 == itile

    !! output flux: write output in certain intervals, see mo_shared_output
    !! 
    !! 
      CALL write_output_flux(tile, model%options, time_day, time_year, 1) ! 1 == itile 

    !! output litterfall and soilphysics: usually this would be output at each timstep
    !! 
    !! 
      IF (base_config%flag_generate_soil_forcing) THEN
        CALL write_soil_forcing(tile, model%options, time_day, time_year, 1) ! 1 == itile 
      END IF

    CASE('plant')
    !! plant --------------------------------------------------------------------------------------
    ! currently same as land !
    !! set vegetation fluxes to zero
    !! 
      CALL set_zero_vegetation_fluxes(tile)   ! mo_veg_util

    !! run Tasks
    !!
      CALL update_radiation       (tile, model%options)
      CALL update_vegsoil_bgc     (tile, model%options)
      CALL update_phenology       (tile, model%options)
      CALL update_canopy_fluxes   (tile, model%options)
      CALL update_growth          (tile, model%options)
      CALL update_veg_turnover    (tile, model%options)
      CALL update_veg_dynamics    (tile, model%options)
      CALL update_sb_pools        (tile, model%options)
      CALL update_soil_physics    (tile, model%options)
      CALL update_veg_pools       (tile, model%options)

    !! accumulation over time - daytime & moving average
    !! 
      CALL update_time_averages(tile, model%options, time_day, timestep_prop_day)    ! see shared

    !! output pool: write output in certain intervals, see mo_shared_output
    !! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
    !! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
    !! 
      CALL write_output_pool(tile, model%options, time_day, time_year, 1) ! 1 == itile

    !! output flux: write output in certain intervals, see mo_shared_output
    !! 
    !! 
      CALL write_output_flux(tile, model%options, time_day, time_year, 1) ! 1 == itile 

    !! output litterfall and soilphysics: usually this would be output at each timstep
    !! 
    !! 
      IF (base_config%flag_generate_soil_forcing) THEN
        CALL write_soil_forcing(tile, model%options, time_day, time_year, 1) ! 1 == itile 
      END IF

    CASE('soil')
    !! soil ---------------------------------------------------------------------------------------

      CALL update_sb_pools        (tile, model%options)
      !! output pool: write output in certain intervals, see mo_shared_output
      !! note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
      !! note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
      !! 
      CALL write_output_pool(tile, model%options, time_day, time_year, 1) ! 1 == itile
  
      !! output flux: write output in certain intervals, see mo_shared_output
      !! 
      !! 
      CALL write_output_flux(tile, model%options, time_day, time_year, 1) ! 1 == itile 

    CASE('canopy')
    !! canopy -------------------------------------------------------------------------------------
    !! set vegetation fluxes to zero
    !! 
      CALL set_zero_vegetation_fluxes(tile)   ! mo_veg_util

    !! run Tasks
    !!
      CALL update_veg_pools       (tile, model%options)
      CALL update_radiation       (tile, model%options)
      CALL update_canopy_fluxes   (tile, model%options)
      CALL update_soil_physics    (tile, model%options)

    !! accumulation over time - daytime & moving average
    !! 
      CALL update_time_averages(tile, model%options, time_day, timestep_prop_day)    ! see shared

    !! output fluxes -> filename: fort.20
    !! output also some atm2land variables
    !! 
      CALL canopy_model_write_data(tile, model%options, time_day, time_year, t_air(1), swpar_srf_down(1))                           

    CASE('test_canopy')
    !! test_canopy --------------------------------------------------------------------------------
      CALL run_test_canopy_model  (tile, model%options)

    CASE DEFAULT
      CALL finish('land_model_interface','Invalid SLM model name. Try either land, plant, soil, canopy or test_canopy.')
    END SELECT  ! select SLM model !

    !> pass atmland variables from land2atm memory back to the driver
    !! 
    IF (l_forcing) THEN       ! l_forcing is true if update_atm2land() is called (depends on the SLM model selected)
      CALL update_land2atm( &
      ! the sequence of input variables is consistent with the sequence of var in t_land2atm_memory !
      tile                                              ,   &
      model%options                                     ,   &
      t_srf                   = t_srf(:)             ,   &
      t_eff_srf               = t_eff_srf(:)         ,   &
      qsat_srf                = qsat_srf(:)          ,   &
      fact_q_air              = fact_q_air(:)        ,   &
      fact_qsat_srf           = fact_qsat_srf(:)     ,   &
      zh_srf                  = zh_srf(:)            ,   &
      zm_srf                  = zm_srf(:)            ,   &
      alb_vis_dif             = alb_vis_dif(:)       ,   &
      alb_nir_dif             = alb_nir_dif(:)       ,   &
      alb_vis_dir             = alb_vis_dir(:)       ,   &
      alb_nir_dir             = alb_nir_dir(:)       ,   &
      albedo                  = albedo(:)            ,   &
      sensible_heat_flx       = sensible_heat_flx(:) ,   &
      latent_heat_flx         = latent_heat_flx(:)   ,   &
      net_radiation           = net_radiation(:)     ,   &
      vpd                     = vpd(:)               ,   &
      ! For lakes:
      latent_hflx_wtr         = latent_hflx_wtr(:)   ,   &
      sensible_hflx_wtr       = sensible_hflx_wtr(:) ,   &
      latent_hflx_ice         = latent_hflx_ice(:)   ,   &
      sensible_hflx_ice       = sensible_hflx_ice(:)     &
      & )
    END IF
  END SUBROUTINE land_model_interface


  !-----------------------------------------------------------------------------------------------------
  !> Write results and data from canopy model
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE canopy_model_write_data(tile, options, time_day, time_year, t_air, swpar_srf_down)

    USE mo_rad_constants,         ONLY: rad2ppfd
    USE mo_base_constants,        ONLY: one_year
    USE mo_base_memory,           ONLY: t_slm_tile, t_model, t_jsb_task_options
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_

    ! Use of process memories
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),           INTENT(inout)     :: tile         !< one tile
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    REAL(wp),        INTENT(in)         :: time_day          !< displays the 'real time' in days; time_day == 1.0 after one day has passed
    INTEGER,         INTENT(in)         :: time_year         !< displays the 'real time' in years; the 1st year of a simulation == 1 (not zero)
    REAL(wp),        INTENT(in)         :: t_air             !< air temp
    REAL(wp),        INTENT(in)         :: swpar_srf_down    !< short wave light at surface level ?
    ! ---------------------------
    ! 0.2 Local
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':canopy_model_write_data'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ASSIMI_)                       !   TYPE(t_assimi_memory), POINTER   :: assimi__mem  
    dsl4jsb_Def_memory(VEG_)                          !   TYPE(t_veg_memory), POINTER      :: veg__mem  
    dsl4jsb_Def_memory(RAD_)                          !     TYPE(t_rad_memory),      POINTER :: rad__mem
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk      :: gross_assimilation        ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: canopy_cond
    dsl4jsb_Real2D_onChunk      :: co2_conc_leaf
    dsl4jsb_Real2D_onChunk      :: appfd                     ! RAD_
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_cl           ! VEG_
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(xx_)
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(RAD_)
    ! Get process variables (Set pointers to variables in memory)
!    model  => Get_model(tile%owner_model_id)
!    lctlib => model%lctlib(tile%lcts%lib_id) 
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, gross_assimilation)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, canopy_cond)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, co2_conc_leaf)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(RAD_, appfd)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_cl)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_cl)



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


      WRITE(20,'(7F16.7)')  time_day+time_year*one_year, t_air, swpar_srf_down * rad2ppfd, &
                            gross_assimilation(1), appfd(1), & 
                            canopy_cond(1), co2_conc_leaf(1) 
      WRITE(21,'(20F16.7)') ppfd_sunlit_cl(1,:)       ! write values of all canopy layers
      WRITE(22,'(20F16.7)') ppfd_shaded_cl(1,:)       ! write values of all canopy layers
      WRITE(23,'(20F16.7)') fleaf_sunlit_cl(1,:)   ! write values of all canopy layers
      WRITE(24,'(20F16.7)') ppfd_sunlit_cl(1,:) * fleaf_sunlit_cl(1,:) + &
                            & ppfd_shaded_cl(1,:) + (1._wp - fleaf_sunlit_cl(1,:))  ! write values of all canopy layers

  END SUBROUTINE canopy_model_write_data


  !-----------------------------------------------------------------------------------------------------
  !> Run the test_canopy model from the land model interface
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE run_test_canopy_model(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model, t_jsb_task_options
    USE mo_base_constants,        ONLY: npft
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_rad_constants,         ONLY: rad2ppfd
    USE mo_atmland_constants,     ONLY: frac_vis_swtot_srf, frac_par_swvis_srf,eps_vpd, standard_press_srf, def_co2_mixing_ratio
    USE mo_atmland_util,          ONLY: calc_spec_humidity_sat
    USE mo_veg_util,              ONLY: set_zero_vegetation_fluxes
    USE mo_shared_constants,      ONLY: molar_mass_N, Tzero 

    ! Use of process memories
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory

    ! Process interfaces:
    USE mo_atmland_interface                ! update_atmland
    USE mo_rad_interface                    ! update_radiation
    USE mo_assimi_interface                 ! update_canopy_fluxes
    USE mo_veg_canopy
    ! get parameters defined in init_slm_test_canopy()

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),           INTENT(inout)     :: tile         !< one tile
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_model),          POINTER           :: model
    TYPE(t_lctlib_element), POINTER           :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp), DIMENSION(options%nc)           :: & 
                  t_air           , &
                  press_srf       , &
                  q_air           , &
                  co2_mixing_ratio, &
                  rain            , &
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
                  cos_zenith_angle, &
                  ga             , &
                  !latitude      , &
                  !longitude     , &
                  vpd
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
    REAL(wp),SAVE :: ref_theta = 1._wp                ! soil water availability in fraction of AWC 
    REAL(wp),SAVE :: ref_lai = 3._wp                  ! LAI
    REAL(wp),SAVE :: ref_nleaf = 2.0_wp               ! N leaf per unit leaf area in gN/m2
    INTEGER       :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':run_test_canopy_model'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ASSIMI_)                       !   TYPE(t_assimi_memory), POINTER      :: assimi__mem  
    dsl4jsb_Def_memory(VEG_)                          !   TYPE(t_veg_memory), POINTER      :: veg__mem  
    dsl4jsb_Def_memory(SOIL_)                         !   TYPE(t_soil_memory), POINTER      :: soil__mem  
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk      :: net_assimilation          ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: canopy_cond
    dsl4jsb_Real2D_onChunk      :: co2_conc_leaf
    dsl4jsb_Real2D_onChunk      :: w_soil_root_theta         ! SOIL_
    dsl4jsb_Real2D_onChunk      :: lai                       ! VEG_
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    dsl4jsb_Real3D_onChunk      :: lai_cl
    dsl4jsb_Real3D_onChunk      :: cumm_lai_cl
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl

! 2 be removed: pointers to veg assimi *slmtemp
    TYPE(t_veg_memory), POINTER         :: veg        ! local pointer to veg__mem for SLM usability
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(PHENO_)
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SOIL_)
    ! Get process variables (Set pointers to variables in memory)
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
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, net_assimilation)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, canopy_cond)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, co2_conc_leaf)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_theta)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, cumm_lai_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_rub_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_et_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_pepc_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_oth_cl)



! 2 be removed  *slmtemp
    veg    => veg__mem 


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 Init
    !! 

    !> 1.1 Initalisation of met-forcing which is assumed constant 
    !! 
    !!  
    t_air(:)            = ref_t_air + TZero 
    press_srf(:)        = standard_press_srf
    q_air(:)            = calc_spec_humidity_sat(t_air(:),press_srf(:)) - ref_vpd / press_srf(1) * eps_vpd
    co2_mixing_ratio(:) = def_co2_mixing_ratio   ! def_co2_mixing_ratio is used instead of ref_co2_mixing_ratio 
    rain(:)             = 0.0_wp
    snow(:)             = 0.0_wp
    swpar_srf_down(:)   = ref_swdown 
    swvis_srf_down(:)   = swpar_srf_down(:) / frac_par_swvis_srf 
    swnir_srf_down(:)   = swvis_srf_down(:) * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf
    swpar_fdiffuse(:)   = ref_fdiffuse
    wind_air(:)         = ref_wind
    wind_10m(:)         = wind_air(:)
    lw_srf_down(:)      = 0.0_wp
    drag_srf(:)         = 1.0_wp
    t_acoef(:)          = 0.0_wp
    t_bcoef(:)          = 0.0_wp
    q_acoef(:)          = 0.0_wp
    q_bcoef(:)          = 0.0_wp
    pch(:)              = 0.0_wp
    cos_zenith_angle(:) = ref_cos_angle
    ga(:)               = 0.4_wp / ( drag_srf(1) * wind_10m(1) )
    !latitude                  = 50.0_wp
    !longitude                 = 0.0_wp

    CALL update_atm2land( &
           tile                                            , &
           options                                         , &
           t_air                   = t_air(:)              , &
           q_air                   = q_air(:)              , &
           wind_air                = wind_air(:)           , &
           wind_10m                = wind_10m(:)           , &
           swvis_srf_down          = swvis_srf_down(:)     , &
           swnir_srf_down          = swnir_srf_down(:)     , &
           swpar_srf_down          = swpar_srf_down(:)     , &
           swpar_fdiffuse          = swpar_fdiffuse(:)     , &
           press_srf               = press_srf(:)          , &
           drag_srf                = drag_srf(:)           , &
           cos_zenith_angle        = cos_zenith_angle(:)   , &
           co2_mixing_ratio        = co2_mixing_ratio(:)   , &
           ga                      = ga(:)                 )

    !> 1.2 Initialisation of soil and vegetation (so far relevant) 
    !!
    !!  
    w_soil_root_theta(:) = ref_theta
    lai(1)   = ref_lai 
    IF(ref_pft==0) THEN
      veg%pool%leaf%nitrogen%ptr     = ref_nleaf / molar_mass_N * lai(1)
    ELSE
      veg%pool%leaf%nitrogen%ptr     = lai(1) / lctlib%sla / lctlib%cn_leaf
      veg%pool%leaf%phosphorus%ptr   = veg%pool%leaf%nitrogen%ptr / lctlib%np_leaf
    ENDIF
    CALL calc_canopy_layers(lctlib = lctlib, & 
                            nc = nc, &
                          leaf_nitrogen = veg%pool%leaf%nitrogen%ptr, &
                          lai = lai(1),&
                          lai_cl = lai_cl(1,:), &
                          cumm_lai_cl = cumm_lai_cl(1,:), &
                          leaf_nitrogen_cl = leaf_nitrogen_cl(1,:), &
                          fn_rub_cl = fn_rub_cl(1,:), &
                          fn_et_cl = fn_et_cl(1,:), &
                          fn_pepc_cl = fn_pepc_cl(1,:), &
                          fn_chl_cl = fn_chl_cl(1,:), &
                          fn_oth_cl = fn_oth_cl(1,:), &
                          inquire_n_fractions=.TRUE.)
    t_jmax_opt_mavg(1) = 20._wp+0.35_wp*(t_air(1) - Tzero)
  !> 2.0 Sensitivity plots
  !!
  !!  

  !> 2.1 PAR (0-1500 µmol m-2 s-1)
  !! 
  !! 
  lun = 10
  WRITE(lun,*) 'PAR A GC CI' 

  swpar_srf_down(:) = 0.0_wp 
  DO i = 1,16

     swpar_srf_down(:) = swpar_srf_down(:) + 100.0_wp / rad2ppfd 
     swvis_srf_down(:) = swpar_srf_down(:) / frac_par_swvis_srf 
     swnir_srf_down(:) = swvis_srf_down(:) * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf
     CALL update_atm2land(tile, &
           options            , &
           swvis_srf_down = swvis_srf_down(:), &
           swnir_srf_down = swnir_srf_down(:), &
           swpar_srf_down = swpar_srf_down(:))

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') swpar_srf_down(1) * rad2ppfd, & 
           net_assimilation(1), canopy_cond(1) *1000._wp, & 
           co2_conc_leaf(1)/co2_mixing_ratio(1) !,&
!           net_assimilation(1),
!           co2_conc_leaf(1)/co2_mixing_ratio(1),veg%pool%leaf%nitrogen%ptr*molar_mass_N,&
!           net_assimilation(1)/(veg%pool%leaf%nitrogen%ptr*molar_mass_N)

  END DO
  swpar_srf_down(:) = ref_swdown 
  swvis_srf_down(:) = swpar_srf_down(:) / frac_par_swvis_srf 
  swnir_srf_down(:) = swvis_srf_down(:) * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf
  CALL update_atm2land(tile, &
           options            , &
           swvis_srf_down = swvis_srf_down, &
           swnir_srf_down = swnir_srf_down, &
           swpar_srf_down = swpar_srf_down)

  !> 2.2 Temperature (-5,60 °C)
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'T A GC CI' 
  t_air(:) = - 5.0_wp + Tzero
  hlp1 = t_jmax_opt_mavg(1)

  DO i = 1,70

     t_air(:) = t_air(:) + 1.0_wp 
     q_air(:) = 0.8_wp * calc_spec_humidity_sat(t_air(:),press_srf(:))
     t_jmax_opt_mavg(1) = 20._wp+0.35_wp*(t_air(1) - Tzero)
     CALL update_atm2land(tile, &
           options            , &
           t_air = t_air, &
           q_air = q_air)

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(5F16.7)') t_air(1), &  
           net_assimilation(1), canopy_cond(1) *1000._wp, &
           co2_conc_leaf(1)/co2_mixing_ratio(1) 

  END DO
  t_air(:)    = ref_t_air + Tzero
  q_air(:)    = calc_spec_humidity_sat(t_air(:),press_srf(:)) - ref_vpd / press_srf(1) * eps_vpd
  CALL update_atm2land(tile, &
        options            , &
        t_air = t_air, &
        q_air = q_air)
  t_jmax_opt_mavg(1) = hlp1 

  !> 2.3 atmospheric CO2 (175 - 1000 ppm) 
  !! 
  !! 
  lun = lun + 1 
  WRITE(lun,*) 'CO2 A GC CI' 
  co2_mixing_ratio(:) = 150._wp
  DO i = 1,34

     co2_mixing_ratio(:) = co2_mixing_ratio(:) + 25._wp
     CALL update_atm2land(tile, &
           options            , &
           co2_mixing_ratio = co2_mixing_ratio)

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') co2_mixing_ratio(1), &
           net_assimilation(1), canopy_cond(1) *1000._wp, &
           co2_conc_leaf(1)/co2_mixing_ratio(1) 

  END DO
  co2_mixing_ratio(:) = ref_co2_mixing_ratio
  CALL update_atm2land(tile    , &
                       options , &
                       co2_mixing_ratio = co2_mixing_ratio)

  !> 2.4 VPD (0-3.3 kPa) 
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'VPD A GC CI' 
  DO i = 1,30
     q_air(:) = calc_spec_humidity_sat(t_air(:),press_srf(:)) - i * 100.0_wp / press_srf(1) * eps_vpd
     CALL update_atm2land(tile    , &
                          options , &
                          q_air = q_air) 

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     vpd(:) = ( calc_spec_humidity_sat(t_air(:),press_srf(:)) - q_air(:) ) / eps_vpd * press_srf(1) / 1000.0_wp

     WRITE(lun,'(4F16.7)') vpd(1), &
           net_assimilation(1), canopy_cond(1) *1000._wp, &
           co2_conc_leaf(1)/co2_mixing_ratio(1) 

  END DO
  q_air(:) = calc_spec_humidity_sat(t_air(:),press_srf(:)) - ref_vpd / press_srf(1) * eps_vpd
  vpd(:) = ( calc_spec_humidity_sat(t_air(:),press_srf(:)) - q_air(:) ) / eps_vpd * press_srf(1) / 1000.0_wp
     CALL update_atm2land(tile    , &
                          options , &
                          q_air = q_air) 

  !> 2.5 Theta (0-100 %)
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'THETA A GC CI' 
  w_soil_root_theta(:) = 0.0_wp
  DO i = 1,100

     w_soil_root_theta(:) = w_soil_root_theta(:) + 0.01_wp 

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') w_soil_root_theta(1) *100._wp, &
           net_assimilation(1), canopy_cond(1) *1000._wp, &
           co2_conc_leaf(1)/co2_mixing_ratio(1) 

  END DO
  w_soil_root_theta(:) = ref_theta

  !> 2.6 LAI (0-11)
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'LAI A GC CI'
  lai(1) = 0.0_wp
! @todo test_canopy_model -> over what are you looping here ? "i" does not change any as far as I can see (Jan 20160927)
  DO i = 1,22 

     lai(1) = lai(1) + 0.5_wp 
     veg%pool%leaf%nitrogen%ptr=ref_nleaf/molar_mass_N*lai(1)
     CALL calc_canopy_layers(lctlib = lctlib, & 
                            nc = nc, &
                          leaf_nitrogen = veg%pool%leaf%nitrogen%ptr, &
                          lai = lai(1),&
                          lai_cl = lai_cl(1,:), &
                          cumm_lai_cl = cumm_lai_cl(1,:), &
                          leaf_nitrogen_cl = leaf_nitrogen_cl(1,:), &
                          fn_rub_cl = fn_rub_cl(1,:), &
                          fn_et_cl = fn_et_cl(1,:), &
                          fn_pepc_cl = fn_pepc_cl(1,:), &
                          fn_chl_cl = fn_chl_cl(1,:), &
                          fn_oth_cl = fn_oth_cl(1,:), &
                          inquire_n_fractions=.TRUE.)

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') lai(1), &
           net_assimilation(1), canopy_cond(1) *1000._wp, &
           co2_conc_leaf(1)/co2_mixing_ratio(1) 

  END DO
  lai(1) = ref_lai 

  !> 2.7 leaf N (0.5-3 g N / m2(LAI))
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'Nleaf A GC CI'
  DO i = 1,24

     veg%pool%leaf%nitrogen%ptr=0.125_wp*(i)/molar_mass_N*lai(1)
     CALL calc_canopy_layers(lctlib = lctlib, & 
                            nc = nc, &
                          leaf_nitrogen = veg%pool%leaf%nitrogen%ptr, &
                          lai = lai(1),&
                          lai_cl = lai_cl(1,:), &
                          cumm_lai_cl = cumm_lai_cl(1,:), &
                          leaf_nitrogen_cl = leaf_nitrogen_cl(1,:), &
                          fn_rub_cl = fn_rub_cl(1,:), &
                          fn_et_cl = fn_et_cl(1,:), &
                          fn_pepc_cl = fn_pepc_cl(1,:), &
                          fn_chl_cl = fn_chl_cl(1,:), &
                          fn_oth_cl = fn_oth_cl(1,:), &
                          inquire_n_fractions=.TRUE.)

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') veg%pool%leaf%nitrogen%ptr/lai(1)*molar_mass_N, & 
           net_assimilation(1), canopy_cond(1) *1000._wp, &
           co2_conc_leaf(1)/co2_mixing_ratio(1) 

  END DO
  veg%pool%leaf%nitrogen%ptr=ref_nleaf/molar_mass_N*lai(1)

  !> 2.8 PFT
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'PFT A GC CI'
  DO i = 1,npft
     veg%pool%leaf%nitrogen%ptr  =lai(1) /model%lctlib(i)%sla /model%lctlib(i)%cn_leaf
     veg%pool%leaf%phosphorus%ptr=veg%pool%leaf%nitrogen%ptr /model%lctlib(i)%np_leaf
     tile%lcts%lib_id = i
     CALL calc_canopy_layers(lctlib = model%lctlib(i), & 
                            nc = nc, &
                          leaf_nitrogen = veg%pool%leaf%nitrogen%ptr, &
                          lai = lai(1),&
                          lai_cl = lai_cl(1,:), &
                          cumm_lai_cl = cumm_lai_cl(1,:), &
                          leaf_nitrogen_cl = leaf_nitrogen_cl(1,:), &
                          fn_rub_cl = fn_rub_cl(1,:), &
                          fn_et_cl = fn_et_cl(1,:), &
                          fn_pepc_cl = fn_pepc_cl(1,:), &
                          fn_chl_cl = fn_chl_cl(1,:), &
                          fn_oth_cl = fn_oth_cl(1,:), &
                          inquire_n_fractions=.TRUE.)

     !! canopy calculation fluxes
     CALL set_zero_vegetation_fluxes(tile)    !< mo_veg_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(I1,3F16.7)') i, &
           net_assimilation(1), canopy_cond(1) *1000._wp, &
           co2_conc_leaf(1)/co2_mixing_ratio(1) 

  END DO 
  !STOP
  END SUBROUTINE run_test_canopy_model

END MODULE mo_land_model_interface
