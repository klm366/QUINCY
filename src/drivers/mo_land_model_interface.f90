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
                                  time_day, time_year, timestep_prop_day &
                                 )

    USE mo_kind
    USE mo_base_memory,             ONLY: t_slm_tile, t_model
    USE mo_base_config_class,       ONLY: base_config
    USE mo_atmland_memory
    USE mo_veg_output,              ONLY: write_output_pool, write_output_flux
    USE mo_veg_config_class,        ONLY: veg_config  ! used to switch on/off the update_soil_pool via flag_base_soil_model

    ! used only for the canopy model:
    USE mo_atmland_constants,       ONLY: rad2par
    USE mo_base_constants,          ONLY: one_year

    ! Process interfaces:
    USE mo_rad_interface
    USE mo_veg_interface
    USE mo_assimi_interface
    !USE mo_soil_interface


    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)      :: tile              !< data structure for dynamic vegetation and soil variables
    TYPE(t_model),INTENT(inout)         :: model             !< model structure
    REAL(wp),INTENT(in)                 :: time_day          !< displays the 'real time' in days; time_day == 1.0 after one day has passed
    INTEGER,INTENT(in)                  :: time_year         !< displays the 'real time' in years; the 1st year of a simulation == 1 (not zero)
    REAL(wp),INTENT(in)                 :: timestep_prop_day !< proportion of one day that is simulated with one timestep (dt/one_day)
    ! 0.2 Local
    ! ...


    !> @todo mo_land_model_interface:land_model_interface put all Tasks into a Queue accordingly to the model selection
    !! (land plant soil canopy test_canopy)
    SELECT CASE(TRIM(base_config%slm_model_name))
    !! land ---------------------------------------------------------------------------------------
    CASE('land')
    !! set vegetation fluxes to zero
    !! 
    !! ==> integrate into vegetation process
      CALL tile%plant%flux%set_zero              ! mo_veg_memory

    !! run Tasks
    !!
      CALL update_radiation       (tile)
      CALL update_vegsoil_bgc     (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_phenology       (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_canopy_fluxes   (tile, model%lctlib(tile%lcts%lib_id))
      CALL calc_growth            (tile, model%lctlib(tile%lcts%lib_id))
      CALL calc_turnover          (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_veg_dynamics    (tile, model%lctlib(tile%lcts%lib_id))
      IF (veg_config%flag_base_soil_model) THEN   ! switch: soil model 
        CALL update_soil_pool       (tile, model%lctlib(tile%lcts%lib_id))
      END IF
      CALL update_vegsoil_physics (tile)
      CALL update_veg_pools       (tile, model%lctlib(tile%lcts%lib_id))

    !! accumulation - daily & moving average
    !! 
      CALL tile%plant%daily_average%calc_avrg(tile%plant, time_day, timestep_prop_day) !< daily avrg of the prev day
      CALL tile%plant%moving_average%update_veg(tile%plant)               !< moving averages - vegetation
      !CALL tile%plant%moving_average%update_soil(tile%plant)             !< moving averages - soil -> plant N & P uptake
      CALL tile%plant%moving_average%update_tile(atm2land%t_air, atm2land%press_srf, atm2land%q_air) !< moving averages - atmland

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
      CALL tile%plant%flux%set_zero              ! mo_veg_memory

    !! run Tasks
    !!
      CALL update_radiation       (tile)
      CALL update_vegsoil_bgc     (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_phenology       (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_canopy_fluxes   (tile, model%lctlib(tile%lcts%lib_id))
      CALL calc_growth            (tile, model%lctlib(tile%lcts%lib_id))
      CALL calc_turnover          (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_veg_dynamics    (tile, model%lctlib(tile%lcts%lib_id))
      IF (veg_config%flag_base_soil_model) THEN   ! switch: soil model 
        CALL update_soil_pool       (tile, model%lctlib(tile%lcts%lib_id))
      END IF
      CALL update_vegsoil_physics (tile)
      CALL update_veg_pools       (tile, model%lctlib(tile%lcts%lib_id))

    !! accumulation - daily & moving average
    !! 
      CALL tile%plant%daily_average%calc_avrg(tile%plant, time_day, timestep_prop_day) !< daily avrg of the prev day
      CALL tile%plant%moving_average%update_veg(tile%plant)               !< moving averages - vegetation
      !CALL tile%plant%moving_average%update_soil(tile%plant)             !< moving averages - soil -> plant N & P uptake
      CALL tile%plant%moving_average%update_tile(atm2land%t_air, atm2land%press_srf, atm2land%q_air) !< moving averages - atmland

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
      CALL tile%plant%flux%set_zero              ! mo_veg_memory

    !! run Tasks
    !!
      CALL update_radiation       (tile)
      CALL update_canopy_fluxes   (tile, model%lctlib(tile%lcts%lib_id))
      CALL update_vegsoil_physics (tile)

    !! accumulation - daily & moving average
    !! 
      CALL tile%plant%daily_average%calc_avrg(tile%plant, time_day, timestep_prop_day) !< daily avrg of the prev day
      CALL tile%plant%moving_average%update_veg(tile%plant)               !< moving averages - vegetation
      !CALL tile%plant%moving_average%update_soil(tile%plant)             !< moving averages - soil -> plant N & P uptake
      CALL tile%plant%moving_average%update_tile(atm2land%t_air, atm2land%press_srf, atm2land%q_air) !< moving averages - atmland

    !! output fluxes -> filename: fort.20
    !! output also some atm2land variables
    !! 
      WRITE(20,'(6F16.7)') time_day+time_year*one_year, atm2land%t_air, atm2land%swpar_srf_down * rad2par, & 
                           tile%plant%flux%gpp, tile%plant%gc, tile%plant%flux%ci

    CASE('test_canopy')
    !! test_canopy --------------------------------------------------------------------------------
      CALL run_test_canopy_model  (tile, model)

    CASE DEFAULT
        CALL finish('land_model_interface','Invalid SLM model name. Try either land, plant, soil, canopy or test_canopy.')
   END SELECT
  END SUBROUTINE land_model_interface


  !-----------------------------------------------------------------------------------------------------
  !> Run the test_canopy model from the land model interface
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE run_test_canopy_model(tile, model)

    USE mo_base_memory,                     ONLY: t_slm_tile, t_model
    USE mo_base_constants,                  ONLY: npft
    USE mo_jsb_lctlib_class,                ONLY: t_lctlib_element
    USE mo_atmland_memory,                  ONLY: atm2land
    USE mo_atmland_constants,               ONLY: rad2par, frac_vis_swtot_srf, frac_par_swvis_srf,eps_vpd
    USE mo_atmland_util,                    ONLY: calc_spec_humidity_sat
    USE mo_shared_constants,                ONLY: molar_mass_N, Tzero 
    ! Process interfaces:
    USE mo_rad_interface                    ! update_radiation
    USE mo_assimi_interface                 ! update_canopy_fluxes
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

    ! parameters (specific to the test_canopy model)
    ! these parameters were also defined in mo_jsb_model_init:init_slm_test_canopy !!
    REAL(wp),SAVE :: ref_t_air = 25.0_wp          ! temperature in deg C 
    REAL(wp),SAVE :: ref_vpd = 1000._wp           ! vpd in Pa
    REAL(wp),SAVE :: ref_swdown = 652.173913_wp   ! PAR in W m-2
    REAL(wp),SAVE :: ref_fdiffuse = 0.2_wp        ! diffuse fraction
    REAL(wp),SAVE :: ref_cos_angle = 1.0_wp       ! cosine of zenith angle
    REAL(wp),SAVE :: ref_co2 = 380.0_wp           ! co2 in ppm
    REAL(wp),SAVE :: ref_wind = 1.0_wp            ! wind in m/s
    REAL(wp),SAVE :: ref_theta = 100._wp          ! soil water availability in %
    REAL(wp),SAVE :: ref_lai = 3._wp              ! LAI
    REAL(wp),SAVE :: ref_nleaf = 2.0_wp           ! N leaf per unit leaf area in gN/m2
    INTEGER ,SAVE :: ref_pft = 1                  ! PFT identifier

    !  bi-directional pointer
    lctlib => model%lctlib(ref_pft)
    tile%plant%moving_average%t_jmax_opt = 20._wp+0.35_wp*(atm2land%t_air - Tzero)

  !> 2.0 Sensitivity plots
  !!
  !!  

  !> 2.1 PAR (0-1500 µmol m-2 s-1)
  !! 
  !! 
  lun = 10
  WRITE(lun,*) 'PAR A GC CI' 
  atm2land%swvis_srf_down = 0.0_wp 
  DO i = 1,16

     atm2land%swpar_srf_down = atm2land%swpar_srf_down + 100.0_wp / rad2par 
     atm2land%swvis_srf_down = atm2land%swpar_srf_down / frac_par_swvis_srf 
     atm2land%swnir_srf_down = atm2land%swvis_srf_down * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf

     !! canopy calculation fluxes
     CALL tile%plant%flux%set_zero
     CALL update_radiation(tile)                       !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, lctlib)        !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') atm2land%swpar_srf_down * rad2par, & 
           tile%plant%flux%an,tile%plant%gc*1000._wp,tile%plant%flux%ci !,&
!           tile%plant%flux%an,
!           tile%plant%flux%ci/atm2land%co2_mixing_ratio,tile%plant%pool%leaf%nitrogen*molar_mass_N,&
!           tile%plant%flux%an/(tile%plant%pool%leaf%nitrogen*molar_mass_N)

  END DO
  atm2land%swpar_srf_down = ref_swdown 
  atm2land%swvis_srf_down = atm2land%swpar_srf_down / frac_par_swvis_srf 
  atm2land%swnir_srf_down = atm2land%swvis_srf_down * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf

  !> 2.2 Temperature (-5,60 °C)
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'T A GC CI' 
  atm2land%t_air = - 5.0_wp + Tzero
  hlp1 = tile%plant%moving_average%t_jmax_opt
  DO i = 1,70

     atm2land%t_air = atm2land%t_air + 1.0_wp 
     atm2land%q_air = 0.8_wp * calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf)
     tile%plant%moving_average%t_jmax_opt = 20._wp+0.35_wp*(atm2land%t_air - Tzero)

     !! canopy calculation fluxes
     CALL tile%plant%flux%set_zero
     CALL update_radiation(tile)                       !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, lctlib)        !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') atm2land%t_air, &  
           tile%plant%flux%an,tile%plant%gc*1000._wp,tile%plant%flux%ci

  END DO
  atm2land%t_air    = ref_t_air + Tzero
  atm2land%q_air    = calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf) - ref_vpd / atm2land%press_srf * eps_vpd
  tile%plant%moving_average%t_jmax_opt = hlp1 

  !> 2.3 atmospheric CO2 (175 - 1000 ppm) 
  !! 
  !! 
  lun = lun + 1 
  WRITE(lun,*) 'CO2 A GC CI' 
  atm2land%co2_mixing_ratio = 150._wp
  DO i = 1,34

     atm2land%co2_mixing_ratio = atm2land%co2_mixing_ratio + 25._wp

     !! canopy calculation fluxes
     CALL tile%plant%flux%set_zero
     CALL update_radiation(tile)                       !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, lctlib)        !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') atm2land%co2_mixing_ratio, &
           tile%plant%flux%an,tile%plant%gc*1000._wp,tile%plant%flux%ci

  END DO
  atm2land%co2_mixing_ratio = ref_co2 

  !> 2.4 VPD (0-3.3 kPa) 
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'VPD A GC CI' 
  DO i = 1,30
     atm2land%q_air = calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf) - i * 100.0_wp / atm2land%press_srf * eps_vpd

     !! canopy calculation fluxes
     CALL tile%plant%flux%set_zero
     CALL update_radiation(tile)                       !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, lctlib)        !< mo_assimi_update_canopy_fluxes

     vpd = ( calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf) - atm2land%q_air ) / eps_vpd * atm2land%press_srf / 1000.0_wp

     WRITE(lun,'(4F16.7)') vpd, &
           tile%plant%flux%an,tile%plant%gc*1000._wp,tile%plant%flux%ci

  END DO
  atm2land%q_air = calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf) - ref_vpd / atm2land%press_srf * eps_vpd
  vpd = ( calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf) - atm2land%q_air ) / eps_vpd * atm2land%press_srf / 1000.0_wp

  !> 2.5 Theta (0-100 %)
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'THETA A GC CI' 
  tile%soil%theta = 0.0_wp
  DO i = 1,100

     tile%soil%theta = tile%soil%theta + 1.0_wp 

     !! canopy calculation fluxes
     CALL tile%plant%flux%set_zero
     CALL update_radiation(tile)                       !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, lctlib)        !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') tile%soil%theta, &
           tile%plant%flux%an,tile%plant%gc*1000._wp,tile%plant%flux%ci

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

     !! canopy calculation fluxes
     CALL tile%plant%flux%set_zero
     CALL update_radiation(tile)                       !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, lctlib)        !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') tile%plant%lai, &
           tile%plant%flux%an,tile%plant%gc*1000._wp,tile%plant%flux%ci

  END DO
  tile%plant%lai = ref_lai 

  !> 2.7 leaf N (0.5-3 g N / m2(LAI))
  !! 
  !! 
  lun = lun + 1
  WRITE(lun,*) 'Nleaf A GC CI'
  DO i = 1,24

     tile%plant%pool%leaf%nitrogen=0.125_wp*(i)/molar_mass_N*tile%plant%lai

     !! canopy calculation fluxes
     CALL tile%plant%flux%set_zero
     CALL update_radiation(tile)                       !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, lctlib)        !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') tile%plant%pool%leaf%nitrogen/tile%plant%lai*molar_mass_N, & 
           tile%plant%flux%an,tile%plant%gc*1000._wp,tile%plant%flux%ci

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

     !! canopy calculation fluxes
     CALL tile%plant%flux%set_zero
     CALL update_radiation(tile)                             !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile, model%lctlib(i))        !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(I1,3F16.7)') tile%plant%pftid, &
           tile%plant%flux%an,tile%plant%gc*1000._wp,tile%plant%flux%ci

  END DO 
  !STOP
  END SUBROUTINE run_test_canopy_model

END MODULE mo_land_model_interface
