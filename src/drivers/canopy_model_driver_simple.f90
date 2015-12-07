!!
!! @brief driver for the stand-alone canopy model (fapar and photosynthesis) 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
PROGRAM canopy_model_driver

  USE mo_kind
  USE mo_base_memory
  USE mo_land_constants
  USE mo_veg_parameters
  USE mo_veg_radiation
  USE mo_veg_photosynthesis

  IMPLICIT NONE

  TYPE(VegetationParameterClass),SAVE :: veg_param     ! data structure for vegetation parameters
  TYPE(TileClass),DIMENSION(ntile), SAVE :: tile       ! data structure for dynamic vegetation and soil variables

  REAL(wp) :: &                                        ! model forcing
      t_air            , &
      q_air            , &
      rain             , &
      snow             , &
      wind_air         , &
      wind_10m         , &
      lw_srf_down      , &
      swvis_srf_down   , &
      swnir_srf_down   , &
      swpar_srf_down   , &
      swpar_fdiffuse   , &
      press_srf        , &
      drag_srf         , &
      t_acoef          , &
      t_bcoef          , &
      q_acoef          , &
      q_bcoef          , &
      pch              , &
      cos_zenith_angle , &
      co2_mixing_ratio
  REAL(wp) :: &                                       ! tentative model output for atmosphere from radiation routine
      alb_vis_dif ,   &
      alb_nir_dif ,   &
      alb_vis_dir ,   &
      alb_nir_dir ,   &  
      ga, &
      vpd
  INTEGER :: i,lun
 
  REAL(wp) :: ref_t_air = 25.0_wp          ! temperature in deg C 
  REAL(wp) :: ref_vpd = 1000._wp           ! vpd in Pa
  REAL(wp) :: ref_swdown = 652.173913_wp   ! PAR in W m-2
  REAL(wp) :: ref_fdiffuse = 0.2_wp        ! diffuse fraction
  REAL(wp) :: ref_cos_angle = 1.0_wp       ! cosine of zenith angle
  REAL(wp) :: ref_co2 = 380.0_wp           ! co2 in ppm
  REAL(wp) :: ref_wind = 1.0_wp            ! wind in m/s
  REAL(wp) :: ref_theta = 100._wp          ! soil water availability in %
  REAL(wp) :: ref_lai = 3._wp              ! LAI
  REAL(wp) :: ref_nleaf = 2.0_wp           ! N leaf per unit leaf area in gN/m2
  !
  ! 1.0 Initialisation of vegetation parameters and constants
  ! 
  CALL init_pft_param_table(veg_param)
  CALL init_veg_parameters
  CALL init_canopy_layer_depth
  !
  ! 1.1 Initalisation of met-forcing which is assumed constant 
  ! 
  t_air = ref_t_air + Tzero
  press_srf = standard_press_srf 
  q_air = calc_spec_humidity_sat(t_air,press_srf) - ref_vpd / press_srf * eps_vpd
  co2_mixing_ratio = ref_co2 
  rain = 0.0_wp
  snow = 0.0_wp
  swvis_srf_down = ref_swdown
  swnir_srf_down = swvis_srf_down * 0.1_wp
  swpar_srf_down = swvis_srf_down 
  swpar_fdiffuse = ref_fdiffuse
  wind_air = ref_wind 
  wind_10m = wind_air
  lw_srf_down = 0.0_wp
  drag_srf = 1.0_wp
  t_acoef = 0.0_wp
  t_bcoef = 0.0_wp
  q_acoef = 0.0_wp
  q_bcoef = 0.0_wp
  pch = 0.0_wp
  cos_zenith_angle = ref_cos_angle 
  ga = .1_wp / ( drag_srf * wind_10m ) 
  !
  ! 1.2 Initialisation of soil and vegetation (so far relevant) 
  !
  veg_param%tjmax_opt(:) = 20._wp+0.35_wp*(t_air - Tzero)
  tile(:)%soil%theta = ref_theta
  tile(1)%plant%lai = ref_lai 
  tile(1)%plant%pool%leaf%nitrogen=ref_nleaf/molar_mass_N14*tile(1)%plant%lai
  ! 
  ! 2.0 Sensitivity plots
  ! 
  ! 2.1 PAR (0-1500 µmol m-2 s-1)
  lun = 10
  WRITE(*,*) 'PAR A GC CI' 
  swvis_srf_down = 0.0_wp
  DO i = 1,15

     swvis_srf_down = swvis_srf_down + 100.0_wp / 2.3_wp
     swnir_srf_down = swvis_srf_down * 0.1_wp
     swpar_srf_down = swvis_srf_down

     CALL update_radiation(tile,swvis_srf_down,swnir_srf_down,swpar_srf_down, &
                              swpar_fdiffuse,cos_zenith_angle,alb_vis_dif,alb_vis_dir, &
                              alb_nir_dif,alb_nir_dir)

     CALL update_canopy_fluxes(tile,t_air,q_air,press_srf,co2_mixing_ratio,ga,veg_param)

     vpd = ( calc_spec_humidity_sat(t_air,press_srf) - q_air ) / eps_vpd * press_srf / 1000.0_wp

     WRITE(*,'(4F16.7)') swpar_srf_down * sw2par, & 
           tile(1)%plant%flux%gpp,tile(1)%plant%flux%gc,tile(1)%plant%flux%ci

  ENDDO
  swvis_srf_down = ref_swdown
  swnir_srf_down = swvis_srf_down * 0.1_wp
  swpar_srf_down = swvis_srf_down 

END PROGRAM canopy_model_driver
