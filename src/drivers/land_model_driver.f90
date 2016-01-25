!!
!! @brief driver for the stand-alone land model
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
PROGRAM land_model_driver
  USE mo_kind
  USE mo_base_constants
  USE mo_atmland_constants
  USE mo_land_constants
  USE mo_land_model_interface, ONLY: land_model_interface 
  USE mo_veg_parameters

  IMPLICIT NONE

   TYPE(VegetationParameterClass),SAVE :: veg_param    ! data structure for vegetation parameter

   REAL(wp) :: &                                       ! model forcing
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

  REAL(wp) :: ref_t_air = 25.0_wp          ! temperature in deg C
  REAL(wp) :: ref_vpd = 1000._wp           ! vpd in Pa
  REAL(wp) :: ref_swdown = 652.173913_wp   ! PAR in W m-2
  REAL(wp) :: ref_fdiffuse = 0.2_wp        ! diffuse fraction
  REAL(wp) :: ref_cos_angle = 1.0_wp       ! cosine of zenith angle
  REAL(wp) :: ref_co2 = 380.0_wp           ! co2 in ppm
  REAL(wp) :: ref_wind = 1.0_wp            ! wind in m/s
  REAL(wp) :: ref_theta = 100._wp          ! soil water availability in %
  REAL(wp) :: ref_lai = 3._wp              ! LAI
  REAL(wp) :: ref_nleaf = 1.2_wp           ! N leaf per unit leaf area in gN/m2

  INTEGER :: i

  dt = 1                   ! model time step (days)
  one_year = 365.25_wp     ! length of one year (days)

  !
  ! 1.0 Initialisation of vegetation parameters and constants
  !
  CALL init_base_constants
  CALL init_atmland_constants
  CALL init_pft_param_table(veg_param)
  CALL init_veg_parameters
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

  DO i = 1,100

    CALL land_model_interface(t_air,q_air,rain,snow,wind_air,wind_10m,lw_srf_down,swvis_srf_down, &
                                  swnir_srf_down,swpar_srf_down,swpar_fdiffuse,press_srf,drag_srf,t_acoef,t_bcoef, &
                                  q_acoef,q_bcoef,pch,cos_zenith_angle,co2_mixing_ratio,veg_param)

  ENDDO

END PROGRAM land_model_driver
