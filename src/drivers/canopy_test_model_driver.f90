!>
!! @brief driver for the stand-alone canopy model (fapar and photosynthesis) 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
PROGRAM canopy_test_model_driver

  USE mo_kind
  USE mo_base_memory
  USE mo_atmland_util
  USE mo_atmland_memory
  USE mo_atmland_forcing
  USE mo_atmland_constants
  USE mo_atmland_config_class,  ONLY: atmland_config
  USE mo_land_constants   ! still necessary ?
  USE mo_veg_parameters
  USE mo_veg_radiation
  USE mo_veg_photosynthesis
  USE mo_veg_config_class,      ONLY: veg_config

  IMPLICIT NONE

  TYPE(VegetationParameterClass),SAVE    :: veg_param     ! data structure for vegetation parameters
  TYPE(TileClass),DIMENSION(ntile), SAVE :: tile       ! data structure for dynamic vegetation and soil variables

  REAL(wp) :: vpd
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
  INTEGER  :: ref_pft = 1                  ! PFT identifier

  ! 
  ! 1.0 Initialisation of vegetation parameters and constants
  ! 
  CALL init_base_constants
  CALL init_atmland_constants
  CALL init_veg_parameters(veg_param)
  !
  ! 1.1 Initalisation of met-forcing which is assumed constant 
  ! 
  atm2land%t_air = def_t_air
  atm2land%press_srf = standard_press_srf
  atm2land%q_air = calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf) - def_vpd / atm2land%press_srf * eps_vpd
  atm2land%co2_mixing_ratio = def_co2
  atm2land%rain = 0.0_wp
  atm2land%snow = 0.0_wp
  atm2land%swvis_srf_down = def_swdown
  atm2land%swnir_srf_down = atm2land%swvis_srf_down * 0.1_wp
  atm2land%swpar_srf_down = atm2land%swvis_srf_down
  atm2land%swpar_fdiffuse = def_fdiffuse
  atm2land%wind_air = def_wind
  atm2land%wind_10m = atm2land%wind_air
  atm2land%lw_srf_down = 0.0_wp
  atm2land%drag_srf = 1.0_wp
  atm2land%t_acoef = 0.0_wp
  atm2land%t_bcoef = 0.0_wp
  atm2land%q_acoef = 0.0_wp
  atm2land%q_bcoef = 0.0_wp
  atm2land%pch = 0.0_wp
  atm2land%cos_zenith_angle = def_cos_angle
  atm2land%ga = .1_wp / ( atm2land%drag_srf * atm2land%wind_10m )
  grid%lat = 50.0
  grid%lon = 0.0
  ! 
  ! 1.2 Initialisation of soil and vegetation (so far relevant) 
  ! 
  tile(:)%soil%theta = ref_theta
  tile(1)%plant%lai = ref_lai 
  IF(ref_pft==0) THEN
    tile(1)%plant%pool%leaf%nitrogen=ref_nleaf/molar_mass_N*tile(1)%plant%lai
  ELSE
    tile(1)%plant%pool%leaf%nitrogen=tile(1)%plant%lai & 
            / veg_param%sla(ref_pft)/veg_param%cn_leaf(ref_pft)
    tile(1)%plant%pool%leaf%phosphorus=tile(1)%plant%pool%leaf%nitrogen & 
             / veg_param%np_leaf(ref_pft)
    tile(1)%plant%pftid=ref_pft
  ENDIF
  ! 
  ! 2.0 Sensitivity plots
  ! 
  ! 2.1 PAR (0-1500 µmol m-2 s-1)
  lun = 10
  WRITE(lun,*) 'PAR A GC CI' 
  atm2land%swvis_srf_down = 0.0_wp
  DO i = 1,15

     atm2land%swvis_srf_down = atm2land%swvis_srf_down + 100.0_wp / 2.3_wp
     atm2land%swnir_srf_down = atm2land%swvis_srf_down * 0.1_wp
     atm2land%swpar_srf_down = atm2land%swvis_srf_down

     CALL update_radiation(tile)
     CALL update_canopy_fluxes(tile(1),veg_param)

     WRITE(lun,'(4F16.7)') atm2land%swpar_srf_down * sw2par, & 
           tile(1)%plant%flux%gpp,tile(1)%plant%gc*1000._wp,tile(1)%plant%flux%ci !,&
!           tile(1)%plant%flux%an,
!           tile(1)%plant%flux%ci/atm2land%co2_mixing_ratio,tile(1)%plant%pool%leaf%nitrogen*molar_mass_N,&
!           tile(1)%plant%flux%an/(tile(1)%plant%pool%leaf%nitrogen*molar_mass_N)

  ENDDO
  atm2land%swvis_srf_down = ref_swdown
  atm2land%swnir_srf_down = atm2land%swvis_srf_down * 0.1_wp
  atm2land%swpar_srf_down = atm2land%swvis_srf_down 

  ! 2.2 Temperature (-5,60 °C)
  lun = lun + 1
  WRITE(lun,*) 'T A GC CI' 
  atm2land%t_air = - 5.0_wp + Tzero
  DO i = 1,70

     atm2land%t_air = atm2land%t_air + 1.0_wp 
     atm2land%q_air = 0.8_wp * calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf)
     veg_param%tjmax_opt(:) = 20._wp+0.35_wp*(atm2land%t_air - Tzero)

     CALL update_radiation(tile)
     CALL update_canopy_fluxes(tile(1),veg_param)

     WRITE(lun,'(4F16.7)') atm2land%t_air, &  
           tile(1)%plant%flux%gpp,tile(1)%plant%gc*1000._wp,tile(1)%plant%flux%ci

  ENDDO
  atm2land%t_air = ref_t_air + Tzero
  atm2land%q_air = calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf) - ref_vpd / atm2land%press_srf * eps_vpd

  ! 2.3 atmospheric CO2 (175 - 1000 ppm) 
  lun = lun + 1 
  WRITE(lun,*) 'CO2 A GC CI' 
  atm2land%co2_mixing_ratio = 150._wp
  DO i = 1,34

     atm2land%co2_mixing_ratio = atm2land%co2_mixing_ratio + 25._wp

     CALL update_radiation(tile)
     CALL update_canopy_fluxes(tile(1),veg_param)

     WRITE(lun,'(4F16.7)') atm2land%co2_mixing_ratio, &
           tile(1)%plant%flux%gpp,tile(1)%plant%gc*1000._wp,tile(1)%plant%flux%ci

  ENDDO
  atm2land%co2_mixing_ratio = ref_co2 

  ! 2.4 VPD (0-3.3 kPa) 
  lun = lun + 1
  WRITE(lun,*) 'VPD A GC CI' 
  DO i = 1,30
     atm2land%q_air = calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf) - i * 100.0_wp / atm2land%press_srf * eps_vpd

     CALL update_radiation(tile)
     CALL update_canopy_fluxes(tile(1),veg_param)

     vpd = ( calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf) - atm2land%q_air ) / eps_vpd * atm2land%press_srf / 1000.0_wp

     WRITE(lun,'(4F16.7)') vpd, &
           tile(1)%plant%flux%gpp,tile(1)%plant%gc*1000._wp,tile(1)%plant%flux%ci

  ENDDO
  atm2land%q_air = calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf) - ref_vpd / atm2land%press_srf * eps_vpd
  vpd = ( calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf) - atm2land%q_air ) / eps_vpd * atm2land%press_srf / 1000.0_wp

  ! 2.5 Theta (0-100 %)
  lun = lun + 1
  WRITE(lun,*) 'THETA A GC CI' 
  tile(1)%soil%theta = 0.0_wp
  DO i = 1,100

     tile(1)%soil%theta = tile(1)%soil%theta + 1.0_wp 

     CALL update_radiation(tile)
     CALL update_canopy_fluxes(tile(1),veg_param)

     WRITE(lun,'(4F16.7)') tile(1)%soil%theta, &
           tile(1)%plant%flux%gpp,tile(1)%plant%gc*1000._wp,tile(1)%plant%flux%ci

  ENDDO
  tile(1)%soil%theta = ref_theta

  ! 2.6 LAI (0-11)
  lun = lun + 1
  WRITE(lun,*) 'LAI A GC CI'
  tile(1)%plant%lai = 0.0_wp
  DO i = 1,22

     tile(1)%plant%lai = tile(1)%plant%lai + 0.5_wp 
     tile(1)%plant%pool%leaf%nitrogen=ref_nleaf/molar_mass_N*tile(1)%plant%lai

     CALL update_radiation(tile)
     CALL update_canopy_fluxes(tile(1),veg_param)

     WRITE(lun,'(4F16.7)') tile(1)%plant%lai, &
           tile(1)%plant%flux%gpp,tile(1)%plant%gc*1000._wp,tile(1)%plant%flux%ci

  ENDDO
  tile(1)%plant%lai = ref_lai 

  ! 2.7 leaf N (0.5-3 g N / m2(LAI))
  lun = lun + 1
  WRITE(lun,*) 'Nleaf A GC CI'
  DO i = 1,24

     tile(1)%plant%pool%leaf%nitrogen=0.125_wp*(i)/molar_mass_N*tile(1)%plant%lai

     CALL update_radiation(tile)
     CALL update_canopy_fluxes(tile(1),veg_param)

     WRITE(lun,'(4F16.7)') tile(1)%plant%pool%leaf%nitrogen/tile(1)%plant%lai*molar_mass_N, & 
           tile(1)%plant%flux%gpp,tile(1)%plant%gc*1000._wp,tile(1)%plant%flux%ci

  ENDDO
  tile(1)%plant%pool%leaf%nitrogen=ref_nleaf/molar_mass_N*tile(1)%plant%lai

  ! 2.8 PFT
  lun = lun + 1
  WRITE(lun,*) 'Nleaf A GC CI'
  DO i = 1,npft
     tile(1)%plant%pool%leaf%nitrogen=tile(1)%plant%lai/veg_param%sla(i)/veg_param%cn_leaf(i)
     tile(1)%plant%pool%leaf%phosphorus=tile(1)%plant%pool%leaf%nitrogen/veg_param%np_leaf(i)
     tile(1)%plant%pftid=i

     CALL update_radiation(tile)
     CALL update_canopy_fluxes(tile(1),veg_param)

     WRITE(lun,'(i,3F16.7)') tile(1)%plant%pftid, &
           tile(1)%plant%flux%gpp,tile(1)%plant%gc*1000._wp,tile(1)%plant%flux%ci

  ENDDO 
STOP
  

END PROGRAM canopy_test_model_driver
