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
  USE mo_atmland_memory
  USE mo_atmland_constants
  USE mo_atmland_util
  USE mo_atmland_forcing, ONLY: init_forcing, update_forcing
  USE mo_land_constants
  USE mo_veg_parameters
  USE mo_veg_radiation
  USE mo_veg_photosynthesis

  IMPLICIT NONE

  TYPE(VegetationParameterClass),SAVE :: veg_param     ! data structure for vegetation parameters
  TYPE(TileClass),DIMENSION(ntile), SAVE :: tile       ! data structure for dynamic vegetation and soil variables
  TYPE(t_slm_grid)                  :: grid             ! data structure for coordinates

  INTEGER :: i,lun
  REAL(wp) :: time
  INTEGER :: ntimestep 
  REAL(wp) :: ref_t_air = 25.0_wp          ! temperature in deg C 
  REAL(wp) :: ref_vpd = 1000._wp           ! vpd in Pa
  REAL(wp) :: ref_rh = 0.8_wp              ! relative humidity in fractions 
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
  CALL init_base_constants
  CALL init_atmland_constants
  CALL init_pft_param_table(veg_param)
  CALL init_veg_parameters
  !
  ! 1.1 Initalisation of parameters for generic model forcing 
  !
  time = 0.0                               ! start at mid-night 1.Jan
  dt = 1800.0_wp                           ! time-step in seconds
  ntimestep = NINT(10.0_wp * one_year * one_day / dt)  ! number of time steps
  !
  ! 1.2 Initalisation of met-forcing to standard values 
  ! 
  CALL init_forcing(grid)
  !
  ! 1.3 Initialisation of soil and vegetation (so far relevant) 
  !
  veg_param%tjmax_opt(:) = 20._wp+0.35_wp*(atm2land%t_air - Tzero)
  tile(:)%soil%theta = ref_theta
  tile(1)%plant%lai = ref_lai 
  tile(1)%plant%pool%leaf%nitrogen=ref_nleaf/molar_mass_N*tile(1)%plant%lai
  ! 
  ! 2.0 Time loop 
  ! 
  DO i = 1,ntimestep
     !
     ! 2.1 advance timer
     !
     time = time + dt/one_day
     ! 
     ! 2.2 update time-variant forcing
     !
     CALL update_forcing(time,grid)
     ! 
     ! 3.0 do canopy calculation
     ! 
     CALL update_radiation(tile)

     CALL update_canopy_fluxes(tile,veg_param)
     ! 
     ! 4.0 write output 
     !
     WRITE(12,'(5F16.7)') atm2land%t_air, atm2land%swpar_srf_down * sw2par, & 
           tile(1)%plant%flux%gpp,tile(1)%plant%flux%gc,tile(1)%plant%flux%ci
     !
     ! 5.0 check timer  
     !
     IF(time > one_year) time = time - one_year
  ENDDO

END PROGRAM canopy_model_driver
