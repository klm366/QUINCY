!>
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
  USE mo_base_constants
!  USE mo_jsb_config_class,      ONLY: t_jsb_config
  ! atmland
  USE mo_atmland_constants
  USE mo_atmland_memory
  USE mo_atmland_util
  USE mo_atmland_forcing
  USE mo_atmland_config_class,  ONLY: atmland_config
  ! land
  USE mo_land_constants         ! still necessary ?
  ! vegetation
  USE mo_veg_parameters
  USE mo_veg_radiation
  USE mo_veg_photosynthesis
  USE mo_veg_config_class,      ONLY: veg_config
  ! shared
  USE mo_namelist,              ONLY: POSITIONED, open_nml, position_nml  ! read model options from file

  IMPLICIT NONE

  TYPE(VegetationParameterClass),SAVE     :: veg_param   ! data structure for vegetation parameters
  TYPE(TileClass),DIMENSION(ntile), SAVE  :: tile        ! data structure for dynamic vegetation and soil variables

  INTEGER :: i
  REAL(wp) :: time
  INTEGER :: ntimestep 
  REAL(wp) :: ref_theta = 100._wp          ! soil water availability in %
  REAL(wp) :: ref_lai = 3._wp              ! LAI
  REAL(wp) :: ref_nleaf = 2.0_wp           ! N leaf per unit leaf area in gN/m2

  ! 
  ! 1.0 Read config paramters from namelist.canopy
  ! 
  veg_config%namelist_filename        = 'namelist.canopy'
  atmland_config%namelist_filename    = 'namelist.canopy'
  CALL veg_config%Init
  CALL atmland_config%Init
  ! 
  ! 1.1 Initialisation of vegetation parameters, constants and forcing
  ! 
  CALL init_base_constants
  CALL init_atmland_constants
  CALL init_veg_parameters(veg_param)
  CALL init_forcing    ! needs atmland_config%Init to run before
  ! 
  ! 1.2 Initialise time step
  ! 
  time = 0.0_wp
  dt = 1800.0_wp     ! time step length in seconds
  ntimestep = NINT(10.0_wp * one_year * one_day / dt)  ! number of time steps
  ! 
  ! 1.3 Initialise tile (Initialisation of soil and vegetation, so far relevant) 
  ! 
  veg_param%tjmax_opt(:) = 20._wp+0.35_wp*(atm2land%t_air - Tzero)  ! conversion from K to °C
  tile(:)%soil%theta = ref_theta
  tile(1)%plant%lai = ref_lai 
  IF(veg_config%pft_id==0) THEN
    tile(1)%plant%pool%leaf%nitrogen=ref_nleaf/molar_mass_N*tile(1)%plant%lai
  ELSE
    tile(1)%plant%pool%leaf%nitrogen=tile(1)%plant%lai & 
             / veg_param%sla(veg_config%pft_id)/veg_param%cn_leaf(veg_config%pft_id)
    tile(1)%plant%pool%leaf%phosphorus=tile(1)%plant%pool%leaf%nitrogen & 
             / veg_param%np_leaf(veg_config%pft_id)
    tile(1)%plant%pftid=veg_config%pft_id
  ENDIF
  ! 
  ! 1.4 open input files (and read the header into dummy variables); see mo_atmland_forcing
  ! 
  CALL open_forcing_input_file  ! only if: options%forcing%option_forcing == 'readfile'
  ! 
  ! 2.0 Time loop 
  ! 
  DO i = 1,ntimestep
     ! 
     ! 2.1 advance timer
     ! 
     time = time + dt/one_day 
     ! 
     ! 2.2 set orbital parameters and update time-variant forcing; 
     ! 
     CALL update_forcing(time)  ! see mo_atmland_forcing
     ! 
     ! 3.0 canopy calculation fluxes
     ! 
     CALL update_radiation(tile)                 ! atm2land land2atm are global TYPES, see mo_atmland_forcing
     CALL update_canopy_fluxes(tile(1), veg_param)  ! atm2land is a global TYPE, see mo_atmland_forcing
     ! 
     ! 4.0 write output  -> filename: fort.20
     !
     WRITE(20,'(5F16.7)') atm2land%t_air, atm2land%swpar_srf_down * sw2par, & 
           tile(1)%plant%flux%gpp,tile(1)%plant%flux%gc,tile(1)%plant%flux%ci
     !
     ! 5.0 check timer  
     !
     IF(time > one_year) time = time - one_year
  ENDDO

END PROGRAM canopy_model_driver
