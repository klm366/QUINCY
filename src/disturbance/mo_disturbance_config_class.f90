!>
!! @brief Contains structures and methods for disturbance config
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (06.04.2017)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_disturb_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_disturb_config

  TYPE, EXTENDS(t_jsb_config) :: t_disturb_config
    INTEGER           :: fire_algorithm      ! Select fire algorithm: 0 : none, 1 : jsbach (default),
                                             ! 2 : Arora and Boer (2005) (not yet), 3 : Thonicke et al.(2001)
    INTEGER           :: windbreak_algorithm ! Select windbreak algorithm: 0 : none, 1 : jsbach (default)
    LOGICAL           :: lburn_pasture
    REAL(wp)          :: wnd_threshold, &    ! factor by which the previous day maximum wind speed
                                             ! must be larger than the climatological daily maximum
                                             ! wind speed to allow any wind damage
                         wnd_damage_scale    ! scaling factor for wind damage

   REAL(wp)           :: fire_litter_threshold,  & ! minimal amount of litter [mol(C)/m^2(grid box)] for fire
                         fire_rel_hum_threshold, & ! maximal relative humidity for fire
                         fire_minimum_woody,     & ! minimal fraction of act_fpc of woody PFT to be burned each year
                         fire_minimum_grass,     & ! minimal fraction of act_fpc of grass PFT to be burned each year
                         fire_tau_woody,         & ! return period of fire for woody PFT [year] at 0% relative humidity
                         fire_tau_grass            ! return period of fire for grass PFT [year] at 0% relative humidity

    INTEGER           :: d_climate_input,        & ! <Rackhun Son>
                         d_climate_hidden,       & ! <Rackhun Son>
                         d_vegetation_input,     & ! <Rackhun Son>
                         d_vegetation_hidden,    & ! <Rackhun Son>
                         d_lnt_input,            & ! <Rackhun Son>
                         d_lnt_hidden,           & ! <Rackhun Son>
                         d_output,               & ! <Rackhun Son>
                         d_emb_input,            & ! <Rackhun Son>
                         d_emb_hidden,           & ! <Rackhun Son>
                         d_emb_output              ! <Rackhun Son>

    REAL(wp), ALLOCATABLE :: climate_i2i_weight(:,:), climate_i2i_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_h2i_weight(:,:), climate_h2i_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_i2f_weight(:,:), climate_i2f_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_h2f_weight(:,:), climate_h2f_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_i2g_weight(:,:), climate_i2g_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_h2g_weight(:,:), climate_h2g_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_i2o_weight(:,:), climate_i2o_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_h2o_weight(:,:), climate_h2o_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_out_weight(:,:), climate_out_bias(:) ! dl fire model parameters <Rackhun Son>

    REAL(wp), ALLOCATABLE :: vegetation_i2i_weight(:,:), vegetation_i2i_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_h2i_weight(:,:), vegetation_h2i_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_i2f_weight(:,:), vegetation_i2f_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_h2f_weight(:,:), vegetation_h2f_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_i2g_weight(:,:), vegetation_i2g_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_h2g_weight(:,:), vegetation_h2g_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_i2o_weight(:,:), vegetation_i2o_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_h2o_weight(:,:), vegetation_h2o_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_out_weight(:,:), vegetation_out_bias(:) ! dl fire model parameters <Rackhun Son>

    REAL(wp), ALLOCATABLE :: lightning_topography_layer1_weight(:,:), lightning_topography_layer1_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: lightning_topography_layer2_weight(:,:), lightning_topography_layer2_bias(:) ! dl fire model parameters <Rackhun Son>

    REAL(wp), ALLOCATABLE :: climate_emb_layer1_weight(:,:), climate_emb_layer1_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_emb_layer2_weight(:,:), climate_emb_layer2_bias(:) ! dl fire model parameters <Rackhun Son>
    
    REAL(wp), ALLOCATABLE :: vegetation_emb_layer1_weight(:,:), vegetation_emb_layer1_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_emb_layer2_weight(:,:), vegetation_emb_layer2_bias(:) ! dl fire model parameters <Rackhun Son>

    REAL(wp), ALLOCATABLE :: climate_min_max(:,:)                ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_min_max(:,:)             ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: embedding_min_max(:,:)              ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: lightning_topography_min_max(:,:)   ! dl fire model parameters <Rackhun Son>
    
    !REAL(wp), ALLOCATABLE :: t_air_minmax(:), t_air_ano_minmax(:), t_air_clim_minmax(:) ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: q_rel_minmax(:), q_rel_ano_minmax(:), q_rel_clim_minmax(:) ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: q_minmax(:),     q_ano_minmax(:),     q_clim_minmax(:)     ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: swl1_minmax(:),  swl1_ano_minmax(:),  swl1_clim_minmax(:)  ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: swl2_minmax(:),  swl2_ano_minmax(:),  swl2_clim_minmax(:)  ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: swl3_minmax(:),  swl3_ano_minmax(:),  swl3_clim_minmax(:)  ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: swl4_minmax(:),  swl4_ano_minmax(:),  swl4_clim_minmax(:)  ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: lai_minmax(:),   lai_ano_minmax(:),   lai_clim_minmax(:)   ! dl fire model parameters <Rackhun Son>

    !REAL(wp), ALLOCATABLE :: wind_10m_minmax(:)                       ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: rain_minmax(:)                           ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: population_density_minmax(:)             ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: lightning_climatology_minmax(:)          ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: elevation_minmax(:)                      ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: agb_minmax(:)                            ! dl fire model parameters <Rackhun Son>
    
    !REAL(wp), ALLOCATABLE :: pft1_minmax(:)                           ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft2_minmax(:)                           ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft3_minmax(:)                           ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft4_minmax(:)                           ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft5_minmax(:)                           ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft6_minmax(:)                           ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft7_minmax(:)                           ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft8_minmax(:)                           ! dl fire model parameters <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft9_minmax(:)                           ! dl fire model parameters <Rackhun Son>

  CHARACTER(len=filename_max) :: bc_lightning_filename         ! lightning frequency (daily) for fire <Rackhun Son>
  CHARACTER(len=filename_max) :: bc_topography_filename        ! topographical complexity for fire    <Rackhun Son>
  !CHARACTER(len=filename_max) :: bc_population_filename        ! population density (daily) for fire  <Rackhun Son>
  CHARACTER(len=filename_max) :: dl_fire_model_params_filename ! dl fire model parameters <Rackhun Son>
  CHARACTER(len=filename_max) :: dl_local_clim_t_air_filename   ! OBS local mean for dl fire model <Rackhun Son>
  CHARACTER(len=filename_max) :: dl_local_clim_rH_cf_filename   ! OBS local mean for dl fire model <Rackhun Son>
  CHARACTER(len=filename_max) :: dl_local_clim_q_filename       ! OBS local mean for dl fire model <Rackhun Son>
  CHARACTER(len=filename_max) :: dl_local_clim_ws10_filename    ! OBS local mean for dl fire model <Rackhun Son>
  CHARACTER(len=filename_max) :: dl_local_clim_swvl1_filename   ! OBS local mean for dl fire model <Rackhun Son>
  CHARACTER(len=filename_max) :: dl_local_clim_swvl2_filename   ! OBS local mean for dl fire model <Rackhun Son>
  CHARACTER(len=filename_max) :: dl_local_clim_swvl3_filename   ! OBS local mean for dl fire model <Rackhun Son>
  CHARACTER(len=filename_max) :: dl_local_clim_swvl4_filename   ! OBS local mean for dl fire model <Rackhun Son>
  CHARACTER(len=filename_max) :: dl_local_clim_lai_filename     ! OBS local mean for dl fire model <Rackhun Son>


   CONTAINS
     PROCEDURE :: Init => Init_disturb_config
  END type t_disturb_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize disturb process
  !!
  !! @param[inout]     config     Configuration type of process (t_disturb_config)
  !!
  SUBROUTINE Init_disturb_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid              ! <Rackhun Son>
    USE mo_jsb_grid,           ONLY: Register_vgrid                      ! <Rackhun Son>
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC                       ! <Rackhun Son>
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input ! <Rackhun Son>

    CLASS(t_disturb_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    CHARACTER(len=filename_max) :: bc_lightning_filename         ! lightning frequency (daily) for fire  <Rackhun Son>
    CHARACTER(len=filename_max) :: bc_topography_filename        ! topographical complexity for fire      <Rackhun Son>
    !CHARACTER(len=filename_max) :: bc_population_filename        ! population density (daily) for fire  <Rackhun Son>
    CHARACTER(len=filename_max) :: dl_fire_model_params_filename ! dl fire model parameters               <Rackhun Son>

    CHARACTER(len=filename_max) :: dl_local_clim_t_air_filename   ! OBS local mean for dl fire model      <Rackhun Son>
    CHARACTER(len=filename_max) :: dl_local_clim_rH_cf_filename   ! OBS local mean for dl fire model      <Rackhun Son>
    CHARACTER(len=filename_max) :: dl_local_clim_q_filename       ! OBS local mean for dl fire model      <Rackhun Son>
    CHARACTER(len=filename_max) :: dl_local_clim_ws10_filename    ! OBS local mean for dl fire model      <Rackhun Son>
    CHARACTER(len=filename_max) :: dl_local_clim_swvl1_filename   ! OBS local mean for dl fire model      <Rackhun Son>
    CHARACTER(len=filename_max) :: dl_local_clim_swvl2_filename   ! OBS local mean for dl fire model      <Rackhun Son>
    CHARACTER(len=filename_max) :: dl_local_clim_swvl3_filename   ! OBS local mean for dl fire model      <Rackhun Son>
    CHARACTER(len=filename_max) :: dl_local_clim_swvl4_filename   ! OBS local mean for dl fire model      <Rackhun Son>
    CHARACTER(len=filename_max) :: dl_local_clim_lai_filename     ! OBS local mean for dl fire model      <Rackhun Son>


    INTEGER                     :: fire_algorithm
    INTEGER                     :: windbreak_algorithm
    LOGICAL                     :: lburn_pasture
    ! Fire namelist parameters:
    REAL(wp)                    :: fire_litter_threshold,  & ! minimal amount of litter [mol(C)/m^2(grid box)] for fire
                                   fire_rel_hum_threshold, & ! maximal relative humidity for fire
                                   fire_minimum_woody,     & ! minimal fraction of act_fpc of woody PFT to be burned each year
                                   fire_minimum_grass,     & ! minimal fraction of act_fpc of grass PFT to be burned each year
                                   fire_tau_woody,         & ! return period of fire for woody PFT [year] at 0% relative humidity
                                   fire_tau_grass            ! return period of fire for grass PFT [year] at 0% relative humidity
    ! Wind break namelist parameters:
    REAL(wp)                    :: wnd_threshold, & ! factor by which the previous day maximum wind speed
                                                    ! must be larger than the climatological daily maximum
                                                    ! wind speed to allow any wind damage
                                   wnd_damage_scale ! scaling factor for wind damage

    NAMELIST /jsb_disturb_nml/     &
         active,                   &
         ic_filename,              &
         bc_filename,              &
         fire_algorithm,           &
         windbreak_algorithm,      &
         lburn_pasture,            &
         fire_litter_threshold,    &
         fire_rel_hum_threshold,   &
         fire_minimum_woody,       &
         fire_minimum_grass,       &
         fire_tau_woody,           &
         fire_tau_grass,           &
         wnd_threshold,            &
         wnd_damage_scale

    INTEGER :: nml_handler, nml_unit, istat

    ! to read fire model parematers
    TYPE(t_input_file) :: input_file   ! <Rackhun Son>
    REAL(wp), POINTER  :: ptr_1D(:)    ! <Rackhun Son>

    INTEGER :: d_climate_input,    d_climate_hidden       ! <Rackhun Son>
    INTEGER :: d_vegetation_input, d_vegetation_hidden    ! <Rackhun Son>
    INTEGER :: d_lnt_input,        d_lnt_hidden           ! <Rackhun Son>
    INTEGER :: d_emb_input,        d_emb_hidden           ! <Rackhun Son>
    INTEGER :: d_output,           d_emb_output           ! <Rackhun Son>

    REAL(wp), ALLOCATABLE :: climate_i2i_weight(:,:), climate_i2i_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_h2i_weight(:,:), climate_h2i_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_i2f_weight(:,:), climate_i2f_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_h2f_weight(:,:), climate_h2f_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_i2g_weight(:,:), climate_i2g_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_h2g_weight(:,:), climate_h2g_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_i2o_weight(:,:), climate_i2o_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_h2o_weight(:,:), climate_h2o_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_out_weight(:,:), climate_out_bias(:) ! <Rackhun Son>

    REAL(wp), ALLOCATABLE :: vegetation_i2i_weight(:,:), vegetation_i2i_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_h2i_weight(:,:), vegetation_h2i_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_i2f_weight(:,:), vegetation_i2f_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_h2f_weight(:,:), vegetation_h2f_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_i2g_weight(:,:), vegetation_i2g_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_h2g_weight(:,:), vegetation_h2g_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_i2o_weight(:,:), vegetation_i2o_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_h2o_weight(:,:), vegetation_h2o_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_out_weight(:,:), vegetation_out_bias(:) ! <Rackhun Son>

    REAL(wp), ALLOCATABLE :: lightning_topography_layer1_weight(:,:), lightning_topography_layer1_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: lightning_topography_layer2_weight(:,:), lightning_topography_layer2_bias(:) ! <Rackhun Son>

    REAL(wp), ALLOCATABLE :: climate_emb_layer1_weight(:,:), climate_emb_layer1_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: climate_emb_layer2_weight(:,:), climate_emb_layer2_bias(:) ! <Rackhun Son>

    REAL(wp), ALLOCATABLE :: vegetation_emb_layer1_weight(:,:), vegetation_emb_layer1_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_emb_layer2_weight(:,:), vegetation_emb_layer2_bias(:) ! <Rackhun Son>

    REAL(wp), ALLOCATABLE :: climate_min_max(:,:)                 ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vegetation_min_max(:,:)              ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: embedding_min_max(:,:)               ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: lightning_topography_min_max(:,:)    ! <Rackhun Son>
    
    !REAL(wp), ALLOCATABLE :: t_air_minmax(:), t_air_ano_minmax(:), t_air_clim_minmax(:) ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: q_rel_minmax(:), q_rel_ano_minmax(:), q_rel_clim_minmax(:) ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: q_minmax(:),     q_ano_minmax(:),     q_clim_minmax(:)     ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: swl1_minmax(:),  swl1_ano_minmax(:),  swl1_clim_minmax(:)  ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: swl2_minmax(:),  swl2_ano_minmax(:),  swl2_clim_minmax(:)  ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: swl3_minmax(:),  swl3_ano_minmax(:),  swl3_clim_minmax(:)  ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: swl4_minmax(:),  swl4_ano_minmax(:),  swl4_clim_minmax(:)  ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: lai_minmax(:),   lai_ano_minmax(:),   lai_clim_minmax(:)   ! <Rackhun Son>

    !REAL(wp), ALLOCATABLE :: wind_10m_minmax(:)                         ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: rain_minmax(:)                             ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: population_density_minmax(:)               ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: lightning_climatology_minmax(:)            ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: elevation_minmax(:)                        ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: agb_minmax(:)                              ! <Rackhun Son>

    !REAL(wp), ALLOCATABLE :: pft1_minmax(:)                             ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft2_minmax(:)                             ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft3_minmax(:)                             ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft4_minmax(:)                             ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft5_minmax(:)                             ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft6_minmax(:)                             ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft7_minmax(:)                             ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft8_minmax(:)                             ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: pft9_minmax(:)                             ! <Rackhun Son>

    TYPE(t_jsb_vgrid), POINTER :: daily            ! New axis for the day of year <Rackhun Son>
    !TYPE(t_jsb_vgrid), POINTER :: monthly           ! New axis for month           <Rackhun Son>
    TYPE(t_jsb_vgrid), POINTER :: d_climate_lstm    ! New axis for climate LSTM    <Rackhun Son>
    TYPE(t_jsb_vgrid), POINTER :: d_vegetation_lstm ! New axis for vegetation LSTM <Rackhun Son>

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_disturb_config'

    CALL message(TRIM(routine), 'Starting disturb configuration')

    ! Set defaults
    active                 = .TRUE.                      ! activate fire model <Rackhun Son>
    bc_filename            = 'bc_land_disturb.nc'
    ic_filename            = 'ic_land_disturb.nc'
    bc_lightning_filename  = 'bc_land_lightning.nc'      ! lightning frequency (daily) for fire <Rackhun Son>
    bc_topography_filename = 'bc_land_topography.nc'     ! topographical complexity for fire    <Rackhun Son>
    !bc_population_filename = 'bc_land_population.nc'     ! population density (daily) for fire  <Rackhun Son>
    dl_fire_model_params_filename = 'dl_model_params.nc' ! dl fire model parameters             <Rackhun Son>

    dl_local_clim_t_air_filename   = 'dl_local_clim_t2m.nc'   ! t2m local mean for fire           <Rackhun Son>
    dl_local_clim_rH_cf_filename   = 'dl_local_clim_rH_cf.nc' ! rel humidity local mean for fire  <Rackhun Son>
    dl_local_clim_q_filename       = 'dl_local_clim_q.nc'     ! humidity local mean for fire      <Rackhun Son>
    dl_local_clim_ws10_filename    = 'dl_local_clim_ws10.nc'  ! wind speed local mean for fire    <Rackhun Son>
    dl_local_clim_swvl1_filename   = 'dl_local_clim_swvl1.nc' ! swvl1 local mean for fire         <Rackhun Son>
    dl_local_clim_swvl2_filename   = 'dl_local_clim_swvl2.nc' ! swvl2 local mean for fire         <Rackhun Son>
    dl_local_clim_swvl3_filename   = 'dl_local_clim_swvl3.nc' ! swvl3 local mean for fire         <Rackhun Son>
    dl_local_clim_swvl4_filename   = 'dl_local_clim_swvl4.nc' ! swvl4 local mean for fire         <Rackhun Son>
    dl_local_clim_lai_filename     = 'dl_local_clim_lai.nc'   ! lai local mean for fire           <Rackhun Son>

    windbreak_algorithm    =  1
    lburn_pasture          = .FALSE.
    !
    fire_litter_threshold  = 16.67_wp   ! minimal amount of litter [mol(C)/m^2(grid box)] for fire
    fire_rel_hum_threshold = 70._wp     ! maximal relative humidity for fire
    fire_minimum_woody     =  0.002_wp  ! minimal fraction of act_fpc of woody PFT to be burned each year
    fire_minimum_grass     =  0.006_wp  ! minimal fraction of act_fpc of grass PFT to be burned each year
    fire_tau_woody         =  6.0_wp    ! return period of fire for woody PFT [year] at 0% relative humidity
    fire_tau_grass         =  2.0_wp    ! return period of fire for grass PFT [year] at 0% relative humidity
    !
    wnd_threshold          =  2.25_wp   ! factor by which the previous day maximum wind speed
                                        ! must be larger than the climatological daily maximum
                                        ! wind speed to allow any wind damage
    wnd_damage_scale       =  5.e-03_wp ! scaling factor for wind damage

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_disturb_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_disturb_nml)

    CALL close_nml(nml_handler)

    ! Read dl fire model parameters  <Rackhun Son>
    input_file = jsb_netcdf_open_input(dl_fire_model_params_filename)

    ! Read model dimensions  <Rackhun Son>
    ptr_1D => input_file%Read_1d(variable_name='d_climate_input')
    d_climate_input = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='d_climate_hidden')
    d_climate_hidden = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)

    ptr_1D => input_file%Read_1d(variable_name='d_vegetation_input')
    d_vegetation_input = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='d_vegetation_hidden')
    d_vegetation_hidden = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)

    ptr_1D => input_file%Read_1d(variable_name='d_lnt_input')
    d_lnt_input = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='d_lnt_hidden')
    d_lnt_hidden = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)

    ptr_1D => input_file%Read_1d(variable_name='d_output')
    d_output = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)

    ptr_1D => input_file%Read_1d(variable_name='d_emb_input')
    d_emb_input = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='d_emb_hidden')
    d_emb_hidden = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='d_emb_output')
    d_emb_output = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)

    ! Read model parameters  <Rackhun Son>
    ptr_1D => input_file%Read_1d(variable_name='climate_i2i_weight')
    climate_i2i_weight = RESHAPE(ptr_1D, (/d_climate_input, d_climate_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_i2i_bias')
    climate_i2i_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_h2i_weight')
    climate_h2i_weight = RESHAPE(ptr_1D, (/d_climate_hidden, d_climate_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_h2i_bias')
    climate_h2i_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_i2f_weight')
    climate_i2f_weight = RESHAPE(ptr_1D, (/d_climate_input, d_climate_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_i2f_bias')
    climate_i2f_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_h2f_weight')
    climate_h2f_weight = RESHAPE(ptr_1D, (/d_climate_hidden, d_climate_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_h2f_bias')
    climate_h2f_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_i2g_weight')
    climate_i2g_weight = RESHAPE(ptr_1D, (/d_climate_input, d_climate_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_i2g_bias')
    climate_i2g_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_h2g_weight')
    climate_h2g_weight = RESHAPE(ptr_1D, (/d_climate_hidden, d_climate_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_h2g_bias')
    climate_h2g_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_i2o_weight')
    climate_i2o_weight = RESHAPE(ptr_1D, (/d_climate_input, d_climate_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_i2o_bias')
    climate_i2o_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_h2o_weight')
    climate_h2o_weight = RESHAPE(ptr_1D, (/d_climate_hidden, d_climate_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_h2o_bias')
    climate_h2o_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_out_weight')
    climate_out_weight = RESHAPE(ptr_1D, (/d_climate_hidden, d_output*d_emb_output/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_out_bias')
    climate_out_bias = ptr_1D
    DEALLOCATE(ptr_1D)

    ptr_1D => input_file%Read_1d(variable_name='vegetation_i2i_weight')
    vegetation_i2i_weight = RESHAPE(ptr_1D, (/d_vegetation_input, d_vegetation_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_i2i_bias')
    vegetation_i2i_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_h2i_weight')
    vegetation_h2i_weight = RESHAPE(ptr_1D, (/d_vegetation_hidden, d_vegetation_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_h2i_bias')
    vegetation_h2i_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_i2f_weight')
    vegetation_i2f_weight = RESHAPE(ptr_1D, (/d_vegetation_input, d_vegetation_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_i2f_bias')
    vegetation_i2f_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_h2f_weight')
    vegetation_h2f_weight = RESHAPE(ptr_1D, (/d_vegetation_hidden, d_vegetation_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_h2f_bias')
    vegetation_h2f_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_i2g_weight')
    vegetation_i2g_weight = RESHAPE(ptr_1D, (/d_vegetation_input, d_vegetation_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_i2g_bias')
    vegetation_i2g_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_h2g_weight')
    vegetation_h2g_weight = RESHAPE(ptr_1D, (/d_vegetation_hidden, d_vegetation_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_h2g_bias')
    vegetation_h2g_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_i2o_weight')
    vegetation_i2o_weight = RESHAPE(ptr_1D, (/d_vegetation_input, d_vegetation_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_i2o_bias')
    vegetation_i2o_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_h2o_weight')
    vegetation_h2o_weight = RESHAPE(ptr_1D, (/d_vegetation_hidden, d_vegetation_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_h2o_bias')
    vegetation_h2o_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_out_weight')
    vegetation_out_weight = RESHAPE(ptr_1D, (/d_vegetation_hidden, d_output*d_emb_output/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_out_bias')
    vegetation_out_bias = ptr_1D
    DEALLOCATE(ptr_1D)

    ptr_1D => input_file%Read_1d(variable_name='lightning_topography_layer1_weight')
    lightning_topography_layer1_weight = RESHAPE(ptr_1D, (/d_lnt_input, d_lnt_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='lightning_topography_layer1_bias')
    lightning_topography_layer1_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='lightning_topography_layer2_weight')
    lightning_topography_layer2_weight = RESHAPE(ptr_1D, (/d_lnt_hidden, d_output/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='lightning_topography_layer2_bias')
    lightning_topography_layer2_bias = ptr_1D
    DEALLOCATE(ptr_1D)

    ptr_1D => input_file%Read_1d(variable_name='climate_emb_layer1_weight')
    climate_emb_layer1_weight = RESHAPE(ptr_1D, (/d_emb_input, d_emb_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_emb_layer1_bias')
    climate_emb_layer1_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_emb_layer2_weight')
    climate_emb_layer2_weight = RESHAPE(ptr_1D, (/d_emb_hidden, d_emb_output/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='climate_emb_layer2_bias')
    climate_emb_layer2_bias = ptr_1D
    DEALLOCATE(ptr_1D)

    ptr_1D => input_file%Read_1d(variable_name='vegetation_emb_layer1_weight')
    vegetation_emb_layer1_weight = RESHAPE(ptr_1D, (/d_emb_input, d_emb_hidden/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_emb_layer1_bias')
    vegetation_emb_layer1_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_emb_layer2_weight')
    vegetation_emb_layer2_weight = RESHAPE(ptr_1D, (/d_emb_hidden, d_emb_output/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_emb_layer2_bias')
    vegetation_emb_layer2_bias = ptr_1D
    DEALLOCATE(ptr_1D)

    ! Read parameters min/max normalization  <Rackhun Son>
    ptr_1D => input_file%Read_1d(variable_name='climate_min_max')
    climate_min_max = RESHAPE(ptr_1D, (/d_climate_input, 2/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='vegetation_min_max')
    vegetation_min_max = RESHAPE(ptr_1D, (/d_vegetation_input, 2/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='embedding_min_max')
    embedding_min_max = RESHAPE(ptr_1D, (/9, 2/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='lightning_topography_min_max')
    lightning_topography_min_max = RESHAPE(ptr_1D, (/d_lnt_input, 2/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='t2m_min_max')
    !t_air_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='t2m_ano_min_max')
    !t_air_ano_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='t2m_clim_min_max')
    !t_air_clim_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='rh_min_max')
    !q_rel_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='rh_ano_min_max')
    !q_rel_ano_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='rh_clim_min_max')
    !q_rel_clim_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='q_min_max')
    !q_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='q_ano_min_max')
    !q_ano_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='q_clim_min_max')
    !q_clim_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='ws_min_max')
    !wind_10m_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='rain_min_max')
    !rain_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='sw1_min_max')
    !swl1_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='sw1_ano_min_max')
    !swl1_ano_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='sw1_clim_min_max')
    !swl1_clim_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='sw2_min_max')
    !swl2_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='sw2_ano_min_max')
    !swl2_ano_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='sw2_clim_min_max')
    !swl2_clim_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='sw3_min_max')
    !swl3_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='sw3_ano_min_max')
    !swl3_ano_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='sw3_clim_min_max')
    !swl3_clim_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='sw4_min_max')
    !swl4_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='sw4_ano_min_max')
    !swl4_ano_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='sw4_clim_min_max')
    !swl4_clim_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='lai_min_max')
    !lai_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='lai_ano_min_max')
    !lai_ano_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='lai_clim_min_max')
    !lai_clim_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='population_min_max')
    !population_density_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='lightning_min_max')
    !lightning_climatology_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='elevation_min_max')
    !elevation_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='agb_min_max')
    !agb_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)

    !ptr_1D => input_file%Read_1d(variable_name='pft1_min_max')
    !pft1_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='pft2_min_max')
    !pft2_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='pft3_min_max')
    !pft3_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='pft4_min_max')
    !pft4_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='pft5_min_max')
    !pft5_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='pft6_min_max')
    !pft6_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='pft7_min_max')
    !pft7_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='pft8_min_max')
    !pft8_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)
    !ptr_1D => input_file%Read_1d(variable_name='pft9_min_max')
    !pft9_minmax = ptr_1D
    !DEALLOCATE(ptr_1D)


    ! Write resulting values in config memory
    config%active                 = active
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename
    config%bc_lightning_filename  = bc_lightning_filename    ! lightning frequency (daily) for fire       <Rackhun Son>
    config%bc_topography_filename = bc_topography_filename   ! topographical complexity (daily) for fire  <Rackhun Son>
    !config%bc_population_filename = bc_population_filename  ! population density (daily) for fire        <Rackhun Son>
    config%dl_fire_model_params_filename = dl_fire_model_params_filename ! dl fire model parameters       <Rackhun Son>

    config%dl_local_clim_t_air_filename   = dl_local_clim_t_air_filename
    config%dl_local_clim_rH_cf_filename   = dl_local_clim_rH_cf_filename
    config%dl_local_clim_q_filename       = dl_local_clim_q_filename
    config%dl_local_clim_ws10_filename    = dl_local_clim_ws10_filename
    config%dl_local_clim_swvl1_filename   = dl_local_clim_swvl1_filename
    config%dl_local_clim_swvl2_filename   = dl_local_clim_swvl2_filename
    config%dl_local_clim_swvl3_filename   = dl_local_clim_swvl3_filename
    config%dl_local_clim_swvl4_filename   = dl_local_clim_swvl4_filename
    config%dl_local_clim_lai_filename     = dl_local_clim_lai_filename

    config%fire_algorithm         = fire_algorithm
    config%windbreak_algorithm    = windbreak_algorithm
    config%lburn_pasture          = lburn_pasture
    config%fire_litter_threshold  = fire_litter_threshold
    config%fire_rel_hum_threshold = fire_rel_hum_threshold
    config%fire_minimum_woody     = fire_minimum_woody
    config%fire_minimum_grass     = fire_minimum_grass
    config%fire_tau_woody         = fire_tau_woody
    config%fire_tau_grass         = fire_tau_grass
    config%wnd_threshold          = wnd_threshold
    config%wnd_damage_scale       = wnd_damage_scale

    config%d_climate_input          = d_climate_input          ! <Rackhun Son>
    config%d_climate_hidden         = d_climate_hidden         ! <Rackhun Son>
    config%d_vegetation_input       = d_vegetation_input       ! <Rackhun Son>
    config%d_vegetation_hidden      = d_vegetation_hidden      ! <Rackhun Son>
    config%d_lnt_input              = d_lnt_input              ! <Rackhun Son>
    config%d_lnt_hidden             = d_lnt_hidden             ! <Rackhun Son>
    config%d_output                 = d_output                 ! <Rackhun Son>
    config%d_emb_input              = d_emb_input              ! <Rackhun Son>
    config%d_emb_hidden             = d_emb_hidden             ! <Rackhun Son>
    config%d_emb_output             = d_emb_output             ! <Rackhun Son>

    config%climate_i2i_weight       = climate_i2i_weight       ! <Rackhun Son> 
    config%climate_i2i_bias         = climate_i2i_bias         ! <Rackhun Son>
    config%climate_h2i_weight       = climate_h2i_weight       ! <Rackhun Son>
    config%climate_h2i_bias         = climate_h2i_bias         ! <Rackhun Son>
    config%climate_i2f_weight       = climate_i2f_weight       ! <Rackhun Son>
    config%climate_i2f_bias         = climate_i2f_bias         ! <Rackhun Son>
    config%climate_h2f_weight       = climate_h2f_weight       ! <Rackhun Son>
    config%climate_h2f_bias         = climate_h2f_bias         ! <Rackhun Son>
    config%climate_i2g_weight       = climate_i2g_weight       ! <Rackhun Son>
    config%climate_i2g_bias         = climate_i2g_bias         ! <Rackhun Son>
    config%climate_h2g_weight       = climate_h2g_weight       ! <Rackhun Son>
    config%climate_h2g_bias         = climate_h2g_bias         ! <Rackhun Son>
    config%climate_i2o_weight       = climate_i2o_weight       ! <Rackhun Son>
    config%climate_i2o_bias         = climate_i2o_bias         ! <Rackhun Son>
    config%climate_h2o_weight       = climate_h2o_weight       ! <Rackhun Son>
    config%climate_h2o_bias         = climate_h2o_bias         ! <Rackhun Son>
    config%climate_out_weight       = climate_out_weight       ! <Rackhun Son>
    config%climate_out_bias         = climate_out_bias         ! <Rackhun Son>

    config%vegetation_i2i_weight    = vegetation_i2i_weight    ! <Rackhun Son>
    config%vegetation_i2i_bias      = vegetation_i2i_bias      ! <Rackhun Son>
    config%vegetation_h2i_weight    = vegetation_h2i_weight    ! <Rackhun Son>
    config%vegetation_h2i_bias      = vegetation_h2i_bias      ! <Rackhun Son>
    config%vegetation_i2f_weight    = vegetation_i2f_weight    ! <Rackhun Son>
    config%vegetation_i2f_bias      = vegetation_i2f_bias      ! <Rackhun Son>
    config%vegetation_h2f_weight    = vegetation_h2f_weight    ! <Rackhun Son>
    config%vegetation_h2f_bias      = vegetation_h2f_bias      ! <Rackhun Son>
    config%vegetation_i2g_weight    = vegetation_i2g_weight    ! <Rackhun Son>
    config%vegetation_i2g_bias      = vegetation_i2g_bias      ! <Rackhun Son>
    config%vegetation_h2g_weight    = vegetation_h2g_weight    ! <Rackhun Son>
    config%vegetation_h2g_bias      = vegetation_h2g_bias      ! <Rackhun Son>
    config%vegetation_i2o_weight    = vegetation_i2o_weight    ! <Rackhun Son>
    config%vegetation_i2o_bias      = vegetation_i2o_bias      ! <Rackhun Son>
    config%vegetation_h2o_weight    = vegetation_h2o_weight    ! <Rackhun Son>
    config%vegetation_h2o_bias      = vegetation_h2o_bias      ! <Rackhun Son>
    config%vegetation_out_weight    = vegetation_out_weight    ! <Rackhun Son>
    config%vegetation_out_bias      = vegetation_out_bias      ! <Rackhun Son>

    config%lightning_topography_layer1_weight = lightning_topography_layer1_weight   ! <Rackhun Son>
    config%lightning_topography_layer1_bias   = lightning_topography_layer1_bias     ! <Rackhun Son>
    config%lightning_topography_layer2_weight = lightning_topography_layer2_weight   ! <Rackhun Son>
    config%lightning_topography_layer2_bias   = lightning_topography_layer2_bias     ! <Rackhun Son>

    config%climate_emb_layer1_weight    = climate_emb_layer1_weight    ! <Rackhun Son>
    config%climate_emb_layer1_bias      = climate_emb_layer1_bias      ! <Rackhun Son>
    config%climate_emb_layer2_weight    = climate_emb_layer2_weight    ! <Rackhun Son>
    config%climate_emb_layer2_bias      = climate_emb_layer2_bias      ! <Rackhun Son>

    config%vegetation_emb_layer1_weight = vegetation_emb_layer1_weight ! <Rackhun Son>
    config%vegetation_emb_layer1_bias   = vegetation_emb_layer1_bias   ! <Rackhun Son>
    config%vegetation_emb_layer2_weight = vegetation_emb_layer2_weight ! <Rackhun Son>
    config%vegetation_emb_layer2_bias   = vegetation_emb_layer2_bias   ! <Rackhun Son>

    config%climate_min_max              = climate_min_max              ! <Rackhun Son>
    config%vegetation_min_max           = vegetation_min_max           ! <Rackhun Son>
    config%embedding_min_max            = embedding_min_max            ! <Rackhun Son>
    config%lightning_topography_min_max = lightning_topography_min_max ! <Rackhun Son>

    !config%t_air_minmax                 = t_air_minmax
    !config%t_air_ano_minmax             = t_air_ano_minmax
    !config%t_air_clim_minmax            = t_air_clim_minmax
    !config%q_rel_minmax                 = q_rel_minmax
    !config%q_rel_ano_minmax             = q_rel_ano_minmax
    !config%q_rel_clim_minmax            = q_rel_clim_minmax
    !config%q_minmax                     = q_minmax
    !config%q_ano_minmax                 = q_ano_minmax
    !config%q_clim_minmax                = q_clim_minmax
    !config%wind_10m_minmax              = wind_10m_minmax
    !config%rain_minmax                  = rain_minmax
    !config%swl1_minmax                  = swl1_minmax
    !config%swl1_ano_minmax              = swl1_ano_minmax
    !config%swl1_clim_minmax             = swl1_clim_minmax
    !config%swl2_minmax                  = swl2_minmax
    !config%swl2_ano_minmax              = swl2_ano_minmax
    !config%swl2_clim_minmax             = swl2_clim_minmax
    !config%swl3_minmax                  = swl3_minmax
    !config%swl3_ano_minmax              = swl3_ano_minmax
    !config%swl3_clim_minmax             = swl3_clim_minmax
    !config%swl4_minmax                  = swl4_minmax
    !config%swl4_ano_minmax              = swl4_ano_minmax
    !config%swl4_clim_minmax             = swl4_clim_minmax
    !config%lai_minmax                   = lai_minmax
    !config%lai_ano_minmax               = lai_ano_minmax
    !config%lai_clim_minmax              = lai_clim_minmax
    !config%population_density_minmax    = population_density_minmax
    !config%lightning_climatology_minmax = lightning_climatology_minmax
    !config%elevation_minmax             = elevation_minmax
    !config%agb_minmax                   = agb_minmax
    !config%pft1_minmax                  = pft1_minmax
    !config%pft2_minmax                  = pft2_minmax
    !config%pft3_minmax                  = pft3_minmax
    !config%pft4_minmax                  = pft4_minmax
    !config%pft5_minmax                  = pft5_minmax
    !config%pft6_minmax                  = pft6_minmax
    !config%pft7_minmax                  = pft7_minmax
    !config%pft8_minmax                  = pft8_minmax
    !config%pft9_minmax                  = pft9_minmax

    ! New axis for the day of year <Rackhun Son>
    daily => new_vgrid('daily', ZAXIS_GENERIC, 366, units='day')
    !monthly => new_vgrid('monthly', ZAXIS_GENERIC, 12, units='month')
    CALL Register_vgrid(daily)
    !CALL Register_vgrid(monthly)

    ! New axis for climate LSTM <Rackhun Son>
    d_climate_lstm => new_vgrid('d_climate_lstm', ZAXIS_GENERIC, d_climate_hidden, units='')
    CALL Register_vgrid(d_climate_lstm)

    ! New axis for vegetation LSTM <Rackhun Son>
    d_vegetation_lstm => new_vgrid('d_vegetation_lstm', ZAXIS_GENERIC, d_vegetation_hidden, units='')
    CALL Register_vgrid(d_vegetation_lstm)

  END SUBROUTINE Init_disturb_config

#endif
END MODULE mo_disturb_config_class
