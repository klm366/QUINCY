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

    INTEGER           :: d_input,                & ! <Rackhun Son>
                         d_hidden,               & ! <Rackhun Son>
                         d_output,               & ! <Rackhun Son>
                         d_region,               & ! <Rackhun Son>
                         d_class,                & ! <Rackhun Son>
                         d_hidden_n_class,       & ! <Rackhun Son>
                         d_output_n_class        

    REAL(wp), ALLOCATABLE :: layer1_weight(:,:), layer1_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer2_weight(:,:), layer2_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer3_weight(:,:), layer3_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer4_weight(:,:), layer4_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer5_weight(:,:), layer5_bias(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer6_weight(:,:), layer6_bias(:) ! dl fire model parameters <Rackhun Son>

    REAL(wp), ALLOCATABLE :: layer1_weight_update(:,:), layer1_bias_update(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer2_weight_update(:,:), layer2_bias_update(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer3_weight_update(:,:), layer3_bias_update(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer4_weight_update(:,:), layer4_bias_update(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer5_weight_update(:,:), layer5_bias_update(:) ! dl fire model parameters <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer6_weight_update(:,:), layer6_bias_update(:) ! dl fire model parameters <Rackhun Son>

    REAL(wp), ALLOCATABLE :: mW1(:,:),           mb1(:)         ! dl fire model 1st momentum <Rackhun Son>
    REAL(wp), ALLOCATABLE :: mW2(:,:),           mb2(:)         ! dl fire model 1st momentum <Rackhun Son>
    REAL(wp), ALLOCATABLE :: mW3(:,:),           mb3(:)         ! dl fire model 1st momentum <Rackhun Son>
    REAL(wp), ALLOCATABLE :: mW4(:,:),           mb4(:)         ! dl fire model 1st momentum <Rackhun Son>
    REAL(wp), ALLOCATABLE :: mW5(:,:),           mb5(:)         ! dl fire model 1st momentum <Rackhun Son>
    REAL(wp), ALLOCATABLE :: mW6(:,:),           mb6(:)         ! dl fire model 1st momentum <Rackhun Son>

    REAL(wp), ALLOCATABLE :: vW1(:,:),           vb1(:)         ! dl fire model 2nd momentum <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vW2(:,:),           vb2(:)         ! dl fire model 2nd momentum <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vW3(:,:),           vb3(:)         ! dl fire model 2nd momentum <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vW4(:,:),           vb4(:)         ! dl fire model 2nd momentum <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vW5(:,:),           vb5(:)         ! dl fire model 2nd momentum <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vW6(:,:),           vb6(:)         ! dl fire model 2nd momentum <Rackhun Son>

    REAL(wp), ALLOCATABLE :: min_max(:,:)                  ! dl fire model parameters <Rackhun Son>

  CHARACTER(len=filename_max) :: bc_lightning_filename          ! lightning frequency (daily) for fire <Rackhun Son>
  CHARACTER(len=filename_max) :: bc_topography_filename         ! topographical complexity for fire    <Rackhun Son>
  CHARACTER(len=filename_max) :: bc_anthropogenic_filename      ! anthropogenic factors for fire       <Rackhun Son>
  CHARACTER(len=filename_max) :: bc_gfed_region_filename        ! GFED4 reference regions for fire     <Rackhun Son>
  CHARACTER(len=filename_max) :: dl_fire_model_params_filename  ! dl fire model parameters             <Rackhun Son>


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
    CHARACTER(len=filename_max) :: bc_lightning_filename         ! lightning frequency (daily) for fire   <Rackhun Son>
    CHARACTER(len=filename_max) :: bc_topography_filename        ! topographical complexity for fire      <Rackhun Son>
    CHARACTER(len=filename_max) :: bc_anthropogenic_filename     ! anthropogenic factor for fire          <Rackhun Son>
    CHARACTER(len=filename_max) :: bc_gfed_region_filename       ! GFED4 reference regions  for fire      <Rackhun Son>
    CHARACTER(len=filename_max) :: dl_fire_model_params_filename ! dl fire model parameters               <Rackhun Son>

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

    INTEGER :: d_input,    d_hidden    ! <Rackhun Son>
    INTEGER :: d_output,   d_region    ! <Rackhun Son>
    INTEGER :: d_hidden_n_class        ! <Rackhun Son>
    INTEGER :: d_output_n_class        ! <Rackhun Son>
    INTEGER :: d_class                 ! <Rackhun Son>

    REAL(wp), ALLOCATABLE :: layer1_weight(:,:), layer1_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer2_weight(:,:), layer2_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer3_weight(:,:), layer3_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer4_weight(:,:), layer4_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer5_weight(:,:), layer5_bias(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer6_weight(:,:), layer6_bias(:) ! <Rackhun Son>

    REAL(wp), ALLOCATABLE :: layer1_weight_update(:,:), layer1_bias_update(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer2_weight_update(:,:), layer2_bias_update(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer3_weight_update(:,:), layer3_bias_update(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer4_weight_update(:,:), layer4_bias_update(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer5_weight_update(:,:), layer5_bias_update(:) ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: layer6_weight_update(:,:), layer6_bias_update(:) ! <Rackhun Son>

    REAL(wp), ALLOCATABLE :: mW1(:,:),           mb1(:)         ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: mW2(:,:),           mb2(:)         ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: mW3(:,:),           mb3(:)         ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: mW4(:,:),           mb4(:)         ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: mW5(:,:),           mb5(:)         ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: mW6(:,:),           mb6(:)         ! <Rackhun Son>

    REAL(wp), ALLOCATABLE :: vW1(:,:),           vb1(:)         ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vW2(:,:),           vb2(:)         ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vW3(:,:),           vb3(:)         ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vW4(:,:),           vb4(:)         ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vW5(:,:),           vb5(:)         ! <Rackhun Son>
    REAL(wp), ALLOCATABLE :: vW6(:,:),           vb6(:)         ! <Rackhun Son>

    !REAL(wp), ALLOCATABLE :: i2i_weight(:,:), i2i_bias(:) ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: h2i_weight(:,:), h2i_bias(:) ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: i2f_weight(:,:), i2f_bias(:) ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: h2f_weight(:,:), h2f_bias(:) ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: i2g_weight(:,:), i2g_bias(:) ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: h2g_weight(:,:), h2g_bias(:) ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: i2o_weight(:,:), i2o_bias(:) ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: h2o_weight(:,:), h2o_bias(:) ! <Rackhun Son>
    !REAL(wp), ALLOCATABLE :: out_weight(:,:), out_bias(:) ! <Rackhun Son>

    REAL(wp), ALLOCATABLE :: min_max(:,:)                 ! <Rackhun Son>
    !REAL(wp)              :: training_memory(35,1000000)  ! <Rackhun Son>
    
    TYPE(t_jsb_vgrid), POINTER :: daily            ! New axis for the day of year <Rackhun Son>
    TYPE(t_jsb_vgrid), POINTER :: memory           ! New axis for the day of year <Rackhun Son>

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_disturb_config'

    CALL message(TRIM(routine), 'Starting disturb configuration')

    ! Set defaults
    active                 = .TRUE.                      ! activate fire model <Rackhun Son>
    bc_filename            = 'bc_land_disturb.nc'
    ic_filename            = 'ic_land_disturb.nc'
    bc_lightning_filename  = 'bc_land_lightning.nc'      ! lightning frequency (daily) for fire <Rackhun Son>
    bc_topography_filename = 'bc_land_topography.nc'     ! topographical complexity for fire    <Rackhun Son>
    bc_anthropogenic_filename     = 'bc_anthropogenic.nc'! anthropogenic factor for fire        <Rackhun Son>
    bc_gfed_region_filename       = 'bc_land_gfed.nc'    ! anthropogenic factor for fire        <Rackhun Son>
    dl_fire_model_params_filename = 'dl_model_params.nc' ! dl fire model parameters             <Rackhun Son>

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
    ptr_1D => input_file%Read_1d(variable_name='d_input')
    d_input = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='d_hidden')
    d_hidden = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='d_output')
    d_output = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='d_region')
    d_region = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='d_class')
    d_class = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='d_hidden_n_class')
    d_hidden_n_class = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='d_output_n_class')
    d_output_n_class = SIZE(ptr_1D)
    DEALLOCATE(ptr_1D)


    ! Read model parameters  <Rackhun Son>
    ptr_1D => input_file%Read_1d(variable_name='layer1_weight')
    layer1_weight = RESHAPE(ptr_1D, (/d_input, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer1_bias')
    layer1_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer2_weight')
    layer2_weight = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer2_bias')
    layer2_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer3_weight')
    layer3_weight = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer3_bias')
    layer3_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer4_weight')
    layer4_weight = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer4_bias')
    layer4_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer5_weight')
    layer5_weight = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer5_bias')
    layer5_bias = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer6_weight')
    layer6_weight = RESHAPE(ptr_1D, (/d_hidden, d_output*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer6_bias')
    layer6_bias = ptr_1D
    DEALLOCATE(ptr_1D)

    ptr_1D => input_file%Read_1d(variable_name='layer1_weight_update')
    layer1_weight_update = RESHAPE(ptr_1D, (/d_input, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer1_bias_update')
    layer1_bias_update = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer2_weight_update')
    layer2_weight_update = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer2_bias_update')
    layer2_bias_update = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer3_weight_update')
    layer3_weight_update = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer3_bias_update')
    layer3_bias_update = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer4_weight_update')
    layer4_weight_update = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer4_bias_update')
    layer4_bias_update = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer5_weight_update')
    layer5_weight_update = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer5_bias_update')
    layer5_bias_update = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer6_weight_update')
    layer6_weight_update = RESHAPE(ptr_1D, (/d_hidden, d_output*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='layer6_bias_update')
    layer6_bias_update = ptr_1D
    DEALLOCATE(ptr_1D)

    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_w_init')
    mW1 = RESHAPE(ptr_1D, (/d_input, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_b_init')
    mb1 = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_w')
    mW2 = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_b')
    mb2 = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_w')
    mW3 = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_b')
    mb3 = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_w')
    mW4 = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_b')
    mb4 = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_w')
    mW5 = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_b')
    mb5 = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_w_last')
    mW6 = RESHAPE(ptr_1D, (/d_hidden, d_output*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_b_last')
    mb6 = ptr_1D
    DEALLOCATE(ptr_1D)

    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_w_init')
    vW1 = RESHAPE(ptr_1D, (/d_input, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_b_init')
    vb1 = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_w')
    vW2 = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_b')
    vb2 = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_w')
    vW3 = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_b')
    vb3 = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_w')
    vW4 = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_b')
    vb4 = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_w')
    vW5 = RESHAPE(ptr_1D, (/d_hidden, d_hidden*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_b')
    vb5 = ptr_1D
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_w_last')
    vW6 = RESHAPE(ptr_1D, (/d_hidden, d_output*d_class/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)
    ptr_1D => input_file%Read_1d(variable_name='zero_momentum_b_last')
    vb6 = ptr_1D
    DEALLOCATE(ptr_1D)


    ! Read parameters min/max normalization  <Rackhun Son>
    ptr_1D => input_file%Read_1d(variable_name='min_max')
    min_max = RESHAPE(ptr_1D, (/d_input, 2/), order=(/2, 1/))
    DEALLOCATE(ptr_1D)

    !training_memory(:,:) = 0.
    
    ! Write resulting values in config memory
    config%active                 = active
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename
    config%bc_lightning_filename  = bc_lightning_filename                ! lightning frequency (daily) for fire      <Rackhun Son>
    config%bc_topography_filename = bc_topography_filename               ! topographical complexity (daily) for fire <Rackhun Son>
    config%bc_anthropogenic_filename = bc_anthropogenic_filename         ! anthropogenic factors for fire            <Rackhun Son>
    config%bc_gfed_region_filename= bc_gfed_region_filename              ! GFED4 reference region for fire           <Rackhun Son>
    config%dl_fire_model_params_filename = dl_fire_model_params_filename ! dl fire model parameters                  <Rackhun Son>

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

    config%d_input                = d_input          ! <Rackhun Son>
    config%d_hidden               = d_hidden         ! <Rackhun Son>
    config%d_output               = d_output         ! <Rackhun Son>
    config%d_region               = d_region         ! <Rackhun Son>
    config%d_class                = d_class          ! <Rackhun Son>
    config%d_hidden_n_class       = d_hidden_n_class ! <Rackhun Son>
    config%d_output_n_class       = d_output_n_class ! <Rackhun Son>

    config%layer1_weight          = layer1_weight    ! <Rackhun Son> 
    config%layer1_bias            = layer1_bias      ! <Rackhun Son>
    config%layer2_weight          = layer2_weight    ! <Rackhun Son> 
    config%layer2_bias            = layer2_bias      ! <Rackhun Son>
    config%layer3_weight          = layer3_weight    ! <Rackhun Son> 
    config%layer3_bias            = layer3_bias      ! <Rackhun Son>
    config%layer4_weight          = layer4_weight    ! <Rackhun Son> 
    config%layer4_bias            = layer4_bias      ! <Rackhun Son>
    config%layer5_weight          = layer5_weight    ! <Rackhun Son> 
    config%layer5_bias            = layer5_bias      ! <Rackhun Son>
    config%layer6_weight          = layer6_weight    ! <Rackhun Son> 
    config%layer6_bias            = layer6_bias      ! <Rackhun Son>

    config%layer1_weight_update   = layer1_weight_update    ! <Rackhun Son> 
    config%layer1_bias_update     = layer1_bias_update      ! <Rackhun Son>
    config%layer2_weight_update   = layer2_weight_update    ! <Rackhun Son> 
    config%layer2_bias_update     = layer2_bias_update      ! <Rackhun Son>
    config%layer3_weight_update   = layer3_weight_update    ! <Rackhun Son> 
    config%layer3_bias_update     = layer3_bias_update      ! <Rackhun Son>
    config%layer4_weight_update   = layer4_weight_update    ! <Rackhun Son> 
    config%layer4_bias_update     = layer4_bias_update      ! <Rackhun Son>
    config%layer5_weight_update   = layer5_weight_update    ! <Rackhun Son> 
    config%layer5_bias_update     = layer5_bias_update      ! <Rackhun Son>
    config%layer6_weight_update   = layer6_weight_update    ! <Rackhun Son> 
    config%layer6_bias_update     = layer6_bias_update      ! <Rackhun Son>

    config%mW1                    = mW1              ! <Rackhun Son> 
    config%mb1                    = mb1              ! <Rackhun Son>
    config%mW2                    = mW2              ! <Rackhun Son> 
    config%mb2                    = mb2              ! <Rackhun Son>
    config%mW3                    = mW3              ! <Rackhun Son> 
    config%mb3                    = mb3              ! <Rackhun Son>
    config%mW4                    = mW4              ! <Rackhun Son> 
    config%mb4                    = mb4              ! <Rackhun Son>
    config%mW5                    = mW5              ! <Rackhun Son> 
    config%mb5                    = mb5              ! <Rackhun Son>
    config%mW6                    = mW6              ! <Rackhun Son> 
    config%mb6                    = mb6              ! <Rackhun Son>

    config%vW1                    = vW1              ! <Rackhun Son> 
    config%vb1                    = vb1              ! <Rackhun Son>
    config%vW2                    = vW2              ! <Rackhun Son> 
    config%vb2                    = vb2              ! <Rackhun Son>
    config%vW3                    = vW3              ! <Rackhun Son> 
    config%vb3                    = vb3              ! <Rackhun Son>
    config%vW4                    = vW4              ! <Rackhun Son> 
    config%vb4                    = vb4              ! <Rackhun Son>
    config%vW5                    = vW5              ! <Rackhun Son> 
    config%vb5                    = vb5              ! <Rackhun Son>
    config%vW6                    = vW6              ! <Rackhun Son> 
    config%vb6                    = vb6              ! <Rackhun Son>

    config%min_max                = min_max          ! <Rackhun Son>


    ! New axis for the day of year <Rackhun Son>
    daily => new_vgrid('daily', ZAXIS_GENERIC, 366, units='day')
    CALL Register_vgrid(daily)
    memory => new_vgrid('memory', ZAXIS_GENERIC, 730, units='day')
    CALL Register_vgrid(memory)

  END SUBROUTINE Init_disturb_config

#endif
END MODULE mo_disturb_config_class
