!>
!! @brief Initialization of the the hydrology memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-09)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_hydro_init
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message
  USE mo_jsb_control,        ONLY: debug_on

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,             ONLY: missval
  USE mo_jsb_impl_constants, ONLY: ifs_nsoil, ifs_soil_depth
  USE mo_util,               ONLY: int2string

  dsl4jsb_Use_processes HYDRO_, PHENO_, SSE_
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(SSE_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: hydro_init

  TYPE t_hydro_init_vars
    REAL(wp), POINTER ::                &
      & elevation      (:,:) => NULL(), &
      & oro_stddev     (:,:) => NULL(), &
      & soil_depth     (:,:) => NULL(), &
      & max_moist      (:,:) => NULL(), &
      & vol_porosity   (:,:) => NULL(), &
      & hyd_cond_sat   (:,:) => NULL(), &
      & matrix_pot     (:,:) => NULL(), &
      & bclapp         (:,:) => NULL(), &
      & pore_size_index(:,:) => NULL(), &
      & vol_field_cap  (:,:) => NULL(), &
      & vol_p_wilt     (:,:) => NULL(), &
      & root_depth     (:,:) => NULL(), &
      & w_snow_soil    (:,:) => NULL(), &
      & w_soil_sl    (:,:,:) => NULL(), &
      & ifs_smi_sl   (:,:,:) => NULL(), &
      & ifs_w_snow     (:,:) => NULL()
  END TYPE t_hydro_init_vars

  TYPE(t_hydro_init_vars) :: hydro_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_init'

CONTAINS

  !
  !> Intialize hydrology process (after memory has been set up)
  !
  SUBROUTINE hydro_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init'

    model => Get_model(tile%owner_model_id)

    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL hydro_read_init_vars(tile)
    END IF

    CALL hydro_init_bc(tile)
    CALL hydro_init_ic(tile)

    ! TODO This is not correct if hydro runs on veg tile!
    IF (tile%Is_last_process_tile(HYDRO_)) THEN
      CALL hydro_finalize_init_vars()
    END IF

  END SUBROUTINE hydro_init

  SUBROUTINE hydro_finalize_init_vars

    DEALLOCATE( &
      & hydro_init_vars%elevation       ,       &
      & hydro_init_vars%oro_stddev      ,       &
      & hydro_init_vars%soil_depth      ,       &
      & hydro_init_vars%max_moist       ,       &
      & hydro_init_vars%vol_porosity    ,       &
      & hydro_init_vars%hyd_cond_sat    ,       &
      & hydro_init_vars%matrix_pot      ,       &
      & hydro_init_vars%bclapp          ,       &
      & hydro_init_vars%pore_size_index ,       &
      & hydro_init_vars%vol_field_cap   ,       &
      & hydro_init_vars%vol_p_wilt      ,       &
      & hydro_init_vars%root_depth      ,       &
      & hydro_init_vars%w_snow_soil)

    IF (ASSOCIATED(hydro_init_vars%w_soil_sl)) DEALLOCATE(hydro_init_vars%w_soil_sl)
    IF (ASSOCIATED(hydro_init_vars%ifs_smi_sl)) DEALLOCATE(hydro_init_vars%ifs_smi_sl)
    IF (ASSOCIATED(hydro_init_vars%ifs_w_snow)) DEALLOCATE(hydro_init_vars%ifs_w_snow)

  END SUBROUTINE hydro_finalize_init_vars

  SUBROUTINE hydro_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(HYDRO_)

    REAL(wp), POINTER :: &
      & ptr_2D(:,  :),   &
      & ptr_3D(:,:,:)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_w
    TYPE(t_input_file) :: input_file
    INTEGER :: i, isoil, nproma, nblks, nsoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_read_init_vars'

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (debug_on()) CALL message(routine, 'Reading/setting hydrology init vars')

    model => Get_model(tile%owner_model_id)

    grid   => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()
    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    ALLOCATE( &
      & hydro_init_vars%elevation       (nproma, nblks),       &
      & hydro_init_vars%oro_stddev      (nproma, nblks),       &
      & hydro_init_vars%soil_depth      (nproma, nblks),       &
      & hydro_init_vars%max_moist       (nproma, nblks),       &
      & hydro_init_vars%vol_porosity    (nproma, nblks),       &
      & hydro_init_vars%hyd_cond_sat    (nproma, nblks),       &
      & hydro_init_vars%matrix_pot      (nproma, nblks),       &
      & hydro_init_vars%bclapp          (nproma, nblks),       &
      & hydro_init_vars%pore_size_index (nproma, nblks),       &
      & hydro_init_vars%vol_field_cap   (nproma, nblks),       &
      & hydro_init_vars%vol_p_wilt      (nproma, nblks),       &
      & hydro_init_vars%root_depth      (nproma, nblks),       &
      & hydro_init_vars%w_snow_soil     (nproma, nblks)        &
      & )

      hydro_init_vars%elevation      (:,:)   = missval
      hydro_init_vars%oro_stddev     (:,:)   = missval
      hydro_init_vars%soil_depth     (:,:)   = missval
      hydro_init_vars%max_moist      (:,:)   = missval
      hydro_init_vars%vol_porosity   (:,:)   = missval
      hydro_init_vars%hyd_cond_sat   (:,:)   = missval
      hydro_init_vars%matrix_pot     (:,:)   = missval
      hydro_init_vars%bclapp         (:,:)   = missval
      hydro_init_vars%pore_size_index(:,:)   = missval
      hydro_init_vars%vol_field_cap  (:,:)   = missval
      hydro_init_vars%vol_p_wilt     (:,:)   = missval
      hydro_init_vars%root_depth     (:,:)   = missval
      hydro_init_vars%w_snow_soil    (:,:)   = missval

    IF (model%config%init_from_ifs) THEN
      ALLOCATE(hydro_init_vars%ifs_smi_sl(nproma,ifs_nsoil,nblks))
      hydro_init_vars%ifs_smi_sl(:,:,:) = missval
      ALLOCATE(hydro_init_vars%ifs_w_snow(nproma,nblks))
      hydro_init_vars%ifs_w_snow(:,:) = missval
    ELSE
      ALLOCATE(hydro_init_vars%w_soil_sl(nproma, nsoil, nblks))
      hydro_init_vars%w_soil_sl(:,:,:) = missval
    END IF

      ! Get hydrology config
    dsl4jsb_Get_config(HYDRO_)

    !----------

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HYDRO_)%bc_sso_filename), model%grid_id)

    ! Orography ...
    ptr_2D => input_file%Read_2d(  &
      & variable_name='elevation', &
      & fill_array = hydro_init_vars%elevation)
    hydro_init_vars%elevation = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)

    ptr_2D => input_file%Read_2d(  &
      & variable_name='orostd',    &
      & fill_array = hydro_init_vars%oro_stddev)
    hydro_init_vars%oro_stddev = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)

    CALL input_file%Close()

    IF (tile%contains_soil) THEN

      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HYDRO_)%bc_filename), model%grid_id)

      ! Total soil depth until bedrock
      ptr_2D => input_file%Read_2d(   &
        & variable_name='soil_depth', &
        & fill_array = hydro_init_vars%soil_depth)
      hydro_init_vars%soil_depth = MAX(0.1_wp, MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp))

      ! Maximum root zone soil moisture
      ptr_2D => input_file%Read_2d(                        &
        & fill_array = hydro_init_vars%max_moist, &
        & variable_name='maxmoist')
      ! max_moist is zero over glaciers ... but since we do this only for tiles with soil, we set max_moist to 0.2 over glaciers
      hydro_init_vars%max_moist = MERGE(ptr_2D, 0.2_wp, ptr_2D > 0._wp)
      IF (dsl4jsb_Config(HYDRO_)%w_soil_limit > 0._wp) THEN
        hydro_init_vars%max_moist = MIN(hydro_init_vars%max_moist, dsl4jsb_Config(HYDRO_)%w_soil_limit)
      END IF

      ! Volumetric soil porosity
      ptr_2D => input_file%Read_2d(                            &
        & fill_array = hydro_init_vars%vol_porosity,  &
        & variable_name='soil_porosity')
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) >= 0._wp)

      ! Saturated hydraulic conductivity
      ptr_2D => input_file%Read_2d(                 &
        & variable_name='hyd_cond_sat',             &
        & fill_array = hydro_init_vars%hyd_cond_sat)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 4.e-6_wp, ptr_2D(:,:) >= 0._wp)

      ! Matrix potential
      ptr_2D => input_file%Read_2d(                 &
        & variable_name='moisture_pot',             &
        & fill_array = hydro_init_vars%matrix_pot)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) >= 0._wp)

      ! Exponent B in Clapp and Hornberger
      ptr_2D => input_file%Read_2d(                 &
        & variable_name='bclapp',                   &
        & fill_array = hydro_init_vars%bclapp)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 6._wp, ptr_2D(:,:) >= 0._wp)

      ! Pore size distribution index
      ptr_2D => input_file%Read_2d(                 &
        & variable_name='pore_size_index',          &
        & fill_array = hydro_init_vars%pore_size_index)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.25_wp, ptr_2D(:,:) >= 0._wp)

      ! Volumetric field capacity [m/m]
      ptr_2D => input_file%Read_2d(       &
        & variable_name='soil_field_cap', & ! TODO change name in infile to w_soil_fcvol or another name that indicates "volumetric"
        & fill_array = hydro_init_vars%vol_field_cap)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

      ! Volumetric permanent wilting point
      ! @todo: Do we need this for bare soil tiles?
      ptr_2D => input_file%Read_2d(                                 &
        & variable_name='wilting_point',                            &
        & fill_array=hydro_init_vars%vol_p_wilt)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

      CALL input_file%Close()

    END IF

    IF (tile%contains_vegetation) THEN

      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HYDRO_)%bc_filename), model%grid_id)

      ptr_2D => input_file%Read_2d(   &
        & variable_name='root_depth', &
        & fill_array = hydro_init_vars%root_depth)
      hydro_init_vars%root_depth = MAX(0.1_wp, MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp))           ! min root_depth is 0.1 m, consistent with soil_depth
      hydro_init_vars%root_depth = MIN(hydro_init_vars%root_depth, hydro_init_vars%soil_depth)  ! Temporary fix (root depth sometimes larger then soil depth)

      CALL input_file%Close()

    END IF

    IF (tile%contains_soil) THEN

      ! Initialize hydrology variables from ic file

      IF (model%config%init_from_ifs) THEN

        DO isoil=1,ifs_nsoil
          ptr_2D => model%config%ifs_input_file%Read_2d_1lev_1time( &
            & variable_name='SMIL'//int2string(isoil))
          hydro_init_vars%ifs_smi_sl(:,isoil,:) = ptr_2D(:,:)
        END DO
        ptr_3D => model%config%ifs_input_file%Read_2d_time( &
          & variable_name='W_SNOW', start_time_step=1, end_time_step=1)
        hydro_init_vars%ifs_w_snow(:,:) = ptr_3D(:,:,1)

      ELSE

        input_file = jsb_netcdf_open_input(dsl4jsb_Config(HYDRO_)%ic_filename, model%grid_id)

        ! Initial snow depth
        IF (input_file%Has_var('snow')) THEN
          ptr_2D => input_file%Read_2d( &
            & variable_name='snow',     &
            & fill_array = hydro_init_vars%w_snow_soil)
          hydro_init_vars%w_snow_soil = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)
        ELSE
          hydro_init_vars%w_snow_soil = 0._wp
          CALL message(TRIM(routine), 'Initializing snow to zero')
        END IF

        DO isoil=1,nsoil
          ptr_3D => input_file%Read_2d_extdim( &
            & variable_name='layer_moist',     &
            & start_extdim=isoil, end_extdim=isoil, extdim_name='soillev')
          hydro_init_vars%w_soil_sl(:,isoil,:) = MERGE(ptr_3D(:,:,1), 0._wp, ptr_3D(:,:,1) >= 0._wp)
        END DO

        CALL input_file%Close()

      END IF

    END IF

  END SUBROUTINE hydro_read_init_vars

  SUBROUTINE hydro_init_bc(tile)

    USE mo_util,                  ONLY: soil_depth_to_layers_2d
    USE mo_hydro_constants,       ONLY: kdiff_sat_max
    
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)

    TYPE(t_jsb_model), POINTER        :: model
    TYPE(t_jsb_vgrid), POINTER        :: soil_w
    TYPE(t_jsb_vgrid), POINTER        :: vgrid_soil_sb        !< Vertical grid
    INTEGER                           :: i, nsoil
    INTEGER                           :: nsoil_sb             !< number of soil layers as used/defined by the SB_ process
    INTEGER                           :: isoil, iblk, nblks
    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init_bc'

    ! SSE_ 3D
    dsl4jsb_Real3D_onDomain    :: soil_lay_width_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_center_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_ubound_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_lbound_sl

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (debug_on()) CALL message(routine, 'Setting  hydrology boundary conditions for tile '//TRIM(tile%name))

    model => Get_model(tile%owner_model_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels
    vgrid_soil_sb  => Get_vgrid('soil_layer_sb')
    nsoil_sb       =  vgrid_soil_sb%n_levels

    ! Get hydrology config
    dsl4jsb_Get_config(HYDRO_)

    ! Get hydrology memory of the tile
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)

    IF (tile%contains_soil .AND. model%config%use_quincy) THEN
      ! SSE_ 3D
      dsl4jsb_Get_var3D_onDomain(SSE_, soil_lay_width_sl)          ! out
      dsl4jsb_Get_var3D_onDomain(SSE_, soil_lay_depth_center_sl)   ! out
      dsl4jsb_Get_var3D_onDomain(SSE_, soil_lay_depth_ubound_sl)   ! out
      dsl4jsb_Get_var3D_onDomain(SSE_, soil_lay_depth_lbound_sl)   ! out
      nblks = SIZE(soil_lay_width_sl, 3)
    END IF

    !----------

    ! Orography ...
    dsl4jsb_var2D_onDomain(HYDRO_,elevation) = hydro_init_vars%elevation

    dsl4jsb_var2D_onDomain(HYDRO_,oro_stddev) = hydro_init_vars%oro_stddev

    !$acc update device(dsl4jsb_var2D_onDomain(HYDRO_,elevation), &
    !$acc dsl4jsb_var2D_onDomain(HYDRO_,oro_stddev))

    IF (tile%contains_soil) THEN

      ! Initialize hydrology variables from bc file

      ! Total soil depth until bedrock
      dsl4jsb_var2D_onDomain(HYDRO_, soil_depth) = hydro_init_vars%soil_depth

      ! Compute actual soil layer depths based on the fixed layer thicknesses from input file and
      ! the bedrock depth just read in.
      dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl) = soil_depth_to_layers_2d( &
        & dsl4jsb_var2D_onDomain(HYDRO_, soil_depth),                          & ! Total soil depth until bedrock (from textures)
        & soil_w%dz(:))                                                          ! Soil layer thicknesses

      ! Maximum root zone soil moisture
      ! max_moist is zero over glaciers ... but since we do this only for tiles with soil, we set max_moist to 0.2 over glaciers
      dsl4jsb_var2D_onDomain(HYDRO_, max_moist) = hydro_init_vars%max_moist

      ! Volumetric soil porosity
      dsl4jsb_var_ptr(HYDRO_, vol_porosity) = hydro_init_vars%vol_porosity

      ! Saturated hydraulic conductivity
      dsl4jsb_var_ptr(HYDRO_, hyd_cond_sat) = hydro_init_vars%hyd_cond_sat

      ! Matrix potential
      dsl4jsb_var_ptr(HYDRO_, matrix_pot) = hydro_init_vars%matrix_pot

      ! Exponent B in Clapp and Hornberger
      dsl4jsb_var_ptr(HYDRO_, bclapp) = hydro_init_vars%bclapp

      ! Pore size distribution index
      dsl4jsb_var_ptr(HYDRO_, pore_size_index) = hydro_init_vars%pore_size_index

      ! Volumetric field capacity [m/m]
      dsl4jsb_var_ptr(HYDRO_, vol_field_cap) = hydro_init_vars%vol_field_cap

      ! Get water content at field capacity by multiplying volumetric field capacity with soil layer thicknesses
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,w_soil_fc_sl) (:,i,:) = hydro_init_vars%vol_field_cap(:,:) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      ! Volumetric permanent wilting point
      ! @todo: Do we need this for bare soil tiles?
      dsl4jsb_var_ptr(HYDRO_, vol_p_wilt) = hydro_init_vars%vol_p_wilt

      ! Get permanent wilting point by multiplying with soil layer thicknesses
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,w_soil_pwp_sl) (:,i,:) = hydro_init_vars%vol_p_wilt(:,:) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      ! Get water holding capacity (water content at saturation) by multiplying volumetric porosity with soil layer thicknesses
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,w_soil_sat_sl) (:,i,:) = dsl4jsb_var2D_onDomain(HYDRO_,vol_porosity) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      !$acc update device(dsl4jsb_var_ptr(HYDRO_,w_soil_sat_sl), &
      !$acc dsl4jsb_var_ptr(HYDRO_,w_soil_pwp_sl), &
      !$acc dsl4jsb_var_ptr(HYDRO_,w_soil_fc_sl), &
      !$acc dsl4jsb_var_ptr(HYDRO_, vol_field_cap), &
      !$acc dsl4jsb_var_ptr(HYDRO_, vol_p_wilt), &
      !$acc dsl4jsb_var_ptr(HYDRO_, pore_size_index), &
      !$acc dsl4jsb_var_ptr(HYDRO_, bclapp), &
      !$acc dsl4jsb_var_ptr(HYDRO_, matrix_pot), &
      !$acc dsl4jsb_var_ptr(HYDRO_, hyd_cond_sat), &
      !$acc dsl4jsb_var_ptr(HYDRO_, vol_porosity), &
      !$acc dsl4jsb_var2D_onDomain(HYDRO_, max_moist), &
      !$acc dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl), &
      !$acc dsl4jsb_var2D_onDomain(HYDRO_, soil_depth) )

    END IF

    IF (tile%contains_vegetation) THEN

      dsl4jsb_var_ptr(HYDRO_, root_depth) = hydro_init_vars%root_depth
      dsl4jsb_var3D_onDomain(HYDRO_, root_depth_sl) = soil_depth_to_layers_2d( &
        & dsl4jsb_var2D_onDomain(HYDRO_, root_depth),                          & ! Total rooting depth
        & soil_w%dz(:))    ! Soil layer thicknesses

      !$acc update device(dsl4jsb_var2D_onDomain(HYDRO_, root_depth), &
      !$acc dsl4jsb_var3D_onDomain(HYDRO_, root_depth_sl))
    END IF

    ! quincy
    IF (tile%contains_soil .AND. model%config%use_quincy) THEN
      ! pass soil-layer values to a SSE_ variable with a more intuitive variable-name
      ! soil_depth_sl is the width of each soil layer until bedrock
      soil_lay_width_sl(:,:,:) = dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)

      ! calc more metrics of the soil layers
      !   i)  lower & upper boundary
      !   ii) depth at the center of each layer
      DO iblk = 1,nblks
        DO isoil = 1,nsoil_sb
          ! lower & upper boundary
          IF (isoil == 1) THEN
            soil_lay_depth_lbound_sl(:,isoil,iblk) = 0.0_wp
            soil_lay_depth_ubound_sl(:,isoil,iblk) = soil_lay_width_sl(:,isoil,iblk)
          ELSE
            soil_lay_depth_lbound_sl(:,isoil,iblk) = soil_lay_depth_lbound_sl(:,isoil-1,iblk) + &
              & soil_lay_width_sl(:,isoil-1,iblk)
            soil_lay_depth_ubound_sl(:,isoil,iblk) = soil_lay_depth_ubound_sl(:,isoil-1,iblk) + &
              & soil_lay_width_sl(:,isoil,iblk)
          END IF                                
          ! soil-layer center ['(lower_bound + upper_bound) * 0.5']
          soil_lay_depth_center_sl(:,isoil,iblk)   = (soil_lay_depth_lbound_sl(:,isoil,iblk) + &
            & soil_lay_depth_ubound_sl(:,isoil,iblk)) * 0.5_wp 
        END DO
      END DO
    END IF  ! tile%contains_soil & quincy

  END SUBROUTINE hydro_init_bc


  SUBROUTINE hydro_init_ic(tile)

    !-----------------------------------------------------------------------
    !  DECLARATIONS
    USE mo_hydro_constants,     ONLY: vol_porosity_org, vol_field_cap_org, vol_p_wilt_org
    USE mo_hydro_util,          ONLY: get_water_in_root_zone
    USE mo_jsb_math_constants,  ONLY: eps8
    USE mo_util,                ONLY: ifs2soil

    !-----------------------------------------------------------------------
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    !-----------------------------------------------------------------------
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(PHENO_)

    TYPE(t_jsb_model), POINTER              :: model
    TYPE(t_jsb_vgrid), POINTER              :: soil_w
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: hlp1
    INTEGER                                 :: i, nsoil, nproma, nblks
    INTEGER                                 :: j,k
    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init_ic'
    !-----------------------------------------------------------------------
    ! HYDRO_ 2D
    dsl4jsb_Real2D_onDomain    :: root_depth
    dsl4jsb_Real2D_onDomain    :: w_soil_root
    dsl4jsb_Real2D_onDomain    :: w_soil_root_fc
    dsl4jsb_Real2D_onDomain    :: w_soil_root_theta
    dsl4jsb_Real2D_onDomain    :: w_soil_root_pot
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onDomain    :: saxtonA 
    dsl4jsb_Real3D_onDomain    :: saxtonB 
    dsl4jsb_Real3D_onDomain    :: soil_depth_sl
    dsl4jsb_Real3D_onDomain    :: w_soil_sl
    dsl4jsb_Real3D_onDomain    :: w_soil_fc_sl
    dsl4jsb_Real3D_onDomain    :: w_soil_pot_sl

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing hydrology memory for tile '//TRIM(tile%name))

    model => Get_model(tile%owner_model_id)
    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    ! Get hydrology config
    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SSE_)
    ! Get memory of the tile
    dsl4jsb_Get_memory(HYDRO_)
    IF (tile%contains_vegetation) THEN
      dsl4jsb_Get_memory(PHENO_)
    END IF

    IF (tile%contains_vegetation .AND. model%config%use_quincy) THEN
      ! HYDRO_ 2D
      dsl4jsb_Get_var2D_onDomain(HYDRO_,  root_depth)          ! in
      dsl4jsb_Get_var2D_onDomain(HYDRO_,  w_soil_root)         ! in
      dsl4jsb_Get_var2D_onDomain(HYDRO_,  w_soil_root_fc)      ! in
      dsl4jsb_Get_var2D_onDomain(HYDRO_,  w_soil_root_theta)   ! out
      dsl4jsb_Get_var2D_onDomain(HYDRO_,  w_soil_root_pot)     ! out
      ! HYDRO_ 3D
      dsl4jsb_Get_var3D_onDomain(HYDRO_,  saxtonA)             ! calculated in init_soil_properties(), used as input here
      dsl4jsb_Get_var3D_onDomain(HYDRO_,  saxtonB)             ! calculated in init_soil_properties(), used as input here
      dsl4jsb_Get_var3D_onDomain(HYDRO_,  soil_depth_sl)       ! in
      dsl4jsb_Get_var3D_onDomain(HYDRO_,  w_soil_sl)           ! in
      dsl4jsb_Get_var3D_onDomain(HYDRO_,  w_soil_fc_sl)        ! in
      dsl4jsb_Get_var3D_onDomain(HYDRO_,  w_soil_pot_sl)       ! out
      ! allocate local var
      nproma = SIZE(saxtonA, 1)
      nblks  = SIZE(saxtonA, 3)
      ALLOCATE(hlp1(nproma, nsoil, nblks))
    END IF

    !----------

    IF (tile%contains_soil) THEN

      IF (model%config%init_from_ifs) THEN

        CALL ifs2soil(hydro_init_vars%ifs_smi_sl(:,:,:), ifs_soil_depth(:), dsl4jsb_var_ptr(HYDRO_,w_soil_sl)(:,:,:), soil_w%ubounds(:))

        ! Convert soil moisture index (SMI) in IFS file to water content
        WHERE (dsl4jsb_var_ptr(HYDRO_,w_soil_sl)(:,:,:) < -999._wp)
          ! Fill missing values
          dsl4jsb_var_ptr(HYDRO_,w_soil_sl)(:,:,:) = &
            &  0.5_wp * (dsl4jsb_var_ptr(HYDRO_,w_soil_pwp_sl)(:,:,:) + dsl4jsb_var_ptr(HYDRO_,w_soil_fc_sl) (:,:,:))
        ELSE WHERE
          dsl4jsb_var_ptr(HYDRO_,w_soil_sl)(:,:,:) = &
            & MIN(dsl4jsb_var_ptr(HYDRO_,w_soil_sat_sl)(:,:,:), &
            &     MAX(0._wp, &
            &         dsl4jsb_var_ptr(HYDRO_,w_soil_sl)(:,:,:) &
            &         * (dsl4jsb_var_ptr(HYDRO_,w_soil_fc_sl) (:,:,:) - dsl4jsb_var_ptr(HYDRO_,w_soil_pwp_sl)(:,:,:)) &
            &         + dsl4jsb_var_ptr(HYDRO_,w_soil_pwp_sl)(:,:,:) &
            &        ) &
            &    )
        END WHERE
        dsl4jsb_var_ptr(HYDRO_, w_snow_soil)(:,:) = hydro_init_vars%ifs_w_snow(:,:)
        dsl4jsb_var2D_onDomain(HYDRO_, w_snow) = dsl4jsb_var2D_onDomain(HYDRO_, w_snow_soil)

      ELSE

        ! Initialize hydrology variables from ic file

        ! Initial snow depth
        dsl4jsb_var_ptr(HYDRO_, w_snow_soil) = hydro_init_vars%w_snow_soil
        dsl4jsb_var2D_onDomain(HYDRO_, w_snow) = dsl4jsb_var2D_onDomain(HYDRO_, w_snow_soil)

        dsl4jsb_var_ptr(HYDRO_,w_soil_sl) (:,:,:) = hydro_init_vars%w_soil_sl(:,:,:)

      END IF

      ! Initial w_soil_overflown water
      dsl4jsb_var2D_onDomain(HYDRO_,w_soil_overflow) = 0._wp

      !$acc update device(dsl4jsb_var2D_onDomain(HYDRO_, w_snow),       &
      !$acc               dsl4jsb_var2D_onDomain(HYDRO_, w_snow_soil),  &
      !$acc               dsl4jsb_var2D_onDomain(HYDRO_, w_soil_overflow))
      ! avoid compiler warnings about pointers not being dereferenced

    END IF

    ! Initialize organic soil layer fractions and modify saturated, field capacity and perm. wilting point capacities accordingly
    ! Note: at the moment, fract_org_sl doesn't change over time, but it could
    IF (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil) THEN
      dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl) = 0._wp
      IF (tile%is_vegetation) THEN
        ! For general VEG tile
        dsl4jsb_var_ptr(HYDRO_, fract_org_sl)(:,1,:) = dsl4jsb_var2D_onDomain(PHENO_, fract_forest)
        IF (tile%lcts(1)%lib_id /= 0)  THEN ! only if the present tile is a PFT (otherwise lctlib not defined)
          ! For PFT tiles
          ! Note: in current setup hydrology does not run on PFTs
          IF (dsl4jsb_Lctlib_param(ForestFlag)) THEN
            dsl4jsb_var_ptr(HYDRO_, fract_org_sl)(:,1,:) = dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max)
          END IF
        END IF
      END IF
      ! The following can change later depending on the organic layers and soil ice content.
      dsl4jsb_var3D_onDomain(HYDRO_, w_soil_sat_sl) =                                                        &
        & (1 - dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl)) * dsl4jsb_var3D_onDomain(HYDRO_, w_soil_sat_sl) &
        & +    dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl)  * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl) * vol_porosity_org
      dsl4jsb_var3D_onDomain(HYDRO_, w_soil_fc_sl) =                                                         &
        & (1 - dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl)) * dsl4jsb_var3D_onDomain(HYDRO_, w_soil_fc_sl)  &
        & +    dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl)  * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl) * vol_field_cap_org
      dsl4jsb_var3D_onDomain(HYDRO_, w_soil_pwp_sl) =                                                        &
        & (1 - dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl)) * dsl4jsb_var3D_onDomain(HYDRO_, w_soil_pwp_sl) &
        & +    dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl)  * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl) * vol_p_wilt_org
      
      !$acc update device(dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl), &
      !$acc dsl4jsb_var3D_onDomain(HYDRO_, w_soil_sat_sl), &
      !$acc dsl4jsb_var3D_onDomain(HYDRO_, w_soil_fc_sl), &
      !$acc dsl4jsb_var3D_onDomain(HYDRO_, w_soil_pwp_sl))
    END IF

    IF (tile%contains_soil) THEN
      ! @todo: The interpolation of the initial layer moisture from T63 to the ICON grid leads to positive values where the
      ! water content should actually be zero. Therefore, make sure that water content is not larger than field capacity.
      dsl4jsb_var3D_onDomain(HYDRO_, w_soil_sl) = &
        & MIN(dsl4jsb_var3D_onDomain(HYDRO_, w_soil_sl), dsl4jsb_var3D_onDomain(HYDRO_, w_soil_fc_sl))
      !$acc update device(dsl4jsb_var3D_onDomain(HYDRO_, w_soil_sl))

      DO i=1,SIZE(dsl4jsb_var_ptr(HYDRO_, w_soil_root), 2)
          CALL get_water_in_root_zone(                       &
          & dsl4jsb_var_ptr(HYDRO_, w_soil_sl)     (:,:,i),  &
          & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl) (:,:,i),  &
          & dsl4jsb_var_ptr(HYDRO_, root_depth_sl) (:,:,i),  &
          & dsl4jsb_var_ptr(HYDRO_, w_soil_column) (:,  i)   )
      END DO
      !$acc update host(dsl4jsb_var_ptr(HYDRO_, w_soil_column))
      WHERE (dsl4jsb_var2D_onDomain(HYDRO_, w_soil_column) > dsl4jsb_var2D_onDomain(HYDRO_, max_moist))
        dsl4jsb_var2D_onDomain(HYDRO_, w_soil_column) = dsl4jsb_var2D_onDomain(HYDRO_, max_moist)
      END WHERE

      WHERE (dsl4jsb_var2D_onDomain(HYDRO_, max_moist) > 0._wp)
        dsl4jsb_var2D_onDomain(HYDRO_, w_soil_rel) = &
          & dsl4jsb_var2D_onDomain(HYDRO_, w_soil_column) / dsl4jsb_var2D_onDomain(HYDRO_, max_moist)
      ELSE WHERE
        dsl4jsb_var2D_onDomain(HYDRO_, w_soil_rel) = 0._wp
      END WHERE

      
      !$acc update device(dsl4jsb_var2D_onDomain(HYDRO_, w_soil_rel),    &
      !$acc               dsl4jsb_var2D_onDomain(HYDRO_, w_soil_column), &
      !$acc               dsl4jsb_var3D_onDomain(HYDRO_, w_soil_sl))

      CALL init_soil_properties(tile)
    END IF

    ! if vegetation (assuming there is also soil)
    IF (tile%contains_vegetation) THEN
      ! Calculate water content and field capacity in the root zone
      DO i=1,SIZE(dsl4jsb_var_ptr(HYDRO_, w_soil_root), 2)
          CALL get_water_in_root_zone( &
          & dsl4jsb_var_ptr(HYDRO_, w_soil_sl)      (:,:,i), & 
          & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, root_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, w_soil_root)    (:,  i)  )
 
          CALL get_water_in_root_zone( & 
          & dsl4jsb_var_ptr(HYDRO_, w_soil_fc_sl)   (:,:,i), & 
          & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, root_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, w_soil_root_fc) (:,  i)  )
          
          CALL get_water_in_root_zone( & 
          & dsl4jsb_var_ptr(HYDRO_, w_soil_pwp_sl)  (:,:,i), & 
          & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, root_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, w_soil_root_pwp)(:,  i)  )
      END DO
      !$acc update host(dsl4jsb_var_ptr(HYDRO_, w_soil_root),     &
      !$acc             dsl4jsb_var_ptr(HYDRO_, w_soil_root_fc),  &
      !$acc             dsl4jsb_var_ptr(HYDRO_, w_soil_root_pwp))
      
      ! quincy
      IF (model%config%use_quincy) THEN
        ! adjust initial water amount to prescribed theta
        ! w_soil_root(:)        = w_soil_root_pwp(:) + & 
        !                         dsl4jsb_Config(SSE_)%soil_theta_prescribe * (w_soil_root_fc(:) - w_soil_root_pwp(:)) 
        w_soil_root_theta(:,:) = dsl4jsb_Config(SSE_)%soil_theta_prescribe
        ! in case there is no water in the root zone
        WHERE (MIN(w_soil_root_fc(:,:), w_soil_root(:,:)) > eps8)
          w_soil_root_pot(:,:) = saxtonA(:,1,:) * &
            & (MIN(w_soil_root_fc(:,:), w_soil_root(:,:)) / root_depth(:,:)) ** saxtonB(:,1,:)
        ELSEWHERE
          w_soil_root_pot(:,:) = 0.0_wp
        ENd WHERE
        
        !> 7.2 soil water potential per layer
        !! 
        ! in case a soil-layer width (soil_depth_sl) is zero AND w_soil_fc_sl & w_soil_sl is zero (hlp1 should not become zero)
        ! (it seem that these two are always zero at the same time, i.e., not useing .OR.)
        WHERE((soil_depth_sl(:,:,:) > eps8) .AND. (MIN(w_soil_fc_sl(:,:,:), w_soil_sl(:,:,:)) > eps8))
          hlp1(:,:,:)          = MIN(w_soil_fc_sl(:,:,:), w_soil_sl(:,:,:)) / soil_depth_sl(:,:,:)
          w_soil_pot_sl(:,:,:) = saxtonA(:,:,:) * hlp1(:,:,:) ** saxtonB(:,:,:)
        ELSEWHERE
          w_soil_pot_sl(:,:,:) = 0.0_wp 
        END WHERE
      END IF ! quincy
    END IF ! vegetation

  END SUBROUTINE hydro_init_ic


  SUBROUTINE init_soil_properties(tile)

    USE mo_hydro_constants, ONLY: vol_porosity_org, hyd_cond_sat_org, bclapp_org, matrix_pot_org, pore_size_index_org, &
      & vol_field_cap_org, vol_p_wilt_org
    USE mo_hydro_constants, ONLY: k_pwp_s, k_pwp_c, k_pwp_sc, k_pwp_a, k_pwp_at, k_pwp_bt, &
      & k_fc_s, k_fc_c, k_fc_sc, k_fc_a, k_fc_at, k_fc_bt, k_fc_ct, &
      & k_sat_s, k_sat_c, k_sat_sc, k_sat_a, k_sat_at, k_sat_bt, k_sat_ct, k_sat_dt, &
      & wpot_pwp, wpot_fc, wpot_fc

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! Local variables
    !
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    INTEGER                                          :: nproma, nblks
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)          :: hlp1
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)          :: theta_pwp_sl, theta_fc_sl

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: &
      & hyd_cond_sat,          &
      & vol_porosity,          &
      & bclapp,                &
      & matrix_pot,            &
      & pore_size_index,       &
      & vol_field_cap,         &
      & vol_p_wilt
    dsl4jsb_Real3D_onDomain :: &
      & fract_org_sl,          &
      & hyd_cond_sat_sl,       &
      & vol_porosity_sl,       &
      & bclapp_sl,             &
      & matrix_pot_sl,         &
      & pore_size_index_sl,    &
      & vol_field_cap_sl,      &
      & vol_p_wilt_sl

    ! SSE_ 3D
    dsl4jsb_Real3D_onDomain    :: sand_sl
    dsl4jsb_Real3D_onDomain    :: silt_sl
    dsl4jsb_Real3D_onDomain    :: clay_sl
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onDomain    :: saxtonA 
    dsl4jsb_Real3D_onDomain    :: saxtonB 
    dsl4jsb_Real3D_onDomain    :: saxtonC 

    ! Locally allocated vectors
    !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_w

    INTEGER  :: nsoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_soil_properties'

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (tile%is_lake) RETURN

    model => Get_model(tile%owner_model_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SSE_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)

    dsl4jsb_Get_var2D_onDomain(HYDRO_, hyd_cond_sat)
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_porosity)
    dsl4jsb_Get_var2D_onDomain(HYDRO_, bclapp)
    dsl4jsb_Get_var2D_onDomain(HYDRO_, matrix_pot)
    dsl4jsb_Get_var2D_onDomain(HYDRO_, pore_size_index)
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_field_cap)
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_p_wilt)
    IF (dsl4jsb_Config(HYDRO_)%l_organic) dsl4jsb_Get_var3D_onDomain(HYDRO_, fract_org_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, hyd_cond_sat_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, bclapp_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, matrix_pot_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, pore_size_index_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_field_cap_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_p_wilt_sl)

    IF (model%config%use_quincy) THEN
      ! SSE_ 3D
      dsl4jsb_Get_var3D_onDomain(SSE_, sand_sl)       ! out  <- run init here (instead of SSE_) to resolve dependency problem (HYDRO_ init needs to be done before SSE_ init)
      dsl4jsb_Get_var3D_onDomain(SSE_, silt_sl)       ! out  <- run init here (instead of SSE_) to resolve dependency problem (HYDRO_ init needs to be done before SSE_ init)
      dsl4jsb_Get_var3D_onDomain(SSE_, clay_sl)       ! out  <- run init here (instead of SSE_) to resolve dependency problem (HYDRO_ init needs to be done before SSE_ init)
      ! HYDRO_ 3D
      dsl4jsb_Get_var3D_onDomain(HYDRO_, saxtonA)     ! out
      dsl4jsb_Get_var3D_onDomain(HYDRO_, saxtonB)     ! out
      dsl4jsb_Get_var3D_onDomain(HYDRO_, saxtonC)     ! out
      ! ----------------------------------------------------------------------------------------------------- !
      nproma = SIZE(sand_sl,1)
      nblks  = SIZE(sand_sl,3)
      ! ----------------------------------------------------------------------------------------------------- !
      ALLOCATE(hlp1(nproma, nsoil, nblks))
      ALLOCATE(theta_pwp_sl(nproma, nsoil, nblks))
      ALLOCATE(theta_fc_sl(nproma, nsoil, nblks))
    END IF


    hyd_cond_sat_sl   (:,:,:) = SPREAD(hyd_cond_sat   (:,:), DIM=2, ncopies=nsoil)
    vol_porosity_sl   (:,:,:) = SPREAD(vol_porosity   (:,:), DIM=2, ncopies=nsoil)
    bclapp_sl         (:,:,:) = SPREAD(bclapp         (:,:), DIM=2, ncopies=nsoil)
    matrix_pot_sl     (:,:,:) = SPREAD(matrix_pot     (:,:), DIM=2, ncopies=nsoil)
    pore_size_index_sl(:,:,:) = SPREAD(pore_size_index(:,:), DIM=2, ncopies=nsoil)
    vol_field_cap_sl  (:,:,:) = SPREAD(vol_field_cap  (:,:), DIM=2, ncopies=nsoil)
    vol_p_wilt_sl     (:,:,:) = SPREAD(vol_p_wilt     (:,:), DIM=2, ncopies=nsoil)
    IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
      hyd_cond_sat_sl   (:,:,:) = (1 - fract_org_sl(:,:,:)) * hyd_cond_sat_sl   (:,:,:) + fract_org_sl(:,:,:) * hyd_cond_sat_org
      vol_porosity_sl   (:,:,:) = (1 - fract_org_sl(:,:,:)) * vol_porosity_sl   (:,:,:) + fract_org_sl(:,:,:) * vol_porosity_org
      bclapp_sl         (:,:,:) = (1 - fract_org_sl(:,:,:)) * bclapp_sl         (:,:,:) + fract_org_sl(:,:,:) * bclapp_org
      matrix_pot_sl     (:,:,:) = (1 - fract_org_sl(:,:,:)) * matrix_pot_sl     (:,:,:) + fract_org_sl(:,:,:) * matrix_pot_org
      pore_size_index_sl(:,:,:) = (1 - fract_org_sl(:,:,:)) * pore_size_index_sl(:,:,:) + fract_org_sl(:,:,:) * pore_size_index_org
      vol_field_cap_sl  (:,:,:) = (1 - fract_org_sl(:,:,:)) * vol_field_cap_sl  (:,:,:) + fract_org_sl(:,:,:) * vol_field_cap_org
      vol_p_wilt_sl     (:,:,:) = (1 - fract_org_sl(:,:,:)) * vol_p_wilt_sl     (:,:,:) + fract_org_sl(:,:,:) * vol_p_wilt_org
    END IF

    !$acc update device(hyd_cond_sat_sl, vol_porosity_sl, bclapp_sl, &
    !$acc matrix_pot_sl, pore_size_index_sl, vol_field_cap_sl, &
    !$acc vol_p_wilt_sl)

    ! quincy
    IF (model%config%use_quincy) THEN
      ! init soil properties
      sand_sl(:,:,:)      = dsl4jsb_Config(SSE_)%soil_sand
      silt_sl(:,:,:)      = dsl4jsb_Config(SSE_)%soil_silt
      clay_sl(:,:,:)      = dsl4jsb_Config(SSE_)%soil_clay
      ! soil hydraulic properties, based on Saxton and Rawls 2006, SSSA
      ! fractional water holding capacity at permanent wilting point
      !
      ! note: this assumes that the soil hydraulic properties are identical between the layers
      hlp1(:,:,:)         = k_pwp_s * sand_sl(:,:,:) + & 
                            k_pwp_c * clay_sl(:,:,:) + &
                            k_pwp_sc * sand_sl(:,:,:) * clay_sl(:,:,:) + &
                            k_pwp_a
      theta_pwp_sl(:,:,:) = k_pwp_at + (1.0_wp + k_pwp_bt) * hlp1(:,:,:) 

      !>  fractional water holding capacity at field capacity 
      !>
      hlp1(:,:,:)         = k_fc_s * sand_sl(:,:,:) + &
                            k_fc_c * clay_sl(:,:,:) + &
                            k_fc_sc * sand_sl(:,:,:) * clay_sl(:,:,:) + &
                            k_fc_a
      theta_fc_sl(:,:,:)  = k_fc_at + (1.0_wp + k_fc_bt) * hlp1(:,:,:) + k_fc_ct * hlp1(:,:,:) ** 2._wp
      
      !>parameters describing the moisture-tension curve                 
      !>
      ! note that for convenience, B has a negative sign and A is in MPa 
      saxtonB(:,:,:)      = -(LOG(wpot_pwp) - LOG(wpot_fc)) / (LOG(theta_fc_sl(:,:,:)) - LOG(theta_pwp_sl(:,:,:)))
      saxtonA(:,:,:)      = -(EXP(LOG(wpot_fc) - saxtonB(:,:,:) * LOG(theta_fc_sl(:,:,:)) )) / 1000.0_wp

      !>parameters describing the moisture-diffusivity curve
      !>
      ! kdiff_sat_sl & k_diff are not used in jsbach4/quincy
      ! kdiff_sat_sl(:,:,:) = kdiff_sat_max * (theta_sat_sl(:,:,:) - theta_fc_sl(:,:,:)) ** (3.0_wp + 1._wp / saxtonB(:,:,:))
      saxtonC(:,:,:)      = 3.0_wp - 2._wp * saxtonB(:,:,:)
    END IF

  END SUBROUTINE init_soil_properties

#endif
END MODULE mo_hydro_init
