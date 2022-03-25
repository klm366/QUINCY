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

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,         ONLY: get_model
#ifdef HAVE_MOFFLEDAC
  USE mo_timer,             ONLY: timer_start, timer_stop
  USE mo_moffledac_iface,   ONLY: timer_moffledac
  USE mo_input,             ONLY: InputVarGroupStart, InputVarGroupEnd, InputFileAdd, InputVarAdd, InputOptAdd, InputGetData, &
                                  InputVarReleaseMem, input_file_list
#else
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,            ONLY: missval
#endif
  USE mo_jsb_impl_constants, ONLY: ifs_nsoil, ifs_soil_depth
  USE mo_util,               ONLY: int2string

  dsl4jsb_Use_processes HYDRO_, PHENO_
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(PHENO_)

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
    !! quincy
    IF (model%config%use_quincy) CALL hydro_init_q_ic(tile)

    ! TODO This is not correct if hydro runs on veg tile!
    IF (tile%Is_last_process_tile(HYDRO_)) THEN
      CALL hydro_finalize_init_vars()
    END IF

  END SUBROUTINE hydro_init

  SUBROUTINE hydro_finalize_init_vars

#ifdef HAVE_MOFFLEDAC
    CALL InputVarReleaseMem('elevation')
    CALL InputVarReleaseMem('orostd')
    CALL InputVarReleaseMem('soil_depth')
    CALL InputVarReleaseMem('maxmoist')
    CALL InputVarReleaseMem('soil_porosity')
    CALL InputVarReleaseMem('hyd_cond_sat')
    CALL InputVarReleaseMem('moisture_pot')
    CALL InputVarReleaseMem('bclapp')
    CALL InputVarReleaseMem('pore_size_index')
    CALL InputVarReleaseMem('soil_field_cap')
    CALL InputVarReleaseMem('wilting_point')
    CALL InputVarReleaseMem('root_depth')
    CALL InputVarReleaseMem('snow')
    CALL InputVarReleaseMem('layer_moist')
! StW: ToDo: Check new ifs-variables with MOFFLEDAC
#else
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

#endif

  END SUBROUTINE hydro_finalize_init_vars

  SUBROUTINE hydro_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(HYDRO_)

    TYPE(t_jsb_model), POINTER :: model

#ifdef HAVE_MOFFLEDAC
    TYPE (input_file_list), POINTER :: BCFile, ICFile, SSOFile
#else

    REAL(wp), POINTER :: &
      & ptr_2D(:,  :),   &
      & ptr_3D(:,:,:)

    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_w
    TYPE(t_input_file) :: input_file
    INTEGER :: i, isoil, nproma, nblks, nsoil
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_read_init_vars'

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (debug_on()) CALL message(routine, 'Reading/setting hydrology init vars')

    model => Get_model(tile%owner_model_id)

#ifndef HAVE_MOFFLEDAC
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
#endif

      ! Get hydrology config
    dsl4jsb_Get_config(HYDRO_)

    !----------
#ifdef HAVE_MOFFLEDAC
    ! Warning: The order of InputVarAdd and InputOptAdd calls below is important: Changing may change the results!
    CALL InputFileAdd(TRIM(dsl4jsb_Config(HYDRO_)%bc_filename),dt_file=0,file_ref=BCFile)
    CALL InputFileAdd(TRIM(dsl4jsb_Config(HYDRO_)%ic_filename),dt_file=0,file_ref=ICFile)
    CALL InputFileAdd(TRIM(dsl4jsb_Config(HYDRO_)%bc_sso_filename),dt_file=0,file_ref=SSOFile)
    CALL InputVarGroupStart('hydro_init',IFile=BCFile,dt_update=0)
    CALL InputVarAdd('nproma,nblock',hydro_init_vars%elevation, 'elevation',IFile=SSOFile)
    CALL InputOptAdd(fill_val=0._wp,valid_range=(/0._wp,1.e300_wp/))
    CALL InputVarAdd('nproma,nblock',hydro_init_vars%oro_stddev,'orostd', IFile=SSOFile)
    CALL InputOptAdd(fill_val=0._wp,valid_range=(/0._wp,1.e300_wp/))
    IF (tile%contains_soil) THEN
      CALL InputVarAdd('nproma,nblock',        hydro_init_vars%soil_depth, 'soil_depth')
      CALL InputOptAdd(fill_val=0.1_wp,  valid_range=(/0.0_wp,1.e300_wp/))
      CALL InputVarAdd('nproma,nblock',        hydro_init_vars%max_moist, 'maxmoist')
      CALL InputOptAdd(fill_val=0.2_wp,  valid_range=(/0.0_wp,1.e300_wp/))
      CALL InputVarAdd('nproma,nblock',        hydro_init_vars%vol_porosity, 'soil_porosity')
      CALL InputOptAdd(fill_val=0.2_wp,  valid_range=(/0.0_wp,1.e300_wp/))
      CALL InputVarAdd('nproma,nblock',        hydro_init_vars%hyd_cond_sat, 'hyd_cond_sat')
      CALL InputOptAdd(fill_val=4.e-6_wp,valid_range=(/0.0_wp,1.e300_wp/))
      CALL InputVarAdd('nproma,nblock',        hydro_init_vars%matrix_pot, 'moisture_pot')
      CALL InputOptAdd(fill_val=0.2_wp,  valid_range=(/0.0_wp,1.e300_wp/))
      CALL InputVarAdd('nproma,nblock',        hydro_init_vars%bclapp, 'bclapp')
      CALL InputOptAdd(fill_val=6.0_wp,  valid_range=(/0.0_wp,1.e300_wp/))
      CALL InputVarAdd('nproma,nblock', hydro_init_vars%pore_size_index,'pore_size_index')
      CALL InputOptAdd(fill_val=0.25_wp, valid_range=(/0.0_wp,1.e300_wp/))
      CALL InputVarAdd('nproma,nblock',        hydro_init_vars%vol_field_cap, 'soil_field_cap')
      CALL InputOptAdd(fill_val=0._wp,   valid_range=(/0._wp,1.e300_wp/))
      CALL InputVarAdd('nproma,nblock',        hydro_init_vars%vol_p_wilt, 'wilting_point')
      CALL InputOptAdd(fill_val=0._wp,   valid_range=(/0._wp,1.e300_wp/))
      CALL InputVarAdd('nproma,nblock',        hydro_init_vars%w_snow_soil, 'snow',       IFile=ICFile)
      CALL InputOptAdd(fill_val=0._wp, valid_range=(/0._wp,1.e300_wp/),def_val=0._wp)
      CALL InputVarAdd('nproma,soillev,nblock',hydro_init_vars%w_soil_sl, 'layer_moist',IFile=ICFIle)
      CALL InputOptAdd(fill_val=0._wp,   valid_range=(/0._wp,1.e300_wp/))
    ENDIF
    IF (tile%contains_vegetation) THEN
      CALL InputVarAdd('nproma,nblock',hydro_init_vars%root_depth,'root_depth')
      CALL InputOptAdd(fill_val=0._wp,valid_range=(/0._wp,1.e300_wp/))
    ENDIF
    CALL InputVarGroupEnd()
    ! Read all hydro variables
    CALL timer_start(timer_moffledac(model%id))
    CALL InputGetData(msg=TRIM(routine))
    CALL timer_stop(timer_moffledac(model%id))

    ! Final adjustments
    IF (tile%contains_soil) THEN
      IF (dsl4jsb_Config(HYDRO_)%w_soil_limit > 0._wp) THEN
        hydro_init_vars%max_moist = MIN(hydro_init_vars%max_moist, dsl4jsb_Config(HYDRO_)%w_soil_limit)
      END IF
    ENDIF
    IF (tile%contains_vegetation) THEN
      hydro_init_vars%root_depth = MIN(hydro_init_vars%root_depth, hydro_init_vars%soil_depth)  ! Temporary fix (root depth sometimes larger then soil depth)
    ENDIF
#else
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
      hydro_init_vars%root_depth = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)
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
#endif

  END SUBROUTINE hydro_read_init_vars

  SUBROUTINE hydro_init_bc(tile)

    USE mo_util,    ONLY: soil_depth_to_layers_2d

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)

    TYPE(t_jsb_vgrid), POINTER :: soil_w
    INTEGER :: i, nsoil

    TYPE(t_jsb_model),       POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init_bc'

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (debug_on()) CALL message(routine, 'Setting  hydrology boundary conditions for tile '//TRIM(tile%name))

    model => Get_model(tile%owner_model_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    ! Get hydrology config
    dsl4jsb_Get_config(HYDRO_)

    ! Get hydrology memory of the tile
    dsl4jsb_Get_memory(HYDRO_)

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
        dsl4jsb_var_ptr(HYDRO_,w_soil_fc_sl) (:,i,:) = hydro_init_vars%vol_field_cap(:,:) * dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      ! Volumetric permanent wilting point
      ! @todo: Do we need this for bare soil tiles?
      dsl4jsb_var_ptr(HYDRO_, vol_p_wilt) = hydro_init_vars%vol_p_wilt

      ! Get permanent wilting point by multiplying with soil layer thicknesses
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,w_soil_pwp_sl) (:,i,:) = hydro_init_vars%vol_p_wilt(:,:) * dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      ! Get water holding capacity (water content at saturation) by multiplying volumetric porosity with soil layer thicknesses
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,w_soil_sat_sl) (:,i,:) = dsl4jsb_var2D_onDomain(HYDRO_,vol_porosity) * dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
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

  END SUBROUTINE hydro_init_bc


  SUBROUTINE hydro_init_ic(tile)

    !-----------------------------------------------------------------------
    !  DECLARATIONS
    USE mo_hydro_constants, ONLY: vol_porosity_org, vol_field_cap_org, vol_p_wilt_org
    USE mo_hydro_util,      ONLY: get_water_in_root_zone
    USE mo_util,            ONLY: ifs2soil

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    dsl4jsb_Def_config(HYDRO_)

    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(PHENO_)

    TYPE(t_jsb_vgrid), POINTER :: soil_w
    INTEGER :: i, nsoil

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init_ic'

    !-----------------------------------------------------------------------
    ! CONTENT

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing hydrology memory for tile '//TRIM(tile%name))

    model => Get_model(tile%owner_model_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    ! Get hydrology config
    dsl4jsb_Get_config(HYDRO_)

    ! Get memory of the tile
    dsl4jsb_Get_memory(HYDRO_)
    IF (tile%contains_vegetation) THEN
      dsl4jsb_Get_memory(PHENO_)
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
    END IF

  END SUBROUTINE hydro_init_ic


  !-----------------------------------------------------------------------------------------------------
  !> Intialize hydro variables initial conditions (ic)
  !! 
  !! quincy
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE hydro_init_q_ic(tile)

    dsl4jsb_Use_memory(HYDRO_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),       POINTER  :: model
    !TYPE(t_input_file)                :: input_file
    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init_q_ic'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    !dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onDomain    :: w_soil_root_theta
    dsl4jsb_Real2D_onDomain    :: w_soil_root_pot
    ! ---------------------------
    ! 0.4 Debug Option
    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of hydro memory (quincy) for tile '// &
      &                          TRIM(tile%name))
    ! ---------------------------
    ! 0.5 Get Memory
    model => get_model(tile%owner_model_id)
    ! Get process config
    !dsl4jsb_Get_config(HYDRO_)
    ! Get process memories
    dsl4jsb_Get_memory(HYDRO_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onDomain(HYDRO_, w_soil_root_theta)
    dsl4jsb_Get_var2D_onDomain(HYDRO_, w_soil_root_pot)
    
    

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 1.0 init constants
    !! 
    !


    !> 2.0 init variables
    !!    
    ! get input file
    !input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HYDRO_)%ic_filename), model%grid_id)

    w_soil_root_theta(:,:) = 1.0_wp   ! this is the default value, maybe better use site-specific value
    w_soil_root_pot(:,:)   = -0.032999999999_wp ! quincy: saxtonA(:,1) * (MIN(w_soil_root_fc,w_soil_root(:)) / root_depth(:)) ** saxtonB(:,1) 


  END SUBROUTINE hydro_init_q_ic


  SUBROUTINE init_soil_properties(tile)

    USE mo_hydro_constants, ONLY: vol_porosity_org, hyd_cond_sat_org, bclapp_org, matrix_pot_org, pore_size_index_org, &
      & vol_field_cap_org, vol_p_wilt_org

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! Local variables
    !
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)

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

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HYDRO_)

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

  END SUBROUTINE init_soil_properties

#endif
END MODULE mo_hydro_init
