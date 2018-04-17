!>
!! @brief Initialization of the the soil and snow memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_sse_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_grid_class,    ONLY: t_jsb_vgrid
  USE mo_jsb_grid,          ONLY: Get_vgrid
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
#ifdef HAVE_MOFFLEDAC
  USE mo_timer,             ONLY: timer_start, timer_stop, timer_input
  USE mo_input,             ONLY: InputVarGroupStart, InputVarGroupEnd, InputFileAdd, InputVarAdd, InputOptAdd, InputDimAdd, &
                                  InputGetData, InputVarReleaseMem, input_file_list, input_var_list
#endif

  dsl4jsb_Use_processes SSE_, SEB_, HYDRO_
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sse_init

  REAL(wp), PARAMETER :: &
    & fao_thermal_diff(0:6) = (/7.4e-7_wp ,8.7e-7_wp ,8.0e-7_wp ,7.4e-7_wp ,7.1e-7_wp ,6.7e-7_wp ,-1._wp/), &
    & fao_vol_hcap    (0:6) = (/2.25e+6_wp,1.93e+6_wp,2.10e+6_wp,2.25e+6_wp,2.36e+6_wp,2.48e+6_wp,-1._wp/)

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_init'

CONTAINS

  !
  !> Intialize soil and snow energy process (after memory has been set up)
  !
  SUBROUTINE sse_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sse_init'

    CALL sse_init_bc(tile)
    CALL sse_init_ic(tile)

  END SUBROUTINE sse_init

  SUBROUTINE sse_init_bc(tile)
#ifdef HAVE_MOFFLEDAC
    TYPE (input_file_list), POINTER :: IFile
    TYPE (input_var_list),  POINTER :: vfao
#endif

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)

    REAL(wp), POINTER :: &
      & fao    (:,  :), & !< FAO soil types
      & ptr_2D (:,  :)    !< temporary pointer

    INTEGER :: iblk

    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sse_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get soil config
    dsl4jsb_Get_config(SSE_)

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting  soil energy boundary conditions for tile '//TRIM(tile%name))

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SSE_)

#ifdef HAVE_MOFFLEDAC
    CALL InputFileAdd(TRIM(dsl4jsb_Config(SSE_)%bc_filename),file_ref=IFile,dt_file=0)
    CALL InputVarGroupStart('sse_'//TRIM(tile%name))
    CALL InputVarAdd('nproma,nblock',fao,'fao_'//TRIM(tile%name),var_ref=vfao,dt_update=0)
    CALL InputVarGroupEnd()
    CALL InputOptAdd(fill_val=0._wp,valid_range=(/REAL(LBOUND(fao_thermal_diff,1),wp),REAL(UBOUND(fao_thermal_diff,1),wp)/), &
                     AVar=vfao,AFile=IFile,var_file_name='fao')
    IF (dsl4jsb_Config(SSE_)%l_heat_cap_map) THEN
      CALL InputOptAdd(var_name='sse_vol_heat_cap_'//TRIM(tile%name),AFile=IFile,var_file_name='heat_capacity', &
                       fill_val=2.e6_wp,valid_range=(/0._wp,1.e300_wp/),missing_val=0._wp)
    ENDIF
    IF (dsl4jsb_Config(SSE_)%l_heat_cond_map) THEN
      CALL InputOptAdd(var_name='sse_heat_cond_'//TRIM(tile%name),AFile=IFile,var_file_name='heat_conductivity', &
                       fill_val=0.2_wp,valid_range=(/0._wp,1.e300_wp/),missing_val=0._wp)
    ENDIF
    CALL timer_start(timer_input)
    CALL InputGetData(IFile=IFile,msg=TRIM(routine))
    CALL timer_stop(timer_input)
#else
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SSE_)%bc_filename), model%grid_id)

    ! Thermal diffusivity of soil
    fao => input_file%Read_2d(                    &
      & variable_name='fao')
    fao(:,:) = MERGE(fao(:,:), 0._wp, fao(:,:) >= 0._wp)
    fao(:,:) = MAX(REAL(LBOUND(fao_thermal_diff,1),wp), MIN(REAL(UBOUND(fao_thermal_diff,1),wp), fao(:,:)))

    IF (dsl4jsb_Config(SSE_)%l_heat_cap_map) THEN
      ptr_2D => input_file%Read_2d(         &
        & variable_name='heat_capacity',            &
        & fill_array = dsl4jsb_var_ptr(SSE_, vol_heat_cap))
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 2.e6_wp, ptr_2D(:,:) > 0._wp)
    ENDIF

    IF (dsl4jsb_Config(SSE_)%l_heat_cond_map) THEN
      ptr_2D => input_file%Read_2d(         &
        & variable_name='heat_conductivity',         &
        & fill_array = dsl4jsb_var_ptr(SSE_, heat_cond))
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) > 0._wp)
    ENDIF 

    CALL input_file%Close()                              ! Temporary
#endif

    DO iblk=1,SIZE(fao,2)
      dsl4jsb_var_ptr(SSE_,thermal_diffusivity)(:,iblk) = fao_thermal_diff(NINT(fao(:,iblk)))
    END DO

    ! Volumetric heat capacity of dry soil [J/(m^3*K)] if not read
    IF (.NOT. dsl4jsb_Config(SSE_)%l_heat_cap_map) THEN
      DO iblk=1,SIZE(fao,2)
        dsl4jsb_var_ptr(SSE_,vol_heat_cap)(:,iblk) = fao_vol_hcap(NINT(fao(:,iblk)))
      END DO
    END IF

    ! Heat conductivity of mineral soil if not read
    IF (.NOT. dsl4jsb_Config(SSE_)%l_heat_cond_map) THEN
      dsl4jsb_var2D_onDomain(SSE_, heat_cond) = &
        & dsl4jsb_var2D_onDomain(SSE_, vol_heat_cap) * dsl4jsb_var2D_onDomain(SSE_, thermal_diffusivity)
    END IF

#ifdef HAVE_MOFFLEDAC
    ! Clean up temporary variables
    CALL InputVarReleaseMem(var_ref=vfao)
#endif

  END SUBROUTINE sse_init_bc

  SUBROUTINE sse_init_ic(tile)

    USE mo_sse_process, ONLY: init_soil_temperature

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

#ifdef HAVE_MOFFLEDAC
    TYPE (input_file_list), POINTER :: IFile
    TYPE (input_var_list),  POINTER :: var
#endif

    TYPE(t_jsb_model),         POINTER :: model
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(SEB_)

    REAL(wp), POINTER :: &
      & t_srf_mon (:,:,:),    & !< Climatological monthly values of surface temperature
      & ptr_2D    (:,  :),    & !< temporary pointer
      & ptr_3D    (:,:,:)

    TYPE(t_input_file) :: input_file
    TYPE(t_jsb_vgrid), POINTER :: soil_e
    INTEGER :: iblk

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sse_init_ic'

    model => get_model(tile%owner_model_id)

    soil_e => Get_vgrid('soil_depth_energy')

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(SEB_)

    ! Get soil config
    dsl4jsb_Get_config(SSE_)

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing soil energy memory for tile '//TRIM(tile%name))

#ifdef HAVE_MOFFLEDAC
    CALL InputDimAdd('clim_monthly_year',Extent=12)
    CALL InputFileAdd(TRIM(dsl4jsb_Config(SSE_)%ic_filename),file_ref=IFile,dt_file=0)
    CALL InputVarGroupStart('sse_'//TRIM(tile%name))
    CALL InputVarAdd('nproma,clim_monthly_year,nblock',t_srf_mon,'surf_temp_'//TRIM(tile%name),IFile,var_ref=var,dt_update=0)
    CALL InputVarGroupEnd()
    CALL InputOptAdd(AVar=var,var_file_name='surf_temp',fill_val=280._wp,valid_range=(/0._wp,1.e300_wp/),missing_val=0._wp)
    CALL InputOptAdd(AFile=IFile,dim_name='time=>clim_monthly_year')
    CALL InputGetData(var=var,msg=TRIM(routine))
    DO iblk=1,SIZE(t_srf_mon, DIM=3)
      CALL init_soil_temperature(                     &
        & soil_e%levels(:),                           &
        & t_srf_mon(:,1:12,iblk),                     &
        & dsl4jsb_var_ptr(SSE_,t_soil)(:,:,iblk),   &
        & dsl4jsb_var_ptr(SEB_,t)(:,iblk))
    END DO

    ! Clean up temporary variables
    CALL InputVarReleaseMem(var_ref=var)
#else
    input_file = jsb_netcdf_open_input(dsl4jsb_Config(SSE_)%ic_filename, model%grid_id)

    ! Soil and surface temperatures
    !
    ! Read climatological surface temperature
    t_srf_mon => input_file%Read_2d_time(         &   ! Temporary
      & variable_name='surf_temp')
    t_srf_mon(:,:,:) = MERGE(t_srf_mon(:,:,:), 280._wp, t_srf_mon(:,:,:) > 0._wp)

    CALL input_file%Close()

    DO iblk=1,SIZE(t_srf_mon, DIM=2)
      CALL init_soil_temperature(                     &
        & soil_e%levels(:),                           &
        & t_srf_mon(:,iblk,1:12),                     &
        & dsl4jsb_var_ptr(SSE_,t_soil)(:,:,iblk),   &
        & dsl4jsb_var_ptr(SEB_,t)(:,iblk))
    END DO
#endif
    dsl4jsb_var2D_onDomain(SEB_,t_unfilt) = dsl4jsb_var2D_onDomain(SEB_,t)

    CALL init_soil_properties(tile)

  END SUBROUTINE sse_init_ic

  SUBROUTINE init_soil_properties(tile)

    USE mo_jsb_physical_constants, ONLY: ci, dens_snow_min
    USE mo_sse_constants,          ONLY: vol_hcap_ice, vol_hcap_snow, vol_hcap_snow_min, &
      &                                  hcond_ice, hcond_org, hcond_snow, hcond_snow_min
    USE mo_sse_process,            ONLY: calc_vol_heat_capacity, calc_thermal_conductivity

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! Local variables
    !
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: &
      & vol_heat_cap, &
      & heat_cond, &
      & vol_heat_cap_snow, &
      & heat_cond_snow
    dsl4jsb_Real3D_onDomain :: &
      & soil_depth_sl, &
      & fract_org_sl, &
      & vol_heat_cap_sl, &
      & heat_cond_sl

    ! Locally allocated vectors
    !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_e

    LOGICAL  :: l_compat401
    INTEGER  :: nsoil, nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_soil_properties'

    model => Get_model(tile%owner_model_id)

    l_compat401 = model%config%l_compat401

    soil_e => Get_vgrid('soil_depth_energy')
    nsoil = soil_e%n_levels

    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(HYDRO_)

    ! Get reference to variables for domain
    !
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var3D_onDomain(SSE_, vol_heat_cap_sl)
    dsl4jsb_Get_var3D_onDomain(SSE_, heat_cond_sl)
    dsl4jsb_Get_var2D_onDomain(SSE_, vol_heat_cap_snow)
    dsl4jsb_Get_var2D_onDomain(SSE_, heat_cond_snow)

    ! On glacier tile, set variables and return
    IF (tile%is_glacier) THEN
      vol_heat_cap_sl  (:,:,:) = vol_hcap_ice
      heat_cond_sl     (:,:,:) = hcond_ice
      vol_heat_cap_snow(:  ,:) = vol_hcap_snow
      heat_cond_snow   (:  ,:) = hcond_snow
      RETURN
    END IF

    dsl4jsb_Get_var3D_onDomain(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_var2D_onDomain(SSE_,   vol_heat_cap)
    dsl4jsb_Get_var2D_onDomain(SSE_,   heat_cond)
    nproma = SIZE(heat_cond,1)
    nblks  = SIZE(heat_cond,2)
    IF (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil) THEN
      dsl4jsb_Get_var3D_onDomain(HYDRO_, fract_org_sl)
    ELSE
      ALLOCATE(fract_org_sl(nproma, nsoil, nblks))
      fract_org_sl(:,:,:) = 0._wp
    END IF

    IF (dsl4jsb_Config(SSE_)%l_heat_cap_dyn) THEN
      vol_heat_cap_sl (:,:,:) = calc_vol_heat_capacity( &
        &                              SPREAD( vol_heat_cap(:,:),                 DIM=2, ncopies=nsoil)               , &
        &       dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl)                                                       , &
        &                                      fract_org_sl(:,:,:)                                                    , &
        & MERGE(dsl4jsb_var3D_onDomain(HYDRO_, w_soil_sl) / soil_depth_sl(:,:,:), 0._wp, soil_depth_sl(:,:,:) > 0._wp), &
        & MERGE(dsl4jsb_var3D_onDomain(HYDRO_, w_ice_sl ) / soil_depth_sl(:,:,:), 0._wp, soil_depth_sl(:,:,:) > 0._wp), &
        &                                      soil_depth_sl(:,:,:)                                                   , &
        &                       SPREAD(SPREAD( soil_e%dz(:),                       DIM=1, ncopies=nproma), DIM=3, ncopies=nblks)  &
        & )
    ELSE
      vol_heat_cap_sl(:,:,:) = SPREAD( vol_heat_cap(:,:), DIM=2, ncopies=nsoil)
    END IF

    IF (dsl4jsb_Config(SSE_)%l_heat_cond_dyn) THEN
      heat_cond_sl (:,:,:) = calc_thermal_conductivity( &
        &                              SPREAD( heat_cond(:,:),                    DIM=2, ncopies=nsoil)               , &
        &       dsl4jsb_var3D_onDomain(HYDRO_, w_soil_sl)                                                             , &
        &       dsl4jsb_var3D_onDomain(HYDRO_, w_ice_sl)                                                              , &
        &       dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl)                                                       , &
        &       dsl4jsb_var3D_onDomain(HYDRO_, vol_field_cap_sl)                                                      , &
        &                                      fract_org_sl(:,:,:)                                                    , &
        &                                      soil_depth_sl(:,:,:)                                                   , &
        &                       SPREAD(SPREAD( soil_e%dz(:),                       DIM=1, ncopies=nproma), DIM=3, ncopies=nblks)  &
        & )
    ELSE
      heat_cond_sl(:,:,:) = SPREAD(heat_cond(:,:), DIM=2, ncopies=nsoil)
    END IF

    IF (.NOT. (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil)) THEN
      DEALLOCATE(fract_org_sl)
    END IF

    IF (dsl4jsb_Config(SSE_)%l_snow) THEN
      IF (model%config%l_compat401) THEN
        vol_heat_cap_snow(:,:) = 634500._wp
        heat_cond_snow   (:,:) = 0.31_wp
      ELSE
        IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
          ! USe snow density values for fresh snow for initialization
          vol_heat_cap_snow(:,:) = ci * dens_snow_min                               ! Use specific heat of ice (Verseghy 1991)
          heat_cond_snow   (:,:) = 2.5e-6_wp*dens_snow_min**2._wp - 0.123e-3_wp*dens_snow_min + 0.024_wp ! Calonne et al. (2011)
          ! vol_heat_cap_snow(:,:) = vol_hcap_snow_min
          ! heat_cond_snow   (:,:) = hcond_snow_min
        ELSE
          vol_heat_cap_snow(:,:) = vol_hcap_snow
          heat_cond_snow   (:,:) = hcond_snow
        END IF
      END IF
    ! ELSE
    !   vol_heat_cap_snow and heat_cond_snow are fixed in calc_soiltemp_old
    END IF

  END SUBROUTINE init_soil_properties

#endif
END MODULE mo_sse_init
