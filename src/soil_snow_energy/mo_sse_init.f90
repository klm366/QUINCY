!>
!! @brief Initialization of the the soil and snow memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_sse_init
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

  dsl4jsb_Use_processes SSE_, SEB_, HYDRO_
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sse_init, init_soil_properties

  REAL(wp), PARAMETER :: &
    & fao_thermal_diff(0:6) = (/7.4e-7_wp ,8.7e-7_wp ,8.0e-7_wp ,7.4e-7_wp ,7.1e-7_wp ,6.7e-7_wp ,-1._wp/), &
    & fao_vol_hcap    (0:6) = (/2.25e+6_wp,1.93e+6_wp,2.10e+6_wp,2.25e+6_wp,2.36e+6_wp,2.48e+6_wp,-1._wp/)

  TYPE t_sse_init_vars
    REAL(wp), POINTER ::     &
      & fao         (:,:  ) => NULL(), &
      & vol_heat_cap(:,:  ) => NULL(), &
      & heat_cond   (:,:  ) => NULL(), &
      & t_srf_mon   (:,:,:) => NULL(), &
      & ifs_tsoil   (:,:,:) => NULL(), &
      & ifs_tskin   (:,:  ) => NULL()
  END TYPE t_sse_init_vars

  TYPE(t_sse_init_vars) :: sse_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_init'

CONTAINS

  !
  !> Intialize soil and snow energy process (after memory has been set up)
  !
  SUBROUTINE sse_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER  :: model
    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init'
    
    model => Get_model(tile%owner_model_id)

    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL sse_read_init_vars(tile)
    END IF

    CALL sse_init_bc(tile)
    CALL sse_init_ic(tile)

    ! IF (tile%Is_last_process_tile(SSE_)) THEN
    !   CALL sse_finalize_init_vars()
    ! END IF

  END SUBROUTINE sse_init

  SUBROUTINE sse_finalize_init_vars

    DEALLOCATE(                     &
      & sse_init_vars%fao,          &
      & sse_init_vars%vol_heat_cap, &
      & sse_init_vars%heat_cond     &
      & )

    IF (ASSOCIATED(sse_init_vars%t_srf_mon)) DEALLOCATE(sse_init_vars%t_srf_mon)

    IF (ASSOCIATED(sse_init_vars%ifs_tsoil)) DEALLOCATE(sse_init_vars%ifs_tsoil)
    IF (ASSOCIATED(sse_init_vars%ifs_tskin)) DEALLOCATE(sse_init_vars%ifs_tskin)

  END SUBROUTINE sse_finalize_init_vars
      
  SUBROUTINE sse_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SSE_)

    REAL(wp), POINTER :: ptr_2D(:,:), ptr_3D(:,:,:) !< temporary pointers

    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e
    INTEGER :: nproma, nblks, isoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()
    soil_e => Get_vgrid('soil_depth_energy')

    ALLOCATE(                                          &
      & sse_init_vars%fao         (nproma, nblks    ), &
      & sse_init_vars%vol_heat_cap(nproma, nblks    ), &
      & sse_init_vars%heat_cond   (nproma, nblks    )  &
      & )

    sse_init_vars%fao         (:,:)   = missval
    sse_init_vars%vol_heat_cap(:,:)   = missval
    sse_init_vars%heat_cond   (:,:)   = missval

    IF (model%config%init_from_ifs) THEN
      ALLOCATE(sse_init_vars%ifs_tsoil(nproma,ifs_nsoil,nblks))
      sse_init_vars%ifs_tsoil(:,:,:) = missval
      ALLOCATE(sse_init_vars%ifs_tskin(nproma,nblks))
      sse_init_vars%ifs_tskin(:,:) = missval
    ELSE
      ALLOCATE(sse_init_vars%t_srf_mon(nproma, nblks, 12))
      sse_init_vars%t_srf_mon(:,:,:) = missval
    END IF


    ! Get soil config
    dsl4jsb_Get_config(SSE_)

    IF (.NOT. tile%Is_process_calculated(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Reading soil energy boundary conditions for tile '//TRIM(tile%name))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SSE_)%bc_filename), model%grid_id)

    ! Thermal diffusivity of soil
    ptr_2D => input_file%Read_2d(                    &
      & variable_name='fao')
    ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)
    sse_init_vars%fao(:,:) = MAX(REAL(LBOUND(fao_thermal_diff,1),wp), MIN(REAL(UBOUND(fao_thermal_diff,1),wp), ptr_2D(:,:)))

    input_file = jsb_netcdf_open_input(dsl4jsb_Config(SSE_)%bc_filename, model%grid_id)

    ! Volumetric heat capacity of dry soil [J/(m^3*K)]
    ptr_2D => input_file%Read_2d(         &
      & variable_name='heat_capacity',            &
      & fill_array = sse_init_vars%vol_heat_cap)
    ptr_2D(:,:) = MERGE(ptr_2D(:,:), 2.e6_wp, ptr_2D(:,:) > 0._wp)

    ! Heat conductivity of mineral soil
    ptr_2D => input_file%Read_2d(         &
      & variable_name='heat_conductivity',         &
      & fill_array = sse_init_vars%heat_cond)
    ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) > 0._wp)

    CALL input_file%Close()

    IF (model%config%init_from_ifs) THEN
      DO isoil=1,ifs_nsoil
        ptr_2D => model%config%ifs_input_file%Read_2d_1lev_1time( &
          & variable_name='STL'//int2string(isoil))
        sse_init_vars%ifs_tsoil(:,isoil,:) = ptr_2D(:,:)
      END DO
      ptr_3D => model%config%ifs_input_file%Read_2d_time( &
        & variable_name='SKT', start_time_step=1, end_time_step=1)
      sse_init_vars%ifs_tskin(:,:) = ptr_3D(:,:,1)

    ELSE

      input_file = jsb_netcdf_open_input(dsl4jsb_Config(SSE_)%ic_filename, model%grid_id)

      ! Read climatological surface temperature
      ptr_3D => input_file%Read_2d_time( &
        & variable_name='surf_temp')
      sse_init_vars%t_srf_mon(:,:,:) = MERGE(ptr_3D, 280._wp, ptr_3D > 0._wp)

      CALL input_file%Close()

    END IF

  END SUBROUTINE sse_read_init_vars

  SUBROUTINE sse_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)

    INTEGER :: iblk

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get soil config
    dsl4jsb_Get_config(SSE_)

    IF (.NOT. tile%Is_process_calculated(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting soil energy boundary conditions for tile '//TRIM(tile%name))

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SSE_)

    ! Thermal diffusivity of soil
    DO iblk=1,SIZE(sse_init_vars%fao,2)
      dsl4jsb_var_ptr(SSE_,thermal_diffusivity)(:,iblk) = fao_thermal_diff(NINT(sse_init_vars%fao(:,iblk)))
    END DO

    ! Volumetric heat capacity of dry soil [J/(m^3*K)]
    IF (dsl4jsb_Config(SSE_)%l_heat_cap_map) THEN
      dsl4jsb_var_ptr(SSE_, vol_heat_cap) = sse_init_vars%vol_heat_cap
    ELSE
      DO iblk=1,SIZE(sse_init_vars%fao,2)
        dsl4jsb_var_ptr(SSE_,vol_heat_cap)(:,iblk) = fao_vol_hcap(NINT(sse_init_vars%fao(:,iblk)))
      END DO
    END IF

    ! Heat conductivity of mineral soil
    IF (dsl4jsb_Config(SSE_)%l_heat_cond_map) THEN
      dsl4jsb_var_ptr(SSE_, heat_cond) = sse_init_vars%heat_cond
    ELSE
      dsl4jsb_var2D_onDomain(SSE_, heat_cond) = &
        & dsl4jsb_var2D_onDomain(SSE_, vol_heat_cap) * dsl4jsb_var2D_onDomain(SSE_, thermal_diffusivity)
    END IF

    !$acc data present(dsl4jsb_var2D_onDomain(SSE_,thermal_diffusivity), &
    !$acc dsl4jsb_var2D_onDomain(SSE_,vol_heat_cap), &
    !$acc dsl4jsb_var2D_onDomain(SSE_,heat_cond))
    !$acc update device(dsl4jsb_var2D_onDomain(SSE_,thermal_diffusivity), &
    !$acc dsl4jsb_var2D_onDomain(SSE_,vol_heat_cap), &
    !$acc dsl4jsb_var2D_onDomain(SSE_,heat_cond))
    !$acc end data

  END SUBROUTINE sse_init_bc

  SUBROUTINE sse_init_ic(tile)

    USE mo_sse_process, ONLY: init_soil_temperature
    USE mo_util,        ONLY: ifs2soil

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_e

    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(SEB_)

    INTEGER :: iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init_ic'

    model => get_model(tile%owner_model_id)
    soil_e => Get_vgrid('soil_depth_energy')

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(SEB_)

    ! Get soil config
    dsl4jsb_Get_config(SSE_)

    IF (.NOT. tile%Is_process_calculated(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing soil energy memory for tile '//TRIM(tile%name))

    ! Soil and surface temperatures
    !
    IF (model%config%init_from_ifs) THEN

      ! CALL ifs2soil(sse_init_vars%ifs_tsoil(:,:,:), ifs_soil_depth(:), dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,:,:), soil_e%ubounds(:), &
      !   & ifs_sfc=sse_init_vars%ifs_tskin(:,:))
      CALL ifs2soil(sse_init_vars%ifs_tsoil(:,:,:), ifs_soil_depth(:), dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,:,:), soil_e%ubounds(:))

    ELSE
      DO iblk=1,SIZE(sse_init_vars%t_srf_mon, DIM=2)
        CALL init_soil_temperature(                   &
          & soil_e%levels(:),                         &
          & sse_init_vars%t_srf_mon(:,iblk,1:12),     &
          & dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,:,iblk))
      END DO
    END IF

    !$acc data present(dsl4jsb_var3D_onDomain(SSE_,t_soil_sl))
    !$acc update device(dsl4jsb_var3D_onDomain(SSE_,t_soil_sl))
    !$acc end data

    ! IF (TRIM(tile%name) /= 'box') THEN
    IF (TRIM(tile%name) == 'land') THEN
      dsl4jsb_var2D_onDomain(SEB_,t)        = dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,1,:)
      dsl4jsb_var2D_onDomain(SEB_,t_unfilt) = dsl4jsb_var2D_onDomain(SEB_,t)

      !$acc data present(dsl4jsb_var2D_onDomain(SEB_,t_unfilt), & 
      !$acc              dsl4jsb_var2D_onDomain(SEB_,t))
      !$acc update device(dsl4jsb_var2D_onDomain(SEB_,t_unfilt), & 
      !$acc               dsl4jsb_var2D_onDomain(SEB_,t))
      !$acc end data
    END IF

    ! Note: called from model_init now
    ! CALL init_soil_properties(tile)

  END SUBROUTINE sse_init_ic

  !-----------------------------------------------------------------------------------------------------
  !> Intialize soil properties - jsbach4
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_soil_properties(tile)

    USE mo_jsb_physical_constants, ONLY: ci, dens_snow_min, tmelt, alf, rhoh2o
    USE mo_sse_constants,          ONLY: vol_hcap_ice, vol_hcap_snow, &   ! vol_hcap_snow_min,
      &                                  vol_hcap_org_top, vol_hcap_org_below, &
      &                                  hcond_ice, hcond_snow, &         !, hcond_snow_min, hcond_org
      &                                  hcond_org_top, hcond_dry_org_top, hcond_org_below, hcond_dry_org_below
    USE mo_hydro_constants,        ONLY: vol_porosity_org_top, vol_porosity_org_below
    USE mo_sse_process,            ONLY: calc_vol_heat_capacity, calc_thermal_conductivity, Get_liquid_max
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! Local variables
    !
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SEB_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: &
      & vol_heat_cap, &
      & heat_cond, &
      & vol_heat_cap_snow, &
      & heat_cond_snow
    dsl4jsb_Real3D_onDomain :: &
      & soil_depth_sl, &
      & wtr_soil_sl, &
      & ice_soil_sl, &
      & t_soil_sl, &
      & matrix_pot_sl, &
      & bclapp_sl, &
      & fract_org_sl, &
      & vol_heat_cap_sl, &
      & heat_cond_sl, &
      & vol_porosity_sl, &
      & vol_field_cap_sl

    ! Locally allocated vectors
    ! 
    dsl4jsb_Real3D_onDomain :: &
      & ws_max_sl,             &
      & liquid_max,            &
      & heat_cap_sl,           &
      & fract_water,           &
      & fract_ice,             &
      & soil_e_dz_tmp,         & 
      & vol_heat_cap_tmp,      &
      & vol_heat_cap_org_tmp,  &
      & heat_cond_tmp,         &
      & vol_porosity_org_tmp,  & 
      & hcond_org_tmp,         & 
      & hcond_dry_org_tmp

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_e

    INTEGER  :: nsoil, nproma, nblks, i, j, b

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_soil_properties'

    model => Get_model(tile%owner_model_id)

    soil_e => Get_vgrid('soil_depth_energy')
    nsoil = soil_e%n_levels

    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(HYDRO_)

    ! Get reference to variables for domain
    !
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)

    dsl4jsb_Get_var3D_onDomain(SSE_, vol_heat_cap_sl)
    dsl4jsb_Get_var3D_onDomain(SSE_, heat_cond_sl)
    dsl4jsb_Get_var2D_onDomain(SSE_, vol_heat_cap_snow)
    dsl4jsb_Get_var2D_onDomain(SSE_, heat_cond_snow)

    ! On glacier tile, set variables and return
    IF (tile%is_glacier) THEN
      !$claw expand parallel update
      vol_heat_cap_sl  (:,:,:) = vol_hcap_ice
      heat_cond_sl     (:,:,:) = hcond_ice
      !$claw end expand
      !$claw expand parallel update
      vol_heat_cap_snow(:  ,:) = vol_hcap_snow
      heat_cond_snow   (:  ,:) = hcond_snow
      !$claw end expand
      !$claw acc update host(vol_heat_cap_sl, heat_cond_sl, vol_heat_cap_snow, heat_cond_snow)
      RETURN
    END IF

    dsl4jsb_Get_var2D_onDomain(SSE_,   vol_heat_cap)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, ice_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)

    nproma = SIZE(vol_heat_cap,1)
    nblks  = SIZE(vol_heat_cap,2)

    IF (model%config%init_from_ifs .AND. dsl4jsb_Config(SSE_)%l_freeze) THEN

      dsl4jsb_Get_var3D_onDomain(SSE_,   t_soil_sl)
      dsl4jsb_Get_var3D_onDomain(HYDRO_, matrix_pot_sl)
      dsl4jsb_Get_var3D_onDomain(HYDRO_, bclapp_sl)
  
      ALLOCATE( &
        & liquid_max (nproma, nsoil, nblks), &
        & ws_max_sl  (nproma, nsoil, nblks), &
        & heat_cap_sl(nproma, nsoil, nblks)  &
        & )
      !$claw acc enter data create(liquid_max, ws_max_sl, heat_cap_sl)

      !$claw expand parallel update
      ws_max_sl(:,:,:) = soil_depth_sl(:,:,:) * vol_porosity_sl(:,:,:)

      IF (dsl4jsb_Config(SSE_)%l_supercool) THEN
        !$claw expand parallel update savepoint=ser_jsb_sse_init_get_liquid_max_1 tile_name=tile%name
        liquid_max(:,:,:) = Get_liquid_max( &
          & t_soil_sl    (:,:,:), &
          & ws_max_sl    (:,:,:), &
          & matrix_pot_sl(:,:,:), &
          & bclapp_sl    (:,:,:)  &
          & )
      ELSE
        !$claw expand parallel update
        liquid_max(:,:,:) = 0._wp
      END IF

      !$claw expand parallel update
      ice_soil_sl(:,:,:) = 0._wp

      DO j = 1, nsoil
        heat_cap_sl(:,j,:) = soil_e%dz(j) * vol_heat_cap(:,:)
      END DO
      !$claw acc update device(heat_cap_sl)

      !$claw acc data present(heat_cap_sl,t_soil_sl,wtr_soil_sl,liquid_max,ice_soil_sl)
      !$claw acc parallel loop gang vector collapse(3)
      DO i = 1, nproma
        DO j = 1, nsoil
          DO b = 1, nblks
            IF (t_soil_sl(i,j,b) < tmelt .AND. wtr_soil_sl(i,j,b) > liquid_max(i,j,b)) THEN
              ice_soil_sl(i,j,b) = MIN(wtr_soil_sl(i,j,b) - liquid_max(i,j,b),                       &
                &                   heat_cap_sl(i,j,b) * (tmelt - t_soil_sl(i,j,b)) / (alf * rhoh2o) &
                &                  )
            END IF
          END DO
        END DO
      END DO
      !$claw acc end data

      !$claw expand parallel update
      wtr_soil_sl(:,:,:) = wtr_soil_sl(:,:,:) - ice_soil_sl(:,:,:)

      !$claw acc update host(wtr_soil_sl,ice_soil_sl)

      ! t_soil_sl(:,:,:) = t_soil_sl(:,:,:) + ice_soil_sl(:,:,:) * alf * rhoh2o / heat_cap_sl(:,:,:)
      ! IF (TRIM(tile%name) == 'land') THEN
      !   dsl4jsb_var_ptr(SEB_,t)(:,:)        = t_soil_sl(:,1,:)
      !   dsl4jsb_var_ptr(SEB_,t_unfilt)(:,:) = t_soil_sl(:,1,:)
      ! END IF

      !$claw acc exit data delete(liquid_max, ws_max_sl, heat_cap_sl)
      DEALLOCATE(liquid_max, ws_max_sl, heat_cap_sl)

      ! dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl) = &
      !   & MIN(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl), dsl4jsb_var3D_onDomain(HYDRO_, w_soil_fc_sl))
      ! !$acc update device(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl))

    END IF

    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_field_cap_sl)
    dsl4jsb_Get_var2D_onDomain(SSE_,   heat_cond)
    nproma = SIZE(heat_cond,1)
    nblks  = SIZE(heat_cond,2)
    
    IF (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil) THEN
      dsl4jsb_Get_var3D_onDomain(HYDRO_, fract_org_sl)
    ELSE
      ALLOCATE(fract_org_sl(nproma, nsoil, nblks))
      !$claw acc enter data create(fract_org_sl)
      !$claw expand parallel update
      fract_org_sl(:,:,:) = 0._wp
    END IF

    ALLOCATE(soil_e_dz_tmp(nproma,nsoil,nblks),vol_heat_cap_tmp(nproma,nsoil,nblks),heat_cond_tmp(nproma,nsoil,nblks))
    !$claw acc enter data pcreate(soil_e_dz_tmp,vol_heat_cap_tmp,heat_cond_tmp)
    soil_e_dz_tmp(:,:,:) = SPREAD(SPREAD(soil_e%dz(:), DIM=1, ncopies=nproma), DIM=3, ncopies=nblks)
    vol_heat_cap_tmp(:,:,:) = SPREAD( vol_heat_cap(:,:), DIM=2, ncopies=nsoil)
    heat_cond_tmp(:,:,:) = SPREAD( heat_cond(:,:),DIM=2, ncopies=nsoil)
    !$claw acc update device(soil_e_dz_tmp,vol_heat_cap_tmp,heat_cond_tmp)

    IF (dsl4jsb_Config(SSE_)%l_heat_cap_dyn) THEN
      ALLOCATE(fract_water(nproma, nsoil, nblks), fract_ice(nproma, nsoil, nblks), vol_heat_cap_org_tmp(nproma, nsoil, nblks))
      !$claw acc enter data pcreate(fract_water,fract_ice,vol_heat_cap_org_tmp)

      !$claw acc data present(soil_depth_sl,fract_water,fract_ice,wtr_soil_sl,ice_soil_sl,vol_heat_cap_org_tmp)
      !$claw acc update device(soil_depth_sl,wtr_soil_sl,ice_soil_sl)
      !$claw acc parallel loop gang vector collapse(3)
      DO i = 1, nproma
        DO j = 1, nsoil
          DO b = 1, nblks
            IF (soil_depth_sl(i,j,b) > 0._wp) THEN
              fract_water(i,j,b) = wtr_soil_sl(i,j,b) / soil_depth_sl(i,j,b)
              fract_ice  (i,j,b) = ice_soil_sl(i,j,b) / soil_depth_sl(i,j,b)
            ELSE
              fract_water(i,j,b) = 0._wp
              fract_ice  (i,j,b) = 0._wp
            END IF
          END DO
        END DO
      END DO
      vol_heat_cap_org_tmp(:,1,:) = vol_hcap_org_top
      DO j = 2, nsoil
        vol_heat_cap_org_tmp(:,j,:) = vol_hcap_org_below
      END DO
      !$claw acc update host(fract_water,fract_ice)
      !$claw acc end data
      
      !$claw expand parallel update savepoint=ser_jsb_sse_init_calc_vol_heat_capacity_1 tile_name=tile%name
      vol_heat_cap_sl (:,:,:) = calc_vol_heat_capacity( &
        !&                         SPREAD( vol_heat_cap     (:,:),   DIM=2, ncopies=nsoil)                        , &
        &                                 vol_heat_cap_tmp (:,:,:)                                               , &
        &                                 vol_heat_cap_org_tmp(:,:,:)                                            , &
        &                                 vol_porosity_sl  (:,:,:)                                               , &
        &                                 fract_org_sl     (:,:,:)                                               , &
        &                                 fract_water      (:,:,:)                                               , &
        &                                 fract_ice        (:,:,:)                                               , &
        &                                 soil_depth_sl    (:,:,:)                                               , &
        &                                 soil_e_dz_tmp    (:,:,:)                                                 &
        !&                   SPREAD(SPREAD(soil_e%dz        (:),     DIM=1, ncopies=nproma), DIM=3, ncopies=nblks)  &
        & )

      !$claw acc exit data delete(fract_water, fract_ice, vol_heat_cap_org_tmp)
      DEALLOCATE(fract_water, fract_ice, vol_heat_cap_org_tmp)
    ELSE

      !$claw acc data present(vol_heat_cap_sl,vol_heat_cap)
      !$claw acc update device(vol_heat_cap)
      !$claw acc parallel loop gang vector collapse(3)
      DO i = 1, nproma
        DO j = 1, nsoil
          DO b=1,nblks
            vol_heat_cap_sl(i,j,b) = vol_heat_cap(i,b)
          END DO
        END DO
      END DO
      !$claw acc end data
    END IF
    !$claw acc update host(vol_heat_cap_sl)

    IF (dsl4jsb_Config(SSE_)%l_heat_cond_dyn) THEN
      ALLOCATE( &
        & vol_porosity_org_tmp(nproma,nsoil,nblks), &
        & hcond_org_tmp(nproma,nsoil,nblks),        &
        & hcond_dry_org_tmp(nproma,nsoil,nblks))

      vol_porosity_org_tmp(:,1,:)   = vol_porosity_org_top
      hcond_org_tmp(:,1,:)          = hcond_org_top
      hcond_dry_org_tmp(:,1,:)      = hcond_dry_org_top
      DO j = 2, nsoil
        vol_porosity_org_tmp(:,j,:) = vol_porosity_org_below
        hcond_org_tmp(:,j,:)        = hcond_org_below
        hcond_dry_org_tmp(:,j,:)    = hcond_dry_org_below
      END DO

      !$claw expand parallel update savepoint=ser_jsb_sse_init_calc_thermal_conductivity_1 tile_name=tile%name
      heat_cond_sl (:,:,:) = calc_thermal_conductivity( &
        !&                              SPREAD( heat_cond(:,:),                    DIM=2, ncopies=nsoil)               , &
        &                                      heat_cond_tmp(:,:,:)                                                   , &
        &                                      wtr_soil_sl(:,:,:)                                                     , &
        &                                      ice_soil_sl(:,:,:)                                                     , &
        &                                      vol_porosity_sl(:,:,:)                                                 , &
        &                                      vol_field_cap_sl(:,:,:)                                                , &
        &                                      fract_org_sl(:,:,:)                                                    , &
        &                                      vol_porosity_org_tmp(:,:,:)                                            , &
        &                                      hcond_org_tmp(:,:,:)                                                   , &
        &                                      hcond_dry_org_tmp(:,:,:)                                               , &
        &                                      soil_depth_sl(:,:,:)                                                   , &
        &                                      soil_e_dz_tmp(:,:,:)                                                     &
        !&                   SPREAD(SPREAD(soil_e%dz        (:),     DIM=1, ncopies=nproma), DIM=3, ncopies=nblks)  &
        & )
      !$claw acc update host(heat_cond_sl)
      DEALLOCATE(vol_porosity_org_tmp,hcond_org_tmp,hcond_dry_org_tmp)
    ELSE
      !$claw acc data present(heat_cond_sl,heat_cond)
      !$claw acc update device(heat_cond)
      !$claw acc parallel loop gang vector collapse(3)
      DO i = 1, nproma
        DO j = 1, nsoil
          DO b=1,nblks
            heat_cond_sl(i,j,b) = heat_cond(i,b)
          END DO
        END DO
      END DO
      !$claw acc update host(heat_cond_sl)
      !$claw acc end data
    END IF

    IF (.NOT. (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil)) THEN
      !$claw acc exit data delete(fract_org_sl)
      DEALLOCATE(fract_org_sl)
    END IF

    IF (dsl4jsb_Config(SSE_)%l_snow) THEN
      IF (model%config%l_compat401) THEN
        !$claw expand parallel update
        vol_heat_cap_snow(:,:) = 634500._wp
        !$claw expand parallel update
        heat_cond_snow   (:,:) = 0.31_wp
        !$cla$w end expand
      ELSE
        IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
          ! USe snow density values for fresh snow for initialization
          !$claw expand parallel update
          vol_heat_cap_snow(:,:) = ci * dens_snow_min                               ! Use specific heat of ice (Verseghy 1991)
          !$claw expand parallel update
          heat_cond_snow   (:,:) = 2.5e-6_wp*dens_snow_min**2._wp - 0.123e-3_wp*dens_snow_min + 0.024_wp ! Calonne et al. (2011)
          !$cla$w end expand
          ! vol_heat_cap_snow(:,:) = vol_hcap_snow_min
          ! heat_cond_snow   (:,:) = hcond_snow_min
        ELSE
          !$claw expand parallel update
          vol_heat_cap_snow(:,:) = vol_hcap_snow
          !$claw expand parallel update
          heat_cond_snow   (:,:) = hcond_snow
          !$cla$w end expand
        END IF
      END IF
      !$claw acc update host(vol_heat_cap_snow, heat_cond_snow)
    ! ELSE
    !   vol_heat_cap_snow and heat_cond_snow are fixed in calc_soiltemp_old
    END IF

    !$acc exit data delete(soil_e_dz_tmp,vol_heat_cap_tmp,heat_cond_tmp)
    DEALLOCATE(soil_e_dz_tmp,vol_heat_cap_tmp,heat_cond_tmp)

  END SUBROUTINE init_soil_properties

#endif
END MODULE mo_sse_init
