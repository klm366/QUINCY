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

  dsl4jsb_Use_processes SSE_, SEB_, HYDRO_
  dsl4jsb_Use_config(SSE_)
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

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init'

    CALL sse_init_bc(tile)
    CALL sse_init_ic(tile)

  END SUBROUTINE sse_init

  SUBROUTINE sse_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)

    REAL(wp), POINTER :: &
      & fao    (:,  :), & !< FAO soil types
      & ptr_2D (:,  :)    !< temporary pointer

    INTEGER :: iblk

    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get soil config
    dsl4jsb_Get_config(SSE_)

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting  soil energy boundary conditions for tile '//TRIM(tile%name))

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SSE_)

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SSE_)%bc_filename), model%grid_id)

    ! Thermal diffusivity of soil
    fao => input_file%Read_2d(                    &
      & variable_name='fao')
    fao(:,:) = MERGE(fao(:,:), 0._wp, fao(:,:) >= 0._wp)
    fao(:,:) = MAX(REAL(LBOUND(fao_thermal_diff,1),wp), MIN(REAL(UBOUND(fao_thermal_diff,1),wp), fao(:,:)))

    DO iblk=1,SIZE(fao,2)
      dsl4jsb_var_ptr(SSE_,thermal_diffusivity)(:,iblk) = fao_thermal_diff(NINT(fao(:,iblk)))
    END DO

    ! Volumetric heat capacity of dry soil [J/(m^3*K)]
    IF (dsl4jsb_Config(SSE_)%l_heat_cap_map) THEN
      ptr_2D => input_file%Read_2d(         &
        & variable_name='heat_capacity',            &
        & fill_array = dsl4jsb_var_ptr(SSE_, vol_heat_cap))
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 2.e6_wp, ptr_2D(:,:) > 0._wp)
    ELSE
      DO iblk=1,SIZE(fao,2)
        dsl4jsb_var_ptr(SSE_,vol_heat_cap)(:,iblk) = fao_vol_hcap(NINT(fao(:,iblk)))
      END DO
    END IF

    ! Heat conductivity of mineral soil
    IF (dsl4jsb_Config(SSE_)%l_heat_cond_map) THEN
      ptr_2D => input_file%Read_2d(         &
        & variable_name='heat_conductivity',         &
        & fill_array = dsl4jsb_var_ptr(SSE_, heat_cond))
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) > 0._wp)
    ELSE
      dsl4jsb_var2D_onDomain(SSE_, heat_cond) = &
        & dsl4jsb_var2D_onDomain(SSE_, vol_heat_cap) * dsl4jsb_var2D_onDomain(SSE_, thermal_diffusivity)
    END IF

    CALL input_file%Close()                              ! Temporary

  END SUBROUTINE sse_init_bc

  SUBROUTINE sse_init_ic(tile)

    USE mo_sse_process, ONLY: init_soil_temperature

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(SEB_)

    REAL(wp), POINTER :: &
      & t_srf_mon (:,:,:) !,    & !< Climatological monthly values of surface temperature
      !& ptr_2D    (:,  :),    & !< temporary pointer
      !& ptr_3D    (:,:,:)

    TYPE(t_input_file) :: input_file
    TYPE(t_jsb_vgrid), POINTER :: soil_e
    INTEGER :: iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init_ic'

    model => get_model(tile%owner_model_id)

    soil_e => Get_vgrid('soil_depth_energy')

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(SEB_)

    ! Get soil config
    dsl4jsb_Get_config(SSE_)

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing soil energy memory for tile '//TRIM(tile%name))

    input_file = jsb_netcdf_open_input(dsl4jsb_Config(SSE_)%ic_filename, model%grid_id)

    ! Soil and surface temperatures
    !
    ! Read climatological surface temperature
    t_srf_mon => input_file%Read_2d_time(         &   ! Temporary
      & variable_name='surf_temp')
    t_srf_mon(:,:,:) = MERGE(t_srf_mon(:,:,:), 280._wp, t_srf_mon(:,:,:) > 0._wp)

    DO iblk=1,SIZE(t_srf_mon, DIM=2)
      CALL init_soil_temperature(                     &
        & soil_e%levels(:),                           &
        & t_srf_mon(:,iblk,1:12),                     &
        & dsl4jsb_var_ptr(SSE_,t_soil)(:,:,iblk),   &
        & dsl4jsb_var_ptr(SEB_,t)(:,iblk))
    END DO
    dsl4jsb_var2D_onDomain(SEB_,t_unfilt) = dsl4jsb_var2D_onDomain(SEB_,t)

    !$acc update device(dsl4jsb_var2D_onDomain(SEB_,t_unfilt), & 
    !$acc dsl4jsb_var2D_onDomain(SEB_,t), &
    !$acc dsl4jsb_var_ptr(SSE_,t_soil))

    CALL input_file%Close()

    CALL init_soil_properties(tile)

  END SUBROUTINE sse_init_ic

  SUBROUTINE init_soil_properties(tile)

    USE mo_jsb_physical_constants, ONLY: ci, dens_snow_min
    USE mo_sse_constants,          ONLY: vol_hcap_ice, vol_hcap_snow, & ! vol_hcap_snow_min,
      &                                  hcond_ice, hcond_snow !, hcond_snow_min, hcond_org
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
      & heat_cond_sl, &
      & w_soil_sl, &
      & w_ice_sl, &
      & vol_porosity_sl, &
      & vol_field_cap_sl

    ! Locally allocated vectors
    ! 
    dsl4jsb_Real3D_onDomain :: &
      & fract_water,           &
      & fract_ice,             &
      & soil_e_dz_tmp,         & 
      & vol_heat_cap_tmp,      &
      & heat_cond_tmp

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_e

    INTEGER  :: nsoil, nproma, nblks, i, j, b, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_soil_properties'

    model => Get_model(tile%owner_model_id)

    soil_e => Get_vgrid('soil_depth_energy')
    nsoil = soil_e%n_levels

    nc = nproma

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

    dsl4jsb_Get_var3D_onDomain(HYDRO_, w_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, w_ice_sl )
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_field_cap_sl)

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
      ALLOCATE(fract_water(nproma, nsoil, nblks), fract_ice(nproma, nsoil, nblks))
      !$claw acc enter data pcreate(fract_water,fract_ice)

      !$claw acc data present(soil_depth_sl,fract_water,fract_ice,w_soil_sl,w_ice_sl)
      !$claw acc update device(soil_depth_sl,w_soil_sl,w_ice_sl)
      !$claw acc parallel loop gang vector collapse(3)
      DO i = 1, nproma
        DO j = 1, nsoil
          DO b = 1, nblks
            IF (soil_depth_sl(i,j,b) > 0._wp) THEN
              fract_water(i,j,b) = w_soil_sl(i,j,b) / soil_depth_sl(i,j,b)
              fract_ice  (i,j,b) = w_ice_sl(i,j,b) / soil_depth_sl(i,j,b)
            ELSE
              fract_water(i,j,b) = 0._wp
              fract_ice  (i,j,b) = 0._wp
            END IF
          END DO
        END DO
      END DO
      !$claw acc update host(fract_water,fract_ice)
      !$claw acc end data
      
      !$claw expand parallel update savepoint=ser_jsb_sse_init_calc_vol_heat_capacity_1 tile_name=tile%name
      vol_heat_cap_sl (:,:,:) = calc_vol_heat_capacity( &
        !&                         SPREAD( vol_heat_cap     (:,:),   DIM=2, ncopies=nsoil)                        , &
        &                                 vol_heat_cap_tmp (:,:,:)                                               , &
        &                                 vol_porosity_sl  (:,:,:)                                               , &
        &                                 fract_org_sl     (:,:,:)                                               , &
        &                                 fract_water      (:,:,:)                                               , &
        &                                 fract_ice        (:,:,:)                                               , &
        &                                 soil_depth_sl    (:,:,:)                                               , &
        &                                 soil_e_dz_tmp    (:,:,:)                                                 &
        !&                   SPREAD(SPREAD(soil_e%dz        (:),     DIM=1, ncopies=nproma), DIM=3, ncopies=nblks)  &
        & )

      !$claw acc exit data delete(fract_water, fract_ice)
      DEALLOCATE(fract_water, fract_ice)
    ELSE

      
      vol_heat_cap_sl(:,:,:) = SPREAD( vol_heat_cap(:,:), DIM=2, ncopies=nsoil)

    END IF

    IF (dsl4jsb_Config(SSE_)%l_heat_cond_dyn) THEN
      !$claw expand parallel update savepoint=ser_jsb_sse_init_calc_thermal_conductivity_1 tile_name=tile%name
      heat_cond_sl (:,:,:) = calc_thermal_conductivity( &
        !&                              SPREAD( heat_cond(:,:),                    DIM=2, ncopies=nsoil)               , &
        &                                      heat_cond_tmp(:,:,:)                                                   , &
        &                                      w_soil_sl(:,:,:)                                                       , &
        &                                      w_ice_sl(:,:,:)                                                        , &
        &                                      vol_porosity_sl(:,:,:)                                                 , &
        &                                      vol_field_cap_sl(:,:,:)                                                , &
        &                                      fract_org_sl(:,:,:)                                                    , &
        &                                      soil_depth_sl(:,:,:)                                                   , &
        &                                      soil_e_dz_tmp(:,:,:)                                                     &
        !&                   SPREAD(SPREAD(soil_e%dz        (:),     DIM=1, ncopies=nproma), DIM=3, ncopies=nblks)  &
        & )
    ELSE
      heat_cond_sl(:,:,:) = SPREAD(heat_cond(:,:), DIM=2, ncopies=nsoil)
    END IF

    IF (.NOT. (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil)) THEN
      !$claw acc exit data delete(fract_org_sl)
      DEALLOCATE(fract_org_sl)
    END IF

    IF (dsl4jsb_Config(SSE_)%l_snow) THEN
      IF (model%config%l_compat401) THEN
        !$claw expand parallel update
        vol_heat_cap_snow(:,:) = 634500._wp
        heat_cond_snow   (:,:) = 0.31_wp
        !$claw end expand
      ELSE
        IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
          ! USe snow density values for fresh snow for initialization
          !$claw expand parallel update
          vol_heat_cap_snow(:,:) = ci * dens_snow_min                               ! Use specific heat of ice (Verseghy 1991)
          heat_cond_snow   (:,:) = 2.5e-6_wp*dens_snow_min**2._wp - 0.123e-3_wp*dens_snow_min + 0.024_wp ! Calonne et al. (2011)
          !$claw end expand
          ! vol_heat_cap_snow(:,:) = vol_hcap_snow_min
          ! heat_cond_snow   (:,:) = hcond_snow_min
        ELSE
          !$claw expand parallel update
          vol_heat_cap_snow(:,:) = vol_hcap_snow
          heat_cond_snow   (:,:) = hcond_snow
          !$claw end expand
        END IF
      END IF
    ! ELSE
    !   vol_heat_cap_snow and heat_cond_snow are fixed in calc_soiltemp_old
    END IF

    !$acc exit data delete(soil_e_dz_tmp,vol_heat_cap_tmp,heat_cond_tmp)
    DEALLOCATE(soil_e_dz_tmp,vol_heat_cap_tmp,heat_cond_tmp)

  END SUBROUTINE init_soil_properties

#endif
END MODULE mo_sse_init
