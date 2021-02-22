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

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message_text, finish, message
  USE mo_jsb_control,       ONLY: debug_on

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,            ONLY: missval

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

  TYPE t_sse_init_vars
    REAL(wp), POINTER ::     &
      & fao         (:,:  ), &
      & vol_heat_cap(:,:  ), &
      & heat_cond   (:,:  ), &
      & t_srf_mon   (:,:,:)
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
    !! quincy
    IF (model%config%use_quincy) THEN
      CALL sse_init_q_(tile)
    ENDIF

    ! TODO This is not correct if sse runs on veg tile!
    IF (.NOT. tile%Has_children() .AND. .NOT. ASSOCIATED(tile%next_sibling_tile)) THEN
      CALL sse_finalize_init_vars()
    END IF

  END SUBROUTINE sse_init

  SUBROUTINE sse_finalize_init_vars

    DEALLOCATE(                     &
      & sse_init_vars%fao,          &
      & sse_init_vars%vol_heat_cap, &
      & sse_init_vars%heat_cond,    &
      & sse_init_vars%t_srf_mon     &
      & )

  END SUBROUTINE sse_finalize_init_vars
      
  SUBROUTINE sse_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SSE_)

    REAL(wp), POINTER :: ptr_2D(:,:), ptr_3D(:,:,:) !< temporary pointers

    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e
    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()
    soil_e => Get_vgrid('soil_depth_energy')

    ALLOCATE(                                          &
      & sse_init_vars%fao         (nproma, nblks    ), &
      & sse_init_vars%vol_heat_cap(nproma, nblks    ), &
      & sse_init_vars%heat_cond   (nproma, nblks    ), &
      & sse_init_vars%t_srf_mon   (nproma, nblks, 12)  &
      & )

    sse_init_vars%fao         (:,:)   = missval
    sse_init_vars%vol_heat_cap(:,:)   = missval
    sse_init_vars%heat_cond   (:,:)   = missval
    sse_init_vars%t_srf_mon   (:,:,:) = missval

    ! Get soil config
    dsl4jsb_Get_config(SSE_)

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

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

    input_file = jsb_netcdf_open_input(dsl4jsb_Config(SSE_)%ic_filename, model%grid_id)

    ! Read climatological surface temperature
    ptr_3D => input_file%Read_2d_time( &
      & variable_name='surf_temp')
    sse_init_vars%t_srf_mon(:,:,:) = MERGE(ptr_3D, 280._wp, ptr_3D > 0._wp)

    CALL input_file%Close()

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

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

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

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing soil energy memory for tile '//TRIM(tile%name))

    ! Soil and surface temperatures
    !
    DO iblk=1,SIZE(sse_init_vars%t_srf_mon, DIM=2)
      CALL init_soil_temperature(                   &
        & soil_e%levels(:),                         &
        & sse_init_vars%t_srf_mon(:,iblk,1:12),     &
        & dsl4jsb_var_ptr(SSE_,t_soil)(:,:,iblk),   &
        & dsl4jsb_var_ptr(SEB_,t)(:,iblk))
    END DO
    dsl4jsb_var2D_onDomain(SEB_,t_unfilt) = dsl4jsb_var2D_onDomain(SEB_,t)

    !$acc data present(dsl4jsb_var2D_onDomain(SEB_,t_unfilt), & 
    !$acc dsl4jsb_var2D_onDomain(SEB_,t), &
    !$acc dsl4jsb_var3D_onDomain(SSE_,t_soil))
    !$acc update device(dsl4jsb_var2D_onDomain(SEB_,t_unfilt), & 
    !$acc dsl4jsb_var2D_onDomain(SEB_,t), &
    !$acc dsl4jsb_var3D_onDomain(SSE_,t_soil))
    !$acc end data

    CALL init_soil_properties(tile)

  END SUBROUTINE sse_init_ic

  !-----------------------------------------------------------------------------------------------------
  !> Intialize sse process for quincy SSE routines
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE sse_init_q_(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_util,                  ONLY: soil_depth_to_layers_2d
    USE mo_jsb_math_constants,    ONLY: eps8
    USE mo_sse_constants,         ONLY: soil_heat_cap, soil_therm_cond, kdiff_sat_max, &
                                        k_pwp_s, k_pwp_c, k_pwp_sc, k_pwp_a, k_pwp_at, k_pwp_bt, &
                                        k_fc_s, k_fc_c, k_fc_sc, k_fc_a, k_fc_at, k_fc_bt, k_fc_ct, &
                                        k_sat_s, k_sat_c, k_sat_sc, k_sat_a, k_sat_at, k_sat_bt, k_sat_ct, k_sat_dt, &
                                        wpot_pwp, wpot_fc, wpot_fc

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SSE_)
    dsl4jsb_Use_config(SB_)

    ! USE PROC_memory
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(VEG_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER                  :: model
    TYPE(t_lctlib_element), POINTER                  :: lctlib        !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER                  :: vgrid_soil_sb !< Vertical grid
    INTEGER                                          :: nsoil_sb      !< number of soil layers as used/defined by the SB_ process
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)            :: hlp1, &
                                                        theta_pwp_sl, theta_fc_sl, theta_sat_sl
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)            :: sand_sse_init_sl, silt_sse_init_sl, clay_sse_init_sl
    REAL(wp), DIMENSION(options%nc)                  :: soil_awc, &
                                                        hlp2
    REAL(wp)                                         :: site_longitude
    LOGICAL                                          :: flag_spp1685
    CHARACTER(len=3)                                 :: site_ID_spp1685
    INTEGER                                          :: isoil, ichunk
    INTEGER                                          :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sse_init_q_'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(VEG_)
    ! Declare pointers to variables in memory
    ! VEG_
    dsl4jsb_Real3D_onChunk    :: root_fraction_sl
    ! SSE_
    dsl4jsb_Real2D_onChunk    :: soil_depth
    dsl4jsb_Real2D_onChunk    :: root_depth
    dsl4jsb_Real2D_onChunk    :: t_srf_old 
    dsl4jsb_Real2D_onChunk    :: w_skin
    dsl4jsb_Real2D_onChunk    :: w_soil_root
    dsl4jsb_Real2D_onChunk    :: w_soil_root_pwp
    dsl4jsb_Real2D_onChunk    :: w_soil_root_fc
    dsl4jsb_Real2D_onChunk    :: w_soil_root_theta
    dsl4jsb_Real2D_onChunk    :: w_soil_root_pot
    dsl4jsb_Real3D_onChunk    :: soil_depth_sl          ! (soil_depth_sl == soil_lay_width_sl)
    dsl4jsb_Real3D_onChunk    :: soil_lay_width_sl
    dsl4jsb_Real3D_onChunk    :: soil_lay_depth_center_sl
    dsl4jsb_Real3D_onChunk    :: soil_lay_depth_ubound_sl
    dsl4jsb_Real3D_onChunk    :: soil_lay_depth_lbound_sl
    dsl4jsb_Real3D_onChunk    :: bulk_dens_sl
    dsl4jsb_Real3D_onChunk    :: heat_capa_sl
    dsl4jsb_Real3D_onChunk    :: therm_cond_sl
    dsl4jsb_Real3D_onChunk    :: w_soil_pwp_sl
    dsl4jsb_Real3D_onChunk    :: w_soil_fc_sl 
    dsl4jsb_Real3D_onChunk    :: w_soil_sat_sl 
    dsl4jsb_Real3D_onChunk    :: w_soil_sl
    dsl4jsb_Real3D_onChunk    :: w_soil_pot_sl
    dsl4jsb_Real3D_onChunk    :: saxtonA 
    dsl4jsb_Real3D_onChunk    :: saxtonB 
    dsl4jsb_Real3D_onChunk    :: saxtonC 
    dsl4jsb_Real3D_onChunk    :: kdiff_sat_sl 
    dsl4jsb_Real3D_onChunk    :: sand_sl
    dsl4jsb_Real3D_onChunk    :: silt_sl
    dsl4jsb_Real3D_onChunk    :: clay_sl
    dsl4jsb_Real3D_onChunk    :: volume_min_sl
    ! Get local variables from options argument
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(SSE_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model          => Get_model(tile%owner_model_id)
    lctlib         => model%lctlib(tile%lcts%lib_id)
    vgrid_soil_sb  => Get_vgrid('soil_layer_sb')
    nsoil_sb       =  vgrid_soil_sb%n_levels
    site_longitude = model%config%lon
    flag_spp1685   = model%config%flag_spp1685
    ! Get process config
    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(SB_)
    ! Get process memories
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(VEG_)
    ! Allocate local allocatable variables (deallocate below at the end of this routine)    
    ALLOCATE(hlp1(nc, nsoil_sb))
    ALLOCATE(theta_pwp_sl(nc, nsoil_sb))
    ALLOCATE(theta_fc_sl(nc, nsoil_sb))
    ALLOCATE(theta_sat_sl(nc, nsoil_sb))
    ALLOCATE(sand_sse_init_sl(nc, nsoil_sb))
    ALLOCATE(clay_sse_init_sl(nc, nsoil_sb))
    ALLOCATE(silt_sse_init_sl(nc, nsoil_sb))
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var3D_onChunk(VEG_,  root_fraction_sl)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SSE_, soil_depth)
    dsl4jsb_Get_var2D_onChunk(SSE_, root_depth)
    dsl4jsb_Get_var2D_onChunk(SSE_, t_srf_old)
    dsl4jsb_Get_var2D_onChunk(SSE_, w_skin)
    dsl4jsb_Get_var2D_onChunk(SSE_, w_soil_root)
    dsl4jsb_Get_var2D_onChunk(SSE_, w_soil_root_pwp)
    dsl4jsb_Get_var2D_onChunk(SSE_, w_soil_root_fc)
    dsl4jsb_Get_var2D_onChunk(SSE_, w_soil_root_theta)
    dsl4jsb_Get_var2D_onChunk(SSE_, w_soil_root_pot)
    dsl4jsb_Get_var3D_onChunk(SSE_, soil_depth_sl)              ! out   (soil_depth_sl == soil_lay_width_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, soil_lay_width_sl)          ! out
    dsl4jsb_Get_var3D_onChunk(SSE_, soil_lay_depth_center_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(SSE_, soil_lay_depth_ubound_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(SSE_, soil_lay_depth_lbound_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(SSE_, bulk_dens_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, heat_capa_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, therm_cond_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, w_soil_pwp_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, w_soil_fc_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, w_soil_sat_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, w_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, w_soil_pot_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, saxtonA)
    dsl4jsb_Get_var3D_onChunk(SSE_, saxtonB)
    dsl4jsb_Get_var3D_onChunk(SSE_, saxtonC)
    dsl4jsb_Get_var3D_onChunk(SSE_, kdiff_sat_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, sand_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, silt_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, clay_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, volume_min_sl)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 1.0 Total soil depth until bedrock
    ! QUINCY soil_depth for one site (read from namelist)
    soil_depth(:)                             = dsl4jsb_Config(SSE_)%soil_depth
    ! code from jsbach4 (hydro init)
    ! for all points of the domain soil_depth(:,:) is at least 0.1
    dsl4jsb_var2D_onDomain(SSE_, soil_depth)  = MAX(0.1_wp, &
                                                    MERGE(dsl4jsb_var2D_onDomain(SSE_, soil_depth), 0._wp, &
                                                          dsl4jsb_var2D_onDomain(SSE_, soil_depth) >= 0._wp))

    !> 2.0 calc actual soil layer depths based on the fixed layer thicknesses from namelist and the bedrock depth
    !!  note, soil_depth_to_layers_2d() returns an actual 3D variable
    dsl4jsb_var3D_onDomain(SSE_, soil_depth_sl) = soil_depth_to_layers_2d( &
                                                    dsl4jsb_var2D_onDomain(SSE_, soil_depth), &   ! Total soil depth until bedrock (from textures), 2D variable
                                                    vgrid_soil_sb%dz(:))                          ! Soil layer thicknesses from namelist
    ! pass these values to a SSE_ variable with a correct name
    dsl4jsb_var3D_onDomain(SSE_, soil_lay_width_sl) = dsl4jsb_var3D_onDomain(SSE_, soil_depth_sl)

    !> 2.1 calc more metrics of the soil layers
    !! i)  lower & upper boundary
    !! ii) depth at the center of each layer
    ! @TODO sse_init must init variables onDomain for jsbach4, write a function similar to soil_depth_to_layers_2d()
    DO isoil = 1,nsoil_sb
      ! lower & upper boundary
      IF (isoil == 1) THEN
        soil_lay_depth_lbound_sl(:,isoil) = 0.0_wp
        soil_lay_depth_ubound_sl(:,isoil) = soil_lay_width_sl(:,isoil)
      ELSE
        soil_lay_depth_lbound_sl(:,isoil) = soil_lay_depth_lbound_sl(:,isoil-1) + soil_lay_width_sl(:,isoil-1)
        soil_lay_depth_ubound_sl(:,isoil) = soil_lay_depth_ubound_sl(:,isoil-1) + soil_lay_width_sl(:,isoil)
      ENDIF                                
      ! soil-layer center
      soil_lay_depth_center_sl(:,isoil)   = (soil_lay_depth_lbound_sl(:,isoil) + soil_lay_depth_ubound_sl(:,isoil)) * 0.5_wp 
    ENDDO

    !> 3.0 check if last soil layer has a thickness larger zero
    !! QUINCY specific code: for jsb4 it would need to work on the domain
    DO ichunk  = 1,nc
      IF (soil_depth_sl(ichunk,nsoil_sb) < eps8) THEN
        WRITE(message_text,*) 'The last soil layer has a thickness < eps8. Increase soil_depth from ', soil_depth, &
                              ' to at least ', SUM(soil_lay_width_sl(ics,1:nsoil_sb-1)), ' but better ', &
                              SUM(soil_lay_width_sl(ics,1:nsoil_sb)), ' and re-start the simulation.'
        CALL finish(TRIM(routine), message_text)
      END IF  
    END DO
    
    !> 4.0 init soil-property variables
    !! 
    !> 4.1 init with one site-specific value for all soil layers
    !! 
    bulk_dens_sl(:,:)   = dsl4jsb_Config(SSE_)%bulk_density 
    heat_capa_sl(:,:)   = soil_heat_cap
    therm_cond_sl(:,:)  = soil_therm_cond
    sand_sl(:,:)        = dsl4jsb_Config(SSE_)%soil_sand
    silt_sl(:,:)        = dsl4jsb_Config(SSE_)%soil_silt
    clay_sl(:,:)        = dsl4jsb_Config(SSE_)%soil_clay
    t_srf_old(:)        = 273.15_wp
    volume_min_sl(:,:)  = 1._wp

    !> 4.2 ensure a certain minimum proportion for sand, clay and silt
    !! this is to ensure applying the empirical functions by Saxton within the range of the texture they were designed for \n
    !! avoiding ranges where calculation would return negative theta_pwp_sl or relatively large kdiff_sat_sl values
    !! 
    !! Saxton equations with too small soil-texture fractions may lead to unrealistic percolation rates and soil water potentials
    !!
    !! for more details see ticket #560 in the quincy svn project system
    !!  
    !! using the modified *_sse_init_sl soil texture for "saxton functions here in sse_init" only
    sand_sse_init_sl(:,:) = sand_sl(:,:)
    clay_sse_init_sl(:,:) = clay_sl(:,:)
    silt_sse_init_sl(:,:) = silt_sl(:,:)
    CALL calc_soil_texture_saxton_compatible(nc, &
                                             nsoil_sb, &
                                             sand_sse_init_sl(:,:), &
                                             silt_sse_init_sl(:,:), &
                                             clay_sse_init_sl(:,:))
    
    !> 5.0 soil hydraulic properties, based on Saxton and Rawls 2006, SSSA
    !! 

    !> 5.1 fractional water holding capacity at permanent wilting point
    !! note that this assumes that the soil hydraulic properties are identical between the layers
    !!
    hlp1(:,:)         = k_pwp_s * sand_sse_init_sl(:,:) + & 
                        k_pwp_c * clay_sse_init_sl(:,:) + &
                        k_pwp_sc * sand_sse_init_sl(:,:) * clay_sse_init_sl(:,:) + &
                        k_pwp_a
    theta_pwp_sl(:,:) = k_pwp_at + (1.0_wp + k_pwp_bt) * hlp1(:,:) 

    !> 5.2 fractional water holding capacity at field capacity 
    !! 
    hlp1(:,:)         = k_fc_s * sand_sse_init_sl(:,:) + &
                        k_fc_c * clay_sse_init_sl(:,:) + &
                        k_fc_sc * sand_sse_init_sl(:,:) * clay_sse_init_sl(:,:) + &
                        k_fc_a
    theta_fc_sl(:,:)  = k_fc_at + (1.0_wp + k_fc_bt) * hlp1(:,:) + k_fc_ct * hlp1(:,:)**2._wp

    !> 5.3 fractional water holding capacity at saturation 
    !! 
    hlp1(:,:)         = k_sat_s * sand_sse_init_sl(:,:) + & 
                        k_sat_c * clay_sse_init_sl(:,:) + &
                        k_sat_sc * sand_sse_init_sl(:,:) * clay_sse_init_sl(:,:) + &
                        k_sat_a
    theta_sat_sl(:,:) = theta_fc_sl(:,:) + k_sat_at + (1.0_wp + k_sat_bt) * hlp1(:,:) - & 
                        k_sat_ct * sand_sse_init_sl(:,:) + k_sat_dt 
    
    !> 5.4 parameters describing the moisture-tension curve                 
    !! note that for convenience, B has a negative sign and A is in MPa 
    saxtonB(:,:)      = - (log(wpot_pwp) - log(wpot_fc))/(log(theta_fc_sl(:,:)) - log(theta_pwp_sl(:,:)))
    saxtonA(:,:)      = - (exp(log(wpot_fc) - saxtonB(:,:) * log(theta_fc_sl(:,:))))/1000.0_wp

    !> 5.5 parameters describing the moisture-diffusivity curve
    !! 
    kdiff_sat_sl(:,:) = kdiff_sat_max * (theta_sat_sl(:,:) - theta_fc_sl(:,:)) ** (3.0_wp + 1._wp/saxtonB(:,:))
    saxtonC(:,:)      = 3.0_wp - 2._wp * saxtonB(:,:)


    !> 6.0 rooting depth given by prescribed AWC
    !! soil water in the rooting zone [m] \n
    !! the init of 'root_fraction_sl(:,:)' (sec 6.3) depends on this init of the root_depth(:)
    !! 
    soil_awc(:)         = 0._wp
    w_soil_root_pwp(:)  = 0._wp
    w_soil_root_fc(:)   = 0._wp
    root_depth(:)       = 0._wp
    DO isoil=1,nsoil_sb
      hlp1(:,isoil)     = soil_lay_width_sl(:,isoil) * (theta_fc_sl(:,isoil) - theta_pwp_sl(:,isoil))
      soil_awc(:)       = soil_awc(:) + hlp1(:,isoil)
      WHERE (soil_awc(:) < (dsl4jsb_Config(SSE_)%soil_awc_prescribe / 1000._wp)) 
        root_depth(:)       = root_depth(:) + soil_lay_width_sl(:,isoil)
        w_soil_root_pwp(:)  = w_soil_root_pwp(:) + theta_pwp_sl(:,isoil) * soil_lay_width_sl(:,isoil)
        w_soil_root_fc(:)   = w_soil_root_fc(:) + theta_fc_sl(:,isoil) * soil_lay_width_sl(:,isoil)
      ELSEWHERE
        hlp2(:)             = MAX(0._wp, dsl4jsb_Config(SSE_)%soil_awc_prescribe / 1000._wp - soil_awc(:) + hlp1(:,isoil))
        root_depth(:)       = root_depth(:) + hlp2(:) / (theta_fc_sl(:,isoil) - theta_pwp_sl(:,isoil))
        w_soil_root_pwp(:)  = w_soil_root_pwp(:) + theta_pwp_sl(:,isoil) * hlp2(:) / (theta_fc_sl(:,isoil) - theta_pwp_sl(:,isoil))
        w_soil_root_fc(:)   = w_soil_root_fc(:) + theta_fc_sl(:,isoil) * hlp2(:) / (theta_fc_sl(:,isoil) - theta_pwp_sl(:,isoil))
      END WHERE
    END DO

    !> 6.1 water in the surface reservoir [m]
    !! 
    w_skin(:) = 0.0_wp

    !> 6.2 adjust initial water amount to prescribed theta 
    !! 
    w_soil_root(:)        = w_soil_root_pwp(:) + & 
                            dsl4jsb_Config(SSE_)%soil_theta_prescribe * (w_soil_root_fc(:) - w_soil_root_pwp(:)) 
    w_soil_root_theta(:)  = dsl4jsb_Config(SSE_)%soil_theta_prescribe
    w_soil_root_pot(:)    = saxtonA(:,1) * (MIN(w_soil_root_fc,w_soil_root(:)) / root_depth(:)) ** saxtonB(:,1) 

    !> 6.3 root fraction 
    !! 
    !! root fractions assuming an exponential distribution of roots \n
    !!   Cr_l = Cr_0 * exp(-lctlib%k_root_dist*depth) \n
    !! integral to bottom of root zone: 1/lctlib%k_root_dist * (1 - exp(-lctlib%k_root_dist*depth) \n
    !!   thus for each layer
    !!
    !! root_fraction_sl(:,:) is a variable of veg_mem and depends on the init of root_depth(:), see section 6.0 and ticket #358
    !! 
    !! note: upper bound is the larger value compared to lower bound !
    !! 
    ! layer 1 (vgrid_soil_sb)
    root_fraction_sl(:,1) = (1._wp - exp(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:,1))) / &
                            (1._wp - exp(-lctlib%k_root_dist * root_depth(:)))
    ! layer 2 (vgrid_soil_sb)
    DO isoil=2,nsoil_sb
       WHERE(root_depth(:) > soil_lay_depth_ubound_sl(:,isoil-1))
          WHERE(root_depth(:) > soil_lay_depth_ubound_sl(:,isoil))
             root_fraction_sl(:,isoil) = (exp(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:,isoil-1)) - &
                                          exp(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:,isoil))) / & 
                                         (1._wp - exp(-lctlib%k_root_dist * root_depth(:)))
          ELSEWHERE
             root_fraction_sl(:,isoil) = (exp(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:,isoil-1)) - &
                                          exp(-lctlib%k_root_dist * root_depth(:))) / & 
                                         (1._wp - exp(-lctlib%k_root_dist * root_depth(:)))
          ENDWHERE
       ELSEWHERE
          root_fraction_sl(:,isoil) = 0.0_wp 
       ENDWHERE
    ENDDO


    !> 7.0 soil water per layer depth
    !! 
    w_soil_pwp_sl(:,:)  = 0.0_wp
    w_soil_fc_sl(:,:)   = 0.0_wp
    w_soil_sat_sl(:,:)  = 0.0_wp
    DO isoil=1,nsoil_sb
       WHERE(soil_depth_sl(:,isoil) > 0.0_wp) 
          w_soil_pwp_sl(:,isoil) = theta_pwp_sl(:,isoil) * soil_depth_sl(:,isoil) 
          w_soil_fc_sl(:,isoil)  = theta_fc_sl(:,isoil) * soil_depth_sl(:,isoil) 
          w_soil_sat_sl(:,isoil) = theta_sat_sl(:,isoil) * soil_depth_sl(:,isoil) 
       ENDWHERE
    ENDDO
    
    !> 7.1 adjust initial water amount to prescribed theta
    !! 
    w_soil_sl(:,:) = w_soil_pwp_sl(:,:) + dsl4jsb_Config(SSE_)%soil_theta_prescribe * (w_soil_fc_sl(:,:) - w_soil_pwp_sl(:,:))

    !> 7.2 soil water potential per layer
    !! 
    WHERE(soil_depth_sl(:,:) > eps8)
       hlp1(:,:)          = MIN(w_soil_fc_sl(:,:),w_soil_sl(:,:)) / soil_depth_sl(:,:)
       w_soil_pot_sl(:,:) = saxtonA(:,:) * hlp1(:,:) ** saxtonB(:,:)
    ELSEWHERE
       w_soil_pot_sl(:,:) = 0.0_wp 
    ENDWHERE


  ! Deallocate local allocatable variables
  DEALLOCATE(hlp1)
  DEALLOCATE(theta_pwp_sl)
  DEALLOCATE(theta_fc_sl)
  DEALLOCATE(theta_sat_sl)

  END SUBROUTINE sse_init_q_

  !-----------------------------------------------------------------------------------------------------
  !> Intialize soil properties - jsbach4
  !! 
  !-----------------------------------------------------------------------------------------------------
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

    dsl4jsb_Get_var3D_onDomain(HYDRO_, w_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, w_ice_sl )
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_field_cap_sl)
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
      !$claw acc update host(heat_cond_sl)
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

  !-----------------------------------------------------------------------------------------------------
  !> ensure a certain minimum proportion for sand, clay and silt - quincy
  !! 
  !! Saxton equations with too small soil-texture fractions may lead to unrealistic percolation rates and soil water potentials
  !!
  !! using these modified soil texture for "saxton functions in sse_init" only
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_soil_texture_saxton_compatible(nc, nsoil_sb, sand_sl, silt_sl, clay_sl)

    USE mo_jsb_math_constants,          ONLY: eps8
    USE mo_kind,                        ONLY: wp
    USE mo_exception,                   ONLY: message


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut 
    INTEGER,                                          INTENT(in)    :: nc, &             !< dimensions
                                                                       nsoil_sb
    REAL(wp), DIMENSION(nc,nsoil_sb),                 INTENT(inout) :: sand_sl, &        !< sand proportion per layer
                                                                       silt_sl, &        !< silt proportion per layer
                                                                       clay_sl           !< clay proportion per layer
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc,nsoil_sb)                :: hlp_sand_sl, &
                                                       hlp_silt_sl, &
                                                       hlp_clay_sl
    REAL(wp), DIMENSION(nc,nsoil_sb)                :: sand_original_sl, &
                                                       silt_original_sl, &
                                                       clay_original_sl
    REAL(wp)                                        :: min_prop_soil_texture, &
                                                       max_prop_soil_texture
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_soil_texture_saxton_compatible'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 0.9 init local variables
    !!  
    min_prop_soil_texture = 0.2_wp
    max_prop_soil_texture = 1.0_wp - (2.0_wp * min_prop_soil_texture)
    hlp_sand_sl(:,:)      = 0.0_wp
    hlp_clay_sl(:,:)      = 0.0_wp
    hlp_silt_sl(:,:)      = 0.0_wp
    sand_original_sl(:,:) = sand_sl(:,:)
    clay_original_sl(:,:) = clay_sl(:,:)
    silt_original_sl(:,:) = silt_sl(:,:)

    !> 1.0 ensure a minimum proportion 'min_prop_soil_texture' for sand, clay and silt
    !! this is to ensure applying the empirical functions by Saxton within the range of the texture they were designed for \n
    !! 
    WHERE (sand_sl(:,:) < min_prop_soil_texture .AND. clay_sl(:,:) < min_prop_soil_texture) 
      sand_sl(:,:)     = min_prop_soil_texture
      clay_sl(:,:)     = min_prop_soil_texture
      silt_sl(:,:)     = max_prop_soil_texture
    ELSEWHERE (sand_sl(:,:) < min_prop_soil_texture .AND. silt_sl(:,:) < min_prop_soil_texture)
      sand_sl(:,:)     = min_prop_soil_texture
      silt_sl(:,:)     = min_prop_soil_texture
      clay_sl(:,:)     = max_prop_soil_texture
    ELSEWHERE (clay_sl(:,:) < min_prop_soil_texture .AND. silt_sl(:,:) < min_prop_soil_texture)
      clay_sl(:,:)     = min_prop_soil_texture
      silt_sl(:,:)     = min_prop_soil_texture
      sand_sl(:,:)     = max_prop_soil_texture
    ELSEWHERE (sand_sl(:,:) < min_prop_soil_texture)
      hlp_sand_sl(:,:) = min_prop_soil_texture - sand_sl(:,:)
      ! proportion of hlp_sand_sl(:,:) to be drawn from clay and silt
      hlp_clay_sl(:,:) = clay_sl(:,:) + (1.0_wp - (clay_sl(:,:) + silt_sl(:,:))) / 2.0_wp
      hlp_silt_sl(:,:) = silt_sl(:,:) + (1.0_wp - (clay_sl(:,:) + silt_sl(:,:))) / 2.0_wp
      ! where new clay proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from silt
      WHERE     ( (clay_sl(:,:) - hlp_sand_sl(:,:) * hlp_clay_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        silt_sl(:,:) = silt_sl(:,:) - hlp_sand_sl(:,:)
      ! where new silt proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from clay
      ELSEWHERE ( (silt_sl(:,:) - hlp_sand_sl(:,:) * hlp_silt_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        clay_sl(:,:) = clay_sl(:,:) - hlp_sand_sl(:,:)
      ELSEWHERE
        clay_sl(:,:) = clay_sl(:,:) - hlp_sand_sl(:,:) * hlp_clay_sl(:,:)
        silt_sl(:,:) = silt_sl(:,:) - hlp_sand_sl(:,:) * hlp_silt_sl(:,:)
      ENDWHERE
      sand_sl(:,:) = min_prop_soil_texture
    ELSEWHERE (clay_sl(:,:) < min_prop_soil_texture)
      hlp_clay_sl(:,:) = min_prop_soil_texture - clay_sl(:,:)
      ! proportion of hlp_clay_sl(:,:) to be drawn from sand and silt
      hlp_sand_sl(:,:) = sand_sl(:,:) + (1.0_wp - (sand_sl(:,:) + silt_sl(:,:))) / 2.0_wp
      hlp_silt_sl(:,:) = silt_sl(:,:) + (1.0_wp - (sand_sl(:,:) + silt_sl(:,:))) / 2.0_wp
      ! where new sand proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from silt
      WHERE     ( (sand_sl(:,:) - hlp_clay_sl(:,:) * hlp_sand_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        silt_sl(:,:) = silt_sl(:,:) - hlp_clay_sl(:,:)
      ! where new silt proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from sand
      ELSEWHERE ( (silt_sl(:,:) - hlp_clay_sl(:,:) * hlp_silt_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        sand_sl(:,:) = sand_sl(:,:) - hlp_clay_sl(:,:)
      ELSEWHERE
        sand_sl(:,:) = sand_sl(:,:) - hlp_clay_sl(:,:) * hlp_sand_sl(:,:)
        silt_sl(:,:) = silt_sl(:,:) - hlp_clay_sl(:,:) * hlp_silt_sl(:,:)
      ENDWHERE
      clay_sl(:,:) = min_prop_soil_texture      
    ELSEWHERE (silt_sl(:,:) < min_prop_soil_texture)
      hlp_silt_sl(:,:) = min_prop_soil_texture - silt_sl(:,:)
      ! proportion of hlp_silt_sl(:,:) to be drawn from sand and clay
      hlp_sand_sl(:,:) = sand_sl(:,:) + (1.0_wp - (sand_sl(:,:) + clay_sl(:,:))) / 2.0_wp
      hlp_clay_sl(:,:) = clay_sl(:,:) + (1.0_wp - (sand_sl(:,:) + clay_sl(:,:))) / 2.0_wp
      ! where new sand proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from clay
      WHERE     ( (sand_sl(:,:) - hlp_silt_sl(:,:) * hlp_sand_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        clay_sl(:,:) = clay_sl(:,:) - hlp_silt_sl(:,:)
      ! where new clay proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from sand
      ELSEWHERE ( (clay_sl(:,:) - hlp_silt_sl(:,:) * hlp_clay_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        sand_sl(:,:) = sand_sl(:,:) - hlp_silt_sl(:,:)
      ELSEWHERE
        sand_sl(:,:) = sand_sl(:,:) - hlp_silt_sl(:,:) * hlp_sand_sl(:,:)
        clay_sl(:,:) = clay_sl(:,:) - hlp_silt_sl(:,:) * hlp_clay_sl(:,:)
      ENDWHERE
      silt_sl(:,:) = min_prop_soil_texture
    ENDWHERE
    
    !> 1.1 ensure all adds up to exactly 1.0_wp
    !! 
    clay_sl(:,:) = 1.0_wp
    clay_sl(:,:) = clay_sl(:,:) - (sand_sl(:,:) + silt_sl(:,:))

    !> 1.2 output any changes in soil texture to stdout or log
    !! 
    IF (ABS(SUM(sand_sl(:,:) - sand_original_sl(:,:))) > eps8) THEN
      CALL message(TRIM(routine), ' Modified the soil texture sand input for saxton calculation.')
      print*,TRIM(routine), " Modified the soil texture sand input for saxton calculation: ", sand_sl(:,:) - sand_original_sl(:,:)
    ENDIF
    IF (ABS(SUM(clay_sl(:,:) - clay_original_sl(:,:))) > eps8) THEN
      CALL message(TRIM(routine), ' Modified the soil texture clay input for saxton calculation.')
      print*,TRIM(routine), " Modified the soil texture clay input for saxton calculation: ", clay_sl(:,:) - clay_original_sl(:,:)
    ENDIF
    IF (ABS(SUM(silt_sl(:,:) - silt_original_sl(:,:))) > eps8) THEN
      CALL message(TRIM(routine), ' Modified the soil texture silt input for saxton calculation.')
      print*,TRIM(routine), " Modified the soil texture silt input for saxton calculation: ", silt_sl(:,:) - silt_original_sl(:,:)
    ENDIF

  END SUBROUTINE calc_soil_texture_saxton_compatible

#endif
END MODULE mo_sse_init
