!>
!! @brief Initialization of the the vegetation memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_veg_init
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_close_file, jsb_netcdf_open_input, &
                                  jsb_netcdf_read_2d_time, jsb_netcdf_read_2d

  dsl4jsb_Use_processes VEG, SOIL
  dsl4jsb_Use_config(veg)
  dsl4jsb_Use_config(soil)
  dsl4jsb_Use_memory(veg)
  dsl4jsb_Use_memory(soil)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: veg_init_bc, veg_init_ic

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_init'

CONTAINS

  SUBROUTINE veg_init_bc(tile)

    TYPE(t_jsb_tile), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model
    dsl4jsb_Def_memory(VEG)

    REAL(wp), POINTER :: lai_mon(:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':veg_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get veg memory of the tile
    dsl4jsb_Get_memory(VEG)

    ! ...

    ! Propagate boundary conditions down the tile structure
    ! TBD

  END SUBROUTINE veg_init_bc

  RECURSIVE SUBROUTINE veg_init_ic(tile)

    USE mo_veg_process, ONLY: foliage_projected_cover
    USE mo_veg_util,    ONLY: water_in_root_zone
    USE mo_soil_util,   ONLY: soil_depth_to_layers_2d
    USE mo_jsb_time,    ONLY: get_time_interpolation_weights

    TYPE(t_jsb_tile),  INTENT(inout) :: tile

    TYPE(t_jsb_model),        POINTER :: model
    dsl4jsb_Def_config(VEG)
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_config(SOIL)
    dsl4jsb_Def_memory(SOIL)

    INTEGER :: i

    REAL(wp), POINTER ::        &
      & lai_mon(:,:,:),         &
      & return_pointer (:,:)       !< temporary pointer

    TYPE(t_input_file) :: input_file

    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':veg_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get veg memory of the tile
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(SOIL)

    ! Get veg config
    dsl4jsb_Get_config(VEG)
    dsl4jsb_Get_config(SOIL)

    IF (.NOT. dsl4jsb_Config(VEG)%active) RETURN

    ! Initialize tile fractions
    ! TBD

    CALL message(TRIM(routine), 'Initializing veg memory for tile '//TRIM(tile%name))

    ! Initialize physical state
    !

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(VEG)%bc_filename), model%id)

    ! Leaf area index
    ! Get monthly LAI climatology
    lai_mon => jsb_netcdf_read_2d_time(    &
      & input_file,                        &
      & variable_name='lai_clim')
    lai_mon(:,:,:) = MERGE(lai_mon(:,:,:), 0._wp, lai_mon(:,:,:) > 0._wp)
    veg__mem%lai_mon(:,:,2:13) = lai_mon(:,:,1:12)
    veg__mem%lai_mon(:,:,1   ) = lai_mon(:,:,12  )
    veg__mem%lai_mon(:,:,14  ) = lai_mon(:,:,1   )
    DEALLOCATE(lai_mon)

    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
    dsl4jsb_var2D_onDomain(VEG, lai) = &
      wgt1 * veg__mem%lai_mon(:,:,nmw1+1) + wgt2 * veg__mem%lai_mon(:,:,nmw2+1)

    ! Forest fraction
    return_pointer => jsb_netcdf_read_2d(        &
      & input_file,                              &
      & variable_name='forest_fract',            &
      & fill_array = dsl4jsb_var2D_ptr(VEG, fract_forest))
    dsl4jsb_var2D_onDomain(VEG, fract_forest) = &
      MERGE(dsl4jsb_var2D_onDomain(VEG, fract_forest), 0._wp, dsl4jsb_var2D_onDomain(VEG, fract_forest) >= 0._wp)

    dsl4jsb_var2D_onDomain(VEG, canopy_cond) = 0._wp

    dsl4jsb_var2D_onDomain(VEG, skin_can_depth) = dsl4jsb_Config(VEG)%wsr_capacity * dsl4jsb_var2D_onDomain(VEG, lai)

    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='albedo_veg_vis',           &
      & fill_array = dsl4jsb_var2D_ptr(VEG, alb_vis_can))
    dsl4jsb_var2D_onDomain(VEG, alb_vis_can) = &
      MERGE(dsl4jsb_var2D_onDomain(VEG, alb_vis_can), 0._wp, dsl4jsb_var2D_onDomain(VEG, alb_vis_can) >= 0._wp)

    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='albedo_veg_nir',           &
      & fill_array = dsl4jsb_var2D_ptr(VEG, alb_nir_can))
    dsl4jsb_var2D_onDomain(VEG, alb_nir_can) = &
      MERGE(dsl4jsb_var2D_onDomain(VEG, alb_nir_can), 0._wp, dsl4jsb_var2D_onDomain(VEG, alb_nir_can) >= 0._wp)

    CALL jsb_netcdf_close_file(input_file)
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SOIL)%bc_filename), model%id)

    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='root_depth',               &
      & fill_array = dsl4jsb_var2D_ptr(VEG, root_depth))
    dsl4jsb_var2D_onDomain(VEG, root_depth) = &
      MERGE(dsl4jsb_var2D_onDomain(VEG, root_depth), 0._wp, dsl4jsb_var2D_onDomain(VEG, root_depth) >= 0._wp)
    dsl4jsb_var2D_onDomain(VEG, root_depth) = &           ! Temporary fix (root depth sometimes larger then soil depth)
      MIN(dsl4jsb_var2D_onDomain(VEG, root_depth), dsl4jsb_var2D_onDomain(SOIL, soil_depth))
!!$    dsl4jsb_var2D_onDomain(VEG, root_depth) = soil_mem%soil_depth%ptr(:,:) * 0.75_wp  ! Temporary

    dsl4jsb_var3D_onDomain(VEG, root_depth_l) = soil_depth_to_layers_2d(               &
      & dsl4jsb_var2D_onDomain(VEG, root_depth),                                         & ! Total rooting depth
      & dsl4jsb_Config(SOIL)%dz_water(1:dsl4jsb_Config(SOIL)%nsoil_water))    ! Soil layer thicknesses from namelist

    CALL jsb_netcdf_close_file(input_file)                              ! Temporary
    input_file = jsb_netcdf_open_input('bc_land_frac.nc', model%id)     ! Temporary

    ! Maximum foliage projected cover
    return_pointer => jsb_netcdf_read_2d(        &
      & input_file,                              &
      & variable_name='veg_ratio_max',           &
      & fill_array = dsl4jsb_var2D_ptr(VEG, fract_fpc_max))
    dsl4jsb_var2D_onDomain(VEG, fract_fpc_max) = &
      MERGE(dsl4jsb_var2D_onDomain(VEG, fract_fpc_max), 0._wp, dsl4jsb_var2D_onDomain(VEG, fract_fpc_max) >= 0._wp)

!!$    veg_mem%fract_fpc(:,:) = foliage_projected_cover(veg_mem%fract_fpc_max%ptr(:,:), veg_mem%lai%ptr(:,:), veg_mem%clump%ptr(:,:))
    dsl4jsb_var2D_onDomain(VEG, fract_fpc) = &
      foliage_projected_cover(dsl4jsb_var2D_onDomain(VEG, fract_fpc_max), dsl4jsb_var2D_onDomain(VEG, lai), &
      SPREAD(SPREAD(2._wp, DIM=1, NCOPIES=SIZE(dsl4jsb_var2D_ptr(VEG, lai),1)), DIM=2, NCOPIES=SIZE(dsl4jsb_var2D_ptr(VEG, lai),2))) ! Temporary

    CALL jsb_netcdf_close_file(input_file)

    ! Initialize water state
    !
    dsl4jsb_var2D_onDomain(VEG, wfract_can) = 0._wp
    dsl4jsb_var2D_onDomain(VEG, wsr_can)    = 0._wp
    dsl4jsb_var2D_onDomain(VEG, sfract_can) = 0._wp
    dsl4jsb_var2D_onDomain(VEG, wsn_can)    = 0._wp

    ! Calculate water content and field capacity in the root zone
    DO i=1,SIZE(dsl4jsb_var2D_ptr(VEG, ws_root), 2)
      dsl4jsb_var2D_ptr(VEG, ws_root)     (:,i) = water_in_root_zone( dsl4jsb_var3D_ptr(SOIL, ws_l)     (:,:,i), &
        & dsl4jsb_var3D_ptr(SOIL, soil_depth_l)(:,:,i), dsl4jsb_var3D_ptr(VEG, root_depth_l)(:,:,i))
      dsl4jsb_var2D_ptr(VEG, ws_fc_root)  (:,i) = water_in_root_zone( dsl4jsb_var3D_ptr(SOIL, ws_fc_l)  (:,:,i), &
        & dsl4jsb_var3D_ptr(SOIL, soil_depth_l)(:,:,i), dsl4jsb_var3D_ptr(VEG, root_depth_l)(:,:,i))
      dsl4jsb_var2D_ptr(VEG, ws_pwp_root) (:,i) = water_in_root_zone( dsl4jsb_var3D_ptr(SOIL, ws_pwp_l) (:,:,i), &
        & dsl4jsb_var3D_ptr(SOIL, soil_depth_l)(:,:,i), dsl4jsb_var3D_ptr(VEG, root_depth_l)(:,:,i))
    END DO
!!$    veg_mem%ws_root%ptr    (:,:) = 0._wp
!!$    veg_mem%ws_fc_root%ptr (:,:) = 0._wp
!!$    veg_mem%ws_pwp_root%ptr(:,:) = 0._wp

    dsl4jsb_var2d_onDomain(VEG, water_stress) = 0._wp

    ! Propagate initial conditions down the tile structure, starting with the tile with which this subroutine was called initially
    ! TBD: sub-tiles should not read the netcdf files again
!!$    IF (tile%no_of_sub_tiles > 0) THEN
!!$      DO i=1,tile%no_of_sub_tiles
!!$        CALL veg_init_ic(tile%sub_tiles(i))
!!$      END DO
!!$    END IF

  END SUBROUTINE veg_init_ic

#endif
END MODULE mo_veg_init
