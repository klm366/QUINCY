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

  dsl$jsb Use_processes VEG, SOIL
  dsl$jsb Use_config(veg)
  dsl$jsb Use_config(soil)
  dsl$jsb Use_memory(veg)
  dsl$jsb Use_memory(soil)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: veg_init_bc, veg_init_ic

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_init'

CONTAINS

  SUBROUTINE veg_init_bc(tile)

    TYPE(t_jsb_tile), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model
    dsl$jsb Def_memory(veg_mem, veg)

    REAL(wp), POINTER :: lai_mon(:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':veg_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get veg memory of the tile
    dsl$jsb Get_memory(veg_mem, VEG)

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
    dsl$jsb Def_config(veg_conf, veg)
    dsl$jsb Def_memory(veg_mem, veg)
    dsl$jsb Def_config(soil_conf, soil)
    dsl$jsb Def_memory(soil_mem, soil)

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
    dsl$jsb Get_memory(veg_mem, VEG)
    dsl$jsb Get_memory(soil_mem, SOIL)

    ! Get veg config
    dsl$jsb Get_config(veg_conf, VEG)
    dsl$jsb Get_config(soil_conf, SOIL)

    IF (.NOT. veg_conf%active) RETURN

    ! Initialize tile fractions
    ! TBD

    CALL message(TRIM(routine), 'Initializing veg memory for tile '//TRIM(tile%name))

    ! Initialize physical state
    !

    input_file = jsb_netcdf_open_input(TRIM(veg_conf%bc_filename), model%id)

    ! Leaf area index
    ! Get monthly LAI climatology
    lai_mon => jsb_netcdf_read_2d_time(    &
      & input_file,                        &
      & variable_name='lai_clim')
    lai_mon(:,:,:) = MERGE(lai_mon(:,:,:), 0._wp, lai_mon(:,:,:) > 0._wp)
    veg_mem%lai_mon(:,:,2:13) = lai_mon(:,:,1:12)
    veg_mem%lai_mon(:,:,1   ) = lai_mon(:,:,12  )
    veg_mem%lai_mon(:,:,14  ) = lai_mon(:,:,1   )
    DEALLOCATE(lai_mon)

    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
    veg_mem%lai%ptr(:,:) = wgt1 * veg_mem%lai_mon(:,:,nmw1+1) + wgt2 * veg_mem%lai_mon(:,:,nmw2+1)

    ! Forest fraction
    return_pointer => jsb_netcdf_read_2d(        &
      & input_file,                              &
      & variable_name='forest_fract',            &
      & fill_array = veg_mem%fract_forest%ptr)
    veg_mem%fract_forest%ptr(:,:) = MERGE(veg_mem%fract_forest%ptr(:,:), 0._wp, veg_mem%fract_forest%ptr(:,:) >= 0._wp)

    veg_mem%canopy_cond%ptr(:,:) = 0._wp

    veg_mem%skin_can_depth%ptr(:,:) = veg_conf%wsr_capacity * veg_mem%lai%ptr(:,:)

    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='albedo_veg_vis',           &
      & fill_array = veg_mem%alb_vis_can%ptr)
    veg_mem%alb_vis_can%ptr(:,:) = MERGE(veg_mem%alb_vis_can%ptr(:,:), 0._wp, veg_mem%alb_vis_can%ptr(:,:) >= 0._wp)

    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='albedo_veg_nir',           &
      & fill_array = veg_mem%alb_nir_can%ptr)
    veg_mem%alb_nir_can%ptr(:,:) = MERGE(veg_mem%alb_nir_can%ptr(:,:), 0._wp, veg_mem%alb_nir_can%ptr(:,:) >= 0._wp)

    CALL jsb_netcdf_close_file(input_file)
    input_file = jsb_netcdf_open_input(TRIM(soil_conf%bc_filename), model%id)

    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='root_depth',               &
      & fill_array = veg_mem%root_depth%ptr)
    veg_mem%root_depth%ptr(:,:) = MERGE(veg_mem%root_depth%ptr(:,:), 0._wp, veg_mem%root_depth%ptr(:,:) >= 0._wp)
    veg_mem%root_depth%ptr(:,:) = MIN(veg_mem%root_depth%ptr(:,:), soil_mem%soil_depth%ptr(:,:))  ! Temporary fix (root depth sometimes larger then soil depth)
!!$    veg_mem%root_depth%ptr(:,:) = soil_mem%soil_depth%ptr(:,:) * 0.75_wp  ! Temporary

    veg_mem%root_depth_l%ptr(:,:,:) = soil_depth_to_layers_2d(               &
      & veg_mem%root_depth%ptr(:,:),                                         & ! Total rooting depth
      & soil_conf%dz_water(1:soil_conf%nsoil_water))    ! Soil layer thicknesses from namelist

    CALL jsb_netcdf_close_file(input_file)                              ! Temporary
    input_file = jsb_netcdf_open_input('bc_land_frac.nc', model%id)     ! Temporary

    ! Maximum foliage projected cover
    return_pointer => jsb_netcdf_read_2d(        &
      & input_file,                              &
      & variable_name='veg_ratio_max',           &
      & fill_array = veg_mem%fract_fpc_max%ptr)
    veg_mem%fract_fpc_max%ptr(:,:) = MERGE(veg_mem%fract_fpc_max%ptr(:,:), 0._wp, veg_mem%fract_fpc_max%ptr(:,:) >= 0._wp)

!!$    veg_mem%fract_fpc(:,:) = foliage_projected_cover(veg_mem%fract_fpc_max%ptr(:,:), veg_mem%lai%ptr(:,:), veg_mem%clump%ptr(:,:))
    veg_mem%fract_fpc%ptr(:,:) = foliage_projected_cover(veg_mem%fract_fpc_max%ptr(:,:), veg_mem%lai%ptr(:,:), &
      SPREAD(SPREAD(2._wp, DIM=1, NCOPIES=SIZE(veg_mem%lai%ptr,1)), DIM=2, NCOPIES=SIZE(veg_mem%lai%ptr,2))) ! Temporary

    CALL jsb_netcdf_close_file(input_file)

    ! Initialize water state
    !
    veg_mem%wfract_can%ptr  (:,:) = 0._wp
    veg_mem%wsr_can%ptr     (:,:) = 0._wp
    veg_mem%sfract_can%ptr  (:,:) = 0._wp
    veg_mem%wsn_can%ptr     (:,:) = 0._wp

    ! Calculate water content and field capacity in the root zone
    DO i=1,SIZE(veg_mem%ws_root%ptr,2)
      veg_mem%ws_root%ptr     (:,i) = water_in_root_zone( soil_mem%ws_l%ptr     (:,:,i), &
        & soil_mem%soil_depth_l%ptr(:,:,i), veg_mem%root_depth_l%ptr(:,:,i))
      veg_mem%ws_fc_root%ptr  (:,i) = water_in_root_zone( soil_mem%ws_fc_l%ptr  (:,:,i), &
        & soil_mem%soil_depth_l%ptr(:,:,i), veg_mem%root_depth_l%ptr(:,:,i))
      veg_mem%ws_pwp_root%ptr (:,i) = water_in_root_zone( soil_mem%ws_pwp_l%ptr (:,:,i), &
        & soil_mem%soil_depth_l%ptr(:,:,i), veg_mem%root_depth_l%ptr(:,:,i))
    END DO
!!$    veg_mem%ws_root%ptr    (:,:) = 0._wp
!!$    veg_mem%ws_fc_root%ptr (:,:) = 0._wp
!!$    veg_mem%ws_pwp_root%ptr(:,:) = 0._wp

    veg_mem%water_stress%ptr(:,:) = 0._wp

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
