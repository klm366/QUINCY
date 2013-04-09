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

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile
  USE mo_veg_memory_class,  ONLY: t_veg_memory, t_veg_phys_state, t_veg_water_state
  USE mo_jsb_base,          ONLY: get_model
  USE mo_netcdf_read,       ONLY: netcdf_read_oncells_2d_time, netcdf_read_oncells_2d
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: veg_init_bc, veg_init_ic

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_init'

CONTAINS

  SUBROUTINE veg_init_bc(tile)

    TYPE(t_jsb_tile), INTENT(in) :: tile

    TYPE(t_jsb_model),       POINTER :: model
    TYPE(t_veg_memory), POINTER :: veg

    REAL(wp), POINTER :: lai_mon(:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':veg_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get veg memory of the tile
    veg => tile%mem%veg

    ! ...

    ! Propagate boundary conditions down the tile structure
    ! TBD

  END SUBROUTINE veg_init_bc

  RECURSIVE SUBROUTINE veg_init_ic(tile)

    USE mo_veg_process, ONLY: foliage_projected_cover
    USE mo_soil_util, ONLY: soil_depth_to_layers_2d

    TYPE(t_jsb_tile),  INTENT(in) :: tile

    TYPE(t_jsb_model),       POINTER :: model
    TYPE(t_veg_memory),      POINTER :: veg
    TYPE(t_veg_phys_state),  POINTER :: phys
    TYPE(t_veg_water_state), POINTER :: water

    INTEGER :: i

    REAL(wp), POINTER ::        &
      & lai_mon(:,:,:),         &
      & return_pointer (:,:)       !< temporary pointer

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':veg_init_ic'

    model => get_model(tile%owner_model_id)

    ! Initialize tile fractions
    ! TBD

    ! Get veg memory of the tile
    veg => tile%mem%veg

    IF (.NOT. ASSOCIATED(veg)) RETURN

    CALL message(TRIM(routine), 'Initializing veg memory for tile '//TRIM(tile%name))

    ! Initialize physical state
    !
    phys => veg%phys_state

    ! TBD - For now set LAI to constant
    phys%lai(:,:) = 5._wp  ! Temporary

    ! Maximum foliage projected cover
!!$    return_pointer => netcdf_read_oncells_2d(        &
!!$      & filename=model%config%veg%bc_filename,       &
!!$      & variable_name='veg_ratio_max',               &
!!$      & fill_array = phys%fract_fpc_max,             &
!!$      & patch=model%grid%patch)
!!$    phys%fract_fpc_max(:,:) = MERGE(phys%fract_fpc_max(:,:), 0._wp, phys%fract_fpc_max(:,:) >= 0._wp)
    phys%fract_fpc_max(:,:) = 1._wp  ! Temporary

!!$    phys%fract_fpc(:,:) = foliage_projected_cover(phys%fract_fpc_max(:,:), phys%lai(:,:), phys%clump(:,:))
    phys%fract_fpc(:,:)     = 1._wp  ! Temporary
    phys%fract_fpc_max(:,:) = MERGE(phys%fract_fpc_max(:,:), 0._wp, phys%fract_fpc_max(:,:) >= 0._wp)

    ! Forest fraction
    return_pointer => netcdf_read_oncells_2d(        &
      & filename=model%config%veg%bc_filename,       &
      & variable_name='forest_fract',                &
      & fill_array = phys%fract_forest,              &
      & patch=model%grid%patch)
    phys%fract_forest(:,:) = MERGE(phys%fract_forest(:,:), 0._wp, phys%fract_forest(:,:) >= 0._wp)

    phys%canopy_cond(:,:) = 0._wp

    phys%skin_can_depth(:,:) = model%config%veg%wsr_capacity * phys%lai(:,:)

    phys%root_depth(:,:) = 2.65355_wp  ! Temporary

    phys%root_depth_l(:,:,:) = soil_depth_to_layers_2d(               &
      & phys%root_depth(:,:),                                         & ! Total rooting depth
      & model%config%soil%dz_water(1:model%config%soil%nsoil_water))    ! Soil layer thicknesses from namelist

    return_pointer => netcdf_read_oncells_2d(         &
      & filename=model%config%veg%bc_filename,       &
      & variable_name='albedo_veg_vis',               &
      & fill_array = phys%alb_vis_can,                &
      & patch=model%grid%patch)
    phys%alb_vis_can(:,:) = MERGE(phys%alb_vis_can(:,:), 0._wp, phys%alb_vis_can(:,:) >= 0._wp)

    return_pointer => netcdf_read_oncells_2d(         &
      & filename=model%config%veg%bc_filename,       &
      & variable_name='albedo_veg_nir',               &
      & fill_array = phys%alb_nir_can,                &
      & patch=model%grid%patch)
    phys%alb_nir_can(:,:) = MERGE(phys%alb_nir_can(:,:), 0._wp, phys%alb_nir_can(:,:) >= 0._wp)

    ! Initialize water state
    !
    water => veg%water_state

    water%wfract_can  (:,:) = 0._wp
    water%wsr_can     (:,:) = 0._wp
    water%sfract_can  (:,:) = 0._wp
    water%wsn_can     (:,:) = 0._wp
    water%ws_root     (:,:) = 0._wp
    water%ws_fc_root  (:,:) = 0._wp
    water%ws_pwp_root (:,:) = 0._wp

    water%water_stress(:,:) = 0._wp

    ! Propagate initial conditions down the tile structure, starting with the tile with which this subroutine was called initially
    ! TBD: sub-tiles should not read the netcdf files again
    IF (tile%no_of_sub_tiles > 0) THEN
      DO i=1,tile%no_of_sub_tiles
        CALL veg_init_ic(tile%sub_tiles(i))
      END DO
    END IF

  END SUBROUTINE veg_init_ic

END MODULE mo_veg_init
