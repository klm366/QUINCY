!>
!! @brief Contains the interfaces to the vegetation processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_veg_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_jsb_control,        ONLY: debug_on
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish
  USE mo_util,               ONLY: int2string

  USE mo_jsb_parallel_iface, ONLY: get_my_global_mpi_id

  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_memory_class,   ONLY: t_jsb_memory, t_jsb_memory_p
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task

  dsl4jsb_Use_processes VEG, SOIL, SRF, ATM2LAND
  dsl4jsb_Use_config(VEG)
  dsl4jsb_Use_config(SOIL)

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl4jsb_Use_memory(atm2land)
  dsl4jsb_Use_memory(soil)
  dsl4jsb_Use_memory(srf)
  dsl4jsb_Use_memory(veg)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_vegetation_1, update_vegetation_2
  PUBLIC :: t_veg_process

  TYPE, EXTENDS(t_jsb_process) :: t_veg_process
  CONTAINS
    PROCEDURE :: Register_tasks => Register_veg_tasks
  END TYPE t_veg_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_interface'

CONTAINS

  SUBROUTINE Register_veg_tasks(this, model_id)

    CLASS(t_veg_process), INTENT(inout) :: this
    INTEGER, INTENT(in) :: model_id

    !CALL this%Register_task(tsk_radiation_surface(name="radiation_surface", process_id=VEG))

  END SUBROUTINE Register_veg_tasks

  SUBROUTINE update_vegetation_1(tile, config, iblk, ics, ice, nc)

    USE mo_veg_process, ONLY: canopy_conductance, water_stress_factor
    USE mo_veg_util,    ONLY: water_in_root_zone
    USE mo_phy_schemes, ONLY: qsat


    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p), INTENT(in)    :: config(:)
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc

    dsl4jsb_Def_config(VEG)
    dsl4jsb_Def_config(SOIL)
    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(ATM2LAND)

    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: swpar_srf_down
    dsl4jsb_Real2D_onChunk :: canopy_cond
    dsl4jsb_Real3D_onChunk :: soil_depth_l
    dsl4jsb_Real3D_onChunk :: root_depth_l
    dsl4jsb_Real2D_onChunk :: ws              !< Water content in soil column
    dsl4jsb_Real3D_onChunk :: ws_l            !< Water content in soil layers
    dsl4jsb_Real3D_onChunk :: ws_fc_l         !< Water content of soil layers at field capacity
    dsl4jsb_Real3D_onChunk :: ws_pwp_l        !< Water content of soil layers at permanent wilting point
    dsl4jsb_Real2D_onChunk :: ws_root         !< Water content in root zone
    dsl4jsb_Real2D_onChunk :: ws_fc_root      !< Water content in root zone at field capacity
    dsl4jsb_Real2D_onChunk :: ws_pwp_root     !< Water content in root zone at permanent wilting point
    dsl4jsb_Real2D_onChunk :: water_stress
    dsl4jsb_Real2D_onChunk :: t_srf
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: press_srf

    INTEGER :: nsoil

    INTEGER :: p_pe

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_vegetation_1'

    IF (.NOT. tile%Is_process_active(VEG)) RETURN

!!$    p_pe = get_my_global_mpi_id()
!!$    IF (iblk == 1) CALL message(TRIM(routine), 'Starting update on PE '//int2string(p_pe), all_print=.TRUE.)
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of vegetation_1 ')

    dsl4jsb_Get_config(VEG)
    dsl4jsb_Get_config(SOIL)

    nsoil = dsl4jsb_Config(SOIL)%nsoil_water

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(ATM2LAND)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swpar_srf_down, IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_air,          IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, press_srf,      IN)

    dsl4jsb_Get_var2D_onChunk(VEG, lai,          OUT)
    dsl4jsb_Get_var2D_onChunk(VEG, canopy_cond,  OUT)
    dsl4jsb_Get_var3D_onChunk(VEG, root_depth_l, OUT)
    dsl4jsb_Get_var2D_onChunk(VEG, ws_root,      OUT)
    dsl4jsb_Get_var2D_onChunk(VEG, ws_fc_root,   OUT)
    dsl4jsb_Get_var2D_onChunk(VEG, ws_pwp_root,  OUT)
    dsl4jsb_Get_var2D_onChunk(VEG, water_stress, OUT)

    dsl4jsb_Get_var3D_onChunk(SOIL, soil_depth_l, OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, ws,           OUT)
    dsl4jsb_Get_var3D_onChunk(SOIL, ws_l,         OUT)
    dsl4jsb_Get_var3D_onChunk(SOIL, ws_fc_l,      OUT)
    dsl4jsb_Get_var3D_onChunk(SOIL, ws_pwp_l,     OUT)
    dsl4jsb_Get_var2D_onChunk(SRF,  t_srf,        OUT)

    ! Compute (max) canopy (stomatal) conductance with no water stress
    SELECT CASE(TRIM(dsl4jsb_Config(VEG)%scheme_stom_cond))
      CASE('echam5')
        canopy_cond(:) = canopy_conductance( lai(:), swpar_srf_down(:) )

!!$      CASE('bethy')
!!$        ! TBD: Call bethy to compute unstressed canopy conductance

      CASE default
        CALL finish( TRIM(routine), 'Canopy conductance scheme '//TRIM(dsl4jsb_Config(VEG)%scheme_stom_cond)//&
                   & 'not supported.')
    END SELECT

    ! Get water stress

    ! Calculate water content and field capacity in the root zone
    ws_root   (:)  = water_in_root_zone(ws_l     (:,:), soil_depth_l(:,:), root_depth_l(:,:))
    ws_fc_root(:)  = water_in_root_zone(ws_fc_l  (:,:), soil_depth_l(:,:), root_depth_l(:,:))
    ws_pwp_root(:) = water_in_root_zone(ws_pwp_l (:,:), soil_depth_l(:,:), root_depth_l(:,:))

    ! Calculate water stress factor
    ! Temporary: later use the distributed permanent wilting point from above and change water_stress_factor function
    water_stress(:) = water_stress_factor(ws_root(:), ws_fc_root(:), ws_pwp_root(:))

    ! Compute (actual) canopy (stomatal) conductance under water stress.
    canopy_cond(:) = &
      canopy_conductance( canopy_cond(:),                                & ! in, unstressed canopy conductance
                        & water_stress(:),                               & ! in, water stress factor
                        & q_air(:) > qsat(t_srf(:), press_srf(:))     & ! in, atmosphere saturated?
                        )

    ! Compute skin reservoir fraction on canopy

  END SUBROUTINE update_vegetation_1

  SUBROUTINE update_vegetation_2(tile, config, iblk, ics, ice, nc)

    USE mo_veg_process, ONLY: canopy_albedo, foliage_projected_cover
    USE mo_jsb_time,    ONLY: get_time_interpolation_weights

    REAL(wp), PARAMETER :: &
      albedo_snow_min = 0.4_wp,  & !< Minimum possible snow albedo   ! Temporary, will depend on vegetation type
      albedo_snow_max = 0.8_wp     !< Maximum possible snow albedo   ! Temporary

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p), INTENT(in)    :: config(:)
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc

    dsl4jsb_Real2D_onChunk :: t_srf
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: fract_fpc
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: fract_forest
!!$   dsl4jsb_Real2D_onChunk :: alb_vis_can
!!$   dsl4jsb_Real2D_onChunk :: alb_nir_can
!!$   dsl4jsb_Real2D_onChunk :: alb_vis_soil
!!$   dsl4jsb_Real2D_onChunk :: alb_nir_soil
    dsl4jsb_Real2D_onChunk :: alb_vis_srf
    dsl4jsb_Real2D_onChunk :: alb_nir_srf
    dsl4jsb_Real2D_onChunk :: albedo_srf
    dsl4jsb_Real2D_onChunk :: sfract_can
    dsl4jsb_Real2D_onChunk :: sfract_soil
    REAL(wp), POINTER :: lai_mon       (:,:)

    INTEGER :: p_pe

    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2

    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_vegetation_2'

    IF (.NOT. tile%Is_process_active(VEG)) RETURN

    p_pe = get_my_global_mpi_id()
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)

    dsl4jsb_Get_var2D_onChunk( VEG,  lai,           OUT)
    dsl4jsb_Get_var2D_onChunk( VEG,  fract_fpc,     OUT)
    dsl4jsb_Get_var2D_onChunk( VEG,  fract_fpc_max, OUT)
    dsl4jsb_Get_var2D_onChunk( VEG,  fract_forest,  OUT)
    dsl4jsb_Get_var2D_onChunk( VEG,  sfract_can,    OUT)
!!$    alb_vis_can  => tile%mem%veg  %alb_vis_can%ptr     (ics:ice,    iblk)
!!$    alb_nir_can  => tile%mem%veg  %alb_nir_can%ptr     (ics:ice,    iblk)

!!$    alb_vis_soil => tile%mem%soil %alb_vis_soil%ptr    (ics:ice,    iblk)
!!$    alb_nir_soil => tile%mem%soil %alb_nir_soil%ptr    (ics:ice,    iblk)
    dsl4jsb_Get_var2D_onChunk( SRF, t_srf,       OUT)
    dsl4jsb_Get_var2D_onChunk( SRF, alb_vis_srf, OUT)
    dsl4jsb_Get_var2D_onChunk( SRF, alb_nir_srf, OUT)
    dsl4jsb_Get_var2D_onChunk( SRF, albedo_srf,  OUT)

    dsl4jsb_Get_var2D_onChunk( SOIL, sfract_soil, OUT)

    lai_mon      => veg__mem%lai_mon (ics:ice,    iblk, :)

    ! Update leaf area index
    !
    ! At the moment by interpolating monthly climatological LAI
    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
    lai(:) = REAL(wgt1,wp) * lai_mon(:,nmw1+1) + REAL(wgt2,wp) * lai_mon(:,nmw2+1)

    ! Update foliage projected cover
    !
    ! TBD: Use clumpiness factor instead of 1._wp
    fract_fpc(:) = foliage_projected_cover(fract_fpc_max(:), lai(:), SPREAD(1._wp, DIM=1, NCOPIES=SIZE(lai,1)))

    ! Update albedo
    ! TBD: use newer albedo scheme
    alb_vis_srf(:) = canopy_albedo(t_srf(:), fract_forest(:), sfract_soil(:), sfract_can(:), lai(:), &
                   & SPREAD(albedo_snow_min, NCOPIES=nc, DIM=1), SPREAD(albedo_snow_max, NCOPIES=nc, DIM=1), albedo_srf(:))
    alb_nir_srf(:) = alb_vis_srf(:)

  END SUBROUTINE update_vegetation_2

#endif
END MODULE mo_veg_interface
