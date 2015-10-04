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

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_parallel_iface, ONLY: get_my_global_mpi_id

  dsl$jsb Use_processes VEG, SOIL, SRF, ATM2LAND
  dsl$jsb Use_config(veg)
  dsl$jsb Use_config(soil)

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl$jsb Use_memory(atm2land)
  dsl$jsb Use_memory(soil)
  dsl$jsb Use_memory(srf)
  dsl$jsb Use_memory(veg)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_vegetation_1, update_vegetation_2

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_interface'

CONTAINS

  SUBROUTINE update_vegetation_1( tile, iblk, ics, ice, nc)

    USE mo_veg_process, ONLY: canopy_conductance, water_stress_factor
    USE mo_veg_util,    ONLY: water_in_root_zone
    USE mo_phy_schemes, ONLY: qsat


    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER    :: model

    dsl$jsb Def_config(veg_conf, veg)
    dsl$jsb Def_config(soil_conf, soil)
    dsl$jsb Def_memory(soil_mem, soil)
    dsl$jsb Def_memory(srf_mem, srf)
    dsl$jsb Def_memory(veg_mem, veg)
    dsl$jsb Def_memory(atm2land_mem, atm2land)

    dsl$jsb onChunk_Real2D :: lai
    dsl$jsb onChunk_Real2D :: swpar_srf_down
    dsl$jsb onChunk_Real2D :: canopy_cond
    dsl$jsb onChunk_Real3D :: soil_depth_l
    dsl$jsb onChunk_Real3D :: root_depth_l
    dsl$jsb onChunk_Real2D :: ws              !< Water content in soil column
    dsl$jsb onChunk_Real3D :: ws_l            !< Water content in soil layers
    dsl$jsb onChunk_Real3D :: ws_fc_l         !< Water content of soil layers at field capacity
    dsl$jsb onChunk_Real3D :: ws_pwp_l        !< Water content of soil layers at permanent wilting point
    dsl$jsb onChunk_Real2D :: ws_root         !< Water content in root zone
    dsl$jsb onChunk_Real2D :: ws_fc_root      !< Water content in root zone at field capacity
    dsl$jsb onChunk_Real2D :: ws_pwp_root     !< Water content in root zone at permanent wilting point
    dsl$jsb onChunk_Real2D :: water_stress
    dsl$jsb onChunk_Real2D :: t_srf
    dsl$jsb onChunk_Real2D :: q_air
    dsl$jsb onChunk_Real2D :: press_srf

    INTEGER :: nsoil

    INTEGER :: p_pe

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_vegetation_1'

    IF (.NOT. tile%Is_process_active(VEG)) RETURN

!!$    p_pe = get_my_global_mpi_id()
!!$    IF (iblk == 1) CALL message(TRIM(routine), 'Starting update on PE '//int2string(p_pe), all_print=.TRUE.)
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of vegetation_1 ')

    model => get_model(tile%owner_model_id)

    dsl$jsb Get_config(veg_conf, VEG)
    dsl$jsb Get_config(soil_conf, SOIL)

    nsoil = soil_conf%nsoil_water

    ! Get reference to variables for current block
    !
    dsl$jsb Get_memory(veg_mem, VEG)
    dsl$jsb Get_memory(srf_mem, SRF)
    dsl$jsb Get_memory(soil_mem, SOIL)
    dsl$jsb Get_memory(atm2land_mem, ATM2LAND)

    lai                => veg_mem  %lai%ptr            (ics:ice,    iblk)
    swpar_srf_down     => atm2land_mem     %swpar_srf_down%ptr (ics:ice,    iblk)
    canopy_cond        => veg_mem  %canopy_cond%ptr    (ics:ice,    iblk)
    soil_depth_l       => soil_mem %soil_depth_l%ptr   (ics:ice, :, iblk)
    root_depth_l       => veg_mem  %root_depth_l%ptr   (ics:ice, :, iblk)
    ws                 => soil_mem %ws%ptr             (ics:ice,    iblk)
    ws_l               => soil_mem %ws_l%ptr           (ics:ice, :, iblk)
    ws_fc_l            => soil_mem %ws_fc_l%ptr        (ics:ice, :, iblk)
    ws_pwp_l           => soil_mem %ws_pwp_l%ptr       (ics:ice, :, iblk)
    ws_root            => veg_mem  %ws_root%ptr        (ics:ice,    iblk)
    ws_fc_root         => veg_mem  %ws_fc_root%ptr     (ics:ice,    iblk)
    ws_pwp_root        => veg_mem  %ws_pwp_root%ptr    (ics:ice,    iblk)
    water_stress       => veg_mem  %water_stress%ptr   (ics:ice,    iblk)
    t_srf              => srf_mem  %t_srf%ptr          (ics:ice,    iblk)
    q_air              => atm2land_mem     %q_air%ptr          (ics:ice,    iblk)
    press_srf          => atm2land_mem     %press_srf%ptr      (ics:ice,    iblk)

    ! Compute (max) canopy (stomatal) conductance with no water stress
    SELECT CASE(TRIM(veg_conf%scheme_stom_cond))
      CASE('echam5')
        canopy_cond(:) = canopy_conductance( lai(:), swpar_srf_down(:) )

!!$      CASE('bethy')
!!$        ! TBD: Call bethy to compute unstressed canopy conductance

      CASE default
        CALL finish( TRIM(routine), 'Canopy conductance scheme '//TRIM(veg_conf%scheme_stom_cond)//&
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

  SUBROUTINE update_vegetation_2( tile, iblk, ics, ice, nc)

    USE mo_veg_process, ONLY: canopy_albedo, foliage_projected_cover
    USE mo_jsb_time,    ONLY: get_time_interpolation_weights

    REAL(wp), PARAMETER :: &
      albedo_snow_min = 0.4_wp,  & !< Minimum possible snow albedo   ! Temporary, will depend on vegetation type
      albedo_snow_max = 0.8_wp     !< Maximum possible snow albedo   ! Temporary

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER    :: model

    dsl$jsb onChunk_Real2D :: t_srf
    dsl$jsb onChunk_Real2D :: lai
    dsl$jsb onChunk_Real2D :: fract_fpc
    dsl$jsb onChunk_Real2D :: fract_fpc_max
    dsl$jsb onChunk_Real2D :: fract_forest
!!$   dsl$jsb onChunk_Real2D :: alb_vis_can
!!$   dsl$jsb onChunk_Real2D :: alb_nir_can
!!$   dsl$jsb onChunk_Real2D :: alb_vis_soil
!!$   dsl$jsb onChunk_Real2D :: alb_nir_soil
    dsl$jsb onChunk_Real2D :: alb_vis_srf
    dsl$jsb onChunk_Real2D :: alb_nir_srf
    dsl$jsb onChunk_Real2D :: albedo_srf
    dsl$jsb onChunk_Real2D :: sfract_can
    dsl$jsb onChunk_Real2D :: sfract_soil
    REAL(wp), POINTER :: lai_mon       (:,:)

    INTEGER :: p_pe

    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2

    dsl$jsb Def_memory(veg_mem, veg)
    dsl$jsb Def_memory(srf_mem, srf)
    dsl$jsb Def_memory(soil_mem, soil)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_vegetation_2'

    IF (.NOT. tile%Is_process_active(VEG)) RETURN

    p_pe = get_my_global_mpi_id()
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl$jsb Get_memory(veg_mem, VEG)
    dsl$jsb Get_memory(srf_mem, SRF)
    dsl$jsb Get_memory(soil_mem, SOIL)

    dsl$jsb Get_var2D_onChunk( veg_mem,  lai,           OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,  fract_fpc,     OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,  fract_fpc_max, OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,  fract_forest,  OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,  sfract_can,    OUT)
!!$    alb_vis_can  => tile%mem%veg  %alb_vis_can%ptr     (ics:ice,    iblk)
!!$    alb_nir_can  => tile%mem%veg  %alb_nir_can%ptr     (ics:ice,    iblk)
!!$    alb_vis_soil => tile%mem%soil %alb_vis_soil%ptr    (ics:ice,    iblk)
!!$    alb_nir_soil => tile%mem%soil %alb_nir_soil%ptr    (ics:ice,    iblk)
    dsl$jsb Get_var2D_onChunk( srf_mem,  t_srf,         OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,  alb_vis_srf,   OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,  alb_nir_srf,   OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,  albedo_srf,    OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem, sfract_soil,   OUT)
    lai_mon      => veg_mem  %lai_mon (ics:ice,    iblk, :)

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
