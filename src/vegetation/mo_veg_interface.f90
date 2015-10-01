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

  Use_processes VEG, SOIL, SRF, ATM2LAND
  Use_config(veg)
  Use_config(soil)

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  Use_memory(atm2land)
  Use_memory(soil)
  Use_memory(srf)
  Use_memory(veg)

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

    Def_config(veg_conf, veg)
    Def_config(soil_conf, soil)
    Def_memory(soil_mem, soil)
    Def_memory(srf_mem, srf)
    Def_memory(veg_mem, veg)
    Def_memory(atm2land_mem, atm2land)

    REAL(wp), POINTER ::                 &
      & lai                    (:),      &
      & swpar_srf_down         (:),      &
      & canopy_cond            (:),      &
      & soil_depth_l           (:,:),    &
      & root_depth_l           (:,:),    &
      & ws                     (:),      & !< Water content in soil column
      & ws_l                   (:,:),    & !< Water content in soil layers
      & ws_fc_l                (:,:),    & !< Water content of soil layers at field capacity
      & ws_pwp_l               (:,:),    & !< Water content of soil layers at permanent wilting point
      & ws_root                (:),      & !< Water content in root zone
      & ws_fc_root             (:),      & !< Water content in root zone at field capacity
      & ws_pwp_root            (:),      & !< Water content in root zone at permanent wilting point
      & water_stress           (:),      &
      & t_srf                  (:),      &
      & q_air                  (:),      &
      & press_srf              (:)

    INTEGER :: nsoil

    INTEGER :: p_pe

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_vegetation_1'

    IF (.NOT. tile%Is_process_active(VEG)) RETURN

!!$    p_pe = get_my_global_mpi_id()
!!$    IF (iblk == 1) CALL message(TRIM(routine), 'Starting update on PE '//int2string(p_pe), all_print=.TRUE.)
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of vegetation_1 ')

    model => get_model(tile%owner_model_id)

    Get_config(veg_conf, VEG)
    Get_config(soil_conf, SOIL)

    nsoil = soil_conf%nsoil_water

    ! Get reference to variables for current block
    !
    Get_memory(veg_mem, VEG)
    Get_memory(srf_mem, SRF)
    Get_memory(soil_mem, SOIL)
    Get_memory(atm2land_mem, ATM2LAND)

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

    REAL(wp), POINTER ::   &
      & t_srf         (:),   &
      & lai           (:),   &
      & lai_mon       (:,:), &
      & fract_fpc     (:),   &
      & fract_fpc_max (:),   &
      & fract_forest  (:),   &
!!$      & alb_vis_can  (:),    &
!!$      & alb_nir_can  (:),    &
!!$      & alb_vis_soil (:),    &
!!$      & alb_nir_soil (:),    &
      & alb_vis_srf   (:),   &
      & alb_nir_srf   (:),   &
      & albedo_srf    (:),   &
      & sfract_can    (:),   &
      & sfract_soil   (:)

    INTEGER :: p_pe

    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2

    Def_memory(veg_mem, veg)
    Def_memory(srf_mem, srf)
    Def_memory(soil_mem, soil)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_vegetation_2'

    IF (.NOT. tile%Is_process_active(VEG)) RETURN

    p_pe = get_my_global_mpi_id()
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    Get_memory(veg_mem, VEG)
    Get_memory(srf_mem, SRF)
    Get_memory(soil_mem, SOIL)

    t_srf        => srf_mem  %t_srf%ptr           (ics:ice,    iblk)
    lai          => veg_mem  %lai%ptr             (ics:ice,    iblk)
    lai_mon      => veg_mem  %lai_mon%ptr         (ics:ice,    iblk, :)
    fract_fpc    => veg_mem  %fract_fpc%ptr       (ics:ice,    iblk)
    fract_fpc_max=> veg_mem  %fract_fpc_max%ptr   (ics:ice,    iblk)
    fract_forest => veg_mem  %fract_forest%ptr    (ics:ice,    iblk)
!!$    alb_vis_can  => tile%mem%veg  %alb_vis_can%ptr     (ics:ice,    iblk)
!!$    alb_nir_can  => tile%mem%veg  %alb_nir_can%ptr     (ics:ice,    iblk)
!!$    alb_vis_soil => tile%mem%soil %alb_vis_soil%ptr    (ics:ice,    iblk)
!!$    alb_nir_soil => tile%mem%soil %alb_nir_soil%ptr    (ics:ice,    iblk)
    alb_vis_srf  => srf_mem  %alb_vis_srf%ptr     (ics:ice,    iblk)
    alb_nir_srf  => srf_mem  %alb_nir_srf%ptr     (ics:ice,    iblk)
    albedo_srf   => srf_mem  %albedo_srf%ptr      (ics:ice,    iblk)
    sfract_can   => veg_mem  %sfract_can%ptr      (ics:ice,    iblk)
    sfract_soil  => soil_mem %sfract_soil%ptr     (ics:ice,    iblk)

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
