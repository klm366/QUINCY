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

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile
  USE mo_jsb_base,           ONLY: get_model
  USE mo_jsb_parallel_iface, ONLY: get_my_global_mpi_id
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish
  USE mo_util,               ONLY: one_of, int2string

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

    REAL(wp), POINTER ::                 &
      & lai                    (:),      &
      & swrad_srf_down         (:),      &
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

    IF (one_of('veg', tile%components) < 1) RETURN

!!$    p_pe = get_my_global_mpi_id()
!!$    IF (iblk == 1) CALL message(TRIM(routine), 'Starting update on PE '//int2string(p_pe), all_print=.TRUE.)
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of vegetation_1 ')

    model => get_model(tile%owner_model_id)

    nsoil = model%config%soil%nsoil_water

    ! Get reference to variables for current block
    !
    lai                => tile%mem%veg  %phys_state   %lai            (ics:ice,    iblk)
    swrad_srf_down     => tile%mem%srf  %energy_flux  %swrad_srf_down (ics:ice,    iblk)
    canopy_cond        => tile%mem%veg  %phys_state   %canopy_cond    (ics:ice,    iblk)
    soil_depth_l       => tile%mem%soil %phys_state   %soil_depth_l   (ics:ice, :, iblk)
    root_depth_l       => tile%mem%veg  %phys_state   %root_depth_l   (ics:ice, :, iblk)
    ws                 => tile%mem%soil %water_state  %ws             (ics:ice,    iblk)
    ws_l               => tile%mem%soil %water_state  %ws_l           (ics:ice, :, iblk)
    ws_fc_l            => tile%mem%soil %water_state  %ws_fc_l        (ics:ice, :, iblk)
    ws_pwp_l           => tile%mem%soil %water_state  %ws_pwp_l       (ics:ice, :, iblk)
    ws_root            => tile%mem%veg  %water_state  %ws_root        (ics:ice,    iblk)
    ws_fc_root         => tile%mem%veg  %water_state  %ws_fc_root     (ics:ice,    iblk)
    ws_pwp_root        => tile%mem%veg  %water_state  %ws_pwp_root    (ics:ice,    iblk)
    water_stress       => tile%mem%veg  %water_state  %water_stress   (ics:ice,    iblk)
    t_srf              => tile%mem%srf  %energy_state %t_srf          (ics:ice,    iblk)
    q_air              => tile%mem      %atm2land     %q_air          (ics:ice,    iblk)
    press_srf          => tile%mem      %atm2land     %press_srf      (ics:ice,    iblk)

    ! Compute (max) canopy (stomatal) conductance with no water stress
    SELECT CASE(TRIM(model%config%veg%scheme_stom_cond))
      CASE('echam5')
        canopy_cond(:) = canopy_conductance( lai(:), swrad_srf_down(:) / 2._wp )

!!$      CASE('bethy')
!!$        ! TBD: Call bethy to compute unstressed canopy conductance

      CASE default
        CALL finish( TRIM(routine), 'Canopy conductance scheme '//TRIM(model%config%veg%scheme_stom_cond)//&
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
    USE mo_jsb_time,    ONLY: wgt1, wgt2, nmw1,nmw2  ! TBD: this should be made more general

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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_vegetation_2'

    IF (one_of('veg', tile%components) < 1) RETURN

    p_pe = get_my_global_mpi_id()
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    t_srf        => tile%mem%srf  %energy_state %t_srf           (ics:ice,    iblk)
    lai          => tile%mem%veg  %phys_state   %lai             (ics:ice,    iblk)
    lai_mon      => tile%mem%veg  %phys_state   %lai_mon         (ics:ice,    iblk, :)
    fract_fpc    => tile%mem%veg  %phys_state   %fract_fpc       (ics:ice,    iblk)
    fract_fpc_max=> tile%mem%veg  %phys_state   %fract_fpc_max   (ics:ice,    iblk)
    fract_forest => tile%mem%veg  %phys_state   %fract_forest    (ics:ice,    iblk)
!!$    alb_vis_can  => tile%mem%veg  %phys_state   %alb_vis_can     (ics:ice,    iblk)
!!$    alb_nir_can  => tile%mem%veg  %phys_state   %alb_nir_can     (ics:ice,    iblk)
!!$    alb_vis_soil => tile%mem%soil %phys_state   %alb_vis_soil    (ics:ice,    iblk)
!!$    alb_nir_soil => tile%mem%soil %phys_state   %alb_nir_soil    (ics:ice,    iblk)
    alb_vis_srf  => tile%mem%srf  %phys_state   %alb_vis_srf     (ics:ice,    iblk)
    alb_nir_srf  => tile%mem%srf  %phys_state   %alb_nir_srf     (ics:ice,    iblk)
    albedo_srf   => tile%mem%srf  %phys_state   %albedo_srf      (ics:ice,    iblk)
    sfract_can   => tile%mem%veg  %water_state  %sfract_can      (ics:ice,    iblk)
    sfract_soil  => tile%mem%soil %water_state  %sfract_soil     (ics:ice,    iblk)

    ! Update leaf area index
    !
    ! At the moment by interpolating monthly climatological LAI
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
