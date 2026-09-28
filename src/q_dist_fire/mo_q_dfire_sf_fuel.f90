!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> Precompute fuel at PFT level
!>
!> Precalc PFT fraction for tree/grass (sum at aggregation)
!> (Oberhagemann only) Calc live grass moisture
!> (Oberhagemann only) Calc cured live grass
!> Convert to desired units required by Spitfire
!> Calc dead and live fuels...
!> Precalc fuel bulk density
!> Move fuel(s) into process memory
!>
MODULE mo_q_dfire_sf_fuel
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message, finish
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_POOL_ID, VEG_BGCM_POOL_ID

  USE mo_veg_constants,           ONLY: ITREE, IGRASS
  USE mo_q_dfire_sf_constants,    ONLY: df_sf_spq_layer_top_sl

  USE mo_q_dfire_sf_constants,    ONLY: df_sf_fc_1hr_idx, df_sf_fc_10hr_idx, df_sf_fc_100hr_idx, df_sf_fc_1000hr_idx
  USE mo_q_dfire_sf_constants,    ONLY: ratio_as_dead_fuel !< variable sb pool ratios to dead fuel

  USE mo_lnd_bgcm_idx,            ONLY: ixC
  USE mo_q_dfire_sf_fuel_class,   ONLY: t_sf_fuel_class
  USE mo_q_dfire_sf_util,         ONLY: calc_ros_fbd
  USE mo_q_dfire_sf_config_class, ONLY: Q_DFIRE_ROS_OBERHAGEMANN, Q_DFIRE_ROS_THONICKE

  !< Live fraction allocation for dead and live fuels
  USE mo_q_dfire_sf_constants,    ONLY: df_sf_live_frac_tree, &
                                        df_sf_live_frac_grass
  USE mo_q_dfire_sf_constants,    ONLY: dead_fuel_frac => df_sf_dead_frac      !< Fuel fraction to allocate from sb_pool

  USE mo_q_dfire_sf_adaptor_class,      ONLY: sf_convert_units
  USE mo_q_dfire_sf_constants,    ONLY: df_sf_units_k_to_c
  USE mo_q_dfire_sf_fmoist_class, ONLY: calc_mflg_ober
  USE mo_q_dfire_sf_classes,      ONLY: t_sf_lfmc


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_q_dfire_sf_fuel

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_fuel'

  !< Livegrass fuel moisture content (Oberhagemann et al. (2025) only)
  !< Group parameters from lctlib
  TYPE(t_sf_lfmc), SAVE, ALLOCATABLE :: lctlib_lfmc

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Compute spitfire for current time step
  !> It links Quincy and Spitfire
  !>
  SUBROUTINE update_q_dfire_sf_fuel(tile, options)
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: VEG_, HYDRO_, Q_DFIRE_FUEL_SF_, Q_DFIRE_SF_, A2L_, Q_PHENO_
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_grid,                ONLY: Get_vgrid, Get_grid
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_grid_class,          ONLY: t_jsb_grid
    USE mo_jsb_math_constants,      ONLY: one_day
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(Q_DFIRE_FUEL_SF_)
    dsl4jsb_Use_memory(Q_DFIRE_SF_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(Q_PHENO_)
    dsl4jsb_Use_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(Q_DFIRE_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_grid),       POINTER       :: grid
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lnd_bgcm_store), POINTER       :: bgcm_store !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc, ic
    INTEGER                               :: no_children
    REAL(wp)                              :: dtime      !< timestep

    INTEGER                               :: lctlib_growthform(options%nc)

    REAL(wp)                              :: tile_fract(options%nc)      !< [%]
    REAL(wp), POINTER                     :: grid_cell_area(:)           ! [m2]
    REAL(wp), POINTER                     :: tile_fract_ref(:)           ! [%]
    REAL(wp), ALLOCATABLE                 :: top_soil_depth_sl(:)        !< Top soil depth layer [m]
    REAL(wp)                              :: cured_lvg_fract(options%nc) !< Cured live grass fraction considered as dead fuel [%]

    REAL(wp), POINTER, CONTIGUOUS         :: live_fuel_frac(:,:)         !< live fuel classes allocation from veg_pool (#fuel classes,#live pools)
    REAL(wp), ALLOCATABLE                 :: sb_pool_slim(:,:)           !< debug only  (#pools, #nc)
    REAL(wp), ALLOCATABLE                 :: dead_fuel_fcl(:,:,:)        !< (#nc,#dead fuel pools,#fuel classes)
    REAL(wp), ALLOCATABLE                 :: live_grass_fuel_fcl(:,:,:)  !< (#nc,#live fuel pools,#fuel classes)
    REAL(wp), ALLOCATABLE                 :: cured_live_grass_fuel_fcl(:,:,:) !< (#nc,#dead fuel pools,#fuel classes)
    REAL(wp), ALLOCATABLE                 :: live_tree_fuel_fcl(:,:,:)   !< (#nc,#live fuel pools,#fuel classes)
    REAL(wp), ALLOCATABLE                 :: dead_fuel_spread_fcl(:,:)   !< (#nc,#dead fuel pools-1,10,100-,#fuel classes)
    REAL(wp), ALLOCATABLE                 :: live_fuel_fcl(:,:,:)        !< Live fuel (#nc,#live fuel pools, #fuel classes) [g DM m-2]
    REAL(wp), ALLOCATABLE                 :: veg_pool_in_gdm(:,:)        !< Live fuel (#nc,#live fuel pools) [g DM m-2]
    REAL(wp), ALLOCATABLE                 :: sb_pool_in_gdm(:,:)         !< Dead fuel (#nc,#dead fuel pools) [g DM m-2]
    TYPE(t_sf_fuel_class)                 :: fuel_cls                    !< Fuel [g DM m-2]

    CHARACTER(len=*), PARAMETER           :: routine = TRIM(modname)//':update_fdire_sf_fuel'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_DFIRE_SF_)
    dsl4jsb_Def_memory(Q_DFIRE_FUEL_SF_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(Q_PHENO_)
    dsl4jsb_Def_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Local vars here
    ! ----------------------------------------------------------------------------------------------------- !
    !< A2L_
    dsl4jsb_Real3D_onChunk      :: soil_depth_sl         ! Soil depth until bedrock [m]
    dsl4jsb_Real2D_onChunk      :: height                ! in [m]
    dsl4jsb_Real2D_onChunk      :: growing_season        ! in [%]

    !< Q_DFIRE_SF_
    dsl4jsb_Real3D_onChunk      :: dead_fuel_1hr_dfp     ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: dead_fuel_10hr_dfp    ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: dead_fuel_100hr_dfp   ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: dead_fuel_1000hr_dfp  ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: cured_live_grass_fuel_1hr_dfp  ! out [g DM m-2]

    dsl4jsb_Real3D_onChunk      :: live_grass_fuel_1hr_lfp    ! out [g DM m-2]

    dsl4jsb_Real3D_onChunk      :: live_tree_fuel_1hr_lfp     ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: live_tree_fuel_10hr_lfp    ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: live_tree_fuel_100hr_lfp   ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: live_tree_fuel_1000hr_lfp  ! out [g DM m-2]

    dsl4jsb_Real2D_onChunk      :: forest_fract       ! out [%]
    dsl4jsb_Real2D_onChunk      :: grass_fract        ! out [%]

    dsl4jsb_Real2D_onChunk      :: fuel_bulk_density    ! out [kg m-3]
    dsl4jsb_Real2D_onChunk      :: moisture_extinction  ! out [0 to 1]
    dsl4jsb_Real2D_onChunk      :: fwd_wind_factor      ! out [0 to 1]
    dsl4jsb_Real2D_onChunk      :: lvg_moisture         ! out [%], note: allowed up to 250%
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_DFIRE_FUEL_SF_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (lctlib%PastureFlag) RETURN !< do not run this routine at tiles like "pasture pft"
    IF (lctlib%CropFlag) RETURN !< do not run this routine at tiles like "crops pft"
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_DFIRE_SF_)
    dsl4jsb_Get_memory(Q_DFIRE_FUEL_SF_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(Q_PHENO_)
    dsl4jsb_Get_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)      ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, height)               ! in
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season)   ! in [%]

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, dead_fuel_1hr_dfp)     ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, dead_fuel_10hr_dfp)    ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, dead_fuel_100hr_dfp)   ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, dead_fuel_1000hr_dfp)  ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, cured_live_grass_fuel_1hr_dfp) ! out [g DM m-2]

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, live_grass_fuel_1hr_lfp)    ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, live_tree_fuel_1hr_lfp)     ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, live_tree_fuel_10hr_lfp)    ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, live_tree_fuel_100hr_lfp)   ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, live_tree_fuel_1000hr_lfp)  ! out [g DM m-2]

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_SF_, forest_fract)          ! out [%]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_SF_, grass_fract)           ! out [%]

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_SF_, fuel_bulk_density)     ! out [kg m-3]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_SF_, moisture_extinction)   ! out [0 to 1]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_SF_, fwd_wind_factor)       ! out [0 to 1]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_SF_, lvg_moisture)          ! out [%]

    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    ! ----------------------------------------------------------------------------------------------------- !
#ifndef __QUINCY_STANDALONE__
    CALL tile%Get_fraction(ics, ice, iblk, tile_fract)
    grid_cell_area => grid % area(ics:ice, iblk)
    tile_fract_ref => tile%fract(ics:ice, iblk)
#endif

    forest_fract=0._wp
    grass_fract=0._wp
    fuel_bulk_density = 0._wp
    top_soil_depth_sl = soil_depth_sl(:,df_sf_spq_layer_top_sl)
    lctlib_growthform(:) = lctlib%growthform

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !

    ! Aggregate fractions by type (Forest, Grass)
    ! todo: move to init (static setup)
#ifdef __QUINCY_STANDALONE__
    IF (lctlib%ForestFlag) THEN
      forest_fract=1._wp
    ELSE IF (lctlib%GrassFlag) THEN
      grass_fract = 1._wp
    !ELSE ! ignore
    END IF
#else
    IF (lctlib%ForestFlag) THEN
      forest_fract=tile_fract
    ELSE IF (lctlib%GrassFlag) THEN
      grass_fract=tile_fract
    !ELSE ! ignore
    END IF
#endif

    !< Oberhagemann related, required to calculate livegrass moisture at PFT level due to PFT parameters
    lvg_moisture = 0._wp
    cured_lvg_fract = 1._wp !< Thonicke's considers all live grass as dead fuel.

    IF (lctlib%growthform == igrass) THEN
      IF (dsl4jsb_Config(Q_DFIRE_SF_)%ros_model_id == Q_DFIRE_ROS_OBERHAGEMANN) THEN
        IF (.NOT. allocated(lctlib_lfmc)) lctlib_lfmc = t_sf_lfmc(lctlib)
        lvg_moisture = calc_mflg_ober(nc, lctlib_lfmc, growing_season)

        !< Oberhagemann et al. (2025), Section "Improvements to the treatment of live grass moisture in SPITFIRE", page 2031-2032
        cured_lvg_fract = max(0._wp, min(1._wp, (lvg_moisture - 0.3_wp) / (1.2_wp - 0.3_wp)))
      END IF
    END IF

    !< PFT related, todo: no need to redefine at every time step
    IF (lctlib%growthform == itree) THEN
      live_fuel_frac => df_sf_live_frac_tree
    ELSE IF (lctlib%growthform == igrass) THEN
      live_fuel_frac => df_sf_live_frac_grass
    END IF

    sb_pool_slim = sb_pool_mt(:,ixC,:,df_sf_spq_layer_top_sl)
    ALLOCATE(sb_pool_in_gdm, mold=sb_pool_slim) ! unsure why is needed, otherwise NAG compiler raises an error. function should be allocated dead_fuel
    sb_pool_in_gdm = sf_convert_units("molC/m3", "gDM/m2", sb_pool_slim, aux = top_soil_depth_sl, nc_dim_opt=2)

    ALLOCATE(veg_pool_in_gdm, mold=veg_pool_mt(:,ixC,:)) ! unsure why is needed, NAG compiler raises an error. sf_convert_units should already allocate live_fuel
    veg_pool_in_gdm = sf_convert_units("molC/m2", "gDM/m2", veg_pool_mt(:,ixC,:))

    fuel_cls = t_sf_fuel_class(nc, lctlib_growthform, &
                            dead_fuel_frac, live_fuel_frac, &
                            sb_pool_in_gdm, veg_pool_in_gdm, &
                            ratio_as_dead_fuel, cured_lvg_fract)

    dead_fuel_fcl = fuel_cls % Calc_dead_fuel() * dsl4jsb_Config(Q_DFIRE_SF_)%dead_fuel_factor

    !> Calculate the weighted average of the Fuel bulk density for the current PFT.
    !> Using the Pfeiffer et al. (2013) paper, equation A17. Remaining equations should be implemented
    !> when accounting for the other pfts.
    !> Thonicke's et al. (2010) is not used as the equation does not seems correct when accounting for 10 and 100 hours
    !> fuel bulk densities (slightly more).
    !> Aggregation average's the PFT values

    dead_fuel_spread_fcl = SUM(dead_fuel_fcl(:,:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx), DIM=2) !< 2 = fuel pools dimension index
    DO ic = 1, nc
      IF (SUM(dead_fuel_spread_fcl(ic,:)) > 0._wp) THEN
        fuel_bulk_density(ic) = calc_ros_fbd(3, lctlib%q_fuel_bulk_density, dead_fuel_spread_fcl(ic,:)) ! 3=fuel classes (1,10 and 100 hours)
      ENDIF
    END DO

    dead_fuel_1hr_dfp = dead_fuel_fcl(:,:,df_sf_fc_1hr_idx)
    dead_fuel_10hr_dfp = dead_fuel_fcl(:,:,df_sf_fc_10hr_idx)
    dead_fuel_100hr_dfp = dead_fuel_fcl(:,:,df_sf_fc_100hr_idx)
    dead_fuel_1000hr_dfp = dead_fuel_fcl(:,:,df_sf_fc_1000hr_idx)

    IF (lctlib%growthform == igrass) THEN
      ! Live grass to add to live fuel (not applies in Thonicke (so 0%), (1-cured)% in oberhagemann)
      live_grass_fuel_fcl = fuel_cls % Calc_live_grass_fuel() * dsl4jsb_Config(Q_DFIRE_SF_)%live_fuel_factor
      live_grass_fuel_1hr_lfp = live_grass_fuel_fcl(:,:,df_sf_fc_1hr_idx)
      ! Live grass to include as dead fuel (not applies to Thonicke (so 100%), x? in oberhagemann)
      cured_live_grass_fuel_fcl = fuel_cls % Calc_cured_live_grass_fuel() * dsl4jsb_Config(Q_DFIRE_SF_)%dead_fuel_factor
      cured_live_grass_fuel_1hr_dfp = cured_live_grass_fuel_fcl(:,:,df_sf_fc_1hr_idx)
    ELSE IF (lctlib%growthform == itree) THEN
      live_tree_fuel_fcl = fuel_cls % Calc_live_tree_fuel() * dsl4jsb_Config(Q_DFIRE_SF_)%live_fuel_factor
      live_tree_fuel_1hr_lfp = live_tree_fuel_fcl(:,:,df_sf_fc_1hr_idx)
      live_tree_fuel_10hr_lfp = live_tree_fuel_fcl(:,:,df_sf_fc_10hr_idx)
      live_tree_fuel_100hr_lfp = live_tree_fuel_fcl(:,:,df_sf_fc_100hr_idx)
      live_tree_fuel_1000hr_lfp = live_tree_fuel_fcl(:,:,df_sf_fc_1000hr_idx)
    END IF

    NULLIFY(live_fuel_frac)
  END SUBROUTINE update_q_dfire_sf_fuel


#endif
END MODULE mo_q_dfire_sf_fuel
