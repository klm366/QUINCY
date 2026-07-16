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
!>#### task to update fire impacts into the model
!>
MODULE mo_q_dfire_sf_update
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_lnd_bgcm_store,        ONLY: t_lnd_bgcm_store
  USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
  USE mo_jsb_grid,              ONLY: Get_grid
  USE mo_jsb_math_constants,    ONLY: eps8

  USE mo_veg_constants,         ONLY: itree
  USE mo_q_dfire_sf_constants,  ONLY: df_sf_spq_layer_top_sl
  USE mo_q_dfire_sf_util,       ONLY: dfire_sf_sb_update_pools, dfire_sf_veg_update_pools
  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store_class,  ONLY: SB_BGCM_POOL_ID, VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID, &
                                      SB_BGCM_FIRE_EMISSIONS_ID, SB_BGCM_FIRE_LITTER_ASH_ID, &
                                      VEG_BGCM_FIRE_EMISSIONS_ID, VEG_BGCM_FIRE_LITTERFALL_ID

  USE mo_q_dfire_sf_constants,  ONLY: df_sf_spq_layer_top_sl

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_q_dfire_sf_update

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_update'

CONTAINS

  ! ======================================================================================================= !
  !>
  !>
  !>
  SUBROUTINE update_q_dfire_sf_update(tile, options)
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_grid_class,          ONLY: t_jsb_grid
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: VEG_, HYDRO_, Q_DFIRE_FUEL_SF_, Q_DFIRE_SF_
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_math_constants,      ONLY: one_day

    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(Q_DFIRE_FUEL_SF_)
    dsl4jsb_Use_memory(Q_DFIRE_SF_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(VEG_)
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
    INTEGER                               :: iblk, ics, ice, nc, ic, ipool, iel, ip
    INTEGER                               :: no_children
    REAL(wp)                              :: dtime      !< timestep

    ! add here local variables
    REAL(wp), POINTER                     :: area(:)
    REAL(wp), DIMENSION(options%nc)       :: cover_fraction             !< current tile area
    REAL(wp), DIMENSION(options%nc)       :: area_burned_fract_pft !< Fraction of the area burned for current PFT [-]
    REAL(wp), DIMENSION(options%nc)       :: fract_area_burned_rel_to_tile  !< area burned fraction relative to the tile area

    CHARACTER(len=*), PARAMETER           :: routine = TRIM(modname)//':update_dfire_sf_fuel'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real2D_onChunk                :: area_burned_fract         ! in [-]
    dsl4jsb_Real2D_onChunk                :: delta_dens_ind            ! out [-]
    dsl4jsb_Real3D_onChunk                :: delta_root_fraction_sl    ! out [-]
    dsl4jsb_Real2D_onChunk                :: dens_ind                  ! in [-]
    dsl4jsb_Real3D_onChunk                :: root_fraction_sl          ! in [-]
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_fire_emissions_mt
    dsl4jsb_Def_mt2L3D :: sb_fire_litter_ash_mt
    dsl4jsb_Def_mt2L2D :: veg_fire_emissions_mt
    dsl4jsb_Def_mt2L2D :: veg_fire_litterfall_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_DFIRE_SF_)
    dsl4jsb_Def_memory(Q_DFIRE_FUEL_SF_)
    dsl4jsb_Def_memory(Q_DFIRE_SF_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Local vars here
    ! ----------------------------------------------------------------------------------------------------- !

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
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (lctlib%PastureFlag) RETURN !< do not run this routine at tiles like "pasture pft"
    IF (lctlib%CropFlag) RETURN !< do not run this routine at tiles like "crops pft"
    ! ----------------------------------------------------------------------------------------------------- !
    grid => Get_grid(model%grid_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_DFIRE_SF_)
    dsl4jsb_Get_memory(Q_DFIRE_FUEL_SF_)
    dsl4jsb_Get_memory(Q_DFIRE_SF_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, area_burned_fract)      ! in [-]
    dsl4jsb_Get_var2D_onChunk(VEG_, delta_dens_ind)                ! out [-]
    dsl4jsb_Get_var3D_onChunk(VEG_, delta_root_fraction_sl)        ! out[-]
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)                      ! in [-]
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)              ! in [-]
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_FIRE_EMISSIONS_ID, sb_fire_emissions_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_FIRE_LITTER_ASH_ID, sb_fire_litter_ash_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_FIRE_EMISSIONS_ID, veg_fire_emissions_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_FIRE_LITTERFALL_ID, veg_fire_litterfall_mt)
    ! ----------------------------------------------------------------------------------------------------- !

    ! init here
#ifdef __QUINCY_STANDALONE__
        ALLOCATE(area(nc))
        area = 1._wp
        cover_fraction = 1._wp
#else
        area => grid %area (ics:ice, iblk)
        CALL tile%Get_fraction(ics, ice, iblk, fract=cover_fraction)
#endif
    area_burned_fract_pft = area_burned_fract * cover_fraction
    delta_dens_ind = 0._wp

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !

    !< Enforce to update pools due to Spitfire impacts. Temporary solution until a permanent solution
    !< is found to properly deal.
    !< Spitfire is now updated here

    !< Update SB pools
    IF (dsl4jsb_Config(Q_DFIRE_SF_)%flag_has_dead_feedback) THEN
      CALL dfire_sf_sb_update_pools(sb_fire_emissions_mt(:,:,:,df_sf_spq_layer_top_sl), &
                                    sb_fire_litter_ash_mt(:,:,:,df_sf_spq_layer_top_sl), &
                                    sb_pool_mt(:,:,:,df_sf_spq_layer_top_sl))
    END IF

    !< Update VEG pools
    IF (dsl4jsb_Config(Q_DFIRE_SF_)%flag_has_live_feedback) THEN
      CALL dfire_sf_veg_update_pools(nc, model, &
                                    veg_fire_emissions_mt, veg_fire_litterfall_mt, &
                                    veg_pool_mt, veg_litterfall_mt)
    END IF

    !>
    !> Update delta_dens_ind to properly run plant_characteristics in mo_q_model_init.f90
    !> Calculate change in individual density from Area burned
    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> - Adapt the change of individual density to account for the harvested fraction (only for tiles with forest type)
    !>
    IF (lctlib%growthform == itree) THEN
      fract_area_burned_rel_to_tile = 0._wp

      !< Calculate area burned fraction relative to tile area
      DO ic = 1,nc
        IF (cover_fraction(ic) > eps8) THEN !< is cover fraction too small to consider?
          fract_area_burned_rel_to_tile(ic) = MIN(1.0_wp, area_burned_fract_pft(ic) / cover_fraction(ic))
        END IF
      END DO

      DO ic = 1,nc
        delta_dens_ind(ic) = - (fract_area_burned_rel_to_tile(ic) * dens_ind(ic))
      END DO
    END IF

  END SUBROUTINE update_q_dfire_sf_update

#endif
END MODULE mo_q_dfire_sf_update
