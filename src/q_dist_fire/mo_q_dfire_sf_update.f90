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
!> Add description here
!>
MODULE mo_q_dfire_sf_update
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_lnd_bgcm_store,        ONLY: t_lnd_bgcm_store
  USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
  USE mo_jsb_grid,              ONLY: Get_grid

  USE mo_q_dfire_util, ONLY: dfire_update_pools
  USE mo_lnd_bgcm_idx

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
    USE mo_jsb_process_class,       ONLY: VEG_, HYDRO_, Q_DFIRE_FUEL_, Q_DFIRE_
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_math_constants,      ONLY: one_day

    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(Q_DFIRE_FUEL_)
    dsl4jsb_Use_memory(Q_DFIRE_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(Q_DFIRE_)
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


    CHARACTER(len=*), PARAMETER           :: routine = TRIM(modname)//':update_dfire_sf_fuel'
    ! ----------------------------------------------------------------------------------------------------- !
    !dsl4jsb_Real3D_onChunk                :: live_tree_fuel_1hr_lfp    ! in [g DM m-2]

    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_fire_emissions_mt
    dsl4jsb_Def_mt2L3D :: sb_fire_litter_ash_mt
    dsl4jsb_Def_mt2L2D :: veg_fire_emissions_mt
    dsl4jsb_Def_mt2L2D :: veg_fire_litterfall_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_DFIRE_)
    dsl4jsb_Def_memory(Q_DFIRE_FUEL_)
    dsl4jsb_Def_memory(Q_DFIRE_)
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
    IF (.NOT. tile%Is_process_calculated(Q_DFIRE_FUEL_)) RETURN
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
    dsl4jsb_Get_config(Q_DFIRE_)
    dsl4jsb_Get_memory(Q_DFIRE_FUEL_)
    dsl4jsb_Get_memory(Q_DFIRE_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    !dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1hr_lfp)    ! in [g DM m-2]
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


    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !

    !< Enforce to update pools due to Spitfire impacts. Temporary solution until a permanent solution
    !< is found to properly deal.
    !< Spitfire is now updated here

    !< Update pools
    CALL dfire_update_pools(sb_fire_emissions_mt(:,:,:,df_sf_spq_layer_top_sl), &
                            sb_fire_litter_ash_mt(:,:,:,df_sf_spq_layer_top_sl), &
                            sb_pool_mt(:,:,:,df_sf_spq_layer_top_sl))

    nr_of_veg_bgcm_comp = SIZE(veg_pool_mt,1)

    !> 3.0 add fire disturbance into litterfall
    DO ic = 1, nc
      DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
        IF (model%config%is_element_used(ielem)) THEN
          ix_elem = model%config%elements_index_map(ielem)    ! get element index in bgcm
          DO ix_comp = 1, nr_of_veg_bgcm_comp
            veg_litterfall_mt(ix_comp, ix_elem, ic) = veg_litterfall_mt(ix_comp, ix_elem, ic) + &
                                                      veg_fire_litterfall_mt(ix_comp, ix_elem, ic)
          END DO ! ix_comp = 1, n_veg_bgcm_comp
        END IF ! model%config%is_element_used(ielem)) THEN
     END DO ! ielem = FIRST_ELEM_ID, LAST_ELEM_ID
    END DO ! ic = 1, nc

    !> update veg fire emissions
    DO ic = 1, nc
      DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
        IF (model%config%is_element_used(ielem)) THEN
          ix_elem = model%config%elements_index_map(ielem)    ! get element index in bgcm
          DO ix_comp = 1, nr_of_veg_bgcm_comp
                veg_pool_mt(ix_comp,ix_elem,ic) = veg_pool_mt(ix_comp,ix_elem,ic) - &
                    veg_fire_emissions_mt(ix_comp,ix_elem,ic) ! dist_fire
          END DO ! ix_comp
        END IF ! is_element_used(ielem)
      END DO ! ielem
    END DO ! ic

  END SUBROUTINE update_q_dfire_sf_update

#endif
END MODULE mo_q_dfire_sf_update
