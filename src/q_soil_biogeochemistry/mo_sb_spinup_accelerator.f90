!> Monthly forcing bookkeeping and replay for QSOIL
!>
!> QUINCY soil biogeochemistry
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: https://doi.org/10.17871/quincy-model-2019
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> A coupled QLAND run accumulates monthly forcing sums in restart memory.
!> A subsequent QSOIL run replays the corresponding monthly means without
!> external forcing files or full-domain temporary caches.
MODULE mo_sb_spinup_accelerator
#ifndef __NO_JSBACH__

  USE mo_kind,            ONLY: lp
  USE mo_exception,       ONLY: finish, message
  USE mo_jsb_class,       ONLY: Get_model
  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_tile_class,  ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class,  ONLY: t_jsb_task_options
  USE mo_sb_config_class, ONLY: sb_spinup_accelerator_n_forcing_bins

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: collect_sb_spinup_accelerator_forcing
  PUBLIC :: update_read_soil_forcing

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_spinup_accelerator'
  LOGICAL, SAVE :: non_daily_timestep_warning_printed = .FALSE.

CONTAINS

  ! ===================================================================================================== !
  !> Accumulate one coupled-model timestep in its calendar-month forcing bin.
  SUBROUTINE collect_sb_spinup_accelerator_forcing(tile, options)

    USE mo_jsb_process_class,       ONLY: SB_, VEG_, HYDRO_, SSE_, A2L_
    USE mo_jsb_time,                ONLY: get_month
    USE mo_quincy_model_config,     ONLY: QLAND
    USE mo_lnd_bgcm_idx,            ONLY: ix_fine_root, ixC
    USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID, &
      &                                   VEG_BGCM_SEED_BED_LITTERFALL_ID
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(A2L_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt1L2D :: seed_bed_litterfall_mt

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store
    CHARACTER(len=*), PARAMETER :: routine = modname//':collect_sb_spinup_accelerator_forcing'
    INTEGER :: ics, ice, iblk
    INTEGER :: forcing_bin, i_source, i_element, i_soil
    INTEGER :: litter_entry, profile_entry, n_veg_sources, n_elements, nsoil
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_forcing_time
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_litterfall_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_fine_root_c_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_f_n_demand_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_f_p_demand_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_nhx_deposition_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_noy_deposition_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_nhx_n15_deposition_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_noy_n15_deposition_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_p_deposition_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_growth_p_limit_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_root_fraction_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_wtr_down_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_wtr_lat_loss_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_t_soil_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_wtr_soil_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_ice_soil_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_wtr_pot_sum
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    dsl4jsb_Real3D_onChunk :: frac_wtr_transp_down_sl
    dsl4jsb_Real3D_onChunk :: frac_w_lat_loss_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl
    dsl4jsb_Real3D_onChunk :: ice_soil_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_sl
    dsl4jsb_Real3D_onChunk :: t_soil_sl
    dsl4jsb_Real2D_onChunk :: f_n_demand
    dsl4jsb_Real2D_onChunk :: f_p_demand
    dsl4jsb_Real2D_onChunk :: growth_p_limit_based_on_n_mavg
    dsl4jsb_Real2D_onChunk :: nhx_deposition
    dsl4jsb_Real2D_onChunk :: noy_deposition
    dsl4jsb_Real2D_onChunk :: nhx_n15_deposition
    dsl4jsb_Real2D_onChunk :: noy_n15_deposition
    dsl4jsb_Real2D_onChunk :: p_deposition

    model => Get_model(tile%owner_model_id)
    IF (model%config%qmodel_id /= QLAND) RETURN
    IF (tile%lcts(1)%lib_id == 0) RETURN

    dsl4jsb_Get_config(SB_)
    IF (.NOT. dsl4jsb_Config(SB_)%flag_sb_model_spinup_accelerator) RETURN

    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    forcing_bin = get_month(options%current_datetime)
    IF (forcing_bin < 1 .OR. forcing_bin > sb_spinup_accelerator_n_forcing_bins) THEN
      CALL finish(TRIM(routine), 'Invalid calendar month for monthly SB forcing bookkeeping.')
    END IF

    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_forcing_time)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_litterfall_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_fine_root_c_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_f_n_demand_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_f_p_demand_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_nhx_deposition_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_noy_deposition_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_nhx_n15_deposition_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_noy_n15_deposition_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_p_deposition_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_growth_p_limit_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_root_fraction_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_wtr_down_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_wtr_lat_loss_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_t_soil_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_wtr_soil_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_ice_soil_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_wtr_pot_sum)
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, frac_wtr_transp_down_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, frac_w_lat_loss_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pot_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)
    dsl4jsb_Get_var2D_onChunk(VEG_, f_n_demand)
    dsl4jsb_Get_var2D_onChunk(VEG_, f_p_demand)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_p_limit_based_on_n_mavg)
    dsl4jsb_Get_var2D_onChunk(A2L_, nhx_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, noy_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, nhx_n15_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, noy_n15_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, p_deposition)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_LITTERFALL_ID, seed_bed_litterfall_mt)

    sb_spinup_accelerator_forcing_time(:, forcing_bin) = &
      & sb_spinup_accelerator_forcing_time(:, forcing_bin) + options%dtime
    sb_spinup_accelerator_fine_root_c_sum(:, forcing_bin) = &
      & sb_spinup_accelerator_fine_root_c_sum(:, forcing_bin) + &
      & veg_pool_mt(ix_fine_root, ixC, :) * options%dtime
    sb_spinup_accelerator_f_n_demand_sum(:, forcing_bin) = &
      & sb_spinup_accelerator_f_n_demand_sum(:, forcing_bin) + f_n_demand * options%dtime
    sb_spinup_accelerator_f_p_demand_sum(:, forcing_bin) = &
      & sb_spinup_accelerator_f_p_demand_sum(:, forcing_bin) + f_p_demand * options%dtime
    sb_spinup_accelerator_nhx_deposition_sum(:, forcing_bin) = &
      & sb_spinup_accelerator_nhx_deposition_sum(:, forcing_bin) + nhx_deposition * options%dtime
    sb_spinup_accelerator_noy_deposition_sum(:, forcing_bin) = &
      & sb_spinup_accelerator_noy_deposition_sum(:, forcing_bin) + noy_deposition * options%dtime
    sb_spinup_accelerator_nhx_n15_deposition_sum(:, forcing_bin) = &
      & sb_spinup_accelerator_nhx_n15_deposition_sum(:, forcing_bin) + nhx_n15_deposition * options%dtime
    sb_spinup_accelerator_noy_n15_deposition_sum(:, forcing_bin) = &
      & sb_spinup_accelerator_noy_n15_deposition_sum(:, forcing_bin) + noy_n15_deposition * options%dtime
    sb_spinup_accelerator_p_deposition_sum(:, forcing_bin) = &
      & sb_spinup_accelerator_p_deposition_sum(:, forcing_bin) + p_deposition * options%dtime
    sb_spinup_accelerator_growth_p_limit_sum(:, forcing_bin) = &
      & sb_spinup_accelerator_growth_p_limit_sum(:, forcing_bin) + &
      & growth_p_limit_based_on_n_mavg * options%dtime

    n_veg_sources = SIZE(veg_litterfall_mt, 1)
    n_elements = SIZE(veg_litterfall_mt, 2)
    IF (SIZE(seed_bed_litterfall_mt, 1) /= n_elements .OR. &
        & sb_spinup_accelerator_n_forcing_bins * (n_veg_sources + 1) * n_elements > &
        & SIZE(sb_spinup_accelerator_litterfall_sum, 2)) THEN
      CALL finish(TRIM(routine), 'Monthly SB litterfall forcing dimensions differ.')
    END IF

    DO i_element = 1, n_elements
      DO i_source = 1, n_veg_sources
        litter_entry = forcing_bin + sb_spinup_accelerator_n_forcing_bins * &
          & ((i_source - 1) + (n_veg_sources + 1) * (i_element - 1))
        sb_spinup_accelerator_litterfall_sum(:, litter_entry) = &
          & sb_spinup_accelerator_litterfall_sum(:, litter_entry) + &
          & veg_litterfall_mt(i_source, i_element, :)
      END DO
      litter_entry = forcing_bin + sb_spinup_accelerator_n_forcing_bins * &
        & (n_veg_sources + (n_veg_sources + 1) * (i_element - 1))
      sb_spinup_accelerator_litterfall_sum(:, litter_entry) = &
        & sb_spinup_accelerator_litterfall_sum(:, litter_entry) + seed_bed_litterfall_mt(i_element, :)
    END DO

    nsoil = SIZE(root_fraction_sl, 2)
    IF (SIZE(frac_wtr_transp_down_sl, 2) /= nsoil .OR. SIZE(frac_w_lat_loss_sl, 2) /= nsoil .OR. &
        & SIZE(wtr_soil_sl, 2) /= nsoil .OR. SIZE(ice_soil_sl, 2) /= nsoil .OR. &
        & SIZE(wtr_soil_pot_sl, 2) /= nsoil .OR. SIZE(t_soil_sl, 2) /= nsoil .OR. &
        & sb_spinup_accelerator_n_forcing_bins * nsoil > SIZE(sb_spinup_accelerator_root_fraction_sum, 2)) THEN
      CALL finish(TRIM(routine), 'Monthly SB soil-profile forcing dimensions differ.')
    END IF

    DO i_soil = 1, nsoil
      profile_entry = forcing_bin + sb_spinup_accelerator_n_forcing_bins * (i_soil - 1)
      sb_spinup_accelerator_root_fraction_sum(:, profile_entry) = &
        & sb_spinup_accelerator_root_fraction_sum(:, profile_entry) + root_fraction_sl(:, i_soil) * options%dtime
      sb_spinup_accelerator_wtr_down_sum(:, profile_entry) = &
        & sb_spinup_accelerator_wtr_down_sum(:, profile_entry) + &
        & frac_wtr_transp_down_sl(:, i_soil) * options%dtime
      sb_spinup_accelerator_wtr_lat_loss_sum(:, profile_entry) = &
        & sb_spinup_accelerator_wtr_lat_loss_sum(:, profile_entry) + frac_w_lat_loss_sl(:, i_soil) * options%dtime
      sb_spinup_accelerator_t_soil_sum(:, profile_entry) = &
        & sb_spinup_accelerator_t_soil_sum(:, profile_entry) + t_soil_sl(:, i_soil) * options%dtime
      sb_spinup_accelerator_wtr_soil_sum(:, profile_entry) = &
        & sb_spinup_accelerator_wtr_soil_sum(:, profile_entry) + wtr_soil_sl(:, i_soil) * options%dtime
      sb_spinup_accelerator_ice_soil_sum(:, profile_entry) = &
        & sb_spinup_accelerator_ice_soil_sum(:, profile_entry) + ice_soil_sl(:, i_soil) * options%dtime
      sb_spinup_accelerator_wtr_pot_sum(:, profile_entry) = &
        & sb_spinup_accelerator_wtr_pot_sum(:, profile_entry) + wtr_soil_pot_sl(:, i_soil) * options%dtime
    END DO

  END SUBROUTINE collect_sb_spinup_accelerator_forcing

  ! ===================================================================================================== !
  !> Replay the current calendar month's forcing means for one QSOIL timestep.
  SUBROUTINE update_read_soil_forcing(tile, options)

    USE mo_jsb_process_class,       ONLY: SB_, VEG_, HYDRO_, SSE_, A2L_
    USE mo_jsb_time,                ONLY: get_month
    USE mo_quincy_model_config,     ONLY: QSOIL
    USE mo_lnd_bgcm_idx,            ONLY: ix_fine_root, ixC
    USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID, &
      &                                   VEG_BGCM_SEED_BED_LITTERFALL_ID
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(A2L_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt1L2D :: seed_bed_litterfall_mt

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_read_soil_forcing'
    INTEGER :: ics, ice, iblk
    INTEGER :: forcing_bin, i_source, i_element, i_soil
    INTEGER :: litter_entry, profile_entry, n_veg_sources, n_elements, nsoil
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_forcing_time
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_litterfall_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_fine_root_c_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_f_n_demand_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_f_p_demand_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_nhx_deposition_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_noy_deposition_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_nhx_n15_deposition_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_noy_n15_deposition_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_p_deposition_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_growth_p_limit_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_root_fraction_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_wtr_down_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_wtr_lat_loss_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_t_soil_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_wtr_soil_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_ice_soil_sum
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_wtr_pot_sum
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    dsl4jsb_Real3D_onChunk :: frac_wtr_transp_down_sl
    dsl4jsb_Real3D_onChunk :: frac_w_lat_loss_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl
    dsl4jsb_Real3D_onChunk :: ice_soil_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_sl
    dsl4jsb_Real3D_onChunk :: t_soil_sl
    dsl4jsb_Real2D_onChunk :: f_n_demand
    dsl4jsb_Real2D_onChunk :: f_p_demand
    dsl4jsb_Real2D_onChunk :: growth_p_limit_based_on_n_mavg
    dsl4jsb_Real2D_onChunk :: nhx_deposition
    dsl4jsb_Real2D_onChunk :: noy_deposition
    dsl4jsb_Real2D_onChunk :: nhx_n15_deposition
    dsl4jsb_Real2D_onChunk :: noy_n15_deposition
    dsl4jsb_Real2D_onChunk :: p_deposition

    model => Get_model(tile%owner_model_id)
    IF (model%config%qmodel_id /= QSOIL) RETURN
    IF (tile%lcts(1)%lib_id == 0) RETURN

    dsl4jsb_Get_config(SB_)
    IF (.NOT. dsl4jsb_Config(SB_)%flag_sb_model_spinup_accelerator) THEN
      CALL finish(TRIM(routine), 'QSOIL requires monthly SB forcing bookkeeping in its restart input.')
    END IF

    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    forcing_bin = get_month(options%current_datetime)
    IF (forcing_bin < 1 .OR. forcing_bin > sb_spinup_accelerator_n_forcing_bins) THEN
      CALL finish(TRIM(routine), 'Invalid calendar month for monthly QSOIL forcing replay.')
    END IF

    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_forcing_time)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_litterfall_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_fine_root_c_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_f_n_demand_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_f_p_demand_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_nhx_deposition_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_noy_deposition_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_nhx_n15_deposition_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_noy_n15_deposition_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_p_deposition_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_growth_p_limit_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_root_fraction_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_wtr_down_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_wtr_lat_loss_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_t_soil_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_wtr_soil_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_ice_soil_sum)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_wtr_pot_sum)
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, frac_wtr_transp_down_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, frac_w_lat_loss_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pot_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)
    dsl4jsb_Get_var2D_onChunk(VEG_, f_n_demand)
    dsl4jsb_Get_var2D_onChunk(VEG_, f_p_demand)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_p_limit_based_on_n_mavg)
    dsl4jsb_Get_var2D_onChunk(A2L_, nhx_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, noy_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, nhx_n15_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, noy_n15_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, p_deposition)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_LITTERFALL_ID, seed_bed_litterfall_mt)

    IF (ANY(sb_spinup_accelerator_forcing_time(:, forcing_bin) <= 0.0_lp)) THEN
      CALL finish(TRIM(routine), 'The QSOIL restart has no sampled forcing for the current calendar month.')
    END IF

    veg_pool_mt(ix_fine_root, ixC, :) = sb_spinup_accelerator_fine_root_c_sum(:, forcing_bin) / &
      & sb_spinup_accelerator_forcing_time(:, forcing_bin)
    f_n_demand = sb_spinup_accelerator_f_n_demand_sum(:, forcing_bin) / &
      & sb_spinup_accelerator_forcing_time(:, forcing_bin)
    f_p_demand = sb_spinup_accelerator_f_p_demand_sum(:, forcing_bin) / &
      & sb_spinup_accelerator_forcing_time(:, forcing_bin)
    nhx_deposition = sb_spinup_accelerator_nhx_deposition_sum(:, forcing_bin) / &
      & sb_spinup_accelerator_forcing_time(:, forcing_bin)
    noy_deposition = sb_spinup_accelerator_noy_deposition_sum(:, forcing_bin) / &
      & sb_spinup_accelerator_forcing_time(:, forcing_bin)
    nhx_n15_deposition = sb_spinup_accelerator_nhx_n15_deposition_sum(:, forcing_bin) / &
      & sb_spinup_accelerator_forcing_time(:, forcing_bin)
    noy_n15_deposition = sb_spinup_accelerator_noy_n15_deposition_sum(:, forcing_bin) / &
      & sb_spinup_accelerator_forcing_time(:, forcing_bin)
    p_deposition = sb_spinup_accelerator_p_deposition_sum(:, forcing_bin) / &
      & sb_spinup_accelerator_forcing_time(:, forcing_bin)
    growth_p_limit_based_on_n_mavg = sb_spinup_accelerator_growth_p_limit_sum(:, forcing_bin) / &
      & sb_spinup_accelerator_forcing_time(:, forcing_bin)

    n_veg_sources = SIZE(veg_litterfall_mt, 1)
    n_elements = SIZE(veg_litterfall_mt, 2)
    IF (SIZE(seed_bed_litterfall_mt, 1) /= n_elements .OR. &
        & sb_spinup_accelerator_n_forcing_bins * (n_veg_sources + 1) * n_elements > &
        & SIZE(sb_spinup_accelerator_litterfall_sum, 2)) THEN
      CALL finish(TRIM(routine), 'Monthly QSOIL litterfall forcing dimensions differ.')
    END IF

    DO i_element = 1, n_elements
      DO i_source = 1, n_veg_sources
        litter_entry = forcing_bin + sb_spinup_accelerator_n_forcing_bins * &
          & ((i_source - 1) + (n_veg_sources + 1) * (i_element - 1))
        veg_litterfall_mt(i_source, i_element, :) = &
          & sb_spinup_accelerator_litterfall_sum(:, litter_entry) / &
          & sb_spinup_accelerator_forcing_time(:, forcing_bin) * options%dtime
      END DO
      litter_entry = forcing_bin + sb_spinup_accelerator_n_forcing_bins * &
        & (n_veg_sources + (n_veg_sources + 1) * (i_element - 1))
      seed_bed_litterfall_mt(i_element, :) = sb_spinup_accelerator_litterfall_sum(:, litter_entry) / &
        & sb_spinup_accelerator_forcing_time(:, forcing_bin) * options%dtime
    END DO

    nsoil = SIZE(root_fraction_sl, 2)
    IF (SIZE(frac_wtr_transp_down_sl, 2) /= nsoil .OR. SIZE(frac_w_lat_loss_sl, 2) /= nsoil .OR. &
        & SIZE(wtr_soil_sl, 2) /= nsoil .OR. SIZE(ice_soil_sl, 2) /= nsoil .OR. &
        & SIZE(wtr_soil_pot_sl, 2) /= nsoil .OR. SIZE(t_soil_sl, 2) /= nsoil .OR. &
        & sb_spinup_accelerator_n_forcing_bins * nsoil > SIZE(sb_spinup_accelerator_root_fraction_sum, 2)) THEN
      CALL finish(TRIM(routine), 'Monthly QSOIL soil-profile forcing dimensions differ.')
    END IF

    DO i_soil = 1, nsoil
      profile_entry = forcing_bin + sb_spinup_accelerator_n_forcing_bins * (i_soil - 1)
      root_fraction_sl(:, i_soil) = sb_spinup_accelerator_root_fraction_sum(:, profile_entry) / &
        & sb_spinup_accelerator_forcing_time(:, forcing_bin)
      frac_wtr_transp_down_sl(:, i_soil) = sb_spinup_accelerator_wtr_down_sum(:, profile_entry) / &
        & sb_spinup_accelerator_forcing_time(:, forcing_bin)
      frac_w_lat_loss_sl(:, i_soil) = sb_spinup_accelerator_wtr_lat_loss_sum(:, profile_entry) / &
        & sb_spinup_accelerator_forcing_time(:, forcing_bin)
      t_soil_sl(:, i_soil) = sb_spinup_accelerator_t_soil_sum(:, profile_entry) / &
        & sb_spinup_accelerator_forcing_time(:, forcing_bin)
      wtr_soil_sl(:, i_soil) = sb_spinup_accelerator_wtr_soil_sum(:, profile_entry) / &
        & sb_spinup_accelerator_forcing_time(:, forcing_bin)
      ice_soil_sl(:, i_soil) = sb_spinup_accelerator_ice_soil_sum(:, profile_entry) / &
        & sb_spinup_accelerator_forcing_time(:, forcing_bin)
      wtr_soil_pot_sl(:, i_soil) = sb_spinup_accelerator_wtr_pot_sum(:, profile_entry) / &
        & sb_spinup_accelerator_forcing_time(:, forcing_bin)
    END DO

    IF (NINT(options%dtime) /= 86400 .AND. .NOT. non_daily_timestep_warning_printed) THEN
      CALL message(TRIM(routine), 'Monthly QSOIL forcing is designed for a daily (P1D) model timestep.')
      non_daily_timestep_warning_printed = .TRUE.
    END IF

  END SUBROUTINE update_read_soil_forcing

#endif
END MODULE mo_sb_spinup_accelerator
