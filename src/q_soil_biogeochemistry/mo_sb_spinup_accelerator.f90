!> Forward soil-biogeochemistry spin-up accelerator
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
!> Structural entry point for replaying bookkept soil forcing during spin-up.
!> The forward relaxation is shared by JSM and simple_1d and starts from the
!> current soil state. It has no dependency on an analytical steady-state solver.
!>
!> Memory contract for the implementation stage:
!> - forcing-bin weighted sums and elapsed-time weights are restart-safe SB memory;
!> - execution count and convergence history are restart-safe SB memory;
!> - live-forcing backups and replay work arrays are temporary local scratch state.
MODULE mo_sb_spinup_accelerator
#ifndef __NO_JSBACH__

  USE mo_kind,            ONLY: lp
  USE mo_exception,       ONLY: finish, message, message_text
  USE mo_jsb_class,       ONLY: Get_model
  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_tile_class,  ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class,  ONLY: t_jsb_task_options
  USE mo_sb_config_class, ONLY: sb_spinup_accelerator_max_forcing_bins

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: apply_sb_spinup_accelerator

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_spinup_accelerator'
  CHARACTER(len=*), PARAMETER :: sb_spinup_debug_filename = 'sb_spinup_accelerator_debug.txt'
  INTEGER, PARAMETER :: sb_spinup_debug_years(10) = [1, 2, 5, 10, 50, 100, 500, 1000, 5000, 10000]

  !> Event-local climatology reconstructed from restart-safe accumulated sums.
  TYPE :: t_sb_spinup_forcing_climatology
    INTEGER :: n_forcing_bins = 0
    INTEGER :: n_soil_layers = 0
    INTEGER :: n_litter_sources = 0
    INTEGER :: n_elements = 0
    REAL(lp), ALLOCATABLE :: forcing_time(:,:)     !< sampled time by cell and bin [s]
    REAL(lp), ALLOCATABLE :: litterfall_rate(:,:)  !< mean litterfall rate, flattened by bin/source/element [mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: fine_root_c(:,:)       !< mean fine-root carbon by cell and bin [mol C m-2]
    REAL(lp), ALLOCATABLE :: f_n_demand(:,:)        !< mean plant nitrogen-demand fraction by cell and bin [-]
    REAL(lp), ALLOCATABLE :: f_p_demand(:,:)        !< mean plant phosphorus-demand fraction by cell and bin [-]
    REAL(lp), ALLOCATABLE :: nhx_deposition(:,:)    !< mean reduced-N deposition by cell and bin [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: noy_deposition(:,:)    !< mean oxidized-N deposition by cell and bin [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: nhx_n15_deposition(:,:) !< mean reduced-15N deposition [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: noy_n15_deposition(:,:) !< mean oxidized-15N deposition [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: p_deposition(:,:)      !< mean phosphorus deposition by cell and bin [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: growth_p_limit(:,:)    !< mean plant P limitation based on N demand by cell and bin [-]
    REAL(lp), ALLOCATABLE :: root_fraction(:,:)    !< mean root fraction, flattened by bin/layer [-]
    REAL(lp), ALLOCATABLE :: wtr_down(:,:)         !< mean downward water-transfer fraction, flattened by bin/layer [-]
    REAL(lp), ALLOCATABLE :: wtr_lat_loss(:,:)     !< mean lateral water-loss fraction, flattened by bin/layer [-]
    REAL(lp), ALLOCATABLE :: t_soil(:,:)           !< mean soil temperature, flattened by bin/layer [K]
    REAL(lp), ALLOCATABLE :: wtr_soil(:,:)         !< mean liquid soil water, flattened by bin/layer [m]
    REAL(lp), ALLOCATABLE :: ice_soil(:,:)         !< mean soil ice, flattened by bin/layer [m]
    REAL(lp), ALLOCATABLE :: wtr_pot(:,:)          !< mean soil water potential, flattened by bin/layer [MPa]
  END TYPE t_sb_spinup_forcing_climatology

  !> Coupled-model forcing values preserved while the soil-only replay is active.
  TYPE :: t_sb_spinup_live_forcing
    REAL(lp), ALLOCATABLE :: litterfall(:,:,:) !< vegetation litterfall by source, element, and cell [mol m-2 timestep-1]
    REAL(lp), ALLOCATABLE :: seed_litterfall(:,:) !< seed-bed litterfall by element and cell [mol m-2 timestep-1]
    REAL(lp), ALLOCATABLE :: fine_root_c(:)       !< fine-root carbon by cell [mol C m-2]
    REAL(lp), ALLOCATABLE :: f_n_demand(:)        !< plant nitrogen-demand fraction by cell [-]
    REAL(lp), ALLOCATABLE :: f_p_demand(:)        !< plant phosphorus-demand fraction by cell [-]
    REAL(lp), ALLOCATABLE :: nhx_deposition(:)    !< reduced-N deposition by cell [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: noy_deposition(:)    !< oxidized-N deposition by cell [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: nhx_n15_deposition(:) !< reduced-15N deposition by cell [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: noy_n15_deposition(:) !< oxidized-15N deposition by cell [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: p_deposition(:)      !< phosphorus deposition by cell [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: growth_p_limit(:)    !< plant P limitation based on N demand by cell [-]
    REAL(lp), ALLOCATABLE :: root_fraction(:,:) !< root fraction by cell and soil layer [-]
    REAL(lp), ALLOCATABLE :: wtr_down(:,:)      !< downward water-transfer fraction by cell and layer [-]
    REAL(lp), ALLOCATABLE :: wtr_lat_loss(:,:)  !< lateral water-loss fraction by cell and layer [-]
    REAL(lp), ALLOCATABLE :: t_soil(:,:)        !< soil temperature by cell and layer [K]
    REAL(lp), ALLOCATABLE :: wtr_soil(:,:)      !< liquid soil water by cell and layer [m]
    REAL(lp), ALLOCATABLE :: ice_soil(:,:)      !< soil ice by cell and layer [m]
    REAL(lp), ALLOCATABLE :: wtr_pot(:,:)       !< soil water potential by cell and layer [MPa]
    REAL(lp), ALLOCATABLE :: net_biosphere_production(:) !< host-timestep net biosphere production [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: net_biosphere_production_c13(:) !< host-timestep 13C NBP [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: net_biosphere_production_c14(:) !< host-timestep 14C NBP [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: biological_n_fixation(:) !< host-timestep biological N fixation [micro-mol m-2 s-1]
  END TYPE t_sb_spinup_live_forcing

CONTAINS

  ! ===================================================================================================== !
  !> Identify the next accelerator event and whether its preceding forcing window is currently active.
  SUBROUTINE get_sb_spinup_accelerator_forcing_window(tile, options, l_collect, event_year, execution_index)

    USE mo_jsb_process_class, ONLY: SB_
    USE mo_jsb_time,          ONLY: get_year
    dsl4jsb_Use_config(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    LOGICAL,                    INTENT(out)   :: l_collect
    INTEGER,                    INTENT(out)   :: event_year
    INTEGER,                    INTENT(out)   :: execution_index

    dsl4jsb_Def_config(SB_)

    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: current_year
    INTEGER :: i_event

    l_collect = .FALSE.
    event_year = 0
    execution_index = 0

    model => Get_model(tile%owner_model_id)
    IF (.NOT. model%config%l_spinup) RETURN
    dsl4jsb_Get_config(SB_)

    current_year = get_year(options%current_datetime)
    DO i_event = 1, dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_max_executions
      event_year = dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_start_year + &
        & (i_event - 1) * dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_frequency

      IF (current_year <= event_year) THEN
        execution_index = i_event
        ! Only complete model years immediately preceding an event contribute to its climatology.
        l_collect = current_year < event_year .AND. &
          & current_year >= event_year - dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_forcing_years
        RETURN
      END IF
    END DO

    event_year = 0
    execution_index = 0

  END SUBROUTINE get_sb_spinup_accelerator_forcing_window

  ! ===================================================================================================== !
  !> Map the current point in the model year to a fixed climatological phase bin.
  SUBROUTINE get_sb_spinup_accelerator_forcing_bin(options, forcing_interval, forcing_bin, n_forcing_bins)

    USE mo_jsb_time, ONLY: get_year, get_year_day, get_year_length

    TYPE(t_jsb_task_options), INTENT(in)  :: options
    CHARACTER(len=*),         INTENT(in)  :: forcing_interval
    INTEGER,                  INTENT(out) :: forcing_bin
    INTEGER,                  INTENT(out) :: n_forcing_bins

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_sb_spinup_accelerator_forcing_bin'
    INTEGER  :: current_year
    INTEGER  :: day_index
    INTEGER  :: year_length
    REAL(lp) :: year_day

    SELECT CASE (TRIM(forcing_interval))
    CASE ('yearly')
      n_forcing_bins = 1
    CASE ('monthly')
      n_forcing_bins = 12
    CASE ('weekly')
      n_forcing_bins = 52
    CASE ('daily')
      n_forcing_bins = sb_spinup_accelerator_max_forcing_bins
    CASE DEFAULT
      CALL finish(TRIM(routine), 'Unsupported forward SB accelerator forcing interval.')
    END SELECT

    current_year = get_year(options%current_datetime)
    year_length = get_year_length(current_year)
    year_day = get_year_day(options%current_datetime)
#ifdef __QUINCY_STANDALONE__
    ! The standalone adapter reports a zero-based day, whereas mtime reports a one-based day of year.
    year_day = year_day + 1.0_lp
#endif

    ! Map whole model days with integer arithmetic. Dividing and multiplying the
    ! floating-point year phase can place exact day boundaries just below an
    ! integer and leave daily climatology bins unsampled.
    day_index = MAX(0, MIN(year_length - 1, FLOOR(year_day - 1.0_lp)))
    forcing_bin = MIN(n_forcing_bins, day_index * n_forcing_bins / year_length + 1)

  END SUBROUTINE get_sb_spinup_accelerator_forcing_bin

  ! ===================================================================================================== !
  !> Accumulate one normal-model timestep into the selected climatological forcing bin.
  SUBROUTINE collect_sb_spinup_accelerator_forcing(tile, options)

    USE mo_jsb_process_class, ONLY: SB_, VEG_, HYDRO_, SSE_, A2L_
    USE mo_lnd_bgcm_idx,       ONLY: ix_fine_root, ixC
    USE mo_lnd_bgcm_store,       ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_store_class, ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID, VEG_BGCM_SEED_BED_LITTERFALL_ID
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
    LOGICAL :: l_collect
    INTEGER :: event_year
    INTEGER :: execution_index
    INTEGER :: forcing_bin
    INTEGER :: n_forcing_bins
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    INTEGER :: i_source
    INTEGER :: i_element
    INTEGER :: i_soil
    INTEGER :: litter_entry
    INTEGER :: profile_entry
    INTEGER :: n_veg_sources
    INTEGER :: n_elements
    INTEGER :: nsoil
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

    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    model => Get_model(tile%owner_model_id)
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_config(SB_)
    CALL get_sb_spinup_accelerator_forcing_window(tile, options, l_collect, event_year, execution_index)
    IF (.NOT. l_collect) RETURN

    CALL get_sb_spinup_accelerator_forcing_bin(options, &
      & dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_forcing_interval, forcing_bin, n_forcing_bins)

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

    IF (event_year <= 0 .OR. execution_index <= 0 .OR. &
        & forcing_bin < 1 .OR. forcing_bin > n_forcing_bins) THEN
      CALL finish(TRIM(routine), 'Invalid forward SB accelerator forcing-window state.')
    END IF

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
      & sb_spinup_accelerator_growth_p_limit_sum(:, forcing_bin) + growth_p_limit_based_on_n_mavg * options%dtime

    n_veg_sources = SIZE(veg_litterfall_mt, 1)
    n_elements = SIZE(veg_litterfall_mt, 2)
    IF (SIZE(seed_bed_litterfall_mt, 1) /= n_elements) THEN
      CALL finish(TRIM(routine), 'Vegetation and seed-bed litterfall element dimensions differ.')
    END IF
    IF (n_forcing_bins * (n_veg_sources + 1) * n_elements > &
        & SIZE(sb_spinup_accelerator_litterfall_sum, 2)) THEN
      CALL finish(TRIM(routine), 'Forward SB accelerator litterfall storage is too small.')
    END IF

    DO i_element = 1, n_elements
      DO i_source = 1, n_veg_sources
        litter_entry = forcing_bin + n_forcing_bins * &
          & ((i_source - 1) + (n_veg_sources + 1) * (i_element - 1))
        sb_spinup_accelerator_litterfall_sum(:, litter_entry) = &
          & sb_spinup_accelerator_litterfall_sum(:, litter_entry) + &
          & veg_litterfall_mt(i_source, i_element, :)
      END DO

      litter_entry = forcing_bin + n_forcing_bins * &
        & (n_veg_sources + (n_veg_sources + 1) * (i_element - 1))
      sb_spinup_accelerator_litterfall_sum(:, litter_entry) = &
        & sb_spinup_accelerator_litterfall_sum(:, litter_entry) + seed_bed_litterfall_mt(i_element, :)
    END DO

    nsoil = SIZE(root_fraction_sl, 2)
    IF (SIZE(frac_wtr_transp_down_sl, 2) /= nsoil .OR. SIZE(frac_w_lat_loss_sl, 2) /= nsoil .OR. &
        & SIZE(wtr_soil_sl, 2) /= nsoil .OR. SIZE(ice_soil_sl, 2) /= nsoil .OR. &
        & SIZE(wtr_soil_pot_sl, 2) /= nsoil .OR. SIZE(t_soil_sl, 2) /= nsoil .OR. &
        & n_forcing_bins * nsoil > SIZE(sb_spinup_accelerator_root_fraction_sum, 2)) THEN
      CALL finish(TRIM(routine), 'Forward SB accelerator soil-profile forcing dimensions differ.')
    END IF

    DO i_soil = 1, nsoil
      profile_entry = forcing_bin + n_forcing_bins * (i_soil - 1)
      sb_spinup_accelerator_root_fraction_sum(:, profile_entry) = &
        & sb_spinup_accelerator_root_fraction_sum(:, profile_entry) + root_fraction_sl(:, i_soil) * options%dtime
      sb_spinup_accelerator_wtr_down_sum(:, profile_entry) = &
        & sb_spinup_accelerator_wtr_down_sum(:, profile_entry) + frac_wtr_transp_down_sl(:, i_soil) * options%dtime
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

    ! Subsequent bookkeeping groups will use the same bin ordering and elapsed-time denominator.

  END SUBROUTINE collect_sb_spinup_accelerator_forcing

  ! ===================================================================================================== !
  !> Decide whether the current spin-up timestep starts a scheduled accelerator event.
  SUBROUTINE check_sb_spinup_accelerator_schedule(tile, options, l_accelerate, execution_index)

    USE mo_jsb_process_class, ONLY: SB_
    USE mo_jsb_time,          ONLY: get_year, is_newyear
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    LOGICAL,                    INTENT(out)   :: l_accelerate
    INTEGER,                    INTENT(out)   :: execution_index

    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)

    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: current_year
    INTEGER :: years_since_start
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_execution_count

    l_accelerate = .FALSE.
    execution_index = 0

    model => Get_model(tile%owner_model_id)
    IF (.NOT. model%config%l_spinup) RETURN
    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_execution_count)
    IF (.NOT. is_newyear(options%current_datetime, options%dtime)) RETURN

    current_year = get_year(options%current_datetime)
    IF (current_year < dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_start_year) RETURN

    years_since_start = current_year - dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_start_year
    IF (MOD(years_since_start, dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_frequency) /= 0) RETURN

    execution_index = years_since_start / dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_frequency + 1
    IF (execution_index > dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_max_executions) THEN
      execution_index = 0
      RETURN
    END IF
    IF (ALL(NINT(sb_spinup_accelerator_execution_count) >= execution_index)) THEN
      execution_index = 0
      RETURN
    END IF

    l_accelerate = .TRUE.

  END SUBROUTINE check_sb_spinup_accelerator_schedule

  ! ===================================================================================================== !
  !> Record successful completion of a scheduled accelerator event in restart memory.
  SUBROUTINE mark_sb_spinup_accelerator_execution(tile, options, execution_index)

    USE mo_jsb_process_class, ONLY: SB_
    dsl4jsb_Use_memory(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    INTEGER,                    INTENT(in)    :: execution_index

    dsl4jsb_Def_memory(SB_)

    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_execution_count

    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_execution_count)

    sb_spinup_accelerator_execution_count(:) = REAL(execution_index, lp)

  END SUBROUTINE mark_sb_spinup_accelerator_execution

  ! ===================================================================================================== !
  !> Clear forcing-bin accumulators after an accelerator execution or forcing-window reset.
  SUBROUTINE reset_sb_spinup_accelerator_forcing(tile, options)

    USE mo_jsb_process_class, ONLY: SB_
    dsl4jsb_Use_memory(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(SB_)

    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
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

    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_memory(SB_)
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

    sb_spinup_accelerator_forcing_time(:,:) = 0.0_lp
    sb_spinup_accelerator_litterfall_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_fine_root_c_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_f_n_demand_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_f_p_demand_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_nhx_deposition_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_noy_deposition_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_nhx_n15_deposition_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_noy_n15_deposition_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_p_deposition_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_growth_p_limit_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_root_fraction_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_wtr_down_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_wtr_lat_loss_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_t_soil_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_wtr_soil_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_ice_soil_sum(:,:) = 0.0_lp
    sb_spinup_accelerator_wtr_pot_sum(:,:) = 0.0_lp

  END SUBROUTINE reset_sb_spinup_accelerator_forcing

  ! ===================================================================================================== !
  !> Apply forward soil-model relaxation at a scheduled spin-up event.
  SUBROUTINE apply_sb_spinup_accelerator(tile, options)

    USE mo_jsb_process_class, ONLY: SB_
    USE mo_jsb_time,          ONLY: executed_sb_model_spinup_accelerator_this_timestep
    USE mo_quincy_model_config, ONLY: QLAND, QSOIL
    dsl4jsb_Use_config(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(SB_)

    TYPE(t_jsb_model), POINTER :: model
    CHARACTER(len=*), PARAMETER :: routine = modname//':apply_sb_spinup_accelerator'
    LOGICAL :: l_accelerate
    LOGICAL :: l_converged
    INTEGER :: execution_index
    INTEGER :: represented_years
    REAL(lp) :: max_relative_change
    REAL(lp) :: max_total_c_absolute_change
    REAL(lp) :: max_total_c_relative_change
    REAL(lp), ALLOCATABLE :: pool_before(:,:,:,:)
    REAL(lp), ALLOCATABLE :: pool_after(:,:,:,:)

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(SB_)
    IF (.NOT. dsl4jsb_Config(SB_)%flag_sb_model_spinup_accelerator) RETURN
    IF (model%config%qmodel_id /= QLAND .AND. model%config%qmodel_id /= QSOIL) RETURN
    IF (.NOT. model%config%l_spinup) RETURN
    IF (model%config%flag_slow_sb_pool_spinup_accelerator) THEN
      CALL finish(TRIM(routine), 'The forward and slow-pool SB spin-up accelerators are mutually exclusive.')
    END IF

    CALL collect_sb_spinup_accelerator_forcing(tile, options)
    CALL check_sb_spinup_accelerator_schedule(tile, options, l_accelerate, execution_index)

    IF (.NOT. l_accelerate) RETURN
    IF (options%iblk == 1) THEN
      WRITE(message_text,*) 'Starting forward SB spin-up accelerator execution ', execution_index, &
        & ' with forcing interval ', TRIM(dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_forcing_interval), '.'
      CALL message(TRIM(routine), TRIM(message_text))
    END IF

    IF (model%config%flag_spinup_debug) THEN
      CALL snapshot_sb_spinup_pools(tile, options, pool_before)
    END IF

    CALL execute_sb_spinup_accelerator(tile, options, dsl4jsb_Config(SB_)%sb_model_scheme, &
      & execution_index, represented_years, max_relative_change, max_total_c_absolute_change, &
      & max_total_c_relative_change, l_converged)
    executed_sb_model_spinup_accelerator_this_timestep = .TRUE.

    IF (model%config%flag_spinup_debug) THEN
      CALL snapshot_sb_spinup_pools(tile, options, pool_after)
      CALL write_sb_spinup_debug(tile, options, execution_index, represented_years, &
        & pool_before, pool_after)
    END IF

    IF (options%iblk == 1) THEN
      IF (l_converged) THEN
        WRITE(message_text,*) 'Finished forward SB spin-up accelerator execution ', execution_index, &
          & '; represented years = ', represented_years, '; converged = ', l_converged, &
          & '; maximum annual relative pool change = ', max_relative_change, &
          & '; maximum annual total-C absolute and relative changes = ', &
          & max_total_c_absolute_change, ' mol C m-2; ', max_total_c_relative_change
      ELSE
        WRITE(message_text,*) 'WARNING: forward SB spin-up accelerator execution ', execution_index, &
          & ' reached sb_model_spinup_accelerator_max_relaxation_years without convergence; continuing normal simulation.', &
          & ' Represented years = ', represented_years, '; maximum annual relative pool change = ', max_relative_change, &
          & '; maximum annual total-C absolute and relative changes = ', &
          & max_total_c_absolute_change, ' mol C m-2; ', max_total_c_relative_change
      END IF
      CALL message(TRIM(routine), TRIM(message_text))
    END IF

  END SUBROUTINE apply_sb_spinup_accelerator

  ! ===================================================================================================== !
  !> Execute one scheduled accelerator event once forcing storage and scheduling are implemented.
  SUBROUTINE execute_sb_spinup_accelerator(tile, options, sb_model_scheme, execution_index, &
      & represented_years, max_relative_change, max_total_c_absolute_change, &
      & max_total_c_relative_change, l_converged)

    USE mo_sb_util, ONLY: reset_sb_fluxes
    USE mo_jsb_time, ONLY: execute_quincy_sb_analytical_accelerator_this_timestep

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    CHARACTER(len=*),           INTENT(in)    :: sb_model_scheme
    INTEGER,                    INTENT(in)    :: execution_index
    INTEGER,                    INTENT(out)   :: represented_years
    REAL(lp),                   INTENT(out)   :: max_relative_change
    REAL(lp),                   INTENT(out)   :: max_total_c_absolute_change
    REAL(lp),                   INTENT(out)   :: max_total_c_relative_change
    LOGICAL,                    INTENT(out)   :: l_converged

    TYPE(t_sb_spinup_forcing_climatology) :: forcing
    TYPE(t_sb_spinup_live_forcing) :: live_forcing
    LOGICAL :: l_analytical_event

    ! This sequence is intentionally visible before the numerical bodies are implemented.
    CALL prepare_sb_spinup_accelerator(tile, options, forcing, live_forcing)
    ! The analytical target has already been applied by update_sb_pools. Mask its
    ! event flag while replay repeatedly calls update_sb_pools in pseudo-time.
    l_analytical_event = execute_quincy_sb_analytical_accelerator_this_timestep
    execute_quincy_sb_analytical_accelerator_this_timestep = .FALSE.
    CALL run_sb_forward_relaxation(tile, options, sb_model_scheme, forcing, execution_index, represented_years, &
      & max_relative_change, max_total_c_absolute_change, max_total_c_relative_change, l_converged)
    execute_quincy_sb_analytical_accelerator_this_timestep = l_analytical_event
    CALL reset_sb_fluxes(tile, options)
    CALL finalize_sb_spinup_accelerator(tile, options, live_forcing)
    CALL reset_sb_spinup_accelerator_forcing(tile, options)
    CALL mark_sb_spinup_accelerator_execution(tile, options, execution_index)

  END SUBROUTINE execute_sb_spinup_accelerator

  ! ===================================================================================================== !
  !> Validate accumulated forcing and reconstruct the event-local climatology.
  SUBROUTINE prepare_sb_spinup_accelerator(tile, options, forcing, live_forcing)

    USE mo_jsb_process_class, ONLY: SB_, VEG_, HYDRO_, SSE_, A2L_
    USE mo_lnd_bgcm_idx,       ONLY: ix_fine_root, ixC
    USE mo_lnd_bgcm_store,       ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_store_class, ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID, VEG_BGCM_SEED_BED_LITTERFALL_ID
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(A2L_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    TYPE(t_sb_spinup_forcing_climatology), INTENT(out) :: forcing
    TYPE(t_sb_spinup_live_forcing), INTENT(out) :: live_forcing

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
    CHARACTER(len=*), PARAMETER :: routine = modname//':prepare_sb_spinup_accelerator'
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    INTEGER :: i_entry
    INTEGER :: forcing_bin
    INTEGER :: n_cells
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
    dsl4jsb_Real2D_onChunk :: net_biosphere_production
    dsl4jsb_Real2D_onChunk :: net_biosphere_production_c13
    dsl4jsb_Real2D_onChunk :: net_biosphere_production_c14
    dsl4jsb_Real2D_onChunk :: biological_n_fixation
    dsl4jsb_Real2D_onChunk :: f_n_demand
    dsl4jsb_Real2D_onChunk :: f_p_demand
    dsl4jsb_Real2D_onChunk :: growth_p_limit_based_on_n_mavg
    dsl4jsb_Real2D_onChunk :: nhx_deposition
    dsl4jsb_Real2D_onChunk :: noy_deposition
    dsl4jsb_Real2D_onChunk :: nhx_n15_deposition
    dsl4jsb_Real2D_onChunk :: noy_n15_deposition
    dsl4jsb_Real2D_onChunk :: p_deposition

    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    model => Get_model(tile%owner_model_id)
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
    dsl4jsb_Get_var2D_onChunk(VEG_, net_biosphere_production)
    dsl4jsb_Get_var2D_onChunk(VEG_, net_biosphere_production_c13)
    dsl4jsb_Get_var2D_onChunk(VEG_, net_biosphere_production_c14)
    dsl4jsb_Get_var2D_onChunk(VEG_, biological_n_fixation)
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

    n_cells = SIZE(sb_spinup_accelerator_forcing_time, 1)
    forcing%n_forcing_bins = SIZE(sb_spinup_accelerator_forcing_time, 2)
    forcing%n_litter_sources = SIZE(veg_litterfall_mt, 1)
    forcing%n_elements = SIZE(veg_litterfall_mt, 2)
    IF (forcing%n_forcing_bins <= 0 .OR. &
        & MOD(SIZE(sb_spinup_accelerator_root_fraction_sum, 2), forcing%n_forcing_bins) /= 0) THEN
      CALL finish(TRIM(routine), 'Invalid stored forcing dimensions.')
    END IF
    forcing%n_soil_layers = SIZE(sb_spinup_accelerator_root_fraction_sum, 2) / forcing%n_forcing_bins

    IF (ANY(sb_spinup_accelerator_forcing_time <= 0.0_lp)) THEN
      CALL finish(TRIM(routine), 'At least one climatological forcing bin has no sampled model time.')
    END IF

    ALLOCATE(forcing%forcing_time(n_cells, forcing%n_forcing_bins), &
      & forcing%litterfall_rate(n_cells, SIZE(sb_spinup_accelerator_litterfall_sum, 2)), &
      & forcing%fine_root_c(n_cells, forcing%n_forcing_bins), &
      & forcing%f_n_demand(n_cells, forcing%n_forcing_bins), &
      & forcing%f_p_demand(n_cells, forcing%n_forcing_bins), &
      & forcing%nhx_deposition(n_cells, forcing%n_forcing_bins), &
      & forcing%noy_deposition(n_cells, forcing%n_forcing_bins), &
      & forcing%nhx_n15_deposition(n_cells, forcing%n_forcing_bins), &
      & forcing%noy_n15_deposition(n_cells, forcing%n_forcing_bins), &
      & forcing%p_deposition(n_cells, forcing%n_forcing_bins), &
      & forcing%growth_p_limit(n_cells, forcing%n_forcing_bins), &
      & forcing%root_fraction(n_cells, SIZE(sb_spinup_accelerator_root_fraction_sum, 2)), &
      & forcing%wtr_down(n_cells, SIZE(sb_spinup_accelerator_wtr_down_sum, 2)), &
      & forcing%wtr_lat_loss(n_cells, SIZE(sb_spinup_accelerator_wtr_lat_loss_sum, 2)), &
      & forcing%t_soil(n_cells, SIZE(sb_spinup_accelerator_t_soil_sum, 2)), &
      & forcing%wtr_soil(n_cells, SIZE(sb_spinup_accelerator_wtr_soil_sum, 2)), &
      & forcing%ice_soil(n_cells, SIZE(sb_spinup_accelerator_ice_soil_sum, 2)), &
      & forcing%wtr_pot(n_cells, SIZE(sb_spinup_accelerator_wtr_pot_sum, 2)))

    forcing%forcing_time = sb_spinup_accelerator_forcing_time
    forcing%fine_root_c = sb_spinup_accelerator_fine_root_c_sum / forcing%forcing_time
    forcing%f_n_demand = sb_spinup_accelerator_f_n_demand_sum / forcing%forcing_time
    forcing%f_p_demand = sb_spinup_accelerator_f_p_demand_sum / forcing%forcing_time
    forcing%nhx_deposition = sb_spinup_accelerator_nhx_deposition_sum / forcing%forcing_time
    forcing%noy_deposition = sb_spinup_accelerator_noy_deposition_sum / forcing%forcing_time
    forcing%nhx_n15_deposition = sb_spinup_accelerator_nhx_n15_deposition_sum / forcing%forcing_time
    forcing%noy_n15_deposition = sb_spinup_accelerator_noy_n15_deposition_sum / forcing%forcing_time
    forcing%p_deposition = sb_spinup_accelerator_p_deposition_sum / forcing%forcing_time
    forcing%growth_p_limit = sb_spinup_accelerator_growth_p_limit_sum / forcing%forcing_time
    DO i_entry = 1, SIZE(forcing%litterfall_rate, 2)
      forcing_bin = MOD(i_entry - 1, forcing%n_forcing_bins) + 1
      forcing%litterfall_rate(:, i_entry) = sb_spinup_accelerator_litterfall_sum(:, i_entry) / &
        & forcing%forcing_time(:, forcing_bin)
    END DO

    DO i_entry = 1, SIZE(forcing%root_fraction, 2)
      forcing_bin = MOD(i_entry - 1, forcing%n_forcing_bins) + 1
      forcing%root_fraction(:, i_entry) = sb_spinup_accelerator_root_fraction_sum(:, i_entry) / &
        & forcing%forcing_time(:, forcing_bin)
      forcing%wtr_down(:, i_entry) = sb_spinup_accelerator_wtr_down_sum(:, i_entry) / &
        & forcing%forcing_time(:, forcing_bin)
      forcing%wtr_lat_loss(:, i_entry) = sb_spinup_accelerator_wtr_lat_loss_sum(:, i_entry) / &
        & forcing%forcing_time(:, forcing_bin)
      forcing%t_soil(:, i_entry) = sb_spinup_accelerator_t_soil_sum(:, i_entry) / &
        & forcing%forcing_time(:, forcing_bin)
      forcing%wtr_soil(:, i_entry) = sb_spinup_accelerator_wtr_soil_sum(:, i_entry) / &
        & forcing%forcing_time(:, forcing_bin)
      forcing%ice_soil(:, i_entry) = sb_spinup_accelerator_ice_soil_sum(:, i_entry) / &
        & forcing%forcing_time(:, forcing_bin)
      forcing%wtr_pot(:, i_entry) = sb_spinup_accelerator_wtr_pot_sum(:, i_entry) / &
        & forcing%forcing_time(:, forcing_bin)
    END DO

    ALLOCATE(live_forcing%litterfall, SOURCE=veg_litterfall_mt)
    ALLOCATE(live_forcing%seed_litterfall, SOURCE=seed_bed_litterfall_mt)
    ALLOCATE(live_forcing%fine_root_c, SOURCE=veg_pool_mt(ix_fine_root, ixC, :))
    ALLOCATE(live_forcing%f_n_demand, SOURCE=f_n_demand)
    ALLOCATE(live_forcing%f_p_demand, SOURCE=f_p_demand)
    ALLOCATE(live_forcing%nhx_deposition, SOURCE=nhx_deposition)
    ALLOCATE(live_forcing%noy_deposition, SOURCE=noy_deposition)
    ALLOCATE(live_forcing%nhx_n15_deposition, SOURCE=nhx_n15_deposition)
    ALLOCATE(live_forcing%noy_n15_deposition, SOURCE=noy_n15_deposition)
    ALLOCATE(live_forcing%p_deposition, SOURCE=p_deposition)
    ALLOCATE(live_forcing%growth_p_limit, SOURCE=growth_p_limit_based_on_n_mavg)
    ALLOCATE(live_forcing%root_fraction, SOURCE=root_fraction_sl)
    ALLOCATE(live_forcing%wtr_down, SOURCE=frac_wtr_transp_down_sl)
    ALLOCATE(live_forcing%wtr_lat_loss, SOURCE=frac_w_lat_loss_sl)
    ALLOCATE(live_forcing%t_soil, SOURCE=t_soil_sl)
    ALLOCATE(live_forcing%wtr_soil, SOURCE=wtr_soil_sl)
    ALLOCATE(live_forcing%ice_soil, SOURCE=ice_soil_sl)
    ALLOCATE(live_forcing%wtr_pot, SOURCE=wtr_soil_pot_sl)
    ALLOCATE(live_forcing%net_biosphere_production, SOURCE=net_biosphere_production)
    ALLOCATE(live_forcing%net_biosphere_production_c13, SOURCE=net_biosphere_production_c13)
    ALLOCATE(live_forcing%net_biosphere_production_c14, SOURCE=net_biosphere_production_c14)
    ALLOCATE(live_forcing%biological_n_fixation, SOURCE=biological_n_fixation)

  END SUBROUTINE prepare_sb_spinup_accelerator

  ! ===================================================================================================== !
  !> Replay the bookkept climatology until the maximum pseudo-time or convergence criterion is reached.
  SUBROUTINE run_sb_forward_relaxation(tile, options, sb_model_scheme, forcing, execution_index, represented_years, &
      & max_relative_change, max_total_c_absolute_change, max_total_c_relative_change, l_converged)

    USE mo_jsb_process_class, ONLY: SB_
    dsl4jsb_Use_config(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    CHARACTER(len=*),           INTENT(in)    :: sb_model_scheme
    TYPE(t_sb_spinup_forcing_climatology), INTENT(in) :: forcing
    INTEGER,                    INTENT(in)    :: execution_index
    INTEGER,                    INTENT(out)   :: represented_years
    REAL(lp),                   INTENT(out)   :: max_relative_change
    REAL(lp),                   INTENT(out)   :: max_total_c_absolute_change
    REAL(lp),                   INTENT(out)   :: max_total_c_relative_change
    LOGICAL,                    INTENT(out)   :: l_converged

    dsl4jsb_Def_config(SB_)

    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: consecutive_converged_years
    INTEGER :: i_year
    INTEGER :: last_debug_year
    REAL(lp), ALLOCATABLE :: previous_pool(:,:,:,:)

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(SB_)
    CALL snapshot_sb_spinup_pools(tile, options, previous_pool)

    last_debug_year = -1
    IF (model%config%flag_spinup_debug) THEN
      CALL write_sb_spinup_forcing_debug(options, execution_index, forcing)
      CALL write_sb_spinup_snapshot_debug(tile, options, execution_index, 0)
      last_debug_year = 0
    END IF

    represented_years = 0
    consecutive_converged_years = 0
    max_relative_change = HUGE(1.0_lp)
    max_total_c_absolute_change = HUGE(1.0_lp)
    max_total_c_relative_change = HUGE(1.0_lp)
    l_converged = .FALSE.

    DO i_year = 1, dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_max_relaxation_years
      represented_years = i_year
      SELECT CASE (TRIM(sb_model_scheme))
      CASE ('jsm')
        CALL run_jsm_forward_relaxation(tile, options, forcing)
      CASE ('simple_1d')
        CALL run_simple_1d_forward_relaxation(tile, options, forcing)
      CASE DEFAULT
        CALL finish(modname//':run_sb_forward_relaxation', 'Unsupported soil-biogeochemistry model scheme.')
      END SELECT

      CALL check_sb_spinup_convergence(tile, options, previous_pool, consecutive_converged_years, &
        & max_relative_change, max_total_c_absolute_change, max_total_c_relative_change, l_converged)
      IF (model%config%flag_spinup_debug .AND. is_sb_spinup_debug_year(i_year)) THEN
        CALL write_sb_spinup_snapshot_debug(tile, options, execution_index, i_year, &
          & max_relative_change, max_total_c_absolute_change, max_total_c_relative_change, l_converged)
        last_debug_year = i_year
      END IF
      IF (l_converged) EXIT
    END DO

    IF (model%config%flag_spinup_debug .AND. represented_years /= last_debug_year) THEN
      CALL write_sb_spinup_snapshot_debug(tile, options, execution_index, represented_years, &
        & max_relative_change, max_total_c_absolute_change, max_total_c_relative_change, l_converged)
    END IF

  END SUBROUTINE run_sb_forward_relaxation

  ! ===================================================================================================== !
  !> Replay forcing through the JSM process sequence without vegetation or atmospheric model updates.
  SUBROUTINE run_jsm_forward_relaxation(tile, options, forcing)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    TYPE(t_sb_spinup_forcing_climatology), INTENT(in) :: forcing

    CALL run_sb_climatological_year(tile, options, forcing, 'jsm')

  END SUBROUTINE run_jsm_forward_relaxation

  ! ===================================================================================================== !
  !> Replay forcing through the simple_1d process sequence without vegetation or atmospheric model updates.
  SUBROUTINE run_simple_1d_forward_relaxation(tile, options, forcing)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    TYPE(t_sb_spinup_forcing_climatology), INTENT(in) :: forcing

    CALL run_sb_climatological_year(tile, options, forcing, 'simple_1d')

  END SUBROUTINE run_simple_1d_forward_relaxation

  ! ===================================================================================================== !
  !> Advance the selected soil model through one 365-day climatological replay year.
  SUBROUTINE run_sb_climatological_year(tile, options, forcing, sb_model_scheme)

    USE mo_jsb_math_constants,   ONLY: one_day, one_year
    USE mo_jsb_process_class,    ONLY: SB_
    USE mo_sb_util,               ONLY: reset_sb_fluxes, calculate_time_average_soilbiogeochemistry
    USE mo_q_sb_rate_modifier,    ONLY: update_sb_rate_modifier
    USE mo_q_sb_jsm_main,         ONLY: update_sb_jsm
    USE mo_q_sb_ssm_main,         ONLY: update_sb_simple_model
    USE mo_q_sb_update_pools,     ONLY: update_sb_pools
    dsl4jsb_Use_config(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    TYPE(t_sb_spinup_forcing_climatology), INTENT(in) :: forcing
    CHARACTER(len=*),           INTENT(in)    :: sb_model_scheme

    dsl4jsb_Def_config(SB_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_task_options) :: replay_options
    CHARACTER(len=*), PARAMETER :: routine = modname//':run_sb_climatological_year'
    INTEGER :: forcing_bin
    INTEGER :: i_step
    INTEGER :: n_steps_per_bin
    REAL(lp) :: bin_duration
    REAL(lp) :: requested_dtime

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(SB_)

    requested_dtime = dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_timestep_days * one_day
    bin_duration = one_year * one_day / REAL(forcing%n_forcing_bins, lp)
    n_steps_per_bin = MAX(1, CEILING(bin_duration / requested_dtime))

    replay_options = options
    replay_options%dtime = bin_duration / REAL(n_steps_per_bin, lp)

    DO forcing_bin = 1, forcing%n_forcing_bins
      DO i_step = 1, n_steps_per_bin
        CALL load_sb_spinup_forcing_bin(tile, replay_options, forcing, forcing_bin)
        CALL reset_sb_fluxes(tile, replay_options)
        CALL update_sb_rate_modifier(tile, replay_options)

        SELECT CASE (TRIM(sb_model_scheme))
        CASE ('jsm')
          CALL update_sb_jsm(tile, replay_options)
        CASE ('simple_1d')
          CALL update_sb_simple_model(tile, replay_options)
        CASE DEFAULT
          CALL finish(TRIM(routine), 'Unsupported soil-biogeochemistry model scheme.')
        END SELECT

        CALL update_sb_pools(tile, replay_options, l_soil_model_forcing=.TRUE.)
        CALL calculate_time_average_soilbiogeochemistry(tile, replay_options)
      END DO
    END DO

  END SUBROUTINE run_sb_climatological_year

  ! ===================================================================================================== !
  !> Load one annual, monthly, weekly, or daily climatological forcing bin into the replay forcing view.
  SUBROUTINE load_sb_spinup_forcing_bin(tile, options, forcing, forcing_bin)

    USE mo_jsb_process_class, ONLY: VEG_, HYDRO_, SSE_, A2L_
    USE mo_lnd_bgcm_idx,       ONLY: ix_fine_root, ixC
    USE mo_lnd_bgcm_store,       ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_store_class, ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID, VEG_BGCM_SEED_BED_LITTERFALL_ID
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(A2L_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    TYPE(t_sb_spinup_forcing_climatology), INTENT(in) :: forcing
    INTEGER,                    INTENT(in)    :: forcing_bin

    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt1L2D :: seed_bed_litterfall_mt

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store
    CHARACTER(len=*), PARAMETER :: routine = modname//':load_sb_spinup_forcing_bin'
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    INTEGER :: i_element
    INTEGER :: i_soil
    INTEGER :: i_source
    INTEGER :: litter_entry
    INTEGER :: profile_entry
    INTEGER :: n_elements
    INTEGER :: n_sources
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

    IF (forcing_bin < 1 .OR. forcing_bin > forcing%n_forcing_bins) THEN
      CALL finish(TRIM(routine), 'Requested forcing bin is outside the prepared climatology.')
    END IF

    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    model => Get_model(tile%owner_model_id)
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(A2L_)
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

    veg_pool_mt(ix_fine_root, ixC, :) = forcing%fine_root_c(:, forcing_bin)
    f_n_demand = forcing%f_n_demand(:, forcing_bin)
    f_p_demand = forcing%f_p_demand(:, forcing_bin)
    growth_p_limit_based_on_n_mavg = forcing%growth_p_limit(:, forcing_bin)
    nhx_deposition = forcing%nhx_deposition(:, forcing_bin)
    noy_deposition = forcing%noy_deposition(:, forcing_bin)
    nhx_n15_deposition = forcing%nhx_n15_deposition(:, forcing_bin)
    noy_n15_deposition = forcing%noy_n15_deposition(:, forcing_bin)
    p_deposition = forcing%p_deposition(:, forcing_bin)

    n_sources = SIZE(veg_litterfall_mt, 1)
    n_elements = SIZE(veg_litterfall_mt, 2)
    DO i_element = 1, n_elements
      DO i_source = 1, n_sources
        litter_entry = forcing_bin + forcing%n_forcing_bins * &
          & ((i_source - 1) + (n_sources + 1) * (i_element - 1))
        veg_litterfall_mt(i_source, i_element, :) = forcing%litterfall_rate(:, litter_entry) * options%dtime
      END DO
      litter_entry = forcing_bin + forcing%n_forcing_bins * &
        & (n_sources + (n_sources + 1) * (i_element - 1))
      seed_bed_litterfall_mt(i_element, :) = forcing%litterfall_rate(:, litter_entry) * options%dtime
    END DO

    DO i_soil = 1, forcing%n_soil_layers
      profile_entry = forcing_bin + forcing%n_forcing_bins * (i_soil - 1)
      root_fraction_sl(:, i_soil) = forcing%root_fraction(:, profile_entry)
      frac_wtr_transp_down_sl(:, i_soil) = forcing%wtr_down(:, profile_entry)
      frac_w_lat_loss_sl(:, i_soil) = forcing%wtr_lat_loss(:, profile_entry)
      t_soil_sl(:, i_soil) = forcing%t_soil(:, profile_entry)
      wtr_soil_sl(:, i_soil) = forcing%wtr_soil(:, profile_entry)
      ice_soil_sl(:, i_soil) = forcing%ice_soil(:, profile_entry)
      wtr_soil_pot_sl(:, i_soil) = forcing%wtr_pot(:, profile_entry)
    END DO

  END SUBROUTINE load_sb_spinup_forcing_bin

  ! ===================================================================================================== !
  !> Return true for the sparse pseudo-years included in the detailed accelerator diagnostics.
  LOGICAL FUNCTION is_sb_spinup_debug_year(represented_year)

    INTEGER, INTENT(in) :: represented_year

    is_sb_spinup_debug_year = ANY(sb_spinup_debug_years == represented_year)

  END FUNCTION is_sb_spinup_debug_year

  ! ===================================================================================================== !
  !> Write the event-local forcing climatology used by one forward SB acceleration.
  SUBROUTINE write_sb_spinup_forcing_debug(options, execution_index, forcing)

    USE mo_mpi, ONLY: p_pe_work_only

    TYPE(t_jsb_task_options), INTENT(in) :: options
    INTEGER, INTENT(in) :: execution_index
    TYPE(t_sb_spinup_forcing_climatology), INTENT(in) :: forcing

    CHARACTER(len=*), PARAMETER :: routine = modname//':write_sb_spinup_forcing_debug'
    INTEGER :: lun
    INTEGER :: ios
    INTEGER :: ic
    INTEGER :: forcing_bin
    INTEGER :: i_element
    INTEGER :: i_source
    INTEGER :: litter_entry
    INTEGER :: i_soil
    INTEGER :: profile_entry

    ! Limit this local-file diagnostic to one work rank and block. For global runs,
    ! it therefore represents the first local chunk rather than a global gather.
    IF (p_pe_work_only /= 0 .OR. options%iblk /= 1) RETURN

    IF (execution_index == 1) THEN
      OPEN(newunit=lun, file=sb_spinup_debug_filename, status='replace', action='write', iostat=ios)
    ELSE
      OPEN(newunit=lun, file=sb_spinup_debug_filename, status='unknown', position='append', &
        & action='write', iostat=ios)
    END IF
    IF (ios /= 0) CALL finish(TRIM(routine), 'Could not open '//sb_spinup_debug_filename)

    WRITE(lun,'(A)') ''
    WRITE(lun,*) 'FORCING_SUMMARY execution, n_bins, n_layers, n_litter_sources, n_elements = ', &
      & execution_index, forcing%n_forcing_bins, forcing%n_soil_layers, &
      & forcing%n_litter_sources, forcing%n_elements
    WRITE(lun,'(A)') 'FORCING_SCALAR_COLUMNS execution ic bin sampled_time_s fine_root_c_mol_m-2 '// &
      & 'f_n_demand f_p_demand nhx_dep_umol_m-2_s-1 noy_dep_umol_m-2_s-1 '// &
      & 'nhx_n15_dep_umol_m-2_s-1 noy_n15_dep_umol_m-2_s-1 p_dep_umol_m-2_s-1 growth_p_limit'
    DO ic = 1, SIZE(forcing%forcing_time, 1)
      DO forcing_bin = 1, forcing%n_forcing_bins
        WRITE(lun,*) 'FORCING_SCALAR ', execution_index, ic, forcing_bin, &
          & forcing%forcing_time(ic,forcing_bin), forcing%fine_root_c(ic,forcing_bin), &
          & forcing%f_n_demand(ic,forcing_bin), forcing%f_p_demand(ic,forcing_bin), &
          & forcing%nhx_deposition(ic,forcing_bin), forcing%noy_deposition(ic,forcing_bin), &
          & forcing%nhx_n15_deposition(ic,forcing_bin), forcing%noy_n15_deposition(ic,forcing_bin), &
          & forcing%p_deposition(ic,forcing_bin), forcing%growth_p_limit(ic,forcing_bin)
      END DO
    END DO

    WRITE(lun,'(A)') 'FORCING_LITTER_COLUMNS execution ic bin source element '// &
      & 'litterfall_rate_mol_m-2_s-1; source n_litter_sources+1 is seed-bed litterfall'
    DO ic = 1, SIZE(forcing%forcing_time, 1)
      DO i_element = 1, forcing%n_elements
        DO i_source = 1, forcing%n_litter_sources + 1
          DO forcing_bin = 1, forcing%n_forcing_bins
            litter_entry = forcing_bin + forcing%n_forcing_bins * &
              & ((i_source - 1) + (forcing%n_litter_sources + 1) * (i_element - 1))
            WRITE(lun,*) 'FORCING_LITTER ', execution_index, ic, forcing_bin, i_source, i_element, &
              & forcing%litterfall_rate(ic,litter_entry)
          END DO
        END DO
      END DO
    END DO

    WRITE(lun,'(A)') 'FORCING_SOIL_COLUMNS execution ic bin layer root_fraction wtr_down_fraction '// &
      & 'wtr_lateral_loss_fraction t_soil_K wtr_soil_m ice_soil_m wtr_potential_MPa'
    DO ic = 1, SIZE(forcing%forcing_time, 1)
      DO i_soil = 1, forcing%n_soil_layers
        DO forcing_bin = 1, forcing%n_forcing_bins
          profile_entry = forcing_bin + forcing%n_forcing_bins * (i_soil - 1)
          WRITE(lun,*) 'FORCING_SOIL ', execution_index, ic, forcing_bin, i_soil, &
            & forcing%root_fraction(ic,profile_entry), forcing%wtr_down(ic,profile_entry), &
            & forcing%wtr_lat_loss(ic,profile_entry), forcing%t_soil(ic,profile_entry), &
            & forcing%wtr_soil(ic,profile_entry), forcing%ice_soil(ic,profile_entry), &
            & forcing%wtr_pot(ic,profile_entry)
        END DO
      END DO
    END DO
    CLOSE(lun)

  END SUBROUTINE write_sb_spinup_forcing_debug

  ! ===================================================================================================== !
  !> Write a sparse pool and physical-state snapshot during forward pseudo-time.
  SUBROUTINE write_sb_spinup_snapshot_debug(tile, options, execution_index, represented_year, &
      & max_relative_change, max_total_c_absolute_change, max_total_c_relative_change, l_converged)

    USE mo_mpi,               ONLY: p_pe_work_only
    USE mo_jsb_process_class, ONLY: SB_, HYDRO_
    USE mo_lnd_bgcm_idx,      ONLY: ixC
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(HYDRO_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options
    INTEGER, INTENT(in) :: execution_index
    INTEGER, INTENT(in) :: represented_year
    REAL(lp), OPTIONAL, INTENT(in) :: max_relative_change
    REAL(lp), OPTIONAL, INTENT(in) :: max_total_c_absolute_change
    REAL(lp), OPTIONAL, INTENT(in) :: max_total_c_relative_change
    LOGICAL, OPTIONAL, INTENT(in) :: l_converged

    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(HYDRO_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':write_sb_spinup_snapshot_debug'
    REAL(lp), ALLOCATABLE :: pool_snapshot(:,:,:,:)
    INTEGER :: lun
    INTEGER :: ios
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    INTEGER :: ic
    INTEGER :: is
    INTEGER :: n_active_layers
    dsl4jsb_Real2D_onChunk :: num_sl_above_bedrock
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onChunk :: bulk_carbon_max
    dsl4jsb_Real3D_onChunk :: bulk_soil_carbon_sl
    dsl4jsb_Real3D_onChunk :: soil_litter_carbon_sl
    dsl4jsb_Real3D_onChunk :: volume_min_sl
    dsl4jsb_Real3D_onChunk :: qmax_org
    dsl4jsb_Real3D_onChunk :: surface_burial_frac
    dsl4jsb_Real3D_onChunk :: soil_lay_width_sl

    IF (p_pe_work_only /= 0 .OR. options%iblk /= 1) RETURN

    ics = options%ics
    ice = options%ice
    iblk = options%iblk

    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, num_sl_above_bedrock)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_lay_width_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, bulk_dens_corr_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, bulk_carbon_max)
    dsl4jsb_Get_var3D_onChunk(SB_, bulk_soil_carbon_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, soil_litter_carbon_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, volume_min_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, qmax_org)
    dsl4jsb_Get_var3D_onChunk(SB_, surface_burial_frac)
    CALL snapshot_sb_spinup_pools(tile, options, pool_snapshot)

    OPEN(newunit=lun, file=sb_spinup_debug_filename, status='unknown', position='append', &
      & action='write', iostat=ios)
    IF (ios /= 0) CALL finish(TRIM(routine), 'Could not open '//sb_spinup_debug_filename)

    WRITE(lun,'(A)') ''
    WRITE(lun,*) 'STATE_SNAPSHOT execution, represented_year = ', execution_index, represented_year
    IF (PRESENT(max_relative_change) .AND. PRESENT(max_total_c_absolute_change) .AND. &
        & PRESENT(max_total_c_relative_change) .AND. PRESENT(l_converged)) THEN
      WRITE(lun,*) 'RELAXATION_STATUS max_pool_relative_change, total_C_absolute_change_mol_m-2, ', &
        & 'total_C_relative_change, converged = ', max_relative_change, max_total_c_absolute_change, &
        & max_total_c_relative_change, l_converged
    END IF
    WRITE(lun,'(A)') 'STATE_COLUMNS execution represented_year ic layer '// &
      & 'C_DOM C_aDOM C_soluble_litter C_polymeric_litter C_woody_litter C_mycorrhiza '// &
      & 'C_microbial C_residue C_aRes C_total_concentration layer_width_m layer_C_stock_mol_m-2 '// &
      & 'bulk_soil_C soil_litter_C bulk_C_max bulk_density_corrected_kg_m-3 mineral_volume '// &
      & 'qmax_org_mol_C_m-3 surface_burial_fraction; C pools and bulk diagnostics are mol C m-3'
    DO ic = 1, SIZE(pool_snapshot, 3)
      n_active_layers = MIN(SIZE(pool_snapshot, 4), NINT(num_sl_above_bedrock(ic)))
      DO is = 1, n_active_layers
        WRITE(lun,*) 'STATE ', execution_index, represented_year, ic, is, &
          & pool_snapshot(:,ixC,ic,is), SUM(pool_snapshot(:,ixC,ic,is)), &
          & soil_lay_width_sl(ic,is), SUM(pool_snapshot(:,ixC,ic,is)) * soil_lay_width_sl(ic,is), &
          & bulk_soil_carbon_sl(ic,is), soil_litter_carbon_sl(ic,is), bulk_carbon_max(ic,is), &
          & bulk_dens_corr_sl(ic,is), volume_min_sl(ic,is), qmax_org(ic,is), surface_burial_frac(ic,is)
      END DO
    END DO
    CLOSE(lun)

  END SUBROUTINE write_sb_spinup_snapshot_debug

  ! ===================================================================================================== !
  !> Copy the complete soil-biogeochemistry pool matrix for annual convergence comparison.
  SUBROUTINE snapshot_sb_spinup_pools(tile, options, pool_snapshot)

    USE mo_lnd_bgcm_store,       ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_store_class, ONLY: SB_BGCM_POOL_ID

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    REAL(lp), ALLOCATABLE,      INTENT(out)   :: pool_snapshot(:,:,:,:)

    dsl4jsb_Def_mt2L3D :: sb_pool_mt

    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk

    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)

    ALLOCATE(pool_snapshot, SOURCE=sb_pool_mt)

  END SUBROUTINE snapshot_sb_spinup_pools

  ! ===================================================================================================== !
  !> Report the total and layer-resolved soil C change caused by one forward accelerator event.
  SUBROUTINE write_sb_spinup_debug(tile, options, execution_index, represented_years, pool_before, pool_after)

    USE mo_jsb_process_class, ONLY: HYDRO_
    USE mo_lnd_bgcm_idx,      ONLY: ixC
    dsl4jsb_Use_memory(HYDRO_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    INTEGER,                    INTENT(in)    :: execution_index
    INTEGER,                    INTENT(in)    :: represented_years
    REAL(lp),                   INTENT(in)    :: pool_before(:,:,:,:)
    REAL(lp),                   INTENT(in)    :: pool_after(:,:,:,:)

    dsl4jsb_Def_memory(HYDRO_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':write_sb_spinup_debug'
    CHARACTER(len=2048) :: debug_message
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    INTEGER :: ic
    INTEGER :: is
    INTEGER :: n_active_layers
    REAL(lp) :: layer_c_before
    REAL(lp) :: layer_c_after
    REAL(lp) :: total_c_before
    REAL(lp) :: total_c_after
    dsl4jsb_Real2D_onChunk :: num_sl_above_bedrock
    dsl4jsb_Real3D_onChunk :: soil_depth_sl

    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, num_sl_above_bedrock)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)

    IF (ANY(SHAPE(pool_before) /= SHAPE(pool_after))) THEN
      CALL finish(TRIM(routine), 'Before and after SB pool dimensions differ.')
    END IF

    total_c_before = 0.0_lp
    total_c_after = 0.0_lp
    DO ic = 1, SIZE(pool_before, 3)
      n_active_layers = MIN(SIZE(pool_before, 4), NINT(num_sl_above_bedrock(ic)))
      DO is = 1, n_active_layers
        total_c_before = total_c_before + &
          & SUM(pool_before(:,ixC,ic,is)) * soil_depth_sl(ic,is)
        total_c_after = total_c_after + &
          & SUM(pool_after(:,ixC,ic,is)) * soil_depth_sl(ic,is)
      END DO
    END DO

    WRITE(debug_message,*) 'DEBUG forward SB acceleration carbon totals execution/represented_years = ', &
      & execution_index, represented_years, '; before/after/delta [mol C m-2] = ', &
      & total_c_before, total_c_after, total_c_after - total_c_before
    CALL message(TRIM(routine), TRIM(debug_message))

    DO ic = 1, SIZE(pool_before, 3)
      n_active_layers = MIN(SIZE(pool_before, 4), NINT(num_sl_above_bedrock(ic)))
      DO is = 1, n_active_layers
        layer_c_before = SUM(pool_before(:,ixC,ic,is))
        layer_c_after = SUM(pool_after(:,ixC,ic,is))
        WRITE(debug_message,*) 'DEBUG forward SB acceleration profile execution/ic/is = ', &
          & execution_index, ic, is, '; total C before/after/delta [mol C m-3] = ', &
          & layer_c_before, layer_c_after, layer_c_after - layer_c_before, &
          & '; pool C before [mol C m-3] = ', pool_before(:,ixC,ic,is), &
          & '; pool C after [mol C m-3] = ', pool_after(:,ixC,ic,is)
        CALL message(TRIM(routine), TRIM(debug_message))
      END DO
    END DO

  END SUBROUTINE write_sb_spinup_debug

  ! ===================================================================================================== !
  !> Test annual layer-, pool-, and element-wise changes and update the convergence state.
  SUBROUTINE check_sb_spinup_convergence(tile, options, previous_pool, consecutive_converged_years, &
      & max_relative_change, max_total_c_absolute_change, max_total_c_relative_change, l_converged)

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    USE mo_jsb_math_constants,   ONLY: eps8
    USE mo_jsb_process_class,    ONLY: SB_, HYDRO_
    USE mo_lnd_bgcm_idx,         ONLY: ixC
    USE mo_lnd_bgcm_store,       ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_store_class, ONLY: SB_BGCM_POOL_ID
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(HYDRO_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    REAL(lp), ALLOCATABLE,      INTENT(inout) :: previous_pool(:,:,:,:)
    INTEGER,                    INTENT(inout) :: consecutive_converged_years
    REAL(lp),                   INTENT(out)   :: max_relative_change
    REAL(lp),                   INTENT(out)   :: max_total_c_absolute_change
    REAL(lp),                   INTENT(out)   :: max_total_c_relative_change
    LOGICAL,                    INTENT(out)   :: l_converged

    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_mt2L3D :: sb_pool_mt

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store
    CHARACTER(len=*), PARAMETER :: routine = modname//':check_sb_spinup_convergence'
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    REAL(lp), ALLOCATABLE :: current_total_c(:)
    REAL(lp), ALLOCATABLE :: previous_total_c(:)
    dsl4jsb_Real3D_onChunk :: soil_depth_sl

    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    model => Get_model(tile%owner_model_id)
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)

    IF (ANY(SHAPE(previous_pool) /= SHAPE(sb_pool_mt))) THEN
      CALL finish(TRIM(routine), 'Stored and current SB pool dimensions differ.')
    END IF
    IF (.NOT. ALL(ieee_is_finite(sb_pool_mt))) THEN
      CALL finish(TRIM(routine), 'Non-finite SB pool encountered during forward relaxation.')
    END IF

    max_relative_change = MAXVAL(ABS(sb_pool_mt - previous_pool) / &
      & MAX(ABS(sb_pool_mt), ABS(previous_pool), eps8))

    ALLOCATE(current_total_c(SIZE(sb_pool_mt, 3)), previous_total_c(SIZE(sb_pool_mt, 3)))
    current_total_c = SUM(SUM(sb_pool_mt(:, ixC, :, :), DIM=1) * soil_depth_sl, DIM=2)
    previous_total_c = SUM(SUM(previous_pool(:, ixC, :, :), DIM=1) * soil_depth_sl, DIM=2)
    max_total_c_absolute_change = MAXVAL(ABS(current_total_c - previous_total_c))
    max_total_c_relative_change = MAXVAL(ABS(current_total_c - previous_total_c) / &
      & MAX(ABS(current_total_c), ABS(previous_total_c), eps8))

    IF (max_relative_change <= dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_convergence_rel_tol) THEN
      consecutive_converged_years = consecutive_converged_years + 1
    ELSE
      consecutive_converged_years = 0
    END IF
    l_converged = consecutive_converged_years >= &
      & dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_convergence_years

    previous_pool = sb_pool_mt

  END SUBROUTINE check_sb_spinup_convergence

  ! ===================================================================================================== !
  !> Commit the relaxed SB state, restore live forcing, and refresh derived soil properties.
  SUBROUTINE finalize_sb_spinup_accelerator(tile, options, live_forcing)

    USE mo_jsb_process_class, ONLY: VEG_, HYDRO_, SSE_, A2L_
    USE mo_lnd_bgcm_idx,       ONLY: ix_fine_root, ixC
    USE mo_lnd_bgcm_store,       ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_store_class, ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID, VEG_BGCM_SEED_BED_LITTERFALL_ID
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(A2L_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    TYPE(t_sb_spinup_live_forcing), INTENT(in) :: live_forcing

    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt1L2D :: seed_bed_litterfall_mt

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    dsl4jsb_Real3D_onChunk :: frac_wtr_transp_down_sl
    dsl4jsb_Real3D_onChunk :: frac_w_lat_loss_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl
    dsl4jsb_Real3D_onChunk :: ice_soil_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_sl
    dsl4jsb_Real3D_onChunk :: t_soil_sl
    dsl4jsb_Real2D_onChunk :: net_biosphere_production
    dsl4jsb_Real2D_onChunk :: net_biosphere_production_c13
    dsl4jsb_Real2D_onChunk :: net_biosphere_production_c14
    dsl4jsb_Real2D_onChunk :: biological_n_fixation
    dsl4jsb_Real2D_onChunk :: f_n_demand
    dsl4jsb_Real2D_onChunk :: f_p_demand
    dsl4jsb_Real2D_onChunk :: growth_p_limit_based_on_n_mavg
    dsl4jsb_Real2D_onChunk :: nhx_deposition
    dsl4jsb_Real2D_onChunk :: noy_deposition
    dsl4jsb_Real2D_onChunk :: nhx_n15_deposition
    dsl4jsb_Real2D_onChunk :: noy_n15_deposition
    dsl4jsb_Real2D_onChunk :: p_deposition

    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    model => Get_model(tile%owner_model_id)
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, frac_wtr_transp_down_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, frac_w_lat_loss_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pot_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)
    dsl4jsb_Get_var2D_onChunk(VEG_, net_biosphere_production)
    dsl4jsb_Get_var2D_onChunk(VEG_, net_biosphere_production_c13)
    dsl4jsb_Get_var2D_onChunk(VEG_, net_biosphere_production_c14)
    dsl4jsb_Get_var2D_onChunk(VEG_, biological_n_fixation)
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

    veg_litterfall_mt = live_forcing%litterfall
    seed_bed_litterfall_mt = live_forcing%seed_litterfall
    veg_pool_mt(ix_fine_root, ixC, :) = live_forcing%fine_root_c
    f_n_demand = live_forcing%f_n_demand
    f_p_demand = live_forcing%f_p_demand
    growth_p_limit_based_on_n_mavg = live_forcing%growth_p_limit
    nhx_deposition = live_forcing%nhx_deposition
    noy_deposition = live_forcing%noy_deposition
    nhx_n15_deposition = live_forcing%nhx_n15_deposition
    noy_n15_deposition = live_forcing%noy_n15_deposition
    p_deposition = live_forcing%p_deposition
    root_fraction_sl = live_forcing%root_fraction
    frac_wtr_transp_down_sl = live_forcing%wtr_down
    frac_w_lat_loss_sl = live_forcing%wtr_lat_loss
    t_soil_sl = live_forcing%t_soil
    wtr_soil_sl = live_forcing%wtr_soil
    ice_soil_sl = live_forcing%ice_soil
    wtr_soil_pot_sl = live_forcing%wtr_pot
    net_biosphere_production = live_forcing%net_biosphere_production
    net_biosphere_production_c13 = live_forcing%net_biosphere_production_c13
    net_biosphere_production_c14 = live_forcing%net_biosphere_production_c14
    biological_n_fixation = live_forcing%biological_n_fixation

  END SUBROUTINE finalize_sb_spinup_accelerator

#endif
END MODULE mo_sb_spinup_accelerator

!> QUINCY IQJ soil-model helper routines
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### helper routines for the IQJ soil-forcing workflow
!>
MODULE mo_sb_spinup_qsoil_forcing
#ifndef __NO_QUINCY__

  USE mo_kind,                  ONLY: lp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_jsb_math_constants,    ONLY: zero
  USE mo_util,                  ONLY: int2string

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store, get_bgcm_idx
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID, VEG_BGCM_SEED_BED_LITTERFALL_ID
#ifndef  __QUINCY_STANDALONE__
  USE mo_jsb_io_netcdf,           ONLY: t_input_file
  USE mo_read_interface,          ONLY: read_2D_extdim, read_3D_extdim, on_cells
#endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_read_soil_forcing

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_spinup_qsoil_forcing'

#ifndef  __QUINCY_STANDALONE__
  TYPE t_soil_forcing_cache
    TYPE(t_input_file) :: soilphysics_input_file
    TYPE(t_input_file) :: veglitterfall_input_file
    LOGICAL :: soilphysics_open = .FALSE.
    LOGICAL :: veglitterfall_open = .FALSE.
    LOGICAL :: arrays_initialized = .FALSE.
    LOGICAL :: preload_loaded = .FALSE.
    LOGICAL :: isotope_forcing_initialized = .FALSE.
    INTEGER :: nproma = 0
    INTEGER :: nsoil = 0
    INTEGER :: nblks = 0
    INTEGER :: npft = 0
    INTEGER :: veg_litter_nr_parts = 0
    INTEGER :: veg_litter_nr_elems = 0
    INTEGER :: seed_bed_nr_elems = 0
    INTEGER :: preload_nsteps = 0
    INTEGER :: cached_year = -HUGE(0)
    INTEGER :: cached_yearday = -HUGE(0)
    INTEGER :: cached_timestep = -HUGE(0)
    CHARACTER(len=:), ALLOCATABLE :: soilphysics_filename
    CHARACTER(len=:), ALLOCATABLE :: veglitterfall_filename
    REAL(lp), ALLOCATABLE :: root_fraction_sl_preload(:,:,:,:,:)      ! nproma, nsoil, nblks, npft, ntime
    REAL(lp), ALLOCATABLE :: demand_uptake_n_preload(:,:,:,:)         ! nproma, nblks, npft, ntime
    REAL(lp), ALLOCATABLE :: demand_uptake_p_preload(:,:,:,:)         ! nproma, nblks, npft, ntime
    REAL(lp), ALLOCATABLE :: frac_wtr_transp_down_sl_preload(:,:,:,:) ! nproma, nsoil, nblks, ntime
    REAL(lp), ALLOCATABLE :: frac_w_lat_loss_sl_preload(:,:,:,:)      ! nproma, nsoil, nblks, ntime
    REAL(lp), ALLOCATABLE :: wtr_soil_sl_preload(:,:,:,:)             ! nproma, nsoil, nblks, ntime
    REAL(lp), ALLOCATABLE :: vol_porosity_sl_preload(:,:,:,:)         ! nproma, nsoil, nblks, ntime
    REAL(lp), ALLOCATABLE :: t_soil_sl_preload(:,:,:,:)               ! nproma, nsoil, nblks, ntime
    REAL(lp), ALLOCATABLE :: nhx_deposition_preload(:,:,:)            ! nproma, nblks, ntime
    REAL(lp), ALLOCATABLE :: noy_deposition_preload(:,:,:)            ! nproma, nblks, ntime
    REAL(lp), ALLOCATABLE :: nhx_n15_deposition_preload(:,:,:)        ! nproma, nblks, ntime
    REAL(lp), ALLOCATABLE :: noy_n15_deposition_preload(:,:,:)        ! nproma, nblks, ntime
    REAL(lp), ALLOCATABLE :: p_deposition_preload(:,:,:)              ! nproma, nblks, ntime
    REAL(lp), ALLOCATABLE :: veg_pool_fine_root_c_preload(:,:,:,:)    ! nproma, nblks, npft, ntime
    REAL(lp), ALLOCATABLE :: veg_litterfall_preload(:,:,:,:,:,:)      ! part, elem, nproma, nblks, npft, ntime
    REAL(lp), ALLOCATABLE :: seed_bed_litterfall_preload(:,:,:,:,:)   ! elem, nproma, nblks, npft, ntime
    REAL(lp), ALLOCATABLE :: scratch_3d(:,:,:)                        ! nproma, nblks, 1
    REAL(lp), ALLOCATABLE :: scratch_4d(:,:,:,:)                      ! nproma, nsoil, nblks, 1
  END TYPE t_soil_forcing_cache

  TYPE(t_soil_forcing_cache), SAVE :: soil_forcing_cache
#endif

CONTAINS

#ifndef __QUINCY_STANDALONE__
  SUBROUTINE log_soil_forcing_timing(caller, label, dt, detail)

    CHARACTER(len=*), INTENT(in) :: caller, label
    REAL(lp),         INTENT(in) :: dt
    CHARACTER(len=*), INTENT(in), OPTIONAL :: detail

    CHARACTER(len=256) :: line

    IF (.NOT. debug_on()) RETURN

    IF (PRESENT(detail)) THEN
      WRITE(line,'(A,": ",F10.6," s (",A,")")') TRIM(label), dt, TRIM(detail)
    ELSE
      WRITE(line,'(A,": ",F10.6," s")') TRIM(label), dt
    END IF
    CALL message(TRIM(caller), TRIM(line))

  END SUBROUTINE log_soil_forcing_timing
#endif

  ! ======================================================================================================= !
  !> read soil forcing from input files (ONLY for IQ)
  !>
#ifndef  __QUINCY_STANDALONE__
  SUBROUTINE update_read_soil_forcing(tile, options)

    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_process_class,       ONLY: SB_, HYDRO_, SSE_, VEG_, A2L_

    USE mo_lnd_bgcm_idx
    USE mo_jsb_time,                ONLY: t_datetime, get_date_components, get_year_day

    dsl4jsb_Use_processes SB_
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_config(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER :: model            !< the model
    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store       !< the bgcm store of this tile
    REAL(lp)                        :: dtime, year_instant
    INTEGER                         :: iblk             !< Dimensions of the horizontal grid
    INTEGER                         :: pft_idx, npft
    INTEGER                         :: year, month, day, hour, minute, second, yearday, forcing_yearday, &
                                    &  ith_timestep, cache_itimestep     !< forcing record index in file / cache
    CHARACTER(len=:), ALLOCATABLE   :: soilphysics_forcing_file, veglitterfall_forcing_file
    LOGICAL                         :: cache_hit, trace_timing, skip_forcing_isotopes_after_init
    REAL(lp)                        :: t_total_start, t_stage_start, t_stage_end

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_read_soil_forcing'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_config(SB_)
    ! A2L_
    dsl4jsb_Real2D_onDomain      :: nhx_deposition
    dsl4jsb_Real2D_onDomain      :: noy_deposition
    dsl4jsb_Real2D_onDomain      :: nhx_n15_deposition
    dsl4jsb_Real2D_onDomain      :: noy_n15_deposition
    dsl4jsb_Real2D_onDomain      :: p_deposition
    ! VEG_ 2D
    dsl4jsb_Real2D_onDomain      :: demand_uptake_n_tuptake_mavg
    dsl4jsb_Real2D_onDomain      :: demand_uptake_p_tuptake_mavg
    ! VEG_ 3D
    dsl4jsb_Real3D_onDomain      :: root_fraction_sl
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onDomain      :: frac_wtr_transp_down_sl
    dsl4jsb_Real3D_onDomain      :: frac_w_lat_loss_sl
    dsl4jsb_Real3D_onDomain      :: wtr_soil_sl
    dsl4jsb_Real3D_onDomain      :: vol_porosity_sl
    ! SSE_ 3D
    dsl4jsb_Real3D_onDomain      :: t_soil_sl
    ! ----------------------------------------------------------------------------------------------------- !
    trace_timing = debug_on() .AND. (TRIM(tile%name) == 'box' .OR. TRIM(tile%name) == 'pft01')
    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    IF (trace_timing) CALL cpu_time(t_total_start)
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    iblk = options%iblk
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SB_)
    skip_forcing_isotopes_after_init = dsl4jsb_Config(SB_)%skip_forcing_isotopes_after_init
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onDomain(A2L_,  nhx_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_,  noy_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_,  nhx_n15_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_,  noy_n15_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_,  p_deposition)
    dsl4jsb_Get_var2D_onDomain(VEG_,  demand_uptake_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_,  demand_uptake_p_tuptake_mavg)
    dsl4jsb_Get_var3D_onDomain(VEG_,  root_fraction_sl)
    ! HYDRO_ 3D
    dsl4jsb_Get_var3D_onDomain(HYDRO_, frac_wtr_transp_down_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, frac_w_lat_loss_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)
    ! SSE_ 3D
    dsl4jsb_Get_var3D_onDomain(SSE_, t_soil_sl)

    !> get the timestep information of current timestep
    dtime   = options%dtime
    CALL get_date_components(options%current_datetime, year, month, day, hour, minute, second)
    year_instant = get_year_day(options%current_datetime)
    yearday = INT(year_instant)
    forcing_yearday = MODULO(yearday - 1, 365) + 1
    CALL get_soil_forcing_info(TRIM(dsl4jsb_Config(SB_)%soilphysics_forcing_filename), &
      & forcing_yearday, hour, minute, second, dtime, soilphysics_forcing_file, ith_timestep)
    CALL get_soil_forcing_info(TRIM(dsl4jsb_Config(SB_)%veglitterfall_forcing_filename), &
      & forcing_yearday, hour, minute, second, dtime, veglitterfall_forcing_file, ith_timestep)
    IF (debug_on()) CALL message(TRIM(routine), 'Reading the forcing record '//int2string(ith_timestep)// &
      & ' at forcing day '//int2string(forcing_yearday)//' (model day '//int2string(yearday)//') from ' // &
      & TRIM(soilphysics_forcing_file)//' ...')
    bgcm_store => tile%bgcm_store
    npft = model%config%npft
    pft_idx = get_pft_tile_index(tile%name)
    IF (trace_timing) CALL cpu_time(t_stage_start)
    CALL ensure_soil_forcing_cache_loaded(model, bgcm_store, root_fraction_sl, soilphysics_forcing_file, &
      & veglitterfall_forcing_file, year, forcing_yearday, ith_timestep, dtime, npft, &
      & skip_forcing_isotopes_after_init, routine, cache_hit)
    IF (trace_timing) THEN
      CALL cpu_time(t_stage_end)
      IF (cache_hit) THEN
        CALL log_soil_forcing_timing(routine, 'cache_ready tile='//TRIM(tile%name), t_stage_end - t_stage_start, 'cache-hit')
      ELSE
        CALL log_soil_forcing_timing(routine, 'cache_ready tile='//TRIM(tile%name), t_stage_end - t_stage_start, 'cache-reload')
      END IF
    END IF
    cache_itimestep = ith_timestep - soil_forcing_cache%cached_timestep + 1
    IF (cache_itimestep < 1 .OR. cache_itimestep > soil_forcing_cache%preload_nsteps) THEN
      CALL finish(TRIM(routine), 'Invalid soil forcing cache index '//int2string(cache_itimestep)// &
        & ' for forcing timestep '//int2string(ith_timestep)//' with cache start '// &
        & int2string(soil_forcing_cache%cached_timestep)//' and preload_nsteps '// &
        & int2string(soil_forcing_cache%preload_nsteps))
    END IF

    IF (pft_idx > 0) THEN
      IF (trace_timing) CALL cpu_time(t_stage_start)
      root_fraction_sl(:,:,:) = soil_forcing_cache%root_fraction_sl_preload(:,:,:,pft_idx,cache_itimestep)
      demand_uptake_n_tuptake_mavg(:,:) = soil_forcing_cache%demand_uptake_n_preload(:,:,pft_idx,cache_itimestep)
      demand_uptake_p_tuptake_mavg(:,:) = soil_forcing_cache%demand_uptake_p_preload(:,:,pft_idx,cache_itimestep)
      CALL apply_litterfall_cache_to_tile(bgcm_store, tile%name, pft_idx, cache_itimestep, routine)
      IF (trace_timing) THEN
        CALL cpu_time(t_stage_end)
        CALL log_soil_forcing_timing(routine, 'apply_pft tile='//TRIM(tile%name), t_stage_end - t_stage_start)
      END IF
    ELSE
      IF (trace_timing) CALL cpu_time(t_stage_start)
      frac_wtr_transp_down_sl(:,:,:) = soil_forcing_cache%frac_wtr_transp_down_sl_preload(:,:,:,cache_itimestep)
      frac_w_lat_loss_sl(:,:,:)      = soil_forcing_cache%frac_w_lat_loss_sl_preload(:,:,:,cache_itimestep)
      wtr_soil_sl(:,:,:)             = soil_forcing_cache%wtr_soil_sl_preload(:,:,:,cache_itimestep)
      vol_porosity_sl(:,:,:)         = soil_forcing_cache%vol_porosity_sl_preload(:,:,:,cache_itimestep)
      t_soil_sl(:,:,:)               = soil_forcing_cache%t_soil_sl_preload(:,:,:,cache_itimestep)
      nhx_deposition(:,:)            = soil_forcing_cache%nhx_deposition_preload(:,:,cache_itimestep)
      noy_deposition(:,:)            = soil_forcing_cache%noy_deposition_preload(:,:,cache_itimestep)
      nhx_n15_deposition(:,:)        = soil_forcing_cache%nhx_n15_deposition_preload(:,:,cache_itimestep)
      noy_n15_deposition(:,:)        = soil_forcing_cache%noy_n15_deposition_preload(:,:,cache_itimestep)
      p_deposition(:,:)              = soil_forcing_cache%p_deposition_preload(:,:,cache_itimestep)
      IF (trace_timing) THEN
        CALL cpu_time(t_stage_end)
        CALL log_soil_forcing_timing(routine, 'apply_box tile='//TRIM(tile%name), t_stage_end - t_stage_start)
      END IF
    END IF

    IF (trace_timing) THEN
      CALL cpu_time(t_stage_end)
      CALL log_soil_forcing_timing(routine, 'total tile='//TRIM(tile%name), t_stage_end - t_total_start)
    END IF
    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_read_soil_forcing

  ! ======================================================================================================= !
  !> Return the number of records to preload for the current forcing file.
  !>
  INTEGER FUNCTION get_soil_forcing_preload_nsteps(dtime)

    REAL(lp), INTENT(in) :: dtime

    IF (NINT(dtime) == 86400) THEN
      get_soil_forcing_preload_nsteps = 1
    ELSE
      get_soil_forcing_preload_nsteps = 86400 / INT(dtime)
    END IF

  END FUNCTION get_soil_forcing_preload_nsteps

  ! ======================================================================================================= !
  ! ======================================================================================================= !
  !> Return the pft index encoded in a tile name like pft01.
  !>
  INTEGER FUNCTION get_pft_tile_index(tile_name)

    CHARACTER(len=*), INTENT(in) :: tile_name
    INTEGER :: ios

    get_pft_tile_index = 0
    IF (LEN_TRIM(tile_name) /= 5) RETURN
    IF (tile_name(1:3) /= 'pft') RETURN

    READ(tile_name(4:5), '(I2)', IOSTAT=ios) get_pft_tile_index
    IF (ios /= 0) get_pft_tile_index = 0

  END FUNCTION get_pft_tile_index

  ! ======================================================================================================= !
  !> Return the canonical pft tile name.
  !>
  FUNCTION get_pft_tile_name(ipft) RESULT(tile_name)

    INTEGER, INTENT(in) :: ipft
    CHARACTER(len=5)    :: tile_name

    WRITE(tile_name, '("pft",I2.2)') ipft

  END FUNCTION get_pft_tile_name

  ! ======================================================================================================= !
  !> Return the litterfall compartment name used in the forcing variable names.
  !>
  FUNCTION get_veg_litter_part_name(i_part) RESULT(part_name)

    INTEGER, INTENT(in)    :: i_part
    CHARACTER(len=:), ALLOCATABLE :: part_name

    SELECT CASE(i_part)
    CASE (1)
      part_name = 'leaf'
    CASE (2)
      part_name = 'fine_root'
    CASE (3)
      part_name = 'coarse_root'
    CASE (4)
      part_name = 'sap_wood'
    CASE (5)
      part_name = 'heart_wood'
    CASE (6)
      part_name = 'labile'
    CASE (7)
      part_name = 'reserve'
    CASE (8)
      part_name = 'fruit'
    CASE DEFAULT
      CALL finish(modname//':get_veg_litter_part_name', 'Unsupported litterfall part index '//int2string(i_part))
    END SELECT

  END FUNCTION get_veg_litter_part_name

  ! ======================================================================================================= !
  !> Return the litterfall element name used in the forcing variable names.
  !>
  FUNCTION get_veg_litter_elem_name(i_elem) RESULT(elem_name)

    INTEGER, INTENT(in)    :: i_elem
    CHARACTER(len=:), ALLOCATABLE :: elem_name

    SELECT CASE(i_elem)
    CASE (1)
      elem_name = 'C'
    CASE (2)
      elem_name = 'N'
    CASE (3)
      elem_name = 'P'
    CASE (4)
      elem_name = 'C13'
    CASE (5)
      elem_name = 'C14'
    CASE (6)
      elem_name = 'N15'
    CASE DEFAULT
      CALL finish(modname//':get_veg_litter_elem_name', 'Unsupported litterfall element index '//int2string(i_elem))
    END SELECT

  END FUNCTION get_veg_litter_elem_name

  ! ======================================================================================================= !
  !> Keep the forcing files open and cache the currently loaded timestep.
  !>
  SUBROUTINE ensure_soil_forcing_cache_loaded(model, bgcm_store, root_fraction_sl, soilphysics_filename, &
    & veglitterfall_filename, year, yearday, itimestep, dtime, npft, skip_forcing_isotopes_after_init, caller, cache_hit_out)

    USE mo_jsb_model_class, ONLY: t_jsb_model

    TYPE(t_jsb_model),       POINTER, INTENT(in) :: model
    TYPE(t_lnd_bgcm_store),  POINTER, INTENT(in) :: bgcm_store
    REAL(lp),                        INTENT(in) :: root_fraction_sl(:,:,:)
    CHARACTER(len=*),                INTENT(in) :: soilphysics_filename, veglitterfall_filename, caller
    INTEGER,                         INTENT(in) :: year, yearday, itimestep, npft
    REAL(lp),                        INTENT(in) :: dtime
    LOGICAL,                         INTENT(in) :: skip_forcing_isotopes_after_init
    LOGICAL, OPTIONAL,               INTENT(out) :: cache_hit_out

    LOGICAL :: cache_hit, skip_isotope_reload
    INTEGER :: preload_nsteps, cache_start_timestep
    REAL(lp) :: t_stage_start, t_stage_end

    preload_nsteps = get_soil_forcing_preload_nsteps(dtime)
    IF (preload_nsteps == 1) THEN
      cache_start_timestep = itimestep
    ELSE
      cache_start_timestep = 1
    END IF
    !$OMP CRITICAL (soil_forcing_netcdf_io)
    IF (debug_on()) CALL cpu_time(t_stage_start)
    CALL ensure_soil_forcing_cache_shape(bgcm_store, root_fraction_sl, npft, preload_nsteps, caller)
    CALL open_cached_soil_forcing_file(model, soil_forcing_cache%soilphysics_input_file, &
      & soil_forcing_cache%soilphysics_open, soil_forcing_cache%soilphysics_filename, soilphysics_filename, caller)
    CALL open_cached_soil_forcing_file(model, soil_forcing_cache%veglitterfall_input_file, &
      & soil_forcing_cache%veglitterfall_open, soil_forcing_cache%veglitterfall_filename, veglitterfall_filename, caller)
    IF (debug_on()) THEN
      CALL cpu_time(t_stage_end)
      CALL log_soil_forcing_timing(caller, 'cache_prepare', t_stage_end - t_stage_start)
    END IF

    cache_hit = .FALSE.
    IF (ALLOCATED(soil_forcing_cache%soilphysics_filename) .AND. ALLOCATED(soil_forcing_cache%veglitterfall_filename)) THEN
      cache_hit = soil_forcing_cache%preload_loaded .AND. soil_forcing_cache%preload_nsteps == preload_nsteps .AND. &
        & TRIM(soil_forcing_cache%soilphysics_filename) == TRIM(soilphysics_filename) .AND. &
        & TRIM(soil_forcing_cache%veglitterfall_filename) == TRIM(veglitterfall_filename) .AND. &
        & soil_forcing_cache%cached_timestep == cache_start_timestep
    END IF

    IF (.NOT. cache_hit) THEN
      IF (debug_on()) CALL cpu_time(t_stage_start)
      skip_isotope_reload = skip_forcing_isotopes_after_init .AND. &
        & soil_forcing_cache%isotope_forcing_initialized

      CALL load_soilphysics_cache(bgcm_store, preload_nsteps, cache_start_timestep, skip_isotope_reload, caller)
      IF (debug_on()) THEN
        CALL cpu_time(t_stage_end)
        CALL log_soil_forcing_timing(caller, 'load_soilphysics_cache', t_stage_end - t_stage_start)
      END IF
      IF (debug_on()) CALL cpu_time(t_stage_start)
      CALL load_litterfall_cache(bgcm_store, preload_nsteps, cache_start_timestep, skip_isotope_reload, caller)
      IF (debug_on()) THEN
        CALL cpu_time(t_stage_end)
        CALL log_soil_forcing_timing(caller, 'load_litterfall_cache', t_stage_end - t_stage_start)
      END IF
      soil_forcing_cache%preload_loaded = .TRUE.
      IF (skip_forcing_isotopes_after_init) soil_forcing_cache%isotope_forcing_initialized = .TRUE.
      soil_forcing_cache%cached_year = year
      soil_forcing_cache%cached_yearday = yearday
      soil_forcing_cache%cached_timestep = cache_start_timestep
    END IF
    !$OMP END CRITICAL (soil_forcing_netcdf_io)
    IF (PRESENT(cache_hit_out)) cache_hit_out = cache_hit

  END SUBROUTINE ensure_soil_forcing_cache_loaded

  ! ======================================================================================================= !
  !> Ensure that the cache arrays match the current model dimensions.
  !>
  SUBROUTINE ensure_soil_forcing_cache_shape(bgcm_store, root_fraction_sl, npft, preload_nsteps, caller)

    TYPE(t_lnd_bgcm_store), POINTER, INTENT(in) :: bgcm_store
    REAL(lp),                      INTENT(in) :: root_fraction_sl(:,:,:)
    INTEGER,                       INTENT(in) :: npft
    INTEGER,                       INTENT(in) :: preload_nsteps
    CHARACTER(len=*),              INTENT(in) :: caller

    INTEGER :: nproma, nsoil, nblks
    INTEGER :: veg_litter_idx, seed_bed_idx
    INTEGER :: veg_litter_store_idx, seed_bed_store_idx
    LOGICAL :: needs_realloc
    CHARACTER(len=5) :: tile_name

    nproma = SIZE(root_fraction_sl, 1)
    nsoil  = SIZE(root_fraction_sl, 2)
    nblks  = SIZE(root_fraction_sl, 3)
    tile_name = get_pft_tile_name(1)

    veg_litter_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_LITTERFALL_ID, tile_name, caller)
    seed_bed_idx   = get_bgcm_idx(bgcm_store, VEG_BGCM_SEED_BED_LITTERFALL_ID, tile_name, caller)
    veg_litter_store_idx = bgcm_store%idx_in_store(veg_litter_idx)
    seed_bed_store_idx   = bgcm_store%idx_in_store(seed_bed_idx)

    needs_realloc = .NOT. soil_forcing_cache%arrays_initialized
    IF (.NOT. needs_realloc) THEN
      needs_realloc = soil_forcing_cache%nproma /= nproma .OR. soil_forcing_cache%nsoil /= nsoil .OR. &
        & soil_forcing_cache%nblks /= nblks .OR. soil_forcing_cache%npft /= npft .OR. &
        & soil_forcing_cache%preload_nsteps /= preload_nsteps .OR. &
        & soil_forcing_cache%veg_litter_nr_parts /= bgcm_store%store_2l_2d_bgcms(veg_litter_store_idx)%nr_of_parts .OR. &
        & soil_forcing_cache%veg_litter_nr_elems /= bgcm_store%store_2l_2d_bgcms(veg_litter_store_idx)%nr_of_elems .OR. &
        & soil_forcing_cache%seed_bed_nr_elems /= bgcm_store%store_1l_2d_bgcms(seed_bed_store_idx)%nr_of_elems
    END IF

    IF (.NOT. needs_realloc) RETURN

    CALL reset_soil_forcing_cache_arrays()

    soil_forcing_cache%nproma = nproma
    soil_forcing_cache%nsoil  = nsoil
    soil_forcing_cache%nblks  = nblks
    soil_forcing_cache%npft   = npft
    soil_forcing_cache%preload_nsteps = preload_nsteps
    soil_forcing_cache%veg_litter_nr_parts = bgcm_store%store_2l_2d_bgcms(veg_litter_store_idx)%nr_of_parts
    soil_forcing_cache%veg_litter_nr_elems = bgcm_store%store_2l_2d_bgcms(veg_litter_store_idx)%nr_of_elems
    soil_forcing_cache%seed_bed_nr_elems   = bgcm_store%store_1l_2d_bgcms(seed_bed_store_idx)%nr_of_elems
    ALLOCATE(soil_forcing_cache%root_fraction_sl_preload(nproma, nsoil, nblks, npft, preload_nsteps))
    ALLOCATE(soil_forcing_cache%demand_uptake_n_preload(nproma, nblks, npft, preload_nsteps))
    ALLOCATE(soil_forcing_cache%demand_uptake_p_preload(nproma, nblks, npft, preload_nsteps))
    ALLOCATE(soil_forcing_cache%frac_wtr_transp_down_sl_preload(nproma, nsoil, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%frac_w_lat_loss_sl_preload(nproma, nsoil, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%wtr_soil_sl_preload(nproma, nsoil, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%vol_porosity_sl_preload(nproma, nsoil, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%t_soil_sl_preload(nproma, nsoil, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%nhx_deposition_preload(nproma, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%noy_deposition_preload(nproma, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%nhx_n15_deposition_preload(nproma, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%noy_n15_deposition_preload(nproma, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%p_deposition_preload(nproma, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%veg_pool_fine_root_c_preload(nproma, nblks, npft, preload_nsteps))
    ALLOCATE(soil_forcing_cache%veg_litterfall_preload(soil_forcing_cache%veg_litter_nr_parts, soil_forcing_cache%veg_litter_nr_elems, &
      & nproma, nblks, npft, preload_nsteps))
    ALLOCATE(soil_forcing_cache%seed_bed_litterfall_preload(soil_forcing_cache%seed_bed_nr_elems, nproma, nblks, npft, preload_nsteps))
    ALLOCATE(soil_forcing_cache%scratch_3d(nproma, nblks, 1))
    ALLOCATE(soil_forcing_cache%scratch_4d(nproma, nsoil, nblks, 1))

    soil_forcing_cache%root_fraction_sl_preload = zero
    soil_forcing_cache%demand_uptake_n_preload = zero
    soil_forcing_cache%demand_uptake_p_preload = zero
    soil_forcing_cache%frac_wtr_transp_down_sl_preload = zero
    soil_forcing_cache%frac_w_lat_loss_sl_preload = zero
    soil_forcing_cache%wtr_soil_sl_preload = zero
    soil_forcing_cache%vol_porosity_sl_preload = zero
    soil_forcing_cache%t_soil_sl_preload = zero
    soil_forcing_cache%nhx_deposition_preload = zero
    soil_forcing_cache%noy_deposition_preload = zero
    soil_forcing_cache%nhx_n15_deposition_preload = zero
    soil_forcing_cache%noy_n15_deposition_preload = zero
    soil_forcing_cache%p_deposition_preload = zero
    soil_forcing_cache%veg_pool_fine_root_c_preload = zero
    soil_forcing_cache%veg_litterfall_preload = zero
    soil_forcing_cache%seed_bed_litterfall_preload = zero
    soil_forcing_cache%scratch_3d = zero
    soil_forcing_cache%scratch_4d = zero

    soil_forcing_cache%arrays_initialized = .TRUE.
    soil_forcing_cache%preload_loaded = .FALSE.
    soil_forcing_cache%cached_year = -HUGE(0)
    soil_forcing_cache%cached_yearday = -HUGE(0)
    soil_forcing_cache%cached_timestep = -HUGE(0)
  END SUBROUTINE ensure_soil_forcing_cache_shape

  ! ======================================================================================================= !
  !> Deallocate the cache arrays while keeping any open files untouched.
  !>
  SUBROUTINE reset_soil_forcing_cache_arrays()

    IF (ALLOCATED(soil_forcing_cache%root_fraction_sl_preload)) DEALLOCATE(soil_forcing_cache%root_fraction_sl_preload)
    IF (ALLOCATED(soil_forcing_cache%demand_uptake_n_preload)) DEALLOCATE(soil_forcing_cache%demand_uptake_n_preload)
    IF (ALLOCATED(soil_forcing_cache%demand_uptake_p_preload)) DEALLOCATE(soil_forcing_cache%demand_uptake_p_preload)
    IF (ALLOCATED(soil_forcing_cache%frac_wtr_transp_down_sl_preload)) DEALLOCATE(soil_forcing_cache%frac_wtr_transp_down_sl_preload)
    IF (ALLOCATED(soil_forcing_cache%frac_w_lat_loss_sl_preload)) DEALLOCATE(soil_forcing_cache%frac_w_lat_loss_sl_preload)
    IF (ALLOCATED(soil_forcing_cache%wtr_soil_sl_preload)) DEALLOCATE(soil_forcing_cache%wtr_soil_sl_preload)
    IF (ALLOCATED(soil_forcing_cache%vol_porosity_sl_preload)) DEALLOCATE(soil_forcing_cache%vol_porosity_sl_preload)
    IF (ALLOCATED(soil_forcing_cache%t_soil_sl_preload)) DEALLOCATE(soil_forcing_cache%t_soil_sl_preload)
    IF (ALLOCATED(soil_forcing_cache%nhx_deposition_preload)) DEALLOCATE(soil_forcing_cache%nhx_deposition_preload)
    IF (ALLOCATED(soil_forcing_cache%noy_deposition_preload)) DEALLOCATE(soil_forcing_cache%noy_deposition_preload)
    IF (ALLOCATED(soil_forcing_cache%nhx_n15_deposition_preload)) DEALLOCATE(soil_forcing_cache%nhx_n15_deposition_preload)
    IF (ALLOCATED(soil_forcing_cache%noy_n15_deposition_preload)) DEALLOCATE(soil_forcing_cache%noy_n15_deposition_preload)
    IF (ALLOCATED(soil_forcing_cache%p_deposition_preload)) DEALLOCATE(soil_forcing_cache%p_deposition_preload)
    IF (ALLOCATED(soil_forcing_cache%veg_pool_fine_root_c_preload)) DEALLOCATE(soil_forcing_cache%veg_pool_fine_root_c_preload)
    IF (ALLOCATED(soil_forcing_cache%veg_litterfall_preload)) DEALLOCATE(soil_forcing_cache%veg_litterfall_preload)
    IF (ALLOCATED(soil_forcing_cache%seed_bed_litterfall_preload)) DEALLOCATE(soil_forcing_cache%seed_bed_litterfall_preload)
    IF (ALLOCATED(soil_forcing_cache%scratch_3d)) DEALLOCATE(soil_forcing_cache%scratch_3d)
    IF (ALLOCATED(soil_forcing_cache%scratch_4d)) DEALLOCATE(soil_forcing_cache%scratch_4d)

    soil_forcing_cache%arrays_initialized = .FALSE.
    soil_forcing_cache%preload_loaded = .FALSE.
    soil_forcing_cache%nproma = 0
    soil_forcing_cache%nsoil = 0
    soil_forcing_cache%nblks = 0
    soil_forcing_cache%npft = 0
    soil_forcing_cache%preload_nsteps = 0
    soil_forcing_cache%veg_litter_nr_parts = 0
    soil_forcing_cache%veg_litter_nr_elems = 0
    soil_forcing_cache%seed_bed_nr_elems = 0
    soil_forcing_cache%cached_year = -HUGE(0)
    soil_forcing_cache%cached_yearday = -HUGE(0)
    soil_forcing_cache%cached_timestep = -HUGE(0)

  END SUBROUTINE reset_soil_forcing_cache_arrays

  ! ======================================================================================================= !
  !> Reuse an open forcing file when possible, otherwise reopen it.
  !>
  SUBROUTINE open_cached_soil_forcing_file(model, input_file, is_open, cached_filename, target_filename, caller)

    USE mo_jsb_model_class, ONLY: t_jsb_model
    USE mo_jsb_io_netcdf,   ONLY: jsb_netcdf_open_input
    USE mo_jsb_parallel,    ONLY: my_process_is_stdio

    TYPE(t_jsb_model), POINTER, INTENT(in) :: model
    TYPE(t_input_file),        INTENT(inout) :: input_file
    LOGICAL,                   INTENT(inout) :: is_open
    CHARACTER(len=:), ALLOCATABLE, INTENT(inout) :: cached_filename
    CHARACTER(len=*),          INTENT(in) :: target_filename, caller

    IF (is_open .AND. ALLOCATED(cached_filename)) THEN
      IF (TRIM(cached_filename) == TRIM(target_filename)) THEN
        RETURN
      END IF
      CALL input_file%Close()
      is_open = .FALSE.
    END IF

    input_file = jsb_netcdf_open_input(TRIM(target_filename), model%grid_id)

    ! The standard JSBACH input wrapper can report is_open=.FALSE. on non-stdio ranks
    ! even though the distributed read path is still valid. Only the stdio rank should
    ! treat this as a hard open failure.
    IF (my_process_is_stdio() .AND. .NOT. input_file%is_open) THEN
      CALL finish(TRIM(caller), 'Failed to open soil forcing file '//TRIM(target_filename))
    END IF

    is_open = .TRUE.
    cached_filename = TRIM(target_filename)
    soil_forcing_cache%preload_loaded = .FALSE.
    soil_forcing_cache%cached_year = -HUGE(0)
    soil_forcing_cache%cached_yearday = -HUGE(0)
    soil_forcing_cache%cached_timestep = -HUGE(0)

  END SUBROUTINE open_cached_soil_forcing_file

  ! ======================================================================================================= !
  !> Read and cache the forcing fields from bc_qsoil_physics.
  !>
  SUBROUTINE load_soilphysics_cache(bgcm_store, nsteps, start_timestep, skip_isotope_reload, caller)

    TYPE(t_lnd_bgcm_store), POINTER, INTENT(in) :: bgcm_store
    INTEGER,                       INTENT(in) :: nsteps
    INTEGER,                       INTENT(in) :: start_timestep
    LOGICAL,                       INTENT(in) :: skip_isotope_reload
    CHARACTER(len=*),              INTENT(in) :: caller

    INTEGER :: ipft, itimestep, source_timestep
    CHARACTER(len=5) :: tile_name
    REAL(lp) :: t_total_start, t_total_end

    IF (debug_on()) CALL cpu_time(t_total_start)
    DO itimestep = 1, nsteps
      source_timestep = start_timestep + itimestep - 1
      DO ipft = 1, soil_forcing_cache%npft
        tile_name = get_pft_tile_name(ipft)

        CALL read_3D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
          & 'veg_root_fraction_sl_'//tile_name, soil_forcing_cache%scratch_4d, &
          & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
        soil_forcing_cache%root_fraction_sl_preload(:,:,:,ipft,itimestep) = MERGE( &
          & soil_forcing_cache%scratch_4d(:,:,:,1), zero, &
          & soil_forcing_cache%scratch_4d(:,:,:,1) >= zero)

        CALL read_2D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
          & 'veg_demand_uptake_n_tuptake_mavg_'//tile_name, soil_forcing_cache%scratch_3d, &
          & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
        soil_forcing_cache%demand_uptake_n_preload(:,:,ipft,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)

        CALL read_2D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
          & 'veg_demand_uptake_p_tuptake_mavg_'//tile_name, soil_forcing_cache%scratch_3d, &
          & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
        soil_forcing_cache%demand_uptake_p_preload(:,:,ipft,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)
      END DO

      CALL read_3D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'hydro_frac_wtr_transp_down_sl_box', soil_forcing_cache%scratch_4d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%frac_wtr_transp_down_sl_preload(:,:,:,itimestep) = soil_forcing_cache%scratch_4d(:,:,:,1)

      CALL read_3D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'hydro_frac_w_lat_loss_sl_box', soil_forcing_cache%scratch_4d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%frac_w_lat_loss_sl_preload(:,:,:,itimestep) = MERGE( &
        & soil_forcing_cache%scratch_4d(:,:,:,1), zero, &
        & soil_forcing_cache%scratch_4d(:,:,:,1) >= zero)

      CALL read_3D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'hydro_wtr_soil_sl_box', soil_forcing_cache%scratch_4d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%wtr_soil_sl_preload(:,:,:,itimestep) = MERGE( &
        & soil_forcing_cache%scratch_4d(:,:,:,1), zero, &
        & soil_forcing_cache%scratch_4d(:,:,:,1) >= zero)

      CALL read_3D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'hydro_vol_porosity_sl_box', soil_forcing_cache%scratch_4d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%vol_porosity_sl_preload(:,:,:,itimestep) = MERGE( &
        & soil_forcing_cache%scratch_4d(:,:,:,1), zero, &
        & soil_forcing_cache%scratch_4d(:,:,:,1) >= zero)

      CALL read_3D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'sse_t_soil_sl_box', soil_forcing_cache%scratch_4d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%t_soil_sl_preload(:,:,:,itimestep) = soil_forcing_cache%scratch_4d(:,:,:,1)

      CALL read_2D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'a2l_nhx_deposition_box', soil_forcing_cache%scratch_3d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%nhx_deposition_preload(:,:,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)

      CALL read_2D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'a2l_noy_deposition_box', soil_forcing_cache%scratch_3d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%noy_deposition_preload(:,:,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)

      IF (.NOT. skip_isotope_reload) THEN
        CALL read_2D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
          & 'a2l_nhx_n15_deposition_box', soil_forcing_cache%scratch_3d, &
          & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
        soil_forcing_cache%nhx_n15_deposition_preload(:,:,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)

        CALL read_2D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
          & 'a2l_noy_n15_deposition_box', soil_forcing_cache%scratch_3d, &
          & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
        soil_forcing_cache%noy_n15_deposition_preload(:,:,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)
      END IF

      CALL read_2D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'a2l_p_deposition_box', soil_forcing_cache%scratch_3d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%p_deposition_preload(:,:,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)
    END DO
    IF (debug_on()) THEN
      CALL cpu_time(t_total_end)
      CALL log_soil_forcing_timing(caller, 'soilphysics_read_total', t_total_end - t_total_start, &
        & 'records='//int2string(nsteps))
    END IF

  END SUBROUTINE load_soilphysics_cache

  ! ======================================================================================================= !
  !> Read and cache the forcing fields from bc_qsoil_litterfall.
  !>
  SUBROUTINE load_litterfall_cache(bgcm_store, nsteps, start_timestep, skip_isotope_reload, caller)

    TYPE(t_lnd_bgcm_store), POINTER, INTENT(in) :: bgcm_store
    INTEGER,                       INTENT(in) :: nsteps
    INTEGER,                       INTENT(in) :: start_timestep
    LOGICAL,                       INTENT(in) :: skip_isotope_reload
    CHARACTER(len=*),              INTENT(in) :: caller

    INTEGER :: ipft, i_part, i_elem, itimestep, source_timestep
    INTEGER :: veg_litter_idx, seed_bed_idx
    INTEGER :: veg_litter_store_idx, seed_bed_store_idx
    CHARACTER(len=5) :: tile_name
    REAL(lp) :: t_total_start, t_total_end

    tile_name = get_pft_tile_name(1)
    veg_litter_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_LITTERFALL_ID, tile_name, caller)
    seed_bed_idx   = get_bgcm_idx(bgcm_store, VEG_BGCM_SEED_BED_LITTERFALL_ID, tile_name, caller)
    veg_litter_store_idx = bgcm_store%idx_in_store(veg_litter_idx)
    seed_bed_store_idx   = bgcm_store%idx_in_store(seed_bed_idx)

    IF (debug_on()) CALL cpu_time(t_total_start)
    DO itimestep = 1, nsteps
      source_timestep = start_timestep + itimestep - 1
      DO ipft = 1, soil_forcing_cache%npft
        tile_name = get_pft_tile_name(ipft)

        CALL read_2D_extdim(soil_forcing_cache%veglitterfall_input_file%stream_id, on_cells, &
          & 'veg_veg_bgcm_pool_fine_root_C_'//tile_name, soil_forcing_cache%scratch_3d, &
          & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
        soil_forcing_cache%veg_pool_fine_root_c_preload(:,:,ipft,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)

        DO i_part = 1, soil_forcing_cache%veg_litter_nr_parts
          DO i_elem = 1, soil_forcing_cache%veg_litter_nr_elems
            IF (skip_isotope_reload .AND. i_elem >= 4) CYCLE
            CALL read_2D_extdim(soil_forcing_cache%veglitterfall_input_file%stream_id, on_cells, &
              & 'veg_veg_bgcm_litterfall_'//get_veg_litter_part_name(i_part)//'_'// &
              & get_veg_litter_elem_name(i_elem)//'_'//tile_name, soil_forcing_cache%scratch_3d, &
              & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
            soil_forcing_cache%veg_litterfall_preload(i_part, i_elem, :, :, ipft, itimestep) = &
              & soil_forcing_cache%scratch_3d(:,:,1)
          END DO
        END DO

        DO i_elem = 1, soil_forcing_cache%seed_bed_nr_elems
          IF (skip_isotope_reload .AND. i_elem >= 4) CYCLE
          CALL read_2D_extdim(soil_forcing_cache%veglitterfall_input_file%stream_id, on_cells, &
            & 'veg_veg_bgcm_seed_bed_litterfall_'//get_veg_litter_elem_name(i_elem)//'_'//tile_name, &
            & soil_forcing_cache%scratch_3d, start_extdim=source_timestep, end_extdim=source_timestep, &
            & extdim_name='time')
          soil_forcing_cache%seed_bed_litterfall_preload(i_elem, :, :, ipft, itimestep) = &
            & soil_forcing_cache%scratch_3d(:,:,1)
        END DO
      END DO
    END DO
    IF (debug_on()) THEN
      CALL cpu_time(t_total_end)
      CALL log_soil_forcing_timing(caller, 'litterfall_read_total', t_total_end - t_total_start, &
        & 'records='//int2string(nsteps))
    END IF

  END SUBROUTINE load_litterfall_cache

  ! ======================================================================================================= !
  !> Copy cached litterfall and fine-root C into the tile bgcm storage.
  !>
  SUBROUTINE apply_litterfall_cache_to_tile(bgcm_store, tile_name, pft_idx, itimestep, caller)

    TYPE(t_lnd_bgcm_store), POINTER, INTENT(inout) :: bgcm_store
    CHARACTER(len=*),              INTENT(in)    :: tile_name, caller
    INTEGER,                       INTENT(in)    :: pft_idx, itimestep

    INTEGER :: veg_pool_idx, veg_litter_idx, seed_bed_idx
    INTEGER :: veg_pool_store_idx, veg_litter_store_idx, seed_bed_store_idx
    INTEGER :: i_part, i_elem
    REAL(lp) :: t_total_start, t_total_end

    IF (debug_on() .AND. TRIM(tile_name) == 'pft01') CALL cpu_time(t_total_start)

    veg_pool_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_POOL_ID, tile_name, caller)
    veg_litter_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_LITTERFALL_ID, tile_name, caller)
    seed_bed_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_SEED_BED_LITTERFALL_ID, tile_name, caller)

    veg_pool_store_idx = bgcm_store%idx_in_store(veg_pool_idx)
    veg_litter_store_idx = bgcm_store%idx_in_store(veg_litter_idx)
    seed_bed_store_idx = bgcm_store%idx_in_store(seed_bed_idx)

    bgcm_store%store_2l_2d_bgcms(veg_pool_store_idx)%collection_2l_p(ix_fine_root, ixC)%p%ptr2d(:,:) = &
      & soil_forcing_cache%veg_pool_fine_root_c_preload(:,:,pft_idx,itimestep)

    DO i_part = 1, soil_forcing_cache%veg_litter_nr_parts
      DO i_elem = 1, soil_forcing_cache%veg_litter_nr_elems
        bgcm_store%store_2l_2d_bgcms(veg_litter_store_idx)%collection_2l_p(i_part, i_elem)%p%ptr2d(:,:) = &
          & soil_forcing_cache%veg_litterfall_preload(i_part, i_elem, :, :, pft_idx, itimestep)
      END DO
    END DO

    DO i_elem = 1, soil_forcing_cache%seed_bed_nr_elems
      bgcm_store%store_1l_2d_bgcms(seed_bed_store_idx)%collection_1l_p(i_elem)%p%ptr2d(:,:) = &
        & soil_forcing_cache%seed_bed_litterfall_preload(i_elem, :, :, pft_idx, itimestep)
    END DO
    IF (debug_on() .AND. TRIM(tile_name) == 'pft01') THEN
      CALL cpu_time(t_total_end)
      CALL log_soil_forcing_timing(caller, 'apply_litterfall tile='//TRIM(tile_name), t_total_end - t_total_start)
    END IF

  END SUBROUTINE apply_litterfall_cache_to_tile

  ! ======================================================================================================= !
  !> Resolve soil forcing filename and record index for the current model time.
  !> The filename argument is interpreted as a full filename for yearly forcing
  !> and as a prefix for daily forcing, yielding <prefix>_<doy>.nc.
  !>
  SUBROUTINE get_soil_forcing_info(base_filename, yearday, hour, minute, second, dtime, &
      & resolved_filename, itimestep)

      CHARACTER(len=*),           INTENT(in)  :: base_filename
      INTEGER,                    INTENT(in)  :: yearday, hour, minute, second
      REAL(lp),                   INTENT(in)  :: dtime
      CHARACTER(len=:), ALLOCATABLE, INTENT(out) :: resolved_filename
      INTEGER,                    INTENT(out) :: itimestep

    INTEGER :: seconds_of_day
    CHARACTER(len=8) :: doy_string
    LOGICAL :: is_daily_timestep

    seconds_of_day = hour * 3600 + minute * 60 + second
    is_daily_timestep = NINT(dtime) == 86400

    IF (is_daily_timestep) THEN
      resolved_filename = TRIM(base_filename)
      itimestep = yearday
    ELSE
      WRITE(doy_string,'(I3.3)') yearday
      resolved_filename = trim_nc_filename(base_filename)//'_'//TRIM(doy_string)//'.nc'
      itimestep = seconds_of_day / INT(dtime) + 1
    END IF

  END SUBROUTINE get_soil_forcing_info

  ! ======================================================================================================= !
  !> Return the filename prefix without a trailing .nc suffix.
  !>
  FUNCTION trim_nc_filename(filename) RESULT(prefix)

    CHARACTER(len=*), INTENT(in) :: filename
    CHARACTER(len=:), ALLOCATABLE :: prefix
    INTEGER :: ilen

    ilen = LEN_TRIM(filename)
    IF (ilen >= 3 .AND. filename(ilen-2:ilen) == '.nc') THEN
      prefix = filename(:ilen-3)
    ELSE
      prefix = TRIM(filename)
    END IF

  END FUNCTION trim_nc_filename

  ! ======================================================================================================= !
  !>
  !> @TODO add docu
  !>
  FUNCTION IMerge_TIME2Sec (kd, kh, km, ks) RESULT (ix)
    !+
    !
    ! IMerge_TIME2Sec [function, integer]
    !    merge day, hour, minute and second to internal time
    !    (
    !    day    [integer] input (day of year, 0...364)
    !    hour   [integer] input (hour of day, 0...23)
    !    minute [integer] input (minute of hour, 0...59)
    !    second [integer] input (second of minute, 0...59)
    !    )
    !
    !-
    INTEGER, INTENT(in) :: kd, kh, km, ks
    INTEGER             :: ix

    ix = (kd-1)*86400 + kh*3600 + km*60 + ks

  END FUNCTION IMerge_TIME2Sec
#endif

#endif
END MODULE mo_sb_spinup_qsoil_forcing
