!> Forward soil-biogeochemistry spin-up accelerator
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
!> Event driver for replaying bookkept soil forcing during spin-up.
!> The forward relaxation is shared by JSM and simple_1d and starts from the
!> current soil state. An optional conditional carbon seed precedes later replays.
!>
!> Memory contract:
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
  USE mo_sb_config_class, ONLY: sb_spinup_accelerator_n_forcing_bins
  USE mo_sb_plant_uptake_adapter, ONLY: t_sb_plant_uptake_drivers

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: apply_sb_spinup_accelerator
  PUBLIC :: resume_sb_spinup_accelerator
  PUBLIC :: apply_sb_coupled_pore_calibration

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_spinup_accelerator'
  INTEGER, PARAMETER :: sb_spinup_debug_years(10) = [1, 2, 5, 10, 50, 100, 500, 1000, 5000, 10000]

  !> Event-local climatology reconstructed from restart-safe accumulated sums.
  TYPE :: t_sb_spinup_forcing_climatology
    INTEGER :: n_forcing_bins = 0
    INTEGER :: n_daily_bins = 0 !< 365, or 360 for the host 360-day calendar [-]
    INTEGER :: n_soil_layers = 0
    INTEGER :: n_litter_sources = 0
    INTEGER :: n_elements = 0
    ! Only a combined-column context supplies these seasonal PFT drivers.
    ! Ordinary per-PFT replay and file-driven QSOIL retain their original path.
    TYPE(t_sb_plant_uptake_drivers), ALLOCATABLE :: uptake_drivers(:)
    REAL(lp), ALLOCATABLE :: partitioned_litter_rate(:,:,:,:,:)
      !< Combined pool/element/cell/layer/month litter [mol m-3 s-1]
    REAL(lp) :: litter_partition_timestep = 0._lp
      !< Substep used for PFT-specific partitioning; must match replay [s]
    REAL(lp), ALLOCATABLE :: litter_to_uptake_n(:), litter_to_uptake_p(:) !< coupled reference ratios [-]
    ! Keep original budget integrals for the combined-column adapter: averaging
    ! litter/uptake ratios would not conserve either reference budget.
    REAL(lp), ALLOCATABLE :: reference_litter_n(:), reference_litter_p(:) !< Window inputs [mol m-2]
    REAL(lp), ALLOCATABLE :: reference_uptake_n(:), reference_uptake_p(:) !< Window uptake [mol m-2]
    REAL(lp), ALLOCATABLE :: daily_time(:,:) !< Actual daily sampling weights [s]
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
    REAL(lp), ALLOCATABLE :: growth_p_limit(:,:)    !< mean transformed P limitation by seasonal bin/layer [-]
    REAL(lp), ALLOCATABLE :: root_fraction(:,:)    !< mean root fraction, flattened by bin/layer [-]
    REAL(lp), ALLOCATABLE :: wtr_down(:,:)         !< mean downward transport rate, flattened by daily bin/layer [s-1]
    REAL(lp), ALLOCATABLE :: wtr_lat_loss(:,:)     !< mean lateral transport rate, flattened by daily bin/layer [s-1]
    REAL(lp), ALLOCATABLE :: t_soil(:,:)           !< mean soil temperature, flattened by bin/layer [K]
    REAL(lp), ALLOCATABLE :: wtr_soil(:,:)         !< mean liquid soil water, flattened by bin/layer [m]
    REAL(lp), ALLOCATABLE :: ice_soil(:,:)         !< mean soil ice, flattened by bin/layer [m]
    REAL(lp), ALLOCATABLE :: wtr_pot(:,:)          !< mean soil water potential, flattened by bin/layer [MPa]
  END TYPE t_sb_spinup_forcing_climatology

  !> Coupled-model forcing values preserved while the soil-only replay is active.
  !! This is the normal per-year forcing restore, NOT an atomic rollback.
  !! The separate replay checkpoint retains process state and both BGCM
  !! representations for all owners before any shared-state reconstruction.
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
    REAL(lp), ALLOCATABLE :: growth_p_limit(:,:)  !< transformed P driver by cell and soil layer [-]
    REAL(lp), ALLOCATABLE :: root_fraction(:,:) !< root fraction by cell and soil layer [-]
    REAL(lp), ALLOCATABLE :: wtr_down(:,:)      !< native downward transfer coefficient [1000 s-1]
    REAL(lp), ALLOCATABLE :: wtr_lat_loss(:,:)  !< native lateral transfer coefficient [1000 s-1]
    REAL(lp), ALLOCATABLE :: t_soil(:,:)        !< soil temperature by cell and layer [K]
    REAL(lp), ALLOCATABLE :: wtr_soil(:,:)      !< liquid soil water by cell and layer [m]
    REAL(lp), ALLOCATABLE :: ice_soil(:,:)      !< soil ice by cell and layer [m]
    REAL(lp), ALLOCATABLE :: wtr_pot(:,:)       !< soil water potential by cell and layer [MPa]
    REAL(lp), ALLOCATABLE :: net_biosphere_production(:) !< host-timestep net biosphere production [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: net_biosphere_production_c13(:) !< host-timestep 13C NBP [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: net_biosphere_production_c14(:) !< host-timestep 14C NBP [micro-mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: biological_n_fixation(:) !< host-timestep biological N fixation [micro-mol m-2 s-1]
  END TYPE t_sb_spinup_live_forcing

  !> Event-local annual nutrient budget and next-year demand correction.
  TYPE :: t_sb_replay_nutrient_balance
    REAL(lp), ALLOCATABLE :: litter_n(:), litter_p(:) !< litter received in one replay year [mol m-2]
    REAL(lp), ALLOCATABLE :: uptake_n(:), uptake_p(:) !< realized plant uptake in one replay year [mol m-2]
    REAL(lp), ALLOCATABLE :: scale_n(:), scale_p(:) !< multipliers of recorded effective demand [-]
  END TYPE t_sb_replay_nutrient_balance

  ! One grid cell at a time, with all calculated tiles at the same artificial year.
  ! Context is event-local; no checkpoint can be written midway through replay.
  TYPE :: t_sb_grid_replay_tile
    LOGICAL :: active = .FALSE. !< Tile has positive cell area and a pending event [-]
    LOGICAL :: pending = .FALSE. !< Event has not yet been committed [-]
    TYPE(t_sb_spinup_forcing_climatology) :: forcing
    TYPE(t_sb_spinup_live_forcing) :: live
    TYPE(t_sb_replay_nutrient_balance) :: balance
    REAL(lp), ALLOCATABLE :: before(:,:,:,:) !< Optional debug snapshot [mol m-3]
  END TYPE t_sb_grid_replay_tile

CONTAINS

  !> Construct one cell's shared forcing without writing any process memory.
  !> Geometry and sampling coverage must agree: mixing unlike columns or
  !> unmatched climatological windows would not conserve their source budgets.
  !> Means are per represented vegetated area; nutrient ratios are calculated
  !> from combined integrals, never averaged from PFT-specific ratios.
  SUBROUTINE prepare_combined_sb_forcing(model, options, context, shared, accepted, reason)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    USE mo_jsb_process_class, ONLY: HYDRO_, VEG_, SB_
    USE mo_hydro_memory_class, ONLY: t_hydro_memory
    USE mo_veg_memory_class, ONLY: t_veg_memory
    USE mo_sb_config_class, ONLY: t_sb_config
    USE mo_lnd_bgcm_store_class, ONLY: SB_BGCM_POOL_ID
    USE mo_q_sb_litter_processes, ONLY: combine_sb_litter_inputs
    USE mo_sb_plant_uptake_adapter, ONLY: t_sb_plant_uptake, prepare_combined_plant_uptake
    USE mo_jsb_math_constants, ONLY: one_day, eps8
    USE mo_jsb_physical_constants, ONLY: Tzero
    USE mo_sb_constants, ONLY: tacclim_uptake_min, tacclim_uptake_max, k_ea_uptake, &
      & ea_decomposition, ed_decomposition, t_ref_decomposition
    USE mo_phy_schemes, ONLY: calc_peaked_arrhenius_function
    TYPE(t_jsb_model), POINTER, INTENT(in) :: model !< Owning model; no mutation
    TYPE(t_jsb_task_options), INTENT(in) :: options !< Single cell and unchanged host clock
    TYPE(t_sb_grid_replay_tile), INTENT(in) :: context(:) !< Prepared PFT snapshots
    TYPE(t_sb_spinup_forcing_climatology), INTENT(out) :: shared !< Scratch shared forcing
    LOGICAL, INTENT(out) :: accepted !< All construction checks passed [-]
    CHARACTER(len=*), INTENT(out) :: reason !< Controlled rejection/ready diagnostic
    TYPE(t_hydro_memory), POINTER :: hydro
    TYPE(t_veg_memory), POINTER :: veg
    TYPE(t_sb_config), POINTER :: conf
    TYPE(t_sb_plant_uptake) :: checked_uptake !< Scratch coefficient validation; recomputed during replay
    REAL(lp), POINTER :: pools(:,:,:,:) !< Read-only BGCM shape reference [mol m-3]
    REAL(lp), ALLOCATABLE :: depth(:,:), fraction(:,:), roots(:,:,:), sla(:), acclim(:)
    REAL(lp), ALLOCATABLE :: litter(:,:,:,:), seed(:,:,:)
    REAL(lp), ALLOCATABLE :: potential(:,:) !< Common seasonal soil water potential [MPa]
    INTEGER, ALLOCATABLE :: members(:), growthform(:)
    REAL(lp) :: nactive(1), area, old_weight, new_weight, dt, root_stock
    REAL(lp) :: reference_litter(2), reference_uptake(2), reference_ratio(2)
    INTEGER :: it, j, first, c, block, ns, nb, ne, nv, np, nt, b, s, e, v, entry, lib
    LOGICAL :: valid, interactive(2)
    CHARACTER(len=16) :: uptake_scheme !< Existing ECA/simple_1d kernel selection [-]

    accepted = .FALSE.
    reason = 'Combined forcing requires a single cell and at least one active PFT.'
    IF (options%nc /= 1 .OR. options%ics /= options%ice) RETURN
    nt = COUNT(context%active)
    IF (nt == 0) RETURN
    members = PACK([(it,it=1,SIZE(context))],context%active)
    first = members(1)
    c = options%ics
    block = options%iblk
    SELECT TYPE (cfg => model%processes(SB_)%p%config)
    TYPE IS (t_sb_config)
      conf => cfg
    END SELECT
    ns = context(first)%forcing%n_soil_layers
    uptake_scheme = conf%sb_model_scheme
    IF (uptake_scheme == 'jsm') uptake_scheme = conf%sb_adsorp_scheme
    nb = context(first)%forcing%n_forcing_bins
    ne = context(first)%forcing%n_elements
    nv = context(first)%forcing%n_litter_sources
    pools => model%tiles(first)%p%bgcm_store%Get_matrix_from_store_2l_3d( &
      & SB_BGCM_POOL_ID,c,c,block,model%tiles(first)%p%name)
    np = SIZE(pools,1)
    ALLOCATE(depth(1,ns),fraction(1,nt),roots(1,ns,nt),sla(nt),growthform(nt),acclim(nt))
    ALLOCATE(litter(nv,ne,1,nt),seed(ne,1,nt),potential(1,ns))
    area = 0._lp
    DO j = 1,nt
      it = members(j)
      fraction(1,j) = model%tiles(it)%p%fract(c,block)
      reason = 'Combined forcing: invalid active PFT area.'
      IF (.NOT. ieee_is_finite(fraction(1,j))) RETURN
      IF (fraction(1,j) <= 0._lp) RETURN
      SELECT TYPE (mem => model%tiles(it)%p%mem(HYDRO_)%p)
      TYPE IS (t_hydro_memory)
        hydro => mem
      END SELECT
      SELECT TYPE (mem => model%tiles(it)%p%mem(VEG_)%p)
      TYPE IS (t_veg_memory)
        veg => mem
      END SELECT
      ASSOCIATE(f => context(it)%forcing)
        reason = 'Combined forcing: incompatible dimensions, soil geometry or sampling windows.'
        IF (f%n_soil_layers /= ns .OR. f%n_forcing_bins /= nb .OR. &
          & f%n_elements /= ne .OR. f%n_litter_sources /= nv) RETURN
        IF (j == 1) THEN
          shared = f
          depth(1,:) = hydro%soil_depth_sl%ptr(c,:,block)
          nactive(1) = hydro%num_sl_above_bedrock%ptr(c,block)
          IF (.NOT. ALL(ieee_is_finite(depth)) .OR. .NOT. ieee_is_finite(nactive(1))) RETURN
          IF (nactive(1) < 1._lp .OR. nactive(1) > REAL(ns,lp)) RETURN
          IF (nactive(1) /= REAL(NINT(nactive(1)),lp)) RETURN
          IF (ANY(depth(1,1:NINT(nactive(1))) <= 0._lp)) RETURN
        ELSE
          IF (f%n_daily_bins /= shared%n_daily_bins) RETURN
          IF (ANY(SHAPE(f%forcing_time) /= SHAPE(shared%forcing_time))) RETURN
          IF (ANY(SHAPE(f%daily_time) /= SHAPE(shared%daily_time))) RETURN
          IF (ANY(f%forcing_time /= shared%forcing_time)) RETURN
          IF (ANY(f%daily_time /= shared%daily_time)) RETURN
          IF (ANY(hydro%soil_depth_sl%ptr(c,:,block) /= depth(1,:))) RETURN
          IF (hydro%num_sl_above_bedrock%ptr(c,block) /= nactive(1)) RETURN
          ! Incremental area mean avoids summing raw absolute-area inventories.
          old_weight = area/(area+fraction(1,j))
          new_weight = fraction(1,j)/(area+fraction(1,j))
          shared%litterfall_rate = old_weight*shared%litterfall_rate + new_weight*f%litterfall_rate
          shared%fine_root_c = old_weight*shared%fine_root_c + new_weight*f%fine_root_c
          shared%f_n_demand = old_weight*shared%f_n_demand + new_weight*f%f_n_demand
          shared%f_p_demand = old_weight*shared%f_p_demand + new_weight*f%f_p_demand
          shared%reference_litter_n = old_weight*shared%reference_litter_n + new_weight*f%reference_litter_n
          shared%reference_litter_p = old_weight*shared%reference_litter_p + new_weight*f%reference_litter_p
          shared%reference_uptake_n = old_weight*shared%reference_uptake_n + new_weight*f%reference_uptake_n
          shared%reference_uptake_p = old_weight*shared%reference_uptake_p + new_weight*f%reference_uptake_p
          shared%nhx_deposition = old_weight*shared%nhx_deposition + new_weight*f%nhx_deposition
          shared%noy_deposition = old_weight*shared%noy_deposition + new_weight*f%noy_deposition
          shared%nhx_n15_deposition = old_weight*shared%nhx_n15_deposition + new_weight*f%nhx_n15_deposition
          shared%noy_n15_deposition = old_weight*shared%noy_n15_deposition + new_weight*f%noy_n15_deposition
          shared%p_deposition = old_weight*shared%p_deposition + new_weight*f%p_deposition
          shared%growth_p_limit = old_weight*shared%growth_p_limit + new_weight*f%growth_p_limit
          shared%wtr_down = old_weight*shared%wtr_down + new_weight*f%wtr_down
          shared%wtr_lat_loss = old_weight*shared%wtr_lat_loss + new_weight*f%wtr_lat_loss
          shared%t_soil = old_weight*shared%t_soil + new_weight*f%t_soil
          shared%wtr_soil = old_weight*shared%wtr_soil + new_weight*f%wtr_soil
          shared%ice_soil = old_weight*shared%ice_soil + new_weight*f%ice_soil
          shared%wtr_pot = old_weight*shared%wtr_pot + new_weight*f%wtr_pot
        END IF
      END ASSOCIATE
      area = area + fraction(1,j)
      lib = model%tiles(it)%p%lcts(1)%lib_id
      sla(j) = model%lctlib(lib)%sla
      growthform(j) = model%lctlib(lib)%growthform
      ! Preserve each PFT's lagged acclimation; no vegetation averaging is advanced.
      acclim(j) = veg%t_soil_root_tacclim_mavg%ptr(c,block)
    END DO
    ! Match the existing inactive-budget policy, but return a recoverable
    ! construction rejection rather than terminating from the ratio helper.
    reference_litter = [shared%reference_litter_n(1),shared%reference_litter_p(1)]
    reference_uptake = [shared%reference_uptake_n(1),shared%reference_uptake_p(1)]
    interactive = [conf%flag_sb_interactive_n,conf%flag_sb_interactive_p]
    reference_ratio = 0._lp
    DO e = 1,2
      IF (.NOT. interactive(e)) CYCLE
      reason = 'Combined forcing: invalid or one-sided inactive litter/uptake budget.'
      IF (.NOT. ieee_is_finite(reference_litter(e)) .OR. .NOT. ieee_is_finite(reference_uptake(e))) RETURN
      IF (MIN(reference_litter(e),reference_uptake(e)) < 0._lp) RETURN
      IF (MAX(reference_litter(e),reference_uptake(e)) <= eps8) CYCLE
      IF (MIN(reference_litter(e),reference_uptake(e)) <= eps8) RETURN
      reference_ratio(e) = reference_litter(e)/reference_uptake(e)
    END DO
    shared%litter_to_uptake_n = reference_ratio(1)
    shared%litter_to_uptake_p = reference_ratio(2)
    dt = one_day/REAL(MAX(1,CEILING(1._lp/conf%sb_model_spinup_accelerator_timestep_days)),lp)
    shared%litter_partition_timestep = dt
    ALLOCATE(shared%partitioned_litter_rate(np,ne,1,ns,nb),shared%uptake_drivers(nb))
    shared%root_fraction = 0._lp
    DO b = 1,nb
      ASSOCIATE(drivers => shared%uptake_drivers(b))
        drivers%fraction = fraction
        ALLOCATE(drivers%roots(1,ns,nt),drivers%temperature(1,ns,nt),drivers%phi_leaf_min(nt), &
          & drivers%vmax(2,nt),drivers%demand(2,1,nt))
        drivers%roots = 0._lp
        drivers%temperature = 0._lp
        DO j = 1,nt
          it = members(j)
          lib = model%tiles(it)%p%lcts(1)%lib_id
          drivers%phi_leaf_min(j) = model%lctlib(lib)%phi_leaf_min
          drivers%vmax(:,j) = [model%lctlib(lib)%vmax_uptake_n,model%lctlib(lib)%vmax_uptake_p]
          ASSOCIATE(f => context(it)%forcing)
            drivers%demand(:,1,j) = [f%f_n_demand(1,b),f%f_p_demand(1,b)]
            DO s = 1,ns
              entry = b+nb*(s-1)
              roots(1,s,j) = f%root_fraction(1,entry)
              root_stock = f%fine_root_c(1,b)*roots(1,s,j)
              shared%root_fraction(1,entry) = shared%root_fraction(1,entry) + fraction(1,j)/area*root_stock
              IF (s > NINT(nactive(1))) CYCLE
              reason = 'Combined forcing: invalid soil depth or root acclimation.'
              IF (depth(1,s) <= 0._lp .OR. .NOT. ieee_is_finite(acclim(j))) RETURN
              drivers%roots(1,s,j) = root_stock/depth(1,s)
              drivers%temperature(1,s,j) = calc_peaked_arrhenius_function(shared%t_soil(1,entry), &
                & ea_decomposition*k_ea_uptake* &
                & MIN(MAX(acclim(j)-Tzero,tacclim_uptake_min),tacclim_uptake_max), &
                & ed_decomposition,t_ref_decomposition)
            END DO
            DO e = 1,ne
              DO v = 1,nv
                entry = b+nb*((v-1)+(nv+1)*(e-1))
                litter(v,e,1,j) = f%litterfall_rate(1,entry)
              END DO
              seed(e,1,j) = f%litterfall_rate(1,b+nb*(nv+(nv+1)*(e-1)))
            END DO
          END ASSOCIATE
        END DO
        DO s = 1,ns
          potential(1,s) = shared%wtr_pot(1,b+nb*(s-1))
        END DO
        CALL prepare_combined_plant_uptake(uptake_scheme,drivers%fraction,drivers%roots, &
          & drivers%temperature,potential,drivers%phi_leaf_min,drivers%vmax,drivers%demand,checked_uptake,valid)
        reason = 'Combined forcing: invalid grouped root uptake drivers.'
        IF (.NOT. valid) RETURN
      END ASSOCIATE
      DO s = 1,ns
        entry = b+nb*(s-1)
        IF (shared%fine_root_c(1,b) > 0._lp) THEN
          shared%root_fraction(1,entry) = shared%root_fraction(1,entry)/shared%fine_root_c(1,b)
        ELSE
          shared%root_fraction(1,entry) = 0._lp
        END IF
      END DO
      CALL combine_sb_litter_inputs(1,ns,np,ne,nv,nt,dt,conf%sb_model_scheme,fraction,nactive,depth, &
        & sla,growthform,roots,litter,seed,shared%partitioned_litter_rate(:,:,:,:,b),valid)
      reason = 'Combined forcing: per-PFT litter partition or elemental budget rejected.'
      IF (.NOT. valid) RETURN
    END DO
    accepted = .TRUE.
    reason = 'Combined forcing prepared for transactional shared-state execution.'
  END SUBROUTINE prepare_combined_sb_forcing

  ! ===================================================================================================== !
  !> Identify the next accelerator event and whether its preceding forcing window is currently active.
  SUBROUTINE get_sb_spinup_accelerator_forcing_window(tile, options, l_collect, event_year, execution_index)

    USE mo_jsb_process_class, ONLY: SB_
    USE mo_jsb_time,          ONLY: get_year
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    LOGICAL,                    INTENT(out)   :: l_collect
    INTEGER,                    INTENT(out)   :: event_year
    INTEGER,                    INTENT(out)   :: execution_index

    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_execution_count

    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: current_year
    INTEGER :: i_event, ics, ice, iblk

    l_collect = .FALSE.
    event_year = 0
    execution_index = 0

    model => Get_model(tile%owner_model_id)
    IF (.NOT. model%config%l_spinup) RETURN
    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_execution_count)

    current_year = get_year(options%current_datetime)
    DO i_event = 1, dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_max_executions
      ! Completed events must not hide the start of the next collection window
      ! in the same calendar year, notably when forcing years equal frequency.
      IF (ALL(NINT(sb_spinup_accelerator_execution_count) >= i_event)) CYCLE
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
  SUBROUTINE get_sb_spinup_accelerator_forcing_bin(options, forcing_bin, n_forcing_bins, daily_bin)
    USE mo_jsb_time, ONLY: get_year, get_year_day, get_year_length
    TYPE(t_jsb_task_options), INTENT(in) :: options
    INTEGER, INTENT(out) :: forcing_bin, n_forcing_bins, daily_bin
    INTEGER :: year_length, day_index

    year_length = get_year_length(get_year(options%current_datetime))
    day_index = MAX(0, MIN(year_length - 1, FLOOR(get_year_day(options%current_datetime) - 1.0_lp)))
    n_forcing_bins = sb_spinup_accelerator_n_forcing_bins
    ! Match QS equal-phase seasonal bins, not calendar months. Compress leap
    ! years onto 365 daily bins, and preserve all 360 days in a 360-day calendar.
    forcing_bin = day_index * n_forcing_bins / year_length + 1
    daily_bin = day_index * MIN(365, year_length) / year_length + 1
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
    INTEGER :: daily_bin, daily_profile_entry, n_daily_bins
    INTEGER :: n_forcing_bins
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    INTEGER :: ic !< Cell index for restart-completeness checks [-]
    INTEGER :: i_source
    INTEGER :: i_element
    INTEGER :: i_soil
    INTEGER :: litter_entry
    INTEGER :: profile_entry
    INTEGER :: n_veg_sources
    INTEGER :: n_elements
    INTEGER :: nsoil
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_litter_n_sum
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_litter_p_sum
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_plant_uptake_n_sum
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_plant_uptake_p_sum
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_forcing_event
    dsl4jsb_Real3D_onChunk :: sb_accel_cue_sum_sl
    dsl4jsb_Real3D_onChunk :: sb_accel_poly_sum_sl
    dsl4jsb_Real3D_onChunk :: sb_accel_ap_sum_sl
    dsl4jsb_Real3D_onChunk :: microbial_cue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly
    dsl4jsb_Real3D_onChunk :: enzyme_frac_AP
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_daily_time
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
    dsl4jsb_Real3D_onChunk :: growth_p_limit_smoothed_sl
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

    CALL get_sb_spinup_accelerator_forcing_bin(options, forcing_bin, n_forcing_bins, daily_bin)
    n_daily_bins = dsl4jsb_Config(SB_)%sb_spinup_accelerator_n_daily_bins

    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_litter_n_sum)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_litter_p_sum)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_plant_uptake_n_sum)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_plant_uptake_p_sum)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_forcing_event)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_accel_cue_sum_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_accel_poly_sum_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_accel_ap_sum_sl)
    IF (dsl4jsb_Config(SB_)%sb_model_scheme == 'jsm') THEN
      dsl4jsb_Get_var3D_onChunk(SB_, microbial_cue_eff_tmic_mavg)
      dsl4jsb_Get_var3D_onChunk(SB_, enzyme_frac_poly)
      dsl4jsb_Get_var3D_onChunk(SB_, enzyme_frac_AP)
    END IF
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_daily_time)
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
    dsl4jsb_Get_var3D_onChunk(SB_, growth_p_limit_smoothed_sl)
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

    ! Never merge stale forcing from another event after a schedule change.
    ! Empty storage can be assigned; populated storage must have the same owner.
    IF (ANY(NINT(sb_spinup_accelerator_forcing_event) /= event_year)) THEN
      IF (ANY(sb_spinup_accelerator_forcing_time > 0.0_lp)) THEN
        CALL finish(TRIM(routine), 'Replay forcing belongs to a different event; recollect the forcing window.')
      END IF
      sb_spinup_accelerator_forcing_event = REAL(event_year, lp)
    END IF
    ! A missing restart integral remains invalid for a partially collected window.
    ! Initialize only genuinely empty windows; otherwise reject the analytical seed later.
    DO ic = 1,options%nc
      IF (SUM(sb_spinup_accelerator_forcing_time(ic,:)) == 0._lp) THEN
        sb_accel_cue_sum_sl(ic,:) = 0._lp
        sb_accel_poly_sum_sl(ic,:) = 0._lp
        sb_accel_ap_sum_sl(ic,:) = 0._lp
      END IF
      IF (ANY(sb_accel_cue_sum_sl(ic,:) < 0._lp) .OR. ANY(sb_accel_poly_sum_sl(ic,:) < 0._lp) .OR. &
        & ANY(sb_accel_ap_sum_sl(ic,:) < 0._lp)) CYCLE
      IF (dsl4jsb_Config(SB_)%sb_model_scheme /= 'jsm') CYCLE
      sb_accel_cue_sum_sl(ic,:) = sb_accel_cue_sum_sl(ic,:) + microbial_cue_eff_tmic_mavg(ic,:)*options%dtime
      sb_accel_poly_sum_sl(ic,:) = sb_accel_poly_sum_sl(ic,:) + enzyme_frac_poly(ic,:)*options%dtime
      sb_accel_ap_sum_sl(ic,:) = sb_accel_ap_sum_sl(ic,:) + enzyme_frac_AP(ic,:)*options%dtime
    END DO
    sb_spinup_accelerator_daily_time(:, daily_bin) = &
      & sb_spinup_accelerator_daily_time(:, daily_bin) + options%dtime
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
      sb_spinup_accelerator_growth_p_limit_sum(:, profile_entry) = &
        & sb_spinup_accelerator_growth_p_limit_sum(:, profile_entry) + &
        & growth_p_limit_smoothed_sl(:, i_soil) * options%dtime
      daily_profile_entry = daily_bin + n_daily_bins * (i_soil - 1)
      sb_spinup_accelerator_root_fraction_sum(:, profile_entry) = &
        & sb_spinup_accelerator_root_fraction_sum(:, profile_entry) + root_fraction_sl(:, i_soil) * options%dtime
      ! HYDRO coefficients are 1000 times the true per-second rate.
      sb_spinup_accelerator_wtr_down_sum(:, daily_profile_entry) = &
        & sb_spinup_accelerator_wtr_down_sum(:, daily_profile_entry) + &
        & MAX(0.0_lp, frac_wtr_transp_down_sl(:, i_soil)) * options%dtime / 1000.0_lp
      ! HYDRO coefficients are 1000 times the true per-second rate.
      sb_spinup_accelerator_wtr_lat_loss_sum(:, daily_profile_entry) = &
        & sb_spinup_accelerator_wtr_lat_loss_sum(:, daily_profile_entry) + &
        & MAX(0.0_lp, frac_w_lat_loss_sl(:, i_soil)) * options%dtime / 1000.0_lp
      sb_spinup_accelerator_t_soil_sum(:, profile_entry) = &
        & sb_spinup_accelerator_t_soil_sum(:, profile_entry) + t_soil_sl(:, i_soil) * options%dtime
      sb_spinup_accelerator_wtr_soil_sum(:, profile_entry) = &
        & sb_spinup_accelerator_wtr_soil_sum(:, profile_entry) + wtr_soil_sl(:, i_soil) * options%dtime
      sb_spinup_accelerator_ice_soil_sum(:, profile_entry) = &
        & sb_spinup_accelerator_ice_soil_sum(:, profile_entry) + ice_soil_sl(:, i_soil) * options%dtime
      sb_spinup_accelerator_wtr_pot_sum(:, profile_entry) = &
        & sb_spinup_accelerator_wtr_pot_sum(:, profile_entry) + wtr_soil_pot_sl(:, i_soil) * options%dtime
    END DO

    ! Post-pool-update uptake is the realized, demand-limited flux. Litter is a
    ! timestep increment, uptake is a volumetric micro-mol/s rate.
    CALL accumulate_replay_nutrient_budget(tile, options, &
      & sb_spinup_accelerator_litter_n_sum, sb_spinup_accelerator_litter_p_sum, &
      & sb_spinup_accelerator_plant_uptake_n_sum, sb_spinup_accelerator_plant_uptake_p_sum)

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
    ! An event is atomic for a tile chunk. A mixed restart must not rerun only
    ! some cells implicitly, nor silently replay already completed cells.
    IF (ANY(NINT(sb_spinup_accelerator_execution_count) /= NINT(sb_spinup_accelerator_execution_count(1)))) THEN
      CALL finish(modname//':check_sb_spinup_accelerator_schedule', 'Inconsistent replay completion state within chunk.')
    END IF
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
    IF (model%config%l_spinup) sb__mem%convergence_event_count%ptr(ics:ice,iblk) = &
      & sb__mem%convergence_event_count%ptr(ics:ice,iblk) + 1._lp

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
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_litter_n_sum
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_litter_p_sum
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_plant_uptake_n_sum
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_plant_uptake_p_sum
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_forcing_event
    dsl4jsb_Real3D_onChunk :: sb_accel_cue_sum_sl
    dsl4jsb_Real3D_onChunk :: sb_accel_poly_sum_sl
    dsl4jsb_Real3D_onChunk :: sb_accel_ap_sum_sl
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_daily_time
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
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_litter_n_sum)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_litter_p_sum)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_plant_uptake_n_sum)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_plant_uptake_p_sum)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_forcing_event)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_accel_cue_sum_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_accel_poly_sum_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_accel_ap_sum_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_daily_time)
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

    sb_spinup_accelerator_litter_n_sum(:) = 0.0_lp
    sb_spinup_accelerator_litter_p_sum(:) = 0.0_lp
    sb_spinup_accelerator_plant_uptake_n_sum(:) = 0.0_lp
    sb_spinup_accelerator_plant_uptake_p_sum(:) = 0.0_lp
    sb_spinup_accelerator_forcing_event(:) = 0.0_lp
    sb_accel_cue_sum_sl = 0._lp
    sb_accel_poly_sum_sl = 0._lp
    sb_accel_ap_sum_sl = 0._lp
    sb_spinup_accelerator_daily_time(:,:) = 0.0_lp
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
  !> Resume an event at the restart boundary before any coupled process runs.
  !! The host supplies the unchanged checkpoint datetime and chunk. Completed
  !! events are skipped by the persistent counter; incomplete forcing is fatal.
  SUBROUTINE resume_sb_spinup_accelerator(model, options)
    USE mo_jsb_model_class, ONLY: MODEL_QUINCY
    USE mo_jsb_process_class, ONLY: SB_
    USE mo_jsb_time, ONLY: is_newyear
    USE mo_sb_config_class, ONLY: t_sb_config
    TYPE(t_jsb_model), POINTER, INTENT(in) :: model
    TYPE(t_jsb_task_options), INTENT(in) :: options
    TYPE(t_jsb_task_options) :: cell_options
    INTEGER :: c
    IF (model%config%model_scheme /= MODEL_QUINCY .OR. .NOT. model%config%l_spinup) RETURN
    IF (.NOT. is_newyear(options%current_datetime,options%dtime)) RETURN
    IF (.NOT. model%processes(SB_)%p%config%active) RETURN
    SELECT TYPE (conf => model%processes(SB_)%p%config)
    TYPE IS (t_sb_config)
      IF (.NOT. conf%flag_sb_model_spinup_accelerator) RETURN
      IF (conf%sb_convergence_cycle_years <= 0) &
        & CALL finish(modname, 'Grid replay convergence requires a positive configured climate forcing cycle.')
    END SELECT
    cell_options=options
    cell_options%nc=1
    DO c=options%ics,options%ice
      cell_options%ics=c
      cell_options%ice=c
      CALL execute_sb_grid_replay(model,cell_options)
    END DO
  END SUBROUTINE resume_sb_spinup_accelerator

  !> Run all positive-area SB tiles in a cell for one artificial year at a time.
  !! Completion depends only on the area-weighted, depth-integrated grid SOC,
  !! not on the number of tiles, MPI decomposition or nproma grouping. Each cell
  !! may finish its local replay independently; its GLOBAL timestamp never moves.
  !! Live forcing is restored after every tile/year because physics memory may
  !! be inherited/shared between PFTs. The driver cannot checkpoint inside this call.
  SUBROUTINE execute_sb_grid_replay(model, options)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    USE mo_sb_replay_state, ONLY: t_sb_replay_checkpoint, capture_sb_replay_state, &
      & restore_sb_replay_state, release_sb_replay_state
    USE mo_jsb_process_class, ONLY: SB_, HYDRO_
    USE mo_hydro_interface, ONLY: update_soil_properties
    USE mo_q_sb_rate_modifier, ONLY: update_sb_rate_modifier
    USE mo_sb_config_class, ONLY: t_sb_config
    USE mo_sb_util, ONLY: reset_sb_fluxes
    USE mo_sb_convergence_kernel, ONLY: update_replay_convergence
    USE mo_jsb_time, ONLY: executed_sb_model_spinup_accelerator_this_timestep, get_year
    TYPE(t_jsb_model), POINTER, INTENT(in) :: model
    TYPE(t_jsb_task_options), INTENT(in) :: options
    TYPE(t_sb_config), POINTER :: conf
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    TYPE(t_sb_grid_replay_tile), ALLOCATABLE :: context(:)
    TYPE(t_sb_spinup_forcing_climatology) :: shared_forcing !< Shared event forcing, per vegetated area
    TYPE(t_sb_replay_checkpoint) :: checkpoint !< Cell-local pre-event state, including shared owners
    LOGICAL, ALLOCATABLE :: selected(:) !< Positive-area pending SB participants [-]
    REAL(lp), ALLOCATABLE :: history(:,:), after(:,:,:,:)
    REAL(lp) :: stocks(5), initial(5), previous(5), relative(5), fraction, tolerance
    REAL(lp) :: clock_start, clock_end
    INTEGER :: it, event, index, n, y, consecutive
    LOGICAL :: pending, converged, balanced, any_pending
    LOGICAL :: shared_ready, combined_done !< Shared preparation and validated execution
    CHARACTER(len=512) :: shared_reason !< One diagnostic per event, including field values and bounds
    CHARACTER(len=1024) :: text
    CHARACTER(len=32) :: tile_action !< Event-level integration decision, not a yearly trace
    SELECT TYPE (c => model%processes(SB_)%p%config)
    TYPE IS (t_sb_config)
      conf => c
    END SELECT
    ALLOCATE(context(SIZE(model%tiles)))
    any_pending=.FALSE.; event=0
    ! First audit completion across ALL tiles before mutating any of them.
    combined_done=.FALSE.
    DO it=1,SIZE(model%tiles)
      tile => model%tiles(it)%p
      IF (.NOT. tile%Is_process_calculated(SB_)) CYCLE
      IF (tile%lcts(1)%lib_id == 0) CYCLE
      CALL check_sb_spinup_accelerator_schedule(tile,options,pending,index)
      context(it)%pending=pending
      IF (.NOT. pending) CYCLE
      IF (event /= 0 .AND. event /= index) CALL finish(modname, 'Inconsistent pending grid replay event indices.')
      event=index
      any_pending=.TRUE.
    END DO
    IF (.NOT. any_pending) RETURN
    DO it=1,SIZE(model%tiles)
      tile => model%tiles(it)%p
      IF (.NOT. tile%Is_process_calculated(SB_)) CYCLE
      IF (tile%lcts(1)%lib_id == 0) CYCLE
      IF (.NOT. context(it)%pending) CALL finish(modname, 'Mixed completed/pending tile restart for grid replay.')
    END DO
    IF (model%config%flag_slow_sb_pool_spinup_accelerator) &
      & CALL finish(modname, 'Forward and conventional slow-pool acceleration cannot run together.')
    ! Capture all participants together BEFORE refreshing shared physics or seeding
    ! any tile. A later candidate rejection must restore the same pre-event state
    ! for every PFT, including forcing accumulators and both BGCM representations.
    ALLOCATE(selected(SIZE(model%tiles)))
    selected = .FALSE.
    DO it = 1, SIZE(model%tiles)
      IF (.NOT. context(it)%pending) CYCLE
      selected(it) = model%tiles(it)%p%fract(options%ics,options%iblk) > 0._lp
    END DO
    CALL capture_sb_replay_state(model,options,selected,checkpoint)
    ! Replay precedes the coupled task queue, including soil_properties. Its derived
    ! field capacity and porosity are not restart variables: after restart they
    ! still contain initialization values (potentially from 3D inputs), whereas a
    ! continuous run has the last coupled-task values. Rebuild them with the SAME
    ! coupled routine before any replay/analytical rate calculation. Visit actual
    ! HYDRO owners, not the PFTs inheriting their memory. No water, pools, forcing
    ! integrals, lagged biology or global datetime are advanced by this refresh.
    DO it=1,SIZE(model%tiles)
      tile => model%tiles(it)%p
      IF (tile%Is_process_calculated(HYDRO_)) CALL update_soil_properties(tile,options)
    END DO
    n=conf%sb_convergence_cycle_years
    tolerance=conf%sb_model_spinup_accelerator_convergence_rel_tol
    ALLOCATE(history(5,0:n))
    CALL CPU_TIME(clock_start)
    CALL grid_replay_stocks(model,options,initial,fraction)
    DO it=1,SIZE(model%tiles)
      IF (.NOT. context(it)%pending) CYCLE
      tile => model%tiles(it)%p
      ! Report once per pending event, before zero-area tiles are skipped.
      ! Fractions are absolute grid-cell fractions, not renormalized PFT weights.
      IF (model%config%flag_spinup_debug) THEN
        tile_action='participating in event'
        IF (tile%fract(options%ics,options%iblk) <= 0._lp) tile_action='skipped: nonpositive area'
        WRITE(text,*) 'Replay tile decision: model/year/event/ic/iblk=',model%id, &
          & get_year(options%current_datetime),event,options%ics,options%iblk, &
          & '; tile=',TRIM(tile%name),'; grid area fraction=',tile%fract(options%ics,options%iblk), &
          & '; action=',TRIM(tile_action)
        CALL message(modname,TRIM(text),all_print=.TRUE.)
      END IF
      IF (tile%fract(options%ics,options%iblk) <= 0._lp) CYCLE
      context(it)%active=.TRUE.
      ! In particular rmm_sorption is diagnostic, not restart state, and is used
      ! when the analytical seed refreshes nutrient sorption capacities. Rebuild
      ! it from live physics before preparing the seed; replay steps subsequently
      ! recompute all modifiers from each loaded climatology bin as before.
      CALL update_sb_rate_modifier(tile,options)
      ! Restart fields have now replaced cold-start geometry. Validate before preparing any
      ! analytical target, including per-PFT events; refresh never rescales inventories.
      IF (conf%sb_model_scheme == 'jsm') &
        & CALL calibrate_sb_pore_volume(tile,options,0,l_refresh_only=.TRUE.)
      IF (model%config%flag_spinup_debug) CALL write_sb_replay_entry_debug(tile,options,event)
      IF (model%config%flag_spinup_debug) CALL snapshot_sb_spinup_pools(tile,options,context(it)%before)
      CALL prepare_sb_spinup_accelerator(tile,options,context(it)%forcing,context(it)%live)
    END DO
    ! Preparation is a distinct read-only phase. All PFT forcing/live-state
    ! snapshots now exist before any analytical target changes soil pools.
    ! Construct combined inputs before any seed changes the snapshot. The
    ! transaction below leases matching shared state, never a PFT state with
    ! just its forcing replaced by a mean.
    ! Explicit per-PFT selection bypasses the shared transaction, not the common
    ! snapshots, event schedule or convergence accounting. No restart state changes.
    shared_ready=.FALSE.
    shared_reason='Per-PFT replay explicitly selected.'
    WRITE(text,*) 'SB replay requested mode: event/cell/block=',event,options%ics,options%iblk, &
      & '; mode=',TRIM(conf%sb_model_spinup_accelerator_replay_mode)
    CALL message(modname,TRIM(text),all_print=.TRUE.)
    IF (conf%sb_model_spinup_accelerator_replay_mode == 'combined') &
      & CALL prepare_combined_sb_forcing(model,options,context,shared_forcing,shared_ready,shared_reason)
    IF (model%config%flag_spinup_debug) THEN
      WRITE(text,*) 'Combined replay forcing: event/cell/block/ready=',event,options%ics,options%iblk, &
        & shared_ready,'; ',TRIM(shared_reason)
      CALL message(modname,TRIM(text),all_print=.TRUE.)
    END IF
    IF (shared_ready) THEN
      CALL try_combined_sb_replay(model,options,context,shared_forcing,checkpoint,event, &
        & stocks,previous,relative,y,converged,combined_done,shared_reason)
      WRITE(text,*) 'Combined replay outcome: event/cell/block/accepted=',event,options%ics,options%iblk, &
        & combined_done,'; ',TRIM(shared_reason)
      CALL message(modname,TRIM(text),all_print=.TRUE.)
      IF (.NOT. combined_done) THEN
        ! Candidate rejection is atomic. Recreate the original derived physics
        ! before the existing per-PFT fallback (including restarted entry).
        CALL restore_sb_replay_state(checkpoint)
        DO it=1,SIZE(model%tiles)
          tile => model%tiles(it)%p
          IF (tile%Is_process_calculated(HYDRO_)) CALL update_soil_properties(tile,options)
        END DO
        DO it=1,SIZE(model%tiles)
          IF (.NOT. context(it)%active) CYCLE
          CALL update_sb_rate_modifier(model%tiles(it)%p,options)
        END DO
      END IF
    END IF
    IF (.NOT. combined_done) THEN
      WRITE(text,*) 'SB replay realized mode: per_pft; event/cell/block=',event,options%ics,options%iblk, &
        & '; ',TRIM(shared_reason)
      CALL message(modname,TRIM(text),all_print=.TRUE.)
      DO it=1,SIZE(model%tiles)
        IF (.NOT. context(it)%active) CYCLE
        tile => model%tiles(it)%p
        CALL seed_sb_replay_analytically(tile,options,context(it)%forcing,event)
        ALLOCATE(context(it)%balance%litter_n(1),context(it)%balance%litter_p(1), &
          & context(it)%balance%uptake_n(1),context(it)%balance%uptake_p(1), &
          & context(it)%balance%scale_n(1),context(it)%balance%scale_p(1))
        context(it)%balance%scale_n=1._lp
        context(it)%balance%scale_p=1._lp
        CALL finalize_sb_spinup_accelerator(tile,options,context(it)%live)
        IF (model%config%flag_spinup_debug) THEN
          CALL write_sb_spinup_forcing_debug(tile,options,event,context(it)%forcing)
          CALL write_sb_spinup_snapshot_debug(tile,options,event,0)
        END IF
      END DO
      ! The optional analytical seed is outside the replay convergence window.
      CALL grid_replay_stocks(model,options,stocks,fraction)
      history=0._lp
      history(:,0)=stocks
      consecutive=0; converged=.FALSE.; y=0
      DO index=1,conf%sb_model_spinup_accelerator_max_relaxation_years
        IF (fraction <= 0._lp) EXIT
        y=index
        DO it=1,SIZE(model%tiles)
          IF (.NOT. context(it)%active) CYCLE
          tile => model%tiles(it)%p
          SELECT CASE(TRIM(conf%sb_model_scheme))
          CASE('jsm')
            CALL run_jsm_forward_relaxation(tile,options,context(it)%forcing,context(it)%balance)
            CALL calibrate_sb_pore_volume(tile,options,y)
          CASE('simple_1d')
            CALL run_simple_1d_forward_relaxation(tile,options,context(it)%forcing,context(it)%balance)
          CASE DEFAULT
            CALL finish(modname, 'Unsupported grid replay soil model.')
          END SELECT
          CALL update_replay_plant_np_demand_balance(tile,options,context(it)%forcing, &
            & context(it)%balance,event,y,balanced)
          ! Nutrient balance still controls next-year demand, but does not veto SOC.
          CALL finalize_sb_spinup_accelerator(tile,options,context(it)%live)
        END DO
        CALL grid_replay_stocks(model,options,stocks,fraction)
        CALL update_replay_convergence(conf%sb_model_spinup_accelerator_convergence_mode, &
          & n,y,tolerance,stocks,history,consecutive,converged,previous,relative)
        IF (model%config%flag_spinup_debug .AND. is_sb_spinup_debug_year(y)) THEN
          DO it=1,SIZE(model%tiles)
            IF (.NOT. context(it)%active) CYCLE
            tile => model%tiles(it)%p
            CALL write_sb_spinup_snapshot_debug(tile,options,event,y)
          END DO
        END IF
        IF (converged) EXIT
      END DO
    END IF
    ! Do not mark a corrupt candidate completed or discard its forcing history.
    ! This restores state before a controlled fatal exit, not after a kernel abort.
    IF (.NOT. ALL(ieee_is_finite(stocks))) THEN
      CALL restore_sb_replay_state(checkpoint)
      CALL release_sb_replay_state(checkpoint)
      CALL finish(modname, 'Non-finite replay inventory; pre-event state restored, event not committed.')
    END IF
    DO it=1,SIZE(model%tiles)
      IF (.NOT. context(it)%pending) CYCLE
      tile => model%tiles(it)%p
      IF (context(it)%active) THEN
        CALL reset_sb_fluxes(tile,options)
        CALL finalize_sb_spinup_accelerator(tile,options,context(it)%live)
        IF (model%config%flag_spinup_debug) THEN
          CALL write_sb_spinup_snapshot_debug(tile,options,event,y)
          CALL snapshot_sb_spinup_pools(tile,options,after)
          CALL write_sb_spinup_debug(tile,options,event,y,context(it)%before,after)
          DEALLOCATE(after)
        END IF
      END IF
      CALL reset_sb_spinup_accelerator_forcing(tile,options)
      CALL mark_sb_spinup_accelerator_execution(tile,options,event)
      CALL tile%bgcm_store%Write_matrices_to_bgcm_vars(tile%name,options%ics,options%ice,options%iblk)
    END DO
    CALL release_sb_replay_state(checkpoint)
    CALL CPU_TIME(clock_end)
    executed_sb_model_spinup_accelerator_this_timestep=.TRUE.
    WRITE(text,*) 'Grid replay: model/year/event/ic/iblk=',model%id,get_year(options%current_datetime), &
      & event,options%ics,options%iblk,'; mode/N/tolerance=',TRIM(conf%sb_model_spinup_accelerator_convergence_mode), &
      & n,tolerance,'; internal years/SOC pass=',y,converged,'; before/after C,N,P,SIN,SIP mol/m2=',initial,stocks, &
      & '; CPU seconds=',clock_end-clock_start,'; returning to unchanged global timestamp.'
    CALL message(modname,TRIM(text),all_print=.TRUE.)
    IF (y > 0) THEN
      WRITE(text,*) 'Grid replay final diagnostic: signed stock changes=',stocks-previous, &
        & '; relative changes C,N,P,SIN,SIP=',relative, &
        & '; cycle-window comparison valid=',y >= n,'; nutrient diagnostics do not gate SOC.'
      CALL message(modname,TRIM(text),all_print=.TRUE.)
    END IF
  END SUBROUTINE execute_sb_grid_replay

  !> Lease one existing SB leaf as scratch execution storage, then restore ALL
  !> owners and stage the complete redistribution. No parent diagnostic tile is
  !> integrated. Fractions, process associations and the host clock never change.
  !> A normal rejection restores the checkpoint; the caller can replay per PFT.
  SUBROUTINE try_combined_sb_replay(model,options,context,forcing,checkpoint,event, &
      & stocks,previous,relative,years,converged,accepted,reason)
    USE mo_sb_constants, ONLY: min_fine_mineral_volume
    USE mo_sb_combined_state, ONLY: t_sb_combined_state, read_sb_combined_state, &
      & write_sb_combined_state, reconstruct_sb_combined_state, &
      & redistribute_sb_combined_state, valid_sb_combined_inventory
    USE mo_sb_replay_state, ONLY: t_sb_replay_checkpoint, restore_sb_replay_state
    USE mo_sb_convergence_kernel, ONLY: update_replay_convergence
    USE mo_jsb_process_class, ONLY: SB_, HYDRO_, SSE_
    USE mo_sb_config_class, ONLY: t_sb_config
    USE mo_sb_memory_class, ONLY: t_sb_memory
    USE mo_hydro_memory_class, ONLY: t_hydro_memory
    USE mo_sse_memory_class, ONLY: t_sse_memory
    USE mo_q_sb_rate_modifier, ONLY: update_sb_rate_modifier
    USE mo_q_sb_jsm_processes, ONLY: solve_pore_volume_for_bulk_density, calc_available_pore_carbon_components, &
      & sb_density_solution_found, sb_density_target_above_maximum, sb_density_target_below_minimum
    USE mo_lnd_bgcm_idx, ONLY: ix_dom_assoc, ix_residue_assoc, ixC
    USE mo_lnd_bgcm_store_class, ONLY: SB_BGCM_POOL_ID
    TYPE(t_jsb_model), POINTER, INTENT(in) :: model
    TYPE(t_jsb_task_options), INTENT(in) :: options
    TYPE(t_sb_grid_replay_tile), INTENT(in) :: context(:)
    TYPE(t_sb_spinup_forcing_climatology), INTENT(in) :: forcing
    TYPE(t_sb_replay_checkpoint), INTENT(in) :: checkpoint
    INTEGER, INTENT(in) :: event
    REAL(lp), INTENT(out) :: stocks(5),previous(5),relative(5)
    INTEGER, INTENT(out) :: years
    LOGICAL, INTENT(out) :: converged,accepted
    CHARACTER(len=*), INTENT(out) :: reason
    TYPE(t_sb_combined_state), ALLOCATABLE :: before(:),after(:)
    TYPE(t_sb_combined_state) :: shared,target,check
    TYPE(t_sb_replay_nutrient_balance) :: balance
    TYPE(t_sb_config), POINTER :: conf
    TYPE(t_sb_memory), POINTER :: soil
    TYPE(t_hydro_memory), POINTER :: hydro, reference_hydro
    TYPE(t_sse_memory), POINTER :: sse, reference_sse
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    REAL(lp), POINTER :: pool(:,:,:,:)
    REAL(lp), ALLOCATABLE :: depth(:,:),fraction(:),lower(:,:),upper(:,:),history(:,:),bulk(:,:)
    REAL(lp) :: area,tolerance,maoc,poc,lo,hi,solved,density,residual,weighted(5),target_stocks(5)
    INTEGER, ALLOCATABLE :: members(:)
    INTEGER :: j,it,rep,ns,nt,c,b,n,consecutive,s,status
    INTEGER :: failed_member !< Input PFT index for reconstruction rejection; zero means shared state
    LOGICAL :: valid,balanced
    CHARACTER(len=256) :: report !< Sparse execution-owner diagnostic

    accepted = .FALSE.
    converged = .FALSE.
    years = 0
    stocks = 0._lp
    previous = 0._lp
    relative = 0._lp
    SELECT TYPE(cfg => model%processes(SB_)%p%config)
    TYPE IS(t_sb_config)
      conf => cfg
    END SELECT
    ! This lease inventory is JSM-specific. Do not silently borrow JSM fields
    ! for simple_1d, legacy volume or unaggregated mycorrhizal traits.
    reason='Combined replay not applicable: retains per-PFT execution for this configuration.'
    IF (conf%sb_model_scheme /= 'jsm' .OR. .NOT. conf%flag_sb_jsm_available_pore_volume) RETURN
    IF (conf%flag_mycorrhiza .OR. conf%flag_mycorrhiza_org .OR. conf%flag_mycorrhiza_prim) RETURN
    members=PACK([(it,it=1,SIZE(context))],context%active)
    nt=SIZE(members)
    IF (nt == 0) RETURN
    rep = members(1)
    c = options%ics
    b = options%iblk
    ALLOCATE(before(nt))
    DO j=1,nt
      CALL read_sb_combined_state(model%tiles(members(j))%p,options,before(j))
    END DO
    ns=SIZE(before(1)%pore)
    SELECT TYPE(mem => model%tiles(rep)%p%mem(HYDRO_)%p)
    TYPE IS(t_hydro_memory)
      reference_hydro => mem
    END SELECT
    SELECT TYPE(mem => model%tiles(rep)%p%mem(SSE_)%p)
    TYPE IS(t_sse_memory)
      reference_sse => mem
    END SELECT
    ALLOCATE(depth(ns,nt),fraction(nt),lower(ns,nt),upper(ns,nt),bulk(ns,nt))
    DO j=1,nt
      it=members(j)
      SELECT TYPE(mem => model%tiles(it)%p%mem(HYDRO_)%p)
      TYPE IS(t_hydro_memory)
        hydro => mem
      END SELECT
      SELECT TYPE(mem => model%tiles(it)%p%mem(SSE_)%p)
      TYPE IS(t_sse_memory)
        sse => mem
      END SELECT
      depth(:,j)=hydro%soil_depth_sl%ptr(c,1:ns,b)
      ! Each owner's geometry was validated at event entry, before the snapshots were constructed.
      bulk(:,j)=sse%bulk_dens_sl%ptr(c,1:ns,b)
      fraction(j)=model%tiles(it)%p%fract(c,b)
      reason='Combined replay requires common soil-physics properties.'
      IF (ANY(bulk(:,j) /= bulk(:,1))) RETURN
      ! Equal soil properties may have different memory owners. Water and
      ! temperature are supplied by the shared forcing, not borrowed here.
      IF (ANY(sse%clay_sl%ptr(c,1:ns,b) /= reference_sse%clay_sl%ptr(c,1:ns,b))) RETURN
      IF (ANY(sse%silt_sl%ptr(c,1:ns,b) /= reference_sse%silt_sl%ptr(c,1:ns,b))) RETURN
      IF (ANY(hydro%vol_field_cap_sl%ptr(c,1:ns,b) /= reference_hydro%vol_field_cap_sl%ptr(c,1:ns,b))) RETURN
    END DO
    CALL reconstruct_sb_combined_state(fraction,depth,before,shared,valid,report,failed_member)
    IF (.NOT. valid) THEN
      IF (failed_member > 0) THEN
        reason='Reconstruction rejected on '//TRIM(model%tiles(members(failed_member))%p%name)//': '//TRIM(report)
      ELSE
        reason='Shared reconstruction rejected: '//TRIM(report)
      END IF
      RETURN
    END IF
    tile => model%tiles(rep)%p
    IF (model%config%flag_spinup_debug) THEN
      WRITE(report,*) 'Combined replay execution: event/cell/block=',event,c,b, &
        & '; temporary owner=',TRIM(tile%name),'; represented PFTs=',nt, &
        & '; area=',SUM(fraction),'; snapshots refer to the combined column until redistribution.'
      CALL message(modname,TRIM(report),all_print=.TRUE.)
    END IF
    CALL write_sb_combined_state(tile,options,shared)
    ALLOCATE(balance%litter_n(1),balance%litter_p(1),balance%uptake_n(1),balance%uptake_p(1), &
      & balance%scale_n(1),balance%scale_p(1))
    balance%scale_n = 1._lp
    balance%scale_p = 1._lp
    ! Load shared physics BEFORE refreshing capacities and seeding. Each normal
    ! replay substep subsequently loads its own climatology bin as usual.
    CALL load_sb_spinup_forcing_bin(tile,options,forcing,1,1,balance)
    CALL update_sb_rate_modifier(tile,options)
    CALL calibrate_sb_pore_volume(tile,options,0,l_refresh_only=.TRUE.)
    CALL seed_sb_replay_analytically(tile,options,forcing,event)
    CALL grid_replay_stocks(model,options,stocks,area,only_tile=rep)
    n=conf%sb_convergence_cycle_years
    tolerance=conf%sb_model_spinup_accelerator_convergence_rel_tol
    ALLOCATE(history(5,0:n),SOURCE=0._lp)
    history(:,0)=stocks
    consecutive=0
    DO years=1,conf%sb_model_spinup_accelerator_max_relaxation_years
      CALL run_jsm_forward_relaxation(tile,options,forcing,balance)
      CALL calibrate_sb_pore_volume(tile,options,years)
      CALL update_replay_plant_np_demand_balance(tile,options,forcing,balance,event,years,balanced)
      CALL grid_replay_stocks(model,options,stocks,area,only_tile=rep)
      CALL update_replay_convergence(conf%sb_model_spinup_accelerator_convergence_mode, &
        & n,years,tolerance,stocks,history,consecutive,converged,previous,relative)
      IF (model%config%flag_spinup_debug .AND. is_sb_spinup_debug_year(years)) &
        & CALL write_sb_spinup_snapshot_debug(tile,options,event,years)
      IF (converged) EXIT
    END DO
    years=MIN(years,conf%sb_model_spinup_accelerator_max_relaxation_years)
    CALL read_sb_combined_state(tile,options,target)
    target_stocks=stocks
    ! Restore every owner BEFORE writing any individual redistribution. This
    ! restores PFT forcing, independent chemistry and shared physics exactly.
    CALL restore_sb_replay_state(checkpoint)
    lower=0._lp
    upper=1._lp-min_fine_mineral_volume
    CALL redistribute_sb_combined_state(fraction,depth,before,target,lower,upper,after,valid)
    reason='Combined redistribution rejected: parent shares, isotopes or mineral structure.'
    IF (.NOT. valid) RETURN
    DO j=1,nt
      tile => model%tiles(members(j))%p
      CALL write_sb_combined_state(tile,options,after(j))
      CALL update_sb_rate_modifier(tile,options)
      CALL calibrate_sb_pore_volume(tile,options,years,l_refresh_only=.TRUE.)
      SELECT TYPE(mem => tile%mem(SB_)%p)
      TYPE IS(t_sb_memory)
        soil => mem
      END SELECT
      pool => tile%bgcm_store%Get_matrix_from_store_2l_3d(SB_BGCM_POOL_ID,c,c,b,tile%name)
      DO s=1,ns
        CALL calc_available_pore_carbon_components(pool,1,s,maoc,poc)
        CALL solve_pore_volume_for_bulk_density(maoc,poc,after(j)%capacity(1,s), &
          & conf%sb_bulk_density_target_fraction*bulk(s,j),1.e-3_lp,lo,hi,solved,density,residual,status)
        valid=status == sb_density_solution_found .OR. status == sb_density_target_above_maximum .OR. &
          & status == sb_density_target_below_minimum
        valid=valid .AND. after(j)%pore(s) >= lo .AND. after(j)%pore(s) <= hi
        valid=valid .AND. soil%volume_pore_maoc_excess_sl%ptr(c,s,b) <= 1.e-12_lp
        IF (s > 1) valid=valid .AND. soil%volume_pore_poc_excess_sl%ptr(c,s,b) <= 1.e-12_lp
        valid=valid .AND. maoc <= soil%qmax_org%ptr(c,s,b)*(1._lp+1.e-12_lp)
        valid=valid .AND. after(j)%inorganic(2,1,s) <= soil%qmax_nh4%ptr(c,s,b)*(1._lp+1.e-12_lp)
        ! Forward plant/microbial competition exposes total qmax_po4 minus
        ! fast-associated PO4. qmax_slow_po4 is an isotherm component, not an
        ! independent inventory bound on the slow/occluded/primary P pools.
        valid=valid .AND. after(j)%inorganic(7,1,s) <= soil%qmax_po4%ptr(c,s,b)*(1._lp+1.e-12_lp)
        IF (.NOT. valid) THEN
          WRITE(reason,'(A,2I6,A,3ES12.4,A,2ES12.4)') 'Redistribution rejected PFT/layer=',members(j),s, &
            & '; pore/bounds=',after(j)%pore(s),lo,hi, &
            & '; MAOC/POC excess=',soil%volume_pore_maoc_excess_sl%ptr(c,s,b), &
            & soil%volume_pore_poc_excess_sl%ptr(c,s,b)
          CALL restore_sb_replay_state(checkpoint)
          RETURN
        END IF
      END DO
      CALL read_sb_combined_state(tile,options,check)
      IF (.NOT. valid_sb_combined_inventory(check,report)) THEN
        reason='Refreshed inventory rejected on '//TRIM(tile%name)//': '//TRIM(report)
        CALL restore_sb_replay_state(checkpoint)
        RETURN
      END IF
    END DO
    CALL grid_replay_stocks(model,options,weighted,area)
    reason='Combined redistribution stock-conservation rejection.'
    IF (ANY(ABS(weighted-target_stocks) > 1.e-12_lp*MAX(ABS(target_stocks),TINY(1._lp)))) THEN
      CALL restore_sb_replay_state(checkpoint)
      RETURN
    END IF
    stocks=weighted
    accepted=.TRUE.
    reason='Combined JSM replay accepted; all PFT redistributions validated.'
  END SUBROUTINE try_combined_sb_replay

  !> Actual BG soil inventory, depth integrated BEFORE weighting calculated tiles.
  !! Reuse the ordinary diagnostic definition; exclude first-layer litter only.
  !! Called after each tile has restored live physics to avoid inherited forcing leakage.
  SUBROUTINE grid_replay_stocks(model,options,stocks,fraction,only_tile)
    USE mo_jsb_process_class, ONLY: SB_, HYDRO_
    USE mo_sb_memory_class, ONLY: t_sb_memory
    USE mo_hydro_memory_class, ONLY: t_hydro_memory
    USE mo_lnd_bgcm_idx, ONLY: ixC, ixN, ixP
    USE mo_lnd_bgcm_store_class, ONLY: SB_BGCM_POOL_ID
    USE mo_q_sb_update_pools, ONLY: calculate_sb_pool_bg_soil_for_element
    TYPE(t_jsb_model), POINTER, INTENT(in) :: model
    TYPE(t_jsb_task_options), INTENT(in) :: options
    INTEGER, OPTIONAL, INTENT(in) :: only_tile !< Temporary shared-column owner; no area/clock mutation
    REAL(lp), INTENT(out) :: stocks(5) !< Active-soil-area C, N, P, SIN, SIP [mol m-2]
    REAL(lp), INTENT(out) :: fraction !< Sum of calculated soil tile fractions [-]
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    TYPE(t_sb_memory), POINTER :: soil
    TYPE(t_hydro_memory), POINTER :: hydro
    REAL(lp), POINTER :: pool(:,:,:,:)
    REAL(lp), ALLOCATABLE :: layer(:,:)
    REAL(lp) :: f
    INTEGER :: it,e,c,b,nz,elem(3)
    stocks=0._lp; fraction=0._lp
    c=options%ics; b=options%iblk
    elem=[ixC,ixN,ixP]
    DO it=1,SIZE(model%tiles)
      IF (PRESENT(only_tile)) THEN
        IF (it /= only_tile) CYCLE
      END IF
      tile => model%tiles(it)%p
      IF (.NOT. tile%Is_process_calculated(SB_)) CYCLE
      IF (tile%lcts(1)%lib_id == 0) CYCLE
      f=tile%fract(c,b)
      IF (f <= 0._lp) CYCLE
      SELECT TYPE (s => tile%mem(SB_)%p)
      TYPE IS (t_sb_memory)
        soil => s
      END SELECT
      SELECT TYPE (h => tile%mem(HYDRO_)%p)
      TYPE IS (t_hydro_memory)
        hydro => h
      END SELECT
      pool => tile%bgcm_store%Get_matrix_from_store_2l_3d(SB_BGCM_POOL_ID,c,c,b,tile%name)
      nz=SIZE(pool,4)
      ALLOCATE(layer(1,nz))
      DO e=1,3
        CALL calculate_sb_pool_bg_soil_for_element(pool,nz,elem(e),layer)
        stocks(e)=stocks(e)+f*SUM(layer(1,:)*hydro%soil_depth_sl%ptr(c,:,b))
      END DO
      stocks(4)=stocks(4)+f*SUM(hydro%soil_depth_sl%ptr(c,:,b)*(soil%nh4_solute%ptr(c,:,b)+ &
        & soil%nh4_assoc%ptr(c,:,b)+soil%no3_solute%ptr(c,:,b)))
      stocks(5)=stocks(5)+f*SUM(hydro%soil_depth_sl%ptr(c,:,b)*(soil%po4_solute%ptr(c,:,b)+ &
        & soil%po4_assoc_fast%ptr(c,:,b)+soil%po4_assoc_slow%ptr(c,:,b)))
      fraction=fraction+f
      DEALLOCATE(layer)
    END DO
    IF (fraction > 0._lp) stocks=stocks/fraction
  END SUBROUTINE grid_replay_stocks

  SUBROUTINE apply_sb_spinup_accelerator(tile, options)
    USE mo_jsb_process_class, ONLY: SB_
    dsl4jsb_Use_config(SB_)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options
    TYPE(t_jsb_model), POINTER :: model
    LOGICAL :: pending
    INTEGER :: execution_index
    dsl4jsb_Def_config(SB_)
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(SB_)
    IF (.NOT. model%config%l_spinup .OR. .NOT. sb__conf%flag_sb_model_spinup_accelerator) RETURN
    ! The driver must execute all pending grid events before coupled tasks.
    ! Never silently fall back to the old tile/chunk-local replay scheduler.
    CALL check_sb_spinup_accelerator_schedule(tile,options,pending,execution_index)
    IF (pending) CALL finish(modname, 'Pending replay reached coupled tasks; grid replay driver was not called.')
    CALL collect_sb_spinup_accelerator_forcing(tile,options)
  END SUBROUTINE apply_sb_spinup_accelerator

  ! ===================================================================================================== !
  !> Seed an eligible replay from its own climatology; no independent clock or forcing collector.
  !> The first COMPLETED replay is structure conditioning. Subsequent columns are transactional:
  !> reject all layers together on invalid forcing, ambiguous composition or failed C equations.
  SUBROUTINE seed_sb_replay_analytically(tile, options, forcing, execution_index)
    USE mo_jsb_process_class, ONLY: SB_, HYDRO_
    USE mo_jsb_math_constants, ONLY: one_day, eps8
    USE mo_jsb_lctlib_class, ONLY: t_lctlib_element
    USE mo_lnd_bgcm_store_class, ONLY: SB_BGCM_POOL_ID
    USE mo_lnd_bgcm_idx
    USE mo_sb_analytical_accelerator, ONLY: t_sb_analytical_inputs, t_sb_local_target, &
      & build_sb_analytical_inputs, solve_simple_carbon_target, reconstruct_analytical_column, &
      & t_sb_column_solution, prepare_analytical_bioturbation_inputs, diagnose_full_profile_carbon_targets
    USE mo_sb_constants, ONLY: tau_metabolic_litter, tau_structural_litter, tau_woody_litter, &
      & tau_fast, tau_slow, frac_woody2dec_litter, frac_litter2fast_som, frac_fast2slow_som, &
      & frac_slow2fast_som, fc_soluable_max, k_fc_soluable, lc_woody_litter
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(HYDRO_)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< Owning SB tile
    TYPE(t_jsb_task_options), INTENT(in) :: options !< Global timestamp and chunk, unchanged
    TYPE(t_sb_spinup_forcing_climatology), INTENT(in) :: forcing !< Validated shared replay means
    INTEGER, INTENT(in) :: execution_index !< Shared schedule event index [-]
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_lctlib_element), POINTER :: lctlib
    TYPE(t_sb_analytical_inputs) :: inputs
    TYPE(t_sb_local_target), ALLOCATABLE :: targets(:,:)
    TYPE(t_sb_column_solution), ALLOCATABLE :: solutions(:) !< Scratch full-column solver diagnostics
    REAL(lp), ALLOCATABLE :: candidate(:,:,:), increments(:)
    REAL(lp), ALLOCATABLE :: partitioned_litter(:,:,:,:,:) !< One-cell shared litter rates [mol m-3 s-1]
    INTEGER, ALLOCATABLE :: pool_indices(:)
    INTEGER :: ics, ice, iblk, nc, ns, c, s, n, failed_layer, accepted_count, rejected_count
    INTEGER :: failed_pool !< Local analytical pool index; zero if rejection is not pool-specific [-]
    INTEGER :: e, element_indices(6) !< Optional element indices for diagnostic increments [-]
    REAL(lp) :: dt, woody_fraction, woody_c, woody_n, rate(5), transfer(4)
    REAL(lp) :: analytical_cpu_start, analytical_cpu_end !< Local process CPU time, not MPI wall time [s]
    REAL(lp) :: increment_total(6) !< Chunk sum for log only, NOT area-weighted global inventory [mol m-2]
    LOGICAL :: accepted, is_jsm
    CHARACTER(len=192) :: reason, first_rejection
    CHARACTER(len=1024) :: report
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    dsl4jsb_Real2D_onChunk :: num_sl_above_bedrock, sb_spinup_accelerator_execution_count
    dsl4jsb_Real3D_onChunk :: soil_depth_sl, vol_field_cap_sl
    dsl4jsb_Real3D_onChunk :: sb_accel_cue_sum_sl, sb_accel_poly_sum_sl, sb_accel_ap_sum_sl

    model => Get_model(tile%owner_model_id)
    IF (.NOT. model%config%flag_analytical_ss_spinup_accelerator) RETURN
    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    nc = options%nc
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_execution_count)
    IF (ALL(sb_spinup_accelerator_execution_count <= 0._lp)) THEN
      CALL message(modname//':analytical_seed', 'Conditioning replay: analytical seed skipped on '//TRIM(tile%name))
      RETURN
    END IF
    ! The analytical equations omit fungal/priming pathways: retain replay rather than a false solution.
    IF (dsl4jsb_Config(SB_)%flag_mycorrhiza .OR. dsl4jsb_Config(SB_)%flag_mycorrhiza_org .OR. &
      & dsl4jsb_Config(SB_)%flag_mycorrhiza_prim) THEN
      CALL message(modname//':analytical_seed', 'Analytical seed rejected: mycorrhizal pathways; replay unchanged')
      RETURN
    END IF
    dsl4jsb_Get_var2D_onChunk(HYDRO_, num_sl_above_bedrock)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_field_cap_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_accel_cue_sum_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_accel_poly_sum_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_accel_ap_sum_sl)
    sb_pool_mt => tile%bgcm_store%Get_matrix_from_store_2l_3d(SB_BGCM_POOL_ID,ics,ice,iblk,tile%name)
    ns = SIZE(sb_pool_mt,4)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    is_jsm = dsl4jsb_Config(SB_)%sb_model_scheme == 'jsm'
    ! Same rounded substep as run_sb_climatological_year, not the host timestep.
    dt = one_day/REAL(MAX(1,CEILING(1._lp/dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_timestep_days)),lp)
    IF (is_jsm) THEN
      pool_indices = [ix_soluable_litter,ix_woody_litter,ix_polymeric_litter,ix_microbial, &
        & ix_residue,ix_dom,ix_dom_assoc,ix_residue_assoc]
    ELSE
      pool_indices = [ix_soluable_litter,ix_woody_litter,ix_polymeric_litter,ix_microbial,ix_residue]
    END IF
    transfer = [frac_woody2dec_litter,frac_litter2fast_som,frac_fast2slow_som,frac_slow2fast_som]
    accepted_count = 0
    rejected_count = 0
    increment_total = 0._lp
    element_indices = [ixC,ixN,ixP,ixC13,ixC14,ixN15]
    first_rejection = ''
    DO c = 1,nc
      IF (sb_spinup_accelerator_execution_count(c) <= 0._lp) CYCLE
      n = NINT(num_sl_above_bedrock(c))
      IF (n <= 0) CYCLE
      ! Build one column at a time: one invalid cell must not reject its MPI neighbours.
      IF (ALLOCATED(forcing%partitioned_litter_rate)) THEN
        IF (forcing%litter_partition_timestep /= dt) &
          & CALL finish(modname//':analytical_seed', 'Combined litter timestep differs from replay.')
        partitioned_litter = forcing%partitioned_litter_rate(:,:,c:c,:,:)
      END IF
      CALL build_sb_analytical_inputs(1,ns,SIZE(sb_pool_mt,1),SIZE(sb_pool_mt,2),forcing%n_litter_sources, &
        & forcing%n_forcing_bins,forcing%n_daily_bins,dt,lctlib%sla,lctlib%growthform, &
        & dsl4jsb_Config(SB_)%sb_model_scheme,num_sl_above_bedrock(c:c),soil_depth_sl(c:c,:), &
        & vol_field_cap_sl(c:c,:),forcing%forcing_time(c:c,:),forcing%litterfall_rate(c:c,:), &
        & forcing%root_fraction(c:c,:),forcing%t_soil(c:c,:),forcing%wtr_soil(c:c,:),forcing%ice_soil(c:c,:), &
        & forcing%wtr_pot(c:c,:),forcing%wtr_down(c:c,:),forcing%wtr_lat_loss(c:c,:), &
        & sb_accel_cue_sum_sl(c:c,:),sb_accel_poly_sum_sl(c:c,:),sb_accel_ap_sum_sl(c:c,:),inputs,reason, &
        & partitioned_litter_rate=partitioned_litter)
      IF (.NOT. inputs%ready) THEN
        rejected_count = rejected_count+1
        IF (LEN_TRIM(first_rejection) == 0) first_rejection = reason
        IF (model%config%flag_spinup_debug) THEN
          WRITE(report,*) 'Rejected analytical forcing: tile/block/cell=',TRIM(tile%name),iblk,ics+c-1, &
            & '; ',TRIM(reason)
          CALL message(modname//':analytical_seed',TRIM(report),all_print=.TRUE.)
        END IF
        CYCLE
      END IF
      inputs%depth = soil_depth_sl(c:c,:)
      IF (is_jsm) THEN
        inputs%pore_volume = sb__mem%volume_pore_total_sl%ptr(ics+c-1:ics+c-1,:,iblk)
        inputs%qmax_org = sb__mem%qmax_org%ptr(ics+c-1:ics+c-1,:,iblk)
        inputs%bulk_density = sb__mem%bulk_dens_corr_sl%ptr(ics+c-1:ics+c-1,:,iblk)
        CALL prepare_analytical_bioturbation_inputs(1,ns,forcing%n_forcing_bins, &
          & num_sl_above_bedrock(c:c),forcing%root_fraction(c:c,:),inputs,reason)
        IF (.NOT. inputs%bioturbation_ready) THEN
          rejected_count = rejected_count+1
          IF (LEN_TRIM(first_rejection) == 0) first_rejection = reason
          CYCLE
        END IF
        ! Always include analytical bioturbation for JSM's available-pore target.
        ! Monthly nonlinear fluxes retain annual CUE/AP. A failed stage never changes
        ! live state; the existing whole-column CNP/isotope reconstruction below is atomic.
        CALL CPU_TIME(analytical_cpu_start)
        CALL diagnose_full_profile_carbon_targets(inputs,sb_pool_mt(pool_indices,ixC,c:c,:), &
          & dsl4jsb_Config(SB_)%sb_spinup_lateral_return_fraction,targets,solutions)
        CALL CPU_TIME(analytical_cpu_end)
        IF (model%config%flag_spinup_debug .OR. .NOT. solutions(1)%admissible) THEN
          WRITE(report,*) 'Analytical profile tile/block/cell=',TRIM(tile%name),iblk,ics+c-1, &
            & '; converged/admissible=',solutions(1)%converged,solutions(1)%admissible, &
            & '; monthly/bioturbation fraction=',solutions(1)%monthly_fraction,solutions(1)%bioturbation_fraction, &
            & '; iterations/evaluations=',solutions(1)%total_iterations,solutions(1)%evaluations, &
            & '; worst layer/equation=',solutions(1)%worst_layer,solutions(1)%worst_equation, &
            & '; local CPU seconds=',analytical_cpu_end-analytical_cpu_start, &
            & '; ',TRIM(solutions(1)%reason)
          CALL message(modname//':analytical_seed',TRIM(report),all_print=.TRUE.)
          IF (solutions(1)%diagnostics%valid) THEN
            WRITE(report,*) 'Analytical profile input/respiration/export/storage/closure [mol C m-2 s-1]=', &
              & solutions(1)%diagnostics%column_budget(:,1)
            CALL message(modname//':analytical_seed',TRIM(report),all_print=.TRUE.)
          END IF
        END IF
      ELSE
        IF (ALLOCATED(targets)) DEALLOCATE(targets)
        ALLOCATE(targets(1,ns))
        DO s = 1,n
          rate = inputs%rate_product_mean(1,s,6) / &
            & [tau_metabolic_litter,tau_woody_litter,tau_structural_litter,tau_fast,tau_slow]
          woody_c = sb_pool_mt(ix_woody_litter,ixC,c,s)
          woody_n = sb_pool_mt(ix_woody_litter,ixN,c,s)
          woody_fraction = 0._lp
          ! Match the forward per-timestep N-loss threshold at the retained woody C:N.
          IF (woody_c > 0._lp .AND. woody_n > 0._lp) THEN
            IF (inputs%litter_mean(ix_woody_litter,ixC,1,s)*dt*(woody_n/woody_c) > eps8) &
              & woody_fraction = MAX(0._lp,fc_soluable_max-k_fc_soluable*lc_woody_litter*woody_c/woody_n)
          END IF
          CALL solve_simple_carbon_target(inputs%litter_mean(pool_indices(1:3),ixC,1,s), &
            & rate,transfer,woody_fraction,dt,targets(1,s))
        END DO
      END IF
      IF (ALLOCATED(candidate)) DEALLOCATE(candidate,increments)
      ALLOCATE(candidate(SIZE(sb_pool_mt,1),SIZE(sb_pool_mt,2),n),increments(SIZE(sb_pool_mt,2)))
      CALL reconstruct_analytical_column(sb_pool_mt(:,:,c,1:n),targets(1,1:n),soil_depth_sl(c,1:n), &
        & pool_indices,ixC,candidate,increments,accepted,reason,failed_layer,failed_pool)
      IF (.NOT. accepted) THEN
        rejected_count = rejected_count+1
        DO s = 1,n
          IF (targets(1,s)%status /= 0) THEN
            failed_layer = s
            reason = targets(1,s)%reason
            EXIT
          END IF
        END DO
        IF (LEN_TRIM(first_rejection) == 0) first_rejection = reason
        IF (model%config%flag_spinup_debug) THEN
          WRITE(report,*) 'Rejected analytical column: tile/block/cell/layer=',TRIM(tile%name),iblk,ics+c-1, &
            & failed_layer,'; ',TRIM(reason),'; layer=0 means not localized'
          CALL message(modname//':analytical_seed',TRIM(report),all_print=.TRUE.)
          IF (failed_layer > 0 .AND. failed_pool > 0) THEN
            WRITE(report,*) 'Rejected reconstruction: local pool/BGCM pool=',failed_pool,pool_indices(failed_pool), &
              & '; C_before/C_target [mol C m-3]=',sb_pool_mt(pool_indices(failed_pool),ixC,c,failed_layer), &
              & targets(1,failed_layer)%carbon(failed_pool)
            CALL message(modname//':analytical_seed',TRIM(report),all_print=.TRUE.)
          END IF
        END IF
        CYCLE
      END IF
      sb_pool_mt(:,:,c,1:n) = candidate
      accepted_count = accepted_count+1
      DO e = 1,SIZE(element_indices)
        IF (element_indices(e) < 1 .OR. element_indices(e) > SIZE(increments)) CYCLE
        increment_total(e) = increment_total(e) + increments(element_indices(e))
      END DO
      IF (is_jsm) THEN
        ! Seed carbon allocation only; preserve AP and nutrient-allocation histories.
        DO s = 1,n
          sb__mem%enzyme_frac_poly%ptr(ics+c-1,s,iblk) = targets(1,s)%alpha_poly
          sb__mem%enzyme_frac_residue%ptr(ics+c-1,s,iblk) = targets(1,s)%alpha_residue
          sb__mem%enzyme_frac_poly_c%ptr(ics+c-1,s,iblk) = targets(1,s)%alpha_poly
          sb__mem%enzyme_frac_poly_c_mavg%ptr(ics+c-1,s,iblk) = targets(1,s)%alpha_poly
        END DO
      END IF
    END DO
    IF (accepted_count > 0) THEN
      IF (is_jsm) CALL calibrate_sb_pore_volume(tile,options,0,l_refresh_only=.TRUE.)
      CALL refresh_seed_pool_diagnostics(tile,options)
      CALL tile%bgcm_store%Write_matrix_to_2l_3d_bgcm_vars(SB_BGCM_POOL_ID,ics,ice,iblk,tile%name,sb_pool_mt)
    END IF
    WRITE(report,*) 'Analytical seed event/tile/block=',execution_index,TRIM(tile%name),iblk, &
      & '; accepted/rejected columns=',accepted_count,rejected_count, &
      & '; chunk-summed artificial C,N,P,C13,C14,N15 increments [mol m-2; C14: pmol m-2]=',increment_total
    CALL message(modname//':analytical_seed',TRIM(report),all_print=.TRUE.)
    IF (rejected_count > 0) CALL message(modname//':analytical_seed',TRIM(first_rejection),all_print=.TRUE.)
    ! This is an artificial carbon seed, not organic N/P or C14 equilibrium. Replay ages C14 normally.
  END SUBROUTINE seed_sb_replay_analytically

  !> Refresh diagnostic views after an artificial seed, without advancing any physical process.
  !> Inorganic pools and lagged nutrient status are unchanged; replay updates these normally.
  SUBROUTINE refresh_seed_pool_diagnostics(tile, options)
    USE mo_jsb_process_class, ONLY: SB_, HYDRO_
    USE mo_lnd_bgcm_store_class, ONLY: SB_BGCM_POOL_ID
    USE mo_lnd_bgcm_idx, ONLY: ixC, ixN, ixP, ixC13, ixC14, ixN15, ix_woody_litter
    USE mo_q_sb_update_pools, ONLY: calculate_sb_pool_bg_soil_for_element, calculate_sb_pool_ag_litter_for_element
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(HYDRO_)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    INTEGER :: ics, ice, iblk, ns, j, e, indices(6)
    REAL(lp), POINTER :: total(:,:), below(:,:), above(:), woody(:) !< Pool diagnostics [mol m-3/m-2]
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)
    sb_pool_mt => tile%bgcm_store%Get_matrix_from_store_2l_3d(SB_BGCM_POOL_ID,ics,ice,iblk,tile%name)
    ns = SIZE(sb_pool_mt,4)
    indices = [ixC,ixN,ixP,ixC13,ixC14,ixN15]
    DO j = 1,6
      e = indices(j)
      IF (e < 1 .OR. e > SIZE(sb_pool_mt,2)) CYCLE
      SELECT CASE (j)
      CASE (1)
        total => sb__mem%sb_pool_total_c%ptr(ics:ice,:,iblk)
        below => sb__mem%sb_pool_total_bg_soil_c%ptr(ics:ice,:,iblk)
        above => sb__mem%sb_pool_total_ag_litter_c%ptr(ics:ice,iblk)
        woody => sb__mem%sb_pool_ag_woody_litter_c%ptr(ics:ice,iblk)
      CASE (2)
        total => sb__mem%sb_pool_total_n%ptr(ics:ice,:,iblk)
        below => sb__mem%sb_pool_total_bg_soil_n%ptr(ics:ice,:,iblk)
        above => sb__mem%sb_pool_total_ag_litter_n%ptr(ics:ice,iblk)
        woody => sb__mem%sb_pool_ag_woody_litter_n%ptr(ics:ice,iblk)
      CASE (3)
        total => sb__mem%sb_pool_total_p%ptr(ics:ice,:,iblk)
        below => sb__mem%sb_pool_total_bg_soil_p%ptr(ics:ice,:,iblk)
        above => sb__mem%sb_pool_total_ag_litter_p%ptr(ics:ice,iblk)
        woody => sb__mem%sb_pool_ag_woody_litter_p%ptr(ics:ice,iblk)
      CASE (4)
        total => sb__mem%sb_pool_total_c13%ptr(ics:ice,:,iblk)
        below => sb__mem%sb_pool_total_bg_soil_c13%ptr(ics:ice,:,iblk)
        above => sb__mem%sb_pool_total_ag_litter_c13%ptr(ics:ice,iblk)
        woody => sb__mem%sb_pool_ag_woody_litter_c13%ptr(ics:ice,iblk)
      CASE (5)
        total => sb__mem%sb_pool_total_c14%ptr(ics:ice,:,iblk)
        below => sb__mem%sb_pool_total_bg_soil_c14%ptr(ics:ice,:,iblk)
        above => sb__mem%sb_pool_total_ag_litter_c14%ptr(ics:ice,iblk)
        woody => sb__mem%sb_pool_ag_woody_litter_c14%ptr(ics:ice,iblk)
      CASE (6)
        total => sb__mem%sb_pool_total_n15%ptr(ics:ice,:,iblk)
        below => sb__mem%sb_pool_total_bg_soil_n15%ptr(ics:ice,:,iblk)
        above => sb__mem%sb_pool_total_ag_litter_n15%ptr(ics:ice,iblk)
        woody => sb__mem%sb_pool_ag_woody_litter_n15%ptr(ics:ice,iblk)
      END SELECT
      total = SUM(sb_pool_mt(:,e,:,:),DIM=1)
      CALL calculate_sb_pool_bg_soil_for_element(sb_pool_mt,ns,e,below)
      CALL calculate_sb_pool_ag_litter_for_element(sb_pool_mt,soil_depth_sl,e,above)
      woody = sb_pool_mt(ix_woody_litter,e,:,1)*soil_depth_sl(:,1)
    END DO
    sb__mem%total_soil_n%ptr(ics:ice,:,iblk) = sb__mem%sb_pool_total_bg_soil_n%ptr(ics:ice,:,iblk) + &
      & sb__mem%nh4_solute%ptr(ics:ice,:,iblk)+sb__mem%nh4_assoc%ptr(ics:ice,:,iblk) + &
      & sb__mem%no3_solute%ptr(ics:ice,:,iblk)
    sb__mem%total_soil_inorg_n%ptr(ics:ice,:,iblk) = sb__mem%nh4_solute%ptr(ics:ice,:,iblk) + &
      & sb__mem%nh4_assoc%ptr(ics:ice,:,iblk)+sb__mem%no3_solute%ptr(ics:ice,:,iblk)
  END SUBROUTINE refresh_seed_pool_diagnostics



  ! ===================================================================================================== !
  !> Validate accumulated forcing and reconstruct the event-local climatology.
  SUBROUTINE prepare_sb_spinup_accelerator(tile, options, forcing, live_forcing)

    USE mo_jsb_process_class, ONLY: SB_, VEG_, HYDRO_, SSE_, A2L_
    USE mo_jsb_time, ONLY: get_year, get_year_length
    USE mo_jsb_math_constants, ONLY: one_day
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    dsl4jsb_Use_config(SB_)
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
    CHARACTER(len=*), PARAMETER :: routine = modname//':prepare_sb_spinup_accelerator'
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    INTEGER :: i_entry
    INTEGER :: forcing_bin
    INTEGER :: n_cells, sampling_year, event_year
    REAL(lp) :: expected_seconds
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_litter_n_sum
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_litter_p_sum
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_plant_uptake_n_sum
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_plant_uptake_p_sum
    dsl4jsb_Real2D_onChunk :: sb_spinup_accelerator_forcing_event
    dsl4jsb_Real3D_onChunk :: sb_spinup_accelerator_daily_time
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
    dsl4jsb_Real3D_onChunk :: growth_p_limit_smoothed_sl
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
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_litter_n_sum)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_litter_p_sum)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_plant_uptake_n_sum)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_plant_uptake_p_sum)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_spinup_accelerator_forcing_event)
    dsl4jsb_Get_var3D_onChunk(SB_, sb_spinup_accelerator_daily_time)
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
    dsl4jsb_Get_var3D_onChunk(SB_, growth_p_limit_smoothed_sl)
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
    forcing%n_daily_bins = SIZE(sb_spinup_accelerator_daily_time, 2)
    forcing%n_litter_sources = SIZE(veg_litterfall_mt, 1)
    forcing%n_elements = SIZE(veg_litterfall_mt, 2)
    IF (forcing%n_forcing_bins <= 0 .OR. &
        & MOD(SIZE(sb_spinup_accelerator_root_fraction_sum, 2), forcing%n_forcing_bins) /= 0) THEN
      CALL finish(TRIM(routine), 'Invalid stored forcing dimensions.')
    END IF
    forcing%n_soil_layers = SIZE(sb_spinup_accelerator_root_fraction_sum, 2) / forcing%n_forcing_bins

    IF (forcing%n_forcing_bins /= sb_spinup_accelerator_n_forcing_bins .OR. &
      & forcing%n_daily_bins /= dsl4jsb_Config(SB_)%sb_spinup_accelerator_n_daily_bins .OR. &
      & SIZE(sb_spinup_accelerator_wtr_down_sum, 2) /= forcing%n_daily_bins * forcing%n_soil_layers .OR. &
      & SIZE(sb_spinup_accelerator_wtr_lat_loss_sum, 2) /= forcing%n_daily_bins * forcing%n_soil_layers) THEN
      CALL finish(TRIM(routine), 'Incompatible replay forcing layout; recollect a complete forcing window.')
    END IF
    IF (ANY(.NOT. ieee_is_finite(sb_spinup_accelerator_forcing_time)) .OR. &
      & ANY(.NOT. ieee_is_finite(sb_spinup_accelerator_daily_time))) THEN
      CALL finish(TRIM(routine), 'Non-finite replay sampling time.')
    END IF
    IF (ANY(sb_spinup_accelerator_forcing_time <= 0.0_lp) .OR. ANY(sb_spinup_accelerator_daily_time <= 0.0_lp)) THEN
      CALL finish(TRIM(routine), 'At least one climatological forcing bin has no sampled model time.')
    END IF

    ! A nonempty bin is not sufficient: a restart halfway through collection
    ! must not masquerade as a complete multi-year climatology. Allow only the
    ! current coupled step as a boundary tolerance, not missing days or years.
    event_year = get_year(options%current_datetime)
    IF (ANY(NINT(sb_spinup_accelerator_forcing_event) /= event_year)) THEN
      CALL finish(TRIM(routine), 'Replay forcing event does not match the scheduled execution year.')
    END IF
    expected_seconds = 0.0_lp
    DO sampling_year = event_year - dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_forcing_years, event_year - 1
      expected_seconds = expected_seconds + REAL(get_year_length(sampling_year), lp) * one_day
    END DO
    IF (ANY(ABS(SUM(sb_spinup_accelerator_forcing_time, DIM=2) - expected_seconds) > options%dtime) .OR. &
      & ANY(ABS(SUM(sb_spinup_accelerator_daily_time, DIM=2) - expected_seconds) > options%dtime)) THEN
      CALL finish(TRIM(routine), 'Incomplete replay forcing window; restart with complete bookkeeping.')
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
      & forcing%growth_p_limit(n_cells, SIZE(sb_spinup_accelerator_growth_p_limit_sum, 2)), &
      & forcing%root_fraction(n_cells, SIZE(sb_spinup_accelerator_root_fraction_sum, 2)), &
      & forcing%wtr_down(n_cells, SIZE(sb_spinup_accelerator_wtr_down_sum, 2)), &
      & forcing%wtr_lat_loss(n_cells, SIZE(sb_spinup_accelerator_wtr_lat_loss_sum, 2)), &
      & forcing%t_soil(n_cells, SIZE(sb_spinup_accelerator_t_soil_sum, 2)), &
      & forcing%wtr_soil(n_cells, SIZE(sb_spinup_accelerator_wtr_soil_sum, 2)), &
      & forcing%ice_soil(n_cells, SIZE(sb_spinup_accelerator_ice_soil_sum, 2)), &
      & forcing%wtr_pot(n_cells, SIZE(sb_spinup_accelerator_wtr_pot_sum, 2)))

    ALLOCATE(forcing%litter_to_uptake_n(n_cells), forcing%litter_to_uptake_p(n_cells))
    forcing%reference_litter_n = sb_spinup_accelerator_litter_n_sum
    forcing%reference_litter_p = sb_spinup_accelerator_litter_p_sum
    forcing%reference_uptake_n = sb_spinup_accelerator_plant_uptake_n_sum
    forcing%reference_uptake_p = sb_spinup_accelerator_plant_uptake_p_sum
    forcing%daily_time = sb_spinup_accelerator_daily_time
    CALL diagnose_replay_reference_ratio(sb_spinup_accelerator_litter_n_sum, &
      & sb_spinup_accelerator_plant_uptake_n_sum, dsl4jsb_Config(SB_)%flag_sb_interactive_n, &
      & forcing%litter_to_uptake_n, 'N')
    CALL diagnose_replay_reference_ratio(sb_spinup_accelerator_litter_p_sum, &
      & sb_spinup_accelerator_plant_uptake_p_sum, dsl4jsb_Config(SB_)%flag_sb_interactive_p, &
      & forcing%litter_to_uptake_p, 'P')

    forcing%forcing_time = sb_spinup_accelerator_forcing_time
    forcing%fine_root_c = sb_spinup_accelerator_fine_root_c_sum / forcing%forcing_time
    forcing%f_n_demand = sb_spinup_accelerator_f_n_demand_sum / forcing%forcing_time
    forcing%f_p_demand = sb_spinup_accelerator_f_p_demand_sum / forcing%forcing_time
    forcing%nhx_deposition = sb_spinup_accelerator_nhx_deposition_sum / forcing%forcing_time
    forcing%noy_deposition = sb_spinup_accelerator_noy_deposition_sum / forcing%forcing_time
    forcing%nhx_n15_deposition = sb_spinup_accelerator_nhx_n15_deposition_sum / forcing%forcing_time
    forcing%noy_n15_deposition = sb_spinup_accelerator_noy_n15_deposition_sum / forcing%forcing_time
    forcing%p_deposition = sb_spinup_accelerator_p_deposition_sum / forcing%forcing_time
    DO i_entry = 1, SIZE(forcing%litterfall_rate, 2)
      forcing_bin = MOD(i_entry - 1, forcing%n_forcing_bins) + 1
      forcing%litterfall_rate(:, i_entry) = sb_spinup_accelerator_litterfall_sum(:, i_entry) / &
        & forcing%forcing_time(:, forcing_bin)
    END DO

    DO i_entry = 1, SIZE(forcing%root_fraction, 2)
      forcing_bin = MOD(i_entry - 1, forcing%n_forcing_bins) + 1
      forcing%growth_p_limit(:, i_entry) = sb_spinup_accelerator_growth_p_limit_sum(:, i_entry) / &
        & forcing%forcing_time(:, forcing_bin)
      forcing%root_fraction(:, i_entry) = sb_spinup_accelerator_root_fraction_sum(:, i_entry) / &
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

    DO i_entry = 1, SIZE(forcing%wtr_down, 2)
      forcing_bin = MOD(i_entry - 1, forcing%n_daily_bins) + 1
      forcing%wtr_down(:, i_entry) = sb_spinup_accelerator_wtr_down_sum(:, i_entry) / &
        & sb_spinup_accelerator_daily_time(:, forcing_bin)
      forcing%wtr_lat_loss(:, i_entry) = sb_spinup_accelerator_wtr_lat_loss_sum(:, i_entry) / &
        & sb_spinup_accelerator_daily_time(:, forcing_bin)
    END DO
    IF (ANY(.NOT. ieee_is_finite(forcing%wtr_down)) .OR. &
      & ANY(.NOT. ieee_is_finite(forcing%wtr_lat_loss))) THEN
      CALL finish(TRIM(routine), 'Non-finite daily transport climatology.')
    END IF

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
    ALLOCATE(live_forcing%growth_p_limit, SOURCE=growth_p_limit_smoothed_sl)
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
  !> One bounded correction after 300 completed coupled experiment years, without replay.
  !> The experiment start survives job restarts. Organic pools and the host clock are unchanged;
  !! mineral-associated inorganic inventories follow the imposed mineral-mass change.
  !> Completion is SB-tile/chunk restart state, never a process-global latch.
  SUBROUTINE apply_sb_coupled_pore_calibration(tile, options)
    USE mo_jsb_process_class, ONLY: SB_
    USE mo_jsb_time, ONLY: get_year, get_year_at_experiment_start, is_newyear
    USE mo_quincy_model_config, ONLY: QLAND
    USE mo_jsb_parallel, ONLY: my_process_is_stdio
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Real2D_onChunk :: sb_coupled_pore_calibration_done
    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: ics, ice, iblk, spinup_age
    INTEGER, PARAMETER :: calibration_age = 300 !< approved one-off conditioning age [yr]
    CHARACTER(len=512) :: text
    CHARACTER(len=*), PARAMETER :: routine = modname//':apply_sb_coupled_pore_calibration'

    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    model => Get_model(tile%owner_model_id)
    IF (.NOT. model%config%l_spinup .OR. model%config%qmodel_id /= QLAND) RETURN
    dsl4jsb_Get_config(SB_)
    IF (dsl4jsb_Config(SB_)%sb_model_scheme /= 'jsm') RETURN
    IF (.NOT. dsl4jsb_Config(SB_)%flag_sb_jsm_available_pore_volume) RETURN
    IF (dsl4jsb_Config(SB_)%flag_sb_model_spinup_accelerator) RETURN
    IF (.NOT. is_newyear(options%current_datetime, options%dtime)) RETURN
    spinup_age = get_year(options%current_datetime) - get_year_at_experiment_start()
    IF (spinup_age /= calibration_age) RETURN
    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_var2D_onChunk(SB_, sb_coupled_pore_calibration_done)
    IF (ALL(sb_coupled_pore_calibration_done == 1.0_lp)) RETURN
    IF (ANY(sb_coupled_pore_calibration_done /= 0.0_lp)) &
      & CALL finish(routine, 'Invalid or mixed coupled pore-calibration completion state.')
    CALL calibrate_sb_pore_volume(tile, options, spinup_age, l_full_correction=.TRUE.)
    sb_coupled_pore_calibration_done = 1.0_lp
    ! Retain history but label windows crossing this one-off structural intervention.
    sb__mem%convergence_event_count%ptr(ics:ice,iblk) = &
      & sb__mem%convergence_event_count%ptr(ics:ice,iblk) + 1._lp
    IF (my_process_is_stdio() .AND. iblk == 1 .AND. ics == 1) THEN
      WRITE(text,*) 'Completed one-off coupled pore calibration: tile, spinup age, calendar year = ', &
        & TRIM(tile%name), spinup_age, get_year(options%current_datetime), &
        & '; bounded correction with mineral-inventory adjustment; completion saved in restart.'
      CALL message(routine, TRIM(text))
    END IF
  END SUBROUTINE apply_sb_coupled_pore_calibration

  !> Bounded annual replay calibration of SB-owned pore structure.
  !!
  !! Input SSE density remains the target reference, never an inout soil property.
  !! MAOC imposes both physical pore and chemical mineral-mass bounds. Only the
  !! structural state changes here; no CNP/isotope mass is created or removed.
  !! Density agreement is diagnostic, not an additional stock-convergence gate.
  SUBROUTINE calibrate_sb_pore_volume(tile, options, represented_year, l_full_correction, l_refresh_only)
    USE mo_jsb_parallel, ONLY: my_process_is_stdio
    USE mo_jsb_process_class, ONLY: SB_, HYDRO_, SSE_
    USE mo_lnd_bgcm_store_class, ONLY: SB_BGCM_POOL_ID
    USE mo_q_sb_jsm_processes, ONLY: solve_pore_volume_for_bulk_density, &
      & calc_available_pore_carbon_components, calc_bulk_soil_carbon_available_pore_volume, &
      & calc_bulk_density_correction_available_pore_volume, calc_qmax_bulk_density_correction, &
      & calc_Psorption_parameter, check_available_pore_geometry, &
      & sb_density_solution_found, sb_density_target_above_maximum, sb_density_target_below_minimum
    USE mo_sb_constants, ONLY: rho_particle_mineral, km_adsorpt_OM_nh4, km_adsorpt_mineral_nh4, &
      & min_fine_mineral_volume
    USE mo_lnd_bgcm_idx, ONLY: ix_woody_litter, ixC
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SSE_)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options
    INTEGER, INTENT(in) :: represented_year !< artificial elapsed years; global clock is unchanged [yr]
    LOGICAL, OPTIONAL, INTENT(in) :: l_full_correction !< full one-off solution rather than replay relaxation [-]
    LOGICAL, OPTIONAL, INTENT(in) :: l_refresh_only !< refresh occupancy/capacity with mineral structure fixed [-]
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: ics, ice, iblk, nc, nsoil_w, ic, is, status
    INTEGER :: bounded_count !< number of layers whose target density is outside admissible bounds [-]
    REAL(lp), PARAMETER :: density_calibration_rel_tol = 1.e-3_lp !< numerical density residual tolerance [-]
    REAL(lp) :: maoc, poc !< carbon concentrations per total soil volume [mol C m-3]
    REAL(lp) :: lower, upper, solved, old_pore !< pore fractions per total soil volume [-]
    REAL(lp) :: density, residual, max_residual, max_change !< density [kg m-3]; normalized residual/fraction [-]
    REAL(lp) :: correction_fraction !< fraction of bounded structural correction applied [-]
    REAL(lp) :: mineral_mass_ratio, new_pore !< new/old mineral mass [-], bounded pore fraction [-]
    REAL(lp) :: old_bound_n, old_bound_p !< pre-calibration mineral-associated inventories [mol m-3]
    REAL(lp) :: inventory_delta_n, inventory_delta_p !< sum of cell-column imposed increments [mol m-2]
    INTEGER :: floor_hold_count !< layers held against further mineral removal [-]
    LOGICAL :: l_report !< sparse replay diagnostics or compulsory one-off summary
    CHARACTER(len=1024) :: calibration_message
    dsl4jsb_Real2D_onChunk :: num_sl_above_bedrock
    dsl4jsb_Real3D_onChunk :: clay_sl, silt_sl
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: nh4_assoc, nh4_n15_assoc
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast, po4_assoc_slow, po4_primary, po4_occluded
    dsl4jsb_Real3D_onChunk :: qmax_nh4_min_sl
    dsl4jsb_Real3D_onChunk :: qmax_nh4
    dsl4jsb_Real3D_onChunk :: km_adsorpt_nh4_sl
    dsl4jsb_Real3D_onChunk :: ph_sl
    dsl4jsb_Real3D_onChunk :: rmm_sorption
    dsl4jsb_Real3D_onChunk :: po4_solute
    dsl4jsb_Real3D_onChunk :: qmax_po4_min_sl
    dsl4jsb_Real3D_onChunk :: qmax_po4_om_sl
    dsl4jsb_Real3D_onChunk :: qmax_po4
    dsl4jsb_Real3D_onChunk :: km_adsorpt_po4_sl
    dsl4jsb_Real3D_onChunk :: qmax_fast_po4
    dsl4jsb_Real3D_onChunk :: qmax_slow_po4
    dsl4jsb_Real3D_onChunk :: km_fast_po4
    dsl4jsb_Real3D_onChunk :: km_slow_po4
    dsl4jsb_Real3D_onChunk :: Qmax_AlFe_cor
    dsl4jsb_Real3D_onChunk :: bulk_dens_sl
    dsl4jsb_Real3D_onChunk :: volume_pore_total_sl
    dsl4jsb_Real3D_onChunk :: volume_pore_maoc_sl
    dsl4jsb_Real3D_onChunk :: volume_pore_poc_sl
    dsl4jsb_Real3D_onChunk :: volume_pore_free_sl
    dsl4jsb_Real3D_onChunk :: volume_pore_maoc_excess_sl
    dsl4jsb_Real3D_onChunk :: volume_pore_poc_excess_sl
    dsl4jsb_Real3D_onChunk :: poc_displaceable_frac_sl
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onChunk :: bulk_soil_carbon_sl
    dsl4jsb_Real3D_onChunk :: soil_litter_carbon_sl
    dsl4jsb_Real3D_onChunk :: volume_min_sl
    dsl4jsb_Real3D_onChunk :: surface_burial_frac
    dsl4jsb_Real3D_onChunk :: qmax_org_min_sl
    dsl4jsb_Real3D_onChunk :: qmax_org

    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    nc = ice - ics + 1
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(SB_)
    IF (.NOT. dsl4jsb_Config(SB_)%flag_sb_jsm_available_pore_volume) RETURN
    correction_fraction = dsl4jsb_Config(SB_)%sb_pore_volume_relaxation
    l_report = model%config%flag_spinup_debug .AND. is_sb_spinup_debug_year(represented_year)
    IF (PRESENT(l_full_correction)) THEN
      IF (l_full_correction) correction_fraction = 1.0_lp
      l_report = l_report .OR. l_full_correction
    END IF
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)
    ! Direct store access keeps ownership on this SB tile, not a shared physics tile.
    sb_pool_mt => tile%bgcm_store%Get_matrix_from_store_2l_3d(SB_BGCM_POOL_ID, ics, ice, iblk, tile%name)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, num_sl_above_bedrock)
    dsl4jsb_Get_var3D_onChunk(SSE_, bulk_dens_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, nh4_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_, nh4_n15_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_, po4_assoc_fast)
    dsl4jsb_Get_var3D_onChunk(SB_, po4_assoc_slow)
    dsl4jsb_Get_var3D_onChunk(SB_, po4_primary)
    dsl4jsb_Get_var3D_onChunk(SB_, po4_occluded)
    dsl4jsb_Get_var3D_onChunk(SB_, volume_pore_total_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, volume_pore_maoc_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, volume_pore_poc_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, volume_pore_free_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, volume_pore_maoc_excess_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, volume_pore_poc_excess_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, poc_displaceable_frac_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, bulk_dens_corr_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, bulk_soil_carbon_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, soil_litter_carbon_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, volume_min_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, surface_burial_frac)
    dsl4jsb_Get_var3D_onChunk(SB_, qmax_org_min_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, qmax_org)
    dsl4jsb_Get_var3D_onChunk(SSE_, clay_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, silt_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, qmax_nh4_min_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, qmax_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_, km_adsorpt_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, ph_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, rmm_sorption)
    dsl4jsb_Get_var3D_onChunk(SB_, po4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_, qmax_po4_min_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, qmax_po4_om_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, qmax_po4)
    dsl4jsb_Get_var3D_onChunk(SB_, km_adsorpt_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, qmax_fast_po4)
    dsl4jsb_Get_var3D_onChunk(SB_, qmax_slow_po4)
    dsl4jsb_Get_var3D_onChunk(SB_, km_fast_po4)
    dsl4jsb_Get_var3D_onChunk(SB_, km_slow_po4)
    dsl4jsb_Get_var3D_onChunk(SB_, Qmax_AlFe_cor)
    nsoil_w = SIZE(sb_pool_mt, 4)
    ! Includes refresh-only calls after restart/reconstruction. Never divide by an invalid old mineral mass.
    CALL check_available_pore_geometry(nc,num_sl_above_bedrock,volume_pore_total_sl,tile%name)
    max_change = 0._lp
    bounded_count = 0
    floor_hold_count = 0
    inventory_delta_n = 0._lp
    inventory_delta_p = 0._lp
    DO ic = 1, nc
      DO is = 1, MIN(nsoil_w, NINT(num_sl_above_bedrock(ic)))
        IF (PRESENT(l_refresh_only)) THEN
          IF (l_refresh_only) CYCLE
        END IF
        CALL calc_available_pore_carbon_components(sb_pool_mt, ic, is, maoc, poc)
        CALL solve_pore_volume_for_bulk_density(maoc, poc, qmax_org_min_sl(ic,is), &
          & dsl4jsb_Config(SB_)%sb_bulk_density_target_fraction * bulk_dens_sl(ic,is), &
          & density_calibration_rel_tol, lower, upper, solved, density, residual, status)
        SELECT CASE (status)
        CASE (sb_density_solution_found, sb_density_target_above_maximum, sb_density_target_below_minimum)
          IF (status /= sb_density_solution_found) bounded_count = bounded_count + 1
          old_pore = volume_pore_total_sl(ic,is)
          new_pore = MIN(upper, MAX(lower, old_pore + correction_fraction * (solved - old_pore)))
          ! No permanent latch: a future target may restore mineral mass (decrease V_pore).
          IF (1._lp-old_pore <= min_fine_mineral_volume+1.e-10_lp .AND. new_pore >= old_pore) THEN
            new_pore = old_pore
            floor_hold_count = floor_hold_count + 1
          END IF
          ! At fixed layer thickness, mineral mass is rho_min*(1-V_pore). Preserve nutrient
          ! concentration per kg mineral during this artificial structure adjustment. This is
          ! NOT weathering/occlusion or biological turnover; solutes and organic pools are unchanged.
          mineral_mass_ratio = (1._lp-new_pore)/(1._lp-old_pore)
          old_bound_n = nh4_assoc(ic,is)
          old_bound_p = po4_assoc_fast(ic,is)+po4_assoc_slow(ic,is)+po4_primary(ic,is)+po4_occluded(ic,is)
          nh4_assoc(ic,is) = nh4_assoc(ic,is)*mineral_mass_ratio
          nh4_n15_assoc(ic,is) = nh4_n15_assoc(ic,is)*mineral_mass_ratio
          po4_assoc_fast(ic,is) = po4_assoc_fast(ic,is)*mineral_mass_ratio
          po4_assoc_slow(ic,is) = po4_assoc_slow(ic,is)*mineral_mass_ratio
          po4_primary(ic,is) = po4_primary(ic,is)*mineral_mass_ratio
          po4_occluded(ic,is) = po4_occluded(ic,is)*mineral_mass_ratio
          inventory_delta_n = inventory_delta_n+old_bound_n*(mineral_mass_ratio-1._lp)*soil_depth_sl(ic,is)
          inventory_delta_p = inventory_delta_p+old_bound_p*(mineral_mass_ratio-1._lp)*soil_depth_sl(ic,is)
          volume_pore_total_sl(ic,is) = new_pore
          max_change = MAX(max_change, ABS(volume_pore_total_sl(ic,is) - old_pore))
        CASE DEFAULT
          WRITE(calibration_message,*) 'Pore calibration has no admissible state: block/cell/layer/status=', &
            & iblk, ic, is, status
          CALL finish(modname//':calibrate_sb_pore_volume', TRIM(calibration_message))
        END SELECT
      END DO
    END DO
    CALL calc_bulk_soil_carbon_available_pore_volume(nc, nsoil_w, num_sl_above_bedrock, 'jsm', sb_pool_mt, &
      & surface_burial_frac, bulk_soil_carbon_sl, soil_litter_carbon_sl, volume_min_sl, &
      & volume_pore_total_sl, volume_pore_maoc_sl, volume_pore_poc_sl, volume_pore_free_sl, &
      & poc_displaceable_frac_sl, volume_pore_maoc_excess_sl, volume_pore_poc_excess_sl)
    bulk_dens_corr_sl = calc_bulk_density_correction_available_pore_volume( &
      & volume_pore_total_sl, volume_pore_maoc_sl, volume_pore_poc_sl)
    CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl, volume_min_sl, bulk_dens_sl, &
      & 0._lp, qmax_org_min_sl, qmax_org, &
      & mineral_mass_density_sl=rho_particle_mineral * (1._lp - volume_pore_total_sl))
    ! Refresh nutrient sorption too: replay must not use previous-year mineral mass.
    ! This also replaces restart-stored phosphate Qmax/Km by the current formulation
    ! before replay entry, without imposing equilibrium on existing inorganic pools.
    CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl, volume_min_sl, bulk_dens_sl, &
      & 0._lp, qmax_nh4_min_sl, qmax_nh4, &
      & mineral_mass_density_sl=rho_particle_mineral * (1._lp - volume_pore_total_sl))
    CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl, volume_min_sl, bulk_dens_sl, &
      & km_adsorpt_OM_nh4, km_adsorpt_mineral_nh4, km_adsorpt_nh4_sl, &
      & mineral_mass_density_sl=rho_particle_mineral * (1._lp - volume_pore_total_sl))
    DO ic = 1, nc
      DO is = 1, MIN(nsoil_w, NINT(num_sl_above_bedrock(ic)))
        CALL calc_Psorption_parameter(dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
          & clay_sl(ic,is), silt_sl(ic,is), rmm_sorption(ic,is) * 1000._lp, ph_sl(ic,is), &
          & bulk_soil_carbon_sl(ic,is), soil_litter_carbon_sl(ic,is), sb_pool_mt(ix_woody_litter,ixC,ic,is), &
          & bulk_dens_sl(ic,is), volume_min_sl(ic,is), po4_solute(ic,is), &
          & qmax_po4_min_sl(ic,is), qmax_po4_om_sl(ic,is), qmax_po4(ic,is), km_adsorpt_po4_sl(ic,is), &
          & qmax_fast_po4(ic,is), qmax_slow_po4(ic,is), km_fast_po4(ic,is), km_slow_po4(ic,is), &
          & Qmax_AlFe_cor(ic,is), mineral_mass_density_sl=rho_particle_mineral * (1._lp - volume_pore_total_sl(ic,is)), &
          & flag_available_pore_volume=.TRUE.)
      END DO
    END DO
    max_residual = 0._lp
    DO ic = 1, nc
      DO is = 1, MIN(nsoil_w, NINT(num_sl_above_bedrock(ic)))
        density = dsl4jsb_Config(SB_)%sb_bulk_density_target_fraction * bulk_dens_sl(ic,is)
        max_residual = MAX(max_residual, ABS(bulk_dens_corr_sl(ic,is) - density) / MAX(1._lp, density))
      END DO
    END DO
    ! Keep output/convergence diagnostics consistent with the artificial mineral-inventory change.
    ! This does not update lagged nutrient status or advance the moving averages.
    CALL refresh_seed_pool_diagnostics(tile,options)
    IF (l_report .AND. my_process_is_stdio() .AND. iblk == 1 .AND. ics == 1) THEN
      WRITE(calibration_message,*) 'Pore calibration tile/block/year=', TRIM(tile%name), iblk, represented_year, &
        & '; local max pore change, density relative residual, bounded layers=', max_change, max_residual, bounded_count
      CALL message(modname//':calibrate_sb_pore_volume', TRIM(calibration_message))
      WRITE(calibration_message,*) 'Calibration imposed inventory increments (sum of local columns, mol m-2): N/P=', &
        & inventory_delta_n,inventory_delta_p,'; mineral-floor holds=',floor_hold_count
      CALL message(modname//':calibrate_sb_pore_volume', TRIM(calibration_message))
    END IF
  END SUBROUTINE calibrate_sb_pore_volume


  !> Replay the bookkept climatology until the maximum pseudo-time or convergence criterion is reached.

  ! ===================================================================================================== !
  !> Replay forcing through the JSM process sequence without vegetation or atmospheric model updates.
  SUBROUTINE run_jsm_forward_relaxation(tile, options, forcing, balance)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    TYPE(t_sb_spinup_forcing_climatology), INTENT(in) :: forcing
    TYPE(t_sb_replay_nutrient_balance), INTENT(inout) :: balance

    CALL run_sb_climatological_year(tile, options, forcing, 'jsm', balance)

  END SUBROUTINE run_jsm_forward_relaxation

  ! ===================================================================================================== !
  !> Replay forcing through the simple_1d process sequence without vegetation or atmospheric model updates.
  SUBROUTINE run_simple_1d_forward_relaxation(tile, options, forcing, balance)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    TYPE(t_sb_spinup_forcing_climatology), INTENT(in) :: forcing
    TYPE(t_sb_replay_nutrient_balance), INTENT(inout) :: balance

    CALL run_sb_climatological_year(tile, options, forcing, 'simple_1d', balance)

  END SUBROUTINE run_simple_1d_forward_relaxation

  ! ===================================================================================================== !
  !> Advance one climatological replay year using the configured host-calendar bins.
  SUBROUTINE run_sb_climatological_year(tile, options, forcing, sb_model_scheme, balance)

    USE mo_jsb_math_constants,   ONLY: one_day
    USE mo_jsb_process_class,    ONLY: SB_, HYDRO_
    USE mo_sb_plant_uptake_adapter, ONLY: t_sb_plant_uptake, prepare_combined_plant_uptake
    USE mo_sb_util,               ONLY: reset_sb_fluxes, calculate_time_average_soilbiogeochemistry
    USE mo_q_sb_rate_modifier,    ONLY: update_sb_rate_modifier
    USE mo_q_sb_jsm_main,         ONLY: update_sb_jsm
    USE mo_q_sb_ssm_main,         ONLY: update_sb_simple_model
    USE mo_q_sb_update_pools,     ONLY: update_sb_pools
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(HYDRO_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    TYPE(t_sb_spinup_forcing_climatology), INTENT(in) :: forcing
    TYPE(t_sb_replay_nutrient_balance), INTENT(inout) :: balance
    CHARACTER(len=*),           INTENT(in)    :: sb_model_scheme

    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_sl

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_task_options) :: replay_options
    CHARACTER(len=*), PARAMETER :: routine = modname//':run_sb_climatological_year'
    INTEGER :: forcing_bin, daily_bin
    INTEGER :: i_step
    INTEGER :: n_steps_per_bin
    REAL(lp) :: bin_duration
    REAL(lp) :: requested_dtime
    TYPE(t_sb_plant_uptake), ALLOCATABLE :: combined_uptake !< Absent on the ordinary PFT path
    REAL(lp), ALLOCATABLE :: combined_litter(:,:,:,:) !< Shared litter increments [mol m-3 timestep-1]
    REAL(lp), ALLOCATABLE :: effective_demand(:,:,:) !< N/P,cell,PFT after annual budget scaling [-]
    CHARACTER(len=16) :: uptake_scheme !< Existing kernel choice, not a new control [-]
    INTEGER :: ics, ice, iblk, it !< Chunk bounds/block and PFT index [-]
    LOGICAL :: uptake_valid !< Adapter input validation [-]
    INTEGER :: uptake_bin !< Last prepared seasonal bin in this internal year [-]

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(SB_)
    IF (ALLOCATED(forcing%uptake_drivers)) THEN
      IF (SIZE(forcing%uptake_drivers) /= forcing%n_forcing_bins) &
        & CALL finish(routine, 'Combined uptake forcing does not cover every seasonal bin.')
      ics = options%ics
      ice = options%ice
      iblk = options%iblk
      dsl4jsb_Get_memory(HYDRO_)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pot_sl)
      uptake_scheme = sb_model_scheme
      IF (sb_model_scheme == 'jsm') uptake_scheme = dsl4jsb_Config(SB_)%sb_adsorp_scheme
      ALLOCATE(combined_uptake)
    END IF

    requested_dtime = dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_timestep_days * one_day
    bin_duration = one_day
    n_steps_per_bin = MAX(1, CEILING(bin_duration / requested_dtime))

    replay_options = options
    replay_options%dtime = bin_duration / REAL(n_steps_per_bin, lp)
    IF (ALLOCATED(forcing%partitioned_litter_rate)) THEN
      IF (forcing%litter_partition_timestep /= replay_options%dtime) &
        & CALL finish(routine, 'Combined litter timestep differs from replay.')
      IF (SIZE(forcing%partitioned_litter_rate,5) /= forcing%n_forcing_bins) &
        & CALL finish(routine, 'Combined litter does not cover every seasonal bin.')
    END IF

    ! Replay a climatological year without advancing the host datetime. Use
    ! actual daily steps (or smaller requested steps), not one step per month.
    balance%litter_n = 0.0_lp
    balance%litter_p = 0.0_lp
    balance%uptake_n = 0.0_lp
    balance%uptake_p = 0.0_lp
    uptake_bin = 0
    DO daily_bin = 1, forcing%n_daily_bins
      forcing_bin = (daily_bin - 1) * forcing%n_forcing_bins / forcing%n_daily_bins + 1
      IF (ALLOCATED(forcing%partitioned_litter_rate)) &
        & combined_litter = forcing%partitioned_litter_rate(:,:,:,:,forcing_bin)*replay_options%dtime
      DO i_step = 1, n_steps_per_bin
        CALL load_sb_spinup_forcing_bin(tile, replay_options, forcing, forcing_bin, daily_bin, balance)
        CALL reset_sb_fluxes(tile, replay_options)
        CALL update_sb_rate_modifier(tile, replay_options)
        ! Root drivers, water potential and demand multipliers are fixed within
        ! a seasonal bin/year. Only transport coefficients vary daily. Prepare
        ! uptake once per bin, not at every replay substep; no solute is cached.
        IF (ALLOCATED(combined_uptake) .AND. forcing_bin /= uptake_bin) THEN
          ASSOCIATE(drivers => forcing%uptake_drivers(forcing_bin))
            IF (.NOT. ALLOCATED(drivers%fraction) .OR. .NOT. ALLOCATED(drivers%roots) .OR. &
              & .NOT. ALLOCATED(drivers%temperature) .OR. .NOT. ALLOCATED(drivers%phi_leaf_min) .OR. &
              & .NOT. ALLOCATED(drivers%vmax) .OR. .NOT. ALLOCATED(drivers%demand)) &
              & CALL finish(routine, 'Incomplete combined plant uptake drivers.')
            IF (ANY(SHAPE(drivers%demand) /= [2,SIZE(balance%scale_n),SIZE(drivers%fraction,2)])) &
              & CALL finish(routine, 'Invalid combined demand dimensions.')
            effective_demand = drivers%demand
            DO it = 1, SIZE(effective_demand,3)
              effective_demand(1,:,it) = MIN(1._lp, MAX(0._lp, drivers%demand(1,:,it)*balance%scale_n))
              effective_demand(2,:,it) = MIN(1._lp, MAX(0._lp, drivers%demand(2,:,it)*balance%scale_p))
            END DO
            CALL prepare_combined_plant_uptake(uptake_scheme,drivers%fraction,drivers%roots, &
              & drivers%temperature,wtr_soil_pot_sl,drivers%phi_leaf_min,drivers%vmax,effective_demand, &
              & combined_uptake,uptake_valid)
          END ASSOCIATE
          IF (.NOT. uptake_valid) CALL finish(routine, 'Invalid combined plant uptake drivers.')
          uptake_bin = forcing_bin
        END IF

        SELECT CASE (TRIM(sb_model_scheme))
        CASE ('jsm')
          CALL update_sb_jsm(tile, replay_options, l_soil_model_forcing=.TRUE., &
            & combined_uptake=combined_uptake, combined_litter=combined_litter)
        CASE ('simple_1d')
          CALL update_sb_simple_model(tile, replay_options, l_balance_weathering=.FALSE., &
            & combined_uptake=combined_uptake, combined_litter=combined_litter)
        CASE DEFAULT
          CALL finish(TRIM(routine), 'Unsupported soil-biogeochemistry model scheme.')
        END SELECT

        ! An unallocated actual denotes an absent OPTIONAL argument (Fortran
        ! 2008). Existing per-PFT calls therefore retain their exact science.
        CALL update_sb_pools(tile, replay_options, l_soil_model_forcing=.TRUE., combined_uptake=combined_uptake)
        CALL accumulate_replay_nutrient_budget(tile, replay_options, balance%litter_n, balance%litter_p, &
          & balance%uptake_n, balance%uptake_p)
        CALL calculate_time_average_soilbiogeochemistry(tile, replay_options)
      END DO
    END DO

  END SUBROUTINE run_sb_climatological_year

  ! ===================================================================================================== !
  !> Load seasonal biology/soil physics and the selected daily water-transport rates.
  SUBROUTINE load_sb_spinup_forcing_bin(tile, options, forcing, forcing_bin, daily_bin, balance)

    USE mo_jsb_math_constants, ONLY: eps8
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
    TYPE(t_sb_spinup_forcing_climatology), INTENT(in) :: forcing
    TYPE(t_sb_replay_nutrient_balance), INTENT(inout) :: balance
    INTEGER,                    INTENT(in)    :: forcing_bin, daily_bin

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
    CHARACTER(len=*), PARAMETER :: routine = modname//':load_sb_spinup_forcing_bin'
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    INTEGER :: i_element
    INTEGER :: i_soil
    INTEGER :: i_source
    INTEGER :: litter_entry
    INTEGER :: profile_entry, daily_profile_entry
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
    dsl4jsb_Real3D_onChunk :: growth_p_limit_smoothed_sl
    dsl4jsb_Real2D_onChunk :: nhx_deposition
    dsl4jsb_Real2D_onChunk :: noy_deposition
    dsl4jsb_Real2D_onChunk :: nhx_n15_deposition
    dsl4jsb_Real2D_onChunk :: noy_n15_deposition
    dsl4jsb_Real2D_onChunk :: p_deposition

    IF (forcing_bin < 1 .OR. forcing_bin > forcing%n_forcing_bins .OR. &
      & daily_bin < 1 .OR. daily_bin > forcing%n_daily_bins) THEN
      CALL finish(TRIM(routine), 'Requested forcing bin is outside the prepared climatology.')
    END IF

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
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, frac_wtr_transp_down_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, frac_w_lat_loss_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pot_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)
    dsl4jsb_Get_var2D_onChunk(VEG_, f_n_demand)
    dsl4jsb_Get_var2D_onChunk(VEG_, f_p_demand)
    dsl4jsb_Get_var3D_onChunk(SB_, growth_p_limit_smoothed_sl)
    dsl4jsb_Get_var2D_onChunk(A2L_, nhx_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, noy_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, nhx_n15_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, noy_n15_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, p_deposition)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_LITTERFALL_ID, seed_bed_litterfall_mt)

    veg_pool_mt(ix_fine_root, ixC, :) = forcing%fine_root_c(:, forcing_bin)
    f_n_demand = MIN(1.0_lp, MAX(0.0_lp, forcing%f_n_demand(:, forcing_bin) * balance%scale_n))
    f_p_demand = MIN(1.0_lp, MAX(0.0_lp, forcing%f_p_demand(:, forcing_bin) * balance%scale_p))
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
      growth_p_limit_smoothed_sl(:, i_soil) = forcing%growth_p_limit(:, profile_entry)
      root_fraction_sl(:, i_soil) = forcing%root_fraction(:, profile_entry)
      daily_profile_entry = daily_bin + forcing%n_daily_bins * (i_soil - 1)
      ! Do not silently clip an unstable daily removal: that would change export.
      ! Transport consumes true rate * dtime; restore HYDRO's native factor 1000.
      IF (ANY((forcing%wtr_down(:, daily_profile_entry) + forcing%wtr_lat_loss(:, daily_profile_entry)) &
        & * options%dtime > 1.0_lp + eps8)) THEN
        WRITE(message_text,*) 'Replay transport removal exceeds one; day/layer/dtime = ', &
          & daily_bin, i_soil, options%dtime, '; reduce sb_model_spinup_accelerator_timestep_days.'
        CALL finish(TRIM(routine), TRIM(message_text))
      END IF
      frac_wtr_transp_down_sl(:, i_soil) = forcing%wtr_down(:, daily_profile_entry) * 1000.0_lp
      frac_w_lat_loss_sl(:, i_soil) = forcing%wtr_lat_loss(:, daily_profile_entry) * 1000.0_lp
      t_soil_sl(:, i_soil) = forcing%t_soil(:, profile_entry)
      wtr_soil_sl(:, i_soil) = forcing%wtr_soil(:, profile_entry)
      ice_soil_sl(:, i_soil) = forcing%ice_soil(:, profile_entry)
      wtr_soil_pot_sl(:, i_soil) = forcing%wtr_pot(:, profile_entry)
    END DO

  END SUBROUTINE load_sb_spinup_forcing_bin

  ! ===================================================================================================== !
  SUBROUTINE accumulate_replay_nutrient_budget(tile, options, annual_litter_n, annual_litter_p, &
      & annual_plant_uptake_n, annual_plant_uptake_p)

    USE mo_jsb_process_class,    ONLY: SB_, HYDRO_
    USE mo_lnd_bgcm_idx,         ONLY: ixN, ixP
    USE mo_lnd_bgcm_store,       ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_store_class, ONLY: VEG_BGCM_LITTERFALL_ID, VEG_BGCM_SEED_BED_LITTERFALL_ID
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(HYDRO_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    REAL(lp),                   INTENT(inout) :: annual_litter_n(:)
    REAL(lp),                   INTENT(inout) :: annual_litter_p(:)
    REAL(lp),                   INTENT(inout) :: annual_plant_uptake_n(:)
    REAL(lp),                   INTENT(inout) :: annual_plant_uptake_p(:)

    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt1L2D :: seed_bed_litterfall_mt

    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl

    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, plant_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, plant_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, plant_uptake_po4_sl)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_LITTERFALL_ID, seed_bed_litterfall_mt)

    ! Litter matrices already contain mol m-2 for this replay timestep; uptake diagnostics are rates.
    annual_litter_n = annual_litter_n + SUM(veg_litterfall_mt(:,ixN,:), DIM=1) + &
      & seed_bed_litterfall_mt(ixN,:)
    annual_litter_p = annual_litter_p + SUM(veg_litterfall_mt(:,ixP,:), DIM=1) + &
      & seed_bed_litterfall_mt(ixP,:)
    annual_plant_uptake_n = annual_plant_uptake_n + &
      & SUM((plant_uptake_nh4_sl + plant_uptake_no3_sl) * soil_depth_sl, DIM=2) * options%dtime / 1.e6_lp
    annual_plant_uptake_p = annual_plant_uptake_p + &
      & SUM(plant_uptake_po4_sl * soil_depth_sl, DIM=2) * options%dtime / 1.e6_lp

  END SUBROUTINE accumulate_replay_nutrient_budget

  !> Diagnose a dimensionless coupled litter/uptake reference from integrated mol m-2 budgets.
  !> Zero denotes an inactive reference, not a replacement ratio. Retain its recorded demand.
  SUBROUTINE diagnose_replay_reference_ratio(litter, uptake, interactive, ratio, element_name)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    USE mo_jsb_math_constants, ONLY: eps8
    REAL(lp), INTENT(in) :: litter(:), uptake(:)
    LOGICAL, INTENT(in) :: interactive
    REAL(lp), INTENT(out) :: ratio(:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: element_name !< nutrient label for rejection diagnostics
    INTEGER :: ic
    CHARACTER(len=512) :: text
    CHARACTER(len=*), PARAMETER :: routine = modname//':diagnose_replay_reference_ratio'

    ratio = 0.0_lp
    IF (.NOT. interactive) RETURN
    DO ic = 1, SIZE(litter)
      IF (.NOT. ieee_is_finite(litter(ic)) .OR. .NOT. ieee_is_finite(uptake(ic))) &
        & CALL finish(routine, 'Non-finite coupled litter/uptake reference.')
      IF (litter(ic) < 0.0_lp .OR. uptake(ic) < 0.0_lp) &
        & CALL finish(routine, 'Negative coupled litter/uptake reference.')
      IF (litter(ic) <= eps8 .AND. uptake(ic) <= eps8) CYCLE
      IF (litter(ic) <= eps8 .OR. uptake(ic) <= eps8) THEN
        WRITE(text,*) 'One-sided zero coupled reference at chunk cell ', ic, '; litter, uptake [mol m-2]: ', &
          & litter(ic), uptake(ic)
        IF (PRESENT(element_name)) text = TRIM(text)//'; element='//element_name
        CALL finish(routine, TRIM(text))
      END IF
      ratio(ic) = litter(ic) / uptake(ic)
    END DO
  END SUBROUTINE diagnose_replay_reference_ratio

  !> Adjust next replay year's plant demand to the coupled litter/uptake ratio, independently for N/P.
  !> This changes a demand multiplier, not pools or nutrient supply. Actual uptake remains process-limited.
  SUBROUTINE update_replay_plant_np_demand_balance(tile, options, forcing, balance, &
      & execution_index, represented_year, l_balanced)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    USE mo_jsb_math_constants, ONLY: eps8
    USE mo_jsb_process_class, ONLY: SB_, HYDRO_
    USE mo_jsb_parallel, ONLY: my_process_is_stdio
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(HYDRO_)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options
    TYPE(t_sb_spinup_forcing_climatology), INTENT(in) :: forcing
    TYPE(t_sb_replay_nutrient_balance), INTENT(inout) :: balance
    INTEGER, INTENT(in) :: execution_index, represented_year
    LOGICAL, INTENT(out) :: l_balanced
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: nh4_solute, nh4_assoc, no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute, po4_assoc_fast, po4_assoc_slow
    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: ics, ice, iblk, ic, element
    LOGICAL :: interactive(2)
    REAL(lp) :: litter(2), uptake(2), reference(2), scale(2), stock(2) ! [mol m-2], ratios [-]
    REAL(lp) :: target, error, max_error, tolerance ! target [mol m-2 yr-1], errors [-]
    CHARACTER(len=1024) :: text
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_replay_plant_np_demand_balance'

    ics = options%ics
    ice = options%ice
    iblk = options%iblk
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SB_, nh4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_, nh4_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_, no3_solute)
    dsl4jsb_Get_var3D_onChunk(SB_, po4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_, po4_assoc_fast)
    dsl4jsb_Get_var3D_onChunk(SB_, po4_assoc_slow)
    interactive = [dsl4jsb_Config(SB_)%flag_sb_interactive_n, dsl4jsb_Config(SB_)%flag_sb_interactive_p]
    tolerance = dsl4jsb_Config(SB_)%sb_model_spinup_accelerator_convergence_rel_tol
    max_error = 0.0_lp
    DO ic = 1, SIZE(balance%scale_n)
      litter = [balance%litter_n(ic), balance%litter_p(ic)]
      uptake = [balance%uptake_n(ic), balance%uptake_p(ic)]
      reference = [forcing%litter_to_uptake_n(ic), forcing%litter_to_uptake_p(ic)]
      scale = [balance%scale_n(ic), balance%scale_p(ic)]
      stock(1) = SUM((nh4_solute(ic,:) + nh4_assoc(ic,:) + no3_solute(ic,:)) * soil_depth_sl(ic,:))
      stock(2) = SUM((po4_solute(ic,:) + po4_assoc_fast(ic,:) + po4_assoc_slow(ic,:)) * soil_depth_sl(ic,:))
      DO element = 1, 2
        IF (.NOT. interactive(element) .OR. reference(element) == 0.0_lp) CYCLE
        IF (.NOT. ieee_is_finite(litter(element)) .OR. .NOT. ieee_is_finite(uptake(element))) &
          & CALL finish(routine, 'Non-finite annual replay nutrient budget.')
        IF (litter(element) < 0.0_lp .OR. uptake(element) < 0.0_lp) &
          & CALL finish(routine, 'Negative annual replay nutrient budget.')
        target = litter(element) / reference(element)
        error = ABS(uptake(element) - target) / MAX(ABS(uptake(element)), ABS(target), eps8)
        max_error = MAX(max_error, error)
        IF (stock(element) <= eps8 .AND. uptake(element) < target * (1.0_lp - tolerance)) THEN
          WRITE(text,*) 'Active nutrient depleted: tile/block/cell/event/year/element(N=1,P=2) ', &
            & TRIM(tile%name), iblk, ic, execution_index, represented_year, element, &
            & '; stock, uptake, target [mol m-2] = ', stock(element), uptake(element), target
          CALL finish(routine, TRIM(text))
        END IF
        IF (uptake(element) > eps8) THEN
          ! The multiplier may exceed one; only effective demand is capped to [0,1] when loaded.
          scale(element) = MAX(0.0_lp, scale(element) * target / uptake(element))
        ELSE
          scale(element) = 1.0_lp
        END IF
      END DO
      IF (ANY(.NOT. ieee_is_finite(scale))) CALL finish(routine, 'Non-finite replay demand multiplier.')
      balance%scale_n(ic) = scale(1)
      balance%scale_p(ic) = scale(2)
    END DO
    l_balanced = max_error <= tolerance
    IF (model%config%flag_spinup_debug .AND. my_process_is_stdio() .AND. iblk == 1 .AND. ics == 1 &
      & .AND. is_sb_spinup_debug_year(represented_year)) THEN
      WRITE(text,*) 'Replay demand balance tile/event/year = ', TRIM(tile%name), execution_index, represented_year, &
        & '; local max relative error, N/P scale ranges = ', max_error, &
        & MINVAL(balance%scale_n), MAXVAL(balance%scale_n), MINVAL(balance%scale_p), MAXVAL(balance%scale_p)
      CALL message(routine, TRIM(text))
    END IF
  END SUBROUTINE update_replay_plant_np_demand_balance


  !> Return true for the sparse pseudo-years included in the detailed accelerator diagnostics.
  LOGICAL FUNCTION is_sb_spinup_debug_year(represented_year)

    INTEGER, INTENT(in) :: represented_year

    is_sb_spinup_debug_year = ANY(sb_spinup_debug_years == represented_year)

  END FUNCTION is_sb_spinup_debug_year

  ! ===================================================================================================== !
  !> Bounded pre-seed diagnostics to distinguish restored state from derived physics.
  !! First owned sample only, once per active tile/event, never per replay step.
  SUBROUTINE write_sb_replay_entry_debug(tile, options, event)
    USE mo_jsb_process_class, ONLY: SB_, HYDRO_
    USE mo_jsb_parallel, ONLY: my_process_is_stdio
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(HYDRO_)
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options
    INTEGER, INTENT(in) :: event !< Scheduled replay index [-]
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(HYDRO_)
    INTEGER :: ics, ice, iblk, is
    CHARACTER(len=2048) :: text
    dsl4jsb_Real2D_onChunk :: num_sl_above_bedrock
    dsl4jsb_Real3D_onChunk :: nh4_solute
    dsl4jsb_Real3D_onChunk :: nh4_assoc
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast
    dsl4jsb_Real3D_onChunk :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: enzyme_frac_AP
    dsl4jsb_Real3D_onChunk :: fact_n_status_mic_c_growth
    dsl4jsb_Real3D_onChunk :: fact_p_status_mic_c_growth
    dsl4jsb_Real3D_onChunk :: microbial_cue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk :: microbial_nue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk :: microbial_pue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk :: vol_field_cap_sl
    dsl4jsb_Real3D_onChunk :: vol_porosity_sl

    IF (.NOT. my_process_is_stdio() .OR. options%iblk /= 1 .OR. options%ics /= 1) RETURN
    ics=options%ics; ice=options%ice; iblk=options%iblk
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, num_sl_above_bedrock)
    dsl4jsb_Get_var3D_onChunk(SB_, nh4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_, nh4_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_, no3_solute)
    dsl4jsb_Get_var3D_onChunk(SB_, po4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_, po4_assoc_fast)
    dsl4jsb_Get_var3D_onChunk(SB_, po4_assoc_slow)
    dsl4jsb_Get_var3D_onChunk(SB_, enzyme_frac_AP)
    dsl4jsb_Get_var3D_onChunk(SB_, fact_n_status_mic_c_growth)
    dsl4jsb_Get_var3D_onChunk(SB_, fact_p_status_mic_c_growth)
    dsl4jsb_Get_var3D_onChunk(SB_, microbial_cue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_, microbial_nue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_, microbial_pue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_field_cap_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl)
    DO is=1,MIN(SIZE(vol_porosity_sl,2),NINT(num_sl_above_bedrock(1)))
      WRITE(text,*) 'Replay entry tile/event/layer=',TRIM(tile%name),event,is, &
        & '; field capacity/porosity [-]=',vol_field_cap_sl(1,is),vol_porosity_sl(1,is), &
        & '; NH4sol/NH4assoc/NO3/PO4sol/PO4fast/PO4slow [mol m-3]=', &
        & nh4_solute(1,is),nh4_assoc(1,is),no3_solute(1,is), &
        & po4_solute(1,is),po4_assoc_fast(1,is),po4_assoc_slow(1,is)
      CALL message(modname,TRIM(text))
      WRITE(text,*) 'Replay entry biology tile/event/layer=',TRIM(tile%name),event,is, &
        & '; enzymeAP/Nstatus/Pstatus/CUEavg/NUEavg/PUEavg [-]=',enzyme_frac_AP(1,is), &
        & fact_n_status_mic_c_growth(1,is),fact_p_status_mic_c_growth(1,is), &
        & microbial_cue_eff_tmic_mavg(1,is),microbial_nue_eff_tmic_mavg(1,is),microbial_pue_eff_tmic_mavg(1,is)
      CALL message(modname,TRIM(text))
    END DO
  END SUBROUTINE write_sb_replay_entry_debug

  !> Write the event-local forcing climatology used by one forward SB acceleration.
  SUBROUTINE write_sb_spinup_forcing_debug(tile, options, execution_index, forcing)
    USE mo_jsb_parallel, ONLY: my_process_is_stdio
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options
    INTEGER, INTENT(in) :: execution_index
    TYPE(t_sb_spinup_forcing_climatology), INTENT(in) :: forcing
    CHARACTER(len=2048) :: text
    CHARACTER(len=*), PARAMETER :: routine = modname//':write_sb_spinup_forcing_debug'
    INTEGER :: bin !< Monthly climatology index; bounds log volume to 12 records per event/tile [-]

    ! Bounded local diagnostics, not a global gather or a per-cell file writer.
    IF (.NOT. my_process_is_stdio() .OR. options%iblk /= 1 .OR. options%ics /= 1) RETURN
    WRITE(text,*) 'Replay forcing tile/event=', TRIM(tile%name), execution_index, &
      & '; monthly/daily bins=', forcing%n_forcing_bins, forcing%n_daily_bins, &
      & '; local sampled seconds min/max=', MINVAL(SUM(forcing%forcing_time, DIM=2)), &
      & MAXVAL(SUM(forcing%forcing_time, DIM=2))
    CALL message(routine, TRIM(text))
    WRITE(text,*) 'Local forcing ranges: litter [mol m-2 s-1]=', &
      & MINVAL(forcing%litterfall_rate), MAXVAL(forcing%litterfall_rate), &
      & '; N/P demand [-]=', MINVAL(forcing%f_n_demand), MAXVAL(forcing%f_n_demand), &
      & MINVAL(forcing%f_p_demand), MAXVAL(forcing%f_p_demand), &
      & '; transformed plant P limitation [-]=', MINVAL(forcing%growth_p_limit), MAXVAL(forcing%growth_p_limit)
    CALL message(routine, TRIM(text))
    WRITE(text,*) 'Local daily transport rate ranges [s-1]: down=', &
      & MINVAL(forcing%wtr_down), MAXVAL(forcing%wtr_down), &
      & '; lateral=', MINVAL(forcing%wtr_lat_loss), MAXVAL(forcing%wtr_lat_loss), &
      & '; monthly soil temperature [K]=', MINVAL(forcing%t_soil), MAXVAL(forcing%t_soil)
    CALL message(routine, TRIM(text))
    ! Min/max alone cannot detect reordered bins or loss of partial-window sums.
    ! Report the actual demand/deposition means and weights used by replay.
    DO bin=1,forcing%n_forcing_bins
      WRITE(text,*) 'Replay monthly forcing tile/event/bin=',TRIM(tile%name),execution_index,bin, &
        & '; seconds/rootC/Ndemand/Pdemand/NHx/NOy/Pdep=',forcing%forcing_time(1,bin), &
        & forcing%fine_root_c(1,bin),forcing%f_n_demand(1,bin),forcing%f_p_demand(1,bin), &
        & forcing%nhx_deposition(1,bin),forcing%noy_deposition(1,bin),forcing%p_deposition(1,bin), &
        & '; units: s, mol C m-2, 1, 1, micro-mol m-2 s-1 (depositions)'
      CALL message(routine,TRIM(text))
    END DO
  END SUBROUTINE write_sb_spinup_forcing_debug

  ! ===================================================================================================== !
  !> Write a sparse pool and physical-state snapshot during forward pseudo-time.
  SUBROUTINE write_sb_spinup_snapshot_debug(tile, options, execution_index, represented_year, &
      & max_relative_change, max_total_c_absolute_change, max_total_c_relative_change, l_converged)

    USE mo_jsb_parallel,      ONLY: my_process_is_stdio
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
    CHARACTER(len=2048) :: text
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

    IF (.NOT. my_process_is_stdio() .OR. options%iblk /= 1 .OR. options%ics /= 1) RETURN

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

    WRITE(text,*) 'Replay sample profile: tile/event/pseudo-year=', TRIM(tile%name), execution_index, represented_year
    CALL message(routine, TRIM(text))
    IF (PRESENT(max_relative_change) .AND. PRESENT(max_total_c_absolute_change) .AND. &
        & PRESENT(max_total_c_relative_change) .AND. PRESENT(l_converged)) THEN
      WRITE(text,*) 'RELAXATION_STATUS max_pool_relative_change, total_C_absolute_change_mol_m-2, ', &
        & 'total_C_relative_change, converged = ', max_relative_change, max_total_c_absolute_change, &
        & max_total_c_relative_change, l_converged
      CALL message(routine, TRIM(text))
    END IF
    CALL message(routine, 'First local cell only: pool C [mol m-3], layer stock [mol m-2], '// &
      & 'bulk/litter C [mol m-3], density [kg m-3], mineral fraction, qmax_org, burial fraction.')
    DO ic = 1, MIN(1, SIZE(pool_snapshot, 3))
      n_active_layers = MIN(SIZE(pool_snapshot, 4), NINT(num_sl_above_bedrock(ic)))
      DO is = 1, n_active_layers
        WRITE(text,*) 'STATE ', execution_index, represented_year, ic, is, &
          & pool_snapshot(:,ixC,ic,is), SUM(pool_snapshot(:,ixC,ic,is)), &
          & soil_lay_width_sl(ic,is), SUM(pool_snapshot(:,ixC,ic,is)) * soil_lay_width_sl(ic,is), &
          & bulk_soil_carbon_sl(ic,is), soil_litter_carbon_sl(ic,is), bulk_carbon_max(ic,is), &
          & bulk_dens_corr_sl(ic,is), volume_min_sl(ic,is), qmax_org(ic,is), surface_burial_frac(ic,is)
        CALL message(routine, TRIM(text))
      END DO
    END DO

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

    USE mo_jsb_parallel, ONLY: my_process_is_stdio
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

    IF (.NOT. my_process_is_stdio() .OR. options%iblk /= 1 .OR. options%ics /= 1) RETURN
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
    DO ic = 1, MIN(1, SIZE(pool_before, 3))
      n_active_layers = MIN(SIZE(pool_before, 4), NINT(num_sl_above_bedrock(ic)))
      DO is = 1, n_active_layers
        total_c_before = total_c_before + &
          & SUM(pool_before(:,ixC,ic,is)) * soil_depth_sl(ic,is)
        total_c_after = total_c_after + &
          & SUM(pool_after(:,ixC,ic,is)) * soil_depth_sl(ic,is)
      END DO
    END DO

    WRITE(debug_message,*) 'DEBUG first local cell SB carbon totals execution/represented_years = ', &
      & execution_index, represented_years, '; before/after/delta [mol C m-2] = ', &
      & total_c_before, total_c_after, total_c_after - total_c_before
    CALL message(TRIM(routine), TRIM(debug_message))

    DO ic = 1, MIN(1, SIZE(pool_before, 3))
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

  ! ===================================================================================================== !
  !> Commit the relaxed SB state, restore live forcing, and refresh derived soil properties.
  SUBROUTINE finalize_sb_spinup_accelerator(tile, options, live_forcing)

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
    TYPE(t_sb_spinup_live_forcing), INTENT(in) :: live_forcing

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
    dsl4jsb_Real3D_onChunk :: growth_p_limit_smoothed_sl
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
    dsl4jsb_Get_var3D_onChunk(SB_, growth_p_limit_smoothed_sl)
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
    growth_p_limit_smoothed_sl = live_forcing%growth_p_limit
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
    TYPE(t_soil_forcing_cache), POINTER :: next => NULL() !< next independent model/grid cache
    INTEGER :: owner_model_id = -1 !< owner of grid-distributed input streams and PFT layout [-]
    INTEGER :: owner_grid_id = -1 !< grid used by the host NetCDF reader [-]
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

  ! One bounded daily cache per model/grid, shared across its PFTs. A linked
  ! registry avoids assuming contiguous model IDs or a fixed number of domains.
  ! Both selecting the active pointer and consuming it occur under the I/O lock.
  TYPE(t_soil_forcing_cache), POINTER, SAVE :: soil_forcing_cache => NULL()
  TYPE(t_soil_forcing_cache), POINTER, SAVE :: soil_forcing_cache_head => NULL()
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
    USE mo_jsb_time,                ONLY: get_date_components, get_year_day

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
    ! The cache is process-local shared scratch. Hold the lock through copying,
    ! not just NetCDF loading: another chunk/model must not replace it meanwhile.
    !$OMP CRITICAL (soil_forcing_netcdf_io)
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
      IF (pft_idx > soil_forcing_cache%npft) CALL finish(routine, 'PFT index exceeds forcing-cache layout.')
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

    !$OMP END CRITICAL (soil_forcing_netcdf_io)
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

    ! File records are equally spaced integer seconds within a daily file.
    ! Validate before division; a fractional/truncated index would load wrong data.
    IF (.NOT. (dtime >= 1.0_lp .AND. dtime <= 86400.0_lp)) &
      & CALL finish(modname//':get_soil_forcing_preload_nsteps', 'QSOIL timestep must be within [1,86400] seconds.')
    IF (REAL(NINT(dtime), lp) /= dtime) &
      & CALL finish(modname//':get_soil_forcing_preload_nsteps', 'QSOIL file timestep must use integer seconds.')
    IF (MOD(86400, NINT(dtime)) /= 0) &
      & CALL finish(modname//':get_soil_forcing_preload_nsteps', 'QSOIL file timestep must divide a day exactly.')
    get_soil_forcing_preload_nsteps = 86400 / NINT(dtime)

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
    ! Caller holds soil_forcing_netcdf_io until all destination copies finish.
    ! Same-shaped grids and identical filenames do not imply the same distributed
    ! ownership. Keep isotope-initialization state independent for each owner too.
    CALL select_soil_forcing_cache(model%id, model%grid_id)
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

      CALL load_soilphysics_cache(preload_nsteps, cache_start_timestep, skip_isotope_reload, caller)
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
    IF (PRESENT(cache_hit_out)) cache_hit_out = cache_hit

  END SUBROUTINE ensure_soil_forcing_cache_loaded

  !> Select/create process-local cache ownership; caller must hold the I/O lock.
  SUBROUTINE select_soil_forcing_cache(model_id, grid_id)
    INTEGER, INTENT(in) :: model_id, grid_id
    soil_forcing_cache => soil_forcing_cache_head
    DO WHILE (ASSOCIATED(soil_forcing_cache))
      IF (soil_forcing_cache%owner_model_id == model_id .AND. &
        & soil_forcing_cache%owner_grid_id == grid_id) RETURN
      soil_forcing_cache => soil_forcing_cache%next
    END DO
    ALLOCATE(soil_forcing_cache)
    soil_forcing_cache%owner_model_id = model_id
    soil_forcing_cache%owner_grid_id = grid_id
    soil_forcing_cache%next => soil_forcing_cache_head
    soil_forcing_cache_head => soil_forcing_cache
  END SUBROUTINE select_soil_forcing_cache

  ! ======================================================================================================= !
  !> Ensure that the cache arrays match the current model dimensions.
  !>
  SUBROUTINE ensure_soil_forcing_cache_shape(bgcm_store, root_fraction_sl, npft, preload_nsteps, caller)

    USE mo_jsb_parallel, ONLY: my_process_is_stdio
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
    REAL(lp) :: cache_bytes !< local payload plus NetCDF scratch, excluding framework/file buffers [byte]
    CHARACTER(len=512) :: cache_message

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
    ! Physics has no PFT dimension; only biological forcing is PFT-specific.
    ! Report before allocation so a failed large global request remains diagnosable.
    cache_bytes = REAL(STORAGE_SIZE(0.0_lp) / 8, lp) * REAL(nproma, lp) * REAL(nblks, lp) * &
      & (REAL(preload_nsteps, lp) * (REAL(nsoil, lp) * REAL(npft, lp) + 3.0_lp * REAL(npft, lp) + &
      & 5.0_lp * REAL(nsoil, lp) + 5.0_lp + REAL(npft, lp) * &
      & (REAL(soil_forcing_cache%veg_litter_nr_parts, lp) * REAL(soil_forcing_cache%veg_litter_nr_elems, lp) + &
      & REAL(soil_forcing_cache%seed_bed_nr_elems, lp))) + REAL(nsoil + 1, lp))
    IF (my_process_is_stdio()) THEN
      WRITE(cache_message,*) 'QSOIL cache model/grid, local cells, PFTs, records, estimated local bytes = ', &
        & soil_forcing_cache%owner_model_id, soil_forcing_cache%owner_grid_id, &
        & REAL(nproma, lp) * REAL(nblks, lp), npft, preload_nsteps, cache_bytes
      CALL message(TRIM(caller), TRIM(cache_message))
    END IF
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
    ! Reallocated isotope arrays contain no first-record forcing yet.
    soil_forcing_cache%isotope_forcing_initialized = .FALSE.
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
  SUBROUTINE load_soilphysics_cache(nsteps, start_timestep, skip_isotope_reload, caller)

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
    CHARACTER(len=5) :: tile_name
    REAL(lp) :: t_total_start, t_total_end

    tile_name = get_pft_tile_name(1)
    ! Retain store validation even though cached dimensions were checked earlier.
    ! get_bgcm_idx can report a missing material; these calls are not dead code.
    veg_litter_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_LITTERFALL_ID, tile_name, caller)
    seed_bed_idx   = get_bgcm_idx(bgcm_store, VEG_BGCM_SEED_BED_LITTERFALL_ID, tile_name, caller)
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
