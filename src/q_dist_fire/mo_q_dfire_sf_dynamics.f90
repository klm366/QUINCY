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
!> Calculate Spitfire fire dynamics
!>
!> -update_dfire_spitfire: integrate Quincy with the spitfire process:
!>    - convert units from Quincy into Thonicke's Spitfire
!>    - gather lctlib parameters
!>    - select precalculated daily precipitation method
!>    - Soil moisture content
!>    - Fraction of the vegetation
!>    - Call calc_spitfire (science)
!>    - Save data into process memory
!>
!> -calc_spitfire: spitfire subroutine which calls
!>    - Initialize t_sf_dynamics, used to output data
!>    - Calc Nesterov index (mo_q_dfire_sf_nesterov.f90)
!>    - Calc Thonicke's / Oberhagemann implementation
!>    -   Calc fuel moisture (mo_q_dfire_sf_fmoist_class.f90)
!>    -   Calc Fire danger index (-)
!>    -   Calc Rate of spread (mo_q_dfire_sf_ros_class.f90)
!>    - Calc Area Burned and its fraction (-)
!>    - Calc ignition events (mo_q_difre_sf_ignition.f90)
!>    - Calc dead fraction of consumption (mo_q_dfire_sf_frac_consump_class.f90)
!>    - Calc dead fuel consumption (mo_q_dfire_sf_fuel_consumption.f90)
!>    - Calc fire intensity (-)
!>
MODULE mo_q_dfire_sf_dynamics
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_lnd_bgcm_store,        ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,  ONLY: SB_BGCM_LOSS_ID, SB_BGCM_POOL_ID
  USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract

  USE mo_q_dfire_sf_constants,  ONLY: sl_ixs => df_sf_sl_indices !< Soil related constants: top soil moisture layers

  USE mo_q_dfire_sf_constants,  ONLY: df_sf_max_density_human_pop
  USE mo_q_dfire_sf_constants,  ONLY: df_sf_fc_spread_nr
  USE mo_q_dfire_sf_constants,  ONLY: df_sf_alpha_nr, df_sf_alphas !< Drying parameters

  USE mo_q_dfire_sf_constants,  ONLY: df_sf_calorific_heat_content, df_sf_isurface_extinc_th

  USE mo_q_dfire_sf_constants,  ONLY: sf_fuel_dead_fuel_nr, sf_fuel_live_fuel_nr, df_sf_fc_full_nr, &
                                      df_sf_fc_1hr_idx, df_sf_fc_10hr_idx, df_sf_fc_100hr_idx, df_sf_fc_1000hr_idx !< fuel classes indices

  USE mo_q_dfire_sf_constants,  ONLY: df_sf_omegas_nr, & !< Fuel moisture dimension size (for omega_o and omega_nl)
                                      df_sf_omega_o_idx, & !< omega_o index from calc_sf_fuel_moisture_for_thonicke output (array)
                                      df_sf_omega_nl_idx

  USE mo_q_dfire_sf_constants,  ONLY: ixs_sf_dead_litter_fuel_pools
  USE mo_q_dfire_sf_constants,  ONLY: df_sf_units_g_to_kg
  USE mo_q_dfire_sf_constants,  ONLY: df_sf_sigma_dd_nr, df_sf_sigma_lv_nr, df_sf_sigma_dd_ros_nr
  USE mo_q_dfire_sf_constants,  ONLY: df_sigmas_dead, df_sigmas_live, df_alphas_dead
  USE mo_q_dfire_sf_config_class,  ONLY: Q_DFIRE_ROS_THONICKE, Q_DFIRE_ROS_OBERHAGEMANN

  USE mo_q_dfire_sf_fuel_class,   ONLY: t_sf_fuel_class, sf_calc_has_any_fuel, sf_calc_aggr_fuels
  USE mo_q_dfire_sf_fmoist_class, ONLY: wrap_calc_sf_fuel_moisture_for_thonicke, wrap_calc_sf_fuel_moisture_for_oberhagemann
  USE mo_q_dfire_sf_ignition_class,     ONLY: wrap_calc_sf_ignitions
  USE mo_q_dfire_sf_nesterov_class,     ONLY: calc_sf_nesterov_index
  USE mo_q_dfire_sf_ros_class,          ONLY: t_sf_veget_cover_frac, sf_calc_rate_of_spread
  USE mo_q_dfire_sf_ros_thon_class,     ONLY: sf_calc_rate_of_spread_thon
  USE mo_q_dfire_sf_ros_ober_class,     ONLY: sf_calc_rate_of_spread_ober, calc_live_me
  USE mo_q_dfire_sf_classes,            ONLY: t_sf_ros, t_sf_dynamics, t_sf_ignitions, t_sf_mf, t_sf_atm, t_sf_fuel_load, t_sf_config
  USE mo_q_dfire_sf_frac_consump_class, ONLY: sf_calc_frac_consump
  USE mo_q_dfire_sf_fuel_consump_class, ONLY: sf_calc_dead_fuel_consumption

  USE mo_q_dfire_sf_adaptor_class,    ONLY: sf_convert_units
  USE mo_q_dfire_sf_util,       ONLY: logical2real
  USE mo_quincy_stats,          ONLY: stats_average_weight

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_q_dfire_sf_dynamics, sf_calc_fire_intensity, calc_ab_moisture

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_dynamics'

CONTAINS

  !>
  !> Compute spitfire for current time step
  !> It links Quincy and Spitfire
  !>
  SUBROUTINE update_q_dfire_sf_dynamics(tile, options)
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_process_class,       ONLY: VEG_, Q_DFIRE_SF_, HYDRO_, Q_DFIRE_FUEL_SF_, A2L_
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_grid,                ONLY: Get_vgrid, Get_grid
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid, t_jsb_grid
    USE mo_jsb_math_constants,      ONLY: one_day

    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(Q_DFIRE_SF_)
    dsl4jsb_Use_memory(Q_DFIRE_FUEL_SF_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(VEG_)
    !dsl4jsb_Use_memory(SB_)
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
    INTEGER                               :: iblk, ics, ice, nc, ic, is
    INTEGER                               :: no_children
    REAL(wp)                              :: dtime      !< timestep
    CHARACTER(len=*), PARAMETER           :: routine = TRIM(modname)//':update_q_dfire_sf_dynamics'

    REAL(wp), POINTER                     :: grid_cell_area(:) ! [m2]
    REAL(wp), POINTER                     :: tile_fract(:) ! [%]
    LOGICAL                               :: mask_is_tile_enabled(options%nc)


    TYPE(t_sf_mf)                         :: sf_fuel_moisture_all
    TYPE(t_sf_veget_cover_frac)           :: veg_cover_frac                               !< Vegetation cover
    TYPE(t_sf_atm)                        :: sf_atm                                       !< Atmosphere related variables / Weather
    TYPE(t_sf_ignitions)                  :: sf_ignitions                                 !< Ignition variables
    TYPE(t_sf_ros)                        :: sf_ros                                       !< Rate of spread
    TYPE(t_sf_fuel_load)                  :: sf_fuel_load                                 !< Fuel load
    TYPE(t_sf_config)                     :: sf_config                                    !< spitfire config

    REAL(wp)                              :: burnable_area_ha(options%nc)                 !< Burnable area [ha]

    TYPE(t_sf_dynamics)                   :: sf_dynamics

    INTEGER                               :: nr_of_children

    INTEGER                               :: ierror
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_DFIRE_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_DFIRE_SF_)
    dsl4jsb_Def_memory(Q_DFIRE_FUEL_SF_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Local vars here
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real3D_onChunk                :: wtr_soil_sl        ! soil water content [m]
    dsl4jsb_Real3D_onChunk                :: wtr_soil_sat_sl    ! Soil water content at saturation (calculated from
                                                                ! soil porosity, reduced by ice)                        [m]
    dsl4jsb_Real3D_onChunk                :: soil_depth_sl      ! Soil depth until bedrock [m]

    !< Ignitions
    dsl4jsb_Real2D_onChunk                :: e_n_ignitions    ! out, Total rate of ignitions events [ha-1 ts-1]
    dsl4jsb_Real2D_onChunk                :: e_n_humans       ! out, Rate of human-events ignitions [ha-1 ts-1]
    dsl4jsb_Real2D_onChunk                :: e_n_lightning    ! out, Rate of lightning-events ignitions [ha-1 ts-1]

    ! DFUEL_ 2D
    ! Note: t_air_max/t_air_min belongs to the previous days due to implementation constraints
    dsl4jsb_Real2D_onChunk                :: t_air_min         ! Min temperature of the day (previous day) [k]
    dsl4jsb_Real2D_onChunk                :: t_air_max         ! Max temperature of the day (previous day) [k]
    dsl4jsb_Real2D_onChunk                :: precip_daily_mmd  ! out [mm day-1]
    dsl4jsb_Real2D_onChunk                :: precip_daily      ! in [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk                :: wind_10m          ! in A2L_ [m s-1]

    dsl4jsb_Real2D_onChunk                :: pot_burn_area     ! Potential burnable area (forest and grass areas) [m2?]

    dsl4jsb_Real2D_onChunk                :: mask_is_fire_intense_enough ! Mask relative to gridcell to define whether the current
                                                               ! fire is intense enough to exist  TRUE=1, FALSE=0 [-]
    dsl4jsb_Real2D_onChunk                :: nesterov_index    ! inout [Celsius 2]
    dsl4jsb_Real2D_onChunk                :: ab_moisture_diag  ! out [-]
    dsl4jsb_Real3D_onChunk                :: fuel_moisture     ! out [-]
    dsl4jsb_Real2D_onChunk                :: fdi               ! out [%]
    dsl4jsb_Real2D_onChunk                :: ros_f             ! out [?]
    dsl4jsb_Real2D_onChunk                :: ros_b             ! out [?]
    dsl4jsb_Real2D_onChunk                :: mean_fire_area    ! out [?]
    dsl4jsb_Real2D_onChunk                :: fire_duration     ! out [?]
    dsl4jsb_Real2D_onChunk                :: distance_travelled ! out [?]
    dsl4jsb_Real2D_onChunk                :: l2b_ratio         ! out [?]

    dsl4jsb_Real2D_onChunk                :: reaction_intensity       ! out
    dsl4jsb_Real2D_onChunk                :: prop_flux_ratio          ! out
    dsl4jsb_Real2D_onChunk                :: wind_factor              ! out
    dsl4jsb_Real2D_onChunk                :: eff_heating_number       ! out
    dsl4jsb_Real2D_onChunk                :: preig_heat               ! out

    dsl4jsb_Real2D_onChunk                :: opt_react_velo           ! out
    dsl4jsb_Real2D_onChunk                :: net_fuel_load            ! out
    dsl4jsb_Real2D_onChunk                :: moist_damp_coeff         ! out
    dsl4jsb_Real2D_onChunk                :: mineral_damp_coeff       ! out

    ! Q_DFIRE_FUEL_ 3D
    dsl4jsb_Real3D_onChunk                :: dead_fuel_1hr_dfp          ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_10hr_dfp         ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_100hr_dfp        ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_1000hr_dfp       ! in [g DM m-2]

    dsl4jsb_Real3D_onChunk                :: live_grass_fuel_1hr_lfp    ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: cured_live_grass_fuel_1hr_dfp    ! in [g DM m-2]

    dsl4jsb_Real3D_onChunk                :: live_tree_fuel_1hr_lfp     ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_tree_fuel_10hr_lfp    ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_tree_fuel_100hr_lfp   ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_tree_fuel_1000hr_lfp  ! in [g DM m-2]

    dsl4jsb_Real2D_onChunk                :: forest_fract      ! in [%]
    dsl4jsb_Real2D_onChunk                :: grass_fract       ! in [%]

    dsl4jsb_Real2D_onChunk                :: fuel_bulk_density          ! in [kg m-3]
    dsl4jsb_Real2D_onChunk                :: moisture_extinction        ! in [-]
    dsl4jsb_Real2D_onChunk                :: fwd_wind_factor            ! in [-]
    dsl4jsb_Real2D_onChunk                :: lvg_moisture               ! in [0 to 2.5]%

    dsl4jsb_Real3D_onChunk                :: dead_fuel_consumption_1hr ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_consumption_10hr ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_consumption_100hr ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_consumption_1000hr ! out [g DM m-2]

    dsl4jsb_Real2D_onChunk                :: area_burned               ! out [ha ts-1]
    dsl4jsb_Real2D_onChunk                :: area_burned_fract         ! out [% ts-1]
    dsl4jsb_Real2D_onChunk                :: area_burned_output        ! out [ha s-1]
    dsl4jsb_Real2D_onChunk                :: area_burned_fract_output  ! out [% s-1]
    dsl4jsb_Real3D_onChunk                :: dead_frac_consumption     ! out [-]
    dsl4jsb_Real2D_onChunk                :: i_surface                 ! out [kW m-1]

    dsl4jsb_Real3D_onChunk                :: frac_consump_fcl          ! out [-]
    dsl4jsb_Real2D_onChunk                :: ic_pop_dens               ! in [#humans km-2]
    dsl4jsb_Real2D_onChunk                :: a2l_lightning_events      ! in [#events ha-1 d-1]
    dsl4jsb_Real2D_onChunk                :: a_nd                      ! in [#ignitions_indiv-1 d-1]
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_DFIRE_SF_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_DFIRE_SF_)
    dsl4jsb_Get_memory(Q_DFIRE_SF_)
    dsl4jsb_Get_memory(Q_DFIRE_FUEL_SF_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)              ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sat_sl)          ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)            ! in

    !< Ignitions
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, e_n_ignitions)          ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, e_n_humans)             ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, e_n_lightning)          ! out

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, t_air_min)              ! in [k]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, t_air_max)              ! in
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, precip_daily_mmd)       ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, precip_daily)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_10m)                   ! in

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, nesterov_index)         ! inout [Celsius 2]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, ab_moisture_diag)       ! out [-]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_SF_, fuel_moisture)          ! out [-]

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, fdi)                    ! out [%]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, ros_f)                  ! out [?]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, ros_b)                  ! out [?]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, mean_fire_area)         ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, fire_duration)          ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, distance_travelled)     ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, l2b_ratio)              ! out

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, pot_burn_area)          ! out [m2?]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, mask_is_fire_intense_enough) ! out [-]

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, dead_fuel_1hr_dfp)     ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, dead_fuel_10hr_dfp)    ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, dead_fuel_100hr_dfp)   ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, dead_fuel_1000hr_dfp)  ! in [g DM m-2]

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, live_grass_fuel_1hr_lfp)    ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, cured_live_grass_fuel_1hr_dfp) ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, live_tree_fuel_1hr_lfp)     ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, live_tree_fuel_10hr_lfp)    ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, live_tree_fuel_100hr_lfp)   ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_SF_, live_tree_fuel_1000hr_lfp)  ! in [g DM m-2]

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_SF_, forest_fract)        ! in [%]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_SF_, grass_fract)         ! in [%]

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_SF_, fuel_bulk_density)   ! in [kg m-3]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_SF_, moisture_extinction) ! in [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_SF_, fwd_wind_factor)     ! in [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_SF_, lvg_moisture)        ! in [0 to 2.5]%

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, reaction_intensity)       ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, prop_flux_ratio)          ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, wind_factor)              ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, eff_heating_number)       ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, preig_heat)               ! out

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, opt_react_velo)           ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, net_fuel_load)            ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, moist_damp_coeff)         ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, mineral_damp_coeff)       ! out

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_SF_, dead_fuel_consumption_1hr) ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_SF_, dead_fuel_consumption_10hr) ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_SF_, dead_fuel_consumption_100hr) ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_SF_, dead_fuel_consumption_1000hr) ! out [g DM m-2]

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, area_burned)               ! out [ha ts-1]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, area_burned_fract)         ! out [% ts-1]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, area_burned_output)        ! out [ha s-1]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, area_burned_fract_output)  ! out [% s-1]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_SF_, dead_frac_consumption)     ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, i_surface)                 ! out [kW m-1]

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_SF_, frac_consump_fcl)          ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, ic_pop_dens)               ! in [#humans km-2]
    dsl4jsb_Get_var2D_onChunk(A2L_, a2l_lightning_events)          ! in [#events ha-1 d-1]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, a_nd)                      ! in [#events ha-1 d-1]
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    ! ----------------------------------------------------------------------------------------------------- !

    no_children = tile%Get_no_of_children()

    !< Init variables (set output variables to 0)
    ab_moisture_diag = 0._wp
    fuel_moisture = 0._wp
    i_surface = 0._wp
    frac_consump_fcl = 0._wp
    area_burned_output = 0._wp
    area_burned_fract_output = 0._wp

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_STANDALONE__
    ALLOCATE(grid_cell_area(options%nc))
    grid_cell_area = 1._wp ! defined as fraction too, QS has site simulations
    burnable_area_ha = 1._wp
    mask_is_tile_enabled = .TRUE.
    pot_burn_area = 1._wp
#else
    ! IQ
    grid_cell_area => grid % area(ics:ice, iblk)
    tile_fract => tile%fract(ics:ice, iblk)
    pot_burn_area = grid_cell_area * (grass_fract + forest_fract)
    burnable_area_ha = sf_convert_units("m2", "ha", pot_burn_area)
    mask_is_tile_enabled = (tile_fract > 0._wp)
#endif

    sf_fuel_load = t_sf_fuel_load(nc, sf_fuel_dead_fuel_nr, sf_fuel_live_fuel_nr, df_sf_fc_full_nr, & ! params
                                  df_sf_fc_1hr_idx, df_sf_fc_100hr_idx, ixs_sf_dead_litter_fuel_pools, &
                    dead_fuel_1hr_dfp, dead_fuel_10hr_dfp, dead_fuel_100hr_dfp, dead_fuel_1000hr_dfp, &
                    live_tree_fuel_1hr_lfp, live_tree_fuel_10hr_lfp, live_tree_fuel_100hr_lfp, live_tree_fuel_1000hr_lfp, &
                    live_grass_fuel_1hr_lfp, &
                    cured_live_grass_fuel_1hr_dfp)

    sf_ros = t_sf_ros(nc, fwd_wind_factor, &
                             fuel_bulk_density, & ! in
                             ros_f, ros_b, & ! inout
                             distance_travelled, l2b_ratio, fire_duration, mean_fire_area, &
                             reaction_intensity, prop_flux_ratio, &
                             wind_factor, eff_heating_number, &
                             preig_heat, &
                             opt_react_velo, net_fuel_load, &
                             moist_damp_coeff, mineral_damp_coeff)

    sf_atm = t_sf_atm(t_air_min, t_air_max, wind_10m, precip_daily, precip_daily_mmd)
    sf_ignitions = t_sf_ignitions(nc, ic_pop_dens, a_nd, a2l_lightning_events, &
                                  e_n_humans, e_n_lightning, e_n_ignitions)

    !< Pack vegetation cover fraction
    veg_cover_frac = t_sf_veget_cover_frac(nc, forest_fract, grass_fract)

    ! todo: change to ab_moisture_diag
    ab_moisture_diag = calc_ab_moisture(soil_depth_sl(:,sl_ixs), wtr_soil_sl(:,sl_ixs), wtr_soil_sat_sl(:,sl_ixs), &
                                    mask=mask_is_tile_enabled)

    sf_fuel_moisture_all = t_sf_mf(nc, df_sf_fc_full_nr, moisture_extinction, & ! core
                                df_sf_alpha_nr, df_sf_alphas, ab_moisture_diag, fuel_moisture, & ! thonicke
                                df_sf_sigma_dd_nr, df_sf_sigma_dd_ros_nr, df_sf_sigma_lv_nr, & ! oberhagemann
                                df_alphas_dead, df_sigmas_dead, df_sigmas_live, & ! oberhagemann
                                lvg_moisture)

    sf_dynamics = t_sf_dynamics(nc, df_sf_fc_full_nr, sf_fuel_dead_fuel_nr, &
                                fdi, &
                                area_burned, area_burned_fract, &
                                dead_frac_consumption, i_surface)

    sf_config = t_sf_config(nc, df_sf_fc_spread_nr, df_sf_fc_full_nr, &
                            sf_fuel_dead_fuel_nr, &
                            dtime, one_day, &
                            dsl4jsb_Config(Q_DFIRE_SF_)%ros_model_id, &
                            df_sf_max_density_human_pop, &
                            df_sf_calorific_heat_content, df_sf_isurface_extinc_th)

    CALL calc_sf_dynamics(sf_config, & ! in...
                        burnable_area_ha, & ! soil properties
                        veg_cover_frac, & ! soil properties
                        sf_atm, & ! weather
                        sf_fuel_load, &
                        sf_ignitions, & ! inout...
                        sf_ros, & ! rate of spread
                        nesterov_index, & ! dryness
                        sf_fuel_moisture_all, &
                        sf_dynamics, &
                        mask=mask_is_tile_enabled)

    IF (dsl4jsb_Config(Q_DFIRE_SF_)%has_diagnostics) THEN
      CALL run_diagnostics_at_tile(nc, iblk, burnable_area_ha, sf_dynamics)
    ENDIF

    !< Ignitions into process memory
    dead_fuel_consumption_1hr = sf_dynamics % fc_dead_cls(:,:,df_sf_fc_1hr_idx)
    dead_fuel_consumption_10hr = sf_dynamics % fc_dead_cls(:,:,df_sf_fc_10hr_idx)
    dead_fuel_consumption_100hr = sf_dynamics % fc_dead_cls(:,:,df_sf_fc_100hr_idx)
    dead_fuel_consumption_1000hr = sf_dynamics % fc_dead_cls(:,:,df_sf_fc_1000hr_idx)

    mask_is_fire_intense_enough = logical2real(sf_dynamics % is_fire_intense_enough)
    area_burned_fract_output = sf_dynamics % area_burned_fract / dtime !< to seconds
    area_burned_output = sf_dynamics % area_burned / dtime !< to seconds
    ! todo: remove due to duplicity with 'dead_frac_consumption'
    frac_consump_fcl = sf_dynamics % dead_frac_consump

#ifdef __QUINCY_STANDALONE__
    DEALLOCATE(grid_cell_area)
#endif
  END SUBROUTINE update_q_dfire_sf_dynamics


  ! ======================================================================================================= !
  !>
  !> Run science
  !>
  SUBROUTINE calc_sf_dynamics(sf_config, &
                  grid_cell_area, & ! fuel properties
                  veg_cover_frac, & ! in
                  sf_atm, & ! in
                  sf_fuel_load, & ! in
                  sf_ignitions, & ! inout
                  sf_ros, & ! inout
                  nesterov_index, & ! inout
                  sf_fuel_moisture, & ! inout
                  spitfire, & ! inout
                  mask) ! in

    TYPE(t_sf_config), INTENT(in) :: sf_config
    REAL(wp), INTENT(in)  :: grid_cell_area(sf_config%nc)                    !< Gridcell area [ha]
    TYPE(t_sf_veget_cover_frac), INTENT(in) :: veg_cover_frac      !< Vegetation cover fraction [0 to 1]
    TYPE(t_sf_atm), INTENT(in) :: sf_atm                           !< Weather
    TYPE(t_sf_fuel_load), INTENT(in) :: sf_fuel_load               !< Fuel load

    TYPE(t_sf_ignitions), INTENT(inout)  :: sf_ignitions
    TYPE(t_sf_ros), INTENT(inout) :: sf_ros                        !< Rate of spread
    REAL(wp), INTENT(INOUT) :: nesterov_index(sf_config%nc)                  !< Dryness index [Celsius^2]
    TYPE(t_sf_mf), INTENT(inout) :: sf_fuel_moisture                  !< Fuel moisture
    TYPE(t_sf_dynamics), INTENT(inout)  :: spitfire                !< Spitfire dynamics
    LOGICAL, INTENT(in), OPTIONAL       :: mask(sf_config%nc)                !< Gridcell can burn

    LOGICAL :: mask_tile_enabled(sf_config%nc)

    ! locals
    LOGICAL                 :: has_fire_conditions(sf_config%nc)          !< 1, 10, 100 and 1000 hours

    mask_tile_enabled = .TRUE.
    IF (PRESENT(mask)) THEN
      mask_tile_enabled = mask
    ENDIF

    !< According to Thonicke et al. (2010), which is based on Rothermel et al.,
    !< fire spread depends on dead fuel availability. If there is no dead fuel,
    !< then fire cannot spread and we stop the calculation. The same applies if FDI is low.
    has_fire_conditions = mask_tile_enabled .AND. sf_fuel_load%has_dead_fuel

    !> Calculate whether there are conditions for fire based on the Nesterov Index,
    !> following Thonicke et al. (2010) with modifications from Pfeiffer et al. (2013)
    !> and own changes for QUINCY.
    nesterov_index = calc_sf_nesterov_index(sf_config%nc, sf_atm%t_air_max_c, sf_atm%t_air_min_c, &
        sf_atm%precip_daily_mmd, nesterov_index, &
        o_tsfactor = (sf_config%dtime / sf_config%one_day), &
        mask= has_fire_conditions)

    !> Fuel moisture and ROS
    !
    !> Calculate wind speed (paragraph before eq. 10 from Thonicke et al. (2010))
    sf_ros % u_f_wind = sf_ros % fwd_wind_factor * sf_atm%wind_sf

    IF (sf_config%ros_model_id == Q_DFIRE_ROS_THONICKE) THEN

      CALL calc_sf_dynamics_for_thonicke(sf_config, sf_fuel_load, &
                                    nesterov_index, &
                                    sf_fuel_moisture, & ! inout
                                    spitfire % fdi, & !< outputs
                                    sf_ros, &
                                    mask=has_fire_conditions) !< optional
    ELSE IF (sf_config%ros_model_id == Q_DFIRE_ROS_OBERHAGEMANN) THEN

      CALL calc_sf_dynamics_for_oberhagemann(sf_config, sf_fuel_load, &
                                    nesterov_index, &
                                    sf_fuel_moisture, & ! inout
                                    spitfire % fdi, & ! out
                                    sf_ros, &
                                    mask = has_fire_conditions) ! optional
    ELSE
      CALL finish("calc_sf_dynamics", "Unexpected fire sub-model found")
    END IF
    CALL sf_calc_rate_of_spread(sf_config%nc, &
                  sf_ros % u_f_wind, &
                  spitfire % fdi, &
                  veg_cover_frac, &
                  has_fire_conditions, &
                  sf_ros)

    !< Calculate fire ignition probability based on population and lightning, following Thonicke et al. (2010)
    CALL wrap_calc_sf_ignitions(sf_config, sf_ignitions)

    !< no ignitions events, no fire conditions
    has_fire_conditions = has_fire_conditions .AND. (sf_ignitions%e_n_ignitions > 0._wp)

    !< Based on ignitions and fire spread conditions, calculate burned area and its fraction.
    CALL calc_area_burned(sf_config%nc, sf_ignitions % e_n_ignitions, & !< in
                          spitfire % fdi, &
                          sf_ros % mean_fire_area, &
                          grid_cell_area, &
                          spitfire % area_burned, &  !< out
                          spitfire % area_burned_fract, & !< out
                          mask=has_fire_conditions) !< in, optional

    !< Fraction of consumption used to calculate dead fuel consumption (Thonicke et al. (2010),
    !< extracted from equations B1, B4, B5 and B6)
    spitfire % dead_frac_consump = sf_calc_frac_consump(sf_config%nc, sf_config%nfcl, sf_fuel_moisture%moisture_extinction, &
        sf_fuel_moisture % fuel_moisture_for_fc, mask=has_fire_conditions)

    !< Calculate dead fuel consumption for the 1h, 10h and 100h fuel classes,
    !< following Thonicke et al. (2010), Appendix B. and well as custom modifications.
    !< See mo_q_dfire_sf_fuel_consumption_class.f90 for furter details.
    spitfire % fc_dead_cls = sf_calc_dead_fuel_consumption(spitfire % area_burned_fract, &
                                                      spitfire % dead_frac_consump, &
                                                      sf_fuel_load%dead_fuel_fcl, mask=has_fire_conditions)

    !< According to Thonicke et al. (2010), fire intensity is calculated based on the fractional fuel combustion,
    !< rate of spread and burned area fraction (Sec. 2.2.5).
    !< Here, it is calculated only for grid cells where fuel is dry enough.
    spitfire % i_surface = sf_calc_fire_intensity(sf_config%nc, sf_config%nddpools, sf_config % nfcl_ros, &
                                                sf_config % const_calorific_heat_content, &
                                                sf_ros %  ros_forward, &
                                                spitfire % area_burned_fract, &
                                                spitfire % fc_dead_cls(:,:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx), &
                                                mask=has_fire_conditions)

    spitfire % is_fire_intense_enough = has_fire_conditions .AND. &
                                          (spitfire % i_surface >= sf_config % const_isurface_extinct_thd)

  END SUBROUTINE calc_sf_dynamics


  ! ======================================================================================================= !
  !>
  !> Calculate spitfire thonicke's implementation
  !>
  SUBROUTINE calc_sf_dynamics_for_thonicke(sf_config, sf_fuel_load, & ! in
                                    nesterov_index, &
                                    sf_mf, & ! inout
                                    fdi, & ! out
                                    sf_ros, &
                                    mask) !< optional

    !USE mo_q_dfire_sf_constants, ONLY:  ix_1hr => df_sf_fc_1hr_idx, &
    !                                   ix_10hr => df_sf_fc_10hr_idx, &
    !                                  ix_100hr => df_sf_fc_100hr_idx, &
    !                                 ix_1000hr => df_sf_fc_1000hr_idx

    TYPE(t_sf_config), INTENT(in) :: sf_config
    TYPE(t_sf_fuel_load), INTENT(in) :: sf_fuel_load
    REAL(wp), INTENT(in)        :: nesterov_index(:)      !< [deg C2]
    TYPE(t_sf_mf), INTENT(inout):: sf_mf
    REAL(wp), INTENT(out)       :: fdi(:)                 !< Fire dange index [%]
    TYPE(t_sf_ros), INTENT(inout) :: sf_ros
    LOGICAL, INTENT(in), OPTIONAL :: mask(:)              !< True: gridcell is enabled

    LOGICAL :: has_fire_conditions(sf_config%nc)

    fdi = 0._wp
    sf_mf % omegas_thon = 0._wp
    has_fire_conditions = .TRUE.
    IF (PRESENT(mask)) has_fire_conditions = mask

    CALL wrap_calc_sf_fuel_moisture_for_thonicke(sf_config, nesterov_index, sf_fuel_load, sf_mf)

    !> Fire danger is defined as zero in the absence of fuel or for wet fuel, and unity for completely dry fuel.
    CALL calc_fire_danger_index(sf_config % nc, sf_mf % omegas_thon(:, df_sf_omega_o_idx), sf_mf % moisture_extinction, fdi, mask=has_fire_conditions)

    !< It accounts for 1, 10 and 100 hours fuel classes
    has_fire_conditions = has_fire_conditions .AND. (fdi > 1.e-6_wp)

    CALL sf_calc_rate_of_spread_thon(sf_config % nc, sf_config % nfcl_ros, & ! in
                      sf_ros % fbd, &
                      sf_ros % u_f_wind, sf_mf % moisture_extinction, &
                      sf_fuel_load % dead_fuel_litter * df_sf_units_g_to_kg, &
                      sf_mf % omegas_thon(:,df_sf_omega_o_idx), &
                      has_fire_conditions, &
                      sf_ros) ! inout
  END SUBROUTINE calc_sf_dynamics_for_thonicke


  ! ======================================================================================================= !
  !>
  !> Calculate spitfire Oberhagemann's implementation
  !>
  SUBROUTINE calc_sf_dynamics_for_oberhagemann(sf_config, sf_fuel_load, & ! in
                  nesterov_index, &
                  fuel_moisture, & ! inout
                  fdi, & ! out
                  ros, &
                  mask) ! optional

    USE mo_q_dfire_sf_constants, ONLY: ixs_ros => df_sf_ix_sigma_dd_ros, & ! 1hr, 10hr, 100hr, cured live grass (dead)
                                    ix_clvg => df_sf_ix_sigma_dd_clvg, &
                                     ix_1hr => df_sf_ix_sigma_dd_1hr, &
                                    ix_10hr => df_sf_ix_sigma_dd_10hr, &
                                   ix_100hr => df_sf_ix_sigma_dd_100hr, &
                                  ix_1000hr => df_sf_ix_sigma_dd_1000hr, &
                                     ix_lvg => df_sf_ix_sigma_lv_lvg, &
                                     ix_wdy => df_sf_ix_sigma_lv_wdy, &
                                 nlitros => df_sf_sigma_dd_ros_litter_nr ! number of dead litter fuel classes involved in fire spread -1 to 100 hours-

    TYPE(t_sf_config), INTENT(in) :: sf_config
    TYPE(t_sf_fuel_load), INTENT(in) :: sf_fuel_load
    REAL(wp), INTENT(in) :: nesterov_index(:) !< Dryness index [deg C2]
    TYPE(t_sf_mf), INTENT(inout) :: fuel_moisture !< Fuel moisture
    REAL(wp), INTENT(out) :: fdi(:) !< Fire danger index [%]
    TYPE(t_sf_ros), INTENT(inout) :: ros !< Rate of spread [-]
    LOGICAL, INTENT(in) :: mask(:) !< True: gridcell enabled

    REAL(wp) :: woi_dd(sf_config%nc,fuel_moisture % sigma_dd_nr)        !< [kg DM m-2] (#nc, #sigma classes-1,10,100hrs,cured lvg-) Note: in kg
    REAL(wp) :: woi_lv(sf_config%nc,fuel_moisture % sigma_lv_nr)        !< [kg DM m-2] (#nc, #sigma classes-lvg, woody-) Note: in kg
    REAL(wp) :: me_dd(sf_config%nc)             !< Dead moisture of extinction [-] (#nc)
    REAL(wp) :: me_lv(sf_config%nc)             !< Live moisture of extinction [-] (#nc)

    fdi = 0._wp

    CALL wrap_calc_sf_fuel_moisture_for_oberhagemann(sf_config, nesterov_index, sf_fuel_load, fuel_moisture, mask)

    woi_dd = 0._wp
    woi_dd(:,ix_1hr) = sf_fuel_load % dead_fuel_litter(:,ix_1hr) * df_sf_units_g_to_kg
    woi_dd(:,ix_10hr) = sf_fuel_load % dead_fuel_litter(:,ix_10hr) * df_sf_units_g_to_kg
    woi_dd(:,ix_100hr) = sf_fuel_load % dead_fuel_litter(:,ix_100hr) * df_sf_units_g_to_kg
    woi_dd(:,ix_clvg) = sf_fuel_load % tot_dead_fuel_grass * df_sf_units_g_to_kg  !< 4 = dry livegrass

    woi_lv = 0._wp
    woi_lv(:,ix_lvg) = 0._wp !< 1 = wet livegrass, not yet implemented
    woi_lv(:,ix_wdy) = sf_fuel_load % tot_live_fuel_tree * df_sf_units_g_to_kg !< 2= woody fuel

    me_dd = fuel_moisture % moisture_extinction
    me_lv = calc_live_me(sf_config%nc, fuel_moisture % sigma_dd_ros_nr, fuel_moisture % sigma_lv_nr, &
                          fuel_moisture%dead(:,ixs_ros), &
                          me_dd, woi_dd(:,ixs_ros), woi_lv, &
                          fuel_moisture % sigmas_dead, fuel_moisture % sigmas_live, &
                          mask)

    !< only involved ROS fuel (and FM)
    CALL sf_calc_rate_of_spread_ober(&
        sf_config%nc, &
        fuel_moisture % sigma_dd_ros_nr, &
        fuel_moisture % sigma_lv_nr, &
        fuel_moisture % sigmas_dead, fuel_moisture % sigmas_live, &
        woi_dd(:,ixs_ros), woi_lv, & ! fuels
        fuel_moisture%dead(:,ixs_ros), fuel_moisture%live, &
        me_dd, me_lv, &
        ros % fbd, & !fuel_bulk_density, &
        ros % u_f_wind, &
        fuel_moisture%val, & ! output too
        mask, &
        ros)

    !> Fire danger is defined as zero in the absence of fuel or for wet fuel, and unity for completely dry fuel.
    CALL calc_fire_danger_index(sf_config%nc, fuel_moisture%val, fuel_moisture%moisture_extinction, fdi, mask=mask)

  END SUBROUTINE calc_sf_dynamics_for_oberhagemann

  ! ======================================================================================================= !
  !>
  !> Calculate Fire Danger Index (FDI)
  !> Eq. 8 from Thonicke et al. (2010)
  !>
  PURE SUBROUTINE calc_fire_danger_index(nc, fuel_moisture, moisture_extinction, fdi, mask)
    INTEGER, INTENT(in)   :: nc
    REAL(wp), INTENT(in)  :: fuel_moisture(nc)
    REAL(wp), INTENT(in)  :: moisture_extinction(nc)
    REAL(wp), INTENT(out) :: fdi(nc)
    LOGICAL, INTENT(in), OPTIONAL :: mask(nc)

    LOGICAL :: active_mask(nc)

    active_mask = .TRUE.
    IF (PRESENT(mask)) active_mask = mask

    fdi = 0._wp
    WHERE (active_mask .AND. moisture_extinction > 0._wp)
      fdi = MAX(0._wp, 1._wp - (fuel_moisture / moisture_extinction))
    END WHERE
  END SUBROUTINE calc_fire_danger_index

  ! ======================================================================================================= !
  !>
  !> Calculate fire intensity (1, 10 and 100 hours dead fuel classes)
  !> Eq. 15 from Thonicke et al. (2010)
  !>
  FUNCTION sf_calc_fire_intensity(nc, ndeadpools, nfuels, calorific_heat_content, ros_forward, area_burned_fract, &
      fc_dead_cls, mask) RESULT (i_surface)
    INTEGER, INTENT(in)  :: nc, ndeadpools, nfuels
    REAL(wp), INTENT(in) :: calorific_heat_content          !< Calorific heat content of fuel [KJ Kg-1]
    REAL(wp), INTENT(in) :: ros_forward(nc)                 !< Forward rate of spread of the fire [m min-1]
    REAL(wp), INTENT(in) :: area_burned_fract(nc)           !< Fraction of the area burned [0 to 1]
    REAL(wp), INTENT(in) :: fc_dead_cls(nc, ndeadpools, nfuels) !< Dead fuel consumption. It includes 1, 10 and 100 hours fuel classes [g DM m-2]
    LOGICAL, INTENT(in)  :: mask(nc)                        !< True = compute gridcell(s)
    REAL(wp)             :: i_surface(nc)                   !< Fire intensity [kW m-1]

    REAL(wp) :: tot_fc_dead(nc)                             !< Total dead fuel consumption aggregated over 1h, 10h, and 100h fuel classes [g DM m-2]
    INTEGER :: inc

    !< Account for fuel classes involved in the spread part
    tot_fc_dead = sf_calc_aggr_fuels(fc_dead_cls)
    i_surface = 0._wp

    DO inc = 1, nc
      IF (mask(inc) .AND. area_burned_fract(inc) > 0._wp) THEN
        !< Calculate the intensity of a surface fire at the flaming front (eq. 15 from Thonicke et al. (2010))
        i_surface(inc) = calorific_heat_content * (tot_fc_dead(inc) * ros_forward(inc)) / &
            (1000._wp * area_burned_fract(inc) * 60._wp)
      END IF
    END DO
  END FUNCTION sf_calc_fire_intensity


  ! ======================================================================================================= !
  !>
  !> Calculate area burned and fraction of the area burned of a gridcell
  !> Eq. 2 from Thonicke et al. (2010)
  !>
  SUBROUTINE calc_area_burned(nc, e_n_ignitions, fdi, mean_fire_area, area, ab, ab_fract, mask)
    INTEGER, INTENT(in)   :: nc
    REAL(wp), INTENT(in)  :: e_n_ignitions(nc)  !< Expected number of ignitions events [ha-1 ts-1]
    REAL(wp), INTENT(in)  :: fdi(nc)            !< Fire danger index [0 to 1]
    REAL(wp), INTENT(in)  :: mean_fire_area(nc) !< [ha]
    REAL(wp), INTENT(in)  :: area(nc)           !< Area of the gridcell [ha]
    REAL(wp), INTENT(out) :: ab(nc)             !< area_burned [ha]
    REAL(wp), INTENT(out) :: ab_fract(nc)       !< area_burned fraction [-]
    LOGICAL, INTENT(in), OPTIONAL :: mask(nc)   !< True=compute gridcell(s)
    LOGICAL :: mask_sel(nc)
    INTEGER :: ic

    ab = 0._wp
    ab_fract = 0._wp

    mask_sel = .TRUE. !< default
    IF (PRESENT(mask)) THEN
      mask_sel = mask
    END IF

    DO ic=1, nc
      IF (mask_sel(ic)) THEN
        !< Calculate Fraction of the Area burned (eq. 2 from Thonicke et al. (2010))
        ab(ic) = e_n_ignitions(ic) * fdi(ic) * mean_fire_area(ic) * area(ic) !< [ha TS-1]
        ab_fract(ic) = ab(ic) / area(ic)
      END IF
    END DO
  END SUBROUTINE calc_area_burned


  ! ======================================================================================================= !
  !>
  !> Calculate the above ground relative soil moisture
  !>
  FUNCTION calc_ab_moisture(ab_soil_depth_sl, ab_wtr_soil_sl, ab_wtr_wsat_sl, mask) RESULT (ab_moisture)
    REAL(wp), INTENT(in)                :: ab_soil_depth_sl(:,:)      !< soil depth [m] (#nc, #above ground soil layers)
    REAL(wp), INTENT(in)                :: ab_wtr_soil_sl(:,:)        !< soil water content [m] (#nc, #above ground soil layers)
    REAL(wp), INTENT(in)                :: ab_wtr_wsat_sl(:,:)        !< soil water saturation [m] (#nc, #above ground soil layers)
    LOGICAL, INTENT(in), OPTIONAL       :: mask(:)                    !< True = gridcell enabled
    REAL(wp), ALLOCATABLE               :: ab_moisture(:)             !< Relative above ground soil moisture [-] (nc)

    !< Select the first 2 layers from the soil moisture and calculate soil moisture content as a fraction of field capacity
    INTEGER :: is, ic
    INTEGER :: nc, nsl
    REAL(wp), ALLOCATABLE :: ab_moisture_sl(:,:)
    LOGICAL, ALLOCATABLE :: mask_in(:)

    nc = SIZE(ab_wtr_wsat_sl, DIM=1)
    nsl = SIZE(ab_wtr_wsat_sl, DIM=2) ! 2=soil layers index
    ALLOCATE(ab_moisture_sl, mold=ab_wtr_soil_sl)

    ALLOCATE(mask_in(nc))
    IF (PRESENT(mask)) THEN
      mask_in=mask
    ELSE
      mask_in = .TRUE.
    ENDIF

    !< Relative soil moisture for each above ground soil layers
    ab_moisture_sl = 0._wp
    DO is=1, nsl
      DO ic=1, nc
        IF (mask_in(ic) .AND. ab_wtr_wsat_sl(ic,is)>0._wp) THEN
          ab_moisture_sl(ic,is) = ab_wtr_soil_sl(ic,is) / ab_wtr_wsat_sl(ic,is)
        END IF
      END DO
    END DO

    !< Average of above ground soil layers
    IF (.NOT. ALLOCATED(ab_moisture)) ALLOCATE(ab_moisture(nc))
    ab_moisture = 0._wp
    DO ic=1, nc
      IF (mask_in(ic)) THEN
        ab_moisture(ic) = stats_average_weight(nsl, ab_moisture_sl(ic,:), ab_soil_depth_sl(ic,:))
      END IF
    ENDDO
  END FUNCTION calc_ab_moisture



  ! ======================================================================================================= !
  !>
  !> Check calculated data is consistent, otherwise report and stop the simulation
  !>
  SUBROUTINE run_diagnostics_at_tile(nc, iblk, grid_cell_area, sf_outputs)
    INTEGER, INTENT(in) :: nc, iblk
    REAL(wp), INTENT(in) :: grid_cell_area(:)
    TYPE(t_sf_dynamics), INTENT(in) :: sf_outputs

    LOGICAL :: test_BAfract_is_outside_limits
    LOGICAL :: test_BA_is_bigger_than_gridcellarea
    LOGICAL :: has_any_test_failed
    INTEGER :: ic
    INTEGER :: nfailedtests

    nfailedtests = 0
    test_BAfract_is_outside_limits = .FALSE.
    test_BA_is_bigger_than_gridcellarea = .FALSE.
    has_any_test_failed = .FALSE.

    ! AB_fract is within fractional limits? 0 to 1
    IF (ANY(sf_outputs % area_burned_fract > 1._wp)) THEN
      DO ic=1, nc
        IF (sf_outputs % area_burned_fract(ic) > 1._wp) THEN
          WRITE(*,*) "spitfire.f90:: ic, iblk, AB_fract=", ic, iblk, &
                  sf_outputs % area_burned_fract(ic)
          test_BAfract_is_outside_limits = .TRUE.
        ENDIF
      ENDDO
      IF (test_BAfract_is_outside_limits) nfailedtests = nfailedtests + 1
      WRITE(*,*) "(spitfire) AB_fract is beyong fractional limits (0 to 1)?", test_BAfract_is_outside_limits
    ENDIF

    !< is Area burnt bigger than grid cell area?
    DO ic=1, nc
      IF (sf_outputs % area_burned(ic) > grid_cell_area(ic)) THEN
        WRITE(*,*) "spitfire.f90:: ic, iblk, AB, grid cell area=", ic, iblk, &
                  sf_outputs % area_burned(ic), grid_cell_area(ic)
        test_BA_is_bigger_than_gridcellarea = .TRUE.
      ENDIF
      ENDDO
    IF (test_BA_is_bigger_than_gridcellarea) THEN
      nfailedtests = nfailedtests + 1
      WRITE(*,*) "(spitfire) AB is bigger than gridcell area?", test_BA_is_bigger_than_gridcellarea
    ENDIF

    ! sum up
    has_any_test_failed = test_BAfract_is_outside_limits .OR. &
                          test_BA_is_bigger_than_gridcellarea

    IF (has_any_test_failed) THEN
      WRITE(*,*) "Failed tests:", nfailedtests
      CALL finish("run_diagnostics", "Test(s) failed")
    ENDIF

  END SUBROUTINE run_diagnostics_at_tile

#endif
END MODULE mo_q_dfire_sf_dynamics
