!>
!>  @par Copyright
!>    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!>    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany
!>    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!>
!>    This file is part of the QUINCY model.
!>    The QUINCY model is free software: you can redistribute it and/or modify
!>    it under the terms of the GNU General Public License version 3,
!>    as published by the Free Software Foundation.
!>
!>    The QUINCY model is distributed in the hope that it will be useful,
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!>    GNU General Public License for more details.
!>
!>    You should have received a copy of the GNU General Public License
!>    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!>
!> @brief module to determine the spitfire fire fuel
!>
!> Main task to compute the spitfire
!>
!> -update_dfire_spitfire: integrate Quincy with the spitfire process:
!>    - convert units from Quincy into spitfire
!>    - gather lctlib parameters
!>    - select precalculated daily precipitation method
!>    - Soil moisture content
!>    - Fraction of the vegetation
!>    - Call calc_spitfire (science)
!>    - Dead feedback
!>    - Live feedback
!>    - Save data into process memory
!>
!> -calc_spitfire: spitfire subroutine which calls
!>    - [mask] method (file)
!>    - [-] Initialize t_sf_spitfire, used to output data
!>    - [-] Calc ignition events (mo_difre_sf_ignition.f90)
!>    - [-] Calc fuel load (mo_q_dfire_sf_fuel_class.f90)
!>    - [-] Calc Nesterov index (mo_q_dfire_sf_nesterov.f90)
!>    - [-] Calc fuel moisture (mo_q_dfire_sf_fmoist_class.f90)
!>    - [-] Calc Fire danger index (-)
!>    - [dead_fuel-all- > 0 & FDI>0] Calc Rate of spread (mo_q_dfire_sf_ros_class.f90)
!>    - [dead_fuel-all- > 0 & FDI>0] Calc Area Burned and its fraction (-)
!>    - [dead_fuel-1,10,100- > 0 & FDI>0] Calc dead fraction of consumption (mo_q_dfire_sf_frac_consump_class.f90)
!>    - [dead_fuel-1,10,100- > 0 & FDI>0]  Calc dead fuel consumption (mo_q_dfire_sf_fuel_consumption.f90)
!>    - [dead_fuel-1,10,100- > 0 & FDI>0]  Calc fire intensity (-)
!>    - [... & I_surface>50] Calc tree mortality (mo_q_dfire_sf_damage.f90)
!>    - [... & I_surface>50]  Calc live fuel consumption (mo_q_dfire_sf_fuel_consumption.f90)
!>    - [... & I_surface>50]  Calc emissions (mo_q_dfire_sf_emissions_class.f90)
!>    - [... & I_surface>50]  Derive fuel consumption for N,P,C13,C14,N15
!>    - [... & I_surface>50]  Calc residue from emissions
!>
!> mo_q_dfire_spitfire.f90 subroutines:
!>    -sf_calc_fire_intensity: calculate fire intensity
!>    -calc_area_burned_and_fract: calculate area burned and its fraction
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_q_dfire_spitfire
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_lnd_bgcm_store,        ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,  ONLY: SB_BGCM_LOSS_ID, SB_BGCM_POOL_ID
  USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract

  USE mo_q_dfire_constants,       ONLY: df_sf_sl_dims_nr, sl_ixs => df_sf_sl_indices !< Soil related constants: top soil moisture layers

  USE mo_q_dfire_constants,       ONLY: t_sf_fuel_dead_fuel_woody_ix
  USE mo_q_dfire_constants,       ONLY: df_sf_max_density_human_pop
  USE mo_q_dfire_constants,       ONLY: df_sf_fc_spread_nr
  USE mo_q_dfire_constants,       ONLY: df_sf_alpha_nr, df_sf_alphas !< Drying parameters

  USE mo_q_dfire_constants,       ONLY: df_sf_calorific_heat_content, df_sf_isurface_extinc_th

  USE mo_q_dfire_constants,       ONLY: sf_fuel_dead_fuel_nr, sf_fuel_live_fuel_nr, df_sf_fc_full_nr, &
                                      df_sf_fc_1hr_idx, df_sf_fc_10hr_idx, df_sf_fc_100hr_idx, df_sf_fc_1000hr_idx !< fuel classes indices

  USE mo_q_dfire_constants,       ONLY: df_sf_omegas_nr, & !< Fuel moisture dimension size (for omega_o and omega_nl)
                                      df_sf_omega_o_idx, & !< omega_o index from sf_calc_fuel_moisture output (array)
                                      df_sf_omega_nl_idx

  USE mo_q_dfire_constants,       ONLY: ixs_sf_dead_litter_fuel_pools, ixs_sf_dead_livegrass_fuel_pools
  USE mo_q_dfire_constants,       ONLY: df_sf_units_g_to_kg
  USE mo_q_dfire_constants,       ONLY: df_sf_sigma_dd_nr, df_sf_sigma_lv_nr, df_sf_sigma_dd_ros_nr
  USE mo_q_dfire_constants,       ONLY: df_sigmas_dead, df_sigmas_live, df_alphas_dead
  USE mo_q_dfire_config_class,    ONLY: Q_DFIRE_ROS_THONICKE, Q_DFIRE_ROS_OBERHAGEMANN

  USE mo_q_dfire_sf_fuel_class,   ONLY: t_sf_fuel_class, sf_calc_has_any_fuel, sf_calc_aggr_fuels, t_sf_fuel_load
  USE mo_q_dfire_sf_fmoist_class, ONLY: sf_calc_fuel_moisture, t_mf, sf_calc_fuel_moisture_ober
  USE mo_q_dfire_sf_ignition_class,     ONLY: calc_q_sf_ignition, t_q_sf_ignitions
  USE mo_q_dfire_sf_nesterov_class,     ONLY: calc_sf_nesterov_index
  USE mo_q_dfire_sf_ros_class,          ONLY: t_sf_veget_cover_frac, sf_calc_rate_of_spread
  USE mo_q_dfire_sf_ros_thon_class,     ONLY: sf_calc_rate_of_spread_thon
  USE mo_q_dfire_sf_ros_ober_class,     ONLY: sf_calc_rate_of_spread_ober
  USE mo_q_dfire_sf_ros_type,           ONLY: t_sf_ros
  USE mo_q_dfire_sf_ros_ober_class,     ONLY: calc_live_me
  USE mo_q_dfire_sf_frac_consump_class, ONLY: sf_calc_frac_consump
  USE mo_q_dfire_sf_fuel_consump_class, ONLY: sf_calc_dead_fuel_consumption

  USE mo_q_dfire_sf_adaptor,    ONLY: sf_convert_units
  USE mo_q_dfire_util,          ONLY: logical2real
  USE mo_quincy_stats,          ONLY: stats_average_weight

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_q_dfire_spitfire, sf_calc_fire_intensity, t_sf_spitfire, calc_ab_moisture

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_spitfire'

  !<
  !< Calculated Spitfire values
  !< Move spitfire data from calc_spitfire <-> update_q_dfire_spitfire subroutines
  !<
  TYPE t_sf_spitfire
    REAL(wp), ALLOCATABLE :: e_n_ignitions(:)
    REAL(wp), ALLOCATABLE :: e_n_humans(:)
    REAL(wp), ALLOCATABLE :: e_n_lightning(:)

    REAL(wp), ALLOCATABLE :: fdi(:) !< fire danger index [0 to 1]
    TYPE(t_sf_ros)        :: ros    !< fire dynamics

    REAL(wp), ALLOCATABLE :: area_burned(:)             !<
    REAL(wp), ALLOCATABLE :: area_burned_fract(:)

    REAL(wp), ALLOCATABLE :: dead_frac_consump(:,:)     !< Dead fraction of consumption (#nc, #fuel classes) [0 to 1]
    REAL(wp), ALLOCATABLE :: fc_dead_cls(:,:,:)         !< dead fuel consumption (#nc, #dead pools, #fuels classes) [g DM m-2]
    REAL(wp), ALLOCATABLE :: i_surface(:)               !< Fire intensity [kW m-1]

    LOGICAL, ALLOCATABLE  :: is_fire_intense_enough(:)  !< mask, TRUE= fire intensity is strong enough to exist
  END TYPE t_sf_spitfire

  INTERFACE t_sf_spitfire
    MODULE PROCEDURE :: create_sf_spitfire
  END INTERFACE t_sf_spitfire

!< Group weather related variables to pass in calc_spitfire
TYPE t_sf_met
    REAL(wp), POINTER :: wind(:)                    !< Wind speed [m s-1]
    REAL(wp), POINTER :: t_air_min_c(:)             !< Daily minimum air temperature [deg C]
    REAL(wp), POINTER :: t_air_max_c(:)             !< Daily maximum air temperature [deg C]
    REAL(wp), POINTER :: precip_daily_mmd(:)        !< Daily precipitation [mm/day]
    REAL(wp), POINTER :: nesterov_index(:)          !< Dryness index [deg C2]
END TYPE t_sf_met

INTERFACE t_sf_met
  MODULE PROCEDURE :: create_sf_met
END INTERFACE t_sf_met

!< Group ignition related variables to pass in calc_spitfire
TYPE :: t_sf_ignitions
    REAL(wp), POINTER :: a_nd(:)!< Propensity of people to produce ignition events [ignitions individual-1 d-1]
    REAL(wp), POINTER :: pop_dens(:) !< Average human population density in the grid cell [individuals km-2]
    REAL(wp), POINTER :: lightning_events(:) !< Lightning frequency [ha-1 d-1]
    REAL(wp)          :: max_pop_dens
END TYPE t_sf_ignitions

INTERFACE t_sf_ignitions
  MODULE PROCEDURE :: create_sf_ignitions
END INTERFACE t_sf_ignitions

!< Group fuel moisture related variables to pass in calc_spitfire
TYPE :: t_sf_fuel_moisture
    REAL(wp), POINTER :: alphas_thon(:) !< Drying parameter for 1, 10 and 100 hours fuel classes
    REAL(wp), POINTER :: alphas_dead_ober(:)
    REAL(wp), POINTER :: sigmas_dead_ober(:)
    REAL(wp), POINTER :: sigmas_live_ober(:)
    REAL(wp), POINTER :: ab_moisture(:) !< above ground relative soil moisture [0 to 1]
    REAL(wp), POINTER :: lvg_moisture(:) !< live grass fuel moisture [%], note: up to 250%
    REAL(wp), POINTER :: moisture_extinction(:) !< Moisture of extinction threshold [0 to 1]
    REAL(wp), POINTER :: fuel_moisture_thon(:)
    TYPE(t_mf), ALLOCATABLE :: fuel_moisture_ober(:)
END TYPE t_sf_fuel_moisture

INTERFACE t_sf_fuel_moisture
  MODULE PROCEDURE :: create_sf_fuel_moisture
END INTERFACE t_sf_fuel_moisture

CONTAINS

  !>
  !> t_sf_spitfire constructor: Allocate and initialize
  !>
  FUNCTION create_sf_fuel_moisture(nc, nddsigmas, nlvsigmas, &
                                    alphas_thon, &
                                    alphas_dead_ober, sigmas_dead_ober, sigmas_live_ober, &
                                    ab_moisture, lvg_moisture, &
                                    moisture_extinction, &
                                    fuel_moisture_thon, fuel_moisture_ober) &
                                    RESULT (sf_fuel_moisture)
    INTEGER, INTENT(in) :: nc
    INTEGER, INTENT(in) :: nddsigmas
    INTEGER, INTENT(in) :: nlvsigmas
    REAL(wp), POINTER, INTENT(in) :: alphas_thon(:)
    REAL(wp), POINTER, INTENT(in) :: alphas_dead_ober(:)
    REAL(wp), POINTER, INTENT(in) :: sigmas_dead_ober(:)
    REAL(wp), POINTER, INTENT(in) :: sigmas_live_ober(:)
    REAL(wp), POINTER, INTENT(in) :: ab_moisture(:)
    REAL(wp), POINTER, INTENT(in) :: lvg_moisture(:)
    REAL(wp), POINTER, INTENT(in) :: moisture_extinction(:)
    REAL(wp), POINTER, INTENT(in) :: fuel_moisture_thon(:)
    TYPE(t_mf), INTENT(in) :: fuel_moisture_ober(:)
    TYPE(t_sf_fuel_moisture) :: sf_fuel_moisture

    sf_fuel_moisture % alphas_thon => alphas_thon
    sf_fuel_moisture % alphas_dead_ober => alphas_dead_ober
    sf_fuel_moisture % sigmas_dead_ober => sigmas_dead_ober
    sf_fuel_moisture % sigmas_live_ober => sigmas_live_ober
    sf_fuel_moisture % ab_moisture => ab_moisture
    sf_fuel_moisture % lvg_moisture => lvg_moisture
    sf_fuel_moisture % moisture_extinction => moisture_extinction
    sf_fuel_moisture % fuel_moisture_thon => fuel_moisture_thon
    ALLOCATE(sf_fuel_moisture % fuel_moisture_ober)
    sf_fuel_moisture % fuel_moisture_ober = t_mf(nc, nddsigmas, nlvsigmas)
  END FUNCTION create_sf_fuel_moisture


  ! ======================================================================================================= !
  !>
  !> t_sf_spitfire constructor: Allocate and initialize
  !>
  FUNCTION create_sf_spitfire(nc, nfuels, ndeadpools) RESULT (sf_spitfire)
    INTEGER, INTENT(in) :: nc          !< #Gridcells
    INTEGER, INTENT(in) :: nfuels      !< # fuel classes
    INTEGER, INTENT(in) :: ndeadpools  !< # dead fuel pools
    TYPE(t_sf_spitfire) :: sf_spitfire !< Output

    ALLOCATE(sf_spitfire % e_n_ignitions(nc))
    ALLOCATE(sf_spitfire % e_n_humans(nc))
    ALLOCATE(sf_spitfire % e_n_lightning(nc))
    ALLOCATE(sf_spitfire % fdi(nc))
    ALLOCATE(sf_spitfire % area_burned(nc))
    ALLOCATE(sf_spitfire % area_burned_fract(nc))
    ALLOCATE(sf_spitfire % dead_frac_consump(nc, nfuels))
    ALLOCATE(sf_spitfire % fc_dead_cls(nc, ndeadpools, nfuels)) ! spitfire dead fuel consumption (fuel classes)
    ALLOCATE(sf_spitfire % i_surface(nc))
    ALLOCATE(sf_spitfire % is_fire_intense_enough(nc))

    sf_spitfire % e_n_ignitions = 0._wp
    sf_spitfire % e_n_humans = 0._wp
    sf_spitfire % e_n_lightning = 0._wp
    sf_spitfire % fdi = 0._wp
    sf_spitfire % ros = t_sf_ros(nc)
    sf_spitfire % area_burned = 0._wp
    sf_spitfire % area_burned_fract = 0._wp
    sf_spitfire % dead_frac_consump = 0._wp
    sf_spitfire % fc_dead_cls = 0._wp
    sf_spitfire % i_surface = 0._wp
    sf_spitfire % is_fire_intense_enough = .FALSE.

  END FUNCTION create_sf_spitfire


  ! ======================================================================================================= !
  !>
  !> Constructor for t_sf_ignitions
  !>
  FUNCTION create_sf_ignitions(a_nd, ic_pop_dens, lightning_events, max_pop_dens) RESULT (sf_ignitions)
    REAL(wp), INTENT(in), POINTER :: a_nd(:)
    REAL(wp), INTENT(in), POINTER :: ic_pop_dens(:)
    REAL(wp), INTENT(in), POINTER :: lightning_events(:)
    REAL(wp), INTENT(in)          :: max_pop_dens
    TYPE(t_sf_ignitions) :: sf_ignitions

    sf_ignitions % a_nd => a_nd
    sf_ignitions % ic_pop_dens => ic_pop_dens
    sf_ignitions % lightning_events => lightning_events
    sf_ignitions % max_pop_dens = max_pop_dens
  END FUNCTION create_sf_ignitions


  ! ======================================================================================================= !
  !>
  !> Constructor for t_sf_met
  !>
  FUNCTION create_sf_met(wind, tmin, tmax, precip, ni) RESULT (sf_met)
    REAL(wp), INTENT(in), POINTER :: wind(:)
    REAL(wp), INTENT(in), POINTER :: tmin(:)
    REAL(wp), INTENT(in), POINTER :: tmax(:)
    REAL(wp), INTENT(in), POINTER :: precip(:)
    REAL(wp), INTENT(in), POINTER :: ni(:)
    TYPE(t_sf_met) :: sf_met

    sf_met % wind => wind
    sf_met % t_air_min_c => tmin
    sf_met % t_air_max_c => tmax
    sf_met % precip_daily_mmd => precip
    sf_met % nesterov_index => ni
  END FUNCTION create_sf_met


  ! ======================================================================================================= !
  !>
  !> Compute spitfire for current time step
  !> It links Quincy and Spitfire
  !>
  SUBROUTINE update_q_dfire_spitfire(tile, options)
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_process_class,       ONLY: VEG_, Q_DFIRE_, HYDRO_, Q_DFIRE_FUEL_, A2L_
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_grid,                ONLY: Get_vgrid, Get_grid
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid, t_jsb_grid
    USE mo_jsb_math_constants,      ONLY: one_day

    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(Q_DFIRE_)
    dsl4jsb_Use_memory(Q_DFIRE_FUEL_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(VEG_)
    !dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(Q_DFIRE_)
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
    CHARACTER(len=*), PARAMETER           :: routine = TRIM(modname)//':update_q_fdire_sf_spitfire'

    REAL(wp), POINTER                     :: grid_cell_area(:) ! [m2]
    REAL(wp), POINTER                     :: tile_fract(:) ! [%]
    LOGICAL                               :: mask_is_tile_enabled(options%nc)

    REAL(wp)                              :: ab_moisture(options%nc)                      !< top relative soil moisture layer [-]

    TYPE(t_sf_veget_cover_frac)           :: veg_cover_frac                               !< Vegetation cover
    TYPE(t_sf_fuel_load)                  :: fuel_load
    TYPE(t_sf_ignitions)                  :: ignitions
    TYPE(t_sf_met)                        :: met
    TYPE(t_sf_fuel_moisture)              :: fuel_moisture_params

    REAL(wp), TARGET                      :: wind_sf(options%nc)                          !< Wind into spitfire units [m min-1]
    REAL(wp), TARGET                      :: t_air_min_c(options%nc)                      !< Minimum daily temperature [Celsius]
    REAL(wp), TARGET                      :: t_air_max_c(options%nc)                      !< Maximum daily temperature [Celsius]
    REAL(wp)                              :: burnable_area_ha(options%nc)                 !< Burnable area [ha]

    TYPE(t_sf_spitfire)                   :: sf_spitfire

    INTEGER                               :: nr_of_children

    REAL(wp), TARGET                      :: a_nd(options%nc)
    INTEGER                               :: ierror
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_DFIRE_)
    dsl4jsb_Def_memory(Q_DFIRE_FUEL_)
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
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_DFIRE_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_DFIRE_)
    dsl4jsb_Get_memory(Q_DFIRE_)
    dsl4jsb_Get_memory(Q_DFIRE_FUEL_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sat_sl)        ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)          ! in

    !< Ignitions
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, e_n_ignitions)          ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, e_n_humans)             ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, e_n_lightning)          ! out

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, t_air_min)              ! in [k]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, t_air_max)              ! in
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, precip_daily_mmd)       ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, precip_daily)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_10m)                 ! in

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, nesterov_index)         ! inout [Celsius 2]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, ab_moisture_diag)       ! out [-]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_, fuel_moisture)          ! out [-]

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, fdi)                    ! out [%]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, ros_f)                  ! out [?]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, ros_b)                  ! out [?]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, mean_fire_area)         ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, fire_duration)          ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, distance_travelled)     ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, l2b_ratio)              ! out

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, pot_burn_area)          ! out [m2?]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, mask_is_fire_intense_enough) ! out [-]

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, dead_fuel_1hr_dfp)     ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, dead_fuel_10hr_dfp)    ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, dead_fuel_100hr_dfp)   ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, dead_fuel_1000hr_dfp)  ! in [g DM m-2]

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_grass_fuel_1hr_lfp)    ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, cured_live_grass_fuel_1hr_dfp) ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1hr_lfp)     ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_10hr_lfp)    ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_100hr_lfp)   ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1000hr_lfp)  ! in [g DM m-2]

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, forest_fract)      ! in [%]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, grass_fract)       ! in [%]

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, fuel_bulk_density)   ! in [kg m-3]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, moisture_extinction) ! in [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, fwd_wind_factor)     ! in [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, lvg_moisture)        ! in [0 to 2.5]%

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, reaction_intensity)     ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, prop_flux_ratio)        ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, wind_factor)            ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, eff_heating_number)     ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, preig_heat)             ! out

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, opt_react_velo)         ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, net_fuel_load)          ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, moist_damp_coeff)       ! out
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, mineral_damp_coeff)     ! out

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_, dead_fuel_consumption_1hr) ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_, dead_fuel_consumption_10hr) ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_, dead_fuel_consumption_100hr) ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_, dead_fuel_consumption_1000hr) ! out [g DM m-2]

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, area_burned)               ! out [ha ts-1]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, area_burned_fract)         ! out [% ts-1]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, area_burned_output)        ! out [ha s-1]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, area_burned_fract_output)  ! out [% s-1]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_, dead_frac_consumption)     ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, i_surface)                 ! out [kW m-1]

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_, frac_consump_fcl)            ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, ic_pop_dens)                 ! in [#humans km-2]
    dsl4jsb_Get_var2D_onChunk(A2L_, a2l_lightning_events) ! in [#events ha-1 d-1]
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    ! ----------------------------------------------------------------------------------------------------- !

    no_children = tile%Get_no_of_children()

    ! !< Init variables
    fuel_load = t_sf_fuel_load(nc, sf_fuel_dead_fuel_nr, sf_fuel_live_fuel_nr, df_sf_fc_full_nr, &
                                 dead_fuel_1hr_dfp, &
                                 dead_fuel_10hr_dfp, &
                                 dead_fuel_100hr_dfp, &
                                 dead_fuel_1000hr_dfp, &
                                 cured_live_grass_fuel_1hr_dfp, &
                                 live_grass_fuel_1hr_lfp, &
                                 live_tree_fuel_1hr_lfp, &
                                 live_tree_fuel_10hr_lfp, &
                                 live_tree_fuel_100hr_lfp, &
                                 live_tree_fuel_1000hr_lfp)

    !< Make sure those values are set to 0, as it means there is no fire.
    ab_moisture_diag = 0._wp
    fuel_moisture = 0._wp
    i_surface = 0._wp
    ab_moisture = 0._wp
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

    !> Units conversion
    t_air_max_c = sf_convert_units("Kelvin", "Celsius", t_air_max)
    t_air_min_c = sf_convert_units("Kelvin", "Celsius", t_air_min)
    precip_daily_mmd = sf_convert_units("kg m-2 s-1", "mm/day", precip_daily)
    wind_sf = sf_convert_units("m/s", "m/min", wind_10m) !< wind to spitfire units
    met = t_sf_met(wind_sf, t_air_min_c, t_air_max_c, precip_daily_mmd, nesterov_index)

    !< Pack vegetation cover fraction
    veg_cover_frac = t_sf_veget_cover_frac(nc, forest_fract, grass_fract)

    ab_moisture = calc_ab_moisture(soil_depth_sl(:,sl_ixs), wtr_soil_sl(:,sl_ixs), wtr_soil_sat_sl(:,sl_ixs), &
                                    mask=mask_is_tile_enabled)

    ! Load input data: rate of human-caused ignitions per person [#ign indiv-1 d-]
    ! Todo: further research on this data, currently taken from a paper.
    a_nd = a2l_lightning_events * 4._wp
    ignitions = t_sf_ignitions(a_nd, ic_pop_dens, a2l_lightning_events, df_sf_max_density_human_pop)

    fuel_moisture_params = t_sf_fuel_moisture(nc, df_sf_sigma_dd_nr, df_sf_sigma_lv_nr, &
                                        df_sf_alphas, & ! thonicke
                                        df_alphas_dead, & ! oberhagemann
                                        df_sigmas_dead, & ! oberhagemann
                                        df_sigmas_live, & ! oberhagemann
                                        ab_moisture, lvg_moisture, &
                                        moisture_extinction, fuel_moisture)

    CALL calc_spitfire(nc, df_sf_alpha_nr, df_sf_fc_full_nr, dtime, one_day, &
                        dsl4jsb_Config(Q_DFIRE_)%ros_model_id, &
                        fuel_bulk_density, fwd_wind_factor, & ! lctlib weighted
                        df_sf_calorific_heat_content, &
                        df_sf_isurface_extinc_th, &
                        fuel_moisture_params, &
                        burnable_area_ha, &
                        ignitions, &
                        fuel_load, &
                        met, &
                        veg_cover_frac, &
                        fuel_moisture, & ! fuel properties
                        sf_spitfire, &
                        mask=mask_is_tile_enabled)

    dead_fuel_consumption_1hr = sf_spitfire % fc_dead_cls(:,:,df_sf_fc_1hr_idx)
    dead_fuel_consumption_10hr = sf_spitfire % fc_dead_cls(:,:,df_sf_fc_10hr_idx)
    dead_fuel_consumption_100hr = sf_spitfire % fc_dead_cls(:,:,df_sf_fc_100hr_idx)
    dead_fuel_consumption_1000hr = sf_spitfire % fc_dead_cls(:,:,df_sf_fc_1000hr_idx)

    !< Ignitions into process memory
    e_n_ignitions = sf_spitfire % e_n_ignitions
    e_n_humans = sf_spitfire % e_n_humans
    e_n_lightning = sf_spitfire % e_n_lightning

    dead_frac_consumption = sf_spitfire % dead_frac_consump
    i_surface = sf_spitfire % i_surface
    mask_is_fire_intense_enough = logical2real(sf_spitfire % is_fire_intense_enough)
    area_burned_fract = sf_spitfire % area_burned_fract
    area_burned = sf_spitfire % area_burned
    area_burned_fract_output = sf_spitfire % area_burned_fract / dtime !< to seconds
    area_burned_output = sf_spitfire % area_burned / dtime !< to seconds

    fdi = sf_spitfire % fdi ! todo: define as pointer

    ab_moisture_diag(:) =  ab_moisture(:)

    ros_f = sf_spitfire % ros % ros_forward
    ros_b = sf_spitfire % ros % ros_backwards

    mean_fire_area = sf_spitfire % ros % mean_fire_area
    fire_duration = sf_spitfire % ros % fire_duration
    distance_travelled = sf_spitfire % ros % distance_travelled
    l2b_ratio = sf_spitfire % ros % l2b_ratio

    reaction_intensity(:) = sf_spitfire % ros % reaction_intensity(:)
    prop_flux_ratio(:) = sf_spitfire % ros % prop_flux_ratio(:)
    wind_factor(:) = sf_spitfire % ros % wind_factor(:)
    eff_heating_number(:) = sf_spitfire % ros % eff_heating_number(:)
    preig_heat(:) = sf_spitfire % ros % pre_ign_heat(:)

    opt_react_velo(:) = sf_spitfire % ros % opt_react_velo(:)
    net_fuel_load(:) = sf_spitfire % ros % net_fuel_load(:)
    moist_damp_coeff(:) = sf_spitfire % ros % moist_damp_coeff(:)
    mineral_damp_coeff(:) = sf_spitfire % ros % mineral_damp_coeff(:)

    frac_consump_fcl = sf_spitfire % dead_frac_consump
#ifdef __QUINCY_STANDALONE__
    DEALLOCATE(grid_cell_area)
#endif
  END SUBROUTINE update_q_dfire_spitfire


  ! ======================================================================================================= !
  !>
  !> Run science
  !>
  SUBROUTINE calc_spitfire(nc, nfuelclasses, dtime, one_day, & ! in
                  config_ros_model_id, &
                  lctlib_fuel_bulk_density, lctlib_fwd_wind_factor, &
                  const_calorific_heat_content, const_i_surface_thd, &
                  ! alphas, &
                  fuel_moisture_params, &
                  grid_cell_area, & ! fuel properties
                  ignitions, &
                  fuel_load, &
                  met, &
                  veg_cover_frac, & ! in
                  fuel_moisture, & ! out
                  spitfire, & ! out
                  mask)

    USE mo_q_dfire_constants, ONLY:  ix_1hr => df_sf_ix_sigma_dd_1hr, &
                                    ix_10hr => df_sf_ix_sigma_dd_10hr, &
                                   ix_100hr => df_sf_ix_sigma_dd_100hr, &
                                  ix_1000hr => df_sf_ix_sigma_dd_1000hr

    INTEGER, INTENT(in)   :: nc, nfuelclasses
    REAL(wp), INTENT(in)  :: dtime, one_day

    INTEGER, INTENT(in)   :: config_ros_model_id                   !< Select ROS implementation

    REAL(wp), INTENT(in)  :: lctlib_fuel_bulk_density(nc)          !< Fuel bulk density [kg m-3]
    REAL(wp), INTENT(in)  :: lctlib_fwd_wind_factor(nc)            !< Forward wind factor [0 to 1]

    REAL(wp), INTENT(in)  :: const_calorific_heat_content          !< Fuel calorific heat content [kJ kg-1]
    REAL(wp), INTENT(in)  :: const_i_surface_thd                   !< Fire intensity threshold [kW m-1]

    TYPE(t_sf_fuel_moisture), INTENT(inout) :: fuel_moisture_params

    REAL(wp), INTENT(in)  :: grid_cell_area(nc)                    !< Gridcell area [ha]

    TYPE(t_sf_ignitions), INTENT(in) :: ignitions !< Ignition related inputs
    TYPE(t_sf_fuel_load), INTENT(in) :: fuel_load !< [g DM m-2]
    TYPE(t_sf_met), INTENT(inout) :: met

    TYPE(t_sf_veget_cover_frac), INTENT(in) :: veg_cover_frac      !< Vegetation cover fraction [0 to 1]

    REAL(wp), INTENT(out) :: fuel_moisture(nc,df_sf_omegas_nr)     !< Fuel moisture (#nc, #omega_nr (omega_o, omega_nl))

    TYPE(t_sf_spitfire), INTENT(out)    :: spitfire                !< Spitfire outputs
    LOGICAL, INTENT(in), OPTIONAL       :: mask(nc)

    LOGICAL :: mask_tile_enabled(nc)

    ! locals
    TYPE(t_q_sf_ignitions)  :: sf_ignitions !< Ignitions outputs
    INTEGER :: ndeadpools
    LOGICAL  :: has_dead_fuel(nc)           !< True = dead fuel
    LOGICAL  :: has_fire_conditions(nc)     !< 1, 10, 100 and 1000 hours

    !< Thonicke related
    REAL(wp), ALLOCATABLE :: fuel_moisture_for_fc(:,:) !< fuel moisture (or omega) [%] (#nc,#fuel classes-1,10,100,1000hrs-)
    !< Oberhagemann related
    TYPE(t_mf), ALLOCATABLE :: fuel_moisture_ober
    !< ROS related
    REAL(wp) :: u_f_wind(nc)  !< Wind speed at vegetation cover [m min-1]

    INTEGER                               :: i, j

    ! todo: deal with tile_frac = 0?
    ndeadpools = fuel_load%count_dead_fuel_pools()
    mask_tile_enabled = .TRUE.
    IF (PRESENT(mask)) THEN
      mask_tile_enabled = mask
    ENDIF

    !< According to Thonicke et al. (2010), which is based on Rothermel et al.,
    !< fire spread depends on dead fuel availability. If there is no dead fuel,
    !< then fire cannot spread and we stop the calculation. The same applies if FDI is low.
    has_dead_fuel = fuel_load%has_dead_fuel(spread=.TRUE.)
    !< It accounts for 1, 10 and 100 hours fuel classes
    has_fire_conditions = mask_tile_enabled .AND. has_dead_fuel

    sf_ignitions = t_q_sf_ignitions(nc)
    spitfire = t_sf_spitfire(nc, ndeadpools, sf_fuel_live_fuel_nr)

    !> Calculate whether there are conditions for fire based on the Nesterov Index,
    !> following Thonicke et al. (2010) with modifications from Pfeiffer et al. (2013)
    !> and own changes for QUINCY.
    met%nesterov_index = calc_sf_nesterov_index(nc, met%t_air_max_c, met%t_air_min_c, &
                                                met%precip_daily_mmd, met%nesterov_index, &
                                                o_tsfactor = (dtime / one_day), &
                                                mask= has_fire_conditions)

    !> Fuel moisture and ROS
    u_f_wind = lctlib_fwd_wind_factor * met%wind_sf
    !< Calculate amount of dead fuel available for combustion following Thonicke et al. (2010)
    IF (config_ros_model_id == Q_DFIRE_ROS_THONICKE) THEN

      CALL calc_spitfire_thonicke(nc, SIZE(fuel_moisture_params%alphas_thon), df_sf_fc_spread_nr, &
                                    fuel_moisture_params%alphas_thon, &
                                    lctlib_fuel_bulk_density, &
                                    u_f_wind, &
                                    met%nesterov_index, &
                                    fuel_moisture_params%ab_moisture, &
                                    fuel_moisture_params%moisture_extinction, &
                                    fuel_load, &
                                    spitfire % fdi, & !< outputs
                                    fuel_moisture_params%fuel_moisture_thon, &
                                    spitfire%ros, &
                                    mask=has_fire_conditions) !< optional
    ELSE IF (config_ros_model_id == Q_DFIRE_ROS_OBERHAGEMANN) THEN

      CALL calc_spitfire_oberhagemann(nc, &
          df_sf_sigma_dd_nr, & ! all dead fuel classes
          df_sf_sigma_dd_ros_nr, & ! dead fuel classes involved in fire spread
          df_sf_sigma_lv_nr, & ! live fuel classes
          fuel_moisture_params%alphas_dead_ober, &
          fuel_moisture_params%sigmas_dead_ober, &
          fuel_moisture_params%sigmas_live_ober, &
          met%nesterov_index, &
          lctlib_fuel_bulk_density, &
          fuel_moisture_params%moisture_extinction, &
          u_f_wind, &
          fuel_load, &
          fuel_moisture_params%lvg_moisture, & ! live grass moisture pre-calculated at fuel.f90
          fuel_moisture_params%fuel_moisture_ober, & ! outputs
          spitfire % fdi, &
          spitfire % ros, &
          mask = has_fire_conditions) ! optional

       !< used in thonicke's implementation
       fuel_moisture_params%fuel_moisture_thon = 0._wp
    ELSE
      CALL finish("calc_spitfire", "Unexpected fire sub-model found")
    END IF
    CALL sf_calc_rate_of_spread(nc, &
                  u_f_wind, &
                  spitfire % fdi, &
                  veg_cover_frac, &
                  has_fire_conditions, &
                  spitfire % ros)

    !< Calculate fire ignition probability based on population and lightning, following Thonicke et al. (2010)
    sf_ignitions = calc_q_sf_ignition(nc, ignitions%pop_dens, &
                                    ignitions%a_nd, ignitions%lightning_events, &
                                    ignitions%max_pop_dens)
    !< From daily to ts
    spitfire%e_n_ignitions = sf_ignitions%total * (dtime / one_day)
    spitfire%e_n_humans = sf_ignitions%human * (dtime / one_day)
    spitfire%e_n_lightning = sf_ignitions%lightning * (dtime / one_day)

    !< no ignitions events, no fire conditions
    has_fire_conditions = has_fire_conditions .AND. (spitfire%e_n_ignitions > 0._wp)

    !< Based on ignitions and fire spread conditions, calculate burned area and its fraction.
    CALL calc_area_burned(nc, spitfire % e_n_ignitions, & !< in
                          spitfire % fdi, &
                          spitfire % ros % mean_fire_area, &
                          grid_cell_area, &
                          spitfire % area_burned, &  !< out
                          spitfire % area_burned_fract, & !< out
                          mask=has_fire_conditions) !< in, optional

    IF (config_ros_model_id == Q_DFIRE_ROS_THONICKE) THEN
      ALLOCATE(fuel_moisture_for_fc(nc,nfuelclasses))
      ASSOCIATE(w_o => fuel_moisture_params%fuel_moisture_thon(:,df_sf_omega_o_idx), &
                w_nl => fuel_moisture_params%fuel_moisture_thon(:,df_sf_omega_nl_idx))
        fuel_moisture_for_fc(:,ix_1hr) = w_nl ! 1hr
        fuel_moisture_for_fc(:,ix_10hr) = w_o ! 10hrs
        fuel_moisture_for_fc(:,ix_100hr) = w_o ! 100hrs
        fuel_moisture_for_fc(:,ix_1000hr) = w_o ! 1000hrs
      END ASSOCIATE
    ELSE IF (config_ros_model_id == Q_DFIRE_ROS_OBERHAGEMANN) THEN
      ALLOCATE(fuel_moisture_for_fc(nc,nfuelclasses))
      fuel_moisture_for_fc(:,ix_1hr) = fuel_moisture_params%fuel_moisture_ober%dead(:,ix_1hr)
      fuel_moisture_for_fc(:,ix_10hr) = fuel_moisture_params%fuel_moisture_ober%dead(:,ix_10hr)
      fuel_moisture_for_fc(:,ix_100hr) = fuel_moisture_params%fuel_moisture_ober%dead(:,ix_100hr)
      fuel_moisture_for_fc(:,ix_1000hr) = fuel_moisture_params%fuel_moisture_ober%dead(:,ix_1000hr)
    ELSE
      CALL finish("calc_spitfire", "Unexpected fire sub-model found")
    END IF

    !< Fraction of consumption used to calculate dead fuel consumption (Thonicke et al. (2010),
    !< extracted from equations B1, B4, B5 and B6)
    spitfire % dead_frac_consump = sf_calc_frac_consump(nc, nfuelclasses, &
        fuel_moisture_params%moisture_extinction, &
        fuel_moisture_for_fc, mask=has_fire_conditions)

    !< Calculate dead fuel consumption for the 1h, 10h and 100h fuel classes,
    !< following Thonicke et al. (2010), Appendix B. and well as custom modifications.
    !< See mo_q_dfire_sf_fuel_consumption_class.f90 for furter details.
    spitfire % fc_dead_cls = sf_calc_dead_fuel_consumption(spitfire % area_burned_fract, &
                                                      spitfire % dead_frac_consump, &
                                                      dead_fuel, mask=has_fire_conditions)

    !< According to Thonicke et al. (2010), fire intensity is calculated based on the fractional fuel combustion,
    !< rate of spread and burned area fraction (Sec. 2.2.5).
    !< Here, it is calculated only for grid cells where fuel is dry enough.
    spitfire % i_surface = sf_calc_fire_intensity(nc, ndeadpools, df_sf_fc_spread_nr, &
                                                const_calorific_heat_content, &
                                                spitfire % ros %  ros_forward, &
                                                spitfire % area_burned_fract, &
                                                spitfire % fc_dead_cls(:,:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx), &
                                                mask=has_fire_conditions)

    spitfire % is_fire_intense_enough = has_fire_conditions .AND. &
                                          (spitfire % i_surface >= const_i_surface_thd)

  END SUBROUTINE calc_spitfire


  ! ======================================================================================================= !
  !>
  !> Calculate spitfire thonicke's implementation
  !>
  SUBROUTINE calc_spitfire_thonicke(nc, nalphas, nspread, &
                                    alphas, &
                                    fuel_bulk_density, &
                                    u_f_wind, &
                                    nesterov_index, ab_moisture, &
                                    moisture_extinction, &
                                    fuel_load, &
                                    fdi, fuel_moisture, & !< output
                                    ros, &
                                    mask) !< optional

    INTEGER, INTENT(in)         :: nc, nalphas, nspread
    REAL(wp), INTENT(in)        :: alphas(:) !< Surface area to volume ratio [cm2/cm3]
    REAL(wp), INTENT(in)        :: fuel_bulk_density(:) !< Fuel bulk density [kg m-3]
    REAL(wp), INTENT(in)        :: u_f_wind(:) !< Wind speed at vegetation cover [m min-1]
    REAL(wp), INTENT(in)        :: nesterov_index(:) !< [deg C2]
    REAL(wp), INTENT(in)        :: ab_moisture(:) !< Aboveground soil moisture [%]
    REAL(wp), INTENT(in)        :: moisture_extinction(:) !< Dead moisture of extinction [%]
    TYPE(t_sf_fuel_load), INTENT(in) :: fuel_load

    REAL(wp), INTENT(out)       :: fdi(:) !< Fire dange index [%]
    REAL(wp), INTENT(out)       :: fuel_moisture(:,:) !< Fuel moisture (#nc,#fuel moisture types) [%]
    TYPE(t_sf_ros), INTENT(out) :: ros
    LOGICAL, INTENT(in), OPTIONAL :: mask(:) !< True: gridcell is enabled
    REAL(wp), ALLOCATABLE       :: dead_litter(:,:) !< (#nc,#fuel classes)

    LOGICAL :: has_fire_conditions(nc)

    fdi = 0._wp
    fuel_moisture = 0._wp
    has_fire_conditions = .TRUE.
    IF (PRESENT(mask)) has_fire_conditions = mask

    !< Compute fuel moisture (omega_o, omega_nl)
    CALL fuel_load%Calc_dead_litter(.TRUE., .TRUE., a2=dead_litter)
    fuel_moisture = sf_calc_fuel_moisture(nc, nalphas, nspread, alphas, &
                                          nesterov_index, ab_moisture, &
                                          dead_litter, &
                                          fuel_load%Calc_total_cured_live_grass())

    !> Fire danger is defined as zero in the absence of fuel or for wet fuel, and unity for completely dry fuel.
    !> Eq. 8 from Thonicke et al. (2010)
    WHERE(moisture_extinction > 0._wp)
      fdi = max(0._wp, 1._wp  - (fuel_moisture(:, df_sf_omega_o_idx) / moisture_extinction))
    ENDWHERE

    !< It accounts for 1, 10 and 100 hours fuel classes
    has_fire_conditions = has_fire_conditions .AND. (fdi > 1.e-6_wp)

    ros = sf_calc_rate_of_spread_thon(nc, nspread, &
                      fuel_bulk_density, &
                      u_f_wind, moisture_extinction, fdi, &
                      dead_litter * df_sf_units_g_to_kg, &
                      fuel_moisture(:,df_sf_omega_o_idx), &
                      has_fire_conditions)
  END SUBROUTINE calc_spitfire_thonicke


  ! ======================================================================================================= !
  !>
  !> Calculate spitfire Oberhagemann's implementation
  !>
  !>
  !>
  SUBROUTINE calc_spitfire_oberhagemann(nc, &
                  ndd, & ! 5=1,10,100,1000 hours, cured live grass
                  nddros, & ! 4 = 1,10,100, cured live grass
                  nlv, & ! 2=live grass, woody
                  alphas_dd, &
                  sigmas_dd, sigmas_lv, &
                  nesterov_index, &
                  fuel_bulk_density, &
                  moisture_extinction, &
                  u_f_wind, &
                  fuel_load, &
                  lvg_moisture, &
                  fuel_moisture, & ! outputs
                  fdi, &
                  ros, &
                  mask) ! optional

    USE mo_q_dfire_constants, ONLY: ixs_ros => df_sf_ix_sigma_dd_ros, &
                                    ix_clvg => df_sf_ix_sigma_dd_clvg, &
                                     ix_1hr => df_sf_ix_sigma_dd_1hr, &
                                    ix_10hr => df_sf_ix_sigma_dd_10hr, &
                                   ix_100hr => df_sf_ix_sigma_dd_100hr, &
                                     ix_lvg => df_sf_ix_sigma_lv_lvg, &
                                     ix_wdy => df_sf_ix_sigma_lv_wdy, &
                                 nlitros => df_sf_sigma_dd_ros_litter_nr ! number of dead litter fuel classes involved in fire spread -1 to 100 hours-

    INTEGER, INTENT(in) :: nc
    INTEGER, INTENT(in) :: ndd !< number of fuel classes (1 to 1000 hours ...)
    INTEGER, INTENT(in) :: nddros !< dead fuel classes involved in fire spread (1 to 100 hours and cured live grass)
    INTEGER, INTENT(in) :: nlv !< live fuel classes (live grass and woody)
    REAL(wp), INTENT(in) :: alphas_dd(:) !< Drying parameter [deg C-2]
    REAL(wp), INTENT(in) :: sigmas_dd(:), sigmas_lv(:) !< Surface area to volume ratio [cm2/cm3]
    REAL(wp), INTENT(in) :: nesterov_index(:) !< Dryness index [deg C2]
    REAL(wp), INTENT(in) :: fuel_bulk_density(:) !< Fuel bulk density [kg m-3]
    REAL(wp), INTENT(in) :: moisture_extinction(:) !< Dead moisture of extinction [%]
    REAL(wp), INTENT(in) :: u_f_wind(:) !< Wind speed at vegetation cover [m min-1]
    TYPE(t_sf_fuel_load), INTENT(in) :: fuel_load(:)
    REAL(wp), INTENT(in) :: lvg_moisture(:) !< [%], Note: up to 250%
    TYPE(t_mf), INTENT(out) :: fuel_moisture !< Fuel moisture
    REAL(wp), INTENT(out) :: fdi(:) !< Fire danger index [%]
    TYPE(t_sf_ros), INTENT(out) :: ros !< Rate of spread [-]
    LOGICAL, OPTIONAL :: mask(:) !< True: gridcell enabled

    REAL(wp), ALLOCATABLE :: dead_fuel_litter(:,:) !< [g DM m-2] (#nc,#Spread fuel classes-1,10,100 hours-)
    REAL(wp) :: woi_dd(nc,ndd)    !< [kg DM m-2] (#nc, #sigma classes-1,10,100hrs,cured lvg-) Note: in kg
    REAL(wp) :: woi_lv(nc,nlv)    !< [kg DM m-2] (#nc, #sigma classes-lvg, woody-) Note: in kg
    REAL(wp) :: me_dd(nc)            !< Dead moisture of extinction [-] (#nc)
    REAL(wp) :: me_lv(nc)            !< Live moisture of extinction [-] (#nc)

    fdi = 0._wp
    CALL fuel_load%Calc_dead_litter(.TRUE., .TRUE., a2=dead_fuel_litter)

    fuel_moisture = sf_calc_fuel_moisture_ober(nc, &
                      ndd, & ! 5 = 1, 10, 100, 1000 hours and cured live grass
                      nlitros, & ! 3=1, 10 and 100 hours
                      nlv, &
                      alphas_dd, &
                      lvg_moisture, &
                      dead_fuel_litter, &
                      fuel_load%Calc_total_cured_live_grass(), &
                      nesterov_index, &
                      mask)

    woi_dd = 0._wp
    woi_dd(:,ix_1hr) = dead_fuel_litter(:,ix_1hr) * df_sf_units_g_to_kg
    woi_dd(:,ix_10hr) = dead_fuel_litter(:,ix_10hr) * df_sf_units_g_to_kg
    woi_dd(:,ix_100hr) = dead_fuel_litter(:,ix_100hr) * df_sf_units_g_to_kg
    woi_dd(:,ix_clvg) = fuel_load%Calc_total_cured_live_grass() * df_sf_units_g_to_kg  !< 4 = dry livegrass

    woi_lv = 0._wp
    woi_lv(:,ix_lvg) = 0._wp !< 1 = wet livegrass, not yet implemented
    woi_lv(:,ix_wdy) = fuel_load%Calc_total_live_tree() * df_sf_units_g_to_kg !< 2= woody fuel

    me_dd = moisture_extinction
    me_lv = calc_live_me(nc, nddros, nlv, fuel_moisture%dead(:,ixs_ros), &
                          me_dd, woi_dd(:,ixs_ros), woi_lv, sigmas_dd, sigmas_lv, mask)

    !< only involved ROS fuel (and FM)
    ros = sf_calc_rate_of_spread_ober(&
        nc, nddros, nlv, &
        sigmas_dd, sigmas_lv, &
        woi_dd(:,ixs_ros), woi_lv, & ! fuels
        fuel_moisture%dead(:,ixs_ros), fuel_moisture%live, &
        me_dd, me_lv, &
        fuel_bulk_density, &
        u_f_wind, &
        fuel_moisture%val, & ! output too
        mask) ! optional

    !> Fire danger is defined as zero in the absence of fuel or for wet fuel, and unity for completely dry fuel.
    !> Eq. 8 from Thonicke et al. (2010)
    WHERE(moisture_extinction > 0._wp)
      fdi = max(0._wp, 1._wp  - (fuel_moisture%val / moisture_extinction))
    ENDWHERE

  END SUBROUTINE calc_spitfire_oberhagemann

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
      IF (mask(ic)) THEN
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

#endif
END MODULE mo_q_dfire_spitfire
