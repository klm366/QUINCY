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

  USE mo_q_dfire_constants,       ONLY: df_sf_sl_dims_nr, df_sf_sl_indices !< Soil related constants: top soil moisture layers 

  USE mo_q_dfire_constants,       ONLY: t_sf_fuel_dead_fuel_woody_ix
  USE mo_q_dfire_constants,       ONLY: df_sf_max_density_human_pop
  USE mo_q_dfire_constants,       ONLY: df_sf_fc_spread_nr
  USE mo_q_dfire_constants,       ONLY: df_sf_alpha_nr, df_sf_alphas !< Drying parameters

  USE mo_q_dfire_constants,       ONLY: df_sf_calorific_heat_content, df_sf_isurface_extinc_th

  USE mo_q_dfire_constants,       ONLY: sf_fuel_dead_fuel_nr, sf_fuel_live_fuel_nr, df_sf_fc_full_nr, &
                                      df_sf_fc_1hr_idx, df_sf_fc_10hr_idx, df_sf_fc_100hr_idx, df_sf_fc_1000hr_idx !< fuel classes indices

  USE mo_q_dfire_constants,       ONLY: df_sf_omegas_nr, & !< Fuel moisture dimension size (for omega_o and omega_nl)
                                      df_sf_omega_o_idx !< omega_o index from sf_calc_fuel_moisture output (array)

  USE mo_q_dfire_sf_fuel_class,   ONLY: t_sf_fuel_class, sf_calc_has_any_fuel, sf_calc_aggr_fuels
  USE mo_q_dfire_sf_fmoist_class, ONLY: sf_calc_fuel_moisture
  USE mo_q_dfire_sf_ignition_class,   ONLY: calc_q_sf_ignition, t_q_sf_ignitions
  USE mo_q_dfire_sf_nesterov_class,   ONLY: calc_sf_nesterov_index
  USE mo_q_dfire_sf_ros_class,    ONLY: t_sf_ros, t_sf_veget_cover_frac, sf_calc_rate_of_spread
  USE mo_q_dfire_sf_frac_consump_class, ONLY: sf_calc_frac_consump
  USE mo_q_dfire_sf_fuel_consump_class, ONLY: sf_calc_dead_fuel_consumption

  USE mo_q_dfire_sf_adaptor,    ONLY: sf_convert_units
  USE mo_q_dfire_util,          ONLY: logical2real
  USE mo_quincy_stats,          ONLY: stats_average_weight

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_q_dfire_spitfire, sf_calc_fire_intensity, t_sf_spitfire

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


  TYPE t_sf_bgcm_2L3D_p ! todo: remove
    REAL(wp), POINTER :: p(:,:,:,:) !< dims: (part,elem,nc,layer)
  END TYPE t_sf_bgcm_2L3D_p

CONTAINS

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

    REAL(wp), POINTER                     :: grid_cell_area(:) ! [m2?]
    REAL(wp), POINTER                     :: tile_fract(:) ! [%]
    REAL(wp), ALLOCATABLE                 :: tile_area(:) ! [m2]
    LOGICAL                               :: mask_is_tile_enabled(options%nc)

    REAL(wp)                              :: ab_soil_depth_sl(options%nc,df_sf_sl_dims_nr)    !< above ground soil depth layers [m]
    REAL(wp)                              :: ab_w_soil_sl(options%nc,df_sf_sl_dims_nr)    !< above ground soil water [m]
    REAL(wp)                              :: ab_w_soil_fc_sl(options%nc,df_sf_sl_dims_nr) !< above ground soil field capacity [m]
    REAL(wp)                              :: ab_moisture(options%nc)                      !< top relative soil moisture layer [-]

    TYPE(t_sf_veget_cover_frac)           :: veg_cover_frac                               !< Vegetation cover

    !< dead fuel (#nc, #dead fuel pools, #fuel classes) [g DM m-2]
    REAL(wp)                              :: dead_fuel_fcl(options%nc,sf_fuel_dead_fuel_nr,df_sf_fc_full_nr)
    !< live grass fuel (#nc, #live fuel pools, #fuel classes) [g DM m-2]
    REAL(wp)                              :: live_grass_fuel_fcl(options%nc,sf_fuel_live_fuel_nr,df_sf_fc_full_nr)
    !< live tree fuel (#nc, #live fuel pools, #fuel classes) [g DM m-2]
    REAL(wp)                              :: live_tree_fuel_fcl(options%nc,sf_fuel_live_fuel_nr,df_sf_fc_full_nr)

    REAL(wp)                              :: wind_sf(options%nc)                          !< Wind into spitfire units [m min-1]
    REAL(wp)                              :: t_air_min_c(options%nc)                      !< Minimum daily temperature [Celsius]
    REAL(wp)                              :: t_air_max_c(options%nc)                      !< Maximum daily temperature [Celsius]
    REAL(wp)                              :: grid_cell_area_ha(options%nc)                !< Grid cell area [ha] 

    TYPE(t_sf_spitfire)                   :: sf_spitfire

    INTEGER                               :: nr_of_children

    REAL(wp)                              :: a_nd(options%nc)
    INTEGER                               :: ierror
    ! ----------------------------------------------------------------------------------------------------- !
    !dsl4jsb_Def_mt2L2D :: veg_pool_mt
    !dsl4jsb_Def_mt2L3D :: sb_loss_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_DFIRE_)
    dsl4jsb_Def_memory(Q_DFIRE_FUEL_)
    dsl4jsb_Def_memory(VEG_)
    !dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Local vars here
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real3D_onChunk                :: wtr_soil_sl        ! soil water content [m]
    dsl4jsb_Real3D_onChunk                :: wtr_soil_fc_sl     ! soil water field capacity [m]
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

    dsl4jsb_Real3D_onChunk                :: live_tree_fuel_1hr_lfp     ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_tree_fuel_10hr_lfp    ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_tree_fuel_100hr_lfp   ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_tree_fuel_1000hr_lfp  ! in [g DM m-2]

    dsl4jsb_Real2D_onChunk                :: forest_fract      ! in [%]
    dsl4jsb_Real2D_onChunk                :: grass_fract       ! in [%]

    dsl4jsb_Real2D_onChunk                :: fuel_bulk_density          ! in [kg m-3]
    dsl4jsb_Real2D_onChunk                :: moisture_extinction        ! in [-]
    dsl4jsb_Real2D_onChunk                :: fwd_wind_factor            ! in [-]

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
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_fc_sl)         ! in
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
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1hr_lfp)     ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_10hr_lfp)    ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_100hr_lfp)   ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1000hr_lfp)  ! in [g DM m-2]

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, forest_fract)      ! in [%]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, grass_fract)       ! in [%]

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, fuel_bulk_density)   ! in [kg m-3]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, moisture_extinction) ! in [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, fwd_wind_factor)     ! in [-]

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

    !< Init variables
    live_grass_fuel_fcl = 0._wp

    ! rebuild fuel variables
    dead_fuel_fcl(:,:,df_sf_fc_1hr_idx) = dead_fuel_1hr_dfp ! dead fuel
    dead_fuel_fcl(:,:,df_sf_fc_10hr_idx) = dead_fuel_1hr_dfp
    dead_fuel_fcl(:,:,df_sf_fc_100hr_idx) = dead_fuel_1hr_dfp
    dead_fuel_fcl(:,:,df_sf_fc_1000hr_idx) = dead_fuel_1hr_dfp
    live_grass_fuel_fcl(:,:,df_sf_fc_1hr_idx) = live_grass_fuel_1hr_lfp ! live grass fuel
    live_tree_fuel_fcl(:,:,df_sf_fc_1hr_idx) = live_tree_fuel_1hr_lfp ! live tree fuel
    live_tree_fuel_fcl(:,:,df_sf_fc_10hr_idx) = live_tree_fuel_10hr_lfp
    live_tree_fuel_fcl(:,:,df_sf_fc_100hr_idx) = live_tree_fuel_100hr_lfp
    live_tree_fuel_fcl(:,:,df_sf_fc_1000hr_idx) = live_tree_fuel_1000hr_lfp

    !< Init variables
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
    ALLOCATE(tile_area(options%nc))
    grid_cell_area = 1._wp ! defined as fraction too, QS has site simulations
    grid_cell_area_ha = 1._wp
    mask_is_tile_enabled = .TRUE.
    tile_area = 1._wp
#else
    ! IQ
    grid_cell_area => grid % area(ics:ice, iblk)
    tile_fract => tile%fract(ics:ice, iblk)
    tile_area = grid_cell_area * tile_fract
    grid_cell_area_ha = sf_convert_units("m2", "ha", tile_area)
    mask_is_tile_enabled = (tile_fract > 0._wp)
#endif
    ! diagnostics only
    pot_burn_area = grid_cell_area * (grass_fract + forest_fract)

    !> Units conversion
    t_air_max_c = sf_convert_units("Kelvin", "Celsius", t_air_max)
    t_air_min_c = sf_convert_units("Kelvin", "Celsius", t_air_min)
    precip_daily_mmd = sf_convert_units("kg m-2 s-1", "mm/day", precip_daily)

    wind_sf = sf_convert_units("m/s", "m/min", wind_10m) !< wind to spitfire units

    !< Pack vegetation cover fraction
    veg_cover_frac = t_sf_veget_cover_frac(nc, forest_fract, grass_fract)

    ab_soil_depth_sl(:,:) = soil_depth_sl(:, [df_sf_sl_indices])
    ab_w_soil_sl(:,:) = wtr_soil_sl(:, [df_sf_sl_indices])
    ab_w_soil_fc_sl(:,:) = wtr_soil_fc_sl(:, [df_sf_sl_indices])
    ab_moisture = calc_ab_moisture(ab_soil_depth_sl, ab_w_soil_sl, ab_w_soil_fc_sl, mask=mask_is_tile_enabled)

    ! Load input data
    a_nd = a2l_lightning_events * 4._wp

    CALL calc_spitfire(nc, df_sf_alpha_nr, df_sf_fc_full_nr, dtime, one_day, &
                        fuel_bulk_density, fwd_wind_factor, & ! lctlib weighted
                        ic_pop_dens, & ! ignition events
                        a_nd, &
                        a2l_lightning_events, &
                        df_sf_calorific_heat_content, &
                        df_sf_isurface_extinc_th, &
                        df_sf_alphas, &
                        moisture_extinction, &
                        ab_moisture, & ! soil properties
                        grid_cell_area_ha, &
                        dead_fuel_fcl, live_grass_fuel_fcl, live_tree_fuel_fcl, & ! fuel
                        t_air_min, t_air_max, precip_daily_mmd, wind_sf, & ! weather
                        veg_cover_frac, &
                        nesterov_index, & ! dryness
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
    DEALLOCATE(tile_area)
#endif
  END SUBROUTINE update_q_dfire_spitfire


  ! ======================================================================================================= !
  !>
  !> Run science
  !>
  SUBROUTINE calc_spitfire(nc, nalphas, nfuelclasses, dtime, one_day, & ! in
                  lctlib_fuel_bulk_density, lctlib_fwd_wind_factor, &
                  density_human_population, a_nd, lightning_frequency, & ! ignition events
                  const_calorific_heat_content, const_i_surface_thd, &
                  alphas, &
                  moisture_extinction, &
                  ab_moisture, & ! soil properties
                  grid_cell_area, & 
                  dead_fuel, live_grass_fuel, live_tree_fuel, & ! fuels
                  t_air_min_c, t_air_max_c, precip_daily_mmd, wind_sf, & ! weather
                  veg_cover_frac, & ! in
                  nesterov_index, & ! inout
                  fuel_moisture, & ! out
                  spitfire, & ! out
                  mask)
    INTEGER, INTENT(in)   :: nc, nalphas, nfuelclasses
    REAL(wp), INTENT(in)  :: dtime, one_day

    REAL(wp), INTENT(in)  :: lctlib_fuel_bulk_density(nc)          !< Fuel bulk density [kg m-3]
    REAL(wp), INTENT(in)  :: lctlib_fwd_wind_factor(nc)            !< Forward wind factor [0 to 1]

    REAL(wp), INTENT(in)  :: density_human_population(nc)          !< Average human population density in the grid cell [individuals km-2]
    REAL(wp), INTENT(in)  :: a_nd(nc)                              !< Propensity of people to produce ignition events [ignitions individual-1 d-1]
    REAL(wp), INTENT(in)  :: lightning_frequency(nc)               !< Lightning frequency [ha-1 d-1]

    REAL(wp), INTENT(in)  :: const_calorific_heat_content          !< Fuel calorific heat content [kJ kg-1]
    REAL(wp), INTENT(in)  :: const_i_surface_thd                   !< Fire intensity threshold [kW m-1] 

    REAL(wp), INTENT(in)  :: alphas(nalphas)                       !< Drying parameter for 1, 10 and 100 hours fuel classes
    REAL(wp), INTENT(in)  :: moisture_extinction(:)                !< Moisture of extinction threshold [0 to 1]
    REAL(wp), INTENT(in)  :: ab_moisture(nc)                       !< above ground relative soil moisture [0 to 1]

    REAL(wp), INTENT(in)  :: grid_cell_area(nc)                    !< Gridcell area [ha]

    REAL(wp), INTENT(in)  :: dead_fuel(:,:,:)                      !< Dead fuel (#nc,#dead fuel pools, #fuel classes) [g DM m-2]
    REAL(wp), INTENT(in)  :: live_grass_fuel(:,:,:)                !< Live grass fuel (#nc,#live fuel pools, #fuel classes) [g DM m-2]
    REAL(wp), INTENT(in)  :: live_tree_fuel(:,:,:)                 !< Live tree fuel (#nc,#live fuel pools, #fuel classes) [g DM m-2]

    REAL(wp), INTENT(in)  :: t_air_min_c(nc), t_air_max_c(nc)      !< [Celsius], [Celsius], todo: refactor to t_weather derive type
    REAL(wp), INTENT(in)  :: precip_daily_mmd(nc)                  !< [mm d-1]
    REAL(wp), INTENT(in)  :: wind_sf(nc)                           !< [m min-1]

    TYPE(t_sf_veget_cover_frac), INTENT(in) :: veg_cover_frac      !< Vegetation cover fraction [0 to 1]

    REAL(wp), INTENT(INOUT) :: nesterov_index(nc)                  !< Dryness index [Celsius^2]

    REAL(wp), INTENT(out) :: fuel_moisture(nc,df_sf_omegas_nr)     !< Fuel moisture (#nc, #omega_nr (omega_o, omega_nl))

    TYPE(t_sf_spitfire), INTENT(out)    :: spitfire                !< Spitfire outputs
    LOGICAL, INTENT(in), OPTIONAL       :: mask(nc)

    LOGICAL :: mask_tile_enabled(nc)

    ! locals
    TYPE(t_q_sf_ignitions)  :: sf_ignitions
    REAL(wp) :: w_oi(nc, df_sf_fc_full_nr) !< Dead fuel by fuel classes [g DM m-2]
    REAL(wp) :: tot_lvg_fuel(nc)           !< Total livegrass fuel [g DM m-2]
    INTEGER :: ndeadpools
    TYPE(t_sf_fuel_class) :: fuel_cls
    LOGICAL  :: has_dead_fuel(nc)          !< True = dead fuel
    LOGICAL  :: has_fire_conditions(nc)    !< 1, 10, 100 and 1000 hours

    INTEGER                               :: i, j

    ! todo: deal with tile_frac = 0?

    ndeadpools = SIZE(dead_fuel, DIM=2) ! 2=dead fuel pools !fuel_cls % Count_dead_fuel_pools()
    mask_tile_enabled = .TRUE.
    IF (PRESENT(mask)) THEN
      mask_tile_enabled = mask
    ENDIF

    !< According to Thonicke et al. (2010), which is based on Rothermel et al.,
    !< fire spread depends on dead fuel availability. If there is no dead fuel,
    !< then fire cannot spread and we stop the calculation. The same applies if FDI is low.
    has_dead_fuel = sf_calc_has_any_fuel(dead_fuel(:,:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx))
    !< It accounts for 1, 10 and 100 hours fuel classes
    has_fire_conditions = mask_tile_enabled .AND. has_dead_fuel

    sf_ignitions = t_q_sf_ignitions(nc)
    spitfire = t_sf_spitfire(nc, ndeadpools, sf_fuel_live_fuel_nr)

    !> Calculate whether there are conditions for fire based on the Nesterov Index,
    !> following Thonicke et al. (2010) with modifications from Pfeiffer et al. (2013)
    !> and own changes for QUINCY.
    nesterov_index = calc_sf_nesterov_index(nc, t_air_max_c, t_air_min_c, &
        precip_daily_mmd, nesterov_index, &
        o_tsfactor = (dtime / one_day), &
        mask= has_fire_conditions)

    !< Calculate amount of dead fuel available for combustion following Thonicke et al. (2010)
    w_oi = SUM(dead_fuel, DIM=2) !< 2= dead fuel pools index, fuel class only

    DO CONCURRENT (i=1:SIZE(live_grass_fuel,2), j=1:SIZE(live_grass_fuel,3)) !2=fuel pools dim, 3=fuel classes dim
      tot_lvg_fuel = SUM(live_grass_fuel(:,i,j)) ! sum multiple dimensions
    END DO
    !< Compute fuel moisture (omega_o, omega_nl)
    fuel_moisture = sf_calc_fuel_moisture(nc, nalphas, df_sf_fc_spread_nr, alphas, &
                                          nesterov_index, ab_moisture, &
                                          w_oi(:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx), &
                                          tot_lvg_fuel)

    !> Fire danger is defined as zero in the absence of fuel or for wet fuel, and unity for completely dry fuel.
    !> Eq. 8 from Thonicke et al. (2010)
    WHERE(moisture_extinction > 0._wp)
      spitfire % fdi = max(0._wp, 1._wp  - (fuel_moisture(:, df_sf_omega_o_idx) / moisture_extinction))
    ENDWHERE

    !< It accounts for 1, 10 and 100 hours fuel classes
    has_fire_conditions = has_fire_conditions .AND. (spitfire % fdi > 1.e-6_wp)

    !> Fire rate of spread (IR, optimum reaction velocity, fuel bulk density, net fuel loading, ...)
    !> calculate the rate of spread following Thonicke et al. (2010)
    spitfire % ros = sf_calc_rate_of_spread(nc, df_sf_fc_spread_nr, & !< 3 = fuel classes for 1, 10 and 100 hours
        lctlib_fuel_bulk_density, lctlib_fwd_wind_factor, &
        wind_sf, moisture_extinction, spitfire % fdi, veg_cover_frac, &
        w_oi(:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx), fuel_moisture(:,df_sf_omega_o_idx), &
        mask=has_fire_conditions)

    !< Calculate fire ignition probability based on population and lightning, following Thonicke et al. (2010)
    sf_ignitions = calc_q_sf_ignition(nc, density_human_population, &
                                    a_nd, lightning_frequency, &
                                    df_sf_max_density_human_pop)
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

    !< Fraction of consumption used to calculate dead fuel consumption (Thonicke et al. (2010),
    !< extracted from equations B1, B4, B5 and B6)
    spitfire % dead_frac_consump = sf_calc_frac_consump(nc, nfuelclasses, moisture_extinction, &
                                                      fuel_moisture, mask=has_fire_conditions)

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
      IF (mask(inc)) THEN
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
  FUNCTION calc_ab_moisture(ab_soil_depth_sl, ab_wtr_soil_sl, ab_wtr_soil_fc_sl, mask) RESULT (ab_moisture)
    REAL(wp), INTENT(in)                :: ab_soil_depth_sl(:,:)      !< soil depth [m] (#nc, #above ground soil layers)
    REAL(wp), INTENT(in)                :: ab_wtr_soil_sl(:,:)        !< soil water content [m] (#nc, #above ground soil layers)
    REAL(wp), INTENT(in)                :: ab_wtr_soil_fc_sl(:,:)     !< soil water field capacity [m] (#nc, #above ground soil layers)
    LOGICAL, INTENT(in)                 :: mask(:)                    !< True = gridcell enabled
    REAL(wp), ALLOCATABLE               :: ab_moisture(:)             !< Relative above ground soil moisture [-] (nc)

    !< Select the first 2 layers from the soil moisture and calculate soil moisture content as a fraction of field capacity
    INTEGER :: is, ic
    INTEGER :: nc, nsl
    REAL(wp), ALLOCATABLE :: ab_moisture_sl(:,:)

    nc = SIZE(ab_wtr_soil_fc_sl, DIM=1)
    nsl = SIZE(ab_wtr_soil_fc_sl, DIM=2) ! 2=soil layers index
    ALLOCATE(ab_moisture_sl, mold=ab_wtr_soil_sl)

    !< Relative soil moisture for each above ground soil layers
    ab_moisture_sl = 0._wp
    DO is=1, nsl
      DO ic=1, nc
        IF (mask(ic) .AND. ab_wtr_soil_fc_sl(ic,is)>0._wp) THEN
          ab_moisture_sl(ic,is) = ab_wtr_soil_sl(ic,is) / ab_wtr_soil_fc_sl(ic,is)
        END IF
      END DO
    END DO

    !< Average of above ground soil layers
    IF (.NOT. ALLOCATED(ab_moisture)) ALLOCATE(ab_moisture(nc))
    ab_moisture = 0._wp
    DO ic=1, nc
      IF (mask(ic)) THEN
        ab_moisture(ic) = stats_average_weight(nsl, ab_moisture_sl(ic,:), ab_soil_depth_sl(ic,:))
      END IF
    ENDDO
  END FUNCTION calc_ab_moisture

#endif
END MODULE mo_q_dfire_spitfire
