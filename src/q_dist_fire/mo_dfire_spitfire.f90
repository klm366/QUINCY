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
!>    - [-] Calc fuel load (mo_dfire_sf_fuel_class.f90)
!>    - [-] Calc Nesterov index (mo_dfire_sf_nesterov.f90)
!>    - [-] Calc fuel moisture (mo_dfire_sf_fmoist_class.f90)
!>    - [-] Calc Fire danger index (-)
!>    - [dead_fuel-all- > 0 & FDI>0] Calc Rate of spread (mo_dfire_sf_ros_class.f90)
!>    - [dead_fuel-all- > 0 & FDI>0] Calc Area Burned and its fraction (-)
!>    - [dead_fuel-1,10,100- > 0 & FDI>0] Calc dead fraction of consumption (mo_dfire_sf_frac_consump_class.f90)
!>    - [dead_fuel-1,10,100- > 0 & FDI>0]  Calc dead fuel consumption (mo_dfire_sf_fuel_consumption.f90)
!>    - [dead_fuel-1,10,100- > 0 & FDI>0]  Calc fire intensity (-)
!>    - [... & I_surface>50] Calc tree mortality (mo_dfire_sf_damage.f90)
!>    - [... & I_surface>50]  Calc live fuel consumption (mo_dfire_sf_fuel_consumption.f90)
!>    - [... & I_surface>50]  Calc emissions (mo_dfire_sf_emissions_class.f90)
!>    - [... & I_surface>50]  Derive fuel consumption for N,P,C13,C14,N15
!>    - [... & I_surface>50]  Calc residue from emissions
!>
!> mo_dfire_spitfire.f90 subroutines:
!>    -sf_calc_fire_intensity: calculate fire intensity
!>    -calc_area_burned_and_fract: calculate area burned and its fraction
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_dfire_spitfire
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_lnd_bgcm_store,        ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,  ONLY: SB_BGCM_LOSS_ID, SB_BGCM_POOL_ID
  USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract

  USE mo_dfire_constants,       ONLY: df_sf_sl_dims_nr, df_sf_sl_indices !< Soil related constants: top soil moisture layers 

  USE mo_dfire_constants,       ONLY: t_sf_fuel_dead_fuel_woody_ix
  USE mo_dfire_constants,       ONLY: df_sf_max_density_human_pop
  USE mo_dfire_constants,       ONLY: df_sf_fc_spread_nr
  USE mo_dfire_constants,       ONLY: df_sf_alpha_nr, df_sf_alphas !< Drying parameters

  USE mo_dfire_constants,       ONLY: sf_fuel_dead_fuel_nr, sf_fuel_live_fuel_nr, df_sf_fc_full_nr, &
                                      df_sf_fc_1hr_idx, df_sf_fc_10hr_idx, df_sf_fc_100hr_idx, df_sf_fc_1000hr_idx !< fuel classes indices

  USE mo_dfire_constants,       ONLY: df_sf_omegas_nr, & !< Fuel moisture dimension size (for omega_o and omega_nl)
                                      df_sf_omega_o_idx !< omega_o index from sf_calc_fuel_moisture output (array)

  USE mo_dfire_sf_fuel_class,   ONLY: t_sf_fuel_class, sf_calc_has_any_fuel
  USE mo_dfire_sf_fmoist_class, ONLY: sf_calc_fuel_moisture
  USE mo_dfire_sf_ignition,     ONLY: calc_sf_ignition, t_sf_ignitions
  USE mo_dfire_sf_nesterov,     ONLY: calc_sf_nesterov_index
  USE mo_dfire_sf_ros_class,    ONLY: t_sf_ros, t_sf_veget_cover_frac, sf_calc_rate_of_spread

  USE mo_dfire_sf_adaptor,      ONLY: sf_convert_units

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_dfire_spitfire

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_spitfire'

  !<
  !< Calculated Spitfire values
  !< Move spitfire data from calc_spitfire <-> update_dfire_spitfire subroutines
  !<
  TYPE t_sf_spitfire
    REAL(wp), ALLOCATABLE :: e_n_ignitions(:)
    REAL(wp), ALLOCATABLE :: e_n_humans(:)
    REAL(wp), ALLOCATABLE :: e_n_lightning(:)

    REAL(wp), ALLOCATABLE :: fdi(:) !< fire danger index [0 to 1]
    TYPE(t_sf_ros)        :: ros    !< fire dynamics
  END TYPE t_sf_spitfire

  INTERFACE t_sf_spitfire
    MODULE PROCEDURE :: create_sf_spitfire
  END INTERFACE t_sf_spitfire


  TYPE t_sf_bgcm_2L3D_p
    REAL(wp), POINTER :: p(:,:,:,:) !< dims: (part,elem,nc,layer)
  END TYPE t_sf_bgcm_2L3D_p

CONTAINS

  !>
  !> t_sf_spitfire constructor: Allocate and initialize
  !>
  FUNCTION create_sf_spitfire(nc) RESULT (sf_spitfire)
    INTEGER, INTENT(in) :: nc          !< #Gridcells
    TYPE(t_sf_spitfire) :: sf_spitfire !< Output

    ALLOCATE(sf_spitfire % e_n_ignitions(nc))
    ALLOCATE(sf_spitfire % e_n_humans(nc))
    ALLOCATE(sf_spitfire % e_n_lightning(nc))
    ALLOCATE(sf_spitfire % fdi(nc))

    sf_spitfire % e_n_ignitions = 0._wp
    sf_spitfire % e_n_humans = 0._wp
    sf_spitfire % e_n_lightning = 0._wp
    sf_spitfire % fdi = 0._wp
    sf_spitfire % ros = t_sf_ros(nc)

  END FUNCTION create_sf_spitfire

  ! ======================================================================================================= !
  !>
  !> Compute spitfire for current time step
  !> It links Quincy and Spitfire
  !>
  SUBROUTINE update_dfire_spitfire(tile, options)
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_process_class,       ONLY: VEG_, DFIRE_, SPQ_, DFIRE_FUEL_, A2L_
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_grid,                ONLY: Get_vgrid, Get_grid
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid, t_jsb_grid
    USE mo_jsb_math_constants,      ONLY: one_day

    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(DFIRE_)
    dsl4jsb_Use_memory(DFIRE_FUEL_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(VEG_)
    !dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_grid),       POINTER       :: grid
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lnd_bgcm_store), POINTER       :: bgcm_store !< the bgcm store of this tile
    INTEGER                               :: iblk, ics, ice, nc, ic
    INTEGER                               :: no_children
    REAL(wp)                              :: dtime      !< timestep
    CHARACTER(len=*), PARAMETER           :: routine = TRIM(modname)//':update_fdire_sf_fuel'

    REAL(wp), POINTER                     :: fract_fpc(:,:) ! in
    REAL(wp), POINTER                     :: grid_cell_area(:) ! [m2?]
    REAL(wp), POINTER                     :: tile_fract(:) ! [%]
    REAL(wp), ALLOCATABLE                 :: tile_area(:) ! [m2]

    REAL(wp)                              :: ab_w_soil_sl(options%nc,df_sf_sl_dims_nr)    !< above ground soil water [m]
    REAL(wp)                              :: ab_w_soil_fc_sl(options%nc,df_sf_sl_dims_nr) !< above ground soil field capacity [m]
    REAL(wp)                              :: ab_moisture_sl(options%nc,df_sf_sl_dims_nr)  !< above ground relative soil moisture layers [-]
    REAL(wp)                              :: ab_moisture(options%nc)                      !< top relative soil moisture layer [-]

    TYPE(t_sf_veget_cover_frac)           :: veg_cover_frac                               !< Vegetation cover

    !< dead fuel (#nc, #dead fuel pools, #fuel classes) [g DM m-2]
    REAL(wp)                              :: dead_fuel_fcl(options%nc,sf_fuel_dead_fuel_nr,df_sf_fc_full_nr)
    !< live grass fuel (#nc, #live fuel pools, #fuel classes) [g DM m-2]
    REAL(wp)                              :: live_grass_fuel_fcl(options%nc,sf_fuel_live_fuel_nr,df_sf_fc_full_nr)
    !< live tree fuel (#nc, #live fuel pools, #fuel classes) [g DM m-2]
    REAL(wp)                              :: live_tree_fuel_fcl(options%nc,sf_fuel_live_fuel_nr,df_sf_fc_full_nr)

    ! todo: move to memory
    REAL(wp)                              :: fuel_moisture(options%nc,df_sf_omegas_nr) ! fuel moisture (#nc, #omegas) [-]
    REAL(wp), ALLOCATABLE, SAVE           :: nesterov_index(:) !< (#nc) [C2]

    REAL(wp)                              :: wind_sf(options%nc)                          !< Wind into spitfire units [m min-1]
    REAL(wp)                              :: t_air_min_c(options%nc)                      !< Minimum daily temperature [Celsius]
    REAL(wp)                              :: t_air_max_c(options%nc)                      !< Maximum daily temperature [Celsius]

    TYPE(t_sf_spitfire)                   :: sf_spitfire

    INTEGER                               :: nr_of_children
    TYPE(t_sf_bgcm_2L3D_p), ALLOCATABLE   :: sb_loss_stores(:) !< Reference to each sb_loss pft pools (#npfts)
    ! ----------------------------------------------------------------------------------------------------- !
    !dsl4jsb_Def_mt2L2D :: veg_pool_mt
    !dsl4jsb_Def_mt2L3D :: sb_loss_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(DFIRE_)
    dsl4jsb_Def_memory(DFIRE_FUEL_)
    dsl4jsb_Def_memory(VEG_)
    !dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Local vars here
    ! ----------------------------------------------------------------------------------------------------- !
    !dsl4jsb_Real2D_onChunk                :: e_n_ignitions    ! out, Total rate of ignitions events [ha-1 ts-1]
    dsl4jsb_Real3D_onChunk                :: w_soil_sl        ! soil water content [m]
    dsl4jsb_Real3D_onChunk                :: w_soil_fc_sl     ! soil water field capacity [m]

    ! DFUEL_ 2D
    ! Note: t_air_max/t_air_min belongs to the previous days due to implementation constraints
    dsl4jsb_Real2D_onChunk                :: t_air_min         ! Min temperature of the day (previous day) [k]
    dsl4jsb_Real2D_onChunk                :: t_air_max         ! Max temperature of the day (previous day) [k]
    dsl4jsb_Real2D_onChunk                :: precip_daily_mmd  ! out [mm day-1]
    dsl4jsb_Real2D_onChunk                :: precip_daily      ! in [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk                :: wind_10m          ! in A2L_ [m s-1]

    dsl4jsb_Real2D_onChunk                :: pot_burn_area     ! Potential burnable area (forest and grass areas) [m2?]
    dsl4jsb_Real2D_onChunk                :: pm_fire           ! Tree probability of mortality by fire [0 to 1]

    ! DFIRE_FUEL_ 3D
    dsl4jsb_Real3D_onChunk                :: dead_fuel                  ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_1hr_dfp          ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_10hr_dfp         ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_100hr_dfp        ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_1000hr_dfp       ! in [g DM m-2]

    dsl4jsb_Real3D_onChunk                :: live_fuel                  ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_grass_fuel_1hr_lfp    ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_grass_fuel_10hr_lfp   ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_grass_fuel_100hr_lfp  ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_grass_fuel_1000hr_lfp ! in [g DM m-2]

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
    dsl4jsb_Real3D_onChunk                :: live_fuel_consumption_1hr ! out [g DM m-2]

    dsl4jsb_Real2D_onChunk                :: area_burned               ! out [ha ts-1]
    dsl4jsb_Real2D_onChunk                :: area_burned_fract         ! out [-]
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(DFIRE_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(DFIRE_)
    dsl4jsb_Get_memory(DFIRE_)
    dsl4jsb_Get_memory(DFIRE_FUEL_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_sl)                ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_fc_sl)             ! in

    dsl4jsb_Get_var2D_onChunk(DFIRE_, t_air_min)              ! in [k]
    dsl4jsb_Get_var2D_onChunk(DFIRE_, t_air_max)              ! in
    dsl4jsb_Get_var2D_onChunk(DFIRE_, precip_daily_mmd)       ! out
    dsl4jsb_Get_var2D_onChunk(DFIRE_, precip_daily)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_10m)                 ! in

    dsl4jsb_Get_var2D_onChunk(DFIRE_, pot_burn_area)          ! out [m2?]
    dsl4jsb_Get_var2D_onChunk(DFIRE_, pm_fire)                ! out [%]

    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel)             ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_1hr_dfp)     ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_10hr_dfp)    ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_100hr_dfp)   ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_1000hr_dfp)  ! in [g DM m-2]

    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_fuel)                  ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_grass_fuel_1hr_lfp)    ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_grass_fuel_10hr_lfp)   ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_grass_fuel_100hr_lfp)  ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_grass_fuel_1000hr_lfp) ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_tree_fuel_1hr_lfp)     ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_tree_fuel_10hr_lfp)    ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_tree_fuel_100hr_lfp)   ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_tree_fuel_1000hr_lfp)  ! in [g DM m-2]

    dsl4jsb_Get_var2D_onChunk(DFIRE_FUEL_, forest_fract)      ! in [%]
    dsl4jsb_Get_var2D_onChunk(DFIRE_FUEL_, grass_fract)       ! in [%]

    dsl4jsb_Get_var2D_onChunk(DFIRE_FUEL_, fuel_bulk_density)   ! in [kg m-3]
    dsl4jsb_Get_var2D_onChunk(DFIRE_FUEL_, moisture_extinction) ! in [-]
    dsl4jsb_Get_var2D_onChunk(DFIRE_FUEL_, fwd_wind_factor)     ! in [-]

    dsl4jsb_Get_var3D_onChunk(DFIRE_, dead_fuel_consumption_1hr) ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_, live_fuel_consumption_1hr) ! out [g DM m-2]

    dsl4jsb_Get_var2D_onChunk(DFIRE_, area_burned)               ! out [ha ts-1]
    dsl4jsb_Get_var2D_onChunk(DFIRE_, area_burned_fract)         ! out [-]
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    !dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    !dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    ! ----------------------------------------------------------------------------------------------------- !
!    dsl4jsb_Get_mt2L3D(SB_BGCM_LOSS_ID, sb_loss_mt)

    no_children = tile%Get_no_of_children()

    !< Init variables
    fract_fpc => veg__mem%collect_var_real2d(ics, ice, iblk, no_children, veg__mem%fract_fpc)

    ! gather children's bgcm stores from specific flux pool
    sb_loss_stores = Get_childrens_bgcm_stores_2l_3d(SB_BGCM_LOSS_ID, ics, ice, iblk, tile)

    ab_moisture_sl = 0._wp
    ab_moisture = 0._wp
    pm_fire = 0._wp
    ! rebuild fuel variables
    dead_fuel_fcl(:,:,df_sf_fc_1hr_idx) = dead_fuel_1hr_dfp ! dead fuel
    dead_fuel_fcl(:,:,df_sf_fc_10hr_idx) = dead_fuel_1hr_dfp
    dead_fuel_fcl(:,:,df_sf_fc_100hr_idx) = dead_fuel_1hr_dfp
    dead_fuel_fcl(:,:,df_sf_fc_1000hr_idx) = dead_fuel_1hr_dfp
    live_grass_fuel_fcl(:,:,df_sf_fc_1hr_idx) = live_grass_fuel_1hr_lfp ! live grass fuel
    live_grass_fuel_fcl(:,:,df_sf_fc_10hr_idx) = live_grass_fuel_10hr_lfp
    live_grass_fuel_fcl(:,:,df_sf_fc_100hr_idx) = live_grass_fuel_100hr_lfp
    live_grass_fuel_fcl(:,:,df_sf_fc_1000hr_idx) = live_grass_fuel_1000hr_lfp
    live_tree_fuel_fcl(:,:,df_sf_fc_1hr_idx) = live_tree_fuel_1hr_lfp ! live tree fuel
    live_tree_fuel_fcl(:,:,df_sf_fc_10hr_idx) = live_tree_fuel_10hr_lfp
    live_tree_fuel_fcl(:,:,df_sf_fc_100hr_idx) = live_tree_fuel_100hr_lfp
    live_tree_fuel_fcl(:,:,df_sf_fc_1000hr_idx) = live_tree_fuel_1000hr_lfp

    write(*,*) "update_spitfire:: fract_fpc shape?", SHAPE(fract_fpc)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !
    grid_cell_area => grid % area(ics:ice, iblk)
    tile_fract => tile%fract(ics:ice, iblk)

    tile_area = grid_cell_area * tile_fract
    pot_burn_area = grid_cell_area * (grass_fract + forest_fract)
    WRITE(*,*) "(spitfire.f90) calc_spitfire:: tile_fract, tile_area, gc area", tile_fract(9), &
            tile_area(9), grid_cell_area(9)

    ! assertion
    IF (ANY(pot_burn_area > tile_area)) THEN
      DO ic=1, nc
        IF (pot_burn_area(ic) > tile_area(ic)) THEN
          WRITE(*,*) "pot_burn_area, tile_fract=", grid_cell_area(ic), tile_fract(ic), &
                                                   (grass_fract(ic) + forest_fract(ic))
          WRITE(*,*) "pot_burn_area, tile_fract=", pot_burn_area(ic), tile_area(ic)
        ENDIF
      ENDDO
      CALL finish(TRIM(routine), "potential_burn_area > tile_fract")
    END IF

    !> Units conversion
    t_air_max_c = sf_convert_units("Kelvin", "Celsius", t_air_max)
    t_air_min_c = sf_convert_units("Kelvin", "Celsius", t_air_min)
    precip_daily_mmd = sf_convert_units("kg m-2 s-1", "mm/day", precip_daily)
    wind_sf = sf_convert_units("m/s", "m/min", wind_10m) !< wind to spitfire units
    !DBH(:) = sf_convert_units("m", "cm", diameter) !< dbh to spitfire units

    !< Select the first 2 layers from the soil moisture and calculate soil moisture content as a fraction of field capacity
    ab_w_soil_sl(:,:) = w_soil_sl(:, [df_sf_sl_indices])
    ab_w_soil_fc_sl(:,:) = w_soil_fc_sl(:, [df_sf_sl_indices])
    WHERE (ab_w_soil_fc_sl > 0._wp) ! should it happen?
      ab_moisture_sl(:,:) = ab_w_soil_sl(:,:) / ab_w_soil_fc_sl(:,:)
    END WHERE

    !< Pack vegetation cover fraction
    veg_cover_frac = t_sf_veget_cover_frac(nc, tree_frac, grass_frac)

    !< moisture average of above ground soil layers
    DO ic=1, nc
      IF (SUM(ab_moisture_sl(ic,:)) > 0._wp) THEN
        ab_moisture(ic) = sum(ab_moisture_sl(ic,:)) / size(ab_moisture_sl, dim=2) !< 2 = soil layers index
      END IF
    ENDDO

    IF (.NOT. ALLOCATED(nesterov_index)) THEN
      ALLOCATE(nesterov_index(nc))
      nesterov_index = 0._wp
    ENDIF

    CALL calc_spitfire(nc, df_sf_alpha_nr, dtime, one_day, &
                        fuel_bulk_density, fwd_wind_factor, & ! lctlib weighted
                        dsl4jsb_Config(DFIRE_)%density_human_population, & ! ignition events
                        dsl4jsb_Config(DFIRE_)%a_nd, &
                        dsl4jsb_Config(DFIRE_)%lightning_frequency, &
                        df_sf_alphas, &
                        moisture_extinction, &
                        ab_moisture, & ! soil properties
                        dead_fuel_fcl, live_grass_fuel_fcl, live_tree_fuel_fcl, & ! fuel
                        t_air_min, t_air_max, precip_daily_mmd, wind_sf, & ! weather
                        nesterov_index, & ! dryness
                        fuel_moisture, & ! fuel properties
                        sf_spitfire)

    !< Test: assume 5% live fuel consumption
    live_fuel_consumption_1hr = live_fuel * 0.045 * 0.05 ! 4.5% 1 hour live fuel, 5%
    !< Test: assume 3% dead fuel consumption
    dead_fuel_consumption_1hr = dead_fuel * 0.045 * 0.03 ! 4.5% 1 hour live fuel, 3%

    area_burned = tile_area * 0.05        !< AB = 5%
    area_burned_fract = tile_fract * 0.05

    write(*,*) "(spitfire) I am burning..., iblk, dead_fuel_X=", iblk, dead_fuel(1,t_sf_fuel_dead_fuel_woody_ix)
    !WRITE(*,*) "(spitfire) t_air_min, t_air_max, forest_fract, grass_fract, area, fbd = ", &
    !        t_air_min(1), t_air_max(1), forest_fract(1), grass_fract(1), pot_burn_area(1), &
    !        fuel_bulk_density(1)
  END SUBROUTINE update_dfire_spitfire


  ! ======================================================================================================= !
  !>
  !> Run science
  !>
  SUBROUTINE calc_spitfire(nc, nalphas, dtime, one_day, & ! in
                  lctlib_fuel_bulk_density, lctlib_fwd_wind_factor, &
                  config_density_human_population, config_a_nd, config_lightning_frequency, & ! ignition events
                  alphas, &
                  moisture_extinction, &
                  ab_moisture, & ! soil properties
                  dead_fuel, live_grass_fuel, live_tree_fuel, & ! fuels
                  t_air_min_c, t_air_max_c, precip_daily_mmd, wind_sf, & ! weather
                  veg_cover_fract, & in
                  nesterov_index, & ! inout
                  fuel_moisture, & ! out
                  spitfire) ! out
    INTEGER, INTENT(in)   :: nc, nalphas
    REAL(wp), INTENT(in)  :: dtime, one_day

    REAL(wp), INTENT(in)  :: lctlib_fuel_bulk_density(nc)          !< Fuel bulk density [kg m-3]
    REAL(wp), INTENT(in)  :: lctlib_fwd_wind_factor(nc)            !< Forward wind factor [0 to 1]

    REAL(wp), INTENT(in)  :: config_density_human_population       !< Average human population density in the grid cell [individuals km-2]
    REAL(wp), INTENT(in)  :: config_a_nd                           !< Propensity of people to produce ignition events [ignitions individual-1 d-1]
    REAL(wp), INTENT(in)  :: config_lightning_frequency            !< Lightning frequency [ha-1 d-1]

    REAL(wp), INTENT(in)  :: alphas(nalphas)                       !< Drying parameter for 1, 10 and 100 hours fuel classes
    REAL(wp), INTENT(in)  :: moisture_extinction(:)                !< Moisture of extinction threshold [0 to 1]
    REAL(wp), INTENT(in)  :: ab_moisture(nc)                       !< above ground relative soil moisture [0 to 1]

    REAL(wp), INTENT(in)  :: dead_fuel(:,:,:)                      !< Dead fuel (#nc,#dead fuel pools, #fuel classes) [g DM m-2]
    REAL(wp), INTENT(in)  :: live_grass_fuel(:,:,:)                !< Live grass fuel (#nc,#live fuel pools, #fuel classes) [g DM m-2]
    REAL(wp), INTENT(in)  :: live_tree_fuel(:,:,:)                 !< Live tree fuel (#nc,#live fuel pools, #fuel classes) [g DM m-2]

    REAL(wp), INTENT(in)  :: t_air_min_c(nc), t_air_max_c(nc)      !< [Celsius], [Celsius], todo: refactor to t_weather derive type
    REAL(wp), INTENT(in)  :: precip_daily_mmd(nc)                  !< [mm d-1]
    REAL(wp), INTENT(in)  :: wind_sf(nc)                           !< [m min-1]

    TYPE(t_sf_veget_cover_frac), INTENT(in) :: veg_cover_frac      !< Vegetation cover fraction [0 to 1]

    REAL(wp), INTENT(INOUT) :: nesterov_index(nc)                  !< Dryness index [Celsius^2]

    REAL(wp), INTENT(out) :: fuel_moisture(nc,df_sf_omegas_nr)     !< Fuel moisture (#nc, #omega_nr (omega_o, omega_nl))

    TYPE(t_sf_spitfire), INTENT(out) :: spitfire                   !< Spitfire outputs

    ! locals
    TYPE(t_sf_ignitions)  :: sf_ignitions
    REAL(wp) :: w_oi(nc, df_sf_fc_full_nr) !< Dead fuel by fuel classes [g DM m-2]
    REAL(wp) :: tot_lvg_fuel(nc)           !< Total livegrass fuel [g DM m-2]
    INTEGER :: ndeadpools
    TYPE(t_sf_fuel_class) :: fuel_cls
    LOGICAL  :: has_dead_fuel(nc)          !< True = dead fuel
    LOGICAL  :: has_fire_conditions(nc)    !< 1, 10, 100 and 1000 hours

    INTEGER                               :: i, j

    ! todo: deal with tile_frac = 0?

    sf_ignitions = t_sf_ignitions(nc)
    spitfire = t_sf_spitfire(nc)

    !> Initialize object by computing dead and live fuels
    !fuel_cls = t_sf_fuel_class(nc, lctlib_growthform, &
    !    dead_fuel_frac, live_fuel_frac, &
    !    dead_fuel_mt, live_fuel_mt, &
    !    ratio_polymeric_litter)
    ndeadpools = SIZE(dead_fuel, DIM=2) ! 2=dead fuel pools !fuel_cls % Count_dead_fuel_pools()

    !> Calculate whether there are conditions for fire based on the Nesterov Index,
    !> following Thonicke et al. (2010) with modifications from Pfeiffer et al. (2013)
    !> and own changes for QUINCY.
    nesterov_index = calc_sf_nesterov_index(nc, t_air_max_c, t_air_min_c, &
        precip_daily_mmd, nesterov_index, &
        o_tsfactor = (dtime / one_day) )

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
    !WRITE(*,*) "(spitfire.f90) calc_spitfire:: nc, me == 0?", nc, COUNT(moisture_extinction == 0._wp)
    !WRITE(*,*) "(spitfire.f90) calc_spitfire:: idx==0?", pack([(ix,ix=1,size(moisture_extinction))],moisture_extinction==0._wp)
    WHERE(moisture_extinction > 0._wp)
      spitfire % fdi = max(0._wp, 1._wp  - (fuel_moisture(:, df_sf_omega_o_idx) / moisture_extinction))
    ENDWHERE

    !< According to Thonicke et al. (2010), which is based on Rothermel et al.,
    !< fire spread depends on dead fuel availability. If there is no dead fuel,
    !< then fire cannot spread and we stop the calculation. The same applies if FDI is low.
    has_dead_fuel = sf_calc_has_any_fuel(dead_fuel(:,:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx))
    !< It accounts for 1, 10 and 100 hours fuel classes
    has_fire_conditions(:) = ANY(has_dead_fuel) .AND. ANY(spitfire % fdi > 1.e-6_wp)

    !< If there are no fire conditions, stop calculation and return burned area and fire emissions as zero.
    !< Previously initialized to zero.
    !< Such conditions accounts for multiple grid cells running at the same time (IQ)
!    IF (ALL(has_fire_conditions .EQV. .FALSE.)) RETURN

    !> Fire rate of spread (IR, optimum reaction velocity, fuel bulk density, net fuel loading, ...)
    !> calculate the rate of spread following Thonicke et al. (2010)
    spitfire % ros = sf_calc_rate_of_spread(nc, df_sf_fc_spread_nr, & !< 3 = fuel classes for 1, 10 and 100 hours
        lctlib_fuel_bulk_density, lctlib_fwd_wind_factor, &
        wind_sf, moisture_extinction, spitfire % fdi, veg_cover_frac, &
        w_oi(:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx), fuel_moisture(:,df_sf_omega_o_idx), &
        mask=has_fire_conditions)

    !< Calculate fire ignition probability based on population and lightning, following Thonicke et al. (2010)
    sf_ignitions = calc_sf_ignition(nc, config_density_human_population, &
                                    config_a_nd, config_lightning_frequency, &
                                    df_sf_max_density_human_pop)
    spitfire%e_n_ignitions = sf_ignitions%total * (dtime / one_day)
    spitfire%e_n_humans = sf_ignitions%human * (dtime / one_day)
    spitfire%e_n_lightning = sf_ignitions%lightning * (dtime / one_day)
  END SUBROUTINE calc_spitfire


  ! ======================================================================================================= !
  !>
  !> Gather children's bgcm stores from specific flux pool
  !>
  FUNCTION Get_childrens_bgcm_stores_2l_3d(pool_id, ics, ice, iblk, tile) RESULT (bgcm_stores)
    INTEGER, INTENT(in) :: pool_id
    INTEGER, INTENT(in) :: ics, ice, iblk
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile         !< one tile with data structure for one lct
    TYPE(t_sf_bgcm_2L3D_p), ALLOCATABLE   :: bgcm_stores(:) ! output

    INTEGER :: i_tile, no_children
    CLASS(t_jsb_tile_abstract),  POINTER  :: current_tile

    no_children = tile%Get_no_of_children()
    ! gather children's bgcm stores from specific flux pool
    ALLOCATE(bgcm_stores(no_children)) !< move to spitfire_init, own function
    current_tile => tile%Get_first_child_tile()
    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile))
      i_tile = i_tile + 1
      !write(*,*) "update_spitfire:: name tile=", current_tile%name
      bgcm_stores(i_tile)%p => current_tile%bgcm_store%Get_matrix_from_store_2l_3d(pool_id, ics, ice, iblk, current_tile%name)
      current_tile => current_tile%Get_next_sibling_tile()
    ENDDO
  END FUNCTION Get_childrens_bgcm_stores_2l_3d

#endif
END MODULE mo_dfire_spitfire
