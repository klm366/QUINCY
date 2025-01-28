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

  USE mo_dfire_sf_ignition,     ONLY: t_sf_ignitions

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

    sf_spitfire % e_n_ignitions = 0._wp
    sf_spitfire % e_n_humans = 0._wp
    sf_spitfire % e_n_lightning = 0._wp

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
    USE mo_jsb_process_class,       ONLY: VEG_, DFIRE_, SPQ_, DFIRE_FUEL_
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
    !dsl4jsb_Use_memory(A2L_)
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
    !dsl4jsb_Def_memory(A2L_)
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

    dsl4jsb_Real2D_onChunk                :: pot_burn_area     ! Potential burnable area (forest and grass areas) [m2?]
    dsl4jsb_Real2D_onChunk                :: pm_fire           ! Tree probability of mortality by fire [0 to 1]

    ! DFIRE_FUEL_ 3D
    dsl4jsb_Real3D_onChunk                :: dead_fuel         ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_fuel         ! in [g DM m-2]

    dsl4jsb_Real2D_onChunk                :: forest_fract      ! [%]
    dsl4jsb_Real2D_onChunk                :: grass_fract       ! [%]

    dsl4jsb_Real2D_onChunk                :: fuel_bulk_density ! [kg m-3]

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
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_sl)                ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_fc_sl)             ! in

    dsl4jsb_Get_var2D_onChunk(DFIRE_, t_air_min)              ! in [k]
    dsl4jsb_Get_var2D_onChunk(DFIRE_, t_air_max)              ! in

    dsl4jsb_Get_var2D_onChunk(DFIRE_, pot_burn_area)          ! out [m2?]
    dsl4jsb_Get_var2D_onChunk(DFIRE_, pm_fire)                ! out [%]

    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel)         ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_fuel)         ! in [g DM m-2]

    dsl4jsb_Get_var2D_onChunk(DFIRE_FUEL_, forest_fract)      ! in [%]
    dsl4jsb_Get_var2D_onChunk(DFIRE_FUEL_, grass_fract)       ! in [%]

    dsl4jsb_Get_var2D_onChunk(DFIRE_FUEL_, fuel_bulk_density) ! in [kg m-3]

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

    ! todo: veg_pool -> change units from molsC/m2 to gDM/m2
    ! todo: sb_pool -> change units from molsC/m3 to gDM/m2

    write(*,*) "update_spitfire:: fract_fpc shape?", SHAPE(fract_fpc)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !
    grid_cell_area => grid % area(ics:ice, iblk)
    tile_fract => tile%fract(ics:ice, iblk)

    tile_area = grid_cell_area * tile_fract
    pot_burn_area = grid_cell_area * (grass_fract + forest_fract)
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

    !< Select the first 2 layers from the soil moisture and calculate soil moisture content as a fraction of field capacity
    ab_w_soil_sl(:,:) = w_soil_sl(:, [df_sf_sl_indices])
    ab_w_soil_fc_sl(:,:) = w_soil_fc_sl(:, [df_sf_sl_indices])
    WHERE (ab_w_soil_fc_sl > 0._wp) ! should it happen?
      ab_moisture_sl(:,:) = ab_w_soil_sl(:,:) / ab_w_soil_fc_sl(:,:)
    END WHERE

    !< moisture average of above ground soil layers
    DO ic=1, nc
      IF (SUM(ab_moisture_sl(ic,:)) > 0._wp) THEN
        ab_moisture(ic) = sum(ab_moisture_sl(ic,:)) / size(ab_moisture_sl, dim=2) !< 2 = soil layers index
      END IF
    ENDDO

    CALL calc_spitfire(nc, dtime, one_day, spitfire)

    !< Test: assume 5% live fuel consumption
    live_fuel_consumption_1hr = live_fuel * 0.045 * 0.05 ! 4.5% 1 hour live fuel, 5%
    !< Test: assume 3% dead fuel consumption
    dead_fuel_consumption_1hr = dead_fuel * 0.045 * 0.03 ! 4.5% 1 hour live fuel, 3%

    area_burned = tile_area * 0.05        !< AB = 5%
    area_burned_fract = tile_fract * 0.05

    write(*,*) "(spitfire) I am burning..., iblk, dead_fuel_X=", iblk, dead_fuel(1,t_sf_fuel_dead_fuel_woody_ix)
    WRITE(*,*) "(spitfire) t_air_min, t_air_max, forest_fract, grass_fract, area, fbd = ", &
            t_air_min(1), t_air_max(1), forest_fract(1), grass_fract(1), pot_burn_area(1), &
            fuel_bulk_density(1)
  END SUBROUTINE update_dfire_spitfire


  ! ======================================================================================================= !
  !>
  !> Run science
  !>
  SUBROUTINE calc_spitfire(nc, dtime, one_day, spitfire)
    INTEGER, INTENT(in)   :: nc
    REAL(wp), INTENT(in)  :: dtime, one_day
    TYPE(t_sf_spitfire)   :: spitfire !< Spitfire outputs

    TYPE(t_sf_ignitions)  :: sf_ignitions

    sf_ignitions = t_sf_ignitions(nc)

    !< Calculate fire ignition probability based on population and lightning, following Thonicke et al. (2010)
    sf_ignitions = calc_sf_ignition(nc, density_human_population, &
                                    a_nd, lightning_frequency, &
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
