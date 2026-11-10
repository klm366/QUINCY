!> Contains the interfaces to the carbon process
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

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_carbon_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on, acc_stream
  USE mo_kind,            ONLY: lp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  USE mo_carbon_constants,   ONLY: molarMassCO2_kg, sec_per_day,           &
                                 & Y_ACID, Y_WATER, Y_ETHANOL,             &
                                 & Y_NONSOLUBLE, Y_HUMUS, WOODLITSIZE,     &
                                 & ROOT_DIST_BETA, permafrost_threshold,   &
                                 & Y_ACID, Y_WATER, Y_ETHANOL, Y_NONSOLUBLE

  ! Use of processes in this module (Get_carbon_memory and Get_carbon_config)
  dsl4jsb_Use_processes CARBON_, ASSIMI_, PHENO_, CLIMBUF_, DISTURB_, NLCC_, L2A_, SSE_, HYDRO_

  ! Use of process configurations
  dsl4jsb_Use_config(PHENO_)
  dsl4jsb_Use_config(CARBON_)
  dsl4jsb_Use_config(HYDRO_)

  ! Use of process memories (t_carbon_memory)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(CLIMBUF_)
  dsl4jsb_Use_memory(L2A_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(HYDRO_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_carbon_tasks
  PUBLIC ::  recalc_carbon_per_tile_vars, recalc_per_tile_vars,                                            &
    &        calculate_current_c_ta_state_sum, check_carbon_conservation, yday_carbon_conservation_test,   &
    &        rescale_carbon_upon_reference_area_change,                                                    &
    &        carbon_transfer_from_active_to_passive_vars_onChunk, global_carbon_diagnostics,               &
    &        rescale_carbon_vars_upon_ageing_induced_ref_area_change, merge_carbon_vars_upon_ageing

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_interface'

  !> Type definition for C_NPP_pot_allocation task
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_C_NPP_pot_allocation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_C_NPP_pot_allocation    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_C_NPP_pot_allocation !< Aggregates computed task variables
  END TYPE   tsk_C_NPP_pot_allocation

  !> Constructor interface for C_NPP_pot_allocation task
  INTERFACE tsk_C_NPP_pot_allocation
    PROCEDURE Create_task_C_NPP_pot_allocation        !< Constructor function for task
  END INTERFACE tsk_C_NPP_pot_allocation

  !> Type definition for ta_vars_after_lcc task
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_ta_vars_after_lcc
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_ta_vars_after_lcc    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_ta_vars_after_lcc !< Aggregates computed task variables
  END TYPE   tsk_ta_vars_after_lcc

  !> Constructor interface for ta_vars_after_lcc task
  INTERFACE tsk_ta_vars_after_lcc
    PROCEDURE Create_task_ta_vars_after_lcc        !< Constructor function for task
  END INTERFACE tsk_ta_vars_after_lcc

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for carbon task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "carbon"
  !!
  FUNCTION Create_task_C_NPP_pot_allocation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_C_NPP_pot_allocation::return_ptr)
    CALL return_ptr%Construct(name='C_NPP_pot_allocation', process_id=CARBON_, owner_model_id=model_id)

  END FUNCTION Create_task_C_NPP_pot_allocation

  ! ====================================================================================================== !
  !
  !> Constructor for ta variables after lcc
  !
  FUNCTION Create_task_ta_vars_after_lcc(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_ta_vars_after_lcc::return_ptr)
    CALL return_ptr%Construct(name='ta_vars_after_lcc', process_id=CARBON_, owner_model_id=model_id)

  END FUNCTION Create_task_ta_vars_after_lcc

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for carbon process
  !!
  !! @param[in,out] this      Instance of carbon process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_carbon_tasks(this, model_id)

    ! in/out
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_C_NPP_pot_allocation(model_id))
    CALL this%Register_task(tsk_ta_vars_after_lcc(model_id))

  END SUBROUTINE Register_carbon_tasks

  ! ================================================================================================================================
  !>
  !> Implementation to calculate the allocation of NPP carbon to the plant and soil carbon pools
  !!        R: Corresponds to update_cbalance_bethy in JSBACH3.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_C_NPP_pot_allocation(tile, options)

    ! Use declarations
    USE mo_carbon_process,           ONLY: calc_plant_c, calc_soil_c
    USE mo_util,                     ONLY: soil_depth_to_layers
    USE mo_jsb_time,                 ONLY: is_newday, is_newyear, timesteps_per_day, &
      &                                    is_time_experiment_start, is_time_restart
    USE mo_carbon_constants,         ONLY: molarMassC_kg, C_2_biomass_ratio, rho_org_bg

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Declare local variables
    TYPE(t_jsb_model), POINTER   :: model
    TYPE(t_jsb_grid),  POINTER   :: grid
    TYPE(t_jsb_vgrid), POINTER   :: soil_c
    LOGICAL                      :: lstart, lresume, new_day, new_year
    LOGICAL                      :: nlcc_is_enabled    ! helper needed on DWD NEC
    INTEGER                      :: iblk, ics, ice, nc, ic, is, nsoil
    INTEGER                      :: tmp_timesteps_per_day
    INTEGER                      :: lctlib_PhenologyType
    REAL(lp)                     :: lctlib_fract_npp_2_woodPool, lctlib_fract_NPP_2_reservePool
    REAL(lp)                     :: lctlib_fract_NPP_2_exudates, lctlib_fract_green_2_herbivory
    REAL(lp)                     :: lctlib_tau_c_woods, lctlib_LAI_shed_constant, lctlib_MaxLAI
    REAL(lp)                     :: lctlib_Max_C_content_woods
    REAL(lp)                     :: lctlib_specificLeafArea_C, lctlib_reserveC2leafC
    CHARACTER(len=*), PARAMETER  :: routine = modname//':update_C_NPP_pot_allocation'
    REAL(lp)                     :: dtime
    REAL(lp)                     :: delta_time_C        !< Carbon time step [s]
    REAL(lp), DIMENSION(options%nc) :: MaxLai
    REAL(lp), DIMENSION(options%nc) :: old_c_state_sum_ta, current_fluxes
    REAL(lp), DIMENSION(options%nc) :: cflux_c_green_2_litter, cflux_c_wood_2_litter, cflux_crop_harvest_2_atmos

    ! Variables required for vertical YASSO
    REAL(lp), ALLOCATABLE  :: layer_depth(:)            !< soil layer depth [m]
    REAL(lp), ALLOCATABLE  :: layer_thick(:)            !< soil layer thickness [m]

    REAL(lp), ALLOCATABLE  :: soil_depth_thaw(:)        !< Minimum of either soil or thaw depth [m]
    REAL(lp), ALLOCATABLE  :: thaw_depth(:)             !< Maximum thaw depth / active layer depth [m]

    REAL(lp) :: mass_organic_soil !< mass of organic soil fraction in a soil layer [kg/m2]

    ! Pool and parameter transfer for YASSO scheme
    REAL(lp) :: LeafLit_coef(5)    !< Array with leaf litter coefficients []
    REAL(lp) :: WoodLit_coef(5)    !< Array with woody litter coefficients []
    REAL(lp) :: MiscPar(2)         !< Array with additional coefficients

    REAL(lp), POINTER            :: tile_fract(:) !< Fraction of the tile (rel. to the grid cell) []

    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_tile_abstract), POINTER :: climbuf   ! Pointer to the climate buffer memory residing on the veg tile
    dsl4jsb_Def_memory_tile(CLIMBUF_, climbuf)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk ::  temp_15d_yDay_climbuf            !< Smoothed air temperature [K]
    dsl4jsb_Real2D_onChunk ::  precip_15d_yDay_climbuf          !< Smoothed precipitation [kg m-2 s-1]
    dsl4jsb_Real3D_onChunk ::  t_soil_sl_15d_yDay_climbuf       !< Smoothed soil temperature [K]
    dsl4jsb_Real3D_onChunk ::  wtr_soil_rel_sl_15d_yDay_climbuf !< Smoothed relative soil water content [-]
    dsl4jsb_Real3D_onChunk ::  fract_org_sl                     !< Vol. soil layer organic matter fraction [-]

    dsl4jsb_Real2D_onChunk ::  cconservation_calcCpools ! C conservation test: Deviation from conservation

    dsl4jsb_Real2D_onChunk ::  cflux_c_greenwood_2_litter           ! Carbon flux from the veget. to the litter pools
                                                      ! [mol(C)/m^2(canopy) s]
    dsl4jsb_Real2D_onChunk ::  c_green                ! C-pool for leaves, fine roots, vegetative organs and
                                                      ! other green (living) parts of vegetation [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_reserve              ! C-pool for carbohydrate reserve (sugars, starches) that
                                                      ! allows plants to survive bad times[mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_woods                ! C-pool for stems, thick roots and other (dead) structural
                                                      !  material of living plants [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_crop_harvest         ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]

    !SIZE CLASS 1
    dsl4jsb_Real2D_onChunk ::  c_acid_ag1             ! Yasso above ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_ag1            ! Yasso above ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag1          ! Yasso above ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag1       ! Yasso above ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_acid_bg1             ! Yasso below ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_bg1            ! Yasso below ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg1          ! Yasso below ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg1       ! Yasso below ground litter-pool for non-soluble litter
                                                      ! [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_humus_1              ! Yasso below ground litter-pool for slow C compartment
                                                      !  [mol(C)/m^2(canopy)]
    !SIZE CLASS 2
    dsl4jsb_Real2D_onChunk ::  c_acid_ag2             ! Yasso above ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_ag2            ! Yasso above ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag2          ! Yasso above ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag2       ! Yasso above ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_acid_bg2             ! Yasso below ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_bg2            ! Yasso below ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg2          ! Yasso below ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg2       ! Yasso below ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_humus_2              ! Yasso below ground litter-pool for slow C compartment
                                                      !  [mol(C)/m^2(canopy)]

    dsl4jsb_Real2D_onChunk ::  c_decomp_humus_1_sum   ! Annual sum of humus decomposed in YASSO (leaf)
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_decomp_humus_2_sum   ! Annual sum of humus decomposed in YASSO (wood)
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_into_humus_1_sum     ! Annual sum of cflux into humus in YASSO (leaf)
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_into_humus_2_sum     ! Annual sum of cflux into humus in YASSO (wood)
                                                      !  [mol(C)/m^2(canopy)]

    dsl4jsb_Real2D_onChunk ::  gross_assimilation_ca  ! Gross assimilation of all canopy layers together
    dsl4jsb_Real2D_onChunk ::  NPP_pot_rate_ca        ! The instantaneous (potential) NPP rate [mol(C)/m^2(canopy) s]
    !dsl4jsb_Real2D_onChunk ::  NPP_pot_Rate_acc      ! averaged NPP rate over one day [mol(C)/m^2(canopy) s]
    dsl4jsb_Real2D_onChunk ::  LAI_sum                ! used to accumulate LAI over a day. Sum of LAI since Midnight.
    dsl4jsb_Real2D_onChunk ::  NPP_pot_sum            ! used to accumulated NPP-Rate over a day
    dsl4jsb_Real2D_onChunk ::  GPP_sum                ! used to accumulated GPP-Rate over a day
    dsl4jsb_Real2D_onChunk ::  LAI_yDayMean           ! previous days mean
    dsl4jsb_Real2D_onChunk ::  LAI_yyDayMean          ! previous previous days mean
    dsl4jsb_Real2D_onChunk ::  NPP_pot_yDayMean       ! mean value of NPP-Rate yesterday (from NPP_pot_sum())
                                                      ! [mol(CO2)/(m^2(canopy) s)]
                                                      ! = cbalance%NPP_pot_sum(...)/time_steps_per_day
    dsl4jsb_Real2D_onChunk ::  NPP_act_yDayMean       ! mean value of actual NPP-Rate yesterday after excess carbon is removed.
                                                      ! [mol(CO2)/(m^2(canopy) s)]
    dsl4jsb_Real2D_onChunk ::  GPP_yDayMean

    dsl4jsb_Real2D_onChunk ::  soil_respiration       ! mean daily rate of heterotrophic (soil) respiration
                                                      ! [mol(CO2)/m^2(ground)].
    dsl4jsb_Real2D_onChunk ::  NPP_flux_correction    ! Daily updated flux correction from yesterdays carbon balance
                                                      ! [mol(CO2)/m^2(canopy) s]
                                                      ! Amount by which the NPP rate entering the routine has to be corrected.
                                                      ! This correction arises either because otherwise the reserve pool would
                                                      ! get negative (positive correction), or the wood pool would exceed its
                                                      ! maximum value (negative correction).
    dsl4jsb_Real2D_onChunk ::  NPP_excess             ! That part of NPP that because of structural limits could not be
                                                      ! allocated in carbon pools [mol(CO2)/m^2(canopy) s]

    dsl4jsb_Real2D_onChunk ::  root_exudates          ! Total root exudates entering to the litter green pools
                                                      ! [mol(C)/m^2(canopy) s]

    dsl4jsb_Real2D_onChunk ::  c_sum_veg_ta
    dsl4jsb_Real2D_onChunk ::  c_sum_litter_ag_ta
    dsl4jsb_Real2D_onChunk ::  c_sum_litter_bg_ta
    dsl4jsb_Real2D_onChunk ::  c_sum_humus_ta
    dsl4jsb_Real2D_onChunk ::  c_sum_natural_ta

    dsl4jsb_Real2D_onChunk ::  cflux_c_green_2_herb     !
    dsl4jsb_Real2D_onChunk ::  cflux_herb_2_littergreen !
    dsl4jsb_Real2D_onChunk ::  cflux_herb_2_atm         !
    dsl4jsb_Real2D_onChunk ::  co2flux_npp_2_atm_yday_ta  ! day CO2 flux from actual NPP, required for cconservation test
    dsl4jsb_Real2D_onChunk ::  co2flux_npp_2_atm_ta       ! grid cell averages of net CO2 fluxes between
    dsl4jsb_Real2D_onChunk ::  co2flux_soilresp_2_atm_ta  ! .. biosphere (due to NPP, soil respiration and
    dsl4jsb_Real2D_onChunk ::  co2flux_herb_2_atm_ta      ! .. grazing) and atmosphere [kg(CO2)/(m^2(ground) s)]


    dsl4jsb_Real2D_onChunk :: veg_fract_correction
    dsl4jsb_Real2D_onChunk :: LAI
    dsl4jsb_Real2D_onChunk :: fract_fpc_max

    ! variables needed for NLCC process
    dsl4jsb_Real2D_onChunk :: max_green_bio
    dsl4jsb_Real2D_onChunk :: sla

    ! Variables needed with forest regrowth
    dsl4jsb_Real2D_onChunk :: maxLAI_allom

    dsl4jsb_Real3D_onChunk :: c_acid_bg1_sl       !< Yasso below ground litter-pool for acid soluble litter
                                                  !< on soil layers [mol(C)/m^2(canopy)]
    dsl4jsb_Real3D_onChunk :: c_water_bg1_sl      !< Yasso below ground litter-pool for water soluble litter
                                                  !< on soil layers [mol(C)/m^2(canopy)]
    dsl4jsb_Real3D_onChunk :: c_ethanol_bg1_sl    !< Yasso below ground litter-pool for ethanol soluble litter
                                                  !< on soil layers [mol(C)/m^2(canopy)]
    dsl4jsb_Real3D_onChunk :: c_nonsoluble_bg1_sl !< Yasso below ground litter-pool for non-soluble litter
                                                  !< on soil layers [mol(C)/m^2(canopy)]
    dsl4jsb_Real3D_onChunk :: c_humus_1_sl        !< Yasso below ground litter-pool for slow C compartment
                                                  !< on soil layers [mol(C)/m^2(canopy)]
    !SIZE CLASS 2
    dsl4jsb_Real3D_onChunk :: c_acid_bg2_sl       !< Yasso below ground litter-pool for acid soluble litter
                                                  !< on soil layers [mol(C)/m^2(canopy)]
    dsl4jsb_Real3D_onChunk :: c_water_bg2_sl      !< Yasso below ground litter-pool for water soluble litter
                                                  !< on soil layers [mol(C)/m^2(canopy)]
    dsl4jsb_Real3D_onChunk :: c_ethanol_bg2_sl    !< Yasso below ground litter-pool for ethanol soluble litter
                                                  !< on soil layers [mol(C)/m^2(canopy)]
    dsl4jsb_Real3D_onChunk :: c_nonsoluble_bg2_sl !< Yasso below ground litter-pool for non-soluble litter
                                                  !< on soil layers [mol(C)/m^2(canopy)]
    dsl4jsb_Real3D_onChunk :: c_humus_2_sl        !< Yasso below ground litter-pool for slow C compartment
                                                  !< on soil layers [mol(C)/m^2(canopy)]
    dsl4jsb_Real3D_onChunk :: root_dist_sl        !< Cumulative root density distribution on soil layers [-]
    dsl4jsb_Real3D_onChunk :: max_layer_sl        !< Maximum soil layer thickness [m]

    ! climatic soil conditions
    dsl4jsb_Real3D_onChunk :: permafrost_period_sl_climbuf !< Period of consecutively froozen soil on layers [yr]

    ! active layer thickness, soil depth
    dsl4jsb_Real2D_onChunk :: max_thaw_prev_year_climbuf !< Last years maximum active layer thickness [m]
    dsl4jsb_Real2D_onChunk :: soil_depth                 !< Soil depth until bedrock [m]
    dsl4jsb_Real3D_onChunk :: soil_depth_sl              !< Soil layer depth until bedock [m]
    dsl4jsb_Real2D_onChunk :: n_thaw_cycle_prev          !< Number of thaw cycles during last year [-]

    tmp_timesteps_per_day = timesteps_per_day(options%dtime)

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(CARBON_)) RETURN

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime
    delta_time_C = sec_per_day

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    nlcc_is_enabled = model%Is_process_enabled(NLCC_)     ! workaround currently needed for DWD NEC

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(CARBON_)
    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)

    ! Prepare soil grid information for vertical yasso
    IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
      soil_c => Get_vgrid('soil_depth_carbon')
      nsoil = soil_c%n_levels

      ALLOCATE(layer_depth(nsoil))
      ALLOCATE(layer_thick(nsoil))
      layer_depth(:) = soil_c%levels
      layer_thick(:) = soil_c%ubounds-soil_c%lbounds
      ALLOCATE(soil_depth_thaw(nc))
      ALLOCATE(thaw_depth(nc))
    END IF

    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(MaxLai, old_c_state_sum_ta, current_fluxes) &
    !$ACC   CREATE(cflux_c_green_2_litter, cflux_c_wood_2_litter, cflux_crop_harvest_2_atmos) &
    !$ACC   CREATE(LeafLit_coef, WoodLit_coef, MiscPar)
    IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
      !$ACC ENTER DATA COPYIN(layer_depth, layer_thick) ASYNC(acc_stream)
    END IF

    ! The climate buffer is calculated on the veg tile, we thus need a pointer to the veg tile. We call it
    ! climbuf, as the pointer name is appended to the climate buffer variables below.
    climbuf => tile%parent_tile
    ! Make sure this is the veg tile
    IF (TRIM(climbuf%name) /= 'veg') CALL finish(routine, 'Setup ERROR: expected parent tile to be the veg tile.')
    dsl4jsb_Get_memory_tile(CLIMBUF_, climbuf)

    ! Set process variables
    dsl4jsb_Get_var2D_onChunk_tile_name(CLIMBUF_, temp_15d_yDay, climbuf)      ! in
    dsl4jsb_Get_var2D_onChunk_tile_name(CLIMBUF_, precip_15d_yDay, climbuf)    ! in
    dsl4jsb_Get_var2D_onChunk_tile_name(CLIMBUF_, max_thaw_prev_year, climbuf) ! in

    IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
      dsl4jsb_Get_var3D_onChunk_tile_name(CLIMBUF_, t_soil_sl_15d_yDay, climbuf)        ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(CLIMBUF_, wtr_soil_rel_sl_15d_yDay, climbuf)  ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(CLIMBUF_, permafrost_period_sl, climbuf)      ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, fract_org_sl)                                   ! inout
    END IF


    !
    ! This is on canopy area:
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cconservation_calcCpools )  ! out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_greenwood_2_litter ) ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest)         ! inout
    !SIZE CLASS 1
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1)            ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1)      ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1)            ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1)      ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1)             ! inout
    !SIZE CLASS 2
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2)            ! inout
    !
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_1_sum) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_2_sum) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_1_sum)   ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_2_sum)   ! inout
    ! xx
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  gross_assimilation_ca )     ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  NPP_pot_rate_ca )           ! inout
    !dsl4jsb_Get_var2D_onChunk(ASSIMI_,  NPP_pot_Rate_acc )         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  LAI_sum )                   ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_sum )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_sum )                   ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  LAI_yDayMean )              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  LAI_yyDayMean )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_yDayMean )          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_act_yDayMean )          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_yDayMean )              ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  soil_respiration )          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_flux_correction)        ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_excess)                 ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  root_exudates )             ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_sum_veg_ta )               ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_sum_litter_ag_ta )         ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_sum_litter_bg_ta )         ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_sum_humus_ta )             ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_sum_natural_ta )           ! out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_green_2_herb )       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_herb_2_littergreen  )  ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_herb_2_atm )           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  co2flux_npp_2_atm_yday_ta )  ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  co2flux_npp_2_atm_ta )       ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  co2flux_soilresp_2_atm_ta )  ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  co2flux_herb_2_atm_ta )      ! out

    ! Variables for vertical YASSO
    IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
      !SIZE CLASS 1
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_acid_bg1_sl)            ! inout
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_water_bg1_sl)           ! inout
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_ethanol_bg1_sl)         ! inout
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_nonsoluble_bg1_sl)      ! inout
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_humus_1_sl)             ! inout
      !SIZE CLASS 2
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_acid_bg2_sl)            ! inout
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_water_bg2_sl)           ! inout
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_ethanol_bg2_sl)         ! inout
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_nonsoluble_bg2_sl)      ! inout
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_humus_2_sl)             ! inout
      !
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   soil_depth)               ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   soil_depth_sl)            ! in

      dsl4jsb_Get_var3D_onChunk(CARBON_,  root_dist_sl)             ! inout
      dsl4jsb_Get_var3D_onChunk(CARBON_,  max_layer_sl)             ! inout
    END IF

    dsl4jsb_Get_var2D_onChunk(PHENO_,   veg_fract_correction )      ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_,   LAI )                       ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_,   fract_fpc_max )             ! in

    ! variables needed for NLCC process
    IF (nlcc_is_enabled) THEN
      dsl4jsb_Get_var2D_onChunk(CARBON_,  max_green_bio)            ! out
      dsl4jsb_Get_var2D_onChunk(CARBON_,  sla)                      ! out
    END IF

    IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
      dsl4jsb_Get_var2D_onChunk(SSE_,     n_thaw_cycle_prev)        ! in
    END IF

    lstart   = is_time_experiment_start(options%current_datetime)
    lresume  = is_time_restart(options%current_datetime)
    new_day  = is_newday(options%current_datetime, dtime)
    new_year = is_newyear(options%current_datetime, dtime)

    ! ---------------------------
    ! Go

    ! Put a comment on the screen...
    IF (debug_on() .AND. iblk==1 .AND. new_day  ) CALL message(TRIM(routine), &
                                                           'First time step of new day')
    IF (debug_on() .AND. iblk==1 .AND. new_year ) CALL message(TRIM(routine), &
                                                           'First time step of new year')

    ! Initializations

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      co2flux_npp_2_atm_yday_ta(ic)  = 0._lp
      co2flux_npp_2_atm_ta(ic)       = 0._lp
      co2flux_soilresp_2_atm_ta(ic)  = 0._lp
      co2flux_herb_2_atm_ta(ic)      = 0._lp
    END DO
    !$ACC END LOOP

    IF( new_year .OR. lstart) THEN
      ! Reset annual sums if newyear or if the simulation just started
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        c_decomp_humus_1_sum(ic) =  0._lp
        c_decomp_humus_2_sum(ic) =  0._lp
        c_into_humus_1_sum(ic)   =  0._lp
        c_into_humus_2_sum(ic)   =  0._lp
      END DO
      !$ACC END LOOP
    END IF
    !$ACC END PARALLEL

    ! R: "IF (read_cpools .AND. (lstart .OR. lresume)) THEN ..." of JS3 is in JS4 now in SUBROUTINE carbon_init_ic(tile)!

    ! R: To keep in mind: In JS3: NPP_pot_Rate is calculated from function NPP_pot_rate_bethy. This function gets
    !    grossAssimilation from theLand%Bethy%gross_assimilation,
    !    which is GPP relative to ground area [MOL(CO2) / M^2(ground) S] (mean value), WHICH MEANS CANOPY AREA!!!
    !    -> cbalance%NPP_pot_sum = cbalance%NPP_pot_sum + cbalance%NPP_pot_rate.
    !    AND: -> cbalance_diag%NPP_pot_yDayMean = cbalance%NPP_pot_sum/time_steps_per_day. -> cbalance_diag%NPP_pot_yDayMean
    !    is used in update_Cpools to calculate the carbon storage of the pools.
    ! R: Should be arbitary now:
    ! Compute net primary production rate
    ! R: In JS4 ist fogende Zeile jetzt im assimi_process.f90: calc_NPP_pot_rate
    !    cbalance%NPP_pot_rate(kidx0:kidx1,:) = NPP_pot_rate_bethy(grossAssimilation(1:nidx,:),darkRespiration(1:nidx,:))
    ! R: NPP_pot_Rate_acc unterscheidet sich in JS3 von NPP_pot_sum darin, dass NPP_pot_rate_acc zu jedem Zeitschritt eine
    !    momentane Aufaddierung vom NPP ist aber NPP_pot_sum immer eine Tagessumme ist.
    !    NPP_pot_Rate_acc wurde nur in den Output geschrieben aber nie weiter verwendet, daher habe ich es raus genommen.
    !    Wenn man es doch drin haben will:
    !NPP_pot_Rate_ca_acc =NPP_pot_Rate_ca_acc  +   NPP_pot_rate_ca / time_steps_per_day
    ! R: The following line should be handeled differently in JS4:
    !    In JS3 the areaWeightingFactor rescales the canopy area up to the whole Gridbox. E.g. for the C pools.
    !    In JS4 this upscaling has normally to be done by the accumulation procedures. However, the canopy area
    !    still has to be scaled to the PFT tile area for JS4. As the C pools are calculated on m^2 canopy also in JS4
    !    => Cpool(on PFT tile) =Cpool(canopy) * veg_fract_correction.
    !
    ! Therefore this gets arbitary:
    ! Prepare area weighting factor to rescale from 1/[m^2(canopy)] to 1/[m^2(grid box)]
    ! areaWeightingFactor(:,:) = veg_fract_correction(:,:) * surface%cover_fract(kidx0:kidx1,:) &
    !                            * SPREAD(surface%veg_ratio_max(kidx0:kidx1),DIM=2,NCOPIES=ntiles)

    ! Get information about the maximum thaw depth on soil layers
    IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
      soil_depth_thaw = MIN(soil_depth, max_thaw_prev_year_climbuf)
      max_layer_sl = soil_depth_to_layers(soil_depth_thaw, layer_thick)
      thaw_depth = max_thaw_prev_year_climbuf
      !$ACC ENTER DATA COPYIN(soil_depth_thaw, thaw_depth) ASYNC(acc_stream)
    END IF

    ! All per tile variables need to be re-calculated in case they are written to output
    CALL recalc_per_tile_vars(tile, options,                                      &
      & c_green = c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:),     &
      & c_crop_harvest = c_crop_harvest(:),                                       &
      & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:),                 &
      & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
      & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:),                 &
      & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
      & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:),                 &
      & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
      & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:),                 &
      & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
      & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:),                       &
      & root_exudates = root_exudates(:),                                         &
      & soil_respiration = soil_respiration(:),                                   &
      & NPP_flux_correction = NPP_flux_correction(:),                             &
      & cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter(:),               &
      & cflux_c_green_2_herb = cflux_c_green_2_herb(:),                           &
      & NPP_act_yDayMean = NPP_act_yDayMean(:),                                   &
      & NPP_pot_yDayMean = NPP_pot_yDayMean(:),                                   &
      & GPP_yDayMean = GPP_yDayMean(:))

    IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
      ! Recalc per tile variables on soil levels
      CALL recalc_per_tile_vars(tile, options,                                                      &
        & c_acid_bg1_sl = c_acid_bg1_sl(:,:), c_water_bg1_sl = c_water_bg1_sl(:,:),                 &
        & c_ethanol_bg1_sl = c_ethanol_bg1_sl(:,:), c_nonsoluble_bg1_sl = c_nonsoluble_bg1_sl(:,:), &
        & c_acid_bg2_sl = c_acid_bg2_sl(:,:), c_water_bg2_sl = c_water_bg2_sl(:,:),                 &
        & c_ethanol_bg2_sl = c_ethanol_bg2_sl(:,:), c_nonsoluble_bg2_sl = c_nonsoluble_bg2_sl(:,:), &
        & c_humus_1_sl = c_humus_1_sl(:,:), c_humus_2_sl = c_humus_2_sl(:,:))
    END IF

    IF( .NOT. new_day .OR. lstart) THEN ! perform daily sums if WE ARE NOT STARTING A NEW DAY
                                            ! or if we JUST STARTED AN TOTALLY NEW EXPERIMENT
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        LAI_sum(ic)      = LAI_sum(ic)      + LAI(ic)
        NPP_pot_sum(ic)  = NPP_pot_sum(ic)  + NPP_pot_rate_ca(ic)
        GPP_sum(ic)      = GPP_sum(ic)      + gross_assimilation_ca(ic)
      END DO
      !$ACC END PARALLEL LOOP

    ELSE                                    ! A NEW DAY BEGINS and we DID NOT START AN
                                            ! TOTALLY NEW EXPERIMENT ==> perform carbon balance
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        ! First save before yesterday's and yesterday's means
        LAI_yyDayMean(ic)  = LAI_yDayMean(ic)
        LAI_yDayMean(ic)   = LAI_sum(ic)/tmp_timesteps_per_day

        NPP_pot_yDayMean(ic)     = NPP_pot_sum(ic)/tmp_timesteps_per_day
        GPP_yDayMean(ic)         = GPP_sum(ic)/tmp_timesteps_per_day

        ! Then start summing up today's values
        LAI_sum(ic)             = LAI(ic)
        NPP_pot_sum(ic)         = NPP_pot_rate_ca(ic)
        GPP_sum(ic)             = gross_assimilation_ca(ic)
        !NPP_pot_rate_ca_acc(ic)        = NPP_pot_rate_ca(ic)/tmp_timesteps_per_day ! s.o.
      END DO
      !$ACC END PARALLEL LOOP

      ! With forest regrowth the maximum LAI depends on location and can change over time.
      IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
        dsl4jsb_Get_var2D_onChunk(PHENO_, maxLAI_allom)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          MaxLai(ic) = maxLAI_allom(ic)
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        lctlib_MaxLAI = dsl4jsb_Lctlib_param(MaxLAI)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          MaxLai(ic) = lctlib_MaxLAI
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      ! Helper parameters needed with GPUs
      lctlib_fract_npp_2_woodPool    = dsl4jsb_Lctlib_param(fract_npp_2_woodPool)
      lctlib_fract_NPP_2_reservePool = dsl4jsb_Lctlib_param(fract_NPP_2_reservePool)
      lctlib_fract_NPP_2_exudates    = dsl4jsb_Lctlib_param(fract_NPP_2_exudates)
      lctlib_fract_green_2_herbivory = dsl4jsb_Lctlib_param(fract_green_2_herbivory)
      lctlib_tau_c_woods             = dsl4jsb_Lctlib_param(tau_c_woods)
      lctlib_LAI_shed_constant       = dsl4jsb_Lctlib_param(LAI_shed_constant)
      lctlib_Max_C_content_woods     = dsl4jsb_Lctlib_param(Max_C_content_woods)
      lctlib_specificLeafArea_C      = dsl4jsb_Lctlib_param(specificLeafArea_C)
      lctlib_reserveC2leafC          = dsl4jsb_Lctlib_param(reserveC2leafC)
      lctlib_PhenologyType           = dsl4jsb_Lctlib_param(PhenologyType)

      ! Set up parameter arrays for cleaner parameter transfer vertical YASSO
      LeafLit_coef(Y_ACID)       = dsl4jsb_Lctlib_param(LeafLit_coef(Y_ACID))
      LeafLit_coef(Y_WATER)      = dsl4jsb_Lctlib_param(LeafLit_coef(Y_WATER))
      LeafLit_coef(Y_ETHANOL)    = dsl4jsb_Lctlib_param(LeafLit_coef(Y_ETHANOL))
      LeafLit_coef(Y_NONSOLUBLE) = dsl4jsb_Lctlib_param(LeafLit_coef(Y_NONSOLUBLE))
      LeafLit_coef(Y_HUMUS)      = dsl4jsb_Lctlib_param(LeafLit_coef(Y_HUMUS))
      WoodLit_coef(Y_ACID)       = dsl4jsb_Lctlib_param(WoodLit_coef(Y_ACID))
      WoodLit_coef(Y_WATER)      = dsl4jsb_Lctlib_param(WoodLit_coef(Y_WATER))
      WoodLit_coef(Y_ETHANOL)    = dsl4jsb_Lctlib_param(WoodLit_coef(Y_ETHANOL))
      WoodLit_coef(Y_NONSOLUBLE) = dsl4jsb_Lctlib_param(WoodLit_coef(Y_NONSOLUBLE))
      WoodLit_coef(Y_HUMUS)      = dsl4jsb_Lctlib_param(WoodLit_coef(Y_HUMUS))
      MiscPar(WOODLITSIZE)         = dsl4jsb_Lctlib_param(WoodLitterSize)
      MiscPar(ROOT_DIST_BETA)      = dsl4jsb_Lctlib_param(root_dist_beta)

      ! Check whether optional values were read from lctlib
      IF (MiscPar(ROOT_DIST_BETA) < -990._lp .AND. dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
        CALL finish(TRIM(routine), 'ROOT_DIST_BETA is required in lctlib file for simulations ' // &
          & 'using vertically stratified yasso scheme')
      END IF

      !$ACC UPDATE DEVICE(LeafLit_coef, WoodLit_coef, MiscPar) ASYNC(acc_stream)

      ! Calculate current carbon sum before the update - for C conservation test
      CALL recalc_carbon_per_tile_vars(tile, options)
      CALL calculate_current_c_ta_state_sum(tile, options, old_c_state_sum_ta(:))

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        ! Now update the C pools for plants
        CALL calc_plant_c( &
          &  delta_time_C,                                      & ! in
          &  LAI_yDayMean(ic),                                  & ! in
          &  LAI_yyDayMean(ic),                                 & ! in
          &  MaxLai(ic),                                        & ! in
          &  NPP_pot_yDayMean(ic),                              & ! in
          !
          &  lctlib_fract_npp_2_woodPool,                       & ! in
          &  lctlib_fract_NPP_2_reservePool,                    & ! in
          &  lctlib_fract_NPP_2_exudates,                       & ! in
          &  lctlib_fract_green_2_herbivory,                    & ! in
          &  lctlib_tau_c_woods,                                & ! in
          &  lctlib_LAI_shed_constant,                          & ! in
          &  lctlib_Max_C_content_woods,                        & ! in
          &  lctlib_specificLeafArea_C,                         & ! in
          &  lctlib_reserveC2leafC,                             & ! in
          &  lctlib_PhenologyType,                              & ! in
          !
          &  c_green(ic),                                       & ! inout
          &  c_woods(ic),                                       & ! inout
          &  c_reserve(ic),                                     & ! inout
          &  c_crop_harvest(ic),                                & ! inout
          !
          &  NPP_flux_correction(ic),                           & ! out
          &  NPP_excess(ic),                                    & ! out
          &  NPP_act_yDayMean(ic),                              & ! out
          &  root_exudates(ic),                                 & ! out
          &  cflux_c_greenwood_2_litter(ic),                    & ! out
          &  cflux_c_green_2_litter(ic),                        & ! out
          &  cflux_c_wood_2_litter(ic),                         & ! out
          &  cflux_c_green_2_herb(ic),                          & ! out
          &  cflux_herb_2_littergreen(ic),                      & ! out
          &  cflux_herb_2_atm(ic),                              & ! out
          &  cflux_crop_harvest_2_atmos(ic)                     & ! out
          & )
      END DO
      !$ACC END PARALLEL LOOP

      ! clear C from frozen pools if requested
      IF ((dsl4jsb_Config(CARBON_)%vertical_yasso) .AND. (dsl4jsb_Config(CARBON_)%remove_permafrost_carbon) &
        &  .AND. (lresume .OR. lstart)) THEN
        CALL message(TRIM(routine), 'Removing carbon from permafrost affected soil layers')
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
        DO ic = 1, nc
          DO is = 1, nsoil
            IF (permafrost_period_sl_climbuf(ic,is) > permafrost_threshold) THEN
              c_acid_bg1_sl(ic,is)       = 0._lp
              c_water_bg1_sl(ic,is)      = 0._lp
              c_ethanol_bg1_sl(ic,is)    = 0._lp
              c_nonsoluble_bg1_sl(ic,is) = 0._lp
              c_acid_bg2_sl(ic,is)       = 0._lp
              c_water_bg2_sl(ic,is)      = 0._lp
              c_ethanol_bg2_sl(ic,is)    = 0._lp
              c_nonsoluble_bg2_sl(ic,is) = 0._lp
              c_humus_1_sl(ic,is)        = 0._lp
              c_humus_2_sl(ic,is)        = 0._lp
            END IF
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      ! Calculate soil carbon pools
      IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
        ! Use multi-layer yasso scheme
        CALL calc_soil_c(                                                                              &
          ! Input (parameter, soil properties, climate forcing, carbon fluxes)
          & nc, nsoil, delta_time_C, LeafLit_coef(:), WoodLit_coef(:), MiscPar(:),                     &
          & soil_depth(:), soil_depth_sl(:,:), soil_depth_thaw(:), max_layer_sl(:,:), layer_depth(:),  &
          & layer_thick(:), thaw_depth(:), permafrost_period_sl_climbuf(:,:),                          &
          & temp_15d_yDay_climbuf(:), precip_15d_yDay_climbuf(:), t_soil_sl_15d_yDay_climbuf(:,:),     &
          & wtr_soil_rel_sl_15d_yDay_climbuf(:,:), n_thaw_cycle_prev(:),                               &
          & cflux_c_wood_2_litter(:), cflux_c_green_2_litter(:), root_exudates(:),                     &
          ! Inout (Carbon pools)
          & c_acid_ag1(:),      c_water_ag1(:),      c_ethanol_ag1(:),      c_nonsoluble_ag1(:),       &
          & c_acid_ag2(:),      c_water_ag2(:),      c_ethanol_ag2(:),      c_nonsoluble_ag2(:),       &
          & c_acid_bg1(:),      c_water_bg1(:),      c_ethanol_bg1(:),      c_nonsoluble_bg1(:),       &
          & c_acid_bg2(:),      c_water_bg2(:),      c_ethanol_bg2(:),      c_nonsoluble_bg2(:),       &
          & c_acid_bg1_sl(:,:), c_water_bg1_sl(:,:), c_ethanol_bg1_sl(:,:), c_nonsoluble_bg1_sl(:,:),  &
          & c_acid_bg2_sl(:,:), c_water_bg2_sl(:,:), c_ethanol_bg2_sl(:,:), c_nonsoluble_bg2_sl(:,:),  &
          & c_humus_1(:),       c_humus_2(:),        c_humus_1_sl(:,:),     c_humus_2_sl(:,:),         &
          ! Output (Soil respiration and humus fluxes)
          & c_into_humus_1_sum(:), c_decomp_humus_1_sum(:),                                            &
          & c_into_humus_2_sum(:), c_decomp_humus_2_sum(:), soil_respiration(:)                        &
          & )
      ELSE
        ! Use single-layer yasso scheme
        CALL calc_soil_c(                                                          &
          ! Input (Parameter, Climate forcing, Carbon fluxes)
          & nc, delta_time_C, LeafLit_coef(:), WoodLit_coef(:), MiscPar(:),        &
          & temp_15d_yDay_climbuf(:), precip_15d_yDay_climbuf(:),                  &
          & cflux_c_wood_2_litter(:), cflux_c_green_2_litter(:), root_exudates(:), &
          ! Inout (Carbon pools)
          & c_acid_ag1(:), c_water_ag1(:), c_ethanol_ag1(:), c_nonsoluble_ag1(:),  &
          & c_acid_ag2(:), c_water_ag2(:), c_ethanol_ag2(:), c_nonsoluble_ag2(:),  &
          & c_acid_bg1(:), c_water_bg1(:), c_ethanol_bg1(:), c_nonsoluble_bg1(:),  &
          & c_acid_bg2(:), c_water_bg2(:), c_ethanol_bg2(:), c_nonsoluble_bg2(:),  &
          & c_humus_1(:),  c_humus_2(:),                                           &
          ! Output (Soil respiration and humus fluxes)
          & c_into_humus_1_sum(:), c_decomp_humus_1_sum(:),                        &
          & c_into_humus_2_sum(:), c_decomp_humus_2_sum(:), soil_respiration(:)    &
          & )
      END IF

      ! Reset humus accumulation diagnostics for too small tiles
      tile_fract => tile%fract(ics:ice, iblk)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        IF (tile_fract(ic) < 1e-10_lp) THEN
          c_decomp_humus_1_sum(ic) = 0._lp
          c_decomp_humus_2_sum(ic) = 0._lp
          c_into_humus_1_sum(ic)   = 0._lp
          c_into_humus_2_sum(ic)   = 0._lp
        END IF
      END DO
      !$ACC END PARALLEL LOOP

      !! ============== Include decay of crop harvest in soil respiration =============
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        soil_respiration(ic) = soil_respiration(ic) + cflux_crop_harvest_2_atmos(ic)
      END DO
      !$ACC END PARALLEL LOOP

      ! determine annual maximum content of green pool for NLCC process
      IF (model%Is_process_enabled(NLCC_)) THEN
        IF (new_year) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
          DO ic = 1, nc
            max_green_bio(ic) = c_green(ic)
          END DO
          !$ACC END PARALLEL LOOP
        ELSE IF (new_day) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
          DO ic = 1, nc
            max_green_bio(ic) = MAX(max_green_bio(ic), c_green(ic))
          END DO
          !$ACC END PARALLEL LOOP
        END IF
      END IF

      ! Compute carbon contents from canopy to tile box area by weighting pools with fractions of grid box covered by vegetation
      ! ------------------------------------------------------------------------------------------------------------------------

      CALL recalc_per_tile_vars(tile, options,                                      &
        & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:),      &
        & c_crop_harvest = c_crop_harvest(:),                                       &
        & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:),                 &
        & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
        & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:),                 &
        & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
        & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:),                 &
        & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
        & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:),                 &
        & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
        & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:),                       &
        & c_decomp_humus_1_sum = c_decomp_humus_1_sum(:),                           &
        & c_decomp_humus_2_sum = c_decomp_humus_2_sum(:),                           &
        & c_into_humus_1_sum = c_into_humus_1_sum(:),                               &
        & c_into_humus_2_sum = c_into_humus_2_sum(:),                               &
        & root_exudates = root_exudates(:),                                         &
        & soil_respiration = soil_respiration(:),                                   &
        & NPP_flux_correction = NPP_flux_correction(:),                             &
        & cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter(:),               &
        & cflux_c_green_2_herb = cflux_c_green_2_herb(:),                           &
        & NPP_act_yDayMean = NPP_act_yDayMean(:),                                   &
        & NPP_pot_yDayMean = NPP_pot_yDayMean(:),                                   &
        & GPP_yDayMean = GPP_yDayMean(:))

      IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
        ! Recalc per tile variables on soil levels
        CALL recalc_per_tile_vars(tile, options,                                                      &
          & c_acid_bg1_sl = c_acid_bg1_sl(:,:), c_water_bg1_sl = c_water_bg1_sl(:,:),                 &
          & c_ethanol_bg1_sl = c_ethanol_bg1_sl(:,:), c_nonsoluble_bg1_sl = c_nonsoluble_bg1_sl(:,:), &
          & c_acid_bg2_sl = c_acid_bg2_sl(:,:), c_water_bg2_sl = c_water_bg2_sl(:,:),                 &
          & c_ethanol_bg2_sl = c_ethanol_bg2_sl(:,:), c_nonsoluble_bg2_sl = c_nonsoluble_bg2_sl(:,:), &
          & c_humus_1_sl = c_humus_1_sl(:,:), c_humus_2_sl = c_humus_2_sl(:,:))
      END IF

      ! cflux_herb_2_atm and soil_respiration are negative (fluxes away from land)
      ! thus + here instead of - to reduce NPP by these fluxes
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        current_fluxes(ic) = (NPP_act_yDayMean(ic) + cflux_herb_2_atm(ic) + soil_respiration(ic)) &
          &                   * sec_per_day * veg_fract_correction(ic) * fract_fpc_max(ic)
      END DO
      !$ACC END PARALLEL LOOP

      CALL check_carbon_conservation(tile, options, old_c_state_sum_ta(:), &
        & current_fluxes(:), cconservation_calcCpools(:))

      ! Update relative soil carbon content on soil layers (if not prescribed by map)
      IF (dsl4jsb_Config(CARBON_)%vertical_yasso .AND. dsl4jsb_Config(HYDRO_)%l_organic &
        & .AND. .NOT. dsl4jsb_Config(HYDRO_)%l_socmap) THEN
        ! Accumulated below ground soil carbon content, convert to organic matter mass and
        ! compute its relative fraction wrt total soil mass.
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream) &
        !$ACC   PRIVATE(mass_organic_soil)
        DO ic = 1, nc
          DO is = 1, nsoil
            mass_organic_soil   = molarMassC_kg / C_2_biomass_ratio * (                   &
              &                   c_humus_1_sl(ic,is)        + c_humus_2_sl(ic,is)        &
              &                 + c_acid_bg1_sl(ic,is)       + c_acid_bg2_sl(ic,is)       &
              &                 + c_ethanol_bg1_sl(ic,is)    + c_ethanol_bg2_sl(ic,is)    &
              &                 + c_water_bg1_sl(ic,is)      + c_water_bg2_sl(ic,is)      &
              &                 + c_nonsoluble_bg1_sl(ic,is) + c_nonsoluble_bg2_sl(ic,is) &
              &                 ) * veg_fract_correction(ic) * fract_fpc_max(ic)
            fract_org_sl(ic,is) = MAX(0._lp, MIN(1.0_lp, &
              & mass_organic_soil / (rho_org_bg * layer_thick(is))))
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END IF

    END IF ! new day

    ! Compute net CO2 fluxes exchanged with atmosphere at each time step
    !-------------------------------------------------------------------
    ! Note: carbon loss of biosphere means a positive CO2 flux to atmosphere (i.e. NEP and net CO2-flux have opposite signs)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      co2flux_npp_2_atm_ta(ic) = molarMassCO2_kg *                    & ! Conversion factor from mol to kg CO2
                 (- veg_fract_correction(ic) * fract_fpc_max(ic)      & ! Minus: atmosphere gain is positive
                      * (NPP_pot_rate_ca(ic)                          & ! current (not actual) NPP rate
                          - (NPP_pot_yDayMean(ic)-NPP_act_yDayMean(ic))) )  ! corrected with yesterdays actual NPP defizit

      co2flux_soilresp_2_atm_ta(ic) = molarMassCO2_kg *               & ! Conversion factor from mol to kg CO2
                 (- veg_fract_correction(ic) *fract_fpc_max(ic)       & ! Minus: atmosphere gain is positive
                      * soil_respiration(ic)      )                     ! .. soil respiration

      co2flux_herb_2_atm_ta(ic) = molarMassCO2_kg *                   & ! Conversion factor from mol to kg CO2
                 (- veg_fract_correction(ic) * fract_fpc_max(ic)      & ! Minus: atmosphere gain is positive
                      * cflux_herb_2_atm(ic) )                          ! .. herbivory

      co2flux_npp_2_atm_yday_ta(ic) = molarMassCO2_kg *               & ! daily C conservation test cannot deal with
                (- veg_fract_correction(ic) * fract_fpc_max(ic)       & ! diurnal cycle of NPP,
                     * (NPP_act_yDayMean(ic)))                          ! therefore here additionaly the day mean
    END DO
    !$ACC END PARALLEL LOOP

    ! Calculate diagnostic carbon sums
    CALL calculate_current_c_ta_state_sum(tile, options, c_sum_natural_ta(:),            &
      & c_sum_veg_ta(:), c_sum_litter_ag_ta(:), c_sum_litter_bg_ta(:), c_sum_humus_ta(:))

    ! Write lct information on variable to make it available on pft-tiles via function collect_var for NLCC process
    IF (model%Is_process_enabled(NLCC_)) THEN
      lctlib_specificLeafArea_C = dsl4jsb_Lctlib_param(specificLeafArea_C)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        sla(ic) = lctlib_specificLeafArea_C
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
      !$ACC EXIT DATA DELETE(layer_depth, layer_thick, soil_depth_thaw, thaw_depth) ASYNC(acc_stream)
    END IF
    !$ACC END DATA
    IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
      DEALLOCATE(layer_depth)
      DEALLOCATE(layer_thick)
      DEALLOCATE(soil_depth_thaw)
      DEALLOCATE(thaw_depth)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_C_NPP_pot_allocation

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "C_NPP_pot_allocation"
  !!
  !! @param[in,out] tile    Tile for which aggregation of child tiles is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_C_NPP_pot_allocation(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_memory(CARBON_)

    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract
    TYPE(t_jsb_model), POINTER   :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_C_NPP_pot_allocation'

    INTEGER  :: iblk , ics, ice

    model => Get_model(tile%owner_model_id)
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_config(CARBON_)
    dsl4jsb_Get_memory(CARBON_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yDayMean,          weighted_by_fract)

    !@todo Those commented out are currently not calculated
    ! non-yasso-pools
    dsl4jsb_Aggregate_onChunk(CARBON_, c_green_ta,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_woods_ta,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_reserve_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_crop_harvest_ta,     weighted_by_fract)
    ! ag1 and bg1
    dsl4jsb_Aggregate_onChunk(CARBON_, c_acid_ag1_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_water_ag1_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_ethanol_ag1_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_nonsoluble_ag1_ta,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_acid_bg1_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_water_bg1_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_ethanol_bg1_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_nonsoluble_bg1_ta,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_humus_1_ta,          weighted_by_fract)
    !ag2 and bg2
    dsl4jsb_Aggregate_onChunk(CARBON_, c_acid_ag2_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_water_ag2_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_ethanol_ag2_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_nonsoluble_ag2_ta,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_acid_bg2_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_water_bg2_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_ethanol_bg2_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_nonsoluble_bg2_ta,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_humus_2_ta,          weighted_by_fract)

    !JN: vars for analytically equilibrating YASSO humus pools outside of jsb4 simulations (compare jsb3 svn ref 9652)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_decomp_humus_1_sum_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_decomp_humus_2_sum_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_into_humus_1_sum_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_into_humus_2_sum_ta,          weighted_by_fract)

    ! sums
    dsl4jsb_Aggregate_onChunk(CARBON_, c_sum_veg_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_sum_litter_ag_ta,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_sum_litter_bg_ta,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_sum_humus_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_sum_natural_ta,      weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_pot_sum_ta,        weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, GPP_sum_ta,            weighted_by_fract)
    ! fluxes
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_c_greenwood_2_litter_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, soil_respiration_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_flux_correction_ta,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_excess_ta,            weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yyDayMean_ta,         weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yDayMean_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_pot_yDayMean_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_act_yDayMean_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, GPP_yDayMean_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, root_exudates_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_c_green_2_herb_ta,  weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_, cflux_herb_2_littergreen_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cconservation_calcCpools, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, co2flux_npp_2_atm_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, co2flux_soilresp_2_atm_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, co2flux_herb_2_atm_ta,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, co2flux_npp_2_atm_yday_ta, weighted_by_fract)

    ! variables for NLCC process
    dsl4jsb_Aggregate_onChunk(CARBON_, max_green_bio,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, sla,                      weighted_by_fract)

    IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
      ! variables on soil levels
      ! bg1
      dsl4jsb_Aggregate_onChunk(CARBON_, c_acid_bg1_sl_ta,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_, c_water_bg1_sl_ta,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_, c_ethanol_bg1_sl_ta,      weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_, c_nonsoluble_bg1_sl_ta,   weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_, c_humus_1_sl_ta,          weighted_by_fract)
      ! bg2
      dsl4jsb_Aggregate_onChunk(CARBON_, c_acid_bg2_sl_ta,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_, c_water_bg2_sl_ta,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_, c_ethanol_bg2_sl_ta,      weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_, c_nonsoluble_bg2_sl_ta,   weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_, c_humus_2_sl_ta,          weighted_by_fract)

      dsl4jsb_Aggregate_onChunk(CARBON_, root_dist_sl,             weighted_by_fract)             ! inout
      dsl4jsb_Aggregate_onChunk(CARBON_, max_layer_sl,             weighted_by_fract)             ! inout
    END IF

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_C_NPP_pot_allocation

  ! ====================================================================================================== !
  !
  !> Implementation of "update" for task "ta_vars_after_lcc"
  !>
  !> Task "ta_vars_after_lcc" updates the ta variables
  !
  SUBROUTINE update_ta_vars_after_lcc(tile, options)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! -------------------------------------------------------------------------------------------------- !

    ! ... currently not required, but could contain something like:
    ! CALL recalc_carbon_per_tile_vars(tile, options)

  END SUBROUTINE update_ta_vars_after_lcc

  ! ====================================================================================================== !
  !
  !> Implementation of "aggregate" for task "ta_vars_after_lcc"
  !>
  !> Aggregate for task "ta_vars_after_lcc" is currently required to aggregate the ta for the co2flux
  !
  ! -------------------------------------------------------------------------------------------------------!
  SUBROUTINE aggregate_ta_vars_after_lcc(tile, options)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! -------------------------------------------------------------------------------------------------- !

    dsl4jsb_Def_memory(CARBON_)

    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_ta_vars_after_lcc'

    INTEGER  :: iblk, ics, ice

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(CARBON_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    ! if running with FLCC (which runs on the veg tile)
    ! the co2flux_fire_all_2_atm_ta variable is changed on the leaves
    ! ... it therefore requires a carbon aggregation task on the leaves to run subsequent to the flcc task
    dsl4jsb_Aggregate_onChunk(CARBON_,  co2flux_fire_all_2_atm_ta, weighted_by_fract)

    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_green_ta,              weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_woods_ta,              weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_reserve_ta,            weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_dist_green_2_soil_ta, weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_dist_woods_2_soil_ta, weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_ag1_ta,           weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_ag1_ta,          weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_ag1_ta,        weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_ag1_ta,     weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_ag2_ta,           weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_ag2_ta,          weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_ag2_ta,        weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_ag2_ta,     weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_humus_1_ta,            weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_,  c_humus_2_ta,            weighted_by_fract)

    CALL weighted_by_fract%EndAggregate()

  END SUBROUTINE aggregate_ta_vars_after_lcc


!##########################################################################################################################
! other subroutines (no further tasks)
!##########################################################################################################################

  ! ================================================================================================================================
  !>
  !> scales all carbon state variables upon reference area change
  !!
  !! @todo discuss how to make more general, e.g. get the info on which state variables to scale from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE rescale_carbon_upon_reference_area_change(tile, options, oldRefArea, newRefArea)

    USE mo_util,                ONLY: real2string

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)            :: options
    REAL(lp),                   INTENT(in)            :: oldRefArea(:),  &
                                                         newRefArea(:)
    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_memory(CARBON_)

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//':rescale_carbon_upon_reference_area_change'
    REAL(lp)                    :: refarea_scaling
    INTEGER                     :: nc, ic, ics, ice, iblk
    TYPE(t_jsb_model), POINTER  :: model

    dsl4jsb_Real2D_onChunk ::  c_green
    dsl4jsb_Real2D_onChunk ::  c_reserve
    dsl4jsb_Real2D_onChunk ::  c_woods
    dsl4jsb_Real2D_onChunk ::  c_crop_harvest

    dsl4jsb_Real2D_onChunk ::  c_acid_ag1
    dsl4jsb_Real2D_onChunk ::  c_water_ag1
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag1
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag1
    dsl4jsb_Real2D_onChunk ::  c_acid_bg1
    dsl4jsb_Real2D_onChunk ::  c_water_bg1
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg1
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg1
    dsl4jsb_Real2D_onChunk ::  c_humus_1

    dsl4jsb_Real2D_onChunk ::  c_acid_ag2
    dsl4jsb_Real2D_onChunk ::  c_water_ag2
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag2
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag2
    dsl4jsb_Real2D_onChunk ::  c_acid_bg2
    dsl4jsb_Real2D_onChunk ::  c_water_bg2
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg2
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg2
    dsl4jsb_Real2D_onChunk ::  c_humus_2

    dsl4jsb_Real2D_onChunk ::  c_decomp_humus_1_sum
    dsl4jsb_Real2D_onChunk ::  c_decomp_humus_2_sum
    dsl4jsb_Real2D_onChunk ::  c_into_humus_1_sum
    dsl4jsb_Real2D_onChunk ::  c_into_humus_2_sum

    dsl4jsb_Real3D_onChunk ::  c_acid_bg1_sl
    dsl4jsb_Real3D_onChunk ::  c_water_bg1_sl
    dsl4jsb_Real3D_onChunk ::  c_ethanol_bg1_sl
    dsl4jsb_Real3D_onChunk ::  c_nonsoluble_bg1_sl
    dsl4jsb_Real3D_onChunk ::  c_humus_1_sl

    dsl4jsb_Real3D_onChunk ::  c_acid_bg2_sl
    dsl4jsb_Real3D_onChunk ::  c_water_bg2_sl
    dsl4jsb_Real3D_onChunk ::  c_ethanol_bg2_sl
    dsl4jsb_Real3D_onChunk ::  c_nonsoluble_bg2_sl
    dsl4jsb_Real3D_onChunk ::  c_humus_2_sl

    ! Assertion: 0 < oldRedArea <= 1 and 0 < newRefArea <= 1
    ! TODO on GPU
#ifndef _OPENACC
    IF ( ANY(oldRefArea .LE. 0) .OR. ANY(oldRefArea .GT. 1) ) THEN
      CALL finish(TRIM(routine), &
        & 'Violation of assertion: Reference areas (here old) need to be >0 and <=1 - min: '&
        & //real2string(MINVAL(oldRefArea)) //'; max: '//real2string(MAXVAL(oldRefArea))  )
    END IF
    IF ( ANY(newRefArea .LE. 0) .OR. ANY(newRefArea .GT. 1) ) THEN
      CALL finish(TRIM(routine), &
        & 'Violation of assertion: Reference areas (here new) need to be >0 and <=1 - min: '&
        & //real2string(MINVAL(newRefArea)) //'; max: '//real2string(MAXVAL(newRefArea))  )
    END IF
#endif

    model => Get_model(tile%owner_model_id)
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    dsl4jsb_Get_config(CARBON_)
    dsl4jsb_Get_memory(CARBON_)

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve )           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest)       ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1)            ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2)            ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_1_sum) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_2_sum) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_1_sum)   ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_2_sum)   ! inout

    IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_acid_bg1_sl )
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_water_bg1_sl )
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_ethanol_bg1_sl )
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_nonsoluble_bg1_sl )
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_humus_1_sl )

      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_acid_bg2_sl )
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_water_bg2_sl )
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_ethanol_bg2_sl )
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_nonsoluble_bg2_sl )
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_humus_2_sl )
    END IF

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(refarea_scaling)
    DO ic = 1, nc
      refarea_scaling = oldRefArea(ic) / newRefArea(ic)
      c_green(ic)          = c_green(ic)          * refarea_scaling
      c_woods(ic)          = c_woods(ic)          * refarea_scaling
      c_reserve(ic)        = c_reserve(ic)        * refarea_scaling
      c_crop_harvest(ic)   = c_crop_harvest(ic)   * refarea_scaling

      c_acid_ag1(ic)       = c_acid_ag1(ic)       * refarea_scaling
      c_water_ag1(ic)      = c_water_ag1(ic)      * refarea_scaling
      c_ethanol_ag1(ic)    = c_ethanol_ag1(ic)    * refarea_scaling
      c_nonsoluble_ag1(ic) = c_nonsoluble_ag1(ic) * refarea_scaling
      c_acid_bg1(ic)       = c_acid_bg1(ic)       * refarea_scaling
      c_water_bg1(ic)      = c_water_bg1(ic)      * refarea_scaling
      c_ethanol_bg1(ic)    = c_ethanol_bg1(ic)    * refarea_scaling
      c_nonsoluble_bg1(ic) = c_nonsoluble_bg1(ic) * refarea_scaling
      c_humus_1(ic)        = c_humus_1(ic)        * refarea_scaling
      c_acid_ag2(ic)       = c_acid_ag2(ic)       * refarea_scaling
      c_water_ag2(ic)      = c_water_ag2(ic)      * refarea_scaling
      c_ethanol_ag2(ic)    = c_ethanol_ag2(ic)    * refarea_scaling
      c_nonsoluble_ag2(ic) = c_nonsoluble_ag2(ic) * refarea_scaling
      c_acid_bg2(ic)       = c_acid_bg2(ic)       * refarea_scaling
      c_water_bg2(ic)      = c_water_bg2(ic)      * refarea_scaling
      c_ethanol_bg2(ic)    = c_ethanol_bg2(ic)    * refarea_scaling
      c_nonsoluble_bg2(ic) = c_nonsoluble_bg2(ic) * refarea_scaling
      c_humus_2(ic)        = c_humus_2(ic)        * refarea_scaling

      c_decomp_humus_1_sum(ic) = c_decomp_humus_1_sum(ic) * refarea_scaling
      c_decomp_humus_2_sum(ic) = c_decomp_humus_2_sum(ic) * refarea_scaling
      c_into_humus_1_sum(ic)   = c_into_humus_1_sum(ic)   * refarea_scaling
      c_into_humus_2_sum(ic)   = c_into_humus_2_sum(ic)   * refarea_scaling
    END DO
    !$ACC END PARALLEL LOOP

    IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
      DO ic = 1, nc
        c_acid_bg1_sl(ic,:)       = c_acid_bg1_sl(ic,:)       * refarea_scaling
        c_water_bg1_sl(ic,:)      = c_water_bg1_sl(ic,:)      * refarea_scaling
        c_ethanol_bg1_sl(ic,:)    = c_ethanol_bg1_sl(ic,:)    * refarea_scaling
        c_nonsoluble_bg1_sl(ic,:) = c_nonsoluble_bg1_sl(ic,:) * refarea_scaling
        c_humus_1_sl(ic,:)        = c_humus_1_sl(ic,:)        * refarea_scaling
        c_acid_bg2_sl(ic,:)       = c_acid_bg2_sl(ic,:)       * refarea_scaling
        c_water_bg2_sl(ic,:)      = c_water_bg2_sl(ic,:)      * refarea_scaling
        c_ethanol_bg2_sl(ic,:)    = c_ethanol_bg2_sl(ic,:)    * refarea_scaling
        c_nonsoluble_bg2_sl(ic,:) = c_nonsoluble_bg2_sl(ic,:) * refarea_scaling
        c_humus_2_sl(ic,:)        = c_humus_2_sl(ic,:)        * refarea_scaling
      END DO
    END IF

    CALL recalc_per_tile_vars(tile, options, &
      & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:), &
      & c_crop_harvest = c_crop_harvest(:), &
      & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:), &
      & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
      & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:), &
      & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
      & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:), &
      & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
      & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:), &
      & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
      & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:),                       &
      & c_decomp_humus_1_sum  = c_decomp_humus_1_sum(:),                          &
      & c_decomp_humus_2_sum  = c_decomp_humus_2_sum(:),                          &
      & c_into_humus_1_sum    = c_into_humus_1_sum(:),                            &
      & c_into_humus_2_sum    = c_into_humus_2_sum(:))

    IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
      CALL recalc_per_tile_vars(tile, options,                                                    &
      & c_acid_bg1_sl = c_acid_bg1_sl(:,:), c_water_bg1_sl = c_water_bg1_sl(:,:),                 &
      & c_ethanol_bg1_sl = c_ethanol_bg1_sl(:,:), c_nonsoluble_bg1_sl = c_nonsoluble_bg1_sl(:,:), &
      & c_acid_bg2_sl = c_acid_bg2_sl(:,:), c_water_bg2_sl = c_water_bg2_sl(:,:),                 &
      & c_ethanol_bg2_sl = c_ethanol_bg2_sl(:,:), c_nonsoluble_bg2_sl = c_nonsoluble_bg2_sl(:,:), &
      & c_humus_1_sl = c_humus_1_sl(:,:), c_humus_2_sl = c_humus_2_sl(:,:))
    END IF

  END SUBROUTINE rescale_carbon_upon_reference_area_change


  ! ====================================================================================================== !
  !
  !> Calculates all carbon tile variables (known to this function) from current canopy states
  !
  SUBROUTINE recalc_carbon_per_tile_vars(tile, options)

  ! -------------------------------------------------------------------------------------------------- !
  CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: tile
  TYPE(t_jsb_task_options),   INTENT(in)            :: options
  ! -------------------------------------------------------------------------------------------------- !

  dsl4jsb_Real2D_onChunk ::  c_green
  dsl4jsb_Real2D_onChunk ::  c_woods
  dsl4jsb_Real2D_onChunk ::  c_reserve
  dsl4jsb_Real2D_onChunk ::  c_crop_harvest

  dsl4jsb_Real2D_onChunk ::  c_acid_ag1
  dsl4jsb_Real2D_onChunk ::  c_water_ag1
  dsl4jsb_Real2D_onChunk ::  c_ethanol_ag1
  dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag1
  dsl4jsb_Real2D_onChunk ::  c_acid_bg1
  dsl4jsb_Real2D_onChunk ::  c_water_bg1
  dsl4jsb_Real2D_onChunk ::  c_ethanol_bg1
  dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg1

  dsl4jsb_Real2D_onChunk ::  c_acid_ag2
  dsl4jsb_Real2D_onChunk ::  c_water_ag2
  dsl4jsb_Real2D_onChunk ::  c_ethanol_ag2
  dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag2
  dsl4jsb_Real2D_onChunk ::  c_acid_bg2
  dsl4jsb_Real2D_onChunk ::  c_water_bg2
  dsl4jsb_Real2D_onChunk ::  c_ethanol_bg2
  dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg2

  dsl4jsb_Real2D_onChunk ::  c_humus_1
  dsl4jsb_Real2D_onChunk ::  c_humus_2

  dsl4jsb_Real2D_onChunk ::  cflux_dist_green_2_soil
  dsl4jsb_Real2D_onChunk ::  cflux_dist_woods_2_soil
  dsl4jsb_Real2D_onChunk ::  co2flux_fire_all_2_atm

  dsl4jsb_Real2D_onChunk ::  soil_respiration
  dsl4jsb_Real2D_onChunk ::  NPP_pot_yDayMean
  dsl4jsb_Real2D_onChunk ::  NPP_act_yDayMean
  dsl4jsb_Real2D_onChunk ::  NPP_flux_correction
  dsl4jsb_Real2D_onChunk ::  GPP_yDayMean
  dsl4jsb_Real2D_onChunk ::  cflux_c_greenwood_2_litter
  dsl4jsb_Real2D_onChunk ::  root_exudates
  dsl4jsb_Real2D_onChunk ::  cflux_c_green_2_herb

  dsl4jsb_Real3D_onChunk ::  c_acid_bg1_sl
  dsl4jsb_Real3D_onChunk ::  c_water_bg1_sl
  dsl4jsb_Real3D_onChunk ::  c_ethanol_bg1_sl
  dsl4jsb_Real3D_onChunk ::  c_nonsoluble_bg1_sl
  dsl4jsb_Real3D_onChunk ::  c_humus_1_sl

  dsl4jsb_Real3D_onChunk ::  c_acid_bg2_sl
  dsl4jsb_Real3D_onChunk ::  c_water_bg2_sl
  dsl4jsb_Real3D_onChunk ::  c_ethanol_bg2_sl
  dsl4jsb_Real3D_onChunk ::  c_nonsoluble_bg2_sl
  dsl4jsb_Real3D_onChunk ::  c_humus_2_sl

  dsl4jsb_Def_config(CARBON_)
  dsl4jsb_Def_memory(CARBON_)

  INTEGER                     :: ics, ice, iblk
  TYPE(t_jsb_model), POINTER  :: model
  ! -------------------------------------------------------------------------------------------------- !
  model => Get_model(tile%owner_model_id)
  dsl4jsb_Get_config(CARBON_)
  dsl4jsb_Get_memory(CARBON_)

  iblk  = options%iblk
  ics   = options%ics
  ice   = options%ice

  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest)

  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1)

  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2)

  dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_dist_green_2_soil)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_dist_woods_2_soil)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  co2flux_fire_all_2_atm)

  dsl4jsb_Get_var2D_onChunk(CARBON_,  soil_respiration)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_yDayMean)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_act_yDayMean)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_flux_correction)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_yDayMean)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_greenwood_2_litter)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  root_exudates)
  dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_green_2_herb)

  IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
    dsl4jsb_Get_var3D_onChunk(CARBON_,  c_acid_bg1_sl )
    dsl4jsb_Get_var3D_onChunk(CARBON_,  c_water_bg1_sl )
    dsl4jsb_Get_var3D_onChunk(CARBON_,  c_ethanol_bg1_sl )
    dsl4jsb_Get_var3D_onChunk(CARBON_,  c_nonsoluble_bg1_sl )
    dsl4jsb_Get_var3D_onChunk(CARBON_,  c_humus_1_sl )

    dsl4jsb_Get_var3D_onChunk(CARBON_,  c_acid_bg2_sl )
    dsl4jsb_Get_var3D_onChunk(CARBON_,  c_water_bg2_sl )
    dsl4jsb_Get_var3D_onChunk(CARBON_,  c_ethanol_bg2_sl )
    dsl4jsb_Get_var3D_onChunk(CARBON_,  c_nonsoluble_bg2_sl )
    dsl4jsb_Get_var3D_onChunk(CARBON_,  c_humus_2_sl )
  END IF


  CALL recalc_per_tile_vars(tile, options,                                    &
  & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:),      &
  & c_crop_harvest = c_crop_harvest(:),                                       &
  & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:),                 &
  & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
  & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:),                 &
  & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
  & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:),                 &
  & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
  & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:),                 &
  & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
  & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:), root_exudates = root_exudates(:),   &
  & soil_respiration = soil_respiration(:), NPP_flux_correction = NPP_flux_correction(:),   &
  & cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter(:),                             &
  & cflux_dist_green_2_soil = cflux_dist_green_2_soil(:),                                   &
  & cflux_dist_woods_2_soil = cflux_dist_woods_2_soil(:),                                   &
  & co2flux_fire_all_2_atm = co2flux_fire_all_2_atm(:),                                     &
  & cflux_c_green_2_herb = cflux_c_green_2_herb(:), NPP_act_yDayMean = NPP_act_yDayMean(:), &
  & NPP_pot_yDayMean = NPP_pot_yDayMean(:), GPP_yDayMean = GPP_yDayMean(:))

  IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
    CALL recalc_per_tile_vars(tile, options,                                                    &
    & c_acid_bg1_sl = c_acid_bg1_sl(:,:), c_water_bg1_sl = c_water_bg1_sl(:,:),                 &
    & c_ethanol_bg1_sl = c_ethanol_bg1_sl(:,:), c_nonsoluble_bg1_sl = c_nonsoluble_bg1_sl(:,:), &
    & c_acid_bg2_sl = c_acid_bg2_sl(:,:), c_water_bg2_sl = c_water_bg2_sl(:,:),                 &
    & c_ethanol_bg2_sl = c_ethanol_bg2_sl(:,:), c_nonsoluble_bg2_sl = c_nonsoluble_bg2_sl(:,:), &
    & c_humus_1_sl = c_humus_1_sl(:,:), c_humus_2_sl = c_humus_2_sl(:,:))
  END IF

END SUBROUTINE recalc_carbon_per_tile_vars


  ! ================================================================================================================================
  !>
  !> Calculates the per tile states from current canopy states
  !!
  !! @todo discuss how to make more general, e.g. get the info on which state variables to work on from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE recalc_per_tile_vars(tile, options,                          &
      & c_green , c_woods, c_reserve, c_crop_harvest,                     &
      & c_acid_ag1, c_water_ag1, c_ethanol_ag1, c_nonsoluble_ag1,         &
      & c_acid_bg1, c_water_bg1, c_ethanol_bg1, c_nonsoluble_bg1,         &
      & c_acid_ag2, c_water_ag2, c_ethanol_ag2, c_nonsoluble_ag2,         &
      & c_acid_bg2, c_water_bg2, c_ethanol_bg2, c_nonsoluble_bg2,         &
      & c_humus_1, c_humus_2,                                             &
      & c_decomp_humus_1_sum, c_decomp_humus_2_sum,                       &
      & c_into_humus_1_sum, c_into_humus_2_sum,                           &
      & soil_respiration,                                                 &
      & NPP_flux_correction, NPP_pot_yDayMean, NPP_act_yDayMean,          &
      & GPP_yDayMean, root_exudates, cflux_c_green_2_herb,                &
      & cflux_dist_green_2_soil, cflux_dist_woods_2_soil,                 &
      & co2flux_fire_all_2_atm,                                           &
      & cflux_c_greenwood_2_litter,                                       &
      & c_acid_bg1_sl, c_water_bg1_sl, c_ethanol_bg1_sl, c_nonsoluble_bg1_sl,         &
      & c_acid_bg2_sl, c_water_bg2_sl, c_ethanol_bg2_sl, c_nonsoluble_bg2_sl,         &
      & c_humus_1_sl, c_humus_2_sl)

    USE mo_carbon_process, ONLY: get_per_tile

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)            :: options
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_green(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_woods(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_reserve(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_crop_harvest(:)

    REAL(lp),                   INTENT(in), OPTIONAL  :: c_acid_ag1(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_water_ag1(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_ethanol_ag1(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_nonsoluble_ag1(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_acid_bg1(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_water_bg1(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_ethanol_bg1(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_nonsoluble_bg1(:)

    REAL(lp),                   INTENT(in), OPTIONAL  :: c_acid_ag2(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_water_ag2(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_ethanol_ag2(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_nonsoluble_ag2(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_acid_bg2(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_water_bg2(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_ethanol_bg2(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_nonsoluble_bg2(:)

    REAL(lp),                   INTENT(in), OPTIONAL  :: c_humus_1(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_humus_2(:)

    REAL(lp),                   INTENT(in), OPTIONAL  :: c_decomp_humus_1_sum(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_decomp_humus_2_sum(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_into_humus_1_sum(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_into_humus_2_sum(:)

    REAL(lp),                   INTENT(in), OPTIONAL  :: cflux_dist_green_2_soil(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: cflux_dist_woods_2_soil(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: co2flux_fire_all_2_atm(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: soil_respiration(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: NPP_pot_yDayMean(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: NPP_act_yDayMean(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: NPP_flux_correction(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: GPP_yDayMean(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: cflux_c_greenwood_2_litter(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: root_exudates(:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: cflux_c_green_2_herb(:)

    REAL(lp),                   INTENT(in), OPTIONAL  :: c_acid_bg1_sl(:,:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_water_bg1_sl(:,:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_ethanol_bg1_sl(:,:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_nonsoluble_bg1_sl(:,:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_acid_bg2_sl(:,:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_water_bg2_sl(:,:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_ethanol_bg2_sl(:,:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_nonsoluble_bg2_sl(:,:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_humus_1_sl(:,:)
    REAL(lp),                   INTENT(in), OPTIONAL  :: c_humus_2_sl(:,:)

    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(PHENO_)

    ! Local variables
    CHARACTER(len=*), PARAMETER     :: routine = modname//':recalc_per_tile_vars'
    INTEGER                         :: ics, ice, iblk
    INTEGER                         :: is, nsoil
    TYPE(t_jsb_model), POINTER      :: model

    dsl4jsb_Real2D_onChunk ::  c_green_ta
    dsl4jsb_Real2D_onChunk ::  c_reserve_ta
    dsl4jsb_Real2D_onChunk ::  c_woods_ta
    dsl4jsb_Real2D_onChunk ::  c_crop_harvest_ta

    dsl4jsb_Real2D_onChunk ::  c_acid_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_water_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_acid_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_water_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_humus_1_ta

    dsl4jsb_Real2D_onChunk ::  c_acid_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_water_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_acid_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_water_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_humus_2_ta

    dsl4jsb_Real2D_onChunk ::  c_decomp_humus_1_sum_ta
    dsl4jsb_Real2D_onChunk ::  c_decomp_humus_2_sum_ta
    dsl4jsb_Real2D_onChunk ::  c_into_humus_1_sum_ta
    dsl4jsb_Real2D_onChunk ::  c_into_humus_2_sum_ta

    dsl4jsb_Real2D_onChunk ::  cflux_dist_green_2_soil_ta
    dsl4jsb_Real2D_onChunk ::  cflux_dist_woods_2_soil_ta
    dsl4jsb_Real2D_onChunk ::  co2flux_fire_all_2_atm_ta
    dsl4jsb_Real2D_onChunk ::  soil_respiration_ta
    dsl4jsb_Real2D_onChunk ::  NPP_pot_yDayMean_ta
    dsl4jsb_Real2D_onChunk ::  NPP_act_yDayMean_ta
    dsl4jsb_Real2D_onChunk ::  NPP_flux_correction_ta
    dsl4jsb_Real2D_onChunk ::  GPP_yDayMean_ta
    dsl4jsb_Real2D_onChunk ::  cflux_c_greenwood_2_litter_ta
    dsl4jsb_Real2D_onChunk ::  root_exudates_ta
    dsl4jsb_Real2D_onChunk ::  cflux_c_green_2_herb_ta

    dsl4jsb_Real2D_onChunk ::  veg_fract_correction
    dsl4jsb_Real2D_onChunk ::  fract_fpc_max

    dsl4jsb_Real3D_onChunk ::  c_acid_bg1_sl_ta
    dsl4jsb_Real3D_onChunk ::  c_water_bg1_sl_ta
    dsl4jsb_Real3D_onChunk ::  c_ethanol_bg1_sl_ta
    dsl4jsb_Real3D_onChunk ::  c_nonsoluble_bg1_sl_ta
    dsl4jsb_Real3D_onChunk ::  c_humus_1_sl_ta

    dsl4jsb_Real3D_onChunk ::  c_acid_bg2_sl_ta
    dsl4jsb_Real3D_onChunk ::  c_water_bg2_sl_ta
    dsl4jsb_Real3D_onChunk ::  c_ethanol_bg2_sl_ta
    dsl4jsb_Real3D_onChunk ::  c_nonsoluble_bg2_sl_ta
    dsl4jsb_Real3D_onChunk ::  c_humus_2_sl_ta

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(CARBON_)
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(PHENO_)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green_ta )             ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve_ta )           ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods_ta )             ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest_ta)       ! opt out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1_ta)           ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1_ta)          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1_ta)        ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1_ta)     ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1_ta)           ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1_ta)          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1_ta)        ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1_ta)     ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1_ta)            ! opt out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2_ta)           ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2_ta)          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2_ta)        ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2_ta)     ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2_ta)           ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2_ta)          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2_ta)        ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2_ta)     ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2_ta)            ! opt out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_1_sum_ta)    ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_2_sum_ta)    ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_1_sum_ta)      ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_2_sum_ta)      ! opt out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_dist_green_2_soil_ta )   ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_dist_woods_2_soil_ta )   ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  co2flux_fire_all_2_atm_ta )    ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  soil_respiration_ta )          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_yDayMean_ta )          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_act_yDayMean_ta )          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_flux_correction_ta)        ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_yDayMean_ta )              ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_greenwood_2_litter_ta )! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  root_exudates_ta )             ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_green_2_herb_ta )      ! opt out

    dsl4jsb_Get_var2D_onChunk(PHENO_,   veg_fract_correction )       ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_,   fract_fpc_max )              ! in

    IF (dsl4jsb_Config(CARBON_)%vertical_yasso) THEN
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_acid_bg1_sl_ta )      ! opt out
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_water_bg1_sl_ta )      ! opt out
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_ethanol_bg1_sl_ta )      ! opt out
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_nonsoluble_bg1_sl_ta )      ! opt out
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_humus_1_sl_ta )      ! opt out

      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_acid_bg2_sl_ta )      ! opt out
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_water_bg2_sl_ta )      ! opt out
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_ethanol_bg2_sl_ta )      ! opt out
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_nonsoluble_bg2_sl_ta )      ! opt out
      dsl4jsb_Get_var3D_onChunk(CARBON_,  c_humus_2_sl_ta )      ! opt out
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    IF (PRESENT(c_green)) CALL get_per_tile(c_green_ta(:), c_green(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_woods)) CALL get_per_tile(c_woods_ta(:), c_woods(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_reserve)) CALL get_per_tile(c_reserve_ta(:), c_reserve(:), veg_fract_correction(:), fract_fpc_max(:))

    IF (PRESENT(c_crop_harvest)) CALL get_per_tile(c_crop_harvest_ta(:), &
      & c_crop_harvest(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_acid_ag1)) CALL get_per_tile(c_acid_ag1_ta(:), &
      & c_acid_ag1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_water_ag1)) CALL get_per_tile(c_water_ag1_ta(:), &
      & c_water_ag1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_ethanol_ag1)) CALL get_per_tile(c_ethanol_ag1_ta(:), &
      & c_ethanol_ag1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_nonsoluble_ag1)) CALL get_per_tile(c_nonsoluble_ag1_ta(:), &
      & c_nonsoluble_ag1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_acid_bg1)) CALL get_per_tile(c_acid_bg1_ta(:), &
      & c_acid_bg1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_water_bg1)) CALL get_per_tile(c_water_bg1_ta(:), &
      & c_water_bg1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_ethanol_bg1)) CALL get_per_tile(c_ethanol_bg1_ta(:), &
      & c_ethanol_bg1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_nonsoluble_bg1)) CALL get_per_tile(c_nonsoluble_bg1_ta(:), &
      & c_nonsoluble_bg1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_humus_1)) CALL get_per_tile(c_humus_1_ta(:), &
      & c_humus_1(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_acid_ag2)) CALL get_per_tile(c_acid_ag2_ta(:), &
      & c_acid_ag2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_water_ag2)) CALL get_per_tile(c_water_ag2_ta(:), &
      & c_water_ag2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_ethanol_ag2)) CALL get_per_tile(c_ethanol_ag2_ta(:), &
      & c_ethanol_ag2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_nonsoluble_ag2)) CALL get_per_tile(c_nonsoluble_ag2_ta(:), &
      & c_nonsoluble_ag2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_acid_bg2)) CALL get_per_tile(c_acid_bg2_ta(:), &
      & c_acid_bg2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_water_bg2)) CALL get_per_tile(c_water_bg2_ta(:), &
      & c_water_bg2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_ethanol_bg2)) CALL get_per_tile(c_ethanol_bg2_ta(:), &
      & c_ethanol_bg2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_nonsoluble_bg2)) CALL get_per_tile(c_nonsoluble_bg2_ta(:), &
      & c_nonsoluble_bg2(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_humus_2)) CALL get_per_tile(c_humus_2_ta(:), &
      & c_humus_2(:), veg_fract_correction(:), fract_fpc_max(:))

    IF (PRESENT(c_decomp_humus_1_sum)) CALL get_per_tile(c_decomp_humus_1_sum_ta(:), &
      & c_decomp_humus_1_sum(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_decomp_humus_2_sum)) CALL get_per_tile(c_decomp_humus_2_sum_ta(:), &
      & c_decomp_humus_2_sum(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_into_humus_1_sum)) CALL get_per_tile(c_into_humus_1_sum_ta(:), &
      & c_into_humus_1_sum(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_into_humus_2_sum)) CALL get_per_tile(c_into_humus_2_sum_ta(:), &
      & c_into_humus_2_sum(:), veg_fract_correction(:), fract_fpc_max(:))

    ! C fluxes
    IF (PRESENT(cflux_dist_green_2_soil)) CALL get_per_tile(cflux_dist_green_2_soil_ta(:), &
      & cflux_dist_green_2_soil(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(cflux_dist_woods_2_soil)) CALL get_per_tile(cflux_dist_woods_2_soil_ta(:), &
      & cflux_dist_woods_2_soil(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(co2flux_fire_all_2_atm))CALL get_per_tile(co2flux_fire_all_2_atm_ta(:), &
      & co2flux_fire_all_2_atm(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(NPP_pot_yDayMean)) CALL get_per_tile(NPP_pot_yDayMean_ta(:), &
      & NPP_pot_yDayMean(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(NPP_act_yDayMean)) CALL get_per_tile(NPP_act_yDayMean_ta(:), &
      & NPP_act_yDayMean(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(soil_respiration)) CALL get_per_tile(soil_respiration_ta(:), &
      & soil_respiration(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(NPP_flux_correction)) CALL get_per_tile(NPP_flux_correction_ta(:), &
      & NPP_flux_correction(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(GPP_yDayMean)) CALL get_per_tile(GPP_yDayMean_ta(:), &
      & GPP_yDayMean(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(cflux_c_greenwood_2_litter)) CALL get_per_tile(cflux_c_greenwood_2_litter_ta(:), &
      & cflux_c_greenwood_2_litter(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(root_exudates)) CALL get_per_tile(root_exudates_ta(:), &
      & root_exudates(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(cflux_c_green_2_herb)) CALL get_per_tile(cflux_c_green_2_herb_ta(:), &
      & cflux_c_green_2_herb(:), veg_fract_correction(:), fract_fpc_max(:))

    IF (PRESENT(c_acid_bg1_sl)) THEN
      nsoil = SIZE(c_acid_bg1_sl, DIM=2)
      DO is = 1,nsoil
        CALL get_per_tile(c_acid_bg1_sl_ta(:,is), &
        & c_acid_bg1_sl(:,is), veg_fract_correction(:), fract_fpc_max(:))
      END DO
    END IF
    IF (PRESENT(c_water_bg1_sl)) THEN
      nsoil = SIZE(c_water_bg1_sl, DIM=2)
      DO is = 1,nsoil
        CALL get_per_tile(c_water_bg1_sl_ta(:,is), &
        & c_water_bg1_sl(:,is), veg_fract_correction(:), fract_fpc_max(:))
      END DO
    END IF
    IF (PRESENT(c_ethanol_bg1_sl)) THEN
      nsoil = SIZE(c_ethanol_bg1_sl, DIM=2)
      DO is = 1,nsoil
        CALL get_per_tile(c_ethanol_bg1_sl_ta(:,is), &
        & c_ethanol_bg1_sl(:,is), veg_fract_correction(:), fract_fpc_max(:))
      END DO
    END IF
    IF (PRESENT(c_nonsoluble_bg1_sl)) THEN
      nsoil = SIZE(c_nonsoluble_bg1_sl, DIM=2)
      DO is = 1,nsoil
        CALL get_per_tile(c_nonsoluble_bg1_sl_ta(:,is), &
        & c_nonsoluble_bg1_sl(:,is), veg_fract_correction(:), fract_fpc_max(:))
      END DO
    END IF
    IF (PRESENT(c_acid_bg2_sl)) THEN
      nsoil = SIZE(c_acid_bg2_sl, DIM=2)
      DO is = 1,nsoil
        CALL get_per_tile(c_acid_bg2_sl_ta(:,is), &
        & c_acid_bg2_sl(:,is), veg_fract_correction(:), fract_fpc_max(:))
      END DO
    END IF
    IF (PRESENT(c_water_bg2_sl)) THEN
      nsoil = SIZE(c_water_bg2_sl, DIM=2)
      DO is = 1,nsoil
        CALL get_per_tile(c_water_bg2_sl_ta(:,is), &
        & c_water_bg2_sl(:,is), veg_fract_correction(:), fract_fpc_max(:))
      END DO
    END IF
    IF (PRESENT(c_ethanol_bg2_sl)) THEN
      nsoil = SIZE(c_ethanol_bg2_sl, DIM=2)
      DO is = 1,nsoil
        CALL get_per_tile(c_ethanol_bg2_sl_ta(:,is), &
        & c_ethanol_bg2_sl(:,is), veg_fract_correction(:), fract_fpc_max(:))
      END DO
    END IF
    IF (PRESENT(c_nonsoluble_bg2_sl)) THEN
      nsoil = SIZE(c_nonsoluble_bg2_sl, DIM=2)
      DO is = 1,nsoil
        CALL get_per_tile(c_nonsoluble_bg2_sl_ta(:,is), &
        & c_nonsoluble_bg2_sl(:,is), veg_fract_correction(:), fract_fpc_max(:))
      END DO
    END IF
    IF (PRESENT(c_humus_1_sl)) THEN
      nsoil = SIZE(c_humus_1_sl, DIM=2)
      DO is = 1,nsoil
        CALL get_per_tile(c_humus_1_sl_ta(:,is), &
        & c_humus_1_sl(:,is), veg_fract_correction(:), fract_fpc_max(:))
      END DO
    END IF
    IF (PRESENT(c_humus_2_sl)) THEN
      nsoil = SIZE(c_humus_2_sl, DIM=2)
      DO is = 1,nsoil
        CALL get_per_tile(c_humus_2_sl_ta(:,is), &
        & c_humus_2_sl(:,is), veg_fract_correction(:), fract_fpc_max(:))
      END DO
    END IF

  END SUBROUTINE recalc_per_tile_vars

  ! ================================================================================================================================
  !>
  !> calculates the current sum of all carbon state variables -> on ta! (per canopy vars are not aggregated!)
  !!
  !! @todo discuss how to make more general, e.g. get the info which are the c state variables from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE calculate_current_c_ta_state_sum(tile, options, current_c_state_sum_ta, &
    & c_sum_veg_ta, c_sum_litter_ag_ta, c_sum_litter_bg_ta, c_sum_humus_ta)

    !USE mo_util,                ONLY: real2string

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)         :: options
    REAL(lp),                   INTENT(out)        :: current_c_state_sum_ta(:)
    REAL(lp), OPTIONAL,         INTENT(out)        :: c_sum_veg_ta(:)
    REAL(lp), OPTIONAL,         INTENT(out)        :: c_sum_litter_ag_ta(:)
    REAL(lp), OPTIONAL,         INTENT(out)        :: c_sum_litter_bg_ta(:)
    REAL(lp), OPTIONAL,         INTENT(out)        :: c_sum_humus_ta(:)

    dsl4jsb_Def_memory(CARBON_)

    ! Local variables
    INTEGER :: nc, ic, ics, ice, iblk
    LOGICAL :: is_present_c_sum_veg_ta, is_present_c_sum_litter_ag_ta, is_present_c_sum_litter_bg_ta, is_present_c_sum_humus_ta

    CHARACTER(len=*), PARAMETER :: routine = modname//':calculate_current_c_ta_state_sum'

    dsl4jsb_Real2D_onChunk ::  c_green_ta
    dsl4jsb_Real2D_onChunk ::  c_reserve_ta
    dsl4jsb_Real2D_onChunk ::  c_woods_ta
    dsl4jsb_Real2D_onChunk ::  c_crop_harvest_ta

    dsl4jsb_Real2D_onChunk ::  c_acid_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_water_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag1_ta
    dsl4jsb_Real2D_onChunk ::  c_acid_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_water_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg1_ta
    dsl4jsb_Real2D_onChunk ::  c_humus_1_ta

    dsl4jsb_Real2D_onChunk ::  c_acid_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_water_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag2_ta
    dsl4jsb_Real2D_onChunk ::  c_acid_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_water_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg2_ta
    dsl4jsb_Real2D_onChunk ::  c_humus_2_ta

    is_present_c_sum_veg_ta       = PRESENT(c_sum_veg_ta)
    is_present_c_sum_litter_ag_ta = PRESENT(c_sum_litter_ag_ta)
    is_present_c_sum_litter_bg_ta = PRESENT(c_sum_litter_bg_ta)
    is_present_c_sum_humus_ta     = PRESENT(c_sum_humus_ta)

    dsl4jsb_Get_memory(CARBON_)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green_ta )               ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve_ta )             ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods_ta )               ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest_ta)         ! in

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1_ta)            ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1_ta)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1_ta)         ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1_ta)      ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1_ta)            ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1_ta)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1_ta)         ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1_ta)      ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1_ta)             ! in

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2_ta)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2_ta)          ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2_ta)        ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2_ta)     ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2_ta)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2_ta)          ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2_ta)        ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2_ta)     ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2_ta)            ! in

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      current_c_state_sum_ta(ic) = c_green_ta(ic) + c_woods_ta(ic) + c_reserve_ta(ic) + c_crop_harvest_ta(ic) + &
        & c_acid_ag1_ta(ic) + c_water_ag1_ta(ic) + c_ethanol_ag1_ta(ic) + c_nonsoluble_ag1_ta(ic) + &
        & c_acid_ag2_ta(ic) + c_water_ag2_ta(ic) + c_ethanol_ag2_ta(ic) + c_nonsoluble_ag2_ta(ic) + &
        & c_acid_bg1_ta(ic) + c_water_bg1_ta(ic) + c_ethanol_bg1_ta(ic) + c_nonsoluble_bg1_ta(ic) + &
        & c_acid_bg2_ta(ic) + c_water_bg2_ta(ic) + c_ethanol_bg2_ta(ic) + c_nonsoluble_bg2_ta(ic) + &
        & c_humus_1_ta(ic)  + c_humus_2_ta(ic)
    END DO
    !$ACC END LOOP

    IF (is_present_c_sum_veg_ta) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        c_sum_veg_ta(ic) = c_green_ta(ic) + c_woods_ta(ic) + c_reserve_ta(ic)
      END DO
      !$ACC END LOOP
    END IF

    IF (is_present_c_sum_litter_ag_ta) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        c_sum_litter_ag_ta(ic) = &
          & c_acid_ag1_ta(ic) + c_water_ag1_ta(ic) + c_ethanol_ag1_ta(ic) + c_nonsoluble_ag1_ta(ic) + &
          & c_acid_ag2_ta(ic) + c_water_ag2_ta(ic) + c_ethanol_ag2_ta(ic) + c_nonsoluble_ag2_ta(ic)
      END DO
      !$ACC END LOOP
    END IF

    IF (is_present_c_sum_litter_bg_ta) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        c_sum_litter_bg_ta(ic) = &
          & c_acid_bg1_ta(ic) + c_water_bg1_ta(ic) + c_ethanol_bg1_ta(ic) + c_nonsoluble_bg1_ta(ic) + &
          & c_acid_bg2_ta(ic) + c_water_bg2_ta(ic) + c_ethanol_bg2_ta(ic) + c_nonsoluble_bg2_ta(ic)
      END DO
      !$ACC END LOOP
    END IF

    IF (is_present_c_sum_humus_ta) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        c_sum_humus_ta(ic) = c_humus_1_ta(ic)  + c_humus_2_ta(ic)
      END DO
      !$ACC END LOOP
    END IF
    !$ACC END PARALLEL

  END SUBROUTINE calculate_current_c_ta_state_sum

  ! ====================================================================================================== !
  !
  !> Relocates content of (expected) own active to own passive vars that are part of lcc_relocations
  !
  SUBROUTINE carbon_transfer_from_active_to_passive_vars_onChunk(lcc_relocations, tile, i_tile, options)

    USE mo_util,                ONLY: one_of
    USE mo_jsb_lcc_class,       ONLY: t_jsb_lcc_proc, collect_matter_of_active_vars_onChunk
    USE mo_jsb_task_class,      ONLY: t_jsb_task_options
    USE mo_carbon_constants,    ONLY: fract_green_aboveGround, fract_wood_aboveGround,           &
      &                               nr_of_yasso_pools, carbon_required_passive_vars, coefficient_ind_of_passive_vars
    USE mo_carbon_process,      ONLY: add_litter_to_yasso_pool
    USE mo_jsb_varlist,         ONLY: VARNAME_LEN

    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(INOUT) :: lcc_relocations
        !< lcc structure of the calling lcc process for distributing matter from an active to a passive var
    CLASS(t_jsb_tile_abstract), INTENT(IN)    :: tile    !< tile for which the relocation is be conducted
    INTEGER,                    INTENT(IN)    :: i_tile  !< index of the tile in lcc structure
    TYPE(t_jsb_task_options),   INTENT(IN)    :: options !< run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':carbon_transfer_from_active_to_passive_vars_onChunk'

    TYPE(t_jsb_model), POINTER   :: model

    INTEGER :: ic, k, i_pool, i_passive, i_litter_coefficient, nc, ics, ice, iblk

    REAL(lp) :: litter_coefficient, fraction, this_fract
    REAL(lp) :: this_fract_green_aboveGround, this_fract_wood_aboveGround, fire_fract_wood_2_atmos, not_burned_AG_fract
    REAL(lp), DIMENSION(options%nc) :: non_woody_litter, woody_litter, this_litter

    dsl4jsb_Def_config(CARBON_)
    ! -------------------------------------------------------------------------------------------------- !

    !JN-TODO: is it less expensive to only do it for tiles with actual area changes?
    nc = options%nc
    ics = options%ics
    ice = options%ice
    iblk = options%iblk

    IF (debug_on() .AND. options%iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(CARBON_)


    !$ACC DATA CREATE(non_woody_litter, woody_litter, this_litter)

    ! Collect matter from potential active variables -- here unfortunately explicit
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      non_woody_litter(ic) = 0._lp
      woody_litter(ic) = 0._lp
    END DO
    !$ACC END PARALLEL LOOP

    CALL collect_matter_of_active_vars_onChunk( non_woody_litter, lcc_relocations, &
      & i_tile, iblk, ics, ice, [character(len=VARNAME_LEN) :: 'c_green', 'c_reserve' ])
    CALL collect_matter_of_active_vars_onChunk( woody_litter, lcc_relocations,     &
      & i_tile, iblk, ics, ice, [character(len=VARNAME_LEN) ::  'c_woods' ])

    ! Usually green, reserve and wood are distributed to ag and bg pools according to given fractions
    this_fract_green_aboveGround = fract_green_aboveGround
    this_fract_wood_aboveGround = fract_wood_aboveGround
    IF (lcc_relocations%name .EQ. 'flcc_lcc') THEN
      ! In case of fire there is no C transferred from green and reserve pools to ag litter pools, since it burned.
      this_fract_green_aboveGround = 0.0_lp
      ! Besides, part of the ag wood C is burned and has already been added to the atmosphere flux (vg: where????)
      fire_fract_wood_2_atmos = dsl4jsb_Config(CARBON_)%fire_fract_wood_2_atmos
      not_burned_AG_fract = (1.0_lp - fire_fract_wood_2_atmos) * fract_wood_aboveGround
      this_fract_wood_aboveGround = not_burned_AG_fract / (not_burned_AG_fract + (1.0_lp - fract_wood_aboveGround))
    END IF

    ! Distribution of the litter to the different yasso carbon pools
    DO i_pool = 1, nr_of_yasso_pools
      i_passive = one_of(carbon_required_passive_vars(i_pool), lcc_relocations%passive_vars_names)
      i_litter_coefficient = coefficient_ind_of_passive_vars(i_pool)

      IF (INDEX(carbon_required_passive_vars(i_pool), '1') > 0) THEN
        ! Pool name contains '1' => non woody C pool
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          this_litter(ic) = non_woody_litter(ic)
        END DO
        !$ACC END PARALLEL LOOP
        ! Fraction of litter allocated to above ground non-woody yasso pools (rather than below ground non-woody pools)
        this_fract = this_fract_green_aboveGround
        ! Fraction of litter allocated to the respective yasso pools (acid vs. water vs. ethanol vs. nonsoluble)
        litter_coefficient = dsl4jsb_Lctlib_param(LeafLit_coef(i_litter_coefficient))
      ELSE IF (INDEX(carbon_required_passive_vars(i_pool), '2') > 0) THEN
        ! Pool name contains '2' => woody C pool
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          this_litter(ic) = woody_litter(ic)
        END DO
        !$ACC END PARALLEL LOOP
        ! Fraction of litter allocated to above ground woody yasso pools (rather than below ground woody pools)
        this_fract = this_fract_wood_aboveGround
        ! Fraction of litter allocated to the respective yasso pools (acid vs. water vs. ethanol vs. nonsoluble)
        litter_coefficient = dsl4jsb_Lctlib_param(WoodLit_coef(i_litter_coefficient))
      ELSE
        CALL finish(TRIM(routine), 'Violation of assertion: ' // carbon_required_passive_vars(i_pool) &
          & // ' did not contain "1" or "2" which are key to identifiy woody vs non-woody yasso pools. Please check!')
      END IF

      IF (INDEX(carbon_required_passive_vars(i_pool), 'ag') > 0) THEN
        ! Pool name contains 'ag' => above ground C pool => fraction to above ground pool
        fraction = this_fract
      ELSE IF (INDEX(carbon_required_passive_vars(i_pool), 'bg') > 0) THEN
        ! Pool name contains 'bg' => below ground C pool => 1 - fraction to above ground pool
        fraction = 1.0_lp - this_fract
      ELSE IF (INDEX(carbon_required_passive_vars(i_pool), 'humus') > 0) THEN
        ! Pool name contains 'humus' => humus C pool
        fraction = 1.0_lp
      ELSE
        CALL finish(TRIM(routine), 'Violation of assertion: ' // carbon_required_passive_vars(i_pool) &
          & // ' did not contain "ag", "bg" or "humus" which are key to identifiy above vs below'     &
          & // ' ground yasso pools. Please check!')
      END IF

      !JN-TODO: better way? Reiner? DSL?
      DO ic = 1, nc
        k = ics + ic - 1
        CALL add_litter_to_yasso_pool( &
          & lcc_relocations%passive_vars(i_passive)%p%relocate_this(k,i_tile,iblk), &
          & this_litter(ic), fraction, litter_coefficient)
      END DO
    END DO
    !$ACC WAIT(acc_stream)
    !$ACC END DATA

  END SUBROUTINE carbon_transfer_from_active_to_passive_vars_onChunk


  ! ====================================================================================================== !
  !
  !> Rescales those carbon variables which need to be merged upon ageing (FAGE process), i.e.
  !> when the fraction with the maximum age of one age class is pushed into the next age class
  !
  SUBROUTINE rescale_carbon_vars_upon_ageing_induced_ref_area_change(tile, options, old_veg_fract_correction)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract),   INTENT(inout), TARGET :: tile
    TYPE(t_jsb_task_options),        INTENT(in)         :: options
    REAL(lp), DIMENSION(options%nc), INTENT(in)         :: old_veg_fract_correction
    ! -------------------------------------------------------------------------------------------------- !
    REAL(lp) :: veg_fract_corr_scaling
    INTEGER  :: iblk, ics, ice, nc, ic

    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(PHENO_)

    dsl4jsb_Real2D_onChunk ::  veg_fract_correction

    dsl4jsb_Real2D_onChunk ::  GPP_sum
    dsl4jsb_Real2D_onChunk ::  NPP_pot_sum
    dsl4jsb_Real2D_onChunk ::  LAI_yDayMean
    dsl4jsb_Real2D_onChunk ::  LAI_sum
    dsl4jsb_Real2D_onChunk ::  current_max_green
    dsl4jsb_Real2D_onChunk ::  veg_carbon_at_max_green
    ! -------------------------------------------------------------------------------------------------- !
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(PHENO_)

    dsl4jsb_Get_var2D_onChunk(CARBON_, GPP_sum)
    dsl4jsb_Get_var2D_onChunk(CARBON_, NPP_pot_sum)
    dsl4jsb_Get_var2D_onChunk(CARBON_, LAI_yDayMean)
    dsl4jsb_Get_var2D_onChunk(CARBON_, LAI_sum)
    dsl4jsb_Get_var2D_onChunk(CARBON_, current_max_green)
    dsl4jsb_Get_var2D_onChunk(CARBON_, veg_carbon_at_max_green)

    dsl4jsb_Get_var2D_onChunk(PHENO_, veg_fract_correction)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(veg_fract_corr_scaling)
    DO ic = 1, nc
      veg_fract_corr_scaling = old_veg_fract_correction(ic) / veg_fract_correction(ic)
      GPP_sum(ic)                 = GPP_sum(ic)                 * veg_fract_corr_scaling
      NPP_pot_sum(ic)             = NPP_pot_sum(ic)             * veg_fract_corr_scaling
      LAI_yDayMean(ic)            = LAI_yDayMean(ic)            * veg_fract_corr_scaling
      LAI_sum(ic)                 = LAI_sum(ic)                 * veg_fract_corr_scaling
      current_max_green(ic)       = current_max_green(ic)       * veg_fract_corr_scaling
      veg_carbon_at_max_green(ic) = veg_carbon_at_max_green(ic) * veg_fract_corr_scaling
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE rescale_carbon_vars_upon_ageing_induced_ref_area_change


  ! ====================================================================================================== !
  !
  !> Induces merging for those carbon variables which need to be merged upon ageing (FAGE process)
  !
  SUBROUTINE merge_carbon_vars_upon_ageing(target, source, options, moved_area)

    USE mo_fage_process,     ONLY : weighted_avg_per_canopy_var_upon_area_movement
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract),      INTENT(inout) :: target     !< target age class
    CLASS(t_jsb_tile_abstract),      INTENT(inout) :: source     !< source age class
    TYPE(t_jsb_task_options),        INTENT(in)    :: options    !< Additional run-time parameters
    REAL(lp), DIMENSION(options%nc), INTENT(in)    :: moved_area !< area moved from source to target ac
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER  :: iblk, ics, ice

    dsl4jsb_Def_memory_tile(PHENO_, target)
    dsl4jsb_Def_memory_tile(PHENO_, source)

    dsl4jsb_Def_memory_tile(CARBON_, target)
    dsl4jsb_Def_memory_tile(CARBON_, source)

    dsl4jsb_Real2D_onChunk :: veg_fract_correction_source
    dsl4jsb_Real2D_onChunk :: veg_fract_correction_target

    dsl4jsb_Real2D_onChunk :: GPP_sum_source
    dsl4jsb_Real2D_onChunk :: GPP_sum_target
    dsl4jsb_Real2D_onChunk :: NPP_pot_sum_source
    dsl4jsb_Real2D_onChunk :: NPP_pot_sum_target
    dsl4jsb_Real2D_onChunk :: LAI_sum_source
    dsl4jsb_Real2D_onChunk :: LAI_sum_target
    dsl4jsb_Real2D_onChunk :: LAI_yDayMean_source
    dsl4jsb_Real2D_onChunk :: LAI_yDayMean_target
    dsl4jsb_Real2D_onChunk :: current_max_green_source
    dsl4jsb_Real2D_onChunk :: current_max_green_target
    dsl4jsb_Real2D_onChunk :: veg_carbon_at_max_green_source
    dsl4jsb_Real2D_onChunk :: veg_carbon_at_max_green_target

    REAL(lp), DIMENSION(options%nc) :: target_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':merge_carbon_vars_upon_ageing'
    ! -------------------------------------------------------------------------------------------------- !
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For target '//TRIM(target%name)//' ...')

    dsl4jsb_Get_memory_tile(PHENO_, source)
    dsl4jsb_Get_memory_tile(PHENO_, target)

    dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, veg_fract_correction, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, veg_fract_correction, target)

    dsl4jsb_Get_memory_tile(CARBON_, source)
    dsl4jsb_Get_memory_tile(CARBON_, target)

    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, GPP_sum, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, GPP_sum, target)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, NPP_pot_sum, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, NPP_pot_sum, target)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, LAI_yDayMean, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, LAI_yDayMean, target)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, LAI_sum, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, LAI_sum, target)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, current_max_green, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, current_max_green, target)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, veg_carbon_at_max_green, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, veg_carbon_at_max_green, target)

    CALL target%Get_fraction(ics, ice, iblk, fract=target_fract(:))

    CALL weighted_avg_per_canopy_var_upon_area_movement(target_fract(:), moved_area(:),   &
      &  veg_fract_correction_source(:), veg_fract_correction_target(:),                  &
      &  GPP_sum_source(:), GPP_sum_target(:))
    CALL weighted_avg_per_canopy_var_upon_area_movement(target_fract(:), moved_area(:),   &
      &  veg_fract_correction_source(:), veg_fract_correction_target(:),                  &
      &  NPP_pot_sum_source(:), NPP_pot_sum_target(:))
    CALL weighted_avg_per_canopy_var_upon_area_movement(target_fract(:), moved_area(:),   &
      &  veg_fract_correction_source(:), veg_fract_correction_target(:),                  &
      &  LAI_sum_source(:), LAI_sum_target(:))
    CALL weighted_avg_per_canopy_var_upon_area_movement(target_fract(:), moved_area(:),   &
      &  veg_fract_correction_source(:), veg_fract_correction_target(:),                  &
      &  LAI_yDayMean_source(:), LAI_yDayMean_target(:))
    CALL weighted_avg_per_canopy_var_upon_area_movement(target_fract(:), moved_area(:),   &
      &  veg_fract_correction_source(:), veg_fract_correction_target(:),                  &
      &  current_max_green_source(:), current_max_green_target(:))
    CALL weighted_avg_per_canopy_var_upon_area_movement(target_fract(:), moved_area(:),   &
      &  veg_fract_correction_source(:), veg_fract_correction_target(:),                  &
      &  veg_carbon_at_max_green_source(:), veg_carbon_at_max_green_target(:))

  END SUBROUTINE merge_carbon_vars_upon_ageing


  ! ================================================================================================================================
  !>
  !> conducts a carbon conservation test on this tile's current states and given old sum + fluxes
  !!        -> on ta! (per canopy vars are not aggregated!)
  !!
  !! @todo discuss how to make more general, e.g. get the info which are the c state variables from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE check_carbon_conservation(tile, options, old_c_state_sum_ta, cflux_ta, Cconserve)

    !USE mo_util,                ONLY: real2string

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)         :: options
    REAL(lp),                   INTENT(in)         :: old_c_state_sum_ta(:)
    REAL(lp),                   INTENT(in)         :: cflux_ta(:)
        ! Note: cflux is not really a flux, but already multiplied with time unit!
    REAL(lp),                   INTENT(inout)      :: Cconserve(:)

    ! Local variables
    INTEGER :: ic, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':check_carbon_conservation'

    nc = options%nc

    CALL calculate_current_c_ta_state_sum(tile, options, Cconserve)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      Cconserve(ic) = Cconserve(ic) - old_c_state_sum_ta(ic) - cflux_ta(ic)
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE check_carbon_conservation

  ! ================================================================================================================================
  !>
  !> conducts a box carbon conservation test for the last day, called each new day -- (jsbach_start_timestep)
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE yday_carbon_conservation_test(tile)

    USE mo_jsb_grid,           ONLY: Get_grid

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile

    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(L2A_)

    dsl4jsb_Real2D_onChunk ::  carbon_conservation_test, yday_c_state_sum
    dsl4jsb_Real2D_onChunk ::  co2flux_npp_2_atm_yday_ta
    dsl4jsb_Real2D_onChunk ::  co2flux_soilresp_2_atm_ta
    dsl4jsb_Real2D_onChunk ::  co2flux_herb_2_atm_ta
    dsl4jsb_Real2D_onChunk ::  co2flux_fire_all_2_atm_ta

    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid
    TYPE(t_jsb_task_options)        :: options

    REAL(lp), ALLOCATABLE :: C_flux_yday(:)
    INTEGER               :: ic, ics, ice, iblk, nc

    dsl4jsb_Get_memory(L2A_)
    dsl4jsb_Get_memory(CARBON_)

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    ! Note: options is defined here as a local variable. It is needed in calculate_current_c_ta_state_sum
    options%ics = 1
    options%ice = grid%nproma

    ics = options%ics
    ice = options%ice

    ALLOCATE(C_flux_yday(grid%nproma))
    !$ACC DATA CREATE(C_flux_yday)

    DO iblk = 1, grid%nblks

      options%iblk = iblk
      IF (iblk == grid%nblks) THEN
        options%ice = grid%npromz
        ice = grid%npromz
        ! TODO: ics and ice should be retrieved for each iblk (are not necessarily constant for each block in ICON)
      END IF
      options%nc = ice - ics + 1
      nc = options%nc

      dsl4jsb_Get_var2D_onChunk(CARBON_, co2flux_npp_2_atm_yday_ta)   ! in
      dsl4jsb_Get_var2D_onChunk(CARBON_, co2flux_soilresp_2_atm_ta)   ! in
      dsl4jsb_Get_var2D_onChunk(CARBON_, co2flux_herb_2_atm_ta)       ! in
      dsl4jsb_Get_var2D_onChunk(CARBON_, co2flux_fire_all_2_atm_ta)   ! in

      dsl4jsb_Get_var2D_onChunk(L2A_,  carbon_conservation_test )     ! inout
      dsl4jsb_Get_var2D_onChunk(L2A_,  yday_c_state_sum )             ! inout

      ! Do not test as long as carbon_conservation_test still carries the initialisation value
      IF (.NOT. ALL(carbon_conservation_test .EQ. -999.0_lp)) THEN

        ! CO2 flux due to npp, soil respiration and herbivory
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          C_flux_yday(ic) = &
            & - ((  co2flux_npp_2_atm_yday_ta(ic)  &
            &     + co2flux_soilresp_2_atm_ta(ic)  &
            &     + co2flux_herb_2_atm_ta(ic)      &
            &   ) * sec_per_day / molarMassCO2_kg)
        END DO
        !$ACC END PARALLEL LOOP

        ! CO2 flux due to fire
        IF (model%Is_process_enabled(DISTURB_)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
          DO ic = 1, nc
            C_flux_yday(ic) =  C_flux_yday(ic) &
              & - co2flux_fire_all_2_atm_ta(ic) * sec_per_day / molarMassCO2_kg
          END DO
          !$ACC END PARALLEL LOOP
        END IF

        CALL calculate_current_c_ta_state_sum(tile, options, carbon_conservation_test)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          carbon_conservation_test(ic) = carbon_conservation_test(ic) - yday_c_state_sum(ic)  - C_flux_yday(ic)
        END DO
        !$ACC END PARALLEL LOOP

      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic = 1, nc
          carbon_conservation_test(ic) = 0._lp
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      CALL calculate_current_c_ta_state_sum(tile, options, yday_c_state_sum(:))

    END DO

    !$ACC WAIT(acc_stream)
    !$ACC END DATA

    DEALLOCATE(C_flux_yday)

  END SUBROUTINE yday_carbon_conservation_test

  ! ================================================================================================================================
  !>
  !> calculations of diagnostic global sums
  !!        called from jsbach_finish_timestep, after the loop over the nproma blocks.
  !! @param[in,out] tile    Tile for which routine is executed.
  !
  SUBROUTINE global_carbon_diagnostics(tile, land_area_scaling)

    USE mo_carbon_constants,        ONLY: molarMassC_kg, sec_per_year
    USE mo_util,                    ONLY: scaled_global_sum_2D
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile                   !< Current tile
    REAL(lp),                   INTENT(in) :: land_area_scaling(:,:) !< Land area scaling factor (accounting for halo cells)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(CARBON_)

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_carbon_diagnostics'

    dsl4jsb_Real2D_onDomain :: C_sum_veg_ta
    dsl4jsb_Real2D_onDomain :: C_sum_litter_ag_ta
    dsl4jsb_Real2D_onDomain :: C_sum_litter_bg_ta
    dsl4jsb_Real2D_onDomain :: C_sum_humus_ta
    dsl4jsb_Real2D_onDomain :: C_sum_natural_ta
    dsl4jsb_Real2D_onDomain :: NPP_act_yDayMean_ta
    dsl4jsb_Real2D_onDomain :: GPP_yDayMean_ta
    dsl4jsb_Real2D_onDomain :: soil_respiration_ta

    REAL(lp), POINTER       :: C_sum_veg_gsum(:)
    REAL(lp), POINTER       :: C_sum_litter_ag_gsum(:)
    REAL(lp), POINTER       :: C_sum_litter_bg_gsum(:)
    REAL(lp), POINTER       :: C_sum_humus_gsum(:)
    REAL(lp), POINTER       :: C_sum_natural_gsum(:)
    REAL(lp), POINTER       :: NPP_act_yDayMean_gsum(:)
    REAL(lp), POINTER       :: GPP_yDayMean_gsum(:)
    REAL(lp), POINTER       :: soil_respiration_gsum(:)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_var2D_onDomain(CARBON_,  C_sum_veg_ta)                ! in
    dsl4jsb_Get_var2D_onDomain(CARBON_,  C_sum_litter_ag_ta)          ! in
    dsl4jsb_Get_var2D_onDomain(CARBON_,  C_sum_litter_bg_ta)          ! in
    dsl4jsb_Get_var2D_onDomain(CARBON_,  C_sum_humus_ta)              ! in
    dsl4jsb_Get_var2D_onDomain(CARBON_,  C_sum_natural_ta)            ! in
    dsl4jsb_Get_var2D_onDomain(CARBON_,  NPP_act_yDayMean_ta)         ! in
    dsl4jsb_Get_var2D_onDomain(CARBON_,  GPP_yDayMean_ta)             ! in
    dsl4jsb_Get_var2D_onDomain(CARBON_,  soil_respiration_ta)         ! in

    C_sum_veg_gsum        => CARBON__mem%C_sum_veg_gsum%ptr(:)        ! out
    C_sum_litter_ag_gsum  => CARBON__mem%C_sum_litter_ag_gsum%ptr(:)  ! out
    C_sum_litter_bg_gsum  => CARBON__mem%C_sum_litter_bg_gsum%ptr(:)  ! out
    C_sum_humus_gsum      => CARBON__mem%C_sum_humus_gsum%ptr(:)      ! out
    C_sum_natural_gsum    => CARBON__mem%C_sum_natural_gsum%ptr(:)    ! out
    NPP_act_yDayMean_gsum => CARBON__mem%NPP_act_yDayMean_gsum%ptr(:) ! out
    GPP_yDayMean_gsum     => CARBON__mem%GPP_yDayMean_gsum%ptr(:)     ! out
    soil_respiration_gsum => CARBON__mem%soil_respiration_gsum%ptr(:) ! out

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')
    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    ! Calculate global carbon inventories, if requested for output
    !  => Conversion from [mol(C)/m^2] to [GtC]
    !     1 mol C = molarMassC_kg kg C   => 1 mol C = molarMassC_kg * e-12 Gt C
    IF (CARBON__mem%C_sum_veg_gsum%is_in_output) THEN
      CALL scaled_global_sum_2D(C_sum_veg_gsum, C_sum_veg_ta, land_area_scaling, molarMassC_kg * 1.e-12_lp, acc_stream)
      !$ACC UPDATE ASYNC(acc_stream) DEVICE(CARBON__mem%C_sum_veg_gsum%ptr(:))
    END IF
    IF (CARBON__mem%C_sum_litter_ag_gsum%is_in_output) THEN
      CALL scaled_global_sum_2D(C_sum_litter_ag_gsum, C_sum_litter_ag_ta, land_area_scaling, molarMassC_kg * 1.e-12_lp, acc_stream)
      !$ACC UPDATE ASYNC(acc_stream) DEVICE(CARBON__mem%C_sum_litter_ag_gsum%ptr(:))
    END IF
    IF (CARBON__mem%C_sum_litter_bg_gsum%is_in_output) THEN
      CALL scaled_global_sum_2D(C_sum_litter_bg_gsum, C_sum_litter_bg_ta, land_area_scaling, molarMassC_kg * 1.e-12_lp, acc_stream)
      !$ACC UPDATE ASYNC(acc_stream) DEVICE(CARBON__mem%C_sum_litter_bg_gsum%ptr(:))
    END IF
    IF (CARBON__mem%C_sum_humus_gsum%is_in_output) THEN
      CALL scaled_global_sum_2D(C_sum_humus_gsum, C_sum_humus_ta, land_area_scaling, molarMassC_kg * 1.e-12_lp, acc_stream)
      !$ACC UPDATE ASYNC(acc_stream) DEVICE(CARBON__mem%C_sum_humus_gsum%ptr(:))
    END IF
    IF (CARBON__mem%C_sum_natural_gsum%is_in_output) THEN
      CALL scaled_global_sum_2D(C_sum_natural_gsum, C_sum_natural_ta, land_area_scaling, molarMassC_kg * 1.e-12_lp, acc_stream)
      !$ACC UPDATE ASYNC(acc_stream) DEVICE(CARBON__mem%C_sum_natural_gsum%ptr(:))
    END IF

    ! ... fluxes * sec_per_year
    IF (CARBON__mem%NPP_act_yDayMean_gsum%is_in_output) THEN
      CALL scaled_global_sum_2D(NPP_act_yDayMean_gsum, NPP_act_yDayMean_ta, land_area_scaling, &
          & molarMassC_kg * 1.e-12_lp * sec_per_year, acc_stream)
      !$ACC UPDATE ASYNC(acc_stream) DEVICE(CARBON__mem%NPP_act_yDayMean_gsum%ptr(:))
    END IF
    IF (CARBON__mem%GPP_yDayMean_gsum%is_in_output) THEN
      CALL scaled_global_sum_2D(GPP_yDayMean_gsum, GPP_yDayMean_ta, land_area_scaling, &
          & molarMassC_kg * 1.e-12_lp * sec_per_year, acc_stream)
      !$ACC UPDATE ASYNC(acc_stream) DEVICE(CARBON__mem%GPP_yDayMean_gsum%ptr(:))
    END IF
    IF (CARBON__mem%soil_respiration_gsum%is_in_output) THEN
      CALL scaled_global_sum_2D(soil_respiration_gsum, soil_respiration_ta, land_area_scaling, &
          & molarMassC_kg * 1.e-12_lp * sec_per_year, acc_stream)
      !$ACC UPDATE ASYNC(acc_stream) DEVICE(CARBON__mem%soil_respiration_gsum%ptr(:))
    END IF

  END SUBROUTINE global_carbon_diagnostics

#endif
END MODULE mo_carbon_interface
