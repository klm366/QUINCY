!>
!! @brief Contains the interfaces to the carbon process
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-04-05)
!! Julia Nabel (2018): l_forestRegrowth flag -> allow for changing MaxLai; forest age cohorts
!!      ++ routines to check carbon conservation 
!!      ++ task to apply c consequences of area changes among forest age cohorts 
!!
MODULE mo_carbon_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, message_text, finish

  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module (Get_carbon_memory and Get_carbon_config)
  dsl4jsb_Use_processes CARBON_, ASSIMI_, HYDRO_, PHENO_, RAD_, A2L_, DISTURB_, FAGE_

  ! Use of process configurations (t_carbon_config)
  dsl4jsb_Use_config(ASSIMI_)
  dsl4jsb_Use_config(CARBON_)
  dsl4jsb_Use_config(DISTURB_)
  dsl4jsb_Use_config(PHENO_)
  dsl4jsb_Use_config(FAGE_)

  ! Use of process memories (t_carbon_memory)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(CARBON_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_carbon_tasks ,                                                       &
               & recalc_per_tile_vars,                                                     &
               & recalc_delta_carbon_rp_of_source_and_sink_upon_area_shift,                &
               & recalc_delta_carbon_rp_upon_disturbance,                                  &
               & apply_delta_carbon_to_carbon_states,                                   &
               & rescale_carbon_upon_reference_area_change, calculate_current_c_state_sum, &
               & check_carbon_conservation_for_operation

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_interface'

  !> Type definition for C_NPP_pot_allocation task
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_C_NPP_pot_allocation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_C_NPP_pot_allocation   !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_C_NPP_pot_allocation !< Aggregates computed task variables
  END TYPE   tsk_C_NPP_pot_allocation

  !> Constructor interface for C_NPP_pot_allocation task
  INTERFACE   tsk_C_NPP_pot_allocation
    PROCEDURE Create_task_C_NPP_pot_allocation        !< Constructor function for task
  END INTERFACE   tsk_C_NPP_pot_allocation

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_apply_age_cohort_changes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_apply_age_cohort_changes    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_apply_age_cohort_changes !< Aggregates computed task variables
  END TYPE   tsk_apply_age_cohort_changes

  INTERFACE   tsk_apply_age_cohort_changes
    PROCEDURE Create_task_apply_age_cohort_changes 
  END INTERFACE   tsk_apply_age_cohort_changes

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

    return_ptr%name = 'C_NPP_pot_allocation'
    return_ptr%process_id = CARBON_
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_C_NPP_pot_allocation
  
  !> Constructor for apply_age_cohort_changes task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "apply_age_cohort_changes"
  !!
  FUNCTION Create_task_apply_age_cohort_changes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_apply_age_cohort_changes::return_ptr)

    return_ptr%name = 'apply_age_cohort_changes'
    return_ptr%process_id = CARBON_
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_apply_age_cohort_changes

  ! =======================================================================================================
  !> Register tasks for carbon process
  !!
  !! @param[in,out] this      Instance of carbon process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_carbon_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_C_NPP_pot_allocation(model_id))
    CALL this%Register_task(tsk_apply_age_cohort_changes(model_id))

  END SUBROUTINE Register_carbon_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation to calculate the allocation of NPP carbon to the plant and soil carbon pools
  !!        R: Corresponds to update_cbalance_bethy in JSBACH3.
  !! 
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_C_NPP_pot_allocation(tile, options)

    ! Use declarations
    USE mo_jsb_grid,          ONLY: Get_grid
    USE mo_carbon_process,    ONLY: calc_Cpools
    USE mo_jsb_time,          ONLY: is_newday, is_newyear, time_steps_per_day, is_time_experiment_start
    USE mo_carbon_constants,  ONLY: molarMassCO2_kg, sec_per_day, days_per_year

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Declare local variables
    TYPE(t_jsb_model), POINTER   :: model
    TYPE(t_jsb_grid),  POINTER   :: grid
    LOGICAL                      :: lstart
    INTEGER                      :: ics, ice, iblk, nc
    REAL(wp)                     :: dtime
    INTEGER                      :: i, itile
    REAL(wp), allocatable        :: precip_total(:)
    REAL(wp), DIMENSION(options%nc) :: MaxLai
    REAL(wp), DIMENSION(options%nc) :: old_c_state_sum, current_fluxes
    CHARACTER(len=*), PARAMETER  :: routine = TRIM(modname)//':update_C_NPP_pot_allocation'

    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_memory(CARBON_)
    !dsl4jsb_Def_config(A2L_)
    dsl4jsb_Def_memory(A2L_)

    ! Declare pointers to variables in memory
    ! R: former cbalance_type variables
    dsl4jsb_Real2D_onChunk ::  pseudo_temp
    dsl4jsb_Real2D_onChunk ::  N_pseudo_temp
    dsl4jsb_Real2D_onChunk ::  F_pseudo_temp
    dsl4jsb_Real2D_onChunk ::  pseudo_temp_yDay
    dsl4jsb_Real2D_onChunk ::  pseudo_precip
    dsl4jsb_Real2D_onChunk ::  pseudo_precip_yDay
    dsl4jsb_Real2D_onChunk ::  N_pseudo_precip
    dsl4jsb_Real2D_onChunk ::  F_pseudo_precip

    dsl4jsb_Real2D_onChunk ::  topSoilTemp           ! was not cbalance_type
    dsl4jsb_Real2D_onChunk ::  alpha                 ! was not cbalance_type. Relative soil moisture (used for
                                                     ! heterotrophic respiration)
    dsl4jsb_Real2D_onChunk :: cconservation_calcCpools

    dsl4jsb_Real2D_onChunk ::  cflux_c_greenwood_2_litter           ! Carbon flux from the veget. to the litter pools
                                                     ! [mol(C)/m^2(canopy) s]
    dsl4jsb_Real2D_onChunk ::  c_green           ! C-pool for leaves, fine roots, vegetative organs and
                                                     ! other green (living) parts of vegetation [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_reserve         ! C-pool for carbohydrate reserve (sugars, starches) that
                                                     ! allows plants to survive bad times[mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_woods           ! C-pool for stems, thick roots and other (dead) structural
                                                     !  material of living plants [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_crop_harvest    ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]
    !SIZE CLASS 1
    dsl4jsb_Real2D_onChunk ::  c_acid_ag1        ! Yasso above ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_ag1       ! Yasso above ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag1     ! Yasso above ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag1  ! Yasso above ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_acid_bg1        ! Yasso below ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_bg1       ! Yasso below ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg1     ! Yasso below ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg1  ! Yasso below ground litter-pool for non-soluble litter
                                                      ! [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_humus_1         ! Yasso below ground litter-pool for slow C compartment
                                                      !  [mol(C)^2/m(canopy)]
    !SIZE CLASS 2
    dsl4jsb_Real2D_onChunk ::  c_acid_ag2        ! Yasso above ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_ag2       ! Yasso above ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag2     ! Yasso above ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag2  ! Yasso above ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_acid_bg2        ! Yasso below ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_bg2       ! Yasso below ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg2     ! Yasso below ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg2  ! Yasso below ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_humus_2         ! Yasso below ground litter-pool for slow C compartment
                                                      !  [mol(C)^2/m(canopy)]
    !
    !
    dsl4jsb_Real2D_onChunk ::  gross_assimilation_ca  ! Gross assimilation of all canopy layers together
    dsl4jsb_Real2D_onChunk ::  NPP_pot_rate_ca        ! The instantaneous (potential) NPP rate [mol(C)/m^2(canopy) s]
    !dsl4jsb_Real2D_onChunk ::  NPP_pot_Rate_acc      ! averaged NPP rate over one day [mol(C)/m^2(canopy) s]
    dsl4jsb_Real2D_onChunk ::  LAI_sum                ! used to accumulate LAI over a day. Sum of LAI since Midnight.
    dsl4jsb_Real2D_onChunk ::  NPP_pot_sum            ! used to accumulated NPP-Rate over a day
    dsl4jsb_Real2D_onChunk ::  GPP_sum                ! used to accumulated GPP-Rate over a day
    dsl4jsb_Real2D_onChunk ::  topSoilTemp_sum        ! used to accumulated topSoilTemp over a day
    dsl4jsb_Real2D_onChunk ::  alpha_sum              ! used to accumulated water stress coefficient alpha over a day.
    dsl4jsb_Real2D_onChunk ::  LAI_yDayMean           ! previous days mean
    dsl4jsb_Real2D_onChunk ::  LAI_yyDayMean          ! previous previous days mean
    dsl4jsb_Real2D_onChunk ::  NPP_pot_yDayMean       ! mean value of NPP-Rate yesterday (from NPP_pot_sum())
                                                      ! [mol(CO2)/(m^2(canopy) s)]
                                                      ! = cbalance%NPP_pot_sum(...)/time_steps_per_day
    dsl4jsb_Real2D_onChunk ::  NPP_act_yDayMean       ! mean value of actual NPP-Rate yesterday, i.e. after N-limitation
                                                      ! Actual NPP after N-limitation and excess carbon drop.
                                                      ! [mol(CO2)/(m^2(canopy) s)]
    dsl4jsb_Real2D_onChunk ::  GPP_yDayMean 
    dsl4jsb_Real2D_onChunk ::  topSoilTemp_yDayMean
    dsl4jsb_Real2D_onChunk ::  alpha_yDayMean

    dsl4jsb_Real2D_onChunk ::  soil_respiration       ! mean daily rate of heterotrophic (soil) respiration
                                                      ! [mol(CO2)/m^2(ground)].
                                                      ! Without N limitation!
    dsl4jsb_Real2D_onChunk ::  NPP_flux_correction    ! Daily updated flux correction from yesterdays carbon balance 
                                                      ! [mol(CO2)/m^2(canopy) s]
                                                      ! Amount by which the NPP rate entering the routine has to be corrected.
                                                      ! This correction arises either because otherwise the reserve pool would
                                                      ! get negative (positive correction), or the wood pool would exceed its
                                                      ! maximum value (negative correction).
    dsl4jsb_Real2D_onChunk ::  excess_NPP             ! That part of NPP that because of structural limits could not be
                                                      ! allocated in carbon pools [mol(CO2)/m^2(canopy) s]



    dsl4jsb_Real2D_onChunk ::  root_exudates          ! Total root exudates entering to the litter green pools
                                                      ! [mol(C)/m^2(canopy) s]

    dsl4jsb_Real2D_onChunk ::  annual_avg_NPP 
    dsl4jsb_Real2D_onChunk ::  biennial_avg_NPP 
    dsl4jsb_Real2D_onChunk ::  decennial_avg_NPP 

    dsl4jsb_Real2D_onChunk ::  cflux_npp_2_atm_yday_ta  ! day CO2 flux from actual NPP, required for cconservation test
    dsl4jsb_Real2D_onChunk ::  cflux_c_green_2_herb        !
    dsl4jsb_Real2D_onChunk ::  cflux_herb_2_littergreen     !
    dsl4jsb_Real2D_onChunk ::  cflux_herb_2_atm  !
    dsl4jsb_Real2D_onChunk ::  cflux_npp_2_atm_ta       ! grid cell averages of net CO2 fluxes between
    dsl4jsb_Real2D_onChunk ::  cflux_soilresp_2_atm_ta  ! .. biosphere (due to NPP, soil respiration and
    dsl4jsb_Real2D_onChunk ::  cflux_herb_2_atm_ta ! .. grazing) and atmosphere [kg(CO2)/(m^2(ground) s)] 

    ! This is on tile area:
    dsl4jsb_Real2D_onChunk ::  soil_respiration_ta      ! out
    dsl4jsb_Real2D_onChunk ::  NPP_pot_yDayMean_ta      ! out
    dsl4jsb_Real2D_onChunk ::  NPP_act_yDayMean_ta      ! out
    dsl4jsb_Real2D_onChunk ::  NPP_flux_correction_ta   ! out
    dsl4jsb_Real2D_onChunk ::  GPP_yDayMean_ta          ! out
!    dsl4jsb_Real2D_onChunk ::  cflux_c_greenwood_2_litter_ta           ! out
    dsl4jsb_Real2D_onChunk ::  root_exudates_ta         ! out
    dsl4jsb_Real2D_onChunk ::  cflux_c_green_2_herb_ta       ! out

    ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
    !dsl4jsb_Real2D_onChunk :: fract_litter_wood_new
    !
    dsl4jsb_Real2D_onChunk :: veg_fract_correction
    dsl4jsb_Real2D_onChunk :: LAI 
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: t_air                       ! Atmosphere temperature (lowest layer) in Kelvin!
    dsl4jsb_Real2D_onChunk :: rain
    dsl4jsb_Real2D_onChunk :: snow

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(CARBON_)) RETURN

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Set model and process config and memory
    model => Get_model(tile%owner_model_id)
    !grid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(ASSIMI_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_config(CARBON_)
    dsl4jsb_Get_memory(A2L_)
    !dsl4jsb_Get_config(A2L_)

    ! Set process variables
    dsl4jsb_Get_var2D_onChunk(CARBON_,  pseudo_temp)                ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  pseudo_temp_yDay)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  N_pseudo_temp)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  F_pseudo_temp)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  pseudo_precip)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  pseudo_precip_yDay)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  N_pseudo_precip)            ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  F_pseudo_precip)            ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  topSoilTemp)                ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  alpha)                      ! in
    !
    ! This is on canopy area:
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cconservation_calcCpools )  ! out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_greenwood_2_litter )               ! out
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
    ! xx
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  gross_assimilation_ca )     ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  NPP_pot_rate_ca )           ! inout 
    !dsl4jsb_Get_var2D_onChunk(ASSIMI_,  NPP_pot_Rate_acc )         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  LAI_sum )                   ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_sum )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_sum )                   ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  topSoilTemp_sum )           ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  alpha_sum )                 ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  LAI_yDayMean )              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  LAI_yyDayMean )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_yDayMean )          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_act_yDayMean )          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_yDayMean )              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  topSoilTemp_yDayMean )      ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  alpha_yDayMean )            ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  soil_respiration )          ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_flux_correction)        ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  excess_NPP)                 ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  root_exudates )             ! inout


    IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
      dsl4jsb_Get_var2D_onChunk(CARBON_,  annual_avg_NPP )             ! out
      dsl4jsb_Get_var2D_onChunk(CARBON_,  biennial_avg_NPP )           ! out
      dsl4jsb_Get_var2D_onChunk(CARBON_,  decennial_avg_NPP )          ! out
    ENDIF

    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_green_2_herb )           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_herb_2_littergreen  )       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_herb_2_atm )     ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_npp_2_atm_yday_ta )    ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_npp_2_atm_ta )         ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_soilresp_2_atm_ta )    ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_herb_2_atm_ta )   ! out

    ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
    !dsl4jsb_Get_var2D_onChunk(CARBON_,  fract_litter_wood_new )
    !
    dsl4jsb_Get_var2D_onChunk(PHENO_,   veg_fract_correction )      ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_,   LAI )                       ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_,   fract_fpc_max )             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,     t_air)                      ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,     rain)                       ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,     snow)                       ! in

    ! ---------------------------
    ! Go

    ! Put a comment on the screen...
    IF (debug_on() .AND. iblk==1 .AND. is_newday()  ) CALL message(TRIM(routine), &
                                                           'update_C_NPP_pot_allocation: First time step of new day')
    IF (debug_on() .AND. iblk==1 .AND. is_newyear() ) CALL message(TRIM(routine), &
                                                           'update_C_NPP_pot_allocation: First time step of new year')

    ! Initializations
    cflux_npp_2_atm_yday_ta  = 0._wp
    cflux_npp_2_atm_ta       = 0._wp
    cflux_soilresp_2_atm_ta  = 0._wp
    cflux_herb_2_atm_ta = 0._wp

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
    ! R: NPP_pot_Rate_acc unterscheidet sich in JS3 von NPP_pot_sum darin, daß NPP_pot_rate_acc zu jedem Zeitschritt eine
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


    ! Update pseudo-15day-mean air temperature and precipitation for pseudo_temp_yDay and pseudo_precip_yDay
    ! These variables are needed for the "CALL calc_Cpools", where they given to yasso to calculte C decomposition.
    ! pseudo_temp is a temperature that depicts the smoth course of the soil temperature following the air temperature.
    pseudo_temp   = t_air * N_pseudo_temp  +  F_pseudo_temp   * pseudo_temp
    ! R: I would like to put this precip_total with another name into hydro memory and calculate this
    !    in HYDRO and here just load it from there...
    precip_total  = rain + snow  ! Precipitation rate [kg/(m^2 s)]
    pseudo_precip = precip_total * N_pseudo_precip  +  F_pseudo_precip * pseudo_precip

    lstart = is_time_experiment_start()

    IF( .NOT. is_newday() .OR. lstart) THEN ! perform daily sums if WE ARE NOT STARTING A NEW DAY 
                                            ! or if we JUST STARTED AN TOTALLY NEW EXPERIMENT

       LAI_sum          = LAI_sum          + LAI
       NPP_pot_sum      = NPP_pot_sum      + NPP_pot_rate_ca
       GPP_sum          = GPP_sum          + gross_assimilation_ca
       topSoilTemp_sum  = topSoilTemp_sum  + topSoilTemp ! 2d
       alpha_sum        = alpha_sum        + alpha

    ELSE                                    ! A NEW DAY BEGINS and we DID NOT START AN 
                                            ! TOTALLY NEW EXPERIMENT ==> perform carbon balance

       ! First save as "previous previous days" and "previous days" means
       LAI_yyDayMean  = LAI_yDayMean
       LAI_yDayMean   = LAI_sum/time_steps_per_day()

       NPP_pot_yDayMean     = NPP_pot_sum/time_steps_per_day()
       GPP_yDayMean         = GPP_sum/time_steps_per_day()
       topSoilTemp_yDayMean = topSoilTemp_sum/time_steps_per_day() ! 2d
       alpha_yDayMean       = alpha_sum/time_steps_per_day()

       ! Then restart summing of this days values
       LAI_sum             = LAI
       NPP_pot_sum         = NPP_pot_rate_ca
       GPP_sum             = gross_assimilation_ca
       topSoilTemp_sum     = topSoilTemp ! 2d
       alpha_sum           = alpha
       !NPP_pot_rate_ca_acc        = NPP_pot_rate_ca/time_steps_per_day() ! s.o.

       ! Save other variables (for yasso) as "previous days"
       pseudo_temp_yDay   = pseudo_temp
       pseudo_precip_yDay = pseudo_precip

       ! In case of l_forestRegrowth = true the maximum LAI can change over time
       IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
          MaxLai(:) = dsl4jsb_var2D_onChunk(PHENO_, maxLAI_allom )
       ELSE
        ! R: Rename MaxLAI in lctlib to LAI_max in JS3 comes from theLand%Vegetation%LAI_max,
        ! which is for each vegetated tile taken from MaxLAI in the lct library.
          MaxLai(:) = dsl4jsb_Lctlib_param(MaxLAI)
       END IF

       CALL calculate_current_c_state_sum(tile, options, old_c_state_sum(:))

       ! Now update the C pools
       CALL calc_Cpools(      LAI_yDayMean,                                     & ! in
                              LAI_yyDayMean,                                    & ! in
                              MaxLai,                                           & ! in
                              NPP_pot_yDayMean,                                 & ! in
                              topSoilTemp_yDayMean,                             & ! in
                              alpha_yDayMean,                                   & ! in
                              !
                              dsl4jsb_Lctlib_param(fract_npp_2_woodPool),        & ! in
                              dsl4jsb_Lctlib_param(fract_NPP_2_reservePool),     & ! in
                              dsl4jsb_Lctlib_param(fract_NPP_2_exudates),        & ! in
                              !
                              dsl4jsb_Lctlib_param(fract_green_2_herbivory),     & ! in
                              dsl4jsb_Lctlib_param(tau_c_litter_leaf),      & ! in
                              dsl4jsb_Lctlib_param(tau_c_litter_wood),      & ! in
                              dsl4jsb_Lctlib_param(tau_c_woods),            & ! in
                              dsl4jsb_Lctlib_param(LAI_shed_constant),          & ! in
                              dsl4jsb_Lctlib_param(fract_C_litter_green2atmos),  & ! in
                              dsl4jsb_Lctlib_param(Max_C_content_woods),        & ! in
                              dsl4jsb_Lctlib_param(specificLeafArea_C),         & ! in
                              dsl4jsb_Lctlib_param(reserveC2leafC),             & ! in
                              dsl4jsb_Lctlib_param(PhenologyType),              & ! in  R: is it crop?
                                                                                  ! If yes then this is 5!
                              time_steps_per_day(),                             & ! in
                              dtime,                                            & ! in
                              !
                              c_green,                                      & ! inout
                              c_woods,                                      & ! inout
                              c_reserve,                                    & ! inout
                              c_crop_harvest,                               & ! inout                                     & ! inout
                              !
                              soil_respiration,                                 & ! out
                              NPP_flux_correction,                              & ! out
                              excess_NPP,                                       & ! out
                              root_exudates,                                    & ! out
                              cflux_c_greenwood_2_litter,                                      & ! out
                              cflux_c_green_2_herb,                                  & ! out
                              cflux_herb_2_littergreen,                               & ! out
                              cflux_herb_2_atm,                            & ! out
                              NPP_act_yDayMean,                                 & ! out
                              !
    !fract_litter_wood_new    = fract_litter_wood_new,       & ! inout, optional. Only for spitfire
                              ! variables only needed with yasso:
                              ! R: Note, in the moment the procedure is only used with Yasso.
   temp2_30d               =  pseudo_temp_yDay,                                 & ! in, optional, necessary with yasso
   precip_30d              =  pseudo_precip_yDay,                               & ! in, optional, necessary with yasso
                              !
   c_acid_ag1         =  c_acid_ag1,                                  & ! inout, optional, necessary with yasso
   c_water_ag1        =  c_water_ag1,                                 & ! inout, optional, necessary with yasso
   c_ethanol_ag1      =  c_ethanol_ag1,                               & ! inout, optional, necessary with yasso
   c_nonsoluble_ag1   =  c_nonsoluble_ag1,                            & ! inout, optional, necessary with yasso
   c_acid_bg1         =  c_acid_bg1,                                  & ! inout, optional, necessary with yasso
   c_water_bg1        =  c_water_bg1,                                 & ! inout, optional, necessary with yasso
   c_ethanol_bg1      =  c_ethanol_bg1,                               & ! inout, optional, necessary with yasso
   c_nonsoluble_bg1   =  c_nonsoluble_bg1,                            & ! inout, optional, necessary with yasso
   c_humus_1          =  c_humus_1,                                   & ! inout, optional, necessary with yasso
   c_acid_ag2         =  c_acid_ag2,                                  & ! inout, optional, necessary with yasso
   c_water_ag2        =  c_water_ag2,                                 & ! inout, optional, necessary with yasso
   c_ethanol_ag2      =  c_ethanol_ag2,                               & ! inout, optional, necessary with yasso
   c_nonsoluble_ag2   =  c_nonsoluble_ag2,                            & ! inout, optional, necessary with yasso
   c_acid_bg2         =  c_acid_bg2,                                  & ! inout, optional, necessary with yasso
   c_water_bg2        =  c_water_bg2,                                 & ! inout, optional, necessary with yasso
   c_ethanol_bg2      =  c_ethanol_bg2,                               & ! inout, optional, necessary with yasso
   c_nonsoluble_bg2   =  c_nonsoluble_bg2,                            & ! inout, optional, necessary with yasso
   c_humus_2          =  c_humus_2,                                   & ! inout, optional, necessary with yasso
                              !
   LeafLit_coef1           =  dsl4jsb_Lctlib_param(LeafLit_coef(1)),            & ! in, optional, necessary with yasso
   LeafLit_coef2           =  dsl4jsb_Lctlib_param(LeafLit_coef(2)),            & ! in, optional, necessary with yasso
   LeafLit_coef3           =  dsl4jsb_Lctlib_param(LeafLit_coef(3)),            & ! in, optional, necessary with yasso
   LeafLit_coef4           =  dsl4jsb_Lctlib_param(LeafLit_coef(4)),            & ! in, optional, necessary with yasso
   LeafLit_coef5           =  dsl4jsb_Lctlib_param(LeafLit_coef(5)),            & ! in, optional, necessary with yasso
   WoodLit_coef1           =  dsl4jsb_Lctlib_param(WoodLit_coef(1)),            & ! in, optional, necessary with yasso
   WoodLit_coef2           =  dsl4jsb_Lctlib_param(WoodLit_coef(2)),            & ! in, optional, necessary with yasso
   WoodLit_coef3           =  dsl4jsb_Lctlib_param(WoodLit_coef(3)),            & ! in, optional, necessary with yasso
   WoodLit_coef4           =  dsl4jsb_Lctlib_param(WoodLit_coef(4)),            & ! in, optional, necessary with yasso
   WoodLit_coef5           =  dsl4jsb_Lctlib_param(WoodLit_coef(5)),            & ! in, optional, necessary with yasso
   WoodLitterSize          =  dsl4jsb_Lctlib_param(WoodLitterSize)              & ! in, optional, necessary with yasso
                              )
       ! R: In JSBACH3 wurden alle C-pools auf canopy Fläche gezogen gerechnet. Ebenso die C Fluxe.
       !    Für den Output wurden alle Variablen auch in grid box umgerechnet und entsprechend doppelt angelegt.
       !    Besser wäre für JSBACH4 eigentlich auf Tile Fläche bezogen zu rechnen um alle variablen von JS4
       !    konsistent auf dieselbe Fläche zu beziehen. Sonst machen wir den Code unnötig unübersichtlich.
       !    Dadurch würde z.B. die Variablen-Verdopplung vermieden.
       !    Ich habe die Assimilation zwar schon auf tile Fläche umgerechnet und könnte daher prinzipiell
       !    die Cpools und Fluxe alle über das NPP-auf-Tile-Fläche rechnen.
       !    ABER: in update_Cpools sind eine ganze Reihe von Berechnungen und Konstanten, die dazu ebenfalls angepaßt 
       !    werden müßten. Z.B: "Decomposition of slow soil pool" und evtl. "decomposition rate of green litter" und
       !    specific_leaf_area_C  (=> c_green_max und c_reserve_optimal), Max_C_content_woods. usw.
       !    Dazu kommt, daß Yasso aufgerufen wird und auch hier alle Parametrisierungen durchgeschaut werden müssten.
       !    Wie auch immer wenn später etwas wie "cbalone" verwendet werden soll, dann ist es besser auf canopy
       !    Fläche zu bleiben.
       !    => daher habe ich das doppelte Flächenkonzept von JSBACH3 übernommen.
       !
       !    Um die Benennungen wenigstens richtig zu machen müsste ich jede dieser 30 Variablen in allen
       !    Files (Konsistenz!) mit _ca versehen. Da der Code hier aber ein Provisorium ist:
       !    => C Pools und Fluxe sind ohne _ca benannt und bekommen stattdessen _ta wenn sie auf Tile Fläche bezogen sind.

       ! Compute carbon contents from canopy to tile box area by weighting pools with fractions of grid box covered by vegetation
       ! ------------------------------------------------------------------------------------------------------------------------

       ! C pools
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
        & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:), root_exudates = root_exudates(:), &
        & soil_respiration = soil_respiration(:), NPP_flux_correction = NPP_flux_correction(:), &
 !      &  cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter(:), &
        & cflux_c_green_2_herb = cflux_c_green_2_herb(:), cflux_herb_2_atm = cflux_herb_2_atm(:), &
        & NPP_act_yDayMean = NPP_act_yDayMean(:), NPP_pot_yDayMean = NPP_pot_yDayMean(:), GPP_yDayMean = GPP_yDayMean(:))


      ! cflux_herb_2_atm and soil_respiration are negative (fluxes away from land)
      ! therefore + here instead of -
      current_fluxes = (NPP_act_yDayMean + cflux_herb_2_atm + soil_respiration) &
        & * sec_per_day * veg_fract_correction * fract_fpc_max
      CALL check_carbon_conservation_for_operation(tile, options, old_c_state_sum(:), &
        & current_fluxes(:), cconservation_calcCpools(:))

      IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
         ! calculate weighted (multi) annual increments
         annual_avg_NPP = (annual_avg_NPP * (days_per_year - 1) / days_per_year) &
           & + ( sec_per_day * NPP_act_yDayMean )
         biennial_avg_NPP = (biennial_avg_NPP * ((2 * days_per_year) - 1) / (2 * days_per_year)) &
           & + ( sec_per_day * NPP_act_yDayMean / 2 )
         decennial_avg_NPP = (decennial_avg_NPP * ((10 * days_per_year) - 1) / (10 * days_per_year)) &
           & + ( sec_per_day * NPP_act_yDayMean / 10 )
      ENDIF
   END IF

   ! Compute net CO2 fluxes exchanged with atmosphere at each time step
   !-------------------------------------------------------------------
   ! Note: carbon loss of biosphere means a positive CO2 flux to atmosphere (i.e. NEP and net CO2-flux have opposite signs)

   cflux_npp_2_atm_ta = molarMassCO2_kg *                     & ! Conversion factor from mol to kg CO2
              (- veg_fract_correction * fract_fpc_max          & ! Minus: atmosphere gain is positive
                   * (NPP_pot_rate_ca                         & ! current (not actual) NPP rate
                       - (NPP_pot_yDayMean-NPP_act_yDayMean)))  ! corrected with yesterdays actual NPP defizit

   cflux_soilresp_2_atm_ta = molarMassCO2_kg *                & ! Conversion factor from mol to kg CO2
              (- veg_fract_correction *fract_fpc_max           & ! Minus: atmosphere gain is positive
                   * soil_respiration      )                    ! .. soil respiration

   cflux_herb_2_atm_ta = molarMassCO2_kg *               & ! Conversion factor from mol to kg CO2
              (- veg_fract_correction * fract_fpc_max          & ! Minus: atmosphere gain is positive
                   * cflux_herb_2_atm )                    ! .. herbivory

   cflux_npp_2_atm_yday_ta = molarMassCO2_kg *                & ! C conservarion test in jsb standalone driver 
              (- veg_fract_correction * fract_fpc_max          & ! cannot deal with diurnal cycle of NPP, 
                   * (NPP_act_yDayMean))                        ! therefore here additionaly the day mean


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

    dsl4jsb_Def_memory(CARBON_)

    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract

    TYPE(t_jsb_model), POINTER                :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_C_NPP_pot_allocation'

    INTEGER  :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(CARBON_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    !@todo Those commented out are currently not used by other processes, 
    !     some of them are even not used at all 
    !     + they are not part of the recalc_per_tile_vars
    !       i.e. anyways not updated! 
    ! ++ I do not need them....
!    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_c_greenwood_2_litter_ta,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_green_ta,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_woods_ta,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_reserve_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_crop_harvest_ta,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_acid_ag1_ta,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_water_ag1_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_ethanol_ag1_ta,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_nonsoluble_ag1_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_acid_bg1_ta,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_water_bg1_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_ethanol_bg1_ta,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_nonsoluble_bg1_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_humus_1_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_acid_ag2_ta,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_water_ag2_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_ethanol_ag2_ta,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_nonsoluble_ag2_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_acid_bg2_ta,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_water_bg2_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_ethanol_bg2_ta,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_nonsoluble_bg2_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_humus_2_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, C_sum_natural_ta,         weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_pot_sum_ta,           weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, GPP_sum_ta,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, soil_respiration_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_flux_correction_ta,   weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, excess_NPP_ta,            weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yyDayMean_ta,         weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yDayMean_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_pot_yDayMean_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_act_yDayMean_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, GPP_yDayMean_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, root_exudates_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_c_green_2_herb_ta,       weighted_by_fract)
   ! dsl4jsb_Aggregate_onChunk(CARBON_, cflux_herb_2_littergreen_ta, weighted_by_fract)
   ! dsl4jsb_Aggregate_onChunk(CARBON_, cflux_herb_2_atm_ta, weighted_by_fract)
   ! dsl4jsb_Aggregate_onChunk(CARBON_, fract_litter_wood_new_ta,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_npp_2_atm_yday_ta,  weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(CARBON_, cconservation_calcCpools, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_npp_2_atm_ta,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_soilresp_2_atm_ta,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_C_NPP_pot_allocation

  ! ================================================================================================================================
  !>
  !! @brief Task "update_apply_age_cohort_changes" applies area and c changes to the age cohorts (using delta variables) 
  !! - TODO-JN-COHORTS area change should probably be a task of another process
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  ! 
  SUBROUTINE update_apply_age_cohort_changes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    LOGICAL :: l_isFirstChild
    REAL(wp), DIMENSION(options%nc)        :: source_fract

    INTEGER  :: iblk, ics, ice, nc
    TYPE(t_jsb_model), POINTER :: model
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_apply_age_cohort_changes'

    ! ---------------------------
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    model => Get_model(tile%owner_model_id)

    IF (.NOT. ASSOCIATED(tile%parent_tile)) RETURN
    IF (.NOT. tile%parent_tile%Is_process_active(FAGE_)) RETURN
    IF (.NOT. dsl4jsb_Lctlib_param(ForestFlag)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    IF(ANY(tile%parent_tile%l_delta_area_change(ics:ice, iblk))) THEN     
      IF((tile%Is_first_child()) .OR. ANY(tile%l_update_c(ics:ice, iblk))) THEN
        ! a) TODO-JN-COHORTS what else to update? 

        ! b) update C states with deltas (+ reset C deltas to 0)
        CALL apply_delta_carbon_to_carbon_states(tile, options)
      ENDIF

      ! c) update area fractions 
      ! TODO-JN-COHORTS: maybe make it a function in mo_jsb_tile_class: update_fraction?
      ! Or make it an own task in a more suited process?
      CALL tile%Get_fraction(ics, ice, iblk, fract=source_fract(:))
      CALL tile%Set_fraction(ics, ice, iblk, fract=source_fract(:) + tile%delta_fract(ics:ice, iblk))
      ! + reset deltas to 0
      tile%delta_fract(ics:ice, iblk) = 0._wp

    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_apply_age_cohort_changes

  SUBROUTINE aggregate_apply_age_cohort_changes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_memory(CARBON_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_apply_age_cohort_changes'

    INTEGER  :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(CARBON_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    IF(ANY(tile%l_update_c(ics:ice, iblk))) THEN 
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_green_ta,               weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_woods_ta,               weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_reserve_ta,             weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_crop_harvest_ta,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_ag1_ta,           weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_ag1_ta,          weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_ag1_ta,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_ag1_ta,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_ag2_ta,           weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_ag2_ta,          weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_ag2_ta,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_ag2_ta,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_bg1_ta,           weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_bg1_ta,          weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_bg1_ta,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_bg1_ta,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_bg2_ta,           weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_bg2_ta,          weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_bg2_ta,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_bg2_ta,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_humus_1_ta,            weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  c_humus_2_ta,            weighted_by_fract)

      ! in case of disturbance also these variables need to be aggregated:
      dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_dist_greenreserve_2_soil_ta, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_dist_woods_2_soil_ta,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_fire_all_2_atm_ta,           weighted_by_fract)

      !TODO-JN-COHORTS: other aggregations required?
    ENDIF

    !NOTE: despite changes, currently no aggregation of fractions because:
    !      for fage fractions are only redistributed on child tiles, i.e. no changes on above tiles

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_apply_age_cohort_changes



!##########################################################################################################################
! other subroutines (no further tasks)
!##########################################################################################################################

  ! ================================================================================================================================
  !>
  !! @brief scales all carbon state variables upon reference area change
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
    REAL(wp),                   INTENT(in)            :: oldRefArea(:),  &
                                                         newRefArea(:)

    dsl4jsb_Def_memory(CARBON_)

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':rescale_carbon_upon_reference_area_change'
    INTEGER                     :: ics, ice, iblk

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

    ! Assertion: 0 < oldRedArea <= 1 and 0 < newRefArea <= 1 
    IF ( ANY(oldRefArea .LE. 0) .OR. ANY(oldRefArea .GT. 1) ) THEN 
      CALL finish(TRIM(routine), &
        & 'Violation of assertion: Reference areas (here old) need to be >0 and <=1 - min: '&
        & //real2string(MINVAL(oldRefArea)) //'; max: '//real2string(MAXVAL(oldRefArea))  )
    ENDIF
    IF ( ANY(newRefArea .LE. 0) .OR. ANY(newRefArea .GT. 1) ) THEN 
      CALL finish(TRIM(routine), &
        & 'Violation of assertion: Reference areas (here new) need to be >0 and <=1 - min: '&
        & //real2string(MINVAL(newRefArea)) //'; max: '//real2string(MAXVAL(newRefArea))  )
    ENDIF

    dsl4jsb_Get_memory(CARBON_)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest)         ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1)            ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1)         ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1)      ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1)            ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1)           ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1)         ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1)      ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1)             ! inout 

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2)           ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2)          ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2)        ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2)          ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2)        ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2)     ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2)            ! inout 

    c_green           = c_green           * (oldRefArea / newRefArea)
    c_woods           = c_woods           * (oldRefArea / newRefArea)
    c_reserve         = c_reserve         * (oldRefArea / newRefArea)
    c_crop_harvest    = c_crop_harvest    * (oldRefArea / newRefArea)

    c_acid_ag1       = c_acid_ag1       * (oldRefArea / newRefArea)
    c_water_ag1      = c_water_ag1      * (oldRefArea / newRefArea)
    c_ethanol_ag1    = c_ethanol_ag1    * (oldRefArea / newRefArea)
    c_nonsoluble_ag1 = c_nonsoluble_ag1 * (oldRefArea / newRefArea)
    c_acid_bg1       = c_acid_bg1       * (oldRefArea / newRefArea)
    c_water_bg1      = c_water_bg1      * (oldRefArea / newRefArea)
    c_ethanol_bg1    = c_ethanol_bg1    * (oldRefArea / newRefArea)
    c_nonsoluble_bg1 = c_nonsoluble_bg1 * (oldRefArea / newRefArea)
    c_humus_1        = c_humus_1        * (oldRefArea / newRefArea)
    c_acid_ag2       = c_acid_ag2       * (oldRefArea / newRefArea)
    c_water_ag2      = c_water_ag2      * (oldRefArea / newRefArea)
    c_ethanol_ag2    = c_ethanol_ag2    * (oldRefArea / newRefArea)
    c_nonsoluble_ag2 = c_nonsoluble_ag2 * (oldRefArea / newRefArea)
    c_acid_bg2       = c_acid_bg2       * (oldRefArea / newRefArea)
    c_water_bg2      = c_water_bg2      * (oldRefArea / newRefArea)
    c_ethanol_bg2    = c_ethanol_bg2    * (oldRefArea / newRefArea)
    c_nonsoluble_bg2 = c_nonsoluble_bg2 * (oldRefArea / newRefArea)
    c_humus_2        = c_humus_2        * (oldRefArea / newRefArea)

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
      & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:))

  END SUBROUTINE rescale_carbon_upon_reference_area_change


  ! ================================================================================================================================
  !>
  !! @brief updates delta of all carbon state variables upon area shift
  !!
  !! @todo discuss how to make more general, e.g. get the info on which state variables to work on from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  ! NOTE: although the density would not change when moving C away the calculation is done symmetrically  
  !       because it makes the transfer from delta to states easier afterwards, particularly since 
  !       area losses and gains happen simultaneously

  SUBROUTINE recalc_delta_carbon_rp_of_source_and_sink_upon_area_shift( &
      sourceT, sinkT, options, delta_fract, prev_fract )

    USE mo_carbon_process,      ONLY: update_delta_carbon_rp

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: sourceT, sinkT
    TYPE(t_jsb_task_options),   INTENT(in)            :: options
    REAL(wp),                   INTENT(in)            :: delta_fract(:), prev_fract(:)

    TYPE(t_jsb_model), POINTER                :: model
    dsl4jsb_Def_config(PHENO_)

    dsl4jsb_Def_memory_tile(CARBON_, sourceT)
    dsl4jsb_Def_memory_tile(CARBON_, sinkT)
    dsl4jsb_Def_memory_tile(PHENO_, sourceT)
    dsl4jsb_Def_memory_tile(PHENO_, sinkT)

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = &
      & TRIM(modname)//':recalc_delta_carbon_rp_of_source_and_sink_upon_area_shift'
    INTEGER                     :: ics, ice, iblk, nc, i

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

    dsl4jsb_Real2D_onChunk ::  GPP_sum
    dsl4jsb_Real2D_onChunk ::  NPP_pot_sum

    dsl4jsb_Real2D_onChunk ::  veg_fract_correction

    dsl4jsb_Real2D_onChunk ::  current_max_green
    dsl4jsb_Real2D_onChunk ::  veg_carbon_at_max_green

    dsl4jsb_Real2D_onChunk ::  day_npp_sum

    dsl4jsb_Real2D_onChunk ::  annual_avg_NPP
    dsl4jsb_Real2D_onChunk ::  biennial_avg_NPP
    dsl4jsb_Real2D_onChunk ::  decennial_avg_NPP


    dsl4jsb_Real2D_onChunk ::  delta_c_green_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_woods_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_reserve_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_crop_harvest_rp_sourceT

    dsl4jsb_Real2D_onChunk ::  delta_c_acid_ag1_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_water_ag1_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_ag1_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_ag1_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_acid_bg1_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_water_bg1_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_bg1_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_bg1_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_humus_1_rp_sourceT

    dsl4jsb_Real2D_onChunk ::  delta_c_acid_ag2_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_water_ag2_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_ag2_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_ag2_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_acid_bg2_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_water_bg2_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_bg2_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_bg2_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_c_humus_2_rp_sourceT

    dsl4jsb_Real2D_onChunk ::  delta_current_max_green_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_veg_carbon_at_max_green_rp_sourceT

    dsl4jsb_Real2D_onChunk ::  delta_day_npp_sum_rp_sourceT

    dsl4jsb_Real2D_onChunk ::  delta_annual_avg_NPP_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_biennial_avg_NPP_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_decennial_avg_NPP_rp_sourceT

    dsl4jsb_Real2D_onChunk ::  delta_GPP_sum_rp_sourceT
    dsl4jsb_Real2D_onChunk ::  delta_NPP_pot_sum_rp_sourceT

    dsl4jsb_Real2D_onChunk ::  delta_c_green_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_woods_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_reserve_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_crop_harvest_rp_sinkT

    dsl4jsb_Real2D_onChunk ::  delta_c_acid_ag1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_water_ag1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_ag1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_ag1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_acid_bg1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_water_bg1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_bg1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_bg1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_humus_1_rp_sinkT

    dsl4jsb_Real2D_onChunk ::  delta_c_acid_ag2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_water_ag2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_ag2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_ag2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_acid_bg2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_water_bg2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_bg2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_bg2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_humus_2_rp_sinkT

    dsl4jsb_Real2D_onChunk ::  delta_current_max_green_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_veg_carbon_at_max_green_rp_sinkT

    dsl4jsb_Real2D_onChunk ::  delta_day_npp_sum_rp_sinkT

    dsl4jsb_Real2D_onChunk ::  delta_annual_avg_NPP_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_biennial_avg_NPP_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_decennial_avg_NPP_rp_sinkT

    dsl4jsb_Real2D_onChunk ::  delta_GPP_sum_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_NPP_pot_sum_rp_sinkT

    !-------------------------
    model => Get_model(sourceT%owner_model_id)

    dsl4jsb_Get_config(PHENO_)

    dsl4jsb_Get_memory_tile(CARBON_, sourceT)
    dsl4jsb_Get_memory_tile(CARBON_, sinkT)

    dsl4jsb_Get_memory_tile(PHENO_, sinkT)
    dsl4jsb_Get_memory_tile(PHENO_, sourceT)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    !Assertion prev_fract + delta_fract .GE. 0._wp
    IF (ANY(prev_fract(:) + delta_fract(:) .LT. 0_wp )) THEN
      DO i = 1,nc
        IF (prev_fract(i) + delta_fract(i) .LT. 0_wp ) THEN
          WRITE (message_text,*) 'Violation of assertion: area to be changed (', delta_fract(i), &
            & ') needs to be .LE. existing area (', prev_fract(i) ,')'
          CALL finish(TRIM(routine), message_text)
        ENDIF
      ENDDO
    ENDIF
    ! TODO-JN-COHORTS: if smallest possible fraction 1e-10 will be re-introduced:  
    !    then the following should not happen and other Assertions and MAX statements will be required!

    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_green, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_woods, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_reserve, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_crop_harvest, sourceT)

    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_acid_ag1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_water_ag1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_ethanol_ag1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_nonsoluble_ag1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_acid_bg1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_water_bg1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_ethanol_bg1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_nonsoluble_bg1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_humus_1, sourceT)

    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_acid_ag2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_water_ag2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_ethanol_ag2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_nonsoluble_ag2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_acid_bg2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_water_bg2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_ethanol_bg2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_nonsoluble_bg2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_humus_2, sourceT)

    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, GPP_sum, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, NPP_pot_sum, sourceT)

    dsl4jsb_Get_var2d_onChunk_tile(PHENO_, veg_fract_correction, sourceT)

    dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, delta_day_npp_sum_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, delta_day_npp_sum_rp, sinkT)

    IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
      dsl4jsb_Get_var2d_onChunk_tile(CARBON_, annual_avg_NPP, sourceT)
      dsl4jsb_Get_var2d_onChunk_tile(CARBON_, biennial_avg_NPP, sourceT)
      dsl4jsb_Get_var2d_onChunk_tile(CARBON_, decennial_avg_NPP, sourceT)
      
      dsl4jsb_Get_var2d_onChunk_tile(CARBON_, current_max_green, sourceT)
      dsl4jsb_Get_var2d_onChunk_tile(CARBON_, veg_carbon_at_max_green, sourceT)

      dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_current_max_green_rp, sourceT)
      dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_veg_carbon_at_max_green_rp, sourceT)

      dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_current_max_green_rp, sinkT)
      dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_veg_carbon_at_max_green_rp, sinkT)

      dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_annual_avg_NPP_rp, sourceT)
      dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_biennial_avg_NPP_rp, sourceT)
      dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_decennial_avg_NPP_rp, sourceT)

      dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_annual_avg_NPP_rp, sinkT)
      dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_biennial_avg_NPP_rp, sinkT)
      dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_decennial_avg_NPP_rp, sinkT)
    ENDIF

    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_green_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_woods_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_reserve_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_crop_harvest_rp, sourceT)

    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_acid_ag1_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_water_ag1_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_ethanol_ag1_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_nonsoluble_ag1_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_acid_bg1_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_water_bg1_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_ethanol_bg1_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_nonsoluble_bg1_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_humus_1_rp, sourceT)

    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_acid_ag2_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_water_ag2_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_ethanol_ag2_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_nonsoluble_ag2_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_acid_bg2_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_water_bg2_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_ethanol_bg2_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_nonsoluble_bg2_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_humus_2_rp, sourceT)

    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_GPP_sum_rp, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_NPP_pot_sum_rp, sourceT)

    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_green_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_woods_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_reserve_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_crop_harvest_rp, sinkT)

    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_acid_ag1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_water_ag1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_ethanol_ag1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_nonsoluble_ag1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_acid_bg1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_water_bg1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_ethanol_bg1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_nonsoluble_bg1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_humus_1_rp, sinkT)

    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_acid_ag2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_water_ag2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_ethanol_ag2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_nonsoluble_ag2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_acid_bg2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_water_bg2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_ethanol_bg2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_nonsoluble_bg2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_humus_2_rp, sinkT)

    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_GPP_sum_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_NPP_pot_sum_rp, sinkT)

    CALL update_delta_carbon_rp(c_green(:), delta_c_green_rp_sourceT(:), &
      & delta_c_green_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_woods(:), delta_c_woods_rp_sourceT(:), &
      & delta_c_woods_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_reserve(:), delta_c_reserve_rp_sourceT(:), &
      & delta_c_reserve_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_crop_harvest(:), delta_c_crop_harvest_rp_sourceT(:), &
      delta_c_crop_harvest_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))

    CALL update_delta_carbon_rp(c_acid_ag1(:), delta_c_acid_ag1_rp_sourceT(:), &
      & delta_c_acid_ag1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_water_ag1(:), delta_c_water_ag1_rp_sourceT(:), &
      & delta_c_water_ag1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_ethanol_ag1(:), delta_c_ethanol_ag1_rp_sourceT(:), &
      & delta_c_ethanol_ag1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_nonsoluble_ag1(:), delta_c_nonsoluble_ag1_rp_sourceT(:), &
      & delta_c_nonsoluble_ag1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_acid_bg1(:), delta_c_acid_bg1_rp_sourceT(:), &
      & delta_c_acid_bg1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_water_bg1(:), delta_c_water_bg1_rp_sourceT(:), &
      & delta_c_water_bg1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_ethanol_bg1(:), delta_c_ethanol_bg1_rp_sourceT(:), &
      & delta_c_ethanol_bg1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_nonsoluble_bg1(:), delta_c_nonsoluble_bg1_rp_sourceT(:), &
      & delta_c_nonsoluble_bg1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_humus_1(:), delta_c_humus_1_rp_sourceT(:), &
      & delta_c_humus_1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))

    CALL update_delta_carbon_rp(c_acid_ag2(:), delta_c_acid_ag2_rp_sourceT(:), &
      & delta_c_acid_ag2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_water_ag2(:), delta_c_water_ag2_rp_sourceT(:), &
      & delta_c_water_ag2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_ethanol_ag2(:), delta_c_ethanol_ag2_rp_sourceT(:), &
      & delta_c_ethanol_ag2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_nonsoluble_ag2(:), delta_c_nonsoluble_ag2_rp_sourceT(:), &
      & delta_c_nonsoluble_ag2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_acid_bg2(:), delta_c_acid_bg2_rp_sourceT(:), &
      & delta_c_acid_bg2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_water_bg2(:), delta_c_water_bg2_rp_sourceT(:), &
      & delta_c_water_bg2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_ethanol_bg2(:), delta_c_ethanol_bg2_rp_sourceT(:), &
      & delta_c_ethanol_bg2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_nonsoluble_bg2(:), delta_c_nonsoluble_bg2_rp_sourceT(:), &
      & delta_c_nonsoluble_bg2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_humus_2(:), delta_c_humus_2_rp_sourceT(:), &
      & delta_c_humus_2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))

    CALL update_delta_carbon_rp(GPP_sum(:), delta_GPP_sum_rp_sourceT(:), &
      & delta_GPP_sum_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(NPP_pot_sum(:), delta_NPP_pot_sum_rp_sourceT(:), &
      & delta_NPP_pot_sum_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))

    !TODO-JN-COHORTS day_npp_sum is not a carbon variable but belongs to pheno 
    ! -- nevertheless averaged for now here in mo_carbon_interface, since for now its the 
    !    only non carbon variable, and its too much effort to make an own task
    CALL update_delta_carbon_rp(day_npp_sum(:), delta_day_npp_sum_rp_sourceT(:), &
      & delta_day_npp_sum_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))

    IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
      CALL update_delta_carbon_rp(current_max_green(:), delta_current_max_green_rp_sourceT(:), &
        & delta_current_max_green_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(veg_carbon_at_max_green(:), delta_veg_carbon_at_max_green_rp_sourceT(:), &
        & delta_veg_carbon_at_max_green_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))

      CALL update_delta_carbon_rp(annual_avg_NPP(:), delta_annual_avg_NPP_rp_sourceT(:), &
        & delta_annual_avg_NPP_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(biennial_avg_NPP(:), delta_biennial_avg_NPP_rp_sourceT(:), &
        & delta_biennial_avg_NPP_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(decennial_avg_NPP(:), delta_decennial_avg_NPP_rp_sourceT(:), &
        & delta_decennial_avg_NPP_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))

    ENDIF

  END SUBROUTINE recalc_delta_carbon_rp_of_source_and_sink_upon_area_shift


  ! ================================================================================================================================
  !>
  !! @brief updates delta of carbon state variables upon disturbance (including clear-cuts!)
  !!        -> source densities stay the same
  !!        -> sink (first age cohort) gets fraction of source litter/soil 
  !!            + fraction of green/reserve/wood are also distributed into litter/soil
  !!            + in case of fire some of the litter pools are redistributed
  !!
  !! @todo discuss how to make more general, e.g. get the info on which state variables to work on from somewhere else?
  !!       or as a vector ?
  !!
  !! @param[in,out] sourceT      Source tile of area and thus carbon - usually no carbon density changes here 
  !! @param[in,out] sinkT        Sink tiles to which area and thus carbon is added (here: always first child!)
  !! @param[in] l_sourceisfirstChild Is first child the source? Then there'll be changes in vegetation carbon densities 
  !! @param[in]     options          Additional run-time parameters
  !! @param[in]     delta_fract      Change in area
  !! @param[in]     prev_fract       Source fract before change
  !! @param[in]     distKind         Kind of the executed disturbance - wind, fire, clearCut 
  !
  SUBROUTINE recalc_delta_carbon_rp_upon_disturbance( &
      sourceT, sinkT, l_sourceisfirstChild, options, delta_fract, prev_fract, distKind )

    USE mo_carbon_process,      ONLY: update_delta_carbon_rp, relocate_carbon_fire, relocate_carbon_damage
    USE mo_carbon_constants,    ONLY: WIND, FIRE, CLEARCUT

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: sourceT, sinkT
    LOGICAL,                    INTENT(in)            :: l_sourceisfirstChild
    TYPE(t_jsb_task_options),   INTENT(in)            :: options
    REAL(wp),                   INTENT(in)            :: delta_fract(:), prev_fract(:)
    INTEGER,                    INTENT(in)            :: distKind

    TYPE(t_jsb_model), POINTER                :: model
    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_config(CARBON_)

    dsl4jsb_Def_memory_tile(CARBON_, sourceT)
    dsl4jsb_Def_memory_tile(CARBON_, sinkT)
    dsl4jsb_Def_memory_tile(PHENO_, sourceT)

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = &
      & TRIM(modname)//':recalc_delta_carbon_rp_upon_disturbance'
    INTEGER                     :: ics, ice, iblk, nc, i
    REAL(wp), DIMENSION(options%nc) :: helper, d_green_new, d_reserve_new, d_wood_new
    REAL(wp), DIMENSION(options%nc) :: d_acid_ag1_new, d_water_ag1_new, d_ethanol_ag1_new, d_nonsoluble_ag1_new
    REAL(wp), DIMENSION(options%nc) :: d_acid_ag2_new, d_water_ag2_new, d_ethanol_ag2_new, d_nonsoluble_ag2_new

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

    dsl4jsb_Real2D_onChunk ::  veg_fract_correction

    dsl4jsb_Real2D_onChunk ::  delta_cflux_dist_greenreserve_2_soil_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_cflux_dist_woods_2_soil_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_cflux_fire_all_2_atm_rp_sinkT

    dsl4jsb_Real2D_onChunk ::  delta_c_green_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_woods_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_reserve_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_crop_harvest_rp_sinkT

    dsl4jsb_Real2D_onChunk ::  delta_c_acid_ag1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_water_ag1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_ag1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_ag1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_acid_bg1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_water_bg1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_bg1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_bg1_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_humus_1_rp_sinkT

    dsl4jsb_Real2D_onChunk ::  delta_c_acid_ag2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_water_ag2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_ag2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_ag2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_acid_bg2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_water_bg2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_bg2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_bg2_rp_sinkT
    dsl4jsb_Real2D_onChunk ::  delta_c_humus_2_rp_sinkT

    ! TODO-JN-COHORTS: if the smallest possible fraction 1e-10 is re-introduced 
    ! then the following should not happen and other Assertions and MAX statements will be required!

    model => Get_model(sourceT%owner_model_id)

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_config(CARBON_)

    dsl4jsb_Get_memory_tile(CARBON_, sourceT)
    dsl4jsb_Get_memory_tile(CARBON_, sinkT)
    dsl4jsb_Get_memory_tile(PHENO_, sourceT)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    ice   = options%nc

    !Assertion prev_fract + delta_fract .GE. 0._wp
    IF (ANY(prev_fract(:) + delta_fract(:) .LT. 0_wp )) THEN
      DO i = 1,nc
        IF (prev_fract(i) + delta_fract(i) .LT. 0_wp ) THEN
          WRITE (message_text,*) 'Violation of assertion: area to be changed (', delta_fract(i), &
            & ') needs to be .LE. existing area (', prev_fract(i) ,')'
          CALL finish(TRIM(routine), message_text)
        ENDIF
      ENDDO
    ENDIF

    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_green, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_woods, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_reserve, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_crop_harvest, sourceT)

    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_acid_ag1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_water_ag1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_ethanol_ag1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_nonsoluble_ag1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_acid_bg1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_water_bg1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_ethanol_bg1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_nonsoluble_bg1, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_humus_1, sourceT)

    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_acid_ag2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_water_ag2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_ethanol_ag2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_nonsoluble_ag2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_acid_bg2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_water_bg2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_ethanol_bg2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_nonsoluble_bg2, sourceT)
    dsl4jsb_Get_var2d_onChunk_tile(CARBON_, c_humus_2, sourceT)

    dsl4jsb_Get_var2d_onChunk_tile(PHENO_, veg_fract_correction, sourceT)

    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_cflux_dist_greenreserve_2_soil_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_cflux_dist_woods_2_soil_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_cflux_fire_all_2_atm_rp, sinkT)

    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_green_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_woods_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_reserve_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_crop_harvest_rp, sinkT)

    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_acid_ag1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_water_ag1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_ethanol_ag1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_nonsoluble_ag1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_acid_bg1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_water_bg1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_ethanol_bg1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_nonsoluble_bg1_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_humus_1_rp, sinkT)

    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_acid_ag2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_water_ag2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_ethanol_ag2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_nonsoluble_ag2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_acid_bg2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_water_bg2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_ethanol_bg2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_nonsoluble_bg2_rp, sinkT)
    dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, delta_c_humus_2_rp, sinkT)

    d_green_new(:) = 0._wp
    d_reserve_new(:) = 0._wp
    d_wood_new(:) = 0._wp

    d_acid_ag1_new(:) = 0._wp 
    d_water_ag1_new(:) = 0._wp
    d_ethanol_ag1_new(:) = 0._wp
    d_nonsoluble_ag1_new(:) = 0._wp
    d_acid_ag2_new(:) = 0._wp
    d_water_ag2_new(:) = 0._wp
    d_ethanol_ag2_new(:) = 0._wp
    d_nonsoluble_ag2_new(:) = 0._wp

    ! For all implemented disturbance kinds the vegetation carbon on the disturbed area 
    ! is entirely distributed to other pools, therefore, additional pools are required 
    CALL update_delta_carbon_rp(c_green(:), helper(:), d_green_new(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_woods(:), helper(:), d_wood_new(:), delta_fract(:), veg_fract_correction(:))
    CALL update_delta_carbon_rp(c_reserve(:), helper(:), d_reserve_new(:), delta_fract(:), veg_fract_correction(:))

    IF ((distKind .EQ. FIRE)) THEN
      ! For fire disturbance, parts of the above ground litter carbon are burned 
      ! to deal with the first and all other childs in the same function call, therefore, additional pools are used
      CALL update_delta_carbon_rp(c_acid_ag1(:), helper(:), d_acid_ag1_new(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_water_ag1(:), helper(:), d_water_ag1_new(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_ethanol_ag1(:), helper(:), &
        & d_ethanol_ag1_new(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_nonsoluble_ag1(:), helper(:), &
        & d_nonsoluble_ag1_new(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_acid_ag2(:), helper(:), d_acid_ag2_new(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_water_ag2(:), helper(:), d_water_ag2_new(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_ethanol_ag2(:), helper(:), &
        & d_ethanol_ag2_new(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_nonsoluble_ag2(:), helper(:), &
        & d_nonsoluble_ag2_new(:), delta_fract(:), veg_fract_correction(:))
    ELSEIF (.NOT. l_sourceisfirstChild) THEN
      ! For other disturbances above ground litter is only redistributed 
      ! However, for the first child the damaged area stays within the first age cohort and so do the ag litter pools
      CALL update_delta_carbon_rp(c_acid_ag1(:), helper(:), &
        & delta_c_acid_ag1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_water_ag1(:), helper(:), &
        & delta_c_water_ag1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_ethanol_ag1(:), helper(:), &
        & delta_c_ethanol_ag1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_nonsoluble_ag1(:), helper(:), &
        & delta_c_nonsoluble_ag1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))

      CALL update_delta_carbon_rp(c_acid_ag2(:), helper(:), &
        & delta_c_acid_ag2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_water_ag2(:), helper(:), &
        & delta_c_water_ag2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_ethanol_ag2(:), helper(:), &
        & delta_c_ethanol_ag2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_nonsoluble_ag2(:), helper(:), &
        & delta_c_nonsoluble_ag2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    ENDIF

    IF (.NOT. l_sourceisfirstChild) THEN
      ! As opposed to the ag litter, the bg litter/soil is not changed by fire, thus here it is always the case that 
      ! when the first child is the source, the damaged area stays within the first age cohort and same for bg soil/litter pools
      CALL update_delta_carbon_rp(c_crop_harvest(:), helper(:), &
        delta_c_crop_harvest_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))

      CALL update_delta_carbon_rp(c_acid_bg1(:), helper(:), &
        & delta_c_acid_bg1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_water_bg1(:), helper(:), &
        & delta_c_water_bg1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_ethanol_bg1(:), helper(:), &
        & delta_c_ethanol_bg1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_nonsoluble_bg1(:), helper(:), &
        & delta_c_nonsoluble_bg1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_humus_1(:), helper(:), &
        & delta_c_humus_1_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))

      CALL update_delta_carbon_rp(c_acid_bg2(:), helper(:), &
        & delta_c_acid_bg2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_water_bg2(:), helper(:), &
        & delta_c_water_bg2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_ethanol_bg2(:), helper(:), &
        & delta_c_ethanol_bg2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_nonsoluble_bg2(:), helper(:), &
        & delta_c_nonsoluble_bg2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
      CALL update_delta_carbon_rp(c_humus_2(:), helper(:), &
        & delta_c_humus_2_rp_sinkT(:), delta_fract(:), veg_fract_correction(:))
    ENDIF

    IF (l_sourceisfirstChild) THEN
      ! If source and sink is the same child - i.e. the first child, 
      !  changes in green, reserve and wood need to be subtracted from the delta carbon pools 
      !  after being distributed to the soil/litter (or probably later in case of harvest to a harvest pool)
      delta_c_green_rp_sinkT(:) = delta_c_green_rp_sinkT(:) - d_green_new(:)
      delta_c_reserve_rp_sinkT(:) = delta_c_reserve_rp_sinkT(:) - d_reserve_new(:)
      delta_c_woods_rp_sinkT(:) = delta_c_woods_rp_sinkT(:) - d_wood_new(:)
    ENDIF

    helper(:) = 1._wp

    IF((distKind .EQ. WIND) .OR. (distKind .EQ. CLEARCUT)) THEN
      !TODO-JN-COHORTS: currently this is used for harvest, too: change as soon as there'll be product pools!

      ! Redistribute the carbon among the pools according to wind damage
      CALL relocate_carbon_damage( helper(:), & ! in
        & model%lctlib(sinkT%lcts(1)%lib_id)%LeafLit_coef(1:5), & ! in
        & model%lctlib(sinkT%lcts(1)%lib_id)%WoodLit_coef(1:5), & ! in
        & d_green_new(:),                     & ! inout
        & d_wood_new(:),                      & ! inout
        & d_reserve_new(:),                   & ! inout
        & delta_cflux_dist_greenreserve_2_soil_rp_sinkT(:),  & ! inout
        & delta_cflux_dist_woods_2_soil_rp_sinkT(:),   & ! inout
        & delta_c_acid_ag1_rp_sinkT(:),            & ! inout
        & delta_c_water_ag1_rp_sinkT(:),           & ! inout
        & delta_c_ethanol_ag1_rp_sinkT(:),         & ! inout
        & delta_c_nonsoluble_ag1_rp_sinkT(:),      & ! inout
        & delta_c_acid_ag2_rp_sinkT(:),            & ! inout
        & delta_c_water_ag2_rp_sinkT(:),           & ! inout
        & delta_c_ethanol_ag2_rp_sinkT(:),         & ! inout
        & delta_c_nonsoluble_ag2_rp_sinkT(:),      & ! inout
        & delta_c_acid_bg1_rp_sinkT(:),            & ! inout
        & delta_c_water_bg1_rp_sinkT(:),           & ! inout
        & delta_c_ethanol_bg1_rp_sinkT(:),         & ! inout
        & delta_c_nonsoluble_bg1_rp_sinkT(:),      & ! inout
        & delta_c_acid_bg2_rp_sinkT(:),            & ! inout
        & delta_c_water_bg2_rp_sinkT(:),           & ! inout
        & delta_c_ethanol_bg2_rp_sinkT(:),         & ! inout
        & delta_c_nonsoluble_bg2_rp_sinkT(:),      & ! inout
        & delta_c_humus_1_rp_sinkT(:),             & ! inout
        & delta_c_humus_2_rp_sinkT(:)              & ! inout
        & )

    ELSEIF (distKind .EQ. FIRE) THEN 

      IF (l_sourceisfirstChild) THEN
        ! in case of fire, the aboveground litter is altered, when the firest child is the source, some of the 
        ! above ground carbon therefore will need to be substracted. 
        ! Here before the function call below, all ag litter carbon on the disturbed area is subtracted
        ! -> below the function call the not-burned part will be then added back   
        delta_c_acid_ag1_rp_sinkT = delta_c_acid_ag1_rp_sinkT - d_acid_ag1_new
        delta_c_water_ag1_rp_sinkT = delta_c_water_ag1_rp_sinkT - d_water_ag1_new
        delta_c_ethanol_ag1_rp_sinkT = delta_c_ethanol_ag1_rp_sinkT - d_ethanol_ag1_new
        delta_c_nonsoluble_ag1_rp_sinkT = delta_c_nonsoluble_ag1_rp_sinkT - d_nonsoluble_ag1_new
        delta_c_acid_ag2_rp_sinkT = delta_c_acid_ag2_rp_sinkT - d_acid_ag2_new
        delta_c_water_ag2_rp_sinkT = delta_c_water_ag2_rp_sinkT - d_water_ag2_new
        delta_c_ethanol_ag2_rp_sinkT = delta_c_ethanol_ag2_rp_sinkT - d_ethanol_ag2_new
        delta_c_nonsoluble_ag2_rp_sinkT = delta_c_nonsoluble_ag2_rp_sinkT - d_nonsoluble_ag2_new
      ENDIF

      CALL relocate_carbon_fire( &
        & helper(:),                     & ! in
        & dsl4jsb_Config(CARBON_)%fire_fract_wood_2_atmos , & ! in
        & model%lctlib(sinkT%lcts(1)%lib_id)%LeafLit_coef(1:5), & ! in
        & model%lctlib(sinkT%lcts(1)%lib_id)%WoodLit_coef(1:5), & ! in
        & d_green_new(:),                & ! inout
        & d_wood_new(:),                 & ! inout
        & d_reserve_new(:),              & ! inout
        & delta_cflux_dist_greenreserve_2_soil_rp_sinkT(:),  & ! inout
        & delta_cflux_dist_woods_2_soil_rp_sinkT(:),   & ! inout
        & d_acid_ag1_new(:),            & ! inout
        & d_water_ag1_new(:),           & ! inout
        & d_ethanol_ag1_new(:),         & ! inout
        & d_nonsoluble_ag1_new(:),      & ! inout
        & d_acid_ag2_new(:),            & ! inout
        & d_water_ag2_new(:),           & ! inout
        & d_ethanol_ag2_new(:),         & ! inout
        & d_nonsoluble_ag2_new(:),      & ! inout
        & delta_c_acid_bg1_rp_sinkT(:),            & ! inout
        & delta_c_water_bg1_rp_sinkT(:),           & ! inout
        & delta_c_ethanol_bg1_rp_sinkT(:),         & ! inout
        & delta_c_nonsoluble_bg1_rp_sinkT(:),      & ! inout
        & delta_c_acid_bg2_rp_sinkT(:),            & ! inout
        & delta_c_water_bg2_rp_sinkT(:),           & ! inout
        & delta_c_ethanol_bg2_rp_sinkT(:),         & ! inout
        & delta_c_nonsoluble_bg2_rp_sinkT(:),      & ! inout
        & delta_c_humus_1_rp_sinkT(:),             & ! inout
        & delta_c_humus_2_rp_sinkT(:),             & ! inout
        & delta_cflux_fire_all_2_atm_rp_sinkT(:)   & ! inout
        & )

      ! Add those parts of the ag litter carbon that remain after the fire
      delta_c_acid_ag1_rp_sinkT = delta_c_acid_ag1_rp_sinkT + d_acid_ag1_new
      delta_c_water_ag1_rp_sinkT = delta_c_water_ag1_rp_sinkT + d_water_ag1_new
      delta_c_ethanol_ag1_rp_sinkT = delta_c_ethanol_ag1_rp_sinkT + d_ethanol_ag1_new
      delta_c_nonsoluble_ag1_rp_sinkT = delta_c_nonsoluble_ag1_rp_sinkT + d_nonsoluble_ag1_new
      delta_c_acid_ag2_rp_sinkT = delta_c_acid_ag2_rp_sinkT + d_acid_ag2_new
      delta_c_water_ag2_rp_sinkT = delta_c_water_ag2_rp_sinkT + d_water_ag2_new
      delta_c_ethanol_ag2_rp_sinkT = delta_c_ethanol_ag2_rp_sinkT + d_ethanol_ag2_new
      delta_c_nonsoluble_ag2_rp_sinkT = delta_c_nonsoluble_ag2_rp_sinkT + d_nonsoluble_ag2_new

    ELSE
      WRITE (message_text,*) 'Unknown disturbance kind (', distKind, ')'
      CALL finish(TRIM(routine), message_text)
    ENDIF

  END SUBROUTINE recalc_delta_carbon_rp_upon_disturbance


  ! ================================================================================================================================
  !>
  !! @brief Calculates the per tile states from current canopy states
  !!
  !! @todo discuss how to make more general, e.g. get the info on which state variables to work on from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE recalc_per_tile_vars(tile, options, &
      & c_green , c_woods, c_reserve, c_crop_harvest, &
      & c_acid_ag1, c_water_ag1, c_ethanol_ag1, c_nonsoluble_ag1, &
      & c_acid_bg1, c_water_bg1, c_ethanol_bg1, c_nonsoluble_bg1, &
      & c_acid_ag2, c_water_ag2, c_ethanol_ag2, c_nonsoluble_ag2, &
      & c_acid_bg2, c_water_bg2, c_ethanol_bg2, c_nonsoluble_bg2, &
      & c_humus_1, c_humus_2, &
      & soil_respiration, NPP_flux_correction, NPP_pot_yDayMean, NPP_act_yDayMean, &
      & GPP_yDayMean, root_exudates, cflux_c_green_2_herb, cflux_herb_2_atm, &
      & cflux_dist_greenreserve_2_soil, cflux_dist_woods_2_soil, cflux_fire_all_2_atm)
!cflux_c_greenwood_2_litter

    USE mo_carbon_process,    ONLY: get_per_tile
    USE mo_carbon_constants,  ONLY: molarMassCO2_kg

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)            :: options
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_green(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_woods(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_reserve(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_crop_harvest(:)

    REAL(wp),                   INTENT(in), OPTIONAL  :: c_acid_ag1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_water_ag1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_ethanol_ag1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_nonsoluble_ag1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_acid_bg1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_water_bg1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_ethanol_bg1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_nonsoluble_bg1(:)

    REAL(wp),                   INTENT(in), OPTIONAL  :: c_acid_ag2(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_water_ag2(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_ethanol_ag2(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_nonsoluble_ag2(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_acid_bg2(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_water_bg2(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_ethanol_bg2(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_nonsoluble_bg2(:)

    REAL(wp),                   INTENT(in), OPTIONAL  :: c_humus_1(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_humus_2(:)

    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_dist_greenreserve_2_soil(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_dist_woods_2_soil(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_fire_all_2_atm(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: soil_respiration(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: NPP_pot_yDayMean(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: NPP_act_yDayMean(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: NPP_flux_correction(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: GPP_yDayMean(:)
!    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_c_greenwood_2_litter(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: root_exudates(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_c_green_2_herb(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_herb_2_atm(:)

    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(PHENO_)

    ! Local variables
    CHARACTER(len=*), PARAMETER     :: routine = TRIM(modname)//':recalc_per_tile_vars'
    INTEGER                         :: ics, ice, iblk

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

    dsl4jsb_Real2D_onChunk ::  cflux_dist_greenreserve_2_soil_ta  ! out
    dsl4jsb_Real2D_onChunk ::  cflux_dist_woods_2_soil_ta   ! out
    dsl4jsb_Real2D_onChunk ::  cflux_fire_all_2_atm_ta  ! out
    dsl4jsb_Real2D_onChunk ::  soil_respiration_ta      ! out
    dsl4jsb_Real2D_onChunk ::  NPP_pot_yDayMean_ta      ! out
    dsl4jsb_Real2D_onChunk ::  NPP_act_yDayMean_ta      ! out
    dsl4jsb_Real2D_onChunk ::  NPP_flux_correction_ta   ! out
    dsl4jsb_Real2D_onChunk ::  GPP_yDayMean_ta          ! out
!    dsl4jsb_Real2D_onChunk ::  cflux_c_greenwood_2_litter_ta           ! out
    dsl4jsb_Real2D_onChunk ::  root_exudates_ta         ! out
    dsl4jsb_Real2D_onChunk ::  cflux_c_green_2_herb_ta       ! out

    dsl4jsb_Real2D_onChunk ::  C_sum_natural_ta

    dsl4jsb_Real2D_onChunk ::  veg_fract_correction
    dsl4jsb_Real2D_onChunk ::  fract_fpc_max

    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(PHENO_)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green_ta )               ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve_ta )             ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods_ta )               ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest_ta)         ! opt out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1_ta)            ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1_ta)           ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1_ta)         ! opt out 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1_ta)      ! opt out 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1_ta)            ! opt out 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1_ta)           ! opt out 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1_ta)         ! opt out 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1_ta)      ! opt out 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1_ta)             ! opt out 

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2_ta)           ! opt out 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2_ta)          ! opt out 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2_ta)        ! opt out 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2_ta)     ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2_ta)           ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2_ta)          ! opt out 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2_ta)        ! opt out 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2_ta)     ! opt out 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2_ta)            ! opt out 

    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_dist_greenreserve_2_soil_ta ) ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_dist_woods_2_soil_ta )  ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_fire_all_2_atm_ta )      ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  soil_respiration_ta )          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_yDayMean_ta )          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_act_yDayMean_ta )          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_flux_correction_ta)        ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_yDayMean_ta )              ! opt out
!    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_greenwood_2_litter_ta )               ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  root_exudates_ta )             ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_green_2_herb_ta )           ! opt out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  C_sum_natural_ta )           ! out

    dsl4jsb_Get_var2D_onChunk(PHENO_,   veg_fract_correction )         ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_,   fract_fpc_max )               ! in

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

    CALL calculate_current_c_state_sum(tile, options, C_sum_natural_ta(:))

    ! C fluxes
    IF (PRESENT(cflux_dist_greenreserve_2_soil)) CALL get_per_tile(cflux_dist_greenreserve_2_soil_ta(:), &
      & cflux_dist_greenreserve_2_soil(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(cflux_dist_woods_2_soil)) CALL get_per_tile(cflux_dist_woods_2_soil_ta(:), &
      & cflux_dist_woods_2_soil(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(cflux_fire_all_2_atm))CALL get_per_tile(cflux_fire_all_2_atm_ta(:), &
      & cflux_fire_all_2_atm(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(NPP_pot_yDayMean)) CALL get_per_tile(NPP_pot_yDayMean_ta(:), &
      & NPP_pot_yDayMean(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(NPP_act_yDayMean)) CALL get_per_tile(NPP_act_yDayMean_ta(:), &
      & NPP_act_yDayMean(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(NPP_flux_correction)) CALL get_per_tile(NPP_flux_correction_ta(:), &
      & NPP_flux_correction(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(GPP_yDayMean)) CALL get_per_tile(GPP_yDayMean_ta(:), &
      & GPP_yDayMean(:), veg_fract_correction(:), fract_fpc_max(:))
!    cflux_c_greenwood_2_litter_ta               = cflux_c_greenwood_2_litter * veg_fract_correction * fract_fpc_max
    IF (PRESENT(root_exudates)) CALL get_per_tile(root_exudates_ta(:), &
      & root_exudates(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(cflux_c_green_2_herb)) CALL get_per_tile(cflux_c_green_2_herb_ta(:), &
      & cflux_c_green_2_herb(:), veg_fract_correction(:), fract_fpc_max(:))

  END SUBROUTINE recalc_per_tile_vars


  ! ================================================================================================================================
  !>
  !! @brief Adds delta carbon to carbon state variables and sets delta to zero
  !!
  !! @todo discuss how to make more general, e.g. get the info on which state variables to work on from somewhere else?
  !! TODO-JN-COHORTS: discuss how c_crop_harvest should be dealt with here - in my opinion it should anyway not be a density!
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE apply_delta_carbon_to_carbon_states(tile, options)

    USE mo_carbon_process,      ONLY: update_state_reset_delta

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)            :: options

    TYPE(t_jsb_model), POINTER                :: model

    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_config(PHENO_)

    ! Local variables
    CHARACTER(len=*), PARAMETER     :: routine = TRIM(modname)//':apply_delta_carbon_to_carbon_states'
    INTEGER                         :: ics, ice, iblk, nc, i
    REAL(wp), DIMENSION(options%nc) :: delta_fract, prev_fract

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

    dsl4jsb_Real2D_onChunk ::  GPP_sum
    dsl4jsb_Real2D_onChunk ::  NPP_pot_sum

    dsl4jsb_Real2D_onChunk ::  cflux_dist_greenreserve_2_soil
    dsl4jsb_Real2D_onChunk ::  cflux_dist_woods_2_soil
    dsl4jsb_Real2D_onChunk ::  cflux_fire_all_2_atm

    dsl4jsb_Real2D_onChunk ::  current_max_green
    dsl4jsb_Real2D_onChunk ::  veg_carbon_at_max_green

    dsl4jsb_Real2D_onChunk ::  day_npp_sum

    dsl4jsb_Real2D_onChunk ::  annual_avg_NPP
    dsl4jsb_Real2D_onChunk ::  biennial_avg_NPP
    dsl4jsb_Real2D_onChunk ::  decennial_avg_NPP

    dsl4jsb_Real2D_onChunk ::  delta_c_green_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_reserve_rp 
    dsl4jsb_Real2D_onChunk ::  delta_c_woods_rp 
    dsl4jsb_Real2D_onChunk ::  delta_c_crop_harvest_rp

    dsl4jsb_Real2D_onChunk ::  delta_c_acid_ag1_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_water_ag1_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_ag1_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_ag1_rp 
    dsl4jsb_Real2D_onChunk ::  delta_c_acid_bg1_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_water_bg1_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_bg1_rp 
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_bg1_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_humus_1_rp 

    dsl4jsb_Real2D_onChunk ::  delta_c_acid_ag2_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_water_ag2_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_ag2_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_ag2_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_acid_bg2_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_water_bg2_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_ethanol_bg2_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_nonsoluble_bg2_rp
    dsl4jsb_Real2D_onChunk ::  delta_c_humus_2_rp

    dsl4jsb_Real2D_onChunk ::  delta_GPP_sum_rp
    dsl4jsb_Real2D_onChunk ::  delta_NPP_pot_sum_rp
    
    dsl4jsb_Real2D_onChunk ::  delta_cflux_dist_greenreserve_2_soil_rp
    dsl4jsb_Real2D_onChunk ::  delta_cflux_dist_woods_2_soil_rp

    dsl4jsb_Real2D_onChunk ::  delta_current_max_green_rp
    dsl4jsb_Real2D_onChunk ::  delta_veg_carbon_at_max_green_rp

    dsl4jsb_Real2D_onChunk ::  delta_day_npp_sum_rp

    dsl4jsb_Real2D_onChunk ::  delta_annual_avg_NPP_rp
    dsl4jsb_Real2D_onChunk ::  delta_biennial_avg_NPP_rp
    dsl4jsb_Real2D_onChunk ::  delta_decennial_avg_NPP_rp

    dsl4jsb_Real2D_onChunk ::  delta_cflux_fire_all_2_atm_rp

    dsl4jsb_Real2D_onChunk ::  veg_fract_correction

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_config(PHENO_)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest)         ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1)            ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1)         ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1)      ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1)            ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1)           ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1)         ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1)      ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1)             ! inout 

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2)           ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2)          ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2)        ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2)          ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2)        ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2)     ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2)            ! inout 

    dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_sum)                ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_sum)            ! inout 

    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_dist_greenreserve_2_soil) ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_dist_woods_2_soil)  ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_fire_all_2_atm)      ! inout 

    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_green_rp )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_reserve_rp )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_woods_rp )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_crop_harvest_rp)         ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_acid_ag1_rp)            ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_water_ag1_rp)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_ethanol_ag1_rp)         ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_nonsoluble_ag1_rp)      ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_acid_bg1_rp)            ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_water_bg1_rp)           ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_ethanol_bg1_rp)         ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_nonsoluble_bg1_rp)      ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_humus_1_rp)             ! inout 

    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_acid_ag2_rp)           ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_water_ag2_rp)          ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_ethanol_ag2_rp)        ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_nonsoluble_ag2_rp)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_acid_bg2_rp)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_water_bg2_rp)          ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_ethanol_bg2_rp)        ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_nonsoluble_bg2_rp)     ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_c_humus_2_rp)            ! inout 

    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_GPP_sum_rp)                ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_NPP_pot_sum_rp)            ! inout 

    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_cflux_dist_greenreserve_2_soil_rp) ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_cflux_dist_woods_2_soil_rp)  ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_cflux_fire_all_2_atm_rp)      ! inout 

    dsl4jsb_Get_var2D_onChunk(PHENO_,  day_npp_sum)          ! inout 
    dsl4jsb_Get_var2D_onChunk(PHENO_,  delta_day_npp_sum_rp) ! inout 

    IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
      dsl4jsb_Get_var2D_onChunk(CARBON_,  current_max_green)                ! inout 
      dsl4jsb_Get_var2D_onChunk(CARBON_,  veg_carbon_at_max_green)          ! inout 

      dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_current_max_green_rp)       ! inout 
      dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_veg_carbon_at_max_green_rp) ! inout 

      dsl4jsb_Get_var2D_onChunk(CARBON_,  annual_avg_NPP)                 ! inout 
      dsl4jsb_Get_var2D_onChunk(CARBON_,  biennial_avg_NPP)           ! inout 
      dsl4jsb_Get_var2D_onChunk(CARBON_,  decennial_avg_NPP)          ! inout 

      dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_annual_avg_NPP_rp)        ! inout 
      dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_biennial_avg_NPP_rp)  ! inout 
      dsl4jsb_Get_var2D_onChunk(CARBON_,  delta_decennial_avg_NPP_rp) ! inout 

    ENDIF

    dsl4jsb_Get_var2D_onChunk(PHENO_,   veg_fract_correction)      ! in

    delta_fract(:)=tile%delta_fract(ics:ice, iblk)
    CALL tile%Get_fraction(ics, ice, iblk, fract=prev_fract(:))

    !Assertion prev_fract + delta_fract .GE. 0._wp
    IF (ANY(prev_fract(:) + delta_fract(:) .LT. 0_wp )) THEN
      DO i = 1,nc
        IF (prev_fract(i) + delta_fract(i) .LT. 0_wp ) THEN
          WRITE (message_text,*) 'Violation of assertion: area to be changed (', delta_fract(i), &
            & ') needs to be .LE. existing area (', prev_fract(i) ,')'
          CALL finish(TRIM(routine), message_text)
        ENDIF
      ENDDO
    ENDIF

    ! TODO-JN-COHORTS: if the smallest possible fraction 1e-10 is re-introduced
    !      then the following should not happen and other Assertions and MAX statements will be required!)

    CALL update_state_reset_delta(c_green(:), delta_c_green_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_reserve(:), delta_c_reserve_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_woods(:), delta_c_woods_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_crop_harvest(:), delta_c_crop_harvest_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))

    CALL update_state_reset_delta(c_acid_ag1(:), delta_c_acid_ag1_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_water_ag1(:), delta_c_water_ag1_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_ethanol_ag1(:), delta_c_ethanol_ag1_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_nonsoluble_ag1(:), delta_c_nonsoluble_ag1_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_acid_bg1(:), delta_c_acid_bg1_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_water_bg1(:), delta_c_water_bg1_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_ethanol_bg1(:), delta_c_ethanol_bg1_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_nonsoluble_bg1(:), delta_c_nonsoluble_bg1_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_humus_1(:), delta_c_humus_1_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))

    CALL update_state_reset_delta(c_acid_ag2(:), delta_c_acid_ag2_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_water_ag2(:), delta_c_water_ag2_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_ethanol_ag2(:), delta_c_ethanol_ag2_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_nonsoluble_ag2(:), delta_c_nonsoluble_ag2_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_acid_bg2(:), delta_c_acid_bg2_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_water_bg2(:), delta_c_water_bg2_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_ethanol_bg2(:), delta_c_ethanol_bg2_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_nonsoluble_bg2(:), delta_c_nonsoluble_bg2_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(c_humus_2(:), delta_c_humus_2_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))

    CALL update_state_reset_delta(GPP_sum(:), delta_GPP_sum_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(NPP_pot_sum(:), delta_NPP_pot_sum_rp(:), delta_fract(:), &
      & veg_fract_correction(:), prev_fract(:))

    ! Following three are only required if with disturbances or harvest, however, they are always defined anyway 
    CALL update_state_reset_delta(cflux_dist_greenreserve_2_soil(:), delta_cflux_dist_greenreserve_2_soil_rp(:), &
      & delta_fract(:), veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(cflux_dist_woods_2_soil(:), delta_cflux_dist_woods_2_soil_rp(:), &
      & delta_fract(:), veg_fract_correction(:), prev_fract(:))
    CALL update_state_reset_delta(cflux_fire_all_2_atm(:), delta_cflux_fire_all_2_atm_rp(:), &
      & delta_fract(:), veg_fract_correction(:), prev_fract(:))

    IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
      CALL update_state_reset_delta(current_max_green(:), delta_current_max_green_rp(:), delta_fract(:), &
        & veg_fract_correction(:), prev_fract(:))
      CALL update_state_reset_delta(veg_carbon_at_max_green(:), delta_veg_carbon_at_max_green_rp(:), delta_fract(:), &
        & veg_fract_correction(:), prev_fract(:))

      CALL update_state_reset_delta(annual_avg_NPP(:), delta_annual_avg_NPP_rp(:), delta_fract(:), &
        & veg_fract_correction(:), prev_fract(:))
      CALL update_state_reset_delta(biennial_avg_NPP(:), delta_biennial_avg_NPP_rp(:), delta_fract(:), &
        & veg_fract_correction(:), prev_fract(:))
      CALL update_state_reset_delta(decennial_avg_NPP(:), delta_decennial_avg_NPP_rp(:), delta_fract(:), &
        & veg_fract_correction(:), prev_fract(:))

    ENDIF

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
      & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:),       &
      & cflux_dist_greenreserve_2_soil = cflux_dist_greenreserve_2_soil(:), &
      & cflux_dist_woods_2_soil = cflux_dist_woods_2_soil(:), cflux_fire_all_2_atm = cflux_fire_all_2_atm(:))

  END SUBROUTINE apply_delta_carbon_to_carbon_states


  ! ================================================================================================================================
  !>
  !! @brief calculates the current sum of all carbon state variables -> on ta! (per canopy vars are not aggregated!)
  !!
  !! @todo discuss how to make more general, e.g. get the info which are the c state variables from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  ! 
  SUBROUTINE calculate_current_c_state_sum(tile, options, current_c_state_sum)

    USE mo_util,                ONLY: real2string

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)         :: options
    REAL(wp),                   INTENT(inout)      :: current_c_state_sum(:)

    dsl4jsb_Def_memory(CARBON_)

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calculate_current_c_state_sum'
    INTEGER                     :: ics, ice, iblk

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

    dsl4jsb_Get_memory(CARBON_)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

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

    !TODO-JN-COHORTS: could be done in mo_carbon_process
    current_c_state_sum = c_green_ta + c_woods_ta + c_reserve_ta + c_crop_harvest_ta + &
        & c_acid_ag1_ta + c_water_ag1_ta + c_ethanol_ag1_ta + c_nonsoluble_ag1_ta + &
        & c_acid_ag2_ta + c_water_ag2_ta + c_ethanol_ag2_ta + c_nonsoluble_ag2_ta + &
        & c_acid_bg1_ta + c_water_bg1_ta + c_ethanol_bg1_ta + c_nonsoluble_bg1_ta + &
        & c_acid_bg2_ta + c_water_bg2_ta + c_ethanol_bg2_ta + c_nonsoluble_bg2_ta + &
        & c_humus_1_ta  + c_humus_2_ta

  END SUBROUTINE calculate_current_c_state_sum


  ! ================================================================================================================================
  !>
  !! @brief conducts a carbon conservation test on current states and given old sum + fluxes
  !!        -> on ta! (per canopy vars are not aggregated!) 
  !!
  !! @todo discuss how to make more general, e.g. get the info which are the c state variables from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  ! 
  SUBROUTINE check_carbon_conservation_for_operation(tile, options, old_c_state_sum, cflux, Cconserve)

    USE mo_util,                ONLY: real2string

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)         :: options
    REAL(wp),                   INTENT(in)         :: old_c_state_sum(:)
    REAL(wp),                   INTENT(in)         :: cflux(:)
    REAL(wp),                   INTENT(inout)      :: Cconserve(:)

    dsl4jsb_Def_memory(CARBON_)

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':check_carbon_conservation_for_operation'
    INTEGER                     :: ics, ice, iblk

    CALL calculate_current_c_state_sum(tile, options, Cconserve)
    !TODO-JN-COHORTS: could be done in mo_carbon_process
    Cconserve = Cconserve - old_c_state_sum - cflux

  END SUBROUTINE check_carbon_conservation_for_operation

#endif
END MODULE mo_carbon_interface
