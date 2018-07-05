!>
!! @brief Contains the interfaces to the carbon process
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-04-05)
!!
MODULE mo_carbon_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module (Get_carbon_memory and Get_carbon_config)
  dsl4jsb_Use_processes CARBON_, ASSIMI_, HYDRO_, PHENO_, RAD_, A2L_, DISTURB_

  ! Use of process configurations (t_carbon_config)
  dsl4jsb_Use_config(ASSIMI_)
  dsl4jsb_Use_config(CARBON_)
  dsl4jsb_Use_config(DISTURB_)

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
  PUBLIC ::  Register_carbon_tasks

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


  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for carbon process
  !!
  !! @param[in,out] this      Instance of carbon process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_carbon_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_C_NPP_pot_allocation(model_id))

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
    USE mo_carbon_constants,  ONLY: molarMassCO2_kg

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
    CHARACTER(len=*), PARAMETER  :: routine = modname//':update_C_NPP_pot_allocation'

    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_memory(ASSIMI_)
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
    dsl4jsb_Real2D_onChunk ::  cflux_c_greenwood_2_litter           ! Carbon flux from the veget. to the litter pools
                                                     ! [mol(C)/m^2(canopy) s]
    dsl4jsb_Real2D_onChunk ::  c_green           ! C-pool for leaves, fine roots, vegetative organs and
                                                     ! other green (living) parts of vegetation [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_reserve         ! C-pool for carbohydrate reserve (sugars, starches) that
                                                     ! allows plants to survive bad times[mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_woods           ! C-pool for stems, thick roots and other (dead) structural
                                                     !  material of living plants [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_crop_harvest    ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]
    dsl4jsb_Real2D_onChunk ::  c_litter_green_ag ! C-pool for above ground non-woody litter especially from leaves
                                                     !     [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_litter_green_bg ! C-pool for below ground non-woody litter especially from
                                                     !     fine roots [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_litter_wood_ag  ! C-pool for litter originating from above ground woody parts of the
                                                     !    plants [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_litter_wood_bg  ! C-pool for litter originating from below ground woody parts of the
                                                     !    plants [mol(C)/m^2(canopy)]
    !CBAL dsl4jsb_Real2D_onChunk ::  c_slow            ! C-pool for below ground organic material (coming from the litter
                                                     !     pools) [mol(C)/m^2(canopy)]
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
    ! TA
    dsl4jsb_Real2D_onChunk ::  C_sum_natural_ta
    dsl4jsb_Real2D_onChunk ::  c_green_ta        ! C-pool for leaves, fine roots, vegetative organs and other green
                                                     ! .. (living) parts of vegetation [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_reserve_ta      ! C-pool for carbohydrate reserve (sugars, starches) that allows plants
                                                     ! .. to survive bad times[mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_woods_ta        ! C-pool for stems, thick roots and other (dead) structural material of
                                                     ! .. living plants [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_crop_harvest_ta ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]     
    dsl4jsb_Real2D_onChunk ::  c_litter_green_ag_ta ! C-pool for above ground non-woody litter especially from leaves
                                                     !     [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_litter_green_bg_ta ! C-pool for below ground non-woody litter especially from fine roots
                                                     !     [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_litter_wood_ag_ta  ! C-pool for litter originating from above ground woody parts of the
                                                     !    plants [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_litter_wood_bg_ta  ! C-pool for litter originating from below ground woody parts of the
                                                     !    plants [mol(C)/m^2(canopy)]
    !CBAL dsl4jsb_Real2D_onChunk ::  c_slow_ta         ! C-pool for below ground organic material (coming from the litter
                                                     !     pools) [mol(C)/m^2(canopy)]

    dsl4jsb_Real2D_onChunk ::  c_bg_sum          ! Sum of belowground Yasso pools  [mol(C)^2/m(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ag_sum_1        ! Sum of aboveground green Yasso pools  [mol(C)^2/m(canopy)]

    !SIZE CLASS 1
    dsl4jsb_Real2D_onChunk ::  c_acid_ag1_ta     ! Yasso above ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_ag1_ta    ! Yasso above ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag1_ta  ! Yasso above ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag1_ta  ! Yasso above ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_acid_bg1_ta     ! Yasso below ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_bg1_ta    ! Yasso below ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg1_ta  ! Yasso below ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg1_ta  ! Yasso below ground litter-pool for non-soluble litter
                                                      ! [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_humus_1_ta      ! Yasso below ground litter-pool for slow C compartment
                                                      !  [mol(C)^2/m(canopy)]
    !SIZE CLASS 2
    dsl4jsb_Real2D_onChunk ::  c_acid_ag2_ta     ! Yasso above ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_ag2_ta    ! Yasso above ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag2_ta  ! Yasso above ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag2_ta  ! Yasso above ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_acid_bg2_ta     ! Yasso below ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_water_bg2_ta    ! Yasso below ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg2_ta  ! Yasso below ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg2_ta  ! Yasso below ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_humus_2_ta      ! Yasso below ground litter-pool for slow C compartment
                                                      !  [mol(C)^2/m(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_bg_sum_ta          ! Sum of belowground Yasso pools  [mol(C)^2/m(canopy)]
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
    dsl4jsb_Real2D_onChunk ::  cflux_c_greenwood_2_litter_ta           ! out
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
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_greenwood_2_litter )               ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_litter_green_ag )     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_litter_green_bg )     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_litter_wood_ag  )     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_litter_wood_bg  )     ! inout
    !CBAL dsl4jsb_Get_var2D_onChunk(CARBON_,  c_slow            )     ! inout
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

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_bg_sum)             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ag_sum_1)           ! inout

    !
    ! This is on tile area:
    dsl4jsb_Get_var2D_onChunk(CARBON_,  C_sum_natural_ta)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green_ta )           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve_ta )         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods_ta )           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest_ta)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_litter_green_ag_ta ) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_litter_green_bg_ta ) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_litter_wood_ag_ta  ) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_litter_wood_bg_ta  ) ! inout
    !CBAL dsl4jsb_Get_var2D_onChunk(CARBON_,  c_slow_ta            ) ! inout
    !SIZE CLASS 1
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1_ta)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1_ta)       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1_ta)     ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1_ta)  ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1_ta)        ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1_ta)       ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1_ta)     ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1_ta)  ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1_ta)         ! inout 
    !SIZE CLASS 2
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2_ta)       ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2_ta)      ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2_ta)    ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2_ta) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2_ta)       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2_ta)      ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2_ta)    ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2_ta) ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2_ta)        ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_bg_sum_ta)             ! inout

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
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_green_2_herb )           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_herb_2_littergreen  )       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_herb_2_atm )     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_npp_2_atm_ta )         ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_soilresp_2_atm_ta )    ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_herb_2_atm_ta )   ! out

    ! This is on tile area:
    dsl4jsb_Get_var2D_onChunk(CARBON_,  soil_respiration_ta )          ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_yDayMean_ta )          ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_act_yDayMean_ta )          ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_flux_correction_ta)        ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_yDayMean_ta )              ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_greenwood_2_litter_ta )               ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  root_exudates_ta )             ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_green_2_herb_ta )           ! out

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

      ! Now update the C pools
      CALL calc_Cpools( &
        &  LAI_yDayMean,                                      & ! in
        &  LAI_yyDayMean,                                     & ! in
        &  NPP_pot_yDayMean,                                  & ! in
        &  topSoilTemp_yDayMean,                              & ! in
        &  alpha_yDayMean,                                    & ! in
        !
        &  dsl4jsb_Lctlib_param(fract_npp_2_woodPool),        & ! in
        &  dsl4jsb_Lctlib_param(fract_NPP_2_reservePool),     & ! in
        &  dsl4jsb_Lctlib_param(fract_NPP_2_exudates),        & ! in
        !
        &  dsl4jsb_Lctlib_param(fract_green_2_herbivory),     & ! in
        &  dsl4jsb_Lctlib_param(tau_c_litter_leaf),           & ! in
        &  dsl4jsb_Lctlib_param(tau_c_litter_wood),           & ! in
        &  dsl4jsb_Lctlib_param(tau_c_woods),                 & ! in
        &  dsl4jsb_Lctlib_param(LAI_shed_constant),           & ! in
        &  dsl4jsb_Lctlib_param(fract_C_litter_green2atmos),  & ! in
        &  dsl4jsb_Lctlib_param(Max_C_content_woods),         & ! in
        &  dsl4jsb_Lctlib_param(specificLeafArea_C),          & ! in
        &  dsl4jsb_Lctlib_param(reserveC2leafC),              & ! in
        &  dsl4jsb_Lctlib_param(MaxLAI),                      & ! in
        &  dsl4jsb_Lctlib_param(PhenologyType),               & ! in  R: is it crop? If yes then this is 5!
        &  time_steps_per_day(),                              & ! in
        &  dtime,                                             & ! in
        !
        &  c_green,                                           & ! inout
        &  c_woods,                                           & ! inout
        &  c_reserve,                                         & ! inout
        &  c_crop_harvest,                                    & ! inout
        &  c_litter_green_ag,                                 & ! inout
        &  c_litter_green_bg,                                 & ! inout
        &  c_litter_wood_ag,                                  & ! inout
        &  c_litter_wood_bg,                                  & ! inout
        ! &  CBAL c_slow,                                       & ! inout
        !
        &  soil_respiration,                                  & ! out
        &  NPP_flux_correction,                               & ! out
        &  excess_NPP,                                        & ! out
        &  root_exudates,                                     & ! out
        &  cflux_c_greenwood_2_litter,                        & ! out
        &  cflux_c_green_2_herb,                              & ! out
        &  cflux_herb_2_littergreen,                          & ! out
        &  cflux_herb_2_atm,                                  & ! out
        &  NPP_act_yDayMean,                                  & ! out
        !
        !&  fract_litter_wood_new    = fract_litter_wood_new,       & ! inout, optional. Only for spitfire
        !
        ! variables only needed with yasso:
        ! R: Note, in the moment the procedure is only used with Yasso.
        &  temp2_30d        =  pseudo_temp_yDay,                      & ! in, optional, necessary with yasso
        &  precip_30d       =  pseudo_precip_yDay,                    & ! in, optional, necessary with yasso
                             !
        &  c_acid_ag1       =  c_acid_ag1,                            & ! inout, optional, necessary with yasso
        &  c_water_ag1      =  c_water_ag1,                           & ! inout, optional, necessary with yasso
        &  c_ethanol_ag1    =  c_ethanol_ag1,                         & ! inout, optional, necessary with yasso
        &  c_nonsoluble_ag1 =  c_nonsoluble_ag1,                      & ! inout, optional, necessary with yasso
        &  c_acid_bg1       =  c_acid_bg1,                            & ! inout, optional, necessary with yasso
        &  c_water_bg1      =  c_water_bg1,                           & ! inout, optional, necessary with yasso
        &  c_ethanol_bg1    =  c_ethanol_bg1,                         & ! inout, optional, necessary with yasso
        &  c_nonsoluble_bg1 =  c_nonsoluble_bg1,                      & ! inout, optional, necessary with yasso
        &  c_humus_1        =  c_humus_1,                             & ! inout, optional, necessary with yasso
        &  c_acid_ag2       =  c_acid_ag2,                            & ! inout, optional, necessary with yasso
        &  c_water_ag2      =  c_water_ag2,                           & ! inout, optional, necessary with yasso
        &  c_ethanol_ag2    =  c_ethanol_ag2,                         & ! inout, optional, necessary with yasso
        &  c_nonsoluble_ag2 =  c_nonsoluble_ag2,                      & ! inout, optional, necessary with yasso
        &  c_acid_bg2       =  c_acid_bg2,                            & ! inout, optional, necessary with yasso
        &  c_water_bg2      =  c_water_bg2,                           & ! inout, optional, necessary with yasso
        &  c_ethanol_bg2    =  c_ethanol_bg2,                         & ! inout, optional, necessary with yasso
        &  c_nonsoluble_bg2 =  c_nonsoluble_bg2,                      & ! inout, optional, necessary with yasso
        &  c_humus_2        =  c_humus_2,                             & ! inout, optional, necessary with yasso
                              !
        &  LeafLit_coef1    =  dsl4jsb_Lctlib_param(LeafLit_coef(1)), & ! in, optional, necessary with yasso
        &  LeafLit_coef2    =  dsl4jsb_Lctlib_param(LeafLit_coef(2)), & ! in, optional, necessary with yasso
        &  LeafLit_coef3    =  dsl4jsb_Lctlib_param(LeafLit_coef(3)), & ! in, optional, necessary with yasso
        &  LeafLit_coef4    =  dsl4jsb_Lctlib_param(LeafLit_coef(4)), & ! in, optional, necessary with yasso
        &  LeafLit_coef5    =  dsl4jsb_Lctlib_param(LeafLit_coef(5)), & ! in, optional, necessary with yasso
        &  WoodLit_coef1    =  dsl4jsb_Lctlib_param(WoodLit_coef(1)), & ! in, optional, necessary with yasso
        &  WoodLit_coef2    =  dsl4jsb_Lctlib_param(WoodLit_coef(2)), & ! in, optional, necessary with yasso
        &  WoodLit_coef3    =  dsl4jsb_Lctlib_param(WoodLit_coef(3)), & ! in, optional, necessary with yasso
        &  WoodLit_coef4    =  dsl4jsb_Lctlib_param(WoodLit_coef(4)), & ! in, optional, necessary with yasso
        &  WoodLit_coef5    =  dsl4jsb_Lctlib_param(WoodLit_coef(5)), & ! in, optional, necessary with yasso
        &  WoodLitterSize   =  dsl4jsb_Lctlib_param(WoodLitterSize)   & ! in, optional, necessary with yasso
        & )

      c_bg_sum =   c_acid_bg1 + c_water_bg1 + c_ethanol_bg1 + c_nonsoluble_bg1 + c_humus_1 &  ! Size class 1; green litter
        &        + c_acid_bg2 + c_water_bg2 + c_ethanol_bg2 + c_nonsoluble_bg2 + c_humus_2    ! Size class 2; woody litter

      c_ag_sum_1 =  c_acid_ag1 + c_water_ag1 + c_ethanol_ag1 + c_nonsoluble_ag1 + c_humus_1

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
      c_green_ta           = c_green           * veg_fract_correction * fract_fpc_max
      c_woods_ta           = c_woods           * veg_fract_correction * fract_fpc_max
      c_reserve_ta         = c_reserve         * veg_fract_correction * fract_fpc_max
      c_crop_harvest_ta    = c_crop_harvest    * veg_fract_correction * fract_fpc_max

      c_acid_ag1_ta       = c_acid_ag1       * veg_fract_correction * fract_fpc_max
      c_water_ag1_ta      = c_water_ag1      * veg_fract_correction * fract_fpc_max
      c_ethanol_ag1_ta    = c_ethanol_ag1    * veg_fract_correction * fract_fpc_max
      c_nonsoluble_ag1_ta = c_nonsoluble_ag1 * veg_fract_correction * fract_fpc_max
      c_acid_bg1_ta       = c_acid_bg1       * veg_fract_correction * fract_fpc_max
      c_water_bg1_ta      = c_water_bg1      * veg_fract_correction * fract_fpc_max
      c_ethanol_bg1_ta    = c_ethanol_bg1    * veg_fract_correction * fract_fpc_max
      c_nonsoluble_bg1_ta = c_nonsoluble_bg1 * veg_fract_correction * fract_fpc_max
      c_humus_1_ta        = c_humus_1        * veg_fract_correction * fract_fpc_max
      c_acid_ag2_ta       = c_acid_ag2       * veg_fract_correction * fract_fpc_max
      c_water_ag2_ta      = c_water_ag2      * veg_fract_correction * fract_fpc_max
      c_ethanol_ag2_ta    = c_ethanol_ag2    * veg_fract_correction * fract_fpc_max
      c_nonsoluble_ag2_ta = c_nonsoluble_ag2 * veg_fract_correction * fract_fpc_max
      c_acid_bg2_ta       = c_acid_bg2       * veg_fract_correction * fract_fpc_max
      c_water_bg2_ta      = c_water_bg2      * veg_fract_correction * fract_fpc_max
      c_ethanol_bg2_ta    = c_ethanol_bg2    * veg_fract_correction * fract_fpc_max
      c_nonsoluble_bg2_ta = c_nonsoluble_bg2 * veg_fract_correction * fract_fpc_max
      c_humus_2_ta        = c_humus_2        * veg_fract_correction * fract_fpc_max

      c_bg_sum_ta         = c_bg_sum         * veg_fract_correction * fract_fpc_max

      C_sum_natural_ta =  &
        &   c_green_ta          + c_woods_ta           &
        & + c_reserve_ta        + c_crop_harvest_ta    &
        & + c_acid_ag1_ta       + c_acid_ag2_ta        &
        & + c_water_ag1_ta      + c_water_ag2_ta       &
        & + c_ethanol_ag1_ta    + c_ethanol_ag2_ta     &
        & + c_nonsoluble_ag1_ta + c_nonsoluble_ag2_ta  &
        & + c_acid_bg1_ta       + c_acid_bg2_ta        &
        & + c_water_bg1_ta      + c_water_bg2_ta       &
        & + c_ethanol_bg1_ta    + c_ethanol_bg2_ta     &
        & + c_nonsoluble_bg1_ta + c_nonsoluble_bg2_ta  &
        & + c_humus_1_ta        + c_humus_2_ta

      ! C fluxes
      soil_respiration_ta           = soil_respiration           * veg_fract_correction * fract_fpc_max
      !soil_respiration_pot_ta       = soil_respiration_pot      * veg_fract_correction * fract_fpc_max
      NPP_pot_yDayMean_ta           = NPP_pot_yDayMean           * veg_fract_correction * fract_fpc_max
      NPP_act_yDayMean_ta           = NPP_act_yDayMean           * veg_fract_correction * fract_fpc_max
      NPP_flux_correction_ta        = NPP_flux_correction        * veg_fract_correction * fract_fpc_max
      GPP_yDayMean_ta               = GPP_yDayMean               * veg_fract_correction * fract_fpc_max
      cflux_c_greenwood_2_litter_ta = cflux_c_greenwood_2_litter * veg_fract_correction * fract_fpc_max
      root_exudates_ta              = root_exudates              * veg_fract_correction * fract_fpc_max
      cflux_c_green_2_herb_ta          = cflux_c_green_2_herb    * veg_fract_correction * fract_fpc_max
      ! cflux_herb_2_littergreen_ta    = cflux_herb_2_littergreen    * veg_fract_correction * fract_fpc_max
      ! cflux_herb_2_atm_ta   = cflux_herb_2_atm * veg_fract_correction * fract_fpc_max

    END IF

    ! Compute net CO2 fluxes exchanged with atmosphere at each time step
    !-------------------------------------------------------------------
    ! Note: carbon loss of biosphere means a positive CO2 flux to atmosphere (i.e. NEP and net CO2-flux have opposite signs)

    cflux_npp_2_atm_ta = molarMassCO2_kg *                      & ! Conversion factor from mol to kg CO2
               (- veg_fract_correction * fract_fpc_max          & ! Minus: atmosphere gain is positive
                    * (NPP_pot_rate_ca                          & ! current (not actual) NPP rate
                        - (NPP_pot_yDayMean-NPP_act_yDayMean)) )  ! corrected with yesterdays actual NPP defizit


    cflux_soilresp_2_atm_ta = molarMassCO2_kg *                 & ! Conversion factor from mol to kg CO2
               (- veg_fract_correction *fract_fpc_max           & ! Minus: atmosphere gain is positive
                    * soil_respiration      )                    ! .. soil respiration

    cflux_herb_2_atm_ta = molarMassCO2_kg *                     & ! Conversion factor from mol to kg CO2
               (- veg_fract_correction * fract_fpc_max          & ! Minus: atmosphere gain is positive
                    * cflux_herb_2_atm )                          ! .. herbivory

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

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_C_NPP_pot_allocation'

    INTEGER  :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(CARBON_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yDayMean,          weighted_by_fract)

    !@todo Those commented out are currently not calculated
    ! same for ta variables
    ! non-yasso-pools
    dsl4jsb_Aggregate_onChunk(CARBON_, c_green_ta,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_woods_ta,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_reserve_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_crop_harvest_ta,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_litter_green_ag_ta,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_litter_green_bg_ta,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_litter_wood_ag_ta,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_litter_wood_bg_ta,   weighted_by_fract)
    !CBAL dsl4jsb_Aggregate_onChunk(CARBON_, c_slow_ta,        weighted_by_fract)
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
    ! sums
    dsl4jsb_Aggregate_onChunk(CARBON_, c_bg_sum_ta,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, C_sum_natural_ta,      weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_pot_sum_ta,        weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, GPP_sum_ta,            weighted_by_fract)
    ! fluxes
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_c_greenwood_2_litter_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, soil_respiration_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_flux_correction_ta,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, excess_NPP_ta,            weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yyDayMean_ta,         weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yDayMean_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_pot_yDayMean_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_act_yDayMean_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, GPP_yDayMean_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, root_exudates_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_c_green_2_herb_ta,  weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_, cflux_herb_2_littergreen_ta, weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_, cflux_herb_2_atm_ta, weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_, fract_litter_wood_new_ta,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_npp_2_atm_ta,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_soilresp_2_atm_ta,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_herb_2_atm_ta,      weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_, cflux_npp_2_atm_test_ta,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_C_NPP_pot_allocation

#endif
END MODULE mo_carbon_interface
