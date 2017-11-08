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
  PUBLIC ::  Register_carbon_tasks , relocate_disturbed_carbon ! R: relocate_disturbed_carbon is not a carbon task but part of
                                                               !    the process disturbance only affecting the C pools.
                                                               !    It is called by tsk_natural_disturbances only if the
                                                               !    CARBON process is active. This logic follows JSBACH3, but
                                                               !    could also be done the other way round...

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

    return_ptr%name = 'C_NPP_pot_allocation'
    return_ptr%process_id = CARBON_
    return_ptr%owner_model_id = model_id

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
    CHARACTER(len=*), PARAMETER  :: routine = TRIM(modname)//':update_C_NPP_pot_allocation'


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
    dsl4jsb_Real2D_onChunk ::  litter_flux           ! Carbon flux from the veget. to the litter pools
                                                     ! [mol(C)/m^2(canopy) s]
    dsl4jsb_Real2D_onChunk ::  Cpool_green           ! C-pool for leaves, fine roots, vegetative organs and
                                                     ! other green (living) parts of vegetation [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  Cpool_reserve         ! C-pool for carbohydrate reserve (sugars, starches) that
                                                     ! allows plants to survive bad times[mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  Cpool_woods           ! C-pool for stems, thick roots and other (dead) structural
                                                     !  material of living plants [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  Cpool_crop_harvest    ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]
    dsl4jsb_Real2D_onChunk ::  Cpool_litter_green_ag ! C-pool for above ground non-woody litter especially from leaves
                                                     !     [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  Cpool_litter_green_bg ! C-pool for below ground non-woody litter especially from
                                                     !     fine roots [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  Cpool_litter_wood_ag  ! C-pool for litter originating from above ground woody parts of the
                                                     !    plants [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  Cpool_litter_wood_bg  ! C-pool for litter originating from below ground woody parts of the
                                                     !    plants [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  Cpool_slow            ! C-pool for below ground organic material (coming from the litter
                                                     !     pools) [mol(C)/m^2(canopy)]
    !SIZE CLASS 1
    dsl4jsb_Real2D_onChunk ::  YCpool_acid_ag1        ! Yasso above ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_water_ag1       ! Yasso above ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_ethanol_ag1     ! Yasso above ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_nonsoluble_ag1  ! Yasso above ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_acid_bg1        ! Yasso below ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_water_bg1       ! Yasso below ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_ethanol_bg1     ! Yasso below ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_nonsoluble_bg1  ! Yasso below ground litter-pool for non-soluble litter
                                                      ! [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_humus_1         ! Yasso below ground litter-pool for slow C compartment
                                                      !  [mol(C)^2/m(canopy)]
    !SIZE CLASS 2
    dsl4jsb_Real2D_onChunk ::  YCpool_acid_ag2        ! Yasso above ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_water_ag2       ! Yasso above ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_ethanol_ag2     ! Yasso above ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_nonsoluble_ag2  ! Yasso above ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_acid_bg2        ! Yasso below ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_water_bg2       ! Yasso below ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_ethanol_bg2     ! Yasso below ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_nonsoluble_bg2  ! Yasso below ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_humus_2         ! Yasso below ground litter-pool for slow C compartment
                                                      !  [mol(C)^2/m(canopy)]
    !
    ! TA
    dsl4jsb_Real2D_onChunk ::  C_sum_natural_ta
    dsl4jsb_Real2D_onChunk ::  Cpool_green_ta        ! C-pool for leaves, fine roots, vegetative organs and other green
                                                     ! .. (living) parts of vegetation [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  Cpool_reserve_ta      ! C-pool for carbohydrate reserve (sugars, starches) that allows plants
                                                     ! .. to survive bad times[mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  Cpool_woods_ta        ! C-pool for stems, thick roots and other (dead) structural material of
                                                     ! .. living plants [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  Cpool_crop_harvest_ta ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]     
    dsl4jsb_Real2D_onChunk ::  Cpool_litter_green_ag_ta ! C-pool for above ground non-woody litter especially from leaves
                                                     !     [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  Cpool_litter_green_bg_ta ! C-pool for below ground non-woody litter especially from fine roots
                                                     !     [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  Cpool_litter_wood_ag_ta  ! C-pool for litter originating from above ground woody parts of the
                                                     !    plants [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  Cpool_litter_wood_bg_ta  ! C-pool for litter originating from below ground woody parts of the
                                                     !    plants [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  Cpool_slow_ta         ! C-pool for below ground organic material (coming from the litter
                                                     !     pools) [mol(C)/m^2(canopy)]
    !SIZE CLASS 1
    dsl4jsb_Real2D_onChunk ::  YCpool_acid_ag1_ta     ! Yasso above ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_water_ag1_ta    ! Yasso above ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_ethanol_ag1_ta  ! Yasso above ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_nonsoluble_ag1_ta  ! Yasso above ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_acid_bg1_ta     ! Yasso below ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_water_bg1_ta    ! Yasso below ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_ethanol_bg1_ta  ! Yasso below ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_nonsoluble_bg1_ta  ! Yasso below ground litter-pool for non-soluble litter
                                                      ! [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_humus_1_ta      ! Yasso below ground litter-pool for slow C compartment
                                                      !  [mol(C)^2/m(canopy)]
    !SIZE CLASS 2
    dsl4jsb_Real2D_onChunk ::  YCpool_acid_ag2_ta     ! Yasso above ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_water_ag2_ta    ! Yasso above ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_ethanol_ag2_ta  ! Yasso above ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_nonsoluble_ag2_ta  ! Yasso above ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_acid_bg2_ta     ! Yasso below ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_water_bg2_ta    ! Yasso below ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_ethanol_bg2_ta  ! Yasso below ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_nonsoluble_bg2_ta  ! Yasso below ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  YCpool_humus_2_ta      ! Yasso below ground litter-pool for slow C compartment
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
    dsl4jsb_Real2D_onChunk ::  Cflux_herbivory        !
    dsl4jsb_Real2D_onChunk ::  Cflux_herbivory_LG     !
    dsl4jsb_Real2D_onChunk ::  Cflux_herbivory_2_atm  !
    dsl4jsb_Real2D_onChunk ::  CO2_flx2atm_npp_ta       ! grid cell averages of net CO2 fluxes between
    dsl4jsb_Real2D_onChunk ::  CO2_flx2atm_soilresp_ta  ! .. biosphere (due to NPP, soil respiration and
    dsl4jsb_Real2D_onChunk ::  CO2_flx2atm_herbivory_ta ! .. grazing) and atmosphere [kg(CO2)/(m^2(ground) s)] 


    ! This is on tile area:
    dsl4jsb_Real2D_onChunk ::  soil_respiration_ta      ! out
    dsl4jsb_Real2D_onChunk ::  NPP_pot_yDayMean_ta      ! out
    dsl4jsb_Real2D_onChunk ::  NPP_act_yDayMean_ta      ! out
    dsl4jsb_Real2D_onChunk ::  NPP_flux_correction_ta   ! out
    dsl4jsb_Real2D_onChunk ::  GPP_yDayMean_ta          ! out
    dsl4jsb_Real2D_onChunk ::  litter_flux_ta           ! out
    dsl4jsb_Real2D_onChunk ::  root_exudates_ta         ! out
    dsl4jsb_Real2D_onChunk ::  Cflux_herbivory_ta       ! out


    ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
    !dsl4jsb_Real2D_onChunk :: frac_litter_wood_new
    !
    dsl4jsb_Real2D_onChunk :: veg_frac_correction
    dsl4jsb_Real2D_onChunk :: LAI 
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: t_air                       ! Atmosphere temperature (lowest layer) in Kelvin!
    dsl4jsb_Real2D_onChunk :: rain
    dsl4jsb_Real2D_onChunk :: snow

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(ASSIMI_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

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
    dsl4jsb_Get_var2D_onChunk(CARBON_,  litter_flux )               ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_green )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_reserve )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_woods )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_crop_harvest)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_litter_green_ag )     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_litter_green_bg )     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_litter_wood_ag  )     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_litter_wood_bg  )     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_slow            )     ! inout
    !SIZE CLASS 1
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_acid_ag1)            ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_water_ag1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_ethanol_ag1)         ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_nonsoluble_ag1)      ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_acid_bg1)            ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_water_bg1)           ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_ethanol_bg1)         ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_nonsoluble_bg1)      ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_humus_1)             ! inout 
    !SIZE CLASS 2
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_acid_ag2)           ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_water_ag2)          ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_ethanol_ag2)        ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_nonsoluble_ag2)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_acid_bg2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_water_bg2)          ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_ethanol_bg2)        ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_nonsoluble_bg2)     ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_humus_2)            ! inout 
    !
    ! This is on tile area:
    dsl4jsb_Get_var2D_onChunk(CARBON_,  C_sum_natural_ta)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_green_ta )           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_reserve_ta )         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_woods_ta )           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_crop_harvest_ta)     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_litter_green_ag_ta ) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_litter_green_bg_ta ) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_litter_wood_ag_ta  ) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_litter_wood_bg_ta  ) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cpool_slow_ta            ) ! inout
    !SIZE CLASS 1
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_acid_ag1_ta)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_water_ag1_ta)       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_ethanol_ag1_ta)     ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_nonsoluble_ag1_ta)  ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_acid_bg1_ta)        ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_water_bg1_ta)       ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_ethanol_bg1_ta)     ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_nonsoluble_bg1_ta)  ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_humus_1_ta)         ! inout 
    !SIZE CLASS 2
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_acid_ag2_ta)       ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_water_ag2_ta)      ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_ethanol_ag2_ta)    ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_nonsoluble_ag2_ta) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_acid_bg2_ta)       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_water_bg2_ta)      ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_ethanol_bg2_ta)    ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_nonsoluble_bg2_ta) ! inout 
    dsl4jsb_Get_var2D_onChunk(CARBON_,  YCpool_humus_2_ta)        ! inout
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
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cflux_herbivory )           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cflux_herbivory_LG  )       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cflux_herbivory_2_atm )     ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  CO2_flx2atm_npp_ta )         ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  CO2_flx2atm_soilresp_ta )    ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  CO2_flx2atm_herbivory_ta )   ! out

    ! This is on tile area:
    dsl4jsb_Get_var2D_onChunk(CARBON_,  soil_respiration_ta )          ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_yDayMean_ta )          ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_act_yDayMean_ta )          ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_flux_correction_ta)        ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_yDayMean_ta )              ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  litter_flux_ta )               ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  root_exudates_ta )             ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  Cflux_herbivory_ta )           ! out

    ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
    !dsl4jsb_Get_var2D_onChunk(CARBON_,  frac_litter_wood_new )
    !
    dsl4jsb_Get_var2D_onChunk(PHENO_,   veg_frac_correction )      ! in
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
    CO2_flx2atm_npp_ta       = 0._wp
    CO2_flx2atm_soilresp_ta  = 0._wp
    CO2_flx2atm_herbivory_ta = 0._wp

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
    !    => Cpool(on PFT tile) =Cpool(canopy) * veg_frac_correction. 
    !
    ! Therefore this gets arbitary:
    ! Prepare area weighting factor to rescale from 1/[m^2(canopy)] to 1/[m^2(grid box)]
    ! areaWeightingFactor(:,:) = veg_frac_correction(:,:) * surface%cover_fract(kidx0:kidx1,:) &
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
       CALL calc_Cpools(      LAI_yDayMean,                                     & ! in
                              LAI_yyDayMean,                                    & ! in
                              NPP_pot_yDayMean,                                 & ! in
                              topSoilTemp_yDayMean,                             & ! in
                              alpha_yDayMean,                                   & ! in
                              !
                              dsl4jsb_Lctlib_param(frac_npp_2_woodPool),        & ! in
                              dsl4jsb_Lctlib_param(frac_NPP_2_reservePool),     & ! in
                              dsl4jsb_Lctlib_param(frac_NPP_2_exudates),        & ! in
                              !
                              dsl4jsb_Lctlib_param(frac_green_2_herbivory),     & ! in
                              dsl4jsb_Lctlib_param(tau_Cpool_litter_leaf),      & ! in
                              dsl4jsb_Lctlib_param(tau_Cpool_litter_wood),      & ! in
                              dsl4jsb_Lctlib_param(tau_Cpool_woods),            & ! in
                              dsl4jsb_Lctlib_param(LAI_shed_constant),          & ! in
                              dsl4jsb_Lctlib_param(frac_C_litter_green2atmos),  & ! in
                              dsl4jsb_Lctlib_param(Max_C_content_woods),        & ! in
                              dsl4jsb_Lctlib_param(specificLeafArea_C),         & ! in
                              dsl4jsb_Lctlib_param(reserveC2leafC),             & ! in
                              dsl4jsb_Lctlib_param(MaxLAI),                     & ! in
                                                                                  ! R: Rename MaxLAI in lctlib to
                                                                                  !    LAI_max in JS3 comes from
                                                                                  !    theLand%Vegetation%LAI_max,
                                                                                  !    which is for each vegetated
                                                                                  !    tile taken from MaxLAI
                                                                                  !    in the lct library.
                              dsl4jsb_Lctlib_param(PhenologyType),              & ! in  R: is it crop?
                                                                                  ! If yes then this is 5!
                              time_steps_per_day(),                             & ! in
                              dtime,                                            & ! in
                              !
                              Cpool_green,                                      & ! inout
                              Cpool_woods,                                      & ! inout
                              Cpool_reserve,                                    & ! inout
                              Cpool_crop_harvest,                               & ! inout
                              Cpool_litter_green_ag,                            & ! inout
                              Cpool_litter_green_bg,                            & ! inout
                              Cpool_litter_wood_ag,                             & ! inout
                              Cpool_litter_wood_bg,                             & ! inout
                              Cpool_slow,                                       & ! inout
                              !
                              soil_respiration,                                 & ! out
                              NPP_flux_correction,                              & ! out
                              excess_NPP,                                       & ! out
                              root_exudates,                                    & ! out
                              litter_flux,                                      & ! out
                              Cflux_herbivory,                                  & ! out
                              Cflux_herbivory_LG,                               & ! out
                              Cflux_herbivory_2_atm,                            & ! out
                              NPP_act_yDayMean,                                 & ! out
                              !
    !frac_litter_wood_new    = frac_litter_wood_new,       & ! inout, optional. Only for spitfire
                              ! variables only needed with yasso:
                              ! R: Note, in the moment the procedure is only used with Yasso.
   temp2_30d               =  pseudo_temp_yDay,                                 & ! in, optional, necessary with yasso
   precip_30d              =  pseudo_precip_yDay,                               & ! in, optional, necessary with yasso
                              !
   YCpool_acid_ag1         =  YCpool_acid_ag1,                                  & ! inout, optional, necessary with yasso
   YCpool_water_ag1        =  YCpool_water_ag1,                                 & ! inout, optional, necessary with yasso
   YCpool_ethanol_ag1      =  YCpool_ethanol_ag1,                               & ! inout, optional, necessary with yasso
   YCpool_nonsoluble_ag1   =  YCpool_nonsoluble_ag1,                            & ! inout, optional, necessary with yasso
   YCpool_acid_bg1         =  YCpool_acid_bg1,                                  & ! inout, optional, necessary with yasso
   YCpool_water_bg1        =  YCpool_water_bg1,                                 & ! inout, optional, necessary with yasso
   YCpool_ethanol_bg1      =  YCpool_ethanol_bg1,                               & ! inout, optional, necessary with yasso
   YCpool_nonsoluble_bg1   =  YCpool_nonsoluble_bg1,                            & ! inout, optional, necessary with yasso
   YCpool_humus_1          =  YCpool_humus_1,                                   & ! inout, optional, necessary with yasso
   YCpool_acid_ag2         =  YCpool_acid_ag2,                                  & ! inout, optional, necessary with yasso
   YCpool_water_ag2        =  YCpool_water_ag2,                                 & ! inout, optional, necessary with yasso
   YCpool_ethanol_ag2      =  YCpool_ethanol_ag2,                               & ! inout, optional, necessary with yasso
   YCpool_nonsoluble_ag2   =  YCpool_nonsoluble_ag2,                            & ! inout, optional, necessary with yasso
   YCpool_acid_bg2         =  YCpool_acid_bg2,                                  & ! inout, optional, necessary with yasso
   YCpool_water_bg2        =  YCpool_water_bg2,                                 & ! inout, optional, necessary with yasso
   YCpool_ethanol_bg2      =  YCpool_ethanol_bg2,                               & ! inout, optional, necessary with yasso
   YCpool_nonsoluble_bg2   =  YCpool_nonsoluble_bg2,                            & ! inout, optional, necessary with yasso
   YCpool_humus_2          =  YCpool_humus_2,                                   & ! inout, optional, necessary with yasso
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
       !    specific_leaf_area_C  (=> Cpool_green_max und Cpool_reserve_optimal), Max_C_content_woods. usw.
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
       Cpool_green_ta           = Cpool_green           * veg_frac_correction * fract_fpc_max
       Cpool_woods_ta           = Cpool_woods           * veg_frac_correction * fract_fpc_max
       Cpool_reserve_ta         = Cpool_reserve         * veg_frac_correction * fract_fpc_max
       Cpool_crop_harvest_ta    = Cpool_crop_harvest    * veg_frac_correction * fract_fpc_max

       YCpool_acid_ag1_ta       = YCpool_acid_ag1       * veg_frac_correction * fract_fpc_max
       YCpool_water_ag1_ta      = YCpool_water_ag1      * veg_frac_correction * fract_fpc_max
       YCpool_ethanol_ag1_ta    = YCpool_ethanol_ag1    * veg_frac_correction * fract_fpc_max
       YCpool_nonsoluble_ag1_ta = YCpool_nonsoluble_ag1 * veg_frac_correction * fract_fpc_max
       YCpool_acid_bg1_ta       = YCpool_acid_bg1       * veg_frac_correction * fract_fpc_max
       YCpool_water_bg1_ta      = YCpool_water_bg1      * veg_frac_correction * fract_fpc_max
       YCpool_ethanol_bg1_ta    = YCpool_ethanol_bg1    * veg_frac_correction * fract_fpc_max
       YCpool_nonsoluble_bg1_ta = YCpool_nonsoluble_bg1 * veg_frac_correction * fract_fpc_max
       YCpool_humus_1_ta        = YCpool_humus_1        * veg_frac_correction * fract_fpc_max
       YCpool_acid_ag2_ta       = YCpool_acid_ag2       * veg_frac_correction * fract_fpc_max
       YCpool_water_ag2_ta      = YCpool_water_ag2      * veg_frac_correction * fract_fpc_max
       YCpool_ethanol_ag2_ta    = YCpool_ethanol_ag2    * veg_frac_correction * fract_fpc_max
       YCpool_nonsoluble_ag2_ta = YCpool_nonsoluble_ag2 * veg_frac_correction * fract_fpc_max
       YCpool_acid_bg2_ta       = YCpool_acid_bg2       * veg_frac_correction * fract_fpc_max
       YCpool_water_bg2_ta      = YCpool_water_bg2      * veg_frac_correction * fract_fpc_max
       YCpool_ethanol_bg2_ta    = YCpool_ethanol_bg2    * veg_frac_correction * fract_fpc_max
       YCpool_nonsoluble_bg2_ta = YCpool_nonsoluble_bg2 * veg_frac_correction * fract_fpc_max
       YCpool_humus_2_ta        = YCpool_humus_2        * veg_frac_correction * fract_fpc_max


       C_sum_natural_ta =   Cpool_green_ta           + Cpool_woods_ta            &
                          + Cpool_reserve_ta         + Cpool_crop_harvest_ta     &
                          + YCpool_acid_ag1_ta       + YCpool_acid_ag2_ta        &
                          + YCpool_water_ag1_ta      + YCpool_water_ag2_ta       &
                          + YCpool_ethanol_ag1_ta    + YCpool_ethanol_ag2_ta     &
                          + YCpool_nonsoluble_ag1_ta + YCpool_nonsoluble_ag2_ta  &
                          + YCpool_acid_bg1_ta       + YCpool_acid_bg2_ta        &
                          + YCpool_water_bg1_ta      + YCpool_water_bg2_ta       &
                          + YCpool_ethanol_bg1_ta    + YCpool_ethanol_bg2_ta     &
                          + YCpool_nonsoluble_bg1_ta + YCpool_nonsoluble_bg2_ta  &
                          + YCpool_humus_1_ta        + YCpool_humus_2_ta

       ! C fluxes
       soil_respiration_ta         = soil_respiration      * veg_frac_correction * fract_fpc_max
       !soil_respiration_pot_ta     = soil_respiration_pot  * veg_frac_correction * fract_fpc_max
       NPP_pot_yDayMean_ta         = NPP_pot_yDayMean      * veg_frac_correction * fract_fpc_max
       NPP_act_yDayMean_ta         = NPP_act_yDayMean      * veg_frac_correction * fract_fpc_max
       NPP_flux_correction_ta      = NPP_flux_correction   * veg_frac_correction * fract_fpc_max
       GPP_yDayMean_ta             = GPP_yDayMean          * veg_frac_correction * fract_fpc_max
       litter_flux_ta              = litter_flux           * veg_frac_correction * fract_fpc_max
       root_exudates_ta            = root_exudates         * veg_frac_correction * fract_fpc_max
       Cflux_herbivory_ta          = Cflux_herbivory       * veg_frac_correction * fract_fpc_max
       ! Cflux_herbivory_LG_ta      = Cflux_herbivory_LG    * veg_frac_correction * fract_fpc_max
       ! Cflux_herbivory_2_atm_ta   = Cflux_herbivory_2_atm * veg_frac_correction * fract_fpc_max

   END IF

   ! Compute net CO2 fluxes exchanged with atmosphere at each time step
   !-------------------------------------------------------------------
   ! Note: carbon loss of biosphere means a positive CO2 flux to atmosphere (i.e. NEP and net CO2-flux have opposite signs)

   CO2_flx2atm_npp_ta = molarMassCO2_kg *                     & ! Conversion factor from mol to kg CO2
              (- veg_frac_correction * fract_fpc_max          & ! Minus: atmosphere gain is positive
                   * (NPP_pot_rate_ca                         & ! current (not actual) NPP rate
                       - (NPP_pot_yDayMean-NPP_act_yDayMean)) ) ! corrected with yesterdays actual NPP defizit


   CO2_flx2atm_soilresp_ta = molarMassCO2_kg *                & ! Conversion factor from mol to kg CO2
              (- veg_frac_correction *fract_fpc_max           & ! Minus: atmosphere gain is positive
                   * soil_respiration      )                    ! .. soil respiration

   CO2_flx2atm_herbivory_ta = molarMassCO2_kg *               & ! Conversion factor from mol to kg CO2
              (- veg_frac_correction * fract_fpc_max          & ! Minus: atmosphere gain is positive
                   * Cflux_herbivory_2_atm )                    ! .. herbivory


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

    dsl4jsb_Aggregate_onChunk(CARBON_, litter_flux_ta,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, Cpool_green_ta,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, Cpool_woods_ta,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, Cpool_reserve_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, Cpool_crop_harvest_ta,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, Cpool_litter_green_ag_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, Cpool_litter_green_bg_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, Cpool_litter_wood_ag_ta,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, Cpool_litter_wood_bg_ta,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, Cpool_slow_ta,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_acid_ag1_ta,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_water_ag1_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_ethanol_ag1_ta,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_nonsoluble_ag1_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_acid_bg1_ta,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_water_bg1_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_ethanol_bg1_ta,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_nonsoluble_bg1_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_humus_1_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_acid_ag2_ta,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_water_ag2_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_ethanol_ag2_ta,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_nonsoluble_ag2_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_acid_bg2_ta,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_water_bg2_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_ethanol_bg2_ta,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_nonsoluble_bg2_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, YCpool_humus_2_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, C_sum_natural_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_pot_sum_ta,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, GPP_sum_ta,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, soil_respiration_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_flux_correction_ta,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, excess_NPP_ta,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yyDayMean_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yDayMean_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_pot_yDayMean_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_act_yDayMean_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, GPP_yDayMean_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, root_exudates_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, Cflux_herbivory_ta,       weighted_by_fract)
   ! dsl4jsb_Aggregate_onChunk(CARBON_, Cflux_herbivory_LG_ta, weighted_by_fract)
   ! dsl4jsb_Aggregate_onChunk(CARBON_, Cflux_herbivory_2_atm_ta, weighted_by_fract)
   ! dsl4jsb_Aggregate_onChunk(CARBON_, frac_litter_wood_new_ta,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, CO2_flx2atm_npp_ta,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, CO2_flx2atm_soilresp_ta,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, CO2_flx2atm_herbivory_ta, weighted_by_fract)
   ! dsl4jsb_Aggregate_onChunk(CARBON_, CO2_flx2atm_npp_test_ta,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_C_NPP_pot_allocation



  ! ================================================================================================================================
  !>
  !! @brief Implementation to relocate carbon caused by disturbances.
  !!        Relocation between the plant and soil pools and the atmosphere
  !! 
  !!
  SUBROUTINE relocate_disturbed_carbon(tile, options,                & ! in
                                       dist_type,                    & ! in (DIST_FIRE or DIST_WINDBREAK)
                                       no_children,                  & ! in
                                       distrubed_frac_st,            & ! in (from burned_frac_st or damaged_frac_st)
                                       fract_st,                     & ! in
                                       veg_frac_correction_st,       & ! in
                                       Cpool_green_st,               & ! inout
                                       Cpool_woods_st,               & ! inout
                                       Cpool_reserve_st,             & ! inout
                                       carbon_2_GreenLitterPools_st, & ! inout
                                       carbon_2_WoodLitterPools_st,  & ! inout
                                       YCpool_acid_ag1_st,           & ! inout
                                       YCpool_water_ag1_st,          & ! inout
                                       YCpool_ethanol_ag1_st,        & ! inout
                                       YCpool_nonsoluble_ag1_st,     & ! inout
                                       YCpool_acid_ag2_st,           & ! inout
                                       YCpool_water_ag2_st,          & ! inout
                                       YCpool_ethanol_ag2_st,        & ! inout
                                       YCpool_nonsoluble_ag2_st,     & ! inout
                                       !
                                       YCpool_acid_bg1_st,           & ! inout
                                       YCpool_water_bg1_st,          & ! inout
                                       YCpool_ethanol_bg1_st,        & ! inout
                                       YCpool_nonsoluble_bg1_st,     & ! inout
                                       YCpool_acid_bg2_st,           & ! inout
                                       YCpool_water_bg2_st,          & ! inout
                                       YCpool_ethanol_bg2_st,        & ! inout
                                       YCpool_nonsoluble_bg2_st,     & ! inout
                                       !
                                       YCpool_humus_1_st,            & ! inout
                                       YCpool_humus_2_st,            & ! inout
                                       carbon_2_atmos_st_ta          & ! out, used only with fire
                                       )

    USE mo_exception,             ONLY: message, finish
    USE mo_hsm_class,             ONLY: t_State
    USE mo_carbon_process,        ONLY: relocate_carbon_fire, relocate_carbon_damage
    USE mo_disturbance_constants, ONLY: DIST_FIRE, DIST_WINDBREAK

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile           ! R: see X1 in mo_disturbance_interface.f90
    TYPE(t_jsb_task_options),   INTENT(in)         :: options

    INTEGER,                    INTENT(in)         :: dist_type,   & ! Sub-disturbance type (one of the DIST_ constants)
                                                      no_children    ! Number of subtiles

    ! The following variables all have the dimensions (ics:ice, no_children). ics:ice values do not change.
    ! They all collect the corresponding variables of the child pfts
    REAL(wp),                   INTENT(in)         :: distrubed_frac_st(:,:),            &
                                                      fract_st(:,:),                     &
                                                      veg_frac_correction_st(:,:)

    REAL(wp),                   INTENT(inout)      :: Cpool_green_st(:,:),               &
                                                      Cpool_woods_st(:,:),               &
                                                      Cpool_reserve_st(:,:),             &
                                                      carbon_2_GreenLitterPools_st(:,:), &
                                                      carbon_2_WoodLitterPools_st(:,:),  &
                                                      YCpool_acid_ag1_st(:,:),           &
                                                      YCpool_water_ag1_st(:,:),          &
                                                      YCpool_ethanol_ag1_st(:,:),        &
                                                      YCpool_nonsoluble_ag1_st(:,:),     &
                                                      YCpool_acid_ag2_st(:,:),           &
                                                      YCpool_water_ag2_st(:,:),          &
                                                      YCpool_ethanol_ag2_st(:,:),        &
                                                      YCpool_nonsoluble_ag2_st(:,:),     &
                                                      !
                                                      YCpool_acid_bg1_st(:,:),           &
                                                      YCpool_water_bg1_st(:,:),          &
                                                      YCpool_ethanol_bg1_st(:,:),        &
                                                      YCpool_nonsoluble_bg1_st(:,:),     &
                                                      YCpool_acid_bg2_st(:,:),           &
                                                      YCpool_water_bg2_st(:,:),          &
                                                      YCpool_ethanol_bg2_st(:,:),        &
                                                      YCpool_nonsoluble_bg2_st(:,:),     &
                                                      !
                                                      YCpool_humus_1_st(:,:),            &
                                                      YCpool_humus_2_st(:,:)
                                                 !cover_frac_pot_st(:,:)

    REAL(wp),                   INTENT(out)   :: carbon_2_atmos_st_ta(:,:)  ! Carbon immedieatly released by fire per tile
                                                                            ! [mol(C)m-2(subtile)]

    ! Local variables
    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_config(DISTURB_)

    TYPE(t_jsb_model), POINTER  :: model
    INTEGER                     :: iblk, ics, ice, nc

    INTEGER                 :: pft
    REAL(wp), allocatable   :: LeafLit_coef(:,:,:) ! Factor to seperate non woody litterfall into yasso pools
                                                   ! (ics:ice,no_children,chemical_classes)
    REAL(wp), allocatable   :: WoodLit_coef(:,:,:) ! Factor to seperate woody litterfall into yasso pools
                                                   ! (ics:ice,no_children,chemical_classes)

    !CLASS(t_jsb_tile_abstract), POINTER :: act_subtile_ptr => null() ! actual pointer to the child tiles of this tile
    CLASS(t_State), POINTER :: act_subtile_ptr => null() ! actual pointer to the child tiles of this tile

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! Get model
    model         => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(CARBON_)
    dsl4jsb_Get_config(DISTURB_)

    ! ---------------------------
    ! Go

    ! Initializations
    carbon_2_atmos_st_ta(:,:)         = 0._wp ! not necessary, just to be sure
    carbon_2_GreenLitterPools_st(:,:) = 0._wp ! not necessary, just to be sure
    carbon_2_WoodLitterPools_st(:,:)  = 0._wp ! not necessary, just to be sure

    ! R: For JSBACH4 I removed the following proof that the values in lctlib are consistent with YASSO.
    !    This was in jsbach_inter_1d and should be moved somewhere else in JS3!
    !DO itile=1,ntiles
    !   LeafLit_coef(1:nidx,itile,:) = &
    !        theLand%LctLibrary%LeafLit_coef(theLand%Surface%cover_type(kidx0:kidx1,itile),:)
    !   WoodLit_coef(1:nidx,itile,:) = &
    !        theLand%LctLibrary%WoodLit_coef(theLand%Surface%cover_type(kidx0:kidx1,itile),:)
    !   DO i = 1,nidx 
    !      IF ( ABS(SUM(LeafLit_coef(i,itile,:)) - 1.0_dp) > EPSILON(1.0_dp) ) THEN
    !         WRITE(*,*) "Leaf litter coefficients have invalid values; sum must be 1, but is:" 
    !         WRITE(*,*) SUM(LeafLit_coef(i,itile,:)), "for tile=", itile
    !         CALL finish('jsbach_interfacei_1d','invalid sum of leaf litter coefficients')
    !      END IF

    !      IF ( ABS(SUM(WoodLit_coef(i,itile,:)) - 1.0_dp) > EPSILON(1.0_dp) ) THEN
    !         WRITE(*,*) "Wood litter coefficients have invalid values; sum must be 1, but is:" 
    !         WRITE(*,*) SUM(WoodLit_coef(i,itile,:)), "for tile=", itile
    !         CALL finish('jsbach_interfacei_1d','invalid sum of wood litter coefficients')
    !      END IF
    !      IF (theOptions%Withnitrogen) THEN
    !         IF ( ABS(LeafLit_coef(i,itile,5)) > EPSILON(1.0_dp) ) THEN
    !            WRITE(*,*) "Leaf litter coefficient has invalid value for use with nitrogen"
    !            WRITE(*,*) "Leaf litter coefficient must be 0 for humus, but is:" 
    !            WRITE(*,*) LeafLit_coef(i,itile,5), "for tile=", itile
    !            CALL finish('jsbach_interfacei_1d','non zero leaf litter coefficient for humus')
    !         ENDIF
    !         IF ( ABS(WoodLit_coef(i,itile,5)) > EPSILON(1.0_dp) ) THEN
    !            WRITE(*,*) "Wood litter coefficient has invalid value for use with nitrogen"
    !            WRITE(*,*) "Wood litter coefficient must be 0 for humus, but is:" 
    !            WRITE(*,*) WoodLit_coef(i,itile,5), "for tile=", itile
    !            CALL finish('jsbach_interfacei_1d','non zero wood litter coefficient for humus') 
    !         ENDIF
    !      ENDIF
    !   END DO
    !END DO

    ! R: We have to prepare the LeafLit_coef(ics:ice,no_children,chemical_classes) from the lct library
    !    dsl4jsb_Lctlib_param(LeafLit_coef(1:5))
    DO pft=1,no_children
       IF ( pft == 1 ) THEN
          act_subtile_ptr => tile%first_child
       ELSE
          !act_subtile_ptr =>  act_subtile_ptr%next_sibling ! R: das geht bestimmt nicht, also wirds umständlich weil wir dann einen
                                                            !    temporäten Pointer brauchen. Aber das macht die Funktion für uns:
          act_subtile_ptr => tile%Get_next_sibling_tile()
       END IF

       SELECT TYPE (act_subtile_ptr)
            CLASS IS (t_jsb_tile_abstract)
            ! R: I'm not sure with DIM=2 ...
            LeafLit_coef(1:nc,pft,1:5) = SPREAD( model%lctlib(act_subtile_ptr%lcts(1)%lib_id)%LeafLit_coef(1:5), DIM=2,NCOPIES=nc)
            WoodLit_coef(1:nc,pft,1:5) = SPREAD( model%lctlib(act_subtile_ptr%lcts(1)%lib_id)%WoodLit_coef(1:5), DIM=2,NCOPIES=nc)
       END SELECT
    ENDDO

    ! call the relocation routines
    IF ( dist_type == DIST_FIRE ) THEN 
      SELECT CASE (dsl4jsb_Config(CARBON_)%fuel_algorithm)
        CASE (0) !! No fuel consumption
        CASE (1) !! jsbach fuel consumption
            CALL relocate_carbon_fire(           &
              no_children,                       & ! in
              distrubed_frac_st(:,:),            & ! in
              fract_st(:,:),                     & ! in
              veg_frac_correction_st(:,:),       & ! in
              dsl4jsb_Config(CARBON_)%fire_frac_wood_2_atmos , & ! in
              Cpool_green_st(:,:),               & ! inout
              Cpool_woods_st(:,:),               & ! inout
              Cpool_reserve_st(:,:),             & ! inout
              carbon_2_atmos_st_ta(:,:),         & ! out
              carbon_2_GreenLitterPools_st(:,:), & ! inout, because relocate_carbon_fire is 
                                                   ! called AFTER relocate_carbon_damage!
              carbon_2_WoodLitterPools_st(:,:),  & ! inout, because relocate_carbon_fire is 
                                                   ! called AFTER relocate_carbon_damage!
                !Cpool_litter_green_ag,          & ! inout, 
                !Cpool_litter_green_bg,          & ! inout, 
                !Cpool_litter_wood_ag,           & ! inout, 
                !Cpool_litter_wood_bg,           & ! inout, 
               ! Yasso
              YCpool_acid_ag1_st      (:,:),     & ! inout, necessary with yasso
              YCpool_acid_ag2_st      (:,:),     & ! inout, necessary with yasso
              YCpool_water_ag1_st     (:,:),     & ! inout, necessary with yasso
              YCpool_water_ag2_st     (:,:),     & ! inout, necessary with yasso
              YCpool_ethanol_ag1_st   (:,:),     & ! inout, necessary with yasso
              YCpool_ethanol_ag2_st   (:,:),     & ! inout, necessary with yasso
              YCpool_nonsoluble_ag1_st(:,:),     & ! inout, necessary with yasso
              YCpool_nonsoluble_ag2_st(:,:),     & ! inout, necessary with yasso
              YCpool_acid_bg1_st      (:,:),     & ! inout, necessary with yasso
              YCpool_acid_bg2_st      (:,:),     & ! inout, necessary with yasso
              YCpool_water_bg1_st     (:,:),     & ! inout, necessary with yasso
              YCpool_water_bg2_st     (:,:),     & ! inout, necessary with yasso
              YCpool_ethanol_bg1_st   (:,:),     & ! inout, necessary with yasso
              YCpool_ethanol_bg2_st   (:,:),     & ! inout, necessary with yasso
              YCpool_nonsoluble_bg1_st(:,:),     & ! inout, necessary with yasso
              YCpool_nonsoluble_bg2_st(:,:),     & ! inout, necessary with yasso
              YCpool_humus_1_st       (:,:),     & ! inout, necessary with yasso
              YCpool_humus_2_st       (:,:),     & ! inout, necessary with yasso
              LeafLit_coef            (:,:,:),   & ! in, optional, necessary with yasso
                                                   ! Coefficient to distribute leaf litter into 5 classes
                                                   !               of chemical composition
              WoodLit_coef            (:,:,:)    & ! in, optional, necessary with yasso
                                     )

        CASE (3) !! relocation of carbon according to the spitfire algorithm

              CALL finish('disturbed_frac','Thornicke spit fire not implemented yet.')

        CASE (4) !! Read fuel consumption
        CASE DEFAULT
          CALL finish('relocate_disturbed_carbon','Unknown fire algorithm')
      END SELECT

      ! R: lemissions Schalter ist für andere Emissionen als C z.B. NOx. Daher den Block hier raus genommen.
      ! Update emissions
      ! IF (dist_opts%lemissions) THEN ... 


    ! Disturbance by windbreak
    ELSEIF ( dist_type == DIST_WINDBREAK ) THEN
      SELECT CASE (dsl4jsb_Config(DISTURB_)%windbreak_algorithm)
        CASE (0) !! No windbreak algorithm
        CASE (1) !! jsbach windbreak algorithm

            ! with yasso, no nitrogen
            CALL relocate_carbon_damage(         &
              no_children,                       & ! in
              distrubed_frac_st(:,:),            & ! in
              fract_st(:,:),                     & ! in
              veg_frac_correction_st(:,:),       & ! in
              Cpool_green_st(:,:),               & ! inout
              Cpool_woods_st(:,:),               & ! inout
              Cpool_reserve_st(:,:),             & ! inout
              carbon_2_GreenLitterPools_st(:,:), & ! out, because relocate_carbon_damage is
                                                   ! called before relocate_carbon_fire!
              carbon_2_WoodLitterPools_st(:,:),  & ! out, because relocate_carbon_damage is 
                                                   ! called before relocate_carbon_fire!
              !Cpool_litter_green_ag,       & ! inout, optional
              !Cpool_litter_green_bg,       & ! inout, optional
              !Cpool_litter_wood_ag,        & ! inout, optional
              !Cpool_litter_wood_bg,        & ! inout, optional
              ! Yasso
              YCpool_acid_ag1_st      (:,:),     & ! inout, necessary with yasso
              YCpool_acid_ag2_st      (:,:),     & ! inout, necessary with yasso
              YCpool_water_ag1_st     (:,:),     & ! inout, necessary with yasso
              YCpool_water_ag2_st     (:,:),     & ! inout, necessary with yasso
              YCpool_ethanol_ag1_st   (:,:),     & ! inout, necessary with yasso
              YCpool_ethanol_ag2_st   (:,:),     & ! inout, necessary with yasso
              YCpool_nonsoluble_ag1_st(:,:),     & ! inout, necessary with yasso
              YCpool_nonsoluble_ag2_st(:,:),     & ! inout, necessary with yasso
              YCpool_acid_bg1_st      (:,:),     & ! inout, necessary with yasso
              YCpool_acid_bg2_st      (:,:),     & ! inout, necessary with yasso
              YCpool_water_bg1_st     (:,:),     & ! inout, necessary with yasso
              YCpool_water_bg2_st     (:,:),     & ! inout, necessary with yasso
              YCpool_ethanol_bg1_st   (:,:),     & ! inout, necessary with yasso
              YCpool_ethanol_bg2_st   (:,:),     & ! inout, necessary with yasso
              YCpool_nonsoluble_bg1_st(:,:),     & ! inout, necessary with yasso
              YCpool_nonsoluble_bg2_st(:,:),     & ! inout, necessary with yasso
              YCpool_humus_1_st       (:,:),     & ! inout, necessary with yasso
              YCpool_humus_2_st       (:,:),     & ! inout, necessary with yasso
              LeafLit_coef            (:,:,:),   & ! in, necessary with yasso
                                                   ! Coefficient to distribute leaf litter into 5 classes 
                                                   !               of chemical composition
              WoodLit_coef            (:,:,:)    & ! in, necessary with yasso
                                      )

          ! Insert other windbreak algorithms here
        CASE DEFAULT
          CALL finish('relocate_disturbed_carbon','Unknown windbreak algorithm')
      END SELECT
    ENDIF

  END SUBROUTINE relocate_disturbed_carbon


#endif
END MODULE mo_carbon_interface
