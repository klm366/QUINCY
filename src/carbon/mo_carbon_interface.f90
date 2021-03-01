!>
!! @brief Contains the interfaces to the carbon process
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-04-05)
!! Julia Nabel (2018): l_forestRegrowth flag -> allow for changing MaxLai; forest age cohorts
!!      ++ subroutines for cpool sums and c conservation test
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_carbon_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  USE mo_carbon_constants,   ONLY: molarMassCO2_kg, sec_per_day

  ! Use of processes in this module (Get_carbon_memory and Get_carbon_config)
  dsl4jsb_Use_processes CARBON_, ASSIMI_, PHENO_, DISTURB_, A2L_, L2A_,HYDRO_,RAD_

  ! Use of process configurations
  dsl4jsb_Use_config(CARBON_)
  dsl4jsb_Use_config(PHENO_)

  ! Use of process memories (t_carbon_memory)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(L2A_)

  ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  ! Use the ELEMENT FUNCTION from HYDRO for JeDi : get_water_stress_factor
  ! 
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(RAD_)
  ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_carbon_tasks
  PUBLIC ::  recalc_per_tile_vars, calculate_current_c_ag_1_and_bg_sums, rescale_carbon_upon_reference_area_change, &
    &        calculate_current_c_ta_state_sum, check_carbon_conservation, yday_carbon_conservation_test
  PUBLIC ::  jedi_biomass_to_fract

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
  
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
 !PHS 24.07.2020 Compute water stress in CARBON at PFT level. 
 !> Type definition for water_stress task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_water_stress
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_water_stress
    PROCEDURE, NOPASS :: Aggregate => aggregate_water_stress
  END TYPE tsk_water_stress

  !> Constructor interface for water_stress task
  INTERFACE tsk_water_stress
    PROCEDURE Create_task_water_stress
  END INTERFACE tsk_water_stress
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

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


! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
 !PHS 24.07.2020 Compute water stress in CARBON at PFT level. 
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for water_stress task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "water_stress"
  !!
  FUNCTION Create_task_water_stress(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_water_stress::return_ptr)
    CALL return_ptr%Construct(name='water_stress', process_id=CARBON_, owner_model_id=model_id)

  END FUNCTION Create_task_water_stress
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 


  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for carbon process
  !!
  !! @param[in,out] this      Instance of carbon process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_carbon_tasks(this, model_id)

    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    ! in/out
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    ! local
    TYPE(t_jsb_model), POINTER :: model

    ! get var / objects 
    model   => Get_model(model_id)


    IF (.NOT. model%config%use_quincy) THEN
      CALL this%Register_task(tsk_C_NPP_pot_allocation(model_id))
      ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
      !PHS 24.07.2020 Compute water stress in CARBON at PFT level. 
       IF (TRIM(model%config%usecase) .eq. 'jsbach_jedi' ) THEN        
          CALL this%Register_task( tsk_water_stress           (model_id))
      END IF 
      ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

    ENDIF

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
    USE mo_jsb_grid,          ONLY: Get_grid,Get_vgrid
    USE mo_util,    ONLY: soil_depth_to_layers_2d
    USE mo_carbon_process,    ONLY: calc_Cpools
    USE mo_carbon_jedi_process, ONLY: calc_Ctradeoff,jedi_water_stress_factor,jedi_soil_depth_to_layers_1d!,drop_expensive_tissue
    USE mo_jsb_physical_constants, ONLY : tmelt  
    USE mo_jsb_control,        ONLY:  jsbach_is_restarted
    USE mo_jsb_time,          ONLY: t_datetime,is_newday, is_newyear, timesteps_per_day, &
                                 &  is_time_experiment_start, get_year_at_model_start,     &
                                 &  get_time_experiment_start, get_year, deallocateDatetime
    USE mo_carbon_constants,  ONLY: max_pipe_transpir
    USE mo_jsb_grid_class,     ONLY:  t_jsb_vgrid !PHS: root_depth
    USE mo_hydro_util,             ONLY: get_water_in_root_zone !PHS: root_depth
    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Declare local variables
    TYPE(t_jsb_model), POINTER   :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_w   !PHS: soil depth
    LOGICAL                      :: lstart
    INTEGER                      :: iblk, ics, ice, nc, i
    REAL(wp)                     :: dtime
    REAL(wp), allocatable        :: precip_total(:)
    REAL(wp), DIMENSION(options%nc) :: MaxLai
    REAL(wp), DIMENSION(options%nc) :: old_c_state_sum, current_fluxes
    CHARACTER(len=*), PARAMETER  :: routine = modname//':update_C_NPP_pot_allocation'
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Declare variavles for JeDi carbon allocation 
    ! 
    TYPE(t_datetime), POINTER    :: inidate
    INTEGER                      :: iniyear
    LOGICAL                      :: lrestart
    LOGICAL                      :: lspinup_jedi
  
    !REAL(wp), parameter  :: rootshoot_opt     = 3.58462104_wp ! median value of root-shoot ratio 50th percentile
    REAL(wp), parameter  :: rootshoot_opt     = 8.473583644190947_wp ! median value of root-shoot ratio 75th percentile
    REAL(wp), parameter  :: c_CN_woody       = 0.00303_wp * (12.01_wp / 14.0067_wp)  ! molC/molN Ratio for Woody Tissue
    REAL(wp), parameter  :: c_sapwood        = 0.05_wp
    REAL(wp), parameter  :: c_Res_N          = 0.218_wp * ( 14.0067_wp / 12.01_wp )  !! specific respiration in molC/molN/d
    REAL(wp), PARAMETER  :: Q10              = 1.8_wp          ! Empirical parameter for temperature dependence of
                                                                !   heterotrophic respiration rate. This value is 
                                                                !   taken from JB4.
    REAL(wp), parameter  :: cTiny            = 1.E-8_wp
    REAL(wp)             ::  dTdt

  

    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
    

    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_config(CARBON_)
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! JeDi needs soil moisture from HYDRO process
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(RAD_)
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    
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

    dsl4jsb_Real2D_onChunk ::  cconservation_calcCpools ! C conservation test: Deviation from conservation

    dsl4jsb_Real2D_onChunk ::  cflux_c_greenwood_2_litter           ! Carbon flux from the veget. to the litter pools
                                                     ! [mol(C)/m^2(canopy) s]
    dsl4jsb_Real2D_onChunk ::  c_green           ! C-pool for leaves, fine roots, vegetative organs and
                                                     ! other green (living) parts of vegetation [mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_reserve         ! C-pool for carbohydrate reserve (sugars, starches) that
                                                     ! allows plants to survive bad times[mol(C)/m^2(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_woods           ! C-pool for stems, thick roots and other (dead) structural
                                                     !  material of living plants [mol(C)/m^2(canopy)]
  
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    !
    dsl4jsb_Real2D_onChunk ::   c_storage               !     Storage carbon pool 
    dsl4jsb_Real2D_onChunk ::   c_leaves                !      Leaves carbon pool             
    dsl4jsb_Real2D_onChunk ::   c_fnroot                !   Fine root carbon pool
    dsl4jsb_Real2D_onChunk ::   c_csroot                ! Coarse root carbon pool
    dsl4jsb_Real2D_onChunk ::   c_seed                  !        Seed carbon pool
    dsl4jsb_Real2D_onChunk ::   c_stem                  !        Stem carbon pool
    dsl4jsb_Real2D_onChunk ::   ff_grow_W    
    dsl4jsb_Real2D_onChunk ::   ff_grow_T
    dsl4jsb_Real2D_onChunk ::   ff_grow_G
    dsl4jsb_Real2D_onChunk ::   Flg_grow
    dsl4jsb_Real2D_onChunk ::   Flg_grow_G
    dsl4jsb_Real2D_onChunk ::   ff_npp
    dsl4jsb_Real2D_onChunk ::   flx_germ
    dsl4jsb_Real2D_onChunk ::   adj_fraction
    dsl4jsb_Real2D_onChunk ::   jedi_light_constraint
    dsl4jsb_Real2D_onChunk ::   Flg_die
    dsl4jsb_Real2D_onChunk ::   Flg_senesc
    dsl4jsb_Real2D_onChunk ::   MaintRes
    dsl4jsb_Real2D_onChunk ::   MaintRes_sum
    dsl4jsb_Real2D_onChunk ::   MaintRes_yDayMean
    dsl4jsb_Real2D_onChunk ::   temp_sum
    dsl4jsb_Real2D_onChunk ::   temp_yDayMean
    dsl4jsb_Real2D_onChunk ::   GrowthResp    
    dsl4jsb_Real2D_onChunk ::   AutoRes
    dsl4jsb_Real2D_onChunk ::  CUE
    dsl4jsb_Real2D_onChunk ::  root_maxmoist
    dsl4jsb_Real2D_onChunk ::  seed_maxmoist
    dsl4jsb_Real2D_onChunk ::  jedi_pot_supply_transp
    dsl4jsb_Real2D_onChunk ::  jedi_root_depth
    dsl4jsb_Real2D_onChunk ::  soil_depth
    dsl4jsb_Real3D_onChunk ::  jedi_root_depth_sl
    dsl4jsb_Real3D_onChunk ::  vol_field_cap_sl
    dsl4jsb_Real3D_onChunk ::  soil_depth_sl
    dsl4jsb_Real3D_onChunk ::  w_soil_sl
    dsl4jsb_Real3D_onChunk ::  seed_depth_sl
    dsl4jsb_Real2D_onChunk ::   C_balan
    dsl4jsb_Real2D_onChunk ::   C_delta
    dsl4jsb_Real2D_onChunk ::   c_total
    dsl4jsb_Real2D_onChunk ::   one_counter
    ! Taking relative soil moisture from the HYDRO process            
    dsl4jsb_Real2D_onChunk ::   alpha_jedi
    dsl4jsb_Real2D_onChunk ::   alpha_jedi_seed
    dsl4jsb_Real2D_onChunk ::   max_moist
    dsl4jsb_Real2D_onChunk ::   w_soil_column_seed
    dsl4jsb_Real2D_onChunk ::   w_soil_column_jedi
    dsl4jsb_Real2D_onChunk ::   ff_T_30d
    dsl4jsb_Real2D_onChunk ::   ff_dTdt        
    dsl4jsb_Real2D_onChunk ::   pseudo_c_total              
    dsl4jsb_Real2D_onChunk ::   ff_grow_T_yDayMean
    dsl4jsb_Real2D_onChunk ::   faPAR_prev_noon
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
    
                                                        
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

    dsl4jsb_Real2D_onChunk ::  c_decomp_humus_1_sum  ! Annual sum of humus decomposed in YASSO (leaf)
                                                      !  [mol(C)^2/m(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_decomp_humus_2_sum  ! Annual sum of humus decomposed in YASSO (wood)
                                                      !  [mol(C)^2/m(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_into_humus_1_sum    ! Annual sum of cflux into humus in YASSO (leaf)
                                                      !  [mol(C)^2/m(canopy)]
    dsl4jsb_Real2D_onChunk ::  c_into_humus_2_sum    ! Annual sum of cflux into humus in YASSO (wood)
                                                      !  [mol(C)^2/m(canopy)]

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
    dsl4jsb_Real2D_onChunk ::  NPP_act_yDayMean       ! mean value of actual NPP-Rate yesterday, i.e. after N-limitation
                                                      ! Actual NPP after N-limitation and excess carbon drop.
                                                      ! [mol(CO2)/(m^2(canopy) s)]
    dsl4jsb_Real2D_onChunk ::  GPP_yDayMean

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

    dsl4jsb_Real2D_onChunk ::  c_sum_litter_ag_ta
    dsl4jsb_Real2D_onChunk ::  c_sum_litter_bg_ta
    dsl4jsb_Real2D_onChunk ::  c_sum_humus_ta
    dsl4jsb_Real2D_onChunk ::  c_sum_natural_ta

    dsl4jsb_Real2D_onChunk ::  cflux_c_green_2_herb     !
    dsl4jsb_Real2D_onChunk ::  cflux_herb_2_littergreen !
    dsl4jsb_Real2D_onChunk ::  cflux_herb_2_atm  !
    dsl4jsb_Real2D_onChunk ::  cflux_npp_2_atm_yday_ta  ! day CO2 flux from actual NPP, required for cconservation test
    dsl4jsb_Real2D_onChunk ::  cflux_npp_2_atm_ta       ! grid cell averages of net CO2 fluxes between
    dsl4jsb_Real2D_onChunk ::  cflux_soilresp_2_atm_ta  ! .. biosphere (due to NPP, soil respiration and
    dsl4jsb_Real2D_onChunk ::  cflux_herb_2_atm_ta      ! .. grazing) and atmosphere [kg(CO2)/(m^2(ground) s)]


    ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
    !dsl4jsb_Real2D_onChunk :: fract_litter_wood_new
    !
    dsl4jsb_Real2D_onChunk :: veg_fract_correction
    dsl4jsb_Real2D_onChunk :: LAI
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: t_air                       ! Atmosphere temperature (lowest layer) in Kelvin!
    dsl4jsb_Real2D_onChunk :: rain
    dsl4jsb_Real2D_onChunk :: snow


    !PHS: root depth
    soil_w => Get_vgrid('soil_depth_water')
    
    
    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(CARBON_)) RETURN

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc   = options%nc
    dtime = options%dtime

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_config(CARBON_)
    dsl4jsb_Get_memory(A2L_)
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! JeDi needs soil moisture from HYDRO process
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_memory(RAD_)
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    
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

    ! This is on canopy area:
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cconservation_calcCpools )  ! out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_greenwood_2_litter ) ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green )               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods )               ! inout
    
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    !
    IF (TRIM(model%config%usecase) == 'jsbach_jedi' ) THEN

    dsl4jsb_Get_var2D_onChunk(CARBON_,   c_storage         )     !     Storage carbon pool 
    dsl4jsb_Get_var2D_onChunk(CARBON_,   c_leaves          )     !      Leaves carbon pool             
    dsl4jsb_Get_var2D_onChunk(CARBON_,   c_fnroot          )     !   Fine root carbon pool
    dsl4jsb_Get_var2D_onChunk(CARBON_,   c_csroot          )     ! Coarse root carbon pool
    dsl4jsb_Get_var2D_onChunk(CARBON_,   c_seed            )     !        Seed carbon pool
    dsl4jsb_Get_var2D_onChunk(CARBON_,   c_stem            )     !        Stem carbon pool
    dsl4jsb_Get_var2D_onChunk(CARBON_,   ff_grow_W         )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   ff_grow_T         )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   ff_grow_G         )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   Flg_grow         )
    dsl4jsb_Get_var2D_onChunk(CARBON_,    Flg_grow_G     )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   ff_npp            )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   flx_germ          )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   adj_fraction      )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   jedi_light_constraint )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   Flg_die           )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   Flg_senesc        )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   alpha_jedi        )
    
    dsl4jsb_Get_var2D_onChunk(CARBON_,   alpha_jedi_seed)
    dsl4jsb_Get_var3D_onChunk(CARBON_,   seed_depth_sl)
    dsl4jsb_Get_var2D_onChunk(CARBON_,   seed_maxmoist)
    dsl4jsb_Get_var2D_onChunk(CARBON_,   w_soil_column_seed)
    dsl4jsb_Get_var2D_onChunk(CARBON_,   w_soil_column_jedi)
    dsl4jsb_Get_var2D_onChunk(CARBON_,  root_maxmoist )
    dsl4jsb_Get_var2D_onChunk(CARBON_,  jedi_pot_supply_transp )
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    max_moist         )

    dsl4jsb_Get_var2D_onChunk(CARBON_,   jedi_root_depth          )
    dsl4jsb_Get_var2D_onChunk(HYDRO_,   soil_depth     )
    ! 3 Dimensional--------------------------------------------------------------------------
    dsl4jsb_Get_var3D_onChunk(HYDRO_,   w_soil_sl                )
    dsl4jsb_Get_var3D_onChunk(CARBON_,   jedi_root_depth_sl ) 
    dsl4jsb_Get_var3D_onChunk(HYDRO_,   vol_field_cap_sl       ) 
    dsl4jsb_Get_var3D_onChunk(HYDRO_,  soil_depth_sl            )
    ! 3 Dimensional--------------------------------------------------------------------------
    dsl4jsb_Get_var2D_onChunk(CARBON_,  one_counter        )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   MaintRes          )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   MaintRes_sum      )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   MaintRes_yDayMean )
    dsl4jsb_Get_var2D_onChunk(CARBON_,    temp_sum         )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   temp_yDayMean  )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   GrowthResp        )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   AutoRes           )
    dsl4jsb_Get_var2D_onChunk(CARBON_,    CUE              )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   C_balan           )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   C_delta           )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   c_total           )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   ff_dTdt           )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   ff_T_30d          )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   pseudo_c_total     )
    dsl4jsb_Get_var2D_onChunk(CARBON_,   ff_grow_T_yDayMean) 
    dsl4jsb_Get_var2D_onChunk(RAD_,        faPAR_prev_noon)

    END IF
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
   
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

    dsl4jsb_Get_var2D_onChunk(CARBON_,  excess_NPP)                 ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  root_exudates )             ! inout

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_sum_litter_ag_ta )         ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_sum_litter_bg_ta )         ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_sum_humus_ta )             ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_sum_natural_ta )           ! out

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
    IF (debug_on() .AND. iblk==1 .AND. is_newday(options%current_datetime, dtime)  ) CALL message(TRIM(routine), &
                                                           'update_C_NPP_pot_allocation: First time step of new day')
    IF (debug_on() .AND. iblk==1 .AND. is_newyear(options%current_datetime, dtime) ) CALL message(TRIM(routine), &
                                                           'update_C_NPP_pot_allocation: First time step of new year')

    ! Initializations
    cflux_npp_2_atm_yday_ta  = 0._wp
    cflux_npp_2_atm_ta       = 0._wp
    cflux_soilresp_2_atm_ta  = 0._wp
    cflux_herb_2_atm_ta = 0._wp

    IF( is_newyear(options%current_datetime,dtime) .OR. lstart) THEN ! reset annual sums if newyear or if just started exp
      c_decomp_humus_1_sum =  0._wp
      c_decomp_humus_2_sum =  0._wp
      c_into_humus_1_sum =  0._wp
      c_into_humus_2_sum =  0._wp
    ENDIF

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

    lstart = is_time_experiment_start(options%current_datetime)
    

    IF (TRIM(model%config%usecase) == 'jsbach_jedi') THEN

    inidate => get_time_experiment_start()
    iniyear = get_year(inidate)
    CALL deallocateDatetime(inidate)
    lrestart= jsbach_is_restarted()      
    lspinup_jedi = dsl4jsb_Config(CARBON_)%l_spinup_jedi
    
    IF (get_year(options%current_datetime) .LE. iniyear + 2 ) THEN    !  PHS: initial carbon
      ff_npp = 0._wp
      Flg_senesc = 0._wp
      Flg_die = 0._wp
      one_counter = 0._wp           
    END IF     !  PHS: initial carbon
    
      seed_depth_sl = jedi_soil_depth_to_layers_1d(  SPREAD(0.05_wp, NCOPIES=nc, DIM=1)  , soil_w%dz(:))         ! Jedi assumes germination of seed depends on the first 50mm of soil wetness
     CALL get_water_in_root_zone(vol_field_cap_sl(:,:) * soil_depth_sl(:,:),  soil_depth_sl(:,:), seed_depth_sl(:,:),     seed_maxmoist(:))    
     CALL get_water_in_root_zone(vol_field_cap_sl(:,:) * soil_depth_sl(:,:), soil_depth_sl(:,:), jedi_root_depth_sl(:,:), root_maxmoist(:))

     CALL get_water_in_root_zone(w_soil_sl(:,:), soil_depth_sl(:,:), seed_depth_sl(:,:), w_soil_column_seed(:))
     CALL get_water_in_root_zone(w_soil_sl(:,:), soil_depth_sl(:,:), jedi_root_depth_sl(:,:), w_soil_column_jedi(:))
      
    alpha_jedi(:)          = jedi_water_stress_factor( w_soil_column_jedi(:),      root_maxmoist(:),  dsl4jsb_Config(HYDRO_)%w_soil_wilt_fract)
    alpha_jedi_seed(:) = jedi_water_stress_factor( w_soil_column_seed(:), seed_maxmoist(:),  dsl4jsb_Config(HYDRO_)%w_soil_wilt_fract)            
    
    ! 2019.11.27 PHS: still at the stage of testing "water stress" 
    WHERE (c_leaves .gt. 0._wp)  
      ! The potential supply rate from soil depends on the root-shoot-ratio. The root-shoot-ratio is the ratio betweem fine-root-carbon and Leaf-area-index.
      !     a(r)=1-exp(- (-log(0.01))/r_optimal   * r  )     a : is a function of root-shoot-ratio  [0...1]
      ! r : root-shoot-ratio [10^-3  ~  10^3]
      ! r_optimal: the median value of the root-shoot-ratio among all PGS
      jedi_pot_supply_transp = (1._wp - exp( -(- log(0.01))/rootshoot_opt * &
        &                         (c_fnroot / (c_leaves*dsl4jsb_Lctlib_param(JeDiSLA)) ) )) * max_pipe_transpir * alpha_jedi
    ELSEWHERE 
      jedi_pot_supply_transp =  0._wp
    END WHERE
    
      ff_grow_T =  ( t_air - tmelt ) / dsl4jsb_Lctlib_param(jedi_N_pseudo_T)      + &
        &          ff_grow_T * dsl4jsb_Lctlib_param(jedi_F_pseudo_T)
      ff_grow_G =  alpha_jedi_seed   / dsl4jsb_Lctlib_param(jedi_N_pseudo_soilM)  + &
        &          ff_grow_G * dsl4jsb_Lctlib_param(jedi_F_pseudo_soilM)
      ff_grow_W =  alpha_jedi        / dsl4jsb_Lctlib_param (jedi_N_pseudo_soilM) + &
        &          ff_grow_W * dsl4jsb_Lctlib_param(jedi_F_pseudo_soilM)
      ff_npp    =  NPP_act_yDayMean  / dsl4jsb_Lctlib_param(jedi_N_pseudo_npp)    + &
        &          ff_npp    * dsl4jsb_Lctlib_param(jedi_F_pseudo_npp)
      
      ff_T_30d  =  ( t_air - tmelt ) / ( 1._wp / (1._wp-EXP(- dtime / (30._wp* 86400._wp )))) + &
        &          ff_T_30d * EXP(- dtime / (30._wp* 86400._wp ))
    END IF

    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!     

    ! All per tile variables need to be re-calculated in case they are written to output
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
      & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:), root_exudates = root_exudates(:),   &
      & soil_respiration = soil_respiration(:), NPP_flux_correction = NPP_flux_correction(:),   &
      & cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter(:),                             &
      & cflux_c_green_2_herb = cflux_c_green_2_herb(:), cflux_herb_2_atm = cflux_herb_2_atm(:), &
      & NPP_act_yDayMean = NPP_act_yDayMean(:), NPP_pot_yDayMean = NPP_pot_yDayMean(:), GPP_yDayMean = GPP_yDayMean(:))

    IF( .NOT. is_newday(options%current_datetime,dtime) .OR. lstart) THEN ! perform daily sums if WE ARE NOT STARTING A NEW DAY
                                            ! or if we JUST STARTED AN TOTALLY NEW EXPERIMENT

      LAI_sum          = LAI_sum          + LAI
      NPP_pot_sum      = NPP_pot_sum      + NPP_pot_rate_ca
      GPP_sum          = GPP_sum          + gross_assimilation_ca
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!     
      ! PHS: Compute Maintenance respration for JeDi. This is computed outside clac_tradeoff to account for 
      !   the average of previous-day value.
      !   MaintRes_sum     = MaintRes_sum + MaintRes
      IF (TRIM(model%config%usecase) == 'jsbach_jedi' ) temp_sum = temp_sum + t_air
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!     

    ELSE                                    ! A NEW DAY BEGINS and we DID NOT START AN
                                            ! TOTALLY NEW EXPERIMENT ==> perform carbon balance

      ! First save as "previous previous days" and "previous days" means
      LAI_yyDayMean  = LAI_yDayMean
      LAI_yDayMean   = LAI_sum/timesteps_per_day(dtime)

      NPP_pot_yDayMean     = NPP_pot_sum/timesteps_per_day(dtime)
      GPP_yDayMean         = GPP_sum/timesteps_per_day(dtime)

      ! Then restart summing of this days values
      LAI_sum             = LAI
      NPP_pot_sum         = NPP_pot_rate_ca
      GPP_sum             = gross_assimilation_ca
      !NPP_pot_rate_ca_acc        = NPP_pot_rate_ca/timesteps_per_day(dtime) ! s.o.

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!     
      ! PHS: Compute Maintenance respration for JeDi. 
      ! PHS try to compute maintenance respiration only once per day
      IF (TRIM(model%config%usecase) == 'jsbach_jedi' ) THEN

       temp_yDayMean  = temp_sum / timesteps_per_day(dtime)
       
       MaintRes = 0._wp
       MaintRes = c_sapwood * c_CN_woody * ( c_stem + c_csroot ) + &
         &        dsl4jsb_Lctlib_param(NitroStatus) * dsl4jsb_Lctlib_param(JeDiSLA) * &
         &             ( c_leaves + c_fnroot )  !NNN
       MaintRes = MaintRes  * c_Res_N    * ( 1.6_wp ** ( ( temp_yDayMean - tmelt - 20._wp ) / 10._wp) ) ! JeDi uses 1.4 for Q10 value
    
       MaintRes = MaintRes / sec_per_day ! The maintanence respiration is daily computed in JeDi.
       MaintRes_yDayMean   = MaintRes
       
       temp_sum = t_air
       
       !PHS 05.08.2019 Tend
       !PHS  I think I made mistake for calculating ff_dTdt 
       !ff_dTdt            = ( ff_T_30d - ff_grow_T_yDayMean ) / ( 1._wp /  ( 1._wp - EXP(- dtime / (30._wp* 86400._wp ))  )) + &
       !  &                  EXP(- dtime / (30._wp* 86400._wp ))  * ff_dTdt  
       ff_dTdt            = ( ff_T_30d - ff_grow_T_yDayMean ) / ( 1._wp /  ( 1._wp - EXP(- 86400._wp / (30._wp* 86400._wp ))  )) + &
         &                  EXP(- 86400._wp / (30._wp* 86400._wp ))  * ff_dTdt   
       ff_grow_T_yDayMean = ff_T_30d

      END IF

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

      !Calculate current sum before operation for c conservation test
      CALL calculate_current_c_ta_state_sum(tile, options, old_c_state_sum(:))

      ! Now update the C pools
      IF (TRIM(model%config%usecase) .eq. 'jsbach_jedi' ) THEN
       
       
       IF (get_year(options%current_datetime) .GT. iniyear + 2 ) THEN     !  PHS: initial carbon
       
        CALL calc_Ctradeoff( &
        & GPP_yDayMean,                                & ! in
        ! JeDi plant trait parameter comes in now <<<<<<<<<<<
        & dsl4jsb_Lctlib_param(GrowRspdMoist),         & ! in
        & dsl4jsb_Lctlib_param(GrowRspdTemp ),         & ! in
        & dsl4jsb_Lctlib_param(CritiGrowTemp),         & ! in
        & dsl4jsb_Lctlib_param(EndGrowTemp),           & ! in
        & dsl4jsb_Lctlib_param(GermFrac     ),         & ! in
        & dsl4jsb_Lctlib_param(Allo2Reproduc),         & ! in
        & dsl4jsb_Lctlib_param(Allo2Abground),         & ! in
        & dsl4jsb_Lctlib_param(Allo2Blground),         & ! in
        & dsl4jsb_Lctlib_param(Allo2Storage ),         & ! in
        & dsl4jsb_Lctlib_param(RelAllo2Abgrd),         & ! in
        & dsl4jsb_Lctlib_param(RelAllo2Blgrd),         & ! in
        & dsl4jsb_Lctlib_param(TurnoverStruc),         & ! in
        & dsl4jsb_Lctlib_param(TurnoverLeaf ),         & ! in
        & dsl4jsb_Lctlib_param(SenesRspdNPP ),         & ! in
        & dsl4jsb_Lctlib_param(NitroStatus  ),         & ! in
        & dsl4jsb_Lctlib_param(JeDiSLA),               &  ! in
        & dsl4jsb_Lctlib_param(LandcoverClass),        &  ! in
        ! JeDi plant trait parameter comes in now >>>>>>>>>>>
        & lspinup_jedi,                                & !in 
        & iniyear,                                     & !in
        & get_year(options%current_datetime),          & !in
        & timesteps_per_day(dtime),                    & ! in
        & faPAR_prev_noon,                                  & !in
        & one_counter,                                 & ! inout
        & c_storage,                                   & ! inout
        & c_leaves,                                    & ! inout
        & c_fnroot,                                    & ! inout
        & c_csroot,                                    & ! inout
        & c_seed,                                      & ! inout
        & c_stem,                                      & ! inout
        & ff_grow_W,                               & ! in
        & ff_grow_T,                                   & ! in
        & ff_grow_G,                                   & ! in
        & ff_dTdt,                                       & !in
        & Flg_grow,                    & ! out
        & Flg_grow_G ,              & ! out
        & ff_npp,                                      & ! inout
        & Flg_die,                                     & ! inout
        & Flg_senesc,                                  & ! out
        & flx_germ,                                    & ! out
        &  adj_fraction,                                  &  !out
        & jedi_light_constraint,                       &  !out
        & MaintRes_yDayMean,                           & ! inout
        & GrowthResp,                                  & ! out
        & AutoRes,                                     & ! out
         & CUE,                          & !out carbon use efficiency
         & jedi_root_depth,                 &  ! out 
        & C_balan,                                     & ! out
        & C_delta,                                     & ! out
        & c_total,                                     & !
        &  c_Res_N,                    & !
        &  c_CN_woody,               & !
        &  c_sapwood,                   & ! 
         &  temp_yDayMean,        & ! in
  !
 !
        & soil_respiration,                            & ! out
 !
        & cflux_c_greenwood_2_litter,                  & ! out
 !
        & NPP_act_yDayMean,                            & ! out
 !
    ! Yasso variables
        &  temp2_30d        =  pseudo_temp_yDay,                      & ! in, optional, necessary with yasso
        &  precip_30d       =  pseudo_precip_yDay,                    & ! in, optional, necessary with yasso
                             !
        &  c_acid_ag1       =  c_acid_ag1,                            & ! inout
        &  c_water_ag1      =  c_water_ag1,                           & ! inout
        &  c_ethanol_ag1    =  c_ethanol_ag1,                         & ! inout
        &  c_nonsoluble_ag1 =  c_nonsoluble_ag1,                      & ! inout
        &  c_acid_bg1       =  c_acid_bg1,                            & ! inout
        &  c_water_bg1      =  c_water_bg1,                           & ! inout
        &  c_ethanol_bg1    =  c_ethanol_bg1,                         & ! inout
        &  c_nonsoluble_bg1 =  c_nonsoluble_bg1,                      & ! inout
        &  c_humus_1        =  c_humus_1,                             & ! inout
        &  c_acid_ag2       =  c_acid_ag2,                            & ! inout
        &  c_water_ag2      =  c_water_ag2,                           & ! inout
        &  c_ethanol_ag2    =  c_ethanol_ag2,                         & ! inout
        &  c_nonsoluble_ag2 =  c_nonsoluble_ag2,                      & ! inout
        &  c_acid_bg2       =  c_acid_bg2,                            & ! inout
        &  c_water_bg2      =  c_water_bg2,                           & ! inout
        &  c_ethanol_bg2    =  c_ethanol_bg2,                         & ! inout
        &  c_nonsoluble_bg2 =  c_nonsoluble_bg2,                      & ! inout
        &  c_humus_2        =  c_humus_2,                             & ! inout
                              !
        &  LeafLit_coef1    =  dsl4jsb_Lctlib_param(LeafLit_coef(1)), & ! in
        &  LeafLit_coef2    =  dsl4jsb_Lctlib_param(LeafLit_coef(2)), & ! in
        &  LeafLit_coef3    =  dsl4jsb_Lctlib_param(LeafLit_coef(3)), & ! in
        &  LeafLit_coef4    =  dsl4jsb_Lctlib_param(LeafLit_coef(4)), & ! in
        &  LeafLit_coef5    =  dsl4jsb_Lctlib_param(LeafLit_coef(5)), & ! in
        &  WoodLit_coef1    =  dsl4jsb_Lctlib_param(WoodLit_coef(1)), & ! in
        &  WoodLit_coef2    =  dsl4jsb_Lctlib_param(WoodLit_coef(2)), & ! in
        &  WoodLit_coef3    =  dsl4jsb_Lctlib_param(WoodLit_coef(3)), & ! in
        &  WoodLit_coef4    =  dsl4jsb_Lctlib_param(WoodLit_coef(4)), & ! in
        &  WoodLit_coef5    =  dsl4jsb_Lctlib_param(WoodLit_coef(5)), & ! in
        &  WoodLitterSize   =  dsl4jsb_Lctlib_param(WoodLitterSize)   & ! in
        & )
      END IF  !  PHS: initial carbon
      jedi_root_depth = MAX( jedi_root_depth, 0.05_wp)
      
      jedi_root_depth = MIN(jedi_root_depth,soil_depth)
      
      jedi_root_depth_sl = jedi_soil_depth_to_layers_1d(  jedi_root_depth, soil_w%dz(:))   

      pseudo_c_total  =  c_total / ( 1._wp / (1._wp-EXP(- 86400._wp / (3650._wp* 86400._wp )))) + &
      &          pseudo_c_total * EXP(- 86400._wp / (3650._wp* 86400._wp ))

      
      ELSE
    
        CALL calc_Cpools( &
          &  LAI_yDayMean,                                      & ! in
          &  LAI_yyDayMean,                                     & ! in
          &  MaxLai,                                            & ! in
          &  NPP_pot_yDayMean,                                  & ! in
          !
          &  dsl4jsb_Lctlib_param(fract_npp_2_woodPool),        & ! in
          &  dsl4jsb_Lctlib_param(fract_NPP_2_reservePool),     & ! in
          &  dsl4jsb_Lctlib_param(fract_NPP_2_exudates),        & ! in
          !
          &  dsl4jsb_Lctlib_param(fract_green_2_herbivory),     & ! in
          &  dsl4jsb_Lctlib_param(tau_c_woods),                 & ! in
          &  dsl4jsb_Lctlib_param(LAI_shed_constant),           & ! in
          &  dsl4jsb_Lctlib_param(Max_C_content_woods),         & ! in
          &  dsl4jsb_Lctlib_param(specificLeafArea_C),          & ! in
          &  dsl4jsb_Lctlib_param(reserveC2leafC),              & ! in
          &  dsl4jsb_Lctlib_param(PhenologyType),               & ! in  R: is it crop? If yes then this is 5!
          !
          &  c_green,                                           & ! inout
          &  c_woods,                                           & ! inout
          &  c_reserve,                                         & ! inout
          &  c_crop_harvest,                                    & ! inout
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
          !
          &  temp2_30d        =  pseudo_temp_yDay,                      & ! in
          &  precip_30d       =  pseudo_precip_yDay,                    & ! in
                              !
          &  c_acid_ag1       =  c_acid_ag1,                            & ! inout
          &  c_water_ag1      =  c_water_ag1,                           & ! inout
          &  c_ethanol_ag1    =  c_ethanol_ag1,                         & ! inout
          &  c_nonsoluble_ag1 =  c_nonsoluble_ag1,                      & ! inout
          &  c_acid_bg1       =  c_acid_bg1,                            & ! inout
          &  c_water_bg1      =  c_water_bg1,                           & ! inout
          &  c_ethanol_bg1    =  c_ethanol_bg1,                         & ! inout
          &  c_nonsoluble_bg1 =  c_nonsoluble_bg1,                      & ! inout
          &  c_humus_1        =  c_humus_1,                             & ! inout
          &  c_acid_ag2       =  c_acid_ag2,                            & ! inout
          &  c_water_ag2      =  c_water_ag2,                           & ! inout
          &  c_ethanol_ag2    =  c_ethanol_ag2,                         & ! inout
          &  c_nonsoluble_ag2 =  c_nonsoluble_ag2,                      & ! inout
          &  c_acid_bg2       =  c_acid_bg2,                            & ! inout
          &  c_water_bg2      =  c_water_bg2,                           & ! inout
          &  c_ethanol_bg2    =  c_ethanol_bg2,                         & ! inout
          &  c_nonsoluble_bg2 =  c_nonsoluble_bg2,                      & ! inout
          &  c_humus_2        =  c_humus_2,                             & ! inout
                                !
          &  LeafLit_coef1    =  dsl4jsb_Lctlib_param(LeafLit_coef(1)), & ! in
          &  LeafLit_coef2    =  dsl4jsb_Lctlib_param(LeafLit_coef(2)), & ! in
          &  LeafLit_coef3    =  dsl4jsb_Lctlib_param(LeafLit_coef(3)), & ! in
          &  LeafLit_coef4    =  dsl4jsb_Lctlib_param(LeafLit_coef(4)), & ! in
          &  LeafLit_coef5    =  dsl4jsb_Lctlib_param(LeafLit_coef(5)), & ! in
          &  WoodLit_coef1    =  dsl4jsb_Lctlib_param(WoodLit_coef(1)), & ! in
          &  WoodLit_coef2    =  dsl4jsb_Lctlib_param(WoodLit_coef(2)), & ! in
          &  WoodLit_coef3    =  dsl4jsb_Lctlib_param(WoodLit_coef(3)), & ! in
          &  WoodLit_coef4    =  dsl4jsb_Lctlib_param(WoodLit_coef(4)), & ! in
          &  WoodLit_coef5    =  dsl4jsb_Lctlib_param(WoodLit_coef(5)), & ! in
          &  WoodLitterSize   =  dsl4jsb_Lctlib_param(WoodLitterSize),  & ! in

          & c_decomp_humus_1_sum  = c_decomp_humus_1_sum,               & ! out
          & c_decomp_humus_2_sum  = c_decomp_humus_2_sum,               & ! out
          & c_into_humus_1_sum    = c_into_humus_1_sum,                 & ! out
          & c_into_humus_2_sum    = c_into_humus_2_sum                  & ! out

          & )
      END IF ! usecase=jedi

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

      CALL calculate_current_c_ag_1_and_bg_sums(tile, options)
      CALL recalc_per_tile_vars(tile, options,                                      &
        & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:),      &
        ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
        ! Variavles for JeDi carbon allocation 
        & c_leaves = c_leaves(:) , c_storage = c_storage(:) , c_fnroot = c_fnroot(:) , &
        & c_csroot = c_csroot(:) , c_stem = c_stem(:)       , c_seed = c_seed(:)     , &
        & c_total = c_total(:),         &
        ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
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
        & c_into_humus_2_sum = c_into_humus_2_sum(:), root_exudates = root_exudates(:),           &
        & soil_respiration = soil_respiration(:), NPP_flux_correction = NPP_flux_correction(:),   &
        & cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter(:),                             &
        & cflux_c_green_2_herb = cflux_c_green_2_herb(:), cflux_herb_2_atm = cflux_herb_2_atm(:), &
        & NPP_act_yDayMean = NPP_act_yDayMean(:), NPP_pot_yDayMean = NPP_pot_yDayMean(:), GPP_yDayMean = GPP_yDayMean(:))

      ! cflux_herb_2_atm and soil_respiration are negative (fluxes away from land)
      ! thus + here instead of - to reduce NPP by these fluxes
      current_fluxes = (NPP_act_yDayMean + cflux_herb_2_atm + soil_respiration) &
        & * sec_per_day * veg_fract_correction * fract_fpc_max
      CALL check_carbon_conservation(tile, options, old_c_state_sum(:), &
        & current_fluxes(:), cconservation_calcCpools(:))

    END IF ! new day

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

    cflux_npp_2_atm_yday_ta = molarMassCO2_kg *                 & ! daily C conservation test cannot deal with 
              (- veg_fract_correction * fract_fpc_max           & ! diurnal cycle of NPP,
                   * (NPP_act_yDayMean))                          ! therefore here additionaly the day mean

    ! Calculate diagnostic carbon sums
    CALL calculate_current_c_ta_state_sum(tile, options, c_sum_natural_ta(:), &
      & c_sum_litter_ag_ta(:), c_sum_litter_bg_ta(:), c_sum_humus_ta(:))

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')    

  END SUBROUTINE update_C_NPP_pot_allocation

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
      !PHS 24.07.2020 Compute water stress in CARBON at PFT level. 
  
! -------------------------------------------------------------------------------------------------------
  !>
  !! @brief Implementation of "update" for task "water_stress"
  !! Task "update_water_stress" calculates water stress
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_water_stress(tile, options)
    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_carbon_constants,  ONLY: max_pipe_transpir

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(HYDRO_)

    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(CARBON_)   !PHS

 
    dsl4jsb_Real2D_onChunk :: water_stress
    dsl4jsb_Real2D_onChunk :: evapopot
    dsl4jsb_Real2D_onChunk :: jedi_pot_supply_transp  !PHS
    dsl4jsb_Real2D_onChunk :: alpha_jedi   !PHS	
   
    INTEGER :: iblk, ics, ice
    
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_water_stress'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (.NOT. tile%Is_process_active(CARBON_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(HYDRO_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(CARBON_)

  
    dsl4jsb_Get_var2D_onChunk(CARBON_,    water_stress)        ! out   !PHS
    
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapopot) !PHS
    dsl4jsb_Get_var2D_onChunk(CARBON_,   jedi_pot_supply_transp) !PHS	
    dsl4jsb_Get_var2D_onChunk(CARBON_,  alpha_jedi)!PHS


    IF (TRIM(model%config%usecase) .eq. 'jsbach_jedi' ) THEN
      
      WHERE (evapopot(:) .GE. 0._wp)  
        water_stress(:) = 0._wp
      ELSEWHERE ( evapopot(:) .LT. 0._wp )   
       
        water_stress(:) = MIN(1._wp , ABS(evapopot(:)) / max_pipe_transpir  * &
        &                (1._wp - exp( - jedi_pot_supply_transp(:) / ABS(evapopot(:)) ) ))

      END WHERE
      
  
    END IF 

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_water_stress
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
    


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

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_C_NPP_pot_allocation'

    INTEGER  :: iblk , ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(CARBON_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yDayMean,          weighted_by_fract)

    !@todo Those commented out are currently not calculated
    ! non-yasso-pools
    dsl4jsb_Aggregate_onChunk(CARBON_, c_green_ta,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_woods_ta,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_reserve_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_crop_harvest_ta,     weighted_by_fract)
    
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    dsl4jsb_Aggregate_onChunk(CARBON_, c_leaves_ta,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_stem_ta,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_seed_ta,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_storage_ta,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_fnroot_ta,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_csroot_ta,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_total_ta,                weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, MaintRes_yDayMean, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_,  GrowthResp,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_,  AutoRes,                  weighted_by_fract)
    !dsl4jsb_Aggregate_onChunk(CARBON_,  root_maxmoist,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_,  jedi_root_depth,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_,  jedi_root_depth_sl,        weighted_by_fract)
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    
    ! Variavles for JeDi carbon allocation 
    
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
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles which are interested to JeDi
    dsl4jsb_Aggregate_onChunk(CARBON_, c_ag_sum_1,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, c_bg_sum,              weighted_by_fract)
    
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    
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
    dsl4jsb_Aggregate_onChunk(CARBON_, excess_NPP_ta,            weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yyDayMean_ta,         weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CARBON_, LAI_yDayMean_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_pot_yDayMean_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, NPP_act_yDayMean_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, GPP_yDayMean_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, root_exudates_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_c_green_2_herb_ta,  weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_, cflux_herb_2_littergreen_ta, weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(CARBON_, fract_litter_wood_new_ta,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cconservation_calcCpools, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_npp_2_atm_ta,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_soilresp_2_atm_ta,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_herb_2_atm_ta,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CARBON_, cflux_npp_2_atm_yday_ta,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_C_NPP_pot_allocation

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
      !PHS 24.07.2020 Compute water stress in CARBON at PFT level. 
  
 ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "water_stress"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_water_stress(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(CARBON_)
    
    TYPE(t_jsb_model),       POINTER :: model
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_water_stress'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(CARBON_)
  
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    ! HYDRO_
    dsl4jsb_Aggregate_onChunk(CARBON_, water_stress,        weighted_by_fract)


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_water_stress
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
  

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
    CHARACTER(len=*), PARAMETER :: routine = modname//':rescale_carbon_upon_reference_area_change'
    INTEGER                     :: ics, ice, iblk

    dsl4jsb_Real2D_onChunk ::  c_green
    dsl4jsb_Real2D_onChunk ::  c_reserve
    dsl4jsb_Real2D_onChunk ::  c_woods
    dsl4jsb_Real2D_onChunk ::  c_crop_harvest
    
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    !
    dsl4jsb_Real2D_onChunk ::  c_storage        
    dsl4jsb_Real2D_onChunk ::  c_leaves             
    dsl4jsb_Real2D_onChunk ::  c_fnroot          
    dsl4jsb_Real2D_onChunk ::  c_csroot       
    dsl4jsb_Real2D_onChunk ::  c_seed          
    dsl4jsb_Real2D_onChunk ::  c_stem
    dsl4jsb_Real2D_onChunk ::  c_total
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
    
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
    
     ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
      ! Variavles for JeDi carbon allocation 
     
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_leaves  )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_storage )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_stem    )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_fnroot  )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_csroot  )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_seed    )             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_total   )             ! inout
     ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    
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

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_1_sum) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_2_sum) ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_1_sum)   ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_2_sum)   ! inout            

    c_green           = c_green           * (oldRefArea / newRefArea)
    c_woods           = c_woods           * (oldRefArea / newRefArea)
    c_reserve         = c_reserve         * (oldRefArea / newRefArea)

    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
         
    c_leaves          = c_leaves          * (oldRefArea / newRefArea)
    c_seed            = c_seed            * (oldRefArea / newRefArea)
    c_stem            = c_stem            * (oldRefArea / newRefArea)
    c_storage         = c_storage         * (oldRefArea / newRefArea)
    c_fnroot          = c_fnroot          * (oldRefArea / newRefArea)
    c_csroot          = c_csroot          * (oldRefArea / newRefArea)
    c_total           = c_total           * (oldRefArea / newRefArea)
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
     
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

    c_decomp_humus_1_sum = c_decomp_humus_1_sum * (oldRefArea / newRefArea)
    c_decomp_humus_2_sum = c_decomp_humus_2_sum * (oldRefArea / newRefArea)
    c_into_humus_1_sum   = c_into_humus_1_sum   * (oldRefArea / newRefArea)
    c_into_humus_2_sum   = c_into_humus_2_sum   * (oldRefArea / newRefArea)

    CALL calculate_current_c_ag_1_and_bg_sums(tile, options)
    CALL recalc_per_tile_vars(tile, options, &
      & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:), &
      ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
      ! Variavles for JeDi carbon allocation 
      & c_leaves = c_leaves(:) , c_storage = c_storage(:) , c_fnroot = c_fnroot(:) , &
      & c_csroot = c_csroot(:) , c_stem = c_stem(:)       , c_seed = c_seed(:)     , &
      & c_total = c_total(:), &
      ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
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

  END SUBROUTINE rescale_carbon_upon_reference_area_change


  ! ================================================================================================================================
  !>
  !! @brief Calculates the per tile states from current canopy states
  !!
  !! @todo discuss how to make more general, e.g. get the info on which state variables to work on from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE recalc_per_tile_vars(tile, options,                  &
      & c_green , c_woods, c_reserve, c_crop_harvest,             &
      ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
      ! Variavles for JeDi carbon allocation 
      & c_leaves , c_storage , c_fnroot , c_csroot , c_stem , c_seed , c_total, &
      ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
      & c_acid_ag1, c_water_ag1, c_ethanol_ag1, c_nonsoluble_ag1, &
      & c_acid_bg1, c_water_bg1, c_ethanol_bg1, c_nonsoluble_bg1, &
      & c_acid_ag2, c_water_ag2, c_ethanol_ag2, c_nonsoluble_ag2, &
      & c_acid_bg2, c_water_bg2, c_ethanol_bg2, c_nonsoluble_bg2, &
      & c_humus_1, c_humus_2,                                     &
      & c_decomp_humus_1_sum, c_decomp_humus_2_sum, c_into_humus_1_sum, c_into_humus_2_sum, &
      & soil_respiration, NPP_flux_correction, NPP_pot_yDayMean, NPP_act_yDayMean,          &
      & GPP_yDayMean, root_exudates, cflux_c_green_2_herb, cflux_herb_2_atm,                &
      & cflux_dist_greenreserve_2_soil, cflux_dist_woods_2_soil, cflux_fire_all_2_atm,      &
      & cflux_c_greenwood_2_litter)

    USE mo_carbon_process,    ONLY: get_per_tile
    !USE mo_carbon_constants,  ONLY: molarMassCO2_kg

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)            :: options
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_green(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_woods(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_reserve(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_crop_harvest(:)
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_leaves(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_stem(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_fnroot(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_csroot(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_seed(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_storage(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_total(:)
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

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

    REAL(wp),                   INTENT(in), OPTIONAL  :: c_decomp_humus_1_sum(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_decomp_humus_2_sum(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_into_humus_1_sum(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_into_humus_2_sum(:)

    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_dist_greenreserve_2_soil(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_dist_woods_2_soil(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_fire_all_2_atm(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: soil_respiration(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: NPP_pot_yDayMean(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: NPP_act_yDayMean(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: NPP_flux_correction(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: GPP_yDayMean(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_c_greenwood_2_litter(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: root_exudates(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_c_green_2_herb(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_herb_2_atm(:)

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(PHENO_)

    ! Local variables
    CHARACTER(len=*), PARAMETER     :: routine = modname//':recalc_per_tile_vars'
    INTEGER                         :: ics, ice, iblk

    dsl4jsb_Real2D_onChunk ::  c_green_ta
    dsl4jsb_Real2D_onChunk ::  c_reserve_ta
    dsl4jsb_Real2D_onChunk ::  c_woods_ta
    dsl4jsb_Real2D_onChunk ::  c_crop_harvest_ta

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    
    dsl4jsb_Real2D_onChunk ::     c_storage_ta
    dsl4jsb_Real2D_onChunk ::     c_seed_ta
    dsl4jsb_Real2D_onChunk ::     c_stem_ta
    dsl4jsb_Real2D_onChunk ::     c_leaves_ta
    dsl4jsb_Real2D_onChunk ::     c_fnroot_ta
    dsl4jsb_Real2D_onChunk ::     c_csroot_ta        
    dsl4jsb_Real2D_onChunk ::     c_total_ta          
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

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

    dsl4jsb_Real2D_onChunk ::  cflux_dist_greenreserve_2_soil_ta
    dsl4jsb_Real2D_onChunk ::  cflux_dist_woods_2_soil_ta
    dsl4jsb_Real2D_onChunk ::  cflux_fire_all_2_atm_ta
    dsl4jsb_Real2D_onChunk ::  soil_respiration_ta
    dsl4jsb_Real2D_onChunk ::  NPP_pot_yDayMean_ta
    dsl4jsb_Real2D_onChunk ::  NPP_act_yDayMean_ta
    dsl4jsb_Real2D_onChunk ::  NPP_flux_correction_ta
    dsl4jsb_Real2D_onChunk ::  GPP_yDayMean_ta
    dsl4jsb_Real2D_onChunk ::  cflux_c_greenwood_2_litter_ta
    dsl4jsb_Real2D_onChunk ::  root_exudates_ta
    dsl4jsb_Real2D_onChunk ::  cflux_c_green_2_herb_ta
    dsl4jsb_Real2D_onChunk ::  cflux_herb_2_atm_ta

    dsl4jsb_Real2D_onChunk ::  veg_fract_correction
    dsl4jsb_Real2D_onChunk ::  fract_fpc_max

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(PHENO_)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_green_ta )               ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_reserve_ta )             ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_woods_ta )               ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_crop_harvest_ta)         ! opt out

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    
    IF (TRIM(model%config%usecase) == 'jsbach_jedi') THEN
      dsl4jsb_Get_var2D_onChunk(CARBON_,  c_leaves_ta  )
      dsl4jsb_Get_var2D_onChunk(CARBON_,  c_fnroot_ta  )
      dsl4jsb_Get_var2D_onChunk(CARBON_,  c_csroot_ta  )
      dsl4jsb_Get_var2D_onChunk(CARBON_,  c_stem_ta    )
      dsl4jsb_Get_var2D_onChunk(CARBON_,  c_storage_ta )
      dsl4jsb_Get_var2D_onChunk(CARBON_,  c_seed_ta    )
      dsl4jsb_Get_var2D_onChunk(CARBON_,  c_total_ta   )
    END IF
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  


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

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_1_sum_ta)    ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_decomp_humus_2_sum_ta)    ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_1_sum_ta)      ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_into_humus_2_sum_ta)      ! opt out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_dist_greenreserve_2_soil_ta ) ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_dist_woods_2_soil_ta )   ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_fire_all_2_atm_ta )      ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  soil_respiration_ta )          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_pot_yDayMean_ta )          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_act_yDayMean_ta )          ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  NPP_flux_correction_ta)        ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  GPP_yDayMean_ta )              ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_greenwood_2_litter_ta )! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  root_exudates_ta )             ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_c_green_2_herb_ta )      ! opt out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  cflux_herb_2_atm_ta )          ! opt out

    dsl4jsb_Get_var2D_onChunk(PHENO_,   veg_fract_correction )       ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_,   fract_fpc_max )              ! in

    IF (PRESENT(c_green)) CALL get_per_tile(c_green_ta(:), c_green(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_woods)) CALL get_per_tile(c_woods_ta(:), c_woods(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_reserve)) CALL get_per_tile(c_reserve_ta(:), c_reserve(:), veg_fract_correction(:), fract_fpc_max(:))

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    
    IF (TRIM(model%config%usecase) == 'jsbach_jedi') THEN
      IF (PRESENT(c_leaves) )  CALL get_per_tile(c_leaves_ta(:) , c_leaves(:) , veg_fract_correction(:), fract_fpc_max(:))    
      IF (PRESENT(c_stem)   )  CALL get_per_tile(c_stem_ta(:)   , c_stem(:)   , veg_fract_correction(:), fract_fpc_max(:))    
      IF (PRESENT(c_fnroot) )  CALL get_per_tile(c_fnroot_ta(:) , c_fnroot(:) , veg_fract_correction(:), fract_fpc_max(:))    
      IF (PRESENT(c_csroot) )  CALL get_per_tile(c_csroot_ta(:) , c_csroot(:) , veg_fract_correction(:), fract_fpc_max(:))
      IF (PRESENT(c_seed)   )  CALL get_per_tile(c_seed_ta(:)   , c_seed(:)   , veg_fract_correction(:), fract_fpc_max(:))   
      IF (PRESENT(c_storage))  CALL get_per_tile(c_storage_ta(:), c_storage(:), veg_fract_correction(:), fract_fpc_max(:))
      IF (PRESENT(c_total)  )  CALL get_per_tile(c_total_ta(:)  , c_total(:)  , veg_fract_correction(:), fract_fpc_max(:))
    END IF
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  


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
    IF (PRESENT(cflux_herb_2_atm)) CALL get_per_tile(cflux_herb_2_atm_ta(:), &
      & cflux_herb_2_atm(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(cflux_c_green_2_herb)) CALL get_per_tile(cflux_c_green_2_herb_ta(:), &
      & cflux_c_green_2_herb(:), veg_fract_correction(:), fract_fpc_max(:))

  END SUBROUTINE recalc_per_tile_vars

  ! ================================================================================================================================
  !>
  !! @brief calculates the current sum of all carbon state variables -> on ta! (per canopy vars are not aggregated!)
  !!
  !! @todo discuss how to make more general, e.g. get the info which are the c state variables from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE calculate_current_c_ta_state_sum(tile, options, current_c_state_sum, &
    & c_sum_litter_ag, c_sum_litter_bg, c_sum_humus)

    !USE mo_util,                ONLY: real2string

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)         :: options
    REAL(wp),                   INTENT(out)        :: current_c_state_sum(:)
    REAL(wp), OPTIONAL,         INTENT(out)        :: c_sum_litter_ag(:)
    REAL(wp), OPTIONAL,         INTENT(out)        :: c_sum_litter_bg(:)
    REAL(wp), OPTIONAL,         INTENT(out)        :: c_sum_humus(:)

    dsl4jsb_Def_memory(CARBON_)

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//':calculate_current_c_ta_state_sum'
    INTEGER                     :: ics, ice, iblk


    dsl4jsb_Real2D_onChunk ::  c_green_ta
    dsl4jsb_Real2D_onChunk ::  c_reserve_ta
    dsl4jsb_Real2D_onChunk ::  c_woods_ta
    dsl4jsb_Real2D_onChunk ::  c_crop_harvest_ta

!! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
!    ! Variavles for JeDi carbon allocation 
!    dsl4jsb_Real2D_onChunk ::     c_storage_ta
!    dsl4jsb_Real2D_onChunk ::     c_seed_ta
!    dsl4jsb_Real2D_onChunk ::     c_stem_ta
!    dsl4jsb_Real2D_onChunk ::     c_leaves_ta
!    dsl4jsb_Real2D_onChunk ::     c_fnroot_ta
!    dsl4jsb_Real2D_onChunk ::     c_csroot_ta                      
!! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

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
    
!! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
!    ! Variavles for JeDi carbon allocation 
!    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_leaves_ta  )
!    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_fnroot_ta  )
!    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_csroot_ta  )
!    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_stem_ta    )
!    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_storage_ta )
!    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_seed_ta    )
!! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

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

    !$claw expand parallel update 
    current_c_state_sum(ics:ice) = c_green_ta(:) + c_woods_ta(:) + c_reserve_ta(:) + c_crop_harvest_ta(:) + &
        & c_acid_ag1_ta(:) + c_water_ag1_ta(:) + c_ethanol_ag1_ta(:) + c_nonsoluble_ag1_ta(:) + &
        & c_acid_ag2_ta(:) + c_water_ag2_ta(:) + c_ethanol_ag2_ta(:) + c_nonsoluble_ag2_ta(:) + &
        & c_acid_bg1_ta(:) + c_water_bg1_ta(:) + c_ethanol_bg1_ta(:) + c_nonsoluble_bg1_ta(:) + &
        & c_acid_bg2_ta(:) + c_water_bg2_ta(:) + c_ethanol_bg2_ta(:) + c_nonsoluble_bg2_ta(:) + &
        & c_humus_1_ta(:)  + c_humus_2_ta(:)

    IF (PRESENT(c_sum_litter_ag)) c_sum_litter_ag(ics:ice) =                                    &
        & c_acid_ag1_ta(:) + c_water_ag1_ta(:) + c_ethanol_ag1_ta(:) + c_nonsoluble_ag1_ta(:) + &
        & c_acid_ag2_ta(:) + c_water_ag2_ta(:) + c_ethanol_ag2_ta(:) + c_nonsoluble_ag2_ta(:)
    IF (PRESENT(c_sum_litter_bg)) c_sum_litter_bg(ics:ice) =                                    &
        & c_acid_bg1_ta(:) + c_water_bg1_ta(:) + c_ethanol_bg1_ta(:) + c_nonsoluble_bg1_ta(:) + &
        & c_acid_bg2_ta(:) + c_water_bg2_ta(:) + c_ethanol_bg2_ta(:) + c_nonsoluble_bg2_ta(:)
    IF (PRESENT(c_sum_humus)) c_sum_humus(ics:ice) = c_humus_1_ta(:)  + c_humus_2_ta(:)

  END SUBROUTINE calculate_current_c_ta_state_sum


  ! ================================================================================================================================
  !>
  !! @brief calculates the current sum of all bg carbon state variables (per canopy)
  !!
  !! @todo discuss how to make more general, e.g. get the info which are the c state variables from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE calculate_current_c_ag_1_and_bg_sums(tile, options)

    !USE mo_util,                ONLY: real2string

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)         :: options

    dsl4jsb_Def_memory(CARBON_)

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//':calculate_current_c_ag_1_and_bg_sums'
    INTEGER                     :: ics, ice, iblk

    dsl4jsb_Real2D_onChunk ::  c_bg_sum
    dsl4jsb_Real2D_onChunk ::  c_ag_sum_1

    dsl4jsb_Real2D_onChunk ::  c_acid_ag1
    dsl4jsb_Real2D_onChunk ::  c_water_ag1
    dsl4jsb_Real2D_onChunk ::  c_ethanol_ag1
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_ag1

    dsl4jsb_Real2D_onChunk ::  c_acid_bg1
    dsl4jsb_Real2D_onChunk ::  c_water_bg1
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg1
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg1
    dsl4jsb_Real2D_onChunk ::  c_humus_1

    dsl4jsb_Real2D_onChunk ::  c_acid_bg2
    dsl4jsb_Real2D_onChunk ::  c_water_bg2
    dsl4jsb_Real2D_onChunk ::  c_ethanol_bg2
    dsl4jsb_Real2D_onChunk ::  c_nonsoluble_bg2
    dsl4jsb_Real2D_onChunk ::  c_humus_2

    dsl4jsb_Get_memory(CARBON_)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_bg_sum )             ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ag_sum_1 )           ! out

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1)            ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1)         ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1)      ! in

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg1)            ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg1)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg1)         ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg1)      ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_1)             ! in

    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_bg2)            ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_bg2)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_bg2)         ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_bg2)      ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_humus_2)             ! in

    !TODO-JN-COHORTS: could be done in mo_carbon_process
    !$claw expand parallel update
    c_bg_sum(:) = c_acid_bg1(:) + c_water_bg1(:) + c_ethanol_bg1(:) + c_nonsoluble_bg1(:) + &
             & c_acid_bg2(:) + c_water_bg2(:) + c_ethanol_bg2(:) + c_nonsoluble_bg2(:) + &
             & c_humus_1(:)  + c_humus_2(:)

    c_ag_sum_1(:) =  c_acid_ag1(:) + c_water_ag1(:) + c_ethanol_ag1(:) + c_nonsoluble_ag1(:) + c_humus_1(:)
    !$claw end expand
  END SUBROUTINE calculate_current_c_ag_1_and_bg_sums


  ! ================================================================================================================================
  !>
  !! @brief conducts a carbon conservation test on this tile's current states and given old sum + fluxes
  !!        -> on ta! (per canopy vars are not aggregated!)
  !!
  !! @todo discuss how to make more general, e.g. get the info which are the c state variables from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE check_carbon_conservation(tile, options, old_c_state_sum, cflux, Cconserve)

    !USE mo_util,                ONLY: real2string

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)         :: options
    REAL(wp),                   INTENT(in)         :: old_c_state_sum(:)
    REAL(wp),                   INTENT(in)         :: cflux(:)
        ! Note: cflux is not really a flux, but already multiplied with time unit!
    REAL(wp),                   INTENT(inout)      :: Cconserve(:)

    !dsl4jsb_Def_memory(CARBON_)

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//':check_carbon_conservation'
    !INTEGER                     :: ics, ice, iblk

    CALL calculate_current_c_ta_state_sum(tile, options, Cconserve)
    
    !$claw expand parallel update
    Cconserve(:) = Cconserve(:) - old_c_state_sum(:) - cflux(:)

  END SUBROUTINE check_carbon_conservation

  ! ================================================================================================================================
  !>
  !! @brief conducts a box carbon conservation test for the last day, called each new day -- (jsbach_start_timestep)
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

    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid
    TYPE(t_jsb_task_options)        :: options

    REAL(wp), ALLOCATABLE :: C_flux_yday(:)
    INTEGER               :: ics, ice, iblk, nc

    dsl4jsb_Get_memory(L2A_)
    dsl4jsb_Get_memory(CARBON_)

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

!    iblk  = options%iblk
    options%ics = 1
    options%ice = grid%nproma

    ics = options%ics
    ice = options%ice
    nc = ice - ics + 1

    ALLOCATE(C_flux_yday(nc))

    DO iblk = 1, grid%nblks

      CALL model%Set_options(iblk=iblk)
      options%iblk = iblk
      IF (iblk == grid%nblks) THEN
        options%ice = grid%npromz
        ice = grid%npromz
        ! TODO: ics and ice should be retrieved for each iblk (are not necessarily constant for each block in ICON)
      END IF

      dsl4jsb_Get_var2D_onChunk(L2A_,  carbon_conservation_test )     ! inout
      dsl4jsb_Get_var2D_onChunk(L2A_,  yday_c_state_sum )             ! inout

      ! Do not test as long as carbon_conservation_test still carries the initialisation value
      IF (.NOT. ALL(carbon_conservation_test .EQ. -999.0_wp)) THEN
      
        !CO2 flux due to npp, soil respiration and herbivory
        C_flux_yday(ics:ice) = - ((dsl4jsb_var2D_onChunk(CARBON_, cflux_npp_2_atm_yday_ta)     &
          & + dsl4jsb_var2D_onChunk(CARBON_, cflux_soilresp_2_atm_ta)                   &
          & + dsl4jsb_var2D_onChunk(CARBON_, cflux_herb_2_atm_ta)) / molarMassCO2_kg) * sec_per_day

        IF (model%processes(DISTURB_)%p%config%active) THEN
          C_flux_yday(ics:ice) =  C_flux_yday(ics:ice) &
            & - (dsl4jsb_var2D_onChunk(CARBON_, cflux_fire_all_2_atm_ta))
        ENDIF

        CALL calculate_current_c_ta_state_sum(tile, options, carbon_conservation_test)
        carbon_conservation_test = carbon_conservation_test - yday_c_state_sum  - C_flux_yday(ics:ice)

      ELSE
        carbon_conservation_test = 0._wp
      ENDIF

      CALL calculate_current_c_ta_state_sum(tile, options, yday_c_state_sum)

    END DO

    DEALLOCATE(C_flux_yday)

  END SUBROUTINE yday_carbon_conservation_test

  SUBROUTINE jedi_biomass_to_fract(tile, options)

    USE mo_jsb_var_class,   ONLY:t_jsb_var_real2d
    USE mo_jsb_grid,           ONLY: Get_grid

    CLASS(t_jsb_tile_abstract), INTENT(inout)         :: tile
    TYPE(t_jsb_task_options)  , INTENT(in)         :: options
    
    TYPE(t_jsb_var_real2d), POINTER :: var
    CLASS(t_jsb_tile_abstract), POINTER :: child

    dsl4jsb_Real2D_onDomain :: pseudo_c_total

    dsl4jsb_Def_memory(CARBON_)

    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid

    REAL(wp), ALLOCATABLE :: c_total_sum(:,:), fract(:,:)
    INTEGER               :: nblks, nproma
    INTEGER               :: i, no_children

    CHARACTER(len=*), PARAMETER :: routine = modname//':jedi_biomass_to_fract'

    IF (TRIM(tile%name) /= 'veg') CALL finish(routine, 'should only be called on veg tile')

    dsl4jsb_Get_memory(CARBON_)

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    nblks = grid%nblks
    nproma = grid%nproma

    var => dsl4jsb_memory(CARBON_)%pseudo_c_total

    ALLOCATE(c_total_sum(nproma,nblks))
    c_total_sum(:,:) = 0._wp

    ALLOCATE(fract(nproma,nblks))
    CALL tile%Get_fraction(fract=fract(:,:)) ! Fractions of veg tile

    no_children = tile%Get_no_of_children()

    IF (ASSOCIATED(dsl4jsb_memory(CARBON_)%children(1)%p%vars(var%child_idx(1))%p%ptr2d)) THEN

      DO i=1,no_children
        IF (var%child_idx(i) > 0) THEN
          ASSOCIATE (data_child => dsl4jsb_memory(CARBON_)%children(i)%p%vars(var%child_idx(i))%p%ptr2d(:,:) )
            c_total_sum(:,:) = c_total_sum(:,:) + data_child(:,:)
          END ASSOCIATE
        END IF
      END DO
      c_total_sum(:,:) = c_total_sum(:,:) * fract(:,:) ! Fractions of children need to sum up to fraction of parent tile in block below

      child => tile%Get_first_child_tile()
      DO i=1,no_children
        ASSOCIATE (data_child => dsl4jsb_memory(CARBON_)%children(i)%p%vars(var%child_idx(i))%p%ptr2d(:,:) )
          WHERE (c_total_sum(:,:) > 0._wp)
            child%fract(:,:) = data_child(:,:) / c_total_sum(:,:)
          END WHERE
        END ASSOCIATE
        child => child%Get_next_sibling_tile()
      END DO

    ELSE IF (ASSOCIATED(var%ptr2d_t)) THEN

      IF (SIZE(var%ptr_t, 2) /= no_children) CALL finish(routine, 'Dimension mismatch in c_total%ptr_t and no_children')

      c_total_sum(:,:) = SUM(var%ptr2d_t(:,:,:), DIM=2)
      child => tile%Get_first_child_tile()
      DO i=1,no_children
        WHERE (c_total_sum(:,:) > 0._wp)
          child%fract(:,:) = var%ptr2d_t(:,i,:) / c_total_sum(:,:) * fract(:,:) ! Fractions of children need to sum up to fraction of parent tile
        END WHERE
        child => child%Get_next_sibling_tile()
      END DO

    ELSE

      CALL finish(routine, 'Could not get children of variable c_total')

    END IF
    
    tile%l_fract_children_changed(:,:) = .TRUE.

    DEALLOCATE(c_total_sum)

  END SUBROUTINE jedi_biomass_to_fract

#endif
END MODULE mo_carbon_interface
