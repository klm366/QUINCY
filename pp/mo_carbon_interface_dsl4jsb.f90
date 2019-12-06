!
! This file was automatically generated from:
!    /scratch/snx3000/rschnur/icon-rc-testing-plus-mc8.test/externals/jsbach/src/carbon/mo_carbon_interface.f90
!    with md5 hash: bc9974fe32f3a4c05a211911dbdf469b
!
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

  ! Use of processes in this module (Get_carbon_memory and Get_carbon_config)
  USE mo_jsb_process_class, ONLY: carbon_, assimi_, pheno_, disturb_, a2l_, l2a_

  ! Use of process configurations
  USE mo_pheno_config_class, ONLY: t_pheno_config

  ! Use of process memories (t_carbon_memory)
  USE mo_a2l_memory_class, ONLY: t_a2l_memory
  USE mo_assimi_memory_class, ONLY: t_assimi_memory
  USE mo_pheno_memory_class, ONLY: t_pheno_memory
  USE mo_carbon_memory_class, ONLY: t_carbon_memory
  USE mo_l2a_memory_class, ONLY: t_l2a_memory

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_carbon_tasks
  PUBLIC ::  recalc_per_tile_vars, calculate_current_c_ag_1_and_bg_sums, rescale_carbon_upon_reference_area_change, &
    &        calculate_current_c_ta_state_sum, check_carbon_conservation, yday_carbon_conservation_test

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
    USE mo_carbon_process,    ONLY: calc_Cpools
    USE mo_jsb_time,          ONLY: is_newday, is_newyear, time_steps_per_day, is_time_experiment_start
    USE mo_carbon_constants,  ONLY: molarMassCO2_kg, sec_per_day

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Declare local variables
    TYPE(t_jsb_model), POINTER   :: model
    LOGICAL                      :: lstart
    INTEGER                      :: iblk, ics, ice
    REAL(wp), allocatable        :: precip_total(:)
    REAL(wp), DIMENSION(options%nc) :: MaxLai
    REAL(wp), DIMENSION(options%nc) :: old_c_state_sum, current_fluxes
    CHARACTER(len=*), PARAMETER  :: routine = modname//':update_C_NPP_pot_allocation'

    ! Declare process configuration and memory Pointers
    TYPE(t_pheno_config), POINTER :: pheno__conf
    TYPE(t_pheno_memory), POINTER :: pheno__mem
    TYPE(t_assimi_memory), POINTER :: assimi__mem
    TYPE(t_carbon_memory), POINTER :: carbon__mem
    TYPE(t_a2l_memory), POINTER :: a2l__mem

    ! Declare pointers to variables in memory
    ! R: former cbalance_type variables
    REAL(wp), POINTER, DIMENSION(:) ::  pseudo_temp
    REAL(wp), POINTER, DIMENSION(:) ::  n_pseudo_temp
    REAL(wp), POINTER, DIMENSION(:) ::  f_pseudo_temp
    REAL(wp), POINTER, DIMENSION(:) ::  pseudo_temp_yDay
    REAL(wp), POINTER, DIMENSION(:) ::  pseudo_precip
    REAL(wp), POINTER, DIMENSION(:) ::  pseudo_precip_yDay
    REAL(wp), POINTER, DIMENSION(:) ::  n_pseudo_precip
    REAL(wp), POINTER, DIMENSION(:) ::  f_pseudo_precip

    REAL(wp), POINTER, DIMENSION(:) ::  topSoilTemp           ! was not cbalance_type
    REAL(wp), POINTER, DIMENSION(:) ::  alpha                 ! was not cbalance_type. Relative soil moisture (used for
                                                     ! heterotrophic respiration)
    REAL(wp), POINTER, DIMENSION(:) :: cconservation_calcCpools

    REAL(wp), POINTER, DIMENSION(:) ::  cflux_c_greenwood_2_litter           ! Carbon flux from the veget. to the litter pools
                                                     ! [mol(C)/m^2(canopy) s]
    REAL(wp), POINTER, DIMENSION(:) ::  c_green           ! C-pool for leaves, fine roots, vegetative organs and
                                                     ! other green (living) parts of vegetation [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_reserve         ! C-pool for carbohydrate reserve (sugars, starches) that
                                                     ! allows plants to survive bad times[mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_woods           ! C-pool for stems, thick roots and other (dead) structural
                                                     !  material of living plants [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_crop_harvest    ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_litter_green_ag ! C-pool for above ground non-woody litter especially from leaves
                                                     !     [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_litter_green_bg ! C-pool for below ground non-woody litter especially from
                                                     !     fine roots [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_litter_wood_ag  ! C-pool for litter originating from above ground woody parts of the
                                                     !    plants [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_litter_wood_bg  ! C-pool for litter originating from below ground woody parts of the
                                                     !    plants [mol(C)/m^2(canopy)]
    !CBAL REAL(wp), POINTER, DIMENSION(:) ::  c_slow            ! C-pool for below ground organic material (coming from the litter
                                                     !     pools) [mol(C)/m^2(canopy)]
    !SIZE CLASS 1
    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_ag1        ! Yasso above ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_ag1       ! Yasso above ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_ag1     ! Yasso above ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_ag1  ! Yasso above ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_bg1        ! Yasso below ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_bg1       ! Yasso below ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_bg1     ! Yasso below ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_bg1  ! Yasso below ground litter-pool for non-soluble litter
                                                      ! [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_humus_1         ! Yasso below ground litter-pool for slow C compartment
                                                      !  [mol(C)^2/m(canopy)]
    !SIZE CLASS 2
    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_ag2        ! Yasso above ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_ag2       ! Yasso above ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_ag2     ! Yasso above ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_ag2  ! Yasso above ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_bg2        ! Yasso below ground litter-pool for acid soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_bg2       ! Yasso below ground litter-pool for water soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_bg2     ! Yasso below ground litter-pool for ethanol soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_bg2  ! Yasso below ground litter-pool for non-soluble litter
                                                      !  [mol(C)/m^2(canopy)]
    REAL(wp), POINTER, DIMENSION(:) ::  c_humus_2         ! Yasso below ground litter-pool for slow C compartment
                                                      !  [mol(C)^2/m(canopy)]

    REAL(wp), POINTER, DIMENSION(:) ::  gross_assimilation_ca  ! Gross assimilation of all canopy layers together
    REAL(wp), POINTER, DIMENSION(:) ::  npp_pot_rate_ca        ! The instantaneous (potential) NPP rate [mol(C)/m^2(canopy) s]
    !REAL(wp), POINTER, DIMENSION(:) ::  npp_pot_Rate_acc      ! averaged NPP rate over one day [mol(C)/m^2(canopy) s]
    REAL(wp), POINTER, DIMENSION(:) ::  lai_sum                ! used to accumulate LAI over a day. Sum of LAI since Midnight.
    REAL(wp), POINTER, DIMENSION(:) ::  npp_pot_sum            ! used to accumulated NPP-Rate over a day
    REAL(wp), POINTER, DIMENSION(:) ::  gpp_sum                ! used to accumulated GPP-Rate over a day
    REAL(wp), POINTER, DIMENSION(:) ::  topSoilTemp_sum        ! used to accumulated topSoilTemp over a day
    REAL(wp), POINTER, DIMENSION(:) ::  alpha_sum              ! used to accumulated water stress coefficient alpha over a day.
    REAL(wp), POINTER, DIMENSION(:) ::  lai_yDayMean           ! previous days mean
    REAL(wp), POINTER, DIMENSION(:) ::  lai_yyDayMean          ! previous previous days mean
    REAL(wp), POINTER, DIMENSION(:) ::  npp_pot_yDayMean       ! mean value of NPP-Rate yesterday (from npp_pot_sum())
                                                      ! [mol(CO2)/(m^2(canopy) s)]
                                                      ! = cbalance%NPP_pot_sum(...)/time_steps_per_day
    REAL(wp), POINTER, DIMENSION(:) ::  npp_act_yDayMean       ! mean value of actual NPP-Rate yesterday, i.e. after N-limitation
                                                      ! Actual NPP after N-limitation and excess carbon drop.
                                                      ! [mol(CO2)/(m^2(canopy) s)]
    REAL(wp), POINTER, DIMENSION(:) ::  gpp_yDayMean
    REAL(wp), POINTER, DIMENSION(:) ::  topSoilTemp_yDayMean
    REAL(wp), POINTER, DIMENSION(:) ::  alpha_yDayMean

    REAL(wp), POINTER, DIMENSION(:) ::  soil_respiration       ! mean daily rate of heterotrophic (soil) respiration
                                                      ! [mol(CO2)/m^2(ground)].
                                                      ! Without N limitation!
    REAL(wp), POINTER, DIMENSION(:) ::  npp_flux_correction    ! Daily updated flux correction from yesterdays carbon balance
                                                      ! [mol(CO2)/m^2(canopy) s]
                                                      ! Amount by which the NPP rate entering the routine has to be corrected.
                                                      ! This correction arises either because otherwise the reserve pool would
                                                      ! get negative (positive correction), or the wood pool would exceed its
                                                      ! maximum value (negative correction).
    REAL(wp), POINTER, DIMENSION(:) ::  excess_NPP             ! That part of NPP that because of structural limits could not be
                                                      ! allocated in carbon pools [mol(CO2)/m^2(canopy) s]

    REAL(wp), POINTER, DIMENSION(:) ::  root_exudates          ! Total root exudates entering to the litter green pools
                                                      ! [mol(C)/m^2(canopy) s]

    REAL(wp), POINTER, DIMENSION(:) ::  cflux_c_green_2_herb     !
    REAL(wp), POINTER, DIMENSION(:) ::  cflux_herb_2_littergreen !
    REAL(wp), POINTER, DIMENSION(:) ::  cflux_herb_2_atm  !
    REAL(wp), POINTER, DIMENSION(:) ::  cflux_npp_2_atm_yday_ta  ! day CO2 flux from actual NPP, required for cconservation test
    REAL(wp), POINTER, DIMENSION(:) ::  cflux_npp_2_atm_ta       ! grid cell averages of net CO2 fluxes between
    REAL(wp), POINTER, DIMENSION(:) ::  cflux_soilresp_2_atm_ta  ! .. biosphere (due to NPP, soil respiration and
    REAL(wp), POINTER, DIMENSION(:) ::  cflux_herb_2_atm_ta      ! .. grazing) and atmosphere [kg(CO2)/(m^2(ground) s)]


    ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
    !REAL(wp), POINTER, DIMENSION(:) :: fract_litter_wood_new
    !
    REAL(wp), POINTER, DIMENSION(:) :: veg_fract_correction
    REAL(wp), POINTER, DIMENSION(:) :: LAI
    REAL(wp), POINTER, DIMENSION(:) :: fract_fpc_max
    REAL(wp), POINTER, DIMENSION(:) :: t_air                       ! Atmosphere temperature (lowest layer) in Kelvin!
    REAL(wp), POINTER, DIMENSION(:) :: rain
    REAL(wp), POINTER, DIMENSION(:) :: snow

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(CARBON_)) RETURN

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    model => Get_model(tile%owner_model_id)
    SELECT TYPE (mm => tile%mem(assimi_)%p); TYPE IS (t_assimi_memory); assimi__mem => mm ; END SELECT
    SELECT TYPE (cc => model%processes(pheno_)%p%config); TYPE IS (t_pheno_config) ; pheno__conf => cc ; END SELECT
    SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT

    ! Set process variables
    pseudo_temp => carbon__mem%pseudo_temp%ptr(ics:ice,iblk)                ! inout
    pseudo_temp_yDay => carbon__mem%pseudo_temp_yDay%ptr(ics:ice,iblk)           ! inout
    n_pseudo_temp => carbon__mem%n_pseudo_temp%ptr(ics:ice,iblk)              ! inout
    f_pseudo_temp => carbon__mem%f_pseudo_temp%ptr(ics:ice,iblk)              ! inout
    pseudo_precip => carbon__mem%pseudo_precip%ptr(ics:ice,iblk)              ! inout
    pseudo_precip_yDay => carbon__mem%pseudo_precip_yDay%ptr(ics:ice,iblk)         ! inout
    n_pseudo_precip => carbon__mem%n_pseudo_precip%ptr(ics:ice,iblk)            ! inout
    f_pseudo_precip => carbon__mem%f_pseudo_precip%ptr(ics:ice,iblk)            ! inout

    topSoilTemp => carbon__mem%topSoilTemp%ptr(ics:ice,iblk)                ! in
    alpha => carbon__mem%alpha%ptr(ics:ice,iblk)                      ! in
    !
    ! This is on canopy area:
    cconservation_calcCpools => carbon__mem%cconservation_calcCpools%ptr(ics:ice,iblk)  ! out

    cflux_c_greenwood_2_litter => carbon__mem%cflux_c_greenwood_2_litter%ptr(ics:ice,iblk) ! out
    c_green => carbon__mem%c_green%ptr(ics:ice,iblk)               ! inout
    c_reserve => carbon__mem%c_reserve%ptr(ics:ice,iblk)             ! inout
    c_woods => carbon__mem%c_woods%ptr(ics:ice,iblk)               ! inout
    c_crop_harvest => carbon__mem%c_crop_harvest%ptr(ics:ice,iblk)         ! inout
    c_litter_green_ag => carbon__mem%c_litter_green_ag%ptr(ics:ice,iblk)     ! inout
    c_litter_green_bg => carbon__mem%c_litter_green_bg%ptr(ics:ice,iblk)     ! inout
    c_litter_wood_ag => carbon__mem%c_litter_wood_ag%ptr(ics:ice,iblk)     ! inout
    c_litter_wood_bg => carbon__mem%c_litter_wood_bg%ptr(ics:ice,iblk)     ! inout
    !CBAL c_slow => carbon__mem%c_slow%ptr(ics:ice,iblk)     ! inout
    !SIZE CLASS 1
    c_acid_ag1 => carbon__mem%c_acid_ag1%ptr(ics:ice,iblk)            ! inout
    c_water_ag1 => carbon__mem%c_water_ag1%ptr(ics:ice,iblk)           ! inout
    c_ethanol_ag1 => carbon__mem%c_ethanol_ag1%ptr(ics:ice,iblk)         ! inout
    c_nonsoluble_ag1 => carbon__mem%c_nonsoluble_ag1%ptr(ics:ice,iblk)      ! inout
    c_acid_bg1 => carbon__mem%c_acid_bg1%ptr(ics:ice,iblk)            ! inout
    c_water_bg1 => carbon__mem%c_water_bg1%ptr(ics:ice,iblk)           ! inout
    c_ethanol_bg1 => carbon__mem%c_ethanol_bg1%ptr(ics:ice,iblk)         ! inout
    c_nonsoluble_bg1 => carbon__mem%c_nonsoluble_bg1%ptr(ics:ice,iblk)      ! inout
    c_humus_1 => carbon__mem%c_humus_1%ptr(ics:ice,iblk)             ! inout
    !SIZE CLASS 2
    c_acid_ag2 => carbon__mem%c_acid_ag2%ptr(ics:ice,iblk)           ! inout
    c_water_ag2 => carbon__mem%c_water_ag2%ptr(ics:ice,iblk)          ! inout
    c_ethanol_ag2 => carbon__mem%c_ethanol_ag2%ptr(ics:ice,iblk)        ! inout
    c_nonsoluble_ag2 => carbon__mem%c_nonsoluble_ag2%ptr(ics:ice,iblk)     ! inout
    c_acid_bg2 => carbon__mem%c_acid_bg2%ptr(ics:ice,iblk)           ! inout
    c_water_bg2 => carbon__mem%c_water_bg2%ptr(ics:ice,iblk)          ! inout
    c_ethanol_bg2 => carbon__mem%c_ethanol_bg2%ptr(ics:ice,iblk)        ! inout
    c_nonsoluble_bg2 => carbon__mem%c_nonsoluble_bg2%ptr(ics:ice,iblk)     ! inout
    c_humus_2 => carbon__mem%c_humus_2%ptr(ics:ice,iblk)            ! inout
    !
    ! xx
    gross_assimilation_ca => assimi__mem%gross_assimilation_ca%ptr(ics:ice,iblk)     ! in
    npp_pot_rate_ca => assimi__mem%npp_pot_rate_ca%ptr(ics:ice,iblk)           ! inout
    !npp_pot_Rate_acc => assimi__mem%npp_pot_Rate_acc%ptr(ics:ice,iblk)         ! inout
    lai_sum => carbon__mem%lai_sum%ptr(ics:ice,iblk)                   ! inout
    npp_pot_sum => carbon__mem%npp_pot_sum%ptr(ics:ice,iblk)               ! inout
    gpp_sum => carbon__mem%gpp_sum%ptr(ics:ice,iblk)                   ! inout
    topSoilTemp_sum => carbon__mem%topSoilTemp_sum%ptr(ics:ice,iblk)           ! inout
    alpha_sum => carbon__mem%alpha_sum%ptr(ics:ice,iblk)                 ! inout
    lai_yDayMean => carbon__mem%lai_yDayMean%ptr(ics:ice,iblk)              ! inout
    lai_yyDayMean => carbon__mem%lai_yyDayMean%ptr(ics:ice,iblk)             ! inout
    npp_pot_yDayMean => carbon__mem%npp_pot_yDayMean%ptr(ics:ice,iblk)          ! inout
    npp_act_yDayMean => carbon__mem%npp_act_yDayMean%ptr(ics:ice,iblk)          ! inout
    gpp_yDayMean => carbon__mem%gpp_yDayMean%ptr(ics:ice,iblk)              ! inout
    topSoilTemp_yDayMean => carbon__mem%topSoilTemp_yDayMean%ptr(ics:ice,iblk)      ! inout
    alpha_yDayMean => carbon__mem%alpha_yDayMean%ptr(ics:ice,iblk)            ! inout

    soil_respiration => carbon__mem%soil_respiration%ptr(ics:ice,iblk)          ! inout
    npp_flux_correction => carbon__mem%npp_flux_correction%ptr(ics:ice,iblk)        ! inout

    excess_NPP => carbon__mem%excess_NPP%ptr(ics:ice,iblk)                 ! inout

    root_exudates => carbon__mem%root_exudates%ptr(ics:ice,iblk)             ! inout
    cflux_c_green_2_herb => carbon__mem%cflux_c_green_2_herb%ptr(ics:ice,iblk)           ! inout
    cflux_herb_2_littergreen => carbon__mem%cflux_herb_2_littergreen%ptr(ics:ice,iblk)       ! inout
    cflux_herb_2_atm => carbon__mem%cflux_herb_2_atm%ptr(ics:ice,iblk)     ! inout
    cflux_npp_2_atm_yday_ta => carbon__mem%cflux_npp_2_atm_yday_ta%ptr(ics:ice,iblk)    ! out
    cflux_npp_2_atm_ta => carbon__mem%cflux_npp_2_atm_ta%ptr(ics:ice,iblk)         ! out
    cflux_soilresp_2_atm_ta => carbon__mem%cflux_soilresp_2_atm_ta%ptr(ics:ice,iblk)    ! out
    cflux_herb_2_atm_ta => carbon__mem%cflux_herb_2_atm_ta%ptr(ics:ice,iblk)   ! out

    ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
    !fract_litter_wood_new => carbon__mem%fract_litter_wood_new%ptr(ics:ice,iblk)
    !
    veg_fract_correction => pheno__mem%veg_fract_correction%ptr(ics:ice,iblk)      ! in
    LAI => pheno__mem%LAI%ptr(ics:ice,iblk)                       ! in
    fract_fpc_max => pheno__mem%fract_fpc_max%ptr(ics:ice,iblk)             ! in
    t_air => a2l__mem%t_air%ptr(ics:ice,iblk)                      ! in
    rain => a2l__mem%rain%ptr(ics:ice,iblk)                       ! in
    snow => a2l__mem%snow%ptr(ics:ice,iblk)                       ! in

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
      IF (pheno__conf%l_forestRegrowth) THEN
         MaxLai(:) = pheno__mem%maxlai_allom%ptr(ics:ice,iblk)
      ELSE
        ! R: Rename MaxLAI in lctlib to LAI_max in JS3 comes from theLand%Vegetation%LAI_max,
        ! which is for each vegetated tile taken from MaxLAI in the lct library.
        MaxLai(:) = model%lctlib(tile%lcts(1)%lib_id)%MaxLAI
      END IF

      !Calculate current sum before operation for c conservation test
      CALL calculate_current_c_ta_state_sum(tile, options, old_c_state_sum(:))

      ! Now update the C pools
      CALL calc_Cpools( &
        &  LAI_yDayMean,                                      & ! in
        &  LAI_yyDayMean,                                     & ! in
        &  MaxLai,                                            & ! in
        &  NPP_pot_yDayMean,                                  & ! in
        &  topSoilTemp_yDayMean,                              & ! in
        &  alpha_yDayMean,                                    & ! in
        !
        &  model%lctlib(tile%lcts(1)%lib_id)%fract_npp_2_woodPool,        & ! in
        &  model%lctlib(tile%lcts(1)%lib_id)%fract_npp_2_reservePool,     & ! in
        &  model%lctlib(tile%lcts(1)%lib_id)%fract_npp_2_exudates,        & ! in
        !
        &  model%lctlib(tile%lcts(1)%lib_id)%fract_green_2_herbivory,     & ! in
        &  model%lctlib(tile%lcts(1)%lib_id)%tau_c_litter_leaf,           & ! in
        &  model%lctlib(tile%lcts(1)%lib_id)%tau_c_litter_wood,           & ! in
        &  model%lctlib(tile%lcts(1)%lib_id)%tau_c_woods,                 & ! in
        &  model%lctlib(tile%lcts(1)%lib_id)%lai_shed_constant,           & ! in
        &  model%lctlib(tile%lcts(1)%lib_id)%fract_c_litter_green2atmos,  & ! in
        &  model%lctlib(tile%lcts(1)%lib_id)%Max_c_content_woods,         & ! in
        &  model%lctlib(tile%lcts(1)%lib_id)%specificLeafArea_C,          & ! in
        &  model%lctlib(tile%lcts(1)%lib_id)%reserveC2leafC,              & ! in
        &  model%lctlib(tile%lcts(1)%lib_id)%PhenologyType,               & ! in  R: is it crop? If yes then this is 5!
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
        &  LeafLit_coef1    =  model%lctlib(tile%lcts(1)%lib_id)%LeafLit_coef(1), & ! in, optional, necessary with yasso
        &  LeafLit_coef2    =  model%lctlib(tile%lcts(1)%lib_id)%LeafLit_coef(2), & ! in, optional, necessary with yasso
        &  LeafLit_coef3    =  model%lctlib(tile%lcts(1)%lib_id)%LeafLit_coef(3), & ! in, optional, necessary with yasso
        &  LeafLit_coef4    =  model%lctlib(tile%lcts(1)%lib_id)%LeafLit_coef(4), & ! in, optional, necessary with yasso
        &  LeafLit_coef5    =  model%lctlib(tile%lcts(1)%lib_id)%LeafLit_coef(5), & ! in, optional, necessary with yasso
        &  WoodLit_coef1    =  model%lctlib(tile%lcts(1)%lib_id)%WoodLit_coef(1), & ! in, optional, necessary with yasso
        &  WoodLit_coef2    =  model%lctlib(tile%lcts(1)%lib_id)%WoodLit_coef(2), & ! in, optional, necessary with yasso
        &  WoodLit_coef3    =  model%lctlib(tile%lcts(1)%lib_id)%WoodLit_coef(3), & ! in, optional, necessary with yasso
        &  WoodLit_coef4    =  model%lctlib(tile%lcts(1)%lib_id)%WoodLit_coef(4), & ! in, optional, necessary with yasso
        &  WoodLit_coef5    =  model%lctlib(tile%lcts(1)%lib_id)%WoodLit_coef(5), & ! in, optional, necessary with yasso
        &  WoodLitterSize   =  model%lctlib(tile%lcts(1)%lib_id)%WoodLitterSize   & ! in, optional, necessary with yasso
        & )

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

      ! cflux_herb_2_atm and soil_respiration are negative (fluxes away from land)
      ! thus + here instead of - to reduce NPP by these fluxes
      current_fluxes = (NPP_act_yDayMean + cflux_herb_2_atm + soil_respiration) &
        & * sec_per_day * veg_fract_correction * fract_fpc_max
      CALL check_carbon_conservation(tile, options, old_c_state_sum(:), &
        & current_fluxes(:), cconservation_calcCpools(:))

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

    cflux_npp_2_atm_yday_ta = molarMassCO2_kg *                 & ! daily C conservation test cannot deal with 
              (- veg_fract_correction * fract_fpc_max           & ! diurnal cycle of NPP,
                   * (NPP_act_yDayMean))                          ! therefore here additionaly the day mean


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

    TYPE(t_carbon_memory), POINTER :: carbon__mem

    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_C_NPP_pot_allocation'

    INTEGER  :: iblk , ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, carbon__mem%lai_yDayMean,ics,ice,iblk,routine)

    !@todo Those commented out are currently not calculated
    ! same for ta variables
    ! non-yasso-pools
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_green_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_woods_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_reserve_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_crop_harvest_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_litter_green_ag_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_litter_green_bg_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_litter_wood_ag_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_litter_wood_bg_ta,ics,ice,iblk,routine)
    !CBAL CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_slow_ta,ics,ice,iblk,routine)
    ! ag1 and bg1
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_acid_ag1_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_water_ag1_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_ethanol_ag1_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_nonsoluble_ag1_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_acid_bg1_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_water_bg1_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_ethanol_bg1_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_nonsoluble_bg1_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_humus_1_ta,ics,ice,iblk,routine)
    !ag2 and bg2
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_acid_ag2_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_water_ag2_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_ethanol_ag2_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_nonsoluble_ag2_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_acid_bg2_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_water_bg2_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_ethanol_bg2_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_nonsoluble_bg2_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_humus_2_ta,ics,ice,iblk,routine)
    ! sums
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_sum_natural_ta,ics,ice,iblk,routine)
!    CALL weighted_by_fract%Aggregate(tile, carbon__mem%npp_pot_sum_ta,ics,ice,iblk,routine)
!    CALL weighted_by_fract%Aggregate(tile, carbon__mem%gpp_sum_ta,ics,ice,iblk,routine)
    ! fluxes
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%cflux_c_greenwood_2_litter_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%soil_respiration_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%npp_flux_correction_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%excess_npp_ta,ics,ice,iblk,routine)
!    CALL weighted_by_fract%Aggregate(tile, carbon__mem%lai_yyDayMean_ta,ics,ice,iblk,routine)
!    CALL weighted_by_fract%Aggregate(tile, carbon__mem%lai_yDayMean_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%npp_pot_yDayMean_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%npp_act_yDayMean_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%gpp_yDayMean_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%root_exudates_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%cflux_c_green_2_herb_ta,ics,ice,iblk,routine)
    ! CALL weighted_by_fract%Aggregate(tile, carbon__mem%cflux_herb_2_littergreen_ta,ics,ice,iblk,routine)
    ! CALL weighted_by_fract%Aggregate(tile, carbon__mem%cflux_herb_2_atm_ta,ics,ice,iblk,routine)
    ! CALL weighted_by_fract%Aggregate(tile, carbon__mem%fract_litter_wood_new_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%cconservation_calcCpools,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%cflux_npp_2_atm_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%cflux_soilresp_2_atm_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%cflux_herb_2_atm_ta,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, carbon__mem%cflux_npp_2_atm_yday_ta,ics,ice,iblk,routine)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_C_NPP_pot_allocation


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

    TYPE(t_carbon_memory), POINTER :: carbon__mem

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//':rescale_carbon_upon_reference_area_change'
    INTEGER                     :: ics, ice, iblk

    REAL(wp), POINTER, DIMENSION(:) ::  c_green
    REAL(wp), POINTER, DIMENSION(:) ::  c_reserve
    REAL(wp), POINTER, DIMENSION(:) ::  c_woods
    REAL(wp), POINTER, DIMENSION(:) ::  c_crop_harvest

    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_ag1
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_ag1
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_ag1
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_ag1
    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_bg1
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_bg1
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_bg1
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_bg1
    REAL(wp), POINTER, DIMENSION(:) ::  c_humus_1

    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_ag2
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_ag2
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_ag2
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_ag2
    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_bg2
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_bg2
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_bg2
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_bg2
    REAL(wp), POINTER, DIMENSION(:) ::  c_humus_2

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

    SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    c_green => carbon__mem%c_green%ptr(ics:ice,iblk)               ! inout
    c_reserve => carbon__mem%c_reserve%ptr(ics:ice,iblk)             ! inout
    c_woods => carbon__mem%c_woods%ptr(ics:ice,iblk)               ! inout
    c_crop_harvest => carbon__mem%c_crop_harvest%ptr(ics:ice,iblk)         ! inout

    c_acid_ag1 => carbon__mem%c_acid_ag1%ptr(ics:ice,iblk)            ! inout
    c_water_ag1 => carbon__mem%c_water_ag1%ptr(ics:ice,iblk)           ! inout
    c_ethanol_ag1 => carbon__mem%c_ethanol_ag1%ptr(ics:ice,iblk)         ! inout
    c_nonsoluble_ag1 => carbon__mem%c_nonsoluble_ag1%ptr(ics:ice,iblk)      ! inout
    c_acid_bg1 => carbon__mem%c_acid_bg1%ptr(ics:ice,iblk)            ! inout
    c_water_bg1 => carbon__mem%c_water_bg1%ptr(ics:ice,iblk)           ! inout
    c_ethanol_bg1 => carbon__mem%c_ethanol_bg1%ptr(ics:ice,iblk)         ! inout
    c_nonsoluble_bg1 => carbon__mem%c_nonsoluble_bg1%ptr(ics:ice,iblk)      ! inout
    c_humus_1 => carbon__mem%c_humus_1%ptr(ics:ice,iblk)             ! inout

    c_acid_ag2 => carbon__mem%c_acid_ag2%ptr(ics:ice,iblk)           ! inout
    c_water_ag2 => carbon__mem%c_water_ag2%ptr(ics:ice,iblk)          ! inout
    c_ethanol_ag2 => carbon__mem%c_ethanol_ag2%ptr(ics:ice,iblk)        ! inout
    c_nonsoluble_ag2 => carbon__mem%c_nonsoluble_ag2%ptr(ics:ice,iblk)     ! inout
    c_acid_bg2 => carbon__mem%c_acid_bg2%ptr(ics:ice,iblk)           ! inout
    c_water_bg2 => carbon__mem%c_water_bg2%ptr(ics:ice,iblk)          ! inout
    c_ethanol_bg2 => carbon__mem%c_ethanol_bg2%ptr(ics:ice,iblk)        ! inout
    c_nonsoluble_bg2 => carbon__mem%c_nonsoluble_bg2%ptr(ics:ice,iblk)     ! inout
    c_humus_2 => carbon__mem%c_humus_2%ptr(ics:ice,iblk)            ! inout

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

    CALL calculate_current_c_ag_1_and_bg_sums(tile, options)
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
  !! @brief Calculates the per tile states from current canopy states
  !!
  !! @todo discuss how to make more general, e.g. get the info on which state variables to work on from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE recalc_per_tile_vars(tile, options,                  &
      & c_green , c_woods, c_reserve, c_crop_harvest,             &
      & c_acid_ag1, c_water_ag1, c_ethanol_ag1, c_nonsoluble_ag1, &
      & c_acid_bg1, c_water_bg1, c_ethanol_bg1, c_nonsoluble_bg1, &
      & c_acid_ag2, c_water_ag2, c_ethanol_ag2, c_nonsoluble_ag2, &
      & c_acid_bg2, c_water_bg2, c_ethanol_bg2, c_nonsoluble_bg2, &
      & c_humus_1, c_humus_2,                                     &
      & soil_respiration, NPP_flux_correction, NPP_pot_yDayMean, NPP_act_yDayMean,     &
      & GPP_yDayMean, root_exudates, cflux_c_green_2_herb, cflux_herb_2_atm,           &
      & cflux_dist_greenreserve_2_soil, cflux_dist_woods_2_soil, cflux_fire_all_2_atm, &
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
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_c_greenwood_2_litter(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: root_exudates(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_c_green_2_herb(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: cflux_herb_2_atm(:)

    TYPE(t_carbon_memory), POINTER :: carbon__mem
    TYPE(t_pheno_memory), POINTER :: pheno__mem

    ! Local variables
    CHARACTER(len=*), PARAMETER     :: routine = modname//':recalc_per_tile_vars'
    INTEGER                         :: ics, ice, iblk

    REAL(wp), POINTER, DIMENSION(:) ::  c_green_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_reserve_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_woods_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_crop_harvest_ta

    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_ag1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_ag1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_ag1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_ag1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_bg1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_bg1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_bg1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_bg1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_humus_1_ta

    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_ag2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_ag2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_ag2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_ag2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_bg2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_bg2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_bg2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_bg2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_humus_2_ta

    REAL(wp), POINTER, DIMENSION(:) ::  cflux_dist_greenreserve_2_soil_ta
    REAL(wp), POINTER, DIMENSION(:) ::  cflux_dist_woods_2_soil_ta
    REAL(wp), POINTER, DIMENSION(:) ::  cflux_fire_all_2_atm_ta
    REAL(wp), POINTER, DIMENSION(:) ::  soil_respiration_ta
    REAL(wp), POINTER, DIMENSION(:) ::  npp_pot_yDayMean_ta
    REAL(wp), POINTER, DIMENSION(:) ::  npp_act_yDayMean_ta
    REAL(wp), POINTER, DIMENSION(:) ::  npp_flux_correction_ta
    REAL(wp), POINTER, DIMENSION(:) ::  gpp_yDayMean_ta
    REAL(wp), POINTER, DIMENSION(:) ::  cflux_c_greenwood_2_litter_ta
    REAL(wp), POINTER, DIMENSION(:) ::  root_exudates_ta
    REAL(wp), POINTER, DIMENSION(:) ::  cflux_c_green_2_herb_ta
    REAL(wp), POINTER, DIMENSION(:) ::  cflux_herb_2_atm_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_sum_natural_ta

    REAL(wp), POINTER, DIMENSION(:) ::  veg_fract_correction
    REAL(wp), POINTER, DIMENSION(:) ::  fract_fpc_max

    SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    c_green_ta => carbon__mem%c_green_ta%ptr(ics:ice,iblk)               ! opt out
    c_reserve_ta => carbon__mem%c_reserve_ta%ptr(ics:ice,iblk)             ! opt out
    c_woods_ta => carbon__mem%c_woods_ta%ptr(ics:ice,iblk)               ! opt out
    c_crop_harvest_ta => carbon__mem%c_crop_harvest_ta%ptr(ics:ice,iblk)         ! opt out

    c_acid_ag1_ta => carbon__mem%c_acid_ag1_ta%ptr(ics:ice,iblk)            ! opt out
    c_water_ag1_ta => carbon__mem%c_water_ag1_ta%ptr(ics:ice,iblk)           ! opt out
    c_ethanol_ag1_ta => carbon__mem%c_ethanol_ag1_ta%ptr(ics:ice,iblk)         ! opt out
    c_nonsoluble_ag1_ta => carbon__mem%c_nonsoluble_ag1_ta%ptr(ics:ice,iblk)      ! opt out
    c_acid_bg1_ta => carbon__mem%c_acid_bg1_ta%ptr(ics:ice,iblk)            ! opt out
    c_water_bg1_ta => carbon__mem%c_water_bg1_ta%ptr(ics:ice,iblk)           ! opt out
    c_ethanol_bg1_ta => carbon__mem%c_ethanol_bg1_ta%ptr(ics:ice,iblk)         ! opt out
    c_nonsoluble_bg1_ta => carbon__mem%c_nonsoluble_bg1_ta%ptr(ics:ice,iblk)      ! opt out
    c_humus_1_ta => carbon__mem%c_humus_1_ta%ptr(ics:ice,iblk)             ! opt out

    c_acid_ag2_ta => carbon__mem%c_acid_ag2_ta%ptr(ics:ice,iblk)           ! opt out
    c_water_ag2_ta => carbon__mem%c_water_ag2_ta%ptr(ics:ice,iblk)          ! opt out
    c_ethanol_ag2_ta => carbon__mem%c_ethanol_ag2_ta%ptr(ics:ice,iblk)        ! opt out
    c_nonsoluble_ag2_ta => carbon__mem%c_nonsoluble_ag2_ta%ptr(ics:ice,iblk)     ! opt out
    c_acid_bg2_ta => carbon__mem%c_acid_bg2_ta%ptr(ics:ice,iblk)           ! opt out
    c_water_bg2_ta => carbon__mem%c_water_bg2_ta%ptr(ics:ice,iblk)          ! opt out
    c_ethanol_bg2_ta => carbon__mem%c_ethanol_bg2_ta%ptr(ics:ice,iblk)        ! opt out
    c_nonsoluble_bg2_ta => carbon__mem%c_nonsoluble_bg2_ta%ptr(ics:ice,iblk)     ! opt out
    c_humus_2_ta => carbon__mem%c_humus_2_ta%ptr(ics:ice,iblk)            ! opt out

    cflux_dist_greenreserve_2_soil_ta => carbon__mem%cflux_dist_greenreserve_2_soil_ta%ptr(ics:ice,iblk) ! opt out
    cflux_dist_woods_2_soil_ta => carbon__mem%cflux_dist_woods_2_soil_ta%ptr(ics:ice,iblk)   ! opt out
    cflux_fire_all_2_atm_ta => carbon__mem%cflux_fire_all_2_atm_ta%ptr(ics:ice,iblk)      ! opt out
    soil_respiration_ta => carbon__mem%soil_respiration_ta%ptr(ics:ice,iblk)          ! opt out
    npp_pot_yDayMean_ta => carbon__mem%npp_pot_yDayMean_ta%ptr(ics:ice,iblk)          ! opt out
    npp_act_yDayMean_ta => carbon__mem%npp_act_yDayMean_ta%ptr(ics:ice,iblk)          ! opt out
    npp_flux_correction_ta => carbon__mem%npp_flux_correction_ta%ptr(ics:ice,iblk)        ! opt out
    gpp_yDayMean_ta => carbon__mem%gpp_yDayMean_ta%ptr(ics:ice,iblk)              ! opt out
    cflux_c_greenwood_2_litter_ta => carbon__mem%cflux_c_greenwood_2_litter_ta%ptr(ics:ice,iblk)! opt out
    root_exudates_ta => carbon__mem%root_exudates_ta%ptr(ics:ice,iblk)             ! opt out
    cflux_c_green_2_herb_ta => carbon__mem%cflux_c_green_2_herb_ta%ptr(ics:ice,iblk)      ! opt out
    cflux_herb_2_atm_ta => carbon__mem%cflux_herb_2_atm_ta%ptr(ics:ice,iblk)          ! opt out

    c_sum_natural_ta => carbon__mem%c_sum_natural_ta%ptr(ics:ice,iblk)           ! out

    veg_fract_correction => pheno__mem%veg_fract_correction%ptr(ics:ice,iblk)       ! in
    fract_fpc_max => pheno__mem%fract_fpc_max%ptr(ics:ice,iblk)              ! in

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

    CALL calculate_current_c_ta_state_sum(tile, options, C_sum_natural_ta(:))

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
  SUBROUTINE calculate_current_c_ta_state_sum(tile, options, current_c_state_sum)

    !USE mo_util,                ONLY: real2string

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)         :: options
    REAL(wp),                   INTENT(inout)      :: current_c_state_sum(:)

    TYPE(t_carbon_memory), POINTER :: carbon__mem

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//':calculate_current_c_ta_state_sum'
    INTEGER                     :: ics, ice, iblk

   

    REAL(wp), POINTER, DIMENSION(:) ::  c_green_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_reserve_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_woods_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_crop_harvest_ta

    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_ag1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_ag1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_ag1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_ag1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_bg1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_bg1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_bg1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_bg1_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_humus_1_ta

    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_ag2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_ag2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_ag2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_ag2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_bg2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_bg2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_bg2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_bg2_ta
    REAL(wp), POINTER, DIMENSION(:) ::  c_humus_2_ta

    SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    c_green_ta => carbon__mem%c_green_ta%ptr(ics:ice,iblk)               ! in
    c_reserve_ta => carbon__mem%c_reserve_ta%ptr(ics:ice,iblk)             ! in
    c_woods_ta => carbon__mem%c_woods_ta%ptr(ics:ice,iblk)               ! in
    c_crop_harvest_ta => carbon__mem%c_crop_harvest_ta%ptr(ics:ice,iblk)         ! in

    c_acid_ag1_ta => carbon__mem%c_acid_ag1_ta%ptr(ics:ice,iblk)            ! in
    c_water_ag1_ta => carbon__mem%c_water_ag1_ta%ptr(ics:ice,iblk)           ! in
    c_ethanol_ag1_ta => carbon__mem%c_ethanol_ag1_ta%ptr(ics:ice,iblk)         ! in
    c_nonsoluble_ag1_ta => carbon__mem%c_nonsoluble_ag1_ta%ptr(ics:ice,iblk)      ! in
    c_acid_bg1_ta => carbon__mem%c_acid_bg1_ta%ptr(ics:ice,iblk)            ! in
    c_water_bg1_ta => carbon__mem%c_water_bg1_ta%ptr(ics:ice,iblk)           ! in
    c_ethanol_bg1_ta => carbon__mem%c_ethanol_bg1_ta%ptr(ics:ice,iblk)         ! in
    c_nonsoluble_bg1_ta => carbon__mem%c_nonsoluble_bg1_ta%ptr(ics:ice,iblk)      ! in
    c_humus_1_ta => carbon__mem%c_humus_1_ta%ptr(ics:ice,iblk)             ! in

    c_acid_ag2_ta => carbon__mem%c_acid_ag2_ta%ptr(ics:ice,iblk)           ! in
    c_water_ag2_ta => carbon__mem%c_water_ag2_ta%ptr(ics:ice,iblk)          ! in
    c_ethanol_ag2_ta => carbon__mem%c_ethanol_ag2_ta%ptr(ics:ice,iblk)        ! in
    c_nonsoluble_ag2_ta => carbon__mem%c_nonsoluble_ag2_ta%ptr(ics:ice,iblk)     ! in
    c_acid_bg2_ta => carbon__mem%c_acid_bg2_ta%ptr(ics:ice,iblk)           ! in
    c_water_bg2_ta => carbon__mem%c_water_bg2_ta%ptr(ics:ice,iblk)          ! in
    c_ethanol_bg2_ta => carbon__mem%c_ethanol_bg2_ta%ptr(ics:ice,iblk)        ! in
    c_nonsoluble_bg2_ta => carbon__mem%c_nonsoluble_bg2_ta%ptr(ics:ice,iblk)     ! in
    c_humus_2_ta => carbon__mem%c_humus_2_ta%ptr(ics:ice,iblk)            ! in

    !$claw expand parallel update 
    current_c_state_sum(ics:ice) = c_green_ta(:) + c_woods_ta(:) + c_reserve_ta(:) + c_crop_harvest_ta(:) + &
        & c_acid_ag1_ta(:) + c_water_ag1_ta(:) + c_ethanol_ag1_ta(:) + c_nonsoluble_ag1_ta(:) + &
        & c_acid_ag2_ta(:) + c_water_ag2_ta(:) + c_ethanol_ag2_ta(:) + c_nonsoluble_ag2_ta(:) + &
        & c_acid_bg1_ta(:) + c_water_bg1_ta(:) + c_ethanol_bg1_ta(:) + c_nonsoluble_bg1_ta(:) + &
        & c_acid_bg2_ta(:) + c_water_bg2_ta(:) + c_ethanol_bg2_ta(:) + c_nonsoluble_bg2_ta(:) + &
        & c_humus_1_ta(:)  + c_humus_2_ta(:)
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

    TYPE(t_carbon_memory), POINTER :: carbon__mem

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//':calculate_current_c_ag_1_and_bg_sums'
    INTEGER                     :: ics, ice, iblk

    REAL(wp), POINTER, DIMENSION(:) ::  c_bg_sum
    REAL(wp), POINTER, DIMENSION(:) ::  c_ag_sum_1

    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_ag1
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_ag1
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_ag1
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_ag1

    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_bg1
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_bg1
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_bg1
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_bg1
    REAL(wp), POINTER, DIMENSION(:) ::  c_humus_1

    REAL(wp), POINTER, DIMENSION(:) ::  c_acid_bg2
    REAL(wp), POINTER, DIMENSION(:) ::  c_water_bg2
    REAL(wp), POINTER, DIMENSION(:) ::  c_ethanol_bg2
    REAL(wp), POINTER, DIMENSION(:) ::  c_nonsoluble_bg2
    REAL(wp), POINTER, DIMENSION(:) ::  c_humus_2

    SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    c_bg_sum => carbon__mem%c_bg_sum%ptr(ics:ice,iblk)             ! out
    c_ag_sum_1 => carbon__mem%c_ag_sum_1%ptr(ics:ice,iblk)           ! out

    c_acid_ag1 => carbon__mem%c_acid_ag1%ptr(ics:ice,iblk)            ! in
    c_water_ag1 => carbon__mem%c_water_ag1%ptr(ics:ice,iblk)           ! in
    c_ethanol_ag1 => carbon__mem%c_ethanol_ag1%ptr(ics:ice,iblk)         ! in
    c_nonsoluble_ag1 => carbon__mem%c_nonsoluble_ag1%ptr(ics:ice,iblk)      ! in

    c_acid_bg1 => carbon__mem%c_acid_bg1%ptr(ics:ice,iblk)            ! in
    c_water_bg1 => carbon__mem%c_water_bg1%ptr(ics:ice,iblk)           ! in
    c_ethanol_bg1 => carbon__mem%c_ethanol_bg1%ptr(ics:ice,iblk)         ! in
    c_nonsoluble_bg1 => carbon__mem%c_nonsoluble_bg1%ptr(ics:ice,iblk)      ! in
    c_humus_1 => carbon__mem%c_humus_1%ptr(ics:ice,iblk)             ! in

    c_acid_bg2 => carbon__mem%c_acid_bg2%ptr(ics:ice,iblk)            ! in
    c_water_bg2 => carbon__mem%c_water_bg2%ptr(ics:ice,iblk)           ! in
    c_ethanol_bg2 => carbon__mem%c_ethanol_bg2%ptr(ics:ice,iblk)         ! in
    c_nonsoluble_bg2 => carbon__mem%c_nonsoluble_bg2%ptr(ics:ice,iblk)      ! in
    c_humus_2 => carbon__mem%c_humus_2%ptr(ics:ice,iblk)             ! in

    !TODO-JN-COHORTS: could be done in mo_carbon_process
    !$claw expand parallel update
    c_bg_sum(ics:ice) = c_acid_bg1(:) + c_water_bg1(:) + c_ethanol_bg1(:) + c_nonsoluble_bg1(:) + &
             & c_acid_bg2(:) + c_water_bg2(:) + c_ethanol_bg2(:) + c_nonsoluble_bg2(:) + &
             & c_humus_1(:)  + c_humus_2(:)

    c_ag_sum_1(ics:ice) =  c_acid_ag1(:) + c_water_ag1(:) + c_ethanol_ag1(:) + c_nonsoluble_ag1(:) + c_humus_1(:)
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
    REAL(wp),                   INTENT(inout)      :: Cconserve(:)

    !TYPE(t_carbon_memory), POINTER :: carbon__mem

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

    USE mo_carbon_constants,   ONLY: molarMassCO2_kg, sec_per_day
    USE mo_jsb_grid,           ONLY: Get_grid

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(in), TARGET :: tile

    TYPE(t_carbon_memory), POINTER :: carbon__mem
    TYPE(t_l2a_memory), POINTER :: l2a__mem
    REAL(wp), POINTER, DIMENSION(:) ::  carbon_conservation_test, yday_c_state_sum

    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid
    TYPE(t_jsb_task_options)        :: options

    REAL(wp), ALLOCATABLE :: C_flux_yday(:)
    INTEGER               :: ics, ice, iblk, nc

    SELECT TYPE (mm => tile%mem(l2a_)%p); TYPE IS (t_l2a_memory); l2a__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT

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

      carbon_conservation_test => l2a__mem%carbon_conservation_test%ptr(ics:ice,iblk)     ! inout
      yday_c_state_sum => l2a__mem%yday_c_state_sum%ptr(ics:ice,iblk)             ! inout

      ! Do not test as long as carbon_conservation_test still carries the initialisation value
      IF (.NOT. ALL(carbon_conservation_test .EQ. -999.0_wp)) THEN
      
        !CO2 flux due to npp and respiration
        c_flux_yday(ics:ice) = - ((carbon__mem%cflux_npp_2_atm_yday_ta%ptr(ics:ice,iblk)     &
          & + carbon__mem%cflux_soilresp_2_atm_ta%ptr(ics:ice,iblk)                   &
          & + carbon__mem%cflux_herb_2_atm_ta%ptr(ics:ice,iblk)) / molarMassco2_kg) * sec_per_day

        IF (model%processes(DISTURB_)%p%config%active) THEN
          C_flux_yday(ics:ice) =  C_flux_yday(ics:ice) &
            & - (carbon__mem%cflux_fire_all_2_atm_ta%ptr(ics:ice,iblk))
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

#endif
END MODULE mo_carbon_interface
