!>
!! @brief Contains the memory class for the carbon process.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-04-6)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_carbon_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: finish
  USE mo_util,             ONLY: One_of
  USE mo_jsb_utils_iface,  ONLY: assign_if_present_allocatable

  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real1d, t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_pool_class,   ONLY: t_jsb_pool
  USE mo_jsb_lct_class,    ONLY: LAND_TYPE, BARE_TYPE, VEG_TYPE
  USE mo_jsb_cqt_class,    ONLY: LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, &
    &                            BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE

  ! Use of processes in this module
  dsl4jsb_Use_processes CARBON_, PHENO_

  ! Use process configurations
  dsl4jsb_Use_config(CARBON_)
  dsl4jsb_Use_config(PHENO_)    ! needed for l_forestRegrowth

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_carbon_memory, max_no_of_vars
  PUBLIC :: t_basic_pool, t_plant_pool, t_abg_pool, t_awenh_pool, t_soil_pool, t_pools, t_my_pools
  PUBLIC :: t_jsb_pool
  PUBLIC :: CARBON, NITROGEN

  INTEGER, PARAMETER :: max_no_of_vars = 300

  TYPE, EXTENDS(t_jsb_pool) :: t_basic_pool
    TYPE(t_jsb_var_real2d), POINTER :: carbon
    TYPE(t_jsb_var_real2d), POINTER :: nitrogen
  CONTAINS
    PROCEDURE :: Init => t_basic_pool_init
    PROCEDURE :: Get_element_name_by_id => t_basic_pool_get_element_name_by_id
  END TYPE t_basic_pool

  ! TYPE, EXTENDS(t_basic_pool) :: t_litter_abg_pool
  !   TYPE(t_basic_pool), POINTER :: ag
  !   TYPE(t_basic_pool), POINTER :: bg
  ! END TYPE t_litter_abg_pool

  ! TYPE, EXTENDS(t_litter_abg_pool) :: t_litter_pool
  !   TYPE(t_litter_abg_pool), POINTER :: green
  !   TYPE(t_litter_abg_pool), POINTER :: wood
  ! END TYPE t_litter_pool

  TYPE, EXTENDS(t_basic_pool) :: t_plant_pool
    TYPE(t_basic_pool),  POINTER :: green
      ! Pool for leaves, fine roots, vegetative organs and other green (living) part of vegetation
    TYPE(t_basic_pool),  POINTER :: wood
      ! Pool for stems, thick roots and other (dead) structural material of living plants
    TYPE(t_basic_pool),  POINTER :: reserve
      ! Pool for carbohydrate reserve (sugars, starches) that allows plants to survive bad times (only carbon)
    TYPE(t_basic_pool),  POINTER :: mobile
      ! Buffer pool for plant internal retranslocation of N (only nitrogen)
  END TYPE t_plant_pool

  ! INTERFACE t_plant_pool
  !   MODULE PROCEDURE :: Create_plant_pool
  ! END INTERFACE t_plant_pool

  TYPE, EXTENDS(t_basic_pool) :: t_abg_pool
    TYPE(t_basic_pool),   POINTER :: acid
    TYPE(t_basic_pool),   POINTER :: water
    TYPE(t_basic_pool),   POINTER :: ethanol
    TYPE(t_basic_pool),   POINTER :: nonsoluble
  END TYPE t_abg_pool

  TYPE, EXTENDS(t_basic_pool) :: t_awenh_pool
    TYPE(t_abg_pool),   POINTER :: ag
    TYPE(t_abg_pool),   POINTER :: bg
    TYPE(t_basic_pool), POINTER :: humus
  END TYPE t_awenh_pool

  TYPE, EXTENDS(t_basic_pool) :: t_soil_pool
    TYPE(t_awenh_pool), POINTER :: woody
    TYPE(t_awenh_pool), POINTER :: nonwoody
    TYPE(t_jsb_var_real2d) ::     &
      & soil_respiration
  END TYPE t_soil_pool

  ! INTERFACE t_soil_pool
  !   MODULE PROCEDURE :: Create_soil_pool
  ! END INTERFACE t_soil_pool

  TYPE, EXTENDS(t_basic_pool) :: t_pools
    TYPE(t_plant_pool), POINTER :: plant      ! Plant pools with unit mol/m^2(canopy)
    ! TYPE(t_plant_pool), POINTER :: plant_ta   ! Plant pools with unit mol/m^2(tile area)
    TYPE(t_basic_pool), POINTER :: crop_harvest
    ! TYPE(t_basic_pool), POINTER :: crop_harvest_ta
    TYPE(t_soil_pool),  POINTER :: soil       ! Soil (yasso) pools with unit mol/m^2(canopy)
    ! TYPE(t_soil_pool),  POINTER :: soil_ta    ! Soil (yasso) pools with unit mol/m^2(tile area)
  END TYPE t_pools

  TYPE, EXTENDS(t_basic_pool) :: t_my_pools
    TYPE(t_pools),      POINTER :: state
    TYPE(t_pools),      POINTER :: state_ta
    TYPE(t_plant_pool), POINTER :: falloc
  END TYPE t_my_pools

  ! INTERFACE t_pools
  !   MODULE PROCEDURE :: Create_pools
  ! END INTERFACE t_pools

  ENUM, BIND(C)
    ENUMERATOR :: CARBON=1, NITROGEN
  END ENUM

  INTEGER, ALLOCATABLE, SAVE :: all_elements(:)

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_carbon_memory

    TYPE(t_my_pools), POINTER :: my_pools

    ! 1d reals:
    !REAL(wp) ::            &
    TYPE(t_jsb_var_real2d) ::     &
       N_pseudo_temp,       & ! Constants used for updating the pseudo-15-day mean temperature
       F_pseudo_temp,       & ! Constants used for updating the pseudo-15-day mean temperature
       N_pseudo_precip,     & ! Constants used for updating the pseudo-15-day mean precipitation
       F_pseudo_precip        ! Constants used for updating the pseudo-15-day mean precipitation

    ! cbalance_type variables
    TYPE(t_jsb_var_real2d) ::     &
      pseudo_temp,                & ! a proxy for 15-day mean air temperature, [K] needed with yasso
      pseudo_temp_yDay,           & !
      pseudo_precip,              & ! a proxy for 15-day mean preciptation rate [Kg/m^2 s] needed with yasso
      pseudo_precip_yDay            !

    TYPE(t_jsb_var_real2d) ::     &
      cconservation_calcCpools,   & ! Carbon conservation test: Deviation from conservation [mol(C)/m^2(grid box)]
      cflux_c_greenwood_2_litter, & ! R: From JS3: cbalance_diag_type. Carbon flux from the vegetation to the litter pools
                                    !    [mol(C)/m^2(canopy) s]
      c_green,                & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      c_woods,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                    ! .. plants [mol(C)/m^2(canopy)]
      c_reserve,              & ! C-pool for carbohydrate reserve (sugars, starches) that allows plants to survive
                                    ! .. bad times[mol(C)/m^2(canopy)]
!@todo why is the c_crop_harvest realised as a density and not as a per tile variable in the first place?
! ++ only required for crop pfts, or?
      c_crop_harvest,         & ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]
     !SIZE CLASS 1
      c_acid_ag1,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      c_water_ag1,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      c_ethanol_ag1,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      c_nonsoluble_ag1,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      c_acid_bg1,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      c_water_bg1,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      c_ethanol_bg1,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      c_nonsoluble_bg1,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      c_humus_1,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)/m^2(canopy)]
     !SIZE CLASS 2
      c_acid_ag2,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      c_water_ag2,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      c_ethanol_ag2,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      c_nonsoluble_ag2,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      c_acid_bg2,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      c_water_bg2,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      c_ethanol_bg2,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      c_nonsoluble_bg2,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      c_humus_2,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)/m^2(canopy)]

      c_bg_sum,              & ! Sum of belowground Yasso pools  [mol(C)/m^2(canopy)]
      c_ag_sum_1,            & ! Sum of aboveground green Yasso pools  [mol(C)/m^2(canopy)]
      max_green_bio,         & ! Annual maximum of C-pool green for NLCC process [mol(C)/m^2(canopy)]
      sla,                   & ! specific leave area for NLCC process [m2(leaf)/mol(C)]

  !JN: vars for analytically equilibrating YASSO humus pools outside of jsb4 simulations (compare jsb3 svn ref 9652)
      c_decomp_humus_1_sum,      & ! annual cflux sum of humus decomposed in YASSO (leaf) [mol(C)/m^2(canopy)]
      c_decomp_humus_2_sum,      & ! annual cflux sum of humus decomposed in YASSO (wood) [mol(C)/m^2(canopy)]
      c_into_humus_1_sum,        & ! annual cflux sum into humus pool in YASSO (leaf) [mol(C)/m^2(canopy)]
      c_into_humus_2_sum,        & ! annual cflux sum into humus pool in YASSO (wood) [mol(C)/m^2(canopy)]

      !
      LAI_sum,                    & ! used to accumulate LAI over a day
      NPP_pot_sum,                    & ! used to accumulated NPP-Rate over a day
      GPP_sum,                    & ! used to accumulated GPP-Rate over a day
      !runoff_sum,                 & ! used to accumulate runoff over a day.  R: Nicht notwendig in JS4.
      soil_respiration,           & ! mean daily rate of heterotrophic (soil) respiration [mol(CO2)/m^2(canopy)].
                                    ! Without N limitation!
      NPP_flux_correction,        & ! Daily updated flux correction from yesterdays carbon balance
                                    ! [mol(CO2)/m^2(canopy) s]
                                    ! Amount by which the NPP rate entering the routine has to be corrected. This
                                    ! correction arises either because otherwise the reserve pool would get
                                    ! negative (positive correction).
      excess_NPP,                 & ! That part of NPP that because of structural limits could not be allocated in
                                    ! carbon pools [mol(CO2)/m^2(canopy) s]
      LAI_yyDayMean,              & ! previous previous day mean
      LAI_yDayMean,               & ! previous day mean
      NPP_pot_yDayMean,           & ! mean value of NPP-Rate yesterday (from NPP_pot_sum()) [mol(CO2)/(m^2(canopy) s)]
      NPP_act_yDayMean,           & ! mean value of actual NPP-Rate after N-limitation and excess carbon drop.
                                    ! "yDayMean" because the actual NPP depends on/is calculated from results of the
                                    ! previous day.
      GPP_yDayMean,               & !
      root_exudates,              & ! Total root exudates entering to the litter green pools [mol(C)/m^2(canopy) s]
      cflux_c_green_2_herb,       & ! Total carbon flux from C pool green to herbivory.
      cflux_herb_2_littergreen,   & ! ..part of cflux_c_green_2_herb that goes into green litter
      cflux_herb_2_atm,           & ! ..part of cflux_c_green_2_herb that goes into atmosphere
      !
      cflux_dist_green_2_soil,    & ! R: should be later renamed to C_2_GreenLitterPools
      cflux_dist_woods_2_soil,    & ! R: should be later renamed to C_2_WoodLitterPools

      co2flux_fire_all_2_atm,     & ! CO2 flux from fire per tile [kg(CO2)/m^2(canopy) s]
      ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
      !fract_litter_wood_new,       & ! new fraction in above ground wood litter pool

  !JN: variables required when l_forestRegrowth = true, i.e. when maxLAI is derived from biomass via allometric relationships
      current_max_green,          & ! Running maximum of the green pool [mol(C) m-2(canopy)]
                                    ! since last maxLAI calculation
      veg_carbon_at_max_green       ! Total vegetation carbon [mol(C) m-2(canopy)]
                                    ! at the time of the green pool maximum (current_max_green)

   ! Output on tile areas:
   TYPE(t_jsb_var_real2d) ::         &
      cflux_c_greenwood_2_litter_ta, & ! R: From JS3: cbalance_diag_type. Carbon flux from the vegetation to the litter pools
                                       !    [mol(C)/m^2(tile area) s]
      c_green_ta,                & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                       ! .. of vegetation [mol(C)/m^2(tile area)]
      c_woods_ta,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                       ! .. plants [mol(C)/m^2(tile area)]
      c_reserve_ta,              & ! C-pool for carbohydrate reserve (sugars, starches) that allows plants to survive
                                       ! .. bad times[mol(C)/m^2(tile area)]
      c_crop_harvest_ta,         & ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]
      !CBAL c_slow_ta,                 & ! C-pool for below ground organic material (coming from the litter pools)
      !CBAL                              !     [mol(C)/m^2(tile area)]
     !SIZE CLASS 1
      c_acid_ag1_ta,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      c_water_ag1_ta,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      c_ethanol_ag1_ta,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      c_nonsoluble_ag1_ta,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      c_acid_bg1_ta,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      c_water_bg1_ta,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      c_ethanol_bg1_ta,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      c_nonsoluble_bg1_ta,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      c_humus_1_ta,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)/m^2(tile area)]
     !SIZE CLASS 2
      c_acid_ag2_ta,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      c_water_ag2_ta,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      c_ethanol_ag2_ta,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      c_nonsoluble_ag2_ta,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      c_acid_bg2_ta,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      c_water_bg2_ta,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      c_ethanol_bg2_ta,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      c_nonsoluble_bg2_ta,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      c_humus_2_ta,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)/m^2(tile area)]

  !JN: vars for analytically equilibrating YASSO humus pools outside of jsb4 simulations (compare jsb3 svn ref 9652)
      c_decomp_humus_1_sum_ta,  & ! annual cflux sum of humus decomposed in YASSO (leaf) [mol(C)/m^2(tile area)]
      c_decomp_humus_2_sum_ta,  & ! annual cflux sum of humus decomposed in YASSO (wood) [mol(C)/m^2(tile area)]
      c_into_humus_1_sum_ta,    & ! annual cflux sum into humus pool in YASSO (leaf) [mol(C)/m^2(tile area)]
      c_into_humus_2_sum_ta,    & ! annual cflux sum into humus pool in YASSO (wood) [mol(C)/m^2(tile area)]

      !
      ! diagnostic carbon sums
      c_sum_veg_ta,             & ! Sum of all vegetation carbon pools [mol(C)/m^2(tile area)]
      c_sum_litter_ag_ta,       & ! Sum of all above ground litter carbon pools [mol(C)/m^2(tile area)]
      c_sum_litter_bg_ta,       & ! Sum of all below ground litter carbon pools [mol(C)/m^2(tile area)]
      c_sum_humus_ta,           & ! Sum of all soil / humus carbon pools [mol(C)/m^2(tile area)]
      c_sum_natural_ta,         & ! Sum of all natural carbon pools [mol(C)/m^2(tile area)]
      !
      !@todo: NPP_pot_sum_ta, GPP_sum_ta, LAI_yyDayMean_ta and LAI_yDayMean_ta are never calculated
!      NPP_pot_sum_ta,                    & ! used to accumulated NPP-Rate over a day
!      GPP_sum_ta,                    & ! used to accumulated GPP-Rate over a day
      soil_respiration_ta,           & ! mean daily rate of heterotrophic (soil) respiration [mol(CO2)/m^2(tile area)].
                                       ! Without N limitation!
      NPP_flux_correction_ta,        & ! Daily updated flux correction from yesterdays carbon balance
                                    ! [mol(CO2)/m^2(tile area) s]
                                    ! Amount by which the NPP rate entering the routine has to be corrected. This
                                    ! correction arises either because otherwise the reserve pool would get
                                    ! negative (positive correction).
      excess_NPP_ta,                 & ! That part of NPP that because of structural limits could not be allocated in
                                    ! carbon pools [mol(CO2)/m^2(tile area) s]
!      LAI_yyDayMean_ta,              & ! previous previous day mean
!      LAI_yDayMean_ta,               & ! previous day mean
      NPP_pot_yDayMean_ta,           & ! mean value of NPP-Rate yesterday (from NPP_pot_sum()) [mol(CO2)/(m^2(tile area) s)]
      NPP_act_yDayMean_ta,           & ! mean value of actual NPP-Rate yesterday, i.e. after N-limitation.
                                       ! Actual NPP (after N-limitation and) excess carbon drop.
                                       ! [mol(CO2)/(m^2(tile area) s)]
      GPP_yDayMean_ta,               & !
      root_exudates_ta,              & ! Total root exudates entering to the litter green pools [mol(C)/m^2(tile area) s]
      cflux_c_green_2_herb_ta,            & ! Total carbon flux from C pool green to herbivory.
      cflux_dist_green_2_soil_ta,  & ! R: should be later renamed to C_2_GreenLitterPools_ta
      cflux_dist_woods_2_soil_ta,   & ! R: should be later renamed to C_2_WoodLitterPools_ta

      ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
      fract_litter_wood_new_ta,      & ! new fraction in above ground wood litter pool

      ! tile averages of net CO2 fluxes between various compartments and atmosphere [kg(CO2)/m^2 s]  
      ! .. positive for emission to atmosphere
      ! .. negative for absorption by biosphere
      co2flux_npp_2_atm_ta,          & ! fluxes due to NPP
      co2flux_soilresp_2_atm_ta,     & ! fluxes due to soil respiration
      co2flux_herb_2_atm_ta,         & ! fluxes due to grazing
      co2flux_fire_all_2_atm_ta,     & ! Carbon immediately released by fire per tile
      co2flux_npp_2_atm_yday_ta        ! previous day CO2 flux from actual NPP, required for C conservation test


    ! Diagnostic global carbon sums e.g. for monitoring (only available with ICON)
    TYPE(t_jsb_var_real1d) ::        &
      c_sum_veg_gsum,                & ! Global sum of all vegetation carbon pools [PgC]
      c_sum_litter_ag_gsum,          & ! Global sum of all above ground litter carbon pools [PgC]
      c_sum_litter_bg_gsum,          & ! Global sum of all below ground litter carbon pools [PgC]
      c_sum_humus_gsum,              & ! Global sum of all soil / humus carbon pools [PgC]
      c_sum_natural_gsum,            & ! Global sum of all natural carbon pools [PgC]
      NPP_act_yDayMean_gsum,         & ! Global actual NPP rate [PgC/yr]
      GPP_yDayMean_gsum,             & ! Global GPP rate [PgC/yr]
      soil_respiration_gsum


  CONTAINS
    PROCEDURE :: Init => Init_carbon_memory
  END TYPE t_carbon_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_memory_class'

CONTAINS

  FUNCTION t_basic_pool_get_element_name_by_id(this, id) RESULT(name)

    CLASS(t_basic_pool), INTENT(in) :: this
    INTEGER,             INTENT(in) :: id
    CHARACTER(LEN=:), ALLOCATABLE   :: name

    SELECT CASE (id)
    CASE (CARBON)
      name='carbon'
    CASE (NITROGEN)
      name='nitrogen'
    END SELECT

  END FUNCTION t_basic_pool_get_element_name_by_id

  SUBROUTINE t_basic_pool_init(this, name, shortname, l_elements, element_list, element_unit)

    CLASS(t_basic_pool),        INTENT(inout) :: this
    CHARACTER(LEN=*),           INTENT(in)    :: name
    CHARACTER(LEN=*),           INTENT(in)    :: shortname
    LOGICAL,          OPTIONAL, INTENT(in)    :: l_elements
    INTEGER,          OPTIONAL, INTENT(in)    :: element_list(:)
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)    :: element_unit

    LOGICAL :: l_add_elements
    INTEGER, ALLOCATABLE :: element_list_loc(:)
    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_basic_pool_init'

    this%name = name
    this%shortname = shortname

    l_add_elements = .TRUE.
    IF (PRESENT(l_elements)) l_add_elements = l_elements

    IF (.NOT. l_add_elements) RETURN

    IF (.NOT. PRESENT(element_unit)) CALL finish(routine, 'Missing unit for elements')

    ALLOCATE(element_list_loc(SIZE(all_elements)))
    element_list_loc(:) = all_elements(:)
    CALL assign_if_present_allocatable(element_list_loc, element_list)
  
    this%element_unit = element_unit
  
    IF (One_of(CARBON, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon', 'C', CARBON, dim=2)

      i = this%Find_element_position(CARBON)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%carbon => selector
        END SELECT
      END IF
    END IF

    IF (One_of(NITROGEN, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen', 'N', NITROGEN, dim=2)
      i = this%Find_element_position(NITROGEN)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%nitrogen => selector
        END SELECT
      END IF
    END IF

  END SUBROUTINE t_basic_pool_init

  FUNCTION Create_pools(name, shortname, unit, lct_ids, pool_prefix, pool_suffix) RESULT(pools)
    ! FUNCTION Create_pools(name, shortname, unit, pool_prefix, pool_suffix) RESULT(pools)

    CHARACTER(LEN=*),           INTENT(in) :: name
    CHARACTER(LEN=*),           INTENT(in) :: shortname
    CHARACTER(LEN=*),           INTENT(in) :: unit
    INTEGER,                    INTENT(in) :: lct_ids(:)
    ! CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: prefix
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: pool_prefix
    ! CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: suffix
    CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: pool_suffix
    TYPE(t_pools),              POINTER    :: pools

    TYPE(t_basic_pool), POINTER :: basic_pool
    TYPE(t_plant_pool), POINTER :: plant
    TYPE(t_soil_pool),  POINTER :: soil

    CHARACTER(LEN=:), ALLOCATABLE :: prefix_loc, suffix_loc

    prefix_loc = ''
    suffix_loc = ''
    CALL assign_if_present_allocatable(prefix_loc, pool_prefix)
    CALL assign_if_present_allocatable(suffix_loc, pool_suffix)
    IF (prefix_loc /= '') prefix_loc = prefix_loc // '_'
    IF (suffix_loc /= '') suffix_loc = '_' // suffix_loc

    ALLOCATE(pools)
    CALL pools%Init(name, shortname, l_elements=.TRUE., element_unit=unit)

    IF ( (     One_of(VEG_TYPE,  lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 &
      &  )) THEN
      plant    => Create_plant_pool(prefix_loc//'plant'//suffix_loc, 'plant', unit=unit)
      CALL pools%Add_pool(plant)
      pools%plant => plant

      ALLOCATE(basic_pool)
      CALL basic_pool%Init(prefix_loc//'crop_harvest'//suffix_loc, 'crop_harvest', l_elements=.TRUE., element_unit=unit)
      CALL pools%Add_pool(basic_pool)
      pools%crop_harvest => basic_pool
      NULLIFY(basic_pool)
    END IF

    IF ( (     One_of(VEG_TYPE,  lct_ids(:)) > 0 &
      &   .OR. One_of(BARE_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 &
      &  )) THEN
      soil    => Create_soil_pool(name=prefix_loc//'soil'//suffix_loc, shortname='soil', unit=unit)
      CALL pools%Add_pool(soil)
      pools%soil => soil
    END IF

    ! CALL pools%Set_paths()

  END FUNCTION Create_pools

  FUNCTION Create_plant_pool(name, shortname, unit, element_list) RESULT(plant)

    CHARACTER(LEN=*),   INTENT(in) :: name
    CHARACTER(LEN=*),   INTENT(in) :: shortname
    CHARACTER(LEN=*),   INTENT(in) :: unit
    INTEGER, OPTIONAL,  INTENT(in) :: element_list(:)
    TYPE(t_plant_pool), POINTER    :: plant

    TYPE(t_basic_pool), POINTER :: basic_pool

    ALLOCATE(plant)
    CALL plant%Init(name=name, shortname=shortname, l_elements=.TRUE., element_unit=unit, element_list=element_list)

      ALLOCATE(basic_pool)
      CALL basic_pool%Init(name='green pool', shortname='green', &
        &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
    CALL plant%Add_pool(basic_pool)
    plant%green => basic_pool
      NULLIFY(basic_pool)
      ALLOCATE(basic_pool)
      CALL basic_pool%Init(name='wood pool', shortname='wood', &
        &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
    CALL plant%Add_pool(basic_pool)
    plant%wood => basic_pool
      NULLIFY(basic_pool)
      ALLOCATE(basic_pool)
      CALL basic_pool%Init(name='reserve pool', shortname='reserve', &
        &                  l_elements=.TRUE., element_unit=unit, element_list=[CARBON])
    CALL plant%Add_pool(basic_pool)
    plant%reserve => basic_pool
      NULLIFY(basic_pool)
      ALLOCATE(basic_pool)
      CALL basic_pool%Init(name='mobile pool', shortname='mobile', &
        &                  l_elements=.TRUE., element_unit=unit, element_list=[NITROGEN])
    CALL plant%Add_pool(basic_pool)
    plant%mobile => basic_pool
    NULLIFY(basic_pool)
    
    ! CALL plant%Add_pool(litter)
    ! plant%litter => litter

  END FUNCTION Create_plant_pool

  FUNCTION Create_soil_pool(name, shortname, unit, element_list) RESULT(soil)

    CHARACTER(LEN=*),   INTENT(in) :: name
    CHARACTER(LEN=*),   INTENT(in) :: shortname
    CHARACTER(LEN=*),   INTENT(in) :: unit
    INTEGER, OPTIONAL,  INTENT(in) :: element_list(:)
    TYPE(t_soil_pool),  POINTER    :: soil

    TYPE(t_basic_pool), POINTER :: basic_pool
    TYPE(t_abg_pool),   POINTER :: abg_pool
    TYPE(t_awenh_pool), POINTER :: awenh_pool

    ALLOCATE(soil)
    CALL soil%Init(name=name, shortname=shortname, l_elements=.TRUE., element_unit=unit, element_list=element_list)

      ALLOCATE(awenh_pool)
      CALL awenh_pool%Init(name='woody pool', shortname='woody', &
        &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        ALLOCATE(abg_pool)
        CALL abg_pool%Init(name='above ground pool', shortname='ag', &
          &                l_elements=.TRUE., element_unit=unit, element_list=element_list)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='acid pool', shortname='acid', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%acid => basic_pool
          NULLIFY(basic_pool)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='water pool', shortname='water', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%water => basic_pool
          NULLIFY(basic_pool)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='ethanol pool', shortname='ethanol', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%ethanol => basic_pool
          NULLIFY(basic_pool)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='nonsoluble pool', shortname='nonsoluble', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%nonsoluble => basic_pool
          NULLIFY(basic_pool)
      CALL awenh_pool%Add_pool(abg_pool)
      awenh_pool%ag => abg_pool
        NULLIFY(abg_pool)
        ALLOCATE(abg_pool)
        CALL abg_pool%Init(name='below ground pool', shortname='bg', &
          &                l_elements=.TRUE., element_unit=unit, element_list=element_list)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='acid pool', shortname='acid', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%acid => basic_pool
          NULLIFY(basic_pool)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='water pool', shortname='water', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%water => basic_pool
          NULLIFY(basic_pool)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='ethanol pool', shortname='ethanol', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%ethanol => basic_pool
          NULLIFY(basic_pool)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='nonsoluble pool', shortname='nonsoluble', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%nonsoluble => basic_pool
          NULLIFY(basic_pool)
      CALL awenh_pool%Add_pool(abg_pool)
      awenh_pool%bg => abg_pool
        ALLOCATE(basic_pool)
        CALL basic_pool%Init(name='humus pool', shortname='humus', &
          &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
      CALL awenh_pool%Add_pool(basic_pool)
      awenh_pool%humus => basic_pool
        NULLIFY(basic_pool)
    CALL soil%Add_pool(awenh_pool)
    soil%woody => awenh_pool
      NULLIFY(awenh_pool)

      ALLOCATE(awenh_pool)
      CALL awenh_pool%Init(name='nonwoody pool', shortname='nonwoody', &
        &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        ALLOCATE(abg_pool)
        CALL abg_pool%Init(name='above ground pool', shortname='ag', &
          &                l_elements=.TRUE., element_unit=unit, element_list=element_list)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='acid pool', shortname='acid', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%acid => basic_pool
          NULLIFY(basic_pool)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='water pool', shortname='water', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%water => basic_pool
          NULLIFY(basic_pool)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='ethanol pool', shortname='ethanol', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%ethanol => basic_pool
          NULLIFY(basic_pool)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='nonsoluble pool', shortname='nonsoluble', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%nonsoluble => basic_pool
          NULLIFY(basic_pool)
      CALL awenh_pool%Add_pool(abg_pool)
      awenh_pool%ag => abg_pool
        NULLIFY(abg_pool)
        ALLOCATE(abg_pool)
        CALL abg_pool%Init(name='below ground pool', shortname='bg', &
          &                l_elements=.TRUE., element_unit=unit, element_list=element_list)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='acid pool', shortname='acid', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%acid => basic_pool
          NULLIFY(basic_pool)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='water pool', shortname='water', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%water => basic_pool
          NULLIFY(basic_pool)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='ethanol pool', shortname='ethanol', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%ethanol => basic_pool
          NULLIFY(basic_pool)
          ALLOCATE(basic_pool)
          CALL basic_pool%Init(name='nonsoluble pool', shortname='nonsoluble', &
            &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
        CALL abg_pool%Add_pool(basic_pool)
        abg_pool%nonsoluble => basic_pool
          NULLIFY(basic_pool)
      CALL awenh_pool%Add_pool(abg_pool)
      awenh_pool%bg => abg_pool
        NULLIFY(abg_pool)
        ALLOCATE(basic_pool)
        CALL basic_pool%Init(name='humus pool', shortname='humus', &
          &                  l_elements=.TRUE., element_unit=unit, element_list=element_list)
      CALL awenh_pool%Add_pool(basic_pool)
      awenh_pool%humus => basic_pool
        NULLIFY(basic_pool)
    CALL soil%Add_pool(awenh_pool)
    soil%nonwoody => awenh_pool
      NULLIFY(awenh_pool)

  END FUNCTION Create_soil_pool

  SUBROUTINE Init_carbon_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_carbon_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_config(PHENO_)    ! needed for l_forestRegrowth

    TYPE(t_jsb_model),  POINTER :: model
    TYPE(t_jsb_grid),   POINTER :: hgrid             ! Horizontal grid
    TYPE(t_jsb_vgrid),  POINTER :: surface, soil_w   ! Vertical grids
    TYPE(t_pools),      POINTER :: pools
    TYPE(t_plant_pool), POINTER :: plant_pool
    ! TYPE(t_my_pools), POINTER :: my_pools
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_carbon_memory'

    IF (lct_ids(1) > 0) CONTINUE ! avoid compiler warning about dummy argument not being used

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    soil_w  => Get_vgrid('soil_depth_water')

    dsl4jsb_Get_config(CARBON_)
    dsl4jsb_Get_config(PHENO_)

    IF (.NOT. ALLOCATED(all_elements)) THEN
      IF (dsl4jsb_Config(CARBON_)%With_nitrogen) THEN
        ALLOCATE(all_elements(2))
        all_elements(:) = [CARBON, NITROGEN]
      ELSE
        ALLOCATE(all_elements(1))
        all_elements(:) = [CARBON]
      END IF
    END IF

    !==============================================================================================================================
    ! Create and populate basic pool structure with densities per canopy area
    !==============================================================================================================================

    ALLOCATE(t_my_pools :: mem%pools)
    CALL mem%pools%Init(name='pools', shortname='pools', l_elements=.FALSE.)

    SELECT TYPE (selected => mem%pools)
    CLASS IS (t_my_pools)
      mem%my_pools => selected
    END SELECT

    !------------------------------------------------------------------------------------------------------------------------------
    ! First the pools
    !------------------------------------------------------------------------------------------------------------------------------
    !
    pools => Create_pools(name='state_pools', shortname='state', unit='mol m-2(canopy)', &
      & lct_ids=lct_ids, pool_suffix='pool')
    ! pools => Create_pools(name='state pools', shortname='state', unit='mol m-2(canopy)', lct_ids=lct_ids)
    CALL mem%my_pools%Add_pool(pools)
    mem%my_pools%state => pools

    pools => Create_pools(name='state_pools_ta', shortname='state_ta', unit='mol m-2(tile area)', &
      & lct_ids=lct_ids, pool_suffix='ta_pool')
    ! pools => Create_pools(name='state pools (ta)', shortname='state_ta', unit='mol m-2(tile area)', lct_ids=lct_ids)
    CALL mem%my_pools%Add_pool(pools)
    mem%my_pools%state_ta => pools

    NULLIFY(pools)

    ! plant_pool => Create_pools(name='falloc_pools', shortname='falloc', unit='1', &
    !   & lct_ids=lct_ids, pool_suffix='pool', l_soil=.FALSE.)
    ! CALL mem%my_pools%Add_pool(plant_pool)
    ! mem%my_pools%falloc => plant_pool

    IF ( (     One_of(VEG_TYPE,  lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 &
      &  )) THEN
      plant_pool => Create_plant_pool('falloc_pool', 'falloc', unit='1', element_list=[CARBON])
      CALL mem%my_pools%Add_pool(plant_pool)
      mem%my_pools%falloc => plant_pool
    END IF

    NULLIFY(plant_pool)
  
    CALL mem%pools%Set_paths()

    CALL mem%Add_var(mem%my_pools%state, &
      & hgrid, surface, soil_w,                                                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & lrestart=.TRUE., initval_r=0.0_wp )

    CALL mem%Add_var(mem%my_pools%state_ta, &
      & hgrid, surface, soil_w,                                                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & lrestart=.FALSE., initval_r=0.0_wp )

    IF (ASSOCIATED(mem%my_pools%falloc)) THEN
      CALL mem%Add_var(mem%my_pools%falloc, &
        & hgrid, surface, soil_w,                                                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )
    END IF

    !------------------------------------------------------------------------------------------------------------------------------
    ! Additional variables for the soil pool
    !------------------------------------------------------------------------------------------------------------------------------
    !
    CALL mem%Add_var(mem%my_pools%state%soil, 'soil_respiration', mem%my_pools%state%soil%soil_respiration,  &
      & hgrid, surface,                                                        &
      & t_cf('soil_respiration', 'mol(C) m-2(canopy) s-1', 'Soil respiration.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    IF (.NOT. ASSOCIATED(mem%parent)) THEN
      ! CALL mem%my_pools%state   %Print()
      ! CALL mem%my_pools%state_ta%Print()
      CALL mem%my_pools%Print()
    END IF

    !==============================================================================================================================
    ! Add other variables to the memory
    !==============================================================================================================================
    !

    ! R: only if we would want to have these variables in 2d (see mo_carbon_init.f90 for their initialization):
    CALL mem%Add_var('n_pseudo_temp', mem%N_pseudo_temp,                                             &
      & hgrid, surface,                                                                              &
      & t_cf('N_pseudo_temp', 'K', 'Normalization for computing pseudo-15-day mean air temperature'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp )

    CALL mem%Add_var('f_pseudo_temp', mem%F_pseudo_temp,                                             &
      & hgrid, surface,                                                                              &
      & t_cf('F_pseudo_temp', '-', 'Exponential factor used for updating the pseudo-15-day mean temperature'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp )

    CALL mem%Add_var('n_pseudo_precip', mem%N_pseudo_precip,                                         &
      & hgrid, surface,                                                                              &
      & t_cf('N_pseudo_precip', 'Kg m-2 s-1', 'Normalization for computing pseudo-15-day mean precipitation rate'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('f_pseudo_precip', mem%F_pseudo_precip,                                         &
      & hgrid, surface,                                                                              &
      & t_cf('F_pseudo_precip', '-', 'Exponential factor used for updating the pseudo-15-day mean precip'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp )


    CALL mem%Add_var('pseudo_temp', mem%pseudo_temp,                                                 &
      & hgrid, surface,                                                                              &
      & t_cf('pseudo_temp', 'K', 'Pseudo-mean air temperature over previous 15 days'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & output_level=MEDIUM,                                                                         &
      & lrestart=.TRUE., initval_r=285.0_wp ) ! initval taken from JSBACH3, arbitary.

    ! R: In the output it should be the same as pseudo_temp, I do not understand why it was written out in JS3.
    CALL mem%Add_var('pseudo_temp_yDay', mem%pseudo_temp_yDay,                                         &
      & hgrid, surface,                                                                                &
      & t_cf('pseudo_temp_yDay', 'K', 'Yesterdays pseudo-mean air temperature over previous 15 days'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & lrestart=.TRUE., initval_r=285.0_wp ) ! initval taken from JSBACH3, arbitary.

    CALL mem%Add_var('pseudo_precip', mem%pseudo_precip,                                             &
      & hgrid, surface,                                                                              &
      & t_cf('pseudo_precip', 'Kg m-2 s-1', 'Pseudo-mean precipitation rate over previous 15 days'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & output_level=MEDIUM,                                                                         &
      & lrestart=.TRUE., initval_r=3.0E-5_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3, arbitary.

    ! R: In the output it should be the same as pseudo_precip, I do not understand why it was written out in JS3.
    CALL mem%Add_var('pseudo_precip_yDay', mem%pseudo_precip_yDay,                                   &
      & hgrid, surface,                                                                              &
      & t_cf('pseudo_precip_yDay', 'Kg m-2 s-1', 'Yesterdays pseudo-mean precipitation rate over previous 15 days'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & lrestart=.TRUE., initval_r=3.0E-5_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3, arbitary.

    ! Carbon variables on canopy area
    CALL mem%Add_var('c_greenwood2litter', mem%cflux_c_greenwood_2_litter,                             &
      & hgrid, surface,                                                                                &
      & t_cf('c_greenwood2litter', 'mol(CO2) m-2(canopy) s-1', 'Total litter flux entering the soil pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=FULL,                                                                             &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_green', mem%c_green,                                                   &
      & hgrid, surface,                                                                        &
      & t_cf('c_green', 'mol(C) m-2(canopy)', 'C-Pool for green parts of vegetation.'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = LIVE_CARBON_CQ_TYPE,    &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_woods', mem%c_woods,                                                              &
      & hgrid, surface,                                                                                   &
      & t_cf('c_woods', 'mol(C) m-2(canopy)', 'C-Pool for structural material/parts of vegetation.'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                &
      & prefix, suffix,                                                                                   &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = LIVE_CARBON_CQ_TYPE,               &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_reserve', mem%c_reserve,                                                                           &
      & hgrid, surface,                                                                                                    &
      & t_cf('c_reserve', 'mol(C) m-2(canopy)', 'C-Pool for reserve carbohydrates (starches, sugars) of vegetation.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = LIVE_CARBON_CQ_TYPE,            &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_crop_harvest', mem%c_crop_harvest,                                             &
      & hgrid, surface,                                                                                &
      & t_cf('c_crop_harvest', 'mol(C) m-2(canopy)', 'C-Pool for Harvested Material from Crops.'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = PRODUCT_CARBON_CQ_TYPE,         &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Yasso 1
    CALL mem%Add_var('c_acid_ag1', mem%c_acid_ag1,                                                               &
      & hgrid, surface,                                                                                          &
      & t_cf('c_acid_ag1', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_ag1', mem%c_water_ag1,                                                               &
      & hgrid, surface,                                                                                            &
      & t_cf('c_water_ag1', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_ag1', mem%c_ethanol_ag1,                                                               &
      & hgrid, surface,                                                                                                &
      & t_cf('c_ethanol_ag1', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix,                                                                                                &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                           &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_ag1', mem%c_nonsoluble_ag1,                                                        &
      & hgrid, surface,                                                                                               &
      & t_cf('c_nonsoluble_ag1', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                           &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_acid_bg1', mem%c_acid_bg1,                                                               &
      & hgrid, surface,                                                                                          &
      & t_cf('c_acid_bg1', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                     &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_bg1', mem%c_water_bg1,                                                               &
      & hgrid, surface,                                                                                            &
      & t_cf('c_water_bg1', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                       &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_bg1', mem%c_ethanol_bg1,                                                              &
      & hgrid, surface,                                                                                               &
      & t_cf('c_ethanol_bg1', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                          &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_bg1', mem%c_nonsoluble_bg1,                                                     &
      & hgrid, surface,                                                                                            &
      & t_cf('c_nonsoluble_bg1', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (belowground).'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_humus_1', mem%c_humus_1,                                     &
      & hgrid, surface,                                                              &
      & t_cf('c_humus_1', 'mol(C) m-2(canopy)', 'C-Pool for humus in Yasso.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                     &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Yasso 2
    CALL mem%Add_var('c_acid_ag2', mem%c_acid_ag2,                                                               &
      & hgrid, surface,                                                                                          &
      & t_cf('c_acid_ag2', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_ag2', mem%c_water_ag2,                                                               &
      & hgrid, surface,                                                                                            &
      & t_cf('c_water_ag2', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_ag2', mem%c_ethanol_ag2,                                                               &
      & hgrid, surface,                                                                                                &
      & t_cf('c_ethanol_ag2', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix,                                                                                                &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                            &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_ag2', mem%c_nonsoluble_ag2,                                                        &
      & hgrid, surface,                                                                                               &
      & t_cf('c_nonsoluble_ag2', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                           &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_acid_bg2', mem%c_acid_bg2,                                                               &
      & hgrid, surface,                                                                                          &
      & t_cf('c_acid_bg2', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_bg2', mem%c_water_bg2,                                                               &
      & hgrid, surface,                                                                                            &
      & t_cf('c_water_bg2', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_bg2', mem%c_ethanol_bg2,                                                              &
      & hgrid, surface,                                                                                               &
      & t_cf('c_ethanol_bg2', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                           &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_bg2', mem%c_nonsoluble_bg2,                                                       &
      & hgrid, surface,                                                                                              &
      & t_cf('c_nonsoluble_bg2', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
      & prefix, suffix,                                                                                              &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                          &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_humus_2', mem%c_humus_2,                                     &
      & hgrid, surface,                                                              &
      & t_cf('c_humus_2', 'mol(C) m-2(canopy)', 'C-Pool for humus in Yasso.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,     &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_bg_sum', mem%c_bg_sum,                                         &
      & hgrid, surface,                                                                &
      & t_cf('c_bg_sum', 'mol(C) m-2(canopy)', 'C-Pool of below ground C in Yasso.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      & prefix, suffix,                                                                &
      & output_level=FULL,                                                             &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('c_ag_sum_1', mem%c_ag_sum_1,                                             &
      & hgrid, surface,                                                                        &
      & t_cf('c_ag_sum_1', 'mol(C) m-2(canopy)', 'C-Pool of above ground green C in Yasso.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=FULL,                                                                     &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('max_green_bio', mem%max_green_bio,                                       &
      & hgrid, surface,                                                                        &
      & t_cf('max_green_bio', 'mol(C) m-2(canopy)', 'annual maximum content of green C-Pool.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=FULL,                                                                     &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('sla', mem%sla,                                                           &
      & hgrid, surface,                                                                        &
      & t_cf('sla', 'm2(leave)/mol(C)', 'specific leave area.'),                               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=FULL,                                                                     &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('lai_sum', mem%LAI_sum,                                         &
      & hgrid, surface,                                                              &
      & t_cf('LAI_sum', '-', 'Used to accumulate LAI over a day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_pot_sum', mem%NPP_pot_sum,                                         &
      & hgrid, surface,                                                                      &
      & t_cf('NPP_pot_sum', 'mol(C) m-2(canopy) s-1', 'Used to accumulate NPP over a day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
      & prefix, suffix,                                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE.) ! initval taken from JSBACH3

    CALL mem%Add_var('gpp_sum', mem%GPP_sum,                                             &
      & hgrid, surface,                                                                  &
      & t_cf('GPP_sum', 'mol(C) m-2(canopy) s-1', 'Used to accumulate GPP over a day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('soil_respiration', mem%soil_respiration,                 &
      & hgrid, surface,                                                        &
      & t_cf('soil_respiration', 'mol(C) m-2(canopy) s-1', 'Soil respiration.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_flux_correction', mem%NPP_flux_correction,                                 &
      & hgrid, surface,                                                                              &
      & t_cf('NPP_flux_correction', 'mol(C) m-2(canopy) s-1',                                        &
      &      'Amount by which the NPP rate has to be corrected to avoid a negative reserve pool.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('excess_npp', mem%excess_NPP,                                                                     &
      & hgrid, surface,                                                                                                &
      & t_cf('excess_NPP', 'mol(C) m-2 (canopy) s-1', 'NPP that could not be not stored in vegetation carbon pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix,                                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('lai_yydaymean', mem%LAI_yyDayMean,                    &
      & hgrid, surface,                                                     &
      & t_cf('LAI_yyDayMean', '-', 'Previous previous day mean.'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var('lai_ydaymean', mem%LAI_yDayMean,                      &
      & hgrid, surface,                                                     &
      & t_cf('LAI_yDayMean', '-', 'Previous day mean.'),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_pot_ydaymean', mem%NPP_pot_yDayMean,                                    &
      & hgrid, surface,                                                                           &
      & t_cf('NPP_pot_yDayMean', 'mol(C) m-2(canopy) s-1', 'Mean NPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
      & prefix, suffix,                                                                           &
      & lrestart=.TRUE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would 
                                                  ! result in a zero divide in the first time step of an initialized run.

    CALL mem%Add_var('npp_act_ydaymean', mem%NPP_act_yDayMean,                                      &
      & hgrid, surface,                                                                             &
      & t_cf('NPP_act_yDayMean', 'mol(C) m-2(canopy) s-1', 'Actual NPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
      & prefix, suffix,                                                                             &
      & lrestart=.TRUE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would 
                                                  ! result in a zero divide in the first time step of an initialized run.

    CALL mem%Add_var('gpp_ydaymean', mem%GPP_yDayMean,                                        &
      & hgrid, surface,                                                                       &
      & t_cf('GPP_yDayMean', 'mol(C) m-2(canopy) s-1', 'Mean GPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
      & prefix, suffix,                                                                       &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('root_exudates', mem%root_exudates,                   &
      & hgrid, surface,                                                    &
      & t_cf('root_exudates', 'mol(C) m-2(canopy) s-1', 'Root exudates'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('cflux_c_green_2_herb', mem%cflux_c_green_2_herb,                                          &
      & hgrid, surface,                                                                                         &
      & t_cf('cflux_c_green_2_herb', 'mol(C) m-2(canopy) s-1', 'Total carbon flux from C pool green to herbivory.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('cflux_c_green_2_herb_lg', mem%cflux_herb_2_littergreen,                                   &
      & hgrid, surface,                                                                                         &
      & t_cf('cflux_herb_2_littergreen', 'mol(C) m-2(canopy) s-1',                                              &
      &      'Carbon flux from cflux_c_green_2_herb into green litter (complementary to cflux_herb_2_atm).'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('cflux_herb_2_atm', mem%cflux_herb_2_atm,                                                  &
      & hgrid, surface,                                                                                         &
      & t_cf('cflux_herb_2_atm', 'mol(C) m-2(canopy) s-1',                                                      &
      &  'Carbon flux from cflux_c_green_2_herb into atmosphere (complementary to cflux_herb_2_littergreen).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('cflux_dist_green_2_soil', mem%cflux_dist_green_2_soil,                            &
      & hgrid, surface,                                                                                 &
      & t_cf('cflux_dist_green_2_soil', 'mol(C) m-2(canopy) s-1',                                       &
      &      'Amount of C relocated by wind damage to the green litter pools.'),                        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
      & prefix, suffix,                                                                                 &
      & output_level=MEDIUM, l_conserve_quan=.TRUE., cons_quan_type_id = FLUX_C_CQ_TYPE,                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('cflux_dist_woods_2_soil', mem%cflux_dist_woods_2_soil,                            &
      & hgrid, surface,                                                                                 &
      & t_cf('cflux_dist_woods_2_soil', 'mol(C) m-2(canopy) s-1',                                           &
      &      'Amount of C relocated by wind and fire damage to the woody litter pools..'),              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
      & prefix, suffix,                                                                                 &
      & output_level=MEDIUM,  l_conserve_quan=.TRUE., cons_quan_type_id = FLUX_C_CQ_TYPE,               &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2flux_fire_all_2_atm', mem%co2flux_fire_all_2_atm,                                        &
      & hgrid, surface,                                                                                           &
      & t_cf('co2flux_fire_all_2_atm', 'kg(CO2) m-2(canopy) s-1', 'CO2 flux into atmosphere due to fire.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                                                           &
      & output_level=MEDIUM, l_conserve_quan=.TRUE., cons_quan_type_id = FLUX_C_CQ_TYPE,                          &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Carbon variables on tile area
    CALL mem%Add_var('c_greenwood2litter_ta', mem%cflux_c_greenwood_2_litter_ta,                       &
      & hgrid, surface,                                                                                &
      & t_cf('c_greenwood2litter_ta', 'mol(CO2) m-2 s-1', 'Total litter flux entering the soil pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=BASIC,                                                                            &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_green_ta', mem%c_green_ta,                                        &
      & hgrid, surface,                                                                   &
      & t_cf('c_green_ta', 'mol(C) m-2', 'C-Pool for green parts of vegetation.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & loutput = .TRUE., output_level=BASIC,                                                                            &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_woods_ta', mem%c_woods_ta,                                                   &
      & hgrid, surface,                                                                              &
      & t_cf('c_woods_ta', 'mol(C) m-2', 'C-Pool for structural material/parts of vegetation.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & loutput = .TRUE., output_level=BASIC,                                                        &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_reserve_ta', mem%c_reserve_ta,                                                                &
      & hgrid, surface,                                                                                               &
      & t_cf('c_reserve_ta', 'mol(C) m-2', 'C-Pool for reserve carbohydrates (starches, sugars) of vegetation.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & loutput = .TRUE., output_level=BASIC,                                                                            &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_crop_harvest_ta', mem%c_crop_harvest_ta,                                  &
      & hgrid, surface,                                                                           &
      & t_cf('c_crop_harvest_ta', 'mol(C) m-2', 'C-Pool for Harvested Material from Crops.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
      & prefix, suffix,                                                                           &
      & output_level=BASIC,                                                                            &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Yasso 1
    CALL mem%Add_var('c_acid_ag1_ta', mem%c_acid_ag1_ta,                                                    &
      & hgrid, surface,                                                                                     &
      & t_cf('c_acid_ag1_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_ag1_ta', mem%c_water_ag1_ta,                                                    &
      & hgrid, surface,                                                                                       &
      & t_cf('c_water_ag1_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
      & prefix, suffix,                                                                                       &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_ag1_ta', mem%c_ethanol_ag1_ta,                                                    &
      & hgrid, surface,                                                                                           &
      & t_cf('c_ethanol_ag1_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                                                           &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_ag1_ta', mem%c_nonsoluble_ag1_ta,                                             &
      & hgrid, surface,                                                                                          &
      & t_cf('c_nonsoluble_ag1_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_acid_bg1_ta', mem%c_acid_bg1_ta,                                                    &
      & hgrid, surface,                                                                                     &
      & t_cf('c_acid_bg1_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_bg1_ta', mem%c_water_bg1_ta,                                                    &
      & hgrid, surface,                                                                                       &
      & t_cf('c_water_bg1_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
      & prefix, suffix,                                                                                       &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_bg1_ta', mem%c_ethanol_bg1_ta,                                                   &
      & hgrid, surface,                                                                                          &
      & t_cf('c_ethanol_bg1_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_bg1_ta', mem%c_nonsoluble_bg1_ta,                                            &
      & hgrid, surface,                                                                                         &
      & t_cf('c_nonsoluble_bg1_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_humus_1_ta', mem%c_humus_1_ta,                          &
      & hgrid, surface,                                                         &
      & t_cf('c_humus_1_ta', 'mol(C) m-2', 'C-Pool for humus in Yasso.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),      &
      & prefix, suffix,                                                         &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Yasso 2
    CALL mem%Add_var('c_acid_ag2_ta', mem%c_acid_ag2_ta,                                                    &
      & hgrid, surface,                                                                                     &
      & t_cf('c_acid_ag2_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_ag2_ta', mem%c_water_ag2_ta,                                                    &
      & hgrid, surface,                                                                                       &
      & t_cf('c_water_ag2_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
      & prefix, suffix,                                                                                       &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_ag2_ta', mem%c_ethanol_ag2_ta,                                                    &
      & hgrid, surface,                                                                                           &
      & t_cf('c_ethanol_ag2_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                                                           &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_ag2_ta', mem%c_nonsoluble_ag2_ta,                                             &
      & hgrid, surface,                                                                                          &
      & t_cf('c_nonsoluble_ag2_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_acid_bg2_ta', mem%c_acid_bg2_ta,                                                    &
      & hgrid, surface,                                                                                     &
      & t_cf('c_acid_bg2_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_bg2_ta', mem%c_water_bg2_ta,                                                    &
      & hgrid, surface,                                                                                       &
      & t_cf('c_water_bg2_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
      & prefix, suffix,                                                                                       &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_bg2_ta', mem%c_ethanol_bg2_ta,                                                   &
      & hgrid, surface,                                                                                          &
      & t_cf('c_ethanol_bg2_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_bg2_ta', mem%c_nonsoluble_bg2_ta,                                            &
      & hgrid, surface,                                                                                         &
      & t_cf('c_nonsoluble_bg2_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_humus_2_ta', mem%c_humus_2_ta,                          &
      & hgrid, surface,                                                         &
      & t_cf('c_humus_2_ta', 'mol(C) m-2', 'C-Pool for humus in Yasso.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),      &
      & prefix, suffix,                                                         &
      & loutput = .TRUE., output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_sum_veg_ta', mem%C_sum_veg_ta,                         &
      & hgrid, surface,                                                        &
      & t_cf('C_sum_veg_ta', 'mol(C) m-2', 'Sum of all vegetation C pools.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & output_level=BASIC,                                                    &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('c_sum_litter_ag_ta', mem%C_sum_litter_ag_ta,                           &
      & hgrid, surface,                                                                      &
      & t_cf('C_sum_litter_ag_ta', 'mol(C) m-2', 'Sum of all above ground litter C pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
      & prefix, suffix,                                                                      &
      & output_level=BASIC,                                                                  &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'c_sum_litter_bg_ta', mem%C_sum_litter_bg_ta,                          &
      & hgrid, surface,                                                                      &
      & t_cf('C_sum_litter_bg_ta', 'mol(C) m-2', 'Sum of all below ground litter C pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
      & prefix, suffix,                                                                      &
      & output_level=BASIC,                                                                  &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('c_sum_humus_ta', mem%C_sum_humus_ta,                     &
      & hgrid, surface,                                                        &
      & t_cf('C_sum_humus_ta', 'mol(C) m-2', 'Sum of all humus C pools.'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & output_level=BASIC,                                                    &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('c_sum_natural_ta', mem%C_sum_natural_ta,                 &
      & hgrid, surface,                                                        &
      & t_cf('C_sum_natural_ta', 'mol(C) m-2', 'Sum of all natural C pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & output_level=BASIC,                                                    &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('soil_respiration_ta', mem%soil_respiration_ta,        &
      & hgrid, surface,                                                     &
      & t_cf('soil_respiration_ta', 'mol(C) m-2 s-1', 'Soil respiration.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & output_level=BASIC,                                                 &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_flux_corr_ta', mem%NPP_flux_correction_ta,                                 &
      & hgrid, surface,                                                                              &
      & t_cf('NPP_flux_correction_ta', 'mol(C) m-2 s-1',                                             &
      &      'Amount by which the NPP rate has to be corrected to avoid a negative reserve pool.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('excess_npp_ta', mem%excess_NPP_ta,                                                          &
      & hgrid, surface,                                                                                           &
      & t_cf('excess_NPP_ta', 'mol(C) m-2  s-1', 'NPP that could not be not stored in vegetation carbon pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                                                           &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    !      CALL mem%Add_var('lai_yydaymean_ta', mem%LAI_yyDayMean_ta,              &
    !        & hgrid, surface,                                                     &
    !        & t_cf('LAI_yyDayMean_ta', '-', 'Previous previous day mean.'),          &
    !        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
    !        & prefix, suffix,                                                     &
    !        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    !     CALL mem%Add_var('lai_ydaymean_ta', mem%LAI_yDayMean_ta,                 &
    !        & hgrid, surface,                                                     &
    !        & t_cf('LAI_yDayMean_ta', '-', 'Previous day mean.'),                    &
    !        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
    !        & prefix, suffix,                                                     &
    !        & output_level=BASIC,                                                                                &
    !        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_pot_ydaymean_ta', mem%NPP_pot_yDayMean_ta,                              &
      & hgrid, surface,                                                                           &
      & t_cf('NPP_pot_yDayMean_ta', 'mol(C) m-2 s-1', 'Mean NPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
      & prefix, suffix,                                                                           &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would result
                                                  ! in a zero divide in the first time step of an initialized run.

    CALL mem%Add_var('npp_act_ydaymean_ta', mem%NPP_act_yDayMean_ta,                                &
      & hgrid, surface,                                                                             &
      & t_cf('NPP_act_yDayMean_ta', 'mol(C) m-2 s-1', 'Actual NPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
      & prefix, suffix,                                                                             &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would result
                                                  ! in a zero divide in the first time step of an initialized run.

    CALL mem%Add_var('gpp_ydaymean_ta', mem%GPP_yDayMean_ta,                             &
      & hgrid, surface,                                                                  &
      & t_cf('GPP_yDayMean_ta', 'mol(C) m-2 s-1', 'Mean GPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3


    CALL mem%Add_var('root_exudates_ta', mem%root_exudates_ta,             &
      & hgrid, surface,                                                    &
      & t_cf('root_exudates_ta', 'mol(C) m-2 s-1', 'Root exudates'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3


    CALL mem%Add_var('cflux_c_green_2_herb_ta', mem%cflux_c_green_2_herb_ta,                                   &
      & hgrid, surface,                                                                                    &
      & t_cf('cflux_c_green_2_herb_ta', 'mol(C) m-2 s-1', 'Total carbon flux from C pool green to herbivory.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                 &
      & prefix, suffix,                                                                                    &
      & output_level=MEDIUM,                                                                               &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c2greenlitter_ta', mem%cflux_dist_green_2_soil_ta,                 &
      & hgrid, surface,                                                                         &
      & t_cf('cflux_dist_green_2_soil_ta', 'mol(C) m-2 s-1', 'Carbon flux to green litter.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      &
      & prefix, suffix,                                                                         &
      & output_level=MEDIUM,                                                                               &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c2woodlitter_ta', mem%cflux_dist_woods_2_soil_ta,                        &
      & hgrid, surface,                                                                        &
      & t_cf('cflux_dist_woods_2_soil_ta', 'mol(C) m-2 s-1', 'Carbon flux to woody litter.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=MEDIUM,loutput=.TRUE.,                                                            &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2_l2a_fire_ta', mem%co2flux_fire_all_2_atm_ta,                                    &
      & hgrid, surface,                                                                                   &
      & t_cf('co2flux_fire_all_2_atm_ta', 'kg(CO2) m-2 s-1', 'CO2 flux from C pools to atmosphere due to fire.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                &
      & prefix, suffix, loutput=.TRUE.,                                                                   &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! R: In JSBACH3 the following variables are not written into the output.
    CALL mem%Add_var('co2_l2a_npp_ta', mem%co2flux_npp_2_atm_ta,                                                     &
      & hgrid, surface,                                                                                              &
      & t_cf('co2flux_npp_2_atm_ta', 'kg(CO2) m-2 s-1', 'Net CO2 fluxes between biosphere and atmosphere due to NPP.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
      & prefix, suffix,                                                                                              &
      & output_level=BASIC,                                                                                &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2_l2a_resp_ta', mem%co2flux_soilresp_2_atm_ta,                                     &
      & hgrid, surface,                                                                                    &
      & t_cf('co2flux_soilresp_2_atm_ta', 'kg(CO2) m-2 s-1',                                               &
      &      'Net CO2 fluxes between biosphere and atmosphere due to soil respiration.'),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                 &
      & prefix, suffix,                                                                                    &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2_l2a_herb_ta', mem%co2flux_herb_2_atm_ta,                                          &
      & hgrid, surface,                                                                                     &
      & t_cf('co2flux_herb_2_atm_ta', 'kg(CO2) m-2 s-1',                                                    &
      &      'Net CO2 fluxes between biosphere and atmosphere due to herbivory.'),                          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2_l2a_npp_yday_ta', mem%co2flux_npp_2_atm_yday_ta,                                           &
      & hgrid, surface,                                                                                              &
      & t_cf('co2flux_npp_2_atm_yday_ta', 'kg(CO2) m-2 s-1', 'Net daily CO2 flux, only needed for conservation test.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
      & prefix, suffix,                                                                                              &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('ccons_calccpools', mem%cconservation_calcCpools,            &
      & hgrid, surface,                                                           &
      & t_cf('cconservation_calcCpools', 'mol(C) m-2(tile area)',                 &
      & 'C conservation test: should be zero'),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
      & prefix, suffix,                                                           &
      & output_level=BASIC,                                                       &
      & loutput=.TRUE.,lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('c_decomp_humus_1_sum', mem%c_decomp_humus_1_sum, hgrid, surface,                                 &
      & t_cf('c_decomp_humus_1_sum', 'mol(C) m-2(canopy)', 'annual sum of carbon decomposed from YASSO humus (leaf)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix, loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('c_decomp_humus_2_sum', mem%c_decomp_humus_2_sum, hgrid, surface,                                 &
      & t_cf('c_decomp_humus_2_sum', 'mol(C) m-2(canopy)', 'annual sum of carbon decomposed from YASSO humus (wood)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix, loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('c_into_humus_1_sum', mem%c_into_humus_1_sum, hgrid, surface,                                     &
      & t_cf('c_into_humus_1_sum', 'mol(C) m-2(canopy)', 'annual sum of carbon flux to humus pool in YASSO (leaf)'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix, loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('c_into_humus_2_sum', mem%c_into_humus_2_sum, hgrid, surface,                                     &
      & t_cf('c_into_humus_2_sum', 'mol(C) m-2(canopy)', 'annual sum of carbon flux to humus pool in YASSO (wood)'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix, loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp)

    CALL mem%Add_var('c_decomp_humus_1_sum_ta', mem%c_decomp_humus_1_sum_ta, hgrid, surface,                      &
      & t_cf('c_decomp_humus_1_sum_ta', 'mol(C) m-2', 'annual sum of carbon decomposed from YASSO humus (leaf)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix, loutput=.FALSE., lrestart=dsl4jsb_Config(CARBON_)%diag_humus_fluxes, initval_r=0.0_wp)
    CALL mem%Add_var('c_decomp_humus_2_sum_ta', mem%c_decomp_humus_2_sum_ta, hgrid, surface,                      &
      & t_cf('c_decomp_humus_2_sum_ta', 'mol(C) m-2', 'annual sum of carbon decomposed from YASSO humus (wood)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix, loutput=.FALSE., lrestart=dsl4jsb_Config(CARBON_)%diag_humus_fluxes, initval_r=0.0_wp)
    CALL mem%Add_var('c_into_humus_1_sum_ta', mem%c_into_humus_1_sum_ta, hgrid, surface,                          &
      & t_cf('c_into_humus_1_sum_ta', 'mol(C) m-2', 'annual sum of carbon flux to humus pool in YASSO (leaf)'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix, loutput=.FALSE., lrestart=dsl4jsb_Config(CARBON_)%diag_humus_fluxes, initval_r=0.0_wp)
    CALL mem%Add_var('c_into_humus_2_sum_ta', mem%c_into_humus_2_sum_ta, hgrid, surface,                          &
      & t_cf('c_into_humus_2_sum_ta', 'mol(C) m-2', 'annual sum of carbon flux to humus pool in YASSO (wood)'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix, loutput=.FALSE., lrestart=dsl4jsb_Config(CARBON_)%diag_humus_fluxes, initval_r=0.0_wp)

    IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
      CALL mem%Add_var('max_green', mem%current_max_green,                   &
        & hgrid, surface,                                                    &
        & t_cf('current_max_green', 'mol(C) m-2(canopy)', 'Running maximum of the green pool'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE.)
      CALL mem%Add_var('max_green_veg_c', mem%veg_carbon_at_max_green,       &
        & hgrid, surface,                                                    &
        & t_cf('veg_carbon_at_max_green', 'mol(C) m-2(canopy)', 'Total vegetation carbon at green pool maximum'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE.)
    ENDIF

#ifdef __ICON__
    ! Diagnostic global carbon sums for experiment monitoring 
    !       (1d stream variables are not supported with echam)
    !
    IF ( TRIM(suffix) == 'box' ) THEN
      CALL mem%Add_var('c_sum_veg_gsum', mem%C_sum_veg_gsum,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('C_sum_veg_gsum', 'PgC', 'Global sum of all vegetation C pools.'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('c_sum_litter_ag_gsum', mem%C_sum_litter_ag_gsum,                         &
        & hgrid, surface,                                                                        &
        & t_cf('C_sum_litter_ag_gsum', 'PgC', 'Global sum of all above ground litter C pools'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('c_sum_litter_bg_gsum', mem%C_sum_litter_bg_gsum,                         &
        & hgrid, surface,                                                                        &
        & t_cf('C_sum_litter_bg_gsum', 'PgC', 'Global sum of all below ground litter C pools.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('c_sum_humus_gsum', mem%C_sum_humus_gsum,                                 &
        & hgrid, surface,                                                                        &
        & t_cf('C_sum_humus_gsum', 'PgC', 'Global sum of all humus C pools.'),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('c_sum_natural_gsum', mem%C_sum_natural_gsum,                             &
        & hgrid, surface,                                                                        &
        & t_cf('C_sum_natural_gsum', 'PgC', 'Global sum of all natural C pools.'),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('npp_act_ydaymean_gsum', mem%NPP_act_yDayMean_gsum,                       &
        & hgrid, surface,                                                                        &
        & t_cf('NPP_act_yDayMean_gsum', 'PgC yr-1', 'Global actual NPP rate of the previous day'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('gpp_ydaymean_gsum', mem%GPP_yDayMean_gsum,                               &
        & hgrid, surface,                                                                        &
        & t_cf('GPP_yDayMean_gsum', 'PgC yr-1', 'Global GPP rate of the previous day'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('soil_respiration_gsum', mem%soil_respiration_gsum,                       &
        & hgrid, surface,                                                                        &
        & t_cf('soil_respiration_gsum', 'PgC yr-1', 'Global soil respiration rate'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )
    END IF
#endif

  END SUBROUTINE Init_carbon_memory

#endif
END MODULE mo_carbon_memory_class
