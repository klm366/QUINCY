!>
!! @brief Contains the memory class for the vegetation processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20160830)
!!
MODULE mo_veg_memory_class

  USE mo_kind,                   ONLY: wp, false
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d         , &           ! *slmjsb4
                                       t_jsb_struct, t_jsb_n_struct                             ! *slmjsb4


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_veg_memory, max_no_of_vars
  PUBLIC :: max_no_of_structs, max_no_of_n_structs            ! *slmjsb4
  PUBLIC :: t_veg_element_a, t_veg_pool_a                     ! *slmtemp
  PUBLIC :: Init_veg_pools_elements_a, &                      ! *slmtemp
            Deallocate_var_veg_pools_elements_a
  PUBLIC :: veg_pool, veg_growth, veg_litterfall, veg_falloc, veg_total_flux      ! *slmtemp
  PUBLIC :: veg_exudation, veg_establishment, veg_reserve_use, veg_total_biomass  ! *slmtemp
  

  INTEGER, PARAMETER :: max_no_of_vars      = 250  !> TYPE t_jsb_var
  INTEGER, PARAMETER :: max_no_of_structs   = 20   !> TYPE t_jsb_struct     *slmjsb4
  INTEGER, PARAMETER :: max_no_of_n_structs = 20   !> TYPE t_jsb_n_struct   *slmjsb4


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for helper elements with ALLOCATABLE variables for elements in veg memory
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_veg_element_a
    REAL(wp), ALLOCATABLE :: carbon(:)      ! dimension of variable in memory is var(ics:ice,iblk)
    REAL(wp), ALLOCATABLE :: nitrogen(:)
    REAL(wp), ALLOCATABLE :: phosphorus(:)
    REAL(wp), ALLOCATABLE :: carbon13(:)
    REAL(wp), ALLOCATABLE :: carbon14(:)
    REAL(wp), ALLOCATABLE :: nitrogen15(:)
  CONTAINS
    PRIVATE
    !! operator overload
    ! +
    PROCEDURE               :: Add_element_a                   ! Add two elements
    PROCEDURE               :: Add_element_a_and_scalar        ! Add scalar to element
    PROCEDURE               :: Add_element_a_and_vector       ! Add REAL 1D vector ( var(ics:ice) ) to element
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_element_a, Add_element_a_and_scalar, Add_element_a_and_vector
    !! operator overload
    ! -
    PROCEDURE               :: Subtract_element_a                   ! Subtract one element from another element
    PROCEDURE               :: Subtract_element_a_and_scalar        ! Subtract scalar from element
    PROCEDURE               :: Subtract_element_a_and_vector        ! Subtract vector ( var(ics:ice) ) from element
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_element_a, Subtract_element_a_and_scalar, Subtract_element_a_and_vector
    !! operator overload
    ! *
    PROCEDURE               :: Multiply_element_a                   ! Multiply two elements
    PROCEDURE               :: Multiply_element_a_and_scalar        ! Multiply element with 
    PROCEDURE, PASS(this2)  :: Multiply_scalar_and_element_a        ! Multiply scalar with element
    PROCEDURE               :: Multiply_element_a_and_vector        ! Multiply element with REAL 1D vector ( var(ics:ice) )
    PROCEDURE, PASS(this2)  :: Multiply_vector_and_element_a        ! Multiply REAL 1D vector ( var(ics:ice) ) with element
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_element_a, Multiply_element_a_and_scalar, Multiply_scalar_and_element_a, &
                                              Multiply_element_a_and_vector, Multiply_vector_and_element_a
  END TYPE t_veg_element_a


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for helper pools with ALLOCATABLE variables for pools in veg memory
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_veg_pool_a
    TYPE(t_veg_element_a)  :: leaf             !< elements of leaf
    TYPE(t_veg_element_a)  :: fine_root        !< elements of fine_root
    TYPE(t_veg_element_a)  :: coarse_root      !< elements of coarse_root
    TYPE(t_veg_element_a)  :: sap_wood         !< elements of sap_wood
    TYPE(t_veg_element_a)  :: heart_wood       !< elements of heart_wood
    TYPE(t_veg_element_a)  :: labile           !< elements of labile
    TYPE(t_veg_element_a)  :: reserve          !< elements of reserve
    TYPE(t_veg_element_a)  :: fruit            !< elements of fruit
    TYPE(t_veg_element_a)  :: seed_bed         !< elements of seed_bed
  CONTAINS
    PRIVATE
    !! operator overload
    ! +
    PROCEDURE               :: Add_pool_a                   ! Add two pools
    PROCEDURE               :: Add_pool_a_and_scalar        ! Add scalar to pool
    PROCEDURE               :: Add_pool_a_and_vector        ! Add REAL 1D vector ( var(ics:ice) ) to pool
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_pool_a, Add_pool_a_and_scalar, Add_pool_a_and_vector
    !! operator overload
    ! -
    PROCEDURE               :: Subtract_pool_a                   ! Subtract one pool from another pool
    PROCEDURE               :: Subtract_pool_a_and_scalar        ! Subtract scalar from pool
    PROCEDURE               :: Subtract_pool_a_and_vector        ! Subtract vector ( var(ics:ice) ) from pool
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_pool_a, Subtract_pool_a_and_scalar, Subtract_pool_a_and_vector
    !! operator overload
    ! *
    PROCEDURE               :: Multiply_pool_a                   ! Multiply two pools
    PROCEDURE               :: Multiply_pool_a_and_scalar        ! Multiply pool with scalar
    PROCEDURE, PASS(this2)  :: Multiply_scalar_and_pool_a        ! Multiply scalar with pool
    PROCEDURE               :: Multiply_pool_a_and_vector        ! Multiply pool with REAL 1D vector ( var(ics:ice) )
    PROCEDURE, PASS(this2)  :: Multiply_vector_and_pool_a        ! Multiply REAL 1D vector ( var(ics:ice) ) with pool
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_pool_a, Multiply_pool_a_and_scalar, Multiply_scalar_and_pool_a, &
                                              Multiply_pool_a_and_vector, Multiply_vector_and_pool_a
  END TYPE t_veg_pool_a



  !-----------------------------------------------------------------------------------------------------
  !> Type definition for elements in veg pool
  !! 
  !! elements in veg would be 2D (while they would be 3D in soil)
  !! 
  !! *slmtemp 1D -> 2D
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_struct) :: t_slm_struct_veg_element
    TYPE(t_jsb_var_real2d)    :: carbon       !< total carbon concentration of pool in mol / m2
    TYPE(t_jsb_var_real2d)    :: nitrogen     !< totol nitrogen concentration of pool in mol / m2
    TYPE(t_jsb_var_real2d)    :: phosphorus   !< phosphorus concentration of pool in mol / m2
    TYPE(t_jsb_var_real2d)    :: carbon13     !< 13carbon concentration of pool in mol / m2
    TYPE(t_jsb_var_real2d)    :: carbon14     !< 14carbon concentration of pool in pmol / m2
    TYPE(t_jsb_var_real2d)    :: nitrogen15   !< 15nitrogen concentration of pool in mol / m2
  CONTAINS
    PROCEDURE :: Init       => Init_t_slm_struct_veg_element
    PROCEDURE :: Add_vars   => Add_element_variables_to_memory
  END TYPE t_slm_struct_veg_element


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for veg pool with nested TYPEs
  !!
  !! used for: pool, and fluxes: growth and litterfall
  !!  
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_n_struct) :: t_slm_struct_veg_pool
    TYPE(t_slm_struct_veg_element)  :: leaf             !< elements of leaf
    TYPE(t_slm_struct_veg_element)  :: fine_root        !< elements of fine_root
    TYPE(t_slm_struct_veg_element)  :: coarse_root      !< elements of coarse_root
    TYPE(t_slm_struct_veg_element)  :: sap_wood         !< elements of sap_wood
    TYPE(t_slm_struct_veg_element)  :: heart_wood       !< elements of heart_wood
    TYPE(t_slm_struct_veg_element)  :: labile           !< elements of labile
    TYPE(t_slm_struct_veg_element)  :: reserve          !< elements of reserve
    TYPE(t_slm_struct_veg_element)  :: fruit            !< elements of fruit
    TYPE(t_slm_struct_veg_element)  :: seed_bed         !< elements of seed_bed
  CONTAINS
    PROCEDURE :: Init           => Init_t_slm_struct_veg_pool
    PROCEDURE :: Add_elements   => Add_pool_elements_to_memory
  END TYPE t_slm_struct_veg_pool


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for veg memory
  !!
  !! @par includes: \n
  !!    veg variables veg\%VAR \n
  !!    pool  \n
  !!    growth flux \n
  !!    litterfall flux \n
  !!    moving averagages of variables VAR_mavg   \n
  !!    daytime averages of variables VAR_daytime \n
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_veg_memory

    !> variables
    !! 
    !! 
    ! basic plant 2D
    TYPE(t_jsb_var_real2d)          :: &   
                                    height            , &            !< plant height [m]
                                    diameter          , &            !< plant diameter [m]
                                    dens_ind          , &            !< density of individuals [# m-2]
                                    lai               , &            !< actual leaf area index [m2 m-2]
                                    target_cn_leaf    , &            !< target foliar C:N ratio [mol C mol-1 N]
                                    target_np_leaf    , &            !< target foliar N:P ratio [mol N mol-1 P]
                                    mean_leaf_age     , &            !< average foliage age in days [days]
                                    cohort_age        , &            !< age of the cohort in years, used with veg_dynamics_scheme cohort [years]
                                    f_n_demand        , &            !< relative plant demand for N
                                    f_p_demand        , &            !< relative plant demand for P
                                    beta_sinklim_ps   , &            !< scaling factor to account for a direct sink limitation contraint on photosynthesis [unitless]
                                    dphi              , &            !< plant-soil water potential gradient [MPa]
                                    growth_req_n      , &            !< moles N required for a unit of C growth under current allocation fractions [mol N m-2 -1]
                                    growth_req_p      , &            !< moles P required for a unit of N growth under current allocation fractions [mol P m-2 -1]
                                    leaf_cn_direction , &            !< indicates direction of change in leaf CN ratio to maximise NPP [unitless]
                                    target_lai        , &            !< LAI implied from sap wood area [m2 m-2]
                                    kstar_labile                     !< labile pool turnover rate assuming only temperature and moisture constraints, no stoichiometry [day-1]

    ! used as logical: 0 = FALSE / 1 = TRUE
    TYPE(t_jsb_var_real2d)          :: &
                                    do_cohort_harvest                !< logical to identify whether the cohort may be harvested, in veg_dynamics_scheme cohort [logical]

    ! basic plant 3D 
    TYPE(t_jsb_var_real3d)          :: &
                                    root_fraction_sl                 !< fractions of root per layer [unitless]  DIMENSION(soil_layer_sb)

    ! target C:N and N:P ratios for various non-leaf tissues (mol/mol)
    TYPE(t_jsb_var_real2d)          :: &
                                    target_cn_fine_root       , &   !< 
                                    target_cn_coarse_root     , &   !< 
                                    target_cn_sap_wood        , &   !< 
                                    target_cn_fruit           , &   !< 
                                    target_np_fine_root       , &   !< 
                                    target_np_coarse_root     , &   !< 
                                    target_np_sap_wood        , &   !< 
                                    target_np_fruit                 !< 

    ! canopy layers
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    leaf_nitrogen_cl  , &           !< N in canopy layer [mmol m-2 LAI]
                                    fn_chl_cl         , &           !< fraction of N in chlorophyll, pigments only [unitless]
                                    fn_et_cl          , &           !< fraction of N in electron transport [unitless]
                                    fn_rub_cl         , &           !< fraction of N in Rubisco [unitless]
                                    fn_pepc_cl        , &           !< fraction of N in PEP carboylase, C4 plants only [unitless]
                                    fn_oth_cl         , &           !< fraction of N not associated with photosynthesis [unitless]
                                    fleaf_sunlit_cl   , &           !< fraction of layer which is sunlit [unitless]
                                    lai_cl            , &           !< leaf area in the canopy layer [m2 m-2]
                                    cumm_lai_cl                     !< cummulative leaf area above the mid-point of the canopy layer [m2 m-2]

    ! fluxes
    TYPE(t_jsb_var_real2d)          :: &   
                                    maint_respiration_pot         , & !< potential maintenance respiration [µmol CO2 m-2 s-1]
                                    maint_respiration             , & !< maintenance respiration [µmol CO2 m-2 s-1]
                                    maint_respiration_c13         , & !< maintenance respiration [µmol 13CO2 m-2 s-1]
                                    maint_respiration_c14         , & !< maintenance respiration [µmol 14CO2 m-2 s-1]
                                    growth_respiration            , & !< growth respiration [µmol CO2 m-2 s-1]
                                    growth_respiration_c13        , & !< growth respiration [µmol 13CO2 m-2 s-1]
                                    growth_respiration_c14        , & !< growth respiration [µmol 14CO2 m-2 s-1]
                                    n_transform_respiration       , & !< respiration associated with N transformation [µmol CO2 m-2 s-1]
                                    n_fixation_respiration        , & !< respiration associated with N fixation [µmol CO2 m-2 s-1]
                                    n_processing_respiration      , & !< respiration associated with N processing [µmol CO2 m-2 s-1]
                                    n_processing_respiration_c13  , & !< respiration associated with N processing [µmol 13CO2 m-2 s-1]
                                    n_processing_respiration_c14  , & !< respiration associated with N processing [µmol 14CO2 m-2 s-1]
                                    npp                           , & !< net primary production [µmol CO2 m-2 s-1]
                                    npp_c13                       , & !< net primary production [µmol 13CO2 m-2 s-1]
                                    npp_c14                       , & !< net primary production [µmol 14CO2 m-2 s-1]
                                    uptake_nh4                    , & !< NH4 uptake from soil [µmol N m-2 s-1]
                                    uptake_nh4_n15                , & !< 15NH4 uptake from soil [µmol N m-2 s-1]
                                    uptake_no3                    , & !< NO3 uptake from soil [µmol N m-2 s-1]
                                    uptake_no3_n15                , & !< 15NO3 uptake from soil [µmol N m-2 s-1]
                                    uptake_po4                    , & !< PO4 uptake from soil [µmol P m-2 s-1]
                                    n_fixation                    , & !< symbiontic N fixation [µmol N m-2 s-1]
                                    n_fixation_n15                , & !< symbiontic N fixation [µmol 15N m-2 s-1]
                                    recycling_leaf_n              , & !< net flux of N from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_leaf_n15            , & !< net flux of 15N from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_leaf_p              , & !< net flux of P from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_fine_root_n         , & !< net flux of N from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_fine_root_n15       , & !< net flux of 15N from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_fine_root_p         , & !< net flux of P from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_heart_wood_n        , & !< net flux of N from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_heart_wood_n15      , & !< net flux of 15N from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_heart_wood_p        , & !< net flux of P from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    net_growth                        !< npp minus turnover [µmol CO2 m-2 s-1]

    ! everything else 2D
    TYPE(t_jsb_var_real2d)          :: &   
                                    fmaint_rate_root          , & !< current N-specific maintentance respiration [~µmol C mol-1 N s-1]
                                    unit_npp                  , & !< incremental change in NPP [µmol CO2 m-2 s-1]
                                    unit_transpiration        , & !< incremental change in transpiration [kg m-2 s-1]
                                    unit_uptake_n_pot         , & !< incremental change in potential N uptake per unit root carbon [µmol N m-2 s-1]
                                    unit_uptake_p_pot         , & !< incremental change in potential P uptake per unit root carbon [µmol P m-2 s-1]
                                    unit_uptake_n_act         , & !< incremental change in actual N uptake per unit root carbon [µmol N m-2 s-1]
                                    unit_uptake_p_act         , & !< incremental change in actual P uptake per unit root carbon [µmol P m-2 s-1]
                                    delta_dens_ind            , & !< change in density of individuals [# m-2 timestep-1]
                                    cost_n_uptake_root        , & !< actual cost for N uptake by roots [mol C mol-1 N]
                                    cost_n_fixation           , & !< actual cost for N uptake through fixation [mol C mol-1 N]
                                    mortality_rate            , & !< background_mort_rate_tree / grass [timestep-1]
                                    k1_opt                    , & !< coefficient for optimal leaf:root ratio
                                    k2_opt                    , & !< coefficient for optimal leaf:root ratio
                                    leaf2sapwood_mass_ratio   , & !< leaf to shoot mass ratio [unitless]
                                    leaf2root_mass_ratio          !< leaf to root mass ratio [unitless]

    !> daytime averages
    !! 
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_daytime_cl          !< fraction of layer which is sunlit average of the previous day [unitless]
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_daytime             , &    !< air temperature average of the previous day [K]
                                    press_srf_daytime         , &    !< air pressure at surface average of the previous day [Pa]
                                    co2_mixing_ratio_daytime  , &    !< co2 mixing ratio of the air average of the previous day [co2 in ppm]
                                    ga_daytime                , &    !< aerodynamic conductance average of the previous day [m s-1]
                                    beta_sinklim_ps_daytime   , &    !< scaling factor to account for soil moisture constraints on photosynthesis average of the previous day [unitless]
                                    t_jmax_opt_daytime               !< optimal temperature for electron transport of photosynthesis average of the previous day [deg C]

    !> accumlated over current daytime
    !! 
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_daytime_dacc_cl     !< fraction of layer which is sunlit value accumulated over current day [unitless]
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_daytime_dacc            , & !< air temperature average of the previous day [K]
                                    press_srf_daytime_dacc        , & !< air pressure at surface average of the previous day [Pa]
                                    co2_mixing_ratio_daytime_dacc , & !< co2 mixing ratio of the air average of the previous day [co2 in ppm]
                                    ga_daytime_dacc               , & !< aerodynamic conductance average of the previous day [m s-1]
                                    beta_sinklim_ps_daytime_dacc  , & !< scaling factor to account for soil moisture constraints on photosynthesis average of the previous day [unitless]
                                    t_jmax_opt_daytime_dacc           !< optimal temperature for electron transport of photosynthesis value accumulated over current day [deg C]
    TYPE(t_jsb_var_real2d)          :: &   
                                    daytime_counter                   !< number of timesteps of current day with light available [#]

    !> moving averages
    !! 
    ! 1.1 memory of conditions for the labile pool dynamics
    TYPE(t_jsb_var_real2d)          :: &   
                                    gpp_tlabile_mavg               , &  !< gpp  [µmol m-2 s-1]
                                    maint_respiration_tlabile_mavg , &  !< maintenance respiration [µmol m-2 s-1]
                                    growth_req_n_tlabile_mavg      , &  !< moles N required for a unit of C growth under current allocation fractions time-averaged [mol N m-2 s-1]
                                    growth_req_p_tlabile_mavg           !< moles P required for a unit of N growth under current allocation fractions time-averaged [mol P m-2 s-1]
    ! 1.2 memory of conditions for phenology calculations
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tphen_mavg                    !< air temperature time-averaged [K]
    ! 1.3 memory for plant uptake (average demand over a couple of days to get rid of diurnal C cycle)
    TYPE(t_jsb_var_real2d)          :: &   
                                    npp_tuptake_mavg               , &  !< npp of the last week [µmol m-2 s-1]
                                    demand_uptake_n_tuptake_mavg   , &  !< scalar of Vmax for root N uptake time-averaged [unitless]
                                    demand_uptake_p_tuptake_mavg   , &  !< scalar of Vmax for root P uptake time-averaged [unitless]
                                    growth_req_n_tuptake_mavg      , &  !< moles N required for a unit of C growth under current allocation fractions time-averaged [mol N m-2 s-1]
                                    growth_req_p_tuptake_mavg      , &  !< moles P required for a unit of N growth under current allocation fractions time-averaged [mol P m-2 s-1]
                                    cost_n_fixation_tuptake_mavg        !< moving average of N fixation cost [mol C mol-1 N]
    ! 1.4 memory of conditions for photosynthesis for the optimisation of foliar N fractions 
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_tfrac_mavg_cl          !< fraction of layer which is sunlit time-averaged [unitless]
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tfrac_mavg               , &  !< air temperature time-averaged [K]
                                    press_srf_tfrac_mavg           , &  !< air pressure at surface time-averaged [Pa]
                                    co2_mixing_ratio_tfrac_mavg    , &  !< co2 mixing ratio of the air time-averaged [co2 in ppm]
                                    ga_tfrac_mavg                  , &  !< aerodynamic conductance time-averaged [m s-1]
                                    beta_sinklim_ps_tfrac_mavg          !< scaling factor to account for soil moisture constraints on photosynthesis time-averaged [unitless]
    ! 1.5 memory of conditions for photosynthesis for the optimisation of foliar N  
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_tcnl_mavg_cl           !< fraction of layer which is sunlit time-averaged [unitless]
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tcnl_mavg                , &  !< air temperature time-averaged [K]
                                    press_srf_tcnl_mavg            , &  !< air pressure at surface time-averaged [Pa]
                                    co2_mixing_ratio_tcnl_mavg     , &  !< co2 mixing ratio of the air time-averaged [co2 in ppm]
                                    ga_tcnl_mavg                   , &  !< aerodynamic conductance time-averaged [m s-1]
                                    beta_sinklim_ps_tcnl_mavg      , &  !< scaling factor to account for soil moisture constraints on conductance time-averaged [unitless]
                                    uptake_n_tcnl_mavg             , &  !< root N uptake time-averaged [µmol N m-2 s-1]
                                    growth_cn_tcnl_mavg            , &  !< new growth CN ratio time-averaged [mol C mol-1 N]
                                    growth_np_tcnl_mavg            , &  !< new growth NP ratio time-averaged [mol N mol-1 P]
                                    npp_tcnl_mavg                  , &  !< NPP time-averaged [µmol CO2 m-2 s-1]
                                    fmaint_rate_tcnl_mavg          , &  !< maintenance respiration rate per unit N for roots
                                    labile_nitrogen_tcnl_mavg      , &  !< labile nitrogen vegetation pool time-averaged [mol N m-2]
                                    labile_carbon_tcnl_mavg        , &  !< labile carbon vegetation pool time-averaged [mol C m-2]
                                    labile_phosphorus_tcnl_mavg         !< labile phosphorus vegetation pool time-averaged [mol P m-2]

    ! 1.6 memory of conditions for photosynthesis for the optimisation of biomass allocation fractions
    TYPE(t_jsb_var_real2d)          :: &   
                                    npp_talloc_mavg                , &  !< NPP [µmol m-2 s-1]
                                    n_fixation_talloc_mavg         , &  !< N fixation time-averaged [µmol N m-2 s-1]
                                    transpiration_talloc_mavg      , &  !< transpiration time-averaged [kg m-2 s-1]
                                    unit_npp_talloc_mavg           , &  !< rate of change in NPP with  increment in leaf mass [µmol C mol-1]
                                    unit_transpiration_talloc_mavg , &  !< incremental change in transpiration time-averaged [kg m-2 s-1]
                                    unit_uptake_n_talloc_mavg      , &  !< actual rate of N uptkae per unit root [µmol N mol-1]
                                    unit_uptake_p_talloc_mavg      , &  !< actual rate of P uptake per unit root [µmol P mol-1]
                                    dphi_talloc_mavg               , &  !< leaf-soil water potential gradient [MPa]
                                    growth_cn_talloc_mavg          , &  !< plant growth cn ratio, equivalent to 1/growth_req_n time-averaged [mol C mol-1 N]
                                    growth_cp_talloc_mavg          , &  !< plant growth cp ratio, equivalent to 1/(growthreq_n*growth_req_p) time-averaged [mol C mol-1 P]
                                    growth_np_talloc_mavg          , &  !< plant growth np ratio, equivalent to 1/(growth_req_p) time-averaged [mol N mol-1 P]
                                    labile_carbon_talloc_mavg      , &  !< labile carbon pool neeeded for allometry calculations time-averaged [mol C m-2]
                                    labile_nitrogen_talloc_mavg    , &  !< labile nitrogen pool neeeded for allometry calculations time-averaged [mol N m-2]
                                    labile_phosphorus_talloc_mavg  , &  !< labile phosphorus pool neeeded for allometry calculations time-averaged [mol P m-2]
                                    reserve_carbon_talloc_mavg     , &  !< reserve carbon needed for optimal allocation time-averaged [mol C m-2]
                                    reserve_nitrogen_talloc_mavg   , &  !< reserve nitrogen needed for optimal allocation time-averaged [mol N m-2]
                                    reserve_phosphorus_talloc_mavg , &  !< reserve phosphorus needed for optimal allocation time-averaged [mol P m-2]
                                    w_root_lim_talloc_mavg              !< root water limitation below which root allocation starts responding time-averaged
                                    
                                    
    ! 1.7 memory over root life time 
    TYPE(t_jsb_var_real2d)          :: &   
                                    leaf2root_troot_mavg           , &  !< leaf to root mass ratio time-averaged [unitless]
                                    unit_npp_troot_mavg            , &  !< rate of change in NPP with increment in leaf mass [µmol C mol-1]
                                    unit_uptake_n_pot_troot_mavg   , &  !< rate of N uptkae per unit root [µmol N mol-1]
                                    unit_uptake_p_pot_troot_mavg   , &  !< rate of P uptake per unit root [µmol P mol-1]
                                    fmaint_rate_troot_mavg              !< moving average of N-specific maintentance respiration [~µmol C mol-1 N s-1]
    ! 1.8 memory over timescale of vegetation dynamics 
    TYPE(t_jsb_var_real2d)          :: &   
                                    an_boc_tvegdyn_mavg            , &  !< net photosynthesis at bottom of canopy, lowest canopy layer existing time-averaged [µmol m-2 s-1]
                                    net_growth_tvegdyn_mavg        , &  !< net growth, npp minus turnover [µmol C m-2 s-1]
                                    lai_tvegdyn_mavg                    !< lai time-averaged [m2 m-2]
    ! 1.9 memory of variables for respiration accclimation 
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tacclim_mavg             , &  !< air temperature for respiration acclimation time-averaged [K]
                                    t_soil_root_tacclim_mavg            !< root temperature for respiration acclimation time-averaged [deg K]
    ! 1.10 memory of certain variables, temperature, ...
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_week_mavg                , &  !< air temperature time-averaged [K]
                                    t_air_month_mavg               , &  !< air temperature time-averaged [K]
                                    t_jmax_opt_mavg                     !< optimal temperature for electron transport of photosynthesis time-averaged [deg C]

    !> structures
    !! 
    !!
    TYPE(t_slm_struct_veg_element)  :: &
                                    exudation       , &   !< flux between fine roots and mycorrhiza (mol m-2 timestep-1)
                                    establishment   , &   !< flux between seed bed and reserve (mol m-2 timestep-1)
                                    reserve_use     , &   !< net flux between reserve and labile pool (mol m-2 timestep-1)
                                    total_biomass         !< per time step sum over veg\%pool\%..

    TYPE(t_slm_struct_veg_pool)     :: &   ! nested pool TYPE t_jsb_n_struct
                                    pool            , &   !< pools of the plant
                                    growth          , &   !< fluxes: growth
                                    litterfall      , &   !< fluxes: litterfall
                                    falloc          , &   !< mo_veg_growth:calc_allocation_fraction -> fractional allocation to pools
                                    total_flux            !< for diagnoses of mass conservation: sum of fluxes


  CONTAINS
    PROCEDURE :: Init => Init_veg_memory
  END TYPE t_veg_memory


  ! definition of veg_ pool and element structures with ALLOCATABLE variables
  ! allocation of variables done with Init_veg_pools_elements_a during model init in init_slm()
  TYPE(t_veg_pool_a)      :: veg_pool           , &
                             veg_growth         , &
                             veg_litterfall     , &
                             veg_falloc         , &
                             veg_total_flux
  TYPE(t_veg_element_a)   :: veg_exudation      , &
                             veg_establishment  , &
                             veg_reserve_use    , &
                             veg_total_biomass


  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class'


CONTAINS

  !-----------------------------------------------------------------------------------------------------!
  ! t_veg_element_a operator overload                                                                   !
  !-----------------------------------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------------------------------
  !> + (t_veg_element_a + t_veg_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_element_a(this,element2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this, element2
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)       ! note: Function results are automatically deallocated after the result is "used" in the outer scope.
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      + element2%carbon(:)
    element3%nitrogen(:)     = this%nitrogen(:)    + element2%nitrogen(:)
    element3%phosphorus(:)   = this%phosphorus(:)  + element2%phosphorus(:)
    element3%carbon13(:)     = this%carbon13(:)    + element2%carbon13(:)
    element3%carbon14(:)     = this%carbon14(:)    + element2%carbon14(:)
    element3%nitrogen15(:)   = this%nitrogen15(:)  + element2%nitrogen15(:)
  END FUNCTION Add_element_a
  !-----------------------------------------------------------------------------------------------------
  !> + (t_veg_element_a + scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_element_a_and_scalar(this,scalar) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this
    REAL(wp),               INTENT(in)  :: scalar
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      + scalar
    element3%nitrogen(:)     = this%nitrogen(:)    + scalar
    element3%phosphorus(:)   = this%phosphorus(:)  + scalar
    element3%carbon13(:)     = this%carbon13(:)    + scalar
    element3%carbon14(:)     = this%carbon14(:)    + scalar
    element3%nitrogen15(:)   = this%nitrogen15(:)  + scalar
  END FUNCTION Add_element_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> + (t_veg_element_a + 1D vector)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_element_a_and_vector(this,vector) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this
    REAL(wp),               INTENT(in)  :: vector(:)
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      + vector(:)
    element3%nitrogen(:)     = this%nitrogen(:)    + vector(:)
    element3%phosphorus(:)   = this%phosphorus(:)  + vector(:)
    element3%carbon13(:)     = this%carbon13(:)    + vector(:)
    element3%carbon14(:)     = this%carbon14(:)    + vector(:)
    element3%nitrogen15(:)   = this%nitrogen15(:)  + vector(:)
  END FUNCTION Add_element_a_and_vector
  !-----------------------------------------------------------------------------------------------------
  !> - (t_veg_element_a - t_veg_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_element_a(this,element2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this, element2
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      - element2%carbon(:)
    element3%nitrogen(:)     = this%nitrogen(:)    - element2%nitrogen(:)
    element3%phosphorus(:)   = this%phosphorus(:)  - element2%phosphorus(:)
    element3%carbon13(:)     = this%carbon13(:)    - element2%carbon13(:)
    element3%carbon14(:)     = this%carbon14(:)    - element2%carbon14(:)
    element3%nitrogen15(:)   = this%nitrogen15(:)  - element2%nitrogen15(:)
  END FUNCTION Subtract_element_a
  !-----------------------------------------------------------------------------------------------------
  !> - (t_veg_element_a - scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_element_a_and_scalar(this,scalar) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this
    REAL(wp),               INTENT(in)  :: scalar
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      - scalar
    element3%nitrogen(:)     = this%nitrogen(:)    - scalar
    element3%phosphorus(:)   = this%phosphorus(:)  - scalar
    element3%carbon13(:)     = this%carbon13(:)    - scalar
    element3%carbon14(:)     = this%carbon14(:)    - scalar
    element3%nitrogen15(:)   = this%nitrogen15(:)  - scalar
  END FUNCTION Subtract_element_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> - (t_veg_element_a - 1D vector)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_element_a_and_vector(this,vector) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this
    REAL(wp),               INTENT(in)  :: vector(:)
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      - vector(:)
    element3%nitrogen(:)     = this%nitrogen(:)    - vector(:)
    element3%phosphorus(:)   = this%phosphorus(:)  - vector(:)
    element3%carbon13(:)     = this%carbon13(:)    - vector(:)
    element3%carbon14(:)     = this%carbon14(:)    - vector(:)
    element3%nitrogen15(:)   = this%nitrogen15(:)  - vector(:)
  END FUNCTION Subtract_element_a_and_vector
  !-----------------------------------------------------------------------------------------------------
  !> * (t_veg_element_a * t_veg_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_element_a(this,element2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this, element2
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      * element2%carbon(:)
    element3%nitrogen(:)     = this%nitrogen(:)    * element2%nitrogen(:)
    element3%phosphorus(:)   = this%phosphorus(:)  * element2%phosphorus(:)
    element3%carbon13(:)     = this%carbon13(:)    * element2%carbon13(:)
    element3%carbon14(:)     = this%carbon14(:)    * element2%carbon14(:)
    element3%nitrogen15(:)   = this%nitrogen15(:)  * element2%nitrogen15(:)
  END FUNCTION Multiply_element_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_veg_element_a * scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_element_a_and_scalar(this,scalar) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this
    REAL(wp),               INTENT(in)  :: scalar
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      * scalar
    element3%nitrogen(:)     = this%nitrogen(:)    * scalar
    element3%phosphorus(:)   = this%phosphorus(:)  * scalar
    element3%carbon13(:)     = this%carbon13(:)    * scalar
    element3%carbon14(:)     = this%carbon14(:)    * scalar
    element3%nitrogen15(:)   = this%nitrogen15(:)  * scalar
  END FUNCTION Multiply_element_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> * (scalar * t_veg_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_scalar_and_element_a(scalar,this2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),               INTENT(in)  :: scalar
    CLASS(t_veg_element_a), INTENT(in)  :: this2
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this2%carbon)
    ALLOCATE(element3%nitrogen,   source=this2%nitrogen)
    ALLOCATE(element3%phosphorus, source=this2%phosphorus)
    ALLOCATE(element3%carbon13,   source=this2%carbon13)
    ALLOCATE(element3%carbon14,   source=this2%carbon14)
    ALLOCATE(element3%nitrogen15, source=this2%nitrogen15)
    element3%carbon(:)       = scalar * this2%carbon(:)
    element3%nitrogen(:)     = scalar * this2%nitrogen(:)
    element3%phosphorus(:)   = scalar * this2%phosphorus(:)
    element3%carbon13(:)     = scalar * this2%carbon13(:)
    element3%carbon14(:)     = scalar * this2%carbon14(:)
    element3%nitrogen15(:)   = scalar * this2%nitrogen15(:)
  END FUNCTION Multiply_scalar_and_element_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_veg_element_a * 1D vector)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_element_a_and_vector(this,vector) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this
    REAL(wp),               INTENT(in)  :: vector(:)
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      * vector(:)
    element3%nitrogen(:)     = this%nitrogen(:)    * vector(:)
    element3%phosphorus(:)   = this%phosphorus(:)  * vector(:)
    element3%carbon13(:)     = this%carbon13(:)    * vector(:)
    element3%carbon14(:)     = this%carbon14(:)    * vector(:)
    element3%nitrogen15(:)   = this%nitrogen15(:)  * vector(:)
  END FUNCTION Multiply_element_a_and_vector
  !-----------------------------------------------------------------------------------------------------
  !> * (vector * t_veg_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_vector_and_element_a(vector,this2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),               INTENT(in)  :: vector(:)
    CLASS(t_veg_element_a), INTENT(in)  :: this2
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this2%carbon)
    ALLOCATE(element3%nitrogen,   source=this2%nitrogen)
    ALLOCATE(element3%phosphorus, source=this2%phosphorus)
    ALLOCATE(element3%carbon13,   source=this2%carbon13)
    ALLOCATE(element3%carbon14,   source=this2%carbon14)
    ALLOCATE(element3%nitrogen15, source=this2%nitrogen15)
    element3%carbon(:)       = vector(:) * this2%carbon(:)
    element3%nitrogen(:)     = vector(:) * this2%nitrogen(:)
    element3%phosphorus(:)   = vector(:) * this2%phosphorus(:)
    element3%carbon13(:)     = vector(:) * this2%carbon13(:)
    element3%carbon14(:)     = vector(:) * this2%carbon14(:)
    element3%nitrogen15(:)   = vector(:) * this2%nitrogen15(:)
  END FUNCTION Multiply_vector_and_element_a
  
  

  !-----------------------------------------------------------------------------------------------------!
  ! t_veg_pool_a operator overload                                                                      !
  !-----------------------------------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------------------------------
  !> + (t_veg_pool_a + t_veg_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_pool_a(this,pool2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this, pool2
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         + pool2%leaf
    pool3%fine_root     = this%fine_root    + pool2%fine_root
    pool3%coarse_root   = this%coarse_root  + pool2%coarse_root
    pool3%sap_wood      = this%sap_wood     + pool2%sap_wood
    pool3%heart_wood    = this%heart_wood   + pool2%heart_wood
    pool3%labile        = this%labile       + pool2%labile
    pool3%reserve       = this%reserve      + pool2%reserve
    pool3%fruit         = this%fruit        + pool2%fruit
    pool3%seed_bed      = this%seed_bed     + pool2%seed_bed
  END FUNCTION Add_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> + (t_veg_pool_a + scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_pool_a_and_scalar(this,scalar) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this
    REAL(wp),            INTENT(in) :: scalar
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         + scalar
    pool3%fine_root     = this%fine_root    + scalar
    pool3%coarse_root   = this%coarse_root  + scalar
    pool3%sap_wood      = this%sap_wood     + scalar
    pool3%heart_wood    = this%heart_wood   + scalar
    pool3%labile        = this%labile       + scalar
    pool3%reserve       = this%reserve      + scalar
    pool3%fruit         = this%fruit        + scalar
    pool3%seed_bed      = this%seed_bed     + scalar
  END FUNCTION Add_pool_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> + (t_veg_pool_a + 1D vector)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_pool_a_and_vector(this,vector) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this
    REAL(wp),            INTENT(in) :: vector(:)
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         + vector(:)
    pool3%fine_root     = this%fine_root    + vector(:)
    pool3%coarse_root   = this%coarse_root  + vector(:)
    pool3%sap_wood      = this%sap_wood     + vector(:)
    pool3%heart_wood    = this%heart_wood   + vector(:)
    pool3%labile        = this%labile       + vector(:)
    pool3%reserve       = this%reserve      + vector(:)
    pool3%fruit         = this%fruit        + vector(:)
    pool3%seed_bed      = this%seed_bed     + vector(:)
  END FUNCTION Add_pool_a_and_vector
  !-----------------------------------------------------------------------------------------------------
  !> - (t_veg_pool_a - t_veg_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_pool_a(this,pool2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this, pool2
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         - pool2%leaf
    pool3%fine_root     = this%fine_root    - pool2%fine_root
    pool3%coarse_root   = this%coarse_root  - pool2%coarse_root
    pool3%sap_wood      = this%sap_wood     - pool2%sap_wood
    pool3%heart_wood    = this%heart_wood   - pool2%heart_wood
    pool3%labile        = this%labile       - pool2%labile
    pool3%reserve       = this%reserve      - pool2%reserve
    pool3%fruit         = this%fruit        - pool2%fruit
    pool3%seed_bed      = this%seed_bed     - pool2%seed_bed
  END FUNCTION Subtract_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> - (t_veg_pool_a - scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_pool_a_and_scalar(this,scalar) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this
    REAL(wp),            INTENT(in) :: scalar
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         - scalar
    pool3%fine_root     = this%fine_root    - scalar
    pool3%coarse_root   = this%coarse_root  - scalar
    pool3%sap_wood      = this%sap_wood     - scalar
    pool3%heart_wood    = this%heart_wood   - scalar
    pool3%labile        = this%labile       - scalar
    pool3%reserve       = this%reserve      - scalar
    pool3%fruit         = this%fruit        - scalar
    pool3%seed_bed      = this%seed_bed     - scalar
  END FUNCTION Subtract_pool_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> - (t_veg_pool_a - 1D vector)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_pool_a_and_vector(this,vector) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this
    REAL(wp),            INTENT(in) :: vector(:)
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         - vector(:)
    pool3%fine_root     = this%fine_root    - vector(:)
    pool3%coarse_root   = this%coarse_root  - vector(:)
    pool3%sap_wood      = this%sap_wood     - vector(:)
    pool3%heart_wood    = this%heart_wood   - vector(:)
    pool3%labile        = this%labile       - vector(:)
    pool3%reserve       = this%reserve      - vector(:)
    pool3%fruit         = this%fruit        - vector(:)
    pool3%seed_bed      = this%seed_bed     - vector(:)
  END FUNCTION Subtract_pool_a_and_vector
  !-----------------------------------------------------------------------------------------------------
  !> * (t_veg_pool_a * t_veg_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_pool_a(this,pool2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this, pool2
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         * pool2%leaf
    pool3%fine_root     = this%fine_root    * pool2%fine_root
    pool3%coarse_root   = this%coarse_root  * pool2%coarse_root
    pool3%sap_wood      = this%sap_wood     * pool2%sap_wood
    pool3%heart_wood    = this%heart_wood   * pool2%heart_wood
    pool3%labile        = this%labile       * pool2%labile
    pool3%reserve       = this%reserve      * pool2%reserve
    pool3%fruit         = this%fruit        * pool2%fruit
    pool3%seed_bed      = this%seed_bed     * pool2%seed_bed
  END FUNCTION Multiply_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_veg_pool_a * scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_pool_a_and_scalar(this,scalar) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this
    REAL(wp),            INTENT(in) :: scalar
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         * scalar
    pool3%fine_root     = this%fine_root    * scalar
    pool3%coarse_root   = this%coarse_root  * scalar
    pool3%sap_wood      = this%sap_wood     * scalar
    pool3%heart_wood    = this%heart_wood   * scalar
    pool3%labile        = this%labile       * scalar
    pool3%reserve       = this%reserve      * scalar
    pool3%fruit         = this%fruit        * scalar
    pool3%seed_bed      = this%seed_bed     * scalar
  END FUNCTION Multiply_pool_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> * (scalar * t_veg_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_scalar_and_pool_a(scalar,this2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),            INTENT(in) :: scalar
    CLASS(t_veg_pool_a), INTENT(in) :: this2
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = scalar * this2%leaf
    pool3%fine_root     = scalar * this2%fine_root
    pool3%coarse_root   = scalar * this2%coarse_root
    pool3%sap_wood      = scalar * this2%sap_wood
    pool3%heart_wood    = scalar * this2%heart_wood
    pool3%labile        = scalar * this2%labile
    pool3%reserve       = scalar * this2%reserve
    pool3%fruit         = scalar * this2%fruit
    pool3%seed_bed      = scalar * this2%seed_bed
  END FUNCTION Multiply_scalar_and_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_veg_pool_a * 1D vector)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_pool_a_and_vector(this,vector) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this
    REAL(wp),            INTENT(in) :: vector(:)
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         * vector(:)
    pool3%fine_root     = this%fine_root    * vector(:)
    pool3%coarse_root   = this%coarse_root  * vector(:)
    pool3%sap_wood      = this%sap_wood     * vector(:)
    pool3%heart_wood    = this%heart_wood   * vector(:)
    pool3%labile        = this%labile       * vector(:)
    pool3%reserve       = this%reserve      * vector(:)
    pool3%fruit         = this%fruit        * vector(:)
    pool3%seed_bed      = this%seed_bed     * vector(:)
  END FUNCTION Multiply_pool_a_and_vector
  !-----------------------------------------------------------------------------------------------------
  !> * (1D vector * t_veg_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_vector_and_pool_a(vector,this2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),            INTENT(in) :: vector(:)
    CLASS(t_veg_pool_a), INTENT(in) :: this2
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = vector(:) * this2%leaf
    pool3%fine_root     = vector(:) * this2%fine_root
    pool3%coarse_root   = vector(:) * this2%coarse_root
    pool3%sap_wood      = vector(:) * this2%sap_wood
    pool3%heart_wood    = vector(:) * this2%heart_wood
    pool3%labile        = vector(:) * this2%labile
    pool3%reserve       = vector(:) * this2%reserve
    pool3%fruit         = vector(:) * this2%fruit
    pool3%seed_bed      = vector(:) * this2%seed_bed
  END FUNCTION Multiply_vector_and_pool_a


  !-----------------------------------------------------------------------------------------------------
  !> Initialisation t_slm_struct_veg_element with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE Init_t_slm_struct_veg_element(this,s)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_veg_element),INTENT(inout) :: this 
    REAL(wp),                       INTENT(in)    :: s
    ! 1.0 Go
    this%carbon%ptr(:,:)     = s
    this%nitrogen%ptr(:,:)   = s
    this%phosphorus%ptr(:,:) = s
    this%carbon13%ptr(:,:)   = s
    this%carbon14%ptr(:,:)   = s
    this%nitrogen15%ptr(:,:) = s
  END SUBROUTINE Init_t_slm_struct_veg_element 

  !-----------------------------------------------------------------------------------------------------
  !> Initialisation t_slm_struct_veg_pool with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE Init_t_slm_struct_veg_pool(this,s)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_veg_pool), INTENT(inout) :: this 
    REAL(wp),                     INTENT(in)    :: s
    ! 1.0 Go
    CALL this%leaf%init(s)
    CALL this%fine_root%init(s)
    CALL this%coarse_root%init(s)
    CALL this%sap_wood%init(s)
    CALL this%heart_wood%init(s)
    CALL this%labile%init(s)
    CALL this%reserve%init(s)
    CALL this%fruit%init(s)
    CALL this%seed_bed%init(s)
  END SUBROUTINE Init_t_slm_struct_veg_pool 

  !-----------------------------------------------------------------------------------------------------
  !> Add all variables of an element TYPE to the process memory
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Add_element_variables_to_memory(this, prefix, suffix, lct_ids, model_id, unit)

    USE mo_slm_varlist,       ONLY: jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits
    USE mo_jsb_io,            ONLY: t_cf, t_grib1, t_grib2
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1.InOut
    CLASS(t_slm_struct_veg_element), INTENT(inout) :: this 
    CHARACTER(len=*),                INTENT(in)    :: prefix          !> ?
    CHARACTER(len=*),                INTENT(in)    :: suffix          !> ?
    INTEGER,                         INTENT(in)    :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,                         INTENT(in)    :: model_id        !> model ID model\%id
    CHARACTER(len=*),                INTENT(in)    :: unit            !> unit of the variables in the elements
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_surface                ! Vertical grid
    INTEGER                     :: table
    CHARACTER(len=50)           :: var_name
    CHARACTER(len=100)          :: var_name_long
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_element_variables_to_memory'
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(model_id)
    table         = 1     ! in jsbach4 this line is: table = tables(1)     *slmtemp
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface => Get_vgrid('surface')


    ! run add_var()
    var_name      = TRIM(this%name)//'_carbon'
    var_name_long = 'Carbon in vegetation '//TRIM(this%name)
    CALL jsb_add_var(this, TRIM(var_name), this%carbon, &
          hgrid, vgrid_surface, &
          t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
          t_grib1(table, 255, grib_bits), &
          t_grib2(255, 255, 255, grib_bits), &
          prefix, suffix, &
          initval_r=0.0_wp)     ! initvalue
    var_name      = TRIM(this%name)//'_nitrogen'
    var_name_long = 'Nitrogen in vegetation '//TRIM(this%name)
    CALL jsb_add_var(this, TRIM(var_name), this%nitrogen, &
          hgrid, vgrid_surface, &
          t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
          t_grib1(table, 255, grib_bits), &
          t_grib2(255, 255, 255, grib_bits), &
          prefix, suffix, &
          initval_r=0.0_wp)     ! initvalue
    var_name      = TRIM(this%name)//'_phosphorus'
    var_name_long = 'Phosphorus in vegetation '//TRIM(this%name)
    CALL jsb_add_var(this, TRIM(var_name), this%phosphorus, &
          hgrid, vgrid_surface, &
          t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
          t_grib1(table, 255, grib_bits), &
          t_grib2(255, 255, 255, grib_bits), &
          prefix, suffix, &
          initval_r=0.0_wp)     ! initvalue
    var_name      = TRIM(this%name)//'_carbon13'
    var_name_long = 'Carbon13 in vegetation '//TRIM(this%name)
    CALL jsb_add_var(this, TRIM(var_name), this%carbon13, &
          hgrid, vgrid_surface, &
          t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
          t_grib1(table, 255, grib_bits), &
          t_grib2(255, 255, 255, grib_bits), &
          prefix, suffix, &
          initval_r=0.0_wp)     ! initvalue
    var_name      = TRIM(this%name)//'_carbon14'
    var_name_long = 'Carbon14 in vegetation '//TRIM(this%name)
    CALL jsb_add_var(this, TRIM(var_name), this%carbon14, &
          hgrid, vgrid_surface, &
          t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
          t_grib1(table, 255, grib_bits), &
          t_grib2(255, 255, 255, grib_bits), &
          prefix, suffix, &
          initval_r=0.0_wp)     ! initvalue
    var_name      = TRIM(this%name)//'_nitrogen15'
    var_name_long = 'Nitrogen15 in vegetation '//TRIM(this%name)
    CALL jsb_add_var(this, TRIM(var_name), this%nitrogen15, &
          hgrid, vgrid_surface, &
          t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
          t_grib1(table, 255, grib_bits), &
          t_grib2(255, 255, 255, grib_bits), &
          prefix, suffix, &
          initval_r=0.0_wp)     ! initvalue
  END SUBROUTINE Add_element_variables_to_memory

  !-----------------------------------------------------------------------------------------------------
  !> Add all elements of a pool TYPE to the process memory
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Add_pool_elements_to_memory(this, prefix, suffix, lct_ids, model_id, unit)

    USE mo_slm_varlist,       ONLY: jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits
    USE mo_jsb_io,            ONLY: t_cf, t_grib1, t_grib2
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1.InOut
    CLASS(t_slm_struct_veg_pool),    INTENT(inout) :: this 
    CHARACTER(len=*),                INTENT(in)    :: prefix          !> ?
    CHARACTER(len=*),                INTENT(in)    :: suffix          !> ?
    INTEGER,                         INTENT(in)    :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,                         INTENT(in)    :: model_id        !> model ID model\%id
    CHARACTER(len=*),                INTENT(in)    :: unit            !> unit of the variables in the elements
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_surface                ! Vertical grid
    INTEGER                     :: table
    CHARACTER(len=40)           :: element_name
    CHARACTER(len=80)           :: element_name_long    
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_pool_elements_to_memory'
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(model_id)
    table         = 1     ! in jsbach4 this line is: table = tables(1)     *slmtemp
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface => Get_vgrid('surface')


    ! run add_var()
    element_name      = TRIM(this%name)//'_leaf'
    element_name_long = 'Leaf element of vegetation '//TRIM(this%name)
      CALL jsb_add_var(this, TRIM(element_name), this%leaf, &
            hgrid, vgrid_surface, &
            t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
    element_name      = TRIM(this%name)//'_fine_root'
    element_name_long = 'Fine_root element of vegetation '//TRIM(this%name)
      CALL jsb_add_var(this, TRIM(element_name), this%fine_root, &
            hgrid, vgrid_surface, &
            t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
    element_name      = TRIM(this%name)//'_coarse_root'
    element_name_long = 'Coarse_root element of vegetation '//TRIM(this%name)
      CALL jsb_add_var(this, TRIM(element_name), this%coarse_root, &
            hgrid, vgrid_surface, &
            t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
    element_name      = TRIM(this%name)//'_sap_wood'
    element_name_long = 'Sap_wood element of vegetation '//TRIM(this%name)
      CALL jsb_add_var(this, TRIM(element_name), this%sap_wood, &
            hgrid, vgrid_surface, &
            t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
    element_name      = TRIM(this%name)//'_heart_wood'
    element_name_long = 'Heart_wood element of vegetation '//TRIM(this%name)
      CALL jsb_add_var(this, TRIM(element_name), this%heart_wood, &
            hgrid, vgrid_surface, &
            t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
    element_name      = TRIM(this%name)//'_labile'
    element_name_long = 'Labile element of vegetation '//TRIM(this%name)
      CALL jsb_add_var(this, TRIM(element_name), this%labile, &
            hgrid, vgrid_surface, &
            t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
    element_name      = TRIM(this%name)//'_reserve'
    element_name_long = 'Reserve element of vegetation '//TRIM(this%name)
      CALL jsb_add_var(this, TRIM(element_name), this%reserve, &
            hgrid, vgrid_surface, &
            t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
    element_name      = TRIM(this%name)//'_fruit'
    element_name_long = 'Fruit element of vegetation '//TRIM(this%name)
      CALL jsb_add_var(this, TRIM(element_name), this%fruit, &
            hgrid, vgrid_surface, &
            t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
    element_name      = TRIM(this%name)//'_seed_bed'
    element_name_long = 'Seed_bed element of vegetation '//TRIM(this%name)
      CALL jsb_add_var(this, TRIM(element_name), this%seed_bed, &
            hgrid, vgrid_surface, &
            t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
  END SUBROUTINE Add_pool_elements_to_memory


  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the VEG_ process
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_veg_memory(mem, prefix, suffix, lct_ids, model_id)

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_slm_varlist,       ONLY: jsb_add_var   ! replace by mo_jsb_varlist later    *slmtemp
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2   ! , TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1.InOut
    CLASS(t_veg_memory),  INTENT(inout), TARGET  :: mem           !> pheno memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_surface                ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_canopy                 ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_soil_sb                ! Vertical grid
    INTEGER                    :: table
    CHARACTER(len=30)          :: unit_veg_pool
    CHARACTER(len=30)          :: unit_veg_flux
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_veg_memory'
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(model_id)
    table         = 1     ! in jsbach4 this line is: table = tables(1)     *slmtemp
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface => Get_vgrid('surface')
    vgrid_canopy  => Get_vgrid('canopy_layer')
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')



    !> 0.9 init local variables
    ! 
    unit_veg_pool = "mol m-2"
    unit_veg_flux = "mol m-2 timestep-1"


    !> 1.0  init memory
    !! 
    IF (.NOT. ASSOCIATED(mem%height%ptr)) THEN

      !> add variables to veg memory veg\%VAR
      !! 
      !! jsb_add_var -> jsb_add_var_list_element_r1d    *slmtemp 1D -> 2D
      !!                jsb_add_var_list_element_l1d
      !!  
      CALL jsb_add_var(mem, 'height', mem%height, &
            hgrid, vgrid_surface, &
            t_cf('height', 'm', 'plant height'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'diameter', mem%diameter, &
            hgrid, vgrid_surface, &
            t_cf('diameter', 'm', 'plant diameter'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'dens_ind', mem%dens_ind, &
            hgrid, vgrid_surface, &
            t_cf('dens_ind', '# m-2', 'density of individuals'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'lai', mem%lai, &
            hgrid, vgrid_surface, &
            t_cf('lai', 'm2 m-2', 'actual leaf area index'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_cn_leaf', mem%target_cn_leaf, &
            hgrid, vgrid_surface, &
            t_cf('target_cn_leaf', '', 'target foliar C:N ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_np_leaf', mem%target_np_leaf, &
            hgrid, vgrid_surface, &
            t_cf('target_np_leaf', '', 'target foliar N:P ratio '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'mean_leaf_age', mem%mean_leaf_age, &
            hgrid, vgrid_surface, &
            t_cf('mean_leaf_age', '', 'average foliage age in days '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'cohort_age', mem%cohort_age, &
            hgrid, vgrid_surface, &
            t_cf('cohort_age', '', 'age of the cohort in years, used with veg_dynamics_scheme cohort'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'f_n_demand', mem%f_n_demand, &
            hgrid, vgrid_surface, &
            t_cf('f_n_demand', '', 'relative plant demand for N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'f_p_demand', mem%f_p_demand, &
            hgrid, vgrid_surface, &
            t_cf('f_p_demand', '', 'relative plant demand for P'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_sinklim_ps', mem%beta_sinklim_ps, &
            hgrid, vgrid_surface, &
            t_cf('beta_sinklim_ps', '', 'scaling factor to account for a direct sink limitation contraint on photosynthesis'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'dphi', mem%dphi, &
            hgrid, vgrid_surface, &
            t_cf('dphi', 'MPa', 'plant-soil water potential gradient'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_req_n', mem%growth_req_n, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_n', '', 'moles N required for a unit of C growth under current allocation fractions '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_req_p', mem%growth_req_p, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_p', '', 'moles P required for a unit of N growth under current allocation fractions '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'leaf_cn_direction', mem%leaf_cn_direction, &
            hgrid, vgrid_surface, &
            t_cf('leaf_cn_direction', '', 'indicates direction of change in leaf CN ratio to maximise NPP'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_lai', mem%target_lai, &
            hgrid, vgrid_surface, &
            t_cf('target_lai', 'm2 m-2', 'LAI implied from sap wood area'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'kstar_labile', mem%kstar_labile, &
            hgrid, vgrid_surface, &
            t_cf('kstar_labile', '', 'labile pool turnover rate assuming only temperature and moisture constraints, no stoichiometry'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'do_cohort_harvest', mem%do_cohort_harvest, &
            hgrid, vgrid_surface, &
            t_cf('do_cohort_harvest', '', 'logical to identify whether the cohort may be harvested, in veg_dynamics_scheme cohort'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=false) ! initvalue

      CALL jsb_add_var(mem, 'root_fraction_sl', mem%root_fraction_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('root_fraction_sl', '', 'fractions of root per layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'target_cn_fine_root', mem%target_cn_fine_root, &
            hgrid, vgrid_surface, &
            t_cf('target_cn_fine_root', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_cn_coarse_root', mem%target_cn_coarse_root, &
            hgrid, vgrid_surface, &
            t_cf('target_cn_coarse_root', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_cn_sap_wood', mem%target_cn_sap_wood, &
            hgrid, vgrid_surface, &
            t_cf('target_cn_sap_wood', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_cn_fruit', mem%target_cn_fruit, &
            hgrid, vgrid_surface, &
            t_cf('target_cn_fruit', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_np_fine_root', mem%target_np_fine_root, &
            hgrid, vgrid_surface, &
            t_cf('target_np_fine_root', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_np_coarse_root', mem%target_np_coarse_root, &
            hgrid, vgrid_surface, &
            t_cf('target_np_coarse_root', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_np_sap_wood', mem%target_np_sap_wood, &
            hgrid, vgrid_surface, &
            t_cf('target_np_sap_wood', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_np_fruit', mem%target_np_fruit, &
            hgrid, vgrid_surface, &
            t_cf('target_np_fruit', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'leaf_nitrogen_cl', mem%leaf_nitrogen_cl, &
            hgrid, vgrid_canopy, &
            t_cf('leaf_nitrogen_cl', 'mmol m-2 LAI', 'N in canopy layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fn_chl_cl', mem%fn_chl_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_chl_cl', '', 'fraction of N in chlorophyll, pigments only'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fn_et_cl', mem%fn_et_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_et_cl', '', 'fraction of N in electron transport '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fn_rub_cl', mem%fn_rub_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_rub_cl', '', 'fraction of N in Rubisco '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fn_pepc_cl', mem%fn_pepc_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_pepc_cl', '', 'fraction of N in PEP carboylase, C4 plants only'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fn_oth_cl', mem%fn_oth_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_oth_cl', '', 'fraction of N not associated with photosynthesis '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fleaf_sunlit_cl', mem%fleaf_sunlit_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_cl', '', 'fraction of layer which is sunlit'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'lai_cl', mem%lai_cl, &
            hgrid, vgrid_canopy, &
            t_cf('lai_cl', 'm2 m-2', 'leaf area in the canopy layer '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'cumm_lai_cl', mem%cumm_lai_cl, &
            hgrid, vgrid_canopy, &
            t_cf('cumm_lai_cl', 'm2 m-2', 'cummulative leaf area above the mid-point of the canopy layer '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'maint_respiration_pot', mem%maint_respiration_pot, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration_pot', 'µmol CO2 m-2 s-1', 'potential maintenance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'maint_respiration', mem%maint_respiration, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration', 'µmol CO2 m-2 s-1', 'maintenance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'maint_respiration_c13', mem%maint_respiration_c13, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration_c13', 'µmol 13CO2 m-2 s-1', 'maintenance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'maint_respiration_c14', mem%maint_respiration_c14, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration_c14', 'µmol 14CO2 m-2 s-1', 'maintenance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_respiration', mem%growth_respiration, &
            hgrid, vgrid_surface, &
            t_cf('growth_respiration', 'µmol CO2 m-2 s-1', 'growth respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_respiration_c13', mem%growth_respiration_c13, &
            hgrid, vgrid_surface, &
            t_cf('growth_respiration_c13', 'µmol 13CO2 m-2 s-1', 'growth respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_respiration_c14', mem%growth_respiration_c14, &
            hgrid, vgrid_surface, &
            t_cf('growth_respiration_c14', 'µmol 14CO2 m-2 s-1', 'growth respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n_transform_respiration', mem%n_transform_respiration, &
            hgrid, vgrid_surface, &
            t_cf('n_transform_respiration', 'µmol CO2 m-2 s-1', 'respiration associated with N transformation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n_fixation_respiration', mem%n_fixation_respiration, &
            hgrid, vgrid_surface, &
            t_cf('n_fixation_respiration', 'µmol CO2 m-2 s-1', 'respiration associated with N fixation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n_processing_respiration', mem%n_processing_respiration, &
            hgrid, vgrid_surface, &
            t_cf('n_processing_respiration', 'µmol CO2 m-2 s-1', 'respiration associated with N processing'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n_processing_respiration_c13', mem%n_processing_respiration_c13, &
            hgrid, vgrid_surface, &
            t_cf('n_processing_respiration_c13', 'µmol 13CO2 m-2 s-1', 'respiration associated with N processing'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n_processing_respiration_c14', mem%n_processing_respiration_c14, &
            hgrid, vgrid_surface, &
            t_cf('n_processing_respiration_c14', 'µmol 14CO2 m-2 s-1', 'respiration associated with N processing'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'npp', mem%npp, &
            hgrid, vgrid_surface, &
            t_cf('npp', 'µmol CO2 m-2 s-1', 'net primary production'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'npp_c13', mem%npp_c13, &
            hgrid, vgrid_surface, &
            t_cf('npp_c13', 'µmol 13CO2 m-2 s-1', 'net primary production'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'npp_c14', mem%npp_c14, &
            hgrid, vgrid_surface, &
            t_cf('npp_c14', 'µmol 14CO2 m-2 s-1', 'net primary production'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'uptake_nh4', mem%uptake_nh4, &
            hgrid, vgrid_surface, &
            t_cf('uptake_nh4', 'µmol N m-2 s-1', 'NH4 uptake from soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'uptake_nh4_n15', mem%uptake_nh4_n15, &
            hgrid, vgrid_surface, &
            t_cf('uptake_nh4_n15', 'µmol N m-2 s-1', '15NH4 uptake from soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'uptake_no3', mem%uptake_no3, &
            hgrid, vgrid_surface, &
            t_cf('uptake_no3', 'µmol N m-2 s-1', 'NO3 uptake from soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'uptake_no3_n15', mem%uptake_no3_n15, &
            hgrid, vgrid_surface, &
            t_cf('uptake_no3_n15', 'µmol N m-2 s-1', '15NO3 uptake from soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'uptake_po4', mem%uptake_po4, &
            hgrid, vgrid_surface, &
            t_cf('uptake_po4', 'µmol P m-2 s-1', 'PO4 uptake from soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n_fixation', mem%n_fixation, &
            hgrid, vgrid_surface, &
            t_cf('n_fixation', 'µmol N m-2 s-1', 'symbiontic N fixation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n_fixation_n15', mem%n_fixation_n15, &
            hgrid, vgrid_surface, &
            t_cf('n_fixation_n15', 'µmol 15N m-2 s-1', 'symbiontic N fixation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_leaf_n', mem%recycling_leaf_n, &
            hgrid, vgrid_surface, &
            t_cf('recycling_leaf_n', 'mol m-2 timestep-1', 'net flux of N from leaf to labile (senescence and maintenance)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_leaf_n15', mem%recycling_leaf_n15, &
            hgrid, vgrid_surface, &
            t_cf('recycling_leaf_n15', 'mol m-2 timestep-1', 'net flux of 15N from leaf to labile (senescence and maintenance)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_leaf_p', mem%recycling_leaf_p, &
            hgrid, vgrid_surface, &
            t_cf('recycling_leaf_p', 'mol m-2 timestep-1', 'net flux of P from leaf to labile (senescence and maintenance)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_fine_root_n', mem%recycling_fine_root_n, &
            hgrid, vgrid_surface, &
            t_cf('recycling_fine_root_n', 'mol m-2 timestep-1', 'net flux of N from fine roots to labile (senescence and maintenance)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_fine_root_n15', mem%recycling_fine_root_n15, &
            hgrid, vgrid_surface, &
            t_cf('recycling_fine_root_n15', 'mol m-2 timestep-1', 'net flux of 15N from fine roots to labile (senescence and maintenance)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_fine_root_p', mem%recycling_fine_root_p, &
            hgrid, vgrid_surface, &
            t_cf('recycling_fine_root_p', 'mol m-2 timestep-1', 'net flux of P from fine roots to labile (senescence and maintenance)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_heart_wood_n', mem%recycling_heart_wood_n, &
            hgrid, vgrid_surface, &
            t_cf('recycling_heart_wood_n', 'mol m-2 timestep-1', 'net flux of N from sapwood to labile (senescence and maintenance)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_heart_wood_n15', mem%recycling_heart_wood_n15, &
            hgrid, vgrid_surface, &
            t_cf('recycling_heart_wood_n15', 'mol m-2 timestep-1', 'net flux of 15N from sapwood to labile (senescence and maintenance)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_heart_wood_p', mem%recycling_heart_wood_p, &
            hgrid, vgrid_surface, &
            t_cf('recycling_heart_wood_p', 'mol m-2 timestep-1', 'net flux of P from sapwood to labile (senescence and maintenance)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_growth', mem%net_growth, &
            hgrid, vgrid_surface, &
            t_cf('net_growth', 'µmol CO2 m-2 s-1', 'npp minus turnover'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fmaint_rate_root', mem%fmaint_rate_root, &
            hgrid, vgrid_surface, &
            t_cf('fmaint_rate_root', '~µmol C mol-1 N s-1', 'current N-specific maintentance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_npp', mem%unit_npp, &
            hgrid, vgrid_surface, &
            t_cf('unit_npp', '', 'incremental change in NPP. Needed for optimal allocation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_transpiration', mem%unit_transpiration, &
            hgrid, vgrid_surface, &
            t_cf('unit_transpiration', '', 'incremental change in transpiration. Needed for optimal allocation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_n_pot', mem%unit_uptake_n_pot, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_n_pot', '', 'incremental change in potential N uptake per unit root carbon, needed for N uptake calculation '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_p_pot', mem%unit_uptake_p_pot, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_p_pot', '', 'incremental change in potential P uptake per unit root carbon, currently unused '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_n_act', mem%unit_uptake_n_act, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_n_act', '', 'incremental change in actual N uptake per unit root carbon, needed for optimal allocation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_p_act', mem%unit_uptake_p_act, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_p_act', '', 'incremental change in actual P uptake per unit root carbon, needed for optimal allocation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'delta_dens_ind', mem%delta_dens_ind, &
            hgrid, vgrid_surface, &
            t_cf('delta_dens_ind', '# m-2 timestep-1', 'change in density of individuals'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'cost_n_uptake_root', mem%cost_n_uptake_root, &
            hgrid, vgrid_surface, &
            t_cf('cost_n_uptake_root', 'mol C mol-1 N', 'actual cost for N uptake by roots'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'cost_n_fixation', mem%cost_n_fixation, &
            hgrid, vgrid_surface, &
            t_cf('cost_n_fixation', 'mol C mol-1 N', 'actual cost for N uptake through fixation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'mortality_rate', mem%mortality_rate, &
            hgrid, vgrid_surface, &
            t_cf('mortality_rate', 'timestep-1', 'background_mort_rate_tree / grass'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'k1_opt', mem%k1_opt, &
            hgrid, vgrid_surface, &
            t_cf('k1_opt', '', 'coefficient for optimal leaf:root ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'k2_opt', mem%k2_opt, &
            hgrid, vgrid_surface, &
            t_cf('k2_opt', '', 'coefficient for optimal leaf:root ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'leaf2sapwood_mass_ratio', mem%leaf2sapwood_mass_ratio, &
            hgrid, vgrid_surface, &
            t_cf('leaf2sapwood_mass_ratio', '', 'leaf to shoot mass ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'leaf2root_mass_ratio', mem%leaf2root_mass_ratio, &
            hgrid, vgrid_surface, &
            t_cf('leaf2root_mass_ratio', '', 'leaf to root mass ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fleaf_sunlit_daytime_cl', mem%fleaf_sunlit_daytime_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_daytime_cl', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_daytime', mem%t_air_daytime, &
            hgrid, vgrid_surface, &
            t_cf('t_air_daytime', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'press_srf_daytime', mem%press_srf_daytime, &
            hgrid, vgrid_surface, &
            t_cf('press_srf_daytime', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'co2_mixing_ratio_daytime', mem%co2_mixing_ratio_daytime, &
            hgrid, vgrid_surface, &
            t_cf('co2_mixing_ratio_daytime', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ga_daytime', mem%ga_daytime, &
            hgrid, vgrid_surface, &
            t_cf('ga_daytime', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_sinklim_ps_daytime', mem%beta_sinklim_ps_daytime, &
            hgrid, vgrid_surface, &
            t_cf('beta_sinklim_ps_daytime', 'unitless', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_jmax_opt_daytime', mem%t_jmax_opt_daytime, &
            hgrid, vgrid_surface, &
            t_cf('t_jmax_opt_daytime', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fleaf_sunlit_daytime_dacc_cl', mem%fleaf_sunlit_daytime_dacc_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_daytime_dacc_cl', '', 'value accumulated over current day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_daytime_dacc', mem%t_air_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('t_air_daytime_dacc', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'press_srf_daytime_dacc', mem%press_srf_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('press_srf_daytime_dacc', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'co2_mixing_ratio_daytime_dacc', mem%co2_mixing_ratio_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('co2_mixing_ratio_daytime_dacc', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ga_daytime_dacc', mem%ga_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('ga_daytime_dacc', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_sinklim_ps_daytime_dacc', mem%beta_sinklim_ps_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('beta_sinklim_ps_daytime_dacc', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_jmax_opt_daytime_dacc', mem%t_jmax_opt_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('t_jmax_opt_daytime_dacc', '', 'value accumulated over current day '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'daytime_counter', mem%daytime_counter, &
            hgrid, vgrid_surface, &
            t_cf('daytime_counter', '', 'number of timesteps of current day with light available'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'gpp_tlabile_mavg', mem%gpp_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('gpp_tlabile_mavg', 'µmol m-2 s-1', 'gpp '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'maint_respiration_tlabile_mavg', mem%maint_respiration_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration_tlabile_mavg', 'µmol m-2 s-1', 'maintenance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_req_n_tlabile_mavg', mem%growth_req_n_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_n_tlabile_mavg', '', 'moles N required for a unit of C growth under current allocation fractions'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_req_p_tlabile_mavg', mem%growth_req_p_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_p_tlabile_mavg', '', 'moles P required for a unit of N growth under current allocation fractions'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_tphen_mavg', mem%t_air_tphen_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_tphen_mavg', '', 'air temp mavg for phenology calc'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'npp_tuptake_mavg', mem%npp_tuptake_mavg, &
            hgrid, vgrid_surface, &
            t_cf('npp_tuptake_mavg', 'µmol m-2 s-1', 'npp of the last week'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'demand_uptake_n_tuptake_mavg', mem%demand_uptake_n_tuptake_mavg, &
            hgrid, vgrid_surface, &
            t_cf('demand_uptake_n_tuptake_mavg', '', 'scalar of Vmax for root N uptake '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'demand_uptake_p_tuptake_mavg', mem%demand_uptake_p_tuptake_mavg, &
            hgrid, vgrid_surface, &
            t_cf('demand_uptake_p_tuptake_mavg', '', 'scalar of Vmax for root P uptake'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_req_n_tuptake_mavg', mem%growth_req_n_tuptake_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_n_tuptake_mavg', '', 'moles N required for a unit of C growth under current allocation fractions'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_req_p_tuptake_mavg', mem%growth_req_p_tuptake_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_p_tuptake_mavg', '', 'moles P required for a unit of N growth under current allocation fractions'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'cost_n_fixation_tuptake_mavg', mem%cost_n_fixation_tuptake_mavg, &
            hgrid, vgrid_surface, &
            t_cf('cost_n_fixation_tuptake_mavg', 'mol C mol-1 N', 'moving average of N fixation cost'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fleaf_sunlit_tfrac_mavg_cl', mem%fleaf_sunlit_tfrac_mavg_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_tfrac_mavg_cl', '', 'fraction of sunlit leaves'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_tfrac_mavg', mem%t_air_tfrac_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_tfrac_mavg', '', 'air temp'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'press_srf_tfrac_mavg', mem%press_srf_tfrac_mavg, &
            hgrid, vgrid_surface, &
            t_cf('press_srf_tfrac_mavg', '', 'air pressure at surface'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'co2_mixing_ratio_tfrac_mavg', mem%co2_mixing_ratio_tfrac_mavg, &
            hgrid, vgrid_surface, &
            t_cf('co2_mixing_ratio_tfrac_mavg', '', 'co2 mixing ratio of the air'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ga_tfrac_mavg', mem%ga_tfrac_mavg, &
            hgrid, vgrid_surface, &
            t_cf('ga_tfrac_mavg', '', 'aerodynamic conductance'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_sinklim_ps_tfrac_mavg', mem%beta_sinklim_ps_tfrac_mavg, &
            hgrid, vgrid_surface, &
            t_cf('beta_sinklim_ps_tfrac_mavg', '', 'scaling factor to account for soil moisture constraints on photosynthesis '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fleaf_sunlit_tcnl_mavg_cl', mem%fleaf_sunlit_tcnl_mavg_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_tcnl_mavg_cl', '', 'fraction of sunlit leaves'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_tcnl_mavg', mem%t_air_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_tcnl_mavg', '', 'air temp'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'press_srf_tcnl_mavg', mem%press_srf_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('press_srf_tcnl_mavg', '', 'air pressure at surface'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'co2_mixing_ratio_tcnl_mavg', mem%co2_mixing_ratio_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('co2_mixing_ratio_tcnl_mavg', '', 'co2 mixing ratio of the air'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ga_tcnl_mavg', mem%ga_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('ga_tcnl_mavg', '', 'aerodynamic conductance'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_sinklim_ps_tcnl_mavg', mem%beta_sinklim_ps_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('beta_sinklim_ps_tcnl_mavg', '', 'scaling factor to account for soil moisture constraints on conductance'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'uptake_n_tcnl_mavg', mem%uptake_n_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('uptake_n_tcnl_mavg', '', 'root N uptake'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'growth_cn_tcnl_mavg', mem%growth_cn_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_cn_tcnl_mavg', '', 'new growth CN ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'growth_np_tcnl_mavg', mem%growth_np_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_np_tcnl_mavg', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'npp_tcnl_mavg', mem%npp_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('npp_tcnl_mavg', '', 'NPP'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'fmaint_rate_tcnl_mavg', mem%fmaint_rate_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('fmaint_rate_tcnl_mavg', '~µmol C mol-1 N s-1', 'N-specific maintentance respiration for optmal leaf N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'labile_nitrogen_tcnl_mavg', mem%labile_nitrogen_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('labile_nitrogen_tcnl_mavg', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'labile_carbon_tcnl_mavg', mem%labile_carbon_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('labile_carbon_tcnl_mavg', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'labile_phosphorus_tcnl_mavg', mem%labile_phosphorus_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('labile_phosphorus_tcnl_mavg', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
  
      CALL jsb_add_var(mem, 'npp_talloc_mavg', mem%npp_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('npp_talloc_mavg', 'µmol m-2 s-1', 'NPP'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'n_fixation_talloc_mavg', mem%n_fixation_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('n_fixation_talloc_mavg', '', 'N fixation averaged at the allocation timescale'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transpiration_talloc_mavg', mem%transpiration_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('transpiration_talloc_mavg', '', 'average transpiration at allocation time-scale'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_npp_talloc_mavg', mem%unit_npp_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('unit_npp_talloc_mavg', 'µmol C mol-1', 'rate of change in NPP with  increment in leaf mass'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_transpiration_talloc_mavg', mem%unit_transpiration_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('unit_transpiration_talloc_mavg', '', 'rate of tranpiration change per unit LAI change?'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_n_talloc_mavg', mem%unit_uptake_n_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_n_talloc_mavg', 'µmol N mol-1', 'actual rate of N uptkae per unit root'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_p_talloc_mavg', mem%unit_uptake_p_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_p_talloc_mavg', 'µmol P mol-1', 'actual rate of P uptake per unit root'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'dphi_talloc_mavg', mem%dphi_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('dphi_talloc_mavg', 'MPa', 'leaf-soil water potential gradient'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_cn_talloc_mavg', mem%growth_cn_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_cn_talloc_mavg', '', 'plant growth cn ratio, equivalent to 1/growth_req_n'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_cp_talloc_mavg', mem%growth_cp_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_cp_talloc_mavg', '', 'plant growth cp ratio, equivalent to 1/(growthreq_n*growth_req_p)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'growth_np_talloc_mavg', mem%growth_np_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_np_talloc_mavg', '', 'plant growth hp ratio, equivalent to 1/(growth_req_p)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'labile_carbon_talloc_mavg', mem%labile_carbon_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('labile_carbon_talloc_mavg', '', 'moving average of labile carbon pool neeeded for allometry calculations'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'labile_nitrogen_talloc_mavg', mem%labile_nitrogen_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('labile_nitrogen_talloc_mavg', '', 'moving average of labile nitrogen pool neeeded for allometry calculations'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'labile_phosphorus_talloc_mavg', mem%labile_phosphorus_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('labile_phosphorus_talloc_mavg', '', 'moving average of labile phosphorus pool neeeded for allometry calculations'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'reserve_carbon_talloc_mavg', mem%reserve_carbon_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('reserve_carbon_talloc_mavg', '', 'moving average of reserve carbon needed for optimal allocation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'reserve_nitrogen_talloc_mavg', mem%reserve_nitrogen_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('reserve_nitrogen_talloc_mavg', '', 'moving average of reserve nitrogen needed for optimal allocation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'reserve_phosphorus_talloc_mavg', mem%reserve_phosphorus_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('reserve_phosphorus_talloc_mavg', '', 'moving average of reserve phosphorus needed for optimal allocation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_root_lim_talloc_mavg', mem%w_root_lim_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('w_root_lim_talloc_mavg', '', 'moving average of ...'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=1.0_wp) ! initvalue
      
      CALL jsb_add_var(mem, 'leaf2root_troot_mavg', mem%leaf2root_troot_mavg, &
            hgrid, vgrid_surface, &
            t_cf('leaf2root_troot_mavg', '', 'leaf to root mass ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_npp_troot_mavg', mem%unit_npp_troot_mavg, &
            hgrid, vgrid_surface, &
            t_cf('unit_npp_troot_mavg', 'µmol C mol-1', 'rate of change in NPP with increment in leaf mass'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_n_pot_troot_mavg', mem%unit_uptake_n_pot_troot_mavg, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_n_pot_troot_mavg', 'µmol N mol-1', 'rate of N uptkae per unit root'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_p_pot_troot_mavg', mem%unit_uptake_p_pot_troot_mavg, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_p_pot_troot_mavg', 'µmol P mol-1', 'rate of P uptake per unit root'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fmaint_rate_troot_mavg', mem%fmaint_rate_troot_mavg, &
            hgrid, vgrid_surface, &
            t_cf('fmaint_rate_troot_mavg', '~µmol C mol-1 N s-1', 'moving average of N-specific maintentance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'an_boc_tvegdyn_mavg', mem%an_boc_tvegdyn_mavg, &
            hgrid, vgrid_surface, &
            t_cf('an_boc_tvegdyn_mavg', '', 'net photosynthesis at bottom of canopy, lowest canopy layer existing'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_growth_tvegdyn_mavg', mem%net_growth_tvegdyn_mavg, &
            hgrid, vgrid_surface, &
            t_cf('net_growth_tvegdyn_mavg', 'µmol C m-2 s-1', 'net growth, npp minus turnover'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'lai_tvegdyn_mavg', mem%lai_tvegdyn_mavg, &
            hgrid, vgrid_surface, &
            t_cf('lai_tvegdyn_mavg', 'm2 m-2', 'lai at vegdynamics timescale'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_tacclim_mavg', mem%t_air_tacclim_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_tacclim_mavg', '', 'average air temperature for respiration acclimation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_soil_root_tacclim_mavg', mem%t_soil_root_tacclim_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_soil_root_tacclim_mavg', '', 'average root temperature for respiration acclimation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_week_mavg', mem%t_air_week_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_week_mavg', '', 'avrg air temperature across 7 days'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_month_mavg', mem%t_air_month_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_month_mavg', '', 'avrg air temperature across 30 days'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_jmax_opt_mavg', mem%t_jmax_opt_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_jmax_opt_mavg', '', 'optimim temperature for electron transport of photosynthesis'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue


      !-----------------------------------------------------------------------------------------------------
      !> veg elements TYPE (t_jsb_struct)
      !! 
      !! jsb_add_var -> slm_add_struct()
      !! 
      !-----------------------------------------------------------------------------------------------------
      CALL jsb_add_var(mem, 'exudation',              mem%exudation, &
            hgrid, vgrid_surface, &
            t_cf('exudation', TRIM(unit_veg_flux), 'Exudation element in vegetation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'establishment',          mem%establishment, &
            hgrid, vgrid_surface, &
            t_cf('establishment', TRIM(unit_veg_flux), 'Establishment element in vegetation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'reserve_use',            mem%reserve_use, &
            hgrid, vgrid_surface, &
            t_cf('reserve_use', TRIM(unit_veg_flux), 'Reserve_use element in vegetation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'total_biomass',          mem%total_biomass, &
            hgrid, vgrid_surface, &
            t_cf('total_biomass', TRIM(unit_veg_pool), 'Total_biomass element in vegetation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)

      !> add variables to the elements
      !! 
      !! veg\%elements\%VAR
      !! 
      CALL mem%exudation%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%establishment%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%reserve_use%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%total_biomass%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_pool))


      !-----------------------------------------------------------------------------------------------------
      !> veg pools (nested) TYPE (t_jsb_n_struct)
      !! 
      !! jsb_add_var -> slm_add_n_struct()
      !! 
      !-----------------------------------------------------------------------------------------------------
      CALL jsb_add_var(mem, 'pool'        , mem%pool, &
            hgrid, vgrid_surface, &
            t_cf('pool', TRIM(unit_veg_pool), 'Pool vegetation pool'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'growth'      , mem%growth, &
            hgrid, vgrid_surface, &
            t_cf('growth', TRIM(unit_veg_flux), 'Growth vegetation flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'litterfall'  , mem%litterfall, &
            hgrid, vgrid_surface, &
            t_cf('litterfall', TRIM(unit_veg_flux), 'Litterfall vegetation flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'falloc'      , mem%falloc, &
            hgrid, vgrid_surface, &
            t_cf('falloc', TRIM(unit_veg_flux), 'Falloc vegetation flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'total_flux'  , mem%total_flux, &
            hgrid, vgrid_surface, &
            t_cf('total_flux', TRIM(unit_veg_flux), 'Total_flux vegetation flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)     ! mass conservation diagnosis
            
      !> add element TYPEs to the pools
      !! 
      !! jsb_add_var -> slm_add_l2_struct()
      !! 
      CALL mem%pool%Add_elements(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_pool))
      CALL mem%growth%Add_elements(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%litterfall%Add_elements(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%falloc%Add_elements(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%total_flux%Add_elements(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))

      !> add variables to the elements in the pools
      !!
      !! veg\%pool\%elements\%VAR
      !! 
      CALL mem%pool%leaf%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_pool))
      CALL mem%pool%fine_root%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_pool))
      CALL mem%pool%coarse_root%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_pool))
      CALL mem%pool%sap_wood%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_pool))
      CALL mem%pool%heart_wood%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_pool))
      CALL mem%pool%labile%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_pool))
      CALL mem%pool%reserve%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_pool))
      CALL mem%pool%fruit%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_pool))
      CALL mem%pool%seed_bed%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_pool))
 
      CALL mem%growth%leaf%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%growth%fine_root%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%growth%coarse_root%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%growth%sap_wood%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%growth%heart_wood%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%growth%labile%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%growth%reserve%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%growth%fruit%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%growth%seed_bed%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
       
      CALL mem%litterfall%leaf%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%litterfall%fine_root%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%litterfall%coarse_root%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%litterfall%sap_wood%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%litterfall%heart_wood%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%litterfall%labile%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%litterfall%reserve%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%litterfall%fruit%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%litterfall%seed_bed%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      
      CALL mem%falloc%leaf%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%falloc%fine_root%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%falloc%coarse_root%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%falloc%sap_wood%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%falloc%heart_wood%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%falloc%labile%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%falloc%reserve%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%falloc%fruit%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%falloc%seed_bed%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      
      CALL mem%total_flux%leaf%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%total_flux%fine_root%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%total_flux%coarse_root%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%total_flux%sap_wood%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%total_flux%heart_wood%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%total_flux%labile%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%total_flux%reserve%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%total_flux%fruit%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%total_flux%seed_bed%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))

    END IF

  END SUBROUTINE Init_veg_memory


  !-----------------------------------------------------------------------------------------------------
  !> ALLOCATE variables in veg_pool_a and veg_element_a structures
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_veg_pools_elements_a(ics, ice, nc, iblk)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER, INTENT(in)     :: ics, ice, nc, iblk             !< information about dimensions of variables in memory
    ! ---------------------------
    ! 0.2 Local
    INTEGER                 :: npoints                        !< number of points in the chunk

    ! set dimension to number of points in this chunk
    npoints = nc


    !< 1.0  ALLOCATE variables
    !! 
    ALLOCATE(veg_pool%leaf%carbon(npoints))
    ALLOCATE(veg_pool%leaf%nitrogen(npoints))
    ALLOCATE(veg_pool%leaf%phosphorus(npoints))
    ALLOCATE(veg_pool%leaf%carbon13(npoints))
    ALLOCATE(veg_pool%leaf%carbon14(npoints))
    ALLOCATE(veg_pool%leaf%nitrogen15(npoints))
    ALLOCATE(veg_pool%fine_root%carbon(npoints))
    ALLOCATE(veg_pool%fine_root%nitrogen(npoints))
    ALLOCATE(veg_pool%fine_root%phosphorus(npoints))
    ALLOCATE(veg_pool%fine_root%carbon13(npoints))
    ALLOCATE(veg_pool%fine_root%carbon14(npoints))
    ALLOCATE(veg_pool%fine_root%nitrogen15(npoints))
    ALLOCATE(veg_pool%coarse_root%carbon(npoints))
    ALLOCATE(veg_pool%coarse_root%nitrogen(npoints))
    ALLOCATE(veg_pool%coarse_root%phosphorus(npoints))
    ALLOCATE(veg_pool%coarse_root%carbon13(npoints))
    ALLOCATE(veg_pool%coarse_root%carbon14(npoints))
    ALLOCATE(veg_pool%coarse_root%nitrogen15(npoints))
    ALLOCATE(veg_pool%sap_wood%carbon(npoints))
    ALLOCATE(veg_pool%sap_wood%nitrogen(npoints))
    ALLOCATE(veg_pool%sap_wood%phosphorus(npoints))
    ALLOCATE(veg_pool%sap_wood%carbon13(npoints))
    ALLOCATE(veg_pool%sap_wood%carbon14(npoints))
    ALLOCATE(veg_pool%sap_wood%nitrogen15(npoints))
    ALLOCATE(veg_pool%heart_wood%carbon(npoints))
    ALLOCATE(veg_pool%heart_wood%nitrogen(npoints))
    ALLOCATE(veg_pool%heart_wood%phosphorus(npoints))
    ALLOCATE(veg_pool%heart_wood%carbon13(npoints))
    ALLOCATE(veg_pool%heart_wood%carbon14(npoints))
    ALLOCATE(veg_pool%heart_wood%nitrogen15(npoints))
    ALLOCATE(veg_pool%labile%carbon(npoints))
    ALLOCATE(veg_pool%labile%nitrogen(npoints))
    ALLOCATE(veg_pool%labile%phosphorus(npoints))
    ALLOCATE(veg_pool%labile%carbon13(npoints))
    ALLOCATE(veg_pool%labile%carbon14(npoints))
    ALLOCATE(veg_pool%labile%nitrogen15(npoints))
    ALLOCATE(veg_pool%reserve%carbon(npoints))
    ALLOCATE(veg_pool%reserve%nitrogen(npoints))
    ALLOCATE(veg_pool%reserve%phosphorus(npoints))
    ALLOCATE(veg_pool%reserve%carbon13(npoints))
    ALLOCATE(veg_pool%reserve%carbon14(npoints))
    ALLOCATE(veg_pool%reserve%nitrogen15(npoints))
    ALLOCATE(veg_pool%fruit%carbon(npoints))
    ALLOCATE(veg_pool%fruit%nitrogen(npoints))
    ALLOCATE(veg_pool%fruit%phosphorus(npoints))
    ALLOCATE(veg_pool%fruit%carbon13(npoints))
    ALLOCATE(veg_pool%fruit%carbon14(npoints))
    ALLOCATE(veg_pool%fruit%nitrogen15(npoints))
    ALLOCATE(veg_pool%seed_bed%carbon(npoints))
    ALLOCATE(veg_pool%seed_bed%nitrogen(npoints))
    ALLOCATE(veg_pool%seed_bed%phosphorus(npoints))
    ALLOCATE(veg_pool%seed_bed%carbon13(npoints))
    ALLOCATE(veg_pool%seed_bed%carbon14(npoints))
    ALLOCATE(veg_pool%seed_bed%nitrogen15(npoints))

    ALLOCATE(veg_growth%leaf%carbon(npoints))
    ALLOCATE(veg_growth%leaf%nitrogen(npoints))
    ALLOCATE(veg_growth%leaf%phosphorus(npoints))
    ALLOCATE(veg_growth%leaf%carbon13(npoints))
    ALLOCATE(veg_growth%leaf%carbon14(npoints))
    ALLOCATE(veg_growth%leaf%nitrogen15(npoints))
    ALLOCATE(veg_growth%fine_root%carbon(npoints))
    ALLOCATE(veg_growth%fine_root%nitrogen(npoints))
    ALLOCATE(veg_growth%fine_root%phosphorus(npoints))
    ALLOCATE(veg_growth%fine_root%carbon13(npoints))
    ALLOCATE(veg_growth%fine_root%carbon14(npoints))
    ALLOCATE(veg_growth%fine_root%nitrogen15(npoints))
    ALLOCATE(veg_growth%coarse_root%carbon(npoints))
    ALLOCATE(veg_growth%coarse_root%nitrogen(npoints))
    ALLOCATE(veg_growth%coarse_root%phosphorus(npoints))
    ALLOCATE(veg_growth%coarse_root%carbon13(npoints))
    ALLOCATE(veg_growth%coarse_root%carbon14(npoints))
    ALLOCATE(veg_growth%coarse_root%nitrogen15(npoints))
    ALLOCATE(veg_growth%sap_wood%carbon(npoints))
    ALLOCATE(veg_growth%sap_wood%nitrogen(npoints))
    ALLOCATE(veg_growth%sap_wood%phosphorus(npoints))
    ALLOCATE(veg_growth%sap_wood%carbon13(npoints))
    ALLOCATE(veg_growth%sap_wood%carbon14(npoints))
    ALLOCATE(veg_growth%sap_wood%nitrogen15(npoints))
    ALLOCATE(veg_growth%heart_wood%carbon(npoints))
    ALLOCATE(veg_growth%heart_wood%nitrogen(npoints))
    ALLOCATE(veg_growth%heart_wood%phosphorus(npoints))
    ALLOCATE(veg_growth%heart_wood%carbon13(npoints))
    ALLOCATE(veg_growth%heart_wood%carbon14(npoints))
    ALLOCATE(veg_growth%heart_wood%nitrogen15(npoints))
    ALLOCATE(veg_growth%labile%carbon(npoints))
    ALLOCATE(veg_growth%labile%nitrogen(npoints))
    ALLOCATE(veg_growth%labile%phosphorus(npoints))
    ALLOCATE(veg_growth%labile%carbon13(npoints))
    ALLOCATE(veg_growth%labile%carbon14(npoints))
    ALLOCATE(veg_growth%labile%nitrogen15(npoints))
    ALLOCATE(veg_growth%reserve%carbon(npoints))
    ALLOCATE(veg_growth%reserve%nitrogen(npoints))
    ALLOCATE(veg_growth%reserve%phosphorus(npoints))
    ALLOCATE(veg_growth%reserve%carbon13(npoints))
    ALLOCATE(veg_growth%reserve%carbon14(npoints))
    ALLOCATE(veg_growth%reserve%nitrogen15(npoints))
    ALLOCATE(veg_growth%fruit%carbon(npoints))
    ALLOCATE(veg_growth%fruit%nitrogen(npoints))
    ALLOCATE(veg_growth%fruit%phosphorus(npoints))
    ALLOCATE(veg_growth%fruit%carbon13(npoints))
    ALLOCATE(veg_growth%fruit%carbon14(npoints))
    ALLOCATE(veg_growth%fruit%nitrogen15(npoints))
    ALLOCATE(veg_growth%seed_bed%carbon(npoints))
    ALLOCATE(veg_growth%seed_bed%nitrogen(npoints))
    ALLOCATE(veg_growth%seed_bed%phosphorus(npoints))
    ALLOCATE(veg_growth%seed_bed%carbon13(npoints))
    ALLOCATE(veg_growth%seed_bed%carbon14(npoints))
    ALLOCATE(veg_growth%seed_bed%nitrogen15(npoints))

    ALLOCATE(veg_litterfall%leaf%carbon(npoints))
    ALLOCATE(veg_litterfall%leaf%nitrogen(npoints))
    ALLOCATE(veg_litterfall%leaf%phosphorus(npoints))
    ALLOCATE(veg_litterfall%leaf%carbon13(npoints))
    ALLOCATE(veg_litterfall%leaf%carbon14(npoints))
    ALLOCATE(veg_litterfall%leaf%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%fine_root%carbon(npoints))
    ALLOCATE(veg_litterfall%fine_root%nitrogen(npoints))
    ALLOCATE(veg_litterfall%fine_root%phosphorus(npoints))
    ALLOCATE(veg_litterfall%fine_root%carbon13(npoints))
    ALLOCATE(veg_litterfall%fine_root%carbon14(npoints))
    ALLOCATE(veg_litterfall%fine_root%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%coarse_root%carbon(npoints))
    ALLOCATE(veg_litterfall%coarse_root%nitrogen(npoints))
    ALLOCATE(veg_litterfall%coarse_root%phosphorus(npoints))
    ALLOCATE(veg_litterfall%coarse_root%carbon13(npoints))
    ALLOCATE(veg_litterfall%coarse_root%carbon14(npoints))
    ALLOCATE(veg_litterfall%coarse_root%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%sap_wood%carbon(npoints))
    ALLOCATE(veg_litterfall%sap_wood%nitrogen(npoints))
    ALLOCATE(veg_litterfall%sap_wood%phosphorus(npoints))
    ALLOCATE(veg_litterfall%sap_wood%carbon13(npoints))
    ALLOCATE(veg_litterfall%sap_wood%carbon14(npoints))
    ALLOCATE(veg_litterfall%sap_wood%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%heart_wood%carbon(npoints))
    ALLOCATE(veg_litterfall%heart_wood%nitrogen(npoints))
    ALLOCATE(veg_litterfall%heart_wood%phosphorus(npoints))
    ALLOCATE(veg_litterfall%heart_wood%carbon13(npoints))
    ALLOCATE(veg_litterfall%heart_wood%carbon14(npoints))
    ALLOCATE(veg_litterfall%heart_wood%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%labile%carbon(npoints))
    ALLOCATE(veg_litterfall%labile%nitrogen(npoints))
    ALLOCATE(veg_litterfall%labile%phosphorus(npoints))
    ALLOCATE(veg_litterfall%labile%carbon13(npoints))
    ALLOCATE(veg_litterfall%labile%carbon14(npoints))
    ALLOCATE(veg_litterfall%labile%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%reserve%carbon(npoints))
    ALLOCATE(veg_litterfall%reserve%nitrogen(npoints))
    ALLOCATE(veg_litterfall%reserve%phosphorus(npoints))
    ALLOCATE(veg_litterfall%reserve%carbon13(npoints))
    ALLOCATE(veg_litterfall%reserve%carbon14(npoints))
    ALLOCATE(veg_litterfall%reserve%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%fruit%carbon(npoints))
    ALLOCATE(veg_litterfall%fruit%nitrogen(npoints))
    ALLOCATE(veg_litterfall%fruit%phosphorus(npoints))
    ALLOCATE(veg_litterfall%fruit%carbon13(npoints))
    ALLOCATE(veg_litterfall%fruit%carbon14(npoints))
    ALLOCATE(veg_litterfall%fruit%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%seed_bed%carbon(npoints))
    ALLOCATE(veg_litterfall%seed_bed%nitrogen(npoints))
    ALLOCATE(veg_litterfall%seed_bed%phosphorus(npoints))
    ALLOCATE(veg_litterfall%seed_bed%carbon13(npoints))
    ALLOCATE(veg_litterfall%seed_bed%carbon14(npoints))
    ALLOCATE(veg_litterfall%seed_bed%nitrogen15(npoints))

    ALLOCATE(veg_falloc%leaf%carbon(npoints))
    ALLOCATE(veg_falloc%leaf%nitrogen(npoints))
    ALLOCATE(veg_falloc%leaf%phosphorus(npoints))
    ALLOCATE(veg_falloc%leaf%carbon13(npoints))
    ALLOCATE(veg_falloc%leaf%carbon14(npoints))
    ALLOCATE(veg_falloc%leaf%nitrogen15(npoints))
    ALLOCATE(veg_falloc%fine_root%carbon(npoints))
    ALLOCATE(veg_falloc%fine_root%nitrogen(npoints))
    ALLOCATE(veg_falloc%fine_root%phosphorus(npoints))
    ALLOCATE(veg_falloc%fine_root%carbon13(npoints))
    ALLOCATE(veg_falloc%fine_root%carbon14(npoints))
    ALLOCATE(veg_falloc%fine_root%nitrogen15(npoints))
    ALLOCATE(veg_falloc%coarse_root%carbon(npoints))
    ALLOCATE(veg_falloc%coarse_root%nitrogen(npoints))
    ALLOCATE(veg_falloc%coarse_root%phosphorus(npoints))
    ALLOCATE(veg_falloc%coarse_root%carbon13(npoints))
    ALLOCATE(veg_falloc%coarse_root%carbon14(npoints))
    ALLOCATE(veg_falloc%coarse_root%nitrogen15(npoints))
    ALLOCATE(veg_falloc%sap_wood%carbon(npoints))
    ALLOCATE(veg_falloc%sap_wood%nitrogen(npoints))
    ALLOCATE(veg_falloc%sap_wood%phosphorus(npoints))
    ALLOCATE(veg_falloc%sap_wood%carbon13(npoints))
    ALLOCATE(veg_falloc%sap_wood%carbon14(npoints))
    ALLOCATE(veg_falloc%sap_wood%nitrogen15(npoints))
    ALLOCATE(veg_falloc%heart_wood%carbon(npoints))
    ALLOCATE(veg_falloc%heart_wood%nitrogen(npoints))
    ALLOCATE(veg_falloc%heart_wood%phosphorus(npoints))
    ALLOCATE(veg_falloc%heart_wood%carbon13(npoints))
    ALLOCATE(veg_falloc%heart_wood%carbon14(npoints))
    ALLOCATE(veg_falloc%heart_wood%nitrogen15(npoints))
    ALLOCATE(veg_falloc%labile%carbon(npoints))
    ALLOCATE(veg_falloc%labile%nitrogen(npoints))
    ALLOCATE(veg_falloc%labile%phosphorus(npoints))
    ALLOCATE(veg_falloc%labile%carbon13(npoints))
    ALLOCATE(veg_falloc%labile%carbon14(npoints))
    ALLOCATE(veg_falloc%labile%nitrogen15(npoints))
    ALLOCATE(veg_falloc%reserve%carbon(npoints))
    ALLOCATE(veg_falloc%reserve%nitrogen(npoints))
    ALLOCATE(veg_falloc%reserve%phosphorus(npoints))
    ALLOCATE(veg_falloc%reserve%carbon13(npoints))
    ALLOCATE(veg_falloc%reserve%carbon14(npoints))
    ALLOCATE(veg_falloc%reserve%nitrogen15(npoints))
    ALLOCATE(veg_falloc%fruit%carbon(npoints))
    ALLOCATE(veg_falloc%fruit%nitrogen(npoints))
    ALLOCATE(veg_falloc%fruit%phosphorus(npoints))
    ALLOCATE(veg_falloc%fruit%carbon13(npoints))
    ALLOCATE(veg_falloc%fruit%carbon14(npoints))
    ALLOCATE(veg_falloc%fruit%nitrogen15(npoints))
    ALLOCATE(veg_falloc%seed_bed%carbon(npoints))
    ALLOCATE(veg_falloc%seed_bed%nitrogen(npoints))
    ALLOCATE(veg_falloc%seed_bed%phosphorus(npoints))
    ALLOCATE(veg_falloc%seed_bed%carbon13(npoints))
    ALLOCATE(veg_falloc%seed_bed%carbon14(npoints))
    ALLOCATE(veg_falloc%seed_bed%nitrogen15(npoints))

    ALLOCATE(veg_total_flux%leaf%carbon(npoints))
    ALLOCATE(veg_total_flux%leaf%nitrogen(npoints))
    ALLOCATE(veg_total_flux%leaf%phosphorus(npoints))
    ALLOCATE(veg_total_flux%leaf%carbon13(npoints))
    ALLOCATE(veg_total_flux%leaf%carbon14(npoints))
    ALLOCATE(veg_total_flux%leaf%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%fine_root%carbon(npoints))
    ALLOCATE(veg_total_flux%fine_root%nitrogen(npoints))
    ALLOCATE(veg_total_flux%fine_root%phosphorus(npoints))
    ALLOCATE(veg_total_flux%fine_root%carbon13(npoints))
    ALLOCATE(veg_total_flux%fine_root%carbon14(npoints))
    ALLOCATE(veg_total_flux%fine_root%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%coarse_root%carbon(npoints))
    ALLOCATE(veg_total_flux%coarse_root%nitrogen(npoints))
    ALLOCATE(veg_total_flux%coarse_root%phosphorus(npoints))
    ALLOCATE(veg_total_flux%coarse_root%carbon13(npoints))
    ALLOCATE(veg_total_flux%coarse_root%carbon14(npoints))
    ALLOCATE(veg_total_flux%coarse_root%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%sap_wood%carbon(npoints))
    ALLOCATE(veg_total_flux%sap_wood%nitrogen(npoints))
    ALLOCATE(veg_total_flux%sap_wood%phosphorus(npoints))
    ALLOCATE(veg_total_flux%sap_wood%carbon13(npoints))
    ALLOCATE(veg_total_flux%sap_wood%carbon14(npoints))
    ALLOCATE(veg_total_flux%sap_wood%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%heart_wood%carbon(npoints))
    ALLOCATE(veg_total_flux%heart_wood%nitrogen(npoints))
    ALLOCATE(veg_total_flux%heart_wood%phosphorus(npoints))
    ALLOCATE(veg_total_flux%heart_wood%carbon13(npoints))
    ALLOCATE(veg_total_flux%heart_wood%carbon14(npoints))
    ALLOCATE(veg_total_flux%heart_wood%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%labile%carbon(npoints))
    ALLOCATE(veg_total_flux%labile%nitrogen(npoints))
    ALLOCATE(veg_total_flux%labile%phosphorus(npoints))
    ALLOCATE(veg_total_flux%labile%carbon13(npoints))
    ALLOCATE(veg_total_flux%labile%carbon14(npoints))
    ALLOCATE(veg_total_flux%labile%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%reserve%carbon(npoints))
    ALLOCATE(veg_total_flux%reserve%nitrogen(npoints))
    ALLOCATE(veg_total_flux%reserve%phosphorus(npoints))
    ALLOCATE(veg_total_flux%reserve%carbon13(npoints))
    ALLOCATE(veg_total_flux%reserve%carbon14(npoints))
    ALLOCATE(veg_total_flux%reserve%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%fruit%carbon(npoints))
    ALLOCATE(veg_total_flux%fruit%nitrogen(npoints))
    ALLOCATE(veg_total_flux%fruit%phosphorus(npoints))
    ALLOCATE(veg_total_flux%fruit%carbon13(npoints))
    ALLOCATE(veg_total_flux%fruit%carbon14(npoints))
    ALLOCATE(veg_total_flux%fruit%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%seed_bed%carbon(npoints))
    ALLOCATE(veg_total_flux%seed_bed%nitrogen(npoints))
    ALLOCATE(veg_total_flux%seed_bed%phosphorus(npoints))
    ALLOCATE(veg_total_flux%seed_bed%carbon13(npoints))
    ALLOCATE(veg_total_flux%seed_bed%carbon14(npoints))
    ALLOCATE(veg_total_flux%seed_bed%nitrogen15(npoints))

    ALLOCATE(veg_exudation%carbon(npoints))
    ALLOCATE(veg_exudation%nitrogen(npoints))
    ALLOCATE(veg_exudation%phosphorus(npoints))
    ALLOCATE(veg_exudation%carbon13(npoints))
    ALLOCATE(veg_exudation%carbon14(npoints))
    ALLOCATE(veg_exudation%nitrogen15(npoints))

    ALLOCATE(veg_establishment%carbon(npoints))
    ALLOCATE(veg_establishment%nitrogen(npoints))
    ALLOCATE(veg_establishment%phosphorus(npoints))
    ALLOCATE(veg_establishment%carbon13(npoints))
    ALLOCATE(veg_establishment%carbon14(npoints))
    ALLOCATE(veg_establishment%nitrogen15(npoints))

    ALLOCATE(veg_reserve_use%carbon(npoints))
    ALLOCATE(veg_reserve_use%nitrogen(npoints))
    ALLOCATE(veg_reserve_use%phosphorus(npoints))
    ALLOCATE(veg_reserve_use%carbon13(npoints))
    ALLOCATE(veg_reserve_use%carbon14(npoints))
    ALLOCATE(veg_reserve_use%nitrogen15(npoints))

    ALLOCATE(veg_total_biomass%carbon(npoints))
    ALLOCATE(veg_total_biomass%nitrogen(npoints))
    ALLOCATE(veg_total_biomass%phosphorus(npoints))
    ALLOCATE(veg_total_biomass%carbon13(npoints))
    ALLOCATE(veg_total_biomass%carbon14(npoints))
    ALLOCATE(veg_total_biomass%nitrogen15(npoints))

  END SUBROUTINE Init_veg_pools_elements_a


  !-----------------------------------------------------------------------------------------------------
  !> DEALLOCATE variables from veg_pool_a and veg_element_a structures at the program end
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Deallocate_var_veg_pools_elements_a

    IMPLICIT NONE

    !< 1.0  DEALLOCATE variables
    !! 
    DEALLOCATE(veg_pool%leaf%carbon)
    DEALLOCATE(veg_pool%leaf%nitrogen)
    DEALLOCATE(veg_pool%leaf%phosphorus)
    DEALLOCATE(veg_pool%leaf%carbon13)
    DEALLOCATE(veg_pool%leaf%carbon14)
    DEALLOCATE(veg_pool%leaf%nitrogen15)
    DEALLOCATE(veg_pool%fine_root%carbon)
    DEALLOCATE(veg_pool%fine_root%nitrogen)
    DEALLOCATE(veg_pool%fine_root%phosphorus)
    DEALLOCATE(veg_pool%fine_root%carbon13)
    DEALLOCATE(veg_pool%fine_root%carbon14)
    DEALLOCATE(veg_pool%fine_root%nitrogen15)
    DEALLOCATE(veg_pool%coarse_root%carbon)
    DEALLOCATE(veg_pool%coarse_root%nitrogen)
    DEALLOCATE(veg_pool%coarse_root%phosphorus)
    DEALLOCATE(veg_pool%coarse_root%carbon13)
    DEALLOCATE(veg_pool%coarse_root%carbon14)
    DEALLOCATE(veg_pool%coarse_root%nitrogen15)
    DEALLOCATE(veg_pool%sap_wood%carbon)
    DEALLOCATE(veg_pool%sap_wood%nitrogen)
    DEALLOCATE(veg_pool%sap_wood%phosphorus)
    DEALLOCATE(veg_pool%sap_wood%carbon13)
    DEALLOCATE(veg_pool%sap_wood%carbon14)
    DEALLOCATE(veg_pool%sap_wood%nitrogen15)
    DEALLOCATE(veg_pool%heart_wood%carbon)
    DEALLOCATE(veg_pool%heart_wood%nitrogen)
    DEALLOCATE(veg_pool%heart_wood%phosphorus)
    DEALLOCATE(veg_pool%heart_wood%carbon13)
    DEALLOCATE(veg_pool%heart_wood%carbon14)
    DEALLOCATE(veg_pool%heart_wood%nitrogen15)
    DEALLOCATE(veg_pool%labile%carbon)
    DEALLOCATE(veg_pool%labile%nitrogen)
    DEALLOCATE(veg_pool%labile%phosphorus)
    DEALLOCATE(veg_pool%labile%carbon13)
    DEALLOCATE(veg_pool%labile%carbon14)
    DEALLOCATE(veg_pool%labile%nitrogen15)
    DEALLOCATE(veg_pool%reserve%carbon)
    DEALLOCATE(veg_pool%reserve%nitrogen)
    DEALLOCATE(veg_pool%reserve%phosphorus)
    DEALLOCATE(veg_pool%reserve%carbon13)
    DEALLOCATE(veg_pool%reserve%carbon14)
    DEALLOCATE(veg_pool%reserve%nitrogen15)
    DEALLOCATE(veg_pool%fruit%carbon)
    DEALLOCATE(veg_pool%fruit%nitrogen)
    DEALLOCATE(veg_pool%fruit%phosphorus)
    DEALLOCATE(veg_pool%fruit%carbon13)
    DEALLOCATE(veg_pool%fruit%carbon14)
    DEALLOCATE(veg_pool%fruit%nitrogen15)
    DEALLOCATE(veg_pool%seed_bed%carbon)
    DEALLOCATE(veg_pool%seed_bed%nitrogen)
    DEALLOCATE(veg_pool%seed_bed%phosphorus)
    DEALLOCATE(veg_pool%seed_bed%carbon13)
    DEALLOCATE(veg_pool%seed_bed%carbon14)
    DEALLOCATE(veg_pool%seed_bed%nitrogen15)

    DEALLOCATE(veg_growth%leaf%carbon)
    DEALLOCATE(veg_growth%leaf%nitrogen)
    DEALLOCATE(veg_growth%leaf%phosphorus)
    DEALLOCATE(veg_growth%leaf%carbon13)
    DEALLOCATE(veg_growth%leaf%carbon14)
    DEALLOCATE(veg_growth%leaf%nitrogen15)
    DEALLOCATE(veg_growth%fine_root%carbon)
    DEALLOCATE(veg_growth%fine_root%nitrogen)
    DEALLOCATE(veg_growth%fine_root%phosphorus)
    DEALLOCATE(veg_growth%fine_root%carbon13)
    DEALLOCATE(veg_growth%fine_root%carbon14)
    DEALLOCATE(veg_growth%fine_root%nitrogen15)
    DEALLOCATE(veg_growth%coarse_root%carbon)
    DEALLOCATE(veg_growth%coarse_root%nitrogen)
    DEALLOCATE(veg_growth%coarse_root%phosphorus)
    DEALLOCATE(veg_growth%coarse_root%carbon13)
    DEALLOCATE(veg_growth%coarse_root%carbon14)
    DEALLOCATE(veg_growth%coarse_root%nitrogen15)
    DEALLOCATE(veg_growth%sap_wood%carbon)
    DEALLOCATE(veg_growth%sap_wood%nitrogen)
    DEALLOCATE(veg_growth%sap_wood%phosphorus)
    DEALLOCATE(veg_growth%sap_wood%carbon13)
    DEALLOCATE(veg_growth%sap_wood%carbon14)
    DEALLOCATE(veg_growth%sap_wood%nitrogen15)
    DEALLOCATE(veg_growth%heart_wood%carbon)
    DEALLOCATE(veg_growth%heart_wood%nitrogen)
    DEALLOCATE(veg_growth%heart_wood%phosphorus)
    DEALLOCATE(veg_growth%heart_wood%carbon13)
    DEALLOCATE(veg_growth%heart_wood%carbon14)
    DEALLOCATE(veg_growth%heart_wood%nitrogen15)
    DEALLOCATE(veg_growth%labile%carbon)
    DEALLOCATE(veg_growth%labile%nitrogen)
    DEALLOCATE(veg_growth%labile%phosphorus)
    DEALLOCATE(veg_growth%labile%carbon13)
    DEALLOCATE(veg_growth%labile%carbon14)
    DEALLOCATE(veg_growth%labile%nitrogen15)
    DEALLOCATE(veg_growth%reserve%carbon)
    DEALLOCATE(veg_growth%reserve%nitrogen)
    DEALLOCATE(veg_growth%reserve%phosphorus)
    DEALLOCATE(veg_growth%reserve%carbon13)
    DEALLOCATE(veg_growth%reserve%carbon14)
    DEALLOCATE(veg_growth%reserve%nitrogen15)
    DEALLOCATE(veg_growth%fruit%carbon)
    DEALLOCATE(veg_growth%fruit%nitrogen)
    DEALLOCATE(veg_growth%fruit%phosphorus)
    DEALLOCATE(veg_growth%fruit%carbon13)
    DEALLOCATE(veg_growth%fruit%carbon14)
    DEALLOCATE(veg_growth%fruit%nitrogen15)
    DEALLOCATE(veg_growth%seed_bed%carbon)
    DEALLOCATE(veg_growth%seed_bed%nitrogen)
    DEALLOCATE(veg_growth%seed_bed%phosphorus)
    DEALLOCATE(veg_growth%seed_bed%carbon13)
    DEALLOCATE(veg_growth%seed_bed%carbon14)
    DEALLOCATE(veg_growth%seed_bed%nitrogen15)

    DEALLOCATE(veg_litterfall%leaf%carbon)
    DEALLOCATE(veg_litterfall%leaf%nitrogen)
    DEALLOCATE(veg_litterfall%leaf%phosphorus)
    DEALLOCATE(veg_litterfall%leaf%carbon13)
    DEALLOCATE(veg_litterfall%leaf%carbon14)
    DEALLOCATE(veg_litterfall%leaf%nitrogen15)
    DEALLOCATE(veg_litterfall%fine_root%carbon)
    DEALLOCATE(veg_litterfall%fine_root%nitrogen)
    DEALLOCATE(veg_litterfall%fine_root%phosphorus)
    DEALLOCATE(veg_litterfall%fine_root%carbon13)
    DEALLOCATE(veg_litterfall%fine_root%carbon14)
    DEALLOCATE(veg_litterfall%fine_root%nitrogen15)
    DEALLOCATE(veg_litterfall%coarse_root%carbon)
    DEALLOCATE(veg_litterfall%coarse_root%nitrogen)
    DEALLOCATE(veg_litterfall%coarse_root%phosphorus)
    DEALLOCATE(veg_litterfall%coarse_root%carbon13)
    DEALLOCATE(veg_litterfall%coarse_root%carbon14)
    DEALLOCATE(veg_litterfall%coarse_root%nitrogen15)
    DEALLOCATE(veg_litterfall%sap_wood%carbon)
    DEALLOCATE(veg_litterfall%sap_wood%nitrogen)
    DEALLOCATE(veg_litterfall%sap_wood%phosphorus)
    DEALLOCATE(veg_litterfall%sap_wood%carbon13)
    DEALLOCATE(veg_litterfall%sap_wood%carbon14)
    DEALLOCATE(veg_litterfall%sap_wood%nitrogen15)
    DEALLOCATE(veg_litterfall%heart_wood%carbon)
    DEALLOCATE(veg_litterfall%heart_wood%nitrogen)
    DEALLOCATE(veg_litterfall%heart_wood%phosphorus)
    DEALLOCATE(veg_litterfall%heart_wood%carbon13)
    DEALLOCATE(veg_litterfall%heart_wood%carbon14)
    DEALLOCATE(veg_litterfall%heart_wood%nitrogen15)
    DEALLOCATE(veg_litterfall%labile%carbon)
    DEALLOCATE(veg_litterfall%labile%nitrogen)
    DEALLOCATE(veg_litterfall%labile%phosphorus)
    DEALLOCATE(veg_litterfall%labile%carbon13)
    DEALLOCATE(veg_litterfall%labile%carbon14)
    DEALLOCATE(veg_litterfall%labile%nitrogen15)
    DEALLOCATE(veg_litterfall%reserve%carbon)
    DEALLOCATE(veg_litterfall%reserve%nitrogen)
    DEALLOCATE(veg_litterfall%reserve%phosphorus)
    DEALLOCATE(veg_litterfall%reserve%carbon13)
    DEALLOCATE(veg_litterfall%reserve%carbon14)
    DEALLOCATE(veg_litterfall%reserve%nitrogen15)
    DEALLOCATE(veg_litterfall%fruit%carbon)
    DEALLOCATE(veg_litterfall%fruit%nitrogen)
    DEALLOCATE(veg_litterfall%fruit%phosphorus)
    DEALLOCATE(veg_litterfall%fruit%carbon13)
    DEALLOCATE(veg_litterfall%fruit%carbon14)
    DEALLOCATE(veg_litterfall%fruit%nitrogen15)
    DEALLOCATE(veg_litterfall%seed_bed%carbon)
    DEALLOCATE(veg_litterfall%seed_bed%nitrogen)
    DEALLOCATE(veg_litterfall%seed_bed%phosphorus)
    DEALLOCATE(veg_litterfall%seed_bed%carbon13)
    DEALLOCATE(veg_litterfall%seed_bed%carbon14)
    DEALLOCATE(veg_litterfall%seed_bed%nitrogen15)

    DEALLOCATE(veg_falloc%leaf%carbon)
    DEALLOCATE(veg_falloc%leaf%nitrogen)
    DEALLOCATE(veg_falloc%leaf%phosphorus)
    DEALLOCATE(veg_falloc%leaf%carbon13)
    DEALLOCATE(veg_falloc%leaf%carbon14)
    DEALLOCATE(veg_falloc%leaf%nitrogen15)
    DEALLOCATE(veg_falloc%fine_root%carbon)
    DEALLOCATE(veg_falloc%fine_root%nitrogen)
    DEALLOCATE(veg_falloc%fine_root%phosphorus)
    DEALLOCATE(veg_falloc%fine_root%carbon13)
    DEALLOCATE(veg_falloc%fine_root%carbon14)
    DEALLOCATE(veg_falloc%fine_root%nitrogen15)
    DEALLOCATE(veg_falloc%coarse_root%carbon)
    DEALLOCATE(veg_falloc%coarse_root%nitrogen)
    DEALLOCATE(veg_falloc%coarse_root%phosphorus)
    DEALLOCATE(veg_falloc%coarse_root%carbon13)
    DEALLOCATE(veg_falloc%coarse_root%carbon14)
    DEALLOCATE(veg_falloc%coarse_root%nitrogen15)
    DEALLOCATE(veg_falloc%sap_wood%carbon)
    DEALLOCATE(veg_falloc%sap_wood%nitrogen)
    DEALLOCATE(veg_falloc%sap_wood%phosphorus)
    DEALLOCATE(veg_falloc%sap_wood%carbon13)
    DEALLOCATE(veg_falloc%sap_wood%carbon14)
    DEALLOCATE(veg_falloc%sap_wood%nitrogen15)
    DEALLOCATE(veg_falloc%heart_wood%carbon)
    DEALLOCATE(veg_falloc%heart_wood%nitrogen)
    DEALLOCATE(veg_falloc%heart_wood%phosphorus)
    DEALLOCATE(veg_falloc%heart_wood%carbon13)
    DEALLOCATE(veg_falloc%heart_wood%carbon14)
    DEALLOCATE(veg_falloc%heart_wood%nitrogen15)
    DEALLOCATE(veg_falloc%labile%carbon)
    DEALLOCATE(veg_falloc%labile%nitrogen)
    DEALLOCATE(veg_falloc%labile%phosphorus)
    DEALLOCATE(veg_falloc%labile%carbon13)
    DEALLOCATE(veg_falloc%labile%carbon14)
    DEALLOCATE(veg_falloc%labile%nitrogen15)
    DEALLOCATE(veg_falloc%reserve%carbon)
    DEALLOCATE(veg_falloc%reserve%nitrogen)
    DEALLOCATE(veg_falloc%reserve%phosphorus)
    DEALLOCATE(veg_falloc%reserve%carbon13)
    DEALLOCATE(veg_falloc%reserve%carbon14)
    DEALLOCATE(veg_falloc%reserve%nitrogen15)
    DEALLOCATE(veg_falloc%fruit%carbon)
    DEALLOCATE(veg_falloc%fruit%nitrogen)
    DEALLOCATE(veg_falloc%fruit%phosphorus)
    DEALLOCATE(veg_falloc%fruit%carbon13)
    DEALLOCATE(veg_falloc%fruit%carbon14)
    DEALLOCATE(veg_falloc%fruit%nitrogen15)
    DEALLOCATE(veg_falloc%seed_bed%carbon)
    DEALLOCATE(veg_falloc%seed_bed%nitrogen)
    DEALLOCATE(veg_falloc%seed_bed%phosphorus)
    DEALLOCATE(veg_falloc%seed_bed%carbon13)
    DEALLOCATE(veg_falloc%seed_bed%carbon14)
    DEALLOCATE(veg_falloc%seed_bed%nitrogen15)

    DEALLOCATE(veg_total_flux%leaf%carbon)
    DEALLOCATE(veg_total_flux%leaf%nitrogen)
    DEALLOCATE(veg_total_flux%leaf%phosphorus)
    DEALLOCATE(veg_total_flux%leaf%carbon13)
    DEALLOCATE(veg_total_flux%leaf%carbon14)
    DEALLOCATE(veg_total_flux%leaf%nitrogen15)
    DEALLOCATE(veg_total_flux%fine_root%carbon)
    DEALLOCATE(veg_total_flux%fine_root%nitrogen)
    DEALLOCATE(veg_total_flux%fine_root%phosphorus)
    DEALLOCATE(veg_total_flux%fine_root%carbon13)
    DEALLOCATE(veg_total_flux%fine_root%carbon14)
    DEALLOCATE(veg_total_flux%fine_root%nitrogen15)
    DEALLOCATE(veg_total_flux%coarse_root%carbon)
    DEALLOCATE(veg_total_flux%coarse_root%nitrogen)
    DEALLOCATE(veg_total_flux%coarse_root%phosphorus)
    DEALLOCATE(veg_total_flux%coarse_root%carbon13)
    DEALLOCATE(veg_total_flux%coarse_root%carbon14)
    DEALLOCATE(veg_total_flux%coarse_root%nitrogen15)
    DEALLOCATE(veg_total_flux%sap_wood%carbon)
    DEALLOCATE(veg_total_flux%sap_wood%nitrogen)
    DEALLOCATE(veg_total_flux%sap_wood%phosphorus)
    DEALLOCATE(veg_total_flux%sap_wood%carbon13)
    DEALLOCATE(veg_total_flux%sap_wood%carbon14)
    DEALLOCATE(veg_total_flux%sap_wood%nitrogen15)
    DEALLOCATE(veg_total_flux%heart_wood%carbon)
    DEALLOCATE(veg_total_flux%heart_wood%nitrogen)
    DEALLOCATE(veg_total_flux%heart_wood%phosphorus)
    DEALLOCATE(veg_total_flux%heart_wood%carbon13)
    DEALLOCATE(veg_total_flux%heart_wood%carbon14)
    DEALLOCATE(veg_total_flux%heart_wood%nitrogen15)
    DEALLOCATE(veg_total_flux%labile%carbon)
    DEALLOCATE(veg_total_flux%labile%nitrogen)
    DEALLOCATE(veg_total_flux%labile%phosphorus)
    DEALLOCATE(veg_total_flux%labile%carbon13)
    DEALLOCATE(veg_total_flux%labile%carbon14)
    DEALLOCATE(veg_total_flux%labile%nitrogen15)
    DEALLOCATE(veg_total_flux%reserve%carbon)
    DEALLOCATE(veg_total_flux%reserve%nitrogen)
    DEALLOCATE(veg_total_flux%reserve%phosphorus)
    DEALLOCATE(veg_total_flux%reserve%carbon13)
    DEALLOCATE(veg_total_flux%reserve%carbon14)
    DEALLOCATE(veg_total_flux%reserve%nitrogen15)
    DEALLOCATE(veg_total_flux%fruit%carbon)
    DEALLOCATE(veg_total_flux%fruit%nitrogen)
    DEALLOCATE(veg_total_flux%fruit%phosphorus)
    DEALLOCATE(veg_total_flux%fruit%carbon13)
    DEALLOCATE(veg_total_flux%fruit%carbon14)
    DEALLOCATE(veg_total_flux%fruit%nitrogen15)
    DEALLOCATE(veg_total_flux%seed_bed%carbon)
    DEALLOCATE(veg_total_flux%seed_bed%nitrogen)
    DEALLOCATE(veg_total_flux%seed_bed%phosphorus)
    DEALLOCATE(veg_total_flux%seed_bed%carbon13)
    DEALLOCATE(veg_total_flux%seed_bed%carbon14)
    DEALLOCATE(veg_total_flux%seed_bed%nitrogen15)

    DEALLOCATE(veg_exudation%carbon)
    DEALLOCATE(veg_exudation%nitrogen)
    DEALLOCATE(veg_exudation%phosphorus)
    DEALLOCATE(veg_exudation%carbon13)
    DEALLOCATE(veg_exudation%carbon14)
    DEALLOCATE(veg_exudation%nitrogen15)

    DEALLOCATE(veg_establishment%carbon)
    DEALLOCATE(veg_establishment%nitrogen)
    DEALLOCATE(veg_establishment%phosphorus)
    DEALLOCATE(veg_establishment%carbon13)
    DEALLOCATE(veg_establishment%carbon14)
    DEALLOCATE(veg_establishment%nitrogen15)

    DEALLOCATE(veg_reserve_use%carbon)
    DEALLOCATE(veg_reserve_use%nitrogen)
    DEALLOCATE(veg_reserve_use%phosphorus)
    DEALLOCATE(veg_reserve_use%carbon13)
    DEALLOCATE(veg_reserve_use%carbon14)
    DEALLOCATE(veg_reserve_use%nitrogen15)

    DEALLOCATE(veg_total_biomass%carbon)
    DEALLOCATE(veg_total_biomass%nitrogen)
    DEALLOCATE(veg_total_biomass%phosphorus)
    DEALLOCATE(veg_total_biomass%carbon13)
    DEALLOCATE(veg_total_biomass%carbon14)
    DEALLOCATE(veg_total_biomass%nitrogen15)

  END SUBROUTINE Deallocate_var_veg_pools_elements_a

END MODULE mo_veg_memory_class
