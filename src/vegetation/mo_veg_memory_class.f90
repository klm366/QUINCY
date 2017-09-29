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

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d         , &
                                       t_jsb_struct, t_jsb_n_struct                                   ! *slmjsb4

  USE mo_slm_varlist,            ONLY: jsb_add_var   ! replace by mo_jsb_varlist later    *slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_veg_memory, max_no_of_vars
  PUBLIC :: t_slm_struct_veg_pool, t_slm_struct_veg_element   ! *slmjsb4 
  PUBLIC :: max_no_of_structs, max_no_of_n_structs            ! *slmjsb4
  PUBLIC :: t_veg_element_a, t_veg_pool_a                     ! *slmtemp

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
    PROCEDURE :: init => init_t_slm_struct_veg_element
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
    PROCEDURE :: init => init_t_slm_struct_veg_pool
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

    ! basic plant
    TYPE(t_jsb_var_real2d)          :: &   
                                    height            , &            !< plant height in m
                                    diameter          , &            !< plant diameter in m
                                    dens_ind          , &            !< density of individuals (#/m2) 
                                    lai               , &            !< actual leaf area index
                                    target_cn_leaf    , &            !< target foliar C:N ratio
                                    target_np_leaf    , &            !< target foliar N:P ratio 
                                    mean_leaf_age     , &            !< average foliage age in days 
                                    f_n_demand        , &            !< relative plant demand for N; see update_vegsoil_bgc
                                    f_p_demand        , &            !< relative plant demand for P; see update_vegsoil_bgc
                                    beta_sinklim_ps   , &            !< scaling factor to account for a direct sink limitation contraint on photosynthesis
                                    dphi              , &            !< plant-soil water potential gradient (MPa)
                                    growth_req_n      , &            !< moles N required for a unit of C growth under current allocation fractions 
                                    growth_req_p      , &            !< moles P required for a unit of N growth under current allocation fractions 
                                    leaf_cn_direction , &            !< indicates direction of change in leaf CN ratio to maximise NPP
                                    target_lai        , &            !< LAI implied from sap wood area
                                    kstar_labile                     !< labile pool turnover rate assuming only temperature and moisture constraints (no stoichiometry)

    ! target C:N and N:P ratios for various non-leaf tissues (mol/mol)
    TYPE(t_jsb_var_real2d)          :: &
                                    target_cn_fine_root       , &
                                    target_cn_coarse_root     , &
                                    target_cn_sap_wood        , &
                                    target_cn_fruit           , &
                                    target_np_fine_root       , &
                                    target_np_coarse_root     , &
                                    target_np_sap_wood        , &
                                    target_np_fruit

    ! canopy layers
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    leaf_nitrogen_cl         , &       !< N in canopy layer (mmol / m2 (LAI))
                                    fn_chl_cl        , &       !< fraction of N in chlorophyll (pigments only!)
                                    fn_et_cl         , &       !< fraction of N in electron transport 
                                    fn_rub_cl        , &       !< fraction of N in Rubisco 
                                    fn_pepc_cl       , &       !< fraction of N in PEP carboylase (C4 plants only)
                                    fn_oth_cl        , &       !< fraction of N not associated with photosynthesis 
                                    fleaf_sunlit_cl   , &       !< fraction of layer which is sunlit
                                    lai_cl            , &       !< leaf area in the canopy layer 
                                    cumm_lai_cl                 !< cummulative leaf area above the mid-point of the canopy layer 

    ! fluxes
    TYPE(t_jsb_var_real2d)          :: &   
                                    maint_respiration_pot     , & !< potential maintenance respiration (µmol CO2 m-2 s-1)
                                    maint_respiration         , & !< maintenance respiration (µmol CO2 m-2 s-1)
                                    maint_respiration_c13     , & !< maintenance respiration (µmol 13CO2 m-2 s-1)
                                    maint_respiration_c14     , & !< maintenance respiration (µmol 14CO2 m-2 s-1)
                                    growth_respiration        , & !< growth respiration (µmol CO2 m-2 s-1)
                                    growth_respiration_c13    , & !< growth respiration (µmol 13CO2 m-2 s-1)
                                    growth_respiration_c14    , & !< growth respiration (µmol 14CO2 m-2 s-1)
                                    n_transform_respiration   , & !< respiration associated with N transformation (µmol CO2 m-2 s-1)
                                    n_fixation_respiration    , & !< respiration associated with N fixation (µmol CO2 m-2 s-1)
                                    n_processing_respiration  , &       !< respiration associated with N processing (µmol CO2 m-2 s-1)
                                    n_processing_respiration_c13 , &  !< respiration associated with N processing (µmol 13CO2 m-2 s-1)
                                    n_processing_respiration_c14 , &  !< respiration associated with N processing (µmol 14CO2 m-2 s-1)
                                    npp                       , & !< net primary production (µmol CO2 m-2 s-1)
                                    npp_c13                   , & !< net primary production (µmol 13CO2 m-2 s-1)
                                    npp_c14                   , & !< net primary production (µmol 14CO2 m-2 s-1)
                                    uptake_nh4                , & !< NH4 uptake from soil (µmol N m-2 s-1) 
                                    uptake_nh4_n15            , & !< 15NH4 uptake from soil (µmol N m-2 s-1)
                                    uptake_no3                , & !< NO3 uptake from soil (µmol N m-2 s-1)
                                    uptake_no3_n15            , & !< 15NO3 uptake from soil (µmol N m-2 s-1)
                                    uptake_po4                , & !< PO4 uptake from soil (µmol P m-2 s-1) 
                                    n_fixation                , & !< symbiontic N fixation (µmol N m-2 s-1)
                                    n_fixation_n15            , & !< symbiontic N fixation (µmol 15N m-2 s-1)
                                    recycling_leaf_n          , & !< net flux of N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
                                    recycling_leaf_n15        , & !< net flux of 15N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1) 
                                    recycling_leaf_p          , & !< net flux of P from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
                                    recycling_fine_root_n     , & !< net flux of N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
                                    recycling_fine_root_n15   , & !< net flux of 15N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1) 
                                    recycling_fine_root_p     , & !< net flux of P from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
                                    recycling_heart_wood_n      , & !< net flux of N from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1)
                                    recycling_heart_wood_n15    , & !< net flux of 15N from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1) 
                                    recycling_heart_wood_p      , & !< net flux of P from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1)
                                    transpiration                 !< transpiration in kg m-2 s-1

    ! everything else...
    TYPE(t_jsb_var_real2d)          :: &   
                                    fmaint_rate               , & !< current N-specific maintentance respiration (~µmol C/mol N / s ) 
                                    unit_npp                  , & !< incremental change in NPP. Needed for optimal allocation
                                    unit_transpiration        , & !< incremental change in transpiration. Needed for optimal allocation
                                    unit_uptake_n_pot         , & !< incremental change in potential N uptake per unit root carbon, needed for N uptake calculation 
                                    unit_uptake_p_pot         , & !< incremental change in potential P uptake per unit root carbon, currently unused 
                                    unit_uptake_n_act         , & !< incremental change in actual N uptake oer unit root carbon, needed for optimal allocation
                                    unit_uptake_p_act         , & !< incremental change in actual P uptake per unit root carbon, needed for optimal allocation
                                    delta_dens_ind            , & !< change in density of individuals (#/m2/timestep)
                                    cost_n_uptake_root        , & !< actual cost for N uptake by roots (mol C/mol N)
                                    cost_n_fixation           , & !< actual cost for N uptake through fixation (mol C/mol N)
                                    mortality_rate            , & !< background_mort_rate_tree / grass converted to 1/timestep
                                    k1_opt                    , & !< coefficient for optimal leaf:root ratio
                                    k2_opt                    , & !< coefficient for optimal leaf:root ratio
                                    leaf2sapwood_mass_ratio   , & !< leaf to shoot mass ratio
                                    leaf2root_mass_ratio          !< leaf to root mass ratio

    ! daytime averages
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_daytime_cl          !<  average of the previous day
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_daytime             , &    !<  average of the previous day
                                    press_srf_daytime         , &    !<  average of the previous day
                                    co2_mixing_ratio_daytime  , &    !<  average of the previous day
                                    ga_daytime                , &    !<  average of the previous day
                                    beta_sinklim_ps_daytime   , &    !<  average of the previous day
                                    t_jmax_opt_daytime               !<  average of the previous day
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_daytime_dacc_cl     !<  value accumulated over current day
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_daytime_dacc            , & !<  average of the previous day
                                    press_srf_daytime_dacc        , & !<  average of the previous day
                                    co2_mixing_ratio_daytime_dacc , & !<  average of the previous day
                                    ga_daytime_dacc               , & !<  average of the previous day
                                    beta_sinklim_ps_daytime_dacc  , & !<  average of the previous day
                                    t_jmax_opt_daytime_dacc           !<  value accumulated over current day 
    TYPE(t_jsb_var_real2d)          :: &   
                                    daytime_counter  !< number of timesteps of current day with light available & daytime_avrg variables were added up

    ! moving averages
    ! 1.1 memory of conditions for the labile pool dynamics
    TYPE(t_jsb_var_real2d)          :: &   
                                    gpp_tlabile_mavg               , &  !< gpp  (µmol m-2 s-1)
                                    maint_respiration_tlabile_mavg , &  !< maintenance respiration (µmol m-2 s-1)
                                    growth_req_n_tlabile_mavg      , &  !< moles N required for a unit of C growth under current allocation fractions
                                    growth_req_p_tlabile_mavg           !< moles P required for a unit of N growth under current allocation fractions
    ! 1.2 memory of conditions for phenology calculations
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tphen_mavg                    !< air temp mavg for phenology calc
    ! 1.3 memory for plant uptake (average demand over a couple of days to get rid of diurnal C cycle)
    TYPE(t_jsb_var_real2d)          :: &   
                                    npp_tuptake_mavg               , &  !< npp of the last week (µmol m-2 s-1)
                                    demand_uptake_n_tuptake_mavg   , &  !< scalar of Vmax for root N uptake 
                                    demand_uptake_p_tuptake_mavg   , &  !< scalar of Vmax for root P uptake
                                    growth_req_n_tuptake_mavg      , &  !< moles N required for a unit of C growth under current allocation fractions
                                    growth_req_p_tuptake_mavg      , &  !< moles P required for a unit of N growth under current allocation fractions
                                    cost_n_fixation_tuptake_mavg        !< moving average of N fixation cost (mol C / mol N)
    ! 1.4 memory of conditions for photosynthesis for the optimisation of foliar N fractions 
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_tfrac_mavg_cl          !< fraction of sunlit leaves
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tfrac_mavg               , &  !< air temp
                                    press_srf_tfrac_mavg           , &  !< air pressure at surface
                                    co2_mixing_ratio_tfrac_mavg    , &  !< co2 mixing ratio of the air
                                    ga_tfrac_mavg                  , &  !< aerodynamic resistance
                                    beta_sinklim_ps_tfrac_mavg          !< scaling factor to account for soil moisture constraints on photosynthesis 
    ! 1.5 memory of conditions for photosynthesis for the optimisation of foliar N  
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_tcnl_mavg_cl           !< fraction of sunlit leaves
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tcnl_mavg                , &  !< air temp
                                    press_srf_tcnl_mavg            , &  !< air pressure at surface
                                    co2_mixing_ratio_tcnl_mavg     , &  !< co2 mixing ratio of the air
                                    ga_tcnl_mavg                   , &  !< aerodynamic resistance
                                    beta_sinklim_ps_tcnl_mavg      , &     !< scaling factor to account for soil moisture constraints on conductance
                                    uptake_n_tcnl_mavg             , &  !< root N uptake
                                    growth_cn_tcnl_mavg            , &  !< new growth CN ratio
                                    npp_tcnl_mavg                       !< NPP

    ! 1.6 memory of conditions for photosynthesis for the optimisation of biomass allocation fractions
    TYPE(t_jsb_var_real2d)          :: &   
                                    npp_talloc_mavg                , &  !< NPP (µmol m-2 s-1) (averaging time-scale: tau_allocation)
                                    n_fixation_talloc_mavg         , &  !< N fixation averaged at the allocation timescale
                                    transpiration_talloc_mavg      , &  !< average transpiration at allocation time-scale
                                    unit_npp_talloc_mavg           , &  !< rate of change in NPP with  increment in leaf mass, mumol C mol-1
                                    unit_transpiration_talloc_mavg , &  !< rate of tranpiration change per unit LAI change?
                                    unit_uptake_n_talloc_mavg  , &  !< actual rate of N uptkae per unit root, mumol N mol-1 (averaging time-scale: tau_allocation)
                                    unit_uptake_p_talloc_mavg  , &  !< actual rate of P uptake per unit root, mumol P mol-1 (time-scale of accumulation: tau_allocation)
                                    dphi_talloc_mavg               , &  !< leaf-soil water potential gradient (MPa)  (averaging time-scale: tau_allocation)
                                    growth_cn_talloc_mavg          , &  !< plant growth cn ratio, equivalent to 1/growth_req_n
                                    growth_cp_talloc_mavg          , &  !< plant growth cp ratio, equivalent to 1/(growthreq_n*growth_req_p)
                                    labile_carbon_talloc_mavg      , &  !< moving average of labile carbon pool neeeded for allometry calculations
                                    labile_nitrogen_talloc_mavg    , &  !< moving average of labile nitrogen pool neeeded for allometry calculations
                                    labile_phosphorus_talloc_mavg   ,&  !< moving average of labile phosphorus pool neeeded for allometry calculations
                                    reserve_carbon_talloc_mavg  ,&  !< moving average of reserve carbon needed for optimal allocation
                                    reserve_nitrogen_talloc_mavg,&  !< moving average of reserve nitrogen needed for optimal allocation
                                    reserve_phosphorus_talloc_mavg  !< moving average of reserve phosphorus needed for optimal allocation
                                    
                                    
    ! 1.7 memory over root life time 
    TYPE(t_jsb_var_real2d)          :: &   
                                    leaf2root_troot_mavg           , &  !< leaf to root mass ratio
                                    unit_npp_troot_mavg            , &  !< rate of change in NPP with  increment in leaf mass, mumol C mol-1
                                    unit_uptake_n_pot_troot_mavg   , &  !< rate of N uptkae per unit root, mumol N mol-1 (time-scale of accumulation: tau_root)
                                    unit_uptake_p_pot_troot_mavg   , &  !< rate of P uptake per unit root, mumol P mol-1 (time-scale of accumulation: tau_root)
                                    fmaint_rate_troot_mavg              !< moving average of N-specific maintentance respiration (~µmol C/mol N / s ) 
    ! 1.8 memory over timescale of vegetation dynamics 
    TYPE(t_jsb_var_real2d)          :: &   
                                    an_boc_tvegdyn_mavg                 !< net photosynthesis at bottom of canopy (lowest canopy layer existing)
    ! 1.9 memory of variables for respiration accclimation 
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tacclim_mavg                  !< average temperature for respiration acclimation
    ! 1.10 memory of certain variables, temperature, ...
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_week_mavg                , &  !< avrg air temperature across 7 days
                                    t_air_month_mavg               , &  !< avrg air temperature across 30 days
                                    t_jmax_opt_mavg                     !< optimim temperature for electron transport of photosynthesis

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
  PURE ELEMENTAL SUBROUTINE init_t_slm_struct_veg_element(this,s)
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
  END SUBROUTINE init_t_slm_struct_veg_element 

  !-----------------------------------------------------------------------------------------------------
  !> Initialisation t_slm_struct_veg_pool with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE init_t_slm_struct_veg_pool(this,s)
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
  END SUBROUTINE init_t_slm_struct_veg_pool 



  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the processes in: vegetation
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_veg_memory(mem, prefix, suffix, lct_ids, model_id)

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits
!    USE mo_jsb_io,            ONLY: TSTEP_CONSTANT, tables
    USE mo_jsb_io,            ONLY: t_cf, t_grib1, t_grib2
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,       ONLY: t_jsb_model

    ! 0.1.InOut
    CLASS(t_veg_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids !(:)
    INTEGER,             INTENT(in)    :: model_id
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
!    TYPE(t_jsb_vgrid), POINTER :: canopy_layer                 ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_veg_memory'

    model => Get_model(model_id)

!    table = tables(1)
    table = 1

!    hgrid   => Get_grid(mem%grid_id)
!    surface => Get_vgrid('surface')
!    canopy_layer => Get_vgrid('canopy_layer')

!    IF (.NOT. ASSOCIATED(mem%height%ptr)) THEN
    IF (.NOT. ALLOCATED(mem%height%ptr)) THEN

      !> add variables to veg memory veg\%VAR
      !! 
      !! jsb_add_var -> jsb_add_var_list_element_r1d    *slmtemp 1D -> 2D
      !!                jsb_add_var_list_element_l1d
      !!  
      CALL jsb_add_var(mem, 'height', mem%height, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'diameter', mem%diameter, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'dens_ind', mem%dens_ind, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'lai', mem%lai, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_cn_leaf', mem%target_cn_leaf, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_np_leaf', mem%target_np_leaf, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'mean_leaf_age', mem%mean_leaf_age, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'f_n_demand', mem%f_n_demand, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'f_p_demand', mem%f_p_demand, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_sinklim_ps', mem%beta_sinklim_ps, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'dphi', mem%dphi, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_req_n', mem%growth_req_n, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_req_p', mem%growth_req_p, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'leaf_cn_direction', mem%leaf_cn_direction, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_lai', mem%target_lai, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'kstar_labile', mem%kstar_labile, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_cn_fine_root', mem%target_cn_fine_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_cn_coarse_root', mem%target_cn_coarse_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_cn_sap_wood', mem%target_cn_sap_wood, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_cn_fruit', mem%target_cn_fruit, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_np_fine_root', mem%target_np_fine_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_np_coarse_root', mem%target_np_coarse_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_np_sap_wood', mem%target_np_sap_wood, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'target_np_fruit', mem%target_np_fruit, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'leaf_nitrogen_cl', mem%leaf_nitrogen_cl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fn_chl_cl', mem%fn_chl_cl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fn_et_cl', mem%fn_et_cl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fn_rub_cl', mem%fn_rub_cl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fn_pepc_cl', mem%fn_pepc_cl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fn_oth_cl', mem%fn_oth_cl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fleaf_sunlit_cl', mem%fleaf_sunlit_cl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'lai_cl', mem%lai_cl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'cumm_lai_cl', mem%cumm_lai_cl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'maint_respiration_pot', mem%maint_respiration_pot, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'maint_respiration', mem%maint_respiration, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'maint_respiration_c13', mem%maint_respiration_c13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'maint_respiration_c14', mem%maint_respiration_c14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_respiration', mem%growth_respiration, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_respiration_c13', mem%growth_respiration_c13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_respiration_c14', mem%growth_respiration_c14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n_transform_respiration', mem%n_transform_respiration, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n_fixation_respiration', mem%n_fixation_respiration, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n_processing_respiration', mem%n_processing_respiration, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n_processing_respiration_c13', mem%n_processing_respiration_c13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n_processing_respiration_c14', mem%n_processing_respiration_c14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'npp', mem%npp, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'npp_c13', mem%npp_c13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'npp_c14', mem%npp_c14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'uptake_nh4', mem%uptake_nh4, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'uptake_nh4_n15', mem%uptake_nh4_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'uptake_no3', mem%uptake_no3, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'uptake_no3_n15', mem%uptake_no3_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'uptake_po4', mem%uptake_po4, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n_fixation', mem%n_fixation, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n_fixation_n15', mem%n_fixation_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_leaf_n', mem%recycling_leaf_n, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_leaf_n15', mem%recycling_leaf_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_leaf_p', mem%recycling_leaf_p, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_fine_root_n', mem%recycling_fine_root_n, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_fine_root_n15', mem%recycling_fine_root_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_fine_root_p', mem%recycling_fine_root_p, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_heart_wood_n', mem%recycling_heart_wood_n, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_heart_wood_n15', mem%recycling_heart_wood_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'recycling_heart_wood_p', mem%recycling_heart_wood_p, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transpiration', mem%transpiration, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fmaint_rate', mem%fmaint_rate, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_npp', mem%unit_npp, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_transpiration', mem%unit_transpiration, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_n_pot', mem%unit_uptake_n_pot, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_p_pot', mem%unit_uptake_p_pot, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_n_act', mem%unit_uptake_n_act, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_p_act', mem%unit_uptake_p_act, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'delta_dens_ind', mem%delta_dens_ind, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'cost_n_uptake_root', mem%cost_n_uptake_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'cost_n_fixation', mem%cost_n_fixation, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'mortality_rate', mem%mortality_rate, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'k1_opt', mem%k1_opt, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'k2_opt', mem%k2_opt, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'leaf2sapwood_mass_ratio', mem%leaf2sapwood_mass_ratio, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'leaf2root_mass_ratio', mem%leaf2root_mass_ratio, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fleaf_sunlit_daytime_cl', mem%fleaf_sunlit_daytime_cl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_daytime', mem%t_air_daytime, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'press_srf_daytime', mem%press_srf_daytime, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'co2_mixing_ratio_daytime', mem%co2_mixing_ratio_daytime, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ga_daytime', mem%ga_daytime, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_sinklim_ps_daytime', mem%beta_sinklim_ps_daytime, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_jmax_opt_daytime', mem%t_jmax_opt_daytime, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fleaf_sunlit_daytime_dacc_cl', mem%fleaf_sunlit_daytime_dacc_cl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_daytime_dacc', mem%t_air_daytime_dacc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'press_srf_daytime_dacc', mem%press_srf_daytime_dacc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'co2_mixing_ratio_daytime_dacc', mem%co2_mixing_ratio_daytime_dacc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ga_daytime_dacc', mem%ga_daytime_dacc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_sinklim_ps_daytime_dacc', mem%beta_sinklim_ps_daytime_dacc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_jmax_opt_daytime_dacc', mem%t_jmax_opt_daytime_dacc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'daytime_counter', mem%daytime_counter, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'gpp_tlabile_mavg', mem%gpp_tlabile_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'maint_respiration_tlabile_mavg', mem%maint_respiration_tlabile_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_req_n_tlabile_mavg', mem%growth_req_n_tlabile_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_req_p_tlabile_mavg', mem%growth_req_p_tlabile_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_tphen_mavg', mem%t_air_tphen_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'npp_tuptake_mavg', mem%npp_tuptake_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'demand_uptake_n_tuptake_mavg', mem%demand_uptake_n_tuptake_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'demand_uptake_p_tuptake_mavg', mem%demand_uptake_p_tuptake_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_req_n_tuptake_mavg', mem%growth_req_n_tuptake_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_req_p_tuptake_mavg', mem%growth_req_p_tuptake_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'cost_n_fixation_tuptake_mavg', mem%cost_n_fixation_tuptake_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fleaf_sunlit_tfrac_mavg_cl', mem%fleaf_sunlit_tfrac_mavg_cl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_tfrac_mavg', mem%t_air_tfrac_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'press_srf_tfrac_mavg', mem%press_srf_tfrac_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'co2_mixing_ratio_tfrac_mavg', mem%co2_mixing_ratio_tfrac_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ga_tfrac_mavg', mem%ga_tfrac_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_sinklim_ps_tfrac_mavg', mem%beta_sinklim_ps_tfrac_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fleaf_sunlit_tcnl_mavg_cl', mem%fleaf_sunlit_tcnl_mavg_cl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_tcnl_mavg', mem%t_air_tcnl_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'press_srf_tcnl_mavg', mem%press_srf_tcnl_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'co2_mixing_ratio_tcnl_mavg', mem%co2_mixing_ratio_tcnl_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ga_tcnl_mavg', mem%ga_tcnl_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_sinklim_ps_tcnl_mavg', mem%beta_sinklim_ps_tcnl_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'uptake_n_tcnl_mavg', mem%uptake_n_tcnl_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'growth_cn_tcnl_mavg', mem%growth_cn_tcnl_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'npp_tcnl_mavg', mem%npp_tcnl_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'npp_talloc_mavg', mem%npp_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'n_fixation_talloc_mavg', mem%n_fixation_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transpiration_talloc_mavg', mem%transpiration_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_npp_talloc_mavg', mem%unit_npp_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_transpiration_talloc_mavg', mem%unit_transpiration_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_n_talloc_mavg', mem%unit_uptake_n_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_p_talloc_mavg', mem%unit_uptake_p_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'dphi_talloc_mavg', mem%dphi_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_cn_talloc_mavg', mem%growth_cn_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'growth_cp_talloc_mavg', mem%growth_cp_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'labile_carbon_talloc_mavg', mem%labile_carbon_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'labile_nitrogen_talloc_mavg', mem%labile_nitrogen_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'labile_phosphorus_talloc_mavg', mem%labile_phosphorus_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'reserve_carbon_talloc_mavg', mem%reserve_carbon_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'reserve_nitrogen_talloc_mavg', mem%reserve_nitrogen_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'reserve_phosphorus_talloc_mavg', mem%reserve_phosphorus_talloc_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'leaf2root_troot_mavg', mem%leaf2root_troot_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_npp_troot_mavg', mem%unit_npp_troot_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_n_pot_troot_mavg', mem%unit_uptake_n_pot_troot_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'unit_uptake_p_pot_troot_mavg', mem%unit_uptake_p_pot_troot_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'fmaint_rate_troot_mavg', mem%fmaint_rate_troot_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'an_boc_tvegdyn_mavg', mem%an_boc_tvegdyn_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_tacclim_mavg', mem%t_air_tacclim_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_week_mavg', mem%t_air_week_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_air_month_mavg', mem%t_air_month_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_jmax_opt_mavg', mem%t_jmax_opt_mavg, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      !> simple TYPE (t_jsb_struct)
      !! 
      !! jsb_add_var -> slm_add_struct()
      !! 
      CALL jsb_add_var(mem, 'veg_exudation',              mem%exudation, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'veg_establishment',          mem%establishment, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'veg_reserve_use',            mem%reserve_use, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'veg_total_biomass',          mem%total_biomass, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)

      !> nested TYPE (t_jsb_n_struct)
      !! 
      !! jsb_add_var -> slm_add_n_struct()
      !! 
      CALL jsb_add_var(mem, 'veg_pool'        , mem%pool, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'veg_growth'      , mem%growth, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'veg_litterfall'  , mem%litterfall, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'veg_falloc'      , mem%falloc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'veg_total_flux'  , mem%total_flux, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)     ! mass conservation diagnosis
      !> level 2 TYPE (l2_struct) nested within the above nested TYPE structures (t_jsb_n_struct)
      !! 
      !! jsb_add_var -> slm_add_l2_struct()
      !! 
      CALL jsb_add_var(mem%pool, 'veg_pool_leaf',         mem%pool%leaf, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)       
      CALL jsb_add_var(mem%pool, 'veg_pool_fine_root',    mem%pool%fine_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%pool, 'veg_pool_coarse_root',  mem%pool%coarse_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%pool, 'veg_pool_sap_wood',     mem%pool%sap_wood, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%pool, 'veg_pool_heart_wood',   mem%pool%heart_wood, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%pool, 'veg_pool_labile',       mem%pool%labile, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%pool, 'veg_pool_reserve',      mem%pool%reserve, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%pool, 'veg_pool_fruit',        mem%pool%fruit, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%pool, 'veg_pool_seed_bed',     mem%pool%seed_bed, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%growth, 'veg_growth_leaf',         mem%growth%leaf, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)       
      CALL jsb_add_var(mem%growth, 'veg_growth_fine_root',    mem%growth%fine_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%growth, 'veg_growth_coarse_root',  mem%growth%coarse_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%growth, 'veg_growth_sap_wood',     mem%growth%sap_wood, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%growth, 'veg_growth_heart_wood',   mem%growth%heart_wood, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%growth, 'veg_growth_labile',       mem%growth%labile, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%growth, 'veg_growth_reserve',      mem%growth%reserve, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%growth, 'veg_growth_fruit',        mem%growth%fruit, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%growth, 'veg_growth_seed_bed',     mem%growth%seed_bed, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%litterfall, 'veg_litterfall_leaf',         mem%litterfall%leaf, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)       
      CALL jsb_add_var(mem%litterfall, 'veg_litterfall_fine_root',    mem%litterfall%fine_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%litterfall, 'veg_litterfall_coarse_root',  mem%litterfall%coarse_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%litterfall, 'veg_litterfall_sap_wood',     mem%litterfall%sap_wood, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%litterfall, 'veg_litterfall_heart_wood',   mem%litterfall%heart_wood, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%litterfall, 'veg_litterfall_labile',       mem%litterfall%labile, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%litterfall, 'veg_litterfall_reserve',      mem%litterfall%reserve, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%litterfall, 'veg_litterfall_fruit',        mem%litterfall%fruit, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%litterfall, 'veg_litterfall_seed_bed',     mem%litterfall%seed_bed, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%falloc, 'veg_falloc_leaf',         mem%falloc%leaf, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)       
      CALL jsb_add_var(mem%falloc, 'veg_falloc_fine_root',    mem%falloc%fine_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%falloc, 'veg_falloc_coarse_root',  mem%falloc%coarse_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%falloc, 'veg_falloc_sap_wood',     mem%falloc%sap_wood, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%falloc, 'veg_falloc_heart_wood',   mem%falloc%heart_wood, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%falloc, 'veg_falloc_labile',       mem%falloc%labile, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%falloc, 'veg_falloc_reserve',      mem%falloc%reserve, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%falloc, 'veg_falloc_fruit',        mem%falloc%fruit, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%falloc, 'veg_falloc_seed_bed',     mem%falloc%seed_bed, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%total_flux, 'veg_total_flux_leaf',         mem%total_flux%leaf, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)       
      CALL jsb_add_var(mem%total_flux, 'veg_total_flux_fine_root',    mem%total_flux%fine_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%total_flux, 'veg_total_flux_coarse_root',  mem%total_flux%coarse_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%total_flux, 'veg_total_flux_sap_wood',     mem%total_flux%sap_wood, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%total_flux, 'veg_total_flux_heart_wood',   mem%total_flux%heart_wood, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%total_flux, 'veg_total_flux_labile',       mem%total_flux%labile, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%total_flux, 'veg_total_flux_reserve',      mem%total_flux%reserve, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%total_flux, 'veg_total_flux_fruit',        mem%total_flux%fruit, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%total_flux, 'veg_total_flux_seed_bed',     mem%total_flux%seed_bed, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)

      !> add variables to the structs and the l2 structs
      !!
      !! veg\%t_jsb_struct\%VAR
      !! 
      !! veg\%t_jsb_n_struct\%t_jsb_struct\%VAR
      !! 
      CALL jsb_add_var(mem%pool%leaf, 'veg_pool_leaf_carbon', mem%pool%leaf%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%leaf, 'veg_pool_leaf_nitrogen', mem%pool%leaf%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%leaf, 'veg_pool_leaf_phosphorus', mem%pool%leaf%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%leaf, 'veg_pool_leaf_carbon13', mem%pool%leaf%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%leaf, 'veg_pool_leaf_carbon14', mem%pool%leaf%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%leaf, 'veg_pool_leaf_nitrogen15', mem%pool%leaf%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%fine_root, 'veg_pool_fine_root_carbon', mem%pool%fine_root%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fine_root, 'veg_pool_fine_root_nitrogen', mem%pool%fine_root%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fine_root, 'veg_pool_fine_root_phosphorus', mem%pool%fine_root%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fine_root, 'veg_pool_fine_root_carbon13', mem%pool%fine_root%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fine_root, 'veg_pool_fine_root_carbon14', mem%pool%fine_root%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fine_root, 'veg_pool_fine_root_nitrogen15', mem%pool%fine_root%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_carbon', mem%pool%coarse_root%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_nitrogen', mem%pool%coarse_root%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_phosphorus', mem%pool%coarse_root%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_carbon13', mem%pool%coarse_root%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_carbon14', mem%pool%coarse_root%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_nitrogen15', mem%pool%coarse_root%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_carbon', mem%pool%sap_wood%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_nitrogen', mem%pool%sap_wood%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_phosphorus', mem%pool%sap_wood%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_carbon13', mem%pool%sap_wood%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_carbon14', mem%pool%sap_wood%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_nitrogen15', mem%pool%sap_wood%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_carbon', mem%pool%heart_wood%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_nitrogen', mem%pool%heart_wood%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_phosphorus', mem%pool%heart_wood%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_carbon13', mem%pool%heart_wood%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_carbon14', mem%pool%heart_wood%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_nitrogen15', mem%pool%heart_wood%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%labile, 'veg_pool_labile_carbon', mem%pool%labile%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%labile, 'veg_pool_labile_nitrogen', mem%pool%labile%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%labile, 'veg_pool_labile_phosphorus', mem%pool%labile%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%labile, 'veg_pool_labile_carbon13', mem%pool%labile%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%labile, 'veg_pool_labile_carbon14', mem%pool%labile%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%labile, 'veg_pool_labile_nitrogen15', mem%pool%labile%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%reserve, 'veg_pool_reserve_carbon', mem%pool%reserve%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%reserve, 'veg_pool_reserve_nitrogen', mem%pool%reserve%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%reserve, 'veg_pool_reserve_phosphorus', mem%pool%reserve%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%reserve, 'veg_pool_reserve_carbon13', mem%pool%reserve%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%reserve, 'veg_pool_reserve_carbon14', mem%pool%reserve%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%reserve, 'veg_pool_reserve_nitrogen15', mem%pool%reserve%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%fruit, 'veg_pool_fruit_carbon', mem%pool%fruit%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fruit, 'veg_pool_fruit_nitrogen', mem%pool%fruit%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fruit, 'veg_pool_fruit_phosphorus', mem%pool%fruit%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fruit, 'veg_pool_fruit_carbon13', mem%pool%fruit%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fruit, 'veg_pool_fruit_carbon14', mem%pool%fruit%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fruit, 'veg_pool_fruit_nitrogen15', mem%pool%fruit%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_carbon', mem%pool%seed_bed%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_nitrogen', mem%pool%seed_bed%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_phosphorus', mem%pool%seed_bed%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_carbon13', mem%pool%seed_bed%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_carbon14', mem%pool%seed_bed%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_nitrogen15', mem%pool%seed_bed%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%growth%leaf, 'veg_growth_leaf_carbon', mem%growth%leaf%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%leaf, 'veg_growth_leaf_nitrogen', mem%growth%leaf%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%leaf, 'veg_growth_leaf_phosphorus', mem%growth%leaf%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%leaf, 'veg_growth_leaf_carbon13', mem%growth%leaf%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%leaf, 'veg_growth_leaf_carbon14', mem%growth%leaf%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%leaf, 'veg_growth_leaf_nitrogen15', mem%growth%leaf%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%growth%fine_root, 'veg_growth_fine_root_carbon', mem%growth%fine_root%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%fine_root, 'veg_growth_fine_root_nitrogen', mem%growth%fine_root%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%fine_root, 'veg_growth_fine_root_phosphorus', mem%growth%fine_root%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%fine_root, 'veg_growth_fine_root_carbon13', mem%growth%fine_root%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%fine_root, 'veg_growth_fine_root_carbon14', mem%growth%fine_root%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%fine_root, 'veg_growth_fine_root_nitrogen15', mem%growth%fine_root%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_carbon', mem%growth%coarse_root%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_nitrogen', mem%growth%coarse_root%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_phosphorus', mem%growth%coarse_root%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_carbon13', mem%growth%coarse_root%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_carbon14', mem%growth%coarse_root%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_nitrogen15', mem%growth%coarse_root%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_carbon', mem%growth%sap_wood%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_nitrogen', mem%growth%sap_wood%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_phosphorus', mem%growth%sap_wood%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_carbon13', mem%growth%sap_wood%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_carbon14', mem%growth%sap_wood%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_nitrogen15', mem%growth%sap_wood%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_carbon', mem%growth%heart_wood%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_nitrogen', mem%growth%heart_wood%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_phosphorus', mem%growth%heart_wood%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_carbon13', mem%growth%heart_wood%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_carbon14', mem%growth%heart_wood%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_nitrogen15', mem%growth%heart_wood%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%growth%labile, 'veg_growth_labile_carbon', mem%growth%labile%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%labile, 'veg_growth_labile_nitrogen', mem%growth%labile%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%labile, 'veg_growth_labile_phosphorus', mem%growth%labile%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%labile, 'veg_growth_labile_carbon13', mem%growth%labile%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%labile, 'veg_growth_labile_carbon14', mem%growth%labile%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%labile, 'veg_growth_labile_nitrogen15', mem%growth%labile%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%growth%reserve, 'veg_growth_reserve_carbon', mem%growth%reserve%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%reserve, 'veg_growth_reserve_nitrogen', mem%growth%reserve%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%reserve, 'veg_growth_reserve_phosphorus', mem%growth%reserve%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%reserve, 'veg_growth_reserve_carbon13', mem%growth%reserve%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%reserve, 'veg_growth_reserve_carbon14', mem%growth%reserve%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%reserve, 'veg_growth_reserve_nitrogen15', mem%growth%reserve%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%growth%fruit, 'veg_growth_fruit_carbon', mem%growth%fruit%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%fruit, 'veg_growth_fruit_nitrogen', mem%growth%fruit%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%fruit, 'veg_growth_fruit_phosphorus', mem%growth%fruit%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%fruit, 'veg_growth_fruit_carbon13', mem%growth%fruit%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%fruit, 'veg_growth_fruit_carbon14', mem%growth%fruit%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%fruit, 'veg_growth_fruit_nitrogen15', mem%growth%fruit%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_carbon', mem%growth%seed_bed%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_nitrogen', mem%growth%seed_bed%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_phosphorus', mem%growth%seed_bed%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_carbon13', mem%growth%seed_bed%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_carbon14', mem%growth%seed_bed%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_nitrogen15', mem%growth%seed_bed%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_carbon', mem%litterfall%leaf%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_nitrogen', mem%litterfall%leaf%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_phosphorus', mem%litterfall%leaf%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_carbon13', mem%litterfall%leaf%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_carbon14', mem%litterfall%leaf%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_nitrogen15', mem%litterfall%leaf%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_carbon', mem%litterfall%fine_root%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_nitrogen', mem%litterfall%fine_root%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_phosphorus', mem%litterfall%fine_root%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_carbon13', mem%litterfall%fine_root%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_carbon14', mem%litterfall%fine_root%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_nitrogen15', mem%litterfall%fine_root%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_carbon', mem%litterfall%coarse_root%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_nitrogen', mem%litterfall%coarse_root%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_phosphorus', mem%litterfall%coarse_root%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_carbon13', mem%litterfall%coarse_root%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_carbon14', mem%litterfall%coarse_root%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_nitrogen15', mem%litterfall%coarse_root%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_carbon', mem%litterfall%sap_wood%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_nitrogen', mem%litterfall%sap_wood%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_phosphorus', mem%litterfall%sap_wood%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_carbon13', mem%litterfall%sap_wood%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_carbon14', mem%litterfall%sap_wood%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_nitrogen15', mem%litterfall%sap_wood%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_carbon', mem%litterfall%heart_wood%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_nitrogen', mem%litterfall%heart_wood%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_phosphorus', mem%litterfall%heart_wood%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_carbon13', mem%litterfall%heart_wood%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_carbon14', mem%litterfall%heart_wood%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_nitrogen15', mem%litterfall%heart_wood%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%litterfall%labile, 'veg_litterfall_labile_carbon', mem%litterfall%labile%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%labile, 'veg_litterfall_labile_nitrogen', mem%litterfall%labile%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%labile, 'veg_litterfall_labile_phosphorus', mem%litterfall%labile%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%labile, 'veg_litterfall_labile_carbon13', mem%litterfall%labile%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%labile, 'veg_litterfall_labile_carbon14', mem%litterfall%labile%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%labile, 'veg_litterfall_labile_nitrogen15', mem%litterfall%labile%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_carbon', mem%litterfall%reserve%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_nitrogen', mem%litterfall%reserve%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_phosphorus', mem%litterfall%reserve%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_carbon13', mem%litterfall%reserve%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_carbon14', mem%litterfall%reserve%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_nitrogen15', mem%litterfall%reserve%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_carbon', mem%litterfall%fruit%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_nitrogen', mem%litterfall%fruit%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_phosphorus', mem%litterfall%fruit%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_carbon13', mem%litterfall%fruit%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_carbon14', mem%litterfall%fruit%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_nitrogen15', mem%litterfall%fruit%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_carbon', mem%litterfall%seed_bed%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_nitrogen', mem%litterfall%seed_bed%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_phosphorus', mem%litterfall%seed_bed%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_carbon13', mem%litterfall%seed_bed%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_carbon14', mem%litterfall%seed_bed%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_nitrogen15', mem%litterfall%seed_bed%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%falloc%leaf, 'veg_falloc_leaf_carbon', mem%falloc%leaf%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%leaf, 'veg_falloc_leaf_nitrogen', mem%falloc%leaf%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%leaf, 'veg_falloc_leaf_phosphorus', mem%falloc%leaf%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%leaf, 'veg_falloc_leaf_carbon13', mem%falloc%leaf%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%leaf, 'veg_falloc_leaf_carbon14', mem%falloc%leaf%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%leaf, 'veg_falloc_leaf_nitrogen15', mem%falloc%leaf%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%falloc%fine_root, 'veg_falloc_fine_root_carbon', mem%falloc%fine_root%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%fine_root, 'veg_falloc_fine_root_nitrogen', mem%falloc%fine_root%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%fine_root, 'veg_falloc_fine_root_phosphorus', mem%falloc%fine_root%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%fine_root, 'veg_falloc_fine_root_carbon13', mem%falloc%fine_root%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%fine_root, 'veg_falloc_fine_root_carbon14', mem%falloc%fine_root%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%fine_root, 'veg_falloc_fine_root_nitrogen15', mem%falloc%fine_root%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%falloc%coarse_root, 'veg_falloc_coarse_root_carbon', mem%falloc%coarse_root%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%coarse_root, 'veg_falloc_coarse_root_nitrogen', mem%falloc%coarse_root%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%coarse_root, 'veg_falloc_coarse_root_phosphorus', mem%falloc%coarse_root%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%coarse_root, 'veg_falloc_coarse_root_carbon13', mem%falloc%coarse_root%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%coarse_root, 'veg_falloc_coarse_root_carbon14', mem%falloc%coarse_root%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%coarse_root, 'veg_falloc_coarse_root_nitrogen15', mem%falloc%coarse_root%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%falloc%sap_wood, 'veg_falloc_sap_wood_carbon', mem%falloc%sap_wood%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%sap_wood, 'veg_falloc_sap_wood_nitrogen', mem%falloc%sap_wood%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%sap_wood, 'veg_falloc_sap_wood_phosphorus', mem%falloc%sap_wood%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%sap_wood, 'veg_falloc_sap_wood_carbon13', mem%falloc%sap_wood%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%sap_wood, 'veg_falloc_sap_wood_carbon14', mem%falloc%sap_wood%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%sap_wood, 'veg_falloc_sap_wood_nitrogen15', mem%falloc%sap_wood%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%falloc%heart_wood, 'veg_falloc_heart_wood_carbon', mem%falloc%heart_wood%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%heart_wood, 'veg_falloc_heart_wood_nitrogen', mem%falloc%heart_wood%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%heart_wood, 'veg_falloc_heart_wood_phosphorus', mem%falloc%heart_wood%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%heart_wood, 'veg_falloc_heart_wood_carbon13', mem%falloc%heart_wood%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%heart_wood, 'veg_falloc_heart_wood_carbon14', mem%falloc%heart_wood%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%heart_wood, 'veg_falloc_heart_wood_nitrogen15', mem%falloc%heart_wood%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%falloc%labile, 'veg_falloc_labile_carbon', mem%falloc%labile%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%labile, 'veg_falloc_labile_nitrogen', mem%falloc%labile%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%labile, 'veg_falloc_labile_phosphorus', mem%falloc%labile%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%labile, 'veg_falloc_labile_carbon13', mem%falloc%labile%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%labile, 'veg_falloc_labile_carbon14', mem%falloc%labile%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%labile, 'veg_falloc_labile_nitrogen15', mem%falloc%labile%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%falloc%reserve, 'veg_falloc_reserve_carbon', mem%falloc%reserve%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%reserve, 'veg_falloc_reserve_nitrogen', mem%falloc%reserve%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%reserve, 'veg_falloc_reserve_phosphorus', mem%falloc%reserve%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%reserve, 'veg_falloc_reserve_carbon13', mem%falloc%reserve%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%reserve, 'veg_falloc_reserve_carbon14', mem%falloc%reserve%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%reserve, 'veg_falloc_reserve_nitrogen15', mem%falloc%reserve%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%falloc%fruit, 'veg_falloc_fruit_carbon', mem%falloc%fruit%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%fruit, 'veg_falloc_fruit_nitrogen', mem%falloc%fruit%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%fruit, 'veg_falloc_fruit_phosphorus', mem%falloc%fruit%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%fruit, 'veg_falloc_fruit_carbon13', mem%falloc%fruit%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%fruit, 'veg_falloc_fruit_carbon14', mem%falloc%fruit%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%fruit, 'veg_falloc_fruit_nitrogen15', mem%falloc%fruit%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%falloc%seed_bed, 'veg_falloc_seed_bed_carbon', mem%falloc%seed_bed%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%seed_bed, 'veg_falloc_seed_bed_nitrogen', mem%falloc%seed_bed%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%seed_bed, 'veg_falloc_seed_bed_phosphorus', mem%falloc%seed_bed%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%seed_bed, 'veg_falloc_seed_bed_carbon13', mem%falloc%seed_bed%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%seed_bed, 'veg_falloc_seed_bed_carbon14', mem%falloc%seed_bed%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%falloc%seed_bed, 'veg_falloc_seed_bed_nitrogen15', mem%falloc%seed_bed%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue

      CALL jsb_add_var(mem%total_flux%leaf, 'veg_total_flux_leaf_carbon', mem%total_flux%leaf%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%leaf, 'veg_total_flux_leaf_nitrogen', mem%total_flux%leaf%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%leaf, 'veg_total_flux_leaf_phosphorus', mem%total_flux%leaf%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%leaf, 'veg_total_flux_leaf_carbon13', mem%total_flux%leaf%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%leaf, 'veg_total_flux_leaf_carbon14', mem%total_flux%leaf%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%leaf, 'veg_total_flux_leaf_nitrogen15', mem%total_flux%leaf%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%total_flux%fine_root, 'veg_total_flux_fine_root_carbon', mem%total_flux%fine_root%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%fine_root, 'veg_total_flux_fine_root_nitrogen', mem%total_flux%fine_root%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%fine_root, 'veg_total_flux_fine_root_phosphorus', mem%total_flux%fine_root%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%fine_root, 'veg_total_flux_fine_root_carbon13', mem%total_flux%fine_root%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%fine_root, 'veg_total_flux_fine_root_carbon14', mem%total_flux%fine_root%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%fine_root, 'veg_total_flux_fine_root_nitrogen15', mem%total_flux%fine_root%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%total_flux%coarse_root, 'veg_total_flux_coarse_root_carbon', mem%total_flux%coarse_root%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%coarse_root, 'veg_total_flux_coarse_root_nitrogen', mem%total_flux%coarse_root%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%coarse_root, 'veg_total_flux_coarse_root_phosphorus', mem%total_flux%coarse_root%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%coarse_root, 'veg_total_flux_coarse_root_carbon13', mem%total_flux%coarse_root%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%coarse_root, 'veg_total_flux_coarse_root_carbon14', mem%total_flux%coarse_root%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%coarse_root, 'veg_total_flux_coarse_root_nitrogen15', mem%total_flux%coarse_root%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%total_flux%sap_wood, 'veg_total_flux_sap_wood_carbon', mem%total_flux%sap_wood%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%sap_wood, 'veg_total_flux_sap_wood_nitrogen', mem%total_flux%sap_wood%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%sap_wood, 'veg_total_flux_sap_wood_phosphorus', mem%total_flux%sap_wood%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%sap_wood, 'veg_total_flux_sap_wood_carbon13', mem%total_flux%sap_wood%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%sap_wood, 'veg_total_flux_sap_wood_carbon14', mem%total_flux%sap_wood%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%sap_wood, 'veg_total_flux_sap_wood_nitrogen15', mem%total_flux%sap_wood%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%total_flux%heart_wood, 'veg_total_flux_heart_wood_carbon', mem%total_flux%heart_wood%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%heart_wood, 'veg_total_flux_heart_wood_nitrogen', mem%total_flux%heart_wood%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%heart_wood, 'veg_total_flux_heart_wood_phosphorus', mem%total_flux%heart_wood%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%heart_wood, 'veg_total_flux_heart_wood_carbon13', mem%total_flux%heart_wood%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%heart_wood, 'veg_total_flux_heart_wood_carbon14', mem%total_flux%heart_wood%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%heart_wood, 'veg_total_flux_heart_wood_nitrogen15', mem%total_flux%heart_wood%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%total_flux%labile, 'veg_total_flux_labile_carbon', mem%total_flux%labile%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%labile, 'veg_total_flux_labile_nitrogen', mem%total_flux%labile%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%labile, 'veg_total_flux_labile_phosphorus', mem%total_flux%labile%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%labile, 'veg_total_flux_labile_carbon13', mem%total_flux%labile%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%labile, 'veg_total_flux_labile_carbon14', mem%total_flux%labile%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%labile, 'veg_total_flux_labile_nitrogen15', mem%total_flux%labile%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%total_flux%reserve, 'veg_total_flux_reserve_carbon', mem%total_flux%reserve%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%reserve, 'veg_total_flux_reserve_nitrogen', mem%total_flux%reserve%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%reserve, 'veg_total_flux_reserve_phosphorus', mem%total_flux%reserve%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%reserve, 'veg_total_flux_reserve_carbon13', mem%total_flux%reserve%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%reserve, 'veg_total_flux_reserve_carbon14', mem%total_flux%reserve%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%reserve, 'veg_total_flux_reserve_nitrogen15', mem%total_flux%reserve%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%total_flux%fruit, 'veg_total_flux_fruit_carbon', mem%total_flux%fruit%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%fruit, 'veg_total_flux_fruit_nitrogen', mem%total_flux%fruit%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%fruit, 'veg_total_flux_fruit_phosphorus', mem%total_flux%fruit%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%fruit, 'veg_total_flux_fruit_carbon13', mem%total_flux%fruit%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%fruit, 'veg_total_flux_fruit_carbon14', mem%total_flux%fruit%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%fruit, 'veg_total_flux_fruit_nitrogen15', mem%total_flux%fruit%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%total_flux%seed_bed, 'veg_total_flux_seed_bed_carbon', mem%total_flux%seed_bed%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%seed_bed, 'veg_total_flux_seed_bed_nitrogen', mem%total_flux%seed_bed%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%seed_bed, 'veg_total_flux_seed_bed_phosphorus', mem%total_flux%seed_bed%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%seed_bed, 'veg_total_flux_seed_bed_carbon13', mem%total_flux%seed_bed%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%seed_bed, 'veg_total_flux_seed_bed_carbon14', mem%total_flux%seed_bed%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_flux%seed_bed, 'veg_total_flux_seed_bed_nitrogen15', mem%total_flux%seed_bed%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue

      CALL jsb_add_var(mem%exudation, 'veg_exudation_carbon', mem%exudation%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%exudation, 'veg_exudation_nitrogen', mem%exudation%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%exudation, 'veg_exudation_phosphorus', mem%exudation%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%exudation, 'veg_exudation_carbon13', mem%exudation%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%exudation, 'veg_exudation_carbon14', mem%exudation%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%exudation, 'veg_exudation_nitrogen15', mem%exudation%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
 
      CALL jsb_add_var(mem%establishment, 'veg_establishment_carbon', mem%establishment%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%establishment, 'veg_establishment_nitrogen', mem%establishment%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%establishment, 'veg_establishment_phosphorus', mem%establishment%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%establishment, 'veg_establishment_carbon13', mem%establishment%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%establishment, 'veg_establishment_carbon14', mem%establishment%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%establishment, 'veg_establishment_nitrogen15', mem%establishment%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%reserve_use, 'veg_reserve_use_carbon', mem%reserve_use%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%reserve_use, 'veg_reserve_use_nitrogen', mem%reserve_use%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%reserve_use, 'veg_reserve_use_phosphorus', mem%reserve_use%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%reserve_use, 'veg_reserve_use_carbon13', mem%reserve_use%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%reserve_use, 'veg_reserve_use_carbon14', mem%reserve_use%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%reserve_use, 'veg_reserve_use_nitrogen15', mem%reserve_use%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue

      CALL jsb_add_var(mem%total_biomass, 'veg_total_biomass_carbon', mem%total_biomass%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_biomass, 'veg_total_biomass_nitrogen', mem%total_biomass%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_biomass, 'veg_total_biomass_phosphorus', mem%total_biomass%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_biomass, 'veg_total_biomass_carbon13', mem%total_biomass%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_biomass, 'veg_total_biomass_carbon14', mem%total_biomass%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%total_biomass, 'veg_total_biomass_nitrogen15', mem%total_biomass%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue

    END IF

  END SUBROUTINE Init_veg_memory

END MODULE mo_veg_memory_class
