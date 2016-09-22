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

  USE mo_kind, ONLY: wp

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d         , &
                                       t_jsb_var_logical2d                        , &                 ! *slmjsb4
                                       t_jsb_struct, t_jsb_n_struct               , &                 ! *slmjsb4
                                       t_jsb_var_real1d, t_jsb_var_logical1d                          ! *slmtemp

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later    *slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_veg_memory, max_no_of_vars
  PUBLIC :: max_no_of_structs, max_no_of_n_structs  ! *slmjsb4

  INTEGER, PARAMETER :: max_no_of_vars      = 160  !> TYPE t_jsb_var
  INTEGER, PARAMETER :: max_no_of_structs   = 20   !> TYPE t_jsb_struct     *slmjsb4
  INTEGER, PARAMETER :: max_no_of_n_structs = 20   !> TYPE t_jsb_n_struct   *slmjsb4



  !-----------------------------------------------------------------------------------------------------
  !> Type definition for elements in veg pool
  !! 
  !! elements in veg would be 2D (while they would be 3D in soil)
  !! 
  !! *slmtemp 1D -> 2D
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_struct) :: t_slm_struct_veg_element
    TYPE(t_jsb_var_real1d)    :: carbon       !< total carbon concentration of pool in mol / m2
    TYPE(t_jsb_var_real1d)    :: nitrogen     !< totol nitrogen concentration of pool in mol / m2
    TYPE(t_jsb_var_real1d)    :: phosphorus   !< phosphorus concentration of pool in mol / m2
    TYPE(t_jsb_var_real1d)    :: carbon13     !< 13carbon concentration of pool in mol / m2
    TYPE(t_jsb_var_real1d)    :: carbon14     !< 14carbon concentration of pool in pmol / m2
    TYPE(t_jsb_var_real1d)    :: nitrogen15   !< 15nitrogen concentration of pool in mol / m2
  CONTAINS
    PRIVATE
    PROCEDURE               :: Add_struct                   ! Add two variables of type t_slm_struct_veg_element
    PROCEDURE               :: Add_struct_and_scalar        ! Add variable of type t_slm_struct_veg_element with REAL scalar
    PROCEDURE, PASS(this2)  :: Add_scalar_and_struct        ! Add variable of type t_slm_struct_veg_element with REAL scalar
!    PROCEDURE               :: Add_struct_and_field2d       ! Add variable of type t_slm_struct_veg_element with REAL 2-d field2d
!    PROCEDURE, PASS(this2)  :: Add_field2d_and_struct       ! Add variable of REAL 2-d field2d with type t_slm_struct_veg_element
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_struct, Add_struct_and_scalar, Add_scalar_and_struct!, &
!                                              Add_struct_and_field2d, Add_field2d_and_struct
    PROCEDURE               :: Subtract_struct              ! Subtract two variables of type t_slm_struct_veg_element
    PROCEDURE               :: Subtract_struct_and_scalar   ! Subtract variable of type t_slm_struct_veg_element with REAL scalar
    PROCEDURE, PASS(this2)  :: Subtract_scalar_and_struct   ! Subtract variable of type t_slm_struct_veg_element with REAL scalar
!    PROCEDURE               :: Subtract_struct_and_field2d  ! Subtract variable of type t_slm_struct_veg_element with REAL 2-d field2d
!    PROCEDURE, PASS(this2)  :: Subtract_field2d_and_struct  ! Subtract variable of REAL 2-d field2d with type t_slm_struct_veg_element
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_struct, Subtract_struct_and_scalar, Subtract_scalar_and_struct!, &
!                                              Subtract_struct_and_field2d, Subtract_field2d_and_struct
    PROCEDURE               :: Multiply_struct              ! Multiply two variables of type t_slm_struct_veg_element
    PROCEDURE               :: Multiply_struct_and_scalar   ! Multiply variable of type t_slm_struct_veg_element with REAL scalar
    PROCEDURE, PASS(this2)  :: Multiply_scalar_and_struct   ! Multiply variable of type t_slm_struct_veg_element with REAL scalar
!    PROCEDURE               :: Multiply_struct_and_field2d  ! Multiply variable of type t_slm_struct_veg_element with REAL 2-d field2d
!    PROCEDURE, PASS(this2)  :: Multiply_field2d_and_struct  ! Multiply variable of REAL 2-d field2d with type t_slm_struct_veg_element
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_struct, Multiply_struct_and_scalar, Multiply_scalar_and_struct!, &
!                                              Multiply_struct_and_field2d, Multiply_field2d_and_struct
    PROCEDURE               :: Divide_struct                ! Divide two variables of type t_slm_struct_veg_element
    PROCEDURE               :: Divide_struct_and_scalar     ! Divide variable of type t_slm_struct_veg_element by REAL scalar
    PROCEDURE, PASS(this2)  :: Divide_scalar_and_struct     ! Divide variable of REAL scalar by type t_slm_struct_veg_element
!    PROCEDURE               :: Divide_struct_and_field2d    ! Divide variable of type t_slm_struct_veg_element by REAL 2-d field2d
!    PROCEDURE, PASS(this2)  :: Divide_field2d_and_struct    ! Divide variable of REAL 2-d field2d by type t_slm_struct_veg_element
    GENERIC, PUBLIC         :: OPERATOR(/) => Divide_struct, Divide_struct_and_scalar, Divide_scalar_and_struct!, &
!                                              Divide_struct_and_field2d, Divide_field2d_and_struct
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
    PRIVATE
    PROCEDURE               :: Add_n_struct                   
    PROCEDURE               :: Add_n_struct_and_scalar        
    PROCEDURE, PASS(this2)  :: Add_scalar_and_n_struct        
!    PROCEDURE               :: Add_n_struct_and_field2d      
!    PROCEDURE, PASS(this2)  :: Add_field2d_and_n_struct      
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_n_struct, Add_n_struct_and_scalar, Add_scalar_and_n_struct!, &
!                                              Add_n_struct_and_field2d, Add_field2d_and_n_struct
    PROCEDURE               :: Subtract_n_struct              
    PROCEDURE               :: Subtract_n_struct_and_scalar   
    PROCEDURE, PASS(this2)  :: Subtract_scalar_and_n_struct   
!    PROCEDURE               :: Subtract_n_struct_and_field2d  
!    PROCEDURE, PASS(this2)  :: Subtract_field2d_and_n_struct  
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_n_struct, Subtract_n_struct_and_scalar, Subtract_scalar_and_n_struct!, &
!                                              Subtract_n_struct_and_field2d, Subtract_field2d_and_n_struct
    PROCEDURE               :: Multiply_n_struct              
    PROCEDURE               :: Multiply_n_struct_and_scalar   
    PROCEDURE, PASS(this2)  :: Multiply_scalar_and_n_struct   
!    PROCEDURE               :: Multiply_n_struct_and_field2d  
!    PROCEDURE, PASS(this2)  :: Multiply_field2d_and_n_struct  
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_n_struct, Multiply_n_struct_and_scalar, Multiply_scalar_and_n_struct!, &
!                                              Multiply_n_struct_and_field2d, Multiply_field2d_and_n_struct
    PROCEDURE               :: Divide_n_struct                
    PROCEDURE               :: Divide_n_struct_and_scalar     
    PROCEDURE, PASS(this2)  :: Divide_scalar_and_n_struct     
!    PROCEDURE               :: Divide_n_struct_and_field2d    
!    PROCEDURE, PASS(this2)  :: Divide_field2d_and_n_struct    
    GENERIC, PUBLIC         :: OPERATOR(/) => Divide_n_struct, Divide_n_struct_and_scalar, Divide_scalar_and_n_struct!, &
!                                              Divide_n_struct_and_field2d, Divide_field2d_and_n_struct
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

! *slmtemp  --> 2D
    !> variables
    !! 
    ! logical
    TYPE(t_jsb_var_logical1d)      :: growing_season    !< logical to identify whether the plant is in the growing season or not
    ! basic plant
    TYPE(t_jsb_var_real1d)          :: &   
                                    height            , &            !< plant height in m
                                    diameter          , &            !< plant diameter in m
                                    dens_ind          , &            !< density of individuals (#/m2) 
                                    lai               , &            !< actual leaf area index
                                    target_cn_leaf    , &            !< target foliar C:N ratio
                                    target_np_leaf    , &            !< target foliar N:P ratio 
                                    mean_leaf_age     , &            !< average foliage age in days 
                                    gc                , &            !< some kind of stomatal conductance ?
                                    f_n_demand        , &            !< relative plant demand for N; see update_vegsoil_bgc
                                    f_p_demand        , &            !< relative plant demand for P; see update_vegsoil_bgc
                                    beta_air          , &            !< scaling factor to account for air humidity effects on conductance
                                    beta_soil_ps      , &            !< scaling factor to account for soil moisture constraints on conductance
                                    dphi              , &            !< plant-soil water potential gradient (MPa)
                                    beta_soil_gs      , &            !< scaling factor to account for soil moisture constraints on conductance
                                    growth_req_n      , &            !< moles N required for a unit of C growth under current allocation fractions 
                                    growth_req_p      , &            !< moles P required for a unit of N growth under current allocation fractions 
                                    gdd               , &            !< growing degree days above 5°C
                                    nd_dormance       , &            !< number of days since last growing season 
                                    t_jmax_opt        , &            !< current optimum temperature for eletron transport of photosynthesis
                                    an_boc            , &            !< net photosynthesis at bottom of canopy (lowest canopy layer existing)
                                    leaf_cn_direction , &            !< indicates direction of change in leaf CN ratio to maximise NPP
                                    apar                             !< ..
    ! canopy layers -> levels, lbounds, ubounds
    TYPE(t_jsb_var_real2d)          :: &   ! DIMENSION(ncanop)
                                    nl_tot_cl         , &       !< N in canopy layer (mmol / m2 (LAI))
                                    fnl_chl_cl        , &       !< fraction of N in chlorophyll (pigments only!)
                                    fnl_et_cl         , &       !< fraction of N in electron transport 
                                    fnl_rub_cl        , &       !< fraction of N in Rubisco 
                                    fnl_pepc_cl       , &       !< fraction of N in PEP carboylase (C4 plants only)
                                    fnl_oth_cl        , &       !< fraction of N not associated with photosynthesis 
                                    par_sunlit_cl     , &       !< PAR on sunlit leaves (µmol m-2 s-1)
                                    par_shaded_cl     , &       !< PAR on shaded leaves (µmol m-2 s-1)
                                    fleaf_sunlit_cl   , &       !< fraction of layer which is sunlit
                                    lai_cl            , &       !< leaf area in the canopy layer 
                                    cumm_lai_cl                 !< cummulative leaf area above the mid-point of the canopy layer 

    ! fluxes
    TYPE(t_jsb_var_real1d)          :: &   
                                    gpp                       , & !< gross photosynthesis (µmol CO2 m-2 s-1)
                                    an                        , & !< net photosynthesis (µmol CO2 m-2 s-1)
                                    gpp_c13                   , & !< gross photosynthesis (µmol 13CO2 m-2 s-1)
                                    gpp_c14                   , & !< gross photosynthesis (µmol 14CO2 m-2 s-1)
                                    maint_respiration_leaf    , & !< foliar maintenance respiration (µmol CO2 m-2 s-1)
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
                                    uptake_n                  , & !< N uptake from soil ( is that correct? )
                                    uptake_n15                , & !< 15N uptake from soil ( is that correct? ) 
                                    uptake_p                  , & !< P uptake from soil ( is that correct? )
                                    n_fixation                , & !< symbiontic N fixation (µmol N m-2 s-1)
                                    n_fixation_n15            , & !< symbiontic N fixation (µmol 15N m-2 s-1)
                                    recycling_leaf_n          , & !< net flux of N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
                                    recycling_leaf_n15        , & !< net flux of 15N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1) 
                                    recycling_leaf_p          , & !< net flux of P from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
                                    recycling_fine_root_n     , & !< net flux of N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
                                    recycling_fine_root_n15   , & !< net flux of 15N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1) 
                                    recycling_fine_root_p     , & !< net flux of P from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
                                    recycling_sap_wood_n      , & !< net flux of N from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1)
                                    recycling_sap_wood_n15    , & !< net flux of 15N from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1) 
                                    recycling_sap_wood_p      , & !< net flux of P from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1)
                                    transpiration                 !< transpiration in kg m-2 s-1
    ! everything else...
    TYPE(t_jsb_var_real1d)          :: &   
                                    fmaint_rate               , & !< current N-specific maintentance respiration (~µmol C/mol N / s ) 
                                    ci                        , & !< leaf internal co2 concentration
                                    unit_npp                  , & !< incremental change in NPP. Needed for optimal allocation
                                    unit_transpiration        , & !< incremental change in transpiration. Needed for optimal allocation
                                    unit_uptake_n             , & !< incremental change in N-based growth, needed for optimal allocation
                                    unit_uptake_p             , & !< incremental change in P-based growth, needed for optimal allocation
                                    delta_dens_ind            , & !< change in density of individuals (#/m2/timestep)
                                    cost_n_uptake_root        , & !< actual cost for N uptake by roots (mol C/mol N)
                                    cost_n_fixation               !< actual cost for N uptake through fixation (mol C/mol N)

    ! daytime averages
    TYPE(t_jsb_var_real2d)          :: &   ! DIMENSION(ncanop)
                                    par_sunlit_daytime_cl     , &    !<  average of the previous day
                                    par_shaded_daytime_cl     , &    !<  average of the previous day
                                    fleaf_sunlit_daytime_cl          !<  average of the previous day
    TYPE(t_jsb_var_real1d)          :: &   
                                    t_air_daytime             , &    !<  average of the previous day
                                    press_srf_daytime         , &    !<  average of the previous day
                                    co2_mixing_ratio_daytime  , &    !<  average of the previous day
                                    ga_daytime                , &    !<  average of the previous day
                                    beta_air_daytime          , &    !<  average of the previous day
                                    beta_soil_ps_daytime      , &    !<  average of the previous day
                                    beta_soil_gs_daytime      , &    !<  average of the previous day
                                    t_jmax_opt_daytime               !<  average of the previous day
    TYPE(t_jsb_var_real2d)          :: &   ! DIMENSION(ncanop)
                                    par_sunlit_daytime_dacc_cl   , & !<  value accumulated over current day
                                    par_shaded_daytime_dacc_cl   , & !<  value accumulated over current day
                                    fleaf_sunlit_daytime_dacc_cl     !<  value accumulated over current day
    TYPE(t_jsb_var_real1d)          :: &   
                                    t_air_daytime_dacc            , & !<  average of the previous day
                                    press_srf_daytime_dacc        , & !<  average of the previous day
                                    co2_mixing_ratio_daytime_dacc , & !<  average of the previous day
                                    ga_daytime_dacc               , & !<  average of the previous day
                                    beta_air_daytime_dacc         , & !<  average of the previous day
                                    beta_soil_ps_daytime_dacc     , & !<  average of the previous day
                                    beta_soil_gs_daytime_dacc     , & !<  average of the previous day
                                    t_jmax_opt_daytime_dacc           !<  value accumulated over current day 
    TYPE(t_jsb_var_real1d)          :: &   
                                    daytime_counter  !< number of timesteps of current day with light available & daytime_avrg variables were added up

    ! moving averages
    TYPE(t_jsb_var_real1d)          :: &   
    ! 1.1 memory of conditions for the labile pool dynamics
                                    gpp_tlabile_mavg               , &  !< gpp  (µmol m-2 s-1)
                                    maint_respiration_tlabile_mavg , &  !< maintenance respiration (µmol m-2 s-1)
                                    growth_req_n_tlabile_mavg      , &  !< moles N required for a unit of C growth under current allocation fractions
                                    growth_req_p_tlabile_mavg           !< moles P required for a unit of N growth under current allocation fractions
    ! 1.2 memory of conditions for phenology calculations
    TYPE(t_jsb_var_real1d)          :: &   
                                    t_air_tphen_mavg               , &  !< air temp mavg for phenology calc
                                    beta_soil_gs_tphen_mavg             !< scaling factor to account for soil moisture constraints on conductance
    ! 1.3 memory for plant uptake (average demand over a couple of days to get rid of diurnal C cycle)
    TYPE(t_jsb_var_real1d)          :: &   
                                    npp_tuptake_mavg               , &  !< npp of the last week (µmol m-2 s-1)
                                    demand_uptake_n_tuptake_mavg   , &  !< scalar of Vmax for root N uptake 
                                    demand_uptake_p_tuptake_mavg   , &  !< scalar of Vmax for root P uptake
                                    growth_req_n_tuptake_mavg      , &  !< moles N required for a unit of C growth under current allocation fractions
                                    growth_req_p_tuptake_mavg      , &  !< moles P required for a unit of N growth under current allocation fractions
                                    cost_n_fixation_tuptake_mavg        !< moving average of N fixation cost (mol C / mol N)
    ! 1.4 memory of conditions for photosynthesis for the optimisation of foliar N fractions 
    TYPE(t_jsb_var_real2d)          :: &   ! DIMENSION(ncanop)
                                    par_sunlit_tfrac_mavg_cl       , &  !< PAR on sunlit leaves
                                    par_shaded_tfrac_mavg_cl       , &  !< PAR on shaded leaves
                                    fleaf_sunlit_tfrac_mavg_cl          !< fraction of sunlit leaves
    TYPE(t_jsb_var_real1d)          :: &   
                                    t_air_tfrac_mavg               , &  !< air temp
                                    press_srf_tfrac_mavg           , &  !< air pressure at surface
                                    co2_mixing_ratio_tfrac_mavg    , &  !< co2 mixing ratio of the air
                                    ga_tfrac_mavg                  , &  !< aerodynamic resistance
                                    beta_air_tfrac_mavg            , &  !< scaling factor to account for air humidity effects on conductance
                                    beta_soil_ps_tfrac_mavg        , &  !< scaling factor to account for soil moisture constraints on photosynthesis 
                                    beta_soil_gs_tfrac_mavg             !< scaling factor to account for soil moisture constraints on conductance
    ! 1.5 memory of conditions for photosynthesis for the optimisation of foliar N  
    TYPE(t_jsb_var_real2d)          :: &   ! DIMENSION(ncanop)
                                    par_sunlit_tcnl_mavg_cl        , &  !< PAR on sunlit leaves
                                    par_shaded_tcnl_mavg_cl        , &  !< PAR on shaded leaves
                                    fleaf_sunlit_tcnl_mavg_cl           !< fraction of sunlit leaves
    TYPE(t_jsb_var_real1d)          :: &   
                                    t_air_tcnl_mavg                , &  !< air temp
                                    press_srf_tcnl_mavg            , &  !< air pressure at surface
                                    co2_mixing_ratio_tcnl_mavg     , &  !< co2 mixing ratio of the air
                                    ga_tcnl_mavg                   , &  !< aerodynamic resistance
                                    beta_air_tcnl_mavg             , &  !< scaling factor to account for air humidity effects on conductance
                                    beta_soil_ps_tcnl_mavg         , &  !< scaling factor to account for soil moisture constraints on conductance
                                    beta_soil_gs_tcnl_mavg              !< scaling factor to account for soil moisture constraints on conductance
    ! 1.6 memory of conditions for photosynthesis for the optimisation of biomass allocation fractions
    TYPE(t_jsb_var_real1d)          :: &   
                                    npp_talloc_mavg                , &  !< NPP (µmol m-2 s-1) (averaging time-scale: tau_allocation)
                                    transpiration_talloc_mavg      , &  !< average transpiration at allocation time-scale
                                    unit_npp_talloc_mavg           , &  !< rate of change in NPP with  increment in leaf mass, mumol C mol-1
                                    unit_transpiration_talloc_mavg , &  !< rate of tranpiration change per unit LAI change?
                                    unit_uptake_n_talloc_mavg      , &  !< rate of N uptkae per unit root, mumol N mol-1 (averaging time-scale: tau_allocation)
                                    unit_uptake_p_talloc_mavg      , &  !< rate of P uptake per unit root, mumol P mol-1 (time-scale of accumulation: tau_allocation)
                                    dphi_talloc_mavg               , &  !< leaf-soil water potential gradient (MPa)  (averaging time-scale: tau_allocation)
                                    growth_cn_talloc_mavg          , &  !< plant growth cn ratio, equivalent to 1/growth_req_n
                                    growth_cp_talloc_mavg          , &  !< plant growth cp ratio, equivalent to 1/(growthreq_n*growth_req_p)
                                    labile_carbon_talloc_mavg      , &  !< moving average of labile carbon pool neeeded for allometry calculations
                                    labile_nitrogen_talloc_mavg    , &  !< moving average of labile nitrogen pool neeeded for allometry calculations
                                    labile_phosphorus_talloc_mavg       !< moving average of labile phosphorus pool neeeded for allometry calculations
    ! 1.7 memory over root life time 
    TYPE(t_jsb_var_real1d)          :: &   
                                    leaf2root_troot_mavg           , &  !< leaf to root mass ratio
                                    unit_npp_troot_mavg            , &  !< rate of change in NPP with  increment in leaf mass, mumol C mol-1
                                    unit_uptake_n_troot_mavg       , &  !< rate of N uptkae per unit root, mumol N mol-1 (time-scale of accumulation: tau_root)
                                    unit_uptake_p_troot_mavg       , &  !< rate of P uptake per unit root, mumol P mol-1 (time-scale of accumulation: tau_root)
                                    fmaint_rate_troot_mavg              !< moving average of N-specific maintentance respiration (~µmol C/mol N / s ) 
    ! 1.8 memory over timescale of vegetation dynamics 
    TYPE(t_jsb_var_real1d)          :: &   
                                    an_boc_tvegdyn_mavg                 !< net photosynthesis at bottom of canopy (lowest canopy layer existing)
    ! 1.6 memory of certain variables, temperature, ...
    TYPE(t_jsb_var_real1d)          :: &   
                                    t_air_week_mavg                , &  !< avrg air temperature across 7 days
                                    t_air_month_mavg               , &  !< avrg air temperature across 30 days
                                    t_jmax_opt_mavg                     !< optimim temperature for electron transport of photosynthesis

    !> structures
    !! 
    TYPE(t_slm_struct_veg_element)  :: &
                                    establishment   , &   !< flux between seed bed and reserve (mol m-2 timestep-1)
                                    reserve_use           !< net flux between reserve and labile pool (mol m-2 timestep-1)

    TYPE(t_slm_struct_veg_pool)     :: &   ! nested pool TYPE t_jsb_n_struct
                                    pool            , &   !< pools of the plant
                                    growth          , &   !< fluxes: growth
                                    litterfall            !< fluxes: litterfall


  CONTAINS
    PROCEDURE :: Init => Init_veg_memory
  END TYPE t_veg_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class'


CONTAINS

  !-----------------------------------------------------------------------------------------------------!
  ! t_slm_struct_veg_element operator overload                                                                                                  !
  !-----------------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------------
  !> + (t_slm_struct_veg_element + t_slm_struct_veg_element)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_veg_element), INTENT(in) :: this, struct2  !< t_slm_struct_veg_element in
    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3        !< t_slm_struct_veg_element out
    ! 0.2 Local
!    INTEGER                                     :: n_elements, & ! number of elements in t_slm_struct_veg_element
!                                                   i_elements    ! iterate over elements

!    TYPE(t_jsb_var_real1d), POINTER             :: struct3var, &
!                                                   thisvar   , &
!                                                   struct2var

    ALLOCATE(struct3)!, source=this)
!    ALLOCATE(struct3%vars(lbound(this%vars,1):ubound(this%vars,1)))
!    n_elements = this%no_of_vars

!    ! loop over elements
!    DO i_elements = 1, n_elements
!      SELECT TYPE (mm => struct3%vars(i_elements)%p)
!        TYPE IS(t_jsb_var_real1d)
!          struct3var => mm
!      END SELECT
!      SELECT TYPE (mm => this%vars(i_elements)%p)
!        TYPE IS(t_jsb_var_real1d)
!          thisvar => mm
!      END SELECT
!      SELECT TYPE (mm => struct2%vars(i_elements)%p)
!        TYPE IS(t_jsb_var_real1d)
!          struct2var => mm
!      END SELECT

!      struct3var = thisvar + struct2var ! struct3%vars(i_elements)%p = this%vars(i_elements)%p + struct2%vars(i_elements)%p
!    END DO

    struct3%carbon       = this%carbon      + struct2%carbon
    struct3%nitrogen     = this%nitrogen    + struct2%nitrogen
    struct3%phosphorus   = this%phosphorus  + struct2%phosphorus
    struct3%carbon13     = this%carbon13    + struct2%carbon13
    struct3%carbon14     = this%carbon14    + struct2%carbon14
    struct3%nitrogen15   = this%nitrogen15  + struct2%nitrogen15
  END FUNCTION Add_struct

  !-----------------------------------------------------------------------------------------------------
  !> + (t_slm_struct_veg_element + scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_struct_and_scalar(this,scalar) RESULT(struct3)

    CLASS(t_slm_struct_veg_element), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this)

    struct3%carbon       = this%carbon      + scalar
    struct3%nitrogen     = this%nitrogen    + scalar
    struct3%phosphorus   = this%phosphorus  + scalar
    struct3%carbon13     = this%carbon13    + scalar
    struct3%carbon14     = this%carbon14    + scalar
    struct3%nitrogen15   = this%nitrogen15  + scalar
  END FUNCTION Add_struct_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> + (scalar + t_slm_struct_veg_element)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_scalar_and_struct(scalar,this2) RESULT(struct3)

    CLASS(t_slm_struct_veg_element), INTENT(in) :: this2
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this2)

    struct3%carbon       = scalar + this2%carbon
    struct3%nitrogen     = scalar + this2%nitrogen
    struct3%phosphorus   = scalar + this2%phosphorus
    struct3%carbon13     = scalar + this2%carbon13
    struct3%carbon14     = scalar + this2%carbon14
    struct3%nitrogen15   = scalar + this2%nitrogen15
  END FUNCTION Add_scalar_and_struct

! *slmdev -> works only with 2D variables (t_jsb_var_real2d)
!  !-----------------------------------------------------------------------------------------------------
!  !> + (t_slm_struct_veg_element + field2d)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Add_struct_and_field2d(this,field2d) RESULT(struct3)

!    CLASS(t_slm_struct_veg_element), INTENT(in) :: this
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this)

!    struct3%carbon       = this%carbon      + field2d
!    struct3%nitrogen     = this%nitrogen    + field2d
!    struct3%phosphorus   = this%phosphorus  + field2d
!    struct3%carbon13     = this%carbon13    + field2d
!    struct3%carbon14     = this%carbon14    + field2d
!    struct3%nitrogen15   = this%nitrogen15  + field2d
!  END FUNCTION Add_struct_and_field2d

!  !-----------------------------------------------------------------------------------------------------
!  !> + (field2d + t_slm_struct_veg_element)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Add_field2d_and_struct(field2d,this2) RESULT(struct3)

!    CLASS(t_slm_struct_veg_element), INTENT(in) :: this2
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this2)

!    struct3%carbon       = field2d + this2%carbon
!    struct3%nitrogen     = field2d + this2%nitrogen
!    struct3%phosphorus   = field2d + this2%phosphorus
!    struct3%carbon13     = field2d + this2%carbon13
!    struct3%carbon14     = field2d + this2%carbon14
!    struct3%nitrogen15   = field2d + this2%nitrogen15
!  END FUNCTION Add_field2d_and_struct

  !-----------------------------------------------------------------------------------------------------
  !> - (t_slm_struct_veg_element - t_slm_struct_veg_element)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_veg_element), INTENT(in) :: this, struct2  !< t_slm_struct_veg_element in
    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3        !< t_slm_struct_veg_element out

    ALLOCATE(struct3)!, source=this)

    struct3%carbon       = this%carbon      - struct2%carbon
    struct3%nitrogen     = this%nitrogen    - struct2%nitrogen
    struct3%phosphorus   = this%phosphorus  - struct2%phosphorus
    struct3%carbon13     = this%carbon13    - struct2%carbon13
    struct3%carbon14     = this%carbon14    - struct2%carbon14
    struct3%nitrogen15   = this%nitrogen15  - struct2%nitrogen15
  END FUNCTION Subtract_struct

  !-----------------------------------------------------------------------------------------------------
  !> - (t_slm_struct_veg_element - scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_struct_and_scalar(this,scalar) RESULT(struct3)

    CLASS(t_slm_struct_veg_element), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this)

    struct3%carbon       = this%carbon      - scalar
    struct3%nitrogen     = this%nitrogen    - scalar
    struct3%phosphorus   = this%phosphorus  - scalar
    struct3%carbon13     = this%carbon13    - scalar
    struct3%carbon14     = this%carbon14    - scalar
    struct3%nitrogen15   = this%nitrogen15  - scalar
  END FUNCTION Subtract_struct_and_scalar

!> @todo operator overload for "scalar - t_slm_struct_veg_element" may not work with 2D variables
!! *slmtmp
  !-----------------------------------------------------------------------------------------------------
  !> - (scalar - t_slm_struct_veg_element)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_scalar_and_struct(scalar,this2) RESULT(struct3)

    CLASS(t_slm_struct_veg_element), INTENT(in) :: this2
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this2)

    struct3%carbon       = scalar - this2%carbon
    struct3%nitrogen     = scalar - this2%nitrogen
    struct3%phosphorus   = scalar - this2%phosphorus
    struct3%carbon13     = scalar - this2%carbon13
    struct3%carbon14     = scalar - this2%carbon14
    struct3%nitrogen15   = scalar - this2%nitrogen15
  END FUNCTION Subtract_scalar_and_struct

! *slmdev -> works only with 2D variables (t_jsb_var_real2d)
!  !-----------------------------------------------------------------------------------------------------
!  !> - (t_slm_struct_veg_element - field2d)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Subtract_struct_and_field2d(this,field2d) RESULT(struct3)

!    CLASS(t_slm_struct_veg_element), INTENT(in) :: this
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this)

!    struct3%carbon       = this%carbon      - field2d
!    struct3%nitrogen     = this%nitrogen    - field2d
!    struct3%phosphorus   = this%phosphorus  - field2d
!    struct3%carbon13     = this%carbon13    - field2d
!    struct3%carbon14     = this%carbon14    - field2d
!    struct3%nitrogen15   = this%nitrogen15  - field2d
!  END FUNCTION Subtract_struct_and_field2d

!  !-----------------------------------------------------------------------------------------------------
!  !> - (field2d - t_slm_struct_veg_element)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Subtract_field2d_and_struct(field2d,this2) RESULT(struct3)

!    CLASS(t_slm_struct_veg_element), INTENT(in) :: this2
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this2)

!    struct3%carbon       = field2d - this2%carbon
!    struct3%nitrogen     = field2d - this2%nitrogen
!    struct3%phosphorus   = field2d - this2%phosphorus
!    struct3%carbon13     = field2d - this2%carbon13
!    struct3%carbon14     = field2d - this2%carbon14
!    struct3%nitrogen15   = field2d - this2%nitrogen15
!  END FUNCTION Subtract_field2d_and_struct

  !-----------------------------------------------------------------------------------------------------
  !> * (t_slm_struct_veg_element * t_slm_struct_veg_element)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_veg_element), INTENT(in) :: this, struct2  !< t_slm_struct_veg_element in
    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3        !< t_slm_struct_veg_element out

    ALLOCATE(struct3)!, source=this)

    struct3%carbon       = this%carbon      * struct2%carbon
    struct3%nitrogen     = this%nitrogen    * struct2%nitrogen
    struct3%phosphorus   = this%phosphorus  * struct2%phosphorus
    struct3%carbon13     = this%carbon13    * struct2%carbon13
    struct3%carbon14     = this%carbon14    * struct2%carbon14
    struct3%nitrogen15   = this%nitrogen15  * struct2%nitrogen15
  END FUNCTION Multiply_struct

  !-----------------------------------------------------------------------------------------------------
  !> * (t_slm_struct_veg_element * scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_struct_and_scalar(this,scalar) RESULT(struct3)

    CLASS(t_slm_struct_veg_element), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this)

    struct3%carbon       = this%carbon      * scalar
    struct3%nitrogen     = this%nitrogen    * scalar
    struct3%phosphorus   = this%phosphorus  * scalar
    struct3%carbon13     = this%carbon13    * scalar
    struct3%carbon14     = this%carbon14    * scalar
    struct3%nitrogen15   = this%nitrogen15  * scalar
  END FUNCTION Multiply_struct_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> * (scalar * t_slm_struct_veg_element)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_scalar_and_struct(scalar,this2) RESULT(struct3)

    CLASS(t_slm_struct_veg_element), INTENT(in) :: this2
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this2)

    struct3%carbon       = scalar * this2%carbon
    struct3%nitrogen     = scalar * this2%nitrogen
    struct3%phosphorus   = scalar * this2%phosphorus
    struct3%carbon13     = scalar * this2%carbon13
    struct3%carbon14     = scalar * this2%carbon14
    struct3%nitrogen15   = scalar * this2%nitrogen15
  END FUNCTION Multiply_scalar_and_struct

! *slmdev -> works only with 2D variables (t_jsb_var_real2d)
!  !-----------------------------------------------------------------------------------------------------
!  !> * (t_slm_struct_veg_element * field2d)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Multiply_struct_and_field2d(this,field2d) RESULT(struct3)

!    CLASS(t_slm_struct_veg_element), INTENT(in) :: this
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this)

!    struct3%carbon       = this%carbon      * field2d
!    struct3%nitrogen     = this%nitrogen    * field2d
!    struct3%phosphorus   = this%phosphorus  * field2d
!    struct3%carbon13     = this%carbon13    * field2d
!    struct3%carbon14     = this%carbon14    * field2d
!    struct3%nitrogen15   = this%nitrogen15  * field2d
!  END FUNCTION Multiply_struct_and_field2d

!  !-----------------------------------------------------------------------------------------------------
!  !> * (field2d * t_slm_struct_veg_element)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Multiply_field2d_and_struct(field2d,this2) RESULT(struct3)

!    CLASS(t_slm_struct_veg_element), INTENT(in) :: this2
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this2)

!    struct3%carbon       = field2d * this2%carbon
!    struct3%nitrogen     = field2d * this2%nitrogen
!    struct3%phosphorus   = field2d * this2%phosphorus
!    struct3%carbon13     = field2d * this2%carbon13
!    struct3%carbon14     = field2d * this2%carbon14
!    struct3%nitrogen15   = field2d * this2%nitrogen15
!  END FUNCTION Multiply_field2d_and_struct

  !-----------------------------------------------------------------------------------------------------
  !> / (t_slm_struct_veg_element / t_slm_struct_veg_element)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_veg_element), INTENT(in) :: this, struct2  !< t_slm_struct_veg_element in
    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3        !< t_slm_struct_veg_element out

    ALLOCATE(struct3)!, source=this)

    struct3%carbon       = this%carbon      / struct2%carbon
    struct3%nitrogen     = this%nitrogen    / struct2%nitrogen
    struct3%phosphorus   = this%phosphorus  / struct2%phosphorus
    struct3%carbon13     = this%carbon13    / struct2%carbon13
    struct3%carbon14     = this%carbon14    / struct2%carbon14
    struct3%nitrogen15   = this%nitrogen15  / struct2%nitrogen15
  END FUNCTION Divide_struct

  !-----------------------------------------------------------------------------------------------------
  !> / (t_slm_struct_veg_element / scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_struct_and_scalar(this,scalar) RESULT(struct3)

    CLASS(t_slm_struct_veg_element), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this)

    struct3%carbon       = this%carbon      / scalar
    struct3%nitrogen     = this%nitrogen    / scalar
    struct3%phosphorus   = this%phosphorus  / scalar
    struct3%carbon13     = this%carbon13    / scalar
    struct3%carbon14     = this%carbon14    / scalar
    struct3%nitrogen15   = this%nitrogen15  / scalar
  END FUNCTION Divide_struct_and_scalar

!> @todo operator overload for "scalar / t_slm_struct_veg_element" may not work with 2D variables
!! *slmtmp
  !-----------------------------------------------------------------------------------------------------
  !> / (scalar / t_slm_struct_veg_element)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_scalar_and_struct(scalar,this2) RESULT(struct3)

    CLASS(t_slm_struct_veg_element), INTENT(in) :: this2
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this2)

    struct3%carbon       = scalar / this2%carbon
    struct3%nitrogen     = scalar / this2%nitrogen
    struct3%phosphorus   = scalar / this2%phosphorus
    struct3%carbon13     = scalar / this2%carbon13
    struct3%carbon14     = scalar / this2%carbon14
    struct3%nitrogen15   = scalar / this2%nitrogen15
  END FUNCTION Divide_scalar_and_struct

! *slmdev -> works only with 2D variables (t_jsb_var_real2d)
!  !-----------------------------------------------------------------------------------------------------
!  !> / (t_slm_struct_veg_element / field2d)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Divide_struct_and_field2d(this,field2d) RESULT(struct3)

!    CLASS(t_slm_struct_veg_element), INTENT(in) :: this
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this)

!    struct3%carbon       = this%carbon      / field2d
!    struct3%nitrogen     = this%nitrogen    / field2d
!    struct3%phosphorus   = this%phosphorus  / field2d
!    struct3%carbon13     = this%carbon13    / field2d
!    struct3%carbon14     = this%carbon14    / field2d
!    struct3%nitrogen15   = this%nitrogen15  / field2d
!  END FUNCTION Divide_struct_and_field2d

!  !-----------------------------------------------------------------------------------------------------
!  !> / (field2d / t_slm_struct_veg_element)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Divide_field2d_and_struct(field2d,this2) RESULT(struct3)

!    CLASS(t_slm_struct_veg_element), INTENT(in) :: this2
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_element), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this2)

!    struct3%carbon       = field2d / this2%carbon
!    struct3%nitrogen     = field2d / this2%nitrogen
!    struct3%phosphorus   = field2d / this2%phosphorus
!    struct3%carbon13     = field2d / this2%carbon13
!    struct3%carbon14     = field2d / this2%carbon14
!    struct3%nitrogen15   = field2d / this2%nitrogen15
!  END FUNCTION Divide_field2d_and_struct


  !-----------------------------------------------------------------------------------------------------!
  ! t_slm_struct_veg_pool operator overload                                                                                                  !
  !-----------------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------------
  !> + (t_slm_struct_veg_pool + t_slm_struct_veg_pool)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_n_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this, struct2  !< t_slm_struct_veg_pool in
    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3        !< t_slm_struct_veg_pool out

    ALLOCATE(struct3)!, source=this)

    struct3%leaf          = this%leaf         + struct2%leaf
    struct3%fine_root     = this%fine_root    + struct2%fine_root
    struct3%coarse_root   = this%coarse_root  + struct2%coarse_root
    struct3%sap_wood      = this%sap_wood     + struct2%sap_wood
    struct3%heart_wood    = this%heart_wood   + struct2%heart_wood
    struct3%labile        = this%labile       + struct2%labile
    struct3%reserve       = this%reserve      + struct2%reserve
    struct3%fruit         = this%fruit        + struct2%fruit
    struct3%seed_bed      = this%seed_bed     + struct2%seed_bed
  END FUNCTION Add_n_struct

  !-----------------------------------------------------------------------------------------------------
  !> + (t_slm_struct_veg_pool + scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_n_struct_and_scalar(this,scalar) RESULT(struct3)

    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this
    REAL(wp),                     INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this)

    struct3%leaf          = this%leaf         + scalar
    struct3%fine_root     = this%fine_root    + scalar
    struct3%coarse_root   = this%coarse_root  + scalar
    struct3%sap_wood      = this%sap_wood     + scalar
    struct3%heart_wood    = this%heart_wood   + scalar
    struct3%labile        = this%labile       + scalar
    struct3%reserve       = this%reserve      + scalar
    struct3%fruit         = this%fruit        + scalar
    struct3%seed_bed      = this%seed_bed     + scalar
  END FUNCTION Add_n_struct_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> + (scalar + t_slm_struct_veg_pool)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_scalar_and_n_struct(scalar,this2) RESULT(struct3)

    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this2
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this2)

    struct3%leaf          = scalar + this2%leaf
    struct3%fine_root     = scalar + this2%fine_root
    struct3%coarse_root   = scalar + this2%coarse_root
    struct3%sap_wood      = scalar + this2%sap_wood
    struct3%heart_wood    = scalar + this2%heart_wood
    struct3%labile        = scalar + this2%labile
    struct3%reserve       = scalar + this2%reserve
    struct3%fruit         = scalar + this2%fruit
    struct3%seed_bed      = scalar + this2%seed_bed
  END FUNCTION Add_scalar_and_n_struct

! *slmdev -> works only with 2D variables (t_jsb_var_real2d)
!  !-----------------------------------------------------------------------------------------------------
!  !> + (t_slm_struct_veg_pool + field2d)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Add_n_struct_and_field2d(this,field2d) RESULT(struct3)

!    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this)

!    struct3%leaf           = this%leaf         + field2d
!    struct3%fine_root      = this%fine_root    + field2d
!    struct3%coarse_root    = this%coarse_root  + field2d
!    struct3%sap_wood       = this%sap_wood     + field2d
!    struct3%heart_wood     = this%heart_wood   + field2d
!    struct3%labile         = this%labile       + field2d
!    struct3%reserve        = this%reserve      + field2d
!    struct3%fruit          = this%fruit        + field2d
!    struct3%seed_bed       = this%seed_bed     + field2d
!  END FUNCTION Add_n_struct_and_field2d

!  !-----------------------------------------------------------------------------------------------------
!  !> + (field2d + t_slm_struct_veg_pool)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Add_field2d_and_n_struct(field2d,this2) RESULT(struct3)

!    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this2
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this2)

!    struct3%leaf           = field2d + this2%leaf
!    struct3%fine_root      = field2d + this2%fine_root
!    struct3%coarse_root    = field2d + this2%coarse_root
!    struct3%sap_wood       = field2d + this2%sap_wood
!    struct3%heart_wood     = field2d + this2%heart_wood
!    struct3%labile         = field2d + this2%labile
!    struct3%reserve        = field2d + this2%reserve
!    struct3%fruit          = field2d + this2%fruit
!    struct3%seed_bed       = field2d + this2%seed_bed
!  END FUNCTION Add_field2d_and_n_struct

  !-----------------------------------------------------------------------------------------------------
  !> - (t_slm_struct_veg_pool - t_slm_struct_veg_pool)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_n_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this, struct2  !< t_slm_struct_veg_pool in
    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3        !< t_slm_struct_veg_pool out

    ALLOCATE(struct3)!, source=this)

    struct3%leaf          = this%leaf         - struct2%leaf
    struct3%fine_root     = this%fine_root    - struct2%fine_root
    struct3%coarse_root   = this%coarse_root  - struct2%coarse_root
    struct3%sap_wood      = this%sap_wood     - struct2%sap_wood
    struct3%heart_wood    = this%heart_wood   - struct2%heart_wood
    struct3%labile        = this%labile       - struct2%labile
    struct3%reserve       = this%reserve      - struct2%reserve
    struct3%fruit         = this%fruit        - struct2%fruit
    struct3%seed_bed      = this%seed_bed     - struct2%seed_bed
  END FUNCTION Subtract_n_struct

  !-----------------------------------------------------------------------------------------------------
  !> - (t_slm_struct_veg_pool - scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_n_struct_and_scalar(this,scalar) RESULT(struct3)

    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this)

    struct3%leaf          = this%leaf         - scalar
    struct3%fine_root     = this%fine_root    - scalar
    struct3%coarse_root   = this%coarse_root  - scalar
    struct3%sap_wood      = this%sap_wood     - scalar
    struct3%heart_wood    = this%heart_wood   - scalar
    struct3%labile        = this%labile       - scalar
    struct3%reserve       = this%reserve      - scalar
    struct3%fruit         = this%fruit        - scalar
    struct3%seed_bed      = this%seed_bed     - scalar
  END FUNCTION Subtract_n_struct_and_scalar

!> @todo operator overload for "scalar - t_slm_struct_veg_pool" may not work with 2D variables
!! *slmtmp
  !-----------------------------------------------------------------------------------------------------
  !> - (scalar - t_slm_struct_veg_pool)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_scalar_and_n_struct(scalar,this2) RESULT(struct3)

    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this2
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this2)

    struct3%leaf          = scalar - this2%leaf
    struct3%fine_root     = scalar - this2%fine_root
    struct3%coarse_root   = scalar - this2%coarse_root
    struct3%sap_wood      = scalar - this2%sap_wood
    struct3%heart_wood    = scalar - this2%heart_wood
    struct3%labile        = scalar - this2%labile
    struct3%reserve       = scalar - this2%reserve
    struct3%fruit         = scalar - this2%fruit
    struct3%seed_bed      = scalar - this2%seed_bed
  END FUNCTION Subtract_scalar_and_n_struct

! *slmdev -> works only with 2D variables (t_jsb_var_real2d)
!  !-----------------------------------------------------------------------------------------------------
!  !> - (t_slm_struct_veg_pool - field2d)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Subtract_n_struct_and_field2d(this,field2d) RESULT(struct3)

!    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this)

!    struct3%leaf           = this%leaf         - field2d
!    struct3%fine_root      = this%fine_root    - field2d
!    struct3%coarse_root    = this%coarse_root  - field2d
!    struct3%sap_wood       = this%sap_wood     - field2d
!    struct3%heart_wood     = this%heart_wood   - field2d
!    struct3%labile         = this%labile       - field2d
!    struct3%reserve        = this%reserve      - field2d
!    struct3%fruit          = this%fruit        - field2d
!    struct3%seed_bed       = this%seed_bed     - field2d
!  END FUNCTION Subtract_n_struct_and_field2d

!  !-----------------------------------------------------------------------------------------------------
!  !> - (field2d - t_slm_struct_veg_pool)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Subtract_field2d_and_n_struct(field2d,this2) RESULT(struct3)

!    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this2
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this2)

!    struct3%leaf           = field2d - this2%leaf
!    struct3%fine_root      = field2d - this2%fine_root
!    struct3%coarse_root    = field2d - this2%coarse_root
!    struct3%sap_wood       = field2d - this2%sap_wood
!    struct3%heart_wood     = field2d - this2%heart_wood
!    struct3%labile         = field2d - this2%labile
!    struct3%reserve        = field2d - this2%reserve
!    struct3%fruit          = field2d - this2%fruit
!    struct3%seed_bed       = field2d - this2%seed_bed
!  END FUNCTION Subtract_field2d_and_n_struct

  !-----------------------------------------------------------------------------------------------------
  !> * (t_slm_struct_veg_pool * t_slm_struct_veg_pool)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_n_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this, struct2  !< t_slm_struct_veg_pool in
    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3        !< t_slm_struct_veg_pool out

    ALLOCATE(struct3)!, source=this)

    struct3%leaf          = this%leaf         * struct2%leaf
    struct3%fine_root     = this%fine_root    * struct2%fine_root
    struct3%coarse_root   = this%coarse_root  * struct2%coarse_root
    struct3%sap_wood      = this%sap_wood     * struct2%sap_wood
    struct3%heart_wood    = this%heart_wood   * struct2%heart_wood
    struct3%labile        = this%labile       * struct2%labile
    struct3%reserve       = this%reserve      * struct2%reserve
    struct3%fruit         = this%fruit        * struct2%fruit
    struct3%seed_bed      = this%seed_bed     * struct2%seed_bed
  END FUNCTION Multiply_n_struct

  !-----------------------------------------------------------------------------------------------------
  !> * (t_slm_struct_veg_pool * scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_n_struct_and_scalar(this,scalar) RESULT(struct3)

    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this)

    struct3%leaf          = this%leaf      * scalar
    struct3%fine_root     = this%fine_root    * scalar
    struct3%coarse_root   = this%coarse_root  * scalar
    struct3%sap_wood      = this%sap_wood    * scalar
    struct3%heart_wood    = this%heart_wood    * scalar
    struct3%labile        = this%labile       * scalar
    struct3%reserve       = this%reserve      * scalar
    struct3%fruit         = this%fruit        * scalar
    struct3%seed_bed      = this%seed_bed     * scalar
  END FUNCTION Multiply_n_struct_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> * (scalar * t_slm_struct_veg_pool)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_scalar_and_n_struct(scalar,this2) RESULT(struct3)

    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this2
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this2)

    struct3%leaf          = scalar * this2%leaf
    struct3%fine_root     = scalar * this2%fine_root
    struct3%coarse_root   = scalar * this2%coarse_root
    struct3%sap_wood      = scalar * this2%sap_wood
    struct3%heart_wood    = scalar * this2%heart_wood
    struct3%labile        = scalar * this2%labile
    struct3%reserve       = scalar * this2%reserve
    struct3%fruit         = scalar * this2%fruit
    struct3%seed_bed      = scalar * this2%seed_bed
  END FUNCTION Multiply_scalar_and_n_struct

! *slmdev -> works only with 2D variables (t_jsb_var_real2d)
!  !-----------------------------------------------------------------------------------------------------
!  !> * (t_slm_struct_veg_pool * field2d)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Multiply_n_struct_and_field2d(this,field2d) RESULT(struct3)

!    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this)

!    struct3%leaf           = this%leaf         * field2d
!    struct3%fine_root      = this%fine_root    * field2d
!    struct3%coarse_root    = this%coarse_root  * field2d
!    struct3%sap_wood       = this%sap_wood     * field2d
!    struct3%heart_wood     = this%heart_wood   * field2d
!    struct3%labile         = this%labile       * field2d
!    struct3%reserve        = this%reserve      * field2d
!    struct3%fruit          = this%fruit        * field2d
!    struct3%seed_bed       = this%seed_bed     * field2d
!  END FUNCTION Multiply_n_struct_and_field2d

!  !-----------------------------------------------------------------------------------------------------
!  !> * (field2d * t_slm_struct_veg_pool)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Multiply_field2d_and_n_struct(field2d,this2) RESULT(struct3)

!    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this2
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this2)

!    struct3%leaf           = field2d * this2%leaf
!    struct3%fine_root      = field2d * this2%fine_root
!    struct3%coarse_root    = field2d * this2%coarse_root
!    struct3%sap_wood       = field2d * this2%sap_wood
!    struct3%heart_wood     = field2d * this2%heart_wood
!    struct3%labile         = field2d * this2%labile
!    struct3%reserve        = field2d * this2%reserve
!    struct3%fruit          = field2d * this2%fruit
!    struct3%seed_bed       = field2d * this2%seed_bed
!  END FUNCTION Multiply_field2d_and_n_struct

  !-----------------------------------------------------------------------------------------------------
  !> / (t_slm_struct_veg_pool / t_slm_struct_veg_pool)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_n_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this, struct2  !< t_slm_struct_veg_pool in
    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3        !< t_slm_struct_veg_pool out

    ALLOCATE(struct3)!, source=this)

    struct3%leaf          = this%leaf         / struct2%leaf
    struct3%fine_root     = this%fine_root    / struct2%fine_root
    struct3%coarse_root   = this%coarse_root  / struct2%coarse_root
    struct3%sap_wood      = this%sap_wood     / struct2%sap_wood
    struct3%heart_wood    = this%heart_wood   / struct2%heart_wood
    struct3%labile        = this%labile       / struct2%labile
    struct3%reserve       = this%reserve      / struct2%reserve
    struct3%fruit         = this%fruit        / struct2%fruit
    struct3%seed_bed      = this%seed_bed     / struct2%seed_bed
  END FUNCTION Divide_n_struct

  !-----------------------------------------------------------------------------------------------------
  !> / (t_slm_struct_veg_pool / scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_n_struct_and_scalar(this,scalar) RESULT(struct3)

    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this)

    struct3%leaf          = this%leaf         / scalar
    struct3%fine_root     = this%fine_root    / scalar
    struct3%coarse_root   = this%coarse_root  / scalar
    struct3%sap_wood      = this%sap_wood     / scalar
    struct3%heart_wood    = this%heart_wood   / scalar
    struct3%labile        = this%labile       / scalar
    struct3%reserve       = this%reserve      / scalar
    struct3%fruit         = this%fruit        / scalar
    struct3%seed_bed      = this%seed_bed     / scalar
  END FUNCTION Divide_n_struct_and_scalar

!> @todo operator overload for "scalar / t_slm_struct_veg_pool" may not work with 2D variables
!! *slmtmp
  !-----------------------------------------------------------------------------------------------------
  !> / (scalar / t_slm_struct_veg_pool)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_scalar_and_n_struct(scalar,this2) RESULT(struct3)

    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this2
    REAL(wp),                     INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3)!, source=this2)

    struct3%leaf          = scalar / this2%leaf
    struct3%fine_root     = scalar / this2%fine_root
    struct3%coarse_root   = scalar / this2%coarse_root
    struct3%sap_wood      = scalar / this2%sap_wood
    struct3%heart_wood    = scalar / this2%heart_wood
    struct3%labile        = scalar / this2%labile
    struct3%reserve       = scalar / this2%reserve
    struct3%fruit         = scalar / this2%fruit
    struct3%seed_bed      = scalar / this2%seed_bed
  END FUNCTION Divide_scalar_and_n_struct

! *slmdev -> works only with 2D variables (t_jsb_var_real2d)
!  !-----------------------------------------------------------------------------------------------------
!  !> / (t_slm_struct_veg_pool / field2d)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Divide_n_struct_and_field2d(this,field2d) RESULT(struct3)

!    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this)

!    struct3%leaf           = this%leaf         / field2d
!    struct3%fine_root      = this%fine_root    / field2d
!    struct3%coarse_root    = this%coarse_root  / field2d
!    struct3%sap_wood       = this%sap_wood     / field2d
!    struct3%heart_wood     = this%heart_wood   / field2d
!    struct3%labile         = this%labile       / field2d
!    struct3%reserve        = this%reserve      / field2d
!    struct3%fruit          = this%fruit        / field2d
!    struct3%seed_bed       = this%seed_bed     / field2d
!  END FUNCTION Divide_n_struct_and_field2d

!  !-----------------------------------------------------------------------------------------------------
!  !> / (field2d / t_slm_struct_veg_pool)
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Divide_field2d_and_n_struct(field2d,this2) RESULT(struct3)

!    CLASS(t_slm_struct_veg_pool), INTENT(in) :: this2
!    REAL(wp),                        INTENT(in) :: field2d(:,:)
!    TYPE(t_slm_struct_veg_pool), ALLOCATABLE :: struct3

!    ALLOCATE(struct3)!, source=this2)

!    struct3%leaf           = field2d / this2%leaf
!    struct3%fine_root      = field2d / this2%fine_root
!    struct3%coarse_root    = field2d / this2%coarse_root
!    struct3%sap_wood       = field2d / this2%sap_wood
!    struct3%heart_wood     = field2d / this2%heart_wood
!    struct3%labile         = field2d / this2%labile
!    struct3%reserve        = field2d / this2%reserve
!    struct3%fruit          = field2d / this2%fruit
!    struct3%seed_bed       = field2d / this2%seed_bed
!  END FUNCTION Divide_field2d_and_n_struct



  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the processes in: vegetation
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_veg_memory(mem, prefix, suffix, lct_ids, model_id, ics, ice, iblk)

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
!    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
!    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_atmland_constants,   ONLY: standard_press_srf, def_co2_mixing_ratio

    CLASS(t_veg_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids !(:)
    INTEGER,             INTENT(in)    :: model_id
    INTEGER,             INTENT(in)    :: ics   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)    :: ice   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)    :: iblk  ! model%options -> for SLM compatibility

!    TYPE(t_jsb_model), POINTER :: model
!    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
!    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer        ! Vertical grid
!    INTEGER :: table

!    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_assimi_memory'

!    model => Get_model(model_id)

!    table = tables(1)

!    hgrid   => Get_grid(mem%grid_id)
!    surface => Get_vgrid('surface')
!    canopy_layer => Get_vgrid('canopy_layer')

    IF (.NOT. ASSOCIATED(mem%growing_season%ptr)) THEN

      !> add variables to veg memory veg\%VAR
      !! 
      !! slm_add_var -> slm_add_var_list_element_r1d    *slmtemp 1D -> 2D
      !!                slm_add_var_list_element_l1d
      !!  
      CALL slm_add_var(mem, 'veg_growing_season', mem%growing_season, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=.TRUE.)     ! initvalue 

      CALL slm_add_var(mem, 'height', mem%height, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'diameter', mem%diameter, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'dens_ind', mem%dens_ind, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'lai', mem%lai, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'target_cn_leaf', mem%target_cn_leaf, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'target_np_leaf', mem%target_np_leaf, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'mean_leaf_age', mem%mean_leaf_age, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gc', mem%gc, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'f_n_demand', mem%f_n_demand, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'f_p_demand', mem%f_p_demand, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air', mem%beta_air, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps', mem%beta_soil_ps, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'dphi', mem%dphi, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs', mem%beta_soil_gs, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_req_n', mem%growth_req_n, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_req_p', mem%growth_req_p, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gdd', mem%gdd, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'nd_dormance', mem%nd_dormance, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_jmax_opt', mem%t_jmax_opt, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'an_boc', mem%an_boc, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'leaf_cn_direction', mem%leaf_cn_direction, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'apar', mem%apar, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'nl_tot_cl', mem%nl_tot_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fnl_chl_cl', mem%fnl_chl_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fnl_et_cl', mem%fnl_et_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fnl_rub_cl', mem%fnl_rub_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fnl_pepc_cl', mem%fnl_pepc_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fnl_oth_cl', mem%fnl_oth_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_sunlit_cl', mem%par_sunlit_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_shaded_cl', mem%par_shaded_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fleaf_sunlit_cl', mem%fleaf_sunlit_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'lai_cl', mem%lai_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'cumm_lai_cl', mem%cumm_lai_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gpp', mem%gpp, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'an', mem%an, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gpp_c13', mem%gpp_c13, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gpp_c14', mem%gpp_c14, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'maint_respiration_leaf', mem%maint_respiration_leaf, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'maint_respiration_pot', mem%maint_respiration_pot, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'maint_respiration', mem%maint_respiration, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'maint_respiration_c13', mem%maint_respiration_c13, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'maint_respiration_c14', mem%maint_respiration_c14, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_respiration', mem%growth_respiration, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_respiration_c13', mem%growth_respiration_c13, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_respiration_c14', mem%growth_respiration_c14, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n_transform_respiration', mem%n_transform_respiration, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n_fixation_respiration', mem%n_fixation_respiration, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n_processing_respiration', mem%n_processing_respiration, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n_processing_respiration_c13', mem%n_processing_respiration_c13, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n_processing_respiration_c14', mem%n_processing_respiration_c14, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'npp', mem%npp, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'npp_c13', mem%npp_c13, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'npp_c14', mem%npp_c14, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'uptake_n', mem%uptake_n, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'uptake_n15', mem%uptake_n15, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'uptake_p', mem%uptake_p, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n_fixation', mem%n_fixation, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n_fixation_n15', mem%n_fixation_n15, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_leaf_n', mem%recycling_leaf_n, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_leaf_n15', mem%recycling_leaf_n15, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_leaf_p', mem%recycling_leaf_p, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_fine_root_n', mem%recycling_fine_root_n, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_fine_root_n15', mem%recycling_fine_root_n15, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_fine_root_p', mem%recycling_fine_root_p, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_sap_wood_n', mem%recycling_sap_wood_n, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_sap_wood_n15', mem%recycling_sap_wood_n15, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_sap_wood_p', mem%recycling_sap_wood_p, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transpiration', mem%transpiration, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fmaint_rate', mem%fmaint_rate, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'ci', mem%ci, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_npp', mem%unit_npp, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_transpiration', mem%unit_transpiration, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_uptake_n', mem%unit_uptake_n, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_uptake_p', mem%unit_uptake_p, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'delta_dens_ind', mem%delta_dens_ind, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'cost_n_uptake_root', mem%cost_n_uptake_root, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'cost_n_fixation', mem%cost_n_fixation, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_sunlit_daytime_cl', mem%par_sunlit_daytime_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_shaded_daytime_cl', mem%par_shaded_daytime_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fleaf_sunlit_daytime_cl', mem%fleaf_sunlit_daytime_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_daytime', mem%t_air_daytime, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 'press_srf_daytime', mem%press_srf_daytime, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=standard_press_srf) ! initvalue
  
      CALL slm_add_var(mem, 'co2_mixing_ratio_daytime', mem%co2_mixing_ratio_daytime, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=def_co2_mixing_ratio) ! initvalue
  
      CALL slm_add_var(mem, 'ga_daytime', mem%ga_daytime, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.2_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air_daytime', mem%beta_air_daytime, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps_daytime', mem%beta_soil_ps_daytime, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_daytime', mem%beta_soil_gs_daytime, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_jmax_opt_daytime', mem%t_jmax_opt_daytime, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=17.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_sunlit_daytime_dacc_cl', mem%par_sunlit_daytime_dacc_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_shaded_daytime_dacc_cl', mem%par_shaded_daytime_dacc_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fleaf_sunlit_daytime_dacc_cl', mem%fleaf_sunlit_daytime_dacc_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_daytime_dacc', mem%t_air_daytime_dacc, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 'press_srf_daytime_dacc', mem%press_srf_daytime_dacc, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=standard_press_srf) ! initvalue
  
      CALL slm_add_var(mem, 'co2_mixing_ratio_daytime_dacc', mem%co2_mixing_ratio_daytime_dacc, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=def_co2_mixing_ratio) ! initvalue
  
      CALL slm_add_var(mem, 'ga_daytime_dacc', mem%ga_daytime_dacc, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.2_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air_daytime_dacc', mem%beta_air_daytime_dacc, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps_daytime_dacc', mem%beta_soil_ps_daytime_dacc, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_daytime_dacc', mem%beta_soil_gs_daytime_dacc, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_jmax_opt_daytime_dacc', mem%t_jmax_opt_daytime_dacc, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'daytime_counter', mem%daytime_counter, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gpp_tlabile_mavg', mem%gpp_tlabile_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'maint_respiration_tlabile_mavg', mem%maint_respiration_tlabile_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_req_n_tlabile_mavg', mem%growth_req_n_tlabile_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1._wp/25._wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_req_p_tlabile_mavg', mem%growth_req_p_tlabile_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1._wp/14._wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_tphen_mavg', mem%t_air_tphen_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_tphen_mavg', mem%beta_soil_gs_tphen_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'npp_tuptake_mavg', mem%npp_tuptake_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'demand_uptake_n_tuptake_mavg', mem%demand_uptake_n_tuptake_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.1_wp) ! initvalue
  
      CALL slm_add_var(mem, 'demand_uptake_p_tuptake_mavg', mem%demand_uptake_p_tuptake_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.1_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_req_n_tuptake_mavg', mem%growth_req_n_tuptake_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1._wp/25._wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_req_p_tuptake_mavg', mem%growth_req_p_tuptake_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1._wp/14._wp) ! initvalue
  
      CALL slm_add_var(mem, 'cost_n_fixation_tuptake_mavg', mem%cost_n_fixation_tuptake_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=10.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_sunlit_tfrac_mavg_cl', mem%par_sunlit_tfrac_mavg_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=150.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_shaded_tfrac_mavg_cl', mem%par_shaded_tfrac_mavg_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=10.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fleaf_sunlit_tfrac_mavg_cl', mem%fleaf_sunlit_tfrac_mavg_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.5_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_tfrac_mavg', mem%t_air_tfrac_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 'press_srf_tfrac_mavg', mem%press_srf_tfrac_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=standard_press_srf) ! initvalue
  
      CALL slm_add_var(mem, 'co2_mixing_ratio_tfrac_mavg', mem%co2_mixing_ratio_tfrac_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=def_co2_mixing_ratio) ! initvalue
  
      CALL slm_add_var(mem, 'ga_tfrac_mavg', mem%ga_tfrac_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.2_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air_tfrac_mavg', mem%beta_air_tfrac_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps_tfrac_mavg', mem%beta_soil_ps_tfrac_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_tfrac_mavg', mem%beta_soil_gs_tfrac_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_sunlit_tcnl_mavg_cl', mem%par_sunlit_tcnl_mavg_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=150.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_shaded_tcnl_mavg_cl', mem%par_shaded_tcnl_mavg_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=10.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fleaf_sunlit_tcnl_mavg_cl', mem%fleaf_sunlit_tcnl_mavg_cl, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.5_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_tcnl_mavg', mem%t_air_tcnl_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 'press_srf_tcnl_mavg', mem%press_srf_tcnl_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=standard_press_srf) ! initvalue
  
      CALL slm_add_var(mem, 'co2_mixing_ratio_tcnl_mavg', mem%co2_mixing_ratio_tcnl_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=def_co2_mixing_ratio) ! initvalue
  
      CALL slm_add_var(mem, 'ga_tcnl_mavg', mem%ga_tcnl_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.2_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air_tcnl_mavg', mem%beta_air_tcnl_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps_tcnl_mavg', mem%beta_soil_ps_tcnl_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_tcnl_mavg', mem%beta_soil_gs_tcnl_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'npp_talloc_mavg', mem%npp_talloc_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transpiration_talloc_mavg', mem%transpiration_talloc_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_npp_talloc_mavg', mem%unit_npp_talloc_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_transpiration_talloc_mavg', mem%unit_transpiration_talloc_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_uptake_n_talloc_mavg', mem%unit_uptake_n_talloc_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_uptake_p_talloc_mavg', mem%unit_uptake_p_talloc_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'dphi_talloc_mavg', mem%dphi_talloc_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=-2.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_cn_talloc_mavg', mem%growth_cn_talloc_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_cp_talloc_mavg', mem%growth_cp_talloc_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'labile_carbon_talloc_mavg', mem%labile_carbon_talloc_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'labile_nitrogen_talloc_mavg', mem%labile_nitrogen_talloc_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'labile_phosphorus_talloc_mavg', mem%labile_phosphorus_talloc_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'leaf2root_troot_mavg', mem%leaf2root_troot_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_npp_troot_mavg', mem%unit_npp_troot_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_uptake_n_troot_mavg', mem%unit_uptake_n_troot_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_uptake_p_troot_mavg', mem%unit_uptake_p_troot_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fmaint_rate_troot_mavg', mem%fmaint_rate_troot_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'an_boc_tvegdyn_mavg', mem%an_boc_tvegdyn_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=5.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_week_mavg', mem%t_air_week_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_month_mavg', mem%t_air_month_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_jmax_opt_mavg', mem%t_jmax_opt_mavg, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=17.0_wp) ! initvalue
      !> simple TYPE (t_jsb_struct)
      !! 
      !! slm_add_var -> slm_add_struct()
      !! 
      CALL slm_add_var(mem, 'veg_establishment',          mem%establishment)
      CALL slm_add_var(mem, 'veg_reserve_use',            mem%reserve_use)
      !> nested TYPE (t_jsb_n_struct)
      !! 
      !! slm_add_var -> slm_add_n_struct()
      !! 
      CALL slm_add_var(mem, 'veg_pool',       mem%pool)
      CALL slm_add_var(mem, 'veg_growth',     mem%growth)
      CALL slm_add_var(mem, 'veg_litterfall', mem%litterfall)
      !> level 2 TYPE (l2_struct) nested within the above nested TYPE structures (t_jsb_n_struct)
      !! 
      !! slm_add_var -> slm_add_l2_struct()
      !! 
      !! veg%pool \n
      !! veg%growth : fluxes  \n
      !! veg%litterfall : fluxes  \n
      CALL slm_add_var(mem%pool, 'veg_pool_leaf',         mem%pool%leaf)       
      CALL slm_add_var(mem%pool, 'veg_pool_fine_root',    mem%pool%fine_root)
      CALL slm_add_var(mem%pool, 'veg_pool_coarse_root',  mem%pool%coarse_root)
      CALL slm_add_var(mem%pool, 'veg_pool_sap_wood',     mem%pool%sap_wood)
      CALL slm_add_var(mem%pool, 'veg_pool_heart_wood',   mem%pool%heart_wood)
      CALL slm_add_var(mem%pool, 'veg_pool_labile',       mem%pool%labile)
      CALL slm_add_var(mem%pool, 'veg_pool_reserve',      mem%pool%reserve)
      CALL slm_add_var(mem%pool, 'veg_pool_fruit',        mem%pool%fruit)
      CALL slm_add_var(mem%pool, 'veg_pool_seed_bed',     mem%pool%seed_bed)
      CALL slm_add_var(mem%growth, 'veg_growth_leaf',         mem%growth%leaf)       
      CALL slm_add_var(mem%growth, 'veg_growth_fine_root',    mem%growth%fine_root)
      CALL slm_add_var(mem%growth, 'veg_growth_coarse_root',  mem%growth%coarse_root)
      CALL slm_add_var(mem%growth, 'veg_growth_sap_wood',     mem%growth%sap_wood)
      CALL slm_add_var(mem%growth, 'veg_growth_heart_wood',   mem%growth%heart_wood)
      CALL slm_add_var(mem%growth, 'veg_growth_labile',       mem%growth%labile)
      CALL slm_add_var(mem%growth, 'veg_growth_reserve',      mem%growth%reserve)
      CALL slm_add_var(mem%growth, 'veg_growth_fruit',        mem%growth%fruit)
      CALL slm_add_var(mem%growth, 'veg_growth_seed_bed',     mem%growth%seed_bed)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_leaf',         mem%litterfall%leaf)       
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_fine_root',    mem%litterfall%fine_root)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_coarse_root',  mem%litterfall%coarse_root)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_sap_wood',     mem%litterfall%sap_wood)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_heart_wood',   mem%litterfall%heart_wood)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_labile',       mem%litterfall%labile)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_reserve',      mem%litterfall%reserve)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_fruit',        mem%litterfall%fruit)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_seed_bed',     mem%litterfall%seed_bed)
      !> add variables to the structs (veg\%t_jsb_struct\%VAR)
      !!
      !! e.g. establishment
      !! slm_add_var -> slm_add_struct_element_r1d    *slmtemp 1D -> 2D
      !! 
      CALL slm_add_var(mem%establishment, 'veg_establishment_carbon', mem%establishment%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_nitrogen', mem%establishment%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_phosphorus', mem%establishment%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_carbon13', mem%establishment%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_carbon14', mem%establishment%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_nitrogen15', mem%establishment%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_carbon', mem%reserve_use%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_nitrogen', mem%reserve_use%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_phosphorus', mem%reserve_use%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_carbon13', mem%reserve_use%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_carbon14', mem%reserve_use%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_nitrogen15', mem%reserve_use%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      !> add variables to the l2 structs (pool and flux variables, veg\%t_jsb_n_struct\%t_jsb_struct\%VAR)
      !! 
      !! slm_add_var -> slm_add_struct_element_r1d    *slmtemp 1D -> 2D
      !! 
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_carbon', mem%pool%leaf%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_nitrogen', mem%pool%leaf%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_phosphorus', mem%pool%leaf%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_carbon13', mem%pool%leaf%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_carbon14', mem%pool%leaf%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_nitrogen15', mem%pool%leaf%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_carbon', mem%pool%fine_root%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_nitrogen', mem%pool%fine_root%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_phosphorus', mem%pool%fine_root%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_carbon13', mem%pool%fine_root%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_carbon14', mem%pool%fine_root%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_nitrogen15', mem%pool%fine_root%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_carbon', mem%pool%coarse_root%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_nitrogen', mem%pool%coarse_root%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_phosphorus', mem%pool%coarse_root%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_carbon13', mem%pool%coarse_root%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_carbon14', mem%pool%coarse_root%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_nitrogen15', mem%pool%coarse_root%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_carbon', mem%pool%sap_wood%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_nitrogen', mem%pool%sap_wood%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_phosphorus', mem%pool%sap_wood%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_carbon13', mem%pool%sap_wood%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_carbon14', mem%pool%sap_wood%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_nitrogen15', mem%pool%sap_wood%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_carbon', mem%pool%heart_wood%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_nitrogen', mem%pool%heart_wood%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_phosphorus', mem%pool%heart_wood%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_carbon13', mem%pool%heart_wood%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_carbon14', mem%pool%heart_wood%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_nitrogen15', mem%pool%heart_wood%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_carbon', mem%pool%labile%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_nitrogen', mem%pool%labile%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_phosphorus', mem%pool%labile%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_carbon13', mem%pool%labile%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_carbon14', mem%pool%labile%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_nitrogen15', mem%pool%labile%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_carbon', mem%pool%reserve%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_nitrogen', mem%pool%reserve%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_phosphorus', mem%pool%reserve%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_carbon13', mem%pool%reserve%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_carbon14', mem%pool%reserve%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_nitrogen15', mem%pool%reserve%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_carbon', mem%pool%fruit%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_nitrogen', mem%pool%fruit%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_phosphorus', mem%pool%fruit%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_carbon13', mem%pool%fruit%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_carbon14', mem%pool%fruit%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_nitrogen15', mem%pool%fruit%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_carbon', mem%pool%seed_bed%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_nitrogen', mem%pool%seed_bed%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_phosphorus', mem%pool%seed_bed%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_carbon13', mem%pool%seed_bed%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_carbon14', mem%pool%seed_bed%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_nitrogen15', mem%pool%seed_bed%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_carbon', mem%growth%leaf%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_nitrogen', mem%growth%leaf%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_phosphorus', mem%growth%leaf%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_carbon13', mem%growth%leaf%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_carbon14', mem%growth%leaf%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_nitrogen15', mem%growth%leaf%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_carbon', mem%growth%fine_root%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_nitrogen', mem%growth%fine_root%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_phosphorus', mem%growth%fine_root%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_carbon13', mem%growth%fine_root%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_carbon14', mem%growth%fine_root%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_nitrogen15', mem%growth%fine_root%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_carbon', mem%growth%coarse_root%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_nitrogen', mem%growth%coarse_root%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_phosphorus', mem%growth%coarse_root%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_carbon13', mem%growth%coarse_root%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_carbon14', mem%growth%coarse_root%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_nitrogen15', mem%growth%coarse_root%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_carbon', mem%growth%sap_wood%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_nitrogen', mem%growth%sap_wood%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_phosphorus', mem%growth%sap_wood%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_carbon13', mem%growth%sap_wood%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_carbon14', mem%growth%sap_wood%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_nitrogen15', mem%growth%sap_wood%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_carbon', mem%growth%heart_wood%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_nitrogen', mem%growth%heart_wood%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_phosphorus', mem%growth%heart_wood%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_carbon13', mem%growth%heart_wood%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_carbon14', mem%growth%heart_wood%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_nitrogen15', mem%growth%heart_wood%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_carbon', mem%growth%labile%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_nitrogen', mem%growth%labile%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_phosphorus', mem%growth%labile%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_carbon13', mem%growth%labile%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_carbon14', mem%growth%labile%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_nitrogen15', mem%growth%labile%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_carbon', mem%growth%reserve%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_nitrogen', mem%growth%reserve%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_phosphorus', mem%growth%reserve%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_carbon13', mem%growth%reserve%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_carbon14', mem%growth%reserve%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_nitrogen15', mem%growth%reserve%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_carbon', mem%growth%fruit%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_nitrogen', mem%growth%fruit%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_phosphorus', mem%growth%fruit%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_carbon13', mem%growth%fruit%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_carbon14', mem%growth%fruit%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_nitrogen15', mem%growth%fruit%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_carbon', mem%growth%seed_bed%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_nitrogen', mem%growth%seed_bed%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_phosphorus', mem%growth%seed_bed%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_carbon13', mem%growth%seed_bed%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_carbon14', mem%growth%seed_bed%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_nitrogen15', mem%growth%seed_bed%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_carbon', mem%litterfall%leaf%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_nitrogen', mem%litterfall%leaf%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_phosphorus', mem%litterfall%leaf%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_carbon13', mem%litterfall%leaf%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_carbon14', mem%litterfall%leaf%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_nitrogen15', mem%litterfall%leaf%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_carbon', mem%litterfall%fine_root%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_nitrogen', mem%litterfall%fine_root%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_phosphorus', mem%litterfall%fine_root%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_carbon13', mem%litterfall%fine_root%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_carbon14', mem%litterfall%fine_root%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_nitrogen15', mem%litterfall%fine_root%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_carbon', mem%litterfall%coarse_root%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_nitrogen', mem%litterfall%coarse_root%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_phosphorus', mem%litterfall%coarse_root%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_carbon13', mem%litterfall%coarse_root%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_carbon14', mem%litterfall%coarse_root%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_nitrogen15', mem%litterfall%coarse_root%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_carbon', mem%litterfall%sap_wood%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_nitrogen', mem%litterfall%sap_wood%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_phosphorus', mem%litterfall%sap_wood%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_carbon13', mem%litterfall%sap_wood%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_carbon14', mem%litterfall%sap_wood%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_nitrogen15', mem%litterfall%sap_wood%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_carbon', mem%litterfall%heart_wood%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_nitrogen', mem%litterfall%heart_wood%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_phosphorus', mem%litterfall%heart_wood%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_carbon13', mem%litterfall%heart_wood%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_carbon14', mem%litterfall%heart_wood%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_nitrogen15', mem%litterfall%heart_wood%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_carbon', mem%litterfall%labile%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_nitrogen', mem%litterfall%labile%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_phosphorus', mem%litterfall%labile%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_carbon13', mem%litterfall%labile%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_carbon14', mem%litterfall%labile%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_nitrogen15', mem%litterfall%labile%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_carbon', mem%litterfall%reserve%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_nitrogen', mem%litterfall%reserve%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_phosphorus', mem%litterfall%reserve%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_carbon13', mem%litterfall%reserve%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_carbon14', mem%litterfall%reserve%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_nitrogen15', mem%litterfall%reserve%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_carbon', mem%litterfall%fruit%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_nitrogen', mem%litterfall%fruit%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_phosphorus', mem%litterfall%fruit%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_carbon13', mem%litterfall%fruit%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_carbon14', mem%litterfall%fruit%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_nitrogen15', mem%litterfall%fruit%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_carbon', mem%litterfall%seed_bed%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_nitrogen', mem%litterfall%seed_bed%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_phosphorus', mem%litterfall%seed_bed%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_carbon13', mem%litterfall%seed_bed%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_carbon14', mem%litterfall%seed_bed%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_nitrogen15', mem%litterfall%seed_bed%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue

    END IF

  END SUBROUTINE Init_veg_memory

END MODULE mo_veg_memory_class
