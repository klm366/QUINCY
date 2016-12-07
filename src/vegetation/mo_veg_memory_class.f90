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
                                       t_jsb_var_logical2d                        , &                 ! *slmjsb4
                                       t_jsb_struct, t_jsb_n_struct               , &                 ! *slmjsb4
                                       t_jsb_var_real1d, t_jsb_var_logical1d                          ! *slmtemp

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later    *slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_veg_memory, max_no_of_vars
  PUBLIC :: t_slm_struct_veg_pool, t_slm_struct_veg_element   ! *slmjsb4 
  PUBLIC :: max_no_of_structs, max_no_of_n_structs            ! *slmjsb4

  INTEGER, PARAMETER :: max_no_of_vars      = 300  !> TYPE t_jsb_var
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
    PROCEDURE, PUBLIC       :: init => init_t_slm_struct_veg_element
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
    PROCEDURE, PUBLIC       :: init => init_t_slm_struct_veg_pool
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
    TYPE(t_jsb_var_logical1d)       :: growing_season     !< logical to identify whether the plant is in the growing season or not
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
                                    beta_soil_ps      , &            !< scaling factor to account for soil moisture constraints on photosynthesis
                                    beta_sinklim_ps   , &            !< scaling factor to account for a direct sink limitation contraint on photosynthesis
                                    dphi              , &            !< plant-soil water potential gradient (MPa)
                                    beta_soil_gs      , &            !< scaling factor to account for soil moisture constraints on conductance
                                    growth_req_n      , &            !< moles N required for a unit of C growth under current allocation fractions 
                                    growth_req_p      , &            !< moles P required for a unit of N growth under current allocation fractions 
                                    gdd               , &            !< growing degree days above 5°C
                                    nd_dormance       , &            !< number of days since last growing season 
                                    t_jmax_opt        , &            !< current optimum temperature for eletron transport of photosynthesis
                                    an_boc            , &            !< net photosynthesis at bottom of canopy (lowest canopy layer existing)
                                    leaf_cn_direction , &            !< indicates direction of change in leaf CN ratio to maximise NPP
                                    apar              , &            !< ..
                                    target_lai        , &            !< LAI implied from sap wood area
                                    kstar_labile                     !< labile pool turnover rate assuming only temperature and moisture constraints (no stoichiometry)

    ! target C:N and N:P ratios for various non-leaf tissues (mol/mol)
    TYPE(t_jsb_var_real1d)          :: &
                                    target_cn_fine_root       , &
                                    target_cn_coarse_root     , &
                                    target_cn_sap_wood        , &
                                    target_cn_fruit           , &
                                    target_np_fine_root       , &
                                    target_np_coarse_root     , &
                                    target_np_sap_wood        , &
                                    target_np_fruit

    ! canopy layers -> levels, lbounds, ubounds
    TYPE(t_jsb_var_real2d)          :: &   ! DIMENSION(ncanopy)
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
    TYPE(t_jsb_var_real2d)          :: &   ! DIMENSION(ncanopy)
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
                                    beta_sinklim_ps_daytime   , &    !<  average of the previous day
                                    beta_soil_gs_daytime      , &    !<  average of the previous day
                                    t_jmax_opt_daytime               !<  average of the previous day
    TYPE(t_jsb_var_real2d)          :: &   ! DIMENSION(ncanopy)
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
                                    beta_sinklim_ps_daytime_dacc  , & !<  average of the previous day
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
    TYPE(t_jsb_var_real2d)          :: &   ! DIMENSION(ncanopy)
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
                                    beta_sinklim_ps_tfrac_mavg     , &  !< scaling factor to account for soil moisture constraints on photosynthesis 
                                    beta_soil_gs_tfrac_mavg             !< scaling factor to account for soil moisture constraints on conductance
    ! 1.5 memory of conditions for photosynthesis for the optimisation of foliar N  
    TYPE(t_jsb_var_real2d)          :: &   ! DIMENSION(ncanopy)
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
                                    beta_sinklim_ps_tcnl_mavg      , &  !< scaling factor to account for soil moisture constraints on conductance
                                    beta_soil_gs_tcnl_mavg              !< scaling factor to account for soil moisture constraints on conductance
    ! 1.6 memory of conditions for photosynthesis for the optimisation of biomass allocation fractions
    TYPE(t_jsb_var_real1d)          :: &   
                                    npp_talloc_mavg                , &  !< NPP (µmol m-2 s-1) (averaging time-scale: tau_allocation)
                                    transpiration_talloc_mavg      , &  !< average transpiration at allocation time-scale
                                    unit_npp_talloc_mavg           , &  !< rate of change in NPP with  increment in leaf mass, mumol C mol-1
                                    unit_transpiration_talloc_mavg , &  !< rate of tranpiration change per unit LAI change?
                                    unit_uptake_n_act_talloc_mavg  , &  !< actual rate of N uptkae per unit root, mumol N mol-1 (averaging time-scale: tau_allocation)
                                    unit_uptake_p_act_talloc_mavg  , &  !< actual rate of P uptake per unit root, mumol P mol-1 (time-scale of accumulation: tau_allocation)
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
                                    unit_uptake_n_pot_troot_mavg   , &  !< rate of N uptkae per unit root, mumol N mol-1 (time-scale of accumulation: tau_root)
                                    unit_uptake_p_pot_troot_mavg   , &  !< rate of P uptake per unit root, mumol P mol-1 (time-scale of accumulation: tau_root)
                                    fmaint_rate_troot_mavg              !< moving average of N-specific maintentance respiration (~µmol C/mol N / s ) 
    ! 1.8 memory over timescale of vegetation dynamics 
    TYPE(t_jsb_var_real1d)          :: &   
                                    an_boc_tvegdyn_mavg                 !< net photosynthesis at bottom of canopy (lowest canopy layer existing)
    ! 1.9 memory of variables for respiration accclimation 
    TYPE(t_jsb_var_real1d)          :: &   
                                    t_air_tacclim_mavg                  !< average temperature for respiration acclimation
    ! 1.10 memory of certain variables, temperature, ...
    TYPE(t_jsb_var_real1d)          :: &   
                                    t_air_week_mavg                , &  !< avrg air temperature across 7 days
                                    t_air_month_mavg               , &  !< avrg air temperature across 30 days
                                    t_jmax_opt_mavg                     !< optimim temperature for electron transport of photosynthesis

    !> structures
    !! 
    TYPE(t_slm_struct_veg_element)  :: &
                                    establishment   , &   !< flux between seed bed and reserve (mol m-2 timestep-1)
                                    reserve_use     , &   !< net flux between reserve and labile pool (mol m-2 timestep-1)
                                    init_pool       , &   !< initialization of some plant variables, considering: a) PFT, and b) flag_tile_veg_dynamics
                                                          !! used only in mo_jsb_model_init:init_slm_plant()
                                    total_biomass         !< per time step sum over veg\%pool\%..

    TYPE(t_slm_struct_veg_pool)     :: &   ! nested pool TYPE t_jsb_n_struct
                                    pool            , &   !< pools of the plant
                                    growth          , &   !< fluxes: growth
                                    litterfall      , &   !< fluxes: litterfall
                                    growth_t        , &   !< fluxes: growth     --> for veg_turnover      *slmdev
                                    litterfall_t    , &   !< fluxes: litterfall --> for veg_turnover      *slmdev
                                    falloc          , &   !< mo_veg_growth:calc_allocation_fraction -> fractional allocation to pools
                                    change_pool     , &   !< for diagnoses of mass conservation: changes in pools per timestep
                                    total_flux            !< for diagnoses of mass conservation: sum of fluxes


  CONTAINS
    PROCEDURE :: Init => Init_veg_memory
  END TYPE t_veg_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class'


CONTAINS

  !-----------------------------------------------------------------------------------------------------!
  ! t_slm_struct_veg_element operator overload                                                          !
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

    ALLOCATE(struct3, source=this)
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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this2)

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

!    ALLOCATE(struct3, source=this)

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

!    ALLOCATE(struct3, source=this2)

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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this2)

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

!    ALLOCATE(struct3, source=this)

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

!    ALLOCATE(struct3, source=this2)

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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this2)

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

!    ALLOCATE(struct3, source=this)

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

!    ALLOCATE(struct3, source=this2)

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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this2)

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

!    ALLOCATE(struct3, source=this)

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

!    ALLOCATE(struct3, source=this2)

!    struct3%carbon       = field2d / this2%carbon
!    struct3%nitrogen     = field2d / this2%nitrogen
!    struct3%phosphorus   = field2d / this2%phosphorus
!    struct3%carbon13     = field2d / this2%carbon13
!    struct3%carbon14     = field2d / this2%carbon14
!    struct3%nitrogen15   = field2d / this2%nitrogen15
!  END FUNCTION Divide_field2d_and_struct

  !-----------------------------------------------------------------------------------------------------
  !> Initialisation t_slm_struct_veg_element with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE init_t_slm_struct_veg_element(this,s)
    IMPLICIT NONE
    CLASS(t_slm_struct_veg_element),INTENT(inout) :: this 
    REAL(wp),                       INTENT(in)    :: s

    CALL this%carbon%init(s) 
    CALL this%nitrogen%init(s) 
    CALL this%phosphorus%init(s) 
    CALL this%carbon13%init(s) 
    CALL this%carbon14%init(s) 
    CALL this%nitrogen15%init(s) 
  END SUBROUTINE init_t_slm_struct_veg_element 


  !-----------------------------------------------------------------------------------------------------!
  ! t_slm_struct_veg_pool operator overload                                                             !
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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this2)

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

!    ALLOCATE(struct3, source=this)

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

!    ALLOCATE(struct3, source=this2)

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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this2)

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

!    ALLOCATE(struct3, source=this)

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

!    ALLOCATE(struct3, source=this2)

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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this2)

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

!    ALLOCATE(struct3, source=this)

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

!    ALLOCATE(struct3, source=this2)

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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this)

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

    ALLOCATE(struct3, source=this2)

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

!    ALLOCATE(struct3, source=this)

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

!    ALLOCATE(struct3, source=this2)

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
  !> Initialisation t_slm_struct_veg_pool with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE init_t_slm_struct_veg_pool(this,s)
    IMPLICIT NONE
    CLASS(t_slm_struct_veg_pool), INTENT(inout) :: this 
    REAL(wp),                     INTENT(in)    :: s

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
!    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
!    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_atmland_constants,   ONLY: standard_press_srf, def_co2_mixing_ratio
    USE mo_base_memory,         ONLY: t_model

    ! 0.1.InOut
    CLASS(t_veg_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids !(:)
    INTEGER,             INTENT(in)    :: model_id
    ! 0.2 Local
    TYPE(t_model), POINTER :: model
!    TYPE(t_jsb_model), POINTER :: model
!    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
!    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer        ! Vertical grid
!    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_veg_memory'

    model => Get_model(model_id)

!    table = tables(1)

!    hgrid   => Get_grid(mem%grid_id)
!    surface => Get_vgrid('surface')
!    canopy_layer => Get_vgrid('canopy_layer')

!    IF (.NOT. ASSOCIATED(mem%growing_season%ptr)) THEN
    IF (.NOT. ALLOCATED(mem%growing_season%ptr)) THEN

      !> add variables to veg memory veg\%VAR
      !! 
      !! slm_add_var -> slm_add_var_list_element_r1d    *slmtemp 1D -> 2D
      !!                slm_add_var_list_element_l1d
      !!  
      CALL slm_add_var(mem, 'growing_season', mem%growing_season, &
            initval_r=.TRUE.) ! initvalue
  
      CALL slm_add_var(mem, 'height', mem%height, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'diameter', mem%diameter, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'dens_ind', mem%dens_ind, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'lai', mem%lai, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'target_cn_leaf', mem%target_cn_leaf, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'target_np_leaf', mem%target_np_leaf, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'mean_leaf_age', mem%mean_leaf_age, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gc', mem%gc, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'f_n_demand', mem%f_n_demand, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'f_p_demand', mem%f_p_demand, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air', mem%beta_air, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps', mem%beta_soil_ps, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_sinklim_ps', mem%beta_sinklim_ps, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'dphi', mem%dphi, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs', mem%beta_soil_gs, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_req_n', mem%growth_req_n, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_req_p', mem%growth_req_p, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gdd', mem%gdd, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'nd_dormance', mem%nd_dormance, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_jmax_opt', mem%t_jmax_opt, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'an_boc', mem%an_boc, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'leaf_cn_direction', mem%leaf_cn_direction, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'apar', mem%apar, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'target_lai', mem%target_lai, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'kstar_labile', mem%kstar_labile, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'target_cn_fine_root', mem%target_cn_fine_root, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'target_cn_coarse_root', mem%target_cn_coarse_root, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'target_cn_sap_wood', mem%target_cn_sap_wood, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'target_cn_fruit', mem%target_cn_fruit, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'target_np_fine_root', mem%target_np_fine_root, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'target_np_coarse_root', mem%target_np_coarse_root, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'target_np_sap_wood', mem%target_np_sap_wood, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'target_np_fruit', mem%target_np_fruit, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'nl_tot_cl', mem%nl_tot_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fnl_chl_cl', mem%fnl_chl_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fnl_et_cl', mem%fnl_et_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fnl_rub_cl', mem%fnl_rub_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fnl_pepc_cl', mem%fnl_pepc_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fnl_oth_cl', mem%fnl_oth_cl, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_sunlit_cl', mem%par_sunlit_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_shaded_cl', mem%par_shaded_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fleaf_sunlit_cl', mem%fleaf_sunlit_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'lai_cl', mem%lai_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'cumm_lai_cl', mem%cumm_lai_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gpp', mem%gpp, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'an', mem%an, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gpp_c13', mem%gpp_c13, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gpp_c14', mem%gpp_c14, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'maint_respiration_leaf', mem%maint_respiration_leaf, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'maint_respiration_pot', mem%maint_respiration_pot, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'maint_respiration', mem%maint_respiration, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'maint_respiration_c13', mem%maint_respiration_c13, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'maint_respiration_c14', mem%maint_respiration_c14, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_respiration', mem%growth_respiration, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_respiration_c13', mem%growth_respiration_c13, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_respiration_c14', mem%growth_respiration_c14, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n_transform_respiration', mem%n_transform_respiration, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n_fixation_respiration', mem%n_fixation_respiration, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n_processing_respiration', mem%n_processing_respiration, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n_processing_respiration_c13', mem%n_processing_respiration_c13, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n_processing_respiration_c14', mem%n_processing_respiration_c14, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'npp', mem%npp, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'npp_c13', mem%npp_c13, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'npp_c14', mem%npp_c14, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'uptake_n', mem%uptake_n, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'uptake_n15', mem%uptake_n15, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'uptake_p', mem%uptake_p, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n_fixation', mem%n_fixation, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n_fixation_n15', mem%n_fixation_n15, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_leaf_n', mem%recycling_leaf_n, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_leaf_n15', mem%recycling_leaf_n15, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_leaf_p', mem%recycling_leaf_p, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_fine_root_n', mem%recycling_fine_root_n, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_fine_root_n15', mem%recycling_fine_root_n15, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_fine_root_p', mem%recycling_fine_root_p, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_sap_wood_n', mem%recycling_sap_wood_n, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_sap_wood_n15', mem%recycling_sap_wood_n15, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'recycling_sap_wood_p', mem%recycling_sap_wood_p, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transpiration', mem%transpiration, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fmaint_rate', mem%fmaint_rate, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'ci', mem%ci, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_npp', mem%unit_npp, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_transpiration', mem%unit_transpiration, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_uptake_n_pot', mem%unit_uptake_n_pot, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_uptake_p_pot', mem%unit_uptake_p_pot, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_uptake_n_act', mem%unit_uptake_n_act, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_uptake_p_act', mem%unit_uptake_p_act, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'delta_dens_ind', mem%delta_dens_ind, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'cost_n_uptake_root', mem%cost_n_uptake_root, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'cost_n_fixation', mem%cost_n_fixation, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'mortality_rate', mem%mortality_rate, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'k1_opt', mem%k1_opt, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'k2_opt', mem%k2_opt, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'leaf2sapwood_mass_ratio', mem%leaf2sapwood_mass_ratio, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'leaf2root_mass_ratio', mem%leaf2root_mass_ratio, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_sunlit_daytime_cl', mem%par_sunlit_daytime_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_shaded_daytime_cl', mem%par_shaded_daytime_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fleaf_sunlit_daytime_cl', mem%fleaf_sunlit_daytime_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_daytime', mem%t_air_daytime, &
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 'press_srf_daytime', mem%press_srf_daytime, &
            initval_r=standard_press_srf) ! initvalue
  
      CALL slm_add_var(mem, 'co2_mixing_ratio_daytime', mem%co2_mixing_ratio_daytime, &
            initval_r=def_co2_mixing_ratio) ! initvalue
  
      CALL slm_add_var(mem, 'ga_daytime', mem%ga_daytime, &
            initval_r=0.2_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air_daytime', mem%beta_air_daytime, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps_daytime', mem%beta_soil_ps_daytime, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_sinklim_ps_daytime', mem%beta_sinklim_ps_daytime, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_daytime', mem%beta_soil_gs_daytime, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_jmax_opt_daytime', mem%t_jmax_opt_daytime, &
            initval_r=17.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_sunlit_daytime_dacc_cl', mem%par_sunlit_daytime_dacc_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_shaded_daytime_dacc_cl', mem%par_shaded_daytime_dacc_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fleaf_sunlit_daytime_dacc_cl', mem%fleaf_sunlit_daytime_dacc_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_daytime_dacc', mem%t_air_daytime_dacc, &
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 'press_srf_daytime_dacc', mem%press_srf_daytime_dacc, &
            initval_r=standard_press_srf) ! initvalue
  
      CALL slm_add_var(mem, 'co2_mixing_ratio_daytime_dacc', mem%co2_mixing_ratio_daytime_dacc, &
            initval_r=def_co2_mixing_ratio) ! initvalue
  
      CALL slm_add_var(mem, 'ga_daytime_dacc', mem%ga_daytime_dacc, &
            initval_r=0.2_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air_daytime_dacc', mem%beta_air_daytime_dacc, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps_daytime_dacc', mem%beta_soil_ps_daytime_dacc, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_sinklim_ps_daytime_dacc', mem%beta_sinklim_ps_daytime_dacc, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_daytime_dacc', mem%beta_soil_gs_daytime_dacc, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_jmax_opt_daytime_dacc', mem%t_jmax_opt_daytime_dacc, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'daytime_counter', mem%daytime_counter, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gpp_tlabile_mavg', mem%gpp_tlabile_mavg, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'maint_respiration_tlabile_mavg', mem%maint_respiration_tlabile_mavg, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_req_n_tlabile_mavg', mem%growth_req_n_tlabile_mavg, &
            initval_r=1._wp/25._wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_req_p_tlabile_mavg', mem%growth_req_p_tlabile_mavg, &
            initval_r=1._wp/14._wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_tphen_mavg', mem%t_air_tphen_mavg, &
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_tphen_mavg', mem%beta_soil_gs_tphen_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'npp_tuptake_mavg', mem%npp_tuptake_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'demand_uptake_n_tuptake_mavg', mem%demand_uptake_n_tuptake_mavg, &
            initval_r=0.1_wp) ! initvalue
  
      CALL slm_add_var(mem, 'demand_uptake_p_tuptake_mavg', mem%demand_uptake_p_tuptake_mavg, &
            initval_r=0.1_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_req_n_tuptake_mavg', mem%growth_req_n_tuptake_mavg, &
            initval_r=1._wp/25._wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_req_p_tuptake_mavg', mem%growth_req_p_tuptake_mavg, &
            initval_r=1._wp/14._wp) ! initvalue
  
      CALL slm_add_var(mem, 'cost_n_fixation_tuptake_mavg', mem%cost_n_fixation_tuptake_mavg, &
            initval_r=10.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_sunlit_tfrac_mavg_cl', mem%par_sunlit_tfrac_mavg_cl, &
            initval_r=150.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_shaded_tfrac_mavg_cl', mem%par_shaded_tfrac_mavg_cl, &
            initval_r=10.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fleaf_sunlit_tfrac_mavg_cl', mem%fleaf_sunlit_tfrac_mavg_cl, &
            initval_r=0.5_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_tfrac_mavg', mem%t_air_tfrac_mavg, &
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 'press_srf_tfrac_mavg', mem%press_srf_tfrac_mavg, &
            initval_r=standard_press_srf) ! initvalue
  
      CALL slm_add_var(mem, 'co2_mixing_ratio_tfrac_mavg', mem%co2_mixing_ratio_tfrac_mavg, &
            initval_r=def_co2_mixing_ratio) ! initvalue
  
      CALL slm_add_var(mem, 'ga_tfrac_mavg', mem%ga_tfrac_mavg, &
            initval_r=0.2_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air_tfrac_mavg', mem%beta_air_tfrac_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps_tfrac_mavg', mem%beta_soil_ps_tfrac_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_sinklim_ps_tfrac_mavg', mem%beta_sinklim_ps_tfrac_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_tfrac_mavg', mem%beta_soil_gs_tfrac_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_sunlit_tcnl_mavg_cl', mem%par_sunlit_tcnl_mavg_cl, &
            initval_r=150.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'par_shaded_tcnl_mavg_cl', mem%par_shaded_tcnl_mavg_cl, &
            initval_r=10.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fleaf_sunlit_tcnl_mavg_cl', mem%fleaf_sunlit_tcnl_mavg_cl, &
            initval_r=0.5_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_tcnl_mavg', mem%t_air_tcnl_mavg, &
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 'press_srf_tcnl_mavg', mem%press_srf_tcnl_mavg, &
            initval_r=standard_press_srf) ! initvalue
  
      CALL slm_add_var(mem, 'co2_mixing_ratio_tcnl_mavg', mem%co2_mixing_ratio_tcnl_mavg, &
            initval_r=def_co2_mixing_ratio) ! initvalue
  
      CALL slm_add_var(mem, 'ga_tcnl_mavg', mem%ga_tcnl_mavg, &
            initval_r=0.2_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air_tcnl_mavg', mem%beta_air_tcnl_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps_tcnl_mavg', mem%beta_soil_ps_tcnl_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_sinklim_ps_tcnl_mavg', mem%beta_sinklim_ps_tcnl_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_tcnl_mavg', mem%beta_soil_gs_tcnl_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'npp_talloc_mavg', mem%npp_talloc_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transpiration_talloc_mavg', mem%transpiration_talloc_mavg, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_npp_talloc_mavg', mem%unit_npp_talloc_mavg, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_transpiration_talloc_mavg', mem%unit_transpiration_talloc_mavg, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_uptake_n_act_talloc_mavg', mem%unit_uptake_n_act_talloc_mavg, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_uptake_p_act_talloc_mavg', mem%unit_uptake_p_act_talloc_mavg, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'dphi_talloc_mavg', mem%dphi_talloc_mavg, &
            initval_r=-2.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_cn_talloc_mavg', mem%growth_cn_talloc_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'growth_cp_talloc_mavg', mem%growth_cp_talloc_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'labile_carbon_talloc_mavg', mem%labile_carbon_talloc_mavg, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'labile_nitrogen_talloc_mavg', mem%labile_nitrogen_talloc_mavg, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'labile_phosphorus_talloc_mavg', mem%labile_phosphorus_talloc_mavg, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'leaf2root_troot_mavg', mem%leaf2root_troot_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_npp_troot_mavg', mem%unit_npp_troot_mavg, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_uptake_n_pot_troot_mavg', mem%unit_uptake_n_pot_troot_mavg, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'unit_uptake_p_pot_troot_mavg', mem%unit_uptake_p_pot_troot_mavg, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'fmaint_rate_troot_mavg', mem%fmaint_rate_troot_mavg, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'an_boc_tvegdyn_mavg', mem%an_boc_tvegdyn_mavg, &
            initval_r=5.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_tacclim_mavg', mem%t_air_tacclim_mavg, &
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_week_mavg', mem%t_air_week_mavg, &
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_air_month_mavg', mem%t_air_month_mavg, &
            initval_r=283.15_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_jmax_opt_mavg', mem%t_jmax_opt_mavg, &
            initval_r=17.0_wp) ! initvalue

      !> simple TYPE (t_jsb_struct)
      !! 
      !! slm_add_var -> slm_add_struct()
      !! 
      CALL slm_add_var(mem, 'veg_establishment',          mem%establishment)
      CALL slm_add_var(mem, 'veg_reserve_use',            mem%reserve_use)
      CALL slm_add_var(mem, 'veg_init_pool',              mem%init_pool)
      CALL slm_add_var(mem, 'veg_total_biomass',          mem%total_biomass)
      !> nested TYPE (t_jsb_n_struct)
      !! 
      !! slm_add_var -> slm_add_n_struct()
      !! 
      CALL slm_add_var(mem, 'veg_pool'        , mem%pool)
      CALL slm_add_var(mem, 'veg_growth'      , mem%growth)
      CALL slm_add_var(mem, 'veg_litterfall'  , mem%litterfall)
      CALL slm_add_var(mem, 'veg_growth_t'    , mem%growth_t)       ! *slmdev
      CALL slm_add_var(mem, 'veg_litterfall_t', mem%litterfall_t)   ! *slmdev
      CALL slm_add_var(mem, 'veg_falloc'      , mem%falloc)
      CALL slm_add_var(mem, 'veg_change_pool' , mem%change_pool)    ! mass conservation diagnosis
      CALL slm_add_var(mem, 'veg_total_flux'  , mem%total_flux)     ! mass conservation diagnosis
      !> level 2 TYPE (l2_struct) nested within the above nested TYPE structures (t_jsb_n_struct)
      !! 
      !! slm_add_var -> slm_add_l2_struct()
      !! 
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
      CALL slm_add_var(mem%falloc, 'veg_falloc_leaf',         mem%falloc%leaf)       
      CALL slm_add_var(mem%falloc, 'veg_falloc_fine_root',    mem%falloc%fine_root)
      CALL slm_add_var(mem%falloc, 'veg_falloc_coarse_root',  mem%falloc%coarse_root)
      CALL slm_add_var(mem%falloc, 'veg_falloc_sap_wood',     mem%falloc%sap_wood)
      CALL slm_add_var(mem%falloc, 'veg_falloc_heart_wood',   mem%falloc%heart_wood)
      CALL slm_add_var(mem%falloc, 'veg_falloc_labile',       mem%falloc%labile)
      CALL slm_add_var(mem%falloc, 'veg_falloc_reserve',      mem%falloc%reserve)
      CALL slm_add_var(mem%falloc, 'veg_falloc_fruit',        mem%falloc%fruit)
      CALL slm_add_var(mem%falloc, 'veg_falloc_seed_bed',     mem%falloc%seed_bed)
      CALL slm_add_var(mem%growth_t, 'veg_growth_t_leaf',         mem%growth_t%leaf)       
      CALL slm_add_var(mem%growth_t, 'veg_growth_t_fine_root',    mem%growth_t%fine_root)
      CALL slm_add_var(mem%growth_t, 'veg_growth_t_coarse_root',  mem%growth_t%coarse_root)
      CALL slm_add_var(mem%growth_t, 'veg_growth_t_sap_wood',     mem%growth_t%sap_wood)
      CALL slm_add_var(mem%growth_t, 'veg_growth_t_heart_wood',   mem%growth_t%heart_wood)
      CALL slm_add_var(mem%growth_t, 'veg_growth_t_labile',       mem%growth_t%labile)
      CALL slm_add_var(mem%growth_t, 'veg_growth_t_reserve',      mem%growth_t%reserve)
      CALL slm_add_var(mem%growth_t, 'veg_growth_t_fruit',        mem%growth_t%fruit)
      CALL slm_add_var(mem%growth_t, 'veg_growth_t_seed_bed',     mem%growth_t%seed_bed)
      CALL slm_add_var(mem%litterfall_t, 'veg_litterfall_t_leaf',         mem%litterfall_t%leaf)       
      CALL slm_add_var(mem%litterfall_t, 'veg_litterfall_t_fine_root',    mem%litterfall_t%fine_root)
      CALL slm_add_var(mem%litterfall_t, 'veg_litterfall_t_coarse_root',  mem%litterfall_t%coarse_root)
      CALL slm_add_var(mem%litterfall_t, 'veg_litterfall_t_sap_wood',     mem%litterfall_t%sap_wood)
      CALL slm_add_var(mem%litterfall_t, 'veg_litterfall_t_heart_wood',   mem%litterfall_t%heart_wood)
      CALL slm_add_var(mem%litterfall_t, 'veg_litterfall_t_labile',       mem%litterfall_t%labile)
      CALL slm_add_var(mem%litterfall_t, 'veg_litterfall_t_reserve',      mem%litterfall_t%reserve)
      CALL slm_add_var(mem%litterfall_t, 'veg_litterfall_t_fruit',        mem%litterfall_t%fruit)
      CALL slm_add_var(mem%litterfall_t, 'veg_litterfall_t_seed_bed',     mem%litterfall_t%seed_bed)
      CALL slm_add_var(mem%change_pool, 'veg_change_pool_leaf',         mem%change_pool%leaf)       
      CALL slm_add_var(mem%change_pool, 'veg_change_pool_fine_root',    mem%change_pool%fine_root)
      CALL slm_add_var(mem%change_pool, 'veg_change_pool_coarse_root',  mem%change_pool%coarse_root)
      CALL slm_add_var(mem%change_pool, 'veg_change_pool_sap_wood',     mem%change_pool%sap_wood)
      CALL slm_add_var(mem%change_pool, 'veg_change_pool_heart_wood',   mem%change_pool%heart_wood)
      CALL slm_add_var(mem%change_pool, 'veg_change_pool_labile',       mem%change_pool%labile)
      CALL slm_add_var(mem%change_pool, 'veg_change_pool_reserve',      mem%change_pool%reserve)
      CALL slm_add_var(mem%change_pool, 'veg_change_pool_fruit',        mem%change_pool%fruit)
      CALL slm_add_var(mem%change_pool, 'veg_change_pool_seed_bed',     mem%change_pool%seed_bed)
      CALL slm_add_var(mem%total_flux, 'veg_total_flux_leaf',         mem%total_flux%leaf)       
      CALL slm_add_var(mem%total_flux, 'veg_total_flux_fine_root',    mem%total_flux%fine_root)
      CALL slm_add_var(mem%total_flux, 'veg_total_flux_coarse_root',  mem%total_flux%coarse_root)
      CALL slm_add_var(mem%total_flux, 'veg_total_flux_sap_wood',     mem%total_flux%sap_wood)
      CALL slm_add_var(mem%total_flux, 'veg_total_flux_heart_wood',   mem%total_flux%heart_wood)
      CALL slm_add_var(mem%total_flux, 'veg_total_flux_labile',       mem%total_flux%labile)
      CALL slm_add_var(mem%total_flux, 'veg_total_flux_reserve',      mem%total_flux%reserve)
      CALL slm_add_var(mem%total_flux, 'veg_total_flux_fruit',        mem%total_flux%fruit)
      CALL slm_add_var(mem%total_flux, 'veg_total_flux_seed_bed',     mem%total_flux%seed_bed)

      !> add variables to the structs and the l2 structs
      !!
      !! veg\%t_jsb_struct\%VAR
      !! 
      !! veg\%t_jsb_n_struct\%t_jsb_struct\%VAR
      !! 
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_carbon', mem%pool%leaf%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_nitrogen', mem%pool%leaf%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_phosphorus', mem%pool%leaf%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_carbon13', mem%pool%leaf%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_carbon14', mem%pool%leaf%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_nitrogen15', mem%pool%leaf%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_carbon', mem%pool%fine_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_nitrogen', mem%pool%fine_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_phosphorus', mem%pool%fine_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_carbon13', mem%pool%fine_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_carbon14', mem%pool%fine_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_nitrogen15', mem%pool%fine_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_carbon', mem%pool%coarse_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_nitrogen', mem%pool%coarse_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_phosphorus', mem%pool%coarse_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_carbon13', mem%pool%coarse_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_carbon14', mem%pool%coarse_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_nitrogen15', mem%pool%coarse_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_carbon', mem%pool%sap_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_nitrogen', mem%pool%sap_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_phosphorus', mem%pool%sap_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_carbon13', mem%pool%sap_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_carbon14', mem%pool%sap_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_nitrogen15', mem%pool%sap_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_carbon', mem%pool%heart_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_nitrogen', mem%pool%heart_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_phosphorus', mem%pool%heart_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_carbon13', mem%pool%heart_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_carbon14', mem%pool%heart_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_nitrogen15', mem%pool%heart_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_carbon', mem%pool%labile%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_nitrogen', mem%pool%labile%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_phosphorus', mem%pool%labile%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_carbon13', mem%pool%labile%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_carbon14', mem%pool%labile%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_nitrogen15', mem%pool%labile%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_carbon', mem%pool%reserve%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_nitrogen', mem%pool%reserve%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_phosphorus', mem%pool%reserve%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_carbon13', mem%pool%reserve%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_carbon14', mem%pool%reserve%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_nitrogen15', mem%pool%reserve%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_carbon', mem%pool%fruit%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_nitrogen', mem%pool%fruit%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_phosphorus', mem%pool%fruit%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_carbon13', mem%pool%fruit%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_carbon14', mem%pool%fruit%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_nitrogen15', mem%pool%fruit%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_carbon', mem%pool%seed_bed%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_nitrogen', mem%pool%seed_bed%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_phosphorus', mem%pool%seed_bed%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_carbon13', mem%pool%seed_bed%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_carbon14', mem%pool%seed_bed%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_nitrogen15', mem%pool%seed_bed%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_carbon', mem%growth%leaf%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_nitrogen', mem%growth%leaf%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_phosphorus', mem%growth%leaf%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_carbon13', mem%growth%leaf%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_carbon14', mem%growth%leaf%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_nitrogen15', mem%growth%leaf%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_carbon', mem%growth%fine_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_nitrogen', mem%growth%fine_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_phosphorus', mem%growth%fine_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_carbon13', mem%growth%fine_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_carbon14', mem%growth%fine_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_nitrogen15', mem%growth%fine_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_carbon', mem%growth%coarse_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_nitrogen', mem%growth%coarse_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_phosphorus', mem%growth%coarse_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_carbon13', mem%growth%coarse_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_carbon14', mem%growth%coarse_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_nitrogen15', mem%growth%coarse_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_carbon', mem%growth%sap_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_nitrogen', mem%growth%sap_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_phosphorus', mem%growth%sap_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_carbon13', mem%growth%sap_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_carbon14', mem%growth%sap_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_nitrogen15', mem%growth%sap_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_carbon', mem%growth%heart_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_nitrogen', mem%growth%heart_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_phosphorus', mem%growth%heart_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_carbon13', mem%growth%heart_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_carbon14', mem%growth%heart_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_nitrogen15', mem%growth%heart_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_carbon', mem%growth%labile%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_nitrogen', mem%growth%labile%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_phosphorus', mem%growth%labile%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_carbon13', mem%growth%labile%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_carbon14', mem%growth%labile%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_nitrogen15', mem%growth%labile%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_carbon', mem%growth%reserve%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_nitrogen', mem%growth%reserve%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_phosphorus', mem%growth%reserve%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_carbon13', mem%growth%reserve%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_carbon14', mem%growth%reserve%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_nitrogen15', mem%growth%reserve%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_carbon', mem%growth%fruit%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_nitrogen', mem%growth%fruit%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_phosphorus', mem%growth%fruit%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_carbon13', mem%growth%fruit%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_carbon14', mem%growth%fruit%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_nitrogen15', mem%growth%fruit%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_carbon', mem%growth%seed_bed%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_nitrogen', mem%growth%seed_bed%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_phosphorus', mem%growth%seed_bed%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_carbon13', mem%growth%seed_bed%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_carbon14', mem%growth%seed_bed%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_nitrogen15', mem%growth%seed_bed%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_carbon', mem%litterfall%leaf%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_nitrogen', mem%litterfall%leaf%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_phosphorus', mem%litterfall%leaf%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_carbon13', mem%litterfall%leaf%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_carbon14', mem%litterfall%leaf%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_nitrogen15', mem%litterfall%leaf%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_carbon', mem%litterfall%fine_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_nitrogen', mem%litterfall%fine_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_phosphorus', mem%litterfall%fine_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_carbon13', mem%litterfall%fine_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_carbon14', mem%litterfall%fine_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_nitrogen15', mem%litterfall%fine_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_carbon', mem%litterfall%coarse_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_nitrogen', mem%litterfall%coarse_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_phosphorus', mem%litterfall%coarse_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_carbon13', mem%litterfall%coarse_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_carbon14', mem%litterfall%coarse_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_nitrogen15', mem%litterfall%coarse_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_carbon', mem%litterfall%sap_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_nitrogen', mem%litterfall%sap_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_phosphorus', mem%litterfall%sap_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_carbon13', mem%litterfall%sap_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_carbon14', mem%litterfall%sap_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_nitrogen15', mem%litterfall%sap_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_carbon', mem%litterfall%heart_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_nitrogen', mem%litterfall%heart_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_phosphorus', mem%litterfall%heart_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_carbon13', mem%litterfall%heart_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_carbon14', mem%litterfall%heart_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_nitrogen15', mem%litterfall%heart_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_carbon', mem%litterfall%labile%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_nitrogen', mem%litterfall%labile%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_phosphorus', mem%litterfall%labile%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_carbon13', mem%litterfall%labile%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_carbon14', mem%litterfall%labile%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_nitrogen15', mem%litterfall%labile%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_carbon', mem%litterfall%reserve%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_nitrogen', mem%litterfall%reserve%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_phosphorus', mem%litterfall%reserve%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_carbon13', mem%litterfall%reserve%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_carbon14', mem%litterfall%reserve%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_nitrogen15', mem%litterfall%reserve%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_carbon', mem%litterfall%fruit%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_nitrogen', mem%litterfall%fruit%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_phosphorus', mem%litterfall%fruit%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_carbon13', mem%litterfall%fruit%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_carbon14', mem%litterfall%fruit%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_nitrogen15', mem%litterfall%fruit%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_carbon', mem%litterfall%seed_bed%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_nitrogen', mem%litterfall%seed_bed%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_phosphorus', mem%litterfall%seed_bed%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_carbon13', mem%litterfall%seed_bed%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_carbon14', mem%litterfall%seed_bed%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_nitrogen15', mem%litterfall%seed_bed%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%falloc%leaf, 'veg_falloc_leaf_carbon', mem%falloc%leaf%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%leaf, 'veg_falloc_leaf_nitrogen', mem%falloc%leaf%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%leaf, 'veg_falloc_leaf_phosphorus', mem%falloc%leaf%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%leaf, 'veg_falloc_leaf_carbon13', mem%falloc%leaf%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%leaf, 'veg_falloc_leaf_carbon14', mem%falloc%leaf%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%leaf, 'veg_falloc_leaf_nitrogen15', mem%falloc%leaf%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%falloc%fine_root, 'veg_falloc_fine_root_carbon', mem%falloc%fine_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%fine_root, 'veg_falloc_fine_root_nitrogen', mem%falloc%fine_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%fine_root, 'veg_falloc_fine_root_phosphorus', mem%falloc%fine_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%fine_root, 'veg_falloc_fine_root_carbon13', mem%falloc%fine_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%fine_root, 'veg_falloc_fine_root_carbon14', mem%falloc%fine_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%fine_root, 'veg_falloc_fine_root_nitrogen15', mem%falloc%fine_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%falloc%coarse_root, 'veg_falloc_coarse_root_carbon', mem%falloc%coarse_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%coarse_root, 'veg_falloc_coarse_root_nitrogen', mem%falloc%coarse_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%coarse_root, 'veg_falloc_coarse_root_phosphorus', mem%falloc%coarse_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%coarse_root, 'veg_falloc_coarse_root_carbon13', mem%falloc%coarse_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%coarse_root, 'veg_falloc_coarse_root_carbon14', mem%falloc%coarse_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%coarse_root, 'veg_falloc_coarse_root_nitrogen15', mem%falloc%coarse_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%falloc%sap_wood, 'veg_falloc_sap_wood_carbon', mem%falloc%sap_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%sap_wood, 'veg_falloc_sap_wood_nitrogen', mem%falloc%sap_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%sap_wood, 'veg_falloc_sap_wood_phosphorus', mem%falloc%sap_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%sap_wood, 'veg_falloc_sap_wood_carbon13', mem%falloc%sap_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%sap_wood, 'veg_falloc_sap_wood_carbon14', mem%falloc%sap_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%sap_wood, 'veg_falloc_sap_wood_nitrogen15', mem%falloc%sap_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%falloc%heart_wood, 'veg_falloc_heart_wood_carbon', mem%falloc%heart_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%heart_wood, 'veg_falloc_heart_wood_nitrogen', mem%falloc%heart_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%heart_wood, 'veg_falloc_heart_wood_phosphorus', mem%falloc%heart_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%heart_wood, 'veg_falloc_heart_wood_carbon13', mem%falloc%heart_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%heart_wood, 'veg_falloc_heart_wood_carbon14', mem%falloc%heart_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%heart_wood, 'veg_falloc_heart_wood_nitrogen15', mem%falloc%heart_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%falloc%labile, 'veg_falloc_labile_carbon', mem%falloc%labile%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%labile, 'veg_falloc_labile_nitrogen', mem%falloc%labile%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%labile, 'veg_falloc_labile_phosphorus', mem%falloc%labile%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%labile, 'veg_falloc_labile_carbon13', mem%falloc%labile%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%labile, 'veg_falloc_labile_carbon14', mem%falloc%labile%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%labile, 'veg_falloc_labile_nitrogen15', mem%falloc%labile%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%falloc%reserve, 'veg_falloc_reserve_carbon', mem%falloc%reserve%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%reserve, 'veg_falloc_reserve_nitrogen', mem%falloc%reserve%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%reserve, 'veg_falloc_reserve_phosphorus', mem%falloc%reserve%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%reserve, 'veg_falloc_reserve_carbon13', mem%falloc%reserve%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%reserve, 'veg_falloc_reserve_carbon14', mem%falloc%reserve%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%reserve, 'veg_falloc_reserve_nitrogen15', mem%falloc%reserve%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%falloc%fruit, 'veg_falloc_fruit_carbon', mem%falloc%fruit%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%fruit, 'veg_falloc_fruit_nitrogen', mem%falloc%fruit%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%fruit, 'veg_falloc_fruit_phosphorus', mem%falloc%fruit%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%fruit, 'veg_falloc_fruit_carbon13', mem%falloc%fruit%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%fruit, 'veg_falloc_fruit_carbon14', mem%falloc%fruit%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%fruit, 'veg_falloc_fruit_nitrogen15', mem%falloc%fruit%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%falloc%seed_bed, 'veg_falloc_seed_bed_carbon', mem%falloc%seed_bed%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%seed_bed, 'veg_falloc_seed_bed_nitrogen', mem%falloc%seed_bed%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%seed_bed, 'veg_falloc_seed_bed_phosphorus', mem%falloc%seed_bed%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%seed_bed, 'veg_falloc_seed_bed_carbon13', mem%falloc%seed_bed%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%seed_bed, 'veg_falloc_seed_bed_carbon14', mem%falloc%seed_bed%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%falloc%seed_bed, 'veg_falloc_seed_bed_nitrogen15', mem%falloc%seed_bed%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth_t%leaf, 'veg_growth_t_leaf_carbon', mem%growth_t%leaf%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%leaf, 'veg_growth_t_leaf_nitrogen', mem%growth_t%leaf%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%leaf, 'veg_growth_t_leaf_phosphorus', mem%growth_t%leaf%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%leaf, 'veg_growth_t_leaf_carbon13', mem%growth_t%leaf%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%leaf, 'veg_growth_t_leaf_carbon14', mem%growth_t%leaf%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%leaf, 'veg_growth_t_leaf_nitrogen15', mem%growth_t%leaf%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth_t%fine_root, 'veg_growth_t_fine_root_carbon', mem%growth_t%fine_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%fine_root, 'veg_growth_t_fine_root_nitrogen', mem%growth_t%fine_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%fine_root, 'veg_growth_t_fine_root_phosphorus', mem%growth_t%fine_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%fine_root, 'veg_growth_t_fine_root_carbon13', mem%growth_t%fine_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%fine_root, 'veg_growth_t_fine_root_carbon14', mem%growth_t%fine_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%fine_root, 'veg_growth_t_fine_root_nitrogen15', mem%growth_t%fine_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth_t%coarse_root, 'veg_growth_t_coarse_root_carbon', mem%growth_t%coarse_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%coarse_root, 'veg_growth_t_coarse_root_nitrogen', mem%growth_t%coarse_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%coarse_root, 'veg_growth_t_coarse_root_phosphorus', mem%growth_t%coarse_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%coarse_root, 'veg_growth_t_coarse_root_carbon13', mem%growth_t%coarse_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%coarse_root, 'veg_growth_t_coarse_root_carbon14', mem%growth_t%coarse_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%coarse_root, 'veg_growth_t_coarse_root_nitrogen15', mem%growth_t%coarse_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth_t%sap_wood, 'veg_growth_t_sap_wood_carbon', mem%growth_t%sap_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%sap_wood, 'veg_growth_t_sap_wood_nitrogen', mem%growth_t%sap_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%sap_wood, 'veg_growth_t_sap_wood_phosphorus', mem%growth_t%sap_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%sap_wood, 'veg_growth_t_sap_wood_carbon13', mem%growth_t%sap_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%sap_wood, 'veg_growth_t_sap_wood_carbon14', mem%growth_t%sap_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%sap_wood, 'veg_growth_t_sap_wood_nitrogen15', mem%growth_t%sap_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth_t%heart_wood, 'veg_growth_t_heart_wood_carbon', mem%growth_t%heart_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%heart_wood, 'veg_growth_t_heart_wood_nitrogen', mem%growth_t%heart_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%heart_wood, 'veg_growth_t_heart_wood_phosphorus', mem%growth_t%heart_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%heart_wood, 'veg_growth_t_heart_wood_carbon13', mem%growth_t%heart_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%heart_wood, 'veg_growth_t_heart_wood_carbon14', mem%growth_t%heart_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%heart_wood, 'veg_growth_t_heart_wood_nitrogen15', mem%growth_t%heart_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth_t%labile, 'veg_growth_t_labile_carbon', mem%growth_t%labile%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%labile, 'veg_growth_t_labile_nitrogen', mem%growth_t%labile%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%labile, 'veg_growth_t_labile_phosphorus', mem%growth_t%labile%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%labile, 'veg_growth_t_labile_carbon13', mem%growth_t%labile%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%labile, 'veg_growth_t_labile_carbon14', mem%growth_t%labile%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%labile, 'veg_growth_t_labile_nitrogen15', mem%growth_t%labile%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth_t%reserve, 'veg_growth_t_reserve_carbon', mem%growth_t%reserve%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%reserve, 'veg_growth_t_reserve_nitrogen', mem%growth_t%reserve%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%reserve, 'veg_growth_t_reserve_phosphorus', mem%growth_t%reserve%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%reserve, 'veg_growth_t_reserve_carbon13', mem%growth_t%reserve%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%reserve, 'veg_growth_t_reserve_carbon14', mem%growth_t%reserve%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%reserve, 'veg_growth_t_reserve_nitrogen15', mem%growth_t%reserve%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth_t%fruit, 'veg_growth_t_fruit_carbon', mem%growth_t%fruit%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%fruit, 'veg_growth_t_fruit_nitrogen', mem%growth_t%fruit%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%fruit, 'veg_growth_t_fruit_phosphorus', mem%growth_t%fruit%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%fruit, 'veg_growth_t_fruit_carbon13', mem%growth_t%fruit%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%fruit, 'veg_growth_t_fruit_carbon14', mem%growth_t%fruit%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%fruit, 'veg_growth_t_fruit_nitrogen15', mem%growth_t%fruit%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth_t%seed_bed, 'veg_growth_t_seed_bed_carbon', mem%growth_t%seed_bed%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%seed_bed, 'veg_growth_t_seed_bed_nitrogen', mem%growth_t%seed_bed%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%seed_bed, 'veg_growth_t_seed_bed_phosphorus', mem%growth_t%seed_bed%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%seed_bed, 'veg_growth_t_seed_bed_carbon13', mem%growth_t%seed_bed%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%seed_bed, 'veg_growth_t_seed_bed_carbon14', mem%growth_t%seed_bed%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth_t%seed_bed, 'veg_growth_t_seed_bed_nitrogen15', mem%growth_t%seed_bed%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall_t%leaf, 'veg_litterfall_t_leaf_carbon', mem%litterfall_t%leaf%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%leaf, 'veg_litterfall_t_leaf_nitrogen', mem%litterfall_t%leaf%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%leaf, 'veg_litterfall_t_leaf_phosphorus', mem%litterfall_t%leaf%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%leaf, 'veg_litterfall_t_leaf_carbon13', mem%litterfall_t%leaf%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%leaf, 'veg_litterfall_t_leaf_carbon14', mem%litterfall_t%leaf%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%leaf, 'veg_litterfall_t_leaf_nitrogen15', mem%litterfall_t%leaf%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall_t%fine_root, 'veg_litterfall_t_fine_root_carbon', mem%litterfall_t%fine_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%fine_root, 'veg_litterfall_t_fine_root_nitrogen', mem%litterfall_t%fine_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%fine_root, 'veg_litterfall_t_fine_root_phosphorus', mem%litterfall_t%fine_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%fine_root, 'veg_litterfall_t_fine_root_carbon13', mem%litterfall_t%fine_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%fine_root, 'veg_litterfall_t_fine_root_carbon14', mem%litterfall_t%fine_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%fine_root, 'veg_litterfall_t_fine_root_nitrogen15', mem%litterfall_t%fine_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall_t%coarse_root, 'veg_litterfall_t_coarse_root_carbon', mem%litterfall_t%coarse_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%coarse_root, 'veg_litterfall_t_coarse_root_nitrogen', mem%litterfall_t%coarse_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%coarse_root, 'veg_litterfall_t_coarse_root_phosphorus', mem%litterfall_t%coarse_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%coarse_root, 'veg_litterfall_t_coarse_root_carbon13', mem%litterfall_t%coarse_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%coarse_root, 'veg_litterfall_t_coarse_root_carbon14', mem%litterfall_t%coarse_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%coarse_root, 'veg_litterfall_t_coarse_root_nitrogen15', mem%litterfall_t%coarse_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall_t%sap_wood, 'veg_litterfall_t_sap_wood_carbon', mem%litterfall_t%sap_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%sap_wood, 'veg_litterfall_t_sap_wood_nitrogen', mem%litterfall_t%sap_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%sap_wood, 'veg_litterfall_t_sap_wood_phosphorus', mem%litterfall_t%sap_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%sap_wood, 'veg_litterfall_t_sap_wood_carbon13', mem%litterfall_t%sap_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%sap_wood, 'veg_litterfall_t_sap_wood_carbon14', mem%litterfall_t%sap_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%sap_wood, 'veg_litterfall_t_sap_wood_nitrogen15', mem%litterfall_t%sap_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall_t%heart_wood, 'veg_litterfall_t_heart_wood_carbon', mem%litterfall_t%heart_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%heart_wood, 'veg_litterfall_t_heart_wood_nitrogen', mem%litterfall_t%heart_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%heart_wood, 'veg_litterfall_t_heart_wood_phosphorus', mem%litterfall_t%heart_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%heart_wood, 'veg_litterfall_t_heart_wood_carbon13', mem%litterfall_t%heart_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%heart_wood, 'veg_litterfall_t_heart_wood_carbon14', mem%litterfall_t%heart_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%heart_wood, 'veg_litterfall_t_heart_wood_nitrogen15', mem%litterfall_t%heart_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall_t%labile, 'veg_litterfall_t_labile_carbon', mem%litterfall_t%labile%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%labile, 'veg_litterfall_t_labile_nitrogen', mem%litterfall_t%labile%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%labile, 'veg_litterfall_t_labile_phosphorus', mem%litterfall_t%labile%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%labile, 'veg_litterfall_t_labile_carbon13', mem%litterfall_t%labile%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%labile, 'veg_litterfall_t_labile_carbon14', mem%litterfall_t%labile%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%labile, 'veg_litterfall_t_labile_nitrogen15', mem%litterfall_t%labile%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall_t%reserve, 'veg_litterfall_t_reserve_carbon', mem%litterfall_t%reserve%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%reserve, 'veg_litterfall_t_reserve_nitrogen', mem%litterfall_t%reserve%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%reserve, 'veg_litterfall_t_reserve_phosphorus', mem%litterfall_t%reserve%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%reserve, 'veg_litterfall_t_reserve_carbon13', mem%litterfall_t%reserve%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%reserve, 'veg_litterfall_t_reserve_carbon14', mem%litterfall_t%reserve%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%reserve, 'veg_litterfall_t_reserve_nitrogen15', mem%litterfall_t%reserve%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall_t%fruit, 'veg_litterfall_t_fruit_carbon', mem%litterfall_t%fruit%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%fruit, 'veg_litterfall_t_fruit_nitrogen', mem%litterfall_t%fruit%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%fruit, 'veg_litterfall_t_fruit_phosphorus', mem%litterfall_t%fruit%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%fruit, 'veg_litterfall_t_fruit_carbon13', mem%litterfall_t%fruit%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%fruit, 'veg_litterfall_t_fruit_carbon14', mem%litterfall_t%fruit%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%fruit, 'veg_litterfall_t_fruit_nitrogen15', mem%litterfall_t%fruit%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall_t%seed_bed, 'veg_litterfall_t_seed_bed_carbon', mem%litterfall_t%seed_bed%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%seed_bed, 'veg_litterfall_t_seed_bed_nitrogen', mem%litterfall_t%seed_bed%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%seed_bed, 'veg_litterfall_t_seed_bed_phosphorus', mem%litterfall_t%seed_bed%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%seed_bed, 'veg_litterfall_t_seed_bed_carbon13', mem%litterfall_t%seed_bed%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%seed_bed, 'veg_litterfall_t_seed_bed_carbon14', mem%litterfall_t%seed_bed%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall_t%seed_bed, 'veg_litterfall_t_seed_bed_nitrogen15', mem%litterfall_t%seed_bed%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%change_pool%leaf, 'veg_change_pool_leaf_carbon', mem%change_pool%leaf%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%leaf, 'veg_change_pool_leaf_nitrogen', mem%change_pool%leaf%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%leaf, 'veg_change_pool_leaf_phosphorus', mem%change_pool%leaf%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%leaf, 'veg_change_pool_leaf_carbon13', mem%change_pool%leaf%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%leaf, 'veg_change_pool_leaf_carbon14', mem%change_pool%leaf%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%leaf, 'veg_change_pool_leaf_nitrogen15', mem%change_pool%leaf%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%change_pool%fine_root, 'veg_change_pool_fine_root_carbon', mem%change_pool%fine_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%fine_root, 'veg_change_pool_fine_root_nitrogen', mem%change_pool%fine_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%fine_root, 'veg_change_pool_fine_root_phosphorus', mem%change_pool%fine_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%fine_root, 'veg_change_pool_fine_root_carbon13', mem%change_pool%fine_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%fine_root, 'veg_change_pool_fine_root_carbon14', mem%change_pool%fine_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%fine_root, 'veg_change_pool_fine_root_nitrogen15', mem%change_pool%fine_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%change_pool%coarse_root, 'veg_change_pool_coarse_root_carbon', mem%change_pool%coarse_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%coarse_root, 'veg_change_pool_coarse_root_nitrogen', mem%change_pool%coarse_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%coarse_root, 'veg_change_pool_coarse_root_phosphorus', mem%change_pool%coarse_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%coarse_root, 'veg_change_pool_coarse_root_carbon13', mem%change_pool%coarse_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%coarse_root, 'veg_change_pool_coarse_root_carbon14', mem%change_pool%coarse_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%coarse_root, 'veg_change_pool_coarse_root_nitrogen15', mem%change_pool%coarse_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%change_pool%sap_wood, 'veg_change_pool_sap_wood_carbon', mem%change_pool%sap_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%sap_wood, 'veg_change_pool_sap_wood_nitrogen', mem%change_pool%sap_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%sap_wood, 'veg_change_pool_sap_wood_phosphorus', mem%change_pool%sap_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%sap_wood, 'veg_change_pool_sap_wood_carbon13', mem%change_pool%sap_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%sap_wood, 'veg_change_pool_sap_wood_carbon14', mem%change_pool%sap_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%sap_wood, 'veg_change_pool_sap_wood_nitrogen15', mem%change_pool%sap_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%change_pool%heart_wood, 'veg_change_pool_heart_wood_carbon', mem%change_pool%heart_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%heart_wood, 'veg_change_pool_heart_wood_nitrogen', mem%change_pool%heart_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%heart_wood, 'veg_change_pool_heart_wood_phosphorus', mem%change_pool%heart_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%heart_wood, 'veg_change_pool_heart_wood_carbon13', mem%change_pool%heart_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%heart_wood, 'veg_change_pool_heart_wood_carbon14', mem%change_pool%heart_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%heart_wood, 'veg_change_pool_heart_wood_nitrogen15', mem%change_pool%heart_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%change_pool%labile, 'veg_change_pool_labile_carbon', mem%change_pool%labile%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%labile, 'veg_change_pool_labile_nitrogen', mem%change_pool%labile%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%labile, 'veg_change_pool_labile_phosphorus', mem%change_pool%labile%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%labile, 'veg_change_pool_labile_carbon13', mem%change_pool%labile%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%labile, 'veg_change_pool_labile_carbon14', mem%change_pool%labile%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%labile, 'veg_change_pool_labile_nitrogen15', mem%change_pool%labile%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%change_pool%reserve, 'veg_change_pool_reserve_carbon', mem%change_pool%reserve%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%reserve, 'veg_change_pool_reserve_nitrogen', mem%change_pool%reserve%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%reserve, 'veg_change_pool_reserve_phosphorus', mem%change_pool%reserve%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%reserve, 'veg_change_pool_reserve_carbon13', mem%change_pool%reserve%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%reserve, 'veg_change_pool_reserve_carbon14', mem%change_pool%reserve%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%reserve, 'veg_change_pool_reserve_nitrogen15', mem%change_pool%reserve%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%change_pool%fruit, 'veg_change_pool_fruit_carbon', mem%change_pool%fruit%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%fruit, 'veg_change_pool_fruit_nitrogen', mem%change_pool%fruit%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%fruit, 'veg_change_pool_fruit_phosphorus', mem%change_pool%fruit%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%fruit, 'veg_change_pool_fruit_carbon13', mem%change_pool%fruit%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%fruit, 'veg_change_pool_fruit_carbon14', mem%change_pool%fruit%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%fruit, 'veg_change_pool_fruit_nitrogen15', mem%change_pool%fruit%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%change_pool%seed_bed, 'veg_change_pool_seed_bed_carbon', mem%change_pool%seed_bed%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%seed_bed, 'veg_change_pool_seed_bed_nitrogen', mem%change_pool%seed_bed%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%seed_bed, 'veg_change_pool_seed_bed_phosphorus', mem%change_pool%seed_bed%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%seed_bed, 'veg_change_pool_seed_bed_carbon13', mem%change_pool%seed_bed%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%seed_bed, 'veg_change_pool_seed_bed_carbon14', mem%change_pool%seed_bed%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%change_pool%seed_bed, 'veg_change_pool_seed_bed_nitrogen15', mem%change_pool%seed_bed%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%total_flux%leaf, 'veg_total_flux_leaf_carbon', mem%total_flux%leaf%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%leaf, 'veg_total_flux_leaf_nitrogen', mem%total_flux%leaf%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%leaf, 'veg_total_flux_leaf_phosphorus', mem%total_flux%leaf%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%leaf, 'veg_total_flux_leaf_carbon13', mem%total_flux%leaf%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%leaf, 'veg_total_flux_leaf_carbon14', mem%total_flux%leaf%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%leaf, 'veg_total_flux_leaf_nitrogen15', mem%total_flux%leaf%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%total_flux%fine_root, 'veg_total_flux_fine_root_carbon', mem%total_flux%fine_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%fine_root, 'veg_total_flux_fine_root_nitrogen', mem%total_flux%fine_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%fine_root, 'veg_total_flux_fine_root_phosphorus', mem%total_flux%fine_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%fine_root, 'veg_total_flux_fine_root_carbon13', mem%total_flux%fine_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%fine_root, 'veg_total_flux_fine_root_carbon14', mem%total_flux%fine_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%fine_root, 'veg_total_flux_fine_root_nitrogen15', mem%total_flux%fine_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%total_flux%coarse_root, 'veg_total_flux_coarse_root_carbon', mem%total_flux%coarse_root%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%coarse_root, 'veg_total_flux_coarse_root_nitrogen', mem%total_flux%coarse_root%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%coarse_root, 'veg_total_flux_coarse_root_phosphorus', mem%total_flux%coarse_root%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%coarse_root, 'veg_total_flux_coarse_root_carbon13', mem%total_flux%coarse_root%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%coarse_root, 'veg_total_flux_coarse_root_carbon14', mem%total_flux%coarse_root%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%coarse_root, 'veg_total_flux_coarse_root_nitrogen15', mem%total_flux%coarse_root%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%total_flux%sap_wood, 'veg_total_flux_sap_wood_carbon', mem%total_flux%sap_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%sap_wood, 'veg_total_flux_sap_wood_nitrogen', mem%total_flux%sap_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%sap_wood, 'veg_total_flux_sap_wood_phosphorus', mem%total_flux%sap_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%sap_wood, 'veg_total_flux_sap_wood_carbon13', mem%total_flux%sap_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%sap_wood, 'veg_total_flux_sap_wood_carbon14', mem%total_flux%sap_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%sap_wood, 'veg_total_flux_sap_wood_nitrogen15', mem%total_flux%sap_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%total_flux%heart_wood, 'veg_total_flux_heart_wood_carbon', mem%total_flux%heart_wood%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%heart_wood, 'veg_total_flux_heart_wood_nitrogen', mem%total_flux%heart_wood%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%heart_wood, 'veg_total_flux_heart_wood_phosphorus', mem%total_flux%heart_wood%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%heart_wood, 'veg_total_flux_heart_wood_carbon13', mem%total_flux%heart_wood%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%heart_wood, 'veg_total_flux_heart_wood_carbon14', mem%total_flux%heart_wood%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%heart_wood, 'veg_total_flux_heart_wood_nitrogen15', mem%total_flux%heart_wood%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%total_flux%labile, 'veg_total_flux_labile_carbon', mem%total_flux%labile%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%labile, 'veg_total_flux_labile_nitrogen', mem%total_flux%labile%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%labile, 'veg_total_flux_labile_phosphorus', mem%total_flux%labile%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%labile, 'veg_total_flux_labile_carbon13', mem%total_flux%labile%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%labile, 'veg_total_flux_labile_carbon14', mem%total_flux%labile%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%labile, 'veg_total_flux_labile_nitrogen15', mem%total_flux%labile%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%total_flux%reserve, 'veg_total_flux_reserve_carbon', mem%total_flux%reserve%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%reserve, 'veg_total_flux_reserve_nitrogen', mem%total_flux%reserve%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%reserve, 'veg_total_flux_reserve_phosphorus', mem%total_flux%reserve%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%reserve, 'veg_total_flux_reserve_carbon13', mem%total_flux%reserve%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%reserve, 'veg_total_flux_reserve_carbon14', mem%total_flux%reserve%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%reserve, 'veg_total_flux_reserve_nitrogen15', mem%total_flux%reserve%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%total_flux%fruit, 'veg_total_flux_fruit_carbon', mem%total_flux%fruit%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%fruit, 'veg_total_flux_fruit_nitrogen', mem%total_flux%fruit%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%fruit, 'veg_total_flux_fruit_phosphorus', mem%total_flux%fruit%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%fruit, 'veg_total_flux_fruit_carbon13', mem%total_flux%fruit%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%fruit, 'veg_total_flux_fruit_carbon14', mem%total_flux%fruit%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%fruit, 'veg_total_flux_fruit_nitrogen15', mem%total_flux%fruit%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%total_flux%seed_bed, 'veg_total_flux_seed_bed_carbon', mem%total_flux%seed_bed%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%seed_bed, 'veg_total_flux_seed_bed_nitrogen', mem%total_flux%seed_bed%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%seed_bed, 'veg_total_flux_seed_bed_phosphorus', mem%total_flux%seed_bed%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%seed_bed, 'veg_total_flux_seed_bed_carbon13', mem%total_flux%seed_bed%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%seed_bed, 'veg_total_flux_seed_bed_carbon14', mem%total_flux%seed_bed%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_flux%seed_bed, 'veg_total_flux_seed_bed_nitrogen15', mem%total_flux%seed_bed%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%establishment, 'veg_establishment_carbon', mem%establishment%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_nitrogen', mem%establishment%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_phosphorus', mem%establishment%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_carbon13', mem%establishment%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_carbon14', mem%establishment%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_nitrogen15', mem%establishment%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_carbon', mem%reserve_use%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_nitrogen', mem%reserve_use%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_phosphorus', mem%reserve_use%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_carbon13', mem%reserve_use%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_carbon14', mem%reserve_use%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_nitrogen15', mem%reserve_use%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%init_pool, 'veg_init_pool_carbon', mem%init_pool%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%init_pool, 'veg_init_pool_nitrogen', mem%init_pool%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%init_pool, 'veg_init_pool_phosphorus', mem%init_pool%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%init_pool, 'veg_init_pool_carbon13', mem%init_pool%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%init_pool, 'veg_init_pool_carbon14', mem%init_pool%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%init_pool, 'veg_init_pool_nitrogen15', mem%init_pool%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%total_biomass, 'veg_total_biomass_carbon', mem%total_biomass%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_biomass, 'veg_total_biomass_nitrogen', mem%total_biomass%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_biomass, 'veg_total_biomass_phosphorus', mem%total_biomass%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_biomass, 'veg_total_biomass_carbon13', mem%total_biomass%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_biomass, 'veg_total_biomass_carbon14', mem%total_biomass%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%total_biomass, 'veg_total_biomass_nitrogen15', mem%total_biomass%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue

    END IF

  END SUBROUTINE Init_veg_memory

END MODULE mo_veg_memory_class
