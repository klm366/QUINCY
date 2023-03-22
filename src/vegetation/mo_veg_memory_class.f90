!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Contains the memory class for the vegetation processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20160830, based on the JSBACH4, MPI-M Hamburg)
!! code additions for bgcm store solution by Julia Nabel, MPI BGC (20230117)
!!
MODULE mo_veg_memory_class
#ifndef __NO_JSBACH__

  USE mo_veg_memory_class_util   ! "TYPE bgc_material" code, e.g., t_veg_bgcm_with_components & t_veg_bgcm_with_elements

  USE mo_kind,                   ONLY: wp
#ifdef __QUINCY_STANDALONE__
#else
  USE mo_exception,              ONLY: finish, message
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_utils_iface,        ONLY: assign_if_present_allocatable
#endif

  USE mo_jsb_process_class,      ONLY: VEG_
  USE mo_lnd_bgcm_class,         ONLY: ELEM_C_ID, ELEM_N_ID, ELEM_P_ID, ELEM_C13_ID, ELEM_C14_ID, ELEM_N15_ID, &
    &                                  LAST_ELEM_ID, get_name_for_element_id
  USE mo_lnd_bgcm_store_class,   ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_TOTAL_BIO_ID
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
#ifdef __QUINCY_STANDALONE__
  USE mo_jsb_var_class,          ONLY: t_jsb_struct, t_jsb_n_struct
#else
  USE mo_jsb_pool_class,         ONLY: ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15, t_jsb_pool
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, VEG_TYPE
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_veg_memory, max_no_of_vars
#ifdef __QUINCY_STANDALONE__
  PUBLIC :: max_no_of_structs, max_no_of_n_structs
  PUBLIC :: t_veg_element_a, t_veg_pool_a                                         ! mo_veg_memory_class_util
  PUBLIC :: Init_veg_pools_elements_a, &                                          ! mo_veg_memory_class_util
            Deallocate_var_veg_pools_elements_a                                   ! mo_veg_memory_class_util
  PUBLIC :: veg_pool, veg_growth, veg_litterfall, veg_falloc, veg_total_flux      ! mo_veg_memory_class_util
  PUBLIC :: veg_exudation, veg_establishment, veg_reserve_use, veg_total_biomass  ! mo_veg_memory_class_util
  INTEGER, PARAMETER :: max_no_of_vars      = 250  !< TYPE t_jsb_var
  INTEGER, PARAMETER :: max_no_of_structs   = 20   !< TYPE t_jsb_struct
  INTEGER, PARAMETER :: max_no_of_n_structs = 20   !< TYPE t_jsb_n_struct
#else
  PUBLIC :: t_veg_main_pool, t_veg_tissues_pool, t_veg_elements_pool
  PUBLIC :: t_veg_bgcm_main, t_veg_bgcm_with_components, t_veg_bgcm_with_elements
  ! INTEGER, PARAMETER :: max_no_of_vars      = 530    !< TYPE t_jsb_var | 'normal' var: ~ 220 | 5 pools: 270 var | 4 tissues: 24
  INTEGER, PARAMETER :: max_no_of_vars      = 999  !< only temporary for adding bgc_material in parallel to "old pools"
#endif

  ! ======================================================================================================= !
  !>Type definition for veg memory
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_veg_memory

#ifdef __QUINCY_STANDALONE__
    !------------------------------------------------------------------------------------------------------ !
    !>QUINCY element-pool structures
    !> 
    TYPE(t_quincy_struct_veg_element) :: &
                                      exudation       , &   !< flux between fine roots and mycorrhiza (mol m-2 timestep-1)
                                      establishment   , &   !< flux between seed bed and reserve (mol m-2 timestep-1)
                                      reserve_use     , &   !< net flux between reserve and labile pool (mol m-2 timestep-1)
                                      total_biomass         !< per time step sum over veg\%pool\%..

    TYPE(t_quincy_struct_veg_pool)    :: &   ! nested pool TYPE t_jsb_n_struct
                                      pool            , &   !< pools of the plant
                                      growth          , &   !< fluxes: growth
                                      litterfall      , &   !< fluxes: litterfall
                                      falloc          , &   !< mo_veg_growth:calc_allocation_fraction -> fractional allocation to pools
                                      total_flux            !< for diagnoses of mass conservation: sum of fluxes
#else
    !------------------------------------------------------------------------------------------------------ !
    !>JSBACH4 element-pool structures
    !> 
    ! the veg_main_pool does contain all other pools & fluxes 
    TYPE(t_veg_main_pool),  POINTER :: veg_main_pool
    TYPE(t_veg_bgcm_main),  POINTER :: veg_bgcm_main
#endif

    !------------------------------------------------------------------------------------------------------ !
    !>veg_pool diagnostics - summary of C, N, P, C13, C14, N15 across the pool
    !> 
    TYPE(t_jsb_var_real2d)          :: &   
                                    veg_pool_total_c      , &         !< sum of all veg_pool%*%carbon      [mol C / m2]
                                    veg_pool_total_n      , &         !< sum of all veg_pool%*%nitrogen    [mol C / m2]
                                    veg_pool_total_p      , &         !< sum of all veg_pool%*%phosphorus  [mol C / m2]
                                    veg_pool_total_c13    , &         !< sum of all veg_pool%*%carbon13
                                    veg_pool_total_c14    , &         !< sum of all veg_pool%*%carbon14
                                    veg_pool_total_n15                !< sum of all veg_pool%*%nitrogen15

    !------------------------------------------------------------------------------------------------------ !
    !>basic plant 2D
    !>
    TYPE(t_jsb_var_real2d)          :: &   
                                    height            , &            !< plant height [m]
                                    diameter          , &            !< plant diameter [m]
                                    dens_ind          , &            !< density of individuals [# m-2]
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
    TYPE(t_jsb_var_real3d)          :: &    ! DIMENSION(soil_layer_sb)
                                    root_fraction_sl                 !< fractions of root per layer [unitless]

    !------------------------------------------------------------------------------------------------------ !
    !>target C:N and N:P ratios for various non-leaf tissues (mol/mol)
    !>
    TYPE(t_jsb_var_real2d)          :: &
                                    target_cn_fine_root       , &   !< 
                                    target_cn_coarse_root     , &   !< 
                                    target_cn_sap_wood        , &   !< 
                                    target_cn_fruit           , &   !< 
                                    target_np_fine_root       , &   !< 
                                    target_np_coarse_root     , &   !< 
                                    target_np_sap_wood        , &   !< 
                                    target_np_fruit                 !< 

    !------------------------------------------------------------------------------------------------------ !
    !>canopy layer information
    !>  need 'lrestart = .TRUE.' because lai and lai_cl (PHENO_) are also 'lrestart = .TRUE.'
    !>
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    leaf_nitrogen_cl  , &           !< N in canopy layer [mmol m-2 LAI]
                                    fn_chl_cl         , &           !< fraction of N in chlorophyll, pigments only [unitless]
                                    fn_et_cl          , &           !< fraction of N in electron transport [unitless]
                                    fn_rub_cl         , &           !< fraction of N in Rubisco [unitless]
                                    fn_pepc_cl        , &           !< fraction of N in PEP carboylase, C4 plants only [unitless]
                                    fn_oth_cl         , &           !< fraction of N not associated with photosynthesis [unitless]
                                    fleaf_sunlit_cl                 !< fraction of layer which is sunlit [unitless]

    !------------------------------------------------------------------------------------------------------ !
    !>fluxes
    !>
    TYPE(t_jsb_var_real2d)          :: &   
                                    maint_respiration_pot         , & !< potential maintenance respiration [micro-mol CO2 m-2 s-1]
                                    maint_respiration             , & !< maintenance respiration [micro-mol CO2 m-2 s-1]
                                    maint_respiration_c13         , & !< maintenance respiration [micro-mol 13CO2 m-2 s-1]
                                    maint_respiration_c14         , & !< maintenance respiration [micro-mol 14CO2 m-2 s-1]
                                    growth_respiration            , & !< growth respiration [micro-mol CO2 m-2 s-1]
                                    growth_respiration_c13        , & !< growth respiration [micro-mol 13CO2 m-2 s-1]
                                    growth_respiration_c14        , & !< growth respiration [micro-mol 14CO2 m-2 s-1]
                                    n_transform_respiration       , & !< respiration associated with N transformation [micro-mol CO2 m-2 s-1]
                                    n_fixation_respiration        , & !< respiration associated with N fixation [micro-mol CO2 m-2 s-1]
                                    n_processing_respiration      , & !< respiration associated with N processing [micro-mol CO2 m-2 s-1]
                                    n_processing_respiration_c13  , & !< respiration associated with N processing [micro-mol 13CO2 m-2 s-1]
                                    n_processing_respiration_c14  , & !< respiration associated with N processing [micro-mol 14CO2 m-2 s-1]
                                    npp                           , & !< net primary production [micro-mol CO2 m-2 s-1]
                                    npp_c13                       , & !< net primary production [micro-mol 13CO2 m-2 s-1]
                                    npp_c14                       , & !< net primary production [micro-mol 14CO2 m-2 s-1]
                                    uptake_nh4                    , & !< NH4 uptake from soil [micro-mol N m-2 s-1]
                                    uptake_nh4_n15                , & !< 15NH4 uptake from soil [micro-mol N m-2 s-1]
                                    uptake_no3                    , & !< NO3 uptake from soil [micro-mol N m-2 s-1]
                                    uptake_no3_n15                , & !< 15NO3 uptake from soil [micro-mol N m-2 s-1]
                                    uptake_po4                    , & !< PO4 uptake from soil [micro-mol P m-2 s-1]
                                    n_fixation                    , & !< symbiontic N fixation [micro-mol N m-2 s-1]
                                    n_fixation_n15                , & !< symbiontic N fixation [micro-mol 15N m-2 s-1]
                                    recycling_leaf_n              , & !< net flux of N from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_leaf_n15            , & !< net flux of 15N from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_leaf_p              , & !< net flux of P from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_fine_root_n         , & !< net flux of N from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_fine_root_n15       , & !< net flux of 15N from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_fine_root_p         , & !< net flux of P from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_heart_wood_n        , & !< net flux of N from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_heart_wood_n15      , & !< net flux of 15N from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_heart_wood_p        , & !< net flux of P from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    net_growth                        !< npp minus turnover [micro-mol CO2 m-2 s-1]

    !------------------------------------------------------------------------------------------------------ !
    !>everything else
    !>
    TYPE(t_jsb_var_real2d)          :: &   
                                    fmaint_rate_root          , & !< current N-specific maintentance respiration [~micro-mol C mol-1 N s-1]
                                    unit_npp                  , & !< incremental change in NPP [micro-mol CO2 m-2 s-1]
                                    unit_transpiration        , & !< incremental change in transpiration [kg m-2 s-1]
                                    unit_uptake_n_pot         , & !< incremental change in potential N uptake per unit root carbon [micro-mol N m-2 s-1]
                                    unit_uptake_p_pot         , & !< incremental change in potential P uptake per unit root carbon [micro-mol P m-2 s-1]
                                    unit_uptake_n_act         , & !< incremental change in actual N uptake per unit root carbon [micro-mol N m-2 s-1]
                                    unit_uptake_p_act         , & !< incremental change in actual P uptake per unit root carbon [micro-mol P m-2 s-1]
                                    delta_dens_ind            , & !< change in density of individuals [# m-2 timestep-1]
                                    cost_n_uptake_root        , & !< actual cost for N uptake by roots [mol C mol-1 N]
                                    cost_n_fixation           , & !< actual cost for N uptake through fixation [mol C mol-1 N]
                                    mortality_rate            , & !< background_mort_rate_tree / grass [timestep-1]
                                    k1_opt                    , & !< coefficient for optimal leaf:root ratio
                                    k2_opt                    , & !< coefficient for optimal leaf:root ratio
                                    leaf2sapwood_mass_ratio   , & !< leaf to shoot mass ratio [unitless]
                                    leaf2root_mass_ratio          !< leaf to root mass ratio [unitless]

    !------------------------------------------------------------------------------------------------------ !
    !>daytime averages
    !> 
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_daytime_cl          !< fraction of layer which is sunlit average of the previous day [unitless]
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_daytime             , &    !< air temperature average of the previous day [K]
                                    press_srf_daytime         , &    !< air pressure at surface average of the previous day [Pa]
                                    co2_mixing_ratio_daytime  , &    !< co2 mixing ratio of the air average of the previous day [co2 in ppm]
                                    ga_daytime                , &    !< aerodynamic conductance average of the previous day [m s-1]
                                    beta_sinklim_ps_daytime   , &    !< scaling factor to account for soil moisture constraints on photosynthesis average of the previous day [unitless]
                                    t_jmax_opt_daytime               !< optimal temperature for electron transport of photosynthesis average of the previous day [deg C]

    !------------------------------------------------------------------------------------------------------ !
    !>accumlated over current daytime
    !> 
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_daytime_dacc_cl     !< fraction of layer which is sunlit value accumulated over current day [unitless]
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_daytime_dacc            , & !< air temperature average of the previous day [K]
                                    press_srf_daytime_dacc        , & !< air pressure at surface average of the previous day [Pa]
                                    co2_mixing_ratio_daytime_dacc , & !< co2 mixing ratio of the air average of the previous day [co2 in ppm]
                                    ga_daytime_dacc               , & !< aerodynamic conductance average of the previous day [m s-1]
                                    beta_sinklim_ps_daytime_dacc  , & !< scaling factor to account for soil moisture constraints on photosynthesis average of the previous day [unitless]
                                    t_jmax_opt_daytime_dacc           !< optimal temperature for electron transport of photosynthesis value accumulated over current day [deg C]

    !------------------------------------------------------------------------------------------------------ !
    !>moving averages
    !> 
    ! 1.1 conditions for the labile pool dynamics
    TYPE(t_jsb_var_real2d)          :: &   
                                    gpp_tlabile_mavg               , &  !< gpp  [micro-mol m-2 s-1]
                                    maint_respiration_tlabile_mavg , &  !< maintenance respiration [micro-mol m-2 s-1]
                                    growth_req_n_tlabile_mavg      , &  !< moles N required for a unit of C growth under current allocation fractions time-averaged [mol N m-2 s-1]
                                    growth_req_p_tlabile_mavg           !< moles P required for a unit of N growth under current allocation fractions time-averaged [mol P m-2 s-1]
    ! 1.2 conditions for phenology calculations
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tphen_mavg                    !< air temperature time-averaged [K]
    ! 1.3 memory for plant uptake (average demand over a couple of days to get rid of diurnal C cycle)
    TYPE(t_jsb_var_real2d)          :: &   
                                    npp_tuptake_mavg               , &  !< npp of the last week [micro-mol m-2 s-1]
                                    demand_uptake_n_tuptake_mavg   , &  !< scalar of Vmax for root N uptake time-averaged [unitless]
                                    demand_uptake_p_tuptake_mavg   , &  !< scalar of Vmax for root P uptake time-averaged [unitless]
                                    growth_req_n_tuptake_mavg      , &  !< moles N required for a unit of C growth under current allocation fractions time-averaged [mol N m-2 s-1]
                                    growth_req_p_tuptake_mavg      , &  !< moles P required for a unit of N growth under current allocation fractions time-averaged [mol P m-2 s-1]
                                    cost_n_fixation_talloc_mavg    , &  !< moving average of N fixation cost [mol C mol-1 N]
                                    exudation_c_tmyc_mavg               !< moving average of C allocation to mycorrhizal fungi [mol C]
    ! 1.4 conditions for photosynthesis for the optimisation of foliar N fractions 
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_tfrac_mavg_cl          !< fraction of layer which is sunlit time-averaged [unitless]
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tfrac_mavg               , &  !< air temperature time-averaged [K]
                                    press_srf_tfrac_mavg           , &  !< air pressure at surface time-averaged [Pa]
                                    co2_mixing_ratio_tfrac_mavg    , &  !< co2 mixing ratio of the air time-averaged [co2 in ppm]
                                    ga_tfrac_mavg                  , &  !< aerodynamic conductance time-averaged [m s-1]
                                    beta_sinklim_ps_tfrac_mavg          !< scaling factor to account for soil moisture constraints on photosynthesis time-averaged [unitless]
    ! 1.5 conditions for photosynthesis for the optimisation of foliar N  
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_tcnl_mavg_cl           !< fraction of layer which is sunlit time-averaged [unitless]
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tcnl_mavg                , &  !< air temperature time-averaged [K]
                                    press_srf_tcnl_mavg            , &  !< air pressure at surface time-averaged [Pa]
                                    co2_mixing_ratio_tcnl_mavg     , &  !< co2 mixing ratio of the air time-averaged [co2 in ppm]
                                    ga_tcnl_mavg                   , &  !< aerodynamic conductance time-averaged [m s-1]
                                    beta_sinklim_ps_tcnl_mavg      , &  !< scaling factor to account for soil moisture constraints on conductance time-averaged [unitless]
                                    uptake_n_tcnl_mavg             , &  !< root N uptake time-averaged [micro-mol N m-2 s-1]
                                    growth_cn_tcnl_mavg            , &  !< new growth CN ratio time-averaged [mol C mol-1 N]
                                    growth_np_tcnl_mavg            , &  !< new growth NP ratio time-averaged [mol N mol-1 P]
                                    npp_tcnl_mavg                  , &  !< NPP time-averaged [micro-mol CO2 m-2 s-1]
                                    fmaint_rate_tcnl_mavg          , &  !< maintenance respiration rate per unit N for roots
                                    labile_nitrogen_tcnl_mavg      , &  !< labile nitrogen vegetation pool time-averaged [mol N m-2]
                                    labile_carbon_tcnl_mavg        , &  !< labile carbon vegetation pool time-averaged [mol C m-2]
                                    labile_phosphorus_tcnl_mavg         !< labile phosphorus vegetation pool time-averaged [mol P m-2]
    ! 1.6 conditions for photosynthesis for the optimisation of biomass allocation fractions
    TYPE(t_jsb_var_real2d)          :: &   
                                    npp_talloc_mavg                , &  !< NPP [micro-mol m-2 s-1]
                                    n_fixation_talloc_mavg         , &  !< N fixation time-averaged [micro-mol N m-2 s-1]
                                    transpiration_talloc_mavg      , &  !< transpiration time-averaged [kg m-2 s-1]
                                    unit_npp_talloc_mavg           , &  !< rate of change in NPP with  increment in leaf mass [micro-mol C mol-1]
                                    unit_transpiration_talloc_mavg , &  !< incremental change in transpiration time-averaged [kg m-2 s-1]
                                    unit_uptake_n_talloc_mavg      , &  !< actual rate of N uptkae per unit root [micro-mol N mol-1]
                                    unit_uptake_p_talloc_mavg      , &  !< actual rate of P uptake per unit root [micro-mol P mol-1]
                                    dphi_talloc_mavg               , &  !< leaf-soil water potential gradient [MPa]
                                    growth_cn_talloc_mavg          , &  !< plant growth cn ratio, equivalent to 1/growth_req_n time-averaged [mol C mol-1 N]
                                    growth_cp_talloc_mavg          , &  !< plant growth cp ratio, equivalent to 1/(growthreq_n*growth_req_p) time-averaged [mol C mol-1 P]
                                    growth_np_talloc_mavg          , &  !< plant growth np ratio, equivalent to 1/(growth_req_p) time-averaged [mol N mol-1 P]
                                    labile_carbon_talloc_mavg      , &  !< labile carbon pool needed for allometry calculations time-averaged [mol C m-2]
                                    labile_nitrogen_talloc_mavg    , &  !< labile nitrogen pool neeeded for allometry calculations time-averaged [mol N m-2]
                                    labile_phosphorus_talloc_mavg  , &  !< labile phosphorus pool neeeded for allometry calculations time-averaged [mol P m-2]
                                    reserve_carbon_talloc_mavg     , &  !< reserve carbon needed for optimal allocation time-averaged [mol C m-2]
                                    reserve_nitrogen_talloc_mavg   , &  !< reserve nitrogen needed for optimal allocation time-averaged [mol N m-2]
                                    reserve_phosphorus_talloc_mavg , &  !< reserve phosphorus needed for optimal allocation time-averaged [mol P m-2]
                                    w_root_lim_talloc_mavg              !< root water limitation below which root allocation starts responding time-averaged
    ! 1.7 memory over root life time 
    TYPE(t_jsb_var_real2d)          :: &   
                                    leaf2root_troot_mavg           , &  !< leaf to root mass ratio time-averaged [unitless]
                                    unit_npp_troot_mavg            , &  !< rate of change in NPP with increment in leaf mass [micro-mol C mol-1]
                                    unit_uptake_n_pot_troot_mavg   , &  !< rate of N uptkae per unit root [micro-mol N mol-1]
                                    unit_uptake_p_pot_troot_mavg   , &  !< rate of P uptake per unit root [micro-mol P mol-1]
                                    fmaint_rate_troot_mavg              !< moving average of N-specific maintentance respiration [~micro-mol C mol-1 N s-1]
    ! 1.8 memory over timescale of vegetation dynamics 
    TYPE(t_jsb_var_real2d)          :: &   
                                    an_boc_tvegdyn_mavg            , &  !< net photosynthesis at bottom of canopy, lowest canopy layer existing time-averaged [micro-mol m-2 s-1]
                                    net_growth_tvegdyn_mavg        , &  !< net growth, npp minus turnover [micro-mol C m-2 s-1]
                                    lai_tvegdyn_mavg                    !< lai time-averaged [m2 m-2]
    ! 1.9 variables for respiration accclimation 
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tacclim_mavg             , &  !< air temperature for respiration acclimation time-averaged [K]
                                    t_soil_root_tacclim_mavg            !< root temperature for respiration acclimation time-averaged [deg K]
    ! 1.10 certain variables, temperature, ...
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_week_mavg                , &  !< air temperature time-averaged [K]
                                    t_air_month_mavg               , &  !< air temperature time-averaged [K]
                                    t_jmax_opt_mavg                     !< optimal temperature for electron transport of photosynthesis time-averaged [deg C]

  CONTAINS
    PROCEDURE :: Init => Init_veg_memory
  END TYPE t_veg_memory

#ifdef __QUINCY_STANDALONE__
#else
  ! ======================================================================================================= !
  !>
  !> Type definition for veg_bgc_material_elements
  !>
  TYPE, EXTENDS(t_jsb_pool) :: t_veg_bgcm_with_elements
    TYPE(t_jsb_var_real2d), POINTER :: carbon    
    TYPE(t_jsb_var_real2d), POINTER :: nitrogen  
    TYPE(t_jsb_var_real2d), POINTER :: phosphorus
    TYPE(t_jsb_var_real2d), POINTER :: carbon13  
    TYPE(t_jsb_var_real2d), POINTER :: carbon14  
    TYPE(t_jsb_var_real2d), POINTER :: nitrogen15
  CONTAINS
    PROCEDURE :: Get_element_name_by_id     => veg_bgcm_with_elements_get_element_name_by_id
    PROCEDURE :: Init                       => veg_bgcm_with_elements_init
  END TYPE t_veg_bgcm_with_elements

  ! ======================================================================================================= !
  !>
  !> Type definition for veg_bgc_material_components
  !>
  TYPE, EXTENDS(t_jsb_pool) :: t_veg_bgcm_with_components
    TYPE(t_veg_bgcm_with_elements), POINTER :: leaf       
    TYPE(t_veg_bgcm_with_elements), POINTER :: fine_root  
    TYPE(t_veg_bgcm_with_elements), POINTER :: coarse_root
    TYPE(t_veg_bgcm_with_elements), POINTER :: sap_wood   
    TYPE(t_veg_bgcm_with_elements), POINTER :: heart_wood 
    TYPE(t_veg_bgcm_with_elements), POINTER :: labile     
    TYPE(t_veg_bgcm_with_elements), POINTER :: reserve    
    TYPE(t_veg_bgcm_with_elements), POINTER :: fruit      
    TYPE(t_veg_bgcm_with_elements), POINTER :: seed_bed   
  CONTAINS
    PROCEDURE :: Get_element_name_by_id     => veg_bgcm_with_components_get_element_name_by_id
    PROCEDURE :: Init                       => veg_bgcm_with_components_init
  END TYPE t_veg_bgcm_with_components
  
  ! ======================================================================================================= !
  !>
  !> Type definition for veg_bgc_material_main
  !> the main bgc_material structure containing all VEG_ bgc_material structures
  !> main structure, added directly into the "TYPE(t_jsb_memory) :: t_veg_memory" as "CLASS(t_jsb_pool) :: bgc_material"
  !>
  TYPE, EXTENDS(t_jsb_pool) :: t_veg_bgcm_main
    !! pools
    TYPE(t_veg_bgcm_with_components),  POINTER :: veg_bgcm_pool
    TYPE(t_veg_bgcm_with_elements),    POINTER :: veg_bgcm_total_biomass
    !! fluxes
    ! TYPE(t_veg_bgcm_with_components),  POINTER :: veg_growth    
    ! TYPE(t_veg_bgcm_with_components),  POINTER :: veg_litterfall
    ! TYPE(t_veg_bgcm_with_components),  POINTER :: veg_total_flux
    ! TYPE(t_veg_bgcm_with_components),  POINTER :: veg_falloc
    ! TYPE(t_veg_bgcm_with_elements),    POINTER :: veg_exudation    
    ! TYPE(t_veg_bgcm_with_elements),    POINTER :: veg_establishment
    ! TYPE(t_veg_bgcm_with_elements),    POINTER :: veg_reserve_use  
  CONTAINS
    PROCEDURE :: Get_element_name_by_id     => veg_bgcm_main_get_element_name_by_id
    PROCEDURE :: Init                       => veg_bgcm_main_init
  END TYPE t_veg_bgcm_main
#endif

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the VEG_ process
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_veg_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
#ifdef __QUINCY_STANDALONE__
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2   ! , tables
#else
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
#endif
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model
    dsl4jsb_Use_config(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_veg_memory),  INTENT(inout), TARGET  :: mem             !< veg memory
    CHARACTER(len=*),     INTENT(in)             :: prefix          !< process name
    CHARACTER(len=*),     INTENT(in)             :: suffix          !< tile name
    INTEGER,              INTENT(in)             :: lct_ids(:)      !< vector of land cover type IDs
    INTEGER,              INTENT(in)             :: model_id        !< model ID model\%id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),          POINTER :: model
    TYPE(t_jsb_grid),           POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid),          POINTER :: surface                      ! Vertical grid
    TYPE(t_jsb_vgrid),          POINTER :: vgrid_canopy                 ! Vertical grid
    TYPE(t_jsb_vgrid),          POINTER :: vgrid_soil_sb                ! Vertical grid
    CHARACTER(len=30)                   :: unit_veg_pool
    CHARACTER(len=30)                   :: unit_veg_flux
    INTEGER                             :: elem_idx_map(LAST_ELEM_ID - 1)  !< element mapper ID -> IND
#ifdef __QUINCY_STANDALONE__
#else
    TYPE(t_veg_elements_pool),  POINTER :: veg_elements_pool    ! e.g. veg_exudation & veg_total_biomass
    TYPE(t_veg_tissues_pool),   POINTER :: veg_tissues_pool     ! e.g. veg_pool & veg_growth
    TYPE(t_veg_bgcm_with_components),  POINTER :: veg_bgcm_components  ! bgc_material components
    TYPE(t_veg_bgcm_with_elements),    POINTER :: veg_bgcm_elements    ! bgc_material elements
#endif
    INTEGER                             :: table
    CHARACTER(len=*), PARAMETER         :: routine = TRIM(modname)//':Init_veg_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
#ifdef __QUINCY_STANDALONE__
    table         = 1
#else
    table         =  tables(1)
#endif
    hgrid         => Get_grid(model%grid_id)
    surface       => Get_vgrid('surface')
#ifdef __QUINCY_STANDALONE__
    vgrid_canopy  => Get_vgrid('canopy_layer')
#else
    vgrid_canopy  => Get_vgrid('canopy_layer')
#endif
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(VEG_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! init local variables
    unit_veg_pool = "mol m-2"
    unit_veg_flux = "mol m-2 timestep-1"

    ! ----------------------------------------------------------------------------------------------------- !
    ! add memory only for LAND & PFT tiles
#ifdef __QUINCY_STANDALONE__
#else
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. &
       & One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN
#endif

      ! --------------------------------------------------------------------------------------------------- !
      !>1.0 element-pool structures
      !> 

#ifdef __QUINCY_STANDALONE__
      ! --------------------------------------------------------------------------------------------------- !
      !>  1.1 QUINCY element-pool structures
      !> 

      ! add pools to memory
      CALL mem%Add_pool('exudation', mem%exudation, &
        & hgrid, surface, &
        & t_cf('exudation', TRIM(unit_veg_flux), 'Exudation element in vegetation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .FALSE., &
        & lrestart = .FALSE.)
        
      CALL mem%Add_pool('establishment', mem%establishment, &
        & hgrid, surface, &
        & t_cf('establishment', TRIM(unit_veg_flux), 'Establishment element in vegetation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .FALSE., &
        & lrestart = .FALSE.)
        
      CALL mem%Add_pool('reserve_use', mem%reserve_use, &
        & hgrid, surface, &
        & t_cf('reserve_use', TRIM(unit_veg_flux), 'Reserve_use element in vegetation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .FALSE., &
        & lrestart = .FALSE.)
        
      CALL mem%Add_pool('total_biomass', mem%total_biomass, &
        & hgrid, surface, &
        & t_cf('total_biomass', TRIM(unit_veg_pool), 'Total_biomass element in vegetation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .FALSE.)

      ! add variables to pools
      CALL mem%exudation%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%establishment%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%reserve_use%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%total_biomass%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_pool))

      ! add pools to memory
      CALL mem%Add_pool('pool', mem%pool, &
        & hgrid, surface, &
        & t_cf('pool', TRIM(unit_veg_pool), 'Pool vegetation pool'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE.)
        
      CALL mem%Add_pool('growth', mem%growth, &
        & hgrid, surface, &
        & t_cf('growth', TRIM(unit_veg_flux), 'Growth vegetation flux'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .FALSE., &
        & lrestart = .FALSE.)
        
      CALL mem%Add_pool('litterfall', mem%litterfall, &
        & hgrid, surface, &
        & t_cf('litterfall', TRIM(unit_veg_flux), 'Litterfall vegetation flux'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .FALSE., &
        & lrestart = .FALSE.)
        
      CALL mem%Add_pool('falloc', mem%falloc, &
        & hgrid, surface, &
        & t_cf('falloc', TRIM(unit_veg_flux), 'Falloc vegetation flux'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .FALSE., &
        & lrestart = .FALSE.)
        
      CALL mem%Add_pool('total_flux', mem%total_flux, &
        & hgrid, surface, &
        & t_cf('total_flux', TRIM(unit_veg_flux), 'Total_flux vegetation flux'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .FALSE., &
        & lrestart = .FALSE.)
            
      ! add pools to pools
      CALL mem%pool%Add_elements(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_pool))
      CALL mem%growth%Add_elements(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%litterfall%Add_elements(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%falloc%Add_elements(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))
      CALL mem%total_flux%Add_elements(prefix, suffix, lct_ids, model_id, TRIM(unit_veg_flux))

      ! add variables to pools
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
#else
      ! --------------------------------------------------------------------------------------------------- !
      !>  1.2 ICON-Land bgc_material structures
      !> 
      ! copy elements_index_map from config to local
      elem_idx_map(:) = model%config%elements_index_map(:)

      ! --------------------------------------------------------------------------------------------------- !
      !> create and populate the main bgc_material structure mem%bgc_material that contains all other bgc_material
      !> note: no spaces in 'name' and/or 'shortname'
      !>
      ALLOCATE(t_veg_bgcm_main :: mem%bgc_material)
      CALL mem%bgc_material%Init(name = 'vegetation_bgc_material_main', shortname = 'veg_bgcm_main', &
        &                        elements_index_map = elem_idx_map(:), l_elements = .FALSE.)
      SELECT TYPE (selected => mem%bgc_material)
      CLASS IS (t_veg_bgcm_main)
        mem%veg_bgcm_main => selected
      END SELECT

      ! --------------------------------------------------------------------------------------------------- !
      !> Add the info to this memory that this memory contains bgc materials
      !>
      mem%has_bgc_materials = .TRUE.

      ! --------------------------------------------------------------------------------------------------- !
      !> create the bgc_material components structures
      !>
      ! "shortname" is used for model output names
      ! veg_pool
      veg_bgcm_components => Create_veg_bgc_material_components( &
        & name = 'vegetation_bgcm_pool', shortname = 'veg_bgcm_pool', unit = unit_veg_pool, &
        & lct_ids = lct_ids, elements_index_map = elem_idx_map(:), &
        & pool_suffix = 'bgcm')                                     ! create the bgcm object (l_elements=.TRUE. by default)
      CALL mem%veg_bgcm_main%Add_pool(veg_bgcm_components, VEG_BGCM_POOL_ID) ! add bgcm object to pool_list(:)
      mem%veg_bgcm_main%veg_bgcm_pool => veg_bgcm_components        ! link to bgcm object in pool_list(:)

      ! --------------------------------------------------------------------------------------------------- !
      !> create the bgc_material elements structures
      !>
      ! "shortname" is used for model output names
      ! veg_total_biomass
      ALLOCATE(veg_bgcm_elements)
      CALL veg_bgcm_elements%Init(name = 'vegetation_bgcm_total_biomass', shortname = 'veg_bgcm_total_biomass', &
        &                        elements_index_map = elem_idx_map(:), &
        &                        l_elements = .TRUE., element_unit = unit_veg_pool)
      CALL mem%veg_bgcm_main%Add_pool(veg_bgcm_elements, VEG_BGCM_TOTAL_BIO_ID)
      mem%veg_bgcm_main%veg_bgcm_total_biomass => veg_bgcm_elements
      NULLIFY(veg_bgcm_elements)
      
      ! create list of CHARACTERs with names veg_bgcm_main to the "lowest childs" of bgcm structures
      CALL mem%bgc_material%Set_paths()

      ! add variables to bgc_material
      CALL mem%Add_var(mem%veg_bgcm_main%veg_bgcm_pool, &
        & hgrid, surface, vgrid_canopy,                                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = BASIC,                                                                        &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .TRUE.,                                                                           &
        ! & initval_r = 0.0_wp)
        & initval_r = 7.2_wp)     ! @TODO ONLY FOR TESTING !!

      CALL mem%Add_var(mem%veg_bgcm_main%veg_bgcm_total_biomass, &
        & hgrid, surface, vgrid_canopy,                                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = MEDIUM,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .TRUE.,                                                                           &
        ! & initval_r = 0.0_wp)
        & initval_r = 8.3_wp)     ! @TODO ONLY FOR TESTING !!


      !! ====================================================================================================
      !! the old "pools" initialization starts here
      !! ====================================================================================================


      ! --------------------------------------------------------------------------------------------------- !
      !> populate all_elements(:) with the element_id of all elements
      !> QUINCY is using always all of the below six elements
      !>
      IF (.NOT. ALLOCATED(all_elements)) THEN
        ALLOCATE(all_elements(6))
        all_elements(:) = [ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15]
      END IF

      ! --------------------------------------------------------------------------------------------------- !
      ! create and populate the main pool structure mem%pools that contains all other pools
      ! note: no spaces in 'name' and/or 'short name'
      ALLOCATE(t_veg_main_pool :: mem%pools)
      CALL mem%pools%Init(name='vegetation_main_pool', shortname='veg_main_pool', &
        &                 elements_index_map=elem_idx_map(:), l_elements=.FALSE.)

      SELECT TYPE (selected => mem%pools)
      CLASS IS (t_veg_main_pool)
        mem%veg_main_pool => selected
      END SELECT

      ! --------------------------------------------------------------------------------------------------- !
      ! create the vegetation tissue pools & fluxes (veg_tissues_pool) 
      ! veg_pool / veg_growth / ...
      !   create pools & fluxes
      !   "shortname" is used for model output names
      ! veg_pool
      veg_tissues_pool => Create_veg_tissues_pool(name='vegetation_pool', shortname='veg_pool', unit='mol m-2', &
        & lct_ids=lct_ids, elements_index_map=elem_idx_map(:), &
        &  pool_suffix='vegpool') ! l_elements=.TRUE. by default
      CALL mem%veg_main_pool%Add_pool(veg_tissues_pool)
      mem%veg_main_pool%veg_pool => veg_tissues_pool
      ! veg_growth
      veg_tissues_pool => Create_veg_tissues_pool(name='veg_growth_flux', shortname='veg_growth', unit='mol m-2', &
        & lct_ids=lct_ids, elements_index_map=elem_idx_map(:), &
        &  pool_suffix='groflux') ! l_elements=.TRUE. by default
      CALL mem%veg_main_pool%Add_pool(veg_tissues_pool)
      mem%veg_main_pool%veg_growth => veg_tissues_pool
      ! veg_litterfall
      veg_tissues_pool => Create_veg_tissues_pool(name='veg_litterfall_flux', shortname='veg_litterfall', unit='mol m-2', &
        & lct_ids=lct_ids, elements_index_map=elem_idx_map(:), &
        &  pool_suffix='litflux') ! l_elements=.TRUE. by default
      CALL mem%veg_main_pool%Add_pool(veg_tissues_pool)
      mem%veg_main_pool%veg_litterfall => veg_tissues_pool
      ! veg_total_flux
      veg_tissues_pool => Create_veg_tissues_pool(name='vegetation_total_flux', shortname='veg_total_flux', unit='mol m-2', &
        & lct_ids=lct_ids, elements_index_map=elem_idx_map(:), &
        &  pool_suffix='totflux') ! l_elements=.TRUE. by default
      CALL mem%veg_main_pool%Add_pool(veg_tissues_pool)
      mem%veg_main_pool%veg_total_flux => veg_tissues_pool
      ! veg_falloc
      veg_tissues_pool => Create_veg_tissues_pool(name='veg_allocation_fraction_flux', shortname='veg_falloc', unit='mol m-2', &
        & lct_ids=lct_ids, elements_index_map=elem_idx_map(:), &
        &  pool_suffix='falflux') ! l_elements=.TRUE. by default
      CALL mem%veg_main_pool%Add_pool(veg_tissues_pool)
      mem%veg_main_pool%veg_falloc => veg_tissues_pool
      NULLIFY(veg_tissues_pool)
      
      ! init and add element pools to the veg_main_pool
      ! (i.e., add them at the same level as the above veg_tissues_pool veg_pool / veg_growth / ...)
      ! "shortname" is used for model output names
      ! veg_total_biomass
      ALLOCATE(veg_elements_pool)
      CALL veg_elements_pool%Init(name='vegetation_total_biomass', shortname='veg_total_biomass', &
        &                         elements_index_map=elem_idx_map(:),                             &
        &                         l_elements=.TRUE., element_unit='mol m-2')
      CALL mem%veg_main_pool%Add_pool(veg_elements_pool)
      mem%veg_main_pool%veg_total_biomass => veg_elements_pool
      NULLIFY(veg_elements_pool)
      ! veg_exudation
      ALLOCATE(veg_elements_pool)
      CALL veg_elements_pool%Init(name='veg_exudation', shortname='veg_exudation', &
        &                         elements_index_map=elem_idx_map(:),              &
        &                         l_elements=.TRUE., element_unit='mol m-2')
      CALL mem%veg_main_pool%Add_pool(veg_elements_pool)
      mem%veg_main_pool%veg_exudation => veg_elements_pool
      NULLIFY(veg_elements_pool)
      ! veg_establishment
      ALLOCATE(veg_elements_pool)
      CALL veg_elements_pool%Init(name='veg_establishment', shortname='veg_establishment', &
        &                         elements_index_map=elem_idx_map(:),                      &
        &                         l_elements=.TRUE., element_unit='mol m-2')
      CALL mem%veg_main_pool%Add_pool(veg_elements_pool)
      mem%veg_main_pool%veg_establishment => veg_elements_pool
      NULLIFY(veg_elements_pool)
      ! veg_reserve_use
      ALLOCATE(veg_elements_pool)
      CALL veg_elements_pool%Init(name='veg_reserve_use', shortname='veg_reserve_use', &
        &                         elements_index_map=elem_idx_map(:),                  &
        &                         l_elements=.TRUE., element_unit='mol m-2')
      CALL mem%veg_main_pool%Add_pool(veg_elements_pool)
      mem%veg_main_pool%veg_reserve_use => veg_elements_pool
      NULLIFY(veg_elements_pool)

      ! set paths of all above initialized & added pools (type-bound procedure of t_jsb_pool)
      CALL mem%pools%Set_paths()

      ! add all variables in the veg_tissues_pool to the VEG_ memory
      CALL mem%Add_var(mem%veg_main_pool%veg_pool, &
        & hgrid, surface, vgrid_canopy,                                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = BASIC,                                                                        &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .TRUE.,                                                                           &
        & initval_r = 0.0_wp )
      CALL mem%Add_var(mem%veg_main_pool%veg_growth, &
        & hgrid, surface, vgrid_canopy,                                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                         &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp )
      CALL mem%Add_var(mem%veg_main_pool%veg_litterfall, &
        & hgrid, surface, vgrid_canopy,                                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                         &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp )
      CALL mem%Add_var(mem%veg_main_pool%veg_total_flux, &
        & hgrid, surface, vgrid_canopy,                                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                         &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp )
      CALL mem%Add_var(mem%veg_main_pool%veg_falloc, &
        & hgrid, surface, vgrid_canopy,                                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                         &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp )

      ! add all variables in the veg_elements_pool to the VEG_ memory
      CALL mem%Add_var(mem%veg_main_pool%veg_total_biomass, &
        & hgrid, surface, vgrid_canopy,                                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                         &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp )
      CALL mem%Add_var(mem%veg_main_pool%veg_exudation, &
        & hgrid, surface, vgrid_canopy,                                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                         &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp )
      CALL mem%Add_var(mem%veg_main_pool%veg_establishment, &
        & hgrid, surface, vgrid_canopy,                                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                         &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp )
      CALL mem%Add_var(mem%veg_main_pool%veg_reserve_use, &
        & hgrid, surface, vgrid_canopy,                                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                         &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp )

      ! output to logfile all pools & their variables at model-init 
      ! (maybe use for debugging; its quite a long list)
      !IF (.NOT. ASSOCIATED(mem%parent)) THEN
      !  CALL mem%veg_main_pool%Print()
      !  ! CALL mem%veg_main_pool%veg_pool%Print()
      !END IF
#endif

      ! --------------------------------------------------------------------------------------------------- !
      !>2.0 2D & 3D variables
      !> 

      ! diagnostic variables specific to QUINCY pools in ICON-Land
      CALL mem%Add_var('veg_pool_total_c', mem%veg_pool_total_c, &
        & hgrid, surface, &
        & t_cf('veg_pool_total_c', '[mol C / m2]', 'total C, sum of all veg_pool tissues'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_total_n', mem%veg_pool_total_n, &
        & hgrid, surface, &
        & t_cf('veg_pool_total_n', '[mol C / m2]', 'total N, sum of all veg_pool tissues'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_total_p', mem%veg_pool_total_p, &
        & hgrid, surface, &
        & t_cf('veg_pool_total_p', '[mol C / m2]', 'total P, sum of all veg_pool tissues'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
        
      CALL mem%Add_var('veg_pool_total_c13', mem%veg_pool_total_c13, &
        & hgrid, surface, &
        & t_cf('veg_pool_total_c13', '..', 'total C13, sum of all veg_pool tissues'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_total_c14', mem%veg_pool_total_c14, &
        & hgrid, surface, &
        & t_cf('veg_pool_total_c14', '..', 'total C14, sum of all veg_pool tissues'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_total_n15', mem%veg_pool_total_n15, &
        & hgrid, surface, &
        & t_cf('veg_pool_total_n15', '..', 'total N15, sum of all veg_pool tissues'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      ! normal VEG_ variables (not specific to QUINCY pools in ICON-Land)
      CALL mem%Add_var('height', mem%height, &
        & hgrid, surface, &
        & t_cf('height', 'm', 'plant height'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('diameter', mem%diameter, &
        & hgrid, surface, &
        & t_cf('diameter', 'm', 'plant diameter'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('dens_ind', mem%dens_ind, &
        & hgrid, surface, &
        & t_cf('dens_ind', '# m-2', 'density of individuals'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('target_cn_leaf', mem%target_cn_leaf, &
        & hgrid, surface, &
        & t_cf('target_cn_leaf', '', 'target foliar C:N ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('target_np_leaf', mem%target_np_leaf, &
        & hgrid, surface, &
        & t_cf('target_np_leaf', '', 'target foliar N:P ratio '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('mean_leaf_age', mem%mean_leaf_age, &
        & hgrid, surface, &
        & t_cf('mean_leaf_age', '', 'average foliage age in days '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('cohort_age', mem%cohort_age, &
        & hgrid, surface, &
        & t_cf('cohort_age', '', 'age of the cohort in years, used with veg_dynamics_scheme cohort'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('f_n_demand', mem%f_n_demand, &
        & hgrid, surface, &
        & t_cf('f_n_demand', '', 'relative plant demand for N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('f_p_demand', mem%f_p_demand, &
        & hgrid, surface, &
        & t_cf('f_p_demand', '', 'relative plant demand for P'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('beta_sinklim_ps', mem%beta_sinklim_ps, &
        & hgrid, surface, &
        & t_cf('beta_sinklim_ps', '', 'scaling factor to account for a direct sink limitation contraint on photosynthesis'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('dphi', mem%dphi, &
        & hgrid, surface, &
        & t_cf('dphi', 'MPa', 'plant-soil water potential gradient'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('growth_req_n', mem%growth_req_n, &
        & hgrid, surface, &
        & t_cf('growth_req_n', '', 'moles N required for a unit of C growth under current allocation fractions '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('growth_req_p', mem%growth_req_p, &
        & hgrid, surface, &
        & t_cf('growth_req_p', '', 'moles P required for a unit of N growth under current allocation fractions '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('leaf_cn_direction', mem%leaf_cn_direction, &
        & hgrid, surface, &
        & t_cf('leaf_cn_direction', '', 'indicates direction of change in leaf CN ratio to maximise NPP'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('target_lai', mem%target_lai, &
        & hgrid, surface, &
        & t_cf('target_lai', 'm2 m-2', 'LAI implied from sap wood area'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('kstar_labile', mem%kstar_labile, &
        & hgrid, surface, &
        & t_cf('kstar_labile', '', 'labile pool turnover rate, only temperature and moisture constraints, no stoichiometry'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('do_cohort_harvest', mem%do_cohort_harvest, &
        & hgrid, surface, &
        & t_cf('do_cohort_harvest', '', 'logical: cohort harvest on/off, in veg_dynamics_scheme cohort'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('root_fraction_sl', mem%root_fraction_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('root_fraction_sl', '', 'fractions of root per layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('target_cn_fine_root', mem%target_cn_fine_root, &
        & hgrid, surface, &
        & t_cf('target_cn_fine_root', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('target_cn_coarse_root', mem%target_cn_coarse_root, &
        & hgrid, surface, &
        & t_cf('target_cn_coarse_root', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('target_cn_sap_wood', mem%target_cn_sap_wood, &
        & hgrid, surface, &
        & t_cf('target_cn_sap_wood', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('target_cn_fruit', mem%target_cn_fruit, &
        & hgrid, surface, &
        & t_cf('target_cn_fruit', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('target_np_fine_root', mem%target_np_fine_root, &
        & hgrid, surface, &
        & t_cf('target_np_fine_root', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('target_np_coarse_root', mem%target_np_coarse_root, &
        & hgrid, surface, &
        & t_cf('target_np_coarse_root', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('target_np_sap_wood', mem%target_np_sap_wood, &
        & hgrid, surface, &
        & t_cf('target_np_sap_wood', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('target_np_fruit', mem%target_np_fruit, &
        & hgrid, surface, &
        & t_cf('target_np_fruit', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('leaf_nitrogen_cl', mem%leaf_nitrogen_cl, &
        & hgrid, vgrid_canopy, &
        & t_cf('leaf_nitrogen_cl', 'mmol m-2 LAI', 'N in canopy layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('fn_chl_cl', mem%fn_chl_cl, &
        & hgrid, vgrid_canopy, &
        & t_cf('fn_chl_cl', '', 'fraction of N in chlorophyll, pigments only'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('fn_et_cl', mem%fn_et_cl, &
        & hgrid, vgrid_canopy, &
        & t_cf('fn_et_cl', '', 'fraction of N in electron transport '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('fn_rub_cl', mem%fn_rub_cl, &
        & hgrid, vgrid_canopy, &
        & t_cf('fn_rub_cl', '', 'fraction of N in Rubisco '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('fn_pepc_cl', mem%fn_pepc_cl, &
        & hgrid, vgrid_canopy, &
        & t_cf('fn_pepc_cl', '', 'fraction of N in PEP carboylase, C4 plants only'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('fn_oth_cl', mem%fn_oth_cl, &
        & hgrid, vgrid_canopy, &
        & t_cf('fn_oth_cl', '', 'fraction of N not associated with photosynthesis '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('fleaf_sunlit_cl', mem%fleaf_sunlit_cl, &
        & hgrid, vgrid_canopy, &
        & t_cf('fleaf_sunlit_cl', '', 'fraction of layer which is sunlit'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('maint_respiration_pot', mem%maint_respiration_pot, &
        & hgrid, surface, &
        & t_cf('maint_respiration_pot', 'micro-mol CO2 m-2 s-1', 'potential maintenance respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('maint_respiration', mem%maint_respiration, &
        & hgrid, surface, &
        & t_cf('maint_respiration', 'micro-mol CO2 m-2 s-1', 'maintenance respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('maint_respiration_c13', mem%maint_respiration_c13, &
        & hgrid, surface, &
        & t_cf('maint_respiration_c13', 'micro-mol 13CO2 m-2 s-1', 'maintenance respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('maint_respiration_c14', mem%maint_respiration_c14, &
        & hgrid, surface, &
        & t_cf('maint_respiration_c14', 'micro-mol 14CO2 m-2 s-1', 'maintenance respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('growth_respiration', mem%growth_respiration, &
        & hgrid, surface, &
        & t_cf('growth_respiration', 'micro-mol CO2 m-2 s-1', 'growth respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('growth_respiration_c13', mem%growth_respiration_c13, &
        & hgrid, surface, &
        & t_cf('growth_respiration_c13', 'micro-mol 13CO2 m-2 s-1', 'growth respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('growth_respiration_c14', mem%growth_respiration_c14, &
        & hgrid, surface, &
        & t_cf('growth_respiration_c14', 'micro-mol 14CO2 m-2 s-1', 'growth respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('n_transform_respiration', mem%n_transform_respiration, &
        & hgrid, surface, &
        & t_cf('n_transform_respiration', 'micro-mol CO2 m-2 s-1', 'respiration associated with N transformation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('n_fixation_respiration', mem%n_fixation_respiration, &
        & hgrid, surface, &
        & t_cf('n_fixation_respiration', 'micro-mol CO2 m-2 s-1', 'respiration associated with N fixation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('n_processing_respiration', mem%n_processing_respiration, &
        & hgrid, surface, &
        & t_cf('n_processing_respiration', 'micro-mol CO2 m-2 s-1', 'respiration associated with N processing'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('n_processing_respiration_c13', mem%n_processing_respiration_c13, &
        & hgrid, surface, &
        & t_cf('n_processing_respiration_c13', 'micro-mol 13CO2 m-2 s-1', 'respiration associated with N processing'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('n_processing_respiration_c14', mem%n_processing_respiration_c14, &
        & hgrid, surface, &
        & t_cf('n_processing_respiration_c14', 'micro-mol 14CO2 m-2 s-1', 'respiration associated with N processing'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('npp', mem%npp, &
        & hgrid, surface, &
        & t_cf('npp', 'micro-mol CO2 m-2 s-1', 'net primary production'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('npp_c13', mem%npp_c13, &
        & hgrid, surface, &
        & t_cf('npp_c13', 'micro-mol 13CO2 m-2 s-1', 'net primary production'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('npp_c14', mem%npp_c14, &
        & hgrid, surface, &
        & t_cf('npp_c14', 'micro-mol 14CO2 m-2 s-1', 'net primary production'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('uptake_nh4', mem%uptake_nh4, &
        & hgrid, surface, &
        & t_cf('uptake_nh4', 'micro-mol N m-2 s-1', 'NH4 uptake from soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('uptake_nh4_n15', mem%uptake_nh4_n15, &
        & hgrid, surface, &
        & t_cf('uptake_nh4_n15', 'micro-mol N m-2 s-1', '15NH4 uptake from soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('uptake_no3', mem%uptake_no3, &
        & hgrid, surface, &
        & t_cf('uptake_no3', 'micro-mol N m-2 s-1', 'NO3 uptake from soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('uptake_no3_n15', mem%uptake_no3_n15, &
        & hgrid, surface, &
        & t_cf('uptake_no3_n15', 'micro-mol N m-2 s-1', '15NO3 uptake from soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('uptake_po4', mem%uptake_po4, &
        & hgrid, surface, &
        & t_cf('uptake_po4', 'micro-mol P m-2 s-1', 'PO4 uptake from soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('n_fixation', mem%n_fixation, &
        & hgrid, surface, &
        & t_cf('n_fixation', 'micro-mol N m-2 s-1', 'symbiontic N fixation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('n_fixation_n15', mem%n_fixation_n15, &
        & hgrid, surface, &
        & t_cf('n_fixation_n15', 'micro-mol 15N m-2 s-1', 'symbiontic N fixation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('recycling_leaf_n', mem%recycling_leaf_n, &
        & hgrid, surface, &
        & t_cf('recycling_leaf_n', 'mol m-2 timestep-1', 'net flux of N from leaf to labile (senescence and maintenance)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('recycling_leaf_n15', mem%recycling_leaf_n15, &
        & hgrid, surface, &
        & t_cf('recycling_leaf_n15', 'mol m-2 timestep-1', 'net flux of 15N from leaf to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('recycling_leaf_p', mem%recycling_leaf_p, &
        & hgrid, surface, &
        & t_cf('recycling_leaf_p', 'mol m-2 timestep-1', 'net flux of P from leaf to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('recycling_fine_root_n', mem%recycling_fine_root_n, &
        & hgrid, surface, &
        & t_cf('recycling_fine_root_n', 'mol m-2 timestep-1', 'net flux of N from fine roots to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('recycling_fine_root_n15', mem%recycling_fine_root_n15, &
        & hgrid, surface, &
        & t_cf('recycling_fine_root_n15', 'mol m-2 timestep-1', 'net flux of 15N from fine roots to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('recycling_fine_root_p', mem%recycling_fine_root_p, &
        & hgrid, surface, &
        & t_cf('recycling_fine_root_p', 'mol m-2 timestep-1', 'net flux of P from fine roots to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('recycling_heart_wood_n', mem%recycling_heart_wood_n, &
        & hgrid, surface, &
        & t_cf('recycling_heart_wood_n', 'mol m-2 timestep-1', 'net flux of N from sapwood to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('recycling_heart_wood_n15', mem%recycling_heart_wood_n15, &
        & hgrid, surface, &
        & t_cf('recycling_heart_wood_n15', 'mol m-2 timestep-1', 'net flux of 15N from sapwood to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('recycling_heart_wood_p', mem%recycling_heart_wood_p, &
        & hgrid, surface, &
        & t_cf('recycling_heart_wood_p', 'mol m-2 timestep-1', 'net flux of P from sapwood to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('net_growth', mem%net_growth, &
        & hgrid, surface, &
        & t_cf('net_growth', 'micro-mol CO2 m-2 s-1', 'npp minus turnover'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('fmaint_rate_root', mem%fmaint_rate_root, &
        & hgrid, surface, &
        & t_cf('fmaint_rate_root', '~micro-mol C mol-1 N s-1', 'current N-specific maintentance respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('unit_npp', mem%unit_npp, &
        & hgrid, surface, &
        & t_cf('unit_npp', '', 'incremental change in NPP. Needed for optimal allocation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('unit_transpiration', mem%unit_transpiration, &
        & hgrid, surface, &
        & t_cf('unit_transpiration', '', 'incremental change in transpiration. Needed for optimal allocation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('unit_uptake_n_pot', mem%unit_uptake_n_pot, &
        & hgrid, surface, &
        & t_cf('unit_uptake_n_pot', '', 'incremental change in potential N uptake per unit root carbon, N uptake calculation '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('unit_uptake_p_pot', mem%unit_uptake_p_pot, &
        & hgrid, surface, &
        & t_cf('unit_uptake_p_pot', '', 'incremental change in potential P uptake per unit root carbon, currently unused '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('unit_uptake_n_act', mem%unit_uptake_n_act, &
        & hgrid, surface, &
        & t_cf('unit_uptake_n_act', '', 'incremental change in actual N uptake per unit root carbon, optimal allocation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('unit_uptake_p_act', mem%unit_uptake_p_act, &
        & hgrid, surface, &
        & t_cf('unit_uptake_p_act', '', 'incremental change in actual P uptake per unit root carbon, optimal allocation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('delta_dens_ind', mem%delta_dens_ind, &
        & hgrid, surface, &
        & t_cf('delta_dens_ind', '# m-2 timestep-1', 'change in density of individuals'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('cost_n_uptake_root', mem%cost_n_uptake_root, &
        & hgrid, surface, &
        & t_cf('cost_n_uptake_root', 'mol C mol-1 N', 'actual cost for N uptake by roots'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('cost_n_fixation', mem%cost_n_fixation, &
        & hgrid, surface, &
        & t_cf('cost_n_fixation', 'mol C mol-1 N', 'actual cost for N uptake through fixation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('mortality_rate', mem%mortality_rate, &
        & hgrid, surface, &
        & t_cf('mortality_rate', 'timestep-1', 'background_mort_rate_tree / grass'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('k1_opt', mem%k1_opt, &
        & hgrid, surface, &
        & t_cf('k1_opt', '', 'coefficient for optimal leaf:root ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('k2_opt', mem%k2_opt, &
        & hgrid, surface, &
        & t_cf('k2_opt', '', 'coefficient for optimal leaf:root ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('leaf2sapwood_mass_ratio', mem%leaf2sapwood_mass_ratio, &
        & hgrid, surface, &
        & t_cf('leaf2sapwood_mass_ratio', '', 'leaf to shoot mass ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('leaf2root_mass_ratio', mem%leaf2root_mass_ratio, &
        & hgrid, surface, &
        & t_cf('leaf2root_mass_ratio', '', 'leaf to root mass ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('fleaf_sunlit_daytime_cl', mem%fleaf_sunlit_daytime_cl, &
        & hgrid, vgrid_canopy, &
        & t_cf('fleaf_sunlit_daytime_cl', '', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('t_air_daytime', mem%t_air_daytime, &
        & hgrid, surface, &
        & t_cf('t_air_daytime', '', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('press_srf_daytime', mem%press_srf_daytime, &
        & hgrid, surface, &
        & t_cf('press_srf_daytime', '', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('co2_mixing_ratio_daytime', mem%co2_mixing_ratio_daytime, &
        & hgrid, surface, &
        & t_cf('co2_mixing_ratio_daytime', '', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('ga_daytime', mem%ga_daytime, &
        & hgrid, surface, &
        & t_cf('ga_daytime', '', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('beta_sinklim_ps_daytime', mem%beta_sinklim_ps_daytime, &
        & hgrid, surface, &
        & t_cf('beta_sinklim_ps_daytime', 'unitless', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('t_jmax_opt_daytime', mem%t_jmax_opt_daytime, &
        & hgrid, surface, &
        & t_cf('t_jmax_opt_daytime', '', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('fleaf_sunlit_daytime_dacc_cl', mem%fleaf_sunlit_daytime_dacc_cl, &
        & hgrid, vgrid_canopy, &
        & t_cf('fleaf_sunlit_daytime_dacc_cl', '', 'value accumulated over current day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('t_air_daytime_dacc', mem%t_air_daytime_dacc, &
        & hgrid, surface, &
        & t_cf('t_air_daytime_dacc', '', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('press_srf_daytime_dacc', mem%press_srf_daytime_dacc, &
        & hgrid, surface, &
        & t_cf('press_srf_daytime_dacc', '', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('co2_mixing_ratio_daytime_dacc', mem%co2_mixing_ratio_daytime_dacc, &
        & hgrid, surface, &
        & t_cf('co2_mixing_ratio_daytime_dacc', '', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('ga_daytime_dacc', mem%ga_daytime_dacc, &
        & hgrid, surface, &
        & t_cf('ga_daytime_dacc', '', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('beta_sinklim_ps_daytime_dacc', mem%beta_sinklim_ps_daytime_dacc, &
        & hgrid, surface, &
        & t_cf('beta_sinklim_ps_daytime_dacc', '', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('t_jmax_opt_daytime_dacc', mem%t_jmax_opt_daytime_dacc, &
        & hgrid, surface, &
        & t_cf('t_jmax_opt_daytime_dacc', '', 'value accumulated over current day '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('gpp_tlabile_mavg', mem%gpp_tlabile_mavg, &
        & hgrid, surface, &
        & t_cf('gpp_tlabile_mavg', 'micro-mol m-2 s-1', 'gpp '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('maint_respiration_tlabile_mavg', mem%maint_respiration_tlabile_mavg, &
        & hgrid, surface, &
        & t_cf('maint_respiration_tlabile_mavg', 'micro-mol m-2 s-1', 'maintenance respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('growth_req_n_tlabile_mavg', mem%growth_req_n_tlabile_mavg, &
        & hgrid, surface, &
        & t_cf('growth_req_n_tlabile_mavg', '', 'moles N required for a unit of C growth under current allocation fractions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('growth_req_p_tlabile_mavg', mem%growth_req_p_tlabile_mavg, &
        & hgrid, surface, &
        & t_cf('growth_req_p_tlabile_mavg', '', 'moles P required for a unit of N growth under current allocation fractions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('t_air_tphen_mavg', mem%t_air_tphen_mavg, &
        & hgrid, surface, &
        & t_cf('t_air_tphen_mavg', '', 'air temp mavg for phenology calc'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('npp_tuptake_mavg', mem%npp_tuptake_mavg, &
        & hgrid, surface, &
        & t_cf('npp_tuptake_mavg', 'micro-mol m-2 s-1', 'npp of the last week'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('demand_uptake_n_tuptake_mavg', mem%demand_uptake_n_tuptake_mavg, &
        & hgrid, surface, &
        & t_cf('demand_uptake_n_tuptake_mavg', '', 'scalar of Vmax for root N uptake '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('demand_uptake_p_tuptake_mavg', mem%demand_uptake_p_tuptake_mavg, &
        & hgrid, surface, &
        & t_cf('demand_uptake_p_tuptake_mavg', '', 'scalar of Vmax for root P uptake'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('growth_req_n_tuptake_mavg', mem%growth_req_n_tuptake_mavg, &
        & hgrid, surface, &
        & t_cf('growth_req_n_tuptake_mavg', '', 'moles N required for a unit of C growth under current allocation fractions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('growth_req_p_tuptake_mavg', mem%growth_req_p_tuptake_mavg, &
        & hgrid, surface, &
        & t_cf('growth_req_p_tuptake_mavg', '', 'moles P required for a unit of N growth under current allocation fractions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('cost_n_fixation_talloc_mavg', mem%cost_n_fixation_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('cost_n_fixation_talloc_mavg', 'mol C mol-1 N', 'moving average of N fixation cost'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('exudation_c_tmyc_mavg', mem%exudation_c_tmyc_mavg, &
        & hgrid, surface, &
        & t_cf('exudation_c_tmyc_mavg', 'mol C', 'moving average of C allocation to mycorrhizae'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fleaf_sunlit_tfrac_mavg_cl', mem%fleaf_sunlit_tfrac_mavg_cl, &
        & hgrid, vgrid_canopy, &
        & t_cf('fleaf_sunlit_tfrac_mavg_cl', '', 'fraction of sunlit leaves'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('t_air_tfrac_mavg', mem%t_air_tfrac_mavg, &
        & hgrid, surface, &
        & t_cf('t_air_tfrac_mavg', '', 'air temp'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('press_srf_tfrac_mavg', mem%press_srf_tfrac_mavg, &
        & hgrid, surface, &
        & t_cf('press_srf_tfrac_mavg', '', 'air pressure at surface'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('co2_mixing_ratio_tfrac_mavg', mem%co2_mixing_ratio_tfrac_mavg, &
        & hgrid, surface, &
        & t_cf('co2_mixing_ratio_tfrac_mavg', '', 'co2 mixing ratio of the air'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('ga_tfrac_mavg', mem%ga_tfrac_mavg, &
        & hgrid, surface, &
        & t_cf('ga_tfrac_mavg', '', 'aerodynamic conductance'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('beta_sinklim_ps_tfrac_mavg', mem%beta_sinklim_ps_tfrac_mavg, &
        & hgrid, surface, &
        & t_cf('beta_sinklim_ps_tfrac_mavg', '', 'scaling factor to account for soil moisture constraints on photosynthesis '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('fleaf_sunlit_tcnl_mavg_cl', mem%fleaf_sunlit_tcnl_mavg_cl, &
        & hgrid, vgrid_canopy, &
        & t_cf('fleaf_sunlit_tcnl_mavg_cl', '', 'fraction of sunlit leaves'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('t_air_tcnl_mavg', mem%t_air_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('t_air_tcnl_mavg', '', 'air temp'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('press_srf_tcnl_mavg', mem%press_srf_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('press_srf_tcnl_mavg', '', 'air pressure at surface'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('co2_mixing_ratio_tcnl_mavg', mem%co2_mixing_ratio_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('co2_mixing_ratio_tcnl_mavg', '', 'co2 mixing ratio of the air'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('ga_tcnl_mavg', mem%ga_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('ga_tcnl_mavg', '', 'aerodynamic conductance'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('beta_sinklim_ps_tcnl_mavg', mem%beta_sinklim_ps_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('beta_sinklim_ps_tcnl_mavg', '', 'scaling factor to account for soil moisture constraints on conductance'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('uptake_n_tcnl_mavg', mem%uptake_n_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('uptake_n_tcnl_mavg', '', 'root N uptake'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_cn_tcnl_mavg', mem%growth_cn_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('growth_cn_tcnl_mavg', '', 'new growth CN ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_np_tcnl_mavg', mem%growth_np_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('growth_np_tcnl_mavg', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('npp_tcnl_mavg', mem%npp_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('npp_tcnl_mavg', '', 'NPP'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fmaint_rate_tcnl_mavg', mem%fmaint_rate_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('fmaint_rate_tcnl_mavg', '~micro-mol C mol-1 N s-1', 'N-specific maintentance respiration for optmal leaf N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('labile_nitrogen_tcnl_mavg', mem%labile_nitrogen_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('labile_nitrogen_tcnl_mavg', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('labile_carbon_tcnl_mavg', mem%labile_carbon_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('labile_carbon_tcnl_mavg', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('labile_phosphorus_tcnl_mavg', mem%labile_phosphorus_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('labile_phosphorus_tcnl_mavg', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
  
      CALL mem%Add_var('npp_talloc_mavg', mem%npp_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('npp_talloc_mavg', 'micro-mol m-2 s-1', 'NPP'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('n_fixation_talloc_mavg', mem%n_fixation_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('n_fixation_talloc_mavg', '', 'N fixation averaged at the allocation timescale'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('transpiration_talloc_mavg', mem%transpiration_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('transpiration_talloc_mavg', '', 'average transpiration at allocation time-scale'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('unit_npp_talloc_mavg', mem%unit_npp_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('unit_npp_talloc_mavg', 'micro-mol C mol-1', 'rate of change in NPP with  increment in leaf mass'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('unit_transpiration_talloc_mavg', mem%unit_transpiration_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('unit_transpiration_talloc_mavg', '', 'rate of tranpiration change per unit LAI change?'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('unit_uptake_n_talloc_mavg', mem%unit_uptake_n_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('unit_uptake_n_talloc_mavg', 'micro-mol N mol-1', 'actual rate of N uptkae per unit root'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('unit_uptake_p_talloc_mavg', mem%unit_uptake_p_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('unit_uptake_p_talloc_mavg', 'micro-mol P mol-1', 'actual rate of P uptake per unit root'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('dphi_talloc_mavg', mem%dphi_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('dphi_talloc_mavg', 'MPa', 'leaf-soil water potential gradient'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('growth_cn_talloc_mavg', mem%growth_cn_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('growth_cn_talloc_mavg', '', 'plant growth cn ratio, equivalent to 1/growth_req_n'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('growth_cp_talloc_mavg', mem%growth_cp_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('growth_cp_talloc_mavg', '', 'plant growth cp ratio, equivalent to 1 / (growthreq_n * growth_req_p)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_np_talloc_mavg', mem%growth_np_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('growth_np_talloc_mavg', '', 'plant growth hp ratio, equivalent to 1/(growth_req_p)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('labile_carbon_talloc_mavg', mem%labile_carbon_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('labile_carbon_talloc_mavg', '', 'moving average of labile carbon pool, allometry calculations'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('labile_nitrogen_talloc_mavg', mem%labile_nitrogen_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('labile_nitrogen_talloc_mavg', '', 'moving average of labile nitrogen pool, allometry calculations'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('labile_phosphorus_talloc_mavg', mem%labile_phosphorus_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('labile_phosphorus_talloc_mavg', '', 'moving average of labile phosphorus pool, allometry calculations'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('reserve_carbon_talloc_mavg', mem%reserve_carbon_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('reserve_carbon_talloc_mavg', '', 'moving average of reserve carbon, optimal allocation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('reserve_nitrogen_talloc_mavg', mem%reserve_nitrogen_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('reserve_nitrogen_talloc_mavg', '', 'moving average of reserve nitrogen, optimal allocation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('reserve_phosphorus_talloc_mavg', mem%reserve_phosphorus_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('reserve_phosphorus_talloc_mavg', '', 'moving average of reserve phosphorus, optimal allocation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('w_root_lim_talloc_mavg', mem%w_root_lim_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('w_root_lim_talloc_mavg', '', 'moving average of ...'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
      
      CALL mem%Add_var('leaf2root_troot_mavg', mem%leaf2root_troot_mavg, &
        & hgrid, surface, &
        & t_cf('leaf2root_troot_mavg', '', 'leaf to root mass ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('unit_npp_troot_mavg', mem%unit_npp_troot_mavg, &
        & hgrid, surface, &
        & t_cf('unit_npp_troot_mavg', 'micro-mol C mol-1', 'rate of change in NPP with increment in leaf mass'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('unit_uptake_n_pot_troot_mavg', mem%unit_uptake_n_pot_troot_mavg, &
        & hgrid, surface, &
        & t_cf('unit_uptake_n_pot_troot_mavg', 'micro-mol N mol-1', 'rate of N uptkae per unit root'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('unit_uptake_p_pot_troot_mavg', mem%unit_uptake_p_pot_troot_mavg, &
        & hgrid, surface, &
        & t_cf('unit_uptake_p_pot_troot_mavg', 'micro-mol P mol-1', 'rate of P uptake per unit root'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('fmaint_rate_troot_mavg', mem%fmaint_rate_troot_mavg, &
        & hgrid, surface, &
        & t_cf('fmaint_rate_troot_mavg', '~micro-mol C mol-1 N s-1', 'moving average of N-specific maintentance respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('an_boc_tvegdyn_mavg', mem%an_boc_tvegdyn_mavg, &
        & hgrid, surface, &
        & t_cf('an_boc_tvegdyn_mavg', '', 'net photosynthesis at bottom of canopy, lowest canopy layer existing'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('net_growth_tvegdyn_mavg', mem%net_growth_tvegdyn_mavg, &
        & hgrid, surface, &
        & t_cf('net_growth_tvegdyn_mavg', 'micro-mol C m-2 s-1', 'net growth, npp minus turnover'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('lai_tvegdyn_mavg', mem%lai_tvegdyn_mavg, &
        & hgrid, surface, &
        & t_cf('lai_tvegdyn_mavg', 'm2 m-2', 'lai at vegdynamics timescale'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('t_air_tacclim_mavg', mem%t_air_tacclim_mavg, &
        & hgrid, surface, &
        & t_cf('t_air_tacclim_mavg', '', 'average air temperature for respiration acclimation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('t_soil_root_tacclim_mavg', mem%t_soil_root_tacclim_mavg, &
        & hgrid, surface, &
        & t_cf('t_soil_root_tacclim_mavg', '', 'average root temperature for respiration acclimation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('t_air_week_mavg', mem%t_air_week_mavg, &
        & hgrid, surface, &
        & t_cf('t_air_week_mavg', '', 'avrg air temperature across 7 days'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('t_air_month_mavg', mem%t_air_month_mavg, &
        & hgrid, surface, &
        & t_cf('t_air_month_mavg', '', 'avrg air temperature across 30 days'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  
      CALL mem%Add_var('t_jmax_opt_mavg', mem%t_jmax_opt_mavg, &
        & hgrid, surface, &
        & t_cf('t_jmax_opt_mavg', '', 'optimim temperature for electron transport of photosynthesis'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

#ifdef __QUINCY_STANDALONE__
#else
    ENDIF  ! IF One_of(LAND_TYPE  OR VEG_TYPE, lct_ids)
#endif

  END SUBROUTINE Init_veg_memory

#ifdef __QUINCY_STANDALONE__
#else
  ! ======================================================================================================= !
  !>
  !> return name of element variable of bgc_material by its ENUM ID
  !> helper function for t_veg_bgcm_with_elements
  !> "bgcm_object%Get_element_name_by_id" => veg_bgcm_with_elements_get_element_name_by_id
  !>
  FUNCTION veg_bgcm_with_elements_get_element_name_by_id(this, id) RESULT(name)
    CLASS(t_veg_bgcm_with_elements), INTENT(in)  :: this
    INTEGER,                         INTENT(in)  :: id
    CHARACTER(LEN=:),                ALLOCATABLE :: name

    ! return name of element variable only if this bgc_material object contains elements
    IF (this%contains_elements) THEN
      name = get_name_for_element_id(id)
    ELSE
      name='no_elements'
    END IF
  END FUNCTION veg_bgcm_with_elements_get_element_name_by_id
  
  ! ======================================================================================================= !
  !>
  !> return name of element variable of bgc_material by its ENUM ID
  !> helper function for t_veg_bgcm_with_components
  !> "bgcm_object%Get_element_name_by_id" => veg_bgcm_with_components_get_element_name_by_id
  !>
  FUNCTION veg_bgcm_with_components_get_element_name_by_id(this, id) RESULT(name)
    CLASS(t_veg_bgcm_with_components), INTENT(in)  :: this
    INTEGER,                           INTENT(in)  :: id
    CHARACTER(LEN=:),                  ALLOCATABLE :: name

    ! return name of element variable only if this bgc_material object contains elements
    IF (this%contains_elements) THEN
      name = get_name_for_element_id(id)
    ELSE
      name='no_elements'
    END IF
  END FUNCTION veg_bgcm_with_components_get_element_name_by_id
  
  ! ======================================================================================================= !
  !>
  !> return name of element variable of bgc_material by its ENUM ID
  !> helper function for t_veg_bgcm_with_components
  !> "bgcm_object%Get_element_name_by_id" => veg_bgcm_main_get_element_name_by_id
  !>
  FUNCTION veg_bgcm_main_get_element_name_by_id(this, id) RESULT(name)
    CLASS(t_veg_bgcm_main), INTENT(in)  :: this
    INTEGER,                INTENT(in)  :: id
    CHARACTER(LEN=:),       ALLOCATABLE :: name

    ! return name of element variable only if this bgc_material object contains elements
    IF (this%contains_elements) THEN
      name = get_name_for_element_id(id)
    ELSE
      name='no_elements'
    END IF
  END FUNCTION veg_bgcm_main_get_element_name_by_id

  ! ======================================================================================================= !
  !>
  !> initializing the t_veg_bgcm_with_elements object
  !> "bgcm_object%Init()"
  !> by default it does contain element variables !
  !>
  SUBROUTINE veg_bgcm_with_elements_init(this, name, shortname, elements_index_map, &
    &                                    l_elements, element_list, element_unit)
    USE mo_lnd_bgcm_class,      ONLY: ELEM_C_ID, ELEM_N_ID, ELEM_P_ID, ELEM_C13_ID, ELEM_C14_ID, ELEM_N15_ID, &
      &                               LAST_ELEM_ID
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_veg_bgcm_with_elements),    INTENT(inout) :: this                   ! object of TYPE t_sb_elements_pool
    CHARACTER(LEN=*),                   INTENT(in)    :: name                   ! pool name       (No spaces !)
    CHARACTER(LEN=*),                   INTENT(in)    :: shortname              ! pool short name (No spaces !)
    INTEGER,                            INTENT(in)    :: elements_index_map(:)  ! elements ID -> IND
    LOGICAL,          OPTIONAL,         INTENT(in)    :: l_elements             ! IF pool has elements 
    INTEGER,          OPTIONAL,         INTENT(in)    :: element_list(:) ! @TODO remove | list of elements; used if only selected elements should be added
    CHARACTER(LEN=*), OPTIONAL,         INTENT(in)    :: element_unit           ! unit of elements
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':veg_bgcm_with_elements_init'

    this%name              = name
    this%shortname         = shortname
    this%contains_elements = l_elements
    this%contains_bgcm     = .FALSE.

    IF (.NOT. PRESENT(element_unit)) CALL finish(routine, 'Missing unit for elements')
    this%element_unit = element_unit

    ! add element variables
    ! carbon
    IF (elements_index_map(ELEM_C_ID) > 0) THEN
      CALL this%Add_element('carbon', 'C', ELEM_C_ID, dim=2)
      SELECT TYPE (selector => this%element_list(elements_index_map(ELEM_C_ID))%p)
      CLASS IS (t_jsb_var_real2d)
        this%carbon => selector
      END SELECT
    END IF
    ! nitrogen
    IF (elements_index_map(ELEM_N_ID) > 0) THEN
      CALL this%Add_element('nitrogen', 'N', ELEM_N_ID, dim=2)
      SELECT TYPE (selector => this%element_list(elements_index_map(ELEM_N_ID))%p)
      CLASS IS (t_jsb_var_real2d)
        this%nitrogen => selector
      END SELECT
    END IF
    ! phosphorus
    IF (elements_index_map(ELEM_P_ID) > 0) THEN
      CALL this%Add_element('phosphorus', 'P', ELEM_P_ID, dim=2)
      SELECT TYPE (selector => this%element_list(elements_index_map(ELEM_P_ID))%p)
      CLASS IS (t_jsb_var_real2d)
        this%phosphorus => selector
      END SELECT
    END IF
    ! carbon13
    IF (elements_index_map(ELEM_C13_ID) > 0) THEN
      CALL this%Add_element('carbon13', 'C13', ELEM_C13_ID, dim=2)
      SELECT TYPE (selector => this%element_list(elements_index_map(ELEM_C13_ID))%p)
      CLASS IS (t_jsb_var_real2d)
        this%carbon13 => selector
      END SELECT
    END IF
    ! carbon14
    IF (elements_index_map(ELEM_C14_ID) > 0) THEN
      CALL this%Add_element('carbon14', 'C14', ELEM_C14_ID, dim=2)
      SELECT TYPE (selector => this%element_list(elements_index_map(ELEM_C14_ID))%p)
      CLASS IS (t_jsb_var_real2d)
        this%carbon14 => selector
      END SELECT
    END IF
    ! nitrogen15
    IF (elements_index_map(ELEM_N15_ID) > 0) THEN
      CALL this%Add_element('nitrogen15', 'N15', ELEM_N15_ID, dim=2)
      SELECT TYPE (selector => this%element_list(elements_index_map(ELEM_N15_ID))%p)
      CLASS IS (t_jsb_var_real2d)
        this%nitrogen15 => selector
      END SELECT
    END IF
  END SUBROUTINE veg_bgcm_with_elements_init

  ! ======================================================================================================= !
  !>
  !> initializing the t_veg_bgcm_with_components object
  !> "bgcm_object%Init()"
  !> by default it does not contain element variables (only t_veg_bgcm_with_elements do)
  !>
  SUBROUTINE veg_bgcm_with_components_init(this, name, shortname, elements_index_map, &
    &                                      l_elements, element_list, element_unit)
    CLASS(t_veg_bgcm_with_components),    INTENT(inout) :: this                   ! object of TYPE t_sb_components_pool
    CHARACTER(LEN=*),                     INTENT(in)    :: name                   ! pool name       (No spaces !)
    CHARACTER(LEN=*),                     INTENT(in)    :: shortname              ! pool short name (No spaces !)
    INTEGER,                              INTENT(in)    :: elements_index_map(:)  ! elements ID -> IND
    LOGICAL,          OPTIONAL,           INTENT(in)    :: l_elements             ! IF pool has elements 
    INTEGER,          OPTIONAL,           INTENT(in)    :: element_list(:)        ! list of elements; used if only selected elements should be added
    CHARACTER(LEN=*), OPTIONAL,           INTENT(in)    :: element_unit           ! unit of elements
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':veg_bgcm_with_components_init'

    this%name              = name
    this%shortname         = shortname
    this%contains_elements = l_elements
    this%contains_bgcm     = .TRUE.

    IF (.NOT. PRESENT(element_unit)) CALL finish(routine, 'Missing unit for elements')
    this%element_unit = element_unit
  END SUBROUTINE veg_bgcm_with_components_init
  
  ! ======================================================================================================= !
  !>
  !> initializing the t_veg_bgcm_main object
  !> "bgcm_object%Init()"
  !> by default it does not contain element variables (only t_veg_bgcm_with_elements do)
  !>
  SUBROUTINE veg_bgcm_main_init(this, name, shortname, elements_index_map, &
    &                           l_elements, element_list, element_unit)
    CLASS(t_veg_bgcm_main),         INTENT(inout) :: this                   ! object of TYPE t_sb_components_pool
    CHARACTER(LEN=*),               INTENT(in)    :: name                   ! pool name       (No spaces !)
    CHARACTER(LEN=*),               INTENT(in)    :: shortname              ! pool short name (No spaces !)
    INTEGER,                        INTENT(in)    :: elements_index_map(:)  ! elements ID -> IND
    LOGICAL,          OPTIONAL,     INTENT(in)    :: l_elements             ! IF pool has elements 
    INTEGER,          OPTIONAL,     INTENT(in)    :: element_list(:)        ! list of elements; used if only selected elements should be added
    CHARACTER(LEN=*), OPTIONAL,     INTENT(in)    :: element_unit           ! unit of elements
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':veg_bgcm_main_init'
  
    this%name              = name
    this%shortname         = shortname
    this%contains_elements = l_elements
    this%contains_bgcm     = .TRUE.
  END SUBROUTINE veg_bgcm_main_init

  ! ======================================================================================================= !
  !>
  !> create and init an object of the VEG_ components bgcm "t_veg_bgcm_with_components"
  !> used for VEG_ memory init
  !>   init and add all its objects of the "t_jsb_pool" with "pool_object%Init()"
  !>   "object => Create_veg_bgc_material_components()"
  !>
  FUNCTION Create_veg_bgc_material_components(name, shortname, unit, lct_ids, elements_index_map, pool_prefix, pool_suffix) &
    &      RESULT(veg_bgcm_components)
    CHARACTER(LEN=*),                     INTENT(in)  :: name
    CHARACTER(LEN=*),                     INTENT(in)  :: shortname
    CHARACTER(LEN=*),                     INTENT(in)  :: unit
    INTEGER,                              INTENT(in)  :: lct_ids(:)
    INTEGER,                              INTENT(in)  :: elements_index_map(:)  ! elements ID -> IND
    CHARACTER(LEN=*), OPTIONAL,           INTENT(in)  :: pool_prefix
    CHARACTER(LEN=*), OPTIONAL,           INTENT(in)  :: pool_suffix
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_veg_bgcm_with_components),  POINTER     :: veg_bgcm_components
    TYPE(t_veg_bgcm_with_elements),    POINTER     :: veg_bgcm_elements
    CHARACTER(LEN=:),                     ALLOCATABLE :: prefix_loc
    CHARACTER(LEN=:),                     ALLOCATABLE :: suffix_loc
    ! ----------------------------------------------------------------------------------------------------- !
    prefix_loc = ''
    suffix_loc = ''
    CALL assign_if_present_allocatable(prefix_loc, pool_prefix)
    CALL assign_if_present_allocatable(suffix_loc, pool_suffix)
    IF (prefix_loc /= '') prefix_loc = prefix_loc // '_'
    IF (suffix_loc /= '') suffix_loc = '_' // suffix_loc

    ! bgcm_components init
    ALLOCATE(veg_bgcm_components)
    CALL veg_bgcm_components%Init(name, shortname, elements_index_map(:), l_elements=.FALSE., element_unit=unit)

    ! bgcm_elements allocate & init & add to bgcm_components
    !   "shortname" is used for model output
    !   "element_list" is optional and can be used to select specific elements (not needed when all elements may be added)
    !! @TODO: maybe extract as a routine given two strings and returning the pointer 
    !!    routine should start with IF(ASSOCIATED) NULLIFY(sb_bgcm_elements)
    !!    i.e. only left: function call with strings and subsequent assignment
    ! leaf
    ALLOCATE(veg_bgcm_elements)
    CALL veg_bgcm_elements%Init(name='leaf_bgcm'//suffix_loc, shortname='leaf', &
      &                        elements_index_map=elements_index_map(:), &
      &                        l_elements=.TRUE., element_unit=unit)
    CALL veg_bgcm_components%Add_pool(veg_bgcm_elements)
    NULLIFY(veg_bgcm_elements)
    ! fine_root
    ALLOCATE(veg_bgcm_elements)
    CALL veg_bgcm_elements%Init(name='fine_root_bgcm'//suffix_loc, shortname='fine_root', &
      &                        elements_index_map=elements_index_map(:), &
      &                        l_elements=.TRUE., element_unit=unit)
    CALL veg_bgcm_components%Add_pool(veg_bgcm_elements)
    NULLIFY(veg_bgcm_elements)
    ! coarse_root
    ALLOCATE(veg_bgcm_elements)
    CALL veg_bgcm_elements%Init(name='coarse_root_bgcm'//suffix_loc, shortname='coarse_root', &
      &                        elements_index_map=elements_index_map(:), &
      &                        l_elements=.TRUE., element_unit=unit)
    CALL veg_bgcm_components%Add_pool(veg_bgcm_elements)
    NULLIFY(veg_bgcm_elements)
    ! sap_wood
    ALLOCATE(veg_bgcm_elements)
    CALL veg_bgcm_elements%Init(name='sap_wood_bgcm'//suffix_loc, shortname='sap_wood', &
      &                        elements_index_map=elements_index_map(:), &
      &                        l_elements=.TRUE., element_unit=unit)
    CALL veg_bgcm_components%Add_pool(veg_bgcm_elements)
    NULLIFY(veg_bgcm_elements)
    ! heart_wood
    ALLOCATE(veg_bgcm_elements)
    CALL veg_bgcm_elements%Init(name='heart_wood_bgcm'//suffix_loc, shortname='heart_wood', &
      &                        elements_index_map=elements_index_map(:), &
      &                        l_elements=.TRUE., element_unit=unit)
    CALL veg_bgcm_components%Add_pool(veg_bgcm_elements)
    NULLIFY(veg_bgcm_elements)
    ! labile
    ALLOCATE(veg_bgcm_elements)
    CALL veg_bgcm_elements%Init(name='labile_bgcm'//suffix_loc, shortname='labile', &
      &                        elements_index_map=elements_index_map(:), &
      &                        l_elements=.TRUE., element_unit=unit)
    CALL veg_bgcm_components%Add_pool(veg_bgcm_elements)
    NULLIFY(veg_bgcm_elements)
    ! reserve
    ALLOCATE(veg_bgcm_elements)
    CALL veg_bgcm_elements%Init(name='reserve_bgcm'//suffix_loc, shortname='reserve', &
      &                        elements_index_map=elements_index_map(:), &
      &                        l_elements=.TRUE., element_unit=unit)
    CALL veg_bgcm_components%Add_pool(veg_bgcm_elements)
    NULLIFY(veg_bgcm_elements)
    ! fruit
    ALLOCATE(veg_bgcm_elements)
    CALL veg_bgcm_elements%Init(name='fruit_bgcm'//suffix_loc, shortname='fruit', &
      &                        elements_index_map=elements_index_map(:), &
      &                        l_elements=.TRUE., element_unit=unit)
    CALL veg_bgcm_components%Add_pool(veg_bgcm_elements)
    NULLIFY(veg_bgcm_elements)
    ! seed_bed
    ALLOCATE(veg_bgcm_elements)
    CALL veg_bgcm_elements%Init(name='seed_bed_bgcm'//suffix_loc, shortname='seed_bed', &
      &                        elements_index_map=elements_index_map(:), &
      &                        l_elements=.TRUE., element_unit=unit)
    CALL veg_bgcm_components%Add_pool(veg_bgcm_elements)
    NULLIFY(veg_bgcm_elements)
  END FUNCTION Create_veg_bgc_material_components
#endif  

#endif
END MODULE mo_veg_memory_class
