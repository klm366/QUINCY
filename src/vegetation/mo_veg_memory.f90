!>
!! @brief Contains the routines for plant memory operations 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_veg_memory

  USE mo_kind,              ONLY: wp, eps8
  USE mo_base_constants,    ONLY: ncanop, dt, one_day
  USE mo_shared_memory,     ONLY: t_element

  IMPLICIT NONE


  !-----------------------------------------------------------------------------------------------------
  !> vegetation pools storing elements in tissues (e.g. leaf, labile, reserve)
  !! 
  !! this type is used for vegetation pools, and fluxes in & out of pools;
  !! pools: there is one pool per plant and tissue
  !! flux: there are two fluxes per pool of a tissue, growth & litterfall (in & out of pool)
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_veg_pool 
    TYPE(t_element) :: leaf           !< foliage
    TYPE(t_element) :: fine_root      !< fine roots (< 2mm diameter)
    TYPE(t_element) :: coarse_root    !< coarse roots (> 2mm diameter)
    TYPE(t_element) :: sap_wood       !< living wood for water and phloem transport
    TYPE(t_element) :: heart_wood     !< dead standing wood
    TYPE(t_element) :: labile         !< meristem carbon
    TYPE(t_element) :: reserve        !< storage
    TYPE(t_element) :: fruit          !< live propagation pool
    TYPE(t_element) :: seed_bed       !< fruit mass on ground
  CONTAINS
    PRIVATE
    PROCEDURE :: add_veg_pools_type_type
    PROCEDURE :: add_veg_pools_type_scalar
    PROCEDURE, PASS(this2) :: add_veg_pools_scalar_type
    GENERIC, PUBLIC :: OPERATOR(+) => add_veg_pools_type_type, add_veg_pools_type_scalar, add_veg_pools_scalar_type
    PROCEDURE :: subtract_veg_pools_type_type
    PROCEDURE :: subtract_veg_pools_type_scalar
    PROCEDURE, PASS(this2) :: subtract_veg_pools_scalar_type
    GENERIC, PUBLIC :: OPERATOR(-) => subtract_veg_pools_type_type, subtract_veg_pools_type_scalar, subtract_veg_pools_scalar_type
    PROCEDURE :: multiply_veg_pools_type_type
    PROCEDURE :: multiply_veg_pools_type_scalar
    PROCEDURE, PASS(this2) :: multiply_veg_pools_scalar_type
    GENERIC, PUBLIC :: OPERATOR(*) => multiply_veg_pools_type_type, multiply_veg_pools_type_scalar, multiply_veg_pools_scalar_type
    PROCEDURE :: divide_veg_pools_type_type
    PROCEDURE :: divide_veg_pools_type_scalar
    PROCEDURE, PASS(this2) :: divide_veg_pools_scalar_type
    GENERIC, PUBLIC :: OPERATOR(/) => divide_veg_pools_type_type, divide_veg_pools_type_scalar, divide_veg_pools_scalar_type
    PROCEDURE, PUBLIC :: init => init_veg_pool
  END TYPE t_slm_veg_pool

!  !-----------------------------------------------------------------------------------------------------
!  !> all fluxes of vegetation (in & out of pools)
!  !! 
!  !! uses the TYPE t_slm_veg_pool to define fluxes: growth & litterfall \n
!  !! plant\%flux\%growth & plant\%flux\%litterfall
!  !! 
!  !! includes additionally: TYPE(t_element) :: reserve_use \n
!  !! plant\%flux\%reserve_use
!  !-----------------------------------------------------------------------------------------------------
!  TYPE t_slm_veg_flux 
!     TYPE(t_element)      :: establishment !< flux between seed bed and reserve (mol m-2 timestep-1)
!     TYPE(t_element)      :: reserve_use  !< net flux between reserve and labile pool (mol m-2 timestep-1)
!                                          !!  (+ve: flux to labile; -ve flux from labile) 
!     TYPE(t_slm_veg_pool) :: growth       !< add to      plant%pool due to growth     per timestep (mol m-2 timestep-1)
!     TYPE(t_slm_veg_pool) :: litterfall   !< remove from plant%pool due to litterfall per timestep (mol m-2 timestep-1)
!     REAL(wp) :: gpp                      !< gross photosynthesis (µmol CO2 m-2 s-1)
!     REAL(wp) :: an                       !< net photosynthesis (µmol CO2 m-2 s-1)
!     REAL(wp) :: gpp_c13                  !< gross photosynthesis (µmol 13CO2 m-2 s-1)
!     REAL(wp) :: gpp_c14                  !< gross photosynthesis (µmol 14CO2 m-2 s-1)
!     REAL(wp) :: maint_respiration_leaf   !< foliar maintenance respiration (µmol CO2 m-2 s-1)
!     REAL(wp) :: maint_respiration_pot    !< potential maintenance respiration (µmol CO2 m-2 s-1)
!     REAL(wp) :: maint_respiration        !< maintenance respiration (µmol CO2 m-2 s-1)
!     REAL(wp) :: maint_respiration_c13    !< maintenance respiration (µmol 13CO2 m-2 s-1)
!     REAL(wp) :: maint_respiration_c14    !< maintenance respiration (µmol 14CO2 m-2 s-1)
!     REAL(wp) :: growth_respiration       !< growth respiration (µmol CO2 m-2 s-1)
!     REAL(wp) :: growth_respiration_c13   !< growth respiration (µmol 13CO2 m-2 s-1)
!     REAL(wp) :: growth_respiration_c14   !< growth respiration (µmol 14CO2 m-2 s-1)
!     REAL(wp) :: npp                      !< net primary production (µmol CO2 m-2 s-1)
!     REAL(wp) :: npp_c13                  !< net primary production (µmol 13CO2 m-2 s-1)
!     REAL(wp) :: npp_c14                  !< net primary production (µmol 14CO2 m-2 s-1)
!     REAL(wp) :: uptake_n                 !< N uptake from soil ( is that correct? )
!     REAL(wp) :: uptake_n15               !< 15N uptake from soil ( is that correct? ) 
!     REAL(wp) :: uptake_p                 !< P uptake from soil ( is that correct? )
!     REAL(wp) :: recycling_leaf_n         !< net flux of N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
!     REAL(wp) :: recycling_leaf_n15       !< net flux of 15N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1) 
!     REAL(wp) :: recycling_leaf_p         !< net flux of P from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
!     REAL(wp) :: recycling_fine_root_n    !< net flux of N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
!     REAL(wp) :: recycling_fine_root_n15  !< net flux of 15N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1) 
!     REAL(wp) :: recycling_fine_root_p    !< net flux of P from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
!     REAL(wp) :: recycling_sap_wood_n     !< net flux of N from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1)
!     REAL(wp) :: recycling_sap_wood_n15   !< net flux of 15N from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1) 
!     REAL(wp) :: recycling_sap_wood_p     !< net flux of P from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1)
!     REAL(wp) :: transpiration            !< transpiration in kg m-2 s-1
!     REAL(wp) :: ci                       !< leaf internal co2 concentration
!     REAL(wp) :: unit_npp                 !< incremental change in NPP. Needed for optimal allocation
!     REAL(wp) :: unit_uptake_n            !< incremental change in N-based growth, needed for optimal allocation
!     REAL(wp) :: unit_uptake_p            !< incremental change in P-based growth, needed for optimal allocation
!     REAL(wp) :: delta_dens_ind           !< change in density of individuals (#/m2/timestep)
!  CONTAINS
!    PRIVATE
!    PROCEDURE, PUBLIC :: set_zero => set_zero_veg_flux
!  END TYPE t_slm_veg_flux


  !-----------------------------------------------------------------------------------------------------
  !> all variables of vegetation per tile
  !! 
  !! including all pools and all fluxes
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_plant 
    TYPE(t_slm_veg_pool)      :: pool             !< the actual vegetation pools
    !TYPE(t_slm_veg_flux)      :: flux             !< contains growth and litter fall fluxes per timestep for each pool
    LOGICAL  :: growing_season                    !< logical to identify whether the plant is in the growing season or not
    REAL(wp) :: height                            !< plant height in m
    REAL(wp) :: diameter                          !< plant diameter in m
    REAL(wp) :: dens_ind                          !< density of individuals (#/m2) 
    REAL(wp) :: lai                               !< actual leaf area index
    REAL(wp) :: target_cn_leaf                    !< target foliar C:N ratio
    REAL(wp) :: target_np_leaf                    !< target foliar N:P ratio 
    REAL(wp) :: mean_leaf_age                     !< average foliage age in days 
    REAL(wp) :: gc                                !< some kind of stomatal conductance ?
    REAL(wp) :: f_n_demand                        !< relative plant demand for N; see update_vegsoil_bgc
    REAL(wp) :: f_p_demand                        !< relative plant demand for P; see update_vegsoil_bgc
    REAL(wp) :: beta_air                          !< scaling factor to account for air humidity effects on conductance
    REAL(wp) :: beta_soil_ps                      !< scaling factor to account for soil moisture constraints on conductance
    REAL(wp) :: dphi                              !< plant-soil water potential gradient (MPa)
    REAL(wp) :: beta_soil_gs                      !< scaling factor to account for soil moisture constraints on conductance
    REAL(wp) :: growth_req_n                      !< moles N required for a unit of C growth under current allocation fractions 
    REAL(wp) :: growth_req_p                      !< moles P required for a unit of N growth under current allocation fractions 
    REAL(wp) :: gdd                               !< growing degree days above 5°C
    REAL(wp) :: nd_dormance                       !< number of days since last growing season 
    REAL(wp) :: t_jmax_opt                        !< current optimum temperature for eletron transport of photosynthesis
    REAL(wp) :: an_boc                            !< net photosynthesis at bottom of canopy (lowest canopy layer existing)
    REAL(wp) :: leaf_cn_direction                 !< indicates direction of change in leaf CN ratio to maximise NPP
    REAL(wp) :: apar 
    !! prior (to r239) the variables ..._cl were stored in t_slm_canopy_layer (without the trailing _cl), which is deleted now
    !!
    REAL(wp),DIMENSION(ncanop) :: nl_tot_cl                !< N in canopy layer (mmol / m2 (LAI))
    REAL(wp),DIMENSION(ncanop) :: fnl_chl_cl               !< fraction of N in chlorophyll (pigments only!)
    REAL(wp),DIMENSION(ncanop) :: fnl_et_cl                !< fraction of N in electron transport 
    REAL(wp),DIMENSION(ncanop) :: fnl_rub_cl               !< fraction of N in Rubisco 
    REAL(wp),DIMENSION(ncanop) :: fnl_pepc_cl              !< fraction of N in PEP carboylase (C4 plants only)
    REAL(wp),DIMENSION(ncanop) :: fnl_oth_cl               !< fraction of N not associated with photosynthesis 
    REAL(wp),DIMENSION(ncanop) :: par_sunlit_cl            !< PAR on sunlit leaves (µmol m-2 s-1)
    REAL(wp),DIMENSION(ncanop) :: par_shaded_cl            !< PAR on shaded leaves (µmol m-2 s-1)
    REAL(wp),DIMENSION(ncanop) :: fleaf_sunlit_cl          !< fraction of layer which is sunlit
    REAL(wp),DIMENSION(ncanop) :: lai_cl                   !< leaf area in the canopy layer 
    REAL(wp),DIMENSION(ncanop) :: cumm_lai_cl              !< cummulative leaf area above the mid-point of the canopy layer 
    !! daytime averages - moved here from t_daily_average at trunk r260 ticket #84
    !! 
    REAL(wp), DIMENSION(ncanop)     :: par_sunlit_daytime_cl             !<  average of the previous day
    REAL(wp), DIMENSION(ncanop)     :: par_shaded_daytime_cl             !<  average of the previous day
    REAL(wp), DIMENSION(ncanop)     :: fleaf_sunlit_daytime_cl           !<  average of the previous day
    REAL(wp)                        :: t_air_daytime                     !<  average of the previous day
    REAL(wp)                        :: press_srf_daytime                 !<  average of the previous day
    REAL(wp)                        :: co2_mixing_ratio_daytime          !<  average of the previous day
    REAL(wp)                        :: ga_daytime                        !<  average of the previous day
    REAL(wp)                        :: beta_air_daytime                  !<  average of the previous day
    REAL(wp)                        :: beta_soil_ps_daytime              !<  average of the previous day
    REAL(wp)                        :: beta_soil_gs_daytime              !<  average of the previous day
    REAL(wp)                        :: t_jmax_opt_daytime                !<  average of the previous day
    REAL(wp), DIMENSION(ncanop)     :: par_sunlit_daytime_dacc_cl        !<  value accumulated over current day
    REAL(wp), DIMENSION(ncanop)     :: par_shaded_daytime_dacc_cl        !<  value accumulated over current day
    REAL(wp), DIMENSION(ncanop)     :: fleaf_sunlit_daytime_dacc_cl      !<  value accumulated over current day
    REAL(wp)                        :: t_air_daytime_dacc                !<  average of the previous day
    REAL(wp)                        :: press_srf_daytime_dacc            !<  average of the previous day
    REAL(wp)                        :: co2_mixing_ratio_daytime_dacc     !<  average of the previous day
    REAL(wp)                        :: ga_daytime_dacc                   !<  average of the previous day
    REAL(wp)                        :: beta_air_daytime_dacc             !<  average of the previous day
    REAL(wp)                        :: beta_soil_ps_daytime_dacc         !<  average of the previous day
    REAL(wp)                        :: beta_soil_gs_daytime_dacc         !<  average of the previous day
    REAL(wp)                        :: t_jmax_opt_daytime_dacc           !<  value accumulated over current day 
    !> @todo t_slm_plant:daytime_counter we may move this counter to a different place, as it is used for multiple processes  *slmdev
    REAL(wp) :: daytime_counter         !< number of timesteps of current day with light available & daytime_avrg variables were added up
    !! moving averages
    !! 
    ! 1.1 memory of conditions for the labile pool dynamics
    REAL(wp) :: gpp_tlabile_mavg                     !< gpp  (µmol m-2 s-1)
    REAL(wp) :: maint_respiration_tlabile_mavg       !< maintenance respiration (µmol m-2 s-1)
    REAL(wp) :: growth_req_n_tlabile_mavg            !< moles N required for a unit of C growth under current allocation fractions
    REAL(wp) :: growth_req_p_tlabile_mavg            !< moles P required for a unit of N growth under current allocation fractions
    ! 1.2 memory of conditions for phenology calculations
    REAL(wp) :: t_air_tphen_mavg                              !< air temp
    REAL(wp) :: beta_soil_gs_tphen_mavg                 !< scaling factor to account for soil moisture constraints on conductance
    ! 1.3 memory for plant uptake (average demand over a couple of days to get rid of diurnal C cycle)
    REAL(wp) :: npp_tuptake_mavg                     !< npp of the last week (µmol m-2 s-1)
    REAL(wp) :: demand_uptake_n_tuptake_mavg         !< scalar of Vmax for root N uptake 
    REAL(wp) :: demand_uptake_p_tuptake_mavg         !< scalar of Vmax for root P uptake
    REAL(wp) :: growth_req_n_tuptake_mavg            !< moles N required for a unit of C growth under current allocation fractions
    REAL(wp) :: growth_req_p_tuptake_mavg            !< moles P required for a unit of N growth under current allocation fractions
    REAL(wp) :: cost_n_fixation_tuptake_mavg         !< moving average of N fixation cost (mol C / mol N)
    ! 1.4 memory of conditions for photosynthesis for the optimisation of foliar N fractions 
    REAL(wp),DIMENSION(ncanop) :: par_sunlit_tfrac_mavg_cl    !< PAR on sunlit leaves
    REAL(wp),DIMENSION(ncanop) :: par_shaded_tfrac_mavg_cl    !< PAR on shaded leaves
    REAL(wp),DIMENSION(ncanop) :: fleaf_sunlit_tfrac_mavg_cl  !< fraction of sunlit leaves
    REAL(wp) :: t_air_tfrac_mavg                              !< air temp
    REAL(wp) :: press_srf_tfrac_mavg                          !< air pressure at surface
    REAL(wp) :: co2_mixing_ratio_tfrac_mavg                   !< co2 mixing ratio of the air
    REAL(wp) :: ga_tfrac_mavg                                 !< aerodynamic resistance
    REAL(wp) :: beta_air_tfrac_mavg                           !< scaling factor to account for air humidity effects on conductance
    REAL(wp) :: beta_soil_ps_tfrac_mavg                       !< scaling factor to account for soil moisture constraints on photosynthesis 
    REAL(wp) :: beta_soil_gs_tfrac_mavg                       !< scaling factor to account for soil moisture constraints on conductance
    ! 1.5 memory of conditions for photosynthesis for the optimisation of foliar N  
    REAL(wp),DIMENSION(ncanop) :: par_sunlit_tcnl_mavg_cl     !< PAR on sunlit leaves
    REAL(wp),DIMENSION(ncanop) :: par_shaded_tcnl_mavg_cl     !< PAR on shaded leaves
    REAL(wp),DIMENSION(ncanop) :: fleaf_sunlit_tcnl_mavg_cl   !< fraction of sunlit leaves
    REAL(wp) :: t_air_tcnl_mavg                               !< air temp
    REAL(wp) :: press_srf_tcnl_mavg                           !< air pressure at surface
    REAL(wp) :: co2_mixing_ratio_tcnl_mavg                    !< co2 mixing ratio of the air
    REAL(wp) :: ga_tcnl_mavg                                  !< aerodynamic resistance
    REAL(wp) :: beta_air_tcnl_mavg                            !< scaling factor to account for air humidity effects on conductance
    REAL(wp) :: beta_soil_ps_tcnl_mavg                        !< scaling factor to account for soil moisture constraints on conductance
    REAL(wp) :: beta_soil_gs_tcnl_mavg                        !< scaling factor to account for soil moisture constraints on conductance
    ! 1.6 memory of conditions for photosynthesis for the optimisation of biomass allocation fractions
    REAL(wp) :: npp_talloc_mavg                 !< NPP (µmol m-2 s-1) (averaging time-scale: tau_allocation)
    REAL(wp) :: transpiration_talloc_mavg       !< average transpiration at allocation time-scale
    REAL(wp) :: unit_npp_talloc_mavg            !< rate of change in NPP with  increment in leaf mass, mumol C mol-1
    REAL(wp) :: unit_transpiration_talloc_mavg  !< rate of tranpiration change per unit LAI change?
    REAL(wp) :: unit_uptake_n_talloc_mavg       !< rate of N uptkae per unit root, mumol N mol-1 (averaging time-scale: tau_allocation)
    REAL(wp) :: unit_uptake_p_talloc_mavg       !< rate of P uptake per unit root, mumol P mol-1 (time-scale of accumulation: tau_allocation)
    REAL(wp) :: dphi_talloc_mavg                !< leaf-soil water potential gradient (MPa)  (averaging time-scale: tau_allocation)
    REAL(wp) :: growth_cn_talloc_mavg           !< plant growth cn ratio, equivalent to 1/growth_req_n
    REAL(wp) :: growth_cp_talloc_mavg           !< plant growth cp ratio, equivalent to 1/(growthreq_n*growth_req_p)
    REAL(wp) :: labile_carbon_talloc_mavg       !< moving average of labile carbon pool neeeded for allometry calculations
    REAL(wp) :: labile_nitrogen_talloc_mavg     !< moving average of labile nitrogen pool neeeded for allometry calculations
    REAL(wp) :: labile_phosphorus_talloc_mavg   !< moving average of labile phosphorus pool neeeded for allometry calculations
    ! 1.7 memory over root life time 
    REAL(wp) :: leaf2root_troot_mavg            !< leaf to root mass ratio
    REAL(wp) :: unit_npp_troot_mavg            !< rate of change in NPP with  increment in leaf mass, mumol C mol-1
    REAL(wp) :: unit_uptake_n_troot_mavg        !< rate of N uptkae per unit root, mumol N mol-1 (time-scale of accumulation: tau_root)
    REAL(wp) :: unit_uptake_p_troot_mavg        !< rate of P uptake per unit root, mumol P mol-1 (time-scale of accumulation: tau_root)
    REAL(wp) :: fmaint_rate_troot_mavg             !< moving average of N-specific maintentance respiration (~µmol C/mol N / s ) 
    ! 1.8 memory over timescale of vegetation dynamics 
    REAL(wp) :: an_boc_tvegdyn_mavg            !< net photosynthesis at bottom of canopy (lowest canopy layer existing)
    ! 1.6 memory of certain variables, temperature, ...
    REAL(wp) :: t_air_week_mavg        !< avrg air temperaure across 7 days
    REAL(wp) :: t_air_month_mavg       !< avrg air temperaure across 30 days
    REAL(wp) :: t_jmax_opt_mavg        !< optimim temperature for electron transport of photosynthesis

    !! fluxes 
    !!
     TYPE(t_element)      :: establishment !< flux between seed bed and reserve (mol m-2 timestep-1)
     TYPE(t_element)      :: reserve_use  !< net flux between reserve and labile pool (mol m-2 timestep-1)
                                          !!  (+ve: flux to labile; -ve flux from labile) 
     TYPE(t_slm_veg_pool) :: growth       !< add to      plant%pool due to growth     per timestep (mol m-2 timestep-1)
     TYPE(t_slm_veg_pool) :: litterfall   !< remove from plant%pool due to litterfall per timestep (mol m-2 timestep-1)
     REAL(wp) :: gpp                      !< gross photosynthesis (µmol CO2 m-2 s-1)
     REAL(wp) :: an                       !< net photosynthesis (µmol CO2 m-2 s-1)
     REAL(wp) :: gpp_c13                  !< gross photosynthesis (µmol 13CO2 m-2 s-1)
     REAL(wp) :: gpp_c14                  !< gross photosynthesis (µmol 14CO2 m-2 s-1)
     REAL(wp) :: maint_respiration_leaf   !< foliar maintenance respiration (µmol CO2 m-2 s-1)
     REAL(wp) :: maint_respiration_pot    !< potential maintenance respiration (µmol CO2 m-2 s-1)
     REAL(wp) :: maint_respiration        !< maintenance respiration (µmol CO2 m-2 s-1)
     REAL(wp) :: maint_respiration_c13    !< maintenance respiration (µmol 13CO2 m-2 s-1)
     REAL(wp) :: maint_respiration_c14    !< maintenance respiration (µmol 14CO2 m-2 s-1)
     REAL(wp) :: growth_respiration       !< growth respiration (µmol CO2 m-2 s-1)
     REAL(wp) :: growth_respiration_c13   !< growth respiration (µmol 13CO2 m-2 s-1)
     REAL(wp) :: growth_respiration_c14   !< growth respiration (µmol 14CO2 m-2 s-1)
     REAL(wp) :: n_transform_respiration  !< respiration associated with N transformation (µmol CO2 m-2 s-1)
     REAL(wp) :: n_fixation_respiration   !< respiration associated with N fixation (µmol CO2 m-2 s-1)
     REAL(wp) :: n_processing_respiration       !< respiration associated with N processing (µmol CO2 m-2 s-1)
     REAL(wp) :: n_processing_respiration_c13   !< respiration associated with N processing (µmol 13CO2 m-2 s-1)
     REAL(wp) :: n_processing_respiration_c14   !< respiration associated with N processing (µmol 14CO2 m-2 s-1)
     REAL(wp) :: npp                      !< net primary production (µmol CO2 m-2 s-1)
     REAL(wp) :: npp_c13                  !< net primary production (µmol 13CO2 m-2 s-1)
     REAL(wp) :: npp_c14                  !< net primary production (µmol 14CO2 m-2 s-1)
     REAL(wp) :: uptake_n                 !< N uptake from soil ( is that correct? )
     REAL(wp) :: uptake_n15               !< 15N uptake from soil ( is that correct? ) 
     REAL(wp) :: uptake_p                 !< P uptake from soil ( is that correct? )
     REAL(wp) :: n_fixation               !< symbiontic N fixation (µmol N m-2 s-1)
     REAL(wp) :: n_fixation_n15           !< symbiontic N fixation (µmol 15N m-2 s-1)
     REAL(wp) :: recycling_leaf_n         !< net flux of N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
     REAL(wp) :: recycling_leaf_n15       !< net flux of 15N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1) 
     REAL(wp) :: recycling_leaf_p         !< net flux of P from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
     REAL(wp) :: recycling_fine_root_n    !< net flux of N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
     REAL(wp) :: recycling_fine_root_n15  !< net flux of 15N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1) 
     REAL(wp) :: recycling_fine_root_p    !< net flux of P from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
     REAL(wp) :: recycling_sap_wood_n     !< net flux of N from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1)
     REAL(wp) :: recycling_sap_wood_n15   !< net flux of 15N from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1) 
     REAL(wp) :: recycling_sap_wood_p     !< net flux of P from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1)
     REAL(wp) :: transpiration            !< transpiration in kg m-2 s-1
     !! 
     !! everything else...
     REAL(wp) :: fmaint_rate              !< current N-specific maintentance respiration (~µmol C/mol N / s ) 
     REAL(wp) :: ci                       !< leaf internal co2 concentration
     REAL(wp) :: unit_npp                 !< incremental change in NPP. Needed for optimal allocation
     REAL(wp) :: unit_transpiration       !< incremental change in transpiration. Needed for optimal allocation
     REAL(wp) :: unit_uptake_n            !< incremental change in N-based growth, needed for optimal allocation
     REAL(wp) :: unit_uptake_p            !< incremental change in P-based growth, needed for optimal allocation
     REAL(wp) :: delta_dens_ind           !< change in density of individuals (#/m2/timestep)
     REAL(wp) :: cost_n_uptake_root       !< actual cost for N uptake by roots (mol C/mol N)
     REAL(wp) :: cost_n_fixation          !< actual cost for N uptake through fixation (mol C/mol N)

  CONTAINS
    PRIVATE
    PROCEDURE, PUBLIC :: init => init_plant
    PROCEDURE, PUBLIC :: set_zero_fluxes => set_zero_veg_flux
   END TYPE t_slm_plant


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> 1.0 addition of two vegetation pools: type + type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION add_veg_pools_type_type(this,t2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this, t2         !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3               !< any t_slm_veg_pool
     t3%leaf          = this%leaf        + t2%leaf         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   + t2%fine_root
     t3%coarse_root   = this%coarse_root + t2%coarse_root
     t3%sap_wood      = this%sap_wood    + t2%sap_wood
     t3%heart_wood    = this%heart_wood  + t2%heart_wood
     t3%labile        = this%labile      + t2%labile
     t3%reserve       = this%reserve     + t2%reserve
     t3%fruit         = this%fruit       + t2%fruit
     t3%seed_bed      = this%seed_bed    + t2%seed_bed
  END FUNCTION add_veg_pools_type_type  

  !-----------------------------------------------------------------------------------------------------
  !> 1.1 addition of two vegetation pools: type + scalar
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION add_veg_pools_type_scalar(this,s2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this        !< any t_slm_veg_pool
     REAL(wp), INTENT(in)              :: s2          !< any scalar
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = this%leaf        + s2         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   + s2
     t3%coarse_root   = this%coarse_root + s2
     t3%sap_wood      = this%sap_wood    + s2
     t3%heart_wood    = this%heart_wood  + s2
     t3%labile        = this%labile      + s2
     t3%reserve       = this%reserve     + s2
     t3%fruit         = this%fruit       + s2
     t3%seed_bed      = this%seed_bed    + s2
  END FUNCTION add_veg_pools_type_scalar  

  !-----------------------------------------------------------------------------------------------------
  !> 1.2 addition of two vegetation pools: scalar + type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION add_veg_pools_scalar_type(s1,this2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     REAL(wp), INTENT(in)              :: s1          !< any scalar
     CLASS(t_slm_veg_pool), INTENT(in) :: this2       !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = s1 + this2%leaf               ! operator overload, see t_element
     t3%fine_root     = s1 + this2%fine_root 
     t3%coarse_root   = s1 + this2%coarse_root
     t3%sap_wood      = s1 + this2%sap_wood  
     t3%heart_wood    = s1 + this2%heart_wood
     t3%labile        = s1 + this2%labile 
     t3%reserve       = s1 + this2%reserve 
     t3%fruit         = s1 + this2%fruit 
     t3%seed_bed      = s1 + this2%seed_bed 
  END FUNCTION add_veg_pools_scalar_type  

  !-----------------------------------------------------------------------------------------------------
  !> 2.0 subtraction of two vegetation pools: type - type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION subtract_veg_pools_type_type(this,t2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this, t2         !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3               !< any t_slm_veg_pool
     t3%leaf          = this%leaf        - t2%leaf         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   - t2%fine_root
     t3%coarse_root   = this%coarse_root - t2%coarse_root
     t3%sap_wood      = this%sap_wood    - t2%sap_wood
     t3%heart_wood    = this%heart_wood  - t2%heart_wood
     t3%labile        = this%labile      - t2%labile
     t3%reserve       = this%reserve     - t2%reserve
     t3%fruit         = this%fruit       - t2%fruit
     t3%seed_bed      = this%seed_bed    - t2%seed_bed
  END FUNCTION subtract_veg_pools_type_type  

  !-----------------------------------------------------------------------------------------------------
  !> 2.1 subtraction of two vegetation pools: type - scalar
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION subtract_veg_pools_type_scalar(this,s2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this        !< any t_slm_veg_pool
     REAL(wp), INTENT(in)              :: s2          !< any scalar
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = this%leaf        - s2         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   - s2
     t3%coarse_root   = this%coarse_root - s2
     t3%sap_wood      = this%sap_wood    - s2
     t3%heart_wood    = this%heart_wood  - s2
     t3%labile        = this%labile      - s2
     t3%reserve       = this%reserve     - s2
     t3%fruit         = this%fruit       - s2
     t3%seed_bed      = this%seed_bed    - s2
  END FUNCTION subtract_veg_pools_type_scalar  

  !-----------------------------------------------------------------------------------------------------
  !> 2.2 subtraction of two vegetation pools: scalar - type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION subtract_veg_pools_scalar_type(s1,this2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     REAL(wp), INTENT(in)              :: s1          !< any scalar
     CLASS(t_slm_veg_pool), INTENT(in) :: this2       !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = s1 - this2%leaf               ! operator overload, see t_element
     t3%fine_root     = s1 - this2%fine_root 
     t3%coarse_root   = s1 - this2%coarse_root
     t3%sap_wood      = s1 - this2%sap_wood  
     t3%heart_wood    = s1 - this2%heart_wood
     t3%labile        = s1 - this2%labile 
     t3%reserve       = s1 - this2%reserve 
     t3%fruit         = s1 - this2%fruit 
     t3%seed_bed      = s1 - this2%seed_bed 
  END FUNCTION subtract_veg_pools_scalar_type  

  !-----------------------------------------------------------------------------------------------------
  !> 3.0 multiplication of two vegetation pools: type * type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION multiply_veg_pools_type_type(this,t2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this, t2         !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3               !< any t_slm_veg_pool
     t3%leaf          = this%leaf        * t2%leaf         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   * t2%fine_root
     t3%coarse_root   = this%coarse_root * t2%coarse_root
     t3%sap_wood      = this%sap_wood    * t2%sap_wood
     t3%heart_wood    = this%heart_wood  * t2%heart_wood
     t3%labile        = this%labile      * t2%labile
     t3%reserve       = this%reserve     * t2%reserve
     t3%fruit         = this%fruit       * t2%fruit
     t3%seed_bed      = this%seed_bed    * t2%seed_bed
  END FUNCTION multiply_veg_pools_type_type  

  !-----------------------------------------------------------------------------------------------------
  !> 3.1 multiplication of two vegetation pools: type * scalar
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION multiply_veg_pools_type_scalar(this,s2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this        !< any t_slm_veg_pool
     REAL(wp), INTENT(in)              :: s2          !< any scalar
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = this%leaf        * s2         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   * s2
     t3%coarse_root   = this%coarse_root * s2
     t3%sap_wood      = this%sap_wood    * s2
     t3%heart_wood    = this%heart_wood  * s2
     t3%labile        = this%labile      * s2
     t3%reserve       = this%reserve     * s2
     t3%fruit         = this%fruit       * s2
     t3%seed_bed      = this%seed_bed    * s2
  END FUNCTION multiply_veg_pools_type_scalar  

  !-----------------------------------------------------------------------------------------------------
  !> 3.2 multiplication of two vegetation pools: scalar * type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION multiply_veg_pools_scalar_type(s1,this2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     REAL(wp), INTENT(in)              :: s1          !< any scalar
     CLASS(t_slm_veg_pool), INTENT(in) :: this2       !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = s1 * this2%leaf               ! operator overload, see t_element
     t3%fine_root     = s1 * this2%fine_root 
     t3%coarse_root   = s1 * this2%coarse_root
     t3%sap_wood      = s1 * this2%sap_wood  
     t3%heart_wood    = s1 * this2%heart_wood
     t3%labile        = s1 * this2%labile 
     t3%reserve       = s1 * this2%reserve 
     t3%fruit         = s1 * this2%fruit 
     t3%seed_bed      = s1 * this2%seed_bed 
  END FUNCTION multiply_veg_pools_scalar_type  

  !-----------------------------------------------------------------------------------------------------
  !> 4.0 division of two vegetation pools: type / type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION divide_veg_pools_type_type(this,t2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this, t2         !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3               !< any t_slm_veg_pool
     t3%leaf          = this%leaf        / t2%leaf         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   / t2%fine_root
     t3%coarse_root   = this%coarse_root / t2%coarse_root
     t3%sap_wood      = this%sap_wood    / t2%sap_wood
     t3%heart_wood    = this%heart_wood  / t2%heart_wood
     t3%labile        = this%labile      / t2%labile
     t3%reserve       = this%reserve     / t2%reserve
     t3%fruit         = this%fruit       / t2%fruit
     t3%seed_bed      = this%seed_bed    / t2%seed_bed
  END FUNCTION divide_veg_pools_type_type  

  !-----------------------------------------------------------------------------------------------------
  !> 4.1 division of two vegetation pools: type / scalar
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION divide_veg_pools_type_scalar(this,s2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this        !< any t_slm_veg_pool
     REAL(wp), INTENT(in)              :: s2          !< any scalar
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = this%leaf        / s2         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   / s2
     t3%coarse_root   = this%coarse_root / s2
     t3%sap_wood      = this%sap_wood    / s2
     t3%heart_wood    = this%heart_wood  / s2
     t3%labile        = this%labile      / s2
     t3%reserve       = this%reserve     / s2
     t3%fruit         = this%fruit       / s2
     t3%seed_bed      = this%seed_bed    / s2
  END FUNCTION divide_veg_pools_type_scalar  

  !-----------------------------------------------------------------------------------------------------
  !> 4.2 division of two vegetation pools: scalar / type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION divide_veg_pools_scalar_type(s1,this2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     REAL(wp), INTENT(in)              :: s1          !< any scalar
     CLASS(t_slm_veg_pool), INTENT(in) :: this2       !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = s1 / this2%leaf               ! operator overload, see t_element
     t3%fine_root     = s1 / this2%fine_root 
     t3%coarse_root   = s1 / this2%coarse_root
     t3%sap_wood      = s1 / this2%sap_wood  
     t3%heart_wood    = s1 / this2%heart_wood
     t3%labile        = s1 / this2%labile 
     t3%reserve       = s1 / this2%reserve 
     t3%fruit         = s1 / this2%fruit 
     t3%seed_bed      = s1 / this2%seed_bed 
  END FUNCTION divide_veg_pools_scalar_type  

  !-----------------------------------------------------------------------------------------------------
  !> 5.0 initialise veg pool with scalar 
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE init_veg_pool(this,s)

    USE mo_shared_memory,   ONLY: t_element

    IMPLICIT NONE
    REAL(wp), INTENT(in)         :: s
    CLASS(t_slm_veg_pool), INTENT(inout) :: this
     CALL this%leaf%init(s)
     CALL this%fine_root%init(s)
     CALL this%coarse_root%init(s)
     CALL this%sap_wood%init(s)
     CALL this%heart_wood%init(s)
     CALL this%labile%init(s)
     CALL this%reserve%init(s)
     CALL this%fruit%init(s)
     CALL this%seed_bed%init(s)
  END SUBROUTINE init_veg_pool 


  !-----------------------------------------------------------------------------------------------------
  !> 5.0 Initialisation all plant variables including pools and canopy layer attributes
  !! 
  !! currently called in the driver \n
  !! for fluxes see mo_veg_memory:set_zero_veg_flux
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE init_plant(this)

    USE mo_atmland_constants,   ONLY: standard_press_srf, def_co2_mixing_ratio

    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_plant), INTENT(inout) :: this 
    ! 0.2 Local
     REAL(wp) :: pool_init, &                ! pools
                 canopy_layer_init, &        ! canopy layers
                 plant_var_init_zero, &      ! plant variables -> growing_season is extra
                 plant_var_init_one          ! plant variables

     pool_init           = 0.0_wp
     canopy_layer_init   = 0.0_wp
     plant_var_init_zero = 0.0_wp
     plant_var_init_one  = 1.0_wp

    ! pools
    CALL this%pool%leaf%init(pool_init)
    CALL this%pool%fine_root%init(pool_init)
    CALL this%pool%coarse_root%init(pool_init)
    CALL this%pool%sap_wood%init(pool_init)
    CALL this%pool%heart_wood%init(pool_init)
    CALL this%pool%labile%init(pool_init)
    CALL this%pool%reserve%init(pool_init)
    CALL this%pool%fruit%init(pool_init)
    CALL this%pool%seed_bed%init(pool_init)
    ! canopy layers
    this%nl_tot_cl(:)       = canopy_layer_init
    this%fnl_chl_cl(:)      = canopy_layer_init
    this%fnl_et_cl(:)       = canopy_layer_init
    this%fnl_rub_cl(:)      = canopy_layer_init
    this%fnl_pepc_cl(:)      = canopy_layer_init
    this%fnl_oth_cl(:)      = canopy_layer_init
    this%par_sunlit_cl(:)   = canopy_layer_init
    this%par_shaded_cl(:)   = canopy_layer_init
    this%fleaf_sunlit_cl(:) = canopy_layer_init
    ! plant variables
    this%growing_season   = .TRUE.     ! LOGICAL
    this%height           = plant_var_init_zero
    this%diameter         = plant_var_init_zero
    this%lai              = plant_var_init_zero
    this%target_cn_leaf   = plant_var_init_zero
    this%target_np_leaf   = plant_var_init_zero
    this%mean_leaf_age    = plant_var_init_zero
    this%gc               = plant_var_init_zero
    this%f_n_demand       = plant_var_init_zero
    this%f_p_demand       = plant_var_init_zero
    this%beta_air         = plant_var_init_one
    this%beta_soil_ps     = plant_var_init_one
    this%beta_soil_gs     = plant_var_init_one
    this%dphi             = plant_var_init_zero
    this%t_jmax_opt       = plant_var_init_zero
    this%growth_req_n     = plant_var_init_zero
    this%growth_req_p     = plant_var_init_zero
    this%gdd                = plant_var_init_one
    this%nd_dormance        = plant_var_init_one
    this%an_boc             = plant_var_init_one
    this%leaf_cn_direction  = plant_var_init_one
    this%apar             = plant_var_init_zero
    !! daytime averages
    this%par_sunlit_daytime_cl    = 0.0_wp    ! this sets the complete array to 0.0_wp
    this%par_shaded_daytime_cl    = 0.0_wp
    this%fleaf_sunlit_daytime_cl  = 0.0_wp
    this%t_air_daytime            = 283.15_wp
    this%press_srf_daytime        = standard_press_srf 
    this%co2_mixing_ratio_daytime = def_co2_mixing_ratio 
    this%ga_daytime               = 0.2_wp
    this%beta_air_daytime         = 1.0_wp
    this%beta_soil_ps_daytime     = 1.0_wp
    this%beta_soil_gs_daytime     = 1.0_wp 
    this%t_jmax_opt_daytime       = 17.0_wp   ! set here to avoid referencing assim_constants
    this%par_sunlit_daytime_dacc_cl     = 0.0_wp
    this%par_shaded_daytime_dacc_cl     = 0.0_wp
    this%fleaf_sunlit_daytime_dacc_cl   = 0.0_wp
    this%t_air_daytime_dacc             = 283.15_wp
    this%press_srf_daytime_dacc         = standard_press_srf 
    this%co2_mixing_ratio_daytime_dacc  = def_co2_mixing_ratio 
    this%ga_daytime_dacc                = 0.2_wp
    this%beta_air_daytime_dacc          = 1.0_wp
    this%beta_soil_ps_daytime_dacc      = 1.0_wp
    this%beta_soil_gs_daytime_dacc      = 1.0_wp 
    this%t_jmax_opt_daytime_dacc        = 0.0_wp
    this%daytime_counter  = plant_var_init_zero
    !! moving averages
    ! 1.1. labile pool 
    this%gpp_tlabile_mavg                   = 0.0_wp
    this%maint_respiration_tlabile_mavg     = 0.0_wp
    this%growth_req_n_tlabile_mavg          = 1._wp / 25._wp 
    this%growth_req_p_tlabile_mavg          = 1._wp / 14._wp
    ! 1.2 phenology
    this%t_air_tphen_mavg                   = 283.15_wp
    this%beta_soil_gs_tphen_mavg            = 1.0_wp
    ! 1.3 nutrient demand for plant uptake
    this%npp_tuptake_mavg                   = 1.0_wp
    this%demand_uptake_n_tuptake_mavg       = 0.1_wp
    this%demand_uptake_p_tuptake_mavg       = 0.1_wp
    this%growth_req_n_tuptake_mavg          = 1._wp / 25._wp 
    this%growth_req_p_tuptake_mavg          = 1._wp / 14._wp
    this%cost_n_fixation_tuptake_mavg       = 10._wp
    ! 1.4 input for optimal fraction calculation 
    this%par_sunlit_tfrac_mavg_cl           = 150.0_wp
    this%par_shaded_tfrac_mavg_cl           = 10.0_wp
    this%fleaf_sunlit_tfrac_mavg_cl         = 0.5_wp
    this%t_air_tfrac_mavg                   = 283.15_wp   ! 10°C; temperature in K
    this%press_srf_tfrac_mavg               = standard_press_srf      ! mo_atmland_constants
    this%co2_mixing_ratio_tfrac_mavg        = def_co2_mixing_ratio    ! mo_atmland_constants
    this%ga_tfrac_mavg                      = 0.2_wp
    this%beta_air_tfrac_mavg                = 1.0_wp
    this%beta_soil_ps_tfrac_mavg            = 1.0_wp
    this%beta_soil_gs_tfrac_mavg            = 1.0_wp
    ! 1.5 input for optimal N calculation  
    this%par_sunlit_tcnl_mavg_cl            = 150.0_wp
    this%par_shaded_tcnl_mavg_cl            = 10.0_wp
    this%fleaf_sunlit_tcnl_mavg_cl          = 0.5_wp
    this%t_air_tcnl_mavg                    = 283.15_wp   ! 10°C; temperature in K
    this%press_srf_tcnl_mavg                = standard_press_srf      ! mo_atmland_constants
    this%co2_mixing_ratio_tcnl_mavg         = def_co2_mixing_ratio    ! mo_atmland_constants
    this%ga_tcnl_mavg                       = 0.2_wp
    this%beta_air_tcnl_mavg                 = 1.0_wp
    this%beta_soil_ps_tcnl_mavg             = 1.0_wp
    this%beta_soil_gs_tcnl_mavg             = 1.0_wp
    ! 1.6 input for optimal allocation
    this%npp_talloc_mavg                    = 1._wp
    this%transpiration_talloc_mavg          = 0._wp
    this%unit_npp_talloc_mavg               = 0._wp
    this%unit_transpiration_talloc_mavg     = 0._wp
    this%unit_uptake_n_talloc_mavg          = 0._wp
    this%unit_uptake_p_talloc_mavg          = 0._wp
    this%dphi_talloc_mavg                   = -2._wp
    this%growth_cn_talloc_mavg              = 1._wp
    this%growth_cp_talloc_mavg              = 1._wp
    this%labile_carbon_talloc_mavg          = 0._wp
    this%labile_nitrogen_talloc_mavg        = 0._wp
    this%labile_phosphorus_talloc_mavg      = 0._wp
    ! 1.7 root life time avearges
    this%leaf2root_troot_mavg               = 1._wp
    this%unit_npp_troot_mavg                = 0._wp
    this%unit_uptake_n_troot_mavg           = 0._wp
    this%unit_uptake_p_troot_mavg           = 0._wp
    this%fmaint_rate_troot_mavg             = 0._wp
    ! 1.8 memory of vegetation dynamics 
    this%an_boc_tvegdyn_mavg                = 5.0_wp     
    ! 1.9 memory of certain variables, temperature [K], ...
    this%t_air_week_mavg                    = 283.15_wp   ! 10°C; 7  days
    this%t_air_month_mavg                   = 283.15_wp   ! 10°C; 30 days
    this%t_jmax_opt_mavg                    = 17._wp   ! 17°C

  END SUBROUTINE init_plant 

  !-----------------------------------------------------------------------------------------------------
  !> initalize or just set zero all fluxes vegetation memory (plant)
  !!
  !! is initially called in the driver, 
  !! and then called every timestep in beginning of mo_veg_interface -> move to 1st Task of veg called
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE set_zero_veg_flux(this)

    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_plant),INTENT(inout)  :: this    !< plant
    ! 0.2 Local
    REAL(wp) :: flux_zero

    flux_zero = 0.0_wp

    CALL this%reserve_use%init(0.0_wp)   ! t_element init
    CALL this%establishment%init(0.0_wp) ! t_element init
    CALL this%growth%init(0.0_wp)        ! t_slm_veg_pool init
    CALL this%litterfall%init(0.0_wp)    ! t_slm_veg_pool init

    this%gpp                     = flux_zero  !< gross photosynthesis (µmol CO2 m-2 s-1)
    this%an                      = flux_zero  !< net photosynthesis (µmol CO2 m-2 s-1)
    this%gpp_c13                 = flux_zero  !< gross photosynthesis (µmol 13CO2 m-2 s-1)
    this%gpp_c14                 = flux_zero  !< gross photosynthesis (µmol 14CO2 m-2 s-1)
    this%maint_respiration_leaf  = flux_zero  !< foliar maintenance respiration (µmol CO2 m-2 s-1)
    this%maint_respiration_pot   = flux_zero  !< maintenance respiration (µmol CO2 m-2 s-1)
    this%maint_respiration       = flux_zero  !< maintenance respiration (µmol CO2 m-2 s-1)
    this%maint_respiration_c13   = flux_zero  !< maintenance respiration (µmol 13CO2 m-2 s-1)
    this%maint_respiration_c14   = flux_zero  !< maintenance respiration (µmol 14CO2 m-2 s-1)
    this%growth_respiration      = flux_zero  !< growth respiration (µmol CO2 m-2 s-1)
    this%growth_respiration_c13  = flux_zero  !< growth respiration (µmol 13CO2 m-2 s-1)
    this%growth_respiration_c14  = flux_zero  !< growth respiration (µmol 14CO2 m-2 s-1)
    this%n_transform_respiration = flux_zero  !< respiration associated with N transformation (µmol CO2 m-2 s-1)
    this%n_fixation_respiration  = flux_zero  !< respiration associated with N fixation (µmol CO2 m-2 s-1)
    this%n_processing_respiration  = flux_zero  !< respiration associated with N fixation (µmol CO2 m-2 s-1)
    this%n_processing_respiration_c13  = flux_zero  !< respiration associated with N fixation (µmol 13CO2 m-2 s-1)
    this%n_processing_respiration_c14  = flux_zero  !< respiration associated with N fixation (µmol 14CO2 m-2 s-1)
    this%npp                     = flux_zero  !< net primary production (µmol CO2 m-2 s-1)
    this%npp_c13                 = flux_zero  !< net primary production (µmol 13CO2 m-2 s-1)
    this%npp_c14                 = flux_zero  !< net primary production (µmol 14CO2 m-2 s-1)
    this%uptake_n                = flux_zero  !< ??
    this%uptake_n15              = flux_zero  !< ??
    this%n_fixation              = flux_zero  !< ??
    this%n_fixation_n15          = flux_zero  !< ??
    this%uptake_p                = flux_zero  !< ??
    this%recycling_leaf_n        = flux_zero  !< net flux of N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
    this%recycling_leaf_n15      = flux_zero  !< net flux of 15N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1) 
    this%recycling_leaf_p        = flux_zero  !< net flux of P from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
    this%recycling_fine_root_n   = flux_zero  !< net flux of N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
    this%recycling_fine_root_n15 = flux_zero  !< net flux of 15N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1) 
    this%recycling_fine_root_p   = flux_zero  !< net flux of P from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
    this%recycling_sap_wood_n    = flux_zero  !< net flux of N from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1)
    this%recycling_sap_wood_n15  = flux_zero  !< net flux of 15N from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1) 
    this%recycling_sap_wood_p    = flux_zero  !< net flux of P from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1)
    this%transpiration           = flux_zero  !< transpiration in kg m-2 s-1
    this%ci                      = flux_zero  !< ??
    this%delta_dens_ind          = flux_zero  !< ??
    this%unit_transpiration      = flux_zero
  END SUBROUTINE set_zero_veg_flux


END MODULE mo_veg_memory
