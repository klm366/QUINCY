!> QUINCY vegetation constants
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### declare and define vegetation constants
!>
MODULE mo_veg_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: jp

  IMPLICIT NONE
  PUBLIC


  ! Respiration parameters
  REAL(jp), SAVE :: &
      ! growth respiration per unit new biomass
      fresp_growth          , &   !< growth respiration per unit new biomass
      ! N-specific foliar maintenance respiration rate (micro-mol CO2 / mol N / s)
      ! Taken from Sprugel et al. Respiration from the organ level to the stand.
      ! In: Smith WK, Hinckley, T.M., editor. Physiological Ecology of Coniferous Forests.
      ! San Diego, CA, USA: Academic Press; 1996.
      ! Note converted from 7.4 e-7 to 1.25 e-6 to account for the different values of the Lloyd&Taylor equation
      ! such that the applied rates a comparable at 15-20degC
      fmaint_rate_base      , &   !< N-specific maintenance respiration rate for roots and leaves (micro-mol CO2 / mmol N / s)
      fmaint_rate_wood      , &   !< N-specific maintenance respiration rate for wood (micro-mol CO2 / mmol N / s)
      ! parameters in the Lloyd and Taylor temperature response function (Lloyd & Taylor, 1996, Funct. Ecol.)
      t_maint_k1            , &   !< Coefficient in the equation for temperature sensitivity of foliar respiration
      t_maint_k2            , &   !< Coefficient in the equation for temperature sensitivity of foliar respiration
      t_maint_k3            , &   !< Coefficient in the equation for temperature sensitivity of foliar respiration
      !
      f_resp_acclim         , &   !< respiration temperature acclimation factor
      t_acclim_zero         , &   !< base temperature for respiration acclimation
      k_maint_resp_suppress       !< factor reduce maintentance respiration [unitless]

  ! Parameters describing the PS ~ An relationship (describing photosynthetic activity )
  REAL(jp), SAVE :: &
      k0_fn_pepc_C3     , &  !< fraction of N in PEP and PPDK (C3 plants)
      k0_fn_pepc_C4     , &  !< fraction of N in PEP and PPDK (C4 plants) (Makino et al. 2003 / Tazoe et al. PCE 2006, 29, 691-700)
      ! in static leaf N mode: prescribed ratio of Jmax25/Vcmax25
      jmax2vcmax_C3     , &  !< text book; ratio of Jmax25/Vcmax25 for prescribed leaf N (C3 plants)
      jmax2vcmax_C4          !< calibrated to have Jmax=vcmax at 25 deg; ratio of Jmax25/Vcmax25 for prescribed leaf N (C4 plants)

  ! Labile pool parameters
  REAL(jp), SAVE :: &
      ! default turnover time of the labile carbo-hydrates (in days)
      tau_labile,           & !< default turnover time of the labile carbo-hydrates (in days)
      ! ratio between the pull of carbon from the reserve pool
      ! and carbon allocation to growth (both w.r.t. labile pool turnover rate)
      k_labphen_tree,       & !< ration between pull from reserve pool and allocation to growth for trees
      k_labphen_grass,      & !< ration between pull from reserve pool and allocation to growth for grass
      ! temperature scalar on turnover time of labile pool (1-exp(-(lambda * T)^k))
      ! tuned to reduce growth between 5degC to zero at 0degC
      lambda_labile_temp,   & !< lambda in temperature response function of labile pool
      k_labile_temp,        & !< k in temperature response function of labile pool
      ! moisture scalar on turnover time of labile pool (1-exp(-(lambda * theta)^k))
      ! tuned to reduce growth from 20% saturation to zero at 0% saturation
      lambda_labile_theta,  & !< lambda in moisture response function of labile pool
      k_labile_theta,       & !< k in moisture response function of labile pool
      ! pull from reserve to labile pool due to phenology (1-exp(-(lambda * lai/target_lai)^k))
      lambda_phiphen,       & !< lambda in storage response function to phenology
      k_phiphen,            & !< k in storage response function to phenology
      ! pull from reserve to labile pool due to maintenance of labile pool (1-exp(-(lambda * labile/target_labile)^k))
      lambda_phimaint,      & !< lambda in storage response function to maintenance needs
      k_phimaint,           & !< k in storage response function to maintenance needs
      lambda_phimaint_nut,  & !< lambda in storage response function to maintenance needs
      k_phimaint_nut,       & !< k in storage response function to maintenance needs
      ! pull from reserve to labile pool due to maintenance of labile pool (1-exp(-(lambda * labile/target_labile)^k))
      lambda_phistore,      & !< lambda in labile response function to storage size
      k_phistore,           & !< k in labile response function to storage size
      ! sink limitation signal strength to photosynthesis (exp(-(lambda * labile/target_labile)^k))
      knut_labile,          & !< rate at which nutrients can retrieved quicker for growth than carbon
      lambda_sinklim_ps,    & !< lambda in photosynthesis response function to labile pool size
      k_sinklim_ps,         & !< k in photosynthesis response function to labile pool size
      k_cnp_sinklim_ps,     & !< scaling factor to limit PS if labile CNP becomes too wide
      beta_sinklim_ps_min,  & !< maximum downacclimation of photosynthetic capacity due to sink limitation
      fstore_leaf_max,      & !< maximum reserve storage in leaves relative to leaf mass
      fstore_fine_root_max, & !< maximum reserve storage in fine roots relative to fine root mass
      fstore_sap_wood_max,  & !< maximum reserve storage in sap wood relative to sap wood mass
      k_phi_interact          !< value of reserve pull phi_X, at which the inverse reserve pull phi_Y is reduced

  ! Resorption parameters
  REAL(jp), SAVE :: &
      resorp_fract_leaf,    & !< default fraction of nutrient resporption before leaf shedding
      resorp_fract_wood,    & !< default fraction of nutrient resporption before wood death
      root2leaf_cn,         & !< relative C:N of fine roots compared to leaves (mol C / mol N)
      wood2leaf_cn,         & !< relative C:N of sapwood biomass compared to leaves (mol C / mol N)
      root2leaf_np,         & !< relative N:P of fine roots compared to leaves (mol N / mol P)
      wood2leaf_np,         & !< relative N:P of sapwood biomass compared to leaves (mol N / mol P)
      tau_nutrient_recycling  !< time scale of enzyme, RNA etc. turnover (as in return to labile pool!, days)

  ! Parameters concerning N uptake calculations
  REAL(jp), SAVE :: &
      kappa_f_demand_nc    , & !< fraction of maximum labile nitrogen concentration at which uptake is reduced to 50%
      kappa_f_demand_pn    , & !< fraction of maximum labile phosphorus  concentration at which uptake is reduced to 50%
      k_f_demand           , & !< k in nutrient uptake response funktion to labile nutrient concentration
      omega_nutrient_demand, & !< the safety ratio margin when comparing the targeted and real labile N:C or P:C (unitless)
      k1_temp_bnf          , & !< k1 .. k4 factors in temperature response function of nitrogenase activity
      k2_temp_bnf          , &
      k3_temp_bnf          , &
      k4_temp_bnf          , &
      transform_cost_nh4   , & !< N transformation costs of converting NH4 to organic N  (mol C / mol N)
      transform_cost_no3   , & !< N transformation costs of converting NO3 to organic N (mol C / mol N)
      transform_cost_org   , & !< N transformation costs of converting organic N to organic N (mol C / mol N)
      f_nacq_no3           , & !< fractions of NO3 and NH4 uptake based on transformation costs
      f_nacq_nh4           , &
      k_cost_symb_bnf      , & !< cost of symbiotic N fixation (mol C / mol N)
      t_opt_symb_bnf       , & !< optimal temperature for BNF (K)
      k_topt_symb_bnf      , & !< slope of optimal temperature for symb. BNF vs. growth temperature
      tacclim_nacq_min     , & !< minimum temperature to which N acquisition acclimates (degC)
      tacclim_nacq_max     , & !< maximum temperature to which N acquisition acclimates (degC)
      ea_symb_bnf          , & !< activation energy for symbiotic N fixation (J mol-1)
      k_ea_symb_bnf        , & !< slope of activation energy for symb. BNF vs. growth temperature
      ed_symb_bnf          , & !< de-activation energy for symbiotic N fixation (J mol-1)
      eta_nfixation        , & !< fractionation of N15 due to N fixation (per mill), parameter value from Robinson, Tree 2001
      f_myc_exudation_max  , & !< maximum fraction of current growth allocated to mycorrhiza
      f_mycorrhization_min , & !< minimum root mycorrhization
      f_mycorrhization_max     !< maximum root mycorrhization

  ! parameters to describe the varying fraction of N in chlorophyll through the canopy
  REAL(jp), SAVE :: &
      k0_fn_chl_C3  ,   &  !< parameter to describe the varying fraction of N in chlorophyll through the canopy
      k0_fn_chl_C4  ,   &  !< parameter to describe the varying fraction of N in chlorophyll through the canopy
      k1_fn_chl_C3  ,   &  !< parameter to describe the varying fraction of N in chlorophyll through the canopy
      k1_fn_chl_C4  ,   &  !< parameter to describe the varying fraction of N in chlorophyll through the canopy
      kfn_chl              !< parameter to describe the varying fraction of N in chlorophyll through the canopy

  ! parameters used for phenology calculations
  REAL(jp), SAVE :: &
      min_lai                , & !< minimum LAI beyond which all leaves are shed
      target_lai_max         , & !< maximum LAI target for reserve use calculations
      max_leaf_shedding_rate     !< maximum rate of leaf shedding (1/days)

  ! vegetation dynamics
  REAL(jp), SAVE :: &
      background_mort_rate_tree   , & !< background mortality rate for trees in static mortality case (1/year)
      background_mort_rate_grass  , & !< background mortality rate for grass in static mortality case (1/year)
      lambda_mort_light           , & !< parameter in the Weibull function controlling self-thinning mortality
      k_mort_light                , & !< parameter in the Weibull function controlling self-thinning mortality
      lambda_est_temp             , & !< parameter in the Weibull function controlling temperature-limited establishment
      k_est_temp                  , & !< parameter in the Weibull function controlling temperature-limited establishment
      lambda_est_moist            , & !< parameter in the Weibull function controlling moisture-limited establishment
      k_est_moist                 , & !< parameter in the Weibull function controlling moisture-limited establishment
      min_greff                   , & !< minimum growth efficiency (mol C / m2 (LAI) / yr)
      k2_mort_greff               , & !< scaling coefficient
      k3_mort_greff               , & !< scaling coefficient in growth efficiency mortality below minimum growth efficiency
      max_crown_area              , & !< maximum crown area (m2)
      min_diameter                , & !< minimum diameter for crown area calculation (m)
      k_crown_area                , & !< scaling parameter in crownarea to diameter relationship
      k_rp                        , & !< scaling exponent in crownarea to diameter relationship
      k_fpc                       , & !< ... see ticket #475
      fpc_max                         !< maximum foliage projective cover

  ! other
  REAL(jp), SAVE :: &
      ! carbon content of leaves per dry matter
      carbon_per_dryweight_leaf   , & !< carbon mass per unit dry weight of leaves (gC / g DW)
      ! minimum permissable height of an individual
      min_height                  , & !< minimum permissable height of an individual
      ! stem area to leaf area ratio for trees
      k_sai2lai                   , & !< stem area to leaf area ratio for trees (--)
      ! stem mass to leaf mass ratio of grasses
      sm2lm_grass                 , & !< stem mass to leaf mass ratio of grasses
      ! Parameters controlling fruit allocation in: exp(-(lambda * reserve_use+offset)^k)
      lambda_fruit_alloc          , & !< parameter in the Weibull function controlling fruit growth
      k4_fruit_alloc              , & !< parameter in the Weibull function controlling fruit growth
      k3_fruit_alloc              , & !< parameter in the Weibull function controlling fruit growth; ~micro-mol C / m2 / s
      ! minimum fraction of allocation going to fruit (flowers etc.)
      k1_fruit_alloc              , & !< minimum fraction of allocation going to fruit (flowers etc.)
      k1_root_alloc               , & !< value of relative root zone water below which root allocation starts responding
      leaf2root_min_ratio         , & !< minimum leaf to root ratio for dynamic allocation
      leaf2root_max_ratio             !< maximum leaf to root ratio for dynamic allocation

  ! parameters for changes in canopy N fractions
  REAL(jp), SAVE :: &
      ! maximum change in photosynthetic N fractions
      delta_n_fraction    , & !< maximum fractional change for optimal adjustment of photosynthetic N fraction
      ! maximum change in CN ratio of canopy
      delta_n_leaf            !< maximum fractional change for optimal adjustment of canopy N allocation

  ! parameters used for product pool decay
  REAL(jp), SAVE :: &
      tau_pp_crop,       & !< exponential time constants for the decay of the crop harvest product pool
      tau_pp_fuelwood,   & !< exponential time constants for the decay of the fuel product pool
      tau_pp_paper,      & !< exponential time constants for the decay of the paper product pool
      tau_pp_fiberboard, & !< exponential time constants for the decay of the fiberboard product pool
      tau_pp_oirw,       & !< exponential time constants for the decay of the other IRW product pool
      tau_pp_pv,         & !< exponential time constants for the decay of the plywood and veneer product pool
      tau_pp_sawnwood      !< exponential time constants for the decay of the sawnwood product pool

  ! parameters used for herbivory
  REAL(jp), SAVE :: &
    & ftroph_loss_herbivory, &            !< carbon respiration associated with herbivory (unitless)
    & k_herbivory_grass,     &            !< asymptotic daily rate of herbivory in natural grasslands (1/day)
    & k_herbivory_pasture,   &            !< asymptotic daily rate of herbivory in pasture (1/day)
    & min_net_growth_herbivory_grass, &   !< minimum long-term average productivity for herbivory to occur in natural grasslands (mol / m2 / yr)
    & min_net_growth_herbivory_pasture, & !< minimum long-term average productivity for herbivory to occur in pasture (mol / m2 / yr)
    & min_lai_herbivory                   !< minimum LAI beyond which herbivory can occur

  ! red to far red ratio calculations
  REAL(jp), SAVE :: &
    & rfr_ratio_toc, &     !< top of the canopy red:farred ratio  (Kull & Kruij 1998)
    & k_r2fr_chl           !< extinction of red:farred ratio on leaf chlorophyll (1/mmol N);
                           !! (Kull & Krujit 1998 for the far red value, tuned to get a R:FR ratio
                           !! of 0.15 at full canopy closure (a common value)

  ! some useful numbers for easy identification of photosynthetic pathways and growth forms
  ! need to be parameters (i.e. constant across model runtime)
  INTEGER, PARAMETER ::    &
      itree         = 1,   &
      igrass        = 2,   &
      icrop         = 3


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_veg_constants'

  !$ACC DECLARE CREATE(t_maint_k1, t_maint_k2, t_maint_k3, min_height, k_fpc, k_sai2lai)
  !$ACC DECLARE CREATE(f_resp_acclim, t_acclim_zero, fmaint_rate_base)

  !$ACC DECLARE CREATE(carbon_per_dryweight_leaf, jmax2vcmax_C4, k0_fn_chl_C4, k1_fn_chl_C4, k0_fn_pepc_C4)
  !$ACC DECLARE CREATE(jmax2vcmax_C3, k0_fn_chl_C3, k1_fn_chl_C3, k0_fn_pepc_C3, kfn_chl)

  !$ACC DECLARE CREATE(rfr_ratio_toc, k_r2fr_chl)

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants used in the vegetation process
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_veg_constants

    USE mo_jsb_math_constants,     ONLY: one_day, one_year
    USE mo_jsb_physical_constants, ONLY: molar_mass_C, molar_mass_N, molar_mass_P, Tzero

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    !
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_veg_constants'


  ! Respiration parameters
      fresp_growth                = 0.25_jp
      fmaint_rate_base            = 1._jp
      fmaint_rate_wood            = 0.25_jp * fmaint_rate_base
      t_maint_k1                  = 308.56_jp
      t_maint_k2                  = 56.02_jp
      t_maint_k3                  = 227.13_jp
      f_resp_acclim               = -0.008_jp
      t_acclim_zero               = 283.15_jp
      k_maint_resp_suppress       = 1.0_jp        ! set to 1.0 on purpose, not reducing maint_respiration; ticket https://projects.bgc-jena.mpg.de/QUINCY/ticket/594

  ! Parameters describing the PS ~ An relationship (describing photosynthetic activity )
      k0_fn_pepc_C3               = 0.0_jp
      k0_fn_pepc_C4               = 0.045_jp
      jmax2vcmax_C3               = 1.92_jp / 4._jp
      jmax2vcmax_C4               = 1.4_jp

  ! Labile pool parameters
      tau_labile                  = 5.0_jp
      k_labphen_grass             = 0.5_jp
      k_labphen_tree              = 0.25_jp
      lambda_labile_temp          = 0.5_jp
      k_labile_temp               = 2.0_jp
      lambda_labile_theta         = 10.0_jp
      k_labile_theta              = 2.0_jp
      lambda_phiphen              = 1.3_jp
      k_phiphen                   = 8.0_jp
      lambda_phimaint             = 4.0_jp
      k_phimaint                  = 1.2_jp
      lambda_phimaint_nut         = 1.6_jp
      k_phimaint_nut              = 3._jp
      lambda_phistore             = 2._jp
      k_phistore                  = 3._jp
      knut_labile                 = 1.2_jp
      lambda_sinklim_ps           = 0.1_jp
      k_sinklim_ps                = 2.0_jp
      k_cnp_sinklim_ps            = 4.0_jp
      beta_sinklim_ps_min         = 0.25_jp
      fstore_leaf_max             = 0.02_jp
      fstore_fine_root_max        = 0.2_jp
      fstore_sap_wood_max         = 0.15_jp
      k_phi_interact              = 0.75_jp


  ! Resorption parameters
      resorp_fract_leaf           = 0.5_jp
      resorp_fract_wood           = 0.2_jp
      root2leaf_cn                = 0.85_jp
      wood2leaf_cn                = 0.145_jp
      root2leaf_np                = 1.0_jp
      wood2leaf_np                = 1.0_jp
      tau_nutrient_recycling      = 10._jp

  ! N uptake parameters
      kappa_f_demand_nc           = 0.75_jp
      kappa_f_demand_pn           = 0.9_jp
      k_f_demand                  = 2.0_jp
      omega_nutrient_demand       = 0.15_jp ! tuned to leave a flexiable space between the targeted and actual labile CN or NP ratios
      k1_temp_bnf                 = 1.25_jp
      k2_temp_bnf                 = -3.62_jp
      k3_temp_bnf                 = 0.27_jp
      k4_temp_bnf                 = 50.3_jp
      transform_cost_nh4          = 1.7_jp    / (molar_mass_C/molar_mass_N)
      transform_cost_no3          = 2.3_jp    / (molar_mass_C/molar_mass_N)
      transform_cost_org          = 1.0_jp    / (molar_mass_C/molar_mass_N)
      f_nacq_nh4                  = transform_cost_no3 / (transform_cost_no3 + transform_cost_nh4)
      f_nacq_no3                  = transform_cost_nh4 / (transform_cost_no3 + transform_cost_nh4)
      k_cost_symb_bnf             = 12.75_jp  / (molar_mass_C/molar_mass_N)     ! Average value from Kern 2021
      t_opt_symb_bnf              = 28.15_jp + Tzero ! calibrated from Bytnerowicz et al. 2022, adjusted for acclimation of growth temperature
      k_topt_symb_bnf             = 0.15_jp          ! calibrated based on Bytnerowicz et al. 2022, mean between measured biomes
      k_ea_symb_bnf               = 0.04_jp          ! calibrated from Bytnerowicz et al. 2022
      ea_symb_bnf                 = 72000._jp        ! calibrated from Bytnerowicz et al. 2022
      ed_symb_bnf                 = 250000._jp       ! calibrated from Bytnerowicz et al. 2022
      tacclim_nacq_min            = 15.0_jp
      tacclim_nacq_max            = 30.0_jp
      eta_nfixation               = 3.0_jp                   ! Robinson, Tree 2001
      f_myc_exudation_max         = 0.3_jp
      f_mycorrhization_min        = 0.1_jp
      f_mycorrhization_max        = 0.5_jp

  ! parameters to describe the varying fraction of N in chlorophyll through the canopy
      k0_fn_chl_C3                = 6.0_jp
      k0_fn_chl_C4                = 15.0_jp
      k1_fn_chl_C3                = 3.6_jp
      k1_fn_chl_C4                = 4.4_jp
      kfn_chl                     = 0.7_jp

  ! parameters used for phenology calculations
      min_lai                     = 0.1_jp
      target_lai_max              = 5.0_jp
      max_leaf_shedding_rate      = 0.05_jp

  ! vegetation dynamics
      background_mort_rate_tree   = 0.01_jp
      background_mort_rate_grass  = 0.05_jp
      lambda_mort_light           = 0.15_jp
      k_mort_light                = 4.0_jp
      lambda_est_temp             = 0.075_jp
      k_est_temp                  = 4.0_jp
      lambda_est_moist            = 10._jp
      k_est_moist                 = 2.0_jp
      k2_mort_greff               = 0.3_jp
      k3_mort_greff               = 25.0_jp
      min_greff                   = 1.0_jp
      max_crown_area              = 100._jp ! adjusted to conform with observations by Pretsch et al. 2015; but cutting off extremes
      min_diameter                = 0.01_jp
      k_crown_area                = 175._jp ! adjusted to medium size trees analysed by Pretch et al. 2015; https://doi.org/10.1016/j.ufug.2015.04.006
      k_rp                        = 1.6_jp
      k_fpc                       = 0.5_jp
      fpc_max                     = 0.95_jp

  ! other
      carbon_per_dryweight_leaf   = 0.48_jp
      min_height                  = 0.1_jp
      k_sai2lai                   = 0.05_jp ! calibrated to get fapar < 0.2 when no leaves at canopy closure
      sm2lm_grass                 = 0.05_jp
      lambda_fruit_alloc          = 10.0_jp
      k1_fruit_alloc              = 0.01_jp
      k3_fruit_alloc              = 0.1_jp
      k4_fruit_alloc              = 2.0_jp
      k1_root_alloc               = 0.8_jp
      leaf2root_min_ratio         = 0.2_jp
      leaf2root_max_ratio         = 15.0_jp

  ! parameters for changes in canopy N fractions
!@todo:should this be 1/20days?
      delta_n_fraction            = 0.048_jp        ! 0.001_jp (changed to account for multiplication with dtime/one_day in the code)
      delta_n_leaf                = 0.0048_jp       ! 0.0001 (changed to account for multiplication with dtime/one_day in the code)

    ! parameters used for product pool decay [years] - these need to be > 0.0_jp
      tau_pp_crop = 1.0_jp ! made up by SZ
      ! following values are all taken from Nuetzel (2021, LMU master thesis)
      tau_pp_fuelwood = 0.53_jp
      tau_pp_paper = 1.50_jp
      tau_pp_fiberboard = 13.2_jp
      tau_pp_oirw = 16.5_jp
      tau_pp_pv = 24.9_jp
      tau_pp_sawnwood = 49.9_jp

    ! parameters used for herbivory
    ! calibrated from data in McNaugton et al. 1989, Nature 341, 142-144, and Cyr and Pace, 1993, Nature 361, 148-150
      ftroph_loss_herbivory             = 0.9_jp   !< assumuming a typical 90% loss of energy per trophic level
      k_herbivory_grass                 = 0.01_jp  !< calibrated from Cyr and Pace, Nature 1993, adjusted to daily rate and fraction of accessible NPP
      k_herbivory_pasture               = 0.015_jp !< calibrated from Cyr and Pace, Nature 1993, adjusted to daily rate and fraction of accessible NPP
      min_net_growth_herbivory_grass    = 20.8_jp  !< 250 gC/m2/yr converted to mol
      min_net_growth_herbivory_pasture  = 4.16_jp  !< 50 gC/m2/yr converted to mol
      min_lai_herbivory                 = 0.1_jp   !< set to min_lai value, needs evaluation

    ! red to far red ratio calculations
      rfr_ratio_toc           = 1.2_jp                !< Kull & Kruijt 1998
      k_r2fr_chl              = 0.0012_jp - 0.012_jp  !< Kull & Kruijt 1998

    !$ACC UPDATE DEVICE(t_maint_k1, t_maint_k2, t_maint_k3, min_height) ASYNC(1)
    !$ACC UPDATE DEVICE(f_resp_acclim, t_acclim_zero, fmaint_rate_base, k_fpc, k_sai2lai) ASYNC(1)
    !$ACC UPDATE DEVICE(carbon_per_dryweight_leaf, jmax2vcmax_C4, k0_fn_chl_C4, k1_fn_chl_C4, k0_fn_pepc_C4) ASYNC(1)
    !$ACC UPDATE DEVICE(jmax2vcmax_C3, k0_fn_chl_C3, k1_fn_chl_C3, k0_fn_pepc_C3, kfn_chl) ASYNC(1)
    !$ACC UPDATE DEVICE(rfr_ratio_toc, k_r2fr_chl) ASYNC(1)
    !$ACC WAIT(1)

  END SUBROUTINE init_veg_constants

#endif
END MODULE mo_veg_constants
