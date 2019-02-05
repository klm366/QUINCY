!>
!! @brief Contains the science routines for the JSM version in SLM 
!!
!! @author
!!  Sönke Zaehle 
!!
!! @par Revision History
!! First version by Sönke Zaehle (20170124): adopted from JSM_interface and JSM_core
!!
MODULE mo_sb_jsm_processes

#include "dsl_jsbach.inc"

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_jsb_math_constants,  ONLY: dtime, eps8
  USE mo_jsb_lctlib_class,    ONLY: t_lctlib_element 


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_litter_turnover, calc_depolymerisation, & 
            calc_microbial_turnover,calc_microbial_growth, &
            calc_bulk_soil_carbon, calc_bulk_density_correction, calc_qmax_bulk_density_correction, &
            calc_om_sorption, calc_depolymerisation_SESAM, &
            calc_inorganic_uptake, calc_sinking_flux_po4, calc_sourcing_flux_po4, calc_fast_po4

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_jsm_processes'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_jsm
  ! 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate turnover of litter pools according to first order decay 
  !! provides a link between the soluable and woody litter pools introduced in SLM, and the
  !! DOM and polymeric pools of JSM
  !!
  !! Input: temperature and moisture dependence of litter decomposition, litter pool size, and its
  !!        turnover time at standard conditions 
  !! Output: transfer of litter to its sink pool (soluable -> DOM; woody -> polymeric litter) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_litter_turnover(nc, nsoil_sb          , &
                                  tau_litter_pool       , &
                                  fc_litter2sink        , &
                                  rtm_depolymerisation  , &
                                  rmm_depolymerisation  , & 
                                  pool_litter           , &
                                  loss_litter           , &
                                  formation_sink        , & 
                                  het_respiration       , & 
                                  het_respiration_c13   , & 
                                  het_respiration_c14)

    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                         !< dimensions
                                                                       nsoil_sb
    REAL(wp),                                         INTENT(in)    :: tau_litter_pool       !< turnover time of the litter pool [1/s]
    REAL(wp),                                         INTENT(in)    :: fc_litter2sink        !< fraction of litter C entering sink pool (rest is respired) 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_depolymerisation  !< temperature rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_depolymerisation  !< moisture rate modifier [-]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_litter           !< litter pool to decompose [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_litter           !< reduction in litter pool [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_sink        !< formation of SOM pool [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: het_respiration       !< heterotrophic respiration [µmol C/m3/s] 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: het_respiration_c13   !< heterotrophic respiration [µmol 13C/m3/s] 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: het_respiration_c14   !< heterotrophic respiration [pmol 14C/m3/s] 
    ! ---------------------------
    ! 0.2 Local 
    REAL(wp), DIMENSION(nc, nsoil_sb)                 :: hlp1 , &
                                                         decay_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_turnover'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ..


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 calulate decay per timestep
    !!
    decay_rate(:,:) = rtm_depolymerisation * rmm_depolymerisation /tau_litter_pool * dtime

    !> 2.0 calculate source and sink terms
    !! assuming to change in composition or mass of the material to decompose
    !! 
    hlp1(:,:)                       = decay_rate(:,:) * pool_litter%carbon(:,:)
    loss_litter%carbon(:,:)         = loss_litter%carbon(:,:) + hlp1(:,:)
    formation_sink%carbon(:,:)      = formation_sink%carbon(:,:) + fc_litter2sink * hlp1(:,:) 
    het_respiration(:,:)            = het_respiration(:,:) + (1._wp - fc_litter2sink) * hlp1(:,:) * 1000000._wp / dtime 

    hlp1(:,:)                       = decay_rate(:,:) * pool_litter%nitrogen(:,:)
    loss_litter%nitrogen(:,:)       = loss_litter%nitrogen(:,:) + hlp1(:,:)
    formation_sink%nitrogen(:,:)    = formation_sink%nitrogen(:,:) + hlp1(:,:) 

    hlp1(:,:)                       = decay_rate(:,:) * pool_litter%phosphorus(:,:)
    loss_litter%phosphorus(:,:)     = loss_litter%phosphorus(:,:) + hlp1(:,:)
    formation_sink%phosphorus(:,:)  = formation_sink%phosphorus(:,:) + hlp1(:,:) 

    hlp1(:,:)                       = decay_rate(:,:) * pool_litter%carbon13(:,:)
    loss_litter%carbon13(:,:)       = loss_litter%carbon13(:,:) + hlp1(:,:)
    formation_sink%carbon13(:,:)    = formation_sink%carbon13(:,:) + fc_litter2sink * hlp1(:,:) 
    het_respiration_c13(:,:)        = het_respiration_c13(:,:) + (1._wp - fc_litter2sink) * hlp1(:,:) * 1000000._wp / dtime 

    hlp1(:,:)                       = decay_rate(:,:) * pool_litter%carbon14(:,:)
    loss_litter%carbon14(:,:)       = loss_litter%carbon14(:,:) + hlp1(:,:)
    formation_sink%carbon14(:,:)    = formation_sink%carbon14(:,:) + fc_litter2sink * hlp1(:,:) 
    het_respiration_c14(:,:)        = het_respiration_c14(:,:) + (1._wp - fc_litter2sink) * hlp1(:,:) * 1000000._wp / dtime 

    hlp1(:,:)                       = decay_rate(:,:) * pool_litter%nitrogen15(:,:)
    loss_litter%nitrogen15(:,:)     = loss_litter%nitrogen15(:,:) + hlp1(:,:)
    formation_sink%nitrogen15(:,:)  = formation_sink%nitrogen15(:,:) + hlp1(:,:) 

  END SUBROUTINE calc_litter_turnover

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_jsm
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate sinks of the bio-available inorganic nutrient pools
  !!  according to Zhu et al. (2017)
  !!
  !! Input: solute concentration, microbial biomass, fine root biomass, rate modifiers 
  !!
  !! Output: microbial and plant inorganic uptake rates
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_inorganic_uptake(sb_adsorp_scheme, lctlib, solute, microbial_carbon, fine_root_carbon, &
                                             rtm_mic_uptake, rmm_mic_uptake, rtm_hsc, rmm_hsc, & 
                                             rtm_plant_uptake, rmm_plant_uptake, km1_upt, km2_upt, &
                                             vmax_mic_uptake, km_mic_uptake, vmax_plant_uptake, km_plant_uptake, & 
                                             k_mic_uptake, k_plant_uptake, microbial_uptake, plant_uptake)

    USE mo_sb_constants, ONLY : f_plant_carrier_n, f_mic_carrier_n

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CHARACTER(len=*), INTENT(in)    :: sb_adsorp_scheme        !< eca_full; 2: eca_none
    TYPE(t_lctlib_element), INTENT(in)  :: lctlib              !< land-cover-type library - parameter across pft's
    REAL(wp), INTENT(in)            :: solute                  !< concentration of the solvent to be taken up (mol / m3)
    REAL(wp), INTENT(in)            :: microbial_carbon        !< microbial carbon density (mol / m3)
    REAL(wp), INTENT(in)            :: fine_root_carbon        !< fine_root density (mol / m3)
    REAL(wp), INTENT(in)            :: rtm_mic_uptake          !< temperature rate modifier for plant uptake
    REAL(wp), INTENT(in)            :: rmm_mic_uptake          !< moisture rate modifier for plant uptake
    REAL(wp), INTENT(in)            :: rtm_hsc                 !< temperature rate modifier for half_saturation constants
    REAL(wp), INTENT(in)            :: rmm_hsc                 !< moisture rate modifier for half_saturation constants
    REAL(wp), INTENT(in)            :: rtm_plant_uptake        !< temperature rate modifier for plant uptake
    REAL(wp), INTENT(in)            :: rmm_plant_uptake        !< moisture rate modifier for plant uptake
    REAL(wp), INTENT(in)            :: km1_upt                 !< low-affinity parameter for plant uptake [m3/mol]
    REAL(wp), INTENT(in)            :: km2_upt                 !< high-affinity parameter for plant uptake [m3/mol]
    REAL(wp), INTENT(in)            :: vmax_mic_uptake         !< maximium uptake rate of microbial biomass (µmol S / mol C / s)
    REAL(wp), INTENT(in)            :: km_mic_uptake           !< half-saturation constant of microbial nutrient uptake (mol S / m3)
    REAL(wp), INTENT(in)            :: vmax_plant_uptake       !< maximium uptake rate of microbial biomass (µmol S / mol C / s)
    REAL(wp), INTENT(in)            :: km_plant_uptake         !< half-saturation constant of microbial nutrient uptake (mol S / m3)
    REAL(wp), INTENT(in)            :: k_mic_uptake            !< reaction rate of microbial nutrient uptake (s -1)
    REAL(wp), INTENT(in)            :: k_plant_uptake          !< reaction rate of plant nutrient uptake (s -1)
    REAL(wp), INTENT(inout)         :: microbial_uptake        !< microbial uptake flux (µmol / m3 / s)
    REAL(wp), INTENT(out)           :: plant_uptake            !< plant uptake flux (µmol / m3 / s)
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)              :: hlp1,hlp2
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_inorganic_uptake'


    !> define uptake scheme to use, 1. ECA; 2. non-ECA

    !> 1.0 uptake
    !! co-limited rates of N uptake for microbes and plants 

    SELECT CASE (sb_adsorp_scheme)
      CASE("eca_none")
      hlp1 = km1_upt / (rtm_hsc * rmm_hsc)
      hlp2 = km2_upt * rtm_hsc * rmm_hsc
      plant_uptake = fine_root_carbon * rtm_plant_uptake * rmm_plant_uptake * lctlib%vmax_uptake_n * & 
                     solute * ( hlp1 + 1._wp / ( hlp2 + solute )) 
      CASE("eca_part")
        hlp1 = rtm_hsc * rmm_hsc * km_mic_uptake
        hlp2 = rtm_hsc * rmm_hsc * km_plant_uptake
        microbial_uptake     = rtm_mic_uptake * rmm_mic_uptake * vmax_mic_uptake * microbial_carbon * solute / &
                             (hlp1 + solute + hlp1/hlp2 * f_plant_carrier_n * fine_root_carbon + f_mic_carrier_n * microbial_carbon)
        plant_uptake         = rtm_plant_uptake * rmm_plant_uptake * vmax_plant_uptake * fine_root_carbon * solute / &
                             (hlp2 + solute + hlp2/hlp1 * f_mic_carrier_n * microbial_carbon + f_plant_carrier_n * fine_root_carbon)
      CASE("eca_full")
        hlp1 = rtm_hsc * rmm_hsc * km_mic_uptake
        hlp2 = rtm_hsc * rmm_hsc * km_plant_uptake
        microbial_uptake     = rtm_mic_uptake * rmm_mic_uptake * k_mic_uptake * f_mic_carrier_n * microbial_carbon * solute / &
                             (hlp1 + solute + hlp1/hlp2 * f_plant_carrier_n * fine_root_carbon + f_mic_carrier_n * microbial_carbon)
        plant_uptake         = rtm_plant_uptake * rmm_plant_uptake * k_plant_uptake * f_plant_carrier_n * fine_root_carbon * solute / &
                             (hlp2 + solute + hlp2/hlp1 * f_mic_carrier_n * microbial_carbon + f_plant_carrier_n * fine_root_carbon)

    END SELECT
    !> 1.1 limit to possible solution
    !!
    hlp1 = microbial_uptake + plant_uptake
    hlp2 = MAX(solute-eps8,0.0_wp) * 1000000._wp / dtime
    IF(hlp1 > hlp2 .AND. hlp1 > eps8)THEN
       microbial_uptake = microbial_uptake / hlp1 * hlp2
       plant_uptake     = plant_uptake     / hlp1 * hlp2 
    ENDIF

  END SUBROUTINE calc_inorganic_uptake

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from sb_jsm_main
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate the flux rates of P processes
  !! weathering and occlusion according to Wang et al. 2010; 
  !! desorption (slow exchange) rate based on the ANIMO model (Groenendijk et al. 2005)
  !! biocehmical mineralization is calculated using MM-kinetics assuming Mic C is the limiting factor
  !! Input: rate modifiers and source pools
  !!
  !! Output: biochemical P mineralization rate and loss of source pools, gross P mineralization, 
  !!         weathering, occlusion and (slow) desorption
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_sourcing_flux_po4(flag_sb_prescribe_po4, microbial_carbon, & 
                                            dom_carbon, dom_nitrogen, dom_phosphorus, &
                                            rtm_mic_uptake, rmm_mic_uptake, rtm_hsc, rmm_hsc, &
                                            rtm_sorption, rmm_sorption, rtm_desorption, rmm_desorption, &
                                            dom_assoc_po4, residue_assoc_po4, SOM_PC, &
                                            loss_dom_assoc_po4, loss_residue_assoc_po4, & 
                                            solute_po4, solute_nh4, solute_no3, &
                                            assoc_fast, assoc_slow, po4_primary, & 
                                            bulk_dens_corr_sl, root_carbon_sl, vmax_weath_mineral_sl, &
                                            biochem_min_po4, gross_mineralisation_po4, weathering_po4, &
                                            mic_upt_nh4, mic_upt_no3, mic_upt_po4, &
                                            fast_desorp_po4, occlusion_po4, slow_exchange_po4)


    USE mo_sb_constants,        ONLY: vmax_biochem_po4, km_mic_biochem_po4, km_mic_uptake, vmax_mic_uptake, microbial_pue, &
                                      k_desorp_f_po4, k_occlude_po4, k_desorpt_po4, k_adsorpt_po4, km_dip_biochem_po4, &
                                      vmax_mic_uptake_nh4, km_mic_uptake_nh4, vmax_mic_uptake_no3, km_mic_uptake_no3, &
                                      microbial_cue_max, microbial_cn, microbial_np, microbial_pue, microbial_nue, &
                                      km_rootc_weath_po4, km_rootc_biochem_po4, f_mic_carrier_po4, f_plant_carrier_po4, biomin_thres_PC

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    LOGICAL,  INTENT(in)            :: flag_sb_prescribe_po4     !< whether PO4 is prescribed or not
    REAL(wp), INTENT(in)            :: microbial_carbon          !< microbial carbon density (mol / m3)
    REAL(wp), INTENT(in)            :: dom_carbon                !< dom carbon density (mol / m3)
    REAL(wp), INTENT(in)            :: dom_nitrogen              !< dom phosphorus density (mol / m3)
    REAL(wp), INTENT(in)            :: dom_phosphorus            !< dom phosphorus density (mol / m3)
    REAL(wp), INTENT(in)            :: rtm_mic_uptake            !< temperature rate modifier for microbial (and plant) uptake
    REAL(wp), INTENT(in)            :: rmm_mic_uptake            !< moisture rate modifier for microbial (and plant) uptake
    REAL(wp), INTENT(in)            :: rtm_hsc                   !< temperature rate modifier for half_saturation constants
    REAL(wp), INTENT(in)            :: rmm_hsc                   !< moisture rate modifier for half_saturation constants
    REAL(wp), INTENT(in)            :: rtm_sorption              !< temperature rate modifier for sorption
    REAL(wp), INTENT(in)            :: rmm_sorption              !< moisture rate modifier for sorption
    REAL(wp), INTENT(in)            :: rtm_desorption            !< temperature rate modifier for desorption
    REAL(wp), INTENT(in)            :: rmm_desorption            !< moisture rate modifier for desorption
    REAL(wp), INTENT(in)            :: dom_assoc_po4             !< PO4 in the mineral associated dom pool [mol/m3]
    REAL(wp), INTENT(in)            :: residue_assoc_po4         !< PO4 in the mineral associated microbial necromass pool [mol/m3]
    REAL(wp), INTENT(in)            :: SOM_PC                    !< P:C ratio in the slow SOM pool [mol/mol]
    REAL(wp), INTENT(in)            :: solute_po4                !< soluble po4 pool [mol/m3]
    REAL(wp), INTENT(in)            :: solute_nh4                !< soluble nh4 pool [mol/m3]
    REAL(wp), INTENT(in)            :: solute_no3                !< soluble no3 pool [mol/m3]
    REAL(wp), INTENT(in)            :: assoc_slow                !< slow minerally associated PO4 pool [mol/m3]
    REAL(wp), INTENT(in)            :: assoc_fast                !< fast minerally associated PO4 pool [mol/m3]
    REAL(wp), INTENT(in)            :: po4_primary               !< primary PO4 pool [mol/m3]
    REAL(wp), INTENT(in)            :: bulk_dens_corr_sl         !< OM corrected soil bulk density [kg/mol]
    REAL(wp), INTENT(in)            :: root_carbon_sl            !< carbon of fine roots per soil layer, [mol/m3]
    REAL(wp), INTENT(in)            :: vmax_weath_mineral_sl     !< weathering coefficient of mineral soil [mol P /m3 /s]
    REAL(wp), INTENT(inout)         :: loss_dom_assoc_po4        !< loss of PO4 in the mineral associated dom pool [mol/m3]
    REAL(wp), INTENT(inout)         :: loss_residue_assoc_po4    !< loss of PO4 in the mineral associated microbial necromass pool [mol/m3]
    REAL(wp), INTENT(out)           :: biochem_min_po4           !< phosphorus biochemical mineralization rate, [µmol P /m3 /s]
    REAL(wp), INTENT(out)           :: gross_mineralisation_po4  !< phosphorus gross(biological) mineralization rate, [µmol P /m3 /s]
    REAL(wp), INTENT(out)           :: weathering_po4            !< phosphorus weathering rate [µmol P / m3 / s]
    REAL(wp), INTENT(out)           :: fast_desorp_po4           !< phosphorus fast desorption rate [µmol P / m3 / s]
    REAL(wp), INTENT(out)           :: occlusion_po4             !< phosphorus occlusion rate, assoc_slow --> occluded [µmol P / m3 / s]
    REAL(wp), INTENT(out)           :: slow_exchange_po4         !< phosphorus desorption rate, assoc_fast --> assoc_slow [µmol P / m3 / s]
    REAL(wp), INTENT(out)           :: mic_upt_nh4               !< microbial nh4 uptake rate [µmol N / m3 / s], for eca_none only
    REAL(wp), INTENT(out)           :: mic_upt_no3               !< microbial no3 uptake rate [µmol N / m3 / s], for eca_none only
    REAL(wp), INTENT(out)           :: mic_upt_po4               !< microbial po4 uptake rate [µmol P / m3 / s], for eca_none only
    ! ---------------------------    
    ! 0.2 Local
    REAL(wp)                        :: hlp1, hlp2, hlp3, hlp4, hlp5, uptake_c, f_uptake, uptake_n, uptake_p, mic_upt_n
    CHARACTER(len=*), PARAMETER     :: routine = TRIM(modname)//':calc_sourcing_flux_po4'
    ! ---------------------------


    !> 1.0 Calculates the biochemical P mineralization, and the loss of the contributing pools.
    !! P is assumed to clevage out from poly litter, dom associated, residue and residue associated 
    IF(.NOT.flag_sb_prescribe_po4)THEN
       hlp1 = (dom_assoc_po4 + residue_assoc_po4) !polymeric_litter_po4 + residue_po4 + 
       IF (hlp1 > eps8) THEN
          hlp2 = rtm_hsc * rmm_hsc * km_dip_biochem_po4 / (rtm_hsc * rmm_hsc * km_dip_biochem_po4 + solute_po4)
          hlp3 = (root_carbon_sl * f_plant_carrier_po4 + microbial_carbon * f_mic_carrier_po4) 
          hlp5 =  hlp3 / (rtm_hsc * rmm_hsc * km_mic_biochem_po4 + hlp3)
          IF (SOM_PC < biomin_thres_PC) THEN
              hlp4 = 0.0_wp
          ELSE 
              hlp4 = SOM_PC / ( 0.005_wp + SOM_PC)
          ENDIF 
          ! hlp5 = microbial_carbon / (rtm_hsc * rmm_hsc * km_mic_biochem_po4 + microbial_carbon)
          ! instantaneous rate of biochemical mineralisation 1/s
          biochem_min_po4 = rtm_mic_uptake * rmm_mic_uptake * vmax_biochem_po4 * hlp2 * hlp4* hlp5
 
          ! loss terms mol/m3/timestep
          ! loss_polymeric_litter_po4 = loss_polymeric_litter_po4 + biochem_min_po4 * polymeric_litter_po4 * dtime
          loss_dom_assoc_po4 = loss_dom_assoc_po4 + biochem_min_po4 * dom_assoc_po4 * dtime
          ! loss_residue_po4 = loss_residue_po4 + biochem_min_po4 * residue_po4 * dtime
          loss_residue_assoc_po4 = loss_residue_assoc_po4 + biochem_min_po4 * residue_assoc_po4 * dtime
          ! biochemical mineralisation in µmol P / m3 / s
          biochem_min_po4 = biochem_min_po4 * hlp1 * 1.e6_wp
       ELSE
          biochem_min_po4 = 0.0_wp
       ENDIF
    ELSE
       biochem_min_po4 = 0.0_wp
    ENDIF
    !> 2.0 calculate the gross P mineralization, based on the potential depolymerization of DOM
    !!
    hlp1 = dom_carbon / (km_mic_uptake * rtm_hsc * rmm_hsc + dom_carbon )
    uptake_c = vmax_mic_uptake * rtm_mic_uptake * rmm_mic_uptake * hlp1 * microbial_carbon
    IF(dom_carbon > eps8) THEN
       f_uptake = MIN(uptake_c / dom_carbon, 1._wp * 1.e6_wp/dtime)
    ELSE
       f_uptake = 0.0_wp
    ENDIF
    uptake_p = f_uptake * dom_phosphorus
    uptake_n = f_uptake * dom_nitrogen
    uptake_c = f_uptake * dom_carbon
    gross_mineralisation_po4 = (1._wp - microbial_pue) * uptake_p
    !> 2.1 calculate the potential microbial N, P uptake
    mic_upt_po4 = MAX(0._wp, uptake_c * microbial_cue_max / microbial_cn / microbial_np - uptake_p * microbial_pue)
    mic_upt_n   = MAX(0._wp, uptake_c * microbial_cue_max / microbial_cn - uptake_n * microbial_nue)
    hlp1        = vmax_mic_uptake_nh4 * microbial_carbon * solute_nh4 / (solute_nh4 + km_mic_uptake_nh4)
    hlp2        = vmax_mic_uptake_no3 * microbial_carbon * solute_no3 / (solute_no3 + km_mic_uptake_no3)
    mic_upt_nh4 = mic_upt_n * hlp1 /(hlp1 + hlp2)
    mic_upt_no3 = mic_upt_n * hlp2 /(hlp1 + hlp2)

    !> 3.0 calculate the fast desorption rate, based on the Langmuir kinetics
    fast_desorp_po4 = k_desorp_f_po4 * assoc_fast * 1.e6_wp

    !> 4.0 Pi fluxes between different mineral P pools
    !! weathering, occlusion and slow exchange between sorbed and labile P pools
    weathering_po4 = rtm_hsc / rmm_hsc * vmax_weath_mineral_sl * bulk_dens_corr_sl * root_carbon_sl / & 
                          (km_rootc_weath_po4 + root_carbon_sl) * 1.e6_wp
    IF (weathering_po4 /1.e6_wp > po4_primary / dtime) weathering_po4 = po4_primary / dtime * 1.e6_wp
    occlusion_po4 = assoc_slow * k_occlude_po4 * 1.e6_wp

    hlp1 = rtm_desorption * k_desorpt_po4 * assoc_slow 
    hlp2 = rtm_sorption * k_adsorpt_po4 * assoc_fast     
    slow_exchange_po4 = (hlp2 - hlp1) * 1.e6_wp

  END SUBROUTINE calc_sourcing_flux_po4

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate sinks of the bio-available inorganic po4 pools
  !! ECA approach modified based on Zhu et al. (2016) and Zhu et al. (2017), fast equilibrium based on Van der Zee et al. 1989
  !!
  !! Input: Pi solute concentration, microbial biomass, fine root biomass, sorption capacity, soil Pi pools, rate modifiers
  !!
  !! Output: microbial and plant inorganic uptake rates, fast adsorption and desorption rates
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_sinking_flux_po4(sb_adsorp_scheme, flag_sb_ECA_upt_po4, lctlib, &
                                             microbial_carbon, fine_root_carbon, solute, &
                                             rtm_mic_uptake, rmm_mic_uptake, rtm_hsc, rmm_hsc, &
                                             rtm_plant_uptake, rmm_plant_uptake, km1_upt, km2_upt, &
                                             qmax_po4, km_adsorpt_po4_sl, assoc_fast, microbial_uptake, plant_uptake, &
                                             fast_desorp_po4, fast_adsorpt_po4, solute_ava, &
                                             eca_divisor_adsorp, eca_divisor_mic, eca_divisor_plant)

    USE mo_sb_constants,    ONLY: f_mic_carrier_po4, km_mic_uptake_po4, km_plant_uptake_po4, f_plant_carrier_po4, &
                                  k_mic_uptake_po4, k_plant_uptake_po4, vmax_mic_uptake_po4, vmax_plant_uptake_po4, k_adsorpt_f_po4, k_desorp_f_po4

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CHARACTER(len=*), INTENT(in)    :: sb_adsorp_scheme          !< eca_full; eca_none; eca_part
    LOGICAL,          INTENT(in)    :: flag_sb_ECA_upt_po4       !< T: eca only applies to uptake; F: eca also applies to adsorption
    TYPE(t_lctlib_element), INTENT(in)  :: lctlib                !< land-cover-type library - parameter across pft's
    REAL(wp)        , INTENT(in)    :: microbial_carbon          !< microbial carbon density (mol / m3)
    REAL(wp)        , INTENT(in)    :: fine_root_carbon          !< fine_root density (mol / m3)
    REAL(wp)        , INTENT(in)    :: rtm_mic_uptake            !< temperature rate modifier for microbial uptake
    REAL(wp)        , INTENT(in)    :: rmm_mic_uptake            !< moisture rate modifier for microbial uptake
    REAL(wp)        , INTENT(in)    :: rtm_hsc                   !< temperature rate modifier for half_saturation constants
    REAL(wp)        , INTENT(in)    :: rmm_hsc                   !< moisture rate modifier for half_saturation constants
    REAL(wp)        , INTENT(in)    :: rtm_plant_uptake          !< temperature rate modifier for plant uptake
    REAL(wp)        , INTENT(in)    :: rmm_plant_uptake          !< moisture rate modifier for plant uptake
    REAL(wp)        , INTENT(in)    :: km1_upt                   !< low-affinity parameter for plant uptake [m3/mol]
    REAL(wp)        , INTENT(in)    :: km2_upt                   !< high-affinity parameter for plant uptake [m3/mol]
    REAL(wp)        , INTENT(in)    :: solute                    !< concentration of the solvent to be taken up (mol / m3)
    REAL(wp)        , INTENT(in)    :: qmax_po4                  !< sorption capacity of po4 to mineral soil [mol P / m3]
    REAL(wp)        , INTENT(in)    :: km_adsorpt_po4_sl         !< half-saturation of po4 sorption to mineral soil [mol P / m3]
    REAL(wp)        , INTENT(in)    :: assoc_fast                !< fast minerally associated PO4 pool [mol/m3]
    REAL(wp)        , INTENT(inout) :: microbial_uptake          !< microbial uptake flux (µmol / m3 / s)
    REAL(wp)        , INTENT(out)   :: plant_uptake              !< plant uptake flux (µmol / m3 / s)
    REAL(wp)        , INTENT(inout) :: solute_ava                !< concentration of the available solvent to be taken up per second(µmol / m3 / s)
    REAL(wp)        , INTENT(inout) :: fast_desorp_po4           !< fast desorption flux (µmol / m3 / s)
    REAL(wp)        , INTENT(out)   :: fast_adsorpt_po4          !< fast adsorption flux (µmol / m3 / s)
    REAL(wp)        , INTENT(inout) :: eca_divisor_adsorp, &     !< denominator of ECA item for po4 adsorption [mol P / m3]
                                       eca_divisor_mic, &
                                       eca_divisor_plant

    ! 0.2 Local
    REAL(wp)                                             :: hlp1, hlp2, hlp3, hlp4, hlp5, avail_sorption_sites
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_sinking_flux_po4'
    ! ---------------------------


    !> define uptake scheme to use, 1. ECA; 2. non-ECA
    SELECT CASE (sb_adsorp_scheme)
      CASE("eca_full", "eca_part")
    !> 1.0 uptake
    !! co-limited rates of P uptake for microbes and plants

      !> 1.1 ECA approach
      !! 
        avail_sorption_sites = MAX(0._wp, qmax_po4 - assoc_fast) 
        hlp5 = solute + assoc_fast
        IF (flag_sb_ECA_upt_po4) THEN
          fast_desorp_po4  = 0._wp
          fast_adsorpt_po4 = 0._wp
          hlp1 = km_mic_uptake_po4 + hlp5 + f_mic_carrier_po4 * microbial_carbon + &
                 km_mic_uptake_po4 / km_plant_uptake_po4 * f_plant_carrier_po4 * fine_root_carbon
          hlp2 = km_plant_uptake_po4 + hlp5 + f_plant_carrier_po4 * fine_root_carbon + &
                 km_plant_uptake_po4 / km_mic_uptake_po4 * f_mic_carrier_po4 * microbial_carbon
        ELSE
          hlp1 = km_mic_uptake_po4 + hlp5 + f_mic_carrier_po4 * microbial_carbon + &
                 km_mic_uptake_po4 / km_plant_uptake_po4 * f_plant_carrier_po4 * fine_root_carbon + &
                 km_mic_uptake_po4 / km_adsorpt_po4_sl * avail_sorption_sites
          hlp2 = km_plant_uptake_po4 + hlp5 + f_plant_carrier_po4 * fine_root_carbon + &
                 km_plant_uptake_po4 / km_mic_uptake_po4 * f_mic_carrier_po4 * microbial_carbon + &
                 km_plant_uptake_po4 / km_adsorpt_po4_sl * avail_sorption_sites

          hlp3 = km_adsorpt_po4_sl + hlp5 + avail_sorption_sites + &
                 km_adsorpt_po4_sl / km_plant_uptake_po4 * f_plant_carrier_po4 * fine_root_carbon + &
                 km_adsorpt_po4_sl / km_mic_uptake_po4 * f_mic_carrier_po4 * microbial_carbon
          hlp4 = k_desorp_f_po4 / km_adsorpt_po4_sl * solute !k_adsorpt_f_po4 * solute
          fast_adsorpt_po4 = MIN(1._wp / dtime, hlp4 * hlp5 / hlp3) * avail_sorption_sites * 1.e6_wp
          eca_divisor_adsorp = hlp3
          eca_divisor_mic    = hlp1
          eca_divisor_plant  = hlp2
        ENDIF
        microbial_uptake = MIN(rtm_mic_uptake * rmm_mic_uptake * k_mic_uptake_po4 * hlp5 * f_mic_carrier_po4 * &
                            microbial_carbon / hlp1, microbial_uptake)
        plant_uptake     = rtm_plant_uptake * rmm_plant_uptake * k_plant_uptake_po4 * hlp5 * f_plant_carrier_po4 * &
                            fine_root_carbon / hlp2 

        IF (sb_adsorp_scheme=="eca_part") THEN
          microbial_uptake = rtm_mic_uptake * rmm_mic_uptake * vmax_mic_uptake_po4 * hlp5 * &
                              microbial_carbon / hlp1!, microbial_uptake)
          plant_uptake     = rtm_plant_uptake * rmm_plant_uptake * vmax_plant_uptake_po4 * hlp5 * &
                              fine_root_carbon / hlp2
        ENDIF
        
      !> 1.2 non-ECA approach
      !! 
      CASE("eca_none")
        plant_uptake       = rtm_plant_uptake * rmm_plant_uptake * fine_root_carbon * lctlib%vmax_uptake_p * &
                             solute * (km1_upt + 1._wp / (km2_upt + solute)) 
        fast_desorp_po4    = 0._wp
        fast_adsorpt_po4   = 0._wp       
    END SELECT

    !! check the availability of soluable po4 considering all sourcing fluxes
    hlp1 = microbial_uptake + plant_uptake + fast_adsorpt_po4
    solute_ava = solute_ava !+ fast_desorp_po4
    IF ((hlp1 > solute_ava)  .AND. hlp1 > 0._wp)THEN
        microbial_uptake = microbial_uptake * solute_ava / hlp1
        plant_uptake     = plant_uptake     * solute_ava / hlp1
        fast_adsorpt_po4 = fast_adsorpt_po4 * solute_ava / hlp1
    ENDIF

  END SUBROUTINE calc_sinking_flux_po4


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate turnover of the microbial pool
  !! corresponds to the term pi * Cb in Ahrens et al. 2015, eq. 7
  !!
  !! Input: microbial biomass
  !!
  !! Output: transfer of microbial biomass to its sink pools (soluable -> DOM; polymeric -> residue)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_microbial_turnover(nc, nsoil_sb, pool_microbial,loss_microbial,formation_dom,formation_residue)

    USE mo_sb_constants,          ONLY: k_mort_microbial, min_mic_biomass, f_soluable_residue, fc_p_res2dom, fc_n_res2dom
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                        INTENT(in)    :: nc, &              !< dimensions
                                                     nsoil_sb
    TYPE(t_sb_element_a),           INTENT(in)    :: pool_microbial     !< microbial biomass [mol/m3]
    TYPE(t_sb_element_a),           INTENT(inout) :: loss_microbial     !< microbial death [mol/m3/timestep]
    TYPE(t_sb_element_a),           INTENT(inout) :: formation_dom      !< DOM formation from microbial death [mol/m3/timestep]
    TYPE(t_sb_element_a),           INTENT(inout) :: formation_residue  !< residue formation from microbial death [mol/m3/timestep]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1 , &
                                                        mortality_rate
    INTEGER                                          :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_microbial_turnover'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 calculate microbial mortality rate, limited to minimum microbial biomass density (-> dormant microbes)
    !!
    mortality_rate(:,:) = k_mort_microbial * dtime
    hlp1(:,:)           = mortality_rate(:,:) * pool_microbial%carbon(:,:) 
    WHERE((pool_microbial%carbon(:,:) - hlp1(:,:)) < min_mic_biomass)
      mortality_rate(:,:) = MAX(0._wp, 1._wp - (min_mic_biomass / pool_microbial%carbon(:,:)))
    ENDWHERE

    !> 2.0 calculate matter fluxes associated with microbial death
    !! a fraction turns to DOM, the remainder becomes residue
    !! the assumption is made that the composition of this material is not changed
    !!@todo: define C:N:P stoichiometries of DOM and residues from microbial death 
    !!
    hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%carbon(:,:)
    loss_microbial%carbon(:,:)        = loss_microbial%carbon(:,:) + hlp1(:,:)
    formation_dom%carbon(:,:)         = formation_dom%carbon(:,:) + f_soluable_residue * hlp1(:,:) 
    formation_residue%carbon(:,:)     = formation_residue%carbon(:,:) + (1._wp - f_soluable_residue) * hlp1(:,:) 

    hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%nitrogen(:,:)
    loss_microbial%nitrogen(:,:)      = loss_microbial%nitrogen(:,:) + hlp1(:,:)
    formation_dom%nitrogen(:,:)       = formation_dom%nitrogen(:,:) + f_soluable_residue * hlp1(:,:) + &
                                        (1._wp - f_soluable_residue) * hlp1(:,:) * fc_n_res2dom
    formation_residue%nitrogen(:,:)   = formation_residue%nitrogen(:,:) + (1._wp - f_soluable_residue) * (1._wp - fc_n_res2dom) * hlp1(:,:) 

    hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%phosphorus(:,:)
    loss_microbial%phosphorus(:,:)    = loss_microbial%phosphorus(:,:) + hlp1(:,:)
    formation_dom%phosphorus(:,:)     = formation_dom%phosphorus(:,:) + f_soluable_residue * hlp1(:,:) + &
                                        (1._wp - f_soluable_residue) * hlp1(:,:) * fc_p_res2dom
    formation_residue%phosphorus(:,:) = formation_residue%phosphorus(:,:) + (1._wp - f_soluable_residue) * (1._wp - fc_p_res2dom) * hlp1(:,:) 

    hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%carbon13(:,:)
    loss_microbial%carbon13(:,:)      = loss_microbial%carbon13(:,:) + hlp1(:,:)
    formation_dom%carbon13(:,:)       = formation_dom%carbon13(:,:) + f_soluable_residue * hlp1(:,:) 
    formation_residue%carbon13(:,:)   = formation_residue%carbon13(:,:) + (1._wp - f_soluable_residue) * hlp1(:,:) 

    hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%carbon14(:,:)
    loss_microbial%carbon14(:,:)      = loss_microbial%carbon14(:,:) + hlp1(:,:)
    formation_dom%carbon14(:,:)       = formation_dom%carbon14(:,:) + f_soluable_residue * hlp1(:,:) 
    formation_residue%carbon14(:,:)   = formation_residue%carbon14(:,:) + (1._wp - f_soluable_residue) * hlp1(:,:) 

    hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%nitrogen15(:,:)
    loss_microbial%nitrogen15(:,:)    = loss_microbial%nitrogen15(:,:) + hlp1(:,:)
    formation_dom%nitrogen15(:,:)     = formation_dom%nitrogen15(:,:) + f_soluable_residue * hlp1(:,:) + &
                                        (1._wp - f_soluable_residue) * hlp1(:,:) * fc_n_res2dom
    formation_residue%nitrogen15(:,:) = formation_residue%nitrogen15(:,:) + (1._wp - f_soluable_residue) * (1._wp - fc_n_res2dom) * hlp1(:,:) 


  END SUBROUTINE calc_microbial_turnover


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_jsm
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate growth of the microbial pool 
  !! corresponds to Ahrens et al. 2015, eq. 2 
  !! extended to account for N and P limitation of microbial growth given a prescribed 
  !! stoichometry
  !!
  !! Input: DOM and microbial pools, microbial nutrient uptake rate, temperature and moisture rate modifiers 
  !!
  !! Output: DOM consumption, microbial growth, hetertrophic respiration and net mineralisation
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_microbial_growth(nc, nsoil_sb, sb_bnf_scheme,pool_dom,pool_microbial, & 
                                   rtm_mic_uptake,rmm_mic_uptake,rtm_hsc,rmm_hsc,microbial_cue_eff_mavg, &
                                   microbial_uptake_nh4_sl,microbial_uptake_nh4_n15_sl, &
                                   microbial_uptake_no3_sl,microbial_uptake_no3_n15_sl,microbial_uptake_po4_sl, &
                                   loss_dom,formation_microbial, &
                                   het_respiration,het_respiration_c13,het_respiration_c14, &
                                   net_mineralisation_nh4,net_mineralisation_nh4_n15, &
                                   net_mineralisation_no3,net_mineralisation_no3_n15,net_mineralisation_po4, &
                                   microbial_cue_eff,asymb_n_fixation)

    USE mo_isotope_util,          ONLY: calc_fractionation
    USE mo_sb_constants,          ONLY: vmax_mic_uptake,km_mic_uptake,microbial_cue,microbial_nue,microbial_pue, & 
                                        microbial_cn,microbial_np,eta_ammonification,microbial_cue_min,microbial_cue_max
    USE mo_jsb_math_constants,    ONLY: dtime
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &               !< dimensions
                                                                       nsoil_sb
    CHARACTER(len=*),                                 INTENT(in)    :: sb_bnf_scheme       !< veg n fixation scheme  
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_dom            !< DOM [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_microbial      !< microbial biomass [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_mic_uptake      !< temperature rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_mic_uptake      !< moisture rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_hsc             !< temperature rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_hsc             !< moisture rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: microbial_cue_eff_mavg       !< microbial effective CUE given NP constraints
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: microbial_uptake_nh4_sl      !< microbial uptake rate [µmol N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: microbial_uptake_nh4_n15_sl  !< microbial uptake rate [µmol 15N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: microbial_uptake_no3_sl      !< microbial uptake rate [µmol N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: microbial_uptake_no3_n15_sl  !< microbial uptake rate [µmol 15N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: microbial_uptake_po4_sl      !< microbial uptake rate [µmol P/m3/s]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_dom            !< microbial consumption of DOM [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_microbial !< microbial growth on DOM [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: het_respiration     !< respiration rate [µmol C/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: het_respiration_c13 !< respiration rate [µmol 13C/m3/s] 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: het_respiration_c14 !< respiration rate [µmol 14C/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: net_mineralisation_nh4     !< net mineralisation rate [µmol N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: net_mineralisation_nh4_n15 !< net mineralisation rate [µmol 15N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: net_mineralisation_no3     !< net mineralisation rate [µmol N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: net_mineralisation_no3_n15 !< net mineralisation rate [µmol 15N/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: net_mineralisation_po4     !< net mineralisation rate [µmol P/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: microbial_cue_eff          !< current effective CUE given NP constraints 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: asymb_n_fixation           !< current effective CUE given NP constraints 
    ! ---------------------------
    ! 0.2 Local 
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1             , &
                                                        f_uptake         , &
                                                        uptake_dom  
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: growth_c,&
                                                        growth_n,&
                                                        growth_p,&
                                                        growth,&
                                                        uptake_c,&
                                                        uptake_n,&
                                                        uptake_p 
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: scal
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: microbial_nue_eff,&
                                                        microbial_pue_eff
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: gross_mineralisation_nh4,&
                                                        gross_mineralisation_nh4_n15,&
                                                        gross_mineralisation_po4 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_microbial_growth'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...
    
    
    

    !> 1.0 potential microbial DOM consumption (µmol / m3 / s) 
    !! 
    hlp1(:,:)     = pool_dom%carbon(:,:) / (km_mic_uptake * rtm_hsc(:,:) * rmm_hsc(:,:) + pool_dom%carbon(:,:) )
    uptake_c(:,:) = vmax_mic_uptake * rtm_mic_uptake(:,:) * rmm_mic_uptake(:,:) * hlp1(:,:) * pool_microbial%carbon(:,:) 
    WHERE(pool_dom%carbon(:,:) > eps8)
       f_uptake(:,:) = MIN(uptake_c(:,:) / pool_dom%carbon(:,:), 1._wp * 1.e6_wp/dtime)
    ELSEWHERE
       f_uptake(:,:) = 0.0_wp
    ENDWHERE
    uptake_n(:,:) = f_uptake(:,:) * pool_dom%nitrogen(:,:) 
    uptake_p(:,:) = f_uptake(:,:) * pool_dom%phosphorus(:,:)

    !> 1.1 microbial growth rates in µmol C / m3 / s, given by C, N and P availability, and the rate-limiting growth rate
    !! 
    growth_c(:,:) = microbial_cue_max * uptake_c(:,:) 
    growth_n(:,:) = (microbial_nue * uptake_n(:,:) + &
                     microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:)) * microbial_cn
    growth_p(:,:) = (microbial_pue * uptake_p(:,:) + &
                      microbial_uptake_po4_sl(:,:)) * microbial_np * microbial_cn
    WHERE(growth_c(:,:) <= growth_n(:,:) .AND. growth_c(:,:) <= growth_p(:,:))
       ! carbon uptake is limiting
       growth(:,:) = growth_c(:,:)
    ELSEWHERE(growth_n(:,:) < growth_c(:,:) .AND. growth_n(:,:) <= growth_p(:,:))
       ! nitrogen uptake is limiting
       growth(:,:) = growth_n(:,:)
    ELSEWHERE(growth_p(:,:) < growth_c(:,:) .AND. growth_p(:,:) < growth_n(:,:))
       ! phosphorus uptake is limiting
       growth(:,:) = growth_p(:,:)
    ENDWHERE

    !> 1.2 Adjust carbon-use efficiency such that the growth is stoichiometrically balanced
    !!
    microbial_cue_eff(:,:) = microbial_cue_max
    WHERE(growth(:,:) < growth_c(:,:) .AND. uptake_c(:,:) > eps8)
       ! excess carbon is going to be respired
       microbial_cue_eff(:,:) = MAX(microbial_cue_min,growth(:,:) / uptake_c(:,:))
    ELSEWHERE(growth(:,:) < growth_c(:,:))
       microbial_cue_eff(:,:) = microbial_cue_min
    ENDWHERE

    !> 2.0 actual growth given actual CUE, adjust uptake rate of DOM to the most 
    !!     limiting nutrient
    !!
    !! growth_n and growth_p are unchanged, only growth_c needs to be adjusted to current 
    !! CUE
    growth_c(:,:) = microbial_cue_eff_mavg(:,:) * uptake_c(:,:) 
    IF(sb_bnf_scheme=='unlimited')THEN
       asymb_n_fixation(:,:) = MAX(0.0_wp, growth_c(:,:) / microbial_cn - &
                                           microbial_nue * f_uptake(:,:) * pool_dom%nitrogen(:,:) - &
                                           microbial_uptake_nh4_sl(:,:) - microbial_uptake_no3_sl(:,:))
       microbial_uptake_nh4_sl(:,:) = microbial_uptake_nh4_sl(:,:) + asymb_n_fixation(:,:)
    ELSE
       WHERE(growth_c(:,:) <= growth_n(:,:) .AND. growth_c(:,:) <= growth_p(:,:))
          ! carbon uptake is limiting, full uptake rate is realised
          f_uptake(:,:) = f_uptake(:,:) 
       ELSEWHERE(growth_n(:,:) < growth_c(:,:) .AND. growth_n(:,:) <= growth_p(:,:).AND.pool_dom%carbon(:,:) > eps8)
          ! nitrogen uptake is limiting
          ! calculate uptake rate that satisfies N constraint and CUE
          f_uptake(:,:) = MAX(0.0_wp,MIN(f_uptake(:,:), & 
                                         ((microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:)) * microbial_cn) / & 
                                         (microbial_cue_eff_mavg(:,:) * pool_dom%carbon(:,:) - & 
                                          microbial_nue * pool_dom%nitrogen(:,:) * microbial_cn)))
       ELSEWHERE(growth_p(:,:) < growth_c(:,:) .AND. growth_p(:,:) < growth_n(:,:) .AND. pool_dom%carbon(:,:) > eps8)
          ! phosphorus uptake is limiting
          ! calculate uptake rate that satisfies P constraint and CUE
          f_uptake(:,:) = MAX(0.0_wp,MIN(f_uptake(:,:), &
                                         (microbial_uptake_po4_sl(:,:) * microbial_np * microbial_cn) / & 
                                         (microbial_cue_eff_mavg(:,:) * pool_dom%carbon(:,:) - & 
                                          microbial_pue * pool_dom%phosphorus(:,:) * microbial_np * microbial_cn)))
       ELSEWHERE
          f_uptake(:,:) = 0.0_wp
       ENDWHERE
   
       growth_c(:,:) = microbial_cue_eff_mavg(:,:) * f_uptake(:,:) * pool_dom%carbon(:,:) 
       uptake_n(:,:) = f_uptake(:,:) * pool_dom%nitrogen(:,:) 
       uptake_p(:,:) = f_uptake(:,:) * pool_dom%phosphorus(:,:)
       growth_n(:,:) = (microbial_nue * uptake_n(:,:) + &
                        microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:)) * microbial_cn
       growth_p(:,:) = (microbial_pue * uptake_p(:,:) + &
                         microbial_uptake_po4_sl(:,:)) * microbial_np * microbial_cn
   
    ENDIF
  
    !> 2.1 Adjust nitrogen-use efficiency and immobilisation flux such that the growth is stoichiometrically balanced
    !!     and calculate gross mineralisation rate (and its isotopic signature)
    microbial_nue_eff(:,:) = microbial_nue
    WHERE(growth_c(:,:) < growth_n(:,:)) 
       ! inorganic N immobilisation flux in agreement with C growth
       hlp1(:,:) = MAX(growth_c(:,:) / microbial_cn - microbial_nue * uptake_n(:,:),0.0_wp)
       WHERE((microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:)) > eps8)
          scal(:,:) = hlp1(:,:) / (microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:)) 
       ELSEWHERE
          scal(:,:) = 0.0_wp
       ENDWHERE
       microbial_uptake_nh4_sl(:,:)     = microbial_uptake_nh4_sl(:,:)     * scal(:,:) 
       microbial_uptake_nh4_n15_sl(:,:) = microbial_uptake_nh4_n15_sl(:,:) * scal(:,:)
       microbial_uptake_no3_sl(:,:)     = microbial_uptake_no3_sl(:,:)     * scal(:,:) 
       microbial_uptake_no3_n15_sl(:,:) = microbial_uptake_no3_n15_sl(:,:) * scal(:,:)
       WHERE(uptake_n(:,:) > eps8)
          microbial_nue_eff(:,:) = MIN(microbial_nue, &
                                      ( growth_c(:,:) / microbial_cn - microbial_uptake_nh4_sl(:,:) - microbial_uptake_no3_sl(:,:) ) / uptake_n(:,:))
       ENDWHERE
    ENDWHERE
    gross_mineralisation_nh4(:,:)     = (1._wp - microbial_nue_eff(:,:)) * uptake_n(:,:) 
    gross_mineralisation_nh4_n15(:,:) = gross_mineralisation_nh4(:,:) * calc_fractionation(pool_dom%nitrogen(:,:), & 
                                                                                           pool_dom%nitrogen15(:,:),& 
                                                                                           eta_ammonification)

    !> 2.2 Adjust phosphorus-use efficiency and immobilisation flux such that the growth is stoichiometrically balanced
    !!     and calculate gross mineralisation rate 
    microbial_pue_eff(:,:) = microbial_pue
    WHERE(growth_c(:,:) < growth_p(:,:))  
       ! inorganic P immobilisation flux, in agreement with C and N uptake
       microbial_uptake_po4_sl(:,:) = MAX(growth_c(:,:) / microbial_cn / microbial_np - microbial_pue * uptake_p(:,:),0.0_wp)
       ! In the case that the source C:N:P is too tight for microbes, cause gross mineralisation
       WHERE(uptake_p(:,:) > eps8)
          microbial_pue_eff(:,:) = ( growth_c(:,:) / microbial_cn / microbial_np - microbial_uptake_po4_sl(:,:) ) / uptake_p(:,:)
       ENDWHERE
    ENDWHERE
    gross_mineralisation_po4(:,:) = (1._wp - microbial_pue_eff(:,:)) * uptake_p(:,:)

    !> 3.0 update fluxes (note that loss and formation of organic matter is in mol/m3/timestep, whereas fluxes are in µmol/m3/s)
    !! carbon
    uptake_dom(:,:)                 = f_uptake(:,:) * pool_dom%carbon(:,:)
    loss_dom%carbon(:,:)            = loss_dom%carbon(:,:) + uptake_dom(:,:) * dtime / 1000000._wp
    formation_microbial%carbon(:,:) = formation_microbial%carbon(:,:) + & 
                                        microbial_cue_eff_mavg(:,:) * uptake_dom(:,:) * dtime / 1000000._wp 
    het_respiration(:,:)            = het_respiration(:,:) + & 
                                        (1._wp - microbial_cue_eff_mavg(:,:)) * uptake_dom(:,:) 

    !! nitrogen
    net_mineralisation_nh4(:,:)       = gross_mineralisation_nh4(:,:) - microbial_uptake_nh4_sl(:,:)
    net_mineralisation_no3(:,:)       = - microbial_uptake_no3_sl(:,:)
    uptake_dom(:,:)                   = f_uptake(:,:) * pool_dom%nitrogen(:,:)
    loss_dom%nitrogen(:,:)            = loss_dom%nitrogen(:,:) + uptake_dom(:,:) * dtime / 1000000.0_wp 
    formation_microbial%nitrogen(:,:) = formation_microbial%nitrogen(:,:) + & 
                                          (microbial_nue_eff(:,:) * uptake_dom(:,:) + & 
                                          microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:)) * dtime / 1000000.0_wp 
   
    !! phosphorus
    net_mineralisation_po4(:,:)         = gross_mineralisation_po4(:,:) - microbial_uptake_po4_sl(:,:)
    uptake_dom(:,:)                     = f_uptake(:,:) * pool_dom%phosphorus(:,:)
    loss_dom%phosphorus(:,:)            = loss_dom%phosphorus(:,:) + uptake_dom(:,:) * dtime / 1000000.0_wp 
    formation_microbial%phosphorus(:,:) = formation_microbial%phosphorus(:,:) + & 
                                            (microbial_pue_eff(:,:) * uptake_dom(:,:) + & 
                                            microbial_uptake_po4_sl(:,:)) * dtime / 1000000.0_wp 
   
    !! carbon13
    uptake_dom(:,:)                   = f_uptake(:,:) * pool_dom%carbon13(:,:)
    loss_dom%carbon13(:,:)            = loss_dom%carbon13(:,:) + uptake_dom(:,:) * dtime / 1000000._wp 
    formation_microbial%carbon13(:,:) = formation_microbial%carbon13(:,:) + & 
                                          microbial_cue_eff_mavg(:,:) * uptake_dom(:,:) * dtime / 1000000.0_wp 
    het_respiration_c13(:,:)          = het_respiration_c13(:,:) + (1._wp - microbial_cue_eff_mavg(:,:) ) * uptake_dom(:,:)

    !! carbon14
    uptake_dom(:,:)                   = f_uptake(:,:) * pool_dom%carbon14(:,:)
    loss_dom%carbon14(:,:)            = loss_dom%carbon14(:,:) + uptake_dom(:,:) * dtime / 1000000.0_wp 
    formation_microbial%carbon14(:,:) = formation_microbial%carbon14(:,:) + & 
                                          microbial_cue_eff_mavg(:,:) * uptake_dom(:,:) * dtime / 1000000.0_wp
    het_respiration_c14(:,:)          = het_respiration_c14(:,:) + (1._wp - microbial_cue_eff_mavg(:,:) ) * uptake_dom(:,:)

    !! nitrogen15
    net_mineralisation_nh4_n15(:,:)     = gross_mineralisation_nh4_n15(:,:) - microbial_uptake_nh4_n15_sl(:,:)
    net_mineralisation_no3_n15(:,:)     = - microbial_uptake_no3_n15_sl(:,:)
    uptake_dom(:,:)                     = f_uptake(:,:) * pool_dom%nitrogen15(:,:)
    loss_dom%nitrogen15(:,:)            = loss_dom%nitrogen15(:,:) + uptake_dom(:,:) * dtime / 1000000.0_wp 
    formation_microbial%nitrogen15(:,:) = formation_microbial%nitrogen15(:,:) + & 
                                          (microbial_nue * uptake_dom(:,:) + & 
                                           microbial_uptake_nh4_n15_sl(:,:) + microbial_uptake_no3_n15_sl(:,:)) * &
                                          dtime / 1000000.0_wp

  END SUBROUTINE calc_microbial_growth


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_jsm
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate depolymerisation rate of polymeric litter and microbial residue 
  !! corresponds Ahrens et al. 2015, eq. 1 
  !! updated to the more recent version of JSM, including microbial residue as alternative polymeric pool
  !!
  !! Input: maximum depolymerisation rate (vmax), temperature and moisture modifiers for
  !!        depolymerisation and half-saturation, polymeric/residue pool and microbial biomass
  !!
  !! Output: transfer of residue/polymeric pool to DOM
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_depolymerisation(nc, nsoil_sb, vmax_depolymerisation,rtm_depolymerisation,rmm_depolymerisation, &
                                   rtm_hsc,rmm_hsc,microbial_carbon,pool_polymeric,loss_polymeric,formation_dom)

    USE mo_sb_constants,          ONLY: km_depolymerisation
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                 !< dimensions
                                                                       nsoil_sb
    REAL(wp),                                         INTENT(in)    :: vmax_depolymerisation !< maximum depolymerisation rate
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_depolymerisation  !< temperature rate modifier for depolymerisation
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_depolymerisation  !< moisture rate modifier for depolymerisation
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_hsc               !< temperature rate modifier for half-saturation constants
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_hsc               !< moisture rate modifier for half-saturation constants
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: microbial_carbon      !< microbial biomass carbon [mol C/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_polymeric        !< polymeric/residue pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_polymeric        !< loss of polymeric/residue pool [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_dom         !< formation of DOM [mol/m3/timestep]
    ! ---------------------------
    ! 0.2 Local 
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1                   , &
                                                        hlp2                   , &
                                                        depolymerisation_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_depolymerisation'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...
    
    

    !> 1.0 determine potential depolymerisation rate [mol / m3 / timestep]
    !!
    hlp1(:,:) = 1._wp   
    hlp1(:,:) = microbial_carbon(:,:) / (km_depolymerisation * rtm_hsc(:,:) * rmm_hsc(:,:) + microbial_carbon(:,:)) 
    hlp2(:,:) = vmax_depolymerisation * rtm_depolymerisation(:,:) * rmm_depolymerisation(:,:) * hlp1(:,:) * dtime * & 
                 pool_polymeric%carbon(:,:) 
    
    !> 2.0 convert to turnover rate [1 / timestep] 
    !!
    WHERE(pool_polymeric%carbon(:,:) > eps8)
        depolymerisation_rate(:,:) = hlp2(:,:) / pool_polymeric%carbon(:,:)
    ELSEWHERE
        depolymerisation_rate(:,:) = 0.0_wp
    ENDWHERE

    !> 3.0 apply to all elements in this pool (should be a simple operator overload, but does currently not work safely)
    !!
    hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%carbon(:,:)
    loss_polymeric%carbon(:,:)      = loss_polymeric%carbon(:,:) + hlp2(:,:) 
    formation_dom%carbon(:,:)       = formation_dom%carbon(:,:) + hlp2(:,:) 

    hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%nitrogen(:,:)
    loss_polymeric%nitrogen(:,:)    = loss_polymeric%nitrogen(:,:) + hlp2(:,:) 
    formation_dom%nitrogen(:,:)     = formation_dom%nitrogen(:,:) + hlp2(:,:) 

    hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%phosphorus(:,:)
    loss_polymeric%phosphorus(:,:)  = loss_polymeric%phosphorus(:,:) + hlp2(:,:) 
    formation_dom%phosphorus(:,:)   = formation_dom%phosphorus(:,:) + hlp2(:,:) 

    hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%carbon13(:,:)
    loss_polymeric%carbon13(:,:)    = loss_polymeric%carbon13(:,:) + hlp2(:,:) 
    formation_dom%carbon13(:,:)     = formation_dom%carbon13(:,:) + hlp2(:,:) 

    hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%carbon14(:,:)
    loss_polymeric%carbon14(:,:)    = loss_polymeric%carbon14(:,:) + hlp2(:,:) 
    formation_dom%carbon14(:,:)     = formation_dom%carbon14(:,:) + hlp2(:,:) 

    hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%nitrogen15(:,:)
    loss_polymeric%nitrogen15(:,:)  = loss_polymeric%nitrogen15(:,:) + hlp2(:,:) 
    formation_dom%nitrogen15(:,:)   = formation_dom%nitrogen15(:,:) + hlp2(:,:) 

  END SUBROUTINE calc_depolymerisation


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate depolymerisation rate of polymeric litter and microbial residue 
  !! corresponds Ahrens et al. 2015, eq. 1 
  !! updated to the more recent version of JSM, including microbial residue as alternative polymeric pool
  !! enzyme partitioning simulated using SESAM approach, by Wutzler et al. 2017
  !!
  !! Input: maximum depolymerisation rate (vmax), temperature and moisture modifiers for
  !!        depolymerisation and half-saturation, polymeric/residue pool and microbial biomass
  !!
  !! Output: transfer of residue/polymeric pool to DOM
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_depolymerisation_SESAM(nc, nsoil_sb, rtm_hsc,rmm_hsc,rtm_depolymerisation,rmm_depolymerisation, &
                                         microbial_carbon, enzyme_frac_poly, enzyme_frac_residue, & 
                                         pool_polymeric,loss_polymeric,pool_residue,loss_residue,formation_dom, &
                                         enzyme_frac_poly_c, enzyme_frac_poly_n, enzyme_frac_poly_p, &
                                         enzyme_frac_poly_c_mavg, enzyme_frac_poly_n_mavg, enzyme_frac_poly_p_mavg)

    USE mo_sb_constants
    USE mo_jsb_math_constants,              ONLY: dtime, one_day
    USE mo_quincy_update_time_averages,     ONLY: mavg_period_tenzyme
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                    !< dimensions
                                                                       nsoil_sb
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_hsc                  !< temperature rate modifier for half-saturation constants
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_hsc                  !< moisture rate modifier for half-saturation constants
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_depolymerisation     !< temperature rate modifier for depolymerisation
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_depolymerisation     !< moisture rate modifier for depolymerisation
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: microbial_carbon         !< microbial biomass carbon [mol C/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_poly         !< fraction of enzyme for poly litter depolymerisation
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_residue      !< fraction of enzyme for microbial residue depolymerisation
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_polymeric           !< polymeric pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_polymeric           !< loss of polymeric pool [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_residue             !< residue pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_residue             !< loss of residue pool [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_dom            !< formation of DOM [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_poly_c       !< temporal fraction of enzyme for poly litter depolymerisation, C
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_poly_n       !< temporal fraction of enzyme for poly litter depolymerisation, N
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_poly_p       !< temporal fraction of enzyme for poly litter depolymerisation, P
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_poly_c_mavg  !< mavg of enzyme for poly litter depolymerisation, C
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_poly_n_mavg  !< mavg of enzyme for poly litter depolymerisation, N
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_poly_p_mavg  !< mavg of enzyme for poly litter depolymerisation, P
    ! ---------------------------
    ! 0.2 Local 
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1, hlp2, hlp3, &
                                                        poly_depolymerisation_rate, &
                                                        residue_depolymerisation_rate, &
                                                        depoly_c, depoly_n, depoly_p, temp_depoly
    LOGICAL, DIMENSION(nsoil_sb)                     :: tmp1,tmp2,tmp3
    INTEGER                                          :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_depolymerisation_SESAM'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...
    


    !> 0.9 init local variables
    !!
    depoly_c(:,:) = 0.0_wp
    depoly_n(:,:) = 0.0_wp
    depoly_p(:,:) = 0.0_wp


    !> 1.0 calculate polymeric litter depolymerisation rate 
    !!
    hlp1(:,:) = 1._wp   
    hlp1(:,:) = microbial_carbon(:,:) * enzyme_frac_poly / &
               (km_depolymerisation * rtm_hsc(:,:) * rmm_hsc(:,:) + microbial_carbon(:,:) * enzyme_frac_poly) 
    hlp2(:,:) = vmax_depolymerisation_litter * rtm_depolymerisation * rmm_depolymerisation * hlp1(:,:) * dtime * & 
                 pool_polymeric%carbon(:,:) 
    
    WHERE(pool_polymeric%carbon(:,:) > eps8)
        poly_depolymerisation_rate(:,:) = hlp2(:,:) / pool_polymeric%carbon(:,:)
    ELSEWHERE
        poly_depolymerisation_rate(:,:) = 0.0_wp
    ENDWHERE

    !> 2.0 calculate microbial residue depolymerisation rate 
    !!
    hlp1(:,:) = 1._wp   
    hlp1(:,:) = microbial_carbon(:,:) * enzyme_frac_residue / &
               (km_depolymerisation * rtm_hsc(:,:) * rmm_hsc(:,:) + microbial_carbon(:,:) * enzyme_frac_residue) 
    hlp2(:,:) = vmax_depolymerisation_residue * rtm_depolymerisation * rmm_depolymerisation * hlp1(:,:) * dtime * & 
                 pool_residue%carbon(:,:) 
    
    WHERE(pool_residue%carbon(:,:) > eps8)
        residue_depolymerisation_rate(:,:) = hlp2(:,:) / pool_residue%carbon(:,:)
    ELSEWHERE
        residue_depolymerisation_rate(:,:) = 0.0_wp
    ENDWHERE

    !> 3.0 apply to all elements in polymeric litter pool (should be a simple operator overload, but does currently not work safely)
    !!
    hlp2(:,:)                       = poly_depolymerisation_rate(:,:) * pool_polymeric%carbon(:,:)
    loss_polymeric%carbon(:,:)      = loss_polymeric%carbon(:,:) + hlp2(:,:) 
    formation_dom%carbon(:,:)       = formation_dom%carbon(:,:) + hlp2(:,:) 
    depoly_c(:,:) = depoly_c(:,:) + hlp2(:,:)

    hlp2(:,:)                       = poly_depolymerisation_rate(:,:) * pool_polymeric%nitrogen(:,:)
    loss_polymeric%nitrogen(:,:)    = loss_polymeric%nitrogen(:,:) + hlp2(:,:) 
    formation_dom%nitrogen(:,:)     = formation_dom%nitrogen(:,:) + hlp2(:,:) 
    depoly_n(:,:) = depoly_n(:,:) + hlp2(:,:)

    hlp2(:,:)                       = poly_depolymerisation_rate(:,:) * pool_polymeric%phosphorus(:,:)
    loss_polymeric%phosphorus(:,:)  = loss_polymeric%phosphorus(:,:) + hlp2(:,:) 
    formation_dom%phosphorus(:,:)   = formation_dom%phosphorus(:,:) + hlp2(:,:) 
    depoly_p(:,:) = depoly_p(:,:) + hlp2(:,:)

    hlp2(:,:)                       = poly_depolymerisation_rate(:,:) * pool_polymeric%carbon13(:,:)
    loss_polymeric%carbon13(:,:)    = loss_polymeric%carbon13(:,:) + hlp2(:,:) 
    formation_dom%carbon13(:,:)     = formation_dom%carbon13(:,:) + hlp2(:,:) 

    hlp2(:,:)                       = poly_depolymerisation_rate(:,:) * pool_polymeric%carbon14(:,:)
    loss_polymeric%carbon14(:,:)    = loss_polymeric%carbon14(:,:) + hlp2(:,:) 
    formation_dom%carbon14(:,:)     = formation_dom%carbon14(:,:) + hlp2(:,:) 

    hlp2(:,:)                       = poly_depolymerisation_rate(:,:) * pool_polymeric%nitrogen15(:,:)
    loss_polymeric%nitrogen15(:,:)  = loss_polymeric%nitrogen15(:,:) + hlp2(:,:) 
    formation_dom%nitrogen15(:,:)   = formation_dom%nitrogen15(:,:) + hlp2(:,:) 

    !> 4.0 apply to all elements in microbial residue pool
    !!
    hlp2(:,:)                       = residue_depolymerisation_rate(:,:) * pool_residue%carbon(:,:)
    loss_residue%carbon(:,:)        = loss_residue%carbon(:,:) + hlp2(:,:) 
    formation_dom%carbon(:,:)       = formation_dom%carbon(:,:) + hlp2(:,:) 
    depoly_c(:,:) = depoly_c(:,:) + hlp2(:,:)

    hlp2(:,:)                       = residue_depolymerisation_rate(:,:) * pool_residue%nitrogen(:,:)
    loss_residue%nitrogen(:,:)      = loss_residue%nitrogen(:,:) + hlp2(:,:) 
    formation_dom%nitrogen(:,:)     = formation_dom%nitrogen(:,:) + hlp2(:,:) 
    depoly_n(:,:) = depoly_n(:,:) + hlp2(:,:)

    hlp2(:,:)                       = residue_depolymerisation_rate(:,:) * pool_residue%phosphorus(:,:)
    loss_residue%phosphorus(:,:)    = loss_residue%phosphorus(:,:) + hlp2(:,:) 
    formation_dom%phosphorus(:,:)   = formation_dom%phosphorus(:,:) + hlp2(:,:) 
    depoly_p(:,:) = depoly_p(:,:) + hlp2(:,:)

    hlp2(:,:)                       = residue_depolymerisation_rate(:,:) * pool_residue%carbon13(:,:)
    loss_residue%carbon13(:,:)      = loss_residue%carbon13(:,:) + hlp2(:,:) 
    formation_dom%carbon13(:,:)     = formation_dom%carbon13(:,:) + hlp2(:,:) 

    hlp2(:,:)                       = residue_depolymerisation_rate(:,:) * pool_residue%carbon14(:,:)
    loss_residue%carbon14(:,:)      = loss_residue%carbon14(:,:) + hlp2(:,:) 
    formation_dom%carbon14(:,:)     = formation_dom%carbon14(:,:) + hlp2(:,:) 

    hlp2(:,:)                       = residue_depolymerisation_rate(:,:) * pool_residue%nitrogen15(:,:)
    loss_residue%nitrogen15(:,:)    = loss_residue%nitrogen15(:,:) + hlp2(:,:) 
    formation_dom%nitrogen15(:,:)   = formation_dom%nitrogen15(:,:) + hlp2(:,:) 

    !> 5.0 calculate the temporal enzyme partitioning for poly litter based on C, N, and P/m3/s

    WHERE (poly_depolymerisation_rate > 0.0_wp .AND. poly_depolymerisation_rate> 0.0_wp)
      !! alpha_C
      hlp1(:,:)               = poly_depolymerisation_rate(:,:) * pool_polymeric%carbon(:,:)
      hlp2(:,:)               = residue_depolymerisation_rate(:,:) * pool_residue%carbon(:,:)
      hlp3(:,:)               = km_depolymerisation * rtm_hsc(:,:) * rmm_hsc(:,:)
      enzyme_frac_poly_c(:,:) = calc_enzyme_allocation(hlp1(:,:), hlp2(:,:), hlp3(:,:), microbial_carbon(:,:))
      !! alpha_N
      hlp1(:,:)               = poly_depolymerisation_rate(:,:) * pool_polymeric%nitrogen(:,:)
      hlp2(:,:)               = residue_depolymerisation_rate(:,:) * pool_residue%nitrogen(:,:)
      enzyme_frac_poly_n(:,:) = calc_enzyme_allocation(hlp1(:,:), hlp2(:,:), hlp3(:,:), microbial_carbon(:,:))
      !! alpha_P
      hlp1(:,:)               = poly_depolymerisation_rate(:,:) * pool_polymeric%phosphorus(:,:)
      hlp2(:,:)               = residue_depolymerisation_rate(:,:) * pool_residue%phosphorus(:,:)
      enzyme_frac_poly_p(:,:) = calc_enzyme_allocation(hlp1(:,:), hlp2(:,:), hlp3(:,:), microbial_carbon(:,:))
    ENDWHERE
    
    !> 6.0 calculate the enzyme partitioning based on depolymerized C, N, P
    hlp1(:,:) = depoly_c(:,:)*microbial_cue
    hlp2(:,:) = depoly_n(:,:)*microbial_cn*microbial_nue
    hlp3(:,:) = depoly_p(:,:)*microbial_cn*microbial_np*microbial_pue

    WHERE (hlp1(:,:)>0.0_wp .AND. hlp2(:,:)>0.0_wp .AND. hlp3(:,:)>0.0_wp)
      WHERE (hlp1(:,:) < hlp2(:,:) .AND. hlp1(:,:)<hlp3(:,:))
        temp_depoly(:,:) = enzyme_frac_poly_c_mavg(:,:)
      ELSEWHERE (hlp2(:,:) < hlp1(:,:) .AND. hlp2(:,:)<hlp3(:,:))
        temp_depoly(:,:) = enzyme_frac_poly_n_mavg(:,:)
      ELSEWHERE (hlp3(:,:) < hlp1(:,:) .AND. hlp3(:,:)<hlp2(:,:))
        temp_depoly(:,:) = enzyme_frac_poly_p_mavg(:,:)
      ENDWHERE
      enzyme_frac_poly(:,:)    = MIN( MAX(enzyme_frac_poly(:,:) + (temp_depoly(:,:) - enzyme_frac_poly(:,:))&
                                 / mavg_period_tenzyme / one_day * dtime, 0.01_wp), 0.99_wp)
      enzyme_frac_residue(:,:) = 1._wp - enzyme_frac_poly(:,:)
    ENDWHERE
  END SUBROUTINE calc_depolymerisation_SESAM


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates the bulk density of the soil corrected by its organic component 
  !! derived from Federer et al. 1993, Can J. For. Res. 23 1026-1032 
  !!
  !! Input: carbon content by layer 
  !!
  !! Output: corrected bulk density of the soil 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_bulk_density_correction(bulk_soil_carbon_sl,bulk_dens_sl) RESULT(bulk_dens_corr_sl)

    USE mo_jsb_phys_constants,    ONLY: molar_mass_C 
    USE mo_sb_constants,          ONLY: rho_bulk_org, carbon_per_dryweight_SOM 

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: bulk_soil_carbon_sl     !< soil carbon contributing to density [mol C / m3] 
    REAL(wp), INTENT(in) :: bulk_dens_sl            !< bulk density of the mineral soil [kg / m3]
    REAL(wp)             :: bulk_dens_corr_sl       !< combined bulk density of the mineral and organic soil [kg / m3]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                    :: om_dens_sl       ! organic matter density [kg / m3] 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bulk_density_correction'

    !> 1.0 calculate organic matter density per soil layer [kg / m3]
    !!
    om_dens_sl = bulk_soil_carbon_sl * molar_mass_C / carbon_per_dryweight_SOM / 1000._wp 

    !> 2.0 actual bulk density of the soil 
    !!
    bulk_dens_corr_sl = MAX(rho_bulk_org, om_dens_sl + bulk_dens_sl - & 
                            om_dens_sl * bulk_dens_sl / rho_bulk_org) 

  END FUNCTION calc_bulk_density_correction 

  !-----------------------------------------------------------------------------------------------------
  ! Called in update_sb_pools and mo_sb_init
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates the maximum sorption capacity for organic material of a soil layer
  !! derived from COMISSION, the implict assumption here is that no association happens on organic material
  !!
  !! Input: bulk density of the soil (mineral and organic combined) 
  !!
  !! Output: sorption capacity for organic material 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_qmax_bulk_density_correction(bulk_soil_carbon_sl, bulk_dens_sl, qmax_om, qmax_min) RESULT(qmax)

    USE mo_jsb_phys_constants,    ONLY: molar_mass_C
    USE mo_sb_constants,          ONLY: carbon_per_dryweight_SOM, rho_bulk_org

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: bulk_soil_carbon_sl   !< organic carbon content of the soil layer [mol C / m3]
    REAL(wp), INTENT(in) :: bulk_dens_sl          !< bulk density of the mineral soil [kg / m3]
    REAL(wp), INTENT(in) :: qmax_om               !< maximium sorption capacity on OM [mol / kg]
    REAL(wp), INTENT(in) :: qmax_min              !< maximium sorption capacity on  mineral soil [mol / kg]
    REAL(wp)             :: qmax                  !< maximum sorption capacity of the soil [mol / m3]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                    :: vol_fraction_om
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_qmax_bulk_density_correction'

    ! calculate volume fraction filled by organic material 
    ! given assumed density and carbon content of OM
    vol_fraction_om = MIN(bulk_soil_carbon_sl * molar_mass_C / 1000._wp / & 
                          carbon_per_dryweight_SOM / rho_bulk_org, 1.0_wp)

    ! maximum sorption capacity per soil volume as the mix of OM and mineral sorption (mol/m3
    qmax = vol_fraction_om           * qmax_om  * rho_bulk_org + & 
           (1._wp - vol_fraction_om) * qmax_min * bulk_dens_sl

  END FUNCTION calc_qmax_bulk_density_correction


  !-----------------------------------------------------------------------------------------------------
  ! Called in update_sb_pools and mo_sb_init
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates bulk soil carbon used for the calculation of soil bulk density and sorption capacity 
  !!
  !! Input: soil pools 
  !!
  !! Output: total_bulk_soil_carbon 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_bulk_soil_carbon(nc, nsoil_sb, sb_pool) RESULT (bulk_soil_carbon_sl)

    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                         INTENT(in)    :: nc, &                   !< dimensions
                                                      nsoil_sb
    TYPE(t_sb_pool_a),               INTENT(in)    :: sb_pool                 !< pools [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb)              :: bulk_soil_carbon_sl
    ! ---------------------------
    ! 0.2 Local
    INTEGER                     :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bulk_soil_carbon'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Get local variables from options argument
    ! ..


    !> 1.0 
    !! 
    DO isoil = 1,nsoil_sb
       bulk_soil_carbon_sl(:,isoil) = &
            sb_pool%microbial%carbon(:,isoil) + &
            sb_pool%residue%carbon(:,isoil) + &
            sb_pool%residue_assoc%carbon(:,isoil) + &
            sb_pool%dom%carbon(:,isoil) + &
            sb_pool%dom_assoc%carbon(:,isoil) + &
            sb_pool%soluable_litter%carbon(:,isoil) + &
            sb_pool%polymeric_litter%carbon(:,isoil)
       IF(isoil > 1 ) THEN       
         bulk_soil_carbon_sl(:,isoil) = &
            bulk_soil_carbon_sl(:,isoil) + &
            sb_pool%woody_litter%carbon(:,isoil)
       ENDIF
    ENDDO

  END FUNCTION calc_bulk_soil_carbon


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from mo_sb_init 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates the fast_assoc_po4 based on Langmuir equilibrium and Labile P pool
  !!
  !! Input: Labile P pool, Qmax, Km, all in the unit [mol m-3]
  !!
  !! Output: fast_assoc_po4
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_fast_po4(Labile_p, Qmax, km, eca_enzyme_hlp) RESULT (fast_po4)
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: Labile_p                    !< labile P pool of the soil layer [mol P / m3]
    REAL(wp), INTENT(in) :: Qmax                        !< maximium sorption capacity of po4 [mol P/ m3]
    REAL(wp), INTENT(in) :: km                          !< sorption affinity of po4 [mol P/ m3]
    REAL(wp)             :: fast_po4                    !< fast associated po4 of the soil layer [mol P / m3]
    REAL(wp), INTENT(in), OPTIONAL :: eca_enzyme_hlp    !< helper item for enzyme abundance in eca approach [mol P / m3]
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fast_po4'
    REAL(wp) :: hlp1, hlp2, hlp3


    !> based on the assumption that labile P is the sum of both,
    !! and they follow the Langmuir equilibrium
    !! Labile_p - fast_po4 = fast_po4 * km / (Labile_p / eca_divisor_adsorp * (Qmax - fast_po4))
    IF (PRESENT(eca_enzyme_hlp)) THEN
      hlp1     = Labile_p + km
      hlp2     = Labile_p * Labile_p + Labile_p * Qmax + (Labile_p + eca_enzyme_hlp) * km
      hlp3     = Labile_p * Labile_p * Qmax
      fast_po4 = (hlp2 - sqrt(hlp2**2 - 4*hlp1*hlp3))/(2*hlp1)
    ELSE
      !! po4_assoc_fast = (p_labile_sl + Smax + K - sqrt((p_labile_sl + Smax + K)**2 - 4*p_labile_sl*Smax)) / 2
      fast_po4 = (Labile_p + Qmax + km- sqrt((Labile_p + Qmax + km)**2 - 4*Labile_p*Qmax)) / 2._wp
    ENDIF
  END FUNCTION calc_fast_po4


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates potential enzyme allocation coefficient based on SESAM concept, by Thomas Wutzler
  !!
  !! Input: potential depolymerisation rate (C, N, P), microbial biomass, half-saturation coefficient
  !!
  !! Output: potential enzyme coefficient for C, N, P
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_enzyme_allocation(dL, dR, km, Cmic) RESULT (alpha)
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: dL          !< organic carbon content of the soil layer [mol C / m3]
    REAL(wp), INTENT(in) :: dR          !< bulk density of the mineral soil [kg / m3]
    REAL(wp), INTENT(in) :: km          !< maximium sorption capacity on OM [mol / kg]
    REAL(wp), INTENT(in) :: Cmic        !< maximium sorption capacity on  mineral soil [mol / kg]
    REAL(wp)             :: alpha       !< maximum sorption capacity of the soil [mol / m3]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                    :: D    !< b**2-4*a*c in the quadratic equation
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_enzyme_allocation'

    !> Solution of the equation:
    !! dL/dR=alpha(km+alpha*Cmic)/((1-alpha)*(km+Cmic*(1-alpha)))
    !! only one of the solution fits in the range [0,1]
    D = 4*Cmic**2*dL*dR + 8*Cmic*km*dL*dR + km**2*dL**2 + 2*km**2*dL*dR + km**2*dR**2
    alpha = (dL*km+2*dL*Cmic+dR*km-sqrt(D))/(2*Cmic*(dL-dR))

  END FUNCTION calc_enzyme_allocation


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates sorption and desorption of organic material
  !! corresponds to the sorption and desorption terms in Ahrens et al. 2015, eq. 5 & 6
  !!
  !! Input: sorption temperature and moisture modifiers, sorption pools
  !!
  !! Output: sorption and desorption rate (mol/m2/timestep) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_om_sorption(nc, nsoil_sb, qmax_org,rtm_sorption,rmm_sorption,rtm_desorption, &
                              rmm_desorption,pool_dom,pool_dom_assoc,pool_residue, & 
                              pool_residue_assoc,formation_dom,formation_dom_assoc,formation_residue, &
                              formation_residue_assoc,loss_dom,loss_dom_assoc,loss_residue,loss_residue_assoc)

    USE mo_sb_constants,          ONLY: k_adsorpt_dom,k_desorpt_dom,k_adsorpt_det,k_desorpt_det 
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                   !< dimensions
                                                                       nsoil_sb
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: qmax_org 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_sorption            !< temperature rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_sorption            !< moisture rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_desorption          !< temperature rate modifier [-]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_desorption          !< moisture rate modifier [-]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_dom                !< DOM pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_residue            !< DOM pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_dom_assoc          !< DOM pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_residue_assoc      !< DOM pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_dom           !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_dom_assoc     !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_residue       !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_residue_assoc !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_dom                !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_dom_assoc          !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_residue            !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_residue_assoc      !< rate of DOM transport [mol/m3/timestep]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1                   , &
                                                        avail_sorption_sites   , &
                                                        sorption_rate          , &
                                                        desorption_rate        
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_om_sorption'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...


    !> 1.0 Determine available sorption sites [mol C m-3]
    !!
    avail_sorption_sites(:,:) = MAX(0._wp, qmax_org(:,:) - (pool_dom_assoc%carbon(:,:) + pool_residue_assoc%carbon(:,:)))

    !> 2.0 mineral association of DOM
    sorption_rate(:,:)   = k_adsorpt_dom * rtm_sorption(:,:)   * rmm_sorption(:,:) * avail_sorption_sites(:,:) * dtime 
    desorption_rate(:,:) = k_desorpt_dom * rtm_desorption(:,:) * rmm_desorption(:,:) * dtime 

    hlp1(:,:) = sorption_rate(:,:) * pool_dom%carbon(:,:)
    formation_dom_assoc%carbon(:,:) = formation_dom_assoc%carbon(:,:) + hlp1(:,:) 
    loss_dom%carbon(:,:)            = loss_dom%carbon(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%carbon(:,:)
    formation_dom%carbon(:,:)       = formation_dom%carbon(:,:) + hlp1(:,:) 
    loss_dom_assoc%carbon(:,:)      = loss_dom_assoc%carbon(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_dom%nitrogen(:,:)
    formation_dom_assoc%nitrogen(:,:) = formation_dom_assoc%nitrogen(:,:) + hlp1(:,:) 
    loss_dom%nitrogen(:,:)            = loss_dom%nitrogen(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%nitrogen(:,:)
    formation_dom%nitrogen(:,:)       = formation_dom%nitrogen(:,:) + hlp1(:,:) 
    loss_dom_assoc%nitrogen(:,:)      = loss_dom_assoc%nitrogen(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_dom%phosphorus(:,:)
    formation_dom_assoc%phosphorus(:,:) = formation_dom_assoc%phosphorus(:,:) + hlp1(:,:) 
    loss_dom%phosphorus(:,:)            = loss_dom%phosphorus(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%phosphorus(:,:)
    formation_dom%phosphorus(:,:)       = formation_dom%phosphorus(:,:) + hlp1(:,:) 
    loss_dom_assoc%phosphorus(:,:)      = loss_dom_assoc%phosphorus(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_dom%carbon13(:,:)
    formation_dom_assoc%carbon13(:,:) = formation_dom_assoc%carbon13(:,:) + hlp1(:,:) 
    loss_dom%carbon13(:,:)            = loss_dom%carbon13(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%carbon13(:,:)
    formation_dom%carbon13(:,:)       = formation_dom%carbon13(:,:) + hlp1(:,:) 
    loss_dom_assoc%carbon13(:,:)      = loss_dom_assoc%carbon13(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_dom%carbon14(:,:)
    formation_dom_assoc%carbon14(:,:) = formation_dom_assoc%carbon14(:,:) + hlp1(:,:) 
    loss_dom%carbon14(:,:)            = loss_dom%carbon14(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%carbon14(:,:)
    formation_dom%carbon14(:,:)       = formation_dom%carbon14(:,:) + hlp1(:,:) 
    loss_dom_assoc%carbon14(:,:)      = loss_dom_assoc%carbon14(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_dom%nitrogen15(:,:)
    formation_dom_assoc%nitrogen15(:,:) = formation_dom_assoc%nitrogen15(:,:) + hlp1(:,:) 
    loss_dom%nitrogen15(:,:)            = loss_dom%nitrogen15(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%nitrogen15(:,:)
    formation_dom%nitrogen15(:,:)       = formation_dom%nitrogen15(:,:) + hlp1(:,:) 
    loss_dom_assoc%nitrogen15(:,:)      = loss_dom_assoc%nitrogen15(:,:) + hlp1(:,:) 

    !> 3.0 mineral association of microbial residue 
    sorption_rate(:,:)   = k_adsorpt_det * rtm_sorption(:,:)   * rmm_sorption(:,:) * avail_sorption_sites(:,:) * dtime 
    desorption_rate(:,:) = k_desorpt_det * rtm_desorption(:,:) * rmm_desorption(:,:) * dtime 

    hlp1(:,:) = sorption_rate(:,:) * pool_residue%carbon(:,:)
    formation_residue_assoc%carbon(:,:) = formation_residue_assoc%carbon(:,:) + hlp1(:,:) 
    loss_residue%carbon(:,:)            = loss_residue%carbon(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%carbon(:,:)
    formation_residue%carbon(:,:)       = formation_residue%carbon(:,:) + hlp1(:,:) 
    loss_residue_assoc%carbon(:,:)      = loss_residue_assoc%carbon(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_residue%nitrogen(:,:)
    formation_residue_assoc%nitrogen(:,:) = formation_residue_assoc%nitrogen(:,:) + hlp1(:,:) 
    loss_residue%nitrogen(:,:)            = loss_residue%nitrogen(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%nitrogen(:,:)
    formation_residue%nitrogen(:,:)       = formation_residue%nitrogen(:,:) + hlp1(:,:) 
    loss_residue_assoc%nitrogen(:,:)      = loss_residue_assoc%nitrogen(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_residue%phosphorus(:,:)
    formation_residue_assoc%phosphorus(:,:) = formation_residue_assoc%phosphorus(:,:) + hlp1(:,:) 
    loss_residue%phosphorus(:,:)            = loss_residue%phosphorus(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%phosphorus(:,:)
    formation_residue%phosphorus(:,:)       = formation_residue%phosphorus(:,:) + hlp1(:,:) 
    loss_residue_assoc%phosphorus(:,:)      = loss_residue_assoc%phosphorus(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_residue%carbon13(:,:)
    formation_residue_assoc%carbon13(:,:) = formation_residue_assoc%carbon13(:,:) + hlp1(:,:) 
    loss_residue%carbon13(:,:)            = loss_residue%carbon13(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%carbon13(:,:)
    formation_residue%carbon13(:,:)       = formation_residue%carbon13(:,:) + hlp1(:,:) 
    loss_residue_assoc%carbon13(:,:)      = loss_residue_assoc%carbon13(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_residue%carbon14(:,:)
    formation_residue_assoc%carbon14(:,:) = formation_residue_assoc%carbon14(:,:) + hlp1(:,:) 
    loss_residue%carbon14(:,:)            = loss_residue%carbon14(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%carbon14(:,:)
    formation_residue%carbon14(:,:)       = formation_residue%carbon14(:,:) + hlp1(:,:) 
    loss_residue_assoc%carbon14(:,:)      = loss_residue_assoc%carbon14(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_residue%nitrogen15(:,:)
    formation_residue_assoc%nitrogen15(:,:) = formation_residue_assoc%nitrogen15(:,:) + hlp1(:,:) 
    loss_residue%nitrogen15(:,:)            = loss_residue%nitrogen15(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%nitrogen15(:,:)
    formation_residue%nitrogen15(:,:)       = formation_residue%nitrogen15(:,:) + hlp1(:,:) 
    loss_residue_assoc%nitrogen15(:,:)      = loss_residue_assoc%nitrogen15(:,:) + hlp1(:,:) 

  END SUBROUTINE calc_om_sorption 

END MODULE mo_sb_jsm_processes
