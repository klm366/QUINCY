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
!! @brief Contains the science routines for the JSM version in SLM 
!!
!! @author
!!  Soenke Zaehle 
!!
!! @par Revision History
!! First version by Soenke Zaehle (20170124): adopted from JSM_interface and JSM_core
!!
MODULE mo_sb_jsm_processes
#ifndef __NO_JSBACH__

#ifdef __QUINCY_STANDALONE__
#include "dsl_jsbach.inc"
#endif

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_jsb_math_constants,  ONLY: dtime, eps8


  IMPLICIT NONE

  PRIVATE
  ! PUBLIC :: calc_litter_turnover, calc_depolymerisation, & 
  !           calc_microbial_turnover, calc_microbial_growth, calc_km_qmax_correction, &
  !           calc_bulk_soil_carbon, calc_bulk_density_correction, calc_qmax_bulk_density_correction, &
  !           calc_om_sorption, calc_depolymerisation_SESAM, &
  !           calc_inorganic_uptake, calc_sinking_flux, calc_sourcing_flux_po4, &
  !           calc_fast_po4, calc_nitri_denitri_partitioning_rate, &
  !           calc_Psorption_parameter
  ! selection of the above functions needed for simple_1d
  PUBLIC :: calc_Psorption_parameter, calc_qmax_bulk_density_correction, &
            calc_bulk_soil_carbon, calc_fast_po4, calc_bulk_density_correction

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_jsm_processes'

CONTAINS


  ! !-----------------------------------------------------------------------------------------------------
  ! ! Sub Task called from update_sb_jsm
  ! ! 
  ! !
  ! !-----------------------------------------------------------------------------------------------------
  ! !> calculate turnover of litter pools according to first order decay 
  ! !! provides a link between the soluable and woody litter pools introduced in SLM, and the
  ! !! DOM and polymeric pools of JSM
  ! !!
  ! !! Input: temperature and moisture dependence of litter decomposition, litter pool size, and its
  ! !!        turnover time at standard conditions 
  ! !! Output: transfer of litter to its sink pool (soluable -> DOM; woody -> polymeric litter) 
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE calc_litter_turnover(nc, nsoil_sb          , &
  !                                 tau_litter_pool       , &
  !                                 fc_litter2sink        , &
  !                                 rtm_depolymerisation  , &
  !                                 rmm_depolymerisation  , & 
  !                                 pool_litter           , &
  !                                 loss_litter           , &
  !                                 formation_sink        , & 
  !                                 het_respiration       , & 
  !                                 het_respiration_c13   , & 
  !                                 het_respiration_c14)
  ! 
  !   dsl4slm_Use_memory_a_sb
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   INTEGER,                                          INTENT(in)    :: nc, &                         !< dimensions
  !                                                                      nsoil_sb
  !   REAL(wp),                                         INTENT(in)    :: tau_litter_pool       !< turnover time of the litter pool [1/s]
  !   REAL(wp),                                         INTENT(in)    :: fc_litter2sink        !< fraction of litter C entering sink pool (rest is respired) 
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_depolymerisation  !< temperature rate modifier [-]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_depolymerisation  !< moisture rate modifier [-]
  !   TYPE(t_sb_element_a),                             INTENT(in)    :: pool_litter           !< litter pool to decompose [mol/m3]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: loss_litter           !< reduction in litter pool [mol/m3/timestep]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: formation_sink        !< formation of SOM pool [mol/m3/timestep]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: het_respiration       !< heterotrophic respiration [micro-mol C/m3/s] 
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: het_respiration_c13   !< heterotrophic respiration [micro-mol 13C/m3/s] 
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: het_respiration_c14   !< heterotrophic respiration [pmol 14C/m3/s] 
  !   ! ---------------------------
  !   ! 0.2 Local 
  !   REAL(wp), DIMENSION(nc, nsoil_sb)                 :: hlp1 , &
  !                                                        decay_rate
  !   CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_turnover'
  !   ! ---------------------------
  !   ! 0.3 Declare Memory
  !   ! ..
  ! 
  ! 
  !   ! ------------------------------------------------------------------------------------------------------------
  !   ! Go science
  !   ! ------------------------------------------------------------------------------------------------------------
  ! 
  ! 
  !   !> 1.0 calulate decay per timestep
  !   !!
  !   decay_rate(:,:) = rtm_depolymerisation * rmm_depolymerisation /tau_litter_pool * dtime
  ! 
  !   !> 2.0 calculate source and sink terms
  !   !! assuming to change in composition or mass of the material to decompose
  !   !! 
  !   hlp1(:,:)                       = decay_rate(:,:) * pool_litter%carbon(:,:)
  !   loss_litter%carbon(:,:)         = loss_litter%carbon(:,:) + hlp1(:,:)
  !   formation_sink%carbon(:,:)      = formation_sink%carbon(:,:) + fc_litter2sink * hlp1(:,:) 
  !   het_respiration(:,:)            = het_respiration(:,:) + (1._wp - fc_litter2sink) * hlp1(:,:) * 1000000._wp / dtime 
  ! 
  !   hlp1(:,:)                       = decay_rate(:,:) * pool_litter%nitrogen(:,:)
  !   loss_litter%nitrogen(:,:)       = loss_litter%nitrogen(:,:) + hlp1(:,:)
  !   formation_sink%nitrogen(:,:)    = formation_sink%nitrogen(:,:) + hlp1(:,:) 
  ! 
  !   hlp1(:,:)                       = decay_rate(:,:) * pool_litter%phosphorus(:,:)
  !   loss_litter%phosphorus(:,:)     = loss_litter%phosphorus(:,:) + hlp1(:,:)
  !   formation_sink%phosphorus(:,:)  = formation_sink%phosphorus(:,:) + hlp1(:,:) 
  ! 
  !   hlp1(:,:)                       = decay_rate(:,:) * pool_litter%carbon13(:,:)
  !   loss_litter%carbon13(:,:)       = loss_litter%carbon13(:,:) + hlp1(:,:)
  !   formation_sink%carbon13(:,:)    = formation_sink%carbon13(:,:) + fc_litter2sink * hlp1(:,:) 
  !   het_respiration_c13(:,:)        = het_respiration_c13(:,:) + (1._wp - fc_litter2sink) * hlp1(:,:) * 1000000._wp / dtime 
  ! 
  !   hlp1(:,:)                       = decay_rate(:,:) * pool_litter%carbon14(:,:)
  !   loss_litter%carbon14(:,:)       = loss_litter%carbon14(:,:) + hlp1(:,:)
  !   formation_sink%carbon14(:,:)    = formation_sink%carbon14(:,:) + fc_litter2sink * hlp1(:,:) 
  !   het_respiration_c14(:,:)        = het_respiration_c14(:,:) + (1._wp - fc_litter2sink) * hlp1(:,:) * 1000000._wp / dtime 
  ! 
  !   hlp1(:,:)                       = decay_rate(:,:) * pool_litter%nitrogen15(:,:)
  !   loss_litter%nitrogen15(:,:)     = loss_litter%nitrogen15(:,:) + hlp1(:,:)
  !   formation_sink%nitrogen15(:,:)  = formation_sink%nitrogen15(:,:) + hlp1(:,:) 
  ! 
  ! END SUBROUTINE calc_litter_turnover
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! ! Sub Task called from update_sb_jsm
  ! !
  ! !-----------------------------------------------------------------------------------------------------
  ! !> calculate sinks of the bio-available inorganic nutrient pools
  ! !!  according to Zhu et al. (2017)
  ! !!
  ! !! Input: solute concentration, microbial biomass, fine root biomass, rate modifiers 
  ! !!
  ! !! Output: microbial and plant inorganic uptake rates
  ! !-----------------------------------------------------------------------------------------------------
  ! ELEMENTAL SUBROUTINE calc_inorganic_uptake(lctlib_vmax_uptake_n, &
  !                                            sb_adsorp_scheme, &
  !                                            solute, microbial_carbon, fine_root_carbon, &
  !                                            rtm_mic_uptake, rmm_mic_uptake, rtm_hsc, rmm_hsc, & 
  !                                            rtm_plant_uptake, rmm_plant_uptake, km1_upt, km2_upt, &
  !                                            vmax_mic_uptake, km_mic_uptake, vmax_plant_uptake, km_plant_uptake, & 
  !                                            k_mic_uptake, k_plant_uptake, microbial_uptake, plant_uptake)
  ! 
  !   USE mo_sb_constants, ONLY : f_plant_carrier_n, f_mic_carrier_n
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   REAL(wp), INTENT(in)            :: lctlib_vmax_uptake_n    !< lctlib paramter
  !   CHARACTER(len=*), INTENT(in)    :: sb_adsorp_scheme        !< eca_full; 2: eca_none
  !   REAL(wp), INTENT(in)            :: solute                  !< concentration of the solvent to be taken up (mol / m3)
  !   REAL(wp), INTENT(in)            :: microbial_carbon        !< microbial carbon density (mol / m3)
  !   REAL(wp), INTENT(in)            :: fine_root_carbon        !< fine_root density (mol / m3)
  !   REAL(wp), INTENT(in)            :: rtm_mic_uptake          !< temperature rate modifier for plant uptake
  !   REAL(wp), INTENT(in)            :: rmm_mic_uptake          !< moisture rate modifier for plant uptake
  !   REAL(wp), INTENT(in)            :: rtm_hsc                 !< temperature rate modifier for half_saturation constants
  !   REAL(wp), INTENT(in)            :: rmm_hsc                 !< moisture rate modifier for half_saturation constants
  !   REAL(wp), INTENT(in)            :: rtm_plant_uptake        !< temperature rate modifier for plant uptake
  !   REAL(wp), INTENT(in)            :: rmm_plant_uptake        !< moisture rate modifier for plant uptake
  !   REAL(wp), INTENT(in)            :: km1_upt                 !< low-affinity parameter for plant uptake [m3/mol]
  !   REAL(wp), INTENT(in)            :: km2_upt                 !< high-affinity parameter for plant uptake [m3/mol]
  !   REAL(wp), INTENT(in)            :: vmax_mic_uptake         !< maximium uptake rate of microbial biomass (micro-mol S / mol C / s)
  !   REAL(wp), INTENT(in)            :: km_mic_uptake           !< half-saturation constant of microbial nutrient uptake (mol S / m3)
  !   REAL(wp), INTENT(in)            :: vmax_plant_uptake       !< maximium uptake rate of microbial biomass (micro-mol S / mol C / s)
  !   REAL(wp), INTENT(in)            :: km_plant_uptake         !< half-saturation constant of microbial nutrient uptake (mol S / m3)
  !   REAL(wp), INTENT(in)            :: k_mic_uptake            !< reaction rate of microbial nutrient uptake (s -1)
  !   REAL(wp), INTENT(in)            :: k_plant_uptake          !< reaction rate of plant nutrient uptake (s -1)
  !   REAL(wp), INTENT(inout)         :: microbial_uptake        !< microbial uptake flux (micro-mol / m3 / s)
  !   REAL(wp), INTENT(out)           :: plant_uptake            !< plant uptake flux (micro-mol / m3 / s)
  !   ! ---------------------------
  !   ! 0.2 Local
  !   REAL(wp)                    :: hlp1,hlp2
  !   CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_inorganic_uptake'
  ! 
  ! 
  !   !> 1.0 uptake
  !   !! co-limited rates of N uptake for microbes and plants 
  !   SELECT CASE (sb_adsorp_scheme)
  !     CASE("eca_none")
  !       hlp1 = km1_upt / (rtm_hsc * rmm_hsc)
  !       hlp2 = km2_upt * rtm_hsc * rmm_hsc
  !       plant_uptake = fine_root_carbon * rtm_plant_uptake * rmm_plant_uptake * lctlib_vmax_uptake_n * & 
  !                      solute * ( hlp1 + 1._wp / ( hlp2 + solute )) 
  !     CASE("eca_part")
  !       hlp1 = rtm_hsc * rmm_hsc * km_mic_uptake
  !       hlp2 = rtm_hsc * rmm_hsc * km_plant_uptake
  !       microbial_uptake = rtm_mic_uptake * rmm_mic_uptake * vmax_mic_uptake * microbial_carbon * solute / &
  !                        (hlp1 + solute + hlp1/hlp2 * f_plant_carrier_n * fine_root_carbon + f_mic_carrier_n * microbial_carbon)
  !       plant_uptake     = rtm_plant_uptake * rmm_plant_uptake * vmax_plant_uptake * fine_root_carbon * solute / &
  !                        (hlp2 + solute + hlp2/hlp1 * f_mic_carrier_n * microbial_carbon + f_plant_carrier_n * fine_root_carbon)
  !     CASE("eca_full")
  !       hlp1 = rtm_hsc * rmm_hsc * km_mic_uptake
  !       hlp2 = rtm_hsc * rmm_hsc * km_plant_uptake
  !       microbial_uptake = rtm_mic_uptake * rmm_mic_uptake * k_mic_uptake * f_mic_carrier_n * microbial_carbon * solute / &
  !                        (hlp1 + solute + hlp1/hlp2 * f_plant_carrier_n * fine_root_carbon + f_mic_carrier_n * microbial_carbon)
  !       plant_uptake     = rtm_plant_uptake * rmm_plant_uptake * k_plant_uptake * f_plant_carrier_n * fine_root_carbon * solute / &
  !                        (hlp2 + solute + hlp2/hlp1 * f_mic_carrier_n * microbial_carbon + f_plant_carrier_n * fine_root_carbon)
  !   END SELECT
  ! 
  !   !> 1.1 limit to possible solution
  !   !!
  !   hlp1 = microbial_uptake + plant_uptake
  !   hlp2 = MAX(solute-eps8,0.0_wp) * 1000000._wp / dtime
  !   IF(hlp1 > hlp2 .AND. hlp1 > eps8)THEN
  !      microbial_uptake = microbial_uptake / hlp1 * hlp2
  !      plant_uptake     = plant_uptake     / hlp1 * hlp2 
  !   ENDIF
  ! 
  ! END SUBROUTINE calc_inorganic_uptake
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! ! Sub Task called from sb_jsm_main
  ! !
  ! !-----------------------------------------------------------------------------------------------------
  ! !> calculate the flux rates of P processes
  ! !! weathering and occlusion according to Wang et al. 2010; 
  ! !! desorption (slow exchange) rate based on the ANIMO model (Groenendijk et al. 2005)
  ! !! biocehmical mineralization is calculated using MM-kinetics assuming Mic C is the limiting factor
  ! !! Input: rate modifiers and source pools
  ! !!
  ! !! Output: biochemical P mineralization rate and loss of source pools, gross P mineralization, 
  ! !!         weathering, occlusion and (slow) desorption, microbial uptake
  ! !-----------------------------------------------------------------------------------------------------
  ! ELEMENTAL SUBROUTINE calc_sourcing_flux_po4(flag_sb_prescribe_po4, microbial_carbon, &                            ! in
  !                                           dom_carbon, dom_nitrogen, dom_phosphorus, &
  !                                           rtm_mic_uptake, rmm_mic_uptake, rtm_hsc, rmm_hsc, &
  !                                           rtm_sorption, rmm_sorption, rtm_desorption, rmm_desorption, &
  !                                           dom_assoc_po4, residue_assoc_po4, residue_po4, &
  !                                           dom_assoc_carbon, residue_assoc_carbon, residue_carbon, &
  !                                           solute_po4, solute_nh4, solute_no3, &
  !                                           assoc_fast, assoc_slow, po4_primary, & 
  !                                           bulk_dens_sl, volume_min_sl, root_carbon_sl, vmax_weath_mineral_sl, &   ! in
  !                                           loss_dom_assoc_po4, loss_residue_assoc_po4, loss_residue_po4, &         ! inout                                            
  !                                           biochem_min_po4, gross_mineralisation_po4, weathering_po4, &            ! out
  !                                           mic_upt_nh4, mic_upt_no3, mic_upt_po4, &
  !                                           occlusion_po4, slow_exchange_po4)
  ! 
  !   USE mo_sb_constants
  ! 
  !   IMPLICIT NONE
  ! 
  !   LOGICAL,  INTENT(in)            :: flag_sb_prescribe_po4     !< whether PO4 is prescribed or not
  !   REAL(wp), INTENT(in)            :: microbial_carbon          !< microbial carbon density (mol / m3)
  !   REAL(wp), INTENT(in)            :: dom_carbon                !< dom carbon density (mol / m3)
  !   REAL(wp), INTENT(in)            :: dom_nitrogen              !< dom phosphorus density (mol / m3)
  !   REAL(wp), INTENT(in)            :: dom_phosphorus            !< dom phosphorus density (mol / m3)
  !   REAL(wp), INTENT(in)            :: rtm_mic_uptake            !< temperature rate modifier for microbial (and plant) uptake
  !   REAL(wp), INTENT(in)            :: rmm_mic_uptake            !< moisture rate modifier for microbial (and plant) uptake
  !   REAL(wp), INTENT(in)            :: rtm_hsc                   !< temperature rate modifier for half_saturation constants
  !   REAL(wp), INTENT(in)            :: rmm_hsc                   !< moisture rate modifier for half_saturation constants
  !   REAL(wp), INTENT(in)            :: rtm_sorption              !< temperature rate modifier for sorption
  !   REAL(wp), INTENT(in)            :: rmm_sorption              !< moisture rate modifier for sorption
  !   REAL(wp), INTENT(in)            :: rtm_desorption            !< temperature rate modifier for desorption
  !   REAL(wp), INTENT(in)            :: rmm_desorption            !< moisture rate modifier for desorption
  !   REAL(wp), INTENT(in)            :: dom_assoc_po4             !< PO4 in the mineral associated dom pool [mol/m3]
  !   REAL(wp), INTENT(in)            :: residue_assoc_po4         !< PO4 in the mineral associated microbial necromass pool [mol/m3]
  !   REAL(wp), INTENT(in)            :: residue_po4               !< PO4 in the microbial necromass pool [mol/m3]
  !   REAL(wp), INTENT(in)            :: dom_assoc_carbon          !< dom carbon density (mol / m3)
  !   REAL(wp), INTENT(in)            :: residue_assoc_carbon      !< dom carbon density (mol / m3)
  !   REAL(wp), INTENT(in)            :: residue_carbon            !< dom carbon density (mol / m3)
  !   REAL(wp), INTENT(in)            :: solute_po4                !< soluble po4 pool [mol/m3]
  !   REAL(wp), INTENT(in)            :: solute_nh4                !< soluble nh4 pool [mol/m3]
  !   REAL(wp), INTENT(in)            :: solute_no3                !< soluble no3 pool [mol/m3]
  !   REAL(wp), INTENT(in)            :: assoc_slow                !< slow minerally associated PO4 pool [mol/m3]
  !   REAL(wp), INTENT(in)            :: assoc_fast                !< fast minerally associated PO4 pool [mol/m3]
  !   REAL(wp), INTENT(in)            :: po4_primary               !< primary PO4 pool [mol/m3]
  !   REAL(wp), INTENT(in)            :: bulk_dens_sl              !< soil bulk density [kg/m3]
  !   REAL(wp), INTENT(in)            :: volume_min_sl             !< mineral soil volume fraction [m3/m3]
  !   REAL(wp), INTENT(in)            :: root_carbon_sl            !< carbon of fine roots per soil layer, [mol/m3]
  !   REAL(wp), INTENT(in)            :: vmax_weath_mineral_sl     !< weathering coefficient of mineral soil [mol P /m3 /s]
  !   REAL(wp), INTENT(inout)         :: loss_dom_assoc_po4        !< loss of PO4 in the mineral associated dom pool [mol/m3]
  !   REAL(wp), INTENT(inout)         :: loss_residue_assoc_po4    !< loss of PO4 in the mineral associated microbial necromass pool [mol/m3]
  !   REAL(wp), INTENT(inout)         :: loss_residue_po4          !< loss of PO4 in the microbial necromass pool [mol/m3]
  !   REAL(wp), INTENT(out)           :: biochem_min_po4           !< phosphorus biochemical mineralization rate, [micro-mol P /m3 /s]
  !   REAL(wp), INTENT(out)           :: gross_mineralisation_po4  !< phosphorus gross(biological) mineralization rate, [micro-mol P /m3 /s]
  !   REAL(wp), INTENT(out)           :: weathering_po4            !< phosphorus weathering rate [micro-mol P / m3 / s]
  !   REAL(wp), INTENT(out)           :: mic_upt_nh4               !< microbial nh4 uptake rate [micro-mol N / m3 / s], for eca_none only
  !   REAL(wp), INTENT(out)           :: mic_upt_no3               !< microbial no3 uptake rate [micro-mol N / m3 / s], for eca_none only
  !   REAL(wp), INTENT(out)           :: mic_upt_po4               !< microbial po4 uptake rate [micro-mol P / m3 / s], for eca_none only
  !   REAL(wp), INTENT(out)           :: occlusion_po4             !< phosphorus occlusion rate, assoc_slow --> occluded [micro-mol P / m3 / s]
  !   REAL(wp), INTENT(out)           :: slow_exchange_po4         !< phosphorus desorption rate, assoc_fast --> assoc_slow [micro-mol P / m3 / s]
  ! 
  !   REAL(wp)                        :: hlp1, hlp2, hlp3, hlp4, hlp5, hlp6, hlp7, &
  !                                      uptake_c, f_uptake, uptake_n, uptake_p, mic_upt_n, &
  !                                      som_pc
  !   CHARACTER(len=*), PARAMETER     :: routine = TRIM(modname)//':calc_sourcing_flux_po4'
  ! 
  !   ! ------------------------------------------------------------------------------------------------------------
  !   ! Go science
  !   ! ------------------------------------------------------------------------------------------------------------
  ! 
  !   !> 1.0 Calculates the biochemical P mineralization, and the loss of the contributing pools.
  !   !! P is assumed to clevage out from dom associated, residue and residue associated 
  !   IF(.NOT. flag_sb_prescribe_po4)THEN
  !      hlp1 = (dom_assoc_po4 + residue_assoc_po4 + residue_po4)
  !      IF (hlp1 > eps8) THEN
  !         ! biomineralisation is regulated by solute_po4 and enzyme concentration (from root and microbe) 
  !         hlp2 = rtm_hsc * rmm_hsc * km_dip_biochem_po4 / (rtm_hsc * rmm_hsc * km_dip_biochem_po4 + solute_po4)
  !         hlp3 = (root_carbon_sl * f_plant_carrier_po4 + microbial_carbon * f_mic_carrier_po4) 
  !         hlp4 =  hlp3 / (rtm_hsc * rmm_hsc * km_mic_biochem_po4 + hlp3)
  !         ! instantaneous rate of biochemical mineralisation 1/s
  !         biochem_min_po4 = rtm_mic_uptake * rmm_mic_uptake * vmax_biochem_po4 * hlp2 * hlp4 * microbial_carbon
  ! 
  !         ! unit of loss terms: mol/m3/timestep
  !         ! loss from mineral associated DOM 
  !         IF (dom_assoc_carbon <eps8) THEN
  !           hlp5 = 0.0_wp
  !         ELSE
  !           som_pc = dom_assoc_po4 / dom_assoc_carbon
  !           IF (som_pc < biomin_thres_pc) THEN
  !             hlp5 = 0.0_wp
  !           ELSE 
  !             hlp5 = som_pc / ( km_pc_biochem_po4 + som_pc)
  !           ENDIF 
  !         ENDIF
  !         loss_dom_assoc_po4     = loss_dom_assoc_po4     + biochem_min_po4 * dom_assoc_po4 * hlp5 * dtime
  ! 
  !         ! loss from mineral associated residue 
  !         IF (residue_assoc_carbon <eps8) THEN
  !           hlp6 = 0.0_wp
  !         ELSE
  !           som_pc = residue_assoc_po4 / residue_assoc_carbon
  !           IF (som_pc < biomin_thres_pc) THEN
  !             hlp6 = 0.0_wp
  !           ELSE 
  !             hlp6 = som_pc / ( km_pc_biochem_po4 + som_pc)
  !           ENDIF 
  !         ENDIF
  !         loss_residue_assoc_po4 = loss_residue_assoc_po4 + biochem_min_po4 * residue_assoc_po4 * hlp6 * dtime
  ! 
  !         ! loss from microbial residue
  !         IF (residue_carbon <eps8) THEN
  !           hlp7 = 0._wp
  !         ELSE
  !           som_pc = residue_po4 / residue_carbon
  !           IF (som_pc < biomin_thres_pc) THEN
  !             hlp7 = 0.0_wp
  !           ELSE 
  !             hlp7 = som_pc / ( km_pc_biochem_po4 + som_pc)
  !           ENDIF
  !         ENDIF
  ! 
  !         loss_residue_po4       = loss_residue_po4       + biochem_min_po4 * residue_po4 * hlp7 * dtime
  ! 
  !         ! biochemical mineralisation in micro-mol P / m3 / s
  !         biochem_min_po4 = biochem_min_po4 * (dom_assoc_po4 * hlp5 + residue_po4 * hlp7 + residue_assoc_po4 * hlp6) * 1.e6_wp
  !      ELSE
  !         biochem_min_po4 = 0.0_wp
  !      ENDIF
  !   ELSE
  !      biochem_min_po4 = 0.0_wp
  !   ENDIF
  ! 
  !   !> 2.0 calculate the gross P mineralization, based on the potential depolymerization of DOM
  !   !!
  !   hlp1 = dom_carbon / (km_mic_uptake * rtm_hsc * rmm_hsc + dom_carbon )
  !   uptake_c = vmax_mic_uptake * rtm_mic_uptake * rmm_mic_uptake * hlp1 * microbial_carbon
  !   IF(dom_carbon > eps8) THEN
  !      f_uptake = MIN(uptake_c / dom_carbon, 1._wp * 1.e6_wp / dtime)
  !   ELSE
  !      f_uptake = 0.0_wp
  !   ENDIF
  !   uptake_p = f_uptake * dom_phosphorus
  !   uptake_n = f_uptake * dom_nitrogen
  !   uptake_c = f_uptake * dom_carbon
  !   gross_mineralisation_po4 = (1._wp - microbial_pue) * uptake_p
  !   !> 2.1 calculate the potential microbial N, P uptake
  !   mic_upt_po4 = MAX(0._wp, uptake_c * microbial_cue_max / microbial_cn / microbial_np - uptake_p * microbial_pue)
  !   mic_upt_n   = MAX(0._wp, uptake_c * microbial_cue_max / microbial_cn - uptake_n * microbial_nue)
  !   hlp1        = vmax_mic_uptake_nh4 * microbial_carbon * solute_nh4 / (solute_nh4 + km_mic_uptake_nh4)
  !   hlp2        = vmax_mic_uptake_no3 * microbial_carbon * solute_no3 / (solute_no3 + km_mic_uptake_no3)
  !   mic_upt_nh4 = mic_upt_n * hlp1 /(hlp1 + hlp2)
  !   mic_upt_no3 = mic_upt_n * hlp2 /(hlp1 + hlp2)
  ! 
  !   !> 3.0 Pi fluxes between different mineral P pools
  !   !! weathering, occlusion and slow exchange between sorbed and labile P pools
  !   hlp3 = (root_carbon_sl * f_plant_carrier_po4 + microbial_carbon * f_mic_carrier_po4) 
  !   hlp4 =  hlp3 / (rtm_hsc * rmm_hsc * km_root_enzyme_weath_po4 + hlp3)
  !   hlp5 = 0.1_wp ! minimum microbial biomass to exudate bio-weathering enzyme
  !   IF (root_carbon_sl < eps8 .AND. microbial_carbon < hlp5) hlp4 = 0._wp
  !   weathering_po4 = rtm_hsc / rmm_hsc * vmax_weath_mineral_sl * bulk_dens_sl * volume_min_sl * hlp4 * 1.e6_wp
  !   IF (weathering_po4 / 1.e6_wp > po4_primary / dtime) weathering_po4 = po4_primary / dtime * 1.e6_wp
  !   occlusion_po4 = assoc_slow * k_occlude_po4 * 1.e6_wp
  ! 
  !   hlp1 = rtm_desorption * k_desorpt_po4 * assoc_slow 
  !   hlp2 = rtm_sorption * k_adsorpt_po4 * assoc_fast     
  !   slow_exchange_po4 = (hlp2 - hlp1) * 1.e6_wp
  ! 
  ! END SUBROUTINE calc_sourcing_flux_po4
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! ! Sub Task called from update_sb_jsm for: po4, nh4, no3
  ! !
  ! !-----------------------------------------------------------------------------------------------------
  ! !> calculate sinks of the bio-available inorganic po4 and nh4 and no3 pools 
  ! !! ECA approach (Equilibrium Chemistry Approximation) modified based on \n
  ! !! Zhu et al. (2016) and Zhu et al. (2017), fast exchange of po4 and nh4 based on Van der Zee et al. (1989)
  ! !!
  ! !! Input: solute concentration, microbial biomass, fine root biomass, rate modifiers, km values, \n
  ! !!        sorption capacity and available sorption sites (po4, nh4)
  ! !! 
  ! !! Output: microbial and plant inorganic uptake rates, fast adsorption rates (po4, nh4), \n
  ! !!         nitrification/denitrification rates (nh4/no3)
  ! !-----------------------------------------------------------------------------------------------------
  ! ELEMENTAL SUBROUTINE calc_sinking_flux(lctlib_vmax_uptake, substrate, sb_adsorp_scheme, &           ! in
  !                                            microbial_carbon, fine_root_carbon, solute, &            ! in
  !                                            rtm_mic_uptake, rmm_mic_uptake, rtm_hsc, rmm_hsc, &      ! in
  !                                            rtm_plant_uptake, rmm_plant_uptake, km1_upt, km2_upt, &  ! in
  !                                            km_nitrification, km_denitrification, &                  ! in
  !                                            rtm_nit_denit_fication, rmm_nit_denit_fication, &        ! in
  !                                            anaerobic_volume_fraction, km_adsorpt_sl, &              ! in
  !                                            avail_sorption_sites, solute_available, &                ! in
  !                                            vmax_nitri_denitri, &      ! OPTIONAL (in)
  !                                            solute_alternative, &      ! OPTIONAL (in)
  !                                            microbial_uptake, &        ! inout
  !                                            plant_uptake, &            ! out
  !                                            fast_adsorpt, &            ! OPTIONAL (out)
  !                                            flux_nitri_denitri)        ! OPTIONAL (out)
  ! 
  !   USE mo_sb_constants
  ! 
  !   IMPLICIT NONE
  ! 
  !   REAL(wp)          , INTENT(in)    :: lctlib_vmax_uptake        !< lctlib uptake paramter, for eca_none
  !   CHARACTER(len=*)  , INTENT(in)    :: substrate                 !< po4, nh4, no3, determine the specific ECA paramaters
  !   CHARACTER(len=*)  , INTENT(in)    :: sb_adsorp_scheme          !< eca_full; eca_none; eca_part
  !   REAL(wp)          , INTENT(in)    :: microbial_carbon          !< microbial carbon density (mol / m3)
  !   REAL(wp)          , INTENT(in)    :: fine_root_carbon          !< fine_root density (mol / m3)
  !   REAL(wp)          , INTENT(in)    :: rtm_mic_uptake            !< temperature rate modifier for microbial uptake
  !   REAL(wp)          , INTENT(in)    :: rmm_mic_uptake            !< moisture rate modifier for microbial uptake
  !   REAL(wp)          , INTENT(in)    :: rtm_hsc                   !< temperature rate modifier for half_saturation constants
  !   REAL(wp)          , INTENT(in)    :: rmm_hsc                   !< moisture rate modifier for half_saturation constants
  !   REAL(wp)          , INTENT(in)    :: rtm_plant_uptake          !< temperature rate modifier for plant uptake
  !   REAL(wp)          , INTENT(in)    :: rmm_plant_uptake          !< moisture rate modifier for plant uptake
  !   REAL(wp)          , INTENT(in)    :: km1_upt                   !< low-affinity parameter for plant uptake [m3/mol]
  !   REAL(wp)          , INTENT(in)    :: km2_upt                   !< high-affinity parameter for plant uptake [m3/mol]
  !   REAL(wp)          , INTENT(in)    :: km_nitrification          !< half-saturation constant for carbon-limited (de-)nitrification [mol C m-3]
  !   REAL(wp)          , INTENT(in)    :: km_denitrification        !< half-saturation constant for nitrate-limited denitrifcation [mol N m-3]
  !   REAL(wp)          , INTENT(in)    :: rtm_nit_denit_fication    !< actual temperature response for de/nitrification
  !   REAL(wp)          , INTENT(in)    :: rmm_nit_denit_fication    !< actual moisture response for de/nitrification
  !   REAL(wp)          , INTENT(in)    :: anaerobic_volume_fraction !< actual anaerobic volume fraction
  !   REAL(wp)          , INTENT(in)    :: km_adsorpt_sl             !< half-saturation of solute sorption to mineral soil [mol / m3]
  !   REAL(wp)          , INTENT(in)    :: solute                    !< concentration of the solvent to be taken up (mol / m3)
  !   REAL(wp)          , INTENT(in)    :: avail_sorption_sites      !< available sorption sites [mol / m3]
  !   REAL(wp)          , INTENT(in)    :: solute_available          !< available solvent for all the sinking processes (mol / m3)
  !   REAL(wp), OPTIONAL, INTENT(in)    :: vmax_nitri_denitri        !< maximum nitrification/denitrification rate, for eca_none
  !   REAL(wp), OPTIONAL, INTENT(in)    :: solute_alternative        !< concentration of the alternative solvent to be taken up by plant (mol / m3)    
  !   REAL(wp)          , INTENT(inout) :: microbial_uptake          !< microbial uptake flux (micro-mol / m3 / s)
  !   REAL(wp)          , INTENT(out)   :: plant_uptake              !< plant uptake flux (micro-mol / m3 / s)
  !   REAL(wp), OPTIONAL, INTENT(out)   :: fast_adsorpt              !< fast adsorption flux (micro-mol / m3 / s)
  !   REAL(wp), OPTIONAL, INTENT(out)   :: flux_nitri_denitri        !< nitrification or denitrification flux (micro-mol / m3 / s)
  ! 
  ! 
  !   REAL(wp)                        :: mic_hlp, plant_hlp, mineral_hlp, nitri_hlp, &
  !                                      hlp1, hlp2, hlp3
  !   REAL(wp)                        :: f_mic_carrier, km1_mic_uptake, km2_mic_uptake, &
  !                                      f_plant_carrier, km_plant_uptake, km1_plant_uptake, &
  !                                      k_mic_uptake, k_plant_uptake, vmax_mic_uptake_sub, vmax_plant_uptake, &
  !                                      k_desorpt_f, vmax_nit_denit, km_nit_denit, & 
  !                                      k_nit_denit, f_nit_carrier, reaction_volume
  !   CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_sinking_flux'
  ! 
  !   ! ------------------------------------------------------------------------------------------------------------
  !   ! Go science
  !   ! ------------------------------------------------------------------------------------------------------------
  ! 
  !   !> 0.9 init local variables
  !   !! assign ECA specific parameter values based on substrate
  !   SELECT CASE (substrate)
  !     CASE ("po4")
  !       f_mic_carrier         = f_mic_carrier_po4
  !       km1_mic_uptake        = km_mic_uptake_po4
  !       km_plant_uptake       = km_plant_uptake_po4
  !       f_plant_carrier       = f_plant_carrier_po4
  !       k_mic_uptake          = k_mic_uptake_po4
  !       k_plant_uptake        = k_plant_uptake_po4
  !       vmax_mic_uptake_sub   = vmax_mic_uptake_po4
  !       vmax_plant_uptake     = vmax_plant_uptake_po4
  !       k_desorpt_f           = k_desorpt_f_po4
  !       km_nit_denit          = 0._wp
  !       k_nit_denit           = 0._wp
  !       f_nit_carrier         = 0._wp
  !       vmax_nit_denit        = 0._wp
  !       reaction_volume       = 0._wp
  !     CASE ("nh4")
  !       f_mic_carrier         = f_mic_carrier_n
  !       km1_mic_uptake        = km_mic_uptake_nh4
  !       km2_mic_uptake        = km_mic_uptake_no3
  !       km_plant_uptake       = km_plant_uptake_nh4
  !       km1_plant_uptake      = km_plant_uptake_no3
  !       f_plant_carrier       = f_plant_carrier_n
  !       k_mic_uptake          = k_mic_uptake_nh4
  !       k_plant_uptake        = k_plant_uptake_nh4
  !       vmax_mic_uptake_sub   = vmax_mic_uptake_nh4
  !       vmax_plant_uptake     = vmax_plant_uptake_nh4
  !       k_desorpt_f           = k_desorpt_nh4
  !       vmax_nit_denit        = vmax_nitrification_nh4
  !       km_nit_denit          = km_nitrification_nh4
  !       f_nit_carrier         = f_nit_carrier_nh4
  !       k_nit_denit           = k_nit_carrier_n
  !       reaction_volume       = 1._wp - anaerobic_volume_fraction
  !     CASE ("no3")
  !       f_mic_carrier         = f_mic_carrier_n
  !       km1_mic_uptake        = km_mic_uptake_no3
  !       km2_mic_uptake        = km_mic_uptake_nh4
  !       km_plant_uptake       = km_plant_uptake_no3
  !       km1_plant_uptake      = km_plant_uptake_nh4
  !       f_plant_carrier       = f_plant_carrier_n
  !       k_mic_uptake          = k_mic_uptake_no3
  !       k_plant_uptake        = k_plant_uptake_no3
  !       vmax_mic_uptake_sub   = vmax_mic_uptake_no3
  !       vmax_plant_uptake     = vmax_plant_uptake_no3
  !       k_desorpt_f           = 0._wp
  !       vmax_nit_denit        = vmax_denitrification_no3
  !       km_nit_denit          = km_denit_no3
  !       f_nit_carrier         = f_denit_carrier_no3
  !       k_nit_denit           = k_nit_carrier_n
  !       reaction_volume       = anaerobic_volume_fraction
  !   END SELECT
  ! 
  !   !> 1.0 uptake
  !   !! co-limited rates of substrate uptake for microbes and plants
  !   !!
  !   !! define uptake scheme to use, 1. ECA; 2. non-ECA
  !   SELECT CASE (sb_adsorp_scheme)
  !     !> 1.1 ECA approach
  !     !! 
  !     CASE("eca_full", "eca_part")
  !       mic_hlp     = km1_mic_uptake + solute_available + f_mic_carrier * microbial_carbon + &
  !                     km1_mic_uptake / km_plant_uptake * f_plant_carrier * fine_root_carbon 
  !       plant_hlp   = km_plant_uptake + solute_available + f_plant_carrier * fine_root_carbon + &
  !                     km_plant_uptake / km1_mic_uptake * f_mic_carrier * microbial_carbon 
  !       mineral_hlp = km_adsorpt_sl + solute_available + avail_sorption_sites + &
  !                     km_adsorpt_sl / km_plant_uptake * f_plant_carrier * fine_root_carbon + &
  !                     km_adsorpt_sl / km1_mic_uptake * f_mic_carrier * microbial_carbon 
  !       nitri_hlp   = km_nit_denit + solute_available + f_nit_carrier * microbial_carbon + &
  !                     km_nit_denit / km1_mic_uptake * f_mic_carrier * microbial_carbon + &
  !                     km_nit_denit / km_plant_uptake * f_plant_carrier * fine_root_carbon
  !       ! nh4 & no3
  !       !! including the effects of alternative solute and nitrification/denitrification process
  !       IF (substrate == "nh4" .OR. substrate == "no3") THEN
  !         plant_hlp = plant_hlp + &
  !                     km_plant_uptake / km1_plant_uptake * solute_alternative + &
  !                     km_plant_uptake / km_nit_denit * f_nit_carrier * microbial_carbon
  !         mic_hlp   = mic_hlp + &
  !                     km1_mic_uptake / km2_mic_uptake * solute_alternative + &
  !                     km1_mic_uptake / km_nit_denit * f_nit_carrier * microbial_carbon
  !         ! nh4
  !         !! including the effect of adsorption
  !         IF (substrate=="nh4") THEN
  !           mic_hlp      = mic_hlp + &
  !                          km1_mic_uptake / km_adsorpt_sl * avail_sorption_sites 
  !           plant_hlp    = plant_hlp + &
  !                          km_plant_uptake / km_adsorpt_sl * avail_sorption_sites 
  !           mineral_hlp  = mineral_hlp + &
  !                          km_adsorpt_sl / km_nit_denit * f_nit_carrier * microbial_carbon
  !           nitri_hlp    = nitri_hlp + &
  !                          km_nit_denit / km_adsorpt_sl * avail_sorption_sites
  !           hlp1         = k_desorpt_f / km_adsorpt_sl * solute 
  !           fast_adsorpt = MIN(1._wp / dtime, hlp1 * solute_available / mineral_hlp) * avail_sorption_sites * 1.e6_wp
  !         ENDIF
  !       ! po4
  !       !! including the effect of adsorption
  !       ELSE
  !         mic_hlp      = mic_hlp + &
  !                        km1_mic_uptake / km_adsorpt_sl * avail_sorption_sites 
  !         plant_hlp    = plant_hlp + &
  !                        km_plant_uptake / km_adsorpt_sl * avail_sorption_sites 
  !         hlp1         = k_desorpt_f / km_adsorpt_sl * solute 
  !         fast_adsorpt = MIN(1._wp / dtime, hlp1 * solute_available / mineral_hlp) * avail_sorption_sites * 1.e6_wp
  !       ENDIF
  !       ! eca_full, rates calculation based on enzyme abundance
  !       !! Calculate plant uptake, microbial uptake and nitrification/denitrification rates
  !       microbial_uptake   = MIN(rtm_mic_uptake * rmm_mic_uptake * k_mic_uptake * solute_available * f_mic_carrier * &
  !                            microbial_carbon / mic_hlp, microbial_uptake)
  !       plant_uptake       = rtm_plant_uptake * rmm_plant_uptake * k_plant_uptake * solute_available * f_plant_carrier_po4 * &
  !                            fine_root_carbon / plant_hlp 
  !       IF (substrate == "nh4" .OR. substrate == "no3") THEN
  !         flux_nitri_denitri = rtm_nit_denit_fication * rmm_nit_denit_fication * reaction_volume * k_nit_denit * &
  !                              solute_available * f_nit_carrier * microbial_carbon / nitri_hlp
  !       ENDIF
  !       ! eca_part, rates calculation based on vmax (prescribed)
  !       !! Calculate plant uptake, microbial uptake and nitrification/denitrification rates
  !       IF (sb_adsorp_scheme == "eca_part") THEN
  !         microbial_uptake   = rtm_mic_uptake * rmm_mic_uptake * vmax_mic_uptake_sub * &
  !                              solute_available * microbial_carbon / mic_hlp
  !         plant_uptake       = rtm_plant_uptake * rmm_plant_uptake * vmax_plant_uptake * &
  !                              solute_available * fine_root_carbon / plant_hlp
  !         IF (substrate == "nh4" .OR. substrate == "no3") THEN
  !           flux_nitri_denitri = rtm_nit_denit_fication * rmm_nit_denit_fication * reaction_volume * vmax_nit_denit * & 
  !                                solute_available * microbial_carbon / nitri_hlp
  !         ENDIF
  !       ENDIF
  ! 
  !     !> 1.2 non-ECA approach
  !     !! same as simple_1d model
  !     CASE("eca_none")
  !       plant_uptake   = rtm_plant_uptake * rmm_plant_uptake * fine_root_carbon * lctlib_vmax_uptake * &
  !                        solute * (km1_upt + 1._wp / (km2_upt + solute)) 
  !       IF (substrate == "po4" .OR. substrate == "nh4") THEN
  !         fast_adsorpt   = 0._wp
  !       ENDIF
  !       IF (substrate == "nh4" .OR. substrate == "no3") THEN
  !         flux_nitri_denitri = rtm_nit_denit_fication * rmm_nit_denit_fication * vmax_nitri_denitri * microbial_carbon / &
  !                              (km_nitrification + microbial_carbon)  * reaction_volume * solute
  !       ENDIF
  !       IF (substrate == "no3") THEN
  !         flux_nitri_denitri = flux_nitri_denitri / (km_denitrification + solute)
  !       ENDIF
  !   END SELECT
  ! 
  !   !> 2.0 
  !   !! check the availability of solute considering all sinking fluxes
  !   hlp2 = solute_available / dtime * 1.e6_wp
  !   SELECT CASE (substrate)
  !     CASE("po4")
  !       hlp1 = microbial_uptake + plant_uptake + fast_adsorpt
  !     CASE("nh4")
  !       hlp1 = microbial_uptake + plant_uptake + fast_adsorpt + flux_nitri_denitri
  !     CASE("no3")
  !       hlp1 = microbial_uptake + plant_uptake + flux_nitri_denitri
  !   END SELECT
  !   IF (hlp1 > hlp2 .AND. hlp1 > 0._wp)THEN
  !     microbial_uptake   = microbial_uptake   * hlp2 / hlp1
  !     plant_uptake       = plant_uptake       * hlp2 / hlp1
  !     SELECT CASE (substrate)
  !       CASE("po4")
  !         fast_adsorpt       = fast_adsorpt       * hlp2 / hlp1
  !       CASE("nh4")
  !         fast_adsorpt       = fast_adsorpt       * hlp2 / hlp1
  !         flux_nitri_denitri = flux_nitri_denitri * hlp2 / hlp1      
  !       CASE("no3")
  !         flux_nitri_denitri = flux_nitri_denitri * hlp2 / hlp1
  !     END SELECT
  !   ENDIF
  ! 
  ! END SUBROUTINE calc_sinking_flux
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! ! Sub Task to update_sb_jsm
  ! !
  ! !-----------------------------------------------------------------------------------------------------
  ! !> For the sb_nloss_scheme = fixed: calculates the amount of N lost in proportion to 
  ! !! net mineralisation \n
  ! !! For the sb_nloss_scheme = dynamic: calculates the partitioning of nitrification and denitrification
  ! !! to NO3, NOy, N2O and N2, and diagnoses the amount of leaching below the rooting zone
  ! !!
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE calc_nitri_denitri_partitioning_rate( &
  !                                   nc, &                         ! in
  !                                   nsoil_sb, &
  !                                   soil_depth_sl, &
  !                                   sb_nloss_scheme, &
  !                                   rtm_denitrification_act, &
  !                                   nitrification_no3, &          ! inout
  !                                   denitrification_n2, &
  !                                   nitrification_noy, &          ! out
  !                                   nitrification_n2o, &
  !                                   denitrification_noy, &
  !                                   denitrification_n2o, &
  !                                   net_mineralisation_nh4, & 
  !                                   net_mineralisation_nh4_n15, & 
  !                                   nitrification_no3_n15, &
  !                                   volatilisation_nh4, &
  !                                   volatilisation_nh4_n15, &
  !                                   emission_n2, &
  !                                   emission_n2_n15)
  ! 
  !   USE mo_sb_constants,            ONLY: f_nit_noy, f_nit_n2o, f_denit_noy, f_denit_n2o,floss_nmin,fnit_nmin
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   INTEGER,                                          INTENT(in)    :: nc, &                           !< dimensions
  !                                                                      nsoil_sb
  !   CHARACTER(len=*),                                 INTENT(in)    :: sb_nloss_scheme                 !< sb n loss scheme 
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: soil_depth_sl, &                !< soil depth per layer [m]
  !                                                                      net_mineralisation_nh4, &       !< net mineralisation rate of NH4 [micro-mol/m2/s]
  !                                                                      net_mineralisation_nh4_n15, &   !< net mineralisation rate of 15NH4 [micro-mol/m2/s]
  !                                                                      rtm_denitrification_act         !< temperature modifier for denitrification
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: nitrification_no3, &            !< nitrification to NO3 [micro-mol/m2/s]
  !                                                                      denitrification_n2              !< denitrification to N2 [micro-mol/m2/s]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(out)   :: nitrification_noy, &            !< nitrification to NOy [micro-mol/m2/s]
  !                                                                      nitrification_n2o, &            !< nitrification to N2O [micro-mol/m2/s]
  !                                                                      denitrification_noy, &          !< denitrification to NOy [micro-mol/m2/s]
  !                                                                      denitrification_n2o, &          !< denitrification to N2O [micro-mol/m2/s]
  !                                                                      nitrification_no3_n15, &        !< nitrification to 15NO3 [micro-mol/m2/s]
  !                                                                      volatilisation_nh4, &           !< volatilisation of NH4 [micro-mol/m2/s]
  !                                                                      volatilisation_nh4_n15          !< volatilisation of 15NH4 [micro-mol/m2/s]
  !   REAL(wp), DIMENSION(nc),                          INTENT(out)   :: emission_n2, &                  !< soil efflux of N2 [micro-mol/m2/s]
  !                                                                      emission_n2_n15                 !< soil efflux of 15N2 [micro-mol/m2/s]
  !   ! ---------------------------
  !   ! 0.2 Local
  !   CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_nitri_denitri_partitioning_rate'
  ! 
  ! 
  !   !> 1.0 case of dynamic N loss calculations, following Zaehle et al. 2011
  !   !!
  !   IF (sb_nloss_scheme=="dynamic") THEN
  !     ! nitrification rate to produce NO3, NOy and N2O from NH4 [mol m-3 timestep-1]
  !     nitrification_noy(:,:) = f_nit_noy * nitrification_no3(:,:)
  !     nitrification_n2o(:,:) = f_nit_n2o * nitrification_no3(:,:)
  !     nitrification_no3(:,:) = nitrification_no3(:,:) - nitrification_noy(:,:) - nitrification_n2o(:,:)
  ! 
  !     ! denitrification rate to produce NOy, N2O and N2 from NO3 [mol m-3 timestep-1]
  !     denitrification_noy(:,:) = f_denit_noy * rtm_denitrification_act(:,:) * denitrification_n2(:,:)
  !     denitrification_n2o(:,:) = f_denit_n2o * rtm_denitrification_act(:,:) * denitrification_n2(:,:)
  !     denitrification_n2(:,:)  = denitrification_n2(:,:) - denitrification_noy(:,:) - denitrification_n2o(:,:)
  ! 
  !   !> 2.0 case of fixed N losses: N loss proportional to net mineralisation
  !   !!
  !   ELSE
  !     WHERE(net_mineralisation_nh4(:,:) > 0.0_wp)
  !       ! N lost to atmosphere during processing of net mineralisation [mol m-3 timestep-1]
  !       volatilisation_nh4(:,:)      = net_mineralisation_nh4(:,:)     * floss_nmin 
  !       volatilisation_nh4_n15(:,:)  = net_mineralisation_nh4_n15(:,:) * floss_nmin 
  !       ! NH4 transformed to NO3 [mol m-3 timestep-1]
  !       nitrification_no3(:,:)     = net_mineralisation_nh4(:,:)     * (1._wp - floss_nmin) * fnit_nmin 
  !       nitrification_no3_n15(:,:) = net_mineralisation_nh4_n15(:,:) * (1._wp - floss_nmin) * fnit_nmin 
  !     ENDWHERE
  !     ! report volatilisation losses as N2 emission [mol m-2 timestep-1]
  !     emission_n2(:)     = SUM(volatilisation_nh4(:,:)     * soil_depth_sl(:,:), DIM=2)
  !     emission_n2_n15(:) = SUM(volatilisation_nh4_n15(:,:) * soil_depth_sl(:,:), DIM=2)
  !   ENDIF
  ! 
  ! END SUBROUTINE calc_nitri_denitri_partitioning_rate
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! ! Sub Task called from update_sb_pools
  ! !
  ! !-----------------------------------------------------------------------------------------------------
  ! !> calculate turnover of the microbial pool
  ! !! corresponds to the term pi * Cb in Ahrens et al. 2015, eq. 7
  ! !!
  ! !! Input: microbial biomass
  ! !!
  ! !! Output: transfer of microbial biomass to its sink pools (soluable -> DOM; polymeric -> residue)
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE calc_microbial_turnover(nc, nsoil_sb, pool_microbial, loss_microbial, formation_dom, formation_residue)
  ! 
  !   USE mo_sb_constants,          ONLY: k_mort_microbial, min_mic_biomass, f_soluable_residue, fc_p_res2dom, fc_n_res2dom
  !   dsl4slm_Use_memory_a_sb
  ! 
  !   IMPLICIT NONE
  ! 
  !   INTEGER,                        INTENT(in)    :: nc, &              !< dimensions
  !                                                    nsoil_sb
  !   TYPE(t_sb_element_a),           INTENT(in)    :: pool_microbial     !< microbial biomass [mol/m3]
  !   TYPE(t_sb_element_a),           INTENT(inout) :: loss_microbial     !< microbial death [mol/m3/timestep]
  !   TYPE(t_sb_element_a),           INTENT(inout) :: formation_dom      !< DOM formation from microbial death [mol/m3/timestep]
  !   TYPE(t_sb_element_a),           INTENT(inout) :: formation_residue  !< residue formation from microbial death [mol/m3/timestep]
  ! 
  !   REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1, &
  !                                                       p_res2dom, &
  !                                                       n_res2dom, &
  !                                                       mortality_rate
  !   REAL(wp)                                         :: hlp2, hlp3
  !   INTEGER                                          :: isoil
  !   CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_microbial_turnover'
  ! 
  !   ! ------------------------------------------------------------------------------------------------------------
  !   ! Go science
  !   ! ------------------------------------------------------------------------------------------------------------
  ! 
  !   !> 0.9 init local variables
  !   !! assuming that p is preferably recycled to dom in top soils, while n is not
  !   !! 
  !   hlp2 = 0.5_wp    ! scaling factor for fc_p_res2dom in the 2nd soil layer
  !   hlp3 = 0.25_wp   ! scaling factor for fc_p_res2dom for 3rd++ soil layer
  !   p_res2dom(:,:) = fc_p_res2dom
  !   n_res2dom(:,:) = fc_n_res2dom
  !   p_res2dom(:,2) = fc_p_res2dom * hlp2
  !   DO isoil = 3,nsoil_sb
  !     p_res2dom(:,isoil) = fc_p_res2dom * hlp3
  !   ENDDO
  ! 
  ! 
  !   !> 1.0 calculate microbial mortality rate, limited to minimum microbial biomass density (-> dormant microbes)
  !   !!
  !   mortality_rate(:,:) = k_mort_microbial * dtime
  !   hlp1(:,:)           = mortality_rate(:,:) * pool_microbial%carbon(:,:) 
  !   WHERE((pool_microbial%carbon(:,:) - hlp1(:,:)) < min_mic_biomass)
  !     mortality_rate(:,:) = MAX(0._wp, 1._wp - (min_mic_biomass / pool_microbial%carbon(:,:)))
  !   ENDWHERE
  ! 
  !   !> 2.0 calculate matter fluxes associated with microbial death
  !   !! a fraction turns to DOM, the remainder becomes residue
  !   !! the assumption is made that the composition of this material is not changed
  !   !!@todo: define C:N:P stoichiometries of DOM and residues from microbial death 
  !   !!
  !   hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%carbon(:,:)
  !   loss_microbial%carbon(:,:)        = loss_microbial%carbon(:,:) + hlp1(:,:)
  !   formation_dom%carbon(:,:)         = formation_dom%carbon(:,:) + f_soluable_residue * hlp1(:,:) 
  !   formation_residue%carbon(:,:)     = formation_residue%carbon(:,:) + (1._wp - f_soluable_residue) * hlp1(:,:) 
  ! 
  !   hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%nitrogen(:,:)
  !   loss_microbial%nitrogen(:,:)      = loss_microbial%nitrogen(:,:) + hlp1(:,:)
  !   formation_dom%nitrogen(:,:)       = formation_dom%nitrogen(:,:) + f_soluable_residue * hlp1(:,:) + &
  !                                       (1._wp - f_soluable_residue) * hlp1(:,:) * n_res2dom(:,:)
  !   formation_residue%nitrogen(:,:)   = formation_residue%nitrogen(:,:) + (1._wp - f_soluable_residue) * &
  !                                       (1._wp - n_res2dom(:,:)) * hlp1(:,:) 
  ! 
  !   hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%phosphorus(:,:)
  !   loss_microbial%phosphorus(:,:)    = loss_microbial%phosphorus(:,:) + hlp1(:,:)
  !   formation_dom%phosphorus(:,:)     = formation_dom%phosphorus(:,:) + f_soluable_residue * hlp1(:,:) + &
  !                                       (1._wp - f_soluable_residue) * hlp1(:,:) * p_res2dom(:,:)
  !   formation_residue%phosphorus(:,:) = formation_residue%phosphorus(:,:) + (1._wp - f_soluable_residue) * &
  !                                       (1._wp - p_res2dom(:,:)) * hlp1(:,:) 
  ! 
  !   hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%carbon13(:,:)
  !   loss_microbial%carbon13(:,:)      = loss_microbial%carbon13(:,:) + hlp1(:,:)
  !   formation_dom%carbon13(:,:)       = formation_dom%carbon13(:,:) + f_soluable_residue * hlp1(:,:) 
  !   formation_residue%carbon13(:,:)   = formation_residue%carbon13(:,:) + (1._wp - f_soluable_residue) * hlp1(:,:) 
  ! 
  !   hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%carbon14(:,:)
  !   loss_microbial%carbon14(:,:)      = loss_microbial%carbon14(:,:) + hlp1(:,:)
  !   formation_dom%carbon14(:,:)       = formation_dom%carbon14(:,:) + f_soluable_residue * hlp1(:,:) 
  !   formation_residue%carbon14(:,:)   = formation_residue%carbon14(:,:) + (1._wp - f_soluable_residue) * hlp1(:,:) 
  ! 
  !   hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%nitrogen15(:,:)
  !   loss_microbial%nitrogen15(:,:)    = loss_microbial%nitrogen15(:,:) + hlp1(:,:)
  !   formation_dom%nitrogen15(:,:)     = formation_dom%nitrogen15(:,:) + f_soluable_residue * hlp1(:,:) + &
  !                                       (1._wp - f_soluable_residue) * hlp1(:,:) * n_res2dom(:,:)
  !   formation_residue%nitrogen15(:,:) = formation_residue%nitrogen15(:,:) + (1._wp - f_soluable_residue) * &
  !                                       (1._wp - n_res2dom(:,:)) * hlp1(:,:) 
  ! 
  ! 
  ! END SUBROUTINE calc_microbial_turnover
  ! 
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! ! Sub Task called from update_sb_jsm
  ! !
  ! !-----------------------------------------------------------------------------------------------------
  ! !> calculate growth of the microbial pool 
  ! !! corresponds to Ahrens et al. 2015, eq. 2 
  ! !! extended to account for N and P limitation of microbial growth given a prescribed 
  ! !! stoichometry
  ! !! enzyme partitioning added according to Wutzler et al. 2017, by Lin
  ! !!
  ! !! Input: DOM and microbial pools, microbial nutrient uptake rate, temperature and moisture rate modifiers 
  ! !!
  ! !! Output: DOM consumption, microbial growth, hetertrophic respiration and net mineralisation
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE calc_microbial_growth(nc, nsoil_sb, sb_bnf_scheme, & 
  !                                  rtm_mic_uptake, rmm_mic_uptake, rtm_hsc, rmm_hsc, microbial_cue_eff_mavg, &
  !                                  pool_dom, pool_microbial, &
  !                                  microbial_uptake_nh4_sl, microbial_uptake_nh4_n15_sl, &
  !                                  microbial_uptake_no3_sl, microbial_uptake_no3_n15_sl, microbial_uptake_po4_sl, &
  !                                  het_respiration, het_respiration_c13, het_respiration_c14, &
  !                                  net_mineralisation_nh4, net_mineralisation_nh4_n15, &
  !                                  net_mineralisation_no3, net_mineralisation_no3_n15, net_mineralisation_po4, &
  !                                  microbial_cue_eff, microbial_nue_eff, microbial_pue_eff, asymb_n_fixation, &
  !                                  loss_dom, formation_microbial)
  ! 
  !   USE mo_isotope_util,          ONLY: calc_fractionation
  !   USE mo_sb_constants,          ONLY: vmax_mic_uptake,km_mic_uptake,microbial_cue,microbial_nue,microbial_pue, & 
  !                                       microbial_cn,microbial_np,eta_ammonification,microbial_cue_min,microbial_cue_max
  !   USE mo_jsb_math_constants,    ONLY: dtime
  !   dsl4slm_Use_memory_a_sb
  ! 
  !   IMPLICIT NONE
  ! 
  !   INTEGER,                                          INTENT(in)    :: nc, &                        !< dimensions
  !                                                                      nsoil_sb
  !   CHARACTER(len=*),                                 INTENT(in)    :: sb_bnf_scheme                !< veg n fixation scheme  
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_mic_uptake               !< temperature rate modifier [-]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_mic_uptake               !< moisture rate modifier [-]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_hsc                      !< temperature rate modifier [-]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_hsc                      !< moisture rate modifier [-]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: microbial_cue_eff_mavg       !< microbial effective CUE given NP constraints
  !   TYPE(t_sb_element_a),                             INTENT(in)    :: pool_dom                     !< DOM [mol/m3]
  !   TYPE(t_sb_element_a),                             INTENT(in)    :: pool_microbial               !< microbial biomass [mol/m3]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: microbial_uptake_nh4_sl      !< microbial uptake rate [micro-mol N/m3/s]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: microbial_uptake_nh4_n15_sl  !< microbial uptake rate [micro-mol 15N/m3/s]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: microbial_uptake_no3_sl      !< microbial uptake rate [micro-mol N/m3/s]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: microbial_uptake_no3_n15_sl  !< microbial uptake rate [micro-mol 15N/m3/s]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: microbial_uptake_po4_sl      !< microbial uptake rate [micro-mol P/m3/s]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: het_respiration              !< respiration rate [micro-mol C/m3/s]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: het_respiration_c13          !< respiration rate [micro-mol 13C/m3/s] 
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: het_respiration_c14          !< respiration rate [micro-mol 14C/m3/s]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: net_mineralisation_nh4       !< net mineralisation rate [micro-mol N/m3/s]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: net_mineralisation_nh4_n15   !< net mineralisation rate [micro-mol 15N/m3/s]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: net_mineralisation_no3       !< net mineralisation rate [micro-mol N/m3/s]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: net_mineralisation_no3_n15   !< net mineralisation rate [micro-mol 15N/m3/s]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: net_mineralisation_po4       !< net mineralisation rate [micro-mol P/m3/s]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: microbial_cue_eff            !< current effective CUE given NP constraints 
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: microbial_nue_eff            !< current effective NUE given NP constraints 
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: microbial_pue_eff            !< current effective PUE given NP constraints 
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: asymb_n_fixation             !< current effective CUE given NP constraints 
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: loss_dom                     !< microbial consumption of DOM [mol/m3/timestep]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: formation_microbial          !< microbial growth on DOM [mol/m3/timestep]
  ! 
  !   REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1, hlp2, hlp3, &
  !                                                       f_uptake, &
  !                                                       uptake_dom  
  !   REAL(wp), DIMENSION(nc, nsoil_sb)                :: growth_c,&
  !                                                       growth_n,&
  !                                                       growth_p,&
  !                                                       growth,&
  !                                                       uptake_c,&
  !                                                       uptake_n,&
  !                                                       uptake_p 
  !   REAL(wp), DIMENSION(nc, nsoil_sb)                :: scal
  !   REAL(wp), DIMENSION(nc, nsoil_sb)                :: gross_mineralisation_nh4,&
  !                                                       gross_mineralisation_nh4_n15,&
  !                                                       gross_mineralisation_po4 
  !   INTEGER                                          :: isoil
  !   CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_microbial_growth'
  !   ! ---------------------------
  !   ! 0.3 Declare Memory
  !   ! ...
  ! 
  ! 
  ! 
  !   !> 1.0 potential microbial DOM consumption (micro-mol / m3 / s) 
  !   !! 
  !   hlp1(:,:)     = pool_dom%carbon(:,:) / (km_mic_uptake * rtm_hsc(:,:) * rmm_hsc(:,:) + pool_dom%carbon(:,:) )
  !   uptake_c(:,:) = vmax_mic_uptake * rtm_mic_uptake(:,:) * rmm_mic_uptake(:,:) * hlp1(:,:) * pool_microbial%carbon(:,:) 
  !   WHERE(pool_dom%carbon(:,:) > eps8)
  !      f_uptake(:,:) = MIN(uptake_c(:,:) / pool_dom%carbon(:,:), 1._wp * 1.e6_wp/dtime)
  !   ELSEWHERE
  !      f_uptake(:,:) = 0.0_wp
  !   ENDWHERE
  !   uptake_n(:,:) = f_uptake(:,:) * pool_dom%nitrogen(:,:) 
  !   uptake_p(:,:) = f_uptake(:,:) * pool_dom%phosphorus(:,:)
  ! 
  !   !> 1.1 microbial growth rates in micro-mol C / m3 / s, given by C, N and P availability, and the rate-limiting growth rate
  !   !! 
  !   growth_c(:,:) = microbial_cue_max * uptake_c(:,:) 
  !   growth_n(:,:) = (microbial_nue * uptake_n(:,:) + &
  !                    microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:)) * microbial_cn
  !   growth_p(:,:) = (microbial_pue * uptake_p(:,:) + &
  !                     microbial_uptake_po4_sl(:,:)) * microbial_np * microbial_cn
  !   WHERE(growth_c(:,:) <= growth_n(:,:) .AND. growth_c(:,:) <= growth_p(:,:))
  !      ! carbon uptake is limiting
  !      growth(:,:) = growth_c(:,:)
  !   ELSEWHERE(growth_n(:,:) < growth_c(:,:) .AND. growth_n(:,:) <= growth_p(:,:))
  !      ! nitrogen uptake is limiting
  !      growth(:,:) = growth_n(:,:)
  !   ELSEWHERE(growth_p(:,:) < growth_c(:,:) .AND. growth_p(:,:) < growth_n(:,:))
  !      ! phosphorus uptake is limiting
  !      growth(:,:) = growth_p(:,:)
  !   ENDWHERE
  ! 
  !   !> 2.0 Adjust carbon-use efficiency such that the growth is stoichiometrically balanced
  !   !!
  !   microbial_cue_eff(:,:) = microbial_cue_max
  !   WHERE(growth(:,:) < growth_c(:,:) .AND. uptake_c(:,:) > eps8)
  !      ! excess carbon is going to be respired
  !      microbial_cue_eff(:,:) = MAX(microbial_cue_min,growth(:,:) / uptake_c(:,:))
  !   ELSEWHERE(growth(:,:) < growth_c(:,:))
  !      microbial_cue_eff(:,:) = microbial_cue_min
  !   ENDWHERE
  ! 
  !   !> 2.1 actual growth given actual CUE, adjust uptake rate of DOM to the most 
  !   !!     limiting nutrient
  !   !!
  !   !! growth_n and growth_p are unchanged, only growth_c needs to be adjusted to current 
  !   !! CUE
  !   growth_c(:,:) = microbial_cue_eff_mavg(:,:) * uptake_c(:,:) 
  !   IF(sb_bnf_scheme=='unlimited')THEN
  !      asymb_n_fixation(:,:) = MAX(0.0_wp, growth_c(:,:) / microbial_cn - &
  !                                          microbial_nue * f_uptake(:,:) * pool_dom%nitrogen(:,:) - &
  !                                          microbial_uptake_nh4_sl(:,:) - microbial_uptake_no3_sl(:,:))
  !      microbial_uptake_nh4_sl(:,:) = microbial_uptake_nh4_sl(:,:) + asymb_n_fixation(:,:)
  !   ELSE
  !      WHERE(growth_c(:,:) <= growth_n(:,:) .AND. growth_c(:,:) <= growth_p(:,:))
  !         ! carbon uptake is limiting, full uptake rate is realised
  !         f_uptake(:,:) = f_uptake(:,:) 
  !      ELSEWHERE(growth_n(:,:) < growth_c(:,:) .AND. growth_n(:,:) <= growth_p(:,:).AND.pool_dom%carbon(:,:) > eps8)
  !         ! nitrogen uptake is limiting
  !         ! calculate uptake rate that satisfies N constraint and CUE
  !         f_uptake(:,:) = MAX(0.0_wp,MIN(f_uptake(:,:), & 
  !                                        ((microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:)) * microbial_cn) / & 
  !                                        (microbial_cue_eff_mavg(:,:) * pool_dom%carbon(:,:) - & 
  !                                         microbial_nue * pool_dom%nitrogen(:,:) * microbial_cn)))
  !      ELSEWHERE(growth_p(:,:) < growth_c(:,:) .AND. growth_p(:,:) < growth_n(:,:) .AND. pool_dom%carbon(:,:) > eps8)
  !         ! phosphorus uptake is limiting
  !         ! calculate uptake rate that satisfies P constraint and CUE
  !         f_uptake(:,:) = MAX(0.0_wp,MIN(f_uptake(:,:), &
  !                                        (microbial_uptake_po4_sl(:,:) * microbial_np * microbial_cn) / & 
  !                                        (microbial_cue_eff_mavg(:,:) * pool_dom%carbon(:,:) - & 
  !                                         microbial_pue * pool_dom%phosphorus(:,:) * microbial_np * microbial_cn)))
  !      ELSEWHERE
  !         f_uptake(:,:) = 0.0_wp
  !      ENDWHERE
  ! 
  !      growth_c(:,:) = microbial_cue_eff_mavg(:,:) * f_uptake(:,:) * pool_dom%carbon(:,:) 
  !      uptake_n(:,:) = f_uptake(:,:) * pool_dom%nitrogen(:,:) 
  !      uptake_p(:,:) = f_uptake(:,:) * pool_dom%phosphorus(:,:)
  !      growth_n(:,:) = (microbial_nue * uptake_n(:,:) + &
  !                       microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:)) * microbial_cn
  !      growth_p(:,:) = (microbial_pue * uptake_p(:,:) + &
  !                        microbial_uptake_po4_sl(:,:)) * microbial_np * microbial_cn
  ! 
  !   ENDIF
  ! 
  !   !> 2.2 Adjust nitrogen-use efficiency and immobilisation flux such that the growth is stoichiometrically balanced
  !   !!     and calculate gross mineralisation rate (and its isotopic signature)
  !   microbial_nue_eff(:,:) = microbial_nue
  !   WHERE(growth_c(:,:) < growth_n(:,:)) 
  !      ! inorganic N immobilisation flux in agreement with C growth
  !      hlp1(:,:) = MAX(growth_c(:,:) / microbial_cn - microbial_nue * uptake_n(:,:),0.0_wp)
  !      WHERE((microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:)) > eps8)
  !         scal(:,:) = hlp1(:,:) / (microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:)) 
  !      ELSEWHERE
  !         scal(:,:) = 0.0_wp
  !      ENDWHERE
  !      microbial_uptake_nh4_sl(:,:)     = microbial_uptake_nh4_sl(:,:)     * scal(:,:) 
  !      microbial_uptake_nh4_n15_sl(:,:) = microbial_uptake_nh4_n15_sl(:,:) * scal(:,:)
  !      microbial_uptake_no3_sl(:,:)     = microbial_uptake_no3_sl(:,:)     * scal(:,:) 
  !      microbial_uptake_no3_n15_sl(:,:) = microbial_uptake_no3_n15_sl(:,:) * scal(:,:)
  !      WHERE(uptake_n(:,:) > eps8)
  !         microbial_nue_eff(:,:) = MIN(microbial_nue, &
  !                                     ( growth_c(:,:) / &
  !                                       microbial_cn - microbial_uptake_nh4_sl(:,:) - microbial_uptake_no3_sl(:,:) ) / &
  !                                     uptake_n(:,:))
  !      ENDWHERE
  !   ENDWHERE
  !   gross_mineralisation_nh4(:,:)     = (1._wp - microbial_nue_eff(:,:)) * uptake_n(:,:) 
  !   gross_mineralisation_nh4_n15(:,:) = gross_mineralisation_nh4(:,:) * calc_fractionation(pool_dom%nitrogen(:,:), & 
  !                                                                                          pool_dom%nitrogen15(:,:),& 
  !                                                                                          eta_ammonification)
  ! 
  !   !> 2.3 Adjust phosphorus-use efficiency and immobilisation flux such that the growth is stoichiometrically balanced
  !   !!     and calculate gross mineralisation rate 
  !   microbial_pue_eff(:,:) = microbial_pue
  !   WHERE(growth_c(:,:) < growth_p(:,:))  
  !      ! inorganic P immobilisation flux, in agreement with C and N uptake
  !      microbial_uptake_po4_sl(:,:) = MAX(growth_c(:,:) / microbial_cn / microbial_np - microbial_pue * uptake_p(:,:),0.0_wp)
  !      ! In the case that the source C:N:P is too tight for microbes, cause gross mineralisation
  !      WHERE(uptake_p(:,:) > eps8)
  !         microbial_pue_eff(:,:) = ( growth_c(:,:) / microbial_cn / microbial_np - microbial_uptake_po4_sl(:,:) ) / uptake_p(:,:)
  !      ENDWHERE
  !   ENDWHERE
  !   gross_mineralisation_po4(:,:) = (1._wp - microbial_pue_eff(:,:)) * uptake_p(:,:)
  ! 
  !   !> 3.0 update fluxes (note that loss and formation of organic matter is in mol/m3/timestep, whereas fluxes are in micro-mol/m3/s)
  !   !! carbon
  !   uptake_dom(:,:)                 = f_uptake(:,:) * pool_dom%carbon(:,:)
  !   loss_dom%carbon(:,:)            = loss_dom%carbon(:,:) + uptake_dom(:,:) * dtime / 1000000._wp
  !   formation_microbial%carbon(:,:) = formation_microbial%carbon(:,:) + & 
  !                                       microbial_cue_eff_mavg(:,:) * uptake_dom(:,:) * dtime / 1000000._wp 
  !   het_respiration(:,:)            = het_respiration(:,:) + & 
  !                                       (1._wp - microbial_cue_eff_mavg(:,:)) * uptake_dom(:,:) 
  ! 
  !   !! nitrogen
  !   net_mineralisation_nh4(:,:)       = gross_mineralisation_nh4(:,:) - microbial_uptake_nh4_sl(:,:)
  !   net_mineralisation_no3(:,:)       = - microbial_uptake_no3_sl(:,:)
  !   uptake_dom(:,:)                   = f_uptake(:,:) * pool_dom%nitrogen(:,:)
  !   loss_dom%nitrogen(:,:)            = loss_dom%nitrogen(:,:) + uptake_dom(:,:) * dtime / 1000000.0_wp 
  !   formation_microbial%nitrogen(:,:) = formation_microbial%nitrogen(:,:) + & 
  !                                         (microbial_nue_eff(:,:) * uptake_dom(:,:) + & 
  !                                         microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:)) * dtime / 1000000.0_wp 
  ! 
  !   !! phosphorus
  !   net_mineralisation_po4(:,:)         = gross_mineralisation_po4(:,:) - microbial_uptake_po4_sl(:,:)
  !   uptake_dom(:,:)                     = f_uptake(:,:) * pool_dom%phosphorus(:,:)
  !   loss_dom%phosphorus(:,:)            = loss_dom%phosphorus(:,:) + uptake_dom(:,:) * dtime / 1000000.0_wp 
  !   formation_microbial%phosphorus(:,:) = formation_microbial%phosphorus(:,:) + & 
  !                                           (microbial_pue_eff(:,:) * uptake_dom(:,:) + & 
  !                                           microbial_uptake_po4_sl(:,:)) * dtime / 1000000.0_wp 
  ! 
  !   !! carbon13
  !   uptake_dom(:,:)                   = f_uptake(:,:) * pool_dom%carbon13(:,:)
  !   loss_dom%carbon13(:,:)            = loss_dom%carbon13(:,:) + uptake_dom(:,:) * dtime / 1000000._wp 
  !   formation_microbial%carbon13(:,:) = formation_microbial%carbon13(:,:) + & 
  !                                         microbial_cue_eff_mavg(:,:) * uptake_dom(:,:) * dtime / 1000000.0_wp 
  !   het_respiration_c13(:,:)          = het_respiration_c13(:,:) + (1._wp - microbial_cue_eff_mavg(:,:) ) * uptake_dom(:,:)
  ! 
  !   !! carbon14
  !   uptake_dom(:,:)                   = f_uptake(:,:) * pool_dom%carbon14(:,:)
  !   loss_dom%carbon14(:,:)            = loss_dom%carbon14(:,:) + uptake_dom(:,:) * dtime / 1000000.0_wp 
  !   formation_microbial%carbon14(:,:) = formation_microbial%carbon14(:,:) + & 
  !                                         microbial_cue_eff_mavg(:,:) * uptake_dom(:,:) * dtime / 1000000.0_wp
  !   het_respiration_c14(:,:)          = het_respiration_c14(:,:) + (1._wp - microbial_cue_eff_mavg(:,:) ) * uptake_dom(:,:)
  ! 
  !   !! nitrogen15
  !   net_mineralisation_nh4_n15(:,:)     = gross_mineralisation_nh4_n15(:,:) - microbial_uptake_nh4_n15_sl(:,:)
  !   net_mineralisation_no3_n15(:,:)     = - microbial_uptake_no3_n15_sl(:,:)
  !   uptake_dom(:,:)                     = f_uptake(:,:) * pool_dom%nitrogen15(:,:)
  !   loss_dom%nitrogen15(:,:)            = loss_dom%nitrogen15(:,:) + uptake_dom(:,:) * dtime / 1000000.0_wp 
  !   formation_microbial%nitrogen15(:,:) = formation_microbial%nitrogen15(:,:) + & 
  !                                         (microbial_nue * uptake_dom(:,:) + & 
  !                                          microbial_uptake_nh4_n15_sl(:,:) + microbial_uptake_no3_n15_sl(:,:)) * &
  !                                         dtime / 1000000.0_wp
  ! 
  ! END SUBROUTINE calc_microbial_growth
  ! 
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! ! Sub Task called from update_sb_jsm
  ! !
  ! !-----------------------------------------------------------------------------------------------------
  ! !> calculate depolymerisation rate of polymeric litter and microbial residue 
  ! !! corresponds Ahrens et al. 2015, eq. 1 
  ! !! updated to the more recent version of JSM, including microbial residue as alternative polymeric pool
  ! !!
  ! !! Input: maximum depolymerisation rate (vmax), temperature and moisture modifiers for
  ! !!        depolymerisation and half-saturation, polymeric/residue pool and microbial biomass
  ! !!
  ! !! Output: transfer of residue/polymeric pool to DOM
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE calc_depolymerisation(nc, nsoil_sb, vmax_depolymerisation,rtm_depolymerisation,rmm_depolymerisation, &
  !                                  rtm_hsc,rmm_hsc,microbial_carbon,pool_polymeric,loss_polymeric,formation_dom)
  ! 
  !   USE mo_sb_constants,          ONLY: km_depolymerisation
  !   dsl4slm_Use_memory_a_sb
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   INTEGER,                                          INTENT(in)    :: nc, &                 !< dimensions
  !                                                                      nsoil_sb
  !   REAL(wp),                                         INTENT(in)    :: vmax_depolymerisation !< maximum depolymerisation rate
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_depolymerisation  !< temperature rate modifier for depolymerisation
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_depolymerisation  !< moisture rate modifier for depolymerisation
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_hsc               !< temperature rate modifier for half-saturation constants
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_hsc               !< moisture rate modifier for half-saturation constants
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: microbial_carbon      !< microbial biomass carbon [mol C/m3]
  !   TYPE(t_sb_element_a),                             INTENT(in)    :: pool_polymeric        !< polymeric/residue pool [mol/m3]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: loss_polymeric        !< loss of polymeric/residue pool [mol/m3/timestep]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: formation_dom         !< formation of DOM [mol/m3/timestep]
  !   ! ---------------------------
  !   ! 0.2 Local 
  !   REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1                   , &
  !                                                       hlp2                   , &
  !                                                       depolymerisation_rate
  !   CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_depolymerisation'
  !   ! ---------------------------
  !   ! 0.3 Declare Memory
  !   ! ...
  ! 
  ! 
  ! 
  !   !> 1.0 determine potential depolymerisation rate [mol / m3 / timestep]
  !   !!
  !   hlp1(:,:) = 1._wp   
  !   hlp1(:,:) = microbial_carbon(:,:) / (km_depolymerisation * rtm_hsc(:,:) * rmm_hsc(:,:) + microbial_carbon(:,:)) 
  !   hlp2(:,:) = vmax_depolymerisation * rtm_depolymerisation(:,:) * rmm_depolymerisation(:,:) * hlp1(:,:) * dtime * & 
  !                pool_polymeric%carbon(:,:) 
  ! 
  !   !> 2.0 convert to turnover rate [1 / timestep] 
  !   !!
  !   WHERE(pool_polymeric%carbon(:,:) > eps8)
  !       depolymerisation_rate(:,:) = hlp2(:,:) / pool_polymeric%carbon(:,:)
  !   ELSEWHERE
  !       depolymerisation_rate(:,:) = 0.0_wp
  !   ENDWHERE
  ! 
  !   !> 3.0 apply to all elements in this pool (should be a simple operator overload, but does currently not work safely)
  !   !!
  !   hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%carbon(:,:)
  !   loss_polymeric%carbon(:,:)      = loss_polymeric%carbon(:,:) + hlp2(:,:) 
  !   formation_dom%carbon(:,:)       = formation_dom%carbon(:,:) + hlp2(:,:) 
  ! 
  !   hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%nitrogen(:,:)
  !   loss_polymeric%nitrogen(:,:)    = loss_polymeric%nitrogen(:,:) + hlp2(:,:) 
  !   formation_dom%nitrogen(:,:)     = formation_dom%nitrogen(:,:) + hlp2(:,:) 
  ! 
  !   hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%phosphorus(:,:)
  !   loss_polymeric%phosphorus(:,:)  = loss_polymeric%phosphorus(:,:) + hlp2(:,:) 
  !   formation_dom%phosphorus(:,:)   = formation_dom%phosphorus(:,:) + hlp2(:,:) 
  ! 
  !   hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%carbon13(:,:)
  !   loss_polymeric%carbon13(:,:)    = loss_polymeric%carbon13(:,:) + hlp2(:,:) 
  !   formation_dom%carbon13(:,:)     = formation_dom%carbon13(:,:) + hlp2(:,:) 
  ! 
  !   hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%carbon14(:,:)
  !   loss_polymeric%carbon14(:,:)    = loss_polymeric%carbon14(:,:) + hlp2(:,:) 
  !   formation_dom%carbon14(:,:)     = formation_dom%carbon14(:,:) + hlp2(:,:) 
  ! 
  !   hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%nitrogen15(:,:)
  !   loss_polymeric%nitrogen15(:,:)  = loss_polymeric%nitrogen15(:,:) + hlp2(:,:) 
  !   formation_dom%nitrogen15(:,:)   = formation_dom%nitrogen15(:,:) + hlp2(:,:) 
  ! 
  ! END SUBROUTINE calc_depolymerisation
  ! 
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! ! Sub Task called from update_sb_pools
  ! !
  ! !-----------------------------------------------------------------------------------------------------
  ! !> calculate depolymerisation rate of polymeric litter and microbial residue 
  ! !! corresponds Ahrens et al. 2015, eq. 1 
  ! !! updated to the more recent version of JSM, including microbial residue as alternative polymeric pool
  ! !! enzyme partitioning simulated using SESAM approach, by Wutzler et al. 2017
  ! !!
  ! !! Input: maximum depolymerisation rate (vmax), temperature and moisture modifiers for
  ! !!        depolymerisation and half-saturation, polymeric/residue pool and microbial biomass
  ! !!
  ! !! Output: transfer of residue/polymeric pool to DOM
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE calc_depolymerisation_SESAM(nc, nsoil_sb,rtm_hsc,rmm_hsc,rtm_depolymerisation,rmm_depolymerisation, &
  !                                        microbial_carbon, microbial_cue_mavg, microbial_nue_mavg, microbial_pue_mavg, &
  !                                        pool_polymeric, pool_residue, &
  !                                        enzyme_frac_poly, enzyme_frac_residue, & 
  !                                        enzyme_frac_poly_c, enzyme_frac_poly_n, enzyme_frac_poly_p, &
  !                                        enzyme_frac_poly_c_mavg, enzyme_frac_poly_n_mavg, enzyme_frac_poly_p_mavg, &
  !                                        loss_polymeric, loss_residue, formation_dom)
  ! 
  !   USE mo_sb_constants
  !   USE mo_jsb_math_constants,              ONLY: dtime, one_day
  !   USE mo_quincy_update_time_averages,     ONLY: mavg_period_tenzyme
  !   dsl4slm_Use_memory_a_sb
  ! 
  !   IMPLICIT NONE
  ! 
  !   INTEGER,                                          INTENT(in)    :: nc, &                    !< dimensions
  !                                                                      nsoil_sb
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_hsc                  !< temperature rate modifier for half-saturation constants
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_hsc                  !< moisture rate modifier for half-saturation constants
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_depolymerisation     !< temperature rate modifier for depolymerisation
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_depolymerisation     !< moisture rate modifier for depolymerisation
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: microbial_carbon         !< microbial biomass carbon [mol C/m3]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: microbial_cue_mavg       !< microbial cue moving average
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: microbial_nue_mavg       !< microbial nue moving average
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: microbial_pue_mavg       !< microbial pue moving average
  !   TYPE(t_sb_element_a),                             INTENT(in)    :: pool_polymeric           !< polymeric pool [mol/m3]
  !   TYPE(t_sb_element_a),                             INTENT(in)    :: pool_residue             !< residue pool [mol/m3]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_poly         !< fraction of enzyme for poly litter depolymerisation
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_residue      !< fraction of enzyme for microbial residue depolymerisation
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_poly_c       !< temporal fraction of enzyme for poly litter depolymerisation, C
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_poly_n       !< temporal fraction of enzyme for poly litter depolymerisation, N
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_poly_p       !< temporal fraction of enzyme for poly litter depolymerisation, P
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_poly_c_mavg  !< mavg of enzyme for poly litter depolymerisation, C
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_poly_n_mavg  !< mavg of enzyme for poly litter depolymerisation, N
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: enzyme_frac_poly_p_mavg  !< mavg of enzyme for poly litter depolymerisation, P
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: loss_polymeric           !< loss of polymeric pool [mol/m3/timestep]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: loss_residue             !< loss of residue pool [mol/m3/timestep]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: formation_dom            !< formation of DOM [mol/m3/timestep]
  ! 
  !   REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1, hlp2, hlp3, &
  !                                                       poly_depolymerisation_rate, &
  !                                                       residue_depolymerisation_rate, &
  !                                                       depoly_c, depoly_n, depoly_p, hlp_depoly
  !   LOGICAL, DIMENSION(nsoil_sb)                     :: tmp1,tmp2,tmp3
  !   INTEGER                                          :: isoil
  !   CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_depolymerisation_SESAM'
  ! 
  !   ! ------------------------------------------------------------------------------------------------------------
  !   ! Go science
  !   ! ------------------------------------------------------------------------------------------------------------
  ! 
  !   !> 0.9 init local variables
  !   !!
  !   depoly_c(:,:) = 0.0_wp
  !   depoly_n(:,:) = 0.0_wp
  !   depoly_p(:,:) = 0.0_wp
  ! 
  ! 
  !   !> 1.0 calculate polymeric litter depolymerisation rate 
  !   !!
  !   hlp1(:,:) = 1._wp   
  !   hlp1(:,:) = microbial_carbon(:,:) * enzyme_frac_poly / &
  !              (km_depolymerisation * rtm_hsc(:,:) * rmm_hsc(:,:) + microbial_carbon(:,:) * enzyme_frac_poly) 
  !   hlp2(:,:) = vmax_depolymerisation_litter * rtm_depolymerisation * rmm_depolymerisation * hlp1(:,:) * dtime * & 
  !                pool_polymeric%carbon(:,:) 
  ! 
  !   WHERE(pool_polymeric%carbon(:,:) > eps8)
  !       poly_depolymerisation_rate(:,:) = hlp2(:,:) / pool_polymeric%carbon(:,:)
  !   ELSEWHERE
  !       poly_depolymerisation_rate(:,:) = 0.0_wp
  !   ENDWHERE
  ! 
  !   !> 2.0 calculate microbial residue depolymerisation rate 
  !   !!
  !   hlp1(:,:) = 1._wp   
  !   hlp1(:,:) = microbial_carbon(:,:) * enzyme_frac_residue / &
  !              (km_depolymerisation * rtm_hsc(:,:) * rmm_hsc(:,:) + microbial_carbon(:,:) * enzyme_frac_residue) 
  !   hlp2(:,:) = vmax_depolymerisation_residue * rtm_depolymerisation * rmm_depolymerisation * hlp1(:,:) * dtime * & 
  !                pool_residue%carbon(:,:) 
  ! 
  !   WHERE(pool_residue%carbon(:,:) > eps8)
  !       residue_depolymerisation_rate(:,:) = hlp2(:,:) / pool_residue%carbon(:,:)
  !   ELSEWHERE
  !       residue_depolymerisation_rate(:,:) = 0.0_wp
  !   ENDWHERE
  ! 
  !   !> 3.0 apply to all elements in polymeric litter pool (should be a simple operator overload, but does currently not work safely)
  !   !!
  !   hlp2(:,:)                       = poly_depolymerisation_rate(:,:) * pool_polymeric%carbon(:,:)
  !   loss_polymeric%carbon(:,:)      = loss_polymeric%carbon(:,:) + hlp2(:,:) 
  !   formation_dom%carbon(:,:)       = formation_dom%carbon(:,:) + hlp2(:,:) 
  !   depoly_c(:,:) = depoly_c(:,:) + hlp2(:,:)
  ! 
  !   hlp2(:,:)                       = poly_depolymerisation_rate(:,:) * pool_polymeric%nitrogen(:,:)
  !   loss_polymeric%nitrogen(:,:)    = loss_polymeric%nitrogen(:,:) + hlp2(:,:) 
  !   formation_dom%nitrogen(:,:)     = formation_dom%nitrogen(:,:) + hlp2(:,:) 
  !   depoly_n(:,:) = depoly_n(:,:) + hlp2(:,:)
  ! 
  !   hlp2(:,:)                       = poly_depolymerisation_rate(:,:) * pool_polymeric%phosphorus(:,:)
  !   loss_polymeric%phosphorus(:,:)  = loss_polymeric%phosphorus(:,:) + hlp2(:,:) 
  !   formation_dom%phosphorus(:,:)   = formation_dom%phosphorus(:,:) + hlp2(:,:) 
  !   depoly_p(:,:) = depoly_p(:,:) + hlp2(:,:)
  ! 
  !   hlp2(:,:)                       = poly_depolymerisation_rate(:,:) * pool_polymeric%carbon13(:,:)
  !   loss_polymeric%carbon13(:,:)    = loss_polymeric%carbon13(:,:) + hlp2(:,:) 
  !   formation_dom%carbon13(:,:)     = formation_dom%carbon13(:,:) + hlp2(:,:) 
  ! 
  !   hlp2(:,:)                       = poly_depolymerisation_rate(:,:) * pool_polymeric%carbon14(:,:)
  !   loss_polymeric%carbon14(:,:)    = loss_polymeric%carbon14(:,:) + hlp2(:,:) 
  !   formation_dom%carbon14(:,:)     = formation_dom%carbon14(:,:) + hlp2(:,:) 
  ! 
  !   hlp2(:,:)                       = poly_depolymerisation_rate(:,:) * pool_polymeric%nitrogen15(:,:)
  !   loss_polymeric%nitrogen15(:,:)  = loss_polymeric%nitrogen15(:,:) + hlp2(:,:) 
  !   formation_dom%nitrogen15(:,:)   = formation_dom%nitrogen15(:,:) + hlp2(:,:) 
  ! 
  !   !> 4.0 apply to all elements in microbial residue pool
  !   !!
  !   hlp2(:,:)                       = residue_depolymerisation_rate(:,:) * pool_residue%carbon(:,:)
  !   loss_residue%carbon(:,:)        = loss_residue%carbon(:,:) + hlp2(:,:) 
  !   formation_dom%carbon(:,:)       = formation_dom%carbon(:,:) + hlp2(:,:) 
  !   depoly_c(:,:) = depoly_c(:,:) + hlp2(:,:)
  ! 
  !   hlp2(:,:)                       = residue_depolymerisation_rate(:,:) * pool_residue%nitrogen(:,:)
  !   loss_residue%nitrogen(:,:)      = loss_residue%nitrogen(:,:) + hlp2(:,:) 
  !   formation_dom%nitrogen(:,:)     = formation_dom%nitrogen(:,:) + hlp2(:,:) 
  !   depoly_n(:,:) = depoly_n(:,:) + hlp2(:,:)
  ! 
  !   hlp2(:,:)                       = residue_depolymerisation_rate(:,:) * pool_residue%phosphorus(:,:)
  !   loss_residue%phosphorus(:,:)    = loss_residue%phosphorus(:,:) + hlp2(:,:) 
  !   formation_dom%phosphorus(:,:)   = formation_dom%phosphorus(:,:) + hlp2(:,:) 
  !   depoly_p(:,:) = depoly_p(:,:) + hlp2(:,:)
  ! 
  !   hlp2(:,:)                       = residue_depolymerisation_rate(:,:) * pool_residue%carbon13(:,:)
  !   loss_residue%carbon13(:,:)      = loss_residue%carbon13(:,:) + hlp2(:,:) 
  !   formation_dom%carbon13(:,:)     = formation_dom%carbon13(:,:) + hlp2(:,:) 
  ! 
  !   hlp2(:,:)                       = residue_depolymerisation_rate(:,:) * pool_residue%carbon14(:,:)
  !   loss_residue%carbon14(:,:)      = loss_residue%carbon14(:,:) + hlp2(:,:) 
  !   formation_dom%carbon14(:,:)     = formation_dom%carbon14(:,:) + hlp2(:,:) 
  ! 
  !   hlp2(:,:)                       = residue_depolymerisation_rate(:,:) * pool_residue%nitrogen15(:,:)
  !   loss_residue%nitrogen15(:,:)    = loss_residue%nitrogen15(:,:) + hlp2(:,:) 
  !   formation_dom%nitrogen15(:,:)   = formation_dom%nitrogen15(:,:) + hlp2(:,:) 
  ! 
  !   !> 5.0 calculate the temporal enzyme partitioning for poly litter based on C, N, and P/m3/s
  !   !!
  !   WHERE (poly_depolymerisation_rate > 0.0_wp .AND. poly_depolymerisation_rate> 0.0_wp)
  !     !! alpha_C
  !     hlp1(:,:)               = poly_depolymerisation_rate(:,:) * pool_polymeric%carbon(:,:)
  !     hlp2(:,:)               = residue_depolymerisation_rate(:,:) * pool_residue%carbon(:,:)
  !     hlp3(:,:)               = km_depolymerisation * rtm_hsc(:,:) * rmm_hsc(:,:)
  !     enzyme_frac_poly_c(:,:) = calc_enzyme_allocation(hlp1(:,:), hlp2(:,:), hlp3(:,:), microbial_carbon(:,:))
  !     !! alpha_N
  !     hlp1(:,:)               = poly_depolymerisation_rate(:,:) * pool_polymeric%nitrogen(:,:)
  !     hlp2(:,:)               = residue_depolymerisation_rate(:,:) * pool_residue%nitrogen(:,:)
  !     enzyme_frac_poly_n(:,:) = calc_enzyme_allocation(hlp1(:,:), hlp2(:,:), hlp3(:,:), microbial_carbon(:,:))
  !     !! alpha_P
  !     hlp1(:,:)               = poly_depolymerisation_rate(:,:) * pool_polymeric%phosphorus(:,:)
  !     hlp2(:,:)               = residue_depolymerisation_rate(:,:) * pool_residue%phosphorus(:,:)
  !     enzyme_frac_poly_p(:,:) = calc_enzyme_allocation(hlp1(:,:), hlp2(:,:), hlp3(:,:), microbial_carbon(:,:))
  !   ENDWHERE
  ! 
  !   !> 6.0 calculate the enzyme partitioning based on depolymerized C, N, P
  !   !!
  !   hlp1(:,:) = depoly_c(:,:)*microbial_cue_mavg
  !   hlp2(:,:) = depoly_n(:,:)*microbial_cn*microbial_nue_mavg
  !   hlp3(:,:) = depoly_p(:,:)*microbial_cn*microbial_np*microbial_pue_mavg
  ! 
  !   WHERE (hlp1(:,:)>0.0_wp .AND. hlp2(:,:)>0.0_wp .AND. hlp3(:,:)>0.0_wp)
  !     WHERE (hlp1(:,:) < hlp2(:,:) .AND. hlp1(:,:)<hlp3(:,:))
  !       hlp_depoly(:,:) = enzyme_frac_poly_c_mavg(:,:)
  !     ELSEWHERE (hlp2(:,:) < hlp1(:,:) .AND. hlp2(:,:)<hlp3(:,:))
  !       hlp_depoly(:,:) = enzyme_frac_poly_n_mavg(:,:)
  !     ELSEWHERE (hlp3(:,:) < hlp1(:,:) .AND. hlp3(:,:)<hlp2(:,:))
  !       hlp_depoly(:,:) = enzyme_frac_poly_p_mavg(:,:)
  !     ENDWHERE
  !     enzyme_frac_poly(:,:)    = MIN( MAX(enzyme_frac_poly(:,:) + (hlp_depoly(:,:) - enzyme_frac_poly(:,:))&
  !                                / mavg_period_tenzyme / one_day * dtime, min_enzyme_fraction), (1._wp - min_enzyme_fraction))
  !     enzyme_frac_residue(:,:) = 1._wp - enzyme_frac_poly(:,:)
  !   ENDWHERE
  ! END SUBROUTINE calc_depolymerisation_SESAM
  ! 

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
  ELEMENTAL FUNCTION calc_bulk_density_correction(bulk_soil_carbon_sl, &
                                                  soil_litter_carbon_sl, &
                                                  bulk_dens_sl ) &
                                                  RESULT(bulk_dens_corr_sl)
  
    USE mo_jsb_physical_constants,    ONLY: molar_mass_C
    USE mo_sb_constants,              ONLY: rho_bulk_org, carbon_per_dryweight_SOM 
  
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in)             :: bulk_soil_carbon_sl     !< soil non-litter carbon contributing to density [mol C / m3] 
    REAL(wp), INTENT(in)             :: soil_litter_carbon_sl   !< soil litter carbon contributing to density [mol C / m3] 
    REAL(wp), INTENT(in)             :: bulk_dens_sl            !< bulk density of the mineral soil [kg / m3]
    REAL(wp)                         :: bulk_dens_corr_sl       !< combined bulk density of the mineral and organic soil [kg / m3]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                    :: om_dens_sl       ! organic matter density [kg / m3] 
    REAL(wp)                    :: lit_dens_sl      ! litter density [kg / m3] 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bulk_density_correction'
  
  
  
    !> 1.0 calculate organic matter density and litter density per soil layer [kg / m3]
    !!
    om_dens_sl  = bulk_soil_carbon_sl * molar_mass_C / carbon_per_dryweight_SOM / 1000._wp 
    lit_dens_sl = soil_litter_carbon_sl * molar_mass_C / carbon_per_dryweight_SOM / 1000._wp 
  
    !> 2.0 actual bulk density of the soil 
    !!
  
    bulk_dens_corr_sl = MAX(rho_bulk_org, om_dens_sl + lit_dens_sl + bulk_dens_sl - & 
                            om_dens_sl * bulk_dens_sl / rho_bulk_org - &
                            lit_dens_sl * bulk_dens_sl / rho_bulk_org) 
  
  
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
  ELEMENTAL SUBROUTINE calc_qmax_bulk_density_correction(bulk_soil_carbon_sl, &
                                                         volume_min_sl, bulk_dens_sl, &
                                                         qmax_om, qmax_min, &
                                                         qmax, qmax1, qmax2)
  
    USE mo_jsb_physical_constants,    ONLY: molar_mass_C
    USE mo_sb_constants,              ONLY: carbon_per_dryweight_SOM, rho_bulk_org
  
    REAL(wp),           INTENT(in)  :: bulk_soil_carbon_sl   !< non-litter organic carbon content of the soil layer [mol C / m3]
    REAL(wp),           INTENT(in)  :: volume_min_sl         !< mineral soil volume fraction the soil layer [m3/ m3]
    REAL(wp),           INTENT(in)  :: bulk_dens_sl          !< bulk density of the mineral soil [kg / m3]
    REAL(wp),           INTENT(in)  :: qmax_om               !< maximium sorption capacity on OM [mol / kg]
    REAL(wp),           INTENT(in)  :: qmax_min              !< maximium sorption capacity on  mineral soil [mol / kg]
    REAL(wp),           INTENT(out) :: qmax                  !< maximum sorption capacity of the soil [mol / m3]
    REAL(wp), OPTIONAL, INTENT(out) :: qmax1                 !< maximum sorption capacity of 1st pool [mol / m3]
    REAL(wp), OPTIONAL, INTENT(out) :: qmax2                 !< maximum sorption capacity of 2nd pool [mol / m3]
  
    REAL(wp)                    :: vol_fraction_om
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_qmax_bulk_density_correction'
  
  
    ! calculate volume fraction filled by organic material 
    ! given assumed density and carbon content of OM and litter
    vol_fraction_om = MIN(bulk_soil_carbon_sl * molar_mass_C / 1000._wp / carbon_per_dryweight_SOM / rho_bulk_org, 1.0_wp)
  
    ! maximum sorption capacity per soil volume as the mix of OM and mineral sorption (mol/m3)
    qmax = vol_fraction_om * qmax_om * rho_bulk_org + volume_min_sl * qmax_min * bulk_dens_sl
    IF (PRESENT(qmax1)) THEN
      qmax1 = vol_fraction_om * qmax_om * rho_bulk_org
    ENDIF
    IF (PRESENT(qmax2)) THEN
      qmax2 = volume_min_sl * qmax_min * bulk_dens_sl
    ENDIF
  
  END SUBROUTINE calc_qmax_bulk_density_correction
  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from calc_Psorption_parameter
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates the km for phosphate of a soil layer
  !! first version by Lin, 2019-08-26
  !!
  !! Input: default km, solute concentration, qmax of different sorption pool (mineral and organic combined) 
  !!
  !! Output: km of the combined sorption pool 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_km_qmax_correction(solute, qmax1, qmax2, km1, km2, w_soil) RESULT(km)
  
    IMPLICIT NONE
  
    REAL(wp),           INTENT(in) :: solute          !< po4 concentration of the soil layer [mol P / m3]
    REAL(wp),           INTENT(in) :: qmax1           !< maximum po4 sorption capacity of the 1st pool [mol P/ m3]
    REAL(wp),           INTENT(in) :: qmax2           !< maximum po4 sorption capacity of the 2nd pool [mol P/ m3]
    REAL(wp),           INTENT(in) :: km1             !< km of the 1st pool [mol P/ m3]
    REAL(wp),           INTENT(in) :: km2             !< km of the 2nd pool [mol P/ m3]
    REAL(wp), OPTIONAL, INTENT(in) :: w_soil          !< water content of the soil layer [L / m3]
    REAL(wp)                       :: km              !< km of the combined pool [mol P/ m3]
  
    REAL(wp)                    :: dq1, dq2, hlp, sol1, sol2
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_km_qmax_correction'
  
  
    !> 1.0 calculate the change of sorbed P pool given the solute concentration
    dq1  = qmax1 * km1 / (km1 + solute) / (km1 + solute)
    dq2  = qmax2 * km2 / (km2 + solute) / (km2 + solute)
    hlp  = dq1 + dq2
  
    !> 2.0 calculate km for the combined P pool
    !! 
    ! Qmax*km/(km+solute)^2 = dq1 + dq2, hlp = dq1 + dq2
    ! km = (Qmax/hlp - 2*solute +/- sqrt((Qmax/hlp)^2-4 * Qmax/hlp))/2
  
    IF (hlp < eps8) THEN
      km = MAX(km1, km2)
    ELSEIF ((km1 + solute) < eps8 .OR. (km2 + solute) < eps8) THEN
      km = MAX(km1, km2)
    ELSE 
      hlp  = (qmax1 + qmax2) / hlp
      sol1 = (hlp - 2._wp * solute + SQRT(hlp * hlp - 4._wp * hlp * solute)) / 2._wp
      sol2 = (hlp - 2._wp * solute - SQRT(hlp * hlp - 4._wp * hlp * solute)) / 2._wp
      IF (sol1 /= sol1 .OR. sol2 /= sol2) THEN
        km = MAX(km1, km2)
      ELSEIF ((km1 - sol1) * (km2 - sol1) <= 0._wp) THEN
        km = sol1
      ELSE
        km = sol2
      ENDIF
    ENDIF
  
    IF (PRESENT(w_soil)) km = km * w_soil
  
  END FUNCTION calc_km_qmax_correction
  
  ! ------------------------------------------------------------------------------------------------------- !
  !>calculate bulk soil carbon used for the calculation of soil bulk density and sorption capacity 
  !>
  SUBROUTINE calc_bulk_soil_carbon( &
    & nc, &
    & nsoil_sb, &
    & sb_model_scheme, &
    & sb_pool_microbial_carbon, &
    & sb_pool_residue_carbon, &
    & sb_pool_residue_assoc_carbon, &
    & sb_pool_dom_carbon, &
    & sb_pool_dom_assoc_carbon, &
    & sb_pool_soluable_litter_carbon, &
    & sb_pool_polymeric_litter_carbon, &
    & sb_pool_woody_litter_carbon, &
    & bulk_soil_carbon_sl, &
    & soil_litter_carbon_sl, &
    & volume_min_sl &
    & )
  
    USE mo_jsb_physical_constants,    ONLY: molar_mass_C
    USE mo_sb_constants,              ONLY: rho_bulk_org, carbon_per_dryweight_SOM
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                          INTENT(in)        :: nc                               !< chunk size
    INTEGER,                          INTENT(in)        :: nsoil_sb                         !< soil-layer number
    CHARACTER(len=*)               ,  INTENT(in)        :: sb_model_scheme                  !< soil model scheme: simple_1d or jsm
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)        :: sb_pool_microbial_carbon         !< sb_pool [mol/m3]
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)        :: sb_pool_residue_carbon           !< sb_pool [mol/m3]
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)        :: sb_pool_residue_assoc_carbon     !< sb_pool [mol/m3]
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)        :: sb_pool_dom_carbon               !< sb_pool [mol/m3]
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)        :: sb_pool_dom_assoc_carbon         !< sb_pool [mol/m3]
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)        :: sb_pool_soluable_litter_carbon   !< sb_pool [mol/m3]
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)        :: sb_pool_polymeric_litter_carbon  !< sb_pool [mol/m3]
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)        :: sb_pool_woody_litter_carbon      !< sb_pool [mol/m3]
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(inout)     :: bulk_soil_carbon_sl              !< mol/m3
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(inout)     :: soil_litter_carbon_sl            !< mol/m3
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(inout)     :: volume_min_sl                    !< m3/m3
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                             :: isoil
    REAL(wp), DIMENSION(nc,nsoil_sb)                    :: om_dens_sl          ! organic matter density [kg / m3] 
    REAL(wp), DIMENSION(nc,nsoil_sb)                    :: lit_dens_sl         ! litter density [kg / m3] 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bulk_soil_carbon'
  
    !>1.0 calculate the content of organic mattter [mol/m3] and litter [mol/m3]
    !>
    ! soil layer 1
    ! woody_litter is not considered for layer 1
    bulk_soil_carbon_sl(:,1) = &
      & sb_pool_microbial_carbon(:,1) + &
      & sb_pool_residue_carbon(:,1) + &
      & sb_pool_residue_assoc_carbon(:,1) + &
      & sb_pool_dom_carbon(:,1) + &
      & sb_pool_dom_assoc_carbon(:,1)     
    soil_litter_carbon_sl(:,1) = &
      & sb_pool_soluable_litter_carbon(:,1) + &
      & sb_pool_polymeric_litter_carbon(:,1)
    ! soil layers 2:x
    DO isoil = 2,nsoil_sb
      bulk_soil_carbon_sl(:,isoil) = &
        & sb_pool_microbial_carbon(:,isoil) + &
        & sb_pool_residue_carbon(:,isoil) + &
        & sb_pool_residue_assoc_carbon(:,isoil) + &
        & sb_pool_dom_carbon(:,isoil) + &
        & sb_pool_dom_assoc_carbon(:,isoil)
      soil_litter_carbon_sl(:,isoil) = &
        & sb_pool_soluable_litter_carbon(:,isoil) + &
        & sb_pool_polymeric_litter_carbon(:,isoil) + &
        & sb_pool_woody_litter_carbon(:,isoil)
    ENDDO
  
    !>  1.1 simple_1d
    !>
    ! only jsm distinguish litter C from SOM C
    IF (sb_model_scheme == 'simple_1d') THEN
      bulk_soil_carbon_sl(:,:)   = bulk_soil_carbon_sl(:,:) + soil_litter_carbon_sl(:,:)
      soil_litter_carbon_sl(:,:) = 0._wp
    ENDIF
  
    !>2.0 calculate the volume of mineral soil [m3/m3]
    !> 
    ! assuming the bulk soil is composed of OM, litter and soil mineral
    om_dens_sl(:,:)  = bulk_soil_carbon_sl(:,:) * molar_mass_C / carbon_per_dryweight_SOM / 1000._wp    ! OM density [kg/m3]
    lit_dens_sl(:,:) = soil_litter_carbon_sl(:,:) * molar_mass_C / carbon_per_dryweight_SOM / 1000._wp  ! litter density [kg/m3]
    DO isoil = 1,nsoil_sb
      volume_min_sl(:,isoil) = MAX(0._wp, 1._wp - om_dens_sl(:,isoil) / rho_bulk_org - lit_dens_sl(:,isoil) / rho_bulk_org)
    ENDDO
  
  END SUBROUTINE calc_bulk_soil_carbon
  
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
  ELEMENTAL FUNCTION calc_fast_po4(Labile_p, Qmax, km) RESULT (fast_po4)
  
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: Labile_p                    !< labile P pool of the soil layer [mol P / m3]
    REAL(wp), INTENT(in) :: Qmax                        !< maximium sorption capacity of po4 [mol P/ m3]
    REAL(wp), INTENT(in) :: km                          !< sorption affinity of po4 [mol P/ m3]
    REAL(wp)             :: fast_po4                    !< fast associated po4 of the soil layer [mol P / m3]
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fast_po4'
  
  
    !> 1.0 
    !! based on the assumption that labile P is the sum of both,
    !! and they follow the Langmuir equilibrium
    !! Labile_p - fast_po4 = fast_po4 * km / (Labile_p * (Qmax - fast_po4))
    fast_po4 = (Labile_p + Qmax + km - sqrt((Labile_p + Qmax + km)**2._wp - 4._wp * Labile_p * Qmax)) / 2._wp
  
  END FUNCTION calc_fast_po4
  
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! ! Sub Task called from update_sb_pools 
  ! !
  ! !-----------------------------------------------------------------------------------------------------
  ! !> calculates potential enzyme allocation coefficient based on SESAM concept, by Thomas Wutzler
  ! !!
  ! !! Input: potential depolymerisation rate (C, N, P), microbial biomass, half-saturation coefficient
  ! !!
  ! !! Output: potential enzyme coefficient for C, N, P
  ! !-----------------------------------------------------------------------------------------------------
  ! ELEMENTAL FUNCTION calc_enzyme_allocation(dL, dR, km, Cmic) RESULT (alpha)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   REAL(wp), INTENT(in) :: dL          !< organic carbon content of the soil layer [mol C / m3]
  !   REAL(wp), INTENT(in) :: dR          !< bulk density of the mineral soil [kg / m3]
  !   REAL(wp), INTENT(in) :: km          !< maximium sorption capacity on OM [mol / kg]
  !   REAL(wp), INTENT(in) :: Cmic        !< maximium sorption capacity on  mineral soil [mol / kg]
  !   REAL(wp)             :: alpha       !< maximum sorption capacity of the soil [mol / m3]
  !   ! ---------------------------
  !   ! 0.2 Local
  !   REAL(wp)                    :: hlp1    !< b**2-4*a*c in the quadratic equation
  !   CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_enzyme_allocation'
  ! 
  !   !> Solution of the equation:
  !   !! dL/dR=alpha(km+alpha*Cmic)/((1-alpha)*(km+Cmic*(1-alpha)))
  !   !! only one of the solution fits in the range [0,1]
  !   hlp1  = 4._wp * Cmic ** 2._wp * dL * dR + &
  !           8._wp * Cmic * km * dL * dR + &
  !           km ** 2._wp * dL ** 2._wp + &
  !           2._wp * km ** 2._wp * dL * dR + &
  !           km ** 2._wp * dR ** 2._wp
  !   alpha = (dL * km + 2._wp * dL * Cmic + dR * km - sqrt(hlp1)) / (2._wp * Cmic * (dL - dR))
  ! 
  ! END FUNCTION calc_enzyme_allocation
  ! 
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! ! Sub Task called from update_sb_jsm
  ! !
  ! !-----------------------------------------------------------------------------------------------------
  ! !> calculates sorption and desorption of organic material
  ! !! corresponds to the sorption and desorption terms in Ahrens et al. 2015, eq. 5 & 6
  ! !!
  ! !! Input: sorption temperature and moisture modifiers, sorption pools
  ! !!
  ! !! Output: sorption and desorption rate (mol/m2/timestep) 
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE calc_om_sorption(nc, nsoil_sb, qmax_org, rtm_sorption, rmm_sorption, rtm_desorption, &
  !                             rmm_desorption, pool_dom, pool_dom_assoc, pool_residue, pool_residue_assoc, &
  !                             formation_dom, formation_dom_assoc, formation_residue, &
  !                             formation_residue_assoc, loss_dom, loss_dom_assoc, loss_residue, loss_residue_assoc)
  ! 
  !   USE mo_sb_constants,          ONLY: k_adsorpt_dom, k_desorpt_dom, k_adsorpt_det, k_desorpt_det 
  !   dsl4slm_Use_memory_a_sb
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   INTEGER,                                          INTENT(in)    :: nc, &                   !< dimensions
  !                                                                      nsoil_sb
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: qmax_org 
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_sorption            !< temperature rate modifier [-]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_sorption            !< moisture rate modifier [-]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rtm_desorption          !< temperature rate modifier [-]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: rmm_desorption          !< moisture rate modifier [-]
  !   TYPE(t_sb_element_a),                             INTENT(in)    :: pool_dom                !< DOM pool [mol/m3]
  !   TYPE(t_sb_element_a),                             INTENT(in)    :: pool_residue            !< DOM pool [mol/m3]
  !   TYPE(t_sb_element_a),                             INTENT(in)    :: pool_dom_assoc          !< DOM pool [mol/m3]
  !   TYPE(t_sb_element_a),                             INTENT(in)    :: pool_residue_assoc      !< DOM pool [mol/m3]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: formation_dom           !< rate of DOM transport [mol/m3/timestep]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: formation_dom_assoc     !< rate of DOM transport [mol/m3/timestep]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: formation_residue       !< rate of DOM transport [mol/m3/timestep]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: formation_residue_assoc !< rate of DOM transport [mol/m3/timestep]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: loss_dom                !< rate of DOM transport [mol/m3/timestep]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: loss_dom_assoc          !< rate of DOM transport [mol/m3/timestep]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: loss_residue            !< rate of DOM transport [mol/m3/timestep]
  !   TYPE(t_sb_element_a),                             INTENT(inout) :: loss_residue_assoc      !< rate of DOM transport [mol/m3/timestep]
  !   ! ---------------------------
  !   ! 0.2 Local
  !   REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1                   , &
  !                                                       avail_sorption_sites   , &
  !                                                       sorption_rate          , &
  !                                                       desorption_rate        
  !   CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_om_sorption'
  !   ! ---------------------------
  !   ! 0.3 Declare Memory
  !   ! ...
  ! 
  ! 
  !   !> 1.0 Determine available sorption sites [mol C m-3]
  !   !!
  !   avail_sorption_sites(:,:) = MAX(0._wp, qmax_org(:,:) - (pool_dom_assoc%carbon(:,:) + pool_residue_assoc%carbon(:,:)))
  ! 
  !   !> 2.0 mineral association of DOM
  !   sorption_rate(:,:)   = k_adsorpt_dom * rtm_sorption(:,:)   * rmm_sorption(:,:) * avail_sorption_sites(:,:) * dtime 
  !   desorption_rate(:,:) = k_desorpt_dom * rtm_desorption(:,:) * rmm_desorption(:,:) * dtime 
  ! 
  !   hlp1(:,:) = sorption_rate(:,:) * pool_dom%carbon(:,:)
  !   formation_dom_assoc%carbon(:,:) = formation_dom_assoc%carbon(:,:) + hlp1(:,:) 
  !   loss_dom%carbon(:,:)            = loss_dom%carbon(:,:) + hlp1(:,:) 
  !   hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%carbon(:,:)
  !   formation_dom%carbon(:,:)       = formation_dom%carbon(:,:) + hlp1(:,:) 
  !   loss_dom_assoc%carbon(:,:)      = loss_dom_assoc%carbon(:,:) + hlp1(:,:) 
  ! 
  !   hlp1(:,:) = sorption_rate(:,:) * pool_dom%nitrogen(:,:)
  !   formation_dom_assoc%nitrogen(:,:) = formation_dom_assoc%nitrogen(:,:) + hlp1(:,:) 
  !   loss_dom%nitrogen(:,:)            = loss_dom%nitrogen(:,:) + hlp1(:,:) 
  !   hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%nitrogen(:,:)
  !   formation_dom%nitrogen(:,:)       = formation_dom%nitrogen(:,:) + hlp1(:,:) 
  !   loss_dom_assoc%nitrogen(:,:)      = loss_dom_assoc%nitrogen(:,:) + hlp1(:,:) 
  ! 
  !   hlp1(:,:) = sorption_rate(:,:) * pool_dom%phosphorus(:,:)
  !   formation_dom_assoc%phosphorus(:,:) = formation_dom_assoc%phosphorus(:,:) + hlp1(:,:) 
  !   loss_dom%phosphorus(:,:)            = loss_dom%phosphorus(:,:) + hlp1(:,:) 
  !   hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%phosphorus(:,:)
  !   formation_dom%phosphorus(:,:)       = formation_dom%phosphorus(:,:) + hlp1(:,:) 
  !   loss_dom_assoc%phosphorus(:,:)      = loss_dom_assoc%phosphorus(:,:) + hlp1(:,:) 
  ! 
  !   hlp1(:,:) = sorption_rate(:,:) * pool_dom%carbon13(:,:)
  !   formation_dom_assoc%carbon13(:,:) = formation_dom_assoc%carbon13(:,:) + hlp1(:,:) 
  !   loss_dom%carbon13(:,:)            = loss_dom%carbon13(:,:) + hlp1(:,:) 
  !   hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%carbon13(:,:)
  !   formation_dom%carbon13(:,:)       = formation_dom%carbon13(:,:) + hlp1(:,:) 
  !   loss_dom_assoc%carbon13(:,:)      = loss_dom_assoc%carbon13(:,:) + hlp1(:,:) 
  ! 
  !   hlp1(:,:) = sorption_rate(:,:) * pool_dom%carbon14(:,:)
  !   formation_dom_assoc%carbon14(:,:) = formation_dom_assoc%carbon14(:,:) + hlp1(:,:) 
  !   loss_dom%carbon14(:,:)            = loss_dom%carbon14(:,:) + hlp1(:,:) 
  !   hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%carbon14(:,:)
  !   formation_dom%carbon14(:,:)       = formation_dom%carbon14(:,:) + hlp1(:,:) 
  !   loss_dom_assoc%carbon14(:,:)      = loss_dom_assoc%carbon14(:,:) + hlp1(:,:) 
  ! 
  !   hlp1(:,:) = sorption_rate(:,:) * pool_dom%nitrogen15(:,:)
  !   formation_dom_assoc%nitrogen15(:,:) = formation_dom_assoc%nitrogen15(:,:) + hlp1(:,:) 
  !   loss_dom%nitrogen15(:,:)            = loss_dom%nitrogen15(:,:) + hlp1(:,:) 
  !   hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%nitrogen15(:,:)
  !   formation_dom%nitrogen15(:,:)       = formation_dom%nitrogen15(:,:) + hlp1(:,:) 
  !   loss_dom_assoc%nitrogen15(:,:)      = loss_dom_assoc%nitrogen15(:,:) + hlp1(:,:) 
  ! 
  !   !> 3.0 mineral association of microbial residue 
  !   sorption_rate(:,:)   = k_adsorpt_det * rtm_sorption(:,:)   * rmm_sorption(:,:) * avail_sorption_sites(:,:) * dtime 
  !   desorption_rate(:,:) = k_desorpt_det * rtm_desorption(:,:) * rmm_desorption(:,:) * dtime 
  ! 
  !   hlp1(:,:) = sorption_rate(:,:) * pool_residue%carbon(:,:)
  !   formation_residue_assoc%carbon(:,:) = formation_residue_assoc%carbon(:,:) + hlp1(:,:) 
  !   loss_residue%carbon(:,:)            = loss_residue%carbon(:,:) + hlp1(:,:) 
  !   hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%carbon(:,:)
  !   formation_residue%carbon(:,:)       = formation_residue%carbon(:,:) + hlp1(:,:) 
  !   loss_residue_assoc%carbon(:,:)      = loss_residue_assoc%carbon(:,:) + hlp1(:,:) 
  ! 
  !   hlp1(:,:) = sorption_rate(:,:) * pool_residue%nitrogen(:,:)
  !   formation_residue_assoc%nitrogen(:,:) = formation_residue_assoc%nitrogen(:,:) + hlp1(:,:) 
  !   loss_residue%nitrogen(:,:)            = loss_residue%nitrogen(:,:) + hlp1(:,:) 
  !   hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%nitrogen(:,:)
  !   formation_residue%nitrogen(:,:)       = formation_residue%nitrogen(:,:) + hlp1(:,:) 
  !   loss_residue_assoc%nitrogen(:,:)      = loss_residue_assoc%nitrogen(:,:) + hlp1(:,:) 
  ! 
  !   hlp1(:,:) = sorption_rate(:,:) * pool_residue%phosphorus(:,:)
  !   formation_residue_assoc%phosphorus(:,:) = formation_residue_assoc%phosphorus(:,:) + hlp1(:,:) 
  !   loss_residue%phosphorus(:,:)            = loss_residue%phosphorus(:,:) + hlp1(:,:) 
  !   hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%phosphorus(:,:)
  !   formation_residue%phosphorus(:,:)       = formation_residue%phosphorus(:,:) + hlp1(:,:) 
  !   loss_residue_assoc%phosphorus(:,:)      = loss_residue_assoc%phosphorus(:,:) + hlp1(:,:) 
  ! 
  !   hlp1(:,:) = sorption_rate(:,:) * pool_residue%carbon13(:,:)
  !   formation_residue_assoc%carbon13(:,:) = formation_residue_assoc%carbon13(:,:) + hlp1(:,:) 
  !   loss_residue%carbon13(:,:)            = loss_residue%carbon13(:,:) + hlp1(:,:) 
  !   hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%carbon13(:,:)
  !   formation_residue%carbon13(:,:)       = formation_residue%carbon13(:,:) + hlp1(:,:) 
  !   loss_residue_assoc%carbon13(:,:)      = loss_residue_assoc%carbon13(:,:) + hlp1(:,:) 
  ! 
  !   hlp1(:,:) = sorption_rate(:,:) * pool_residue%carbon14(:,:)
  !   formation_residue_assoc%carbon14(:,:) = formation_residue_assoc%carbon14(:,:) + hlp1(:,:) 
  !   loss_residue%carbon14(:,:)            = loss_residue%carbon14(:,:) + hlp1(:,:) 
  !   hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%carbon14(:,:)
  !   formation_residue%carbon14(:,:)       = formation_residue%carbon14(:,:) + hlp1(:,:) 
  !   loss_residue_assoc%carbon14(:,:)      = loss_residue_assoc%carbon14(:,:) + hlp1(:,:) 
  ! 
  !   hlp1(:,:) = sorption_rate(:,:) * pool_residue%nitrogen15(:,:)
  !   formation_residue_assoc%nitrogen15(:,:) = formation_residue_assoc%nitrogen15(:,:) + hlp1(:,:) 
  !   loss_residue%nitrogen15(:,:)            = loss_residue%nitrogen15(:,:) + hlp1(:,:) 
  !   hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%nitrogen15(:,:)
  !   formation_residue%nitrogen15(:,:)       = formation_residue%nitrogen15(:,:) + hlp1(:,:) 
  !   loss_residue_assoc%nitrogen15(:,:)      = loss_residue_assoc%nitrogen15(:,:) + hlp1(:,:) 
  ! 
  ! END SUBROUTINE calc_om_sorption 

  ! ======================================================================================================= !
  !>calculate P sorption parameters, including both single and double Langmuir
  !>
  ELEMENTAL SUBROUTINE calc_Psorption_parameter( &
    & soil_depth_sl, &
    & flag_sb_double_langmuir, &
    & clay_sl, silt_sl, soil_water, ph_sl, & 
    & bulk_soil_carbon_sl, soil_litter_carbon_sl, woody_litter_carbon_sl, &
    & bulk_dens_sl, volume_min_sl, po4_solute, &
    & qmax_po4_min_sl, qmax_po4_om_sl, qmax_po4, km_adsorpt_po4_sl, &             ! inout
    & qmax_fast_po4, qmax_slow_po4, km_fast_po4, km_slow_po4, Qmax_AlFe_cor)      ! inout

    USE mo_sb_constants,        ONLY: qmax_po4_silt, qmax_po4_sand, qmax_po4_clay, qmax_po4_OM, qmax_po4_mineral, &
      &                               km_adsorpt_silt_po4, km_adsorpt_clay_po4, km_adsorpt_sand_po4, km_adsorpt_OM_po4, &
      &                               km_adsorpt_mineral_po4, km_po4_ph, min_fraction_default_sb
    USE mo_spq_util,            ONLY: calc_qmax_texture
    USE mo_jsb_math_constants,  ONLY: eps8

    REAL(wp),     INTENT(in)    :: soil_depth_sl                !< soil layer width    
    LOGICAL,      INTENT(in)    :: flag_sb_double_langmuir      !< conventional or double Langmuir
    REAL(wp),     INTENT(in)    :: clay_sl                      !< clay content of soil layer [-]
    REAL(wp),     INTENT(in)    :: silt_sl                      !< silt content of soil layer [-]
    REAL(wp),     INTENT(in)    :: soil_water                   !< soil water content of soil layer [L/m3]
    REAL(wp),     INTENT(in)    :: ph_sl                        !< ph of soil layer [-]
    REAL(wp),     INTENT(in)    :: bulk_soil_carbon_sl          !< bulk soil carbon content of soil layer [mol C/m3]
    REAL(wp),     INTENT(in)    :: soil_litter_carbon_sl        !< soil litter carbon content of soil layer [mol C/m3]
    REAL(wp),     INTENT(in)    :: woody_litter_carbon_sl       !< woody litter carbon content of soil layer [mol C/m3]
    REAL(wp),     INTENT(in)    :: bulk_dens_sl                 !< soil bulk density [kg/m3]
    REAL(wp),     INTENT(in)    :: volume_min_sl                !< soil mineral volumetric content [m3/m3]
    REAL(wp),     INTENT(in)    :: po4_solute                   !< soluble P concentration of soil layer [mol P/m3]
    REAL(wp),     INTENT(inout) :: qmax_po4_min_sl              !< maximium po4 sorption capacity for mineral soil [mol P/m3]
    REAL(wp),     INTENT(inout) :: qmax_po4_om_sl               !< maximium po4 sorption capacity for OM soil [mol P/m3]
    REAL(wp),     INTENT(inout) :: qmax_po4                     !< maximium po4 sorption capacity for soil layer [mol P/m3]
    REAL(wp),     INTENT(inout) :: km_adsorpt_po4_sl            !< half-saturation po4 concentration [mol P/m3]
    REAL(wp),     INTENT(inout) :: qmax_fast_po4                !< maximium po4 sorption capacity of fast associated pool [mol P/m3]
    REAL(wp),     INTENT(inout) :: qmax_slow_po4                !< maximium po4 sorption capacity of slow associated pool [mol P/m3]
    REAL(wp),     INTENT(inout) :: km_fast_po4                  !< half-saturation po4 concentration of fast associated pool [mol P/m3]
    REAL(wp),     INTENT(inout) :: km_slow_po4                  !< half-saturation po4 concentration of slow associated pool [mol P/m3]
    REAL(wp),     INTENT(inout) :: Qmax_AlFe_cor                !< qmax correlation factor of Al/Fe content [-]

    REAL(wp)  :: hlp1, hlp2, hlp3, hlp4
    REAL(wp)  :: sand_sl
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_Psorption_parameter'
    

    ! run only for soil layers with width larger zero
    !
    IF (soil_depth_sl > eps8) THEN

      !>0.8 initialize local variable
      !>
      sand_sl = 1._wp - clay_sl - silt_sl

      !>0.9 set zero if lower eps8
      !>
      IF (Qmax_AlFe_cor < eps8) THEN
        Qmax_AlFe_cor = 0._wp
      END IF

      !>1.0 double-surface Langmuir
      !>
      IF (flag_sb_double_langmuir) THEN
        !>  1.1 calculate Qmax based on soil texture
        !> 
        ! fine soil + sandy soil [mol P /kg soil]
        qmax_po4_min_sl = (calc_qmax_texture(qmax_po4_silt, clay_sl, silt_sl) + & 
          & calc_qmax_texture(qmax_po4_sand, 0._wp, sand_sl)) * Qmax_AlFe_cor 
        ! soil OM [mol P /kg OM]
        qmax_po4_om_sl = qmax_po4_OM
        CALL calc_qmax_bulk_density_correction( &
          & bulk_soil_carbon_sl, &
          & volume_min_sl, &
          & bulk_dens_sl, & 
          & qmax_po4_om_sl, &
          & qmax_po4_min_sl, &
          & qmax_po4, &
          & qmax_fast_po4, &
          & qmax_slow_po4)
        ! Fe/Al in OM [mol P/m3]
        hlp1 = calc_qmax_texture(qmax_po4_clay, 0._wp, clay_sl) * volume_min_sl * bulk_dens_sl * Qmax_AlFe_cor
        ! update qmax_fast_po4 and qmax_po4
        qmax_fast_po4 = qmax_fast_po4 + hlp1
        qmax_po4      = qmax_po4 + hlp1

        !>  1.2 calculate km based on soil texture
        !>
        ! km for the OM-Fe/Al pool, fast pool
        hlp2 = qmax_fast_po4 - hlp1                        ! qmax for non-adsorbed Fe/Al [mol P/m3]
        hlp3 = km_po4_ph * ph_sl * km_adsorpt_clay_po4     ! km for Fe/Al, mol P/m3
        hlp4 = km_po4_ph * ph_sl * km_adsorpt_OM_po4       ! km for OM, mol P/m3
        km_fast_po4 = calc_km_qmax_correction(po4_solute, hlp1, hlp2, hlp3, hlp4, soil_water)
        ! km for the clay-silt-sand pool, slow pool
        ! catch "qmax_po4_min_sl < eps8" to avoid division by zero
        IF (qmax_slow_po4 < eps8 .OR. qmax_po4_min_sl < eps8) THEN
          qmax_slow_po4 = 0._wp
          hlp3        = km_po4_ph * ph_sl * km_adsorpt_silt_po4    ! km for clay/silt [mol P/L]
          hlp4        = km_po4_ph * ph_sl * km_adsorpt_sand_po4    ! km for sand [mol P/L]
          km_slow_po4 = MAX(hlp3, hlp4)
        ELSE
          hlp1        = qmax_slow_po4 * MIN(calc_qmax_texture(qmax_po4_silt, clay_sl, silt_sl) * & 
            & Qmax_AlFe_cor / qmax_po4_min_sl, (1._wp - min_fraction_default_sb)) ! qmax for clay/silt
          hlp2        = qmax_slow_po4 - hlp1                                      ! qmax for sand
          hlp3        = km_po4_ph * ph_sl * km_adsorpt_silt_po4                   ! km for clay/silt [mol P/L]
          hlp4        = km_po4_ph * ph_sl * km_adsorpt_sand_po4                   ! km for sand [mol P/L]
          km_slow_po4 = calc_km_qmax_correction(po4_solute, hlp1, hlp2, hlp3, hlp4, soil_water)
        END IF
        ! km for the whole soil layer
        km_adsorpt_po4_sl = calc_km_qmax_correction( &
          & po4_solute, qmax_fast_po4,qmax_slow_po4, km_fast_po4, km_slow_po4)

      !>2.0 single-surface Langmuir
      !> 
      ! Qmax and km are based on texture and calculated both as the sum of OM and mineral soil
      ELSE
        qmax_po4_min_sl   = calc_qmax_texture(qmax_po4_mineral, 0.0_wp, 1._wp)
        qmax_po4_om_sl    = qmax_po4_min_sl * min_fraction_default_sb
        CALL calc_qmax_bulk_density_correction( &
          & bulk_soil_carbon_sl, &    ! in 
          & volume_min_sl, &          ! in 
          & bulk_dens_sl, &           ! in 
          & qmax_po4_om_sl, &         ! in 
          & qmax_po4_min_sl, &        ! in 
          & qmax_po4, &               ! out
          & qmax_fast_po4, &          ! out
          & qmax_slow_po4)            ! out
        qmax_po4          = qmax_po4 * soil_water / 1000._wp ! corrected by water content
        CALL calc_qmax_bulk_density_correction( &
          & bulk_soil_carbon_sl, &    ! in 
          & volume_min_sl, &          ! in 
          & bulk_dens_sl, &           ! in 
          & km_adsorpt_OM_po4, &      ! in 
          & km_adsorpt_mineral_po4, & ! in 
          & km_adsorpt_po4_sl, &      ! out
          & km_fast_po4, &            ! out
          & km_slow_po4)              ! out
      END IF  ! double/single Langmuir

    END IF ! IF (soil_depth_sl > eps8)
           
  END SUBROUTINE calc_Psorption_parameter

#endif
END MODULE mo_sb_jsm_processes
