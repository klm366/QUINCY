!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., and Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Contains the science routines for vertical transport along the soil column 
!!
!! @author
!!  Sönke Zaehle 
!!
!! @par Revision History
!! First version by Sönke Zaehle (20170124): adopted from JSM_interface and JSM_core
!!
MODULE mo_sb_jsm_transport

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: message, finish
  USE mo_jsb_math_constants,    ONLY: dtime


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_particle_fluxrate, calc_bioturbation_rate, &
            calc_bioturbation_transport_wrapper, calc_liquid_phase_transport_wrapper , & 
            calc_particle_transport_wrapper, & 
            calc_bioturbation_transport, calc_liquid_phase_transport

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_jsm_processes'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> Interface to the liquid transport routine, owing to the lack of a suitable operator overload 
  !!
  !! Input: leachable pools percolation rate, and soil depth 
  !!
  !! Output: vertical transport rate (mol/m2/timestep) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_liquid_phase_transport_wrapper(nc, nsoil_sb, percolation_sl,soil_depth_sl, &
                                         pool_dom, nh4_solute, no3_solute, po4_solute, nh4_n15_solute, no3_n15_solute, & 
                                         transport_dom, transport_nh4_solute, transport_no3_solute, & 
                                         transport_po4_solute, transport_nh4_n15_solute, transport_no3_n15_solute, &
                                         leaching_dom_carbon, leaching_dom_nitrogen, leaching_dom_phosphorus, & 
                                         leaching_dom_carbon13, leaching_dom_carbon14, leaching_dom_nitrogen15, &
                                         leaching_nh4_solute, leaching_no3_solute, leaching_po4_solute, &
                                         leaching_nh4_n15_solute, leaching_no3_n15_solute)
                                         
    dsl4slm_Use_memory_a_sb
                                         
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                  !< dimensions
                                                                       nsoil_sb
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: percolation_sl         !< fraction of pool transported up or down per second
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: soil_depth_sl          !< depth of each soil layer [m] 
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_dom               !< DOM pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: nh4_solute             !< NH4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: no3_solute             !< NO3 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: po4_solute             !< PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: nh4_n15_solute         !< N15H4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: no3_n15_solute         !< N15O3 pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(inout) :: transport_dom          !< rate of DOM transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_nh4_solute   !< rate of NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_no3_solute   !< rate of NO3 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_po4_solute   !< rate of PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_nh4_n15_solute !< rate of N15H4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_no3_n15_solute !< rate of N15O3 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: leaching_dom_carbon      !< rate of DOM-C leaching [µmol/m2/s]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: leaching_dom_nitrogen    !< rate of DOM-N leaching [µmol/m2/s]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: leaching_dom_phosphorus  !< rate of DOM-P leaching [µmol/m2/s]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: leaching_dom_carbon13    !< rate of DOM-C13 leaching [µmol/m2/s]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: leaching_dom_carbon14    !< rate of DOM-C14 leaching [µmol/m2/s]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: leaching_dom_nitrogen15  !< rate of DOM-N15 leaching [µmol/m2/s]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: leaching_nh4_solute      !< rate of NH4 leaching [µmol/m2/s]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: leaching_no3_solute      !< rate of NO3 leaching [µmol/m2/s]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: leaching_po4_solute      !< rate of PO4 leaching [µmol/m2/s]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: leaching_nh4_n15_solute  !< rate of 14NH4 leaching [µmol/m2/s]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: leaching_no3_n15_solute  !< rate of 15NO3 leaching [µmol/m2/s]
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_liquid_phase_transport_wrapper'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 transport of DOM
    !!
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), soil_depth_sl(:,:), &
                                        pool_dom%carbon(:,:),transport_dom%carbon(:,:), &
                                        leaching_dom_carbon(:)) 
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), soil_depth_sl(:,:), &
                                        pool_dom%nitrogen(:,:),transport_dom%nitrogen(:,:), &
                                        leaching_dom_nitrogen(:)) 
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), soil_depth_sl(:,:), &
                                        pool_dom%phosphorus(:,:),transport_dom%phosphorus(:,:), & 
                                        leaching_dom_phosphorus(:)) 
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), soil_depth_sl(:,:), &
                                        pool_dom%carbon13(:,:),transport_dom%carbon13(:,:), &
                                        leaching_dom_carbon13(:)) 
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), soil_depth_sl(:,:), &
                                        pool_dom%carbon14(:,:),transport_dom%carbon14(:,:), &
                                        leaching_dom_carbon14(:)) 
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), soil_depth_sl(:,:), &
                                        pool_dom%nitrogen15(:,:),transport_dom%nitrogen15(:,:), & 
                                        leaching_dom_nitrogen15(:)) 

    !> 2.0 transport of other solutes (NH4, NO3, PO4, and their isotopes)
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), soil_depth_sl(:,:), &
                                        nh4_solute(:,:),transport_nh4_solute(:,:), &
                                        leaching_nh4_solute(:)) 
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), soil_depth_sl(:,:), &
                                        no3_solute(:,:),transport_no3_solute(:,:), &
                                        leaching_no3_solute(:)) 
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), soil_depth_sl(:,:), &
                                        po4_solute(:,:),transport_po4_solute(:,:), &
                                        leaching_po4_solute(:)) 
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), soil_depth_sl(:,:), &
                                        nh4_n15_solute(:,:),transport_nh4_n15_solute(:,:), &
                                        leaching_nh4_n15_solute(:)) 
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), soil_depth_sl(:,:), &
                                        no3_n15_solute(:,:),transport_no3_n15_solute(:,:), & 
                                        leaching_no3_n15_solute(:)) 

  END SUBROUTINE calc_liquid_phase_transport_wrapper


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates matter transport by bioturbation, by calling calc_bioturbation_transport
  !! partly necessary because of a missing operator overload...
  !!
  !! Input:  pools subject to bioturbation, percolation rate, and soil depth 
  !!         This explicitly excludes mycorrhiza, as these die when eaten 
  !!
  !! Output: vertical transport rate (mol/m2/timestep) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_bioturbation_transport_wrapper(nc, nsoil_sb, k_bioturb, soil_depth_sl, &
                                          pool_polymeric_litter, pool_dom_assoc, & 
                                          pool_fungi, pool_microbial, pool_residue, pool_residue_assoc, & 
                                          nh4_assoc, po4_assoc_fast, po4_assoc_slow, & 
                                          po4_occluded, nh4_n15_assoc, & 
                                          transport_polymeric_litter, transport_dom_assoc, & 
                                          transport_fungi, transport_microbial, transport_residue, & 
                                          transport_residue_assoc, transport_nh4_assoc, & 
                                          transport_po4_assoc_fast, transport_po4_assoc_slow, & 
                                          transport_po4_occluded, transport_nh4_n15_assoc)

    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                  !< dimensions
                                                                       nsoil_sb
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: k_bioturb              !< diffusion factor for bioturbation 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: soil_depth_sl          !< depth of each soil layer [m] 
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_polymeric_litter  !< polymeric litter pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_dom_assoc         !< minerally associated DOM pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_fungi             !< fungi pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_microbial         !< microbial pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_residue           !< microbial residue pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_residue_assoc     !< minerally associated microbial residue pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: nh4_assoc              !< minerally associated NH4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: po4_assoc_fast         !< fast minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: po4_assoc_slow         !< slow minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: po4_occluded           !< occluded PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: nh4_n15_assoc          !< N15H4 pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(inout) :: transport_polymeric_litter !< rate of polymeric litter transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: transport_dom_assoc        !< rate of minerally associated DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: transport_fungi            !< rate of fungi transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: transport_microbial        !< rate of microbial transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: transport_residue          !< rate of microbial residue transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: transport_residue_assoc    !< rate of minerally associated microbial residue transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_nh4_assoc        !< rate of NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_po4_assoc_fast   !< rate of fast associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_po4_assoc_slow   !< rate of slow associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_po4_occluded     !< rate of occluded PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_nh4_n15_assoc    !< rate of N15H4 transport [mol/m3/timestep]
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bioturbation_transport_wrapper'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 transport of organic pools
    !! 
    !> 1.1 polymeric litter
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_polymeric_litter%carbon(:,:), transport_polymeric_litter%carbon(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_polymeric_litter%nitrogen(:,:), transport_polymeric_litter%nitrogen(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_polymeric_litter%phosphorus(:,:), transport_polymeric_litter%phosphorus(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_polymeric_litter%carbon13(:,:), transport_polymeric_litter%carbon13(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_polymeric_litter%carbon14(:,:), transport_polymeric_litter%carbon14(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_polymeric_litter%nitrogen15(:,:), transport_polymeric_litter%nitrogen15(:,:))

    !> 1.2 minerally associated DOM 
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_dom_assoc%carbon(:,:), transport_dom_assoc%carbon(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_dom_assoc%nitrogen(:,:), transport_dom_assoc%nitrogen(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_dom_assoc%phosphorus(:,:), transport_dom_assoc%phosphorus(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_dom_assoc%carbon13(:,:), transport_dom_assoc%carbon13(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_dom_assoc%carbon14(:,:), transport_dom_assoc%carbon14(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_dom_assoc%nitrogen15(:,:), transport_dom_assoc%nitrogen15(:,:))

    !> 1.3 fungi biomass 
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_fungi%carbon(:,:), transport_fungi%carbon(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_fungi%nitrogen(:,:), transport_fungi%nitrogen(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_fungi%phosphorus(:,:), transport_fungi%phosphorus(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_fungi%carbon13(:,:), transport_fungi%carbon13(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_fungi%carbon14(:,:), transport_fungi%carbon14(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_fungi%nitrogen15(:,:), transport_fungi%nitrogen15(:,:))

    !> 1.4 microbial biomass 
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_microbial%carbon(:,:), transport_microbial%carbon(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_microbial%nitrogen(:,:), transport_microbial%nitrogen(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_microbial%phosphorus(:,:), transport_microbial%phosphorus(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_microbial%carbon13(:,:), transport_microbial%carbon13(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_microbial%carbon14(:,:), transport_microbial%carbon14(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_microbial%nitrogen15(:,:), transport_microbial%nitrogen15(:,:))

    !> 1.5 microbial residue 
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_residue%carbon(:,:), transport_residue%carbon(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_residue%nitrogen(:,:), transport_residue%nitrogen(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_residue%phosphorus(:,:), transport_residue%phosphorus(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_residue%carbon13(:,:), transport_residue%carbon13(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_residue%carbon14(:,:), transport_residue%carbon14(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_residue%nitrogen15(:,:), transport_residue%nitrogen15(:,:))

    !> 1.6 minerally associated microbial residue 
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_residue_assoc%carbon(:,:), transport_residue_assoc%carbon(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_residue_assoc%nitrogen(:,:), transport_residue_assoc%nitrogen(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_residue_assoc%phosphorus(:,:), transport_residue_assoc%phosphorus(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_residue_assoc%carbon13(:,:), transport_residue_assoc%carbon13(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_residue_assoc%carbon14(:,:), transport_residue_assoc%carbon14(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:), & 
                                        pool_residue_assoc%nitrogen15(:,:), transport_residue_assoc%nitrogen15(:,:))

    !> 2.0 inorganic pools 
    ! CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:),nh4_assoc(:,:),transport_nh4_assoc(:,:)) 
    ! CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:),po4_assoc_fast(:,:),transport_po4_assoc_fast(:,:)) 
    ! CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:),po4_assoc_slow(:,:),transport_po4_assoc_slow(:,:)) 
    ! CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:),po4_occluded(:,:),transport_po4_occluded(:,:)) 
    ! CALL calc_bioturbation_transport(nc, nsoil_sb, k_bioturb(:,:),soil_depth_sl(:,:),nh4_n15_assoc(:,:),transport_nh4_n15_assoc(:,:)) 

  END SUBROUTINE calc_bioturbation_transport_wrapper


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates particile flux rate given the change in organic pools through formation 
  !! and loss
  !! corresponds to the sum of the omega term in Ahrens et al. 2015, eq. S1-S3
  !!
  !! Input: change of soil pools through formation and loss
  !!
  !! Output: omega_transport_term 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_particle_fluxrate(nc, nsoil_sb, soil_depth_sl,bulk_dens_corr_sl,sb_formation,sb_loss) RESULT (particle_fluxrate)

    USE mo_sb_constants,          ONLY: carbon_per_dryweight_SOM, rho_bulk_org
    USE mo_jsb_phys_constants,    ONLY: molar_mass_C
    USE mo_jsb_math_constants,    ONLY: eps8
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in) :: nc, &                  !< dimensions
                                                                    nsoil_sb
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in) :: soil_depth_sl          !< depth of soil layer [m] 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in) :: bulk_dens_corr_sl      !< bulk density of mineral+organic soil [kg/m3]
    TYPE(t_sb_pool_a),                                INTENT(in) :: sb_formation           !< formation of pools [mol/m3/timestep]
    TYPE(t_sb_pool_a),                                INTENT(in) :: sb_loss                !< loss from pools [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb)                            :: particle_fluxrate 
    ! ---------------------------
    ! 0.2 Local
    INTEGER                     :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_particle_fluxrate'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 calculate net change of organic matter cummulative down the profile (kg/m2/dtime)
    !! note the omission of woody litter in the first layer, which is assumed not to be part of the profile
    !!
    particle_fluxrate(:,1) = soil_depth_sl(:,1) * &
           (sb_formation%soluable_litter%carbon(:,1)  - sb_loss%soluable_litter%carbon(:,1) + &
           sb_formation%polymeric_litter%carbon(:,1) - sb_loss%polymeric_litter%carbon(:,1) + &
           sb_formation%dom%carbon(:,1)              - sb_loss%dom%carbon(:,1) + &
           sb_formation%dom_assoc%carbon(:,1)        - sb_loss%dom_assoc%carbon(:,1) + &
           sb_formation%fungi%carbon(:,1)            - sb_loss%fungi%carbon(:,1) + &
           sb_formation%mycorrhiza%carbon(:,1)       - sb_loss%mycorrhiza%carbon(:,1) + &
           sb_formation%microbial%carbon(:,1)        - sb_loss%microbial%carbon(:,1) + &
           sb_formation%residue%carbon(:,1)          - sb_loss%residue%carbon(:,1) + &
           sb_formation%residue_assoc%carbon(:,1)    - sb_loss%residue_assoc%carbon(:,1)) / &
           (carbon_per_dryweight_SOM / molar_mass_C * 1000._wp)

    DO isoil = 2,nsoil_sb
       particle_fluxrate(:,isoil) = soil_depth_sl(:,isoil) * &
          (sb_formation%soluable_litter%carbon(:,isoil)  - sb_loss%soluable_litter%carbon(:,isoil) + &
           sb_formation%polymeric_litter%carbon(:,isoil) - sb_loss%polymeric_litter%carbon(:,isoil) + &
           sb_formation%woody_litter%carbon(:,isoil)     - sb_loss%woody_litter%carbon(:,isoil) + &
           sb_formation%dom%carbon(:,isoil)              - sb_loss%dom%carbon(:,isoil) + &
           sb_formation%dom_assoc%carbon(:,isoil)        - sb_loss%dom_assoc%carbon(:,isoil) + &
           sb_formation%fungi%carbon(:,isoil)            - sb_loss%fungi%carbon(:,isoil) + &
           sb_formation%mycorrhiza%carbon(:,isoil)       - sb_loss%mycorrhiza%carbon(:,isoil) + &
           sb_formation%microbial%carbon(:,isoil)        - sb_loss%microbial%carbon(:,isoil) + &
           sb_formation%residue%carbon(:,isoil)          - sb_loss%residue%carbon(:,isoil) + &
           sb_formation%residue_assoc%carbon(:,isoil)    - sb_loss%residue_assoc%carbon(:,isoil)) / &
           (carbon_per_dryweight_SOM / molar_mass_C * 1000._wp)
   ENDDO

   !> 2.0 convert change in pools in density change and thus transport term (m per dtime)
   !!
   WHERE(bulk_dens_corr_sl(:,:) > eps8)
      particle_fluxrate(:,:) = particle_fluxrate(:,:) / rho_bulk_org
   ENDWHERE

  END FUNCTION calc_particle_fluxrate


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates matter transport due to particle flux, by calling calc_particle_transport
  !! partly necessary because of a missing operator overload...
  !!
  !! Input:  pools subject to particle flux, and soil depth 
  !!
  !! Output: vertical transport rate (mol/m2/timestep) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_particle_transport_wrapper(nc, nsoil_sb, particle_fluxrate, soil_depth_sl, &
                                          sb_pool, nh4_solute, nh4_assoc, no3_solute, po4_solute, po4_assoc_fast, & 
                                          po4_assoc_slow, po4_occluded, po4_primary, nh4_n15_solute, nh4_n15_assoc, no3_n15_solute, & 
                                          sb_transport, transport_nh4_solute, transport_nh4_assoc, transport_no3_solute, & 
                                          transport_po4_solute, transport_po4_assoc_fast, transport_po4_assoc_slow, & 
                                          transport_po4_occluded, transport_po4_primary, transport_nh4_n15_solute, &
                                          transport_nh4_n15_assoc, transport_no3_n15_solute)
                                          
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                      !< dimensions
                                                                       nsoil_sb
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: particle_fluxrate          !< particle fluxrate (1/timestep) 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: soil_depth_sl              !< depth of each soil layer [m] 
    TYPE(t_sb_pool_a),                                INTENT(in)    :: sb_pool                    !< organic pools [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: nh4_solute                 !< NH4 pool in soil solution [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: nh4_assoc                  !< minerally associated NH4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: no3_solute                 !< NO3 pool in soil solution [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: po4_solute                 !< PO4 pool in soil solution [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: po4_assoc_fast             !< fast minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: po4_assoc_slow             !< slow minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: po4_occluded               !< occluded PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: po4_primary                !< occluded PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: nh4_n15_solute             !< N15H4 pool in soil solution [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: nh4_n15_assoc              !< minerally associated N15H4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: no3_n15_solute             !< N15O3 pool in soil solution [mol/m3]
    TYPE(t_sb_pool_a),                                INTENT(inout) :: sb_transport               !< xminerally associated DOM pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_nh4_solute       !< rate of solute NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_nh4_assoc        !< rate of associated NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_no3_solute       !< rate of solute NO3 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_po4_solute       !< rate of solute PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_po4_assoc_fast   !< rate of associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_po4_assoc_slow   !< rate of associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_po4_occluded     !< rate of associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_po4_primary      !< rate of associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_nh4_n15_solute   !< rate of solute N15H4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_no3_n15_solute   !< rate of solute N15O3 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport_nh4_n15_assoc    !< rate of associated N15H4 transport [mol/m3/timestep]
    ! ---------------------------
    ! 0.2 Local
    INTEGER                                          :: isoil
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: particle_fluxrate_woody
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_particle_transport_wrapper'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 transport of organic pools
    !! 
    !> 1.1 soluable litter
    CALL calc_particle_transport(nc, nsoil_sb,particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%soluable_litter%carbon(:,:), sb_transport%soluable_litter%carbon(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%soluable_litter%nitrogen(:,:), sb_transport%soluable_litter%nitrogen(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%soluable_litter%phosphorus(:,:), sb_transport%soluable_litter%phosphorus(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%soluable_litter%carbon13(:,:), sb_transport%soluable_litter%carbon13(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%soluable_litter%carbon14(:,:), sb_transport%soluable_litter%carbon14(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%soluable_litter%nitrogen15(:,:), sb_transport%soluable_litter%nitrogen15(:,:))

    !> 1.2 polymeric litter
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%polymeric_litter%carbon(:,:), sb_transport%polymeric_litter%carbon(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%polymeric_litter%nitrogen(:,:), sb_transport%polymeric_litter%nitrogen(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%polymeric_litter%phosphorus(:,:), sb_transport%polymeric_litter%phosphorus(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%polymeric_litter%carbon13(:,:), sb_transport%polymeric_litter%carbon13(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%polymeric_litter%carbon14(:,:), sb_transport%polymeric_litter%carbon14(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%polymeric_litter%nitrogen15(:,:), sb_transport%polymeric_litter%nitrogen15(:,:))

    !> 1.3 woody litter
    !! first layer of woody litter does not belong to profile 
    DO isoil = 1,nsoil_sb
      particle_fluxrate_woody(:,isoil) = particle_fluxrate(:,isoil) - particle_fluxrate(:,1)
    ENDDO
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate_woody(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%woody_litter%carbon(:,:), sb_transport%woody_litter%carbon(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate_woody(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%woody_litter%nitrogen(:,:), sb_transport%woody_litter%nitrogen(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate_woody(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%woody_litter%phosphorus(:,:), sb_transport%woody_litter%phosphorus(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate_woody(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%woody_litter%carbon13(:,:), sb_transport%woody_litter%carbon13(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate_woody(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%woody_litter%carbon14(:,:), sb_transport%woody_litter%carbon14(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate_woody(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%woody_litter%nitrogen15(:,:), sb_transport%woody_litter%nitrogen15(:,:))

    !> 1.5 DOM 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%dom%carbon(:,:), sb_transport%dom%carbon(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%dom%nitrogen(:,:), sb_transport%dom%nitrogen(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%dom%phosphorus(:,:), sb_transport%dom%phosphorus(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%dom%carbon13(:,:), sb_transport%dom%carbon13(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%dom%carbon14(:,:), sb_transport%dom%carbon14(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%dom%nitrogen15(:,:), sb_transport%dom%nitrogen15(:,:))

    !> 1.6 minerally associated DOM 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%dom_assoc%carbon(:,:), sb_transport%dom_assoc%carbon(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%dom_assoc%nitrogen(:,:), sb_transport%dom_assoc%nitrogen(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%dom_assoc%phosphorus(:,:), sb_transport%dom_assoc%phosphorus(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%dom_assoc%carbon13(:,:), sb_transport%dom_assoc%carbon13(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%dom_assoc%carbon14(:,:), sb_transport%dom_assoc%carbon14(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%dom_assoc%nitrogen15(:,:), sb_transport%dom_assoc%nitrogen15(:,:))

    !> 1.7 fungi biomass 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%fungi%carbon(:,:), sb_transport%fungi%carbon(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%fungi%nitrogen(:,:), sb_transport%fungi%nitrogen(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%fungi%phosphorus(:,:), sb_transport%fungi%phosphorus(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%fungi%carbon13(:,:), sb_transport%fungi%carbon13(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%fungi%carbon14(:,:), sb_transport%fungi%carbon14(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%fungi%nitrogen15(:,:), sb_transport%fungi%nitrogen15(:,:))

    !> 1.8 mycorrhiza biomass 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%mycorrhiza%carbon(:,:), sb_transport%mycorrhiza%carbon(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%mycorrhiza%nitrogen(:,:), sb_transport%mycorrhiza%nitrogen(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%mycorrhiza%phosphorus(:,:), sb_transport%mycorrhiza%phosphorus(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%mycorrhiza%carbon13(:,:), sb_transport%mycorrhiza%carbon13(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%mycorrhiza%carbon14(:,:), sb_transport%mycorrhiza%carbon14(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%mycorrhiza%nitrogen15(:,:), sb_transport%mycorrhiza%nitrogen15(:,:))

    !> 1.9 microbial biomass 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%microbial%carbon(:,:), sb_transport%microbial%carbon(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%microbial%nitrogen(:,:), sb_transport%microbial%nitrogen(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%microbial%phosphorus(:,:), sb_transport%microbial%phosphorus(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%microbial%carbon13(:,:), sb_transport%microbial%carbon13(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%microbial%carbon14(:,:), sb_transport%microbial%carbon14(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%microbial%nitrogen15(:,:), sb_transport%microbial%nitrogen15(:,:))

    !> 1.10 microbial residue 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%residue%carbon(:,:), sb_transport%residue%carbon(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%residue%nitrogen(:,:), sb_transport%residue%nitrogen(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%residue%phosphorus(:,:), sb_transport%residue%phosphorus(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%residue%carbon13(:,:), sb_transport%residue%carbon13(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%residue%carbon14(:,:), sb_transport%residue%carbon14(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%residue%nitrogen15(:,:), sb_transport%residue%nitrogen15(:,:))

    !> 1.11 minerally associated microbial residue 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%residue_assoc%carbon(:,:), sb_transport%residue_assoc%carbon(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%residue_assoc%nitrogen(:,:), sb_transport%residue_assoc%nitrogen(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%residue_assoc%phosphorus(:,:), sb_transport%residue_assoc%phosphorus(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%residue_assoc%carbon13(:,:), sb_transport%residue_assoc%carbon13(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%residue_assoc%carbon14(:,:), sb_transport%residue_assoc%carbon14(:,:))
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:),soil_depth_sl(:,:), & 
                                        sb_pool%residue_assoc%nitrogen15(:,:), sb_transport%residue_assoc%nitrogen15(:,:))

    !> 2.0 inorganic pools 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:), soil_depth_sl(:,:), &
                                 nh4_solute(:,:),transport_nh4_solute(:,:)) 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:), soil_depth_sl(:,:), &
                                 nh4_assoc(:,:),transport_nh4_assoc(:,:)) 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:), soil_depth_sl(:,:), &
                                 no3_solute(:,:),transport_no3_solute(:,:)) 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:), soil_depth_sl(:,:), &
                                 po4_solute(:,:),transport_po4_solute(:,:)) 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:), soil_depth_sl(:,:), &
                                 po4_assoc_fast(:,:),transport_po4_assoc_fast(:,:)) 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:), soil_depth_sl(:,:), &
                                 po4_assoc_slow(:,:),transport_po4_assoc_slow(:,:)) 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:), soil_depth_sl(:,:), &
                                 po4_occluded(:,:),transport_po4_occluded(:,:)) 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:), soil_depth_sl(:,:), &
                                 po4_primary(:,:),transport_po4_primary(:,:)) 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:), soil_depth_sl(:,:), &
                                 nh4_n15_solute(:,:),transport_nh4_n15_solute(:,:)) 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:), soil_depth_sl(:,:), &
                                 nh4_n15_assoc(:,:),transport_nh4_n15_assoc(:,:)) 
    CALL calc_particle_transport(nc, nsoil_sb, particle_fluxrate(:,:), soil_depth_sl(:,:), &
                                 no3_n15_solute(:,:),transport_no3_n15_solute(:,:)) 

  END SUBROUTINE calc_particle_transport_wrapper


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from calc_particle_transport 
  !
  !-----------------------------------------------------------------------------------------------------
  !> Advective Transport in the vertical dimension 
  !! corresponds to the particule flux transport term in Ahrens et al. 2015, eq 4-7 
  !!
  !! Input: pool to be transported, percolation rate, and soil depth 
  !!
  !! Output: vertical transport rate (mol/m2/s) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_particle_transport(nc, nsoil_sb, particle_fluxrate,soil_depth_sl,pool,transport)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                  !< dimensions
                                                                       nsoil_sb
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: particle_fluxrate      !< fraction of pool transported up or down per timestep 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: soil_depth_sl          !< depth of each soil layer [m] 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: pool                   !< current state of pool to be transported
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport              !< transport between layers per timestep 
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)         :: hlp1 
    INTEGER                         :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_particle_transport'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 advective transport in the vertical axis of the soil (i.e. down)
    !! I can do 1,n-1 here, because by definition, there cannot be any upward transport from the uppermost layer
    DO isoil = 1,nsoil_sb-1
       WHERE(particle_fluxrate(:,isoil) > 0.0_wp)
          hlp1(:) = particle_fluxrate(:,isoil) * pool(:,isoil)
       ELSEWHERE
          hlp1(:) = particle_fluxrate(:,isoil) * pool(:,isoil+1)
       ENDWHERE
       transport(:,isoil)   = transport(:,isoil)   - hlp1(:) / soil_depth_sl(:,isoil)
       transport(:,isoil+1) = transport(:,isoil+1) + hlp1(:) / soil_depth_sl(:,isoil+1)
       
   ENDDO

   !> 2.0 lower boundary condition,??
   !!  based on the assumption that there is no upward flux, and that all other 
   !!  material is lost
!   WHERE(particle_fluxrate(:,nsoil_sb) > 0.0_wp)
!      hlp1(:) = particle_fluxrate(:,nsoil_sb) * pool(:,nsoil_sb) * soil_depth_sl(:,nsoil_sb)
!   ELSEWHERE
      hlp1(:) = 0.0_wp 
!   ENDWHERE
   transport(:,nsoil_sb) = transport(:,nsoil_sb) - hlp1(:) / soil_depth_sl(:,nsoil_sb)

  END SUBROUTINE calc_particle_transport


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from calc_liquid_phase_transport_wrapper
  !
  !-----------------------------------------------------------------------------------------------------
  !> Advective Transport in the vertical dimension 
  !! corresponds to the water flux transport term in Ahrens et al. 2015, eq 5 
  !!
  !! Input: pool to be transported, percolation rate, and soil depth 
  !!
  !! Output: vertical transport rate (mol/m2/s) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl,soil_depth_sl,pool,transport,leaching)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &              !< dimensions
                                                                       nsoil_sb
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: percolation_sl     !< fraction of pool transported up or down per second
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: soil_depth_sl      !< depth of each soil layer [m] 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: pool               !< current state of pool to be transported [mol m-3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport          !< transport between layers per timestep [mol m-3]
    REAL(wp), DIMENSION(nc),                          INTENT(out)   :: leaching           !< transport through the lower surface of the lowest 
                                                                                          !! layer per timestep [mol m-3 timestep]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)         :: hlp1 
    INTEGER                         :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_liquid_phase_transport'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 advective transport in the vertical axis of the soil (i.e. down)
    !! I can do 1,n-1 here, because by definition, there cannot be any upward transport from the uppermost layer
    DO isoil = 1,nsoil_sb-1
       WHERE(percolation_sl(:,isoil) > 0.0_wp)
          hlp1(:) = percolation_sl(:,isoil) * dtime/1000._wp * pool(:,isoil) * soil_depth_sl(:,isoil)
!! BLARPP for now no upward transport - need to properly think about this in update_sse_physics
       ELSEWHERE
          hlp1(:) = 0.0_wp ! percolation_sl(:,isoil) * dtime/1000._wp * pool(:,isoil+1) * soil_depth_sl(:,isoil+1)
       ENDWHERE
       transport(:,isoil)   = transport(:,isoil)   - hlp1(:) / soil_depth_sl(:,isoil)
       transport(:,isoil+1) = transport(:,isoil+1) + hlp1(:) / soil_depth_sl(:,isoil+1)
    ENDDO

    !> 2.0 lower boundary condition,
    !!  based on the assumption that there is no upward flux from below, 
    !!  and that all other material is lost
    WHERE(percolation_sl(:,nsoil_sb) > 0.0_wp)
      hlp1(:) = percolation_sl(:,nsoil_sb) * dtime/1000._wp * pool(:,nsoil_sb) * soil_depth_sl(:,nsoil_sb)
    ELSEWHERE
      hlp1(:) = 0.0_wp 
    ENDWHERE
    transport(:,nsoil_sb) = transport(:,nsoil_sb) - hlp1(:) / soil_depth_sl(:,nsoil_sb)
    leaching(:) = hlp1(:)

  END SUBROUTINE calc_liquid_phase_transport


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates rate of bioturbation 
  !!
  !! Input: bulk density 
  !!
  !! Output: vertical transport rate (1/s) 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_bioturbation_rate(nc, nsoil_sb, soil_depth_sl,bulk_dens_corr_sl,root_fraction_sl) RESULT(k_bioturb)

    USE mo_sb_constants,          ONLY: k_diff_org
    USE mo_jsb_math_constants,    ONLY: eps8

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                  !< dimensions
                                                                       nsoil_sb
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: soil_depth_sl          !< thickness of soil layers [m] 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: bulk_dens_corr_sl      !< OM corrected bulk density [kg m-3] 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: root_fraction_sl       !< root fraction in soil layer 
    REAL(wp), DIMENSION(nc, nsoil_sb)                               :: k_bioturb              !< diffusion factor for bioturbation 
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc, nsoil_sb)               :: hlp1
    INTEGER                                         :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bioturbation_rate'



    !> 1.0 calculate root density relative to top-soil root density  
    !!
    WHERE(soil_depth_sl(:,:) > eps8)
       hlp1(:,:) = root_fraction_sl(:,:) / soil_depth_sl(:,:)
    ELSEWHERE
       hlp1(:,:) = 0.0_wp
    ENDWHERE
    DO isoil = nsoil_sb,1,-1
       WHERE(hlp1(:,1) > eps8)
          hlp1(:,isoil) = hlp1(:,isoil)/hlp1(:,1)
       ELSEWHERE
          hlp1(:,isoil) = 0.0_wp
       ENDWHERE
    ENDDO

    !> 2.0 adjust bioturbation coefficient by bulk density and root density
    !!
    WHERE(bulk_dens_corr_sl(:,:) > eps8) 
       k_bioturb(:,:) = hlp1(:,:) * k_diff_org / bulk_dens_corr_sl(:,:)
    ELSEWHERE
       k_bioturb(:,:) = 0.0_wp
    ENDWHERE

  END FUNCTION calc_bioturbation_rate


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from calc_bioturbation_transport 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates matter transport by bioturbation, by calling calc_bioturbation_transport
  !! corresponds to the diffusion via bioturbation term in Ahrens et al. 2015, eq. 4, 6 and 7
  !!
  !! Input:  pools subject to bioturbation, percolation rate, and soil depth 
  !!
  !! Output: vertical transport rate (mol/m2/timestep) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_bioturbation_transport(nc, nsoil_sb,  k_bioturb, soil_depth_sl, pool,transport)

    USE mo_jsb_math_constants,    ONLY: eps8

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                  !< dimensions
                                                                       nsoil_sb
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: k_bioturb              !< diffusion factor for bioturbation,[m2/s] 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: soil_depth_sl          !< depth of each soil layer [m] 
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: pool                   !< current state of pool to be transported
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport              !< transport between layers per timestep 
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)         :: hlp1,hlp2 
    INTEGER                         :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bioturbation_transport'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 diffusion transport in the vertical axis of the soil (i.e. down)
    !! I can do 1,n-1 here, because by definition, there cannot be any upward transport from the uppermost layer
    !! I don't need nsoil_sb here, because the assumption is made that the gradient at the bottom of the
    !! soil column is zero
    DO isoil = 1,nsoil_sb-1

       ! distance between the centres of the soil layers
       hlp1(:) = (soil_depth_sl(:,isoil)/2.0_wp + soil_depth_sl(:,isoil+1)/2.0_wp)

       WHERE(hlp1(:) > eps8)  
          hlp2(:) = k_bioturb(:,isoil) * dtime * (pool(:,isoil) - pool(:,isoil+1)) / hlp1(:) 
       ELSEWHERE
          hlp2(:) = 0.0_wp
       ENDWHERE
       transport(:,isoil)   = transport(:,isoil)   - hlp2(:) / soil_depth_sl(:,isoil)
       transport(:,isoil+1) = transport(:,isoil+1) + hlp2(:) / soil_depth_sl(:,isoil+1)

   ENDDO

  END SUBROUTINE calc_bioturbation_transport

END MODULE mo_sb_jsm_transport
