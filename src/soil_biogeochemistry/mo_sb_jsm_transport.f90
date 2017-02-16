!>
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

  USE mo_kind,            ONLY: wp, eps4, eps8
  USE mo_exception,       ONLY: message, finish
  USE mo_base_memory,     ONLY: t_jsb_task_options
  USE mo_base_constants,  ONLY: dtime
  USE mo_sb_memory_class, ONLY: t_slm_struct_sb_pool,t_slm_struct_sb_element

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_particle_fluxrate, calc_bioturbation_rate, &
            calc_bioturbation_transport_wrapper, calc_liquid_phase_transport_wrapper , & 
            calc_particle_transport_wrapper

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
  SUBROUTINE calc_liquid_phase_transport_wrapper(options,nsoil,percolation_sl,soil_depth_l, &
                                         pool_dom, nh4_solute, no3_solute, po4_solute, nh4_n15_solute, no3_n15_solute, & 
                                         transport_dom, transport_nh4_solute, transport_no3_solute, & 
                                         transport_po4_solute, transport_nh4_n15_solute, transport_no3_n15_solute)
    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options
    INTEGER,                               INTENT(in)    :: nsoil
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: percolation_sl         !< fraction of pool transported up or down per second
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: soil_depth_l           !< depth of each soil layer [m] 
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_dom               !< DOM pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: nh4_solute             !< NH4 pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: no3_solute             !< NO3 pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: po4_solute             !< PO4 pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: nh4_n15_solute         !< N15H4 pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: no3_n15_solute         !< N15O3 pool [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: transport_dom          !< rate of DOM transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_nh4_solute   !< rate of NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_no3_solute   !< rate of NO3 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_po4_solute   !< rate of PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_nh4_n15_solute !< rate of N15H4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_no3_n15_solute !< rate of N15O3 transport [mol/m3/timestep]
    ! 0.2 Local
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_liquid_phase_transport_wrapper'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    !> 1.0 transport of DOM
    !!
    CALL calc_liquid_phase_transport(options, nsoil, percolation_sl(:,:), soil_depth_l(:,:), &
                                        pool_dom%carbon%ptr(ics:ice,:,iblk),transport_dom%carbon%ptr(ics:ice,:,iblk)) 
    CALL calc_liquid_phase_transport(options, nsoil, percolation_sl(:,:), soil_depth_l(:,:), &
                                        pool_dom%nitrogen%ptr(ics:ice,:,iblk),transport_dom%nitrogen%ptr(ics:ice,:,iblk)) 
    CALL calc_liquid_phase_transport(options, nsoil, percolation_sl(:,:), soil_depth_l(:,:), &
                                        pool_dom%phosphorus%ptr(ics:ice,:,iblk),transport_dom%phosphorus%ptr(ics:ice,:,iblk)) 
    CALL calc_liquid_phase_transport(options, nsoil, percolation_sl(:,:), soil_depth_l(:,:), &
                                        pool_dom%carbon13%ptr(ics:ice,:,iblk),transport_dom%carbon13%ptr(ics:ice,:,iblk)) 
    CALL calc_liquid_phase_transport(options, nsoil, percolation_sl(:,:), soil_depth_l(:,:), &
                                        pool_dom%carbon14%ptr(ics:ice,:,iblk),transport_dom%carbon14%ptr(ics:ice,:,iblk)) 
    CALL calc_liquid_phase_transport(options, nsoil, percolation_sl(:,:), soil_depth_l(:,:), &
                                        pool_dom%nitrogen15%ptr(ics:ice,:,iblk),transport_dom%nitrogen15%ptr(ics:ice,:,iblk)) 

    !> 2.0 transport of other solutes (NH4, NO3, PO4, and their isotopes)
    CALL calc_liquid_phase_transport(options, nsoil, percolation_sl(:,:), soil_depth_l(:,:), &
                                        nh4_solute(:,:),transport_nh4_solute(:,:)) 
    CALL calc_liquid_phase_transport(options, nsoil, percolation_sl(:,:), soil_depth_l(:,:), &
                                        no3_solute(:,:),transport_no3_solute(:,:)) 
    CALL calc_liquid_phase_transport(options, nsoil, percolation_sl(:,:), soil_depth_l(:,:), &
                                        po4_solute(:,:),transport_po4_solute(:,:)) 
    CALL calc_liquid_phase_transport(options, nsoil, percolation_sl(:,:), soil_depth_l(:,:), &
                                        nh4_n15_solute(:,:),transport_nh4_n15_solute(:,:)) 
    CALL calc_liquid_phase_transport(options, nsoil, percolation_sl(:,:), soil_depth_l(:,:), &
                                        no3_n15_solute(:,:),transport_no3_n15_solute(:,:)) 

  END SUBROUTINE calc_liquid_phase_transport_wrapper

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates matter transport by bioturbation, by calling calc_bioturbation_transport
  !! partly necessary because of a missing operator overload...
  !!
  !! Input:  pools subject to bioturbation, percolation rate, and soil depth 
  !!
  !! Output: vertical transport rate (mol/m2/timestep) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_bioturbation_transport_wrapper(options, nsoil, k_bioturb, soil_depth_l, &
                                          pool_polymeric_litter, pool_dom_assoc, & 
                                          pool_microbial, pool_residue, pool_residue_assoc, & 
                                          nh4_assoc, po4_assoc, nh4_n15_assoc, & 
                                          transport_polymeric_litter, transport_dom_assoc, & 
                                          transport_microbial, transport_residue, & 
                                          transport_residue_assoc, transport_nh4_assoc, & 
                                          transport_po4_assoc, transport_nh4_n15_assoc)

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options
    INTEGER,                               INTENT(in)    :: nsoil
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: k_bioturb              !< diffusion factor for bioturbation 
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: soil_depth_l           !< depth of each soil layer [m] 
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_polymeric_litter  !< polymeric litter pool [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_dom_assoc         !< minerally associated DOM pool [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_microbial         !< microbial pool [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_residue           !< microbial residue pool [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_residue_assoc     !< minerally associated microbial residue pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: nh4_assoc              !< NH4 pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: po4_assoc              !< PO4 pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: nh4_n15_assoc          !< N15H4 pool [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: transport_polymeric_litter !< rate of polymeric litter transport [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: transport_dom_assoc        !< rate of minerally associated DOM transport [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: transport_microbial        !< rate of microbial transport [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: transport_residue          !< rate of microbial residue transport [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: transport_residue_assoc    !< rate of minerally associated microbial residue transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_nh4_assoc        !< rate of NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_po4_assoc        !< rate of PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_nh4_n15_assoc    !< rate of N15H4 transport [mol/m3/timestep]
    ! 0.2 Local
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bioturbation_transport_wrapper'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    !> 1.0 transport of organic pools
    !! 
    !> 1.1 polymeric litter
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_polymeric_litter%carbon%ptr(ics:ice,:,iblk), &
                                        transport_polymeric_litter%carbon%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_polymeric_litter%nitrogen%ptr(ics:ice,:,iblk), &
                                        transport_polymeric_litter%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_polymeric_litter%phosphorus%ptr(ics:ice,:,iblk), &
                                        transport_polymeric_litter%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_polymeric_litter%carbon13%ptr(ics:ice,:,iblk), &
                                        transport_polymeric_litter%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_polymeric_litter%carbon14%ptr(ics:ice,:,iblk), &
                                        transport_polymeric_litter%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk), &
                                        transport_polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk))

    !> 1.2 minerally associated DOM 
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_dom_assoc%carbon%ptr(ics:ice,:,iblk), &
                                        transport_dom_assoc%carbon%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_dom_assoc%nitrogen%ptr(ics:ice,:,iblk), &
                                        transport_dom_assoc%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_dom_assoc%phosphorus%ptr(ics:ice,:,iblk), &
                                        transport_dom_assoc%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_dom_assoc%carbon13%ptr(ics:ice,:,iblk), &
                                        transport_dom_assoc%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_dom_assoc%carbon14%ptr(ics:ice,:,iblk), &
                                        transport_dom_assoc%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_dom_assoc%nitrogen15%ptr(ics:ice,:,iblk), &
                                        transport_dom_assoc%nitrogen15%ptr(ics:ice,:,iblk))

    !> 1.3 microbial biomass 
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_microbial%carbon%ptr(ics:ice,:,iblk), &
                                        transport_microbial%carbon%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_microbial%nitrogen%ptr(ics:ice,:,iblk), &
                                        transport_microbial%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_microbial%phosphorus%ptr(ics:ice,:,iblk), &
                                        transport_microbial%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_microbial%carbon13%ptr(ics:ice,:,iblk), &
                                        transport_microbial%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_microbial%carbon14%ptr(ics:ice,:,iblk), &
                                        transport_microbial%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_microbial%nitrogen15%ptr(ics:ice,:,iblk), &
                                        transport_microbial%nitrogen15%ptr(ics:ice,:,iblk))

    !> 1.4 microbial residue 
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_residue%carbon%ptr(ics:ice,:,iblk), &
                                        transport_residue%carbon%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_residue%nitrogen%ptr(ics:ice,:,iblk), &
                                        transport_residue%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_residue%phosphorus%ptr(ics:ice,:,iblk), &
                                        transport_residue%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_residue%carbon13%ptr(ics:ice,:,iblk), &
                                        transport_residue%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_residue%carbon14%ptr(ics:ice,:,iblk), &
                                        transport_residue%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_residue%nitrogen15%ptr(ics:ice,:,iblk), &
                                        transport_residue%nitrogen15%ptr(ics:ice,:,iblk))

    !> 1.5 minerally associated microbial residue 
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_residue_assoc%carbon%ptr(ics:ice,:,iblk), &
                                        transport_residue_assoc%carbon%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_residue_assoc%nitrogen%ptr(ics:ice,:,iblk), &
                                        transport_residue_assoc%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_residue_assoc%phosphorus%ptr(ics:ice,:,iblk), &
                                        transport_residue_assoc%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_residue_assoc%carbon13%ptr(ics:ice,:,iblk), &
                                        transport_residue_assoc%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_residue_assoc%carbon14%ptr(ics:ice,:,iblk), &
                                        transport_residue_assoc%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        pool_residue_assoc%nitrogen15%ptr(ics:ice,:,iblk), &
                                        transport_residue_assoc%nitrogen15%ptr(ics:ice,:,iblk))

    !> 2.0 inorganic pools 
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        nh4_assoc,transport_nh4_assoc) 
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        po4_assoc,transport_po4_assoc) 
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l, & 
                                        nh4_n15_assoc,transport_nh4_n15_assoc) 


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
  FUNCTION calc_particle_fluxrate(options,nsoil,soil_depth_l,bulk_dens_corr_sl,sb_formation,sb_loss) & 
     RESULT (particle_fluxrate)

    USE mo_sb_constants,          ONLY: carbon_per_dryweight_SOM 
    USE mo_shared_constants,      ONLY: molar_mass_C

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in) :: options
    INTEGER,                               INTENT(in) :: nsoil
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in) :: soil_depth_l           !< depth of soil layer [m] 
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in) :: bulk_dens_corr_sl      !< bulk density of mineral+organic soil [kg/m3]
    TYPE(t_slm_struct_sb_pool),            INTENT(in) :: sb_formation           !< formation of pools [mol/m3/timestep]
    TYPE(t_slm_struct_sb_pool),            INTENT(in) :: sb_loss                !< loss from pools [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil)             :: particle_fluxrate 
    ! 0.2 Local
    ! 0.2 Local
    INTEGER                     :: iblk, ics, ice, nc
    INTEGER                     :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_particle_fluxrate'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    !> 1.0 calculate net change of organic matter cummulative down the profile 
    !! note the omission of woody litter in the first layer, which is assumed not to be part of the profile 
    !! 
    particle_fluxrate(:,1) = soil_depth_l(:,1) * & 
          (sb_formation%soluable_litter%carbon%ptr(ics:ice,1,iblk)  - sb_loss%soluable_litter%carbon%ptr(ics:ice,1,iblk) + & 
           sb_formation%polymeric_litter%carbon%ptr(ics:ice,1,iblk) - sb_loss%polymeric_litter%carbon%ptr(ics:ice,1,iblk) + & 
           sb_formation%dom%carbon%ptr(ics:ice,1,iblk)              - sb_loss%dom%carbon%ptr(ics:ice,1,iblk) + & 
           sb_formation%dom_assoc%carbon%ptr(ics:ice,1,iblk)        - sb_loss%dom_assoc%carbon%ptr(ics:ice,1,iblk) + & 
           sb_formation%microbial%carbon%ptr(ics:ice,1,iblk)        - sb_loss%microbial%carbon%ptr(ics:ice,1,iblk) + & 
           sb_formation%residue%carbon%ptr(ics:ice,1,iblk)          - sb_loss%residue%carbon%ptr(ics:ice,1,iblk) + & 
           sb_formation%residue_assoc%carbon%ptr(ics:ice,1,iblk)    - sb_loss%residue_assoc%carbon%ptr(ics:ice,1,iblk))
          
    DO isoil = 2,nsoil
       particle_fluxrate(:,isoil) = particle_fluxrate(:,isoil-1) + soil_depth_l(:,isoil) * & 
          (sb_formation%soluable_litter%carbon%ptr(ics:ice,isoil,iblk)  - sb_loss%soluable_litter%carbon%ptr(ics:ice,isoil,iblk) + & 
           sb_formation%polymeric_litter%carbon%ptr(ics:ice,isoil,iblk) - sb_loss%polymeric_litter%carbon%ptr(ics:ice,isoil,iblk) + & 
           sb_formation%woody_litter%carbon%ptr(ics:ice,isoil,iblk)     - sb_loss%woody_litter%carbon%ptr(ics:ice,isoil,iblk) + & 
           sb_formation%dom%carbon%ptr(ics:ice,isoil,iblk)              - sb_loss%dom%carbon%ptr(ics:ice,isoil,iblk) + & 
           sb_formation%dom_assoc%carbon%ptr(ics:ice,isoil,iblk)        - sb_loss%dom_assoc%carbon%ptr(ics:ice,isoil,iblk) + & 
           sb_formation%microbial%carbon%ptr(ics:ice,isoil,iblk)        - sb_loss%microbial%carbon%ptr(ics:ice,isoil,iblk) + & 
           sb_formation%residue%carbon%ptr(ics:ice,isoil,iblk)          - sb_loss%residue%carbon%ptr(ics:ice,isoil,iblk) + & 
           sb_formation%residue_assoc%carbon%ptr(ics:ice,isoil,iblk)    - sb_loss%residue_assoc%carbon%ptr(ics:ice,isoil,iblk))
   ENDDO  
  
   !> 2.0 convert change in pools in density change and thus transport term 
   !! bulk-density is converted to mol C
   WHERE(bulk_dens_corr_sl(:,:) > eps8)
      particle_fluxrate(:,:) = particle_fluxrate(:,:) / & 
            (bulk_dens_corr_sl(:,:) * carbon_per_dryweight_SOM * molar_mass_C * 1000._wp) 
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
  SUBROUTINE calc_particle_transport_wrapper(options, nsoil, particle_fluxrate, soil_depth_l, &
                                          sb_pool, nh4_solute, nh4_assoc, no3_solute, po4_solute, po4_assoc, & 
                                          nh4_n15_solute, nh4_n15_assoc, no3_n15_solute, & 
                                          sb_transport, transport_nh4_solute, transport_nh4_assoc, transport_no3_solute, & 
                                          transport_po4_solute, transport_po4_assoc, transport_nh4_n15_solute, &
                                          transport_nh4_n15_assoc, transport_no3_n15_solute)

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options
    INTEGER,                               INTENT(in)    :: nsoil
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: particle_fluxrate          !< particle fluxrate (1/timestep) 
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: soil_depth_l               !< depth of each soil layer [m] 
    TYPE(t_slm_struct_sb_pool),            INTENT(in)    :: sb_pool                    !< organic pools [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: nh4_solute                 !< NH4 pool in soil solution [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: nh4_assoc                  !< minerally associated NH4 pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: no3_solute                 !< NO3 pool in soil solution [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: po4_solute                 !< PO4 pool in soil solution [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: po4_assoc                  !< minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: nh4_n15_solute             !< N15H4 pool in soil solution [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: nh4_n15_assoc              !< minerally associated N15H4 pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: no3_n15_solute             !< N15O3 pool in soil solution [mol/m3]
    TYPE(t_slm_struct_sb_pool),            INTENT(inout) :: sb_transport               !< xminerally associated DOM pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_nh4_solute       !< rate of solute NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_nh4_assoc        !< rate of associated NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_no3_solute       !< rate of solute NO3 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_po4_solute       !< rate of solute PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_po4_assoc        !< rate of associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_nh4_n15_solute   !< rate of solute N15H4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_no3_n15_solute   !< rate of solute N15O3 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport_nh4_n15_assoc    !< rate of associated N15H4 transport [mol/m3/timestep]
    ! 0.2 Local
    INTEGER                               :: iblk, ics, ice, nc
    INTEGER                               :: isoil
    REAL(wp), DIMENSION(options%nc,nsoil) :: particle_fluxrate_woody
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_particle_transport_wrapper'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    !> 1.0 transport of organic pools
    !! 
    !> 1.1 soluable litter
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%soluable_litter%carbon%ptr(ics:ice,:,iblk), &
                                        sb_transport%soluable_litter%carbon%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%soluable_litter%nitrogen%ptr(ics:ice,:,iblk), &
                                        sb_transport%soluable_litter%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%soluable_litter%phosphorus%ptr(ics:ice,:,iblk), &
                                        sb_transport%soluable_litter%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%soluable_litter%carbon13%ptr(ics:ice,:,iblk), &
                                        sb_transport%soluable_litter%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%soluable_litter%carbon14%ptr(ics:ice,:,iblk), &
                                        sb_transport%soluable_litter%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk), &
                                        sb_transport%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk))

    !> 1.2 polymeric litter
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk), &
                                        sb_transport%polymeric_litter%carbon%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk), &
                                        sb_transport%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk), &
                                        sb_transport%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk), &
                                        sb_transport%polymeric_litter%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk), &
                                        sb_transport%polymeric_litter%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk), &
                                        sb_transport%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk))

    !> 1.3 woody litter
    !! first layer of woody litter does not belong to profile 
    DO isoil = 1,nsoil
      particle_fluxrate_woody(:,isoil) = particle_fluxrate(:,isoil) - particle_fluxrate(:,1)
    ENDDO
    CALL calc_particle_transport(options,nsoil,particle_fluxrate_woody,soil_depth_l, & 
                                        sb_pool%woody_litter%carbon%ptr(ics:ice,:,iblk), &
                                        sb_transport%woody_litter%carbon%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate_woody,soil_depth_l, & 
                                        sb_pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk), &
                                        sb_transport%woody_litter%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate_woody,soil_depth_l, & 
                                        sb_pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk), &
                                        sb_transport%woody_litter%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate_woody,soil_depth_l, & 
                                        sb_pool%woody_litter%carbon13%ptr(ics:ice,:,iblk), &
                                        sb_transport%woody_litter%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate_woody,soil_depth_l, & 
                                        sb_pool%woody_litter%carbon14%ptr(ics:ice,:,iblk), &
                                        sb_transport%woody_litter%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate_woody,soil_depth_l, & 
                                        sb_pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk), &
                                        sb_transport%woody_litter%nitrogen15%ptr(ics:ice,:,iblk))

    !> 1.5 DOM 
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%dom%carbon%ptr(ics:ice,:,iblk), &
                                        sb_transport%dom%carbon%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%dom%nitrogen%ptr(ics:ice,:,iblk), &
                                        sb_transport%dom%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%dom%phosphorus%ptr(ics:ice,:,iblk), &
                                        sb_transport%dom%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%dom%carbon13%ptr(ics:ice,:,iblk), &
                                        sb_transport%dom%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%dom%carbon14%ptr(ics:ice,:,iblk), &
                                        sb_transport%dom%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%dom%nitrogen15%ptr(ics:ice,:,iblk), &
                                        sb_transport%dom%nitrogen15%ptr(ics:ice,:,iblk))

    !> 1.6 minerally associated DOM 
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%dom_assoc%carbon%ptr(ics:ice,:,iblk), &
                                        sb_transport%dom_assoc%carbon%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%dom_assoc%nitrogen%ptr(ics:ice,:,iblk), &
                                        sb_transport%dom_assoc%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%dom_assoc%phosphorus%ptr(ics:ice,:,iblk), &
                                        sb_transport%dom_assoc%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%dom_assoc%carbon13%ptr(ics:ice,:,iblk), &
                                        sb_transport%dom_assoc%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%dom_assoc%carbon14%ptr(ics:ice,:,iblk), &
                                        sb_transport%dom_assoc%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk), &
                                        sb_transport%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk))

    !> 1.7 microbial biomass 
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%microbial%carbon%ptr(ics:ice,:,iblk), &
                                        sb_transport%microbial%carbon%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%microbial%nitrogen%ptr(ics:ice,:,iblk), &
                                        sb_transport%microbial%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%microbial%phosphorus%ptr(ics:ice,:,iblk), &
                                        sb_transport%microbial%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%microbial%carbon13%ptr(ics:ice,:,iblk), &
                                        sb_transport%microbial%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%microbial%carbon14%ptr(ics:ice,:,iblk), &
                                        sb_transport%microbial%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%microbial%nitrogen15%ptr(ics:ice,:,iblk), &
                                        sb_transport%microbial%nitrogen15%ptr(ics:ice,:,iblk))

    !> 1.8 microbial residue 
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%residue%carbon%ptr(ics:ice,:,iblk), &
                                        sb_transport%residue%carbon%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%residue%nitrogen%ptr(ics:ice,:,iblk), &
                                        sb_transport%residue%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%residue%phosphorus%ptr(ics:ice,:,iblk), &
                                        sb_transport%residue%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%residue%carbon13%ptr(ics:ice,:,iblk), &
                                        sb_transport%residue%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%residue%carbon14%ptr(ics:ice,:,iblk), &
                                        sb_transport%residue%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%residue%nitrogen15%ptr(ics:ice,:,iblk), &
                                        sb_transport%residue%nitrogen15%ptr(ics:ice,:,iblk))

    !> 1.9 minerally associated microbial residue 
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%residue_assoc%carbon%ptr(ics:ice,:,iblk), &
                                        sb_transport%residue_assoc%carbon%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%residue_assoc%nitrogen%ptr(ics:ice,:,iblk), &
                                        sb_transport%residue_assoc%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%residue_assoc%phosphorus%ptr(ics:ice,:,iblk), &
                                        sb_transport%residue_assoc%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%residue_assoc%carbon13%ptr(ics:ice,:,iblk), &
                                        sb_transport%residue_assoc%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%residue_assoc%carbon14%ptr(ics:ice,:,iblk), &
                                        sb_transport%residue_assoc%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        sb_pool%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk), &
                                        sb_transport%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk))

    !> 2.0 inorganic pools 
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        nh4_solute,transport_nh4_solute) 
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        nh4_assoc,transport_nh4_assoc) 
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        no3_solute,transport_no3_solute) 
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        po4_solute,transport_po4_solute) 
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        po4_assoc,transport_po4_assoc) 
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        nh4_n15_solute,transport_nh4_n15_solute) 
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        nh4_n15_assoc,transport_nh4_n15_assoc) 
    CALL calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l, & 
                                        no3_n15_solute,transport_no3_n15_solute) 


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
  SUBROUTINE calc_particle_transport(options,nsoil,particle_fluxrate,soil_depth_l,pool,transport)

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options
    INTEGER,                               INTENT(in)    :: nsoil
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: particle_fluxrate      !< fraction of pool transported up or down per timestep 
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: soil_depth_l           !< depth of each soil layer [m] 
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: pool                   !< current state of pool to be transported
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport              !< transport between layers per timestep 
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc) :: hlp1 
    INTEGER                         :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_particle_transport'

    !> 1.0 advective transport in the vertical axis of the soil (i.e. down)
    !! I can do 1,n-1 here, because by definition, there cannot be any upward transport from the uppermost layer
    DO isoil = 1,nsoil-1
       WHERE(particle_fluxrate(:,isoil)>0.0_wp)
          hlp1 = particle_fluxrate(:,isoil) * pool(:,isoil) * soil_depth_l(:,isoil)
       ELSEWHERE
          hlp1 = particle_fluxrate(:,isoil) * pool(:,isoil+1) * soil_depth_l(:,isoil+1)
       ENDWHERE
       transport(:,isoil)   = transport(:,isoil)   - hlp1 / soil_depth_l(:,isoil)
       transport(:,isoil+1) = transport(:,isoil+1) + hlp1 / soil_depth_l(:,isoil+1)
       
   ENDDO

   !> 2.0 lower boundary condition,??
   !!  based on the assumption that there is no upward flux, and that all other 
   !!  material is lost
!   WHERE(particle_fluxrate(:,nsoil)>0.0_wp)
!      hlp1 = particle_fluxrate(:,nsoil) * pool(:,nsoil) * soil_depth_l(:,nsoil)
!   ELSEWHERE
      hlp1 = 0.0_wp 
!   ENDWHERE
   transport(:,nsoil) = transport(:,nsoil) - hlp1 / soil_depth_l(:,nsoil)

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
  SUBROUTINE calc_liquid_phase_transport(options,nsoil,percolation_sl,soil_depth_l,pool,transport)

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options
    INTEGER,                               INTENT(in)    :: nsoil
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: percolation_sl         !< fraction of pool transported up or down per second
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: soil_depth_l           !< depth of each soil layer [m] 
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: pool                   !< current state of pool to be transported
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport              !< transport between layers per timestep 
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc) :: hlp1 
    INTEGER                         :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_liquid_transport'

    !> 1.0 advective transport in the vertical axis of the soil (i.e. down)
    !! I can do 1,n-1 here, because by definition, there cannot be any upward transport from the uppermost layer
    DO isoil = 1,nsoil-1
       WHERE(percolation_sl(:,isoil)>0.0_wp)
          hlp1 = percolation_sl(:,isoil) * dtime/1000._wp * pool(:,isoil) * soil_depth_l(:,isoil)
!! BLARPP for now no upward transport - need to properly think about this in update_soil_physics
       ELSEWHERE
          hlp1 = 0.0_wp ! percolation_sl(:,isoil) * dtime/1000._wp * pool(:,isoil+1) * soil_depth_l(:,isoil+1)
       ENDWHERE
       transport(:,isoil)   = transport(:,isoil)   - hlp1 / soil_depth_l(:,isoil)
       transport(:,isoil+1) = transport(:,isoil+1) + hlp1 / soil_depth_l(:,isoil+1)
       
   ENDDO

   !> 2.0 lower boundary condition,
   !!  based on the assumption that there is no upward flux, and that all other 
   !!  material is lost
   WHERE(percolation_sl(:,nsoil)>0.0_wp)
      hlp1 = percolation_sl(:,nsoil) * dtime/1000._wp * pool(:,nsoil) * soil_depth_l(:,nsoil)
   ELSEWHERE
      hlp1 = 0.0_wp 
   ENDWHERE
   transport(:,nsoil) = transport(:,nsoil) - hlp1 / soil_depth_l(:,nsoil)

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
  FUNCTION calc_bioturbation_rate(options,nsoil,bulk_dens_corr_sl) RESULT(k_bioturb)

    USE mo_sb_constants, ONLY: k_diff_org

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options
    INTEGER,                               INTENT(in)    :: nsoil
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: bulk_dens_corr_sl      !< OM corrected bulk density [kg m-3] 
    REAL(wp), DIMENSION(options%nc,nsoil)                :: k_bioturb              !< diffusion factor for bioturbation 
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bioturbation_rate'

    !> 1.0 
    WHERE(bulk_dens_corr_sl(:,:) > eps8) 
       k_bioturb(:,:) = k_diff_org / bulk_dens_corr_sl(:,:)
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
  SUBROUTINE calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_l,pool,transport)

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options
    INTEGER,                               INTENT(in)    :: nsoil
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: k_bioturb              !< diffusion factor for bioturbation 
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: soil_depth_l           !< depth of each soil layer [m] 
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: pool                   !< current state of pool to be transported
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: transport              !< transport between layers per timestep 
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc) :: hlp1,hlp2 
    INTEGER                         :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bioturbation_transport'

    !> 1.0 diffusion transport in the vertical axis of the soil (i.e. down)
    !! I can do 1,n-1 here, because by definition, there cannot be any upward transport from the uppermost layer
    !! I don't need nsoil here, because the assumption is made that the gradient at the bottom of the
    !! soil column is zero
    DO isoil = 1,nsoil-1

       ! distance between the centres of the soil layers
       hlp1 = (soil_depth_l(:,isoil)/2.0_wp + soil_depth_l(:,isoil+1)/2.0_wp)

       WHERE(hlp1 > eps8)  
          hlp2 = k_bioturb(:,isoil) * dtime * (pool(:,isoil) - pool(:,isoil+1)) / hlp1 
       ELSEWHERE
          hlp2 = 0.0_wp
       ENDWHERE
       transport(:,isoil)   = transport(:,isoil)   - hlp2 / soil_depth_l(:,isoil)
       transport(:,isoil+1) = transport(:,isoil+1) + hlp2 / soil_depth_l(:,isoil+1)

   ENDDO

  END SUBROUTINE calc_bioturbation_transport

END MODULE mo_sb_jsm_transport
