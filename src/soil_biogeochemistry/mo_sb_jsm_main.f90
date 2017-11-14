!>
!! @brief Contains the main routines for the jsm soil model
!!
!! @author
!!  Sönke Zaehle 
!!
!! @par Revision History
!! First version by Jan Engel (20170508): copyied from mo_sb_update_pools.f90
!!
MODULE mo_sb_jsm_main

#include "dsl_jsbach.inc"

  USE mo_kind,            ONLY: wp, eps8, eps4
  USE mo_exception,       ONLY: message, message_text, finish

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_jsm


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_jsm_main'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil biogeochemical pools 
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sb_jsm(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_base_constants,        ONLY: dtime
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_shared_isotope,        ONLY: calc_fractionation
    USE mo_sb_constants
    USE mo_sb_jsm_processes,      ONLY: calc_inorganic_uptake, calc_po4_uptake_ECA, calc_litter_turnover, calc_depolymerisation, &
                                        calc_transfer_litterfall_pools, calc_microbial_growth, calc_flux_po4, &
                                        calc_microbial_turnover, calc_om_sorption, calc_adsorption_po4
    USE mo_sb_jsm_transport,      ONLY: calc_liquid_phase_transport_wrapper, calc_bioturbation_transport_wrapper, &
                                        calc_bioturbation_rate, calc_particle_fluxrate, calc_particle_transport_wrapper

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)

    ! Use of process memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4slm_Use_memory_a_sb
    USE mo_sb_memory_class,       ONLY: sb_pool, sb_formation, sb_loss, sb_transport
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_soil_memory
    dsl4slm_Use_memory_a_veg        !  USE mo_veg_memory_class,      ONLY: t_veg_element_a, t_veg_pool_a
    USE mo_veg_memory_class,      ONLY: veg_pool, veg_litterfall


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout), TARGET  :: tile         !< one tile
    TYPE(t_jsb_task_options), INTENT(in)             :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: nsoil_sb
    INTEGER                               :: iblk, ics, ice, nc
!!
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: fine_root_carbon
!!
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_jsm'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)                               ! TYPE(t_sb_memory),       POINTER :: sb__mem
    dsl4jsb_Def_memory(SOIL_)                             ! TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(VEG_)                              ! TYPE(t_veg_memory),      POINTER :: veg__mem
    ! Declare pointers to variables in memory
    ! SOIL_
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    dsl4jsb_Real3D_onChunk :: percolation_sl 
    ! SB_
    dsl4jsb_Real3D_onChunk :: rtm_depolymerisation 
    dsl4jsb_Real3D_onChunk :: rtm_mic_uptake
    dsl4jsb_Real3D_onChunk :: rtm_sorption
    dsl4jsb_Real3D_onChunk :: rtm_desorption
    dsl4jsb_Real3D_onChunk :: rtm_hsc
    dsl4jsb_Real3D_onChunk :: rmm_depolymerisation 
    dsl4jsb_Real3D_onChunk :: rmm_mic_uptake
    dsl4jsb_Real3D_onChunk :: rmm_sorption
    dsl4jsb_Real3D_onChunk :: rmm_desorption
    dsl4jsb_Real3D_onChunk :: rmm_hsc
    dsl4jsb_Real3D_onChunk :: het_respiration
    dsl4jsb_Real3D_onChunk :: het_respiration_c13
    dsl4jsb_Real3D_onChunk :: het_respiration_c14
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: weathering_po4
    dsl4jsb_Real3D_onChunk :: occlusion_po4
    dsl4jsb_Real3D_onChunk :: slow_sorption_po4
    dsl4jsb_Real3D_onChunk :: fast_adsorption_po4
    dsl4jsb_Real3D_onChunk :: biochem_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: gross_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: nh4_solute 
    dsl4jsb_Real3D_onChunk :: nh4_assoc 
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute 
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast 
    dsl4jsb_Real3D_onChunk :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: po4_occluded
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: nh4_n15_assoc 
    dsl4jsb_Real3D_onChunk :: no3_n15_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_assoc 
    dsl4jsb_Real3D_onChunk :: transport_no3_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_solute 
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_fast 
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_slow 
    dsl4jsb_Real3D_onChunk :: transport_po4_occluded
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_assoc 
    dsl4jsb_Real3D_onChunk :: transport_no3_n15_solute 
    dsl4jsb_Real3D_onChunk :: k_bioturb 
    dsl4jsb_Real3D_onChunk :: qmax_org 
    dsl4jsb_Real3D_onChunk :: qmax_po4
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl 
    dsl4jsb_Real3D_onChunk :: particle_fluxrate 
    ! Get local variables from options argument
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    nsoil_sb  = options%nsoil_sb
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id) 
    ! Get process config
    dsl4jsb_Get_config(SB_)
    ! Get process memories
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(VEG_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     soil_depth_sl)               ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_,     root_fraction_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_,     percolation_sl)              ! in
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_depolymerisation)        ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_mic_uptake)              ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_sorption)                ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_desorption)              ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_hsc)                     ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_depolymerisation)        ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_mic_uptake)              ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_sorption)                ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_desorption)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_hsc)                     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration)             ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c13)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c14)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4_n15)  ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3_n15)  ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_po4)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       weathering_po4)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       occlusion_po4)               ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       slow_sorption_po4)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       fast_adsorption_po4)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       biochem_mineralisation_po4)  ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       gross_mineralisation_po4)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_assoc)                   ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_fast)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_slow)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_occluded)                ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_assoc)               ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_solute)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_assoc)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_solute)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_solute)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_fast)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_slow)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_occluded)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_solute)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_assoc)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_n15_solute)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       k_bioturb)                   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_org)                    ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_po4)                    ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       particle_fluxrate)           ! out
    ! Get process memory pool structures
    dsl4slm_get_vegpool_var(pool,fine_root,carbon)
#define Get_veg_litterfall
#define Get_sb_pool
#define Get_sb_formation
#define Get_sb_loss
#define Get_sb_transport
#include "slm_pools_a_incl.f90"
#undef Get_veg_litterfall
#undef Get_sb_pool
#undef Get_sb_formation
#undef Get_sb_loss
#undef Get_sb_transport


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

!BLARPP
!percolation_sl(:,options%nsoil_sb) = 0.0_wp
rmm_mic_uptake = 1.0_wp
rmm_depolymerisation = 1.0_wp
rmm_hsc = 1.0_wp
rmm_sorption = 1.0_wp
rmm_desorption = 1.0_wp

    fine_root_carbon(:,:) = root_fraction_sl(:,:) * 600._wp 
    !fine_root_carbon(:,:) = root_fraction_sl(:,:) * veg_pool%fine_root%carbon(:)

    !> X.0 Inorganic N assimilation
    !!
    !> X.1 NH4
    CALL calc_inorganic_uptake(nh4_solute(:,:), sb_pool%microbial%carbon(:,:), fine_root_carbon(:,:), &
                               rtm_mic_uptake(:,:), rmm_mic_uptake(:,:), rtm_hsc(:,:), rmm_hsc(:,:), & 
                               vmax_mic_uptake_nh4, km_mic_uptake_nh4, vmax_plant_uptake_nh4, km_plant_uptake_nh4, & 
                               microbial_uptake_nh4_sl(:,:), plant_uptake_nh4_sl(:,:))
    microbial_uptake_nh4_n15_sl(:,:)  = microbial_uptake_nh4_sl(:,:) * &
                                          calc_fractionation(nh4_solute,nh4_n15_solute,eta_mic_uptake_nh4)
    plant_uptake_nh4_n15_sl(:,:)      = plant_uptake_nh4_sl(:,:) * &
                                          calc_fractionation(nh4_solute,nh4_n15_solute,eta_plant_uptake_nh4)

    !> X.2 NO3
    CALL calc_inorganic_uptake(no3_solute(:,:), sb_pool%microbial%carbon(:,:), fine_root_carbon(:,:), &
                               rtm_mic_uptake(:,:), rmm_mic_uptake(:,:), rtm_hsc(:,:), rmm_hsc(:,:), & 
                               vmax_mic_uptake_no3, km_mic_uptake_no3, vmax_plant_uptake_no3, km_plant_uptake_no3, & 
                               microbial_uptake_no3_sl(:,:), plant_uptake_no3_sl(:,:))
    microbial_uptake_no3_n15_sl(:,:)  = microbial_uptake_no3_sl(:,:) * &
                                          calc_fractionation(no3_solute,no3_n15_solute,eta_mic_uptake_no3)
    plant_uptake_no3_n15_sl(:,:)      = plant_uptake_no3_sl(:,:) * &
                                          calc_fractionation(no3_solute,no3_n15_solute,eta_plant_uptake_no3)

    !> X.3 PO4

    !! Calculate the P flux rates
    !!
    CALL calc_flux_po4(sb_pool%microbial%carbon(:,:), sb_pool%dom%carbon(:,:),sb_pool%dom%phosphorus(:,:),&
                            rtm_mic_uptake(:,:), rmm_mic_uptake(:,:), rtm_hsc(:,:), rmm_hsc(:,:), &
                            sb_pool%polymeric_litter%phosphorus(:,:), sb_pool%dom_assoc%phosphorus(:,:), &
                            sb_pool%residue%phosphorus(:,:), sb_pool%residue_assoc%phosphorus(:,:), &
                            sb_loss%polymeric_litter%phosphorus(:,:), sb_loss%dom_assoc%phosphorus(:,:), &
                            sb_loss%residue%phosphorus(:,:), sb_loss%residue_assoc%phosphorus(:,:), po4_assoc_slow(:,:), &
                            biochem_mineralisation_po4(:,:), gross_mineralisation_po4(:,:), &
                            weathering_po4(:,:), occlusion_po4(:,:), slow_sorption_po4(:,:))

    !! Microbial and plant P uptake
    !!
    CALL calc_po4_uptake_ECA(sb_pool%microbial%carbon(:,:), fine_root_carbon(:,:), po4_solute(:,:), &
                                rtm_mic_uptake(:,:), rmm_mic_uptake(:,:), rtm_hsc(:,:), rmm_hsc(:,:), &
                                qmax_po4(:,:), po4_assoc_fast(:,:), microbial_uptake_po4_sl(:,:), &
                                plant_uptake_po4_sl(:,:), dsl4jsb_Config(SB_)%flag_sb_ECA_upt_po4)

    !! P adsorption
    !!
    CALL calc_adsorption_po4(options, po4_solute(:,:), po4_assoc_fast(:,:), po4_assoc_slow(:,:), po4_occluded(:,:), &
                               qmax_po4(:,:), weathering_po4(:,:), occlusion_po4(:,:), slow_sorption_po4(:,:), &
                               percolation_sl(:,:), bulk_dens_corr_sl(:,:), &
                               gross_mineralisation_po4(:,:), biochem_mineralisation_po4(:,:), &
                               sb_pool%microbial%carbon(:,:), fine_root_carbon(:,:), &
                               microbial_uptake_po4_sl(:,:), plant_uptake_po4_sl(:,:), fast_adsorption_po4(:,:), &
                               TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &
                               dsl4jsb_Config(SB_)%flag_sb_CLE_po4, &
                               dsl4jsb_Config(SB_)%flag_sb_LK_po4)

    !> 2.0 Litter-related fluxes  
    !!
    !> 2.1 decay of soluable litter
    !!
    CALL calc_litter_turnover(options, tau_soluable_litter, &
                              rtm_depolymerisation(:,:),rmm_depolymerisation(:,:), &
                              sb_pool%soluable_litter,sb_loss%soluable_litter,sb_formation%dom)

    !> 2.2 decay of woody litter
    !!
    CALL calc_litter_turnover(options, tau_woody_litter, &
                              rtm_depolymerisation(:,:),rmm_depolymerisation(:,:), &
                              sb_pool%woody_litter,sb_loss%woody_litter,sb_formation%polymeric_litter)

    !> 2.3 co-metabolisation of polymeric litter 
    !! 
    CALL calc_depolymerisation(options,vmax_depolymerisation_litter, &
                               rtm_depolymerisation(:,:),rmm_depolymerisation(:,:), &
                               rtm_hsc(:,:),rmm_hsc(:,:), &
                               sb_pool%microbial%carbon(:,:), & 
                               sb_pool%polymeric_litter,sb_loss%polymeric_litter,sb_formation%dom)

    !> 2.4 co-metabolisation of polymeric residues
    !! 
    CALL calc_depolymerisation(options,vmax_depolymerisation_residue, &
                               rtm_depolymerisation(:,:),rmm_depolymerisation(:,:), &
                               rtm_hsc(:,:),rmm_hsc(:,:), &
                               sb_pool%microbial%carbon(:,:), & 
                               sb_pool%residue,sb_loss%residue,sb_formation%dom)

    !> 2.5 update litter pools with recent litter fall
    CALL calc_transfer_litterfall_pools(options,soil_depth_sl, & 
                                        root_fraction_sl,veg_litterfall,sb_formation%soluable_litter, & 
                                        sb_formation%polymeric_litter,sb_formation%woody_litter)
                                             

    !> 3.0 Microbial fluxes
    !!
    !> 3.1 microbial growth
    !! 
    CALL calc_microbial_growth(options,sb_pool%dom,sb_pool%microbial, & 
                               rtm_mic_uptake(:,:),rmm_mic_uptake(:,:),rtm_hsc(:,:),rmm_hsc(:,:), &
                               microbial_uptake_nh4_sl(:,:),microbial_uptake_nh4_n15_sl(:,:), & 
                               microbial_uptake_no3_sl(:,:),microbial_uptake_no3_n15_sl(:,:), microbial_uptake_po4_sl(:,:), &
                               sb_loss%dom,sb_formation%microbial, &
                               het_respiration(:,:),het_respiration_c13(:,:),het_respiration_c14(:,:), &
                               net_mineralisation_nh4(:,:),net_mineralisation_nh4_n15(:,:), & 
                               net_mineralisation_no3(:,:),net_mineralisation_no3_n15(:,:), & 
                               net_mineralisation_po4(:,:))

    !> 3.2 microbial death
    !! determine mortality rate, constrained by to avoiding a complete death of the community
    CALL calc_microbial_turnover(options,sb_pool%microbial,sb_loss%microbial,& 
                                 sb_formation%dom,sb_formation%residue)

    !> 4.0 mineral stabilisation
    !!
    !> 4.1 mineral stabilisation of organic matter
    !!
    CALL calc_om_sorption(options,qmax_org,rtm_sorption,rmm_sorption,rtm_desorption, &
                          rmm_desorption,sb_pool%dom,sb_pool%dom_assoc,sb_pool%residue, & 
                          sb_pool%residue_assoc,sb_formation%dom,sb_formation%dom_assoc,sb_formation%residue, &
                          sb_formation%residue_assoc,sb_loss%dom,sb_loss%dom_assoc,sb_loss%residue,sb_loss%residue_assoc)

    !> 5.0 vertical transport
    !!
    !> 5.1 liquid phase transport
    !!
    CALL calc_liquid_phase_transport_wrapper(options, & 
                                   percolation_sl(:,:), soil_depth_sl(:,:), &
                                   sb_pool%dom, nh4_solute(:,:), no3_solute(:,:), po4_solute(:,:), & 
                                   nh4_n15_solute(:,:), no3_n15_solute(:,:), & 
                                   sb_transport%dom, transport_nh4_solute(:,:), & 
                                   transport_no3_solute(:,:),transport_po4_solute(:,:), & 
                                   transport_nh4_n15_solute(:,:),transport_no3_n15_solute(:,:))

    !> 5.2 transport due to bioturbation
    !!
    k_bioturb(:,:) = calc_bioturbation_rate(options, bulk_dens_corr_sl(:,:))
    CALL calc_bioturbation_transport_wrapper(options, k_bioturb(:,:), soil_depth_sl(:,:), &
                                    sb_pool%polymeric_litter, sb_pool%dom_assoc, sb_pool%fungi, & 
                                    sb_pool%microbial, sb_pool%residue, sb_pool%residue_assoc, & 
                                    nh4_assoc(:,:), po4_assoc_fast(:,:), po4_assoc_slow(:,:), & 
                                    po4_occluded(:,:), nh4_n15_assoc(:,:), & 
                                    sb_transport%polymeric_litter, sb_transport%dom_assoc, & 
                                    sb_transport%fungi, sb_transport%microbial, sb_transport%residue, & 
                                    sb_transport%residue_assoc, transport_nh4_assoc(:,:), & 
                                    transport_po4_assoc_fast(:,:), transport_po4_assoc_slow(:,:), & 
                                    transport_po4_occluded(:,:), transport_nh4_n15_assoc(:,:))

    !> 5.3 particle transport
    !!
    particle_fluxrate(:,:) = calc_particle_fluxrate(options, & 
                                    soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), & 
                                    sb_formation, sb_loss) 
    CALL calc_particle_transport_wrapper(options, particle_fluxrate(:,:), soil_depth_sl(:,:), &
                                    sb_pool, &
                                    nh4_solute(:,:), nh4_assoc(:,:), no3_solute(:,:), po4_solute(:,:), & 
                                    po4_assoc_fast(:,:), po4_assoc_slow(:,:), po4_occluded(:,:), & 
                                    nh4_n15_solute(:,:), nh4_n15_assoc(:,:), no3_n15_solute(:,:), & 
                                    sb_transport, &
                                    transport_nh4_solute(:,:), transport_nh4_assoc(:,:), & 
                                    transport_no3_solute(:,:), transport_po4_solute(:,:), & 
                                    transport_po4_assoc_fast(:,:), transport_po4_assoc_slow(:,:), transport_po4_occluded(:,:), & 
                                    transport_nh4_n15_solute(:,:), transport_nh4_n15_assoc(:,:), transport_no3_n15_solute(:,:))

    !> 6.0 convert inorganic transport fluxes to µmol m3 s-1
    !!
    transport_nh4_solute(:,:)     = transport_nh4_solute(:,:) * 1.e6_wp / dtime
    transport_nh4_assoc(:,:)      = transport_nh4_assoc(:,:) * 1.e6_wp / dtime
    transport_no3_solute(:,:)     = transport_no3_solute(:,:) * 1.e6_wp / dtime
    transport_po4_solute(:,:)     = transport_po4_solute(:,:) * 1.e6_wp / dtime
    transport_po4_assoc_fast(:,:) = transport_po4_assoc_fast(:,:) * 1.e6_wp / dtime
    transport_po4_assoc_slow(:,:) = transport_po4_assoc_slow(:,:) * 1.e6_wp / dtime
    transport_po4_occluded(:,:)   = transport_po4_occluded(:,:) * 1.e6_wp / dtime
    transport_nh4_n15_solute(:,:) = transport_nh4_n15_solute(:,:) * 1.e6_wp / dtime
    transport_nh4_n15_assoc(:,:)  = transport_nh4_n15_assoc(:,:) * 1.e6_wp / dtime
    transport_no3_n15_solute(:,:) = transport_no3_n15_solute(:,:)  * 1.e6_wp / dtime


! Pass back values to process memory
dsl4slm_passback_vegpool_var(pool,fine_root,carbon)
#define Passback_veg_litterfall
#define Passback_sb_pool
#define Passback_sb_formation
#define Passback_sb_loss
#define Passback_sb_transport
#include "slm_pools_a_incl.f90"
#undef Passback_veg_litterfall
#undef Passback_sb_pool
#undef Passback_sb_formation
#undef Passback_sb_loss
#undef Passback_sb_transport

  END SUBROUTINE update_sb_jsm

END MODULE mo_sb_jsm_main
