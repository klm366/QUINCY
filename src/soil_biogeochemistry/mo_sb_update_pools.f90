!>
!! @brief Contains the routines for the Task update_sb_pools
!!
!! @author
!!  Silvia Caldararu, Sönke Zaehle 
!!
!! @par Revision History
!! First version by Sönke Zaehle (20170105): copyied from vegetation/mo_veg_soil.f90
!!
MODULE mo_sb_update_pools

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp
  USE mo_exception,       ONLY: message, finish
  USE mo_base_memory,     ONLY: t_jsb_task_options

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_pools


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_update_pools'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil biogeochemical pools 
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sb_pools(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_base_memory,           ONLY: t_slm_tile, t_model, t_jsb_task_options
    USE mo_base_constants,        ONLY: dtime
    USE mo_base_config_class,     ONLY: base_config
    USE mo_soil_config_class,     ONLY: soil_config
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element    
    USE mo_shared_constants,      ONLY: lambda_C14
    USE mo_sb_constants
    USE mo_sb_jsm_processes
    USE mo_sb_jsm_transport
 
    ! Use of process memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_soil_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout), TARGET  :: tile         !< one tile
    TYPE(t_jsb_task_options), INTENT(in)             :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_pools'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(SB_)                               ! TYPE(t_sb_memory),       POINTER :: sb__mem
    dsl4jsb_Def_memory(SOIL_)                             ! TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(VEG_)                             ! TYPE(t_soil_memory),     POINTER :: soil__mem
    ! Declare pointers to variables in memory
    dsl4jsb_Real3D_onChunk :: t_soil_sl
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    dsl4jsb_Real3D_onChunk :: w_soil_fc_sl
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    dsl4jsb_Real3D_onChunk :: percolation_sl 
    dsl4jsb_Real3D_onChunk :: bulk_dens_sl 
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
    dsl4jsb_Real3D_onChunk :: net_mineralisation_n
    dsl4jsb_Real3D_onChunk :: net_mineralisation_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_p
    dsl4jsb_Real3D_onChunk :: nh4_solute 
    dsl4jsb_Real3D_onChunk :: nh4_assoc 
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute 
    dsl4jsb_Real3D_onChunk :: po4_assoc 
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: nh4_n15_assoc 
    dsl4jsb_Real3D_onChunk :: no3_n15_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_assoc 
    dsl4jsb_Real3D_onChunk :: transport_no3_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_solute 
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc 
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_assoc 
    dsl4jsb_Real3D_onChunk :: transport_no3_n15_solute 
    dsl4jsb_Real3D_onChunk :: k_bioturb 
    dsl4jsb_Real3D_onChunk :: qmax_org 
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl 
    dsl4jsb_Real3D_onChunk :: bulk_soil_carbon_sl
    dsl4jsb_Real3D_onChunk :: particle_fluxrate 
! 2 be removed?: pointers to veg assimi *slmtemp
    TYPE(t_sb_memory), POINTER              :: sb         ! local pointer to sb__mem for SLM usability
    TYPE(t_veg_memory), POINTER             :: veg       ! local pointer to soil__mem for SLM usability
    ! Get local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(SB_)
    ! Get process memories
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(VEG_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id) 
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SOIL_,     t_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     w_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     w_soil_fc_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     root_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     percolation_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     bulk_dens_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_depolymerisation)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_mic_uptake)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_sorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_desorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_hsc)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_depolymerisation)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_mic_uptake)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_sorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_desorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_hsc)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c13)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c14)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_n)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_p)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_assoc) 
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_assoc) 
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_assoc) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_assoc) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       k_bioturb) 
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_org) 
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl) 
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_soil_carbon_sl) 
    dsl4jsb_Get_var3D_onChunk(SB_,       particle_fluxrate) 


! 2 be removed  *slmtemp
    sb   => sb__mem
    veg   => veg__mem


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 Calculate biophysical rate modifiers
    !!
    !! for the full and plant only model, update soil modifiers. In the soil model, these are read from a file 
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant')

         rtm_depolymerisation = calc_temperature_rate_modifier(t_soil_sl,ea_depolymerisation) 
         rtm_mic_uptake       = calc_temperature_rate_modifier(t_soil_sl,ea_mic_uptake)
         rtm_sorption         = calc_temperature_rate_modifier(t_soil_sl,ea_sorption) 
         rtm_desorption       = calc_temperature_rate_modifier(t_soil_sl,ea_desorption) 
         rtm_hsc              = calc_temperature_rate_modifier(t_soil_sl,ea_hsc) 
     
         rmm_depolymerisation = calc_oxygen_rate_modifier(w_soil_sl,w_soil_fc_sl,kmr_depolymerisation) 
         rmm_mic_uptake       = calc_oxygen_rate_modifier(w_soil_sl,w_soil_fc_sl,kmr_mic_uptake)
         rmm_sorption         = calc_convert_bulk_to_liquid(w_soil_sl)
         rmm_desorption       = calc_convert_bulk_to_liquid(w_soil_sl) 
         rmm_hsc              = calc_substrate_rate_modifier(w_soil_sl,w_soil_fc_sl,kmr_hsc) 

    ENDSELECT
 
    !> 2.0 Calculate soil pools 
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant')

         CALL update_sb_pools_simple(tile, options)

      CASE('soil')

!BLARPP
!percolation_sl(:,soil_config%nsoil_water) = 0.0_wp
!rmm_mic_uptake = 1.0_wp
!rmm_depolymerisation = 1.0_wp
!rmm_hsc = 1.0_wp

rmm_sorption = 1.0_wp
rmm_desorption = 1.0_wp
         !> 1.0 initialise fluxes
         !!
         CALL sb%loss%init(0.0_wp)
         CALL sb%formation%init(0.0_wp)
         CALL sb%transport%init(0.0_wp)
         het_respiration(:,:) = 0.0_wp
         het_respiration_c13(:,:) = 0.0_wp
         het_respiration_c14(:,:) = 0.0_wp
         net_mineralisation_n(:,:) = 0.0_wp
         net_mineralisation_n15(:,:) = 0.0_wp
         net_mineralisation_p(:,:) = 0.0_wp
         transport_nh4_solute(:,:) = 0.0_wp 
         transport_nh4_assoc(:,:) = 0.0_wp 
         transport_no3_solute(:,:) = 0.0_wp 
         transport_po4_solute(:,:) = 0.0_wp 
         transport_po4_assoc(:,:) = 0.0_wp 
         transport_nh4_n15_solute(:,:) = 0.0_wp 
         transport_nh4_n15_assoc(:,:) = 0.0_wp 
         transport_no3_n15_solute(:,:) = 0.0_wp 

         !> 2.0 Litter-related fluxes  
         !!
         !> 2.1 decay of soluable litter
         !!
         CALL calc_litter_turnover(options,soil_config%nsoil_water,tau_soluable_litter, &
                                   rtm_depolymerisation(:,:),rmm_depolymerisation(:,:), &
                                   sb%pool%soluable_litter,sb%loss%soluable_litter,sb%formation%dom)

         !> 2.1 decay of woody litter
         !!
         CALL calc_litter_turnover(options,soil_config%nsoil_water,tau_woody_litter, &
                                   rtm_depolymerisation(:,:),rmm_depolymerisation(:,:), &
                                   sb%pool%woody_litter,sb%loss%woody_litter,sb%formation%polymeric_litter)

         !> 2.3 co-metabolisation of polymeric litter 
         !! 
         CALL calc_depolymerisation(options,soil_config%nsoil_water,vmax_depolymerisation_litter, &
                                    rtm_depolymerisation(:,:),rmm_depolymerisation(:,:), &
                                    rtm_hsc(:,:),rmm_hsc(:,:), &
                                    sb%pool%microbial%carbon%ptr(ics:ice,:,iblk), & 
                                    sb%pool%polymeric_litter,sb%loss%polymeric_litter,sb%formation%dom)

         !> 2.4 update litter pools with recent litter fall
         CALL calc_transfer_litterfall_pools(options,soil_config%nsoil_water,soil_depth_sl, & 
                                             root_fraction_sl,veg%litterfall,sb%formation%soluable_litter, & 
                                             sb%formation%polymeric_litter,sb%formation%woody_litter)
                                             

         !> 3.0 Microbial fluxes
         !!
         !> 3.1 microbial growth
         !! currently unconstrained by nutrients 
         CALL calc_microbial_growth(options,soil_config%nsoil_water,sb%pool%dom,sb%pool%microbial, & 
                                    rtm_mic_uptake(:,:),rmm_mic_uptake(:,:),rtm_hsc(:,:),rmm_hsc(:,:), &
                                    sb%loss%dom,sb%formation%microbial, &
                                    het_respiration(:,:),het_respiration_c13(:,:),het_respiration_c14(:,:), &
                                    net_mineralisation_n(:,:),net_mineralisation_n15(:,:),net_mineralisation_p(:,:))

         !> 3.2 microbial death
         !! determine mortality rate, constrained by to avoiding a complete death of the community
         CALL calc_microbial_turnover(options,soil_config%nsoil_water,sb%pool%microbial,sb%loss%microbial,& 
                                      sb%formation%dom,sb%formation%residue)

         !> 3.3 co-metabolisation of polymeric residues
         !! 
         CALL calc_depolymerisation(options,soil_config%nsoil_water,vmax_depolymerisation_residue, &
                                    rtm_depolymerisation(:,:),rmm_depolymerisation(:,:), &
                                    rtm_hsc(:,:),rmm_hsc(:,:), &
                                    sb%pool%microbial%carbon%ptr(ics:ice,:,iblk), & 
                                    sb%pool%residue,sb%loss%residue,sb%formation%dom)

         !> 4.0 mineral stabilisation
         !!
         !> 4.1 mineral stabilisation of organic matter
         !!
         CALL calc_om_sorption(options,soil_config%nsoil_water,qmax_org,rtm_sorption,rmm_sorption,rtm_desorption, &
                              rmm_desorption,sb%pool%dom,sb%pool%dom_assoc,sb%pool%residue, & 
                              sb%pool%residue_assoc,sb%formation%dom,sb%formation%dom_assoc,sb%formation%residue, &
                              sb%formation%residue_assoc,sb%loss%dom,sb%loss%dom_assoc,sb%loss%residue,sb%loss%residue_assoc)

         !> 5.0 vertical transport
         !!
         !> 5.1 liquid phase transport
         !!
         CALL calc_liquid_phase_transport_wrapper(options,soil_config%nsoil_water, & 
                                          percolation_sl(:,:), soil_depth_sl(:,:), &
                                          sb%pool%dom, nh4_solute(:,:), no3_solute(:,:), po4_solute(:,:), & 
                                          nh4_n15_solute(:,:), no3_n15_solute(:,:), & 
                                          sb%transport%dom, transport_nh4_solute(:,:), & 
                                          transport_no3_solute(:,:),transport_po4_solute(:,:), & 
                                          transport_nh4_n15_solute(:,:),transport_no3_n15_solute(:,:))

         !> 5.2 transport due to bioturbation
         !!
         k_bioturb(:,:) = calc_bioturbation_rate(options, soil_config%nsoil_water, bulk_dens_corr_sl(:,:))
         CALL calc_bioturbation_transport_wrapper(options, soil_config%nsoil_water, k_bioturb(:,:), soil_depth_sl(:,:), &
                                          sb%pool%polymeric_litter, sb%pool%dom_assoc, & 
                                          sb%pool%microbial, sb%pool%residue, sb%pool%residue_assoc, & 
                                          nh4_assoc(:,:), po4_assoc(:,:), nh4_n15_assoc(:,:), & 
                                          sb%transport%polymeric_litter, sb%transport%dom_assoc, & 
                                          sb%transport%microbial, sb%transport%residue, & 
                                          sb%transport%residue_assoc, transport_nh4_assoc(:,:), & 
                                          transport_po4_assoc(:,:), transport_nh4_n15_assoc(:,:))

         !> 5.3 particle transport
         !!
         particle_fluxrate(:,:) = calc_particle_fluxrate(options,soil_config%nsoil_water, & 
                                          soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), & 
                                          sb%formation,sb%loss) 
         CALL calc_particle_transport_wrapper(options, soil_config%nsoil_water, particle_fluxrate(:,:), soil_depth_sl(:,:), &
                                          sb%pool, nh4_solute(:,:), nh4_assoc(:,:), no3_solute(:,:), po4_solute(:,:), & 
                                          po4_assoc(:,:), nh4_n15_solute(:,:), nh4_n15_assoc(:,:), no3_n15_solute(:,:), & 
                                          sb%transport, transport_nh4_solute(:,:), transport_nh4_assoc(:,:), & 
                                          transport_no3_solute(:,:), transport_po4_solute(:,:), transport_po4_assoc(:,:), & 
                                          transport_nh4_n15_solute(:,:), transport_nh4_n15_assoc(:,:), transport_no3_n15_solute(:,:))

         !> 6.0 radioactive decay
         !!
         sb%loss%soluable_litter%carbon14%ptr(ics:ice,:,iblk)  = sb%loss%soluable_litter%carbon14%ptr(ics:ice,:,iblk) + & 
            sb%pool%soluable_litter%carbon14%ptr(ics:ice,:,iblk)  * lambda_C14 * dtime
         sb%loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = sb%loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + &
            sb%pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) * lambda_C14 * dtime
         sb%loss%woody_litter%carbon14%ptr(ics:ice,:,iblk)     = sb%loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) + &
            sb%pool%woody_litter%carbon14%ptr(ics:ice,:,iblk)     * lambda_C14 * dtime
         sb%loss%dom%carbon14%ptr(ics:ice,:,iblk)              = sb%loss%dom%carbon14%ptr(ics:ice,:,iblk) + &
            sb%pool%dom%carbon14%ptr(ics:ice,:,iblk)              * lambda_C14 * dtime
         sb%loss%dom_assoc%carbon14%ptr(ics:ice,:,iblk)        = sb%loss%dom_assoc%carbon14%ptr(ics:ice,:,iblk) + &
            sb%pool%dom_assoc%carbon14%ptr(ics:ice,:,iblk)        * lambda_C14 * dtime
         sb%loss%microbial%carbon14%ptr(ics:ice,:,iblk)        = sb%loss%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
            sb%pool%microbial%carbon14%ptr(ics:ice,:,iblk)        * lambda_C14 * dtime
         sb%loss%residue%carbon14%ptr(ics:ice,:,iblk)          = sb%loss%residue%carbon14%ptr(ics:ice,:,iblk) + &
            sb%pool%residue%carbon14%ptr(ics:ice,:,iblk)          * lambda_C14 * dtime
         sb%loss%residue_assoc%carbon14%ptr(ics:ice,:,iblk)    = sb%loss%residue_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
            sb%pool%residue_assoc%carbon14%ptr(ics:ice,:,iblk)    * lambda_C14 * dtime 

         !> 7.0 update pools
         !!
         CALL calc_jsm_pool_change(options,sb%pool,sb%formation,sb%loss,sb%transport)
         !nh4_solute(:,:) = nh4_solute(:,:) + net_mineralisation_n(:,:) + transport_nh4_solute(:,:)
         nh4_solute(:,:) = nh4_solute(:,:) + transport_nh4_solute(:,:)
         nh4_assoc(:,:) = nh4_assoc(:,:) + transport_nh4_assoc(:,:)
         no3_solute(:,:) = no3_solute(:,:) + transport_no3_solute(:,:)
         !po4_solute(:,:) = po4_solute(:,:) + net_mineralisation_p(:,:) + transport_po4_solute(:,:)
         po4_solute(:,:) = po4_solute(:,:) + transport_po4_solute(:,:)
         po4_assoc(:,:) = po4_assoc(:,:) + transport_po4_assoc(:,:)
         !nh4_n15_solute(:,:) = nh4_n15_solute(:,:) + net_mineralisation_n15(:,:) + transport_nh4_n15_solute(:,:)
         nh4_n15_solute(:,:) = nh4_n15_solute(:,:) + transport_nh4_n15_solute(:,:)
         nh4_n15_assoc(:,:) = nh4_n15_assoc(:,:) + transport_nh4_n15_assoc(:,:)
         no3_n15_solute(:,:) = no3_n15_solute(:,:) + transport_no3_n15_solute(:,:)

         !> 7.1 update diagnostics and other state variables 
         !!
         bulk_soil_carbon_sl = calc_bulk_soil_carbon(options,soil_config%nsoil_water, sb%pool)
         bulk_dens_corr_sl = calc_bulk_density_correction(bulk_soil_carbon_sl,bulk_dens_sl)
         qmax_org = calc_qmax_org(bulk_dens_corr_sl)

      CASE('canopy','test_canopy')

         CALL message(TRIM(routine), 'invalid model-type: '//TRIM(base_config%slm_model_name))

    ENDSELECT

  END SUBROUTINE update_sb_pools


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task
  !
  !-----------------------------------------------------------------------------------------------------
  !> routine to update soil litter and mineral pools (simple version for vegetation model testing)
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sb_pools_simple(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model, t_jsb_task_options
    USE mo_base_constants,        ONLY: dtime,one_day,one_year
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_veg_config_class,      ONLY: veg_config
    USE mo_sb_constants

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout), TARGET    :: tile          !< one tile
    TYPE(t_jsb_task_options), INTENT(in)               :: options       !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib               !< land-cover-type library - parameter across pft's
    INTEGER                               :: isoil
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_pools_simple'
    ! ---------------------------
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(SOIL_)                         ! TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(SB_)                           ! TYPE(t_sb_memory),       POINTER :: sb__mem
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk      :: uptake_nh4
    dsl4jsb_Real2D_onChunk      :: uptake_no3
    dsl4jsb_Real2D_onChunk      :: uptake_p
!! *slmtemp    
    TYPE(t_veg_memory), POINTER             :: veg        ! local pointer to veg__mem for SLM usability
    TYPE(t_soil_memory), POINTER            :: soil       ! local pointer to soil__mem for SLM usability
    TYPE(t_sb_memory), POINTER              :: sb         ! local pointer to sb__mem for SLM usability
    ! Get local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(xx_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(SB_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_p)


!! *slmtemp
    veg  => veg__mem
    soil => soil__mem
    sb   => sb__mem



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 Calculate soil litter pool decomposition.
    !!
    sb%het_respiration%ptr(ics:ice,:,iblk) = sb%rtm_depolymerisation%ptr(ics:ice,:,iblk)/tau_litter/one_day/one_year * & 
                                             sb%pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) 
    sb%net_mineralisation_n%ptr(ics:ice,:,iblk) = sb%rtm_depolymerisation%ptr(ics:ice,:,iblk)/tau_litter/one_day/one_year * & 
                                                  sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) 
    sb%net_mineralisation_p%ptr(ics:ice,:,iblk) = sb%rtm_depolymerisation%ptr(ics:ice,:,iblk)/tau_litter/one_day/one_year * & 
                                                  sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) 
    sb%het_respiration_c13%ptr(ics:ice,:,iblk) = sb%rtm_depolymerisation%ptr(ics:ice,:,iblk)/tau_litter/one_day/one_year * & 
                                                 sb%pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) 
    sb%het_respiration_c14%ptr(ics:ice,:,iblk) = sb%rtm_depolymerisation%ptr(ics:ice,:,iblk)/tau_litter/one_day/one_year * & 
                                                 sb%pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) 
    sb%net_mineralisation_n15%ptr(ics:ice,:,iblk) = sb%rtm_depolymerisation%ptr(ics:ice,:,iblk)/tau_litter/one_day/one_year * & 
                                                    sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) 

    ! > 1.1 subtract losses
    !!
    sb%pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) = sb%pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) - & 
                                                          sb%het_respiration%ptr(ics:ice,:,iblk) * dtime
    sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) = sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                            sb%net_mineralisation_n%ptr(ics:ice,:,iblk) * dtime
    sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                              sb%net_mineralisation_p%ptr(ics:ice,:,iblk) * dtime
    sb%pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) = sb%pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) - & 
                                                            sb%het_respiration_c13%ptr(ics:ice,:,iblk) * dtime
    sb%pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = sb%pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) - & 
                                                            sb%het_respiration_c14%ptr(ics:ice,:,iblk) * dtime
    sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                              sb%net_mineralisation_n15%ptr(ics:ice,:,iblk) * dtime

    !> 2 Update soil litter pool with litter fall
    !!
    sb%pool%polymeric_litter%carbon%ptr(ics:ice,1,iblk) = sb%pool%polymeric_litter%carbon%ptr(ics:ice,1,iblk) + & 
       (veg%litterfall%leaf%carbon%ptr + veg%litterfall%sap_wood%carbon%ptr + & 
        veg%litterfall%heart_wood%carbon%ptr + veg%litterfall%fruit%carbon%ptr + & 
        veg%litterfall%labile%carbon%ptr + veg%litterfall%reserve%carbon%ptr + & 
        veg%litterfall%seed_bed%carbon%ptr ) / soil%soil_depth_sl%ptr(ics:ice,1,iblk)
    sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) = sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) +  & 
       (veg%litterfall%leaf%nitrogen%ptr + veg%litterfall%sap_wood%nitrogen%ptr + & 
        veg%litterfall%heart_wood%nitrogen%ptr + veg%litterfall%fruit%nitrogen%ptr + & 
        veg%litterfall%labile%nitrogen%ptr + veg%litterfall%reserve%nitrogen%ptr + & 
        veg%litterfall%seed_bed%nitrogen%ptr ) / soil%soil_depth_sl%ptr(ics:ice,1,iblk)
    sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) = sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) + & 
       (veg%litterfall%leaf%phosphorus%ptr + veg%litterfall%sap_wood%phosphorus%ptr + & 
        veg%litterfall%heart_wood%phosphorus%ptr + veg%litterfall%fruit%phosphorus%ptr + & 
        veg%litterfall%labile%phosphorus%ptr + veg%litterfall%reserve%phosphorus%ptr + & 
        veg%litterfall%seed_bed%phosphorus%ptr ) / soil%soil_depth_sl%ptr(ics:ice,1,iblk)
    sb%pool%polymeric_litter%carbon13%ptr(ics:ice,1,iblk) = sb%pool%polymeric_litter%carbon13%ptr(ics:ice,1,iblk) + & 
       (veg%litterfall%leaf%carbon13%ptr + veg%litterfall%sap_wood%carbon13%ptr + & 
        veg%litterfall%heart_wood%carbon13%ptr + veg%litterfall%fruit%carbon13%ptr + & 
        veg%litterfall%labile%carbon13%ptr + veg%litterfall%reserve%carbon13%ptr + & 
        veg%litterfall%seed_bed%carbon13%ptr ) / soil%soil_depth_sl%ptr(ics:ice,1,iblk)
    sb%pool%polymeric_litter%carbon14%ptr(ics:ice,1,iblk) = sb%pool%polymeric_litter%carbon14%ptr(ics:ice,1,iblk) + & 
       (veg%litterfall%leaf%carbon14%ptr + veg%litterfall%sap_wood%carbon14%ptr + & 
        veg%litterfall%heart_wood%carbon14%ptr + veg%litterfall%fruit%carbon14%ptr + & 
        veg%litterfall%labile%carbon14%ptr + veg%litterfall%reserve%carbon14%ptr + & 
        veg%litterfall%seed_bed%carbon14%ptr ) / soil%soil_depth_sl%ptr(ics:ice,1,iblk)
    sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) = sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) + &
       (veg%litterfall%leaf%nitrogen15%ptr + veg%litterfall%sap_wood%nitrogen15%ptr + & 
        veg%litterfall%heart_wood%nitrogen15%ptr + veg%litterfall%fruit%nitrogen15%ptr + & 
        veg%litterfall%labile%nitrogen15%ptr + veg%litterfall%reserve%nitrogen15%ptr + & 
        veg%litterfall%seed_bed%nitrogen15%ptr ) / soil%soil_depth_sl%ptr(ics:ice,1,iblk)

    DO isoil = 1,1

       sb%pool%polymeric_litter%carbon%ptr(ics:ice,isoil,iblk) = sb%pool%polymeric_litter%carbon%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%carbon%ptr + veg%litterfall%coarse_root%carbon%ptr ) / & ! * & 
          !soil%root_fraction_sl%ptr(ics:ice,isoil,iblk) / & 
            soil%soil_depth_sl%ptr(ics:ice,isoil,iblk) 
       sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,isoil,iblk) = sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%nitrogen%ptr + veg%litterfall%coarse_root%nitrogen%ptr ) / & ! * & 
          !soil%root_fraction_sl%ptr(ics:ice,isoil,iblk) / & 
            soil%soil_depth_sl%ptr(ics:ice,isoil,iblk) 
       sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,isoil,iblk) = sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%phosphorus%ptr + veg%litterfall%coarse_root%phosphorus%ptr ) / & ! * & 
          !soil%root_fraction_sl%ptr(ics:ice,isoil,iblk) / & 
            soil%soil_depth_sl%ptr(ics:ice,isoil,iblk) 
       sb%pool%polymeric_litter%carbon13%ptr(ics:ice,isoil,iblk) = sb%pool%polymeric_litter%carbon13%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%carbon13%ptr + veg%litterfall%coarse_root%carbon13%ptr ) / & ! * & 
          !soil%root_fraction_sl%ptr(ics:ice,isoil,iblk) / & 
            soil%soil_depth_sl%ptr(ics:ice,isoil,iblk) 
       sb%pool%polymeric_litter%carbon14%ptr(ics:ice,isoil,iblk) = sb%pool%polymeric_litter%carbon14%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%carbon14%ptr + veg%litterfall%coarse_root%carbon14%ptr ) / & ! * & 
          !soil%root_fraction_sl%ptr(ics:ice,isoil,iblk) / & 
            soil%soil_depth_sl%ptr(ics:ice,isoil,iblk) 
       sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,isoil,iblk) = sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%nitrogen15%ptr + veg%litterfall%coarse_root%nitrogen15%ptr ) / & ! * & 
          !soil%root_fraction_sl%ptr(ics:ice,isoil,iblk) / & 
            soil%soil_depth_sl%ptr(ics:ice,isoil,iblk) 

    ENDDO

    !> 3 Update mineral pool
    !!
    IF(veg_config%n_uptake_scheme /= 'prescribed') THEN
       sb%volatilisation_n%ptr(ics:ice,:,iblk) = sb%net_mineralisation_n%ptr(ics:ice,:,iblk) * floss_nmin * 1.e6_wp 
       sb%nitrification%ptr(ics:ice,:,iblk) = sb%net_mineralisation_n%ptr(ics:ice,:,iblk) * (1._wp - floss_nmin) * fnit_nmin * 1.e6_wp 
       sb%nh4_solute%ptr(ics:ice,:,iblk) = sb%nh4_solute%ptr(ics:ice,:,iblk) + & 
                                           sb%net_mineralisation_n%ptr(ics:ice,:,iblk) * (1._wp - floss_nmin) * (1._wp - fnit_nmin) * dtime
       sb%no3_solute%ptr(ics:ice,:,iblk) = sb%no3_solute%ptr(ics:ice,:,iblk) + & 
                                           sb%net_mineralisation_n%ptr(ics:ice,:,iblk) * (1._wp - floss_nmin) * fnit_nmin * dtime
       sb%nh4_solute%ptr(ics:ice,1,iblk) = sb%nh4_solute%ptr(ics:ice,1,iblk) - & 
                                           uptake_nh4(1:nc) / soil%soil_depth_sl%ptr(ics:ice,1,iblk) * dtime/1.e6_wp
       sb%no3_solute%ptr(ics:ice,1,iblk) = sb%no3_solute%ptr(ics:ice,1,iblk) - & 
                                           uptake_no3(1:nc) / soil%soil_depth_sl%ptr(ics:ice,1,iblk) * dtime/1.e6_wp
    ENDIF

    IF(veg_config%p_uptake_scheme /= 'prescribed') THEN
       sb%po4_solute%ptr(ics:ice,:,iblk) = sb%po4_solute%ptr(ics:ice,:,iblk) + sb%net_mineralisation_p%ptr(ics:ice,:,iblk) * dtime
       sb%po4_solute%ptr(ics:ice,1,iblk) = sb%po4_solute%ptr(ics:ice,1,iblk) - uptake_p(1:nc) / soil%soil_depth_sl%ptr(ics:ice,1,iblk) * dtime/1.e6_wp
    ENDIF

    END SUBROUTINE update_sb_pools_simple

END MODULE mo_sb_update_pools
