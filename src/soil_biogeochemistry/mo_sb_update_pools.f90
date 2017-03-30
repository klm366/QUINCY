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

  USE mo_kind,                  ONLY: wp, eps8, eps4
  USE mo_exception,       ONLY: message, message_text, finish
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
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_base_constants,        ONLY: dtime
    USE mo_base_config_class,     ONLY: base_config
    USE mo_soil_config_class,     ONLY: soil_config
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element    
    USE mo_shared_constants,      ONLY: lambda_C14
    USE mo_sb_config_class,       ONLY: sb_config
    USE mo_sb_constants
    USE mo_sb_jsm_processes
    USE mo_sb_jsm_transport
    USE mo_sb_simple_model 
    USE mo_shared_isotope
 
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
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc
    INTEGER                               :: isoil
!!
    REAL(wp), DIMENSION(options%nc,soil_config%nsoil_water) :: fine_root_carbon,hlp1
!!
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
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_po4
!    dsl4jsb_Real3D_onChunk :: nitrification 
!    dsl4jsb_Real3D_onChunk :: denitrification 
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl
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
    dsl4jsb_Get_var3D_onChunk(SOIL_,     t_soil_sl)                   ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_,     w_soil_sl)                   ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_,     w_soil_fc_sl)                ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_,     soil_depth_sl)               ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_,     root_fraction_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_,     percolation_sl)              ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_,     bulk_dens_sl)                ! in
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_depolymerisation)      ! soil: in / land+plant: out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_mic_uptake)            ! soil: in / land+plant: out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_sorption)              ! soil: in / land+plant: out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_desorption)            ! soil: in / land+plant: out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_hsc)                   ! soil: in / land+plant: out
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
!    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification)
!    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_assoc)                   ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc)                   ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_assoc)               ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_solute)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_assoc)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_solute)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_solute)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_solute)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_assoc)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_n15_solute)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       k_bioturb)                   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_org)                    ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_soil_carbon_sl)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       particle_fluxrate)           ! out


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
    SELECT CASE(TRIM(sb_config%sb_model_scheme))
      CASE('one_pool','four_pool')

         CALL update_sb_simple_model_process(tile,options)

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
            sb%pool%residue_assoc%carbon14%ptr(ics:ice,:,iblk)          * lambda_C14 * dtime

      CASE('jsm')

!BLARPP
!percolation_sl(:,soil_config%nsoil_water) = 0.0_wp
rmm_mic_uptake = 1.0_wp
rmm_depolymerisation = 1.0_wp
rmm_hsc = 1.0_wp
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
         net_mineralisation_nh4(:,:) = 0.0_wp
         net_mineralisation_nh4_n15(:,:) = 0.0_wp
         net_mineralisation_no3(:,:) = 0.0_wp
         net_mineralisation_no3_n15(:,:) = 0.0_wp
         net_mineralisation_po4(:,:) = 0.0_wp
         plant_uptake_nh4_sl(:,:) = 0.0_wp
         plant_uptake_nh4_n15_sl(:,:) = 0.0_wp
         plant_uptake_no3_sl(:,:) = 0.0_wp
         plant_uptake_no3_n15_sl(:,:) = 0.0_wp
         plant_uptake_po4_sl(:,:) = 0.0_wp
         transport_nh4_solute(:,:) = 0.0_wp 
         transport_nh4_assoc(:,:) = 0.0_wp 
         transport_no3_solute(:,:) = 0.0_wp 
         transport_po4_solute(:,:) = 0.0_wp 
         transport_po4_assoc(:,:) = 0.0_wp 
         transport_nh4_n15_solute(:,:) = 0.0_wp 
         transport_nh4_n15_assoc(:,:) = 0.0_wp 
         transport_no3_n15_solute(:,:) = 0.0_wp 

         !fine_root_carbon = root_fraction_sl * 6._wp !veg%pool%fine_root%carbon%ptr
         fine_root_carbon = root_fraction_sl * veg%pool%fine_root%carbon%ptr

         !> X.0 Inorganic N assimilation
         !!
         !> X.1 NH4
         CALL calc_inorganic_uptake(nh4_solute(:,:), sb%pool%microbial%carbon%ptr(ics:ice,:,iblk), fine_root_carbon(:,:), &
                                    rtm_mic_uptake(:,:), rmm_mic_uptake(:,:), rtm_hsc(:,:), rmm_hsc(:,:), & 
                                    vmax_mic_uptake_nh4, km_mic_uptake_nh4, vmax_plant_uptake_nh4, km_plant_uptake_nh4, & 
                                    net_mineralisation_nh4(:,:), plant_uptake_nh4_sl(:,:), & 
                                    nh4_n15_solute(:,:), eta_mic_uptake_nh4, eta_plant_uptake_nh4, & 
                                    net_mineralisation_nh4_n15(:,:), plant_uptake_nh4_n15_sl(:,:))
         !> X.2 NO3
         CALL calc_inorganic_uptake(no3_solute(:,:), sb%pool%microbial%carbon%ptr(ics:ice,:,iblk), fine_root_carbon(:,:), &
                                    rtm_mic_uptake(:,:), rmm_mic_uptake(:,:), rtm_hsc(:,:), rmm_hsc(:,:), & 
                                    vmax_mic_uptake_no3, km_mic_uptake_no3, vmax_plant_uptake_no3, km_plant_uptake_no3, & 
                                    net_mineralisation_no3(:,:), plant_uptake_no3_sl(:,:), & 
                                    no3_n15_solute(:,:), eta_mic_uptake_no3, eta_plant_uptake_no3, & 
                                    net_mineralisation_no3_n15(:,:), plant_uptake_no3_n15_sl(:,:))
         !> X.3 PO4
         CALL calc_inorganic_uptake(po4_solute(:,:), sb%pool%microbial%carbon%ptr(ics:ice,:,iblk), fine_root_carbon(:,:), &
                                    rtm_mic_uptake(:,:), rmm_mic_uptake(:,:), rtm_hsc(:,:), rmm_hsc(:,:), & 
                                    vmax_mic_uptake_po4, km_mic_uptake_po4, vmax_plant_uptake_po4, km_plant_uptake_po4, & 
                                    net_mineralisation_po4(:,:), plant_uptake_po4_sl(:,:))

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
                                    net_mineralisation_nh4(:,:),net_mineralisation_nh4_n15(:,:), & 
                                    net_mineralisation_no3(:,:),net_mineralisation_no3_n15(:,:), & 
                                    net_mineralisation_po4(:,:))

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
         nh4_solute(:,:)     = nh4_solute(:,:)     + (net_mineralisation_nh4(:,:) - plant_uptake_nh4_sl(:,:))/ 1000000.0_wp &
                                                   + transport_nh4_solute(:,:)
         nh4_n15_solute(:,:) = nh4_n15_solute(:,:) + (net_mineralisation_nh4_n15(:,:) - plant_uptake_nh4_n15_sl(:,:)) / 1000000.0_wp &
                                                   + transport_nh4_n15_solute(:,:)
         no3_solute(:,:)     = no3_solute(:,:)     + (net_mineralisation_no3(:,:) - plant_uptake_no3_sl(:,:)) / 1000000.0_wp &
                                                   + transport_no3_solute(:,:)
         no3_n15_solute(:,:) = no3_n15_solute(:,:) + (net_mineralisation_no3_n15(:,:) - plant_uptake_no3_n15_sl(:,:)) / 1000000.0_wp &
                                                   + transport_no3_n15_solute(:,:)
         po4_solute(:,:)     = po4_solute(:,:)     + (net_mineralisation_po4(:,:) - plant_uptake_po4_sl(:,:)) / 1000000.0_wp &
                                                   + transport_po4_solute(:,:)

         IF(sb_config%flag_sb_prescribe_nh4)THEN
            hlp1(:,:) = nh4_solute_prescribe - nh4_solute(:,:)
            nh4_solute(:,:) = nh4_solute(:,:) + hlp1(:,:) 
            nh4_n15_solute(:,:) = nh4_solute(:,:) * calc_mixing_ratio_N15N14(0.0_wp) 
         ENDIF
         IF(sb_config%flag_sb_prescribe_no3)THEN
            hlp1(:,:) = no3_solute_prescribe - no3_solute(:,:)
            no3_solute(:,:) = no3_solute(:,:) + hlp1(:,:) 
            no3_n15_solute(:,:) = no3_solute(:,:) * calc_mixing_ratio_N15N14(0.0_wp) 
         ENDIF
         IF(sb_config%flag_sb_prescribe_po4)THEN
            po4_solute(:,:) = po4_solute_prescribe 
         ENDIF
         !nh4_assoc(:,:) = nh4_assoc(:,:) + transport_nh4_assoc(:,:)
         !po4_assoc(:,:) = po4_assoc(:,:) + transport_po4_assoc(:,:)
         !nh4_n15_assoc(:,:) = nh4_n15_assoc(:,:) + transport_nh4_n15_assoc(:,:)

         !> 7.1 update diagnostics and other state variables 
         !!
         bulk_soil_carbon_sl = calc_bulk_soil_carbon(options,soil_config%nsoil_water, sb%pool)
         bulk_dens_corr_sl = calc_bulk_density_correction(bulk_soil_carbon_sl,bulk_dens_sl)
         qmax_org = calc_qmax_org(bulk_dens_corr_sl)

      CASE('canopy','test_canopy')

         CALL message(TRIM(routine), 'invalid model-type: '//TRIM(base_config%slm_model_name))

    ENDSELECT

  END SUBROUTINE update_sb_pools

END MODULE mo_sb_update_pools
