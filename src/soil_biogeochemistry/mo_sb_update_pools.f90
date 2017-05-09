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

  USE mo_kind,            ONLY: wp, eps8, eps4
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
    USE mo_sb_jsm_main,           ONLY: update_sb_jsm
 
    ! Use of process memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_soil_memory
    dsl4jsb_Use_memory(A2L_)

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
    dsl4jsb_Def_memory(A2L_)
    ! Declare pointers to variables in memory
    ! A2L
    dsl4jsb_Real2D_onChunk :: nhx_deposition
    dsl4jsb_Real2D_onChunk :: noy_deposition
    ! ..
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
    dsl4jsb_Get_memory(A2L_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id) 
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(A2L_,     nhx_deposition)               ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,     noy_deposition)               ! in
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
      CASE('simple')

        CALL update_sb_simple_model_process(tile,options)

        CALL calc_organic_pool_change(options,sb%pool,sb%formation,sb%loss,sb%transport)

        IF(sb_config%flag_n_deposition) THEN
          nh4_solute(:,1)     = nh4_solute(:,1)   + nhx_deposition(:) / soil_depth_sl(:,1) * dtime
          no3_solute(:,1)     = no3_solute(:,1)   + noy_deposition(:) / soil_depth_sl(:,1) * dtime
        END IF

        veg%uptake_nh4%ptr(ice:ice,iblk) = plant_uptake_nh4_sl(:,1) * soil_depth_sl(:,1)
        veg%uptake_no3%ptr(ice:ice,iblk) = plant_uptake_no3_sl(:,1) * soil_depth_sl(:,1)
        veg%uptake_p%ptr(ice:ice,iblk)   = plant_uptake_po4_sl(:,1) * soil_depth_sl(:,1)
        IF(.NOT.sb_config%flag_sb_prescribe_nh4) THEN
          nh4_solute(:,:)     = nh4_solute(:,:)     + (net_mineralisation_nh4(:,:) - sb%nitrification%ptr(ice:ice,:,iblk) & 
                                          - sb%volatilisation_n%ptr(ice:ice,:,iblk)- plant_uptake_nh4_sl(:,:)) * dtime/1.e6_wp
          nh4_n15_solute(:,:) = nh4_n15_solute(:,:) + net_mineralisation_nh4_n15(:,:) / 1000000.0_wp * dtime
        ENDIF
        IF(.NOT.sb_config%flag_sb_prescribe_no3) THEN
          no3_solute(:,:)     = no3_solute(:,:)     + (sb%nitrification%ptr(ice:ice,:,iblk) - plant_uptake_no3_sl(:,:)) / dtime/1.e6_wp 
          no3_n15_solute(:,:) = no3_n15_solute(:,:) + net_mineralisation_no3_n15(:,:) / 1000000.0_wp * dtime
        ENDIF
        IF(.NOT.sb_config%flag_sb_prescribe_po4) THEN
          po4_solute(:,:)     = po4_solute(:,:) + (net_mineralisation_po4(:,:) - plant_uptake_po4_sl(:,:)) * dtime/1.e6_wp 
        ENDIF

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

         CALL update_sb_jsm(tile, options)


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
         sb%loss%fungi%carbon14%ptr(ics:ice,:,iblk)            = sb%loss%fungi%carbon14%ptr(ics:ice,:,iblk) + & 
            sb%pool%fungi%carbon14%ptr(ics:ice,:,iblk)            * lambda_C14 * dtime
         sb%loss%mycorrhiza%carbon14%ptr(ics:ice,:,iblk)       = sb%loss%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) + & 
            sb%pool%mycorrhiza%carbon14%ptr(ics:ice,:,iblk)       * lambda_C14 * dtime
         sb%loss%microbial%carbon14%ptr(ics:ice,:,iblk)        = sb%loss%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
            sb%pool%microbial%carbon14%ptr(ics:ice,:,iblk)        * lambda_C14 * dtime
         sb%loss%residue%carbon14%ptr(ics:ice,:,iblk)          = sb%loss%residue%carbon14%ptr(ics:ice,:,iblk) + &
            sb%pool%residue%carbon14%ptr(ics:ice,:,iblk)          * lambda_C14 * dtime
         sb%loss%residue_assoc%carbon14%ptr(ics:ice,:,iblk)    = sb%loss%residue_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
            sb%pool%residue_assoc%carbon14%ptr(ics:ice,:,iblk)    * lambda_C14 * dtime 

         !> 7.0 update pools
         !!
         CALL calc_organic_pool_change(options,sb%pool,sb%formation,sb%loss,sb%transport)
         IF(sb_config%flag_n_deposition) THEN
           nh4_solute(:,1)     = nh4_solute(:,1)   + nhx_deposition(:) / soil_depth_sl(:,1) * dtime
           no3_solute(:,1)     = no3_solute(:,1)   + noy_deposition(:) / soil_depth_sl(:,1) * dtime
         END IF

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


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> updates all SB pools given the JSM equations for formation, loss and transport 
  !! should be a simple pool operator overload function, but that doesn't work yet
  !! 
  !! Input: formation, loss and transport rate, pool sizes 
  !!
  !! Output: updated pool sizes 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_organic_pool_change(options,pool,formation,loss,transport)

    USE mo_sb_memory_class,       ONLY: t_slm_struct_sb_pool, t_slm_struct_sb_element

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options), INTENT(in)          :: options      !< model options
    TYPE(t_slm_struct_sb_pool),  INTENT(in)       :: formation    !< rate of formation of a pool [mol/m3/timestep]
    TYPE(t_slm_struct_sb_pool),  INTENT(in)       :: loss         !< rate of loss of a pool [mol/m3/timestep]
    TYPE(t_slm_struct_sb_pool),  INTENT(in)       :: transport    !< rate of transport of a pool [mol/m3/timestep]
    TYPE(t_slm_struct_sb_pool),  INTENT(inout)    :: pool         !< the pool [mol/m3]
    ! ---------------------------
    ! 0.2 Local 
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_organic_pool_change'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Get local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! carbon
    pool%dom%carbon%ptr(ics:ice,:,iblk) = pool%dom%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%carbon%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%carbon%ptr(ics:ice,:,iblk) = pool%dom_assoc%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%carbon%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%carbon%ptr(ics:ice,:,iblk) = pool%soluable_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%carbon%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) = pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%carbon%ptr(ics:ice,:,iblk)  
    pool%woody_litter%carbon%ptr(ics:ice,:,iblk) = pool%woody_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%carbon%ptr(ics:ice,:,iblk)  
    pool%fungi%carbon%ptr(ics:ice,:,iblk) = pool%fungi%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%fungi%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%fungi%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%fungi%carbon%ptr(ics:ice,:,iblk)  
    pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk) = pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%mycorrhiza%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%mycorrhiza%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%mycorrhiza%carbon%ptr(ics:ice,:,iblk)  
    pool%microbial%carbon%ptr(ics:ice,:,iblk) = pool%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%carbon%ptr(ics:ice,:,iblk)  
    pool%residue%carbon%ptr(ics:ice,:,iblk) = pool%residue%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%carbon%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%carbon%ptr(ics:ice,:,iblk) = pool%residue_assoc%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%carbon%ptr(ics:ice,:,iblk)  

    ! nitrogen
    pool%dom%nitrogen%ptr(ics:ice,:,iblk) = pool%dom%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%nitrogen%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%nitrogen%ptr(ics:ice,:,iblk) = pool%dom_assoc%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%nitrogen%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) = pool%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%nitrogen%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) = pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)  
    pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk) = pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%nitrogen%ptr(ics:ice,:,iblk)  
    pool%fungi%nitrogen%ptr(ics:ice,:,iblk) = pool%fungi%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%fungi%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%fungi%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%fungi%nitrogen%ptr(ics:ice,:,iblk)  
    pool%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) = pool%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk)  
    pool%microbial%nitrogen%ptr(ics:ice,:,iblk) = pool%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%nitrogen%ptr(ics:ice,:,iblk)  
    pool%residue%nitrogen%ptr(ics:ice,:,iblk) = pool%residue%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%nitrogen%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%nitrogen%ptr(ics:ice,:,iblk) = pool%residue_assoc%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%nitrogen%ptr(ics:ice,:,iblk)  

    ! phosphorus
    pool%dom%phosphorus%ptr(ics:ice,:,iblk) = pool%dom%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%phosphorus%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%phosphorus%ptr(ics:ice,:,iblk) = pool%dom_assoc%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%phosphorus%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) = pool%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%phosphorus%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk)  
    pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk) = pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%phosphorus%ptr(ics:ice,:,iblk)  
    pool%fungi%phosphorus%ptr(ics:ice,:,iblk) = pool%fungi%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%fungi%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%fungi%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%fungi%phosphorus%ptr(ics:ice,:,iblk)  
    pool%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) = pool%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk)  
    pool%microbial%phosphorus%ptr(ics:ice,:,iblk) = pool%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%phosphorus%ptr(ics:ice,:,iblk)  
    pool%residue%phosphorus%ptr(ics:ice,:,iblk) = pool%residue%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%phosphorus%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%phosphorus%ptr(ics:ice,:,iblk) = pool%residue_assoc%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%phosphorus%ptr(ics:ice,:,iblk)  

    ! C13
    pool%dom%carbon13%ptr(ics:ice,:,iblk) = pool%dom%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%carbon13%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%carbon13%ptr(ics:ice,:,iblk) = pool%dom_assoc%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%carbon13%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%carbon13%ptr(ics:ice,:,iblk) = pool%soluable_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%carbon13%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) = pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%carbon13%ptr(ics:ice,:,iblk)  
    pool%woody_litter%carbon13%ptr(ics:ice,:,iblk) = pool%woody_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%carbon13%ptr(ics:ice,:,iblk)  
    pool%fungi%carbon13%ptr(ics:ice,:,iblk) = pool%fungi%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%fungi%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%fungi%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%fungi%carbon13%ptr(ics:ice,:,iblk)  
    pool%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) = pool%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%mycorrhiza%carbon13%ptr(ics:ice,:,iblk)  
    pool%microbial%carbon13%ptr(ics:ice,:,iblk) = pool%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%carbon13%ptr(ics:ice,:,iblk)  
    pool%residue%carbon13%ptr(ics:ice,:,iblk) = pool%residue%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%carbon13%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%carbon13%ptr(ics:ice,:,iblk) = pool%residue_assoc%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%carbon13%ptr(ics:ice,:,iblk)  
    ! C14
    pool%dom%carbon14%ptr(ics:ice,:,iblk) = pool%dom%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%carbon14%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%carbon14%ptr(ics:ice,:,iblk) = pool%dom_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%carbon14%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%carbon14%ptr(ics:ice,:,iblk) = pool%soluable_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%carbon14%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%carbon14%ptr(ics:ice,:,iblk)  
    pool%woody_litter%carbon14%ptr(ics:ice,:,iblk) = pool%woody_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%carbon14%ptr(ics:ice,:,iblk)  
    pool%fungi%carbon14%ptr(ics:ice,:,iblk) = pool%fungi%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%fungi%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%fungi%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%fungi%carbon14%ptr(ics:ice,:,iblk)  
    pool%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) = pool%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%mycorrhiza%carbon14%ptr(ics:ice,:,iblk)  
    pool%microbial%carbon14%ptr(ics:ice,:,iblk) = pool%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%carbon14%ptr(ics:ice,:,iblk)  
    pool%residue%carbon14%ptr(ics:ice,:,iblk) = pool%residue%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%carbon14%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%carbon14%ptr(ics:ice,:,iblk) = pool%residue_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%carbon14%ptr(ics:ice,:,iblk)  
    ! N15
    pool%dom%nitrogen15%ptr(ics:ice,:,iblk) = pool%dom%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) = pool%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) = pool%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) = pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%fungi%nitrogen15%ptr(ics:ice,:,iblk) = pool%fungi%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%fungi%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%fungi%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%fungi%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk) = pool%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%microbial%nitrogen15%ptr(ics:ice,:,iblk) = pool%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%residue%nitrogen15%ptr(ics:ice,:,iblk) = pool%residue%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) = pool%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk)  

  END SUBROUTINE calc_organic_pool_change


END MODULE mo_sb_update_pools
