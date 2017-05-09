!>
!! @brief Contains the routines for the Task sb_simple_model 
!!
!! @author
!!  Silvia Caldararu, Sönke Zaehle 
!!
!! @par Revision History
!! First version by Sönke Zaehle (20170105): copyied from vegetation/mo_veg_soil.f90
!! Second version by Sönke Zaehle 20170331): added a 4 pool model 
!!
MODULE mo_sb_simple_model

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp, eps8, eps4
  USE mo_exception,       ONLY: message, message_text, finish
  USE mo_base_memory,     ONLY: t_jsb_task_options

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_simple_model_process 


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_simple_model'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil biogeochemical pools 
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sb_simple_model_process(tile, options)

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
    USE mo_shared_isotope
    USE mo_sb_jsm_processes

 
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
    dsl4jsb_Real2D_onChunk      :: nhx_deposition
    dsl4jsb_Real2D_onChunk      :: noy_deposition
    ! ..
    dsl4jsb_Real3D_onChunk :: t_soil_sl
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    dsl4jsb_Real3D_onChunk :: w_soil_fc_sl
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
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
    dsl4jsb_Real3D_onChunk :: nitrification 
    dsl4jsb_Real3D_onChunk :: denitrification 
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: nh4_solute 
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute 
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: no3_n15_solute 
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
    dsl4jsb_Get_var2D_onChunk(A2L_, nhx_deposition)                       ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, noy_deposition)                       ! in
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SOIL_,     t_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     w_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     w_soil_fc_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     soil_depth_sl)
    ! ---------------------------
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
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute) 


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
         !!BLARPP
         rmm_depolymerisation(:,:) = 1.0_wp
         rmm_mic_uptake       = calc_oxygen_rate_modifier(w_soil_sl,w_soil_fc_sl,kmr_mic_uptake)
         rmm_sorption         = calc_convert_bulk_to_liquid(w_soil_sl)
         rmm_desorption       = calc_convert_bulk_to_liquid(w_soil_sl) 
         rmm_hsc              = calc_substrate_rate_modifier(w_soil_sl,w_soil_fc_sl,kmr_hsc) 

    END SELECT

    CALL calc_sb_fluxes_simple(tile, options)

!    CALL calc_jsm_pool_change(options,sb%pool,sb%formation,sb%loss,sb%transport)

!    IF(sb_config%flag_n_deposition) THEN
!      nh4_solute(:,1)     = nh4_solute(:,1)   + nhx_deposition(:) / soil_depth_sl(:,1) * dtime
!      no3_solute(:,1)     = no3_solute(:,1)   + noy_deposition(:) / soil_depth_sl(:,1) * dtime
!    END IF

!    veg%uptake_nh4%ptr(ice:ice,iblk) = plant_uptake_nh4_sl(:,1) * soil_depth_sl(:,1)
!    veg%uptake_no3%ptr(ice:ice,iblk) = plant_uptake_no3_sl(:,1) * soil_depth_sl(:,1)
!    veg%uptake_p%ptr(ice:ice,iblk)   = plant_uptake_po4_sl(:,1) * soil_depth_sl(:,1)
!    IF(.NOT.sb_config%flag_sb_prescribe_nh4) THEN
!      nh4_solute(:,:)     = nh4_solute(:,:)     + (net_mineralisation_nh4(:,:) - sb%nitrification%ptr(ice:ice,:,iblk) & 
!                                                       - sb%volatilisation_n%ptr(ice:ice,:,iblk)- plant_uptake_nh4_sl(:,:)) * dtime/1.e6_wp
!      nh4_n15_solute(:,:) = nh4_n15_solute(:,:) + net_mineralisation_nh4_n15(:,:) / 1000000.0_wp * dtime
!    ENDIF
!    IF(.NOT.sb_config%flag_sb_prescribe_no3) THEN
!      no3_solute(:,:)     = no3_solute(:,:)     + (sb%nitrification%ptr(ice:ice,:,iblk) - plant_uptake_no3_sl(:,:)) / dtime/1.e6_wp 
!      no3_n15_solute(:,:) = no3_n15_solute(:,:) + net_mineralisation_no3_n15(:,:) / 1000000.0_wp * dtime
!    ENDIF
!    IF(.NOT.sb_config%flag_sb_prescribe_po4) THEN
!      po4_solute(:,:)     = po4_solute(:,:) + (net_mineralisation_po4(:,:) - plant_uptake_po4_sl(:,:)) * dtime/1.e6_wp 
!    ENDIF

  END SUBROUTINE update_sb_simple_model_process


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task
  !
  !-----------------------------------------------------------------------------------------------------
  !> routine to update soil litter and mineral pools (simple version for vegetation model testing)
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_sb_fluxes_simple(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_base_constants,        ONLY: dtime
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_soil_config_class,     ONLY: soil_config
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
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib               !< land-cover-type library - parameter across pft's
    INTEGER                               :: isoil
    INTEGER                               :: iblk, ics, ice, nc
    REAL(wp), DIMENSION(options%nc,soil_config%nsoil_water) :: apparent_mic_cue,growth_c,growth_n,hlp1
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_pools_lesssimple'
    ! ---------------------------
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(SOIL_)                         ! TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(SB_)                           ! TYPE(t_sb_memory),       POINTER :: sb__mem
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: uptake_nh4
    dsl4jsb_Real2D_onChunk :: uptake_no3
    dsl4jsb_Real2D_onChunk :: uptake_p
    dsl4jsb_Real3D_onChunk :: t_soil_sl
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    dsl4jsb_Real3D_onChunk :: w_soil_fc_sl
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
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
    dsl4jsb_Real3D_onChunk :: nitrification 
    dsl4jsb_Real3D_onChunk :: volatilisation_n
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: nh4_solute 
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute 
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: no3_n15_solute 
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
    dsl4jsb_Get_var3D_onChunk(SOIL_,     t_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     w_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     w_soil_fc_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     root_fraction_sl)
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
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_n)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute) 


!! *slmtemp
    veg  => veg__mem
    sb   => sb__mem

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------
    !> 1.0 Initialise output 
    !!
    CALL sb%loss%init(0.0_wp)
    CALL sb%formation%init(0.0_wp)
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
    volatilisation_n(:,:) = 0.0_wp
    nitrification(:,:) = 0.0_wp


    !> 1.1 Communicate vegetation uptake to soil layers (BLARPPP - temporary for now)
    plant_uptake_nh4_sl(:,1) = uptake_nh4(:) / soil_depth_sl(:,1) * dtime/1.e6 ! convert to mol/timestep
    !plant_uptake_no3_sl(:,1) = uptake_no3(:) / soil_depth_sl(:,1)
    plant_uptake_po4_sl(:,1) = uptake_p(:) / soil_depth_sl(:,1) * dtime/1.e6 ! convert to mol/timestep

    !> 2.0 Calculate decomposition of pools: dP/dt = - f(T)*f(M)/tau_pool * P
    !!
    !> 2.1 non-woody litter (= polymeric pool)
    !! 
    sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) = sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + &
                                                          rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_non_woody_litter * dtime * &
                                                          sb%pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk)
    sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) = sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + &
                                                            rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_non_woody_litter * dtime * &
                                                            sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)
    sb%loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = sb%loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + &
                                                              rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_non_woody_litter * dtime * &
                                                              sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk)
    sb%loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) = sb%loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + &
                                                            rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_non_woody_litter * dtime * &
                                                            sb%pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk)
    sb%loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = sb%loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + &
                                                            rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_non_woody_litter * dtime * &
                                                            sb%pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk)
    sb%loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb%loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + &
                                                              rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_non_woody_litter * dtime * &
                                                              sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)

    !> 2.2 woody litter
    !! 
    sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk) = sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk) + &
                                                      rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_woody_litter * dtime * &
                                                      sb%pool%woody_litter%carbon%ptr(ics:ice,:,iblk)
    sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) = sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + &
                                                        rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_woody_litter * dtime * &
                                                        sb%pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk)
    sb%loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) = sb%loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) + &
                                                          rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_woody_litter * dtime * &
                                                          sb%pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk)
    sb%loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) = sb%loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) + &
                                                        rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_woody_litter * dtime * &
                                                        sb%pool%woody_litter%carbon13%ptr(ics:ice,:,iblk)
    sb%loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) = sb%loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) + &
                                                        rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_woody_litter * dtime * &
                                                        sb%pool%woody_litter%carbon14%ptr(ics:ice,:,iblk)
    sb%loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb%loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) + &
                                                          rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_woody_litter * dtime * &
                                                          sb%pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk)

    !> 2.3 fast SOM pool (=microbial pool)
    !!
    sb%loss%microbial%carbon%ptr(ics:ice,:,iblk) = sb%loss%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                   rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_fast * dtime * &
                                                   sb%pool%microbial%carbon%ptr(ics:ice,:,iblk)
    sb%loss%microbial%nitrogen%ptr(ics:ice,:,iblk) = sb%loss%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_fast * dtime * &
                                                     sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk)
    sb%loss%microbial%phosphorus%ptr(ics:ice,:,iblk) = sb%loss%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                       rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_fast * dtime * &
                                                       sb%pool%microbial%phosphorus%ptr(ics:ice,:,iblk)
    sb%loss%microbial%carbon13%ptr(ics:ice,:,iblk) = sb%loss%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_fast * dtime * &
                                                     sb%pool%microbial%carbon13%ptr(ics:ice,:,iblk)
    sb%loss%microbial%carbon14%ptr(ics:ice,:,iblk) = sb%loss%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_fast * dtime * &
                                                     sb%pool%microbial%carbon14%ptr(ics:ice,:,iblk)
    sb%loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) = sb%loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                       rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_fast * dtime * &
                                                       sb%pool%microbial%nitrogen15%ptr(ics:ice,:,iblk)

    !> 2.4 slow SOM pool (=residue pool)
    !!
    sb%loss%residue%carbon%ptr(ics:ice,:,iblk) = sb%loss%residue%carbon%ptr(ics:ice,:,iblk) + & 
                                                 rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_slow * dtime * &
                                                 sb%pool%residue%carbon%ptr(ics:ice,:,iblk)
    sb%loss%residue%nitrogen%ptr(ics:ice,:,iblk) = sb%loss%residue%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                   rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_slow * dtime * &
                                                   sb%pool%residue%nitrogen%ptr(ics:ice,:,iblk)
    sb%loss%residue%phosphorus%ptr(ics:ice,:,iblk) = sb%loss%residue%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_slow * dtime * &
                                                     sb%pool%residue%phosphorus%ptr(ics:ice,:,iblk)
    sb%loss%residue%carbon13%ptr(ics:ice,:,iblk) = sb%loss%residue%carbon13%ptr(ics:ice,:,iblk) + & 
                                                   rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_slow * dtime * &
                                                   sb%pool%residue%carbon13%ptr(ics:ice,:,iblk)
    sb%loss%residue%carbon14%ptr(ics:ice,:,iblk) = sb%loss%residue%carbon14%ptr(ics:ice,:,iblk) + & 
                                                   rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_slow * dtime * &
                                                   sb%pool%residue%carbon14%ptr(ics:ice,:,iblk)
    sb%loss%residue%nitrogen15%ptr(ics:ice,:,iblk) = sb%loss%residue%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation(:,:)*rmm_depolymerisation(:,:)/tau_slow * dtime * &
                                                     sb%pool%residue%nitrogen15%ptr(ics:ice,:,iblk)

    !> 3.0 formation of SOM resulting from the decomposition dPrec/dt = CUE_Pdon * dPdon/dt
    !! net mineralisation follows the stoichiometric difference of the donor and recipient pool, taking account for CUE
    !! for fast -> slow -> fast this only causes positive mineralisation because of the C:N:P and CUE
    !! for litter -> fast this causes in most circumstances negative mineralisation, which is then taken from the 
    !!               mineral pools (or organic matter for P)

    !> 3.1 formation of fast SOM from litter 
    !! 
    ! calculate C growth of fast pool
    growth_c(:,:) = microbial_cue * (sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk))
    
    ! calculate potential immobilisation rate 
    hlp1(:,:) = growth_c(:,:) * microbial_cn - (sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk))
    ! limit to available N
    WHERE(hlp1(:,:) + plant_uptake_nh4_sl(:,:) > nh4_solute(:,:) .AND. hlp1(:,:) + plant_uptake_nh4_sl(:,:) > 1.e-13)
        hlp1(:,:) = nh4_solute(:,:) * hlp1(:,:) / (hlp1(:,:) + plant_uptake_nh4_sl(:,:)) 
        plant_uptake_nh4_sl(:,:) = nh4_solute(:,:) * plant_uptake_nh4_sl(:,:) / (hlp1(:,:) + plant_uptake_nh4_sl(:,:)) 
    ENDWHERE 
    
    ! calculate C growth based on N availability
    growth_n(:,:) = microbial_cn * (sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                    sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + hlp1(:,:))

    ! adjust carbon-use efficiency to match N availability by unchanged litter decomposition
    WHERE(growth_c(:,:) <= growth_n(:,:)) ! enough N is available -> constant CUE
       apparent_mic_cue(:,:) = microbial_cue
    ELSEWHERE(growth_c(:,:) > eps8) ! not enough N is available -> changed CUE 
       apparent_mic_cue(:,:) = growth_n * microbial_cue / growth_c
    ELSEWHERE ! growth C is too small -> CUE does not matter
       apparent_mic_cue(:,:) = 0.0_wp
    ENDWHERE

    sb%formation%microbial%carbon%ptr(ics:ice,:,iblk) = sb%formation%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                        apparent_mic_cue(:,:) * & 
                                                        (sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) 
    sb%formation%microbial%nitrogen%ptr(ics:ice,:,iblk) = sb%formation%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                          apparent_mic_cue(:,:) / microbial_cn * & 
                                                          (sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) 
    sb%formation%microbial%phosphorus%ptr(ics:ice,:,iblk) = sb%formation%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                            apparent_mic_cue(:,:) / microbial_cn / microbial_np * & 
                                                            (sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) 
    sb%formation%microbial%carbon13%ptr(ics:ice,:,iblk) = sb%formation%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                          apparent_mic_cue(:,:) * & 
                                                          (sb%loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%carbon13%ptr(ics:ice,:,iblk)) 
    sb%formation%microbial%carbon14%ptr(ics:ice,:,iblk) = sb%formation%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                          apparent_mic_cue(:,:) * & 
                                                          (sb%loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%carbon14%ptr(ics:ice,:,iblk)) 
    WHERE(sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)+sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk)>eps8)
       sb%formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) = sb%formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                              (sb%loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)+sb%loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk)) / & 
                              (sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)+sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk)) * & 
                               apparent_mic_cue(:,:) / microbial_cn * & 
                              (sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) 
    ENDWHERE

    het_respiration(:,:)     = het_respiration(:,:) + (1._wp - apparent_mic_cue(:,:)) * & 
                               (sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) 
    het_respiration_c13(:,:) = het_respiration_c13(:,:) + (1._wp - apparent_mic_cue(:,:)) * & 
                               (sb%loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%carbon13%ptr(ics:ice,:,iblk)) 
    het_respiration_c14(:,:) = het_respiration_c14(:,:) + (1._wp - apparent_mic_cue(:,:)) * & 
                               (sb%loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%carbon14%ptr(ics:ice,:,iblk)) 
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) + & 
                                  (sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk)) - & 
                                  apparent_mic_cue(:,:) / microbial_cn * & 
                                  (sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) 
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) + & 
                                  (sb%loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk)) - & 
                                  apparent_mic_cue(:,:) / microbial_cn / microbial_np * & 
                                  (sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) 
    WHERE(sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)+sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk)>eps8)
       net_mineralisation_nh4_n15(:,:)     = net_mineralisation_nh4_n15(:,:) + & 
                              sb%loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) - & 
                              (sb%loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)+sb%loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk)) / & 
                              (sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)+sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk)) * & 
                               apparent_mic_cue(:,:) / microbial_cn * & 
                              (sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) 
    ENDWHERE

    !> 3.2 formation of slow SOM (=residue pool) from fast SOM decay
    !! 
    sb%formation%residue%carbon%ptr(ics:ice,:,iblk) = sb%formation%residue%carbon%ptr(ics:ice,:,iblk) + & 
                                                      residue_cue * sb%loss%microbial%carbon%ptr(ics:ice,:,iblk)
    sb%formation%residue%nitrogen%ptr(ics:ice,:,iblk) = sb%formation%residue%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                        residue_cue / microbial_cn * sb%loss%microbial%carbon%ptr(ics:ice,:,iblk)
    sb%formation%residue%phosphorus%ptr(ics:ice,:,iblk) = sb%formation%residue%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                          residue_cue / microbial_cn / microbial_np * sb%loss%microbial%carbon%ptr(ics:ice,:,iblk)
    sb%formation%residue%carbon13%ptr(ics:ice,:,iblk) = sb%formation%residue%carbon13%ptr(ics:ice,:,iblk) + & 
                                                        residue_cue * sb%loss%microbial%carbon13%ptr(ics:ice,:,iblk)
    sb%formation%residue%carbon14%ptr(ics:ice,:,iblk) = sb%formation%residue%carbon14%ptr(ics:ice,:,iblk) + & 
                                                        residue_cue * sb%loss%microbial%carbon14%ptr(ics:ice,:,iblk)
    WHERE(sb%loss%microbial%nitrogen%ptr(ics:ice,:,iblk) > eps8)
       sb%formation%residue%nitrogen15%ptr(ics:ice,:,iblk) = sb%formation%residue%nitrogen15%ptr(ics:ice,:,iblk) + &
                                                             sb%loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) / &  
                                                             sb%loss%microbial%nitrogen%ptr(ics:ice,:,iblk) * & 
                                                             residue_cue / microbial_cn * sb%loss%microbial%carbon%ptr(ics:ice,:,iblk)
    ENDWHERE
    het_respiration(:,:)     = het_respiration(:,:)     + (1._wp - residue_cue) * sb%loss%microbial%carbon%ptr(ics:ice,:,iblk)
    het_respiration_c13(:,:) = het_respiration_c13(:,:) + (1._wp - residue_cue) * sb%loss%microbial%carbon13%ptr(ics:ice,:,iblk)
    het_respiration_c14(:,:) = het_respiration_c14(:,:) + (1._wp - residue_cue) * sb%loss%microbial%carbon14%ptr(ics:ice,:,iblk)
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) + sb%loss%microbial%nitrogen%ptr(ics:ice,:,iblk) - &
                                                                residue_cue / microbial_cn * sb%loss%microbial%carbon%ptr(ics:ice,:,iblk)
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) + sb%loss%microbial%phosphorus%ptr(ics:ice,:,iblk) - &
                                                                residue_cue / microbial_cn / microbial_np * sb%loss%microbial%carbon%ptr(ics:ice,:,iblk)
    WHERE(sb%loss%microbial%nitrogen%ptr(ics:ice,:,iblk) > eps8)
       net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) + & 
                                    sb%loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) - &
                                    sb%loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) / sb%loss%microbial%nitrogen%ptr(ics:ice,:,iblk) * & 
                                    residue_cue / microbial_cn * sb%loss%microbial%carbon%ptr(ics:ice,:,iblk)
    ENDWHERE

    !> 3.3 formation of fast SOM (=microbial pool) from slow SOM decay 
    !! 
    sb%formation%microbial%carbon%ptr(ics:ice,:,iblk) = sb%formation%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                        residue_cue * sb%loss%residue%carbon%ptr(ics:ice,:,iblk)
    sb%formation%microbial%nitrogen%ptr(ics:ice,:,iblk) = sb%formation%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                          residue_cue / microbial_cn * sb%loss%residue%carbon%ptr(ics:ice,:,iblk)
    sb%formation%microbial%phosphorus%ptr(ics:ice,:,iblk) = sb%formation%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                            residue_cue / microbial_cn / microbial_np * sb%loss%residue%carbon%ptr(ics:ice,:,iblk)
    sb%formation%microbial%carbon13%ptr(ics:ice,:,iblk) = sb%formation%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                          residue_cue * sb%loss%residue%carbon13%ptr(ics:ice,:,iblk)
    sb%formation%microbial%carbon14%ptr(ics:ice,:,iblk) = sb%formation%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                          residue_cue * sb%loss%residue%carbon14%ptr(ics:ice,:,iblk)
    WHERE(sb%loss%microbial%nitrogen%ptr(ics:ice,:,iblk)>eps8)
       sb%formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) = sb%formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                               sb%loss%residue%nitrogen15%ptr(ics:ice,:,iblk) / &  
                                                               sb%loss%residue%nitrogen%ptr(ics:ice,:,iblk) * &  
                                                               residue_cue / microbial_cn * sb%loss%residue%carbon%ptr(ics:ice,:,iblk)
    ENDWHERE

    het_respiration(:,:)      = het_respiration(:,:)     + (1._wp - residue_cue) * sb%loss%residue%carbon%ptr(ics:ice,:,iblk)
    het_respiration_c13(:,:)  = het_respiration_c13(:,:) + (1._wp - residue_cue) * sb%loss%residue%carbon13%ptr(ics:ice,:,iblk)
    het_respiration_c14(:,:)  = het_respiration_c14(:,:) + (1._wp - residue_cue) * sb%loss%residue%carbon14%ptr(ics:ice,:,iblk)
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) + sb%loss%residue%nitrogen%ptr(ics:ice,:,iblk) - &
                                                                residue_cue / microbial_cn * sb%loss%residue%carbon%ptr(ics:ice,:,iblk)
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) + sb%loss%residue%phosphorus%ptr(ics:ice,:,iblk) - &
                                                                residue_cue / microbial_cn / microbial_np * sb%loss%residue%carbon%ptr(ics:ice,:,iblk)
    WHERE(sb%loss%microbial%nitrogen%ptr(ics:ice,:,iblk)>eps8)
       net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) + & 
                                    sb%loss%residue%nitrogen15%ptr(ics:ice,:,iblk) - &
                                    sb%loss%residue%nitrogen15%ptr(ics:ice,:,iblk) / sb%loss%residue%nitrogen%ptr(ics:ice,:,iblk) * &  
                                    residue_cue / microbial_cn * sb%loss%residue%carbon%ptr(ics:ice,:,iblk)
    ENDWHERE

    !> 4.0 calculate biomineralisation
    !! calculate potential P immobilisation rate 
    hlp1(:,:) = ABS(net_mineralisation_po4(:,:))
    ! limit plant uptake to avaialble P 
    WHERE(net_mineralisation_po4(:,:) < 0.0_wp .AND. & 
          hlp1(:,:) + plant_uptake_po4_sl(:,:) > po4_solute(:,:) .AND. hlp1(:,:) + plant_uptake_po4_sl(:,:) > 1.e-13)
        plant_uptake_po4_sl(:,:) = po4_solute(:,:) * plant_uptake_po4_sl(:,:) / (hlp1(:,:) + plant_uptake_po4_sl(:,:)) 
    ENDWHERE 
    
    ! limit to available P 
    hlp1(:,:) = ABS(net_mineralisation_po4(:,:)) - (po4_solute(:,:) - plant_uptake_po4_sl(:,:))
    WHERE(net_mineralisation_po4(:,:) < 0.0_wp .AND. hlp1(:,:) > 0.0_wp) 
          net_mineralisation_po4(:,:) = - (po4_solute(:,:) - plant_uptake_po4_sl(:,:))
          sb%loss%residue%phosphorus%ptr(ics:ice,:,iblk) = sb%loss%residue%phosphorus%ptr(ics:ice,:,iblk) + hlp1
    ENDWHERE

    !> 4.1 calculate N losses and nitrification
    !! 
    WHERE(net_mineralisation_nh4(:,:)>0.0_wp)
       volatilisation_n(:,:) = net_mineralisation_nh4(:,:) * floss_nmin 
       nitrification(:,:)    = net_mineralisation_nh4(:,:) * (1._wp - floss_nmin) * fnit_nmin 
    ENDWHERE

    !> 5.0 Update soil litter pools with litter fall
    !!
    sb%formation%polymeric_litter%carbon%ptr(ics:ice,1,iblk) = sb%formation%polymeric_litter%carbon%ptr(ics:ice,1,iblk) + & 
       (veg%litterfall%leaf%carbon%ptr + veg%litterfall%fruit%carbon%ptr + & 
        veg%litterfall%labile%carbon%ptr + veg%litterfall%reserve%carbon%ptr + & 
        veg%litterfall%seed_bed%carbon%ptr ) / soil_depth_sl(:,1)
    sb%formation%polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) = sb%formation%polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) +  & 
       (veg%litterfall%leaf%nitrogen%ptr + veg%litterfall%fruit%nitrogen%ptr + & 
        veg%litterfall%labile%nitrogen%ptr + veg%litterfall%reserve%nitrogen%ptr + & 
        veg%litterfall%seed_bed%nitrogen%ptr ) / soil_depth_sl(:,1)
    sb%formation%polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) = sb%formation%polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) + & 
       (veg%litterfall%leaf%phosphorus%ptr + veg%litterfall%fruit%phosphorus%ptr + & 
        veg%litterfall%labile%phosphorus%ptr + veg%litterfall%reserve%phosphorus%ptr + & 
        veg%litterfall%seed_bed%phosphorus%ptr ) / soil_depth_sl(:,1)
    sb%formation%polymeric_litter%carbon13%ptr(ics:ice,1,iblk) = sb%formation%polymeric_litter%carbon13%ptr(ics:ice,1,iblk) + & 
       (veg%litterfall%leaf%carbon13%ptr + veg%litterfall%fruit%carbon13%ptr + & 
        veg%litterfall%labile%carbon13%ptr + veg%litterfall%reserve%carbon13%ptr + & 
        veg%litterfall%seed_bed%carbon13%ptr ) / soil_depth_sl(:,1)
    sb%formation%polymeric_litter%carbon14%ptr(ics:ice,1,iblk) = sb%formation%polymeric_litter%carbon14%ptr(ics:ice,1,iblk) + & 
       (veg%litterfall%leaf%carbon14%ptr + veg%litterfall%fruit%carbon14%ptr + & 
        veg%litterfall%labile%carbon14%ptr + veg%litterfall%reserve%carbon14%ptr + & 
        veg%litterfall%seed_bed%carbon14%ptr ) / soil_depth_sl(:,1)
    sb%formation%polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) = sb%formation%polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) + &
       (veg%litterfall%leaf%nitrogen15%ptr + veg%litterfall%fruit%nitrogen15%ptr + & 
        veg%litterfall%labile%nitrogen15%ptr + veg%litterfall%reserve%nitrogen15%ptr + & 
        veg%litterfall%seed_bed%nitrogen15%ptr ) / soil_depth_sl(:,1)

    sb%formation%woody_litter%carbon%ptr(ics:ice,1,iblk) = sb%formation%woody_litter%carbon%ptr(ics:ice,1,iblk) + & 
       (veg%litterfall%sap_wood%carbon%ptr + veg%litterfall%heart_wood%carbon%ptr) / & 
       soil_depth_sl(:,1)
    sb%formation%woody_litter%nitrogen%ptr(ics:ice,1,iblk) = sb%formation%woody_litter%nitrogen%ptr(ics:ice,1,iblk) +  & 
       (veg%litterfall%sap_wood%nitrogen%ptr + veg%litterfall%heart_wood%nitrogen%ptr) / & 
       soil_depth_sl(:,1)
    sb%formation%woody_litter%phosphorus%ptr(ics:ice,1,iblk) = sb%formation%woody_litter%phosphorus%ptr(ics:ice,1,iblk) + & 
       (veg%litterfall%sap_wood%phosphorus%ptr + veg%litterfall%heart_wood%phosphorus%ptr) / & 
       soil_depth_sl(:,1)
    sb%formation%woody_litter%carbon13%ptr(ics:ice,1,iblk) = sb%formation%woody_litter%carbon13%ptr(ics:ice,1,iblk) + & 
       (veg%litterfall%sap_wood%carbon13%ptr + veg%litterfall%heart_wood%carbon13%ptr) / & 
       soil_depth_sl(:,1)
    sb%formation%woody_litter%carbon14%ptr(ics:ice,1,iblk) = sb%formation%woody_litter%carbon14%ptr(ics:ice,1,iblk) + & 
       (veg%litterfall%sap_wood%carbon14%ptr + veg%litterfall%heart_wood%carbon14%ptr) / & 
       soil_depth_sl(:,1)
    sb%formation%woody_litter%nitrogen15%ptr(ics:ice,1,iblk) = sb%formation%woody_litter%nitrogen15%ptr(ics:ice,1,iblk) + &
       (veg%litterfall%sap_wood%nitrogen15%ptr + veg%litterfall%heart_wood%nitrogen15%ptr) / &
       soil_depth_sl(:,1)

    DO isoil = 1,1

       sb%pool%polymeric_litter%carbon%ptr(ics:ice,isoil,iblk) = sb%pool%polymeric_litter%carbon%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%carbon%ptr + veg%litterfall%coarse_root%carbon%ptr ) / & ! * & 
          !root_fraction_sl(:,isoil) / & 
            soil_depth_sl(:,isoil) 
       sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,isoil,iblk) = sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%nitrogen%ptr + veg%litterfall%coarse_root%nitrogen%ptr ) / & ! * & 
          !root_fraction_sl(:,isoil) / & 
            soil_depth_sl(:,isoil) 
       sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,isoil,iblk) = sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%phosphorus%ptr + veg%litterfall%coarse_root%phosphorus%ptr ) / & ! * & 
          !root_fraction_sl(:,isoil) / & 
            soil_depth_sl(:,isoil) 
       sb%pool%polymeric_litter%carbon13%ptr(ics:ice,isoil,iblk) = sb%pool%polymeric_litter%carbon13%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%carbon13%ptr + veg%litterfall%coarse_root%carbon13%ptr ) / & ! * & 
          !root_fraction_sl(:,isoil) / & 
            soil_depth_sl(:,isoil) 
       sb%pool%polymeric_litter%carbon14%ptr(ics:ice,isoil,iblk) = sb%pool%polymeric_litter%carbon14%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%carbon14%ptr + veg%litterfall%coarse_root%carbon14%ptr ) / & ! * & 
          !root_fraction_sl(:,isoil) / & 
            soil_depth_sl(:,isoil) 
       sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,isoil,iblk) = sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%nitrogen15%ptr + veg%litterfall%coarse_root%nitrogen15%ptr ) / & ! * & 
          !root_fraction_sl(:,isoil) / & 
            soil_depth_sl(:,isoil) 

    ENDDO

    !> 6.0 convert rates to µmol/ m3 / s
    !!
    het_respiration(:,:) = het_respiration(:,:) / dtime * 1000000._wp
    het_respiration_c13(:,:) = het_respiration_c13(:,:) / dtime * 1000000._wp
    het_respiration_c14(:,:) = het_respiration_c14(:,:) / dtime * 1000000._wp
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) / dtime * 1000000._wp
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) / dtime * 1000000._wp
    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) / dtime * 1000000._wp
    plant_uptake_nh4_sl(:,:) = plant_uptake_nh4_sl(:,:) / dtime * 1000000._wp 
    plant_uptake_po4_sl(:,:) = plant_uptake_po4_sl(:,:) / dtime * 1000000._wp 
    volatilisation_n(:,:) = volatilisation_n(:,:) / dtime * 1000000._wp
    nitrification(:,:) = nitrification(:,:) / dtime * 1000000._wp 


  END SUBROUTINE calc_sb_fluxes_simple 

END MODULE mo_sb_simple_model
