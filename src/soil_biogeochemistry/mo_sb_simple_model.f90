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

  USE mo_kind,            ONLY: wp, eps8, eps4
  USE mo_exception,       ONLY: message, message_text, finish

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_simple_model 


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
  SUBROUTINE update_sb_simple_model(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout), TARGET  :: tile         !< one tile
    TYPE(t_jsb_task_options), INTENT(in)             :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    !TYPE(t_jsb_model),      POINTER       :: model
    !TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    !INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_simple_model'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...





    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    CALL calc_sb_fluxes_simple(tile, options)


  END SUBROUTINE update_sb_simple_model


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
    USE mo_base_constants,        ONLY: dtime, one_day
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_soil_config_class,     ONLY: soil_config
    USE mo_sb_config_class,       ONLY: sb_config
    USE mo_sb_constants
    USE mo_sb_jsm_transport

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
    INTEGER                               :: isoil                !< loop over soil_config%nsoil_water
    INTEGER                               :: iblk, ics, ice, nc
    REAL(wp), DIMENSION(options%nc,soil_config%nsoil_water) ::      &
                                             apparent_mic_cue     , &
                                             growth_c             , &
                                             growth_n             , &
                                             hlp1                 , &
                                             hlp2                 , &
                                             hlp3                 , &
                                             mic_uptake_sl        , &
                                             root_fraction_act_sl , &      ! and below needed for switch between 0d and 1d
                                             rtm_depolymerisation_act, &
                                             rmm_depolymerisation_act, & 
                                             km1_plant_uptake_n, &
                                             km2_plant_uptake_n, &
                                             km1_plant_uptake_p, &
                                             km2_plant_uptake_p, &
                                             km_denitrification_c_act, &
                                             km_denitrification_no3_act, &
                                             rtm_nitrification_act, &
                                             rtm_denitrification_act, &
                                             rmm_nitrification_act, &
                                             rtm_gasdiffusion_act, &
                                             rmm_gasdiffusion_act, &
                                             anaerobic_volume_fraction_act
    REAL(wp), DIMENSION(options%nc) :: root_depth_avg
    
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_pools_lesssimple'
    ! ---------------------------
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(SOIL_)                         ! TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(SB_)                           ! TYPE(t_sb_memory),       POINTER :: sb__mem
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: unit_uptake_n_pot
    dsl4jsb_Real2D_onChunk :: unit_uptake_p_pot
    dsl4jsb_Real2D_onChunk :: drainage_fraction 
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    dsl4jsb_Real3D_onChunk :: percolation_sl
    dsl4jsb_Real3D_onChunk :: rtm_depolymerisation 
    dsl4jsb_Real3D_onChunk :: rtm_nitrification
    dsl4jsb_Real3D_onChunk :: rtm_denitrification
    dsl4jsb_Real3D_onChunk :: rtm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: rmm_depolymerisation 
    dsl4jsb_Real3D_onChunk :: rmm_mic_uptake
    dsl4jsb_Real3D_onChunk :: rmm_sorption
    dsl4jsb_Real3D_onChunk :: rmm_desorption
    dsl4jsb_Real3D_onChunk :: rmm_hsc
    dsl4jsb_Real3D_onChunk :: rmm_nitrification
    dsl4jsb_Real3D_onChunk :: rmm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: anaerobic_volume_fraction_sl 
    dsl4jsb_Real3D_onChunk :: het_respiration
    dsl4jsb_Real3D_onChunk :: het_respiration_c13
    dsl4jsb_Real3D_onChunk :: het_respiration_c14
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: nitrification_no3 
    dsl4jsb_Real3D_onChunk :: nitrification_noy 
    dsl4jsb_Real3D_onChunk :: nitrification_n2o
    dsl4jsb_Real3D_onChunk :: denitrification_noy 
    dsl4jsb_Real3D_onChunk :: denitrification_n2o 
    dsl4jsb_Real3D_onChunk :: denitrification_n2 
    dsl4jsb_Real3D_onChunk :: volatilisation_n
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: transport_nh4_solute
    dsl4jsb_Real3D_onChunk :: transport_no3_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_solute
    dsl4jsb_Real2D_onChunk :: leaching_nh4_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_solute
    dsl4jsb_Real2D_onChunk :: leaching_po4_solute
    dsl4jsb_Real3D_onChunk :: transport_noy
    dsl4jsb_Real3D_onChunk :: transport_n2o
    dsl4jsb_Real3D_onChunk :: transport_n2 
    dsl4jsb_Real3D_onChunk :: noy
    dsl4jsb_Real3D_onChunk :: n2o
    dsl4jsb_Real3D_onChunk :: n2 
    dsl4jsb_Real3D_onChunk :: nh4_solute 
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute 
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: no3_n15_solute 
    dsl4jsb_Real2D_onChunk :: emission_noy
    dsl4jsb_Real2D_onChunk :: emission_n2o
    dsl4jsb_Real2D_onChunk :: emission_n2
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
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_pot)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_pot)
    dsl4jsb_Get_var2D_onChunk(SOIL_,     drainage_fraction)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     root_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     percolation_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_depolymerisation)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_nitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_denitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_gasdiffusion)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_depolymerisation)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_mic_uptake)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_sorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_desorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_hsc)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_nitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_gasdiffusion)
    dsl4jsb_Get_var3D_onChunk(SB_,       anaerobic_volume_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c13)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c14)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2)
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_n)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       noy) 
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o) 
    dsl4jsb_Get_var3D_onChunk(SB_,       n2) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_po4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy)          
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o)         
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2)         


!! *slmtemp
    veg  => veg__mem
    sb   => sb__mem

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> differentiate between simple_0d & simple_1d (1 & many soil layers)
    !!
    !! for the OD model, all processes happen only in the first layer of the soil
    !! therefore the soil temperature and moisture modifiers are aggregated 
    !! across the profile (assuming that the weighting follows the prescribed root profile)
    !! The km's are expressed as mol/m3 or m3/mol and therefore need to be adjusted to avoid
    !! the implicit assumption that all material would be concentrated in the top layer 
    !! this is done here by assuming the material is distributed within the average rooting depth
    IF(sb_config%sb_model_scheme == 'simple_0d') THEN
      root_fraction_act_sl(:,1)                         = 1.0_wp
      root_fraction_act_sl(:,2:soil_config%nsoil_water) = 0.0_wp
      ! set to default values to avoid division by zero below (2:nsoil_water are unused in OD)
      km1_plant_uptake_n(:,:) = lctlib%soil_n_aff1 
      km1_plant_uptake_p(:,:) = lctlib%soil_p_aff1 
      km2_plant_uptake_n(:,:) = lctlib%soil_n_aff2 
      km2_plant_uptake_p(:,:) = lctlib%soil_p_aff2
      km_denitrification_c_act(:,:) = km_denitrification_c
      km_denitrification_no3_act(:,:) = km_denitrification_no3
      root_depth_avg(:) = sum(root_fraction_sl(:,:)*soil_depth_sl(:,:),DIM=2)
      WHERE(root_depth_avg(:)>eps8) 
         km1_plant_uptake_n(:,1) = lctlib%soil_n_aff1 * soil_depth_sl(:,1) / root_depth_avg(:)
         km1_plant_uptake_p(:,1) = lctlib%soil_p_aff1 * soil_depth_sl(:,1) / root_depth_avg(:)
         km2_plant_uptake_n(:,1) = lctlib%soil_n_aff2 / soil_depth_sl(:,1) * root_depth_avg(:)
         km2_plant_uptake_p(:,1) = lctlib%soil_p_aff2 / soil_depth_sl(:,1) * root_depth_avg(:)
         km_denitrification_c_act(:,1) = km_denitrification_c / soil_depth_sl(:,1) * root_depth_avg(:)
         km_denitrification_no3_act(:,1) = km_denitrification_no3 / soil_depth_sl(:,1) * root_depth_avg(:)
      ENDWHERE
      rtm_depolymerisation_act(:,:) = 0.0_wp 
      rmm_depolymerisation_act(:,:) = 0.0_wp 
      rtm_nitrification_act(:,:) = 0.0_wp 
      rmm_nitrification_act(:,:) = 0.0_wp 
      rtm_denitrification_act(:,:) = 0.0_wp 
      rtm_gasdiffusion_act(:,:) = 0.0_wp 
      rmm_gasdiffusion_act(:,:) = 0.0_wp 
      anaerobic_volume_fraction_act(:,:) = 0.0_wp 
      DO isoil = 1,soil_config%nsoil_water
         rtm_depolymerisation_act(:,1) = rtm_depolymerisation_act(:,1) + & 
                                         rtm_depolymerisation(:,isoil) * root_fraction_sl(:,isoil)
         rmm_depolymerisation_act(:,1) = rmm_depolymerisation_act(:,1) + & 
                                         rmm_depolymerisation(:,isoil) * root_fraction_sl(:,isoil)
         rtm_nitrification_act(:,1) = rtm_nitrification_act(:,1) + & 
                                         rtm_nitrification(:,isoil) * root_fraction_sl(:,isoil)
         rmm_nitrification_act(:,1) = rmm_nitrification_act(:,1) + & 
                                         rmm_nitrification(:,isoil) * root_fraction_sl(:,isoil)
         rtm_denitrification_act(:,1) = rtm_denitrification_act(:,1) + & 
                                         rtm_denitrification(:,isoil) * root_fraction_sl(:,isoil)
         rtm_gasdiffusion_act(:,1) = rtm_gasdiffusion_act(:,1) + & 
                                         rtm_gasdiffusion(:,isoil) * root_fraction_sl(:,isoil)
         rmm_gasdiffusion_act(:,1) = rmm_gasdiffusion_act(:,1) + & 
                                         rmm_gasdiffusion(:,isoil) * root_fraction_sl(:,isoil)
         anaerobic_volume_fraction_act(:,1) = anaerobic_volume_fraction_act(:,1) + & 
                                         anaerobic_volume_fraction_sl(:,isoil) * root_fraction_sl(:,isoil)
      ENDDO
    ELSE ! simple_1d, simply copy values
      root_fraction_act_sl(:,:) = root_fraction_sl(:,:)
      rtm_depolymerisation_act(:,:) = rtm_depolymerisation(:,:)
      rmm_depolymerisation_act(:,:) = rmm_depolymerisation(:,:) 
      rtm_nitrification_act(:,:) = rtm_nitrification(:,:)
      rmm_nitrification_act(:,:) = rmm_nitrification(:,:) 
      rtm_denitrification_act(:,:) = rtm_denitrification(:,:)
      rtm_gasdiffusion_act(:,:) = rtm_gasdiffusion(:,:)
      rmm_gasdiffusion_act(:,:) = rmm_gasdiffusion(:,:) 
      anaerobic_volume_fraction_act(:,:) = anaerobic_volume_fraction_sl(:,:) 
      km1_plant_uptake_n(:,:) = lctlib%soil_n_aff1 
      km1_plant_uptake_p(:,:) = lctlib%soil_p_aff1 
      km2_plant_uptake_n(:,:) = lctlib%soil_n_aff2 
      km2_plant_uptake_p(:,:) = lctlib%soil_p_aff2 
      km_denitrification_c_act(:,:) = km_denitrification_c
      km_denitrification_no3_act(:,:) = km_denitrification_no3
    ENDIF

    !> 1.0 Plant nutrient uptake rates
    !!
    unit_uptake_n_pot(:) = 0.0_wp
    unit_uptake_p_pot(:) = 0.0_wp
    DO isoil = 1,soil_config%nsoil_water
      !> 1.1 inorganic N uptake 
      !!
      !! potential N uptake rate (µmol N / mol C fine-root / s), weighted by root distribution
      hlp1(:,isoil) = root_fraction_act_sl(:,isoil) * lctlib%vmax_uptake_n * & 
                  nh4_solute(:,isoil) * ( km1_plant_uptake_n(:,isoil) + & 
                  1._wp / ( km2_plant_uptake_n(:,isoil) + nh4_solute(:,isoil) )) 
      hlp2(:,isoil) = root_fraction_act_sl(:,isoil) * lctlib%vmax_uptake_n * & 
                  no3_solute(:,isoil) * ( km1_plant_uptake_n(:,isoil) + & 
                  1._wp / ( km2_plant_uptake_n(:,isoil) + no3_solute(:,isoil) )) 
      unit_uptake_n_pot(:) = unit_uptake_n_pot(:) + hlp1(:,isoil) + hlp2(:,isoil) 
  
      ! potential uptake rate in mol N / m3 / timiestep
      plant_uptake_nh4_sl(:,isoil) = hlp1(:,isoil) * veg%pool%fine_root%carbon%ptr / & 
                                     soil_depth_sl(:,isoil) * dtime/1.e6
      plant_uptake_no3_sl(:,isoil) = hlp2(:,isoil) * veg%pool%fine_root%carbon%ptr / & 
                                     soil_depth_sl(:,isoil) * dtime/1.e6
  
      !> 1.2 inorganic P uptake 
      !!
      !! potential P uptake rate (µmol P / mol C fine-root / s), weighted by root distribution
      hlp1(:,isoil) = root_fraction_act_sl(:,isoil) * lctlib%vmax_uptake_p * &
                             po4_solute(:,isoil) * ( km1_plant_uptake_p(:,isoil) + &
                             1._wp / ( km2_plant_uptake_p(:,isoil) + po4_solute(:,isoil) ))
      unit_uptake_p_pot(:) = unit_uptake_p_pot(:) + hlp1(:,isoil)
      ! potential uptake rate in mol P / m3 / timestep
      plant_uptake_po4_sl(:,isoil) = hlp1(:,isoil) * & 
                                     veg%pool%fine_root%carbon%ptr / soil_depth_sl(:,isoil) * dtime/1.e6
    END DO

    !> 2.0 Loss pathways for the dynamic loss calculation scheme 
    !!
    IF(sb_config%sb_nloss_scheme == 'dynamic') THEN
       !> 2.1 Potential nitrification and denitrification rates (according to Zaehle et al. 2011)
       !! potential nitrification rate from NH4 to NO3 [mol/m3/timestep]
       nitrification_no3(:,:)   = vmax_nitrification * rtm_nitrification_act(:,:) * rmm_nitrification_act(:,:) * & 
                                  (1._wp - anaerobic_volume_fraction_act(:,:)) * nh4_solute(:,:) * dtime/1.e6 
       !! potential denitrification rate from NO3 to N2 [mol/m3/timestep]
       denitrification_n2(:,:)  = vmax_denitrification * rtm_denitrification_act(:,:) * &
                                  sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) / ( km_denitrification_c_act(:,:) + sb%pool%microbial%carbon%ptr(ics:ice,:,iblk))  * &
                                  anaerobic_volume_fraction_act(:,:) * no3_solute(:,:) / ( km_denitrification_no3_act(:,:) + no3_solute(:,:)) * dtime/1.e6 

       !> 2.2. Leaching rate
       !!
       IF(sb_config%sb_model_scheme == 'simple_0d')THEN
          ! For the 0D soil model, loss is proportional to the fractional waterloss of the soil column
          ! reduction of layer concentration in mol m-3 timestep-1
          transport_nh4_solute(:,1) = - fleach_nh4 * nh4_solute(:,1) * drainage_fraction(:) * dtime
          transport_no3_solute(:,1) = - no3_solute(:,1) * drainage_fraction(:) * dtime
          transport_po4_solute(:,1) = - fleach_po4 * po4_solute(:,1) * drainage_fraction(:) * dtime

       ELSE 
          IF(soil_config%soil_moisture_scheme=='layer')THEN
            ! For the 1D scheme with multiple soil hydrology layers, use adjection-diffusion scheme of JSM
            ! reduction of layer concentration in mol m-3 timestep-1, leaching rate in mol m-2 timestep-1
            CALL calc_liquid_phase_transport(options, soil_config%nsoil_water, percolation_sl(:,:), soil_depth_sl(:,:), &
                                             fleach_nh4 * nh4_solute(:,:),transport_nh4_solute(:,:),leaching_nh4_solute(:))
            CALL calc_liquid_phase_transport(options, soil_config%nsoil_water, percolation_sl(:,:), soil_depth_sl(:,:), &
                                             no3_solute(:,:),transport_no3_solute(:,:),leaching_no3_solute(:))
            CALL calc_liquid_phase_transport(options, soil_config%nsoil_water, percolation_sl(:,:), soil_depth_sl(:,:), &
                                             fleach_po4 * po4_solute(:,:),transport_po4_solute(:,:),leaching_po4_solute(:))
          ENDIF

       ENDIF

    ENDIF 

    !> 3.0 Calculate decomposition of pools: dP/dt = - f(T)*f(M)/tau_pool * P
    !!
    !> 3.1 non-woody litter (= polymeric pool)
    !! 
    sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) = sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + &
                                                          rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_non_woody_litter * dtime * &
                                                          sb%pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk)
    sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) = sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + &
                                                            rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_non_woody_litter * dtime * &
                                                            sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)
    sb%loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = sb%loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + &
                                                              rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_non_woody_litter * dtime * &
                                                              sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk)
    sb%loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) = sb%loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + &
                                                            rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_non_woody_litter * dtime * &
                                                            sb%pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk)
    sb%loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = sb%loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + &
                                                            rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_non_woody_litter * dtime * &
                                                            sb%pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk)
    sb%loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb%loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + &
                                                              rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_non_woody_litter * dtime * &
                                                              sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)

    !> 3.2 woody litter
    !! 
    sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk) = sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk) + &
                                                      rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_woody_litter * dtime * &
                                                      sb%pool%woody_litter%carbon%ptr(ics:ice,:,iblk)
    sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) = sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + &
                                                        rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_woody_litter * dtime * &
                                                        sb%pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk)
    sb%loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) = sb%loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) + &
                                                          rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_woody_litter * dtime * &
                                                          sb%pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk)
    sb%loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) = sb%loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) + &
                                                        rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_woody_litter * dtime * &
                                                        sb%pool%woody_litter%carbon13%ptr(ics:ice,:,iblk)
    sb%loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) = sb%loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) + &
                                                        rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_woody_litter * dtime * &
                                                        sb%pool%woody_litter%carbon14%ptr(ics:ice,:,iblk)
    sb%loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb%loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) + &
                                                          rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_woody_litter * dtime * &
                                                          sb%pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk)

    !> 3.3 fast SOM pool (=microbial pool)
    !!
    sb%loss%microbial%carbon%ptr(ics:ice,:,iblk) = sb%loss%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                   rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_fast * dtime * &
                                                   sb%pool%microbial%carbon%ptr(ics:ice,:,iblk)
    sb%loss%microbial%nitrogen%ptr(ics:ice,:,iblk) = sb%loss%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_fast * dtime * &
                                                     sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk)
    sb%loss%microbial%phosphorus%ptr(ics:ice,:,iblk) = sb%loss%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                       rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_fast * dtime * &
                                                       sb%pool%microbial%phosphorus%ptr(ics:ice,:,iblk)
    sb%loss%microbial%carbon13%ptr(ics:ice,:,iblk) = sb%loss%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_fast * dtime * &
                                                     sb%pool%microbial%carbon13%ptr(ics:ice,:,iblk)
    sb%loss%microbial%carbon14%ptr(ics:ice,:,iblk) = sb%loss%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_fast * dtime * &
                                                     sb%pool%microbial%carbon14%ptr(ics:ice,:,iblk)
    sb%loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) = sb%loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                       rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_fast * dtime * &
                                                       sb%pool%microbial%nitrogen15%ptr(ics:ice,:,iblk)

    !> 3.4 slow SOM pool (=residue pool)
    !!
    sb%loss%residue%carbon%ptr(ics:ice,:,iblk) = sb%loss%residue%carbon%ptr(ics:ice,:,iblk) + & 
                                                 rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_slow * dtime * &
                                                 sb%pool%residue%carbon%ptr(ics:ice,:,iblk)
    sb%loss%residue%nitrogen%ptr(ics:ice,:,iblk) = sb%loss%residue%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                   rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_slow * dtime * &
                                                   sb%pool%residue%nitrogen%ptr(ics:ice,:,iblk)
    sb%loss%residue%phosphorus%ptr(ics:ice,:,iblk) = sb%loss%residue%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_slow * dtime * &
                                                     sb%pool%residue%phosphorus%ptr(ics:ice,:,iblk)
    sb%loss%residue%carbon13%ptr(ics:ice,:,iblk) = sb%loss%residue%carbon13%ptr(ics:ice,:,iblk) + & 
                                                   rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_slow * dtime * &
                                                   sb%pool%residue%carbon13%ptr(ics:ice,:,iblk)
    sb%loss%residue%carbon14%ptr(ics:ice,:,iblk) = sb%loss%residue%carbon14%ptr(ics:ice,:,iblk) + & 
                                                   rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_slow * dtime * &
                                                   sb%pool%residue%carbon14%ptr(ics:ice,:,iblk)
    sb%loss%residue%nitrogen15%ptr(ics:ice,:,iblk) = sb%loss%residue%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_slow * dtime * &
                                                     sb%pool%residue%nitrogen15%ptr(ics:ice,:,iblk)

    !> 4.0 formation of SOM resulting from the decomposition dPrec/dt = CUE_Pdon * dPdon/dt
    !! net mineralisation follows the stoichiometric difference of the donor and recipient pool, taking account for CUE
    !! for fast -> slow -> fast this only causes positive mineralisation because of the C:N:P and CUE
    !! for litter -> fast this causes in most circumstances negative mineralisation, which is then taken from the 
    !!               mineral pools (or organic matter for P)

    !> 4.1 formation of fast SOM from litter 
    !! 
    ! calculate C growth of fast pool
    growth_c(:,:) = microbial_cue_max * & 
                    (sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + & 
                     sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk))
    
    ! calculate potential microbial immobilisation rate as difference between maximal growth rate and N in litter 
    mic_uptake_sl(:,:) = growth_c(:,:) / microbial_cn - & 
                         (sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                          sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk))

    ! limit microbial and plant uptake, as well as nitrification to available NH4 minus transport losses
    hlp1(:,:) = mic_uptake_sl(:,:) + plant_uptake_nh4_sl(:,:) + nitrification_no3(:,:) 
    hlp2(:,:) = nh4_solute(:,:) + transport_nh4_solute(:,:)
    WHERE(hlp2(:,:) < hlp1(:,:)) 
        WHERE(hlp1(:,:) > eps8)
           mic_uptake_sl(:,:) = hlp2(:,:) * mic_uptake_sl(:,:) / hlp1(:,:) 
           plant_uptake_nh4_sl(:,:) = hlp2(:,:) * plant_uptake_nh4_sl(:,:) / hlp1(:,:)
           nitrification_no3(:,:) = hlp2(:,:) * nitrification_no3(:,:) / hlp1(:,:)
        ELSEWHERE
           mic_uptake_sl(:,:) = 0.0_wp 
           plant_uptake_nh4_sl(:,:) = 0.0_wp 
           nitrification_no3(:,:) = 0.0_wp 
        ENDWHERE
    ENDWHERE 
   
    ! limit plant uptake and denitrification to available NO3 minus transport losses
    hlp1(:,:) = plant_uptake_no3_sl(:,:) + denitrification_n2(:,:) 
    hlp2(:,:) = no3_solute(:,:) + transport_no3_solute(:,:)
    WHERE(hlp2(:,:) < hlp1(:,:))
        WHERE(hlp1(:,:) > eps8)
           plant_uptake_no3_sl(:,:) = hlp2(:,:) * plant_uptake_no3_sl(:,:) / hlp1(:,:)
           denitrification_n2(:,:) = hlp2(:,:) * denitrification_n2(:,:) / hlp1(:,:)
        ELSEWHERE
           plant_uptake_no3_sl(:,:) = 0.0_wp 
           denitrification_n2(:,:) = 0.0_wp 
        ENDWHERE
    ENDWHERE 
    
    ! calculate C growth based on N availability
    growth_n(:,:) = microbial_cn * (sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                    sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + mic_uptake_sl(:,:))

    ! adjust carbon-use efficiency to match N availability by unchanged litter decomposition
    WHERE(growth_c(:,:) <= growth_n(:,:)) ! enough N is available -> constant CUE
       apparent_mic_cue(:,:) = microbial_cue_max
    ELSEWHERE(growth_c(:,:) > eps8) ! not enough N is available -> changed CUE (within min/max CUE)
       apparent_mic_cue(:,:) = MAX(microbial_cue_min, & 
                                   growth_n(:,:) / ( growth_c(:,:) / microbial_cue_max ))
    ELSEWHERE ! growth C is too small -> CUE does not matter
       apparent_mic_cue(:,:) = microbial_cue_min 
    ENDWHERE

    ! calculate actual litter loss rate, calculated such that litter decomposition
    ! matches the actual immobilisation flux and CN uptake from litter given actual CUE
    hlp1(:,:) = (apparent_mic_cue(:,:) / microbial_cn * &
                      (sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) - &
                      (sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk)))
    hlp2(:,:) = 1.0_wp
    WHERE(mic_uptake_sl(:,:)<hlp1(:,:))
       WHERE(hlp1(:,:)>eps8)
           hlp2(:,:) = mic_uptake_sl(:,:) / hlp1(:,:)
       ELSEWHERE
           hlp2(:,:) = 0.0_wp
       ENDWHERE
       sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb%loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk)
       sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb%loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)
       sb%loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb%loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk)
       sb%loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb%loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk)
       sb%loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb%loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk)
       sb%loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb%loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)

       sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb%loss%woody_litter%carbon%ptr(ics:ice,:,iblk)
       sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb%loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk)
       sb%loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb%loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk)
       sb%loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb%loss%woody_litter%carbon13%ptr(ics:ice,:,iblk)
       sb%loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb%loss%woody_litter%carbon14%ptr(ics:ice,:,iblk)
       sb%loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb%loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk)
    ENDWHERE 

    !> 5.0 formation of SOM and associated respiration, net mineralisation fluxes
    !!
    !> 5.1 formation of fast SOM from litter 
    !! 
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

    !> 5.2 formation of slow SOM (=residue pool) from fast SOM decay
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

    !> 5.3 formation of fast SOM (=microbial pool) from slow SOM decay 
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

    !> 5.4 calculate biomineralisation
    !! calculate potential P immobilisation rate 
    hlp1(:,:) = ABS(net_mineralisation_po4(:,:))
    hlp2(:,:) = hlp1(:,:) + plant_uptake_po4_sl(:,:) 
    hlp3(:,:) = po4_solute + transport_po4_solute(:,:)
    ! limit plant uptake to avaialble P 
    WHERE(net_mineralisation_po4(:,:) < 0.0_wp .AND. & 
          hlp2(:,:) > hlp3(:,:) .AND. hlp2(:,:) > 1.e-13)
        plant_uptake_po4_sl(:,:) = hlp3(:,:) * plant_uptake_po4_sl(:,:) / hlp2(:,:)
    ENDWHERE 
    
    ! limit to available P 
    hlp1(:,:) = ABS(net_mineralisation_po4(:,:)) - (po4_solute(:,:) - plant_uptake_po4_sl(:,:))
    WHERE(net_mineralisation_po4(:,:) < 0.0_wp .AND. hlp1(:,:) > 0.0_wp) 
          net_mineralisation_po4(:,:) = - (po4_solute(:,:) - plant_uptake_po4_sl(:,:))
          sb%loss%residue%phosphorus%ptr(ics:ice,:,iblk) = sb%loss%residue%phosphorus%ptr(ics:ice,:,iblk) + hlp1
    ENDWHERE

    !> 6.0 calculate N losses and nitrification
    !! 
    IF(sb_config%sb_nloss_scheme == 'fixed') THEN
       !> 6.1 case of fixed N losses: N loss proportional to net mineralisation
       WHERE(net_mineralisation_nh4(:,:)>0.0_wp)
          ! N lost to atmosphere during processing of net mineralisation [mol m-3 timestep-1]
          volatilisation_n(:,:)  = net_mineralisation_nh4(:,:) * floss_nmin 
          ! NH4 transformed to NO3 [mol m-3 timestep-1]
          nitrification_no3(:,:) = net_mineralisation_nh4(:,:) * (1._wp - floss_nmin) * fnit_nmin 
       ENDWHERE
       ! report volatilisation losses as N2 emission [mol m-2 timestep-1]
       emission_n2(:) = SUM(volatilisation_n(:,:) * soil_depth_sl(:,:),DIM=2)
    ELSE
       !> 6.2 case of dynamic N loss calculations, following Zaehle et al. 2011
       !!
       ! nitrification rate to produce NO3, NOy and N2O from NH4 [mol m-3 timestep-1]
       nitrification_noy(:,:) = f_nit_noy * nitrification_no3(:,:)
       nitrification_n2o(:,:) = f_nit_n2o * nitrification_no3(:,:)
       nitrification_no3(:,:) = nitrification_no3(:,:) - nitrification_noy(:,:) - nitrification_n2o(:,:)
 
       ! denitrification rate to produce NOy, N2O and N2 from NO3 [mol m-3 timestep-1]
       denitrification_noy(:,:) = f_denit_noy * rtm_denitrification_act(:,:) * denitrification_n2(:,:)
       denitrification_n2o(:,:) = f_denit_n2o * rtm_denitrification_act(:,:) * denitrification_n2(:,:)
       denitrification_n2(:,:) = denitrification_n2(:,:) - denitrification_noy(:,:) - denitrification_n2o(:,:)

       ! gaseous loss of NOy, N2O, and N2, for now following Xu-Ri and Prentice 2008, thus ignoring any actual gas-transport
       ! and instead assume that from each layer a certain fraction excapes directly to the atmosphere.
       ! [mol m-3 timestep-1] 
       transport_noy(:,:) = - rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * noy(:,:) / one_day * dtime
       transport_n2o(:,:) = - rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * n2o(:,:) / one_day * dtime 
       transport_n2(:,:)  = - rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * n2(:,:) / one_day * dtime 

       ! emission rates as positive flux to the atmosphere, integrated over the soil column [mol m-2 timestep-1]
       emission_noy(:) = SUM(-transport_noy(:,:) * soil_depth_sl(:,:),DIM=2)
       emission_n2o(:) = SUM(-transport_n2o(:,:) * soil_depth_sl(:,:),DIM=2)
       emission_n2(:) = SUM(-transport_n2(:,:) * soil_depth_sl(:,:),DIM=2)

       ! record leaching losses for simple_0d model as [mol m-2 timestep-1] (simple_1d has already been calculated in 2.)
       IF(sb_config%sb_model_scheme=='simple_0d')THEN
          leaching_nh4_solute(:) = -transport_nh4_solute(:,1) * soil_depth_sl(:,1)
          leaching_no3_solute(:) = -transport_no3_solute(:,1) * soil_depth_sl(:,1)
          leaching_po4_solute(:) = -transport_po4_solute(:,1) * soil_depth_sl(:,1)
       ENDIF

    ENDIF

    !> 7.0 Update soil litter pools with litter fall
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

    DO isoil = 1,soil_config%nsoil_water

       sb%formation%polymeric_litter%carbon%ptr(ics:ice,isoil,iblk)     = sb%formation%polymeric_litter%carbon%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%carbon%ptr     + veg%litterfall%coarse_root%carbon%ptr )     * & 
          root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb%formation%polymeric_litter%nitrogen%ptr(ics:ice,isoil,iblk)   = sb%formation%polymeric_litter%nitrogen%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%nitrogen%ptr   + veg%litterfall%coarse_root%nitrogen%ptr )   * & 
          root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb%formation%polymeric_litter%phosphorus%ptr(ics:ice,isoil,iblk) = sb%formation%polymeric_litter%phosphorus%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%phosphorus%ptr + veg%litterfall%coarse_root%phosphorus%ptr ) * & 
          root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb%formation%polymeric_litter%carbon13%ptr(ics:ice,isoil,iblk)   = sb%formation%polymeric_litter%carbon13%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%carbon13%ptr   + veg%litterfall%coarse_root%carbon13%ptr )   * & 
          root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb%formation%polymeric_litter%carbon14%ptr(ics:ice,isoil,iblk)   = sb%formation%polymeric_litter%carbon14%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%carbon14%ptr   + veg%litterfall%coarse_root%carbon14%ptr )   * & 
          root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb%formation%polymeric_litter%nitrogen15%ptr(ics:ice,isoil,iblk) = sb%formation%polymeric_litter%nitrogen15%ptr(ics:ice,isoil,iblk) + & 
          ( veg%litterfall%fine_root%nitrogen15%ptr + veg%litterfall%coarse_root%nitrogen15%ptr ) * & 
          root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
    END DO

    !> 8.0 convert rates to µmol/ m3 / s (µmol / m2 / s for surface fluxes)
    !!
    het_respiration(:,:) = het_respiration(:,:) / dtime * 1000000._wp
    het_respiration_c13(:,:) = het_respiration_c13(:,:) / dtime * 1000000._wp
    het_respiration_c14(:,:) = het_respiration_c14(:,:) / dtime * 1000000._wp
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) / dtime * 1000000._wp
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) / dtime * 1000000._wp
    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) / dtime * 1000000._wp
    plant_uptake_nh4_sl(:,:) = plant_uptake_nh4_sl(:,:) / dtime * 1000000._wp 
    plant_uptake_no3_sl(:,:) = plant_uptake_no3_sl(:,:) / dtime * 1000000._wp 
    plant_uptake_po4_sl(:,:) = plant_uptake_po4_sl(:,:) / dtime * 1000000._wp 
    volatilisation_n(:,:) = volatilisation_n(:,:) / dtime * 1000000._wp
    nitrification_no3(:,:) = nitrification_no3(:,:) / dtime * 1000000._wp 
    nitrification_noy(:,:) = nitrification_noy(:,:) / dtime * 1000000._wp 
    nitrification_n2o(:,:) = nitrification_n2o(:,:) / dtime * 1000000._wp 
    denitrification_noy(:,:) = denitrification_noy(:,:) / dtime * 1000000._wp 
    denitrification_n2o(:,:) = denitrification_n2o(:,:) / dtime * 1000000._wp 
    denitrification_n2(:,:) = denitrification_n2(:,:) / dtime * 1000000._wp 
    transport_nh4_solute(:,:) = transport_nh4_solute(:,:) / dtime * 1000000._wp
    transport_no3_solute(:,:) = transport_no3_solute(:,:) / dtime * 1000000._wp
    transport_po4_solute(:,:) = transport_po4_solute(:,:) / dtime * 1000000._wp
    leaching_nh4_solute(:) = leaching_nh4_solute(:) / dtime * 1000000._wp
    leaching_no3_solute(:) = leaching_no3_solute(:) / dtime * 1000000._wp
    leaching_po4_solute(:) = leaching_po4_solute(:) / dtime * 1000000._wp
    transport_noy(:,:) = transport_noy(:,:) / dtime * 1000000._wp
    transport_n2o(:,:) = transport_n2o(:,:) / dtime * 1000000._wp
    transport_n2(:,:) = transport_n2(:,:) / dtime * 1000000._wp
    emission_noy(:) = emission_noy(:) / dtime * 1000000._wp
    emission_n2o(:) = emission_n2o(:) / dtime * 1000000._wp
    emission_n2(:) = emission_n2(:) / dtime * 1000000._wp


  END SUBROUTINE calc_sb_fluxes_simple 

END MODULE mo_sb_simple_model
