!>
!! @brief Contains the routines for the plant growth 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_veg_growth

#include "dsl_jsbach.inc"

  USE mo_kind
  USE mo_exception,             ONLY: message_text, finish
  USE mo_base_constants
  USE mo_veg_constants


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_growth
  PUBLIC :: calc_diameter_from_woody_biomass, calc_height_from_diameter

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_growth'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update plant growth (flux\%growth)
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_growth(tile, options)

    USE mo_jsb_class,               ONLY: Get_model
    USE mo_base_memory,             ONLY: t_slm_tile, t_model, t_jsb_task_options
    USE mo_veg_config_class,        ONLY: veg_config
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_veg_respiration,         ONLY: temperature_response_respiration
    USE mo_veg_canopy,              ONLY: calc_dir_optimal_cn_leaf, calc_marginal_canopy_flux_increment
    USE mo_assimi_config_class,     ONLY: assimi_config

    ! Use of process memories
    dsl4jsb_Use_memory(A2L_)   !  USE mo_a2l_memory_class, ONLY: t_a2l_memory
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,    ONLY: t_pheno_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout) :: tile             !< one tile
    TYPE(t_jsb_task_options), INTENT(in)    :: options          !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp), DIMENSION(options%nc)       :: dheight    !< delta height
    INTEGER                               :: ichunk     !< loop over chunks ics:ice / 1:nc
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_growth'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(A2L_)                   !     TYPE(t_a2l_memory), POINTER :: a2l__mem
    dsl4jsb_Def_memory(ASSIMI_)                     !     TYPE(t_assimi_memory),   POINTER :: assimi__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(PHENO_)                      !     TYPE(t_pheno_memory),    POINTER :: pheno__mem
    dsl4jsb_Def_memory(RAD_)                        !     TYPE(t_rad_memory),      POINTER :: rad__mem
    dsl4jsb_Def_memory(SOIL_)                       !     TYPE(t_soil_memory),     POINTER :: soil__mem
    ! Declare pointers to variables in memory
    dsl4jsb_Logical2D_onChunk   :: growing_season             ! PHENO_
    dsl4jsb_Real2D_onChunk      :: t_air                      ! A2L_
    dsl4jsb_Real2D_onChunk      :: q_air
    dsl4jsb_Real2D_onChunk      :: press_srf
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio
    dsl4jsb_Real2D_onChunk      :: ga
    dsl4jsb_Real2D_onChunk      :: swpar_srf_down
    dsl4jsb_Real2D_onChunk      :: swpar_fdiffuse
    dsl4jsb_Real2D_onChunk      :: beta_air                   ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs
    dsl4jsb_Real2D_onChunk      :: maint_respiration_leaf
    dsl4jsb_Real2D_onChunk      :: net_assimilation
    dsl4jsb_Real2D_onChunk      :: canopy_cond
    dsl4jsb_Real2D_onChunk      :: beta_air_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: w_soil_root_theta          ! SOIL_
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_cl             ! RAD_
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tcnl_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tcnl_mavg_cl
    dsl4jsb_Real2D_onChunk      :: target_cn_leaf             ! VEG_
    dsl4jsb_Real2D_onChunk      :: lai
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: press_srf_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: ga_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tacclim_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: leaf_cn_direction
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps
    dsl4jsb_Real2D_onChunk      :: unit_npp
    dsl4jsb_Real2D_onChunk      :: unit_transpiration
    dsl4jsb_Real2D_onChunk      :: kstar_labile
    dsl4jsb_Real2D_onChunk      :: fmaint_rate
    dsl4jsb_Real2D_onChunk      :: fmaint_rate_troot_mavg
    dsl4jsb_Real2D_onChunk      :: unit_npp_troot_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_pot_troot_mavg
    dsl4jsb_Real2D_onChunk      :: cost_n_uptake_root
    dsl4jsb_Real2D_onChunk      :: cost_n_fixation
    dsl4jsb_Real2D_onChunk      :: target_np_leaf
    dsl4jsb_Real2D_onChunk      :: target_cn_fine_root
    dsl4jsb_Real2D_onChunk      :: target_cn_coarse_root
    dsl4jsb_Real2D_onChunk      :: target_cn_sap_wood
    dsl4jsb_Real2D_onChunk      :: target_cn_fruit
    dsl4jsb_Real2D_onChunk      :: target_np_fine_root
    dsl4jsb_Real2D_onChunk      :: target_np_coarse_root
    dsl4jsb_Real2D_onChunk      :: target_np_sap_wood
    dsl4jsb_Real2D_onChunk      :: target_np_fruit
    dsl4jsb_Real2D_onChunk      :: height
    dsl4jsb_Real2D_onChunk      :: dens_ind
    dsl4jsb_Real2D_onChunk      :: uptake_n_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cn_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: npp_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: labile_carbon_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_carbon_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cn_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cp_talloc_mavg 
    dsl4jsb_Real2D_onChunk      :: npp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_npp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: transpiration_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_transpiration_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: n_fixation_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_p_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: dphi_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: leaf2sapwood_mass_ratio 
    dsl4jsb_Real2D_onChunk      :: leaf2root_mass_ratio
    dsl4jsb_Real2D_onChunk      :: gpp_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: maint_respiration_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_p_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: leaf2root_troot_mavg
    dsl4jsb_Real2D_onChunk      :: target_lai
    dsl4jsb_Real2D_onChunk      :: k1_opt
    dsl4jsb_Real2D_onChunk      :: k2_opt
    dsl4jsb_Real2D_onChunk      :: growth_req_p
    dsl4jsb_Real2D_onChunk      :: growth_req_n
    dsl4jsb_Real2D_onChunk      :: uptake_nh4
    dsl4jsb_Real2D_onChunk      :: uptake_no3
    dsl4jsb_Real2D_onChunk      :: n_fixation
    dsl4jsb_Real2D_onChunk      :: n_transform_respiration
    dsl4jsb_Real2D_onChunk      :: n_fixation_respiration
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c13
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c14
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c13
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c14
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration_c13
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration_c14
    dsl4jsb_Real2D_onChunk      :: maint_respiration_pot
    dsl4jsb_Real2D_onChunk      :: maint_respiration
    dsl4jsb_Real2D_onChunk      :: growth_respiration
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tcnl_mavg_cl
!! *slmtemp    
    TYPE(t_veg_memory),      POINTER :: veg         ! local pointer to veg__mem for SLM usability
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(VEG_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(xx_)
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(PHENO_)          ! veg__mem => veg_memory
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SOIL_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(PHENO_, growing_season)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)
    dsl4jsb_Get_var2D_onChunk(A2L_, co2_mixing_ratio)
    dsl4jsb_Get_var2D_onChunk(A2L_, ga)
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)
     dsl4jsb_Get_var2D_onChunk(A2L_, swpar_fdiffuse)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, maint_respiration_leaf)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, net_assimilation)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, canopy_cond)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tcnl_mavg)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_theta)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tcnl_mavg_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tcnl_mavg_cl)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_leaf)
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tacclim_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, leaf_cn_direction)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_npp)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_transpiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, kstar_labile)
    dsl4jsb_Get_var2D_onChunk(VEG_, fmaint_rate)
    dsl4jsb_Get_var2D_onChunk(VEG_, fmaint_rate_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_npp_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_pot_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, cost_n_uptake_root)
    dsl4jsb_Get_var2D_onChunk(VEG_, cost_n_fixation)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_leaf)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_fine_root)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_coarse_root)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_sap_wood)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_fruit)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_fine_root)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_coarse_root)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_sap_wood)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_fruit)
    dsl4jsb_Get_var2D_onChunk(VEG_, height)
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_n_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cn_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_carbon_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_nitrogen_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_phosphorus_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_carbon_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_nitrogen_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_phosphorus_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cn_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_npp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, dphi_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, leaf2sapwood_mass_ratio)
    dsl4jsb_Get_var2D_onChunk(VEG_, leaf2root_mass_ratio)
    dsl4jsb_Get_var2D_onChunk(VEG_, gpp_tlabile_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_tlabile_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tlabile_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tlabile_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, leaf2root_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_lai)
    dsl4jsb_Get_var2D_onChunk(VEG_, k1_opt)
    dsl4jsb_Get_var2D_onChunk(VEG_, k2_opt)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_transform_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c13)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c14)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c13)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c14)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration_c13)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration_c14)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_pot)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_et_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_rub_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_pepc_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_oth_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tcnl_mavg_cl)



!! *slmtemp
    veg     => veg__mem 


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------



    !> 1.0 calculate potential growth rate relative to labile pool size
    !! 
    !! 
    kstar_labile(:) = calc_meristem_activity(t_air(:), w_soil_root_theta(:), growing_season(:))

    !> 1.1 calculate costs of N uptake and fixation
    !!
    veg%fmaint_rate%ptr = fmaint_rate_base * temperature_response_respiration(t_air(1), veg%t_air_tacclim_mavg%ptr)
    CALL calc_nitrogen_uptake_costs(lctlib                               , &
                                    t_air(:)                             , &
                                    veg%pool%fine_root%carbon%ptr        , &
                                    veg%pool%fine_root%nitrogen%ptr      , &
!                                    veg%pool%fine_root%carbon(:)      , &
!                                    veg%pool%fine_root%nitrogen(:)    , &
                                    fmaint_rate_troot_mavg(:)        , &
                                    unit_npp_troot_mavg(:)           , &
                                    unit_uptake_n_pot_troot_mavg(:)  , &
                                    cost_n_uptake_root(:)            , &
                                    cost_n_fixation(:)               , &
                                    uptake_no3(:)                    , &
                                    uptake_nh4(:))           
  
    !> 2.0 Calculate target C:N:P for each tissue type
    !! 
    !! runs only IF (growing_season(ichunk) .AND. veg_pool_leaf_nitrogen > eps4 ) \n
    !! loops internally over: ichunk=1,nc
    !! 
    IF (veg_config%flag_optimal_Nleaf) THEN
      CALL calc_dir_optimal_cn_leaf( &
                                        nc                               , & ! in
                                        assimi_config%ncanopy            , &
                                        lctlib                           , &
                                        growing_season(:)                , &
                                        veg%pool%leaf%carbon%ptr         , &  
                                        veg%pool%leaf%nitrogen%ptr       , & 
!                                        veg%pool%leaf%carbon(ics:ice)    , & 
!                                        veg%pool%leaf%nitrogen(ics:ice)  , & 
                                        lai(:)                           , &
                                        t_jmax_opt_mavg(:)               , &
                                        target_cn_leaf(:)                , &
                                        t_air_tacclim_mavg(:)            , &
                                        t_air_tcnl_mavg(:)               , &
                                        ga_tcnl_mavg(:)                  , &
                                        press_srf_tcnl_mavg(:)           , &
                                        co2_mixing_ratio_tcnl_mavg(:)    , &
                                        beta_sinklim_ps_tcnl_mavg(:)     , &
                                        beta_air_tcnl_mavg(:)            , &
                                        beta_soil_ps_tcnl_mavg(:)        , &
                                        beta_soil_gs_tcnl_mavg(:)        , &
                                        fn_chl_cl(:,:)                   , &
                                        fn_et_cl(:,:)                    , &
                                        fn_rub_cl(:,:)                   , &
                                        fn_pepc_cl(:,:)                  , &
                                        fn_oth_cl(:,:)                   , &
                                        fleaf_sunlit_tcnl_mavg_cl(:,:)   , &
                                        ppfd_sunlit_tcnl_mavg_cl(:,:)    , &
                                        ppfd_shaded_tcnl_mavg_cl(:,:)    , &
                                        uptake_n_tcnl_mavg               , &
                                        growth_cn_tcnl_mavg                 , &
                                        npp_tcnl_mavg                    , &
                                        leaf_cn_direction(:)             )    ! out
    END IF ! veg_config%flag_optimal_Nleaf

    CALL calc_stoichiometry_changes(lctlib                        , & ! in
                                      growing_season(:)           , &
                                      leaf_cn_direction(:)        , &
                                      target_np_leaf(:)           , &
                                      target_cn_leaf(:)           , & ! inout
                                      target_cn_fine_root(:)      , & ! out
                                      target_cn_coarse_root(:)    , &
                                      target_cn_sap_wood(:)       , &
                                      target_cn_fruit(:)          , &
                                      target_np_fine_root(:)      , &
                                      target_np_coarse_root(:)    , &
                                      target_np_sap_wood(:)       , &
                                      target_np_fruit(:)          )

    !> 3.0 calculate current allometry 
    !! 
    !! 
    !> 3.1 calculate marginal NPP increment per unit C growth, needed for optimal allocation scheme
    !!
    CALL calc_marginal_canopy_flux_increment( &
                                              options%nc                  , & ! in
                                              assimi_config%ncanopy       , &
                                              lctlib                      , &
                                              t_air(:)                    , &
                                              press_srf(:)                , &
                                              q_air(:)                    , &
                                              co2_mixing_ratio(:)         , &
                                              ga(:)                       , & 
                                              net_assimilation(:)         , & 
                                              canopy_cond(:)              , & 
                                              beta_air(:)                 , & 
                                              beta_soil_ps(:)             , & 
                                              beta_soil_gs(:)             , & 
                                              veg%pool%leaf%carbon%ptr    , &  
                                              veg%pool%leaf%nitrogen%ptr  , &
!                                              veg%pool%leaf%carbon(:)     , &  
!                                              veg%pool%leaf%nitrogen(:)   , &  
                                              lai(:)                      , &
                                              beta_sinklim_ps(:)          , &
                                              t_air_tacclim_mavg(:)       , &
                                              t_jmax_opt_mavg(:)          , &
                                              fn_chl_cl(:,:)              , &
                                              fn_et_cl(:,:)               , &
                                              fn_rub_cl(:,:)              , &
                                              fn_pepc_cl(:,:)             , &
                                              fn_oth_cl(:,:)              , &
                                              fleaf_sunlit_cl(:,:)        , &
                                              ppfd_sunlit_cl(:,:)         , &
                                              ppfd_shaded_cl(:,:)         , &
                                              swpar_srf_down              , &
                                              swpar_fdiffuse              , &
                                              unit_npp(:)                 , & ! out
                                              unit_transpiration(:)       )   
    
    !> 3.2 calculate allometry    
    !!
    CALL calc_allometry(lctlib                           , &
                        height(:)                        , &
                        veg%pool%sap_wood%carbon%ptr     , &
                        veg%pool%heart_wood%carbon%ptr   , &
!                        veg%pool%sap_wood%carbon(:)     , &
!                        veg%pool%heart_wood%carbon(:)   , &
                        dens_ind(:)                      , & 
                        labile_carbon_talloc_mavg(:)     , &
                        labile_nitrogen_talloc_mavg(:)   , & 
                        labile_phosphorus_talloc_mavg(:) , &
                        reserve_carbon_talloc_mavg(:)    , &
                        reserve_nitrogen_talloc_mavg(:)  , &
                        reserve_phosphorus_talloc_mavg(:), &
                        growth_cn_talloc_mavg(:)         , &
                        growth_cp_talloc_mavg(:)         , & 
                        npp_talloc_mavg(:)               , &
                        unit_npp_talloc_mavg(:)          , &
                        transpiration_talloc_mavg(:)     , &
                        unit_transpiration_talloc_mavg(:), &
                        n_fixation_talloc_mavg(:)        , &
                        unit_uptake_n_talloc_mavg(:) , &
                        unit_uptake_p_talloc_mavg(:) , &
                        dphi_talloc_mavg(:)              , &
                        kstar_labile(:)                  , &
                        dheight(:)                       , &
                        leaf2sapwood_mass_ratio(:)       , & 
                        leaf2root_mass_ratio(:)          , &
                        target_lai(:)                    , &
                        k1_opt(:)                        , &
                        k2_opt(:)) 

    !> 4.0 calculate net flux between reserve and labile pool 
    !! 
    !! Info: veg%reserve_use is of TYPE(t_slm_struct_veg_element)
    call calc_reserve_useage(lctlib                               , &
                             growing_season(1)                    , &
                             kstar_labile(1)                      , &
                             gpp_tlabile_mavg(1)                  , &
                             maint_respiration_tlabile_mavg(1)    , &
                             growth_req_n_tlabile_mavg(1)         , &
                             growth_req_p_tlabile_mavg(1)         , &
                             target_cn_leaf(1)                    , &
                             target_cn_fine_root(1)               , & 
                             target_np_leaf(1)                    , &
                             target_np_fine_root(1)               , &
                             leaf2root_troot_mavg(1)              , &
                             lai(1)                               , &
                             target_lai(1)                        , &
                             veg%pool%labile                      , &
                             veg%pool%reserve                     , &
                             veg%reserve_use                      , &
                             beta_sinklim_ps(1))


    !> 5.0 calculate allocation coefficients for structural plant pools (leaves, wood, fine and coarse root, fruits)
    !! 
    !! 
    CALL calc_allocation_fraction(lctlib                          , & ! in
                                  dheight(1)                      , &
                                  kstar_labile(1)                 , &
                                  veg%pool%leaf%carbon%ptr        , &
                                  veg%pool%sap_wood%carbon%ptr    , &
                                  veg%pool%coarse_root%carbon%ptr , &
                                  veg%pool%fine_root%carbon%ptr   , &
                                  veg%pool%labile%carbon%ptr      , &
                                  height(1)                       , &
                                  leaf2sapwood_mass_ratio(1)      , &
                                  leaf2root_mass_ratio(1)         , & 
                                  k1_opt(1)                       , &
                                  k2_opt(1)                       , &
                                  target_cn_leaf(1)               , &
                                  target_cn_fine_root(1)          , &
                                  target_cn_coarse_root(1)        , &
                                  target_cn_sap_wood(1)           , &
                                  target_cn_fruit(1)              , &
                                  target_np_leaf(1)               , &
                                  target_np_fine_root(1)          , &
                                  target_np_coarse_root(1)        , &
                                  target_np_sap_wood(1)           , &
                                  target_np_fruit(1)              , &
                                  veg%reserve_use%carbon%ptr      , & 
                                  veg%falloc                      , & ! inout  TYPE(t_slm_struct_veg_pool)
                                  growth_req_n(1)                 , & ! out
                                  growth_req_p(1)                 )

    !> 6.0 partition labile carbon into growth, growth_respiration, and maintenance 
    !! 
    !! 
    CALL calc_partitioning(t_air(1)                               , & ! in
                           t_air_tacclim_mavg(1)            , &
                           kstar_labile(1)                  , &
                           maint_respiration_leaf(1)        , &
                           cost_n_fixation(1)               , &
                           growth_req_p(1)                  , &
                           growth_req_n(1)                  , &
                           veg%falloc                       , &
                           veg%pool                         , &
                           uptake_nh4(1)                    , & ! inout
                           uptake_no3(1)                    , & ! inout
                           n_fixation(1)                    , &
                           n_transform_respiration(1)       , &
                           n_fixation_respiration(1)        , &
                           n_processing_respiration(1)      , &
                           maint_respiration_c13(1)         , &
                           maint_respiration_c14(1)         , &
                           growth_respiration_c13(1)        , &
                           growth_respiration_c14(1)        , &
                           n_processing_respiration_c13(1)  , &
                           n_processing_respiration_c14(1)  , &
                           veg%growth                       , &
                           maint_respiration_pot(1)         , & ! out
                           maint_respiration(1)             , &
                           growth_respiration(1)  )  

  END SUBROUTINE update_growth 


  !-----------------------------------------------------------------------------------------------------
  ! 1.0 Sub Task of update_growth
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculates the activity of the meristem for growth 
  !!
  !! As this is proportional to the size of the labile pool, there is no need for a direct temperature 
  !! function of the turnover, as this will result from the dynamics of the labile pool. \nThis function
  !!  therefore only takes account of the impact of low temperatures (<5°C) and low moisture 
  !! (<20\% saturation). \nAdditionally, phenological dormancy is accounted for.
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_meristem_activity(t_air,w_soil_root_theta,growing_season) RESULT(kstar_labile)

    USE mo_shared_constants,  ONLY: Tzero
    USE mo_veg_constants,     ONLY: lambda_labile_temp  , &
                                    lambda_labile_theta , &
                                    k_labile_temp       , & 
                                    k_labile_theta      , &
                                    tau_labile
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),   INTENT(in)    :: t_air                  !< air temperature (K) 
    REAL(wp),   INTENT(in)    :: w_soil_root_theta      !< soil moisture scalar (fraction of available water holding capacity) 
    LOGICAL,    INTENT(in)    :: growing_season         !< is the plant growing? 
    REAL(wp)                  :: kstar_labile           !< current potential labile pool turnover (1/day)
    ! 0.2 Local
    REAL(wp)                :: ftemp                  ! temperature limit on cell construction rate
    REAL(wp)                :: fmoist                 ! moisture limit on cell construction rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_meristem_activity'
  
    IF(t_air > Tzero) THEN 
        ftemp = 1.0_wp - exp(-(lambda_labile_temp * (t_air - Tzero))**k_labile_temp)
      ELSE 
        ftemp = 0.0_wp
    ENDIF
    fmoist = 1.0_wp - exp(-(lambda_labile_theta * w_soil_root_theta)**k_labile_theta) 
    IF(growing_season)THEN
        kstar_labile = ftemp * fmoist / tau_labile
      ELSE 
        kstar_labile = 0.0_wp
    ENDIF

  END FUNCTION calc_meristem_activity 


  !-----------------------------------------------------------------------------------------------------
  ! 1.1 Sub Task of update_growth
  ! 
  !----------------------------------------------------------------------------------------------------------
  !> Routine to calculate the C costs for (inorganic) N uptake and BNF in order to optimalize N acquisition
  !!   
  !! Output: C costs of root uptake and N fixation
  !! 
  !! @todo: reconsider the time scale of these arguments. Since cost_n_fixation is used for respiration
  !! calculation, it should be based on actual temperature, not its average
  !----------------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_nitrogen_uptake_costs(lctlib                    , &
                                                  t_air                     , &
                                                  fine_root_carbon          , &
                                                  fine_root_nitrogen        , &
                                                  fmaint_rate_troot_mavg    , &
                                                  unit_npp_troot_mavg       , &
                                                  unit_uptake_n_troot_mavg  , &
                                                  cost_n_uptake_root        , &
                                                  cost_n_fixation           , &
                                                  no3                       , &
                                                  nh4)         

    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element    
    USE mo_shared_constants,        ONLY: Tzero
    USE mo_veg_constants,           ONLY: transform_cost_nh4,  transform_cost_no3, fresp_growth, & 
                                          cost_n_fixation_max, cost_n_fixation_min, & 
                                          k1_temp_bnf, k2_temp_bnf,  k3_temp_bnf,  k4_temp_bnf
    
    IMPLICIT NONE
    
    ! 0.1 InOut
    TYPE(t_lctlib_element),   INTENT(in)  :: lctlib                   !< land-cover-type library - parameter across pft's
    REAL(wp),                 INTENT(in)  :: t_air                    !< air temperature (K) !! should be soil temperature!
    REAL(wp),                 INTENT(in)  :: fine_root_carbon         !< fine root carbon mass (mol C / m2)
    REAL(wp),                 INTENT(in)  :: fine_root_nitrogen       !< fine root nitrogen mass (mol N / m2) 
    REAL(wp),                 INTENT(in)  :: fmaint_rate_troot_mavg   !< N-speficic root maintentance respiration rate (µmol C/mol N / s)
    REAL(wp),                 INTENT(in)  :: unit_npp_troot_mavg      !< increase of NPP with increased leaf C investment (??)
    REAL(wp),                 INTENT(in)  :: unit_uptake_n_troot_mavg !< increase of root N uptake with increased root C inventment (??)
    REAL(wp),                 INTENT(in)  :: no3                      !< previous uptake of NO3
    REAL(wp),                 INTENT(in)  :: nh4                      !< previous uptake of NH4
    REAL(wp),                 INTENT(out) :: cost_n_uptake_root       !< cost of root N uptake (mol C / mol N)
    REAL(wp),                 INTENT(out) :: cost_n_fixation          !< cost of N fixation (mol C / mol N)
    ! 0.2 Local 
    REAL(wp) :: fresp_maint
    REAL(wp) :: gc
    REAL(wp) :: gn
    REAL(wp) :: tf                                    !! BNF-temperature-function after Houlton (2008)
    REAL(wp) :: f_nacq_nh4
    REAL(wp) :: f_nacq_no3
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_nitrogen_uptake_costs'


    !> 1.0 C costs for N acquisition
    !!

    ! maintenance respiration cost of a unit fine root carbon
    IF(fine_root_carbon > eps8)THEN
        fresp_maint = fmaint_rate_troot_mavg * fine_root_nitrogen / fine_root_carbon  * &
                      lctlib%tau_fine_root * one_day * one_year / 1.e6_wp
    ELSE
        fresp_maint = 0.0_wp
    ENDIF

    !> 1.1 C costs for (inorganic) N uptake
    !! based on Rastetter (2001) and Meyerholt (2016, OPT)

    ! photosynthesis gain relative to carbon investment into leafs
    gc = MAX(unit_npp_troot_mavg,0.0_wp)
        
    ! n-uptake gain relative to carbon investment into roots
    ! multiplier is Croot / Croot * (1 + fresp_growth + fresp_maint
    gn = unit_uptake_n_troot_mavg /(1._wp+fresp_growth+fresp_maint)
    
    f_nacq_nh4 = nh4 / (nh4 + no3)
    f_nacq_no3 = no3 / (nh4 + no3)
    IF(gn > eps8)THEN
        cost_n_uptake_root = gc / gn + &
                             (f_nacq_nh4 * transform_cost_nh4 + f_nacq_no3 * transform_cost_no3)
    ELSE
        cost_n_uptake_root = 1000._wp ! an arbitrarily large value !
    ENDIF
                       
    !> 1.2 C costs for BNF
    !! 

    ! Temperature function based on Houlton (2008)
    tf = k1_temp_bnf * exp(k2_temp_bnf + k3_temp_bnf * (t_air - Tzero) * &
                           (1.0_wp-(t_air-Tzero)/k4_temp_bnf))
       
    cost_n_fixation = (cost_n_fixation_max-(cost_n_fixation_max-cost_n_fixation_min)*tf) 

  END SUBROUTINE calc_nitrogen_uptake_costs 


  !-----------------------------------------------------------------------------------------------------
  ! 2.0 Sub Task of update_growth
  ! 
  ! ------------------------------------------------------------------------------------------
  !> Routine to calculate the target stoichiometry of the newly growth tissue
  !! 
  !! Includes optimal CN ratio calculation \n
  !! Output: target stoichiometries
  ! --------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_stoichiometry_changes(lctlib                  , &   ! in
                                                  growing_season          , &
                                                  leaf_cn_direction       , &
                                                  target_np_leaf          , &
                                                  target_cn_leaf          , &   ! inout
                                                  target_cn_fine_root     , &   ! out
                                                  target_cn_coarse_root   , &
                                                  target_cn_sap_wood      , &
                                                  target_cn_fruit         , &
                                                  target_np_fine_root     , &
                                                  target_np_coarse_root   , &
                                                  target_np_sap_wood      , &
                                                  target_np_fruit  )         

    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_veg_config_class,      ONLY: veg_config
    USE mo_base_constants,        ONLY: dtime, one_day
    USE mo_veg_constants,         ONLY: delta_n_leaf

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),   INTENT(in)    :: lctlib                    !< land-cover-type library - parameter across pft's
    LOGICAL,                  INTENT(in)    :: growing_season            !< 
    REAL(wp),                 INTENT(in)    :: leaf_cn_direction     , & !< 
                                               target_np_leaf            !< 
    REAL(wp),                 INTENT(inout) :: target_cn_leaf            !< 
    REAL(wp),                 INTENT(out)   :: target_cn_fine_root   , & !< fine root target C:N
                                               target_cn_coarse_root , & !< coarse root target C:N
                                               target_cn_sap_wood    , & !< sap wood target C:N
                                               target_cn_fruit       , & !< fruit target C:N
                                               target_np_fine_root   , & !< fine root target N:P
                                               target_np_coarse_root , & !< coarse root target N:P
                                               target_np_sap_wood    , & !< sap wood target N:P
                                               target_np_fruit           !< fruit target N:P
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_stoichiometry_changes'


    !> 1.0 calculate optimal allocation fractions
    !! 
    IF (veg_config%flag_optimal_Nleaf .AND. growing_season) THEN
       ! update leaf CN target. leaf_cn_direction has values of -1, 0 or 1
       target_cn_leaf = target_cn_leaf * (1.0_wp - leaf_cn_direction * delta_n_leaf * dtime/one_day)
    ENDIF

    !> 2.0 determine C:N of other plant tissues
    !!
    target_cn_fine_root = lctlib%cn_leaf / root2leaf_cn
    SELECT CASE(lctlib%growthform)
      CASE (itree) ! trees
         target_cn_coarse_root  = lctlib%cn_leaf / wood2leaf_cn
         target_cn_sap_wood     = lctlib%cn_leaf / wood2leaf_cn 
      CASE (igrass) ! grasses 
         target_cn_coarse_root  = lctlib%cn_leaf / root2leaf_cn
         target_cn_sap_wood     = lctlib%cn_leaf / root2leaf_cn 
    END SELECT
    target_cn_fruit = lctlib%cn_leaf / root2leaf_cn

    !> 3.0 determine N:P of other plant tissues
    !!
    target_np_fine_root = target_np_leaf / root2leaf_np
    SELECT CASE(lctlib%growthform)
      CASE (itree) ! trees
         target_np_coarse_root  = target_np_leaf / wood2leaf_np
         target_np_sap_wood     = target_np_leaf / wood2leaf_np
      CASE (igrass) ! grasses 
         target_np_coarse_root  = target_np_leaf / root2leaf_np
         target_np_sap_wood     = target_np_leaf / root2leaf_np
    END SELECT
    target_np_fruit = target_np_leaf / root2leaf_np

  END SUBROUTINE calc_stoichiometry_changes


  !-----------------------------------------------------------------------------------------------------
  ! 3.2 Sub Task of update_growth  (3.1 calc_marginal_canopy_flux_increment:mo_veg_canopy)
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculates the current allometry of a tree 
  !!
  !! diagnoses diameter and height from woody biomass by assuming that the trunk is a cylinder following
  !!  a given allometric relationship: H = k_allom1 * D ** k_allom2
  !!
  !! calculates leaf to sapwood and root mass ratios
  !!
  !! calculates maximum LAI given current root and sapwood mass
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE calc_allometry(lctlib                         , &
                                            height                        , &
                                            veg_pool_sap_wood_carbon      , &
                                            veg_pool_heart_wood_carbon    , &
                                            dens_ind                      , &
                                            labile_carbon                 , &
                                            labile_nitrogen               , &
                                            labile_phosphorus             , &
                                            reserve_carbon            , &
                                            reserve_nitrogen          , &
                                            reserve_phosphorus        , &
                                            growth_cn                     , &
                                            growth_cp                     , &
                                            npp                           , &
                                            unit_npp                      , &
                                            transpiration                 , &
                                            unit_transpiration            , &
                                            n_fixation                    , &
                                            unit_uptake_n                 , &
                                            unit_uptake_p                 , &
                                            dphi                          , &
                                            kstar_labile                  , &
                                            dheight                       , &
                                            leaf2sapwood_mass_ratio       , &
                                            leaf2root_mass_ratio          , &
                                            target_lai                    , &
                                            k1_opt                        , &
                                            k2_opt )                    

    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_veg_config_class,        ONLY: veg_config    
    USE mo_base_constants,          ONLY: pi
    USE mo_shared_constants,        ONLY: rhoh2o, mavg_period_tgrowth
    USE mo_veg_constants,           ONLY: min_height, itree, igrass, sm2lm_grass, target_lai_max

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),INTENT(in) :: lctlib           !< land-cover-type library - parameters for one pft
    REAL(wp), INTENT(in)  :: height                       !< height (m)
    REAL(wp), INTENT(in)  :: veg_pool_sap_wood_carbon     !< woody biomass of an individual tree (mol C/m2)
    REAL(wp), INTENT(in)  :: veg_pool_heart_wood_carbon   !< woody biomass of an individual tree (mol C/m2)
    REAL(wp), INTENT(in)  :: dens_ind                     !< density of populations (number of indivuduals /m2) 
    REAl(wp), INTENT(in)  :: labile_carbon                !< labile carbon (mol C/m2) 
    REAL(wp), INTENT(in)  :: labile_nitrogen              !< labile nitrogen (mol N/m2)
    REAL(wp), INTENT(in)  :: labile_phosphorus            !< labile phosphorus (mol P/m2)
    REAL(wp), INTENT(in)  :: reserve_carbon               !< reserve carbon (mol C)
    REAL(wp), INTENT(in)  :: reserve_nitrogen             !< reserve nitrogen (mol N)
    REAL(wp), INTENT(in)  :: reserve_phosphorus           !< reserve phosphorus (mol P)
    REAL(wp), INTENT(in)  :: growth_cn                    !< C:N of whole plant growth increment
    REAL(wp), INTENT(in)  :: growth_cp                    !< C:P of whole plant growth increment
    REAL(wp), INTENT(in)  :: npp                          !< current average NPP (µmol C m-2 s-1)
    REAL(wp), INTENT(in)  :: unit_npp                     !< increase in NPP for unit increase in leaf mass (µmol C m-2 s-1 mol-1 leaf C)
    REAL(wp), INTENT(in)  :: transpiration                !< average transpiration flux (kg m-2 s-1)
    REAL(wp), INTENT(in)  :: unit_transpiration           !< increase in transpiration with unit increase in leaf mass (kg m-1 s-2 mol-1 leaf C)
    REAL(wp), INTENT(in)  :: n_fixation                   !< N fixation (µmol N m-2 s-1)
    REAL(wp), INTENT(in)  :: unit_uptake_n                !< N uptake per unit root mass (µmol N m-2 s-1 mol-1 root C)
    REAL(wp), INTENT(in)  :: unit_uptake_p                !< P uptake per unit root mass (µmol P m-2 s-1 mol-1 root C)
    REAL(wp), INTENT(in)  :: dphi                         !< leaf-soil water potential gradient (MPa) 
    REAL(wp), INTENT(in)  :: kstar_labile                 !< labile pool turnover rate
    REAL(wp), INTENT(out) :: dheight                      !< change in height (m) resulting from a small change in woody biomass
    REAL(wp), INTENT(out) :: leaf2sapwood_mass_ratio      !< leaf to sapwood mass ratio 
    REAL(wp), INTENT(out) :: leaf2root_mass_ratio         !< leaf to root mass ratio
    REAL(wp), INTENT(out) :: target_lai                   !< lai implied by current available supporting structure 
    REAL(wp), INTENT(out) :: k1_opt                       !< coefficient in optimized leaf:root ratio
    REAL(wp), INTENT(out) :: k2_opt                       !< coefficient in optimized leaf:root ratio
    ! 0.2 Local
    REAL(wp)             :: diameter2                     ! diameter (m) resulting from a infitisemal small increase in woody biomas 
    REAL(wp)             :: k1_opt_w  , &
                            k2_opt_w
    REAL(wp)             :: k1_opt_n      , &
                            k1_opt_p
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_allometry'


    !> 1.0 case growthform
    !! 
    SELECT CASE (lctlib%growthform)
      !> 1.1 trees
      !! 
      CASE (itree)
        ! Diameter and Height calculated assuming that 
        !   all woody biomass forms a cylinder: diameter = SQRT(4 * wood volume / height / pi)
        !   height = allom_k1 * diameter ^ allom_k2
        IF(dens_ind > eps8) THEN
           ! Calculate change of height given a small change in woody biomass
           diameter2 = calc_diameter_from_woody_biomass(lctlib,veg_pool_sap_wood_carbon+eps4,veg_pool_heart_wood_carbon,dens_ind)
           dheight   = calc_height_from_diameter(lctlib,diameter2) - height
        ELSE 
           dheight   = 0.0_wp
        ENDIF
 
        ! scaling factor for leaf to wood mass, derived from pipe-model (Shinozaki et al. 1964, J J Ecol)
        ! leaf area = k_latosa * sap_wood area
        ! => Cl * sla = k * Cs / wood_density / height
        ! => Cl = k_latosa / sla / wood_density / height * Cs = leaf2sapwood_mass_ratio * Cs / height
        leaf2sapwood_mass_ratio = lctlib%k_latosa / lctlib%sla / lctlib%wood_density

        ! full leaf area implied by current sapwood area
        IF (height > eps8) THEN
          target_lai = veg_pool_sap_wood_carbon * lctlib%k_latosa / lctlib%wood_density / height
        ELSE
          target_lai = 0.0_wp
        END IF

        IF (veg_config%flag_optimal_allocation) THEN
           !> Calculate allocation fractions assuming optimal relationship between leaf and root allocation
           !! potential carbon growth rate = potential nittrogen growth rate * plant cn ratio
           !! NPP = gN * R * growth_cn
           !! Or similar relationships for P
           !! Difference between BNF schemes - for optimal scheme fixation depends on root mass
           SELECT CASE (veg_config%n_uptake_scheme)
           CASE ("prescribed", "standard","resistance")
              k1_opt_n = (npp*mavg_period_tgrowth*one_day +&
                   (labile_carbon + reserve_carbon - &
                   (labile_nitrogen + reserve_nitrogen) * growth_cn) * 1.e6_wp - &
                   n_fixation * growth_cn*mavg_period_tgrowth*one_day) / & 
                   MAX(unit_uptake_n * growth_cn*mavg_period_tgrowth*one_day, eps8)
               
              k1_opt_p = (npp*mavg_period_tgrowth*one_day +&
                   (labile_carbon + reserve_carbon - &
                   (labile_phosphorus + reserve_phosphorus) * growth_cp) * 1.e6_wp) / & 
                   MAX(unit_uptake_p * growth_cp * mavg_period_tgrowth*one_day,eps8)
              IF ( k1_opt_n >= k1_opt_p) THEN 
                 ! if N is more limiting than P
                 k1_opt = k1_opt_n 
                 k2_opt = unit_npp/MAX(unit_uptake_n * growth_cn, eps8)
              ELSE
                 k1_opt = k1_opt_p 
                 k2_opt = unit_npp/MAX(unit_uptake_p * growth_cp, eps8)
              ENDIF
           CASE("optimal")
              k1_opt_n = (npp*mavg_period_tgrowth*one_day +&
                   (labile_carbon + reserve_carbon - &
                   (labile_nitrogen + reserve_nitrogen) * growth_cn) * 1.e6_wp) / & 
                   MAX((unit_uptake_n + n_fixation) * growth_cn*mavg_period_tgrowth*one_day, eps8)
              k1_opt_p = (npp*mavg_period_tgrowth*one_day +&
                   (labile_carbon + reserve_carbon - &
                   (labile_phosphorus + reserve_phosphorus) * growth_cp) * 1.e6_wp) / & 
                   MAX(unit_uptake_p * growth_cp * mavg_period_tgrowth*one_day,eps8)
              IF ( k1_opt_n >= k1_opt_p) THEN 
                 ! if N is more limiting than P
                 k1_opt = k1_opt_n 
                 k2_opt = unit_npp/MAX((unit_uptake_n+n_fixation) * growth_cn, eps8)
              ELSE
                 k1_opt = k1_opt_p 
                 k2_opt = unit_npp/MAX(unit_uptake_p * growth_cp, eps8)
              ENDIF
           END SELECT

           !> water based optimality, following Magnani 2000 
           !! R = k1 + k2*dL +k3*dS 
           k1_opt_w = transpiration / rhoh2o * (1._wp + lctlib%c0_allom * height) / (lctlib%k_root * dphi) 
           k2_opt_w = unit_transpiration/ rhoh2o * (1._wp + lctlib%c0_allom * height) / (lctlib%k_root * dphi) 

           ! Compare root increments for a hypothetical unit of leaf growth 
           IF (k1_opt + k2_opt < k1_opt_w + k2_opt_w) THEN 
              ! water is limiting 
              k1_opt = k1_opt_w 
              k2_opt = k2_opt_w 
           ENDIF 
        ELSE             
           ! scaling factor for root to leaf mass
           ! ad 1: derived from hydraulic homeostatis of root and wood hydraulic resistance (Magnani et al. 2000, PCE)
           ! Note: currently still lacks a moisture response!
           ! sapwood mass / height =  k_rtos * fine_root_mass
           ! Cr = Cs / k_rtos / height
           ! let Cs = Cl / leaf2sapwood_mass_ratio * height
           ! => Cl = Cr * leaf2sapwood_mass_ratio * k_rtos = Cl = Cr * leaf2root_mass_ratio 
           ! ad 2: nitrogen/phosphorus sensitivity still needed
           leaf2root_mass_ratio = lctlib%k_rtos * leaf2sapwood_mass_ratio
        ENDIF 
      !> 1.1 grass
      !! 
      CASE (igrass)
         ! leaf to stem mass ratio
         leaf2sapwood_mass_ratio = 1._wp / sm2lm_grass
         ! leaf:root ratio
         IF (veg_config%flag_optimal_allocation) THEN
            !> Calculate allocation fractions assuming optimal relationship between leaf and root allocation
            !! potential carbon growth rate = potential nitrogen growth rate * plant cn ratio
            !! NPP = gN * R * growth_cn
            !! Or similar relationships for P
               !! Difference between BNF schemes - for optimal scheme fixation depends on root mass
           SELECT CASE (veg_config%n_uptake_scheme)
           CASE ("prescribed", "standard","resistance")
              k1_opt_n = (npp*mavg_period_tgrowth*one_day +&
                   (labile_carbon + reserve_carbon - &
                   (labile_nitrogen + reserve_nitrogen) * growth_cn) * 1.e6_wp - &
                   n_fixation * growth_cn*mavg_period_tgrowth*one_day) / & 
                   MAX(unit_uptake_n * growth_cn*mavg_period_tgrowth*one_day, eps8)
           CASE("optimal")
              k1_opt_n = (npp*mavg_period_tgrowth*one_day +&
                   (labile_carbon + reserve_carbon - &
                   (labile_nitrogen + reserve_nitrogen) * growth_cn) * 1.e6_wp) / & 
                   MAX((unit_uptake_n + n_fixation) * growth_cn*mavg_period_tgrowth*one_day, eps8)
           END SELECT
               
           k1_opt_p = (npp*mavg_period_tgrowth*one_day +&
                (labile_carbon + reserve_carbon - &
                (labile_phosphorus + reserve_phosphorus) * growth_cp) * 1.e6_wp) / & 
                MAX(unit_uptake_p * growth_cp * mavg_period_tgrowth*one_day,eps8)
           IF ( k1_opt_n >= k1_opt_p) THEN 
              ! if N is more limiting than P
              k1_opt = k1_opt_n 
              k2_opt = unit_npp/MAX(unit_uptake_n * growth_cn, eps8)
           ELSE
              k1_opt = k1_opt_p 
              k2_opt = unit_npp/MAX(unit_uptake_p * growth_cp, eps8)
           ENDIF

            !> water based optimality, following Magnani 2000 
          
           k1_opt_w = transpiration / (lctlib%k_root * dphi) 
           k2_opt_w = unit_transpiration / (lctlib%k_root * dphi) 

           ! Compare root increments for a hypothetical unit of leaf growth 
           IF (k1_opt + k2_opt < k1_opt_w + k2_opt_w) THEN 
              ! water is limiting 
              k1_opt = k1_opt_w 
              k2_opt = k2_opt_w 
           ENDIF 

         ELSE      
            leaf2root_mass_ratio = sm2lm_grass * lctlib%k_rtos 
         ENDIF
         ! full leaf area - should depend on carbon balance of lowest leaf - works in principle but needs some filtering
         target_lai = 6.0_wp
    ENDSELECT

  END SUBROUTINE calc_allometry 


  !-----------------------------------------------------------------------------------------------------
  ! 4.0 Sub Task of update_growth
  ! 
  !----------------------------------------------------------------------------------------------------------
  !> Calculate flux between reserve and labile pools for each element 
  !! 
  !! Input:\n
  !!    (1) weekly GPP and maintenance respiration, and N,P growth requirements to calculate
  !!        target labile pool size;\n
  !!    (2) target LAI, leaf-to-root mass ratio and target C:N:P of foliage and fine roots to calculate
  !!        target reserve pool size;\n
  !!    (3) current LAI, labile and reserve pools, growing season status and potential turnover rate
  !!          of the labile pool;\n    
  !!  
  !! Output: net flux between reserve and labile pools 
  !----------------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE calc_reserve_useage(lctlib                  , &
                                                  growing_season        , &
                                                  kstar_labile          , &
                                                  gpp_mavg              , &
                                                  rm_mavg               , & 
                                                  growth_req_n_mavg     , &
                                                  growth_req_p_mavg     , &
                                                  target_cn_leaf        , &
                                                  target_cn_fine_root   , &
                                                  target_np_leaf        , &
                                                  target_np_fine_root   , & 
                                                  leaf2root_mass_ratio  , &
                                                  lai                   , &
                                                  target_lai            , &
                                                  veg_pool_labile       , &
                                                  veg_pool_reserve      , &
                                                  veg_reserve_use       , & 
                                                  beta_sinklim_ps )              

    USE mo_veg_memory_class,      ONLY: t_slm_struct_veg_element    ! *slmdev -> dsl4jsb macro [like dsl4jsb_Use_memory(VEG_)]
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),        INTENT(in)     :: lctlib                    !< land-cover-type library - parameters for one pft
    LOGICAL,                       INTENT(in)     :: growing_season            !< does the plant intend to grow?
    REAL(wp),                      INTENT(in)     :: kstar_labile              !< potential labile pool turnover (1/day) 
    REAL(wp),                      INTENT(in)     :: gpp_mavg                  !< last week's GPP (µmol CO2/m2/s) 
    REAL(wp),                      INTENT(in)     :: rm_mavg                   !< last week's Rm (µmol CO2/m2/s)
    REAL(wp),                      INTENT(in)     :: growth_req_n_mavg         !< last week's requirement for unit growth in C and N
    REAL(wp),                      INTENT(in)     :: growth_req_p_mavg         !< last week's rP requirement for unit growth in C and N
    REAL(wp),                      INTENT(in)     :: target_cn_leaf            !< target C:N and N:P rations for tissue type 
    REAL(wp),                      INTENT(in)     :: target_cn_fine_root       !< target C:N and N:P rations for tissue type 
    REAL(wp),                      INTENT(in)     :: target_np_leaf            !< target C:N and N:P rations for tissue type 
    REAL(wp),                      INTENT(in)     :: target_np_fine_root       !< target C:N and N:P rations for tissue type 
    REAL(wp),                      INTENT(in)     :: leaf2root_mass_ratio      !< leaf to root mass ratio
    REAL(wp),                      INTENT(in)     :: lai                       !< current LAI 
    REAL(wp),                      INTENT(in)     :: target_lai                !< target LAI 
    TYPE(t_slm_struct_veg_element),INTENT(in)     :: veg_pool_labile           !< the plant's labile pool size
    TYPE(t_slm_struct_veg_element),INTENT(in)     :: veg_pool_reserve          !< the plant's reserve pool size
    TYPE(t_slm_struct_veg_element),INTENT(inout)  :: veg_reserve_use           !< current resource use rate (mol/time step)
                                                                               !!  (+ve: net flux to labile / -ve: net flux to reserve)
    REAL(wp),                      INTENT(out)    :: beta_sinklim_ps           !< sink limitation factor for photosynthesis 
    ! 0.2 Local
    REAL(wp)             :: phi_phenology                    ! relative pull from reserve due to phenology
    REAL(wp)             :: phi_maintenance                  ! relative pull from reserve due to maintenance resp. 
    REAL(wp)             :: phi_storage                      ! relative pull from labile due to low reserve
    REAL(wp)             :: target_labile_pool_carbon        ! desirable size (mol C/m2) of the labile pool
    REAL(wp)             :: target_labile_pool_nitrogen      ! desirable size (mol N/m2) of the labile pool
    REAL(wp)             :: target_labile_pool_phosphorus    ! desirable size (mol P/m2) of the labile pool
    REAL(wp)             :: target_reserve_pool_carbon       ! desirable size (mol C/m2) of the reserve pool
    REAL(wp)             :: target_reserve_pool_nitrogen     ! desirable size (mol N/m2) of the reserve pool
    REAL(wp)             :: target_reserve_pool_phosphorus   ! desirable size (mol P/m2) of the reserve pool
    REAL(wp)             :: hlp1                             ! helper variable
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_reserve_useage'


    !> 1.0 target size of the labile and reserve pools
    !!

    !> 1.1 labile pool
    !! 
    !! This is taken as the average idaily input + maintenance requirements of the last week 
    !!  scaled to the default turnover time of the labile pool
    !! 
    !! for N and P this corresponds to the current growth C:N:P  
    target_labile_pool_carbon     = MAX( gpp_mavg , rm_mavg ) * one_day / 1.e6_wp * tau_labile
    target_labile_pool_nitrogen   = veg_pool_labile%carbon%ptr    * growth_req_n_mavg
    target_labile_pool_phosphorus = veg_pool_labile%nitrogen%ptr  * growth_req_p_mavg

    !> 1.2 reserve pool
    !!
    !! For the reserve pool, this is the C:N:P needed to regrow one year worth of folige and roots,
    !! accounting for losses to respiration and the nutrients currently required to grow these
    target_reserve_pool_carbon = target_lai / lctlib%sla * & 
           ( MIN(1.0_wp,1.0_wp/lctlib%tau_leaf) + leaf2root_mass_ratio/lctlib%tau_fine_root ) * (1._wp + fresp_growth)
    target_reserve_pool_nitrogen = target_lai / lctlib%sla * & 
           ( MIN(1.0_wp,1.0_wp/lctlib%tau_leaf) / target_cn_leaf + leaf2root_mass_ratio/lctlib%tau_fine_root / target_cn_fine_root) 
    target_reserve_pool_phosphorus = target_lai / lctlib%sla * & 
           ( MIN(1.0_wp,1.0_wp/lctlib%tau_leaf) / target_cn_leaf / target_np_leaf & 
             + leaf2root_mass_ratio/lctlib%tau_fine_root / target_cn_fine_root / target_np_fine_root)  

    !> 2.0 carbon reserve flow
    !! 
    !! calculate the pull from reserve to labile (phi_phenology,phi_maintenance) or
    !! 
    !! push from labile to reserve (phi_storage) depending on the difference of the LAI/labile/reserve from their targets
    !! 
    !! phi_storage is adjusted by phi_phenology to avoid storage building to interfere with flushing (push phi_storage towards 
    !!  reserve at times of low LAI or labile pool to give priority to growth and maintentance)
    phi_phenology   = 0._wp
    phi_maintenance = 0._wp
    phi_storage     = 0._wp
    hlp1 = MIN(MAX(min_lai,target_lai),target_lai_max)
    IF(lai < hlp1 .AND. growing_season) & 
       phi_phenology = exp(-(lambda_phiphen * lai / hlp1 )**k_phiphen) & 
                       - exp(-(lambda_phiphen)**k_phiphen) 

    hlp1 = MAX(0.0_wp,veg_pool_labile%carbon%ptr)
    IF(target_labile_pool_carbon > eps8 .AND. hlp1 < target_labile_pool_carbon) & 
       phi_maintenance = exp(-(lambda_phimaint * hlp1 / target_labile_pool_carbon)**k_phimaint) & 
                       - exp(-(lambda_phimaint)**k_phimaint)

    hlp1 = MAX(0.0_wp,veg_pool_reserve%carbon%ptr)
    IF(target_reserve_pool_carbon > eps8 .AND. hlp1 < target_reserve_pool_carbon) & 
       phi_storage = exp(-(lambda_phistore * hlp1 / target_reserve_pool_carbon)**k_phistore) & 
                       - exp(-(lambda_phistore)**k_phistore)

    phi_storage = phi_storage * ( 1.0_wp - phi_phenology) 
    IF(phi_maintenance > 0.75_wp) phi_storage = phi_storage * 4.0_wp * ( 1.0_wp - phi_maintenance) 

    !> 2.2 convert phi to turnover times (1/timestep)
    !! 
    !! taking account that flushing can only happen if the meristem is active
    phi_phenology = phi_phenology * kstar_labile/one_day * dtime
    phi_maintenance = phi_maintenance / tau_labile/one_day * dtime
    phi_storage = phi_storage / tau_labile/one_day * dtime

    veg_reserve_use%carbon%ptr   = ((phi_phenology + phi_maintenance) * veg_pool_reserve%carbon%ptr   - & 
                                 phi_storage   * veg_pool_labile%carbon%ptr) 
    veg_reserve_use%carbon13%ptr = ((phi_phenology + phi_maintenance) * veg_pool_reserve%carbon13%ptr - & 
                                 phi_storage   * veg_pool_labile%carbon13%ptr) 
    veg_reserve_use%carbon14%ptr = ((phi_phenology + phi_maintenance) * veg_pool_reserve%carbon14%ptr - & 
                                 phi_storage   * veg_pool_labile%carbon14%ptr) 

    !> 3.0 net transfer of N from labile to reserve
    !! 
    !! 
    phi_maintenance = 0._wp
    phi_storage = 0._wp
    hlp1 = MAX(0.0_wp,veg_pool_labile%nitrogen%ptr)
    IF(target_labile_pool_nitrogen > eps8 .AND. hlp1 < target_labile_pool_nitrogen) & 
       phi_maintenance = exp(-(lambda_phimaint * hlp1 / target_labile_pool_nitrogen)**k_phimaint) & 
                       - exp(-(lambda_phimaint)**k_phimaint)
    IF(target_reserve_pool_nitrogen > eps8 .AND. veg_pool_reserve%nitrogen%ptr < target_reserve_pool_nitrogen) & 
       phi_storage = exp(-(lambda_phistore * veg_pool_reserve%nitrogen%ptr / target_reserve_pool_nitrogen)**k_phistore) & 
                       - exp(-(lambda_phistore)**k_phistore)
       
    veg_reserve_use%nitrogen%ptr    = (phi_maintenance * veg_pool_reserve%nitrogen%ptr - & 
                                    phi_storage * veg_pool_labile%nitrogen%ptr) / tau_labile/one_day * dtime
    veg_reserve_use%nitrogen15%ptr  = (phi_maintenance * veg_pool_reserve%nitrogen15%ptr - & 
                                    phi_storage * veg_pool_labile%nitrogen15%ptr) / tau_labile/one_day * dtime

    !> 4.0 net transfer of P from labile to reserve
    !! 
    !! 
    phi_maintenance = 0._wp
    phi_storage = 0._wp
    hlp1 = MAX(0.0_wp,veg_pool_labile%phosphorus%ptr)
    IF(target_labile_pool_phosphorus > eps8 .AND. hlp1 < target_labile_pool_phosphorus) & 
       phi_maintenance = exp(-(lambda_phimaint * hlp1 / target_labile_pool_phosphorus)**k_phimaint) & 
                       - exp(-(lambda_phimaint)**k_phimaint)
    IF(target_reserve_pool_phosphorus > eps8 .AND. veg_pool_reserve%phosphorus%ptr < target_reserve_pool_phosphorus) & 
       phi_storage = exp(-(lambda_phistore * veg_pool_reserve%phosphorus%ptr / target_reserve_pool_phosphorus)**k_phistore) & 
                       - exp(-(lambda_phistore)**k_phistore)

    veg_reserve_use%phosphorus%ptr = (phi_maintenance * veg_pool_reserve%phosphorus%ptr - & 
                                  phi_storage * veg_pool_labile%phosphorus%ptr) / tau_labile/one_day * dtime

    !> 5.0 signal of carbon saturation to photosynthesis
    !!
    !!
    IF(target_labile_pool_carbon>eps8 .AND. veg_pool_labile%carbon%ptr>target_labile_pool_carbon)THEN
        hlp1=(veg_pool_labile%carbon%ptr-target_labile_pool_carbon)/target_labile_pool_carbon
        beta_sinklim_ps=beta_sinklim_ps_min + (1._wp - beta_sinklim_ps_min) * exp(-(lambda_sinklim_ps*hlp1)**k_sinklim_ps)
    ELSE 
        beta_sinklim_ps=1.0_wp
    ENDIF


  END SUBROUTINE calc_reserve_useage



  !-----------------------------------------------------------------------------------------------------
  ! 5.0 Sub Task of update_growth
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> Calculates allocation fractions for each tissue type and element
  !! 
  !! Input: target stoichiometry of the tissue types, plant allometry, and current pool sizes
  !! 
  !! Output:\n 
  !!         (1) Allocation coefficients for the distribution of biomass for growth;\n
  !!         (2) total N and P required to sustain a unit of growth (in mol); \n
  !!
  !! @todo: calc_allocation_fraction: check calculation of b_inc. Currently a fixed value for stability,
  !!        Needs to be adjusted to better match prescriebed allometry and to account for falloc_fruit
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_allocation_fraction(lctlib                      , &
                                      dheight                     , &
                                      kstar_labile                , &
                                      leaf_carbon                 , &
                                      sap_wood_carbon             , &
                                      coarse_root_carbon          , &
                                      root_carbon                 , &
                                      labile_carbon               , &
                                      height                      , &
                                      leaf2sapwood_mass_ratio     , &
                                      leaf2root_mass_ratio        , &
                                      k1_opt                      , &
                                      k2_opt                      , &
                                      target_cn_leaf              , &
                                      target_cn_fine_root         , &
                                      target_cn_coarse_root       , & 
                                      target_cn_sap_wood          , &
                                      target_cn_fruit             , &
                                      target_np_leaf              , &
                                      target_np_fine_root         , & 
                                      target_np_coarse_root       , &
                                      target_np_sap_wood          , &
                                      target_np_fruit             , &
                                      veg_reserve_use_carbon      , &
                                      veg_falloc                  , &       
                                      growth_req_n                , &     
                                      growth_req_p )       

    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_veg_memory_class,        ONLY: t_slm_struct_veg_pool     ! *slmdev -> dsl4jsb macro [like dsl4jsb_Use_memory(VEG_)]
    USE mo_veg_config_class,        ONLY: veg_config
    USE mo_veg_constants,           ONLY: itree, igrass, lambda_fruit_alloc, k1_fruit_alloc, k3_fruit_alloc, k4_fruit_alloc

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),       INTENT(in)    :: lctlib                     !< land-cover-type library - parameter across pft's
    REAL(wp),                     INTENT(in)    :: dheight                    !< derivative of height per unit sapwood mass change 
    REAL(wp),                     INTENT(in)    :: kstar_labile               !< potential labile pool turnover rate (1/day) without stoichiometric limit
    REAL(wp),                     INTENT(in)    :: leaf_carbon                !< carbon in leaves (mol/individual)
    REAL(wp),                     INTENT(in)    :: sap_wood_carbon            !< carbon in sap wood (mol/individual)
    REAL(wp),                     INTENT(in)    :: coarse_root_carbon         !< carbon in heart wood (mol/individual)
    REAL(wp),                     INTENT(in)    :: root_carbon                !< carbon in roots (mol/individual) 
    REAL(wp),                     INTENT(in)    :: labile_carbon              !< carbon in labile pool (mol/individual) 
    REAL(wp),                     INTENT(in)    :: height                     !< plant height  
    REAL(wp),                     INTENT(in)    :: leaf2sapwood_mass_ratio    !< leaf to shoot mass ratio
    REAL(wp),                     INTENT(in)    :: leaf2root_mass_ratio       !< leaf to root mass ratio
    REAL(wp),                     INTENT(in)    :: k1_opt                     !< coefficient for optimal leaf:root ratio
    REAL(wp),                     INTENT(in)    :: k2_opt                     !< coefficient for optimal leaf:root ratio
    REAL(wp),                     INTENT(in)    :: target_cn_leaf             !< target C:N ratio for tissue type 
    REAL(wp),                     INTENT(in)    :: target_cn_fine_root        !< target C:N ratio for tissue type 
    REAL(wp),                     INTENT(in)    :: target_cn_coarse_root      !< target C:N ratio for tissue type 
    REAL(wp),                     INTENT(in)    :: target_cn_sap_wood         !< target C:N ratio for tissue type 
    REAL(wp),                     INTENT(in)    :: target_cn_fruit            !< target C:N ratio for tissue type 
    REAL(wp),                     INTENT(in)    :: target_np_leaf             !< target N:P ratio for tissue type 
    REAL(wp),                     INTENT(in)    :: target_np_fine_root        !< target N:P ratio for tissue type 
    REAL(wp),                     INTENT(in)    :: target_np_coarse_root      !< target N:P ratio for tissue type 
    REAL(wp),                     INTENT(in)    :: target_np_sap_wood         !< target N:P ratio for tissue type 
    REAL(wp),                     INTENT(in)    :: target_np_fruit            !< target N:P ratio for tissue type 
    REAL(wp),                     INTENT(in)    :: veg_reserve_use_carbon     !< current reserve usage (mol C / m2 / timestep)
    TYPE(t_slm_struct_veg_pool),  INTENT(inout) :: veg_falloc                 !< allocation coefficients 
    REAL(wp),                     INTENT(out)   :: growth_req_n               !< N requirement for unit growth in C and N
    REAL(wp),                     INTENT(out)   :: growth_req_p               !< P requirement for unit growth in C and N
    ! 0.2 Local
    REAL(wp) :: b_inc,                   &             ! total and compartment-wise carbon increment (mol/dtime)
                all_inc,                 &
                leaf_carbon_inc,         &
                root_carbon_inc,         &
                sap_wood_carbon_inc,     &
                coarse_root_carbon_inc 
    REAL(wp) :: leaf2wood_mass_ratio,    &
                wood_carbon
    REAL(wp) :: hlp1

    !> 1.0 Allocation factors for C (this should be calculated based on various principles of plant growth and allometry)
    !! 
    !! 

    ! SZ: uncomment this to revert to fixed constant allocation
    ! veg_falloc%leaf%carbon = 0.3_wp
    ! veg_falloc%fine_root%carbon = 0.3_wp
    ! veg_falloc%coarse_root%carbon = 0.1_wp
    ! veg_falloc%sap_wood%carbon = 0.2_wp
    ! veg_falloc%fruit%carbon = 0.1_wp 

    ! 
    ! seed for allocation coefficients is the potential growth rate per time step.
    ! given by the average turnover of the labile pool 
    b_inc = kstar_labile * labile_carbon * dtime/one_day

    !> 1.1 Calculate allocation fraction to reproduction
    !!  Fruit allocation depends on whether or not the plant is able to build reserves
    !!  as soon as the carbon balance exceeds a PFT specific offset, the maximal fruit growth rate is achieved
    hlp1                    = MAX(0.0_wp, veg_reserve_use_carbon / dtime * 1000000.0_wp + k3_fruit_alloc)
    veg_falloc%fruit%carbon%ptr = k1_fruit_alloc + ( lctlib%k2_fruit_alloc - k1_fruit_alloc ) * & 
                              exp(-(lambda_fruit_alloc * hlp1)**k4_fruit_alloc)

    !> 1.2 Calculate allocation fractions to functional organs
    !!
    !! 
    SELECT CASE(lctlib%growthform)
      !> 1.2 trees
      !!
      CASE (itree)
       ! carbon increment
       IF(b_inc > eps8) THEN

          !! needed for stability - to be checked BLARPP! 
          !! @todo calc_allocation_fraction() -> 1.2 Calculate allocation fractions -> check stability issues with b_inc
          b_inc = 0.2_wp * dtime / one_day 

          leaf2wood_mass_ratio = leaf2sapwood_mass_ratio / (1._wp + lctlib%k_crtos)
          wood_carbon = sap_wood_carbon + coarse_root_carbon

          !! @todo calc_allocation_fraction() -> 1.2 Calculate allocation fractions -> 
          !! rewrite "IF (veg_config\%flag_optimal_allocation)" such that instead of kopt1 and 2 leaf2root_mass ratio can be given
          IF (veg_config%flag_optimal_allocation) THEN

             !> 1.2a Calculate allocation fractions assuming optimal relationship between leaf and root allocation
             !! 
             leaf_carbon_inc = MAX((b_inc - k1_opt + root_carbon + sap_wood_carbon - dheight/leaf2sapwood_mass_ratio * &
                  (b_inc-k1_opt + root_carbon) - height/leaf2sapwood_mass_ratio*leaf_carbon) / &
                  (1.0_wp + k2_opt + dheight/leaf2sapwood_mass_ratio * &
                  (b_inc - (1.0_wp + k2_opt) * leaf_carbon - k1_opt + root_carbon) + &
                  height/leaf2sapwood_mass_ratio), 0.0_wp)
             sap_wood_carbon_inc = MAX((( leaf_carbon + leaf_carbon_inc ) / leaf2sapwood_mass_ratio * height - sap_wood_carbon) / & 
                  (1._wp - ( leaf_carbon + leaf_carbon_inc ) / leaf2sapwood_mass_ratio * dheight ), 0.0_wp)        
             IF (sap_wood_carbon_inc > 0.0_wp) THEN
                coarse_root_carbon_inc = MAX((sap_wood_carbon + sap_wood_carbon_inc) * lctlib%k_crtos & 
                     - coarse_root_carbon, 0.0_wp)
                root_carbon_inc = MAX(k1_opt + k2_opt * leaf_carbon_inc - root_carbon, 0.0_wp)
             ELSE
                leaf_carbon_inc = MAX((b_inc - k1_opt + root_carbon)/(1.0_wp+k2_opt), 0.0_wp)
                root_carbon_inc = MAX(k1_opt + k2_opt * leaf_carbon_inc - root_carbon, 0.0_wp)
                sap_wood_carbon_inc = 0.0_wp
                coarse_root_carbon_inc = 0.0_wp
             ENDIF
          ELSE 

             !> 1.2b solving the allometric rations between leaf and root, leaf and sapwood and sapwood and coarse roots 
             !! simulataneously with the height-to-diameter relationship 
             leaf_carbon_inc = MAX(( b_inc + leaf_carbon * dheight / leaf2wood_mass_ratio * & 
                            ( leaf_carbon / leaf2root_mass_ratio - b_inc - root_carbon ) + wood_carbon + root_carbon - & 
                            leaf_carbon * (height / leaf2wood_mass_ratio + 1.0_wp / leaf2root_mass_ratio )) / & 
                            (1.0_wp + ( b_inc - (1.0_wp + 2.0_wp/leaf2root_mass_ratio) * leaf_carbon + root_carbon ) & 
                            * dheight / leaf2wood_mass_ratio + height / leaf2wood_mass_ratio + 1._wp / leaf2root_mass_ratio ), & 
                            0.0_wp) 
   
             sap_wood_carbon_inc = MAX((( leaf_carbon + leaf_carbon_inc ) / leaf2sapwood_mass_ratio * height - sap_wood_carbon) / & 
                                (1._wp - ( leaf_carbon + leaf_carbon_inc ) / leaf2sapwood_mass_ratio * dheight ), 0.0_wp)          
             IF(sap_wood_carbon_inc > 0.0_wp) THEN
                coarse_root_carbon_inc = MAX((sap_wood_carbon + sap_wood_carbon_inc) * lctlib%k_crtos & 
                                           - coarse_root_carbon, 0.0_wp)
                root_carbon_inc = MAX(( leaf_carbon + leaf_carbon_inc ) / leaf2root_mass_ratio - root_carbon, 0.0_wp)
             ELSE ! special case of no wood growth (phenology etc.): redistribute available carbon 
                leaf_carbon_inc = MAX(( b_inc - leaf_carbon / leaf2root_mass_ratio + root_carbon ) / & 
                                    (1.0_wp + 1._wp/ leaf2root_mass_ratio),0.0_wp)
                root_carbon_inc = MAX(0.0_wp,MIN(b_inc - leaf_carbon_inc, & 
                                       ( leaf_carbon_inc + leaf_carbon ) / leaf2root_mass_ratio - root_carbon))
                sap_wood_carbon_inc = 0.0_wp
                coarse_root_carbon_inc = 0.0_wp
             ENDIF
          ENDIF
        ENDIF
      !> 1.2 grass
      !!
      CASE(igrass)
         IF (veg_config%flag_optimal_allocation) THEN
            leaf_carbon_inc = MAX((b_inc - k1_opt - leaf_carbon/leaf2sapwood_mass_ratio + root_carbon + &
                 sap_wood_carbon) / (1.0_wp + k2_opt + 1.0_wp/leaf2sapwood_mass_ratio), 0.0_wp)
            root_carbon_inc = MAX(k1_opt + k2_opt * leaf_carbon_inc - root_carbon, 0.0_wp)
            sap_wood_carbon_inc = ( leaf_carbon_inc + leaf_carbon ) / leaf2sapwood_mass_ratio - sap_wood_carbon
            coarse_root_carbon_inc = 0.0_wp
         ELSE
            leaf_carbon_inc = MAX( ( b_inc - leaf_carbon / leaf2root_mass_ratio + sap_wood_carbon + root_carbon ) / & 
                 (1.0_wp + 1._wp / leaf2sapwood_mass_ratio + 1._wp / leaf2root_mass_ratio ) , 0.0_wp) 
            root_carbon_inc = MAX(MIN( b_inc, ( leaf_carbon + leaf_carbon_inc ) / leaf2root_mass_ratio - root_carbon), 0.0_wp)
            sap_wood_carbon_inc = ( leaf_carbon_inc + leaf_carbon ) / leaf2sapwood_mass_ratio - sap_wood_carbon
            coarse_root_carbon_inc = 0.0_wp
         ENDIF
      !> 1.2 default - write and finish
      !!
      CASE DEFAULT
       !! Write statement
       !! 
       WRITE(message_text,'(a,i4)') 'invalid growthform type ',lctlib%growthform
       CALL finish("mo_veg_growth",message_text)
    END SELECT

    !> 1.3 Convert biomass increments to allocation fractions
    !! taking account of the previously calculated fraction of growth that goes to fruits
    !!
    IF(b_inc > eps8) THEN
       all_inc = leaf_carbon_inc + sap_wood_carbon_inc + root_carbon_inc + coarse_root_carbon_inc
       IF(all_inc > eps8) THEN 
          veg_falloc%leaf%carbon%ptr        = (1._wp - veg_falloc%fruit%carbon%ptr ) * leaf_carbon_inc / all_inc
          veg_falloc%fine_root%carbon%ptr   = (1._wp - veg_falloc%fruit%carbon%ptr ) * root_carbon_inc / all_inc
          veg_falloc%coarse_root%carbon%ptr = (1._wp - veg_falloc%fruit%carbon%ptr ) * coarse_root_carbon_inc / all_inc
          veg_falloc%sap_wood%carbon%ptr    = (1._wp - veg_falloc%fruit%carbon%ptr ) * sap_wood_carbon_inc/ all_inc
       ENDIF
    ENDIF

    IF(b_inc > eps8) THEN

       !> 2.0 Allocation factors for N, derived from the target C:N ratios and the C partitioning coefficients
       !!  
 
       !> 2.1 N required to sustain a unit of C growth
       !!
       growth_req_n = veg_falloc%leaf%carbon%ptr        / target_cn_leaf        + & 
                      veg_falloc%fine_root%carbon%ptr   / target_cn_fine_root   + &
                      veg_falloc%coarse_root%carbon%ptr / target_cn_coarse_root + &
                      veg_falloc%sap_wood%carbon%ptr    / target_cn_sap_wood    + & 
                      veg_falloc%fruit%carbon%ptr       / target_cn_fruit 
   
       !> 2.2 Adjustment of N allocation factors
       !!
       IF(growth_req_n > eps8)THEN
          veg_falloc%leaf%nitrogen%ptr        = veg_falloc%leaf%carbon%ptr        / target_cn_leaf        / growth_req_n
          veg_falloc%fine_root%nitrogen%ptr   = veg_falloc%fine_root%carbon%ptr   / target_cn_fine_root   / growth_req_n 
          veg_falloc%coarse_root%nitrogen%ptr = veg_falloc%coarse_root%carbon%ptr / target_cn_coarse_root / growth_req_n
          veg_falloc%sap_wood%nitrogen%ptr    = veg_falloc%sap_wood%carbon%ptr    / target_cn_sap_wood    / growth_req_n
          veg_falloc%fruit%nitrogen%ptr       = veg_falloc%fruit%carbon%ptr       / target_cn_fruit       / growth_req_n
       ENDIF

       !> 3.0 Allocation factors for P, derived from the target N:P ratios and the N partitioning coefficients
       !! 
       !! 
       ! P required to sustain a unit of N growth
       growth_req_p = veg_falloc%leaf%nitrogen%ptr        / target_np_leaf        + & 
                      veg_falloc%fine_root%nitrogen%ptr   / target_np_fine_root   + &
                      veg_falloc%coarse_root%nitrogen%ptr / target_np_coarse_root + &
                      veg_falloc%sap_wood%nitrogen%ptr    / target_np_sap_wood    + & 
                      veg_falloc%fruit%nitrogen%ptr       / target_np_fruit 
   
       ! Adjustment of P allocation factors
       IF(growth_req_p > eps8)THEN
          veg_falloc%leaf%phosphorus%ptr        = veg_falloc%leaf%nitrogen%ptr        / target_np_leaf        / growth_req_p
          veg_falloc%fine_root%phosphorus%ptr   = veg_falloc%fine_root%nitrogen%ptr   / target_np_fine_root   / growth_req_p 
          veg_falloc%coarse_root%phosphorus%ptr = veg_falloc%coarse_root%nitrogen%ptr / target_np_coarse_root / growth_req_p
          veg_falloc%sap_wood%phosphorus%ptr    = veg_falloc%sap_wood%nitrogen%ptr    / target_np_sap_wood    / growth_req_p
          veg_falloc%fruit%phosphorus%ptr       = veg_falloc%fruit%nitrogen%ptr       / target_np_fruit       / growth_req_p
       ENDIF
   
       !> 4.0 Allocation factors for isotopes follow the main isotope 
       !! 
       !! 
       veg_falloc%leaf%carbon13%ptr         = veg_falloc%leaf%carbon%ptr
       veg_falloc%fine_root%carbon13%ptr    = veg_falloc%fine_root%carbon%ptr
       veg_falloc%coarse_root%carbon13%ptr  = veg_falloc%coarse_root%carbon%ptr
       veg_falloc%sap_wood%carbon13%ptr     = veg_falloc%sap_wood%carbon%ptr
       veg_falloc%fruit%carbon13%ptr        = veg_falloc%fruit%carbon%ptr
   
       veg_falloc%leaf%carbon14%ptr         = veg_falloc%leaf%carbon%ptr
       veg_falloc%fine_root%carbon14%ptr    = veg_falloc%fine_root%carbon%ptr
       veg_falloc%coarse_root%carbon14%ptr  = veg_falloc%coarse_root%carbon%ptr
       veg_falloc%sap_wood%carbon14%ptr     = veg_falloc%sap_wood%carbon%ptr
       veg_falloc%fruit%carbon14%ptr        = veg_falloc%fruit%carbon%ptr
   
       veg_falloc%leaf%nitrogen15%ptr         = veg_falloc%leaf%nitrogen%ptr
       veg_falloc%fine_root%nitrogen15%ptr    = veg_falloc%fine_root%nitrogen%ptr
       veg_falloc%coarse_root%nitrogen15%ptr  = veg_falloc%coarse_root%nitrogen%ptr
       veg_falloc%sap_wood%nitrogen15%ptr     = veg_falloc%sap_wood%nitrogen%ptr
       veg_falloc%fruit%nitrogen15%ptr        = veg_falloc%fruit%nitrogen%ptr
    ENDIF

  END SUBROUTINE calc_allocation_fraction


  !-----------------------------------------------------------------------------------------------------
  ! 6.0 Sub Task of update_growth
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> Calculate partitioning of labile pool outflow into maintenance, growth respiration and growth 
  !! Follows the concept of Amthor: Growth = Yg * (P - Rm), where\n
  !!   Yg = growth efficiency of biomass production (1-fresp_growth)\n
  !!   P  = current turnover of labile pool (Amthor's production)\n
  !!   Rm = maintenance respiration
  !!
  !! In the current implementation, maintenance is given preference over growth, but 
  !! nevertheless, maintenance is limited to the available C   
  !!
  !! Input: veg
  !! 
  !! Output: growth of tissue types, maintenance respiration and growth respiration 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_partitioning(t_air                          , &  
                                t_air_tacclim_mavg            , &
                                kstar_labile                  , &
                                maint_respiration_leaf        , &
                                cost_n_fixation               , &
                                growth_req_p                  , &
                                growth_req_n                  , &
                                veg_falloc                    , &
                                veg_pool                      , &
                                uptake_nh4                    , &
                                uptake_no3                    , &
                                n_fixation                    , &
                                n_transform_respiration       , &
                                n_fixation_respiration        , &
                                n_processing_respiration      , &
                                maint_respiration_c13         , &
                                maint_respiration_c14         , &
                                growth_respiration_c13        , &
                                growth_respiration_c14        , &
                                n_processing_respiration_c13  , &
                                n_processing_respiration_c14  , &
                                veg_growth                    , &
                                maint_respiration_pot         , &
                                maint_respiration             , &
                                growth_respiration  )  

    USE mo_veg_respiration,         ONLY: calc_maintenance_respiration
    USE mo_veg_config_class,        ONLY: veg_config
    USE mo_veg_memory_class,        ONLY: t_slm_struct_veg_pool

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),                     INTENT(in)    :: t_air                        , & !< air temperature
                                                   t_air_tacclim_mavg           , & !< .. 
                                                   kstar_labile                 , & !< potential labile pool turnover (1/day) 
                                                   maint_respiration_leaf       , & !< foliar maintenance respiration
                                                   cost_n_fixation              , & !< actual cost for N uptake through fixation
                                                   growth_req_p                 , & !< moles P required for a unit of N growth under current allocation fractions
                                                   growth_req_n                     !< moles N required for a unit of C growth under current allocation fractions
    TYPE(t_slm_struct_veg_pool),  INTENT(in)    :: veg_falloc                   , & !< fractional allocation to pools
                                                   veg_pool                         !< pools of the plant
    REAL(wp),                     INTENT(inout) :: uptake_nh4                   , & !< NH4 uptake from soil
                                                   uptake_no3                   , & !< NO3 uptake from soil
                                                   n_fixation                   , & !< symbiontic N fixation
                                                   n_transform_respiration      , & !< respiration associated with N transformation
                                                   n_fixation_respiration       , & !< respiration associated with N fixation
                                                   n_processing_respiration     , & !< respiration associated with N processing
                                                   maint_respiration_c13        , & !< maintenance respiration
                                                   maint_respiration_c14        , & !< maintenance respiration
                                                   growth_respiration_c13       , & !< growth respiration
                                                   growth_respiration_c14       , & !< growth respiration
                                                   n_processing_respiration_c13 , & !< respiration associated with N processing
                                                   n_processing_respiration_c14     !< respiration associated with N processing
    TYPE(t_slm_struct_veg_pool),  INTENT(inout) :: veg_growth                       !< fluxes: growth
    REAL(wp),                     INTENT(out)   :: maint_respiration_pot        , & !< potential maintenance respiration
                                                   maint_respiration            , & !< maintenance respiration
                                                   growth_respiration               !< growth respiration
    ! 0.2 Local
    REAL(wp)                                    :: n_processing_respiration_pot , & ! ..
                                                   n_processing_respiration_act     ! ..
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_partitioning' 


    !> 1.0 Determine C use for maintenance 
    !! 
    !! based on tissue N and Temperature, constrained by default turnover of the labile pool
    !! to avoid that maitenance respiration consumes all C in times of low GPP 

    ! potential rate of maintenance respiration (µmol C / m2 / s)
    maint_respiration_pot = calc_maintenance_respiration(               &
                                     &      t_air                             , &  ! t_stem
                                     &      t_air                             , &  ! t_root
                                     &      t_air_tacclim_mavg                , &  ! average temperature for acclimation (K) 
                                     &      veg_pool%sap_wood%nitrogen%ptr    , & 
                                     &      veg_pool%fine_root%nitrogen%ptr   , &
                                     &      veg_pool%coarse_root%nitrogen%ptr)  & 
                                     & + maint_respiration_leaf
    ! actual rate of maintentance respiration (µmol C / m2 / s), constrained to 
    ! support leaf respiration, but not more than the current labile pool turnover
    maint_respiration = MAX(maint_respiration_leaf                                      , & 
                                      &   MIN(1.0_wp/tau_labile/one_day * veg_pool%labile%carbon%ptr * 1.e6_wp, &
                                      &       maint_respiration_pot)                              &
                                      &   )

    !> 1.5 Determine C use for N uptake
    IF(veg_config%n_uptake_scheme == 'optimal' .OR. veg_config%n_uptake_scheme == 'resistance' .OR. veg_config%n_uptake_scheme == 'std_plus_costs') THEN
      n_processing_respiration_pot = uptake_nh4 * transform_cost_nh4 + uptake_no3 * transform_cost_no3 + & 
                                     n_fixation * cost_n_fixation
      n_processing_respiration_act = MIN(1.0_wp/tau_labile/one_day * veg_pool%labile%carbon%ptr * 1.e6_wp, & 
                                          n_processing_respiration_pot)
       ! adjust N uptake to C expenditure
       IF(n_processing_respiration_pot > eps8 .AND. n_processing_respiration_pot > n_processing_respiration_act) THEN 
          uptake_nh4    = uptake_nh4  * n_processing_respiration_act / n_processing_respiration_pot
          uptake_no3    = uptake_no3  * n_processing_respiration_act / n_processing_respiration_pot
          n_fixation  = n_fixation * n_processing_respiration_act / n_processing_respiration_pot
       ENDIF
       n_transform_respiration  = uptake_nh4 * transform_cost_nh4 + uptake_no3 * transform_cost_no3
       n_fixation_respiration   = n_fixation * cost_n_fixation
       n_processing_respiration = n_fixation_respiration + n_transform_respiration
!    ELSEIF(veg_config%n_uptake_scheme == 'std_plus_costs') THEN
!      n_processing_respiration_pot = uptake_nh4 * transform_cost_nh4 + uptake_no3 * transform_cost_no3 
!      n_processing_respiration_act = MIN(1.0_wp/tau_labile/one_day * veg_pool%labile%carbon%ptr * 1.e6_wp, & 
!                                         n_processing_respiration_pot)
       ! adjust N uptake to C expenditure
       IF(n_processing_respiration_pot > eps8 .AND. n_processing_respiration_pot > n_processing_respiration_act) THEN 
          uptake_nh4    = uptake_nh4  * n_processing_respiration_act / n_processing_respiration_pot
          uptake_no3    = uptake_no3  * n_processing_respiration_act / n_processing_respiration_pot
       ENDIF
       n_transform_respiration  = uptake_nh4 * transform_cost_nh4 + uptake_no3 * transform_cost_no3
       n_processing_respiration = n_fixation_respiration + n_transform_respiration
    ENDIF

    !> 2.0 Determine C, N, P use for growth
    !!
    !> 2.1 potential C growth rate (mol C / m2 / timestep) 
    !! is determined from the labile turnover minus the carbon used for maintenance, i
    !! accounting for the fraction of growth respiration required to achieve this growth, to ensure that
    !! growth does not empty the labile pool
    veg_growth%labile%carbon%ptr = (veg_pool%labile%carbon%ptr * kstar_labile / one_day - &
                                    &  (maint_respiration + n_processing_respiration) / 1.e6_wp ) * & 
                                    &   dtime /  (1._wp + fresp_growth) 
    IF(veg_growth%labile%carbon%ptr < eps8) veg_growth%labile%carbon%ptr = 0.0_wp

    !> 2.2 Elemental limitation of growth by N & P (mol [N|P] / m2 / timestep)
    !! is determined by the minimum turnover of the labile N and P pools and the stoichiometric growth requirement
    !! given the potential growth rate of C
    veg_growth%labile%nitrogen%ptr   = MIN(growth_req_n * veg_growth%labile%carbon%ptr, & 
                                            veg_pool%labile%nitrogen%ptr * kstar_labile * dtime/one_day) 
    veg_growth%labile%phosphorus%ptr = MIN(growth_req_p * veg_growth%labile%nitrogen%ptr, & 
                                            veg_pool%labile%phosphorus%ptr * kstar_labile * dtime/one_day) 
    ! adjust N growth rate to P growth rate if P is limiting 
    IF(growth_req_p > eps8) THEN
      IF(veg_growth%labile%nitrogen%ptr > veg_growth%labile%phosphorus%ptr / growth_req_p) & 
        veg_growth%labile%nitrogen%ptr = veg_growth%labile%phosphorus%ptr / growth_req_p
    END IF

    ! adjust C growth rate to N growth rate if N is limiting 
    IF(growth_req_n > eps8) THEN
      IF(veg_growth%labile%carbon%ptr > veg_growth%labile%nitrogen%ptr / growth_req_n) & 
        veg_growth%labile%carbon%ptr = veg_growth%labile%nitrogen%ptr / growth_req_n
    END IF
 
    ! Diagnose growth respiration (µmol C / m2 / s)
    growth_respiration = fresp_growth * veg_growth%labile%carbon%ptr * 1.e6_wp / dtime

    !> 3.0 Calculate isotopic signature of respiration fluxes, and growth
    !!
    IF(veg_pool%labile%carbon%ptr > eps8) THEN
       ! growth (in mol C / m2 / timestep)
       veg_growth%labile%carbon13%ptr = veg_growth%labile%carbon%ptr * veg_pool%labile%carbon13%ptr/veg_pool%labile%carbon%ptr
       veg_growth%labile%carbon14%ptr = veg_growth%labile%carbon%ptr * veg_pool%labile%carbon14%ptr/veg_pool%labile%carbon%ptr
       ! respiration (in µmol C / m2 / s)
       maint_respiration_c13  = maint_respiration    * veg_pool%labile%carbon13%ptr/veg_pool%labile%carbon%ptr
       maint_respiration_c14  = maint_respiration    * veg_pool%labile%carbon14%ptr/veg_pool%labile%carbon%ptr
       growth_respiration_c13 = growth_respiration   * veg_pool%labile%carbon13%ptr/veg_pool%labile%carbon%ptr
       growth_respiration_c14 = growth_respiration   * veg_pool%labile%carbon14%ptr/veg_pool%labile%carbon%ptr
       n_processing_respiration_c13 = n_processing_respiration * & 
                                                  veg_pool%labile%carbon13%ptr / veg_pool%labile%carbon%ptr 
       n_processing_respiration_c14 = n_processing_respiration * & 
                                                  veg_pool%labile%carbon14%ptr / veg_pool%labile%carbon%ptr 
    ENDIF
    IF(veg_pool%labile%nitrogen%ptr > eps8) THEN
       ! growth (in mol  / m2 / timestep)
       veg_growth%labile%nitrogen15%ptr = veg_growth%labile%nitrogen%ptr * &
                                              veg_pool%labile%nitrogen15%ptr/veg_pool%labile%nitrogen%ptr
    ENDIF

    !> 4.0 partition growth to tissue types
    !! 
    veg_growth%leaf        = veg_falloc%leaf        * veg_growth%labile
    veg_growth%fine_root   = veg_falloc%fine_root   * veg_growth%labile 
    veg_growth%coarse_root = veg_falloc%coarse_root * veg_growth%labile
    veg_growth%sap_wood    = veg_falloc%sap_wood    * veg_growth%labile
    veg_growth%fruit       = veg_falloc%fruit       * veg_growth%labile

    !> 4.1 invert flux of labile to be out of the labile pool
    !! 
    veg_growth%labile%carbon     = veg_growth%labile%carbon     * (-1.0_wp)
    veg_growth%labile%nitrogen   = veg_growth%labile%nitrogen   * (-1.0_wp)
    veg_growth%labile%phosphorus = veg_growth%labile%phosphorus * (-1.0_wp)
    veg_growth%labile%carbon13   = veg_growth%labile%carbon13   * (-1.0_wp)
    veg_growth%labile%carbon14   = veg_growth%labile%carbon14   * (-1.0_wp)
    veg_growth%labile%nitrogen15 = veg_growth%labile%nitrogen15 * (-1.0_wp)

  END SUBROUTINE calc_partitioning 



  !-----------------------------------------------------------------------------------------------------
  ! Additonal function used in:
  !   calc_allometry
  !   update_veg_pools
  ! 
  !----------------------------------------------------------------------------------------------------------
  !> Function to calculate diameter of a tree given its biomass
  !! 
  !! Diameter and Height calculated assuming that: \n
  !!   all woody biomass forms a cylinder: diameter = SQRT(4 * wood volume / height / pi) \n
  !!   height = allom_k1 * diameter ^ allom_k2
  !----------------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_diameter_from_woody_biomass( lctlib            , &
                                                            sap_wood_carbon   , &
                                                            heart_wood_carbon , &
                                                            dens_ind)           &
                                                            RESULT(diameter)

    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_base_constants,          ONLY: pi

    IMPLICIT NONE

    !0.1 InOut
    TYPE(t_lctlib_element), INTENT(in)  :: lctlib                  !< land-cover-type library - parameters for one pft
    REAL(wp),               INTENT(in)  :: sap_wood_carbon         !< woody biomass of an individual tree (mol C/m2)
    REAL(wp),               INTENT(in)  :: heart_wood_carbon       !< woody biomass of an individual tree (mol C/m2)
    REAL(wp),               INTENT(in)  :: dens_ind                !< density of populations (#indivuduals /m2) 
    REAL(wp)                            :: diameter                !< diameter of an individual (m2)
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_diameter_from_woody_biomass'


    diameter = ( ( sap_wood_carbon + heart_wood_carbon ) / dens_ind / & 
               (lctlib%wood_density * pi/4._wp * lctlib%allom_k1)) ** (1._wp/(2._wp + lctlib%allom_k2))

  END FUNCTION calc_diameter_from_woody_biomass


  !-----------------------------------------------------------------------------------------------------
  ! Additonal function used in:
  !   calc_allometry
  !   update_veg_pools
  ! 
  !----------------------------------------------------------------------------------------------------------
  !> Function to calculate height of a tree given its diameter 
  !!
  !! height = allom_k1 * diameter ^ allom_k2
  !----------------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_height_from_diameter(lctlib, diameter) RESULT(height)

    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_veg_constants,           ONLY: min_height

    IMPLICIT NONE

    !0.1 InOut
    TYPE(t_lctlib_element), INTENT(in)  :: lctlib      !< land-cover-type library - parameters for one pft
    REAL(wp),               INTENT(in)  :: diameter    !< diameter of an individual (m2) 
    REAL(wp)                            :: height      !< height of an individual (m)
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_height_from_diameter'


    height = MAX(min_height,lctlib%allom_k1 * (diameter ** lctlib%allom_k2))

  END FUNCTION calc_height_from_diameter

END MODULE mo_veg_growth
