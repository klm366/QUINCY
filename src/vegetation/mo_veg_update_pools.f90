!>
!! @brief Contains the routines for the Task update_canopy_fluxes (photosynthesis related processes)
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!! 2nd version by Jan Engel (20160719): svn move vegetation/mo_veg_photosynthesis.f90 assimilation/mo_assimi_update_canopy_fluxes.f90
!!
MODULE mo_veg_update_pools

#include "dsl_jsbach.inc"

  USE mo_kind
  USE mo_jsb_lctlib_class,    ONLY: t_lctlib_element


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_veg_pools

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_update_pools'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update vegetation pools
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_veg_pools(tile, options)

    USE mo_jsb_class,           ONLY: Get_model
    USE mo_base_memory,         ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,     ONLY: t_jsb_model
    USE mo_base_constants,      ONLY: dtime, one_day
    USE mo_base_config_class,   ONLY: base_config
    USE mo_jsb_process_class,   ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_shared_constants,    ONLY: lambda_C14, rhoh2o, grav
    USE mo_veg_canopy,          ONLY: calc_canopy_layers
    USE mo_veg_growth,          ONLY: calc_diameter_from_woody_biomass, calc_height_from_diameter
    USE mo_veg_constants,       ONLY: itree

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(ASSIMI_)    
    
    ! Use of process memories
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4slm_Use_memory_a_veg
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4slm_Use_memory_a_sb
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout) :: tile             !< one tile
    TYPE(t_jsb_task_options), INTENT(in)    :: options          !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: nsoil_sb
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_pools'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_memory(ASSIMI_)                       ! TYPE(t_assimi_memory),      POINTER :: assimi__mem
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(SB_)                           ! TYPE(t_sb_memory),      POINTER :: sb__mem
    dsl4jsb_Def_memory(RAD_)                          !     TYPE(t_rad_memory),      POINTER :: rad__mem
    dsl4jsb_Def_memory(SOIL_)                         !     TYPE(t_soil_memory),     POINTER :: soil__mem
    ! Declare pointers to variables in memory
    ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: gross_assimilation
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C13
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C14
    dsl4jsb_Real2D_onChunk      :: beta_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tfrac_mavg
    ! SOIL_
    dsl4jsb_Real2D_onChunk      :: w_soil_root_pot
    ! SOIL_ 3D
    dsl4jsb_Real3D_onChunk      :: soil_depth_sl
    ! RAD_
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tfrac_mavg_cl
    ! VEG_
    dsl4jsb_Real2D_onChunk      :: lai
    dsl4jsb_Real2D_onChunk      :: dens_ind
    dsl4jsb_Real2D_onChunk      :: delta_dens_ind
    dsl4jsb_Real2D_onChunk      :: diameter
    dsl4jsb_Real2D_onChunk      :: height
    dsl4jsb_Real2D_onChunk      :: mean_leaf_age
    dsl4jsb_Real2D_onChunk      :: dphi
    dsl4jsb_Real2D_onChunk      :: npp
    dsl4jsb_Real2D_onChunk      :: growth_respiration
    dsl4jsb_Real2D_onChunk      :: maint_respiration
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration
    dsl4jsb_Real2D_onChunk      :: npp_c13
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c13
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c13
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration_c13
    dsl4jsb_Real2D_onChunk      :: npp_c14
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c14
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c14
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration_c14
    dsl4jsb_Real2D_onChunk      :: uptake_nh4
    dsl4jsb_Real2D_onChunk      :: uptake_no3
    dsl4jsb_Real2D_onChunk      :: n_fixation
    dsl4jsb_Real2D_onChunk      :: uptake_nh4_n15
    dsl4jsb_Real2D_onChunk      :: uptake_no3_n15
    dsl4jsb_Real2D_onChunk      :: n_fixation_n15
    dsl4jsb_Real2D_onChunk      :: uptake_po4
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_p
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n15
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_p
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_p
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n15
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n15
    dsl4jsb_Real2D_onChunk      :: t_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tacclim_mavg
    dsl4jsb_Real2D_onChunk      :: press_srf_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: ga_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: lai_cl
    dsl4jsb_Real3D_onChunk      :: cumm_lai_cl
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl
    ! Declare pool structures in memory
#define Def_veg_pool
#define Def_veg_litterfall
#define Def_veg_growth
#define Def_veg_establishment
#define Def_veg_exudation
#define Def_veg_total_biomass
#define Def_veg_reserve_use
#define Def_sb_mycorrhiza_export
#include "slm_pools_a_incl.f90"
#undef Def_veg_pool
#undef Def_veg_litterfall
#undef Def_veg_growth
#undef Def_veg_establishment
#undef Def_veg_exudation
#undef Def_veg_total_biomass
#undef Def_veg_reserve_use
#undef Def_sb_mycorrhiza_export
    ! Get local variables from options argument
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    nsoil_sb  = options%nsoil_sb
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(VEG_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! Get process config
    dsl4jsb_Get_config(ASSIMI_)
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SOIL_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, gross_assimilation)          ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, gross_assimilation_C13)      ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, gross_assimilation_C14)      ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_tfrac_mavg)         ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_tfrac_mavg)     ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tfrac_mavg)     ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_pot)               ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_, soil_depth_sl)                 ! in
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tfrac_mavg_cl)      ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tfrac_mavg_cl)      ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)                            ! CASE "plant" & "land": out / other Cases: in
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)                       ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, delta_dens_ind)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, diameter)                       ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, height)                         ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, mean_leaf_age)                  ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, dphi)                           ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, npp)                            ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration)       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_c13)                        ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c13)         ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c13)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration_c13)   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_c14)                        ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c14)         ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c14)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration_c14)   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4_n15)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3_n15)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_n15)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_po4)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_p)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n15)        ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n)         ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_p)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_p)         ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n15)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n15)       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tfrac_mavg)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tacclim_mavg)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_tfrac_mavg)           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_tfrac_mavg)    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_tfrac_mavg)                  ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tfrac_mavg)     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)                ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl)     ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)                         ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, cumm_lai_cl)                    ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)               ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_rub_cl)                      ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_et_cl)                       ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_pepc_cl)                     ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)                      ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_oth_cl)                      ! inout
    ! Get process memory pool structures
#define Get_veg_pool
#define Get_veg_litterfall
#define Get_veg_growth
#define Get_veg_establishment
#define Get_veg_exudation
#define Get_veg_total_biomass
#define Get_veg_reserve_use
#define Get_sb_mycorrhiza_export
#include "slm_pools_a_incl.f90"
#undef Get_veg_pool
#undef Get_veg_litterfall
#undef Get_veg_growth
#undef Get_veg_establishment
#undef Get_veg_exudation
#undef Get_veg_total_biomass
#undef Get_veg_reserve_use
#undef Get_sb_mycorrhiza_export



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 for plant and land model, update vegetation pools with fluxes
    !! 
    !!
    SELECT CASE(TRIM(base_config%slm_model_name))
    CASE('plant','land')
   
       !> 1.1 net balance of the labile carbon pool due to photosynthesis and respiration
       !! 
       !! 
       npp(:)      = gross_assimilation(:)     - & 
           growth_respiration(:) - maint_respiration(:) - n_processing_respiration(:) 
       npp_c13(:)  = gross_assimilation_C13(:) - &
           growth_respiration_c13(:) - maint_respiration_c13(:) - n_processing_respiration_c13(:) 
       npp_c14(:)  = gross_assimilation_C14(:) - & 
           growth_respiration_c14(:) - maint_respiration_c14(:) - n_processing_respiration_c14(:)

       veg_growth%labile%carbon(:)   = veg_growth%labile%carbon(:)   + npp(:)     * dtime / 1000000.0_wp
       veg_growth%labile%carbon13(:) = veg_growth%labile%carbon13(:) + npp_c13(:) * dtime / 1000000.0_wp
       veg_growth%labile%carbon14(:) = veg_growth%labile%carbon14(:) + npp_c14(:) * dtime / 1000000.0_wp
  
       !> 1.2 Add nutrients from uptake calculations
       !! 
       !! 
       veg_growth%labile%nitrogen(:)   = veg_growth%labile%nitrogen(:) + & 
                                          ( uptake_nh4(:) + uptake_no3(:) + n_fixation(:)) * dtime / 1000000.0_wp
       veg_growth%labile%nitrogen15(:) = veg_growth%labile%nitrogen15(:) + & 
                                          ( uptake_nh4_n15(:) + uptake_no3_n15(:) + n_fixation_n15(:)) * dtime / 1000000.0_wp
       veg_growth%labile%phosphorus(:) = veg_growth%labile%phosphorus(:) + & 
                                            uptake_po4(:) * dtime / 1000000.0_wp

       !> 1.2a Add mycorrhiza export
       !!
       veg_growth%labile%carbon(:)     = veg_growth%labile%carbon(:)   + & 
                                        SUM(sb_mycorrhiza_export%carbon(:,:) * soil_depth_sl(:,:),DIM=2)
       veg_growth%labile%carbon13(:)   = veg_growth%labile%carbon13(:) + & 
                                        SUM(sb_mycorrhiza_export%carbon13(:,:) * soil_depth_sl(:,:),DIM=2)
       veg_growth%labile%carbon14(:)   = veg_growth%labile%carbon14(:) + & 
                                        SUM(sb_mycorrhiza_export%carbon14(:,:) * soil_depth_sl(:,:),DIM=2)
       veg_growth%labile%nitrogen(:)   = veg_growth%labile%nitrogen(:) + & 
                                        SUM(sb_mycorrhiza_export%nitrogen(:,:) * soil_depth_sl(:,:),DIM=2)
       veg_growth%labile%nitrogen15(:) = veg_growth%labile%nitrogen15(:) + & 
                                        SUM(sb_mycorrhiza_export%nitrogen15(:,:) * soil_depth_sl(:,:),DIM=2)
       veg_growth%labile%phosphorus(:) = veg_growth%labile%phosphorus(:) + & 
                                        SUM(sb_mycorrhiza_export%phosphorus(:,:) * soil_depth_sl(:,:),DIM=2)



       !> 1.3 update labile pool and reserve 
       !! 
       veg_pool%labile  = veg_pool%labile  + veg_reserve_use
       veg_pool%reserve = veg_pool%reserve - veg_reserve_use
   
       !> 1.4 update fluxes from continuous nutrient turnover and resorption 
       !!
       veg_pool%labile%nitrogen(:)       = veg_pool%labile%nitrogen(:)       + recycling_leaf_n(:) + & 
                                         recycling_fine_root_n(:) + recycling_heart_wood_n(:)
       veg_pool%leaf%nitrogen(:)         = veg_pool%leaf%nitrogen(:)         - recycling_leaf_n(:)
       veg_pool%fine_root%nitrogen(:)    = veg_pool%fine_root%nitrogen(:)    - recycling_fine_root_n(:)
       veg_pool%heart_wood%nitrogen(:)   = veg_pool%heart_wood%nitrogen(:)   - recycling_heart_wood_n(:)
   
       veg_pool%labile%phosphorus(:)     = veg_pool%labile%phosphorus(:)     + recycling_leaf_p(:) + & 
                                         recycling_fine_root_p(:) + recycling_heart_wood_p(:)
       veg_pool%leaf%phosphorus(:)       = veg_pool%leaf%phosphorus(:)       - recycling_leaf_p(:)
       veg_pool%fine_root%phosphorus(:)  = veg_pool%fine_root%phosphorus(:)  - recycling_fine_root_p(:)
       veg_pool%heart_wood%phosphorus(:) = veg_pool%heart_wood%phosphorus(:) - recycling_heart_wood_p(:)
   
       veg_pool%labile%nitrogen15(:)     = veg_pool%labile%nitrogen15(:)     + recycling_leaf_n15(:) + & 
                                         recycling_fine_root_n15(:) + recycling_heart_wood_n15(:)
       veg_pool%leaf%nitrogen15(:)       = veg_pool%leaf%nitrogen15(:)       - recycling_leaf_n15(:)
       veg_pool%fine_root%nitrogen15(:)  = veg_pool%fine_root%nitrogen15(:)  - recycling_fine_root_n15(:)
       veg_pool%heart_wood%nitrogen15(:) = veg_pool%heart_wood%nitrogen15(:) - recycling_heart_wood_n15(:)

       !> 1.5 update labile pool with fluxes to mycorrhiza
       veg_pool%labile%carbon     = veg_pool%labile%carbon       - veg_exudation%carbon 
       veg_pool%labile%nitrogen   = veg_pool%labile%nitrogen     - veg_exudation%nitrogen 
       veg_pool%labile%phosphorus = veg_pool%labile%phosphorus   - veg_exudation%phosphorus 
       veg_pool%labile%carbon13   = veg_pool%labile%carbon13     - veg_exudation%carbon13 
       veg_pool%labile%carbon14   = veg_pool%labile%carbon14     - veg_exudation%carbon14 
       veg_pool%labile%nitrogen15 = veg_pool%labile%nitrogen15   - veg_exudation%nitrogen15 

   
       !> 1.5 Update plant biogeochemical pools with growth, litterfall and establishment
       !! 
       veg_pool           = veg_pool - veg_litterfall
       veg_pool           = veg_pool + veg_growth
       veg_pool%labile    = veg_pool%labile   + veg_establishment
       veg_pool%seed_bed  = veg_pool%seed_bed - veg_establishment
   
       !> 1.6 Radioactive decay of C14 in living plant tissue
       !!
       veg_pool%leaf%carbon14         = veg_pool%leaf%carbon14        * (1._wp - lambda_C14 * dtime )
       veg_pool%fine_root%carbon14    = veg_pool%fine_root%carbon14   * (1._wp - lambda_C14 * dtime )
       veg_pool%sap_wood%carbon14     = veg_pool%sap_wood%carbon14    * (1._wp - lambda_C14 * dtime )
       veg_pool%heart_wood%carbon14   = veg_pool%heart_wood%carbon14  * (1._wp - lambda_C14 * dtime )
       veg_pool%coarse_root%carbon14  = veg_pool%coarse_root%carbon14 * (1._wp - lambda_C14 * dtime )
       veg_pool%fruit%carbon14        = veg_pool%fruit%carbon14       * (1._wp - lambda_C14 * dtime )
       veg_pool%labile%carbon14       = veg_pool%labile%carbon14      * (1._wp - lambda_C14 * dtime )
       veg_pool%reserve%carbon14      = veg_pool%reserve%carbon14     * (1._wp - lambda_C14 * dtime )
       veg_pool%seed_bed%carbon14     = veg_pool%seed_bed%carbon14    * (1._wp - lambda_C14 * dtime )
   
       !> 2.0 Update plant diagnostics given the above fluxes
       !!
       !! 
       lai(:)      = veg_pool%leaf%carbon(:) * lctlib%sla
       dens_ind(:) = dens_ind(:) + delta_dens_ind(:)
       WHERE(dens_ind(:) > 0.0_wp .AND. lctlib%growthform == itree)
          diameter(:)  = calc_diameter_from_woody_biomass(lctlib                        , &
                                                          veg_pool%sap_wood%carbon(:)   , &
                                                          veg_pool%heart_wood%carbon(:) , &
                                                          dens_ind(:))
          height(:)    = calc_height_from_diameter(lctlib, diameter(:)) 
       END WHERE

       WHERE(veg_pool%leaf%carbon(:) > eps8)
           mean_leaf_age(:) = (mean_leaf_age(:) + dtime/one_day) * (1._wp - veg_growth%leaf%carbon(:) / veg_pool%leaf%carbon(:))
         ELSE WHERE
           mean_leaf_age(:) = 0.0_wp
       END WHERE

       dphi(:) = w_soil_root_pot(:) - lctlib%phi_leaf_min - grav * rhoh2o * height(:) * 1.e-6_wp

    END SELECT

    !> 3.0 Given current vegetation pools, update canopy layers 
    !!
    !!
    CALL calc_canopy_layers( &
                             nc                                             , & ! in
                             lctlib                                         , &
                             dsl4jsb_Config(ASSIMI_)%ncanopy                , &
                             dsl4jsb_Config(ASSIMI_)%flag_optimal_Nfraction , &
                             dsl4jsb_Config(ASSIMI_)%canopy_layer_depth(:)  , &
                             veg_pool%leaf%nitrogen(:)          , &
                             lai(:)                             , &
                             ppfd_sunlit_tfrac_mavg_cl(:,:)     , &
                             ppfd_shaded_tfrac_mavg_cl(:,:)     , & 
                             fleaf_sunlit_tfrac_mavg_cl(:,:)    , &
                             t_air_tfrac_mavg(:)                , &
                             t_air_tacclim_mavg(:)              , &
                             press_srf_tfrac_mavg(:)            , &
                             co2_mixing_ratio_tfrac_mavg(:)     , &
                             ga_tfrac_mavg(:)                   , &
                             beta_air_tfrac_mavg(:)             , &
                             beta_soil_ps_tfrac_mavg(:)         , &
                             beta_sinklim_ps_tfrac_mavg(:)      , &
                             beta_soil_gs_tfrac_mavg(:)         , & 
                             t_jmax_opt_mavg(:)                 , &
                             lai_cl(:,:)                        , & ! out
                             cumm_lai_cl(:,:)                   , & ! out
                             leaf_nitrogen_cl(:,:)              , & ! inout
                             fn_rub_cl(:,:)                     , &
                             fn_et_cl(:,:)                      , &
                             fn_pepc_cl(:,:)                    , &
                             fn_chl_cl(:,:)                     , &
                             fn_oth_cl(:,:)                     )

    !> 4.0 analysis of mass conservation: sum elements across pools into veg\%total_biomass
    !!
    !!
    veg_total_biomass =   veg_pool%leaf        + &     !< elements of leaf
                          veg_pool%fine_root   + &     !< elements of fine_root
                          veg_pool%coarse_root + &     !< elements of coarse_root
                          veg_pool%sap_wood    + &     !< elements of sap_wood
                          veg_pool%heart_wood  + &     !< elements of heart_wood
                          veg_pool%labile      + &     !< elements of labile
                          veg_pool%reserve     + &     !< elements of reserve
                          veg_pool%fruit       + &     !< elements of fruit
                          veg_pool%seed_bed            !< elements of seed_bed

! Pass back values to process memory & deallocate pool structures 
#define Dealloc_veg_pool
#define Dealloc_veg_litterfall
#define Dealloc_veg_growth
#define Dealloc_veg_establishment
#define Dealloc_veg_exudation
#define Dealloc_veg_total_biomass
#define Dealloc_veg_reserve_use
#define Dealloc_sb_mycorrhiza_export
#include "slm_pools_a_incl.f90"
#undef Dealloc_veg_pool
#undef Dealloc_veg_litterfall
#undef Dealloc_veg_growth
#undef Dealloc_veg_establishment
#undef Dealloc_veg_exudation
#undef Dealloc_veg_total_biomass
#undef Dealloc_veg_reserve_use
#undef Dealloc_sb_mycorrhiza_export

  END SUBROUTINE update_veg_pools

END MODULE mo_veg_update_pools
