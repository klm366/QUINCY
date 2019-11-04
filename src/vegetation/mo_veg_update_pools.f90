!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
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

  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message

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

    USE mo_kind,                            ONLY: wp
    USE mo_jsb_impl_constants,              ONLY: test_false_true
    USE mo_jsb_class,                       ONLY: Get_model
    USE mo_jsb_tile_class,                  ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,                  ONLY: t_jsb_task_options
    USE mo_jsb_lctlib_class,                ONLY: t_lctlib_element
    USE mo_jsb_model_class,                 ONLY: t_jsb_model
    USE mo_jsb_process_class,               ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_grid_class,                  ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                        ONLY: Get_vgrid
    USE mo_jsb_phys_constants,              ONLY: lambda_C14, rhoh2o, grav
    USE mo_jsb_math_constants,              ONLY: dtime, one_day, eps4, eps8
    USE mo_isotope_util,                    ONLY: calc_mixing_ratio_N15N14, calc_fractionation, calc_mixing_ratio_C14C
    USE mo_rad_constants,                   ONLY: rfr_ratio_toc 
    USE mo_veg_canopy,                      ONLY: calc_canopy_layers
    USE mo_veg_growth,                      ONLY: calc_diameter_from_woody_biomass, calc_height_from_diameter
    USE mo_veg_constants,                   ONLY: itree,eta_nfixation, max_leaf_shedding_rate, k_leafon_canopy, fresp_growth
    USE mo_assimi_constants,                ONLY: CiCa_default_C3, CiCa_default_C4, ic3phot, ic4phot
    USE mo_atmland_constants,               ONLY: def_co2_mixing_ratio, def_co2_mixing_ratio_C13, def_co2_mixing_ratio_C14, &
                                                  def_co2_deltaC13, def_co2_deltaC14
    USE mo_assimi_update_canopy_fluxes,     ONLY: discrimination_ps

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(ASSIMI_)    
    
    ! Use of process memories
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,    ONLY: t_pheno_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4slm_Use_memory_a_veg
    USE mo_veg_memory_class,    ONLY: veg_pool, veg_growth, veg_litterfall, &
                                      veg_exudation, veg_establishment, veg_reserve_use, veg_total_biomass
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4slm_Use_memory_a_sb
    USE mo_sb_memory_class,     ONLY: sb_mycorrhiza_export
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    dsl4jsb_Use_memory(SSE_)       !  USE mo_sse_memory_class,     ONLY: t_sse_memory


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),        POINTER               :: model
    TYPE(t_lctlib_element),   POINTER               :: lctlib           !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),        POINTER               :: vgrid_canopy     !< Vertical grid
    TYPE(t_jsb_vgrid),        POINTER               :: vgrid_soil_sb    !< Vertical grid
    INTEGER                                         :: ncanopy          !< number of canopy layers
    INTEGER                                         :: nsoil_sb         !< number of soil layers as used/defined by the SB_ process
    REAL(wp),                 DIMENSION(options%nc) :: hlp1             !< dummy variable
    REAL(wp), ALLOCATABLE                           :: canopy_layer_depth(:)  !< thickness of each single canopy layer [m2 m-2 LAI]
    INTEGER                                         :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_pools'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_memory(ASSIMI_)                       ! TYPE(t_assimi_memory),      POINTER :: assimi__mem
    dsl4jsb_Def_memory(PHENO_)                        ! TYPE(t_pheno_memory),       POINTER :: pheno__mem
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),         POINTER :: veg__mem
    dsl4jsb_Def_memory(SB_)                           ! TYPE(t_sb_memory),          POINTER :: sb__mem
    dsl4jsb_Def_memory(RAD_)                          ! TYPE(t_rad_memory),         POINTER :: rad__mem
    dsl4jsb_Def_memory(SSE_)                         ! TYPE(t_sse_memory),        POINTER :: sse__mem
    ! Declare pointers to variables in memory
    ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: gross_assimilation
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C13
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C14
    dsl4jsb_Real2D_onChunk      :: beta_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soa_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tfrac_mavg
    ! PHENO_
    dsl4jsb_Real2D_onChunk      :: growing_season 
    dsl4jsb_Real2D_onChunk      :: lai_max
    ! SSE_ 2D
    dsl4jsb_Real2D_onChunk      :: w_soil_root_pot
    ! SSE_ 3D
    dsl4jsb_Real3D_onChunk      :: soil_depth_sl
    ! RAD_ 2D
    dsl4jsb_Real2D_onChunk      :: rfr_ratio_boc_tvegdyn_mavg
    ! RAD_ 3D
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tfrac_mavg_cl
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: do_cohort_harvest
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
    dsl4jsb_Real2D_onChunk      :: net_growth
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
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: lai_cl
    dsl4jsb_Real3D_onChunk      :: cumm_lai_cl
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl
    ! Get local variables from options argument
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(VEG_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts%lib_id)
    vgrid_canopy  => Get_vgrid('canopy_layer')
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')    
    ncanopy       =  vgrid_canopy%n_levels
    nsoil_sb      =  vgrid_soil_sb%n_levels           ! needed for sb_mycorrhiza_export
    ALLOCATE(canopy_layer_depth(ncanopy))
    canopy_layer_depth(:) = vgrid_canopy%dz(:)
    ! Get process config
    dsl4jsb_Get_config(ASSIMI_)
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SSE_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, gross_assimilation)          ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, gross_assimilation_C13)      ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, gross_assimilation_C14)      ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_tfrac_mavg)         ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soa_tphen_mavg)         ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_tfrac_mavg)     ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tfrac_mavg)     ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(PHENO_, growing_season)               ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_, lai_max)                      ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SSE_, w_soil_root_pot)                ! in
    dsl4jsb_Get_var3D_onChunk(SSE_, soil_depth_sl)                  ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(RAD_, rfr_ratio_boc_tvegdyn_mavg)     ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tfrac_mavg_cl)      ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tfrac_mavg_cl)      ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, do_cohort_harvest)              ! in
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
    dsl4jsb_Get_var2D_onChunk(VEG_, net_growth)                     ! out
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
    SELECT CASE(TRIM(model%config%slm_model_name))
    CASE('plant','land')
   
       !> 1.1 net balance of the labile carbon pool due to photosynthesis and respiration
       !! 
       !! 
       npp(:)        = gross_assimilation(:)     - & 
                       growth_respiration(:)     - maint_respiration(:)     - n_processing_respiration(:) 
       npp_c13(:)    = gross_assimilation_C13(:) - &
                       growth_respiration_c13(:) - maint_respiration_c13(:) - n_processing_respiration_c13(:) 
       npp_c14(:)    = gross_assimilation_C14(:) - & 
                       growth_respiration_c14(:) - maint_respiration_c14(:) - n_processing_respiration_c14(:)
       net_growth(:) = net_growth(:) + npp(:)

       veg_growth%labile%carbon(:)   = veg_growth%labile%carbon(:)   + npp(:)     * dtime / 1000000.0_wp
       veg_growth%labile%carbon13(:) = veg_growth%labile%carbon13(:) + npp_c13(:) * dtime / 1000000.0_wp
       veg_growth%labile%carbon14(:) = veg_growth%labile%carbon14(:) + npp_c14(:) * dtime / 1000000.0_wp
  
       !> 1.2 Add nutrients from uptake calculations
       !! 
       !! 
       n_fixation_n15(:)               = n_fixation(:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(-eta_nfixation))
       veg_growth%labile%nitrogen(:)   = veg_growth%labile%nitrogen(:) + & 
                                          ( uptake_nh4(:) + uptake_no3(:) + n_fixation(:)) * dtime / 1000000.0_wp
       veg_growth%labile%nitrogen15(:) = veg_growth%labile%nitrogen15(:) + & 
                                          ( uptake_nh4_n15(:) + uptake_no3_n15(:) + n_fixation_n15(:)) * dtime / 1000000.0_wp
       veg_growth%labile%phosphorus(:) = veg_growth%labile%phosphorus(:) + & 
                                            uptake_po4(:) * dtime / 1000000.0_wp

       !> 1.2a Add mycorrhiza export
       !!
       veg_growth%labile%carbon(:)     = veg_growth%labile%carbon(:)   + & 
                                         SUM(sb_mycorrhiza_export%carbon(:,:)     * soil_depth_sl(:,:),DIM=2)
       veg_growth%labile%carbon13(:)   = veg_growth%labile%carbon13(:) + & 
                                         SUM(sb_mycorrhiza_export%carbon13(:,:)   * soil_depth_sl(:,:),DIM=2)
       veg_growth%labile%carbon14(:)   = veg_growth%labile%carbon14(:) + & 
                                         SUM(sb_mycorrhiza_export%carbon14(:,:)   * soil_depth_sl(:,:),DIM=2)
       veg_growth%labile%nitrogen(:)   = veg_growth%labile%nitrogen(:) + & 
                                         SUM(sb_mycorrhiza_export%nitrogen(:,:)   * soil_depth_sl(:,:),DIM=2)
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
                                           recycling_fine_root_n(:)          + recycling_heart_wood_n(:)
       veg_pool%leaf%nitrogen(:)         = veg_pool%leaf%nitrogen(:)         - recycling_leaf_n(:)
       veg_pool%fine_root%nitrogen(:)    = veg_pool%fine_root%nitrogen(:)    - recycling_fine_root_n(:)
       veg_pool%heart_wood%nitrogen(:)   = veg_pool%heart_wood%nitrogen(:)   - recycling_heart_wood_n(:)
   
       veg_pool%labile%phosphorus(:)     = veg_pool%labile%phosphorus(:)     + recycling_leaf_p(:) + & 
                                           recycling_fine_root_p(:)          + recycling_heart_wood_p(:)
       veg_pool%leaf%phosphorus(:)       = veg_pool%leaf%phosphorus(:)       - recycling_leaf_p(:)
       veg_pool%fine_root%phosphorus(:)  = veg_pool%fine_root%phosphorus(:)  - recycling_fine_root_p(:)
       veg_pool%heart_wood%phosphorus(:) = veg_pool%heart_wood%phosphorus(:) - recycling_heart_wood_p(:)
   
       veg_pool%labile%nitrogen15(:)     = veg_pool%labile%nitrogen15(:)     + recycling_leaf_n15(:) + & 
                                           recycling_fine_root_n15(:)        + recycling_heart_wood_n15(:)
       veg_pool%leaf%nitrogen15(:)       = veg_pool%leaf%nitrogen15(:)       - recycling_leaf_n15(:)
       veg_pool%fine_root%nitrogen15(:)  = veg_pool%fine_root%nitrogen15(:)  - recycling_fine_root_n15(:)
       veg_pool%heart_wood%nitrogen15(:) = veg_pool%heart_wood%nitrogen15(:) - recycling_heart_wood_n15(:)

       !> 1.5 update labile pool with fluxes to mycorrhiza
       veg_pool%labile%carbon(:)     = veg_pool%labile%carbon(:)       - veg_exudation%carbon(:) 
       veg_pool%labile%nitrogen(:)   = veg_pool%labile%nitrogen(:)     - veg_exudation%nitrogen(:) 
       veg_pool%labile%phosphorus(:) = veg_pool%labile%phosphorus(:)   - veg_exudation%phosphorus(:) 
       veg_pool%labile%carbon13(:)   = veg_pool%labile%carbon13(:)     - veg_exudation%carbon13(:) 
       veg_pool%labile%carbon14(:)   = veg_pool%labile%carbon14(:)     - veg_exudation%carbon14(:) 
       veg_pool%labile%nitrogen15(:) = veg_pool%labile%nitrogen15(:)   - veg_exudation%nitrogen15(:) 

   
       !> 1.5 Update plant biogeochemical pools with growth, litterfall and establishment
       !! 
       veg_pool           = veg_pool          - veg_litterfall
       veg_pool           = veg_pool          + veg_growth
       veg_pool%labile    = veg_pool%labile   + veg_establishment
       veg_pool%seed_bed  = veg_pool%seed_bed - veg_establishment
   
       !> 1.6 Radioactive decay of C14 in living plant tissue
       !!
       veg_pool%leaf%carbon14(:)         = veg_pool%leaf%carbon14(:)        * (1._wp - lambda_C14 * dtime )
       veg_pool%fine_root%carbon14(:)    = veg_pool%fine_root%carbon14(:)   * (1._wp - lambda_C14 * dtime )
       veg_pool%sap_wood%carbon14(:)     = veg_pool%sap_wood%carbon14(:)    * (1._wp - lambda_C14 * dtime )
       veg_pool%heart_wood%carbon14(:)   = veg_pool%heart_wood%carbon14(:)  * (1._wp - lambda_C14 * dtime )
       veg_pool%coarse_root%carbon14(:)  = veg_pool%coarse_root%carbon14(:) * (1._wp - lambda_C14 * dtime )
       veg_pool%fruit%carbon14(:)        = veg_pool%fruit%carbon14(:)       * (1._wp - lambda_C14 * dtime )
       veg_pool%labile%carbon14(:)       = veg_pool%labile%carbon14(:)      * (1._wp - lambda_C14 * dtime )
       veg_pool%reserve%carbon14(:)      = veg_pool%reserve%carbon14(:)     * (1._wp - lambda_C14 * dtime )
       veg_pool%seed_bed%carbon14(:)     = veg_pool%seed_bed%carbon14(:)    * (1._wp - lambda_C14 * dtime )
   
    CASE ('canopy')
   
       ! diagnose required change in leaf C given phenology to attain prescribed LAI max
       WHERE(growing_season(:) > test_false_true)
          veg_growth%leaf%carbon(:) = lai_max(:) / lctlib%sla * &
                                      (1._wp - exp(-k_leafon_canopy * &
                                       (lai_max(:) - lai(:)) / lai_max(:))) * &
                                      dtime / one_day 
       ELSEWHERE
          veg_litterfall%leaf%carbon(:) = MIN(veg_pool%leaf%carbon(:), & 
                                              lai(:)/lctlib%sla * max_leaf_shedding_rate * dtime / one_day * & 
                                              lai_max(:) / &
                                              MAX(eps4, lai(:)))
       ENDWHERE
 
       ! update leaf C:N:P pools 
       veg_pool%leaf%carbon(:)     = veg_pool%leaf%carbon(:)   + veg_growth%leaf%carbon(:) - veg_litterfall%leaf%carbon(:) 
       veg_pool%leaf%nitrogen(:)   = veg_pool%leaf%carbon(:)   / lctlib%cn_leaf 
       veg_pool%leaf%phosphorus(:) = veg_pool%leaf%nitrogen(:) / lctlib%np_leaf 

    END SELECT


    !> 2.0 Update plant diagnostics given the above fluxes
    !!
    !! 
    lai(:)      = veg_pool%leaf%carbon(:) * lctlib%sla
    dens_ind(:) = dens_ind(:) + delta_dens_ind(:)
    WHERE(dens_ind(:) > 0.0_wp .AND. lctlib%growthform == itree)
       diameter(:)  = calc_diameter_from_woody_biomass(lctlib%allom_k1               , &
                                                       lctlib%allom_k2               , &
                                                       lctlib%wood_density           , &
                                                       veg_pool%sap_wood%carbon(:)   , &
                                                       veg_pool%heart_wood%carbon(:) , &
                                                       dens_ind(:))
       height(:)    = calc_height_from_diameter(lctlib%allom_k1, lctlib%allom_k2, diameter(:)) 
    END WHERE

    WHERE(veg_pool%leaf%carbon(:) > eps8)
        mean_leaf_age(:) = (mean_leaf_age(:) + dtime / one_day) * (1._wp - veg_growth%leaf%carbon(:) / veg_pool%leaf%carbon(:))
      ELSEWHERE
        mean_leaf_age(:) = 0.0_wp
    END WHERE

    dphi(:) = w_soil_root_pot(:) - lctlib%phi_leaf_min - grav * rhoh2o * height(:) * 1.e-6_wp


    !> 3.0 stand/cohort harvest
    !!
    !! in case of complete harvest, set all pools to zero
    !! (to account for wood extraction, which would otherwise leave some wood, and reset the key plant diagnostics)
    IF(model%config%l_do_stand_replacing_harvest .OR. do_cohort_harvest(ics) > test_false_true) THEN 
      !! differ between cohort & population scheme
      IF (do_cohort_harvest(ics) > test_false_true) THEN    ! -> cohort scheme
        dens_ind(:)                    = 1.0_wp
        veg_pool%reserve%carbon(:)     = 1.0_wp
      ELSE ! model%config%l_do_stand_replacing_harvest -> population scheme
        dens_ind(:)                    = 0.1_wp
        veg_pool%reserve%carbon(:)     = 0.1_wp
      ENDIF
      !! zero
      veg_pool%leaf                  = veg_pool%leaf        * 0.0_wp
      veg_pool%fine_root             = veg_pool%fine_root   * 0.0_wp
      veg_pool%coarse_root           = veg_pool%coarse_root * 0.0_wp
      veg_pool%sap_wood              = veg_pool%sap_wood    * 0.0_wp
      veg_pool%heart_wood            = veg_pool%heart_wood  * 0.0_wp
      veg_pool%fruit                 = veg_pool%fruit       * 0.0_wp
      veg_pool%labile                = veg_pool%labile      * 0.0_wp
      height(:)                      = 0.0_wp
      diameter(:)                    = 0.0_wp
      lai(:)                         = 0.0_wp
      mean_leaf_age(:)               = 0.0_wp
      !! init
      dphi(:)                        = w_soil_root_pot(:) - lctlib%phi_leaf_min 
      rfr_ratio_boc_tvegdyn_mavg(:)  = rfr_ratio_toc 
      ! pool%reserve% C N P
      veg_pool%reserve%nitrogen(:)   = veg_pool%reserve%carbon(:)    / (1._wp + fresp_growth) / lctlib%cn_leaf
      veg_pool%reserve%phosphorus(:) = veg_pool%reserve%nitrogen(:)  / lctlib%np_leaf
      ! pool%reserve%carbon13
      IF(lctlib%ps_pathway == ic3phot)THEN
         hlp1(:) = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
      ELSE
         hlp1(:) = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
      ENDIF
      veg_pool%reserve%carbon13(:)   = calc_fractionation(def_co2_mixing_ratio, def_co2_mixing_ratio_C13, hlp1(:)) * &
                                       veg_pool%reserve%carbon(:)
      ! pool%reserve%carbon14
      IF(lctlib%ps_pathway == ic3phot)THEN
         hlp1(:) = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
      ELSE
         hlp1(:) = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
      ENDIF
      veg_pool%reserve%carbon14(:)   = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14-hlp1(:)) * & 
                                       veg_pool%reserve%carbon(:)
      ! pool%reserve%nitrogen15
      veg_pool%reserve%nitrogen15(:) = veg_pool%reserve%nitrogen(:) / (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    ENDIF


    !> 4.0 Given current vegetation pools, update canopy layers 
    !!
    !!
    CALL calc_canopy_layers( nc                                             , & ! in
                             ncanopy                                        , &
                             vgrid_canopy%dz(:)                             , &
                             vgrid_canopy%lbounds(:)                        , &
                             vgrid_canopy%ubounds(:)                        , &
                             lctlib%ps_pathway                              , &
                             lctlib%k0_fn_struc                             , &
                             lctlib%fn_oth_min                              , &
                             lctlib%sla                                     , &
                             lctlib%np_leaf                                 , &
                             lctlib%gmin                                    , &
                             lctlib%g0                                      , &
                             lctlib%g1                                      , &
                             lctlib%t_jmax_omega                            , &
                             dsl4jsb_Config(ASSIMI_)%flag_optimal_Nfraction , &
                             canopy_layer_depth(:)                          , &
                             veg_pool%leaf%nitrogen(:)                      , &
                             lai(:)                                         , &
                             ppfd_sunlit_tfrac_mavg_cl(:,:)                 , &
                             ppfd_shaded_tfrac_mavg_cl(:,:)                 , & 
                             fleaf_sunlit_tfrac_mavg_cl(:,:)                , &
                             t_air_tfrac_mavg(:)                            , &
                             t_air_tacclim_mavg(:)                          , &
                             press_srf_tfrac_mavg(:)                        , &
                             co2_mixing_ratio_tfrac_mavg(:)                 , &
                             ga_tfrac_mavg(:)                               , &
                             beta_air_tfrac_mavg(:)                         , &
                             beta_soa_tphen_mavg(:)                         , &
                             beta_soil_ps_tfrac_mavg(:)                     , &
                             beta_sinklim_ps_tfrac_mavg(:)                  , &
                             beta_soil_gs_tfrac_mavg(:)                     , & 
                             t_jmax_opt_mavg(:)                             , &
                             lai_cl(:,:)                                    , & ! out
                             cumm_lai_cl(:,:)                               , & ! out
                             leaf_nitrogen_cl(:,:)                          , & ! inout
                             fn_rub_cl(:,:)                                 , &
                             fn_et_cl(:,:)                                  , &
                             fn_pepc_cl(:,:)                                , &
                             fn_chl_cl(:,:)                                 , &
                             fn_oth_cl(:,:)                                 )


    !> 5.0 analysis of mass conservation: sum elements across pools into veg\%total_biomass
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

! Pass back values to process memory
#define Passback_veg_pool
#define Passback_veg_litterfall
#define Passback_veg_growth
#define Passback_veg_establishment
#define Passback_veg_exudation
#define Passback_veg_total_biomass
#define Passback_veg_reserve_use
#define Passback_sb_mycorrhiza_export
#include "slm_pools_a_incl.f90"
#undef Passback_veg_pool
#undef Passback_veg_litterfall
#undef Passback_veg_growth
#undef Passback_veg_establishment
#undef Passback_veg_exudation
#undef Passback_veg_total_biomass
#undef Passback_veg_reserve_use
#undef Passback_sb_mycorrhiza_export

  END SUBROUTINE update_veg_pools

END MODULE mo_veg_update_pools
