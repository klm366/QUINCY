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
!! @brief Initialization of the the vegetation memory.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20161213)
!!
MODULE mo_veg_init
#ifndef __NO_JSBACH__

#ifdef __QUINCY_NO_JSB__
#include "dsl_jsbach.inc"
#endif

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message, finish, message_text
  USE mo_jsb_control,             ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: veg_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize vegetation process (after memory has been set up)
  !!
  !! differtiates between the QUINCY models: land, plant, soil, canopy
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE veg_init(tile)

    USE mo_jsb_class,                       ONLY: Get_model
    USE mo_jsb_tile_class,                  ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,                  ONLY: t_jsb_task_options
    USE mo_jsb_model_class,                 ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,                ONLY: t_lctlib_element
    USE mo_jsb_process_class,               ONLY: VEG_, PHENO_

    USE mo_jsb_math_constants,              ONLY: eps8
    USE mo_assimi_constants,                ONLY: CiCa_default_C3, CiCa_default_C4, ic3phot, ic4phot
    USE mo_atmland_constants,               ONLY: def_co2_mixing_ratio, def_co2_mixing_ratio_C13, def_co2_mixing_ratio_C14, &
                                                  standard_press_srf, def_co2_deltaC13, def_co2_deltaC14
    USE mo_jsb_physical_constants,          ONLY: molar_mass_N
    USE mo_veg_constants,                   ONLY: fresp_growth, background_mort_rate_tree, background_mort_rate_grass, w_root_lim_max
    USE mo_isotope_util,                    ONLY: calc_fractionation,calc_mixing_ratio_N15N14,calc_mixing_ratio_C14C
    USE mo_assimi_update_canopy_fluxes,     ONLY: discrimination_ps
    USE mo_veg_constants,                   ONLY: root2leaf_cn, root2leaf_np, itree

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(VEG_)
    ! Use of process memories
    dsl4jsb_Use_memory(PHENO_)
    dsl4jsb_Use_memory(VEG_)                ! USE mo_veg_memory_class,      ONLY: t_veg_memory
    ! Use of pools
#ifdef __QUINCY_NO_JSB__
    dsl4slm_Use_memory_a_veg                ! USE mo_veg_memory_class,      ONLY: t_veg_element_a, t_veg_pool_a
    USE mo_veg_memory_class,                ONLY: veg_pool
#else
    dsl4jsb_Use_pools(VEG_) t_veg_tissues_pool
#endif

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp), ALLOCATABLE                        :: hlp1(:,:)  !< dummy variable
    REAL(wp)                                     :: ref_lai    !< LAI                                 | for canopy model init
    REAL(wp)                                     :: ref_nleaf  !< N leaf per unit leaf area in gN/m2  | for canopy model init
    INTEGER                                      :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':veg_init'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration, process memory and pools
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)
#ifdef __QUINCY_NO_JSB__
#else
    dsl4jsb_Def_pool(t_veg_tissues_pool)  :: veg_pool
#endif
    ! Declare pointers to variables in memory
    ! PHENO_
    dsl4jsb_Real2D_onDomain    :: lai
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: target_cn_leaf
    dsl4jsb_Real2D_onChunk      :: target_np_leaf
    dsl4jsb_Real2D_onChunk      :: target_cn_fine_root
    dsl4jsb_Real2D_onChunk      :: target_np_fine_root
    dsl4jsb_Real2D_onChunk      :: dens_ind
    dsl4jsb_Real2D_onChunk      :: height 
    dsl4jsb_Real2D_onChunk      :: dphi
    dsl4jsb_Real2D_onChunk      :: npp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: n_fixation_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: transpiration_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_npp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_transpiration_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_p_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: dphi_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cn_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_np_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_carbon_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_carbon_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: press_srf_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: ga_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: press_srf_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: uptake_n_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cn_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: growth_np_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: npp_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps
    dsl4jsb_Real2D_onChunk      :: leaf_cn_direction
    dsl4jsb_Real2D_onChunk      :: t_air_daytime
    dsl4jsb_Real2D_onChunk      :: press_srf_daytime
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_daytime
    dsl4jsb_Real2D_onChunk      :: ga_daytime
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_daytime
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_daytime
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: t_air_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: press_srf_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: ga_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_p_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: npp_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: demand_uptake_n_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: demand_uptake_p_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_p_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: cost_n_fixation_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: ga_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: leaf2root_troot_mavg
    dsl4jsb_Real2D_onChunk      :: an_boc_tvegdyn_mavg
    dsl4jsb_Real2D_onChunk      :: net_growth_tvegdyn_mavg
    dsl4jsb_Real2D_onChunk      :: lai_tvegdyn_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tacclim_mavg
    dsl4jsb_Real2D_onChunk      :: t_soil_root_tacclim_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_week_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_month_mavg
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    dsl4jsb_Real2D_onChunk      :: w_root_lim_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: exudation_c_tmyc_mavg
    ! VEG_ 2D pool helper variables
    dsl4jsb_Real2D_onDomain     :: veg_pool_leaf_carbon
    dsl4jsb_Real2D_onDomain     :: veg_growth_leaf_carbon
    dsl4jsb_Real2D_onDomain     :: veg_litterfall_leaf_carbon
    dsl4jsb_Real2D_onDomain     :: veg_pool_leaf_nitrogen
    dsl4jsb_Real2D_onDomain     :: veg_pool_leaf_phosphorus
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tcnl_mavg_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl
    dsl4jsb_Real3D_onChunk      :: lai_cl
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(VEG_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    IF (tile%lcts(1)%lib_id /= 0) THEN            ! work with lctlib only if the present tile is a pft
      lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ELSE
      CALL finish(TRIM(routine), ' Cannot be executed on a non-PFT tile, as the lctlib is needed.')
    ENDIF
    ! Get process config
    dsl4jsb_Get_config(VEG_)
    ! Get process memories
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(VEG_)
#ifdef __QUINCY_NO_JSB__
#else
    ! Get pools
    dsl4jsb_Get_pool(VEG_, veg_main_pool%veg_pool,          veg_pool)
#endif
    ! Get process variables (Set pointers to variables in memory)
    ! PHENO_
    dsl4jsb_Get_var2D_onDomain(PHENO_, lai)                           ! in
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_leaf)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_leaf)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_fine_root)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_fine_root)
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)
    dsl4jsb_Get_var2D_onChunk(VEG_, height)
    dsl4jsb_Get_var2D_onChunk(VEG_, dphi)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_npp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, dphi_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cn_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_np_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_carbon_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_nitrogen_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_phosphorus_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_carbon_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_nitrogen_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_phosphorus_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_n_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cn_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_np_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps)
    dsl4jsb_Get_var2D_onChunk(VEG_, leaf_cn_direction)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tlabile_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tlabile_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tphen_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, demand_uptake_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, demand_uptake_p_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, cost_n_fixation_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, leaf2root_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, an_boc_tvegdyn_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, net_growth_tvegdyn_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, lai_tvegdyn_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tacclim_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_soil_root_tacclim_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_week_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_month_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, w_root_lim_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, exudation_c_tmyc_mavg)
    ! VEG_ 2D pool helper variables
    dsl4jsb_Get_var2D_onDomain(VEG_, veg_pool_leaf_carbon)
    dsl4jsb_Get_var2D_onDomain(VEG_, veg_growth_leaf_carbon)
    dsl4jsb_Get_var2D_onDomain(VEG_, veg_litterfall_leaf_carbon)
    dsl4jsb_Get_var2D_onDomain(VEG_, veg_pool_leaf_nitrogen)
    dsl4jsb_Get_var2D_onDomain(VEG_, veg_pool_leaf_phosphorus)
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tcnl_mavg_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_oth_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)
#ifdef __QUINCY_NO_JSB__
    ! Get pools
#define Get_veg_pool
#include "slm_pools_a_incl.f90"
#undef Get_veg_pool
#endif


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 1.0 vegetation init specific to the SLM model
    !!
    !! SLM models: land, plant,soil (case 1) and canopy (case 2)
    SELECT CASE(TRIM(model%config%slm_model_name))

      !> 1.1 Land & Plant & Soil
      !!
      CASE('land','plant','soil')
        target_cn_leaf(:)             = lctlib%cn_leaf
        target_np_leaf(:)             = lctlib%np_leaf
        target_cn_fine_root(:)        = target_cn_leaf(:) / root2leaf_cn
        target_np_fine_root(:)        = target_np_leaf(:) / root2leaf_np

        !> differ between vegetation dynamics schemes
        !! 
        SELECT CASE (TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme))
          !> start bareground with vegetation dynamics switched on: cohort scheme
          !! 
          CASE ("cohort")
            background_mort_rate_tree     = 0.0_wp    ! default: 0.01_wp
            background_mort_rate_grass    = 0.0_wp    ! default: 0.05_wp 
            dens_ind(:)                   = 1.0_wp

            ! pool%reserve% C N P
            veg_pool%reserve%carbon(:)     = 1.0_wp
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
            veg_pool%reserve%nitrogen15(:) = veg_pool%reserve%nitrogen(:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          !> start bareground with vegetation dynamics switched on: population scheme
          !! 
          CASE ("population")
            veg_pool%seed_bed%carbon(:)      = lctlib%seed_size * 10.0_wp
            dens_ind(:)                      = 0.0_wp

            veg_pool%seed_bed%nitrogen(:)    = veg_pool%seed_bed%carbon(:)    / (1._wp + fresp_growth) / lctlib%cn_leaf
            veg_pool%seed_bed%phosphorus(:)  = veg_pool%seed_bed%nitrogen(:)  / lctlib%np_leaf

            ! pool%seed_bed%carbon13
            IF(lctlib%ps_pathway == ic3phot)THEN
               hlp1(:) = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
            ELSE
               hlp1(:) = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
            ENDIF
            veg_pool%seed_bed%carbon13(:) = veg_pool%seed_bed%carbon(:) * &
                                            calc_fractionation(def_co2_mixing_ratio, def_co2_mixing_ratio_C13, hlp1(:))

            ! pool%seed_bed%carbon14
            IF(lctlib%ps_pathway == ic3phot)THEN
               hlp1(:) = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
            ELSE
               hlp1(:) = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
            ENDIF
            veg_pool%seed_bed%carbon14(:)    = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14-hlp1(:)) * & 
                                               veg_pool%seed_bed%carbon(:) 

            veg_pool%seed_bed%nitrogen15(:)  = veg_pool%seed_bed%nitrogen(:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          !> no vegetation dynamics and not starting with bare ground
          !! 
          CASE ("none")
            veg_pool%reserve%carbon(:)      = 5.0_wp
            dens_ind(:)                     = 0.4_wp

            veg_pool%reserve%nitrogen(:)    = veg_pool%reserve%carbon(:)    / (1._wp + fresp_growth) / lctlib%cn_leaf
            veg_pool%reserve%phosphorus(:)  = veg_pool%reserve%nitrogen(:)  / lctlib%np_leaf

            ! pool%reserve%carbon13
            IF(lctlib%ps_pathway == ic3phot)THEN
               hlp1(:) = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
            ELSE
               hlp1(:) = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
            ENDIF
            veg_pool%reserve%carbon13(:) = veg_pool%reserve%carbon(:) * &
                                           calc_fractionation(def_co2_mixing_ratio, def_co2_mixing_ratio_C13, hlp1(:))

            ! pool%reserve%carbon14
            IF(lctlib%ps_pathway == ic3phot)THEN
               hlp1(:) = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
            ELSE
               hlp1(:) = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
            ENDIF
            veg_pool%reserve%carbon14(:)   = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14-hlp1(:)) * & 
                                             veg_pool%reserve%carbon(:) 

            veg_pool%reserve%nitrogen15(:) = veg_pool%reserve%nitrogen(:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          CASE DEFAULT
            WRITE(message_text,'(2a)') 'not a valid veg_dynamics_scheme: ',TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme)
            CALL finish("mo_veg_config_class:init_veg_config",message_text) 
        ENDSELECT


      !> 1.2 canopy model
      !!
      CASE('canopy')
        ref_lai   = 6.0_wp  !! the ref_lai of 6.0 works great for all sites that are not too dry; but use 0.5_wp for dry sites like ZA-Kru
        ref_nleaf = 2.0_wp
        lai(:)    = ref_lai
        !leaf_nitrogen
        IF(dsl4jsb_Config(VEG_)%pft_id==0) THEN
          veg_pool%leaf%nitrogen(:)   = ref_nleaf / molar_mass_N * lai(:)
          height                      = 10.0_wp
        ELSE
          veg_pool%leaf%nitrogen(:)   = lai(:) / lctlib%sla / lctlib%cn_leaf
          veg_pool%leaf%phosphorus(:) = veg_pool%leaf%nitrogen(:) / lctlib%np_leaf
          IF(lctlib%growthform==itree)THEN
             height(:)                = 10.0_wp
          ELSE 
             height(:)                = 1.0_wp
          ENDIF
        ENDIF
        leaf_nitrogen_cl(:,:) = 0.0_wp  ! this is intentionally wrong; i.e., it will not work at the 1st timestep (which is probably night anyway)
        lai_cl(:,:)           = 0.0_wp  ! this is intentionally wrong; i.e., it will not work at the 1st timestep (which is probably night anyway)

        ! code for initial jsb4 canopy-model version
        ! 2D var used for vegetation Pool/Flux replacement
        veg_pool_leaf_carbon(:,:)           = 0.0_wp
        veg_growth_leaf_carbon(:,:)         = 0.0_wp
        veg_litterfall_leaf_carbon(:,:)     = 0.0_wp
        IF (tile%lcts(1)%lib_id /= 0) THEN      ! work with lctlib only if the present tile is a pft
          veg_pool_leaf_carbon(:,:)     = lai(:,:) / lctlib%sla
          veg_pool_leaf_nitrogen(:,:)   = lai(:,:) / lctlib%sla / lctlib%cn_leaf
          veg_pool_leaf_phosphorus(:,:) = veg_pool_leaf_nitrogen(:,:) / lctlib%np_leaf
          IF(lctlib%growthform == itree)THEN
             height(:,:)                = 10.0_wp
          ELSE 
             height(:,:)                = 1.0_wp
          ENDIF
        ENDIF

      !> 1.3 Default (simply for error catching)
      !!
      CASE DEFAULT
        CALL finish('veg_init','Invalid QUINCY model name. Use either land, plant, soil or canopy.')
    ENDSELECT


    !> 2.0 vegetation init for all SLM models
    !!
    ! plant variables
    dphi(:)                             = -1.0_wp * lctlib%phi_leaf_min
    beta_sinklim_ps(:)                  = 1.0_wp
    leaf_cn_direction(:)                = 0.0_wp
    fn_oth_cl(:,:)                      = 1.0_wp  ! it is necessary to start with fn_oth = 1.0 leaving the rest of the fraction 0.0
                                                  ! fn_chl_cl, fn_rub_cl, fn_et_cl, fn_pepc_cl are set to 0.0_wp by default with the add_var()

    ! initialising moving averages
    ! allocation
    npp_talloc_mavg(:)                  = 0.0_wp
    n_fixation_talloc_mavg(:)           = 0.0_wp
    transpiration_talloc_mavg(:)        = eps8
    unit_transpiration_talloc_mavg(:)   = 1.e-6_wp
    unit_uptake_n_talloc_mavg(:)        = eps8
    unit_uptake_p_talloc_mavg(:)        = 4._wp * 1.e-10_wp
    dphi_talloc_mavg(:)                 = -1.0_wp * lctlib%phi_leaf_min
    growth_cn_talloc_mavg(:)            = 60.0_wp
    unit_npp_talloc_mavg(:)             = unit_uptake_n_talloc_mavg(:) * growth_cn_talloc_mavg(:)
    growth_cp_talloc_mavg(:)            = 2000.0_wp
    growth_np_talloc_mavg(:)            = growth_cp_talloc_mavg(:) / growth_cn_talloc_mavg(:)
    labile_carbon_talloc_mavg(:)        = 0.0_wp
    labile_nitrogen_talloc_mavg(:)      = 0.0_wp
    labile_phosphorus_talloc_mavg(:)    = 0.0_wp
    reserve_carbon_talloc_mavg(:)       = 0.0_wp
    reserve_nitrogen_talloc_mavg(:)     = 0.0_wp
    reserve_phosphorus_talloc_mavg(:)   = 0.0_wp
    w_root_lim_talloc_mavg(:)           = w_root_lim_max
    exudation_c_tmyc_mavg(:)           = 0.0_wp

    ! N fraction optimisation
    fleaf_sunlit_tfrac_mavg_cl(:,:)     = 0.5_wp
    t_air_tfrac_mavg(:)                 = 283.15_wp
    press_srf_tfrac_mavg(:)             = standard_press_srf
    co2_mixing_ratio_tfrac_mavg(:)      = def_co2_mixing_ratio
    ga_tfrac_mavg(:)                    = 0.2_wp
    beta_sinklim_ps_tfrac_mavg(:)       = 1.0_wp

    ! leaf N optimisation
    fleaf_sunlit_tcnl_mavg_cl(:,:)      = 0.5_wp
    t_air_tcnl_mavg(:)                  = 283.15_wp
    press_srf_tcnl_mavg(:)              = standard_press_srf
    co2_mixing_ratio_tcnl_mavg(:)       = def_co2_mixing_ratio
    beta_sinklim_ps_tcnl_mavg(:)        = 1.0_wp
    growth_cn_tcnl_mavg(:)              = 60.0_wp
    growth_np_tcnl_mavg(:)              = 30.0_wp
    uptake_n_tcnl_mavg(:)               = eps8
    npp_tcnl_mavg(:)                    = uptake_n_tcnl_mavg(:) * growth_cn_tcnl_mavg(:)

    ! daytime averages plus accumulation var
    t_air_daytime(:)                    = 283.15_wp
    press_srf_daytime(:)                = standard_press_srf
    co2_mixing_ratio_daytime(:)         = def_co2_mixing_ratio
    ga_daytime(:)                       = 0.2_wp
    beta_sinklim_ps_daytime(:)          = 1.0_wp
    t_jmax_opt_daytime(:)               = 17.0_wp
    t_jmax_opt_daytime_dacc(:)          = 17.0_wp
    t_air_daytime_dacc(:)               = 283.15_wp
    press_srf_daytime_dacc(:)           = standard_press_srf
    co2_mixing_ratio_daytime_dacc(:)    = def_co2_mixing_ratio
    ga_daytime_dacc(:)                  = 0.2_wp
    beta_sinklim_ps_daytime_dacc(:)     = 1.0_wp

    ! moving averages
    growth_req_n_tlabile_mavg(:)        = 1._wp / 25._wp
    growth_req_p_tlabile_mavg(:)        = 1._wp / 14._wp
    t_air_tphen_mavg(:)                 = 283.15_wp
    npp_tuptake_mavg(:)                 = 1.0_wp
    demand_uptake_n_tuptake_mavg(:)     = 0.1_wp
    demand_uptake_p_tuptake_mavg(:)     = 0.1_wp
    growth_req_n_tuptake_mavg(:)        = 1._wp / 25._wp
    growth_req_p_tuptake_mavg(:)        = 1._wp / 14._wp
    cost_n_fixation_talloc_mavg(:)      = 10.0_wp
    ga_tcnl_mavg(:)                     = 0.2_wp
    leaf2root_troot_mavg(:)             = 1.0_wp
    an_boc_tvegdyn_mavg(:)              = 5.0_wp
    net_growth_tvegdyn_mavg(:)          = 5.0_wp
    lai_tvegdyn_mavg(:)                 = 0.1_wp
    t_air_tacclim_mavg(:)               = 283.15_wp     ! this init value may equal t_acclim_zero (283.15) to make work: mo_veg_respiration:temperature_response_respiration
    t_soil_root_tacclim_mavg(:)         = 283.15_wp     ! this init value may equal t_acclim_zero (283.15) to make work: mo_veg_respiration:temperature_response_respiration
    t_air_week_mavg(:)                  = 283.15_wp
    t_air_month_mavg(:)                 = 283.15_wp
    t_jmax_opt_mavg(:)                  = 17.0_wp

#ifdef __QUINCY_NO_JSB__
! Pass back values to process memory
#define Passback_veg_pool
#include "slm_pools_a_incl.f90"
#undef Passback_veg_pool
#endif

  END SUBROUTINE veg_init

#endif
END MODULE mo_veg_init
