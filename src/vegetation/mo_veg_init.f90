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

#ifdef __QUINCY_STANDALONE__
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

  ! ======================================================================================================= !
  !>Intialize VEG_ process memory for QUINCY
  !>
  !>  differtiates between the QUINCY models: land, plant, soil, canopy
  !>
  SUBROUTINE veg_init(tile)

    USE mo_jsb_class,                       ONLY: Get_model
    USE mo_jsb_tile_class,                  ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,                 ONLY: t_jsb_model, MODEL_JSBACH
    USE mo_jsb_lctlib_class,                ONLY: t_lctlib_element
    USE mo_jsb_grid_class,                  ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                        ONLY: Get_vgrid
    ! USE mo_jsb_process_class,               ONLY: VEG_, PHENO_, SSE_
    USE mo_jsb_process_class,               ONLY: VEG_, PHENO_, SPQ_

    USE mo_jsb_math_constants,              ONLY: eps8
    USE mo_assimi_constants,                ONLY: CiCa_default_C3, CiCa_default_C4, ic3phot, ic4phot
    USE mo_atmland_constants,               ONLY: def_co2_mixing_ratio, def_co2_mixing_ratio_C13, def_co2_mixing_ratio_C14, &
                                                  standard_press_srf, def_co2_deltaC13, def_co2_deltaC14
    USE mo_jsb_physical_constants,          ONLY: molar_mass_N
    USE mo_veg_constants,                   ONLY: fresp_growth, background_mort_rate_tree, background_mort_rate_grass, &
                                                  w_root_lim_max
    USE mo_isotope_util,                    ONLY: calc_fractionation, calc_mixing_ratio_N15N14, calc_mixing_ratio_C14C
    USE mo_assimi_update_canopy_fluxes,     ONLY: discrimination_ps
    USE mo_veg_constants,                   ONLY: root2leaf_cn, root2leaf_np, itree
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_memory(PHENO_)
    ! dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_STANDALONE__
    dsl4slm_Use_memory_a_veg                ! USE mo_veg_memory_class, ONLY: t_veg_element_a, t_veg_pool_a
    USE mo_veg_memory_class,                ONLY: veg_pool
#else
    dsl4jsb_Use_pools(VEG_) t_veg_tissues_pool
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)    :: tile         !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER              :: model             !
    TYPE(t_lctlib_element), POINTER              :: lctlib            ! land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER              :: vgrid_soil_sb     ! Vertical grid
    REAL(wp), ALLOCATABLE                        :: hlp1(:,:)         ! dummy variable
    REAL(wp)                                     :: ref_lai           ! LAI | for canopy model init
    INTEGER                                      :: isoil, nsoil_sb   ! looping
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':veg_init'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(PHENO_)
    ! dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_STANDALONE__
#else
    dsl4jsb_Def_pool(t_veg_tissues_pool)  :: veg_pool
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    ! PHENO_ 2D
    dsl4jsb_Real2D_onDomain     :: lai
    ! PHENO_ 3D
    dsl4jsb_Real3D_onDomain     :: lai_cl
    ! ! SSE_ 3D
    ! dsl4jsb_Real3D_onDomain     :: soil_lay_depth_ubound_sl
    ! SPQ_ 3D
    dsl4jsb_Real3D_onDomain     :: soil_lay_depth_ubound_sl
    ! VEG_ 2D
    dsl4jsb_Real2D_onDomain     :: target_cn_leaf
    dsl4jsb_Real2D_onDomain     :: target_np_leaf
    dsl4jsb_Real2D_onDomain     :: target_cn_fine_root
    dsl4jsb_Real2D_onDomain     :: target_np_fine_root
    dsl4jsb_Real2D_onDomain     :: dens_ind
    dsl4jsb_Real2D_onDomain     :: height 
    dsl4jsb_Real2D_onDomain     :: dphi
    dsl4jsb_Real2D_onDomain     :: npp_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: n_fixation_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: transpiration_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: unit_npp_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: unit_transpiration_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: unit_uptake_n_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: unit_uptake_p_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: dphi_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: growth_cn_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: growth_cp_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: growth_np_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: labile_carbon_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: labile_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: labile_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: reserve_carbon_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: reserve_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: reserve_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: t_air_tfrac_mavg
    dsl4jsb_Real2D_onDomain     :: press_srf_tfrac_mavg
    dsl4jsb_Real2D_onDomain     :: co2_mixing_ratio_tfrac_mavg
    dsl4jsb_Real2D_onDomain     :: ga_tfrac_mavg
    dsl4jsb_Real2D_onDomain     :: beta_sinklim_ps_tfrac_mavg
    dsl4jsb_Real2D_onDomain     :: t_air_tcnl_mavg
    dsl4jsb_Real2D_onDomain     :: press_srf_tcnl_mavg
    dsl4jsb_Real2D_onDomain     :: co2_mixing_ratio_tcnl_mavg
    dsl4jsb_Real2D_onDomain     :: beta_sinklim_ps_tcnl_mavg
    dsl4jsb_Real2D_onDomain     :: uptake_n_tcnl_mavg
    dsl4jsb_Real2D_onDomain     :: growth_cn_tcnl_mavg
    dsl4jsb_Real2D_onDomain     :: growth_np_tcnl_mavg
    dsl4jsb_Real2D_onDomain     :: npp_tcnl_mavg
    dsl4jsb_Real2D_onDomain     :: beta_sinklim_ps
    dsl4jsb_Real2D_onDomain     :: leaf_cn_direction
    dsl4jsb_Real2D_onDomain     :: t_air_daytime
    dsl4jsb_Real2D_onDomain     :: press_srf_daytime
    dsl4jsb_Real2D_onDomain     :: co2_mixing_ratio_daytime
    dsl4jsb_Real2D_onDomain     :: ga_daytime
    dsl4jsb_Real2D_onDomain     :: beta_sinklim_ps_daytime
    dsl4jsb_Real2D_onDomain     :: t_jmax_opt_daytime
    dsl4jsb_Real2D_onDomain     :: t_jmax_opt_daytime_dacc
    dsl4jsb_Real2D_onDomain     :: t_air_daytime_dacc
    dsl4jsb_Real2D_onDomain     :: press_srf_daytime_dacc
    dsl4jsb_Real2D_onDomain     :: co2_mixing_ratio_daytime_dacc
    dsl4jsb_Real2D_onDomain     :: ga_daytime_dacc
    dsl4jsb_Real2D_onDomain     :: beta_sinklim_ps_daytime_dacc
    dsl4jsb_Real2D_onDomain     :: growth_req_n_tlabile_mavg
    dsl4jsb_Real2D_onDomain     :: growth_req_p_tlabile_mavg
    dsl4jsb_Real2D_onDomain     :: t_air_tphen_mavg
    dsl4jsb_Real2D_onDomain     :: npp_tuptake_mavg
    dsl4jsb_Real2D_onDomain     :: demand_uptake_n_tuptake_mavg
    dsl4jsb_Real2D_onDomain     :: demand_uptake_p_tuptake_mavg
    dsl4jsb_Real2D_onDomain     :: growth_req_n_tuptake_mavg
    dsl4jsb_Real2D_onDomain     :: growth_req_p_tuptake_mavg
    dsl4jsb_Real2D_onDomain     :: cost_n_fixation_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: ga_tcnl_mavg
    dsl4jsb_Real2D_onDomain     :: leaf2root_troot_mavg
    dsl4jsb_Real2D_onDomain     :: an_boc_tvegdyn_mavg
    dsl4jsb_Real2D_onDomain     :: net_growth_tvegdyn_mavg
    dsl4jsb_Real2D_onDomain     :: lai_tvegdyn_mavg
    dsl4jsb_Real2D_onDomain     :: t_air_tacclim_mavg
    dsl4jsb_Real2D_onDomain     :: t_soil_root_tacclim_mavg
    dsl4jsb_Real2D_onDomain     :: t_air_week_mavg
    dsl4jsb_Real2D_onDomain     :: t_air_month_mavg
    dsl4jsb_Real2D_onDomain     :: t_jmax_opt_mavg
    dsl4jsb_Real2D_onDomain     :: w_root_lim_talloc_mavg
    dsl4jsb_Real2D_onDomain     :: exudation_c_tmyc_mavg
    ! VEG_ 3D
    dsl4jsb_Real3D_onDomain     :: fleaf_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onDomain     :: fleaf_sunlit_tcnl_mavg_cl
    dsl4jsb_Real3D_onDomain     :: fn_oth_cl
    dsl4jsb_Real3D_onDomain     :: leaf_nitrogen_cl
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(VEG_)) RETURN
    IF (debug_on()) CALL message(TRIM(routine), 'Initializing QUINCY VEG_ memory for tile '//TRIM(tile%name))
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    IF (model%config%model_scheme == MODEL_JSBACH) RETURN     ! only with QUINCY model
    IF (tile%lcts(1)%lib_id == 0) RETURN                      ! run this init only if the present tile is a pft
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_memory(PHENO_)
    ! dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_STANDALONE__
#else
    ! Get pools
    dsl4jsb_Get_pool(VEG_, veg_main_pool%veg_pool,          veg_pool)
    ! associate domain variables
    CALL veg_pool%Associate_domain_var_ptr()
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    ! PHENO_ 2D
    dsl4jsb_Get_var2D_onDomain(PHENO_, lai)                           ! in
    ! PHENO_ 3D
    dsl4jsb_Get_var3D_onDomain(PHENO_, lai_cl)                        ! in
    ! ! SSE_ 3D
    ! dsl4jsb_Get_var3D_onDomain(SSE_,   soil_lay_depth_ubound_sl)      ! in
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onDomain(SPQ_,   soil_lay_depth_ubound_sl)      ! in
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onDomain(VEG_, target_cn_leaf)
    dsl4jsb_Get_var2D_onDomain(VEG_, target_np_leaf)
    dsl4jsb_Get_var2D_onDomain(VEG_, target_cn_fine_root)
    dsl4jsb_Get_var2D_onDomain(VEG_, target_np_fine_root)
    dsl4jsb_Get_var2D_onDomain(VEG_, dens_ind)
    dsl4jsb_Get_var2D_onDomain(VEG_, height)
    dsl4jsb_Get_var2D_onDomain(VEG_, dphi)
    dsl4jsb_Get_var2D_onDomain(VEG_, npp_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, n_fixation_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, unit_npp_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, unit_transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, unit_uptake_n_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, unit_uptake_p_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, dphi_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_cn_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_cp_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_np_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, labile_carbon_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, labile_nitrogen_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, labile_phosphorus_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, reserve_carbon_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, reserve_nitrogen_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, reserve_phosphorus_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, press_srf_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, co2_mixing_ratio_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, ga_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, beta_sinklim_ps_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, press_srf_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, co2_mixing_ratio_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, beta_sinklim_ps_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, uptake_n_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_cn_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_np_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, npp_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, beta_sinklim_ps)
    dsl4jsb_Get_var2D_onDomain(VEG_, leaf_cn_direction)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_daytime)
    dsl4jsb_Get_var2D_onDomain(VEG_, press_srf_daytime)
    dsl4jsb_Get_var2D_onDomain(VEG_, co2_mixing_ratio_daytime)
    dsl4jsb_Get_var2D_onDomain(VEG_, ga_daytime)
    dsl4jsb_Get_var2D_onDomain(VEG_, beta_sinklim_ps_daytime)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_jmax_opt_daytime)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_jmax_opt_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(VEG_, press_srf_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(VEG_, co2_mixing_ratio_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(VEG_, ga_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(VEG_, beta_sinklim_ps_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_req_n_tlabile_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_req_p_tlabile_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_tphen_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, npp_tuptake_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, demand_uptake_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, demand_uptake_p_tuptake_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_req_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_req_p_tuptake_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, cost_n_fixation_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, ga_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, leaf2root_troot_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, an_boc_tvegdyn_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, net_growth_tvegdyn_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, lai_tvegdyn_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_tacclim_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_soil_root_tacclim_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_week_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_month_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_jmax_opt_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, w_root_lim_talloc_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, exudation_c_tmyc_mavg)
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onDomain(VEG_, fleaf_sunlit_tfrac_mavg_cl)
    dsl4jsb_Get_var3D_onDomain(VEG_, fleaf_sunlit_tcnl_mavg_cl)
    dsl4jsb_Get_var3D_onDomain(VEG_, fn_oth_cl)
    dsl4jsb_Get_var3D_onDomain(VEG_, leaf_nitrogen_cl)
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(hlp1(SIZE(lai,1), SIZE(lai,2)))
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_STANDALONE__
#define Get_veg_pool
#include "slm_pools_a_incl.f90"
#undef Get_veg_pool
#endif


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    ! ----------------------------------------------------------------------------------------------------- !
    !>1.0 vegetation init specific to the QUINCY model
    !>
    !> QUINCY models: land, plant,soil (case 1) and canopy (case 2)
    !>
    SELECT CASE(TRIM(model%config%quincy_model))

      ! --------------------------------------------------------------------------------------------------- !
      !>  1.1 land & plant & soil model
      !>
      CASE('land','plant','soil')
        target_cn_leaf(:,:)           = lctlib%cn_leaf
        target_np_leaf(:,:)           = lctlib%np_leaf
        target_cn_fine_root(:,:)      = target_cn_leaf(:,:) / root2leaf_cn
        target_np_fine_root(:,:)      = target_np_leaf(:,:) / root2leaf_np

        ! ------------------------------------------------------------------------------------------------- !
        !>  differ between vegetation dynamics schemes: cohort, population, none
        !> 
        SELECT CASE (TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme))
          ! ----------------------------------------------------------------------------------------------- !
          ! start bareground with vegetation dynamics switched on: cohort scheme
          CASE ("cohort")
            ! modify two parameters
            background_mort_rate_tree  = 0.0_wp    ! default: 0.01_wp
            background_mort_rate_grass = 0.0_wp    ! default: 0.05_wp 
            ! individual density
            dens_ind(:,:) = 1.0_wp
            ! reserve C N P
            veg_pool%reserve%domain_c(:,:) = 1.0_wp
            veg_pool%reserve%domain_n(:,:) = veg_pool%reserve%domain_c(:,:) / (1._wp + fresp_growth) / lctlib%cn_leaf
            veg_pool%reserve%domain_p(:,:) = veg_pool%reserve%domain_n(:,:) / lctlib%np_leaf
            ! reserve C13
            IF(lctlib%ps_pathway == ic3phot) THEN
              hlp1(:,:) = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13', lctlib%ps_pathway)
            ELSE
              hlp1(:,:) = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13', lctlib%ps_pathway)
            ENDIF
            veg_pool%reserve%domain_c13(:,:) = calc_fractionation(def_co2_mixing_ratio, def_co2_mixing_ratio_C13, hlp1(:,:)) * &
              & veg_pool%reserve%domain_c(:,:)
            ! reserve C14
            IF(lctlib%ps_pathway == ic3phot) THEN
              hlp1(:,:) = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14', lctlib%ps_pathway)
            ELSE
              hlp1(:,:) = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14', lctlib%ps_pathway)
            ENDIF
            veg_pool%reserve%domain_c14(:,:) = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14-hlp1(:,:)) * & 
              & veg_pool%reserve%domain_c(:,:)
            ! reserve N15
            veg_pool%reserve%domain_n15(:,:) = veg_pool%reserve%domain_n(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          ! ----------------------------------------------------------------------------------------------- !
          ! start bareground with vegetation dynamics switched on: population scheme
          CASE ("population")
            ! individual density
            dens_ind(:,:) = 0.0_wp
            ! seed_bed C N P
            veg_pool%seed_bed%domain_c(:,:) = lctlib%seed_size * 10.0_wp
            veg_pool%seed_bed%domain_n(:,:) = veg_pool%seed_bed%domain_c(:,:) / (1._wp + fresp_growth) / lctlib%cn_leaf
            veg_pool%seed_bed%domain_p(:,:) = veg_pool%seed_bed%domain_n(:,:) / lctlib%np_leaf
            ! seed_bed C13
            IF(lctlib%ps_pathway == ic3phot) THEN
              hlp1(:,:) = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13', lctlib%ps_pathway)
            ELSE
              hlp1(:,:) = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13', lctlib%ps_pathway)
            ENDIF
            veg_pool%seed_bed%domain_c13(:,:) = veg_pool%seed_bed%domain_c(:,:) * &
              & calc_fractionation(def_co2_mixing_ratio, def_co2_mixing_ratio_C13, hlp1(:,:))
            ! seed_bed C14
            IF(lctlib%ps_pathway == ic3phot) THEN
              hlp1(:,:) = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14', lctlib%ps_pathway)
            ELSE
              hlp1(:,:) = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14', lctlib%ps_pathway)
            ENDIF
            veg_pool%seed_bed%domain_c14(:,:) = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14-hlp1(:,:)) * & 
              & veg_pool%seed_bed%domain_c(:,:) 
            ! seed_bed N15
            veg_pool%seed_bed%domain_n15(:,:) = veg_pool%seed_bed%domain_n(:,:) / &
              & ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          ! ----------------------------------------------------------------------------------------------- !
          ! no vegetation dynamics and not starting with bare ground
          CASE ("none")
            ! individual density
            dens_ind(:,:) = 0.4_wp
            ! reserve C N P
            veg_pool%reserve%domain_c(:,:) = 5.0_wp
            veg_pool%reserve%domain_n(:,:) = veg_pool%reserve%domain_c(:,:) / (1._wp + fresp_growth) / lctlib%cn_leaf
            veg_pool%reserve%domain_p(:,:) = veg_pool%reserve%domain_n(:,:) / lctlib%np_leaf
            ! reserve C13
            IF(lctlib%ps_pathway == ic3phot) THEN
              hlp1(:,:) = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13', lctlib%ps_pathway)
            ELSE
              hlp1(:,:) = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13', lctlib%ps_pathway)
            ENDIF
            veg_pool%reserve%domain_c13(:,:) = veg_pool%reserve%domain_c(:,:) * &
              & calc_fractionation(def_co2_mixing_ratio, def_co2_mixing_ratio_C13, hlp1(:,:))
            ! reserve C14
            IF(lctlib%ps_pathway == ic3phot) THEN
              hlp1(:,:) = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14', lctlib%ps_pathway)
            ELSE
              hlp1(:,:) = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14', lctlib%ps_pathway)
            ENDIF
            veg_pool%reserve%domain_c14(:,:) = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14 - hlp1(:,:)) * & 
              & veg_pool%reserve%domain_c(:,:) 
            ! reserve N15
            veg_pool%reserve%domain_n15(:,:) = veg_pool%reserve%domain_n(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          ! ----------------------------------------------------------------------------------------------- !
          ! error catching
          CASE DEFAULT
            WRITE(message_text,'(2a)') 'not a valid veg_dynamics_scheme: ', TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme)
            CALL finish(routine, message_text) 
        ENDSELECT

      ! --------------------------------------------------------------------------------------------------- !
      !>  1.2 canopy model
      !>
      CASE('canopy')
        ! the ref_lai of 6.0 works great for all sites that are not too dry; but use 0.5_wp for dry sites like ZA-Kru
        ref_lai       = 6.0_wp
        lai(:,:)      = ref_lai
        lai_cl(:,:,:) = 0.0_wp  ! this is intentionally wrong; i.e., it will not work at the 1st timestep (which is probably night anyway)

        ! leaf C
        veg_pool%leaf%domain_c(:,:) = lai(:,:) / lctlib%sla
        ! leaf N
        veg_pool%leaf%domain_n(:,:) = lai(:,:) / lctlib%sla / lctlib%cn_leaf
        veg_pool%leaf%domain_p(:,:) = veg_pool%leaf%domain_n(:,:) / lctlib%np_leaf
        IF(lctlib%growthform == itree) THEN
          height(:,:) = 10.0_wp
        ELSE 
          height(:,:) = 1.0_wp
        ENDIF
        leaf_nitrogen_cl(:,:,:) = 0.0_wp  ! this is intentionally wrong; i.e., it will not work at the 1st timestep (which is probably night anyway)

      ! --------------------------------------------------------------------------------------------------- !
      !>  1.3 default - error catching in model selection
      !>
      CASE DEFAULT
        CALL finish(routine, 'Invalid QUINCY model name. Use either land, plant, soil or canopy.')
    ENDSELECT

    ! ----------------------------------------------------------------------------------------------------- !
    !>2.0 vegetation init for all QUINCY models
    !>
    ! plant variables
    dphi(:,:)                             = -1.0_wp * lctlib%phi_leaf_min
    beta_sinklim_ps(:,:)                  = 1.0_wp
    leaf_cn_direction(:,:)                = 0.0_wp
    fn_oth_cl(:,:,:)                      = 1.0_wp  ! it is necessary to start with fn_oth = 1.0 leaving the rest of the fraction 0.0
                                                    ! fn_chl_cl, fn_rub_cl, fn_et_cl, fn_pepc_cl are set to 0.0_wp by default with the add_var()

    ! various moving averages
    growth_req_n_tlabile_mavg(:,:)        = 1._wp / 25._wp
    growth_req_p_tlabile_mavg(:,:)        = 1._wp / 14._wp
    t_air_tphen_mavg(:,:)                 = 283.15_wp
    npp_tuptake_mavg(:,:)                 = 1.0_wp
    demand_uptake_n_tuptake_mavg(:,:)     = 0.1_wp
    demand_uptake_p_tuptake_mavg(:,:)     = 0.1_wp
    growth_req_n_tuptake_mavg(:,:)        = 1._wp / 25._wp
    growth_req_p_tuptake_mavg(:,:)        = 1._wp / 14._wp
    cost_n_fixation_talloc_mavg(:,:)      = 10.0_wp
    ga_tcnl_mavg(:,:)                     = 0.2_wp
    leaf2root_troot_mavg(:,:)             = 1.0_wp
    an_boc_tvegdyn_mavg(:,:)              = 5.0_wp
    net_growth_tvegdyn_mavg(:,:)          = 5.0_wp
    lai_tvegdyn_mavg(:,:)                 = 0.1_wp
    t_air_tacclim_mavg(:,:)               = 283.15_wp     ! this init value may equal t_acclim_zero (283.15) to make work: mo_veg_respiration:temperature_response_respiration
    t_soil_root_tacclim_mavg(:,:)         = 283.15_wp     ! this init value may equal t_acclim_zero (283.15) to make work: mo_veg_respiration:temperature_response_respiration
    t_air_week_mavg(:,:)                  = 283.15_wp
    t_air_month_mavg(:,:)                 = 283.15_wp
    t_jmax_opt_mavg(:,:)                  = 17.0_wp
    
    ! allocation moving averages
    npp_talloc_mavg(:,:)                  = 0.0_wp
    n_fixation_talloc_mavg(:,:)           = 0.0_wp
    transpiration_talloc_mavg(:,:)        = eps8
    unit_transpiration_talloc_mavg(:,:)   = 1.e-6_wp
    unit_uptake_n_talloc_mavg(:,:)        = eps8
    unit_uptake_p_talloc_mavg(:,:)        = 4._wp * 1.e-10_wp
    dphi_talloc_mavg(:,:)                 = -1.0_wp * lctlib%phi_leaf_min
    growth_cn_talloc_mavg(:,:)            = 60.0_wp
    unit_npp_talloc_mavg(:,:)             = unit_uptake_n_talloc_mavg(:,:) * growth_cn_talloc_mavg(:,:)
    growth_cp_talloc_mavg(:,:)            = 2000.0_wp
    growth_np_talloc_mavg(:,:)            = growth_cp_talloc_mavg(:,:) / growth_cn_talloc_mavg(:,:)
    labile_carbon_talloc_mavg(:,:)        = 0.0_wp
    labile_nitrogen_talloc_mavg(:,:)      = 0.0_wp
    labile_phosphorus_talloc_mavg(:,:)    = 0.0_wp
    reserve_carbon_talloc_mavg(:,:)       = 0.0_wp
    reserve_nitrogen_talloc_mavg(:,:)     = 0.0_wp
    reserve_phosphorus_talloc_mavg(:,:)   = 0.0_wp
    w_root_lim_talloc_mavg(:,:)           = w_root_lim_max
    exudation_c_tmyc_mavg(:,:)            = 0.0_wp

    ! N fraction optimisation moving averages
    fleaf_sunlit_tfrac_mavg_cl(:,:,:)     = 0.5_wp
    t_air_tfrac_mavg(:,:)                 = 283.15_wp
    press_srf_tfrac_mavg(:,:)             = standard_press_srf
    co2_mixing_ratio_tfrac_mavg(:,:)      = def_co2_mixing_ratio
    ga_tfrac_mavg(:,:)                    = 0.2_wp
    beta_sinklim_ps_tfrac_mavg(:,:)       = 1.0_wp

    ! leaf N optimisation moving averages
    fleaf_sunlit_tcnl_mavg_cl(:,:,:)      = 0.5_wp
    t_air_tcnl_mavg(:,:)                  = 283.15_wp
    press_srf_tcnl_mavg(:,:)              = standard_press_srf
    co2_mixing_ratio_tcnl_mavg(:,:)       = def_co2_mixing_ratio
    beta_sinklim_ps_tcnl_mavg(:,:)        = 1.0_wp
    growth_cn_tcnl_mavg(:,:)              = 60.0_wp
    growth_np_tcnl_mavg(:,:)              = 30.0_wp
    uptake_n_tcnl_mavg(:,:)               = eps8
    npp_tcnl_mavg(:,:)                    = uptake_n_tcnl_mavg(:,:) * growth_cn_tcnl_mavg(:,:)

    ! daytime averages and their accumulation var
    t_air_daytime(:,:)                    = 283.15_wp
    press_srf_daytime(:,:)                = standard_press_srf
    co2_mixing_ratio_daytime(:,:)         = def_co2_mixing_ratio
    ga_daytime(:,:)                       = 0.2_wp
    beta_sinklim_ps_daytime(:,:)          = 1.0_wp
    t_jmax_opt_daytime(:,:)               = 17.0_wp
    t_jmax_opt_daytime_dacc(:,:)          = 17.0_wp
    t_air_daytime_dacc(:,:)               = 283.15_wp
    press_srf_daytime_dacc(:,:)           = standard_press_srf
    co2_mixing_ratio_daytime_dacc(:,:)    = def_co2_mixing_ratio
    ga_daytime_dacc(:,:)                  = 0.2_wp
    beta_sinklim_ps_daytime_dacc(:,:)     = 1.0_wp

#ifdef __QUINCY_STANDALONE__
! Pass back values to process memory
#define Passback_veg_pool
#include "slm_pools_a_incl.f90"
#undef Passback_veg_pool
#endif

  END SUBROUTINE veg_init

#endif
END MODULE mo_veg_init
