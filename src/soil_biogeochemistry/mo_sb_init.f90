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
!! @brief Initialization of the the soil_biogeochemistry memory.
!!
!! @author
!!  Soenke Zaehle, MPI-BGC Jena
!!
MODULE mo_sb_init
#ifndef __NO_JSBACH__

#ifdef __QUINCY_STANDALONE__
#include "dsl_jsbach.inc"
#endif

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message, finish, message_text
  USE mo_jsb_control,             ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sb_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_init'

CONTAINS

  ! ======================================================================================================= !
  !>Intialize SB_ process memory for QUINCY
  !>
  SUBROUTINE sb_init(tile)

    USE mo_jsb_class,                   ONLY: Get_model
    USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,             ONLY: t_jsb_model, MODEL_JSBACH
    USE mo_jsb_lctlib_class,            ONLY: t_lctlib_element
    USE mo_jsb_grid_class,              ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,                    ONLY: Get_grid, Get_vgrid
    USE mo_jsb_process_class,           ONLY: SB_, VEG_, SPQ_
    ! USE mo_jsb_process_class,           ONLY: SB_, SSE_, VEG_, HYDRO_, SPQ_
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_lnd_bgcm_store,        ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_store_class,  ONLY: SB_BGCM_POOL_ID
    USE mo_lnd_bgcm_class,        ONLY: ELEM_C_ID, ELEM_N_ID, ELEM_P_ID, ELEM_C13_ID, ELEM_C14_ID, ELEM_N15_ID, &
      &                                 LAST_ELEM_ID
    USE mo_sb_config_class,       ONLY: get_number_of_sb_compartments, &
      &                                 SB_PART_DOM_IDX, SB_PART_DOM_ASSOC_IDX, SB_PART_SOLUBLE_LITTER_IDX, &
      &                                 SB_PART_POLYMERIC_LITTER_IDX, SB_PART_WOODY_LITTER_IDX, SB_PART_FUNGI_IDX, &
      &                                 SB_PART_MYCORRHIZA_IDX, SB_PART_MICROBIAL_IDX, SB_PART_RESIDUE_IDX, &
      &                                 SB_PART_RESIDUE_ASSOC_IDX

    USE mo_isotope_util,                ONLY: calc_mixing_ratio_C13C12, calc_mixing_ratio_C14C, &
      & calc_mixing_ratio_N15N14, calc_fractionation
    USE mo_assimi_update_canopy_fluxes, ONLY: discrimination_ps
    USE mo_sb_jsm_processes,            ONLY: calc_bulk_soil_carbon, calc_Psorption_parameter, &
      & calc_fast_po4, calc_qmax_bulk_density_correction, calc_bulk_density_correction
!    USE mo_sb_jsm_transport,            ONLY: calc_bioturbation_rate, ...
    ! USE mo_sse_util,                    ONLY: calc_qmax_texture
    USE mo_spq_util,                    ONLY: calc_qmax_texture

    USE mo_jsb_math_constants,          ONLY: eps4, eps8
    USE mo_jsb_physical_constants,      ONLY: molar_mass_C, molar_mass_N, molar_mass_P
    USE mo_atmland_constants,           ONLY: def_co2_mixing_ratio, def_co2_mixing_ratio_C13, &
      & def_co2_mixing_ratio_C14, def_co2_deltaC13, def_co2_deltaC14
    USE mo_assimi_parameters,           ONLY: CiCa_default_C3, CiCa_default_C4
    USE mo_assimi_constants,            ONLY: ic3phot, ic4phot
    USE mo_veg_constants,               ONLY: itree
    USE mo_sb_constants                 ! e.g., ONLY: nh4_solute_prescribe, no3_solute_prescribe
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(VEG_)
    ! dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_STANDALONE__
    dsl4slm_Use_memory_a_sb
    USE mo_sb_memory_class,       ONLY: sb_pool
#else
    dsl4jsb_Use_pools(SB_) t_sb_components_pool
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile          !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER               :: model                !< model
    TYPE(t_lctlib_element), POINTER               :: lctlib               !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_grid),       POINTER               :: hgrid                !< Horizontal grid
    TYPE(t_jsb_vgrid),      POINTER               :: vgrid_soil_sb        !< Vertical grid
    INTEGER                                       :: nsoil_sb             !< number of soil layers as used/defined by the SB_ process
    REAL(wp)                                      :: frac_c13                         
    REAL(wp)                                      :: frac_c14                       
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)       :: p_labile                       
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)       :: p_labile_slow_site                           
    REAL(wp)                                      :: hlp1                           
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)       :: hlp2                           
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)       :: hlp3                           
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)       :: hlp4                           
    INTEGER                                       :: isoil, iblk, nblks
    INTEGER                                       :: nproma               ! number of grid points per iblk
    TYPE(t_lnd_bgcm_store), POINTER               :: bgcm_store                     !< the bgcm store of this tile
    INTEGER                                       :: elem_idx_map(LAST_ELEM_ID - 1) !< elements ID -> IDX
    INTEGER                                       :: nr_elems                       !< number of elements 
    INTEGER                                       :: nr_sb_parts                    !< number of bgcm compartments
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init'
#ifdef __QUINCY_STANDALONE__
    REAL(wp)                                     :: site_longitude
    LOGICAL                                      :: flag_spp1685
    CHARACTER(len=3)                             :: site_ID_spp1685
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    ! dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(VEG_)
    ! dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_STANDALONE__
#else
    dsl4jsb_Def_pool(t_sb_components_pool)  :: sb_pool
#endif
    REAL(wp), ALLOCATABLE                 :: sb_bgcm_pool_mt(:,:,:,:,:)   !< bgcm sb_pool matrix onDomain
    ! ----------------------------------------------------------------------------------------------------- !
    ! VEG_
    dsl4jsb_Real3D_onDomain :: root_fraction_sl
    ! ! HYDRO_ 2D
    ! dsl4jsb_Real2D_onDomain :: soil_depth
    ! dsl4jsb_Real2D_onDomain :: root_depth
    ! SPQ_ 2D
    dsl4jsb_Real2D_onDomain :: soil_depth
    dsl4jsb_Real2D_onDomain :: root_depth
    ! ! HYDRO_ 3D
    ! dsl4jsb_Real3D_onDomain :: w_soil_sl
    ! ! SSE_ 3D
    ! dsl4jsb_Real3D_onDomain :: bulk_dens_sl
    ! dsl4jsb_Real3D_onDomain :: soil_lay_width_sl
    ! dsl4jsb_Real3D_onDomain :: soil_lay_depth_center_sl
    ! dsl4jsb_Real3D_onDomain :: soil_lay_depth_ubound_sl
    ! dsl4jsb_Real3D_onDomain :: soil_lay_depth_lbound_sl
    ! dsl4jsb_Real3D_onDomain :: qmax_org_min_sl
    ! dsl4jsb_Real3D_onDomain :: qmax_po4_min_sl
    ! dsl4jsb_Real3D_onDomain :: qmax_nh4_min_sl
    ! dsl4jsb_Real3D_onDomain :: qmax_po4_om_sl
    ! dsl4jsb_Real3D_onDomain :: sand_sl
    ! dsl4jsb_Real3D_onDomain :: silt_sl
    ! dsl4jsb_Real3D_onDomain :: clay_sl
    ! dsl4jsb_Real3D_onDomain :: volume_min_sl
    ! SPQ_ 3D
    dsl4jsb_Real3D_onDomain :: w_soil_sl
    dsl4jsb_Real3D_onDomain :: bulk_dens_sl
    dsl4jsb_Real3D_onDomain :: soil_lay_width_sl
    dsl4jsb_Real3D_onDomain :: soil_lay_depth_center_sl
    dsl4jsb_Real3D_onDomain :: soil_lay_depth_ubound_sl
    dsl4jsb_Real3D_onDomain :: soil_lay_depth_lbound_sl
    dsl4jsb_Real3D_onDomain :: qmax_org_min_sl
    dsl4jsb_Real3D_onDomain :: qmax_po4_min_sl
    dsl4jsb_Real3D_onDomain :: qmax_nh4_min_sl
    dsl4jsb_Real3D_onDomain :: qmax_po4_om_sl
    dsl4jsb_Real3D_onDomain :: sand_sl
    dsl4jsb_Real3D_onDomain :: silt_sl
    dsl4jsb_Real3D_onDomain :: clay_sl
    dsl4jsb_Real3D_onDomain :: volume_min_sl
    ! SB_ 3D
    dsl4jsb_Real3D_onDomain :: nh4_assoc
    dsl4jsb_Real3D_onDomain :: nh4_solute
    dsl4jsb_Real3D_onDomain :: no3_solute
    dsl4jsb_Real3D_onDomain :: po4_solute
    dsl4jsb_Real3D_onDomain :: noy
    dsl4jsb_Real3D_onDomain :: n2o
    dsl4jsb_Real3D_onDomain :: n2
    dsl4jsb_Real3D_onDomain :: nh4_n15_assoc
    dsl4jsb_Real3D_onDomain :: nh4_n15_solute
    dsl4jsb_Real3D_onDomain :: no3_n15_solute
    dsl4jsb_Real3D_onDomain :: bulk_soil_carbon_sl
    dsl4jsb_Real3D_onDomain :: soil_litter_carbon_sl
    dsl4jsb_Real3D_onDomain :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onDomain :: qmax_org
    dsl4jsb_Real3D_onDomain :: qmax_po4
    dsl4jsb_Real3D_onDomain :: qmax_nh4
    dsl4jsb_Real3D_onDomain :: qmax_fast_po4
    dsl4jsb_Real3D_onDomain :: qmax_slow_po4
    dsl4jsb_Real3D_onDomain :: km_fast_po4
    dsl4jsb_Real3D_onDomain :: km_slow_po4
    dsl4jsb_Real3D_onDomain :: km_adsorpt_po4_sl
    dsl4jsb_Real3D_onDomain :: km_adsorpt_nh4_sl
    dsl4jsb_Real3D_onDomain :: k_bioturb
    dsl4jsb_Real3D_onDomain :: po4_assoc_fast
    dsl4jsb_Real3D_onDomain :: po4_assoc_slow
    dsl4jsb_Real3D_onDomain :: po4_occluded
    dsl4jsb_Real3D_onDomain :: ph_sl
    dsl4jsb_Real3D_onDomain :: Qmax_AlFe_cor
    dsl4jsb_Real3D_onDomain :: lithology_class_sl
    dsl4jsb_Real3D_onDomain :: po4_primary
    dsl4jsb_Real3D_onDomain :: microbial_cue_eff_tmic_mavg
    dsl4jsb_Real3D_onDomain :: microbial_nue_eff_tmic_mavg
    dsl4jsb_Real3D_onDomain :: microbial_pue_eff_tmic_mavg
    dsl4jsb_Real3D_onDomain :: vmax_weath_mineral_sl
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly
    dsl4jsb_Real3D_onDomain :: enzyme_frac_residue
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly_c
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly_n
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly_p
    dsl4jsb_Real3D_onDomain :: dom_cn
    dsl4jsb_Real3D_onDomain :: dom_cp
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (.NOT. tile%contains_soil) THEN
      IF (debug_on()) CALL message(TRIM(routine), 'Not initializing QUINCY SB_ memory for tile w/o soil: '//TRIM(tile%name))
      RETURN
    END IF
    IF (debug_on()) CALL message(TRIM(routine), 'Initializing QUINCY SB_ memory for tile '//TRIM(tile%name))
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    IF (model%config%model_scheme == MODEL_JSBACH) RETURN     ! only with QUINCY model
    IF (tile%lcts(1)%lib_id == 0) RETURN                      ! run this init only if the present tile is a pft
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    hgrid         => Get_grid(model%grid_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    nblks         =  hgrid%nblks
    nproma        =  hgrid%nproma
    bgcm_store    => tile%bgcm_store
#ifdef __QUINCY_STANDALONE__
    site_longitude = model%config%lon
    flag_spp1685  = model%config%flag_spp1685
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)
    ! dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(VEG_)
    ! dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_STANDALONE__
#define Get_sb_pool
#include "slm_pools_a_incl.f90"
#undef Get_sb_pool
#else
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_pool, sb_pool)
    CALL sb_pool%Associate_domain_var_ptr()
#endif
    elem_idx_map(:) = model%config%elements_index_map(:)
    nr_elems        = model%config%nr_of_used_elements
    nr_sb_parts     = get_number_of_sb_compartments()
    ALLOCATE(sb_bgcm_pool_mt(nr_sb_parts, nr_elems, nproma, nsoil_sb, nblks))
    ! ----------------------------------------------------------------------------------------------------- !
    ! VEG_ 2D
    dsl4jsb_Get_var3D_onDomain(VEG_,   root_fraction_sl)           ! out
    ! ! HYDRO_ 2D
    ! dsl4jsb_Get_var2D_onDomain(HYDRO_, soil_depth)                 ! in
    ! dsl4jsb_Get_var2D_onDomain(HYDRO_, root_depth)                 ! in
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onDomain(SPQ_, soil_depth)                 ! in
    dsl4jsb_Get_var2D_onDomain(SPQ_, root_depth)                 ! in
    ! ! HYDRO_ 3D
    ! dsl4jsb_Get_var3D_onDomain(HYDRO_, w_soil_sl)
    ! ! SSE_ 3D
    ! dsl4jsb_Get_var3D_onDomain(SSE_,   bulk_dens_sl)
    ! dsl4jsb_Get_var3D_onDomain(SSE_,   soil_lay_width_sl)          ! in
    ! dsl4jsb_Get_var3D_onDomain(SSE_,   soil_lay_depth_center_sl)   ! in
    ! dsl4jsb_Get_var3D_onDomain(SSE_,   soil_lay_depth_ubound_sl)   ! in
    ! dsl4jsb_Get_var3D_onDomain(SSE_,   soil_lay_depth_lbound_sl)   ! in
    ! dsl4jsb_Get_var3D_onDomain(SSE_,   qmax_org_min_sl)
    ! dsl4jsb_Get_var3D_onDomain(SSE_,   qmax_po4_min_sl)
    ! dsl4jsb_Get_var3D_onDomain(SSE_,   qmax_nh4_min_sl)
    ! dsl4jsb_Get_var3D_onDomain(SSE_,   qmax_po4_om_sl)
    ! dsl4jsb_Get_var3D_onDomain(SSE_,   sand_sl)
    ! dsl4jsb_Get_var3D_onDomain(SSE_,   silt_sl)
    ! dsl4jsb_Get_var3D_onDomain(SSE_,   clay_sl)
    ! dsl4jsb_Get_var3D_onDomain(SSE_,   volume_min_sl)
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onDomain(SPQ_,   w_soil_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   bulk_dens_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   soil_lay_width_sl)          ! in
    dsl4jsb_Get_var3D_onDomain(SPQ_,   soil_lay_depth_center_sl)   ! in
    dsl4jsb_Get_var3D_onDomain(SPQ_,   soil_lay_depth_ubound_sl)   ! in
    dsl4jsb_Get_var3D_onDomain(SPQ_,   soil_lay_depth_lbound_sl)   ! in
    dsl4jsb_Get_var3D_onDomain(SPQ_,   qmax_org_min_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   qmax_po4_min_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   qmax_nh4_min_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   qmax_po4_om_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   sand_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   silt_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   clay_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   volume_min_sl)
    ! SB_ 3D
    dsl4jsb_Get_var3D_onDomain(SB_,    nh4_assoc)
    dsl4jsb_Get_var3D_onDomain(SB_,    nh4_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    no3_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    noy)
    dsl4jsb_Get_var3D_onDomain(SB_,    n2o)
    dsl4jsb_Get_var3D_onDomain(SB_,    n2)
    dsl4jsb_Get_var3D_onDomain(SB_,    nh4_n15_assoc)
    dsl4jsb_Get_var3D_onDomain(SB_,    nh4_n15_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    no3_n15_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    bulk_soil_carbon_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    soil_litter_carbon_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    bulk_dens_corr_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_org)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_nh4)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_fast_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_slow_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    km_fast_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    km_slow_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    km_adsorpt_po4_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    km_adsorpt_nh4_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    k_bioturb)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_assoc_fast)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_assoc_slow)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_occluded) 
    dsl4jsb_Get_var3D_onDomain(SB_,    ph_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    Qmax_AlFe_cor)
    dsl4jsb_Get_var3D_onDomain(SB_,    lithology_class_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_primary)
    dsl4jsb_Get_var3D_onDomain(SB_,    microbial_cue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onDomain(SB_,    microbial_nue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onDomain(SB_,    microbial_pue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onDomain(SB_,    vmax_weath_mineral_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_poly)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_residue)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_poly_c)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_poly_n)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_poly_p)
    dsl4jsb_Get_var3D_onDomain(SB_,    dom_cn)
    dsl4jsb_Get_var3D_onDomain(SB_,    dom_cp)
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(p_labile(nproma, nsoil_sb, nblks))
    ALLOCATE(p_labile_slow_site(nproma, nsoil_sb, nblks))
    ALLOCATE(hlp2(nproma, nsoil_sb, nblks))
    ALLOCATE(hlp3(nproma, nsoil_sb, nblks))
    ALLOCATE(hlp4(nproma, nsoil_sb, nblks))

    ! -------------------------------------------------------------------------------------------------------
    ! Go science
    ! -------------------------------------------------------------------------------------------------------

    !> 0.9 init bgcm with zero
    !>   note: reading the restart files is done after init
    !>   note: not using a "Get" function in process init for bgcm
    !> 
    sb_bgcm_pool_mt = 0.0_wp

    !>1.0 all soil models (simple_1d & jsm)
    !> 

    !>  1.1 calc frac_c13 & frac_c14 - PFT specific
    !> 
    IF (lctlib%ps_pathway == ic3phot) THEN
      hlp1 = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
    ELSE
      hlp1 = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
    ENDIF
    frac_c13 = calc_fractionation(def_co2_mixing_ratio, def_co2_mixing_ratio_C13, hlp1)

    IF (lctlib%ps_pathway == ic3phot) THEN
      hlp1 = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
    ELSE
      hlp1 = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
    ENDIF
    frac_c14 = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14 - hlp1)

    !>  1.2 site-specific soil properties
    !> 
    ph_sl(:,:,:)              = dsl4jsb_Config(SB_)%soil_ph
    lithology_class_sl(:,:,:) = dsl4jsb_Config(SB_)%lithology_class
    Qmax_AlFe_cor(:,:,:)      = 1._wp

    !>  1.3 init SB_ variables from constants
    !> 
    microbial_cue_eff_tmic_mavg(:,:,:)  = microbial_cue_max
    microbial_nue_eff_tmic_mavg(:,:,:)  = microbial_nue
    microbial_pue_eff_tmic_mavg(:,:,:)  = microbial_pue
    vmax_weath_mineral_sl(:,:,:)        = k_weath_mineral 
    dom_cn(:,:,:)                       = microbial_cn
    dom_cp(:,:,:)                       = microbial_cn * microbial_np

    !>  1.4 root fraction
    !> 
    ! root fractions assuming an exponential distribution of roots:
    !
    !   `Cr_l = Cr_0 * EXP(-lctlib%k_root_dist*depth)`
    !
    ! integral to bottom of root zone: `1 / lctlib%k_root_dist * (1 - EXP(-lctlib%k_root_dist * depth)`
    !
    ! thus for each layer:
    ! root_fraction_sl(:,:) is a variable of veg_mem and depends on the init of root_depth(:), see section 6.0 and ticket #358
    ! 
    ! NOTE: upper bound is the larger value compared to lower bound !
    ! 
    ! layer 1 (vgrid_soil_sb)
    root_fraction_sl(:,1,:) = (1._wp - EXP(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:,1,:))) / &
      &                       (1._wp - EXP(-lctlib%k_root_dist * root_depth(:,:)))
    ! layers 2:x (vgrid_soil_sb)
    DO isoil = 2,nsoil_sb
      WHERE(root_depth(:,:) > soil_lay_depth_ubound_sl(:,isoil-1,:))
        WHERE(root_depth(:,:) > soil_lay_depth_ubound_sl(:,isoil,:))
          root_fraction_sl(:,isoil,:) = (EXP(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:,isoil-1,:)) - &
            &                            EXP(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:,isoil,:))) / & 
            &                           (1._wp - EXP(-lctlib%k_root_dist * root_depth(:,:)))
        ELSEWHERE
          root_fraction_sl(:,isoil,:) = (EXP(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:,isoil-1,:)) - &
            &                            EXP(-lctlib%k_root_dist * root_depth(:,:))) / & 
            &                           (1._wp - EXP(-lctlib%k_root_dist * root_depth(:,:)))
        ENDWHERE
      ELSEWHERE
        root_fraction_sl(:,isoil,:)   = 0.0_wp 
      ENDWHERE
    ENDDO

    ! !------------------------------------------------------------------------------------------------------------------------------
    ! !>  1.x SPP site-set: identify site based on longitude
    ! !>
    ! ! five sites from the "SPP 1685 Project: Forest Strategies for limited Phosphorus Resources"
    ! site_ID_spp1685 = "OTH"   ! default: other
    ! IF (flag_spp1685) THEN
    !   IF (ABS(site_longitude - 9.75_wp) < eps8) THEN
    !     site_ID_spp1685 = "BBR"
    !   ELSEIF (ABS(site_longitude - 10.25_wp) < eps8) THEN
    !     site_ID_spp1685 = "LUE"
    !   ELSEIF (ABS(site_longitude - 10.75_wp) < eps8) THEN
    !     site_ID_spp1685 = "VES"
    !   ELSEIF (ABS(site_longitude - 12.75_wp) < eps8) THEN
    !     site_ID_spp1685 = "MIT"
    !   ELSEIF (ABS(site_longitude - 7.75_wp)  < eps8) THEN
    !     site_ID_spp1685 = "CON"
    !   ENDIF
    ! ENDIF
    ! !------------------------------------------------------------------------------------------------------------------------------

    !>2.0 soil model init
    !>
    ! i)  simple_1d
    ! ii) jsm
    SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
    !!==========================================================================================================================
    ! simple_1d soil model
    ! a four pool soil biogeochemical model
    CASE("simple_1d")
      ! sb_pool carbon
      ! soil layer 1
      ! catch potential division by zero
      ! WHERE ((soil_lay_width_sl(:,1,:) > eps8) .AND. (root_fraction_sl(:,1,:) > eps8))
      WHERE (soil_lay_width_sl(:,1,:) > eps8)
        sb_bgcm_pool_mt(SB_PART_SOLUBLE_LITTER_IDX, elem_idx_map(ELEM_C_ID), :, 1, :)           &
          & = def_sb_pool_metabol_litter /  molar_mass_C / soil_lay_width_sl(:,1,:)            &
          &   * (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1,:))

        sb_pool%soluable_litter%domain_c(:,1,:) = def_sb_pool_metabol_litter /  molar_mass_C / soil_lay_width_sl(:,1,:) * &
          & (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1,:)) 
        sb_pool%polymeric_litter%domain_c(:,1,:) = def_sb_pool_simple_struct_litter / molar_mass_C / soil_lay_width_sl(:,1,:) * &
          & (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1,:))
      END WHERE
      IF (lctlib%growthform == itree) THEN
        ! catch potential division by zero
        WHERE ((soil_lay_width_sl(:,1,:) > eps8))
          sb_pool%woody_litter%domain_c(:,1,:) = def_sb_pool_woody_litter / molar_mass_C / soil_lay_width_sl(:,1,:)
        END WHERE          
      END IF
      ! catch potential division by zero
      WHERE (soil_lay_width_sl(:,1,:) > eps8)
        sb_pool%microbial%domain_c(:,1,:) = def_sb_pool_fast_som / molar_mass_C / soil_lay_width_sl(:,1,:) * &
          & (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1,:)) 
        sb_pool%residue%domain_c(:,1,:)   = def_sb_pool_slow_som / molar_mass_C / soil_lay_width_sl(:,1,:) * & 
          & (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1,:)) 
      END WHERE
      ! soil layers 2:x
      DO isoil = 2,nsoil_sb
        ! catch potential division by zero
        WHERE (soil_lay_width_sl(:,isoil,:) > eps8)
          ! example on how to initialise the bgcm via the matrix written back at the end of this init routine
          sb_bgcm_pool_mt(SB_PART_SOLUBLE_LITTER_IDX, elem_idx_map(ELEM_C_ID), :, isoil, :)           &
            & = def_sb_pool_metabol_litter /  molar_mass_C &
            &   / soil_lay_width_sl(:,isoil,:) * (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil,:)

          sb_pool%soluable_litter%domain_c(:,isoil,:) = def_sb_pool_metabol_litter /  molar_mass_C / &
            & soil_lay_width_sl(:,isoil,:) * (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil,:) 
          sb_pool%polymeric_litter%domain_c(:,isoil,:) = def_sb_pool_simple_struct_litter / molar_mass_C / &
            & soil_lay_width_sl(:,isoil,:) * (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil,:) 
          sb_pool%microbial%domain_c(:,isoil,:) = def_sb_pool_fast_som / molar_mass_C / soil_lay_width_sl(:,isoil,:) * &
            & (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil,:) 
          sb_pool%residue%domain_c(:,isoil,:) = def_sb_pool_slow_som / molar_mass_C / soil_lay_width_sl(:,isoil,:) * & 
            & (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil,:)
        END WHERE
      END DO

      ! note: this differs on purpose from the calculation in mo_sb_update_pools for the prescribed option
      !       this is because here I'm trying to ensure that the models are initialised with the same amount of nutrient
      !       whereas in the prescribed option the target is unlimiting concentrations across the profile. The latter leads
      !       to more nutrient in the soil column in 1d compared to 0d.
      DO isoil = 1,nsoil_sb
        ! avoid division by zero
        WHERE (soil_lay_depth_ubound_sl(:,nsoil_sb,:) > eps8)
          nh4_solute(:,isoil,:) = nh4_solute_prescribe / soil_lay_depth_ubound_sl(:,nsoil_sb,:)
          no3_solute(:,isoil,:) = no3_solute_prescribe / soil_lay_depth_ubound_sl(:,nsoil_sb,:)
        END WHERE
      END DO
      
      ! inital bulk density and qmax_po4 based on current soil carbon content
      DO iblk = 1,nblks
        CALL calc_bulk_soil_carbon( &
          & nproma, &
          & nsoil_sb, &
          & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
          & sb_pool%microbial%domain_c(:,:,iblk), &
          & sb_pool%residue%domain_c(:,:,iblk), &
          & sb_pool%residue_assoc%domain_c(:,:,iblk), &
          & sb_pool%dom%domain_c(:,:,iblk), &
          & sb_pool%dom_assoc%domain_c(:,:,iblk), &
          & sb_pool%soluable_litter%domain_c(:,:,iblk), &
          & sb_pool%polymeric_litter%domain_c(:,:,iblk), &
          & sb_pool%woody_litter%domain_c(:,:,iblk), &
          & bulk_soil_carbon_sl(:,:,iblk), &                ! inout
          & soil_litter_carbon_sl(:,:,iblk), &              ! inout
          & volume_min_sl(:,:,iblk))                        ! inout
      END DO
      bulk_dens_corr_sl(:,:,:) = calc_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,:), &
        & soil_litter_carbon_sl(:,:,:), &
        & bulk_dens_sl(:,:,:))  
      qmax_org_min_sl(:,:,:) = calc_qmax_texture( &
        & dsl4jsb_Config(SB_)%qmax_org_fine_particle, &
        & silt_sl(:,:,:), &
        & clay_sl(:,:,:)) 
      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,:), &   ! in
        & volume_min_sl(:,:,:), &         ! in
        & bulk_dens_sl(:,:,:), &          ! in
        & 0._wp, &                        ! in
        & qmax_org_min_sl(:,:,:), &       ! in
        & qmax_org(:,:,:))                ! out

      ! po4 adsorption parameter intialisation 
      hlp2(:,:,:) = po4_solute_prescribe
      ! ! specific init values for SPP sites for: ph_sl(:,:) & Qmax_AlFe_cor(:,:)
      ! IF (flag_spp1685) THEN
      !   CALL init_soil_properties_spp_sites(nc, nsoil_sb, site_ID_spp1685, ph_sl, Qmax_AlFe_cor)
      ! ENDIF
      ! Correct the first layer Qmax_AlFe_cor with the woody litter
      ! catch runtime error (potentially a division by zero)
      WHERE (soil_lay_width_sl(:,1,:) > eps8)
        Qmax_AlFe_cor(:,1,:) = Qmax_AlFe_cor(:,1,:) * 1._wp / ((1._wp - volume_min_sl(:,1,:)) * &
          & (bulk_soil_carbon_sl(:,1,:) + soil_litter_carbon_sl(:,1,:)  + sb_pool%woody_litter%domain_c(:,1,:)) / &
          & (bulk_soil_carbon_sl(:,1,:) + soil_litter_carbon_sl(:,1,:)) + volume_min_sl(:,1,:))
      END WHERE
      ! catch potential division by zero        
      WHERE (soil_lay_width_sl(:,:,:) > eps8)
        hlp3(:,:,:) = w_soil_sl(:,:,:) / soil_lay_width_sl(:,:,:) * 1000._wp
      END WHERE
      CALL calc_Psorption_parameter( &
        & soil_lay_width_sl(:,:,:), &
        & dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
        & clay_sl(:,:,:), &
        & silt_sl(:,:,:), &
        & hlp3(:,:,:), &
        & ph_sl(:,:,:), &
        & bulk_soil_carbon_sl(:,:,:), &
        & soil_litter_carbon_sl(:,:,:), &
        & sb_pool%woody_litter%domain_c(:,:,:), &
        & bulk_dens_sl(:,:,:), &
        & volume_min_sl(:,:,:), &
        & hlp2(:,:,:), &                            ! in    (po4_solute_prescribe)
        & qmax_po4_min_sl(:,:,:), &                 ! inout            
        & qmax_po4_om_sl(:,:,:), &                  ! inout            
        & qmax_po4(:,:,:), &                        ! inout      
        & km_adsorpt_po4_sl(:,:,:), &               ! inout              
        & qmax_fast_po4(:,:,:), &                   ! inout          
        & qmax_slow_po4(:,:,:), &                   ! inout          
        & km_fast_po4(:,:,:), &                     ! inout        
        & km_slow_po4(:,:,:), &                     ! inout        
        & Qmax_AlFe_cor(:,:,:))                     ! inout        
      ! calculate the partition coefficient between labile and slow P pool based on 
      ! their Qmax (double-Langmuir) or global average values (single-Langmuir)
      IF (dsl4jsb_Config(SB_)%flag_sb_double_langmuir) THEN
        WHERE (qmax_po4(:,:,:) > eps8)
          p_labile_slow_site(:,:,:) = MAX(qmax_fast_po4(:,:,:) / qmax_po4(:,:,:), eps4)
        ELSEWHERE
          p_labile_slow_site(:,:,:) = eps4
        ENDWHERE
      ELSE
        p_labile_slow_site(:,:,:)   = p_labile_slow_global_avg
      ENDIF
      
      ! nh4 adsorption parameter intialisation 
      qmax_nh4_min_sl(:,:,:) = calc_qmax_texture( &
        & qmax_nh4_clay, &
        & 0.0_wp, &
        & clay_sl(:,:,:))
      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,:), &   ! in
        & volume_min_sl(:,:,:), &         ! in
        & bulk_dens_sl(:,:,:), &          ! in
        & 0.0_wp, &                       ! in
        & qmax_nh4_min_sl(:,:,:), &       ! in
        & qmax_nh4(:,:,:))                ! out
      ! avoid potential division by zero
      WHERE (soil_lay_width_sl(:,:,:) > eps8)
        qmax_nh4(:,:,:) = qmax_nh4(:,:,:) * w_soil_sl(:,:,:) / soil_lay_width_sl(:,:,:)
      END WHERE
      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,:), &   ! in
        & volume_min_sl(:,:,:), &         ! in
        & bulk_dens_sl(:,:,:), &          ! in
        & km_adsorpt_OM_nh4, &            ! in
        & km_adsorpt_mineral_nh4, &       ! in
        & km_adsorpt_nh4_sl(:,:,:))       ! out
      ! avoid potential division by zero
      WHERE ((nh4_solute(:,:,:) + km_adsorpt_nh4_sl(:,:,:)) > eps8)
        nh4_assoc(:,:,:)     = qmax_nh4(:,:,:) * nh4_solute (:,:,:) / (nh4_solute(:,:,:) + km_adsorpt_nh4_sl(:,:,:))
      END WHERE
      nh4_n15_assoc(:,:,:) = nh4_assoc(:,:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 

      ! all other elements N, P, C13, C14, N15
      ! metabolic litter
      sb_pool%soluable_litter%domain_n(:,:,:)     = sb_pool%soluable_litter%domain_c(:,:,:) / k_init_soluable_litter_cn 
      sb_pool%soluable_litter%domain_p(:,:,:)     = sb_pool%soluable_litter%domain_n(:,:,:) / k_init_litter_np 
      sb_pool%soluable_litter%domain_c13(:,:,:)   = sb_pool%soluable_litter%domain_c(:,:,:) * frac_c13
      sb_pool%soluable_litter%domain_c14(:,:,:)   = sb_pool%soluable_litter%domain_c(:,:,:) * frac_c14 
      sb_pool%soluable_litter%domain_n15(:,:,:)   = sb_pool%soluable_litter%domain_n(:,:,:) / & 
        & (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
      ! non-woody litter
      sb_pool%polymeric_litter%domain_n(:,:,:)    = sb_pool%polymeric_litter%domain_c(:,:,:) / k_init_polymeric_litter_cn 
      sb_pool%polymeric_litter%domain_p(:,:,:)    = sb_pool%polymeric_litter%domain_n(:,:,:) / k_init_litter_np 
      sb_pool%polymeric_litter%domain_c13(:,:,:)  = sb_pool%polymeric_litter%domain_c(:,:,:) * frac_c13
      sb_pool%polymeric_litter%domain_c14(:,:,:)  = sb_pool%polymeric_litter%domain_c(:,:,:) * frac_c14
      sb_pool%polymeric_litter%domain_n15(:,:,:)  = sb_pool%polymeric_litter%domain_n(:,:,:) / & 
        & (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
      ! woody litter all
      sb_pool%woody_litter%domain_n(:,:,:)    = sb_pool%woody_litter%domain_c(:,:,:) / k_init_woody_litter_cn
      sb_pool%woody_litter%domain_p(:,:,:)    = sb_pool%woody_litter%domain_n(:,:,:) / k_init_litter_np 
      sb_pool%woody_litter%domain_c13(:,:,:)  = sb_pool%woody_litter%domain_c(:,:,:) * frac_c13
      sb_pool%woody_litter%domain_c14(:,:,:)  = sb_pool%woody_litter%domain_c(:,:,:) * frac_c14
      sb_pool%woody_litter%domain_n15(:,:,:)  = sb_pool%woody_litter%domain_n(:,:,:) / & 
        & (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
      ! fast pool (= microbial)
      sb_pool%microbial%domain_n(:,:,:)   = sb_pool%microbial%domain_c(:,:,:) / (k_fast_som_cn_min + k_fast_som_cn_max) * 2.0_wp 
      sb_pool%microbial%domain_p(:,:,:)   = sb_pool%microbial%domain_n(:,:,:) / k_fast_som_np
      sb_pool%microbial%domain_c13(:,:,:) = sb_pool%microbial%domain_c(:,:,:) * frac_c13
      sb_pool%microbial%domain_c14(:,:,:) = sb_pool%microbial%domain_c(:,:,:) * frac_c14
      sb_pool%microbial%domain_n15(:,:,:) = sb_pool%microbial%domain_n(:,:,:) / & 
        & (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
      ! residue
      sb_pool%residue%domain_n(:,:,:)     = sb_pool%residue%domain_c(:,:,:) / k_slow_som_cn
      sb_pool%residue%domain_p(:,:,:)     = sb_pool%residue%domain_n(:,:,:) / k_slow_som_np 
      sb_pool%residue%domain_c13(:,:,:)   = sb_pool%residue%domain_c(:,:,:) * frac_c13
      sb_pool%residue%domain_c14(:,:,:)   = sb_pool%residue%domain_c(:,:,:) * frac_c14 
      sb_pool%residue%domain_n15(:,:,:)   = sb_pool%residue%domain_n(:,:,:) / & 
        & (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

      ! ..
      nh4_n15_solute(:,:,:) = nh4_solute(:,:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
      no3_n15_solute(:,:,:) = no3_solute(:,:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
      
      ! init Soil P pools
      !
      !! SPP site-set
      ! IF (flag_spp1685) THEN
      !   ! for SPP sites, read in initial state from site-specific input data
      !   CALL init_soil_p_pools_spp_sites(nc, nsoil_sb, soil_lay_depth_center_sl(:,:), & 
      !                 dsl4jsb_Config(SB_)%soil_p_labile,   dsl4jsb_Config(SB_)%soil_p_slow, &
      !                 dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
      !                 site_ID_spp1685, p_labile_slow_site(:,:), &
      !                 p_labile(:,:), po4_assoc_slow(:,:), po4_occluded(:,:), po4_primary(:,:))
      !   IF (dsl4jsb_Config(SB_)%flag_sb_double_langmuir) THEN
      !     po4_assoc_fast(:,:) = calc_fast_po4(p_labile(:,:),qmax_fast_po4,km_fast_po4)
      !   ELSE
      !     po4_assoc_fast(:,:) = calc_fast_po4(p_labile(:,:),qmax_po4(:,:),km_adsorpt_po4_sl(:,:))
      !   ENDIF
      !   po4_solute(:,:)     = p_labile(:,:) - po4_assoc_fast(:,:)
      ! 
      !! all sites but the SPP site-set
      ! ELSE
        ! initialise mineral soil P pools
        DO iblk = 1,nblks
          CALL init_soil_p_pools( &
            & nproma, &
            & nsoil_sb, &
            & soil_depth(:,iblk), &
            & soil_lay_width_sl(:,:,iblk), &
            & soil_lay_depth_center_sl(:,:,iblk), &
            & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), & 
            & TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &
            & dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
            & bulk_dens_corr_sl(:,:,iblk), &
            & dsl4jsb_Config(SB_)%soil_p_labile, &
            & dsl4jsb_Config(SB_)%soil_p_slow, &
            & dsl4jsb_Config(SB_)%soil_p_occluded, &
            & dsl4jsb_Config(SB_)%soil_p_primary, &
            & dsl4jsb_Config(SB_)%soil_p_depth, &
            & p_labile_slow_site(:,:,iblk), & 
            & po4_solute(:,:,iblk), &                     ! inout
            & po4_assoc_fast(:,:,iblk), &                 ! inout
            & po4_assoc_slow(:,:,iblk), &                 ! inout
            & po4_occluded(:,:,iblk), &                   ! inout
            & po4_primary(:,:,iblk), &                    ! inout
            & qmax_po4(:,:,iblk), &                       ! inout
            & km_adsorpt_po4_sl(:,:,iblk), &              ! inout  argument name within routin: km_po4
            & qmax_fast_po4(:,:,iblk), &                  ! inout
            & qmax_slow_po4(:,:,iblk), &                  ! inout
            & km_fast_po4(:,:,iblk), &                    ! inout
            & km_slow_po4(:,:,iblk))                      ! inout
          END DO
      ! ENDIF
    !!==========================================================================================================================
    !> JSM (Jena Soil Model)
    !! published with: https://doi.org/10.5194/gmd-13-783-2020
    ! CASE('jsm')
    !   ! initialisation of litter pools
    !   sb_pool%soluable_litter%carbon(:,1)  = def_sb_pool_metabol_litter /  molar_mass_C / soil_lay_width_sl(:,1) * &
    !                                          (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1)) 
    !   sb_pool%polymeric_litter%carbon(:,1) = def_sb_pool_jsm_struct_litter / molar_mass_C / soil_lay_width_sl(:,1) * &
    !                                          (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1)) 
    !   IF(lctlib%growthform==itree) THEN
    !     sb_pool%woody_litter%carbon(:,1)  = def_sb_pool_woody_litter / molar_mass_C / soil_lay_width_sl(:,1)
    !   ENDIF
    !   DO isoil=2,nsoil_sb
    !     sb_pool%soluable_litter%carbon(:,isoil)  = def_sb_pool_metabol_litter /  molar_mass_C / soil_lay_width_sl(:,isoil) * &
    !                                                (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil) 
    !     sb_pool%polymeric_litter%carbon(:,isoil) = def_sb_pool_jsm_struct_litter / molar_mass_C / soil_lay_width_sl(:,isoil) * &
    !                                                (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil) 
    !   ENDDO
    !   ! for each soil layer
    !   DO isoil=1,nsoil_sb
    !     sb_pool%microbial%carbon(:,isoil)  = MAX(def_sb_pool_microbial_biomass / soil_lay_width_sl(:,isoil) * &
    !                                              root_fraction_sl(:,isoil), min_mic_biomass)
    !     sb_pool%dom%carbon(:,isoil)        = def_sb_pool_dom / soil_lay_width_sl(:,isoil) * root_fraction_sl(:,isoil)
    !     sb_pool%dom_assoc%carbon(:,isoil)  = k_sb_pool_dom_assoc_c_corr * sb_pool%dom%carbon(:,isoil)
    !     sb_pool%residue%carbon(:,isoil)    = def_sb_pool_microbial_necromass / soil_lay_width_sl(:,isoil) * &
    !                                          root_fraction_sl(:,isoil)
    !     frac_prod_not_root                 = 1._wp / isoil
    !     ! deep soil is initialised with very low MOC, corrected by k_deep_soil_moc_corr
    !     IF (frac_prod_not_root < k_deep_soil_moc_corr) THEN 
    !       frac_prod_not_root = frac_prod_not_root * k_deep_soil_moc_corr
    !     ENDIF
    !     sb_pool%residue_assoc%carbon(:,isoil) = dsl4jsb_Config(SB_)%qmax_org_fine_particle * &
    !                                             (clay_sl(:,isoil) + silt_sl(:,isoil)) * & 
    !                                             bulk_dens_sl(:,isoil) * frac_prod_not_root
    !     ! To avoid over-saturation of OM in soil layer
    !     ! calculate the fast C and stable C content (saturation at 1._wp) for each layer
    !     IF (isoil==1) THEN
    !       hlp2(:,isoil) = (sb_pool%microbial%carbon(:,isoil) + sb_pool%dom%carbon(:,isoil) + &
    !                           sb_pool%dom_assoc%carbon(:,isoil)) * &
    !                           molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
    !       hlp3(:,isoil) = (sb_pool%residue%carbon(:,isoil) + sb_pool%residue_assoc%carbon(:,isoil)) * &
    !                           molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
    !     ELSE
    !       hlp2(:,isoil) = (sb_pool%soluable_litter%carbon(:,isoil) + sb_pool%polymeric_litter%carbon(:,isoil) + & 
    !                           sb_pool%microbial%carbon(:,isoil) + sb_pool%dom%carbon(:,isoil) + &
    !                           sb_pool%dom_assoc%carbon(:,isoil)) * &
    !                           molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
    !       hlp3(:,isoil) = (sb_pool%residue%carbon(:,isoil) + sb_pool%residue_assoc%carbon(:,isoil)) * &
    !                           molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
    !     ENDIF
    !     ! check if soil is over-saturated with fast C
    !     WHERE (hlp2(:,isoil) > 1._wp )
    !       sb_pool%soluable_litter%carbon(:,isoil)  = sb_pool%soluable_litter%carbon(:,isoil)  * 1._wp / hlp2(:,isoil)
    !       sb_pool%polymeric_litter%carbon(:,isoil) = sb_pool%polymeric_litter%carbon(:,isoil) * 1._wp / hlp2(:,isoil)
    !       sb_pool%microbial%carbon(:,isoil)        = sb_pool%microbial%carbon(:,isoil)        * 1._wp / hlp2(:,isoil)
    !       sb_pool%dom%carbon(:,isoil)              = sb_pool%dom%carbon(:,isoil)              * 1._wp / hlp2(:,isoil)
    !       sb_pool%dom_assoc%carbon(:,isoil)        = sb_pool%dom_assoc%carbon(:,isoil)        * 1._wp / hlp2(:,isoil)
    !       sb_pool%residue%carbon(:,isoil)          = 0._wp
    !       sb_pool%residue_assoc%carbon(:,isoil)    = 0._wp
    !     ! check if soil is over-saturated with fast and stable C
    !     ELSEWHERE ((hlp2(:,isoil) + hlp3(:,isoil)) > 1._wp )
    !       hlp4(:,isoil)                     = (1._wp - hlp2(:,isoil)) * f1_repartition_om_init * &
    !                                               dsl4jsb_Config(SB_)%qmax_org_fine_particle * &
    !                                               (clay_sl(:,isoil) + silt_sl(:,isoil)) * bulk_dens_sl(:,isoil)
    !       hlp3(:,isoil)                     = (1._wp - hlp2(:,isoil)) * (1._wp - f1_repartition_om_init) * &
    !                                               rho_bulk_org * 1000._wp * carbon_per_dryweight_SOM / molar_mass_C
    !       sb_pool%residue_assoc%carbon(:,isoil) = calc_fast_po4(hlp3(:,isoil), hlp4(:,isoil), &
    !                                                             k_desorpt_det / k_adsorpt_det)
    !       sb_pool%residue%carbon(:,isoil)       = hlp3(:,isoil) - sb_pool%residue_assoc%carbon(:,isoil)
    !     ! check if the available mineral surface is enough for sorbed OM
    !     ELSEWHERE
    !       hlp4(:,isoil) = (1._wp - hlp2(:,isoil) - hlp3(:,isoil)) * dsl4jsb_Config(SB_)%qmax_org_fine_particle * &
    !                           (clay_sl(:,isoil) + silt_sl(:,isoil)) * bulk_dens_sl(:,isoil)
    !       hlp3(:,isoil) = (sb_pool%residue%carbon(:,isoil) + sb_pool%residue_assoc%carbon(:,isoil))
    !       WHERE (hlp3(:,isoil) > hlp4(:,isoil))
    !         sb_pool%soluable_litter%carbon(:,isoil)  = sb_pool%soluable_litter%carbon(:,isoil) + &
    !                                                    (hlp3(:,isoil)-hlp4(:,isoil)) * f2_repartition_om_init
    !         sb_pool%polymeric_litter%carbon(:,isoil) = sb_pool%polymeric_litter%carbon(:,isoil) + &
    !                                                    (hlp3(:,isoil)-hlp4(:,isoil)) * (1._wp - f2_repartition_om_init)
    !         hlp3(:,isoil)                        = hlp4(:,isoil)
    !       ENDWHERE
    !       sb_pool%residue_assoc%carbon(:,isoil) = calc_fast_po4(hlp3(:,isoil), hlp4(:,isoil), &
    !                                                             k_desorpt_det / k_adsorpt_det)
    !       sb_pool%residue%carbon(:,isoil)       = hlp3(:,isoil) - sb_pool%residue_assoc%carbon(:,isoil)
    !     ENDWHERE
    !   ENDDO
    ! 
    !   ! initial values of soluble nh4 and no3
    !   DO ichunk = 1,nc
    !     nh4_solute(ichunk,:)      = nh4_solute_prescribe / soil_lay_depth_ubound_sl(ichunk,nsoil_sb)
    !     nh4_n15_solute(ichunk,:)  = nh4_solute(ichunk,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
    !     no3_solute(ichunk,:)      = no3_solute_prescribe / soil_lay_depth_ubound_sl(ichunk,nsoil_sb)
    !     no3_n15_solute(ichunk,:)  = no3_solute(ichunk,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
    !   ENDDO
    !   ! inital bulk density, Qmax_org & Qmax_po4 based on current soil carbon content
    !   CALL calc_bulk_soil_carbon( &
    !     & nc, &
    !     & nsoil_sb, &
    !     & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
    !     & sb_pool%microbial%carbon(:,:), &
    !     & sb_pool%residue%carbon(:,:), &
    !     & sb_pool%residue_assoc%carbon(:,:), &
    !     & sb_pool%dom%carbon(:,:), &
    !     & sb_pool%dom_assoc%carbon(:,:), &
    !     & sb_pool%soluable_litter%carbon(:,:), &
    !     & sb_pool%polymeric_litter%carbon(:,:), &
    !     & sb_pool%woody_litter%carbon(:,:), &
    !     & bulk_soil_carbon_sl(:,:), &                ! inout
    !     & soil_litter_carbon_sl(:,:), &              ! inout
    !     & volume_min_sl(:,:))                        ! inout
    !   bulk_dens_corr_sl(:,:) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:), soil_litter_carbon_sl(:,:), &
    !                                                         bulk_dens_sl(:,:))
    !   !> po4/OC adsorption parameter intialisation 
    !   !! OC adsorption
    !   qmax_org_min_sl(:,:)   = calc_qmax_texture(dsl4jsb_Config(SB_)%qmax_org_fine_particle, silt_sl(:,:), clay_sl(:,:))
    !   CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), &
    !                                          bulk_dens_sl(:,:), 0._wp, qmax_org_min_sl(:,:), qmax_org(:,:))
    !   qmax_org_min_sl(:,:)   = calc_qmax_texture(dsl4jsb_Config(SB_)%qmax_org_fine_particle, silt_sl(:,:), clay_sl(:,:)) 
    !   CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), bulk_dens_sl(:,:), 0._wp, &
    !                                          qmax_org_min_sl(:,:), qmax_org(:,:))   
    !   !! po4 adsorption
    !   hlp2(:,:) = po4_solute_prescribe         
    !   ! specific init values for SPP sites for: ph_sl(:,:) & Qmax_AlFe_cor(:,:)
    !   IF (flag_spp1685) THEN
    !     CALL init_soil_properties_spp_sites(nc, nsoil_sb, site_ID_spp1685, ph_sl, Qmax_AlFe_cor)
    !   ENDIF
    !   ! Correct the first layer Qmax_AlFe_cor with the woody litter
    !   Qmax_AlFe_cor(:,1) = Qmax_AlFe_cor(:,1) * 1._wp / ((1._wp - volume_min_sl(:,1)) * &
    !                        (bulk_soil_carbon_sl(:,1) + soil_litter_carbon_sl(:,1)  + sb_pool%woody_litter%carbon(:,1)) / &
    !                        (bulk_soil_carbon_sl(:,1) + soil_litter_carbon_sl(:,1)) + volume_min_sl(:,1))
    !   CALL calc_Psorption_parameter(dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
    !                                 clay_sl(:,:), silt_sl(:,:), w_soil_sl(:,:) / soil_lay_width_sl(:,:) * 1000._wp, ph_sl(:,:), &
    !                                 bulk_soil_carbon_sl(:,:), soil_litter_carbon_sl(:,:), sb_pool%woody_litter%carbon(:,:), &
    !                                 bulk_dens_sl(:,:), volume_min_sl(:,:), hlp2(:,:), &
    !                                 qmax_po4_min_sl(:,:), qmax_po4_om_sl(:,:), qmax_po4(:,:), km_adsorpt_po4_sl(:,:), &
    !                                 qmax_fast_po4(:,:), qmax_slow_po4(:,:), km_fast_po4(:,:), km_slow_po4(:,:), &
    !                                 Qmax_AlFe_cor(:,:))
    !   ! calculate the partition coefficient between labile and slow P pool
    !   IF (dsl4jsb_Config(SB_)%flag_sb_double_langmuir) THEN
    !     WHERE (qmax_po4(:,:) > eps8)
    !       p_labile_slow_site(:,:)   = MAX(qmax_fast_po4(:,:) / qmax_po4(:,:),eps4)
    !     ELSEWHERE
    !       p_labile_slow_site(:,:)   = eps4
    !     ENDWHERE
    !   ELSE
    !     p_labile_slow_site(:,:)     = p_labile_slow_global_avg
    !   ENDIF
    !   !> nh4 adsorption parameter intialisation
    !   !!
    !   qmax_nh4_min_sl(:,:)     = calc_qmax_texture(qmax_nh4_clay, 0.0_wp, clay_sl(:,:))
    !   CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), bulk_dens_sl(:,:), &
    !                                          0.0_wp,qmax_nh4_min_sl(:,:), qmax_nh4(:,:))
    !   qmax_nh4(:,:)            = qmax_nh4(:,:) * w_soil_sl(:,:) / soil_lay_width_sl(:,:)
    !   CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), bulk_dens_sl(:,:), &
    !                                          km_adsorpt_OM_nh4, km_adsorpt_mineral_nh4, km_adsorpt_nh4_sl(:,:))
    !   nh4_assoc(:,:)           = qmax_nh4(:,:) * nh4_solute (:,:) / (nh4_solute(:,:) + km_adsorpt_nh4_sl(:,:))
    !   nh4_n15_assoc(:,:)       = nh4_assoc(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
    !   !> Soil P pools intialisation
    !   !! 
    !   !! SPP site-set
    !   IF (flag_spp1685) THEN
    !     ! for SPP sites, read in initial state from site-specific input data
    !     CALL init_soil_p_pools_spp_sites(nc, nsoil_sb, soil_lay_depth_center_sl(:,:), & 
    !                      dsl4jsb_Config(SB_)%soil_p_labile,   dsl4jsb_Config(SB_)%soil_p_slow, &
    !                      dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
    !                      site_ID_spp1685, p_labile_slow_site(:,:), &
    !                      p_labile(:,:), po4_assoc_slow(:,:), po4_occluded(:,:), po4_primary(:,:))
    !     ! calculate the soluble and adsorbed P pool based on Langmuir isotherm
    !     IF (dsl4jsb_Config(SB_)%flag_sb_double_langmuir) THEN
    !       po4_assoc_fast(:,:) = calc_fast_po4(p_labile(:,:),qmax_fast_po4,km_fast_po4)
    !     ELSE
    !       po4_assoc_fast(:,:) = calc_fast_po4(p_labile(:,:),qmax_po4(:,:),km_adsorpt_po4_sl(:,:))
    !     ENDIF
    !     po4_solute(:,:)     = p_labile(:,:) - po4_assoc_fast(:,:)
    !   !! all sites but the SPP site-set
    !   ELSE
    !     ! mineral soil P pools intialisation
    !     CALL init_soil_p_pools(nc, nsoil_sb, soil_lay_depth_center_sl(:,:), TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), & 
    !                           TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
    !                           soil_lay_width_sl(:,:), bulk_dens_corr_sl(:,:), root_fraction_sl(:,:), &
    !                           dsl4jsb_Config(SB_)%soil_p_labile, dsl4jsb_Config(SB_)%soil_p_slow, &
    !                           dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
    !                           dsl4jsb_Config(SB_)%soil_p_depth, p_labile_slow_site(:,:), & 
    !                           po4_solute(:,:), po4_assoc_fast(:,:), po4_assoc_slow(:,:), &
    !                           po4_occluded(:,:), po4_primary(:,:), qmax_po4(:,:), km_adsorpt_po4_sl(:,:), &
    !                           qmax_fast_po4(:,:), qmax_slow_po4(:,:), km_fast_po4(:,:), km_slow_po4(:,:))
    !   ENDIF ! (flag_spp1685)
    !   ! initialise diffusion by bioturbation
    !   k_bioturb(:,:) = calc_bioturbation_rate(nc, nsoil_sb,bulk_dens_corr_sl(:,:), root_fraction_sl(:,:), soil_lay_width_sl(:,:)) 
    ! 
    !   ! initialisation of fractional enzyme allocation (range [0,1]) 
    !   enzyme_frac_poly(:,:)     = init_enzyme_fraction
    !   enzyme_frac_residue(:,:)  = init_enzyme_fraction 
    !   enzyme_frac_poly_c(:,:)   = init_enzyme_fraction
    !   enzyme_frac_poly_n(:,:)   = init_enzyme_fraction
    !   enzyme_frac_poly_p(:,:)   = init_enzyme_fraction
    !   ! isotopic signals in all pools
    !   ! microbial, DOM, aDOM, residue, and aRes pools
    !   sb_pool%soluable_litter%nitrogen(:,:)     = sb_pool%soluable_litter%carbon(:,:) / k_init_soluable_litter_cn 
    !   sb_pool%soluable_litter%phosphorus(:,:)   = sb_pool%soluable_litter%nitrogen(:,:) / k_init_litter_np 
    !   sb_pool%soluable_litter%carbon13(:,:)     = sb_pool%soluable_litter%carbon(:,:) * frac_c13 
    !   sb_pool%soluable_litter%carbon14(:,:)     = sb_pool%soluable_litter%carbon(:,:) * frac_c14 
    !   sb_pool%soluable_litter%nitrogen15(:,:)   = sb_pool%soluable_litter%nitrogen(:,:) / &
    !                                               (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    !   sb_pool%polymeric_litter%nitrogen(:,:)    = sb_pool%polymeric_litter%carbon(:,:) / k_init_polymeric_litter_cn
    !   sb_pool%polymeric_litter%phosphorus(:,:)  = sb_pool%polymeric_litter%nitrogen(:,:) / k_init_litter_np 
    !   sb_pool%polymeric_litter%carbon13(:,:)    = sb_pool%polymeric_litter%carbon(:,:) * frac_c13 
    !   sb_pool%polymeric_litter%carbon14(:,:)    = sb_pool%polymeric_litter%carbon(:,:) * frac_c14 
    !   sb_pool%polymeric_litter%nitrogen15(:,:)  = sb_pool%polymeric_litter%nitrogen(:,:) / &
    !                                               (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    !   sb_pool%woody_litter%nitrogen(:,:)        = sb_pool%woody_litter%carbon(:,:) / k_init_woody_litter_cn
    !   sb_pool%woody_litter%phosphorus(:,:)      = sb_pool%woody_litter%nitrogen(:,:) / k_init_litter_np
    !   sb_pool%woody_litter%carbon13(:,:)        = sb_pool%woody_litter%carbon(:,:) * frac_c13 
    !   sb_pool%woody_litter%carbon14(:,:)        = sb_pool%woody_litter%carbon(:,:) * frac_c14 
    !   sb_pool%woody_litter%nitrogen15(:,:)      = sb_pool%woody_litter%nitrogen(:,:) / &
    !                                               (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    !   sb_pool%microbial%nitrogen(:,:)           = sb_pool%microbial%carbon(:,:) / microbial_cn 
    !   sb_pool%microbial%phosphorus(:,:)         = sb_pool%microbial%nitrogen(:,:) / microbial_np
    !   sb_pool%microbial%carbon13(:,:)           = sb_pool%microbial%carbon(:,:) * frac_c13 
    !   sb_pool%microbial%carbon14(:,:)           = sb_pool%microbial%carbon(:,:) * frac_c14 
    !   sb_pool%microbial%nitrogen15(:,:)         = sb_pool%microbial%nitrogen(:,:) / &
    !                                               (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    !   sb_pool%residue%nitrogen(:,:)             = sb_pool%residue%carbon(:,:) / microbial_cn 
    !   sb_pool%residue%phosphorus(:,:)           = sb_pool%residue%nitrogen(:,:) / k_slow_som_np
    !   sb_pool%residue%carbon13(:,:)             = sb_pool%residue%carbon(:,:) * frac_c13 
    !   sb_pool%residue%carbon14(:,:)             = sb_pool%residue%carbon(:,:) * frac_c14 
    !   sb_pool%residue%nitrogen15(:,:)           = sb_pool%residue%nitrogen(:,:) / &
    !                                               (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    !   sb_pool%residue_assoc%nitrogen(:,:)       = sb_pool%residue_assoc%carbon(:,:) / microbial_cn
    !   sb_pool%residue_assoc%phosphorus(:,:)     = sb_pool%residue_assoc%nitrogen(:,:) / k_slow_som_np
    !   sb_pool%residue_assoc%carbon13(:,:)       = sb_pool%residue_assoc%carbon(:,:) * frac_c13 
    !   sb_pool%residue_assoc%carbon14(:,:)       = sb_pool%residue_assoc%carbon(:,:) * frac_c14 
    !   sb_pool%residue_assoc%nitrogen15(:,:)     = sb_pool%residue_assoc%nitrogen(:,:) / &
    !                                               (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    !   sb_pool%dom%nitrogen(:,:)                 = sb_pool%dom%carbon(:,:) / microbial_cn
    !   sb_pool%dom%phosphorus(:,:)               = sb_pool%dom%nitrogen(:,:) / microbial_np
    !   sb_pool%dom%carbon13(:,:)                 = sb_pool%dom%carbon(:,:) * frac_c13
    !   sb_pool%dom%carbon14(:,:)                 = sb_pool%dom%carbon(:,:) * frac_c14
    !   sb_pool%dom%nitrogen15(:,:)               = sb_pool%dom%nitrogen(:,:) / &
    !                                               (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    !   sb_pool%dom_assoc%nitrogen(:,:)           = sb_pool%dom_assoc%carbon(:,:) / microbial_cn
    !   sb_pool%dom_assoc%phosphorus(:,:)         = sb_pool%dom_assoc%nitrogen(:,:) / microbial_np
    !   sb_pool%dom_assoc%carbon13(:,:)           = sb_pool%dom_assoc%carbon(:,:) * frac_c13
    !   sb_pool%dom_assoc%carbon14(:,:)           = sb_pool%dom_assoc%carbon(:,:) * frac_c14
    !   sb_pool%dom_assoc%nitrogen15(:,:)         = sb_pool%dom_assoc%nitrogen(:,:) / & 
    !                                               (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    !   ! fungi and mycorrhiza      
    !   sb_pool%fungi%carbon13(:,:)               = sb_pool%fungi%carbon(:,:) * frac_c13
    !   sb_pool%fungi%carbon14(:,:)               = sb_pool%fungi%carbon(:,:) * frac_c14 
    !   sb_pool%fungi%nitrogen15(:,:)             = sb_pool%fungi%nitrogen(:,:) / & 
    !                                               (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    !   sb_pool%mycorrhiza%carbon13(:,:)          = sb_pool%mycorrhiza%carbon(:,:) * frac_c13
    !   sb_pool%mycorrhiza%carbon14(:,:)          = sb_pool%mycorrhiza%carbon(:,:) * frac_c14 
    !   sb_pool%mycorrhiza%nitrogen15(:,:)        = sb_pool%mycorrhiza%nitrogen(:,:) / & 
    !                                              (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    CASE DEFAULT
      CALL finish(routine,'Invalid setting of sb_model_scheme.')
    END SELECT  ! simple_1d | jsm

  ! Deallocate local allocatable variables
  DEALLOCATE(hlp2)
  DEALLOCATE(hlp3)
  DEALLOCATE(hlp4)
  DEALLOCATE(p_labile)
  DEALLOCATE(p_labile_slow_site)


#ifdef __QUINCY_STANDALONE__
! Pass back values to process memory
#define Passback_sb_pool
#include "slm_pools_a_incl.f90"
#undef Passback_sb_pool
#endif

    ! write back data from bgcm_matrix to bgc_material, and de-allocate local bgcm_matrix
    DO iblk = 1,nblks
      CALL bgcm_store%Write_matrix_to_bgcm_2l( &
        & SB_BGCM_POOL_ID, 1, nproma, iblk, tile%name, sb_bgcm_pool_mt(:,:,:,:,iblk))
    ENDDO
    DEALLOCATE(sb_bgcm_pool_mt)

  END SUBROUTINE sb_init

  ! ======================================================================================================= !
  !>Intialize soil P pools for all layers from the soil P maps (Yang et al. 2013)
  !>assuming that the soil P pool sizes are correlated with soil weight 
  !> 
  SUBROUTINE init_soil_p_pools( &
    & nc, &
    & nsoil_sb, &
    & soil_depth, &
    & soil_lay_width_sl, &
    & soil_lay_depth_center_sl, &
    & sb_model_scheme, &
    & sb_adsorp_scheme, &
    & flag_sb_double_langmuir, & 
    & bulk_dens_sl, &
    & soil_p_labile, &
    & soil_p_slow, &
    & soil_p_occluded, &
    & soil_p_primary, &
    & soil_p_depth, &
    & p_labile_slow_ratio, &
    & po4_solute, &
    & po4_assoc_fast, &
    & po4_assoc_slow, &
    & po4_occluded, &
    & po4_primary, &
    & qmax_po4, &
    & km_po4, &
    & qmax_fast_po4, &
    & qmax_slow_po4, &
    & km_fast_po4, &
    & km_slow_po4)

    USE mo_sb_jsm_processes,          ONLY: calc_fast_po4
    USE mo_jsb_math_constants,        ONLY: eps8
    USE mo_jsb_physical_constants,    ONLY: molar_mass_P
    USE mo_sb_constants               ! e.g.: k1_p_pool_depth_corr, k2_p_pool_depth_corr, k3_p_pool_depth_corr, &
                                      !       k4_p_pool_depth_corr, k_p_pool_depth_corr_min_occlud, p_labile_slow_global_avg
    !------------------------------------------------------------------------------------------------------ !
    INTEGER,                            INTENT(in)     :: nc                        !< number of points at the chunk
    INTEGER,                            INTENT(in)     :: nsoil_sb                  !< number of soil layers as used/defined by the SB_ process
    REAL(wp), DIMENSION(nc),            INTENT(in)     :: soil_depth                !< total soil depth until bedrock (from textures)
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)     :: soil_lay_width_sl         !< soil width per soil layer, [m]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)     :: soil_lay_depth_center_sl  !< depth at the center of each soil layer
    CHARACTER(len=*),                   INTENT(in)     :: sb_model_scheme           !< soil model scheme 
    CHARACTER(len=*),                   INTENT(in)     :: sb_adsorp_scheme          !< soil model po4 adsorption scheme
    LOGICAL,                            INTENT(in)     :: flag_sb_double_langmuir   !< TRUE: double Langmuir; FALSE: traditional (single) Langmuir    
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)     :: bulk_dens_sl              !< corrected soil bulk density per soil layer, [kg m-3]
    REAL(wp),                           INTENT(in)     :: soil_p_labile             !< total labile P pool (solute + assoc_fast) of the soil P maps, [g m-2]
    REAL(wp),                           INTENT(in)     :: soil_p_slow               !< total slow P pool (secondary P) of the soil P maps, [g m-2]
    REAL(wp),                           INTENT(in)     :: soil_p_occluded           !< total occluded P pool of the soil P maps, [g m-2]
    REAL(wp),                           INTENT(in)     :: soil_p_primary            !< total primary P pool of the soil P maps, [g m-2]
    REAL(wp),                           INTENT(in)     :: soil_p_depth              !< soil depth of the soil P maps, i.e., P related init data until this soil depth
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)     :: p_labile_slow_ratio       !< ratio between labile and labile+slow P pools
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout)  :: po4_solute                !< soluble Pi pool per soil layer, [mol m-3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout)  :: po4_assoc_fast            !< fast minerally associated PO4 pool per soil layer, [mol m-3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout)  :: po4_assoc_slow            !< slow minerally associated PO4 pool per soil layer, [mol m-3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout)  :: po4_occluded              !< occluded PO4 pool per soil layer, [mol m-3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout)  :: po4_primary               !< primary PO4 pool per soil layer, [mol m-3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout)  :: qmax_po4                  !< maximum PO4 sorption capacity per soil layer, [mol m-3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout)  :: km_po4                    !< half-saturation PO4 concentration for sorption per soil layer, [mol m-3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout)  :: qmax_fast_po4             !< maximum PO4 sorption capacity for po4_assoc_fast, [mol m-3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout)  :: qmax_slow_po4             !< maximum PO4 sorption capacity for po4_assoc_slow, [mol m-3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout)  :: km_fast_po4               !< half-saturation PO4 concentration for po4_assoc_fast, [mol m-3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout)  :: km_slow_po4               !< half-saturation PO4 concentration for po4_assoc_slow, [mol m-3]
    !------------------------------------------------------------------------------------------------------ !
    INTEGER                               :: isoil                              ! looping
    REAL(wp)                              :: total_p_at_p_depth                 ! = soil_p_slow + soil_p_labile + soil_p_occluded + soil_p_primary [g P / m2]
    REAL(wp), DIMENSION(nc)               :: density_p_at_p_depth               ! density of P within the soil until 'soil_p_depth' (or 'soil_depth', whatever is lower) [g P / kg soil]
    REAL(wp), DIMENSION(nc)               :: soil_weight_at_p_depth             ! total soil weight up to the soil_p_depth, [kg]
    INTEGER,  DIMENSION(nc)               :: soil_layer_at_p_depth              ! up to soil_layer_at_p_depth use the data from the P init map info
    REAL(wp), DIMENSION(nc)               :: hlp1, hlp2, hlp3                   ! helper
    LOGICAL,  DIMENSION(nc)               :: l_calc_soil_properties_at_p_depth  ! flag for looping
    REAL(wp), DIMENSION(nc,nsoil_sb)      :: soil_weight_sl                     ! soil weight per soil layer
    REAL(wp), DIMENSION(nc,nsoil_sb)      :: hlp5
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_soil_p_pools'


    !>1.0 calculate the total soil weight up to the soil_p_depth, 
    !>and the soil weight per soil layer
    !>
    hlp1(:)                               = 0._wp
    hlp2(:)                               = 0._wp
    total_p_at_p_depth                    = soil_p_slow + soil_p_labile + soil_p_occluded + soil_p_primary
    l_calc_soil_properties_at_p_depth(:)  = .TRUE.

    ! if 'total soil depth until bedrock' < .OR. eq 'soil_p_depth'
    WHERE ((soil_depth(:) < soil_p_depth) .OR. (ABS((soil_depth(:) - soil_p_depth)) < eps8))
      ! calc soil weight until soil_depth (instead of until 'soil_p_depth' which is considered below for cases of 'total soil depth until bedrock' > 'soil_p_depth'
      soil_weight_at_p_depth(:) = (soil_depth(:) * SUM(bulk_dens_sl(:,:), DIM=2)) * (soil_depth(:) / soil_p_depth)
      ! set flag to F, as calculations are done for this grid point
      l_calc_soil_properties_at_p_depth(:) = .FALSE.
    END WHERE
    ! save the soil layer that includes the 'soil_depth'
    DO isoil = 1,nsoil_sb
      ! all grid points with 'total soil depth until bedrock' < .OR. eq 'soil_p_depth'
      WHERE (.NOT. l_calc_soil_properties_at_p_depth(:))
        ! if the width of the current soil layer is larger than zero
        ! then it can include the 'soil_depth'
        ! the last soil layer for which this is TRUE actually includes the 'soil_depth'
        WHERE (soil_lay_width_sl(:,isoil) > eps8)
          soil_layer_at_p_depth(:)  = isoil
        END WHERE
      END WHERE
    END DO

    ! if 'total soil depth until bedrock' > 'soil_p_depth'
    DO isoil = 1,nsoil_sb
      soil_weight_sl(:,isoil) = soil_lay_width_sl(:,isoil) * bulk_dens_sl(:,isoil)
      ! if the bottom of the soil layer (upper bound in vgrid) is within the 'soil_p_depth' (which is 0.5 m by default), i.e., l_calc_soil_properties_at_p_depth=T
      WHERE (l_calc_soil_properties_at_p_depth(:)) 
        ! sum up the width of each soil layer and the soil weight
        hlp1(:) = hlp1(:) + soil_lay_width_sl(:,isoil)
        hlp2(:) = hlp2(:) + soil_weight_sl(:,isoil)
        ! if the bottom of the soil layer is deeper than the 'soil_p_depth', set l_calc_soil_properties_at_p_depth=F
        WHERE ((hlp1(:) > soil_p_depth))
          ! calc soil weight until 'soil_p_depth'
          soil_weight_at_p_depth(:) = hlp2(:) - soil_weight_sl(:,isoil) * (hlp1(:) - soil_p_depth) / soil_lay_width_sl(:,isoil)
          ! save the soil layer that includes the 'soil_p_depth'
          soil_layer_at_p_depth(:)  = isoil
          ! set flag to F, as calculations are done for this grid point
          l_calc_soil_properties_at_p_depth(:) = .FALSE.
        END WHERE
      END WHERE
    END DO

    ! calc the P density in soil [g P / kg soil]
    WHERE (soil_depth(:) > soil_p_depth)
      ! until 'soil_p_depth'
      density_p_at_p_depth(:) = total_p_at_p_depth / soil_weight_at_p_depth(:)
    ELSEWHERE
      ! until 'soil_depth'
      density_p_at_p_depth(:) = total_p_at_p_depth * (soil_depth(:) / soil_p_depth) / soil_weight_at_p_depth(:)
    END WHERE

    !>2.0 calculate the initial slow P, occluded P, primary P pool, and p_labile
    !>
    !>  based on the assumption that soil P pool are linearly corrected with soil weight
    !>  the depth correction of initial P pool distributions are fitted against the observed pattern of a young cambisol soil 
    !>  with basalt as parent material (BBR, Lang et al. 2017) that, 
    !>  the proportion of: labile+sorbed P decreases with increasing depth, and primary P increases with increasing depth (Yu et al., 2020)
    !>
    DO isoil = 1,nsoil_sb
      WHERE (soil_lay_width_sl(:,isoil) > eps8)
        ! primary P, corrected with depth, [g P] 
        hlp1(:) = MIN(k2_p_pool_depth_corr * soil_lay_depth_center_sl(:,isoil) ** k4_p_pool_depth_corr, &
                      (k1_p_pool_depth_corr - k_p_pool_depth_corr_min_occlud)) / k1_p_pool_depth_corr
        hlp2(:) = k3_p_pool_depth_corr ** (k2_p_pool_depth_corr * soil_lay_depth_center_sl(:,isoil)) / k1_p_pool_depth_corr
        ! occluded P, corrected with depth
        hlp3(:) = MAX((1.0_wp - hlp1(:) - hlp2(:)), k_p_pool_depth_corr_min_occlud / k1_p_pool_depth_corr) 

        WHERE ((hlp1(:) + hlp2(:) + hlp3(:)) > 1._wp)
          hlp1(:) = hlp1(:) / (hlp1(:) + hlp2(:) + hlp3(:))
          hlp2(:) = hlp2(:) / (hlp1(:) + hlp2(:) + hlp3(:))
          hlp3(:) = hlp3(:) / (hlp1(:) + hlp2(:) + hlp3(:))
        END WHERE
        WHERE (soil_layer_at_p_depth(:) > isoil) 
          po4_assoc_slow(:,isoil) = bulk_dens_sl(:,isoil) * density_p_at_p_depth(:) * (soil_p_slow + soil_p_labile) / &
            & total_p_at_p_depth / molar_mass_P * (1._wp - p_labile_slow_ratio(:,isoil))
          po4_occluded(:,isoil)   = bulk_dens_sl(:,isoil) * density_p_at_p_depth(:) * soil_p_occluded / &
            & total_p_at_p_depth / molar_mass_P
          po4_primary(:,isoil)    = bulk_dens_sl(:,isoil) * density_p_at_p_depth(:) * soil_p_primary / &
            & total_p_at_p_depth / molar_mass_P
          hlp5(:,isoil)           = bulk_dens_sl(:,isoil) * density_p_at_p_depth(:) * (soil_p_slow + soil_p_labile) / &
            & total_p_at_p_depth / molar_mass_P * p_labile_slow_ratio(:,isoil) 
        ELSEWHERE
          po4_assoc_slow(:,isoil) = bulk_dens_sl(:,isoil) * density_p_at_p_depth(:) * hlp2(:) / molar_mass_P * &
            &                       (1._wp - p_labile_slow_ratio(:,isoil))
          po4_occluded(:,isoil)   = bulk_dens_sl(:,isoil) * density_p_at_p_depth(:) * hlp3(:) / molar_mass_P
          po4_primary(:,isoil)    = bulk_dens_sl(:,isoil) * density_p_at_p_depth(:) * hlp1(:) / molar_mass_P
          hlp5(:,isoil)           = bulk_dens_sl(:,isoil) * density_p_at_p_depth(:) * hlp2(:) / molar_mass_P * &
            &                       p_labile_slow_ratio(:,isoil) 
        END WHERE
      END WHERE ! soil_lay_width_sl(:,isoil) > eps8
    END DO

    !>  2.1 calculate the initial solute P and fast P,
    !>
    !>  based on the assumption that labile P is the sum of both,
    !>  and they follow the Langmuir equilibrium
    !>  po4_assoc_fast = (p_labile_sl + Smax + K - sqrt((p_labile_sl + Smax + K)**2._wp - 4*p_labile_sl*Smax)) / 2._wp
    !>
    IF (flag_sb_double_langmuir) THEN
      DO isoil = 1,nsoil_sb
        WHERE (soil_lay_width_sl(:,isoil) > eps8)
          po4_assoc_fast(:,isoil) = calc_fast_po4(hlp5(:,isoil), qmax_fast_po4(:,isoil), km_fast_po4(:,isoil))
          po4_solute(:,isoil)     = hlp5(:,isoil) - po4_assoc_fast(:,isoil)
        END WHERE ! soil_lay_width_sl(:,isoil) > eps8
      END DO
    ELSE
      DO isoil = 1,nsoil_sb
        WHERE (soil_lay_width_sl(:,isoil) > eps8)
          po4_assoc_fast(:,isoil) = calc_fast_po4(hlp5(:,isoil), qmax_po4(:,isoil), km_po4(:,isoil))
          po4_solute(:,isoil)     = hlp5(:,isoil) - po4_assoc_fast(:,isoil)
        END WHERE ! soil_lay_width_sl(:,isoil) > eps8
      END DO
    ENDIF

  END SUBROUTINE init_soil_p_pools

  ! !-----------------------------------------------------------------------------------------------------
  ! !> init soil properties for SPP sites only
  ! !! 
  ! !! Output: ph and Qmax_AlFe_cor prescribed for spp sites
  ! !!
  ! !! Qmax_AlFe_cor is further modified with calc_Psorption_parameter()
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE init_soil_properties_spp_sites(nc, nsoil_sb, site_ID_spp1685, ph_sl, Qmax_AlFe_cor)
  ! 
  !   IMPLICIT NONE
  ! 
  !   INTEGER,                                          INTENT(in)    :: nc, &                        !< dimensions
  !                                                                      nsoil_sb
  !   CHARACTER(len=*), OPTIONAL,                       INTENT(in)    :: site_ID_spp1685              !< site ID for spp1685  
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(out)   :: ph_sl                        !< ph of soil layer [-]
  !   REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(out)   :: Qmax_AlFe_cor                !< qmax correlation factor of Al/Fe content [-]  
  ! 
  !   INTEGER  :: i
  ! 
  !   CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_soil_properties_spp_sites'
  ! 
  !   ! ------------------------------------------------------------------------------------------------------------
  !   ! Go science
  !   ! ------------------------------------------------------------------------------------------------------------
  ! 
  !   !> 1.0 prescribe ph and Qmax_AlFe_cor from measurements
  !   !! 
  !   SELECT CASE(site_ID_spp1685)
  !   CASE ("BBR")
  !     DO i=1,nc
  !       ph_sl(i,1:nsoil_sb)         = (/3.84_wp, 4.14_wp, 4.45_wp, 4.66_wp, 4.85_wp, &
  !                                       4.92_wp,5.19_wp,5.64_wp,5.8_wp,6._wp, &
  !                                       6.2_wp, 6.5_wp, 6.5_wp, 6.5_wp, 6.5_wp/)
  !       Qmax_AlFe_cor(i,1:nsoil_sb) = (/3.04_wp, 3.60_wp, 3.60_wp, 3.295_wp, 2.92_wp, &
  !                                       2.75_wp,2.40_wp,2.04_wp,2._wp,2._wp, &
  !                                       2._wp, 2._wp, 2._wp, 2._wp, 2._wp/)
  !     ENDDO
  !   CASE ("LUE")
  !     DO i=1,nc
  !       ph_sl(i,1:nsoil_sb)         = (/3.51_wp,3.78_wp,3.84_wp,4.14_wp,4.20_wp, &
  !                                       4.46_wp,4.45_wp,4.62_wp,4.63_wp,4.8_wp, &
  !                                       5._wp, 5.2_wp, 5.5_wp, 6._wp, 6.5_wp/)
  !       Qmax_AlFe_cor(i,1:nsoil_sb) = (/0.17_wp, 0.29_wp, 0.41_wp, 0.41_wp, 0.396_wp, &
  !                                       0.3_wp,0.289_wp,0.236_wp,0.2_wp,0.2_wp, &
  !                                       0.2_wp, 0.2_wp, 0.2_wp, 0.2_wp, 0.2_wp/)
  !     ENDDO
  !   CASE ("VES")
  !     DO i=1,nc
  !       ph_sl(i,1:nsoil_sb)         = (/3.47_wp,3.91_wp,4.25_wp,4.46_wp,4.52_wp, &
  !                                       4.66_wp,4.78_wp,4.77_wp,5._wp,5.2_wp, &
  !                                       5.5_wp, 5.8_wp, 6.0_wp, 6.5_wp, 6.5_wp/)
  !       Qmax_AlFe_cor(i,1:nsoil_sb) = (/1.275_wp, 1.515_wp, 1.378_wp, 1.459_wp, 1.90_wp, &
  !                                       2.469_wp,3.428_wp,5.33_wp,6._wp,6._wp, &
  !                                       6._wp, 6._wp, 6._wp, 6._wp, 6._wp/)
  !     ENDDO
  !   CASE ("MIT")
  !     DO i=1,nc
  !       ph_sl(i,1:nsoil_sb)         = (/3.64_wp,4.02_wp,4.53_wp,4.46_wp,4.46_wp, &
  !                                       4.58_wp,4.56_wp,4.64_wp,4.7_wp,4.8_wp, &
  !                                       5._wp, 5.2_wp, 5.5_wp, 5.8_wp, 6._wp/)
  !       Qmax_AlFe_cor(i,1:nsoil_sb) = (/1.326_wp, 1.472_wp, 1.573_wp, 1.5_wp, 1.426_wp, &
  !                                       1.3_wp,1.258_wp,1.2_wp,1.2_wp,1.2_wp, &
  !                                       1.1_wp, 1.1_wp, 1._wp, 1._wp, 1._wp/)
  !     ENDDO
  !   CASE ("CON")
  !     DO i=1,nc
  !       ph_sl(i,1:nsoil_sb)         = (/4.03_wp,4.12_wp,4.32_wp,4.37_wp,4.42_wp, &
  !                                       4.54_wp,4.60_wp,4.61_wp,4.7_wp,4.8_wp, &
  !                                       5._wp, 5.2_wp, 5.5_wp, 5.8_wp, 6._wp/)
  !       Qmax_AlFe_cor(i,1:nsoil_sb) = (/1.9_wp, 1.891_wp, 1.87_wp, 1.8_wp, 1.73_wp, &
  !                                       1.588_wp,1.436_wp,1.24_wp,1.2_wp,1.2_wp, &
  !                                       1.1_wp, 1.1_wp, 1._wp, 1._wp, 1._wp/)
  !     ENDDO
  !   END SELECT
  ! 
  ! END SUBROUTINE init_soil_properties_spp_sites
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! ! SPP site-set
  ! !
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Intialize soil_biogeochemistry process specifically for all SPP sites (after memory has been set up)
  ! !! site of the SPP site-set: BBR, MIT, VES, CON, LUE
  ! !! 
  ! !! @TODO this routine works only at one point! because soil properties for exactly one site are selected
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE init_soil_p_pools_spp_sites (nc, nsoil_sb, soil_lay_depth_center_sl, &
  !                                       soil_p_labile, soil_p_slow, soil_p_occluded, soil_p_primary, &
  !                                       sitename, p_labile_slow_ratio, &
  !                                       p_labile, po4_assoc_slow, po4_occluded, po4_primary)
  ! 
  !   USE mo_sb_constants !,         ONLY: 
  ! 
  !   IMPLICIT NONE
  ! 
  !   INTEGER                         , INTENT(in)     :: nc                        !< dimensions
  !   INTEGER                         , INTENT(in)     :: nsoil_sb                  !< number of soil layers as used/defined by the SB_ process
  !   REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)     :: soil_lay_depth_center_sl  !< depth at the center of each soil layer
  !   REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)     :: p_labile_slow_ratio   !< ratio between labile and labile+slow P pools
  !   REAL(wp)                        , INTENT(in)     :: soil_p_labile, &      !< labile P from Yang et al. 2013, g/m2
  !                                                       soil_p_slow, &        !< secondary P from Yang et al. 2013, g/m2
  !                                                       soil_p_occluded, &    !< Occluded P from Yang et al. 2013, g/m2
  !                                                       soil_p_primary        !< Primary from Yang et al. 2013, g/m2
  !   CHARACTER(len=*)                , INTENT(in)     :: sitename              !< spp site name: BBR or LUE (MIT, VES, CON)
  !   REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(inout)  :: p_labile, &           !< PO4 in solution per soil layer [mol / m3]
  !                                                       po4_assoc_slow, &     !< slow minerally associated PO4 pool per soil layer [mol/m3]
  !                                                       po4_occluded, &       !< occluded PO4 pool per soil layer [mol/m3]
  !                                                       po4_primary           !< primary PO4 pool per soil layer [mol/m3]
  ! 
  !   INTEGER                                      :: isoil, ichunk
  !   REAL(wp)                                     :: hlp1, hlp2, hlp3, hlp4, hlp5, hlp6
  !   REAL(wp), DIMENSION(nc)                      :: arr_hlp2, arr_hlp3, arr_hlp4, soil_layer_at_p_depth
  !   REAL(wp), DIMENSION(nsoil_sb)                :: arr_hlp5
  !   CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_soil_p_pools_spp_sites'
  ! 
  !   !> 1.0 set site-level soil properties
  !   !!
  !   SELECT CASE (sitename)
  !     CASE("BBR")
  !       arr_hlp5(:) = (/9.24_wp,36.92_wp,39.44_wp,40.36_wp,42.75_wp, &  
  !                      45.72_wp,32.45_wp,95.89_wp,124.14_wp,133.93_wp, &
  !                      100._wp,126._wp,122.18_wp,145._wp,90._wp/)       ! observed total SIP, mol P m-3
  !       hlp1            = soil_p_labile     ! p_labile from Yang et al. 2013, g/m2
  !       hlp2            = soil_p_slow       ! p_slow from Yang et al. 2013, g/m2
  !       hlp3            = soil_p_primary    ! p_occluded from Yang et al. 2013, g/m2
  !       hlp4            = soil_p_occluded   ! p_primary from Yang et al. 2013, g/m2
  !       soil_layer_at_p_depth(:) = 6._wp             ! field measurement depth corresponded model soil layers
  !       ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
  !       !hlp5        = 1750._wp / 1000._wp ! g P / kg soil
  !       !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
  !     CASE("LUE")
  !       arr_hlp5(:) = (/0.8867_wp,3.187_wp,2.662_wp,4.887_wp,4.792_wp, &  
  !                      2.987_wp,3.313_wp,6.0616_wp,7.56_wp,5.787_wp, &
  !                      5.932_wp,6.845_wp,8.294_wp,8.413_wp,8.532_wp/)   ! observed total SIP, mol P m-3
  !       hlp1            = soil_p_labile     ! p_labile from Yang et al. 2013, g/m2
  !       hlp2            = soil_p_slow       ! p_slow from Yang et al. 2013, g/m2
  !       hlp3            = soil_p_occluded   ! p_occluded from Yang et al. 2013, g/m2
  !       hlp4            = soil_p_primary    ! p_primary from Yang et al. 2013, g/m2
  !       soil_layer_at_p_depth(:) = 6._wp             ! field measurement depth corresponded model soil layers
  !       ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
  !       !hlp5        = 120._wp / 1000._wp ! g P / kg soil
  !       !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
  !     CASE("VES")
  !       arr_hlp5(:) = (/2.388_wp,8.509_wp,9.033_wp,10.811_wp,15.365_wp, &  
  !                      20.661_wp,22.721_wp,29.321_wp,30.952_wp,32.821_wp, &
  !                      40.673_wp,41.880_wp,41.111_wp,48.438_wp,54._wp/)   ! observed total SIP, mol P m-3
  !       hlp1            = soil_p_labile     ! p_labile from Yang et al. 2013, g/m2
  !       hlp2            = soil_p_slow       ! p_slow from Yang et al. 2013, g/m2
  !       hlp3            = soil_p_primary    ! p_occluded from Yang et al. 2013, g/m2
  !       hlp4            = soil_p_occluded   ! p_primary from Yang et al. 2013, g/m2
  !       soil_layer_at_p_depth(:) = 6._wp             ! field measurement depth corresponded model soil layers
  !       ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
  !       !hlp5        = 900._wp / 1000._wp ! g P / kg soil
  !       !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
  !     CASE("MIT")
  !       arr_hlp5(:) = (/3.723_wp,6.272_wp,7.556_wp,8.292_wp,6.93_wp, &  
  !                      13.731_wp,13.296_wp,16.511_wp,12.890_wp,15._wp, &
  !                      16.667_wp,16.875_wp,17.5_wp,15.104_wp,15.625_wp/)  ! observed total SIP, mol P m-3
  !       hlp1            = soil_p_labile     ! p_labile from Yang et al. 2013, g/m2
  !       hlp2            = soil_p_slow       ! p_slow from Yang et al. 2013, g/m2
  !       hlp3            = soil_p_primary    ! p_occluded from Yang et al. 2013, g/m2
  !       hlp4            = soil_p_occluded   ! p_primary from Yang et al. 2013, g/m2
  !       soil_layer_at_p_depth(:) = 9._wp             ! field measurement depth corresponded model soil layers
  !       ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
  !       !hlp5        = 900._wp / 1000._wp ! g P / kg soil
  !       !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
  !     CASE("CON")
  !       arr_hlp5(:) = (/1.862_wp,2.924_wp,4.027_wp,4.048_wp,4.069_wp, &  
  !                      7.167_wp,10.855_wp,7.719_wp,18.75_wp,18.33_wp, &
  !                      23.44_wp,25._wp,21.093_wp,21.875_wp,15.104_wp/)   ! observed total SIP, mol P m-3
  !       hlp1            = soil_p_labile     ! p_labile from Yang et al. 2013, g/m2
  !       hlp2            = soil_p_slow       ! p_slow from Yang et al. 2013, g/m2
  !       hlp3            = soil_p_primary    ! p_occluded from Yang et al. 2013, g/m2
  !       hlp4            = soil_p_occluded   ! p_primary from Yang et al. 2013, g/m2
  !       soil_layer_at_p_depth(:) = 9._wp             ! field measurement depth corresponded model soil layers
  !       ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
  !       !hlp5        = 900._wp / 1000._wp ! g P / kg soil
  !       !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
  !   END SELECT
  ! 
  !   !> 2.0 calculate the initial slow P, occluded P, primary P pool, and p_labile
  !   !! @TODO this is inconsistent with section 1.0 of this routine, as section 1.0 does not work for multiple points in a chunk
  !   !! based on the assumption that soil P pool are linearly corrected with soil weight
  !   !! the depth correction of initial P pool distributions are fitted against the observed pattern of a young cambisol soil 
  !   !! with basalt as parent material (BBR, Lang et al. 2017) that, 
  !   !! the proportion of: labile+sorbed P decreases with increasing depth, and primary P increases with increasing depth (Yu et al., 2020)
  !   DO isoil = 1, nsoil_sb
  !     ! primary corrected for depth, [g P] 
  !     arr_hlp4(:) = MIN(k2_p_pool_depth_corr * soil_lay_depth_center_sl(:,isoil) ** k4_p_pool_depth_corr, &
  !                   (k1_p_pool_depth_corr - k_p_pool_depth_corr_min_occlud)) / k1_p_pool_depth_corr
  !     ! labile + slow
  !     arr_hlp2(:) = k3_p_pool_depth_corr ** (k2_p_pool_depth_corr * soil_lay_depth_center_sl(:,isoil)) / k1_p_pool_depth_corr
  !     ! occluded corrected for depth
  !     arr_hlp3(:) = MAX((1.0_wp - arr_hlp4(:) - arr_hlp2(:)), k_p_pool_depth_corr_min_occlud / k1_p_pool_depth_corr) 
  ! 
  !     WHERE ((arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:)) > 1._wp)
  !       arr_hlp4(:) = arr_hlp4(:) / (arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:))
  !       arr_hlp2(:) = arr_hlp2(:) / (arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:))
  !       arr_hlp3(:) = arr_hlp3(:) / (arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:))
  !     ENDWHERE
  ! 
  !     WHERE (isoil < soil_layer_at_p_depth(:)) 
  !       po4_assoc_slow(:,isoil) = arr_hlp5(isoil) * (hlp1 + hlp2) / (hlp1 + hlp2 + hlp3 + hlp4) * &
  !                                 (1._wp - p_labile_slow_ratio(:,isoil))
  !       po4_occluded(:,isoil)   = arr_hlp5(isoil) * hlp3 / (hlp1 + hlp2 + hlp3 + hlp4)
  !       po4_primary(:,isoil)    = arr_hlp5(isoil) * hlp4 / (hlp1 + hlp2 + hlp3 + hlp4)
  !       p_labile(:,isoil)       = arr_hlp5(isoil) * (hlp1 + hlp2) / (hlp1 + hlp2 + hlp3 + hlp4) * p_labile_slow_ratio(:,isoil)
  !     ELSEWHERE
  !       po4_assoc_slow(:,isoil) = arr_hlp5(isoil) * arr_hlp2(:) * (1._wp - p_labile_slow_ratio(:,isoil))
  !       po4_occluded(:,isoil)   = arr_hlp5(isoil) * arr_hlp3(:)
  !       po4_primary(:,isoil)    = arr_hlp5(isoil) * arr_hlp4(:)
  !       p_labile(:,isoil)       = arr_hlp5(isoil) * arr_hlp2(:) * p_labile_slow_ratio(:,isoil)
  !     ENDWHERE
  !   ENDDO
  ! 
  ! END SUBROUTINE init_soil_p_pools_spp_sites

#endif
END MODULE mo_sb_init
