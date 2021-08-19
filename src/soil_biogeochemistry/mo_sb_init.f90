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
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version for SLM by Jan Engel (20161208) based on the JBACH4 version by Reiner Schnur (2013-04-04)
!!
MODULE mo_sb_init

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish,message
  USE mo_jsb_control,           ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sb_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize soil_biogeochemistry process (after memory has been set up)
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE sb_init(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_atmland_constants,     ONLY: def_co2_mixing_ratio, def_co2_mixing_ratio_C13, def_co2_mixing_ratio_C14, &
                                        def_co2_deltaC13, def_co2_deltaC14
    USE mo_isotope_util,          ONLY: calc_mixing_ratio_C13C12, calc_mixing_ratio_C14C, calc_mixing_ratio_N15N14, &
                                        calc_fractionation
    USE mo_jsb_phys_constants,              ONLY: molar_mass_C, molar_mass_N, molar_mass_P
    USE mo_jsb_math_constants,              ONLY: eps4, eps8
    USE mo_assimi_constants,                ONLY: CiCa_default_C3, CiCa_default_C4, ic3phot, ic4phot
    USE mo_assimi_update_canopy_fluxes,     ONLY: discrimination_ps
    USE mo_sb_jsm_processes                 ! e.g. calc_bulk_soil_carbon(), calc_Psorption_parameter()
    USE mo_sse_util,                        ONLY: calc_qmax_texture
    USE mo_sb_jsm_transport                 ! e.g. calc_bioturbation_rate()
    USE mo_veg_constants,                   ONLY: itree
    USE mo_sb_constants                     ! e.g.  nh4_solute_prescribe, no3_solute_prescribe, po4_solute_prescribe

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)

    ! USE PROC_memory
    dsl4jsb_Use_memory(SB_)
    dsl4slm_Use_memory_a_sb
    USE mo_sb_memory_class,       ONLY: sb_pool
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(VEG_)
    dsl4slm_Use_memory_a_veg
    USE mo_veg_memory_class,      ONLY: veg_pool

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib        !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER              :: vgrid_soil_sb !< Vertical grid
    INTEGER                                      :: nsoil_sb      !< number of soil layers as used/defined by the SB_ process
    REAL(wp)                                     :: site_longitude
    LOGICAL                                      :: flag_spp1685
    CHARACTER(len=3)                             :: site_ID_spp1685
    INTEGER                                      :: isoil, ichunk
    INTEGER                                      :: iblk, ics, ice, nc
    REAL(wp)                                     :: hlp1
    REAL(wp)                                     :: polyLit_cn_ratio, woodyLit_cn_ratio
    REAL(wp)                                     :: frac_c13, frac_c14
    REAL(wp), ALLOCATABLE, DIMENSION(:)          :: arr_hlp6, root_depth_avg          ! dim: nc
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)        :: arr_hlp7, arr_hlp8, arr_hlp9
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)        :: p_labile, p_labile_slow_site
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init'
    ! ---------------------------
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(VEG_)
    ! Declare pointers to variables in memory
    ! VEG_
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    ! SSE_ 3D
    dsl4jsb_Real3D_onChunk :: bulk_dens_sl
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    dsl4jsb_Real3D_onChunk :: soil_depth_sl       ! (soil_depth_sl == soil_lay_width_sl)
    dsl4jsb_Real3D_onChunk :: soil_lay_width_sl
    dsl4jsb_Real3D_onChunk :: soil_lay_depth_center_sl
    dsl4jsb_Real3D_onChunk :: soil_lay_depth_ubound_sl
    dsl4jsb_Real3D_onChunk :: soil_lay_depth_lbound_sl
    dsl4jsb_Real3D_onChunk :: qmax_org_min_sl
    dsl4jsb_Real3D_onChunk :: qmax_po4_min_sl
    dsl4jsb_Real3D_onChunk :: qmax_nh4_min_sl
    dsl4jsb_Real3D_onChunk :: qmax_po4_om_sl
    dsl4jsb_Real3D_onChunk :: sand_sl
    dsl4jsb_Real3D_onChunk :: silt_sl
    dsl4jsb_Real3D_onChunk :: clay_sl
    dsl4jsb_Real3D_onChunk :: volume_min_sl
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk :: nh4_assoc
    dsl4jsb_Real3D_onChunk :: nh4_solute
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute
    dsl4jsb_Real3D_onChunk :: noy
    dsl4jsb_Real3D_onChunk :: n2o
    dsl4jsb_Real3D_onChunk :: n2
    dsl4jsb_Real3D_onChunk :: nh4_n15_assoc
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute
    dsl4jsb_Real3D_onChunk :: no3_n15_solute
    dsl4jsb_Real3D_onChunk :: bulk_soil_carbon_sl
    dsl4jsb_Real3D_onChunk :: soil_litter_carbon_sl
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onChunk :: qmax_org
    dsl4jsb_Real3D_onChunk :: qmax_po4
    dsl4jsb_Real3D_onChunk :: qmax_nh4
    dsl4jsb_Real3D_onChunk :: qmax_fast_po4
    dsl4jsb_Real3D_onChunk :: qmax_slow_po4
    dsl4jsb_Real3D_onChunk :: km_fast_po4
    dsl4jsb_Real3D_onChunk :: km_slow_po4
    dsl4jsb_Real3D_onChunk :: km_adsorpt_po4_sl
    dsl4jsb_Real3D_onChunk :: km_adsorpt_nh4_sl
    dsl4jsb_Real3D_onChunk :: k_bioturb
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast
    dsl4jsb_Real3D_onChunk :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: po4_occluded
    dsl4jsb_Real3D_onChunk :: ph_sl
    dsl4jsb_Real3D_onChunk :: Qmax_AlFe_cor
    dsl4jsb_Real3D_onChunk :: lithology_class_sl
    dsl4jsb_Real3D_onChunk :: po4_primary
    dsl4jsb_Real3D_onChunk :: microbial_cue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk :: microbial_nue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk :: microbial_pue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk :: vmax_weath_mineral_sl
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly
    dsl4jsb_Real3D_onChunk :: enzyme_frac_residue
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_c
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_n
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_p
    dsl4jsb_Real3D_onChunk :: dom_cn
    dsl4jsb_Real3D_onChunk :: dom_cp
    ! Get local variables from options argument
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(SB_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts%lib_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    site_longitude = model%config%lon
    flag_spp1685  = model%config%flag_spp1685
    ! Get process config
    dsl4jsb_Get_config(SB_)
    ! Get process memories
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(VEG_)
    ! Allocate local allocatable variables (deallocate below at the end of this routine)    
    ALLOCATE(arr_hlp6(nc))
    ALLOCATE(root_depth_avg(nc))
    ALLOCATE(arr_hlp7(nc, nsoil_sb))
    ALLOCATE(arr_hlp8(nc, nsoil_sb))
    ALLOCATE(arr_hlp9(nc, nsoil_sb))
    ALLOCATE(p_labile(nc, nsoil_sb))
    ALLOCATE(p_labile_slow_site(nc, nsoil_sb))
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var3D_onChunk(VEG_,   root_fraction_sl)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SSE_,   bulk_dens_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,   w_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,   soil_depth_sl)            ! in  (soil_depth_sl == soil_lay_width_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,   soil_lay_width_sl)        ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,   soil_lay_depth_center_sl) ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,   soil_lay_depth_ubound_sl) ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,   soil_lay_depth_lbound_sl) ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,   qmax_org_min_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,   qmax_po4_min_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,   qmax_nh4_min_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,   qmax_po4_om_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,   sand_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,   silt_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,   clay_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,   volume_min_sl)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SB_,    nh4_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_,    nh4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,    no3_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,    po4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,    noy)
    dsl4jsb_Get_var3D_onChunk(SB_,    n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,    n2)
    dsl4jsb_Get_var3D_onChunk(SB_,    nh4_n15_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_,    nh4_n15_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,    no3_n15_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,    bulk_soil_carbon_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,    soil_litter_carbon_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,    bulk_dens_corr_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,    qmax_org)
    dsl4jsb_Get_var3D_onChunk(SB_,    qmax_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,    qmax_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,    qmax_fast_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,    qmax_slow_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,    km_fast_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,    km_slow_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,    km_adsorpt_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,    km_adsorpt_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,    k_bioturb)
    dsl4jsb_Get_var3D_onChunk(SB_,    po4_assoc_fast)
    dsl4jsb_Get_var3D_onChunk(SB_,    po4_assoc_slow)
    dsl4jsb_Get_var3D_onChunk(SB_,    po4_occluded) 
    dsl4jsb_Get_var3D_onChunk(SB_,    ph_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,    Qmax_AlFe_cor)
    dsl4jsb_Get_var3D_onChunk(SB_,    lithology_class_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,    po4_primary)
    dsl4jsb_Get_var3D_onChunk(SB_,    microbial_cue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,    microbial_nue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,    microbial_pue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,    vmax_weath_mineral_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,    enzyme_frac_poly)
    dsl4jsb_Get_var3D_onChunk(SB_,    enzyme_frac_residue)
    dsl4jsb_Get_var3D_onChunk(SB_,    enzyme_frac_poly_c)
    dsl4jsb_Get_var3D_onChunk(SB_,    enzyme_frac_poly_n)
    dsl4jsb_Get_var3D_onChunk(SB_,    enzyme_frac_poly_p)
    dsl4jsb_Get_var3D_onChunk(SB_,    dom_cn)
    dsl4jsb_Get_var3D_onChunk(SB_,    dom_cp)
    ! Get process memory pool structures
    dsl4slm_get_vegpool_var(pool,fine_root,carbon)      ! needs veg_pool
#define Get_sb_pool
#include "slm_pools_a_incl.f90"
#undef Get_sb_pool

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 1.0 all soil models
    !! 

    !> 1.1 calc frac_c13 & frac_c14
    !! 
    IF(lctlib%ps_pathway == ic3phot)THEN
       hlp1 = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
    ELSE
       hlp1 = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
    ENDIF
    frac_c13 = calc_fractionation(def_co2_mixing_ratio, def_co2_mixing_ratio_C13, hlp1)

    IF(lctlib%ps_pathway == ic3phot)THEN
       hlp1 = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
    ELSE
       hlp1 = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
    ENDIF
    frac_c14 = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14-hlp1)

    !> 1.2 init site-specific soil properties
    !! 
    ph_sl(:,:)              = dsl4jsb_Config(SB_)%soil_ph
    lithology_class_sl(:,:) = dsl4jsb_Config(SB_)%lithology_class
    Qmax_AlFe_cor(:,:)      = 1._wp

    !> 1.3 init SB_ variables from constants
    !! 
    microbial_cue_eff_tmic_mavg(:,:)  = microbial_cue_max
    microbial_nue_eff_tmic_mavg(:,:)  = microbial_nue
    microbial_pue_eff_tmic_mavg(:,:)  = microbial_pue
    vmax_weath_mineral_sl(:,:)        = k_weath_mineral 
    dom_cn(:,:)                       = microbial_cn
    dom_cp(:,:)                       = microbial_cn * microbial_np

    !------------------------------------------------------------------------------------------------------------------------------
    !> 1.4 SPP site-set: identify site based on longitude
    !! five sites from the "SPP 1685 Project: Forest Strategies for limited Phosphorus Resources"
    site_ID_spp1685 = "OTH"   ! default: other
    IF (flag_spp1685) THEN
      IF (ABS(site_longitude - 9.75_wp) < eps8) THEN
        site_ID_spp1685 = "BBR"
      ELSEIF (ABS(site_longitude - 10.25_wp) < eps8) THEN
        site_ID_spp1685 = "LUE"
      ELSEIF (ABS(site_longitude - 10.75_wp) < eps8) THEN
        site_ID_spp1685 = "VES"
      ELSEIF (ABS(site_longitude - 12.75_wp) < eps8) THEN
        site_ID_spp1685 = "MIT"
      ELSEIF (ABS(site_longitude - 7.75_wp)  < eps8) THEN
        site_ID_spp1685 = "CON"
      ENDIF
    ENDIF
    !------------------------------------------------------------------------------------------------------------------------------

    !> 2.0 soil model init
    !! i)  simple_1d
    !! ii) JSM
    SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
       !!==========================================================================================================================
       !> simple_1d soil models
       !! a four pool soil biogeochemical model
       CASE("simple_1d")
          sb_pool%soluable_litter%carbon(:,1)  = def_sb_pool_metabol_litter /  molar_mass_C / soil_lay_width_sl(:,1) * &
                                                 (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1)) 
          sb_pool%polymeric_litter%carbon(:,1) = def_sb_pool_simple_struct_litter / molar_mass_C / soil_lay_width_sl(:,1) * &
                                                 (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1)) 
          IF(lctlib%growthform==itree) THEN
            sb_pool%woody_litter%carbon(:,1) = def_sb_pool_woody_litter / molar_mass_C / soil_lay_width_sl(:,1)
          ENDIF
          sb_pool%microbial%carbon(:,1) = def_sb_pool_fast_som / molar_mass_C / soil_lay_width_sl(:,1) * &
                                          (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1)) 
          sb_pool%residue%carbon(:,1)   = def_sb_pool_slow_som / molar_mass_C / soil_lay_width_sl(:,1) * & 
                                          (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1)) 
          DO isoil=2,nsoil_sb
            sb_pool%soluable_litter%carbon(:,isoil)  = def_sb_pool_metabol_litter /  molar_mass_C / &
                                                       soil_lay_width_sl(:,isoil) * &
                                                       (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil) 
            sb_pool%polymeric_litter%carbon(:,isoil) = def_sb_pool_simple_struct_litter / molar_mass_C / &
                                                       soil_lay_width_sl(:,isoil) * &
                                                       (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil) 
            sb_pool%microbial%carbon(:,isoil)        = def_sb_pool_fast_som / molar_mass_C / soil_lay_width_sl(:,isoil) * &
                                                       (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil) 
            sb_pool%residue%carbon(:,isoil)          = def_sb_pool_slow_som / molar_mass_C / soil_lay_width_sl(:,isoil) * & 
                                                       (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil)
          ENDDO
          ! note: this differs on purpose from the calculation in mo_sb_update_pools for the prescribed option
          !       this is because here I'm trying to ensure that the models are initialised with the same amount of nutrient
          !       whereas in the prescribed option the target is unlimiting concentrations across the profile. The latter leads
          !       to more nutrient in the soil column in 1d compared to 0d.
          DO ichunk = 1,nc
            nh4_solute(ichunk,:) = nh4_solute_prescribe / soil_lay_depth_ubound_sl(ichunk,nsoil_sb)
            no3_solute(ichunk,:) = no3_solute_prescribe / soil_lay_depth_ubound_sl(ichunk,nsoil_sb)
          ENDDO

          ! inital bulk density and qmax_po4 based on current soil carbon content
          CALL calc_bulk_soil_carbon(nc, nsoil_sb, sb_pool, bulk_soil_carbon_sl, soil_litter_carbon_sl, &
                                        volume_min_sl, TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
          bulk_dens_corr_sl(:,:) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:),soil_litter_carbon_sl(:,:), &
                                                                   bulk_dens_sl(:,:))  
          qmax_org_min_sl(:,:)   = calc_qmax_texture(dsl4jsb_Config(SB_)%qmax_org_fine_particle, silt_sl(:,:), clay_sl(:,:)) 
          CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), bulk_dens_sl(:,:), 0._wp, &
                                                 qmax_org_min_sl(:,:), qmax_org(:,:))                                                                  
          !> po4 adsorption parameter intialisation 
          !!
          arr_hlp7(:,:) = po4_solute_prescribe
          ! specific init values for SPP sites for: ph_sl(:,:) & Qmax_AlFe_cor(:,:)
          IF (flag_spp1685) THEN
            CALL init_soil_properties_spp_sites(nc, nsoil_sb, site_ID_spp1685, ph_sl, Qmax_AlFe_cor)
          ENDIF
          ! Correct the first layer Qmax_AlFe_cor with the woody litter
          Qmax_AlFe_cor(:,1) = Qmax_AlFe_cor(:,1) * 1._wp / ((1._wp - volume_min_sl(:,1)) * &
                               (bulk_soil_carbon_sl(:,1) + soil_litter_carbon_sl(:,1)  + sb_pool%woody_litter%carbon(:,1)) / &
                               (bulk_soil_carbon_sl(:,1) + soil_litter_carbon_sl(:,1)) + volume_min_sl(:,1))
          CALL calc_Psorption_parameter(dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
                                       clay_sl(:,:), silt_sl(:,:), w_soil_sl(:,:) / soil_depth_sl(:,:) * 1000._wp, ph_sl(:,:), &
                                       bulk_soil_carbon_sl(:,:), soil_litter_carbon_sl(:,:), sb_pool%woody_litter%carbon(:,:), &
                                       bulk_dens_sl(:,:), volume_min_sl(:,:), arr_hlp7(:,:), &
                                       qmax_po4_min_sl(:,:), qmax_po4_om_sl(:,:), qmax_po4(:,:), km_adsorpt_po4_sl(:,:), &
                                       qmax_fast_po4(:,:), qmax_slow_po4(:,:), km_fast_po4(:,:), km_slow_po4(:,:), &
                                       Qmax_AlFe_cor(:,:))
          ! calculate the partition coefficient between labile and slow P pool based on 
          ! their Qmax (double-Langmuir) or global average values (single-Langmuir)
          IF (dsl4jsb_Config(SB_)%flag_sb_double_langmuir) THEN
            WHERE (qmax_po4(:,:) > eps8)
              p_labile_slow_site(:,:)   = MAX(qmax_fast_po4(:,:) / qmax_po4(:,:), eps4)
            ELSEWHERE
              p_labile_slow_site(:,:)   = eps4
            ENDWHERE
          ELSE
            p_labile_slow_site(:,:)     = p_labile_slow_global_avg
          ENDIF

          !> nh4 adsorption parameter intialisation 
          qmax_nh4_min_sl(:,:)     = calc_qmax_texture(qmax_nh4_clay, 0.0_wp, clay_sl(:,:))
          CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), &
                                     bulk_dens_sl(:,:), 0.0_wp, qmax_nh4_min_sl(:,:), qmax_nh4(:,:))
          qmax_nh4(:,:)            = qmax_nh4(:,:) * w_soil_sl(:,:) / soil_depth_sl(:,:)
          CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), &
                                     bulk_dens_sl(:,:), km_adsorpt_OM_nh4, km_adsorpt_mineral_nh4, km_adsorpt_nh4_sl(:,:))
          nh4_assoc(:,:)           = qmax_nh4(:,:) * nh4_solute (:,:) / (nh4_solute(:,:) + km_adsorpt_nh4_sl(:,:))
          nh4_n15_assoc(:,:)       = nh4_assoc(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
 
          ! metabolic litter all other elements
          sb_pool%soluable_litter%nitrogen(:,:)    = sb_pool%soluable_litter%carbon(:,:) / k_init_soluable_litter_cn 
          sb_pool%soluable_litter%phosphorus(:,:)  = sb_pool%soluable_litter%nitrogen(:,:) / k_init_litter_np 
          sb_pool%soluable_litter%carbon13(:,:)    = sb_pool%soluable_litter%carbon(:,:) * frac_c13
          sb_pool%soluable_litter%carbon14(:,:)    = sb_pool%soluable_litter%carbon(:,:) * frac_c14 
          sb_pool%soluable_litter%nitrogen15(:,:)  = sb_pool%soluable_litter%nitrogen(:,:) / & 
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! non-woody litter all other elements
          sb_pool%polymeric_litter%nitrogen(:,:)   = sb_pool%polymeric_litter%carbon(:,:) / k_init_polymeric_litter_cn 
          sb_pool%polymeric_litter%phosphorus(:,:) = sb_pool%polymeric_litter%nitrogen(:,:) / k_init_litter_np 
          sb_pool%polymeric_litter%carbon13(:,:)   = sb_pool%polymeric_litter%carbon(:,:) * frac_c13
          sb_pool%polymeric_litter%carbon14(:,:)   = sb_pool%polymeric_litter%carbon(:,:) * frac_c14
          sb_pool%polymeric_litter%nitrogen15(:,:) = sb_pool%polymeric_litter%nitrogen(:,:) / & 
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! woody litter all other elements
          sb_pool%woody_litter%nitrogen(:,:)       = sb_pool%woody_litter%carbon(:,:) / k_init_woody_litter_cn
          sb_pool%woody_litter%phosphorus(:,:)     = sb_pool%woody_litter%nitrogen(:,:) / k_init_litter_np 
          sb_pool%woody_litter%carbon13(:,:)       = sb_pool%woody_litter%carbon(:,:) * frac_c13
          sb_pool%woody_litter%carbon14(:,:)       = sb_pool%woody_litter%carbon(:,:) * frac_c14
          sb_pool%woody_litter%nitrogen15(:,:)     = sb_pool%woody_litter%nitrogen(:,:) / & 
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          ! fast pool (= microbial) all other elements
          sb_pool%microbial%nitrogen(:,:)          = sb_pool%microbial%carbon(:,:) / (k_fast_som_cn_min + k_fast_som_cn_max) * &
                                                     2.0_wp 
          sb_pool%microbial%phosphorus(:,:)        = sb_pool%microbial%nitrogen(:,:) / k_fast_som_np
          sb_pool%microbial%carbon13(:,:)          = sb_pool%microbial%carbon(:,:) * frac_c13
          sb_pool%microbial%carbon14(:,:)          = sb_pool%microbial%carbon(:,:) * frac_c14
          sb_pool%microbial%nitrogen15(:,:)        = sb_pool%microbial%nitrogen(:,:) / & 
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          sb_pool%residue%nitrogen(:,:)            = sb_pool%residue%carbon(:,:) / k_slow_som_cn
          sb_pool%residue%phosphorus(:,:)          = sb_pool%residue%nitrogen(:,:) / k_slow_som_np 
          sb_pool%residue%carbon13(:,:)            = sb_pool%residue%carbon(:,:) * frac_c13
          sb_pool%residue%carbon14(:,:)            = sb_pool%residue%carbon(:,:) * frac_c14 
          sb_pool%residue%nitrogen15(:,:)          = sb_pool%residue%nitrogen(:,:) / & 
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          nh4_n15_solute(:,:) = nh4_solute(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
          no3_n15_solute(:,:) = no3_solute(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
          
          !> init Soil P pools
          !! 
          !! SPP site-set
          IF (flag_spp1685) THEN
            ! for SPP sites, read in initial state from site-specific input data
            CALL init_soil_p_pools_spp_sites(nc, nsoil_sb, soil_lay_depth_center_sl(:,:), & 
                          dsl4jsb_Config(SB_)%soil_p_labile,   dsl4jsb_Config(SB_)%soil_p_slow, &
                          dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
                          site_ID_spp1685, p_labile_slow_site(:,:), &
                          p_labile(:,:), po4_assoc_slow(:,:), po4_occluded(:,:), po4_primary(:,:))

            IF (dsl4jsb_Config(SB_)%flag_sb_double_langmuir) THEN
              po4_assoc_fast(:,:) = calc_fast_po4(p_labile(:,:),qmax_fast_po4,km_fast_po4)
            ELSE
              po4_assoc_fast(:,:) = calc_fast_po4(p_labile(:,:),qmax_po4(:,:),km_adsorpt_po4_sl(:,:))
            ENDIF
            po4_solute(:,:)     = p_labile(:,:) - po4_assoc_fast(:,:)
            
          !! all sites but the SPP site-set
          ELSE
            ! initialise mineral soil P pools
            CALL init_soil_p_pools(nc, nsoil_sb, soil_lay_depth_center_sl(:,:), TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), & 
                                   TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
                                   soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), root_fraction_sl(:,:), &
                                   dsl4jsb_Config(SB_)%soil_p_labile, dsl4jsb_Config(SB_)%soil_p_slow, &
                                   dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
                                   dsl4jsb_Config(SB_)%soil_p_depth, p_labile_slow_site(:,:), & 
                                   po4_solute(:,:), po4_assoc_fast(:,:), po4_assoc_slow(:,:), &
                                   po4_occluded(:,:), po4_primary(:,:), qmax_po4(:,:), km_adsorpt_po4_sl(:,:), &
                                   qmax_fast_po4(:,:), qmax_slow_po4(:,:), km_fast_po4(:,:), km_slow_po4(:,:))
          ENDIF

       !!==========================================================================================================================
       !> JSM (Jena Soil Model)
       !! published with: https://doi.org/10.5194/gmd-13-783-2020
       CASE('jsm')
         ! initialisation of litter pools
         sb_pool%soluable_litter%carbon(:,1)  = def_sb_pool_metabol_litter /  molar_mass_C / soil_lay_width_sl(:,1) * &
                                                (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1)) 
         sb_pool%polymeric_litter%carbon(:,1) = def_sb_pool_jsm_struct_litter / molar_mass_C / soil_lay_width_sl(:,1) * &
                                                (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1)) 
         IF(lctlib%growthform==itree) THEN
           sb_pool%woody_litter%carbon(:,1)  = def_sb_pool_woody_litter / molar_mass_C / soil_lay_width_sl(:,1)
         ENDIF
         DO isoil=2,nsoil_sb
           sb_pool%soluable_litter%carbon(:,isoil)  = def_sb_pool_metabol_litter /  molar_mass_C / soil_lay_width_sl(:,isoil) * &
                                                      (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil) 
           sb_pool%polymeric_litter%carbon(:,isoil) = def_sb_pool_jsm_struct_litter / molar_mass_C / soil_lay_width_sl(:,isoil) * &
                                                      (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil) 
         ENDDO

         ! for each soil layer
         DO isoil=1,nsoil_sb
           sb_pool%microbial%carbon(:,isoil)  = MAX(def_sb_pool_microbial_biomass / soil_lay_width_sl(:,isoil) * &
                                                    root_fraction_sl(:,isoil), min_mic_biomass)
           sb_pool%dom%carbon(:,isoil)        = def_sb_pool_dom / soil_lay_width_sl(:,isoil) * root_fraction_sl(:,isoil)
           sb_pool%dom_assoc%carbon(:,isoil)  = k_sb_pool_dom_assoc_c_corr * sb_pool%dom%carbon(:,isoil)
           sb_pool%residue%carbon(:,isoil)    = def_sb_pool_microbial_necromass / soil_lay_width_sl(:,isoil) * &
                                                root_fraction_sl(:,isoil)
           frac_prod_not_root                 = 1._wp / isoil
           ! deep soil is initialised with very low MOC, corrected by k_deep_soil_moc_corr
           IF (frac_prod_not_root < k_deep_soil_moc_corr) THEN 
             frac_prod_not_root = frac_prod_not_root * k_deep_soil_moc_corr
           ENDIF
           sb_pool%residue_assoc%carbon(:,isoil) = dsl4jsb_Config(SB_)%qmax_org_fine_particle * &
                                                   (clay_sl(:,isoil) + silt_sl(:,isoil)) * & 
                                                   bulk_dens_sl(:,isoil) * frac_prod_not_root
           ! To avoid over-saturation of OM in soil layer
           ! calculate the fast C and stable C content (saturation at 1._wp) for each layer
           IF (isoil==1) THEN
             arr_hlp7(:,isoil) = (sb_pool%microbial%carbon(:,isoil) + sb_pool%dom%carbon(:,isoil) + &
                                 sb_pool%dom_assoc%carbon(:,isoil)) * &
                                 molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
             arr_hlp8(:,isoil) = (sb_pool%residue%carbon(:,isoil) + sb_pool%residue_assoc%carbon(:,isoil)) * &
                                 molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
           ELSE
             arr_hlp7(:,isoil) = (sb_pool%soluable_litter%carbon(:,isoil) + sb_pool%polymeric_litter%carbon(:,isoil) + & 
                                 sb_pool%microbial%carbon(:,isoil) + sb_pool%dom%carbon(:,isoil) + &
                                 sb_pool%dom_assoc%carbon(:,isoil)) * &
                                 molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
             arr_hlp8(:,isoil) = (sb_pool%residue%carbon(:,isoil) + sb_pool%residue_assoc%carbon(:,isoil)) * &
                                 molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
           ENDIF
           ! check if soil is over-saturated with fast C
           WHERE (arr_hlp7(:,isoil) > 1._wp )
             sb_pool%soluable_litter%carbon(:,isoil)  = sb_pool%soluable_litter%carbon(:,isoil)  * 1._wp / arr_hlp7(:,isoil)
             sb_pool%polymeric_litter%carbon(:,isoil) = sb_pool%polymeric_litter%carbon(:,isoil) * 1._wp / arr_hlp7(:,isoil)
             sb_pool%microbial%carbon(:,isoil)        = sb_pool%microbial%carbon(:,isoil)        * 1._wp / arr_hlp7(:,isoil)
             sb_pool%dom%carbon(:,isoil)              = sb_pool%dom%carbon(:,isoil)              * 1._wp / arr_hlp7(:,isoil)
             sb_pool%dom_assoc%carbon(:,isoil)        = sb_pool%dom_assoc%carbon(:,isoil)        * 1._wp / arr_hlp7(:,isoil)
             sb_pool%residue%carbon(:,isoil)          = 0._wp
             sb_pool%residue_assoc%carbon(:,isoil)    = 0._wp
           ! check if soil is over-saturated with fast and stable C
           ELSEWHERE ((arr_hlp7(:,isoil) + arr_hlp8(:,isoil)) > 1._wp )
             arr_hlp9(:,isoil)                     = (1._wp - arr_hlp7(:,isoil)) * f1_repartition_om_init * &
                                                     dsl4jsb_Config(SB_)%qmax_org_fine_particle * &
                                                     (clay_sl(:,isoil) + silt_sl(:,isoil)) * bulk_dens_sl(:,isoil)
             arr_hlp8(:,isoil)                     = (1._wp - arr_hlp7(:,isoil)) * (1._wp - f1_repartition_om_init) * &
                                                     rho_bulk_org * 1000._wp * carbon_per_dryweight_SOM / molar_mass_C
             sb_pool%residue_assoc%carbon(:,isoil) = calc_fast_po4(arr_hlp8(:,isoil), arr_hlp9(:,isoil), &
                                                                   k_desorpt_det / k_adsorpt_det)
             sb_pool%residue%carbon(:,isoil)       = arr_hlp8(:,isoil) - sb_pool%residue_assoc%carbon(:,isoil)
           ! check if the available mineral surface is enough for sorbed OM
           ELSEWHERE
             arr_hlp9(:,isoil) = (1._wp - arr_hlp7(:,isoil) - arr_hlp8(:,isoil)) * dsl4jsb_Config(SB_)%qmax_org_fine_particle * &
                                 (clay_sl(:,isoil) + silt_sl(:,isoil)) * bulk_dens_sl(:,isoil)
             arr_hlp8(:,isoil) = (sb_pool%residue%carbon(:,isoil) + sb_pool%residue_assoc%carbon(:,isoil))
             WHERE (arr_hlp8(:,isoil) > arr_hlp9(:,isoil))
               sb_pool%soluable_litter%carbon(:,isoil)  = sb_pool%soluable_litter%carbon(:,isoil) + &
                                                          (arr_hlp8(:,isoil)-arr_hlp9(:,isoil)) * f2_repartition_om_init
               sb_pool%polymeric_litter%carbon(:,isoil) = sb_pool%polymeric_litter%carbon(:,isoil) + &
                                                          (arr_hlp8(:,isoil)-arr_hlp9(:,isoil)) * (1._wp - f2_repartition_om_init)
               arr_hlp8(:,isoil)                        = arr_hlp9(:,isoil)
             ENDWHERE
             sb_pool%residue_assoc%carbon(:,isoil) = calc_fast_po4(arr_hlp8(:,isoil), arr_hlp9(:,isoil), &
                                                                   k_desorpt_det / k_adsorpt_det)
             sb_pool%residue%carbon(:,isoil)       = arr_hlp8(:,isoil) - sb_pool%residue_assoc%carbon(:,isoil)
           ENDWHERE
         ENDDO
     
         ! initial values of soluble nh4 and no3
         DO ichunk = 1,nc
           nh4_solute(ichunk,:)      = nh4_solute_prescribe / soil_lay_depth_ubound_sl(ichunk,nsoil_sb)
           nh4_n15_solute(ichunk,:)  = nh4_solute(ichunk,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
           no3_solute(ichunk,:)      = no3_solute_prescribe / soil_lay_depth_ubound_sl(ichunk,nsoil_sb)
           no3_n15_solute(ichunk,:)  = no3_solute(ichunk,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
         ENDDO
 
         ! inital bulk density, Qmax_org & Qmax_po4 based on current soil carbon content
         CALL calc_bulk_soil_carbon(nc, nsoil_sb, sb_pool, bulk_soil_carbon_sl, soil_litter_carbon_sl, &
                                    volume_min_sl, TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
         bulk_dens_corr_sl(:,:) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:), soil_litter_carbon_sl(:,:), &
                                                               bulk_dens_sl(:,:))
         !> po4/OC adsorption parameter intialisation 
         !! OC adsorption
         qmax_org_min_sl(:,:)   = calc_qmax_texture(dsl4jsb_Config(SB_)%qmax_org_fine_particle, silt_sl(:,:), clay_sl(:,:))
         CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), &
                                                bulk_dens_sl(:,:), 0._wp, qmax_org_min_sl(:,:), qmax_org(:,:))
         qmax_org_min_sl(:,:)   = calc_qmax_texture(dsl4jsb_Config(SB_)%qmax_org_fine_particle, silt_sl(:,:), clay_sl(:,:)) 
         CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), bulk_dens_sl(:,:), 0._wp, &
                                                qmax_org_min_sl(:,:), qmax_org(:,:))   
         !! po4 adsorption
         arr_hlp7(:,:) = po4_solute_prescribe         
         ! specific init values for SPP sites for: ph_sl(:,:) & Qmax_AlFe_cor(:,:)
         IF (flag_spp1685) THEN
           CALL init_soil_properties_spp_sites(nc, nsoil_sb, site_ID_spp1685, ph_sl, Qmax_AlFe_cor)
         ENDIF
         ! Correct the first layer Qmax_AlFe_cor with the woody litter
         Qmax_AlFe_cor(:,1) = Qmax_AlFe_cor(:,1) * 1._wp / ((1._wp - volume_min_sl(:,1)) * &
                              (bulk_soil_carbon_sl(:,1) + soil_litter_carbon_sl(:,1)  + sb_pool%woody_litter%carbon(:,1)) / &
                              (bulk_soil_carbon_sl(:,1) + soil_litter_carbon_sl(:,1)) + volume_min_sl(:,1))
         CALL calc_Psorption_parameter(dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
                                       clay_sl(:,:), silt_sl(:,:), w_soil_sl(:,:) / soil_depth_sl(:,:) * 1000._wp, ph_sl(:,:), &
                                       bulk_soil_carbon_sl(:,:), soil_litter_carbon_sl(:,:), sb_pool%woody_litter%carbon(:,:), &
                                       bulk_dens_sl(:,:), volume_min_sl(:,:), arr_hlp7(:,:), &
                                       qmax_po4_min_sl(:,:), qmax_po4_om_sl(:,:), qmax_po4(:,:), km_adsorpt_po4_sl(:,:), &
                                       qmax_fast_po4(:,:), qmax_slow_po4(:,:), km_fast_po4(:,:), km_slow_po4(:,:), &
                                       Qmax_AlFe_cor(:,:))
         ! calculate the partition coefficient between labile and slow P pool
         IF (dsl4jsb_Config(SB_)%flag_sb_double_langmuir) THEN
           WHERE (qmax_po4(:,:) > eps8)
             p_labile_slow_site(:,:)   = MAX(qmax_fast_po4(:,:) / qmax_po4(:,:),eps4)
           ELSEWHERE
             p_labile_slow_site(:,:)   = eps4
           ENDWHERE
         ELSE
           p_labile_slow_site(:,:)     = p_labile_slow_global_avg
         ENDIF

         !> nh4 adsorption parameter intialisation
         !!
         qmax_nh4_min_sl(:,:)     = calc_qmax_texture(qmax_nh4_clay, 0.0_wp, clay_sl(:,:))
         CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), bulk_dens_sl(:,:), &
                                                0.0_wp,qmax_nh4_min_sl(:,:), qmax_nh4(:,:))
         qmax_nh4(:,:)            = qmax_nh4(:,:) * w_soil_sl(:,:) / soil_depth_sl(:,:)
         CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), bulk_dens_sl(:,:), &
                                                km_adsorpt_OM_nh4, km_adsorpt_mineral_nh4, km_adsorpt_nh4_sl(:,:))
         nh4_assoc(:,:)           = qmax_nh4(:,:) * nh4_solute (:,:) / (nh4_solute(:,:) + km_adsorpt_nh4_sl(:,:))
         nh4_n15_assoc(:,:)       = nh4_assoc(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 

         !> Soil P pools intialisation
         !! 
         !! SPP site-set
         IF (flag_spp1685) THEN
           ! for SPP sites, read in initial state from site-specific input data
           CALL init_soil_p_pools_spp_sites(nc, nsoil_sb, soil_lay_depth_center_sl(:,:), & 
                            dsl4jsb_Config(SB_)%soil_p_labile,   dsl4jsb_Config(SB_)%soil_p_slow, &
                            dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
                            site_ID_spp1685, p_labile_slow_site(:,:), &
                            p_labile(:,:), po4_assoc_slow(:,:), po4_occluded(:,:), po4_primary(:,:))
           ! calculate the soluble and adsorbed P pool based on Langmuir isotherm
           IF (dsl4jsb_Config(SB_)%flag_sb_double_langmuir) THEN
             po4_assoc_fast(:,:) = calc_fast_po4(p_labile(:,:),qmax_fast_po4,km_fast_po4)
           ELSE
             po4_assoc_fast(:,:) = calc_fast_po4(p_labile(:,:),qmax_po4(:,:),km_adsorpt_po4_sl(:,:))
           ENDIF
           po4_solute(:,:)     = p_labile(:,:) - po4_assoc_fast(:,:)
         !! all sites but the SPP site-set
         ELSE
           ! mineral soil P pools intialisation
           CALL init_soil_p_pools(nc, nsoil_sb, soil_lay_depth_center_sl(:,:), TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), & 
                                 TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
                                 soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), root_fraction_sl(:,:), &
                                 dsl4jsb_Config(SB_)%soil_p_labile, dsl4jsb_Config(SB_)%soil_p_slow, &
                                 dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
                                 dsl4jsb_Config(SB_)%soil_p_depth, p_labile_slow_site(:,:), & 
                                 po4_solute(:,:), po4_assoc_fast(:,:), po4_assoc_slow(:,:), &
                                 po4_occluded(:,:), po4_primary(:,:), qmax_po4(:,:), km_adsorpt_po4_sl(:,:), &
                                 qmax_fast_po4(:,:), qmax_slow_po4(:,:), km_fast_po4(:,:), km_slow_po4(:,:))
         ENDIF ! (flag_spp1685)

         ! initialise diffusion by bioturbation
         k_bioturb(:,:) = calc_bioturbation_rate(nc, nsoil_sb,bulk_dens_corr_sl(:,:), root_fraction_sl(:,:), soil_depth_sl(:,:)) 
         
         ! initialisation of fractional enzyme allocation (range [0,1]) 
         enzyme_frac_poly(:,:)     = init_enzyme_fraction
         enzyme_frac_residue(:,:)  = init_enzyme_fraction 
         enzyme_frac_poly_c(:,:)   = init_enzyme_fraction
         enzyme_frac_poly_n(:,:)   = init_enzyme_fraction
         enzyme_frac_poly_p(:,:)   = init_enzyme_fraction

         ! isotopic signals in all pools
         ! microbial, DOM, aDOM, residue, and aRes pools
         sb_pool%soluable_litter%nitrogen(:,:)     = sb_pool%soluable_litter%carbon(:,:) / k_init_soluable_litter_cn 
         sb_pool%soluable_litter%phosphorus(:,:)   = sb_pool%soluable_litter%nitrogen(:,:) / k_init_litter_np 
         sb_pool%soluable_litter%carbon13(:,:)     = sb_pool%soluable_litter%carbon(:,:) * frac_c13 
         sb_pool%soluable_litter%carbon14(:,:)     = sb_pool%soluable_litter%carbon(:,:) * frac_c14 
         sb_pool%soluable_litter%nitrogen15(:,:)   = sb_pool%soluable_litter%nitrogen(:,:) / &
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
         sb_pool%polymeric_litter%nitrogen(:,:)    = sb_pool%polymeric_litter%carbon(:,:) / k_init_polymeric_litter_cn
         sb_pool%polymeric_litter%phosphorus(:,:)  = sb_pool%polymeric_litter%nitrogen(:,:) / k_init_litter_np 
         sb_pool%polymeric_litter%carbon13(:,:)    = sb_pool%polymeric_litter%carbon(:,:) * frac_c13 
         sb_pool%polymeric_litter%carbon14(:,:)    = sb_pool%polymeric_litter%carbon(:,:) * frac_c14 
         sb_pool%polymeric_litter%nitrogen15(:,:)  = sb_pool%polymeric_litter%nitrogen(:,:) / &
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
         sb_pool%woody_litter%nitrogen(:,:)        = sb_pool%woody_litter%carbon(:,:) / k_init_woody_litter_cn
         sb_pool%woody_litter%phosphorus(:,:)      = sb_pool%woody_litter%nitrogen(:,:) / k_init_litter_np
         sb_pool%woody_litter%carbon13(:,:)        = sb_pool%woody_litter%carbon(:,:) * frac_c13 
         sb_pool%woody_litter%carbon14(:,:)        = sb_pool%woody_litter%carbon(:,:) * frac_c14 
         sb_pool%woody_litter%nitrogen15(:,:)      = sb_pool%woody_litter%nitrogen(:,:) / &
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
         sb_pool%microbial%nitrogen(:,:)           = sb_pool%microbial%carbon(:,:) / microbial_cn 
         sb_pool%microbial%phosphorus(:,:)         = sb_pool%microbial%nitrogen(:,:) / microbial_np
         sb_pool%microbial%carbon13(:,:)           = sb_pool%microbial%carbon(:,:) * frac_c13 
         sb_pool%microbial%carbon14(:,:)           = sb_pool%microbial%carbon(:,:) * frac_c14 
         sb_pool%microbial%nitrogen15(:,:)         = sb_pool%microbial%nitrogen(:,:) / &
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
         sb_pool%residue%nitrogen(:,:)             = sb_pool%residue%carbon(:,:) / microbial_cn 
         sb_pool%residue%phosphorus(:,:)           = sb_pool%residue%nitrogen(:,:) / k_slow_som_np
         sb_pool%residue%carbon13(:,:)             = sb_pool%residue%carbon(:,:) * frac_c13 
         sb_pool%residue%carbon14(:,:)             = sb_pool%residue%carbon(:,:) * frac_c14 
         sb_pool%residue%nitrogen15(:,:)           = sb_pool%residue%nitrogen(:,:) / &
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
         sb_pool%residue_assoc%nitrogen(:,:)       = sb_pool%residue_assoc%carbon(:,:) / microbial_cn
         sb_pool%residue_assoc%phosphorus(:,:)     = sb_pool%residue_assoc%nitrogen(:,:) / k_slow_som_np
         sb_pool%residue_assoc%carbon13(:,:)       = sb_pool%residue_assoc%carbon(:,:) * frac_c13 
         sb_pool%residue_assoc%carbon14(:,:)       = sb_pool%residue_assoc%carbon(:,:) * frac_c14 
         sb_pool%residue_assoc%nitrogen15(:,:)     = sb_pool%residue_assoc%nitrogen(:,:) / &
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
         sb_pool%dom%nitrogen(:,:)                 = sb_pool%dom%carbon(:,:) / microbial_cn
         sb_pool%dom%phosphorus(:,:)               = sb_pool%dom%nitrogen(:,:) / microbial_np
         sb_pool%dom%carbon13(:,:)                 = sb_pool%dom%carbon(:,:) * frac_c13
         sb_pool%dom%carbon14(:,:)                 = sb_pool%dom%carbon(:,:) * frac_c14
         sb_pool%dom%nitrogen15(:,:)               = sb_pool%dom%nitrogen(:,:) / &
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
         sb_pool%dom_assoc%nitrogen(:,:)           = sb_pool%dom_assoc%carbon(:,:) / microbial_cn
         sb_pool%dom_assoc%phosphorus(:,:)         = sb_pool%dom_assoc%nitrogen(:,:) / microbial_np
         sb_pool%dom_assoc%carbon13(:,:)           = sb_pool%dom_assoc%carbon(:,:) * frac_c13
         sb_pool%dom_assoc%carbon14(:,:)           = sb_pool%dom_assoc%carbon(:,:) * frac_c14
         sb_pool%dom_assoc%nitrogen15(:,:)         = sb_pool%dom_assoc%nitrogen(:,:) / & 
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
         ! fungi and mycorrhiza      
         sb_pool%fungi%carbon13(:,:)               = sb_pool%fungi%carbon(:,:) * frac_c13
         sb_pool%fungi%carbon14(:,:)               = sb_pool%fungi%carbon(:,:) * frac_c14 
         sb_pool%fungi%nitrogen15(:,:)             = sb_pool%fungi%nitrogen(:,:) / & 
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
         sb_pool%mycorrhiza%carbon13(:,:)          = sb_pool%mycorrhiza%carbon(:,:) * frac_c13
         sb_pool%mycorrhiza%carbon14(:,:)          = sb_pool%mycorrhiza%carbon(:,:) * frac_c14 
         sb_pool%mycorrhiza%nitrogen15(:,:)        = sb_pool%mycorrhiza%nitrogen(:,:) / & 
                                                    (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

       CASE DEFAULT
         CALL finish(routine,'Invalid setting of sb_model_scheme.')
     END SELECT  ! simple_1d | jsm

  ! Deallocate local allocatable variables
  DEALLOCATE(p_labile)
  DEALLOCATE(arr_hlp6)
  DEALLOCATE(arr_hlp7)
  DEALLOCATE(arr_hlp8)
  DEALLOCATE(arr_hlp9)
  DEALLOCATE(root_depth_avg)
  DEALLOCATE(p_labile_slow_site)


! Pass back values to process memory
dsl4slm_passback_vegpool_var(pool,fine_root,carbon)
#define Passback_sb_pool
#include "slm_pools_a_incl.f90"
#undef Passback_sb_pool

  END SUBROUTINE sb_init

  !-----------------------------------------------------------------------------------------------------
  !> Intialize soil P pools for all layers from the soil P maps (Yang et al. 2013)
  !! assuming that the soil P pool sizes are correlated with soil weight
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_soil_p_pools(nc, nsoil_sb, soil_lay_depth_center_sl, &
                               sb_model_scheme, sb_adsorp_scheme, flag_sb_double_langmuir, & 
                               soil_depth_sl, bulk_dens_sl, root_fraction_sl,&
                               soil_p_labile, soil_p_slow, soil_p_occluded, soil_p_primary, &
                               soil_p_depth, p_labile_slow_ratio, &
                               po4_solute, po4_assoc_fast, po4_assoc_slow, po4_occluded, po4_primary, &
                               qmax_po4,km_po4, qmax_fast_po4, km_fast_po4, qmax_slow_po4, km_slow_po4)

    USE mo_jsb_phys_constants,    ONLY: molar_mass_P
    USE mo_sb_constants           ! e.g.: k1_p_pool_depth_corr, k2_p_pool_depth_corr, k3_p_pool_depth_corr, &
                                  !       k4_p_pool_depth_corr, k_p_pool_depth_corr_min_occlud, p_labile_slow_global_avg
    USE mo_sb_jsm_processes,      ONLY: calc_fast_po4

    IMPLICIT NONE

    INTEGER,                          INTENT(in)     :: nc                        !< dimensions
    INTEGER                         , INTENT(in)     :: nsoil_sb                  !< number of soil layers as used/defined by the SB_ process
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)     :: soil_lay_depth_center_sl  !< depth at the center of each soil layer
    CHARACTER(len=*)                , INTENT(in)     :: sb_model_scheme       !< soil model scheme 
    CHARACTER(len=*)                , INTENT(in)     :: sb_adsorp_scheme      !< soil model po4 adsorption scheme
    LOGICAL                         , INTENT(in)     :: flag_sb_double_langmuir   !< T: double Langmuir; F: traditional Langmuir    
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)     :: soil_depth_sl, &      !< soil depth per soil layer, [m]
                                                        root_fraction_sl, &   !< 
                                                        bulk_dens_sl, &       !< corrected soil bulk density per soil layer, [kg m-3]
                                                        p_labile_slow_ratio   !< ratio between labile and labile+slow P pools
    REAL(wp)                        , INTENT(in)     :: soil_p_depth, &       !< soil depth of the soil P maps
                                                        soil_p_labile, &      !< total labile P pool (solute + assoc_fast) of the soil P maps, [g m-2]
                                                        soil_p_slow, &        !< total slow P pool (secondary P) of the soil P maps, [g m-2]
                                                        soil_p_occluded, &    !< total occluded P pool of the soil P maps, [g m-2]
                                                        soil_p_primary        !< total primary P pool of the soil P maps, [g m-2]
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(inout)  :: po4_solute, &         !< soluble Pi pool per soil layer, [mol m-3]
                                                        po4_assoc_fast, &     !< fast minerally associated PO4 pool per soil layer, [mol m-3]
                                                        po4_assoc_slow, &     !< slow minerally associated PO4 pool per soil layer, [mol m-3]
                                                        po4_occluded, &       !< occluded PO4 pool per soil layer, [mol m-3]
                                                        po4_primary, &        !< primary PO4 pool per soil layer, [mol m-3]
                                                        qmax_po4, &           !< maximum PO4 sorption capacity per soil layer, [mol m-3]
                                                        km_po4, &             !< half-saturation PO4 concentration for sorption per soil layer, [mol m-3]
                                                        qmax_fast_po4, &      !< maximum PO4 sorption capacity for po4_assoc_fast, [mol m-3]
                                                        km_fast_po4, &        !< half-saturation PO4 concentration for po4_assoc_fast, [mol m-3]
                                                        qmax_slow_po4, &      !< maximum PO4 sorption capacity for po4_assoc_slow, [mol m-3]
                                                        km_slow_po4           !< half-saturation PO4 concentration for po4_assoc_slow, [mol m-3]
    INTEGER                                          :: i
    REAL(wp)                                         :: rhlp1
    REAL(wp), DIMENSION(nc)                          :: soil_p_weight, soil_weight_0d, &  !< total soil weight up to the soil_p_depth/soil_depth_0d, [kg]
                                                        soil_p_layer                      !< up to soil_p_layer keep the P init map info
    REAL(wp), DIMENSION(nc)                          :: hlp1, hlp2, hlp3, hlp4
    REAL(wp), DIMENSION(nc)                          :: soil_depth_0d, root_depth_avg
    LOGICAL,  DIMENSION(nc)                          :: aflag  
    REAL(wp), DIMENSION(nc,nsoil_sb)                 :: soil_weight_sl, &                 !< soil weight per soil layer
                                                        arr_hlp1
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_soil_p_pools'


    !> 0.9 init output variables
    !! 
    po4_solute(:,:)     = 0._wp
    po4_assoc_fast(:,:) = 0._wp
    po4_assoc_slow(:,:) = 0._wp 
    po4_occluded(:,:)   = 0._wp
    po4_primary(:,:)    = 0._wp

    !> 1.0 calculate the total soil weight up to the soil_p_depth,
    !! and the soil weight per soil layer
    hlp1(:)  = 0._wp
    hlp2(:)  = 0._wp
    aflag(:) = .TRUE.
    rhlp1    = soil_p_slow + soil_p_labile + soil_p_occluded + soil_p_primary
    DO i = 1,nsoil_sb
      soil_weight_sl(:,i) = soil_depth_sl(:,i) * bulk_dens_sl(:,i)
      WHERE (aflag) 
        hlp1(:) = hlp1(:) + soil_depth_sl(:,i)
        hlp2(:) = hlp2(:) + soil_weight_sl(:,i)
        WHERE (hlp1(:) > soil_p_depth)
          soil_p_weight(:) = hlp2(:) - soil_weight_sl(:,i) * (hlp1(:) - soil_p_depth) / soil_depth_sl(:,i)
          aflag(:)         = .FALSE.
          soil_p_layer(:)  = i
        ENDWHERE
      ENDWHERE
    END DO
    hlp4(:) = rhlp1 / soil_p_weight(:) ! the P density from init map, [g P / kg soil]
    
    !> 2.0 calculate the initial slow P, occluded P, primary P pool, and p_labile
    !! based on the assumption that soil P pool are linearly corrected with soil weight
    !! the depth correction of initial P pool distributions are fitted against the observed pattern of a young cambisol soil 
    !! with basalt as parent material (BBR, Lang et al. 2017) that, 
    !! the proportion of: labile+sorbed P decreases with increasing depth, and primary P increases with increasing depth (Yu et al., 2020)
    DO i=1,nsoil_sb
      ! primary P, corrected with depth, [g P] 
      hlp1(:) = MIN(k2_p_pool_depth_corr * soil_lay_depth_center_sl(:,i) ** k4_p_pool_depth_corr, &
                    (k1_p_pool_depth_corr - k_p_pool_depth_corr_min_occlud)) / k1_p_pool_depth_corr
      hlp2(:) = k3_p_pool_depth_corr ** (k2_p_pool_depth_corr * soil_lay_depth_center_sl(:,i)) / k1_p_pool_depth_corr
      ! occluded P, corrected with depth
      hlp3(:) = MAX((1.0_wp - hlp1(:) - hlp2(:)), k_p_pool_depth_corr_min_occlud / k1_p_pool_depth_corr) 

      WHERE ((hlp1(:) + hlp2(:) + hlp3(:)) > 1._wp)
        hlp1(:) = hlp1(:) / (hlp1(:) + hlp2(:) + hlp3(:))
        hlp2(:) = hlp2(:) / (hlp1(:) + hlp2(:) + hlp3(:))
        hlp3(:) = hlp3(:) / (hlp1(:) + hlp2(:) + hlp3(:))
      ENDWHERE
      WHERE (soil_p_layer(:) > i) 
        po4_assoc_slow(:,i) = bulk_dens_sl(:,i) * hlp4(:) * (soil_p_slow + soil_p_labile) / rhlp1 / molar_mass_P * &
                              (1._wp - p_labile_slow_ratio(:,i))
        po4_occluded(:,i)   = bulk_dens_sl(:,i) * hlp4(:) * soil_p_occluded               / rhlp1 / molar_mass_P
        po4_primary(:,i)    = bulk_dens_sl(:,i) * hlp4(:) * soil_p_primary                / rhlp1 / molar_mass_P
        arr_hlp1(:,i)       = bulk_dens_sl(:,i) * hlp4(:) * (soil_p_slow + soil_p_labile) / rhlp1 / molar_mass_P * &
                              p_labile_slow_ratio(:,i) 
      ELSEWHERE
        po4_assoc_slow(:,i) = bulk_dens_sl(:,i) * hlp4(:) * hlp2(:) / molar_mass_P * (1._wp - p_labile_slow_ratio(:,i))
        po4_occluded(:,i)   = bulk_dens_sl(:,i) * hlp4(:) * hlp3(:) / molar_mass_P
        po4_primary(:,i)    = bulk_dens_sl(:,i) * hlp4(:) * hlp1(:) / molar_mass_P
        arr_hlp1(:,i)       = bulk_dens_sl(:,i) * hlp4(:) * hlp2(:) / molar_mass_P * p_labile_slow_ratio(:,i) 
      ENDWHERE
    END DO

    !> 2.1 calculate the initial solute P and fast P,
    !! based on the assumption that labile P is the sum of both,
    !! and they follow the Langmuir equilibrium
    !! po4_assoc_fast = (p_labile_sl + Smax + K - sqrt((p_labile_sl + Smax + K)**2._wp - 4*p_labile_sl*Smax)) / 2._wp
    IF (flag_sb_double_langmuir) THEN
      DO i=1,nsoil_sb
        po4_assoc_fast(:,i) = calc_fast_po4(arr_hlp1(:,i), qmax_fast_po4(:,i), km_fast_po4(:,i))
        po4_solute(:,i)     = arr_hlp1(:,i) - po4_assoc_fast(:,i)
      END DO
    ELSE
      DO i=1,nsoil_sb
        po4_assoc_fast(:,i) = calc_fast_po4(arr_hlp1(:,i), qmax_po4(:,i), km_po4(:,i))
        po4_solute(:,i)     = arr_hlp1(:,i) - po4_assoc_fast(:,i)
      END DO
    ENDIF

  END SUBROUTINE init_soil_p_pools

  !-----------------------------------------------------------------------------------------------------
  !> init soil properties for SPP sites only
  !! 
  !! Output: ph and Qmax_AlFe_cor prescribed for spp sites
  !!
  !! Qmax_AlFe_cor is further modified with calc_Psorption_parameter()
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_soil_properties_spp_sites(nc, nsoil_sb, site_ID_spp1685, ph_sl, Qmax_AlFe_cor)

    IMPLICIT NONE

    INTEGER,                                          INTENT(in)    :: nc, &                        !< dimensions
                                                                       nsoil_sb
    CHARACTER(len=*), OPTIONAL,                       INTENT(in)    :: site_ID_spp1685              !< site ID for spp1685  
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(out)   :: ph_sl                        !< ph of soil layer [-]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(out)   :: Qmax_AlFe_cor                !< qmax correlation factor of Al/Fe content [-]  

    INTEGER  :: i

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_soil_properties_spp_sites'
    
    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 1.0 prescribe ph and Qmax_AlFe_cor from measurements
    !! 
    SELECT CASE(site_ID_spp1685)
    CASE ("BBR")
      DO i=1,nc
        ph_sl(i,1:nsoil_sb)         = (/3.84_wp, 4.14_wp, 4.45_wp, 4.66_wp, 4.85_wp, &
                                        4.92_wp,5.19_wp,5.64_wp,5.8_wp,6._wp, &
                                        6.2_wp, 6.5_wp, 6.5_wp, 6.5_wp, 6.5_wp/)
        Qmax_AlFe_cor(i,1:nsoil_sb) = (/3.04_wp, 3.60_wp, 3.60_wp, 3.295_wp, 2.92_wp, &
                                        2.75_wp,2.40_wp,2.04_wp,2._wp,2._wp, &
                                        2._wp, 2._wp, 2._wp, 2._wp, 2._wp/)
      ENDDO
    CASE ("LUE")
      DO i=1,nc
        ph_sl(i,1:nsoil_sb)         = (/3.51_wp,3.78_wp,3.84_wp,4.14_wp,4.20_wp, &
                                        4.46_wp,4.45_wp,4.62_wp,4.63_wp,4.8_wp, &
                                        5._wp, 5.2_wp, 5.5_wp, 6._wp, 6.5_wp/)
        Qmax_AlFe_cor(i,1:nsoil_sb) = (/0.17_wp, 0.29_wp, 0.41_wp, 0.41_wp, 0.396_wp, &
                                        0.3_wp,0.289_wp,0.236_wp,0.2_wp,0.2_wp, &
                                        0.2_wp, 0.2_wp, 0.2_wp, 0.2_wp, 0.2_wp/)
      ENDDO
    CASE ("VES")
      DO i=1,nc
        ph_sl(i,1:nsoil_sb)         = (/3.47_wp,3.91_wp,4.25_wp,4.46_wp,4.52_wp, &
                                        4.66_wp,4.78_wp,4.77_wp,5._wp,5.2_wp, &
                                        5.5_wp, 5.8_wp, 6.0_wp, 6.5_wp, 6.5_wp/)
        Qmax_AlFe_cor(i,1:nsoil_sb) = (/1.275_wp, 1.515_wp, 1.378_wp, 1.459_wp, 1.90_wp, &
                                        2.469_wp,3.428_wp,5.33_wp,6._wp,6._wp, &
                                        6._wp, 6._wp, 6._wp, 6._wp, 6._wp/)
      ENDDO
    CASE ("MIT")
      DO i=1,nc
        ph_sl(i,1:nsoil_sb)         = (/3.64_wp,4.02_wp,4.53_wp,4.46_wp,4.46_wp, &
                                        4.58_wp,4.56_wp,4.64_wp,4.7_wp,4.8_wp, &
                                        5._wp, 5.2_wp, 5.5_wp, 5.8_wp, 6._wp/)
        Qmax_AlFe_cor(i,1:nsoil_sb) = (/1.326_wp, 1.472_wp, 1.573_wp, 1.5_wp, 1.426_wp, &
                                        1.3_wp,1.258_wp,1.2_wp,1.2_wp,1.2_wp, &
                                        1.1_wp, 1.1_wp, 1._wp, 1._wp, 1._wp/)
      ENDDO
    CASE ("CON")
      DO i=1,nc
        ph_sl(i,1:nsoil_sb)         = (/4.03_wp,4.12_wp,4.32_wp,4.37_wp,4.42_wp, &
                                        4.54_wp,4.60_wp,4.61_wp,4.7_wp,4.8_wp, &
                                        5._wp, 5.2_wp, 5.5_wp, 5.8_wp, 6._wp/)
        Qmax_AlFe_cor(i,1:nsoil_sb) = (/1.9_wp, 1.891_wp, 1.87_wp, 1.8_wp, 1.73_wp, &
                                        1.588_wp,1.436_wp,1.24_wp,1.2_wp,1.2_wp, &
                                        1.1_wp, 1.1_wp, 1._wp, 1._wp, 1._wp/)
      ENDDO
    END SELECT
                                      
  END SUBROUTINE init_soil_properties_spp_sites

  !-----------------------------------------------------------------------------------------------------
  ! SPP site-set
  !
  !-----------------------------------------------------------------------------------------------------
  !> Intialize soil_biogeochemistry process specifically for all SPP sites (after memory has been set up)
  !! site of the SPP site-set: BBR, MIT, VES, CON, LUE
  !! 
  !! @TODO this routine works only at one point! because soil properties for exactly one site are selected
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_soil_p_pools_spp_sites (nc, nsoil_sb, soil_lay_depth_center_sl, &
                                        soil_p_labile, soil_p_slow, soil_p_occluded, soil_p_primary, &
                                        sitename, p_labile_slow_ratio, &
                                        p_labile, po4_assoc_slow, po4_occluded, po4_primary)

    USE mo_sb_constants !,         ONLY: 

    IMPLICIT NONE

    INTEGER                         , INTENT(in)     :: nc                        !< dimensions
    INTEGER                         , INTENT(in)     :: nsoil_sb                  !< number of soil layers as used/defined by the SB_ process
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)     :: soil_lay_depth_center_sl  !< depth at the center of each soil layer
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)     :: p_labile_slow_ratio   !< ratio between labile and labile+slow P pools
    REAL(wp)                        , INTENT(in)     :: soil_p_labile, &      !< labile P from Yang et al. 2013, g/m2
                                                        soil_p_slow, &        !< secondary P from Yang et al. 2013, g/m2
                                                        soil_p_occluded, &    !< Occluded P from Yang et al. 2013, g/m2
                                                        soil_p_primary        !< Primary from Yang et al. 2013, g/m2
    CHARACTER(len=*)                , INTENT(in)     :: sitename              !< spp site name: BBR or LUE (MIT, VES, CON)
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(inout)  :: p_labile, &           !< PO4 in solution per soil layer [mol / m3]
                                                        po4_assoc_slow, &     !< slow minerally associated PO4 pool per soil layer [mol/m3]
                                                        po4_occluded, &       !< occluded PO4 pool per soil layer [mol/m3]
                                                        po4_primary           !< primary PO4 pool per soil layer [mol/m3]

    INTEGER                                      :: isoil, ichunk
    REAL(wp)                                     :: hlp1, hlp2, hlp3, hlp4, hlp5, hlp6
    REAL(wp), DIMENSION(nc)                      :: arr_hlp2, arr_hlp3, arr_hlp4, soil_p_layer
    REAL(wp), DIMENSION(nsoil_sb)                :: arr_hlp5
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_soil_p_pools_spp_sites'

    !> 1.0 set site-level soil properties
    !!
    SELECT CASE (sitename)
      CASE("BBR")
        arr_hlp5(:) = (/9.24_wp,36.92_wp,39.44_wp,40.36_wp,42.75_wp, &  
                       45.72_wp,32.45_wp,95.89_wp,124.14_wp,133.93_wp, &
                       100._wp,126._wp,122.18_wp,145._wp,90._wp/)       ! observed total SIP, mol P m-3
        hlp1            = soil_p_labile     ! p_labile from Yang et al. 2013, g/m2
        hlp2            = soil_p_slow       ! p_slow from Yang et al. 2013, g/m2
        hlp3            = soil_p_primary    ! p_occluded from Yang et al. 2013, g/m2
        hlp4            = soil_p_occluded   ! p_primary from Yang et al. 2013, g/m2
        soil_p_layer(:) = 6._wp             ! field measurement depth corresponded model soil layers
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        !hlp5        = 1750._wp / 1000._wp ! g P / kg soil
        !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
      CASE("LUE")
        arr_hlp5(:) = (/0.8867_wp,3.187_wp,2.662_wp,4.887_wp,4.792_wp, &  
                       2.987_wp,3.313_wp,6.0616_wp,7.56_wp,5.787_wp, &
                       5.932_wp,6.845_wp,8.294_wp,8.413_wp,8.532_wp/)   ! observed total SIP, mol P m-3
        hlp1            = soil_p_labile     ! p_labile from Yang et al. 2013, g/m2
        hlp2            = soil_p_slow       ! p_slow from Yang et al. 2013, g/m2
        hlp3            = soil_p_occluded   ! p_occluded from Yang et al. 2013, g/m2
        hlp4            = soil_p_primary    ! p_primary from Yang et al. 2013, g/m2
        soil_p_layer(:) = 6._wp             ! field measurement depth corresponded model soil layers
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        !hlp5        = 120._wp / 1000._wp ! g P / kg soil
        !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
      CASE("VES")
        arr_hlp5(:) = (/2.388_wp,8.509_wp,9.033_wp,10.811_wp,15.365_wp, &  
                       20.661_wp,22.721_wp,29.321_wp,30.952_wp,32.821_wp, &
                       40.673_wp,41.880_wp,41.111_wp,48.438_wp,54._wp/)   ! observed total SIP, mol P m-3
        hlp1            = soil_p_labile     ! p_labile from Yang et al. 2013, g/m2
        hlp2            = soil_p_slow       ! p_slow from Yang et al. 2013, g/m2
        hlp3            = soil_p_primary    ! p_occluded from Yang et al. 2013, g/m2
        hlp4            = soil_p_occluded   ! p_primary from Yang et al. 2013, g/m2
        soil_p_layer(:) = 6._wp             ! field measurement depth corresponded model soil layers
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        !hlp5        = 900._wp / 1000._wp ! g P / kg soil
        !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
      CASE("MIT")
        arr_hlp5(:) = (/3.723_wp,6.272_wp,7.556_wp,8.292_wp,6.93_wp, &  
                       13.731_wp,13.296_wp,16.511_wp,12.890_wp,15._wp, &
                       16.667_wp,16.875_wp,17.5_wp,15.104_wp,15.625_wp/)  ! observed total SIP, mol P m-3
        hlp1            = soil_p_labile     ! p_labile from Yang et al. 2013, g/m2
        hlp2            = soil_p_slow       ! p_slow from Yang et al. 2013, g/m2
        hlp3            = soil_p_primary    ! p_occluded from Yang et al. 2013, g/m2
        hlp4            = soil_p_occluded   ! p_primary from Yang et al. 2013, g/m2
        soil_p_layer(:) = 9._wp             ! field measurement depth corresponded model soil layers
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        !hlp5        = 900._wp / 1000._wp ! g P / kg soil
        !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
      CASE("CON")
        arr_hlp5(:) = (/1.862_wp,2.924_wp,4.027_wp,4.048_wp,4.069_wp, &  
                       7.167_wp,10.855_wp,7.719_wp,18.75_wp,18.33_wp, &
                       23.44_wp,25._wp,21.093_wp,21.875_wp,15.104_wp/)   ! observed total SIP, mol P m-3
        hlp1            = soil_p_labile     ! p_labile from Yang et al. 2013, g/m2
        hlp2            = soil_p_slow       ! p_slow from Yang et al. 2013, g/m2
        hlp3            = soil_p_primary    ! p_occluded from Yang et al. 2013, g/m2
        hlp4            = soil_p_occluded   ! p_primary from Yang et al. 2013, g/m2
        soil_p_layer(:) = 9._wp             ! field measurement depth corresponded model soil layers
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        !hlp5        = 900._wp / 1000._wp ! g P / kg soil
        !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
    END SELECT

    !> 2.0 calculate the initial slow P, occluded P, primary P pool, and p_labile
    !! @TODO this is inconsistent with section 1.0 of this routine, as section 1.0 does not work for multiple points in a chunk
    !! based on the assumption that soil P pool are linearly corrected with soil weight
    !! the depth correction of initial P pool distributions are fitted against the observed pattern of a young cambisol soil 
    !! with basalt as parent material (BBR, Lang et al. 2017) that, 
    !! the proportion of: labile+sorbed P decreases with increasing depth, and primary P increases with increasing depth (Yu et al., 2020)
    DO isoil = 1, nsoil_sb
      ! primary corrected for depth, [g P] 
      arr_hlp4(:) = MIN(k2_p_pool_depth_corr * soil_lay_depth_center_sl(:,isoil) ** k4_p_pool_depth_corr, &
                    (k1_p_pool_depth_corr - k_p_pool_depth_corr_min_occlud)) / k1_p_pool_depth_corr
      ! labile + slow
      arr_hlp2(:) = k3_p_pool_depth_corr ** (k2_p_pool_depth_corr * soil_lay_depth_center_sl(:,isoil)) / k1_p_pool_depth_corr
      ! occluded corrected for depth
      arr_hlp3(:) = MAX((1.0_wp - arr_hlp4(:) - arr_hlp2(:)), k_p_pool_depth_corr_min_occlud / k1_p_pool_depth_corr) 

      WHERE ((arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:)) > 1._wp)
        arr_hlp4(:) = arr_hlp4(:) / (arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:))
        arr_hlp2(:) = arr_hlp2(:) / (arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:))
        arr_hlp3(:) = arr_hlp3(:) / (arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:))
      ENDWHERE
      
      WHERE (isoil < soil_p_layer(:)) 
        po4_assoc_slow(:,isoil) = arr_hlp5(isoil) * (hlp1 + hlp2) / (hlp1 + hlp2 + hlp3 + hlp4) * &
                                  (1._wp - p_labile_slow_ratio(:,isoil))
        po4_occluded(:,isoil)   = arr_hlp5(isoil) * hlp3 / (hlp1 + hlp2 + hlp3 + hlp4)
        po4_primary(:,isoil)    = arr_hlp5(isoil) * hlp4 / (hlp1 + hlp2 + hlp3 + hlp4)
        p_labile(:,isoil)       = arr_hlp5(isoil) * (hlp1 + hlp2) / (hlp1 + hlp2 + hlp3 + hlp4) * p_labile_slow_ratio(:,isoil)
      ELSEWHERE
        po4_assoc_slow(:,isoil) = arr_hlp5(isoil) * arr_hlp2(:) * (1._wp - p_labile_slow_ratio(:,isoil))
        po4_occluded(:,isoil)   = arr_hlp5(isoil) * arr_hlp3(:)
        po4_primary(:,isoil)    = arr_hlp5(isoil) * arr_hlp4(:)
        p_labile(:,isoil)       = arr_hlp5(isoil) * arr_hlp2(:) * p_labile_slow_ratio(:,isoil)
      ENDWHERE
    ENDDO

  END SUBROUTINE init_soil_p_pools_spp_sites

END MODULE mo_sb_init
