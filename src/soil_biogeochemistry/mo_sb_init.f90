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
    USE mo_jsb_math_constants,              ONLY: eps8
    USE mo_assimi_constants,                ONLY: CiCa_default_C3, CiCa_default_C4, ic3phot, ic4phot
    USE mo_assimi_update_canopy_fluxes,     ONLY: discrimination_ps
    USE mo_sb_jsm_processes                 ! e.g. calc_bulk_soil_carbon(), convert_s1d_to_s0d_var_soilproperty()
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
    REAL(wp)                                     :: hlp1, hlp2, hlp3, hlp4, hlp5, &
                                                    hlp7, hlp8, hlp9, hlp10, hlp11, &
                                                    polyLit_cn_ratio, woodyLit_cn_ratio
    REAL(wp)                                     :: frac_prod_not_root
    REAL(wp)                                     :: frac_c13,frac_c14
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)        :: fine_root_carbon_sl, labile_p, soil_depth_acc_sl, arr_hlp7, arr_hlp8, arr_hlp9
    REAL(wp), ALLOCATABLE, DIMENSION(:)          :: arr_hlp6, root_depth_avg  ! dim: nc
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
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
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
    dsl4jsb_Real3D_onChunk :: km_adsorpt_po4_sl
    dsl4jsb_Real3D_onChunk :: km_adsorpt_nh4_sl
    dsl4jsb_Real3D_onChunk :: k_bioturb
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast
    dsl4jsb_Real3D_onChunk :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: po4_occluded
    dsl4jsb_Real3D_onChunk :: ph_sl
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
    ALLOCATE(fine_root_carbon_sl(nc, nsoil_sb))
    ALLOCATE(labile_p(nc, nsoil_sb))
    ALLOCATE(soil_depth_acc_sl(nc, nsoil_sb))
    ALLOCATE(arr_hlp6(nc))
    ALLOCATE(arr_hlp7(nc, nsoil_sb))
    ALLOCATE(arr_hlp8(nc, nsoil_sb))
    ALLOCATE(arr_hlp9(nc, nsoil_sb))
    ALLOCATE(root_depth_avg(nc))
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var3D_onChunk(VEG_,      root_fraction_sl)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SSE_,     bulk_dens_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     w_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     qmax_org_min_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     qmax_po4_min_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     qmax_nh4_min_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     qmax_po4_om_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     sand_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     silt_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     clay_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     volume_min_sl)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       noy)
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,       n2)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_soil_carbon_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       soil_litter_carbon_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_org)
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,       km_adsorpt_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       km_adsorpt_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       k_bioturb)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_fast)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_slow)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_occluded) 
    dsl4jsb_Get_var3D_onChunk(SB_,       ph_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       lithology_class_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_primary)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_cue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_nue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_pue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,       vmax_weath_mineral_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly)
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_residue)
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly_c)
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly_n)
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly_p)
    dsl4jsb_Get_var3D_onChunk(SB_,       dom_cn)
    dsl4jsb_Get_var3D_onChunk(SB_,       dom_cp)
    ! Get process memory pool structures
    dsl4slm_get_vegpool_var(pool,fine_root,carbon)      ! needs veg_pool
#define Get_sb_pool
#include "slm_pools_a_incl.f90"
#undef Get_sb_pool



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !>   C3 / C4 plants
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


    !> 1.0 all soil models
    !! 
    ph_sl(:,:)              = dsl4jsb_Config(SB_)%soil_ph
    lithology_class_sl(:,:) = dsl4jsb_Config(SB_)%lithology_class
    ! SPP site-set
    IF (flag_spp1685) THEN
        DO isoil = 1,nsoil_sb
          IF (isoil == 1) THEN
            soil_depth_acc_sl(:,isoil) = soil_depth_sl(:,isoil) / 2._wp
          ELSE
            soil_depth_acc_sl(:,isoil) = soil_depth_acc_sl(:,isoil-1) + (soil_depth_sl(:,isoil) + soil_depth_sl(:,isoil-1)) / 2._wp  
          ENDIF  
        ENDDO 

        ! set site ID of SPP site-set from longitude
        IF (ABS(site_longitude - 9.75_wp) < eps8) THEN
          site_ID_spp1685 = "BBR"
        ELSEIF (ABS(site_longitude - 10.25_wp) < eps8) THEN
          site_ID_spp1685 = "LUE"
        ELSEIF (ABS(site_longitude - 10.75_wp) < eps8) THEN
          site_ID_spp1685 = "VES"
        ELSEIF (ABS(site_longitude - 12.75_wp) < eps8) THEN
          site_ID_spp1685 = "MIT"
        ELSEIF (ABS(site_longitude - 7.75_wp) < eps8) THEN
          site_ID_spp1685 = "CON"
        ELSE
          site_ID_spp1685 = "ERR"
        ENDIF
    ENDIF

    !> 2.0 soil model init
    !! differ between simple_0d/1d and JSM
    !! 
    frac_prod_not_root = 0.6_wp   ! fraction of production not roots 
    SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
       !!==========================================================================================================================
       !> simple_0d & simple_1d soil models
       !! Initial values for Soenke's four pool soil biogeochemical model
       CASE("simple_0d", "simple_1d")

          hlp1 = 100._wp ! g C m-2 metabolic litter
          hlp2 = 100._wp ! g C m-2 structural litter
          hlp3 = 1000._wp ! g C m-2 woody litter
          hlp4 = 1500._wp ! g C m-2 fast SOM
          hlp5 = 10000._wp ! g C m-2 slow SOM

          IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_0d") THEN
             sb_pool%soluable_litter%carbon(:,1)  = hlp1 / molar_mass_C / vgrid_soil_sb%dz(1)
             sb_pool%polymeric_litter%carbon(:,1) = hlp2 / molar_mass_C / vgrid_soil_sb%dz(1)
             IF(lctlib%growthform==itree) THEN
                sb_pool%woody_litter%carbon(:,1) = hlp3 / molar_mass_C / vgrid_soil_sb%dz(1)
             ENDIF
             sb_pool%microbial%carbon(:,1) = hlp4 / molar_mass_C / vgrid_soil_sb%dz(1)
             sb_pool%residue%carbon(:,1)   = hlp5 / molar_mass_C / vgrid_soil_sb%dz(1)

             nh4_solute(:,1) = nh4_solute_prescribe / vgrid_soil_sb%dz(1)
             no3_solute(:,1) = no3_solute_prescribe / vgrid_soil_sb%dz(1)

             ! inital Qmax based on current soil carbon content
             CALL calc_bulk_soil_carbon(nc, nsoil_sb, sb_pool, bulk_soil_carbon_sl, soil_litter_carbon_sl, &
                                        volume_min_sl, TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
             bulk_dens_corr_sl(:,:) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:),soil_litter_carbon_sl(:,:), &
                                                                   bulk_dens_sl(:,:))
             !> po4 adsorption parameter intialisation 
             qmax_po4_min_sl(:,:)     = calc_qmax_texture(qmax_po4_mineral, 0.0_wp, 1._wp)
             qmax_po4_om_sl(:,:)      = qmax_po4_min_sl(:,:) / 100._wp
             qmax_po4(:,:)            = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), &
                                                                          bulk_dens_sl(:,:), &
                                                                          qmax_po4_om_sl(:,:),qmax_po4_min_sl(:,:))
             qmax_po4(:,:)            = qmax_po4(:,:) * w_soil_sl(:,:) / soil_depth_sl(:,:)
             km_adsorpt_po4_sl(:,:)   = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), &
                                                                          bulk_dens_sl(:,:), &
                                                                          km_adsorpt_OM_po4,km_adsorpt_mineral_po4)
             !> nh4 adsorption parameter intialisation 
             qmax_nh4_min_sl(:,:)     = calc_qmax_texture(qmax_nh4_clay, 0.0_wp, clay_sl(:,:))
             qmax_nh4(:,:)            = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), &
                                                                          bulk_dens_sl(:,:), &
                                                                          0.0_wp, qmax_nh4_min_sl(:,:))
             qmax_nh4(:,:)            = qmax_nh4(:,:) * w_soil_sl(:,:) / soil_depth_sl(:,:)
             km_adsorpt_nh4_sl(:,:)   = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), &
                                                                          bulk_dens_sl(:,:), &
                                                                          km_adsorpt_OM_nh4, km_adsorpt_mineral_nh4)

             root_depth_avg(:) = SUM(root_fraction_sl(:,:) * soil_depth_sl(:,:), DIM=2)
             arr_hlp6(:)       = SUM(km_adsorpt_nh4_sl(:,:) * root_fraction_sl(:,:), DIM=2) / root_depth_avg(:) 
             CALL convert_s1d_to_s0d_var_soilproperty(nc, nsoil_sb, soil_depth_sl(:,:), qmax_nh4(:,:))

             nh4_assoc(:,1)           = qmax_nh4(:,1) * nh4_solute(:,1) / (nh4_solute(:,1) + arr_hlp6(:))
             nh4_n15_assoc(:,1)       = nh4_assoc(:,1) / (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 

             ! initialise mineral soil P pools
             CALL init_soil_pi_pools(nc, nsoil_sb, & 
                                    TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &
                                    soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), root_fraction_sl(:,:), &
                                    dsl4jsb_Config(SB_)%soil_p_labile, dsl4jsb_Config(SB_)%soil_p_slow, &
                                    dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
                                    dsl4jsb_Config(SB_)%soil_p_depth, & 
                                    po4_solute(:,:), po4_assoc_fast(:,:), po4_assoc_slow(:,:), &
                                    po4_occluded(:,:), po4_primary(:,:), &
                                    qmax_po4(:,:),km_adsorpt_po4_sl(:,:))
          ELSE 
             sb_pool%soluable_litter%carbon(:,1)  = hlp1 /  molar_mass_C / vgrid_soil_sb%dz(1) * &
                                                    (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1)) 
             sb_pool%polymeric_litter%carbon(:,1) = hlp2 / molar_mass_C / vgrid_soil_sb%dz(1) * &
                                                    (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1)) 
             IF(lctlib%growthform==itree) THEN
               sb_pool%woody_litter%carbon(:,1) = hlp3 / molar_mass_C / vgrid_soil_sb%dz(1)
             ENDIF
             sb_pool%microbial%carbon(:,1) = hlp4 / molar_mass_C / vgrid_soil_sb%dz(1) * &
                                             (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1)) 
             sb_pool%residue%carbon(:,1)   = hlp5 / molar_mass_C / vgrid_soil_sb%dz(1) * & 
                                             (frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1)) 
             DO isoil=2,nsoil_sb
                sb_pool%soluable_litter%carbon(:,isoil)  = hlp1 /  molar_mass_C / vgrid_soil_sb%dz(isoil) * &
                                                           (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil) 
                sb_pool%polymeric_litter%carbon(:,isoil) = hlp2 / molar_mass_C / vgrid_soil_sb%dz(isoil) * &
                                                           (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil) 
                sb_pool%microbial%carbon(:,isoil)        = hlp4 / molar_mass_C / vgrid_soil_sb%dz(isoil) * &
                                                           (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil) 
                sb_pool%residue%carbon(:,isoil)          = hlp5 / molar_mass_C / vgrid_soil_sb%dz(isoil) * & 
                                                           (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil)
             ENDDO
             ! note: this differs on purpose from the calculation in mo_sb_update_pools for the prescribed option
             !       this is because here I'm trying to enure that the models are initialised with the same amount of nutrient
             !       whereas in the prescribed option the target is unlimiting concentrations across the profile. The latter leads
             !       to more nutrient in the soil column in 1d compared to 0d
             nh4_solute(:,:) = nh4_solute_prescribe / SUM(vgrid_soil_sb%dz(:))
             no3_solute(:,:) = no3_solute_prescribe / SUM(vgrid_soil_sb%dz(:))

             ! inital bulk density and qmax_po4 based on current soil carbon content
             CALL calc_bulk_soil_carbon(nc, nsoil_sb, sb_pool, bulk_soil_carbon_sl, soil_litter_carbon_sl, &
                                        volume_min_sl, TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
             bulk_dens_corr_sl(:,:) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:),soil_litter_carbon_sl(:,:), &
                                                                   bulk_dens_sl(:,:))      
             !> po4 adsorption parameter intialisation 
             qmax_po4_min_sl(:,:) = calc_qmax_texture(qmax_po4_mineral, 0.0_wp, 1._wp)
             qmax_po4_om_sl(:,:)  = qmax_po4_min_sl(:,:) / 100._wp
             qmax_po4(:,:) = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:),bulk_dens_sl(:,:), & 
                                                               qmax_po4_om_sl(:,:),qmax_po4_min_sl(:,:))
             qmax_po4(:,:) = qmax_po4(:,:) * w_soil_sl(:,:) / soil_depth_sl(:,:)
             km_adsorpt_po4_sl(:,:) = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), &
                                                                        bulk_dens_sl(:,:), &
                                                                        km_adsorpt_OM_po4,km_adsorpt_mineral_po4)

             !> nh4 adsorption parameter intialisation 
             qmax_nh4_min_sl(:,:)     = calc_qmax_texture(qmax_nh4_clay, 0.0_wp, clay_sl(:,:))
             qmax_nh4(:,:)            = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), &
                                                                          bulk_dens_sl(:,:), &
                                                                          0.0_wp,qmax_nh4_min_sl(:,:))
             qmax_nh4(:,:)            = qmax_nh4(:,:) * w_soil_sl(:,:) / soil_depth_sl(:,:)
             km_adsorpt_nh4_sl(:,:)   = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), &
                                                                          bulk_dens_sl(:,:), &
                                                                          km_adsorpt_OM_nh4,km_adsorpt_mineral_nh4)
             nh4_assoc(:,:)           = qmax_nh4(:,:) * nh4_solute (:,:) / (nh4_solute(:,:) + km_adsorpt_nh4_sl(:,:))
             nh4_n15_assoc(:,:)       = nh4_assoc(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
          ENDIF
 
          ! metabolic litter all other elements
          sb_pool%soluable_litter%nitrogen(:,:)   = sb_pool%soluable_litter%carbon(:,:) / k_init_soluable_litter_cn 
          sb_pool%soluable_litter%phosphorus(:,:) = sb_pool%soluable_litter%nitrogen(:,:) / k_init_litter_np 
          sb_pool%soluable_litter%carbon13(:,:)   = sb_pool%soluable_litter%carbon(:,:) * frac_c13
          sb_pool%soluable_litter%carbon14(:,:)   = sb_pool%soluable_litter%carbon(:,:) * frac_c14 
          sb_pool%soluable_litter%nitrogen15(:,:) = sb_pool%soluable_litter%nitrogen(:,:) / & 
                                                    (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! non-woody litter all other elements
          sb_pool%polymeric_litter%nitrogen(:,:)   = sb_pool%polymeric_litter%carbon(:,:) / k_init_polymeric_litter_cn 
          sb_pool%polymeric_litter%phosphorus(:,:) = sb_pool%polymeric_litter%nitrogen(:,:) / k_init_litter_np 
          sb_pool%polymeric_litter%carbon13(:,:)   = sb_pool%polymeric_litter%carbon(:,:) * frac_c13
          sb_pool%polymeric_litter%carbon14(:,:)   = sb_pool%polymeric_litter%carbon(:,:) * frac_c14
          sb_pool%polymeric_litter%nitrogen15(:,:) = sb_pool%polymeric_litter%nitrogen(:,:) / & 
                                                     (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! woody litter all other elements
          sb_pool%woody_litter%nitrogen(:,:)   = sb_pool%woody_litter%carbon(:,:) / k_init_woody_litter_cn
          sb_pool%woody_litter%phosphorus(:,:) = sb_pool%woody_litter%nitrogen(:,:) / k_init_litter_np 
          sb_pool%woody_litter%carbon13(:,:)   = sb_pool%woody_litter%carbon(:,:) * frac_c13
          sb_pool%woody_litter%carbon14(:,:)   = sb_pool%woody_litter%carbon(:,:) * frac_c14
          sb_pool%woody_litter%nitrogen15(:,:) = sb_pool%woody_litter%nitrogen(:,:) / & 
                                                 (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          ! fast pool (=microbial) all other elements
          sb_pool%microbial%nitrogen(:,:)   = sb_pool%microbial%carbon(:,:) / (k_fast_som_cn_min + k_fast_som_cn_max) * 2.0_wp 
          sb_pool%microbial%phosphorus(:,:) = sb_pool%microbial%nitrogen(:,:) / k_fast_som_np
          sb_pool%microbial%carbon13(:,:)   = sb_pool%microbial%carbon(:,:) * frac_c13
          sb_pool%microbial%carbon14(:,:)   = sb_pool%microbial%carbon(:,:) * frac_c14
          sb_pool%microbial%nitrogen15(:,:) = sb_pool%microbial%nitrogen(:,:) / & 
                                              (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          sb_pool%residue%nitrogen(:,:)   = sb_pool%residue%carbon(:,:) / k_slow_som_cn
          sb_pool%residue%phosphorus(:,:) = sb_pool%residue%nitrogen(:,:) / k_slow_som_np 
          sb_pool%residue%carbon13(:,:)   = sb_pool%residue%carbon(:,:) * frac_c13
          sb_pool%residue%carbon14(:,:)   = sb_pool%residue%carbon(:,:) * frac_c14 
          sb_pool%residue%nitrogen15(:,:) = sb_pool%residue%nitrogen(:,:) / & 
                                            (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          nh4_n15_solute(:,:) = nh4_solute(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
          no3_n15_solute(:,:) = no3_solute(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
          
          ! SPP site-set
          IF (flag_spp1685) THEN
            ! for SPP sites, read in initial state from data
            CALL sb_init_SPP(nc, nsoil_sb, soil_depth_acc_sl(:,:), & 
                          dsl4jsb_Config(SB_)%soil_p_labile, dsl4jsb_Config(SB_)%soil_p_slow, &
                          dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
                          labile_p(:,:),po4_assoc_slow(:,:),po4_occluded(:,:), po4_primary(:,:), &
                          site_ID_spp1685)

            po4_assoc_fast(:,:) = calc_fast_po4(labile_P(:,:),qmax_po4(:,:),km_adsorpt_po4_sl(:,:))
            po4_solute(:,:) = labile_P(:,:) - po4_assoc_fast(:,:)
            IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_0d") THEN
              po4_assoc_slow(:,1) = SUM(po4_assoc_slow(:,:) * soil_depth_sl(:,:)* root_fraction_sl(:,:),DIM=2) / &
                                    root_depth_avg(:) / soil_depth_sl(:,1)
              po4_occluded(:,1)   = SUM(po4_occluded(:,:)   * soil_depth_sl(:,:)* root_fraction_sl(:,:),DIM=2) / &
                                    root_depth_avg(:) / soil_depth_sl(:,1)
              po4_primary(:,1)    = SUM(po4_primary(:,:)    * soil_depth_sl(:,:)* root_fraction_sl(:,:),DIM=2) / &
                                    root_depth_avg(:) / soil_depth_sl(:,1)
              po4_assoc_fast(:,1) = SUM(po4_assoc_fast(:,:) * soil_depth_sl(:,:)* root_fraction_sl(:,:),DIM=2) / &
                                    root_depth_avg(:) / soil_depth_sl(:,1)
              po4_solute(:,1)     = SUM(po4_solute(:,:)     * soil_depth_sl(:,:)* root_fraction_sl(:,:),DIM=2) / &
                                    root_depth_avg(:) / soil_depth_sl(:,1)
            ENDIF
          ELSE
          ! all sites but the SPP site-set
            ! initialise mineral soil P pools
            CALL init_soil_pi_pools(nc, nsoil_sb, & 
                                   TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &
                                   soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), root_fraction_sl(:,:), &
                                   dsl4jsb_Config(SB_)%soil_p_labile, dsl4jsb_Config(SB_)%soil_p_slow, &
                                   dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
                                   dsl4jsb_Config(SB_)%soil_p_depth, & 
                                   po4_solute(:,:), po4_assoc_fast(:,:), po4_assoc_slow(:,:), &
                                   po4_occluded(:,:), po4_primary(:,:), &
                                   qmax_po4(:,:),km_adsorpt_po4_sl(:,:))
          ENDIF

       !!==========================================================================================================================
       !> JSM (Jena Soil Model)
       !!
       CASE('jsm')

          hlp1 = 100._wp   ! g C m-2 metabolic litter
          hlp2 = 1000._wp  ! g C m-2 structural litter
          hlp3 = 1000._wp  ! g C m-2 woody litter
          hlp4 = 30._wp    ! g C m-2 microbial biomass, 30
          hlp5 = 150._wp   ! g C m-2 microbial necromass, 150
          hlp7 = 1._wp     ! g C m-2 DOM, 1

          hlp8 = 0.1_wp           ! numerical coefficient for correction, 0.1_wp
          hlp9 = 5._wp            ! numerical coefficient for correction, 5._wp
          hlp10 = 1._wp / 3._wp    ! re-partition factor when OM got over-saturation in intialisation
          hlp11 = 1._wp / 11._wp  ! re-partition factor when OM got over-saturation in intialisation

          sb_pool%soluable_litter%carbon(:,1)  = hlp1 /  molar_mass_C / vgrid_soil_sb%dz(1) * &
                                                   ( frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1)) 
          sb_pool%polymeric_litter%carbon(:,1) = hlp2 / molar_mass_C / vgrid_soil_sb%dz(1) * &
                                                   ( frac_prod_not_root + (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,1)) 
         IF(lctlib%growthform==itree) THEN
           sb_pool%woody_litter%carbon(:,1)  = hlp3 / molar_mass_C / vgrid_soil_sb%dz(1)
         ENDIF
         DO isoil=2,nsoil_sb
             sb_pool%soluable_litter%carbon(:,isoil)  = hlp1 /  molar_mass_C / vgrid_soil_sb%dz(isoil) * &
                                                          (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil) 
             sb_pool%polymeric_litter%carbon(:,isoil) = hlp2 / molar_mass_C / vgrid_soil_sb%dz(isoil) * &
                                                          (1.0_wp - frac_prod_not_root) * root_fraction_sl(:,isoil) 
         ENDDO

         DO isoil=1,nsoil_sb
             sb_pool%microbial%carbon(:,isoil)  = MAX(hlp4/vgrid_soil_sb%dz(isoil) * root_fraction_sl(:,isoil), &
                                                      min_mic_biomass)
             sb_pool%dom%carbon(:,isoil)        = hlp7 /vgrid_soil_sb%dz(isoil) * root_fraction_sl(:,isoil)
             sb_pool%dom_assoc%carbon(:,isoil)  = hlp9 * sb_pool%dom%carbon(:,isoil)
             sb_pool%residue%carbon(:,isoil)    = hlp5/vgrid_soil_sb%dz(isoil) * root_fraction_sl(:,isoil)
             frac_prod_not_root = 1._wp / isoil
             ! deep soil is initialised with very low MOC, corrected by hlp8
             IF (frac_prod_not_root < hlp8) THEN 
               frac_prod_not_root = frac_prod_not_root * hlp8
             ENDIF
             sb_pool%residue_assoc%carbon(:,isoil) = dsl4jsb_Config(SB_)%qmax_org_fine_particle * &
                                                     (clay_sl(:,isoil) + silt_sl(:,isoil)) * & 
                                                     bulk_dens_sl(:,isoil) * frac_prod_not_root
             ! To avoid over-saturation of OM in soil layer
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

             WHERE (arr_hlp7(:,isoil) > 1._wp )
                sb_pool%soluable_litter%carbon(:,isoil)  = sb_pool%soluable_litter%carbon(:,isoil)  * 1._wp / arr_hlp7(:,isoil)
                sb_pool%polymeric_litter%carbon(:,isoil) = sb_pool%polymeric_litter%carbon(:,isoil) * 1._wp / arr_hlp7(:,isoil)
                sb_pool%microbial%carbon(:,isoil)        = sb_pool%microbial%carbon(:,isoil)        * 1._wp / arr_hlp7(:,isoil)
                sb_pool%dom%carbon(:,isoil)              = sb_pool%dom%carbon(:,isoil)              * 1._wp / arr_hlp7(:,isoil)
                sb_pool%dom_assoc%carbon(:,isoil)        = sb_pool%dom_assoc%carbon(:,isoil)        * 1._wp / arr_hlp7(:,isoil)
                sb_pool%residue%carbon(:,isoil)          = 0._wp
                sb_pool%residue_assoc%carbon(:,isoil)    = 0._wp
             ELSEWHERE ((arr_hlp7(:,isoil) + arr_hlp8(:,isoil)) > 1._wp )
                arr_hlp9(:,isoil) = (1._wp - arr_hlp7(:,isoil)) * hlp10 * dsl4jsb_Config(SB_)%qmax_org_fine_particle * &
                                    (clay_sl(:,isoil) + silt_sl(:,isoil)) * bulk_dens_sl(:,isoil)
                arr_hlp8(:,isoil) = (1._wp - arr_hlp7(:,isoil)) * (1._wp - hlp10) * &
                                    rho_bulk_org * 1000._wp * carbon_per_dryweight_SOM / molar_mass_C
                sb_pool%residue_assoc%carbon(:,isoil) = calc_fast_po4(arr_hlp8(:,isoil), arr_hlp9(:,isoil), &
                                                                      k_desorpt_det / k_adsorpt_det)
                sb_pool%residue%carbon(:,isoil)       = arr_hlp8(:,isoil) - sb_pool%residue_assoc%carbon(:,isoil)
             ELSEWHERE
                arr_hlp9(:,isoil) = (1._wp - arr_hlp7(:,isoil) - arr_hlp8(:,isoil)) * dsl4jsb_Config(SB_)%qmax_org_fine_particle * &
                                    (clay_sl(:,isoil) + silt_sl(:,isoil)) * bulk_dens_sl(:,isoil)
                arr_hlp8(:,isoil) = (sb_pool%residue%carbon(:,isoil) + sb_pool%residue_assoc%carbon(:,isoil))
                WHERE (arr_hlp8(:,isoil) > arr_hlp9(:,isoil))
                  sb_pool%soluable_litter%carbon(:,isoil) = sb_pool%soluable_litter%carbon(:,isoil) + &
                                                            (arr_hlp8(:,isoil)-arr_hlp9(:,isoil)) * hlp11
                  sb_pool%polymeric_litter%carbon(:,isoil)= sb_pool%polymeric_litter%carbon(:,isoil) + &
                                                            (arr_hlp8(:,isoil)-arr_hlp9(:,isoil)) * (1._wp - hlp11)
                  arr_hlp8(:,isoil) = arr_hlp9(:,isoil)
                ENDWHERE
                sb_pool%residue_assoc%carbon(:,isoil) = calc_fast_po4(arr_hlp8(:,isoil), arr_hlp9(:,isoil), &
                                                                      k_desorpt_det / k_adsorpt_det)
                sb_pool%residue%carbon(:,isoil)       = arr_hlp8(:,isoil) - sb_pool%residue_assoc%carbon(:,isoil)
             ENDWHERE
         ENDDO
     
         ! initial values
         nh4_solute(:,:)      = nh4_solute_prescribe / SUM(vgrid_soil_sb%dz(:))
         nh4_n15_solute(:,:)  = nh4_solute(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
         no3_solute(:,:)      = no3_solute_prescribe / SUM(vgrid_soil_sb%dz(:))
         no3_n15_solute(:,:)  = no3_solute(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
 
         ! inital bulk density, Qmax_org & Qmax_po4 based on current soil carbon content
         CALL calc_bulk_soil_carbon(nc, nsoil_sb, sb_pool, bulk_soil_carbon_sl, soil_litter_carbon_sl, &
                                    volume_min_sl, TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
         bulk_dens_corr_sl(:,:) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:), soil_litter_carbon_sl(:,:), &
                                                               bulk_dens_sl(:,:))
         !> po4/OC adsorption parameter intialisation 
         qmax_org_min_sl(:,:) = calc_qmax_texture(dsl4jsb_Config(SB_)%qmax_org_fine_particle, silt_sl(:,:), clay_sl(:,:))
         qmax_po4_min_sl(:,:) = calc_qmax_texture(qmax_po4_mineral, 0.0_wp, 1._wp)
         qmax_po4_om_sl(:,:) = qmax_po4_min_sl(:,:)* hlp8  ! assuming the P sorption capacity of OM is one-tenth of mineral soil
         qmax_org(:,:) = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:),bulk_dens_sl(:,:), &
                                                           0._wp, qmax_org_min_sl(:,:))

         qmax_po4(:,:) = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:),bulk_dens_sl(:,:), &
                                                           qmax_po4_om_sl(:,:),qmax_po4_min_sl(:,:))
         qmax_po4(:,:) = qmax_po4(:,:) * w_soil_sl(:,:) / soil_depth_sl(:,:)
         km_adsorpt_po4_sl(:,:) = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), &
                                                                    bulk_dens_sl(:,:), km_adsorpt_OM_po4, km_adsorpt_mineral_po4)

         !> nh4 adsorption parameter intialisation 
         qmax_nh4_min_sl(:,:)     = calc_qmax_texture(qmax_nh4_clay, 0.0_wp, clay_sl(:,:))
         qmax_nh4(:,:)            = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), &
                                                                      bulk_dens_sl(:,:), &
                                                                      0.0_wp,qmax_nh4_min_sl(:,:))
         qmax_nh4(:,:)            = qmax_nh4(:,:) * w_soil_sl(:,:) / soil_depth_sl(:,:)
         km_adsorpt_nh4_sl(:,:)   = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), &
                                                                      bulk_dens_sl(:,:), &
                                                                      km_adsorpt_OM_nh4,km_adsorpt_mineral_nh4)
         nh4_assoc(:,:)           = qmax_nh4(:,:) * nh4_solute (:,:) / (nh4_solute(:,:) + km_adsorpt_nh4_sl(:,:))
         nh4_n15_assoc(:,:)       = nh4_assoc(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
        ! initialise mineral soil P pools
        DO isoil = 1,nsoil_sb
          WHERE(soil_depth_sl(:,isoil) > eps8)
            fine_root_carbon_sl(:,isoil) = root_fraction_sl(:,isoil) / soil_depth_sl(:,isoil) * & 
                                           veg_pool%fine_root%carbon(:)
          ENDWHERE
          IF (isoil == 1) THEN
            soil_depth_acc_sl(:,isoil) = soil_depth_sl(:,isoil) / 2._wp
          ELSE
            soil_depth_acc_sl(:,isoil) = soil_depth_acc_sl(:,isoil-1) + (soil_depth_sl(:,isoil) + soil_depth_sl(:,isoil-1)) / 2._wp  
          ENDIF   
        ENDDO

       IF (flag_spp1685) THEN
         ! for SPP sites, read in initial state from data
         CALL sb_init_SPP(nc, nsoil_sb, soil_depth_acc_sl(:,:), & 
                          dsl4jsb_Config(SB_)%soil_p_labile, dsl4jsb_Config(SB_)%soil_p_slow, &
                          dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
                          labile_p(:,:),po4_assoc_slow(:,:),po4_occluded(:,:), po4_primary(:,:), &
                          site_ID_spp1685)

          po4_assoc_fast(:,:) = calc_fast_po4(labile_P(:,:),qmax_po4(:,:),km_adsorpt_po4_sl(:,:))
          po4_solute(:,:) = labile_P(:,:) - po4_assoc_fast(:,:)
        ELSE
          CALL init_soil_pi_pools(nc, nsoil_sb, &
                                  TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &
                                  soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), root_fraction_sl(:,:), &
                                  dsl4jsb_Config(SB_)%soil_p_labile, dsl4jsb_Config(SB_)%soil_p_slow, &
                                  dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
                                  dsl4jsb_Config(SB_)%soil_p_depth, & 
                                  po4_solute(:,:), po4_assoc_fast(:,:), po4_assoc_slow(:,:), &
                                  po4_occluded(:,:), po4_primary(:,:), &
                                  qmax_po4(:,:),km_adsorpt_po4_sl(:,:))
        ENDIF

        ! initialise diffusion by bioturbation
        k_bioturb(:,:) = calc_bioturbation_rate(nc, nsoil_sb,bulk_dens_corr_sl(:,:), root_fraction_sl(:,:), soil_depth_sl(:,:)) 
        
        ! enzyme allocation initialise 
        enzyme_frac_poly(:,:)     = 0.5_wp 
        enzyme_frac_residue(:,:)  = 0.5_wp 
        enzyme_frac_poly_c(:,:)   = 0.5_wp 
        enzyme_frac_poly_n(:,:)   = 0.5_wp 
        enzyme_frac_poly_p(:,:)   = 0.5_wp 

        ! isotopic signals in all pools
        ! microbial, DOM, aDOM, residue, and aRes pools
        sb_pool%soluable_litter%nitrogen(:,:)     = sb_pool%soluable_litter%carbon(:,:) / k_init_soluable_litter_cn 
        sb_pool%soluable_litter%phosphorus(:,:)   = sb_pool%soluable_litter%nitrogen(:,:) / k_init_litter_np 
        sb_pool%soluable_litter%carbon13(:,:)     = sb_pool%soluable_litter%carbon(:,:) * frac_c13 
        sb_pool%soluable_litter%carbon14(:,:)     = sb_pool%soluable_litter%carbon(:,:) * frac_c14 
        sb_pool%soluable_litter%nitrogen15(:,:)   = sb_pool%soluable_litter%nitrogen(:,:) / &
                                                      ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        sb_pool%polymeric_litter%nitrogen(:,:)    = sb_pool%polymeric_litter%carbon(:,:) / k_init_polymeric_litter_cn
        sb_pool%polymeric_litter%phosphorus(:,:)  = sb_pool%polymeric_litter%nitrogen(:,:) / k_init_litter_np 
        sb_pool%polymeric_litter%carbon13(:,:)    = sb_pool%polymeric_litter%carbon(:,:) * frac_c13 
        sb_pool%polymeric_litter%carbon14(:,:)    = sb_pool%polymeric_litter%carbon(:,:) * frac_c14 
        sb_pool%polymeric_litter%nitrogen15(:,:)  = sb_pool%polymeric_litter%nitrogen(:,:) / &
                                                      ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        sb_pool%woody_litter%nitrogen(:,:)        = sb_pool%woody_litter%carbon(:,:) / k_init_woody_litter_cn
        sb_pool%woody_litter%phosphorus(:,:)      = sb_pool%woody_litter%nitrogen(:,:) / k_init_litter_np
        sb_pool%woody_litter%carbon13(:,:)        = sb_pool%woody_litter%carbon(:,:) * frac_c13 
        sb_pool%woody_litter%carbon14(:,:)        = sb_pool%woody_litter%carbon(:,:) * frac_c14 
        sb_pool%woody_litter%nitrogen15(:,:)      = sb_pool%woody_litter%nitrogen(:,:) / &
                                                      ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        sb_pool%microbial%nitrogen(:,:)           = sb_pool%microbial%carbon(:,:) / microbial_cn 
        sb_pool%microbial%phosphorus(:,:)         = sb_pool%microbial%nitrogen(:,:) / microbial_np
        sb_pool%microbial%carbon13(:,:)           = sb_pool%microbial%carbon(:,:) * frac_c13 
        sb_pool%microbial%carbon14(:,:)           = sb_pool%microbial%carbon(:,:) * frac_c14 
        sb_pool%microbial%nitrogen15(:,:)         = sb_pool%microbial%nitrogen(:,:) / &
                                                      ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        sb_pool%residue%nitrogen(:,:)             = sb_pool%residue%carbon(:,:) / microbial_cn 
        sb_pool%residue%phosphorus(:,:)           = sb_pool%residue%nitrogen(:,:) / k_slow_som_np
        sb_pool%residue%carbon13(:,:)             = sb_pool%residue%carbon(:,:) * frac_c13 
        sb_pool%residue%carbon14(:,:)             = sb_pool%residue%carbon(:,:) * frac_c14 
        sb_pool%residue%nitrogen15(:,:)           = sb_pool%residue%nitrogen(:,:) / &
                                                      ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        sb_pool%residue_assoc%nitrogen(:,:)       = sb_pool%residue_assoc%carbon(:,:) / microbial_cn
        sb_pool%residue_assoc%phosphorus(:,:)     = sb_pool%residue_assoc%nitrogen(:,:) / k_slow_som_np
        sb_pool%residue_assoc%carbon13(:,:)       = sb_pool%residue_assoc%carbon(:,:) * frac_c13 
        sb_pool%residue_assoc%carbon14(:,:)       = sb_pool%residue_assoc%carbon(:,:) * frac_c14 
        sb_pool%residue_assoc%nitrogen15(:,:)     = sb_pool%residue_assoc%nitrogen(:,:) / &
                                                      ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

        sb_pool%dom%nitrogen(:,:)           = sb_pool%dom%carbon(:,:) / microbial_cn
        sb_pool%dom%phosphorus(:,:)         = sb_pool%dom%nitrogen(:,:) / microbial_np
        sb_pool%dom%carbon13(:,:)           = sb_pool%dom%carbon(:,:) * frac_c13
        sb_pool%dom%carbon14(:,:)           = sb_pool%dom%carbon(:,:) * frac_c14
        sb_pool%dom%nitrogen15(:,:)         = sb_pool%dom%nitrogen(:,:) / &
                                                (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

        sb_pool%dom_assoc%nitrogen(:,:)     = sb_pool%dom_assoc%carbon(:,:) / microbial_cn
        sb_pool%dom_assoc%phosphorus(:,:)   = sb_pool%dom_assoc%nitrogen(:,:) / microbial_np
        sb_pool%dom_assoc%carbon13(:,:)     = sb_pool%dom_assoc%carbon(:,:) * frac_c13
        sb_pool%dom_assoc%carbon14(:,:)     = sb_pool%dom_assoc%carbon(:,:) * frac_c14
        sb_pool%dom_assoc%nitrogen15(:,:)   = sb_pool%dom_assoc%nitrogen(:,:) / & 
                                                (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        ! fungi and mycorrhiza
        sb_pool%fungi%carbon13(:,:)         = sb_pool%fungi%carbon(:,:) * frac_c13
        sb_pool%fungi%carbon14(:,:)         = sb_pool%fungi%carbon(:,:) * frac_c14 
        sb_pool%fungi%nitrogen15(:,:)       = sb_pool%fungi%nitrogen(:,:) / & 
                                                (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        sb_pool%mycorrhiza%carbon13(:,:)    = sb_pool%mycorrhiza%carbon(:,:) * frac_c13
        sb_pool%mycorrhiza%carbon14(:,:)    = sb_pool%mycorrhiza%carbon(:,:) * frac_c14 
        sb_pool%mycorrhiza%nitrogen15(:,:)  = sb_pool%mycorrhiza%nitrogen(:,:) / & 
                                                (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
       CASE DEFAULT
         CALL finish(routine,'Invalid setting of sb_model_scheme.')
 
    END SELECT


    !> 3.0 other SB_ variables - all soil models
    !! 
    microbial_cue_eff_tmic_mavg(:,:)  = microbial_cue_max
    microbial_nue_eff_tmic_mavg(:,:)  = microbial_nue
    microbial_pue_eff_tmic_mavg(:,:)  = microbial_pue
    vmax_weath_mineral_sl(:,:)        = k_weath_mineral 
    dom_cn(:,:)                       = microbial_cn
    dom_cp(:,:)                       = microbial_cn * microbial_np

  ! Deallocate local allocatable variables
  DEALLOCATE(fine_root_carbon_sl)
  DEALLOCATE(labile_p)
  DEALLOCATE(soil_depth_acc_sl)
  DEALLOCATE(arr_hlp6)
  DEALLOCATE(arr_hlp7)
  DEALLOCATE(arr_hlp8)
  DEALLOCATE(arr_hlp9)
  DEALLOCATE(root_depth_avg)


! Pass back values to process memory
dsl4slm_passback_vegpool_var(pool,fine_root,carbon)
#define Passback_sb_pool
#include "slm_pools_a_incl.f90"
#undef Passback_sb_pool

  END SUBROUTINE sb_init

  !-----------------------------------------------------------------------------------------------------
  !> Intialize soil pi pools for all layers from the soil P maps (Yang et al. 2013)
  !! assuming that the soil pi pool sizes are correlated with soil weight
  !! 
  !! also used for simple_0d
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_soil_pi_pools(nc, nsoil_sb, sb_model_scheme, sb_adsorp_scheme, & 
                               soil_depth_sl, bulk_dens_sl, root_fraction_sl,&
                               soil_p_labile, soil_p_slow, soil_p_occluded, soil_p_primary, soil_p_depth, &
                               po4_solute, po4_assoc_fast, po4_assoc_slow, po4_occluded, po4_primary, &
                               qmax_po4,km)

    USE mo_jsb_phys_constants,    ONLY: molar_mass_P
    USE mo_sb_constants
    USE mo_sb_jsm_processes,      ONLY: calc_fast_po4, convert_s1d_to_s0d_var_soilproperty

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)     :: nc                 !< dimensions
    INTEGER                                         , INTENT(in)     :: nsoil_sb           ! number of soil layers as used/defined by the SB_ process
    CHARACTER(len=*)                                , INTENT(in)     :: sb_model_scheme    !< soil model scheme 
    CHARACTER(len=*)                                , INTENT(in)     :: sb_adsorp_scheme   !< soil model po4 adsorption scheme 
    REAL(wp), DIMENSION(nc,nsoil_sb)                , INTENT(in)     :: soil_depth_sl, &   !< soil depth per soil layer, [m]
                                                                        root_fraction_sl, &
                                                                        bulk_dens_sl, &    !< corrected soil bulk density per soil layer, [kg m-3]
                                                                        qmax_po4, &
                                                                        km
    REAL(wp)                                        , INTENT(in)     :: soil_p_depth, &    !< soil depth of the soil P maps
                                                                        soil_p_labile, &   !< total labile P pool (solute + assoc_fast) of the soil P maps, [g m-2]
                                                                        soil_p_slow, &     !< total slow P pool (secondary P) of the soil P maps, [g m-2]
                                                                        soil_p_occluded, & !< total occluded P pool of the soil P maps, [g m-2]
                                                                        soil_p_primary     !< total primary P pool of the soil P maps, [g m-2]
    REAL(wp), DIMENSION(nc,nsoil_sb)                , INTENT(inout)  :: po4_solute, &      !< soluble Pi pool per soil layer, [mol m-3]
                                                                        po4_assoc_fast, &  !< fast minerally associated PO4 pool per soil layer, [mol m-3]
                                                                        po4_assoc_slow, &  !< slow minerally associated PO4 pool per soil layer, [mol m-3]
                                                                        po4_occluded, &    !< occluded PO4 pool per soil layer, [mol m-3]
                                                                        po4_primary        !< primary PO4 pool per soil layer, [mol m-3]
    ! ---------------------------
    ! 0.2 Local
    INTEGER                                          :: i
    REAL(wp)                                         :: rhlp1,rhlp2,rhlp3,rhlp4,rhlp5,rhlp6                           !< logical helper for input Pi density calculation
    REAL(wp), DIMENSION(nc)                          :: soil_p_weight, soil_weight_0d, &  !< total soil weight up to the soil_p_depth/soil_depth_0d, [kg]
                                                        soil_p_layer                      !< up to soil_p_layer keep the P init map info
    REAL(wp), DIMENSION(nc)                          :: hlp1, hlp2, hlp3, hlp4
    REAL(wp), DIMENSION(nc)                          :: soil_depth_0d
    LOGICAL,  DIMENSION(nc)                          :: aflag  
    REAL(wp), DIMENSION(nc,nsoil_sb)                 :: soil_weight_sl, &                 !< soil weight per soil layer
                                                        arr_hlp1, arr_hlp2, soil_depth_acc_sl
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_soil_pi_pools'


    !> 0.9
    !! initialise output
    po4_solute(:,:)     = 0._wp
    po4_assoc_fast(:,:) = 0._wp
    po4_assoc_slow(:,:) = 0._wp 
    po4_occluded(:,:)   = 0._wp
    po4_primary(:,:)    = 0._wp

    !> 1.0 calculate the total soil weight up to the soil_p_depth,
    !! and the soil weight per soil layer
    hlp1(:) = 0._wp
    hlp2(:) = 0._wp
    aflag(:) = .TRUE.
    rhlp1 = soil_p_slow + soil_p_labile + soil_p_occluded + soil_p_primary
    DO i=1,nsoil_sb
      soil_weight_sl(:,i) = soil_depth_sl(:,i) * bulk_dens_sl(:,i)
      WHERE (aflag) 
        hlp1(:) = hlp1(:) + soil_depth_sl(:,i)
        hlp2(:) = hlp2(:) + soil_weight_sl(:,i)
        WHERE (hlp1(:) > soil_p_depth)
          soil_p_weight(:) = hlp2(:) - soil_weight_sl(:,i) * (hlp1(:) - soil_p_depth) / soil_depth_sl(:,i)
          aflag(:) = .FALSE.
          soil_p_layer(:) = i
        ENDWHERE
      ENDWHERE
    END DO
    hlp4(:) = rhlp1 / soil_p_weight(:) ! the Pi density from init map, [g P / kg soil]
    
    !> 2.0 calculate the initial slow P, occluded P, primary P pool, and labile_P
    !! based on the assumption that soil P pool are linearly corrected with soil weight
    DO i=1,nsoil_sb
      IF (i==1) THEN
        soil_depth_acc_sl(:,i)        = soil_depth_sl(:,i) / 2._wp
      ELSE
        soil_depth_acc_sl(:,i)          = soil_depth_acc_sl(:,i-1) + (soil_depth_sl(:,i) + soil_depth_sl(:,i-1)) / 2._wp  
      ENDIF   
    ENDDO
    !> the depth correction of initial Pi pool distributions are fitted against the observed pattern of a young cambisol soil 
    !! with basalt as parent material (BBR, Lang et al. 2017) that, 
    !! the proportion of: labile+sorbed P decreases with increasing depth, and primary P increases with increasing depth (Yu et al., 2020)
    rhlp2 = 3.5_wp            ! denominator coefficient for P pool depth correction
    rhlp3 = 2._wp             ! numerical coefficient for P pool depth correction
    rhlp4 = 0.1_wp            ! numerical coefficient of minimum occluded P for P pool depth correction
    rhlp5 = 0.8_wp            ! numerical coefficient for P pool depth correction
    rhlp6 = 9._wp / 17._wp   ! global average value of labile_P/(labile_P + slow_P)
    DO i=1,nsoil_sb
      hlp1(:)             = MIN(rhlp3 * soil_depth_acc_sl(:,i) ** 0.5_wp, (rhlp2 - rhlp4)) / rhlp2 ! primary P, corrected with depth, [g P] 
      hlp2(:)             = rhlp5 ** (rhlp3 * soil_depth_acc_sl(:,i)) / rhlp2
      hlp3(:)             = MAX((1.0_wp - hlp1(:) - hlp2(:)), rhlp4 / rhlp2) ! occluded P, corrected with depth

      WHERE ((hlp1(:) + hlp2(:) + hlp3(:)) > 1._wp)
        hlp1(:) = hlp1(:) / (hlp1(:) + hlp2(:) + hlp3(:))
        hlp2(:) = hlp2(:) / (hlp1(:) + hlp2(:) + hlp3(:))
        hlp3(:) = hlp3(:) / (hlp1(:) + hlp2(:) + hlp3(:))
      ENDWHERE
      WHERE (soil_p_layer(:) > i) 
        po4_assoc_slow(:,i) = bulk_dens_sl(:,i) * hlp4(:) * (soil_p_slow + soil_p_labile) / rhlp1 / molar_mass_P * (1._wp - rhlp6)
        po4_occluded(:,i)   = bulk_dens_sl(:,i) * hlp4(:) * soil_p_occluded               / rhlp1 / molar_mass_P
        po4_primary(:,i)    = bulk_dens_sl(:,i) * hlp4(:) * soil_p_primary                / rhlp1 / molar_mass_P
        arr_hlp1(:,i)       = bulk_dens_sl(:,i) * hlp4(:) * (soil_p_slow + soil_p_labile) / rhlp1 / molar_mass_P * rhlp6 !p_labile_sl in the follwing Eqn, [mol / m3]
      ELSEWHERE
        po4_assoc_slow(:,i) = bulk_dens_sl(:,i) * hlp4(:) * hlp2(:) / molar_mass_P * (1._wp - rhlp6)
        po4_occluded(:,i)   = bulk_dens_sl(:,i) * hlp4(:) * hlp3(:) / molar_mass_P
        po4_primary(:,i)    = bulk_dens_sl(:,i) * hlp4(:) * hlp1(:) / molar_mass_P
        arr_hlp1(:,i)       = bulk_dens_sl(:,i) * hlp4(:) * hlp2(:) / molar_mass_P * rhlp6 !p_labile_sl in the follwing Eqn, [mol / m3]
      ENDWHERE
    END DO

    !> 2.1 calculate the initial solute P and fast P,
    !! based on the assumption that labile P is the sum of both,
    !! and they follow the Langmuir equilibrium
    !! po4_assoc_fast = (p_labile_sl + Smax + K - sqrt((p_labile_sl + Smax + K)**2 - 4*p_labile_sl*Smax)) / 2
    DO i=1,nsoil_sb
      po4_assoc_fast(:,i) = calc_fast_po4(arr_hlp1(:,i),qmax_po4(:,i),km(:,i))
      po4_solute(:,i)     = arr_hlp1(:,i) - po4_assoc_fast(:,i)
    END DO

    !> 3.0 simple_0d should represent the first meter from the simple_1d
    !!
    IF (sb_model_scheme=='simple_0d') THEN
      CALL convert_s1d_to_s0d_var_soilproperty (nc, nsoil_sb, soil_depth_sl, po4_assoc_fast)
      CALL convert_s1d_to_s0d_var_soilproperty (nc, nsoil_sb, soil_depth_sl, po4_solute)
      CALL convert_s1d_to_s0d_var_soilproperty (nc, nsoil_sb, soil_depth_sl, po4_assoc_slow)
      CALL convert_s1d_to_s0d_var_soilproperty (nc, nsoil_sb, soil_depth_sl, po4_occluded)
      CALL convert_s1d_to_s0d_var_soilproperty (nc, nsoil_sb, soil_depth_sl, po4_primary)
    END IF

  END SUBROUTINE init_soil_pi_pools

  !-----------------------------------------------------------------------------------------------------
  !> Intialize soil_biogeochemistry process specifically for all SPP sites (after memory has been set up)
  !! site of the SPP site-set: BBR, MIT, VES, CON, LUE
  !! 
  !! @TODO this routine works only at one point! because soil properties for exactly one site are selected
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE sb_init_SPP (nc, nsoil_sb, soil_depth_acc_sl, &
                          soil_p_labile, soil_p_slow, soil_p_occluded, soil_p_primary, &
                          labile_p, po4_assoc_slow, po4_occluded, po4_primary, sitename)

    USE mo_sb_constants !,         ONLY: 

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER                         , INTENT(in)     :: nc                 !< dimensions
    INTEGER                         , INTENT(in)     :: nsoil_sb           ! number of soil layers as used/defined by the SB_ process
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)     :: soil_depth_acc_sl  !< accumulated depth array of soil profile [m]
    REAL(wp)                        , INTENT(in)     :: soil_p_labile, &   !< labile P from Yang et al. 2013, g/m2
                                                        soil_p_slow, &     !< secondary P from Yang et al. 2013, g/m2
                                                        soil_p_occluded, & !< Occluded P from Yang et al. 2013, g/m2
                                                        soil_p_primary     !< Primary from Yang et al. 2013, g/m2
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(inout)  :: po4_assoc_slow, &  !< slow minerally associated PO4 pool per soil layer [mol/m3]
                                                        po4_occluded, &    !< occluded PO4 pool per soil layer [mol/m3]
                                                        po4_primary, &     !< primary PO4 pool per soil layer [mol/m3]
                                                        labile_p           !< PO4 in solution per soil layer [mol / m3]
    CHARACTER(len=*)                , INTENT(in)     :: sitename           !< spp site name: BBR or LUE (MIT, VES, CON)
    ! ---------------------------
    ! 0.2 Local
    INTEGER                                      :: isoil, ichunk
    REAL(wp)                                     :: hlp1, hlp2, hlp3, hlp4, hlp5, hlp6
    REAL(wp),DIMENSION(nsoil_sb)                 :: arr_hlp5
    REAL(wp),DIMENSION(nc)                       :: arr_hlp2, arr_hlp3, arr_hlp4, soil_p_layer
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init_SPP'


    !> 1.0 set site-level soil properties
    !!
    SELECT CASE (sitename)
      CASE("BBR")
        arr_hlp5(:) = (/9.24_wp,36.92_wp,39.44_wp,40.36_wp,42.75_wp, &	
                       45.72_wp,32.45_wp,95.89_wp,124.14_wp,133.93_wp, &
                       100._wp,126._wp,122.18_wp,145._wp,90._wp/) ! SIP, mol P m-3
        hlp1 = soil_p_labile     !labile_P from Yang et al. 2013, g/m2
        hlp2 = soil_p_slow       !slow_P from Yang et al. 2013, g/m2
        hlp3 = soil_p_primary    !occluded_P from Yang et al. 2013, g/m2
        hlp4 = soil_p_occluded   !primary_P from Yang et al. 2013, g/m2
        soil_p_layer(:) = 6._wp 
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        hlp5        = 1750._wp / 1000._wp ! g P / kg soil
      CASE("LUE")
        arr_hlp5(:) = (/0.8867_wp,3.187_wp,2.662_wp,4.887_wp,4.792_wp, &	
                       2.987_wp,3.313_wp,6.0616_wp,7.56_wp,5.787_wp, &
                       5.932_wp,6.845_wp,8.294_wp,8.413_wp,8.532_wp/) !,1.676_wp
        hlp1 = soil_p_labile     !labile_P from Yang et al. 2013, g/m2
        hlp2 = soil_p_slow       !slow_P from Yang et al. 2013, g/m2
        hlp3 = soil_p_occluded   !occluded_P from Yang et al. 2013, g/m2
        hlp4 = soil_p_primary	   !primary_P from Yang et al. 2013, g/m2
        soil_p_layer(:) = 6._wp 
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        hlp5        = 120._wp / 1000._wp ! g P / kg soil
      CASE("VES")
        arr_hlp5(:) = (/2.388_wp,8.509_wp,9.033_wp,10.811_wp,15.365_wp, &	
                       20.661_wp,22.721_wp,29.321_wp,30.952_wp,32.821_wp, &
                       40.673_wp,41.880_wp,41.111_wp,48.438_wp,54._wp/) !,6.760_wp
        hlp1 = soil_p_labile     !labile_P from Yang et al. 2013, g/m2
        hlp2 = soil_p_slow       !slow_P from Yang et al. 2013, g/m2
        hlp3 = soil_p_primary    !occluded_P from Yang et al. 2013, g/m2
        hlp4 = soil_p_occluded   !primary_P from Yang et al. 2013, g/m2
        soil_p_layer(:) = 6._wp 
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        hlp5        = 900._wp / 1000._wp ! g P / kg soil
      CASE("MIT")
        arr_hlp5(:) = (/3.723_wp,6.272_wp,7.556_wp,8.292_wp,6.93_wp, &	
                       13.731_wp,13.296_wp,16.511_wp,12.890_wp,15._wp, &
                       16.667_wp,16.875_wp,17.5_wp,15.104_wp,15.625_wp/) !,6.760_wp
        hlp1 = soil_p_labile     !labile_P from Yang et al. 2013, g/m2
        hlp2 = soil_p_slow       !slow_P from Yang et al. 2013, g/m2
        hlp3 = soil_p_primary    !occluded_P from Yang et al. 2013, g/m2
        hlp4 = soil_p_occluded   !primary_P from Yang et al. 2013, g/m2
        soil_p_layer(:) = 9._wp 
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        hlp5        = 900._wp / 1000._wp ! g P / kg soil
        ! arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
      CASE("CON")
        arr_hlp5(:) = (/1.862_wp,2.924_wp,4.027_wp,4.048_wp,4.069_wp, &	
                       7.167_wp,10.855_wp,7.719_wp,18.75_wp,18.33_wp, &
                       23.44_wp,25._wp,21.093_wp,21.875_wp,15.104_wp/) !,6.760_wp
        hlp1 = soil_p_labile     !labile_P from Yang et al. 2013, g/m2
        hlp2 = soil_p_slow       !slow_P from Yang et al. 2013, g/m2
        hlp3 = soil_p_primary    !occluded_P from Yang et al. 2013, g/m2
        hlp4 = soil_p_occluded   !primary_P from Yang et al. 2013, g/m2
        soil_p_layer(:) = 9._wp 
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        hlp5        = 900._wp / 1000._wp ! g P / kg soil
        ! arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
    END SELECT

    !> 2.0 run calculations of soil-properties for all sites, i.e, all points of the chunk
    !! @TODO this is inconsistent with section 1.0 of this routine
    DO isoil = 1, nsoil_sb
      arr_hlp4(:) = MIN(2._wp * soil_depth_acc_sl(:,isoil) ** 0.5_wp, 3.4_wp) / 3.5_wp ! primary corrected for depth, [g P] 
      arr_hlp2(:) = 0.8_wp ** (2._wp * soil_depth_acc_sl(:,isoil)) / 3.5_wp ! labile + slow
      arr_hlp3(:) = MAX((1.0_wp - arr_hlp4(:) - arr_hlp2(:)), 0.2_wp / 3.5_wp) ! occluded corrected for depth

      WHERE ((arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:)) > 1._wp)
        arr_hlp4(:) = arr_hlp4(:) / (arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:))
        arr_hlp2(:) = arr_hlp2(:) / (arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:))
        arr_hlp3(:) = arr_hlp3(:) / (arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:))
      ENDWHERE
      
      WHERE (isoil < soil_p_layer(:)) 
        po4_assoc_slow(:,isoil) = arr_hlp5(isoil) * (hlp1 + hlp2) / (hlp1 + hlp2 + hlp3 + hlp4) * 8._wp / 17._wp
        po4_occluded(:,isoil)   = arr_hlp5(isoil) * hlp3 / (hlp1 + hlp2 + hlp3 + hlp4)
        po4_primary(:,isoil)    = arr_hlp5(isoil) * hlp4 / (hlp1 + hlp2 + hlp3 + hlp4)
        labile_P(:,isoil)       = arr_hlp5(isoil) * (hlp1 + hlp2) / (hlp1 + hlp2 + hlp3 + hlp4) * 9._wp / 17._wp !p_labile_sl in the follwing Eqn, [mol / m3]
      ELSEWHERE
        po4_assoc_slow(:,isoil) = arr_hlp5(isoil) * arr_hlp2(:) * 8._wp / 17._wp
        po4_occluded(:,isoil)   = arr_hlp5(isoil) * arr_hlp3(:)
        po4_primary(:,isoil)    = arr_hlp5(isoil) * arr_hlp4(:)
        labile_P(:,isoil)       = arr_hlp5(isoil) * arr_hlp2(:) * 9._wp / 17._wp !p_labile_sl in the follwing Eqn, [mol / m3]
      ENDWHERE
    ENDDO

  END SUBROUTINE sb_init_SPP

END MODULE mo_sb_init
