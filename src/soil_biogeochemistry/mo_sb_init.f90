!>
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
    USE mo_sb_jsm_processes                 ! e.g. calc_bulk_soil_carbon()
    USE mo_sse_util,                        ONLY: calc_qmax_texture
    USE mo_sb_jsm_transport                 ! e.g. calc_bioturbation_rate()
    USE mo_veg_constants,                   ONLY: itree
    USE mo_sb_constants                     ! e.g.  nh4_solute_prescribe, no3_solute_prescribe, po4_solute_prescribe

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_config(SSE_)

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
    INTEGER                                      :: isoil, ichunk
    INTEGER                                      :: iblk, ics, ice, nc
    REAL(wp)                                     :: hlp1,hlp2,hlp3,hlp4,hlp5,hlp6,frac_c13,frac_c14
    REAL(wp),DIMENSION(options%nc, options%nsoil_sb) :: fine_root_carbon_sl, labile_p, arr_hlp1, arr_z
    REAL(wp),DIMENSION(options%nsoil_sb)         :: arr_hlp5, arr_hlp6
    REAL(wp),DIMENSION(options%nc)               :: soil_p_layer, arr_hlp4, arr_hlp2, arr_hlp3
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init'
    ! ---------------------------
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_config(SSE_)
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
    dsl4jsb_Real3D_onChunk :: qmax_po4_om_sl
    dsl4jsb_Real3D_onChunk :: silt_sl
    dsl4jsb_Real3D_onChunk :: clay_sl
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk :: nh4_solute
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute
    dsl4jsb_Real3D_onChunk :: no3_n15_solute
    dsl4jsb_Real3D_onChunk :: bulk_soil_carbon_sl
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onChunk :: qmax_org
    dsl4jsb_Real3D_onChunk :: qmax_po4
    dsl4jsb_Real3D_onChunk :: km_adsorpt_po4_sl
    dsl4jsb_Real3D_onChunk :: eca_divisor_adsorp
    dsl4jsb_Real3D_onChunk :: k_bioturb
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast
    dsl4jsb_Real3D_onChunk :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: po4_occluded
    dsl4jsb_Real3D_onChunk :: ph_sl
    dsl4jsb_Real3D_onChunk :: lithology_class_sl
    dsl4jsb_Real3D_onChunk :: po4_primary
    dsl4jsb_Real3D_onChunk :: microbial_cue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk :: vmax_weath_mineral_sl
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly
    dsl4jsb_Real3D_onChunk :: enzyme_frac_residue
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_c
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_n
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_p
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
    ! Get process config
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_config(SSE_)
    ! Get process memories
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(VEG_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var3D_onChunk(VEG_,      root_fraction_sl)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SSE_,     bulk_dens_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     w_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     qmax_org_min_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     qmax_po4_min_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     qmax_po4_om_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     silt_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_,     clay_sl)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_soil_carbon_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_org)
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       km_adsorpt_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       eca_divisor_adsorp)
    dsl4jsb_Get_var3D_onChunk(SB_,       k_bioturb)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_fast)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_slow)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_occluded) 
    dsl4jsb_Get_var3D_onChunk(SB_,       ph_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       lithology_class_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_primary)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_cue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,       vmax_weath_mineral_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly)
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_residue)
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly_c)
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly_n)
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly_p)
    ! Get process memory pool structures
    dsl4slm_get_vegpool_var(pool,fine_root,carbon)
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


    !> 2.0 differentiate soil model init
    !! 
    hlp6 = 0.6_wp   ! fraction of production not roots 

    SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))

       ! Initial values for Soenke's four pool soil biogeochemical model
       CASE("simple_0d", "simple_1d")

          hlp1 = 100._wp ! g C m-2 metabolic litter
          hlp2 = 100._wp ! g C m-2 structural litter
          hlp3 = 1000._wp ! g C m-2 woody litter
          hlp4 = 1500._wp ! g C m-2 fast SOM
          hlp5 = 10000._wp ! g C m-2 slow SOM

          IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_0d") THEN
             sb_pool%soluable_litter%carbon(:,1) = hlp1 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(1)
             sb_pool%polymeric_litter%carbon(:,1) = hlp2 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(1)
             IF(lctlib%growthform==itree) &
                sb_pool%woody_litter%carbon(:,1) = hlp3 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(1)
             sb_pool%microbial%carbon(:,1) = hlp4 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(1)
             sb_pool%residue%carbon(:,1) = hlp5 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(1)

             nh4_solute(:,1) = nh4_solute_prescribe / dsl4jsb_Config(SSE_)%dz_water(1)
             no3_solute(:,1) = no3_solute_prescribe / dsl4jsb_Config(SSE_)%dz_water(1)


             ! inital Qmax based on current soil carbon content
             bulk_soil_carbon_sl(:,:) = calc_bulk_soil_carbon(nc, nsoil_sb, sb_pool)
             bulk_dens_corr_sl(:,:)   = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:), bulk_dens_sl(:,:))
             qmax_po4_min_sl(:,:)     = calc_qmax_texture(qmax_po4_mineral, 0.0_wp, 1._wp)
             qmax_po4_om_sl(:,:)      = qmax_po4_min_sl(:,:)/100._wp
             qmax_po4(:,:)            = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:),bulk_dens_sl(:,:), & 
                                        qmax_po4_om_sl(:,:),qmax_po4_min_sl(:,:))
             qmax_po4(:,:)            = qmax_po4(:,:) * w_soil_sl(:,:) / soil_depth_sl(:,:)
             km_adsorpt_po4_sl(:,:) = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:),bulk_dens_sl(:,:), &
                                      km_adsorpt_OM_po4,km_adsorpt_mineral_po4)
             ! initialise mineral soil P pools
             CALL init_soil_pi_pools(nc, nsoil_sb, &
                                    TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &
                                    soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), &
                                    dsl4jsb_Config(SB_)%soil_p_labile, dsl4jsb_Config(SB_)%soil_p_slow, &
                                    dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
                                    dsl4jsb_Config(SB_)%soil_p_depth, dsl4jsb_Config(SB_)%flag_sb_CLE_po4, & 
                                    po4_solute(:,:), po4_assoc_fast(:,:), po4_assoc_slow(:,:), &
                                    po4_occluded(:,:), po4_primary(:,:), &
                                    qmax_po4(:,:),km_adsorpt_po4_sl(:,:))
          ELSE 
             sb_pool%soluable_litter%carbon(:,1)  = hlp1 /  molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(1) * &
                                                                       ( hlp6 + (1.0_wp - hlp6) * root_fraction_sl(:,1)) 
             sb_pool%polymeric_litter%carbon(:,1) = hlp2 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(1) * &
                                                                       ( hlp6 + (1.0_wp - hlp6) * root_fraction_sl(:,1)) 
             IF(lctlib%growthform==itree) &
                sb_pool%woody_litter%carbon(:,1)  = hlp3 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(1)
             sb_pool%microbial%carbon(:,1)        = hlp4 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(1) * &
                                                                       ( hlp6 + (1.0_wp - hlp6) * root_fraction_sl(:,1)) 
             sb_pool%residue%carbon(:,1)          = hlp5 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(1) * & 
                                                                       ( hlp6 + (1.0_wp - hlp6) * root_fraction_sl(:,1)) 
             DO isoil=2,nsoil_sb
                sb_pool%soluable_litter%carbon(:,isoil)  = hlp1 /  molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(isoil) * &
                                                                          (1.0_wp - hlp6) * root_fraction_sl(:,isoil) 
                sb_pool%polymeric_litter%carbon(:,isoil) = hlp2 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(isoil) * &
                                                                          (1.0_wp - hlp6) * root_fraction_sl(:,isoil) 
                sb_pool%microbial%carbon(:,isoil)        = hlp4 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(isoil) * &
                                                                          (1.0_wp - hlp6) * root_fraction_sl(:,isoil) 
                sb_pool%residue%carbon(:,isoil)          = hlp5 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(isoil) * & 
                                                                          (1.0_wp - hlp6) * root_fraction_sl(:,isoil)

             ENDDO
! note: this differs on purpose from the calculation in mo_sb_update_pools for the prescribed option
!       this is because here I'm trying to enure that the models are initialised with the same amount of nutrient
!       whereas in the prescribed option the target is unlimiting concentrations across the profile. The latter leads
!       to more nutrient in the soil column in 1d compared to 0d
             nh4_solute(:,:) = nh4_solute_prescribe / SUM(dsl4jsb_Config(SSE_)%dz_water(:))
             no3_solute(:,:) = no3_solute_prescribe / SUM(dsl4jsb_Config(SSE_)%dz_water(:))

             ! inital bulk density and qmax_po4 based on current soil carbon content
             bulk_soil_carbon_sl(:,:) = calc_bulk_soil_carbon(nc, nsoil_sb, sb_pool)
             bulk_dens_corr_sl(:,:) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:), bulk_dens_sl(:,:))            
             qmax_po4_min_sl(:,:) = calc_qmax_texture(qmax_po4_mineral, 0.0_wp, 1._wp)
             qmax_po4_om_sl(:,:) = qmax_po4_min_sl(:,:)/100._wp
             qmax_po4(:,:) = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:),bulk_dens_sl(:,:), & 
                                      qmax_po4_om_sl(:,:),qmax_po4_min_sl(:,:))
             qmax_po4(:,:)            = qmax_po4(:,:) * w_soil_sl(:,:) / soil_depth_sl(:,:)
             km_adsorpt_po4_sl(:,:) = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:),bulk_dens_sl(:,:), &
                                      km_adsorpt_OM_po4,km_adsorpt_mineral_po4)
             ! initialise mineral soil P pools
             CALL init_soil_pi_pools(nc, nsoil_sb, &
                                  TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &
                                  soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), &
                                  dsl4jsb_Config(SB_)%soil_p_labile, dsl4jsb_Config(SB_)%soil_p_slow, &
                                  dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
                                  dsl4jsb_Config(SB_)%soil_p_depth, dsl4jsb_Config(SB_)%flag_sb_CLE_po4, & 
                                  po4_solute(:,:), po4_assoc_fast(:,:), po4_assoc_slow(:,:), &
                                  po4_occluded(:,:), po4_primary(:,:), &
                                  qmax_po4(:,:),km_adsorpt_po4_sl(:,:))
          ENDIF
 
          ! metabolic litter all other elements
          sb_pool%soluable_litter%nitrogen(:,:)   = sb_pool%soluable_litter%carbon(:,:) / & 
                                                                    (25._wp * molar_mass_N / molar_mass_C)
          sb_pool%soluable_litter%phosphorus(:,:) = sb_pool%soluable_litter%nitrogen(:,:) / & 
                                                                    (14._wp * molar_mass_P / molar_mass_N) 
          sb_pool%soluable_litter%carbon13(:,:)   = sb_pool%soluable_litter%carbon(:,:) * frac_c13
          sb_pool%soluable_litter%carbon14(:,:)   = sb_pool%soluable_litter%carbon(:,:) * frac_c14 
          sb_pool%soluable_litter%nitrogen15(:,:) = sb_pool%soluable_litter%nitrogen(:,:) / & 
                                                                    (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! non-woody litter all other elements
          sb_pool%polymeric_litter%nitrogen(:,:)   = sb_pool%polymeric_litter%carbon(:,:) / & 
                                                                    (100._wp * molar_mass_N / molar_mass_C)
          sb_pool%polymeric_litter%phosphorus(:,:) = sb_pool%polymeric_litter%nitrogen(:,:) / & 
                                                                    (14._wp * molar_mass_P / molar_mass_N) 
          sb_pool%polymeric_litter%carbon13(:,:)   = sb_pool%polymeric_litter%carbon(:,:) * frac_c13
          sb_pool%polymeric_litter%carbon14(:,:)   = sb_pool%polymeric_litter%carbon(:,:) * frac_c14
          sb_pool%polymeric_litter%nitrogen15(:,:) = sb_pool%polymeric_litter%nitrogen(:,:) / & 
                                                                    (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! woody litter all other elements
          sb_pool%woody_litter%nitrogen(:,:)   = sb_pool%woody_litter%carbon(:,:) / & 
                                                                (260._wp * molar_mass_N / molar_mass_C)
          sb_pool%woody_litter%phosphorus(:,:) = sb_pool%woody_litter%nitrogen(:,:) / & 
                                                                (14._wp * molar_mass_P / molar_mass_N) 
          sb_pool%woody_litter%carbon13(:,:)   = sb_pool%woody_litter%carbon(:,:) * frac_c13
          sb_pool%woody_litter%carbon14(:,:)   = sb_pool%woody_litter%carbon(:,:) * frac_c14
          sb_pool%woody_litter%nitrogen15(:,:) = sb_pool%woody_litter%nitrogen(:,:) / & 
                                                                (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          ! fast pool (=microbial) all other elements
          sb_pool%microbial%nitrogen(:,:)   = sb_pool%microbial%carbon(:,:) / SOM_cn
          sb_pool%microbial%phosphorus(:,:) = sb_pool%microbial%nitrogen(:,:) / SOM_np 
          sb_pool%microbial%carbon13(:,:)   = sb_pool%microbial%carbon(:,:) * frac_c13
          sb_pool%microbial%carbon14(:,:)   = sb_pool%microbial%carbon(:,:) * frac_c14
          sb_pool%microbial%nitrogen15(:,:) = sb_pool%microbial%nitrogen(:,:) / & 
                                                             (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          sb_pool%residue%nitrogen(:,:)   = sb_pool%residue%carbon(:,:) / SOM_cn
          sb_pool%residue%phosphorus(:,:) = sb_pool%residue%nitrogen(:,:) / SOM_np 
          sb_pool%residue%carbon13(:,:)   = sb_pool%residue%carbon(:,:) * frac_c13
          sb_pool%residue%carbon14(:,:)   = sb_pool%residue%carbon(:,:) * frac_c14 
          sb_pool%residue%nitrogen15(:,:) = sb_pool%residue%nitrogen(:,:) / & 
                                                             (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          nh4_n15_solute(:,:) = nh4_solute(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
          no3_n15_solute(:,:) = no3_solute(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 

       CASE('jsm')

          hlp1 = 100._wp  ! g C m-2 metabolic litter
          hlp2 = 100._wp  ! g C m-2 structural litter
          hlp3 = 1000._wp ! g C m-2 woody litter
          hlp4 = 30._wp   ! g C m-2 microbial biomass
          hlp5 = 150._wp  ! g C m-2 microbial necromass

          sb_pool%soluable_litter%carbon(:,1)  = hlp1 /  molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(1) * &
                                                                    ( hlp6 + (1.0_wp - hlp6) * root_fraction_sl(:,1)) 
          sb_pool%polymeric_litter%carbon(:,1) = hlp2 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(1) * &
                                                                    ( hlp6 + (1.0_wp - hlp6) * root_fraction_sl(:,1)) 
          IF(lctlib%growthform==itree) &
             sb_pool%woody_litter%carbon(:,1)  = hlp3 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(1)
          DO isoil=2,nsoil_sb
             sb_pool%soluable_litter%carbon(:,isoil)  = hlp1 /  molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(isoil) * &
                                                                       (1.0_wp - hlp6) * root_fraction_sl(:,isoil) 
             sb_pool%polymeric_litter%carbon(:,isoil) = hlp2 / molar_mass_C / dsl4jsb_Config(SSE_)%dz_water(isoil) * &
                                                                       (1.0_wp - hlp6) * root_fraction_sl(:,isoil) 
          ENDDO

          DO isoil=1,nsoil_sb
             sb_pool%microbial%carbon(:,isoil) = MAX(hlp4/dsl4jsb_Config(SSE_)%dz_water(isoil) * root_fraction_sl(:,isoil), &
                                                     min_mic_biomass)
             sb_pool%dom%carbon(:,isoil) = sb_pool%microbial%carbon(:,isoil) * 0.1_wp
             sb_pool%residue%carbon(:,isoil) = hlp5/dsl4jsb_Config(SSE_)%dz_water(isoil) * root_fraction_sl(:,isoil)
             hlp6 = 0.3_wp !MIN(0.3_wp, isoil*0.1_wp)
             sb_pool%residue_assoc%carbon(:,isoil) = dsl4jsb_Config(SB_)%qmax_org_fine_particle * (clay_sl(:,isoil) + silt_sl(:,isoil)) * & 
                                                     bulk_dens_sl(:,isoil) * hlp6 / &
                                                     dsl4jsb_Config(SSE_)%dz_water(isoil) * root_fraction_sl(:,isoil)
          ENDDO
          sb_pool%soluable_litter%nitrogen(:,:)   = sb_pool%soluable_litter%carbon(:,:) / microbial_cn 
          sb_pool%soluable_litter%phosphorus(:,:) = sb_pool%soluable_litter%nitrogen(:,:) / SOM_np
          sb_pool%soluable_litter%carbon13(:,:)   = sb_pool%soluable_litter%carbon(:,:) * frac_c13 
          sb_pool%soluable_litter%carbon14(:,:)   = sb_pool%soluable_litter%carbon(:,:) * frac_c14 
          sb_pool%soluable_litter%nitrogen15(:,:) = sb_pool%soluable_litter%nitrogen(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          sb_pool%polymeric_litter%nitrogen(:,:)   = sb_pool%polymeric_litter%carbon(:,:) / 50._wp 
          sb_pool%polymeric_litter%phosphorus(:,:) = sb_pool%polymeric_litter%nitrogen(:,:) / SOM_np
          sb_pool%polymeric_litter%carbon13(:,:)   = sb_pool%polymeric_litter%carbon(:,:) * frac_c13 
          sb_pool%polymeric_litter%carbon14(:,:)   = sb_pool%polymeric_litter%carbon(:,:) * frac_c14 
          sb_pool%polymeric_litter%nitrogen15(:,:) = sb_pool%polymeric_litter%nitrogen(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          sb_pool%woody_litter%nitrogen(:,:)   = sb_pool%woody_litter%carbon(:,:) / 200._wp 
          sb_pool%woody_litter%phosphorus(:,:) = sb_pool%woody_litter%nitrogen(:,:) / SOM_np
          sb_pool%woody_litter%carbon13(:,:)   = sb_pool%woody_litter%carbon(:,:) * frac_c13 
          sb_pool%woody_litter%carbon14(:,:)   = sb_pool%woody_litter%carbon(:,:) * frac_c14 
          sb_pool%woody_litter%nitrogen15(:,:) = sb_pool%woody_litter%nitrogen(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          sb_pool%microbial%nitrogen(:,:)   = sb_pool%microbial%carbon(:,:) / microbial_cn 
          sb_pool%microbial%phosphorus(:,:) = sb_pool%microbial%nitrogen(:,:) / microbial_np
          sb_pool%microbial%carbon13(:,:)   = sb_pool%microbial%carbon(:,:) * frac_c13 
          sb_pool%microbial%carbon14(:,:)   = sb_pool%microbial%carbon(:,:) * frac_c14 
          sb_pool%microbial%nitrogen15(:,:) = sb_pool%microbial%nitrogen(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          sb_pool%residue%nitrogen(:,:)   = sb_pool%residue%carbon(:,:) / microbial_cn 
          sb_pool%residue%phosphorus(:,:) = sb_pool%residue%nitrogen(:,:) / microbial_np
          sb_pool%residue%carbon13(:,:)   = sb_pool%residue%carbon(:,:) * frac_c13 
          sb_pool%residue%carbon14(:,:)   = sb_pool%residue%carbon(:,:) * frac_c14 
          sb_pool%residue%nitrogen15(:,:) = sb_pool%residue%nitrogen(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          sb_pool%residue_assoc%nitrogen(:,:)   = sb_pool%residue_assoc%carbon(:,:) / SOM_cn 
          sb_pool%residue_assoc%phosphorus(:,:) = sb_pool%residue_assoc%nitrogen(:,:) / SOM_np
          sb_pool%residue_assoc%carbon13(:,:)   = sb_pool%residue_assoc%carbon(:,:) * frac_c13 
          sb_pool%residue_assoc%carbon14(:,:)   = sb_pool%residue_assoc%carbon(:,:) * frac_c14 
          sb_pool%residue_assoc%nitrogen15(:,:) = sb_pool%residue_assoc%nitrogen(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
     
         ! some initial values for testing
         nh4_solute(:,:) = nh4_solute_prescribe / SUM(dsl4jsb_Config(SSE_)%dz_water(:))
         nh4_n15_solute(:,:) = nh4_solute(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
         no3_solute(:,:) = no3_solute_prescribe / SUM(dsl4jsb_Config(SSE_)%dz_water(:))
         no3_n15_solute(:,:) = no3_solute(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
 
         ! inital bulk density, Qmax_org & Qmax_po4 based on current soil carbon content
         bulk_soil_carbon_sl(:,:) = calc_bulk_soil_carbon(nc, nsoil_sb, sb_pool)
         bulk_dens_corr_sl(:,:) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:), bulk_dens_sl(:,:))
         qmax_org_min_sl(:,:) = calc_qmax_texture(dsl4jsb_Config(SB_)%qmax_org_fine_particle, silt_sl(:,:), clay_sl(:,:))
         qmax_po4_min_sl(:,:) = calc_qmax_texture(qmax_po4_mineral, 0.0_wp, 1._wp)
         qmax_po4_om_sl(:,:) = qmax_po4_min_sl(:,:)/10._wp
         qmax_org(:,:) = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:),bulk_dens_sl(:,:), &
                         0._wp,qmax_org_min_sl(:,:))
         qmax_po4(:,:) = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:),bulk_dens_sl(:,:), &
                         qmax_po4_om_sl(:,:),qmax_po4_min_sl(:,:))
         qmax_po4(:,:) = qmax_po4(:,:) * w_soil_sl(:,:) / soil_depth_sl(:,:)
         km_adsorpt_po4_sl(:,:) = calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:),bulk_dens_sl(:,:), &
                         km_adsorpt_OM_po4,km_adsorpt_mineral_po4)

        ! initialise mineral soil P pools
        DO isoil = 1,nsoil_sb
          WHERE(soil_depth_sl(:,isoil) > eps8)
            fine_root_carbon_sl(:,isoil) = root_fraction_sl(:,isoil) / soil_depth_sl(:,isoil) * & 
                                        veg_pool%fine_root%carbon(:)
          ENDWHERE
          IF (isoil == 1) THEN
            arr_z(:,isoil) = soil_depth_sl(:,isoil) / 2._wp
          ELSE
            arr_z(:,isoil) = arr_z(:,isoil-1) + (soil_depth_sl(:,isoil) + soil_depth_sl(:,isoil-1)) / 2._wp  
          ENDIF   
        ENDDO
        arr_hlp1(:,:) = km_adsorpt_po4_sl + qmax_po4 + &
                        km_adsorpt_po4_sl / km_plant_uptake_po4 * f_plant_carrier_po4 * fine_root_carbon_sl + &
                        km_adsorpt_po4_sl / km_mic_uptake_po4 * f_mic_carrier_po4 * sb_pool%microbial%carbon(:,:)

        CALL init_soil_pi_pools(nc, nsoil_sb, &
                                TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &
                                soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), &
                                dsl4jsb_Config(SB_)%soil_p_labile, dsl4jsb_Config(SB_)%soil_p_slow, &
                                dsl4jsb_Config(SB_)%soil_p_occluded, dsl4jsb_Config(SB_)%soil_p_primary, &
                                dsl4jsb_Config(SB_)%soil_p_depth, dsl4jsb_Config(SB_)%flag_sb_CLE_po4, & 
                                po4_solute(:,:), po4_assoc_fast(:,:), po4_assoc_slow(:,:), &
                                po4_occluded(:,:), po4_primary(:,:), &
                                qmax_po4(:,:),km_adsorpt_po4_sl(:,:),arr_hlp1(:,:),eca_divisor_adsorp(:,:))

        ! initialise diffusion by bioturbation
        k_bioturb(:,:) = calc_bioturbation_rate(nc, nsoil_sb, bulk_dens_corr_sl(:,:), root_fraction_sl(:,:), soil_depth_sl(:,:)) 
        ! enzyme allocation initialise 
        enzyme_frac_poly(:,:)     = 0.5_wp 
        enzyme_frac_residue(:,:)  = 0.5_wp 
        enzyme_frac_poly_c(:,:)   = 0.5_wp 
        enzyme_frac_poly_n(:,:)   = 0.5_wp 
        enzyme_frac_poly_p(:,:)   = 0.5_wp 

        ! isotopic signals in all pools
        ! microbial, DOM, aDOM, residue, and aRes pools
        sb_pool%dom%nitrogen(:,:)     = sb_pool%dom%carbon(:,:) / microbial_cn
        sb_pool%dom%phosphorus(:,:)   = sb_pool%dom%nitrogen(:,:) / microbial_np
        sb_pool%dom%carbon13(:,:)     = sb_pool%dom%carbon(:,:) * frac_c13
        sb_pool%dom%carbon14(:,:)     = sb_pool%dom%carbon(:,:) * frac_c14
        sb_pool%dom%nitrogen15(:,:)   = sb_pool%dom%nitrogen(:,:) / (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

        sb_pool%dom_assoc%nitrogen(:,:)     = sb_pool%dom_assoc%carbon(:,:) / SOM_cn
        sb_pool%dom_assoc%phosphorus(:,:)   = sb_pool%dom_assoc%nitrogen(:,:) / SOM_np
        sb_pool%dom_assoc%carbon13(:,:)     = sb_pool%dom_assoc%carbon(:,:) * frac_c13
        sb_pool%dom_assoc%carbon14(:,:)     = sb_pool%dom_assoc%carbon(:,:) * frac_c14
        sb_pool%dom_assoc%nitrogen15(:,:)   = sb_pool%dom_assoc%nitrogen(:,:) / & 
                                                                  (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

        ! sb_pool%residue_assoc%carbon13(:,:)   = sb_pool%residue_assoc%carbon(:,:) * frac_c13
        ! sb_pool%residue_assoc%carbon14(:,:)   = sb_pool%residue_assoc%carbon(:,:) * frac_c14 
        ! sb_pool%residue_assoc%nitrogen15(:,:) = sb_pool%residue_assoc%nitrogen(:,:) / & 
                                                                  ! (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        ! fungi and mycorrhiza
        sb_pool%fungi%carbon13(:,:)   = sb_pool%fungi%carbon(:,:) * frac_c13
        sb_pool%fungi%carbon14(:,:)   = sb_pool%fungi%carbon(:,:) * frac_c14 
        sb_pool%fungi%nitrogen15(:,:) = sb_pool%fungi%nitrogen(:,:) / & 
                                                                  (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        sb_pool%mycorrhiza%carbon13(:,:)   = sb_pool%mycorrhiza%carbon(:,:) * frac_c13
        sb_pool%mycorrhiza%carbon14(:,:)   = sb_pool%mycorrhiza%carbon(:,:) * frac_c14 
        sb_pool%mycorrhiza%nitrogen15(:,:) = sb_pool%mycorrhiza%nitrogen(:,:) / & 
                                                                  (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
       CASE DEFAULT
         CALL finish(routine,'Invalid setting of sb_model_scheme.')
 
    END SELECT


    !> 3.0 other SB_ variables - all soil models
    !! 
    microbial_cue_eff_tmic_mavg(:,:)  = microbial_cue_max
    vmax_weath_mineral_sl(:,:)        = k_weath_mineral                       ! (parameter from mo_sb_constants)

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
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_soil_pi_pools(nc, nsoil_sb, sb_model_scheme, sb_adsorp_scheme, soil_depth_sl, bulk_dens_sl, &
                               soil_p_labile, soil_p_slow, soil_p_occluded, soil_p_primary, soil_p_depth, flag_sb_CLE_po4, &
                               po4_solute, po4_assoc_fast, po4_assoc_slow, po4_occluded, po4_primary, &
                               qmax_po4,km, eca_enzyme_hlp, eca_divisor_adsorp)

    USE mo_jsb_phys_constants,    ONLY: molar_mass_P
    USE mo_sb_constants
    USE mo_sb_jsm_processes,      ONLY: calc_fast_po4

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)     :: nc                 !< dimensions
    INTEGER                                         , INTENT(in)     :: nsoil_sb           ! number of soil layers as used/defined by the SB_ process
    CHARACTER(len=*)                                , INTENT(in)     :: sb_model_scheme    !< soil model scheme 
    CHARACTER(len=*)                                , INTENT(in)     :: sb_adsorp_scheme   !< soil model po4 adsorption scheme 
    LOGICAL                                         , INTENT(in)     :: flag_sb_CLE_po4    !< TRUE: equilibrium is affected by plant and microbe;
    REAL(wp), DIMENSION(nc,nsoil_sb)                , INTENT(in)     :: soil_depth_sl, &   !< soil depth per soil layer, [m]
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
    REAL(wp), DIMENSION(nc,nsoil_sb),         INTENT(in),  OPTIONAL  :: eca_enzyme_hlp     !< helper item for enzyme abundance in eca approach [mol P / m3]
    REAL(wp), DIMENSION(nc,nsoil_sb),         INTENT(out), OPTIONAL  :: eca_divisor_adsorp !< denominator of ECA item for po4 adsorption [mol P / m3]
    ! ---------------------------
    ! 0.2 Local
    INTEGER                                          :: i
    REAL(wp)                                         :: rhlp1, rhlp2
    LOGICAL,  DIMENSION(nc)                          :: aflag
    REAL(wp), DIMENSION(nc)                          :: soil_p_weight, soil_weight_0d, &  !< total soil weight up to the soil_p_depth/soil_depth_0d, [kg]
                                                        soil_p_layer                      !< up to soil_p_layer keep the P init map info
    REAL(wp), DIMENSION(nc)                          :: hlp1, hlp2, hlp3, hlp4
    REAL(wp), DIMENSION(nc)                          :: soil_depth_0d
    REAL(wp), DIMENSION(nc,nsoil_sb)                 :: soil_weight_sl, &                 !< soil weight per soil layer
                                                        arr_hlp1, arr_hlp2, arr_z
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_soil_pi_pools'




    !> 0.9
    !! init local variables
    soil_depth_0d(:) = 1._wp ! soil depth in simple_0d

    !> 0.9.1
    !! initialise output
    po4_solute(:,:)     = 0._wp
    po4_assoc_fast(:,:) = 0._wp
    po4_assoc_slow(:,:) = 0._wp 
    po4_occluded(:,:)   = 0._wp
    po4_primary(:,:)    = 0._wp


    !> 1.0 calculate the total soil weight up to the soil_p_depth,
    !! and the soil weight per soil layer
    hlp1(:) = 0._wp; hlp2(:) = 0._wp; aflag(:) = .TRUE.
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
    
    !> 2.0 simple_0d should represent the first meter, 
    !! simple_1d and jsm use the soil profile depth
    IF (sb_model_scheme=='simple_0d') THEN
      rhlp1 = 9._wp / 17._wp * (soil_p_slow + soil_p_labile)  ! temp correction for global average
      rhlp2 = 8._wp / 17._wp * (soil_p_slow + soil_p_labile)  ! to get equilibrium in fast and slow pool
      hlp1(:) = 0._wp; hlp2(:) = 0._wp; aflag(:) = .TRUE.
      DO i=1,nsoil_sb
        WHERE (aflag) 
          hlp1(:) = hlp1(:) + soil_depth_sl(:,i)
          hlp2(:) = hlp2(:) + soil_weight_sl(:,i)
          WHERE (hlp1(:) > soil_depth_0d(:))
            soil_weight_0d(:) = hlp2(:) - soil_weight_sl(:,i) * (hlp1(:) - soil_depth_0d(:)) / soil_depth_sl(:,i)
            aflag(:) = .FALSE.
          ENDWHERE
        ENDWHERE
      END DO
      po4_assoc_slow(:,1) = rhlp2 * soil_p_weight(:) / soil_weight_0d(:) / molar_mass_P
      po4_occluded(:,1)   = soil_p_occluded * soil_p_weight(:) / soil_weight_0d(:) / molar_mass_P
      po4_primary(:,1)    = soil_p_primary * soil_p_weight(:) / soil_weight_0d(:) / molar_mass_P
      arr_hlp1(:,1)       = rhlp1 * soil_p_weight(:) / soil_weight_0d(:) / molar_mass_P
      hlp1(:)             = soil_weight_0d(:) * qmax_po4_mineral !Smax in Eqn, [mol / m3]
      hlp2(:)             = soil_weight_0d(:) * km_adsorpt_mineral_po4 !K in Eqn, [mol / m3]
      po4_assoc_fast(:,1) = calc_fast_po4(arr_hlp1(:,1),hlp1(:),hlp2(:))
      po4_solute(:,1)     = arr_hlp1(:,1) - po4_assoc_fast(:,1)
    ELSE
      !> 2.1 calculate the initial slow P, occluded P, primary P pool, and labile_P
      !! based on the assumption that soil P pool are linearly corrected with soil weight
      DO i=1,nsoil_sb
        IF (i==1) THEN
          arr_z(:,i)        = soil_depth_sl(:,i) / 2._wp
        ELSE
          arr_z(:,i)          = arr_z(:,i-1) + (soil_depth_sl(:,i) + soil_depth_sl(:,i-1)) / 2._wp  
        ENDIF   
      ENDDO

      DO i=1,nsoil_sb
        hlp1(:)             = MIN(2._wp * arr_z(:,i) ** 0.5_wp, 3.4_wp) / 3.5_wp ! primary corrected for depth, [g P] 
        ! hlp2(:)             = MIN(3._wp, arr_z(:,i) ** (-0.75_wp)) / 3.5_wp ! labile + slow
        hlp2(:)             = 0.8_wp ** (2._wp * arr_z(:,i)) / 3.5_wp
        hlp3(:)             = MAX((1 - hlp1(:) - hlp2(:)), 0.2_wp / 3.5_wp) ! occluded corrected for depth

        WHERE ((hlp1(:) + hlp2(:) + hlp3(:)) > 1._wp)
          hlp1(:) = hlp1(:) / (hlp1(:) + hlp2(:) + hlp3(:))
          hlp2(:) = hlp2(:) / (hlp1(:) + hlp2(:) + hlp3(:))
          hlp3(:) = hlp3(:) / (hlp1(:) + hlp2(:) + hlp3(:))
        ENDWHERE
        WHERE (soil_p_layer(:) > i) 
          po4_assoc_slow(:,i) = bulk_dens_sl(:,i) * hlp4(:) * (soil_p_slow + soil_p_labile) / rhlp1 / molar_mass_P * 8._wp / 17._wp
          po4_occluded(:,i)   = bulk_dens_sl(:,i) * hlp4(:) * soil_p_occluded               / rhlp1 / molar_mass_P
          po4_primary(:,i)    = bulk_dens_sl(:,i) * hlp4(:) * soil_p_primary                / rhlp1 / molar_mass_P
          arr_hlp1(:,i)       = bulk_dens_sl(:,i) * hlp4(:) * (soil_p_slow + soil_p_labile) / rhlp1 / molar_mass_P * 9._wp / 17._wp !p_labile_sl in the follwing Eqn, [mol / m3]
        ELSEWHERE
          po4_assoc_slow(:,i) = bulk_dens_sl(:,i) * hlp4(:) * hlp2(:) / molar_mass_P * 8._wp / 17._wp
          po4_occluded(:,i)   = bulk_dens_sl(:,i) * hlp4(:) * hlp3(:) / molar_mass_P
          po4_primary(:,i)    = bulk_dens_sl(:,i) * hlp4(:) * hlp1(:) / molar_mass_P
          arr_hlp1(:,i)       = bulk_dens_sl(:,i) * hlp4(:) * hlp2(:) / molar_mass_P * 9._wp / 17._wp !p_labile_sl in the follwing Eqn, [mol / m3]
        ENDWHERE
      END DO

      !> 2.2 calculate the initial solute P and fast P,
      !! based on the assumption that labile P is the sum of both,
      !! and they follow the Langmuir equilibrium
      !! po4_assoc_fast = (p_labile_sl + Smax + K - sqrt((p_labile_sl + Smax + K)**2 - 4*p_labile_sl*Smax)) / 2
    
      DO i=1,nsoil_sb
        IF (sb_model_scheme == 'jsm' .AND. sb_adsorp_scheme /= 'eca_none' ) THEN
            po4_assoc_fast(:,i)        = calc_fast_po4(arr_hlp1(:,i),qmax_po4(:,i),km(:,i), eca_enzyme_hlp(:,i))
            po4_solute(:,i)            = arr_hlp1(:,i) - po4_assoc_fast(:,i)
            eca_divisor_adsorp(:,i)    = arr_hlp1(:,i) - po4_assoc_fast(:,i) + eca_enzyme_hlp(:,i)
        ELSE
          po4_assoc_fast(:,i) = calc_fast_po4(arr_hlp1(:,i),qmax_po4(:,i),km(:,i))
          po4_solute(:,i)     = arr_hlp1(:,i) - po4_assoc_fast(:,i)
        ENDIF

      END DO
    END IF

  END SUBROUTINE init_soil_pi_pools

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from mo_sb_init 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates the fast_assoc_po4 based on Langmuir equilibrium and Labile P pool
  !!
  !! Input: Labile P pool, Qmax, Km, all in the unit [mol m-3]
  !!
  !! Output: fast_assoc_po4
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_fast_po4(Labile_p, Qmax, km) RESULT (fast_po4)
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: Labile_p    !< labile P pool of the soil layer [mol P / m3]
    REAL(wp), INTENT(in) :: Qmax        !< maximium sorption capacity of po4 [mol P/ m3]
    REAL(wp), INTENT(in) :: km          !< sorption affinity of po4 [mol P/ m3]
    REAL(wp)             :: fast_po4    !< fast associated po4 of the soil layer [mol P / m3]
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fast_po4'

    !> based on the assumption that labile P is the sum of both,
    !! and they follow the Langmuir equilibrium
    !! po4_assoc_fast = (p_labile_sl + Smax + K - sqrt((p_labile_sl + Smax + K)**2 - 4*p_labile_sl*Smax)) / 2
    fast_po4 = (Labile_p + Qmax + km- sqrt((Labile_p + Qmax + km)**2.0_wp - 4.0_wp*Labile_p*Qmax)) / 2._wp

  END FUNCTION calc_fast_po4
END MODULE mo_sb_init
