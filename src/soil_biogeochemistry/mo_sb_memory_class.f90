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
!! @brief Contains the memory class for the soil_biogeochemistry processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20170103, based on the JSBACH4, MPI-M Hamburg)
!!
MODULE mo_sb_memory_class

  USE mo_sb_memory_class_util   ! pools code for quincy & jsb4-quincy 

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message
#ifdef __QUINCY_NO_JSB__
#else
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_utils_iface,        ONLY: assign_if_present_allocatable
#endif
  USE mo_jsb_impl_constants,     ONLY: false

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
#ifdef __QUINCY_NO_JSB__
  USE mo_jsb_var_class,          ONLY: t_jsb_struct, t_jsb_n_struct
#else
  USE mo_jsb_pool_class,         ONLY: t_jsb_pool, ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, VEG_TYPE
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_sb_memory, max_no_of_vars
#ifdef __QUINCY_NO_JSB__
  PUBLIC :: max_no_of_structs, max_no_of_n_structs
  PUBLIC :: t_sb_element_a, t_sb_pool_a                     ! mo_sb_memory_class_util
  PUBLIC :: Init_sb_pools_elements_a, &                     ! mo_sb_memory_class_util
            Deallocate_var_sb_pools_elements_a              ! mo_sb_memory_class_util
  PUBLIC :: sb_pool, sb_formation, sb_loss, sb_transport    ! mo_sb_memory_class_util
  PUBLIC :: sb_delta_pools, sb_mycorrhiza_export            ! mo_sb_memory_class_util
  INTEGER, PARAMETER :: max_no_of_vars      = 200  !> TYPE t_jsb_var
  INTEGER, PARAMETER :: max_no_of_structs   = 20   !> TYPE t_jsb_struct
  INTEGER, PARAMETER :: max_no_of_n_structs = 20   !> TYPE t_jsb_n_struct
#else
  PUBLIC :: t_sb_main_pool, t_sb_components_pool, t_sb_elements_pool
  INTEGER, PARAMETER :: max_no_of_vars      = 480  !> TYPE t_jsb_var | 'normal' var: ~ 200 | 4 pools: 240 var | 2 tissues: 12
#endif

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for sb memory
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_sb_memory

    !! ####################################################
    !> pools 3D
    !!
    TYPE(t_jsb_var_real3d)  :: &
                            nh4_solute                , &  !< ammonium in soil solution [mol m-3]
                            nh4_assoc                 , &  !< minerally associated ammonium [mol m-3]
                            no3_solute                , &  !< nitrate in soil solution [mol m-3]
                            noy                       , &  !< NO_y [mol m-3]
                            n2o                       , &  !< N2O [mol m-3]
                            n2                        , &  !< N2 [mol m-3]
                            nh4_n15_solute            , &  !< N15-ammonium in soil solution [mol m-3]
                            nh4_n15_assoc             , &  !< minerally associated N15-ammonium [mol m-3]
                            no3_n15_solute            , &  !< N15-nitrate in soil solution [mol m-3]
                            noy_n15                   , &  !< N15-NO_y [mol m-3]
                            n2o_n15                   , &  !< N15-N2O [mol m-3]
                            n2_n15                    , &  !< N15-N2 [mol m-3]
                            po4_assoc_fast            , &  !< fast minerally associated phosphate  [mol m-3]
                            po4_assoc_slow            , &  !< slow minerally associated phosphate  [mol m-3] 
                            po4_occluded              , &  !< occluded phosphate [mol m-3] 
                            po4_primary               , &  !< primary phosphate mineral [mol m-3] 
                            po4_solute                , &  !< phophate in soil solution [mol m-3] 
                            soil_litter_carbon_sl          !< (bulk) soil density of litter carbon [mol C m-3 soil]

    !! ####################################################
    !> transport fluxes 3D
    !!
    TYPE(t_jsb_var_real3d)  :: &
                            transport_nh4_solute         , &  !< transport of ammonium in soil solution [mol m-3 timestep-1]
                            transport_nh4_assoc          , &  !< transport of minerally associated ammonium [mol m-3 timestep-1]
                            transport_no3_solute         , &  !< transport of nitrate in soil solution [mol m-3 timestep-1]
                            transport_noy                , &  !< transport of NO_y [mol m-3 timestep-1]
                            transport_n2o                , &  !< transport of N2O [mol m-3 timestep-1]
                            transport_n2                 , &  !< transport of N2 [mol m-3 timestep-1]
                            transport_nh4_n15_solute     , &  !< transport of N15-ammonium in soil solution [mol m-3 timestep-1]
                            transport_nh4_n15_assoc      , &  !< transport of minerally associated N15-ammonium [mol m-3 timestep-1]
                            transport_no3_n15_solute     , &  !< transport of N15-nitrate in soil solution [mol m-3 timestep-1]
                            transport_noy_n15            , &  !< transport of N15-NO_y [mol m-3 timestep-1]
                            transport_n2o_n15            , &  !< transport of N15-N2O [mol m-3 timestep-1]
                            transport_n2_n15             , &  !< transport of N15-N2 [mol m-3 timestep-1]
                            transport_po4_assoc_fast     , &  !< transport of fast minerally associated phosphate  [mol m-3 timestep-1]
                            transport_po4_assoc_slow     , &  !< transport of slow minerally associated phosphate  [mol m-3 timestep-1] 
                            transport_po4_occluded       , &  !< transport of occluded phosphate [mol m-3 timestep-1] 
                            transport_po4_primary        , &  !< transport of primary phosphate [mol m-3 timestep-1] 
                            transport_po4_solute              !< transport of minerally associated phosphate  [mol m-3 timestep-1]

    !! ####################################################
    !> response functions 3D
    !!
    ! temperature response functions
    TYPE(t_jsb_var_real3d)  :: &
                            rtm_decomposition, &              !< rate modifier for decomposition processes [--]
                            rtm_depolymerisation, &           !< rate modifier for depolymerisation processes [--]
                            rtm_mic_uptake, &                 !< rate modifier for microbial uptake [--] 
                            rtm_plant_uptake, &               !< rate modifier for plant uptake [--] 
                            rtm_sorption, &                   !< rate modifier for sorption [--]
                            rtm_desorption, &                 !< rate modifier for desorption [--]
                            rtm_hsc, &                        !< rate modifier for half-saturation constants [--] 
                            rtm_nitrification, &              !< rate modifier for nitrification [--] 
                            rtm_denitrification, &            !< rate modifier for denitrification [--] 
                            rtm_gasdiffusion                  !< rate modifier for gasdiffusion [--] 

    ! moisture response functions
    TYPE(t_jsb_var_real3d)  :: &
                            rmm_decomposition, &              !< rate modifier for decomposition processes [--]
                            rmm_depolymerisation, &           !< rate modifier for depolymerisation processes [--]
                            rmm_mic_uptake, &                 !< rate modifier for microbial uptake [--] 
                            rmm_plant_uptake, &               !< rate modifier for plant uptake [--] 
                            rmm_sorption, &                   !< rate modifier for sorption [--]
                            rmm_desorption, &                 !< rate modifier for desorption [--]
                            rmm_hsc, &                        !< rate modifier for half-saturation constants [--] 
                            rmm_nitrification, &              !< rate modifier for nitrification [--] 
                            rmm_gasdiffusion                  !< rate modifier for gasdiffusion [--] 

    !! ####################################################
    !> fluxes 3D | respiration, mineralisation, uptake, nitrification, denitrification, ...
    !!
    TYPE(t_jsb_var_real3d)  :: &
                            het_respiration, &                    !< heterotrophic respiration  [µmol m-3 s-1]
                            het_respiration_c13, &                !< c13 heterotrophic respiration   [µmol m-3 s-1]
                            het_respiration_c14, &                !< c14 heterotrophic respiration   [µmol m-3 s-1]
                            myc_respiration, &                    !< mycorrhizae respiration  [µmol m-3 s-1]
                            myc_respiration_c13, &                !< c13 mycorrhizae respiration   [µmol m-3 s-1]
                            myc_respiration_c14, &                !< c14 mycorrhizae respiration   [µmol m-3 s-1]
                            net_mineralisation_nh4, &             !< net mineralisation of N    [µmol m-3 s-1]
                            net_mineralisation_nh4_n15, &         !< net mineralisation of N15  [µmol m-3 s-1] 
                            net_mineralisation_no3, &             !< net mineralisation of N    [µmol m-3 s-1]
                            net_mineralisation_no3_n15, &         !< net mineralisation of N15  [µmol m-3 s-1] 
                            net_mineralisation_po4, &             !< net mineralisation of P  [µmol m-3 s-1]
                            plant_uptake_nh4_sl, &                !< plant uptake of N    [µmol m-3 s-1]
                            plant_uptake_nh4_n15_sl, &            !< plant uptake of N15  [µmol m-3 s-1] 
                            plant_uptake_no3_sl, &                !< plant uptake of N    [µmol m-3 s-1]
                            plant_uptake_no3_n15_sl, &            !< plant uptake of N15  [µmol m-3 s-1] 
                            plant_uptake_po4_sl, &                !< plant uptake of P  [µmol m-3 s-1]
                            microbial_uptake_nh4_sl, &            !< microbial uptake of N    [µmol m-3 s-1]
                            microbial_uptake_nh4_n15_sl, &        !< microbial uptake of N15  [µmol m-3 s-1] 
                            microbial_uptake_no3_sl, &            !< microbial uptake of N    [µmol m-3 s-1]
                            microbial_uptake_no3_n15_sl, &        !< microbial uptake of N15  [µmol m-3 s-1] 
                            microbial_uptake_po4_sl, &            !< microbial uptake of P  [µmol m-3 s-1]
                            mycorrhiza_uptake_nh4_sl, &           !< mycorrhiza uptake of N    [µmol m-3 s-1]
                            mycorrhiza_uptake_nh4_n15_sl, &       !< mycorrhiza uptake of N15  [µmol m-3 s-1] 
                            mycorrhiza_uptake_no3_sl, &           !< mycorrhiza uptake of N    [µmol m-3 s-1]
                            mycorrhiza_uptake_no3_n15_sl, &       !< mycorrhiza uptake of N15  [µmol m-3 s-1] 
                            mycorrhiza_uptake_norg_sl, &          !< mycorrhiza uptake of N    [µmol m-3 s-1]
                            mycorrhiza_uptake_norg_n15_sl, &      !< mycorrhiza uptake of N15  [µmol m-3 s-1] 
                            mycorrhiza_uptake_po4_sl, &           !< mycorrhiza uptake of P  [µmol m-3 s-1] 
                            myc_export_n_tlabile_mavg_sl, &       !< time averaged export from mycorrhizae to plants of N [mol m-3]
                            myc_export_p_tlabile_mavg_sl, &       !< time averaged export from mycorrhizae to plants of P [mol m-3]
                            myc_export_c_tmyc_mavg_sl, &          !< time averaged export from mycorrhizae to plants of C [mol m-3]
                            myc_export_n_tmyc_mavg_sl, &          !< time averaged export from mycorrhizae to plants of N [mol m-3]
                            volatilisation_nh4, &                 !< volatilisation of NH4  [µmol m-3 s-1]
                            volatilisation_nh4_n15, &             !< volatilisation of 15NH4  [µmol m-3 s-1]
                            nitrification_no3, &                  !< nitrification NO3 [µmol m-3 s-1]
                            nitrification_no3_n15, &              !< nitrification 15NO3 [µmol m-3 s-1]
                            nitrification_noy, &                  !< nitrification to NOy  [µmol m-3 s-1]
                            nitrification_noy_n15, &              !< nitrification to 15NOy  [µmol m-3 s-1]
                            nitrification_n2o, &                  !< nitrification to N2O [µmol m-3 s-1]
                            nitrification_n2o_n15, &              !< nitrification to 15N2O [µmol m-3 s-1]
                            denitrification_noy, &                !< nitrification NOy [µmol m-3 s-1]
                            denitrification_n2o_n15, &            !< nitrification to 15N2O  [µmol m-3 s-1]
                            denitrification_n2o, &                !< nitrification to N2O  [µmol m-3 s-1]
                            denitrification_noy_n15, &            !< nitrification 15NOy [µmol m-3 s-1]
                            denitrification_n2, &                 !< denitrification to N2  [µmol m-3 s-1]
                            denitrification_n2_n15, &             !< denitrification to 15N2  [µmol m-3 s-1]
                            asymb_n_fixation, &                   !< asymboitic N fixation [µmol m-3 s-1]
                            asymb_n_fixation_n15, &               !< asymboitic 15N fixation [µmol m-3 s-1]
                            weathering_po4, &                     !< phosphorus weathering rate [µmol P m-3 s-1]
                            occlusion_po4, &                      !< phosphorus occlusion rate, assoc_slow --> occluded [µmol P m-3 s-1]
                            slow_exchange_po4, &                  !< phosphorus slow exchange rate, assoc_fast <--> assoc_slow [µmol P m-3 s-1]
                            fast_exchange_po4, &                  !< phosphorus fast exchange rate, solute <--> assoc_fast [µmol P m-3 s-1]
                            fast_exchange_nh4, &                  !< nh4 fast exchange rate, solute <--> nh4_assoc [µmol N m-3 s-1]
                            fast_adsorpt_po4, &                   !< phosphorus fast exchange rate, solute <--> assoc_fast [µmol P m-3 s-1]
                            fast_adsorpt_nh4, &                   !< nh4 fast adsorption rate, solute --> nh4_assoc [µmol P m-3 s-1]
                            biochem_mineralisation_po4, &         !< phosphorus biochemical mineralization rate [µmol P m-3 s-1]
                            gross_mineralisation_po4              !< phosphorus gross (biological) mineralization rate [µmol P m-3 s-1]

    !! ####################################################
    !> 2D transport fluxes for communication with groundwater and atmosphere
    !!
    TYPE(t_jsb_var_real2d)  :: & 
                            leaching_dom_carbon, &                !< land to groundwater flux of DOM-C [µmol m-2 s-1]
                            leaching_dom_nitrogen, &              !< land to groundwater flux of DOM-N [µmol m-2 s-1]
                            leaching_dom_phosphorus, &            !< land to groundwater flux of DOM-P [µmol m-2 s-1]
                            leaching_dom_carbon13, &              !< land to groundwater flux of DOM-C13 [µmol m-2 s-1]
                            leaching_dom_carbon14, &              !< land to groundwater flux of DOM-C14 [µmol m-2 s-1]
                            leaching_dom_nitrogen15, &            !< land to groundwater flux of DOM-N15 [µmol m-2 s-1]
                            leaching_nh4_solute, &                !< land to groundwater flux of NH4 [µmol m-2 s-1]
                            leaching_no3_solute, &                !< land to groundwater flux of NO3 [µmol m-2 s-1]
                            leaching_po4_solute, &                !< land to groundwater flux of PO4 [µmol m-2 s-1]
                            leaching_nh4_n15_solute, &            !< land to groundwater flux of 15NH4 [µmol m-2 s-1]
                            leaching_no3_n15_solute, &            !< land to groundwater flux of 15NO3 [µmol m-2 s-1]
                            emission_nh3, &                       !< land to atmosphere flux of NH3 [µmol m-2 s-1] 
                            emission_noy, &                       !< land to atmosphere flux of NOy [µmol m-2 s-1] 
                            emission_n2o, &                       !< land to atmosphere flux of N2O [µmol m-2 s-1]
                            emission_n2, &                        !< land to atmosphere flux of N2 [µmol m-2 s-1]
                            emission_nh3_n15, &                   !< land to atmosphere flux of 15NH3 [µmol m-2 s-1]
                            emission_noy_n15, &                   !< land to atmosphere flux of 15NOy [µmol m-2 s-1]
                            emission_n2o_n15, &                   !< land to atmosphere flux of 15N2O [µmol m-2 s-1]
                            emission_n2_n15                       !< land to atmosphere flux of 15N2 [µmol m-2 s-1]

    !! ####################################################
    !> 3D                       
    !!
    TYPE(t_jsb_var_real3d)  :: &
                            qmax_org, &                           !< sorption capacity of dom of mineral soil [kg C kg-1 mineral soil]
                            qmax_po4, &                           !< sorption capacity of po4 of mineral soil [mol P m-3 soil]
                            qmax_nh4, &                           !< sorption capacity of nh4 of mineral soil [mol N m-3 soil]
                            qmax_fast_po4, &                      !< sorption capacity of po4 of the Fe/Al pool [mol P m-3 soil]
                            qmax_slow_po4, &                      !< sorption capacity of po4 of the clay/silt pool [mol P m-3 soil]
                            km_fast_po4, &                        !< half-saturation concentration of po4 for the Fe/Al pool [mol P m-3 soil]
                            km_slow_po4, &                        !< half-saturation concentration of po4 for the clay/silt pool [mol P m-3 soil]
                            km_adsorpt_po4_sl, &                  !< half saturation of po4 sorption isotherm [mol P m-3 soil]
                            km_adsorpt_nh4_sl, &                  !< half saturation of nh4 sorption isotherm [mol N m-3 soil]
                            k_bioturb, &                          !< diffusion rate resulting from bioturbation
                            bulk_dens_corr_sl, &                  !< organic fraction corrected bulk-density of the soil [kg m-3]
                            bulk_soil_carbon_sl, &                !< bulk soil density of carbon [kg m-3] 
                            particle_fluxrate, &                  !< burial rate due to organic matter fall and decay [m s-1] 
                            anaerobic_volume_fraction_sl, &       !< anaerobic volume fraction [unitless]
                            microbial_cue_eff, &                  !< instantaneous microbial CUE given NP constraints [unitless]
                            microbial_cue_eff_tmic_mavg, &        !< time-averaged microbial CUE given NP constraints [unitless]
                            microbial_nue_eff, &                  !< instantaneous microbial NUE given NP constraints [unitless]
                            microbial_nue_eff_tmic_mavg, &        !< time-averaged microbial NUE given NP constraints [unitless]
                            microbial_pue_eff, &                  !< instantaneous microbial PUE given NP constraints [unitless]
                            microbial_pue_eff_tmic_mavg, &        !< time-averaged microbial PUE given NP constraints [unitless]
                            dom_cn, &                             !< dom CN ratio
                            dom_cp, &                             !< dom CP ratio
                            dom_cn_mavg, &                        !< time-averaged dom CN ratio
                            dom_cp_mavg                           !< time-averaged dom CP ratio
 

    !! ####################################################
    !> 3D soil properties
    !!
    TYPE(t_jsb_var_real3d)  :: &
                            ph_sl, &                              !< ph per soil layer
                            lithology_class_sl, &                 !< ..  [unitless]
                            vmax_weath_mineral_sl, &              !< weathering coefficient of mineral soil [mol P m-3 s-1]
                            Qmax_AlFe_cor, &                      !< correlation factor for phosphate qmax due to the unknown Al/Fe content
                            enzyme_frac_poly, &                   !< fraction of depolymerisation enzyme which depolymerises poly litter [unitless]
                            enzyme_frac_residue, &                !< fraction of depolymerisation enzyme which depolymerises microbial residue [unitless]
                            enzyme_frac_poly_c, &                 !< fraction of depolymerisation enzyme which depolymerises poly litter, determined by C [unitless]
                            enzyme_frac_poly_n, &                 !< fraction of depolymerisation enzyme which depolymerises poly litter, determined by N [unitless]
                            enzyme_frac_poly_p,&                  !< fraction of depolymerisation enzyme which depolymerises poly litter, determined by P [unitless]
                            enzyme_frac_poly_c_mavg, &            !< moving average of enzyme_frac_poly_c [unitless]
                            enzyme_frac_poly_n_mavg, &            !< moving average of enzyme_frac_poly_n [unitless]
                            enzyme_frac_poly_p_mavg               !< moving average of enzyme_frac_poly_p [unitless]

    !! ####################################################
    !> 2D total fluxes for mass balance checks in soil_biogeochemistry model
    !!
    TYPE(t_jsb_var_real2d)  :: & 
                            total_flux_carbon, &                  !< carbon mass balance checks in soil_biogeochemistry model [mol C m-2 s-1]
                            total_flux_nitrogen, &                !< nitrogen mass balance checks in soil_biogeochemistry model [mol N m-2 s-1]
                            total_flux_phosphorus, &              !< phosphorus mass balance checks in soil_biogeochemistry model [mol P m-2 s-1]
                            total_flux_carbon13, &                !< carbon13 mass balance checks in soil_biogeochemistry model [mol C13 m-2 s-1]
                            total_flux_carbon14, &                !< carbon14 mass balance checks in soil_biogeochemistry model [mol C14 m-2 s-1]
                            total_flux_nitrogen15                 !< nitrogen15 mass balance checks in soil_biogeochemistry model [mol N15 m-2 s-1]

    !! ####################################################
    !> pool structures
    !! 
    TYPE(t_quincy_struct_sb_element) :: &   ! simple structure for elements: t_jsb_struct
                                        delta_pools, &          !< for mass balance checks in soil_biogeochemistry model
                                        mycorrhiza_export       !< nutrient flux from mycorrhiza to host plant


    TYPE(t_quincy_struct_sb_pool)    :: &   ! nested pool TYPE t_jsb_n_struct
                                        pool, &                 !< pools of the sb
                                        formation, &            !< 
                                        loss, &                 !< 
                                        transport               !< 

  CONTAINS
    PROCEDURE :: Init => Init_sb_memory
  END TYPE t_sb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_memory_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the SB_ process
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_sb_memory(mem, prefix, suffix, lct_ids, model_id)

#ifdef __QUINCY_NO_JSB__
    USE mo_slm_varlist,       ONLY: jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2   ! , tables
#else
    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
#endif
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_sb_memory),   INTENT(inout), TARGET :: mem             !> SB_ memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids !(:)      !> vector of land cover type IDs
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ---------------------------
    TYPE(t_jsb_model),          POINTER :: model
    TYPE(t_jsb_grid),           POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid),          POINTER :: vgrid_surface                ! Vertical grid
    TYPE(t_jsb_vgrid),          POINTER :: vgrid_soil_sb                ! Vertical grid
#ifdef __QUINCY_NO_JSB__
    CHARACTER(len=30)                   :: unit_sb_pool
    CHARACTER(len=30)                   :: unit_sb_flux    
#else
    TYPE(t_sb_elements_pool),   POINTER :: sb_elements_pool             ! pools/fluxes like sb_delta_pools & sb_mycorrhiza_export
    TYPE(t_sb_components_pool), POINTER :: sb_components_pool           ! pools/fluxes like sb_pool & sb_formation
#endif
    INTEGER                             :: table
    CHARACTER(len=*), PARAMETER         :: routine = TRIM(modname)//':Init_sb_memory'
    ! ---------------------------
    model         => Get_model(model_id)
#ifdef __QUINCY_NO_JSB__
    table         = 1
#else
    table         =  tables(1)
#endif
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface => Get_vgrid('surface')
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')

#ifdef __QUINCY_NO_JSB__
    !> 0.9 init local variables
    !! 
    unit_sb_pool = "mol m-3"
    unit_sb_flux = "mol m-3 timestep-1"
#endif

    !> 1.0 sb memory init
    !! 
    !! 
    IF (.NOT. ASSOCIATED(mem%nh4_solute%ptr)) THEN

      !> add variables to soil memory soil\%VAR
      !! 
      !! soil\%VAR
      !!  
      CALL jsb_add_var(mem, 'nh4_solute', mem%nh4_solute, &
            hgrid, vgrid_soil_sb, &
            t_cf('nh4_solute', 'mol m-3', 'ammonium in soil solution'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'nh4_assoc', mem%nh4_assoc, &
            hgrid, vgrid_soil_sb, &
            t_cf('nh4_assoc', 'mol m-3', 'minerally associated ammonium'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'no3_solute', mem%no3_solute, &
            hgrid, vgrid_soil_sb, &
            t_cf('no3_solute', 'mol m-3', 'nitrate in soil solution'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'noy', mem%noy, &
            hgrid, vgrid_soil_sb, &
            t_cf('noy', 'mol m-3', 'NO_y'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n2o', mem%n2o, &
            hgrid, vgrid_soil_sb, &
            t_cf('n2o', 'mol m-3', 'N2O'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n2', mem%n2, &
            hgrid, vgrid_soil_sb, &
            t_cf('n2', 'mol m-3', 'N2'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'nh4_n15_solute', mem%nh4_n15_solute, &
            hgrid, vgrid_soil_sb, &
            t_cf('nh4_n15_solute', 'mol m-3', 'N15-ammonium in soil solution'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'nh4_n15_assoc', mem%nh4_n15_assoc, &
            hgrid, vgrid_soil_sb, &
            t_cf('nh4_n15_assoc', 'mol m-3', 'minerally associated N15-ammonium'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'no3_n15_solute', mem%no3_n15_solute, &
            hgrid, vgrid_soil_sb, &
            t_cf('no3_n15_solute', 'mol m-3', 'N15-nitrate in soil solution'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'noy_n15', mem%noy_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('noy_n15', 'mol m-3', 'N15-NO_y'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n2o_n15', mem%n2o_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('n2o_n15', 'mol m-3', 'N15-N2O'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n2_n15', mem%n2_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('n2_n15', 'mol m-3', 'N15-N2'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'po4_assoc_fast', mem%po4_assoc_fast, &
            hgrid, vgrid_soil_sb, &
            t_cf('po4_assoc_fast', 'mol m-3', 'fast minerally associated phosphate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'po4_assoc_slow', mem%po4_assoc_slow, &
            hgrid, vgrid_soil_sb, &
            t_cf('po4_assoc_slow', 'mol m-3', 'slow minerally associated phosphate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'po4_occluded', mem%po4_occluded, &
            hgrid, vgrid_soil_sb, &
            t_cf('po4_occluded', 'mol m-3', 'occluded phosphate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'po4_primary', mem%po4_primary, &
            hgrid, vgrid_soil_sb, &
            t_cf('po4_primary', 'mol m-3', 'primary phosphate mineral'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'po4_solute', mem%po4_solute, &
            hgrid, vgrid_soil_sb, &
            t_cf('po4_solute', 'mol m-3', 'phophate in soil solution'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_nh4_solute', mem%transport_nh4_solute, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_nh4_solute', 'mol m-3 timestep-1', 'transport of ammonium in soil solution'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_nh4_assoc', mem%transport_nh4_assoc, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_nh4_assoc', 'mol m-3 timestep-1', 'transport of minerally associated ammonium'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_no3_solute', mem%transport_no3_solute, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_no3_solute', 'mol m-3 timestep-1', 'transport of nitrate in soil solution'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_noy', mem%transport_noy, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_noy', 'mol m-3 timestep-1', 'transport of NO_y'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_n2o', mem%transport_n2o, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_n2o', 'mol m-3 timestep-1', 'transport of N2O'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_n2', mem%transport_n2, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_n2', 'mol m-3 timestep-1', 'transport of N2'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_nh4_n15_solute', mem%transport_nh4_n15_solute, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_nh4_n15_solute', 'mol m-3 timestep-1', 'transport of N15-ammonium in soil solution'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_nh4_n15_assoc', mem%transport_nh4_n15_assoc, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_nh4_n15_assoc', 'mol m-3 timestep-1', 'transport of minerally associated N15-ammonium'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_no3_n15_solute', mem%transport_no3_n15_solute, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_no3_n15_solute', 'mol m-3 timestep-1', 'transport of N15-nitrate in soil solution'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_noy_n15', mem%transport_noy_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_noy_n15', 'mol m-3 timestep-1', 'transport of N15-NO_y'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_n2o_n15', mem%transport_n2o_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_n2o_n15', 'mol m-3 timestep-1', 'transport of N15-N2O'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_n2_n15', mem%transport_n2_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_n2_n15', 'mol m-3 timestep-1', 'transport of N15-N2'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_po4_assoc_fast', mem%transport_po4_assoc_fast, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_po4_assoc_fast', 'mol m-3 timestep-1', 'transport of fast minerally associated phosphate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_po4_assoc_slow', mem%transport_po4_assoc_slow, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_po4_assoc_slow', 'mol m-3 timestep-1', 'transport of slow minerally associated phosphate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_po4_occluded', mem%transport_po4_occluded, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_po4_occluded', 'mol m-3 timestep-1', 'transport of occluded phosphate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_po4_primary', mem%transport_po4_primary, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_po4_primary', 'mol m-3 timestep-1', 'transport of primary phosphate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue  

      CALL jsb_add_var(mem, 'transport_po4_solute', mem%transport_po4_solute, &
            hgrid, vgrid_soil_sb, &
            t_cf('transport_po4_solute', 'mol m-3 timestep-1', 'transport of minerally associated phosphate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_decomposition', mem%rtm_decomposition, &
            hgrid, vgrid_soil_sb, &
            t_cf('rtm_decomposition', '--', 'rate modifier for decomposition processes'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_depolymerisation', mem%rtm_depolymerisation, &
            hgrid, vgrid_soil_sb, &
            t_cf('rtm_depolymerisation', '--', 'rate modifier for depolymerisation processes'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_mic_uptake', mem%rtm_mic_uptake, &
            hgrid, vgrid_soil_sb, &
            t_cf('rtm_mic_uptake', '--', 'rate modifier for microbial uptake'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_plant_uptake', mem%rtm_plant_uptake, &
            hgrid, vgrid_soil_sb, &
            t_cf('rtm_plant_uptake', '--', 'rate modifier for plant uptake'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_sorption', mem%rtm_sorption, &
            hgrid, vgrid_soil_sb, &
            t_cf('rtm_sorption', '--', 'rate modifier for sorption'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_desorption', mem%rtm_desorption, &
            hgrid, vgrid_soil_sb, &
            t_cf('rtm_desorption', '--', 'rate modifier for desorption'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_hsc', mem%rtm_hsc, &
            hgrid, vgrid_soil_sb, &
            t_cf('rtm_hsc', '--', 'rate modifier for half-saturation constants'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_nitrification', mem%rtm_nitrification, &
            hgrid, vgrid_soil_sb, &
            t_cf('rtm_nitrification', '--', 'rate modifier for nitrification'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_denitrification', mem%rtm_denitrification, &
            hgrid, vgrid_soil_sb, &
            t_cf('rtm_denitrification', '--', 'rate modifier for denitrification'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_gasdiffusion', mem%rtm_gasdiffusion, &
            hgrid, vgrid_soil_sb, &
            t_cf('rtm_gasdiffusion', '--', 'rate modifier for gasdiffusion'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_decomposition', mem%rmm_decomposition, &
            hgrid, vgrid_soil_sb, &
            t_cf('rmm_decomposition', '--', 'rate modifier for decomposition processes'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_depolymerisation', mem%rmm_depolymerisation, &
            hgrid, vgrid_soil_sb, &
            t_cf('rmm_depolymerisation', '--', 'rate modifier for depolymerisation processes'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_mic_uptake', mem%rmm_mic_uptake, &
            hgrid, vgrid_soil_sb, &
            t_cf('rmm_mic_uptake', '--', 'rate modifier for microbial uptake'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_plant_uptake', mem%rmm_plant_uptake, &
            hgrid, vgrid_soil_sb, &
            t_cf('rmm_plant_uptake', '--', 'rate modifier for plant uptake'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_sorption', mem%rmm_sorption, &
            hgrid, vgrid_soil_sb, &
            t_cf('rmm_sorption', '--', 'rate modifier for sorption'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_desorption', mem%rmm_desorption, &
            hgrid, vgrid_soil_sb, &
            t_cf('rmm_desorption', '--', 'rate modifier for desorption'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_hsc', mem%rmm_hsc, &
            hgrid, vgrid_soil_sb, &
            t_cf('rmm_hsc', '--', 'rate modifier for half-saturation constants'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_nitrification', mem%rmm_nitrification, &
            hgrid, vgrid_soil_sb, &
            t_cf('rmm_nitrification', '--', 'rate modifier for nitrification'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_gasdiffusion', mem%rmm_gasdiffusion, &
            hgrid, vgrid_soil_sb, &
            t_cf('rmm_gasdiffusion', '--', 'rate modifier for gasdiffusion'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'het_respiration', mem%het_respiration, &
            hgrid, vgrid_soil_sb, &
            t_cf('het_respiration', 'µmol m-3 s-1', 'heterotrophic respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'het_respiration_c13', mem%het_respiration_c13, &
            hgrid, vgrid_soil_sb, &
            t_cf('het_respiration_c13', 'µmol m-3 s-1', 'c13 heterotrophic respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'het_respiration_c14', mem%het_respiration_c14, &
            hgrid, vgrid_soil_sb, &
            t_cf('het_respiration_c14', 'µmol m-3 s-1', 'c14 heterotrophic respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'myc_respiration', mem%myc_respiration, &
            hgrid, vgrid_soil_sb, &
            t_cf('myc_respiration', 'µmol m-3 s-1', 'mycorrhizae respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'myc_respiration_c13', mem%myc_respiration_c13, &
            hgrid, vgrid_soil_sb, &
            t_cf('myc_respiration_c13', 'µmol m-3 s-1', 'c13 mycorrhizae respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'myc_respiration_c14', mem%myc_respiration_c14, &
            hgrid, vgrid_soil_sb, &
            t_cf('myc_respiration_c14', 'µmol m-3 s-1', 'c14 mycorrhizae respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_mineralisation_nh4', mem%net_mineralisation_nh4, &
            hgrid, vgrid_soil_sb, &
            t_cf('net_mineralisation_nh4', 'µmol m-3 s-1', 'net mineralisation of N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_mineralisation_nh4_n15', mem%net_mineralisation_nh4_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('net_mineralisation_nh4_n15', 'µmol m-3 s-1', 'net mineralisation of N15'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_mineralisation_no3', mem%net_mineralisation_no3, &
            hgrid, vgrid_soil_sb, &
            t_cf('net_mineralisation_no3', 'µmol m-3 s-1', 'net mineralisation of N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_mineralisation_no3_n15', mem%net_mineralisation_no3_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('net_mineralisation_no3_n15', 'µmol m-3 s-1', 'net mineralisation of N15'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_mineralisation_po4', mem%net_mineralisation_po4, &
            hgrid, vgrid_soil_sb, &
            t_cf('net_mineralisation_po4', 'µmol m-3 s-1', 'net mineralisation of P'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'plant_uptake_nh4_sl', mem%plant_uptake_nh4_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('plant_uptake_nh4_sl', 'µmol m-3 s-1', 'plant uptake of N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'plant_uptake_nh4_n15_sl', mem%plant_uptake_nh4_n15_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('plant_uptake_nh4_n15_sl', 'µmol m-3 s-1', 'plant uptake of N15'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'plant_uptake_no3_sl', mem%plant_uptake_no3_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('plant_uptake_no3_sl', 'µmol m-3 s-1', 'plant uptake of N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'plant_uptake_no3_n15_sl', mem%plant_uptake_no3_n15_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('plant_uptake_no3_n15_sl', 'µmol m-3 s-1', 'plant uptake of N15'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'plant_uptake_po4_sl', mem%plant_uptake_po4_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('plant_uptake_po4_sl', 'µmol m-3 s-1', 'plant uptake of P'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'microbial_uptake_nh4_sl', mem%microbial_uptake_nh4_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('microbial_uptake_nh4_sl', 'µmol m-3 s-1', 'microbial uptake of N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'microbial_uptake_nh4_n15_sl', mem%microbial_uptake_nh4_n15_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('microbial_uptake_nh4_n15_sl', 'µmol m-3 s-1', 'microbial uptake of N15'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'microbial_uptake_no3_sl', mem%microbial_uptake_no3_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('microbial_uptake_no3_sl', 'µmol m-3 s-1', 'microbial uptake of N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'microbial_uptake_no3_n15_sl', mem%microbial_uptake_no3_n15_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('microbial_uptake_no3_n15_sl', 'µmol m-3 s-1', 'microbial uptake of N15'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'microbial_uptake_po4_sl', mem%microbial_uptake_po4_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('microbial_uptake_po4_sl', 'µmol m-3 s-1', 'microbial uptake of P'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'mycorrhiza_uptake_nh4_sl', mem%mycorrhiza_uptake_nh4_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('mycorrhiza_uptake_nh4_sl', 'µmol m-3 s-1', 'mycorrhiza uptake of N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'mycorrhiza_uptake_nh4_n15_sl', mem%mycorrhiza_uptake_nh4_n15_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('mycorrhiza_uptake_nh4_n15_sl', 'µmol m-3 s-1', 'mycorrhiza uptake of N15'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'mycorrhiza_uptake_no3_sl', mem%mycorrhiza_uptake_no3_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('mycorrhiza_uptake_no3_sl', 'µmol m-3 s-1', 'mycorrhiza uptake of N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'mycorrhiza_uptake_no3_n15_sl', mem%mycorrhiza_uptake_no3_n15_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('mycorrhiza_uptake_no3_n15_sl', 'µmol m-3 s-1', 'mycorrhiza uptake of N15'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'mycorrhiza_uptake_norg_sl', mem%mycorrhiza_uptake_norg_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('mycorrhiza_uptake_norg_sl', 'µmol m-3 s-1', 'mycorrhiza uptake of N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'mycorrhiza_uptake_norg_n15_sl', mem%mycorrhiza_uptake_norg_n15_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('mycorrhiza_uptake_norg_n15_sl', 'µmol m-3 s-1', 'mycorrhiza uptake of N15'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'mycorrhiza_uptake_po4_sl', mem%mycorrhiza_uptake_po4_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('mycorrhiza_uptake_po4_sl', 'µmol m-3 s-1', 'mycorrhiza uptake of P'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'myc_export_n_tlabile_mavg_sl', mem%myc_export_n_tlabile_mavg_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('myc_export_n_tlabile_mavg_sl', 'mol m-3', 'time averaged mycorrhizae export of N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'myc_export_p_tlabile_mavg_sl', mem%myc_export_p_tlabile_mavg_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('myc_export_p_tlabile_mavg_sl', 'mol m-3', 'time averaged mycorrhizae export of P'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'myc_export_c_tmyc_mavg_sl', mem%myc_export_c_tmyc_mavg_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('myc_export_c_tmyc_mavg_sl', 'mol m-3', 'time averaged mycorrhizae export of C'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'myc_export_n_tmyc_mavg_sl', mem%myc_export_n_tmyc_mavg_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('myc_export_n_tmyc_mavg_sl', 'mol m-3', 'time averaged mycorrhizae export of N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'volatilisation_nh4', mem%volatilisation_nh4, &
            hgrid, vgrid_soil_sb, &
            t_cf('volatilisation_nh4', 'µmol m-3 s-1', 'volatilisation of NH4'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'volatilisation_nh4_n15', mem%volatilisation_nh4_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('volatilisation_nh4_n15', 'µmol m-3 s-1', 'volatilisation of 15NH4'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'nitrification_no3', mem%nitrification_no3, &
            hgrid, vgrid_soil_sb, &
            t_cf('nitrification_no3', 'µmol m-3 s-1', 'nitrification NO3'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'nitrification_no3_n15', mem%nitrification_no3_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('nitrification_no3_n15', 'µmol m-3 s-1', 'nitrification 15NO3'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'nitrification_noy', mem%nitrification_noy, &
            hgrid, vgrid_soil_sb, &
            t_cf('nitrification_noy', 'µmol m-3 s-1', 'nitrification to NOy'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'nitrification_noy_n15', mem%nitrification_noy_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('nitrification_noy_n15', 'µmol m-3 s-1', 'nitrification to 15NOy'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'nitrification_n2o', mem%nitrification_n2o, &
            hgrid, vgrid_soil_sb, &
            t_cf('nitrification_n2o', 'µmol m-3 s-1', 'nitrification to N2O'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'nitrification_n2o_n15', mem%nitrification_n2o_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('nitrification_n2o_n15', 'µmol m-3 s-1', 'nitrification to 15N2O'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'denitrification_noy', mem%denitrification_noy, &
            hgrid, vgrid_soil_sb, &
            t_cf('denitrification_noy', 'µmol m-3 s-1', 'nitrification NOy'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'denitrification_n2o_n15', mem%denitrification_n2o_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('denitrification_n2o_n15', 'µmol m-3 s-1', 'nitrification to 15N2O'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'denitrification_n2o', mem%denitrification_n2o, &
            hgrid, vgrid_soil_sb, &
            t_cf('denitrification_n2o', 'µmol m-3 s-1', 'nitrification to N2O'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'denitrification_noy_n15', mem%denitrification_noy_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('denitrification_noy_n15', 'µmol m-3 s-1', 'nitrification 15NOy'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'denitrification_n2', mem%denitrification_n2, &
            hgrid, vgrid_soil_sb, &
            t_cf('denitrification_n2', 'µmol m-3 s-1', 'denitrification to N2'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'denitrification_n2_n15', mem%denitrification_n2_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('denitrification_n2_n15', 'µmol m-3 s-1', 'denitrification to 15N2'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'asymb_n_fixation', mem%asymb_n_fixation, &
            hgrid, vgrid_soil_sb, &
            t_cf('asymb_n_fixation', 'µmol m-3 s-1', 'asymboitic N fixation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'asymb_n_fixation_n15', mem%asymb_n_fixation_n15, &
            hgrid, vgrid_soil_sb, &
            t_cf('asymb_n_fixation_n15', 'µmol m-3 s-1', 'asymboitic 15N fixation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'weathering_po4', mem%weathering_po4, &
            hgrid, vgrid_soil_sb, &
            t_cf('weathering_po4', 'µmol m-3 s-1', 'phosphorus weathering rate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'occlusion_po4', mem%occlusion_po4, &
            hgrid, vgrid_soil_sb, &
            t_cf('occlusion_po4', 'µmol m-3 s-1', 'phosphorus occlusion rate, assoc_slow --> occluded'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'slow_exchange_po4', mem%slow_exchange_po4, &
            hgrid, vgrid_soil_sb, &
            t_cf('slow_exchange_po4', 'µmol P m-3 s-1', 'phosphorus slow exchange rate, assoc_fast <--> assoc_slow'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'fast_exchange_po4', mem%fast_exchange_po4, &
            hgrid, vgrid_soil_sb, &
            t_cf('fast_exchange_po4', 'µmol P m-3 s-1', 'phosphorus fast exchange rate, solute <--> assoc_fast'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'fast_exchange_nh4', mem%fast_exchange_nh4, &
            hgrid, vgrid_soil_sb, &
            t_cf('fast_exchange_nh4', 'µmol N m-3 s-1', 'nh4 fast exchange rate, solute <--> nh4_assoc'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'fast_adsorpt_po4', mem%fast_adsorpt_po4, &
            hgrid, vgrid_soil_sb, &
            t_cf('fast_adsorpt_po4', 'µmol P m-3 s-1', 'phosphorus fast exchange rate, solute <--> assoc_fast'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'fast_adsorpt_nh4', mem%fast_adsorpt_nh4, &
            hgrid, vgrid_soil_sb, &
            t_cf('fast_adsorpt_nh4', 'µmol N m-3 s-1', 'nh4 fast adsorption rate, solute --> nh4_assoc'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'biochem_mineralisation_po4', mem%biochem_mineralisation_po4, &
            hgrid, vgrid_soil_sb, &
            t_cf('biochem_mineralisation_po4', 'µmol m-3 s-1', 'phosphorus biochemical mineralization rate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'gross_mineralisation_po4', mem%gross_mineralisation_po4, &
            hgrid, vgrid_soil_sb, &
            t_cf('gross_mineralisation_po4', 'µmol m-3 s-1', 'phosphorus gross (biological) mineralization rate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_dom_carbon', mem%leaching_dom_carbon, &
            hgrid, vgrid_surface, &
            t_cf('leaching_dom_carbon', 'µmol m-2 s-1', 'land to groundwater flux of DOM-C'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_dom_nitrogen', mem%leaching_dom_nitrogen, &
            hgrid, vgrid_surface, &
            t_cf('leaching_dom_nitrogen', 'µmol m-2 s-1', 'land to groundwater flux of DOM-N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_dom_phosphorus', mem%leaching_dom_phosphorus, &
            hgrid, vgrid_surface, &
            t_cf('leaching_dom_phosphorus', 'µmol m-2 s-1', 'land to groundwater flux of DOM-P'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_dom_carbon13', mem%leaching_dom_carbon13, &
            hgrid, vgrid_surface, &
            t_cf('leaching_dom_carbon13', 'µmol m-2 s-1', 'land to groundwater flux of DOM-C13'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_dom_carbon14', mem%leaching_dom_carbon14, &
            hgrid, vgrid_surface, &
            t_cf('leaching_dom_carbon14', 'µmol m-2 s-1', 'land to groundwater flux of DOM-C14'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_dom_nitrogen15', mem%leaching_dom_nitrogen15, &
            hgrid, vgrid_surface, &
            t_cf('leaching_dom_nitrogen15', 'µmol m-2 s-1', 'land to groundwater flux of DOM-N15'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_nh4_solute', mem%leaching_nh4_solute, &
            hgrid, vgrid_surface, &
            t_cf('leaching_nh4_solute', 'µmol m-2 s-1', 'land to groundwater flux of NH4'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue


      CALL jsb_add_var(mem, 'leaching_no3_solute', mem%leaching_no3_solute, &
            hgrid, vgrid_surface, &
            t_cf('leaching_no3_solute', 'µmol m-2 s-1', 'land to groundwater flux of NO3'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_po4_solute', mem%leaching_po4_solute, &
            hgrid, vgrid_surface, &
            t_cf('leaching_po4_solute', 'µmol m-2 s-1', 'land to groundwater flux of PO4'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_nh4_n15_solute', mem%leaching_nh4_n15_solute, &
            hgrid, vgrid_surface, &
            t_cf('leaching_nh4_n15_solute', 'µmol m-2 s-1', 'land to groundwater flux of 15NH4'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_no3_n15_solute', mem%leaching_no3_n15_solute, &
            hgrid, vgrid_surface, &
            t_cf('leaching_no3_n15_solute', 'µmol m-2 s-1', 'land to groundwater flux of 15NO3'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_nh3', mem%emission_nh3, &
            hgrid, vgrid_surface, &
            t_cf('emission_nh3', 'µmol m-2 s-1', 'land to atmosphere flux of NH3'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_noy', mem%emission_noy, &
            hgrid, vgrid_surface, &
            t_cf('emission_noy', 'µmol m-2 s-1', 'land to atmosphere flux of NOy'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_n2o', mem%emission_n2o, &
            hgrid, vgrid_surface, &
            t_cf('emission_n2o', 'µmol m-2 s-1', 'land to atmosphere flux of N2O'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_n2', mem%emission_n2, &
            hgrid, vgrid_surface, &
            t_cf('emission_n2', 'µmol m-2 s-1', 'land to atmosphere flux of N2'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_nh3_n15', mem%emission_nh3_n15, &
            hgrid, vgrid_surface, &
            t_cf('emission_nh3_n15', 'µmol m-2 s-1', 'land to atmosphere flux of 15NH3'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_noy_n15', mem%emission_noy_n15, &
            hgrid, vgrid_surface, &
            t_cf('emission_noy_n15', 'µmol m-2 s-1', 'land to atmosphere flux of 15NOy'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_n2o_n15', mem%emission_n2o_n15, &
            hgrid, vgrid_surface, &
            t_cf('emission_n2o_n15', 'µmol m-2 s-1', 'land to atmosphere flux of 15N2O'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_n2_n15', mem%emission_n2_n15, &
            hgrid, vgrid_surface, &
            t_cf('emission_n2_n15', 'µmol m-2 s-1', 'land to atmosphere flux of 15N2'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'qmax_org', mem%qmax_org, &
            hgrid, vgrid_soil_sb, &
            t_cf('qmax_org', 'kg C kg-1 mineral soil', 'sorption capacity of dom of mineral soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'qmax_po4', mem%qmax_po4, &
            hgrid, vgrid_soil_sb, &
            t_cf('qmax_po4', 'mol P m-3 soil', 'sorption capacity of po4 of mineral soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'km_adsorpt_po4_sl', mem%km_adsorpt_po4_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('km_adsorpt_po4_sl', 'mol P m-3 soil', 'half saturation of po4 sorption isotherm'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'qmax_nh4', mem%qmax_nh4, &
            hgrid, vgrid_soil_sb, &
            t_cf('qmax_nh4', 'mol N m-3 soil', 'sorption capacity of nh4 of mineral soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'km_adsorpt_nh4_sl', mem%km_adsorpt_nh4_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('km_adsorpt_nh4_sl', 'mol N m-3 soil', 'half saturation of nh4 sorption isotherm'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
			
      CALL jsb_add_var(mem, 'qmax_fast_po4', mem%qmax_fast_po4, &
            hgrid, vgrid_soil_sb, &
            t_cf('qmax_fast_po4', 'mol P m-3 soil', 'sorption capacity of po4 of the clay(Fe/Al) and SOM pool'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'qmax_slow_po4', mem%qmax_slow_po4, &
            hgrid, vgrid_soil_sb, &
            t_cf('qmax_slow_po4', 'mol P m-3 soil', 'sorption capacity of po4 of the silt/sand pool'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'km_fast_po4', mem%km_fast_po4, &
            hgrid, vgrid_soil_sb, &
            t_cf('km_fast_po4', 'mol P m-3 soil', 'half-saturation concentration of po4 for the clay(Fe/Al) and SOM pool'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'km_slow_po4', mem%km_slow_po4, &
            hgrid, vgrid_soil_sb, &
            t_cf('km_slow_po4', 'mol P m-3 soil', 'half-saturation concentration of po4 for the silt/sand pool'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'vmax_weath_mineral_sl', mem%vmax_weath_mineral_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('vmax_weath_mineral_sl', 'mol P m-3 s-1', 'weathering coefficient of mineral soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'k_bioturb', mem%k_bioturb, &
            hgrid, vgrid_soil_sb, &
            t_cf('k_bioturb', '', 'diffusion rate resulting from bioturbation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'bulk_dens_corr_sl', mem%bulk_dens_corr_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('bulk_dens_corr_sl', 'kg m-3', 'organic fraction corrected bulk-density of the soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'bulk_soil_carbon_sl', mem%bulk_soil_carbon_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('bulk_soil_carbon_sl', 'mol m-3', 'bulk soil density of non-litter carbon'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'soil_litter_carbon_sl', mem%soil_litter_carbon_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('soil_litter_carbon_sl', 'mol m-3', 'soil density of litter carbon'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'particle_fluxrate', mem%particle_fluxrate, &
            hgrid, vgrid_soil_sb, &
            t_cf('particle_fluxrate', 'm s-1', 'burial rate due to organic matter fall and decay'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'anaerobic_volume_fraction_sl', mem%anaerobic_volume_fraction_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('anaerobic_volume_fraction_sl', 'unitless', 'anaerobic volume fraction'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'microbial_cue_eff', mem%microbial_cue_eff, &
            hgrid, vgrid_soil_sb, &
            t_cf('microbial_cue_eff', 'unitless', 'instantaneous microbial CUE given NP constraints'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'microbial_cue_eff_tmic_mavg', mem%microbial_cue_eff_tmic_mavg, &
            hgrid, vgrid_soil_sb, &
            t_cf('microbial_cue_eff_tmic_mavg', 'unitless', 'time-averaged microbial CUE given NP constraints'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'microbial_nue_eff', mem%microbial_nue_eff, &
            hgrid, vgrid_soil_sb, &
            t_cf('microbial_nue_eff', 'unitless', 'instantaneous microbial NUE given NP constraints'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'microbial_nue_eff_tmic_mavg', mem%microbial_nue_eff_tmic_mavg, &
            hgrid, vgrid_soil_sb, &
            t_cf('microbial_nue_eff_tmic_mavg', 'unitless', 'time-averaged microbial NUE given NP constraints'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'microbial_pue_eff', mem%microbial_pue_eff, &
            hgrid, vgrid_soil_sb, &
            t_cf('microbial_pue_eff', 'unitless', 'instantaneous microbial PUE given NP constraints'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'microbial_pue_eff_tmic_mavg', mem%microbial_pue_eff_tmic_mavg, &
            hgrid, vgrid_soil_sb, &
            t_cf('microbial_pue_eff_tmic_mavg', 'unitless', 'time-averaged microbial CUE given NP constraints'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'dom_cn', mem%dom_cn, &
            hgrid, vgrid_soil_sb, &
            t_cf('dom_cn', '', 'dom CN ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'dom_cp', mem%dom_cp, &
            hgrid, vgrid_soil_sb, &
            t_cf('dom_cp', '', 'dom CP ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'dom_cn_mavg', mem%dom_cn_mavg, &
            hgrid, vgrid_soil_sb, &
            t_cf('dom_cn_mavg', '', 'time-averaged dom CN ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'dom_cp_mavg', mem%dom_cp_mavg, &
            hgrid, vgrid_soil_sb, &
            t_cf('dom_cp_mavg', '', 'time-averaged dom CP ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'ph_sl', mem%ph_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('ph_sl', '', 'ph per soil layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'lithology_class_sl', mem%lithology_class_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('lithology_class_sl', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
			
      CALL jsb_add_var(mem, 'Qmax_AlFe_cor', mem%Qmax_AlFe_cor, &
            hgrid, vgrid_soil_sb, &
            t_cf('Qmax_AlFe_cor', '', 'PO4 Qmax correction factor for clay/silt soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=1.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'enzyme_frac_poly', mem%enzyme_frac_poly, &
            hgrid, vgrid_soil_sb, &
            t_cf('enzyme_frac_poly', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'enzyme_frac_residue', mem%enzyme_frac_residue, &
            hgrid, vgrid_soil_sb, &
            t_cf('enzyme_frac_residue', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'enzyme_frac_poly_c', mem%enzyme_frac_poly_c, &
            hgrid, vgrid_soil_sb, &
            t_cf('enzyme_frac_poly_c', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'enzyme_frac_poly_n', mem%enzyme_frac_poly_n, &
            hgrid, vgrid_soil_sb, &
            t_cf('enzyme_frac_poly_n', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'enzyme_frac_poly_p', mem%enzyme_frac_poly_p, &
            hgrid, vgrid_soil_sb, &
            t_cf('enzyme_frac_poly_p', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'enzyme_frac_poly_c_mavg', mem%enzyme_frac_poly_c_mavg, &
            hgrid, vgrid_soil_sb, &
            t_cf('enzyme_frac_poly_c_mavg', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'enzyme_frac_poly_n_mavg', mem%enzyme_frac_poly_n_mavg, &
            hgrid, vgrid_soil_sb, &
            t_cf('enzyme_frac_poly_n_mavg', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'enzyme_frac_poly_p_mavg', mem%enzyme_frac_poly_p_mavg, &
            hgrid, vgrid_soil_sb, &
            t_cf('enzyme_frac_poly_p_mavg', 'unitless', 'moving average of enzyme_frac_poly_p'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'total_flux_carbon', mem%total_flux_carbon, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon', 'mol C m-2 s-1', 'carbon mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'total_flux_nitrogen', mem%total_flux_nitrogen, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_nitrogen', 'mol N m-2 s-1', 'nitrogen mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'total_flux_phosphorus', mem%total_flux_phosphorus, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_phosphorus', 'mol P m-2 s-1', 'phosphorus mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'total_flux_carbon13', mem%total_flux_carbon13, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon13', 'mol C13 m-2 s-1', 'carbon13 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'total_flux_carbon14', mem%total_flux_carbon14, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon14', 'mol C14 m-2 s-1', 'carbon14 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'total_flux_nitrogen15', mem%total_flux_nitrogen15, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_nitrogen15', 'mol N15 m-2 s-1', 'nitrogen15 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue


      !-----------------------------------------------------------------------------------------------------
      !> sb elements TYPE (t_jsb_struct)
      !! 
      !! jsb_add_var -> slm_add_struct()
      !! 
      !-----------------------------------------------------------------------------------------------------
      CALL jsb_add_var(mem, 'delta_pools', mem%delta_pools, &
            hgrid, vgrid_surface, &
            t_cf('delta_pools', TRIM(unit_sb_pool), 'Delta_pools element in soil biogeochemistry'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)

      CALL jsb_add_var(mem, 'mycorrhiza_export', mem%mycorrhiza_export, &
            hgrid, vgrid_surface, &
            t_cf('mycorrhiza_export', TRIM(unit_sb_flux), 'Mycorrhiza_export element in soil biogeochemistry'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)

      !> add variables to the elements
      !! 
      !! sb\%elements\%VAR
      !! 
      CALL mem%delta_pools%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_pool))
      CALL mem%mycorrhiza_export%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))


      !-----------------------------------------------------------------------------------------------------
      !> sb pools (nested) TYPE (t_jsb_n_struct)
      !! 
      !! jsb_add_var -> slm_add_n_struct()
      !! 
      !-----------------------------------------------------------------------------------------------------
      CALL jsb_add_var(mem, 'pool'           , mem%pool, &
            hgrid, vgrid_surface, &
            t_cf('pool', TRIM(unit_sb_pool), 'Pool soil biogeochemistry pool'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'formation'      , mem%formation, &
            hgrid, vgrid_surface, &
            t_cf('formation', TRIM(unit_sb_flux), 'Formation soil biogeochemistry flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'loss'           , mem%loss, &
            hgrid, vgrid_surface, &
            t_cf('loss', TRIM(unit_sb_flux), 'Loss soil biogeochemistry flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'transport'      , mem%transport, &
            hgrid, vgrid_surface, &
            t_cf('transport', TRIM(unit_sb_flux), 'Transport soil biogeochemistry flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)

      !> add element TYPEs to the pools
      !! 
      !! jsb_add_var -> slm_add_l2_struct()
      !! 
      CALL mem%pool%Add_elements(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_pool))
      CALL mem%formation%Add_elements(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%loss%Add_elements(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%transport%Add_elements(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))

      !> add variables to the elements in the pools
      !!
      !! sb\%pool\%elements\%VAR
      !! 
      CALL mem%pool%dom%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_pool))
      CALL mem%pool%dom_assoc%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_pool))
      CALL mem%pool%soluable_litter%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_pool))
      CALL mem%pool%polymeric_litter%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_pool))
      CALL mem%pool%woody_litter%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_pool))
      CALL mem%pool%fungi%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_pool))
      CALL mem%pool%mycorrhiza%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_pool))
      CALL mem%pool%microbial%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_pool))
      CALL mem%pool%residue%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_pool))
      CALL mem%pool%residue_assoc%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_pool))

      CALL mem%formation%dom%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%formation%dom_assoc%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%formation%soluable_litter%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%formation%polymeric_litter%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%formation%woody_litter%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%formation%fungi%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%formation%mycorrhiza%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%formation%microbial%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%formation%residue%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%formation%residue_assoc%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      
      CALL mem%loss%dom%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%loss%dom_assoc%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%loss%soluable_litter%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%loss%polymeric_litter%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%loss%woody_litter%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%loss%fungi%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%loss%mycorrhiza%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%loss%microbial%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%loss%residue%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%loss%residue_assoc%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      
      CALL mem%transport%dom%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%transport%dom_assoc%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%transport%soluable_litter%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%transport%polymeric_litter%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%transport%woody_litter%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%transport%fungi%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%transport%mycorrhiza%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%transport%microbial%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%transport%residue%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))
      CALL mem%transport%residue_assoc%Add_vars(prefix, suffix, lct_ids, model_id, TRIM(unit_sb_flux))

    END IF

  END SUBROUTINE Init_sb_memory

END MODULE mo_sb_memory_class
