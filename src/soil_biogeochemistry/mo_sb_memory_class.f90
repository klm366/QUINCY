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
#ifdef __QUINCY_NO_JSB__
#else
  USE mo_util,                   ONLY: One_of
#endif

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
#ifdef __QUINCY_NO_JSB__
  USE mo_jsb_var_class,          ONLY: t_jsb_struct, t_jsb_n_struct
#else
  USE mo_jsb_pool_class,         ONLY: ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15
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

    !> pools & fluxes
    !! 
    ! the sb_main_pool does contain all other pools & fluxes 
    TYPE(t_sb_main_pool), POINTER :: sb_main_pool


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
  !! pools and 
  !! 2d/3d variables
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


    !> add memory only for LAND & PFT tiles
    !!
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. &
       & One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN

      !=============================================================================================================================
      !> pools & fluxes
      !=============================================================================================================================

      !-----------------------------------------------------------------------------------------------------------------------------
      !> populate all_elements(:) with the element_id of all elements
      !! QUINCY is using always all of the below six elements
      IF (.NOT. ALLOCATED(all_elements)) THEN
        ALLOCATE(all_elements(6))
        all_elements(:) = [ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15]
      END IF

      !-----------------------------------------------------------------------------------------------------------------------------
      !> create and populate the main pool structure mem%pools that contains all other pools
      !! note: no spaces in 'name' and/or 'short name'
      ALLOCATE(t_sb_main_pool :: mem%pools)
      CALL mem%pools%Init(name='soil_biogeochemistry_main_pool', shortname='sb_main_pool', l_elements=.FALSE.)

      SELECT TYPE (selected => mem%pools)
      CLASS IS (t_sb_main_pool)
        mem%sb_main_pool => selected
      END SELECT

      !-----------------------------------------------------------------------------------------------------------------------------
      !> create the soil biogeochemistry components pools & fluxes (sb_components_pool) 
      !! sb_pool / sb_formation / ...
      ! create pools & fluxes
      ! "shortname" is used for model output names
      ! sb_pool
      sb_components_pool => Create_sb_components_pool(name='soil_biogeochemistry_pool', shortname='sb_pool', unit='mol m-3', &
        & lct_ids=lct_ids, pool_suffix='sbpool') ! l_elements=.TRUE. by default
      CALL mem%sb_main_pool%Add_pool(sb_components_pool)
      mem%sb_main_pool%sb_pool => sb_components_pool
      ! sb_formation
      sb_components_pool => Create_sb_components_pool(name='sb_formation_flux', shortname='sb_formation', unit='mol m-3', &
        & lct_ids=lct_ids, pool_suffix='sbform') ! l_elements=.TRUE. by default
      CALL mem%sb_main_pool%Add_pool(sb_components_pool)
      mem%sb_main_pool%sb_formation => sb_components_pool
      ! sb_loss
      sb_components_pool => Create_sb_components_pool(name='sb_loss_flux', shortname='sb_loss', unit='mol m-3', &
        & lct_ids=lct_ids, pool_suffix='sbloss') ! l_elements=.TRUE. by default
      CALL mem%sb_main_pool%Add_pool(sb_components_pool)
      mem%sb_main_pool%sb_loss => sb_components_pool
      ! sb_transport
      sb_components_pool => Create_sb_components_pool(name='sb_transport_flux', shortname='sb_transport', unit='mol m-3', &
        & lct_ids=lct_ids, pool_suffix='sbtran') ! l_elements=.TRUE. by default
      CALL mem%sb_main_pool%Add_pool(sb_components_pool)
      mem%sb_main_pool%sb_transport => sb_components_pool
      NULLIFY(sb_components_pool)

      ! init and add element pools to the sb_main_pool
      ! (i.e., add them at the same level as the above sb_components_pool sb_pool / sb_formation / ...)
      ! "shortname" is used for model output names
      ! sb_delta_pools
      ALLOCATE(sb_elements_pool)
      CALL sb_elements_pool%Init(name='sb_delta_pools_elem_pool', shortname='sb_delta_pools', &
        &                        l_elements=.TRUE., element_unit='mol m-3')
      CALL mem%sb_main_pool%Add_pool(sb_elements_pool)
      mem%sb_main_pool%sb_delta_pools => sb_elements_pool
      NULLIFY(sb_elements_pool)
      ! sb_mycorrhiza_export
      ALLOCATE(sb_elements_pool)
      CALL sb_elements_pool%Init(name='sb_mycorrhiza_export_elem_flux', shortname='sb_mycorrhiza_export', &
        &                        l_elements=.TRUE., element_unit='mol m-3')
      CALL mem%sb_main_pool%Add_pool(sb_elements_pool)
      mem%sb_main_pool%sb_mycorrhiza_export => sb_elements_pool
      NULLIFY(sb_elements_pool)

      ! set paths of all above initialized & added pools (type-bound procedure of t_jsb_pool)
      CALL mem%pools%Set_paths()

      ! add all variables in the sb_components_pool to the SB_ memory
      CALL mem%Add_var(mem%sb_main_pool%sb_pool, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )
      CALL mem%Add_var(mem%sb_main_pool%sb_formation, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )                                  ! no restart files needed for fluxes
      CALL mem%Add_var(mem%sb_main_pool%sb_loss, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )                                  ! no restart files needed for fluxes
      CALL mem%Add_var(mem%sb_main_pool%sb_transport, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )                                  ! no restart files needed for fluxes

      ! add all variables in the sb_elements_pool to the SB_ memory
      CALL mem%Add_var(mem%sb_main_pool%sb_delta_pools, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )
      CALL mem%Add_var(mem%sb_main_pool%sb_mycorrhiza_export, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )                                  ! no restart files needed for fluxes

      ! output to logfile all pools & their variables at model-init 
      ! (maybe use for debugging; its quite a long list)
      !IF (.NOT. ASSOCIATED(mem%parent)) THEN
      !  CALL mem%sb_main_pool%Print()
      !  ! CALL mem%sb_main_pool%sb_pool%Print()
      !END IF

      !=============================================================================================================================
      !> 2D & 3D variables
      !=============================================================================================================================

      CALL mem%Add_var('total_flux_carbon', mem%total_flux_carbon, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon', 'mol C m-2 s-1', 'carbon mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_nitrogen', mem%total_flux_nitrogen, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_nitrogen', 'mol N m-2 s-1', 'nitrogen mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_phosphorus', mem%total_flux_phosphorus, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_phosphorus', 'mol P m-2 s-1', 'phosphorus mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_carbon13', mem%total_flux_carbon13, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon13', 'mol C13 m-2 s-1', 'carbon13 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_carbon14', mem%total_flux_carbon14, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon14', 'mol C14 m-2 s-1', 'carbon14 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_nitrogen15', mem%total_flux_nitrogen15, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_nitrogen15', 'mol N15 m-2 s-1', 'nitrogen15 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('ph_sl', mem%ph_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('ph_sl', '', 'ph per soil layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue
    ENDIF

  END SUBROUTINE Init_sb_memory

END MODULE mo_sb_memory_class
