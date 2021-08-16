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
!! @brief Contains the routines for the Task update_sb_pools
!!
!! @author
!!  Silvia Caldararu, Sönke Zaehle 
!!
!! @par Revision History
!! First version by Sönke Zaehle (20170105): copyied from vegetation/mo_veg_soil.f90
!!
MODULE mo_sb_update_pools

#include "dsl_jsbach.inc"

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message, message_text, finish
  USE mo_jsb_math_constants,      ONLY: eps8

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_pools


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_update_pools'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil biogeochemical pools 
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sb_pools(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element    
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_jsb_math_constants,    ONLY: dtime, eps8
    USE mo_jsb_phys_constants,    ONLY: lambda_C14
    USE mo_sb_constants
    USE mo_sb_jsm_processes,      ONLY: calc_bulk_soil_carbon, calc_bulk_density_correction, calc_Psorption_parameter, &
                                        calc_qmax_bulk_density_correction, calc_fast_po4, convert_s1d_to_s0d_var_soilproperty
    USE mo_sb_jsm_transport,      ONLY: calc_particle_fluxrate, calc_particle_transport_wrapper
    USE mo_isotope_util,          ONLY: calc_mixing_ratio_N15N14
    USE mo_sse_util,              ONLY: calc_qmax_texture

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)
 
    ! Use of process memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4slm_Use_memory_a_sb
    USE mo_sb_memory_class,       ONLY: sb_pool, sb_formation, sb_loss, sb_transport, sb_delta_pools, sb_mycorrhiza_export
    dsl4jsb_Use_memory(SSE_)       !  USE mo_sse_memory_class,     ONLY: t_sse_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_sse_memory
    dsl4slm_Use_memory_a_veg
    USE mo_veg_memory_class,      ONLY: veg_litterfall, veg_exudation, veg_pool
    dsl4jsb_Use_memory(A2L_) 

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib               !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER       :: vgrid_soil_sb        !< Vertical grid
    INTEGER                               :: nsoil_sb             !< number of soil layers as used/defined by the SB_ process
    TYPE(t_sb_pool_a)                     :: par_transport, sb_last_timestep
    INTEGER                               :: iblk, ics, ice, nc
    INTEGER                               :: isoil
    REAL(wp), ALLOCATABLE, DIMENSION(:)   :: root_depth_avg   ! Dim: nc
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: p_deposition_sl, &
                                             f_myc_sl, &
                                             hlp1, hlp2, hlp3
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: km_adsorpt_po4_act, qmax_po4_mineral_act, & 
                                             km_fast_po4_act, qmax_fast_po4_act, & 
                                             km_slow_po4_act, qmax_slow_po4_act, & 
                                             km_adsorpt_nh4_act, qmax_nh4_mineral_act, & 
                                             k_desorpt_po4_act, k_adsorpt_po4_act, &
                                             rtm_sorption_act, rmm_sorption_act, &
                                             rtm_desorption_act, rmm_desorption_act, &
                                             partition_coef
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_pools'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)                               ! TYPE(t_sb_memory),       POINTER :: sb__mem
    dsl4jsb_Def_memory(SSE_)                             ! TYPE(t_sse_memory),     POINTER :: sse__mem
    dsl4jsb_Def_memory(VEG_)                              ! TYPE(t_sse_memory),     POINTER :: sse__mem
    dsl4jsb_Def_memory(A2L_) 
    ! Declare pointers to variables in memory
    ! A2L
    dsl4jsb_Real2D_onChunk :: nhx_deposition
    dsl4jsb_Real2D_onChunk :: noy_deposition
    dsl4jsb_Real2D_onChunk :: nhx_n15_deposition
    dsl4jsb_Real2D_onChunk :: noy_n15_deposition
    dsl4jsb_Real2D_onChunk :: p_deposition
    ! VEG_
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    dsl4jsb_Real2D_onChunk :: f_n_demand
    dsl4jsb_Real2D_onChunk :: f_p_demand
    ! SOIL
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: bulk_dens_sl 
    dsl4jsb_Real3D_onChunk :: percolation_sl
    dsl4jsb_Real3D_onChunk :: qmax_org_min_sl
    dsl4jsb_Real3D_onChunk :: qmax_po4_min_sl
    dsl4jsb_Real3D_onChunk :: qmax_nh4_min_sl
    dsl4jsb_Real3D_onChunk :: qmax_po4_om_sl
    dsl4jsb_Real3D_onChunk :: sand_sl
    dsl4jsb_Real3D_onChunk :: silt_sl
    dsl4jsb_Real3D_onChunk :: clay_sl
    dsl4jsb_Real3D_onChunk :: volume_min_sl
    ! SB_ 2D
    dsl4jsb_Real2D_onChunk :: total_flux_carbon
    dsl4jsb_Real2D_onChunk :: total_flux_nitrogen
    dsl4jsb_Real2D_onChunk :: total_flux_phosphorus
    dsl4jsb_Real2D_onChunk :: total_flux_carbon13
    dsl4jsb_Real2D_onChunk :: total_flux_carbon14
    dsl4jsb_Real2D_onChunk :: total_flux_nitrogen15
    dsl4jsb_Real2D_onChunk :: leaching_nh4_solute
    dsl4jsb_Real2D_onChunk :: leaching_nh4_n15_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_n15_solute
    dsl4jsb_Real2D_onChunk :: leaching_po4_solute
    dsl4jsb_Real2D_onChunk :: leaching_dom_carbon
    dsl4jsb_Real2D_onChunk :: leaching_dom_nitrogen
    dsl4jsb_Real2D_onChunk :: leaching_dom_phosphorus
    dsl4jsb_Real2D_onChunk :: leaching_dom_carbon13
    dsl4jsb_Real2D_onChunk :: leaching_dom_carbon14
    dsl4jsb_Real2D_onChunk :: leaching_dom_nitrogen15
    dsl4jsb_Real2D_onChunk :: emission_nh3
    dsl4jsb_Real2D_onChunk :: emission_nh3_n15
    dsl4jsb_Real2D_onChunk :: emission_noy
    dsl4jsb_Real2D_onChunk :: emission_noy_n15
    dsl4jsb_Real2D_onChunk :: emission_n2o
    dsl4jsb_Real2D_onChunk :: emission_n2o_n15
    dsl4jsb_Real2D_onChunk :: emission_n2
    dsl4jsb_Real2D_onChunk :: emission_n2_n15
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk :: het_respiration
    dsl4jsb_Real3D_onChunk :: het_respiration_c13
    dsl4jsb_Real3D_onChunk :: het_respiration_c14
    dsl4jsb_Real3D_onChunk :: myc_respiration
    dsl4jsb_Real3D_onChunk :: myc_respiration_c13
    dsl4jsb_Real3D_onChunk :: myc_respiration_c14
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4 
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: nitrification_no3 
    dsl4jsb_Real3D_onChunk :: nitrification_no3_n15 
    dsl4jsb_Real3D_onChunk :: nitrification_noy 
    dsl4jsb_Real3D_onChunk :: nitrification_noy_n15 
    dsl4jsb_Real3D_onChunk :: nitrification_n2o 
    dsl4jsb_Real3D_onChunk :: nitrification_n2o_n15
    dsl4jsb_Real3D_onChunk :: denitrification_noy 
    dsl4jsb_Real3D_onChunk :: denitrification_noy_n15 
    dsl4jsb_Real3D_onChunk :: denitrification_n2o 
    dsl4jsb_Real3D_onChunk :: denitrification_n2o_n15 
    dsl4jsb_Real3D_onChunk :: denitrification_n2 
    dsl4jsb_Real3D_onChunk :: denitrification_n2_n15 
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation 
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation_n15 
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: weathering_po4
    dsl4jsb_Real3D_onChunk :: occlusion_po4
    dsl4jsb_Real3D_onChunk :: slow_exchange_po4
    dsl4jsb_Real3D_onChunk :: fast_exchange_po4
    dsl4jsb_Real3D_onChunk :: fast_exchange_nh4
    dsl4jsb_Real3D_onChunk :: biochem_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: gross_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: particle_fluxrate 
    dsl4jsb_Real3D_onChunk :: transport_nh4_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_assoc 
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_assoc 
    dsl4jsb_Real3D_onChunk :: transport_no3_solute
    dsl4jsb_Real3D_onChunk :: transport_no3_n15_solute 
    dsl4jsb_Real3D_onChunk :: transport_po4_solute 
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_fast 
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_slow 
    dsl4jsb_Real3D_onChunk :: transport_po4_occluded
    dsl4jsb_Real3D_onChunk :: transport_po4_primary
    dsl4jsb_Real3D_onChunk :: transport_noy
    dsl4jsb_Real3D_onChunk :: transport_noy_n15
    dsl4jsb_Real3D_onChunk :: transport_n2o
    dsl4jsb_Real3D_onChunk :: transport_n2o_n15
    dsl4jsb_Real3D_onChunk :: transport_n2
    dsl4jsb_Real3D_onChunk :: transport_n2_n15 
    dsl4jsb_Real3D_onChunk :: nh4_solute 
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: nh4_assoc 
    dsl4jsb_Real3D_onChunk :: nh4_n15_assoc
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: no3_n15_solute 
    dsl4jsb_Real3D_onChunk :: po4_solute 
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast 
    dsl4jsb_Real3D_onChunk :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: po4_occluded
    dsl4jsb_Real3D_onChunk :: po4_primary
    dsl4jsb_Real3D_onChunk :: noy 
    dsl4jsb_Real3D_onChunk :: noy_n15
    dsl4jsb_Real3D_onChunk :: n2o 
    dsl4jsb_Real3D_onChunk :: n2o_n15
    dsl4jsb_Real3D_onChunk :: n2
    dsl4jsb_Real3D_onChunk :: n2_n15
    dsl4jsb_Real3D_onChunk :: rtm_sorption
    dsl4jsb_Real3D_onChunk :: rtm_desorption
    dsl4jsb_Real3D_onChunk :: rmm_sorption
    dsl4jsb_Real3D_onChunk :: rmm_desorption
    dsl4jsb_Real3D_onChunk :: qmax_org
    dsl4jsb_Real3D_onChunk :: qmax_po4
    dsl4jsb_Real3D_onChunk :: qmax_fast_po4
    dsl4jsb_Real3D_onChunk :: qmax_slow_po4
    dsl4jsb_Real3D_onChunk :: km_fast_po4
    dsl4jsb_Real3D_onChunk :: km_slow_po4
    dsl4jsb_Real3D_onChunk :: km_adsorpt_po4_sl
    dsl4jsb_Real3D_onChunk :: qmax_nh4
    dsl4jsb_Real3D_onChunk :: km_adsorpt_nh4_sl
    dsl4jsb_Real3D_onChunk :: dom_cn
    dsl4jsb_Real3D_onChunk :: dom_cp
    dsl4jsb_Real3D_onChunk :: ph_sl
    dsl4jsb_Real3D_onChunk :: Qmax_AlFe_cor
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onChunk :: bulk_soil_carbon_sl
    dsl4jsb_Real3D_onChunk :: soil_litter_carbon_sl
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
    ! Get process memories
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(A2L_) 
    ! Allocate local allocatable variables (deallocate below at the end of this routine)    
    ALLOCATE(p_deposition_sl(nc, nsoil_sb))
    ALLOCATE(f_myc_sl(nc, nsoil_sb))
    ALLOCATE(hlp1(nc, nsoil_sb))
    ALLOCATE(hlp2(nc, nsoil_sb))
    ALLOCATE(hlp3(nc, nsoil_sb))
    ALLOCATE(km_adsorpt_po4_act(nc, nsoil_sb))
    ALLOCATE(qmax_po4_mineral_act(nc, nsoil_sb))
    ALLOCATE(km_fast_po4_act(nc, nsoil_sb))
    ALLOCATE(qmax_fast_po4_act(nc, nsoil_sb))
    ALLOCATE(km_slow_po4_act(nc, nsoil_sb))
    ALLOCATE(qmax_slow_po4_act(nc, nsoil_sb))
    ALLOCATE(km_adsorpt_nh4_act(nc, nsoil_sb))
    ALLOCATE(qmax_nh4_mineral_act(nc, nsoil_sb))
    ALLOCATE(k_desorpt_po4_act(nc, nsoil_sb))
    ALLOCATE(k_adsorpt_po4_act(nc, nsoil_sb))
    ALLOCATE(rtm_sorption_act(nc, nsoil_sb))
    ALLOCATE(rmm_sorption_act(nc, nsoil_sb))
    ALLOCATE(rtm_desorption_act(nc, nsoil_sb))
    ALLOCATE(rmm_desorption_act(nc, nsoil_sb))
    ALLOCATE(partition_coef(nc, nsoil_sb))
    ALLOCATE(root_depth_avg(nc))
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(A2L_,      nhx_deposition)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      noy_deposition)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      nhx_n15_deposition)          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      noy_n15_deposition)          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      p_deposition)                ! in
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(VEG_,      root_fraction_sl)            ! in
    dsl4jsb_Get_var2D_onChunk(VEG_,      f_n_demand)                  ! in
    dsl4jsb_Get_var2D_onChunk(VEG_,      f_p_demand)                  ! in
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SSE_,     soil_depth_sl)               ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,     bulk_dens_sl)                ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,     percolation_sl)              ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,     qmax_org_min_sl)             ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,     qmax_po4_min_sl)             ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,     qmax_nh4_min_sl)             ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,     qmax_po4_om_sl)              ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,     sand_sl)                     ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,     silt_sl)                     ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,     clay_sl)                     ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,     volume_min_sl)               ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon)           ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_phosphorus)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon13)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon14)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen15)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_solute)         ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_n15_solute)     ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_solute)         ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_n15_solute)     ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_po4_solute)         ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_carbon)         ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_nitrogen)       ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_phosphorus)     ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_carbon13)       ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_carbon14)       ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_nitrogen15)     ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_nh3)                ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_nh3_n15)            ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy)                ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy_n15)            ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o)                ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o_n15)            ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2)                 ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2_n15)             ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration)             ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c13)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c14)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       myc_respiration)             ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       myc_respiration_c13)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       myc_respiration_c14)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4_n15)  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3_n15)  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_po4)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3)           ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3_n15)       ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy)           ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy_n15)       ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o)           ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o_n15)       ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy)         ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy_n15)     ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o_n15)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2)          ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2_n15)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_nh4)          ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_nh4_n15)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation)            ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation_n15)        ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_nh4_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_no3_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_po4_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_nh4_n15_sl) ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_no3_n15_sl) ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_sl)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_n15_sl)! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_sl)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_n15_sl)! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_norg_sl)   ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_norg_n15_sl) ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_po4_sl)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       particle_fluxrate)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_solute)        ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_solute)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_assoc)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_assoc)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_solute)        ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_n15_solute)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_solute)        ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_fast)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_slow)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_occluded)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_primary)       ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy)               ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy_n15)           ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o)               ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o_n15)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2)                ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2_n15)            ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_assoc)                   ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_assoc)               ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_fast)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_slow)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_occluded)                ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_primary)                 ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       noy)                         ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       noy_n15)                     ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o)                         ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o_n15)                     ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       n2)                          ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       n2_n15)                      ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       weathering_po4)              ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       occlusion_po4)               ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       slow_exchange_po4)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       fast_exchange_po4)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       fast_exchange_nh4)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       biochem_mineralisation_po4)  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       gross_mineralisation_po4)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_sorption)                ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_desorption)              ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_sorption)                ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_desorption)              ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_org)                    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_po4)                    ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_fast_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_slow_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       km_fast_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       km_slow_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       km_adsorpt_po4_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_nh4)                    ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       km_adsorpt_nh4_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       dom_cn)                      ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       dom_cp)                      ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       ph_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       Qmax_AlFe_cor)
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_soil_carbon_sl)         ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       soil_litter_carbon_sl)       ! inout
    ! Get process memory pool structures
    dsl4slm_get_vegpool_var(pool,fine_root,carbon)
#define Get_veg_litterfall
#define Get_veg_exudation
#define Get_sb_pool
#define Get_sb_formation
#define Get_sb_loss
#define Get_sb_transport
#define Get_sb_delta_pools
#define Get_sb_mycorrhiza_export
#include "slm_pools_a_incl.f90"
#undef Get_veg_litterfall
#undef Get_veg_exudation
#undef Get_sb_pool
#undef Get_sb_formation
#undef Get_sb_loss
#undef Get_sb_transport
#undef Get_sb_delta_pools
#undef Get_sb_mycorrhiza_export



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 1.0 correct the plant uptake of nh4, no3 and po4 for soil model
    !! 
    IF (TRIM(model%config%slm_model_name)=='soil') THEN
        DO isoil = 1,nsoil_sb
          plant_uptake_nh4_sl(:,isoil)     = plant_uptake_nh4_sl(:,isoil)     * f_n_demand(:)
          plant_uptake_nh4_n15_sl(:,isoil) = plant_uptake_nh4_n15_sl(:,isoil) * f_n_demand(:)
          plant_uptake_no3_sl(:,isoil)     = plant_uptake_no3_sl(:,isoil)     * f_n_demand(:)
          plant_uptake_no3_n15_sl(:,isoil) = plant_uptake_no3_n15_sl(:,isoil) * f_n_demand(:)
          plant_uptake_po4_sl(:,isoil)     = plant_uptake_po4_sl(:,isoil)     * f_p_demand(:)
        ENDDO
    ENDIF


    !> 2.0 infrastructure for mass balance checks
    !!
    !! fluxes are calculated at the end of this routine (to include het_respiration)
    !!   
    !! changes in pools per time step are calculated here AND at the end of this routine
    sb_delta_pools%carbon(:,:) = &
       sb_pool%dom%carbon(:,:) + sb_pool%dom_assoc%carbon(:,:) + &
       sb_pool%soluable_litter%carbon(:,:) + sb_pool%polymeric_litter%carbon(:,:) + &
       sb_pool%woody_litter%carbon(:,:) + sb_pool%fungi%carbon(:,:) + &
       sb_pool%mycorrhiza%carbon(:,:) + sb_pool%microbial%carbon(:,:) + &
       sb_pool%residue%carbon(:,:) + sb_pool%residue_assoc%carbon(:,:)
    sb_delta_pools%nitrogen(:,:) = &
       sb_pool%dom%nitrogen(:,:) + sb_pool%dom_assoc%nitrogen(:,:) + &
       sb_pool%soluable_litter%nitrogen(:,:) + sb_pool%polymeric_litter%nitrogen(:,:) + &
       sb_pool%woody_litter%nitrogen(:,:) + sb_pool%fungi%nitrogen(:,:) + &
       sb_pool%mycorrhiza%nitrogen(:,:) + sb_pool%microbial%nitrogen(:,:) + &
       sb_pool%residue%nitrogen(:,:) + sb_pool%residue_assoc%nitrogen(:,:) + &
       nh4_solute(:,:) + nh4_assoc(:,:) + no3_solute(:,:) + noy(:,:) + n2o(:,:) + n2(:,:)
    sb_delta_pools%phosphorus(:,:) = &
       sb_pool%dom%phosphorus(:,:) + sb_pool%dom_assoc%phosphorus(:,:) + &
       sb_pool%soluable_litter%phosphorus(:,:) + sb_pool%polymeric_litter%phosphorus(:,:) + &
       sb_pool%woody_litter%phosphorus(:,:) + sb_pool%fungi%phosphorus(:,:) + &
       sb_pool%mycorrhiza%phosphorus(:,:) + sb_pool%microbial%phosphorus(:,:) + &
       sb_pool%residue%phosphorus(:,:) + sb_pool%residue_assoc%phosphorus(:,:) + &
       po4_solute(:,:) + po4_occluded(:,:) + po4_assoc_slow(:,:) + po4_assoc_fast(:,:) + po4_primary(:,:)
    sb_delta_pools%carbon13(:,:) = &
       sb_pool%dom%carbon13(:,:) + sb_pool%dom_assoc%carbon13(:,:) + &
       sb_pool%soluable_litter%carbon13(:,:) + sb_pool%polymeric_litter%carbon13(:,:) + &
       sb_pool%woody_litter%carbon13(:,:) + sb_pool%fungi%carbon13(:,:) + &
       sb_pool%mycorrhiza%carbon13(:,:) + sb_pool%microbial%carbon13(:,:) + &
       sb_pool%residue%carbon13(:,:) + sb_pool%residue_assoc%carbon13(:,:)
    ! for mass balance purposes, already reduce by radioactive decay, since this is not a reported flux
    sb_delta_pools%carbon14(:,:) = &
       (sb_pool%dom%carbon14(:,:) + sb_pool%dom_assoc%carbon14(:,:) + &
       sb_pool%soluable_litter%carbon14(:,:) + sb_pool%polymeric_litter%carbon14(:,:) + &
       sb_pool%woody_litter%carbon14(:,:) + sb_pool%fungi%carbon14(:,:) + &
       sb_pool%mycorrhiza%carbon14(:,:) + sb_pool%microbial%carbon14(:,:) + &
       sb_pool%residue%carbon14(:,:) + sb_pool%residue_assoc%carbon14(:,:)) * &
       (1._wp - lambda_C14 * dtime)
    sb_delta_pools%nitrogen15(:,:) = &
       sb_pool%dom%nitrogen15(:,:) + sb_pool%dom_assoc%nitrogen15(:,:) + &
       sb_pool%soluable_litter%nitrogen15(:,:) + sb_pool%polymeric_litter%nitrogen15(:,:) + &
       sb_pool%woody_litter%nitrogen15(:,:) + sb_pool%fungi%nitrogen15(:,:) + &
       sb_pool%mycorrhiza%nitrogen15(:,:) + sb_pool%microbial%nitrogen15(:,:) + &
       sb_pool%residue%nitrogen15(:,:) + sb_pool%residue_assoc%nitrogen15(:,:) + &
       nh4_n15_solute(:,:) + nh4_n15_assoc(:,:) + no3_n15_solute(:,:) + noy_n15(:,:) + n2o_n15(:,:) + n2_n15(:,:)


    !> 3.0 loss of C14 due to radioactive decay
    !!
    sb_loss%soluable_litter%carbon14(:,:)  = sb_loss%soluable_litter%carbon14(:,:) + & 
       sb_pool%soluable_litter%carbon14(:,:)  * lambda_C14 * dtime
    sb_loss%polymeric_litter%carbon14(:,:) = sb_loss%polymeric_litter%carbon14(:,:) + &
       sb_pool%polymeric_litter%carbon14(:,:) * lambda_C14 * dtime
    sb_loss%woody_litter%carbon14(:,:)     = sb_loss%woody_litter%carbon14(:,:) + &
       sb_pool%woody_litter%carbon14(:,:)     * lambda_C14 * dtime
    sb_loss%dom%carbon14(:,:)              = sb_loss%dom%carbon14(:,:) + &
       sb_pool%dom%carbon14(:,:)              * lambda_C14 * dtime
    sb_loss%dom_assoc%carbon14(:,:)        = sb_loss%dom_assoc%carbon14(:,:) + &
       sb_pool%dom_assoc%carbon14(:,:)        * lambda_C14 * dtime
    sb_loss%fungi%carbon14(:,:)            = sb_loss%fungi%carbon14(:,:) + & 
       sb_pool%fungi%carbon14(:,:)            * lambda_C14 * dtime
    sb_loss%mycorrhiza%carbon14(:,:)       = sb_loss%mycorrhiza%carbon14(:,:) + & 
       sb_pool%mycorrhiza%carbon14(:,:)       * lambda_C14 * dtime
    sb_loss%microbial%carbon14(:,:)        = sb_loss%microbial%carbon14(:,:) + & 
       sb_pool%microbial%carbon14(:,:)        * lambda_C14 * dtime
    sb_loss%residue%carbon14(:,:)          = sb_loss%residue%carbon14(:,:) + &
       sb_pool%residue%carbon14(:,:)          * lambda_C14 * dtime
    sb_loss%residue_assoc%carbon14(:,:)    = sb_loss%residue_assoc%carbon14(:,:) + & 
       sb_pool%residue_assoc%carbon14(:,:)    * lambda_C14 * dtime 


    !> 4.0 update mycorrhiza with flux due to exudation from vegetation and N uptake from soil AND correct flux due to export to plants
    !!
    IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_0d") THEN
       f_myc_sl(:,1)          = 1.0_wp
       f_myc_sl(:,2:nsoil_sb) = 0.0_wp
    ELSE
       f_myc_sl(:,:)          = root_fraction_sl(:,:) 
    ENDIF

    DO isoil = 1,nsoil_sb
       sb_formation%mycorrhiza%carbon(:,isoil)     = sb_formation%mycorrhiza%carbon(:,isoil) + &
                                                     f_myc_sl(:,isoil) * mycorrhiza_cue * veg_exudation%carbon(:) / soil_depth_sl(:,isoil)
       sb_formation%mycorrhiza%nitrogen(:,isoil)   = sb_formation%mycorrhiza%nitrogen(:,isoil) + &
                                                     f_myc_sl(:,isoil) * veg_exudation%nitrogen(:) / soil_depth_sl(:,isoil) + &
                                                     (mycorrhiza_uptake_nh4_sl(:,isoil) + mycorrhiza_uptake_no3_sl(:,isoil)) * &
                                                     dtime / 1000000._wp
       sb_formation%mycorrhiza%phosphorus(:,isoil) = sb_formation%mycorrhiza%phosphorus(:,isoil) + &
                                                     f_myc_sl(:,isoil) * veg_exudation%phosphorus(:) / soil_depth_sl(:,isoil) + &
                                                     mycorrhiza_uptake_po4_sl(:,isoil) * dtime / 1000000._wp
       sb_formation%mycorrhiza%carbon13(:,isoil)   = sb_formation%mycorrhiza%carbon13(:,isoil) + &
                                                     f_myc_sl(:,isoil) * mycorrhiza_cue * veg_exudation%carbon13(:) / soil_depth_sl(:,isoil)
       sb_formation%mycorrhiza%carbon14(:,isoil)   = sb_formation%mycorrhiza%carbon14(:,isoil) + &
                                                     f_myc_sl(:,isoil) * mycorrhiza_cue * veg_exudation%carbon14(:) / soil_depth_sl(:,isoil)
       sb_formation%mycorrhiza%nitrogen15(:,isoil) = sb_formation%mycorrhiza%nitrogen15(:,isoil) + &
                                                     f_myc_sl(:,isoil) * veg_exudation%nitrogen15(:) / soil_depth_sl(:,isoil) + &
                                                     (mycorrhiza_uptake_nh4_n15_sl(:,isoil) + mycorrhiza_uptake_no3_n15_sl(:,isoil)) * &
                                                     dtime / 1000000._wp

       het_respiration(:,isoil)     = het_respiration(:,isoil)     + f_myc_sl(:,isoil) * (1._wp - mycorrhiza_cue) * & 
                                      veg_exudation%carbon(:) / soil_depth_sl(:,isoil) / dtime * 1000000._wp                          
       myc_respiration(:,isoil)     = myc_respiration(:,isoil)     + f_myc_sl(:,isoil) * (1._wp - mycorrhiza_cue) * &
                                      veg_exudation%carbon(:) / soil_depth_sl(:,isoil) / dtime * 1000000._wp
       het_respiration_c13(:,isoil) = het_respiration_c13(:,isoil) + f_myc_sl(:,isoil) * (1._wp - mycorrhiza_cue) * & 
                                      veg_exudation%carbon13(:) / soil_depth_sl(:,isoil) / dtime * 1000000._wp                          
       myc_respiration_c13(:,isoil) = myc_respiration_c13(:,isoil) + f_myc_sl(:,isoil) * (1._wp - mycorrhiza_cue) * &
                                      veg_exudation%carbon13(:) / soil_depth_sl(:,isoil) / dtime * 1000000._wp
       het_respiration_c14(:,isoil) = het_respiration_c14(:,isoil) + f_myc_sl(:,isoil) * (1._wp - mycorrhiza_cue) * & 
                                      veg_exudation%carbon14(:) / soil_depth_sl(:,isoil) / dtime * 1000000._wp                          
       myc_respiration_c14(:,isoil) = myc_respiration_c14(:,isoil) + f_myc_sl(:,isoil) * (1._wp - mycorrhiza_cue) * &
                                      veg_exudation%carbon14(:) / soil_depth_sl(:,isoil) / dtime * 1000000._wp
    END DO

    sb_loss%mycorrhiza%carbon(:,:)     = sb_loss%mycorrhiza%carbon(:,:)     + sb_mycorrhiza_export%carbon(:,:) 
    sb_loss%mycorrhiza%nitrogen(:,:)   = sb_loss%mycorrhiza%nitrogen(:,:)   + sb_mycorrhiza_export%nitrogen(:,:)
    sb_loss%mycorrhiza%phosphorus(:,:) = sb_loss%mycorrhiza%phosphorus(:,:) + sb_mycorrhiza_export%phosphorus(:,:)
    sb_loss%mycorrhiza%carbon13(:,:)   = sb_loss%mycorrhiza%carbon13(:,:)   + sb_mycorrhiza_export%carbon13(:,:)
    sb_loss%mycorrhiza%carbon14(:,:)   = sb_loss%mycorrhiza%carbon14(:,:)   + sb_mycorrhiza_export%carbon14(:,:)
    sb_loss%mycorrhiza%nitrogen15(:,:) = sb_loss%mycorrhiza%nitrogen15(:,:) + sb_mycorrhiza_export%nitrogen15(:,:)


    !> 5.0 update of the organic pool: simple_0d / simple_1d
    !!
    !! update all components of the sb_pool applying the JSM equation for formation, loss and transport
    !! 
    SELECT CASE (TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
      CASE ("simple_0d", "simple_1d")
        sb_pool = sb_pool + sb_formation - sb_loss + sb_transport 
    END SELECT

    !> 6.0 JSM: 
    !!  a) update of the organic pool
    !!  b) update particle transport related diagnostics and state variables
    !!
    !! docu: 
    !!  a) sb_last_timestep pool is used only for calc_particle_fluxrate() to calc (pseudo code):
    !!     particle_fluxrate = sb_pool - sb_last_timestep  (i.e., current timestep - previous timestep)
    !!  b) par_transport pool is used only as hlp to collect the output of calc_particle_transport() 
    !!     within calc_particle_transport_wrapper() and apply these temporary values to sb_transport & sb_pool
    IF (TRIM(dsl4jsb_Config(SB_)%sb_model_scheme)=="jsm") THEN
      !> 6.1 update all components of the sb_pool applying the JSM equation for formation, loss and transport
      !!
      sb_last_timestep = sb_pool
      sb_pool          = sb_pool + sb_formation - sb_loss + sb_transport 
      !> 6.2 particle transport
      !!
      par_transport          = sb_transport * 0._wp
      particle_fluxrate(:,:) = calc_particle_fluxrate(nc, nsoil_sb, & 
                                  soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), & 
                                  sb_pool, sb_last_timestep) 
      CALL calc_particle_transport_wrapper(nc, nsoil_sb, particle_fluxrate(:,:), soil_depth_sl(:,:), &                        ! in
                                  sb_pool, nh4_solute(:,:), nh4_assoc(:,:), no3_solute(:,:), po4_solute(:,:), & 
                                  po4_assoc_fast(:,:), po4_assoc_slow(:,:), po4_occluded(:,:), po4_primary(:,:), & 
                                  nh4_n15_solute(:,:), nh4_n15_assoc(:,:), no3_n15_solute(:,:), & 
                                  par_transport, transport_nh4_solute(:,:), transport_nh4_assoc(:,:), &                       ! inout
                                  transport_no3_solute(:,:), transport_po4_solute(:,:), & 
                                  transport_po4_assoc_fast(:,:), transport_po4_assoc_slow(:,:), transport_po4_occluded(:,:), &
                                  transport_po4_primary(:,:), & 
                                  transport_nh4_n15_solute(:,:), transport_nh4_n15_assoc(:,:), transport_no3_n15_solute(:,:))
      sb_transport = sb_transport + par_transport
      sb_pool      = sb_pool      + par_transport 
      dom_cn(:,:)  = sb_pool%dom%carbon(:,:) / sb_pool%dom%nitrogen(:,:)
      dom_cp(:,:)  = sb_pool%dom%carbon(:,:) / sb_pool%dom%phosphorus(:,:)
      !! convert inorganic transport fluxes to µmol m3 s-1 (for compatibility with simple_0d & simple_1d)
      !!
      transport_nh4_solute(:,:)     = transport_nh4_solute(:,:)     * 1.e6_wp / dtime
      transport_nh4_assoc(:,:)      = transport_nh4_assoc(:,:)      * 1.e6_wp / dtime
      transport_no3_solute(:,:)     = transport_no3_solute(:,:)     * 1.e6_wp / dtime
      transport_po4_solute(:,:)     = transport_po4_solute(:,:)     * 1.e6_wp / dtime
      transport_po4_assoc_fast(:,:) = transport_po4_assoc_fast(:,:) * 1.e6_wp / dtime
      transport_po4_assoc_slow(:,:) = transport_po4_assoc_slow(:,:) * 1.e6_wp / dtime
      transport_po4_occluded(:,:)   = transport_po4_occluded(:,:)   * 1.e6_wp / dtime
      transport_po4_primary(:,:)    = transport_po4_primary(:,:)    * 1.e6_wp / dtime
      transport_nh4_n15_solute(:,:) = transport_nh4_n15_solute(:,:) * 1.e6_wp / dtime
      transport_nh4_n15_assoc(:,:)  = transport_nh4_n15_assoc(:,:)  * 1.e6_wp / dtime
      transport_no3_n15_solute(:,:) = transport_no3_n15_solute(:,:) * 1.e6_wp / dtime
    ENDIF

    !> 7.0 update of inorganic pools
    !!
    !> 7.1 update of N pools  
    !!
    no3_solute(:,:)     = no3_solute(:,:)     + &
                          (net_mineralisation_no3(:,:) + nitrification_no3(:,:) + transport_no3_solute(:,:) - & 
                           plant_uptake_no3_sl(:,:) - mycorrhiza_uptake_no3_sl(:,:) - &
                           denitrification_noy(:,:) - denitrification_n2o(:,:) - denitrification_n2(:,:)) * dtime/1.e6_wp
    no3_n15_solute(:,:) = no3_n15_solute(:,:) + &
                          (net_mineralisation_no3_n15(:,:) + nitrification_no3_n15(:,:) + transport_no3_n15_solute(:,:) - & 
                           plant_uptake_no3_n15_sl(:,:) - mycorrhiza_uptake_no3_n15_sl(:,:) + &
                           denitrification_noy_n15(:,:) - denitrification_n2o_n15(:,:) - denitrification_n2_n15(:,:)) * dtime/1.e6_wp


    noy(:,:)            = noy(:,:)            + (nitrification_noy(:,:) + denitrification_noy(:,:) + & 
                                                 transport_noy(:,:)) * dtime/1.e6_wp
    noy_n15(:,:)        = noy_n15(:,:)        + (nitrification_noy_n15(:,:) + denitrification_noy_n15(:,:) + & 
                                                 transport_noy_n15(:,:)) * dtime/1.e6_wp
    n2o(:,:)            = n2o(:,:)            + (nitrification_n2o(:,:) + denitrification_n2o(:,:) + & 
                                                 transport_n2o(:,:)) * dtime/1.e6_wp
    n2o_n15(:,:)        = n2o_n15(:,:)        + (nitrification_n2o_n15(:,:) + denitrification_n2o_n15(:,:) + & 
                                                 transport_n2o_n15(:,:)) * dtime/1.e6_wp
    n2(:,:)             = n2(:,:)             + (denitrification_n2(:,:) + transport_n2(:,:)) * dtime/1.e6_wp
    n2_n15(:,:)         = n2_n15(:,:)         + (denitrification_n2_n15(:,:) + transport_n2_n15(:,:)) * dtime/1.e6_wp


    !> 8.0 update diagnostics and other state variables depending on organic matter content 
    !! OC
    CALL calc_bulk_soil_carbon(nc, nsoil_sb, sb_pool, bulk_soil_carbon_sl, soil_litter_carbon_sl, &
                               volume_min_sl, TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
    bulk_dens_corr_sl(:,:) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:), soil_litter_carbon_sl(:,:), bulk_dens_sl(:,:))
    qmax_org_min_sl(:,:)   = calc_qmax_texture(dsl4jsb_Config(SB_)%qmax_org_fine_particle, silt_sl(:,:), clay_sl(:,:))
    CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), bulk_dens_sl(:,:), &
                                           0._wp, qmax_org_min_sl(:,:), qmax_org(:,:))
    !! po4                         
    CALL calc_Psorption_parameter(nc, nsoil_sb, dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
                                  clay_sl(:,:), silt_sl(:,:), rmm_sorption(:,:) * 1000._wp, ph_sl(:,:), &
                                  bulk_soil_carbon_sl(:,:), soil_litter_carbon_sl(:,:), sb_pool%woody_litter%carbon(:,:), &
                                  bulk_dens_sl(:,:), volume_min_sl(:,:), po4_solute(:,:), &
                                  qmax_po4_min_sl(:,:), qmax_po4_om_sl(:,:), qmax_po4(:,:), km_adsorpt_po4_sl(:,:), &
                                  qmax_fast_po4(:,:), qmax_slow_po4(:,:), km_fast_po4(:,:), km_slow_po4(:,:), &
                                  Qmax_AlFe_cor(:,:))
    !! nh4                              
    CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), bulk_dens_sl(:,:), &
                                           0.0_wp,qmax_nh4_min_sl(:,:), qmax_nh4(:,:))
    qmax_nh4(:,:)          = qmax_nh4(:,:) * rmm_sorption(:,:)
    CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), bulk_dens_sl(:,:), &
                                           km_adsorpt_OM_nh4,km_adsorpt_mineral_nh4, km_adsorpt_nh4_sl(:,:))

    !> 9.0 atmospheric deposition
    !!

    no3_solute(:,1)      = no3_solute(:,1)     + noy_deposition(:) / 1.e6_wp / soil_depth_sl(:,1) * dtime
    no3_n15_solute(:,1)  = no3_n15_solute(:,1) + noy_n15_deposition(:) / 1.e6_wp / soil_depth_sl(:,1) * dtime

    !> 10.0 update of mineral NH4 & P pools
    !!
    p_deposition_sl(:,:) = 0._wp
    p_deposition_sl(:,1) = p_deposition(:) / soil_depth_sl(:,1) ! µmol P / m3 / s
    !> 10.1 init vertical variables
    !! 
    IF (TRIM(dsl4jsb_Config(SB_)%sb_model_scheme)=="simple_0d") THEN
      rtm_sorption_act(:,:)     = 0.0_wp 
      rmm_sorption_act(:,:)     = 0.0_wp
      rtm_desorption_act(:,:)   = 0.0_wp 
      rmm_desorption_act(:,:)   = 0.0_wp
      k_adsorpt_po4_act(:,:)    = 0.0_wp 
      k_desorpt_po4_act(:,:)    = 0.0_wp
      km_adsorpt_po4_act(:,:)   = 0.0_wp
      qmax_po4_mineral_act(:,:) = 0.0_wp
      km_fast_po4_act(:,:)      = 0.0_wp
      qmax_fast_po4_act(:,:)    = 0.0_wp
      km_slow_po4_act(:,:)      = 0.0_wp
      qmax_slow_po4_act(:,:)    = 0.0_wp
      km_adsorpt_nh4_act(:,:)   = 0.0_wp
      qmax_nh4_mineral_act(:,:) = 0.0_wp

      DO isoil = 1,nsoil_sb
         rtm_sorption_act(:,1)   = rtm_sorption_act(:,1)   + rtm_sorption(:,isoil)      * root_fraction_sl(:,isoil) 
         rmm_sorption_act(:,1)   = rmm_sorption_act(:,1)   + rmm_sorption(:,isoil)      * root_fraction_sl(:,isoil) 
         rtm_desorption_act(:,1) = rtm_desorption_act(:,1) + rtm_desorption(:,isoil)    * root_fraction_sl(:,isoil) 
         rmm_desorption_act(:,1) = rmm_desorption_act(:,1) + rmm_desorption(:,isoil)    * root_fraction_sl(:,isoil) 
         km_adsorpt_po4_act(:,1) = km_adsorpt_po4_act(:,1) + km_adsorpt_po4_sl(:,isoil) * root_fraction_sl(:,isoil) 
         km_adsorpt_nh4_act(:,1) = km_adsorpt_nh4_act(:,1) + km_adsorpt_nh4_sl(:,isoil) * root_fraction_sl(:,isoil) 
         km_fast_po4_act(:,1)    = km_fast_po4_act(:,1)    + km_fast_po4(:,isoil)       * root_fraction_sl(:,isoil) 
         km_slow_po4_act(:,1)    = km_slow_po4_act(:,1)    + km_slow_po4(:,isoil)       * root_fraction_sl(:,isoil)          
      ENDDO  

      root_depth_avg(:)         = SUM(root_fraction_sl(:,:) * soil_depth_sl(:,:), DIM=2)
      qmax_nh4_mineral_act(:,:) = qmax_nh4(:,:)
      CALL convert_s1d_to_s0d_var_soilproperty (nc, nsoil_sb, soil_depth_sl(:,:), qmax_nh4_mineral_act(:,:))
      km_adsorpt_nh4_act(:,1)   = km_adsorpt_nh4_act(:,1) / root_depth_avg(:)
      qmax_po4_mineral_act(:,:) = qmax_po4(:,:)
      CALL convert_s1d_to_s0d_var_soilproperty (nc, nsoil_sb, soil_depth_sl(:,:), qmax_po4_mineral_act(:,:))
      qmax_fast_po4_act(:,:)    = qmax_fast_po4(:,:)
      CALL convert_s1d_to_s0d_var_soilproperty (nc, nsoil_sb, soil_depth_sl(:,:), qmax_fast_po4_act(:,:))
      qmax_slow_po4_act(:,:)    = qmax_slow_po4(:,:)
      CALL convert_s1d_to_s0d_var_soilproperty (nc, nsoil_sb, soil_depth_sl(:,:), qmax_slow_po4_act(:,:))
      km_adsorpt_po4_act(:,1)   = km_adsorpt_po4_act(:,1) / root_depth_avg(:)
      km_fast_po4_act(:,1)      = km_fast_po4_act(:,1) / root_depth_avg(:)
      km_slow_po4_act(:,1)      = km_slow_po4_act(:,1) / root_depth_avg(:)
      k_adsorpt_po4_act(:,1)    = rtm_sorption_act(:,1)   * rmm_sorption_act(:,1)   * k_adsorpt_po4
      k_desorpt_po4_act(:,1)    = rtm_desorption_act(:,1) * rmm_desorption_act(:,1) * k_desorpt_po4
    ELSE
      k_adsorpt_po4_act(:,:)    = rtm_sorption(:,:)   * rmm_sorption(:,:)   * k_adsorpt_po4
      k_desorpt_po4_act(:,:)    = rtm_desorption(:,:) * rmm_desorption(:,:) * k_desorpt_po4
      km_adsorpt_po4_act(:,:)   = km_adsorpt_po4_sl(:,:)
      qmax_po4_mineral_act(:,:) = qmax_po4(:,:)
      km_fast_po4_act(:,:)      = km_fast_po4(:,:)
      qmax_fast_po4_act(:,:)    = qmax_fast_po4(:,:)
      km_slow_po4_act(:,:)      = km_slow_po4(:,:)
      qmax_slow_po4_act(:,:)    = qmax_slow_po4(:,:)
      km_adsorpt_nh4_act(:,:)   = km_adsorpt_nh4_sl(:,:)
      qmax_nh4_mineral_act(:,:) = qmax_nh4(:,:)
    ENDIF

    !> 10.2 update NH4 pools if not prescribed
    !! 
    IF(.NOT. dsl4jsb_Config(SB_)%flag_sb_prescribe_nh4) THEN
      hlp1(:,:) = (net_mineralisation_nh4(:,:) + asymb_n_fixation(:,:) + transport_nh4_solute(:,:) - & 
                   plant_uptake_nh4_sl(:,:) - mycorrhiza_uptake_nh4_sl(:,:) - & 
                   nitrification_no3(:,:) - nitrification_noy(:,:) - nitrification_n2o(:,:) - &
                   volatilisation_nh4(:,:) + transport_nh4_assoc(:,:)) * dtime / 1.e6_wp 
      hlp1(:,1) = hlp1(:,1) + nhx_deposition(:) / 1.e6_wp / soil_depth_sl(:,1) * dtime
      hlp2(:,:) = -transport_nh4_assoc(:,:) * dtime / 1.e6_wp
      
      IF (TRIM(dsl4jsb_Config(SB_)%sb_model_scheme)=="simple_0d") THEN
        CALL calc_langmuir_kinetics(hlp1(:,1), nh4_solute(:,1),nh4_assoc(:,1), & 
                                                     km_adsorpt_nh4_act(:,1), qmax_nh4_mineral_act(:,1), &
                                                     partition_coef(:,1))
      ELSE
        CALL calc_langmuir_kinetics(hlp1(:,:), nh4_solute(:,:),nh4_assoc(:,:), & 
                                                     km_adsorpt_nh4_act(:,:), qmax_nh4_mineral_act(:,:), &
                                                     partition_coef(:,:))
      ENDIF
      
      nh4_solute(:,:)        = nh4_solute(:,:) + partition_coef(:,:) * hlp1(:,:)
      nh4_assoc(:,:)         = nh4_assoc(:,:)  + (1._wp - partition_coef(:,:)) * hlp1(:,:) 
      fast_exchange_nh4(:,:) = (partition_coef(:,:) * hlp1(:,:)+hlp2(:,:)) * 1.e6_wp / dtime
      hlp2(:,:) = (net_mineralisation_nh4_n15(:,:) + asymb_n_fixation_n15(:,:) + transport_nh4_n15_solute(:,:) - & 
                   plant_uptake_nh4_n15_sl(:,:) - mycorrhiza_uptake_nh4_n15_sl(:,:) - &
                   nitrification_no3_n15(:,:) - nitrification_noy_n15(:,:) - nitrification_n2o_n15(:,:) - &
                   volatilisation_nh4_n15(:,:) + transport_nh4_n15_assoc(:,:)) * dtime / 1.e6_wp 
      hlp2(:,1) = hlp2(:,1) + nhx_n15_deposition(:) / 1.e6_wp / soil_depth_sl(:,1) * dtime

      WHERE ((nh4_n15_solute(:,:) + nh4_n15_assoc(:,:)) < -hlp2(:,:))
        hlp2(:,:) = -(nh4_n15_solute(:,:) + nh4_n15_assoc(:,:))
      ENDWHERE 
      WHERE (nh4_solute(:,:) + nh4_assoc(:,:) > eps8)
        partition_coef(:,:) = nh4_solute(:,:) / (nh4_solute(:,:) + nh4_assoc(:,:))
      ENDWHERE

      nh4_n15_solute(:,:) = partition_coef(:,:) * (nh4_n15_solute(:,:) + nh4_n15_assoc(:,:) + hlp2(:,:))
      nh4_n15_assoc(:,:)  = (1._wp - partition_coef(:,:)) * (nh4_n15_solute(:,:) + nh4_n15_assoc(:,:) + hlp2(:,:))
    ENDIF

    !> 10.3 update mineral P pools if not prescibed
    !! 
    IF(.NOT. dsl4jsb_Config(SB_)%flag_sb_prescribe_po4) THEN
      CALL calc_p_inorg_dynamics(nc, nsoil_sb, &
                                       TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
                                       TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &
                                       dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
                                       po4_solute(:,:), &
                                       po4_assoc_fast(:,:), &
                                       po4_assoc_slow(:,:), &
                                       qmax_po4_mineral_act(:,:), &
                                       km_adsorpt_po4_act(:,:), &
                                       qmax_fast_po4_act(:,:), &
                                       km_fast_po4_act(:,:), &
                                       qmax_slow_po4_act(:,:), &
                                       km_slow_po4_act(:,:), &
                                       k_adsorpt_po4_act(:,:), &
                                       k_desorpt_po4_act(:,:), &
                                       p_deposition_sl(:,:), &
                                       biochem_mineralisation_po4(:,:), &
                                       mycorrhiza_uptake_po4_sl(:,:), &
                                       plant_uptake_po4_sl(:,:), &
                                       microbial_uptake_po4_sl(:,:), &
                                       gross_mineralisation_po4(:,:), &
                                       net_mineralisation_po4(:,:), &
                                       transport_po4_solute(:,:), &
                                       transport_po4_assoc_fast(:,:), &
                                       transport_po4_assoc_slow(:,:), &
                                       weathering_po4(:,:), &
                                       occlusion_po4(:,:), &
                                       slow_exchange_po4(:,:), &
                                       fast_exchange_po4(:,:))

      !! not update the primary and the occlusion pool during transient-forcing spinup (i.e., before the actual transient simulation)
      IF (model%config%l_transient_spinup) THEN
        po4_occluded(:,:)   = po4_occluded(:,:)   + transport_po4_occluded(:,:) * dtime / 1.e6_wp
        po4_primary(:,:)    = po4_primary(:,:)    + transport_po4_primary(:,:)  * dtime / 1.e6_wp
      ELSE
        po4_occluded(:,:)   = po4_occluded(:,:)   + (occlusion_po4(:,:)   + transport_po4_occluded(:,:)) * dtime / 1.e6_wp
        po4_primary(:,:)    = po4_primary(:,:)    + (-weathering_po4(:,:) + transport_po4_primary(:,:))  * dtime / 1.e6_wp
      ENDIF
    ENDIF
   
    !> 11.0 Ensure that mineral nutrients are constant when prescibed
    !!
    !> 11.1 NH4
    !!
    IF(dsl4jsb_Config(SB_)%flag_sb_prescribe_nh4)THEN
      IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_0d") THEN
        hlp1(:,:) = 0.0_wp
        hlp1(:,1) = nh4_solute_prescribe / soil_depth_sl(:,1) - nh4_solute(:,1)
      ELSE 
        hlp1(:,:) = nh4_solute_prescribe - nh4_solute(:,:)
      END IF
      nh4_solute(:,:)     = nh4_solute(:,:) + hlp1(:,:) 
      nh4_n15_solute(:,:) = nh4_n15_solute(:,:) + hlp1(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
    ENDIF

    !> 11.2 NO3
    !!
    IF(dsl4jsb_Config(SB_)%flag_sb_prescribe_no3)THEN
      IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_0d") THEN
        hlp1(:,:) = 0.0_wp
        hlp1(:,1) = no3_solute_prescribe / soil_depth_sl(:,1) - no3_solute(:,1)
      ELSE
        hlp1(:,:) = no3_solute_prescribe - no3_solute(:,:)
      END IF
      no3_solute(:,:)     = no3_solute(:,:) + hlp1(:,:) 
      no3_n15_solute(:,:) = no3_n15_solute(:,:) + hlp1(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
    ENDIF
    
    !> 11.3 PO4
    !!
    IF(dsl4jsb_Config(SB_)%flag_sb_prescribe_po4)THEN
      IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_0d") THEN
        po4_solute(:,:) = 0.0_wp
        po4_solute(:,1) = po4_solute_prescribe / soil_depth_sl(:,1)
      ELSE
        po4_solute(:,:) = po4_solute_prescribe 
      END IF
    ENDIF


    !> 12.0 infrastructure for mass balance checks (see also section 0.9)
    !!
    !! pools: calc the difference to what was calculated at the beginning of this routine
    sb_delta_pools%carbon(:,:) = &
       (sb_pool%dom%carbon(:,:) + sb_pool%dom_assoc%carbon(:,:) + &
       sb_pool%soluable_litter%carbon(:,:) + sb_pool%polymeric_litter%carbon(:,:) + &
       sb_pool%woody_litter%carbon(:,:) + sb_pool%fungi%carbon(:,:) + &
       sb_pool%mycorrhiza%carbon(:,:) + sb_pool%microbial%carbon(:,:) + &
       sb_pool%residue%carbon(:,:) + sb_pool%residue_assoc%carbon(:,:)) - & 
       sb_delta_pools%carbon(:,:)
    sb_delta_pools%nitrogen(:,:) = &
       (sb_pool%dom%nitrogen(:,:) + sb_pool%dom_assoc%nitrogen(:,:) + &
       sb_pool%soluable_litter%nitrogen(:,:) + sb_pool%polymeric_litter%nitrogen(:,:) + &
       sb_pool%woody_litter%nitrogen(:,:) + sb_pool%fungi%nitrogen(:,:) + &
       sb_pool%mycorrhiza%nitrogen(:,:) + sb_pool%microbial%nitrogen(:,:) + &
       sb_pool%residue%nitrogen(:,:) + sb_pool%residue_assoc%nitrogen(:,:) + &
       nh4_solute(:,:) + nh4_assoc(:,:) + no3_solute(:,:) + noy(:,:) + n2o(:,:) + n2(:,:)) - &
       sb_delta_pools%nitrogen(:,:)
    sb_delta_pools%phosphorus(:,:) = &
       (sb_pool%dom%phosphorus(:,:) + sb_pool%dom_assoc%phosphorus(:,:) + &
       sb_pool%soluable_litter%phosphorus(:,:) + sb_pool%polymeric_litter%phosphorus(:,:) + &
       sb_pool%woody_litter%phosphorus(:,:) + sb_pool%fungi%phosphorus(:,:) + &
       sb_pool%mycorrhiza%phosphorus(:,:) + sb_pool%microbial%phosphorus(:,:) + &
       sb_pool%residue%phosphorus(:,:) + sb_pool%residue_assoc%phosphorus(:,:) + &
       po4_solute(:,:) + po4_occluded(:,:) + po4_assoc_slow(:,:) + po4_assoc_fast(:,:) + po4_primary(:,:)) - &
       sb_delta_pools%phosphorus(:,:)
    sb_delta_pools%carbon13(:,:) = &
       (sb_pool%dom%carbon13(:,:) + sb_pool%dom_assoc%carbon13(:,:) + &
       sb_pool%soluable_litter%carbon13(:,:) + sb_pool%polymeric_litter%carbon13(:,:) + &
       sb_pool%woody_litter%carbon13(:,:) + sb_pool%fungi%carbon13(:,:) + &
       sb_pool%mycorrhiza%carbon13(:,:) + sb_pool%microbial%carbon13(:,:) + &
       sb_pool%residue%carbon13(:,:) + sb_pool%residue_assoc%carbon13(:,:)) - &
       sb_delta_pools%carbon13(:,:)
    sb_delta_pools%carbon14(:,:) = &
       (sb_pool%dom%carbon14(:,:) + sb_pool%dom_assoc%carbon14(:,:) + &
       sb_pool%soluable_litter%carbon14(:,:) + sb_pool%polymeric_litter%carbon14(:,:) + &
       sb_pool%woody_litter%carbon14(:,:) + sb_pool%fungi%carbon14(:,:) + &
       sb_pool%mycorrhiza%carbon14(:,:) + sb_pool%microbial%carbon14(:,:) + &
       sb_pool%residue%carbon14(:,:) + sb_pool%residue_assoc%carbon14(:,:)) - &
       sb_delta_pools%carbon14(:,:)
    sb_delta_pools%nitrogen15(:,:) = &
       (sb_pool%dom%nitrogen15(:,:) + sb_pool%dom_assoc%nitrogen15(:,:) + &
       sb_pool%soluable_litter%nitrogen15(:,:) + sb_pool%polymeric_litter%nitrogen15(:,:) + &
       sb_pool%woody_litter%nitrogen15(:,:) + sb_pool%fungi%nitrogen15(:,:) + &
       sb_pool%mycorrhiza%nitrogen15(:,:) + sb_pool%microbial%nitrogen15(:,:) + &
       sb_pool%residue%nitrogen15(:,:) + sb_pool%residue_assoc%nitrogen15(:,:) + &
       nh4_n15_solute(:,:) + nh4_n15_assoc(:,:) + no3_n15_solute(:,:) + noy_n15(:,:) + n2o_n15(:,:) + n2_n15(:,:)) - &
       sb_delta_pools%nitrogen15(:,:)

    !! fluxes
    total_flux_carbon(:) =                  veg_litterfall%leaf%carbon(:)        + &
      veg_litterfall%fine_root%carbon(:) + veg_litterfall%coarse_root%carbon(:) + &
      veg_litterfall%sap_wood%carbon(:)  + veg_litterfall%heart_wood%carbon(:)  + &
      veg_litterfall%labile%carbon(:)    + veg_litterfall%reserve%carbon(:)     + &
      veg_litterfall%fruit%carbon(:)     + veg_litterfall%seed_bed%carbon(:)    + &
      veg_exudation%carbon(:)                                                    - &
      SUM(sb_mycorrhiza_export%carbon(:,:) * soil_depth_sl(:,:), DIM=2)  - &
      (leaching_dom_carbon(:)                                            + &
      SUM(het_respiration(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime / 1.e6_wp   ! DIM=2: sum only over dimension 2, i.e., soil layers
    total_flux_nitrogen(:) =                  veg_litterfall%leaf%nitrogen(:)    + &
      veg_litterfall%fine_root%nitrogen(:) + veg_litterfall%coarse_root%nitrogen(:) + &
      veg_litterfall%sap_wood%nitrogen(:)  + veg_litterfall%heart_wood%nitrogen(:)  + &
      veg_litterfall%labile%nitrogen(:)    + veg_litterfall%reserve%nitrogen(:)     + &
      veg_litterfall%fruit%nitrogen(:)     + veg_litterfall%seed_bed%nitrogen(:)    + &
      SUM(asymb_n_fixation(:,:) * soil_depth_sl(:,:), DIM=2) *dtime/1.e6_wp   + &
      (nhx_deposition(:) + noy_deposition(:)) * dtime/1.e6_wp                 + & 
      veg_exudation%nitrogen(:)                                                  - &
      SUM(sb_mycorrhiza_export%nitrogen(:,:) * soil_depth_sl(:,:), DIM=2)     - &
      (leaching_nh4_solute(:) + leaching_no3_solute(:) + leaching_dom_nitrogen(:) + &
       emission_nh3(:) + emission_noy(:) + emission_n2o(:) + emission_n2(:)   + &
       SUM(plant_uptake_nh4_sl(:,:) * soil_depth_sl(:,:), DIM=2)              + &
       SUM(plant_uptake_no3_sl(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime/1.e6_wp                       
    total_flux_phosphorus(:) =                  veg_litterfall%leaf%phosphorus(:)    + &
      veg_litterfall%fine_root%phosphorus(:) + veg_litterfall%coarse_root%phosphorus(:) + &
      veg_litterfall%sap_wood%phosphorus(:)  + veg_litterfall%heart_wood%phosphorus(:)  + &
      veg_litterfall%labile%phosphorus(:)    + veg_litterfall%reserve%phosphorus(:)     + &
      veg_litterfall%fruit%phosphorus(:)     + veg_litterfall%seed_bed%phosphorus(:)    + &
      p_deposition(:) * dtime/1.e6_wp                                             + &
      veg_exudation%phosphorus(:)                                                    - &
      SUM(sb_mycorrhiza_export%phosphorus(:,:) * soil_depth_sl(:,:), DIM=2)       - &
      (leaching_po4_solute(:) + leaching_dom_phosphorus(:)                        + &
       SUM(plant_uptake_po4_sl(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime/1.e6_wp       
    total_flux_carbon13(:) =                  veg_litterfall%leaf%carbon13(:)        + &
      veg_litterfall%fine_root%carbon13(:) + veg_litterfall%coarse_root%carbon13(:) + &
      veg_litterfall%sap_wood%carbon13(:)  + veg_litterfall%heart_wood%carbon13(:)  + &
      veg_litterfall%labile%carbon13(:)    + veg_litterfall%reserve%carbon13(:)     + &
      veg_litterfall%fruit%carbon13(:)     + veg_litterfall%seed_bed%carbon13(:)    + &
      veg_exudation%carbon13(:)                                                      - &
      SUM(sb_mycorrhiza_export%carbon13(:,:) * soil_depth_sl(:,:), DIM=2)    - &
      (leaching_dom_carbon13(:)                                              + &
      SUM(het_respiration_c13(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime/1.e6_wp
    total_flux_carbon14(:) =                  veg_litterfall%leaf%carbon14(:)        + &
      veg_litterfall%fine_root%carbon14(:) + veg_litterfall%coarse_root%carbon14(:) + &
      veg_litterfall%sap_wood%carbon14(:)  + veg_litterfall%heart_wood%carbon14(:)  + &
      veg_litterfall%labile%carbon14(:)    + veg_litterfall%reserve%carbon14(:)     + &
      veg_litterfall%fruit%carbon14(:)     + veg_litterfall%seed_bed%carbon14(:)    + &
      veg_exudation%carbon14(:)                                                      - &
      SUM(sb_mycorrhiza_export%carbon14(:,:) * soil_depth_sl(:,:), DIM=2)    - &
      (leaching_dom_carbon14(:)                                              + &
      SUM(het_respiration_c14(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime/1.e6_wp
    total_flux_nitrogen15(:) =                  veg_litterfall%leaf%nitrogen15(:)      + &
      veg_litterfall%fine_root%nitrogen15(:) + veg_litterfall%coarse_root%nitrogen15(:)   + &
      veg_litterfall%sap_wood%nitrogen15(:)  + veg_litterfall%heart_wood%nitrogen15(:)    + &
      veg_litterfall%labile%nitrogen15(:)    + veg_litterfall%reserve%nitrogen15(:)       + &
      veg_litterfall%fruit%nitrogen15(:)     + veg_litterfall%seed_bed%nitrogen15(:)      + &
      SUM(asymb_n_fixation_n15(:,:) * soil_depth_sl(:,:), DIM=2) *dtime/1.e6_wp     + &
      (nhx_deposition(:) + noy_deposition(:)) /                                       & 
         ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) *dtime/1.e6_wp         + &
      veg_exudation%nitrogen15(:)                                                   - &
      SUM(sb_mycorrhiza_export%nitrogen15(:,:) * soil_depth_sl(:,:), DIM=2)         - &
      (leaching_nh4_n15_solute(:) + leaching_no3_n15_solute(:)                      + &
       emission_nh3_n15(:) + emission_noy_n15(:) + emission_n2o_n15(:)              + & 
       emission_n2_n15(:) + leaching_dom_nitrogen15(:)                              + &
       SUM(plant_uptake_nh4_n15_sl(:,:) * soil_depth_sl(:,:), DIM=2)                + &
       SUM(plant_uptake_no3_n15_sl(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime/1.e6_wp  


  ! Deallocate local allocatable variables
  DEALLOCATE(p_deposition_sl)
  DEALLOCATE(f_myc_sl)
  DEALLOCATE(hlp1)
  DEALLOCATE(hlp2)
  DEALLOCATE(hlp3)
  DEALLOCATE(km_adsorpt_po4_act)
  DEALLOCATE(qmax_po4_mineral_act)
  DEALLOCATE(km_fast_po4_act)
  DEALLOCATE(qmax_fast_po4_act)
  DEALLOCATE(km_slow_po4_act)
  DEALLOCATE(qmax_slow_po4_act)
  DEALLOCATE(k_desorpt_po4_act)
  DEALLOCATE(k_adsorpt_po4_act)
  DEALLOCATE(rtm_sorption_act)
  DEALLOCATE(rmm_sorption_act)
  DEALLOCATE(rtm_desorption_act)
  DEALLOCATE(rmm_desorption_act)
  DEALLOCATE(partition_coef)
  DEALLOCATE(root_depth_avg)

! Pass back values to process memory
dsl4slm_passback_vegpool_var(pool,fine_root,carbon)
#define Passback_veg_litterfall
#define Passback_veg_exudation
#define Passback_sb_pool
#define Passback_sb_formation
#define Passback_sb_loss
#define Passback_sb_transport
#define Passback_sb_delta_pools
#define Passback_sb_mycorrhiza_export
#include "slm_pools_a_incl.f90"
#undef Passback_veg_litterfall
#undef Passback_veg_exudation
#undef Passback_sb_pool
#undef Passback_sb_formation
#undef Passback_sb_loss
#undef Passback_sb_transport
#undef Passback_sb_delta_pools
#undef Passback_sb_mycorrhiza_export

  END SUBROUTINE update_sb_pools



  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to calculate fast_exchange_po4 based on Langmuir equilibrium (Yang et al. 2014)
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_p_inorg_dynamics(nc, nsoil_sb, &
                                   sb_model_scheme, &
                                   sb_adsorp_scheme, &
                                   flag_sb_double_langmuir, &
                                   po4_solute, &
                                   po4_assoc_fast, &
                                   po4_assoc_slow, &
                                   qmax_po4_mineral_act, &
                                   km_adsorpt_po4_act, &
                                   qmax_fast_po4_act, &
                                   km_fast_po4_act, &
                                   qmax_slow_po4_act, &
                                   km_slow_po4_act, &
                                   k_adsorpt_po4_act, &
                                   k_desorpt_po4_act, &
                                   p_deposition_sl, &
                                   biochem_mineralisation_po4, &
                                   mycorrhiza_uptake_po4_sl, &
                                   plant_uptake_po4_sl, &
                                   microbial_uptake_po4_sl, &
                                   gross_mineralisation_po4, &
                                   net_mineralisation_po4, &
                                   transport_po4_solute, &
                                   transport_po4_assoc_fast, &
                                   transport_po4_assoc_slow, &
                                   weathering_po4, &
                                   occlusion_po4, &
                                   slow_exchange_po4, &
                                   fast_exchange_po4)

    USE mo_sb_constants
    USE mo_sb_jsm_processes,    ONLY: calc_fast_po4
    USE mo_jsb_math_constants,  ONLY: eps4, eps8, dtime

    IMPLICIT NONE

    INTEGER,                                          INTENT(in)    :: nc, &                         !< dimensions
                                                                       nsoil_sb
    CHARACTER(len=*),                                 INTENT(in)    :: sb_model_scheme               !< from sb config
    CHARACTER(len=*),                                 INTENT(in)    :: sb_adsorp_scheme              !< eca_full or eca_none
    LOGICAL ,                                         INTENT(in)    :: flag_sb_double_langmuir       !< T: double Langmuir; F: traditional Langmuir
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: qmax_po4_mineral_act, &       !< maximum sorption capacity of po4, [mol P/ m3]
                                                                       km_adsorpt_po4_act, &         !< half saturation of po4 sorption capacity, [mol P/ m3]
                                                                       qmax_fast_po4_act, &          !< maximum sorption capacity of fast_assoc_po4, [mol P/ m3]
                                                                       km_fast_po4_act, &            !< half saturation of po4 sorption in fast_assoc_po4, [mol P/ m3]
                                                                       qmax_slow_po4_act, &          !< maximum sorption capacity of slow_assoc_po4, [mol P/ m3]
                                                                       km_slow_po4_act, &            !< half saturation of po4 sorption in slow_assoc_po4, [mol P/ m3]
                                                                       k_adsorpt_po4_act, &          !< slow adsorption rate of po4, [1 / s]
                                                                       k_desorpt_po4_act, &          !< slow desorption rate of po4, [1 / s]
                                                                       p_deposition_sl, &            !< actual P deposition rate, [µmol P / m3 / s]
                                                                       biochem_mineralisation_po4, & !< biomineralisation rate [µmol P / m3 / s]
                                                                       transport_po4_solute, &       !< transport of soluble PO4 [µmol P / m3 / s]
                                                                       transport_po4_assoc_fast, &   !< transport of fast_assoc_po4 [µmol P / m3 / s]
                                                                       transport_po4_assoc_slow, &   !< transport of slow_assoc_po4 [µmol P / m3 / s]
                                                                       weathering_po4, &             !< phosphorus weathering rate [µmol P / m3 / s]   
                                                                       mycorrhiza_uptake_po4_sl, &   !< mycorrhizal PO4 uptake [µmol P / m3 / s]
                                                                       plant_uptake_po4_sl, &        !< plant PO4 uptake [µmol P / m3 / s]
                                                                       microbial_uptake_po4_sl, &    !< microbial PO4 uptake [µmol P / m3 / s]
                                                                       gross_mineralisation_po4, &   !< gross mineralisation of PO4 [µmol P / m3 / s]
                                                                       net_mineralisation_po4        !< net mineralisation of PO4 [µmol P / m3 / s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: po4_solute, &                 !< PO4 in solution [mol / m3]
                                                                       po4_assoc_fast, &             !< fast minerally associated PO4 pool [mol/m3]
                                                                       po4_assoc_slow, &             !< slow minerally associated PO4 pool [mol/m3]
                                                                       occlusion_po4, &              !< phosphorus occlusion rate [µmol P / m3 / s]
                                                                       slow_exchange_po4             !< phosphorus desorption rate, assoc_fast --> assoc_slow [µmol P / m3 / s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(out)   :: fast_exchange_po4             !< phosphorus adsorption rate, solute --> assoc_fast [µmol P / m3 / s]

    REAL(wp), DIMENSION(nc, nsoil_sb)           :: hlp1, hlp2, hlp3, hlp4, hlp5, &
                                                   partition_coef, & !< partition coefficient to soluble po4
                                                   partition_coef_1  !< partition coefficient to absorbed po4
    INTEGER                                     :: isoil, ichunk     !< loop over nsoil_sb and nc
    LOGICAL , DIMENSION(nc, nsoil_sb)           :: flag_arr          !< logical array helper
    CHARACTER(len=*), PARAMETER                 :: routine = TRIM(modname)//':calc_p_inorg_dynamics'

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------
    !> 0.9 initilize local variables
    partition_coef(:,:)   = 0._wp
    partition_coef_1(:,:) = 0._wp
    flag_arr(:,:) = .TRUE.
    !> 1.0 limit to possible assoc_slow, simple_0d model
    !!
    IF (sb_model_scheme=="simple_0d") THEN ! if_model
      hlp1(:,1) = k_desorpt_po4_act(:,1) * po4_assoc_slow(:,1) * dtime
      hlp2(:,1) = k_adsorpt_po4_act(:,1) * po4_assoc_fast(:,1) * dtime
      hlp3(:,1) = po4_assoc_slow(:,1) + hlp2(:,1) - hlp1(:,1)  - occlusion_po4(:,1) * dtime / 1000000._wp
      WHERE (hlp3(:,1) < eps8)
        hlp3(:,1) = po4_assoc_slow(:,1) + hlp2(:,1)
        hlp4(:,1) = hlp1(:,1) + occlusion_po4(:,1) * dtime / 1000000._wp
        hlp1(:,1) = hlp1(:,1) * hlp3(:,1) / hlp4(:,1)
        occlusion_po4(:,1) = occlusion_po4(:,1) * hlp3(:,1) / hlp4(:,1)
        slow_exchange_po4(:,1) = (hlp2(:,1) - hlp1(:,1)) / dtime * 1000000._wp
      ENDWHERE
      
      !> 1.1 Langmuir equilibrium between soluble P and labile P, 0d model
      !!
      hlp1(:,:) = 0.0_wp
      hlp2(:,:) = 0.0_wp
      !> 1.1.1 traditional Langmuir isotherm
      !!
      IF (.NOT. flag_sb_double_langmuir) THEN ! if_iso1
        hlp1(:,1) = (weathering_po4(:,1) + biochem_mineralisation_po4(:,1) + net_mineralisation_po4(:,1) + &
                    transport_po4_solute(:,1) + p_deposition_sl(:,1) - &
                    (plant_uptake_po4_sl(:,1) + mycorrhiza_uptake_po4_sl(:,1) + slow_exchange_po4(:,1)) + &
                    transport_po4_assoc_fast(:,1)) * dtime / 1000000._wp
        hlp2(:,1) = (slow_exchange_po4(:,1) - transport_po4_assoc_fast(:,1)) * dtime / 1000000._wp
        CALL calc_langmuir_kinetics(hlp1(:,1), & 
                                    po4_solute(:,1),po4_assoc_fast(:,1),  & 
                                    km_adsorpt_po4_act(:,1), qmax_po4_mineral_act(:,1), &
                                    partition_coef(:,1))

      !> 1.1.2 double-surface Langmuir isotherm
      !!
      ELSE
        flag_arr(:,1) = ((qmax_slow_po4_act(:,1) - po4_assoc_slow(:,1))<eps4 .AND. po4_assoc_slow(:,1) > eps4)
        WHERE (flag_arr(:,1))
          slow_exchange_po4(:,1) = - k_desorpt_po4_act(:,1) * po4_assoc_slow(:,1) * 1000000._wp
        ELSEWHERE
          slow_exchange_po4(:,1) = 0._wp
        ENDWHERE      
        hlp1(:,1) = (weathering_po4(:,1) + biochem_mineralisation_po4(:,1) + net_mineralisation_po4(:,1) + &
                      transport_po4_solute(:,1) + p_deposition_sl(:,1) - &
                      (plant_uptake_po4_sl(:,1) + mycorrhiza_uptake_po4_sl(:,1) + slow_exchange_po4(:,1)) + &
                      transport_po4_assoc_fast(:,1)) * dtime / 1000000._wp
        hlp2(:,1) = (- transport_po4_assoc_fast(:,1)) * dtime / 1000000._wp
        
        DO ichunk = 1, nc
          IF (flag_arr(ichunk,1)) THEN
            CALL calc_langmuir_kinetics(hlp1(:,1), & 
                                        po4_solute(:,1),po4_assoc_fast(:,1), & 
                                        km_fast_po4_act(:,1), qmax_fast_po4_act(:,1), &
                                        partition_coef(:,1))
          ELSE                            
            CALL calc_langmuir_kinetics(hlp1(:,1), & 
                                        po4_solute(:,1),po4_assoc_fast(:,1),  & 
                                        km_adsorpt_po4_act(:,1), qmax_po4_mineral_act(:,1), &
                                        partition_coef(:,1), &
                                        po4_assoc_slow(:,1), &
                                        km_fast_po4_act(:,1), qmax_fast_po4_act(:,1), &
                                        km_slow_po4_act(:,1), qmax_slow_po4_act(:,1), &
                                        partition_coef_1(:,1), slow_exchange_po4(:,1))
          ENDIF                      
        ENDDO
      ENDIF ! if_iso1
    !> 2.0 limit to possible assoc_slow, 1d model and jsm
    !!
    ELSE ! if_model
      hlp1(:,:) = k_desorpt_po4_act(:,:) * po4_assoc_slow(:,:) * dtime
      hlp2(:,:) = k_adsorpt_po4_act(:,:) * po4_assoc_fast(:,:) * dtime
      hlp3(:,:) = po4_assoc_slow(:,:) + hlp2(:,:) - hlp1(:,:) + (transport_po4_assoc_slow(:,:) - occlusion_po4(:,:)) * &
                  dtime / 1000000._wp
      WHERE (hlp3(:,:) < eps8)
        hlp3(:,:) = po4_assoc_slow(:,:) + transport_po4_assoc_slow(:,:) * dtime / 1000000._wp + hlp2(:,:)
        hlp4(:,:) = hlp1(:,:) + occlusion_po4(:,:) * dtime / 1000000._wp
        WHERE (hlp4(:,:) < eps8)
          hlp1(:,:) = 0._wp
          occlusion_po4(:,:) = 0._wp
        ELSEWHERE
          hlp1(:,:) = hlp1(:,:) * hlp3(:,:) / hlp4(:,:)
          occlusion_po4(:,:) = occlusion_po4(:,:) * hlp3(:,:) / hlp4(:,:)            
        ENDWHERE
        slow_exchange_po4(:,:) = (hlp2(:,:) - hlp1(:,:)) / dtime * 1000000._wp
      ENDWHERE

      !> 2.1 Langmuir equilibrium between soluble P and labile P, simple_1d and jsm model
      !!
      
      ! !> 2.1.1 jsm
      ! !!
      ! IF (sb_adsorp_scheme /= 'eca_none' .AND. sb_model_scheme == "jsm") THEN
        ! hlp2(:,:) = (weathering_po4(:,:) + biochem_mineralisation_po4(:,:) + net_mineralisation_po4(:,:) + &
                     ! transport_po4_solute(:,:) + p_deposition_sl(:,:) - &
                     ! (plant_uptake_po4_sl(:,:) + mycorrhiza_uptake_po4_sl(:,:) + slow_exchange_po4(:,:)) + &
                     ! transport_po4_assoc_fast(:,:)) * dtime / 1000000._wp
        ! hlp3(:,:) = (po4_solute(:,:) + po4_assoc_fast(:,:) + hlp2(:,:))
        ! hlp5(:,:) = po4_assoc_fast(:,:)

        ! !! update pools
        ! !! 
        ! po4_assoc_fast(:,:)    = calc_fast_po4(hlp3, qmax_po4_mineral_act, km_adsorpt_po4_act)
        ! po4_solute(:,:)        = hlp3(:,:) - po4_assoc_fast(:,:)
        ! fast_exchange_po4(:,:) = po4_assoc_fast(:,:) - hlp5(:,:) + (slow_exchange_po4(:,:) - transport_po4_assoc_fast(:,:)) * &
                                 ! dtime / 1000000._wp
        ! fast_exchange_po4(:,:) = fast_exchange_po4(:,:) / dtime * 1000000._wp 
        ! WHERE (po4_assoc_fast(:,:) > qmax_po4_mineral_act(:,:))
          ! hlp4(:,:)              = po4_assoc_fast(:,:) - qmax_po4_mineral_act(:,:)
          ! po4_solute(:,:)        = po4_solute(:,:) + hlp4(:,:)
          ! po4_assoc_fast(:,:)    = qmax_po4_mineral_act(:,:)
          ! fast_exchange_po4(:,:) = fast_exchange_po4(:,:) - hlp4(:,:)
        ! ENDWHERE

      !> 2.1.1 traditional Langmuir isotherm
      !!
      IF (.NOT. flag_sb_double_langmuir) THEN ! if_iso2
        hlp1(:,:) = (weathering_po4(:,:) + biochem_mineralisation_po4(:,:) + net_mineralisation_po4(:,:) + &
                     transport_po4_solute(:,:) + p_deposition_sl(:,:) - &
                     (plant_uptake_po4_sl(:,:) + mycorrhiza_uptake_po4_sl(:,:) + slow_exchange_po4(:,:)) + &
                     transport_po4_assoc_fast(:,:)) * dtime / 1000000._wp
        hlp2(:,:) = (slow_exchange_po4(:,:) - transport_po4_assoc_fast(:,:)) * dtime / 1000000._wp

        CALL calc_langmuir_kinetics(hlp1(:,:), & 
                                    po4_solute(:,:),po4_assoc_fast(:,:), & 
                                    km_adsorpt_po4_act(:,:), qmax_po4_mineral_act(:,:), &
                                    partition_coef(:,:))
      !> 2.1.2 double-surface Langmuir isotherm
      !!
      ELSE
        flag_arr(:,:) = ((qmax_slow_po4_act(:,:) - po4_assoc_slow(:,:))<eps4 .AND. po4_assoc_slow(:,:) > eps4)
        WHERE (flag_arr(:,:))
          slow_exchange_po4(:,:) = - k_desorpt_po4_act * po4_assoc_slow(:,:) * 1000000._wp
        ELSEWHERE
          slow_exchange_po4(:,:) = 0._wp
        ENDWHERE
        hlp1(:,:) = (weathering_po4(:,:) + biochem_mineralisation_po4(:,:) + net_mineralisation_po4(:,:) + &
                     transport_po4_solute(:,:) + p_deposition_sl(:,:) - &
                     (plant_uptake_po4_sl(:,:) + mycorrhiza_uptake_po4_sl(:,:) + slow_exchange_po4(:,:)) + &
                     transport_po4_assoc_fast(:,:)) * dtime / 1000000._wp
        hlp2(:,:) = ( - transport_po4_assoc_fast(:,:)) * dtime / 1000000._wp  
        
        DO ichunk = 1, nc
          DO isoil = 1, nsoil_sb
            IF (flag_arr(ichunk,isoil)) THEN
              CALL calc_langmuir_kinetics(hlp1(ichunk,isoil), & 
                                          po4_solute(ichunk,isoil),po4_assoc_fast(ichunk,isoil), & 
                                          km_fast_po4_act(ichunk,isoil), qmax_fast_po4_act(ichunk,isoil), &
                                          partition_coef(ichunk,isoil))
            ELSE                            
              CALL calc_langmuir_kinetics(hlp1(ichunk,isoil), & 
                                          po4_solute(ichunk,isoil),po4_assoc_fast(ichunk,isoil), & 
                                          km_adsorpt_po4_act(ichunk,isoil), qmax_po4_mineral_act(ichunk,isoil), &
                                          partition_coef(ichunk,isoil), &
                                          po4_assoc_slow(ichunk,isoil), &
                                          km_fast_po4_act(ichunk,isoil), qmax_fast_po4_act(ichunk,isoil), &
                                          km_slow_po4_act(ichunk,isoil), qmax_slow_po4_act(ichunk,isoil), &
                                          partition_coef_1(ichunk,isoil), slow_exchange_po4(ichunk,isoil))
            ENDIF
          ENDDO                     
        ENDDO
      ENDIF ! if_iso2
    ENDIF ! if_model
    
    !> 3.0 update po4_solute, fast_assoc_po4, and slow_assoc_po4
    !!
    po4_solute(:,:)        = po4_solute(:,:) + partition_coef(:,:) * hlp1(:,:)
    po4_assoc_fast(:,:)    = po4_assoc_fast(:,:)  + (1._wp - partition_coef(:,:) - partition_coef_1(:,:)) * hlp1(:,:) 
    fast_exchange_po4(:,:) = ((1._wp-partition_coef(:,:)-partition_coef_1(:,:)) * hlp1(:,:)+hlp2(:,:)) * 1.e6_wp / dtime
    po4_assoc_slow(:,:)    = po4_assoc_slow(:,:) + &
                            (slow_exchange_po4(:,:) - occlusion_po4(:,:) + transport_po4_assoc_slow(:,:)) * dtime / 1.e6_wp
  END SUBROUTINE calc_p_inorg_dynamics


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to calculate Langmuir kinetics (modified from Yang et al. 2014)
  !!
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_langmuir_kinetics(solute_flux, solute, assoc_fast, & 
                                              km, qmax, solute_partition_coef, &
                                              assoc_slow, km_fast, qmax_fast, km_slow, qmax_slow, &
                                              slow_partition_coef, slow_exchange)
    USE mo_sb_constants
    USE mo_jsb_math_constants,  ONLY: eps4, eps8, dtime

    IMPLICIT NONE
    ! 0.1 inout
    REAL(wp)          , INTENT(in)    :: solute_flux            !< total solute flux (mol/m3/dtime)
    REAL(wp)          , INTENT(in)    :: solute                 !< solute pool (mol/m3)
    REAL(wp)          , INTENT(in)    :: assoc_fast             !< adsorbed (fast) pool (mol/m3)
    REAL(wp)          , INTENT(in)    :: km                     !< half-saturation concentration of Langmuir isotherm (mol/m3)
    REAL(wp)          , INTENT(in)    :: qmax                   !< maximum sorption capacity of adsorbed pool (mol/m3)
    REAL(wp)          , INTENT(inout) :: solute_partition_coef  !< partition coefficient to soluble pool
    REAL(wp), OPTIONAL, INTENT(in)    :: assoc_slow             !< absorbed (slow) pool (mol/m3)
    REAL(wp), OPTIONAL, INTENT(in)    :: km_fast                !< half-saturation concentration of Langmuir isotherm (mol/m3)
    REAL(wp), OPTIONAL, INTENT(in)    :: qmax_fast              !< maximum sorption capacity of adsorbed pool (mol/m3)
    REAL(wp), OPTIONAL, INTENT(in)    :: km_slow                !< half-saturation concentration of Langmuir isotherm (mol/m3)
    REAL(wp), OPTIONAL, INTENT(in)    :: qmax_slow              !< maximum sorption capacity of adsorbed pool (mol/m3)
    REAL(wp), OPTIONAL, INTENT(inout) :: slow_partition_coef    !< partition coefficient to absorbed pool (po4 only)
    REAL(wp), OPTIONAL, INTENT(inout) :: slow_exchange          !< slow exchange rate between soluble and absorbed po4 pool [µmol P / m3 / s], double Langmuir only

    REAL(wp)                          :: hlp1, hlp2, hlp3, hlp4   
    CHARACTER(len=*), PARAMETER       :: routine = TRIM(modname)//':calc_langmuir_kinetics'

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------
    !> calculate partition coefficient of soluble pool and absorbed pool (po4 double Langmuir only)
    !!
    
    
    !> 1.0 double Langmuir isotherm (po4 only)
    !!
    IF (PRESENT(assoc_slow)) THEN
      hlp1 = qmax_fast * km_fast / ((solute + km_fast)**2.0_wp )
      hlp2 = qmax_slow * km_slow / ((solute + km_slow)**2.0_wp )
      hlp3 = 1._wp
      IF ((qmax_fast * km_fast) < eps8 .AND. (qmax_slow * km_slow) < eps8) THEN
        hlp1 = 0._wp
        hlp2 = 0._wp
      ELSEIF ((qmax_fast * km_fast) < eps8) THEN
        hlp1 = 0._wp
      ELSEIF ((qmax_slow * km_slow) < eps8) THEN
        hlp2 = 0._wp
      ELSE
        IF ( assoc_fast < (-solute_flux) .OR. assoc_slow < (-solute_flux)) THEN
          hlp1 = assoc_fast / (assoc_fast + assoc_slow + solute)
          hlp2 = assoc_slow / (assoc_fast + assoc_slow + solute)
          hlp3 = solute     / (assoc_fast + assoc_slow + solute)
        ELSEIF ( (solute < po4_solute_min) .AND. solute_flux<0._wp ) THEN
          hlp3 = 0._wp
        ENDIF
      ENDIF
      solute_partition_coef = hlp3 / (hlp1 + hlp2 + hlp3)
      slow_partition_coef   = hlp2 / (hlp1 + hlp2 + hlp3)
      slow_exchange = slow_exchange + &
                      hlp2 / (hlp1 + hlp2 + hlp3) * solute_flux / dtime * 1000000._wp

    !> 1.1 check fast pool size
      hlp4 = assoc_fast + (1._wp - solute_partition_coef - slow_partition_coef) * solute_flux      
      IF ((hlp4 > qmax) .AND. (qmax > eps8)) THEN 
         IF (solute_flux < eps8) THEN
           solute_partition_coef = 0._wp
           slow_partition_coef   = 0._wp
         ELSEIF (solute_flux > eps8) THEN
           solute_partition_coef = solute_partition_coef / (solute_partition_coef + slow_partition_coef)
           slow_partition_coef   = slow_partition_coef   / (solute_partition_coef + slow_partition_coef)
         ENDIF
      ENDIF
    
    !> 2.0 traditional Langmuir isotherm
    ELSE
      IF((km > eps8) .AND. (qmax > eps8)) THEN
        !! calculate partition coefficient based on derivatives of solute
        hlp1 = qmax * km / ((solute + km) ** 2.0_wp )
        hlp3 = 1._wp
        solute_partition_coef = hlp3 / (hlp1 + hlp3)
        !! check for extreme low concentrations of soluble and adsorbed pools
        IF (assoc_fast < (-solute_flux)) THEN
          hlp1 = assoc_fast / (assoc_fast + solute)
          !! avoid assoc_fast to be lower than po4_solute_min
          IF ((assoc_fast + hlp1 * solute_flux) < po4_solute_min) THEN
            hlp1 = MAX(0._wp, (assoc_fast - po4_solute_min)) / (-solute_flux)
          ENDIF
          hlp3                  = 1._wp - hlp1
          solute_partition_coef = hlp3 / (hlp1 + hlp3)
        ELSEIF ((solute < po4_solute_min) .AND. solute_flux < 0._wp) THEN
          hlp3                  = 0._wp
          solute_partition_coef = hlp3 / (hlp1 + hlp3)
        ELSEIF (solute_flux < (-eps8)) THEN
          solute_partition_coef = MAX(hlp3 / (hlp1 + hlp3), solute / (solute + assoc_fast))
        ENDIF
      ELSE
        hlp1                  = 0._wp
        hlp3                  = 1._wp
        solute_partition_coef = hlp3 / (hlp1 + hlp3)
      ENDIF

      !> 2.1 check fast pool size against Qmax
      !! 
      hlp4 = assoc_fast + (1._wp - solute_partition_coef) * solute_flux    
      IF ((solute + assoc_fast) > eps8 ) THEN      
         IF ((hlp4 > qmax) .AND. (qmax > eps8)) THEN 
           IF (solute_flux < 0._wp) THEN
             solute_partition_coef = 1._wp - MIN(1._wp, MAX(0._wp, (qmax-po4_solute_min) / (-solute_flux)))
           ELSEIF (solute_flux > 0._wp) THEN
             solute_partition_coef = 1._wp
           ENDIF
         ENDIF
      ENDIF
    ENDIF  ! IF (PRESENT(assoc_slow))

  END SUBROUTINE calc_langmuir_kinetics

END MODULE mo_sb_update_pools
