!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., and Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Contains the main routines for the jsm soil model
!!
!! @author
!!  Sönke Zaehle 
!!
!! @par Revision History
!! First version by Jan Engel (20170508): copyied from mo_sb_update_pools.f90
!!
MODULE mo_sb_jsm_main

#include "dsl_jsbach.inc"

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, message_text, finish
  USE mo_jsb_control,     ONLY: debug_on

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_jsm

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_jsm_main'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil biogeochemical pools 
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sb_jsm(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_jsb_math_constants,    ONLY: dtime, eps8
    USE mo_isotope_util,          ONLY: calc_fractionation
    USE mo_sb_constants
    USE mo_sb_litter_processes,   ONLY: calc_litter_partitioning
    USE mo_sb_jsm_processes,      ONLY: calc_inorganic_uptake, calc_sinking_flux_po4, calc_litter_turnover, &
                                        calc_microbial_growth, calc_sourcing_flux_po4, &
                                        calc_microbial_turnover, calc_om_sorption, &
                                        calc_depolymerisation, calc_depolymerisation_SESAM
    USE mo_sb_jsm_transport,      ONLY: calc_liquid_phase_transport_wrapper, calc_bioturbation_transport_wrapper, &
                                        calc_bioturbation_rate, calc_particle_fluxrate, calc_particle_transport_wrapper

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)

    ! Use of process memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4slm_Use_memory_a_sb
    USE mo_sb_memory_class,       ONLY: sb_pool, sb_formation, sb_loss, sb_transport
    dsl4jsb_Use_memory(SSE_)       !  USE mo_sse_memory_class,     ONLY: t_sse_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_sse_memory
    dsl4slm_Use_memory_a_veg        !  USE mo_veg_memory_class,      ONLY: t_veg_element_a, t_veg_pool_a
    USE mo_veg_memory_class,      ONLY: veg_pool, veg_litterfall
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
    INTEGER                               :: i, isoil
    INTEGER                               :: iblk, ics, ice, nc
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fine_root_carbon_sl, p_deposition_sl, & 
                                             km1_uptake_po4_act, km2_uptake_po4_act, &
                                             eca_divisor_mic, eca_divisor_plant, &
                                             arr_hlp1, arr_hlp2, arr_hlp3, arr_hlp4, arr_hlp5
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_jsm'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)                               ! TYPE(t_sb_memory),       POINTER :: sb__mem
    dsl4jsb_Def_memory(SSE_)                             ! TYPE(t_sse_memory),     POINTER :: sse__mem
    dsl4jsb_Def_memory(VEG_) 
    dsl4jsb_Def_memory(A2L_)                              ! TYPE(t_veg_memory),      POINTER :: a2l__mem
    ! Declare pointers to variables in memory
    ! SSE_
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: percolation_sl 
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    ! A2L
    dsl4jsb_Real2D_onChunk :: p_deposition
    ! VEG_
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    ! SB_ 2D
    dsl4jsb_Real2D_onChunk :: leaching_dom_carbon
    dsl4jsb_Real2D_onChunk :: leaching_dom_nitrogen
    dsl4jsb_Real2D_onChunk :: leaching_dom_phosphorus
    dsl4jsb_Real2D_onChunk :: leaching_dom_carbon13
    dsl4jsb_Real2D_onChunk :: leaching_dom_carbon14
    dsl4jsb_Real2D_onChunk :: leaching_dom_nitrogen15
    dsl4jsb_Real2D_onChunk :: leaching_nh4_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_solute
    dsl4jsb_Real2D_onChunk :: leaching_po4_solute
    dsl4jsb_Real2D_onChunk :: leaching_nh4_n15_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_n15_solute
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk :: rtm_depolymerisation 
    dsl4jsb_Real3D_onChunk :: rtm_mic_uptake
    dsl4jsb_Real3D_onChunk :: rtm_sorption
    dsl4jsb_Real3D_onChunk :: rtm_desorption
    dsl4jsb_Real3D_onChunk :: rtm_hsc
    dsl4jsb_Real3D_onChunk :: rtm_plant_uptake
    dsl4jsb_Real3D_onChunk :: rmm_depolymerisation 
    dsl4jsb_Real3D_onChunk :: rmm_mic_uptake
    dsl4jsb_Real3D_onChunk :: rmm_sorption
    dsl4jsb_Real3D_onChunk :: rmm_desorption
    dsl4jsb_Real3D_onChunk :: rmm_hsc
    dsl4jsb_Real3D_onChunk :: rmm_plant_uptake
    dsl4jsb_Real3D_onChunk :: het_respiration
    dsl4jsb_Real3D_onChunk :: het_respiration_c13
    dsl4jsb_Real3D_onChunk :: het_respiration_c14
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation 
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
    dsl4jsb_Real3D_onChunk :: microbial_cue_eff
    dsl4jsb_Real3D_onChunk :: microbial_cue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk :: weathering_po4
    dsl4jsb_Real3D_onChunk :: occlusion_po4
    dsl4jsb_Real3D_onChunk :: slow_exchange_po4
    dsl4jsb_Real3D_onChunk :: fast_exchange_po4
    dsl4jsb_Real3D_onChunk :: fast_desorp_po4
    dsl4jsb_Real3D_onChunk :: fast_adsorpt_po4
    dsl4jsb_Real3D_onChunk :: biochem_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: gross_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: nh4_solute 
    dsl4jsb_Real3D_onChunk :: nh4_assoc 
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute 
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast 
    dsl4jsb_Real3D_onChunk :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: po4_occluded
    dsl4jsb_Real3D_onChunk :: po4_primary
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: nh4_n15_assoc 
    dsl4jsb_Real3D_onChunk :: no3_n15_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_assoc 
    dsl4jsb_Real3D_onChunk :: transport_no3_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_solute 
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_fast 
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_slow 
    dsl4jsb_Real3D_onChunk :: transport_po4_occluded
    dsl4jsb_Real3D_onChunk :: transport_po4_primary
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_assoc 
    dsl4jsb_Real3D_onChunk :: transport_no3_n15_solute 
    dsl4jsb_Real3D_onChunk :: k_bioturb 
    dsl4jsb_Real3D_onChunk :: qmax_org 
    dsl4jsb_Real3D_onChunk :: qmax_po4
    dsl4jsb_Real3D_onChunk :: km_adsorpt_po4_sl
    dsl4jsb_Real3D_onChunk :: eca_divisor_adsorp
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl 
    dsl4jsb_Real3D_onChunk :: vmax_weath_mineral_sl
    dsl4jsb_Real3D_onChunk :: particle_fluxrate 
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly
    dsl4jsb_Real3D_onChunk :: enzyme_frac_residue
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_c
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_n
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_p
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_c_mavg
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_n_mavg
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_p_mavg
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
    ALLOCATE(fine_root_carbon_sl(nc, nsoil_sb))
    ALLOCATE(p_deposition_sl(nc, nsoil_sb))
    ALLOCATE(km1_uptake_po4_act(nc, nsoil_sb))
    ALLOCATE(km2_uptake_po4_act(nc, nsoil_sb))
    ALLOCATE(eca_divisor_mic(nc, nsoil_sb))
    ALLOCATE(eca_divisor_plant(nc, nsoil_sb))
    ALLOCATE(arr_hlp1(nc, nsoil_sb))
    ALLOCATE(arr_hlp2(nc, nsoil_sb))
    ALLOCATE(arr_hlp3(nc, nsoil_sb))
    ALLOCATE(arr_hlp4(nc, nsoil_sb))
    ALLOCATE(arr_hlp5(nc, nsoil_sb))
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(A2L_,      p_deposition)                ! in
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(VEG_,      root_fraction_sl)            ! in
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SSE_,     soil_depth_sl)               ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,     percolation_sl)              ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,     w_soil_sl)                   ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_carbon)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_nitrogen)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_phosphorus)     ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_carbon13)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_carbon14)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_nitrogen15)     ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_solute)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_solute)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_po4_solute)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_n15_solute)     ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_n15_solute)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_depolymerisation)        ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_mic_uptake)              ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_sorption)                ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_desorption)              ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_hsc)                     ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_plant_uptake)            ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_depolymerisation)        ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_mic_uptake)              ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_sorption)                ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_desorption)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_hsc)                     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_plant_uptake)            ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration)             ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c13)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c14)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4_n15)  ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3_n15)  ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_po4)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_cue_eff)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_cue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,       weathering_po4)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       occlusion_po4)               ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       slow_exchange_po4)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       fast_exchange_po4)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       fast_desorp_po4)             ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       fast_adsorpt_po4)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       biochem_mineralisation_po4)  ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       gross_mineralisation_po4)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_assoc)                   ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_fast)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_slow)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_occluded)                ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_primary)                 ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_assoc)               ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_solute)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_assoc)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_solute)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_solute)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_fast)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_slow)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_occluded)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_primary)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_solute)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_assoc)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_n15_solute)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       k_bioturb)                   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_org)                    ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_po4)                    ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       km_adsorpt_po4_sl)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       eca_divisor_adsorp)          ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       vmax_weath_mineral_sl)       ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       particle_fluxrate)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly)            ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_residue)         ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly_c)          ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly_n)          ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly_p)          ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly_c_mavg)     ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly_n_mavg)     ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       enzyme_frac_poly_p_mavg)     ! inout
    ! Get process memory pool structures
    dsl4slm_get_vegpool_var(pool,fine_root,carbon)
#define Get_veg_litterfall
#define Get_sb_pool
#define Get_sb_formation
#define Get_sb_loss
#define Get_sb_transport
#include "slm_pools_a_incl.f90"
#undef Get_veg_litterfall
#undef Get_sb_pool
#undef Get_sb_formation
#undef Get_sb_loss
#undef Get_sb_transport


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

!BLARPP
!percolation_sl(:,nsoil_sb) = 0.0_wp
! rmm_mic_uptake = 1.0_wp
! rmm_depolymerisation = 1.0_wp
! rmm_hsc = 1.0_wp
! rmm_sorption = 1.0_wp
! rmm_desorption = 1.0_wp

!    fine_root_carbon_sl(:,:) = root_fraction_sl(:,:) * 600._wp 
    DO isoil = 1,nsoil_sb
       WHERE(soil_depth_sl(:,isoil) > eps8)
         fine_root_carbon_sl(:,isoil) = root_fraction_sl(:,isoil) / soil_depth_sl(:,isoil) * veg_pool%fine_root%carbon(:)
       ELSEWHERE
         fine_root_carbon_sl(:,isoil) = 0.0_wp
       ENDWHERE
    ENDDO
    p_deposition_sl(:,:) = 0._wp
    p_deposition_sl(:,1) = p_deposition(:) / soil_depth_sl(:,1) !µmol P / m3 / s

    !> X.0 Inorganic N assimilation
    !!
    !> X.1 PO4

    !! Calculate the P flux rates, and potential microbial N,P uptake for eca_none
    !!
    WHERE ((sb_pool%dom_assoc%carbon(:,:) + sb_pool%residue_assoc%carbon(:,:)) > eps8)
        arr_hlp2(:,:) = (sb_pool%dom_assoc%phosphorus(:,:) + sb_pool%residue_assoc%phosphorus(:,:))/ &
                        (sb_pool%dom_assoc%carbon(:,:) + sb_pool%residue_assoc%carbon(:,:))
    ELSEWHERE
        arr_hlp2(:,:) = 0._wp
    ENDWHERE
    CALL calc_sourcing_flux_po4(dsl4jsb_Config(SB_)%flag_sb_prescribe_po4, sb_pool%microbial%carbon(:,:), &
                            sb_pool%dom%carbon(:,:), sb_pool%dom%nitrogen(:,:), sb_pool%dom%phosphorus(:,:), &
                            rtm_mic_uptake(:,:), rmm_mic_uptake(:,:), rtm_hsc(:,:), rmm_hsc(:,:), &
                            rtm_sorption(:,:), rmm_sorption(:,:), rtm_desorption(:,:), rmm_desorption(:,:), &
                            sb_pool%dom_assoc%phosphorus(:,:), sb_pool%residue_assoc%phosphorus(:,:), arr_hlp2(:,:), &
                            sb_loss%dom_assoc%phosphorus(:,:), sb_loss%residue_assoc%phosphorus(:,:), & 
                            po4_solute(:,:), nh4_solute(:,:), no3_solute(:,:), & 
                            po4_assoc_fast(:,:), po4_assoc_slow(:,:), po4_primary(:,:), &
                            bulk_dens_corr_sl(:,:), fine_root_carbon_sl(:,:), vmax_weath_mineral_sl, &
                            biochem_mineralisation_po4(:,:), gross_mineralisation_po4(:,:), weathering_po4(:,:), &
                            microbial_uptake_nh4_sl(:,:), microbial_uptake_no3_sl(:,:), microbial_uptake_po4_sl(:,:), &
                            fast_desorp_po4(:,:), occlusion_po4(:,:), slow_exchange_po4(:,:))
    !> X.2 NH4
    CALL calc_inorganic_uptake(lctlib%vmax_uptake_n, &
                               TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &
                               nh4_solute(:,:), sb_pool%microbial%carbon(:,:), fine_root_carbon_sl(:,:), &
                               rtm_mic_uptake(:,:), rmm_mic_uptake(:,:), rtm_hsc(:,:), rmm_hsc(:,:), & 
                               rtm_plant_uptake(:,:), rmm_plant_uptake(:,:), km1_up_nh4, km2_up_nh4, &
                               vmax_mic_uptake_nh4, km_mic_uptake_nh4, vmax_plant_uptake_nh4, km_plant_uptake_nh4, & 
                               k_mic_uptake_nh4, k_plant_uptake_nh4, microbial_uptake_nh4_sl(:,:), plant_uptake_nh4_sl(:,:))
    microbial_uptake_nh4_n15_sl(:,:)  = microbial_uptake_nh4_sl(:,:) * &
                                          calc_fractionation(nh4_solute(:,:), nh4_n15_solute(:,:), eta_mic_uptake_nh4)
    plant_uptake_nh4_n15_sl(:,:)      = plant_uptake_nh4_sl(:,:) * &
                                          calc_fractionation(nh4_solute(:,:) ,nh4_n15_solute(:,:) ,eta_plant_uptake_nh4)

    !> X.3 NO3
    CALL calc_inorganic_uptake(lctlib%vmax_uptake_n, &
                               TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &
                               no3_solute(:,:), sb_pool%microbial%carbon(:,:), fine_root_carbon_sl(:,:), &
                               rtm_mic_uptake(:,:), rmm_mic_uptake(:,:), rtm_hsc(:,:), rmm_hsc(:,:), & 
                               rtm_plant_uptake(:,:), rmm_plant_uptake(:,:), km1_up_no3, km2_up_no3, &
                               vmax_mic_uptake_no3, km_mic_uptake_no3, vmax_plant_uptake_no3, km_plant_uptake_no3, & 
                               k_mic_uptake_no3, k_plant_uptake_no3, microbial_uptake_no3_sl(:,:), plant_uptake_no3_sl(:,:))
    microbial_uptake_no3_n15_sl(:,:)  = microbial_uptake_no3_sl(:,:) * &
                                          calc_fractionation(no3_solute(:,:), no3_n15_solute(:,:), eta_mic_uptake_no3)
    plant_uptake_no3_n15_sl(:,:)      = plant_uptake_no3_sl(:,:) * &
                                          calc_fractionation(no3_solute(:,:), no3_n15_solute(:,:), eta_plant_uptake_no3)

    !> 2.0 Litter-related fluxes  
    !!
    !> 2.1 decay of soluable litter
    !!
    CALL calc_litter_turnover(nc, nsoil_sb, tau_soluable_litter, fc_soluable2dom, &
                              rtm_depolymerisation(:,:),rmm_depolymerisation(:,:), &
                              sb_pool%soluable_litter,sb_loss%soluable_litter,sb_formation%dom, &
                              het_respiration(:,:),het_respiration_c13(:,:),het_respiration_c14(:,:))

    !> 2.2 decay of woody litter
    !!
    CALL calc_litter_turnover(nc, nsoil_sb, tau_woody_litter, fc_woody2polymeric, &
                              rtm_depolymerisation(:,:),rmm_depolymerisation(:,:), &
                              sb_pool%woody_litter,sb_loss%woody_litter,sb_formation%polymeric_litter, &
                              het_respiration(:,:),het_respiration_c13(:,:),het_respiration_c14(:,:))

    !> 2.4 co-metabolisation of polymeric litter & residues
    !!

    CALL calc_depolymerisation_SESAM(nc, nsoil_sb, rtm_hsc(:,:),rmm_hsc(:,:), &
                                    rtm_depolymerisation(:,:),rmm_depolymerisation(:,:), &
                                    sb_pool%microbial%carbon(:,:), &
                                    enzyme_frac_poly(:,:), enzyme_frac_residue(:,:), &
                                    sb_pool%polymeric_litter,sb_loss%polymeric_litter, &
                                    sb_pool%residue,sb_loss%residue,sb_formation%dom, &
                                    enzyme_frac_poly_c, enzyme_frac_poly_n, enzyme_frac_poly_p,&
                                    enzyme_frac_poly_c_mavg, enzyme_frac_poly_n_mavg, enzyme_frac_poly_p_mavg)

    !> 2.5 Update soil litter pools with litter fall
    !!
    CALL calc_litter_partitioning(nc                      , &
                                  nsoil_sb                , &
                                  lctlib%sla              , &
                                  lctlib%growthform       , &
                                  TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
                                  soil_depth_sl(:,:)      , &
                                  root_fraction_sl(:,:)   , &
                                  veg_litterfall          , &
                                  sb_formation            )
                                             
    !> 3.0 Microbial fluxes
    !!
    !> 3.1 microbial death
    !! determine mortality rate, constrained by to avoiding a complete death of the community
    CALL calc_microbial_turnover(nc, nsoil_sb, sb_pool%microbial,sb_loss%microbial,& 
                                 sb_formation%dom,sb_formation%residue)

    !> 4.0 mineral stabilisation
    !!
    !> 4.1 mineral stabilisation of organic matter
    !!
    CALL calc_om_sorption(nc, nsoil_sb, qmax_org(:,:), rtm_sorption(:,:), rmm_sorption(:,:), rtm_desorption(:,:), &
                          rmm_desorption(:,:), sb_pool%dom,sb_pool%dom_assoc,sb_pool%residue, & 
                          sb_pool%residue_assoc,sb_formation%dom,sb_formation%dom_assoc,sb_formation%residue, &
                          sb_formation%residue_assoc,sb_loss%dom,sb_loss%dom_assoc,sb_loss%residue,sb_loss%residue_assoc)

    !> 5.0 vertical transport
    !!
    !> 5.1 liquid phase transport
    !!
    CALL calc_liquid_phase_transport_wrapper(nc, nsoil_sb, & 
                                   percolation_sl(:,:), soil_depth_sl(:,:), &
                                   sb_pool%dom, nh4_solute(:,:), no3_solute(:,:), po4_solute(:,:), & 
                                   nh4_n15_solute(:,:), no3_n15_solute(:,:), & 
                                   sb_transport%dom, transport_nh4_solute(:,:), & 
                                   transport_no3_solute(:,:),transport_po4_solute(:,:), & 
                                   transport_nh4_n15_solute(:,:),transport_no3_n15_solute(:,:), &
                                   leaching_dom_carbon(:), leaching_dom_nitrogen(:), leaching_dom_phosphorus(:), & 
                                   leaching_dom_carbon13(:), leaching_dom_carbon14(:), leaching_dom_nitrogen15(:), &
                                   leaching_nh4_solute(:), leaching_no3_solute(:), leaching_po4_solute(:), &
                                   leaching_nh4_n15_solute(:), leaching_no3_n15_solute(:))

    !> 5.2 transport due to bioturbation
    !!
    k_bioturb(:,:) = calc_bioturbation_rate(nc, nsoil_sb, soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), root_fraction_sl(:,:))
    CALL calc_bioturbation_transport_wrapper(nc, nsoil_sb, k_bioturb(:,:), soil_depth_sl(:,:), &
                                    sb_pool%polymeric_litter, sb_pool%dom_assoc, sb_pool%fungi, & 
                                    sb_pool%microbial, sb_pool%residue, sb_pool%residue_assoc, & 
                                    nh4_assoc(:,:), po4_assoc_fast(:,:), po4_assoc_slow(:,:), & 
                                    po4_occluded(:,:), nh4_n15_assoc(:,:), & 
                                    sb_transport%polymeric_litter, sb_transport%dom_assoc, & 
                                    sb_transport%fungi, sb_transport%microbial, sb_transport%residue, & 
                                    sb_transport%residue_assoc, transport_nh4_assoc(:,:), & 
                                    transport_po4_assoc_fast(:,:), transport_po4_assoc_slow(:,:), & 
                                    transport_po4_occluded(:,:), transport_nh4_n15_assoc(:,:))

    !> 5.3 particle transport
    !!
    particle_fluxrate(:,:) = calc_particle_fluxrate(nc, nsoil_sb, & 
                                    soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), & 
                                    sb_formation, sb_loss) 
    CALL calc_particle_transport_wrapper(nc, nsoil_sb, particle_fluxrate(:,:), soil_depth_sl(:,:), &
                                    sb_pool, nh4_solute(:,:), nh4_assoc(:,:), no3_solute(:,:), po4_solute(:,:), & 
                                    po4_assoc_fast(:,:), po4_assoc_slow(:,:), po4_occluded(:,:), po4_primary(:,:), & 
                                    nh4_n15_solute(:,:), nh4_n15_assoc(:,:), no3_n15_solute(:,:), & 
                                    sb_transport, transport_nh4_solute(:,:), transport_nh4_assoc(:,:), & 
                                    transport_no3_solute(:,:), transport_po4_solute(:,:), & 
                                    transport_po4_assoc_fast(:,:), transport_po4_assoc_slow(:,:), transport_po4_occluded(:,:), &
                                    transport_po4_primary(:,:), & 
                                    transport_nh4_n15_solute(:,:), transport_nh4_n15_assoc(:,:), transport_no3_n15_solute(:,:))

    !! Microbial and plant P uptake
    !!
    arr_hlp1(:,:) = (po4_solute(:,:)+ transport_po4_solute(:,:))/ dtime * 1.e6_wp
    km1_uptake_po4_act(:,:) = km1_up_po4 / (rtm_hsc(:,:) * rmm_hsc(:,:)) / 1000._wp / (w_soil_sl(:,:) / soil_depth_sl(:,:))
    km2_uptake_po4_act(:,:) = km2_up_po4 * rtm_hsc(:,:) * rmm_hsc(:,:) * 1000._wp * (w_soil_sl(:,:) / soil_depth_sl(:,:))
    CALL calc_sinking_flux_po4(lctlib%vmax_uptake_p, &
                                TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), dsl4jsb_Config(SB_)%flag_sb_ECA_upt_po4, &
                                sb_pool%microbial%carbon(:,:), fine_root_carbon_sl(:,:), po4_solute(:,:), &
                                rtm_mic_uptake(:,:), rmm_mic_uptake(:,:), rtm_hsc(:,:), rmm_hsc(:,:), &
                                rtm_plant_uptake(:,:), rmm_plant_uptake(:,:), km1_uptake_po4_act(:,:), km2_uptake_po4_act(:,:), &
                                qmax_po4(:,:), km_adsorpt_po4_sl(:,:), po4_assoc_fast(:,:), microbial_uptake_po4_sl(:,:), &
                                plant_uptake_po4_sl(:,:), fast_desorp_po4(:,:), fast_adsorpt_po4(:,:), arr_hlp1(:,:), &
                                eca_divisor_adsorp(:,:), eca_divisor_mic, eca_divisor_plant)

    CALL calc_microbial_growth(nc, nsoil_sb, TRIM(dsl4jsb_Config(SB_)%sb_bnf_scheme), &
                               sb_pool%dom,sb_pool%microbial, & 
                               rtm_mic_uptake(:,:),rmm_mic_uptake(:,:),rtm_hsc(:,:),rmm_hsc(:,:),microbial_cue_eff_tmic_mavg(:,:), &
                               microbial_uptake_nh4_sl(:,:),microbial_uptake_nh4_n15_sl(:,:), & 
                               microbial_uptake_no3_sl(:,:),microbial_uptake_no3_n15_sl(:,:), microbial_uptake_po4_sl(:,:), &
                               sb_loss%dom,sb_formation%microbial, &
                               het_respiration(:,:),het_respiration_c13(:,:),het_respiration_c14(:,:), &
                               net_mineralisation_nh4(:,:),net_mineralisation_nh4_n15(:,:), & 
                               net_mineralisation_no3(:,:),net_mineralisation_no3_n15(:,:), & 
                               net_mineralisation_po4(:,:),microbial_cue_eff(:,:),asymb_n_fixation(:,:))


    !> 6.0 convert inorganic transport fluxes to µmol m3 s-1
    !!
    transport_nh4_solute(:,:)     = transport_nh4_solute(:,:) * 1.e6_wp / dtime
    transport_nh4_assoc(:,:)      = transport_nh4_assoc(:,:) * 1.e6_wp / dtime
    transport_no3_solute(:,:)     = transport_no3_solute(:,:) * 1.e6_wp / dtime
    transport_po4_solute(:,:)     = transport_po4_solute(:,:) * 1.e6_wp / dtime
    transport_po4_assoc_fast(:,:) = transport_po4_assoc_fast(:,:) * 1.e6_wp / dtime
    transport_po4_assoc_slow(:,:) = transport_po4_assoc_slow(:,:) * 1.e6_wp / dtime
    transport_po4_occluded(:,:)   = transport_po4_occluded(:,:) * 1.e6_wp / dtime
    transport_po4_primary(:,:)    = transport_po4_primary(:,:) * 1.e6_wp / dtime
    transport_nh4_n15_solute(:,:) = transport_nh4_n15_solute(:,:) * 1.e6_wp / dtime
    transport_nh4_n15_assoc(:,:)  = transport_nh4_n15_assoc(:,:) * 1.e6_wp / dtime
    transport_no3_n15_solute(:,:) = transport_no3_n15_solute(:,:)  * 1.e6_wp / dtime


  ! Deallocate local allocatable variables
  DEALLOCATE(fine_root_carbon_sl)
  DEALLOCATE(p_deposition_sl)
  DEALLOCATE(km1_uptake_po4_act)
  DEALLOCATE(km2_uptake_po4_act)
  DEALLOCATE(eca_divisor_mic)
  DEALLOCATE(eca_divisor_plant)
  DEALLOCATE(arr_hlp1)
  DEALLOCATE(arr_hlp2)
  DEALLOCATE(arr_hlp3)
  DEALLOCATE(arr_hlp4)
  DEALLOCATE(arr_hlp5)

! Pass back values to process memory
dsl4slm_passback_vegpool_var(pool,fine_root,carbon)
#define Passback_veg_litterfall
#define Passback_sb_pool
#define Passback_sb_formation
#define Passback_sb_loss
#define Passback_sb_transport
#include "slm_pools_a_incl.f90"
#undef Passback_veg_litterfall
#undef Passback_sb_pool
#undef Passback_sb_formation
#undef Passback_sb_loss
#undef Passback_sb_transport

  END SUBROUTINE update_sb_jsm

END MODULE mo_sb_jsm_main
