!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> Calculate post-processed Spitfire fire impacts and model feedback
!>
!> Gather data from memory, rebuild variables, ...
!> (PFT is tree?) calc prob of tree mortality
!> (PFT is tree?) calc live fuel consumption
!> Calc dead fuel consumption
!> Allocate fuel consumption for remaining nutrients and isotopes
!> Calc fire emissions
!> Calc fire residue
!> (online?) Apply impacts on the model
!> Calc diagnostics varibles
!> Adjust flux units
!> (enabled?) Run diagnostics
!>
MODULE mo_q_dfire_sf_impacts
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_lnd_bgcm_store,        ONLY: t_lnd_bgcm_store
  USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract

  USE mo_veg_constants,         ONLY: itree, igrass
  USE mo_q_dfire_sf_adaptor,    ONLY: sf_convert_units
  USE mo_jsb_grid,              ONLY: Get_grid

  USE mo_q_dfire_sf_constants,  ONLY: df_sf_spq_layer_top_sl
  USE mo_q_dfire_sf_constants,  ONLY: df_sf_fc_spread_nr, df_sf_fc_full_nr
  USE mo_q_dfire_sf_constants,  ONLY: df_sf_fc_1hr_idx, df_sf_fc_10hr_idx, df_sf_fc_100hr_idx, df_sf_fc_1000hr_idx
  USE mo_q_dfire_sf_constants,  ONLY: sf_fuel_dead_fuel_nr, sf_fuel_live_fuel_nr
  USE mo_q_dfire_sf_constants,  ONLY: df_sf_elements_nr
  USE mo_q_dfire_sf_constants,  ONLY: df_sf_species_nr

  USE mo_q_dfire_sf_constants,  ONLY: ixs_sf_dead_litter_fuel_pools

  USE mo_q_dfire_sf_constants,  ONLY: df_sf_species_nr, df_sf_species_co2_idx, df_sf_species_co_idx, &
                                        df_sf_species_ch4_idx, &
                                        !df_sf_species_voc_idx, df_sf_species_tpm_idx, & df_sf_species_nox_idx, &
                                        df_sf_species_nh3_idx, df_sf_species_n2o_idx, &
                                        df_sf_species_n2_idx

  USE mo_q_dfire_sf_damage_class, ONLY: t_sf_prob_mort, calc_prob_mortality
  USE mo_q_dfire_sf_util,         ONLY: real2logical, dfire_update_pools

  USE mo_q_dfire_sf_fuel_consump_class, ONLY: sf_calc_live_fuel_consumption, &
                                              sf_calc_dead_fuel_consumption, &
                                            sf_calc_dead_fuel_consumption_populate_elements, &
                                            sf_calc_live_fuel_consumption_populate_elements

  USE mo_q_dfire_sf_feedback_class,   ONLY: calc_sf_feedback_dead, calc_sf_feedback_live

  USE mo_q_dfire_sf_emissions_class,  ONLY: sf_calc_emissions

  USE mo_lnd_bgcm_store_class,  ONLY: SB_BGCM_POOL_ID, VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID, &
                                      SB_BGCM_FIRE_EMISSIONS_ID, SB_BGCM_FIRE_LITTER_ASH_ID, &
                                      VEG_BGCM_FIRE_EMISSIONS_ID, VEG_BGCM_FIRE_LITTERFALL_ID

  USE mo_lnd_bgcm_idx

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_q_dfire_sf_impacts, sf_calc_residue

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_impacts'

CONTAINS

  ! ======================================================================================================= !
  !>
  !>
  !>
  SUBROUTINE update_q_dfire_sf_impacts(tile, options)
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_grid_class,          ONLY: t_jsb_grid
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: VEG_, HYDRO_, Q_DFIRE_FUEL_, Q_DFIRE_SF_
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_math_constants,      ONLY: one_day

    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(Q_DFIRE_FUEL_)
    dsl4jsb_Use_memory(Q_DFIRE_SF_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(Q_DFIRE_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_grid),       POINTER       :: grid
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lnd_bgcm_store), POINTER       :: bgcm_store !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc, ic, ipool, iel, ip
    INTEGER                               :: no_children
    REAL(wp)                              :: dtime      !< timestep

    ! add here local variables
    REAL(wp), POINTER                     :: area(:), tile_fract(:)

    REAL(wp)                              :: dbh(options%nc) !< Tree diameter breast at height [cm]

    REAL(wp)                              :: w_oi(options%nc, df_sf_fc_full_nr) !< Dead fuel by fuel classes [g DM m-2]
    LOGICAL                               :: is_fire_intense_enough(options%nc) !< TRUE=1, FALSE=0
    REAL(wp)                              :: fc_dead_cls(options%nc, sf_fuel_dead_fuel_nr, df_sf_fc_full_nr) !< Dead fuel consumption for current PFT [g DM m-2]
    REAL(wp)                              :: fc_live(options%nc,sf_fuel_live_fuel_nr)   !< Live fuel consumption (#nc, #live pools) for current PFT [g DM m-2]
    REAL(wp)                              :: emissions_live_spc(options%nc,sf_fuel_live_fuel_nr,df_sf_species_nr) !< emissions (#nc, #live pools, #df_sf_species_nr) [mol X Kg-1 DM-1 m-2]
    REAL(wp)                              :: emissions_dead_spc(options%nc,sf_fuel_dead_fuel_nr,df_sf_species_nr) !< emissions (#nc, #dead pools, #df_sf_species_nr) [mol X Kg-1 DM-1 m-2]

    REAL(wp)                              :: emissions_live_el(options%nc,sf_fuel_live_fuel_nr,df_sf_elements_nr) !< emissions by elements (#nc, #live pools, #elements) [mol X m-2]
    REAL(wp)                              :: emissions_dead_el(options%nc,sf_fuel_dead_fuel_nr,df_sf_elements_nr) !< emissions by elements (#nc, #dead pools, #elements) [mol X m-2]

    REAL(wp)                              :: area_burned_pft(options%nc) !< Area burned for current PFT [ha]
    REAL(wp)                              :: area_burned_fract_pft(options%nc) !< Fraction of the area burned for current PFT [-]

    TYPE(t_sf_prob_mort)                  :: prob_mort                 !< Tree mortality
    REAL(wp), ALLOCATABLE                 :: live_tree_fuel_cls(:,:,:) !< live fuel for 1, 10, 100 and 1000 hours fuel classes [g DM m-2]

    REAL(wp), ALLOCATABLE                 :: fc_live_mol_m2(:,:)  !< Total live fuel consumption (nc, live fuels) [molX/m2]
    REAL(wp), ALLOCATABLE                 :: fc_dead_mol_m2(:,:)  !< Total dead fuel consumption (nc, t_sf_fuel_dead_fuel_nr) [molX/m2]

    REAL(wp), ALLOCATABLE                 :: fc_dead_el(:,:,:)   !< dead fuel consumption (#nc, #dead pools, #Elements) [molX m-2]
    REAL(wp), ALLOCATABLE                 :: fc_live_el(:,:,:)   !< live fuel consumption (#nc, #live pools, #Elemens) [molX m-2]

    REAL(wp), ALLOCATABLE                 :: dead_fuel_fcl(:,:,:) !< Dead fuel [g MD m-2]

    REAL(wp)                              :: residue_dead_el(options%nc,sf_fuel_dead_fuel_nr,df_sf_elements_nr) !< emissions (#nc, #dead pools, #elements) [mol m-2]
    REAL(wp)                              :: residue_live_el(options%nc,sf_fuel_live_fuel_nr,df_sf_elements_nr) !< emissions (#nc, #live pools, #elements) [mol m-2]

    REAL(wp)                              :: mass_imbalance_emitted_dead_el(options%nc,sf_fuel_dead_fuel_nr,df_sf_elements_nr) !< emissions (#nc, #dead pools, #elems) [mol m-2]
    REAL(wp)                              :: mass_imbalance_emitted_live_el(options%nc,sf_fuel_live_fuel_nr,df_sf_elements_nr) !< emissions (#nc, #live pools, #elems) [mol m-2]

    REAL(wp), ALLOCATABLE                 :: sb_pool_in_m2(:,:,:)               !< BGCM sb_pool (#pools,#elements,#nc) [molX/m2]
    REAL(wp), ALLOCATABLE                 :: emissions_dead_in_m3(:,:,:)        !< (#nc,#pools,#elements) sb_pool fuels in molX/m3, veg_pool fuels in molX/m2
    REAL(wp), ALLOCATABLE                 :: residue_dead_in_m3(:,:,:)          !< (#nc,#pools,#elements) sb_pool fuels in molX/m3, veg_pool fuels in molX/m2

    REAL(wp), ALLOCATABLE                 :: lctlib_emission_factors(:)         !< Array packed with emission factors from lctlib (#species)

    INTEGER :: ix_elem, ielem, ix_comp, nr_of_veg_bgcm_comp

    CHARACTER(len=*), PARAMETER           :: routine = TRIM(modname)//':update_fdire_sf_fuel'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real3D_onChunk                :: live_tree_fuel_1hr_lfp    ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_tree_fuel_10hr_lfp   ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_tree_fuel_100hr_lfp  ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_tree_fuel_1000hr_lfp ! in [g DM m-2]

    dsl4jsb_Real3D_onChunk                :: dead_fuel_1hr_dfp         ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_10hr_dfp        ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_100hr_dfp       ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_1000hr_dfp      ! in [g DM m-2]

    dsl4jsb_Real2D_onChunk                :: area_burned               ! in [ha ts-1]
    dsl4jsb_Real2D_onChunk                :: area_burned_fract         ! in [-]
    dsl4jsb_Real3D_onChunk                :: dead_frac_consumption     ! in [-]
    dsl4jsb_Real2D_onChunk                :: i_surface                 ! in [kW m-1]
    dsl4jsb_Real2D_onChunk                :: mask_is_fire_intense_enough ! in [-]
    dsl4jsb_Real2D_onChunk                :: pm_prop_ck                ! out [-]
    dsl4jsb_Real2D_onChunk                :: pm_fire                   ! out, Tree probability of mortality by fire [0 to 1 ts-1]
    dsl4jsb_Real2D_onChunk                :: pm_fire_output            ! out, Tree probability of mortality by fire [0 to 1 s-1]
    dsl4jsb_Real2D_onChunk                :: pm_tau                    ! out [-]
    dsl4jsb_Real2D_onChunk                :: tau_l                     ! out [-]
    dsl4jsb_Real2D_onChunk                :: tau_c                     ! out [-]
    dsl4jsb_Real2D_onChunk                :: bark_thickness            ! out [?]
    dsl4jsb_Real2D_onChunk                :: DBH_diag                  ! out [?]
    dsl4jsb_Real2D_onChunk                :: par1_diag                 ! out [?]
    dsl4jsb_Real2D_onChunk                :: par2_diag                 ! out [?]

    dsl4jsb_Real2D_onChunk                :: ck                        ! out [?]
    dsl4jsb_Real2D_onChunk                :: r_ck                      ! out [?]
    dsl4jsb_Real2D_onChunk                :: ck_p                      ! out [?]
    dsl4jsb_Real2D_onChunk                :: scorch_height             ! out [?]
    dsl4jsb_Real2D_onChunk                :: F_scorch_height           ! out [?]
    dsl4jsb_Real2D_onChunk                :: crown_length              ! out [?]
    dsl4jsb_Real2D_onChunk                :: tree_height               ! out [m]

    dsl4jsb_Real3D_onChunk                :: dead_fuel_consumption_1hr ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_consumption_10hr ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_consumption_100hr ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_consumption_1000hr ! in [g DM m-2]

    dsl4jsb_Real3D_onChunk                :: fuel_consumption_els     ! out [mol m-2]
    dsl4jsb_Real3D_onChunk                :: emissions_spc            ! out
    dsl4jsb_Real3D_onChunk                :: emissions_els            ! out
    dsl4jsb_Real3D_onChunk                :: residue_els              ! out
    dsl4jsb_Real3D_onChunk                :: overflow_els             ! out

    dsl4jsb_Real2D_onChunk                :: diameter                   ! in [?]
    dsl4jsb_Real2D_onChunk                :: height                     ! in [?]
    dsl4jsb_Real3D_onChunk                :: soil_depth_sl              ! in [m]

    dsl4jsb_Real3D_onChunk                :: veg_fire_litterfall_total_els ! out [umol X m-2 s-1]
    dsl4jsb_Real3D_onChunk                :: veg_fire_emissions_total_els ! out [umol X m-2 s-1]
    dsl4jsb_Real3D_onChunk                :: sb_fire_emissions_total_els ! out [umol X m-2 s-1]
    dsl4jsb_Real3D_onChunk                :: sb_fire_litter_ash_total_els ! out [umol X m-2 s-1]
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_fire_emissions_mt
    dsl4jsb_Def_mt2L3D :: sb_fire_litter_ash_mt
    dsl4jsb_Def_mt2L2D :: veg_fire_emissions_mt
    dsl4jsb_Def_mt2L2D :: veg_fire_litterfall_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_DFIRE_SF_)
    dsl4jsb_Def_memory(Q_DFIRE_FUEL_)
    dsl4jsb_Def_memory(Q_DFIRE_SF_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Local vars here
    ! ----------------------------------------------------------------------------------------------------- !

    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_DFIRE_FUEL_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (lctlib%PastureFlag) RETURN !< do not run this routine at tiles like "pasture pft"
    IF (lctlib%CropFlag) RETURN !< do not run this routine at tiles like "crops pft"
    ! ----------------------------------------------------------------------------------------------------- !
    grid => Get_grid(model%grid_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_DFIRE_SF_)
    dsl4jsb_Get_memory(Q_DFIRE_FUEL_)
    dsl4jsb_Get_memory(Q_DFIRE_SF_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)               ! in [?]

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1hr_lfp)    ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_10hr_lfp)   ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_100hr_lfp)  ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1000hr_lfp) ! in [g DM m-2]

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, dead_fuel_1hr_dfp)     ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, dead_fuel_10hr_dfp)    ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, dead_fuel_100hr_dfp)   ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, dead_fuel_1000hr_dfp)  ! in [g DM m-2]

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, area_burned)               ! in [ha ts-1]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, area_burned_fract)         ! in [-]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_SF_, dead_frac_consumption)     ! in [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, i_surface)                 ! in [kW m-1]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_SF_, mask_is_fire_intense_enough)! in [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, pm_prop_ck)           ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, pm_fire)              ! out [- ts-1]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, pm_fire_output)       ! out [- s-1]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, pm_tau)               ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, tau_l)                ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, tau_c)                ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, bark_thickness)       ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, DBH_diag)             ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, par1_diag)            ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, par2_diag)            ! out [-]

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, ck)                   ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, r_ck)                 ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, ck_p)                 ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, scorch_height)        ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, F_scorch_height)      ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, crown_length)         ! out [-]
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_FUEL_, tree_height)          ! out [m]

    dsl4jsb_Get_var2D_onChunk(VEG_, diameter)                    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, height)                      ! in

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_SF_, dead_fuel_consumption_1hr) ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_SF_, dead_fuel_consumption_10hr) ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_SF_, dead_fuel_consumption_100hr) ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_SF_, dead_fuel_consumption_1000hr) ! in [g DM m-2]

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, fuel_consumption_els) ! out
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, emissions_spc)        ! out
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, emissions_els)        ! out
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, residue_els)          ! out
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, overflow_els)         ! out

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, veg_fire_litterfall_total_els) ! out [umol X m-2 s-1]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, veg_fire_emissions_total_els) ! out [umol X m-2 s-1]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, sb_fire_emissions_total_els) ! out [umol X m-2 s-1]
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, sb_fire_litter_ash_total_els) ! out [umol X m-2 s-1]
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_FIRE_EMISSIONS_ID, sb_fire_emissions_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_FIRE_LITTER_ASH_ID, sb_fire_litter_ash_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_FIRE_EMISSIONS_ID, veg_fire_emissions_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_FIRE_LITTERFALL_ID, veg_fire_litterfall_mt)
    ! ----------------------------------------------------------------------------------------------------- !

    ! init here
#ifdef __QUINCY_STANDALONE__
    ALLOCATE(area(1), tile_fract(1))
    area = 1._wp
    tile_fract = 1._wp
#else
    area => grid %area (ics:ice, iblk)
    tile_fract => tile%fract(ics:ice, iblk)
#endif
    is_fire_intense_enough = real2logical(mask_is_fire_intense_enough)
    fc_live = 0._wp
    pm_fire = 0._wp
    pm_fire_output = 0._wp
    pm_prop_ck = 0._wp
    pm_tau = 0._wp
    tau_l = 0._wp
    tau_c = 0._wp
    bark_thickness = 0._wp
    DBH_diag = 0._wp
    par1_diag = 0._wp
    par2_diag = 0._wp

    ck = 0._wp
    r_ck = 0._wp
    ck_p = 0._wp
    scorch_height = 0._wp
    F_scorch_height = 0._wp
    crown_length = 0._wp
    tree_height = 0._wp

    fuel_consumption_els = 0._wp
    emissions_spc = 0._wp
    emissions_els = 0._wp
    residue_els = 0._wp
    overflow_els = 0._wp
    veg_fire_litterfall_total_els = 0._wp
    veg_fire_emissions_total_els = 0._wp
    sb_fire_emissions_total_els = 0._wp
    sb_fire_litter_ash_total_els = 0._wp

    ALLOCATE( dead_fuel_fcl(options%nc,sf_fuel_dead_fuel_nr,df_sf_fc_full_nr) )
    dead_fuel_fcl(:,:,df_sf_fc_1hr_idx) = dead_fuel_1hr_dfp
    dead_fuel_fcl(:,:,df_sf_fc_10hr_idx) = dead_fuel_10hr_dfp
    dead_fuel_fcl(:,:,df_sf_fc_100hr_idx) = dead_fuel_100hr_dfp
    dead_fuel_fcl(:,:,df_sf_fc_1000hr_idx) = dead_fuel_1000hr_dfp
    area_burned_pft = area_burned * tile_fract
    area_burned_fract_pft = area_burned_fract * tile_fract

    IF (.NOT. ALLOCATED(lctlib_emission_factors)) THEN ! only 1st timestep
      ALLOCATE(lctlib_emission_factors(df_sf_species_nr))
      !> Pack emission factors in an array
      lctlib_emission_factors(df_sf_species_co2_idx) = lctlib%q_emissions_co2_mol
      lctlib_emission_factors(df_sf_species_co_idx) = lctlib%q_emissions_co_mol
      lctlib_emission_factors(df_sf_species_ch4_idx) = lctlib%q_emissions_ch4_mol
      lctlib_emission_factors(df_sf_species_nh3_idx) = lctlib%q_emissions_nh3_mol
      lctlib_emission_factors(df_sf_species_n2o_idx) = lctlib%q_emissions_n2o_mol
      lctlib_emission_factors(df_sf_species_n2_idx) = lctlib%q_emissions_n2_mol
      !< Decide what to do later
      !lctlib_emission_factors(df_sf_species_voc_idx) = lctlib%q_emissions_voc
      !lctlib_emission_factors(df_sf_species_tpm_idx) = lctlib%q_emissions_tpm
      !lctlib_emission_factors(df_sf_species_nrox_idx) = lctlib%q_emissions_nox
    ENDIF

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !

    ! calc tree mortality
    IF (lctlib%growthform == itree) THEN
      ALLOCATE(live_tree_fuel_cls(nc,sf_fuel_live_fuel_nr, df_sf_fc_full_nr))
      live_tree_fuel_cls(:,:,df_sf_fc_1hr_idx) = live_tree_fuel_1hr_lfp
      live_tree_fuel_cls(:,:,df_sf_fc_10hr_idx) = live_tree_fuel_10hr_lfp
      live_tree_fuel_cls(:,:,df_sf_fc_100hr_idx) = live_tree_fuel_100hr_lfp
      live_tree_fuel_cls(:,:,df_sf_fc_1000hr_idx) = live_tree_fuel_1000hr_lfp

      ! ref vars and init to 0
      prob_mort = t_sf_prob_mort(nc, pm_fire, pm_tau, tau_l, tau_c, &
                               bark_thickness, pm_prop_ck, scorch_height, &
                               crown_length )

      !< Calculate amount of dead fuel available for combustion following Thonicke et al. (2010)
      w_oi = SUM(dead_fuel_fcl, DIM=2) !< 2= dead fuel pools index, fuel class only

      dbh(:) = sf_convert_units("m", "cm", diameter) !< dbh to spitfire units

      is_fire_intense_enough = is_fire_intense_enough .AND. &
                                (tile_fract > 0._wp) .AND. & ! tile is empty, skip
                                (height > 0._wp) ! height is not yet initialized, skip

      CALL calc_prob_mortality(nc, df_sf_fc_spread_nr, &
                                    lctlib%q_crown_length, lctlib%q_crown_damage_p, lctlib%q_crown_damage_r_ck, &
                                    lctlib%q_scorch_height_f, &
                                    lctlib%q_bark_thickness_par1, lctlib%q_bark_thickness_par2, &
                                    height * dsl4jsb_Config(Q_DFIRE_SF_)%tree_height_factor, & ! Tree characteristics
                                    dbh * dsl4jsb_Config(Q_DFIRE_SF_)%tree_diameter_factor, &
                                    i_surface, & ! Fire intensity
                                    w_oi(:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx), &
                                    dead_frac_consumption(:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx), &
                                    prob_mort, & ! output
                                    mask_opt=is_fire_intense_enough)

      fc_live = sf_calc_live_fuel_consumption(nc, sf_fuel_live_fuel_nr, df_sf_fc_full_nr, &
                                    dsl4jsb_Config(Q_DFIRE_SF_)%fract_crown_scorch, &
                                    prob_mort % prop_ck, &
                                    area_burned_fract_pft, &
                                    live_tree_fuel_cls, &
                                    mask=is_fire_intense_enough)

      !< output calculated values
      pm_prop_ck = prob_mort % prop_ck
      bark_thickness = prob_mort % bark_thickness
      par1_diag = lctlib%q_bark_thickness_par1
      par2_diag = lctlib%q_bark_thickness_par2
      dbh_diag = dbh * dsl4jsb_Config(Q_DFIRE_SF_)%tree_diameter_factor

      ck = prob_mort % prop_ck
      scorch_height = prob_mort % scorch_height
      crown_length = prob_mort % crown_length
      !< Diagnostics
      r_ck = lctlib%q_crown_damage_r_ck
      ck_p = lctlib%q_crown_damage_p
      F_scorch_height = lctlib%q_scorch_height_f
      tree_height = height * dsl4jsb_Config(Q_DFIRE_SF_)%tree_height_factor

      If (dsl4jsb_Config(Q_DFIRE_SF_)%flag_has_prob_mort_feedback) THEN
        ! todo: add pm_fire_diag -> always output (even when false)
        pm_fire = pm_fire * area_burned_fract_pft
        pm_fire_output = pm_fire / dtime !< from ts-1 to seconds-1
      End If
    END IF

    !< Recalculate dead fuel consumption at a PFT level
    fc_dead_cls = sf_calc_dead_fuel_consumption(area_burned_fract_pft, dead_frac_consumption, &
                                                dead_fuel_fcl, mask=is_fire_intense_enough)

    !< Units conversion of Dead fuel consumption into molX/m2
    fc_dead_mol_m2 = sf_convert_units("gDM/m2", "molC/m2", SUM(fc_dead_cls, dim=3)) !< 3=fuel classes

    !< Units conversion of Live fuel consumption into Quincy units
    fc_live_mol_m2 = sf_convert_units("gDM/m2", "molC/m2", fc_live) !< (#nc,#live_pools)

    !< To use in fuel consumption, emissions, residue, ...
    sb_pool_in_m2 =  sf_convert_units("-/m3", "-/m2", sb_pool_mt(:,:,:,df_sf_spq_layer_top_sl), nc_dim_opt=3, &
                                      aux=soil_depth_sl(:,df_sf_spq_layer_top_sl)) !< 3=nc dimension

    !< Derive dead fuel consumption elements and isotopes from dead fuel consumption, sb_pool and veg_pool
    fc_dead_el = sf_calc_dead_fuel_consumption_populate_elements(nc, sf_fuel_dead_fuel_nr, df_sf_elements_nr, & ! [molX/m2]
                                                              veg_pool_mt, sb_pool_in_m2, & ! [molX/m2], [molX/m2]
                                                              fc_dead_mol_m2, & ! [molX/m2]
                                                              mask=is_fire_intense_enough)

    !< Derive live fuel consumption elements and isotopes from live fuel consumption and veg_pool
    fc_live_el = sf_calc_live_fuel_consumption_populate_elements(nc, sf_fuel_live_fuel_nr, & ! [molX/m2]
                                                              df_sf_elements_nr, &
                                                              veg_pool_mt, fc_live_mol_m2, & ! [molX/m2], [molX/m2]
                                                              mask=is_fire_intense_enough)

    !< Calculate fire damage to plants (total biomass burned and tree mortality) based on total dead and
    !< live fuel, as in Peiffer et al. (2013)
    !< Also calculate emissions by elements/nutrients so it can be removed from pools
    CALL sf_calc_emissions(lctlib_emission_factors, & ! in
                          SUM(fc_dead_cls, dim=3), & ! 3 = fuel classes index
                          fc_dead_el, &
                          emissions_dead_spc, emissions_dead_el, & ! out
                          mask=is_fire_intense_enough)
    CALL sf_calc_emissions(lctlib_emission_factors, & !in
                          fc_live, &
                          fc_live_el, &
                          emissions_live_spc, emissions_live_el, & ! out
                          mask=is_fire_intense_enough)

    !< Is there any C, N and/or P left after Emissions? Calculate residue
    CALL sf_calc_residue(fc_dead_el, & ! in [molX/m2]
                          emissions_dead_el, & ! [molX/m2]
                          residue_dead_el, & ! out [molX/m2]
                          mass_imbalance_emitted_dead_el) ! [molX/m2]

    CALL sf_calc_residue(fc_live_el, & ! in [molX/m2]
                          emissions_live_el, & ! [molX/m2]
                          residue_live_el, & ! out [molX/m2]
                          mass_imbalance_emitted_live_el) ! [molX/m2]

    ! Feedbacks
    IF (dsl4jsb_Config(Q_DFIRE_SF_)%flag_has_dead_feedback) THEN
      ! Convert sb_pool related dead fuel pools into m3 to operate with sb_pool related fluxes
      ALLOCATE(emissions_dead_in_m3, source=emissions_dead_el)
      emissions_dead_in_m3 = 0._wp
      emissions_dead_in_m3(:,ixs_sf_dead_litter_fuel_pools,:) = sf_convert_units("-/m2", "-/m3", &
                            emissions_dead_el(:,ixs_sf_dead_litter_fuel_pools,:), &
                            soil_depth_sl(:,df_sf_spq_layer_top_sl))

      ALLOCATE(residue_dead_in_m3, source=residue_dead_el)
      residue_dead_in_m3 = 0._wp
      residue_dead_in_m3(:,ixs_sf_dead_litter_fuel_pools,:) = sf_convert_units("-/m2", "-/m3", &
                            residue_dead_el(:,ixs_sf_dead_litter_fuel_pools,:), &
                            soil_depth_sl(:,df_sf_spq_layer_top_sl))

      !< Feedback: dead fuel consumption is removed through veg_fire_emissions, sb_fire_emissions_mt for emissions
      !< Residue is transfered through sb_fire_litter_ash and veg_fire_litterfall
      CALL calc_sf_feedback_dead(emissions_dead_in_m3, &
                                residue_dead_in_m3, &
                                sb_fire_emissions_mt(:,:,:,df_sf_spq_layer_top_sl), & ! todo: create flux pool with top layer only
                                sb_fire_litter_ash_mt(:,:,:,df_sf_spq_layer_top_sl), &
                                veg_fire_emissions_mt, veg_fire_litterfall_mt, &
                                mask=is_fire_intense_enough)

    END IF

    IF (dsl4jsb_Config(Q_DFIRE_SF_)%flag_has_live_feedback) THEN
      !< Feedback: live fuel consumption is removed using veg_fire_emissions (combusted)
      !< Residue is transfered through veg_litterfall
      CALL calc_sf_feedback_live(emissions_live_el, &
                                residue_live_el, &
                                veg_fire_emissions_mt, &
                                veg_fire_litterfall_mt, &
                                mask=is_fire_intense_enough)

    END IF

    !< Enforce to update pools due to Spitfire impacts. Temporary solution until a permanent solution
    !< is found to properly deal.
    !< Spitfire is now updated here

    !< Update pools
    CALL dfire_update_pools(sb_fire_emissions_mt(:,:,:,df_sf_spq_layer_top_sl), &
                            sb_fire_litter_ash_mt(:,:,:,df_sf_spq_layer_top_sl), &
                            sb_pool_mt(:,:,:,df_sf_spq_layer_top_sl))

    nr_of_veg_bgcm_comp = SIZE(veg_pool_mt,1)

    !> 3.0 add fire disturbance into litterfall
    DO ic = 1, nc
      DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
        IF (model%config%is_element_used(ielem)) THEN
          ix_elem = model%config%elements_index_map(ielem)    ! get element index in bgcm
          DO ix_comp = 1, nr_of_veg_bgcm_comp
            veg_litterfall_mt(ix_comp, ix_elem, ic) = veg_litterfall_mt(ix_comp, ix_elem, ic) + &
                                                      veg_fire_litterfall_mt(ix_comp, ix_elem, ic)
          END DO ! ix_comp = 1, n_veg_bgcm_comp
        END IF ! model%config%is_element_used(ielem)) THEN
     END DO ! ielem = FIRST_ELEM_ID, LAST_ELEM_ID
    END DO ! ic = 1, nc

    !> update veg fire emissions
    DO ic = 1, nc
      DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
        IF (model%config%is_element_used(ielem)) THEN
          ix_elem = model%config%elements_index_map(ielem)    ! get element index in bgcm
          DO ix_comp = 1, nr_of_veg_bgcm_comp
                veg_pool_mt(ix_comp,ix_elem,ic) = veg_pool_mt(ix_comp,ix_elem,ic) - &
                    veg_fire_emissions_mt(ix_comp,ix_elem,ic) ! dist_fire
          END DO ! ix_comp
        END IF ! is_element_used(ielem)
      END DO ! ielem
    END DO ! ic


    ! for diagnostics:
    fuel_consumption_els = (SUM(fc_dead_el, DIM=2) + SUM(fc_live_el, DIM=2)) / dtime * 1000000._wp !< 2 = fuel pool index
    emissions_spc = (SUM(emissions_dead_spc, DIM=2) + SUM(emissions_live_spc, DIM=2)) / dtime !< 2 = fuel pool index
    emissions_els = (SUM(emissions_dead_el, DIM=2) + SUM(emissions_live_el, DIM=2)) / dtime * 1000000._wp !< 2 = fuel pool index
    residue_els = (SUM(residue_dead_el, DIM=2) + SUM(residue_live_el, DIM=2)) / dtime * 1000000._wp !< 2 = fuel pool index
    overflow_els = (SUM(mass_imbalance_emitted_dead_el, DIM=2) + & !< 2 = fuel pool index
                        SUM(mass_imbalance_emitted_live_el, DIM=2)) / dtime * 1000000._wp
    ! for diagnostics (fluxes):
    DO ipool=1, SIZE(veg_fire_litterfall_mt,DIM=1)
      DO iel=1, df_sf_elements_nr
        DO ic=1, nc ! (#pool, #elements, #nc)
          veg_fire_litterfall_total_els(ic,iel) = veg_fire_litterfall_total_els(ic,iel) + veg_fire_litterfall_mt(ipool,iel,ic)
          sb_fire_emissions_total_els(ic,iel) = sb_fire_emissions_total_els(ic,iel) + sb_fire_emissions_mt(ipool,iel,ic,1)
          veg_fire_emissions_total_els(ic,iel) = veg_fire_emissions_total_els(ic,iel) + veg_fire_emissions_mt(ipool,iel,ic)
          sb_fire_litter_ash_total_els(ic,iel) = sb_fire_litter_ash_total_els(ic,iel) + sb_fire_litter_ash_mt(ipool,iel,ic,1)
        END DO
      END DO
    END DO
    veg_fire_litterfall_total_els = veg_fire_litterfall_total_els / dtime * 1000000._wp ! to umols s-1
    sb_fire_emissions_total_els = sb_fire_emissions_total_els / dtime * 1000000._wp
    veg_fire_emissions_total_els = veg_fire_emissions_total_els / dtime * 1000000._wp
    sb_fire_litter_ash_total_els = sb_fire_litter_ash_total_els / dtime * 1000000._wp

    IF (dsl4jsb_Config(Q_DFIRE_SF_)%has_diagnostics) THEN
      CALL run_diagnostics_at_pft(nc, sb_pool_mt, sb_fire_emissions_mt, sb_fire_litter_ash_mt)
    ENDIF
  END SUBROUTINE update_q_dfire_sf_impacts


  ! ======================================================================================================= !
  !>
  !> Calculate residue and any mass imbalance
  !>
  SUBROUTINE run_diagnostics_at_pft(nc, sb_pool_mt, sb_fire_emissions_mt, sb_fire_litter_ash_mt)
    INTEGER, INTENT(in) :: nc
    REAL(wp), INTENT(in) :: sb_pool_mt(:,:,:,:)
    REAL(wp), INTENT(in) :: sb_fire_emissions_mt(:,:,:,:)
    REAL(wp), INTENT(in) :: sb_fire_litter_ash_mt(:,:,:,:)

    INTEGER :: ic, ip, iel
    REAL(wp) :: diff
    LOGICAL :: test_is_sb_pool_negative
    LOGICAL :: has_any_test_failed

    test_is_sb_pool_negative = .FALSE.
    has_any_test_failed = .FALSE.

    IF (ANY((sb_pool_mt(:,:,:,df_sf_spq_layer_top_sl) - &
            (sb_fire_emissions_mt(:,:,:,df_sf_spq_layer_top_sl) + &
             sb_fire_litter_ash_mt(:,:,:,df_sf_spq_layer_top_sl))) < 0._wp)) THEN

      DO ic=1, nc
        DO ip=1, SIZE(sb_pool_mt, DIM=1) !ndeadpools
          DO iel=1, SIZE(sb_pool_mt, DIM=2) ! n elements
            diff = sb_pool_mt(ip,iel,ic,df_sf_spq_layer_top_sl) - &
                   (sb_fire_emissions_mt(ip,iel,ic,df_sf_spq_layer_top_sl) + &
                   sb_fire_litter_ash_mt(ip,iel,ic,df_sf_spq_layer_top_sl))
            IF (diff < 0._wp) THEN
              WRITE(*,*) "post.f90:: ic, pool, el, sb(ic,pool,el), ab_fract =", &
                      ic, ip, iel, diff
              WRITE(*,*) "post.f90:: sb_pool, sb_fire_em, sb_fire_ash =", &
                       sb_pool_mt(ip,iel,ic,df_sf_spq_layer_top_sl), &
                       sb_fire_emissions_mt(ip,iel,ic,df_sf_spq_layer_top_sl), &
                       sb_fire_litter_ash_mt(ip,iel,ic,df_sf_spq_layer_top_sl)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      test_is_sb_pool_negative = .TRUE.
      WRITE(*,*) "(post) is too much fuel removed from sb_pool?" , test_is_sb_pool_negative
    ENDIF

    has_any_test_failed = test_is_sb_pool_negative !< add here more tests

    IF (has_any_test_failed) THEN
      CALL finish("post.f90", "A test(s) failed")
    ENDIF

  END SUBROUTINE run_diagnostics_at_pft
  ! ======================================================================================================= !
  !>
  !> Calculate residue and any mass imbalance
  !>
  ELEMENTAL SUBROUTINE sf_calc_residue(fuel_consumption, emitted_elements, residue, mass_overflow)
    REAL(wp), INTENT(in) :: fuel_consumption !< (#nc,#pools,#elements) [molX/m2]
    REAL(wp), INTENT(in) :: emitted_elements !< (#nc,#pools,#elements) [molX/m2]
    REAL(wp), INTENT(out) :: residue         !< (#nc,#pools,#elements) [molX/m2]
    REAL(wp), INTENT(out) :: mass_overflow   !< (#nc,#pools,#elements) [molX/m2]

    !< Compute residue, where Emitted elements < FC
    residue = fuel_consumption - emitted_elements
    mass_overflow = 0._wp

    IF (residue < 0._wp) THEN
      mass_overflow = ABS(residue)
      residue = 0._wp
    END IF
    ! Use for debugging
    !DO ic=1, SIZE(residue, DIM=1) ! 1= nc dimension index
    !  DO ipool=1, SIZE(residue, DIM=2) ! 2=pool dimension index
    !    DO iel=1, SIZE(residue, DIM=3) ! 3=elements dimension index
    !      IF (residue(ic,ipool,iel) < 0._wp .AND. mask(ic) THEN
    !        mass_overflow(ic,ipool,iel) = ABS(residue(ic,ipool,iel))
    !        residue(ic,ipool,iel) = 0._wp
    !      END IF
    !    END DO
    !  END DO
    !END DO
  END SUBROUTINE sf_calc_residue

#endif
END MODULE mo_q_dfire_sf_impacts
