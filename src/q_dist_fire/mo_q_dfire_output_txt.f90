!> Extended functions to complement QUINCY model output as ascii text for dfire
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
MODULE mo_q_dfire_output_txt

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message_text, message, finish
  USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
  USE mo_qs_forcing,            ONLY: forcing_options
  USE mo_q_dfire_sf_fuel_class, ONLY: t_sf_fuel_class
  USE mo_quincy_model_config,   ONLY: QLAND, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION
!  USE mo_qs_output_txt_util

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,        ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,  ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_GROWTH_ID, &
    &                                   VEG_BGCM_LITTERFALL_ID, VEG_BGCM_FRAC_ALLOC_ID,              &
    &                                   VEG_BGCM_EXUDATION_ID, VEG_BGCM_ESTABLISHMENT_ID, VEG_BGCM_RESERVE_USE_ID, &
    &                                   VEG_BGCM_FIRE_EMISSIONS_ID,  VEG_BGCM_FIRE_LITTERFALL_ID
  USE mo_lnd_bgcm_store_class,  ONLY: SB_BGCM_POOL_ID, SB_BGCM_DELTA_POOLS_ID, SB_BGCM_FORMATION_ID, SB_BGCM_LOSS_ID, &
    &                                   SB_BGCM_TRANSPORT_ID, SB_BGCM_MYCO_EXPORT_ID, &
    &                                   SB_BGCM_FIRE_LITTER_ASH_ID, SB_BGCM_FIRE_EMISSIONS_ID

  USE mo_q_dfire_constants,   ONLY: sf_fuel_live_fuel_nr, df_sf_species_nr
  USE mo_q_dfire_constants,   ONLY: df_sf_fc_full_nr, df_sf_fc_fast_nr,  &
                                      df_sf_fc_1hr_idx, df_sf_fc_10hr_idx, df_sf_fc_spread_nr, &
                                      df_sf_fc_100hr_idx, df_sf_fc_1000hr_idx
  USE mo_q_dfire_constants,   ONLY: df_sf_species_co2_idx, df_sf_species_co_idx, &
                                      df_sf_species_ch4_idx,&
                                      !df_sf_species_voc_idx, &
                                      !df_sf_species_tpm_idx, df_sf_species_nrox_idx, &
                                      df_sf_species_nh3_idx, df_sf_species_n2o_idx, &
                                      df_sf_species_n2_idx
  USE mo_q_dfire_constants,   ONLY: df_sf_omega_o_idx, df_sf_omega_nl_idx
  USE mo_q_dfire_constants,   ONLY: df_sf_elements_nr
  USE mo_jsb_process_class,   ONLY: SB_, Q_DFIRE_, Q_DFIRE_FUEL_

  USE mo_linked_list,             ONLY: t_q_linkedlist, t_q_linkedlist_node
  USE mo_q_dfire_output_txt_util, ONLY: t_qs_output_file, t_qs_output_2d, t_qs_output_3d, add_to_list, t_qs_soutputs
  USE mo_q_dfire_output_txt_util, ONLY: EN_OP_AVERAGE, EN_PRE_OP_SUM

  dsl4jsb_Use_config(SB_)

  IMPLICIT NONE

  PUBLIC  ! all is public

  INTEGER, PARAMETER, PRIVATE   :: ic = 1                             !< the point in ics:ice of which data will be written to file; per default = 1

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_dfire_output_txt'

  TYPE(t_qs_soutputs), SAVE :: sout !< simple output library
  !TYPE(t_q_linkedlist), SAVE :: output_files
  !< params
  CHARACTER(:), SAVE, ALLOCATABLE :: spinup_id_param                !< add ID to filename if file contains data from spinup simulation
  CHARACTER(:), SAVE, ALLOCATABLE :: output_interval_flux_str_param !< output interval for "flux" variables
CONTAINS

  ! ======================================================================================================= !
  !>
  !> Provide dependencies used accros its methods
  !>
  SUBROUTINE q_dfire_setparams_outputs(spinup_id, output_interval_flux_str)
    CHARACTER(*), INTENT(in) :: spinup_id
    CHARACTER(*), INTENT(in) :: output_interval_flux_str
    spinup_id_param = spinup_id
    output_interval_flux_str_param = output_interval_flux_str
  END SUBROUTINE q_dfire_setparams_outputs

  ! ======================================================================================================= !
  !>
  !> Write data to outputs
  !>
  SUBROUTINE q_dfire_init_outputs(tile, options)
    USE mo_jsb_class,              ONLY: Get_model
    USE mo_jsb_model_class,        ONLY: t_jsb_model
    USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,         ONLY: t_jsb_task_options
    USE mo_jsb_process_class,      ONLY: Q_DFIRE_, Q_DFIRE_FUEL_
    USE mo_jsb_grid_class,         ONLY: t_jsb_vgrid
    USE mo_jsb_grid,               ONLY: Get_vgrid
    USE mo_jsb_math_constants,     ONLY: eps8
    USE mo_jsb_physical_constants, ONLY: molar_mass_C, molar_mass_N
    USE mo_sb_constants,           ONLY: organic_cn
    USE mo_isotope_util,           ONLY: calc_delta_C13, calc_delta_C14, calc_delta_N15
    USE mo_spq_config_class,       ONLY: max_soil_layers
    USE mo_phy_schemes,            ONLY: qsat_mixed
    USE mo_atmland_constants,      ONLY: eps_vpd
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(Q_DFIRE_)
    dsl4jsb_Use_memory(Q_DFIRE_FUEL_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(in)  :: tile                     !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)  :: options                  !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER   :: model                        !< the model
    TYPE(t_lnd_bgcm_store), POINTER   :: bgcm_store                   !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER   :: lctlib                       !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid), POINTER        :: vgrid_soil_w                 !< Vertical grid
    INTEGER                           :: nsoil_w                      !< number of soil layers as used/defined by the SB_ process
    LOGICAL                           :: set_zero_flux_acc_var_local  !< local variable of input variable: set_zero_flux_acc_variables
    INTEGER                           :: iblk, ics, ice, nc, dtime
    CHARACTER(:), ALLOCATABLE         :: fname
    CHARACTER(:), ALLOCATABLE         :: units
    INTEGER :: fid_em, fid_weather, fid_ems
    INTEGER :: fid_veg_fire_fall, fid_veg_fire_ems
    INTEGER :: fid_sb_fire_ash, fid_sb_fire_ems
    INTEGER :: fid_fc, fid_ofw, fid_res, fid_fire
    INTEGER :: fid_livefuel, fid_deadfuel
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_fire_emissions_mt
    dsl4jsb_Def_mt2L2D :: veg_fire_litterfall_mt
    dsl4jsb_Def_mt2L3D :: sb_fire_emissions_mt
    dsl4jsb_Def_mt2L3D :: sb_fire_litter_ash_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_DFIRE_)
    dsl4jsb_Def_memory(Q_DFIRE_FUEL_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_DFIRE_
    dsl4jsb_Real2D_onChunk      :: nesterov_index
    dsl4jsb_Real2D_onChunk      :: t_air_min
    dsl4jsb_Real2D_onChunk      :: t_air_max
    dsl4jsb_Real2D_onChunk      :: precip_daily_mmd

    dsl4jsb_Real2D_onChunk      :: area_burned_output
    dsl4jsb_Real2D_onChunk      :: area_burned_fract_output
    dsl4jsb_Real2D_onChunk      :: fdi
    dsl4jsb_Real2D_onChunk      :: mean_fire_area
    dsl4jsb_Real2D_onChunk      :: e_n_ignitions

    ! Q_DFIRE_FUEL_
    dsl4jsb_Real3D_onChunk      :: emissions_spc
    dsl4jsb_Real3D_onChunk      :: emissions_els
    dsl4jsb_Real3D_onChunk      :: residue_els
    dsl4jsb_Real3D_onChunk      :: fuel_consumption_els
    dsl4jsb_Real3D_onChunk      :: overflow_els

    dsl4jsb_Real3D_onChunk      :: live_grass_fuel_1hr_lfp
    dsl4jsb_Real3D_onChunk      :: live_tree_fuel_1hr_lfp
    dsl4jsb_Real3D_onChunk      :: live_tree_fuel_10hr_lfp
    dsl4jsb_Real3D_onChunk      :: live_tree_fuel_100hr_lfp
    dsl4jsb_Real3D_onChunk      :: live_tree_fuel_1000hr_lfp

    dsl4jsb_Real3D_onChunk      :: dead_fuel_1hr_dfp
    dsl4jsb_Real3D_onChunk      :: dead_fuel_10hr_dfp
    dsl4jsb_Real3D_onChunk      :: dead_fuel_100hr_dfp
    dsl4jsb_Real3D_onChunk      :: dead_fuel_1000hr_dfp
    ! ----------------------------------------------------------------------------------------------------- !
    iblk     = options%iblk
    ics      = options%ics
    ice      = options%ice
    nc       = options%nc
    dtime    = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       =  vgrid_soil_w%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_DFIRE_)
    dsl4jsb_Get_memory(Q_DFIRE_FUEL_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_FIRE_EMISSIONS_ID, veg_fire_emissions_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_FIRE_LITTERFALL_ID, veg_fire_litterfall_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_FIRE_EMISSIONS_ID, sb_fire_emissions_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_FIRE_LITTER_ASH_ID, sb_fire_litter_ash_mt)
    ! ----------------------------------------------------------------------------------------------------- !

    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_DFIRE_
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, nesterov_index)
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, t_air_min)
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, t_air_max)
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, precip_daily_mmd)

    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, area_burned_output)
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, area_burned_fract_output)
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, fdi)
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, mean_fire_area)
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, e_n_ignitions)

    ! DFIRE_FUEL_
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, emissions_spc)
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, emissions_els)
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, residue_els)
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, fuel_consumption_els)
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, overflow_els)

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_grass_fuel_1hr_lfp)
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1hr_lfp)
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_10hr_lfp)
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_100hr_lfp)
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1000hr_lfp)

    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, dead_fuel_1hr_dfp)
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, dead_fuel_10hr_dfp)
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, dead_fuel_100hr_dfp)
    dsl4jsb_Get_var3D_onChunk(Q_DFIRE_FUEL_, dead_fuel_1000hr_dfp)
    ! ----------------------------------------------------------------------------------------------------- !
    sout = t_qs_soutputs()

    fname = 'dfire_emissions_spc'//'_'//TRIM(spinup_id_param)//TRIM(output_interval_flux_str_param)//'.txt'
    fid_em = sout % add_file(fname)
    ! todo: * units_mols_to_micromols / dtime) ! to flux
    units="g_X_m-2_ts-1" ! format with _ so R script can properly parse names
    CALL sout % add_var(fid_em, "CO["//units//"]", units, op=EN_OP_AVERAGE, &
                                  pdata2d=emissions_spc(ic,df_sf_species_co_idx))
    CALL sout % add_var(fid_em, "CO2["//units//"]", units, op=EN_OP_AVERAGE, &
                                  pdata2d=emissions_spc(ic,df_sf_species_co2_idx))
    CALL sout % add_var(fid_em, "CH4["//units//"]", units, op=EN_OP_AVERAGE, &
                                  pdata2d=emissions_spc(ic,df_sf_species_ch4_idx))
    CALL sout % add_var(fid_em, "NH3["//units//"]", units, op=EN_OP_AVERAGE, &
                                  pdata2d=emissions_spc(ic,df_sf_species_nh3_idx))
    CALL sout % add_var(fid_em, "N2O["//units//"]", units, op=EN_OP_AVERAGE, &
                                  pdata2d=emissions_spc(ic,df_sf_species_n2o_idx))
    CALL sout % add_var(fid_em, "N2["//units//"]", units, op=EN_OP_AVERAGE, &
                                  pdata2d=emissions_spc(ic,df_sf_species_n2_idx))

    fname = 'dfire_weather'//'_'//TRIM(spinup_id_param)//TRIM(output_interval_flux_str_param)//'.txt'
    fid_weather = sout % add_file(fname)
    CALL sout % add_var(fid_weather, "nesterov_index[C2]", "C2", op=EN_OP_AVERAGE, pdata2d=nesterov_index(ic))
    CALL sout % add_var(fid_weather, "t_air_min[C]", "C", op=EN_OP_AVERAGE, pdata2d=t_air_min(ic))
    CALL sout % add_var(fid_weather, "t_air_max[C]", "C", op=EN_OP_AVERAGE, pdata2d=t_air_max(ic))
    CALL sout % add_var(fid_weather, "precip_daily_mmd[mm_d-1]", "mm_d-1", op=EN_OP_AVERAGE, &
                                pdata2d=precip_daily_mmd(ic))

    ! dead fuels
    units = 'gDM_m-2'
    fname = 'dfire_deadfuel'//'_'//TRIM(spinup_id_param)//TRIM(output_interval_flux_str_param)//'.txt'
    fid_deadfuel = sout % add_file(fname)
    CALL sout % add_var(fid_deadfuel, "1hr["//units//"]", units, op=EN_OP_AVERAGE, &
                                pdata3d=dead_fuel_1hr_dfp(ic,:), preop=EN_PRE_OP_SUM)
    CALL sout % add_var(fid_deadfuel, "10hr["//units//"]", units, op=EN_OP_AVERAGE, &
                                pdata3d=dead_fuel_10hr_dfp(ic,:), preop=EN_PRE_OP_SUM)
    CALL sout % add_var(fid_deadfuel, "100hr["//units//"]", units, op=EN_OP_AVERAGE, &
                                pdata3d=dead_fuel_100hr_dfp(ic,:), preop=EN_PRE_OP_SUM)
    CALL sout % add_var(fid_deadfuel, "1000hr["//units//"]", units, op=EN_OP_AVERAGE, &
                                pdata3d=dead_fuel_1000hr_dfp(ic,:), preop=EN_PRE_OP_SUM)

    ! live fuels
    units = 'gDM_m-2'
    fname = 'dfire_livefuel'//'_'//TRIM(spinup_id_param)//TRIM(output_interval_flux_str_param)//'.txt'
    fid_livefuel = sout % add_file(fname)
    CALL sout % add_var(fid_livefuel, "grass_1hr["//units//"]", units, op=EN_OP_AVERAGE, &
                                  pdata3d=live_grass_fuel_1hr_lfp(ic,:), preop=EN_PRE_OP_SUM)
    CALL sout % add_var(fid_livefuel, "tree_1hr["//units//"]", units, op=EN_OP_AVERAGE, &
                                  pdata3d=live_tree_fuel_1hr_lfp(ic,:), preop=EN_PRE_OP_SUM)
    CALL sout % add_var(fid_livefuel, "tree_10hr["//units//"]", units, op=EN_OP_AVERAGE, &
                                  pdata3d=live_tree_fuel_10hr_lfp(ic,:), preop=EN_PRE_OP_SUM)
    CALL sout % add_var(fid_livefuel, "tree_100hr["//units//"]", units, op=EN_OP_AVERAGE, &
                                  pdata3d=live_tree_fuel_100hr_lfp(ic,:), preop=EN_PRE_OP_SUM)
    CALL sout % add_var(fid_livefuel, "tree_1000hr["//units//"]", units, op=EN_OP_AVERAGE, &
                                  pdata3d=live_tree_fuel_1000hr_lfp(ic,:), preop=EN_PRE_OP_SUM)


    ! diagnostics
    units = 'umols_m-2_s-1'
    fname = 'dfire_emissions_els'//'_'//TRIM(spinup_id_param)//TRIM(output_interval_flux_str_param)//'.txt'
    CALL create_output_txt_as_var_els(sout, fname, units, emissions_els)
    fname = 'dfire_residue_els'//'_'//TRIM(spinup_id_param)//TRIM(output_interval_flux_str_param)//'.txt'
    CALL create_output_txt_as_var_els(sout, fname, units, residue_els)
    fname = 'dfire_fuel_consumption_els'//'_'//TRIM(spinup_id_param)//TRIM(output_interval_flux_str_param)//'.txt'
    CALL create_output_txt_as_var_els(sout, fname, units, fuel_consumption_els)
    fname = 'dfire_overflow_els'//'_'//TRIM(spinup_id_param)//TRIM(output_interval_flux_str_param)//'.txt'
    CALL create_output_txt_as_var_els(sout, fname, units, overflow_els)

    ! fluxes
    units = 'umols_m-2_s-1'
    fname = 'dfire_veg_fire_litterfall_flux'//'_'//TRIM(spinup_id_param)//TRIM(output_interval_flux_str_param)//'.txt'
    CALL create_output_txt_as_flux(sout, fname, units, veg_fire_litterfall_mt)

    fname = 'dfire_veg_fire_emissions_flux'//'_'//TRIM(spinup_id_param)//TRIM(output_interval_flux_str_param)//'.txt'
    CALL create_output_txt_as_flux(sout, fname, units, veg_fire_emissions_mt)

    fname = 'dfire_sb_fire_emissions_flux'//'_'//TRIM(spinup_id_param)//TRIM(output_interval_flux_str_param)//'.txt'
    CALL create_output_txt_as_flux(sout, fname, units, sb_fire_emissions_mt(:,:,:,1)) !< 1 = soil top layer

    fname = 'dfire_sb_fire_ash_flux'//'_'//TRIM(spinup_id_param)//TRIM(output_interval_flux_str_param)//'.txt'
    CALL create_output_txt_as_flux(sout, fname, units, sb_fire_litter_ash_mt(:,:,:,1)) !< = top soil layer

    ! main outputs
    fname = 'dfire_areaburned'//'_'//TRIM(spinup_id_param)//TRIM(output_interval_flux_str_param)//'.txt'
    fid_fire = sout % add_file(fname)
    CALL sout % add_var(fid_fire, "area_burned[fract]", "fract", op=EN_OP_AVERAGE, pdata2d=area_burned_output(ic))
    CALL sout % add_var(fid_fire, "area_burned_fract[fract]", "fract", op=EN_OP_AVERAGE, pdata2d=area_burned_fract_output(ic))
    CALL sout % add_var(fid_fire, "fdi[-]", "-", op=EN_OP_AVERAGE, pdata2d=fdi(ic))
    CALL sout % add_var(fid_fire, "mean_fire_area[ha]", "ha", op=EN_OP_AVERAGE, pdata2d=mean_fire_area(ic))
    CALL sout % add_var(fid_fire, "e_n_ignitions[ha_ts-1]", "ha_ts-1", op=EN_OP_AVERAGE, pdata2d=e_n_ignitions(ic))

    CALL sout % open()
    CALL sout % init()
  CONTAINS
    SUBROUTINE initAllFiles(node)
      TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

      SELECT TYPE(v=>node%value)
      CLASS IS (t_qs_output_file)
        CALL v % init()
      CLASS DEFAULT
        ! unexpected type found, raise error
        CALL finish("initAllFiles", "Unexpected class found")
      END SELECT
    END SUBROUTINE initAllFiles

    SUBROUTINE openAllFiles(node)
      TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

      SELECT TYPE(v=>node%value)
      CLASS IS (t_qs_output_file)
        CALL v % open()
      CLASS DEFAULT
        ! unexpected type found, raise error
        CALL finish("openAllFiles", "Unexpected class found")
      END SELECT
    END SUBROUTINE openAllFiles
  END SUBROUTINE q_dfire_init_outputs

  ! ======================================================================================================= !
  !>
  !> Add new output for a flux variable
  !>
  SUBROUTINE create_output_txt_as_flux(sout, name, units, var)
    TYPE(t_qs_soutputs), INTENT(inout) :: sout
    CHARACTER(*), INTENT(in)        :: name, units
    REAL(wp), POINTER, INTENT(in)   :: var(:,:,:)  !< Quincy flux type variable (#pool, #elems, #nc)
    INTEGER :: fid

    fid = sout % add_file(name)

    !< Sum pools when writting data
    CALL sout % add_var(fid, "C["//units//"]", units, op=EN_OP_AVERAGE, &
                pdata3d=var(:,ixC,ic), preop=EN_PRE_OP_SUM)
    CALL sout % add_var(fid, "C13["//units//"]", units, op=EN_OP_AVERAGE, &
                pdata3d=var(:,ixC13,ic), preop=EN_PRE_OP_SUM)
    CALL sout % add_var(fid, "C14["//units//"]", units, op=EN_OP_AVERAGE, &
                pdata3d=var(:,ixC14,ic), preop=EN_PRE_OP_SUM)
    CALL sout % add_var(fid, "N["//units//"]", units, op=EN_OP_AVERAGE, &
                pdata3d=var(:,ixN,ic), preop=EN_PRE_OP_SUM)
    CALL sout % add_var(fid, "N15["//units//"]", units, op=EN_OP_AVERAGE, &
                pdata3d=var(:,ixN15,ic), preop=EN_PRE_OP_SUM)
    CALL sout % add_var(fid, "P["//units//"]", units, op=EN_OP_AVERAGE, &
                pdata3d=var(:,ixP,ic), preop=EN_PRE_OP_SUM)
  END SUBROUTINE create_output_txt_as_flux

  ! ======================================================================================================= !
  !>
  !> Add new output and print each element within (C,C13,C14,N,N15,P)
  !>
  SUBROUTINE create_output_txt_as_var_els(sout, name, units, var)
    TYPE(t_qs_soutputs), INTENT(inout) :: sout
    CHARACTER(*), INTENT(in)        :: name, units
    REAL(wp), POINTER, INTENT(in)   :: var(:,:)  !< Quincy flux type variable (#nc,#elements)
    INTEGER :: fid

    fid = sout % add_file(name)

    !< Output for each element
    CALL sout % add_var(fid, "C["//units//"]", units, op=EN_OP_AVERAGE, pdata2d=var(ic,ixC))
    CALL sout % add_var(fid, "C13["//units//"]", units, op=EN_OP_AVERAGE, pdata2d=var(ic,ixC13))
    CALL sout % add_var(fid, "C14["//units//"]", units, op=EN_OP_AVERAGE, pdata2d=var(ic,ixC14))
    CALL sout % add_var(fid, "N["//units//"]", units, op=EN_OP_AVERAGE, pdata2d=var(ic,ixN))
    CALL sout % add_var(fid, "N15["//units//"]", units, op=EN_OP_AVERAGE, pdata2d=var(ic,ixN15))
    CALL sout % add_var(fid, "P["//units//"]", units, op=EN_OP_AVERAGE, pdata2d=var(ic,ixP))
  END SUBROUTINE create_output_txt_as_var_els


  ! ======================================================================================================= !
  !>
  !> Update all files with new data
  !>
  SUBROUTINE q_dfire_update_outputs()
    CALL sout % update()
    !CALL output_files % traverse(runUpdates)
  CONTAINS
    SUBROUTINE runUpdates(node)
      TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

      SELECT TYPE(v=>node%value)
      CLASS IS (t_qs_output_file)
        CALL v % update()
      CLASS DEFAULT
        ! unexpected type found, raise error
        CALL finish("runUpdates", "Unexpected class found")
      END SELECT
    END SUBROUTINE runUpdates
  END SUBROUTINE q_dfire_update_outputs

  ! ======================================================================================================= !
  !>
  !> Write data to outputs
  !>
  SUBROUTINE q_dfire_write_outputs(time_d)
    REAL(wp), INTENT(in) :: time_d !< Time in days
    CALL sout % write(time_d)
    !CALL output_files % traverse(writeAllFiles)
  CONTAINS
    SUBROUTINE writeAllFiles(node)
      TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

      SELECT TYPE(v=>node%value)
      CLASS IS (t_qs_output_file)
        CALL v % write(time_d)
      CLASS DEFAULT
        ! unexpected type found, raise error
        CALL finish("writeAllFiles", "Unexpected class found")
      END SELECT
    END SUBROUTINE writeAllFiles
  END SUBROUTINE q_dfire_write_outputs

  ! ======================================================================================================= !
  !>
  !> Reset all outputs (set to 0 all counter on average variables)
  !>
  SUBROUTINE q_dfire_reset_outputs()
    CALL sout % reset()
    !CALL output_files % traverse(resetAllFiles)
  CONTAINS
    SUBROUTINE resetAllFiles(node)
      TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

      SELECT TYPE(v=>node%value)
      CLASS IS (t_qs_output_file)
        CALL v % reset()
      CLASS DEFAULT
        ! unexpected type found, raise error
        CALL finish("resetAllFiles", "Unexpected class found")
      END SELECT
    END SUBROUTINE resetAllFiles
  END SUBROUTINE q_dfire_reset_outputs

  ! ======================================================================================================= !
  !!
  !! Close output files
  !!
  ! ======================================================================================================= !


  ! ======================================================================================================= !
  !>Pools and Fluxes: close files
  !>
  !>  at this point of the simulation all data has been written
  !>
  SUBROUTINE q_dfire_close_outputs
    CALL sout % close()
    !CALL output_files%traverse(closeAllFiles)
  CONTAINS
    ! ======================================================================================================= !
    !>
    !> Execute runOperation subroutine for each element
    !>
    SUBROUTINE closeAllFiles(node)
      TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

      SELECT TYPE(v=>node%value)
      CLASS IS (t_qs_output_file)
        CALL v % close()
      CLASS DEFAULT
        ! unexpected type found, raise error
        CALL finish("closeAllFiles", "Unexpected class found")
      END SELECT
    END SUBROUTINE closeAllFiles
  END SUBROUTINE q_dfire_close_outputs

END MODULE mo_q_dfire_output_txt

