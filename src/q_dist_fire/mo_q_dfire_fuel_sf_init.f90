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
MODULE mo_q_dfire_fuel_sf_init
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message
  USE mo_jsb_control,             ONLY: debug_on

  USE mo_jsb_process_class,       ONLY: Q_DFIRE_SF_
  USE mo_veg_constants,           ONLY: itree, igrass   !< lctlib%growthform

  dsl4jsb_Use_config(Q_DFIRE_SF_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: q_dfire_fuel_sf_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_fuel_sf_init'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Intialize disturbance fire process (after memory has been set up)
  !>
  SUBROUTINE q_dfire_fuel_sf_init(tile)
    USE mo_jsb_class,               ONLY: get_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: Q_DFIRE_FUEL_SF_, Q_DFIRE_SF_
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,         ONLY: t_jsb_model

   USE mo_q_dfire_sf_config_class,  ONLY: Q_DFIRE_ROS_THONICKE, Q_DFIRE_ROS_OBERHAGEMANN

    dsl4jsb_Use_memory(Q_DFIRE_FUEL_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib !< land-cover-type library - parameter across pfts
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':q_dfire_fuel_sf_init'
    ! ----------------------------------------------------------------------------------------------------- !
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(Q_DFIRE_SF_)
    dsl4jsb_Def_memory(Q_DFIRE_FUEL_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real3D_onDomain      :: dead_fuel_1hr_dfp
    dsl4jsb_Real3D_onDomain      :: dead_fuel_10hr_dfp
    dsl4jsb_Real3D_onDomain      :: dead_fuel_100hr_dfp
    dsl4jsb_Real3D_onDomain      :: dead_fuel_1000hr_dfp

    dsl4jsb_Real3D_onDomain      :: live_grass_fuel_1hr_lfp
    dsl4jsb_Real3D_onDomain      :: cured_live_grass_fuel_1hr_dfp

    dsl4jsb_Real3D_onDomain      :: live_tree_fuel_1hr_lfp
    dsl4jsb_Real3D_onDomain      :: live_tree_fuel_10hr_lfp
    dsl4jsb_Real3D_onDomain      :: live_tree_fuel_100hr_lfp
    dsl4jsb_Real3D_onDomain      :: live_tree_fuel_1000hr_lfp

    dsl4jsb_Real2D_onDomain      :: forest_fract
    dsl4jsb_Real2D_onDomain      :: grass_fract

    dsl4jsb_Real2D_onDomain      :: fuel_bulk_density
    dsl4jsb_Real2D_onDomain      :: moisture_extinction
    dsl4jsb_Real2D_onDomain      :: fwd_wind_factor

    dsl4jsb_Real2D_onDomain      :: pm_fire
    dsl4jsb_Real2D_onDomain      :: pm_fire_output
    dsl4jsb_Real2D_onDomain      :: pm_prop_ck
    dsl4jsb_Real2D_onDomain      :: pm_tau
    dsl4jsb_Real2D_onDomain      :: tau_l
    dsl4jsb_Real2D_onDomain      :: tau_c
    dsl4jsb_Real2D_onDomain      :: par1_diag
    dsl4jsb_Real2D_onDomain      :: par2_diag
    dsl4jsb_Real2D_onDomain      :: bark_thickness
    dsl4jsb_Real2D_onDomain      :: DBH_diag

    dsl4jsb_Real3D_onDomain      :: fuel_consumption_els
    dsl4jsb_Real3D_onDomain      :: emissions_spc
    dsl4jsb_Real3D_onDomain      :: emissions_els
    dsl4jsb_Real3D_onDomain      :: residue_els
    dsl4jsb_Real3D_onDomain      :: overflow_els

    dsl4jsb_Real2D_onDomain      :: ck
    dsl4jsb_Real2D_onDomain      :: r_ck
    dsl4jsb_Real2D_onDomain      :: ck_p
    dsl4jsb_Real2D_onDomain      :: scorch_height
    dsl4jsb_Real2D_onDomain      :: F_scorch_height
    dsl4jsb_Real2D_onDomain      :: crown_length
    dsl4jsb_Real2D_onDomain      :: tree_height
    dsl4jsb_Real2D_onDomain      :: lvg_moisture

    dsl4jsb_Real3D_onDomain      :: veg_fire_litterfall_total_els
    dsl4jsb_Real3D_onDomain      :: veg_fire_emissions_total_els
    dsl4jsb_Real3D_onDomain      :: sb_fire_emissions_total_els
    dsl4jsb_Real3D_onDomain      :: sb_fire_litter_ash_total_els
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_active(Q_DFIRE_FUEL_SF_)) RETURN
    IF (tile%lcts(1)%lib_id == 0) RETURN                !< run this init only if the present tile is a pft
    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model  => get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_DFIRE_SF_)
    dsl4jsb_Get_memory(Q_DFIRE_FUEL_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, dead_fuel_1hr_dfp)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, dead_fuel_10hr_dfp)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, dead_fuel_100hr_dfp)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, dead_fuel_1000hr_dfp)

    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, live_grass_fuel_1hr_lfp)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, cured_live_grass_fuel_1hr_dfp)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, live_tree_fuel_1hr_lfp)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, live_tree_fuel_10hr_lfp)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, live_tree_fuel_100hr_lfp)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, live_tree_fuel_1000hr_lfp)

    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, forest_fract)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, grass_fract)

    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, fuel_bulk_density)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, moisture_extinction)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, fwd_wind_factor)

    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, pm_fire)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, pm_fire_output)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, pm_prop_ck)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, pm_tau)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, tau_l)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, tau_c)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, par1_diag)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, par2_diag)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, bark_thickness)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, DBH_diag)

    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, ck)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, r_ck)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, ck_p)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, scorch_height)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, F_scorch_height)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, crown_length)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, tree_height)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_FUEL_SF_, lvg_moisture)

    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, fuel_consumption_els)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, emissions_spc)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, emissions_els)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, residue_els)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, overflow_els)

    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, veg_fire_litterfall_total_els)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, veg_fire_emissions_total_els)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, sb_fire_emissions_total_els)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_FUEL_SF_, sb_fire_litter_ash_total_els)
    ! ----------------------------------------------------------------------------------------------------- !
    dead_fuel_1hr_dfp = 0._wp
    dead_fuel_10hr_dfp = 0._wp
    dead_fuel_100hr_dfp = 0._wp
    dead_fuel_1000hr_dfp = 0._wp

    live_grass_fuel_1hr_lfp = 0._wp
    cured_live_grass_fuel_1hr_dfp = 0._wp
    live_tree_fuel_1hr_lfp = 0._wp
    live_tree_fuel_10hr_lfp = 0._wp
    live_tree_fuel_100hr_lfp = 0._wp
    live_tree_fuel_1000hr_lfp = 0._wp

    forest_fract = 0._wp
    grass_fract = 0._wp

    fuel_bulk_density = 0._wp

    pm_fire = 0._wp
    pm_fire_output = 0._wp
    pm_prop_ck = 0._wp
    pm_tau = 0._wp
    tau_l = 0._wp
    tau_c = 0._wp
    par1_diag = 0._wp
    par2_diag = 0._wp
    bark_thickness = 0._wp
    DBH_diag = 0._wp

    ck = 0._wp
    r_ck = 0._wp
    ck_p = 0._wp
    scorch_height = 0._wp
    F_scorch_height = 0._wp
    crown_length = 0._wp
    tree_height = 0._wp
    lvg_moisture = 0._wp

    fuel_consumption_els = 0._wp
    emissions_spc = 0._wp
    emissions_els = 0._wp
    residue_els = 0._wp
    overflow_els = 0._wp

    veg_fire_litterfall_total_els = 0._wp
    veg_fire_emissions_total_els = 0._wp
    sb_fire_emissions_total_els = 0._wp
    sb_fire_litter_ash_total_els = 0._wp

    !< All PFTs values are loaded, it is useful for aggregation to account for all of them
    !< Initialization at the process only loads those values that are specific to the PFT,
    !< so the aggregation could be less than expected
    IF (dsl4jsb_Config(Q_DFIRE_SF_)%ros_model_id == Q_DFIRE_ROS_THONICKE) THEN
      moisture_extinction = dsl4jsb_Lctlib_param(q_moisture_extinction_thon)
    ELSE IF (dsl4jsb_Config(Q_DFIRE_SF_)%ros_model_id == Q_DFIRE_ROS_OBERHAGEMANN) THEN
      moisture_extinction = dsl4jsb_Lctlib_param(q_moisture_extinction_ober)
    END IF
    fwd_wind_factor = dsl4jsb_Lctlib_param(q_fwd_wind_factor)

    ! Mark process as successfully initialized, MR 428
    !CALL model%processes(Q_DFIRE_FUEL_SF_)%p%Set_init_status_true
  END SUBROUTINE q_dfire_fuel_sf_init

#endif
END MODULE mo_q_dfire_fuel_sf_init
