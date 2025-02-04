!> QUINCY 'jena soil model' calculate soil-biogeochemical pools
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!> Compute spitfire thonicke effects and feedbacks
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_dfire_sf_post
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_lnd_bgcm_store,        ONLY: t_lnd_bgcm_store
  USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract

  USE mo_veg_constants,         ONLY: itree, igrass
  USE mo_dfire_sf_adaptor,      ONLY: sf_convert_units
  USE mo_jsb_grid,              ONLY: Get_grid
  USE mo_dfire_constants,       ONLY: df_sf_fc_spread_nr, df_sf_fc_full_nr, df_sf_fc_1hr_idx, df_sf_fc_100hr_idx
  USE mo_dfire_sf_damage,       ONLY: t_sf_prob_mort, calc_prob_mortality
  USE mo_dfire_util,            ONLY: real2logical
!  USE mo_lnd_bgcm_store_class,  ONLY: SB_BGCM_POOL_ID, VEG_BGCM_POOL_ID

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_dfire_sf_post

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_sf_post'

CONTAINS

  ! ======================================================================================================= !
  !>
  !>
  !>
  SUBROUTINE update_dfire_sf_post(tile, options)
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_grid_class,          ONLY: t_jsb_grid
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: VEG_, SPQ_, DFIRE_FUEL_, DFIRE_
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_math_constants,      ONLY: one_day

    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(DFIRE_FUEL_)
    dsl4jsb_Use_memory(DFIRE_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_grid),       POINTER       :: grid
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lnd_bgcm_store), POINTER       :: bgcm_store !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc, ic
    INTEGER                               :: no_children
    REAL(wp)                              :: dtime      !< timestep

    ! add here local variables
    REAL(wp), POINTER                     :: area(:), tile_fract(:)
    REAL(wp), ALLOCATABLE                 :: dead_fc_1hr_tile(:,:) !< 1 hour dead fuel consumption for current tile (#nc, #dead fuel pools) [g DM m-2]
    REAL(wp), ALLOCATABLE                 :: live_fc_1hr_tile(:,:) !< 1 hour live fuel consumptino for current tile (#nc, #live fuel pools) [g DM m-2]

    REAL(wp)                              :: dbh(options%nc) !< Tree diameter breast at height [cm]
    REAL(wp)                              :: lctlib_growthform(options%nc)

    REAL(wp)                              :: w_oi(options%nc, df_sf_fc_full_nr) !< Dead fuel by fuel classes [g DM m-2]
    LOGICAL                               :: is_fire_intense_enough(options%nc) !< TRUE=1, FALSE=0

    TYPE(t_sf_prob_mort)                  :: prob_mort     !< Tree mortality

    CHARACTER(len=*), PARAMETER           :: routine = TRIM(modname)//':update_fdire_sf_fuel'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real2D_onChunk                :: area_burned               ! in [ha ts-1]
    dsl4jsb_Real2D_onChunk                :: area_burned_fract         ! in [-]
    dsl4jsb_Real3D_onChunk                :: dead_frac_consumption     ! in [-]
    dsl4jsb_Real2D_onChunk                :: i_surface                 ! in [kW m-1]
    dsl4jsb_Real2D_onChunk                :: mask_is_fire_intense_enough ! in [-]

    dsl4jsb_Real3D_onChunk                :: dead_fuel_consumption_1hr ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_fuel_consumption_1hr ! in [g DM m-2]

    dsl4jsb_Real2D_onChunk                :: diameter
    dsl4jsb_Real2D_onChunk                :: height
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(DFIRE_FUEL_)
    dsl4jsb_Def_memory(DFIRE_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Local vars here
    ! ----------------------------------------------------------------------------------------------------- !
!    dsl4jsb_Real3D_onChunk      :: soil_depth_sl      ! in [m]

    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(DFIRE_FUEL_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id) 
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    ! ----------------------------------------------------------------------------------------------------- !
    grid => Get_grid(model%grid_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(DFIRE_FUEL_)
    dsl4jsb_Get_memory(DFIRE_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
!    dsl4jsb_Get_var3D_onChunk(SPQ_, soil_depth_sl)              ! in

    dsl4jsb_Get_var2D_onChunk(DFIRE_, area_burned)               ! in [ha ts-1]
    dsl4jsb_Get_var2D_onChunk(DFIRE_, area_burned_fract)         ! in [-]
    dsl4jsb_Get_var3D_onChunk(DFIRE_, dead_frac_consumption)     ! in [-]
    dsl4jsb_Get_var2D_onChunk(DFIRE_, i_surface)                 ! in [kW m-1]
    dsl4jsb_Get_var2D_onChunk(DFIRE_, mask_is_fire_intense_enough)! in [-]

    dsl4jsb_Get_var2D_onChunk(VEG_, diameter)                    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, height)                      ! in

    dsl4jsb_Get_var3D_onChunk(DFIRE_, dead_fuel_consumption_1hr) ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_, live_fuel_consumption_1hr) ! in [g DM m-2]

    ! ----------------------------------------------------------------------------------------------------- !
!    bgcm_store => tile%bgcm_store
!    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
!    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(dead_fc_1hr_tile, mold=dead_fuel_consumption_1hr)
    ALLOCATE(live_fc_1hr_tile, mold=live_fuel_consumption_1hr)

    ! init here
    area => grid %area (ics:ice, iblk)
    tile_fract => tile%fract(ics:ice, iblk)
    is_fire_intense_enough = .FALSE.

    lctlib_growthform(:) = lctlib%growthform
    prob_mort = t_sf_prob_mort(nc)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !

    dbh(:) = sf_convert_units("m", "cm", diameter) !< dbh to spitfire units

    ! calc tree mortality
    IF (lctlib%growthform == itree) THEN
      is_fire_intense_enough = real2logical(mask_is_fire_intense_enough)

      prob_mort = calc_prob_mortality(nc, df_sf_fc_spread_nr, &
                                      lctlib%q_crown_length, lctlib%q_crown_damage_p, lctlib%q_crown_damage_r_ck, &
                                      lctlib%q_F_scorch_height, & 
                                      lctlib%q_bark_thickness_par1, lctlib%q_bark_thickness_par2, &
                                      height, dbh, & ! Tree characteristics
                                      i_surface, & ! Fire intensity
                                      w_oi(:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx), &
                                      dead_frac_consumption(:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx), & 
                                      mask_opt=is_fire_intense_enough)

!      live_fuel = fuel_cls % Calc_live_fuel()
!      fc_live = sf_calc_live_fuel_consumption(nc, t_sf_fuel_live_fuel_nr, nfuels, &
!                                              pm_prop_ck, &
!                                              config_fract_crown_scorch, &
!                                              sf_area_burned_fract, &
!                                              live_fuel, &
!                                              mask=spitfire % is_fire_intense_enough)
    END IF

    ! derive dead FC (DM * AB)
    ! calc live FC (DM)
    ! calc EM (CN)
    ! backcalculate EM and FC (DM -> CNP, CN -> CNP131415)
    ! calc Residue
    ! Feedbacks

    WRITE(*,*) "update_dfire_post:: area, tile_fract=", area(9), tile_fract(9)

  END SUBROUTINE update_dfire_sf_post

#endif
END MODULE mo_dfire_sf_post
