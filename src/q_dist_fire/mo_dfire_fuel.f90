!>
!>  @par Copyright
!>    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!>    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany
!>    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!>
!>    This file is part of the QUINCY model.
!>    The QUINCY model is free software: you can redistribute it and/or modify
!>    it under the terms of the GNU General Public License version 3,
!>    as published by the Free Software Foundation.
!>
!>    The QUINCY model is distributed in the hope that it will be useful,
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!>    GNU General Public License for more details.
!>
!>    You should have received a copy of the GNU General Public License
!>    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!>
!> @brief module to determine the spitfire fire fuel
!>
!> Main task to compute the spitfire fuel
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_dfire_fuel
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_lnd_bgcm_store,        ONLY: t_lnd_bgcm_store
  USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
  USE mo_lnd_bgcm_store_class,  ONLY: SB_BGCM_POOL_ID, VEG_BGCM_POOL_ID

  USE mo_veg_constants,         ONLY: ITREE, IGRASS
  USE mo_dfire_constants,       ONLY: t_sf_fuel_dead_fuel_woody_ix, t_sf_fuel_dead_fuel_polymeric_ix ! dead fuel indices from sb_pool
  USE mo_dfire_constants,       ONLY: t_sf_fuel_dead_fuel_grass_fruit_ix, t_sf_fuel_dead_fuel_grass_labile_ix, &  ! dead fuel indices from veg_pool
                                      t_sf_fuel_dead_fuel_grass_reserve_ix, t_sf_fuel_dead_fuel_grass_leaf_ix, &
                                      t_sf_fuel_dead_fuel_grass_sap_ix
  USE mo_dfire_constants,       ONLY: df_sf_spq_layer_top_sl

  USE mo_dfire_constants,       ONLY: t_sf_fuel_live_fuel_fruit_ix, t_sf_fuel_live_fuel_leaf_ix, & ! live fuel indices from veg_pool
                                      t_sf_fuel_live_fuel_labile_ix, t_sf_fuel_live_fuel_reserve_ix, &
                                      t_sf_fuel_live_fuel_sap_wood_ix, t_sf_fuel_live_fuel_heart_wood_ix

  USE mo_dfire_constants,       ONLY: sf_fuel_dead_fuel_nr, sf_fuel_live_fuel_nr ! total number of fuel pools
  USE mo_dfire_constants,       ONLY: sf_fuel_classes_nr ! total number of fuel classes
  USE mo_dfire_constants,       ONLY: sf_slow_fuel, sf_fast_fuel ! Fraction for each fuel class [%]
  USE mo_dfire_constants,       ONLY: df_sf_fuel_class_factors !< Fuel bulk density factors for each spreading fuel class (1 to 100 hours)
  USE mo_dfire_constants,       ONLY: df_sf_fc_1hr_idx, df_sf_fc_100hr_idx, df_sf_fc_spread_nr !< Spread fuel classes: 1, 10 and 100 hours

  USE mo_lnd_bgcm_idx,          ONLY: ixC
  USE mo_lnd_bgcm_idx,          ONLY: ix_woody_litter, ix_polymeric_litter
  USE mo_lnd_bgcm_idx,          ONLY: ix_leaf, ix_fruit, ix_labile, ix_reserve, ix_sap_wood, ix_heart_wood
  USE mo_dfire_sf_fuel_class,   ONLY: t_sf_fuel_class

  !USE mo_dfire_constants,       ONLY: live_fuel_frac => live_frac_allocation !< Live fraction ratio for live fuels

  !< Live fraction allocation for dead and live fuels
  USE mo_dfire_constants,       ONLY: df_sf_live_frac_tree, &
                                        df_sf_live_frac_grass
  USE mo_dfire_constants,       ONLY: dead_fuel_frac => df_sf_dead_frac      !< Fuel fraction to allocate from sb_pool

  USE mo_dfire_sf_adaptor,      ONLY: sf_convert_units

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_dfire_fuel

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_fuel'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Compute spitfire for current time step
  !> It links Quincy and Spitfire
  !>
  SUBROUTINE update_dfire_fuel(tile, options)
    USE mo_jsb_class,               ONLY: Get_model
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
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lnd_bgcm_store), POINTER       :: bgcm_store !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc, ic
    INTEGER                               :: no_children
    REAL(wp)                              :: dtime      !< timestep

    INTEGER                               :: lctlib_growthform(options%nc)

    REAL(wp)                              :: tile_fract(options%nc) !<  [%]
    REAL(wp), ALLOCATABLE                 :: top_soil_depth_sl(:) !< Top soil depth layer [m]
    REAL(wp)                              :: config_polymeric_litter_fract !< how much polymeric litter account as dead fuel (0 to 1)

    REAL(wp), POINTER, CONTIGUOUS         :: live_fuel_frac(:,:) !< live fuel classes allocation from veg_pool (#fuel classes,#live pools)
    REAL(wp), ALLOCATABLE                 :: sb_pool_slim(:,:) !< debug only  (#pools, #nc)
    REAL(wp), ALLOCATABLE                 :: dead_fuel_fcl(:,:,:) !< (#nc,#dead fuel pools,#fuel classes)
    REAL(wp), ALLOCATABLE                 :: live_grass_fuel_fcl(:,:,:) !< (#nc,#live fuel pools,#fuel classes)
    REAL(wp), ALLOCATABLE                 :: live_tree_fuel_fcl(:,:,:) !< (#nc,#live fuel pools,#fuel classes)
    REAL(wp), ALLOCATABLE                 :: dead_fuel_spread_fcl(:,:) !< (#nc,#dead fuel pools-1,10,100-,#fuel classes)
    TYPE(t_sf_fuel_class)                 :: fuel_cls

    CHARACTER(len=*), PARAMETER           :: routine = TRIM(modname)//':update_fdire_sf_fuel'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(DFIRE_FUEL_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Local vars here
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real3D_onChunk      :: soil_depth_sl      ! in [m]

    dsl4jsb_Real3D_onChunk      :: dead_fuel             ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: dead_fuel_1hr_dfp     ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: dead_fuel_10hr_dfp    ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: dead_fuel_100hr_dfp   ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: dead_fuel_1000hr_dfp  ! out [g DM m-2]

    dsl4jsb_Real3D_onChunk      :: live_fuel                  ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: live_grass_fuel_1hr_lfp    ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: live_grass_fuel_10hr_lfp   ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: live_grass_fuel_100hr_lfp  ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: live_grass_fuel_1000hr_lfp ! out [g DM m-2]

    dsl4jsb_Real3D_onChunk      :: live_tree_fuel_1hr_lfp     ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: live_tree_fuel_10hr_lfp    ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: live_tree_fuel_100hr_lfp   ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: live_tree_fuel_1000hr_lfp  ! out [g DM m-2]

    dsl4jsb_Real2D_onChunk      :: forest_fract       ! out [%]
    dsl4jsb_Real2D_onChunk      :: grass_fract        ! out [%]

    dsl4jsb_Real2D_onChunk      :: fuel_bulk_density  ! out [kg m-3]
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
    dsl4jsb_Get_config(DFIRE_)
    dsl4jsb_Get_memory(DFIRE_FUEL_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onChunk(SPQ_, soil_depth_sl)                ! in

    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel)             ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_1hr_dfp)     ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_10hr_dfp)    ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_100hr_dfp)   ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_1000hr_dfp)  ! out [g DM m-2]

    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_fuel)             ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_grass_fuel_1hr_lfp)    ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_grass_fuel_10hr_lfp)   ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_grass_fuel_100hr_lfp)  ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_grass_fuel_1000hr_lfp) ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_tree_fuel_1hr_lfp)     ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_tree_fuel_10hr_lfp)    ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_tree_fuel_100hr_lfp)   ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_tree_fuel_1000hr_lfp)  ! out [g DM m-2]

    dsl4jsb_Get_var2D_onChunk(DFIRE_FUEL_, forest_fract)          ! out [%]
    dsl4jsb_Get_var2D_onChunk(DFIRE_FUEL_, grass_fract)           ! out [%]

    dsl4jsb_Get_var2D_onChunk(DFIRE_FUEL_, fuel_bulk_density)     ! out [kg m-3]
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    CALL tile%Get_fraction(ics, ice, iblk, tile_fract)
    !WRITE(*,*) "update_dfire_fuel:: computing fuel... (growth form=", lctlib%growthform ,")", ", fraction=", &
    !             fract(1)
    dead_fuel = 0._wp
    live_fuel = 0._wp
    forest_fract=0._wp
    grass_fract=0._wp
    fuel_bulk_density = 0._wp
    top_soil_depth_sl = soil_depth_sl(:,df_sf_spq_layer_top_sl)
    lctlib_growthform(:) = lctlib%growthform
    config_polymeric_litter_fract = dsl4jsb_Config(DFIRE_)%ratio_polymeric_litter_to_dead_fuel

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !

    ! Aggregate fractions by type (Forest, Grass)
    ! todo: move to init (static setup)
    IF (lctlib%ForestFlag) THEN
      forest_fract=tile_fract
    ELSE IF (lctlib%GrassFlag) THEN
      grass_fract=tile_fract
    !ELSE ! ignore
    END IF

    !< PFT related, todo: no need to redefine at every time step
    IF (lctlib%growthform == itree) THEN
      live_fuel_frac => df_sf_live_frac_tree
    ELSE !< %growthform == igrass
      live_fuel_frac => df_sf_live_frac_grass
    END IF

    sb_pool_slim = sb_pool_mt(:,ixC,:,df_sf_spq_layer_top_sl)
    ALLOCATE(dead_fuel, mold=sb_pool_slim) ! unsure why is needed, otherwise NAG compiler raises an error. function should be allocated dead_fuel
    dead_fuel = sf_convert_units("molC/m3", "gDM/m2", sb_pool_slim, aux = top_soil_depth_sl, nc_dim_opt=2)

    ALLOCATE(live_fuel, mold=veg_pool_mt(:,ixC,:)) ! unsure why is needed, NAG compiler raises an error. sf_convert_units should already allocate live_fuel
    live_fuel = sf_convert_units("molC/m2", "gDM/m2", veg_pool_mt(:,ixC,:))

    fuel_cls = t_sf_fuel_class(nc, lctlib_growthform, &
                            dead_fuel_frac, live_fuel_frac, &
                            dead_fuel, live_fuel, &
                            config_polymeric_litter_fract)

    live_grass_fuel_fcl = fuel_cls % calc_live_grass_fuel()
    live_tree_fuel_fcl = fuel_cls % calc_live_tree_fuel()
    dead_fuel_fcl = fuel_cls % calc_dead_fuel()

    !> Calculate the weighted average of the Fuel bulk density for the current PFT.
    !> Using the Pfeiffer et al. (2013) paper, equation A17. Remaining equations should be implemented
    !> when accounting for the other pfts.
    !> Thonicke's et al. (2010) is not used as the equation does not seems correct when accounting for 10 and 100 hours
    !> fuel bulk densities (slightly more).
    !> Aggregation average's the PFT values

    dead_fuel_spread_fcl = SUM(dead_fuel_fcl(:,:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx), DIM=2) !< 2 = fuel pools dimension index
    DO ic = 1, nc
      fuel_bulk_density(:) = lctlib%fbd_pb * SUM(dead_fuel_spread_fcl(ic,:) * df_sf_fuel_class_factors)
    END DO

    dead_fuel_1hr_dfp = dead_fuel_fcl(:,:,df_sf_fc_1hr_idx)
    dead_fuel_10hr_dfp = dead_fuel_fcl(:,:,df_sf_fc_10hr_idx)
    dead_fuel_100hr_dfp = dead_fuel_fcl(:,:,df_sf_fc_100hr_idx)
    dead_fuel_1000hr_dfp = dead_fuel_fcl(:,:,df_sf_fc_1000hr_idx)
    
    live_grass_fuel_1hr_lfp = live_grass_fuel_fcl(:,:,df_sf_fc_1hr_idx)
    live_grass_fuel_10hr_lfp = live_grass_fuel_fcl(:,:,df_sf_fc_10hr_idx)
    live_grass_fuel_100hr_lfp = live_grass_fuel_fcl(:,:,df_sf_fc_100hr_idx)
    live_grass_fuel_1000hr_lfp = live_grass_fuel_fcl(:,:,df_sf_fc_1000hr_idx)

    live_tree_fuel_1hr_lfp = live_tree_fuel_fcl(:,:,df_sf_fc_1hr_idx)
    live_tree_fuel_10hr_lfp = live_tree_fuel_fcl(:,:,df_sf_fc_10hr_idx)
    live_tree_fuel_100hr_lfp = live_tree_fuel_fcl(:,:,df_sf_fc_100hr_idx)
    live_tree_fuel_1000hr_lfp = live_tree_fuel_fcl(:,:,df_sf_fc_1000hr_idx)

    WRITE(*,*) "update_dfire_fuel:: computing fuel... (growth form=", lctlib%growthform ,")", &
                        ", fraction, forest_fract=", &
                        tile_fract(1), forest_fract(1)

    NULLIFY(live_fuel_frac)
  END SUBROUTINE update_dfire_fuel


#endif
END MODULE mo_dfire_fuel

