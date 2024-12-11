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

  USE mo_dfire_constants,       ONLY: t_sf_fuel_live_fuel_fruit_ix, t_sf_fuel_live_fuel_leaf_ix, & ! live fuel indices from veg_pool
                                      t_sf_fuel_live_fuel_labile_ix, t_sf_fuel_live_fuel_reserve_ix, &
                                      t_sf_fuel_live_fuel_sap_wood_ix, t_sf_fuel_live_fuel_heart_wood_ix

  USE mo_dfire_constants,       ONLY: sf_fuel_dead_fuel_nr, sf_fuel_live_fuel_nr ! total number of fuel pools
  USE mo_dfire_constants,       ONLY: sf_fuel_classes_nr ! total number of fuel classes

  USE mo_lnd_bgcm_idx,          ONLY: ixC
  USE mo_lnd_bgcm_idx,          ONLY: ix_woody_litter, ix_polymeric_litter
  USE mo_lnd_bgcm_idx,          ONLY: ix_leaf, ix_fruit, ix_labile, ix_reserve, ix_sap_wood, ix_heart_wood

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
    USE mo_jsb_process_class,       ONLY: VEG_, SPQ_, DFIRE_FUEL_
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_math_constants,      ONLY: one_day

    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(DFIRE_FUEL_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(VEG_)
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

    REAL(wp) :: dead_fuel(options%nc,sf_fuel_dead_fuel_nr,sf_fuel_classes_n) !< (#nc,#dead pools,#fuel classes) [g DM m-2]
    REAL(wp) :: live_fuel(options%nc,sf_fuel_live_fuel_nr,sf_fuel_classes_n) !< (#nc,#live pools,#fuel classes) [g DM m-2]

    CHARACTER(len=*), PARAMETER           :: routine = TRIM(modname)//':update_fdire_sf_fuel'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(DFIRE_FUEL_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Local vars here
    ! ----------------------------------------------------------------------------------------------------- !

    dsl4jsb_Real3D_onChunk      :: dead_fuel_1hr    ! out [g DM m-2]
    dsl4jsb_Real3D_onChunk      :: dead_fuel_10hr
    dsl4jsb_Real3D_onChunk      :: dead_fuel_100hr
    dsl4jsb_Real3D_onChunk      :: dead_fuel_1000hr
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
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(DFIRE_FUEL_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_1hr)         ! out [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_10hr)
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_100hr)
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_1000hr)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    WRITE(*,*) "update_dfire_fuel:: computing fuel... (growth form=", lctlib%growthform ,")"
    dead_fuel = 0._wp
    live_fuel = 0._wp

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !

    IF (lctlib%growthform == ITREE) THEN
      dead_fuel(:,t_sf_fuel_dead_fuel_woody_ix,:) = sb_pool_mt(ix_woody_litter,ixC,:,1) * sf_slow_fuel
      dead_fuel(:,t_sf_fuel_dead_fuel_polymeric_ix,:) = sb_pool_mt(ix_polymeric_litter,ixC,:,1) * sf_slow_fuel
      dead_fuel(:,t_sf_fuel_dead_fuel_grass_leaf_ix,:) = veg_pool_mt(ix_leaf,ixC,:) * sf_fast_fuel
      dead_fuel(:,t_sf_fuel_dead_fuel_grass_fruit_ix,:) = veg_pool_mt(ix_fruit,ixC,:) * sf_fast_fuel
      dead_fuel(:,t_sf_fuel_dead_fuel_grass_labile_ix,:) = veg_pool_mt(ix_labile,ixC,:) * sf_fast_fuel
      dead_fuel(:,t_sf_fuel_dead_fuel_grass_reserve_ix,:) = veg_pool_mt(ix_reserve,ixC,:) * sf_fast_fuel
      dead_fuel(:,t_sf_fuel_dead_fuel_grass_sap_ix,:) = veg_pool_mt(ix_sap_wood,ixC,:) * sf_fast_fuel
    ELSE IF (lctlib%growthform == IGRASS) THEN
      live_fuel(:,t_sf_fuel_live_fuel_leaf_ix,:) = veg_pool_mt(ix_leaf,ixC,:) * sf_fast_fuel
      live_fuel(:,t_sf_fuel_live_fuel_fruit_ix,:) = veg_pool_mt(ix_fruit,ixC,:) * sf_fast_fuel
      live_fuel(:,t_sf_fuel_live_fuel_labile_ix,:) = veg_pool_mt(ix_labile,ixC,:) * sf_slow_fuel
      live_fuel(:,t_sf_fuel_live_fuel_reserve_ix,:) = veg_pool_mt(ix_reserve,ixC,:) * sf_slow_fuel
      live_fuel(:,t_sf_fuel_live_fuel_sap_wood_ix,:) = veg_pool_mt(ix_sap_wood,ixC,:) * sf_slow_fuel
      live_fuel(:,t_sf_fuel_live_fuel_heart_wood_ix,:) = veg_pool_mt(ix_heart_wood,ixC,:) * sf_slow_fuel
    ELSE ! something unexpected
      CALL finish(routine, "Unexpected growth form found")
    END IF

    dead_fuel_1hr = dead_fuel(:,:,1)
    dead_fuel_10hr = dead_fuel(:,:,2)
    dead_fuel_100hr = dead_fuel(:,:,3)
    dead_fuel_1000hr = dead_fuel(:,:,4)

    live_fuel_1hr = live_fuel(:,:,1)
    live_fuel_10hr = live_fuel(:,:,2)
    live_fuel_100hr = live_fuel(:,:,3)
    live_fuel_1000hr = live_fuel(:,:,4)

  END SUBROUTINE update_dfire_fuel


#endif
END MODULE mo_dfire_fuel
