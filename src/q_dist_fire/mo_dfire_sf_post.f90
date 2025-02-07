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

  USE mo_dfire_constants,       ONLY: df_sf_spq_layer_top_sl
  USE mo_dfire_constants,       ONLY: df_sf_fc_spread_nr, df_sf_fc_full_nr
  USE mo_dfire_constants,       ONLY: df_sf_fc_1hr_idx, df_sf_fc_10hr_idx, df_sf_fc_100hr_idx, df_sf_fc_1000hr_idx
  USE mo_dfire_constants,       ONLY: sf_fuel_dead_fuel_nr, sf_fuel_live_fuel_nr
  USE mo_dfire_constants,       ONLY: df_sf_species_nr
  USE mo_dfire_constants,       ONLY: df_sf_elements_nr

  USE mo_dfire_constants,       ONLY: lctlib_emission_factors !< todo: check values depends on PFT

  USE mo_dfire_sf_damage,       ONLY: t_sf_prob_mort, calc_prob_mortality
  USE mo_dfire_util,            ONLY: real2logical

  USE mo_dfire_sf_fuel_consump_class, ONLY: sf_calc_live_fuel_consumption, &
                                            sf_calc_dead_fuel_consumption_populate_elements, &
                                            sf_calc_live_fuel_consumption_populate_elements

  USE mo_dfire_sf_emissions_class,    ONLY: sf_calc_emissions, sf_calc_emissions_elements
  USE mo_lnd_bgcm_store_class,  ONLY: SB_BGCM_POOL_ID, VEG_BGCM_POOL_ID

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
    dsl4jsb_Use_config(DFIRE_)
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

    TYPE(t_sf_prob_mort)                  :: prob_mort     !< Tree mortality
    REAL(wp)                              :: live_fuel_cls(options%nc, sf_fuel_live_fuel_nr, df_sf_fc_full_nr) !< live fuel for 1, 10, 100 and 1000 hours fuel classes [g DM m-2]

    REAL(wp), ALLOCATABLE                 :: fc_live_mol_m2(:,:)  !< Total live fuel consumption (nc, live fuels) [molX/m2]
    REAL(wp), ALLOCATABLE                 :: fc_dead_mol_m2(:,:)  !< Total dead fuel consumption (nc, t_sf_fuel_dead_fuel_nr) [molX/m2]

    REAL(wp), ALLOCATABLE                 :: fc_dead_el(:,:,:) !< dead fuel consumption (#nc, #dead pools, #Elements) [molX m-2]
    REAL(wp), ALLOCATABLE                 :: fc_live_el(:,:,:)   !< live fuel consumption (#nc, #live pools, #Elemens) [molX m-2]

    REAL(wp)                              :: residue_dead_el(options%nc,sf_fuel_dead_fuel_nr,df_sf_elements_nr) !< emissions (#nc, #dead pools, #elements) [mol m-2]
    REAL(wp)                              :: residue_live_el(options%nc,sf_fuel_live_fuel_nr,df_sf_elements_nr) !< emissions (#nc, #live pools, #elements) [mol m-2]
                                          
    REAL(wp)                              :: mass_imbalance_emitted_dead_el(options%nc,sf_fuel_dead_fuel_nr,df_sf_elements_nr) !< emissions (#nc, #dead pools, #elems) [mol m-2]
    REAL(wp)                              :: mass_imbalance_emitted_live_el(options%nc,sf_fuel_live_fuel_nr,df_sf_elements_nr) !< emissions (#nc, #live pools, #elems) [mol m-2]

    REAL(wp), ALLOCATABLE                 :: sb_pool_in_m2(:,:,:) !< BGCM sb_pool (#pools,#elements,#nc) [molX/m2]

    INTEGER                               :: ip ! iterate

    CHARACTER(len=*), PARAMETER           :: routine = TRIM(modname)//':update_fdire_sf_fuel'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real3D_onChunk                :: live_fuel_1hr_lfp         ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_fuel_10hr_lfp        ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_fuel_100hr_lfp       ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: live_fuel_1000hr_lfp      ! in [g DM m-2]

    dsl4jsb_Real2D_onChunk                :: area_burned               ! in [ha ts-1]
    dsl4jsb_Real2D_onChunk                :: area_burned_fract         ! in [-]
    dsl4jsb_Real3D_onChunk                :: dead_frac_consumption     ! in [-]
    dsl4jsb_Real2D_onChunk                :: i_surface                 ! in [kW m-1]
    dsl4jsb_Real2D_onChunk                :: mask_is_fire_intense_enough ! in [-]
    dsl4jsb_Real2D_onChunk                :: pm_prop_ck                ! out [-]

    dsl4jsb_Real3D_onChunk                :: dead_fuel_consumption_1hr ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_consumption_10hr ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_consumption_100hr ! in [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_consumption_1000hr ! in [g DM m-2]

    dsl4jsb_Real2D_onChunk                :: diameter                   ! in [?]
    dsl4jsb_Real2D_onChunk                :: height                     ! in [?]
    dsl4jsb_Real3D_onChunk                :: soil_depth_sl              ! in [m]
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(DFIRE_)
    dsl4jsb_Def_memory(DFIRE_FUEL_)
    dsl4jsb_Def_memory(DFIRE_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(SPQ_)
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
    dsl4jsb_Get_config(DFIRE_)
    dsl4jsb_Get_memory(DFIRE_FUEL_)
    dsl4jsb_Get_memory(DFIRE_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onChunk(SPQ_, soil_depth_sl)               ! in [?]

    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_fuel_1hr_lfp)    ! in [?]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_fuel_10hr_lfp)   ! in [?]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_fuel_100hr_lfp)  ! in [?]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, live_fuel_1000hr_lfp) ! in [?]

    dsl4jsb_Get_var2D_onChunk(DFIRE_, area_burned)               ! in [ha ts-1]
    dsl4jsb_Get_var2D_onChunk(DFIRE_, area_burned_fract)         ! in [-]
    dsl4jsb_Get_var3D_onChunk(DFIRE_, dead_frac_consumption)     ! in [-]
    dsl4jsb_Get_var2D_onChunk(DFIRE_, i_surface)                 ! in [kW m-1]
    dsl4jsb_Get_var2D_onChunk(DFIRE_, mask_is_fire_intense_enough)! in [-]
    dsl4jsb_Get_var2D_onChunk(DFIRE_, pm_prop_ck)                ! out [-]

    dsl4jsb_Get_var2D_onChunk(VEG_, diameter)                    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, height)                      ! in

    dsl4jsb_Get_var3D_onChunk(DFIRE_, dead_fuel_consumption_1hr) ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_, dead_fuel_consumption_10hr) ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_, dead_fuel_consumption_100hr) ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_, dead_fuel_consumption_1000hr) ! in [g DM m-2]

    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    ! ----------------------------------------------------------------------------------------------------- !

    ! init here
    area => grid %area (ics:ice, iblk)
    tile_fract => tile%fract(ics:ice, iblk)
    is_fire_intense_enough = .FALSE.

    prob_mort = t_sf_prob_mort(nc)
    !< Rebuild dead fuel consumption
    !< * Derive dead FC (DM * AB) *
    !< Calculate from Tile to PFT (* tile_fract)
    DO ip=1, sf_fuel_dead_fuel_nr ! iterate over dead fuel pools
      fc_dead_cls(:,ip,df_sf_fc_1hr_idx) = dead_fuel_consumption_1hr(:,ip) * tile_fract
      fc_dead_cls(:,ip,df_sf_fc_10hr_idx) = dead_fuel_consumption_10hr(:,ip) * tile_fract
      fc_dead_cls(:,ip,df_sf_fc_100hr_idx) = dead_fuel_consumption_100hr(:,ip) * tile_fract
      fc_dead_cls(:,ip,df_sf_fc_1000hr_idx) = dead_fuel_consumption_1000hr(:,ip) * tile_fract
    ENDDO
    live_fuel_cls(:,:,df_sf_fc_1hr_idx) = live_fuel_1hr_lfp
    live_fuel_cls(:,:,df_sf_fc_10hr_idx) = live_fuel_10hr_lfp
    live_fuel_cls(:,:,df_sf_fc_100hr_idx) = live_fuel_100hr_lfp
    live_fuel_cls(:,:,df_sf_fc_1000hr_idx) = live_fuel_1000hr_lfp
    area_burned_pft = area_burned * tile_fract
    area_burned_fract_pft = area_burned_fract * tile_fract

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !


    ! calc tree mortality
    IF (lctlib%growthform == itree) THEN
      dbh(:) = sf_convert_units("m", "cm", diameter) !< dbh to spitfire units
      !WRITE(*,*) "update_dfire_post:: area, tile_fract, height, dbh=", &
      !      area(1), tile_fract(1), (tile_fract(1) > 0._wp), height(1), dbh(1)

      is_fire_intense_enough = real2logical(mask_is_fire_intense_enough) .AND. &
                                (tile_fract > 0._wp) .AND. & ! tile is empty, skip 
                                (height > 0._wp) ! height is not yet initialized, skip

      prob_mort = calc_prob_mortality(nc, df_sf_fc_spread_nr, &
                                      lctlib%q_crown_length, lctlib%q_crown_damage_p, lctlib%q_crown_damage_r_ck, &
                                      lctlib%q_F_scorch_height, & 
                                      lctlib%q_bark_thickness_par1, lctlib%q_bark_thickness_par2, &
                                      height, dbh, & ! Tree characteristics
                                      i_surface, & ! Fire intensity
                                      w_oi(:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx), &
                                      dead_frac_consumption(:,df_sf_fc_1hr_idx:df_sf_fc_100hr_idx), & 
                                      mask_opt=is_fire_intense_enough)

      fc_live = sf_calc_live_fuel_consumption(nc, sf_fuel_live_fuel_nr, df_sf_fc_full_nr, &
                                              prob_mort % prop_ck, &
                                              dsl4jsb_Config(DFIRE_)%fract_crown_scorch, &
                                              area_burned_fract, &
                                              live_fuel_cls, &
                                              mask=is_fire_intense_enough)

      pm_prop_ck = prob_mort % prop_ck
    END IF

    ! * Calc Emissions (by Species) *
    !< Calculate fire damage to plants (total biomass burned and tree mortality) based on total dead and
    !< live fuel, as in Peiffer et al. (2013)
    emissions_dead_spc = sf_calc_emissions(nc, sf_fuel_dead_fuel_nr, df_sf_species_nr, lctlib_emission_factors, &
                                                        SUM(fc_dead_cls, dim=3), & !< 3 = fuel classes dim index
                                                        mask=is_fire_intense_enough)
    emissions_live_spc = sf_calc_emissions(nc, sf_fuel_live_fuel_nr, df_sf_species_nr, lctlib_emission_factors, &
                                                        fc_live, &
                                                        mask=is_fire_intense_enough)

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

    !< Derive missing nutrients/isotopes from emissions
    CALL sf_calc_emissions_elements(nc,sf_fuel_dead_fuel_nr, 1, & ! in, 1=dead fuel type, todo: weird approach
                                    fc_dead_el, & ! [molX/m2]
                                    emissions_dead_spc, & ! [molX/m2]
                                    emissions_dead_el, & ! out [molX/m2]
                                    mask=is_fire_intense_enough)

    CALL sf_calc_emissions_elements(nc,sf_fuel_live_fuel_nr, 2, & ! in, 2 = live fuel type
                                    fc_live_el, & ! [molX/m2]
                                    emissions_live_spc, & ! [molX/m2]
                                    emissions_live_el, & ! out [molX/m2]
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

!    WRITE(*,*) "update_dfire_post:: area, tile_fract=", area(9), tile_fract(9)

  END SUBROUTINE update_dfire_sf_post


  ! ======================================================================================================= !
  !>
  !> Calculate residue and any mass imbalance
  !>
  SUBROUTINE sf_calc_residue(fuel_consumption, emitted_elements, residue, mass_overflow)
    REAL(wp), INTENT(in) :: fuel_consumption(:,:,:) !< (#nc,#pools,#elements) [molX/m2]
    REAL(wp), INTENT(in) :: emitted_elements(:,:,:) !< (#nc,#pools,#elements) [molX/m2]
    REAL(wp), INTENT(out) :: residue(:,:,:)         !< (#nc,#pools,#elements) [molX/m2]
    REAL(wp), INTENT(out) :: mass_overflow(:,:,:)   !< (#nc,#pools,#elements) [molX/m2]

    !< Compute residue, where Emitted elements < FC
    residue = fuel_consumption - emitted_elements
    mass_overflow = 0._wp

    WHERE (residue < 0._wp)
      mass_overflow = ABS(residue)
      residue = fuel_consumption ! set to max value
    END WHERE
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
END MODULE mo_dfire_sf_post
