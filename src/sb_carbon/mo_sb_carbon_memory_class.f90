!> Contains the memory class for the carbon process.
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
MODULE mo_sb_carbon_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_cqt_class,          ONLY: LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, &
    &                                  BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real1d, t_jsb_var_real2d， t_jsb_var_real3d

  ! Use of prcesses in this module
  dsl4jsb_Use_processes PHENO_, SB_CARBON_, NLCC_

  ! Use process configurations
  dsl4jsb_Use_config(PHENO_)    ! needed for l_forestRegrowth
  dsl4jsb_Use_config(SB_CARBON_)    ! needed for l_forestRegrowth

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_sb_carbon_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 155

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_sb_carbon_memory

    TYPE(t_jsb_var_real2d) ::     &
      c_leaf,                 & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      c_fineroot,             & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      c_fruit,             & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      c_seedbed,             & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      c_coarseroot,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                    ! .. plants [mol(C)/m^2(canopy)]
      c_sapwood,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                    ! .. plants [mol(C)/m^2(canopy)]
      c_heartwood,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                    ! .. plants [mol(C)/m^2(canopy)]
      tree_height,            & ! for quincy biomass allocation calculation
      f_lai                  & ! growth stage indicator for height estimation


    TYPE(t_jsb_var_real2d) ::     &
      litterfall_leaf,                 & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      litterfall_fineroot,             & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      litterfall_fruit,             & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      litterfall_seedbed,             & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      litterfall_sapwood,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                    ! .. plants [mol(C)/m^2(canopy)]
      litterfall_heartwood,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                    ! .. plants [mol(C)/m^2(canopy)]
      litterfall_coarseroot               & ! C-pool for stems, thick roots and other (dead) structural material of living
                                    ! .. plants [mol(C)/m^2(canopy)]

    TYPE(t_jsb_var_real2d) ::     &
      sb_pool_ag_litter,                 & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      sb_pool_bg_litter,             & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      sb_pool_bg_humus,             & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      sb_pool_bg_mic,             & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      sb_pool_bg_res            & ! C-pool for stems, thick roots and other (dead) structural material of living
                                    ! .. plants [mol(C)/m^2(canopy)]
    

   ! Output on tile areas:
   TYPE(t_jsb_var_real2d) ::         &
      c_leaf_ta,                & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                       ! .. of vegetation [mol(C)/m^2(tile area)]
      c_fineroot_ta,                & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                       ! .. of vegetation [mol(C)/m^2(tile area)]
      c_fruit_ta,                & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                       ! .. of vegetation [mol(C)/m^2(tile area)]
      c_seedbed_ta,                & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                       ! .. of vegetation [mol(C)/m^2(tile area)]
      c_coarseroot_ta,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                       ! .. plants [mol(C)/m^2(tile area)]
      c_sapwood_ta,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                       ! .. plants [mol(C)/m^2(tile area)]
      c_heartwood_ta,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                       ! .. plants [mol(C)/m^2(tile area)]
      tree_height_ta,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                       ! .. plants [mol(C)/m^2(tile area)]
      f_lai_ta,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                       ! .. plants [mol(C)/m^2(tile area)]
      litterfall_leaf_ta,                 & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(tile area)]
      litterfall_fineroot_ta,             & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(tile area)]
      litterfall_fruit_ta,             & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(tile area)]
      litterfall_seedbed_ta,             & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(tile area)]
      litterfall_sapwood_ta,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                    ! .. plants [mol(C)/m^2(tile area)]
      litterfall_heartwood_ta,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                    ! .. plants [mol(C)/m^2(tile area)]
      litterfall_coarseroot_ta               & ! C-pool for stems, thick roots and other (dead) structural material of living
                                    ! .. plants [mol(C)/m^2(tile area)]
      sb_pool_ag_litter_ta,
      sb_pool_bg_litter_ta,
      sb_pool_bg_humus_ta
      

    TYPE(t_jsb_var_real3d)          :: &    ! DIMENSION(soil_depth_water)
                                    root_fraction_sl_fixed, &               !< fixed fractions of root per layer from QUINCY veg [unitless]
                                    fineroot_carbon_sl_fixed, &
                                    sb_pool_bg_humus_sl, &
                                    sb_pool_bg_mic_sl, &
                                    sb_pool_bg_res_sl
                                    

  CONTAINS
    PROCEDURE :: Init => Init_sb_carbon_memory
  END TYPE t_sb_carbon_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_carbon_memory_class'

CONTAINS

  SUBROUTINE Init_sb_carbon_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_sb_carbon_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: lib_id     !< id of primary lct in lctlib
    INTEGER,             INTENT(in)    :: model_id

    dsl4jsb_Def_config(PHENO_)    ! needed for l_forestRegrowth
    dsl4jsb_Def_config(SB_CARBON_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid     ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface   ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_carbon_memory'

    IF (lct_ids(1) > 0) CONTINUE ! avoid compiler warning about dummy argument not being used

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_config(SB_CARBON_)

    CALL mem%Add_var('pseudo_temp', mem%pseudo_temp,                                                 &
      & hgrid, surface,                                                                              &
      & t_cf('pseudo_temp', 'K', 'Pseudo-mean air temperature over previous 15 days'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & output_level=MEDIUM,                                                                         &
      & lrestart=.TRUE., initval_r=285.0_wp ) ! initval taken from JSBACH3, arbitary.

    ! R: In the output it should be the same as pseudo_temp, I do not understand why it was written out in JS3.
    CALL mem%Add_var('pseudo_temp_yDay', mem%pseudo_temp_yDay,                                         &
      & hgrid, surface,                                                                                &
      & t_cf('pseudo_temp_yDay', 'K', 'Yesterdays pseudo-mean air temperature over previous 15 days'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & lrestart=.TRUE., initval_r=285.0_wp ) ! initval taken from JSBACH3, arbitary.

    CALL mem%Add_var('pseudo_precip', mem%pseudo_precip,                                             &
      & hgrid, surface,                                                                              &
      & t_cf('pseudo_precip', 'Kg m-2 s-1', 'Pseudo-mean precipitation rate over previous 15 days'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & output_level=MEDIUM,                                                                         &
      & lrestart=.TRUE., initval_r=3.0E-5_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3, arbitary.

    ! R: In the output it should be the same as pseudo_precip, I do not understand why it was written out in JS3.
    CALL mem%Add_var('pseudo_precip_yDay', mem%pseudo_precip_yDay,                                   &
      & hgrid, surface,                                                                              &
      & t_cf('pseudo_precip_yDay', 'Kg m-2 s-1', 'Yesterdays pseudo-mean precipitation rate over previous 15 days'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & lrestart=.TRUE., initval_r=3.0E-5_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3, arbitary.

    ! Carbon variables on canopy area
    CALL mem%Add_var('c_greenwood2litter', mem%cflux_c_greenwood_2_litter,                             &
      & hgrid, surface,                                                                                &
      & t_cf('c_greenwood2litter', 'mol(CO2) m-2(canopy) s-1', 'Total litter flux entering the soil pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=FULL,                                                                             &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_green', mem%c_green,                                                   &
      & hgrid, surface,                                                                        &
      & t_cf('c_green', 'mol(C) m-2(canopy)', 'C-Pool for green parts of vegetation.'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = LIVE_CARBON_CQ_TYPE,    &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_woods', mem%c_woods,                                                              &
      & hgrid, surface,                                                                                   &
      & t_cf('c_woods', 'mol(C) m-2(canopy)', 'C-Pool for structural material/parts of vegetation.'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                &
      & prefix, suffix,                                                                                   &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = LIVE_CARBON_CQ_TYPE,               &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_reserve', mem%c_reserve,                                                                           &
      & hgrid, surface,                                                                                                    &
      & t_cf('c_reserve', 'mol(C) m-2(canopy)', 'C-Pool for reserve carbohydrates (starches, sugars) of vegetation.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = LIVE_CARBON_CQ_TYPE,            &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_crop_harvest', mem%c_crop_harvest,                                             &
      & hgrid, surface,                                                                                &
      & t_cf('c_crop_harvest', 'mol(C) m-2(canopy)', 'C-Pool for Harvested Material from Crops.'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = PRODUCT_CARBON_CQ_TYPE,         &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Yasso 1: Green carbon pools
    CALL mem%Add_var('c_acid_ag1', mem%c_acid_ag1,                                                               &
      & hgrid, surface,                                                                                          &
      & t_cf('c_acid_ag1', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_ag1', mem%c_water_ag1,                                                               &
      & hgrid, surface,                                                                                            &
      & t_cf('c_water_ag1', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                          &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_ag1', mem%c_ethanol_ag1,                                                               &
      & hgrid, surface,                                                                                                &
      & t_cf('c_ethanol_ag1', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix,                                                                                                &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                             &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_ag1', mem%c_nonsoluble_ag1,                                                        &
      & hgrid, surface,                                                                                               &
      & t_cf('c_nonsoluble_ag1', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                             &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_acid_bg1', mem%c_acid_bg1,                                                               &
      & hgrid, surface,                                                                                          &
      & t_cf('c_acid_bg1', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                       &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_bg1', mem%c_water_bg1,                                                               &
      & hgrid, surface,                                                                                            &
      & t_cf('c_water_bg1', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                         &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_bg1', mem%c_ethanol_bg1,                                                              &
      & hgrid, surface,                                                                                               &
      & t_cf('c_ethanol_bg1', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                            &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_bg1', mem%c_nonsoluble_bg1,                                                     &
      & hgrid, surface,                                                                                            &
      & t_cf('c_nonsoluble_bg1', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (belowground).'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                          &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_humus_1', mem%c_humus_1,                                     &
      & hgrid, surface,                                                              &
      & t_cf('c_humus_1', 'mol(C) m-2(canopy)', 'C-Pool for humus in Yasso.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=FULL,  l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,  &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Yasso 2: Woody carbon pools
    CALL mem%Add_var('c_acid_ag2', mem%c_acid_ag2,                                                               &
      & hgrid, surface,                                                                                          &
      & t_cf('c_acid_ag2', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_ag2', mem%c_water_ag2,                                                               &
      & hgrid, surface,                                                                                            &
      & t_cf('c_water_ag2', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                          &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_ag2', mem%c_ethanol_ag2,                                                               &
      & hgrid, surface,                                                                                                &
      & t_cf('c_ethanol_ag2', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix,                                                                                                &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                              &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_ag2', mem%c_nonsoluble_ag2,                                                        &
      & hgrid, surface,                                                                                               &
      & t_cf('c_nonsoluble_ag2', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = AG_DEAD_C_CQ_TYPE,                             &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_acid_bg2', mem%c_acid_bg2,                                                               &
      & hgrid, surface,                                                                                          &
      & t_cf('c_acid_bg2', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_bg2', mem%c_water_bg2,                                                               &
      & hgrid, surface,                                                                                            &
      & t_cf('c_water_bg2', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
      & prefix, suffix,                                                                                            &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                          &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_bg2', mem%c_ethanol_bg2,                                                              &
      & hgrid, surface,                                                                                               &
      & t_cf('c_ethanol_bg2', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                             &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_bg2', mem%c_nonsoluble_bg2,                                                       &
      & hgrid, surface,                                                                                              &
      & t_cf('c_nonsoluble_bg2', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (belowground).'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
      & prefix, suffix,                                                                                              &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,                            &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_humus_2', mem%c_humus_2,                                     &
      & hgrid, surface,                                                              &
      & t_cf('c_humus_2', 'mol(C) m-2(canopy)', 'C-Pool for humus in Yasso.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=FULL, l_conserve_quan=.TRUE., cons_quan_type_id = BG_DEAD_C_CQ_TYPE,  &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Carbon variables needed to calculate natural land cover changes (NLCC)
    IF (model%Is_process_enabled(NLCC_)) THEN
      CALL mem%Add_var('max_green_bio', mem%max_green_bio,                                     &
        & hgrid, surface,                                                                      &
        & t_cf('max_green_bio', 'mol(C) m-2(canopy)', 'annual maximum content of green C-Pool.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
        & prefix, suffix,                                                                      &
        & output_level=FULL,                                                                   &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('sla', mem%sla,                                                         &
        & hgrid, surface,                                                                      &
        & t_cf('sla', 'm2(leave)/mol(C)', 'specific leave area.'),                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
        & prefix, suffix,                                                                      &
        & output_level=FULL,                                                                   &
       & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )
    END IF

    CALL mem%Add_var('lai_sum', mem%LAI_sum,                                         &
      & hgrid, surface,                                                              &
      & t_cf('LAI_sum', '-', 'Used to accumulate LAI over a day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_pot_sum', mem%NPP_pot_sum,                                         &
      & hgrid, surface,                                                                      &
      & t_cf('NPP_pot_sum', 'mol(C) m-2(canopy) s-1', 'Used to accumulate NPP over a day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
      & prefix, suffix,                                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE.) ! initval taken from JSBACH3

    CALL mem%Add_var('gpp_sum', mem%GPP_sum,                                             &
      & hgrid, surface,                                                                  &
      & t_cf('GPP_sum', 'mol(C) m-2(canopy) s-1', 'Used to accumulate GPP over a day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('soil_respiration', mem%soil_respiration,                 &
      & hgrid, surface,                                                        &
      & t_cf('soil_respiration', 'mol(C) m-2(canopy) s-1', 'Soil respiration.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_flux_correction', mem%NPP_flux_correction,                                 &
      & hgrid, surface,                                                                              &
      & t_cf('NPP_flux_correction', 'mol(C) m-2(canopy) s-1',                                        &
      &      'Amount by which the NPP rate has to be corrected to avoid a negative reserve pool.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('excess_npp', mem%excess_NPP,                                                                     &
      & hgrid, surface,                                                                                                &
      & t_cf('excess_NPP', 'mol(C) m-2 (canopy) s-1', 'NPP that could not be not stored in vegetation carbon pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix,                                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('lai_yydaymean', mem%LAI_yyDayMean,                    &
      & hgrid, surface,                                                     &
      & t_cf('LAI_yyDayMean', '-', 'Previous previous day mean.'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var('lai_ydaymean', mem%LAI_yDayMean,                      &
      & hgrid, surface,                                                     &
      & t_cf('LAI_yDayMean', '-', 'Previous day mean.'),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_pot_ydaymean', mem%NPP_pot_yDayMean,                                    &
      & hgrid, surface,                                                                           &
      & t_cf('NPP_pot_yDayMean', 'mol(C) m-2(canopy) s-1', 'Mean NPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
      & prefix, suffix,                                                                           &
      & lrestart=.TRUE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would
                                                  ! result in a zero divide in the first time step of an initialized run.

    CALL mem%Add_var('npp_act_ydaymean', mem%NPP_act_yDayMean,                                      &
      & hgrid, surface,                                                                             &
      & t_cf('NPP_act_yDayMean', 'mol(C) m-2(canopy) s-1', 'Actual NPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
      & prefix, suffix,                                                                             &
      & lrestart=.TRUE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would
                                                  ! result in a zero divide in the first time step of an initialized run.

    CALL mem%Add_var('gpp_ydaymean', mem%GPP_yDayMean,                                        &
      & hgrid, surface,                                                                       &
      & t_cf('GPP_yDayMean', 'mol(C) m-2(canopy) s-1', 'Mean GPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
      & prefix, suffix,                                                                       &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('root_exudates', mem%root_exudates,                   &
      & hgrid, surface,                                                    &
      & t_cf('root_exudates', 'mol(C) m-2(canopy) s-1', 'Root exudates'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('cflux_c_green_2_herb', mem%cflux_c_green_2_herb,                                          &
      & hgrid, surface,                                                                                         &
      & t_cf('cflux_c_green_2_herb', 'mol(C) m-2(canopy) s-1', 'Total carbon flux from C pool green to herbivory.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('cflux_c_green_2_herb_lg', mem%cflux_herb_2_littergreen,                                   &
      & hgrid, surface,                                                                                         &
      & t_cf('cflux_herb_2_littergreen', 'mol(C) m-2(canopy) s-1',                                              &
      &      'Carbon flux from cflux_c_green_2_herb into green litter (complementary to cflux_herb_2_atm).'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('cflux_herb_2_atm', mem%cflux_herb_2_atm,                                                  &
      & hgrid, surface,                                                                                         &
      & t_cf('cflux_herb_2_atm', 'mol(C) m-2(canopy) s-1',                                                      &
      &  'Carbon flux from cflux_c_green_2_herb into atmosphere (complementary to cflux_herb_2_littergreen).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('cflux_dist_green_2_soil', mem%cflux_dist_green_2_soil,                            &
      & hgrid, surface,                                                                                 &
      & t_cf('cflux_dist_green_2_soil', 'mol(C) m-2(canopy) s-1',                                       &
      &      'Amount of C relocated by wind damage to the green litter pools.'),                        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
      & prefix, suffix,                                                                                 &
      & output_level=MEDIUM, l_conserve_quan=.TRUE., cons_quan_type_id = FLUX_C_CQ_TYPE,                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('cflux_dist_woods_2_soil', mem%cflux_dist_woods_2_soil,                            &
      & hgrid, surface,                                                                                 &
      & t_cf('cflux_dist_woods_2_soil', 'mol(C) m-2(canopy) s-1',                                       &
      &      'Amount of C relocated by wind and fire damage to the woody litter pools..'),              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
      & prefix, suffix,                                                                                 &
      & output_level=MEDIUM,  l_conserve_quan=.TRUE., cons_quan_type_id = FLUX_C_CQ_TYPE,               &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2flux_fire_all_2_atm', mem%co2flux_fire_all_2_atm,                                        &
      & hgrid, surface,                                                                                           &
      & t_cf('co2flux_fire_all_2_atm', 'kg(CO2) m-2(canopy) s-1', 'CO2 flux into atmosphere due to fire.'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                                                           &
      & output_level=MEDIUM, l_conserve_quan=.TRUE., cons_quan_type_id = FLUX_C_CQ_TYPE,                          &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Carbon variables on tile area (_ta)
    CALL mem%Add_var('c_greenwood2litter_ta', mem%cflux_c_greenwood_2_litter_ta,                       &
      & hgrid, surface,                                                                                &
      & t_cf('c_greenwood2litter_ta', 'mol(CO2) m-2 s-1', 'Total litter flux entering the soil pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=BASIC,                                                                            &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_green_ta', mem%c_green_ta,                                        &
      & hgrid, surface,                                                                   &
      & t_cf('c_green_ta', 'mol(C) m-2', 'C-Pool for green parts of vegetation.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & loutput = .TRUE., output_level=BASIC,                                             &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_woods_ta', mem%c_woods_ta,                                                   &
      & hgrid, surface,                                                                              &
      & t_cf('c_woods_ta', 'mol(C) m-2', 'C-Pool for structural material/parts of vegetation.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & loutput = .TRUE., output_level=BASIC,                                                        &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_reserve_ta', mem%c_reserve_ta,                                                                &
      & hgrid, surface,                                                                                               &
      & t_cf('c_reserve_ta', 'mol(C) m-2', 'C-Pool for reserve carbohydrates (starches, sugars) of vegetation.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
      & prefix, suffix,                                                                                               &
      & loutput = .TRUE., output_level=BASIC,                                                                         &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_crop_harvest_ta', mem%c_crop_harvest_ta,                                  &
      & hgrid, surface,                                                                           &
      & t_cf('c_crop_harvest_ta', 'mol(C) m-2', 'C-Pool for Harvested Material from Crops.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
      & prefix, suffix,                                                                           &
      & output_level=BASIC,                                                                            &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Yasso 1: Green carbon pools
    CALL mem%Add_var('c_acid_ag1_ta', mem%c_acid_ag1_ta,                                                    &
      & hgrid, surface,                                                                                     &
      & t_cf('c_acid_ag1_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_ag1_ta', mem%c_water_ag1_ta,                                                    &
      & hgrid, surface,                                                                                       &
      & t_cf('c_water_ag1_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
      & prefix, suffix,                                                                                       &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_ag1_ta', mem%c_ethanol_ag1_ta,                                                    &
      & hgrid, surface,                                                                                           &
      & t_cf('c_ethanol_ag1_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                                                           &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_ag1_ta', mem%c_nonsoluble_ag1_ta,                                             &
      & hgrid, surface,                                                                                          &
      & t_cf('c_nonsoluble_ag1_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_acid_bg1_ta', mem%c_acid_bg1_ta,                                                    &
      & hgrid, surface,                                                                                     &
      & t_cf('c_acid_bg1_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_bg1_ta', mem%c_water_bg1_ta,                                                    &
      & hgrid, surface,                                                                                       &
      & t_cf('c_water_bg1_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
      & prefix, suffix,                                                                                       &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_bg1_ta', mem%c_ethanol_bg1_ta,                                                   &
      & hgrid, surface,                                                                                          &
      & t_cf('c_ethanol_bg1_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_bg1_ta', mem%c_nonsoluble_bg1_ta,                                            &
      & hgrid, surface,                                                                                         &
      & t_cf('c_nonsoluble_bg1_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_humus_1_ta', mem%c_humus_1_ta,                          &
      & hgrid, surface,                                                         &
      & t_cf('c_humus_1_ta', 'mol(C) m-2', 'C-Pool for humus in Yasso.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),      &
      & prefix, suffix,                                                         &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Yasso 2: Woody carbon pools
    CALL mem%Add_var('c_acid_ag2_ta', mem%c_acid_ag2_ta,                                                    &
      & hgrid, surface,                                                                                     &
      & t_cf('c_acid_ag2_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_ag2_ta', mem%c_water_ag2_ta,                                                    &
      & hgrid, surface,                                                                                       &
      & t_cf('c_water_ag2_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
      & prefix, suffix,                                                                                       &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_ag2_ta', mem%c_ethanol_ag2_ta,                                                    &
      & hgrid, surface,                                                                                           &
      & t_cf('c_ethanol_ag2_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                                                           &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_ag2_ta', mem%c_nonsoluble_ag2_ta,                                             &
      & hgrid, surface,                                                                                          &
      & t_cf('c_nonsoluble_ag2_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_acid_bg2_ta', mem%c_acid_bg2_ta,                                                    &
      & hgrid, surface,                                                                                     &
      & t_cf('c_acid_bg2_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_water_bg2_ta', mem%c_water_bg2_ta,                                                    &
      & hgrid, surface,                                                                                       &
      & t_cf('c_water_bg2_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
      & prefix, suffix,                                                                                       &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_ethanol_bg2_ta', mem%c_ethanol_bg2_ta,                                                   &
      & hgrid, surface,                                                                                          &
      & t_cf('c_ethanol_bg2_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
      & prefix, suffix,                                                                                          &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_nonsoluble_bg2_ta', mem%c_nonsoluble_bg2_ta,                                            &
      & hgrid, surface,                                                                                         &
      & t_cf('c_nonsoluble_bg2_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (belowground).'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                                                         &
      & output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_humus_2_ta', mem%c_humus_2_ta,                          &
      & hgrid, surface,                                                         &
      & t_cf('c_humus_2_ta', 'mol(C) m-2', 'C-Pool for humus in Yasso.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),      &
      & prefix, suffix,                                                         &
      & loutput = .TRUE., output_level=MEDIUM,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c_sum_veg_ta', mem%C_sum_veg_ta,                         &
      & hgrid, surface,                                                        &
      & t_cf('C_sum_veg_ta', 'mol(C) m-2', 'Sum of all vegetation C pools.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & output_level=BASIC,                                                    &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('c_sum_litter_ag_ta', mem%C_sum_litter_ag_ta,                           &
      & hgrid, surface,                                                                      &
      & t_cf('C_sum_litter_ag_ta', 'mol(C) m-2', 'Sum of all above ground litter C pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
      & prefix, suffix,                                                                      &
      & output_level=BASIC,                                                                  &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'c_sum_litter_bg_ta', mem%C_sum_litter_bg_ta,                          &
      & hgrid, surface,                                                                      &
      & t_cf('C_sum_litter_bg_ta', 'mol(C) m-2', 'Sum of all below ground litter C pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
      & prefix, suffix,                                                                      &
      & output_level=BASIC,                                                                  &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('c_sum_humus_ta', mem%C_sum_humus_ta,                     &
      & hgrid, surface,                                                        &
      & t_cf('C_sum_humus_ta', 'mol(C) m-2', 'Sum of all humus C pools.'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & output_level=BASIC,                                                    &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('c_sum_natural_ta', mem%C_sum_natural_ta,                 &
      & hgrid, surface,                                                        &
      & t_cf('C_sum_natural_ta', 'mol(C) m-2', 'Sum of all natural C pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & output_level=BASIC,                                                    &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('soil_respiration_ta', mem%soil_respiration_ta,        &
      & hgrid, surface,                                                     &
      & t_cf('soil_respiration_ta', 'mol(C) m-2 s-1', 'Soil respiration.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & output_level=BASIC,                                                 &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_flux_corr_ta', mem%NPP_flux_correction_ta,                                 &
      & hgrid, surface,                                                                              &
      & t_cf('NPP_flux_correction_ta', 'mol(C) m-2 s-1',                                             &
      &      'Amount by which the NPP rate has to be corrected to avoid a negative reserve pool.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('excess_npp_ta', mem%excess_NPP_ta,                                                          &
      & hgrid, surface,                                                                                           &
      & t_cf('excess_NPP_ta', 'mol(C) m-2  s-1', 'NPP that could not be not stored in vegetation carbon pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                                                           &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    !      CALL mem%Add_var('lai_yydaymean_ta', mem%LAI_yyDayMean_ta,              &
    !        & hgrid, surface,                                                     &
    !        & t_cf('LAI_yyDayMean_ta', '-', 'Previous previous day mean.'),          &
    !        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
    !        & prefix, suffix,                                                     &
    !        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    !     CALL mem%Add_var('lai_ydaymean_ta', mem%LAI_yDayMean_ta,                 &
    !        & hgrid, surface,                                                     &
    !        & t_cf('LAI_yDayMean_ta', '-', 'Previous day mean.'),                    &
    !        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
    !        & prefix, suffix,                                                     &
    !        & output_level=BASIC,                                                                                &
    !        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_pot_ydaymean_ta', mem%NPP_pot_yDayMean_ta,                              &
      & hgrid, surface,                                                                           &
      & t_cf('NPP_pot_yDayMean_ta', 'mol(C) m-2 s-1', 'Mean NPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
      & prefix, suffix,                                                                           &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would result
                                                  ! in a zero divide in the first time step of an initialized run.

    CALL mem%Add_var('npp_act_ydaymean_ta', mem%NPP_act_yDayMean_ta,                                &
      & hgrid, surface,                                                                             &
      & t_cf('NPP_act_yDayMean_ta', 'mol(C) m-2 s-1', 'Actual NPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
      & prefix, suffix,                                                                             &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would result
                                                  ! in a zero divide in the first time step of an initialized run.

    CALL mem%Add_var('gpp_ydaymean_ta', mem%GPP_yDayMean_ta,                             &
      & hgrid, surface,                                                                  &
      & t_cf('GPP_yDayMean_ta', 'mol(C) m-2 s-1', 'Mean GPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3


    CALL mem%Add_var('root_exudates_ta', mem%root_exudates_ta,             &
      & hgrid, surface,                                                    &
      & t_cf('root_exudates_ta', 'mol(C) m-2 s-1', 'Root exudates'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3


    CALL mem%Add_var('cflux_c_green_2_herb_ta', mem%cflux_c_green_2_herb_ta,                                   &
      & hgrid, surface,                                                                                    &
      & t_cf('cflux_c_green_2_herb_ta', 'mol(C) m-2 s-1', 'Total carbon flux from C pool green to herbivory.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                 &
      & prefix, suffix,                                                                                    &
      & output_level=MEDIUM,                                                                               &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c2greenlitter_ta', mem%cflux_dist_green_2_soil_ta,                        &
      & hgrid, surface,                                                                         &
      & t_cf('cflux_dist_green_2_soil_ta', 'mol(C) m-2 s-1', 'Carbon flux to green litter.'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      &
      & prefix, suffix,                                                                         &
      & output_level=MEDIUM,                                                                               &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('c2woodlitter_ta', mem%cflux_dist_woods_2_soil_ta,                        &
      & hgrid, surface,                                                                        &
      & t_cf('cflux_dist_woods_2_soil_ta', 'mol(C) m-2 s-1', 'Carbon flux to woody litter.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=MEDIUM,loutput=.TRUE.,                                                    &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2_l2a_fire_ta', mem%co2flux_fire_all_2_atm_ta,                                    &
      & hgrid, surface,                                                                                   &
      & t_cf('co2flux_fire_all_2_atm_ta', 'kg(CO2) m-2 s-1', 'CO2 flux from C pools to atmosphere due to fire.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                &
      & prefix, suffix, loutput=.TRUE.,                                                                   &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2_l2a_npp_ta', mem%co2flux_npp_2_atm_ta,                                                     &
      & hgrid, surface,                                                                                              &
      & t_cf('co2flux_npp_2_atm_ta', 'kg(CO2) m-2 s-1', 'Net CO2 fluxes between biosphere and atmosphere due to NPP.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
      & prefix, suffix,                                                                                              &
      & output_level=BASIC,                                                                                &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2_l2a_resp_ta', mem%co2flux_soilresp_2_atm_ta,                                     &
      & hgrid, surface,                                                                                    &
      & t_cf('co2flux_soilresp_2_atm_ta', 'kg(CO2) m-2 s-1',                                               &
      &      'Net CO2 fluxes between biosphere and atmosphere due to soil respiration.'),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                 &
      & prefix, suffix,                                                                                    &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2_l2a_herb_ta', mem%co2flux_herb_2_atm_ta,                                          &
      & hgrid, surface,                                                                                     &
      & t_cf('co2flux_herb_2_atm_ta', 'kg(CO2) m-2 s-1',                                                    &
      &      'Net CO2 fluxes between biosphere and atmosphere due to herbivory.'),                          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
      & prefix, suffix,                                                                                     &
      & output_level=BASIC,                                                                                &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('co2_l2a_npp_yday_ta', mem%co2flux_npp_2_atm_yday_ta,                                           &
      & hgrid, surface,                                                                                              &
      & t_cf('co2flux_npp_2_atm_yday_ta', 'kg(CO2) m-2 s-1', 'Net daily CO2 flux, only needed for conservation test.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
      & prefix, suffix,                                                                                              &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('ccons_calccpools', mem%cconservation_calcCpools,            &
      & hgrid, surface,                                                           &
      & t_cf('cconservation_calcCpools', 'mol(C) m-2(tile area)',                 &
      & 'C conservation test: should be zero'),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
      & prefix, suffix,                                                           &
      & output_level=BASIC,                                                       &
      & loutput=.TRUE.,lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var('c_decomp_humus_1_sum', mem%c_decomp_humus_1_sum, hgrid, surface,                                 &
      & t_cf('c_decomp_humus_1_sum', 'mol(C) m-2(canopy)', 'annual sum of carbon decomposed from YASSO humus (leaf)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix, loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('c_decomp_humus_2_sum', mem%c_decomp_humus_2_sum, hgrid, surface,                                 &
      & t_cf('c_decomp_humus_2_sum', 'mol(C) m-2(canopy)', 'annual sum of carbon decomposed from YASSO humus (wood)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix, loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('c_into_humus_1_sum', mem%c_into_humus_1_sum, hgrid, surface,                                     &
      & t_cf('c_into_humus_1_sum', 'mol(C) m-2(canopy)', 'annual sum of carbon flux to humus pool in YASSO (leaf)'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix, loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('c_into_humus_2_sum', mem%c_into_humus_2_sum, hgrid, surface,                                     &
      & t_cf('c_into_humus_2_sum', 'mol(C) m-2(canopy)', 'annual sum of carbon flux to humus pool in YASSO (wood)'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
      & prefix, suffix, loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp)

    CALL mem%Add_var('c_decomp_humus_1_sum_ta', mem%c_decomp_humus_1_sum_ta, hgrid, surface,                      &
      & t_cf('c_decomp_humus_1_sum_ta', 'mol(C) m-2', 'annual sum of carbon decomposed from YASSO humus (leaf)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix, loutput=.FALSE., lrestart=dsl4jsb_Config(CARBON_)%diag_humus_fluxes, initval_r=0.0_wp)
    CALL mem%Add_var('c_decomp_humus_2_sum_ta', mem%c_decomp_humus_2_sum_ta, hgrid, surface,                      &
      & t_cf('c_decomp_humus_2_sum_ta', 'mol(C) m-2', 'annual sum of carbon decomposed from YASSO humus (wood)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix, loutput=.FALSE., lrestart=dsl4jsb_Config(CARBON_)%diag_humus_fluxes, initval_r=0.0_wp)
    CALL mem%Add_var('c_into_humus_1_sum_ta', mem%c_into_humus_1_sum_ta, hgrid, surface,                          &
      & t_cf('c_into_humus_1_sum_ta', 'mol(C) m-2', 'annual sum of carbon flux to humus pool in YASSO (leaf)'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix, loutput=.FALSE., lrestart=dsl4jsb_Config(CARBON_)%diag_humus_fluxes, initval_r=0.0_wp)
    CALL mem%Add_var('c_into_humus_2_sum_ta', mem%c_into_humus_2_sum_ta, hgrid, surface,                          &
      & t_cf('c_into_humus_2_sum_ta', 'mol(C) m-2', 'annual sum of carbon flux to humus pool in YASSO (wood)'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix, loutput=.FALSE., lrestart=dsl4jsb_Config(CARBON_)%diag_humus_fluxes, initval_r=0.0_wp)

    ! Used in phenology if forest regrowth is active
    IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
      CALL mem%Add_var('max_green', mem%current_max_green,                   &
        & hgrid, surface,                                                    &
        & t_cf('current_max_green', 'mol(C) m-2(canopy)', 'Running maximum of the green pool'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE.)
      CALL mem%Add_var('max_green_veg_c', mem%veg_carbon_at_max_green,       &
        & hgrid, surface,                                                    &
        & t_cf('veg_carbon_at_max_green', 'mol(C) m-2(canopy)', 'Total vegetation carbon at green pool maximum'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE.)
    ENDIF

    ! Diagnostic global carbon sums for experiment monitoring - only calculated on top (i.e. box) tile
    IF (ASSOCIATED(mem%parent)) RETURN

    CALL mem%Add_var('c_sum_veg_gsum', mem%C_sum_veg_gsum,                                     &
      & hgrid, surface,                                                                        &
      & t_cf('C_sum_veg_gsum', 'PgC', 'Global sum of all vegetation C pools.'),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
    CALL mem%Add_var('c_sum_litter_ag_gsum', mem%C_sum_litter_ag_gsum,                         &
      & hgrid, surface,                                                                        &
      & t_cf('C_sum_litter_ag_gsum', 'PgC', 'Global sum of all above ground litter C pools'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
    CALL mem%Add_var('c_sum_litter_bg_gsum', mem%C_sum_litter_bg_gsum,                         &
      & hgrid, surface,                                                                        &
      & t_cf('C_sum_litter_bg_gsum', 'PgC', 'Global sum of all below ground litter C pools.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
    CALL mem%Add_var('c_sum_humus_gsum', mem%C_sum_humus_gsum,                                 &
      & hgrid, surface,                                                                        &
      & t_cf('C_sum_humus_gsum', 'PgC', 'Global sum of all humus C pools.'),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
    CALL mem%Add_var('c_sum_natural_gsum', mem%C_sum_natural_gsum,                             &
      & hgrid, surface,                                                                        &
      & t_cf('C_sum_natural_gsum', 'PgC', 'Global sum of all natural C pools.'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
    CALL mem%Add_var('npp_act_ydaymean_gsum', mem%NPP_act_yDayMean_gsum,                       &
      & hgrid, surface,                                                                        &
      & t_cf('NPP_act_yDayMean_gsum', 'PgC yr-1', 'Global actual NPP rate of the previous day'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
    CALL mem%Add_var('gpp_ydaymean_gsum', mem%GPP_yDayMean_gsum,                               &
      & hgrid, surface,                                                                        &
      & t_cf('GPP_yDayMean_gsum', 'PgC yr-1', 'Global GPP rate of the previous day'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )
    CALL mem%Add_var('soil_respiration_gsum', mem%soil_respiration_gsum,                       &
      & hgrid, surface,                                                                        &
      & t_cf('soil_respiration_gsum', 'PgC yr-1', 'Global soil respiration rate'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_groups=['mon'], lrestart=.FALSE., initval_r=0.0_wp )

  END SUBROUTINE Init_sb_carbon_memory

#endif
END MODULE mo_sb_carbon_memory_class
