!> Contains the memory class for the lex process.
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
!>#### Contains memory definitions for lex
!>
MODULE mo_lex_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist,            ONLY: BASIC, MEDIUM, FULL

  ! Use of prcesses in this module
  dsl4jsb_Use_processes LEX_

  ! Use process configurations
  dsl4jsb_Use_config(LEX_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lex_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_lex_memory

    ! Common variables independent of lct_type
    TYPE(t_jsb_var_real2d) :: &
      & slope_eff,            &  !< effective slope (micro- and macro-scale) [/]
      & fact_lag_srf,         &  !< effective lag for lateral surface fluxes [s-1]
      & fact_lag_blg,         &  !< effective lag for lateral sub-surface fluxes [s-1]
      & wtr_runoff_in,        &  !< surface runoff received from adjacent tiles [kg m-2 s-1]
      & wtr_infil_in,         &  !< infiltration from adjacent tiles [kg m-2 s-1]  
      & weq_snowblow_in,      &  !< snowblow from adjacent tiles [kg m-2 s-1]  
      & wtr_runoff_out,       &  !< surface runoff received by adjacent tiles [kg m-2 s-1] 
      & wtr_infil_out,        &  !< infiltration to adjacent tiles [kg m-2 s-1]  
      & weq_snowblow_out         !< snowblow to adjacent tiles [kg m-2 s-1]

    TYPE(t_jsb_var_real3d) :: &
      & h_offset_ct,          &  !< difference in surface elevation between adjacent tiles [m]
      & length_scale_ct,      &  !< char. length scale of srf. heterogeneity [m]
      & length_iface_ct,      &  !< char. interface length between adjacent tiles [m m-2]
      & microcon_grvt_ct,     &  !< micro-scale connection fraction with surface offset [/]
      & microcon_dfsn_ct,     &  !< micro-scale connection fraction [/]
      & macrocon_grvt_ct,     &  !< macro-scale connection fraction [/]
      & fract_latflow_srf_ct, &  !< fraction of lateral (srf) flow to adjacent tile [/]
      & fract_latflow_blg_ct, &  !< fraction of lateral (below ground) flow to adjacent tile [/]  
      & micro_slope_ct,       &  !< micro-scale slope between connected tiles [/] 
      & wtr_drainage_in_sl,   &  !< drainage received from adjacent tiles [kg m-2 s-1] 
      & wtr_diffusion_in_sl,  &  !< diffusion (water) received from adjacent tile [kg m-2 s-1]
      & wtr_drainage_out_sl,  &  !< drainage received by adjacent tiles [kg m-2 s-1] 
      & wtr_diffusion_out_sl, &  !< diffusion (water) received by adjacent tile  [kg m-2 s-1] 
      & nrg_diffusion_in_sl,  &  !< diffusion (energy) received from adjacent tile [J m-2 s-1]
      & nrg_diffusion_out_sl, &  !< diffusion (energy) received by adjacent tile [J m-2 s-1]
      & nrg_convection_in_sl, &  !< convection (energy contained in water) received from adjacent tile [J m-2 s-1]
      & nrg_convection_out_sl    !< convection (energy contained in water) received by adjacent tile [J m-2 s-1]

  CONTAINS
    PROCEDURE :: Init => Init_lex_memory
  END TYPE t_lex_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lex_memory_class'

CONTAINS

  SUBROUTINE Init_lex_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_lex_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    dsl4jsb_Def_config(LEX_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: tiles                        ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: soil_w                       ! Vertical grid

    INTEGER :: table
    TYPE(t_grib2) :: grib2_desc

    LOGICAL  :: l_lex

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_lex_memory'

    model => Get_model(model_id)

    dsl4jsb_Get_Config(LEX_)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    tiles   => Get_vgrid('connecttiles')
    soil_w  => Get_vgrid('soil_depth_water')

    ! Common variables

    IF (TRIM(mem%owner_tile_name) == 'box') THEN
      grib2_desc = t_grib2(1,0,202, grib_bits)
    ELSE
      grib2_desc = t_grib2(255, 255, 255, grib_bits)
    END IF

    l_lex = dsl4jsb_Config(LEX_)%active

    ! add your common variables here e.g:
    CALL mem%Add_var( 'slope_eff', mem%slope_eff,                                     &
        & hgrid, surface,                                                             &
        & t_cf('slope_eff', '/', 'effective slope (micro- and macro-scale)'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),          &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                     &
        & loutput=l_lex, output_level=FULL,                                           &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'fact_lag_srf', mem%fact_lag_srf,                               &
        & hgrid, surface,                                                             &
        & t_cf('fact_lag_srf', 's-1', 'effective lag for lateral surface fluxes'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),          &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                     &
        & loutput=l_lex, output_level=FULL,                                           &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'fact_lag_blg', mem%fact_lag_blg,                               &
        & hgrid, surface,                                                             &
        & t_cf('fact_lag_blg', '/', 'effective lag for lateral sub-surface fluxes'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),          &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                     &
        & loutput=l_lex, output_level=FULL,                                           &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'wtr_runoff_out', mem%wtr_runoff_out,                                   &
        & hgrid, surface,                                                                     &
        & t_cf('wtr_runoff_out', 'kg m-2 s-1', 'surface runoff received by adjacent tiles'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                             &
        & loutput=l_lex, output_level=BASIC,                                                  &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'wtr_runoff_in', mem%wtr_runoff_in,                                      &
        & hgrid, surface,                                                                      &
        & t_cf('wtr_runoff_in', 'kg m-2 s-1', 'surface runoff received from adjacent tiles'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                              &
        & loutput=l_lex, output_level=BASIC,                                                   &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'wtr_infil_out', mem%wtr_infil_out,                              &
        & hgrid, surface,                                                              &
        & t_cf('wtr_infil_out', 'kg m-2 s-1', 'infiltration to adjacent tiles'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                      &
        & loutput=l_lex, output_level=BASIC,                                           &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'wtr_infil_in', mem%wtr_infil_in,                                &
        & hgrid, surface,                                                              &
        & t_cf('wtr_infil_in', 'kg m-2 s-1', 'infiltration from adjacent tiles'),      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                      &
        & loutput=l_lex, output_level=BASIC,                                           &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'weq_snowblow_out', mem%weq_snowblow_out,                        &
        & hgrid, surface,                                                              &
        & t_cf('weq_snowblow_out', 'kg m-2 s-1', 'snowblow to adjacent tiles'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                      &
        & loutput=l_lex, output_level=BASIC,                                           &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'weq_snowblow_in', mem%weq_snowblow_in,                          &
        & hgrid, surface,                                                              &
        & t_cf('weq_snowblow_in', 'kg m-2 s-1', 'snowblow from adjacent tiles'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                      &
        & loutput=l_lex, output_level=BASIC,                                           &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'h_offset_ct', mem%h_offset_ct,                                          &
        & hgrid, tiles,                                                                        &
        & t_cf('h_offset_ct', '/', 'difference in surface elevation between adjacent tiles'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                              &
        & loutput=l_lex, output_level=FULL,                                                    &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'length_scale_ct', mem%length_scale_ct,                        &
        & hgrid, tiles,                                                              &
        & t_cf('length_scale_ct', 'm', 'char. length scale of srf. heterogeneity'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                    &
        & loutput=l_lex, output_level=FULL,                                          &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'length_iface_ct', mem%length_iface_ct,                                 &
        & hgrid, tiles,                                                                       &
        & t_cf('length_iface_ct', 'm m-2', 'char. interface length between adjacent tiles'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                             &
        & loutput=l_lex, output_level=FULL,                                                   &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'microcon_grvt_ct', mem%microcon_grvt_ct,                                  &
        & hgrid, tiles,                                                                          &
        & t_cf('microcon_grvt_ct', '/', 'micro-scale connection fraction with surface offset'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                                &
        & loutput=l_lex, output_level=FULL,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'microcon_dfsn_ct', mem%microcon_dfsn_ct,               &
        & hgrid, tiles,                                                       &
        & t_cf('microcon_dfsn_ct', '/', 'micro-scale connection fraction'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        & prefix, suffix, l_aggregate_all=.TRUE.,                             &
        & loutput=l_lex, output_level=FULL,                                   &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'macrocon_grvt_ct', mem%macrocon_grvt_ct,               &
        & hgrid, tiles,                                                       &
        & t_cf('macrocon_grvt_ct', '/', 'macro-scale connection fraction'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        & prefix, suffix, l_aggregate_all=.TRUE.,                             &
        & loutput=l_lex, output_level=FULL,                                   &
        & lrestart=.TRUE., initval_r=0.0_wp)


    CALL mem%Add_var( 'fract_latflow_srf_ct', mem%fract_latflow_srf_ct,                          &
        & hgrid, tiles,                                                                          &
        & t_cf('fract_latflow_srf_ct', '/', 'fraction of lateral (srf) flow to adjacent tile'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                                &
        & loutput=l_lex, output_level=FULL,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'fract_latflow_blg_ct', mem%fract_latflow_blg_ct,                          &
        & hgrid, tiles,                                                                          &
        & t_cf('fract_latflow_blg_ct', '/',                                                      &
        & 'fraction of lateral (below ground) flow to adjacent tile'),                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                                &
        & loutput=l_lex, output_level=FULL,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'micro_slope_ct', mem%micro_slope_ct,                          &
        & hgrid, tiles,                                                              &
        & t_cf('micro_slope_ct', '/', 'micro-scale slope between connected tiles'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                    &
        & loutput=l_lex, output_level=FULL,                                          &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'wtr_drainage_out_sl', mem%wtr_drainage_out_sl,                        &
        & hgrid, soil_w,                                                                     &
        & t_cf('wtr_drainage_out_sl', 'kg m-2 s-1', 'drainage received by adjacent tiles'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                            &
        & loutput=l_lex, output_level=BASIC,                                                 &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'wtr_drainage_in_sl', mem%wtr_drainage_in_sl,                           &
        & hgrid, soil_w,                                                                      &
        & t_cf('wtr_drainage_in_sl', 'kg m-2 s-1', 'drainage received from adjacent tiles'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                             &
        & loutput=l_lex, output_level=BASIC,                                                  &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'wtr_diffusion_out_sl', mem%wtr_diffusion_out_sl,                          &
        & hgrid, soil_w,                                                                         &
        & t_cf('wtr_diffusion_out_sl', 'kg m-2 s-1',                                             &
        & 'diffusion (water) received by adjacent tile'),                                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                                &
        & loutput=l_lex, output_level=BASIC,                                                     &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'wtr_diffusion_in_sl', mem%wtr_diffusion_in_sl,                            &
        & hgrid, soil_w,                                                                         &
        & t_cf('wtr_diffusion_in_sl', 'kg m-2 s-1',                                              &
        & 'diffusion (water) received from adjacent tile'),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                                &
        & loutput=l_lex, output_level=BASIC,                                                     &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'nrg_diffusion_out_sl', mem%nrg_diffusion_out_sl,                          &
        & hgrid, soil_w,                                                                         &
        & t_cf('nrg_diffusion_out_sl', 'J m-2 s-1',                                              &
        & 'diffusion (energy) received by adjacent tile'),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                                &
        & loutput=l_lex, output_level=BASIC,                                                     &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'nrg_diffusion_in_sl', mem%nrg_diffusion_in_sl,                            &
        & hgrid, soil_w,                                                                         &
        & t_cf('nrg_diffusion_in_sl', 'J m-2 s-1',                                               &
        & 'diffusion (energy) received from adjacent tile'),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                                &
        & loutput=l_lex, output_level=BASIC,                                                     &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'nrg_convection_out_sl', mem%nrg_convection_out_sl,                        &
        & hgrid, soil_w,                                                                         &
        & t_cf('nrg_convection_out_sl', 'J m-2 s-1',                                             &
        & 'convection (energy contained in water) received by adjacent tile'),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                                &
        & loutput=l_lex, output_level=BASIC,                                                     &
        & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'nrg_convection_in_sl', mem%nrg_convection_in_sl,                          &
        & hgrid, soil_w,                                                                         &
        & t_cf('nrg_convection_in_sl', 'J m-2 s-1',                                              &
        & 'convection (energy contained in water) received from adjacent tile'),                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix, l_aggregate_all=.TRUE.,                                                &
        & loutput=l_lex, output_level=BASIC,                                                     &
        & lrestart=.TRUE., initval_r=0.0_wp)

  END SUBROUTINE Init_lex_memory

#endif
END MODULE mo_lex_memory_class
