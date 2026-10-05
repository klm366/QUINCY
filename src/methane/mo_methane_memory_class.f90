!> Contains the memory class for the methane production and transport
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
!>
!> author
!>  Stiig Wilkenskjeld, MPI-M Hamburg
!>
!> Revision History
!> Adopted from JSBACH3 code by Thomas Kleinen, 2023-03-17
!>

MODULE mo_ch4_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: lp
  USE mo_exception,              ONLY: finish
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE, BARE_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_process_class,      ONLY: TOP_, CH4_

  ! Use of prcesses in this module
  dsl4jsb_Use_processes TOP_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_ch4_memory, max_no_of_vars, t_chem_specie, t_chem_ptrs, t_timeseries

  INTEGER, PARAMETER :: max_no_of_vars = 130

  ! ======================================================================================================= !
  !>Type definitions for methane production and transport memory
  !>
  TYPE t_chem_specie
    TYPE(t_jsb_var_real2d) ::   &
      & prod_out,               &
      & oxid_out,               &
      & diff_out,               &
      & other_out,              & ! Ebulition or "other"
      & plant_out,              &
      & in_out,                 &
      & excess,                 &
      & total_out,              &
      & total_initial,          & ! StW: Initial for ts or run? In latter case ic?
      & total_final
    TYPE(t_jsb_var_real2d) ::   &
      & test2d

    ! Per layer
    TYPE(t_jsb_var_real3d) ::   &
      & dens_sl                   ! Current layer density [mol m-2]

    ! Debug per layer
    TYPE(t_jsb_var_real3d) ::   &
      & prod_out_sl,            &
      & oxid_out_sl,            &
      & other_out_sl,           & ! Ebulition or "other"
      & plant_out_sl
    TYPE(t_jsb_var_real3d) ::   &
      & test3d

    CHARACTER(len=8) :: specie_name

  END TYPE t_chem_specie

  TYPE t_timeseries
    REAL(lp), POINTER :: atm_vmr(:)
    INTEGER           :: atm_year_min
    INTEGER           :: atm_year_max
  END TYPE t_timeseries

  TYPE t_chem_ptrs
    REAL(lp), DIMENSION(:), POINTER ::     &
      & prod_out,                          &
      & oxid_out,                          &
      & diff_out,                          &
      & other_out,                         & ! Ebulition or "other"
      & plant_out,                         &
      & in_out,                            &
      & excess,                            &
      & total_out,                         &
      & total_initial,                     & ! StW: Initial for ts or run? In latter case ic?
      & total_final
    REAL(lp), DIMENSION(:), POINTER ::     &
      test2d

    ! Per layer
    REAL(lp), DIMENSION(:,:), POINTER ::   &
      & dens_sl                              ! Current layer density [mol m-2]

    ! Debug per layer
    REAL(lp), DIMENSION(:,:), POINTER ::   &
      & prod_out_sl,            &
      & oxid_out_sl,            &
      & other_out_sl,           & ! Ebulition or "other"
      & plant_out_sl
    REAL(lp), DIMENSION(:,:), POINTER ::   &
      test3d

    CHARACTER(len=8) :: specie_name

  END TYPE t_chem_ptrs

  TYPE, EXTENDS(t_jsb_memory) :: t_ch4_memory

    ! General types for chemical species in oxic/anoxic environments
    TYPE(t_chem_specie) ::      &
      & CH4_dry,                &
      & CH4_wet,                &
!      & CH4_tot,                &
      & CO2_dry,                &
      & CO2_wet,                &
      & O2_dry,                 &
      & O2_wet

    ! Basic variables and constants
    TYPE(t_jsb_var_real2d) ::   &
      & CH4_conc,               & ! StW: Collect for all tiles?
      & CO2_conc,               &
      & O2_conc,                &
      & fire_em_NO,             &
      & fire_em_CO,             &
      & fire_em_iso,            &
      & mxlay_above_rock,       & ! Contant throughout run
      & tillers,                &
      & termite_biomass,        &
!
      & CH4_prod_ratio_out_wet, &
!
      & CH4_out_net,            &
      & CO2_out_net,            &
      & CH4_out_fire,           &
      & CH4_out_termite,        &
      & CH4_out_herbivore,      &
!
      & anox_prod_in,           & ! StW: Production of what?
      & oxic_prod_in              ! StW: Production of what?

    ! Variables for coupling to the atmospheric model
    TYPE (t_jsb_var_real2d) ::  &
      & ch4_2_atm,              & ! Total methane surface emissions
      & co_2_atm,               & ! Total carbon monoxide emissions
      & no_2_atm,               & ! Total nitric emissions
      & iso_2_atm                 ! Total isoprene emissions

    ! Data for obtaining time dependent atmospheric GHG from file
    TYPE(t_timeseries) ::       &
      & GHG_atm_CH4,            &
      & GHG_atm_CO2,            &
      & GHG_atm_O2

    ! per layer
    TYPE(t_jsb_var_real3d) :: &
      & root_dist_sl

    REAL(lp), POINTER ::                &
      & layer_midpoint_dz(:) => NULL()

  CONTAINS

    PROCEDURE :: Init => Init_methane_memory

  END TYPE t_ch4_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_ch4_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !>initialize memory (variables) for the process: methane
  !>
  SUBROUTINE Init_methane_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, TSTEP_INSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model
    dsl4jsb_Use_config(CH4_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_ch4_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)            :: prefix
    CHARACTER(len=*),    INTENT(in)            :: suffix
    INTEGER,             INTENT(in)            :: lct_ids(:)
    INTEGER,             INTENT(in)            :: lib_id        !< ID of the tile's primary land cover type (e.g. LAND_TYPE)
    INTEGER,             INTENT(in)            :: model_id
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(CH4_)
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: surface, soil_w              ! Vertical grids
    INTEGER                     :: table
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_methane_memory'
    INTEGER :: default_output_level = BASIC
    ! ----------------------------------------------------------------------------------------------------- !
    IF (SIZE(lct_ids) < 1) CALL finish(routine, 'lcd_ids not specified')
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
    table         = tables(1)
    hgrid         => Get_grid(mem%grid_id)
    surface       => Get_vgrid('surface')
    soil_w        => Get_vgrid('soil_depth_water')

    dsl4jsb_Get_Config(CH4_)
    !> !$acc enter data pcreate(mem)

    ! ----------------------------------------------------------------------------------------------------- !
    ! jsbach4 memory
    !

    ALLOCATE(mem%layer_midpoint_dz(soil_w%n_levels))
    !$acc enter data create(mem%layer_midpoint_dz)
    !$acc enter data attach(mem%layer_midpoint_dz)

    NULLIFY(mem%GHG_atm_CH4%atm_vmr)
    NULLIFY(mem%GHG_atm_CO2%atm_vmr)
    NULLIFY(mem%GHG_atm_O2%atm_vmr)

    CALL Init_chem_specie('CH4_dry',mem%CH4_dry,.FALSE.)
    CALL Init_chem_specie('CH4_wet',mem%CH4_wet,.TRUE.,'ebul')
    CALL Init_chem_specie('CO2_dry',mem%CO2_dry,.FALSE.)
    CALL Init_chem_specie('CO2_wet',mem%CO2_wet,.FALSE.)
    CALL Init_chem_specie('O2_dry', mem%O2_dry, .TRUE.)
    CALL Init_chem_specie('O2_wet', mem%O2_wet, .TRUE.)
!    CALL Init_chem_specie('CH4_tot',mem%CH4_tot,.TRUE.,'ebul')

    CALL mem%Add_var( 'CH4_conc', mem%CH4_conc,                              &
      & hgrid, surface,                                                      &
      & t_cf('CH4_conc', 'mol m-3', ''),                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix, lrestart = .FALSE.,                                  &
      & output_level = FULL, l_aggregate_all=.FALSE.,                        &
      & initval_r=7.706281e-5_lp, isteptype=TSTEP_CONSTANT ) ! May be variable in future

    CALL mem%Add_var( 'CO2_conc', mem%CO2_conc,                              &
      & hgrid, surface,                                                      &
      & t_cf('CO2_conc', 'mol m-3', ''),                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix, lrestart = .FALSE.,                                  &
      & output_level = FULL, l_aggregate_all=.FALSE.,                        &
      & initval_r=1.6349598227e-2_lp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'O2_conc', mem%O2_conc,                                &
      & hgrid, surface,                                                      &
      & t_cf('O2_conc', 'mol m-3', ''),                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix, lrestart = .FALSE.,                                  &
      & output_level = FULL, l_aggregate_all=.FALSE.,                        &
      & initval_r=1.6349598227e-2_lp, isteptype=TSTEP_CONSTANT ) ! May be variable in future

    CALL mem%Add_var( 'mxlay_above_rock', mem%mxlay_above_rock,              &
      & hgrid, surface,                                                      &
      & t_cf('mxlay_above_rock', '', ''),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix, lrestart = .FALSE.,                                  &
      & output_level = default_output_level, l_aggregate_all=.FALSE.,        &
      & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'tillers', mem%tillers,                                &
      & hgrid, surface,                                                      &
      & t_cf('tillers', '', ''),                                             & ! unit?
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix, lrestart = .FALSE.,                                  &
      & output_level = default_output_level, l_aggregate_all=.TRUE.,         & ! .FALSE.?
      & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'CH4_prod_ratio_out_wet', mem%CH4_prod_ratio_out_wet,  &
      & hgrid, surface,                                                      &
      & t_cf('CH4_prod_ratio_out_wet', '', ''),                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix, lrestart = .FALSE.,                                  &
      & output_level = default_output_level, l_aggregate_all=.FALSE.,        &
      & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'CH4_out_net', mem%CH4_out_net,                        &
      & hgrid, surface,                                                      &
      & t_cf('CH4_out_net', '', ''),                                         & ! unit?
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix, lrestart = .FALSE.,                                  &
      & output_level = FULL, l_aggregate_all=.TRUE.,                         &
      & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'CO2_out_net', mem%CO2_out_net,                        &
      & hgrid, surface,                                                      &
      & t_cf('CO2_out_net', '', ''),                                         & ! unit?
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix, lrestart = .FALSE.,                                  &
      & output_level = FULL, l_aggregate_all=.TRUE.,                         &
      & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

    IF (dsl4jsb_Config(CH4_)%do_fire) THEN
      CALL mem%Add_var( 'CH4_out_fire', mem%CH4_out_fire,                      &
        & hgrid, surface,                                                      &
        & t_cf('CH4_out_fire', 'mol m-2 s-1', ''),                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'fire_em_NO', mem%fire_em_NO,                          &
        & hgrid, surface,                                                      &
        & t_cf('fire_em_NO', 'kg m-2 s-1', ''),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'fire_em_CO', mem%fire_em_CO,                          &
        & hgrid, surface,                                                      &
        & t_cf('fire_em_CO', 'kg m-2 s-1', ''),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'fire_em_iso', mem%fire_em_iso,                        &
        & hgrid, surface,                                                      &
        & t_cf('fire_em_iso', 'g m-2 s-1', ''),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )
    ENDIF

    IF (dsl4jsb_Config(CH4_)%do_termite) THEN
      CALL mem%Add_var( 'CH4_out_termite', mem%CH4_out_termite,                &
        & hgrid, surface,                                                      &
        & t_cf('CH4_out_termite', 'mol m-2 s-1', ''),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'termite_biomass', mem%termite_biomass,                &
        & hgrid, surface,                                                      &
        & t_cf('termite_biomass', 'kg m-2', ''),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )
    ENDIF

    IF (dsl4jsb_Config(CH4_)%do_herbivore) THEN
      CALL mem%Add_var( 'CH4_out_herbivore', mem%CH4_out_herbivore,            &
        & hgrid, surface,                                                      &
        & t_cf('CH4_out_herbivore', 'mol m-2 s-1', ''),                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )
    ENDIF

    IF (dsl4jsb_Config(CH4_)%atm_couple) THEN
      CALL mem%Add_var( 'ch4_2_atm', mem%ch4_2_atm,                            &
        & hgrid, surface,                                                      &
        & t_cf('Total terrestrial methane emisions', 'kg(CH4) m-2 s-1', ''),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .TRUE., lrestart_cont = .TRUE.,           &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,                         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'co_2_atm', mem%co_2_atm,                              &
        & hgrid, surface,                                                      &
        & t_cf('Total terrestrial carbon monoxide emisions', 'kg(CO) m-2 s-1', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .TRUE., lrestart_cont = .TRUE.,           &
        & output_level =default_output_level, l_aggregate_all=.TRUE.,                         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'no_2_atm', mem%no_2_atm,                              &
        & hgrid, surface,                                                      &
        & t_cf('Total terrestrial nitric emisions', 'kg(NO) m-2 s-1', ''),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .TRUE., lrestart_cont = .TRUE.,           &
        & output_level =default_output_level, l_aggregate_all=.TRUE.,                         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'iso_2_atm', mem%iso_2_atm,                            &
        & hgrid, surface,                                                      &
        & t_cf('Total terrestrial isoprene emisions', 'kg m-2 s-1', ''),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .TRUE., lrestart_cont = .TRUE.,           &
        & output_level =default_output_level, l_aggregate_all=.TRUE.,                         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

    ENDIF

    CALL mem%Add_var( 'root_dist_sl', mem%root_dist_sl,                      &
      & hgrid, soil_w,                                                       &
      & t_cf('root_dist_sl', '', ''),                                        & ! unit?
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix, lrestart = .FALSE.,                                  &
      & output_level = default_output_level, l_aggregate_all=.FALSE.,        &
      & initval_r=0.0_lp, isteptype=TSTEP_CONSTANT ) ! By dynamic roots: TSTEP_INSTANT

  CONTAINS

    SUBROUTINE Init_chem_specie(specie_name, specie, use_other, other_name)

      USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
      USE mo_jsb_grid,          ONLY: Get_vgrid
      USE mo_exception,         ONLY: message_text, message
      USE mo_jsb_io,            ONLY: TSTEP_CONSTANT, TSTEP_INSTANT
      dsl4jsb_Use_config(CH4_)

      CHARACTER(len=*),             INTENT(in)     :: specie_name
      TYPE (t_chem_specie), TARGET, INTENT(inout)  :: specie
      LOGICAL,          OPTIONAL,   INTENT(in)     :: use_other
      CHARACTER(len=*), OPTIONAL,   INTENT(in)     :: other_name

      TYPE(t_jsb_vgrid), POINTER  :: layer               ! Vertical grids
      CHARACTER(len=256) :: str
      LOGICAL :: uother

      dsl4jsb_Def_config(CH4_)
      dsl4jsb_Get_config(CH4_)

      uother = .TRUE.; IF (PRESENT(use_other)) uother = use_other

      surface       => Get_vgrid('surface')
      layer         => Get_vgrid('soil_depth_water')

      specie%specie_name = TRIM(specie_name)

      CALL mem%Add_var( specie_name//'_prod_out', specie%prod_out,             &
        & hgrid, surface,                                                      &
        & t_cf(specie_name//'_prod_out', 'mol m-2 s-1', ''),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( specie_name//'_oxid_out', specie%oxid_out,             &
        & hgrid, surface,                                                      &
        & t_cf(specie_name//'_oxid_out', 'mol m-2 s-1', ''),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( specie_name//'_diff_out', specie%diff_out,             &
        & hgrid, surface,                                                      &
        & t_cf(specie_name//'_diff_out', 'mol m-2 s-1', ''),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( specie_name//'_plant_out', specie%plant_out,           &
        & hgrid, surface,                                                      &
        & t_cf(specie_name//'_plant_out', 'mol m-2 s-1', ''),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level,                                 &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( specie_name//'_in_out', specie%in_out,                 &
        & hgrid, surface,                                                      &
        & t_cf(specie_name//'_in_out', 'mol m-2 s-1', ''),                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level,                                 &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

!      CALL mem%Add_var( specie_name//'_excess_transp', specie%excess,          &
!        & hgrid, surface,                                                      &
!        & t_cf(specie_name//'_excess_transp', 'mol m-2 s-1', ''),              &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
!        & prefix, suffix, lrestart = .FALSE.,                                  &
!        & output_level = FULL, l_aggregate_all=.TRUE.,                         &
!        & lrestart = .FALSE.,                                                  &
!        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      IF (dsl4jsb_Config(CH4_)%dbg_process) THEN

        CALL mem%Add_var( specie_name//'_prod_out_sl', specie%prod_out_sl,       &
          & hgrid, layer,                                                        &
          & t_cf(specie_name//'_prod_out_sl', 'mol m-2 s-1', ''),                &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix, lrestart = .FALSE.,                                  &
          & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
          & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

        CALL mem%Add_var( specie_name//'_oxid_out_sl', specie%oxid_out_sl,       &
          & hgrid, layer,                                                        &
          & t_cf(specie_name//'_oxid_out_sl', 'mol m-2 s-1', ''),                &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix, lrestart = .FALSE.,                                  &
          & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
          & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

        CALL mem%Add_var( specie_name//'_plant_out_sl', specie%plant_out_sl,     &
          & hgrid, layer,                                                        &
          & t_cf(specie_name//'_plant_out_sl', 'mol m-2 s-1', ''),               &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix, lrestart = .FALSE.,                                  &
          & output_level = default_output_level,                                 &
          & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

        CALL mem%Add_var( specie_name//'_test3d', specie%test3d,                 &
          & hgrid, layer,                                                        &
          & t_cf(specie_name//'_test3d', '???', ''),                             &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix, lrestart = .FALSE.,                                  &
          & output_level = default_output_level,                                 &
          & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

        CALL mem%Add_var( specie_name//'_test2d', specie%test2d,                 &
          & hgrid, surface,                                                      &
          & t_cf(specie_name//'_test2d', '???', ''),                             &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix, lrestart = .FALSE.,                                  &
          & output_level = default_output_level,                                 &
          & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      ENDIF

      IF (uother) THEN
        IF (PRESENT(other_name)) THEN
          str = specie_name//'_'//other_name//'_out'
        ELSE
          str = specie_name//'_other_out'
        ENDIF
        CALL mem%Add_var( TRIM(str), specie%other_out,                           &
          & hgrid, surface,                                                      &
          & t_cf(TRIM(str), 'mol m-2 s-1', ''),                                  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix, lrestart = .FALSE.,                                  &
          & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
          & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

        IF (dsl4jsb_Config(CH4_)%dbg_process) THEN
          CALL mem%Add_var( TRIM(str)//'_sl', specie%other_out_sl,                 &
            & hgrid, layer,                                                        &
            & t_cf(TRIM(str)//'_sl', 'mol m-2 s-1', ''),                           &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
            & prefix, suffix, lrestart = .FALSE.,                                  &
            & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
            & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )
        ENDIF
      ENDIF

      CALL mem%Add_var( specie_name//'_total_out', specie%total_out,           &
        & hgrid, surface,                                                      &
        & t_cf(specie_name//'_total_out', 'mol m-2 s-1', ''),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( specie_name//'_total_initial', specie%total_initial,   &
        & hgrid, surface,                                                      &
        & t_cf(specie_name//'_total_initial', 'mol m-2', ''),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = FULL, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( specie_name//'_total_final', specie%total_final,       &
        & hgrid, surface,                                                      &
        & t_cf(specie_name//'_total_final', 'mol m-2', ''),                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( specie_name//'_sl', specie%dens_sl,                    &
        & hgrid, layer,                                                        &
        & t_cf(specie_name//'_sl', 'mol m-2', ''),                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_lp, isteptype=TSTEP_INSTANT )

    END SUBROUTINE Init_chem_specie

  END SUBROUTINE Init_methane_memory

#endif
END MODULE mo_ch4_memory_class
