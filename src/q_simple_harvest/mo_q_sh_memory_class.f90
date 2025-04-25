!> QUINCY simple harvest process memory
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
!>#### definition and init of (memory) variables for the quincy simple harvest process
!>
MODULE mo_q_sh_memory_class
#ifndef __NO_QUINCY__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  dsl4jsb_Use_processes VEG_
  dsl4jsb_Use_config(VEG_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_q_sh_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 10

  ! ======================================================================================================= !
  !> Type definition for q_sh memory
  !>
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_q_sh_memory

    TYPE(t_jsb_var_real2d) ::        &
      & fract_harvest,               & !< Area harvested from tile as fraction of the grid-cell (unitless)
      & fract_wood_to_slash,         & !< Fraction of harvested wood allocated to slash (unitless)
      & fract_wood_to_pp_fuel,       & !< Fraction of harvested wood allocated to the fuel product pool (unitless)
      & fract_wood_to_pp_paper,      & !< Fraction of harvested wood allocated to the paper pp (unitless)
      & fract_wood_to_pp_fiberboard, & !< Fraction of harvested wood allocated to the fiberboard pp (unitless)
      & fract_wood_to_pp_oirw,       & !< Fraction of harvested wood allocated to other industrial roundwood (unitless)
      & fract_wood_to_pp_pv,         & !< Fraction of harvested wood allocated to the plywood and veneer pp (unitless)
      & fract_wood_to_pp_sawnwood      !< Fraction of harvested wood allocated to the sawnwood pp (unitless)

  CONTAINS
    PROCEDURE :: Init => Init_q_sh_memory
  END TYPE t_q_sh_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_sh_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !> initialize memory (variables) for the process: q_sh
  !>
  !>
  SUBROUTINE Init_q_sh_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)
    USE mo_jsb_model_class,     ONLY: t_jsb_model
    USE mo_jsb_class,           ONLY: Get_model
    USE mo_jsb_varlist,         ONLY: BASIC, MEDIUM, FULL, NONE
    USE mo_jsb_io,              ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
    USE mo_quincy_output_class, ONLY: unitless
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_sh_memory),    INTENT(inout), TARGET :: mem
    CHARACTER(len=*),        INTENT(in)            :: prefix
    CHARACTER(len=*),        INTENT(in)            :: suffix
    INTEGER,                 INTENT(in)            :: lct_ids(:)
    INTEGER,                 INTENT(in)            :: lib_id     !< id of primary lct in lctlib
    INTEGER,                 INTENT(in)            :: model_id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),  POINTER :: model                      !< model
    TYPE(t_jsb_grid),   POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid),  POINTER :: surface                      ! Vertical grid
    INTEGER                     :: table                        ! ...
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_q_sh_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    model                 => Get_model(model_id)
    table                 = tables(1)
    hgrid                 => Get_grid(mem%grid_id)
    surface               => Get_vgrid('surface')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(VEG_)

    ! create memory at tiles of LAND_TYPE or VEG_TYPE
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN

      CALL mem%Add_var('fract_harvest', mem%fract_harvest, &
        & hgrid, surface, &
        & t_cf('fract_harvest', unitless, 'Area harvested from tile as fraction of the grid-cell'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & lrestart = .TRUE., &
        & lrestart_cont = .TRUE., &
        & initval_r = 0.0_wp)

      IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
        CALL mem%Add_var('fract_wood_to_slash', mem%fract_wood_to_slash, &
          & hgrid, surface, &
          & t_cf('fract_wood_to_slash', unitless, 'Fraction of harvested wood allocated to slash'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = NONE, &
          & lrestart = .TRUE., &
          & lrestart_cont = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('fract_wood_to_pp_fuel', mem%fract_wood_to_pp_fuel, &
          & hgrid, surface, &
          & t_cf('fract_wood_to_pp_fuel', unitless, 'Fraction of harvested wood allocated to the fuel product pool'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = NONE, &
          & lrestart = .TRUE., &
          & lrestart_cont = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('fract_wood_to_pp_paper', mem%fract_wood_to_pp_paper, &
          & hgrid, surface, &
          & t_cf('fract_wood_to_pp_paper', unitless, 'Fraction of harvested wood allocated to the paper product pool'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = NONE, &
          & lrestart = .TRUE., &
          & lrestart_cont = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('fract_wood_to_pp_fiberboard', mem%fract_wood_to_pp_fiberboard, &
          & hgrid, surface, &
          & t_cf('fract_wood_to_pp_fiberboard', unitless, 'Fraction of harvested wood allocated to the fiberboard product pool'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = NONE, &
          & lrestart = .TRUE., &
          & lrestart_cont = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('fract_wood_to_pp_oirw', mem%fract_wood_to_pp_oirw, &
          & hgrid, surface, &
          & t_cf('fract_wood_to_pp_oirw', unitless, 'Fraction of harvested wood allocated to the other industrial roundwood pp'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = NONE, &
          & lrestart = .TRUE., &
          & lrestart_cont = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('fract_wood_to_pp_pv', mem%fract_wood_to_pp_pv, &
          & hgrid, surface, &
          & t_cf('fract_wood_to_pp_pv', unitless, 'Fraction of harvested wood allocated to the plywood and veneer product pool'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = NONE, &
          & lrestart = .TRUE., &
          & lrestart_cont = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('fract_wood_to_pp_sawnwood', mem%fract_wood_to_pp_sawnwood, &
          & hgrid, surface, &
          & t_cf('fract_wood_to_pp_sawnwood', unitless, 'Fraction of harvested wood allocated to the sawnwood product pool'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = NONE, &
          & lrestart = .TRUE., &
          & lrestart_cont = .TRUE., &
          & initval_r = 0.0_wp)
      END IF ! IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
    END IF
  END SUBROUTINE Init_q_sh_memory

#endif
END MODULE mo_q_sh_memory_class
