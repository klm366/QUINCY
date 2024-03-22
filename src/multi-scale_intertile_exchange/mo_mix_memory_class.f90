!> Contains the memory class for the multi-scale intertile exchange process.
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
MODULE mo_mix_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: !X e.g: LAND_TYPE, LAKE_TYPE, VEG_TYPE, BARE_TYPE, GLACIER_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist,            ONLY: BASIC, MEDIUM, FULL
  ! USE mo_jsb_physical_constants, ONLY: !X e.g: tmelt

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_mix_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 10

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_mix_memory

    ! Common variables independent of lct_type
    TYPE(t_jsb_var_real2d) ::   &
      & wtr_latflow,            & !< Lateral water flow between tiles               [kg m-2 s-1]
      & runoff_retention          !< Runoff retention reservoir                     [m water equivalent]
    TYPE(t_jsb_var_real3d) ::   &
      & wtr_diflat_sl,          & !< Lateral diffusion flux between tiles           [kg m-2 s-1]
      & drainage_retention_sl     !< Drainage retention reservoir                   [m water equivalent]

  CONTAINS
    PROCEDURE :: Init => Init_mix_memory
  END TYPE t_mix_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_mix_memory_class'

CONTAINS

  SUBROUTINE Init_mix_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_mix_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, soil_w              ! Vertical grids

    INTEGER :: table
    TYPE(t_grib2) :: grib2_desc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_mix_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    soil_w  => Get_vgrid('soil_depth_water')

    IF (TRIM(mem%owner_tile_name) == 'box') THEN
      grib2_desc = t_grib2(1,0,202, grib_bits)
    ELSE
      grib2_desc = t_grib2(255, 255, 255, grib_bits)
    END IF


    ! Common surface variables

    CALL mem%Add_var('wtr_latflow', mem%wtr_latflow,                                 &
      & hgrid, surface,                                                              &
      & t_cf('wtr_latflow', 'kg m-2 s-1', 'Lateral water flow between tiles'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix, l_aggregate_all=.TRUE.,                                      &
      & output_level=BASIC,                                                          &
      & lrestart=.FALSE., initval_r=0.0_wp )

    CALL mem%Add_var( 'runoff_retention', mem%runoff_retention,                      &
      & hgrid, surface,                                                              &
      & t_cf('runoff_retention', 'm (water equivalent)',                             &
      &      'Lateral overland flow retention reservoir'),                           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix, l_aggregate_all=.TRUE.,                                      &
      & output_level=BASIC,                                                          &
      & lrestart=.TRUE., initval_r=0.0_wp)

    ! Common soil layer variables

    CALL mem%Add_var( 'wtr_diflat_sl', mem%wtr_diflat_sl,                            &
      & hgrid, soil_w,                                                               &
      & t_cf('wtr_diflat_sl', 'kg m-2 s-1',                                          &
      &      'Diffusion flux from neighbouring tile'),                               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=BASIC,                                                          &
      & lrestart=.FALSE., initval_r=0.0_wp)

    CALL mem%Add_var( 'drainage_retention_sl', mem%drainage_retention_sl,            &
      & hgrid, soil_w,                                                               &
      & t_cf('drainage_retention_sl', 'm (water equivalent)',                        &
      &      'Retention reservoir for drainage within soil layers'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix, l_aggregate_all=.TRUE.,                                      &
      & output_level=BASIC,                                                          &
      & lrestart=.TRUE., initval_r=0.0_wp)


  END SUBROUTINE Init_mix_memory

#endif
END MODULE mo_mix_memory_class
