!> Contains the memory class for the land ice (lice) process
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
!>#### Memory of the land ice (lice) process
!>
!> The variables are defined on the glacier tile (where the process runs) and on all
!> tiles that contain the glacier land cover type (land, box), where they are filled
!> by aggregation.
!>
MODULE mo_lice_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: lp
  USE mo_util,             ONLY: One_of
  USE mo_jsb_lct_class,    ONLY: GLACIER_TYPE
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d
  USE mo_jsb_varlist,      ONLY: BASIC

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lice_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 10

  TYPE, EXTENDS(t_jsb_memory) :: t_lice_memory

    TYPE(t_jsb_var_real2d) :: &
      & icefract,             & !< Ice-covered fraction of the grid cell received from the external
                                !< land ice model (EBFM); -1: no data from external model for this cell
      & icefract_valid,       & !< 1 where icefract has been received (>= 0), 0 elsewhere
      & glac_fract_diff         !< icefract minus glacier fraction of the grid cell in ICON-Land
                                !< (0 where no data from the external model)

  CONTAINS
    PROCEDURE :: Init => Init_lice_memory
  END TYPE t_lice_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lice_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !> Initialize memory (variables) for the lice process
  !>
  SUBROUTINE Init_lice_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_lice_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),     INTENT(in)    :: prefix
    CHARACTER(len=*),     INTENT(in)    :: suffix
    INTEGER,              INTENT(in)    :: lct_ids(:)
    INTEGER,              INTENT(in)    :: lib_id     !< id of primary lct in lctlib
    INTEGER,              INTENT(in)    :: model_id

    TYPE(t_jsb_grid),  POINTER :: hgrid     ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface   ! Vertical grid

    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_lice_memory'

    IF (lib_id > 0)   CONTINUE ! avoid compiler warning about dummy argument not being used
    IF (model_id > 0) CONTINUE ! avoid compiler warning about dummy argument not being used

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    ! Variables exist on the glacier tile and on all tiles containing the glacier lct (aggregation)
    IF (One_of(GLACIER_TYPE, lct_ids(:)) > 0) THEN

      CALL mem%Add_var('icefract', mem%icefract,                                                 &
        & hgrid, surface,                                                                       &
        & t_cf('icefract', '1', 'Ice-covered fraction of grid cell from external land ice model (-1: no data)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
        & prefix, suffix,                                                                       &
        & output_level=BASIC,                                                                   &
        & lrestart=.TRUE.,                                                                      &
        & initval_r=-1.0_lp )

      CALL mem%Add_var('icefract_valid', mem%icefract_valid,                                     &
        & hgrid, surface,                                                                       &
        & t_cf('icefract_valid', '1', 'Flag: ice fraction received from external land ice model'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
        & prefix, suffix,                                                                       &
        & output_level=BASIC,                                                                   &
        & lrestart=.FALSE.,                                                                     &
        & initval_r=0.0_lp )

      CALL mem%Add_var('glac_fract_diff', mem%glac_fract_diff,                                   &
        & hgrid, surface,                                                                       &
        & t_cf('glac_fract_diff', '1', 'Ice fraction from external land ice model minus ICON-Land glacier fraction'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
        & prefix, suffix,                                                                       &
        & output_level=BASIC,                                                                   &
        & lrestart=.FALSE.,                                                                     &
        & initval_r=0.0_lp )

    END IF

  END SUBROUTINE Init_lice_memory

#endif
END MODULE mo_lice_memory_class
