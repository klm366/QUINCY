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
!> The variables are defined on the lice tile (the glacier tile coupled to EBFM, where the
!> process runs) and on its parent tiles (land, box), where they are filled
!> by aggregation. Where the external model provides no data, the state variables carry the
!> CDI missing value (mo_jsb_io::missval) and the mass fluxes/amounts are zero.
!>
MODULE mo_lice_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: lp
  USE mo_util,             ONLY: One_of
  USE mo_jsb_lct_class,    ONLY: GLACIER_TYPE
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d
  USE mo_jsb_varlist,      ONLY: BASIC
  USE mo_jsb_io,           ONLY: missval

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lice_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 20

  TYPE, EXTENDS(t_jsb_memory) :: t_lice_memory

    TYPE(t_jsb_var_real2d) :: &
      & icefract,             & !< Ice-covered fraction of the grid cell received from the external
                                !< land ice model (EBFM); missing value where the external model has no data
      & icefract_valid,       & !< 1 where icefract has been received (>= 0), 0 elsewhere
      & glac_fract_diff,      & !< icefract minus fraction of the lice tile (ICON-Land's EBFM-coupled ice)
                                !< (missing value where no data from the external model)
      & albedo,               & !< Surface albedo received from the external land ice model (EBFM);
                                !< missing value where the external model has no data
      & t_sub,                & !< Temperature of the first subsurface firn layer of EBFM [K]; missing value where no data
      & ghf_cond,             & !< Conductance between surface and first firn layer of EBFM [W m-2 K-1]; missing value where no data
      & hcap_sub,             & !< Heat capacity of the firn surface layer of EBFM [J m-2 K-1]; missing value where no data
      & runoff,               & !< Runoff from the EBFM firn column [kg m-2 s-1] (per unit ice area; 0 where no data)
      & smb,                  & !< Climatic surface mass balance of EBFM [kg m-2 s-1] (0 where no data)
      & snowmass                !< Snow mass on top of the ice in EBFM [kg m-2] (0 where no data)

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

    ! Variables exist on tiles with the glacier lct where the process is calculated or aggregated
    ! (lice tile and its parents; the sibling glacier tile 'glac' has no lice memory)
    IF (One_of(GLACIER_TYPE, lct_ids(:)) > 0) THEN

      CALL mem%Add_var('icefract', mem%icefract,                                                 &
        & hgrid, surface,                                                                       &
        & t_cf('icefract', '1', 'Ice-covered fraction of grid cell from external land ice model'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
        & prefix, suffix,                                                                       &
        & output_level=BASIC,                                                                   &
        & lrestart=.TRUE.,                                                                      &
        & initval_r=missval )

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
        & t_cf('glac_fract_diff', '1', 'Ice fraction from external land ice model minus lice tile fraction'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
        & prefix, suffix,                                                                       &
        & output_level=BASIC,                                                                   &
        & lrestart=.FALSE.,                                                                     &
        & initval_r=missval )

      CALL mem%Add_var('albedo', mem%albedo,                                                     &
        & hgrid, surface,                                                                       &
        & t_cf('albedo', '1', 'Surface albedo from external land ice model'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
        & prefix, suffix,                                                                       &
        & output_level=BASIC,                                                                   &
        & lrestart=.TRUE.,                                                                      &
        & initval_r=missval )

      CALL mem%Add_var('t_sub', mem%t_sub,                                                       &
        & hgrid, surface,                                                                       &
        & t_cf('t_sub', 'K', 'Temperature of first subsurface firn layer from external land ice model'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
        & prefix, suffix,                                                                       &
        & output_level=BASIC,                                                                   &
        & lrestart=.TRUE.,                                                                      &
        & initval_r=missval )

      CALL mem%Add_var('ghf_cond', mem%ghf_cond,                                                 &
        & hgrid, surface,                                                                       &
        & t_cf('ghf_cond', 'W m-2 K-1', 'Conductance surface - first firn layer from external land ice model'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
        & prefix, suffix,                                                                       &
        & output_level=BASIC,                                                                   &
        & lrestart=.TRUE.,                                                                      &
        & initval_r=missval )

      CALL mem%Add_var('hcap_sub', mem%hcap_sub,                                                 &
        & hgrid, surface,                                                                       &
        & t_cf('hcap_sub', 'J m-2 K-1', 'Heat capacity of the firn surface layer from external land ice model'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
        & prefix, suffix,                                                                       &
        & output_level=BASIC,                                                                   &
        & lrestart=.TRUE.,                                                                      &
        & initval_r=missval )

      CALL mem%Add_var('runoff', mem%runoff,                                                     &
        & hgrid, surface,                                                                       &
        & t_cf('runoff', 'kg m-2 s-1', 'Runoff from external land ice model (per unit ice area)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
        & prefix, suffix,                                                                       &
        & output_level=BASIC,                                                                   &
        & lrestart=.TRUE.,                                                                      &
        & initval_r=0.0_lp )

      CALL mem%Add_var('smb', mem%smb,                                                           &
        & hgrid, surface,                                                                       &
        & t_cf('smb', 'kg m-2 s-1', 'Climatic surface mass balance from external land ice model'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
        & prefix, suffix,                                                                       &
        & output_level=BASIC,                                                                   &
        & lrestart=.FALSE.,                                                                     &
        & initval_r=0.0_lp )

      CALL mem%Add_var('snowmass', mem%snowmass,                                                 &
        & hgrid, surface,                                                                       &
        & t_cf('snowmass', 'kg m-2', 'Snow mass on top of the ice from external land ice model'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
        & prefix, suffix,                                                                       &
        & output_level=BASIC,                                                                   &
        & lrestart=.FALSE.,                                                                     &
        & initval_r=0.0_lp )

    END IF

  END SUBROUTINE Init_lice_memory

#endif
END MODULE mo_lice_memory_class
