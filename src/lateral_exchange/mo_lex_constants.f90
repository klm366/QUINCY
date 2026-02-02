!> Contains constants used in the lateral exchange processes
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
!>#### Connection matrixes used in the lateral exchange process
!>
!> This module contains base connection matrixes for the lateral exchange processes between tiles.
!>
!> Note! The 1-D arrays defined below are being reshaped to ntiles x ntiles matrixes
!> in the lex_init routine: M = reshape(1-D_array, ntiles x ntiles, order=(/2,1/))
!
!> Note! Not all connectivities are directly used in the model. Some are adjusted to
!> be consistent with the topographic information provided.
!>
!> Base connection matrixes:
!> micro-vs. macro-scale: The separation between micro- and macro scale is somewhat ambiguous.
!>                    A basic feature of micro-scale connectivities is a comparativly short distance
!>                    from the center of a tile (or rather of a patch/cluster represented by a tile)
!>                    to the edge of the tile. For lack of a better definition, "short" means that
!>                    the distance is so small that diffusive and conductive fluxes become relevant.
!>                    Examples of micro-scale connections are small-scale depressions / ponds,
!>                    patches of organic soils or groves. Whereas macro-scale connections could refer
!>                    to those between patches located on uplands and slopes to those in the lowlands.
!>
!> cmx_micro_ctc:     Connectivity based on direct contact between two interfacing tiles.
!>                    The matrix elements constitute the relative (geometrical) contact lengths,
!>                    i.e. the fraction of the edge of a tile that interfaces with a given sibling.
!>                    The entries need to be given for the tile (i.e. on the row)
!>                    which is partially encompassed by the sibling (indicated by the column).
!>                    If this can not be specified, then the entries should be given for the tile
!>                    with the smaller characteristic surface area.
!>
!> cmx_micro_hfp_srf: Connectivity for interactions via hydrological flow paths at the surface.
!>                    The matrix elements give the preferential fraction of fluxes originating
!>                    on the tile given by the row to the tile indicated by the column.
!>                    Note that the preferential fraction is not being used directly,
!>                    but only in combination with the cover-fraction of the receiving tile
!>                    relative to the sum of the cover fractions of all receiving tiles.
!>
!> cmx_micro_hfp_blg: Same as cmx_micro_hfp_srf but for below-ground hydrological flow paths
!> cmx_macro_hfp_srf: Same as cmx_micro_hfp_srf but for macro-scale interactions
!> cmx_macro_hfp_blg: Same as cmx_micro_hfp_blg but for macro-scale interactions
!>
!> For cmx_micro_ctc the sum of entries on a row may not exeed 1.
!> Furthermore, the model does not (yet) properly consider nested configurations.
!> The possibilities for such setups are limited to a case in which a single surface feature SF1
!> is fully encompassed by a second surface feature SF2, which, in turn,
!> is encompassed by another (other) surface feature(s).
!> If the model is set up for a cascade of connected tiles (l_nested = FALSE),
!> the geometric relationships can only be determined for each pair of tiles individually.
!>
!> For cmx_micro_hfp_x and cmx_macro_hfp_x, indiviual matrix elements should not exceed 1.
!> However, it is intended that the sum of the elements on a row can be larger than 1.
!> The flow lags assumed for surface runoff and drainage are linearly dependend on the matrix elements
!> of cmx_micro_hfp_x and cmx_macro_hfp_x. Hence, the respective matrix entries determine the overall
!> retention time hence the outgoing fluxes. But the distribution amongst all receiving tiles
!> is done relative to the weighted sum of all entries on a row.
!> Consequently, the distribution only depends on the values relative to each other.
!> If, however, the sum of all entries on a row is smaller than 1, a fraction (1 - SUM(row elements))
!> of the outgoing fluxes is directly given to the streamflow.
!>
MODULE mo_lex_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  !--------------------------------------------------------------------------------------------------------
  ! Connection matrixes for lateral exchange between two tiles
  !
  REAL(wp), PARAMETER, DIMENSION(4) ::  &
    cmx_micro_ctc_2tilelex =            &
      & (/ 0.0, 0.0,                    &
      &    1.0, 0.0 /)

  REAL(wp), PARAMETER, DIMENSION(4) ::  &
    cmx_micro_hfp_srf_2tilelex =        &
      & (/ 0.0, 1.0,                    &
      &    0.0, 0.0 /)

  REAL(wp), PARAMETER, DIMENSION(4) ::  &
    cmx_micro_hfp_blg_2tilelex =        &
      & (/ 0.0, 1.0,                    &
      &    0.0, 0.0 /)

  REAL(wp), PARAMETER, DIMENSION(4) ::  &
    cmx_macro_hfp_srf_2tilelex =        &
      & (/ 0.0, 0.0,                    &
      &    0.0 ,0.0 /)

  REAL(wp), PARAMETER, DIMENSION(4) ::  &
    cmx_macro_hfp_blg_2tilelex =        &
      & (/ 0.0, 0.0,                    &
      &    0.0 ,0.0 /)

  !--------------------------------------------------------------------------------------------------------
  ! Connection matrixes for lateral exchange between three tiles
  !
  REAL(wp), PARAMETER, DIMENSION(9) ::  &
    cmx_micro_ctc_3tilelex =            &
      & (/ 0.0, 1.0, 0.0,               &
      &    0.0, 0.0, 1.0,               &
      &    0.0, 0.0, 0.0 /)

  REAL(wp), PARAMETER, DIMENSION(9) ::  &
    cmx_micro_hfp_srf_3tilelex =        &
      & (/ 0.0, 0.0, 0.0,               &
      &    0.5, 0.0, 0.5,               &
      &    0.0, 0.0, 0.0 /)

  REAL(wp), PARAMETER, DIMENSION(9) ::  &
    cmx_micro_hfp_blg_3tilelex =        &
      & (/ 0.0, 0.0, 0.0,               &
      &    0.5, 0.0, 0.5,               &
      &    0.0, 0.0, 0.0 /)

  REAL(wp), PARAMETER, DIMENSION(9) ::  &
    cmx_macro_hfp_srf_3tilelex =        &
      & (/ 0.0, 0.0, 0.0,               &
      &    0.0, 0.0, 0.0,               &
      &    0.0, 0.0, 0.0 /)

  REAL(wp), PARAMETER, DIMENSION(9) ::  &
    cmx_macro_hfp_blg_3tilelex =        &
      & (/ 0.0, 0.0, 0.0,               &
      &    0.0, 0.0, 0.0,               &
      &    0.0, 0.0, 0.0 /)

  !--------------------------------------------------------------------------------------------------------
  ! Connection matrixes for lateral exchange between four tiles (called topolex setup)
  !
  REAL(wp), PARAMETER, DIMENSION(16) :: &
    cmx_micro_ctc_topolex =             &
      & (/ 0.0, 0.0, 0.0 ,0.0,          &
      &    1.0, 0.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 1.0 ,0.0 /)

  REAL(wp), PARAMETER, DIMENSION(16) :: &
    cmx_micro_hfp_srf_topolex =         &
      & (/ 0.0, 1.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,1.0,          &
      &    0.0, 0.0, 0.0 ,0.0 /)

  REAL(wp), PARAMETER, DIMENSION(16) :: &
    cmx_micro_hfp_blg_topolex =         &
      & (/ 0.0, 1.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,1.0,          &
      &    0.0, 0.0, 0.0 ,0.0 /)

  REAL(wp), PARAMETER, DIMENSION(16) :: &
    cmx_macro_hfp_srf_topolex =         &
      & (/ 0.0, 0.0, 0.0 ,1.0,          &
      &    0.0, 0.0, 0.0 ,1.0,          &
      &    0.0, 0.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,0.0 /)

  REAL(wp), PARAMETER, DIMENSION(16) :: &
    cmx_macro_hfp_blg_topolex =         &
      & (/ 0.0, 0.0, 0.5 ,0.5,          &
      &    0.0, 0.0, 0.5 ,0.5,          &
      &    0.0, 0.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,0.0 /)

#endif
END MODULE mo_lex_constants