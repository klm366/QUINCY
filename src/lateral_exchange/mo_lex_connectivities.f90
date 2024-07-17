!> Contains base connection matrixes for the later exchange processes between tiles
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
MODULE mo_lex_connectivities
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  ! Base connection matrixes:
  !
  ! Note! The 1-D arrays that need to be provided below, are being reshaped to ntiles x ntiles matrixes in the lex_init routine!
  !
  ! In case of cmx_micro_hfp and cmx_macro_hfp, the base connection matrixes give the fraction of fluxes (runoff and drainage) 
  ! from a given tile to a specific, connected (downstream) tile assuming (preffered) hydrological flow paths.
  ! In case of cmx_micro_ctc they indicate the direct connections, i.e. the fraction of the outer edge of a tile 
  ! that is assumed to be connected to a specific sibling tile. 
  ! Note that connectivties are not directly used in the model, but adjusted based on their consistance
  ! with the topographic information provided. 
  !
  ! cmx_micro_ctc: Connectivity based on direct contact between two connected tiles. 
  !                Assumed relevant only for the micro-scale. Here, micro scale connections reflect landscapes where
  !                the tiles are well mixed and have comparativly small characteristic surface areas. 
  !                such as ponds and small scale depressions. 
  !                The matrix elements need to be provided so that entries are given for the tile which is encompassed 
  !                by the other tile, if this can not be specified, then the intries should be given for the tile with the 
  !                smaller characteristic surface area (given by the row)  
  !                which is assumed to be partially encompassed by a larger tile (indicated by the column)
  !
  ! cmx_micro_hfp: Connectivity for connections via hydrological flow paths (surface runoff & lateral drainage)
  !                between two connected tiles. Here, micro scale connections reflect landscapes where
  !                the tiles are well mixed and have comparativly small characteristic surface areas. 
  !                such as ponds and small scale depressions. 
  !                The matrix elements need to be provided to indicate the fractions of fluxes originating
  !                on the tile given by the row to the tile indicated by the column.
  !
  ! cmx_macro_hfp: Connectivity for connections via hydrological flow paths (surface runoff & lateral drainage)
  !                between two tiles that are spatially seperated across larger distances, such as uplands, slopes and lowlands.
  !                The matrix elements need to be provided to so that entries are given for the tile with the 
  !                The matrix elements need to be provided to indicate the fractions of fluxes originating
  !                on the tile given by the row to the tile indicated by the column.
  !
  ! For cmx_micro_ctc the sum of entries on a row may not exeed 1. Furthermore, the model does not (yet) properly considder nested configurations.
  ! These possibility for such setups are limited to a case in which a single surface feature SF1 is fully encompassed by a second surface feature SF2, 
  ! which ,in turn, is encompassed by another (other) surface feature(s). If the model is set up for a cascade of connected tiles (l_nested = FALSE),
  ! the geometric relashionships can only be determined for each pair of tiles individually.
  ! For cmx_micro_hfp and cmx_macro_hfp, indiviual matrix entries should not exceed 1. 
  ! However, it is intended that the sum of a row can be larger than 1. The flow lags assumed for surface runoof and drainage
  ! are linearly dependend on the matrix elements of cmx_micro_hfp and cmx_macro_hfp. Hence, the respective matrix entries 
  ! determine the overall retention time hence the outgoing fluxes. But the distribution amongst all receiving tiles
  ! is done relative to the weighted sum of all entries on a row, so the latter only depends on the values relatve to each other.
  ! If, however, the sum of all entries on a row is smaller than 1, a fraction (1 - SUM(row elements)) of the outgoing fluxes
  ! is directily given to the streamflow.


! -----------------------------------------------------------------------------------
! 2tilelex

  REAL(wp), PARAMETER, DIMENSION(4) ::  &
    cmx_micro_ctc_2tilelex =           &
      & (/ 0.0, 0.0,                    &
      &    1.0, 0.0 /)

  REAL(wp), PARAMETER, DIMENSION(4) ::  &
    cmx_micro_hfp_2tilelex =           &
      & (/ 0.0, 1.0,                    &
      &    0.0, 0.0 /)

  REAL(wp), PARAMETER, DIMENSION(4) ::  &
    cmx_macro_hfp_2tilelex =           &
      & (/ 0.0, 0.0,                    &
      &    0.0 ,0.0 /)

! -----------------------------------------------------------------------------------
! topolex

  REAL(wp), PARAMETER, DIMENSION(16) :: &
    cmx_micro_ctc_topolex =            &
      & (/ 0.0, 0.0, 0.0 ,0.0,          &
      &    1.0, 0.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 1.0 ,0.0 /)

  REAL(wp), PARAMETER, DIMENSION(16) :: &
    cmx_micro_hfp_topolex =            &
      & (/ 0.0, 1.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,1.0,          &
      &    0.0, 0.0, 0.0 ,0.0 /)

  REAL(wp), PARAMETER, DIMENSION(16) :: &
    cmx_macro_hfp_topolex =            &
      & (/ 0.0, 0.0, 0.0 ,1.0,          &
      &    0.0, 0.0, 0.0 ,1.0,          &
      &    0.0, 0.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,0.0 /)

! -----------------------------------------------------------------------------------

#endif
END MODULE mo_lex_connectivities
