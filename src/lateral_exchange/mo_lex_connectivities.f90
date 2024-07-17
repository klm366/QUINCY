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
  ! In case of micro_cnct_grvt and macro_cnct_grvt, the base connection matrixes give the fraction of fluxes (runoff and drainage) 
  ! from a given tile to a specific, connected (downstream) tile. In case of micro_cnct_dfsn they indicate
  ! the fraction of the outer border of a tile that is assumed to be connected to a specific tile. 
  ! Note that connectivties are not directly used in the model, but adjusted based on their consitance
  ! with the topographic information provided. 
  !
  ! micro_cnct_dfsn: Connectivity for diffusive fluxes between two interfacing tiles. 
  !                  Here, micro scale connections reflect landscapes where
  !                  the tiles are well mixed and have comparativly small characteristic surface areas. 
  !                  such as ponds and small scale depressions. 
  !                  The matrix elements need to be provided so that entries are given for the tile which is encompassed 
  !                  by the other tile, if this can not be specified, then the intries should be given for the tile with the 
  !                  smaller characteristic surface area (given by the row)  
  !                  which is assumed to be partially encompassed by a larger tile (indicated by the column)
  !
  ! micro_cnct_grvt: Connectivity for gravity driven fluxes (surface runoff & lateral drainage)
  !                  between two interfacing tiles. Here, micro scale connections reflect landscapes where
  !                  the tiles are well mixed and have comparativly small characteristic surface areas. 
  !                  such as ponds and small scale depressions. 
  !                  The matrix elements need to be provided to indicate the fractions of fluxes originating
  !                  on the tile given by the row to the tile indicated by the column.
  !
  !
  ! macro_cnct_grvt: Connectivity for gravity driven fluxes (surface runoff & lateral drainage)
  !                  between two tiles that are spatially seperated across larger distances, such as uplands, slopes and lowlands.
  !                  The matrix elements need to be provided to so that entries are given for the tile with the 
  !                  The matrix elements need to be provided to indicate the fractions of fluxes originating
  !                  on the tile given by the row to the tile indicated by the column.
  !
  ! For micro_cnct_dfsn the sum of entries on a row may not exeed 1. Furthermore, the model does not (yet) properly considder nested configurations.
  ! These possibility for such setups are limited to a case in which a single surface feature SF1 is fully encompassed by a second surface feature SF2, 
  ! which ,in turn, is encompassed by another (other) surface feature(s). If the model is set up for a cascade of connected tiles (l_nested = FALSE),
  ! the geometric relashionships can only be determined for each pair of tiles individually.
  ! For micro_cnct_grvt and macro_cnct_grvt, indiviual matrix entries should not exceed 1. 
  ! However, it is intended that the sum of a row can be larger than 1. The flow lags assumed for surface runoof and drainage
  ! are linearly dependend on the matrix elements of micro_cnct_grvt and macro_cnct_grvt. Hence, the respective matrix entries 
  ! determine the overall retention time hence the outgoing fluxes. But the distribution amongst all receiving tiles
  ! is done relative to the weighted sum of all entries on a row, so the latter only depends on the values relatve to each other.
  ! If, however, the sum of all entries on a row is smaller than 1, a fraction (1 - SUM(row elements)) of the outgoing fluxes
  ! is directily given to the streamflow.


! -----------------------------------------------------------------------------------
! 2tilelex

  REAL(wp), PARAMETER, DIMENSION(4) ::  &
    micro_cnct_dfsn_2tilelex =          &
      & (/ 0.0, 0.0,                    &
      &    1.0, 0.0 /)

  REAL(wp), PARAMETER, DIMENSION(4) ::  &
    micro_cnct_grvt_2tilelex =          &
      & (/ 0.0, 1.0,                    &
      &    0.0, 0.0 /)

  REAL(wp), PARAMETER, DIMENSION(4) ::  &
    macro_cnct_grvt_2tilelex =          &
      & (/ 0.0, 0.0,                    &
      &    0.0 ,0.0 /)

! -----------------------------------------------------------------------------------
! topolex

  REAL(wp), PARAMETER, DIMENSION(16) :: &
    micro_cnct_dfsn_topolex =           &
      & (/ 0.0, 0.0, 0.0 ,0.0,          &
      &    1.0, 0.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 1.0 ,0.0 /)

  REAL(wp), PARAMETER, DIMENSION(16) :: &
    micro_cnct_grvt_topolex =           &
      & (/ 0.0, 1.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,1.0,          &
      &    0.0, 0.0, 0.0 ,0.0 /)

  REAL(wp), PARAMETER, DIMENSION(16) :: &
    macro_cnct_grvt_topolex =           &
      & (/ 0.0, 0.0, 0.0 ,1.0,          &
      &    0.0, 0.0, 0.0 ,1.0,          &
      &    0.0, 0.0, 0.0 ,0.0,          &
      &    0.0, 0.0, 0.0 ,0.0 /)

! -----------------------------------------------------------------------------------

#endif
END MODULE mo_lex_connectivities
