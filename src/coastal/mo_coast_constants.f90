! Author: Dr. Sri Nandini-Weiß
! Date: 2026-09-23
!> Constants of the coastal-constraint process
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
!>#### Fixed parameters of the salinity and tidal-inundation constraints
!>
!> Parameters that a user may wish to change without recompiling (the minima of the two
!> control factors, the ventilation thresholds and weights, the stand density) are
!> namelist parameters of mo_coast_config_class. The values below are fixed properties
!> of the published equations or of the tidal constituents.
!>
!> Source: Zhao et al. (2025), ORCHIDEE-MAN, J. Adv. Model. Earth Syst.,
!> doi:10.1029/2025MS005185, Equation 1 and Table 1; tidal speeds from the Doodson
!> numbers of the constituents as used by TPXO9-atlas (Egbert and Erofeeva, 2002).
!>
MODULE mo_coast_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  !> lctlib id of the mangrove land cover type [-]
  !> In the npft=13 usecase the tiles pft01 to pft13 are created with the lctlib ids
  !> 1,2,3,4,9,10,11,12,14,15,19,20,21 (mo_jsb_model_usecases), so the mangrove clone of
  !> tropical evergreen trees, lctlib column 22, has id 21.
  INTEGER,  PARAMETER :: mangrove_lct_id = 21

  !> Coefficients of the salinity control factor, Cs = cs_a2*S**2 + cs_a1*S + cs_a0
  !> [PSS-2, PSS-1, -]; Zhao et al. (2025), Equation 1, fitted to 26 observations of
  !> mangrove gross primary productivity against the 2000-2014 mean surface salinity of
  !> Bio-ORACLE v2.0.
  REAL(wp), PARAMETER :: cs_a2 = -0.0214_wp
  REAL(wp), PARAMETER :: cs_a1 =  1.28_wp
  REAL(wp), PARAMETER :: cs_a0 = -18.3_wp

  !> Maximum of the salinity quadratic, reached at S = -cs_a1/(2*cs_a2) = 29.907 PSS [-]
  !> The publication describes the maximum of the factor as one; the published
  !> coefficients give 0.8402. It is used by the normalisation switch l_normalise_cs.
  REAL(wp), PARAMETER :: cs_quadratic_maximum = cs_a0 - cs_a1 * cs_a1 / (4._wp * cs_a2)

  !> Number of tidal constituents, in the order M2, S2, N2, K2, K1, O1, P1, Q1 [-]
  INTEGER,  PARAMETER :: ntide = 8

  !> Angular speeds of the tidal constituents [degree h-1]
  REAL(wp), PARAMETER :: tide_speed(ntide) = [ &
    & 28.9841042_wp, 30.0000000_wp, 28.4397295_wp, 30.0821373_wp, &
    & 15.0410686_wp, 13.9430356_wp, 14.9589314_wp, 13.3986609_wp ]

  !> Maximum height of stilt roots [m]; Zhao et al. (2025), Table 1 (H_AGR_MAX_ST)
  REAL(wp), PARAMETER :: h_agrs_max = 2.0_wp

  !> Maximum height of pneumatophores [m]; Zhao et al. (2025), Table 1 (H_AGR_MAX_PN)
  REAL(wp), PARAMETER :: h_agrp_max = 0.4_wp

  !> Pneumatophore height per unit pneumatophore biomass [m (Mg ha-1)-1];
  !> Zhao et al. (2025), Equation 6, from a mean height of about 0.3 m at 4.55 Mg ha-1
  REAL(wp), PARAMETER :: h_agrp_per_biomass = 0.066_wp

#endif
END MODULE mo_coast_constants
