!> QUINCY calculate fuel moisture
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
!> @brief Module to calculate fuel moisture based on Thonicke et al. (2010)
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_q_dfire_sf_fmoist_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message
  USE mo_q_dfire_sf_fuel_class, ONLY: sf_calc_has_any_fuel
  USE mo_quincy_stats,        ONLY: stats_average_weight
  !< Total number of fuel moisture(omega_o and omega_nl), omega_o index, omega_nl index
  USE mo_q_dfire_constants,   ONLY: df_sf_omegas_nr, df_sf_omega_o_idx, df_sf_omega_nl_idx

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sf_calc_fuel_moisture

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_fmoist_class'

CONTAINS

  !>
  !> Fuel moisture calculation:
  !> - omega_o: dead fuel for 1, 10 and 100 hours + live grass (accounted for 1 hour)
  !> - omega_n,l: 1 hour dead fuel + live grass
  !>
  FUNCTION sf_calc_fuel_moisture(nc, nalphas, nfuels, alphas, nesterov_index, litter_moisture, w_oi, live_fuel_lvg) &
                                      RESULT (fuel_moisture)
    INTEGER, INTENT(in)                 :: nc, nalphas, nfuels        !< Number of gridcells, number of alpha params, number of fuels used to calc fuel moisture
    REAL(wp), INTENT(in)                :: alphas(nalphas)            !< Drying parameter for 1, 10 and 100 hours fuel class [Celsius-2]
    REAL(wp), INTENT(in)                :: nesterov_index(nc)         !< Nesterov index = Dryness index [Celsius^2]
    REAL(wp), INTENT(in)                :: litter_moisture(nc)        !< Relative litter moisture content, approximated by relative top-layer soil moisture content (== w_s,1) [0 to 1]
    REAL(wp), INTENT(in)                :: w_oi(nc,nfuels)            !< Dead fuel involved in the calculation of fuel moisture [g DM m-2]
    REAL(wp), INTENT(in)                :: live_fuel_lvg(nc)          !< Total live grass fuel [g DM m-2]
    REAL(wp)                            :: fuel_moisture(nc, df_sf_omegas_nr) !< output [0 to 1]

    LOGICAL :: has_live_fuel_lvg(nc), has_dead_fuel(nc)
    INTEGER :: inc
    REAL(wp) :: omegas(nc,df_sf_omegas_nr) !< 2 = omega_o, omega_n_l

    omegas = 0._wp
    has_dead_fuel = sf_calc_has_any_fuel(w_oi)
    has_live_fuel_lvg = (live_fuel_lvg > 0._wp)

    !< init outputs
    fuel_moisture = 0._wp

    DO inc = 1, nc
      IF (has_live_fuel_lvg(inc) .OR. has_dead_fuel(inc)) THEN
        omegas(inc,:) = calc_all_omegas(nalphas, nfuels, alphas, nesterov_index(inc), &
            litter_moisture(inc), live_fuel_lvg(inc), w_oi(inc,:))
      END IF
    END DO

    fuel_moisture = omegas
  END FUNCTION sf_calc_fuel_moisture


  ! ======================================================================================================= !
  !>
  !> Calculate relative fuel moisture content (omega_o and omega_nl) for a single gridcell,
  !> following Thonicke et al. (2010) Eq. 6 with adaptations based on Pfeiffer et al. (2013).
  !>
  PURE FUNCTION calc_all_omegas(nalphas, nfuels, alphas, nesterov_index, litter_moisture, live_fuel_lvg, w_oi) RESULT(omegas)
    INTEGER, INTENT(in)  :: nalphas, nfuels         !< Number of alpha params, number of fuel classes
    REAL(wp), INTENT(in) :: alphas(nalphas)         !< Drying parameter for 1, 10, 100 hours fuel classes [Celsius^2]
    REAL(wp), INTENT(in) :: nesterov_index          !< Nesterov index [Celsius^2]
    REAL(wp), INTENT(in) :: litter_moisture         !< Relative litter moisture content, approximated by relative top-layer soil moisture content (== w_s,1) [0 to 1]
    REAL(wp), INTENT(in) :: live_fuel_lvg           !< Livegrass fuel [g DM m-2]
    REAL(wp), INTENT(in) :: w_oi(nfuels)            !< Total dead fuel (1, 10 and 100 hours) [g DM m-2]
    REAL(wp)             :: omegas(df_sf_omegas_nr) !< output 1=omega_o, 2=omega_n_l

    INTEGER, PARAMETER :: ix_woi_1hr = 1
    INTEGER, PARAMETER :: nr_fuels = 4              !< 4 = (3 = total dead fuels (1,10,100)) + (1 = livegrass)
    INTEGER, PARAMETER :: ix_fuel_lvg = 1, ix_fuel_1hr = 2, ix_fuel_100hr = 4
    INTEGER, PARAMETER :: ix_alpha_1hr = 1
    INTEGER, PARAMETER :: nr_fine_fuels = 2         !< 1 hour dead fuel + livegrass

    REAL(wp) :: omega_lvg, omega_1hr

    REAL(wp) :: alpha_lvg, alpha_1hr
    REAL(wp) :: alphas_o(nr_fuels), alpha_o_weighted, weighted_alphas

    REAL(wp) :: dead_fuel_1hr
    REAL(wp) :: fuels_o(nr_fuels) !< fuels used to calculate omega_o (dead 1, 10 and 100 hours + live grass)
    LOGICAL :: has_live_fuel_lvg, has_dead_fuel, has_dead_fuel_1hr

    has_dead_fuel = (SUM(w_oi) > 0._wp)
    has_dead_fuel_1hr = (w_oi(ix_fuel_1hr) > 0._wp)
    has_live_fuel_lvg = (live_fuel_lvg > 0._wp)

    !< init outputs
    omegas(df_sf_omega_o_idx) = 0._wp   !< omega_o: combine all fuels (1,10,100 and 1000 hours dead fuel + live grass)
    omegas(df_sf_omega_nl_idx) = 0._wp  !< omega_n,l: combine 1 hour dead fuel and live grass

    !< Calculate omega_n,l (dead 1hr + livegrass)
    !< According to Thonicke et al. (2010) eq. B3, the omega_1hr and omega_lvg must be weighted in order to compute omega_n,l

    !< Compute omega_1hr
    weighted_alphas = 0._wp
    omega_1hr = 0._wp
    alpha_1hr = 0._wp
    IF (has_dead_fuel_1hr) THEN
      alpha_1hr = alphas(ix_alpha_1hr)
      omega_1hr = calc_omega(alpha_1hr, nesterov_index)
    END IF

    !< Compute alpha livegrass (Pfeiffer et al. (2013), Eq. A6 and A7)
    omega_lvg = 0._wp
    alpha_lvg = 0._wp
    IF (has_live_fuel_lvg .AND. nesterov_index > 0._wp) THEN
      !< Thonicke et al. (2010) Eq. B2
      omega_lvg = max(0._wp, (10._wp/9._wp) * litter_moisture - (1._wp / 9._wp) )
      !< Derived from Eq. 6 from Thonicke et al. (2010)
      alpha_lvg = - (log(omega_lvg) / nesterov_index)
    END IF

    !< Calculate 1hr fuel moisture
    dead_fuel_1hr = w_oi(ix_woi_1hr)

    !< Mixed 1 hour fuel class and livegrass
    omegas(df_sf_omega_nl_idx) = stats_average_weight(nr_fine_fuels, &
                                                        [omega_lvg, omega_1hr], &
                                                        [live_fuel_lvg, dead_fuel_1hr] )

    !< Calculate omega_o (1, 10, 100 and 1000 hours dead fuel + live grass)
    !< According to Thonicke et al. (2010), Eq.6 the ALPHAS must be weighted in order to compule omega_o

    !< Incorporate livegrass fuel (1 hour)
    fuels_o(ix_fuel_lvg) = live_fuel_lvg
    fuels_o(ix_fuel_1hr:ix_fuel_100hr) = w_oi

    alphas_o(ix_fuel_lvg) = alpha_lvg
    alphas_o(ix_fuel_1hr:ix_fuel_100hr) = alphas

    alpha_o_weighted = stats_average_weight(nr_fuels, alphas_o, fuels_o)
    omegas(df_sf_omega_o_idx) = calc_omega(alpha_o_weighted, nesterov_index)
  END FUNCTION calc_all_omegas


  ! ======================================================================================================= !
  !>
  !> calc_omega
  !> Calculate the fuel moisture
  !>
  !> Eq. 6 from Thonicke et al. (2010)
  !>
  ELEMENTAL FUNCTION calc_omega(alpha, nesterov_index) RESULT (omega)
    REAL(wp), INTENT(in) :: alpha          !< Drying parameter for fuel class [Celsius-2]
    REAL(wp), INTENT(in) :: nesterov_index !< Nesterov index [Celsius^2]
    REAL(wp)             :: omega          !< Relative fuel moisture [-]

    omega = exp(-alpha * nesterov_index)
  END FUNCTION calc_omega

#endif
END MODULE mo_q_dfire_sf_fmoist_class
