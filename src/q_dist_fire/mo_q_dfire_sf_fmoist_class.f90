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
!>#### routine to calculate fuel moisture
!>
!> There are 2 implementations found in the code
!>  - Thonicke et al. (2010)
!>  - Oberhagemann et al. (2025)
!>
MODULE mo_q_dfire_sf_fmoist_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message, finish
  USE mo_jsb_math_constants,    ONLY: eps8
  USE mo_q_dfire_sf_fuel_class, ONLY: sf_calc_has_any_fuel
  USE mo_quincy_stats,        ONLY: stats_average_weight
  USE mo_jsb_lctlib_class,    ONLY: t_lctlib_element
  !< Total number of fuel moisture(omega_o and omega_nl), omega_o index, omega_nl index
  USE mo_q_dfire_sf_constants, ONLY: df_sf_omegas_nr, df_sf_omega_o_idx, df_sf_omega_nl_idx
  USE mo_q_dfire_sf_classes,   ONLY: t_sf_lfmc, t_sf_mf

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_sf_fuel_moisture_for_thonicke, calc_sf_fuel_moisture_for_oberhagemann, &
            calc_mflg_ober

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_fmoist_class'
CONTAINS
  ! ======================================================================================================= !
  !>
  !> Fuel moisture calculation:
  !> - omega_o: dead fuel for 1, 10 and 100 hours + live grass (accounted for 1 hour)
  !> - omega_n,l: 1 hour dead fuel + live grass
  !>
  FUNCTION calc_sf_fuel_moisture_for_thonicke(nc, nalphas, nfuels, alphas, nesterov_index, litter_moisture, w_oi, live_fuel_lvg) &
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
  END FUNCTION calc_sf_fuel_moisture_for_thonicke


  ! ======================================================================================================= !
  !>
  !> Calculate relative fuel moisture content (omega_o and omega_nl) for a single gridcell,
  !> following Thonicke et al. (2010) Eq. 6 with adaptations based on Pfeiffer et al. (2013).
  !>
  PURE FUNCTION calc_all_omegas(nalphas, nfuels, alphas, nesterov_index, &
                                litter_moisture, &
                                live_fuel_lvg, w_oi) RESULT(omegas)
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
    !< According to Thonicke et al. (2010) eq. B3, the omega_1hr and lvg_moisture must be weighted in order to compute omega_n,l

    !< Compute omega_1hr
    weighted_alphas = 0._wp
    omega_1hr = 0._wp
    alpha_1hr = 0._wp
    omega_lvg = 0._wp
    IF (has_dead_fuel_1hr) THEN
      alpha_1hr = alphas(ix_alpha_1hr)
      omega_1hr = calc_omega(alpha_1hr, nesterov_index)
    END IF

    !< Compute alpha livegrass (Pfeiffer et al. (2013), Eq. A6 and A7)
    alpha_lvg = 0._wp
    IF (has_live_fuel_lvg .AND. nesterov_index > 0._wp) THEN
      !< Thonicke et al. (2010) Eq. B2
      omega_lvg = max(0._wp, (10._wp/9._wp) * litter_moisture - (1._wp / 9._wp) )

      !< Derived from Eq. 6 from Thonicke et al. (2010)
      IF (omega_lvg > 0._wp) alpha_lvg = - (log(omega_lvg) / nesterov_index)
    END IF

    !< Calculate 1hr fuel moisture
    dead_fuel_1hr = w_oi(ix_woi_1hr)

    !< Mixed 1 hour fuel class and livegrass
    omegas(df_sf_omega_nl_idx) = stats_average_weight(nr_fine_fuels, &
                                                        [omega_lvg, omega_1hr], &
                                                        [live_fuel_lvg, dead_fuel_1hr] )

    !< Calculate omega_o (1, 10, 100 and 1000 hours dead fuel + live grass)
    !< According to Thonicke et al. (2010), Eq.6 the ALPHAS must be weighted in order to compute omega_o

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
  !> Calculate relative fuel moisture content (omega)
  !>
  !> Eq. 6 from Thonicke et al. (2010)
  !>
  ELEMENTAL FUNCTION calc_omega(alpha, nesterov_index) RESULT (omega)
    REAL(wp), INTENT(in) :: alpha          !< Drying parameter for fuel class [Celsius-2]
    REAL(wp), INTENT(in) :: nesterov_index !< Nesterov index [Celsius^2]
    REAL(wp)             :: omega          !< Relative fuel moisture [-]

    omega = exp(-alpha * nesterov_index)
  END FUNCTION calc_omega


  ! ======================================================================================================= !
  !>
  !> calc_sf_fuel_moisture_for_oberhagemann
  !> Calculate the fuel moisture as in Pfeiffer, Andrews PL and Oberhagemann
  !>
  !> Fuel moisture is calculated for dead fuel classes and live fuel classes separately:
  !>  - Dead fuel classes (1, 10, 100 and 1000 hours) and dry cured grass
  !>  - Live fuel classes (woody-not yet implemented- and live grass)
  !>
  SUBROUTINE calc_sf_fuel_moisture_for_oberhagemann(nc, ndd, nlitros, nlv, &
                                        alphas_dd, &
                                        lvg_moisture, &
                                        dead_fuel_litter, &
                                        dead_fuel_grass, &
                                        nesterov_index, &
                                        mask, &
                                        mf)
    USE mo_q_dfire_sf_constants, ONLY: ix_1hr => df_sf_ix_sigma_dd_1hr, & !< dead fuel moisture indices
                                        ix_10hr => df_sf_ix_sigma_dd_10hr, &
                                        ix_100hr => df_sf_ix_sigma_dd_100hr, &
                                        ix_1000hr => df_sf_ix_sigma_dd_1000hr, &
                                        ix_clvg => df_sf_ix_sigma_dd_clvg
    USE mo_q_dfire_sf_constants, ONLY: ix_lvg => df_sf_ix_sigma_lv_lvg, & ! live fuel moisture indices
                                        ix_wdy => df_sf_ix_sigma_lv_wdy

    INTEGER, INTENT(in)       :: nc
    INTEGER, INTENT(in)       :: ndd                     !< Dead fuel classes
    INTEGER, INTENT(in)       :: nlitros                 !< Dead litter fuel classes involved in fire spread -1 to 100 hours-
    INTEGER, INTENT(in)       :: nlv                     !< Live fuel classes
    REAL(wp), INTENT(in)      :: alphas_dd(4)            !< Drying parameter for dead fuels [cm2/cm3] (1, 10, 100 hours, cured livegrass)
    REAL(wp), INTENT(in)      :: lvg_moisture(nc)        !< live grass fuel moisture [%]
    REAL(wp), INTENT(in)      :: dead_fuel_litter(nc,nlitros) !< dead fuel [g DM m-2] (#nc, #dead litter ros fuel classes -1 to 100 hours-)
    REAL(wp), INTENT(in)      :: dead_fuel_grass(nc)     !< dead fuel grass [g DM m-2] (#nc)
    REAL(wp), INTENT(in)      :: nesterov_index(nc)      !< Nesterov index = Dryness index [Celsius^2]
    LOGICAL, INTENT(in)       :: mask(nc)
    TYPE(t_sf_mf), INTENT(inout) :: mf

    REAL(wp) :: dead_fuel(nc,4) !< [g DM m-2]
    REAL(wp) :: alphas_weighted_ros(nc) !<  [deg C-2]
    REAL(wp) :: alphas_weighted_dd(nc) !<  [deg C-2]
    REAL(wp) :: mf_10hr
    INTEGER :: ic

    dead_fuel(:,ix_1hr:ix_100hr) = dead_fuel_litter
    dead_fuel(:,4)               = dead_fuel_grass

    DO ic = 1, nc
      IF (.NOT. mask(ic)) CYCLE
      mf_10hr = 0._wp

      !< Ensure there is fuel
      IF (SUM(dead_fuel_litter(ic,:)) > eps8) THEN
        alphas_weighted_ros(ic) = stats_average_weight(nlitros, alphas_dd(ix_1hr:ix_100hr), dead_fuel_litter(ic,:))
        mf_10hr = calc_omega(alphas_weighted_ros(ic), nesterov_index(ic))
      ENDIF

      !< Ensure there is fuel
      IF (SUM(dead_fuel(ic,:)) > eps8) THEN
        alphas_weighted_dd(ic) = stats_average_weight(4, alphas_dd, dead_fuel(ic,:))
        mf%dead(ic,ix_1hr) = calc_omega(alphas_weighted_dd(ic), nesterov_index(ic)) !< Pfeiffer
      ENDIF

      mf%dead(ic,ix_10hr) = mf_10hr * 1.01_wp !< Andrews PL (2018), section 5.3.1, 2nd paragraph
      mf%dead(ic,ix_100hr) = mf_10hr * 1.02_wp
      mf%dead(ic,ix_1000hr) = mf_10hr * 1.02_wp
      mf%dead(ic,ix_clvg) = mf%dead(ic,ix_1hr) ! Same as in Oberhagemann et al. (2025), fuelload.c

      !< in case  mf_10hr > 1.0_wp, limit to 1.0_wp
      mf%dead(ic,ix_10hr) = MIN(mf%dead(ic,ix_10hr), 1.0_wp)
      mf%dead(ic,ix_100hr) = MIN(mf%dead(ic,ix_100hr), 1.0_wp)
      mf%dead(ic,ix_1000hr) = MIN(mf%dead(ic,ix_1000hr), 1.0_wp)
    ENDDO
    mf%live(:,ix_lvg) = lvg_moisture
    mf%live(:,ix_wdy) = 0._wp !  For later, for shrubs
  END SUBROUTINE calc_sf_fuel_moisture_for_oberhagemann


  ! ======================================================================================================= !
  !>
  !> calc_mflg_ober
  !> Calculate livegrass moisture based on Oberhagemann et al. (2025), eq. 4.
  !>
  FUNCTION calc_mflg_ober(nc, lfmc_par, phen) RESULT (mf_lg)
    INTEGER, INTENT(in) :: nc
    TYPE(t_sf_lfmc), INTENT(in) :: lfmc_par  !< livegrass moisture parameters
    REAL(wp), INTENT(in) :: phen(nc)      !< Livegrass phenological status [0 to 1]
    REAL(wp), ALLOCATABLE :: mf_lg(:)     ! Livegrass moisture [0.3 to 2.5]

    ! As paper's source code
    mf_lg = lfmc_par%b + (lfmc_par%max - lfmc_par%b) * phen
    WHERE (mf_lg < lfmc_par%min) mf_lg = lfmc_par%min !< set min value
    WHERE (mf_lg > lfmc_par%max) mf_lg = lfmc_par%max !< up to 250% livegrass moisture
  END FUNCTION calc_mflg_ober

#endif
END MODULE mo_q_dfire_sf_fmoist_class
