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

  !< GSI parameters, Forkel et al. (2014), Figure 2
  !< "slope": slope of X limiting function tmax_sl */
  !< "base" : inflection point of X limiting function (deg C) */
  !< "tau"  : change rate of actual to previous day X limiting fct */
  REAL(wp), PARAMETER :: gsi_heat(3) = [41.51_wp, 1.74_wp, 0.2_wp] ! base, slope, tau
  REAL(wp), PARAMETER :: gsi_cold(3) = [7._wp, 0.2591_wp, 0.2_wp] ! base, slope, tau
  REAL(wp), PARAMETER :: gsi_light(3) = [59.78_wp, 58._wp, 0.2_wp] ! base, slope, tau
  REAL(wp), PARAMETER :: gsi_wtr(3) = [20.96_wp, 5.24_wp, 0.2_wp] ! base, slope, tau
  INTEGER, PARAMETER :: ix_base = 1, ix_sl = 2, ix_tau=3
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


  ! ======================================================================================================= !
  !>
  !> calc_omega
  !> Calculate the livegrass moisture from Oberhagemann et al. (2026), eq 4
  !>
  !> Eq. 6 from Thonicke et al. (2010)
  !>
  FUNCTION calc_livegrass_moisture_ober(nc) RESULT (lvg)

  END FUNCTION calc_livegrass_moisture_ober


  ! ======================================================================================================= !
  !>
  !> calc_phenolgy_ober
  !> Calculate livegrass moisture based on Oberhagemann et al. (2025), related to equation 4.
  !>
  FUNCTION calc_mflg_ober(nc, lfmc_p, comp_p, tair, sw_rad, wscal) RESULT (mf_lg)
    INTEGER, INTENT(in) :: nc
    REAL(wp), INTENT(in) :: tair(nc)     !< Air temperature [deg C]
    REAL(wp), INTENT(in) :: sw_rad(nc)   !< Shortwave-downward radiation [W m-2]
    REAL(wp), INTENT(in) :: wscal(nc)    !< Daily water limiting function value [0 to 1]
    REAL(wp), INTENT(in) :: lfmc_p(3)    !< 3 = number of scaling factors parameters
                                         !< 1=a (min value), 2=b, 3=c
    REAL(wp), INTENT(in) :: comp_p(4, 3) !< 4=number of components(tmin, tmax, light, water), 3=number of function paremeters (slope,base,tau)
    REAL(wp) :: mf_lg(nc) ! Livegrass moisture [0.3 to 2.5]

    REAL(wp) :: f_cold(nc), f_heat(nc), f_water(nc), f_light(nc) !< phenology components
    REAL(wp), ALLOCATABLE, SAVE :: f_cold_old(:), f_heat_old(:), f_water_old(:), f_light_old(:) !< from previous timestep
    REAL(wp) :: klg(nc) !< scaling factor
    REAL(wp), PARAMETER :: klg_max = 2.5_wp

    IF (.NOT. ALLOCATED(f_cold_old)) THEN
      ALLOCATE(f_cold_old(nc), f_heat_old(nc), f_water_old(nc), f_light_old(nc))
      f_cold_old = f_cold
      f_heat_old = f_heat
      f_water_old = f_water
      f_light_old = f_light
    END IF

    !< Calc phenology
    f_cold = calc_f_comp_PFT(nc, comp_p(1,:), f_cold_old, tair)
    f_heat = calc_f_comp_PFT(nc, comp_p(2,:), f_head_old, tair)
    f_water = calc_f_comp_PFT(nc, comp_p(3,:), f_water_old, wscal)
    f_light = calc_f_comp_PFT(nc, comp_p(4,:), f_light_old, sw_rad)
    !< Forkel et al. (2014), eq. 5
    phen = f_cold * f_heat * f_water * f_light

    !< save for next TS
    f_cold_old = f_cold
    f_heat_old = f_heat
    f_water_old = f_water
    f_light_old = f_light

    ! klg scaling factor
    klg = lfmc_p(2) + lfmc_p(3) * phen
    WHERE (klg < lfmc_p(1)) klg = lfmc_p(1) !< set min value
    WHERE (klg > klg_max)   klg = klg_max !< up to 250% livegrass moisture
  END FUNCTION calc_mflg_ober


  ! ======================================================================================================= !
  !>
  !> calc_limit_functions
  !>
  !> Paper: https://bg.copernicus.org/articles/11/7025/2014/bg-11-7025-2014.pdf, eqs. 6, 7, 8, 9
  !> Component limiting function
  !>
  FUNCTION calc_f_comp_PFT(nc, params, f_old, comp) RESULT (f)
    INTEGER, INTENT(in) :: nc
    REAL(wp), INTENT(in) :: params(3) !< slope, base, tau
    REAL(wp), INTENT(in) :: comp(nc) !< physical variable
    REAL(wp) :: f(nc)

    REAL(wp) :: term(nc)

    term = 1.0_wp / (1._wp + exp(-params(1) * (comp - parmas(2))) ) !< 1=slope index, 2=base index
    f = f_old + (term - f_old)  * params(3) !< 3 = tau index

  END FUNCTION calc_f_comp_PFT
#endif
END MODULE mo_q_dfire_sf_fmoist_class
