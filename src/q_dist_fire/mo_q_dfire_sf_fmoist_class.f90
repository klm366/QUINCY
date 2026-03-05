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
  USE mo_exception,           ONLY: message, finish
  USE mo_q_dfire_sf_fuel_class, ONLY: sf_calc_has_any_fuel
  USE mo_quincy_stats,        ONLY: stats_average_weight
  USE mo_jsb_lctlib_class,    ONLY: t_lctlib_element
  !< Total number of fuel moisture(omega_o and omega_nl), omega_o index, omega_nl index
  USE mo_q_dfire_constants,   ONLY: df_sf_omegas_nr, df_sf_omega_o_idx, df_sf_omega_nl_idx

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sf_calc_fuel_moisture, sf_calc_fuel_moisture_ober, t_gsi_param, t_gsi, t_lfmc, t_mf, calc_mflg_ober

  !> GSI parameters
  TYPE :: t_gsi_param
    REAL(wp) :: slope
    REAL(wp) :: base
    REAL(wp) :: tau
  END TYPE t_gsi_param

  !> GSI phenology components, from Oberhagemann et al. (2025)
  TYPE :: t_gsi
    TYPE(t_gsi_param) :: cold
    TYPE(t_gsi_param) :: heat
    TYPE(t_gsi_param) :: light
    TYPE(t_gsi_param) :: water
  END TYPE t_gsi

  INTERFACE t_gsi
    MODULE PROCEDURE :: create_gsi_from_lctlib
  END INTERFACE t_gsi

  !> Livegrass fuel moisture content paremeters from Oberhagemann et al. (2025)
  TYPE :: t_lfmc
    REAL(wp) :: min, max, b, c
  END TYPE t_lfmc

  INTERFACE t_lfmc
    MODULE PROCEDURE :: create_lfmc_from_lctlib
  END INTERFACE

  TYPE :: t_mf !< fuel moisture
    !< oberhagemann et al. (2025)
    REAL(wp), ALLOCATABLE :: dead(:,:) !< dead fuel moisture (#nc,#dead fuel moisture -1hr, 10hrs, 100hrs, dry livegrass- )
    REAL(wp), ALLOCATABLE :: live(:,:) !< live fuel moisture (wet livegrass, woody)
    !< Weighted fuel moisture
    REAL(wp), POINTER :: val(:) !< weighted fuel moisture (mixed-1,10,100 and lvg-, fine fuels-1hr and lvg-)
  END TYPE t_mf

  INTERFACE t_mf
    MODULE PROCEDURE :: create_t_mf
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_fmoist_class'
CONTAINS

  !>
  !> t_mf constructor
  !>
  FUNCTION create_t_mf(nc, ndd, nlv) RESULT (mf)
    INTEGER, INTENT(in) :: nc, ndd, nlv
    TYPE(t_mf) :: mf

    ALLOCATE(mf%dead(nc,ndd))
    ALLOCATE(mf%live(nc,nlv))
    mf%dead = 0._wp
    mf%live = 0._wp
  END FUNCTION create_t_mf


  ! ======================================================================================================= !
  !>
  !> t_lfmc constructor
  !>
  FUNCTION create_lfmc_from_lctlib(lctlib) RESULT (lfmc)
    TYPE(t_lctlib_element) :: lctlib
    TYPE(t_lfmc) :: lfmc

    lfmc % min = lctlib%q_lfmc_min
    lfmc % max = lctlib%q_lfmc_max
    lfmc % b = lctlib%q_lfmc_b
    lfmc % c = lctlib%q_lfmc_c
  END FUNCTION create_lfmc_from_lctlib


  ! ======================================================================================================= !
  !>
  !> t_gsi constructor
  !>
  FUNCTION create_gsi_from_lctlib(lctlib) RESULT (gsi)
    TYPE(t_lctlib_element) :: lctlib
    TYPE(t_gsi) :: gsi

    gsi % cold = t_gsi_param( lctlib%q_tmin_slope, lctlib%q_tmin_base, lctlib%q_tmin_tau)
    gsi % heat = t_gsi_param( lctlib%q_tmax_slope, lctlib%q_tmax_base, lctlib%q_tmax_tau)
    gsi % light = t_gsi_param( lctlib%q_light_slope, lctlib%q_light_base, lctlib%q_light_tau)
    gsi % water = t_gsi_param( lctlib%q_water_slope, lctlib%q_water_base, lctlib%q_water_tau)
  END FUNCTION create_gsi_from_lctlib


  ! ======================================================================================================= !
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
  !> Calculate the dead fuel moisture
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
  !> sf_calc_fuel_moisture_ober
  !> Calculate the fuel moisture as in Oberhagemann et al. (2025)
  !>
  FUNCTION sf_calc_fuel_moisture_ober(nc, ndd, nlv, &
                                        sigmas_dd, &
                                        lvg_moisture, &
                                        nesterov_index) &
!                                        lfmc_par, phen) &
                                        RESULT (mf)
    INTEGER, INTENT(in)       :: nc, ndd, nlv
    REAL(wp), INTENT(in)      :: sigmas_dd(:)            !< surface area to voluem ratio for dead and live fuels [cm2/cm3]
    REAL(wp), INTENT(in)      :: lvg_moisture(nc)        !< live grass fuel moisture [%]
    REAL(wp), INTENT(in)      :: nesterov_index(nc)      !< Nesterov index = Dryness index [Celsius^2]
    TYPE(t_mf) :: mf

    REAL(wp), ALLOCATABLE :: mf_lvg(:) !< Live grass fuel moisture [0 to 2.5]%
    REAL(wp), ALLOCATABLE :: alphasdd(:) !<  [deg C-2]
    INTEGER, PARAMETER :: ix_1hr=1, ix_10hr=2, ix_100hr=3, ix_clvg=4  ! litter 1hr, 10 hrs, 100 hrs and cured livegrass
    INTEGER, PARAMETER :: ix_lvg=1, ix_wdy=2                        ! livegrass, woody

!    mf_lvg = calc_mflg_ober(nc, lfmc_par, phen)
    mf = t_mf(nc, ndd, nlv)
    alphasdd = calc_alphas(sigmas_dd)

    mf%dead(:,ix_1hr) = calc_omega(alphasdd(ix_1hr), nesterov_index) ! Same method as Thonicke et al. (2010), but declared obsolete though
    mf%dead(:,ix_10hr) = calc_omega(alphasdd(ix_10hr), nesterov_index) !
    mf%dead(:,ix_100hr) = calc_omega(alphasdd(ix_100hr), nesterov_index) !
    mf%dead(:,ix_clvg) = mf%dead(:,ix_1hr) ! Same as in Oberhagemann et al. (2025), fuelload.c

    mf%live(:,ix_lvg) = lvg_moisture
    mf%live(:,ix_wdy) = 0._wp  ! Same method as Oberhagemann, not yet implemented.
  END FUNCTION sf_calc_fuel_moisture_ober


  ! ======================================================================================================= !
  !>
  !> calc_alphas
  !> Calculate alphas from Thonicke et. al. (2010).
  !> Alphas are the inverse proportion to their Surface-area-to-volume ratios. [deg Celsius-2]
  !>
  FUNCTION calc_alphas(sigmas) RESULT (alphas)
    REAL(wp), INTENT(in)  :: sigmas(:) !< Surface-area-to-volume ratio [cm2/cm3]
    REAL(wp), ALLOCATABLE :: alphas(:) !< Inverse proportionto to their SAV [deg Celsius-2]
    REAL(wp), PARAMETER :: init_alpha = 0.001_wp
    INTEGER :: isig
    REAL(wp) :: ratio

    ratio = 1.
    ALLOCATE(alphas, mold=sigmas)

    alphas(1) = init_alpha
    DO isig = 2, SIZE(sigmas)
      ratio = ratio * sigmas(isig-1) / sigmas(isig)
      alphas(isig) = alphas(1) / ratio
    END DO
  END FUNCTION calc_alphas


  ! ======================================================================================================= !
  !>
  !> calc_mflg_with_gsi_ober
  !> Calculate livegrass moisture based on Oberhagemann et al. (2025), related to equation 4.
  !>
  FUNCTION calc_mflg_with_gsi_ober(nc, lfmc_par, gsi_par, tair, sw_rad, water) RESULT (mf_lg)
    INTEGER, INTENT(in) :: nc
    TYPE(t_lfmc), INTENT(in) :: lfmc_par !< livegrass moisture parameters
    TYPE(t_gsi), INTENT(in) :: gsi_par   !< Parameters used to calculate gsi phenology
    REAL(wp), INTENT(in) :: tair(nc)     !< Air temperature [deg C]
    REAL(wp), INTENT(in) :: sw_rad(nc)   !< Shortwave-downward radiation [W m-2]
    REAL(wp), INTENT(in) :: water(nc)    !< Daily water availability: is a ratio between water
                                         !< supply from soil moisture and atmospheric water demand [0 to 100]%
    REAL(wp), ALLOCATABLE :: mf_lg(:) ! Livegrass moisture [0.3 to 2.5]

    REAL(wp) :: f_cold(nc), f_heat(nc), f_water(nc), f_light(nc) !< phenology components
    REAL(wp), ALLOCATABLE, SAVE :: f_cold_old(:), f_heat_old(:), f_water_old(:), f_light_old(:) !< from previous timestep
    !REAL(wp), PARAMETER :: mflg_max = 2.5_wp
    REAL(wp) :: phen(nc)

    IF (.NOT. ALLOCATED(f_cold_old)) THEN ! initilize
      ALLOCATE(f_cold_old(nc), f_heat_old(nc), f_water_old(nc), f_light_old(nc))
      f_cold_old = calc_term(tair, - gsi_par%cold%slope, gsi_par%cold%base)! slope, base
      f_heat_old = calc_term(tair, gsi_par%heat%slope, gsi_par%heat%base)
      f_light_old = calc_term(sw_rad, -gsi_par%light%slope, gsi_par%light%base)
      f_water_old = calc_term(water, -gsi_par%water%slope, gsi_par%water%base)
    END IF
    ALLOCATE(mf_lg(nc))

    !< Calc phenology
    f_cold = calc_f_comp_PFT(nc, gsi_par%cold, f_cold_old, tair, .FALSE.)
    f_heat = calc_f_comp_PFT(nc, gsi_par%heat, f_heat_old, tair, .TRUE.)
    f_light = calc_f_comp_PFT(nc, gsi_par%light, f_light_old, sw_rad, .FALSE.)
    f_water = calc_f_comp_PFT(nc, gsi_par%water, f_water_old, water, .FALSE.)
    !< Forkel et al. (2014), eq. 5
    phen = f_cold * f_heat * f_water * f_light

    !< save for next TS
    f_cold_old = f_cold
    f_heat_old = f_heat
    f_water_old = f_water
    f_light_old = f_light

    ! As paper's source code, eq. 4, todo: remove and rename to calc_gsi_phen
    mf_lg = lfmc_par%b + lfmc_par%c * phen
    WHERE (mf_lg < lfmc_par%min) mf_lg = lfmc_par%min !< set min value
    WHERE (mf_lg > lfmc_par%max) mf_lg = lfmc_par%max !< up to 250% livegrass moisture
  END FUNCTION calc_mflg_with_gsi_ober


  ! ======================================================================================================= !
  !>
  !> calc_mflg_ober
  !> Calculate livegrass moisture based on Oberhagemann et al. (2025), eq. 4.
  !>
  FUNCTION calc_mflg_ober(nc, lfmc_par, phen) RESULT (mf_lg)
    INTEGER, INTENT(in) :: nc
    TYPE(t_lfmc), INTENT(in) :: lfmc_par  !< livegrass moisture parameters
    REAL(wp), INTENT(in) :: phen(nc)      !< Livegrass phenological status [0 to 1]
    REAL(wp), ALLOCATABLE :: mf_lg(:)     ! Livegrass moisture [0.3 to 2.5]

    ! As paper's source code
    mf_lg = lfmc_par%b + lfmc_par%c * phen
    WHERE (mf_lg < lfmc_par%min) mf_lg = lfmc_par%min !< set min value
    WHERE (mf_lg > lfmc_par%max) mf_lg = lfmc_par%max !< up to 250% livegrass moisture
  END FUNCTION calc_mflg_ober


  ! ======================================================================================================= !
  !>
  !> calc_f_comp_PFT
  !>
  !> Paper: https://bg.copernicus.org/articles/11/7025/2014/bg-11-7025-2014.pdf, eqs. 6, 7, 8, 9
  !> Component limiting function
  !>
  FUNCTION calc_f_comp_PFT(nc, params, f_old, comp, ispositive) RESULT (f)
    INTEGER, INTENT(in)           :: nc
    TYPE(t_gsi_param), INTENT(in) :: params
    REAL(wp), INTENT(in)          :: f_old(nc) !< previous TS value
    REAL(wp), INTENT(in)          :: comp(nc) !< physical variable
    LOGICAL, INTENT(in)           :: ispositive !< true: increase, false=decrease
    REAL(wp) :: f(nc)

    REAL(wp) :: term(nc)
    INTEGER :: ic
    REAL(wp) :: slope
    REAL(wp) :: tau_30min

    IF (ispositive) THEN
      slope = params%slope ! slope increase
    ELSE
      slope = -params%slope ! slope decrease
    ENDIF

    ASSOCIATE(base => params%base, tau_daily => params%tau)
      tau_30min = 1._wp - (1._wp - tau_daily)**(1.0_wp / 48.0_wp)

      DO ic=1, nc
        term(ic) = calc_term(comp(ic), slope, base)
        f(ic) = f_old(ic) + (term(ic) - f_old(ic))  * tau_30min
      ENDDO
    END ASSOCIATE

  END FUNCTION calc_f_comp_PFT


  ! ======================================================================================================= !
  !>
  !> calc_term
  !>
  ELEMENTAL PURE FUNCTION calc_term(comp, slope, base) RESULT (f)
     REAL(wp), INTENT(in) :: comp
     REAL(wp), INTENT(in) :: slope !< sign included here
     REAL(wp), INTENT(in) :: base
     REAL(wp) :: f
     REAL(wp), PARAMETER :: overflow_thd = 200._wp

     IF (slope * (comp - base) < 200) THEN
       f = 1.0_wp / (1._wp + exp(slope * (comp - base))) !< 1=slope index, 2=base index
     ELSE
       f = 0._wp
     ENDIF
  END FUNCTION calc_term
#endif
END MODULE mo_q_dfire_sf_fmoist_class
