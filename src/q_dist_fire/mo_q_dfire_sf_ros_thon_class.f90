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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### routine to calculate fire rate of spread from Thonicke et al. (2010)
!>
!>  Most of the equations are based on Thonicke et al. (2010). Nonetheless, some equations are also
!>  using Pfeiffer et al. (2013). Furter improvements are applied using Oberhagemann et al. (2025)
!>  paper. Those includes weighting fuels accordingly.
!>
!>  This involves equations from:
!>    - Section 2.2.4, Equations 9 to 14
!>    - Table A1
!>    - Appendix A (equations A1 to A9)
!>
MODULE mo_q_dfire_sf_ros_thon_class
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: lp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  !< Total number of sigmas (1,10,100) hours,  Pfeiffer el al. (2013), Eq. A18 - Used to calculate fuel bulk density
  USE mo_q_dfire_sf_constants,  ONLY: df_sf_sigma_tho_nr
  !<
  !<  Rate of spread constants
  USE mo_q_dfire_sf_constants,  ONLY: df_sf_oven_dry_part_dens, df_sigmas_tho, &
                                      df_sf_mineral_fract, df_sf_mineral_damp_coeff, &
                                      df_sf_fuel_heat_content

  USE mo_quincy_stats,          ONLY: stats_average_weight
  USE mo_q_dfire_sf_util,       ONLY: calc_ros_opt_react_velocity, &
                                       calc_ros_prop_flux_ratio, &
                                       calc_ros_wind_factor
  USE mo_q_dfire_sf_classes,    ONLY: t_sf_ros


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sf_calc_rate_of_spread_thon
  PUBLIC :: calc_sigma, calc_moist_damp_coeff

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_ros_thon_class'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> calc
  !>
  !> Do the science
  !>
  SUBROUTINE sf_calc_rate_of_spread_thon(nc, nfuels, fuel_bulk_density, u_f_wind, moist_extinction, &
                                        w_oi, omega_o, mask, ros)
    INTEGER, INTENT(in)                  :: nc
    INTEGER, INTENT(in)                  :: nfuels                      !< 1, 10 and 100 hours fuel classes
    REAL(lp), INTENT(in)                 :: fuel_bulk_density(nc)       !< [kg m-3]
    REAL(lp), INTENT(in)                 :: u_f_wind(nc)                !< [m min-1]
    REAL(lp), INTENT(in)                 :: moist_extinction(nc)        !< [0 to 1]
    REAL(lp), INTENT(in)                 :: w_oi(nc, nfuels)            !< [kg DM m-2] Note in kg
    !< According to Pfeiffer et al. (2013), omega_n must be treated as omega_o
    REAL(lp), INTENT(in)                 :: omega_o(nc)                 !< [0 to 1] (#nc)
    LOGICAL, INTENT(in)                  :: mask(nc)                    !< True =  compute gridcell(s)
    TYPE(t_sf_ros), INTENT(INOUT)        :: ros                         !< output [-]

    ! ----------------------------------------------------------------------------------------------------- !
    REAL(lp) :: beta(nc), beta_op(nc)
    REAL(lp) :: sigma(nc)                 !< Surface-area-to-volume-ratio [cm-1]
    REAL(lp) :: A_var(nc)
    REAL(lp) :: max_reaction_velocity(nc)

    INTEGER :: inc
    ! ----------------------------------------------------------------------------------------------------- !
    ros%mean_fire_area = 0._lp

    DO inc = 1, nc
      IF (mask(inc)) THEN
        !<  (Table A1, 5th from Thonicke et al. (2010))
        ros%net_fuel_load(inc) = (1._lp - df_sf_mineral_fract) * SUM(w_oi(inc,:))
        ros%fbd(inc) = fuel_bulk_density(inc) !< fuel bulk density [kg m-3], precomputed at mo_q_dfire_fuel.f90 + aggregation to         average

        sigma(inc) = calc_sigma(df_sf_sigma_tho_nr, nfuels, df_sigmas_tho, w_oi(inc,:)) !< Table A1, 3rd equation [cm-1] from Thonicke et al. (2010)
        !< (Table A1, 3rd from Thonicke et al. (2010))
        beta(inc) =  ros%fbd(inc) / df_sf_oven_dry_part_dens !< [-]
        !< Calculate beta_op (Eq. A6 from Thonicke et al. (2010))
        beta_op(inc) = 0.200395_lp * (sigma(inc) ** (-0.8189_lp)) !< [-]

        A_var(inc) = 8.9033_lp * (sigma(inc) ** (-0.7913)) !< Table A1, 4th equation from Thonicke et al. (2010)
        max_reaction_velocity(inc) = 1._lp / (0.0591_lp  + 2.926_lp * sigma(inc) ** (-1.5_lp)) !< Table A1, 2nd equation from Thonicke et al. (2010) [min-1]
        ros%opt_react_velo(inc) = calc_ros_opt_react_velocity(max_reaction_velocity(inc), beta(inc), &
            beta_op(inc), A_var(inc)) ! Table A1, 1st equation from Thonicke et al. (2010) [min-1]

        ros%moist_damp_coeff(inc) = calc_moist_damp_coeff(omega_o(inc), moist_extinction(inc)) !< Table A1, 6th equation from Thonicke et al. (2010)
        ros%mineral_damp_coeff(inc) = df_sf_mineral_damp_coeff

        !< (Eq. A1 from Thonicke et al. (2010))
        ros%reaction_intensity(inc) = ros%opt_react_velo(inc) * &
                                      ros%net_fuel_load(inc) * df_sf_fuel_heat_content * &
                                      ros%moist_damp_coeff(inc) * df_sf_mineral_damp_coeff !< [kJ m-2 min-1]

        ros%prop_flux_ratio(inc) = calc_ros_prop_flux_ratio(sigma(inc), beta(inc)) !< Eq. A2 from Thonicke et al. (2010) [-]
        ros%wind_factor(inc) = calc_ros_wind_factor(u_f_wind(inc), sigma(inc), beta(inc), beta_op(inc)) !< Eq. A5 from Thonicke et al. (2010) [-]
        ros%eff_heating_number(inc) = exp(-4.528_lp / sigma(inc)) !< (Eq. A3 from Thonicke et al. (2010)) [-]
        ros%pre_ign_heat(inc) = 581._lp + 2594._lp * omega_o(inc) !< (Eq. A4 from Thonicke et al. (2010)) [kJ kg-1]

        !< ( Eq. 9 from Thonicke et al. (2010))
        ros%ros_forward(inc) = ( ros%reaction_intensity(inc)  * ros%prop_flux_ratio(inc) * &
                                (1 + ros%wind_factor(inc)) ) / &
                                ( ros%fbd(inc) * ros%eff_heating_number(inc) * ros%pre_ign_heat(inc) ) !< [m min-1]
      END IF ! IF (mask(inc)) THEN
    ENDDO ! DO inc = 1, nc

  END SUBROUTINE sf_calc_rate_of_spread_thon


  ! ======================================================================================================= !
  !>
  !> Calculate the moisture dampening coefficient (Table A1, 6th from Thonicke et al. (2010))
  !>
  !> when equal to 1, the fuel is moisture-free
  !>
  ELEMENTAL PURE FUNCTION calc_moist_damp_coeff(omega_o, moist_extinction) RESULT (moist_damp_coeff)
    REAL(lp), INTENT(in) :: omega_o           !< Fuel moisture [0 to 1]
    REAL(lp), INTENT(in) :: moist_extinction  !< Moisture of extinction [0 to 1]
    REAL(lp)             :: moist_damp_coeff  !< output [-]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(lp) :: ratio
    ! ----------------------------------------------------------------------------------------------------- !
    ratio = omega_o / moist_extinction
    moist_damp_coeff = 1._lp - (2.59_lp * ratio) + (5.11_lp * (ratio ** 2_lp)) - (3.52_lp * (ratio ** 3_lp))
  END FUNCTION calc_moist_damp_coeff


  ! ======================================================================================================= !
  !>
  !> Calculate sigma (Table A1, 2nd form Thonicke et al. (2010))
  !>
  FUNCTION calc_sigma(nsigma, nfuels, sigmas, woi) RESULT (sigma)
    INTEGER, INTENT(in)  :: nsigma, nfuels                    !<
    REAL(lp), INTENT(in) :: sigmas(nsigma)                    !<
    REAL(lp), INTENT(in) :: woi(nfuels)                       !< dead fuel [kg m-2]
    REAL(lp)             :: sigma                             !< output [cm-1]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(lp)             :: total_woi
    ! ----------------------------------------------------------------------------------------------------- !
    total_woi = sum(woi(:))
    IF (total_woi < EPSILON(0._lp)) THEN ! woi == 0.
      sigma = 0.0001_lp
    ELSE
      sigma = stats_average_weight(nsigma, sigmas(:), woi(:))
    END IF
  END FUNCTION calc_sigma


#endif
END MODULE mo_q_dfire_sf_ros_thon_class
