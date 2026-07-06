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
MODULE mo_q_dfire_sf_frac_consump_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message, finish

  USE mo_q_dfire_sf_constants, ONLY: df_sf_omega_o_idx, df_sf_omega_nl_idx !< omega_o and omega_nl indices
  USE mo_q_dfire_sf_constants, ONLY: df_sf_fc_1hr_idx, df_sf_fc_10hr_idx, & !< fuel consumption indices for each fuel class
                                    df_sf_fc_100hr_idx, df_sf_fc_1000hr_idx

  IMPLICIT NONE

  PRIVATE
  PUBLIC ::  calc_1hr, calc_10hr, calc_100hr, calc_1000hr
  PUBLIC ::  sf_calc_frac_consump

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_frac_consump_class'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Compute fraction of consumption of dead fuel for each fuel class
  !>
  FUNCTION sf_calc_frac_consump(nc, nfuelclasses, moist_extinction, fuel_moisture, mask) RESULT (frac_consump)
    INTEGER, INTENT(in)                   :: nc                            !< Total number of gridcell(s)
    INTEGER, INTENT(in)                   :: nfuelclasses                  !< Total number of fuel classes
    REAL(wp), INTENT(in)                  :: moist_extinction(nc)          !< Moisture of extinction [0 to 1]
    REAL(wp), INTENT(in)                  :: fuel_moisture(nc,nfuelclasses)!< Fuel moisture [%] (#nc,#fuel classes-1,10,100,1000 hours-)
    LOGICAL, INTENT(in)                   :: mask(nc)                      !< True = gridcells to compute
    REAL(wp)                              :: frac_consump(nc,nfuelclasses) !< Fraction of consumption [0 to 1]

    INTEGER :: ic

    frac_consump = 0._wp

    DO ic=1, nc
      IF (mask(ic)) THEN
        frac_consump(ic,df_sf_fc_1hr_idx) = calc_1hr(fuel_moisture(ic,df_sf_fc_1hr_idx), moist_extinction(ic))
        frac_consump(ic,df_sf_fc_10hr_idx) = calc_10hr(fuel_moisture(ic,df_sf_fc_10hr_idx), moist_extinction(ic))
        frac_consump(ic,df_sf_fc_100hr_idx) = calc_100hr(fuel_moisture(ic,df_sf_fc_100hr_idx), moist_extinction(ic))
        frac_consump(ic,df_sf_fc_1000hr_idx) = calc_1000hr(fuel_moisture(ic,df_sf_fc_1000hr_idx), moist_extinction(ic))
      END IF
    END DO
  END FUNCTION sf_calc_frac_consump


  ! ======================================================================================================= !
  !>
  !> Calculate the moisture ratio of fuel consume for 1hr
  !>
  !> Originally using the values  of eq. B1 for 1hr from Thonicke et al. (2010)
  !> Then moved to the Pfeiffer constants values due to a wrong constant.
  !> Finally, the Pfeiffer values have been adjusted to avoid values > 1.
  !> The transition between equations do not properly overlap.
  !>
  !> It it properly tested at the unit tests part.
  !>
  !> Check the script below to visualize its differences:
  !> https://git.bgc-jena.mpg.de/forexd/wp3-modelling-quincy-output-analysis/spitfire/-/blob/main/fraction_consumption.py?ref_type=heads
  !>
  ELEMENTAL FUNCTION calc_1hr(omega_n_l, moist_ext) RESULT (fract_consumption)
    REAL(wp), INTENT(in) :: omega_n_l, moist_ext
    REAL(wp) :: fract_consumption !< output

    REAL(wp) :: ratio
    ratio = omega_n_l / moist_ext
    ratio = min(ratio, 1.0_wp) !< is the fire extinguished?

    IF (ratio <= 0.18_wp) THEN
      fract_consumption = 1.0_wp
      !< Use 0.732 to ensure continuity
    ELSE IF  (ratio > 0.732_wp) THEN !< Pfeiffer et al. (2013): 0.73
      fract_consumption = 2.45_wp - 2.45_wp * ratio
    ELSE ! (0.18_wp > ratio .AND. ratio <= 0.732_wp)
      !< Thonicke et al. (2010): 1.2 (not correct),  Pfeiffer et al. (2013): 1.10 (not perfect)
      fract_consumption = 1.111_wp - 0.62_wp * ratio
    END IF
  END FUNCTION calc_1hr


  ! ======================================================================================================= !
  !>
  !> Calculate the moisture ratio of fuel consume for 10hr
  !>
  !> Par of eq. B4, from Thonicke et al. (2010)
  !> However a constant value has been tuned to have a propper overlap between equations.
  !> Otherwise, small bumps can occur.
  !>
  !> This is properly tested at the unit tests part.
  !>
  ELEMENTAL FUNCTION calc_10hr(omega_o, moist_ext) RESULT (fract_consumption)
    REAL(wp), INTENT(in) :: omega_o, moist_ext
    REAL(wp) :: fract_consumption !< output

    REAL(wp) :: ratio
    ratio = omega_o / moist_ext
    ratio = min(ratio, 1._wp) !< fire extinguished?

    !< Use 0.125 to ensure continuity
    IF (ratio <= 0.125_wp) THEN !<  Thonicke et al. (2010): 0.12, Pfeiffer et al. (2013): 0.12
      fract_consumption = 1.0_wp
    ELSE IF (ratio > 0.51_wp) THEN
      fract_consumption = 1.47_wp - 1.47_wp * ratio
    ELSE !IF (ratio > 0.125_wp .AND. ratio <= 0.51_wp) THEN
      fract_consumption = 1.09_wp - 0.72_wp * ratio
    END IF
  END FUNCTION calc_10hr


  ! ======================================================================================================= !
  !>
  !> Calculate the moisture ratio of fuel consume for 10hr
  !>
  !> Par of eq. B5 from Thonicke et al. (2010)
  !> However a constant value has been tuned to have a propper overlap between equations.
  !> Otherwise, small bumps take place.
  !>
  !> This is properly tested at the unit tests part.
  !>
  ELEMENTAL FUNCTION calc_100hr(omega_o, moist_ext) RESULT (fract_consumption)
    REAL(wp), INTENT(in) :: omega_o, moist_ext
    REAL(wp) :: fract_consumption !< output

    REAL(wp) :: ratio
    ratio = omega_o / moist_ext
    ratio = min(ratio, 1.0_wp)

    !< Use 0.381 to ensure continuity
    IF (ratio <= 0.381_wp) THEN !< Thonicke et al. (2010): 0.38, Pfeiffer et al. (2013): 0.38
      fract_consumption = 0.98_wp - 0.85_wp * ratio
    ELSE !IF (ratio > 0.381_wp) THEN
      fract_consumption = 1.06_wp - 1.06_wp * ratio
    END IF
  END FUNCTION calc_100hr


  ! ======================================================================================================= !
  !>
  !> Calculate the moisture ratio of fuel consume for 10hr
  !> Par of eq. B6 from Thonicke et al. (2010)
  !>
  ELEMENTAL PURE FUNCTION calc_1000hr( omega_o, moist_ext) RESULT (fract_consumption)
    REAL(wp), INTENT(in) :: omega_o, moist_ext
    REAL(wp) :: fract_consumption !< output

    REAL(wp) :: ratio
    ratio = omega_o / moist_ext
    ratio = min(ratio, 1.0_wp)

    fract_consumption = -0.8_wp * ratio + 0.8_wp
  END FUNCTION calc_1000hr

#endif
END MODULE mo_q_dfire_sf_frac_consump_class
