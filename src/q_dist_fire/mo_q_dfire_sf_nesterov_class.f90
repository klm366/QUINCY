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
MODULE mo_q_dfire_sf_nesterov_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message, finish

  !< Frost and precipitation thresholds for Nesterov index
  USE mo_q_dfire_sf_constants, ONLY: df_sf_frost_trhd, df_sf_precip_trhd

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_sf_nesterov_index

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_nesterov_class'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Calculate Nesterov index
  !>
  !> Nesterov Index in Thonicke et al. (2010) is calculated at daily timestep.
  !> We follow the same approach since calculating NI at sub-daily time-step could lead to instabilities.
  !>
  !> Use o_tsfactor to redefine the timestep (e.g to minutes).
  !> Constants can be redefined as well with o_precip_thd and o_ice_thd.
  !> (Section 2.2.2, equation 5 from Thonicke et al. (2010))
  !>
  !> Set default time-step as one day, except if given by optional argument in the function call
  !>
  !> Note: Take into account?
  !> from jsb constants: rain_as_snow_temp_thresh = Tzero + 2.5_wp     ! source: OCN-weather generator
  !>
  !> https://wikifire.wsl.ch/tiki-indexea83.html?page=Nesterov+ignition+index
  !>
  FUNCTION calc_sf_nesterov_index(nc, tmax_daily, tmin_daily, &
                               &  precipitation_daily, nesterov_index_prev, &
                               &  o_precip_thd, o_frost_thd, o_tsfactor, &
                               &  mask) &
                               &  RESULT (nesterov_index)

    INTEGER, INTENT(in)                 :: nc !< dimenion
    REAL(wp), INTENT(in)                :: tmax_daily(nc)           !< daily max temperature [Celsius]
    REAL(wp), INTENT(in)                :: tmin_daily(nc)           !< daily min temperature [Celsius]
    REAL(wp), INTENT(in)                :: precipitation_daily(nc)  !< daily precipitation [mm day-1]
    REAL(wp), DIMENSION(nc), INTENT(in) :: nesterov_index_prev      !< cumulative function [Celsius^2]
    REAL(wp), DIMENSION(nc)             :: nesterov_index           !< output [Celsius^2]
    REAL(wp), OPTIONAL                  :: o_precip_thd             !< redefine precipitation threshold [mm day-1]
    REAL(wp), OPTIONAL                  :: o_frost_thd              !< redefine frosting threshold [Celsius]
    REAL(wp), OPTIONAL                  :: o_tsfactor               !< timestep adjust (daily * tsfactor)
    LOGICAL, OPTIONAL                   :: mask(nc)                 !< True if gridcell is enabled
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: tdew_daily(nc)          !< dew temperature
    LOGICAL  :: is_frost(nc)
    LOGICAL  :: is_wet_enough(nc)
    REAL(wp) :: ni_inc(nc)              !< nesterov index increase
    REAL(wp) :: precip_thd              !< precipitation threshold
    REAL(wp) :: frost_thd               !< rain as frost threshold
    REAL(wp) :: tsfactor
    INTEGER  :: ic
    LOGICAL  :: mask_in(nc)
    ! ----------------------------------------------------------------------------------------------------- !
    precip_thd = df_sf_precip_trhd
    IF (PRESENT(o_precip_thd)) precip_thd = o_precip_thd

    frost_thd = df_sf_frost_trhd
    IF (PRESENT(o_frost_thd)) frost_thd = o_frost_thd

    tsfactor = 1._wp
    IF (PRESENT(o_tsfactor)) tsfactor = o_tsfactor

    mask_in = .TRUE.
    IF (PRESENT(mask)) mask_in = mask

    is_frost = (tmin_daily <= frost_thd)
    is_wet_enough = (precipitation_daily > precip_thd)

    DO ic=1, nc
      IF (is_wet_enough(ic) .OR. is_frost(ic) .OR. .NOT. mask_in(ic)) THEN
        ! reset to 0
        nesterov_index(ic) = 0._wp
      ELSE
        tdew_daily(ic) = tmin_daily(ic) - frost_thd
        ni_inc(ic) = (tmax_daily(ic) * (tmax_daily(ic) - tdew_daily(ic))) * tsfactor
        nesterov_index(ic) = nesterov_index_prev(ic) + ni_inc(ic)
      END IF
    END DO
  END FUNCTION calc_sf_nesterov_index

#endif
END MODULE mo_q_dfire_sf_nesterov_class
