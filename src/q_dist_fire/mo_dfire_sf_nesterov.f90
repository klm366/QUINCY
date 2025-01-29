!> QUINCY Calculate nesterov index
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
!> @brief module to determine the nesterov index (dryness)
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_dfire_sf_nesterov
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message, finish

  !< Frost and precipitation thresholds for Nesterov index
  USE mo_dfire_constants,     ONLY: df_sf_frost_trhd, df_sf_precip_trhd

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_sf_nesterov_index

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_sf_fuel'

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
      &  o_precip_thd, o_frost_thd, o_tsfactor) &
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
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: tdew_daily(nc)          !< dew temperature
    LOGICAL  :: is_frost(nc)
    LOGICAL  :: is_wet_enough(nc)
    REAL(wp) :: ni_inc(nc)              !< nesterov index increase
    REAL(wp) :: precip_thd              !< precipitation threshold
    REAL(wp) :: frost_thd               !< rain as frost threshold
    REAL(wp) :: tsfactor
    ! ----------------------------------------------------------------------------------------------------- !
    precip_thd = df_sf_precip_trhd
    IF (PRESENT(o_precip_thd)) precip_thd = o_precip_thd

    frost_thd = df_sf_frost_trhd
    IF (PRESENT(o_frost_thd)) frost_thd = o_frost_thd

    tsfactor = 1._wp
    IF (PRESENT(o_tsfactor)) tsfactor = o_tsfactor

    is_frost = (tmin_daily <= frost_thd)
    is_wet_enough = (precipitation_daily > precip_thd)

    WHERE (is_wet_enough(:) .OR. is_frost(:))
      ! reset to 0
      nesterov_index(:) = 0._wp
    ELSEWHERE
      tdew_daily(:) = tmin_daily(:) - frost_thd
      ni_inc(:) = (tmax_daily(:) * (tmax_daily(:) - tdew_daily(:))) * tsfactor
      nesterov_index(:) = nesterov_index_prev(:) + ni_inc
    ENDWHERE
  END FUNCTION calc_sf_nesterov_index

#endif
END MODULE mo_dfire_sf_nesterov
