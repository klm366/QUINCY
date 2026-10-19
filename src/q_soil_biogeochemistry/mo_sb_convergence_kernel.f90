!> SOC-first convergence arithmetic shared by coupled monitoring and replay diagnostics
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
MODULE mo_sb_convergence_kernel
#ifndef __NO_JSBACH__
  USE mo_kind, ONLY: lp
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: stock_relative_change, soc_window_pass
  ! User-approved provisional SOC criterion (2026-10-19), not an ecological constant.
  REAL(lp), PARAMETER, PUBLIC :: soc_relative_tolerance = 0.001_lp !< Relative cycle change [-]
  REAL(lp), PARAMETER, PUBLIC :: soc_passing_fraction = 0.98_lp !< Required fraction of active cells [-]
CONTAINS
  !> Relative endpoint change against the OLD stock, without annualization.
  !! Stocks use identical units/area bases. Zero-to-zero passes; a positive stock
  !! emerging from zero is not assigned an artificial denominator. Negative or
  !! nonfinite stocks are invalid. HUGE is a reportable failure sentinel, not a flux.
  PURE ELEMENTAL FUNCTION stock_relative_change(previous, current) RESULT(change)
    REAL(lp), INTENT(in) :: previous, current !< Endpoint stocks [same mass/area units]
    REAL(lp) :: change !< Absolute relative change, or HUGE for invalid/unbounded [-]
    change = HUGE(1._lp)
    IF (.NOT. ieee_is_finite(previous) .OR. .NOT. ieee_is_finite(current)) RETURN
    IF (MIN(previous,current) < 0._lp) RETURN
    IF (previous == 0._lp) THEN
      IF (current == 0._lp) change = 0._lp
      RETURN
    END IF
    ! Avoid overflow when a tiny reference is followed by a large positive stock.
    IF (previous < 1._lp) THEN
      IF (ABS(current-previous) > HUGE(1._lp)*previous) RETURN
    END IF
    change = ABS(current-previous)/previous
  END FUNCTION stock_relative_change

  PURE ELEMENTAL FUNCTION soc_window_pass(previous, current, tolerance) RESULT(passed)
    REAL(lp), INTENT(in) :: previous, current !< Endpoint SOC stocks [same mass/area units]
    REAL(lp), INTENT(in) :: tolerance !< Strict relative threshold [-]
    LOGICAL :: passed !< This window passes, not a permanently latched decision [-]
    passed = .FALSE.
    IF (.NOT. ieee_is_finite(tolerance) .OR. tolerance <= 0._lp) RETURN
    passed = stock_relative_change(previous,current) < tolerance
  END FUNCTION soc_window_pass
#endif
END MODULE mo_sb_convergence_kernel
