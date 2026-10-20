!> SOC-first convergence arithmetic shared by coupled monitoring and replay diagnostics
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
MODULE mo_sb_convergence_kernel
#ifndef __NO_JSBACH__
  USE mo_kind, ONLY: lp
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: stock_relative_change, soc_window_pass
  PUBLIC :: update_replay_convergence
  ! User-approved provisional SOC criterion (2026-10-19), not an ecological constant.
  REAL(lp), PARAMETER, PUBLIC :: soc_relative_tolerance = 0.001_lp !< Relative cycle change [-]
  REAL(lp), PARAMETER, PUBLIC :: soc_passing_fraction = 0.98_lp !< Required fraction of active cells [-]
CONTAINS
  !> Advance the event-local rolling stock history AFTER one completed replay year.
  !! The caller seeds history(:,0) after any analytical jump and initializes
  !! consecutive=0 and converged=.FALSE. History has N+1 slots so that writing
  !! this year's stock cannot overwrite the N-year comparison endpoint first.
  !! Before N years, cycle_window does not change the convergence decision.
  !! Nutrient changes are diagnostic only; SOC is component 1. Configuration
  !! validation guarantees N>0 and a supported mode before entering this kernel.
  PURE SUBROUTINE update_replay_convergence(mode, n, year, tolerance, stocks, history, &
    & consecutive, converged, previous, relative)
    CHARACTER(len=*), INTENT(in) :: mode !< cycle_window or consecutive_years
    INTEGER, INTENT(in) :: n, year !< Climate-cycle length and completed internal years [-]
    REAL(lp), INTENT(in) :: tolerance !< Strict relative SOC threshold [-]
    REAL(lp), INTENT(in) :: stocks(:) !< Current depth/area-integrated stocks [mol m-2]
    REAL(lp), INTENT(inout) :: history(:,0:) !< Same stocks at circular year endpoints [mol m-2]
    INTEGER, INTENT(inout) :: consecutive !< Number of consecutive passing annual comparisons [-]
    LOGICAL, INTENT(inout) :: converged !< Local replay completion, not a global stop request
    REAL(lp), INTENT(out) :: previous(:), relative(:) !< Comparison endpoint [mol m-2], relative change [-]
    INTEGER :: slot, old_slot !< Circular history slots [-]

    slot = MODULO(year,n+1)
    IF (mode == 'cycle_window') THEN
      old_slot = MODULO(year-n,n+1)
      IF (year >= n) converged = soc_window_pass(history(1,old_slot),stocks(1),tolerance)
    ELSE
      old_slot = MODULO(year-1,n+1)
      IF (soc_window_pass(history(1,old_slot),stocks(1),tolerance)) THEN
        consecutive = consecutive+1
      ELSE
        consecutive = 0
      END IF
      converged = consecutive >= n
    END IF
    previous = history(:,old_slot)
    relative = stock_relative_change(previous,stocks)
    history(:,slot) = stocks
  END SUBROUTINE update_replay_convergence

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
