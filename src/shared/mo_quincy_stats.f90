!> QUINCY functions for statistics
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
!> @brief
!>  Contains common statistics function
!>
!> @author
!>  Albert Jornet Puig, MPI-BGC Jena
!>
MODULE mo_quincy_stats

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish

  IMPLICIT NONE
  PUBLIC :: stats_average_weight, stats_average_weight_1d
  PUBLIC :: stats_average
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_quincy_stats'

  INTERFACE stats_average_weight
    MODULE PROCEDURE :: stats_average_weight_1d
  END INTERFACE stats_average_weight

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Calculates the average weight for 1d arrays
  !>
  !> weights are absolute values. Those are normalized within the function
  !>
  PURE FUNCTION stats_average_weight_1d(sz, values, weights) RESULT (average)
    INTEGER, INTENT(in)   :: sz !<
    REAL(wp), INTENT(in)  :: values(sz)
    REAL(wp), INTENT(in)  :: weights(sz)
    REAL(wp)              :: average
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: total_weights
    REAL(wp) :: weighted(sz), norm_weights(sz)
    ! ----------------------------------------------------------------------------------------------------- !
    weighted = 0._wp
    norm_weights = weights

    !< normalize
    total_weights = sum(weights)
    IF (total_weights /= 1._wp) THEN
      norm_weights = weights / total_weights
    END IF

    weighted(:) = values(:) * norm_weights(:)
    average = SUM(weighted)
  END FUNCTION stats_average_weight_1d


  ! ======================================================================================================= !
  !>
  !> Calculates the average for reals
  !>
  PURE FUNCTION stats_average_r1d(values) RESULT (average)
    REAL(wp), INTENT(in)  :: values(:)
    REAL(wp)              :: average
    ! ----------------------------------------------------------------------------------------------------- !
    IF (SIZE(values) == 0) THEN
      CALL finish("stats_average_r1d", "No data given in values")
    END IF

    IF (SUM(values) == 0._wp) THEN
      average = 0._wp
      RETURN
    END IF

    average = SUM(values) / SIZE(values)
  END FUNCTION stats_average_r1d

END MODULE mo_quincy_stats
