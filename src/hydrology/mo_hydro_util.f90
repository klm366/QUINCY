!>
!! @brief Contains utilities for the hydrology process
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_hydro_util
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_water_in_root_zone

  INTERFACE get_water_in_root_zone
    MODULE PROCEDURE get_water_in_root_zone_2d
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_util'

CONTAINS

  FUNCTION get_water_in_root_zone_2d(ws, dz, dz_root) RESULT(ws_root)

    REAL(wp), INTENT(in) :: &
      & ws     (:,:),       &        !< Water content of soil layers
      & dz     (:,:),       &        !< Thicknesses of soil layers
      & dz_root(:,:)                 !< Thicknesses of root layers
    REAL(wp) :: ws_root(SIZE(ws,1))  !< Water content in root zone

    REAL(wp) ::                       &
      & bottom     (SIZE(ws     ,1)), & !< Cumulative depth of soil down to a layer (bottom depth of layer)
      & bottom_m1  (SIZE(ws_root,1)), & !< Cumulative depth of soil down to previous layer
      & root_depth (SIZE(ws,1))         !< Total depth of root_zone

    INTEGER :: n, i, j

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_water_in_root_zone_2d'

    !$claw acc routine seq

    n  = SIZE(dz,2)       ! Number of soil layers

    DO j = 1, SIZE(ws, 1)
      bottom (j)    = SUM(dz (j,:))
      root_depth(j) = SUM(dz_root(j,:))
      IF (root_depth(j) > bottom(j) + EPSILON(1._wp)) THEN
        !$claw remove
        CALL finish(TRIM(routine), 'Root zone larger than soil depth.')
        !$claw end remove
      ELSE
        bottom(j)    = 0._wp
        bottom_m1(j) = 0._wp
        ws_root  (j) = 0._wp
        DO i=1,n
          bottom(j)    = bottom(j) + dz(j,i)
          IF (root_depth(j) >= bottom(j)) THEN                ! Root zone extends to below current layer
            ws_root(j) = ws_root(j) + ws(j,i)
          ELSE IF (root_depth(j) > bottom_m1(j)) THEN         ! Root zone above current layer but below previous layer, i.e.
            ws_root(j) = ws_root(j) + ws(j,i) * (root_depth(j) - bottom_m1(j)) / dz(j,i) ! Root z. partially extends into current layer
          END IF
          bottom_m1(j) = bottom_m1(j) + dz(j,i)
        END DO
      END IF
    END DO

  END FUNCTION get_water_in_root_zone_2d

#endif
END MODULE mo_hydro_util
