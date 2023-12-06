!>
!! @brief Contains utilities for the hydrology process
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_hydro_util
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_water_in_root_zone

  INTERFACE get_water_in_root_zone
    ! MODULE PROCEDURE get_water_in_root_zone_1d
    MODULE PROCEDURE get_water_in_root_zone_2d
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_util'

CONTAINS

  SUBROUTINE get_water_in_root_zone_1d(ws, dz, dz_root, ws_root)

    REAL(wp), INTENT(in) :: &
      & ws     (:),       &        !< Water content of soil layers
      & dz     (:),       &        !< Thicknesses of soil layers
      & dz_root(:)                 !< Thicknesses of root layers
    REAL(wp) :: ws_root            !< Water content in root zone

    REAL(wp) ::                       &
      & bottom    , &   !< Cumulative depth of soil down to a layer (bottom depth of layer)
      & bottom_m1  , &  !< Cumulative depth of soil down to previous layer
      & root_depth      !< Total depth of root_zone

    INTEGER :: n, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_water_in_root_zone_d'

    !$ACC ROUTINE SEQ

    n  = SIZE(dz)       ! Number of soil layers
    bottom    = SUM(dz (:))
    root_depth = SUM(dz_root(:))
#ifndef _OPENACC
    IF (root_depth > bottom + EPSILON(1._wp)) THEN
      CALL finish(TRIM(routine), 'Root zone larger than soil depth.')
    ELSE
#endif
      bottom    = 0._wp
      bottom_m1 = 0._wp
      ws_root   = 0._wp
      DO i=1,n
        bottom    = bottom + dz(i)
        IF (root_depth >= bottom) THEN                ! Root zone extends to below current layer
          ws_root = ws_root + ws(i)
        ELSE IF (root_depth > bottom_m1) THEN         ! Root zone above current layer but below previous layer, i.e.
          ws_root = ws_root + ws(i) * (root_depth - bottom_m1) / dz(i) ! Root z. partially extends into current layer
        END IF
        bottom_m1 = bottom_m1 + dz(i)
      END DO
#ifndef _OPENACC
    END IF
#endif

  END SUBROUTINE get_water_in_root_zone_1d



  SUBROUTINE get_water_in_root_zone_2d(ws, dz, dz_root, ws_root)

    REAL(wp), INTENT(in) :: &
      & ws     (:,:),       &        !< Water content of soil layers
      & dz     (:,:),       &        !< Thicknesses of soil layers
      & dz_root(:,:)                 !< Thicknesses of root layers

    REAL(wp), INTENT(inout) :: ws_root(SIZE(ws,1))  !< Water content in root zone

    REAL(wp) ::                       &
      & bottom     (SIZE(ws     ,1)), & !< Cumulative depth of soil down to a layer (bottom depth of layer)
      & bottom_m1  (SIZE(ws_root,1)), & !< Cumulative depth of soil down to previous layer
      & root_depth (SIZE(ws,1))         !< Total depth of root_zone

    INTEGER :: n, nc, i, j

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_water_in_root_zone_2d'

    nc = SIZE(ws,1)
    n  = SIZE(dz,2)       ! Number of soil layers

    !$ACC DATA CREATE(bottom, bottom_m1, root_depth)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO j = 1, nc
      root_depth(j) = 0._wp
        bottom(j)    = 0._wp
        bottom_m1(j) = 0._wp
        ws_root  (j) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

    !$ACC LOOP SEQ
    DO i=1,n
      !$ACC LOOP GANG(STATIC:1) VECTOR
      DO j = 1, nc
        root_depth(j) = root_depth(j) + dz_root(j,i)
      END DO
    END DO

    !$ACC LOOP SEQ
    DO i=1,n
      !$ACC LOOP GANG(STATIC:1) VECTOR
      DO j = 1, nc
          bottom(j)    = bottom(j) + dz(j,i)
          IF (root_depth(j) >= bottom(j)) THEN                ! Root zone extends to below current layer
            ws_root(j) = ws_root(j) + ws(j,i)
          ELSE IF (root_depth(j) > bottom_m1(j)) THEN         ! Root zone above current layer but below previous layer, i.e.
            ws_root(j) = ws_root(j) + ws(j,i) * (root_depth(j) - bottom_m1(j)) / dz(j,i) ! Root z. partially extends into current layer
          END IF
          bottom_m1(j) = bottom_m1(j) + dz(j,i)
        END DO
    END DO

    !$ACC END PARALLEL

#ifndef _OPENACC
    DO j = 1, nc
      IF (root_depth(j) > bottom(j) + EPSILON(1._wp)) THEN
        CALL finish(TRIM(routine), 'Root zone larger than soil depth.')
      END IF
    END DO
#endif

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE get_water_in_root_zone_2d

#endif
END MODULE mo_hydro_util
