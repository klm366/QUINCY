!>
!! @brief Contains utilities for the veg component
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2012-09-27)
!!
MODULE mo_veg_util

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: water_in_root_zone

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_util'

CONTAINS

  FUNCTION water_in_root_zone(w, z, zr) RESULT(wr)

    REAL(wp), INTENT(in) :: &
      & w (:,:), &         !< Water content of soil layers
      & z (:,:), &         !< Depths of soil layers
      & zr(:,:)            !< Depths of root layers
    REAL(wp) :: wr(SIZE(w,1))      !< Water content in root zone

    REAL(wp) :: &
      & z_sum    (SIZE(w,1)), &  !< Cumulative depth of soil down to a layer
      & z_sum_m1 (SIZE(w,1)), &  !< Cumulative depth of soil down to previous layer
      & zr_total (SIZE(w,1))     !< Total depth of root_zone

    INTEGER :: &
      & n  ,   &
      & nr
    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':water_in_root_zone'

    n  = SIZE(z,2)    ! Number of soil layers
    nr = SIZE(zr,2)  ! Number of root layers

    z_sum (:) = SUM(z (:,:), DIM=2)
    zr_total(:) = SUM(zr(:,:), DIM=2)
    IF (ANY(zr_total > z_sum)) &
      & CALL finish(TRIM(routine), 'Root zone larger than soil depth.')

    z_sum(:)    = 0._wp
    z_sum_m1(:) = 0._wp
    wr   (:)    = 0._wp
    DO i=1,n
      z_sum(:)    = z_sum(:) + z(:,i)
      WHERE(zr_total(:) > z_sum(:))                                    ! Root zone below current layer
        wr(:) = wr(:) + w(:,i)
      ELSEWHERE (zr_total(:) > z_sum_m1(:))                            ! Root zone above current layer but below previous layer, i.e.
        wr(:) = wr(:) + w(:,i) * (zr_total(:) - z_sum_m1(:)) / z(:,i)  !   root zone partially extends into current layer
      END WHERE
      z_sum_m1(:) = z_sum_m1(:) + z(:,i)
    END DO

  END FUNCTION water_in_root_zone

END MODULE mo_veg_util
