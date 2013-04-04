!>
!! @brief Contains the utilities for the soil component
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2012-09-27)
!!
MODULE mo_soil_util

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: soil_depth_to_layers_2d

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_util'

CONTAINS

  FUNCTION soil_depth_to_layers_2d(soil_depth, dz) RESULT(soil_depth_l)

    REAL(wp), INTENT(in) :: soil_depth(:,:)
    REAL(wp), INTENT(in) :: dz(:)
    REAL(wp)             :: soil_depth_l(SIZE(soil_depth,1), SIZE(dz), SIZE(soil_depth,2))

    REAL(wp) :: dz_cumsum(SIZE(dz)), dz_cumsum_m1(SIZE(dz))
    INTEGER  :: nsoil, i

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_depth_to_layer_2d'

    nsoil = SIZE(dz)

    DO i=2,nsoil
      dz_cumsum(i)    = SUM(dz(1:i))
      dz_cumsum_m1(i) = SUM(dz(1:i-1))
    END DO
    dz_cumsum(1)    = dz(1)
    dz_cumsum_m1(1) = 0._wp

    soil_depth_l(:,:,:) = 0._wp

    DO i=1,nsoil
      WHERE (soil_depth(:,:) > dz_cumsum(i))
        soil_depth_l(:,i,:) = dz(i)
      ELSEWHERE (soil_depth(:,:) > dz_cumsum_m1(i))
        soil_depth_l(:,i,:) = soil_depth(:,:) - dz_cumsum_m1(i)
      END WHERE
    END DO

  END FUNCTION soil_depth_to_layers_2d

END MODULE mo_soil_util
