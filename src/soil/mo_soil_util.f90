!>
!! @brief Contains the utilities for the soil process
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! Second version by Jan Engel (20161208 jsbach4:git) - integration into SLM, no code changes though
!!
MODULE mo_soil_util
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: soil_depth_to_layers_2d

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_util'

CONTAINS

  FUNCTION soil_depth_to_layers_2d(soil_depth, dz) RESULT(soil_depth_sl)

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: soil_depth(:,:)
    REAL(wp), INTENT(in) :: dz(:)
    REAL(wp)             :: soil_depth_sl(SIZE(soil_depth,1), SIZE(dz), SIZE(soil_depth,2))
    ! 0.2 Local
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

    soil_depth_sl(:,:,:) = 0._wp

    DO i=1,nsoil
      WHERE (soil_depth(:,:) > dz_cumsum(i))
        soil_depth_sl(:,i,:) = dz(i)
      ELSEWHERE (soil_depth(:,:) > dz_cumsum_m1(i))
        soil_depth_sl(:,i,:) = soil_depth(:,:) - dz_cumsum_m1(i)
      END WHERE
    END DO

  END FUNCTION soil_depth_to_layers_2d

#endif
END MODULE mo_soil_util
