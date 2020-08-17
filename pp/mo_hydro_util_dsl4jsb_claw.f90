MODULE mo_hydro_util
 USE mo_kind , ONLY: wp
 USE mo_exception , ONLY: finish
 PRIVATE :: get_water_in_root_zone_1d
 PRIVATE :: get_water_in_root_zone_2d
 PUBLIC :: get_water_in_root_zone
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_hydro_util"
 INTERFACE get_water_in_root_zone
  MODULE PROCEDURE get_water_in_root_zone_2d
 END INTERFACE

CONTAINS
 SUBROUTINE get_water_in_root_zone_1d ( ws , dz , dz_root , ws_root )

  REAL ( KIND= wp ) , INTENT(IN) :: ws ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: dz ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: dz_root ( : )
  REAL ( KIND= wp ) :: ws_root
  REAL ( KIND= wp ) :: bottom
  REAL ( KIND= wp ) :: bottom_m1
  REAL ( KIND= wp ) :: root_depth
  INTEGER :: n
  INTEGER :: i
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_util" //&
  ":get_water_in_root_zone_d"

!$acc routine seq
  n = size ( dz )
  bottom = sum ( dz ( : ) )
  root_depth = sum ( dz_root ( : ) )
  IF ( root_depth > bottom + epsilon ( 1._wp ) ) THEN
  ELSE
   bottom = 0._wp
   bottom_m1 = 0._wp
   ws_root = 0._wp
   DO i = 1 , n , 1
    bottom = bottom + dz ( i )
    IF ( root_depth >= bottom ) THEN
     ws_root = ws_root + ws ( i )
    ELSE
     IF ( root_depth > bottom_m1 ) THEN
      ws_root = ws_root + ws ( i ) * ( root_depth - bottom_m1 ) / dz ( i )
     END IF
    END IF
    bottom_m1 = bottom_m1 + dz ( i )
   END DO
  END IF
 END SUBROUTINE get_water_in_root_zone_1d

 SUBROUTINE get_water_in_root_zone_2d ( ws , dz , dz_root , ws_root )

  REAL ( KIND= wp ) , INTENT(IN) :: ws ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: dz ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: dz_root ( : , : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: ws_root ( 1 : size ( ws , 1 ) )
  REAL ( KIND= wp ) :: bottom ( 1 : size ( ws , 1 ) )
  REAL ( KIND= wp ) :: bottom_m1 ( 1 : size ( ws_root , 1 ) )
  REAL ( KIND= wp ) :: root_depth ( 1 : size ( ws , 1 ) )
  INTEGER :: n
  INTEGER :: i
  INTEGER :: j
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_util" //&
  ":get_water_in_root_zone_2d"

  n = size ( dz , 2 )
!$acc data present(ws, dz, dz_root, ws_root) create(bottom, bottom_m1, root_depth)
!$acc parallel default(none)
!$acc loop gang vector
  DO j = 1 , size ( ws , 1 ) , 1
   bottom ( j ) = sum ( dz ( j , : ) )
   root_depth ( j ) = sum ( dz_root ( j , : ) )
   IF ( root_depth ( j ) > bottom ( j ) + epsilon ( 1._wp ) ) THEN
   ELSE
    bottom ( j ) = 0._wp
    bottom_m1 ( j ) = 0._wp
    ws_root ( j ) = 0._wp
!$acc loop seq
    DO i = 1 , n , 1
     bottom ( j ) = bottom ( j ) + dz ( j , i )
     IF ( root_depth ( j ) >= bottom ( j ) ) THEN
      ws_root ( j ) = ws_root ( j ) + ws ( j , i )
     ELSE
      IF ( root_depth ( j ) > bottom_m1 ( j ) ) THEN
       ws_root ( j ) = ws_root ( j ) + ws ( j , i ) * ( root_depth ( j ) -&
        bottom_m1 ( j ) ) / dz ( j , i )
      END IF
     END IF
     bottom_m1 ( j ) = bottom_m1 ( j ) + dz ( j , i )
    END DO
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE get_water_in_root_zone_2d

END MODULE mo_hydro_util

