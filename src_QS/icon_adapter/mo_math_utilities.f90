! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

!>
!! Contains the implementation of various mathematical algorithms
!! used by the shallow water model.
!!

! ---------------------------------------------------------------
! LL: xlc has trouble optimizing routines with implicit shaped parameters
! #ifdef __xlC__
! @process nohot
! ! @PROCESS NOOPTIMIZE
! #endif
MODULE mo_math_utilities

  USE mo_kind,    ONLY: wp
  ! USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64
  ! USE mo_math_constants, ONLY: pi, pi_2, dbl_eps
  ! USE mo_gridman_constants, ONLY: SUCCESS, TORUS_MAX_LAT
  ! USE mo_exception, ONLY: finish
  ! USE mo_lib_grid_geometry_info, ONLY: t_grid_geometry_info, planar_torus_geometry, sphere_geometry
  ! USE mo_math_types, ONLY: t_cartesian_coordinates, t_geographical_coordinates, &
  !   &                                   t_line, t_tangent_vectors
  ! USE mo_util_sort, ONLY: quicksort
  IMPLICIT NONE

  PRIVATE

!   PUBLIC :: cc2gc
!   PUBLIC :: gc2cc
!   PUBLIC :: cc2tv
!   PUBLIC :: gvec2cvec
!   PUBLIC :: cvec2gvec
!   PUBLIC :: arc_length
!   PUBLIC :: arc_length_v
!   PUBLIC :: cartesian_to_geographical, geographical_to_cartesian
!   PUBLIC :: arc_length_on_unitsphere
!   PUBLIC :: cos_arc_length
!   PUBLIC :: plane_torus_distance, plane_torus_closest_coordinates
!   PUBLIC :: plane_angle_of_three_points
!   PUBLIC :: sphere_cartesian_midpoint
!   PUBLIC :: norma, normalize
!   PUBLIC :: angle_of_vectors, sin_cc
!   PUBLIC :: triangle_area
!   PUBLIC :: circum_center
!   PUBLIC :: bary_center
!   PUBLIC :: middle
!   PUBLIC :: sphere_tangent_coordinates
!   PUBLIC :: normal_vector
!   PUBLIC :: spherical_intersection, spherical_intersection2, spherical_intersection2_v
!   PUBLIC :: rotate_z, rotate_x, rotate_y
!   PUBLIC :: norma_of_vector_product
!   PUBLIC :: project_point_to_plane

!   PUBLIC :: disp_new
!   PUBLIC :: disp_new_vect

!   PUBLIC :: cc_dot_product
!   PUBLIC :: cc_norm
!   PUBLIC :: vector_product
!   PUBLIC :: integral_over_triangle
!   PUBLIC :: rotate_latlon
!   PUBLIC :: rotate_latlon_vec

!   PUBLIC :: gnomonic_proj
!   PUBLIC :: orthogr_proj
!   PUBLIC :: az_eqdist_proj
! !   PUBLIC :: sphere_cell_mean_char_length
!   PUBLIC :: ccw
!   PUBLIC :: line_intersect
!   PUBLIC :: lintersect
  PUBLIC :: tdma_solver
  PUBLIC :: tdma_solver_vec
!   PUBLIC :: check_orientation

!   !  vertical coordinates routines
!   PUBLIC :: set_zlev

!   ! subroutines for value sets
!   PUBLIC :: t_value_set
!   PUBLIC :: merge_values_into_set
!   PUBLIC :: find_values_in_set
!   PUBLIC :: deallocate_set

!   PUBLIC :: OPERATOR(+)
!   PUBLIC :: OPERATOR(-)
!   PUBLIC :: OPERATOR(*)

!   INTERFACE OPERATOR(+)
!     MODULE PROCEDURE cartesian_coordinates_plus
!   END INTERFACE
!   INTERFACE OPERATOR(-)
!     MODULE PROCEDURE cartesian_coordinates_minus
!   END INTERFACE
!   INTERFACE OPERATOR(*)
!     MODULE PROCEDURE cartesian_coordinates_prod
!   END INTERFACE

!   INTERFACE arc_length
!     MODULE PROCEDURE arc_length_sphere
!     MODULE PROCEDURE arc_length_generic
!   END INTERFACE
!   INTERFACE arc_length_v
!     MODULE PROCEDURE arc_length_v_sphere
!     MODULE PROCEDURE arc_length_v_generic
!   END INTERFACE
!   INTERFACE gc2cc
!     MODULE PROCEDURE gc2cc_sphere
!     MODULE PROCEDURE gc2cc_generic
!   END INTERFACE
!   INTERFACE cc2gc
!     MODULE PROCEDURE cc2gc_sphere
!     MODULE PROCEDURE cc2gc_generic
!   END INTERFACE

!   !> Derived type specifying a set of REAL(wp) values.
!   !
!   TYPE t_value_set
!     INTEGER               :: nvalues !< no. of values in set
!     LOGICAL               :: sort_smallest_first !< Flag. If true, set is sorted smallest-to-largest
!     REAL(wp), ALLOCATABLE :: values(:) !< sorted set of values.

!     !> Everytime when new values are inserted into the set, the total
!     !  list of values is re-sorted.
!   END TYPE t_value_set

!   REAL(wp), PARAMETER :: fcirc = 2.0_wp*pi
!   REAL(wp), PARAMETER :: hcirc = fcirc*0.5_wp, qcirc = fcirc*0.25_wp
!   ! constants for fixpoint<->float conversion of latitude/longitude angles
!   REAL(wp), PARAMETER :: lat_spread = 0.25_wp*REAL(HUGE(1), wp)/qcirc, &
!                          lat_contract = 4.0_wp*qcirc/REAL(HUGE(1), wp), &
!                          lon_spread = 0.25_wp*REAL(HUGE(1), wp)/hcirc, &
!                          lon_contract = 4.0_wp*hcirc/REAL(HUGE(1), wp)

!   PUBLIC :: normalized_lat, fxp_lat, flp_lat
!   PUBLIC :: normalized_lon, fxp_lon, flp_lon

! !-----------------------------------------------------------------------
! ! Basic geometry functions definitions
! #define d_norma_3d(v) SQRT(DOT_PRODUCT(v%x,v%x))
! #define d_norma_3d_x(x) SQRT(DOT_PRODUCT(x,x))
! #define d_normalize(v) v%x=v%x/d_norma_3d(v)
! #define d_normalize_x(x) x=x/d_norma_3d_x(x)
! #define d_arc_of_hord_normalsphere(l) (2.0_wp*ASIN((l)*0.5_wp))
! #define d_normalize_f(v) v%x/d_norma_3d(v)
! #define d_sqrdistance_3d(v1,v2) DOT_PRODUCT((v1%x-v2%x),(v1%x-v2%x))

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_math_utilities'

CONTAINS

  !! 
  !! quincy: code parts deleted because some docu created issues with dsl python script
  !!

  !-------------------------------------------------------------------------
  !>
  !! TDMA tridiagonal matrix solver for a_i*x_(i-1) + b_i*x_i + c_i*x_(i+1) = d_i
  !!
  !!       a - sub-diagonal (means it is the diagonal below the main diagonal)
  !!       b - the main diagonal
  !!       c - sup-diagonal (means it is the diagonal above the main diagonal)
  !!       d - right part
  !!       x - the answer
  !!       n - number of equations
  SUBROUTINE tdma_solver(a, b, c, d, n, varout)

    !$ACC ROUTINE SEQ

    INTEGER, INTENT(in) :: n
    REAL(wp), DIMENSION(n), INTENT(in)  :: a, b, c, d
    REAL(wp), DIMENSION(n), INTENT(out) :: varout

    REAL(wp):: m, cp(n), dp(n)
    INTEGER :: i

! initialize c-prime and d-prime
    cp(1) = c(1)/b(1)
    dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
    DO i = 2, n
      m = b(i) - cp(i - 1)*a(i)
      cp(i) = c(i)/m
      dp(i) = (d(i) - dp(i - 1)*a(i))/m
    END DO
! initialize varout
    varout(n) = dp(n)
! solve for varout from the vectors c-prime and d-prime
    DO i = n - 1, 1, -1
      varout(i) = dp(i) - cp(i)*varout(i + 1)
    END DO

  END SUBROUTINE tdma_solver

  !-------------------------------------------------------------------------
  !>
  !! TDMA tridiagonal matrix solver for a_i*x_(i-1) + b_i*x_i + c_i*x_(i+1) = d_i
  !!
  !!       a - sub-diagonal (means it is the diagonal below the main diagonal)
  !!       b - the main diagonal
  !!       c - sup-diagonal (means it is the diagonal above the main diagonal)
  !!       d - right part
  !!  varout - the answer (identical to x in description above)
  !!    slev - start level (top)
  !!    elev - end level (bottom)
  SUBROUTINE tdma_solver_vec(a, b, c, d, slev, elev, startidx, endidx, varout, opt_acc_queue)
    INTEGER, INTENT(IN) :: slev, elev
    INTEGER, INTENT(IN) :: startidx, endidx
    REAL(wp), INTENT(IN) :: a(:, :), b(:, :), c(:, :), d(:, :)
    REAL(wp), INTENT(OUT) :: varout(:, :)
    INTEGER, OPTIONAL, INTENT(IN) :: opt_acc_queue

    !
    ! local
    REAL(wp):: m, cp(SIZE(a, 1), SIZE(a, 2)), dp(SIZE(a, 1), SIZE(a, 2))
    INTEGER :: i
    INTEGER :: jc
    INTEGER :: acc_queue

    IF (PRESENT(opt_acc_queue)) THEN
      acc_queue = opt_acc_queue
    ELSE
      acc_queue = 1
    END IF

    ! initialize c-prime and d-prime
    !$ACC PARALLEL DEFAULT(PRESENT) CREATE(cp, dp) ASYNC(acc_queue)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jc = startidx, endidx
      cp(jc, slev) = c(jc, slev)/b(jc, slev)
      dp(jc, slev) = d(jc, slev)/b(jc, slev)
    END DO
    ! solve for vectors c-prime and d-prime
    !$ACC LOOP SEQ
!NEC$ outerloop_unroll(4)
    DO i = slev + 1, elev
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(m)
      DO jc = startidx, endidx
        m = 1._wp/(b(jc, i) - cp(jc, i - 1)*a(jc, i))
        cp(jc, i) = c(jc, i)*m
        dp(jc, i) = (d(jc, i) - dp(jc, i - 1)*a(jc, i))*m
      END DO
    END DO
    ! initialize varout
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jc = startidx, endidx
      varout(jc, elev) = dp(jc, elev)
    END DO
    ! solve for varout from the vectors c-prime and d-prime
    !$ACC LOOP SEQ
!NEC$ outerloop_unroll(4)
    DO i = elev - 1, slev, -1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = startidx, endidx
        varout(jc, i) = dp(jc, i) - cp(jc, i)*varout(jc, i + 1)
      END DO
    END DO
    !$ACC END PARALLEL

    IF (.NOT. PRESENT(opt_acc_queue)) THEN
      !$ACC WAIT(acc_queue)
    END IF

  END SUBROUTINE tdma_solver_vec
  !-------------------------------------------------------------------------

  ! !-------------------------------------------------------------------------
  ! !
  ! !!! Helper functions for computing the vertical layer structure
  ! SUBROUTINE set_zlev(zlev_i, zlev_m, n_zlev, dzlev_m)
  !   INTEGER, INTENT(IN)  :: n_zlev
  !   REAL(wp), INTENT(OUT) :: zlev_i(n_zlev + 1), zlev_m(n_zlev)
  !   REAL(wp), INTENT(IN)  :: dzlev_m(:) ! namelist input of layer thickness

  !   INTEGER :: jk
  !   REAL(wp) :: accum

  !   accum = 0.0_wp
  !   ! zlev_i    : upper border surface of vertical cells
  !   ! zlev_m    : position of coordinate surfaces in meters below zero surface.
  !   DO jk = 1, n_zlev
  !     zlev_i(jk) = accum
  !     zlev_m(jk) = 0.5_wp*(2.0_wp*accum + dzlev_m(jk))
  !     accum = accum + dzlev_m(jk)
  !   END DO
  !   zlev_i(n_zlev + 1) = accum

  ! END SUBROUTINE set_zlev

  ! !> normalize latitude to range [-pi/2,pi/2]
  ! ELEMENTAL FUNCTION normalized_lat(lat) RESULT(normalized)
  !   REAL(wp), INTENT(in) :: lat
  !   REAL(wp) :: normalized
  !   REAL(wp) :: norm_hcirc, qcirc_div, offset, sign_mult
  !   offset = SIGN(qcirc, lat)
  !   qcirc_div = AINT((lat + offset)/qcirc)
  !   norm_hcirc = MOD(lat + offset, hcirc) - offset
  !   sign_mult = -2.0_wp*ABS(MOD(AINT(qcirc_div*.5_wp), 2.0_wp)) + 1.0_wp
  !   normalized = sign_mult*norm_hcirc
  ! END FUNCTION normalized_lat

  ! !> convert latitude to normalized fixed-point integer
  ! ELEMENTAL FUNCTION fxp_lat(lat)
  !   REAL(wp), INTENT(in) :: lat
  !   INTEGER :: fxp_lat
  !   fxp_lat = NINT(lat_spread*normalized_lat(lat))
  ! END FUNCTION fxp_lat

  ! !> convert fixed-point latitude to floating-point value in range [-pi/2,pi/2]
  ! ELEMENTAL FUNCTION flp_lat(lat)
  !   INTEGER, INTENT(in) :: lat
  !   REAL(wp) :: flp_lat
  !   flp_lat = lat_contract*REAL(lat, wp)
  ! END FUNCTION flp_lat

  ! !> normalize latitude to range [-pi,pi]
  ! ELEMENTAL FUNCTION normalized_lon(lon) RESULT(normalized)
  !   REAL(wp), INTENT(in) :: lon
  !   REAL(wp) :: normalized, offset
  !   offset = SIGN(hcirc, lon)
  !   normalized = MOD(lon + offset, fcirc) - offset
  ! END FUNCTION normalized_lon

  ! !> convert longitude to normalized fixed-point integer
  ! ELEMENTAL FUNCTION fxp_lon(lon)
  !   REAL(wp), INTENT(in) :: lon
  !   INTEGER :: fxp_lon
  !   fxp_lon = NINT(lon_spread*normalized_lon(lon))
  ! END FUNCTION fxp_lon

  ! !> convert fixed-point longitude to floating-point value in range [-pi,pi]
  ! ELEMENTAL FUNCTION flp_lon(lon)
  !   INTEGER, INTENT(in) :: lon
  !   REAL(wp) :: flp_lon
  !   flp_lon = lon_contract*REAL(lon, wp)
  ! END FUNCTION flp_lon

  ! !-------------------------------------------------------------------------
  ! !> Merges a list of REAL(wp) numbers into another list, yielding the
  ! !  union set (without duplicates).
  ! !
  ! SUBROUTINE merge_values_into_set(nin_values, in_values, value_set, opt_tol)
  !   INTEGER, INTENT(IN)           :: nin_values !< no. of input values
  !   REAL(wp), INTENT(IN)           :: in_values(:) !< list of input values
  !   TYPE(t_value_set), INTENT(INOUT)        :: value_set !< result values (union set)
  !   REAL(wp), INTENT(IN), OPTIONAL :: opt_tol !< (optional:) tolerance for floating-point equality
  !   ! local variables
  !   CHARACTER(LEN=*), PARAMETER :: routine = '::merge_values_into_set'
  !   REAL(wp), PARAMETER :: default_tol = 1.e-12_wp
  !   REAL(wp) :: tol
  !   INTEGER  :: i, k, n0, nn0, ierrstat
  !   INTEGER, ALLOCATABLE :: idx(:)

  !   ! set tolerance for floating-point equality test
  !   tol = default_tol
  !   IF (PRESENT(opt_tol)) tol = opt_tol
  !   ! make sure that we have enough space in our "t_value_set" object:
  !   CALL resize_set(value_set, nin_values)

  !   n0 = value_set%nvalues
  !   nn0 = n0 + nin_values ! size of new set (with duplicates)
  !   ! first step: simply append the value to the end of the list:
  !   value_set%values((n0 + 1):nn0) = in_values(1:nin_values)
  !   ! second step: sort the whole set (including duplicates):
  !   CALL quicksort(value_set%values(1:nn0))
  !   ! revert levels largest-to-smallest if required
  !   IF (.NOT. value_set%sort_smallest_first) THEN
  !     value_set%values(1:nn0) = value_set%values(nn0:1:-1)
  !   END IF
  !   ! allocate temporary index buffer
  !   ALLOCATE (idx(nn0), STAT=ierrstat)
  !   IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE of temporary failed!")
  !   ! third step: loop over the sorted list and mark duplicates
  !   k = 0
  !   DO i = 1, nn0
  !     IF (k == 0) THEN
  !       k = k + 1
  !     ELSE
  !       IF (ABS(value_set%values(i - 1) - value_set%values(i)) > tol) k = k + 1
  !     END IF
  !     idx(i) = k
  !   END DO
  !   value_set%nvalues = k
  !   ! update value and index list (remove duplicates):
  !   DO i = 1, nn0
  !     value_set%values(idx(i)) = value_set%values(i)
  !   END DO
  !   ! clean up
  !   DEALLOCATE (idx, STAT=ierrstat)
  !   IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE of temporary failed!")
  ! END SUBROUTINE merge_values_into_set

  ! !-------------------------------------------------------------------------
  ! !> Compares a list of REAL(wp) numbers with the contents of a
  ! !  "t_value_set" data object. For a "value_set" with N values, this
  ! !  subroutine returns a LOGICAL array "out_flag" of size 1:N, where
  ! !  "out_flag(i) == .TRUE.", if the value "value_set%values(i)"
  ! !  matches one of the input numbers.
  ! !
  ! !  @todo For the sake of simplicity, this algorithm has quadratic
  ! !        complexity, though we could exploit the fact the the value
  ! !        set is sorted.
  ! !
  ! SUBROUTINE find_values_in_set(nin_values, in_values, value_set, out_flag, opt_tol)
  !   INTEGER, INTENT(IN)           :: nin_values !< no. of input values
  !   REAL(wp), INTENT(IN)           :: in_values(:) !< list of input values
  !   TYPE(t_value_set), INTENT(IN)           :: value_set !< set of values
  !   LOGICAL, INTENT(OUT)          :: out_flag(:) !<
  !   REAL(wp), INTENT(IN), OPTIONAL :: opt_tol !< (optional:) tolerance for floating-point equality
  !   ! local variables
  !   CHARACTER(LEN=*), PARAMETER :: routine = '::find_values_in_set'
  !   REAL(wp), PARAMETER :: default_tol = 1.e-12_wp
  !   REAL(wp) :: tol
  !   INTEGER  :: i, j

  !   ! set tolerance for floating-point equality test
  !   tol = default_tol
  !   IF (PRESENT(opt_tol)) tol = opt_tol
  !   ! make sure that our result array has the required size:
  !   IF (SIZE(out_flag) < value_set%nvalues) THEN
  !     CALL finish(routine, "Wrong dimension of output argument!")
  !   END IF
  !   out_flag(:) = .FALSE.
  !   DO i = 1, value_set%nvalues
  !     inner_loop: DO j = 1, nin_values
  !       IF (ABS(value_set%values(i) - in_values(j)) <= tol) THEN
  !         out_flag(i) = .TRUE.
  !         EXIT inner_loop
  !       END IF
  !     END DO inner_loop
  !   END DO

  ! END SUBROUTINE find_values_in_set

  ! !> Resizes the data structures of a "t_value_set" derived type object.
  ! !
  ! SUBROUTINE resize_set(value_set, nadditional_values)
  !   TYPE(t_value_set), INTENT(INOUT) :: value_set !< result values (union set)
  !   INTEGER, INTENT(IN)    :: nadditional_values
  !   ! local variables
  !   CHARACTER(LEN=*), PARAMETER :: routine = '::deallocate_set'
  !   INTEGER               :: ierrstat, oldsize, newsize
  !   REAL(wp), ALLOCATABLE :: tmp_rbuf(:)

  !   ! increase size of REAL(wp) buffer
  !   oldsize = 0
  !   IF (ALLOCATED(value_set%values)) THEN
  !     oldsize = SIZE(value_set%values)
  !     ALLOCATE (tmp_rbuf(oldsize), STAT=ierrstat)
  !     IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
  !     tmp_rbuf(1:oldsize) = value_set%values(1:oldsize)
  !     DEALLOCATE (value_set%values, STAT=ierrstat)
  !     IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  !   ELSE
  !     value_set%nvalues = 0
  !   END IF
  !   newsize = oldsize + nadditional_values
  !   ALLOCATE (value_set%values(newsize), STAT=ierrstat)
  !   IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
  !   IF (ALLOCATED(tmp_rbuf)) THEN
  !     value_set%values(1:oldsize) = tmp_rbuf(1:oldsize)
  !     ! clean up
  !     DEALLOCATE (tmp_rbuf, STAT=ierrstat)
  !     IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  !   END IF
  ! END SUBROUTINE resize_set

  ! !> Frees a "t_value_set" data type.
  ! !
  ! SUBROUTINE deallocate_set(value_set)
  !   TYPE(t_value_set), INTENT(INOUT) :: value_set !< result values (union set)
  !   ! local variables
  !   CHARACTER(LEN=*), PARAMETER :: routine = '::deallocate_set'
  !   INTEGER :: ierrstat

  !   value_set%nvalues = 0
  !   DEALLOCATE (value_set%values, STAT=ierrstat)
  !   IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  ! END SUBROUTINE deallocate_set

END MODULE mo_math_utilities
