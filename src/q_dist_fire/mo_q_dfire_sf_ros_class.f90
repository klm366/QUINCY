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
!>
!> Available implementations to compute Rate of Spread:
!> - Thonicke et al. (2010)
!> - Oberhagemann et al. (2025)
!>
MODULE mo_q_dfire_sf_ros_class
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish

  USE mo_q_dfire_sf_classes,        ONLY: t_sf_ros
  USE mo_q_dfire_config_class,      ONLY: Q_DFIRE_ROS_THONICKE, Q_DFIRE_ROS_OBERHAGEMANN
  USE mo_q_dfire_sf_ros_ober_class, ONLY: sf_calc_rate_of_spread_ober, calc_live_me
  USE mo_q_dfire_sf_ros_thon_class, ONLY: sf_calc_rate_of_spread_thon

  USE mo_q_dfire_constants,         ONLY: df_sf_fire_duration_thd
  USE mo_q_dfire_constants,         ONLY: df_sf_sigma_dd_nr, df_sf_sigma_lv_nr
  USE mo_q_dfire_constants,         ONLY: df_sigmas_dead, df_sigmas_live

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sf_calc_rate_of_spread, t_sf_veget_cover_frac, calc_mean_fire_area

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_ros_class'


  !< Vegetation cover fraction
  TYPE t_sf_veget_cover_frac
    REAL(wp), ALLOCATABLE :: tree(:) !< How much of tree covers the tile [0 to 1]
    REAL(wp), ALLOCATABLE :: grass(:) !< How much grass covers the tile [0 to 1]
  END TYPE t_sf_veget_cover_frac

  INTERFACE t_sf_veget_cover_frac
    MODULE PROCEDURE :: create_veget_cover_frac
  END INTERFACE t_sf_veget_cover_frac


CONTAINS

  !>
  !> Calculate rate of spread common all sub fire models
  !>
  !>
  SUBROUTINE sf_calc_rate_of_spread(nc, &
                  fwd_wind_factor, wind, &
                  fdi, &
                  veg_cover_frac, &
                  mask, ros)
    INTEGER, INTENT(in)                  :: nc
    REAL(wp), INTENT(in)                 :: fwd_wind_factor(nc)         !< [0 to 1]
    REAL(wp), INTENT(in)                 :: wind(nc)                    !< [m min-1]
    REAL(wp), INTENT(in)                 :: fdi(nc)                     !< Fire dange index [0 to 1]
    TYPE(t_sf_veget_cover_frac), INTENT(in) :: veg_cover_frac           !< Vegetation cover fraction [0 to 1]
    LOGICAL, INTENT(in)                  :: mask(nc)                    !< True =  compute gridcell(s)
    TYPE(t_sf_ros), INTENT(inout)        :: ros                         !< output [-]

    REAL(wp) :: u_f_wind(nc)  !< Wind speed at vegetation cover [m min-1]
    INTEGER  :: inc

    !< Calculate wind speed (paragraph before eq. 10 from Thonicke et al. (2010))
    u_f_wind = fwd_wind_factor * wind

    DO inc=1, nc
      IF (.NOT. mask(inc)) CYCLE !< skip gridcell

      !< ( Eq. 10 from Thonicke et al. (2010))
      ros%ros_backwards(inc) = ros%ros_forward(inc) * exp(-0.012_wp * u_f_wind(inc)) !< [m min-1]

      !< (Eqs. A38, A39, A40 from Pfeiffer et al. (2013))
      ros%l2b_ratio(inc) = calc_length_to_breadth_ratio(u_f_wind(inc), &
          veg_cover_frac%tree(inc), &
          veg_cover_frac%grass(inc))

      ros%fire_duration(inc) = 241._wp / (1._wp + 240._wp * exp(-11.06_wp * fdi(inc))) !< Eq. 14 from Thonicke et al. (2010) [min]
      IF (ros%fire_duration(inc) > df_sf_fire_duration_thd) THEN
        ros%fire_duration(inc) = df_sf_fire_duration_thd
      END IF
      ros%distance_travelled(inc) = (ros%ros_forward(inc) + ros%ros_backwards(inc)) * &
          ros%fire_duration(inc) !< (Eq. A42 from Pfeiffer et al. (2013)) [m]
      !< mean fire area (Eq. 11 from Thonicke et al. (2010))
      ros%mean_fire_area(inc) = calc_mean_fire_area(ros%l2b_ratio(inc), ros%distance_travelled(inc)) !< [ha]
    ENDDO

  END SUBROUTINE sf_calc_rate_of_spread


  ! ======================================================================================================= !
  !>
  !> Calculate length to breadth grass ratio
  !> Eq. A39 from Pfeiffer et al. (2013)
  !>
  PURE ELEMENTAL FUNCTION calc_l2b_grass_ratio(u_f_wind) RESULT (l2b_ratio)
    REAL(wp), INTENT(in)        :: u_f_wind   !< [m min-1]
    REAL(wp)                    :: l2b_ratio  !< output [-]
    l2b_ratio = 1.1 + (u_f_wind ** 0.464)
  END FUNCTION calc_l2b_grass_ratio


  ! ======================================================================================================= !
  !>
  !> Calculate length to breadth ratio (Eqs. A38, A39, A40 from Pfeiffer et al. (2013))
  !> Eqs A38 and A39 used 0.06 is the conversion factor from km/h to m/min,
  !> It is not included as units are already in m/min
  !>
  ELEMENTAL PURE FUNCTION calc_length_to_breadth_ratio(u_f_wind, tree_cover, grass_cover) RESULT (l2b_ratio)
    USE mo_quincy_stats,          ONLY: stats_average_weight

    REAL(wp), INTENT(in) :: u_f_wind    !< [m min-1]
    REAL(wp), INTENT(in) :: tree_cover  !< [0 to 1]
    REAL(wp), INTENT(in) :: grass_cover !< [0 to 1]
    REAL(wp)             :: l2b_ratio   !< output [-]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)            :: lb_tree, lb_grass
    INTEGER, PARAMETER  :: nr_elems = 2
    REAL(wp), PARAMETER :: wind_speed_thd = 16.67 !< Threshold for low wind speeds proposed by Pfeiffer et al. (2013) after Eq. A40 [m min-1] == 1 km / h
    REAL(wp), PARAMETER :: max_l2b_ratio_thd = 8. !< Threshold for maximum ratio of a fire proposed by (After eq. A40, from Pfeiffer et al. (2013)) [-]
    LOGICAL             :: is_wind_poor, has_litter_only
    ! -----------------------------------------------------------------------------this%------------------------ !
    is_wind_poor = (u_f_wind < wind_speed_thd)
    has_litter_only = (tree_cover + grass_cover < 1.e-6_wp) !< Close to 0
    IF (is_wind_poor) THEN
      l2b_ratio = 1.0
    ELSE IF (has_litter_only) THEN !< Litter behaves as Grass (Pfeiffer's et al. (2013) source code)
      l2b_ratio = calc_l2b_grass_ratio(u_f_wind)
    ELSE
      lb_tree = 1.0_wp + 8.729_wp * ((1._wp - exp(-0.03_wp * u_f_wind))**2.155_wp) !< (Eq. A38)
      lb_grass =  calc_l2b_grass_ratio(u_f_wind)
      l2b_ratio = stats_average_weight(nr_elems, [lb_tree, lb_grass], [tree_cover, grass_cover])
    END IF

    l2b_ratio = min(l2b_ratio, max_l2b_ratio_thd)
  END FUNCTION calc_length_to_breadth_ratio


  ! ======================================================================================================= !
  !>
  !> Calculate the mean fire area (Eq. 11 from Thonicke et al. (2010) )
  !>
  ELEMENTAL PURE FUNCTION calc_mean_fire_area(l2b_ratio, distance_travelled) RESULT (mean_fire_area)
    USE mo_jsb_math_constants,    ONLY: pi

    REAL(wp), INTENT(in)        :: l2b_ratio          !< length to breadth ratio [-]
    REAL(wp), INTENT(in)        :: distance_travelled !< Total distance travelled by fire [m]
    REAL(wp)                    :: mean_fire_area     !< Output [ha]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: lb_term
    ! ----------------------------------------------------------------------------------------------------- !
    lb_term = pi / (4.0_wp * l2b_ratio)
    mean_fire_area = lb_term * (distance_travelled ** 2) / 10000._wp
  END FUNCTION calc_mean_fire_area


  ! ======================================================================================================= !
  !>
  !> t_sf_veget_cover_frac constructor
  !>
  FUNCTION create_veget_cover_frac(nc, tree_cover_frac, grass_cover_frac) RESULT  (veg_cover_frac)
    INTEGER, INTENT(in) :: nc
    REAL(wp), INTENT(in) :: tree_cover_frac(nc)   !< [0 to 1]
    REAL(wp), INTENT(in) :: grass_cover_frac(nc)  !< [0 to 1]
    TYPE(t_sf_veget_cover_frac) :: veg_cover_frac

    ALLOCATE(veg_cover_frac%tree, source = tree_cover_frac)
    ALLOCATE(veg_cover_frac%grass, source = grass_cover_frac)
  END FUNCTION create_veget_cover_frac


#endif
END MODULE mo_q_dfire_sf_ros_class
