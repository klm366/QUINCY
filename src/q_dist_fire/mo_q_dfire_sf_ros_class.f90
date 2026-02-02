!> QUINCY calculate fire rate of spread
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
!> @brief Main interface to calculate rate of spread and related variables (mean fire area, ...)
!>
!> Available implementations to compute Rate of Spread:
!> - Thonicke et al. (2010)
!> - Oberhagemann et al. (2025)
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_q_dfire_sf_ros_class
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish

  USE mo_q_dfire_config_class,  ONLY: Q_DFIRE_ROS_THONICKE, Q_DFIRE_ROS_OBERHAGEMANN
  USE mo_q_dfire_sf_ros_ober_class, ONLY: sf_calc_rate_of_spread_ober
  USE mo_q_dfire_sf_ros_thon_class, ONLY: sf_calc_rate_of_spread_thon

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_rate_of_spread

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_ros_class'


  !< Vegetation cover fraction
  TYPE t_sf_veget_cover_frac
    REAL(wp), ALLOCATABLE :: tree(:) !< How much of tree covers the tile [0 to 1]
    REAL(wp), ALLOCATABLE :: grass(:) !< How much grass covers the tile [0 to 1]
  END TYPE t_sf_veget_cover_frac

  INTERFACE t_sf_veget_cover_frac
    MODULE PROCEDURE :: create_veget_cover_frac
  END INTERFACE t_sf_veget_cover_frac


  !< Rate of spread characteristics type
  TYPE t_sf_ros
    INTEGER               :: nc
    REAL(wp), ALLOCATABLE :: ros_forward(:)        !< Forward rate of spread of the fire [m min-1]
    REAL(wp), ALLOCATABLE :: ros_backwards(:)      !< Backwards rate of spread of the fire [m min-1]
    REAL(wp), ALLOCATABLE :: distance_travelled(:) !< Total distance travelled [m]
    REAL(wp), ALLOCATABLE :: l2b_ratio(:)          !< length to breadth ratio [0 to 1]
    REAL(wp), ALLOCATABLE :: fire_duration(:)      !< [min]
    REAL(wp), ALLOCATABLE :: mean_fire_area(:)     !< [ha]

    REAL(wp), ALLOCATABLE :: reaction_intensity(:) !< I_R: the energy release rate per unit area of fire front [kJ m−2 min−1]
    REAL(wp), ALLOCATABLE :: prop_flux_ratio(:)    !< the proportion of IR that heats adjacent fuel particles to ignition [-]
    REAL(wp), ALLOCATABLE :: wind_factor(:)        !< wind_factor is a multiplier that accounts for the effect of wind in increasing the effective value of prop_flux_ratio[-]
    REAL(wp), ALLOCATABLE :: fbd(:)                !< Weighted average across PFT for fuel bulk density [kg m−3]
    REAL(wp), ALLOCATABLE :: eff_heating_number(:) !< the proportion of a fuel particle that is heated to ...
    !< the proportion of a fuel particle that is heated to [-]
    REAL(wp), ALLOCATABLE :: pre_ign_heat(:)       !< the amount of heat required to ignite a given mass of fuel [kJ kg−1]

    REAL(wp), ALLOCATABLE :: opt_react_velo(:)     !< Optimum reaction velocity [min-1]
    REAL(wp), ALLOCATABLE :: net_fuel_load(:)      !< Amount of fuel after substraction of its mineral content [kg m-2]
    REAL(wp), ALLOCATABLE :: moist_damp_coeff(:)   !< Moisture dampening coefficient [-]
    REAL(wp), ALLOCATABLE :: mineral_damp_coeff(:) !< Mineral dampening coefficient [-]
  END TYPE t_sf_ros

  INTERFACE t_sf_ros
    MODULE PROCEDURE :: create_sf_ros_init
  END INTERFACE t_sf_ros


CONTAINS

  FUNCTION calc_rate_of_spread(nc, nddlv, ndd, nlv, fuel_bulk_density, fwd_wind_factor, wind, moist_extinction, fdi, &
      veg_cover_frac, woi_dd, woi_dd, mf_dd, mf_lv, me_dd, me_lv, omega_o, mask) RESULT (ros)
    INTEGER, INTENT(in)                  :: nc
    INTEGER, INTENT(in)                  :: nddlv                       !< Dead fuel classes with livegrass in 1hr
    INTEGER, INTENT(in)                  :: ndd                         !< Dead fuel classes (1, 10, 100, dry livegrass)
    INTEGER, INTENT(in)                  :: nlv                         !< Live fuel classes (wet livegrass, woody)
    INTEGER, INTENT(in)                  :: ros_model                   !< ROS model implementation
    REAL(wp), INTENT(in)                 :: fuel_bulk_density(nc)       !< [kg m-3]
    REAL(wp), INTENT(in)                 :: fwd_wind_factor(nc)         !< [0 to 1]
    REAL(wp), INTENT(in)                 :: wind(nc)                    !< [m min-1]
    REAL(wp), INTENT(in)                 :: moist_extinction(nc)        !< [0 to 1]
    REAL(wp), INTENT(in)                 :: fdi(nc)                     !< Fire dange index [0 to 1]
    TYPE(t_sf_veget_cover_frac), INTENT(in) :: veg_cover_frac           !< [0 to 1]
    REAL(wp), INTENT(in)                 :: woi_dd_lvg(nc, nddlv)       !< Dead fuel as live grass in 1 hour (1hr+dry livegrass, 10, 100 hours)
    REAL(wp), INTENT(in)                 :: woi_dd(nc, ndd)             !< Dead fuel (1hr, 10 hors, 100hrs, dry livegress) dead fuel [kg DM m-2], note in kg
    REAL(wp), INTENT(in)                 :: woi_lv(nc, nlv)             !< Live fuel (wet livegrass, woody) live fuel [kg DM m-2], note in kg
    REAL(wp), INTENT(in)                 :: mf_dd(nc), mf_lv(nc)        !< Fuel moisture [0 to 1] (#nc)
    REAL(wp), INTENT(in)                 :: me_dd(nc), me_lv(nc)        !< Moisture of extinction [0 to 1] (#nc)
    LOGICAL, INTENT(in)                  :: mask(nc)                    !< True =  compute gridcell(s)
    TYPE(t_sf_ros)                       :: ros                         !< output [-]

    IF (ros_model == Q_DFIRE_ROS_THONICKE) THEN
      !> Fire rate of spread (IR, optimum reaction velocity, fuel bulk density, net fuel loading, ...)
      !> calculate the rate of spread following Thonicke et al. (2010)
      ros = sf_calc_rate_of_spread_thon(nc, nddlv, &
                                      lctlib_fuel_bulk_density, lctlib_fwd_wind_factor, &
                                      wind_sf, me_dd, fdi, veg_cover_frac, &
                                      woi_dd_lvg, &
                                      mf_dd, &
                                      mask=has_fire_conditions)
    ELSE IF (ros_model == Q_DFIRE_ROS_OBERHAGEMANN) THEN

      ros = sf_calc_rate_of_spread_ober(&
          nc, df_sf_sigma_dd_nr, df_sf_sigma_lv_nr, &
          df_sigmas_dead, df_sigmas_live, &
          woi_dd, woi_lv, &
          mf_dd, mf_lv, & ! todo
          me_dd, me_lv, & ! todo
          lctlib_fuel_bulk_density, & ! todo: rename to fbd
          lctlib_fwd_wind_factor, wind_sf, & !< todo: rename to fwd_wind_factor
          mask=has_fire_conditions)

    ELSE
      CALL finish("calc_spitfire", "Unexpected Rate of spread model found ('thonicke' or 'oberhagemann')")
    END IF

    !< ( Eq. 10 from Thonicke et al. (2010))
    ros%ros_backwards(inc) = ros_forward(inc) * exp(-0.012_wp * u_f_wind(inc)) !< [m min-1]

    !< (Eqs. A38, A39, A40 from Pfeiffer et al. (2013))
    ros%l2b_ratio(inc) = calc_length_to_breadth_ratio(u_f_wind(inc), &
        veg_cover_frac%tree(inc), &
        veg_cover_frac%grass(inc))

    ros%fire_duration(inc) = 241._wp / (1._wp + 240._wp * exp(-11.06_wp * fdi(inc))) !< Eq. 14 from Thonicke et al. (2010) [min]
    IF (fire_duration(inc) > df_sf_fire_duration_thd) THEN
      fire_duration(inc) = df_sf_fire_duration_thd
    END IF
    ros%distance_travelled(inc) = (ros_forward(inc) + ros_backwards(inc)) * &
        fire_duration(inc) !< (Eq. A42 from Pfeiffer et al. (2013)) [m]
    !< mean fire area (Eq. 11 from Thonicke et al. (2010))
    ros%mean_fire_area(inc) = calc_mean_fire_area(l2b_ratio(inc), distance_travelled(inc)) !< [ha]

  END FUNCTION calc_rate_of_spread


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
  !> Calculate wind factor (Eq. A5 from Thonicke et al. (2010))
  !>
  ELEMENTAL PURE FUNCTION calc_wind_factor(u_f_wind, sigma, beta, beta_op) RESULT (wind_factor)
    REAL(wp), INTENT(in)        :: u_f_wind     !< Forward wind [m min-1]
    REAL(wp), INTENT(in)        :: sigma        !< Surface area to volume ratio of the fuel [cm-1]
    REAL(wp), INTENT(in)        :: beta         !< Packing ratio [-]
    REAL(wp), INTENT(in)        :: beta_op      !< [-]
    REAL(wp)                    :: wind_factor  !< Output [-]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                    :: C_var, B_var, E_var
    REAL(wp)                    :: beta_ratio
    ! ----------------------------------------------------------------------------------------------------- !
    C_var = 7.47_wp * exp(-0.8711_wp * (sigma**0.55_wp)) !< Eq. A8
    B_var = 0.15988_wp * (sigma ** 0.54_wp) !< Eq. A7
    E_var = 0.7515_wp * exp(-0.01094_wp * sigma) !< Eq. A9
    beta_ratio = beta / beta_op
    wind_factor = C_var * ((3.281_wp * u_f_wind) ** B_var) * (beta_ratio ** (-E_var))
  END FUNCTION calc_wind_factor


  ! ======================================================================================================= !
  !>
  !> Calculate propagation flux ratio (Eq. A2 from Thonicke et al. (2010))
  !>
  ELEMENTAL PURE FUNCTION calc_prop_flux_ratio(sigma, beta) RESULT (prop_flux_ratio)
    REAL(wp), INTENT(in) :: sigma           !< surface-area to volume ratio of the fuel [cm-1]
    REAL(wp), INTENT(in) :: beta            !< packing ratio [-]
    REAL(wp)             :: prop_flux_ratio !< output [-]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: e_calc
    ! ----------------------------------------------------------------------------------------------------- !
    e_calc = (0.792_wp + 3.7597_wp * sqrt(sigma)) * (beta + 0.1_wp)
    prop_flux_ratio = exp(e_calc) / (192._wp + 7.9095_wp * sigma)
  END FUNCTION calc_prop_flux_ratio


  ! ======================================================================================================= !
  !>
  !> Calculate The optimum reaction velocity (Table A1, 1st from Thonicke et al. (2010))
  !>
  ELEMENTAL PURE FUNCTION calc_opt_react_velocity(max_react_vel, beta, beta_op, A_var) RESULT (opt_react_vel)
    REAL(wp), INTENT(in)        :: max_react_vel   !< Maximum reaction velocity [min-1]
    REAL(wp), INTENT(in)        :: beta            !< packing ratio [kg m-3]
    REAL(wp), INTENT(in)        :: beta_op         !< Eq. A6 [-]
    REAL(wp), INTENT(in)        :: A_var           !< [-]
    REAL(wp)                    :: opt_react_vel   !< output [min-1]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: beta_ratio
    ! ----------------------------------------------------------------------------------------------------- !
    beta_ratio = beta / beta_op
    opt_react_vel = max_react_vel * (beta_ratio ** A_var) * exp(A_var * (1 - beta_ratio))
  END FUNCTION calc_opt_react_velocity


  ! ======================================================================================================= !
  !>
  !> Calculate the weighted average of the Fuel bulk density for the current PFT.
  !> Using the Pfeiffer et al. (2013) paper, equation A17. Remaining equations should be implemented
  !> when accounting for the other pfts.
  !> Thonicke's et al. (2010) is not used as the equation does not seems correct when accounting for
  !> 10 and 100 hours fuel bulk densities (slightly more).
  !>
  FUNCTION calc_fbd(nfclasses, fuel_bulk_density, dead_fuel) RESULT (fbd)
    INTEGER, INTENT(in)         :: nfclasses !< The number of used fuel classes
    REAL(wp), INTENT(in)        :: fuel_bulk_density !< [kg m-3]
    REAL(wp), INTENT(in)        :: dead_fuel(nfclasses) !< 1, 10, 100 hours [gDM m-2]
    REAL(wp)                    :: fbd !< output [kg m-3]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                    :: total_dead_fuel
    REAL(wp)                    :: fbd_accounting_fc
    ! ----------------------------------------------------------------------------------------------------- !
    total_dead_fuel = sum(dead_fuel)
    fbd_accounting_fc = fuel_bulk_density * SUM(dead_fuel * df_sf_fuel_class_factors)
    fbd = fbd_accounting_fc / total_dead_fuel
  END FUNCTION calc_fbd
#endif
END MODULE mo_q_dfire_sf_ros_class
