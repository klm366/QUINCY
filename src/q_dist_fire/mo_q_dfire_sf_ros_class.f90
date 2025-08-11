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
!> @brief Module to determine the fire Rate of spread
!>    This involves equations from:
!>    - Section 2.2.4, Equations 9 to 14
!>    - Table A1
!>    - Appendix A (equations A1 to A9)
!>
!>  Most of the equations are based on Thonicke et al. (2010). Nonetheless, some equations are also
!>  using Pfeiffer et al. (2013). In all cases, It is noted above each equation.
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
  !< Total number of sigmas (1,10,100) hours,  Pfeiffer el al. (2013), Eq. A18 - Used to calculate fuel bulk density
  USE mo_q_dfire_constants,     ONLY: df_sf_sigma_nr, df_sf_fuel_class_factors
  !<
  USE mo_q_dfire_constants,     ONLY: df_sf_fire_duration_thd
  !<  Rate of spread constants
  USE mo_q_dfire_constants,     ONLY: df_sf_oven_dry_part_dens, df_sigmas, &
                                      df_sf_mineral_fract, df_sf_mineral_damp_coeff, &
                                      df_sf_fuel_heat_content
  USE mo_q_dfire_sf_adaptor,    ONLY: sf_convert_units

  USE mo_quincy_stats,          ONLY: stats_average_weight
  USE mo_jsb_math_constants,    ONLY: pi

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_sf_ros, t_sf_veget_cover_frac, sf_calc_rate_of_spread

  PUBLIC :: calc_fbd, calc_sigma, &
            calc_opt_react_velocity, calc_moist_damp_coeff, &
            calc_prop_flux_ratio, calc_wind_factor, &
            calc_length_to_breadth_ratio, calc_mean_fire_area

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

  ! ======================================================================================================= !
  !>
  !> calc
  !>
  !> Do the science
  !>
  FUNCTION sf_calc_rate_of_spread(nc, nfuels, fuel_bulk_density, fwd_wind_factor, wind, moist_extinction, fdi, &
                                        veg_cover_frac, w_oi, omega_o, mask) RESULT(ros)
    INTEGER, INTENT(in)                  :: nc
    INTEGER, INTENT(in)                  :: nfuels                      !< 1, 10 and 100 hours fuel classes
    REAL(wp), INTENT(in)                 :: fuel_bulk_density(nc)       !< [kg m-3]
    REAL(wp), INTENT(in)                 :: fwd_wind_factor(nc)         !< [0 to 1]
    REAL(wp), INTENT(in)                 :: wind(nc)                    !< [m min-1]
    REAL(wp), INTENT(in)                 :: moist_extinction(nc)        !< [0 to 1]
    REAL(wp), INTENT(in)                 :: fdi(nc)                     !< Fire dange index [0 to 1]
    TYPE(t_sf_veget_cover_frac), INTENT(in) :: veg_cover_frac           !< [0 to 1]
    REAL(wp), INTENT(in)                 :: w_oi(nc, nfuels)            !< [g DM m-2]
    !< According to Pfeiffer et al. (2013), omega_n must be treated as omega_o
    REAL(wp), INTENT(in)                 :: omega_o(nc)                 !< [0 to 1] (#nc)
    LOGICAL, INTENT(in)                  :: mask(nc)                    !< True =  compute gridcell(s)
    TYPE(t_sf_ros)                       :: ros                         !< output [-]

    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: beta(nc), beta_op(nc)
    REAL(wp) :: sigma(nc)                 !< Surface-area-to-volume-ratio [cm-1]
    REAL(wp) :: A_var(nc)
    REAL(wp) :: max_reaction_velocity(nc)
    REAL(wp) :: u_f_wind(nc)
    REAL(wp) :: cover_fract_pft_ncs(nc) !< Foliar projective cover

    INTEGER :: inc
    ! ----------------------------------------------------------------------------------------------------- !
    ros = t_sf_ros(nc)
    cover_fract_pft_ncs(:) = 1._wp !TODO: properly calculate the fpc for IQ
    ros%mean_fire_area = 0._wp

    DO inc = 1, nc
      IF (mask(inc)) THEN
        !< Calculate wind speed (paragraph before eq. 10 from Thonicke et al. (2010))
        u_f_wind(inc) = fwd_wind_factor(inc) * cover_fract_pft_ncs(inc) * wind(inc) !< [m min-1]

        !<  (Table A1, 5th from Thonicke et al. (2010))
        ros%net_fuel_load(inc) = (1._wp - df_sf_mineral_fract) * SUM(w_oi(inc,:))
        ros%net_fuel_load(inc) = sf_convert_units("g", "kg", ros%net_fuel_load(inc))
        ros%fbd(inc) = fuel_bulk_density(inc) !< fuel bulk density [kg m-3], precomputed at mo_dfire_fuel.f90 + aggregation to         average

        sigma(inc) = calc_sigma(df_sf_sigma_nr, nfuels, df_sigmas, w_oi(inc,:)) !< Table A1, 3rd equation [cm-1] from Thonicke et al. (2010)
        !< (Table A1, 3rd from Thonicke et al. (2010))
        beta(inc) =  ros%fbd(inc) / df_sf_oven_dry_part_dens !< [-]
        !< Calculate beta_op (Eq. A6 from Thonicke et al. (2010))
        beta_op(inc) = 0.200395_wp * (sigma(inc) ** (-0.8189_wp)) !< [-]

        A_var(inc) = 8.9033_wp * (sigma(inc) ** (-0.7913)) !< Table A1, 4th equation from Thonicke et al. (2010)
        max_reaction_velocity(inc) = 1._wp / (0.0591_wp  + 2.926_wp * sigma(inc) ** (-1.5_wp)) !< Table A1, 2nd equation from Thonicke et al. (2010) [min-1]
        ros%opt_react_velo(inc) = calc_opt_react_velocity(max_reaction_velocity(inc), beta(inc), &
            beta_op(inc), A_var(inc)) ! Table A1, 1st equation from Thonicke et al. (2010) [min-1]

        ros%moist_damp_coeff(inc) = calc_moist_damp_coeff(omega_o(inc), moist_extinction(inc)) !< Table A1, 6th equation from Thonicke et al. (2010)
        ros%mineral_damp_coeff(inc) = df_sf_mineral_damp_coeff

        !< (Eq. A1 from Thonicke et al. (2010))
        ros%reaction_intensity(inc) = ros%opt_react_velo(inc) * &
                                      ros%net_fuel_load(inc) * df_sf_fuel_heat_content * &
                                      ros%moist_damp_coeff(inc) * df_sf_mineral_damp_coeff !< [kJ m−2 min−1]

        ros%prop_flux_ratio(inc) = calc_prop_flux_ratio(sigma(inc), beta(inc)) !< Eq. A2 from Thonicke et al. (2010) [-]
        ros%wind_factor(inc) = calc_wind_factor(u_f_wind(inc), sigma(inc), beta(inc), beta_op(inc)) !< Eq. A5 from Thonicke et al. (2010) [-]
        ros%eff_heating_number(inc) = exp(-4.528_wp / sigma(inc)) !< (Eq. A3 from Thonicke et al. (2010)) [-]
        ros%pre_ign_heat(inc) = 581._wp + 2594._wp * omega_o(inc) !< (Eq. A4 from Thonicke et al. (2010)) [kJ kg-1]

        !< ( Eq. 9 from Thonicke et al. (2010))
        ros%ros_forward(inc) = ( ros%reaction_intensity(inc)  * ros%prop_flux_ratio(inc) * &
                                (1 + ros%wind_factor(inc)) ) / &
                                ( ros%fbd(inc) * ros%eff_heating_number(inc) * ros%pre_ign_heat(inc) ) !< [m min-1]
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
      END IF ! IF (mask(inc)) THEN
    ENDDO ! DO inc = 1, nc

  END FUNCTION sf_calc_rate_of_spread


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


  ! ======================================================================================================= !
  !>
  !> Initialize t_sf_ros to 0
  !>
  FUNCTION create_sf_ros_init(nc) RESULT(sf_ros)
    INTEGER, INTENT(in)                  :: nc        !< Total number of gridcells
    TYPE(t_sf_ros)                       :: sf_ros    !< output [-]
    ! ----------------------------------------------------------------------------------------------------- !
    sf_ros%nc = nc
    ALLOCATE(sf_ros%ros_forward(nc))
    ALLOCATE(sf_ros%ros_backwards(nc))
    ALLOCATE(sf_ros%distance_travelled(nc))
    ALLOCATE(sf_ros%l2b_ratio(nc))
    ALLOCATE(sf_ros%fire_duration(nc))
    ALLOCATE(sf_ros%mean_fire_area(nc))

    ALLOCATE(sf_ros%reaction_intensity(nc))
    ALLOCATE(sf_ros%prop_flux_ratio(nc))
    ALLOCATE(sf_ros%wind_factor(nc))
    ALLOCATE(sf_ros%fbd(nc))
    ALLOCATE(sf_ros%eff_heating_number(nc))
    ALLOCATE(sf_ros%pre_ign_heat(nc))

    ALLOCATE(sf_ros%opt_react_velo(nc))
    ALLOCATE(sf_ros%net_fuel_load(nc))
    ALLOCATE(sf_ros%moist_damp_coeff(nc))
    ALLOCATE(sf_ros%mineral_damp_coeff(nc))

    sf_ros%ros_forward(:) = 0._wp
    sf_ros%ros_backwards(:) = 0._wp
    sf_ros%distance_travelled(:) = 0._wp
    sf_ros%l2b_ratio(:) = 0._wp
    sf_ros%fire_duration(:) = 0._wp
    sf_ros%mean_fire_area = 0._wp

    sf_ros%reaction_intensity(:) = 0._wp
    sf_ros%prop_flux_ratio(:) = 0._wp
    sf_ros%wind_factor(:) = 0._wp
    sf_ros%fbd(:) = 0._wp
    sf_ros%eff_heating_number(:) = 0._wp
    sf_ros%pre_ign_heat(:) = 0._wp

    sf_ros%opt_react_velo(:) = 0._wp
    sf_ros%net_fuel_load(:) = 0._wp
    sf_ros%moist_damp_coeff(:) = 0._wp
    sf_ros%mineral_damp_coeff(:) = 0._wp
  END FUNCTION create_sf_ros_init


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
  !> Calculate the moisture dampening coefficient (Table A1, 6th from Thonicke et al. (2010))
  !>
  !> when equal to 1, the fuel is moisture-free
  !> Todo: apply Pfeiffer approach on weighting moisture of extinction (A25)
  !>
  ELEMENTAL PURE FUNCTION calc_moist_damp_coeff(omega_o, moist_extinction) RESULT (moist_damp_coeff)
    REAL(wp), INTENT(in) :: omega_o           !< Fuel moisture [0 to 1]
    REAL(wp), INTENT(in) :: moist_extinction  !< Moisture of extinction [0 to 1]
    REAL(wp)             :: moist_damp_coeff  !< output [-]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: ratio
    ! ----------------------------------------------------------------------------------------------------- !
    ratio = omega_o / moist_extinction
    moist_damp_coeff = 1._wp - (2.59_wp * ratio) + (5.11_wp * (ratio ** 2_wp)) - (3.52_wp * (ratio ** 3_wp))
  END FUNCTION calc_moist_damp_coeff


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
  !> Calculate sigma (Table A1, 2nd form Thonicke et al. (2010))
  !>
  FUNCTION calc_sigma(nsigma, nfuels, sigmas, woi) RESULT (sigma)
    INTEGER, INTENT(in)  :: nsigma, nfuels                    !<
    REAL(wp), INTENT(in) :: sigmas(nsigma)                    !<
    REAL(wp), INTENT(in) :: woi(nfuels)             !< dead fuel [kg m-2]
    REAL(wp)             :: sigma                             !< output [cm-1]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)             :: total_woi
    ! ----------------------------------------------------------------------------------------------------- !
    total_woi = sum(woi(:))
    IF (total_woi < EPSILON(0._wp)) THEN ! woi == 0.
      sigma = 0.0001_wp
    ELSE
      sigma = stats_average_weight(nsigma, sigmas(:), woi(:))
    END IF
  END FUNCTION calc_sigma


  ! ======================================================================================================= !
  !>
  !> Calculate the weighted average of the Fuel bulk density for the current PFT.
  !> Using the Pfeiffer et al. (2013) paper, equation A17. Remaining equations should be implemented
  !> when accounting for the other pfts.
  !> Thonicke's et al. (2010) is not used as the equation does not seems correct when accounting for 10 and 100 hours
  !> fuel bulk densities (slightly more).
  !>
  !> TODO: currently done for 1 PFT. Later on, it will be changed for each PFT.
  !>
  FUNCTION calc_fbd(nfclasses, fuel_bulk_density, dead_fuel) RESULT (fbd)
    INTEGER, INTENT(in)         :: nfclasses !< The number of used fuel classes
    REAL(wp), INTENT(in)        :: fuel_bulk_density !< [kg m-3]
    REAL(wp), INTENT(in)        :: dead_fuel(nfclasses) !< 1, 10, 100 hours [gDM m-2]
    REAL(wp)                    :: fbd !< output
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
