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
!> @brief Rate of spread based on Oberhagemann et al. (2025) and Thonicke et al. (2010)
!>
!>  Oberhagemann et al. (2025) applies weighting factors on the fuel classes.
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_q_dfire_sf_ros_ober_class
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

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_ros_ober_class'

  INTEGER, PARAMETER :: ix_dead = 1
  INTEGER, PARAMETER :: ix_live = 2

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
    REAL(wp), ALLOCATABLE :: net_fuel_load(:,:)    !< Amount of fuel after substraction of its mineral content [kg m-2] (#dead/live, #nc)
    REAL(wp), ALLOCATABLE :: moist_damp_coeff(:)   !< Moisture dampening coefficient [-]
    REAL(wp), ALLOCATABLE :: mineral_damp_coeff(:) !< Mineral dampening coefficient [-]
  END TYPE t_sf_ros

  INTERFACE t_sf_ros
    MODULE PROCEDURE :: create_sf_ros_init
  END INTERFACE t_sf_ros

  !< weight factor variable
  TYPE :: t_fwvar_r1d
    REAL(wp), ALLOCATABLE :: d(:) !< dead (#nc)
    REAL(wp), ALLOCATABLE :: l(:) !< live (#nc)
    REAL(wp), ALLOCATABLE :: total(:) !< total, or dead + live (#nc)
  END TYPE t_fwvar_r1d

  INTERFACE t_fwvar_r1d
    MODULE PROCEDURE :: create_fwvar_r1d
  END INTERFACE t_fwvar_r1d

  !< weight factor variable
  TYPE :: t_fwvar_r2d
    REAL(wp), ALLOCATABLE :: d(:,:) !< dead (#nc)
    REAL(wp), ALLOCATABLE :: l(:,:) !< live (#nc)
    REAL(wp), ALLOCATABLE :: total(:,:) !< total, or dead + live (#nc)
  END TYPE t_fwvar_r2d

  INTERFACE t_fwvar_r2d
    MODULE PROCEDURE :: create_fwvar_r2d
  END INTERFACE t_fwvar_r2d


  TYPE :: t_sf_weighting_factors
    TYPE(t_fwvar_r2d)     :: Aij    !< (#fuel classes-dead, live-, #nc)
    TYPE(t_fwvar_r1d)     :: Ai     !< (#nc)
    REAL(wp), ALLOCATABLE :: At(:)  !< A_dead + A_live (#nc)
    TYPE(t_fwvar_r2d)     :: fij    !< (#fuel classes-dead, live-, #nc)
    TYPE(t_fwvar_r1d)     :: fi     !< (#nc)
    TYPE(t_fwvar_r2d)     :: gij    !< (#nfuelclasses-dead, live-, #nc)
  END TYPE t_sf_weighting_factors

  INTERFACE t_sf_weighting_factors
    MODULE PROCEDURE :: create_sf_weighting_factors
  END INTERFACE t_sf_weighting_factors

CONTAINS

  !>
  !> Constructor for weighting factors derived type
  !>
  FUNCTION create_sf_weighting_factors(nc, ndeadfc, nlivefc) RESULT (wfactors)
    INTEGER, INTENT(in) :: nc
    INTEGER, INTENT(in) :: ndeadfc !< 1 hr, 10 hr, ..., live grass, ...
    INTEGER, INTENT(in) :: nlivefc !< 1 hr, 10 hr, ..., live grass, ...
    TYPE(t_sf_weighting_factors) :: wfactors

    ALLOCATE(wfactors % At(nc))

    wfactors % Aij = t_fwvar_r2d(nc, ndeadfc, nlivefc)
    wfactors % Ai = t_fwvar_r1d(nc)
    wfactors % At = 0._wp

    wfactors % fij = t_fwvar_r2d(nc, ndeadfc, nlivefc)
    wfactors % fi = t_fwvar_r1d(nc)

    wfactors % gij = t_fwvar_r2d(nc, ndeadfc, nlivefc)
  END FUNCTION create_sf_weighting_factors


  ! ======================================================================================================= !
  !>
  !> Constructor for t_fwvar
  !>
  FUNCTION create_fwvar_r1d(nc, has_total) RESULT (var)
    INTEGER, INTENT(in) :: nc
    LOGICAL, INTENT(in), OPTIONAL :: has_total
    TYPE(t_fwvar), INTENT(in) :: var
    ALLOCATE(var%l(nc))
    ALLOCATE(var%d(nc))

    var%l = 0._wp
    var%d = 0._wp

    IF (PRESENT(has_total) .AND. has_total) THEN
      ALLOCATE(var%t(nc))
      var%t = 0._wp
    END IF
  END FUNCTION create_fwvar_r1d


  ! ======================================================================================================= !
  !>
  !> Constructor for t_fwvar_r2d
  !>
  FUNCTION create_fwvar_r2d(nc, nddfc, nlvfc) RESULT (var)
    INTEGER, INTENT(in) :: nc
    INTEGER, INTENT(in) :: nddfc !< n dead fuel classes
    INTEGER, INTENT(in) :: nlvfc !< n live fuel classes
    TYPE(t_fwvar_r2d), INTENT(in) :: var
    ALLOCATE(var%l(nddfc, nc))
    ALLOCATE(var%d(nlvfc, nc))
    !    ALLOCATE(var%t(nc))
    var%l = 0._wp
    var%d = 0._wp
    !    var%t = 0._wp
  END FUNCTION create_fwvar_r2d


  ! ======================================================================================================= !
  !>
  !> calc
  !>
  !> Do the science
  !>
  FUNCTION sf_calc_rate_of_spread(nc, ndd, nlv, woi_dd, woi_lv, sigmas_dd, sigmas_lv, mask) RESULT(ros)
    INTEGER, INTENT(in) :: nc
    INTEGER, INTENT(in) :: ndd, nlv  !< number of dead fuel classes, number of live fuel classes
    REAL(wp), INTENT(in) :: woi_dd(ndd), woi_lv(nlv) !< dead/live fuel load [kg DM m-1] -> notes units in kg (instead of g)
    REAL(wp), INTENT(in) :: sigmas_dd(ndd), sigmas_lv(nlv) !< surface area to volume ratio [m2/m3]

    LOGICAL, INTENT(in) :: mask(nc)
    TYPE(t_fwvar_r1d) :: wn !< net fuel load [kg/m2]
    ! ----------------------------------------------------------------------------------------------------- !
    ros = t_sf_ros(nc)
    wn = t_fwvar_r1d(nc)

    !< Corrections from Oberhagemann et al. (2025)
    wfactors = calc_weighting_factors(sigmas_dd, woi_dd, &
                                      sigmas_lv, woi_lv, &
                                      df_sf_oven_dry_part_dens)

    DO inc = 1, nc
      IF (mask(inc)) THEN
        ! CHARACTERISTICS DEAD/LIVE CATEGORIES
        wn%d = SUM(wfactors%gij%d * woi_dd * (1._wp-df_sf_mineral_fract))
        wn%l = SUM(wfactors%gij%l * woi_lv * (1._wp-df_sf_mineral_fract))

        ! FUEL BED CHARACTERISTICS
        ! calc:
        ! - surface are to volume ratio (sigmas)
        ! - mean fuel bulk density (pb)
        ! - mean packing ratio (beta)
        ! - optimum packing ratio (beta_op)
        ! - relative packing ratio (beta/beta_op)
      END IF ! IF (mask(inc)) THEN
    ENDDO ! DO inc = 1, nc

  END FUNCTION sf_calc_rate_of_spread


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
  !> Calculate f weighting factors
  !>
  FUNCTION calc_weighting_factors(sigmas_dead, fuel_load_dead, &
                                  sigmas_live, fuel_load_live, part_dens) &
                                RESULT (wfactors)
    REAL(wp), INTENT(in) :: sigmas_dead(:)      !< surface are to volume ratio (#fuel classes) [m2/m3]
    REAL(wp), INTENT(in) :: fuel_load_dead(:,:) !< oven-dry fuel load (#fuel classes) [kg/m2]
    REAL(wp), INTENT(in) :: sigmas_live(:)      !< surface are to volume ratio (#fuel classes) [m2/m3]
    REAL(wp), INTENT(in) :: fuel_load_live(:,:) !< oven-dry fuel load (#fuel classes) [kg/m2]
    REAL(wp), INTENT(in) :: part_dens           !< particle density [kg/m3]
    TYPE(t_sf_weighting_factors) :: wfactors

    INTEGER :: nc, ndeadfc, nlivefc

    nc = SIZE(fuel_load_dead, DIM=2)
    ndeadfc = SIZE(fuel_load_dead, DIM=1)
    nlivefc = SIZE(fuel_load_live, DIM=1)

    wfactors = t_sf_weighting_factors(nc, ndeadfc, nlivefc)

    DO ic = 1, nc
      wfactors % Aij % d(:,ic) = sigma_dead * fuel_load_dead(:,ic) / part_dens
      wfactors % Aij % l(:,ic) = sigma_live * fuel_load_live(:,ic) / part_dens
      wfactors % Ai % d(ic) = SUM(wfactors % Aij % d(:,ic))
      wfactors % Ai % l(ic) = SUM(wfactors % Aij % l(:,ic))
      wfactors % At(ic) = wfactors % Ai % l(ic) + wfactors % Ai % d(ic)

      wfactors % fij % d(:,ic) = wfactors % Aij % d(:,ic) / wfactors % Ai % d(ic)
      wfactors % fij % l(:,ic) = wfactors % Aij % d(:,ic) / wfactors % Ai % l(ic)
      wfactors % fi % d(:,ic) = SUM(wfactors % fij % d(:,ic))
      wfactors % fi % l(:,ic) = SUM(wfactors % fij % l(:,ic))

      wfactors % gij % d(:,ic) = calc_g_weighting_factors(wfactors % fij % d(:,ic), sigmas_dead)
      wfactors % gij % l(:,ic) = calc_g_weighting_factors(wfactors % fij % l(:,ic), sigmas_live)
    END DO
  END FUNCTION calc_weighting_factors


  ! ======================================================================================================= !
  !>
  !> Calculate g weighting factors
  !>
  FUNCTION calc_g_weighting_factors(fs, sigmas) RESULT (gs)
    REAL(wp), INTENT(in)  :: fs(:)     !< f weighting factor (#fuel classes) [-]
    REAL(wp), INTENT(in)  :: sigmas(:) !< surface area to volume ratio (#fuel classes) [m2/m3]
    REAL(wp), ALLOCATABLE :: gs(:)     !< g weighting factors (#fuel classes) [-]

    INTEGER, PARAMETER :: NGLIM = 5
    INTEGER :: nfuelclass
    !< g partitions, Andrews PL et al. (2018), page 17
    REAL(wp), PARAMETER :: glim(NGLIM) = [16._wp,48._wp,96._wp,1200._wp, HUGE(4.0_wp)]
    REAL(wp) :: fsum(NGLIM)
    INTEGER :: ifuel, ibin

    fsum = 0.
    gs = 0.
    nfuelclass = SIZE(fs)

    ! find f value for each j bin
    do ifuel = 1, NFUELCLASS
      do ibin = 1, NGLIM
        if (sigma_dead(ifuel) < glim(ibin)) exit
      end do
      if (ibin == 1) cycle ! ignore 1st bin, so leave to 0
      fsum(ibin) = fsum(ibin) + fs(ifuel)
    end do

    ! assign g value to each fuel class
    do ifuel = 1, NFUELCLASS
      do ibin = 1, NGLIM
        if (sigma_dead(ifuel) < glim(ibin)) exit
      end do
      gs(ifuel) = fsum(ibin)
    end do
  END FUNCTION calc_g_weighting_factors

#endif
END MODULE mo_q_dfire_sf_ros_ober_class
