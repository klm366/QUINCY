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
!>#### routine to calculate fire rate of spread including Oberhagemann improvements from Thonicke's paper
!>
MODULE mo_q_dfire_sf_ros_ober_class
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: lp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  !< Total number of sigmas (1,10,100) hours,  Pfeiffer el al. (2013), Eq. A18 - Used to calculate fuel bulk density
  !
  USE mo_q_dfire_sf_constants,  ONLY: df_sf_fire_duration_thd
  !<  Rate of spread constants
  USE mo_q_dfire_sf_constants,  ONLY: df_sf_oven_dry_part_dens, &
                                      df_sf_mineral_fract, df_sf_mineral_damp_coeff, &
                                      df_sf_fuel_heat_content
  USE mo_q_dfire_sf_adaptor_class,    ONLY: sf_convert_units

  USE mo_quincy_stats,          ONLY: stats_average_weight
  USE mo_jsb_math_constants,    ONLY: pi

  USE mo_q_dfire_sf_util,       ONLY: calc_ros_opt_react_velocity, calc_ros_prop_flux_ratio, calc_ros_wind_factor
  USE mo_q_dfire_sf_classes,    ONLY: t_sf_ros

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sf_calc_rate_of_spread_ober, calc_live_me, calc_weighting_factors, t_fwvar_r1d, t_sf_weighting_factors

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_ros_ober_class'

  INTEGER, PARAMETER :: ix_dead = 1
  INTEGER, PARAMETER :: ix_live = 2

  !< weight factor variable
  TYPE :: t_fwvar_r1d
    REAL(lp), ALLOCATABLE :: d(:) !< dead (#nc)
    REAL(lp), ALLOCATABLE :: l(:) !< live (#nc)
    REAL(lp), ALLOCATABLE :: total(:) !< total, or dead + live (#nc)
  END TYPE t_fwvar_r1d

  INTERFACE t_fwvar_r1d
    MODULE PROCEDURE :: create_fwvar_r1d
  END INTERFACE t_fwvar_r1d

  !< weight factor variable
  TYPE :: t_fwvar_r2d
    REAL(lp), ALLOCATABLE :: d(:,:) !< dead (#nc)
    REAL(lp), ALLOCATABLE :: l(:,:) !< live (#nc)
    REAL(lp), ALLOCATABLE :: total(:,:) !< total, or dead + live (#nc)
  END TYPE t_fwvar_r2d

  INTERFACE t_fwvar_r2d
    MODULE PROCEDURE :: create_fwvar_r2d
  END INTERFACE t_fwvar_r2d


  TYPE :: t_sf_weighting_factors
    TYPE(t_fwvar_r2d)     :: Aij    !< (#fuel classes-dead, live-, #nc)
    TYPE(t_fwvar_r1d)     :: Ai     !< (#nc)
    REAL(lp), ALLOCATABLE :: At(:)  !< A_dead + A_live (#nc)
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
    wfactors % At = 0._lp

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
    TYPE(t_fwvar_r1d) :: var
    ALLOCATE(var%l(nc))
    ALLOCATE(var%d(nc))

    var%l = 0._lp
    var%d = 0._lp

    IF (PRESENT(has_total)) THEN
    IF (has_total) THEN
      ALLOCATE(var%total(nc))
      var%total = 0._lp
    ENDIF
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
    TYPE(t_fwvar_r2d) :: var
    ALLOCATE(var%d(nc, nddfc))
    ALLOCATE(var%l(nc, nlvfc))
    !    ALLOCATE(var%t(nc))
    var%d = 0._lp
    var%l = 0._lp
    !    var%t = 0._lp
  END FUNCTION create_fwvar_r2d


  ! ======================================================================================================= !
  !>
  !> calc
  !>
  !> Do the science
  !>
  SUBROUTINE sf_calc_rate_of_spread_ober(nc, ndd, nlv, &
                                  sigmas_dd, sigmas_lv, &
                                  woi_dd, woi_lv, &
                                  mf_dd, mf_lv, &
                                  me_dd, me_lv, &
                                  fbd, &
                                  u_f_wind, &
                                  mask, &
                                  ros)
    INTEGER, INTENT(in) :: nc
    !< number of dead fuel classes (1h, 10hrs, 100 hrs, dry livegrass)
    !< number of live fuel classes (wet livegrass, woody)
    INTEGER, INTENT(in) :: ndd, nlv
    REAL(lp), INTENT(in) :: woi_dd(nc,ndd), woi_lv(nc,nlv) !< dead/live fuel load [kg DM m-1] -> notes units in kg (instead of g)
    REAL(lp), INTENT(in) :: sigmas_dd(ndd), sigmas_lv(nlv) !< surface area to volume ratio [m2/m3]
    REAL(lp), INTENT(in) :: mf_dd(nc,ndd), mf_lv(nc,nlv) !< Fuel moisture content [-]
    REAL(lp), INTENT(in) :: me_dd(nc), me_lv(nc) !< Moisture of extinction [-]
    REAL(lp), INTENT(in) :: fbd(nc) !< Fuel bulk density [kg/m3]
    REAL(lp), INTENT(in) :: u_f_wind(nc) !< Wind speed at midflame height [m min-1]
    LOGICAL, INTENT(in)  :: mask(nc)
    TYPE(t_sf_ros), INTENT(inout) :: ros !< Output: rate of spread [m min-1]

    INTEGER :: ic
    REAL(lp) :: u_forward(nc)
    REAL(lp) :: beta(nc), beta_op(nc), beta_ratio(nc) !< [-], [-], [-]
    REAL(lp) :: gamma(nc), gamma_max(nc) !< [min-1], [min-1]
    REAL(lp) :: A_var(nc) !< [-]
    REAL(lp) :: U_var(nc) !< wind speed at midflame height [m min-1]
    REAL(lp) :: xi(nc) !< propagating flux ratio [-]
    REAL(lp) :: phi_w(nc), phi_s(nc) !< wind factor [-], wind limit [m min-1]
    REAL(lp) :: heat_source(nc) !< Equations nominator [?]
    REAL(lp) :: ros_forward(nc) !< Forward rate of spread [m min-1]
    TYPE(t_fwvar_r1d) :: wn !< net fuel load [kg/m2]
    TYPE(t_fwvar_r1d) :: mf !< moisture content [-]
    TYPE(t_fwvar_r1d) :: rm !< fuel moisture ratio [0 to 1]
    TYPE(t_fwvar_r1d) :: eta_m !< moisture dampening coefficient [?]
    TYPE(t_fwvar_r1d) :: sigma !< surface area to volume ratio [m2/m3]
    TYPE(t_fwvar_r1d) :: ir !< reaction intensity [KJ m-2 min-1]
    TYPE(t_fwvar_r2d) :: qig_ij !< Heat of preignition [kJ kg-1]
    TYPE(t_fwvar_r1d) :: heat_sink !< Equations denominator [?]
    TYPE(t_sf_weighting_factors) :: wfactors
    ! ----------------------------------------------------------------------------------------------------- !
    wn = t_fwvar_r1d(nc)
    mf = t_fwvar_r1d(nc, has_total=.TRUE.)
    rm = t_fwvar_r1d(nc)
    eta_m = t_fwvar_r1d(nc)
    sigma = t_fwvar_r1d(nc, has_total=.TRUE.)
    ir = t_fwvar_r1d(nc, has_total=.TRUE.)
    qig_ij = t_fwvar_r2d(nc, ndd, nlv)
    heat_sink = t_fwvar_r1d(nc, has_total=.TRUE.)
    u_forward = 0._lp
    beta = 0._lp
    beta_op = 0._lp
    beta_ratio = 0._lp
    gamma = 0._lp
    gamma_max = 0._lp
    A_var = 0._lp
    phi_w = 0._lp
    phi_s = 0._lp
    heat_source = 0._lp
    ros_forward = 0._lp
    xi = 0._lp

    !< Corrections from Oberhagemann et al. (2025)
    wfactors = calc_weighting_factors(sigmas_dd, woi_dd, &
                                      sigmas_lv, woi_lv, &
                                      df_sf_oven_dry_part_dens)
    U_var = u_f_wind

    DO ic = 1, nc
      IF (mask(ic)) THEN
        ! CHARACTERISTICS DEAD/LIVE CATEGORIES
        ! - net fuel load
        wn%d(ic) = SUM(wfactors%gij%d(ic,:) * woi_dd(ic,:) * (1._lp-df_sf_mineral_fract))
        wn%l(ic) = SUM(wfactors%gij%l(ic,:) * woi_lv(ic,:) * (1._lp-df_sf_mineral_fract))
        ! - heat content, constant
        ! - effective mineral content, constant (thonicke=0.01, from Scott and Reinhardt, 2001)
        ! - mineral dampening coefficient, constant (thonicke=0.41739)
        ! - moisture content
        mf%d(ic) = SUM(wfactors%fij%d(ic,:) * mf_dd(ic,:))
        mf%l(ic) = SUM(wfactors%fij%l(ic,:) * mf_lv(ic,:))
        mf%total(ic) = wfactors%fi%l(ic) * mf%l(ic) + wfactors%fi%d(ic) * mf%d(ic)
        ! - moisture dampening coefficient (oberhagemann et al. (2025))
        IF (me_dd(ic) > 0._lp) rm%d(ic) = MIN(1._lp, mf%d(ic) / me_dd(ic))
        IF (me_lv(ic) > 0._lp) rm%l(ic) = MIN(1._lp, mf%l(ic) / me_lv(ic))
        eta_m%d(ic) = 1._lp - (2.59_lp * rm%d(ic)) + (5.11_lp * (rm%d(ic)**2)) - (3.52_lp * (rm%d(ic)**3))
        eta_m%l(ic) = 1._lp - (2.59_lp * rm%l(ic)) + (5.11_lp * (rm%l(ic)**2)) - (3.52_lp * (rm%l(ic)**3))
        ! - surface area to volume ratio (sigma_i)
        sigma%d(ic) = SUM(wfactors%fij%d(ic,:) * sigmas_dd)
        sigma%l(ic) = SUM(wfactors%fij%l(ic,:) * sigmas_lv)

        ! FUEL BED CHARACTERISTICS
        ! calc:
        ! - surface are to volume ratio (sigmas)
        sigma%total(ic) = wfactors%fi%d(ic) * sigma%d(ic) + wfactors%fi%l(ic) * sigma%l(ic)
        ! - mean fuel bulk density (pb), nothing to calculate
        ! - mean packing ratio (beta)
        beta(ic) = fbd(ic) / df_sf_oven_dry_part_dens
        ! - optimum packing ratio (beta_op)
        beta_op(ic) = 0.200395_lp * (sigma%total(ic) ** (-0.8189_lp))
        ! - relative packing ratio (beta/beta_op)
        beta_ratio(ic) = beta(ic) / beta_op(ic)


        ! HEAT SOURCE (RELATED)
        ! - max reaction velocity
        gamma_max(ic) = 1._lp / (0.0591_lp  + 2.926_lp * sigma%total(ic) ** (-1.5_lp)) !< Table A1, 2nd equation from Thonicke et al. (2010) [min-1]
        ! - optimum reaction velocity
        A_var(ic) = 8.9033_lp * (sigma%total(ic) ** (-0.7913)) !< Table A1, 4th equation from Thonicke et al. (2010)
        gamma(ic) = calc_ros_opt_react_velocity(gamma_max(ic), beta(ic), beta_op(ic), A_var(ic)) ! Table A1, 1st equation from Thonicke et al. (2010) [min-1]
        ! - reaction intensity
        ir%d(ic) =  wn%d(ic) * df_sf_fuel_heat_content * eta_m%d(ic) * df_sf_mineral_damp_coeff !< [kJ m-2 min-1]
        ir%l(ic) =  wn%l(ic) * df_sf_fuel_heat_content * eta_m%l(ic) * df_sf_mineral_damp_coeff !< [kJ m-2 min-1]
        ir%total(ic) = gamma(ic) * (ir%d(ic) + ir%l(ic))
        ! - propagating flux ratio
        xi(ic) = calc_ros_prop_flux_ratio(sigma%total(ic), beta(ic))
        ! - no-wind, no-slope propagating flux

        ! WIND AND SLOPE
        ! - slope factor (defined in Rothermel but not implemented in paper)
        phi_s(ic) = 0._lp
        ! - wind limit
        !< From Rothermel (1972), in Oberhagemann they use ...
        !< wind speed limit according to Patricia L. Andrews et al. 2012.
        IF (U_var(ic) / ir%total(ic) > 0.9_lp) U_var(ic) = 0.9_lp * ir%total(ic)
        ! - wind factor
        phi_w(ic) = calc_ros_wind_factor(U_var(ic), sigma%total(ic), beta(ic), beta_op(ic)) !< Eq. A5 from Thonicke et al. (2010) [-]

        ! HEAT SOURCE
        IF (phi_w(ic) + phi_s(ic) > 0._lp) THEN !< if there is wind speed or slope then account for it
          heat_source(ic) = ir%total(ic) * xi(ic) * (1._lp + phi_w(ic) + phi_s(ic))
        ELSE
          heat_source(ic) = ir%total(ic) * xi(ic)
        END IF

        ! HEAT SINK
        qig_ij%d(ic,:) = 581._lp + 2594._lp * mf_dd(ic,:)
        qig_ij%l(ic,:) = 581._lp + 2594._lp * mf_lv(ic,:)
        heat_sink%d(ic) = SUM(wfactors%fij%d(ic,:) * exp(-4.528_lp / sigmas_dd(:)) * qig_ij%d(ic,:))
        heat_sink%l(ic) = SUM(wfactors%fij%l(ic,:) * exp(-4.528_lp / sigmas_lv(:)) * qig_ij%l(ic,:))
        heat_sink%total(ic) = fbd(ic) * (wfactors%fi%d(ic) * heat_sink%d(ic) + wfactors%fi%l(ic) * heat_sink%l(ic))

        ros%ros_forward(ic) = heat_source(ic) / heat_sink%total(ic)
      END IF ! IF (mask(ic)) THEN
    ENDDO ! DO ic = 1, nc

!    ros%ros_forward =  ros_forward
    ros%reaction_intensity = ir%total
    ros%prop_flux_ratio = xi
    ros%wind_factor=phi_w
    !ros%eff_heting_number = epsilon ! todo: exp(-4.528 / sigmas)
    !ros%pre_ign_heat = qig ! dead, live
    ros%opt_react_velo = gamma
    !ros%max_react_velo = gamma_max
    !ros%net_fuel_load = wn ! dead, live
    !ros%moist_damp_coeff  = eta_m  ! dead, live
    ros%mineral_damp_coeff = df_sf_mineral_damp_coeff

  END SUBROUTINE sf_calc_rate_of_spread_ober


  ! ======================================================================================================= !
  !>
  !> Calculate f weighting factors
  !>
  FUNCTION calc_weighting_factors(sigmas_dead, fuel_load_dead, &
                                  sigmas_live, fuel_load_live, part_dens) &
                                RESULT (wfactors)
    REAL(lp), INTENT(in) :: sigmas_dead(:)      !< surface are to volume ratio (#fuel classes) [m2/m3]
    REAL(lp), INTENT(in) :: fuel_load_dead(:,:) !< oven-dry fuel load (#nc, #fuel classes) [kg/m2]
    REAL(lp), INTENT(in) :: sigmas_live(:)      !< surface are to volume ratio (#fuel classes) [m2/m3]
    REAL(lp), INTENT(in) :: fuel_load_live(:,:) !< oven-dry fuel load (#nc, #fuel classesc) [kg/m2]
    REAL(lp), INTENT(in) :: part_dens           !< particle density [kg/m3]
    TYPE(t_sf_weighting_factors) :: wfactors

    INTEGER :: nc, ndeadfc, nlivefc, ic

    nc = SIZE(fuel_load_dead, DIM=1)
    ndeadfc = SIZE(fuel_load_dead, DIM=2)
    nlivefc = SIZE(fuel_load_live, DIM=2)

    wfactors = t_sf_weighting_factors(nc, ndeadfc, nlivefc)

    DO ic = 1, nc
      wfactors % Aij % d(ic,:) = sigmas_dead * fuel_load_dead(ic,:) / part_dens
      wfactors % Aij % l(ic,:) = sigmas_live * fuel_load_live(ic,:) / part_dens
      wfactors % Ai % d(ic) = SUM(wfactors % Aij % d(ic,:))
      wfactors % Ai % l(ic) = SUM(wfactors % Aij % l(ic,:))
      wfactors % At(ic) = wfactors % Ai % l(ic) + wfactors % Ai % d(ic)

      IF (wfactors % Ai % d(ic) > 0._lp)  &
        wfactors % fij % d(ic,:) = wfactors % Aij % d(ic,:) / wfactors % Ai % d(ic)
      IF (wfactors % Ai % l(ic) > 0._lp)  &
        wfactors % fij % l(ic,:) = wfactors % Aij % l(ic,:) / wfactors % Ai % l(ic)
      wfactors % fi % d(ic) = SUM(wfactors % fij % d(ic,:))
      wfactors % fi % l(ic) = SUM(wfactors % fij % l(ic,:))

      wfactors % gij % d(ic,:) = calc_g_weighting_factors(wfactors % fij % d(ic,:), sigmas_dead)
      wfactors % gij % l(ic,:) = calc_g_weighting_factors(wfactors % fij % l(ic,:), sigmas_live)
    END DO

  END FUNCTION calc_weighting_factors


  ! ======================================================================================================= !
  !>
  !> Calculate g weighting factors
  !>
  FUNCTION calc_g_weighting_factors(fs, sigmas) RESULT (gs)
    USE mo_q_dfire_sf_constants, ONLY: feet_to_cm => df_sf_units_ft_to_cm

    REAL(lp), INTENT(in)  :: fs(:)     !< f weighting factor (#fuel classes) [-]
    REAL(lp), INTENT(in)  :: sigmas(:) !< surface area to volume ratio (#fuel classes) [m2/m3]
    REAL(lp), ALLOCATABLE :: gs(:)     !< g weighting factors (#fuel classes) [-]

    INTEGER, PARAMETER :: NGLIM = 5
    INTEGER :: nfuelclass
    !< g partitions, Andrews PL et al. (2018), page 17
    REAL(lp), PARAMETER :: glim(NGLIM) = [16._lp,48._lp,96._lp,1200._lp, HUGE(4.0_lp)] / feet_to_cm
    REAL(lp) :: fsum(NGLIM)
    INTEGER :: ifuel, ibin

    ALLOCATE(gs(SIZE(fs)))
    fsum = 0.
    gs = 0.
    nfuelclass = SIZE(fs)

    ! find f value for each j bin
    do ifuel = 1, NFUELCLASS
      do ibin = 1, NGLIM
        if (sigmas(ifuel) < glim(ibin)) exit
      end do
      if (ibin == 1) cycle ! ignore 1st bin, so leave to 0
      fsum(ibin) = fsum(ibin) + fs(ifuel)
    end do

    ! assign g value to each fuel class
    do ifuel = 1, NFUELCLASS
      do ibin = 1, NGLIM
        if (sigmas(ifuel) < glim(ibin)) exit
      end do
      gs(ifuel) = fsum(ibin)
    end do
  END FUNCTION calc_g_weighting_factors


  ! ======================================================================================================= !
  !>
  !> Calculate live fuel moisture of extinction
  !> USDA Forest Service RMRS-GTR-371. 2018 Table 6a
  !>
  FUNCTION calc_live_me(nc, ndd, nlv, mf_dd, me_dd, woi_dd, woi_lv, sigmas_dd, sigmas_lv, mask) RESULT (me)
    USE mo_q_dfire_sf_constants, ONLY: df_sf_units_ft_to_cm

    INTEGER, INTENT(in) :: nc, ndd, nlv
    REAL(lp), INTENT(in) :: mf_dd(nc,ndd) !< Dead fuel moisture [-]
    REAL(lp), INTENT(in) :: me_dd(nc) !< Dead fuel moisture of extinction [-]
    REAL(lp), INTENT(in) :: woi_dd(nc,ndd), woi_lv(nc,nlv) !< dead/live fuel load [kg DM m-2]
    REAL(lp), INTENT(in) :: sigmas_dd(ndd), sigmas_lv(nlv) !< Surface area to volume ratio [-]
    LOGICAL, INTENT(in)  :: mask(nc)
    REAL(lp) :: me(nc)

    INTEGER :: is, ic
    REAL(lp) :: w(nc) !< dead-to-live load ratio [-]
    REAL(lp) :: wdead(nc,ndd), wlive(nc,nlv) !< dead/live fuel load [kg DM m-2]
    REAL(lp) :: wdead_sum(nc), wlive_sum(nc) !< dead/live fuel load [kg DM m-2]
    REAL(lp) :: mfdead(nc) !< fine dead fuel moisture [-]
    INTEGER, PARAMETER :: ndd_ff = 2 !< number of fine dead fuel moisture (1hr, livegrass)
    INTEGER, PARAMETER :: ix_ff(ndd_ff) = [1,4] !< 1=1hr, 4=dead live grass

    w = 0._lp
    mfdead = 0._lp
    me = 0._lp
    wdead = 0._lp
    wlive = 0._lp

    ! calc dead-to-live load ratio
    DO is=1, ndd
      DO ic=1, nc
        IF (mask(ic)) THEN
          wdead(ic,is) = woi_dd(ic,is) * exp((-138._lp/df_sf_units_ft_to_cm)/sigmas_dd(is))
        ENDIF
      END DO
    END DO

    DO is=1, nlv
      DO ic=1, nc
        IF (mask(ic)) THEN
          wlive(ic,is) = woi_lv(ic,is) * exp((-500._lp/df_sf_units_ft_to_cm)/sigmas_lv(is))
        ENDIF
      ENDDO
    END DO

    wdead_sum = SUM(wdead, DIM=2)
    wlive_sum = SUM(wlive, DIM=2)

    DO ic=1, nc
      IF (mask(ic) .AND. wlive_sum(ic) > 0._lp) THEN
        w(ic) = wdead_sum(ic) / wlive_sum(ic)
      ENDIF
    ENDDO

    !< Note: only account for 1hr and dry livegrass
    mfdead = calc_fine_dead_fuel_moisture(nc, ndd_ff, mf_dd(:,ix_ff), woi_dd(:,ix_ff), sigmas_dd(ix_ff), mask)

    DO ic=1, nc
      IF (mask(ic) .AND. wlive_sum(ic) > 0._lp) THEN
        me(ic) = 2.9_lp * w(ic) * (1._lp - mfdead(ic) / me_dd(ic)) - 0.226_lp
        me(ic) = MAX(me(ic), me_dd(ic))
      ENDIF
    ENDDO
  END FUNCTION calc_live_me


  ! ======================================================================================================= !
  !>
  !> Calculate fine dead fuel moisture
  !> USDA Forest Service RMRS-GTR-371. 2018 Table 6a
  !>
  FUNCTION calc_fine_dead_fuel_moisture(nc, ndd, mf_dd, woi_dd, sigmas_dd, mask) RESULT (mfdead)
    USE mo_q_dfire_sf_constants, ONLY: df_sf_units_ft_to_cm
    INTEGER, INTENT(in) :: nc, ndd
    REAL(lp), INTENT(in) :: mf_dd(nc,ndd) !< Dead fuel moisture [-]
    REAL(lp), INTENT(in) :: woi_dd(nc,ndd) !< Dead fuel load [kg DM m-2]
    REAL(lp), INTENT(in) :: sigmas_dd(ndd) !< Surface area to volume ratio [-]
    LOGICAL, INTENT(in)  :: mask(nc)
    REAL(lp) :: mfdead(nc) !< fine dead fuel moisture [-]

    REAL(lp) :: mffm_dd(nc,ndd)
    REAL(lp) :: woi(nc,ndd)
    REAL(lp) :: woi_sum(nc)
    REAL(lp) :: mffm_sum(nc)
    INTEGER :: is, ic

    mffm_dd = 0._lp
    woi = 0._lp
    mfdead = 0._lp

    DO is=1, ndd
      DO ic=1, nc
        IF (mask(ic)) THEN
          mffm_dd(ic,is) = mf_dd(ic,is) * woi_dd(ic,is) * exp(-138._lp/df_sf_units_ft_to_cm/sigmas_dd(is))
          woi(ic,is) = woi_dd(ic,is) * exp(-138._lp/df_sf_units_ft_to_cm/sigmas_dd(is))
        ENDIF
      ENDDO
    END DO

    mffm_sum = SUM(mffm_dd, DIM=2)
    woi_sum = SUM(woi, DIM=2)

    DO ic=1, nc
      IF (mask(ic) .AND. woi_sum(ic) > 0._lp) THEN
        mfdead(ic) = mffm_sum(ic) / woi_sum(ic)
      ENDIF
    ENDDO
  END FUNCTION calc_fine_dead_fuel_moisture
#endif
END MODULE mo_q_dfire_sf_ros_ober_class
