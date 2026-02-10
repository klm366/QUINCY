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
  !
  USE mo_q_dfire_constants,     ONLY: df_sf_fire_duration_thd
  !<  Rate of spread constants
  USE mo_q_dfire_constants,     ONLY: df_sf_oven_dry_part_dens, &
                                      df_sf_mineral_fract, df_sf_mineral_damp_coeff, &
                                      df_sf_fuel_heat_content
  USE mo_q_dfire_sf_adaptor,    ONLY: sf_convert_units

  USE mo_quincy_stats,          ONLY: stats_average_weight
  USE mo_jsb_math_constants,    ONLY: pi

  USE mo_q_dfire_sf_ros_eqs,    ONLY: calc_opt_react_velocity, calc_prop_flux_ratio, calc_wind_factor
  USE mo_q_dfire_sf_ros_type,   ONLY: t_sf_ros

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sf_calc_rate_of_spread_ober, calc_live_me

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_ros_ober_class'

  INTEGER, PARAMETER :: ix_dead = 1
  INTEGER, PARAMETER :: ix_live = 2

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
    TYPE(t_fwvar_r1d) :: var
    ALLOCATE(var%l(nc))
    ALLOCATE(var%d(nc))

    var%l = 0._wp
    var%d = 0._wp

    IF (PRESENT(has_total)) THEN
    IF (has_total) THEN
      ALLOCATE(var%total(nc))
      var%total = 0._wp
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
    var%d = 0._wp
    var%l = 0._wp
    !    var%t = 0._wp
  END FUNCTION create_fwvar_r2d


  ! ======================================================================================================= !
  !>
  !> calc
  !>
  !> Do the science
  !>
  FUNCTION sf_calc_rate_of_spread_ober(nc, ndd, nlv, &
                                  sigmas_dd, sigmas_lv, &
                                  woi_dd, woi_lv, &
                                  mf_dd, mf_lv, &
                                  me_dd, me_lv, &
                                  fbd, &
                                  u_f_wind, &
                                  mask) RESULT(ros)
    INTEGER, INTENT(in) :: nc
    !< number of dead fuel classes (1h, 10hrs, 100 hrs, dry livegrass)
    !< number of live fuel classes (wet livegrass, woody)
    INTEGER, INTENT(in) :: ndd, nlv
    REAL(wp), INTENT(in) :: woi_dd(nc,ndd), woi_lv(nc,nlv) !< dead/live fuel load [kg DM m-1] -> notes units in kg (instead of g)
    REAL(wp), INTENT(in) :: sigmas_dd(ndd), sigmas_lv(nlv) !< surface area to volume ratio [m2/m3]
    REAL(wp), INTENT(in) :: mf_dd(nc,ndd), mf_lv(nc,nlv) !< Fuel moisture content [-]
    REAL(wp), INTENT(in) :: me_dd(nc), me_lv(nc) !< Moisture of extinction [-]
    REAL(wp), INTENT(in) :: fbd(nc) !< Fuel bulk density [kg/m3]
    REAL(wp), INTENT(in) :: u_f_wind(nc) !< Wind speed at midflame height [m min-1]
    LOGICAL, INTENT(in)  :: mask(nc)
    TYPE(t_sf_ros)       :: ros !< Output: rate of spread [m min-1]

    INTEGER :: ic
    REAL(wp) :: u_forward(nc)
    REAL(wp) :: beta(nc), beta_op(nc), beta_ratio(nc) !< [-], [-], [-]
    REAL(wp) :: gamma(nc), gamma_max(nc) !< [min-1], [min-1]
    REAL(wp) :: A_var(nc) !< [-]
    REAL(wp) :: U_var(nc) !< wind speed at midflame height [m min-1]
    REAL(wp) :: xi(nc) !< propagating flux ratio [-]
    REAL(wp) :: phi_w(nc), phi_s(nc) !< wind factor [-], wind limit [m min-1]
    REAL(wp) :: heat_source(nc) !< Equations nominator [?]
    REAL(wp) :: ros_forward(nc) !< Forward rate of spread [m min-1]
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
    mf = t_fwvar_r1d(nc)
    rm = t_fwvar_r1d(nc)
    eta_m = t_fwvar_r1d(nc)
    sigma = t_fwvar_r1d(nc, has_total=.TRUE.)
    ir = t_fwvar_r1d(nc, has_total=.TRUE.)
    qig_ij = t_fwvar_r2d(nc, ndd, nlv)
    heat_sink = t_fwvar_r1d(nc, has_total=.TRUE.)
    u_forward = 0._wp
    beta = 0._wp
    beta_op = 0._wp
    beta_ratio = 0._wp
    gamma = 0._wp
    gamma_max = 0._wp
    A_var = 0._wp
    phi_w = 0._wp
    phi_s = 0._wp
    heat_source = 0._wp
    ros_forward = 0._wp
    ros = t_sf_ros(nc)
    xi = 0._wp

    !< Corrections from Oberhagemann et al. (2025)
    wfactors = calc_weighting_factors(sigmas_dd, woi_dd, &
                                      sigmas_lv, woi_lv, &
                                      df_sf_oven_dry_part_dens)
    U_var = u_f_wind 

    DO ic = 1, nc
      IF (mask(ic)) THEN
        ! CHARACTERISTICS DEAD/LIVE CATEGORIES
        ! - net fuel load
        wn%d(ic) = SUM(wfactors%gij%d(ic,:) * woi_dd(ic,:) * (1._wp-df_sf_mineral_fract))
        wn%l(ic) = SUM(wfactors%gij%l(ic,:) * woi_lv(ic,:) * (1._wp-df_sf_mineral_fract))
        ! - heat content, constant
        ! - effective mineral content, constant (thonicke=0.01, from Scott and Reinhardt, 2001)
        ! - mineral dampening coefficient, constant (thonicke=0.41739)
        ! - moisture content
        mf%d(ic) = SUM(wfactors%fij%d(ic,:) * mf_dd(ic,:))
        mf%l(ic) = SUM(wfactors%fij%l(ic,:) * mf_lv(ic,:))
        ! - moisture dampening coefficient (oberhagemann et al. (2025))
        rm%d(ic) = MIN(1._wp, mf%d(ic) / me_dd(ic))
        rm%l(ic) = MIN(1._wp, mf%l(ic) / me_lv(ic))
        eta_m%d(ic) = 1._wp - (2.59_wp * rm%d(ic)) + (5.11_wp * (rm%d(ic)**2)) - (3.52_wp * (rm%d(ic)**3))
        eta_m%l(ic) = 1._wp - (2.59_wp * rm%l(ic)) + (5.11_wp * (rm%l(ic)**2)) - (3.52_wp * (rm%l(ic)**3))
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
        beta_op(ic) = 0.200395_wp * (sigma%total(ic) ** (-0.8189_wp))
        ! - relative packing ratio (beta/beta_op)
        beta_ratio(ic) = beta(ic) / beta_op(ic)


        ! HEAT SOURCE (RELATED)
        ! - max reaction velocity
        gamma_max(ic) = 1._wp / (0.0591_wp  + 2.926_wp * sigma%total(ic) ** (-1.5_wp)) !< Table A1, 2nd equation from Thonicke et al. (2010) [min-1]
        ! - optimum reaction velocity
        A_var(ic) = 8.9033_wp * (sigma%total(ic) ** (-0.7913)) !< Table A1, 4th equation from Thonicke et al. (2010)
        gamma(ic) = calc_opt_react_velocity(gamma_max(ic), beta(ic), beta_op(ic), A_var(ic)) ! Table A1, 1st equation from Thonicke et al. (2010) [min-1]
        ! - reaction intensity
        ir%d(ic) =  wn%d(ic) * df_sf_fuel_heat_content * eta_m%d(ic) * df_sf_mineral_damp_coeff !< [kJ m−2 min−1]
        ir%l(ic) =  wn%l(ic) * df_sf_fuel_heat_content * eta_m%l(ic) * df_sf_mineral_damp_coeff !< [kJ m−2 min−1]
        ir%total(ic) = gamma(ic) * (ir%d(ic) + ir%l(ic))
        ! - propagating flux ratio
        xi(ic) = calc_prop_flux_ratio(sigma%total(ic), beta(ic))
        ! - no-wind, no-slope propagating flux

        ! WIND AND SLOPE
        ! - slope factor (defined in Rothermel but not implemented in paper)
        phi_s(ic) = 0._wp
        ! - wind limit
        !< From Rothermel (1972), in Oberhagemann they use ...
        !< wind speed limit according to Patricia L. Andrews et al. 2012.
        IF (U_var(ic) / ir%total(ic) > 0.9_wp) U_var(ic) = 0.9_wp * ir%total(ic)
        ! - wind factor
        phi_w(ic) = calc_wind_factor(U_var(ic), sigma%total(ic), beta(ic), beta_op(ic)) !< Eq. A5 from Thonicke et al. (2010) [-]

        ! HEAT SOURCE
        IF (phi_w(ic) + phi_s(ic) > 0._wp) THEN !< if there is wind speed or slope then account for it
          heat_source(ic) = ir%total(ic) * xi(ic) * (1._wp + phi_w(ic) + phi_s(ic))
        ELSE
          heat_source(ic) = ir%total(ic) * xi(ic)
        END IF

        ! HEAT SINK
        qig_ij%d(ic,:) = 581._wp + 2594._wp * mf_dd(ic,:)
        qig_ij%l(ic,:) = 581._wp + 2594._wp * mf_lv(ic,:)
        heat_sink%d(ic) = SUM(wfactors%fij%d(ic,:) * exp(-4.528_wp / sigmas_dd(:)) * qig_ij%d(ic,:))
        heat_sink%l(ic) = SUM(wfactors%fij%l(ic,:) * exp(-4.528_wp / sigmas_lv(:)) * qig_ij%l(ic,:))
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

  END FUNCTION sf_calc_rate_of_spread_ober


  ! ======================================================================================================= !
  !>
  !> Calculate f weighting factors
  !>
  FUNCTION calc_weighting_factors(sigmas_dead, fuel_load_dead, &
                                  sigmas_live, fuel_load_live, part_dens) &
                                RESULT (wfactors)
    REAL(wp), INTENT(in) :: sigmas_dead(:)      !< surface are to volume ratio (#fuel classes) [m2/m3]
    REAL(wp), INTENT(in) :: fuel_load_dead(:,:) !< oven-dry fuel load (#nc, #fuel classes) [kg/m2]
    REAL(wp), INTENT(in) :: sigmas_live(:)      !< surface are to volume ratio (#fuel classes) [m2/m3]
    REAL(wp), INTENT(in) :: fuel_load_live(:,:) !< oven-dry fuel load (#nc, #fuel classesc) [kg/m2]
    REAL(wp), INTENT(in) :: part_dens           !< particle density [kg/m3]
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

      IF (wfactors % Ai % d(ic) > 0._wp)  &
        wfactors % fij % d(ic,:) = wfactors % Aij % d(ic,:) / wfactors % Ai % d(ic)
      IF (wfactors % Ai % l(ic) > 0._wp)  &
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
    USE mo_q_dfire_constants, ONLY: feet_to_cm => df_sf_units_ft_to_cm 

    REAL(wp), INTENT(in)  :: fs(:)     !< f weighting factor (#fuel classes) [-]
    REAL(wp), INTENT(in)  :: sigmas(:) !< surface area to volume ratio (#fuel classes) [m2/m3]
    REAL(wp), ALLOCATABLE :: gs(:)     !< g weighting factors (#fuel classes) [-]

    INTEGER, PARAMETER :: NGLIM = 5
    INTEGER :: nfuelclass
    !< g partitions, Andrews PL et al. (2018), page 17
    REAL(wp), PARAMETER :: glim(NGLIM) = [16._wp,48._wp,96._wp,1200._wp, HUGE(4.0_wp)] / feet_to_cm 
    REAL(wp) :: fsum(NGLIM)
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
  FUNCTION calc_live_me(nc, ndd, nlv, mf_dd, me_dd, woi_dd, woi_lv, sigmas_dd, sigmas_lv) RESULT (me)
    USE mo_q_dfire_constants, ONLY: df_sf_units_ft_to_cm

    INTEGER, INTENT(in) :: nc, ndd, nlv
    REAL(wp), INTENT(in) :: mf_dd(nc,ndd) !< Dead fuel moisture [-]
    REAL(wp), INTENT(in) :: me_dd(nc) !< Dead fuel moisture of extinction [-]
    REAL(wp), INTENT(in) :: woi_dd(nc,ndd), woi_lv(nc,nlv) !< dead/live fuel load [kg DM m-2]
    REAL(wp), INTENT(in) :: sigmas_dd(ndd), sigmas_lv(nlv) !< Surface area to volume ratio [-]
    REAL(wp) :: me(nc)

    INTEGER :: is
    REAL(wp) :: w(nc) !< dead-to-live load ratio [-]
    REAL(wp) :: wdead(nc,ndd), wlive(nc,nlv) !< dead/live fuel load [kg DM m-2]
    REAL(wp) :: mfdead(nc) !< fine dead fuel moisture [-]
    INTEGER, PARAMETER :: ndd_ff = 2 !< number of fine dead fuel moisture (1hr, livegrass)
    INTEGER, PARAMETER :: ix_ff(ndd_ff) = [1,4] !< 1=1hr, 4=live grass

    w = 0._wp
    mfdead = 0._wp

    ! calc dead-to-live load ratio
    DO is=1, ndd
      wdead(:,is) = woi_dd(:,is) * exp((-138._wp/df_sf_units_ft_to_cm)/sigmas_dd(is))
    END DO
    DO is=1, nlv
      wlive(:,is) = woi_lv(:,is) * exp((-500._wp/df_sf_units_ft_to_cm)/sigmas_lv(is))
    END DO
    WHERE (SUM(wlive, DIM=2) > 0._wp) w = SUM(wdead, DIM=2) / SUM(wlive, DIM=2)

    !< Note: only account for 1hr and dry livegrass
    mfdead = calc_fine_dead_fuel_moisture(nc, ndd_ff, mf_dd(:,ix_ff), woi_dd(:,ix_ff), sigmas_dd(ix_ff))

    me = 2.9_wp * w * (1._wp - mfdead / me_dd) - 0.226_wp
    me = MAX(me, me_dd)
  END FUNCTION calc_live_me


  ! ======================================================================================================= !
  !>
  !> Calculate fine dead fuel moisture
  !> USDA Forest Service RMRS-GTR-371. 2018 Table 6a
  !>
  FUNCTION calc_fine_dead_fuel_moisture(nc, ndd, mf_dd, woi_dd, sigmas_dd) RESULT (mfdead)
    USE mo_q_dfire_constants, ONLY: df_sf_units_ft_to_cm
    INTEGER, INTENT(in) :: nc, ndd
    REAL(wp), INTENT(in) :: mf_dd(nc,ndd) !< Dead fuel moisture [-]
    REAL(wp), INTENT(in) :: woi_dd(nc,ndd) !< Dead fuel load [kg DM m-2]
    REAL(wp), INTENT(in) :: sigmas_dd(ndd) !< Surface area to volume ratio [-]
    REAL(wp) :: mfdead(nc) !< fine dead fuel moisture [-]

    REAL(wp) :: mffm_dd(nc,ndd)
    REAL(wp) :: woi_sum(nc,ndd)
    INTEGER :: is

    DO is=1, ndd
      mffm_dd(:,is) = mf_dd(:,is) * woi_dd(:,is) * exp(-138._wp/df_sf_units_ft_to_cm/sigmas_dd(is))
      woi_sum(:,is) = woi_dd(:,is) * exp(-138._wp/df_sf_units_ft_to_cm/sigmas_dd(is))
    END DO
    mfdead = SUM(mffm_dd) / SUM(woi_sum)
  END FUNCTION calc_fine_dead_fuel_moisture
#endif
END MODULE mo_q_dfire_sf_ros_ober_class
