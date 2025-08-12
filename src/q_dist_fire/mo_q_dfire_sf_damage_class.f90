!> QUINCY Calculate tree mortality
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
!> @brief module to determine tree mortality
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_q_dfire_sf_damage_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message, finish
  USE mo_jsb_math_constants,  ONLY: eps8, one_day

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_tau_l, calc_prop_ck, t_sf_prob_mort
  PUBLIC :: calc_prob_mortality

  TYPE t_sf_prob_mort
    REAL(wp), ALLOCATABLE :: pm_fire(:) !< mix of tau and ck [-]
    REAL(wp), ALLOCATABLE :: pm_tau(:)  !< probability of post fire mortality due to cambial damage [-]
    REAL(wp), ALLOCATABLE :: pm_ck(:)   !< probability of post fire mortality due to crown damage [-]
    REAL(wp), ALLOCATABLE :: tau_l(:)!< Reside:e time of the fire [min]
    REAL(wp), ALLOCATABLE :: tau_c(:)!< Critical time for cambial damage [min]
    REAL(wp), ALLOCATABLE :: bark_thickness(:) !< [cm]
    REAL(wp), ALLOCATABLE :: prop_ck(:) !< Fraction of the crown scorch [-]
    REAL(wp), ALLOCATABLE :: r_ck !< Resistance factor against crown damage [-]
    REAL(wp), ALLOCATABLE :: ck_p !< Parameter for woody PFTs used in Pm (CK) equation [-]
    REAL(wp), ALLOCATABLE :: scorch_height(:)   !< Scorch height [m]
    REAL(wp), ALLOCATABLE :: F_scorch_height    !< PFT-parameter in crown scorch equation [-]
    REAL(wp), ALLOCATABLE :: crown_length(:)    !< Crown length of woody PFT [m]
  CONTAINS
    PROCEDURE :: copy_from_prob_mort
    GENERIC, PUBLIC :: assignment(=) => copy_from_prob_mort
  END TYPE t_sf_prob_mort


  INTERFACE t_sf_prob_mort
    MODULE PROCEDURE :: create_t_sf_prob_mort
  END INTERFACE t_sf_prob_mort

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_effects'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Constructor for t_sf_prob_mort
  !>
  FUNCTION create_t_sf_prob_mort(nc) RESULT (sf_pm)
    INTEGER, INTENT(in) :: nc
    TYPE(t_sf_prob_mort) :: sf_pm

    ALLOCATE(sf_pm%pm_fire(nc))
    ALLOCATE(sf_pm%pm_tau(nc))
    ALLOCATE(sf_pm%pm_ck(nc))
    ALLOCATE(sf_pm%tau_l(nc))
    ALLOCATE(sf_pm%tau_c(nc))
    ALLOCATE(sf_pm%bark_thickness(nc))
    ALLOCATE(sf_pm%prop_ck(nc))
    ALLOCATE(sf_pm%r_ck)
    ALLOCATE(sf_pm%ck_p)
    ALLOCATE(sf_pm%scorch_height(nc))
    ALLOCATE(sf_pm%F_scorch_height)
    ALLOCATE(sf_pm%crown_length(nc))

    sf_pm%pm_fire = 0._wp
    sf_pm%pm_tau = 0._wp
    sf_pm%pm_ck = 0._wp
    sf_pm%tau_l = 0._wp
    sf_pm%tau_c = 0._wp
    sf_pm%bark_thickness = 0._wp
    sf_pm%prop_ck = 0._wp
    sf_pm%r_ck = 0._wp
    sf_pm%ck_p = 0._wp
    sf_pm%scorch_height = 0._wp
    sf_pm%F_scorch_height = 0._wp
    sf_pm%crown_length = 0._wp
  END FUNCTION create_t_sf_prob_mort


  ! ======================================================================================================= !
  !>
  !> Copy from "other" into "self"
  !>
  PURE SUBROUTINE copy_from_prob_mort(self, other)
    CLASS(t_sf_prob_mort), INTENT(inout) :: self
    CLASS(t_sf_prob_mort), INTENT(in)    :: other

    self % pm_fire = other % pm_fire
    self % pm_tau = other % pm_tau
    self % pm_ck = other % pm_ck
    self % tau_l = other % tau_l
    self % tau_c = other % tau_c
    self % bark_thickness = other % bark_thickness
    self % prop_ck = other % prop_ck
    self % r_ck = other % r_ck
    self % ck_p = other % ck_p
    self % scorch_height = other % scorch_height
    self % F_scorch_height = other % F_scorch_height
    self % crown_length = other % crown_length
  END SUBROUTINE copy_from_prob_mort


  ! ======================================================================================================= !
  !>
  !> Calculate probability of mortality (Section 2.2.6 from Thonicke et al. (2010))
  !> (Tree only)
  !>
  FUNCTION calc_prob_mortality(nc, nfuels, lctlib_crown_length, lctlib_r_ck, lctlib_ck_p, &
                              lctlib_F_scorch_height, lctlib_par1, lctlib_par2, &
                              height, dbh, i_surface, dead_fuel, frac_consump, &
                              mask_opt) RESULT (prob_mort)
    INTEGER, INTENT(in)  :: nc
    INTEGER, INTENT(in)  :: nfuels                          !< Dead fuel classes
    REAL(wp), INTENT(in) :: lctlib_crown_length             !< Crown length of woody PFT [m]
    REAL(wp), INTENT(in) :: height(nc)                      !< Tree height [m]
    REAL(wp), INTENT(in) :: dbh(nc)                         !< Diameter at breast height [cm]
    REAL(wp), INTENT(in) :: i_surface(nc)                   !< Fire intensity [kW m-1]
    REAL(wp), INTENT(in) :: dead_fuel(nc,nfuels)            !< Dead fuel involved in the spread of the fire (1,10,100) [gDM m-2]
    REAL(wp), INTENT(in) :: frac_consump(nc,nfuels)         !< fraction of fuel consumption for spread fuel class(1,10,100) [0 to 1]
    REAL(wp), INTENT(in) :: lctlib_r_ck                     !< Resistant factor of the crown (damage) [0 to 1]
    REAL(wp), INTENT(in) :: lctlib_ck_p                     !< Defoliation factor of the crown scorch [3 to 4]
    REAL(wp), INTENT(in) :: lctlib_F_scorch_height          !< Scorch height parameter [-]
    REAL(wp), INTENT(in) :: lctlib_par1                     !< Bark thickness par1 parameter [-]
    REAL(wp), INTENT(in) :: lctlib_par2                     !< Bark thickness par2 parameter [-]
    LOGICAL, INTENT(in), OPTIONAL :: mask_opt(nc)           !< True = Compute gridcell(s) (default = true)
    TYPE(t_sf_prob_mort) :: prob_mort                       !< output

    LOGICAL :: mask(nc)
    INTEGER :: inc

    mask = .TRUE. !< default
    IF (PRESENT(mask_opt)) THEN
      mask = mask_opt
    END IF
    prob_mort = t_sf_prob_mort(nc) !< Initialize outputs to 0

    DO inc = 1, nc
      IF (mask(inc)) THEN
        !< probability of mortality
        prob_mort % scorch_height(inc) = lctlib_F_scorch_height * (i_surface(inc) ** 0.667) ! (Eq. 16 from Thonicke et al. (2010))
        prob_mort % crown_length(inc) = lctlib_crown_length * height(inc)
        prob_mort % prop_ck(inc) = calc_prop_ck(prob_mort % scorch_height(inc), &
                                                   height(inc), &
                                                   prob_mort % crown_length(inc))

        !< Probability of mortality as a result of the crown scorching (eq. 22 from Thonicke et al. (2010))
        prob_mort % pm_ck(inc) = lctlib_r_ck * (prob_mort % prop_ck(inc) ** lctlib_ck_p)

        !< Calculate Pm(tau) from Eq. 19 in Thonicke et al. (2010).
        prob_mort % tau_l(inc) = calc_tau_l(nfuels, dead_fuel(inc,:), frac_consump(inc,:))
        !< Calculate bark thickness (Eq. 21 from Thonicke et al. (2010))
        prob_mort % bark_thickness(inc) = lctlib_par1 * dbh(inc) + lctlib_par2
        !< Calculate critical time for cambial damage (tau_c - eq. 20 from Thonicke et al. (2010))
        prob_mort % tau_c(inc) = 2.9_wp * (prob_mort % bark_thickness(inc) ** 2)
        prob_mort % pm_tau(inc) = calc_prob_mortality_tau(prob_mort % tau_l(inc), prob_mort % tau_c(inc))

        !< Eq. 18 from Thonicke et al. (2010)
        prob_mort%pm_fire(:) = prob_mort % pm_tau(inc) + prob_mort % pm_ck(inc) - &
                                    prob_mort % pm_tau(inc) * prob_mort % pm_ck(inc)
        prob_mort%r_ck = lctlib_r_ck
        prob_mort%ck_p = lctlib_ck_p
        prob_mort%F_scorch_height = lctlib_F_scorch_height
      END IF
    END DO
  END FUNCTION calc_prob_mortality


  ! ======================================================================================================= !
  !>
  !> Calculate the proportion of the crown scorch (Eq. 17 from Thonicke et al. (2010))
  !>
  !>              ________________________________________________________
  !>           _-_            ^                    ^
  !>        /~~   ~~\         |                    |
  !>     /~~         ~~\      |                    |
  !>    {               }     | Crown length       |                ^
  !>     \  _-     -_  /      |                    | Tree height    |
  !>       ~  \\ //  ~        _                    |                | Scorch height
  !>           | |                                 |                | (Flames height)
  !>           | |                                 |                |
  !>          // \\                                |                |
  !>            ----------------------------------------------------------
  !>
  FUNCTION calc_prop_ck(scorch_height, height, crown_length) RESULT (prop_ck)
    REAL(wp), INTENT(in) :: scorch_height !< [m]
    REAL(wp), INTENT(in) :: height !< [m]
    REAL(wp), INTENT(in) :: crown_length !< [m]
    REAL(wp)             :: prop_ck !< output

    LOGICAL :: are_flames_taller_than_tree
    LOGICAL :: are_flames_shorter_than_crown_length

    are_flames_taller_than_tree = (scorch_height > height)
    !< If flames are not reaching high enough, no crown scorch.
    are_flames_shorter_than_crown_length = (height > (scorch_height + crown_length))
    IF (are_flames_taller_than_tree) THEN
      prop_ck = 1._wp ! Full crown scorch
    ELSE IF (are_flames_shorter_than_crown_length) THEN
      prop_ck = 0._wp ! No crown scorch
    ELSE !< Otheriwse, let's find out how much [%]
      prop_ck = (scorch_height - height + crown_length) / crown_length ! (Eq. 17 from Thonicke et al. (2010))
    ENDIF
  END FUNCTION calc_prop_ck


  ! ======================================================================================================= !
  !>
  !> Calculate probability of mortality due to cambial damage
  !> Eq. 19 from Thonicke et al. (2010)
  !>
  !> In the original Thonicke et al (2010), there were discontinuities and peaks between
  !> the three mortality regimes in Eq 19, and Pm(tau) was sometimes above 1.
  !> Here we adjust the thresholds in order to avoid these small peaks when the three terms of the equation meet.
  !>
  FUNCTION calc_prob_mortality_tau(tau_l, tau_c) RESULT (pm_tau)
    REAL(wp), INTENT(in) :: tau_l, tau_c
    REAL(wp) :: pm_tau !< output [0 ot 1]
    REAL(wp) :: residence_time_ratio !< [0 ot 1]

    residence_time_ratio = tau_l / tau_c

    !< Changed to 0.224 for continuity
    IF (residence_time_ratio <= 0.224_wp) THEN !< Thonicke et al. (2010): 0.22
      pm_tau = 0._wp
      !< Changed to 1.996 for continuity
    ELSE IF (residence_time_ratio >= 1.996_wp) THEN !< Thonicke et al. (2010): 2.0
      pm_tau = 1._wp
    ELSE
      pm_tau = 0.563_wp * residence_time_ratio - 0.125_wp
    ENDIF
  END FUNCTION calc_prob_mortality_tau


  ! ======================================================================================================= !
  !>
  !> Calculate  (tau_l) (Peterson and Ryan 1986 Eq.8 -> Eq. 34 from Pfeiffer et al. (2013))
  !>    The Peterson and Ryan 1986 has the units expressed in cm2. On the other hand, Pfeiffer has in m2.
  !>    A units conversion is required
  !>
  FUNCTION calc_tau_l(nfuels, dead_fuel, frac_consump) RESULT (tau_l)
    USE mo_q_dfire_constants, ONLY: df_sf_units_m2_to_cm2
    INTEGER, INTENT(in)  :: nfuels                    !< Number of fuel classes
    REAL(wp), INTENT(in) :: dead_fuel(nfuels)         !< [gDM/m2]
    REAL(wp), INTENT(in) :: frac_consump(nfuels)      !< Fraction of consumption [0 to 1]
    REAL(wp)             :: tau_l                     !< output [min]

    REAL(wp)             :: tmp(nfuels)               !< Intermediate values
    REAL(wp)             :: frac_consump_sqrd(nfuels) !< [0 to 1]

    frac_consump_sqrd(:) = 1._wp -  ((1._wp - frac_consump(:)) ** 0.5_wp)
    tmp(:) = (dead_fuel(:) * 1._wp/df_sf_units_m2_to_cm2) * frac_consump_sqrd(:)
    tau_l =  39.4_wp * sum(tmp)  ! 2=fuel index

  END FUNCTION calc_tau_l

#endif
END MODULE mo_q_dfire_sf_damage_class
