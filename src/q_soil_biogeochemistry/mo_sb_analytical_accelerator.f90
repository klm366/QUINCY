!> Analytical steady-state accelerator for QUINCY soil biogeochemistry
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
MODULE mo_sb_analytical_accelerator
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message
  USE mo_jsb_math_constants, ONLY: eps4, eps8, one_day, one_year
  USE mo_lnd_bgcm_idx

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_analytical_ss_jsm, apply_analytical_ss_jsm_accelerator

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_analytical_accelerator' !< module name used in log contexts

CONTAINS

  ! ======================================================================================================= !
  !> evaluate one JSM analytical steady-state layer for a prescribed microbial biomass
  !> The helper projects the biogeochemical solution onto an admissible soil-volume state
  !> without changing the forward model equations.
  !>
  SUBROUTINE calc_analytical_ss_jsm_layer_state( &
    & l_surface_layer, &
    & input_sol, &
    & input_poly, &
    & input_wood, &
    & dom_carbon, &
    & microbial_carbon, &
    & l_microbial_dormant, &
    & enzyme_frac_ap_layer, &
    & enzyme_frac_total_scale, &
    & particle_poly_c, &
    & particle_residue_c, &
    & particle_dom_assoc_c, &
    & particle_residue_assoc_c, &
    & km_depoly_eff, &
    & vmax_poly_eff, &
    & vmax_res_eff, &
    & qmax_org_min_local, &
    & bulk_dens_local, &
    & rtm_sorb_eff, &
    & rmm_sorb_eff, &
    & rtm_desorb_eff, &
    & rmm_desorb_eff, &
    & alpha_poly, &
    & alpha_residue, &
    & polymeric_carbon, &
    & residue_carbon, &
    & dom_assoc_carbon, &
    & residue_assoc_carbon, &
    & qmax_org_layer, &
    & bulk_carbon_layer)

    USE mo_q_sb_jsm_processes,      ONLY: calc_enzyme_allocation
    USE mo_sb_constants
    USE mo_jsb_physical_constants,  ONLY: molar_mass_C

    LOGICAL,  INTENT(in)  :: l_surface_layer             !< TRUE for the uppermost mineral soil layer
    REAL(wp), INTENT(in)  :: input_sol                   !< soluble-litter C input rate [mol C m-3 s-1]
    REAL(wp), INTENT(in)  :: input_poly                  !< polymeric-litter C input rate [mol C m-3 s-1]
    REAL(wp), INTENT(in)  :: input_wood                  !< woody-litter C input rate [mol C m-3 s-1]
    REAL(wp), INTENT(in)  :: dom_carbon                  !< dissolved organic matter C density [mol C m-3]
    REAL(wp), INTENT(in)  :: microbial_carbon            !< prescribed microbial biomass C density [mol C m-3]
    LOGICAL,  INTENT(in)  :: l_microbial_dormant         !< TRUE when biomass is held at its dormant minimum
    REAL(wp), INTENT(in)  :: enzyme_frac_ap_layer        !< phosphatase fraction of total enzyme investment [-]
    REAL(wp), INTENT(in)  :: enzyme_frac_total_scale     !< physical-constraint multiplier of total enzyme fraction [-]
    REAL(wp), INTENT(in)  :: particle_poly_c             !< net polymeric-litter transport tendency [mol C m-3 s-1]
    REAL(wp), INTENT(in)  :: particle_residue_c          !< net residue transport tendency [mol C m-3 s-1]
    REAL(wp), INTENT(in)  :: particle_dom_assoc_c        !< net mineral-associated DOM tendency [mol C m-3 s-1]
    REAL(wp), INTENT(in)  :: particle_residue_assoc_c    !< net mineral-associated residue tendency [mol C m-3 s-1]
    REAL(wp), INTENT(in)  :: km_depoly_eff               !< modified depolymerisation half-saturation [mol C m-3]
    REAL(wp), INTENT(in)  :: vmax_poly_eff               !< modified polymeric-litter depolymerisation rate [s-1]
    REAL(wp), INTENT(in)  :: vmax_res_eff                !< modified residue depolymerisation rate [s-1]
    REAL(wp), INTENT(in)  :: qmax_org_min_local          !< texture-dependent sorption capacity [mol C kg-1]
    REAL(wp), INTENT(in)  :: bulk_dens_local             !< mineral soil bulk density [kg m-3]
    REAL(wp), INTENT(in)  :: rtm_sorb_eff                !< temperature modifier of sorption [-]
    REAL(wp), INTENT(in)  :: rmm_sorb_eff                !< moisture modifier of sorption [-]
    REAL(wp), INTENT(in)  :: rtm_desorb_eff              !< temperature modifier of desorption [-]
    REAL(wp), INTENT(in)  :: rmm_desorb_eff              !< moisture modifier of desorption [-]
    REAL(wp), INTENT(out) :: alpha_poly                  !< enzyme allocation fraction to polymeric litter [-]
    REAL(wp), INTENT(out) :: alpha_residue               !< enzyme allocation fraction to residue [-]
    REAL(wp), INTENT(out) :: polymeric_carbon            !< analytical polymeric-litter C density [mol C m-3]
    REAL(wp), INTENT(out) :: residue_carbon              !< analytical microbial-residue C density [mol C m-3]
    REAL(wp), INTENT(out) :: dom_assoc_carbon            !< analytical mineral-associated DOM C [mol C m-3]
    REAL(wp), INTENT(out) :: residue_assoc_carbon        !< analytical mineral-associated residue C [mol C m-3]
    REAL(wp), INTENT(out) :: qmax_org_layer              !< available organic-matter sorption capacity [mol C m-3]
    REAL(wp), INTENT(out) :: bulk_carbon_layer           !< C occupying the mineral-soil layer volume [mol C m-3]

    REAL(wp) :: soluable_carbon              !< analytical soluble-litter C density [mol C m-3]
    REAL(wp) :: woody_carbon                 !< analytical woody-litter C density [mol C m-3]
    REAL(wp) :: enzyme_depoly                !< depolymerisation-enzyme C density [mol C m-3]
    REAL(wp) :: enzyme_residue               !< effective residue-depolymerisation enzyme C density [mol C m-3]
    REAL(wp) :: enzyme_frac_total_local      !< scaled microbial enzyme investment fraction [-]
    REAL(wp) :: d_litter                     !< polymeric-substrate source including transport [mol C m-3 s-1]
    REAL(wp) :: d_residue                    !< residue source including transport [mol C m-3 s-1]
    REAL(wp) :: qmax_a                       !< texture capacity expressed per soil volume [mol C m-3]
    REAL(wp) :: qmax_b                       !< C density of pure organic matter [mol C m-3]
    REAL(wp) :: c0_star                      !< non-associated SOM C density occupying soil volume [mol C m-3]
    REAL(wp) :: ratio_dom                    !< DOM adsorption-to-desorption coefficient ratio [m3 mol-1]
    REAL(wp) :: ratio_res                    !< residue adsorption-to-desorption coefficient ratio [m3 mol-1]
    REAL(wp) :: qmax_f                       !< occupied fraction in the sorption-capacity closure [-]
    REAL(wp) :: sorption_denom               !< competitive sorption denominator [-]
    REAL(wp) :: beta_dom                     !< transported associated-DOM offset [mol C m-3]
    REAL(wp) :: beta_residue                 !< transported associated-residue offset [mol C m-3]
    REAL(wp) :: beta_total                   !< sum of associated-pool transport offsets [mol C m-3]
    REAL(wp) :: free_sorption_capacity       !< unoccupied equilibrium sorption capacity [mol C m-3]
    LOGICAL  :: l_particle_tendency          !< TRUE when any particle-transport tendency is nonzero

    !> C_wl*  = I_wl [mol C m-3 s-1] * tau_wl [s]
    woody_carbon    = input_wood * tau_woody_litter
    !> C_sol* = I_sol [mol C m-3 s-1] * tau_sol [s]
    soluable_carbon = input_sol  * tau_soluable_litter

    enzyme_frac_total_local = MAX(0.0_wp, enzyme_frac_total * enzyme_frac_total_scale)
    !> Retain a background extracellular-enzyme capacity when living microbial biomass is very small.
    !> E_depoly* = max(C_B* * enzyme_frac_total * (1 - enzyme_frac_AP), E_persistent,min) [mol C m-3]
    enzyme_depoly = MAX(microbial_carbon * enzyme_frac_total_local * (1._wp - enzyme_frac_ap_layer), &
      &                 min_persistent_enzyme)
    !> d_L = I_poly + eta_poly^wl * I_wl  [mol C m-3 s-1]
    d_litter        = MAX(0.0_wp, input_poly + fc_woody2polymeric * input_wood + particle_poly_c)
    !> d_R = (1 - f_soluable_residue) * k_B * C_B*  [mol C m-3 s-1]
    IF (l_microbial_dormant) THEN
      d_residue = MAX(0.0_wp, particle_residue_c + particle_residue_assoc_c)
    ELSE
      d_residue = MAX(0.0_wp, &
        & (1._wp - f_soluable_residue) * k_mort_microbial * microbial_carbon + &
        & particle_residue_c + particle_residue_assoc_c)
    END IF

    alpha_poly = calc_enzyme_allocation(d_litter, d_residue, km_depoly_eff, enzyme_depoly)

    IF (alpha_poly /= alpha_poly) alpha_poly = init_enzyme_fraction
    alpha_poly    = MIN(MAX(alpha_poly, min_enzyme_fraction), 1._wp - min_enzyme_fraction)
    alpha_residue = 1._wp - alpha_poly
    enzyme_residue = MAX(alpha_residue * enzyme_depoly, min_persistent_enzyme)

    !> C_poly*=(I_poly + eta_poly^wl * I_wl) [mol C m-3 s-1]
    !>         / (vmax_poly [s-1] * alpha_poly [-] * E_depoly* [mol C m-3])
    !>         * (Km_depoly [mol C m-3] + alpha_poly * E_depoly*)
    IF (d_litter > 0.0_wp .AND. alpha_poly * enzyme_depoly > 0.0_wp .AND. vmax_poly_eff > 0.0_wp) THEN
      polymeric_carbon = d_litter * (km_depoly_eff + alpha_poly * enzyme_depoly) / &
        &                (vmax_poly_eff * alpha_poly * enzyme_depoly)
    ELSE
      polymeric_carbon = eps8
    END IF

    !> C_res*=((1 - f_soluable_residue) * k_B * C_B*) [mol C m-3 s-1]
    !>        / (vmax_res [s-1] * E_res* [mol C m-3])
    !>        * (Km_depoly [mol C m-3] + E_res*)
    !> where E_res*=max(alpha_residue * E_depoly*, E_persistent,min).
    IF (d_residue > 0.0_wp .AND. enzyme_residue > 0.0_wp .AND. vmax_res_eff > 0.0_wp) THEN
      residue_carbon = d_residue * (km_depoly_eff + enzyme_residue) / &
        &             (vmax_res_eff * enzyme_residue)
    ELSE
      residue_carbon = eps8
    END IF

    !> A = qmax_min [mol C kg-1] * rho_bulk [kg m-3]
    qmax_a = qmax_org_min_local * bulk_dens_local
    !> B = f_C^SOM [-] * rho_org [kg SOM m-3], converted to [mol C m-3]
    qmax_b = carbon_per_dryweight_SOM * rho_bulk_org * 1000._wp / molar_mass_C

    IF (l_surface_layer) THEN
      !> Keep the first-layer physical constraint aligned with the mineral-volume treatment:
      !> woody litter stays above the mineral soil matrix, while soluable litter, polymeric litter,
      !> microbial biomass, DOM, and residue are the non-associated surface SOM that can be buried into
      !> the first layer.
      !> C0*_surface = C_sol* + C_poly* + C_D* + C_B* + C_res*
      c0_star = soluable_carbon + polymeric_carbon + dom_carbon + microbial_carbon + residue_carbon
    ELSE
      !> C0* = C_sol* + C_poly* + C_D* + C_B* + C_res*
      c0_star = soluable_carbon + polymeric_carbon + dom_carbon + microbial_carbon + residue_carbon
    END IF

    !> r_D = k_ads,D / k_des,D
    ratio_dom = k_adsorpt_dom * rtm_sorb_eff * rmm_sorb_eff / (k_desorpt_dom * rtm_desorb_eff * rmm_desorb_eff)
    !> r_R = k_ads,R / k_des,R
    ratio_res = k_adsorpt_det * rtm_sorb_eff * rmm_sorb_eff / (k_desorpt_det * rtm_desorb_eff * rmm_desorb_eff)
    l_particle_tendency = ABS(particle_poly_c) + ABS(particle_residue_c) + ABS(particle_dom_assoc_c) &
      &              + ABS(particle_residue_assoc_c) > 0.0_wp
    beta_dom = 0.0_wp
    beta_residue = 0.0_wp
    IF (l_particle_tendency) THEN
      IF (k_desorpt_dom * rtm_desorb_eff * rmm_desorb_eff > eps8) THEN
        beta_dom = particle_dom_assoc_c / (k_desorpt_dom * rtm_desorb_eff * rmm_desorb_eff)
      END IF
      IF (k_desorpt_det * rtm_desorb_eff * rmm_desorb_eff > eps8) THEN
        beta_residue = particle_residue_assoc_c / (k_desorpt_det * rtm_desorb_eff * rmm_desorb_eff)
      END IF
    END IF
    beta_total = beta_dom + beta_residue

    !> F = (r_D * C_D* + r_R * C_res*) / (1 + r_D * C_D* + r_R * C_res*)
    qmax_f = (ratio_dom * dom_carbon + ratio_res * residue_carbon) / &
      &      (1._wp + ratio_dom * dom_carbon + ratio_res * residue_carbon)

    !> Qmax* = (A * (1 - C0* / B)) / (1 + A * F / B)
    IF (qmax_b > eps8) THEN
      IF (l_particle_tendency) THEN
        qmax_org_layer = MAX(0.0_wp, qmax_a * &
          & (1._wp - (c0_star + (1._wp - qmax_f) * beta_total) / qmax_b) / &
          & (1._wp + qmax_a * qmax_f / qmax_b))
      ELSE
        qmax_org_layer = MAX(0.0_wp, qmax_a * (1._wp - c0_star / qmax_b) / (1._wp + qmax_a * qmax_f / qmax_b))
      END IF
    ELSE
      qmax_org_layer = eps8
    END IF

    sorption_denom = 1._wp + ratio_res * residue_carbon + ratio_dom * dom_carbon
    IF (sorption_denom > eps8 .AND. l_particle_tendency) THEN
      free_sorption_capacity = MAX(0.0_wp, (qmax_org_layer - beta_total) / sorption_denom)
      dom_assoc_carbon = MAX(0.0_wp, ratio_dom * dom_carbon * free_sorption_capacity + beta_dom)
      residue_assoc_carbon = MAX(0.0_wp, ratio_res * residue_carbon * free_sorption_capacity + beta_residue)
    ELSE IF (sorption_denom > eps8) THEN
      !> C_aDOM* = (r_D * C_D*) * Qmax* / (1 + r_R * C_res* + r_D * C_D*)
      dom_assoc_carbon     = ratio_dom * dom_carbon     * qmax_org_layer / sorption_denom
      !> C_aRes* = (r_R * C_res*) * Qmax* / (1 + r_R * C_res* + r_D * C_D*)
      residue_assoc_carbon = ratio_res * residue_carbon * qmax_org_layer / sorption_denom
    ELSE
      dom_assoc_carbon     = eps8
      residue_assoc_carbon = eps8
    END IF

    IF (l_surface_layer) THEN
      !> C_bulk* for the first mineral layer includes minerally associated SOM fully,
      !> and all other carbon pools as potentially buried surface SOM.
      bulk_carbon_layer = soluable_carbon + polymeric_carbon + microbial_carbon + dom_carbon &
                          + residue_carbon + dom_assoc_carbon + residue_assoc_carbon
    ELSE
      !> C_bulk* = C_wl* + C_sol* + C_poly* + C_D* + C_B* + C_res* + C_aDOM* + C_aRes*
      bulk_carbon_layer = woody_carbon + soluable_carbon + polymeric_carbon + microbial_carbon + dom_carbon &
        &                 + residue_carbon + dom_assoc_carbon + residue_assoc_carbon
    END IF
  END SUBROUTINE calc_analytical_ss_jsm_layer_state

  ! ======================================================================================================= !
  !> Solve the coupled microbial-biomass and DOM steady state with a prescribed microbial
  !> transport tendency. Carbon source excludes the outgoing DOC term lambda_doc * C_DOM.
  SUBROUTINE calc_microbial_dom_ss_transport( &
    & carbon_source, lambda_doc, microbial_transport_c, microbial_cue, microbial_mortality, &
    & km_uptake, vmax_uptake, microbial_carbon, dom_carbon, l_solution)

    REAL(wp), INTENT(in)  :: carbon_source          !< net C source excluding DOC export [mol C m-3 s-1]
    REAL(wp), INTENT(in)  :: lambda_doc             !< first-order downward DOC transport coefficient [s-1]
    REAL(wp), INTENT(in)  :: microbial_transport_c  !< net microbial C transport tendency [mol C m-3 s-1]
    REAL(wp), INTENT(in)  :: microbial_cue          !< microbial carbon-use efficiency [-]
    REAL(wp), INTENT(in)  :: microbial_mortality    !< microbial mortality rate coefficient [s-1]
    REAL(wp), INTENT(in)  :: km_uptake              !< modified microbial uptake half-saturation [mol C m-3]
    REAL(wp), INTENT(in)  :: vmax_uptake            !< modified maximum microbial uptake rate [s-1]
    REAL(wp), INTENT(out) :: microbial_carbon       !< positive microbial steady-state C density [mol C m-3]
    REAL(wp), INTENT(out) :: dom_carbon              !< positive DOM steady-state C density [mol C m-3]
    LOGICAL,  INTENT(out) :: l_solution              !< TRUE when a physically admissible root was found

    REAL(wp) :: biomass_intercept       !< biomass at zero DOM implied by the layer C balance [mol C m-3]
    REAL(wp) :: doc_slope               !< biomass decrease per unit DOM from DOC export [-]
    REAL(wp) :: coef_a                  !< quadratic coefficient of DOM^2 [s-1]
    REAL(wp) :: coef_b                  !< quadratic coefficient of DOM [mol C m-3 s-1]
    REAL(wp) :: coef_c                  !< quadratic constant [(mol C m-3)^2 s-1]
    REAL(wp) :: discriminant            !< quadratic discriminant [(mol C m-3 s-1)^2]
    REAL(wp) :: dom_candidates(2)       !< candidate roots for DOM C density [mol C m-3]
    REAL(wp) :: biomass_candidate       !< microbial C associated with a candidate DOM root [mol C m-3]
    REAL(wp) :: microbial_loss          !< net microbial loss after transport [mol C m-3 s-1]
    REAL(wp) :: microbial_uptake        !< candidate gross microbial DOC uptake [mol C m-3 s-1]
    REAL(wp) :: residual                !< microbial steady-state equation residual [mol C m-3 s-1]
    REAL(wp) :: residual_scale          !< magnitude used to normalize the residual [mol C m-3 s-1]
    INTEGER  :: iroot                   !< index of the two quadratic roots [-]

    microbial_carbon = 0.0_wp
    dom_carbon = 0.0_wp
    l_solution = .FALSE.
    IF (microbial_cue <= 0.0_wp .OR. microbial_cue >= 1.0_wp .OR. microbial_mortality <= 0.0_wp) RETURN

    !> Combining the whole-layer C balance with
    !>   CUE * U - k_B * C_B + T_B = 0
    !> gives C_B = (T_B + CUE * (S - lambda * C_DOM)) / (k_B * (1 - CUE)).
    biomass_intercept = (microbial_transport_c + microbial_cue * carbon_source) / &
      &                 (microbial_mortality * (1.0_wp - microbial_cue))
    doc_slope = microbial_cue * MAX(0.0_wp, lambda_doc) / &
      &         (microbial_mortality * (1.0_wp - microbial_cue))

    coef_a = -doc_slope * (microbial_cue * vmax_uptake - microbial_mortality)
    coef_b = (microbial_cue * vmax_uptake - microbial_mortality) * biomass_intercept &
      &      + microbial_transport_c + km_uptake * microbial_mortality * doc_slope
    coef_c = -km_uptake * (microbial_mortality * biomass_intercept - microbial_transport_c)

    dom_candidates(:) = -1.0_wp
    IF (ABS(coef_a) <= TINY(1.0_wp)) THEN
      IF (ABS(coef_b) > TINY(1.0_wp)) dom_candidates(1) = -coef_c / coef_b
    ELSE
      discriminant = coef_b * coef_b - 4.0_wp * coef_a * coef_c
      IF (discriminant >= 0.0_wp) THEN
        discriminant = SQRT(discriminant)
        dom_candidates(1) = (-coef_b + discriminant) / (2.0_wp * coef_a)
        dom_candidates(2) = (-coef_b - discriminant) / (2.0_wp * coef_a)
      END IF
    END IF

    DO iroot = 1, 2
      IF (dom_candidates(iroot) < 0.0_wp) CYCLE
      biomass_candidate = biomass_intercept - doc_slope * dom_candidates(iroot)
      microbial_loss = microbial_mortality * biomass_candidate - microbial_transport_c
      IF (km_uptake + dom_candidates(iroot) > TINY(1.0_wp)) THEN
        microbial_uptake = vmax_uptake * biomass_candidate * dom_candidates(iroot) / &
          &                (km_uptake + dom_candidates(iroot))
      ELSE
        microbial_uptake = 0.0_wp
      END IF
      residual = microbial_cue * microbial_uptake - microbial_loss
      residual_scale = MAX(TINY(1.0_wp), ABS(microbial_cue * microbial_uptake), ABS(microbial_loss))
      ! Exclude the algebraic B=0 root that appears when T_B=0; it does not
      ! represent the positive microbial steady state selected by the legacy solution.
      IF (biomass_candidate > eps8 .AND. microbial_loss >= -eps8 .AND. &
        & ABS(residual) <= 1000.0_wp * EPSILON(1.0_wp) * residual_scale) THEN
        microbial_carbon = biomass_candidate
        dom_carbon = dom_candidates(iroot)
        l_solution = .TRUE.
        EXIT
      END IF
    END DO

  END SUBROUTINE calc_microbial_dom_ss_transport

  ! ======================================================================================================= !
  !> calculate analytical steady-state carbon pools for the Jena Soil Model
  !> The current C-only analytical draft uses
  !>   C_wl*   = I_wl tau_wl
  !>   C_sol*  = I_sol tau_sol
  !>   C_D*    = (k_B * K_m,U) / (CUE * vmax_U - k_B)
  !> and the layer-wise DOC transport closure
  !>   C_B,sl*  = CUE_sl / (k_B * (1 - CUE_sl)) * C_eff_input_sl
  !>   C_eff_input_sl = (eta_DOM^sol * I_sol,sl + I_poly,sl + eta_poly^wl * I_wl,sl
  !>                - lambda_sl * C_D,sl* + lambda_(sl-1) * C_D,sl-1*)
  !> together with the associated-pool and sorption-capacity closure
  !>   C_aDOM* = ((k_ads,D [s-1] * C_D* [mol C m-3]) / k_des,D [s-1]) * Qmax* [mol C m-3]
  !>             / (1 + (k_ads,R / k_des,R) * C_res* + (k_ads,D / k_des,D) * C_D*)
  !>   C_aRes* = ((k_ads,R [s-1] * C_res* [mol C m-3]) / k_des,R [s-1]) * Qmax* [mol C m-3]
  !>             / (1 + (k_ads,R / k_des,R) * C_res* + (k_ads,D / k_des,D) * C_D*)
  !>   Qmax*   = (A * (1 - C0* / B)) / (1 + A * F / B)
  !>   C0*     = C_sol* + C_poly* + C_D* + C_B* + C_res*
  !>   A       = qmax_min [mol C kg-1] * rho_bulk [kg m-3]
  !>   B       = f_C^SOM [-] * rho_org [kg SOM m-3], converted below to [mol C m-3]
  !>   F       = (r_D * C_D* + r_R * C_res*) / (1 + r_D * C_D* + r_R * C_res*)
  !> before solving the polymeric-litter and residue pools from the enzyme-allocation closure.
  !>   C_poly*=(I_poly + eta_poly^wl * I_wl)/(vmax_poly * enzyme_frac_poly * C_B* )(Km_depoly + enzyme_frac_poly * C_B*)
  !>   C_res*=((1-f_soluable_residue) * k_B * C_B*)/(vmax_res * enzyme_frac_residue * C_B*)(Km_depoly + enzyme_frac_residue * C_B*)
  !>
  !> Unit convention of the explicit analytical forcing fields used here:
  !>   I_sol, I_poly, I_wl        [mol C m-3 s-1]
  !>   lambda_doc                 [s-1]
  !>   tau_*                      [s]
  !>   k_B, vmax_*, k_ads, k_des  [s-1]
  !>   K_m,*                      [mol C m-3]
  !>   qmax_org_fine_particle     [mol C kg-1]
  !>   rho_bulk                   [kg m-3]
  !>   rtm_*, rmm_*, enzyme_frac* [-]
  !>   C_*, Qmax*                 [mol C m-3]
  !> Litter-formation and DOC-transport inputs must be mean rates, not timestep sums.
  SUBROUTINE calc_analytical_ss_jsm( &
    & nproma, &
    & nsoil_w, &
    & nblks, &
    & num_sl_above_bedrock, &
    & qmax_org_fine_particle, &
    & soil_lay_depth_center_sl, &
    & soil_lay_depth_ubound_sl, &
    & formation_soluable_litter_c, &
    & formation_polymeric_litter_c, &
    & formation_woody_litter_c, &
    & frac_wtr_transp_down_sl, &
    & enzyme_frac_ap_forcing, &
    & silt_sl, &
    & clay_sl, &
    & bulk_dens_sl, &
    & bulk_carbon_max, &
    & qmax_org, &
    & rtm_depolymerisation, &
    & rmm_depolymerisation, &
    & rtm_mic_uptake, &
    & rmm_mic_uptake, &
    & rtm_sorption, &
    & rmm_sorption, &
    & rtm_desorption, &
    & rmm_desorption, &
    & rtm_hsc, &
    & rmm_hsc, &
    & sb_pool_mt_domain, &
    & surface_burial_frac, &
    & enzyme_frac_poly, &
    & enzyme_frac_residue, &
    & enzyme_frac_poly_c, &
    & enzyme_frac_poly_n, &
    & enzyme_frac_poly_p, &
    & enzyme_frac_AP, &
    & l_constraint_debug, &
    & qmax_org_min_forcing, &
    & transport_soluable_c, &
    & transport_poly_c, &
    & transport_woody_c, &
    & transport_dom_c, &
    & transport_microbial_c, &
    & transport_residue_c, &
    & transport_dom_assoc_c, &
    & transport_residue_assoc_c)

    USE mo_sb_util,                       ONLY: calc_qmax_texture
    USE mo_sb_constants
    USE mo_jsb_physical_constants, ONLY: molar_mass_C

    INTEGER,  INTENT(in)    :: nproma                         !< horizontal points per block [-]
    INTEGER,  INTENT(in)    :: nsoil_w                        !< number of soil layers [-]
    INTEGER,  INTENT(in)    :: nblks                          !< number of horizontal blocks [-]
    REAL(wp), INTENT(in)    :: num_sl_above_bedrock(:,:)      !< active soil-layer count per point [-]
    REAL(wp), INTENT(in)    :: qmax_org_fine_particle(:,:,:)  !< fine-particle sorption capacity [mol C kg-1]
    REAL(wp), INTENT(in)    :: soil_lay_depth_center_sl(:,:,:)!< soil-layer midpoint depth [m]
    REAL(wp), INTENT(in)    :: soil_lay_depth_ubound_sl(:,:,:) !< soil-layer lower-bound depth [m]
    REAL(wp), INTENT(in)    :: formation_soluable_litter_c(:,:,:) !< soluble-litter input [mol C m-3 s-1]
    REAL(wp), INTENT(in)    :: formation_polymeric_litter_c(:,:,:) !< polymeric-litter input [mol C m-3 s-1]
    REAL(wp), INTENT(in)    :: formation_woody_litter_c(:,:,:) !< woody-litter input [mol C m-3 s-1]
    REAL(wp), INTENT(in)    :: frac_wtr_transp_down_sl(:,:,:) !< downward DOC transport coefficient [s-1]
    REAL(wp), INTENT(in)    :: enzyme_frac_ap_forcing(:,:,:)  !< mean phosphatase enzyme fraction [-]
    REAL(wp), INTENT(in)    :: silt_sl(:,:,:)                 !< soil silt fraction [-]
    REAL(wp), INTENT(in)    :: clay_sl(:,:,:)                 !< soil clay fraction [-]
    REAL(wp), INTENT(in)    :: bulk_dens_sl(:,:,:)            !< mineral soil bulk density [kg m-3]
    REAL(wp), INTENT(inout) :: bulk_carbon_max(:,:,:)         !< maximum admissible soil C density [mol C m-3]
    REAL(wp), INTENT(inout) :: qmax_org(:,:,:)                 !< analytical organic sorption capacity [mol C m-3]
    REAL(wp), INTENT(in)    :: rtm_depolymerisation(:,:,:)    !< temperature modifier of depolymerisation [-]
    REAL(wp), INTENT(in)    :: rmm_depolymerisation(:,:,:)    !< moisture modifier of depolymerisation [-]
    REAL(wp), INTENT(in)    :: rtm_mic_uptake(:,:,:)          !< temperature modifier of microbial uptake [-]
    REAL(wp), INTENT(in)    :: rmm_mic_uptake(:,:,:)          !< moisture modifier of microbial uptake [-]
    REAL(wp), INTENT(in)    :: rtm_sorption(:,:,:)            !< temperature modifier of sorption [-]
    REAL(wp), INTENT(in)    :: rmm_sorption(:,:,:)            !< moisture modifier of sorption [-]
    REAL(wp), INTENT(in)    :: rtm_desorption(:,:,:)          !< temperature modifier of desorption [-]
    REAL(wp), INTENT(in)    :: rmm_desorption(:,:,:)          !< moisture modifier of desorption [-]
    REAL(wp), INTENT(in)    :: rtm_hsc(:,:,:)                 !< temperature modifier of half-saturation constants [-]
    REAL(wp), INTENT(in)    :: rmm_hsc(:,:,:)                 !< moisture modifier of half-saturation constants [-]
    REAL(wp), INTENT(inout) :: sb_pool_mt_domain(:,:,:,:,:)   !< pool by element, point, layer and block [mol m-3]
    REAL(wp), INTENT(inout) :: surface_burial_frac(:,:,:)     !< fraction of surface SOM buried in layer one [-]
    REAL(wp), INTENT(out)   :: enzyme_frac_poly(:,:,:)        !< enzyme allocation to polymeric litter [-]
    REAL(wp), INTENT(out)   :: enzyme_frac_residue(:,:,:)     !< enzyme allocation to microbial residue [-]
    REAL(wp), INTENT(out)   :: enzyme_frac_poly_c(:,:,:)      !< C-based polymeric enzyme allocation [-]
    REAL(wp), INTENT(out)   :: enzyme_frac_poly_n(:,:,:)      !< N-based polymeric enzyme allocation [-]
    REAL(wp), INTENT(out)   :: enzyme_frac_poly_p(:,:,:)      !< P-based polymeric enzyme allocation [-]
    REAL(wp), INTENT(out)   :: enzyme_frac_AP(:,:,:)          !< phosphatase fraction of total enzyme [-]
    LOGICAL,  INTENT(in), OPTIONAL :: l_constraint_debug      !< emit physical-constraint diagnostics when TRUE
    REAL(wp), INTENT(in), OPTIONAL :: qmax_org_min_forcing(:,:,:) !< prescribed capacity [mol C kg-1]
    REAL(wp), INTENT(in), OPTIONAL :: transport_soluable_c(:,:,:) !< soluble-C transport tendency [mol C m-3 s-1]
    REAL(wp), INTENT(in), OPTIONAL :: transport_poly_c(:,:,:) !< polymeric-C transport tendency [mol C m-3 s-1]
    REAL(wp), INTENT(in), OPTIONAL :: transport_woody_c(:,:,:) !< woody-C transport tendency [mol C m-3 s-1]
    REAL(wp), INTENT(in), OPTIONAL :: transport_dom_c(:,:,:)  !< DOM-C transport tendency [mol C m-3 s-1]
    REAL(wp), INTENT(in), OPTIONAL :: transport_microbial_c(:,:,:) !< microbial-C tendency [mol C m-3 s-1]
    REAL(wp), INTENT(in), OPTIONAL :: transport_residue_c(:,:,:) !< residue-C tendency [mol C m-3 s-1]
    REAL(wp), INTENT(in), OPTIONAL :: transport_dom_assoc_c(:,:,:) !< associated-DOM tendency [mol C m-3 s-1]
    REAL(wp), INTENT(in), OPTIONAL :: transport_residue_assoc_c(:,:,:) !< associated-residue tendency [mol C m-3 s-1]

    REAL(wp) :: input_sol                       !< local soluble-litter input [mol C m-3 s-1]
    REAL(wp) :: input_poly                      !< local polymeric-litter input [mol C m-3 s-1]
    REAL(wp) :: input_wood                      !< local woody-litter input [mol C m-3 s-1]
    REAL(wp) :: eff_input_c                     !< net C input available to microbial growth [mol C m-3 s-1]
    REAL(wp) :: carbon_source                   !< layer C source excluding downward DOC loss [mol C m-3 s-1]
    REAL(wp) :: transport_soluable_local        !< local soluble-C transport tendency [mol C m-3 s-1]
    REAL(wp) :: transport_poly_local            !< local polymeric-C transport tendency [mol C m-3 s-1]
    REAL(wp) :: transport_woody_local           !< local woody-C transport tendency [mol C m-3 s-1]
    REAL(wp) :: transport_dom_local             !< local DOM-C transport tendency [mol C m-3 s-1]
    REAL(wp) :: transport_microbial_local       !< local microbial-C transport tendency [mol C m-3 s-1]
    REAL(wp) :: transport_residue_local         !< local residue-C transport tendency [mol C m-3 s-1]
    REAL(wp) :: transport_dom_assoc_local       !< local associated-DOM tendency [mol C m-3 s-1]
    REAL(wp) :: transport_residue_assoc_local   !< local associated-residue tendency [mol C m-3 s-1]
    REAL(wp) :: alpha_poly                      !< enzyme allocation to polymeric litter [-]
    REAL(wp) :: alpha_residue                   !< enzyme allocation to microbial residue [-]
    REAL(wp) :: denominator                     !< positive microbial uptake-minus-mortality rate [s-1]
    REAL(wp) :: km_uptake_eff                   !< environmentally modified uptake half-saturation [mol C m-3]
    REAL(wp) :: vmax_uptake_eff                 !< environmentally modified uptake rate [s-1]
    REAL(wp) :: km_depoly_eff                   !< environmentally modified depolymerisation half-saturation [mol C m-3]
    REAL(wp) :: vmax_poly_eff                   !< modified polymeric-litter depolymerisation rate [s-1]
    REAL(wp) :: vmax_res_eff                    !< modified residue depolymerisation rate [s-1]
    REAL(wp) :: qmax_org_min_local              !< local texture-dependent sorption capacity [mol C kg-1]
    REAL(wp) :: qmax_b                          !< C density of pure organic matter [mol C m-3]
    REAL(wp) :: lambda_doc                      !< downward DOC transport coefficient of this layer [s-1]
    REAL(wp) :: lambda_doc_above                !< downward DOC transport coefficient of layer above [s-1]
    REAL(wp) :: rtm_dep_eff                    !< valid local temperature modifier of depolymerisation [-]
    REAL(wp) :: rmm_dep_eff                    !< valid local moisture modifier of depolymerisation [-]
    REAL(wp) :: rtm_upt_eff                    !< valid local temperature modifier of microbial uptake [-]
    REAL(wp) :: rmm_upt_eff                    !< valid local moisture modifier of microbial uptake [-]
    REAL(wp) :: rtm_sorb_eff                   !< valid local temperature modifier of sorption [-]
    REAL(wp) :: rmm_sorb_eff                   !< valid local moisture modifier of sorption [-]
    REAL(wp) :: rtm_desorb_eff                 !< valid local temperature modifier of desorption [-]
    REAL(wp) :: rmm_desorb_eff                 !< valid local moisture modifier of desorption [-]
    REAL(wp) :: rtm_hsc_eff                    !< valid local temperature modifier of half-saturation constants [-]
    REAL(wp) :: rmm_hsc_eff                    !< valid local moisture modifier of half-saturation constants [-]
    REAL(wp) :: bulk_dens_local                !< local mineral bulk density [kg m-3]
    REAL(wp) :: silt_local                     !< local silt fraction [-]
    REAL(wp) :: clay_local                     !< local clay fraction [-]
    REAL(wp) :: microbial_c_bio                !< unconstrained biological microbial C solution [mol C m-3]
    REAL(wp) :: microbial_c_bio_unclamped      !< microbial C solution before applying the numerical minimum [mol C m-3]
    REAL(wp) :: polymeric_carbon               !< analytical polymeric-litter C density [mol C m-3]
    REAL(wp) :: residue_carbon                 !< analytical microbial-residue C density [mol C m-3]
    REAL(wp) :: dom_assoc_carbon               !< analytical mineral-associated DOM C [mol C m-3]
    REAL(wp) :: residue_assoc_carbon           !< analytical mineral-associated residue C [mol C m-3]
    REAL(wp) :: bulk_carbon_bio                !< bulk C from the current biological solution [mol C m-3]
    REAL(wp) :: bulk_carbon_unconstrained      !< bulk C before applying physical constraints [mol C m-3]
    REAL(wp) :: bulk_carbon_mid                !< bulk C at the current enzyme-scale midpoint [mol C m-3]
    REAL(wp) :: min_volume_min_layer           !< prescribed minimum mineral volume fraction [-]
    REAL(wp) :: cb_bracket_tol                 !< enzyme-scale bisection bracket tolerance [-]
    REAL(wp) :: surface_bulk_nonsorbed         !< non-associated surface C eligible for burial [mol C m-3]
    REAL(wp) :: surface_bulk_assoc             !< fully buried associated surface C [mol C m-3]
    REAL(wp) :: qmax_org_target                !< capacity implied by minimum mineral volume [mol C m-3]
    REAL(wp) :: assoc_scale                    !< multiplier mapping unconstrained to target associated pools [-]
    REAL(wp) :: enzyme_frac_total_scale        !< current multiplier of total enzyme investment [-]
    REAL(wp) :: enzyme_scale_lower             !< lower enzyme-scale bisection bound [-]
    REAL(wp) :: enzyme_scale_upper             !< upper enzyme-scale bisection bound [-]
    REAL(wp) :: enzyme_scale_mid               !< current enzyme-scale bisection midpoint [-]
    REAL(wp) :: bulk_carbon_scale_upper        !< bulk C evaluated at upper enzyme-scale bound [mol C m-3]
    REAL(wp) :: max_enzyme_frac_total_scale    !< largest enzyme multiplier used in this solve [-]
    REAL(wp), PARAMETER :: enzyme_frac_total_phys_scale = 10.0_wp !< factor used to expand enzyme bracket [-]
    LOGICAL  :: l_surface_layer                !< TRUE while solving the first mineral layer
    LOGICAL  :: l_transport_tendencies_present !< TRUE when all optional pool tendencies were supplied
    LOGICAL  :: l_microbial_solution           !< TRUE when an admissible positive microbial solution exists
    LOGICAL  :: l_microbial_dormant            !< TRUE when microbial C is held at its minimum
    INTEGER  :: iblk                           !< horizontal block index [-]
    INTEGER  :: is                             !< soil-layer index [-]
    INTEGER  :: ic                             !< horizontal point index [-]
    INTEGER  :: iter_phys_constraint           !< current physical-constraint bisection iteration [-]
    INTEGER  :: constraint_iterations          !< iterations used for the current constrained layer [-]
    INTEGER  :: n_surface_constraints          !< number of surface burial constraints applied [-]
    INTEGER  :: n_enzyme_constraints           !< number of deeper-layer enzyme constraints applied [-]
    INTEGER, PARAMETER :: max_iter_phys_constraint = 200 !< maximum enzyme-scale bisection iterations [-]
    LOGICAL  :: l_print_constraint_debug       !< local switch for constraint diagnostics
    CHARACTER(len=1024) :: message_text         !< formatted model-log message
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_analytical_ss_jsm' !< log context
    ! This flag describes this individual analytical solve. It is independent of
    ! the model option that enables the outer one-pass transport correction.
    l_transport_tendencies_present = PRESENT(transport_soluable_c) .AND. PRESENT(transport_poly_c) .AND. &
      & PRESENT(transport_woody_c) .AND. PRESENT(transport_dom_c) .AND. PRESENT(transport_microbial_c) .AND. &
      & PRESENT(transport_residue_c) .AND. PRESENT(transport_dom_assoc_c) .AND. &
      & PRESENT(transport_residue_assoc_c)
    IF (ANY([PRESENT(transport_soluable_c), PRESENT(transport_poly_c), PRESENT(transport_woody_c), &
      & PRESENT(transport_dom_c), PRESENT(transport_microbial_c), PRESENT(transport_residue_c), &
      & PRESENT(transport_dom_assoc_c), PRESENT(transport_residue_assoc_c)]) .AND. &
      & .NOT. l_transport_tendencies_present) THEN
      CALL finish(TRIM(routine), 'All eight analytically solved transport tendencies must be supplied together.')
    END IF
    l_print_constraint_debug = .FALSE.
    IF (PRESENT(l_constraint_debug)) l_print_constraint_debug = l_constraint_debug
    n_surface_constraints = 0
    n_enzyme_constraints = 0
    max_enzyme_frac_total_scale = 1.0_wp
    !> init out var
    !> The current analytical steady-state solution only solves for the carbon pools, so we can safely
    !> initialise the other elements with zeros here. The actual NetCDF reader will need to fill
    sb_pool_mt_domain(ix_soluable_litter, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_polymeric_litter, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_woody_litter, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_microbial, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_dom, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_dom_assoc, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_residue, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_residue_assoc, ixC, :, :, :) = 0.0_wp
    qmax_org(:,:,:) = 0.0_wp
    surface_burial_frac(:,:,:) = 1.0_wp
    enzyme_frac_poly(:,:,:) = 0.0_wp
    enzyme_frac_residue(:,:,:) = 0.0_wp
    enzyme_frac_poly_c(:,:,:) = 0.0_wp
    enzyme_frac_poly_n(:,:,:) = 0.0_wp
    enzyme_frac_poly_p(:,:,:) = 0.0_wp
    enzyme_frac_AP(:,:,:) = 0.0_wp

    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))

          input_sol  = MAX(0.0_wp, formation_soluable_litter_c(ic, is, iblk))
          input_poly = MAX(0.0_wp, formation_polymeric_litter_c(ic, is, iblk))
          input_wood = MAX(0.0_wp, formation_woody_litter_c(ic, is, iblk))
          transport_soluable_local = 0.0_wp
          transport_poly_local = 0.0_wp
          transport_woody_local = 0.0_wp
          transport_dom_local = 0.0_wp
          transport_microbial_local = 0.0_wp
          transport_residue_local = 0.0_wp
          transport_dom_assoc_local = 0.0_wp
          transport_residue_assoc_local = 0.0_wp
          IF (l_transport_tendencies_present) THEN
            transport_soluable_local = transport_soluable_c(ic,is,iblk)
            transport_poly_local = transport_poly_c(ic,is,iblk)
            transport_woody_local = transport_woody_c(ic,is,iblk)
            transport_dom_local = transport_dom_c(ic,is,iblk)
            transport_microbial_local = transport_microbial_c(ic,is,iblk)
            transport_residue_local = transport_residue_c(ic,is,iblk)
            transport_dom_assoc_local = transport_dom_assoc_c(ic,is,iblk)
            transport_residue_assoc_local = transport_residue_assoc_c(ic,is,iblk)
          END IF

          ! Soluble and woody litter have linear analytical lifetimes. Their net
          ! particle tendencies therefore enter as additional signed source rates.
          input_sol = MAX(0.0_wp, input_sol + transport_soluable_local)
          input_wood = MAX(0.0_wp, input_wood + transport_woody_local)

          bulk_dens_local = bulk_dens_sl(ic, is, iblk)
          silt_local = silt_sl(ic, is, iblk)
          clay_local = clay_sl(ic, is, iblk)

          rtm_dep_eff    = rtm_depolymerisation(ic, is, iblk)
          rmm_dep_eff    = rmm_depolymerisation(ic, is, iblk)
          rtm_upt_eff    = rtm_mic_uptake(ic, is, iblk)
          rmm_upt_eff    = rmm_mic_uptake(ic, is, iblk)
          rtm_sorb_eff   = rtm_sorption(ic, is, iblk)
          rmm_sorb_eff   = rmm_sorption(ic, is, iblk)
          rtm_desorb_eff = rtm_desorption(ic, is, iblk)
          rmm_desorb_eff = rmm_desorption(ic, is, iblk)
          rtm_hsc_eff    = rtm_hsc(ic, is, iblk)
          rmm_hsc_eff    = rmm_hsc(ic, is, iblk)

          rtm_dep_eff    = MERGE(rtm_dep_eff,    1.0_wp, rtm_dep_eff    > eps8)
          rmm_dep_eff    = MERGE(rmm_dep_eff,    1.0_wp, rmm_dep_eff    > eps8)
          rtm_upt_eff    = MERGE(rtm_upt_eff,    1.0_wp, rtm_upt_eff    > eps8)
          rmm_upt_eff    = MERGE(rmm_upt_eff,    1.0_wp, rmm_upt_eff    > eps8)
          rtm_sorb_eff   = MERGE(rtm_sorb_eff,   1.0_wp, rtm_sorb_eff   > eps8)
          rmm_sorb_eff   = MERGE(rmm_sorb_eff,   1.0_wp, rmm_sorb_eff   > eps8)
          rtm_desorb_eff = MERGE(rtm_desorb_eff, 1.0_wp, rtm_desorb_eff > eps8)
          rmm_desorb_eff = MERGE(rmm_desorb_eff, 1.0_wp, rmm_desorb_eff > eps8)
          rtm_hsc_eff    = MERGE(rtm_hsc_eff,    1.0_wp, rtm_hsc_eff    > eps8)
          rmm_hsc_eff    = MERGE(rmm_hsc_eff,    1.0_wp, rmm_hsc_eff    > eps8)

          km_uptake_eff = km_mic_uptake * rtm_hsc_eff * rmm_hsc_eff
          vmax_uptake_eff = vmax_mic_uptake * rtm_upt_eff * rmm_upt_eff
          km_depoly_eff = km_depolymerisation * rtm_hsc_eff * rmm_hsc_eff
          vmax_poly_eff = vmax_depolymerisation_litter * rtm_dep_eff * rmm_dep_eff
          vmax_res_eff = vmax_depolymerisation_residue * rtm_dep_eff * rmm_dep_eff
          l_surface_layer = (is == 1)

          !> C_wl*   = I_wl [mol C m-3 s-1] * tau_wl [s]
          sb_pool_mt_domain(ix_woody_litter, ixC, ic, is, iblk)    = input_wood * tau_woody_litter
          !> C_sol*  = I_sol [mol C m-3 s-1] * tau_sol [s]
          sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) = input_sol * tau_soluable_litter

          lambda_doc = MAX(0.0_wp, frac_wtr_transp_down_sl(ic, is, iblk))
          IF (is > 1) THEN
            lambda_doc_above = MAX(0.0_wp, frac_wtr_transp_down_sl(ic, is - 1, iblk))
          ELSE
            lambda_doc_above = 0.0_wp
          END IF

          carbon_source = 0.0_wp
          eff_input_c = 0.0_wp
          l_microbial_solution = .FALSE.
          IF (l_transport_tendencies_present) THEN
            carbon_source = fc_soluable2dom * input_sol + input_poly + fc_woody2polymeric * input_wood &
              &             + transport_poly_local + transport_residue_local &
              &             + transport_dom_local + transport_dom_assoc_local + transport_residue_assoc_local
            IF (is > 1) carbon_source = carbon_source + lambda_doc_above * &
              & sb_pool_mt_domain(ix_dom, ixC, ic, is - 1, iblk)
            CALL calc_microbial_dom_ss_transport( &
              & carbon_source, lambda_doc, transport_microbial_local, microbial_cue_max, k_mort_microbial, &
              & km_uptake_eff, vmax_uptake_eff, microbial_c_bio, &
              & sb_pool_mt_domain(ix_dom, ixC, ic,is,iblk), l_microbial_solution)
            IF (.NOT. l_microbial_solution) THEN
              microbial_c_bio_unclamped = 0.0_wp
              microbial_c_bio = min_mic_biomass
              l_microbial_dormant = .TRUE.
              sb_pool_mt_domain(ix_dom, ixC, ic,is,iblk) = 0.0_wp
            ELSE
              microbial_c_bio_unclamped = microbial_c_bio
              l_microbial_dormant = microbial_c_bio <= min_mic_biomass
              microbial_c_bio = MAX(min_mic_biomass, microbial_c_bio)
            END IF
          ELSE
            denominator = microbial_cue_max * vmax_uptake_eff - k_mort_microbial
            IF (denominator <= eps8) THEN
              CALL finish(TRIM(routine), &
                & 'Analytical JSM initialisation requires modified microbial uptake to exceed microbial mortality.')
            END IF
            sb_pool_mt_domain(ix_dom, ixC, ic,is,iblk) = k_mort_microbial * km_uptake_eff / denominator
            eff_input_c = fc_soluable2dom * input_sol + input_poly + fc_woody2polymeric * input_wood &
              &            - lambda_doc * sb_pool_mt_domain(ix_dom, ixC, ic,is,iblk)
            IF (is > 1) THEN
              eff_input_c = eff_input_c + lambda_doc_above * sb_pool_mt_domain(ix_dom, ixC, ic,is - 1,iblk)
            END IF
            IF (eff_input_c > 0.0_wp) THEN
              microbial_c_bio = microbial_cue_max / &
                & (k_mort_microbial * (1._wp - microbial_cue_max)) * eff_input_c
              microbial_c_bio_unclamped = microbial_c_bio
              l_microbial_solution = .TRUE.
              l_microbial_dormant = microbial_c_bio <= min_mic_biomass
              microbial_c_bio = MAX(min_mic_biomass, microbial_c_bio)
            ELSE
              microbial_c_bio_unclamped = 0.0_wp
              microbial_c_bio = min_mic_biomass
              l_microbial_dormant = .TRUE.
            END IF
          END IF

          IF (l_microbial_dormant .AND. l_print_constraint_debug) THEN
            WRITE(message_text,*) 'DEBUG analytical microbial dormancy activated at ic/is/iblk = ', &
              & ic, is, iblk, &
              & '; transport-coupled solve, microbial solution available = ', &
              & l_transport_tendencies_present, l_microbial_solution, &
              & '; C_mic_raw, C_mic_min, C_mic_used [mol C m-3] = ', &
              & microbial_c_bio_unclamped, min_mic_biomass, microbial_c_bio, &
              & '; effective/local C source, DOM [mol C m-3 s-1; mol C m-3] = ', &
              & MERGE(carbon_source, eff_input_c, l_transport_tendencies_present), &
              & sb_pool_mt_domain(ix_dom, ixC, ic,is,iblk)
            CALL message(TRIM(routine), TRIM(message_text))
          END IF

          enzyme_frac_AP(ic, is, iblk) = MIN(MAX(enzyme_frac_ap_forcing(ic, is, iblk), &
            & min_enzyme_fraction), max_enzyme_fraction)
          IF (PRESENT(qmax_org_min_forcing)) THEN
            qmax_org_min_local = qmax_org_min_forcing(ic, is, iblk)
          ELSE
            qmax_org_min_local = calc_qmax_texture(qmax_org_fine_particle(ic, is, iblk), silt_local, clay_local)
            qmax_org_min_local = qmax_org_min_local &
              & * EXP(-is / (INT(num_sl_above_bedrock(ic, iblk)) * k_sl_conv_moc_corr)) * 0.5_wp
            IF (soil_lay_depth_ubound_sl(ic, is, iblk) > z_deep_soil_moc_corr) THEN
              qmax_org_min_local = qmax_org_min_local &
                & * MAX(0.0001_wp, 10._wp**(-soil_lay_depth_ubound_sl(ic, is, iblk)))
            END IF
          END IF

          !> Steady-state-only physical correction:
          !> phi_min^min(1) = 0.05
          !> phi_min^min(k) = phi_min^min(k-1) + 0.01 * (z_mid,k - z_mid,k-1)
          !> C_bulk^max(k) = B * (1 - phi_min^min(k))
          !> For the first mineral layer, if C_bulk^bio > C_bulk^max, solve the geometric
          !> surface-soil partition directly with
          !>   f_burial = (C_bulk^max - (C_aDOM* + C_aRes*))
          !>              / (C_sol* + C_poly* + C_B* + C_D* + C_res*)
          !> so only the buried fraction of the non-woody, non-associated surface SOM occupies soil volume.
          !> For deeper layers, keep the current bisection search on C_B.
          qmax_b = carbon_per_dryweight_SOM * rho_bulk_org * 1000._wp / molar_mass_C
          bulk_carbon_max(ic, is, iblk) = MAX(0.0_wp, bulk_carbon_max(ic, is, iblk))
          min_volume_min_layer = MAX(0.0_wp, 1.0_wp - bulk_carbon_max(ic, is, iblk) / qmax_b)
          enzyme_frac_total_scale = 1.0_wp

          CALL calc_analytical_ss_jsm_layer_state( &
            & l_surface_layer, &
            & input_sol, &
            & input_poly, &
            & input_wood, &
            & sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk), &
            & microbial_c_bio, &
            & l_microbial_dormant, &
            & enzyme_frac_AP(ic, is, iblk), &
            & enzyme_frac_total_scale, &
            & transport_poly_local, transport_residue_local, &
            & transport_dom_assoc_local, transport_residue_assoc_local, &
            & km_depoly_eff, &
            & vmax_poly_eff, &
            & vmax_res_eff, &
            & qmax_org_min_local, &
            & bulk_dens_local, &
            & rtm_sorb_eff, &
            & rmm_sorb_eff, &
            & rtm_desorb_eff, &
            & rmm_desorb_eff, &
            & alpha_poly, &
            & alpha_residue, &
            & polymeric_carbon, &
            & residue_carbon, &
            & dom_assoc_carbon, &
            & residue_assoc_carbon, &
            & qmax_org(ic, is, iblk), &
            & bulk_carbon_bio)
          bulk_carbon_unconstrained = bulk_carbon_bio

          surface_burial_frac(ic, is, iblk) = 1.0_wp
          IF (l_surface_layer .AND. bulk_carbon_bio > bulk_carbon_max(ic, is, iblk)) THEN
            n_surface_constraints = n_surface_constraints + 1
            !> For the capped surface layer, qmax_org is reset by the prescribed minimum
            !> mineral-soil volume rather than by the unconstrained analytical C_B solution:
            !>   Qmax_org = phi_min^min * qmax_org_min * rho_bulk
            qmax_org_target = min_volume_min_layer * qmax_org_min_local * bulk_dens_local
            IF (qmax_org(ic, is, iblk) > eps8) THEN
              assoc_scale = qmax_org_target / qmax_org(ic, is, iblk)
            ELSE
              assoc_scale = 0.0_wp
            END IF
            qmax_org(ic, is, iblk) = qmax_org_target
            dom_assoc_carbon       = dom_assoc_carbon * assoc_scale
            residue_assoc_carbon   = residue_assoc_carbon * assoc_scale
            !> C_bulk,1^soil = f_burial * (C_sol* + C_poly* + C_B* + C_D* + C_res*) + C_aDOM* + C_aRes*
            !> f_burial = (C_bulk,max - (C_aDOM* + C_aRes*))
            !>              / (C_sol* + C_poly* + C_B* + C_D* + C_res*)
            surface_bulk_nonsorbed = input_sol * tau_soluable_litter + polymeric_carbon + microbial_c_bio &
              &                     + sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk) + residue_carbon
            surface_bulk_assoc     = dom_assoc_carbon + residue_assoc_carbon
            IF (surface_bulk_nonsorbed > eps8) THEN
              surface_burial_frac(ic, is, iblk) = (bulk_carbon_max(ic, is, iblk) - surface_bulk_assoc) / surface_bulk_nonsorbed
              surface_burial_frac(ic, is, iblk) = MIN(1.0_wp, MAX(0.0_wp, surface_burial_frac(ic, is, iblk)))
            ELSE
              surface_burial_frac(ic, is, iblk) = 0.0_wp
            END IF
            IF (l_print_constraint_debug) THEN
              WRITE(message_text,*) 'DEBUG analytical physical constraint surface_burial ic/is/iblk = ', &
                & ic, is, iblk, '; C_bulk_unconstrained, C_bulk_max, burial_frac = ', &
                & bulk_carbon_unconstrained, bulk_carbon_max(ic, is, iblk), surface_burial_frac(ic, is, iblk), &
                & '; input_sol/poly/wood, C_microbial/poly/residue = ', input_sol, input_poly, input_wood, &
                & microbial_c_bio, polymeric_carbon, residue_carbon
              CALL message(TRIM(routine), TRIM(message_text))
            END IF
          ELSE IF (bulk_carbon_bio > bulk_carbon_max(ic, is, iblk)) THEN
            n_enzyme_constraints = n_enzyme_constraints + 1
            constraint_iterations = 0
            enzyme_scale_lower = 1.0_wp
            enzyme_scale_upper = MAX(enzyme_frac_total_phys_scale, bulk_carbon_bio / bulk_carbon_max(ic, is, iblk))
            enzyme_frac_total_scale = enzyme_scale_upper
            CALL calc_analytical_ss_jsm_layer_state( &
              & l_surface_layer, &
              & input_sol, &
              & input_poly, &
              & input_wood, &
              & sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk), &
              & microbial_c_bio, &
              & l_microbial_dormant, &
              & enzyme_frac_AP(ic, is, iblk), &
              & enzyme_frac_total_scale, &
              & transport_poly_local, transport_residue_local, &
              & transport_dom_assoc_local, transport_residue_assoc_local, &
              & km_depoly_eff, &
              & vmax_poly_eff, &
              & vmax_res_eff, &
              & qmax_org_min_local, &
              & bulk_dens_local, &
              & rtm_sorb_eff, &
              & rmm_sorb_eff, &
              & rtm_desorb_eff, &
              & rmm_desorb_eff, &
              & alpha_poly, &
              & alpha_residue, &
              & polymeric_carbon, &
              & residue_carbon, &
              & dom_assoc_carbon, &
              & residue_assoc_carbon, &
              & qmax_org(ic, is, iblk), &
              & bulk_carbon_scale_upper)
            IF (bulk_carbon_scale_upper <= bulk_carbon_max(ic, is, iblk)) THEN
              cb_bracket_tol = eps4 * MAX(1.0_wp, enzyme_scale_upper)
              !> Bisection solve for enzyme_frac_total_scale such that
              !> C_bulk(C_B^bio, enzyme_frac_total_scale) = C_bulk^max
              DO iter_phys_constraint = 1, max_iter_phys_constraint
                constraint_iterations = iter_phys_constraint
                enzyme_scale_mid = 0.5_wp * (enzyme_scale_lower + enzyme_scale_upper)
                CALL calc_analytical_ss_jsm_layer_state( &
                  & l_surface_layer, &
                  & input_sol, &
                  & input_poly, &
                  & input_wood, &
                  & sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk), &
                  & microbial_c_bio, &
                  & l_microbial_dormant, &
                  & enzyme_frac_AP(ic, is, iblk), &
                  & enzyme_scale_mid, &
                  & transport_poly_local, transport_residue_local, &
                  & transport_dom_assoc_local, transport_residue_assoc_local, &
                  & km_depoly_eff, &
                  & vmax_poly_eff, &
                  & vmax_res_eff, &
                  & qmax_org_min_local, &
                  & bulk_dens_local, &
                  & rtm_sorb_eff, &
                  & rmm_sorb_eff, &
                  & rtm_desorb_eff, &
                  & rmm_desorb_eff, &
                  & alpha_poly, &
                  & alpha_residue, &
                  & polymeric_carbon, &
                  & residue_carbon, &
                  & dom_assoc_carbon, &
                  & residue_assoc_carbon, &
                  & qmax_org(ic, is, iblk), &
                  & bulk_carbon_mid)
                IF (bulk_carbon_mid > bulk_carbon_max(ic, is, iblk)) THEN
                  enzyme_scale_lower = enzyme_scale_mid
                ELSE
                  enzyme_scale_upper = enzyme_scale_mid
                  bulk_carbon_scale_upper = bulk_carbon_mid
                  IF ((bulk_carbon_max(ic, is, iblk) - bulk_carbon_mid) <= 0.01_wp * bulk_carbon_max(ic, is, iblk)) EXIT
                END IF
                IF ((enzyme_scale_upper - enzyme_scale_lower) <= cb_bracket_tol) EXIT
              END DO
              enzyme_frac_total_scale = enzyme_scale_upper
            ELSE IF (bulk_carbon_bio > bulk_carbon_max(ic, is, iblk)) THEN
              enzyme_frac_total_scale = enzyme_scale_upper
            END IF

            CALL calc_analytical_ss_jsm_layer_state( &
              & l_surface_layer, &
              & input_sol, &
              & input_poly, &
              & input_wood, &
              & sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk), &
              & microbial_c_bio, &
              & l_microbial_dormant, &
              & enzyme_frac_AP(ic, is, iblk), &
              & enzyme_frac_total_scale, &
              & transport_poly_local, transport_residue_local, &
              & transport_dom_assoc_local, transport_residue_assoc_local, &
              & km_depoly_eff, &
              & vmax_poly_eff, &
              & vmax_res_eff, &
              & qmax_org_min_local, &
              & bulk_dens_local, &
              & rtm_sorb_eff, &
              & rmm_sorb_eff, &
              & rtm_desorb_eff, &
              & rmm_desorb_eff, &
              & alpha_poly, &
              & alpha_residue, &
              & polymeric_carbon, &
              & residue_carbon, &
              & dom_assoc_carbon, &
              & residue_assoc_carbon, &
              & qmax_org(ic, is, iblk), &
              & bulk_carbon_bio)
            max_enzyme_frac_total_scale = MAX(max_enzyme_frac_total_scale, enzyme_frac_total_scale)
            IF (l_print_constraint_debug) THEN
              WRITE(message_text,*) 'DEBUG analytical physical constraint enzyme_scale ic/is/iblk = ', &
                & ic, is, iblk, '; C_bulk_unconstrained, C_bulk_max, C_bulk_corrected, scale, iterations = ', &
                & bulk_carbon_unconstrained, bulk_carbon_max(ic, is, iblk), bulk_carbon_bio, &
                & enzyme_frac_total_scale, constraint_iterations, &
                & '; enzyme_frac_base/effective, input_sol/poly/wood, C_microbial/poly/residue = ', &
                & enzyme_frac_total, enzyme_frac_total * enzyme_frac_total_scale, &
                & input_sol, input_poly, input_wood, microbial_c_bio, polymeric_carbon, residue_carbon
              CALL message(TRIM(routine), TRIM(message_text))
            END IF
          END IF

          sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk) = microbial_c_bio
          enzyme_frac_poly(ic, is, iblk)    = alpha_poly
          enzyme_frac_residue(ic, is, iblk) = alpha_residue
          !> In the C-only draft, the carbon-based analytical allocation is the only closed solution.
          !> Keep the N and P diagnostic fractions aligned with the carbon solution until the
          !> nutrient-limited analytical system is added.
          enzyme_frac_poly_c(ic, is, iblk)  = alpha_poly
          enzyme_frac_poly_n(ic, is, iblk)  = alpha_poly
          enzyme_frac_poly_p(ic, is, iblk)  = alpha_poly
          sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) = polymeric_carbon
          sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk)          = residue_carbon
          sb_pool_mt_domain(ix_dom_assoc, ixC, ic, is, iblk)        = dom_assoc_carbon
          sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk)    = residue_assoc_carbon
        END DO
      END DO
    END DO

    IF (l_print_constraint_debug) THEN
      WRITE(message_text,*) 'DEBUG analytical physical constraint summary', &
        & '; surface_burial_count, deeper_enzyme_scale_count, max_enzyme_scale = ', &
        & n_surface_constraints, n_enzyme_constraints, max_enzyme_frac_total_scale
      CALL message(TRIM(routine), TRIM(message_text))
    END IF

  END SUBROUTINE calc_analytical_ss_jsm

  ! ======================================================================================================= !
  !> Calculate and apply one analytical steady-state target on a model chunk.
  !> Carbon follows the analytical solution. C14 follows its own steady-state tracer solution including
  !> radioactive decay; other elements retain their current pool-specific ratios.
  !>
  !> Current accelerator workflow:
  !>
  !>   time-mean forcing and pre-jump state
  !>                  |
  !>                  v
  !>   solve transport-free analytical C target
  !>                  |
  !>                  v
  !>   enforce physical constraints
  !>   (minimum mineral volume, surface burial, qmax, enzyme scaling)
  !>                  |
  !>                  +---------------- bioturbation smoothing disabled ------------+
  !>                  |                                                             |
  !>                  v                                                             |
  !>   integrate bioturbation directly over transport_smoothing_years                |
  !>   while updating bulk density and mineral volume annually                       |
  !>                  |                                                             |
  !>                  v                                                             |
  !>   use the mass-conserving bioturbation-smoothed target                          |
  !>                  |                                                             |
  !>                  +-------------------------------------------------------------+
  !>                  |
  !>                  v
  !>   reconstruct N, P, and isotope pools from target C
  !>   (C14 uses a decay-aware steady-state solution)
  !>                  |
  !>                  v
  !>   replace model pools and associated qmax, burial, and enzyme diagnostics
  !>
  !> The jump is instantaneous in model time. Particle transport is deliberately
  !> excluded here and subsequently evolves from actual forward-model formation and loss rates.
  SUBROUTINE apply_analytical_ss_jsm_accelerator( &
    & nc, nsoil_w, num_sl_above_bedrock, soil_depth_sl, &
    & soil_lay_depth_center_sl, soil_lay_depth_ubound_sl, &
    & input_sol_c, input_poly_c, input_wood_c, input_sol_c14, input_poly_c14, input_wood_c14, &
    & l_c14_enabled, doc_transport, enzyme_ap_forcing, &
    & bulk_dens_sl, bulk_carbon_max, qmax_org_min_sl, depoly_mod, uptake_mod, &
    & sorption_mod, desorption_mod, hsc_mod, l_use_bioturbation, l_transport_debug, &
    & transport_smoothing_years, bioturb_numerator, &
    & sb_pool_mt, qmax_org, surface_burial_frac, &
    & enzyme_frac_poly, enzyme_frac_residue, enzyme_frac_poly_c, enzyme_frac_poly_n, &
    & enzyme_frac_poly_p, enzyme_frac_AP)

    USE mo_sb_constants, ONLY: microbial_cn, microbial_np, k_slow_som_np, &
      & k_init_soluable_litter_cn, k_init_polymeric_litter_cn, &
      & k_init_woody_litter_cn, k_init_litter_np, fc_soluable2dom, &
      & fc_woody2polymeric, f_soluable_residue, k_mort_microbial, &
      & k_desorpt_dom, k_desorpt_det
    USE mo_jsb_physical_constants, ONLY: lambda_C14
    USE mo_q_sb_jsm_processes, ONLY: calc_bulk_soil_carbon, calc_bulk_density_correction

    INTEGER,  INTENT(in)    :: nc                  !< number of horizontal points in this chunk [-]
    INTEGER,  INTENT(in)    :: nsoil_w             !< number of soil layers [-]
    REAL(wp), INTENT(in)    :: num_sl_above_bedrock(:) !< active soil-layer count per point [-]
    REAL(wp), INTENT(in)    :: soil_depth_sl(:,:)  !< soil-layer thickness [m]
    REAL(wp), INTENT(in)    :: soil_lay_depth_center_sl(:,:) !< soil-layer midpoint depth [m]
    REAL(wp), INTENT(in)    :: soil_lay_depth_ubound_sl(:,:) !< soil-layer lower-bound depth [m]
    REAL(wp), INTENT(in)    :: input_sol_c(:,:)    !< mean soluble-litter C input [mol C m-3 s-1]
    REAL(wp), INTENT(in)    :: input_poly_c(:,:)   !< mean polymeric-litter C input [mol C m-3 s-1]
    REAL(wp), INTENT(in)    :: input_wood_c(:,:)   !< mean woody-litter C input [mol C m-3 s-1]
    REAL(wp), INTENT(in)    :: input_sol_c14(:,:)  !< mean soluble-litter C14 input [mol C14 m-3 s-1]
    REAL(wp), INTENT(in)    :: input_poly_c14(:,:) !< mean polymeric-litter C14 input [mol C14 m-3 s-1]
    REAL(wp), INTENT(in)    :: input_wood_c14(:,:) !< mean woody-litter C14 input [mol C14 m-3 s-1]
    LOGICAL,  INTENT(in)    :: l_c14_enabled       !< TRUE when C14 pools are active
    REAL(wp), INTENT(in)    :: doc_transport(:,:)  !< downward DOC transport coefficient [s-1]
    REAL(wp), INTENT(in)    :: enzyme_ap_forcing(:,:) !< mean phosphatase enzyme fraction [-]
    REAL(wp), INTENT(in)    :: bulk_dens_sl(:,:)   !< mineral soil bulk density [kg m-3]
    REAL(wp), INTENT(inout) :: bulk_carbon_max(:,:)!< maximum admissible soil C density [mol C m-3]
    REAL(wp), INTENT(in)    :: qmax_org_min_sl(:,:)!< texture-dependent sorption capacity [mol C kg-1]
    REAL(wp), INTENT(in)    :: depoly_mod(:,:)     !< combined environmental depolymerisation modifier [-]
    REAL(wp), INTENT(in)    :: uptake_mod(:,:)     !< combined environmental microbial-uptake modifier [-]
    REAL(wp), INTENT(in)    :: sorption_mod(:,:)   !< combined environmental sorption modifier [-]
    REAL(wp), INTENT(in)    :: desorption_mod(:,:) !< combined environmental desorption modifier [-]
    REAL(wp), INTENT(in)    :: hsc_mod(:,:)        !< combined half-saturation-constant modifier [-]
    LOGICAL,  INTENT(in)    :: l_use_bioturbation  !< smooth the analytical target with bioturbation
    LOGICAL,  INTENT(in)    :: l_transport_debug   !< print jump summaries and profiles
    INTEGER,  INTENT(in)    :: transport_smoothing_years !< direct bioturbation integration period [yr]
    REAL(wp), INTENT(in)    :: bioturb_numerator(:,:) !< numerator of bioturbation diffusivity [m2 s-1 kg m-3]
    REAL(wp), INTENT(inout) :: sb_pool_mt(:,:,:,:) !< pool by element, point and layer [mol m-3]
    REAL(wp), INTENT(inout) :: qmax_org(:,:)       !< organic-matter sorption capacity [mol C m-3]
    REAL(wp), INTENT(inout) :: surface_burial_frac(:,:) !< buried fraction of first-layer surface SOM [-]
    REAL(wp), INTENT(inout) :: enzyme_frac_poly(:,:) !< polymeric enzyme allocation [-]
    REAL(wp), INTENT(inout) :: enzyme_frac_residue(:,:) !< residue enzyme allocation [-]
    REAL(wp), INTENT(inout) :: enzyme_frac_poly_c(:,:) !< C-based polymeric enzyme allocation [-]
    REAL(wp), INTENT(inout) :: enzyme_frac_poly_n(:,:) !< N-based polymeric enzyme allocation [-]
    REAL(wp), INTENT(inout) :: enzyme_frac_poly_p(:,:) !< P-based polymeric enzyme allocation [-]
    REAL(wp), INTENT(inout) :: enzyme_frac_AP(:,:) !< phosphatase fraction of total enzyme [-]

    REAL(wp), ALLOCATABLE :: pool_before(:,:,:,:)      !< pre-jump element pools [mol m-3]
    REAL(wp), ALLOCATABLE :: pool_domain(:,:,:,:,:)    !< accepted analytical element pools [mol m-3]
    REAL(wp), ALLOCATABLE :: candidate_domain(:,:,:,:,:) !< bioturbation-smoothed target pools [mol m-3]
    REAL(wp), ALLOCATABLE :: num_sl_domain(:,:)        !< active soil-layer count with block dimension [-]
    REAL(wp), ALLOCATABLE :: zero_domain(:,:,:)        !< reusable zero modifier/forcing array [-]
    REAL(wp), ALLOCATABLE :: one_domain(:,:,:)         !< reusable unit modifier array [-]
    REAL(wp), ALLOCATABLE :: bulk_carbon_max_domain(:,:,:) !< accepted maximum bulk C [mol C m-3]
    REAL(wp), ALLOCATABLE :: qmax_org_domain(:,:,:)    !< accepted sorption capacity [mol C m-3]
    REAL(wp), ALLOCATABLE :: surface_burial_domain(:,:,:) !< accepted surface burial fraction [-]
    REAL(wp), ALLOCATABLE :: enzyme_poly_domain(:,:,:)   !< accepted polymeric enzyme fraction [-]
    REAL(wp), ALLOCATABLE :: enzyme_res_domain(:,:,:)    !< accepted residue enzyme fraction [-]
    REAL(wp), ALLOCATABLE :: enzyme_poly_c_domain(:,:,:) !< accepted C-based enzyme fraction [-]
    REAL(wp), ALLOCATABLE :: enzyme_poly_n_domain(:,:,:) !< accepted N-based enzyme fraction [-]
    REAL(wp), ALLOCATABLE :: enzyme_poly_p_domain(:,:,:) !< accepted P-based enzyme fraction [-]
    REAL(wp), ALLOCATABLE :: enzyme_ap_domain(:,:,:)     !< accepted phosphatase fraction [-]
    REAL(wp), ALLOCATABLE :: qmax_org_min_domain(:,:,:)  !< prescribed texture capacity [mol C kg-1]
    REAL(wp), ALLOCATABLE :: candidate_bulk_soil_c(:,:) !< candidate bulk soil C density [mol C m-3]
    REAL(wp), ALLOCATABLE :: candidate_litter_c(:,:)    !< candidate soil litter C density [mol C m-3]
    REAL(wp), ALLOCATABLE :: candidate_volume_min(:,:)  !< candidate mineral soil volume fraction [-]
    REAL(wp), ALLOCATABLE :: candidate_bulk_density(:,:)!< candidate corrected bulk density [kg m-3]
    REAL(wp), ALLOCATABLE :: candidate_k_bioturb(:,:)   !< candidate bioturbation diffusivity [m2 s-1]
    REAL(wp), ALLOCATABLE :: candidate_surface_burial(:,:) !< mutable candidate burial fraction [-]
    REAL(wp), ALLOCATABLE :: bioturbation_pool(:,:) !< buried pool undergoing bioturbation [mol C m-3]
    REAL(wp), ALLOCATABLE :: transported_pool(:,:) !< pool after one implicit bioturbation year [mol C m-3]
    REAL(wp), ALLOCATABLE :: debug_bulk_c_before(:,:)   !< pre-jump bulk C density [mol C m-3]
    REAL(wp), ALLOCATABLE :: debug_bulk_c_after(:,:)    !< post-jump bulk C density [mol C m-3]
    REAL(wp), ALLOCATABLE :: debug_litter_c_before(:,:) !< pre-jump litter C density [mol C m-3]
    REAL(wp), ALLOCATABLE :: debug_litter_c_after(:,:)  !< post-jump litter C density [mol C m-3]
    REAL(wp), ALLOCATABLE :: debug_density_before(:,:)  !< pre-jump corrected bulk density [kg m-3]
    REAL(wp), ALLOCATABLE :: debug_density_after(:,:)   !< post-jump corrected bulk density [kg m-3]
    REAL(wp), ALLOCATABLE :: debug_volume_min_before(:,:) !< pre-jump mineral volume fraction [-]
    REAL(wp), ALLOCATABLE :: debug_volume_min_after(:,:) !< post-jump mineral volume fraction [-]
    REAL(wp), ALLOCATABLE :: debug_burial_before(:,:)   !< pre-jump surface burial fraction [-]
    REAL(wp), ALLOCATABLE :: debug_burial_after(:,:)    !< post-jump surface burial fraction [-]
    INTEGER :: ic                                      !< horizontal point index [-]
    INTEGER :: is                                      !< soil-layer index [-]
    CHARACTER(len=1024) :: message_text                 !< formatted model-log message
    INTEGER, PARAMETER :: bioturb_pools(5) = [ix_polymeric_litter, ix_microbial, ix_residue, & !< pool IDs mixed by bioturbation
      & ix_dom_assoc, ix_residue_assoc]

    ALLOCATE(pool_before(SIZE(sb_pool_mt,1), SIZE(sb_pool_mt,2), nc, nsoil_w))
    ALLOCATE(pool_domain(SIZE(sb_pool_mt,1), SIZE(sb_pool_mt,2), nc, nsoil_w, 1))
    ALLOCATE(num_sl_domain(nc,1), zero_domain(nc,nsoil_w,1), one_domain(nc,nsoil_w,1))
    ALLOCATE(bulk_carbon_max_domain(nc,nsoil_w,1), qmax_org_domain(nc,nsoil_w,1))
    ALLOCATE(surface_burial_domain(nc,nsoil_w,1), enzyme_poly_domain(nc,nsoil_w,1))
    ALLOCATE(enzyme_res_domain(nc,nsoil_w,1), enzyme_poly_c_domain(nc,nsoil_w,1))
    ALLOCATE(enzyme_poly_n_domain(nc,nsoil_w,1), enzyme_poly_p_domain(nc,nsoil_w,1))
    ALLOCATE(enzyme_ap_domain(nc,nsoil_w,1), qmax_org_min_domain(nc,nsoil_w,1))
    IF (l_use_bioturbation .OR. l_transport_debug) THEN
      ALLOCATE(candidate_domain(SIZE(sb_pool_mt,1), SIZE(sb_pool_mt,2), nc, nsoil_w, 1))
      ALLOCATE(candidate_bulk_soil_c(nc,nsoil_w), candidate_litter_c(nc,nsoil_w))
      ALLOCATE(candidate_volume_min(nc,nsoil_w), candidate_bulk_density(nc,nsoil_w))
      ALLOCATE(candidate_k_bioturb(nc,nsoil_w), candidate_surface_burial(nc,nsoil_w))
      ALLOCATE(bioturbation_pool(nc,nsoil_w), transported_pool(nc,nsoil_w))
      IF (l_transport_debug) THEN
        ALLOCATE(debug_bulk_c_before(nc,nsoil_w), debug_bulk_c_after(nc,nsoil_w))
        ALLOCATE(debug_litter_c_before(nc,nsoil_w), debug_litter_c_after(nc,nsoil_w))
        ALLOCATE(debug_density_before(nc,nsoil_w), debug_density_after(nc,nsoil_w))
        ALLOCATE(debug_volume_min_before(nc,nsoil_w), debug_volume_min_after(nc,nsoil_w))
        ALLOCATE(debug_burial_before(nc,nsoil_w), debug_burial_after(nc,nsoil_w))
      END IF
    END IF

    pool_before = sb_pool_mt
    pool_domain(:,:,:,:,1) = sb_pool_mt
    num_sl_domain(:,1) = num_sl_above_bedrock
    zero_domain = 0.0_wp
    one_domain = 1.0_wp
    bulk_carbon_max_domain(:,:,1) = bulk_carbon_max
    qmax_org_domain(:,:,1) = qmax_org
    surface_burial_domain(:,:,1) = surface_burial_frac
    qmax_org_min_domain(:,:,1) = qmax_org_min_sl

    ! Stage 1: solve and physically constrain the local analytical target.
    CALL solve_particle_free_constrained_target()

    ! Stage 2: optionally smooth the target directly with mass-conserving bioturbation.
    IF (l_use_bioturbation) CALL apply_bioturbation_target_smoothing()

    IF (l_transport_debug) CALL write_analytical_target_debug()

    CALL reconstruct_target_elements(nc, nsoil_w, num_sl_above_bedrock, pool_before, pool_domain(:,:,:,:,1))
    IF (l_c14_enabled) THEN
      IF (ixC14 > SIZE(sb_pool_mt,2)) THEN
        CALL finish(TRIM(modname)//':apply_analytical_ss_jsm_accelerator', &
          & 'C14 is enabled but absent from the SB element matrix.')
      END IF
      CALL reconstruct_target_c14(nc, nsoil_w, num_sl_above_bedrock, &
        & input_sol_c, input_poly_c, input_wood_c, input_sol_c14, input_poly_c14, input_wood_c14, &
        & doc_transport, desorption_mod, pool_domain(:,:,:,:,1))
    END IF

    sb_pool_mt = pool_domain(:,:,:,:,1)
    bulk_carbon_max = bulk_carbon_max_domain(:,:,1)
    qmax_org = qmax_org_domain(:,:,1)
    surface_burial_frac = surface_burial_domain(:,:,1)
    enzyme_frac_poly = enzyme_poly_domain(:,:,1)
    enzyme_frac_residue = enzyme_res_domain(:,:,1)
    enzyme_frac_poly_c = enzyme_poly_c_domain(:,:,1)
    enzyme_frac_poly_n = enzyme_poly_n_domain(:,:,1)
    enzyme_frac_poly_p = enzyme_poly_p_domain(:,:,1)
    enzyme_frac_AP = enzyme_ap_domain(:,:,1)

  CONTAINS

    SUBROUTINE solve_particle_free_constrained_target()
      ! No particle-tendency arguments are supplied here. Consequently the
      ! analytical solve uses only biological forcing and then applies burial,
      ! mineral-volume, qmax, and enzyme-fraction constraints.
      CALL solve_analytical_target(pool_domain, bulk_carbon_max_domain, qmax_org_domain, &
        & surface_burial_domain, enzyme_poly_domain, enzyme_res_domain, enzyme_poly_c_domain, &
        & enzyme_poly_n_domain, enzyme_poly_p_domain, enzyme_ap_domain)
    END SUBROUTINE solve_particle_free_constrained_target

    SUBROUTINE write_analytical_target_debug()
      candidate_domain(:,:,:,:,1) = pool_before
      debug_burial_before = surface_burial_frac
      CALL update_candidate_soil_state(candidate_domain, debug_burial_before)
      debug_bulk_c_before = candidate_bulk_soil_c
      debug_litter_c_before = candidate_litter_c
      debug_density_before = candidate_bulk_density
      debug_volume_min_before = candidate_volume_min

      candidate_domain = pool_domain
      debug_burial_after = surface_burial_domain(:,:,1)
      CALL update_candidate_soil_state(candidate_domain, debug_burial_after)
      debug_bulk_c_after = candidate_bulk_soil_c
      debug_litter_c_after = candidate_litter_c
      debug_density_after = candidate_bulk_density
      debug_volume_min_after = candidate_volume_min

      IF (l_use_bioturbation) THEN
        WRITE(message_text,*) 'DEBUG analytical jump with direct bioturbation-smoothed target'
      ELSE
        WRITE(message_text,*) 'DEBUG analytical direct jump without transport'
      END IF
      CALL message(TRIM(modname)//':apply_analytical_ss_jsm_accelerator', TRIM(message_text))
      CALL write_debug_carbon_totals('analytical target jump', &
        & pool_before, candidate_domain(:,:,:,:,1), &
        & debug_bulk_c_before, debug_bulk_c_after, debug_litter_c_before, debug_litter_c_after)

      DO ic = 1, nc
        DO is = 1, NINT(num_sl_above_bedrock(ic))
          WRITE(message_text,*) 'DEBUG analytical target state ic/is = ', ic, is, &
            & '; C_bulk before/after/delta = ', debug_bulk_c_before(ic,is), &
            & debug_bulk_c_after(ic,is), debug_bulk_c_after(ic,is) - debug_bulk_c_before(ic,is), &
            & '; C_litter before/after/delta = ', debug_litter_c_before(ic,is), &
            & debug_litter_c_after(ic,is), debug_litter_c_after(ic,is) - debug_litter_c_before(ic,is), &
            & '; bulk_density before/after/delta = ', debug_density_before(ic,is), &
            & debug_density_after(ic,is), debug_density_after(ic,is) - debug_density_before(ic,is), &
            & '; mineral_volume before/after/delta = ', debug_volume_min_before(ic,is), &
            & debug_volume_min_after(ic,is), debug_volume_min_after(ic,is) - debug_volume_min_before(ic,is), &
            & '; burial before/after/delta = ', debug_burial_before(ic,is), &
            & debug_burial_after(ic,is), debug_burial_after(ic,is) - debug_burial_before(ic,is)
          CALL message(TRIM(modname)//':apply_analytical_ss_jsm_accelerator', TRIM(message_text))
        END DO
      END DO
    END SUBROUTINE write_analytical_target_debug

    SUBROUTINE write_debug_carbon_totals(label, state_before, state_after, &
      & bulk_before, bulk_after, litter_before, litter_after)

      CHARACTER(len=*), INTENT(in) :: label !< diagnostic label identifying the compared states
      REAL(wp), INTENT(in) :: state_before(:,:,:,:) !< pre-jump element pools [mol m-3]
      REAL(wp), INTENT(in) :: state_after(:,:,:,:)  !< post-jump element pools [mol m-3]
      REAL(wp), INTENT(in) :: bulk_before(:,:)      !< pre-jump bulk soil C density [mol C m-3]
      REAL(wp), INTENT(in) :: bulk_after(:,:)       !< post-jump bulk soil C density [mol C m-3]
      REAL(wp), INTENT(in) :: litter_before(:,:)    !< pre-jump litter C density [mol C m-3]
      REAL(wp), INTENT(in) :: litter_after(:,:)     !< post-jump litter C density [mol C m-3]

      REAL(wp) :: total_bulk_before    !< depth-integrated pre-jump bulk C [mol C m-2]
      REAL(wp) :: total_bulk_after     !< depth-integrated post-jump bulk C [mol C m-2]
      REAL(wp) :: total_litter_before  !< depth-integrated pre-jump litter C [mol C m-2]
      REAL(wp) :: total_litter_after   !< depth-integrated post-jump litter C [mol C m-2]
      REAL(wp) :: total_soc_before     !< depth-integrated pre-jump total soil C [mol C m-2]
      REAL(wp) :: total_soc_after      !< depth-integrated post-jump total soil C [mol C m-2]
      INTEGER  :: ic_local             !< horizontal point index [-]
      INTEGER  :: is_local             !< soil-layer index [-]
      CHARACTER(len=1024) :: totals_message !< formatted model-log message

      total_bulk_before = 0.0_wp
      total_bulk_after = 0.0_wp
      total_litter_before = 0.0_wp
      total_litter_after = 0.0_wp
      total_soc_before = 0.0_wp
      total_soc_after = 0.0_wp
      DO ic_local = 1, nc
        DO is_local = 1, NINT(num_sl_above_bedrock(ic_local))
          total_bulk_before = total_bulk_before + &
            & bulk_before(ic_local,is_local) * soil_depth_sl(ic_local,is_local)
          total_bulk_after = total_bulk_after + &
            & bulk_after(ic_local,is_local) * soil_depth_sl(ic_local,is_local)
          total_litter_before = total_litter_before + &
            & litter_before(ic_local,is_local) * soil_depth_sl(ic_local,is_local)
          total_litter_after = total_litter_after + &
            & litter_after(ic_local,is_local) * soil_depth_sl(ic_local,is_local)
          total_soc_before = total_soc_before + &
            & SUM(state_before(:,ixC,ic_local,is_local)) * soil_depth_sl(ic_local,is_local)
          total_soc_after = total_soc_after + &
            & SUM(state_after(:,ixC,ic_local,is_local)) * soil_depth_sl(ic_local,is_local)
        END DO
      END DO

      WRITE(totals_message,*) 'DEBUG ', TRIM(label), ' carbon totals [mol C m-2]', &
        & '; C_bulk before/after/delta = ', total_bulk_before, total_bulk_after, &
        & total_bulk_after - total_bulk_before, &
        & '; C_litter before/after/delta = ', total_litter_before, total_litter_after, &
        & total_litter_after - total_litter_before, &
        & '; SOC before/after/delta = ', total_soc_before, total_soc_after, &
        & total_soc_after - total_soc_before
      CALL message(TRIM(modname)//':apply_analytical_ss_jsm_accelerator', TRIM(totals_message))
    END SUBROUTINE write_debug_carbon_totals

    SUBROUTINE solve_analytical_target(target_pool, target_bulk_carbon_max, target_qmax_org, &
      & target_surface_burial, target_enzyme_poly, target_enzyme_res, target_enzyme_poly_c, &
      & target_enzyme_poly_n, target_enzyme_poly_p, target_enzyme_ap)

      REAL(wp), INTENT(out) :: target_pool(:,:,:,:,:) !< solved pool by element, point, layer, block [mol m-3]
      REAL(wp), INTENT(out) :: target_bulk_carbon_max(:,:,:) !< solved maximum bulk C [mol C m-3]
      REAL(wp), INTENT(out) :: target_qmax_org(:,:,:) !< solved sorption capacity [mol C m-3]
      REAL(wp), INTENT(out) :: target_surface_burial(:,:,:) !< solved surface burial fraction [-]
      REAL(wp), INTENT(out) :: target_enzyme_poly(:,:,:) !< solved polymeric enzyme allocation [-]
      REAL(wp), INTENT(out) :: target_enzyme_res(:,:,:) !< solved residue enzyme allocation [-]
      REAL(wp), INTENT(out) :: target_enzyme_poly_c(:,:,:) !< solved C-based allocation [-]
      REAL(wp), INTENT(out) :: target_enzyme_poly_n(:,:,:) !< solved N-based allocation [-]
      REAL(wp), INTENT(out) :: target_enzyme_poly_p(:,:,:) !< solved P-based allocation [-]
      REAL(wp), INTENT(out) :: target_enzyme_ap(:,:,:) !< solved phosphatase fraction [-]

      target_pool(:,:,:,:,1) = pool_before
      target_bulk_carbon_max(:,:,1) = bulk_carbon_max
      target_qmax_org(:,:,1) = qmax_org
      target_surface_burial(:,:,1) = surface_burial_frac
      target_enzyme_poly(:,:,1) = enzyme_frac_poly
      target_enzyme_res(:,:,1) = enzyme_frac_residue
      target_enzyme_poly_c(:,:,1) = enzyme_frac_poly_c
      target_enzyme_poly_n(:,:,1) = enzyme_frac_poly_n
      target_enzyme_poly_p(:,:,1) = enzyme_frac_poly_p
      target_enzyme_ap(:,:,1) = enzyme_frac_AP

      CALL calc_analytical_ss_jsm( &
        & nc, nsoil_w, 1, num_sl_domain, zero_domain, &
        & RESHAPE(soil_lay_depth_center_sl, [nc,nsoil_w,1]), &
        & RESHAPE(soil_lay_depth_ubound_sl, [nc,nsoil_w,1]), &
        & RESHAPE(input_sol_c, [nc,nsoil_w,1]), RESHAPE(input_poly_c, [nc,nsoil_w,1]), &
        & RESHAPE(input_wood_c, [nc,nsoil_w,1]), RESHAPE(doc_transport, [nc,nsoil_w,1]), &
        & RESHAPE(enzyme_ap_forcing, [nc,nsoil_w,1]), zero_domain, zero_domain, &
        & RESHAPE(bulk_dens_sl, [nc,nsoil_w,1]), target_bulk_carbon_max, target_qmax_org, &
        & RESHAPE(depoly_mod, [nc,nsoil_w,1]), one_domain, RESHAPE(uptake_mod, [nc,nsoil_w,1]), one_domain, &
        & RESHAPE(sorption_mod, [nc,nsoil_w,1]), one_domain, &
        & RESHAPE(desorption_mod, [nc,nsoil_w,1]), one_domain, RESHAPE(hsc_mod, [nc,nsoil_w,1]), one_domain, &
        & target_pool, target_surface_burial, target_enzyme_poly, target_enzyme_res, &
        & target_enzyme_poly_c, target_enzyme_poly_n, target_enzyme_poly_p, target_enzyme_ap, &
        & l_constraint_debug=l_transport_debug, qmax_org_min_forcing=qmax_org_min_domain)
    END SUBROUTINE solve_analytical_target

    SUBROUTINE apply_bioturbation_target_smoothing()
      REAL(wp) :: unburied_surface(nc) !< surface material excluded from mineral-soil mixing [mol C m-3]
      REAL(wp) :: max_pool_change      !< maximum absolute C-density change caused by smoothing [mol C m-3]
      REAL(wp) :: max_bulk_excess      !< maximum post-smoothing excess above physical bulk-C limit [mol C m-3]
      INTEGER  :: smoothing_year       !< annual bioturbation integration step [-]
      INTEGER  :: transport_pool_index !< index in bioturb_pools [-]
      INTEGER  :: local_ic             !< horizontal point index [-]
      INTEGER  :: local_is             !< soil-layer index [-]

      IF (transport_smoothing_years <= 0) THEN
        CALL finish(TRIM(modname)//':apply_bioturbation_target_smoothing', &
          & 'The analytical bioturbation smoothing period must be greater than zero.')
      END IF

      candidate_domain = pool_domain
      candidate_surface_burial = surface_burial_domain(:,:,1)
      IF (l_transport_debug) THEN
        CALL write_bioturbation_target_profile('transport-free target', candidate_domain, &
          & bulk_carbon_max_domain(:,:,1), qmax_org_domain(:,:,1), candidate_surface_burial, &
          & enzyme_poly_domain(:,:,1), enzyme_res_domain(:,:,1), enzyme_ap_domain(:,:,1))
      END IF

      DO smoothing_year = 1, transport_smoothing_years
        ! Recalculate organic volume and bulk density before each annual mixing
        ! step because the bioturbation coefficient depends on bulk density.
        CALL update_candidate_soil_state(candidate_domain, candidate_surface_burial)
        candidate_k_bioturb = 0.0_wp
        WHERE (candidate_bulk_density > eps8)
          candidate_k_bioturb = bioturb_numerator / candidate_bulk_density
        END WHERE

        DO transport_pool_index = 1, SIZE(bioturb_pools)
          bioturbation_pool = candidate_domain(bioturb_pools(transport_pool_index),ixC,:,:,1)
          unburied_surface = 0.0_wp
          SELECT CASE (bioturb_pools(transport_pool_index))
          CASE (ix_polymeric_litter, ix_microbial, ix_residue)
            unburied_surface = bioturbation_pool(:,1) * (1.0_wp - candidate_surface_burial(:,1))
            bioturbation_pool(:,1) = bioturbation_pool(:,1) * candidate_surface_burial(:,1)
          END SELECT

          CALL implicit_bioturbation_step(nc, nsoil_w, num_sl_above_bedrock, soil_depth_sl, &
            & candidate_k_bioturb, one_day * one_year, bioturbation_pool, transported_pool)
          candidate_domain(bioturb_pools(transport_pool_index),ixC,:,:,1) = transported_pool
          candidate_domain(bioturb_pools(transport_pool_index),ixC,:,1,1) = &
            & candidate_domain(bioturb_pools(transport_pool_index),ixC,:,1,1) + unburied_surface
        END DO
      END DO

      CALL update_candidate_soil_state(candidate_domain, candidate_surface_burial)
      max_pool_change = MAXVAL(ABS(candidate_domain(:,ixC,:,:,1) - pool_domain(:,ixC,:,:,1)))
      max_bulk_excess = 0.0_wp
      DO local_ic = 1, nc
        DO local_is = 1, NINT(num_sl_above_bedrock(local_ic))
          max_bulk_excess = MAX(max_bulk_excess, candidate_bulk_soil_c(local_ic,local_is) - &
            & bulk_carbon_max_domain(local_ic,local_is,1))
        END DO
      END DO

      pool_domain = candidate_domain
      surface_burial_domain(:,:,1) = candidate_surface_burial

      IF (l_transport_debug) THEN
        CALL write_bioturbation_target_profile('bioturbation-smoothed target', candidate_domain, &
          & bulk_carbon_max_domain(:,:,1), qmax_org_domain(:,:,1), candidate_surface_burial, &
          & enzyme_poly_domain(:,:,1), enzyme_res_domain(:,:,1), enzyme_ap_domain(:,:,1))
        WRITE(message_text,*) 'DEBUG analytical bioturbation target smoothing; years, ', &
          & 'max pool change [mol C m-3], max bulk-C excess [mol C m-3] = ', &
          & transport_smoothing_years, max_pool_change, max_bulk_excess
        CALL message(TRIM(modname)//':apply_analytical_ss_jsm_accelerator', TRIM(message_text))
      END IF
    END SUBROUTINE apply_bioturbation_target_smoothing

    SUBROUTINE write_bioturbation_target_profile(label, target_pool, target_bulk_max, target_qmax, &
      & target_burial, target_enzyme_poly, target_enzyme_res, target_enzyme_ap)
      CHARACTER(len=*), INTENT(in) :: label !< name of the analytical target being reported
      REAL(wp), INTENT(in) :: target_pool(:,:,:,:,:) !< target element pools [mol m-3]
      REAL(wp), INTENT(in) :: target_bulk_max(:,:) !< maximum admissible bulk C [mol C m-3]
      REAL(wp), INTENT(in) :: target_qmax(:,:) !< target organic sorption capacity [mol C m-3]
      REAL(wp), INTENT(in) :: target_burial(:,:) !< target first-layer burial fraction [-]
      REAL(wp), INTENT(in) :: target_enzyme_poly(:,:) !< polymeric enzyme allocation fraction [-]
      REAL(wp), INTENT(in) :: target_enzyme_res(:,:) !< residue enzyme allocation fraction [-]
      REAL(wp), INTENT(in) :: target_enzyme_ap(:,:) !< phosphatase enzyme allocation fraction [-]

      REAL(wp) :: total_carbon !< total target organic C density in one layer [mol C m-3]
      REAL(wp) :: profile_bulk_soil_c(nc,nsoil_w) !< diagnostic bulk soil C density [mol C m-3]
      REAL(wp) :: profile_litter_c(nc,nsoil_w) !< diagnostic litter C density [mol C m-3]
      REAL(wp) :: profile_volume_min(nc,nsoil_w) !< diagnostic mineral soil volume fraction [-]
      REAL(wp) :: profile_bulk_density(nc,nsoil_w) !< diagnostic corrected bulk density [kg m-3]
      REAL(wp) :: profile_burial(nc,nsoil_w) !< diagnostic burial fraction after constraint update [-]
      INTEGER  :: profile_ic   !< horizontal point index [-]
      INTEGER  :: profile_is   !< soil-layer index [-]

      profile_burial = target_burial
      profile_bulk_soil_c = 0.0_wp
      profile_litter_c = 0.0_wp
      profile_volume_min = 0.0_wp
      CALL calc_bulk_soil_carbon(nc, nsoil_w, num_sl_above_bedrock, 'jsm', &
        & target_pool(:,:,:,:,1), target_bulk_max, profile_burial, &
        & profile_bulk_soil_c, profile_litter_c, profile_volume_min)
      profile_bulk_density = calc_bulk_density_correction( &
        & profile_bulk_soil_c, profile_litter_c, bulk_dens_sl)

      DO profile_ic = 1, nc
        DO profile_is = 1, NINT(num_sl_above_bedrock(profile_ic))
          total_carbon = SUM(target_pool(:,ixC,profile_ic,profile_is,1))
          WRITE(message_text,*) 'DEBUG analytical ', TRIM(label), ' profile ic/is = ', profile_ic, profile_is, &
            & '; C_sol, C_poly, C_wood, C_DOM, C_aDOM, C_mic, C_res, C_aRes [mol C m-3] = ', &
            & target_pool(ix_soluable_litter,ixC,profile_ic,profile_is,1), &
            & target_pool(ix_polymeric_litter,ixC,profile_ic,profile_is,1), &
            & target_pool(ix_woody_litter,ixC,profile_ic,profile_is,1), &
            & target_pool(ix_dom,ixC,profile_ic,profile_is,1), &
            & target_pool(ix_dom_assoc,ixC,profile_ic,profile_is,1), &
            & target_pool(ix_microbial,ixC,profile_ic,profile_is,1), &
            & target_pool(ix_residue,ixC,profile_ic,profile_is,1), &
            & target_pool(ix_residue_assoc,ixC,profile_ic,profile_is,1), &
            & '; C_total, C_bulk, C_litter [mol C m-3] = ', total_carbon, &
            & profile_bulk_soil_c(profile_ic,profile_is), profile_litter_c(profile_ic,profile_is), &
            & '; C_bulk_max, qmax [mol C m-3], bulk_density [kg m-3], mineral_volume, burial = ', &
            & target_bulk_max(profile_ic,profile_is), target_qmax(profile_ic,profile_is), &
            & profile_bulk_density(profile_ic,profile_is), profile_volume_min(profile_ic,profile_is), &
            & profile_burial(profile_ic,profile_is), &
            & '; enzyme_poly, enzyme_res, enzyme_AP = ', target_enzyme_poly(profile_ic,profile_is), &
            & target_enzyme_res(profile_ic,profile_is), target_enzyme_ap(profile_ic,profile_is)
          CALL message(TRIM(modname)//':apply_analytical_ss_jsm_accelerator', TRIM(message_text))
        END DO
      END DO
    END SUBROUTINE write_bioturbation_target_profile


    SUBROUTINE update_candidate_soil_state(candidate_pool, candidate_burial)

      REAL(wp), INTENT(in)    :: candidate_pool(:,:,:,:,:) !< candidate element pools [mol m-3]
      REAL(wp), INTENT(inout) :: candidate_burial(:,:) !< burial fraction updated by bulk-C constraints [-]

      ! These host-associated scratch arrays are deliberately updated here so
      ! subsequent bioturbation smoothing uses a mutually consistent C density,
      ! mineral volume, corrected bulk density, and burial fraction.

      candidate_bulk_soil_c = 0.0_wp
      candidate_litter_c = 0.0_wp
      candidate_volume_min = 0.0_wp
      CALL calc_bulk_soil_carbon(nc, nsoil_w, num_sl_above_bedrock, 'jsm', &
        & candidate_pool(:,:,:,:,1), bulk_carbon_max_domain(:,:,1), candidate_burial, &
        & candidate_bulk_soil_c, candidate_litter_c, candidate_volume_min)
      candidate_bulk_density = calc_bulk_density_correction( &
        & candidate_bulk_soil_c, candidate_litter_c, bulk_dens_sl)
    END SUBROUTINE update_candidate_soil_state

    SUBROUTINE implicit_bioturbation_step(nc_local, nsoil_local, nsl, layer_thickness, &
      & diffusivity, step_seconds, old_pool, new_pool)

      INTEGER,  INTENT(in)  :: nc_local          !< number of horizontal points [-]
      INTEGER,  INTENT(in)  :: nsoil_local       !< number of soil layers [-]
      REAL(wp), INTENT(in)  :: nsl(:)            !< active soil-layer count per point [-]
      REAL(wp), INTENT(in)  :: layer_thickness(:,:) !< soil-layer thickness [m]
      REAL(wp), INTENT(in)  :: diffusivity(:,:)  !< bioturbation diffusion coefficient [m2 s-1]
      REAL(wp), INTENT(in)  :: step_seconds      !< implicit transport integration period [s]
      REAL(wp), INTENT(in)  :: old_pool(:,:)     !< pre-transport pool density [mol C m-3]
      REAL(wp), INTENT(out) :: new_pool(:,:)     !< post-transport pool density [mol C m-3]

      REAL(wp) :: lower(nsoil_local)      !< lower diagonal of dimensionless transport matrix [-]
      REAL(wp) :: diagonal(nsoil_local)   !< main diagonal of dimensionless transport matrix [-]
      REAL(wp) :: upper(nsoil_local)      !< upper diagonal of dimensionless transport matrix [-]
      REAL(wp) :: rhs(nsoil_local)        !< tridiagonal right-hand-side pool density [mol C m-3]
      REAL(wp) :: solution(nsoil_local)   !< tridiagonal solution pool density [mol C m-3]
      REAL(wp) :: interface_conductance   !< diffusive conductance across a layer interface [m s-1]
      REAL(wp) :: midpoint_distance       !< distance between adjacent layer midpoints [m]
      REAL(wp) :: elimination_factor      !< Thomas-algorithm elimination multiplier [-]
      INTEGER  :: ic_local                !< horizontal point index [-]
      INTEGER  :: is_local                !< soil-layer index [-]
      INTEGER  :: n_active                !< active soil-layer count for one point [-]

      new_pool = old_pool
      DO ic_local = 1, nc_local
        n_active = NINT(nsl(ic_local))
        IF (n_active <= 1) CYCLE

        lower = 0.0_wp
        diagonal = 1.0_wp
        upper = 0.0_wp
        rhs = 0.0_wp
        solution = 0.0_wp
        rhs(1:n_active) = old_pool(ic_local,1:n_active)

        DO is_local = 1, n_active - 1
          midpoint_distance = 0.5_wp * (layer_thickness(ic_local,is_local) + &
            & layer_thickness(ic_local,is_local + 1))
          IF (midpoint_distance <= eps8 .OR. layer_thickness(ic_local,is_local) <= eps8 .OR. &
            & layer_thickness(ic_local,is_local + 1) <= eps8) CYCLE
          interface_conductance = MAX(0.0_wp, diffusivity(ic_local,is_local)) / midpoint_distance
          diagonal(is_local) = diagonal(is_local) + step_seconds * interface_conductance / &
            & layer_thickness(ic_local,is_local)
          upper(is_local) = upper(is_local) - step_seconds * interface_conductance / &
            & layer_thickness(ic_local,is_local)
          lower(is_local + 1) = lower(is_local + 1) - step_seconds * interface_conductance / &
            & layer_thickness(ic_local,is_local + 1)
          diagonal(is_local + 1) = diagonal(is_local + 1) + step_seconds * interface_conductance / &
            & layer_thickness(ic_local,is_local + 1)
        END DO

        DO is_local = 2, n_active
          elimination_factor = lower(is_local) / diagonal(is_local - 1)
          diagonal(is_local) = diagonal(is_local) - elimination_factor * upper(is_local - 1)
          rhs(is_local) = rhs(is_local) - elimination_factor * rhs(is_local - 1)
        END DO

        solution(n_active) = rhs(n_active) / diagonal(n_active)
        DO is_local = n_active - 1, 1, -1
          solution(is_local) = (rhs(is_local) - upper(is_local) * solution(is_local + 1)) / diagonal(is_local)
        END DO
        new_pool(ic_local,1:n_active) = MAX(0.0_wp, solution(1:n_active))
      END DO
    END SUBROUTINE implicit_bioturbation_step

    SUBROUTINE reconstruct_target_elements(nc_local, nsoil_local, nsl, old_pool, target_pool)
      INTEGER,  INTENT(in)    :: nc_local       !< number of horizontal points [-]
      INTEGER,  INTENT(in)    :: nsoil_local    !< number of soil layers [-]
      REAL(wp), INTENT(in)    :: nsl(:)         !< active soil-layer count per point [-]
      REAL(wp), INTENT(in)    :: old_pool(:,:,:,:) !< pre-jump element pools [mol m-3]
      REAL(wp), INTENT(inout) :: target_pool(:,:,:,:) !< analytical C and reconstructed elements [mol m-3]
      REAL(wp) :: old_c       !< pre-jump pool C density [mol C m-3]
      REAL(wp) :: new_c       !< analytical target pool C density [mol C m-3]
      REAL(wp) :: layer_c     !< pre-jump total layer C density [mol C m-3]
      REAL(wp) :: layer_ratio(SIZE(old_pool,2)) !< layer-average element-to-C molar ratios [-]
      REAL(wp) :: ratio       !< selected pool or layer element-to-C molar ratio [-]
      INTEGER  :: ic          !< horizontal point index [-]
      INTEGER  :: is          !< soil-layer index [-]
      INTEGER  :: ipool       !< biogeochemical pool index [-]
      INTEGER  :: ielem       !< element index [-]

      ! Preserve each pool's existing stoichiometry where possible; fall back
      ! first to the layer-average ratio and finally to model-default ratios.

      DO ic = 1, nc_local
        DO is = NINT(nsl(ic)) + 1, nsoil_local
          target_pool(:,:,ic,is) = old_pool(:,:,ic,is)
        END DO
        DO is = 1, NINT(nsl(ic))
          layer_c = SUM(old_pool(:,ixC,ic,is))
          layer_ratio(:) = 0.0_wp
          IF (layer_c > eps8) THEN
            DO ielem = 1, SIZE(old_pool,2)
              IF (ielem /= ixC) layer_ratio(ielem) = SUM(old_pool(:,ielem,ic,is)) / layer_c
            END DO
          END IF

          DO ipool = 1, SIZE(old_pool,1)
            old_c = old_pool(ipool,ixC,ic,is)
            new_c = MAX(0.0_wp, target_pool(ipool,ixC,ic,is))
            target_pool(ipool,ixC,ic,is) = new_c
            DO ielem = 1, SIZE(old_pool,2)
              IF (ielem == ixC) CYCLE
              IF (old_c > eps8) THEN
                ratio = old_pool(ipool,ielem,ic,is) / old_c
              ELSE IF (layer_c > eps8) THEN
                ratio = layer_ratio(ielem)
              ELSE
                ratio = default_element_ratio(ipool, ielem)
              END IF
              target_pool(ipool,ielem,ic,is) = MAX(0.0_wp, new_c * ratio)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE reconstruct_target_elements

    SUBROUTINE reconstruct_target_c14(nc_local, nsoil_local, nsl, &
      & sol_c_in, poly_c_in, wood_c_in, sol_c14_in, poly_c14_in, wood_c14_in, &
      & doc_transport_in, desorption_mod_in, target_pool)
      INTEGER,  INTENT(in)    :: nc_local       !< number of horizontal points [-]
      INTEGER,  INTENT(in)    :: nsoil_local    !< number of soil layers [-]
      REAL(wp), INTENT(in)    :: nsl(:)         !< active soil-layer count per point [-]
      REAL(wp), INTENT(in)    :: sol_c_in(:,:)  !< soluble-litter C input [mol C m-3 s-1]
      REAL(wp), INTENT(in)    :: poly_c_in(:,:) !< polymeric-litter C input [mol C m-3 s-1]
      REAL(wp), INTENT(in)    :: wood_c_in(:,:) !< woody-litter C input [mol C m-3 s-1]
      REAL(wp), INTENT(in)    :: sol_c14_in(:,:)  !< soluble-litter C14 input [mol C14 m-3 s-1]
      REAL(wp), INTENT(in)    :: poly_c14_in(:,:) !< polymeric-litter C14 input [mol C14 m-3 s-1]
      REAL(wp), INTENT(in)    :: wood_c14_in(:,:) !< woody-litter C14 input [mol C14 m-3 s-1]
      REAL(wp), INTENT(in)    :: doc_transport_in(:,:) !< downward DOC transport coefficient [s-1]
      REAL(wp), INTENT(in)    :: desorption_mod_in(:,:) !< environmental desorption modifier [-]
      REAL(wp), INTENT(inout) :: target_pool(:,:,:,:) !< target element pools including C14 [mol m-3]

      REAL(wp) :: c_dom, c_mic, c_res, c_adom, c_ares !< target DOM/microbial/residue C pools [mol C m-3]
      REAL(wp) :: c_poly, c_sol, c_wood !< target litter C pools [mol C m-3]
      REAL(wp) :: d_litter              !< source rate to polymeric litter [mol C m-3 s-1]
      REAL(wp) :: d_residue             !< source rate to microbial residue [mol C m-3 s-1]
      REAL(wp) :: mortality_flux        !< microbial mortality C flux [mol C m-3 s-1]
      REAL(wp) :: doc_transport_above   !< DOC transport coefficient of layer above [s-1]
      REAL(wp) :: kdes_dom_eff          !< modified associated-DOM desorption rate [s-1]
      REAL(wp) :: kdes_res_eff          !< modified associated-residue desorption rate [s-1]
      REAL(wp) :: source_c              !< total source C rate for a target pool [mol C m-3 s-1]
      REAL(wp) :: source_c14            !< total source C14 rate for a target pool [mol C14 m-3 s-1]
      REAL(wp) :: ratio_dom, ratio_mic, ratio_res !< C14-to-C molar ratios of DOM, microbes and residue [-]
      REAL(wp) :: ratio_poly, ratio_sol, ratio_wood !< C14-to-C molar ratios of litter pools [-]
      REAL(wp) :: new_dom_c14           !< updated steady-state DOM C14 density [mol C14 m-3]
      REAL(wp) :: old_dom_c14           !< previous DOM C14 iterate [mol C14 m-3]
      REAL(wp) :: max_change            !< absolute DOM C14 iteration change [mol C14 m-3]
      INTEGER  :: ic                    !< horizontal point index [-]
      INTEGER  :: is                    !< soil-layer index [-]
      INTEGER  :: iter                  !< C14 fixed-point iteration index [-]
      INTEGER, PARAMETER :: max_c14_iterations = 100 !< maximum C14 closure iterations [-]

      IF (SUM(MAX(0.0_wp, sol_c14_in)) + SUM(MAX(0.0_wp, poly_c14_in)) + &
        & SUM(MAX(0.0_wp, wood_c14_in)) <= eps8) RETURN

      DO ic = 1, nc_local
        DO is = 1, NINT(nsl(ic))
          c_sol  = target_pool(ix_soluable_litter,ixC,ic,is)
          c_wood = target_pool(ix_woody_litter,ixC,ic,is)
          c_poly = target_pool(ix_polymeric_litter,ixC,ic,is)
          c_dom  = target_pool(ix_dom,ixC,ic,is)
          c_mic  = target_pool(ix_microbial,ixC,ic,is)
          c_res  = target_pool(ix_residue,ixC,ic,is)
          c_adom = target_pool(ix_dom_assoc,ixC,ic,is)
          c_ares = target_pool(ix_residue_assoc,ixC,ic,is)

          target_pool(ix_soluable_litter,ixC14,ic,is) = tracer_steady_state( &
            & MAX(0.0_wp, sol_c14_in(ic,is)), MAX(0.0_wp, sol_c_in(ic,is)), c_sol)
          target_pool(ix_woody_litter,ixC14,ic,is) = tracer_steady_state( &
            & MAX(0.0_wp, wood_c14_in(ic,is)), MAX(0.0_wp, wood_c_in(ic,is)), c_wood)

          d_litter = MAX(0.0_wp, poly_c_in(ic,is) + fc_woody2polymeric * wood_c_in(ic,is))
          ratio_wood = safe_tracer_ratio(target_pool(ix_woody_litter,ixC14,ic,is), c_wood)
          source_c14 = MAX(0.0_wp, poly_c14_in(ic,is)) + &
            & fc_woody2polymeric * MAX(0.0_wp, wood_c_in(ic,is)) * ratio_wood
          target_pool(ix_polymeric_litter,ixC14,ic,is) = tracer_steady_state(source_c14, d_litter, c_poly)

          mortality_flux = k_mort_microbial * c_mic
          d_residue = (1.0_wp - f_soluable_residue) * mortality_flux
          kdes_dom_eff = MAX(0.0_wp, k_desorpt_dom * desorption_mod_in(ic,is))
          kdes_res_eff = MAX(0.0_wp, k_desorpt_det * desorption_mod_in(ic,is))

          DO iter = 1, max_c14_iterations
            old_dom_c14 = target_pool(ix_dom,ixC14,ic,is)
            ratio_dom = safe_tracer_ratio(old_dom_c14, c_dom)

            target_pool(ix_microbial,ixC14,ic,is) = tracer_steady_state( &
              & mortality_flux * ratio_dom, mortality_flux, c_mic)
            ratio_mic = safe_tracer_ratio(target_pool(ix_microbial,ixC14,ic,is), c_mic)

            target_pool(ix_residue,ixC14,ic,is) = tracer_steady_state( &
              & d_residue * ratio_mic, d_residue, c_res)
            ratio_res = safe_tracer_ratio(target_pool(ix_residue,ixC14,ic,is), c_res)

            target_pool(ix_dom_assoc,ixC14,ic,is) = tracer_steady_state( &
              & kdes_dom_eff * c_adom * ratio_dom, kdes_dom_eff * c_adom, c_adom)
            target_pool(ix_residue_assoc,ixC14,ic,is) = tracer_steady_state( &
              & kdes_res_eff * c_ares * ratio_res, kdes_res_eff * c_ares, c_ares)

            ratio_sol  = safe_tracer_ratio(target_pool(ix_soluable_litter,ixC14,ic,is), c_sol)
            ratio_poly = safe_tracer_ratio(target_pool(ix_polymeric_litter,ixC14,ic,is), c_poly)
            source_c = fc_soluable2dom * MAX(0.0_wp, sol_c_in(ic,is)) + d_litter &
              & + f_soluable_residue * mortality_flux + d_residue &
              & + kdes_dom_eff * c_adom + kdes_res_eff * c_ares
            source_c14 = fc_soluable2dom * MAX(0.0_wp, sol_c_in(ic,is)) * ratio_sol &
              & + d_litter * ratio_poly + f_soluable_residue * mortality_flux * ratio_mic &
              & + d_residue * ratio_res &
              & + kdes_dom_eff * target_pool(ix_dom_assoc,ixC14,ic,is) &
              & + kdes_res_eff * target_pool(ix_residue_assoc,ixC14,ic,is)

            IF (is > 1) THEN
              doc_transport_above = MAX(0.0_wp, doc_transport_in(ic,is-1))
              source_c = source_c + doc_transport_above * target_pool(ix_dom,ixC,ic,is-1)
              source_c14 = source_c14 + doc_transport_above * target_pool(ix_dom,ixC14,ic,is-1)
            END IF

            new_dom_c14 = tracer_steady_state(source_c14, source_c, c_dom)
            target_pool(ix_dom,ixC14,ic,is) = new_dom_c14
            max_change = ABS(new_dom_c14 - old_dom_c14)
            IF (max_change <= eps8 * MAX(1.0_wp, ABS(new_dom_c14))) EXIT
          END DO
        END DO
      END DO
    END SUBROUTINE reconstruct_target_c14

    PURE REAL(wp) FUNCTION tracer_steady_state(source_tracer, source_carbon, pool_carbon) RESULT(pool_tracer)
      REAL(wp), INTENT(in) :: source_tracer !< tracer source rate [mol tracer m-3 s-1]
      REAL(wp), INTENT(in) :: source_carbon !< total C source rate [mol C m-3 s-1]
      REAL(wp), INTENT(in) :: pool_carbon   !< prescribed target C density [mol C m-3]
      ! pool_tracer is the decay-corrected steady-state tracer density [mol tracer m-3].

      ! The carbon target defines the implied gross turnover rate source_carbon / pool_carbon.
      IF (source_carbon > eps8 .AND. pool_carbon > eps8) THEN
        pool_tracer = MAX(0.0_wp, source_tracer) * pool_carbon / &
          & (source_carbon + lambda_C14 * pool_carbon)
      ELSE
        pool_tracer = 0.0_wp
      END IF
    END FUNCTION tracer_steady_state

    PURE REAL(wp) FUNCTION safe_tracer_ratio(pool_tracer, pool_carbon) RESULT(ratio)
      REAL(wp), INTENT(in) :: pool_tracer !< tracer pool density [mol tracer m-3]
      REAL(wp), INTENT(in) :: pool_carbon !< carbon pool density [mol C m-3]
      ! ratio is the non-negative tracer-to-carbon molar ratio [-].

      IF (pool_carbon > eps8) THEN
        ratio = MAX(0.0_wp, pool_tracer) / pool_carbon
      ELSE
        ratio = 0.0_wp
      END IF
    END FUNCTION safe_tracer_ratio

    REAL(wp) FUNCTION default_element_ratio(ipool, ielem) RESULT(ratio)
      INTEGER, INTENT(in) :: ipool !< biogeochemical pool identifier [-]
      INTEGER, INTENT(in) :: ielem !< requested element identifier [-]
      REAL(wp) :: cn_ratio         !< default pool C:N molar ratio [-]
      REAL(wp) :: np_ratio         !< default pool N:P molar ratio [-]
      ! ratio is the requested element-to-C molar ratio [-].

      SELECT CASE(ipool)
      CASE(ix_soluable_litter)
        cn_ratio = k_init_soluable_litter_cn
        np_ratio = k_init_litter_np
      CASE(ix_polymeric_litter)
        cn_ratio = k_init_polymeric_litter_cn
        np_ratio = k_init_litter_np
      CASE(ix_woody_litter)
        cn_ratio = k_init_woody_litter_cn
        np_ratio = k_init_litter_np
      CASE(ix_residue, ix_residue_assoc)
        cn_ratio = microbial_cn
        np_ratio = k_slow_som_np
      CASE DEFAULT
        cn_ratio = microbial_cn
        np_ratio = microbial_np
      END SELECT

      ratio = 0.0_wp
      IF (ielem == ixN)  ratio = 1.0_wp / cn_ratio
      IF (ielem == ixP)  ratio = 1.0_wp / cn_ratio / np_ratio
      IF (ielem == ixN15) ratio = 0.0_wp
    END FUNCTION default_element_ratio

  END SUBROUTINE apply_analytical_ss_jsm_accelerator
#endif
END MODULE mo_sb_analytical_accelerator
