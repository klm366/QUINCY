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

  USE mo_kind,               ONLY: lp
  USE mo_lnd_bgcm_idx
  USE mo_sb_constants, ONLY: min_fine_mineral_volume

  IMPLICIT NONE

  PRIVATE
  ! Scientific targets used only by the shared replay event driver.
  PUBLIC :: t_sb_local_coefficients, t_sb_local_target, solve_local_carbon_target
  PUBLIC :: solve_analytical_maoc_partition, reconstruct_analytical_column
  PUBLIC :: t_sb_analytical_inputs, build_sb_analytical_inputs, diagnose_local_carbon_targets
  PUBLIC :: solve_density_dependent_microbial_target
  PUBLIC :: solve_simple_carbon_target
  PUBLIC :: t_sb_monthly_coefficients, t_sb_column_residual, t_sb_column_solution
  PUBLIC :: evaluate_analytical_bioturbation, prepare_analytical_bioturbation_inputs
  PUBLIC :: evaluate_analytical_column_residual, evaluate_analytical_monthly_residual
  PUBLIC :: solve_analytical_carbon_column, solve_analytical_carbon_profile
  PUBLIC :: diagnose_full_profile_carbon_targets

  !> Local pool order is S,W,L,B,R,D,AD,AR; concentrations refer to TOTAL soil volume.
  !> Coefficients are prescribed, not additional fitted parameters or nutrient-equilibrium closures.
  TYPE :: t_sb_local_coefficients
    REAL(lp) :: dtime = 0.0_lp !< Forward replay timestep [s]
    REAL(lp) :: litter_input_carbon(3) = 0.0_lp !< Soluble, woody, polymeric litter inputs [mol C m-3 s-1]
    REAL(lp) :: transport_dom_inflow = 0.0_lp, dom_export_rate = 0.0_lp !< Target-derived DOM influx [mol C m-3 s-1], export [s-1]
    REAL(lp) :: litter_turnover_rate(2) = 0.0_lp !< Effective soluble and woody decay constants [s-1]
    REAL(lp) :: carbon_transfer_fraction(3) = 0.0_lp !< Soluble-to-DOM, woody-to-poly, mortality-to-DOM fractions [-]
    REAL(lp) :: microbial_decay_rate = 0.0_lp, minimum_microbial_carbon = 0.0_lp !< Mortality [s-1], floor [mol C m-3]
    REAL(lp) :: mortality_exponent = 1._lp !< beta=1 linear regression; beta>1 density dependence [-]
    REAL(lp) :: mortality_reference = 1._lp !< Reference biomass per TOTAL soil volume [mol C m-3]
    REAL(lp) :: microbial_cue = 0.0_lp, enzyme_frac_ap_layer = 0.0_lp !< Prescribed CUE and AP allocation [-]
    REAL(lp) :: alpha_poly = 0.0_lp !< Historical allocation for input validation only; target solves its own [-]
    REAL(lp) :: vmax_uptake_eff = 0.0_lp, km_uptake_eff = 0.0_lp !< Molar uptake capacity [s-1], HSC [mol C m-3]
    REAL(lp) :: vmax_depoly_eff(2) = 0.0_lp, km_depoly_eff = 0.0_lp !< Poly/residue capacity [s-1], HSC [mol C m-3]
    REAL(lp) :: enzyme_frac_total_local = 0.0_lp, persistent_enzyme_carbon = 0.0_lp !< Allocation [-], persistent C [mol C m-3]
    REAL(lp) :: adsorption_rate(2) = 0.0_lp !< DOM/residue adsorption [m3 mol C-1 s-1]
    REAL(lp) :: desorption_rate(2) = 0.0_lp !< Associated DOM/residue desorption [s-1]
    REAL(lp) :: qmax_org_layer = 0.0_lp, volume_pore_total = 0.0_lp !< Chemical capacity [mol C m-3], fixed pore fraction [-]
    REAL(lp) :: volume_per_carbon_maoc = 0.0_lp, volume_per_carbon_poc = 0.0_lp !< Solid volume per C concentration [m3 mol C-1]
    LOGICAL :: l_surface_layer = .FALSE. !< Surface excess POC can remain above the mineral soil [-]
  END TYPE t_sb_local_coefficients

  TYPE :: t_sb_monthly_coefficients
    REAL(lp), ALLOCATABLE :: weight(:) !< Fractions of replay year, nonnegative and summing to one [-]
    TYPE(t_sb_local_coefficients), ALLOCATABLE :: layer(:,:) !< Layer/bin conditional coefficients
  END TYPE t_sb_monthly_coefficients

  !> A valid carbon result still requires column-wide organic element reconstruction before a jump.
  TYPE :: t_sb_local_target
    LOGICAL :: l_minimum_microbial_seed = .FALSE. !< Quasi-steady floor branch; microbial growth retained [-]
    INTEGER :: status = 1 !< 0 conditional seed; 1 invalid; 2 no root; 3 limiter; 4 geometry; 5 residual [-]
    CHARACTER(len=192) :: reason = 'Not evaluated' !< Rejection or conditional-success explanation [-]
    REAL(lp) :: carbon(8) = 0.0_lp !< S,W,L,B,R,D,AD,AR candidate pools [mol C m-3]
    REAL(lp) :: residual(8) = 0.0_lp !< Same-order local carbon balance residuals [mol C m-3 s-1]
    REAL(lp) :: burial = 1.0_lp, maoc_volume = 0.0_lp, poc_volume = 0.0_lp !< Candidate geometry [-]
    REAL(lp) :: poc_excess = 0.0_lp !< Solid POC volume outside fixed pores [-]
    REAL(lp) :: transport(3) = 0.0_lp !< DOM influx, downward export, net lateral export [mol C m-3 s-1]
    REAL(lp) :: uptake_step = 0.0_lp !< Unrestricted DOC uptake over one replay timestep [mol C m-3]
    REAL(lp) :: microbial_mortality = 0._lp !< Candidate mortality flux, zero at biomass floor [mol C m-3 s-1]
    REAL(lp) :: alpha_poly = 0.0_lp, alpha_residue = 0.0_lp !< Solved carbon allocation fractions [-]
    REAL(lp) :: enzyme_depoly = 0.0_lp !< Total effective depolymerizing enzyme C [mol C m-3]
    REAL(lp) :: enzyme_poly = 0.0_lp, enzyme_residue = 0.0_lp !< Effective substrate-specific enzyme C [mol C m-3]
  END TYPE t_sb_local_target



  !> Read-only representative inputs, rebuilt from replay forcing at each eligible event.
  !> Mean products of modifiers are retained: product(mean(rt),mean(rm)) is NOT equivalent.
  TYPE :: t_sb_analytical_inputs
    LOGICAL :: ready = .FALSE. !< Complete and valid forcing available [-]
    REAL(lp) :: timestep = 0.0_lp !< Replay timestep used to partition litter [s]
    REAL(lp), ALLOCATABLE :: month_weight(:) !< Replay day-count weights, sum equals one [-]
    REAL(lp), ALLOCATABLE :: litter_mean(:,:,:,:) !< pool/element/cell/layer mean input [mol m-3 s-1]
    REAL(lp), ALLOCATABLE :: litter_month(:,:,:,:,:) !< pool/element/cell/layer/month input [mol m-3 s-1]
    REAL(lp), ALLOCATABLE :: modifier_month(:,:,:,:) !< cell/layer/modifier/month kinetics [-]
    REAL(lp), ALLOCATABLE :: rate_product_mean(:,:,:) !< depoly/uptake/sorption/desorption/HSC/SSM rt*rm [-]
    REAL(lp), ALLOCATABLE :: cue(:,:), alpha_poly(:,:), alpha_ap(:,:) !< Coupled-window means [-]
    REAL(lp), ALLOCATABLE :: water_rate_mean(:,:,:) !< Downward/lateral hydrological export rates [s-1]
    REAL(lp), ALLOCATABLE :: litter_budget_residual(:,:) !< Partitioned minus source element rate [mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: depth(:,:), pore_volume(:,:) !< Thickness [m], fixed pore fraction [-]
    REAL(lp), ALLOCATABLE :: qmax_org(:,:) !< Chemical sorption capacity [mol C m-3 soil]
    LOGICAL :: bioturbation_ready = .FALSE. !< Validated roots and geometry, not permission to commit [-]
    REAL(lp), ALLOCATABLE :: nactive(:) !< Active-layer count including partial bottom layer [-]
    REAL(lp), ALLOCATABLE :: root_fraction_mean(:,:) !< Annual mean root fractions, not renormalized [-]
    REAL(lp), ALLOCATABLE :: bulk_density(:,:) !< Entry bulk density, diagnostic only [kg m-3]
    REAL(lp), ALLOCATABLE :: bioturbation_entry(:,:) !< Entry diffusivity, not frozen target K [m2 s-1]
    REAL(lp), ALLOCATABLE :: interface_distance(:,:) !< Layer-centre spacing [m]
  END TYPE t_sb_analytical_inputs

  TYPE :: t_sb_column_residual
    LOGICAL :: valid = .FALSE. !< All active candidates evaluated [-]
    CHARACTER(len=192) :: reason = 'Not evaluated' !< Failure or approximation description [-]
    REAL(lp), ALLOCATABLE :: reaction(:,:,:), bioturbation(:,:,:), total(:,:,:) !< pool/cell/layer [mol C/m3/s]
    REAL(lp), ALLOCATABLE :: dom_transport(:,:), dom_interface(:,:) !< [mol C/m3/s], cell/interface [mol C/m2/s]
    REAL(lp), ALLOCATABLE :: bio_interface(:,:,:) !< pool/cell/interface, positive downward [mol C/m2/s]
    REAL(lp), ALLOCATABLE :: density(:,:), diffusivity(:,:), burial(:,:), poc_excess(:,:) !< [kg/m3], [m2/s], [-], [-]
    REAL(lp), ALLOCATABLE :: allocation_residual(:,:) !< alpha minus carbon allocation update [-]
    REAL(lp), ALLOCATABLE :: process_rate(:,:,:) !< M,U,dL,dR,JD,JR / cell / layer [mol C/m3/s]
    REAL(lp), ALLOCATABLE :: throughput(:,:,:) !< Gross per-pool source+sink for unchanged tolerances [mol C/m3/s]
    REAL(lp), ALLOCATABLE :: process_month(:,:,:), reaction_month(:,:,:) !< Process or pool/layer/bin [mol C/m3/s]
    REAL(lp), ALLOCATABLE :: allocation_month(:,:) !< Monthly carbon allocation residual [layer/bin, -]
    LOGICAL, ALLOCATABLE :: allocation_active_month(:,:) !< Monthly carbon allocation closure active [-]
    LOGICAL, ALLOCATABLE :: caps_month(:,:,:) !< layer/bin/warning category [-]
    REAL(lp), ALLOCATABLE :: mean_process(:,:,:), mean_total(:,:,:) !< Annual-coefficient counterfactual rates
    LOGICAL, ALLOCATABLE :: allocation_defined(:,:), geometry_admissible(:,:) !< Closure active; no deep overflow [-]
    LOGICAL, ALLOCATABLE :: timestep_warning(:,:,:) !< cell/layer: mortality, DOC uptake, other donor limits [-]
    REAL(lp), ALLOCATABLE :: column_budget(:,:) !< input, respiration, export, storage, closure error [mol C/m2/s]
  END TYPE t_sb_column_residual

  !> One attempted continuation stage; failed trials never replace the last successful seed.
  TYPE :: t_sb_continuation_stage
    REAL(lp) :: fraction = 0._lp, increment = 0._lp !< Bioturbation multiplier and attempted increment [-]
    REAL(lp) :: monthly_fraction = 0._lp !< Annual-coefficient to monthly-flux continuation multiplier [-]
    LOGICAL :: converged = .FALSE., admissible = .FALSE. !< Equation and physical checks [-]
    INTEGER :: iterations = 0, evaluations = 0 !< Accepted updates and residual calls [-]
    INTEGER :: worst_layer = 0, worst_equation = 0 !< Equation index S,W,L,B,R,D,AD,AR,alpha [-]
    REAL(lp), ALLOCATABLE :: norm_history(:), step_history(:) !< Normalized residual and accepted step [-]
    CHARACTER(len=192) :: reason = 'Not attempted' !< Stage outcome [-]
  END TYPE t_sb_continuation_stage

  !> Conditional stationary profile; no caller may infer permission to commit from finite pools alone.
  TYPE :: t_sb_column_solution
    LOGICAL :: converged = .FALSE. !< All pool and allocation equations passed numerical tolerances [-]
    LOGICAL :: admissible = .FALSE. !< Converged, geometry and selected finite-step guards passed [-]
    CHARACTER(len=192) :: reason = 'Not solved' !< Outcome; failed candidates remain diagnostic only [-]
    INTEGER :: iterations = 0, evaluations = 0 !< Accepted updates and residual calls [-]
    REAL(lp), ALLOCATABLE :: carbon(:,:), alpha_poly(:) !< Last candidate [mol C/m3], allocation [-]
    REAL(lp), ALLOCATABLE :: norm_history(:), step_history(:) !< Max normalized residual, accepted step [-]
    INTEGER :: worst_layer = 0, worst_equation = 0 !< Largest final normalized equation [-]
    REAL(lp) :: bioturbation_fraction = 0._lp !< Multiplier of diagnostic candidate; one required when enabled [-]
    REAL(lp) :: monthly_fraction = 0._lp !< Monthly nonlinear evaluation must reach one when supplied [-]
    INTEGER :: total_iterations = 0 !< Accepted updates across all continuation attempts [-]
    REAL(lp), ALLOCATABLE :: seed_carbon(:,:), seed_alpha(:) !< Previous analytical seed [mol C/m3], [-]
    INTEGER, ALLOCATABLE :: seed_status(:) !< Original local-solver status, not final acceptance [-]
    TYPE(t_sb_continuation_stage), ALLOCATABLE :: continuation(:) !< All attempted stages, including failures
    TYPE(t_sb_column_residual) :: diagnostics !< Final residual, fluxes, geometry and cap warnings
  END TYPE t_sb_column_solution

CONTAINS

  !> Solve nine equations per active layer: eight carbon balances and enzyme allocation closure.
  !! Newton steps use scaled finite differences, bounded line searches and QR fallback.
  !! Trial states are scratch only. Pools are mol C/m3 TOTAL soil; rates are mol C/m3/s.
  !! Residual tolerance is 1e-18 + 1e-8*gross throughput, allocation tolerance 1e-8.
  !! Candidate density changes diffusivity at every residual evaluation. Finite-step
  !! mortality/other-donor caps and deep pore overflow reject the final target;
  !! the DOC uptake cap is diagnostic. CUE and AP are prescribed annual means.
  SUBROUTINE solve_analytical_carbon_column(ns, thickness, root_fraction_mean, water_rate, lateral_return, &
    & coefficients, initial_carbon, initial_alpha, solution, max_iterations, include_bioturbation, &
    & bioturbation_fraction, monthly, monthly_fraction)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    USE mo_sb_constants, ONLY: min_enzyme_fraction
    INTEGER, INTENT(in) :: ns !< Number of active layers, including partial bottom layer [-]
    REAL(lp), INTENT(in) :: thickness(ns), root_fraction_mean(ns), water_rate(ns,2), lateral_return
    !< Thickness [m], root fractions [-], effective water export [s-1], return fraction [-]
    TYPE(t_sb_local_coefficients), INTENT(in) :: coefficients(ns) !< Fixed conditional kinetics/geometry
    REAL(lp), INTENT(in) :: initial_carbon(8,ns), initial_alpha(ns) !< [mol C/m3 total soil], [-]
    TYPE(t_sb_column_solution), INTENT(out) :: solution !< Diagnostic candidate; no live-state transaction
    INTEGER, OPTIONAL, INTENT(in) :: max_iterations !< Default 80 accepted iterations; zero evaluates seed only [-]
    LOGICAL, OPTIONAL, INTENT(in) :: include_bioturbation !< Default true; analytical transport only [-]
    REAL(lp), OPTIONAL, INTENT(in) :: bioturbation_fraction !< Continuation lambda in [0,1], default one [-]
    TYPE(t_sb_monthly_coefficients), OPTIONAL, INTENT(in) :: monthly !< Monthly nonlinear rate coefficients
    REAL(lp), OPTIONAL, INTENT(in) :: monthly_fraction !< Monthly-flux continuation in [0,1] [-]
    TYPE(t_sb_column_residual) :: current, trial !< Scratch residual evaluations
    REAL(lp), ALLOCATABLE :: x(:), trial_x(:), variable_scale(:), lower(:), upper(:)
    REAL(lp), ALLOCATABLE :: residual(:), trial_residual(:), jacobian(:,:), direction(:)
    REAL(lp), ALLOCATABLE :: matrix(:,:), rhs(:), equation_scale(:)
    REAL(lp) :: delta, step, merit, trial_merit, damping, max_norm !< Dimensionless solve quantities
    INTEGER :: limit, size_system, iteration, i, j, s, offset, attempt, backtrack, worst
    LOGICAL :: valid, linear_ok, improved !< Evaluation, linear solve and line-search status [-]

    solution%reason = 'Invalid full-profile initial state or iteration budget'
    IF (PRESENT(monthly)) solution%monthly_fraction = 1._lp
    IF (PRESENT(monthly_fraction)) solution%monthly_fraction = monthly_fraction
    IF (.NOT. ieee_is_finite(solution%monthly_fraction)) RETURN
    IF (solution%monthly_fraction < 0._lp .OR. solution%monthly_fraction > 1._lp) RETURN
    IF (solution%monthly_fraction > 0._lp .AND. .NOT. PRESENT(monthly)) RETURN
    solution%bioturbation_fraction = 1._lp
    IF (PRESENT(bioturbation_fraction)) solution%bioturbation_fraction = bioturbation_fraction
    IF (PRESENT(include_bioturbation)) THEN
      IF (.NOT. include_bioturbation) solution%bioturbation_fraction = 0._lp
    END IF
    IF (.NOT. ieee_is_finite(solution%bioturbation_fraction)) RETURN
    IF (solution%bioturbation_fraction < 0._lp .OR. solution%bioturbation_fraction > 1._lp) RETURN
    limit = 80
    IF (PRESENT(max_iterations)) limit = max_iterations
    IF (ns <= 0 .OR. limit < 0) RETURN
    IF (.NOT. ALL(ieee_is_finite(initial_carbon)) .OR. .NOT. ALL(ieee_is_finite(initial_alpha))) RETURN
    IF (ANY(initial_carbon < 0._lp)) RETURN
    IF (ANY(initial_alpha < min_enzyme_fraction) .OR. ANY(initial_alpha > 1._lp-min_enzyme_fraction)) RETURN
    DO s = 1,ns
      IF (.NOT. ieee_is_finite(coefficients(s)%minimum_microbial_carbon)) RETURN
      IF (initial_carbon(4,s) < coefficients(s)%minimum_microbial_carbon) RETURN
    END DO
    size_system = 9*ns
    ALLOCATE(x(size_system), trial_x(size_system), variable_scale(size_system), &
      & lower(size_system), upper(size_system), residual(size_system), trial_residual(size_system), &
      & jacobian(size_system,size_system), direction(size_system), matrix(size_system,size_system), &
      & rhs(size_system), equation_scale(size_system))
    ALLOCATE(solution%norm_history(0:limit), SOURCE=-1._lp)
    ALLOCATE(solution%step_history(limit), SOURCE=0._lp)
    DO s = 1,ns
      offset = 9*(s-1)
      ! Scale to the supplied profile without raising the microbial floor or altering any pool.
      variable_scale(offset+1:offset+8) = MAX(initial_carbon(:,s),1.e-12_lp)
      variable_scale(offset+9) = 1._lp
      x(offset+1:offset+8) = initial_carbon(:,s)/variable_scale(offset+1:offset+8)
      x(offset+9) = initial_alpha(s)
      lower(offset+1:offset+8) = 0._lp
      lower(offset+4) = coefficients(s)%minimum_microbial_carbon/variable_scale(offset+4)
      lower(offset+9) = min_enzyme_fraction
      ! Only a floating-point safety bound, not an ecological maximum pool size.
      upper(offset+1:offset+8) = SQRT(HUGE(1._lp))/MAX(1._lp,variable_scale(offset+1:offset+8))
      upper(offset+9) = 1._lp-min_enzyme_fraction
    END DO
    CALL evaluate(x,current,residual,valid,natural_scale=equation_scale)
    IF (.NOT. valid) THEN
      solution%reason = 'Invalid initial column residual: '//TRIM(current%reason)
      RETURN
    END IF
    DO iteration = 0,limit
      max_norm = MAXVAL(ABS(residual))
      solution%norm_history(iteration) = max_norm
      solution%iterations = iteration
      solution%total_iterations = iteration
      worst = MAXLOC(ABS(residual),DIM=1)
      solution%worst_layer = (worst-1)/9+1
      solution%worst_equation = MOD(worst-1,9)+1
      IF (max_norm <= 1._lp) THEN
        solution%converged = .TRUE.
        EXIT
      END IF
      IF (iteration == limit) THEN
        solution%reason = 'Full-profile iteration limit; no jump'
        EXIT
      END IF
      merit = SUM(residual**2)/REAL(size_system,lp)
      ! One-sided perturbation, reversed at an upper bound; the original iterate is never modified.
      DO j = 1,size_system
        delta = 1.e-5_lp*MAX(ABS(x(j)),1.e-4_lp)
        IF (x(j)+delta > upper(j)) delta = -delta
        delta = MAX(lower(j)-x(j),MIN(upper(j)-x(j),delta))
        IF (delta == 0._lp) EXIT
        trial_x = x
        trial_x(j) = x(j)+delta
        CALL evaluate(trial_x,trial,trial_residual,valid,fixed_scale=equation_scale)
        IF (.NOT. valid) EXIT
        jacobian(:,j) = (trial_residual-residual)/delta
      END DO
      IF (j <= size_system) THEN
        solution%reason = 'Invalid finite-difference Jacobian; no jump'
        EXIT
      END IF
      IF (.NOT. ALL(ieee_is_finite(jacobian))) THEN
        solution%reason = 'Nonfinite finite-difference Jacobian; no jump'
        EXIT
      END IF
      improved = .FALSE.
      DO attempt = 0,7
        IF (attempt == 0) THEN
          matrix = jacobian
          rhs = -residual
        ELSE
          ! Solve the augmented least-squares system by QR; J^T J would square its condition number.
          damping = 10._lp**(attempt-4)
          CALL solve_column_regularized_step(jacobian,residual,damping,direction,linear_ok)
        END IF
        IF (attempt == 0) CALL solve_column_linear_system(matrix,rhs,direction,linear_ok)
        IF (.NOT. linear_ok) CYCLE
        ! Project directions at active bounds; do not project a whole candidate onto a different root.
        WHERE (x <= lower .AND. direction < 0._lp) direction = 0._lp
        WHERE (x >= upper .AND. direction > 0._lp) direction = 0._lp
        ! Bound/trust-limit each coordinate BEFORE backtracking. A near-floor deep pool must not
        ! multiply every upper-layer update by its tiny distance to a bound. The merit test below
        ! still decides whether the modified direction improves the unmodified scientific equations.
        DO i = 1,size_system
          direction(i) = MAX(-2._lp*MAX(ABS(x(i)),1._lp), &
            & MIN(2._lp*MAX(ABS(x(i)),1._lp),direction(i)))
          IF (direction(i) < 0._lp) direction(i) = MAX(direction(i),0.995_lp*(lower(i)-x(i)))
          IF (direction(i) > 0._lp) direction(i) = MIN(direction(i),0.995_lp*(upper(i)-x(i)))
        END DO
        step = 1._lp
        DO backtrack = 1,24
          trial_x = MAX(lower,MIN(upper,x+step*direction)) ! Roundoff at a feasible endpoint only.
          CALL evaluate(trial_x,trial,trial_residual,valid,fixed_scale=equation_scale)
          IF (valid) THEN
            trial_merit = SUM(trial_residual**2)/REAL(size_system,lp)
            IF (trial_merit < merit*(1._lp-1.e-4_lp*step)) THEN
              improved = .TRUE.
              EXIT
            END IF
          END IF
          step = step*0.5_lp
        END DO
        IF (improved) EXIT
      END DO
      IF (.NOT. improved) THEN
        solution%reason = 'Full-profile line search stalled; no jump'
        EXIT
      END IF
      x = trial_x
      ! Acceptance always uses the original candidate-dependent throughput tolerances.
      ! Freeze denominators only within one Jacobian/line search, not across accepted iterations.
      CALL evaluate(x,current,residual,valid,natural_scale=equation_scale)
      IF (.NOT. valid) THEN
        solution%reason = 'Invalid accepted column residual; no jump'
        EXIT
      END IF
      solution%step_history(iteration+1) = step
    END DO
    ALLOCATE(solution%carbon(8,ns),solution%alpha_poly(ns))
    DO s = 1,ns
      offset = 9*(s-1)
      solution%carbon(:,s) = x(offset+1:offset+8)*variable_scale(offset+1:offset+8)
      solution%alpha_poly(s) = x(offset+9)
    END DO
    solution%diagnostics = current
    IF (.NOT. solution%converged) RETURN
    solution%admissible = ALL(current%geometry_admissible) .AND. &
      & .NOT. ANY(current%timestep_warning(:,:,1)) .AND. .NOT. ANY(current%timestep_warning(:,:,3))
    solution%admissible = solution%admissible .AND. ABS(current%column_budget(5,1)) <= &
      & 1.e-18_lp+1.e-10_lp*SUM(ABS(current%column_budget(1:3,1)))
    solution%reason = 'Conditional full-profile root; CUE/AP prescribed; DOC cap warnings remain diagnostic'
    IF (.NOT. ALL(current%allocation_defined)) &
      & solution%reason = 'Conditional root with inactive allocation retained from initial history'
    IF (.NOT. solution%admissible) &
      & solution%reason = 'Residual converged but geometry, finite-step or column-closure guards failed; no jump'

  CONTAINS
    !> All coefficients/densities/fluxes are reevaluated, including during Jacobian perturbations.
    SUBROUTINE evaluate(candidate, diagnostics, normalized, ok, fixed_scale, natural_scale)
      REAL(lp), INTENT(in) :: candidate(size_system) !< Scaled unknowns [-]
      TYPE(t_sb_column_residual), INTENT(out) :: diagnostics
      REAL(lp), INTENT(out) :: normalized(size_system) !< Residual divided by numerical allowance [-]
      LOGICAL, INTENT(out) :: ok
      REAL(lp), OPTIONAL, INTENT(in) :: fixed_scale(size_system) !< Frozen search denominators, not tolerances
      REAL(lp), OPTIONAL, INTENT(out) :: natural_scale(size_system) !< Original candidate-based allowances
      REAL(lp) :: pools(8,1,ns), alpha(1,ns), allowance(8)
      INTEGER :: layer, base
      DO layer = 1,ns
        base = 9*(layer-1)
        pools(:,1,layer) = candidate(base+1:base+8)*variable_scale(base+1:base+8)
        alpha(1,layer) = candidate(base+9)
      END DO
      IF (PRESENT(monthly)) THEN
        CALL evaluate_analytical_monthly_residual(ns,thickness,root_fraction_mean,water_rate,lateral_return, &
          & coefficients,pools,alpha,monthly,solution%monthly_fraction,solution%bioturbation_fraction,diagnostics)
      ELSE
        CALL evaluate_analytical_column_residual(1,ns,[REAL(ns,lp)],RESHAPE(thickness,[1,ns]), &
          & RESHAPE(root_fraction_mean,[1,ns]),RESHAPE(water_rate,[1,ns,2]),lateral_return, &
          & RESHAPE(coefficients,[1,ns]),pools,alpha,diagnostics,include_bioturbation, &
          & solution%bioturbation_fraction)
      END IF
      solution%evaluations = solution%evaluations+1
      ok = diagnostics%valid
      normalized = 0._lp
      IF (.NOT. ok) RETURN
      DO layer = 1,ns
        base = 9*(layer-1)
        allowance = 1.e-18_lp+1.e-8_lp*diagnostics%throughput(:,1,layer)
        IF (PRESENT(natural_scale)) THEN
          natural_scale(base+1:base+8) = allowance
          natural_scale(base+9) = 1.e-8_lp
        END IF
        IF (PRESENT(fixed_scale)) allowance = fixed_scale(base+1:base+8)
        normalized(base+1:base+8) = diagnostics%total(:,1,layer)/allowance
        IF (diagnostics%allocation_defined(1,layer)) THEN
          normalized(base+9) = diagnostics%allocation_residual(1,layer)/1.e-8_lp
        ELSE
          normalized(base+9) = (alpha(1,layer)-initial_alpha(layer))/1.e-8_lp
        END IF
      END DO
      ok = ALL(ieee_is_finite(normalized))
    END SUBROUTINE evaluate
  END SUBROUTINE solve_analytical_carbon_column

  !> Pivoted dense linear solve on the scaled Newton system; singular systems return failure.
  SUBROUTINE solve_column_linear_system(matrix, rhs, solution, ok)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    REAL(lp), INTENT(inout) :: matrix(:,:), rhs(:) !< Scaled Jacobian/system and right hand side [-]
    REAL(lp), INTENT(out) :: solution(:) !< Newton or regularized least-squares direction [-]
    LOGICAL, INTENT(out) :: ok !< Finite nonsingular solve [-]
    REAL(lp) :: row(SIZE(rhs)), value, factor, row_scale !< Elimination scratch [-]
    INTEGER :: n, i, j, pivot
    ok = .FALSE.
    solution = 0._lp
    n = SIZE(rhs)
    ! Row scaling improves pivot detection when different equations have different units/scales.
    DO i = 1,n
      row_scale = MAXVAL(ABS(matrix(i,:)))
      IF (row_scale <= TINY(1._lp)) RETURN
      matrix(i,:) = matrix(i,:)/row_scale
      rhs(i) = rhs(i)/row_scale
    END DO
    DO j = 1,n
      pivot = j-1+MAXLOC(ABS(matrix(j:n,j)),DIM=1)
      IF (ABS(matrix(pivot,j)) < 1.e-14_lp) RETURN
      row = matrix(j,:); matrix(j,:) = matrix(pivot,:); matrix(pivot,:) = row
      value = rhs(j); rhs(j) = rhs(pivot); rhs(pivot) = value
      DO i = j+1,n
        factor = matrix(i,j)/matrix(j,j)
        matrix(i,j:n) = matrix(i,j:n)-factor*matrix(j,j:n)
        rhs(i) = rhs(i)-factor*rhs(j)
      END DO
    END DO
    DO i = n,1,-1
      solution(i) = (rhs(i)-SUM(matrix(i,i+1:n)*solution(i+1:n)))/matrix(i,i)
    END DO
    ok = ALL(ieee_is_finite(solution))
  END SUBROUTINE solve_column_linear_system

  !> Regularized least-squares direction using QR, avoiding the squared conditioning of J-transpose-J.
  SUBROUTINE solve_column_regularized_step(jacobian, residual, damping, direction, ok)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    REAL(lp), INTENT(in) :: jacobian(:,:), residual(:), damping !< Scaled J, residual and damping [-]
    REAL(lp), INTENT(out) :: direction(:) !< Regularized search direction [-]
    LOGICAL, INTENT(out) :: ok !< Finite QR solve [-]
    REAL(lp) :: augmented(2*SIZE(residual),SIZE(residual)), rhs(2*SIZE(residual))
    REAL(lp) :: reflector(2*SIZE(residual)), length, diagonal, projection
    INTEGER :: n, j, k
    n = SIZE(residual)
    ok = .FALSE.
    direction = 0._lp
    augmented = 0._lp
    augmented(1:n,:) = jacobian
    rhs = 0._lp
    rhs(1:n) = -residual
    DO j = 1,n
      augmented(n+j,j) = SQRT(damping)*MAX(NORM2(jacobian(:,j)),1.e-12_lp)
    END DO
    DO j = 1,n
      length = NORM2(augmented(j:2*n,j))
      IF (length <= TINY(1._lp)) RETURN
      diagonal = -SIGN(length,augmented(j,j))
      reflector(j:2*n) = augmented(j:2*n,j)
      reflector(j) = reflector(j)-diagonal
      length = NORM2(reflector(j:2*n))
      IF (length <= TINY(1._lp)) RETURN
      reflector(j:2*n) = reflector(j:2*n)/length
      DO k = j,n
        projection = 2._lp*DOT_PRODUCT(reflector(j:2*n),augmented(j:2*n,k))
        augmented(j:2*n,k) = augmented(j:2*n,k)-projection*reflector(j:2*n)
      END DO
      projection = 2._lp*DOT_PRODUCT(reflector(j:2*n),rhs(j:2*n))
      rhs(j:2*n) = rhs(j:2*n)-projection*reflector(j:2*n)
      augmented(j,j) = diagonal
      augmented(j+1:2*n,j) = 0._lp
    END DO
    DO j = n,1,-1
      direction(j) = (rhs(j)-DOT_PRODUCT(augmented(j,j+1:n),direction(j+1:n)))/augmented(j,j)
    END DO
    ok = ALL(ieee_is_finite(direction))
  END SUBROUTINE solve_column_regularized_step

  !> Continue from the local analytical seed to monthly kinetics, then to full bioturbation.
  !! Each successful intermediate is the next seed. Failed stages halve the increment;
  !! no intermediate is accepted as the final target. Both requested fractions must reach one.
  SUBROUTINE solve_analytical_carbon_profile(ns, thickness, root_fraction_mean, water_rate, lateral_return, &
    & coefficients, solution, include_bioturbation, max_iterations, monthly)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    USE mo_sb_constants, ONLY: min_enzyme_fraction
    INTEGER, INTENT(in) :: ns !< Complete active profile [-]
    REAL(lp), INTENT(in) :: thickness(ns), root_fraction_mean(ns), water_rate(ns,2), lateral_return
    !< Layer thickness [m], annual root fraction [-], water rates [s-1], return fraction [-]
    TYPE(t_sb_local_coefficients), INTENT(in) :: coefficients(ns) !< Fixed conditional kinetics and geometry
    TYPE(t_sb_column_solution), INTENT(out) :: solution !< Final candidate, seed and all stage traces
    LOGICAL, OPTIONAL, INTENT(in) :: include_bioturbation !< Default true; final lambda one or zero [-]
    INTEGER, OPTIONAL, INTENT(in) :: max_iterations !< Per-stage accepted-update budget; default 80 [-]
    TYPE(t_sb_monthly_coefficients), OPTIONAL, INTENT(in) :: monthly !< Required final monthly-flux problem
    INTEGER, PARAMETER :: max_stages = 64 !< Bounded numerical work, including retries [-]
    REAL(lp), PARAMETER :: minimum_increment = 1.e-4_lp !< Abort continuation below this step [-]
    TYPE(t_sb_local_coefficients) :: local_coefficients
    TYPE(t_sb_local_target) :: local_target
    TYPE(t_sb_column_solution) :: stage
    TYPE(t_sb_continuation_stage) :: trace(max_stages)
    REAL(lp) :: seed(8,ns), alpha(ns), original_seed(8,ns), original_alpha(ns)
    REAL(lp) :: inflow, fraction, requested_fraction, trial_fraction, increment
    REAL(lp) :: seasonal_fraction, trial_monthly, requested_monthly
    INTEGER :: seed_status(ns), layer, count, evaluations, iterations
    LOGICAL :: completed_monthly !< Reset the step when switching from monthly to bioturbation continuation

    solution%reason = 'Invalid analytical continuation inputs'
    IF (ns <= 0) RETURN
    IF (.NOT. ALL(ieee_is_finite(thickness)) .OR. ANY(thickness <= 0._lp)) RETURN
    IF (.NOT. ALL(ieee_is_finite(water_rate)) .OR. ANY(water_rate < 0._lp)) RETURN
    IF (.NOT. ieee_is_finite(lateral_return)) RETURN
    IF (lateral_return < 0._lp .OR. lateral_return > 1._lp) RETURN
    seed = 0._lp
    alpha = min_enzyme_fraction
    seed_status = -1
    inflow = 0._lp
    DO layer = 1,ns
      local_coefficients = coefficients(layer)
      local_coefficients%transport_dom_inflow = inflow/thickness(layer)
      local_coefficients%dom_export_rate = water_rate(layer,1)+(1._lp-lateral_return)*water_rate(layer,2)
      CALL solve_local_carbon_target(local_coefficients,local_target)
      seed_status(layer) = local_target%status
      ! Local status 3/4/5 may contain a complete algebraic candidate rejected by a local guard.
      ! It is only a numerical seed. Stage zero rechecks every equation; the final stage rechecks
      ! all existing physical guards. Invalid/incomplete local results do not fall back to live pools.
      IF (local_target%status <= 2 .AND. local_target%status /= 0) THEN
        WRITE(solution%reason,'(A,I0,2A)') 'Analytical seed failed at layer ',layer,': ',TRIM(local_target%reason)
        EXIT
      END IF
      IF (.NOT. ALL(ieee_is_finite(local_target%carbon))) EXIT
      IF (ANY(local_target%carbon < 0._lp)) EXIT
      IF (local_target%carbon(4) < coefficients(layer)%minimum_microbial_carbon) EXIT
      IF (.NOT. ieee_is_finite(local_target%alpha_poly)) EXIT
      IF (local_target%alpha_poly < min_enzyme_fraction .OR. &
        & local_target%alpha_poly > 1._lp-min_enzyme_fraction) EXIT
      seed(:,layer) = local_target%carbon
      alpha(layer) = local_target%alpha_poly
      inflow = thickness(layer)*water_rate(layer,1)*seed(6,layer)
    END DO
    IF (layer <= ns) THEN
      solution%seed_carbon = seed
      solution%seed_alpha = alpha
      solution%seed_status = seed_status
      RETURN
    END IF
    original_seed = seed
    original_alpha = alpha
    requested_fraction = 1._lp
    IF (PRESENT(include_bioturbation)) THEN
      IF (.NOT. include_bioturbation) requested_fraction = 0._lp
    END IF
    fraction = 0._lp
    trial_fraction = 0._lp
    seasonal_fraction = 0._lp
    trial_monthly = 0._lp
    requested_monthly = 0._lp
    IF (PRESENT(monthly)) requested_monthly = 1._lp
    increment = 0.1_lp
    evaluations = 0
    iterations = 0
    DO count = 1,max_stages
      CALL solve_analytical_carbon_column(ns,thickness,root_fraction_mean,water_rate,lateral_return, &
        & coefficients,seed,alpha,stage,max_iterations=max_iterations,bioturbation_fraction=trial_fraction, &
        & monthly=monthly,monthly_fraction=trial_monthly)
      evaluations = evaluations+stage%evaluations
      iterations = iterations+stage%iterations
      trace(count)%fraction = trial_fraction
      trace(count)%increment = trial_fraction-fraction
      IF (seasonal_fraction < requested_monthly) trace(count)%increment = trial_monthly-seasonal_fraction
      trace(count)%monthly_fraction = trial_monthly
      trace(count)%converged = stage%converged
      trace(count)%admissible = stage%admissible
      trace(count)%iterations = stage%iterations
      trace(count)%evaluations = stage%evaluations
      trace(count)%worst_layer = stage%worst_layer
      trace(count)%worst_equation = stage%worst_equation
      trace(count)%reason = stage%reason
      IF (ALLOCATED(stage%norm_history)) trace(count)%norm_history = stage%norm_history
      IF (ALLOCATED(stage%step_history)) trace(count)%step_history = stage%step_history
      solution = stage
      IF (stage%converged) THEN
        completed_monthly = seasonal_fraction < requested_monthly .AND. trial_monthly == requested_monthly
        fraction = trial_fraction
        seasonal_fraction = trial_monthly
        seed = stage%carbon
        alpha = stage%alpha_poly
        IF (fraction == requested_fraction .AND. seasonal_fraction == requested_monthly) EXIT
        IF (completed_monthly) THEN
          increment = 0.1_lp
        ELSE IF (stage%iterations <= 8 .AND. count > 1) THEN
          increment = MIN(0.25_lp,1.5_lp*increment)
        END IF
      ELSE
        IF (count == 1) THEN
          solution%reason = 'Zero-bioturbation seed validation failed: '//TRIM(stage%reason)
          EXIT
        END IF
        increment = 0.5_lp*increment
        IF (increment < minimum_increment) THEN
          solution%reason = 'Analytical continuation minimum step reached: '//TRIM(stage%reason)
          EXIT
        END IF
      END IF
      IF (seasonal_fraction < requested_monthly) THEN
        trial_monthly = MIN(requested_monthly,seasonal_fraction+increment)
        trial_fraction = fraction
      ELSE
        trial_monthly = seasonal_fraction
        trial_fraction = MIN(requested_fraction,fraction+increment)
      END IF
    END DO
    IF (count > max_stages) THEN
      count = max_stages
      solution%reason = 'Analytical continuation attempt limit; no jump'
    END IF
    ! A converged intermediate is never substituted for the requested transport operator.
    IF (fraction /= requested_fraction .OR. seasonal_fraction /= requested_monthly .OR. .NOT. stage%converged) THEN
      solution%converged = .FALSE.
      solution%admissible = .FALSE.
    END IF
    solution%continuation = trace(1:count)
    solution%seed_carbon = original_seed
    solution%seed_alpha = original_alpha
    solution%seed_status = seed_status
    solution%evaluations = evaluations
    solution%total_iterations = iterations
  END SUBROUTINE solve_analytical_carbon_profile

  !> Evaluate closed-boundary particulate diffusion and carbon reactions on one candidate profile.
  !! Interface flux uses unequal-layer centre spacing; divergence divides by each receiving
  !! layer thickness. Only polymeric litter, microbes and residue diffuse in the fixed-pore
  !! scheme. DOC follows prescribed water export; MAOC remains fixed to the mineral matrix.
  !! Density includes mineral solids and occupied organic pore volume, not unburied surface POC.
  SUBROUTINE evaluate_analytical_column_residual(nc, ns, nactive, thickness, root_fraction_mean, &
    & water_rate, lateral_return, coefficients, carbon, alpha_poly, result, include_bioturbation, &
    & bioturbation_fraction)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    USE mo_jsb_math_constants, ONLY: eps8
    USE mo_sb_constants, ONLY: min_enzyme_fraction, min_dom_uptake_concentration
    USE mo_q_sb_jsm_processes, ONLY: calc_bulk_density_correction_available_pore_volume, &
      & calc_microbial_mortality_fraction, calc_enzyme_allocation
    USE mo_q_sb_jsm_transport, ONLY: calc_bioturbation_rate
    INTEGER, INTENT(in) :: nc, ns !< Cells and layers [-]
    REAL(lp), INTENT(in) :: nactive(nc), thickness(nc,ns), root_fraction_mean(nc,ns) !< [-], [m], [-]
    REAL(lp), INTENT(in) :: water_rate(nc,ns,2), lateral_return !< Down/lateral effective [s-1], returned fraction [-]
    TYPE(t_sb_local_coefficients), INTENT(in) :: coefficients(nc,ns) !< Effective layer kinetics; local DOM terms ignored
    REAL(lp), INTENT(in) :: carbon(8,nc,ns), alpha_poly(nc,ns) !< Candidate pools [mol C/m3], allocation [-]
    TYPE(t_sb_column_residual), INTENT(out) :: result !< Fresh scratch diagnostics and residuals
    LOGICAL, OPTIONAL, INTENT(in) :: include_bioturbation !< Default true; never changes forward transport [-]
    REAL(lp), OPTIONAL, INTENT(in) :: bioturbation_fraction !< Scratch continuation multiplier in [0,1] [-]
    TYPE(t_sb_local_coefficients) :: k !< Current layer coefficients
    REAL(lp) :: pools(8), dt, mortality, uptake, enzyme, enzyme_floor, enzymes(2), depoly(2)
    REAL(lp) :: maoc_volume, poc_volume, occupied_poc, free_pore, chemical_free, physical_free
    REAL(lp) :: adsorption(2), desorption(2), net_sorption(2), scale, turnover(2), respiration, lateral
    REAL(lp) :: allocation, log_mortality, forward_mortality !< Allocation [-], logarithm of rate, rate [mol C/m3/s]
    INTEGER :: c, s, p, n !< Cell, layer, pool and active count [-]

    result%reason = 'Invalid column geometry, forcing, coefficients or candidate.'
    IF (MIN(nc,ns) <= 0) RETURN
    IF (.NOT. ALL(ieee_is_finite(nactive))) RETURN
    IF (ANY(nactive < 1._lp) .OR. ANY(nactive > REAL(ns,lp))) RETURN
    IF (ANY(nactive /= REAL(NINT(nactive),lp))) RETURN
    IF (.NOT. ieee_is_finite(lateral_return)) RETURN
    IF (lateral_return < 0._lp .OR. lateral_return > 1._lp) RETURN
    IF (PRESENT(bioturbation_fraction)) THEN
      IF (.NOT. ieee_is_finite(bioturbation_fraction)) RETURN
      IF (bioturbation_fraction < 0._lp .OR. bioturbation_fraction > 1._lp) RETURN
    END IF
    ALLOCATE(result%reaction(8,nc,ns), result%bioturbation(8,nc,ns), result%total(8,nc,ns), SOURCE=0._lp)
    ALLOCATE(result%dom_transport(nc,ns), result%dom_interface(nc,0:ns), SOURCE=0._lp)
    ALLOCATE(result%bio_interface(8,nc,0:ns), result%column_budget(5,nc), SOURCE=0._lp)
    ALLOCATE(result%density(nc,ns), result%diffusivity(nc,ns), result%poc_excess(nc,ns), SOURCE=0._lp)
    ALLOCATE(result%burial(nc,ns), SOURCE=1._lp)
    ALLOCATE(result%allocation_residual(nc,ns), SOURCE=0._lp)
    ALLOCATE(result%process_rate(6,nc,ns), result%throughput(8,nc,ns), SOURCE=0._lp)
    ALLOCATE(result%allocation_defined(nc,ns), result%geometry_admissible(nc,ns), SOURCE=.FALSE.)
    ALLOCATE(result%timestep_warning(nc,ns,3), SOURCE=.FALSE.)
    DO c = 1,nc
      n = NINT(nactive(c))
      DO s = 1,n
        k = coefficients(c,s)
        pools = carbon(:,c,s)
        dt = k%dtime
        IF (.NOT. ALL(ieee_is_finite([pools,alpha_poly(c,s),thickness(c,s),root_fraction_mean(c,s), &
          & water_rate(c,s,:),dt,k%litter_input_carbon,k%litter_turnover_rate,k%carbon_transfer_fraction, &
          & k%microbial_decay_rate,k%minimum_microbial_carbon,k%mortality_exponent,k%mortality_reference, &
          & k%microbial_cue,k%enzyme_frac_ap_layer,k%vmax_uptake_eff,k%km_uptake_eff,k%vmax_depoly_eff, &
          & k%km_depoly_eff,k%enzyme_frac_total_local,k%persistent_enzyme_carbon,k%adsorption_rate, &
          & k%desorption_rate,k%qmax_org_layer,k%volume_pore_total,k%volume_per_carbon_maoc, &
          & k%volume_per_carbon_poc]))) RETURN
        IF (MINVAL([pools,root_fraction_mean(c,s),water_rate(c,s,:),k%litter_input_carbon, &
          & k%litter_turnover_rate,k%carbon_transfer_fraction,k%microbial_decay_rate,k%minimum_microbial_carbon, &
          & k%microbial_cue,k%enzyme_frac_ap_layer,k%vmax_uptake_eff,k%vmax_depoly_eff, &
          & k%enzyme_frac_total_local,k%persistent_enzyme_carbon,k%adsorption_rate,k%desorption_rate, &
          & k%qmax_org_layer,k%volume_pore_total]) < 0._lp) RETURN
        IF (MINVAL([dt,thickness(c,s),k%mortality_reference,k%km_uptake_eff,k%km_depoly_eff, &
          & k%volume_per_carbon_maoc,k%volume_per_carbon_poc]) <= 0._lp) RETURN
        IF (k%mortality_exponent < 1._lp .OR. &
          & MAXVAL([k%carbon_transfer_fraction,k%microbial_cue,k%enzyme_frac_ap_layer, &
          & k%enzyme_frac_total_local,k%volume_pore_total]) > 1._lp) RETURN
        IF (k%volume_pore_total > 1._lp-min_fine_mineral_volume+1.e-10_lp) RETURN
        IF (alpha_poly(c,s) < min_enzyme_fraction .OR. alpha_poly(c,s) > 1._lp-min_enzyme_fraction) RETURN

        maoc_volume = k%volume_per_carbon_maoc*SUM(pools(7:8))
        poc_volume = k%volume_per_carbon_poc*SUM(pools([1,3,4,5]))
        free_pore = MAX(0._lp,k%volume_pore_total-maoc_volume)
        occupied_poc = MIN(poc_volume,free_pore)
        result%poc_excess(c,s) = MAX(0._lp,poc_volume-free_pore)
        result%geometry_admissible(c,s) = maoc_volume <= k%volume_pore_total+1.e-12_lp .AND. &
          & (s == 1 .OR. result%poc_excess(c,s) <= 1.e-12_lp)
        IF (s == 1 .AND. poc_volume > eps8) result%burial(c,s) = MIN(1._lp,free_pore/poc_volume)
        result%density(c,s) = calc_bulk_density_correction_available_pore_volume( &
          & k%volume_pore_total,maoc_volume,occupied_poc)

        ! Floor biomass remains capable of uptake and depolymerization; it is not a dormant closed layer.
        mortality = 0._lp
        IF (pools(4) > k%minimum_microbial_carbon .AND. k%microbial_decay_rate > 0._lp) THEN
          log_mortality = LOG(k%microbial_decay_rate)+LOG(pools(4)) + &
            & (k%mortality_exponent-1._lp)*(LOG(pools(4))-LOG(k%mortality_reference))
          IF (log_mortality >= LOG(HUGE(1._lp))) RETURN
          mortality = EXP(log_mortality)
        END IF
        forward_mortality = pools(4)*calc_microbial_mortality_fraction(pools(4),k%microbial_decay_rate, &
          & k%mortality_reference,k%mortality_exponent,dt,k%minimum_microbial_carbon)/dt
        result%timestep_warning(c,s,1) = ABS(mortality-forward_mortality) > &
          & 1.e-12_lp*MAX(mortality,forward_mortality,TINY(1._lp))
        uptake = 0._lp
        IF (pools(6) > min_dom_uptake_concentration) &
          & uptake = k%vmax_uptake_eff*pools(4)*pools(6)/(k%km_uptake_eff+pools(6))
        result%timestep_warning(c,s,2) = uptake*dt > pools(6)
        enzyme_floor = 0._lp
        IF (k%mortality_exponent == 1._lp) enzyme_floor = k%persistent_enzyme_carbon
        enzyme = MAX(k%enzyme_frac_total_local*(1._lp-k%enzyme_frac_ap_layer)*pools(4),enzyme_floor)
        enzymes = [alpha_poly(c,s)*enzyme,MAX((1._lp-alpha_poly(c,s))*enzyme,enzyme_floor)]
        depoly = k%vmax_depoly_eff*enzymes/(k%km_depoly_eff+enzymes)*pools([3,5])
        ! Match ICON forward SESAM: the eps8 substrate cutoff belongs to beta=1 only.
        IF (k%mortality_exponent == 1._lp) THEN
          WHERE (pools([3,5]) <= eps8) depoly = 0._lp
        END IF
        ! Forward SESAM updates carbon allocation only when both substrate removal rates are positive.
        ! An inactive branch is underdetermined, not proof of a unique allocation equilibrium.
        result%allocation_defined(c,s) = ALL(depoly > 0._lp)
        IF (result%allocation_defined(c,s)) THEN
          allocation = calc_enzyme_allocation(depoly(1),depoly(2),k%km_depoly_eff,enzyme)
          allocation = MIN(1._lp-min_enzyme_fraction,MAX(min_enzyme_fraction,allocation))
          result%allocation_residual(c,s) = alpha_poly(c,s)-allocation
        END IF

        ! Evaluate exchange at the candidate AD/AR, never solve/reassign these pools inside a residual.
        chemical_free = MAX(0._lp,k%qmax_org_layer-SUM(pools(7:8)))
        physical_free = MAX(0._lp,k%volume_pore_total/k%volume_per_carbon_maoc-SUM(pools(7:8)))
        adsorption = pools([6,5])*MIN(1._lp,k%adsorption_rate*chemical_free*dt)
        scale = 1._lp
        IF (SUM(adsorption) > 0._lp) scale = MIN(1._lp,chemical_free/SUM(adsorption), &
          & physical_free/SUM(adsorption))
        adsorption = scale*adsorption
        desorption = MIN(1._lp,k%desorption_rate*dt)*pools(7:8)
        net_sorption = (adsorption-desorption)/dt
        result%process_rate(:,c,s) = [mortality,uptake,depoly,net_sorption]
        turnover = k%litter_turnover_rate*pools(1:2)
        result%timestep_warning(c,s,3) = ANY(turnover*dt > pools(1:2)) .OR. &
          & depoly(1)*dt > pools(3) .OR. depoly(2)*dt+adsorption(2) > pools(5) .OR. &
          & SUM(water_rate(c,s,:)*[1._lp,1._lp-lateral_return])*dt > 1._lp
        result%reaction(1:2,c,s) = k%litter_input_carbon(1:2)-turnover
        result%reaction(3,c,s) = k%litter_input_carbon(3)+k%carbon_transfer_fraction(2)*turnover(2)-depoly(1)
        result%reaction(4,c,s) = k%microbial_cue*uptake-mortality
        result%reaction(5,c,s) = (1._lp-k%carbon_transfer_fraction(3))*mortality-depoly(2)-net_sorption(2)
        result%reaction(6,c,s) = k%carbon_transfer_fraction(1)*turnover(1)+SUM(depoly) + &
          & k%carbon_transfer_fraction(3)*mortality-uptake-net_sorption(1)
        result%reaction(7:8,c,s) = net_sorption
        result%dom_interface(c,s) = thickness(c,s)*water_rate(c,s,1)*pools(6)
        lateral = (1._lp-lateral_return)*water_rate(c,s,2)*pools(6)
        result%dom_transport(c,s) = (result%dom_interface(c,s-1)-result%dom_interface(c,s))/thickness(c,s)-lateral
        ! Gross source+sink, not net exchange: preserve the original acceptance allowances.
        ! Store it with each monthly evaluation so nonlinear/capped sorption is averaged correctly.
        result%throughput(1:2,c,s) = k%litter_input_carbon(1:2)+turnover
        result%throughput(3,c,s) = k%litter_input_carbon(3)+k%carbon_transfer_fraction(2)*turnover(2)+depoly(1)
        result%throughput(4,c,s) = k%microbial_cue*uptake+mortality
        result%throughput(5,c,s) = (1._lp-k%carbon_transfer_fraction(3))*mortality+depoly(2) + &
          & (adsorption(2)+desorption(2))/dt
        result%throughput(6,c,s) = k%carbon_transfer_fraction(1)*turnover(1)+SUM(depoly) + &
          & k%carbon_transfer_fraction(3)*mortality+uptake+(adsorption(1)+desorption(1))/dt + &
          & (result%dom_interface(c,s-1)+result%dom_interface(c,s))/thickness(c,s)+lateral
        result%throughput(7:8,c,s) = (adsorption+desorption)/dt
        respiration = SUM((1._lp-k%carbon_transfer_fraction(1:2))*turnover)+(1._lp-k%microbial_cue)*uptake
        result%column_budget(1,c) = result%column_budget(1,c)+thickness(c,s)*SUM(k%litter_input_carbon)
        result%column_budget(2,c) = result%column_budget(2,c)+thickness(c,s)*respiration
        result%column_budget(3,c) = result%column_budget(3,c)+thickness(c,s)*lateral
      END DO
      result%column_budget(3,c) = result%column_budget(3,c)+result%dom_interface(c,n)
    END DO
    result%diffusivity = calc_bioturbation_rate(nc,ns,nactive,thickness,result%density,root_fraction_mean)
    IF (PRESENT(bioturbation_fraction)) result%diffusivity = bioturbation_fraction*result%diffusivity
    ! Controlled experiment: omit only analytical particulate exchange, retaining density and DOM transport.
    IF (PRESENT(include_bioturbation)) THEN
      IF (.NOT. include_bioturbation) result%diffusivity = 0._lp
    END IF
    DO p = 3,5
      CALL evaluate_analytical_bioturbation(nc,ns,nactive,thickness,result%diffusivity,carbon(p,:,:), &
        & result%bioturbation(p,:,:),result%bio_interface(p,:,:))
    END DO
    result%total = result%reaction+result%bioturbation
    result%total(6,:,:) = result%total(6,:,:)+result%dom_transport
    DO c = 1,nc
      DO s = 1,NINT(nactive(c))
        DO p = 3,5
          result%throughput(p,c,s) = result%throughput(p,c,s) + &
            & (ABS(result%bio_interface(p,c,s-1))+ABS(result%bio_interface(p,c,s)))/thickness(c,s)
        END DO
        result%column_budget(4,c) = result%column_budget(4,c)+thickness(c,s)*SUM(result%total(:,c,s))
      END DO
      result%column_budget(5,c) = result%column_budget(4,c) - &
        & (result%column_budget(1,c)-result%column_budget(2,c)-result%column_budget(3,c))
    END DO
    IF (.NOT. ALL(ieee_is_finite(result%total)) .OR. &
      & .NOT. ALL(ieee_is_finite(result%column_budget)) .OR. &
      & .NOT. ALL(ieee_is_finite(result%density)) .OR. &
      & .NOT. ALL(ieee_is_finite(result%diffusivity)) .OR. &
      & .NOT. ALL(ieee_is_finite(result%process_rate)) .OR. &
      & .NOT. ALL(ieee_is_finite(result%throughput)) .OR. &
      & .NOT. ALL(ieee_is_finite(result%allocation_residual))) RETURN
    result%valid = .TRUE.
    result%reason = 'Conditional carbon residual evaluated; inspect geometry, allocation and timestep warnings.'
  END SUBROUTINE evaluate_analytical_column_residual

  !> Average monthly nonlinear fluxes at one candidate, with annual CUE/AP and root fractions.
  !! This is not evaluation at mean kinetics. Monthly biology experiments are not used.
  !! Closure averages allocation updates only over months where both substrates are active.
  SUBROUTINE evaluate_analytical_monthly_residual(ns, thickness, roots, water_rate, lateral_return, &
    & coefficients, carbon, alpha, monthly, monthly_fraction, bioturbation_fraction, result)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    INTEGER, INTENT(in) :: ns !< Active layers [-]
    REAL(lp), INTENT(in) :: thickness(ns), roots(ns), water_rate(ns,2), lateral_return
    !< Thickness [m], annual fractions [-], mean water rates [s-1], lateral return [-]
    TYPE(t_sb_local_coefficients), INTENT(in) :: coefficients(ns) !< Former annual-coefficient problem
    REAL(lp), INTENT(in) :: carbon(8,1,ns), alpha(1,ns) !< Frozen candidate [mol C/m3], [-]
    TYPE(t_sb_monthly_coefficients), INTENT(in) :: monthly !< Monthly coefficients and replay-day weights
    REAL(lp), INTENT(in) :: monthly_fraction, bioturbation_fraction !< Continuation parameters in [0,1] [-]
    TYPE(t_sb_column_residual), INTENT(out) :: result !< Weighted residual and per-month diagnostics
    TYPE(t_sb_column_residual) :: baseline, average, bin
    REAL(lp) :: weight, active_weight(ns), fraction
    INTEGER :: nbin, b, s

    result%reason = 'Invalid monthly analytical forcing or continuation fraction'
    IF (.NOT. ieee_is_finite(monthly_fraction)) RETURN
    IF (monthly_fraction < 0._lp .OR. monthly_fraction > 1._lp) RETURN
    IF (.NOT. ALLOCATED(monthly%weight) .OR. .NOT. ALLOCATED(monthly%layer)) RETURN
    nbin = SIZE(monthly%weight)
    IF (nbin < 1) RETURN
    IF (SIZE(monthly%layer,1) /= ns .OR. SIZE(monthly%layer,2) /= nbin) RETURN
    IF (.NOT. ALL(ieee_is_finite(monthly%weight)) .OR. ANY(monthly%weight < 0._lp)) RETURN
    IF (ABS(SUM(monthly%weight)-1._lp) > 1.e-12_lp) RETURN
    DO b = 1,nbin
      DO s = 1,ns
        ! Only environmental kinetics/litter vary by month. Biology, dt and geometry stay fixed.
        IF (monthly%layer(s,b)%microbial_cue /= coefficients(s)%microbial_cue .OR. &
          & monthly%layer(s,b)%enzyme_frac_ap_layer /= coefficients(s)%enzyme_frac_ap_layer .OR. &
          & monthly%layer(s,b)%dtime /= coefficients(s)%dtime .OR. &
          & monthly%layer(s,b)%volume_pore_total /= coefficients(s)%volume_pore_total .OR. &
          & monthly%layer(s,b)%qmax_org_layer /= coefficients(s)%qmax_org_layer) RETURN
      END DO
    END DO
    CALL evaluate_analytical_column_residual(1,ns,[REAL(ns,lp)],RESHAPE(thickness,[1,ns]), &
      & RESHAPE(roots,[1,ns]),RESHAPE(water_rate,[1,ns,2]),lateral_return, &
      & RESHAPE(coefficients,[1,ns]),carbon,alpha,baseline,bioturbation_fraction=bioturbation_fraction)
    IF (.NOT. baseline%valid) THEN
      result%reason = baseline%reason
      RETURN
    END IF
    IF (monthly_fraction == 0._lp) THEN
      result = baseline
      RETURN
    END IF
    average = baseline
    average%reaction = 0._lp
    average%process_rate = 0._lp
    average%throughput = 0._lp
    average%column_budget = 0._lp
    average%allocation_residual = 0._lp
    average%allocation_defined = .FALSE.
    average%timestep_warning = .FALSE.
    active_weight = 0._lp
    ALLOCATE(average%process_month(6,ns,nbin),average%reaction_month(8,ns,nbin),SOURCE=0._lp)
    ALLOCATE(average%allocation_month(ns,nbin),SOURCE=0._lp)
    ALLOCATE(average%allocation_active_month(ns,nbin),SOURCE=.FALSE.)
    ALLOCATE(average%caps_month(ns,nbin,3),SOURCE=.FALSE.)
    DO b = 1,nbin
      weight = monthly%weight(b)
      IF (weight == 0._lp) CYCLE
      CALL evaluate_analytical_column_residual(1,ns,[REAL(ns,lp)],RESHAPE(thickness,[1,ns]), &
        & RESHAPE(roots,[1,ns]),RESHAPE(water_rate,[1,ns,2]),lateral_return, &
        & RESHAPE(monthly%layer(:,b),[1,ns]),carbon,alpha,bin, &
        & bioturbation_fraction=bioturbation_fraction)
      IF (.NOT. bin%valid) THEN
        WRITE(result%reason,'(A,I0,2A)') 'Invalid monthly analytical residual, bin ',b,': ',TRIM(bin%reason)
        RETURN
      END IF
      average%process_month(:,:,b) = bin%process_rate(:,1,:)
      average%reaction_month(:,:,b) = bin%reaction(:,1,:)
      average%allocation_month(:,b) = bin%allocation_residual(1,:)
      average%allocation_active_month(:,b) = bin%allocation_defined(1,:)
      average%caps_month(:,b,:) = bin%timestep_warning(1,:,:)
      average%reaction = average%reaction+weight*bin%reaction
      average%process_rate = average%process_rate+weight*bin%process_rate
      average%throughput = average%throughput+weight*bin%throughput
      average%column_budget(1:3,:) = average%column_budget(1:3,:)+weight*bin%column_budget(1:3,:)
      average%timestep_warning = average%timestep_warning .OR. bin%timestep_warning
      DO s = 1,ns
        IF (.NOT. bin%allocation_defined(1,s)) CYCLE
        active_weight(s) = active_weight(s)+weight
        average%allocation_residual(1,s) = average%allocation_residual(1,s)+weight*bin%allocation_residual(1,s)
      END DO
    END DO
    DO s = 1,ns
      IF (active_weight(s) > 0._lp) THEN
        average%allocation_defined(1,s) = .TRUE.
        average%allocation_residual(1,s) = average%allocation_residual(1,s)/active_weight(s)
      END IF
    END DO
    ! Geometry, annual-root bioturbation and fixed-DOC water transport are identical in all bins.
    ! Apply them once, rather than introducing a fictitious monthly timestep or extra transport.
    average%total = average%reaction+baseline%bioturbation
    average%total(6,:,:) = average%total(6,:,:)+baseline%dom_transport
    DO s = 1,ns
      average%column_budget(4,1) = average%column_budget(4,1)+thickness(s)*SUM(average%total(:,1,s))
    END DO
    average%column_budget(5,:) = average%column_budget(4,:) - &
      & (average%column_budget(1,:)-average%column_budget(2,:)-average%column_budget(3,:))
    average%mean_process = baseline%process_rate
    average%mean_total = baseline%total
    result = average
    fraction = monthly_fraction
    IF (fraction < 1._lp) THEN
      result%reaction = (1._lp-fraction)*baseline%reaction+fraction*average%reaction
      result%total = (1._lp-fraction)*baseline%total+fraction*average%total
      result%process_rate = (1._lp-fraction)*baseline%process_rate+fraction*average%process_rate
      result%throughput = (1._lp-fraction)*baseline%throughput+fraction*average%throughput
      result%column_budget = (1._lp-fraction)*baseline%column_budget+fraction*average%column_budget
      result%timestep_warning = baseline%timestep_warning .OR. average%timestep_warning
      result%allocation_defined = baseline%allocation_defined .OR. average%allocation_defined
      result%allocation_residual = (1._lp-fraction)*baseline%allocation_residual+fraction*average%allocation_residual
    END IF
    result%valid = ALL(ieee_is_finite(result%total)) .AND. ALL(ieee_is_finite(result%throughput)) .AND. &
      & ALL(ieee_is_finite(result%column_budget)) .AND. ALL(ieee_is_finite(result%allocation_residual))
    result%reason = 'Monthly nonlinear carbon fluxes at frozen pools; annual CUE/AP and roots prescribed'
  END SUBROUTINE evaluate_analytical_monthly_residual

  !> Reuse the forward conservative operator with a one-second interval to obtain tendencies.
  PURE SUBROUTINE evaluate_analytical_bioturbation(nc, ns, nactive, thickness, diffusivity, carbon, &
    & tendency, interface_flux)
    USE mo_q_sb_jsm_transport, ONLY: calc_bioturbation_transport
    INTEGER, INTENT(in) :: nc, ns !< Cells and soil layers [-]
    REAL(lp), INTENT(in) :: nactive(nc) !< Active-layer counts, as in forward memory [-]
    REAL(lp), INTENT(in) :: thickness(nc,ns), diffusivity(nc,ns) !< Thickness [m], K [m2/s]
    REAL(lp), INTENT(in) :: carbon(nc,ns) !< One transported carbon pool [mol C/m3 total soil]
    REAL(lp), INTENT(out) :: tendency(nc,ns) !< Net bioturbation rate [mol C/m3/s]
    REAL(lp), INTENT(out) :: interface_flux(nc,0:ns) !< Downward interface rate [mol C/m2/s]

    tendency = 0._lp
    CALL calc_bioturbation_transport(nc,ns,1._lp,nactive,diffusivity,thickness,carbon, &
      & tendency,interface_flux)
  END SUBROUTINE evaluate_analytical_bioturbation

  !> Prepare annual root fractions using replay-day weights, without depth renormalization.
  !! Entry density is diagnostic only; the nonlinear solver recomputes density and K.
  SUBROUTINE prepare_analytical_bioturbation_inputs(nc, ns, nbin, nactive, roots, inputs, reason)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    USE mo_q_sb_jsm_transport, ONLY: calc_bioturbation_rate
    INTEGER, INTENT(in) :: nc, ns, nbin !< Cells, layers, replay bins [-]
    REAL(lp), INTENT(in) :: nactive(nc), roots(nc,nbin*ns) !< Active counts and monthly fractions [-]
    TYPE(t_sb_analytical_inputs), INTENT(inout) :: inputs !< Prepared forcing and entry geometry
    CHARACTER(len=*), INTENT(out) :: reason !< Preparation status [-]
    INTEGER :: c, s, b !< Cell, layer, bin [-]

    inputs%bioturbation_ready = .FALSE.
    reason = 'Bioturbation inputs rejected: invalid roots, weights or entry geometry.'
    IF (.NOT. inputs%ready .OR. MIN(nc,ns,nbin) <= 0) RETURN
    IF (.NOT. ALLOCATED(inputs%depth) .OR. .NOT. ALLOCATED(inputs%bulk_density) .OR. &
      & .NOT. ALLOCATED(inputs%month_weight)) RETURN
    IF (ANY(SHAPE(inputs%depth) /= [nc,ns]) .OR. &
      & ANY(SHAPE(inputs%bulk_density) /= [nc,ns]) .OR. SIZE(inputs%month_weight) /= nbin) RETURN
    IF (.NOT. ALL(ieee_is_finite(nactive)) .OR. .NOT. ALL(ieee_is_finite(roots))) RETURN
    IF (ANY(roots < 0._lp) .OR. ANY(nactive < 1._lp) .OR. ANY(nactive > REAL(ns,lp))) RETURN
    IF (ANY(nactive /= REAL(NINT(nactive),lp))) RETURN
    IF (.NOT. ALL(ieee_is_finite(inputs%month_weight))) RETURN
    IF (ANY(inputs%month_weight < 0._lp) .OR. ABS(SUM(inputs%month_weight)-1._lp) > 1.e-12_lp) RETURN
    DO c = 1,nc
      DO s = 1,NINT(nactive(c))
        IF (.NOT. ieee_is_finite(inputs%depth(c,s)) .OR. &
          & .NOT. ieee_is_finite(inputs%bulk_density(c,s))) RETURN
        IF (inputs%depth(c,s) <= 0._lp .OR. inputs%bulk_density(c,s) < 0._lp) RETURN
      END DO
    END DO
    inputs%nactive = nactive
    IF (ALLOCATED(inputs%root_fraction_mean)) DEALLOCATE(inputs%root_fraction_mean)
    IF (ALLOCATED(inputs%interface_distance)) DEALLOCATE(inputs%interface_distance)
    ALLOCATE(inputs%root_fraction_mean(nc,ns), SOURCE=0._lp)
    ALLOCATE(inputs%interface_distance(nc,0:ns), SOURCE=0._lp)
    DO c = 1,nc
      DO s = 1,NINT(nactive(c))
        DO b = 1,nbin
          inputs%root_fraction_mean(c,s) = inputs%root_fraction_mean(c,s) + &
            & inputs%month_weight(b)*roots(c,b+nbin*(s-1))
        END DO
      END DO
      DO s = 1,NINT(nactive(c))-1
        inputs%interface_distance(c,s) = inputs%depth(c,s)/2._lp + inputs%depth(c,s+1)/2._lp
      END DO
    END DO
    inputs%bioturbation_entry = calc_bioturbation_rate(nc,ns,nactive,inputs%depth, &
      & inputs%bulk_density,inputs%root_fraction_mean)
    IF (.NOT. ALL(ieee_is_finite(inputs%bioturbation_entry))) RETURN
    inputs%bioturbation_ready = .TRUE.
    reason = 'Annual-root bioturbation inputs prepared; no live pools modified.'
  END SUBROUTINE prepare_analytical_bioturbation_inputs

  !> Build monthly conditional kinetics and solve each complete active soil column.
  !! All layers share acceptance; non-carbon reconstruction is a separate atomic transaction.
  SUBROUTINE diagnose_full_profile_carbon_targets(inputs, carbon_before, lateral_return, targets, solutions, &
    & include_bioturbation)
    USE mo_lnd_bgcm_idx, ONLY: ixC, ix_soluable_litter, ix_woody_litter, ix_polymeric_litter
    USE mo_jsb_physical_constants, ONLY: molar_mass_C
    USE mo_sb_constants, ONLY: tau_soluable_litter, tau_woody_litter, fc_soluable2dom, &
      & fc_woody2polymeric, f_soluable_residue, k_mort_microbial, min_mic_biomass, &
      & microbial_mortality_exponent, microbial_mortality_reference, calc_persistent_enzyme_floor, &
      & vmax_mic_uptake, km_mic_uptake, vmax_depolymerisation_litter, vmax_depolymerisation_residue, &
      & km_depolymerisation, enzyme_frac_total, min_persistent_enzyme, min_enzyme_fraction, &
      & k_adsorpt_dom, k_adsorpt_det, k_desorpt_dom, k_desorpt_det, carbon_per_dryweight_SOM, &
      & rho_solid_maoc, rho_bulk_org
    TYPE(t_sb_analytical_inputs), INTENT(in) :: inputs !< Validated forcing and fixed geometry
    REAL(lp), INTENT(in) :: carbon_before(:,:,:), lateral_return !< S,W,L,B,R,D,AD,AR/cell/layer [mol C/m3], [-]
    TYPE(t_sb_local_target), ALLOCATABLE, INTENT(out) :: targets(:,:) !< Mapped full-column targets, including failures
    TYPE(t_sb_column_solution), ALLOCATABLE, INTENT(out) :: solutions(:) !< Detailed solver outcomes per cell
    LOGICAL, OPTIONAL, INTENT(in) :: include_bioturbation !< Default true; analytical residual only [-]
    TYPE(t_sb_local_coefficients), ALLOCATABLE :: coefficients(:) !< Conditional layer kinetics
    TYPE(t_sb_monthly_coefficients) :: monthly !< Annual biology with monthly nonlinear kinetics
    REAL(lp), ALLOCATABLE :: alpha(:) !< Historical allocation for local coefficient validation [-]
    REAL(lp) :: products(5) !< Monthly rt*rm for depoly/uptake/sorption/desorption/HSC [-]
    REAL(lp) :: floor_enzyme !< Persistent enzyme C, linear mortality only [mol C/m3]
    INTEGER :: c, s, n, nc, ns, b, nbin !< Cell/layer/bin dimensions [-]

    nc = SIZE(carbon_before,2)
    ns = SIZE(carbon_before,3)
    ALLOCATE(targets(nc,ns),solutions(nc))
    solutions%reason = 'Incomplete or invalid monthly analytical input arrays'
    IF (.NOT. inputs%ready .OR. .NOT. inputs%bioturbation_ready .OR. SIZE(carbon_before,1) /= 8) RETURN
    IF (.NOT. ALLOCATED(inputs%month_weight) .OR. .NOT. ALLOCATED(inputs%modifier_month) .OR. &
      & .NOT. ALLOCATED(inputs%litter_month)) RETURN
    nbin = SIZE(inputs%month_weight)
    IF (nbin < 1) RETURN
    IF (SIZE(inputs%modifier_month,1) /= nc .OR. SIZE(inputs%modifier_month,2) /= ns .OR. &
      & SIZE(inputs%modifier_month,3) /= 10 .OR. SIZE(inputs%modifier_month,4) /= nbin) RETURN
    IF (SIZE(inputs%litter_month,1) < MAX(ix_soluable_litter,ix_woody_litter,ix_polymeric_litter) .OR. &
      & SIZE(inputs%litter_month,2) < ixC .OR. SIZE(inputs%litter_month,3) /= nc .OR. &
      & SIZE(inputs%litter_month,4) /= ns .OR. SIZE(inputs%litter_month,5) /= nbin) RETURN
    floor_enzyme = calc_persistent_enzyme_floor(microbial_mortality_exponent,min_persistent_enzyme)
    DO c = 1,nc
      n = NINT(inputs%nactive(c))
      ALLOCATE(coefficients(n),alpha(n))
      alpha = MIN(1._lp-min_enzyme_fraction,MAX(min_enzyme_fraction,inputs%alpha_poly(c,1:n)))
      DO s = 1,n
        coefficients(s)%dtime = inputs%timestep
        coefficients(s)%litter_input_carbon = &
          & inputs%litter_mean([ix_soluable_litter,ix_woody_litter,ix_polymeric_litter],ixC,c,s)
        coefficients(s)%litter_turnover_rate = inputs%rate_product_mean(c,s,1)/ &
          & [tau_soluable_litter,tau_woody_litter]
        coefficients(s)%carbon_transfer_fraction = [fc_soluable2dom,fc_woody2polymeric,f_soluable_residue]
        coefficients(s)%microbial_decay_rate = k_mort_microbial
        coefficients(s)%minimum_microbial_carbon = min_mic_biomass
        coefficients(s)%mortality_exponent = microbial_mortality_exponent
        coefficients(s)%mortality_reference = microbial_mortality_reference
        coefficients(s)%microbial_cue = inputs%cue(c,s)
        coefficients(s)%alpha_poly = alpha(s)
        coefficients(s)%enzyme_frac_ap_layer = inputs%alpha_ap(c,s)
        ! Uptake converts the forward micromolar coefficient to molar units; products are mean(rt*rm).
        coefficients(s)%vmax_uptake_eff = 1.e-6_lp*vmax_mic_uptake*inputs%rate_product_mean(c,s,2)
        coefficients(s)%km_uptake_eff = km_mic_uptake*inputs%rate_product_mean(c,s,5)
        coefficients(s)%vmax_depoly_eff = [vmax_depolymerisation_litter,vmax_depolymerisation_residue]* &
          & inputs%rate_product_mean(c,s,1)
        coefficients(s)%km_depoly_eff = km_depolymerisation*inputs%rate_product_mean(c,s,5)
        coefficients(s)%enzyme_frac_total_local = enzyme_frac_total
        coefficients(s)%persistent_enzyme_carbon = floor_enzyme
        coefficients(s)%adsorption_rate = [k_adsorpt_dom,k_adsorpt_det]*inputs%rate_product_mean(c,s,3)
        coefficients(s)%desorption_rate = [k_desorpt_dom,k_desorpt_det]*inputs%rate_product_mean(c,s,4)
        coefficients(s)%qmax_org_layer = inputs%qmax_org(c,s)
        coefficients(s)%volume_pore_total = inputs%pore_volume(c,s)
        coefficients(s)%volume_per_carbon_maoc = molar_mass_C/1000._lp/carbon_per_dryweight_SOM/rho_solid_maoc
        coefficients(s)%volume_per_carbon_poc = molar_mass_C/1000._lp/carbon_per_dryweight_SOM/rho_bulk_org
        coefficients(s)%l_surface_layer = s == 1
      END DO
      monthly%weight = inputs%month_weight
      ALLOCATE(monthly%layer(n,nbin))
      DO b = 1,nbin
        DO s = 1,n
          ! Preserve annual CUE/AP and fixed geometry/mortality; only kinetic rates/litter vary.
          monthly%layer(s,b) = coefficients(s)
          products = inputs%modifier_month(c,s,1:5,b)*inputs%modifier_month(c,s,6:10,b)
          monthly%layer(s,b)%litter_input_carbon = &
            & inputs%litter_month([ix_soluable_litter,ix_woody_litter,ix_polymeric_litter],ixC,c,s,b)
          monthly%layer(s,b)%litter_turnover_rate = products(1)/[tau_soluable_litter,tau_woody_litter]
          monthly%layer(s,b)%vmax_uptake_eff = 1.e-6_lp*vmax_mic_uptake*products(2)
          monthly%layer(s,b)%km_uptake_eff = km_mic_uptake*products(5)
          monthly%layer(s,b)%vmax_depoly_eff = &
            & [vmax_depolymerisation_litter,vmax_depolymerisation_residue]*products(1)
          monthly%layer(s,b)%km_depoly_eff = km_depolymerisation*products(5)
          monthly%layer(s,b)%adsorption_rate = [k_adsorpt_dom,k_adsorpt_det]*products(3)
          monthly%layer(s,b)%desorption_rate = [k_desorpt_dom,k_desorpt_det]*products(4)
        END DO
      END DO
      CALL solve_analytical_carbon_profile(n,inputs%depth(c,1:n),inputs%root_fraction_mean(c,1:n), &
        & inputs%water_rate_mean(c,1:n,:),lateral_return,coefficients,solutions(c), &
        & include_bioturbation=include_bioturbation,monthly=monthly)
      DO s = 1,n
        targets(c,s)%status = 5
        targets(c,s)%reason = solutions(c)%reason
        ! Status zero is granted only to the whole converged/admissible column, never to isolated layers.
        IF (solutions(c)%converged .AND. solutions(c)%admissible) targets(c,s)%status = 0
        IF (.NOT. ALLOCATED(solutions(c)%carbon)) CYCLE
        targets(c,s)%carbon = solutions(c)%carbon(:,s)
        targets(c,s)%alpha_poly = solutions(c)%alpha_poly(s)
        targets(c,s)%alpha_residue = 1._lp-targets(c,s)%alpha_poly
        targets(c,s)%enzyme_depoly = MAX(enzyme_frac_total*(1._lp-inputs%alpha_ap(c,s))* &
          & targets(c,s)%carbon(4),floor_enzyme)
        targets(c,s)%enzyme_poly = targets(c,s)%alpha_poly*targets(c,s)%enzyme_depoly
        targets(c,s)%enzyme_residue = MAX(targets(c,s)%alpha_residue*targets(c,s)%enzyme_depoly,floor_enzyme)
        targets(c,s)%maoc_volume = coefficients(s)%volume_per_carbon_maoc*SUM(targets(c,s)%carbon(7:8))
        targets(c,s)%poc_volume = coefficients(s)%volume_per_carbon_poc*SUM(targets(c,s)%carbon([1,3,4,5]))
        IF (.NOT. solutions(c)%diagnostics%valid) CYCLE
        targets(c,s)%residual = solutions(c)%diagnostics%total(:,1,s)
        targets(c,s)%burial = solutions(c)%diagnostics%burial(1,s)
        targets(c,s)%poc_excess = solutions(c)%diagnostics%poc_excess(1,s)
        targets(c,s)%microbial_mortality = solutions(c)%diagnostics%process_rate(1,1,s)
        targets(c,s)%uptake_step = solutions(c)%diagnostics%process_rate(2,1,s)*inputs%timestep
        targets(c,s)%transport = [solutions(c)%diagnostics%dom_interface(1,s-1)/inputs%depth(c,s), &
          & solutions(c)%diagnostics%dom_interface(1,s)/inputs%depth(c,s), &
          & (1._lp-lateral_return)*inputs%water_rate_mean(c,s,2)*targets(c,s)%carbon(6)]
        targets(c,s)%l_minimum_microbial_seed = targets(c,s)%carbon(4) <= min_mic_biomass
      END DO
      DEALLOCATE(coefficients,alpha,monthly%layer)
    END DO
  END SUBROUTINE diagnose_full_profile_carbon_targets

  !> Conditional five-pool SSM carbon target, without nutrient or solid-transport equilibrium.
  !> S,W,L,F,R are concentrations per TOTAL soil volume. Fluxes are mol C m-3 s-1.
  !> kF*F = fLF*(IS+IL+fWD*IW)/(1-fFR*fRF); kR*R=fFR*kF*F.
  !> Woody splitting is supplied from the forward C:N-dependent rule at the retained composition.
  !> Zero removal makes a zero-input pool nonunique: reject, never manufacture an empty target.
  PURE SUBROUTINE solve_simple_carbon_target(input, rate, transfer, woody_soluble, dt, target)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    REAL(lp), INTENT(in) :: input(3) !< S,W,L input [mol C m-3 s-1]
    REAL(lp), INTENT(in) :: rate(5) !< S,W,L,F,R first-order loss rates [s-1]
    REAL(lp), INTENT(in) :: transfer(4) !< Woody-to-litter, litter-to-fast, fast-to-slow, slow-to-fast [-]
    REAL(lp), INTENT(in) :: woody_soluble, dt !< Woody soluble fraction [-], timestep [s]
    TYPE(t_sb_local_target), INTENT(out) :: target !< S,W,L,F,R target, unused JSM entries zero
    REAL(lp) :: supply(5), loss(5), tolerance(5), respiration !< Carbon rates [mol C m-3 s-1]
    target = t_sb_local_target()
    target%reason = 'Rejected: invalid SSM forcing, rates or transfer fractions'
    IF (.NOT. ALL(ieee_is_finite(input)) .OR. .NOT. ALL(ieee_is_finite(rate)) .OR. &
      & .NOT. ALL(ieee_is_finite(transfer)) .OR. .NOT. ieee_is_finite(woody_soluble)) RETURN
    IF (ANY(input < 0._lp) .OR. ANY(rate < 0._lp) .OR. ANY(transfer < 0._lp) .OR. &
      & ANY(transfer > 1._lp) .OR. woody_soluble < 0._lp .OR. woody_soluble > 1._lp .OR. dt <= 0._lp) RETURN
    target%reason = 'Rejected: SSM zero-removal or non-dissipative recycling has no unique target'
    IF (ANY(rate == 0._lp) .OR. transfer(3)*transfer(4) >= 1._lp) RETURN
    target%reason = 'Rejected: SSM target activates the forward timestep donor limit'
    IF (ANY(rate*dt > 1._lp)) RETURN
    supply(1) = input(1) + woody_soluble*transfer(1)*input(2)
    supply(2) = input(2)
    supply(3) = input(3) + (1._lp-woody_soluble)*transfer(1)*input(2)
    supply(4) = transfer(2)*(supply(1)+supply(3))/(1._lp-transfer(3)*transfer(4))
    supply(5) = transfer(3)*supply(4)
    WHERE (rate < 1._lp)
      tolerance = HUGE(1._lp)*rate
    ELSEWHERE
      tolerance = HUGE(1._lp)
    END WHERE
    IF (ANY(supply > tolerance)) RETURN
    target%carbon(1:5) = supply/rate
    loss = rate*target%carbon(1:5)
    target%residual(1:5) = [input(1)+woody_soluble*transfer(1)*loss(2)-loss(1), &
      & input(2)-loss(2), input(3)+(1._lp-woody_soluble)*transfer(1)*loss(2)-loss(3), &
      & transfer(2)*(loss(1)+loss(3))+transfer(4)*loss(5)-loss(4), transfer(3)*loss(4)-loss(5)]
    tolerance = 1.e-18_lp+1.e-9_lp*MAX(supply,loss)
    target%status = 5
    target%reason = 'Rejected: SSM local or column carbon residual'
    IF (ANY(ABS(target%residual(1:5)) > tolerance)) RETURN
    respiration = (1._lp-transfer(1))*loss(2)+(1._lp-transfer(2))*(loss(1)+loss(3)) + &
      & (1._lp-transfer(3))*loss(4)+(1._lp-transfer(4))*loss(5)
    IF (ABS(SUM(input)-respiration) > 1.e-16_lp+1.e-9_lp*MAX(SUM(input),respiration)) RETURN
    target%status = 0
    target%reason = 'Accepted: conditional SSM carbon target; nutrient limitation and solid transport not solved'
  END SUBROUTINE solve_simple_carbon_target


  !> Same integer day-to-month mapping in forward replay and analytical input preparation.
  PURE INTEGER FUNCTION sb_replay_month(day, n_days, n_months) RESULT(month)
    INTEGER, INTENT(in) :: day, n_days, n_months !< One-based day and bin counts [-]
    month = MIN(n_months, (day - 1) * n_months / n_days + 1)
  END FUNCTION sb_replay_month

  !> Prepare coefficients and litter on scratch arrays only; no target solve or live-state update.
  SUBROUTINE build_sb_analytical_inputs(nc, ns, npool, nelem, nsource, nbin, nday, dt, sla, growthform, scheme, &
      & nactive, depth, field_cap, forcing_time, litter, roots, temperature, water, ice, water_potential, &
      & down_daily, lateral_daily, cue_sum, poly_sum, ap_sum, inputs, reason, partitioned_litter_rate)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    USE mo_q_sb_litter_processes, ONLY: calc_litter_partitioning
    USE mo_q_sb_rate_modifier, ONLY: calc_jsm_carbon_rate_modifiers, calc_ssm_carbon_rate_modifiers
    INTEGER, INTENT(in) :: nc, ns, npool, nelem, nsource, nbin, nday !< Array dimensions [-]
    CHARACTER(len=*), INTENT(in) :: scheme !< JSM or simple_1d partitioning and rate equations [-]
    REAL(lp), INTENT(in) :: water_potential(nc,nbin*ns) !< Monthly water potential [forward HYDRO unit]
    REAL(lp), INTENT(in) :: dt !< Replay timestep [s]
    REAL(lp), INTENT(in) :: sla !< Leaf specific area, as used by forward litter partitioning [m2 mol C-1]
    INTEGER, INTENT(in) :: growthform !< Vegetation growth form identifier [-]
    REAL(lp), INTENT(in) :: nactive(nc), depth(nc,ns), field_cap(nc,ns) !< Layers [-], thickness [m], FC [m3 m-3]
    REAL(lp), INTENT(in) :: forcing_time(nc,nbin) !< Coupled-window sampled time by bin [s]
    REAL(lp), INTENT(in) :: litter(:,:) !< Flattened bin/source/element litter [mol m-2 s-1]
    REAL(lp), INTENT(in) :: roots(nc,nbin*ns) !< Monthly root fractions [-]
    REAL(lp), INTENT(in) :: temperature(nc,nbin*ns) !< Monthly soil temperature [K]
    REAL(lp), INTENT(in) :: water(nc,nbin*ns), ice(nc,nbin*ns) !< Monthly liquid and ice depths [m]
    REAL(lp), INTENT(in) :: down_daily(nc,nday*ns), lateral_daily(nc,nday*ns) !< Daily water rates [s-1]
    REAL(lp), INTENT(in) :: cue_sum(nc,ns), poly_sum(nc,ns), ap_sum(nc,ns) !< Coupled biological integrals [s]
    TYPE(t_sb_analytical_inputs), INTENT(out) :: inputs !< Prepared read-only inputs; ready false on rejection
    CHARACTER(len=*), INTENT(out) :: reason !< Preparation outcome [-]
    REAL(lp), OPTIONAL, INTENT(in) :: partitioned_litter_rate(:,:,:,:,:)
      !< Prepartitioned pool/element/cell/layer/month input [mol m-3 s-1]
    REAL(lp) :: veg(nsource,nelem,nc), seed(nelem,nc) !< Scratch litter per replay step [mol m-2]
    REAL(lp) :: formation(npool,nelem,nc,ns), root(nc,ns) !< Scratch increments [mol m-3], roots [-]
    REAL(lp) :: coverage(nc), reference(nelem,nc), reconstructed, tolerance !< Time [s], rates [mol m-2 s-1]
    REAL(lp) :: modifiers(nc,ns,10), ssm_rt, ssm_rm !< Monthly scratch rate multipliers [-]
    INTEGER :: b, d, s, c, e, source, entry, k !< Bin, day, layer, cell, element, source and field indices [-]

    reason = 'Input preparation rejected: invalid dimensions, timestep or incomplete bookkeeping.'
    IF (MIN(nc,ns,npool,nelem,nsource,nbin,nday) <= 0 .OR. dt <= 0.0_lp) RETURN
    IF (SIZE(litter,1) /= nc .OR. SIZE(litter,2) /= nbin*(nsource+1)*nelem) RETURN
    IF (.NOT. ALL(ieee_is_finite(forcing_time)) .OR. ANY(forcing_time <= 0.0_lp)) RETURN
    coverage = SUM(forcing_time, DIM=2)
    IF (.NOT. ALL(ieee_is_finite(cue_sum)) .OR. .NOT. ALL(ieee_is_finite(poly_sum)) .OR. &
        & .NOT. ALL(ieee_is_finite(ap_sum))) RETURN
    IF (.NOT. ALL(ieee_is_finite(nactive)) .OR. .NOT. ALL(ieee_is_finite(depth)) .OR. &
        & .NOT. ALL(ieee_is_finite(field_cap))) RETURN
    DO c = 1,nc
      IF (nactive(c) < 1.0_lp .OR. nactive(c) > REAL(ns,lp)) RETURN
      DO s = 1,NINT(nactive(c))
        IF (MIN(cue_sum(c,s),poly_sum(c,s),ap_sum(c,s)) < 0.0_lp) RETURN
        IF (MAX(cue_sum(c,s),poly_sum(c,s),ap_sum(c,s)) > coverage(c)*(1.0_lp+1.e-12_lp)) RETURN
        IF (depth(c,s) <= 0.0_lp .OR. field_cap(c,s) <= 0.0_lp) RETURN
      END DO
    END DO
    IF (.NOT. ALL(ieee_is_finite(water_potential))) RETURN
    IF (.NOT. ALL(ieee_is_finite(litter)) .OR. ANY(litter < 0.0_lp)) RETURN
    IF (.NOT. ALL(ieee_is_finite(roots)) .OR. ANY(roots < 0.0_lp)) RETURN
    IF (.NOT. ALL(ieee_is_finite(temperature)) .OR. .NOT. ALL(ieee_is_finite(water)) .OR. &
        & .NOT. ALL(ieee_is_finite(ice))) RETURN
    IF (.NOT. ALL(ieee_is_finite(down_daily)) .OR. .NOT. ALL(ieee_is_finite(lateral_daily))) RETURN
    IF (ANY(down_daily < 0.0_lp) .OR. ANY(lateral_daily < 0.0_lp)) RETURN
    IF (PRESENT(partitioned_litter_rate)) THEN
      reason = 'Input preparation rejected: invalid prepartitioned litter.'
      IF (ANY(SHAPE(partitioned_litter_rate) /= [npool,nelem,nc,ns,nbin])) RETURN
      IF (.NOT. ALL(ieee_is_finite(partitioned_litter_rate))) RETURN
      IF (ANY(partitioned_litter_rate < 0._lp)) RETURN
      DO k = 1,npool
        IF (k == ix_soluable_litter .OR. k == ix_polymeric_litter .OR. k == ix_woody_litter) CYCLE
        IF (ANY(partitioned_litter_rate(k,:,:,:,:) /= 0._lp)) RETURN
      END DO
      DO c = 1,nc
        IF (ANY(partitioned_litter_rate(:,:,c,NINT(nactive(c))+1:ns,:) /= 0._lp)) RETURN
      END DO
    END IF

    inputs%timestep = dt
    ALLOCATE(inputs%month_weight(nbin), SOURCE=0.0_lp)
    DO d = 1,nday
      b = sb_replay_month(d, nday, nbin)
      inputs%month_weight(b) = inputs%month_weight(b) + 1.0_lp / REAL(nday,lp)
    END DO
    ALLOCATE(inputs%litter_mean(npool,nelem,nc,ns), SOURCE=0.0_lp)
    IF (scheme == 'jsm') THEN
      ! Scratch arrays for the current column only; no new forcing cache or restart fields.
      ALLOCATE(inputs%litter_month(npool,nelem,nc,ns,nbin), SOURCE=0._lp)
      ALLOCATE(inputs%modifier_month(nc,ns,10,nbin), SOURCE=0._lp)
    END IF
    ALLOCATE(inputs%rate_product_mean(nc,ns,6), SOURCE=0.0_lp)
    ALLOCATE(inputs%cue(nc,ns), inputs%alpha_poly(nc,ns), inputs%alpha_ap(nc,ns))
    ALLOCATE(inputs%water_rate_mean(nc,ns,2), SOURCE=0.0_lp)
    ALLOCATE(inputs%litter_budget_residual(nelem,nc), SOURCE=0.0_lp)
    inputs%cue = 0.0_lp
    inputs%alpha_poly = 0.0_lp
    inputs%alpha_ap = 0.0_lp
    DO c = 1,nc
      DO s = 1,NINT(nactive(c))
        inputs%cue(c,s) = cue_sum(c,s) / coverage(c)
        inputs%alpha_poly(c,s) = poly_sum(c,s) / coverage(c)
        inputs%alpha_ap(c,s) = ap_sum(c,s) / coverage(c)
        entry = (s-1)*nday
        inputs%water_rate_mean(c,s,1) = SUM(REAL(down_daily(c,entry+1:entry+nday),lp)) / REAL(nday,lp)
        inputs%water_rate_mean(c,s,2) = SUM(REAL(lateral_daily(c,entry+1:entry+nday),lp)) / REAL(nday,lp)
      END DO
    END DO
    reference = 0.0_lp
    DO b = 1,nbin
      modifiers = 0._lp
      DO e = 1,nelem
        DO source = 1,nsource
          entry = b + nbin*((source-1)+(nsource+1)*(e-1))
          veg(source,e,:) = litter(:,entry) * dt
        END DO
        entry = b + nbin*(nsource+(nsource+1)*(e-1))
        seed(e,:) = litter(:,entry) * dt
        reference(e,:) = reference(e,:) + &
          & inputs%month_weight(b)*(SUM(veg(:,e,:),DIM=1)+seed(e,:))/dt
      END DO
      DO s = 1,ns
        entry = b+nbin*(s-1)
        root(:,s) = roots(:,entry)
        DO c = 1,nc
          IF (s > NINT(nactive(c))) CYCLE
          IF (temperature(c,entry) <= 0.0_lp .OR. water(c,entry) < 0.0_lp .OR. ice(c,entry) < 0.0_lp) RETURN
          IF (scheme == 'jsm') THEN
            CALL calc_jsm_carbon_rate_modifiers(temperature(c,entry), water(c,entry), ice(c,entry), &
              & field_cap(c,s)*depth(c,s), depth(c,s), &
              & modifiers(c,s,1), modifiers(c,s,2), &
              & modifiers(c,s,3), modifiers(c,s,4), &
              & modifiers(c,s,5), modifiers(c,s,6), &
              & modifiers(c,s,7), modifiers(c,s,8), &
              & modifiers(c,s,9), modifiers(c,s,10))
          ELSE
            CALL calc_ssm_carbon_rate_modifiers(temperature(c,entry), water_potential(c,entry), ssm_rt, ssm_rm)
            inputs%rate_product_mean(c,s,6) = inputs%rate_product_mean(c,s,6) + &
              & inputs%month_weight(b)*ssm_rt*ssm_rm
          END IF
        END DO
      END DO
      formation = 0.0_lp
      ! Preserve tissue CNP, seed input, coarse-root depth distribution and the forward eps8 branches.
      IF (PRESENT(partitioned_litter_rate)) THEN
        ! Same per-PFT partition as replay, not a new partition of mixed tissues.
        formation = partitioned_litter_rate(:,:,:,:,b)*dt
      ELSE
        CALL calc_litter_partitioning(nc,ns,nactive,sla,growthform,scheme,depth,root,veg,seed,formation)
      END IF
      inputs%litter_mean = inputs%litter_mean + inputs%month_weight(b)*formation/dt
      IF (scheme == 'jsm') THEN
        inputs%litter_month(:,:,:,:,b) = formation/dt
        inputs%modifier_month(:,:,:,b) = modifiers
      END IF
      DO k = 1,5
        inputs%rate_product_mean(:,:,k) = inputs%rate_product_mean(:,:,k) + inputs%month_weight(b) * &
          & modifiers(:,:,k)*modifiers(:,:,k+5)
      END DO
    END DO
    IF (.NOT. ALL(ieee_is_finite(inputs%rate_product_mean)) .OR. &
        & .NOT. ALL(ieee_is_finite(inputs%litter_mean))) RETURN
    DO c = 1,nc
      DO e = 1,nelem
        reconstructed = 0.0_lp
        DO s = 1,NINT(nactive(c))
          reconstructed = reconstructed + SUM(inputs%litter_mean(:,e,c,s))*depth(c,s)
        END DO
        inputs%litter_budget_residual(e,c) = reconstructed - reference(e,c)
        tolerance = 1.e-10_lp*ABS(reference(e,c)) + 1.e-18_lp
        IF (ABS(inputs%litter_budget_residual(e,c)) > tolerance) THEN
          reason = 'Input preparation rejected: replay litter partition does not conserve the source input.'
          RETURN
        END IF
      END DO
    END DO
    inputs%ready = .TRUE.
    reason = 'Inputs prepared; mean-state approximation only; no analytical target applied.'
  END SUBROUTINE build_sb_analytical_inputs

  ! ======================================================================================================= !
  !> Diagnose independent layer targets from prepared means, with mineral structure held fixed.
  !>
  !> This adapter supplies current model constants to the infrastructure-free scalar solver. Effective
  !> rates use mean(rt*rm), not mean(rt)*mean(rm). In particular vmax_mic_uptake produces micromoles
  !> in calc_microbial_growth, so its conversion to the molar carbon equations requires 1.e-6.
  !> CUE and allocation are prescribed coupled-window means. No N/P feasibility, transport equilibrium,
  !> seasonal fixed point is claimed by a successful result. Optional coupling includes DOM advection only.
  !> Equation references below refer to JSM_pore_volume_analytical_CNP_derivation.md, Sections 3-11.
  !> Fixed entry geometry is a local-solver assumption only: annual calibration during forward replay
  !> can subsequently change V and Q, so this target is not a jointly calibrated equilibrium.
  SUBROUTINE diagnose_local_carbon_targets(inputs, nactive, targets, lateral_return)
    USE mo_lnd_bgcm_idx, ONLY: ixC, ix_soluable_litter, ix_woody_litter, ix_polymeric_litter
    USE mo_jsb_physical_constants, ONLY: molar_mass_C
    USE mo_sb_constants, ONLY: tau_soluable_litter, tau_woody_litter, fc_soluable2dom, &
      & fc_woody2polymeric, f_soluable_residue, k_mort_microbial, min_mic_biomass, &
      & vmax_mic_uptake, km_mic_uptake, vmax_depolymerisation_litter, vmax_depolymerisation_residue, &
      & km_depolymerisation, enzyme_frac_total, min_persistent_enzyme, k_adsorpt_dom, k_adsorpt_det, &
      & microbial_mortality_exponent, microbial_mortality_reference, calc_persistent_enzyme_floor, &
      & k_desorpt_dom, k_desorpt_det, carbon_per_dryweight_SOM, rho_solid_maoc, rho_bulk_org
    TYPE(t_sb_analytical_inputs), INTENT(in) :: inputs !< Prepared forcing and fixed geometry [various]
    REAL(lp), INTENT(in) :: nactive(:) !< Active soil-layer count per cell [-]
    TYPE(t_sb_local_target), ALLOCATABLE, INTENT(out) :: targets(:,:) !< Diagnostic-only cell/layer results
    REAL(lp), OPTIONAL, INTENT(in) :: lateral_return !< Enable downward/dom_lateral_export_rate coupling; returned fraction [0,1]
    TYPE(t_sb_local_coefficients) :: coefficients !< Layer effective coefficients [units as declared in type]
    INTEGER :: ic, is !< Cell and soil-layer indices [-]
    !> Interface influx [mol C m-2 s-1], export coefficients [s-1]
    REAL(lp) :: dom_inflow_flux, dom_vertical_export_rate, dom_lateral_export_rate
    REAL(lp) :: column_carbon_input, column_carbon_loss !< Independent whole-column input and external loss [mol C m-2 s-1]
    REAL(lp) :: column_carbon_storage !< Integrated seed storage tendency, including floor-biomass growth [mol C m-2 s-1]

    IF (.NOT. inputs%ready) RETURN
    ALLOCATE(targets(SIZE(inputs%cue,1),SIZE(inputs%cue,2)))
    coefficients%dtime = inputs%timestep
    coefficients%carbon_transfer_fraction = [fc_soluable2dom, fc_woody2polymeric, f_soluable_residue]
    coefficients%microbial_decay_rate = k_mort_microbial
    coefficients%mortality_exponent = microbial_mortality_exponent
    coefficients%mortality_reference = microbial_mortality_reference
    coefficients%minimum_microbial_carbon = min_mic_biomass
    coefficients%enzyme_frac_total_local = enzyme_frac_total
    coefficients%persistent_enzyme_carbon = &
      & calc_persistent_enzyme_floor(microbial_mortality_exponent, min_persistent_enzyme)
    !> Eq. (1): gamma_A=M_C/(1000*f_C*rho_A), gamma_P=M_C/(1000*f_C*rho_P) [m3/mol C].
    !> Multiplication by a concentration per TOTAL soil volume gives an occupied volume fraction.
    !> These are dry-solid occupancies; neither dissolved DOM nor woody litter occupies these pores.
    coefficients%volume_per_carbon_maoc = molar_mass_C / 1000.0_lp / carbon_per_dryweight_SOM / rho_solid_maoc
    coefficients%volume_per_carbon_poc = molar_mass_C / 1000.0_lp / carbon_per_dryweight_SOM / rho_bulk_org
    DO ic = 1,SIZE(nactive)
      dom_inflow_flux = 0.0_lp
      column_carbon_input = 0.0_lp
      column_carbon_loss = 0.0_lp
      column_carbon_storage = 0.0_lp
      DO is = 1,NINT(nactive(ic))
        coefficients%transport_dom_inflow = 0.0_lp
        coefficients%dom_export_rate = 0.0_lp
        dom_vertical_export_rate = 0.0_lp
        dom_lateral_export_rate = 0.0_lp
        IF (PRESENT(lateral_return)) THEN
          ! Eq. (24): j_i=F_(i-1)/h_i, ell_i=q_i+(1-r)*lambda_i, F_i=q_i*h_i*D_i.
          ! Bookkeeping ALREADY divided native hydrological coefficients by 1000: these are effective s-1.
          ! Only dom_lateral_export_rate loss receives (1-r). Bottom F_n remains an external leaching loss.
          IF (lateral_return < 0.0_lp .OR. lateral_return > 1.0_lp) RETURN
          IF (inputs%depth(ic,is) <= 0.0_lp) EXIT
          dom_vertical_export_rate = inputs%water_rate_mean(ic,is,1)
          dom_lateral_export_rate = (1.0_lp-lateral_return)*inputs%water_rate_mean(ic,is,2)
          coefficients%transport_dom_inflow = dom_inflow_flux/inputs%depth(ic,is)
          coefficients%dom_export_rate = dom_vertical_export_rate+dom_lateral_export_rate
        END IF
        coefficients%litter_input_carbon = inputs%litter_mean([ix_soluable_litter,ix_woody_litter,ix_polymeric_litter],ixC,ic,is)
        !> Section 4.1: kS=<rt_d*rm_d>/tau_S and kW=<rt_d*rm_d>/tau_W [s-1].
        !> The mean-state approximation freezes coefficients, not the full seasonally forced trajectory:
        !> Eq. (36), <v*B*D/(K+D)> /= <v>*<B>*<D>/(<K>+<D>), even with correct modifier averaging.
        coefficients%litter_turnover_rate = inputs%rate_product_mean(ic,is,1) / [tau_soluable_litter,tau_woody_litter]
        coefficients%microbial_cue = inputs%cue(ic,is)
        coefficients%alpha_poly = inputs%alpha_poly(ic,is)
        coefficients%enzyme_frac_ap_layer = inputs%alpha_ap(ic,is)
        !> Eq. (5): U=vU*B*D/(KU+D); vU=1.e-6*vmax_code*<rt_U*rm_U>, KU=km_U*<rt_h*rm_h>.
        !> This is POTENTIAL molar uptake. Prescribed CUE does not guarantee that N/P availability
        !> permits this uptake in the forward model; a state-dependent uptake cap cannot be absorbed
        !> into CUE while retaining the same DOM equation. No such nutrient feasibility test occurs here.
        coefficients%vmax_uptake_eff = 1.e-6_lp * vmax_mic_uptake * inputs%rate_product_mean(ic,is,2)
        coefficients%km_uptake_eff = km_mic_uptake * inputs%rate_product_mean(ic,is,5)
        coefficients%vmax_depoly_eff = [vmax_depolymerisation_litter,vmax_depolymerisation_residue] * &
          & inputs%rate_product_mean(ic,is,1)
        coefficients%km_depoly_eff = km_depolymerisation * inputs%rate_product_mean(ic,is,5)
        coefficients%adsorption_rate = [k_adsorpt_dom,k_adsorpt_det] * inputs%rate_product_mean(ic,is,3)
        coefficients%desorption_rate = [k_desorpt_dom,k_desorpt_det] * inputs%rate_product_mean(ic,is,4)
        !> Eq. (2): Q=q_min*rho_m*(1-V) is supplied by the entry-state mineral geometry [mol C/m3].
        !> Hold Q and V fixed in this solve. Do not replace Q by min(Q,V/gamma_A) in sorption kinetics:
        !> chemical sites set potential sorption; physical space only limits the resulting increments.
        coefficients%qmax_org_layer = inputs%qmax_org(ic,is)
        coefficients%volume_pore_total = inputs%pore_volume(ic,is)
        coefficients%l_surface_layer = is == 1
        CALL solve_local_carbon_target(coefficients,targets(ic,is))
        targets(ic,is)%transport = [coefficients%transport_dom_inflow,dom_vertical_export_rate*targets(ic, &
      & is)%carbon(6),dom_lateral_export_rate*targets(ic,is)%carbon(6)]
        IF (PRESENT(lateral_return)) THEN
          ! No upward advection or bioturbation in this target: the coupling is triangular, not iterative.
          ! If an upstream root is inadmissible, no downstream source can be certified. Reject the column.
          IF (targets(ic,is)%status /= 0) EXIT
          dom_inflow_flux = dom_vertical_export_rate*inputs%depth(ic,is)*targets(ic,is)%carbon(6)
          column_carbon_input = column_carbon_input+SUM(coefficients%litter_input_carbon)*inputs%depth(ic,is)
          column_carbon_loss = column_carbon_loss + ((1.0_lp-coefficients%carbon_transfer_fraction(1))* &
      & coefficients%litter_input_carbon(1) + &
            & (1.0_lp-coefficients%carbon_transfer_fraction(2))*coefficients%litter_input_carbon(2) + &
      & (1.0_lp-coefficients%microbial_cue)*targets(ic,is)%uptake_step/coefficients%dtime + &
            & dom_lateral_export_rate*targets(ic,is)%carbon(6))*inputs%depth(ic,is)
          column_carbon_storage = column_carbon_storage + SUM(targets(ic,is)%residual)*inputs%depth(ic,is)
        END IF
      END DO
      IF (PRESENT(lateral_return)) THEN
        IF (ALL(targets(ic,1:NINT(nactive(ic)))%status == 0)) THEN
          ! Eq. (12): internal interface transfers cancel in thickness-integrated balances.
          ! Input = external loss + storage tendency. The minimum-biomass seed retains positive B growth;
          ! its residual is storage, NOT respiration or an external loss.
          column_carbon_loss = column_carbon_loss+dom_inflow_flux
          IF (ABS(column_carbon_input-column_carbon_loss-column_carbon_storage) > &
            & 1.e-16_lp+1.e-9_lp*MAX(column_carbon_input,column_carbon_loss,ABS(column_carbon_storage))) THEN
            targets(ic,1:NINT(nactive(ic)))%status = 5
            targets(ic,1:NINT(nactive(ic)))%reason = 'Rejected: independent column C budget does not close'
          END IF
        END IF
      END IF
    END DO
  END SUBROUTINE diagnose_local_carbon_targets

  !> Fixed-pore local C kernel, ported from QS d79f601 (linear mortality path).
  !> Pool concentrations are per TOTAL soil volume; prescribed rates are per second.
  !> With U=vU*B*D/(KU+D) and M=kB*B, e*U=M and
  !> J=(1-e)*U+ellD*D give D=kB*KU/(e*vU-kB), B=e*(J-ellD*D)/(kB*(1-e)).
  !> L=(IL+fW*IW)/aL, R=(1-fR)*M/aR. MAOC uses the finite-step capped exchange,
  !> not an unlimited Langmuir extrapolation. Pore volume and chemical capacity stay fixed.
  !> A biomass-floor seed closes donor balances but retains positive microbial growth:
  !> it is explicitly NOT full steady state. Failed, nonunique or inadmissible roots reject.
  !> No solid transport, nutrient equilibrium or density calibration is solved here.
  PURE SUBROUTINE solve_local_carbon_target(coefficients, target)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    USE mo_jsb_math_constants, ONLY: eps8
    USE mo_q_sb_jsm_processes, ONLY: calc_enzyme_allocation
    USE mo_sb_constants, ONLY: min_enzyme_fraction, min_dom_uptake_concentration
    TYPE(t_sb_local_coefficients), INTENT(in) :: coefficients !< Fixed coefficients [see type declaration]
    TYPE(t_sb_local_target), INTENT(out) :: target !< Candidate, validity and balance diagnostics
    REAL(lp) :: carbon_pools(8) !< Candidate concentrations in S,W,L,B,R,D,AD,AR order [mol C m-3]
    !> Net litter supply, microbial death and gross DOC uptake [mol C m-3 s-1]
    REAL(lp) :: carbon_source, microbial_decay, microbial_uptake
    !> Total depolymerizing, polymeric and residue enzyme C [mol C m-3]
    REAL(lp) :: enzyme_depoly, enzyme_poly, enzyme_residue
    REAL(lp) :: alpha_poly, alpha_residue !< Carbon-based enzyme allocation fractions [-]
    REAL(lp) :: depolymerisation_rate(2) !< Effective polymeric/residue first-order depolymerization rates [s-1]
    REAL(lp) :: substrate_input_rate(2), net_sorption_rate(2) !< Substrate supply and net adsorption [mol C m-3 s-1]
    REAL(lp) :: sorption_increment(2) !< Gross DOM/residue sorption increments [mol C m-3 timestep-1]
    REAL(lp) :: residual_scale(8) !< Per-pool numerical balance tolerance [mol C m-3 s-1]
    REAL(lp) :: volume_pore_available !< Pore fraction available after MAOC occupancy [-]
    LOGICAL :: l_maoc_solution !< Finite-timestep MAOC solve passed [-]
    LOGICAL :: l_microbial_solution !< Nonlinear unrestricted root is resolved and admissible [-]
    REAL(lp) :: dom_lower_bound, dom_upper_bound, dom_removal !< DOC bracket [mol C/m3], removal [mol C/m3/s]
    INTEGER :: iter_dom !< Bounded monotonic DOC root iteration [-]

    ! Named pool aliases follow the reference layer-state solver; vector storage is unchanged.
    ASSOCIATE( &
      & soluable_carbon => carbon_pools(1), &
      & woody_carbon => carbon_pools(2), &
      & polymeric_carbon => carbon_pools(3), &
      & microbial_carbon => carbon_pools(4), &
      & residue_carbon => carbon_pools(5), &
      & dom_carbon => carbon_pools(6), &
      & dom_assoc_carbon => carbon_pools(7), &
      & residue_assoc_carbon => carbon_pools(8), &
      & input_sol => coefficients%litter_input_carbon(1), &
      & input_wood => coefficients%litter_input_carbon(2), &
      & input_poly => coefficients%litter_input_carbon(3), &
      & d_litter => substrate_input_rate(1), &
      & d_residue => substrate_input_rate(2))

    target = t_sb_local_target()
    target%reason = 'Invalid local coefficients'
    IF (.NOT. ALL(ieee_is_finite([coefficients%dtime,coefficients%litter_input_carbon, &
      & coefficients%litter_turnover_rate,coefficients%carbon_transfer_fraction, &
      & coefficients%microbial_decay_rate,coefficients%minimum_microbial_carbon, &
      & coefficients%mortality_exponent,coefficients%mortality_reference, &
      & coefficients%microbial_cue,coefficients%alpha_poly,coefficients%enzyme_frac_ap_layer, &
      & coefficients%vmax_uptake_eff,coefficients%km_uptake_eff,coefficients%vmax_depoly_eff, &
      & coefficients%km_depoly_eff,coefficients%enzyme_frac_total_local,coefficients%persistent_enzyme_carbon, &
      & coefficients%adsorption_rate,coefficients%desorption_rate,coefficients%qmax_org_layer, &
      & coefficients%volume_pore_total,coefficients%volume_per_carbon_maoc, &
      & coefficients%volume_per_carbon_poc,coefficients%transport_dom_inflow,coefficients%dom_export_rate]))) RETURN
    IF (coefficients%dtime <= 0.0_lp .OR. coefficients%microbial_decay_rate <= 0.0_lp .OR. &
      & coefficients%microbial_cue <= 0.0_lp .OR. coefficients%microbial_cue >= 1.0_lp) RETURN
    IF (coefficients%mortality_exponent < 1._lp .OR. coefficients%mortality_reference <= 0._lp) RETURN
    IF (MINVAL([coefficients%litter_input_carbon,coefficients%litter_turnover_rate, &
      & coefficients%minimum_microbial_carbon,coefficients%vmax_uptake_eff,coefficients%vmax_depoly_eff, &
      & coefficients%enzyme_frac_total_local, &
      & coefficients%persistent_enzyme_carbon,coefficients%adsorption_rate,coefficients%desorption_rate, &
      & coefficients%qmax_org_layer,coefficients%volume_pore_total,coefficients%transport_dom_inflow, &
      & coefficients%dom_export_rate]) < 0.0_lp) RETURN
    IF (MINVAL([coefficients%carbon_transfer_fraction,coefficients%alpha_poly, &
      & coefficients%enzyme_frac_ap_layer]) < 0.0_lp .OR. MAXVAL([coefficients%carbon_transfer_fraction, &
      & coefficients%alpha_poly,coefficients%enzyme_frac_ap_layer]) > 1.0_lp) RETURN
    IF (coefficients%volume_pore_total > 1._lp-min_fine_mineral_volume+1.e-10_lp) THEN
      target%reason = 'Fixed geometry violates minimum fine-mineral solid volume'
      RETURN
    END IF
    IF (coefficients%volume_pore_total > 1.0_lp .OR. MIN(coefficients%volume_per_carbon_maoc, &
      & coefficients%volume_per_carbon_poc,coefficients%km_uptake_eff,coefficients%km_depoly_eff) <= 0.0_lp) RETURN
    target%status = 2
    target%reason = 'No unique positive-growth root: zero removal/source or insufficient uptake capacity'
    !> Eqs. (13),(15): kS,kW>0 and e*vU-kB>0 are required for the divisions on this branch.
    !> A zero decay coefficient with positive input has no finite stationary litter stock; with zero
    !> input its stock is undetermined. Reject both rather than silently prescribing a zero pool.
    IF (ANY(coefficients%litter_turnover_rate <= 0.0_lp)) RETURN
    carbon_source = coefficients%carbon_transfer_fraction(1)*input_sol + input_poly + &
      & coefficients%carbon_transfer_fraction(2)*input_wood + coefficients%transport_dom_inflow
    !> Eqs. (16)-(17): J=fS*IS+IL+fW*IW+jD; net growth supply is J-ellD*D>0.
    !> jD comes from the solved upstream TARGET, not a recorded historical material flux.
    carbon_pools = 0.0_lp
    soluable_carbon = input_sol / coefficients%litter_turnover_rate(1)
    woody_carbon = input_wood / coefficients%litter_turnover_rate(2)
    !> Eq. (17) follows from 0=J+kB*B-U-ellD*D and e*U=kB*B. It relies on stationary MAOC net_sorption_rate
    !> and no microbial/residue transport; it cannot be retained after simply clipping B to B_min.
    !> Eq. (15): cancel B>0 in e*vU*B*D/(KU+D)=kB*B to obtain D=kB*KU/(e*vU-kB).
    !> This is not a DOC prescription for a dormant layer or a layer where actual U is capped.
    target%l_minimum_microbial_seed = .TRUE.
    IF (coefficients%mortality_exponent > 1._lp) THEN
      ! q(B)=k_ref*(B/B_ref)**(beta-1); e*U=M and J=(1-e)*M/e+ell*D.
      ! A failed nonlinear solve is NOT a reason to prescribe a dormant layer. Only a resolved
      ! root at/below B_min (or J=0) enters the explicit zero-mortality quasi-steady branch.
      CALL solve_density_dependent_microbial_target(carbon_source, coefficients%dom_export_rate, &
        & coefficients%microbial_cue, coefficients%vmax_uptake_eff, coefficients%km_uptake_eff, &
        & coefficients%microbial_decay_rate, coefficients%mortality_reference, coefficients%mortality_exponent, &
        & microbial_carbon, dom_carbon, microbial_decay, l_microbial_solution)
      target%reason = 'Nonlinear microbial/DOC root absent or unresolved; no analytical jump'
      IF (.NOT. l_microbial_solution) RETURN
      target%l_minimum_microbial_seed = microbial_carbon <= coefficients%minimum_microbial_carbon
    ELSE IF (coefficients%microbial_cue*coefficients%vmax_uptake_eff > coefficients%microbial_decay_rate) THEN
      dom_carbon = coefficients%microbial_decay_rate*coefficients%km_uptake_eff/(coefficients%microbial_cue* &
        & coefficients%vmax_uptake_eff-coefficients%microbial_decay_rate)
      microbial_carbon = coefficients%microbial_cue*(carbon_source-coefficients%dom_export_rate*dom_carbon)/ &
        & (coefficients%microbial_decay_rate*(1.0_lp-coefficients%microbial_cue))
      target%l_minimum_microbial_seed = microbial_carbon <= coefficients%minimum_microbial_carbon
    END IF
    IF (coefficients%mortality_exponent == 1._lp) &
      & microbial_decay = coefficients%microbial_decay_rate*microbial_carbon
    IF (target%l_minimum_microbial_seed) THEN
      microbial_carbon = coefficients%minimum_microbial_carbon
      microbial_decay = 0.0_lp
      !> At the floor forward mortality is zero, but uptake remains active. Solve the donor balance,
      !> not the active-growth D=kB*KU/(e*vU-kB). Removal is monotone in D for fixed B_min.
      !> If ellD>0, [0,J/ellD] brackets the root. With ellD=0 require J<vU*B_min;
      !> then D=J*KU/(vU*B_min-J). Zero supply/removal together is history dependent and rejected.
      target%reason = 'Minimum microbial seed has no finite unique DOC balance'
      IF (coefficients%dom_export_rate > 0.0_lp) THEN
        dom_lower_bound = 0.0_lp
        dom_upper_bound = carbon_source/coefficients%dom_export_rate
        DO iter_dom = 1,160
          dom_carbon = dom_lower_bound+0.5_lp*(dom_upper_bound-dom_lower_bound)
          dom_removal = coefficients%dom_export_rate*dom_carbon + &
            & coefficients%vmax_uptake_eff*microbial_carbon*dom_carbon/(coefficients%km_uptake_eff+dom_carbon)
          IF (dom_removal > carbon_source) THEN
            dom_upper_bound = dom_carbon
          ELSE
            dom_lower_bound = dom_carbon
          END IF
        END DO
      ELSE IF (carbon_source < coefficients%vmax_uptake_eff*microbial_carbon) THEN
        dom_carbon = carbon_source*coefficients%km_uptake_eff / &
          & (coefficients%vmax_uptake_eff*microbial_carbon-carbon_source)
      ELSE
        RETURN
      END IF
    END IF
    target%microbial_mortality = microbial_decay
    !> Eq. (8): E=max(fE*(1-AP)*B,E0), EL=alpha*E, ER=max((1-alpha)*E,E0).
    !> E0 is retained only for beta=1; nonlinear mortality uses E0=0 in both forward and target.
    !> Only ER has the additional independent floor. For small B with ER=E0, Eq. (18) gives R -> 0
    !> as B -> 0 while linear microbial_decay holds; dropping that floor can leave a spurious finite R.
    !> Eqs. (21)-(22): derive allocation from stationary carbon fluxes, not the historical allocation.
    !> d_litter/d_residue = alpha_poly*(Km+alpha_poly*E) /
    !>                     (alpha_residue*(Km+alpha_residue*E)).
    !> The shared forward function handles zero/equal demands; bound its result as in forward JSM.
    !> This is the model's carbon allocation rule, not a new optimum for the asymmetric residue floor.
    !> AP and CUE remain prescribed; nutrient-dependent allocation and its fixed point are not solved.
    d_litter = input_poly + coefficients%carbon_transfer_fraction(2) * input_wood
    d_residue = (1.0_lp - coefficients%carbon_transfer_fraction(3)) * microbial_decay
    enzyme_depoly = MAX(coefficients%enzyme_frac_total_local*(1.0_lp-coefficients%enzyme_frac_ap_layer)* &
      & microbial_carbon,coefficients%persistent_enzyme_carbon)
    alpha_poly = calc_enzyme_allocation(d_litter, d_residue, coefficients%km_depoly_eff, enzyme_depoly)
    IF (.NOT. ieee_is_finite(alpha_poly)) THEN
      target%reason = 'Nonfinite carbon enzyme allocation'
      RETURN
    END IF
    alpha_poly = MIN(MAX(alpha_poly, min_enzyme_fraction), 1.0_lp - min_enzyme_fraction)
    alpha_residue = 1.0_lp - alpha_poly
    enzyme_poly = alpha_poly * enzyme_depoly
    enzyme_residue = MAX(alpha_residue * enzyme_depoly, coefficients%persistent_enzyme_carbon)
    target%alpha_poly = alpha_poly
    target%alpha_residue = alpha_residue
    target%enzyme_depoly = enzyme_depoly
    target%enzyme_poly = enzyme_poly
    target%enzyme_residue = enzyme_residue
    depolymerisation_rate = coefficients%vmax_depoly_eff * [enzyme_poly/(coefficients%km_depoly_eff+ &
      & enzyme_poly),enzyme_residue/(coefficients%km_depoly_eff+enzyme_residue)]
    IF (ANY(depolymerisation_rate <= 0.0_lp)) THEN
      ! No arbitrary stock is assigned to an undetermined zero-source/zero-removal pool.
      ! Keep the live column and continue forward replay, which safely handles zero enzymes.
      target%reason = 'Zero enzyme removal: no unique stationary polymeric/residue target'
      RETURN
    END IF
    !> Eqs. (13),(14),(18): dL=IL+fW*IW, dR=(1-fR)*kB*B, L=dL/aL, R=dR/aR,
    !> with aL=vL*EL/(KE+EL), aR=vR*ER/(KE+ER). Require positive removal coefficients;
    !> Positive-source candidates must also remain above the forward eps8 depolymerization thresholds.
    !> The residue relation requires zero NET residue sorption, not absence of gross sorption/desorption.
    !> d_litter and d_residue are source rates [mol C m-3 s-1], not annual pool increments.
    polymeric_carbon = d_litter/depolymerisation_rate(1)
    residue_carbon = d_residue/depolymerisation_rate(2)

    ! Solve the actual finite-timestep MAOC exchange, not an unrestricted Langmuir extrapolation.
    CALL solve_analytical_maoc_partition(coefficients%dtime,carbon_pools([6,5]), &
      & coefficients%adsorption_rate,coefficients%desorption_rate, &
      & coefficients%qmax_org_layer,coefficients%volume_pore_total/coefficients%volume_per_carbon_maoc, &
      & carbon_pools(7:8),net_sorption_rate,sorption_increment,l_maoc_solution)
    target%carbon = carbon_pools
    IF (.NOT. l_maoc_solution) THEN
      target%reason = 'MAOC partition is degenerate or bounded solve did not close'
      RETURN
    END IF
    microbial_uptake = coefficients%vmax_uptake_eff*microbial_carbon*dom_carbon/(coefficients%km_uptake_eff+dom_carbon)
    target%uptake_step = microbial_uptake*coefficients%dtime
    !> Eq. (11), T_D=jD-ellD*D and other T_X=0; JD/JR are net adsorption [mol C/m3/s]:
    !> rS=IS-kS*S; rW=IW-kW*W; rL=IL+fW*kW*W-aL*L; rB=e*U-M;
    !> rR=(1-fR)*M-aR*R-JR; rD=fS*kS*S+aL*L+aR*R+fR*M-U-JD+jD-ellD*D; rAD=JD; rAR=JR.
    !> Their sum is Eq. (12): IS+IW+IL-(1-fS)*kS*S-(1-fW)*kW*W-(1-e)*U+jD-ellD*D.
    !> Internal transfers cancel only with these signs. Each pool must close, not merely their sum.
    target%residual(1:2) = coefficients%litter_input_carbon(1:2)-coefficients%litter_turnover_rate*carbon_pools(1:2)
    target%residual(3) = input_poly+coefficients%carbon_transfer_fraction(2)* &
      & coefficients%litter_turnover_rate(2)*woody_carbon-depolymerisation_rate(1)*polymeric_carbon
    target%residual(4) = coefficients%microbial_cue*microbial_uptake-microbial_decay
    target%residual(5) = (1.0_lp-coefficients%carbon_transfer_fraction(3))*microbial_decay- &
      & depolymerisation_rate(2)*residue_carbon-net_sorption_rate(2)
    target%residual(6) = coefficients%carbon_transfer_fraction(1)*coefficients%litter_turnover_rate(1)* &
      & soluable_carbon+SUM(depolymerisation_rate*carbon_pools([3,5])) + &
      & coefficients%carbon_transfer_fraction(3)*microbial_decay-microbial_uptake-net_sorption_rate(1)+ &
      & coefficients%transport_dom_inflow-coefficients%dom_export_rate*dom_carbon
    target%residual(7:8) = net_sorption_rate
    !> Eqs. (3)-(4): VA=gamma_A*(AD+AR), VP_raw=gamma_P*(S+L+B+R), free=max(0,V-VA).
    !> For admissible MAOC (VA<=V), surface b=min(1,free/VP_raw), or b=1 for negligible VP_raw.
    !> Excess=max(0,VP_raw-free) is NOT deleted C. At the surface it remains above mineral soil;
    !> at depth this solver rejects it rather than solving the displacement map of Eq. (25).
    !> Burial does not multiply B,R,L in the biochemical equations: they use full pool concentrations.
    target%maoc_volume = coefficients%volume_per_carbon_maoc*SUM(carbon_pools(7:8))
    target%poc_volume = coefficients%volume_per_carbon_poc*SUM(carbon_pools([1,3,4,5]))
    volume_pore_available = MAX(0.0_lp,coefficients%volume_pore_total-target%maoc_volume)
    target%poc_excess = MAX(0.0_lp,target%poc_volume-volume_pore_available)
    IF (coefficients%l_surface_layer .AND. target%poc_volume > eps8) &
      & target%burial = MIN(1.0_lp,volume_pore_available/target%poc_volume)

    !> Eq. (7): unrestricted M(B) requires B-M*dt>=B_min for either mortality law.
    !> U*dt>D is diagnostic only: continuous substrate_input_rate can replenish DOC within the interval.
    !> Forward replay caps actual microbial_uptake against DOC availability. This can cause adjustment away
    !> from the continuous carbon root, but is not itself a reason to reject that root.
    !> Additional C/N/P availability caps in calc_microbial_growth are NOT proved inactive here.
    !> The DOC activation cutoff is shared with forward uptake (1.e-16 mol C/m3), not generic eps8.
    !> This admits low-DOC deep-soil roots without imposing a DOC floor or changing their carbon balance.
    !> The conservative donor checks kS*dt,kW*dt,aL*dt<=1 and dR*dt+ads_R_step<=R avoid
    !> requiring more old substrate than exists; the latter is a sufficient guard, not a complete
    !> reconstruction of the staged forward CNP update including concurrent formation and transport.
    target%status = 3
    ! Report the first failing guard in the original order. Internal formatted
    ! writes keep this pure solver independent of logging; the caller supplies
    ! tile/cell/layer context. Each message fits the existing 192-character field.
    IF (.NOT. target%l_minimum_microbial_seed) THEN
      IF (microbial_carbon-microbial_decay*coefficients%dtime < coefficients%minimum_microbial_carbon) THEN
        WRITE(target%reason,'(A,5(1X,ES16.8))') &
          & 'Microbial floor: B,B-M*dt,Bmin [mol C/m3], M [mol C/m3/s], dt [s]=', &
          & microbial_carbon,microbial_carbon-microbial_decay*coefficients%dtime, &
          & coefficients%minimum_microbial_carbon,microbial_decay,coefficients%dtime
        RETURN
      END IF
      IF (dom_carbon <= min_dom_uptake_concentration) THEN
        WRITE(target%reason,'(A,2(1X,ES16.8))') &
          & 'DOC activation: D <= cutoff [mol C/m3]; D,cutoff=',dom_carbon,min_dom_uptake_concentration
        RETURN
      END IF
    END IF
    ! Match forward depolymerisation: only the linear regression path uses eps8
    ! as a pool activation cutoff. Density-dependent mortality permits tiny positive pools.
    IF (coefficients%mortality_exponent == 1._lp .AND. &
      & d_litter > 0.0_lp .AND. polymeric_carbon <= eps8) THEN
      WRITE(target%reason,'(A,3(1X,ES16.8))') &
        & 'Poly pool cutoff: positive dL [mol C/m3/s], L <= eps8 [mol C/m3]; dL,L,eps8=', &
        & d_litter,polymeric_carbon,eps8
      RETURN
    END IF
    IF (coefficients%mortality_exponent == 1._lp .AND. &
      & d_residue > 0.0_lp .AND. residue_carbon <= eps8) THEN
      WRITE(target%reason,'(A,3(1X,ES16.8))') &
        & 'Residue pool cutoff: positive dR [mol C/m3/s], R <= eps8 [mol C/m3]; dR,R,eps8=', &
        & d_residue,residue_carbon,eps8
      RETURN
    END IF
    IF (ANY(coefficients%litter_turnover_rate*coefficients%dtime > 1.0_lp)) THEN
      WRITE(target%reason,'(A,3(1X,ES16.8))') &
        & 'Litter turnover: kS*dt or kW*dt > 1; kS*dt,kW*dt [-], dt [s]=', &
        & coefficients%litter_turnover_rate*coefficients%dtime,coefficients%dtime
      RETURN
    END IF
    IF (coefficients%dom_export_rate*coefficients%dtime > 1.0_lp) THEN
      WRITE(target%reason,'(A,3(1X,ES16.8))') &
        & 'DOM export: rate*dt > 1; rate [s-1], dt [s], rate*dt [-]=', &
        & coefficients%dom_export_rate,coefficients%dtime,coefficients%dom_export_rate*coefficients%dtime
      RETURN
    END IF
    IF (depolymerisation_rate(1)*coefficients%dtime > 1.0_lp) THEN
      WRITE(target%reason,'(A,3(1X,ES16.8))') &
        & 'Poly removal: aL*dt > 1; aL [s-1], dt [s], aL*dt [-]=', &
        & depolymerisation_rate(1),coefficients%dtime,depolymerisation_rate(1)*coefficients%dtime
      RETURN
    END IF
    IF (d_residue*coefficients%dtime+sorption_increment(2) > residue_carbon) THEN
      WRITE(target%reason,'(A,5(1X,ES16.8))') &
        & 'Residue donor: dR*dt+ads > R; dR*dt,ads,R,total [mol C/m3], dt [s]=', &
        & d_residue*coefficients%dtime,sorption_increment(2),residue_carbon, &
        & d_residue*coefficients%dtime+sorption_increment(2),coefficients%dtime
      RETURN
    END IF
    target%status = 4
    target%reason = 'Deep solid POC requires displacement: no local transport-free equilibrium'
    IF (.NOT. coefficients%l_surface_layer .AND. target%poc_excess > 1.e-12_lp) RETURN
    target%status = 5
    target%reason = 'Nonfinite, negative or nonclosing conditional carbon target'
    IF (.NOT. ALL(ieee_is_finite(carbon_pools)) .OR. ANY(carbon_pools < 0.0_lp)) RETURN
    IF (.NOT. ALL(ieee_is_finite(target%residual))) RETURN
    ! Numerical verification tolerance, not the much looser ecological spinup-convergence threshold.
    ! Use each pool's own throughput so large DOM cycling cannot conceal a small MAOC imbalance.
    residual_scale = 1.e-18_lp + 1.e-9_lp * &
      & [coefficients%litter_input_carbon(1:2),d_litter,microbial_decay,d_residue,MAX(carbon_source,microbial_uptake), &
      & sorption_increment/coefficients%dtime]
    ! The floor seed closes all donor/associated pools, but deliberately retains e*U in microbial storage.
    ! Do not zero that diagnostic or include it in respiration; the column budget accounts for storage.
    IF (target%l_minimum_microbial_seed) residual_scale(4) = ABS(target%residual(4)) + 1.e-18_lp
    IF (ANY(ABS(target%residual) > residual_scale)) RETURN
    target%status = 0
    target%reason = 'Conditional C root; prescribed biology, no solid transport; nutrient equilibrium not solved'
    IF (target%l_minimum_microbial_seed) &
      & target%reason = 'Minimum microbial quasi-steady seed: zero mortality/residue/aRes; positive growth retained'
    END ASSOCIATE
  END SUBROUTINE solve_local_carbon_target

  !> Solve stationary MAOC partition using the forward finite-timestep shared caps.
  PURE SUBROUTINE solve_analytical_maoc_partition(dtime, donor_carbon, adsorption_rate, desorption_rate, &
      & qmax_org_layer, maoc_pore_capacity, &
      & associated_carbon, sorption_residual, sorption_increment, l_maoc_solution)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    REAL(lp), INTENT(in) :: dtime, donor_carbon(2), adsorption_rate(2), desorption_rate(2) !< [s], [mol C/m3], [m3/mol C/s], [s-1]
    !> Chemical and physical carbon capacities [mol C/m3]
    REAL(lp), INTENT(in) :: qmax_org_layer, maoc_pore_capacity
    REAL(lp), INTENT(out) :: associated_carbon(2) !< Stationary associated_carbon DOM/residue [mol C m-3]
    REAL(lp), INTENT(out) :: sorption_residual(2) !< Sorption minus desorption tendency [mol C m-3 s-1]
    REAL(lp), INTENT(out) :: sorption_increment(2) !< Gross sorption increments [mol C m-3 timestep-1]
    LOGICAL, INTENT(out) :: l_maoc_solution !< Both original equations and shared capacity bounds passed [-]
    REAL(lp) :: maoc_lower_bound, maoc_upper_bound, maoc_carbon !< Bracket endpoints and total MAOC [mol C m-3]
    !> Chemical and physical free capacities [mol C m-3]
    REAL(lp) :: chemical_capacity_available, physical_capacity_available
    !> Shared sorption limiter and desorbed fractions per timestep [-]
    REAL(lp) :: sorption_capacity_fraction, desorption_fraction(2)
    !> Potential increments and implied stationary partitions [mol C m-3]
    REAL(lp) :: potential_sorption(2), associated_carbon_candidate(2)
    INTEGER :: iter_maoc !< Bounded bisection counter [-]

    l_maoc_solution = .FALSE.
    associated_carbon = 0.0_lp
    sorption_residual = 0.0_lp
    sorption_increment = 0.0_lp
    IF (.NOT. ALL(ieee_is_finite([dtime,donor_carbon,adsorption_rate,desorption_rate,qmax_org_layer,maoc_pore_capacity]))) RETURN
    IF (dtime <= 0.0_lp .OR. MINVAL([donor_carbon,adsorption_rate,qmax_org_layer,maoc_pore_capacity]) < 0.0_lp) RETURN
    IF (ANY(desorption_rate <= 0.0_lp)) RETURN
    !> Eq. (31): delta_j=min(1,bj*dtime)>0 allows Aj=ads_j_step/delta_j.
    !> Do not divide by a zero desorption coefficient: an immobile associated_carbon stock can be history dependent.
    desorption_fraction = MIN(1.0_lp,desorption_rate*dtime)
    maoc_lower_bound = 0.0_lp
    maoc_upper_bound = MIN(qmax_org_layer,maoc_pore_capacity)
    !> Scalar equation g(A)=sum[s(A)*pj(A)/delta_j]-A=0 with 0<=A<=min(Q,KV).
    !> At A=0, g>=0; at the upper bound at least one shared capacity is zero, so g=-A<=0.
    !> The cap belongs in the root interval and shared limiter, NOT in the chemical potential sorption.
    DO iter_maoc = 1,160
      maoc_carbon = maoc_lower_bound+0.5_lp*(maoc_upper_bound-maoc_lower_bound)
      !> Eq. (29): HQ=max(0,Q-A), HV=max(0,KV-A), pj=Xj*min(1,aj*HQ*dtime).
      !> Xj is old-state D or R. Each pj is an increment [mol C/m3], not a rate [mol C/m3/s].
      chemical_capacity_available = MAX(0.0_lp,qmax_org_layer-maoc_carbon)
      physical_capacity_available = MAX(0.0_lp,maoc_pore_capacity-maoc_carbon)
      potential_sorption = donor_carbon*MIN(1.0_lp,adsorption_rate*chemical_capacity_available*dtime)
      sorption_capacity_fraction = 1.0_lp
      !> Eq. (30): s=min(1,HQ/p,HV/p), p=sum(pj). If p=0 keep s=1 and adsorption=0;
      !> this branch prevents 0/0 without inventing a nonzero sorption flux.
      IF (SUM(potential_sorption) > 0.0_lp) sorption_capacity_fraction = MIN(1.0_lp, &
      & chemical_capacity_available/SUM(potential_sorption),physical_capacity_available/SUM(potential_sorption))
      associated_carbon_candidate = sorption_capacity_fraction*potential_sorption/desorption_fraction
      IF (SUM(associated_carbon_candidate) > maoc_carbon) THEN
        maoc_lower_bound = maoc_carbon
      ELSE
        maoc_upper_bound = maoc_carbon
      END IF
    END DO
    !> Eqs. (19)-(20) would give Aj=Q*(aj/bj)*Xj/(1+sum[(aj/bj)*Xj]) ONLY when
    !> aj*HQ*dtime<=1, bj*dtime<=1 and p<=HQ,HV (no source, desorption or shared-chemical_capacity_available cap changes fluxes).
    !> Even Q<=KV does not guarantee inactive chemical capacity limitation. With one donor and an
    !> active chemical_capacity_available cap, sufficient potential_sorption sorption and b*dtime<1, Eq. (32) instead gives
    !> A=min(Q,KV)/(1+b*dtime). Do not use that one-pool special case for a two-pool partition.
    ! Normalize the partition to the bounded scalar root, then reevaluate the ORIGINAL operator there.
    ! This prevents near-saturation floating-point error from producing an over-capacity candidate.
    IF (SUM(associated_carbon_candidate) > 0.0_lp) associated_carbon = maoc_carbon* &
      & (associated_carbon_candidate/SUM(associated_carbon_candidate))
    maoc_carbon = SUM(associated_carbon)
    chemical_capacity_available = MAX(0.0_lp,qmax_org_layer-maoc_carbon)
    physical_capacity_available = MAX(0.0_lp,maoc_pore_capacity-maoc_carbon)
    potential_sorption = donor_carbon*MIN(1.0_lp,adsorption_rate*chemical_capacity_available*dtime)
    sorption_capacity_fraction = 1.0_lp
    IF (SUM(potential_sorption) > 0.0_lp) sorption_capacity_fraction = MIN(1.0_lp, &
      & chemical_capacity_available/SUM(potential_sorption),physical_capacity_available/SUM(potential_sorption))
    sorption_increment = sorption_capacity_fraction*potential_sorption
    !> Eq. (31): rj=(s*pj-delta_j*Aj)/dtime=0 [mol C/m3/s]. Reevaluate both original equations:
    !> a small scalar bracket alone cannot certify the partition when rates/headrooms are very small.
    sorption_residual = (sorption_increment-desorption_fraction*associated_carbon)/dtime
    l_maoc_solution = ALL(ieee_is_finite(sorption_residual)) .AND. &
      & ALL(ABS(sorption_residual) <= 1.e-18_lp+1.e-9_lp*MAX(sorption_increment,desorption_fraction* &
      & associated_carbon)/dtime) .AND. &
      & maoc_carbon <= MIN(qmax_org_layer,maoc_pore_capacity)+1.e-12_lp
  END SUBROUTINE solve_analytical_maoc_partition

  !> Reconstruct organic elements on scratch arrays after every layer passes.
  !> Scale each existing pool by Cnew/Cold, retaining N/C, P/C and isotope ratios.
  !> C14/C preservation is a NON-equilibrium seed; radioactive decay remains in replay.
  !> Undefined positive-target composition rejects the column without changing live state.
  !> Element increments [mol m-2] are artificial jump bookkeeping, never physical fluxes.
  PURE SUBROUTINE reconstruct_analytical_column(pool_before, targets, depth, pool_index, carbon_index, &
      & pool_candidate, element_increment, l_column_accepted, reason, failed_layer, failed_pool)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    REAL(lp), INTENT(in) :: pool_before(:,:,:), depth(:) !< pool/element/layer [mol m-3], thickness [m]
    TYPE(t_sb_local_target), INTENT(in) :: targets(:) !< Active-layer carbon targets
    INTEGER, INTENT(in) :: pool_index(:), carbon_index !< S,W,L,B,R[,D,AD,AR] and C matrix indices [-]
    REAL(lp), INTENT(out) :: pool_candidate(:,:,:), element_increment(:) !< Candidate [mol m-3], artificial jump [mol m-2]
    LOGICAL, INTENT(out) :: l_column_accepted !< All layers and elements validated [-]
    CHARACTER(len=*), INTENT(out) :: reason !< Column acceptance reason
    INTEGER, OPTIONAL, INTENT(out) :: failed_layer, failed_pool !< Layer/local pool index; 0 if not localized [-]
    INTEGER :: is, ipool_local, ipool, ielem !< Layer, local pool, matrix pool and element indices [-]
    REAL(lp) :: carbon_scaling_factor, carbon_before, carbon_target !< Scaling [-], pre/post carbon concentrations [mol C m-3]

    pool_candidate = pool_before
    element_increment = 0.0_lp
    l_column_accepted = .FALSE.
    IF (PRESENT(failed_layer)) failed_layer = 0
    IF (PRESENT(failed_pool)) failed_pool = 0
    reason = 'Rejected: at least one layer lacks an admissible transport-coupled carbon target'
    IF (ANY(targets%status /= 0)) RETURN
    reason = 'Rejected: invalid pre-jump composition or geometry'
    IF (.NOT. ALL(ieee_is_finite(pool_before)) .OR. .NOT. ALL(ieee_is_finite(depth))) RETURN
    IF (ANY(pool_before < 0.0_lp) .OR. ANY(depth <= 0.0_lp)) RETURN
    DO is = 1,SIZE(targets)
      DO ipool_local = 1,SIZE(pool_index)
        ipool = pool_index(ipool_local)
        carbon_before = pool_before(ipool,carbon_index,is)
        carbon_target = targets(is)%carbon(ipool_local)
        IF (carbon_before == 0.0_lp) THEN
          IF (carbon_target > 0.0_lp .OR. ANY(pool_before(ipool,:,is) /= 0.0_lp)) THEN
            pool_candidate = pool_before
            reason = 'Rejected: zero pool_before carbon leaves target pool stoichiometry undefined'
            IF (PRESENT(failed_layer)) failed_layer = is
            IF (PRESENT(failed_pool)) failed_pool = ipool_local
            RETURN
          END IF
          CYCLE
        END IF
        ! Guard the ratio before division under debug floating-point traps.
        IF (carbon_before < 1.0_lp) THEN
          IF (carbon_target > HUGE(1.0_lp)*carbon_before) THEN
            pool_candidate = pool_before
            reason = 'Rejected: unrepresentable pool scaling'
            IF (PRESENT(failed_layer)) failed_layer = is
            IF (PRESENT(failed_pool)) failed_pool = ipool_local
            RETURN
          END IF
        END IF
        carbon_scaling_factor = carbon_target/carbon_before
        IF (carbon_scaling_factor > 1.0_lp) THEN
          IF (ANY(pool_before(ipool,:,is) > HUGE(1.0_lp)/carbon_scaling_factor)) THEN
            pool_candidate = pool_before
            reason = 'Rejected: element scaling would overflow'
            IF (PRESENT(failed_layer)) failed_layer = is
            IF (PRESENT(failed_pool)) failed_pool = ipool_local
            RETURN
          END IF
        END IF
        pool_candidate(ipool,:,is) = pool_before(ipool,:,is)*carbon_scaling_factor
        pool_candidate(ipool,carbon_index,is) = carbon_target
      END DO
    END DO
    IF (.NOT. ALL(ieee_is_finite(pool_candidate)) .OR. ANY(pool_candidate < 0.0_lp)) THEN
      pool_candidate = pool_before
      reason = 'Rejected: reconstructed inventory is nonfinite or negative'
      RETURN
    END IF
    DO ielem = 1,SIZE(pool_before,2)
      element_increment(ielem) = SUM(SUM(pool_candidate(:,ielem,:)-pool_before(:,ielem,:),DIM=1)*depth)
    END DO
    l_column_accepted = .TRUE.
    reason = 'Accepted: advective C target; old organic ratios retained, including non-equilibrium 14C/C'
  END SUBROUTINE reconstruct_analytical_column

  !> Nonlinear microbial-DOC balances: e*U=M and J=(1-e)*M/e+ell*D.
  !> Per-capita mortality q=k_ref*(B/B_ref)**(beta-1)=e*v*D/(K+D).
  !> For ell>0 bracket D in (0,J/ell) in log space and recover B from q.
  !> For ell=0 use M=e*J/(1-e), then B from the mortality law and D=K*q/(e*v-q).
  !> Reject the uptake saturation pole, unrepresentable roots and failed balances.
  !> The caller handles the B_min boundary and finite-step donor/pore restrictions.
  !> Uses the approved outer 1e-18+1e-9*throughput tolerance, not a spinup criterion.
  PURE SUBROUTINE solve_density_dependent_microbial_target(carbon_source, export_rate, cue, uptake_rate, &
      & uptake_half_saturation, reference_rate, reference_biomass, exponent, biomass, dom, mortality, valid)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    REAL(lp), INTENT(in) :: carbon_source !< Effective litter plus upstream DOC [mol C m-3 s-1]
    REAL(lp), INTENT(in) :: export_rate, uptake_rate, reference_rate !< DOC export, uptake and reference death [s-1]
    REAL(lp), INTENT(in) :: cue, exponent !< Effective CUE and mortality exponent [-]
    REAL(lp), INTENT(in) :: uptake_half_saturation, reference_biomass !< Concentrations [mol C m-3]
    REAL(lp), INTENT(out) :: biomass, dom !< Unrestricted stationary concentrations [mol C m-3]
    REAL(lp), INTENT(out) :: mortality !< Stationary mortality flux [mol C m-3 s-1]
    LOGICAL, INTENT(out) :: valid !< Both balances close and the finite positive root is resolved [-]
    REAL(lp) :: log_lower, log_upper, log_dom, log_biomass !< Logarithmic concentration bracket/work values [-]
    REAL(lp) :: per_capita_death, growth_capacity !< q=M/B and e*v [s-1]
    REAL(lp) :: removal, uptake, tolerance !< Total removal, DOC uptake and balance allowance [mol C m-3 s-1]
    INTEGER :: iteration !< Scalar root iteration [-]

    valid = .FALSE.
    biomass = 0._lp
    dom = 0._lp
    mortality = 0._lp
    IF (.NOT. ALL(ieee_is_finite([carbon_source,export_rate,cue,uptake_rate, &
        & uptake_half_saturation,reference_rate,reference_biomass,exponent]))) RETURN
    IF (MIN(carbon_source,export_rate) < 0._lp .OR. cue <= 0._lp .OR. cue >= 1._lp .OR. &
        & MIN(uptake_rate,uptake_half_saturation,reference_rate,reference_biomass) <= 0._lp .OR. &
        & exponent <= 1._lp) RETURN
    IF (carbon_source == 0._lp) THEN
      valid = .TRUE.
      RETURN
    END IF
    growth_capacity = cue * uptake_rate
    IF (export_rate == 0._lp) THEN
      mortality = cue * carbon_source / (1._lp-cue)
      log_biomass = (LOG(mortality)-LOG(reference_rate)+ &
        & (exponent-1._lp)*LOG(reference_biomass))/exponent
      IF (log_biomass >= LOG(HUGE(1._lp)) .OR. log_biomass <= LOG(TINY(1._lp))) RETURN
      biomass = EXP(log_biomass)
      per_capita_death = mortality / biomass
      IF (per_capita_death >= growth_capacity) RETURN
      dom = uptake_half_saturation * per_capita_death / (growth_capacity-per_capita_death)
    ELSE
      log_lower = LOG(TINY(1._lp))
      log_upper = LOG(carbon_source)-LOG(export_rate)
      IF (log_upper >= LOG(HUGE(1._lp)) .OR. log_upper <= log_lower) RETURN
      DO iteration = 1,160
        log_dom = log_lower + 0.5_lp*(log_upper-log_lower)
        dom = EXP(log_dom)
        per_capita_death = growth_capacity * (dom/(uptake_half_saturation+dom))
        IF (per_capita_death <= 0._lp) THEN
          log_lower = log_dom
          CYCLE
        END IF
        ! q rounding to e*v loses the finite-root information; shrink the upper bracket.
        IF (per_capita_death >= growth_capacity) THEN
          log_upper = log_dom
          CYCLE
        END IF
        log_biomass = LOG(reference_biomass) + &
          & (LOG(per_capita_death)-LOG(reference_rate))/(exponent-1._lp)
        IF (log_biomass >= LOG(HUGE(1._lp))) THEN
          log_upper = log_dom
          CYCLE
        END IF
        biomass = 0._lp
        IF (log_biomass > LOG(TINY(1._lp))) biomass = EXP(log_biomass)
        mortality = per_capita_death * biomass
        removal = (1._lp-cue)*mortality/cue + export_rate*dom
        IF (removal > carbon_source) THEN
          log_upper = log_dom
        ELSE
          log_lower = log_dom
        END IF
      END DO
    END IF
    IF (.NOT. ALL(ieee_is_finite([biomass,dom,mortality]))) RETURN
    IF (MIN(biomass,dom,mortality) <= 0._lp) RETURN
    ! Independently verify the mortality law, not only uptake and total carbon balance.
    per_capita_death = reference_rate * EXP((exponent-1._lp)*(LOG(biomass)-LOG(reference_biomass)))
    IF (ABS(per_capita_death-mortality/biomass) > &
        & 1.e-9_lp*MAX(per_capita_death,mortality/biomass)) RETURN
    uptake = uptake_rate * biomass * (dom/(uptake_half_saturation+dom))
    removal = (1._lp-cue)*mortality/cue + export_rate*dom
    tolerance = 1.e-18_lp + 1.e-9_lp*MAX(carbon_source,removal)
    valid = ABS(removal-carbon_source) <= tolerance .AND. &
      & ABS(cue*uptake-mortality) <= 1.e-18_lp+1.e-9_lp*MAX(cue*uptake,mortality)
  END SUBROUTINE solve_density_dependent_microbial_target

#endif
END MODULE mo_sb_analytical_accelerator
