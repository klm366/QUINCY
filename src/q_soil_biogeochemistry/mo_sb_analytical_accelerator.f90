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

  IMPLICIT NONE

  PRIVATE
  ! Scientific targets used only by the shared replay event driver.
  PUBLIC :: t_sb_local_coefficients, t_sb_local_target, solve_local_carbon_target
  PUBLIC :: solve_analytical_maoc_partition, reconstruct_analytical_column
  PUBLIC :: t_sb_analytical_inputs, build_sb_analytical_inputs, diagnose_local_carbon_targets
  PUBLIC :: solve_density_dependent_microbial_target
  PUBLIC :: solve_simple_carbon_target

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
    REAL(lp), ALLOCATABLE :: rate_product_mean(:,:,:) !< depoly/uptake/sorption/desorption/HSC/SSM rt*rm [-]
    REAL(lp), ALLOCATABLE :: cue(:,:), alpha_poly(:,:), alpha_ap(:,:) !< Coupled-window means [-]
    REAL(lp), ALLOCATABLE :: water_rate_mean(:,:,:) !< Downward/lateral hydrological export rates [s-1]
    REAL(lp), ALLOCATABLE :: litter_budget_residual(:,:) !< Partitioned minus source element rate [mol m-2 s-1]
    REAL(lp), ALLOCATABLE :: depth(:,:), pore_volume(:,:) !< Thickness [m], fixed pore fraction [-]
    REAL(lp), ALLOCATABLE :: qmax_org(:,:) !< Chemical sorption capacity [mol C m-3 soil]
  END TYPE t_sb_analytical_inputs

CONTAINS

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
