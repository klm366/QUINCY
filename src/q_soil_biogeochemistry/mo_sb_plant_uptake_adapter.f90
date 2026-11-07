!> Root uptake aggregation for a shared replay soil column
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
MODULE mo_sb_plant_uptake_adapter
#ifndef __NO_JSBACH__
  USE mo_kind, ONLY: lp
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_sb_plant_uptake, prepare_combined_plant_uptake, valid_sb_plant_uptake
  PUBLIC :: t_sb_plant_uptake_drivers

  !> One seasonal bin of event-local PFT drivers, before shared soil evaluation.
  !! The combined-column constructor fills these from participating contexts.
  !! No soil solute pools belong here: nutrient supply is always the live shared
  !! column. Demand contains recorded fractions BEFORE annual replay correction.
  TYPE :: t_sb_plant_uptake_drivers
    REAL(lp), ALLOCATABLE :: fraction(:,:) !< Cell,PFT participating grid fractions [-]
    REAL(lp), ALLOCATABLE :: roots(:,:,:) !< Cell,layer,PFT roots [mol C m-3 total soil]
    REAL(lp), ALLOCATABLE :: temperature(:,:,:) !< Cell,layer,PFT temperature response [-]
    REAL(lp), ALLOCATABLE :: phi_leaf_min(:) !< PFT minimum leaf water potential [MPa]
    REAL(lp), ALLOCATABLE :: vmax(:,:) !< N/P,PFT uptake capacity [umol X mol-1 C s-1]
    REAL(lp), ALLOCATABLE :: demand(:,:,:) !< N/P,cell,PFT recorded demand fractions [-]
  END TYPE t_sb_plant_uptake_drivers

  !> Temporary uptake coefficients; never registered state or a second soil model.
  !! All concentrations are per total soil volume of the combined vegetated area.
  !! Demand is NOT included in potential capacity. Apply it only after the shared
  !! soil kernel has limited uptake by supply, matching ordinary replay ordering.
  TYPE :: t_sb_plant_uptake
    CHARACTER(len=16) :: scheme = '' !< simple_1d, eca_none, eca_part or eca_full [-]
    REAL(lp), ALLOCATABLE :: root(:,:) !< Sum(w_i R_i), cell,layer [mol C m-3]
    REAL(lp), ALLOCATABLE :: activity(:,:) !< Sum(w_i R_i t_i m_i), cell,layer [mol C m-3]
    REAL(lp), ALLOCATABLE :: capacity(:,:,:) !< Sum(w_i R_i t_i m_i V_Xi), N/P,cell,layer [umol X m-3 s-1]
    REAL(lp), ALLOCATABLE :: demand(:,:,:) !< Capacity-weighted effective N/P demand, N/P,cell,layer [-]
  END TYPE t_sb_plant_uptake
CONTAINS
  !> Form grouped environmental response and nutrient-specific uptake coefficients.
  !! Inputs are already climatological drivers; this does not recover covariances
  !! lost during time averaging. Roots must first be converted from layer stocks
  !! to concentrations using the common positive layer thickness.
  !!
  !! m_i=max(1-psi/phi_i,0). Summing PFT contributions is equivalent to grouping
  !! the two current phi values, without hard-coding their values or PFT indices.
  !! Unlike choosing the dominant PFT, it varies continuously with participating
  !! root inventories. Temperature acclimation can differ among PFTs.
  !!
  !! ECA: total roots enter ALL competition denominators; activity enters only
  !! the plant numerator. Non-ECA: capacity includes PFT-specific Vmax. Each
  !! kernel retains its own substrate response and units conversion. Do not
  !! apply separate supply limits to PFTs/groups sharing one nutrient reservoir.
  !!
  !! Preserve the correlation of demand with potential uptake coefficients:
  !! D_X=sum(a_Xi D_Xi)/sum(a_Xi). Shared concentration factors cancel here,
  !! including a common supply limiter. For ECA a_Xi=w_i R_i t_i m_i;
  !! for non-ECA additionally multiply by V_Xi. Zero capacity gives zero demand.
  !! Annual litter/uptake feedback must be applied to each input demand BEFORE
  !! aggregation and bounded to [0,1], not multiplied into potential capacity.
  PURE SUBROUTINE prepare_combined_plant_uptake(scheme, fraction, roots, temperature, water_potential, &
      & phi_leaf_min, vmax, demand, uptake, accepted, detail)
    CHARACTER(len=*), OPTIONAL, INTENT(out) :: detail !< First rejected grouped uptake driver
    CHARACTER(len=*), INTENT(in) :: scheme !< Existing soil/adsorption scheme [-]
    REAL(lp), INTENT(in) :: fraction(:,:) !< Cell,PFT absolute participating area [-]
    REAL(lp), INTENT(in) :: roots(:,:,:) !< Cell,layer,PFT root C [mol C m-3]
    REAL(lp), INTENT(in) :: temperature(:,:,:) !< Cell,layer,PFT uptake temperature modifier [-]
    REAL(lp), INTENT(in) :: water_potential(:,:) !< Shared cell,layer soil water potential [MPa]
    REAL(lp), INTENT(in) :: phi_leaf_min(:) !< PFT minimum leaf water potential [MPa], strictly negative
    REAL(lp), INTENT(in) :: vmax(:,:) !< N/P,PFT non-ECA uptake capacity [umol X mol-1 C s-1]
    REAL(lp), INTENT(in) :: demand(:,:,:) !< N/P,cell,PFT effective seasonal demand [0,1]
    TYPE(t_sb_plant_uptake), INTENT(out) :: uptake !< Candidate coefficients, usable only when accepted
    LOGICAL, INTENT(out) :: accepted !< Finite, nonnegative and shape-consistent inputs [-]
    INTEGER :: nc, nl, nt, ic, il, it, ix !< Dimensions and indices [-]
    REAL(lp) :: area, weight, response, contribution, coefficient !< Intermediate weights and capacities
    LOGICAL :: eca !< Shared-kinetic rather than PFT-Vmax formulation [-]

    accepted = .FALSE.
    IF (PRESENT(detail)) detail='unsupported scheme, incompatible dimensions or invalid area fractions'
    SELECT CASE (scheme)
    CASE ('eca_part','eca_full')
      eca = .TRUE.
    CASE ('eca_none','simple_1d')
      eca = .FALSE.
    CASE DEFAULT
      RETURN
    END SELECT
    nc = SIZE(roots,1)
    nl = SIZE(roots,2)
    nt = SIZE(roots,3)
    IF (ANY(SHAPE(fraction) /= [nc,nt])) RETURN
    IF (ANY(SHAPE(temperature) /= SHAPE(roots))) RETURN
    IF (ANY(SHAPE(water_potential) /= [nc,nl])) RETURN
    IF (SIZE(phi_leaf_min) /= nt) RETURN
    IF (ANY(SHAPE(vmax) /= [2,nt]) .OR. ANY(SHAPE(demand) /= [2,nc,nt])) RETURN
    IF (.NOT. ALL(ieee_is_finite(fraction))) RETURN
    IF (ANY(fraction < 0._lp) .OR. ANY(fraction > 1._lp)) RETURN
    ALLOCATE(uptake%root(nc,nl), uptake%activity(nc,nl), &
      & uptake%capacity(2,nc,nl), uptake%demand(2,nc,nl))
    uptake%scheme = scheme
    uptake%root = 0._lp
    uptake%activity = 0._lp
    uptake%capacity = 0._lp
    uptake%demand = 0._lp
    DO ic = 1, nc
      area = SUM(fraction(ic,:))
      IF (area <= 0._lp) CYCLE
      IF (.NOT. ALL(ieee_is_finite(water_potential(ic,:)))) THEN
        IF (PRESENT(detail)) THEN
          il=FINDLOC(.NOT.ieee_is_finite(water_potential(ic,:)),.TRUE.,DIM=1)
          WRITE(detail,'(A,2I5,A,ES22.14)') 'water potential member/layer=',0,il,' value=',water_potential(ic,il)
        END IF
        RETURN
      END IF
      DO it = 1, nt
        IF (fraction(ic,it) == 0._lp) CYCLE
        IF (.NOT. ALL(ieee_is_finite(roots(ic,:,it)))) THEN
          IF (PRESENT(detail)) THEN
            il=FINDLOC(.NOT.ieee_is_finite(roots(ic,:,it)),.TRUE.,DIM=1)
            WRITE(detail,'(A,2I5,A,ES22.14)') 'root C member/layer=',it,il,' value=',roots(ic,il,it)
          END IF
          RETURN
        END IF
        IF (ANY(roots(ic,:,it) < 0._lp)) THEN
          IF (PRESENT(detail)) THEN
            il=FINDLOC(roots(ic,:,it)<0._lp,.TRUE.,DIM=1)
            WRITE(detail,'(A,2I5,A,ES22.14)') 'root C member/layer=',it,il,' value=',roots(ic,il,it)
          END IF
          RETURN
        END IF
        IF (.NOT. ALL(ieee_is_finite(temperature(ic,:,it)))) THEN
          IF (PRESENT(detail)) THEN
            il=FINDLOC(.NOT.ieee_is_finite(temperature(ic,:,it)),.TRUE.,DIM=1)
            WRITE(detail,'(A,2I5,A,ES22.14)') 'temperature modifier member/layer=',it,il,' value=',temperature(ic,il,it)
          END IF
          RETURN
        END IF
        IF (ANY(temperature(ic,:,it) < 0._lp)) THEN
          IF (PRESENT(detail)) THEN
            il=FINDLOC(temperature(ic,:,it)<0._lp,.TRUE.,DIM=1)
            WRITE(detail,'(A,2I5,A,ES22.14)') 'temperature modifier member/layer=',it,il,' value=',temperature(ic,il,it)
          END IF
          RETURN
        END IF
        IF (.NOT. ieee_is_finite(phi_leaf_min(it))) THEN
          IF (PRESENT(detail)) WRITE(detail,'(A,I5,A,2ES22.14)') 'phi_leaf_min member=',it,' value=',phi_leaf_min(it)
          RETURN
        END IF
        IF (phi_leaf_min(it) >= 0._lp) THEN
          IF (PRESENT(detail)) WRITE(detail,'(A,I5,A,2ES22.14)') 'phi_leaf_min member=',it,' value=',phi_leaf_min(it)
          RETURN
        END IF
        IF (.NOT. ALL(ieee_is_finite(demand(:,ic,it)))) THEN
          IF (PRESENT(detail)) WRITE(detail,'(A,I5,A,2ES22.14)') 'N/P demand member=',it,' value=',demand(:,ic,it)
          RETURN
        END IF
        IF (ANY(demand(:,ic,it) < 0._lp) .OR. ANY(demand(:,ic,it) > 1._lp)) THEN
          IF (PRESENT(detail)) WRITE(detail,'(A,I5,A,2ES22.14)') 'N/P demand member=',it,' value=',demand(:,ic,it)
          RETURN
        END IF
        IF (.NOT. eca) THEN
          IF (.NOT. ALL(ieee_is_finite(vmax(:,it)))) THEN
            IF (PRESENT(detail)) WRITE(detail,'(A,I5,A,2ES22.14)') 'N/P Vmax member=',it,' value=',vmax(:,it)
            RETURN
          END IF
          IF (ANY(vmax(:,it) < 0._lp)) THEN
            IF (PRESENT(detail)) WRITE(detail,'(A,I5,A,2ES22.14)') 'N/P Vmax member=',it,' value=',vmax(:,it)
            RETURN
          END IF
        END IF
        weight = fraction(ic,it) / area
        DO il = 1, nl
          response = MAX(1._lp - water_potential(ic,il) / phi_leaf_min(it), 0._lp)
          contribution = weight * roots(ic,il,it) * temperature(ic,il,it) * response
          uptake%root(ic,il) = uptake%root(ic,il) + weight * roots(ic,il,it)
          uptake%activity(ic,il) = uptake%activity(ic,il) + contribution
          DO ix = 1, 2
            coefficient = contribution
            IF (.NOT. eca) coefficient = coefficient * vmax(ix,it)
            IF (.NOT. eca) uptake%capacity(ix,ic,il) = uptake%capacity(ix,ic,il) + coefficient
            uptake%demand(ix,ic,il) = uptake%demand(ix,ic,il) + coefficient * demand(ix,ic,it)
          END DO
        END DO
      END DO
      DO il = 1, nl
        DO ix = 1, 2
          coefficient = uptake%activity(ic,il)
          IF (.NOT. eca) coefficient = uptake%capacity(ix,ic,il)
          IF (coefficient > 0._lp) THEN
            uptake%demand(ix,ic,il) = uptake%demand(ix,ic,il) / coefficient
          ELSE
            uptake%demand(ix,ic,il) = 0._lp
          END IF
        END DO
      END DO
    END DO
    accepted = valid_sb_plant_uptake(uptake,nc,nl,scheme)
    IF (.NOT. accepted .AND. PRESENT(detail)) THEN
      WRITE(detail,'(A,4ES23.15)') 'grouped output: min root/activity, min/max demand=', &
        & MINVAL(uptake%root),MINVAL(uptake%activity),MINVAL(uptake%demand),MAXVAL(uptake%demand)
    END IF
  END SUBROUTINE prepare_combined_plant_uptake

  !> Validate an optional adapter before any model memory is modified.
  PURE LOGICAL FUNCTION valid_sb_plant_uptake(uptake,nc,nl,scheme) RESULT(valid)
    TYPE(t_sb_plant_uptake), INTENT(in) :: uptake !< Temporary shared uptake coefficients
    INTEGER, INTENT(in) :: nc, nl !< Expected cells and soil layers [-]
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: scheme !< Expected kinetic scheme [-]
    valid = .FALSE.
    IF (.NOT. ALLOCATED(uptake%root) .OR. .NOT. ALLOCATED(uptake%activity)) RETURN
    IF (.NOT. ALLOCATED(uptake%capacity) .OR. .NOT. ALLOCATED(uptake%demand)) RETURN
    IF (ANY(SHAPE(uptake%root) /= [nc,nl])) RETURN
    IF (ANY(SHAPE(uptake%activity) /= [nc,nl])) RETURN
    IF (ANY(SHAPE(uptake%capacity) /= [2,nc,nl])) RETURN
    IF (ANY(SHAPE(uptake%demand) /= [2,nc,nl])) RETURN
    IF (PRESENT(scheme)) THEN
      IF (uptake%scheme /= scheme) RETURN
    END IF
    IF (.NOT. ALL(ieee_is_finite(uptake%root))) RETURN
    IF (.NOT. ALL(ieee_is_finite(uptake%activity))) RETURN
    IF (.NOT. ALL(ieee_is_finite(uptake%capacity))) RETURN
    IF (.NOT. ALL(ieee_is_finite(uptake%demand))) RETURN
    IF (ANY(uptake%root < 0._lp) .OR. ANY(uptake%activity < 0._lp)) RETURN
    IF (ANY(uptake%capacity < 0._lp)) RETURN
    IF (ANY(uptake%demand < 0._lp) .OR. ANY(uptake%demand > 1._lp)) RETURN
    valid = .TRUE.
  END FUNCTION valid_sb_plant_uptake
#endif
END MODULE mo_sb_plant_uptake_adapter
