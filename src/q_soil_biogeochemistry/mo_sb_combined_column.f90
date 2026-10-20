!> Conservative event-local column arithmetic for combined soil acceleration
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
MODULE mo_sb_combined_column
#ifndef __NO_JSBACH__
  USE mo_kind, ONLY: lp
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: combine_sb_columns, redistribute_sb_columns
  PUBLIC :: combine_sb_structure, redistribute_sb_structure, combine_sb_physiology
CONTAINS
  !> Combine concentrations over positive-area owners with identical active geometry.
  !! Caller selects vegetated SB owners only. Fractions are absolute grid fractions;
  !! output concentrations refer to their combined area, not the whole grid box.
  !! Equal layer thickness makes this concentration average exactly an inventory
  !! average. Different geometries are rejected rather than silently interpolated.
  !! This routine neither calibrates structure nor changes live tile memory.
  PURE SUBROUTINE combine_sb_columns(fraction, thickness, pools, weights, combined, accepted)
    REAL(lp), INTENT(in) :: fraction(:) !< Absolute participating tile areas [-]
    REAL(lp), INTENT(in) :: thickness(:,:) !< Active layer thickness: layer,tile [m]
    REAL(lp), INTENT(in) :: pools(:,:,:,:) !< Pool,element,layer,tile [mol element m-3]
    REAL(lp), INTENT(out) :: weights(SIZE(fraction)) !< Normalized participating areas [-]
    REAL(lp), INTENT(out) :: combined(SIZE(pools,1),SIZE(pools,2),SIZE(pools,3))
      !< Combined pool concentrations [mol element m-3]
    LOGICAL, INTENT(out) :: accepted !< Valid finite inputs and common geometry [-]
    REAL(lp) :: area !< Total represented grid fraction [-]
    INTEGER :: it, first !< Tile indices [-]

    accepted = .FALSE.
    weights = 0._lp
    combined = 0._lp
    IF (SIZE(pools,4) /= SIZE(fraction)) RETURN
    IF (SIZE(thickness,2) /= SIZE(fraction) .OR. SIZE(thickness,1) /= SIZE(pools,3)) RETURN
    IF (.NOT. ALL(ieee_is_finite(fraction))) RETURN
    IF (ANY(fraction < 0._lp) .OR. ANY(fraction > 1._lp)) RETURN
    area = SUM(fraction)
    IF (area <= 0._lp) RETURN
    first = 0
    DO it = 1, SIZE(fraction)
      IF (fraction(it) == 0._lp) CYCLE
      IF (.NOT. ALL(ieee_is_finite(thickness(:,it)))) RETURN
      IF (ANY(thickness(:,it) <= 0._lp)) RETURN
      IF (.NOT. ALL(ieee_is_finite(pools(:,:,:,it)))) RETURN
      IF (ANY(pools(:,:,:,it) < 0._lp)) RETURN
      IF (first == 0) first = it
      ! Geometry comes from the same registered grid; no approximate remapping.
      IF (ANY(thickness(:,it) /= thickness(:,first))) RETURN
    END DO
    weights = fraction / area
    DO it = 1, SIZE(fraction)
      IF (weights(it) == 0._lp) CYCLE
      combined = combined + weights(it) * pools(:,:,:,it)
    END DO
    accepted = ALL(ieee_is_finite(combined))
  END SUBROUTINE combine_sb_columns

  !> Build, but never commit, a conservative redistribution candidate.
  !! For positive pre-event parent stock, every element uses old_parent/mean_parent.
  !! The parent is carbon for organic pools, or the parent nutrient for an inorganic
  !! pool plus its isotope. Consequently final composition is shared across PFTs,
  !! while original parent-stock shares are retained. Completely empty pools use
  !! uniform concentrations. Parent-free nonzero elements are rejected explicitly.
  !! This is only the arithmetic step: caller MUST check pore/sorption/isotope
  !! admissibility on the entire candidate before committing any live state.
  !! Nonparticipants are copied unchanged. On rejection, discard the candidate.
  PURE SUBROUTINE redistribute_sb_columns(fraction, thickness, before, target, parent, candidate, accepted)
    REAL(lp), INTENT(in) :: fraction(:) !< Absolute participating tile areas [-]
    REAL(lp), INTENT(in) :: thickness(:,:) !< Active layer thickness: layer,tile [m]
    REAL(lp), INTENT(in) :: before(:,:,:,:) !< Pre-event pool,element,layer,tile [mol element m-3]
    REAL(lp), INTENT(in) :: target(:,:,:) !< Final combined pool,element,layer [mol element m-3]
    INTEGER, INTENT(in) :: parent !< Parent-element index common to the supplied pool group [-]
    REAL(lp), INTENT(out) :: candidate(SIZE(before,1),SIZE(before,2),SIZE(before,3),SIZE(before,4))
      !< Proposed PFT concentrations, not yet physically validated [mol element m-3]
    LOGICAL, INTENT(out) :: accepted !< Arithmetic redistribution is admissible [-]
    REAL(lp) :: weights(SIZE(fraction)) !< Normalized tile areas [-]
    REAL(lp) :: mean_before(SIZE(before,1),SIZE(before,2),SIZE(before,3))
      !< Pre-event area mean [mol element m-3]
    REAL(lp) :: scale, old_parent, mean_parent !< Share multiplier [-]; parent concentrations [mol m-3]
    INTEGER :: it, layer, pool !< Indices [-]
    LOGICAL :: valid !< Aggregation input validation [-]

    candidate = before
    accepted = .FALSE.
    IF (parent < 1 .OR. parent > SIZE(before,2)) RETURN
    IF (ANY(SHAPE(target) /= [SIZE(before,1), SIZE(before,2), SIZE(before,3)])) RETURN
    CALL combine_sb_columns(fraction, thickness, before, weights, mean_before, valid)
    IF (.NOT. valid) RETURN
    IF (.NOT. ALL(ieee_is_finite(target))) RETURN
    IF (ANY(target < 0._lp)) RETURN
    DO layer = 1, SIZE(before,3)
      DO pool = 1, SIZE(before,1)
        mean_parent = mean_before(pool,parent,layer)
        DO it = 1, SIZE(fraction)
          IF (weights(it) == 0._lp) CYCLE
          old_parent = before(pool,parent,layer,it)
          ! Do not silently erase nutrients/isotopes from carbon-free material.
          IF (old_parent == 0._lp .AND. ANY(before(pool,:,layer,it) > 0._lp)) RETURN
          IF (target(pool,parent,layer) == 0._lp .AND. ANY(target(pool,:,layer) > 0._lp)) RETURN
          scale = 1._lp
          IF (mean_parent > 0._lp) THEN
            ! Guard division/multiplication before they can trigger trapping FP.
            IF (mean_parent < 1._lp) THEN
              IF (old_parent > HUGE(1._lp) * mean_parent) RETURN
            END IF
            scale = old_parent / mean_parent
          END IF
          IF (scale > 1._lp) THEN
            IF (ANY(target(pool,:,layer) > HUGE(1._lp) / scale)) RETURN
          END IF
          candidate(pool,:,layer,it) = target(pool,:,layer) * scale
        END DO
      END DO
    END DO
    accepted = .TRUE.
  END SUBROUTINE redistribute_sb_columns

  !> Reconstruct the independent structure, not its dependent density/capacity diagnostics.
  !! Pore fraction is volume per TOTAL soil volume. With common layer geometry
  !! and mineral particle density, the mineral-mass weight is area*(1-pore).
  !! Thus mean mineral mass times mean capacity per mineral mass reproduces the
  !! original area-mean chemical capacity. Averaging specific capacity by area
  !! alone would not conserve that capacity when PFT structures differ.
  !! The caller must subsequently recompute occupancy, bulk density, burial and
  !! effective sorption capacities with the ordinary forward-model routines.
  PURE SUBROUTINE combine_sb_structure(weights, pore, specific_capacity, shared_pore, &
    & shared_specific_capacity, accepted)
    REAL(lp), INTENT(in) :: weights(:) !< Normalized participating areas from combine_sb_columns [-]
    REAL(lp), INTENT(in) :: pore(:,:) !< Total pore fraction: layer,tile [m3 m-3]
    REAL(lp), INTENT(in) :: specific_capacity(:,:,:) !< Capacity,layer,tile [mol kg-1 mineral]
    REAL(lp), INTENT(out) :: shared_pore(SIZE(pore,1)) !< Combined pore fraction [m3 m-3]
    REAL(lp), INTENT(out) :: shared_specific_capacity(SIZE(specific_capacity,1),SIZE(pore,1))
      !< Combined specific chemical capacities [mol kg-1 mineral]
    LOGICAL, INTENT(out) :: accepted !< Independent state is finite and physically defined [-]
    REAL(lp) :: area, solid, mineral_weight !< Area sum, mean solid fraction, weighted solid fraction [-]
    INTEGER :: it, layer !< Tile and layer indices [-]

    accepted = .FALSE.
    shared_pore = 0._lp
    shared_specific_capacity = 0._lp
    IF (SIZE(pore,2) /= SIZE(weights)) RETURN
    IF (SIZE(specific_capacity,2) /= SIZE(pore,1)) RETURN
    IF (SIZE(specific_capacity,3) /= SIZE(weights)) RETURN
    IF (.NOT. ALL(ieee_is_finite(weights))) RETURN
    IF (ANY(weights < 0._lp) .OR. ANY(weights > 1._lp)) RETURN
    area = SUM(weights)
    IF (area <= 0._lp) RETURN
    DO it = 1, SIZE(weights)
      IF (weights(it) == 0._lp) CYCLE
      IF (.NOT. ALL(ieee_is_finite(pore(:,it)))) RETURN
      ! A mineral-free column cannot define capacities per mineral mass.
      IF (ANY(pore(:,it) < 0._lp) .OR. ANY(pore(:,it) >= 1._lp)) RETURN
      IF (.NOT. ALL(ieee_is_finite(specific_capacity(:,:,it)))) RETURN
      IF (ANY(specific_capacity(:,:,it) < 0._lp)) RETURN
      shared_pore = shared_pore + (weights(it) / area) * pore(:,it)
    END DO
    DO layer = 1, SIZE(pore,1)
      solid = 0._lp
      DO it = 1, SIZE(weights)
        IF (weights(it) == 0._lp) CYCLE
        solid = solid + (weights(it) / area) * (1._lp - pore(layer,it))
      END DO
      IF (solid <= 0._lp) RETURN
      DO it = 1, SIZE(weights)
        IF (weights(it) == 0._lp) CYCLE
        mineral_weight = ((weights(it) / area) * (1._lp - pore(layer,it))) / solid
        shared_specific_capacity(:,layer) = shared_specific_capacity(:,layer) + &
          & mineral_weight * specific_capacity(:,layer,it)
      END DO
    END DO
    accepted = ALL(ieee_is_finite(shared_specific_capacity))
  END SUBROUTINE combine_sb_structure

  !> Map shared calibration to PFT mineral solid fractions without clipping.
  !! Apply the same relative change to every participating mineral solid
  !! fraction: s_i,new=s_i,old*s_shared,new/mean(s_old), where s=1-pore.
  !! This preserves pre-event mineral shares and the final area-mean pore volume.
  !! Bounds are supplied by the existing physical scheme, not new tolerances.
  !! Reject the WHOLE candidate if any PFT violates them. Do not redistribute
  !! only the feasible PFTs or clip independently: either loses the shared target.
  PURE SUBROUTINE redistribute_sb_structure(weights, pore_before, shared_pore, lower, upper, &
    & pore_candidate, accepted)
    REAL(lp), INTENT(in) :: weights(:) !< Participating areas, normalized internally [-]
    REAL(lp), INTENT(in) :: pore_before(:,:) !< Pre-event layer,tile pore fraction [m3 m-3]
    REAL(lp), INTENT(in) :: shared_pore(:) !< Final combined layer pore fraction [m3 m-3]
    REAL(lp), INTENT(in) :: lower(:,:), upper(:,:) !< Per-PFT physical pore bounds [m3 m-3]
    REAL(lp), INTENT(out) :: pore_candidate(SIZE(pore_before,1),SIZE(pore_before,2))
      !< Proposed structure, valid only when accepted [m3 m-3]
    LOGICAL, INTENT(out) :: accepted !< All participating structures satisfy physical bounds [-]
    REAL(lp) :: area, solid, scale, new_solid !< Area, mean solid fraction, multiplier, new fraction [-]
    INTEGER :: layer, it !< Layer and tile indices [-]

    accepted = .FALSE.
    pore_candidate = pore_before
    IF (SIZE(pore_before,2) /= SIZE(weights)) RETURN
    IF (SIZE(shared_pore) /= SIZE(pore_before,1)) RETURN
    IF (ANY(SHAPE(lower) /= SHAPE(pore_before))) RETURN
    IF (ANY(SHAPE(upper) /= SHAPE(pore_before))) RETURN
    IF (.NOT. ALL(ieee_is_finite(weights))) RETURN
    IF (ANY(weights < 0._lp) .OR. ANY(weights > 1._lp)) RETURN
    area = SUM(weights)
    IF (area <= 0._lp) RETURN
    IF (.NOT. ALL(ieee_is_finite(shared_pore))) RETURN
    IF (ANY(shared_pore < 0._lp) .OR. ANY(shared_pore >= 1._lp)) RETURN
    DO it = 1, SIZE(weights)
      IF (weights(it) == 0._lp) CYCLE
      IF (.NOT. ALL(ieee_is_finite(pore_before(:,it)))) RETURN
      IF (ANY(pore_before(:,it) < 0._lp) .OR. ANY(pore_before(:,it) >= 1._lp)) RETURN
      IF (.NOT. ALL(ieee_is_finite(lower(:,it)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(upper(:,it)))) RETURN
      IF (ANY(lower(:,it) < 0._lp) .OR. ANY(upper(:,it) >= 1._lp)) RETURN
      IF (ANY(lower(:,it) > upper(:,it))) RETURN
    END DO
    DO layer = 1, SIZE(pore_before,1)
      solid = 0._lp
      DO it = 1, SIZE(weights)
        IF (weights(it) == 0._lp) CYCLE
        solid = solid + (weights(it) / area) * (1._lp - pore_before(layer,it))
      END DO
      IF (solid <= 0._lp) RETURN
      IF (solid < 1._lp) THEN
        IF (1._lp - shared_pore(layer) > HUGE(1._lp) * solid) RETURN
      END IF
      scale = (1._lp - shared_pore(layer)) / solid
      DO it = 1, SIZE(weights)
        IF (weights(it) == 0._lp) CYCLE
        new_solid = scale * (1._lp - pore_before(layer,it))
        IF (new_solid > 1._lp) RETURN
        pore_candidate(layer,it) = 1._lp - new_solid
        IF (pore_candidate(layer,it) < lower(layer,it)) RETURN
        IF (pore_candidate(layer,it) > upper(layer,it)) RETURN
      END DO
    END DO
    accepted = .TRUE.
  END SUBROUTINE redistribute_sb_structure

  !> Combine physiological state using microbial C inventory weights at each layer.
  !! Geometry has already been checked by combine_sb_columns. Normalize biomass
  !! by its maximum BEFORE weighting to avoid overflow or underflow of an area
  !! times tiny biomass. If all active microbial stocks are exactly zero, use
  !! area weights. No artificial microbial floor is introduced by aggregation.
  !! The caller selects only physiological states/histories, never pool stocks,
  !! cumulative forcing budgets, event counters or soil structural quantities.
  PURE SUBROUTINE combine_sb_physiology(weights, microbial_c, state, shared_state, accepted)
    REAL(lp), INTENT(in) :: weights(:) !< Participating areas [-]
    REAL(lp), INTENT(in) :: microbial_c(:,:) !< Layer,tile biomass [mol C m-3 total soil]
    REAL(lp), INTENT(in) :: state(:,:,:) !< State,layer,tile; units retained from each input field
    REAL(lp), INTENT(out) :: shared_state(SIZE(state,1),SIZE(state,2)) !< Combined state; input units
    LOGICAL, INTENT(out) :: accepted !< Finite inputs and positive represented area [-]
    REAL(lp) :: local_weights(SIZE(weights)), maximum, denominator !< Normalized weights, C, weight sum
    INTEGER :: it, layer !< Tile and layer indices [-]

    accepted = .FALSE.
    shared_state = 0._lp
    IF (SIZE(microbial_c,2) /= SIZE(weights) .OR. SIZE(state,3) /= SIZE(weights)) RETURN
    IF (SIZE(state,2) /= SIZE(microbial_c,1)) RETURN
    IF (.NOT. ALL(ieee_is_finite(weights))) RETURN
    IF (ANY(weights < 0._lp) .OR. ANY(weights > 1._lp)) RETURN
    IF (SUM(weights) <= 0._lp) RETURN
    DO it = 1, SIZE(weights)
      IF (weights(it) == 0._lp) CYCLE
      IF (.NOT. ALL(ieee_is_finite(microbial_c(:,it)))) RETURN
      IF (ANY(microbial_c(:,it) < 0._lp)) RETURN
      IF (.NOT. ALL(ieee_is_finite(state(:,:,it)))) RETURN
    END DO
    DO layer = 1, SIZE(state,2)
      maximum = 0._lp
      DO it = 1, SIZE(weights)
        IF (weights(it) > 0._lp) maximum = MAX(maximum, microbial_c(layer,it))
      END DO
      local_weights = weights
      IF (maximum > 0._lp) THEN
        DO it = 1, SIZE(weights)
          IF (weights(it) == 0._lp) CYCLE
          local_weights(it) = weights(it) * (microbial_c(layer,it) / maximum)
        END DO
      END IF
      denominator = SUM(local_weights)
      IF (denominator <= 0._lp) RETURN
      DO it = 1, SIZE(weights)
        IF (local_weights(it) == 0._lp) CYCLE
        shared_state(:,layer) = shared_state(:,layer) + &
          & (local_weights(it) / denominator) * state(:,layer,it)
      END DO
    END DO
    accepted = ALL(ieee_is_finite(shared_state))
  END SUBROUTINE combine_sb_physiology
#endif
END MODULE mo_sb_combined_column
