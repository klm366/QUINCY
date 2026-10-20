!> QUINCY litter transport and partitioning
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### routines for litter transport and partitioning for the both JSM and SSM (soil models in QUINCY)
!>
MODULE mo_q_sb_litter_processes
#ifndef __NO_JSBACH__

  USE mo_kind,            ONLY: lp
  USE mo_lnd_bgcm_idx

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_litter_partitioning
  PUBLIC :: combine_sb_litter_inputs

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_sb_litter_processes'

CONTAINS

  !> Partition individual PFT tissues before taking their area mean.
  !> Shared columns require common active layer geometry. Fractions need not
  !> sum to one: outputs are per represented vegetated area, not per grid box.
  !> Use the replay timestep to retain the forward partitioner's numerical
  !> branches. Conservation is checked in rates using the approved adapter
  !> tolerance [mol element m-2 s-1], independently for every element/cell.
  SUBROUTINE combine_sb_litter_inputs(nc, ns, npool, nelem, nsource, ntile, dt, scheme, &
      & fraction, nactive, depth, sla, growthform, roots, veg_rate, seed_rate, formation_rate, accepted)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    INTEGER, INTENT(in) :: nc, ns, npool, nelem, nsource, ntile !< Dimensions [-]
    REAL(lp), INTENT(in) :: dt !< Actual replay substep [s]
    CHARACTER(len=*), INTENT(in) :: scheme !< Existing JSM/simple_1d scheme [-]
    REAL(lp), INTENT(in) :: fraction(nc,ntile) !< Absolute PFT areas [-]
    REAL(lp), INTENT(in) :: nactive(nc), depth(nc,ns) !< Common active layers [-], thickness [m]
    REAL(lp), INTENT(in) :: sla(ntile) !< PFT specific leaf area [m2 mol C-1]
    INTEGER, INTENT(in) :: growthform(ntile) !< PFT growth form [-]
    REAL(lp), INTENT(in) :: roots(nc,ns,ntile) !< Individual root fractions [-]
    REAL(lp), INTENT(in) :: veg_rate(nsource,nelem,nc,ntile), seed_rate(nelem,nc,ntile)
      !< Tissue and seed-bed litter rates [mol element m-2 s-1]
    REAL(lp), INTENT(out) :: formation_rate(npool,nelem,nc,ns)
      !< Combined pool litter sources [mol element m-3 s-1]; zero if rejected
    LOGICAL, INTENT(out) :: accepted !< All geometry and elemental budget checks passed [-]
    REAL(lp) :: trial(npool,nelem,nc,ns), single(npool,nelem,1,ns)
    REAL(lp) :: reference, recovered, weight, area !< Rates [mol m-2 s-1], weights [-]
    INTEGER :: c, t, e, s !< Cell, tile, element, layer [-]

    accepted = .FALSE.
    formation_rate = 0._lp
    trial = 0._lp
    IF (.NOT. ieee_is_finite(dt)) RETURN
    IF (dt <= 0._lp) RETURN
    IF (.NOT. ALL(ieee_is_finite(fraction))) RETURN
    IF (ANY(fraction < 0._lp)) RETURN
    IF (.NOT. ALL(ieee_is_finite(nactive)) .OR. .NOT. ALL(ieee_is_finite(depth))) RETURN
    DO c = 1,nc
      IF (nactive(c) < 1._lp .OR. nactive(c) > REAL(ns,lp)) RETURN
      IF (nactive(c) /= REAL(NINT(nactive(c)),lp)) RETURN
      IF (ANY(depth(c,1:NINT(nactive(c))) <= 0._lp)) RETURN
      area = SUM(fraction(c,:))
      IF (.NOT. ieee_is_finite(area) .OR. area <= 0._lp) RETURN
      DO t = 1,ntile
        IF (fraction(c,t) == 0._lp) CYCLE
        IF (.NOT. ieee_is_finite(sla(t))) RETURN
        IF (.NOT. ALL(ieee_is_finite(roots(c,:,t)))) RETURN
        IF (ANY(roots(c,:,t) < 0._lp)) RETURN
        IF (.NOT. ALL(ieee_is_finite(veg_rate(:,:,c,t)))) RETURN
        IF (.NOT. ALL(ieee_is_finite(seed_rate(:,c,t)))) RETURN
        IF (ANY(veg_rate(:,:,c,t) < 0._lp) .OR. ANY(seed_rate(:,c,t) < 0._lp)) RETURN
        single = 0._lp
        CALL calc_litter_partitioning(1,ns,nactive(c:c),sla(t),growthform(t),scheme,depth(c:c,:), &
          & roots(c:c,:,t),veg_rate(:,:,c:c,t)*dt,seed_rate(:,c:c,t)*dt,single)
        IF (.NOT. ALL(ieee_is_finite(single)) .OR. ANY(single < 0._lp)) RETURN
        weight = fraction(c,t)/area
        DO e = 1,nelem
          reference = SUM(veg_rate(:,e,c,t)) + seed_rate(e,c,t)
          recovered = 0._lp
          DO s = 1,NINT(nactive(c))
            recovered = recovered + SUM(single(:,e,1,s))/dt*depth(c,s)
          END DO
          IF (ABS(recovered-reference) > 1.e-18_lp + 1.e-10_lp*reference) RETURN
        END DO
        trial(:,:,c,:) = trial(:,:,c,:) + weight*single(:,:,1,:)/dt
      END DO
    END DO
    IF (.NOT. ALL(ieee_is_finite(trial))) RETURN
    formation_rate = trial
    accepted = .TRUE.
  END SUBROUTINE combine_sb_litter_inputs

  ! ======================================================================================================= !
  !>Calculates the partitioning of litter fall from vegetation into metabolic, structural,
  !>and woody litter given
  !>  a) tissue specifc preferences (labile + reserve = metabolic; sapwood and heartwood = woody);
  !>  b) the tissue-specific lignin concentration and the actual tissue C:N ratio.
  !>All above-ground litter enters the first layer, only fine and coarse roots are distributed
  !>according to their root profile
  !>
  SUBROUTINE calc_litter_partitioning( &
    & nc, &                            ! in
    & nsoil_w, &
    & num_sl_above_bedrock, &
    & lctlib_sla, &
    & lctlib_growthform, &
    & sb_model_scheme, &
    & soil_depth_sl, &
    & root_fraction_sl, &
    & veg_litterfall_mt, &             ! in
    & seed_bed_litterfall_mt, &        ! in
    & sb_formation_mt, &               ! inout
    & prescribed_formation )           ! optional in: already partitioned combined replay litter

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    USE mo_exception, ONLY: finish
    USE mo_veg_constants,           ONLY: ITREE, IGRASS, ICROP
    USE mo_jsb_math_constants,      ONLY: eps8
    USE mo_sb_constants,            ONLY: fc_soluable_max, k_fc_soluable, k_fn_soluable, k_fp_soluable, &
                                          lc_leaf_max, lc_leaf2sla, lc_fine_root, lc_coarse_root, lc_fruit, lc_seed_bed
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc                        !< dimension
    INTEGER,                           INTENT(in)    :: nsoil_w                   !< nr soil layer
    REAL(lp), DIMENSION(nc),           INTENT(in)    :: num_sl_above_bedrock      !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(lp),                          INTENT(in)    :: lctlib_sla                !< lctlib parameter
    INTEGER,                           INTENT(in)    :: lctlib_growthform         !< lctlib parameter
    CHARACTER(len=*),                  INTENT(in)    :: sb_model_scheme           !< config option: soil model
    REAL(lp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: soil_depth_sl             !< soil layer thickness
    REAL(lp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: root_fraction_sl          !< root fraction per soil layer
    REAL(lp),                          INTENT(in)    :: veg_litterfall_mt(:,:,:)  !< bgcm flux: vegetation litterfall
    REAL(lp),                          INTENT(in)    :: seed_bed_litterfall_mt(:,:)  !< bgcm flux: litterfall from seed bed turnover
    REAL(lp),                          INTENT(inout) :: sb_formation_mt(:,:,:,:)  !< bgcm flux: soil biogeochemistry formation flux
    REAL(lp), OPTIONAL, INTENT(in) :: prescribed_formation(:,:,:,:)
      !< PFT-partitioned pool,element,cell,layer litter increments [mol m-3 timestep-1]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(lp), DIMENSION(nc) :: fc_soluable_leaf           !< fc: ...
    REAL(lp), DIMENSION(nc) :: fc_soluable_fine_root
    REAL(lp), DIMENSION(nc) :: fc_soluable_coarse_root
    REAL(lp), DIMENSION(nc) :: fc_soluable_fruit
    REAL(lp), DIMENSION(nc) :: fc_soluable_seed_bed
    REAL(lp), DIMENSION(nc) :: fc_soluable_labile
    REAL(lp), DIMENSION(nc) :: fc_soluable_reserve
    REAL(lp), DIMENSION(nc) :: fn_soluable_leaf           !< fn: ...
    REAL(lp), DIMENSION(nc) :: fn_soluable_fine_root
    REAL(lp), DIMENSION(nc) :: fn_soluable_coarse_root
    REAL(lp), DIMENSION(nc) :: fn_soluable_fruit
    REAL(lp), DIMENSION(nc) :: fn_soluable_seed_bed
    REAL(lp), DIMENSION(nc) :: fn_soluable_labile
    REAL(lp), DIMENSION(nc) :: fn_soluable_reserve
    REAL(lp), DIMENSION(nc) :: fp_soluable_leaf           !< fp: ...
    REAL(lp), DIMENSION(nc) :: fp_soluable_fine_root
    REAL(lp), DIMENSION(nc) :: fp_soluable_coarse_root
    REAL(lp), DIMENSION(nc) :: fp_soluable_fruit
    REAL(lp), DIMENSION(nc) :: fp_soluable_seed_bed
    REAL(lp), DIMENSION(nc) :: fp_soluable_labile
    REAL(lp), DIMENSION(nc) :: fp_soluable_reserve
    INTEGER                 :: ic, is                     !< loop over dimensions
    INTEGER                 :: isoil                      !< loop over soil layers
    INTEGER :: pool !< Validation index [-]
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_partitioning'
    ! ----------------------------------------------------------------------------------------------------- !

    ! A shared replay column must not partition mixed tissue with a representative
    ! PFT's SLA/growth form. Partition each PFT first, then add the conservative
    ! area mean here. Preserve other formation terms already calculated this step.
    ! Source budgets remain the combined raw litter matrices, also used by replay
    ! diagnostics. Exact shape/finite checks precede mutation of the formation.
    IF (PRESENT(prescribed_formation)) THEN
      IF (ANY(SHAPE(prescribed_formation) /= SHAPE(sb_formation_mt))) &
        & CALL finish(routine, 'Invalid prescribed litter shape.')
      IF (.NOT. ALL(ieee_is_finite(prescribed_formation))) &
        & CALL finish(routine, 'Non-finite prescribed litter.')
      IF (ANY(prescribed_formation < 0._lp)) CALL finish(routine, 'Negative prescribed litter.')
      DO pool = 1, SIZE(sb_formation_mt,1)
        IF (pool == ix_soluable_litter .OR. pool == ix_polymeric_litter .OR. pool == ix_woody_litter) CYCLE
        IF (ANY(prescribed_formation(pool,:,:,:) /= 0._lp)) &
          & CALL finish(routine, 'Prescribed litter contains a non-litter source.')
      END DO
      DO ic = 1, nc
        DO isoil = NINT(num_sl_above_bedrock(ic))+1, nsoil_w
          IF (ANY(prescribed_formation(:,:,ic,isoil) /= 0._lp)) &
            & CALL finish(routine, 'Prescribed litter below active soil.')
        END DO
      END DO
      ! The constructor checks flux-rate conservation with the approved rate
      ! tolerance. Do not introduce a timestep-dependent absolute tolerance here.
      sb_formation_mt = sb_formation_mt + prescribed_formation
      RETURN
    END IF

    !>0.9 init local variables
    !>  labile and reserve compartments are init with 1.0
    !>
    fc_soluable_leaf        = 0.0_lp
    fc_soluable_fine_root   = 0.0_lp
    fc_soluable_coarse_root = 0.0_lp
    fc_soluable_fruit       = 0.0_lp
    fc_soluable_seed_bed    = 0.0_lp
    fn_soluable_leaf        = 0.0_lp
    fn_soluable_fine_root   = 0.0_lp
    fn_soluable_coarse_root = 0.0_lp
    fn_soluable_fruit       = 0.0_lp
    fn_soluable_seed_bed    = 0.0_lp
    fp_soluable_leaf        = 0.0_lp
    fp_soluable_fine_root   = 0.0_lp
    fp_soluable_coarse_root = 0.0_lp
    fp_soluable_fruit       = 0.0_lp
    fp_soluable_seed_bed    = 0.0_lp
    ! labile & reserve
    fc_soluable_labile      = 1.0_lp
    fc_soluable_reserve     = 1.0_lp
    fn_soluable_labile      = 1.0_lp
    fn_soluable_reserve     = 1.0_lp
    fp_soluable_labile      = 1.0_lp
    fp_soluable_reserve     = 1.0_lp

    !>1.0 determine allocation of fresh litter fall to litter pools
    !>
    DO ic = 1,nc
      ! leaf
      IF(veg_litterfall_mt(ix_leaf, ixN, ic) > eps8) THEN
        fc_soluable_leaf(ic) = MAX(0.0_lp, fc_soluable_max - k_fc_soluable * (lc_leaf_max + lc_leaf2sla * lctlib_sla) &
          &                    * veg_litterfall_mt(ix_leaf, ixC, ic) / veg_litterfall_mt(ix_leaf, ixN, ic))
        IF(fc_soluable_leaf(ic) > eps8) THEN
          fn_soluable_leaf(ic) = 1._lp / (1._lp + (1._lp - fc_soluable_leaf(ic)) / (k_fn_soluable * fc_soluable_leaf(ic)))
          fp_soluable_leaf(ic) = 1._lp / (1._lp + (1._lp - fc_soluable_leaf(ic)) / (k_fp_soluable * fc_soluable_leaf(ic)))
        END IF
      END IF
      ! fine_root
      IF(veg_litterfall_mt(ix_fine_root, ixN, ic) > eps8) THEN
        fc_soluable_fine_root(ic) = MAX(0.0_lp, fc_soluable_max - k_fc_soluable * lc_fine_root &
          &                         * veg_litterfall_mt(ix_fine_root, ixC, ic) / veg_litterfall_mt(ix_fine_root, ixN, ic))
        IF(fc_soluable_fine_root(ic) > eps8) THEN
          fn_soluable_fine_root(ic) = 1._lp &
            &                         / (1._lp + (1._lp - fc_soluable_fine_root(ic)) / (k_fn_soluable * fc_soluable_fine_root(ic)))
          fp_soluable_fine_root(ic) = 1._lp &
            &                         / (1._lp + (1._lp - fc_soluable_fine_root(ic)) / (k_fp_soluable * fc_soluable_fine_root(ic)))
        END IF
      END IF
      ! coarse_root
      IF(veg_litterfall_mt(ix_coarse_root, ixN, ic) > eps8) THEN
        fc_soluable_coarse_root(ic) = MAX(0.0_lp, fc_soluable_max - k_fc_soluable * lc_coarse_root &
          &                           * veg_litterfall_mt(ix_coarse_root, ixC, ic) / veg_litterfall_mt(ix_coarse_root, ixN, ic))
        IF(fc_soluable_coarse_root(ic) > eps8) THEN
          fn_soluable_coarse_root(ic) = 1._lp / (1._lp + (1._lp - fc_soluable_coarse_root(ic)) &
            &                           / (k_fn_soluable*fc_soluable_coarse_root(ic)))
          fp_soluable_coarse_root(ic) = 1._lp / (1._lp + (1._lp - fc_soluable_coarse_root(ic)) &
            &                           / (k_fp_soluable*fc_soluable_coarse_root(ic)))
         END IF
      END IF
      ! fruit
      IF(veg_litterfall_mt(ix_fruit, ixN, ic) > eps8) THEN
        fc_soluable_fruit(ic) = MAX(0.0_lp, fc_soluable_max - k_fc_soluable * lc_fruit &
          &                     * veg_litterfall_mt(ix_fruit, ixC, ic) / veg_litterfall_mt(ix_fruit, ixN, ic))
        IF(fc_soluable_fruit(ic) > eps8) THEN
          fn_soluable_fruit(ic) = 1._lp / (1._lp + (1._lp - fc_soluable_fruit(ic)) / (k_fn_soluable * fc_soluable_fruit(ic)))
          fp_soluable_fruit(ic) = 1._lp / (1._lp + (1._lp - fc_soluable_fruit(ic)) / (k_fp_soluable * fc_soluable_fruit(ic)))
        END IF
      END IF
      ! seed_bed
      IF(seed_bed_litterfall_mt(ixN, ic) > eps8) THEN
        fc_soluable_seed_bed(ic) = MAX(0.0_lp, fc_soluable_max - k_fc_soluable * lc_seed_bed &
          &                        * seed_bed_litterfall_mt(ixC, ic) / seed_bed_litterfall_mt(ixN, ic))
        IF(fc_soluable_seed_bed(ic) > eps8) THEN
          fn_soluable_seed_bed(ic) = 1._lp &
            &                        / (1._lp + (1._lp - fc_soluable_seed_bed(ic)) / (k_fn_soluable * fc_soluable_seed_bed(ic)))
          fp_soluable_seed_bed(ic) = 1._lp &
            &                        / (1._lp + (1._lp - fc_soluable_seed_bed(ic)) / (k_fp_soluable * fc_soluable_seed_bed(ic)))
        END IF
      END IF
    END DO ! loop over nc

    !>2.0 Update surface metabolic and soluable litter pools with litter fall
    !>
    !>  2.1 soluable litter
    !>
    sb_formation_mt(ix_soluable_litter, ixC, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixC, :, 1) &
      &                              + (fc_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixC, :)      &
      &                                 + fc_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixC, :)     &
      &                                 + fc_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixC, :)    &
      &                                 + fc_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixC, :)   &
      &                                 + fc_soluable_seed_bed(:) * seed_bed_litterfall_mt(ixC, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixN, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixN, :, 1) &
      &                              + (fn_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixN, :)      &
      &                                 + fn_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixN, :)     &
      &                                 + fn_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixN, :)    &
      &                                 + fn_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixN, :)   &
      &                                 + fn_soluable_seed_bed(:) * seed_bed_litterfall_mt(ixN, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixP, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixP, :, 1) &
      &                              + (fp_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixP, :)      &
      &                                 + fp_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixP, :)     &
      &                                 + fp_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixP, :)    &
      &                                 + fp_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixP, :)   &
      &                                 + fp_soluable_seed_bed(:) * seed_bed_litterfall_mt(ixP, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixC13, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixC13, :, 1) &
      &                              + (fc_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixC13, :)      &
      &                                 + fc_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixC13, :)     &
      &                                 + fc_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixC13, :)    &
      &                                 + fc_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixC13, :)   &
      &                                 + fc_soluable_seed_bed(:) * seed_bed_litterfall_mt(ixC13, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixC14, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixC14, :, 1) &
      &                              + (fc_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixC14, :)      &
      &                                 + fc_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixC14, :)     &
      &                                 + fc_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixC14, :)    &
      &                                 + fc_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixC14, :)   &
      &                                 + fc_soluable_seed_bed(:) * seed_bed_litterfall_mt(ixC14, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixN15, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixN15, :, 1) &
      &                              + (fn_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixN15, :)      &
      &                                 + fn_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixN15, :)     &
      &                                 + fn_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixN15, :)    &
      &                                 + fn_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixN15, :)   &
      &                                 + fn_soluable_seed_bed(:) * seed_bed_litterfall_mt(ixN15, :)) &
      &                              / soil_depth_sl(:,1)
    !>  2.2 polymeric litter
    !>
    sb_formation_mt(ix_polymeric_litter, ixC, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixC, :, 1) &
      &                              + ((1.0_lp - fc_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixC, :)      &
      &                                 + (1.0_lp - fc_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixC, :)     &
      &                                 + (1.0_lp - fc_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixC, :)    &
      &                                 + (1.0_lp - fc_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixC, :)   &
      &                                 + (1.0_lp - fc_soluable_seed_bed(:)) * seed_bed_litterfall_mt(ixC, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixN, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixN, :, 1) &
      &                              + ((1.0_lp - fn_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixN, :)      &
      &                                 + (1.0_lp - fn_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixN, :)     &
      &                                 + (1.0_lp - fn_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixN, :)    &
      &                                 + (1.0_lp - fn_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixN, :)   &
      &                                 + (1.0_lp - fn_soluable_seed_bed(:)) * seed_bed_litterfall_mt(ixN, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixP, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixP, :, 1) &
      &                              + ((1.0_lp - fp_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixP, :)      &
      &                                 + (1.0_lp - fp_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixP, :)     &
      &                                 + (1.0_lp - fp_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixP, :)    &
      &                                 + (1.0_lp - fp_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixP, :)   &
      &                                 + (1.0_lp - fp_soluable_seed_bed(:)) * seed_bed_litterfall_mt(ixP, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) &
      &                              + ((1.0_lp - fc_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixC13, :)      &
      &                                 + (1.0_lp - fc_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixC13, :)     &
      &                                 + (1.0_lp - fc_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixC13, :)    &
      &                                 + (1.0_lp - fc_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixC13, :)   &
      &                                 + (1.0_lp - fc_soluable_seed_bed(:)) * seed_bed_litterfall_mt(ixC13, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) &
      &                              + ((1.0_lp - fc_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixC14, :)      &
      &                                 + (1.0_lp - fc_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixC14, :)     &
      &                                 + (1.0_lp - fc_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixC14, :)    &
      &                                 + (1.0_lp - fc_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixC14, :)   &
      &                                 + (1.0_lp - fc_soluable_seed_bed(:)) * seed_bed_litterfall_mt(ixC14, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) &
      &                              + ((1.0_lp - fn_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixN15, :)      &
      &                                 + (1.0_lp - fn_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixN15, :)     &
      &                                 + (1.0_lp - fn_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixN15, :)    &
      &                                 + (1.0_lp - fn_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixN15, :)   &
      &                                 + (1.0_lp - fn_soluable_seed_bed(:)) * seed_bed_litterfall_mt(ixN15, :)) &
      &                              / soil_depth_sl(:,1)

    !>3.0 Add surface litter fall from sapwood and heartwood pools to surface woody litter pool for trees
    !>    for grasses, stems and halms (=sapwood) are not actually wood, and therefore paritioned to
    !>    metabolic and structural litter as if their were fine roots
    SELECT CASE(lctlib_growthform)
    ! tree PFT
    CASE (ITREE)
      ! woody litter
      ! @TODO may be re-written to a one-liner using bgcm matrices
      sb_formation_mt(ix_woody_litter, ixC, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixC, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixC, :) + veg_litterfall_mt(ix_heart_wood, ixC, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixN, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixN, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixN, :) + veg_litterfall_mt(ix_heart_wood, ixN, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixP, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixP, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixP, :) + veg_litterfall_mt(ix_heart_wood, ixP, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixC13, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixC13, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixC13, :) + veg_litterfall_mt(ix_heart_wood, ixC13, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixC14, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixC14, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixC14, :) + veg_litterfall_mt(ix_heart_wood, ixC14, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixN15, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixN15, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixN15, :) + veg_litterfall_mt(ix_heart_wood, ixN15, :)) / soil_depth_sl(:,1)
    ! grass PFT
    CASE (IGRASS,ICROP)
      ! soluable litter
      sb_formation_mt(ix_soluable_litter, ixC, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixC, :, 1) &
        &             + fc_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC, :) + veg_litterfall_mt(ix_heart_wood, ixC, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixN, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixN, :, 1) &
        &             + fn_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixN, :) + veg_litterfall_mt(ix_heart_wood, ixN, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixP, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixP, :, 1) &
        &             + fp_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixP, :) + veg_litterfall_mt(ix_heart_wood, ixP, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixC13, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixC13, :, 1) &
        &             + fc_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC13, :) + veg_litterfall_mt(ix_heart_wood, ixC13, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixC14, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixC14, :, 1) &
        &             + fc_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC14, :) + veg_litterfall_mt(ix_heart_wood, ixC14, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixN15, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixN15, :, 1) &
        &             + fn_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixN15, :) + veg_litterfall_mt(ix_heart_wood, ixN15, :)) &
        &             / soil_depth_sl(:,1)
      ! polymeric litter
      sb_formation_mt(ix_polymeric_litter, ixC, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixC, :, 1) &
        &             + (1.0_lp - fc_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC, :) + veg_litterfall_mt(ix_heart_wood, ixC, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixN, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixN, :, 1) &
        &             + (1.0_lp - fn_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixN, :) + veg_litterfall_mt(ix_heart_wood, ixN, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixP, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixP, :, 1) &
        &             + (1.0_lp - fp_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixP, :) + veg_litterfall_mt(ix_heart_wood, ixP, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) &
        &             + (1.0_lp - fc_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC13, :) + veg_litterfall_mt(ix_heart_wood, ixC13, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) &
        &             + (1.0_lp - fc_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC14, :) + veg_litterfall_mt(ix_heart_wood, ixC14, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) &
        &             + (1.0_lp - fn_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixN15, :) + veg_litterfall_mt(ix_heart_wood, ixN15, :)) &
        &             / soil_depth_sl(:,1)
    END SELECT

    !>4.0 Distribute fine and coarse root litter to soil depth according to the root profile
    !>
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        ! soluable litter (fine root & coarse root)
        sb_formation_mt(ix_soluable_litter, ixC, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixC, ic, is) &
          &             + (fc_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixC, ic) &
          &             +  fc_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixC, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixN, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixN, ic, is) &
          &             + (fn_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixN, ic) &
          &             +  fn_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixN, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixP, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixP, ic, is) &
          &             + (fp_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixP, ic) &
          &             +  fp_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixP, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixC13, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixC13, ic, is) &
          &             + (fc_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixC13, ic) &
          &             +  fc_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixC13, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixC14, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixC14, ic, is) &
          &             + (fc_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixC14, ic) &
          &             +  fc_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixC14, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixN15, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixN15, ic, is) &
          &             + (fn_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixN15, ic) &
          &             +  fn_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixN15, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        ! simple soil model / jsm
        SELECT CASE(TRIM(sb_model_scheme))
        CASE ("simple_1d")
          ! polymeric litter (fine root & coarse root)
          sb_formation_mt(ix_polymeric_litter, ixC, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixC, ic, is) &
            &           + ((1.0_lp - fc_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixC, ic) &
            &           +  (1.0_lp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixN, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixN, ic, is) &
            &           + ((1.0_lp - fn_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixN, ic) &
            &           +  (1.0_lp - fn_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixN, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixP, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixP, ic, is) &
            &           + ((1.0_lp - fp_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixP, ic) &
            &           +  (1.0_lp - fp_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixP, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixC13, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixC13, ic, is) &
            &           + ((1.0_lp - fc_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixC13, ic) &
            &           +  (1.0_lp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC13, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixC14, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixC14, ic, is) &
            &           + ((1.0_lp - fc_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixC14, ic) &
            &           +  (1.0_lp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC14, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixN15, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixN15, ic, is) &
            &           + ((1.0_lp - fn_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixN15, ic) &
            &           +  (1.0_lp - fn_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixN15, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
        CASE ("jsm")
          ! polymeric litter (fine root)
          sb_formation_mt(ix_polymeric_litter, ixC, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixC, ic, is) &
            &           + ((1.0_lp - fc_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixC, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixN, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixN, ic, is) &
            &           + ((1.0_lp - fn_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixN, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixP, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixP, ic, is) &
            &           + ((1.0_lp - fp_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixP, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixC13, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixC13, ic, is) &
            &           + ((1.0_lp - fc_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixC13, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixC14, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixC14, ic, is) &
            &           + ((1.0_lp - fc_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixC14, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixN15, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixN15, ic, is) &
            &           + ((1.0_lp - fn_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixN15, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          ! polymeric litter (coarse root)
          sb_formation_mt(ix_woody_litter, ixC, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixC, ic, is) &
            &           + ((1.0_lp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixN, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixN, ic, is) &
            &           + ((1.0_lp - fn_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixN, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixP, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixP, ic, is) &
            &           + ((1.0_lp - fp_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixP, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixC13, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixC13, ic, is) &
            &           + ((1.0_lp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC13, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixC14, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixC14, ic, is) &
            &           + ((1.0_lp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC14, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixN15, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixN15, ic, is) &
            &           + ((1.0_lp - fn_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixN15, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
        END SELECT
      END DO
    END DO
  END SUBROUTINE calc_litter_partitioning

#endif
END MODULE mo_q_sb_litter_processes
