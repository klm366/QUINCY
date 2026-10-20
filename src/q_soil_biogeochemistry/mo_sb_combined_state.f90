!> Explicit persistent state inventory for transactional combined soil replay
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
MODULE mo_sb_combined_state
#ifndef __NO_JSBACH__
  USE mo_kind, ONLY: lp
  USE mo_jsb_physical_constants, ONLY: scale_correction_c14
  USE mo_jsb_tile_class, ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class, ONLY: t_jsb_task_options
  USE mo_sb_memory_class, ONLY: t_sb_memory
  USE mo_jsb_process_class, ONLY: SB_, HYDRO_
  USE mo_hydro_memory_class, ONLY: t_hydro_memory
  USE mo_lnd_bgcm_store_class, ONLY: SB_BGCM_POOL_ID
  USE mo_lnd_bgcm_idx, ONLY: ixC, ixN, ixC13, ixC14, ixN15, ix_microbial
  USE mo_sb_combined_column, ONLY: combine_sb_columns, combine_sb_structure, &
    & combine_sb_physiology, redistribute_sb_columns, redistribute_sb_structure
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_sb_combined_state, read_sb_combined_state, write_sb_combined_state
  PUBLIC :: reconstruct_sb_combined_state, redistribute_sb_combined_state
  PUBLIC :: valid_sb_combined_inventory

  ! Inorganic groups 1:6 are NH4(solute,associated), NO3, NOy, N2O, N2.
  ! Groups 7:11 are PO4(fast,slow,occluded,primary,solute).
  ! Component 1 is parent N/P, component 2 is N15 (zero for P).
  TYPE :: t_sb_combined_state
    REAL(lp), ALLOCATABLE :: organic(:,:,:) !< Pool,element,layer [mol m-3 TOTAL soil; C14 in pmol m-3]
    REAL(lp), ALLOCATABLE :: inorganic(:,:,:) !< Group,parent/isotope,layer [mol m-3]
    REAL(lp), ALLOCATABLE :: biology(:,:) !< 21 physiological/lagged states,layer [native units]
    REAL(lp), ALLOCATABLE :: analytical(:,:) !< CUE/poly/AP sampled integrals,layer [s]
    REAL(lp), ALLOCATABLE :: pore(:) !< Independent total pore fraction,layer [m3 m-3]
    REAL(lp), ALLOCATABLE :: capacity(:,:) !< OM/NH4/PO4 specific mineral capacities [mol kg-1 mineral]
    REAL(lp), ALLOCATABLE :: chemistry(:,:) !< Organic PO4 capacity, pH, weathering [native units],layer
  END TYPE
  TYPE :: t_layer_reference
    REAL(lp), POINTER :: p(:) => NULL() !< Non-owning one-cell process field
  END TYPE
CONTAINS
  !> One explicit registry for both directions prevents mismatched read/write lists.
  SUBROUTINE state_references(tile, options, refs)
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options
    TYPE(t_layer_reference), INTENT(out) :: refs(48)
    TYPE(t_sb_memory), POINTER :: soil
    INTEGER :: c,b
    c=options%ics; b=options%iblk
    SELECT TYPE(m => tile%mem(SB_)%p)
    TYPE IS(t_sb_memory)
      soil => m
    END SELECT
    refs(1)%p => soil%nh4_solute%ptr(c,:,b)
    refs(2)%p => soil%nh4_assoc%ptr(c,:,b)
    refs(3)%p => soil%no3_solute%ptr(c,:,b)
    refs(4)%p => soil%noy%ptr(c,:,b)
    refs(5)%p => soil%n2o%ptr(c,:,b)
    refs(6)%p => soil%n2%ptr(c,:,b)
    refs(7)%p => soil%po4_assoc_fast%ptr(c,:,b)
    refs(8)%p => soil%po4_assoc_slow%ptr(c,:,b)
    refs(9)%p => soil%po4_occluded%ptr(c,:,b)
    refs(10)%p => soil%po4_primary%ptr(c,:,b)
    refs(11)%p => soil%po4_solute%ptr(c,:,b)
    refs(12)%p => soil%nh4_n15_solute%ptr(c,:,b)
    refs(13)%p => soil%nh4_n15_assoc%ptr(c,:,b)
    refs(14)%p => soil%no3_n15_solute%ptr(c,:,b)
    refs(15)%p => soil%noy_n15%ptr(c,:,b)
    refs(16)%p => soil%n2o_n15%ptr(c,:,b)
    refs(17)%p => soil%n2_n15%ptr(c,:,b)
    refs(18)%p => soil%microbial_cue_eff%ptr(c,:,b)
    refs(19)%p => soil%microbial_cue_eff_tmic_mavg%ptr(c,:,b)
    refs(20)%p => soil%microbial_nue_eff%ptr(c,:,b)
    refs(21)%p => soil%microbial_nue_eff_tmic_mavg%ptr(c,:,b)
    refs(22)%p => soil%microbial_pue_eff%ptr(c,:,b)
    refs(23)%p => soil%microbial_pue_eff_tmic_mavg%ptr(c,:,b)
    refs(24)%p => soil%fact_n_status_mic_c_growth%ptr(c,:,b)
    refs(25)%p => soil%fact_p_status_mic_c_growth%ptr(c,:,b)
    refs(26)%p => soil%enzyme_frac_AP%ptr(c,:,b)
    refs(27)%p => soil%enzyme_frac_poly%ptr(c,:,b)
    refs(28)%p => soil%enzyme_frac_residue%ptr(c,:,b)
    refs(29)%p => soil%enzyme_frac_poly_c%ptr(c,:,b)
    refs(30)%p => soil%enzyme_frac_poly_n%ptr(c,:,b)
    refs(31)%p => soil%enzyme_frac_poly_p%ptr(c,:,b)
    refs(32)%p => soil%enzyme_frac_poly_c_mavg%ptr(c,:,b)
    refs(33)%p => soil%enzyme_frac_poly_n_mavg%ptr(c,:,b)
    refs(34)%p => soil%enzyme_frac_poly_p_mavg%ptr(c,:,b)
    refs(35)%p => soil%dom_cn%ptr(c,:,b)
    refs(36)%p => soil%dom_cp%ptr(c,:,b)
    refs(37)%p => soil%dom_cn_mavg%ptr(c,:,b)
    refs(38)%p => soil%dom_cp_mavg%ptr(c,:,b)
    refs(39)%p => soil%sb_accel_cue_sum_sl%ptr(c,:,b)
    refs(40)%p => soil%sb_accel_poly_sum_sl%ptr(c,:,b)
    refs(41)%p => soil%sb_accel_ap_sum_sl%ptr(c,:,b)
    refs(42)%p => soil%volume_pore_total_sl%ptr(c,:,b)
    refs(43)%p => soil%qmax_org_min_sl%ptr(c,:,b)
    refs(44)%p => soil%qmax_nh4_min_sl%ptr(c,:,b)
    refs(45)%p => soil%qmax_po4_min_sl%ptr(c,:,b)
    refs(46)%p => soil%qmax_po4_om_sl%ptr(c,:,b)
    refs(47)%p => soil%ph_sl%ptr(c,:,b)
    refs(48)%p => soil%vmax_weath_mineral_sl%ptr(c,:,b)
  END SUBROUTINE

  !> Read actual BGCM matrices, not potentially stale output variables.
  SUBROUTINE read_sb_combined_state(tile,options,state)
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options
    TYPE(t_sb_combined_state), INTENT(out) :: state
    TYPE(t_layer_reference) :: refs(48)
    REAL(lp), POINTER :: pool(:,:,:,:)
    TYPE(t_hydro_memory), POINTER :: hydro
    INTEGER :: k,ns
    pool => tile%bgcm_store%Get_matrix_from_store_2l_3d( &
      & SB_BGCM_POOL_ID,options%ics,options%ics,options%iblk,tile%name)
    SELECT TYPE(m => tile%mem(HYDRO_)%p)
    TYPE IS(t_hydro_memory)
      hydro => m
    END SELECT
    ns=NINT(hydro%num_sl_above_bedrock%ptr(options%ics,options%iblk))
    state%organic=pool(:,:,1,1:ns)
    CALL state_references(tile,options,refs)
    ALLOCATE(state%inorganic(11,2,ns),SOURCE=0._lp)
    ALLOCATE(state%biology(21,ns),state%analytical(3,ns),state%pore(ns), &
      & state%capacity(3,ns),state%chemistry(3,ns))
    DO k=1,11
      state%inorganic(k,1,:)=refs(k)%p(1:ns)
    END DO
    DO k=1,6
      state%inorganic(k,2,:)=refs(k+11)%p(1:ns)
    END DO
    DO k=1,21
      state%biology(k,:)=refs(k+17)%p(1:ns)
    END DO
    DO k=1,3
      state%analytical(k,:)=refs(k+38)%p(1:ns)
      state%capacity(k,:)=refs(k+42)%p(1:ns)
    END DO
    state%pore=refs(42)%p(1:ns)
    state%chemistry(1,:)=refs(46)%p(1:ns)
    state%chemistry(2,:)=refs(47)%p(1:ns)
    state%chemistry(3,:)=refs(48)%p(1:ns)
  END SUBROUTINE

  !> Write only an explicitly staged state; caller owns checkpoint and validation.
  !> No synchronization of unrelated matrices, no tile fraction or clock changes.
  SUBROUTINE write_sb_combined_state(tile,options,state)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options
    TYPE(t_sb_combined_state), INTENT(in) :: state
    TYPE(t_layer_reference) :: refs(48)
    REAL(lp), POINTER :: pool(:,:,:,:)
    INTEGER :: k,ns
    pool => tile%bgcm_store%Get_matrix_from_store_2l_3d( &
      & SB_BGCM_POOL_ID,options%ics,options%ics,options%iblk,tile%name)
    ns=SIZE(state%organic,3)
    pool(:,:,1,1:ns)=state%organic
    CALL state_references(tile,options,refs)
    DO k=1,11
      refs(k)%p(1:ns)=state%inorganic(k,1,:)
    END DO
    DO k=1,6
      refs(k+11)%p(1:ns)=state%inorganic(k,2,:)
    END DO
    DO k=1,21
      refs(k+17)%p(1:ns)=state%biology(k,:)
    END DO
    DO k=1,3
      refs(k+38)%p(1:ns)=state%analytical(k,:)
      refs(k+42)%p(1:ns)=state%capacity(k,:)
    END DO
    refs(42)%p(1:ns)=state%pore
    refs(46)%p(1:ns)=state%chemistry(1,:)
    refs(47)%p(1:ns)=state%chemistry(2,:)
    refs(48)%p(1:ns)=state%chemistry(3,:)
  END SUBROUTINE

  !> Aggregate inventories by area, biology by microbial biomass, and capacities
  !> by mineral mass. Derived physics is deliberately absent from this inventory.
  SUBROUTINE reconstruct_sb_combined_state(fraction,depth,before,shared,accepted,reason,failed_member)
    REAL(lp), INTENT(in) :: fraction(:),depth(:,:) !< Tile area [-], layer/tile thickness [m]
    TYPE(t_sb_combined_state), INTENT(in) :: before(:)
    TYPE(t_sb_combined_state), INTENT(out) :: shared
    LOGICAL, INTENT(out) :: accepted
    CHARACTER(len=*), OPTIONAL, INTENT(out) :: reason !< First failed reconstruction check
    INTEGER, OPTIONAL, INTENT(out) :: failed_member !< Input PFT index; zero denotes shared state
    CHARACTER(len=256) :: detail !< Inventory failure with field indices and native value
    REAL(lp), ALLOCATABLE :: organic(:,:,:,:), inorg(:,:,:,:), pore(:,:), capacity(:,:,:)
    REAL(lp), ALLOCATABLE :: biology(:,:,:), microbial(:,:), weights(:)
    INTEGER :: nt,ns,np,ne,j
    accepted=.FALSE.
    IF (PRESENT(reason)) reason='Invalid participating tile count or fraction-array size.'
    IF (PRESENT(failed_member)) failed_member=0
    nt=SIZE(before)
    IF (nt == 0 .OR. SIZE(fraction) /= nt) RETURN
    shared=before(1)
    np=SIZE(shared%organic,1); ne=SIZE(shared%organic,2); ns=SIZE(shared%organic,3)
    ALLOCATE(organic(np,ne,ns,nt),inorg(11,2,ns,nt),pore(ns,nt),capacity(3,ns,nt), &
      & biology(21,ns,nt),microbial(ns,nt),weights(nt))
    DO j=1,nt
      IF (.NOT. valid_sb_combined_inventory(before(j),detail)) THEN
        IF (PRESENT(reason)) reason=detail
        IF (PRESENT(failed_member)) failed_member=j
        RETURN
      END IF
      organic(:,:,:,j)=before(j)%organic
      inorg(:,:,:,j)=before(j)%inorganic
      pore(:,j)=before(j)%pore
      capacity(:,:,j)=before(j)%capacity
      biology(:,:,j)=before(j)%biology
      microbial(:,j)=before(j)%organic(ix_microbial,ixC,:)
    END DO
    IF (PRESENT(reason)) reason='Organic aggregation: invalid area weights, layer geometry or inventory.'
    CALL combine_sb_columns(fraction,depth,organic,weights,shared%organic,accepted)
    IF (.NOT. accepted) RETURN
    IF (PRESENT(reason)) reason='Inorganic aggregation: invalid area weights, layer geometry or inventory.'
    CALL combine_sb_columns(fraction,depth,inorg,weights,shared%inorganic,accepted)
    IF (.NOT. accepted) RETURN
    IF (PRESENT(reason)) reason='Structure aggregation: invalid pore volume or mineral-weighted capacity.'
    CALL combine_sb_structure(weights,pore,capacity,shared%pore,shared%capacity,accepted)
    IF (.NOT. accepted) RETURN
    IF (PRESENT(reason)) reason='Physiology aggregation: invalid microbial weights or physiological state.'
    CALL combine_sb_physiology(weights,microbial,biology,shared%biology,accepted)
    IF (.NOT. accepted) RETURN
    shared%analytical=0._lp
    shared%chemistry=0._lp
    DO j=1,nt
      shared%analytical=shared%analytical+weights(j)*before(j)%analytical
      shared%chemistry=shared%chemistry+weights(j)*before(j)%chemistry
    END DO
    accepted=valid_sb_combined_inventory(shared,detail)
    IF (PRESENT(reason)) reason=detail
  END SUBROUTINE

  !> Propose all PFT states without mutation. Caller must recompute forward
  !> physics and check capacity/volume admissibility before accepting ANY PFT.
  SUBROUTINE redistribute_sb_combined_state(fraction,depth,before,target,lower,upper,after,accepted)
    REAL(lp), INTENT(in) :: fraction(:),depth(:,:),lower(:,:),upper(:,:)
      !< Areas [-], thickness [m], existing per-PFT pore bounds [-]
    TYPE(t_sb_combined_state), INTENT(in) :: before(:),target
    TYPE(t_sb_combined_state), ALLOCATABLE, INTENT(out) :: after(:)
    LOGICAL, INTENT(out) :: accepted
    REAL(lp), ALLOCATABLE :: old(:,:,:,:), new(:,:,:,:), pore(:,:), result(:,:)
    REAL(lp), ALLOCATABLE :: old_n(:,:,:,:),new_n(:,:,:,:)
    REAL(lp), ALLOCATABLE :: mean(:,:,:),mean_n(:,:,:),weights(:)
    INTEGER :: nt,ns,np,ne,j
    accepted=.FALSE.
    IF (.NOT. valid_sb_combined_inventory(target)) RETURN
    nt=SIZE(before); ns=SIZE(target%organic,3)
    np=SIZE(target%organic,1); ne=SIZE(target%organic,2)
    ALLOCATE(old(np,ne,ns,nt),new(np,ne,ns,nt),old_n(11,2,ns,nt),new_n(11,2,ns,nt))
    ALLOCATE(pore(ns,nt),result(ns,nt))
    ALLOCATE(mean(np,ne,ns),mean_n(11,2,ns),weights(nt))
    DO j=1,nt
      old(:,:,:,j)=before(j)%organic
      old_n(:,:,:,j)=before(j)%inorganic
      pore(:,j)=before(j)%pore
    END DO
    CALL redistribute_sb_columns(fraction,depth,old,target%organic,ixC,new,accepted)
    IF (.NOT. accepted) RETURN
    CALL redistribute_sb_columns(fraction,depth,old_n,target%inorganic,1,new_n,accepted)
    IF (.NOT. accepted) RETURN
    CALL redistribute_sb_structure(fraction/SUM(fraction),pore,target%pore,lower,upper,result,accepted)
    IF (.NOT. accepted) RETURN
    after=before
    DO j=1,nt
      after(j)%organic=new(:,:,:,j)
      after(j)%inorganic=new_n(:,:,:,j)
      after(j)%pore=result(:,j)
      ! Shared physiological condition, not biomass-scaled efficiencies/ratios.
      after(j)%biology=target%biology
      ! Restore individual forcing histories and specific mineral properties.
      ! Their post-event reset is performed only by the scheduler on commit.
      IF (.NOT. valid_sb_combined_inventory(after(j))) THEN
        accepted=.FALSE.
        RETURN
      END IF
    END DO
    ! Every pool, layer and isotope must close, not merely the five aggregated
    ! convergence diagnostics. Reuse the approved normalized roundoff allowance.
    CALL combine_sb_columns(fraction,depth,new,weights,mean,accepted)
    IF (.NOT. accepted) RETURN
    IF (ANY(ABS(mean-target%organic) > 1.e-12_lp*MAX(ABS(target%organic),TINY(1._lp)))) THEN
      accepted=.FALSE.
      RETURN
    END IF
    CALL combine_sb_columns(fraction,depth,new_n,weights,mean_n,accepted)
    IF (.NOT. accepted) RETURN
    IF (ANY(ABS(mean_n-target%inorganic) > 1.e-12_lp*MAX(ABS(target%inorganic),TINY(1._lp)))) &
      & accepted=.FALSE.
  END SUBROUTINE

  !> Inventory-only checks. Capacity and occupancy require forward refresh.
  LOGICAL FUNCTION valid_sb_combined_inventory(state,reason) RESULT(valid)
    TYPE(t_sb_combined_state), INTENT(in) :: state
    CHARACTER(len=*), OPTIONAL, INTENT(out) :: reason !< First failure, indices and value/bound
    INTEGER :: s,p !< Layer and pool/group indices
    REAL(lp) :: isotope_c !< Sum of carbon isotopes in parent-carbon units [mol m-3]
    valid=.FALSE.
    IF (PRESENT(reason)) reason=''
    ! Check finiteness before arithmetic, including when FP exception trapping is enabled.
    DO s=1,SIZE(state%organic,3)
      IF (.NOT. check_fields('organic',state%organic(:,:,s),s,.TRUE.)) RETURN
    END DO
    DO s=1,SIZE(state%inorganic,3)
      IF (.NOT. check_fields('inorganic',state%inorganic(:,:,s),s,.TRUE.)) RETURN
    END DO
    DO s=1,SIZE(state%pore)
      IF (.NOT. check_fields('biology',state%biology(:,s:s),s,.FALSE.)) RETURN
      IF (.NOT. check_fields('capacity',state%capacity(:,s:s),s,.TRUE.)) RETURN
      IF (.NOT. check_fields('chemistry',state%chemistry(:,s:s),s,.FALSE.)) RETURN
      IF (.NOT. check_fields('analytical',state%analytical(:,s:s),s,.FALSE.)) RETURN
      IF (.NOT. ieee_is_finite(state%pore(s))) THEN
        CALL reject('pore nonfinite',s,0,0,state%pore(s),1._lp)
        RETURN
      END IF
      IF (state%pore(s) < 0._lp .OR. state%pore(s) >= 1._lp) THEN
        CALL reject('pore outside [0,1)',s,0,0,state%pore(s),1._lp)
        RETURN
      END IF
    END DO
    DO s=1,SIZE(state%inorganic,3)
      DO p=1,SIZE(state%inorganic,1)
        IF (state%inorganic(p,2,s) > state%inorganic(p,1,s)) THEN
          CALL reject('inorganic isotope > parent',s,p,2,state%inorganic(p,2,s),state%inorganic(p,1,s))
          RETURN
        END IF
      END DO
    END DO
    DO s=1,SIZE(state%organic,3)
      DO p=1,SIZE(state%organic,1)
        ! C and C13 are mol; C14 is stored as pmol for numerical precision.
        ! Convert only for validation, never rescale the persistent tracer itself.
        isotope_c=state%organic(p,ixC13,s)+state%organic(p,ixC14,s)/scale_correction_c14
        IF (isotope_c > state%organic(p,ixC,s)) THEN
          CALL reject('C13 + C14/scale > C [mol m-3]',s,p,ixC,isotope_c,state%organic(p,ixC,s))
          RETURN
        END IF
        IF (state%organic(p,ixN15,s) > state%organic(p,ixN,s)) THEN
          CALL reject('organic N15 > N [mol m-3]',s,p,ixN15,state%organic(p,ixN15,s),state%organic(p,ixN,s))
          RETURN
        END IF
      END DO
    END DO
    valid=.TRUE.
  CONTAINS
    !> First offending native-unit field; rank-2 physiological fields use component=1.
    LOGICAL FUNCTION check_fields(label,values,layer,nonnegative) RESULT(ok)
      CHARACTER(len=*), INTENT(in) :: label
      REAL(lp), INTENT(in) :: values(:,:)
      INTEGER, INTENT(in) :: layer
      LOGICAL, INTENT(in) :: nonnegative
      INTEGER :: i,j
      ok=.FALSE.
      DO j=1,SIZE(values,2)
        DO i=1,SIZE(values,1)
          IF (.NOT. ieee_is_finite(values(i,j))) THEN
            CALL reject(label//' nonfinite',layer,i,j,values(i,j),0._lp)
            RETURN
          END IF
          IF (nonnegative .AND. values(i,j) < 0._lp) THEN
            CALL reject(label//' negative',layer,i,j,values(i,j),0._lp)
            RETURN
          END IF
        END DO
      END DO
      ok=.TRUE.
    END FUNCTION check_fields

    SUBROUTINE reject(label,layer,field,component,value,bound)
      CHARACTER(len=*), INTENT(in) :: label
      INTEGER, INTENT(in) :: layer,field,component
      REAL(lp), INTENT(in) :: value,bound
      IF (PRESENT(reason)) WRITE(reason,'(A,A,3(I0,1X),A,2(ES24.16E3,1X))') &
        & TRIM(label),'; layer/field/component=',layer,field,component,'; value/bound=',value,bound
    END SUBROUTINE reject
  END FUNCTION
#endif
END MODULE mo_sb_combined_state
