!> Cell-local transactional backups for soil replay and combined-column staging
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
MODULE mo_sb_replay_state
#ifndef __NO_JSBACH__
  USE mo_kind, ONLY: lp
  USE mo_exception, ONLY: finish
  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_tile_class, ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class, ONLY: t_jsb_task_options
  USE mo_jsb_memory_class, ONLY: t_jsb_memory_p
  USE mo_jsb_var_class, ONLY: REAL1D, REAL2D, REAL3D
  USE mo_jsb_process_class, ONLY: SB_, VEG_, HYDRO_, SSE_, A2L_
  USE mo_lnd_bgcm_store_class, ONLY: t_lnd_bgcm_store_abstract
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_sb_replay_checkpoint, capture_sb_replay_state, restore_sb_replay_state
  PUBLIC :: release_sb_replay_state

  ! A backup owns VALUES, not model memory. Pointers are restore destinations only;
  ! they must not be deallocated or rebound. All horizontal arrays are sliced to
  ! the current cell before saving; no whole-domain forcing/state copy is made.
  TYPE :: t_replay_field
    REAL(lp), POINTER :: p1(:) => NULL(), p2(:,:) => NULL()
    REAL(lp), POINTER :: p3(:,:,:) => NULL(), p4(:,:,:,:) => NULL()
    REAL(lp), ALLOCATABLE :: v1(:), v2(:,:), v3(:,:,:), v4(:,:,:,:)
      !< Saved values, units identical to the original registered variable/matrix
  END TYPE
  TYPE :: t_store_reference
    CLASS(t_lnd_bgcm_store_abstract), POINTER :: p => NULL() !< Non-owning store reference
  END TYPE
  TYPE :: t_sb_replay_checkpoint
    PRIVATE
    TYPE(t_replay_field), ALLOCATABLE :: fields(:) !< Owned value snapshots
    INTEGER :: count = 0 !< Number of populated field snapshots [-]
  END TYPE
CONTAINS
  !> Capture all touched process owners BEFORE changing any participating tile.
  !! Process registries include restart state, lagged biology, forcing sums and
  !! non-restart diagnostics. BGCM matrix arrays are saved separately: registered
  !! variables and matrices need not be synchronized at this entry point, and
  !! rollback must restore BOTH representations exactly, without imposing a sync.
  !! The caller supplies pending positive-area SB tiles. Inherited HYDRO/SSE/A2L
  !! memory and shared BGCM stores are captured once by pointer identity, not name.
  !! No calendar, tile fraction, task schedule, allocation or registration changes
  !! are allowed inside the transaction. This is a CPU, one-cell driver facility;
  !! it must not be called concurrently for aliases of the same process memory.
  SUBROUTINE capture_sb_replay_state(model, options, selected, checkpoint)
    TYPE(t_jsb_model), INTENT(in) :: model
    TYPE(t_jsb_task_options), INTENT(in) :: options
    LOGICAL, INTENT(in) :: selected(:) !< Participating model%tiles entries [-]
    TYPE(t_sb_replay_checkpoint), INTENT(inout) :: checkpoint
    TYPE(t_jsb_memory_p), ALLOCATABLE :: memories(:)
    TYPE(t_store_reference), ALLOCATABLE :: stores(:)
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    INTEGER, PARAMETER :: processes(5) = [SB_, VEG_, HYDRO_, SSE_, A2L_]
    INTEGER :: it, ip, proc, j, nm, ns, capacity, k, v, c, b
    LOGICAL :: found

    IF (ALLOCATED(checkpoint%fields)) CALL finish('capture_sb_replay_state', 'Checkpoint already active.')
    IF (SIZE(selected) /= SIZE(model%tiles)) CALL finish('capture_sb_replay_state', 'Invalid tile mask.')
    IF (options%ics /= options%ice) CALL finish('capture_sb_replay_state', 'Expected one owned grid cell.')
    c = options%ics
    b = options%iblk
    ALLOCATE(memories(SIZE(selected)*SIZE(processes)), stores(SIZE(selected)))
    nm = 0
    ns = 0
    DO it = 1, SIZE(selected)
      IF (.NOT. selected(it)) CYCLE
      tile => model%tiles(it)%p
      DO ip = 1, SIZE(processes)
        proc = processes(ip)
        IF (.NOT. ASSOCIATED(tile%mem(proc)%p)) CYCLE
        found = .FALSE.
        DO j = 1, nm
          IF (ASSOCIATED(memories(j)%p, tile%mem(proc)%p)) found = .TRUE.
        END DO
        IF (found) CYCLE
        nm = nm + 1
        memories(nm)%p => tile%mem(proc)%p
      END DO
      IF (.NOT. ASSOCIATED(tile%bgcm_store)) CYCLE
      found = .FALSE.
      DO j = 1, ns
        IF (ASSOCIATED(stores(j)%p, tile%bgcm_store)) found = .TRUE.
      END DO
      IF (.NOT. found) THEN
        ns = ns + 1
        stores(ns)%p => tile%bgcm_store
      END IF
    END DO
    capacity = 0
    DO j = 1, nm
      capacity = capacity + memories(j)%p%no_of_vars
    END DO
    DO j = 1, ns
      capacity = capacity + stores(j)%p%nr_of_bgc_materials
    END DO
    ALLOCATE(checkpoint%fields(capacity))
    k = 0
    DO j = 1, nm
      DO v = 1, memories(j)%p%no_of_vars
        ASSOCIATE(var => memories(j)%p%vars(v)%p)
          k = k + 1
          SELECT CASE(var%type)
          CASE(REAL1D)
            IF (.NOT. ASSOCIATED(var%ptr1d)) CYCLE
            checkpoint%fields(k)%p1 => var%ptr1d
            checkpoint%fields(k)%v1 = var%ptr1d
          CASE(REAL2D)
            IF (.NOT. ASSOCIATED(var%ptr2d)) CYCLE
            checkpoint%fields(k)%p1 => var%ptr2d(c:c,b)
            checkpoint%fields(k)%v1 = var%ptr2d(c:c,b)
          CASE(REAL3D)
            IF (.NOT. ASSOCIATED(var%ptr3d)) CYCLE
            checkpoint%fields(k)%p2 => var%ptr3d(c:c,:,b)
            checkpoint%fields(k)%v2 = var%ptr3d(c:c,:,b)
          CASE DEFAULT
            CALL finish('capture_sb_replay_state', 'Unsupported variable rank in replay process memory.')
          END SELECT
        END ASSOCIATE
      END DO
    END DO
    DO j = 1, ns
      DO v = 1, stores(j)%p%nr_of_1l_2d_bgc_materials
        k = k + 1
        checkpoint%fields(k)%p2 => stores(j)%p%store_1l_2d_bgcms(v)%mt_1l_2d_bgcm(:,c:c,b)
        checkpoint%fields(k)%v2 = checkpoint%fields(k)%p2
      END DO
      DO v = 1, stores(j)%p%nr_of_1l_3d_bgc_materials
        k = k + 1
        checkpoint%fields(k)%p3 => stores(j)%p%store_1l_3d_bgcms(v)%mt_1l_3d_bgcm(:,c:c,:,b)
        checkpoint%fields(k)%v3 = checkpoint%fields(k)%p3
      END DO
      DO v = 1, stores(j)%p%nr_of_2l_2d_bgc_materials
        k = k + 1
        checkpoint%fields(k)%p3 => stores(j)%p%store_2l_2d_bgcms(v)%mt_2l_2d_bgcm(:,:,c:c,b)
        checkpoint%fields(k)%v3 = checkpoint%fields(k)%p3
      END DO
      DO v = 1, stores(j)%p%nr_of_2l_3d_bgc_materials
        k = k + 1
        checkpoint%fields(k)%p4 => stores(j)%p%store_2l_3d_bgcms(v)%mt_2l_3d_bgcm(:,:,c:c,:,b)
        checkpoint%fields(k)%v4 = checkpoint%fields(k)%p4
      END DO
    END DO
    checkpoint%count = k
  END SUBROUTINE capture_sb_replay_state

  !> Restore pre-event values without resetting or synchronizing either BGCM view.
  !! A normal rejected candidate can then be discarded transactionally. This is
  !! not recovery from finish(), MPI abort, or an arithmetic trap inside a kernel.
  SUBROUTINE restore_sb_replay_state(checkpoint)
    TYPE(t_sb_replay_checkpoint), INTENT(in) :: checkpoint
    INTEGER :: k
    IF (.NOT. ALLOCATED(checkpoint%fields)) CALL finish('restore_sb_replay_state', 'No active checkpoint.')
    DO k = 1, checkpoint%count
      ASSOCIATE(f => checkpoint%fields(k))
        IF (ALLOCATED(f%v1)) f%p1 = f%v1
        IF (ALLOCATED(f%v2)) f%p2 = f%v2
        IF (ALLOCATED(f%v3)) f%p3 = f%v3
        IF (ALLOCATED(f%v4)) f%p4 = f%v4
      END ASSOCIATE
    END DO
  END SUBROUTINE restore_sb_replay_state

  !> Release owned scratch values only; leave model targets and their pointers intact.
  SUBROUTINE release_sb_replay_state(checkpoint)
    TYPE(t_sb_replay_checkpoint), INTENT(inout) :: checkpoint
    IF (ALLOCATED(checkpoint%fields)) DEALLOCATE(checkpoint%fields)
    checkpoint%count = 0
  END SUBROUTINE release_sb_replay_state
#endif
END MODULE mo_sb_replay_state
