!> Methods to deal with conserved quantities in land cover change (lcc) processes with active or passive relocations
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
!>#### Provides methods to deal with active or passive relocation of conserved matter in land cover change (lcc) processes
!>
!> Each lcc process that actively (from one variable of a shrinking tile to another variable of an increasing tile)
!> or passively (from a variable on a shrinking tile to the same variable on another tile) relocates any conserved quantity (CQ)
!> needs an own lcc structure ([[t_jsb_lcc_proc]]).
!> Which is located on the tile for which the process is calculated, i.e. on a parent or grandparent of the
!> involved tiles in between which matter is relocated.
!> This module provides methods to cerate and modify an instance of such an lcc structure.
!>
!> NOTE: This is currently only implemented to deal with leaves and down to two layers,
!>       if something more generic is required this needs to be re-engineered.

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_jsb_lcc
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_util,                ONLY: int2string, one_of, tolower
  USE mo_kind,                ONLY: lp

  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_varlist,         ONLY: VARNAME_LEN
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN
  USE mo_jsb_var_class,       ONLY: REAL2D, REAL3D

  USE mo_jsb_lcc_class,       ONLY: t_jsb_lcc_proc, t_jsb_lcc_var_p, max_tolerated_cf_mismatch
  USE mo_jsb_task_class,      ONLY: t_jsb_task_options
  USE mo_jsb_cqt_class,       ONLY: t_jsb_consQuan_p
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract

  USE mo_jsb_process_class,   ONLY: Get_process_id, Get_process_name

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: end_lcc_reloc, start_lcc_reloc, init_lcc_reloc, init_lcc, transfer_active_to_passive_onChunk
  PUBLIC :: count_descendants_for_up_to_2layers, collect_names_of_descendant_leaves_up_to_2layers

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_lcc'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initializes the lcc structure instance for the given process on the given tile on which the process is calculated.
  !> The routine is called by the "init lcc" routines of the lcc processes that need to relocate conserved quantities,
  !> e.g. from [[mo_alcc_init_lcc:alcc_init_lcc]] of the ALCC (anthropogenic land cover change) process.
  !
  SUBROUTINE init_lcc(process_name, tile, active_cqts, passive_cqts, involved_descendant_tiles)

    USE mo_jsb_class,         ONLY: Get_model
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid
    USE mo_jsb_grid,          ONLY: Get_grid
    USE mo_carbon_constants,  ONLY: carbon_potential_active_vars, carbon_required_passive_vars

    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*),  INTENT(in) :: process_name                 !< Name of the encompassing lcc process
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile             !< Tile on which the lcc process is calculated
    INTEGER,           INTENT(in) :: active_cqts(:)               !< CQTs that are actively relocated by this lcc process
    INTEGER,           INTENT(in) :: passive_cqts(:)              !< CQTs that are passively relocated by this lcc process
    CHARACTER(len=*),  INTENT(in) :: involved_descendant_tiles(:) !< Descendants tiles involved in this lcc process
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_lcc' !< Name of this routine

    TYPE(t_jsb_model), POINTER          :: model                  !< Current instance of the model
    TYPE(t_jsb_grid),  POINTER          :: hgrid                  !< Pointer to a horizontal grid
    TYPE(t_jsb_lcc_proc), POINTER       :: lcc_relocations        !< The new instance of the lcc relocation (lcc structure)
    CLASS(t_jsb_tile_abstract), POINTER :: child, &               !< Pointer to the current child of the main tile
      &                                    leave                  !< Pointer to the current leave: either the current child or
                                                                  !< a child of the current child - depending on the tile tree

    INTEGER :: nr_of_children, &                       !< Number of children in the tile tree below the main tile
      &        nr_of_leaves, &                         !< Number of leaves in the tile tree below a child of the main tile
      &        nr_of_active_vars, &                    !< Number of active variables
      &        nr_of_passive_vars, &                   !< Number of passive variables
      &        nr_of_tiles, &                          !< Number of tiles involved with this lcc process
      &        nr_of_vars, &                           !< Total number of to be conserved quantities
      &        nr_of_procs_with_active_vars, &         !< Number of processes with active variables
      &        nr_of_active_vars_this_proc, &          !< Number of active variables on a certain process
      &        nr_of_passive_vars_this_proc, &         !< Number of passive variables on a certain process
      &        lcc_process_id, &                       !< ID of this lcc process
      &        this_process_id, &                      !< ID of the current process with a conserved quantity
      &        ind, &                                  !< Generic helper variable used for index locations in different arrays
      &        i_child, &                              !< TODO: not required can be replaced
      &        i_leave, &                              !< TODO: not required can be replaced
      &        i_cqt, &                                !< Index of the current CQT
      &        i_var, &                                !< Index of the current variable
      &        i_active, &                             !< Index of an active variable
      &        i_active_this_tile, &                   !< TODO: could probably simply reuse i_active since different loop
      &        i_passive, &                            !< Index of a passive variable
      &        i_passive_this_tile, &                  !< TODO: could probably simply reuse i_passive since different loop
      &        i_tile, &                               !< Tile index
      &        i_proc, &                               !< Index of the current active process
      &        nproma, &                               !< Block length
      &        nblks                                   !< Number of blocks
    LOGICAL :: is_involved_tile                        !< Flag indicating if a child tile is involved with the lcc process
    LOGICAL, ALLOCATABLE :: is_descendant(:)           !< Check if all involved tiles are descendants of the main tile
    INTEGER, ALLOCATABLE :: active_var_process_ids(:)  !< Array of process IDs for processes holding active variables
    CHARACTER(len=:), ALLOCATABLE :: this_proc_name    !< Name of the current process
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: &
      &        active_vars_this_proc(:), &             !< Names of active variables on this process
      &        passive_vars_this_proc(:)               !< Names of passive variables on this process

    ! -------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    hgrid => Get_grid(model%grid_id)
    lcc_process_id = Get_process_id(process_name)

    ! Allocate the new lcc structure for this process on this tile
    ALLOCATE(t_jsb_lcc_proc::lcc_relocations)
    tile%lcc_processes(lcc_process_id)%p => lcc_relocations

    IF (debug_on()) CALL message(TRIM(routine), 'For proc '//TRIM(process_name)//' ...')

    nproma = hgrid%nproma
    nblks  = hgrid%nblks

    nr_of_tiles = SIZE(involved_descendant_tiles(:))

    lcc_relocations%nr_of_tiles = nr_of_tiles
    ALLOCATE(is_descendant(nr_of_tiles))
    is_descendant = .FALSE.

    lcc_relocations%name = TRIM(process_name) // "_lcc"
    lcc_relocations%tile_names = involved_descendant_tiles(:)

    nr_of_active_vars = 0
    nr_of_passive_vars = 0
    nr_of_procs_with_active_vars = 0

    ! Get the number of descendants of the tile in order to iterate over these and to collect the CQ variables
    ! of all CQTs which are actively or passively relocated by the calling lcc process
    IF (tile%Has_children()) THEN
      nr_of_children = tile%Get_no_of_children()
      child => tile%Get_first_child_tile()
    ELSE
      nr_of_children = 1
      child => tile%self_tile
    END IF

    i_tile = 0
    DO i_child = 1,nr_of_children !@TODO: no do loop with counter required: do while!

      ! Check if this is a leave (currently lcc can only run on leaves!) and if not go to the children of this node
      IF (child%Has_children()) THEN
        leave => child%Get_first_child_tile()
        nr_of_leaves = child%Get_no_of_children()
      ELSE
        leave => child
        nr_of_leaves = 1
      END IF

      DO i_leave = 1,nr_of_leaves !@TODO: no do loop with counter required: do while!

        ! Only descendant tiles listed as involved are used
        is_involved_tile = .FALSE.
        ind = one_of(leave%name, involved_descendant_tiles)
        IF (ind > 0) THEN
          i_tile = i_tile + 1
          is_involved_tile = .TRUE.
          is_descendant(ind) = .TRUE.
        END IF

        ! The lcc structure is created according to active and passive cqs on the first involved tile
        IF (is_involved_tile .AND. (i_tile == 1)) THEN

          ! In the first step only collect the number of active and passive variables
          DO i_cqt = 1, leave%nr_of_cqts
            IF (ANY(leave%conserved_quantities(i_cqt)%p%type_id == active_cqts(:) )) THEN
              nr_of_active_vars = nr_of_active_vars &
                & + leave%conserved_quantities(i_cqt)%p%no_of_vars
            END IF
            IF (ANY(leave%conserved_quantities(i_cqt)%p%type_id == passive_cqts(:) )) THEN
              nr_of_passive_vars = nr_of_passive_vars &
                & + leave%conserved_quantities(i_cqt)%p%no_of_vars
            END IF

            ! Note: For FLCC (used in JSBACH configurations running forest age classes (FAGE)) a CQT can be active and passive
            !       at the same time: ag woody litter is burning upon fire, but at the same time filled with new woody litter
            !       which then needs to be passively allocated, therefore two separate ifs and no assertion on disjunct cqt groups

          END DO

          lcc_relocations%nr_of_active_vars = nr_of_active_vars
          lcc_relocations%nr_of_passive_vars = nr_of_passive_vars

          ! Allocate all collector arrays on the lcc structure
          ALLOCATE(lcc_relocations%active_vars(nr_of_active_vars))
          ALLOCATE(lcc_relocations%active_vars_names(nr_of_active_vars))
          ALLOCATE(lcc_relocations%active_vars_cqt(nr_of_active_vars))
          ALLOCATE(lcc_relocations%active_vars_process_id(nr_of_active_vars))
          ALLOCATE(lcc_relocations%passive_vars(nr_of_passive_vars))
          ALLOCATE(lcc_relocations%passive_vars_names(nr_of_passive_vars))
          ALLOCATE(lcc_relocations%passive_vars_cqt(nr_of_passive_vars))
          ALLOCATE(lcc_relocations%passive_vars_process_id(nr_of_passive_vars))

          ALLOCATE(lcc_relocations%initial_area(nproma, nr_of_tiles, nblks))
          ALLOCATE(lcc_relocations%gained_area(nproma, nr_of_tiles, nblks))
          ALLOCATE(lcc_relocations%lost_area(nproma, nr_of_tiles, nblks))

          ALLOCATE(active_var_process_ids(nr_of_active_vars))
          active_var_process_ids = 0

          ! ... and then collect vars and pointers for the first involved tile
          i_active = 0
          i_passive = 0
          DO i_cqt = 1, leave%nr_of_cqts
            ! For the active cqts
            IF (ANY(leave%conserved_quantities(i_cqt)%p%type_id == active_cqts(:) )) THEN
              nr_of_vars = leave%conserved_quantities(i_cqt)%p%no_of_vars
              DO i_var = 1,nr_of_vars
                i_active = i_active + 1

                this_process_id = leave%conserved_quantities(i_cqt)%p%associated_process(i_var)
                IF (.NOT. ANY(active_var_process_ids == this_process_id)) THEN
                  nr_of_procs_with_active_vars = nr_of_procs_with_active_vars + 1
                  active_var_process_ids(nr_of_procs_with_active_vars) = this_process_id
                END IF

                CALL add_cq_to_lcc_structure(nr_of_tiles, nproma, nblks, i_tile, i_cqt, i_var, i_active,                 &
                  & leave%conserved_quantities, lcc_relocations%active_vars_cqt, lcc_relocations%active_vars_process_id, &
                  & lcc_relocations%active_vars, lcc_relocations%active_vars_names)

              END DO
            END IF
            ! ... and separately for the passive cqts
            IF (ANY(leave%conserved_quantities(i_cqt)%p%type_id == passive_cqts(:) )) THEN
              nr_of_vars = leave%conserved_quantities(i_cqt)%p%no_of_vars
              DO i_var = 1,nr_of_vars
                i_passive = i_passive + 1

                CALL add_cq_to_lcc_structure(nr_of_tiles, nproma, nblks, i_tile, i_cqt, i_var, i_passive,                   &
                  & leave%conserved_quantities, lcc_relocations%passive_vars_cqt,  lcc_relocations%passive_vars_process_id, &
                  & lcc_relocations%passive_vars, lcc_relocations%passive_vars_names )

              END DO
            END IF !... for active or passive CQTs

          END DO !... for all CQTs

          ! Collect the names of the processes carrying variables that are actively relocated by this lcc process
          lcc_relocations%nr_of_procs_with_active_vars = nr_of_procs_with_active_vars
          ALLOCATE(lcc_relocations%unique_proc_names_active_vars(nr_of_procs_with_active_vars))
          i_proc = 0
          DO ind = 1,nr_of_active_vars
            IF (active_var_process_ids(ind) /= 0) THEN
              i_proc = i_proc + 1
              lcc_relocations%unique_proc_names_active_vars(i_proc) &
                & = Get_process_name(active_var_process_ids(ind))
            ELSE
              EXIT
            END IF
          END DO
        END IF !... for the first involved tile

        ! For all other involved tiles (i.e. other than the first) only the pointers need to be collected
        ! + assert that the demanded active and passive variables are part of all involved tiles
        IF (is_involved_tile .AND. (i_tile > 1)) THEN

          i_active_this_tile = 0
          i_passive_this_tile = 0

          DO i_cqt = 1, leave%nr_of_cqts
            IF (ANY(leave%conserved_quantities(i_cqt)%p%type_id == active_cqts(:) )) THEN
              CALL collect_pointers_for_cqs(i_cqt, i_tile, i_active_this_tile, leave%name, &
                & TRIM(involved_descendant_tiles(1)), leave%conserved_quantities,          &
                & lcc_relocations%active_vars, lcc_relocations%nr_of_active_vars, 'active')
            END IF
            IF (ANY(leave%conserved_quantities(i_cqt)%p%type_id == passive_cqts(:) )) THEN
              CALL collect_pointers_for_cqs(i_cqt, i_tile, i_passive_this_tile, leave%name, &
                & TRIM(involved_descendant_tiles(1)), leave%conserved_quantities,           &
                & lcc_relocations%passive_vars, lcc_relocations%nr_of_passive_vars, 'passive')
            END IF
          END DO

          !>
          !> Assertion: all tiles involved in an lcc process need to have the same number of conserved quantities.
          !>
          IF (.NOT. i_active_this_tile == nr_of_active_vars ) THEN
            CALL finish(TRIM(routine), &
              & 'Violation of assertion: Not all active variables found on first tile '//TRIM(involved_descendant_tiles(1)) &
              & //' were also active cqs of this tile '//TRIM(involved_descendant_tiles(i_tile))//'. Found '                &
              & //int2string(nr_of_active_vars)//' on first and '//int2string(i_active_this_tile)//' on this. Please check.')
          END IF
          IF (.NOT. i_passive_this_tile == nr_of_passive_vars ) THEN
            CALL finish(TRIM(routine), &
              & 'Violation of assertion: Not all passive variables found on first tile '//TRIM(involved_descendant_tiles(1)) &
              & //' were also passive cqs of this tile '//TRIM(involved_descendant_tiles(i_tile))//'. Found ' &
              & //int2string(nr_of_passive_vars)//' on first and '//int2string(i_passive_this_tile)//' on this. Please check.')
          END IF
        END IF

        leave => leave%Get_next_sibling_tile()
      END DO ! i_leave = 1,nr_of_leaves

      child => child%Get_next_sibling_tile()
    END DO ! i_child = 1,nr_of_children

    !>
    !> Assertion: all tiles involved in an lcc process need to be descendant tiles of the tile on which the lcc process is running.
    !>
    IF (.NOT. ALL(is_descendant)) THEN
      DO i_tile = 1,nr_of_tiles
        IF (.NOT. is_descendant(i_tile)) THEN
          CALL finish(TRIM(routine), &
            & 'Violation of assertion: Conceptually each involved tile needs to be a leave and a descendant tile '&
            & // 'of the tile on which the lcc proc is running -- did not found ' &
            & // TRIM(involved_descendant_tiles(i_tile))//'. Please check.' )
        END IF
      END DO
    END IF

    DEALLOCATE(is_descendant, active_var_process_ids)

    IF (debug_on()) THEN
      ! Following code is only a check! Therefore, its surrounded by debug flag, but could also have an own enable/disable check
      !
      ! Assert: check against potential active and required passive variables for each process with active vars
      !
      ! ... to enable active relocation, a process needs to 'expect' the variables assumed as active by a lcc process
      !     in order to provide a method to relocate them to passive variables. Passive variables used in an active
      !     relocation are therefore required to be passive on the calling lcc process.
      DO i_proc = 1, lcc_relocations%nr_of_procs_with_active_vars
        this_proc_name = TRIM(lcc_relocations%unique_proc_names_active_vars(i_proc))
        this_process_id = Get_process_id(this_proc_name)

        ! Count and collect active vars for this process
        nr_of_active_vars_this_proc = COUNT( lcc_relocations%active_vars_process_id == this_process_id )
        ALLOCATE(active_vars_this_proc(nr_of_active_vars_this_proc))
        ind = 0
        DO i_var = 1, lcc_relocations%nr_of_active_vars
          IF (lcc_relocations%active_vars_process_id(i_var) == this_process_id) THEN
            ind = ind + 1
            active_vars_this_proc(ind) = TRIM(lcc_relocations%active_vars(i_var)%p%name)
          END IF
        END DO

        ! ... and of passive vars of this proc
        nr_of_passive_vars_this_proc = COUNT(lcc_relocations%passive_vars_process_id == this_process_id)
        IF (nr_of_passive_vars_this_proc == 0) THEN
          ! Assertion: if there is an active var, at least one passive var is required, too.
          CALL finish(TRIM(routine), TRIM(lcc_relocations%name) &
            & // ' contains actively relocated variables for ' // TRIM(this_proc_name) &
            & // ' but no passively relocated sink variable for this proc - please check!')
        END IF
        ALLOCATE(passive_vars_this_proc(nr_of_passive_vars_this_proc))
        ind = 0
        DO i_var = 1, lcc_relocations%nr_of_passive_vars
          IF (lcc_relocations%passive_vars_process_id(i_var) == this_process_id) THEN
            ind = ind + 1
            passive_vars_this_proc(ind) = TRIM(lcc_relocations%passive_vars(i_var)%p%name)
          END IF
        END DO

        !... and check them
        SELECT CASE(tolower(this_proc_name))

        CASE ('carbon')
          CALL check_active_lcc_vars(lcc_relocations, this_proc_name, active_vars_this_proc, carbon_potential_active_vars)
          CALL check_passive_lcc_vars(lcc_relocations, this_proc_name, passive_vars_this_proc, carbon_required_passive_vars)

        CASE DEFAULT

          ! Assertion: Each process that has actively relocated conserved quantities needs
          !            to provide lists with potential active and required passive variables
          CALL finish(TRIM(routine), this_proc_name &
            & // ' contains variables that are actively relocated by ' // TRIM(lcc_relocations%name) &
            & // ' but did not specify lists with potential active and required passive vars - please check!')
        END SELECT

        DEALLOCATE(active_vars_this_proc, passive_vars_this_proc)
      END DO ! i_proc
    END IF

  END SUBROUTINE init_lcc

  ! ====================================================================================================== !
  !
  !> Initialize the lcc structure [[t_jsb_lcc_proc]] which is needed to relocate conservative quantities in
  !> the "update" routines of land cover change (lcc) processes. The routine collects the initial areas (cover
  !> fractions) of all involved (child) tiles and initializes/resets the relocation arrays for active or passive
  !> matter relocation.
  !>
  SUBROUTINE init_lcc_reloc(lcc_relocations, options, tile, initial_area)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations   !< LCC process structure for relocations
    TYPE(t_jsb_task_options),   INTENT(in)    :: options           !< Runtime parameters
    CLASS(t_jsb_tile_abstract), INTENT(in)    :: tile              !< Tile on which the LCC process runs
    REAL(lp),                   INTENT(inout) :: initial_area(:,:) !< Current area of involved tiles
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_lcc_reloc' !< Name of the routine

    CLASS(t_jsb_tile_abstract), POINTER :: current_tile, &         !< Pointer to the current tile
      &                                    parent                  !< Pointer to the parent tile
    LOGICAL                             :: on_second_level         !< Flag indicating if the child of a
                                                                   !< tile has children
    REAL(lp), DIMENSION(options%nc)     :: current_fract           !< Current cover fractions of a tile
    INTEGER :: i_tile, &    !< Loop index over tiles
      &         i_var, &    !< Loop index over variables
      &          iblk, &    !< Block index
      &           ics, &    !< Start index
      &           ice, &    !< End index
      &            nc, &    !< Number of cells
      &   nr_of_tiles       !< Number of tiles
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    nr_of_tiles = lcc_relocations%nr_of_tiles

    ! Assertion: initial area size needs to match number of tiles and vector length
    IF ((SIZE(initial_area,DIM=2) /= nr_of_tiles) .OR. (SIZE(initial_area,DIM=1) /= nc)) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Array "initial_area" has unexpected size. Please check!')
    ENDIF

    ! Collect initial areas of all involved tiles
    on_second_level = .FALSE.
    IF (tile%Has_children()) THEN
      current_tile => tile%Get_first_child_tile()
    ELSE
      current_tile => tile%self_tile
    END IF

    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile) .AND. i_tile < lcc_relocations%nr_of_tiles)

      IF (current_tile%Has_children()) THEN
        on_second_level = .TRUE.
        parent => current_tile
        current_tile => parent%Get_first_child_tile()
      END IF

      IF (one_of(current_tile%name, lcc_relocations%tile_names) > 0) THEN
        i_tile = i_tile + 1

        CALL current_tile%Get_fraction(ics, ice, iblk, fract=current_fract(:))

        initial_area(:, i_tile) = current_fract(:)
      END IF

      IF (.NOT. ASSOCIATED(current_tile%next_sibling_tile) .AND. on_second_level) THEN
        current_tile => parent
        on_second_level = .FALSE.
      END IF

      current_tile => current_tile%Get_next_sibling_tile()
    END DO

    ! Initialise / reset relocation arrays
    DO i_var = 1,lcc_relocations%nr_of_active_vars
      IF (lcc_relocations%active_vars(i_var)%p%var_type == REAL2D) THEN
        lcc_relocations%active_vars(i_var)%p%relocate_this(ics:ice,:,iblk) = 0.0_lp
      ELSE IF (lcc_relocations%active_vars(i_var)%p%var_type == REAL3D) THEN
        lcc_relocations%active_vars(i_var)%p%relocate_this_3D(ics:ice,:,:,iblk) = 0.0_lp
      ELSE
        CALL finish(TRIM(routine), 'Not supported variable type: '// &
            & int2string(lcc_relocations%active_vars(i_var)%p%var_type)//'. Please check!')
      END IF
    END DO
    DO i_var = 1,lcc_relocations%nr_of_passive_vars
      IF (lcc_relocations%passive_vars(i_var)%p%var_type == REAL2D) THEN
        lcc_relocations%passive_vars(i_var)%p%relocate_this(ics:ice,:,iblk) = 0.0_lp
      ELSE IF (lcc_relocations%passive_vars(i_var)%p%var_type == REAL3D) THEN
        lcc_relocations%passive_vars(i_var)%p%relocate_this_3D(ics:ice,:,:,iblk) = 0.0_lp
      ELSE
        CALL finish(TRIM(routine), 'Not supported variable type: '// &
            & int2string(lcc_relocations%passive_vars(i_var)%p%var_type)// '. Please check!')
      END IF
    END DO

  END SUBROUTINE init_lcc_reloc

  ! ====================================================================================================== !
  !>
  !> Starts the relocation of conserved quantities (CQs) by collecting the matter that will be transferred by
  !> active relocation (from one CQ variable on one tile to another CQ variable on another tile) and/ or
  !> passive relocation (to the same CQ variable on another tile).
  !>
  !
  SUBROUTINE start_lcc_reloc(lcc_relocations, options, lost_area, initial_area)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations   !< LCC process structure for relocations
    TYPE(t_jsb_task_options),   INTENT(in)    :: options           !< Runtime parameters
    REAL(lp),                   INTENT(in)    :: lost_area(:,:)    !< Area lost from involved tiles
    REAL(lp),                   INTENT(in)    :: initial_area(:,:) !< Initial area of involved tiles
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':start_lcc_reloc' !< Name of the routine

    INTEGER :: &
      & iblk,           & !< Block index
      & nc,             & !< Number of cells
      & nr_of_tiles       !< Number of tiles
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    iblk    = options%iblk
    nc      = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    nr_of_tiles = lcc_relocations%nr_of_tiles

    ! Assertion: dimensions of lost_area and initial_area need to match
    IF ((SIZE(lost_area,DIM=2) /= nr_of_tiles) .OR. (SIZE(initial_area,DIM=2) /= nr_of_tiles) &
      & .OR. (SIZE(lost_area,DIM=1) /= nc) .OR. (SIZE(initial_area,DIM=1) /= nc)) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Area arrays have unexpected size. Please check!')
    END IF

    ! Assertion: lost area is >= 0.0_lp and <= 1.0_lp
    IF (ANY(lost_area < 0.0_lp) .OR. ANY(lost_area > 1.0_lp + max_tolerated_cf_mismatch) &
      & .OR. ANY(initial_area < 0.0_lp) .OR. ANY(initial_area > 1.0_lp + max_tolerated_cf_mismatch)) THEN
      WRITE (message_text,*) 'Violation of assertion: Area arrays have values out of range ([0.0, 1.0]),', &
        & ' Minimum lost_area: ', MINVAL(lost_area), ' at ', MINLOC(lost_area),           &
        & ' Maximum lost_area: ', MAXVAL(lost_area), ' at ', MAXLOC(lost_area),           &
        & ' Minimum initial_area: ', MINVAL(initial_area), ' at ', MINLOC(initial_area),  &
        & ' Maximum initial_area: ', MAXVAL(initial_area), ' at ', MAXLOC(initial_area)
      CALL finish(routine, message_text)
    END IF

    ! Assertion: lost area <= initial_area
    ! -> JN-TODO: we might want/need to discuss the concept of the e-10 areas here
    IF (ANY((initial_area - lost_area) < -EPSILON(1._lp) )) THEN
      WRITE (message_text,*) 'Violation of assertion: Lost area exceeds available area' &
        & //' for lcc process '// TRIM(lcc_relocations%name) &
        & //' -- MINVAL: ', MINVAL(initial_area-lost_area), ' at ', MINLOC(initial_area-lost_area), &
        & ' (nc,tile)'
      CALL finish(routine, message_text)
    END IF

    ! NOTE: in the current implementation it is crucial that the active vars are filled first, because
    !       implementation of FLCC made it necessary that variables can be active and passive at the same time
    !       and currently active vars are relocated first and then passive vars are relocated
    CALL collect_to_be_transferred_matter(options, lost_area, initial_area, &
      & nr_of_tiles, lcc_relocations%nr_of_active_vars, lcc_relocations%active_vars)

    CALL collect_to_be_transferred_matter(options, lost_area, initial_area, &
      & nr_of_tiles, lcc_relocations%nr_of_passive_vars, lcc_relocations%passive_vars)

  END SUBROUTINE start_lcc_reloc

  ! ====================================================================================================== !
  !>
  !> Finalizes the relocation of conserved quantities (CQ) by distributing all matter that had been collected
  !> from shrinking tiles beforehand to extending tiles.
  !> At this point, matter from active relocation - i.e. matter to be transferred from an (active) variable on
  !> a shrinking tile to another variable on an expanding tile - has already been transferred to the receiving
  !> variable (e.g. via [[transfer_active_to_passive_onChunk]]). Thus only passive relocation remains to be
  !> done here, i.e. the relocation of matter from a (passive) variable on a shrinking tile to the same
  !> variable on an expanding tile.
  !>
  !> Note: the concept here is that the matter to be distributed has been collected from ALL shrinking tiles,
  !>       thus in case there are several shrinking tiles, there is no directed movement from one tile to
  !>       another tile but the sum of the matter from the shrinking tiles is distributed together to the
  !>       extending tiles proportional to the extending tile areas.
  !>
  !
  SUBROUTINE end_lcc_reloc(lcc_relocations, options, gained_area)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations   !< LCC process structure for relocations
    TYPE(t_jsb_task_options),   INTENT(in)    :: options           !< Runtime parameters
    REAL(lp),                   INTENT(in)    :: gained_area(:,:)  !< Area gained on involved tiles
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':end_lcc_reloc'  !< Name of the routine

    LOGICAL :: not_zero     !< Boolean to check if any active variable still carries any not relocated matter

    INTEGER :: i_tile, &    !< Loop index over tiles
      &         i_var, &    !< Loop index over variables
      &          iblk, &    !< Block index
      &           ics, &    !< Start index
      &           ice, &    !< End index
      &            nc, &    !< Number of cells
      &   nr_of_tiles, &    !< Number of tiles
      &         i_dim, &    !< Index of additional dimension in case of REAL3D typed passive variables
      &         n_dim       !< Number of additional dimension in case of REAL3D typed passive variables
    REAL(lp) :: passive_var_sum(options%nc), &  !< Sum of a passively relocated variable over all involved tiles
      &         gained_area_sum(options%nc)     !< Sum of all gained area for all involved tiles
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    nr_of_tiles = lcc_relocations%nr_of_tiles
    gained_area_sum = SUM(gained_area(:, :), DIM=2)

    ! Assertion: dimensions of gained_area need to be (nr_of_tiles, nc)
    IF ((SIZE(gained_area,DIM=2) /= nr_of_tiles) .OR. (SIZE(gained_area,DIM=1) /= nc)) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Area array has unexpected size. Please check!')
    END IF

    ! Assertion: relocated matter of active variables is expected to be zero here (as it should have been relocated
    !            to a passive variable before the routine call)
    DO i_var = 1,lcc_relocations%nr_of_active_vars
      DO i_tile = 1,nr_of_tiles
        not_zero = .TRUE.

        IF (lcc_relocations%active_vars(i_var)%p%var_type == REAL2D) THEN
          not_zero = ANY(lcc_relocations%active_vars(i_var)%p%relocate_this(ics:ice,i_tile,iblk) /= 0.0_lp)
        ELSE IF (lcc_relocations%active_vars(i_var)%p%var_type == REAL3D) THEN
          not_zero = ANY(lcc_relocations%active_vars(i_var)%p%relocate_this_3D(ics:ice,:,i_tile,iblk) /= 0.0_lp)
        ELSE
          CALL finish(TRIM(routine), 'Not supported variable type: '// &
              & int2string(lcc_relocations%active_vars(i_var)%p%var_type)// '. Please check!')
        END IF

        IF (not_zero) THEN
          CALL finish(TRIM(routine), 'Violation of assertion: All active variables need to already have been ' &
            & // 'relocated upon call of this routine. However, found != 0.0 values for active var '           &
            & // TRIM(lcc_relocations%active_vars(i_var)%p%name) // ' on tile '                                &
            & // TRIM(lcc_relocations%tile_names(i_tile)) // ' for lcc process '                               &
            & // TRIM(lcc_relocations%name) // '. Please check!')
        END IF
      END DO
    END DO

    DO i_var = 1,lcc_relocations%nr_of_passive_vars
      ! Matter from passive variables is proportionally distributed to these passive variables on extending tiles
      IF (lcc_relocations%passive_vars(i_var)%p%var_type == REAL2D) THEN
        passive_var_sum(:) = SUM(lcc_relocations%passive_vars(i_var)%p%relocate_this(ics:ice,:,iblk), DIM=2)
        DO i_tile = 1,nr_of_tiles
          WHERE (gained_area(:, i_tile) > 0.0_lp)
            lcc_relocations%passive_vars(i_var)%p%var_on_tile(i_tile)%p%ptr2d(ics:ice,iblk) &
              & = lcc_relocations%passive_vars(i_var)%p%var_on_tile(i_tile)%p%ptr2d(ics:ice,iblk) &
              & + passive_var_sum(:) * (gained_area(:, i_tile) / gained_area_sum(:))
          END WHERE
        END DO
        lcc_relocations%passive_vars(i_var)%p%relocate_this(ics:ice,:,iblk) = 0.0_lp

      ELSE IF (lcc_relocations%passive_vars(i_var)%p%var_type == REAL3D) THEN
        n_dim = SIZE(lcc_relocations%passive_vars(i_var)%p%relocate_this_3D(ics:ice,:,:,iblk), DIM=2)

        DO i_dim = 1,n_dim
          passive_var_sum(:) = SUM(lcc_relocations%passive_vars(i_var)%p%relocate_this_3D(ics:ice,i_dim,:,iblk), DIM=2)
          DO i_tile = 1,nr_of_tiles
            WHERE (gained_area(:, i_tile) > 0.0_lp)
              lcc_relocations%passive_vars(i_var)%p%var_on_tile(i_tile)%p%ptr3d(ics:ice,i_dim,iblk) &
                & = lcc_relocations%passive_vars(i_var)%p%var_on_tile(i_tile)%p%ptr3d(ics:ice,i_dim,iblk) &
                & + passive_var_sum(:) * (gained_area(:, i_tile) / gained_area_sum(:))
            END WHERE
          END DO
          lcc_relocations%passive_vars(i_var)%p%relocate_this_3D(ics:ice,i_dim,:,iblk) = 0.0_lp
        END DO
      ELSE
        CALL finish(TRIM(routine), 'Not supported variable type: '// &
            & int2string(lcc_relocations%passive_vars(i_var)%p%var_type)// '. Please check!')
      END IF
    END DO

  END SUBROUTINE end_lcc_reloc


  ! ====================================================================================================== !
  !>
  !> The routine handles the relocation of matter from active variables (source variables of active matter
  !> relocation) to passive variables (receiving variables of the active matter relocation). As the treatment
  !> of different conserved quantity variables is very specific we are calling relocation functions of
  !> the respective land cover change (lcc) processes.
  !>
  !> The routine is only needed with active matter relocation. It only handles the matter relocation from
  !> the active to the passive variable, the relocation from shrinking to expanding tiles is however handled
  !> later in routine [[end_lcc_reloc]].
  !>
  ! Note: the tile and tile index passed into the routine do not correspond to the tile on which the lcc
  !       process is calculated but to the current descendant tile for which the matter transfer is executed.
  !
  SUBROUTINE transfer_active_to_passive_onChunk(lcc_relocations, tile, i_tile, options)
    USE mo_carbon_interface,    ONLY: carbon_transfer_from_active_to_passive_vars_onChunk

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations !< LCC structure of the calling LCC process
    CLASS(t_jsb_tile_abstract), INTENT(in)    :: tile            !< Tile on which relocation should be conducted
    INTEGER,                    INTENT(in)    :: i_tile          !< Index of the tile in LCC structure
    TYPE(t_jsb_task_options),   INTENT(in)    :: options         !< Runtime parameters
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':transfer_active_to_passive_onChunk' !< Name of the routine

    CHARACTER(len=:), ALLOCATABLE :: this_proc_name !< Name of the process of the currently treated active variable
    INTEGER :: i_proc,  &                           !< Process index
      &        i_var,   &                           !< Variable index
      &        proc_id, &                           !< Process ID
      &        ics,     &                           !< Start index
      &        ice,     &                           !< End index
      &        iblk                                 !< Block index
    LOGICAL :: has_non_zero_active_var              !< Flag indicating if any active variable has non-zero content
    ! -------------------------------------------------------------------------------------------------- !
    ics = options%ics
    ice = options%ice
    iblk = options%iblk

    IF (debug_on() .AND. options%iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    DO i_proc = 1, lcc_relocations%nr_of_procs_with_active_vars

      this_proc_name = TRIM(lcc_relocations%unique_proc_names_active_vars(i_proc))
      proc_id = Get_process_id(this_proc_name)

      !JN-TODO: required?
      ! Count number of active vars with non-zero relocations for this process
      has_non_zero_active_var = .FALSE.
      DO i_var = 1, lcc_relocations%nr_of_active_vars
        IF (lcc_relocations%active_vars_process_id(i_var) == proc_id) THEN
          IF (ANY(lcc_relocations%active_vars(i_var)%p%relocate_this(ics:ice,i_tile,iblk) .NE. 0.0_lp ) ) THEN
            has_non_zero_active_var = .TRUE.
          END IF
        END IF
      END DO

      IF (has_non_zero_active_var) THEN

        SELECT CASE(tolower(this_proc_name))

        CASE ('carbon')
          CALL carbon_transfer_from_active_to_passive_vars_onChunk(lcc_relocations, tile, i_tile, options)

        CASE DEFAULT
          !>
          !> Assertion: Each process that has conserved quantities used as actively relocated CQs needs
          !>            to provide a specific routine for the active relocation
          !>
          CALL finish(TRIM(routine), this_proc_name &
            & // ' contains variables that are actively relocated by ' // TRIM(lcc_relocations%name) &
            & // ' but no specific routine is defined - please check!')
        END SELECT

      END IF ! has_non_zero_active_var
    END DO ! i_proc

  END SUBROUTINE transfer_active_to_passive_onChunk

  ! ====================================================================================================== !
  !
  !> Helper routine to count the descendants of a tile down to two layers, i.e. the number of children and
  !> grandchildren. More than two layers are not supported.
  !
  SUBROUTINE count_descendants_for_up_to_2layers(tile, nr_of_descendants)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(INOUT) :: tile
    INTEGER,                    INTENT(OUT)   :: nr_of_descendants
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), POINTER   :: current_tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':count_descendants_for_up_to_2layers'
    ! -------------------------------------------------------------------------------------------------- !

    ! count leaves
    nr_of_descendants = 0
    current_tile => tile%Get_first_child_tile()
    DO WHILE (ASSOCIATED(current_tile))
      IF (current_tile%Has_children()) THEN
        nr_of_descendants = nr_of_descendants + current_tile%Get_no_of_children()
      ELSE
        nr_of_descendants = nr_of_descendants + 1
      END IF
      current_tile => current_tile%Get_next_sibling_tile()
    END DO

  END SUBROUTINE count_descendants_for_up_to_2layers

  ! ====================================================================================================== !
  !
  !> Helper routine to collect the names of descendant tiles down to two layers, i.e. the names of the
  !> children and grandchildren. More than two layers are not supported.
  !
  SUBROUTINE collect_names_of_descendant_leaves_up_to_2layers(tile, descendant_leave_tiles)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(INOUT)  :: tile
    CHARACTER(len=SHORT_NAME_LEN), INTENT(OUT) :: descendant_leave_tiles(:)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: ind_list, ind_end
    CLASS(t_jsb_tile_abstract), POINTER   :: current_tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':collect_names_of_descendant_leaves_up_to_2layers'
    ! -------------------------------------------------------------------------------------------------- !

    current_tile => tile%Get_first_child_tile()
    ind_list = 1
    DO WHILE (ASSOCIATED(current_tile))
      IF (current_tile%Has_children()) THEN
        ! Assertion: Only implemented for 1 or 2 hierarchical layers!
        IF (current_tile%first_child%Has_children()) THEN
          CALL finish(TRIM(routine), &
            & 'Violation of assertion: number of involved layers exceeds two below given tile '&
            & //TRIM(tile%name)//'. Which is not implemented. Please check.')
        END IF

        ind_end = ind_list+current_tile%Get_no_of_children()-1
        CALL current_tile%Get_children_names(descendant_leave_tiles(ind_list:ind_end))
        ind_list = ind_end + 1
      ELSE
        ! Add name of the tile itself
        descendant_leave_tiles(ind_list) = current_tile%name
        ind_list = ind_list + 1
      END IF
      current_tile => current_tile%Get_next_sibling_tile()
    END DO

  END SUBROUTINE collect_names_of_descendant_leaves_up_to_2layers

  ! ====================================================================================================== !
  !>
  !> Helper routine initializing the structures required to collect information and memory pointers for a
  !> specific conserved quantity (CQ).
  !>
  !> This routine creates a new collector for memory pointers for the given CQ variable which is appended
  !> to a given list of CQ variables (var_list). The calling routine thereby can pass the list of active or
  !> passive variables contained in the encompassing land cover change (lcc) structure.
  !> In addition to creating the new collector, meta data such as the name and the variable type of the CQ is
  !> added using the information stored about the variable in the passed conserved_quantities list from the
  !> first tile involved in the relocation for this lcc process.
  !>
  !> The routine is called during the initialization of the respective lcc process.
  !>
  SUBROUTINE add_cq_to_lcc_structure(nr_of_tiles, nproma, nblks, i_tile, i_cqt, i_var, i_of_kind, &
    & conserved_quantities, cqt_of_vars, vars_process_id, var_list, var_names )

    USE mo_jsb_lcc_class,       ONLY: t_jsb_lcc_var

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in) :: &
      & nr_of_tiles, &          !< Number of descendant tiles involved in this lcc process
      & nproma, &               !< Vector length for the blocks
      & nblks, &                !< Number of blocks
      & i_tile, &               !< TODO: remove -- this is by definition 1 because it is the first decendant tile
      & i_cqt, &                !< Index of the conserved quantity type (CQT) to which this CQ belongs to in
                                !< the collection of conserved variables on the first decendant tile
      & i_var, &                !< Index of the CQ variable in the collection of conserved variables
                                !< on the first decendant tile for the given CQT
      & i_of_kind               !< Index in the variable list for the newly created collection for this CQ
    TYPE(t_jsb_consQuan_p), INTENT(inout) :: &
      & conserved_quantities(:) !< Array of conserved quantities on the first decendant tile involved in the
                                !< lcc process for which the lcc structure is initalised here
    INTEGER,                INTENT(inout) :: &
      & cqt_of_vars(:),  &      !< List collecting the CQT ID of this CQ
      & vars_process_id(:)      !< List collecting the process ID of this CQ
    TYPE(t_jsb_lcc_var_p),  INTENT(inout) :: &
      & var_list(:)             !< List to add this CQ variable to (either list of active or passive variables
                                !< on a lcc structure)
    CHARACTER(len=*),       INTENT(inout) :: var_names(:) !< List collecting variable names
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_var), POINTER  :: a_lcc_var !< A new list to collect memory pointers for the current CQ
    INTEGER :: extra_dim                       !< Length of the extra dimension in case of REAL3D variables
    CHARACTER(len=*), PARAMETER :: routine = modname//':add_cq_to_lcc_structure' !< Name of the routine
    ! -------------------------------------------------------------------------------------------------- !

    cqt_of_vars(i_of_kind) = conserved_quantities(i_cqt)%p%type_id
    vars_process_id(i_of_kind) = conserved_quantities(i_cqt)%p%associated_process(i_var)

    ALLOCATE(a_lcc_var)
    var_list(i_of_kind)%p => a_lcc_var
    var_list(i_of_kind)%p%name = conserved_quantities(i_cqt)%p%cq_vars(i_var)%p%name
    var_list(i_of_kind)%p%standard_name = conserved_quantities(i_cqt)%p%cq_vars(i_var)%p%standard_name
    var_list(i_of_kind)%p%var_type = conserved_quantities(i_cqt)%p%var_type
    var_list(i_of_kind)%p%transfer_all = conserved_quantities(i_cqt)%p%transfer_all
    var_names(i_of_kind) = var_list(i_of_kind)%p%name

    ALLOCATE(var_list(i_of_kind)%p%var_on_tile(nr_of_tiles))
    var_list(i_of_kind)%p%var_on_tile(i_tile)%p => conserved_quantities(i_cqt)%p%cq_vars(i_var)%p

    IF (conserved_quantities(i_cqt)%p%var_type == REAL2D) THEN
      ALLOCATE(var_list(i_of_kind)%p%relocate_this(nproma,nr_of_tiles,nblks))
      var_list(i_of_kind)%p%relocate_this(:,:,:) = 0.0_lp
    ELSE IF (conserved_quantities(i_cqt)%p%var_type == REAL3D) THEN
      extra_dim = SIZE(conserved_quantities(i_cqt)%p%cq_vars(i_var)%p%ptr3d,DIM=2)
      ALLOCATE(var_list(i_of_kind)%p%relocate_this_3D(nproma,extra_dim,nr_of_tiles,nblks))
      var_list(i_of_kind)%p%relocate_this_3D(:,:,:,:) = 0.0_lp
    ELSE
      CALL finish(TRIM(routine), 'Not supported variable type: '// &
          & int2string(conserved_quantities(i_cqt)%p%var_type)// '. Please check!')
    END IF

  END SUBROUTINE add_cq_to_lcc_structure


  ! ====================================================================================================== !
  !>
  !> Helper routine to collect memory pointers on a given tile for a specific conserved quantity type (CQT).
  !>
  !> Collects the pointers to the memory of all conserved quantity (CQ) variables of a given conserved
  !> quantity type (CQT) on a given tile. This routine thereby only needs to set the pointers, the
  !> required structure has already been allocated in the routine [[add_cq_to_lcc_structure]].
  !>
  !> The routine is only called upon init of a land cover change (lcc) process.
  !>
  SUBROUTINE collect_pointers_for_cqs(i_cqt, i_tile, i_var_this_tile, this_tile_name, &
    & first_involved_child_name, conserved_quantities, var_list, nr_of_vars, this_case)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                INTENT(in)    :: i_cqt, &        !< Index of the current CQT
      &                                      i_tile          !< Index of the current tile
    INTEGER,                INTENT(inout) :: i_var_this_tile !< Counter for variables on the current tile
    CHARACTER(len=*),       INTENT(in)    :: this_tile_name  !< Name of the current tile
    CHARACTER(len=*),       INTENT(in)    :: first_involved_child_name !< Name of first involved child tile
    TYPE(t_jsb_consQuan_p), INTENT(inout) :: conserved_quantities(:)   !< Array CQs sorted by CQTs on the
                                                                       !< current tile
    TYPE(t_jsb_lcc_var_p),  INTENT(inout) :: var_list(:) !< List of all (active or passive) CQs for the
                                                         !< lcc process - part of the lcc structure
    INTEGER,                INTENT(in)    :: nr_of_vars  !< Number of variables in the var_list
    CHARACTER(len=*),       INTENT(in)    :: this_case   !< Case identifier ('active' or 'passive')
                                                         !< only required in case of an error message
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: nr_of_vars_on_cq, &     !< Number of variables on current CQT -- TODO: rename -> of_cqt
      &        i_var, &                !< Looping index for variables in CQT list
      &        i                       !< Looping index for variables in the given lcc variable list
    LOGICAL :: has_searched_var        !< Flag indicating whether a searched variable was found
    CHARACTER(len=VARNAME_LEN) :: thisName !< Temporary variable name buffer
    CHARACTER(len=*), PARAMETER :: routine = modname//':collect_pointers_for_cqs' !< Name of the routine
    ! -------------------------------------------------------------------------------------------------- !
    nr_of_vars_on_cq = conserved_quantities(i_cqt)%p%no_of_vars
    DO i_var = 1,nr_of_vars_on_cq
      i_var_this_tile = i_var_this_tile + 1
      thisName = conserved_quantities(i_cqt)%p%cq_vars(i_var)%p%standard_name

      has_searched_var = .FALSE.
      DO i = 1,nr_of_vars
        IF (TRIM(var_list(i)%p%standard_name) .EQ. TRIM(thisName)) THEN
          has_searched_var = .TRUE.
          var_list(i)%p%var_on_tile(i_tile)%p => conserved_quantities(i_cqt)%p%cq_vars(i_var)%p
          EXIT
        END IF
      END DO
      ! Assertion: all passive and active variables need to exist on all tiles
      IF (.NOT. has_searched_var) THEN
        CALL finish(TRIM(routine), &
          & 'Violation of assertion: Conceptually each tile needs to contain the same active and passive vars. ' &
          & //this_case //' var ' // TRIM(thisName) // ' of tile ' // TRIM(this_tile_name) // &
          & ' does not exist on tile ' // first_involved_child_name //'. Please check.' )
      END IF
    END DO

  END SUBROUTINE collect_pointers_for_cqs


  ! ====================================================================================================== !
  !>
  !> Collect to be transferred matter from a given list of conserved quantity (CQ) variables
  !>
  !> The matter to be relocated is removed from all shrinking tiles and collected in the lcc structure to be
  !> later distributed to the expanding tiles.
  !>
  SUBROUTINE collect_to_be_transferred_matter(options, lost_area, initial_area, &
    & nr_of_tiles, nr_of_vars, var_list)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_task_options),   INTENT(in)    :: options           !< Runtime parameters
    REAL(lp),                   INTENT(in)    :: lost_area(:,:)    !< Area lost on involved tiles
    REAL(lp),                   INTENT(in)    :: initial_area(:,:) !< Initial area of involved tiles
    INTEGER,                    INTENT(in)    :: nr_of_tiles, &    !< Number of tiles
      &                                          nr_of_vars        !< Number of variables
    TYPE(t_jsb_lcc_var_p),      INTENT(inout) :: var_list(:)       !< List of LCC variable pointers
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: i_dim, &           !< Index of additional dimension (for REAL3D typed variables)
      &        n_dim              !< Size of additional dimension (for REAL3D typed variables)
    INTEGER :: nc, &              !< Number of cells
      &        ic, &              !< Cell index
      &        ics, &             !< Start index
      &        ice, &             !< End index
      &        iblk, &            !< Block index
      &        i_var, &           !< Variable index
      &        i_tile             !< Tile index
    REAL(lp) :: to_be_transferred !< Matter to be transferred
    CHARACTER(len=*), PARAMETER :: routine = modname//':collect_to_be_transferred_matter' !< Name of the routine
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    DO i_var = 1,nr_of_vars
      DO i_tile = 1,nr_of_tiles
        DO ic = 1,nc
          IF (lost_area(ic, i_tile) > 0.0_lp) THEN
            ! Collect to be transferred matter
            IF (var_list(i_var)%p%var_type == REAL2D) THEN
              ! Depending on the transfer_all flag either all matter from the variable is transferred
              ! or only a share proportionally to the lost area
              IF (var_list(i_var)%p%transfer_all) THEN
                to_be_transferred = var_list(i_var)%p%var_on_tile(i_tile)%p%ptr2d(ics+ic-1,iblk)
              ELSE
                to_be_transferred = var_list(i_var)%p%var_on_tile(i_tile)%p%ptr2d(ics+ic-1,iblk) &
                  &   * (lost_area(ic, i_tile) / initial_area(ic, i_tile))
              END IF
              var_list(i_var)%p%relocate_this(ics+ic-1,i_tile,iblk) = to_be_transferred
              ! ... and subtract it
              var_list(i_var)%p%var_on_tile(i_tile)%p%ptr2d(ics+ic-1,iblk) &
                & = var_list(i_var)%p%var_on_tile(i_tile)%p%ptr2d(ics+ic-1,iblk) - to_be_transferred

            ELSE IF (var_list(i_var)%p%var_type == REAL3D) THEN
              n_dim = SIZE(var_list(i_var)%p%var_on_tile(i_tile)%p%ptr3d(ics+ic-1,:,iblk))
              DO i_dim = 1, n_dim
                ! Depending on the transfer_all flag either all matter from the variable is transferred
                ! or only a share proportionally to the lost area
                IF (var_list(i_var)%p%transfer_all) THEN
                  to_be_transferred = var_list(i_var)%p%var_on_tile(i_tile)%p%ptr3d(ics+ic-1,i_dim,iblk)
                ELSE
                  to_be_transferred = var_list(i_var)%p%var_on_tile(i_tile)%p%ptr3d(ics+ic-1,i_dim,iblk) &
                    &   * (lost_area(ic, i_tile) / initial_area(ic, i_tile))
                END IF
                var_list(i_var)%p%relocate_this_3D(ics+ic-1,i_dim,i_tile,iblk) = to_be_transferred
                ! ... and subtract it
                var_list(i_var)%p%var_on_tile(i_tile)%p%ptr3d(ics+ic-1,i_dim,iblk) &
                  & = var_list(i_var)%p%var_on_tile(i_tile)%p%ptr3d(ics+ic-1,i_dim,iblk) - to_be_transferred
              END DO
            ELSE
              CALL finish(TRIM(routine), 'Not supported variable type: '// &
                  & int2string(var_list(i_var)%p%var_type)// '. Please check!')
            END IF
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE collect_to_be_transferred_matter


  ! ====================================================================================================== !
  !>
  !> Test routine asserting that declared active variables are listed as potentially active variables.
  !>
  !> The processes to which the conserved quantity (CQ) variables belong (e.g. CARBON) define a list of all
  !> potential variables for active matter relocation of this process. This routine ensures that all variables
  !> defined as active CQ variables in the calling land cover change (LCC) process are in this list.
  !>
  SUBROUTINE check_active_lcc_vars(lcc_relocations, this_proc_name, active_vars_on_this_lcc_proc, potential_active_vars)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc), INTENT(INOUT) :: &
      & lcc_relocations                    !< The lcc structure of the calling lcc process
    CHARACTER(len=*), INTENT(IN) :: &
      & this_proc_name,                  & !< Name of the process whose active variables are checked
      & active_vars_on_this_lcc_proc(:), & !< Active variables of the calling lcc process that are
                                           !< defined in the tested process
      & potential_active_vars(:)           !< List of potentially active variables defined in the tested
                                           !< process
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':check_active_lcc_vars'  !< Name of the routine
    INTEGER                     :: i_var   !< Looping index for variables
    ! -------------------------------------------------------------------------------------------------- !

    ! Assertion: calling lcc process should not call this routine with active vars
    ! that are not contained in the specified (and handled) potential active vars
    DO i_var = 1,SIZE(active_vars_on_this_lcc_proc)
      IF (.NOT. ANY(active_vars_on_this_lcc_proc(i_var) .EQ. potential_active_vars)) THEN
        CALL finish(TRIM(routine), 'Violation of assertion: ' // TRIM(lcc_relocations%name)       &
          & // ' called with unexpected active var ' // TRIM(active_vars_on_this_lcc_proc(i_var)) &
          & // ' of process ' // this_proc_name)
      END IF
    END DO

  END SUBROUTINE check_active_lcc_vars

  ! ====================================================================================================== !
  !>
  !> Test routine asserting that passive variables required for active relocation are listed as passive
  !> variables in the lcc process.
  !>
  !> The processes to which the conserved quantity (CQ) variables belong (e.g. CARBON) define a list of the
  !> required passive variables for active matter relocation. This routine ensures that all variables in this
  !> list are defined as passive CQ variables in the calling land cover change (LCC) process.
  !>
  SUBROUTINE check_passive_lcc_vars(lcc_relocations, this_proc_name, passive_vars_on_this_lcc_proc , required_passive_vars)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc), INTENT(INOUT) :: &
      & lcc_relocations                     !< The lcc structure of the calling lcc process
    CHARACTER(len=*), INTENT(IN) :: &
      & this_proc_name,                   & !< Name of the process whose active variables are checked
      & passive_vars_on_this_lcc_proc(:), & !< Passive variables of the calling lcc process that are
                                            !< defined in the tested process
      & required_passive_vars(:)            !< List of required passive variables defined in the
                                            !< tested process
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':check_passive_lcc_vars' !< Name of the routine
    INTEGER                     :: i_var    !< Looping index for variables
    ! -------------------------------------------------------------------------------------------------- !

    ! Assertion: calling lcc process needs to call this routine with all required passive vars listed
    DO i_var = 1,SIZE(required_passive_vars)
      IF (.NOT. ANY(required_passive_vars(i_var) .EQ. passive_vars_on_this_lcc_proc)) THEN
        CALL finish(TRIM(routine), 'Violation of assertion: ' // TRIM(lcc_relocations%name) &
          & // ' did not list required passive var ' // TRIM(required_passive_vars(i_var))  &
          & // ' of process ' // this_proc_name)
      END IF
    END DO

  END SUBROUTINE check_passive_lcc_vars

#endif
END MODULE mo_jsb_lcc
