!> Types and helper subroutines for lcc processes with actively or passively relocated conserved quantities.
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
!>#### Contains types required for dealing with actively or passively relocated conserved quantities upon lcc.
!>
!> In addition this module provides some helper subroutines dealing with special purpose objectives related to lcc.
MODULE mo_jsb_lcc_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: jp
  USE mo_util,                ONLY: one_of, int2string
  USE mo_exception,           ONLY: finish, message
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_var_class,       ONLY: t_jsb_var_p
  USE mo_jsb_varlist,         ONLY: VARNAME_LEN
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN
  USE mo_jsb_math_constants,  ONLY: one_year
  USE mo_jsb_var_class,       ONLY: REAL2D, REAL3D

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: min_daily_cf_change, min_annual_cf_change, max_tolerated_fract_mismatch, max_tolerated_cf_mismatch,       &
    & t_jsb_lcc_proc, t_jsb_lcc_proc_p, t_jsb_lcc_var, t_jsb_lcc_var_p,                                               &
    & transfer_from_active_to_passive_var_onChunk, collect_matter_of_active_vars_onChunk,                             &
    & copy_from_active_to_passive_var_onChunk

  REAL(jp), PARAMETER :: min_daily_cf_change = 1.0E-12_jp
    !< Min change in cover fractions for daily lcc, i.e. if the change is smaller than this threshold the lcc is not applied
  REAL(jp), PARAMETER :: min_annual_cf_change = min_daily_cf_change * one_year
    !< Min change in cover fractions for annual lcc, i.e. if the change is smaller than this threshold the lcc is not applied
  REAL(jp), PARAMETER :: max_tolerated_cf_mismatch = 1.0E-08_jp
    !< Max mismatch in cover fractions which is tolerated when comparing different cover fraction sums
  REAL(jp), PARAMETER :: max_tolerated_fract_mismatch = 1.0E-14_jp
    !< Maximum tolerated mismatch in forest age fraction sums (only used when running with FAGE)

  !> Type to collect variable pointers for one specific conserved quantity (CQ) on all tiles involved in
  !> a certain lcc process which actively or passively relocates this CQ.
  !> In addition, the type is used to gather the amount of to be relocated matter for this CQ.
  TYPE :: t_jsb_lcc_var
    CHARACTER(len=VARNAME_LEN)    :: name              !< Name of the conserved quantity (CQ; variable name)
    CHARACTER(len=VARNAME_LEN)    :: standard_name     !< Standard name of the variable
    INTEGER                       :: var_type          !< Variable type (supported: REAL2D or REAL3D)
    LOGICAL                       :: transfer_all      !< True: The complete matter from area loosing tiles
                                                       !< should be collected for relocation for this variable
                                                       !< False: Only the fraction of matter proportional to
                                                       !< the lost area should be collected (default)
!    INTEGER :: type_id = -1     !< one of the CQ_TYPEs ( mo_jsb_cqt_class ) -- JN-TODO: required? -> also in encompassing type...
    TYPE(t_jsb_var_p), ALLOCATABLE :: var_on_tile(:)   !< Pointer to the variable in the memory on each tile
                                                       !< involved in the lcc process
    REAL(jp), ALLOCATABLE :: relocate_this(:,:,:)      !< Collector for the amount of matter to be relocated
                                                       !< for a REAL2D variable from each involved tile.
                                                       !< Dimensions: (nproma, nr_of_tiles, nblks)
    REAL(jp), ALLOCATABLE :: relocate_this_3D(:,:,:,:) !< Collector for the amount of matter to be relocated
                                                       !< for a REAL3D variable from each involved tile.
                                                       !< Dimensions: (nproma, extra_dim, nr_of_tiles, nblks)
  END TYPE t_jsb_lcc_var

  !> Helper type used to collect pointers to actively or passively relocated conserved quantities.
  TYPE :: t_jsb_lcc_var_p
    TYPE(t_jsb_lcc_var), POINTER :: p => NULL()
  END TYPE t_jsb_lcc_var_p

  !> Main type used to organise land cover changes (lcc) for an encompassing lcc process with active
  !> and/or passive relocations (also called "lcc structure").
  !> It is used to collect information on the process, conserved quantities, involved tiles, as well as area changes.
  TYPE :: t_jsb_lcc_proc
    CHARACTER(len=:), ALLOCATABLE :: name                       !< Name of the structure: "process ID + _lcc"
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: tile_names(:) !< Names of the tiles involved with this lcc process
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: unique_proc_names_active_vars(:)
        !< Collection of unique names of the processes carrying variables that are actively relocated by this lcc process
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: active_vars_names(:)
        !< List of active var names - same order as other lists for active vars (size: nr_of_active_vars)
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: passive_vars_names(:)
        !< List of passive var names - same order as other lists for passive vars (size: nr_of_passive_vars)
    INTEGER, ALLOCATABLE :: active_vars_cqt(:)  ! JN-TODO: remove - CQTs are not relevant in this dictionary (only for construction)
    INTEGER, ALLOCATABLE :: passive_vars_cqt(:) ! JN-TODO: remove
    INTEGER, ALLOCATABLE :: active_vars_process_id(:)  !< ID of the process carrying the active var (size: nr_of_active_vars)
    INTEGER, ALLOCATABLE :: passive_vars_process_id(:) !< ID of the process carrying the passive var (size: nr_of_passive_vars)
    INTEGER :: nr_of_tiles = 0          !< Number of tiles involved with this lcc process
    INTEGER :: nr_of_procs_with_active_vars = 0
        !< Number of processes carrying variables that are actively relocated by this lcc process
    INTEGER :: nr_of_active_vars = 0    !< Number of variables with active relocation in the encompassing lcc process
    INTEGER :: nr_of_passive_vars = 0   !< Number of variables with passive relocation in the encompassing lcc process
    TYPE(t_jsb_lcc_var_p), ALLOCATABLE :: active_vars(:)
        !< Collection of lcc var types that are actively relocated for the encompassing lcc process (size: nr_of_active_vars)
    TYPE(t_jsb_lcc_var_p), ALLOCATABLE :: passive_vars(:)
        !< Collection of lcc var types that are passively relocated for the encompassing lcc process  (size: nr_of_passive_vars)
    REAL(jp), ALLOCATABLE :: initial_area(:,:,:) !< Initial area before lcc event (nproma,nr_of_tiles,nblks)
    REAL(jp), ALLOCATABLE :: gained_area(:,:,:)  !< Area gained during lcc event (nproma,nr_of_tiles,nblks)
    REAL(jp), ALLOCATABLE :: lost_area(:,:,:)    !< Area lost during lcc event (nproma,nr_of_tiles,nblks)
  END TYPE t_jsb_lcc_proc

  !> Helper type used to collect pointers to lcc structures on different lcc processes
  TYPE t_jsb_lcc_proc_p
    TYPE(t_jsb_lcc_proc), POINTER :: p
  END TYPE t_jsb_lcc_proc_p

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_lcc_class'

CONTAINS


  ! ====================================================================================================== !
  !
  !> Relocates matter from the given active to the given passive variable.
  !> An optional passed fraction can be used to determine the share of the matter on the active variable to be
  !> relocated.
  !>
  !
  SUBROUTINE transfer_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice,           &
      &                                                  active_var_name, passive_var_name, transfer_fract, &
      &                                                  conversion_fact)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc), INTENT(inout) :: lcc_relocations !< The lcc structure used to distribute matter
                                                           !< from an active to a passive variable
    INTEGER,              INTENT(in)    :: i_tile  !< Index of the tile in the lcc structure for which the
                                                   !< relocation should be conducted
    INTEGER,              INTENT(in)    :: iblk    !< Number of current block (chunk)
    INTEGER,              INTENT(in)    :: ics     !< Index of chunk start
    INTEGER,              INTENT(in)    :: ice     !< Index of chunk end
    CHARACTER(len=*),     INTENT(in)    :: active_var_name  !< Active (source) variable name
    CHARACTER(len=*),     INTENT(in)    :: passive_var_name !< Passive (sink) variable name
    REAL(jp),      INTENT(in), OPTIONAL :: transfer_fract   !< Fraction of matter transferred from the active
                                                            !< variable to the passive variable.
                                                            !< If not passed, 1.0 is assumed
    REAL(jp),      INTENT(in), OPTIONAL :: conversion_fact  !< Conversion fraction in case that the active
                                                            !< variable has a different unit then the passive
                                                            !< variable. If not passed, 1.0 is assumed
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':transfer_from_active_to_passive_var_onChunk'

    REAL(jp) :: &
      & this_transfer_fract, &  !< Fraction of matter transferred from the active variable
      & this_conversion_fact    !< Conversion factor in case active and passive variables have different units
    INTEGER  :: &
      & i_active, &             !< Index of this active variable in the lcc structure
      & i_passive               !< Index of this passive variable in the lcc structure
    ! -------------------------------------------------------------------------------------------------- !

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    i_active = one_of( active_var_name , lcc_relocations%active_vars_names)
    i_passive = one_of( passive_var_name , lcc_relocations%passive_vars_names)

    !>
    !> Assertion source variable needs to be active variable for this lcc process
    !>
    IF (i_active < 0) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Passed active var ' // TRIM(active_var_name) &
        & // ' is not an active var in lcc structure for ' // lcc_relocations%name )
    ENDIF

    !>
    !> Assertion sink variable needs to be passive variable for this lcc process
    !>
    IF (i_passive < 0) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Passed passive var ' // TRIM(passive_var_name) &
        & // ' is not a passive var in lcc structure for ' // lcc_relocations%name )
    ENDIF

    IF (PRESENT(transfer_fract)) THEN
      this_transfer_fract = transfer_fract
    ELSE
      this_transfer_fract = 1.0
    ENDIF

    IF (PRESENT(conversion_fact)) THEN
      this_conversion_fact = conversion_fact
    ELSE
      this_conversion_fact = 1.0
    ENDIF

    ! Currently it is assumed that the passive and the active variable are of the same type (REAL2D),
    ! if relocation between REAL3D or even between REAL2D and REAL3D is required it needs to be implemented here
    IF (.NOT. ((lcc_relocations%passive_vars(i_passive)%p%var_type == REAL2D) &
        & .AND. (lcc_relocations%active_vars(i_active)%p%var_type == REAL2D))) THEN
      CALL finish(TRIM(routine), 'Transfer from active to passive var is currently only implemented for REAL2D variables,' &
          & // ' - Active variable is of type ' // int2string(lcc_relocations%active_vars(i_active)%p%var_type)            &
          & // ' - Passive variable is of type ' // int2string(lcc_relocations%passive_vars(i_passive)%p%var_type))
    END IF

    ! Add fraction of active content to be transferred to passive var,
    ! potentially converting to the unit of the passive var
    lcc_relocations%passive_vars(i_passive)%p%relocate_this(ics:ice,i_tile,iblk)       &
      & = lcc_relocations%passive_vars(i_passive)%p%relocate_this(ics:ice,i_tile,iblk) &
      & + (lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk)  &
      &    * this_transfer_fract * this_conversion_fact)

    ! And subtract from active - of course using the unit of the active var (i.e. no conversion)
    lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk)        &
      & = lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk)  &
      & - (lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk) &
      &    * this_transfer_fract)

  END SUBROUTINE transfer_from_active_to_passive_var_onChunk


  ! ====================================================================================================== !
  !>
  !> Copies (a fraction of) matter from an active (providing) to a passive (receiving) variable. In contrast
  !> to the routine [[transfer_from_active_to_passive_var_onChunk]], the matter is not transferred but copied,
  !> since in contrast here it is not subtracted from the active variable.
  !> This can be required with diagnostic variables.
  !>
  SUBROUTINE copy_from_active_to_passive_var_onChunk(lcc_relocations, i_tile, iblk, ics, ice, &
      &                         active_var_name, passive_var_name, copy_fract, conversion_fact)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc), INTENT(inout) :: lcc_relocations !< The lcc structure of the calling process
    INTEGER,              INTENT(in)    :: i_tile  !< Index of the tile in the lcc structure for which the
                                                   !< copy should be conducted
    INTEGER,              INTENT(in)    :: iblk    !< Number of current block (chunk)
    INTEGER,              INTENT(in)    :: ics     !< Index of chunk start
    INTEGER,              INTENT(in)    :: ice     !< Index of chunk end
    CHARACTER(len=*), INTENT(in)    :: active_var_name  !< Active (source) variable
    CHARACTER(len=*), INTENT(in)    :: passive_var_name !< Passive (sink) variable
    REAL(jp),  INTENT(in), OPTIONAL :: copy_fract       !< Fraction of matter copied from the active to the
                                                        !< passive variable - If not passed, 1.0 is assumed
    REAL(jp),  INTENT(in), OPTIONAL :: conversion_fact  !< Conversion factor in case that the if active
                                                        !< variable has a different unit then the passive
                                                        !< variable - If not passed, 1.0 is assumed
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':copy_from_active_to_passive_var_onChunk'

    REAL(jp) :: &
      & this_copy_fract, &      !< Fraction of the matter on the active variable which should be copied
      & this_conversion_fact    !< Conversion factor in case active and passive variables have different units
    INTEGER  :: &
      & i_active, &             !< Index of this active variable in the lcc structure
      & i_passive               !< Index of this passive variable in the lcc structure

    ! -------------------------------------------------------------------------------------------------- !

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    i_active = one_of( active_var_name , lcc_relocations%active_vars_names)
    i_passive = one_of( passive_var_name , lcc_relocations%passive_vars_names)

    !>
    !> Assertion source variable needs to be active variable for this lcc process
    !>
    IF (i_active < 0) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Passed active var ' // TRIM(active_var_name) &
        & // ' is not an active var in lcc structure for ' // lcc_relocations%name )
    ENDIF

    !>
    !> Assertion sink variable needs to be passive variable for this lcc process
    !>
    IF (i_passive < 0) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Passed passive var ' // TRIM(passive_var_name) &
        & // ' is not a passive var in lcc structure for ' // lcc_relocations%name )
    ENDIF

    IF (PRESENT(copy_fract)) THEN
      this_copy_fract = copy_fract
    ELSE
      this_copy_fract = 1.0
    ENDIF

    IF (PRESENT(conversion_fact)) THEN
      this_conversion_fact = conversion_fact
    ELSE
      this_conversion_fact = 1.0
    ENDIF

    ! Currently it is assumed that the passive and the active variable are of the same type (REAL2D),
    ! if a copy between REAL3D or even between REAL2D and REAL3D is required it needs to be implemented here
    IF (.NOT. ((lcc_relocations%passive_vars(i_passive)%p%var_type == REAL2D) &
        & .AND. (lcc_relocations%active_vars(i_active)%p%var_type == REAL2D))) THEN
      CALL finish(TRIM(routine), 'Copy from active to passive var is currently only implemented for REAL2D variables,' &
          & // ' - Active variable is of type ' // int2string(lcc_relocations%active_vars(i_active)%p%var_type)        &
          & // ' - Passive variable is of type ' // int2string(lcc_relocations%passive_vars(i_passive)%p%var_type))
    END IF

    ! Add active content to passive var
    lcc_relocations%passive_vars(i_passive)%p%relocate_this(ics:ice,i_tile,iblk) &
      & = lcc_relocations%passive_vars(i_passive)%p%relocate_this(ics:ice,i_tile,iblk) &
      & + (lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk) * this_copy_fract * this_conversion_fact)

  END SUBROUTINE copy_from_active_to_passive_var_onChunk

  ! ====================================================================================================== !
  !> This routine sums the amount of matter provided by all active (providing) variables.
  !>
  SUBROUTINE collect_matter_of_active_vars_onChunk(collected_matter, lcc_relocations, &
      &                                       i_tile, iblk, ics, ice, active_var_names)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    REAL(jp),                   INTENT(INOUT) :: collected_matter(:) !< Array for the collected matter
    TYPE(t_jsb_lcc_proc),       INTENT(INOUT) :: lcc_relocations !< Lcc structure of the lcc process
    INTEGER,                    INTENT(IN)    :: i_tile  !< Index of the tile in lcc structure
    INTEGER,                    INTENT(in)    :: iblk    !< Number of current block (chunk)
    INTEGER,                    INTENT(in)    :: ics     !< Index of chunk start
    INTEGER,                    INTENT(in)    :: ice     !< Index of chunk end
    CHARACTER(VARNAME_LEN),     INTENT(IN)    :: active_var_names(:) !< List of the active (source) variables
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':collect_matter_of_active_vars_onChunk'
    INTEGER ::   &
      & i,       & !< Looping index for active variables in the passed array of variable names
      & i_active   !< Index of a certain active variable in the lcc structure
    ! -------------------------------------------------------------------------------------------------- !

    DO i = 1,SIZE(active_var_names)
      i_active = one_of( TRIM(active_var_names(i)) , lcc_relocations%active_vars_names)
      IF (i_active > 0) THEN
        collected_matter(:) = collected_matter(:) + lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk)
        lcc_relocations%active_vars(i_active)%p%relocate_this(ics:ice,i_tile,iblk) = 0.0_jp
      END IF
    END DO

  END SUBROUTINE collect_matter_of_active_vars_onChunk

#endif
END MODULE mo_jsb_lcc_class
