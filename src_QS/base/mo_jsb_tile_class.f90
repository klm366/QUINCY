!> Contains abstract tile class that contains the surface structure and memory state
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
!>### Definition of abstract class and methods for tiles
!> The most important components of a tile instance are
!>
!> - the fraction of a grid box that the tile occupies
!> - the memories of all processes (i.e. variables) owned by this tile
!> - a list of LCT (land cover type) id's that are represented on this tile
!> - a store pointing to the element variables of all BGC materials across all processes
!> - a list of actions to be applied on this tile for the different processes
!> - a list of aggregators used to aggregate variables of the child tiles to this tile

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_jsb_tile_class
#ifndef __NO_JSBACH__

  USE mo_jsb_control,         ONLY: debug_on, acc_stream, acc_async_sync
  USE mo_kind,                ONLY: wp, dp
  USE mo_exception,           ONLY: message, finish
  USE mo_io_units,            ONLY: filename_max
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN

  USE mo_hsm_class,           ONLY: t_Hsm, t_State, t_Message
  USE mo_jsb_memory_class,    ONLY: t_jsb_memory_p, t_jsb_memory
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct, &
    &                               LAND_TYPE, VEG_TYPE, BARE_TYPE, GLACIER_TYPE, LAKE_TYPE
  USE mo_jsb_var_class,       ONLY: t_jsb_var_p, t_jsb_var, t_jsb_var_real2d, t_jsb_var_real3d
  ! USE mo_jsb_cqt_class,       ONLY: t_jsb_consQuan_p
  ! USE mo_jsb_lcc_class,       ONLY: t_jsb_lcc_proc_p
#ifndef __NO_QUINCY__
  USE mo_lnd_bgcm_store,      ONLY: t_lnd_bgcm_store
#endif


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile_abstract
  PUBLIC :: t_jsb_aggregator, t_jsb_aggregator_p, t_jsb_aggregator_weighted_by_fract


!!$  INTEGER, PARAMETER :: max_depth_of_tile_tree = 3

  ! Declarations for aggregators
  TYPE, ABSTRACT :: t_jsb_aggregator
    CHARACTER(len=25) :: name    !< Aggregators are distinguished by using unique names
  CONTAINS
    ! Begin/End routines enable multi-stream OpenACC aggregation
    PROCEDURE(begin_interface),        DEFERRED             :: BeginAggregate !< Start  multi-stream aggregation (useful for OpenACC)
    PROCEDURE(end_interface),          DEFERRED             :: EndAggregate   !< Finish multi-stream aggregation (useful for OpenACC)
    PROCEDURE(aggregate_interface_2d), DEFERRED             :: Aggregate_2d !< Aggregate operator for 2d variables
    PROCEDURE(aggregate_interface_3d), DEFERRED             :: Aggregate_3d !< Aggregate operator for 3d variables
    GENERIC                                                 :: Aggregate => Aggregate_2d, Aggregate_3d
  END TYPE t_jsb_aggregator

  TYPE t_jsb_aggregator_p
    CLASS(t_jsb_aggregator), POINTER :: p => NULL()
  END TYPE t_jsb_aggregator_p

  ! Declarations for tile
  TYPE, ABSTRACT, EXTENDS(t_State) :: t_jsb_tile_abstract
    TYPE(t_jsb_lct),        ALLOCATABLE    :: lcts(:)              !< List of land cover types for tile
                                                                   !< The first entry is always the primary LCT for this tile;
                                                                   !< the remaining lcts are those inherited for aggregation
    ! Convinience settings to easily check primary LCT types and those of the child tiles
    LOGICAL                                :: is_vegetation = .FALSE.         !< whether primary lct is VEG_TYPE
    LOGICAL                                :: is_bare = .FALSE.               !< whether primary lct is BARE_TYPE
    LOGICAL                                :: is_land = .FALSE.               !< whether primary lct is LAND_TYPE
    LOGICAL                                :: is_glacier = .FALSE.            !< whether primary lct is GLACIER_TYPE
    LOGICAL                                :: is_lake = .FALSE.               !< whether primary lct is LAKE_TYPE
    LOGICAL                                :: contains_vegetation = .FALSE.   !< whether any lct of tile is VEG_TYPE
    LOGICAL                                :: contains_bare = .FALSE.         !< whether any lct of tile is BARE_TYPE
    LOGICAL                                :: contains_land = .FALSE.         !< whether any lct of tile is LAND_TYPE
    LOGICAL                                :: contains_glacier = .FALSE.      !< whether any lct of tile is GLACIER_TYPE
    LOGICAL                                :: contains_lake = .FALSE.         !< whether any lct of tile is LAKE_TYPE
    LOGICAL                                :: contains_soil = .FALSE.         !< whether any lct of tile is one that has soil
    !
    INTEGER,                ALLOCATABLE    :: process_action(:)
    REAL(wp),                     POINTER  :: fract     (:,:)          !< Fraction of this tile (relative to grid box)
    REAL(wp),                     POINTER  :: fract_old (:,:)          !< Old fractions of this tile (relative to grid box)
    CHARACTER(LEN=filename_max)            :: fract_filename           !< Name of file to read tile fractions from. The variable
                                                                       !< names must be either "fract[_max]_"//tile%name,
                                                                       !< or specified by fract_varname (in the use case).
    CHARACTER(LEN=20)                      :: fract_varname                 !< Name of variable containing fractions for this tile
    LOGICAL,                ALLOCATABLE    :: l_fract_children_changed(:,:)
    TYPE(t_jsb_aggregator_p), ALLOCATABLE  :: aggregators(:)
    TYPE(t_lnd_bgcm_store), POINTER        :: bgcm_store => NULL()    !< the bgcm store collecting pools and their elements on this tile
    !CLASS(t_jsb_param_p),       POINTER    :: param(:)
    CLASS(t_jsb_memory_p),      POINTER    :: mem(:)
    ! The next three pointers basically repeat what's already in the base type, but they save time-consuming "SELECT TYPE"s
    CLASS(t_jsb_tile_abstract), POINTER    :: parent_tile => NULL()
    CLASS(t_jsb_tile_abstract), POINTER    :: first_child_tile => NULL()
    CLASS(t_jsb_tile_abstract), POINTER    :: next_sibling_tile => NULL()
    !
    INTEGER                                :: grid_id = -1
    INTEGER                                :: owner_model_id = 0
  CONTAINS
    PROCEDURE                              :: Add_lct
    PROCEDURE :: Set_process_action => Set_process_action_tile
    PROCEDURE :: Set_process_action_parents
    PROCEDURE                              :: Register_aggregator
    PROCEDURE                              :: Get_aggregator
    PROCEDURE                              :: Get_fraction_1d
    PROCEDURE                              :: Get_fraction_2d
    GENERIC                                :: Get_fraction => Get_fraction_1d, Get_fraction_2d
    PROCEDURE(Init_interface),           DEFERRED :: Init
    PROCEDURE(Init_vars_interface),      DEFERRED :: Count_and_classify_bgc_materials
    PROCEDURE(Init_vars_interface),      DEFERRED :: Collect_bgc_materials
    PROCEDURE(Check_on_process),         DEFERRED :: Is_process_active
    PROCEDURE(Check_on_process),         DEFERRED :: Is_process_calculated
    PROCEDURE(Check_on_process),         DEFERRED :: Has_process_memory
  END type t_jsb_tile_abstract

  !< Deferred subroutine interfaces (must be provided as part of t_jsb_tile_abstract).
  !< The subroutines themselves are located in the jsb_tile module.

  ABSTRACT INTERFACE
    SUBROUTINE Init_interface(this, varlist_name, prefix, suffix, grid_id, in_var_groups)
      IMPORT :: t_jsb_tile_abstract
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
      CHARACTER(len=*),           INTENT(in)    :: varlist_name
      CHARACTER(len=*),           INTENT(in)    :: prefix
      CHARACTER(len=*),           INTENT(in)    :: suffix
      INTEGER,                    INTENT(in)    :: grid_id
      LOGICAL,                    INTENT(in)    :: in_var_groups
    END SUBROUTINE Init_interface
  END INTERFACE

  ! ABSTRACT INTERFACE
  !   SUBROUTINE Init_fractions_interface(this, varlist_name, prefix, suffix, l_fixed_fractions, l_rel_fractions)
  !     IMPORT :: t_jsb_tile_abstract
  !     CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
  !     CHARACTER(len=*),           INTENT(in)    :: varlist_name
  !     CHARACTER(len=*),           INTENT(in)    :: prefix
  !     CHARACTER(len=*),           INTENT(in)    :: suffix
  !     LOGICAL,                    INTENT(in)    :: l_fixed_fractions
  !     LOGICAL,                    INTENT(in)    :: l_rel_fractions
  !   END SUBROUTINE Init_fractions_interface
  ! END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE Init_vars_interface(this)
      IMPORT :: t_jsb_tile_abstract
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    END SUBROUTINE Init_vars_interface
  END INTERFACE

  ABSTRACT INTERFACE
    LOGICAL FUNCTION Check_on_process(this, iproc)
      IMPORT :: t_jsb_tile_abstract
      CLASS(t_jsb_tile_abstract), INTENT(in) :: this
      INTEGER, INTENT(in) :: iproc
    END FUNCTION Check_on_process
  END INTERFACE
  
  ! ABSTRACT INTERFACE
  !   FUNCTION Handler_interface(this, msg_in) RESULT(return_ptr)
  !     IMPORT :: t_jsb_tile_abstract, t_Hsm, t_Message
  !     CLASS(t_jsb_tile_abstract),  INTENT(inout) :: this
  !     CLASS(t_Message),            INTENT(in)    :: msg_in
  !     CLASS(t_Message),            POINTER       :: return_ptr
  !   END FUNCTION Handler_interface
  ! END INTERFACE
  ! 
  ! ABSTRACT INTERFACE
  !   SUBROUTINE Print_interface(this)
  !     IMPORT :: t_jsb_tile_abstract
  !     CLASS(t_jsb_tile_abstract), INTENT(in) :: this
  !   END SUBROUTINE Print_interface
  ! END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE begin_interface(this)
      IMPORT :: t_jsb_aggregator
      CLASS(t_jsb_aggregator),    INTENT(inout), TARGET :: this
    END SUBROUTINE begin_interface
    SUBROUTINE end_interface(this)
      IMPORT :: t_jsb_aggregator
      CLASS(t_jsb_aggregator),    INTENT(inout), TARGET :: this
    END SUBROUTINE end_interface

    !> Interface for deferred subroutine `Aggregate_2d` of `t_jsb_aggregator`
    SUBROUTINE aggregate_interface_2d(this, tile, var, ics, ice, iblk, caller)
      IMPORT :: t_jsb_aggregator, t_jsb_tile_abstract, t_jsb_var_real2d
      CLASS(t_jsb_aggregator),    INTENT(inout) :: this
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
      TYPE(t_jsb_var_real2d),     INTENT(inout) :: var
      INTEGER,                    INTENT(in)    :: ics, ice, iblk
      CHARACTER(len=*),           INTENT(in)    :: caller
    END SUBROUTINE aggregate_interface_2d
    SUBROUTINE aggregate_interface_3d(this, tile, var, ics, ice, iblk,caller)
      IMPORT :: t_jsb_aggregator, t_jsb_tile_abstract, t_jsb_var_real3d
      CLASS(t_jsb_aggregator),    INTENT(inout) :: this
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
      TYPE(t_jsb_var_real3d),     INTENT(inout) :: var
      INTEGER,                    INTENT(in)    :: ics, ice, iblk
      CHARACTER(len=*),           INTENT(in)    :: caller
    END SUBROUTINE aggregate_interface_3d
    SUBROUTINE aggregate_interface_setopt_logical(this, name, value)
      IMPORT :: t_jsb_aggregator
      CLASS(t_jsb_aggregator), INTENT(inout) :: this
      CHARACTER(len=*)       , INTENT(in)    :: name
      LOGICAL,                 INTENT(in)    :: value
    END SUBROUTINE aggregate_interface_setopt_logical
  END INTERFACE

  TYPE, EXTENDS(t_jsb_aggregator) :: t_jsb_aggregator_weighted_by_fract
    REAL(wp), ALLOCATABLE :: fractions(:,:,:) ! R: (grid%nproma, grid%nblks, no_children))
    REAL(wp), ALLOCATABLE :: fract_sum(:,:)
    REAL(wp), ALLOCATABLE :: data_sum(:,:)

    INTEGER :: stream, original_stream !< OpenACC variables for multi-stream aggregation
  CONTAINS
    !> Functions to make use of the aggregator
    PROCEDURE :: BeginAggregate   => BeginAggregate_weighted_by_fract !< Start  multi-stream aggregation (useful for OpenACC)
    PROCEDURE :: EndAggregate     => EndAggregate_weighted_by_fract   !< Finish multi-stream aggregation (useful for OpenACC)
    PROCEDURE :: Set_fractions_1d => Set_fractions_aggregator_1d  !< Set fraction field for chunk. TODO: Rename "1d" to "chunk"
    PROCEDURE :: Set_fractions_2d => Set_fractions_aggregator_2d  !< Set fraction field for entire domain. TODO: "2d"=>"domain"
    GENERIC   :: Set_fractions    => Set_fractions_1d, Set_fractions_2d !< Generic function to set tile fractions
    PROCEDURE :: Aggregate_2d     => Aggregate_weighted_by_fract_2d     !< Specific implementation of agg. of 2d vars
    PROCEDURE :: Aggregate_3d     => Aggregate_weighted_by_fract_3d     !< Specific implementation of agg. of 3d vars
    FINAL     :: Finalize_aggregator_weighted_by_fract            !< Destructor to deallocate fields when aggregator is deallocated
  END TYPE t_jsb_aggregator_weighted_by_fract

  INTERFACE t_jsb_aggregator_weighted_by_fract
    PROCEDURE Create_aggregator_weighted_by_fract                 !< Constructor of aggregator (call to initialize)
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_tile_class'

CONTAINS

  SUBROUTINE Set_process_action_tile(this, iproc, action, action_to_parents)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    INTEGER,                    INTENT(in)    :: iproc
    INTEGER,                    INTENT(in)    :: action
    INTEGER, OPTIONAL,          INTENT(in)    :: action_to_parents

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_process_action_tile'

    this%process_action(iproc) = action

    ! Propagate action action_to_parents to all parents of this tile
    IF (PRESENT(action_to_parents)) CALL this%Set_process_action_parents(iproc, action_to_parents)

  END SUBROUTINE Set_process_action_tile

  SUBROUTINE Set_process_action_parents(this, iproc, action)

    CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: this
    INTEGER,                    INTENT(in)            :: iproc
    INTEGER,                    INTENT(in)            :: action

    CLASS(t_jsb_tile_abstract), POINTER :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_process_action_parents'

    tile => this
    DO WHILE (ASSOCIATED(tile%parent))
      SELECT TYPE (parent=>tile%parent)
      CLASS IS (t_jsb_tile_abstract)
        parent%process_action(iproc) = action
        tile => parent
      END SELECT
    END DO

  END SUBROUTINE Set_process_action_parents

  ! ====================================================================================================== !
  !
  !> adds given lct to this tile and sets the "is" and/or "contains" properties on the tile accordingly
  !
  SUBROUTINE Add_lct(this, lct)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this  !< this tile
    TYPE(t_jsb_lct),            INTENT(in)    :: lct   !< the landcover type to add for this tile
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lct), ALLOCATABLE :: temp_lcts(:)
    INTEGER :: ilct, nlct

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_lct'
    ! -------------------------------------------------------------------------------------------------- !
    IF (.NOT. ALLOCATED(this%lcts)) THEN

      ALLOCATE(this%lcts(1))
      this%lcts(1) = lct

      ! the "is" properties are only set for the primary lct
      IF (ASSOCIATED(this%parent)) THEN    ! Root tile has no primary LCT
        SELECT CASE (lct%id)
        CASE (VEG_TYPE)
          this%is_vegetation = .TRUE.
          this%is_land       = .TRUE.
        CASE (BARE_TYPE)
          this%is_bare = .TRUE.
        ! CASE (LAND_TYPE)
        !   this%is_land = .TRUE.
        CASE (GLACIER_TYPE)
          this%is_glacier = .TRUE.
        CASE (LAKE_TYPE)
          this%is_lake = .TRUE.
        END SELECT
      END IF

    ELSE

      nlct = SIZE(this%lcts)
      DO ilct=1,nlct
        ! If this lct is already registered in the tile%lcts, do nothing and return
        IF (this%lcts(ilct)%id == lct%id) RETURN
      END DO
      ALLOCATE(temp_lcts(nlct+1))
      temp_lcts(1:nlct) = this%lcts
      temp_lcts(nlct+1) = lct
      CALL MOVE_ALLOC(temp_lcts, this%lcts)

    END IF !(.NOT. ALLOCATED(this%lcts))

    ! the "contains" properties are set according to all lcts of this tile
    SELECT CASE (lct%id)
    CASE (VEG_TYPE)
      this%contains_vegetation = .TRUE.
      this%contains_land       = .TRUE.
    CASE (BARE_TYPE)
      this%contains_bare = .TRUE.
    ! CASE (LAND_TYPE)
    !   this%contains_land = .TRUE.
    CASE (GLACIER_TYPE)
      this%contains_glacier = .TRUE.
    CASE (LAKE_TYPE)
      this%contains_lake = .TRUE.
    END SELECT

    this%contains_soil       = this%contains_vegetation .OR. this%contains_bare .OR. this%contains_land

  END SUBROUTINE Add_lct

  SUBROUTINE Get_fraction_1d(this, ics, ice, iblk, fract, fract_old)

    CLASS(t_jsb_tile_abstract), INTENT(in)    :: this
    INTEGER,                    INTENT(in)    :: ics, ice, iblk
    REAL(wp), OPTIONAL,         INTENT(inout) :: &
      & fract(:),        &
      & fract_old(:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_fraction_1d'

    IF (PRESENT(fract)) THEN
       IF (SIZE(fract) /= ice-ics+1) CALL finish(TRIM(routine), 'Dimension mismatch')
       fract(:) = 1.0_wp
    !   fract(:) = this%fract(ics:ice, iblk)
    END IF

    IF (PRESENT(fract_old)) THEN
      IF (SIZE(fract_old) /= ice-ics+1) CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_old(:) = this%fract_old(ics:ice, iblk)
    END IF

  END SUBROUTINE Get_fraction_1d

  SUBROUTINE Get_fraction_2d(this, fract, fract_old)

    CLASS(t_jsb_tile_abstract), INTENT(in)    :: this
    REAL(wp), OPTIONAL,         INTENT(inout) :: &
      & fract(:,:),        &
      & fract_old(:,:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_fraction_2d'

    IF (PRESENT(fract)) THEN
      IF (SIZE(fract,1) /= SIZE(this%fract,1) .OR. SIZE(fract,2) /= SIZE(this%fract,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      fract(:,:) = this%fract(:,:)
    END IF

    IF (PRESENT(fract_old)) THEN
      IF (SIZE(fract_old,1) /= SIZE(this%fract_old,1) .OR. SIZE(fract_old,2) /= SIZE(this%fract_old,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_old(:,:) = this%fract_old(:,:)
    END IF

  END SUBROUTINE Get_fraction_2d

  SUBROUTINE Register_aggregator(this, aggregator)

    CLASS(t_jsb_tile_abstract), INTENT(inout)        :: this
    CLASS(t_jsb_aggregator),    INTENT(in),   TARGET :: aggregator

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Register_aggregator'


    IF (.NOT. ALLOCATED(this%aggregators)) THEN
      this%aggregators = (/ t_jsb_aggregator_p(aggregator) /)
    ELSE
      this%aggregators = (/ this%aggregators, t_jsb_aggregator_p(aggregator) /)
    END IF

    IF (debug_on()) CALL message(TRIM(routine), 'Registered aggregator '//TRIM(aggregator%name)//' for tile '//TRIM(this%name))

  END SUBROUTINE Register_aggregator

  FUNCTION Get_aggregator(this, name, copy) RESULT(aggregator)

    CLASS(t_jsb_tile_abstract), INTENT(in) :: this
    CHARACTER(len=*),           INTENT(in) :: name
    LOGICAL, OPTIONAL,          INTENT(in) :: copy
    CLASS(t_jsb_aggregator),    POINTER    :: aggregator

    INTEGER :: i
    LOGICAL :: l_copy

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_aggregator'

    l_copy = .FALSE.
    IF (PRESENT(copy)) l_copy = copy

    aggregator => NULL()
    IF (ALLOCATED(this%aggregators)) THEN
      DO i=1,SIZE(this%aggregators)
        IF (TRIM(this%aggregators(i)%p%name) == TRIM(name)) THEN
          IF (l_copy) THEN
            ALLOCATE(aggregator, source=this%aggregators(i)%p)
          ELSE
            aggregator => this%aggregators(i)%p
          END IF
          EXIT
        END IF
      END DO
    END IF

    IF (.NOT.(ASSOCIATED(aggregator))) &
      & CALL finish(TRIM(routine), 'Aggregator "'//TRIM(name)//'" not found on tile '//TRIM(this%name))

  END FUNCTION Get_aggregator

  FUNCTION Create_aggregator_weighted_by_fract(tile) RESULT(aggregator)

    CLASS(t_jsb_tile_abstract), INTENT(inout)          :: tile
    CLASS(t_jsb_aggregator_weighted_by_fract), POINTER :: aggregator

    INTEGER                   :: no_children

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Create_aggregator_weighted_by_fract'

    ALLOCATE(t_jsb_aggregator_weighted_by_fract::aggregator)
    aggregator%name = 'weighted_by_fract'

    IF (debug_on()) CALL message(TRIM(routine), 'Created aggregator '//TRIM(aggregator%name)//' on tile '//TRIM(tile%name))

  END FUNCTION Create_aggregator_weighted_by_fract

  SUBROUTINE Finalize_aggregator_weighted_by_fract(aggregator)

    TYPE(t_jsb_aggregator_weighted_by_fract) :: aggregator

    CHARACTER(len=*), PARAMETER :: routine = modname//':Finalize_aggregator_weighted_by_fract'

    !CALL message(TRIM(routine), 'Finalizing aggregator '//TRIM(aggregator%name))

    IF (ALLOCATED(aggregator%fractions)) THEN
      !$acc exit data delete(aggregator%fractions)
      DEALLOCATE(aggregator%fractions)
    END IF
    IF (ALLOCATED(aggregator%data_sum)) THEN
      !$acc exit data delete(aggregator%data_sum)
      DEALLOCATE(aggregator%data_sum)
    END IF
    IF (ALLOCATED(aggregator%fract_sum)) THEN
      !$acc exit data delete(aggregator%fract_sum)
      DEALLOCATE(aggregator%fract_sum)
    END IF

  END SUBROUTINE Finalize_aggregator_weighted_by_fract

  SUBROUTINE BeginAggregate_weighted_by_fract(this)
    CLASS(t_jsb_aggregator_weighted_by_fract),    INTENT(inout), TARGET :: this

    IF (this%stream /= acc_async_sync) THEN
      this%original_stream = acc_stream
      this%stream = this%original_stream + 1
    END IF
  END SUBROUTINE BeginAggregate_weighted_by_fract
  SUBROUTINE EndAggregate_weighted_by_fract(this)
    CLASS(t_jsb_aggregator_weighted_by_fract),    INTENT(inout), TARGET :: this
    INTEGER :: stream, dummy

    IF (this%original_stream /= -1) THEN
      DO stream = this%original_stream + 1, this%stream - 1
        !$ACC WAIT(stream) ASYNC(this%original_stream)
      END DO
      ! Dummy kernel for simpler graph dependency resolution
      !$ACC PARALLEL ASYNC(this%original_stream)
      dummy = 1
      !$ACC END PARALLEL
    END IF
    this%original_stream = -1
    this%stream = acc_stream
  END SUBROUTINE EndAggregate_weighted_by_fract

  ! ====================================================================================================== !
  !
  !> Updates the part of the fraction field of an aggregator belonging to the specified chunk for all child tiles
  !
  ! TODO: Change name to Set_fractions_aggregator_chunk
  SUBROUTINE Set_fractions_aggregator_1d(this, tile, ics, ice, iblk)

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile
    INTEGER,                                   INTENT(in)    :: ics, ice, iblk

    INTEGER :: no_children, i
    CLASS(t_jsb_tile_abstract), POINTER :: current
    ! REAL(wp), ALLOCATABLE :: fract(:), fract_sum(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_fractions_aggregator_1d'

    ! Check if child tile fractions have changed and, if necessary, update fractions of aggregator
    ! by looping over children of tile
    IF (ANY(tile%l_fract_children_changed(ics:ice,iblk))) THEN
      no_children = tile%Get_no_of_children()
      current => tile%first_child_tile
      DO i=1,no_children
        CALL current%Get_fraction(ics, ice, iblk, fract=this%fractions(ics:ice,iblk,i))
        current => current%next_sibling_tile
      END DO

      !$acc update device( this%fractions(ics:ice,iblk,:) )

      ! ALLOCATE(fract(ice-ics+1))
      ! ALLOCATE(fract_sum(ice-ics+1))
      ! fract(:) = tile%fract(ics:ice,iblk)
      ! fract_sum(:) = SUM(this%fractions(ics:ice,iblk,:), DIM=2)
      ! ! IF (ASSOCIATED(tile%parent) .AND. ANY( fract_sum(:) < 1._wp - EPSILON(1._wp) .AND. fract(:) > 0._wp)) THEN
      ! !   CALL message(TRIM(routine), 'Sum of child tile fractions on tile '//TRIM(tile%name)//': '//&
      ! !     &                         real2string(MINVAL(fract_sum(:)))//' - '//real2string(MAXVAL(fract_sum(:))), &
      ! !     &          all_print=.TRUE.)
      ! !   CALL finish(TRIM(routine), 'Sum of tile fractions is not 1 on tile')
      ! ! END IF
      ! DEALLOCATE(fract, fract_sum)

      tile%l_fract_children_changed(ics:ice,iblk) = .FALSE.
    END IF

  END SUBROUTINE Set_fractions_aggregator_1d

  SUBROUTINE Set_fractions_aggregator_2d(this, tile)

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile

    INTEGER :: no_children, i
    CLASS(t_jsb_tile_abstract), POINTER :: current
    ! REAL(wp), ALLOCATABLE :: fract(:), fract_sum(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_fractions_aggregator_2d'

    ! print*,'EEE ',iblk,ics,ice,ANY(tile%l_fract_children_changed(ics:ice,iblk))
    ! Check if child tile fractions have changed and, if necessary, update fractions of aggregator
    ! by looping over children of tile
    IF (ANY(tile%l_fract_children_changed(:,:))) THEN
      no_children = tile%Get_no_of_children()
      current => tile%first_child_tile
      DO i=1,no_children
        CALL current%Get_fraction_2d(fract=this%fractions(:,:,i))
        current => current%next_sibling_tile
      END DO

      !$acc update device( this%fractions(:,:,:) )

      ! ALLOCATE(fract(ice-ics+1))
      ! ALLOCATE(fract_sum(ice-ics+1))
      ! fract(:) = tile%fract(ics:ice,iblk)
      ! fract_sum(:) = SUM(this%fractions(ics:ice,iblk,:), DIM=2)
      ! ! IF (ASSOCIATED(tile%parent) .AND. ANY( fract_sum(:) < 1._wp - EPSILON(1._wp) .AND. fract(:) > 0._wp)) THEN
      ! !   CALL message(TRIM(routine), 'Sum of child tile fractions on tile '//TRIM(tile%name)//': '//&
      ! !     &                         real2string(MINVAL(fract_sum(:)))//' - '//real2string(MAXVAL(fract_sum(:))), &
      ! !     &          all_print=.TRUE.)
      ! !   CALL finish(TRIM(routine), 'Sum of tile fractions is not 1 on tile')
      ! ! END IF
      ! DEALLOCATE(fract, fract_sum)

      tile%l_fract_children_changed(:,:) = .FALSE.
    END IF

  END SUBROUTINE Set_fractions_aggregator_2d

  SUBROUTINE Aggregate_weighted_by_fract_2d(this, tile, var, ics, ice, iblk, caller)

    USE mo_hsm_class,      ONLY: t_State

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile
    TYPE(t_jsb_var_real2d),                    INTENT(inout) :: var
    INTEGER,                                   INTENT(in)    :: ics, ice, iblk
    CHARACTER(len=*),                          INTENT(in)    :: caller

    INTEGER :: no_children, i, j, iproc, nc
    ! CLASS(t_jsb_tile_abstract), POINTER :: current
    ! TYPE(t_jsb_var_real2d), POINTER :: var_tmp
    CLASS(t_jsb_var), POINTER            :: var_tmp
    REAL(wp), ALLOCATABLE :: data(:,:), fractions(:,:)    ! Collected data and fractions from children
    REAL(wp)              :: fract_sum(ice-ics+1)
    REAL(dp) :: missval

    CHARACTER(len=256) :: routine

    IF (debug_on()) THEN
      routine = TRIM(modname)//':Aggregate_weighted_by_fract_2d'//'(called by '//TRIM(caller)//')'
    ELSE
      routine = TRIM(modname)//':Aggregate_weighted_by_fract_2d'
    END IF

    IF (.NOT. ASSOCIATED(var%ptr)) RETURN

    nc = ice - ics + 1

    missval = var%missval

    iproc = var%owner_proc_id
    IF (iproc < 0) CALL finish(TRIM(routine), 'Unknown process for variable '//TRIM(var%name)//' on tile '//TRIM(tile%name))

    no_children = tile%Get_no_of_children()
    if (no_children == 0) return

    IF (.NOT. ALLOCATED(var%child_idx)) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' not allocated')
    ELSE IF (SIZE(var%child_idx) /= no_children) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' has wrong size')
    END IF

    ALLOCATE(data     (nc,no_children))
    ALLOCATE(fractions(nc,no_children))

    j = 0
    fractions(:,:) = 0._wp
    fract_sum(:)   = 0._wp
    data     (:,:) = 0._wp

    ! current => tile%first_child_tile
    DO i=1,no_children
      ! IF (ASSOCIATED(current%mem(iproc)%p)) THEN  ! Process is active on child tile
      !   SELECT TYPE (v => current%mem(iproc)%p%Get_var_real2d(TRIM(var%name)))
      !   TYPE IS (t_jsb_var_real2d)
      !     var_tmp => v
      !   END SELECT
      ! IF (var%child_idx(i) > 0) THEN
      !   var_tmp => tile%mem(iproc)%p%children(i)%p%vars(var%child_idx(i))%p
      ! ELSE
      !   var_tmp => NULL()
      ! END IF
      ! IF (ASSOCIATED(var_tmp)) THEN
      IF (var%child_idx(i) > 0 .AND. ASSOCIATED(tile%mem(iproc)%p)) THEN
        IF (.NOT. ASSOCIATED(tile%mem(iproc)%p%children(i)%p%vars(var%child_idx(i))%p)) THEN
          CALL finish(TRIM(routine), 'ERROR - this should not happen!')
        ELSE
          var_tmp => tile%mem(iproc)%p%children(i)%p%vars(var%child_idx(i))%p
        END IF
        fractions(:,i) = this%fractions(ics:ice,iblk,i)
        IF (.NOT. var%l_aggregate_all) fract_sum(:)   = fract_sum(:) + fractions(:,i)
        data(:,i) = var_tmp%field(ics:ice,iblk,1,1,1)
        ! data(:,i) = var_tmp%ptr(ics:ice,iblk)
        j = j + 1
      END IF
      ! END IF
      ! current => current%next_sibling_tile
    END DO
    j=no_children
    IF (j == 0) CALL finish(TRIM(routine), 'Variable '//TRIM(var%name)//' not found in any of children of tile '//&
      &                                    TRIM(tile%name)//' - this should not happen!')

    ! IF (ANY(fract_sum(:) > 1._wp + EPSILON(1._wp))) THEN
    !   ! print*, fract_sum
    !   CALL message(TRIM(routine), 'On tile '//TRIM(tile%name)//': fract_sum > 1 - '// &
    !     &                         real2string(MINVAL(fract_sum(:)))//' - '//real2string(MAXVAL(fract_sum(:))), &
    !     &          all_print=.TRUE.)
    !   CALL finish(TRIM(routine),  '    ... this should not happen')
    ! END IF

    IF (j == 1 .AND. .NOT. var%l_aggregate_all) THEN
      ! Only one tile to aggregate, avoid multiplying and dividing by a potentially small fraction 
      var%ptr(ics:ice,iblk) = SUM(data(:,:), DIM=2)
    ELSE IF (var%l_aggregate_all) THEN              ! Sum of fractions is 1
      ! For tiles where variable is not defined data(:,tile) is 0
      var%ptr(ics:ice,iblk) = SUM(fractions(:,:) * data(:,:), DIM=2)
    ELSE
      WHERE (fract_sum(:) > 0._wp)
        var%ptr(ics:ice,iblk) = SUM(fractions(:,:) * data(:,:), DIM=2) / fract_sum(:)
      !ELSE WHERE
      !  var%ptr(ics:ice,iblk) = missval
      END WHERE
    END IF

    DEALLOCATE(data, fractions)

    IF (var%l_aggregate_all) THEN
      IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//'; all_tiles)')
    ELSE
      IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//')')
    END IF

  END SUBROUTINE Aggregate_weighted_by_fract_2d

  SUBROUTINE Aggregate_weighted_by_fract_3d(this, tile, var, ics, ice, iblk, caller)

    USE mo_hsm_class,      ONLY: t_State

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile
    TYPE(t_jsb_var_real3d),                    INTENT(inout) :: var
    INTEGER,                                   INTENT(in)    :: ics, ice, iblk
    CHARACTER(len=*),                          INTENT(in)    :: caller

    INTEGER :: no_children, i, j, iproc, nc, ilev
    ! CLASS(t_jsb_tile_abstract), POINTER :: current
    ! TYPE(t_jsb_var_real3d), POINTER :: var_tmp
    CLASS(t_jsb_var), POINTER :: var_tmp
    REAL(wp), ALLOCATABLE :: data(:,:,:), fractions(:,:)    ! Collected data and fractions from children
    REAL(wp)              :: fract_sum(ice-ics+1)
    REAL(dp) :: missval

    CHARACTER(len=256) :: routine

    IF (debug_on()) THEN
      routine = TRIM(modname)//':Aggregate_weighted_by_fract_2d'//'(called by '//TRIM(caller)//')'
    ELSE
      routine = TRIM(modname)//':Aggregate_weighted_by_fract_2d'
    END IF

    IF (.NOT. ASSOCIATED(var%ptr)) RETURN

    CALL this%Set_fractions(tile, ics, ice, iblk)

    nc = ice - ics + 1

    missval = var%missval

    iproc = var%owner_proc_id
    IF (iproc < 0) CALL finish(TRIM(routine), 'Unknown process for variable '//TRIM(var%name)//' on tile '//TRIM(tile%name))

    no_children = tile%Get_no_of_children()

    IF (.NOT. ALLOCATED(var%child_idx)) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' not allocated')
    ELSE IF (SIZE(var%child_idx) /= no_children) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' has wrong size')
    END IF

    ALLOCATE(data     (nc, SIZE(var%ptr,2), no_children))
    ALLOCATE(fractions(nc,                  no_children))

    j = 0
    fractions(:,:)   = 0._wp
    fract_sum(:  )   = 0._wp
    data     (:,:,:) = 0._wp
    ! current => tile%first_child_tile
    DO i=1,no_children
      ! IF (ASSOCIATED(current%mem(iproc)%p)) THEN   ! Process is active on child tile
      !   SELECT TYPE (v => current%mem(iproc)%p%Get_var_real2d(TRIM(var%name)))
      !   TYPE IS (t_jsb_var_real3d)
      !     var_tmp => v
      !   END SELECT
        ! var_tmp => current%mem(iproc)%p%Get_var_real3d(TRIM(var%name))
        ! IF (ASSOCIATED(var_tmp)) THEN
      IF (var%child_idx(i) > 0 .AND. ASSOCIATED(tile%mem(iproc)%p)) THEN
        IF (.NOT. ASSOCIATED(tile%mem(iproc)%p%children(i)%p%vars(var%child_idx(i))%p)) THEN
          CALL finish(TRIM(routine), 'ERROR - this should not happen!')
        ELSE
          var_tmp => tile%mem(iproc)%p%children(i)%p%vars(var%child_idx(i))%p
        END IF
        fractions(:,i) = this%fractions(ics:ice,iblk,i)
        IF (.NOT. var%l_aggregate_all) fract_sum(:)   = fract_sum(:) + fractions(:,i)
        data(:,:,i) = var_tmp%field(ics:ice,:,iblk,1,1)
        ! data(:,:,i) = var_tmp%ptr(ics:ice,:,iblk)
        j = j + 1
      END IF
      ! END IF
      ! current => current%next_sibling_tile
    END DO
    j=no_children
    IF (j == 0) CALL finish(TRIM(routine), 'Variable '//TRIM(var%name)//' not found in any of children of tile '//&
      &                                    TRIM(tile%name)//' - this should not happen!')

    ! IF (ANY(fract_sum(:) > 1._wp + EPSILON(1._wp))) THEN
    !   CALL message(TRIM(routine), 'On tile '//TRIM(tile%name)//': fract_sum > 1 - '// &
    !     &                         real2string(MINVAL(fract_sum(:)))//' - '//real2string(MAXVAL(fract_sum(:))), &
    !     &          all_print=.TRUE.)
    !   CALL finish(TRIM(routine),  '    ... this should not happen')
    ! END IF

    DO ilev=1,SIZE(data,2)
      IF (j == 1 .AND. .NOT. var%l_aggregate_all) THEN
        ! Only one tile to aggregate, avoid multiplying and dividing by a potentially small fraction 
        var%ptr(ics:ice,ilev,iblk) = SUM(data(:,ilev,:), DIM=2)
      ELSE IF (var%l_aggregate_all) THEN              ! Sum of fractions is 1
        ! For tiles where variable is not defined data(:,ilev,tile) is 0
        var%ptr(ics:ice,ilev,iblk) = SUM(fractions(:,:) * data(:,ilev,:), DIM=2)
      ELSE
        WHERE (fract_sum(:) > 0._wp)
          var%ptr(ics:ice,ilev,iblk) = SUM(fractions(:,:) * data(:,ilev,:), DIM=2) / fract_sum(:)
        ! ELSE WHERE
          ! var%ptr(ics:ice,ilev,iblk) = missval
        END WHERE
      END IF
    END DO

    DEALLOCATE(data, fractions)

    IF (var%l_aggregate_all) THEN
      IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//'; all_tiles)')
    ELSE
      IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//')')
    END IF

  END SUBROUTINE Aggregate_weighted_by_fract_3d

#endif
END MODULE mo_jsb_tile_class
