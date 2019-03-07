!>
!! @brief Contains abstract tile class that contains the surface structure and memory state
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_tile_class
#ifndef __NO_JSBACH__

  USE mo_jsb_control,         ONLY: debug_on
  USE mo_kind,                ONLY: wp, dp
  USE mo_util,                ONLY: real2string
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_io_units,            ONLY: filename_max

  USE mo_hsm_class,           ONLY: t_Hsm, t_State, t_Message
  USE mo_jsb_memory_class,    ONLY: t_jsb_memory_p, t_jsb_memory
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct, &
    &                               LAND_TYPE, VEG_TYPE, BARE_TYPE, GLACIER_TYPE, LAKE_TYPE
  USE mo_jsb_var_class,       ONLY: t_jsb_var, t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile_abstract
  PUBLIC :: t_jsb_aggregator, t_jsb_aggregator_p, t_jsb_aggregator_weighted_by_fract

!!$  INTEGER, PARAMETER :: max_depth_of_tile_tree = 3

  ! Declarations for aggregators
  TYPE, ABSTRACT :: t_jsb_aggregator
    CHARACTER(len=25) :: name
  CONTAINS
    PROCEDURE(aggregate_interface_2d), DEFERRED             :: Aggregate_2d
    PROCEDURE(aggregate_interface_3d), DEFERRED             :: Aggregate_3d
    GENERIC                                                 :: Aggregate => Aggregate_2d, Aggregate_3d
    PROCEDURE(aggregate_interface_setopt_logical), DEFERRED :: Set_aggregator_option_logical
    GENERIC                                                 :: Set_option => Set_aggregator_option_logical
  END TYPE t_jsb_aggregator

  TYPE t_jsb_aggregator_p
    CLASS(t_jsb_aggregator), POINTER :: p => NULL()
  END TYPE t_jsb_aggregator_p

  ! Declarations for tile
  TYPE, ABSTRACT, EXTENDS(t_State) :: t_jsb_tile_abstract
    TYPE(t_jsb_lct),        ALLOCATABLE    :: lcts(:)              !< List of land cover types for tile
                                                                   !! The first entry is always the primary LCT for this tile;
                                                                   !! the remaining lcts are those inherited for aggregation
    LOGICAL                                :: is_vegetation = .FALSE.
    LOGICAL                                :: is_bare = .FALSE.
    LOGICAL                                :: is_land = .FALSE.
    LOGICAL                                :: is_glacier = .FALSE.
    LOGICAL                                :: is_lake = .FALSE.
    LOGICAL                                :: contains_vegetation = .FALSE.
    LOGICAL                                :: contains_bare = .FALSE.
    LOGICAL                                :: contains_land = .FALSE.
    LOGICAL                                :: contains_glacier = .FALSE.
    LOGICAL                                :: contains_lake = .FALSE.
    LOGICAL                                :: contains_soil = .FALSE.
    INTEGER,                ALLOCATABLE    :: process_action(:)
    REAL(wp), PRIVATE,          POINTER    :: fract     (:,:)      !< Fraction of this tile (relative to parent tile)
    REAL(wp),                   POINTER    :: delta_fract(:,:)     !< Current operation delta fract of this tile (rel. to parent) 
    REAL(wp),                   POINTER    :: fract_max (:,:)      !< Maximum potential fraction for this tile (rel. to parent)
    REAL(wp), PRIVATE,          POINTER    :: fract_old     (:,:)  !< Old fract
    REAL(wp),                   POINTER    :: fract_max_old (:,:)  !< Old fract_max
    CHARACTER(LEN=filename_max)            :: fract_filename       !< Name of file to read tile fractions from. The variable names
                                                                   !! must be either "fract[_max]_"//tile%name, or specified by
                                                                   !! fract_varname.
    CHARACTER(LEN=20)                      :: fract_varname        !< Name of variable containing fractions for this tile
    LOGICAL,                ALLOCATABLE    :: l_fract_children_changed(:,:)
    LOGICAL,                ALLOCATABLE    :: l_delta_area_change(:,:) !TODO-JN: probably nblks would be sufficient
    LOGICAL,                ALLOCATABLE    :: l_update_c(:,:)          !         required dim only depends on parallelisation
    TYPE(t_jsb_aggregator_p), ALLOCATABLE  :: aggregators(:)
    !CLASS(t_jsb_param_p),       POINTER    :: param(:)
    TYPE(t_jsb_memory_p),       POINTER    :: mem(:)
    ! The next three pointers basically repeat what's already in the base type, but they save time-consuming "SELECT TYPE"s
    CLASS(t_jsb_tile_abstract), POINTER    :: parent_tile => NULL()
    CLASS(t_jsb_tile_abstract), POINTER    :: first_child_tile => NULL()
    CLASS(t_jsb_tile_abstract), POINTER    :: next_sibling_tile => NULL()
    !
    INTEGER                                :: grid_id = -1
    INTEGER                                :: owner_model_id = 0
  CONTAINS
    PROCEDURE                              :: Add_lct
    PROCEDURE                              :: Add_lct_to_parents
    PROCEDURE                              :: Get_parent_memory
    PROCEDURE                              :: Is_process_active
    PROCEDURE                              :: Set_process_action => Set_process_action_tile
    PROCEDURE                              :: Set_process_action_parents
    PROCEDURE                              :: Register_aggregator
    PROCEDURE                              :: Get_aggregator
    PROCEDURE                              :: Set_fraction_1d
    PROCEDURE                              :: Set_fraction_2d
    GENERIC                                :: Set_fraction => Set_fraction_1d, Set_fraction_2d
    !PROCEDURE                              :: Set_fractions_max
    PROCEDURE                              :: Get_fraction_1d
    PROCEDURE                              :: Get_fraction_2d
    GENERIC                                :: Get_fraction => Get_fraction_1d, Get_fraction_2d
    PROCEDURE                              :: Set_fract_ptr
    PROCEDURE                              :: Set_delta_fract_ptr
    PROCEDURE                              :: Set_fract_old_ptr
    PROCEDURE                              :: Set_fract_max_ptr
    PROCEDURE                              :: Set_fract_max_old_ptr
    !PROCEDURE                              :: Get_fractions_max
    PROCEDURE                              :: Get_first_child_tile
    PROCEDURE                              :: Get_next_sibling_tile
    PROCEDURE                              :: Get_prev_sibling_tile
    PROCEDURE(Init_interface),      DEFERRED :: Init
    PROCEDURE(Init_vars_interface), DEFERRED :: Init_vars
    PROCEDURE(Print_interface),     DEFERRED :: Print
    PROCEDURE(Handler_interface),   DEFERRED :: Handler
  END type t_jsb_tile_abstract

  ABSTRACT INTERFACE
    SUBROUTINE Init_interface(this, varlist_name, prefix, suffix, l_fixed_fractions, in_var_groups)
      IMPORT :: t_jsb_tile_abstract
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
      CHARACTER(len=*),           INTENT(in)    :: varlist_name
      CHARACTER(len=*),           INTENT(in)    :: prefix
      CHARACTER(len=*),           INTENT(in)    :: suffix
      LOGICAL,                    INTENT(in)    :: l_fixed_fractions
      LOGICAL,                    INTENT(in)    :: in_var_groups
    END SUBROUTINE Init_interface
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE Init_vars_interface(this)
      IMPORT :: t_jsb_tile_abstract
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    END SUBROUTINE Init_vars_interface
  END INTERFACE

  ABSTRACT INTERFACE
    ! FUNCTION Handler_interface(this, hsm, msg_in) RESULT(return_ptr)
    FUNCTION Handler_interface(this, msg_in) RESULT(return_ptr)
      IMPORT :: t_jsb_tile_abstract, t_Hsm, t_Message
      CLASS(t_jsb_tile_abstract),  INTENT(inout) :: this
      ! CLASS(t_Hsm),                INTENT(inout) :: hsm
      CLASS(t_Message),            INTENT(in)    :: msg_in
      CLASS(t_Message),            POINTER       :: return_ptr
    END FUNCTION Handler_interface
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE Print_interface(this)
      IMPORT :: t_jsb_tile_abstract
      CLASS(t_jsb_tile_abstract), INTENT(in) :: this
    END SUBROUTINE Print_interface
  END INTERFACE

  ABSTRACT INTERFACE
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
  CONTAINS
    PROCEDURE :: Set_fractions => Set_fractions_aggregator
!    PROCEDURE :: Aggregate => Aggregate_weighted_by_fract_2d
!    PROCEDURE :: Aggregate => Aggregate_weighted_by_fract_3d
    PROCEDURE :: Aggregate_2d                  => Aggregate_weighted_by_fract_2d
    PROCEDURE :: Aggregate_3d                  => Aggregate_weighted_by_fract_3d
    PROCEDURE :: Set_aggregator_option_logical => Aggregate_weighted_by_fract_set_option_logical
    FINAL     :: Finalize_aggregator_weighted_by_fract
  END TYPE t_jsb_aggregator_weighted_by_fract

  INTERFACE t_jsb_aggregator_weighted_by_fract
    PROCEDURE Create_aggregator_weighted_by_fract
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_tile_class'

CONTAINS

  FUNCTION Get_parent_memory(this, iproc) RESULT(return_ptr)

    CLASS(t_jsb_tile_abstract),  INTENT(in) :: this
    INTEGER,                     INTENT(in) :: iproc
    CLASS(t_jsb_memory),         POINTER   :: return_ptr

    CLASS(*), POINTER :: parent

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_parent_memory'

    parent => this%Get_parent()
    IF (ASSOCIATED(parent)) THEN
      SELECT TYPE (parent)
      CLASS IS (t_jsb_tile_abstract)
        return_ptr => parent%mem(iproc)%p
      END SELECT
    ELSE
      return_ptr => NULL()
    END IF

  END FUNCTION Get_parent_memory

  FUNCTION Is_process_active(this, iproc) RESULT(return_value)

    CLASS(t_jsb_tile_abstract), INTENT(in) :: this
    INTEGER,                    INTENT(in) :: iproc
    LOGICAL                                :: return_value

    CHARACTER(len=*), PARAMETER :: routine = modname//':Is_process_active'

    return_value = .FALSE.

    SELECT TYPE (this)
    CLASS IS (t_jsb_tile_abstract)
      return_value = ASSOCIATED(this%mem(iproc)%p)
    END SELECT

  END FUNCTION Is_process_active

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

  SUBROUTINE Add_lct(this, lct)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    TYPE(t_jsb_lct),            INTENT(in)    :: lct

    TYPE(t_jsb_lct), ALLOCATABLE :: temp_lcts(:)
    INTEGER :: ilct, nlct

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_lct'

    IF (.NOT. ALLOCATED(this%lcts)) THEN

      ALLOCATE(this%lcts(1))
      this%lcts(1) = lct

      IF (ASSOCIATED(this%parent)) THEN    ! Root tile has no primary LCT
        SELECT CASE (lct%id)
        CASE (VEG_TYPE)
          this%is_vegetation = .TRUE.
        CASE (BARE_TYPE)
          this%is_bare = .TRUE.
        CASE (LAND_TYPE)
          this%is_land = .TRUE.
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

    END IF

    SELECT CASE (lct%id)
    CASE (VEG_TYPE)
      this%contains_vegetation = .TRUE.
    CASE (BARE_TYPE)
      this%contains_bare = .TRUE.
    CASE (LAND_TYPE)
      this%contains_land = .TRUE.
    CASE (GLACIER_TYPE)
      this%contains_glacier = .TRUE.
    CASE (LAKE_TYPE)
      this%contains_lake = .TRUE.
    END SELECT

    this%contains_soil       = this%contains_vegetation .OR. this%contains_bare .OR. this%contains_land

  END SUBROUTINE Add_lct

  SUBROUTINE Add_lct_to_parents(this, lct)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    TYPE(t_jsb_lct),            INTENT(in)    :: lct

    CLASS(t_State), POINTER :: current

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_lct_to_parents'

    current => this%parent

    DO WHILE (ASSOCIATED(current))
      SELECT TYPE (current)
      CLASS IS (t_jsb_tile_abstract)
        CALL current%Add_lct(lct)
      END SELECT
      current => current%parent
    END DO

  END SUBROUTINE Add_lct_to_parents

  SUBROUTINE Set_fraction_1d(this, ics, ice, iblk, fract, fract_max, fract_old, fract_max_old)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    INTEGER,                    INTENT(in)    :: ics, ice, iblk
    REAL(wp), OPTIONAL,         INTENT(in)    :: &
      & fract(:),        &
      & fract_max(:),    &
      & fract_old(:),    &
      & fract_max_old(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_fraction_1d'

    IF (PRESENT(fract)) THEN
      IF (SIZE(fract) /= ice-ics+1) &
       & CALL finish(TRIM(routine), 'Dimension mismatch')

      this%fract_old(ics:ice, iblk) = this%fract(ics:ice, iblk)
      this%fract    (ics:ice, iblk) = fract(:)

      IF (ASSOCIATED(this%parent)) THEN
        SELECT TYPE (parent=>this%parent)
        CLASS IS (t_jsb_tile_abstract)
          WHERE (this%fract(ics:ice, iblk) /= this%fract_old(ics:ice, iblk))
            parent%l_fract_children_changed(ics:ice, iblk) = .TRUE.
          END WHERE
        END SELECT
      END IF
    END IF

    IF (PRESENT(fract_max)) THEN
      IF (SIZE(fract_max) /= ice-ics+1) &
       & CALL finish(TRIM(routine), 'Dimension mismatch')

      this%fract_max_old(ics:ice, iblk) = this%fract_max(ics:ice, iblk)
      this%fract_max    (ics:ice, iblk) = fract_max(:)
    END IF

    IF (PRESENT(fract_old)) THEN
      IF (SIZE(fract_old) /= ice-ics+1) &
       & CALL finish(TRIM(routine), 'Dimension mismatch')

      this%fract_old(ics:ice, iblk) = fract_old(:)
    END IF

    IF (PRESENT(fract_max_old)) THEN
      IF (SIZE(fract_max_old) /= ice-ics+1) &
       & CALL finish(TRIM(routine), 'Dimension mismatch')

      this%fract_max_old(ics:ice, iblk) = fract_max_old(:)
    END IF

  END SUBROUTINE Set_fraction_1d

  SUBROUTINE Set_fraction_2d(this, fract, fract_max, fract_old, fract_max_old)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    REAL(wp), OPTIONAL,         INTENT(in)    :: &
      & fract(:,:),        &
      & fract_max(:,:),    &
      & fract_old(:,:),    &
      & fract_max_old(:,:)

    INTEGER :: ilct

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_fraction_2d'

    IF (PRESENT(fract)) THEN
      IF (SIZE(fract,1) /= SIZE(this%fract,1) .OR. SIZE(fract,2) /= SIZE(this%fract,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')

      IF (ANY(fract(:,:) /= this%fract(:,:))) THEN

        this%fract_old(:,:) = this%fract(:,:)
        this%fract    (:,:) = fract(:,:)

        ! IF (ASSOCIATED(this%parent)) THEN
        !   SELECT TYPE (parent=>this%parent)
        !   CLASS IS (t_jsb_tile_abstract)
        !     WHERE (this%fract(:,:) /= this%fract_old(:,:))
        !       parent%l_fract_children_changed(:,:) = .TRUE.
        !     END WHERE
        !     IF (ANY(parent%l_fract_children_changed(:,:))) THEN
        !       DO ilct=1,SIZE(parent%lcts)
        !         IF (this%lcts(1)%id == parent%lcts(ilct)%id) THEN
        !           parent%lcts(ilct)%fract(:,:) = parent%lcts(ilct)%fract(:,:) - this%fract_old(:,:) + this%fract(:,:)
        !           EXIT
        !         END IF
        !       END DO
        !     END IF
        !   END SELECT
        ! END IF
      END IF
    END IF

    IF (PRESENT(fract_max)) THEN
      IF (SIZE(fract_max,1) /= SIZE(this%fract_max,1) .OR. SIZE(fract_max,2) /= SIZE(this%fract_max,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')

      IF (ANY(fract_max(:,:) /= this%fract_max(:,:))) THEN
        this%fract_max_old(:,:) = this%fract_max(:,:)
        this%fract_max    (:,:) = fract_max(:,:)
      END IF

    END IF

    IF (PRESENT(fract_old)) THEN
      IF (SIZE(fract_old,1) /= SIZE(this%fract_old,1) .OR. SIZE(fract_old,2) /= SIZE(this%fract_old,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')

      IF (ANY(fract_old(:,:) /= this%fract_old(:,:))) THEN
        this%fract_old(:,:) = fract_old(:,:)
      END IF

    END IF

    IF (PRESENT(fract_max_old)) THEN
      IF (SIZE(fract_max_old,1) /= SIZE(this%fract_max_old,1) .OR. SIZE(fract_max_old,2) /= SIZE(this%fract_max_old,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')

      IF (ANY(fract_max_old(:,:) /= this%fract_max_old(:,:))) THEN
        this%fract_max_old(:,:) = fract_max_old(:,:)
      END IF

    END IF

  END SUBROUTINE Set_fraction_2d

  SUBROUTINE Get_fraction_1d(this, ics, ice, iblk, fract, fract_max, fract_old, fract_max_old)

    CLASS(t_jsb_tile_abstract), INTENT(in)    :: this
    INTEGER,                    INTENT(in)    :: ics, ice, iblk
    REAL(wp), OPTIONAL,         INTENT(inout) :: &
      & fract(:),        &
      & fract_max(:),    &
      & fract_old(:),    &
      & fract_max_old(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_fraction_1d'

    IF (PRESENT(fract)) THEN
      IF (SIZE(fract) /= ice-ics+1) CALL finish(TRIM(routine), 'Dimension mismatch')
      fract(:) = this%fract(ics:ice, iblk)
    END IF

    IF (PRESENT(fract_max)) THEN
      IF (SIZE(fract_max) /= ice-ics+1) CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_max(:) = this%fract_max(ics:ice, iblk)
    END IF

    IF (PRESENT(fract_old)) THEN
      IF (SIZE(fract_old) /= ice-ics+1) CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_old(:) = this%fract_old(ics:ice, iblk)
    END IF

    IF (PRESENT(fract_max_old)) THEN
      IF (SIZE(fract_max_old) /= ice-ics+1) CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_max_old(:) = this%fract_max_old(ics:ice, iblk)
    END IF

  END SUBROUTINE Get_fraction_1d

  SUBROUTINE Get_fraction_2d(this, fract, fract_max, fract_old, fract_max_old)

    CLASS(t_jsb_tile_abstract), INTENT(in)    :: this
    REAL(wp), OPTIONAL,         INTENT(inout) :: &
      & fract(:,:),        &
      & fract_max(:,:),    &
      & fract_old(:,:),    &
      & fract_max_old(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_fraction_2d'

    IF (PRESENT(fract)) THEN
      IF (SIZE(fract,1) /= SIZE(this%fract,1) .OR. SIZE(fract,2) /= SIZE(this%fract,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      fract(:,:) = this%fract(:,:)
    END IF

    IF (PRESENT(fract_max)) THEN
      IF (SIZE(fract_max,1) /= SIZE(this%fract_max,1) .OR. SIZE(fract_max,2) /= SIZE(this%fract_max,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_max(:,:) = this%fract_max(:,:)
    END IF

    IF (PRESENT(fract_old)) THEN
      IF (SIZE(fract_old,1) /= SIZE(this%fract_old,1) .OR. SIZE(fract_old,2) /= SIZE(this%fract_old,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_old(:,:) = this%fract_old(:,:)
    END IF

    IF (PRESENT(fract_max_old)) THEN
      IF (SIZE(fract_max_old,1) /= SIZE(this%fract_max_old,1) .OR. SIZE(fract_max_old,2) /= SIZE(this%fract_max_old,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_max_old(:,:) = this%fract_max_old(:,:)
    END IF

  END SUBROUTINE Get_fraction_2d

  SUBROUTINE Set_fract_ptr(this, ptr)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    REAL(wp), POINTER                         :: ptr(:,:)

    this%fract => ptr

  END SUBROUTINE Set_fract_ptr

  SUBROUTINE Set_delta_fract_ptr(this, ptr)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    REAL(wp), POINTER                         :: ptr(:,:)

    this%delta_fract => ptr

  END SUBROUTINE Set_delta_fract_ptr

  SUBROUTINE Set_fract_old_ptr(this, ptr)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    REAL(wp), POINTER                         :: ptr(:,:)

    this%fract_old => ptr

  END SUBROUTINE Set_fract_old_ptr

  SUBROUTINE Set_fract_max_ptr(this, ptr)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    REAL(wp), POINTER                         :: ptr(:,:)

    this%fract_max => ptr

  END SUBROUTINE Set_fract_max_ptr

  SUBROUTINE Set_fract_max_old_ptr(this, ptr)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    REAL(wp), POINTER                         :: ptr(:,:)

    this%fract_max_old => ptr

  END SUBROUTINE Set_fract_max_old_ptr

  SUBROUTINE Register_aggregator(this, aggregator)

    CLASS(t_jsb_tile_abstract), INTENT(inout)        :: this
    CLASS(t_jsb_aggregator),    INTENT(in),   TARGET :: aggregator

    TYPE(t_jsb_aggregator_p), ALLOCATABLE :: temp_aggregators(:)
    INTEGER :: n
    CHARACTER(len=*), PARAMETER :: routine = modname//':Register_aggregator'

    IF (ALLOCATED(this%aggregators)) THEN
      n = SIZE(this%aggregators)
    ELSE
      n = 0
    END IF
    ALLOCATE(temp_aggregators(n+1))
    IF (ALLOCATED(this%aggregators)) temp_aggregators(1:n) = this%aggregators
    temp_aggregators(n+1) = t_jsb_aggregator_p(aggregator)
    CALL move_ALLOC(temp_aggregators, this%aggregators)

    IF (debug_on()) CALL message(TRIM(routine), 'Registered aggregator '//TRIM(aggregator%name)//' for tile '//TRIM(this%name))

  END SUBROUTINE Register_aggregator

  FUNCTION Get_aggregator(this, name, copy) RESULT(aggregator)

    CLASS(t_jsb_tile_abstract), INTENT(in) :: this
    CHARACTER(len=*),           INTENT(in) :: name
    LOGICAL, OPTIONAL,          INTENT(in) :: copy
    CLASS(t_jsb_aggregator),    POINTER    :: aggregator

    INTEGER :: i
    LOGICAL :: l_copy

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_aggregator'

    l_copy = .FALSE.
    IF (PRESENT(copy)) l_copy = copy

    aggregator => NULL()
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

    IF (.NOT.(ASSOCIATED(aggregator))) &
      & CALL finish(TRIM(routine), 'Aggregator "'//TRIM(name)//'" not found on tile '//TRIM(this%name))

  END FUNCTION Get_aggregator

  FUNCTION Create_aggregator_weighted_by_fract(grid_id, no_of_children) RESULT(aggregator)

    USE mo_jsb_grid_class, ONLY: t_jsb_grid
    USE mo_jsb_grid,       ONLY: Get_grid

    INTEGER, INTENT(in) :: grid_id
    INTEGER, INTENT(in) :: no_of_children
    TYPE(t_jsb_aggregator_weighted_by_fract), POINTER :: aggregator

    TYPE(t_jsb_grid), POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':Create_aggregator_weighted_by_fract'

    IF (no_of_children < 1) CALL finish(TRIM(routine), 'This should not happen!')

    ALLOCATE(aggregator)
    aggregator%name = 'weighted_by_fract'

    grid => Get_grid(grid_id)
    ALLOCATE(aggregator%fractions(grid%nproma, grid%nblks, no_of_children))
    ALLOCATE(aggregator%data_sum (grid%nproma, grid%nblks))
    ALLOCATE(aggregator%fract_sum(grid%nproma, grid%nblks))

  END FUNCTION Create_aggregator_weighted_by_fract

  SUBROUTINE Finalize_aggregator_weighted_by_fract(aggregator)

    TYPE(t_jsb_aggregator_weighted_by_fract) :: aggregator

    CHARACTER(len=*), PARAMETER :: routine = modname//':Finalize_aggregator_weighted_by_fract'

    !CALL message(TRIM(routine), 'Finalizing aggregator '//TRIM(aggregator%name))

    IF (ALLOCATED(aggregator%fractions)) DEALLOCATE(aggregator%fractions)
    IF (ALLOCATED(aggregator%data_sum))  DEALLOCATE(aggregator%data_sum)
    IF (ALLOCATED(aggregator%fract_sum)) DEALLOCATE(aggregator%fract_sum)

  END SUBROUTINE Finalize_aggregator_weighted_by_fract

  SUBROUTINE Set_fractions_aggregator(this, tile, ics, ice, iblk)

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile
    INTEGER,                                   INTENT(in)    :: ics, ice, iblk

    INTEGER :: no_children, i
    CLASS(t_State), POINTER :: current
    REAL(wp), ALLOCATABLE :: fract(:), fract_sum(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_fractions_aggregator'

    ! Check if child tile fractions have changed and, if necessary, update fractions of aggregator
    ! by looping over children of tile
    IF (ANY(tile%l_fract_children_changed(ics:ice,iblk))) THEN
      no_children = tile%Get_no_of_children()
      current => tile%first_child
      DO i=1,no_children
        SELECT TYPE (current)
        CLASS IS (t_jsb_tile_abstract)
          CALL current%Get_fraction(ics, ice, iblk, fract=this%fractions(ics:ice,iblk,i))
        END SELECT
        current => current%next_sibling
      END DO

      ALLOCATE(fract(ice-ics+1))
      ALLOCATE(fract_sum(ice-ics+1))
      fract(:) = tile%fract(ics:ice,iblk)
      fract_sum(:) = SUM(this%fractions(ics:ice,iblk,:), DIM=2)
      ! IF (ASSOCIATED(tile%parent) .AND. ANY( fract_sum(:) < 1._wp - EPSILON(1._wp) .AND. fract(:) > 0._wp)) THEN
      !   CALL message(TRIM(routine), 'Sum of child tile fractions on tile '//TRIM(tile%name)//': '//&
      !     &                         real2string(MINVAL(fract_sum(:)))//' - '//real2string(MAXVAL(fract_sum(:))), &
      !     &          all_print=.TRUE.)
      !   CALL finish(TRIM(routine), 'Sum of tile fractions is not 1 on tile')
      ! END IF
      DEALLOCATE(fract, fract_sum)

      tile%l_fract_children_changed(ics:ice,iblk) = .FALSE.
    END IF

  END SUBROUTINE Set_fractions_aggregator

  SUBROUTINE Aggregate_weighted_by_fract_2d(this, tile, var, ics, ice, iblk, caller)

    USE mo_hsm_class,      ONLY: t_State

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile
    TYPE(t_jsb_var_real2d),                    INTENT(inout) :: var
    INTEGER,                                   INTENT(in)    :: ics, ice, iblk
    CHARACTER(len=*),                          INTENT(in)    :: caller

    INTEGER :: no_children, i, iproc, nc
    REAL(dp) :: missval

    CHARACTER(len=256) :: routine

    IF (debug_on()) THEN
      routine = modname//':Aggregate_weighted_by_fract_2d'//'(called by '//TRIM(caller)//')'
    ELSE
      routine = modname//':Aggregate_weighted_by_fract_2d'
    END IF

    IF (.NOT. ASSOCIATED(var%ptr)) RETURN

    CALL this%Set_fractions(tile, ics, ice, iblk)

    nc = ice - ics + 1

    missval = var%missval

    iproc = var%owner_proc_id
    IF (iproc < 0) CALL finish(TRIM(routine), 'Unknown process for variable '//TRIM(var%name)//' on tile '//TRIM(tile%name))

    IF (.NOT. ASSOCIATED(tile%mem(iproc)%p)) RETURN

    no_children = tile%Get_no_of_children()

    IF (.NOT. ALLOCATED(var%child_idx)) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' not allocated')
    ELSE IF (SIZE(var%child_idx) /= no_children) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' has wrong size')
    END IF

    this%data_sum (:,iblk) = 0._wp
    this%fract_sum(:,iblk) = 0._wp

    DO i=1,no_children
      IF (var%child_idx(i) > 0) THEN
        ASSOCIATE (data_child => tile%mem(iproc)%p%children(i)%p%vars(var%child_idx(i))%p%field(ics:ice,iblk,1,1,1), &
          &        fractions  => this%fractions(ics:ice,iblk,i))
          this%data_sum(ics:ice,iblk) = this%data_sum(ics:ice,iblk) + fractions * data_child
          this%fract_sum(ics:ice,iblk) = this%fract_sum(ics:ice,iblk) + fractions
        END ASSOCIATE
      END IF
    END DO
    ! print*,'AAA ',ics,ice,MINVAL(this%data_sum(ics:ice)), MAXVAL(this%data_sum(ics:ice)),&
      ! & MINVAL(this%fract_sum(ics:ice)), MAXVAL(this%fract_sum(ics:ice))

    ! IF (ANY(fract_sum(:) > 1._wp + EPSILON(1._wp))) THEN
    !   ! print*, fract_sum
    !   CALL message(TRIM(routine), 'On tile '//TRIM(tile%name)//': fract_sum > 1 - '// &
    !     &                         real2string(MINVAL(fract_sum(:)))//' - '//real2string(MAXVAL(fract_sum(:))), &
    !     &          all_print=.TRUE.)
    !   CALL finish(TRIM(routine),  '    ... this should not happen')
    ! END IF

    IF (no_children == 1 .AND. .NOT. var%l_aggregate_all) THEN
      ! Only one tile to aggregate, avoid multiplying and dividing by a potentially small fraction
      var%ptr(ics:ice,iblk) = this%data_sum(ics:ice,iblk)
    ELSE IF (var%l_aggregate_all) THEN              ! Sum of fractions is 1
      ! For tiles where variable is not defined data(:,tile) is 0
      var%ptr(ics:ice,iblk) = this%data_sum(ics:ice,iblk)
    ELSE
      WHERE (this%fract_sum(ics:ice,iblk) > 0._wp)
        var%ptr(ics:ice,iblk) = this%data_sum(ics:ice,iblk) / this%fract_sum(ics:ice,iblk)
      !ELSE WHERE
      !  var%ptr(ics:ice,iblk) = missval
      END WHERE
    END IF

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

    INTEGER :: no_children, i, iproc, nc, ilev
    REAL(dp) :: missval

    CHARACTER(len=256) :: routine

    IF (debug_on()) THEN
      routine = modname//':Aggregate_weighted_by_fract_3d'//'(called by '//TRIM(caller)//')'
    ELSE
      routine = modname//':Aggregate_weighted_by_fract_3d'
    END IF

    IF (.NOT. ASSOCIATED(var%ptr)) RETURN

    CALL this%Set_fractions(tile, ics, ice, iblk)

    nc = ice - ics + 1

    missval = var%missval

    iproc = var%owner_proc_id
    IF (iproc < 0) CALL finish(TRIM(routine), 'Unknown process for variable '//TRIM(var%name)//' on tile '//TRIM(tile%name))

    IF (.NOT. ASSOCIATED(tile%mem(iproc)%p)) RETURN

    no_children = tile%Get_no_of_children()

    IF (.NOT. ALLOCATED(var%child_idx)) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' not allocated')
    ELSE IF (SIZE(var%child_idx) /= no_children) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' has wrong size')
    END IF

    DO ilev=1,SIZE(var%ptr,2)
      this%fract_sum(:,iblk) = 0._wp
      this%data_sum (:,iblk) = 0._wp
      DO i=1,no_children
        IF (var%child_idx(i) > 0) THEN
          ASSOCIATE (data_child => tile%mem(iproc)%p%children(i)%p%vars(var%child_idx(i))%p%field(ics:ice,ilev,iblk,1,1), &
            &        fractions  => this%fractions(ics:ice,iblk,i))
            this%data_sum(ics:ice,iblk) = this%data_sum(ics:ice,iblk) + fractions * data_child
            this%fract_sum(ics:ice,iblk) = this%fract_sum(ics:ice,iblk) + fractions
          END ASSOCIATE
        END IF
      END DO
      IF (no_children == 1 .AND. .NOT. var%l_aggregate_all) THEN
        ! Only one tile to aggregate, avoid multiplying and dividing by a potentially small fraction
        var%ptr(ics:ice,ilev,iblk) = this%data_sum(ics:ice,iblk)
      ELSE IF (var%l_aggregate_all) THEN              ! Sum of fractions is 1
        ! For tiles where variable is not defined data(:,ilev,tile) is 0
        var%ptr(ics:ice,ilev,iblk) = this%data_sum(ics:ice,iblk)
      ELSE
        WHERE (this%fract_sum(ics:ice,iblk) > 0._wp)
          var%ptr(ics:ice,ilev,iblk) = this%data_sum(ics:ice,iblk) / this%fract_sum(ics:ice,iblk)
        ! ELSE WHERE
          ! var%ptr(ics:ice,ilev,iblk) = missval
        END WHERE
      END IF
    END DO

    ! IF (ANY(fract_sum(:) > 1._wp + EPSILON(1._wp))) THEN
    !   CALL message(TRIM(routine), 'On tile '//TRIM(tile%name)//': fract_sum > 1 - '// &
    !     &                         real2string(MINVAL(fract_sum(:)))//' - '//real2string(MAXVAL(fract_sum(:))), &
    !     &          all_print=.TRUE.)
    !   CALL finish(TRIM(routine),  '    ... this should not happen')
    ! END IF

    IF (var%l_aggregate_all) THEN
      IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//'; all_tiles)')
    ELSE
      IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//')')
    END IF

  END SUBROUTINE Aggregate_weighted_by_fract_3d

  SUBROUTINE Aggregate_weighted_by_fract_set_option_logical(this, name, value)

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this
    CHARACTER(len=*),                          INTENT(in)    :: name
    LOGICAL,                                   INTENT(in)    :: value

    ! SELECT CASE (TRIM(name))
    ! CASE ('option')
    !   this%option = value
    ! END SELECT

  END SUBROUTINE Aggregate_weighted_by_fract_set_option_logical

  FUNCTION Get_first_child_tile(this) RESULT(return_ptr)

    CLASS(t_jsb_tile_abstract), INTENT(in)  :: this
    CLASS(t_jsb_tile_abstract), POINTER     :: return_ptr

    CLASS(t_State), POINTER :: next

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_tile_class:Get_first_child_tile'

    next => this%first_child
    IF (.NOT. ASSOCIATED(next)) THEN
      return_ptr => NULL()
      RETURN
    END IF

    SELECT TYPE (next)
    CLASS IS (t_jsb_tile_abstract)
      return_ptr => next
    CLASS IS (t_State)
      CALL finish(TRIM(routine), 'current tile is of type t_State, should be t_jsb_tile_abstract')
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

    NULLIFY(next)

    IF (.NOT. ASSOCIATED(return_ptr)) &
      & CALL finish(TRIM(routine), 'Could not find first child tile')

  END FUNCTION Get_first_child_tile

  FUNCTION Get_next_sibling_tile(this) RESULT(return_ptr)

    CLASS(t_jsb_tile_abstract), INTENT(in)  :: this
    CLASS(t_jsb_tile_abstract), POINTER :: return_ptr

    CLASS(t_State), POINTER :: next

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_tile_class:Get_next_sibling_tile'

    next => this%next_sibling
    IF (.NOT. ASSOCIATED(next)) THEN
      return_ptr => NULL()
      RETURN
    END IF

    SELECT TYPE (next)
    CLASS IS (t_jsb_tile_abstract)
      return_ptr => next
    CLASS IS (t_State)
      CALL finish(TRIM(routine), 'current tile is of type t_State, should be t_jsb_tile_abstract')
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

    NULLIFY(next)

    IF (.NOT. ASSOCIATED(return_ptr)) &
      & CALL finish(TRIM(routine), 'Could not find next sibling tile')

  END FUNCTION Get_next_sibling_tile

  FUNCTION Get_prev_sibling_tile(this) RESULT(return_ptr)

    CLASS(t_jsb_tile_abstract), INTENT(in)  :: this
    CLASS(t_jsb_tile_abstract), POINTER :: return_ptr

    CLASS(t_State), POINTER :: prev

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_tile_class:Get_prev_sibling_tile'

    prev => this%prev_sibling
    IF (.NOT. ASSOCIATED(prev)) THEN
      return_ptr => NULL()
      RETURN
    END IF

    SELECT TYPE (prev)
    CLASS IS (t_jsb_tile_abstract)
      return_ptr => prev
    CLASS IS (t_State)
      CALL finish(TRIM(routine), 'current tile is of type t_State, should be t_jsb_tile_abstract')
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

    NULLIFY(prev)

    IF (.NOT. ASSOCIATED(return_ptr)) &
      & CALL finish(TRIM(routine), 'Could not find prev sibling tile')

  END FUNCTION Get_prev_sibling_tile

#endif
END MODULE mo_jsb_tile_class
