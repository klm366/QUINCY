!>
!! @brief Contains abstract tile class that contains the surface structure and memory state
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_tile_class
#ifndef __NO_JSBACH__

  USE mo_jsb_control,         ONLY: debug_on
  USE mo_kind,                ONLY: wp, dp
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_io_units,            ONLY: filename_max

  USE mo_hsm_class,           ONLY: t_Hsm, t_State, t_Message
  USE mo_jsb_memory_class,    ONLY: t_jsb_memory_p, t_jsb_memory
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct
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
    PROCEDURE(aggregate_interface_2d), DEFERRED :: Aggregate_2d
    PROCEDURE(aggregate_interface_3d), DEFERRED :: Aggregate_3d
    GENERIC                                     :: Aggregate => Aggregate_2d, Aggregate_3d
  END TYPE t_jsb_aggregator

  TYPE t_jsb_aggregator_p
    CLASS(t_jsb_aggregator), POINTER :: p => NULL()
  END TYPE t_jsb_aggregator_p

  ! Declarations for tile
  TYPE, ABSTRACT, EXTENDS(t_State) :: t_jsb_tile_abstract
    TYPE(t_jsb_lct),        ALLOCATABLE    :: lcts(:)              !< List of land cover types for tile
                                                                   !! The first entry is always the primary LCT for this tile;
                                                                   !! the remaining lcts are those inherited for aggregation
    INTEGER,                ALLOCATABLE    :: process_action(:)
    REAL(wp), PRIVATE,          POINTER    :: fract     (:,:)      !< Fraction of this tile (relative to parent tile)
    REAL(wp),                   POINTER    :: fract_max (:,:)      !< Maximum potential fraction for this tile (rel. to parent)
    REAL(wp), PRIVATE,          POINTER    :: fract_old     (:,:)  !< Old fract
    REAL(wp),                   POINTER    :: fract_max_old (:,:)  !< Old fract_max
    CHARACTER(LEN=filename_max)            :: fract_filename       !< Name of file to read tile fractions from (all tiles except
                                                                   !! root). The variable names must be "fract[_max]_"//tile%name
    LOGICAL,                ALLOCATABLE    :: l_fract_children_changed(:,:)
    TYPE(t_jsb_aggregator_p), ALLOCATABLE  :: aggregators(:)
    !CLASS(t_jsb_param_p),       POINTER    :: param(:)
    CLASS(t_jsb_memory_p),      POINTER    :: mem(:)
    CLASS(t_jsb_tile_abstract), POINTER    :: parent_tile => NULL()
    INTEGER                                :: grid_id = -1
    INTEGER                                :: owner_model_id = 0
  CONTAINS
    PROCEDURE                              :: Add_lct
    PROCEDURE                              :: Add_lct_to_parents
    PROCEDURE                              :: Get_parent_memory
    PROCEDURE                              :: Is_process_active
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
    PROCEDURE                              :: Set_fract_old_ptr
    PROCEDURE                              :: Set_fract_max_ptr
    PROCEDURE                              :: Set_fract_max_old_ptr
    !PROCEDURE                              :: Get_fractions_max
    PROCEDURE(Init_interface),    DEFERRED :: Init
    PROCEDURE(Print_interface),   DEFERRED :: Print
    PROCEDURE(Handler_interface), DEFERRED :: Handler
  END type t_jsb_tile_abstract

  ABSTRACT INTERFACE
    SUBROUTINE Init_interface(this, varlist_name, prefix, suffix)
      IMPORT :: t_jsb_tile_abstract
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
      CHARACTER(len=*),           INTENT(in)    :: varlist_name
      CHARACTER(len=*),           INTENT(in)    :: prefix
      CHARACTER(len=*),           INTENT(in)    :: suffix
    END SUBROUTINE Init_interface
  END INTERFACE

  ABSTRACT INTERFACE
    FUNCTION Handler_interface(this, hsm, msg_in) RESULT(return_ptr)
      IMPORT :: t_jsb_tile_abstract, t_Hsm, t_Message
      CLASS(t_jsb_tile_abstract),  INTENT(inout) :: this
      CLASS(t_Hsm),                INTENT(in)    :: hsm
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
    SUBROUTINE aggregate_interface_2d(this, tile, var, ics, ice, iblk)
      IMPORT :: t_jsb_aggregator, t_jsb_tile_abstract, t_jsb_var_real2d
      CLASS(t_jsb_aggregator),    INTENT(inout) :: this
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
      TYPE(t_jsb_var_real2d),     INTENT(inout) :: var
      INTEGER,                    INTENT(in)    :: ics, ice, iblk
    END SUBROUTINE aggregate_interface_2d
    SUBROUTINE aggregate_interface_3d(this, tile, var, ics, ice, iblk)
      IMPORT :: t_jsb_aggregator, t_jsb_tile_abstract, t_jsb_var_real3d
      CLASS(t_jsb_aggregator),    INTENT(inout) :: this
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
      TYPE(t_jsb_var_real3d),     INTENT(inout) :: var
      INTEGER,                    INTENT(in)    :: ics, ice, iblk
    END SUBROUTINE aggregate_interface_3d
  END INTERFACE

  TYPE, EXTENDS(t_jsb_aggregator) :: t_jsb_aggregator_weighted_by_fract
    REAL(wp), ALLOCATABLE :: fractions(:,:,:)
  CONTAINS
    PROCEDURE :: Set_fractions => Set_fractions_aggregator
!    PROCEDURE :: Aggregate => Aggregate_weighted_by_fract_2d
!    PROCEDURE :: Aggregate => Aggregate_weighted_by_fract_3d
    PROCEDURE :: Aggregate_2d => Aggregate_weighted_by_fract_2d
    PROCEDURE :: Aggregate_3d => Aggregate_weighted_by_fract_3d
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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_parent_memory'

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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Is_process_active'

    return_value = .FALSE.

    SELECT TYPE (this)
    CLASS IS (t_jsb_tile_abstract)
      return_value = ASSOCIATED(this%mem(iproc)%p)
    END SELECT

  END FUNCTION Is_process_active

  SUBROUTINE Set_process_action_parents(this, iproc, action)

    CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: this
    INTEGER,                    INTENT(in)            :: iproc
    INTEGER,                    INTENT(in)            :: action

    CLASS(t_jsb_tile_abstract), POINTER :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Set_process_action_parents'

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

    INTEGER :: ilct

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_lct'

    IF (.NOT. ALLOCATED(this%lcts)) THEN
      ALLOCATE(this%lcts(1))
      this%lcts(1) = lct
    ELSE
      DO ilct=1,SIZE(this%lcts)
        ! If this lct is already registered in the tile%lcts, do nothing and return
        IF (this%lcts(ilct)%id == lct%id) RETURN
      END DO
      this%lcts = (/ this%lcts, lct /)
    END IF

  END SUBROUTINE Add_lct

  SUBROUTINE Add_lct_to_parents(this, lct)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    TYPE(t_jsb_lct),            INTENT(in)    :: lct

    CLASS(t_State), POINTER :: current

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_lct_to_parents'

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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Set_fraction_1d'

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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Set_fraction_2d'

    IF (PRESENT(fract)) THEN
      IF (SIZE(fract,1) /= SIZE(this%fract,1) .OR. SIZE(fract,2) /= SIZE(this%fract,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')

      this%fract_old(:,:) = this%fract(:,:)
      this%fract    (:,:) = fract(:,:)

      IF (ASSOCIATED(this%parent)) THEN
        SELECT TYPE (parent=>this%parent)
        CLASS IS (t_jsb_tile_abstract)
          WHERE (this%fract(:,:) /= this%fract_old(:,:))
            parent%l_fract_children_changed(:,:) = .TRUE.
          END WHERE
        END SELECT
      END IF
    END IF

    IF (PRESENT(fract_max)) THEN
      IF (SIZE(fract_max,1) /= SIZE(this%fract_max,1) .OR. SIZE(fract_max,2) /= SIZE(this%fract_max,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')

      this%fract_max_old(:,:) = this%fract_max(:,:)
      this%fract_max    (:,:) = fract_max(:,:)
    END IF

    IF (PRESENT(fract_old)) THEN
      IF (SIZE(fract_old,1) /= SIZE(this%fract_old,1) .OR. SIZE(fract_old,2) /= SIZE(this%fract_old,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')

      this%fract_old(:,:) = fract_old(:,:)
    END IF

    IF (PRESENT(fract_max_old)) THEN
      IF (SIZE(fract_max_old,1) /= SIZE(this%fract_max_old,1) .OR. SIZE(fract_max_old,2) /= SIZE(this%fract_max_old,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')

      this%fract_max_old(:,:) = fract_max_old(:,:)
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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_fraction_1d'

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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_fraction_2d'

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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Register_aggregator'


    IF (.NOT. ALLOCATED(this%aggregators)) THEN
      this%aggregators = (/ t_jsb_aggregator_p(aggregator) /)
    ELSE
      this%aggregators = (/ this%aggregators, t_jsb_aggregator_p(aggregator) /)
    END IF

    CALL message(TRIM(routine), 'Registered aggregator '//TRIM(aggregator%name)//' for tile '//TRIM(this%name))

  END SUBROUTINE Register_aggregator

  FUNCTION Get_aggregator(this, name) RESULT(aggregator)

    USE mo_hsm_class,      ONLY: t_State

    CLASS(t_jsb_tile_abstract), INTENT(in) :: this
    CHARACTER(len=*),           INTENT(in) :: name
    CLASS(t_jsb_aggregator),    POINTER    :: aggregator

    CLASS(t_State),   POINTER :: current
    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_aggregator'

    aggregator => NULL()
    DO i=1,SIZE(this%aggregators)
      IF (TRIM(this%aggregators(i)%p%name) == TRIM(name)) THEN
        aggregator => this%aggregators(i)%p
        EXIT
      END IF
    END DO

    IF (.NOT.(ASSOCIATED(aggregator))) &
      & CALL finish(TRIM(routine), 'Aggregator "'//TRIM(name)//'" not found')

  END FUNCTION Get_aggregator

  FUNCTION Create_aggregator_weighted_by_fract(tile) RESULT(aggregator)

    USE mo_jsb_grid_class, ONLY: t_jsb_grid
    USE mo_jsb_grid,       ONLY: Get_grid

    CLASS(t_jsb_tile_abstract), INTENT(inout)          :: tile
    CLASS(t_jsb_aggregator_weighted_by_fract), POINTER :: aggregator

    INTEGER                   :: no_children
    TYPE(t_jsb_grid), POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Create_aggregator_weighted_by_fract'

    IF (.NOT. tile%Has_children()) CALL finish(TRIM(routine), 'This should not happen!')

    ALLOCATE(t_jsb_aggregator_weighted_by_fract::aggregator)
    aggregator%name = 'weighted_by_fract'

    grid => Get_grid(tile%grid_id)
    no_children = tile%Get_no_of_children()
    ALLOCATE(aggregator%fractions(grid%nproma, grid%nblks, no_children))

  END FUNCTION Create_aggregator_weighted_by_fract

  SUBROUTINE Set_fractions_aggregator(this, tile, ics, ice, iblk)

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile
    INTEGER,                                   INTENT(in)    :: ics, ice, iblk

    INTEGER :: no_children, i
    CLASS(t_State), POINTER :: current

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Set_fractions_aggregator'

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

      IF (ASSOCIATED(tile%parent) .AND. ANY(SUM(this%fractions(ics:ice,iblk,:), DIM=2) < 1._wp - EPSILON(1._wp))) &
        & CALL finish(TRIM(routine), 'Sum of tile fractions is not 1')

      tile%l_fract_children_changed(ics:ice,iblk) = .FALSE.
    END IF

  END SUBROUTINE Set_fractions_aggregator

  SUBROUTINE Aggregate_weighted_by_fract_2d(this, tile, var, ics, ice, iblk)

    USE mo_hsm_class,      ONLY: t_State

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile
    TYPE(t_jsb_var_real2d),                    INTENT(inout) :: var
    INTEGER,                                   INTENT(in)    :: ics, ice, iblk

    INTEGER :: no_children, i, j, iproc, nc
    CLASS(t_State), POINTER :: current
    TYPE(t_jsb_var_real2d), POINTER :: var_tmp
    REAL(wp), ALLOCATABLE :: data(:,:), fractions(:,:)    ! Collected data and fractions from children
    REAL(wp), ALLOCATABLE :: data_sum(:), fract_sum(:)
    REAL(dp) :: missval

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Aggregate_weighted_by_fract_2d'

    IF (.NOT. ASSOCIATED(var%ptr)) RETURN

    CALL this%Set_fractions(tile, ics, ice, iblk)

    nc = ice - ics + 1

    missval = var%missval

    iproc = var%owner_proc_id
    IF (iproc < 0) CALL finish(TRIM(routine), 'Unknown process for variable '//TRIM(var%name)//' on tile '//TRIM(tile%name))

    no_children = tile%Get_no_of_children()

    ALLOCATE(data     (nc,no_children))
    ALLOCATE(fractions(nc,no_children))
    ALLOCATE(data_sum (nc))
    ALLOCATE(fract_sum(nc))

    j = 0
    fractions(:,:) = 0._wp
    fract_sum(:) = 0._wp
    current => tile%first_child
    DO i=1,no_children
      SELECT TYPE (current)
      CLASS IS (t_jsb_tile_abstract)
        var_tmp => current%mem(iproc)%p%Get_var_real2d(TRIM(var%name))
        IF (ASSOCIATED(var_tmp)) THEN
          fractions(:,i) = this%fractions(ics:ice,iblk,i)
          fract_sum(:)   = fract_sum(:) + fractions(:,i)
          data     (:,i) = var_tmp%ptr(ics:ice,iblk)
          j = j + 1
        ELSE
          data(:,i) = 0._wp
        END IF
      END SELECT
      current => current%next_sibling
    END DO

    IF (j == 0) CALL finish(TRIM(routine), 'Variable '//TRIM(var%name)//' not found in any of children of tile '//&
      &                                    TRIM(tile%name)//' - this should not happen!')

    IF (ANY(fract_sum(:) > 1._wp)) CALL finish(TRIM(routine), 'fract_sum > 1 - this should not happen')

    data_sum(:) = SUM(fractions(:,:) * data(:,:), DIM=2)
    WHERE (fract_sum(:) > EPSILON(1._wp))
      var%ptr(ics:ice,iblk) = data_sum(:) / fract_sum(:)
    ELSE WHERE
      ! var%ptr(ics:ice,iblk) = SPREAD(var%missval, DIM=1, NCOPIES=nc) 
      var%ptr(ics:ice,iblk) = missval
    END WHERE

    IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//')')

  END SUBROUTINE Aggregate_weighted_by_fract_2d

  SUBROUTINE Aggregate_weighted_by_fract_3d(this, tile, var, ics, ice, iblk)

    USE mo_hsm_class,      ONLY: t_State

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile
    TYPE(t_jsb_var_real3d),                    INTENT(inout) :: var
    INTEGER,                                   INTENT(in)    :: ics, ice, iblk

    INTEGER :: no_children, i, j, iproc, nc, ilev
    CLASS(t_State), POINTER :: current
    TYPE(t_jsb_var_real3d), POINTER :: var_tmp
    REAL(wp), ALLOCATABLE :: data(:,:,:), fractions(:,:)    ! Collected data and fractions from children
    REAL(wp), ALLOCATABLE :: data_sum(:), fract_sum(:)
    REAL(dp) :: missval

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Aggregate_weighted_by_fract_3d'

    IF (.NOT. ASSOCIATED(var%ptr)) RETURN

    CALL this%Set_fractions(tile, ics, ice, iblk)

    nc = ice - ics + 1

    missval = var%missval

    iproc = var%owner_proc_id
    IF (iproc < 0) CALL finish(TRIM(routine), 'Unknown process for variable '//TRIM(var%name)//' on tile '//TRIM(tile%name))

    no_children = tile%Get_no_of_children()

    ALLOCATE(data     (nc, SIZE(var%ptr,2), no_children))
    ALLOCATE(fractions(nc,                  no_children))
    ALLOCATE(data_sum (nc                              ))
    ALLOCATE(fract_sum(nc                              ))

    j = 0
    fractions(:,:) = 0._wp
    fract_sum(:) = 0._wp
    current => tile%first_child
    DO i=1,no_children
      SELECT TYPE (current)
      CLASS IS (t_jsb_tile_abstract)
        var_tmp => current%mem(iproc)%p%Get_var_real3d(TRIM(var%name))
        IF (ASSOCIATED(var_tmp)) THEN
          fractions(:,i) = this%fractions(ics:ice,iblk,i)
          fract_sum(:)   = fract_sum(:) + fractions(:,i)
          data     (:,:,i) = var_tmp%ptr(ics:ice,:,iblk)
          j = j + 1
        ELSE
          data(:,:,i) = 0._wp
        END IF
      END SELECT
      current => current%next_sibling
    END DO

    IF (j == 0) CALL finish(TRIM(routine), 'Variable '//TRIM(var%name)//' not found in any of children of tile '//&
      &                                    TRIM(tile%name)//' - this should not happen!')

    IF (ANY(fract_sum(:) > 1._wp)) CALL finish(TRIM(routine), 'fract_sum > 1 - this should not happen')

    DO ilev=1,SIZE(data,2)
      data_sum(:) = SUM(fractions(:,:) * data(:,ilev,:), DIM=2)
      WHERE (fract_sum(:) > EPSILON(1._wp))
        var%ptr(ics:ice,ilev,iblk) = data_sum(:) / fract_sum(:)
      ELSE WHERE
        var%ptr(ics:ice,ilev,iblk) = missval
      END WHERE
    END DO

    IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//')')

  END SUBROUTINE Aggregate_weighted_by_fract_3d

#endif
END MODULE mo_jsb_tile_class
