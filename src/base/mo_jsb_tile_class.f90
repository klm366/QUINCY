!>
!! @brief Contains tile class that contains the surface structure and memory state
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_tile_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  !USE mo_jsb_varlist_class, ONLY: t_jsb_varlist
  USE mo_hsm_class,           ONLY: t_Hsm, t_State, t_Message, Handler_interface
  USE mo_jsb_memory_class,    ONLY: t_jsb_memory_p, t_jsb_memory
  USE mo_jsb_process_factory,       ONLY: max_no_of_processes
  USE mo_jsb_process_class,   ONLY: SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_, &
                                    t_jsb_process_task
  USE mo_exception,           ONLY: message_text, message, finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile

  INTEGER, PARAMETER :: max_depth_of_tile_tree = 3

  TYPE, EXTENDS(t_State) :: t_jsb_tile
    INTEGER                          :: lct_type
    INTEGER                          :: process_action(max_no_of_processes)
    REAL(wp),                POINTER :: fract     (:,:)  !< Cover fraction of this tile (relative to parent tile)
    REAL(wp),                POINTER :: fract_max (:,:)  !< Maximum potential cover fraction for this tile (rel. to parent)
    !CLASS(t_jsb_param_p),    POINTER :: param(:)
    CLASS(t_jsb_memory_p),   POINTER :: mem(:)
    TYPE(t_jsb_tile),        POINTER :: parent_tile => NULL()
    INTEGER                          :: grid_id = -1
    INTEGER                          :: owner_model_id = 0
  CONTAINS
!!$    PROCEDURE :: Get_memory
    PROCEDURE :: Get_parent_memory
    PROCEDURE :: Is_process_active
    PROCEDURE :: Set_process_action_parents
    PROCEDURE :: Print => Print_tile
    PROCEDURE :: Handler => Process_task_on_tile
  END type t_jsb_tile

  ! Derived-type constructor
  INTERFACE t_jsb_tile
    PROCEDURE Construct_tile
  END INTERFACE t_jsb_tile

  !>
  !! Container CLASS for tiles
  !!
  !! Contains pointer to an instance of type t_jsb_tile in order to allow arrays
  !! of pointers to tiles.
  !!
  TYPE t_jsb_tile_p
    TYPE(t_jsb_tile), POINTER :: p
  END TYPE t_jsb_tile_p

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_tile_class'

CONTAINS

  FUNCTION Construct_tile(name, description, processes, varlist_name, process_actions, &
    & parent, model_id, grid_id) RESULT(return_ptr)

!!$    USE mo_hsm_class, ONLY: Handler_interface
    USE mo_jsb_process_factory,   ONLY: Create_process_memory, Get_process_name, SRF, ATM2LAND, LAND2ATM
    USE mo_jsb_memory_class, ONLY: t_jsb_memory
    USE mo_util,            ONLY: one_of

    CHARACTER(len=*), INTENT(in)           :: name
    CHARACTER(len=*), INTENT(in)           :: description
    INTEGER,          INTENT(in)           :: processes(:)
    CHARACTER(len=*), INTENT(in)           :: varlist_name
    INTEGER,          INTENT(in), OPTIONAL :: process_actions(:)
    CLASS(t_State),   POINTER,    OPTIONAL :: parent
    INTEGER,          INTENT(in), OPTIONAL :: model_id   ! Only needed if parent is not present
    INTEGER,          INTENT(in), OPTIONAL :: grid_id    ! Only needed if parent is not present
    TYPE(t_jsb_tile), POINTER              :: return_ptr

    CLASS(t_jsb_memory), POINTER :: tmp_ptr
    TYPE(t_jsb_tile),    POINTER :: tmp_tile
    CHARACTER(LEN=15) :: process_name
    LOGICAL :: must_create_memory_parents
    INTEGER :: i, j, iproc
!!$    PROCEDURE(Handler_interface), POINTER :: task_handler

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Construct_tile'

    ! Allocate and initialize base structure (=state) for tile
    ALLOCATE(return_ptr)
!!$    task_handler => Process_task_on_tile
!!$    CALL return_ptr%Init(name, description, parent, task_handler)
    CALL return_ptr%Init_state(name, description, parent)
!!$    NULLIFY(task_handler)

    ! Initialize further components of t_jsb_tile extended type

    IF (PRESENT(process_actions)) THEN
      IF (SIZE(processes) /= SIZE(process_actions)) &
        & CALL finish(TRIM(routine), 'Size of processes and process_actions must be equal!')
    END IF

    IF (PRESENT(parent)) THEN
      SELECT TYPE (parent)
      TYPE IS (t_jsb_tile)
        return_ptr%parent_tile => parent
        return_ptr%grid_id = parent%grid_id
        return_ptr%owner_model_id = parent%owner_model_id
      END SELECT
    ELSE
      IF (PRESENT(model_id)) return_ptr%owner_model_id = model_id
      IF (PRESENT(grid_id))  return_ptr%grid_id = grid_id
    END IF

    IF (ASSOCIATED(return_ptr%parent_tile)) THEN
      CALL Init_tile(return_ptr, TRIM(varlist_name), '', TRIM(name), return_ptr%grid_id)
    ELSE
      CALL Init_tile(return_ptr, TRIM(varlist_name), '', '', return_ptr%grid_id)
    END IF

    ! Allocate all memory components for new tile and initialize their data part with NULL()
    ALLOCATE(return_ptr%mem(max_no_of_processes))
    DO iproc=1,max_no_of_processes
      return_ptr%mem(iproc)%p => NULL()
    END DO

    ! Loop over all processes and set actions for this tile
    DO i=1,max_no_of_processes
      return_ptr%process_action(i) = SKIP_
      must_create_memory_parents = .FALSE.
      iproc = -1
      j = one_of(i, processes)
      IF (j > 0) THEN
        iproc = processes(j)
        IF (PRESENT(process_actions)) THEN
          return_ptr%process_action(iproc) = process_actions(iproc)
        ELSE
          IF (iproc == ATM2LAND) THEN
            return_ptr%process_action(iproc) = ON_TILE_
          ELSE IF (iproc == LAND2ATM) THEN
            return_ptr%process_action(iproc) = ON_TILE_
          ELSE
            return_ptr%process_action(iproc) = ON_LEAFS_
          END IF
        END IF
        ! All parents of this tile aggregate memory for this new process
        CALL return_ptr%Set_process_action_parents(iproc, AGGREGATE_)
        must_create_memory_parents = .TRUE.
      ELSE IF (ASSOCIATED(return_ptr%parent_tile)) THEN
        IF (return_ptr%parent_tile%process_action(i) > AGGREGATE_) THEN
          iproc = i
          SELECT CASE (return_ptr%parent_tile%process_action(iproc))
          CASE (SKIP_)
            return_ptr%process_action(iproc) = SKIP_
          CASE (INHERIT_)
            return_ptr%process_action(iproc) = INHERIT_
          CASE (ON_SUBTREE_)
            return_ptr%process_action(iproc) = ON_SUBTREE_
          CASE (ON_LEAFS_)
            return_ptr%process_action(iproc) = ON_LEAFS_
            return_ptr%parent_tile%process_action(iproc) = AGGREGATE_
          CASE (ON_TILE_)
            IF (iproc == LAND2ATM) THEN
              return_ptr%process_action(iproc) = SKIP_
            ELSE
              return_ptr%process_action(iproc) = INHERIT_
            END IF
          END SELECT
        END IF
      END IF

      IF (iproc < 1) CYCLE
      IF (return_ptr%process_action(iproc) == SKIP_) CYCLE

      process_name = Get_process_name(iproc)

      IF (return_ptr%process_action(iproc) == INHERIT_) THEN
        return_ptr%mem(iproc)%p => return_ptr%parent_tile%mem(iproc)%p
      ELSE
        !return_ptr%param(iproc)%p => t_jsb_param(iproc)
        tmp_ptr => Create_process_memory(iproc)
        SELECT TYPE (tmp_ptr)
        CLASS IS (t_jsb_memory)
          return_ptr%mem(iproc)%p => tmp_ptr
          return_ptr%mem(iproc)%p%id = iproc
          return_ptr%mem(iproc)%p%varlist_name = TRIM(varlist_name)//'_'//TRIM(process_name)
          return_ptr%mem(iproc)%p%grid_id = return_ptr%grid_id
          return_ptr%mem(iproc)%p%owner_model_id = return_ptr%owner_model_id
          return_ptr%mem(iproc)%p%owner_proc_id = iproc
          ALLOCATE(return_ptr%mem(iproc)%p%owner_tile_path(SIZE(return_ptr%path)))
          return_ptr%mem(iproc)%p%owner_tile_path(:) = return_ptr%path(:)
          IF (ASSOCIATED(return_ptr%parent_tile)) THEN
            CALL return_ptr%mem(iproc)%p%Init(TRIM(process_name), TRIM(name))
          ELSE
            CALL return_ptr%mem(iproc)%p%Init(TRIM(process_name), '')
          END IF
        END SELECT
        NULLIFY(tmp_ptr)
      END IF

      IF (must_create_memory_parents) THEN
        tmp_tile => return_ptr%parent_tile
        DO WHILE (ASSOCIATED(tmp_tile))
          IF (ASSOCIATED(tmp_tile%mem(iproc)%p)) CYCLE
          !return_ptr%param(iproc)%p => t_jsb_param(iproc)
          tmp_ptr => Create_process_memory(iproc)
          SELECT TYPE (tmp_ptr)
            CLASS IS (t_jsb_memory)
            tmp_tile%mem(iproc)%p => tmp_ptr
            tmp_tile%mem(iproc)%p%id = iproc
            tmp_tile%mem(iproc)%p%varlist_name = TRIM(varlist_name)//'_'//TRIM(process_name)
            tmp_tile%mem(iproc)%p%grid_id = return_ptr%grid_id
            tmp_tile%mem(iproc)%p%owner_model_id = return_ptr%owner_model_id
            IF (ASSOCIATED(tmp_tile%parent_tile)) THEN
              CALL tmp_tile%mem(iproc)%p%Init(TRIM(process_name), TRIM(tmp_tile%name))
            ELSE
              CALL tmp_tile%mem(iproc)%p%Init(TRIM(process_name), '')
            END IF
          END SELECT
          NULLIFY(tmp_ptr)
          tmp_tile => tmp_tile%parent_tile
        END DO
        NULLIFY(tmp_tile)
      END IF
    END DO

  END FUNCTION Construct_tile

  SUBROUTINE Init_tile(tile, varlist_name, prefix, suffix, grid_id)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, Get_varlist, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    TYPE(t_jsb_tile),   POINTER,  INTENT(inout) :: tile
    CHARACTER(len=*),             INTENT(in)    :: varlist_name
    CHARACTER(len=*),             INTENT(in)    :: prefix        ! Prefix for variable names
    CHARACTER(len=*),             INTENT(in)    :: suffix
    INTEGER,                      INTENT(in)    :: grid_id

    CHARACTER(len=10) :: parent_tile_name

    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_tile'

    hgrid   => Get_grid(grid_id)
    surface => Get_vgrid('surface')

    table = tables(1)

    parent_tile_name = tile%Get_parent_name()
    IF (parent_tile_name == '') parent_tile_name = 'grid box  '

    CALL jsb_add_var(TRIM(varlist_name), prefix//'fract'//suffix, tile%fract,                                 &
      & hgrid, surface,                                                                                         &
      & t_cf('land_area_fraction', '',                                                                          &
      &      'fraction of '//TRIM(parent_tile_name)//' covered by '//TRIM(tile%name)),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp)

    CALL jsb_add_var(TRIM(varlist_name), prefix//'fract_max'//suffix, tile%fract_max,                         &
      & hgrid, surface,                                                                                         &
      & t_cf('land_area_fraction', '',                                                                          &
      &      'fraction of '//TRIM(parent_tile_name)//' covered by '//TRIM(tile%name)),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp)

    tile%fract = 1._wp
    tile%fract_max = 1._wp

  END SUBROUTINE Init_tile

  ! Task handler function for tile
  FUNCTION Process_task_on_tile(this, hsm, msg_in) RESULT(return_ptr)

    CLASS(t_jsb_tile),     INTENT(inout) :: this
    CLASS(t_Hsm),       POINTER       :: hsm
    CLASS(t_Message),   POINTER       :: msg_in
    CLASS(t_Message),   POINTER       :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Process_task_on_tile'

    CALL message(TRIM(routine),'Started')

    ALLOCATE(t_jsb_process_task::return_ptr)

    CALL this%Print()

    CALL message(TRIM(routine),'Finished')

  END FUNCTION Process_task_on_tile

  FUNCTION Get_parent_memory(this, iproc) RESULT(return_ptr)

    CLASS(t_jsb_tile),  INTENT(in) :: this
    INTEGER,            INTENT(in) :: iproc
    CLASS(t_jsb_memory), POINTER   :: return_ptr

    CLASS(*), POINTER :: parent

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_parent_memory'

    parent => this%Get_parent()
    IF (ASSOCIATED(parent)) THEN
      SELECT TYPE (parent)
      TYPE IS (t_jsb_tile)
        return_ptr => parent%mem(iproc)%p
      END SELECT
    ELSE
      return_ptr => NULL()
    END IF

  END FUNCTION Get_parent_memory

  FUNCTION Is_process_active(this, iproc) RESULT(return_value)

    CLASS(t_jsb_tile), INTENT(in) :: this
    INTEGER,           INTENT(in) :: iproc
    LOGICAL                       :: return_value

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Is_process_active'

    return_value = .FALSE.

    SELECT TYPE (this)
    TYPE IS (t_jsb_tile)
      return_value = ASSOCIATED(this%mem(iproc)%p)
    END SELECT

  END FUNCTION Is_process_active

  RECURSIVE SUBROUTINE Set_process_action_parents(this, iproc, action)

    CLASS(t_jsb_tile), INTENT(inout) :: this
    INTEGER,           INTENT(in)    :: iproc
    INTEGER,           INTENT(in)    :: action

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Set_process_action_parents'

    IF (ASSOCIATED(this%parent_tile)) THEN
      this%parent_tile%process_action(iproc) = action
      CALL this%parent_tile%Set_process_action_parents(iproc, action)
    END IF

  END SUBROUTINE Set_process_action_parents

  SUBROUTINE Print_tile(this)

    USE mo_jsb_process_factory, ONLY: Get_process_name
    USE mo_jsb_process_class,   ONLY: Get_action_name

    CLASS(t_jsb_tile), INTENT(in) :: this

    INTEGER :: iproc

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Print_tile'

    ! Print for basic type
    CALL this%Print_state()

    ! Additional print for tile type
    DO iproc=1,max_no_of_processes
      IF (this%process_action(iproc) > 0) THEN
        CALL message('     Process '//Get_process_name(iproc), 'Action '//Get_action_name(this%process_action(iproc)))
      END IF
    END DO

  END SUBROUTINE Print_tile

#endif
END MODULE mo_jsb_tile_class
