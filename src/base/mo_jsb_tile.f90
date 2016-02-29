!>
!! @brief Contains concrete tile class that contains the surface structure and memory state.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! This had to be separated from mo_jsb_tile_class in order to avoid circular dependencies.
!!
MODULE mo_jsb_tile
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message_text, message, finish

  USE mo_hsm_class,           ONLY: t_Hsm, t_State, t_Message
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_process_class,   ONLY: t_jsb_process_p, &
    &                               SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_
  USE mo_jsb_process_factory, ONLY: max_no_of_processes
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task
  USE mo_jsb_config_class,    ONLY: t_jsb_config_p

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile, t_jsb_tile_p

  TYPE, EXTENDS(t_jsb_tile_abstract) :: t_jsb_tile
  CONTAINS
    PROCEDURE :: Print   => Print_tile
    PROCEDURE :: Handler => Process_task_on_tile
  END TYPE t_jsb_tile

  TYPE t_jsb_tile_p
    TYPE(t_jsb_tile), POINTER :: p
  END TYPE t_jsb_tile_p

  ! Derived-type constructor
  INTERFACE t_jsb_tile
    PROCEDURE Construct_tile
  END INTERFACE t_jsb_tile

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_tile'

CONTAINS

    FUNCTION Construct_tile(name, description, processes, varlist_name, process_actions, &
    & parent, model_id, grid_id) RESULT(return_ptr)

    USE mo_jsb_process_factory,   ONLY: Create_process_memory, Get_process_name, SRF, ATM2LAND, LAND2ATM
    USE mo_jsb_memory_class, ONLY: t_jsb_memory
    USE mo_util,            ONLY: one_of

    CHARACTER(len=*),  INTENT(in)           :: name
    CHARACTER(len=*),  INTENT(in)           :: description
    INTEGER,           INTENT(in)           :: processes(:)
    CHARACTER(len=*),  INTENT(in)           :: varlist_name
    INTEGER,           INTENT(in), OPTIONAL :: process_actions(:)
    CLASS(t_jsb_tile), POINTER,    OPTIONAL :: parent
    INTEGER,           INTENT(in), OPTIONAL :: model_id   ! Only needed if parent is not present
    INTEGER,           INTENT(in), OPTIONAL :: grid_id    ! Only needed if parent is not present
    CLASS(t_jsb_tile), POINTER              :: return_ptr

    CHARACTER(LEN=15)                     :: process_name
    LOGICAL                               :: must_create_memory_parents
    INTEGER                               :: i, j, iproc
    CLASS(t_jsb_memory),          POINTER :: tmp_ptr
    CLASS(t_jsb_tile_abstract),   POINTER :: tmp_tile
    CLASS(t_State),               POINTER :: parent_tmp

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Construct_tile'

    ! Allocate and initialize base structure (=state) for tile
    ALLOCATE(t_jsb_tile::return_ptr)
    IF (PRESENT(parent)) THEN
      parent_tmp => parent
      CALL return_ptr%Init_state(name, description, parent_tmp)
    ELSE
      CALL return_ptr%Init_state(name, description)
    END IF

    ! Initialize further components of t_jsb_tile_abstract extended type

    IF (PRESENT(process_actions)) THEN
      IF (SIZE(processes) /= SIZE(process_actions)) &
        & CALL finish(TRIM(routine), 'Size of processes and process_actions must be equal!')
    END IF

    IF (PRESENT(parent)) THEN
      SELECT TYPE (parent)
      CLASS IS (t_jsb_tile)
        return_ptr%parent_tile => parent
        return_ptr%grid_id = parent%grid_id
        return_ptr%owner_model_id = parent%owner_model_id
      END SELECT
    ELSE
      IF (PRESENT(model_id)) THEN
        return_ptr%owner_model_id = model_id
      END IF
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
          return_ptr%process_action(iproc) = process_actions(j)
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
!!$      ELSE IF (ASSOCIATED(return_ptr%prev_sibling)) THEN
!!$        iproc = i
!!$        SELECT TYPE (sibling => return_ptr%prev_sibling)
!!$        CLASS IS (t_jsb_tile_abstract)
!!$          return_ptr%process_action(iproc) = sibling%process_action(iproc)
!!$        END SELECT
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
            !return_ptr%parent_tile%process_action(iproc) = AGGREGATE_
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

    CLASS(t_jsb_tile), POINTER,  INTENT(inout) :: tile
    CHARACTER(len=*),            INTENT(in)    :: varlist_name
    CHARACTER(len=*),            INTENT(in)    :: prefix        ! Prefix for variable names
    CHARACTER(len=*),            INTENT(in)    :: suffix
    INTEGER,                     INTENT(in)    :: grid_id

    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    CHARACTER(len=10)          :: parent_tile_name
    INTEGER                    :: table

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

    CLASS(t_jsb_tile),  INTENT(inout) :: this
    CLASS(t_Hsm),       INTENT(in)    :: hsm
    CLASS(t_Message),   INTENT(in)    :: msg_in
    CLASS(t_Message),   POINTER       :: return_ptr

    CLASS(t_Message), POINTER :: msg
    INTEGER :: iproc

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Process_task_on_tile'

    SELECT TYPE(model=>hsm)
    CLASS IS (t_jsb_model)
      IF (model%options%iblk == 1) THEN
        CALL message(TRIM(routine),'Started')
        CALL message('    ', 'Processing task '//TRIM(msg_in%name)//'/action '//&
          &                   TRIM(msg_in%action%name)//' on tile '//TRIM(this%name))
        CALL this%PRINT()
      END IF
    END SELECT

    ALLOCATE(msg, source=msg_in)
    return_ptr => NULL()

    SELECT TYPE(model=>hsm)
    CLASS IS (t_jsb_model)
      SELECT TYPE(task=>msg)
      CLASS IS (t_jsb_process_task)
        SELECT CASE (TRIM(task%action%name))
        CASE ('START')
          IF (this%Has_children()) THEN

          END IF
        CASE ('INTEGRATE')
          IF (this%process_action(task%process_id) == ON_LEAFS_) THEN
            IF (.NOT. this%Has_children()) THEN
              CALL task%Do_it(this, model%Get_process_configs(), model%options)
            END IF
            CALL model%Goto(model%Get_next_tile())
            !model%next => model%Get_next_tile()
            IF (.NOT. ASSOCIATED(model%next)) THEN
              return_ptr => NULL()
            ELSE
              return_ptr => task
            END IF
          ELSE IF (this%process_action(task%process_id) == AGGREGATE_) THEN
            IF (this%Has_children()) THEN
              CALL model%Goto(model%Get_next_tile())
              !model%next => model%Get_next_tile()
              IF (.NOT. ASSOCIATED(model%next)) THEN
                return_ptr => NULL()
              ELSE
                return_ptr => task
              END IF
            END IF
          END IF
        END SELECT
      END SELECT
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Error running task')
    END SELECT

    SELECT TYPE(model=>hsm)
    CLASS IS (t_jsb_model)
      IF (model%options%iblk == 1) THEN
        CALL message(TRIM(routine),'Finished')
      END IF
    END SELECT

  END FUNCTION Process_task_on_tile

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
END MODULE mo_jsb_tile
