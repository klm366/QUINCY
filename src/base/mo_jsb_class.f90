MODULE mo_jsb_class
#ifndef __NO_JSBACH__

  USE mo_exception,       ONLY: finish

  USE mo_jsb_model_class, ONLY: t_jsb_model, t_jsb_model_m
!!$  USE mo_jsb_tile,  ONLY: t_jsb_tile

  IMPLICIT NONE
  PRIVATE

  ! t_jsb needs to be public for the NEC compiler
  PUBLIC :: t_jsb, jsbach, Get_model !, Get_top_tile, Get_next_tile

  TYPE t_jsb
     INTEGER :: no_of_models
     TYPE(t_jsb_model_m), POINTER :: models(:)
  END TYPE t_jsb

  TYPE(t_jsb) :: jsbach

  INTERFACE get_model
     MODULE PROCEDURE get_model_by_id
     !MODULE PROCEDURE get_model_by_name
  END INTERFACE get_model

CONTAINS

!!$  FUNCTION get_grid_by_id(id) RESULT(grid)
!!$
!!$    USE mo_jsb_grid_class, ONLY: t_jsb_grid
!!$
!!$    INTEGER, INTENT(in)        :: id
!!$    TYPE(t_jsb_grid), POINTER :: grid
!!$
!!$    grid => jsbach%models(id)%m%grid
!!$
!!$  END FUNCTION get_grid_by_id

  FUNCTION get_model_by_id(id) RESULT(model)

    INTEGER, INTENT(in)        :: id
    TYPE(t_jsb_model), POINTER :: model

    model => jsbach%models(id)%m

  END FUNCTION get_model_by_id

!!$  FUNCTION Get_top_tile(model) RESULT(return_ptr)
!!$
!!$    CLASS(t_jsb_model), INTENT(in) :: model
!!$    TYPE(t_jsb_tile),   POINTER    :: return_ptr
!!$
!!$    return_ptr => NULL()
!!$
!!$    SELECT TYPE (tile => model%top)
!!$    CLASS IS (t_jsb_tile)
!!$      return_ptr => tile
!!$    END SELECT
!!$
!!$  END FUNCTION Get_top_tile
!!$
!!$  FUNCTION Get_next_tile(model) RESULT(return_ptr)
!!$
!!$    USE mo_hsm_class, ONLY: t_State
!!$
!!$    CLASS(t_jsb_model), INTENT(inout) :: model
!!$    TYPE(t_jsb_tile), POINTER     :: return_ptr
!!$
!!$    CLASS(t_State), POINTER :: next
!!$
!!$    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_class:Get_next_tile'
!!$
!!$    next => model%Goto_next()
!!$    IF (.NOT. ASSOCIATED(next)) THEN
!!$      return_ptr => NULL()
!!$      RETURN
!!$    END IF
!!$
!!$    SELECT TYPE (next)
!!$    TYPE IS (t_jsb_tile)
!!$      return_ptr => next
!!$    CLASS IS (t_State)
!!$      CALL finish(TRIM(routine), 'current tile is of type t_State, should be t_jsb_tile')
!!$    CLASS DEFAULT
!!$      CALL finish(TRIM(routine), 'Unkown type for tile')
!!$    END SELECT
!!$
!!$    NULLIFY(next)
!!$
!!$  END FUNCTION Get_next_tile

#endif
END MODULE mo_jsb_class
