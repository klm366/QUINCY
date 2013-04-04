MODULE mo_jsb_class

  USE mo_jsb_model_class, ONLY: t_jsb_model_m

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach!, get_model

  TYPE t_jsb
     INTEGER :: no_of_models
     TYPE(t_jsb_model_m), POINTER :: models(:)
  END TYPE t_jsb

  TYPE(t_jsb) :: jsbach

!!$  INTERFACE get_model
!!$     MODULE PROCEDURE get_model_by_id
!!$     !MODULE PROCEDURE get_model_by_name
!!$  END INTERFACE get_model

!!$CONTAINS

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

!!$  FUNCTION get_model_by_id(id) RESULT(model)
!!$
!!$    USE mo_jsb_model_class, ONLY: t_jsb_model
!!$
!!$    INTEGER, INTENT(in)        :: id
!!$    TYPE(t_jsb_model), POINTER :: model
!!$
!!$    model => jsbach%models(id)%m
!!$
!!$  END FUNCTION get_model_by_id

END MODULE mo_jsb_class
