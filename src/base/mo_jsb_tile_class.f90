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

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_io_units,            ONLY: filename_max

  USE mo_hsm_class,           ONLY: t_Hsm, t_State, t_Message
  USE mo_jsb_memory_class,    ONLY: t_jsb_memory_p, t_jsb_memory
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct
  USE mo_jsb_config_class,    ONLY: t_jsb_config_p

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile_abstract

!!$  INTEGER, PARAMETER :: max_depth_of_tile_tree = 3

  TYPE, ABSTRACT, EXTENDS(t_State) :: t_jsb_tile_abstract
    TYPE(t_jsb_lct),        ALLOCATABLE    :: lcts(:)          !< List of land cover types for tile
    INTEGER,                ALLOCATABLE    :: process_action(:)
    REAL(wp),                   POINTER    :: fract     (:,:)  !< Cover fraction of this tile (relative to parent tile)
    REAL(wp),                   POINTER    :: fract_max (:,:)  !< Maximum potential cover fraction for this tile (rel. to parent)
    CHARACTER(LEN=filename_max)            :: fract_filename   !< Name of file to read tile fractions from (all tiles except root)
                                                               !! The variable names must be "fract[_max]_"//tile%name
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

#endif
END MODULE mo_jsb_tile_class
