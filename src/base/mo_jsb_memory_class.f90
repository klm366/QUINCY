!>
!! @brief Contains the class and methods for the JSBACH4 memory structure
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! Copied to SLM by Jan Engel (20160815)
!!
MODULE mo_jsb_memory_class
!#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_jsb_config_class,    ONLY: t_jsb_config
  USE mo_jsb_var_class,       ONLY: t_jsb_var_p, t_jsb_var, t_jsb_var_real2d, t_jsb_var_real3d, &
                                    t_jsb_struct_p, t_jsb_struct, t_jsb_n_struct_p, t_jsb_n_struct        ! *slmjsb4


!  USE mo_jsb_process_factory, ONLY: Create_process_memory

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_memory_p, t_jsb_memory, t_jsb_param

  TYPE, ABSTRACT :: t_jsb_param
  END TYPE t_jsb_param

  ! This type contains the memory state for each tile for one process
  TYPE, ABSTRACT :: t_jsb_memory
    INTEGER :: id = 0
    INTEGER :: grid_id = 0
    INTEGER :: owner_model_id = 0
    INTEGER :: owner_proc_id = 0
    INTEGER, ALLOCATABLE :: owner_tile_path(:)
    CHARACTER(len=10)    :: owner_tile_name
    INTEGER :: max_no_of_vars = 0
    INTEGER :: no_of_vars = 0
    INTEGER :: max_no_of_structs = 0                  ! *slmjsb4
    INTEGER :: no_of_structs = 0                      ! *slmjsb4
    INTEGER :: max_no_of_n_structs = 0                ! *slmjsb4    --> "n" is for "nested"
    INTEGER :: no_of_n_structs = 0                    ! *slmjsb4    --> "n" is for "nested"
    CHARACTER(LEN=50) :: varlist_name
    CHARACTER(LEN=50) :: structlist_name              ! *slmjsb4
    CHARACTER(LEN=50) :: n_structlist_name            ! *slmjsb4    --> "n" is for "nested"
    TYPE(t_jsb_var_p), ALLOCATABLE :: vars(:)
    TYPE(t_jsb_struct_p), ALLOCATABLE :: structs(:)       ! *slmjsb4
    TYPE(t_jsb_n_struct_p), ALLOCATABLE :: n_structs(:)   ! *slmjsb4  --> "n" is for "nested"
    CLASS(t_jsb_param), POINTER :: param => NULL()
  CONTAINS
    PROCEDURE (Init_memory), DEFERRED, PASS(mem) :: Init
    PROCEDURE                                    :: Get_var_real2d
    PROCEDURE                                    :: Get_var_real3d
  END type t_jsb_memory

  TYPE t_jsb_memory_p
    CLASS(t_jsb_memory), POINTER :: p => NULL()
  END TYPE t_jsb_memory_p

  ! Note: argument "mem" needs to have the TARGET attribute so that we can point to individual variables
  !       (of class t_jsb_var) contained in the memory structure.
  ABSTRACT INTERFACE
    SUBROUTINE Init_memory(mem, prefix, suffix, lct_ids, model_id, ics, ice, iblk)
      IMPORT :: t_jsb_memory, t_jsb_config
      CLASS(t_jsb_memory), INTENT(inout), TARGET :: mem
      CHARACTER(LEN=*),    INTENT(in)    :: prefix, suffix
      INTEGER,             INTENT(in)    :: lct_ids!(:)
      INTEGER,             INTENT(in)    :: model_id
      INTEGER,             INTENT(in)    :: ics   ! model%options -> for SLM compatibility  ! *slmtemp
      INTEGER,             INTENT(in)    :: ice   ! model%options -> for SLM compatibility  ! *slmtemp
      INTEGER,             INTENT(in)    :: iblk  ! model%options -> for SLM compatibility  ! *slmtemp
    END SUBROUTINE Init_memory
  END INTERFACE

!  INTERFACE t_jsb_memory
!    PROCEDURE Create_process_memory
!  END INTERFACE t_jsb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_memory_class'

CONTAINS

  FUNCTION Get_var_real2d(this, name) RESULT(var_ptr)

    CLASS(t_jsb_memory),     INTENT(in) :: this
    CHARACTER(len=*),        INTENT(in) :: name
    CLASS(t_jsb_var_real2d), POINTER    :: var_ptr

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_var_real2d'

    IF (.NOT. ASSOCIATED(this%vars(1)%p)) CALL finish(TRIM(routine), 'Invalid vector of vars')

    var_ptr => NULL()
    DO i=1,this%no_of_vars
      IF (TRIM(this%vars(i)%p%name) == TRIM(name)) THEN
        SELECT TYPE (v => this%vars(i)%p)
        TYPE IS (t_jsb_var_real2d)
          var_ptr => v
        END SELECT
        EXIT
      END IF
    END DO
      
    ! Note: Not finding the variable in this memory, i.e. returning NULL(), should only happen because
    !       the variable belongs to an LCT that's not handled by this tile.

    !IF (.NOT. ASSOCIATED(var_ptr)) &
    !  & CALL finish(TRIM(routine), 'Variable '//TRIM(name)//' on tile '//TRIM(this%owner_tile_name)//' not found')

  END FUNCTION Get_var_real2d

  FUNCTION Get_var_real3d(this, name) RESULT(var_ptr)

    CLASS(t_jsb_memory),     INTENT(in) :: this
    CHARACTER(len=*),        INTENT(in) :: name
    CLASS(t_jsb_var_real3d), POINTER    :: var_ptr

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_var_real3d'

    var_ptr => NULL()
    DO i=1,this%no_of_vars
      IF (TRIM(this%vars(i)%p%name) == TRIM(name)) THEN
        SELECT TYPE (v => this%vars(i)%p)
        TYPE IS (t_jsb_var_real3d)
          var_ptr => v
        END SELECT
        EXIT
      END IF
    END DO
      
    ! Note: Not finding the variable in this memory, i.e. returning NULL(), should only happen because
    !       the variable belongs to an LCT that's not handled by this tile.

    !IF (.NOT. ASSOCIATED(var_ptr)) &
    !  & CALL finish(TRIM(routine), 'Variable '//TRIM(name)//' on tile '//TRIM(this%owner_tile_name)//' not found')

  END FUNCTION Get_var_real3d

!#endif
END MODULE mo_jsb_memory_class
