!>
!! @brief Contains the class and methods for the JSBACH memory structure
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_memory_class
#ifndef __NO_JSBACH__

!!$  USE mo_cdi_constants,     ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, &
!!$                                  ZA_SURFACE, &
!!$                                  FILETYPE_NC2, &
!!$                                  GRID_CELL, GRID_REFERENCE, GRID_UNSTRUCTURED_CELL
  USE mo_kind,              ONLY: wp
!!$  USE mo_exception,         ONLY: message
  USE mo_jsb_var_class, ONLY: t_jsb_var_p

!  USE mo_jsb_process_factory, ONLY: Create_process_memory

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_memory_p, t_jsb_memory, t_jsb_param

  TYPE, ABSTRACT :: t_jsb_param
  END TYPE t_jsb_param

  ! This type contains the memory state for each tile
  ! It includes all memory types from all sub-modules
  TYPE, ABSTRACT :: t_jsb_memory
    INTEGER :: id = 0
    INTEGER :: grid_id = 0
    INTEGER :: owner_model_id = 0
    INTEGER :: max_no_of_vars = 0
    INTEGER :: no_of_vars = 0
    CHARACTER(LEN=50) :: varlist_name
    TYPE(t_jsb_var_p), ALLOCATABLE :: vars(:)
    CLASS(t_jsb_param), POINTER :: param => NULL()
  CONTAINS
    PROCEDURE (Init_memory), DEFERRED, PASS(mem) :: Init
  END type t_jsb_memory

  TYPE t_jsb_memory_p
    CLASS(t_jsb_memory), POINTER :: p => NULL()
  END TYPE t_jsb_memory_p

  ABSTRACT INTERFACE
    SUBROUTINE Init_memory(mem, prefix, suffix)
      IMPORT :: t_jsb_memory
      CLASS(t_jsb_memory), INTENT(inout) :: mem
      CHARACTER(LEN=*), INTENT(in) :: prefix, suffix
    END SUBROUTINE Init_memory
  END INTERFACE

!  INTERFACE t_jsb_memory
!    PROCEDURE Create_process_memory
!  END INTERFACE t_jsb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_memory_class'

#endif
END MODULE mo_jsb_memory_class
