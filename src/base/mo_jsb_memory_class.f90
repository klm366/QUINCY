!>
!! @brief Contains the class and methods for the JSBACH memory structure
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

MODULE mo_jsb_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: finish
  USE mo_jsb_config_class,  ONLY: t_jsb_config
  USE mo_jsb_var_class,     ONLY: t_jsb_var_p, t_jsb_var, t_jsb_var_real2d, t_jsb_var_real3d, REAL2D, REAL3D
  USE mo_jsb_grid_class,    ONLY: t_jsb_vgrid
  USE mo_util,              ONLY: int2string

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
    CHARACTER(LEN=:), ALLOCATABLE :: owner_proc_name
    INTEGER :: no_of_children = 0
    TYPE(t_jsb_memory_p), POINTER :: children(:) => NULL()
    CLASS(t_jsb_memory) , POINTER :: parent => NULL()
    INTEGER, ALLOCATABLE :: owner_tile_path(:)
    CHARACTER(len=10)    :: owner_tile_name
    INTEGER :: max_no_of_vars = 0
    INTEGER :: no_of_vars = 0
    INTEGER :: no_of_leaves = 0
    INTEGER :: pos_in_leaves = 0
    CHARACTER(LEN=:), ALLOCATABLE :: varlist_name
    TYPE(t_jsb_var_p), ALLOCATABLE :: vars(:)
    CLASS(t_jsb_param), POINTER :: param => NULL()
    LOGICAL :: in_var_groups = .FALSE.
    LOGICAL :: single_tile_output = .TRUE.
    TYPE(t_jsb_vgrid), POINTER :: vgrid_t       => NULL()         ! "vertical" grid type for tiled output, i.e. single_tile_output=.FALSE.
  CONTAINS
    PROCEDURE (Init_memory), DEFERRED, PASS(mem) :: Init
    PROCEDURE                                    :: Get_var
    PROCEDURE                                    :: Get_var_position
    PROCEDURE                                    :: Collect_var_real2d  ! Collect var from its child tiles
                                                                        ! R: Note, we do not want to keep a "var collection" in the memory
                                                                        !    because we already have all child values in the memory.
    PROCEDURE                                    :: Collect_var_real3d  ! Collect var from its child tiles
    PROCEDURE                                    :: HandDown_var_real2d ! Hand back var to child tiles
    PROCEDURE                                    :: Average_over_children_real2d ! Average var over the child tiles
  END TYPE t_jsb_memory

  TYPE t_jsb_memory_p
    CLASS(t_jsb_memory), POINTER :: p => NULL()
  END TYPE t_jsb_memory_p

  ! Note: argument "mem" needs to have the TARGET attribute so that we can point to individual variables
  !       (of class t_jsb_var) contained in the memory structure.
  ABSTRACT INTERFACE
    SUBROUTINE Init_memory(mem, prefix, suffix, lct_ids, model_id)
      IMPORT :: t_jsb_memory, t_jsb_config
      CLASS(t_jsb_memory), INTENT(inout), TARGET :: mem
      CHARACTER(LEN=*),    INTENT(in)    :: prefix, suffix
      INTEGER,             INTENT(in)    :: lct_ids(:)
      INTEGER,             INTENT(in)    :: model_id
    END SUBROUTINE Init_memory
  END INTERFACE

!  INTERFACE t_jsb_memory
!    PROCEDURE Create_process_memory
!  END INTERFACE t_jsb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_memory_class'

CONTAINS

  FUNCTION Get_var(this, name) RESULT(var_ptr)

    CLASS(t_jsb_memory),     INTENT(in) :: this
    CHARACTER(len=*),        INTENT(in) :: name
    CLASS(t_jsb_var),        POINTER    :: var_ptr

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_var'

    IF (.NOT. ASSOCIATED(this%vars(1)%p)) CALL finish(TRIM(routine), 'Invalid vector of vars')

    var_ptr => NULL()
    DO i=1,this%no_of_vars
      IF (TRIM(this%vars(i)%p%name) == TRIM(name)) THEN
        var_ptr => this%vars(i)%p
        EXIT
      END IF
    END DO

    ! Note: Not finding the variable in this memory, i.e. returning NULL(), should only happen because
    !       the variable belongs to an LCT that's not handled by this tile.

    !IF (.NOT. ASSOCIATED(var_ptr)) &
    !  & CALL finish(TRIM(routine), 'Variable '//TRIM(name)//' on tile '//TRIM(this%owner_tile_name)//' not found')

  END FUNCTION Get_var

  FUNCTION Get_var_position(this, name) RESULT(pos)

    CLASS(t_jsb_memory),     INTENT(in) :: this
    CHARACTER(len=*),        INTENT(in) :: name
    INTEGER                             :: pos

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_var_position'

    IF (.NOT. ASSOCIATED(this%vars(1)%p)) CALL finish(TRIM(routine), 'Invalid vector of vars')

    pos = 0
    DO i=1,this%no_of_vars
      IF (TRIM(this%vars(i)%p%name) == TRIM(name)) THEN
        pos = i
        EXIT
      END IF
    END DO

    ! Note: Not finding the variable in this memory, i.e. returning 0, should only happen because
    !       the variable belongs to an LCT that's not handled by this tile.

    !IF (.NOT. ASSOCIATED(var_ptr)) &
    !  & CALL finish(TRIM(routine), 'Variable '//TRIM(name)//' on tile '//TRIM(this%owner_tile_name)//' not found')

  END FUNCTION Get_var_position

  FUNCTION Collect_var_real2d(this, ics, ice, iblk, no_children,  var) RESULT(var_collected)

    CLASS(t_jsb_memory),    INTENT(in)     :: this        ! mem_process of the parent tile
    INTEGER,                INTENT(in)     :: ics, ice, iblk
    INTEGER,                INTENT(in)     :: no_children ! of the parent tile
    TYPE(t_jsb_var_real2d), INTENT(in)     :: var         ! of the parent tile: tile%mem(iproc)%p%var
    REAL(wp), POINTER                      :: var_collected(:,:)

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Collect_var_real2d'

    ALLOCATE(var_collected(ice-ics+1,no_children))

    DO i=1,no_children
       var_collected(:,i) = this%children(i)%p%vars(var%child_idx(i))%p%ptr2d(ics:ice,iblk)
    END DO

  END FUNCTION Collect_var_real2d


  FUNCTION Collect_var_real3d(this, ics, ice, iblk, no_children,  var) RESULT(var_collected)

    CLASS(t_jsb_memory),    INTENT(in)     :: this        ! mem_process of the parent tile
    INTEGER,                INTENT(in)     :: ics, ice, iblk
    INTEGER,                INTENT(in)     :: no_children ! of the parent tile
    TYPE(t_jsb_var_real3d), INTENT(in)     :: var         ! of the parent tile: tile%mem(iproc)%p%var
    REAL(wp), POINTER                      :: var_collected(:,:,:)

    INTEGER :: i !, iproc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Collect_var_real3d'

    ALLOCATE(var_collected(ice - ics + 1, &
      &                    SIZE(this%children(1)%p%vars(var%child_idx(1))%p%ptr3d, DIM=2), &
      &                    no_children))

    DO i=1,no_children
       ! R: Unterschied zu %ptr(ics:ice,iblk) ???
       var_collected(:,:,i) = this%children(i)%p%vars(var%child_idx(i))%p%ptr3d(ics:ice,:, iblk)
    END DO

  END FUNCTION Collect_var_real3d


  SUBROUTINE HandDown_var_real2d(this, ics, ice, iblk, no_children,  var, var_handdown)

    CLASS(t_jsb_memory), INTENT(inout)  :: this        ! mem_process of the parent tile: tile%mem(iproc)%p
    INTEGER,             INTENT(in)     :: ics, ice, iblk
    INTEGER,             INTENT(in)     :: no_children ! of the parent tile
    CLASS(t_jsb_var),    INTENT(in)     :: var         ! of the parent tile only to find child memories
    REAL(wp),            INTENT(in)     :: var_handdown(:,:) ! e.g. corresponds to (ics:ice,no_children)

    INTEGER :: i !, iproc

    CHARACTER(len=*), PARAMETER :: routine = modname//':HandDown_var_real2d'

    DO i=1,no_children
      this%children(i)%p%vars(var%child_idx(i))%p%ptr2d(ics:ice, iblk) = var_handdown(:,i)
    END DO

  END SUBROUTINE HandDown_var_real2d


  FUNCTION Average_over_children_real2d(this, ics, ice, iblk, no_children,  var, fract_st) RESULT(var_averaged)

    CLASS(t_jsb_memory), INTENT(in)     :: this          ! mem_process of the parent tile: tile%mem(iproc)%p
    INTEGER,             INTENT(in)     :: ics, ice, iblk
    INTEGER,             INTENT(in)     :: no_children   ! of the parent tile
    CLASS(t_jsb_var),    INTENT(in)     :: var           ! of the parent tile: tile%mem(iproc)%p%var
                                                         ! R: kann man auf TYPE(t_jsb_var_real2d) setzen
    REAL(wp),            INTENT(in)     :: fract_st(:,:) ! fract_st(:,:) =tile%aggregators(1)%p%fractions(ics:ice,iblk,:)
    REAL(wp), allocatable               :: var_averaged(:)

    ! local variables
    INTEGER :: i
    REAL(wp), allocatable               :: var_collected(:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':Average_over_children'

    ALLOCATE(var_collected(ice-ics+1,no_children))
    ! Collect corresponding vars of subtiles
    DO i=1,no_children
      IF ( .NOT. ASSOCIATED(this%children(i)%p%vars(var%child_idx(i))%p) ) THEN
        CALL finish(TRIM(routine), 'This var does not exist on the child tile '//TRIM(int2string(i))//' .')
      END IF
      var_collected(:,i) = this%children(i)%p%vars(var%child_idx(i))%p%ptr2d(ics:ice, iblk)
    END DO
    ! Note, this works because the fract_st sum over DIM=2 is always = 1 !
    var_averaged(:) = SUM( var_collected(:,:) * fract_st(:,:) , DIM=2)
    DEALLOCATE(var_collected)

  END FUNCTION Average_over_children_real2d

#endif
END MODULE mo_jsb_memory_class
