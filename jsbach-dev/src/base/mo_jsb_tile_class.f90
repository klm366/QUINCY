!>
!! @brief Contains tile class that contains the surface structure and memory state
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2012-09-24)
!!
MODULE mo_jsb_tile_class

  !USE mo_jsb_varlist_class, ONLY: t_jsb_varlist
  USE mo_jsb_memory_class,   ONLY: t_jsb_memory
  USE mo_exception,          ONLY: message_text, message, finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile
  PUBLIC :: new_tile, register_varlist

  INTEGER, PARAMETER :: max_no_of_sub_tiles = 30
  INTEGER, PARAMETER :: max_depth_of_tile_tree = 3
  INTEGER, PARAMETER :: max_no_of_components = 20

  TYPE t_jsb_tile
     CHARACTER(len=10)           :: name
     CHARACTER(len=50)           :: description
     INTEGER                     :: tree_path(1:max_depth_of_tile_tree)
     CHARACTER(len=20)           :: components(max_no_of_components)
     !INTEGER                     :: no_of_varlists = 0
     !TYPE(t_jsb_varlist), TARGET :: varlists(1:max_no_of_varlists)
     TYPE(t_jsb_memory), POINTER :: mem => NULL()
     TYPE(t_jsb_tile), POINTER   :: parent_tile => NULL()
     INTEGER                     :: no_of_sub_tiles = 0
     INTEGER                     :: idx_of_sub_tile = 0
     TYPE(t_jsb_tile), POINTER   :: sub_tiles(:) => NULL()
     INTEGER                     :: owner_model_id = 0
  END type t_jsb_tile

  INTERFACE new_tile
     MODULE PROCEDURE new_root_tile
     MODULE PROCEDURE new_sub_tile
  END INTERFACE new_tile

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_tile_class'

CONTAINS

  FUNCTION new_sub_tile(parent, name, description, components) RESULT(tile)

    TYPE(t_jsb_tile), POINTER, INTENT(inout) :: parent
    CHARACTER(len=*),          INTENT(in)    :: name
    CHARACTER(len=*),          INTENT(in)    :: description
    CHARACTER(len=*),          INTENT(in)    :: components(:)
    TYPE(t_jsb_tile), POINTER                :: tile

    INTEGER :: level
    CHARACTER(len=2*max_depth_of_tile_tree-1) :: path

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_sub_tile'

    WRITE(message_text,'(A)') 'Adding sub-tile "'//TRIM(name)//'" to tile "'//TRIM(parent%name)//'":'
    CALL message(TRIM(routine), message_text)
    CALL message(' -- Description', TRIM(description))

    ! If this is the first sub-tile for this parent, allocate memory for the sub-tiles
    IF (parent%no_of_sub_tiles == 0) ALLOCATE(parent%sub_tiles(max_no_of_sub_tiles))

    parent%no_of_sub_tiles = parent%no_of_sub_tiles + 1
    IF (parent%no_of_sub_tiles > max_no_of_sub_tiles) &
       CALL finish(TRIM(routine), 'Maximum number of subtiles reached.')

    tile => parent%sub_tiles(parent%no_of_sub_tiles)
    tile%idx_of_sub_tile = parent%no_of_sub_tiles
    tile%parent_tile => parent

    tile%name        = TRIM(name)
    tile%description = TRIM(description)
    tile%components(1:max_no_of_components) = '                    '
    tile%components(1:SIZE(components)) = components(:)

    path = ''
    DO level=1,max_depth_of_tile_tree
       IF (parent%tree_path(level) > 0) THEN
          tile%tree_path(level) = parent%tree_path(level)
          WRITE(path(2*level-1:2*level), '(I1,A1)') level, '-'
       ELSE
          tile%tree_path(level) = parent%no_of_sub_tiles
          WRITE(path(2*level-1:2*level-1), '(I1)') parent%no_of_sub_tiles
          EXIT
       END IF
    END DO
    CALL message(' -- Breadcrumb', TRIM(path))

    tile%owner_model_id = parent%owner_model_id

  END FUNCTION new_sub_tile

  FUNCTION new_root_tile(model_id, name, description, components) RESULT(tile)

    INTEGER,                  INTENT(in)    :: model_id
    CHARACTER(len=*),         INTENT(in)    :: name
    CHARACTER(len=*),         INTENT(in)    :: description
    CHARACTER(len=*),         INTENT(in)    :: components(:)
    TYPE(t_jsb_tile), POINTER               :: tile

    TYPE(t_jsb_tile), POINTER :: root

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_root_tile'

    CALL message(TRIM(routine), 'Adding root tile "'//TRIM(name)//'":')
    CALL message(' -- ', TRIM(description))

    ALLOCATE(root)
    !ALLOCATE(root%sub_tiles(max_no_of_sub_tiles))

    root%name           = TRIM(name)
    root%description    = TRIM(description)
    root%components(1:max_no_of_components) = '                    '
    root%components(1:SIZE(components)) = components(:)
    root%tree_path(:)   = 0
    root%tree_path(1)   = 1
    root%owner_model_id = model_id

    tile => root

  END FUNCTION new_root_tile

  SUBROUTINE register_varlist()

  END SUBROUTINE register_varlist

END MODULE mo_jsb_tile_class
