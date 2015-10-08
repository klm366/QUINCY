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
  USE mo_hsm_class,           ONLY: t_State
  USE mo_jsb_memory_class,    ONLY: t_jsb_memory_p, t_jsb_memory
  USE mo_jsb_process_factory,       ONLY: max_no_of_processes
  USE mo_exception,           ONLY: message_text, message, finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile

  INTEGER, PARAMETER :: max_depth_of_tile_tree = 3

  TYPE, EXTENDS(t_State) :: t_jsb_tile
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

  FUNCTION Construct_tile(name, description, processes, varlist_name, parent, model_id, grid_id) RESULT(return_ptr)

    USE mo_jsb_process_factory,   ONLY: Create_process_memory, Get_process_name
    USE mo_jsb_memory_class, ONLY: t_jsb_memory

    CHARACTER(len=*), INTENT(in)           :: name
    CHARACTER(len=*), INTENT(in)           :: description
    INTEGER,          INTENT(in)           :: processes(:)
    CHARACTER(len=*), INTENT(in)           :: varlist_name
    CLASS(t_State),   POINTER,    OPTIONAL :: parent
    INTEGER,          INTENT(in), OPTIONAL :: model_id   ! Only needed if parent is not present
    INTEGER,          INTENT(in), OPTIONAL :: grid_id    ! Only needed if parent is not present
    TYPE(t_jsb_tile), POINTER              :: return_ptr

    CLASS(t_jsb_memory), POINTER :: tmp_ptr
    INTEGER :: i, iproc

    ! Allocate and initialize base structure (=state) for tile
    ALLOCATE(return_ptr)
    CALL return_ptr%Init(name, description, parent)

    ! Initialize further components of t_jsb_tile extended type

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

    ALLOCATE(return_ptr%mem(max_no_of_processes))
    DO i=1,SIZE(processes)
      iproc = processes(i)
      !return_ptr%param(iproc)%p => t_jsb_param(iproc)
      tmp_ptr => Create_process_memory(iproc)
      SELECT TYPE (tmp_ptr)
      CLASS IS (t_jsb_memory)
        return_ptr%mem(iproc)%p => tmp_ptr
        return_ptr%mem(iproc)%p%id = iproc
        return_ptr%mem(iproc)%p%varlist_name = TRIM(varlist_name)//'_'//TRIM(Get_process_name(iproc))
        return_ptr%mem(iproc)%p%grid_id = return_ptr%grid_id
        return_ptr%mem(iproc)%p%owner_model_id = return_ptr%owner_model_id
        IF (ASSOCIATED(return_ptr%parent_tile)) THEN
          CALL return_ptr%mem(iproc)%p%Init(TRIM(Get_process_name(iproc)), TRIM(name))
        ELSE
!!$          CALL return_ptr%mem(iproc)%p%Init(TRIM(Get_process_name(iproc)), '')
          CALL return_ptr%mem(iproc)%p%Init('', '')
        END IF
      END SELECT
      NULLIFY(tmp_ptr)
!!$      return_ptr%mem(iproc)%p => Create_process_memory(iproc)
    END DO

  END FUNCTION Construct_tile

  SUBROUTINE Init_tile(tile, varlist_name, prefix, suffix, grid_id)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, Get_varlist, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, nc_bits, t_cf, t_grib1, t_grib2, &
                                    FILETYPE_NC2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    TYPE(t_jsb_tile),   POINTER,  INTENT(inout) :: tile
    CHARACTER(len=*),             INTENT(in)    :: varlist_name
    CHARACTER(len=*),             INTENT(in)    :: prefix        ! Prefix for variable names
    CHARACTER(len=*),             INTENT(in)    :: suffix
    INTEGER,                      INTENT(in)    :: grid_id

    TYPE(t_jsb_varlist), POINTER :: this_list

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
      hgrid, surface,                                                                                         &
      t_cf('land_area_fraction', '',                                                                          &
               'fraction of '//TRIM(parent_tile_name)//' covered by '//TRIM(tile%name), nc_bits),             &
      t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      initval_r=0.0_wp)

    CALL jsb_add_var(TRIM(varlist_name), prefix//'fract_max'//suffix, tile%fract_max,                         &
      hgrid, surface,                                                                                         &
      t_cf('land_area_fraction', '',                                                                          &
               'fraction of '//TRIM(parent_tile_name)//' covered by '//TRIM(tile%name), nc_bits),             &
      t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      initval_r=0.0_wp)

    tile%fract = 1._wp
    tile%fract_max = 1._wp

  END SUBROUTINE Init_tile

  FUNCTION Get_parent_memory(this, iproc) RESULT(return_ptr)

    CLASS(t_jsb_tile),  INTENT(in) :: this
    INTEGER,            INTENT(in) :: iproc
    CLASS(t_jsb_memory), POINTER   :: return_ptr

    CLASS(*), POINTER :: parent

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

    return_value = ASSOCIATED(this%mem(iproc)%p)

  END FUNCTION Is_process_active

#endif
END MODULE mo_jsb_tile_class
