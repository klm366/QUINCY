MODULE mo_jsb_grid
#ifndef __NO_JSBACH__

  USE mo_exception,      ONLY: message, finish
  USE mo_jsb_control,    ONLY: debug_on
  USE mo_jsb_grid_class, ONLY: t_jsb_grid, t_jsb_grid_ptr, t_jsb_vgrid, t_jsb_vgrid_ptr
  USE mo_jsb_control,    ONLY: max_no_of_grids => max_no_of_models

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_grid, Register_vgrid, Get_grid, Get_vgrid

  ! INTEGER, PARAMETER :: max_no_of_grids = 1
  INTEGER, PARAMETER :: max_no_of_vgrids = 10

  TYPE(t_jsb_grid_ptr),  SAVE :: jsbach_grids(max_no_of_grids)
  TYPE(t_jsb_vgrid_ptr), SAVE :: jsbach_vgrids(max_no_of_vgrids)

CONTAINS

  SUBROUTINE Register_grid(grid)

    TYPE(t_jsb_grid),  POINTER, INTENT(in)    :: grid

    INTEGER :: i, idx
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid:Register_grid'

    IF (debug_on()) CALL message(TRIM(routine), 'going to register new JSBACH grid instance')

    idx = 0

    DO i=1,max_no_of_grids
      IF (ASSOCIATED(jsbach_grids(i)%p)) THEN
        IF (jsbach_grids(i)%p%id == grid%id) EXIT
      END IF
    END DO
    IF (idx /= 0) RETURN        ! axis already exists and is already registered for model

    ! Find free slot for vgrid
    idx = 0
    DO i=1,max_no_of_grids
      IF (.NOT. ASSOCIATED(jsbach_grids(i)%p)) THEN
        idx = i
        EXIT
      END IF
    END DO
    IF(idx == 0) CALL finish(TRIM(routine), 'Max no of grids reached')

    jsbach_grids(idx)%p => grid

    CALL message(TRIM(routine), 'registered new JSBACH grid instance '//TRIM(grid%name))

  END SUBROUTINE Register_grid

  SUBROUTINE Register_vgrid(vgrid)

    TYPE(t_jsb_vgrid), POINTER, INTENT(in)    :: vgrid

    INTEGER :: i, idx
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid:Register_vgrid'

    idx = 0
    DO i=1,max_no_of_vgrids
      IF (ASSOCIATED(jsbach_vgrids(i)%p)) THEN
        IF (jsbach_vgrids(i)%p%cdiZaxisID == vgrid%cdiZaxisID) EXIT
      END IF
    END DO
    IF (idx /= 0) RETURN        ! axis already exists and is already registered for model

    ! Find free slot for vgrid
    idx = 0
    DO i=1,max_no_of_vgrids
      IF (.NOT. ASSOCIATED(jsbach_vgrids(i)%p)) THEN
        idx = i
        EXIT
      END IF
    END DO
    IF(idx == 0) CALL finish(TRIM(routine), 'Max no of vgrids reached')

    jsbach_vgrids(idx)%p => vgrid

!!$    self%vgrids(idx)%id         = vgrid%id
!!$    self%vgrids(idx)%name       = TRIM(vgrid%name)
!!$    self%vgrids(idx)%longname   = TRIM(vgrid%longname)
!!$    self%vgrids(idx)%units      = TRIM(vgrid%units)
!!$    self%vgrids(idx)%n_levels   = vgrid%n_levels
!!$    self%vgrids(idx)%cdiZaxisID = vgrid%cdiZaxisID
!!$    ALLOCATE(self%vgrids(idx)%levels (vgrid%n_levels))
!!$    ALLOCATE(self%vgrids(idx)%lbounds(vgrid%n_levels))
!!$    ALLOCATE(self%vgrids(idx)%ubounds(vgrid%n_levels))
!!$    self%vgrids(idx)%levels(:)  = vgrid%levels(:)
!!$    self%vgrids(idx)%lbounds(:) = vgrid%lbounds(:)
!!$    self%vgrids(idx)%ubounds(:) = vgrid%ubounds(:)

  END SUBROUTINE Register_vgrid

  FUNCTION Get_grid(id) RESULT(grid)

    INTEGER, INTENT(in)        :: id
    TYPE(t_jsb_grid), POINTER  :: grid

    INTEGER :: i, idx
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid:Get_grid'

    grid => NULL()
    DO i=1,SIZE(jsbach_grids)
      IF (ASSOCIATED(jsbach_grids(i)%p)) THEN
        IF (jsbach_grids(i)%p%id == id) THEN
          grid => jsbach_grids(i)%p
          EXIT
        END IF
      END IF
    END DO

    IF (.NOT. ASSOCIATED(grid)) THEN
      CALL finish(TRIM(routine), 'Grid not found.')
    END IF

  END FUNCTION Get_grid

  FUNCTION Get_vgrid(name) RESULT(vgrid)

    CHARACTER(len=*),           INTENT(in)    :: name
    TYPE(t_jsb_vgrid), POINTER                :: vgrid

    INTEGER :: i, idx
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid:Get_vgrid'

    vgrid => NULL()
    DO i=1,SIZE(jsbach_vgrids)
      IF (ASSOCIATED(jsbach_vgrids(i)%p)) THEN
        IF (jsbach_vgrids(i)%p%id > 0 .AND. TRIM(jsbach_vgrids(i)%p%name) == TRIM(name)) THEN
          vgrid => jsbach_vgrids(i)%p
          EXIT
        END IF
      END IF
    END DO

    IF (.NOT. ASSOCIATED(vgrid)) THEN
      CALL finish(TRIM(routine), 'Vgrid '//TRIM(name)//' not found.')
    END IF

  END FUNCTION Get_vgrid

#endif
END MODULE mo_jsb_grid
