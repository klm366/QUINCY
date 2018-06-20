!>
!! @brief mo_jsb_grid from JSBACH4
!! 
MODULE mo_jsb_grid
#ifndef __NO_JSBACH__

  USE mo_exception,      ONLY: message, finish
  USE mo_jsb_grid_class, ONLY: t_jsb_grid, t_jsb_grid_ptr
  USE mo_jsb_control,    ONLY: debug_on, max_no_of_grids => max_no_of_models

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_grid, Get_grid

  ! INTEGER, PARAMETER :: max_no_of_grids = 1

  TYPE(t_jsb_grid_ptr),  SAVE :: jsbach_grids(max_no_of_grids)

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

#endif
END MODULE mo_jsb_grid
