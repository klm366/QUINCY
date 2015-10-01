!>
!! @brief Contains main JSBACH model class.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_model_class
#ifndef __NO_JSBACH__

  USE mo_hsm_class,        ONLY: t_hsm, t_state
  USE mo_jsb_config_class, ONLY: t_jsb_model_config
  USE mo_jsb_process_class, ONLY: t_jsb_process

  USE mo_exception,        ONLY: message_text, message, finish
  USE mo_io_units,         ONLY: filename_max

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_model_m, t_jsb_model
  PUBLIC :: new_model, print_model, delete_model

  !>
  !! Type for model structure
  !!
  !! This derived type defines one instance of the JSBACH model
  !!
  TYPE t_jsb_model
    INTEGER                           :: id                !< Model id
    CHARACTER(len=30)                 :: name              !< Model name
    CHARACTER(len=10)                 :: shortname         !< Model short name
    CHARACTER(len=132)                :: description       !< Model description
    CHARACTER(len=filename_max)       :: namelist_filename !< Model namelist file
    TYPE(t_jsb_model_config), POINTER :: config            !< Configuration
    INTEGER                           :: grid_id
    TYPE(t_jsb_process), ALLOCATABLE  :: processes(:)
    TYPE(t_hsm),              POINTER :: hsm
  CONTAINS
    PROCEDURE :: Configure_processes => Configure_model_processes
    PROCEDURE :: Get_shortname       => Get_model_shortname
    PROCEDURE :: Get_top_tile        => Get_model_top_tile
  END TYPE t_jsb_model

  !>
  !! Container CLASS for JSBACH model
  !!
  !! Contains pointer to an instance of type t_jsb_model in order to allow arrays
  !! of pointers to models.
  !!
  TYPE t_jsb_model_m
    TYPE(t_jsb_model), POINTER :: m
  END TYPE t_jsb_model_m

CONTAINS

  !>
  !! Constructor method for model instance.
  !!
  !! Subroutine to construct new instance of basic model structure.
  !!
  FUNCTION new_model(name, shortname, description, namelist_filename) RESULT(return_ptr)

    CHARACTER(len=*),           INTENT(in) :: name
    CHARACTER(len=*),           INTENT(in) :: shortname
    CHARACTER(len=*),           INTENT(in) :: description
    CHARACTER(len=*),           INTENT(in) :: namelist_filename !< Model namelist file
    TYPE(t_jsb_model), POINTER             :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:new_model'

    CALL message(TRIM(routine), 'starting construction of new JSBACH model instance.')

    ALLOCATE(return_ptr)

    return_ptr%name              = TRIM(name)
    return_ptr%shortname         = TRIM(shortname)
    return_ptr%description       = TRIM(description)
    return_ptr%namelist_filename = TRIM(namelist_filename)
    return_ptr%config => NULL()

    CALL message(TRIM(routine), 'construction of JSBACH model completed.')

  END FUNCTION new_model

  SUBROUTINE Configure_model_processes(this, processes)

    CLASS(t_jsb_model), INTENT(inout) :: this
    INTEGER,            INTENT(in)    :: processes(:)

    INTEGER :: i, iproc

    DO i=1,SIZE(processes)
      iproc = processes(i)
      CALL this%processes(iproc)%Configure(iproc, this%namelist_filename)
    END DO

  END SUBROUTINE Configure_model_processes

  FUNCTION Get_model_name(this) RESULT(return_value)

    CLASS(t_jsb_model), INTENT(in) :: this
    CHARACTER(LEN=:), ALLOCATABLE :: return_value

    return_value = this%name

  END FUNCTION Get_model_name

  FUNCTION Get_model_shortname(this) RESULT(return_value)

    CLASS(t_jsb_model), INTENT(in) :: this
    CHARACTER(LEN=:), ALLOCATABLE :: return_value

    return_value = this%shortname

  END FUNCTION Get_model_shortname

  FUNCTION Get_model_top_tile(this) RESULT(return_ptr)

    USE mo_jsb_tile_class,   ONLY: t_jsb_tile

    CLASS(t_jsb_model), INTENT(in) :: this
    TYPE(t_jsb_tile), POINTER :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Get_model_top_tile'

    return_ptr => NULL()

    SELECT TYPE (top => this%Hsm%top)
    TYPE IS (t_jsb_tile)
      return_ptr => top
    CLASS IS (t_State)
      CALL finish(TRIM(routine), 'tile is of type t_State, should be t_jsb_tile')
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

  END FUNCTION Get_model_top_tile

  !>
  !! Destructor method for model instance.
  !!
  !! Subroutine to delete model instance from memory.
  !!
  SUBROUTINE delete_model(self)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: self

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:delete_model'

    CALL message(TRIM(routine), 'starting destruction of JSBACH model instance')

    DEALLOCATE(self)

    CALL message(TRIM(routine), 'destruction of JSBACH model instance completed.')

  END SUBROUTINE delete_model

  !>
  !! Print method for model instance.
  !!
  !! Subroutine to print description of model instance.
  !!
  SUBROUTINE print_model(self)

    USE mo_jsb_control,    ONLY: model_base_dir
    USE mo_jsb_grid_class, ONLY: t_jsb_grid
    USE mo_jsb_grid,       ONLY: Get_grid

    TYPE(t_jsb_model), INTENT(in) :: self

    TYPE(t_jsb_grid), POINTER :: grid
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:print_model'

    CALL message(TRIM(routine), 'print JSBACH model description')

    WRITE(message_text,'(A,I2,A,A,A,A,A)') 'Model id: ', self%id, ', model name: ', TRIM(self%name), &
                                     ' (', TRIM(self%shortname), ')'
    CALL message('', message_text)
    CALL message('', '... '//TRIM(self%description))

    CALL message('', 'Namelist file: '//TRIM(model_base_dir)//'/'//TRIM(self%namelist_filename))

    grid => Get_grid(self%grid_id)
    IF (ASSOCIATED(grid)) CALL grid%Print

  END SUBROUTINE print_model

#endif
END MODULE mo_jsb_model_class
