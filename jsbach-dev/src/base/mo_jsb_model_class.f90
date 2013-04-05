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

  USE mo_jsb_config_class, ONLY: t_jsb_model_config
  USE mo_jsb_grid_class,   ONLY: t_jsb_grid
  USE mo_jsb_tile_class,   ONLY: t_jsb_tile

  USE mo_exception,        ONLY: message_text, message, finish
  USE mo_io_units,         ONLY: filename_max

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_model_m, t_jsb_model
  PUBLIC :: new_model, print_model, delete_model
  PUBLIC :: register_grid

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
    TYPE(t_jsb_grid),         POINTER :: grid              !< Horizontal grid
    TYPE(t_jsb_tile),         POINTER :: land              !< Land structure and model state
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
  FUNCTION new_model(name, shortname, description, namelist_filename) RESULT(model)

    CHARACTER(len=*),           INTENT(in) :: name
    CHARACTER(len=*),           INTENT(in) :: shortname
    CHARACTER(len=*),           INTENT(in) :: description
    CHARACTER(len=*),           INTENT(in) :: namelist_filename !< Model namelist file
    TYPE(t_jsb_model), POINTER             :: model

    TYPE(t_jsb_model), POINTER :: new

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:new_model'

    CALL message(TRIM(routine), 'starting construction of new JSBACH model instance.')

    ALLOCATE(new)

    new%name              = TRIM(name)
    new%shortname         = TRIM(shortname)
    new%description       = TRIM(description)
    new%namelist_filename = TRIM(namelist_filename)
    new%config => NULL()
    new%grid => NULL()

    CALL message(TRIM(routine), 'construction of JSBACH model completed.')

    model => new

  END FUNCTION new_model

  !>
  !! Destructor method for model instance.
  !!
  !! Subroutine to delete model instance from memory.
  !!
  SUBROUTINE delete_model(self)

    USE mo_jsb_grid_class, ONLY: delete_grid

#ifdef HAVE_F2003
    TYPE(t_jsb_model), POINTER, INTENT(inout) :: self
#else
    TYPE(t_jsb_model), POINTER                :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:delete_model'

    CALL message(TRIM(routine), 'starting destruction of JSBACH model instance')

    IF (ASSOCIATED(self%grid)) CALL delete_grid(self%grid)
    DEALLOCATE(self)

    CALL message(TRIM(routine), 'destruction of JSBACH model instance completed.')

  END SUBROUTINE delete_model

  !>
  !! Print method for model instance.
  !!
  !! Subroutine to print description of model instance.
  !!
  SUBROUTINE print_model(self)

    USE mo_jsb_grid_class, ONLY: print_grid
    USE mo_jsb_control,    ONLY: model_base_dir

    TYPE(t_jsb_model), INTENT(in) :: self

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:print_model'

    CALL message(TRIM(routine), 'print JSBACH model description')

    WRITE(message_text,'(A,I2,A,A,A,A,A)') 'Model id: ', self%id, ', model name: ', TRIM(self%name), &
                                     ' (', TRIM(self%shortname), ')'
    CALL message('', message_text)
    CALL message('', '... '//TRIM(self%description))

    CALL message('', 'Namelist file: '//TRIM(model_base_dir)//'/'//TRIM(self%namelist_filename))

    IF (ASSOCIATED(self%grid)) CALL print_grid(self%grid)

  END SUBROUTINE print_model

  SUBROUTINE register_grid(self, grid)

    USE mo_jsb_grid_class, ONLY: t_jsb_grid

#ifdef HAVE_F2003
    TYPE(t_jsb_model), POINTER, INTENT(inout) :: self
#else
    TYPE(t_jsb_model), POINTER                :: self
#endif
    TYPE(t_jsb_grid),  POINTER, INTENT(in) :: grid

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:register_grid'

    CALL message(TRIM(routine), 'going to register new JSBACH grid instance')

    NULLIFY(self%grid)
    self%grid => grid

    CALL message(TRIM(routine), 'registered new JSBACH grid instance '//TRIM(grid%name)//&
                                ' for model '//TRIM(self%name))

  END SUBROUTINE register_grid

END MODULE mo_jsb_model_class
