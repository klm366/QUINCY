!>
!! @brief Initialization of the the jedi memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-04-12)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jedi_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_time,          ONLY: get_time_dt
  USE mo_jsb_control,       ONLY: debug_on

  dsl4jsb_Use_processes JEDI_
  dsl4jsb_Use_config(JEDI_)
  dsl4jsb_Use_memory(JEDI_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jedi_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jedi_init'

CONTAINS

  !
  !> Intialize jedi process (after memory has been set up)
  !
  SUBROUTINE jedi_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! TYPE(t_jsb_model), POINTER                :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':jedi_init'

    ! dsl4jsb_Def_memory(JEDI_)
    ! dsl4jsb_Def_config(JEDI_)

    ! model => get_model(tile%owner_model_id)

    ! dsl4jsb_Get_memory(JEDI_)
    ! dsl4jsb_Get_config(JEDI_)

    ! CALL jedi_init_ic(tile)
    ! CALL jedi_init_bc(tile)

  END SUBROUTINE jedi_init

  SUBROUTINE jedi_init_bc(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    ! TYPE(t_jsb_model),       POINTER :: model
    ! dsl4jsb_Def_memory(JEDI_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':jedi_init_bc'

    !-----------------------------------------------------------------------
    ! CONTENT

    ! model => get_model(tile%owner_model_id)

    ! Get jedi memory of the tile
    ! dsl4jsb_Get_memory(JEDI_)

    IF (.NOT. tile%Is_process_active(JEDI_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting  jedi boundary conditions for tile '//TRIM(tile%name))

    ! R: nothing to do in the moment ...

  END SUBROUTINE jedi_init_bc

  SUBROUTINE jedi_init_ic(tile)

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    ! TYPE(t_jsb_model),         POINTER :: model
    ! dsl4jsb_Def_config(JEDI_)
    ! dsl4jsb_Def_memory(JEDI_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':jedi_init_ic'

    !-----------------------------------------------------------------------
    ! CONTENT

    IF (.NOT. tile%Is_process_active(JEDI_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing jedi memory for tile '//TRIM(tile%name))

    ! model => Get_model(tile%owner_model_id)

    ! Get jedi memory of the tile
    ! dsl4jsb_Get_memory(JEDI_)

    ! R: nothing to do in the moment ...

  END SUBROUTINE jedi_init_ic

#endif
END MODULE mo_jedi_init
