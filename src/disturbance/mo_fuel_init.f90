!>
!! @brief Initialization of the fuel memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-06-20)
!!
MODULE mo_fuel_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,  ONLY: t_jsb_config_p
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_time_iface,    ONLY: is_time_experiment_start, is_time_restart

  dsl4jsb_Use_processes FUEL_
  dsl4jsb_Use_config(FUEL_)
  dsl4jsb_Use_memory(FUEL_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: fuel_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fuel_init'

CONTAINS

  !
  !> Intialize fuel process (after memory has been set up)
  !
  SUBROUTINE fuel_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':fuel_init'

    dsl4jsb_Def_memory(FUEL_)
    dsl4jsb_Get_memory(FUEL_)

    ! No initialization needed at this time
    RETURN

    CALL fuel_init_ic(tile)

    CALL fuel_init_bc(tile)



  END SUBROUTINE fuel_init

  SUBROUTINE fuel_init_bc(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    TYPE(t_jsb_model),       POINTER :: model
    dsl4jsb_Def_memory(FUEL_)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':fuel_init_bc'

    !-----------------------------------------------------------------------
    ! CONTENT

    model => get_model(tile%owner_model_id)

    ! Get fuel memory of the tile
    dsl4jsb_Get_memory(FUEL_)

    IF (.NOT. tile%Is_process_active(FUEL_)) RETURN

    CALL message(TRIM(routine), 'Setting fuel boundary conditions for tile '//TRIM(tile%name))

    ! R: nothing to do in the moment ...

  END SUBROUTINE fuel_init_bc

  RECURSIVE SUBROUTINE fuel_init_ic(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model

    dsl4jsb_Def_config(FUEL_)
    dsl4jsb_Def_memory(FUEL_)

    REAL(wp), POINTER :: &
      & return_pointer  (:,  :) !, & !< temporary pointer
      ! & return_pointer3d(:,:,:)

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':fuel_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get fuel memory and config  of the tile
    dsl4jsb_Get_memory(FUEL_)
    dsl4jsb_Get_config(FUEL_)

    IF (.NOT. tile%Is_process_active(FUEL_)) RETURN

    CALL message(TRIM(routine), 'Setting fuel initial conditions for tile '//TRIM(tile%name))

    ! Initialize physical state
    !
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(FUEL_)%ic_filename), model%grid_id)

    CALL input_file%Close()

  END SUBROUTINE fuel_init_ic

#endif
END MODULE mo_fuel_init
