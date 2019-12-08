!
! This file was automatically generated from:
!    /users/rschnur/icon-les-dev-staging.pgi_gpu/externals/jsbach/src/disturbance/mo_fuel_init.f90
!    with md5 hash: 7c37bc43d869b4aaa7821e20ba3a5768
!
!>
!! @brief Initialization of the fuel memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-06-20)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_fuel_init
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input

  USE mo_jsb_process_class, ONLY: fuel_
  USE mo_fuel_config_class, ONLY: t_fuel_config
  !USE mo_fuel_memory_class, ONLY: t_fuel_memory

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

    CHARACTER(len=*), PARAMETER :: routine = modname//':fuel_init'

    !TYPE(t_fuel_memory), POINTER :: fuel__mem
    !SELECT TYPE (mm => tile%mem(fuel_)%p); TYPE IS (t_fuel_memory); fuel__mem => mm ; END SELECT

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
    !TYPE(t_jsb_model),       POINTER :: model
    !TYPE(t_fuel_memory), POINTER :: fuel__mem

    CHARACTER(len=*), PARAMETER :: routine = modname//':fuel_init_bc'

    !-----------------------------------------------------------------------
    ! CONTENT

    !model => get_model(tile%owner_model_id)

    ! Get fuel memory of the tile
    !SELECT TYPE (mm => tile%mem(fuel_)%p); TYPE IS (t_fuel_memory); fuel__mem => mm ; END SELECT

    IF (.NOT. tile%Is_process_active(FUEL_)) RETURN

    CALL message(TRIM(routine), 'Setting fuel boundary conditions for tile '//TRIM(tile%name))

    ! R: nothing to do in the moment ...

  END SUBROUTINE fuel_init_bc

  RECURSIVE SUBROUTINE fuel_init_ic(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model

    TYPE(t_fuel_config), POINTER :: fuel__conf
    !TYPE(t_fuel_memory), POINTER :: fuel__mem

    !REAL(wp), POINTER :: &
    !  & return_pointer  (:,  :) !, & !< temporary pointer
    !  & return_pointer3d(:,:,:)

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = modname//':fuel_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get fuel memory and config  of the tile
    !SELECT TYPE (mm => tile%mem(fuel_)%p); TYPE IS (t_fuel_memory); fuel__mem => mm ; END SELECT
    SELECT TYPE (cc => model%processes(fuel_)%p%config); TYPE IS (t_fuel_config) ; fuel__conf => cc ; END SELECT

    IF (.NOT. tile%Is_process_active(FUEL_)) RETURN

    CALL message(TRIM(routine), 'Setting fuel initial conditions for tile '//TRIM(tile%name))

    ! Initialize physical state
    !
    input_file = jsb_netcdf_open_input(TRIM(fuel__conf%ic_filename), model%grid_id)

    CALL input_file%Close()

  END SUBROUTINE fuel_init_ic

#endif
END MODULE mo_fuel_init
