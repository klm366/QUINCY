!>
!! @brief Initialization of the harv memory.
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2018-03-19)
!!
MODULE mo_harv_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,  ONLY: t_jsb_config_p
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_time_iface,    ONLY: is_time_experiment_start, is_time_restart

  dsl4jsb_Use_processes HARV_
  dsl4jsb_Use_config(HARV_)
  dsl4jsb_Use_memory(HARV_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: harv_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_harv_init'

CONTAINS

  !
  !> Intialize harv process (after memory has been set up)
  !
  SUBROUTINE harv_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':harv_init'

    dsl4jsb_Def_memory(HARV_)
    dsl4jsb_Get_memory(HARV_)

    ! No initialization needed at this time
    RETURN

    CALL harv_init_ic(tile)

    CALL harv_init_bc(tile)



  END SUBROUTINE harv_init

  SUBROUTINE harv_init_bc(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    TYPE(t_jsb_model),       POINTER :: model
    dsl4jsb_Def_memory(HARV_)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':harv_init_bc'

    !-----------------------------------------------------------------------
    ! CONTENT

    model => get_model(tile%owner_model_id)

    ! Get harv memory of the tile
    dsl4jsb_Get_memory(HARV_)

    IF (.NOT. tile%Is_process_active(HARV_)) RETURN

    CALL message(TRIM(routine), 'Setting harv boundary conditions for tile '//TRIM(tile%name))

    ! JN: nothing to do in the moment ...

  END SUBROUTINE harv_init_bc

  RECURSIVE SUBROUTINE harv_init_ic(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model

    dsl4jsb_Def_config(HARV_)
    dsl4jsb_Def_memory(HARV_)

    REAL(wp), POINTER :: &
      & return_pointer  (:,  :) !, & !< temporary pointer
      ! & return_pointer3d(:,:,:)

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':harv_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get harv memory and config  of the tile
    dsl4jsb_Get_memory(HARV_)
    dsl4jsb_Get_config(HARV_)

    IF (.NOT. tile%Is_process_active(HARV_)) RETURN

    CALL message(TRIM(routine), 'Setting harv initial conditions for tile '//TRIM(tile%name))

    ! Initialize physical state
    !
    ! JN: currently nothing to do   
!    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HARV_)%ic_filename), model%grid_id)

!    CALL input_file%Close()

  END SUBROUTINE harv_init_ic

#endif
END MODULE mo_harv_init
