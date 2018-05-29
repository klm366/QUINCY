!>
!! @brief Initialization of the fm memory.
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2018-03-19)
!!
MODULE mo_fm_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,  ONLY: t_jsb_config_p
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_time_iface,    ONLY: is_time_experiment_start, is_time_restart

  dsl4jsb_Use_processes FM_
  dsl4jsb_Use_config(FM_)
  dsl4jsb_Use_memory(FM_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: fm_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fm_init'

CONTAINS

  !
  !> Intialize fm process (after memory has been set up)
  !
  SUBROUTINE fm_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':fm_init'

    dsl4jsb_Def_memory(FM_)
    dsl4jsb_Get_memory(FM_)

    ! No initialization needed at this time -- dummy 
    RETURN

    CALL fm_init_ic(tile)

    CALL fm_init_bc(tile)



  END SUBROUTINE fm_init

  SUBROUTINE fm_init_bc(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    TYPE(t_jsb_model),       POINTER :: model
    dsl4jsb_Def_memory(FM_)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':fm_init_bc'

    !-----------------------------------------------------------------------
    ! CONTENT

    model => get_model(tile%owner_model_id)

    ! Get fm memory of the tile
    dsl4jsb_Get_memory(FM_)

    IF (.NOT. tile%Is_process_active(FM_)) RETURN

    CALL message(TRIM(routine), 'Setting fm boundary conditions for tile '//TRIM(tile%name))

    ! JN: nothing to do in the moment ...

  END SUBROUTINE fm_init_bc

  RECURSIVE SUBROUTINE fm_init_ic(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model

    dsl4jsb_Def_config(FM_)
    dsl4jsb_Def_memory(FM_)

    REAL(wp), POINTER :: &
      & return_pointer  (:,  :) !, & !< temporary pointer
      ! & return_pointer3d(:,:,:)

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':fm_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get fm memory and config  of the tile
    dsl4jsb_Get_memory(FM_)
    dsl4jsb_Get_config(FM_)

    IF (.NOT. tile%Is_process_active(FM_)) RETURN

    CALL message(TRIM(routine), 'Setting fm initial conditions for tile '//TRIM(tile%name))

  END SUBROUTINE fm_init_ic

#endif
END MODULE mo_fm_init
