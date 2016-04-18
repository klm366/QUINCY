!>
!! @brief Initialization of the the seb memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_seb_init
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input

  dsl4jsb_Use_processes SEB
  dsl4jsb_Use_config(SEB)
  dsl4jsb_Use_memory(SEB)
  ! Hack until createMakefile handles cpp includes (TBD)
  USE mo_seb_memory_class, ONLY: t_seb_memory

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: seb_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_init'

CONTAINS

  !
  !> Intialize soil process (after memory has been set up)
  !
  SUBROUTINE seb_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':seb_init'

    CALL seb_init_bc(tile)
    CALL seb_init_ic(tile)

  END SUBROUTINE seb_init

  SUBROUTINE seb_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':seb_init_bc'

    model => get_model(tile%owner_model_id)

  END SUBROUTINE seb_init_bc

  RECURSIVE SUBROUTINE seb_init_ic(tile)

    USE mo_jsb_time, ONLY: get_time_interpolation_weights

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SEB)
    dsl4jsb_Def_memory(SEB)

    REAL(wp), POINTER :: t(:,:), t_mon(:,:,:)
    REAL(wp), POINTER :: &
      & return_pointer (:,  :)       !< temporary pointer

    TYPE(t_input_file) :: input_file

    INTEGER  :: i
    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2, n1, n2

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':seb_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get seb config
    dsl4jsb_Get_config(SEB)

    IF (.NOT. dsl4jsb_Config(SEB)%active) RETURN

    ! Get seb memory of the tile
    dsl4jsb_Get_memory(SEB)

    ! Initialize parameters
    !
    CALL message(TRIM(routine), 'Setting initial state of seb memory for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(SEB)%ic_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SEB)%ic_filename), model%grid_id)

    ! Initial surface temperature
!!$    return_pointer => input_file%Read_2d(         &
!!$      & filename=seb_conf%ic_filename,        &
!!$      & variable_name='surface_temperature',          &
!!$      fill_array=seb__mem%t%ptr,                        &
!!$      & model_id=model%id)
    t_mon => input_file%Read_2d_time(         &   ! Temporary
      &        variable_name='surf_temp')
    t_mon(:,:,:) = MERGE(t_mon(:,:,:), 280._wp, t_mon(:,:,:) > 0._wp)
    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
    IF (nmw1 == 0) THEN
      n1 = 12
    ELSE IF (nmw1 == 13) THEN
      n1 = 1
    ELSE
      n1 = nmw1
    END IF
    IF (nmw2 == 0) THEN
      n2 = 12
    ELSE IF (nmw2 == 13) THEN
      n2 = 1
    ELSE
      n2 = nmw2
    END IF
    dsl4jsb_var2D_onDomain(SEB,t) = wgt1 * t_mon(:,:,n1) + wgt2 * t_mon(:,:,n2)

    DEALLOCATE(t_mon)

    CALL input_file%Close()

  END SUBROUTINE seb_init_ic

#endif
END MODULE mo_seb_init
