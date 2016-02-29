!>
!! @brief Initialization of the the lake memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_lake_init
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile,          ONLY: t_jsb_tile

  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_config_class,  ONLY: t_jsb_config_p
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_close_file, jsb_netcdf_open_input, &
    jsb_netcdf_read_2d, jsb_netcdf_read_2d_time
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
#ifndef __NO_JSBACH_HD__
  USE mo_impl_constants,    ONLY: LAND_BOUNDARY, SEA_BOUNDARY, SEA
#endif

  dsl4jsb_Use_processes LAKE
  dsl4jsb_Use_config(LAKE)
  dsl4jsb_Use_memory(LAKE)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: lake_init_ic, lake_init_bc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lake_init'

CONTAINS

  SUBROUTINE lake_init_bc(tile, config)

    TYPE(t_jsb_tile), INTENT(inout) :: tile
    CLASS(t_jsb_config_p), INTENT(in)    :: config(:)

    TYPE(t_jsb_model),       POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':lake_init_bc'

    model => get_model(tile%owner_model_id)

  END SUBROUTINE lake_init_bc

  RECURSIVE SUBROUTINE lake_init_ic(tile, config)

    TYPE(t_jsb_tile), INTENT(inout) :: tile
    CLASS(t_jsb_config_p), INTENT(in)    :: config(:)

    TYPE(t_jsb_model),        POINTER :: model
    dsl4jsb_Def_config(LAKE)
    dsl4jsb_Def_memory(LAKE)

    REAL(wp), POINTER :: return_pointer (:,:) !< temporary pointer

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':lake_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get lake memory of the tile
    dsl4jsb_Get_memory(LAKE)

    ! Get lake config
    dsl4jsb_Get_config(LAKE)

    IF (.NOT. dsl4jsb_Config(LAKE)%active) RETURN

    ! Initialize parameters
    !
    CALL message(TRIM(routine), 'Initializing lake memory for tile '//TRIM(tile%name)//' from bc_land_frac.nc')
    input_file = jsb_netcdf_open_input('bc_land_frac.nc', model%id)     ! Temporary
!!$    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(LAKE)%bc_filename), model%id)

    ! Land mask (including glacies and lakes)
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='notsea',                   &        !< ocean 0., land 1.
      & fill_array = dsl4jsb_var2D_ptr(LAKE,land))

    ! Lake mask
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='lake',                     &
      & fill_array = dsl4jsb_var2D_ptr(LAKE,lake))
    ! set lake mask to zero on ocean grid cells (instead of missing value)
    WHERE (dsl4jsb_var2D_onDomain(LAKE,land) == 0._wp)
      dsl4jsb_var2D_onDomain(LAKE,lake) = 0._wp
    END WHERE

    CALL jsb_netcdf_close_file(input_file)

    ! Propagate initial conditions down the tile structure, starting with the tile with which this subroutine was called initially
    ! TBD: sub-tiles should not read the netcdf files again
!!$    IF (tile%no_of_sub_tiles > 0) THEN
!!$      DO i=1,tile%no_of_sub_tiles
!!$        CALL lake_init_ic(tile%sub_tiles(i))
!!$      END DO
!!$    END IF

  END SUBROUTINE lake_init_ic

#endif
END MODULE mo_lake_init
