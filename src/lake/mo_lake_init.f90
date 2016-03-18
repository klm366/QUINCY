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

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_config_class,  ONLY: t_jsb_config_p
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_close_file, jsb_netcdf_open_input, &
    jsb_netcdf_read_2d, jsb_netcdf_read_2d_time, jsb_netcdf_read_2d_extdim
#ifndef __NO_JSBACH_HD__
  USE mo_impl_constants,    ONLY: LAND_BOUNDARY, SEA_BOUNDARY, SEA
#endif

  dsl4jsb_Use_processes LAKE
  dsl4jsb_Use_config(LAKE)
  dsl4jsb_Use_memory(LAKE)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: lake_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lake_init'

CONTAINS
  !
  !> Intialize lake process (after memory has been set up)
  !
  SUBROUTINE lake_init(tile, config)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p), INTENT(in)    :: config(:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':lake_init'

    CALL lake_init_bc(tile, config)
    CALL lake_init_ic(tile, config)

  END SUBROUTINE lake_init

  SUBROUTINE lake_init_bc(tile, config)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':lake_init_bc'

    model => get_model(tile%owner_model_id)

  END SUBROUTINE lake_init_bc

  RECURSIVE SUBROUTINE lake_init_ic(tile, config)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)

    TYPE(t_jsb_model), POINTER :: model

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
      & fill_array = dsl4jsb_var2D_ptr(LAKE,land_fract))

    ! Lake mask
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='lake',                     &
      & fill_array = tile%fract)
    ! set lake mask to zero on ocean grid cells (instead of missing value)
    WHERE (dsl4jsb_var2D_onDomain(LAKE,land_fract) == 0._wp)
      tile%fract = 0._wp
    END WHERE

    CALL jsb_netcdf_close_file(input_file)

  END SUBROUTINE lake_init_ic

#endif
END MODULE mo_lake_init
