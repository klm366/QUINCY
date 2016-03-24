!>
!! @brief Initialization of the the glacier memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_glacier_init
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

  dsl4jsb_Use_processes GLACIER
  dsl4jsb_Use_config(GLACIER)
  dsl4jsb_Use_memory(GLACIER)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: glacier_init_ic, glacier_init_bc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_glacier_init'

CONTAINS

  SUBROUTINE glacier_init_bc(tile, config)

    TYPE(t_jsb_tile), INTENT(inout) :: tile
    CLASS(t_jsb_config_p), INTENT(in)    :: config(:)

    TYPE(t_jsb_model),       POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':glacier_init_bc'

    model => get_model(tile%owner_model_id)

  END SUBROUTINE glacier_init_bc

  RECURSIVE SUBROUTINE glacier_init_ic(tile, config)

    TYPE(t_jsb_tile), INTENT(inout) :: tile
    CLASS(t_jsb_config_p), INTENT(in)    :: config(:)

    TYPE(t_jsb_model),        POINTER :: model
    dsl4jsb_Def_config(GLACIER)
    dsl4jsb_Def_memory(GLACIER)

    REAL(wp), POINTER :: return_pointer (:,:) !< temporary pointer

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':glacier_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get glacier memory of the tile
    dsl4jsb_Get_memory(GLACIER)

    ! Get glacier config
    dsl4jsb_Get_config(GLACIER)

    IF (.NOT. dsl4jsb_Config(GLACIER)%active) RETURN

    ! Initialize parameters
    !
    CALL message(TRIM(routine), 'Initializing glacier memory for tile '//TRIM(tile%name)//' from bc_land_frac.nc')
    input_file = jsb_netcdf_open_input('bc_land_frac.nc', model%id)     ! Temporary
!!$    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(GLACIER)%bc_filename), model%id)

    ! Land mask (including glacies and lakes)
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='notsea',                   &        !< ocean 0., land 1.
      & fill_array = dsl4jsb_var2D_ptr(GLACIER,land))

    ! Glacier mask
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='glac',                     &
      & fill_array = dsl4jsb_var2D_ptr(GLACIER,glacier))
    ! set glacier mask to zero on ocean grid cells (instead of missing value)
    WHERE (dsl4jsb_var2D_onDomain(GLACIER,land) == 0._wp)
      dsl4jsb_var2D_onDomain(GLACIER,glacier) = 0._wp
    END WHERE

    CALL jsb_netcdf_close_file(input_file)

    ! Propagate initial conditions down the tile structure, starting with the tile with which this subroutine was called initially
    ! TBD: sub-tiles should not read the netcdf files again
!!$    IF (tile%no_of_sub_tiles > 0) THEN
!!$      DO i=1,tile%no_of_sub_tiles
!!$        CALL glac_init_ic(tile%sub_tiles(i))
!!$      END DO
!!$    END IF

  END SUBROUTINE glacier_init_ic

#endif
END MODULE mo_glacier_init
