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

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract

  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_config_class,  ONLY: t_jsb_config_p
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_close_file, jsb_netcdf_open_input, &
    jsb_netcdf_read_2d, jsb_netcdf_read_2d_time
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
#ifndef __NO_JSBACH_HD__
  USE mo_impl_constants,    ONLY: LAND_BOUNDARY, SEA_BOUNDARY, SEA
#endif

  dsl4jsb_Use_processes SEB, HD
  dsl4jsb_Use_config(SEB)
  dsl4jsb_Use_config(HD)
  dsl4jsb_Use_memory(SEB)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: seb_init_ic, seb_init_bc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_init'

CONTAINS

  SUBROUTINE seb_init_bc(tile, config)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p), INTENT(in)    :: config(:)

    TYPE(t_jsb_model),       POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':seb_init_bc'

    model => get_model(tile%owner_model_id)

  END SUBROUTINE seb_init_bc

  RECURSIVE SUBROUTINE seb_init_ic(tile, config)

    USE mo_jsb_time, ONLY: get_time_interpolation_weights

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p), INTENT(in)    :: config(:)

    TYPE(t_jsb_model),        POINTER :: model
    dsl4jsb_Def_config(SEB)
    dsl4jsb_Def_config(HD)
    dsl4jsb_Def_memory(SEB)

    REAL(wp), POINTER :: t(:,:), t_mon(:,:,:), cell_sea_land_mask(:,:)
    REAL(wp), POINTER :: &
      & return_pointer (:,  :)       !< temporary pointer

    TYPE(t_input_file) :: input_file

    INTEGER  :: i
    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2, n1, n2

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':seb_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get seb config
    dsl4jsb_Get_config(SEB)
    dsl4jsb_Get_config(HD)

    IF (.NOT. dsl4jsb_Config(SEB)%active) RETURN

    ! Get seb memory of the tile
    dsl4jsb_Get_memory(SEB)

    ! Initialize parameters
    !
    CALL message(TRIM(routine), &
         'Initializing seb memory for tile '//TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(SEB)%bc_sso_filename))
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SEB)%bc_sso_filename), model%id)

    ! Orography
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='orostd',                   &
      & fill_array = dsl4jsb_var2D_ptr(SEB,oro_stddev))
    dsl4jsb_var2D_onDomain(SEB,oro_stddev) = &
      MERGE(dsl4jsb_var2D_onDomain(SEB,oro_stddev), 0._wp, dsl4jsb_var2D_onDomain(SEB,oro_stddev) >= 0._wp)

    CALL jsb_netcdf_close_file(input_file)

    CALL message(TRIM(routine), 'Initializing seb memory for tile '//TRIM(tile%name)//' from bc_land_frac.nc')
    input_file = jsb_netcdf_open_input('bc_land_frac.nc', model%id)     ! Temporary
!!$    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SEB)%bc_filename), model%id)

    ! Land mask (including glacies and lakes)
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='notsea',                   &        !< ocean 0., land 1.
      & fill_array = dsl4jsb_var2D_ptr(SEB,land))

    ! Glacier mask
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='glac',                     &
      & fill_array = dsl4jsb_var2D_ptr(SEB,glacier))
    ! set glacier mask to zero on ocean grid cells (instead of missing value)
    WHERE (dsl4jsb_var2D_onDomain(SEB,land) == 0._wp)
      dsl4jsb_var2D_onDomain(SEB,glacier) = 0._wp
    END WHERE

    ! Lake mask
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='lake',                     &
      & fill_array = dsl4jsb_var2D_ptr(SEB,lake))

    CALL jsb_netcdf_close_file(input_file)

#ifndef __NO_JSBACH_HD__
    IF (dsl4jsb_Config(HD)%active) THEN

      input_file = jsb_netcdf_open_input('r2b4_amip.nc', model%id)    !vg Name !!!

      ! coastal cells of the land and ocean
      cell_sea_land_mask => jsb_netcdf_read_2d( &
        & input_file,                                &
        & variable_name='cell_sea_land_mask')
      dsl4jsb_var2D_onDomain(SEB,coast_land)  = &
        MERGE(1._wp, 0._wp, cell_sea_land_mask(:,:) == LAND_BOUNDARY)
      dsl4jsb_var2D_onDomain(SEB,coast_ocean) = &
        MERGE(1._wp, 0._wp, cell_sea_land_mask(:,:) == SEA_BOUNDARY)
      dsl4jsb_var2D_onDomain(SEB,ocean)       = MERGE(1._wp, 0._wp, &
        cell_sea_land_mask(:,:) == SEA .OR. cell_sea_land_mask(:,:) == SEA_BOUNDARY)
      DEALLOCATE(cell_sea_land_mask)

      CALL jsb_netcdf_close_file(input_file)

    END IF
#endif

    ! Initialize energy state
    !

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SEB)%ic_filename), model%id)

    ! Initial surface temperature
!!$    return_pointer => jsb_netcdf_read_2d(         &
!!$      & filename=seb_conf%ic_filename,        &
!!$      & variable_name='surface_temperature',          &
!!$      fill_array=seb_mem%t_seb%ptr,                        &
!!$      & model_id=model%id)
!!$    t_seb => seb_mem%t_seb
    t_mon => jsb_netcdf_read_2d_time(         &   ! Temporary
      & input_file,                               &
      & variable_name='surf_temp')
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

    CALL jsb_netcdf_close_file(input_file)

    t => dsl4jsb_var2D_ptr(SEB,t)

    ! Initialize energy fluxes
    !
    dsl4jsb_var2D_onDomain(SEB,latent_hflx)     = 0._wp
    dsl4jsb_var2D_onDomain(SEB,sensible_hflx)   = 0._wp

  END SUBROUTINE seb_init_ic

#endif
END MODULE mo_seb_init
