!>
!! @brief Initialization of the the srf memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_srf_init
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile

  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_close_file, jsb_netcdf_open_input, &
    jsb_netcdf_read_2d, jsb_netcdf_read_2d_time
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
#ifndef __NO_JSBACH_HD__
  USE mo_impl_constants,    ONLY: LAND_BOUNDARY, SEA_BOUNDARY, SEA
#endif

  Use_processes SRF
  Use_config(srf)
  Use_memory(srf)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: srf_init_ic, srf_init_bc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_srf_init'

CONTAINS

  SUBROUTINE srf_init_bc(tile)

    TYPE(t_jsb_tile), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':srf_init_bc'

    model => get_model(tile%owner_model_id)

  END SUBROUTINE srf_init_bc

  RECURSIVE SUBROUTINE srf_init_ic(tile)

    USE mo_jsb_time, ONLY: get_time_interpolation_weights

    TYPE(t_jsb_tile), INTENT(inout) :: tile

    TYPE(t_jsb_model),        POINTER :: model
    Def_config(srf_conf, srf)
    Def_memory(srf_mem, srf)

    REAL(wp), POINTER :: t_srf(:,:), t_srf_mon(:,:,:), cell_sea_land_mask(:,:)
    REAL(wp), POINTER :: &
      & return_pointer (:,  :)       !< temporary pointer

    TYPE(t_input_file) :: input_file

    INTEGER  :: i
    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2, n1, n2

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':srf_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get srf memory of the tile
    Get_memory(srf_mem, SRF)

    ! Get srf config
    Get_config(srf_conf, SRF)

    IF (.NOT. srf_conf%active) RETURN

    ! Initialize parameters
    !
    CALL message(TRIM(routine), &
         'Initializing srf memory for tile '//TRIM(tile%name)//' from '//TRIM(srf_conf%bc_sso_filename))
    input_file = jsb_netcdf_open_input(TRIM(srf_conf%bc_sso_filename), model%id)

    ! Orography
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='orostd',                   &
      & fill_array = srf_mem%oro_stddev%ptr)
    srf_mem%oro_stddev%ptr(:,:) = MERGE(srf_mem%oro_stddev%ptr(:,:), 0._wp, srf_mem%oro_stddev%ptr(:,:) >= 0._wp)

    CALL jsb_netcdf_close_file(input_file)

    CALL message(TRIM(routine), 'Initializing srf memory for tile '//TRIM(tile%name)//' from bc_land_frac.nc')
    input_file = jsb_netcdf_open_input('bc_land_frac.nc', model%id)     ! Temporary
!!$    input_file = jsb_netcdf_open_input(TRIM(srf_conf%bc_filename), model%id)

    ! Land mask (including glacies and lakes)
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='notsea',                   &        !< ocean 0., land 1.
      & fill_array = srf_mem%land%ptr)

    ! Glacier mask
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='glac',                     &
      & fill_array = srf_mem%glacier%ptr)
    ! set glacier mask to zero on ocean grid cells (instead of missing value)
    WHERE (srf_mem%land%ptr(:,:) == 0._wp)
      srf_mem%glacier%ptr(:,:) = 0._wp
    END WHERE

    ! Lake mask
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='lake',                     &
      & fill_array = srf_mem%lake%ptr)

    CALL jsb_netcdf_close_file(input_file)

#ifndef __NO_JSBACH_HD__

    input_file = jsb_netcdf_open_input('r2b4_amip.nc', model%id)    !vg Name !!!

    ! coastal cells of the land and ocean
    cell_sea_land_mask => jsb_netcdf_read_2d( &
      & input_file,                                &
      & variable_name='cell_sea_land_mask')
    srf_mem%coast_land%ptr(:,:)  = MERGE(1._wp, 0._wp, cell_sea_land_mask(:,:) == LAND_BOUNDARY)
    srf_mem%coast_ocean%ptr(:,:) = MERGE(1._wp, 0._wp, cell_sea_land_mask(:,:) == SEA_BOUNDARY)
    srf_mem%ocean%ptr(:,:)       = MERGE(1._wp, 0._wp, &
      cell_sea_land_mask(:,:) == SEA .OR. cell_sea_land_mask(:,:) == SEA_BOUNDARY)
    DEALLOCATE(cell_sea_land_mask)

    CALL jsb_netcdf_close_file(input_file)

#endif

    input_file = jsb_netcdf_open_input(TRIM(srf_conf%bc_filename), model%id)

    ! Surface roughness length
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='roughness_length',         &
      & fill_array = srf_mem%rough_m_srf%ptr)
    srf_mem%rough_m_srf%ptr(:,:) = MERGE(srf_mem%rough_m_srf%ptr(:,:), 1._wp, srf_mem%rough_m_srf%ptr(:,:) >= 0._wp)

    ! Albedo
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='albedo',                   &
      & fill_array = srf_mem%albedo_srf%ptr)
    srf_mem%albedo_srf%ptr(:,:) = MERGE(srf_mem%albedo_srf%ptr(:,:), 0._wp, srf_mem%albedo_srf%ptr(:,:) >= 0._wp)

    srf_mem%alb_vis_srf%ptr(:,:) = srf_mem%albedo_srf%ptr(:,:)
    srf_mem%alb_nir_srf%ptr(:,:) = srf_mem%albedo_srf%ptr(:,:)

    CALL jsb_netcdf_close_file(input_file)

    ! Initialize energy state
    !

    input_file = jsb_netcdf_open_input(TRIM(srf_conf%ic_filename), model%id)

    ! Initial surface temperature
!!$    return_pointer => jsb_netcdf_read_2d(         &
!!$      & filename=srf_conf%ic_filename,        &
!!$      & variable_name='surface_temperature',          &
!!$      fill_array=srf_mem%t_srf%ptr,                        &
!!$      & model_id=model%id)
!!$    t_srf => srf_mem%t_srf
    t_srf_mon => jsb_netcdf_read_2d_time(         &   ! Temporary
      & input_file,                               &
      & variable_name='surf_temp')
    t_srf_mon(:,:,:) = MERGE(t_srf_mon(:,:,:), 280._wp, t_srf_mon(:,:,:) > 0._wp)
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
    srf_mem%t_srf%ptr(:,:) = wgt1 * t_srf_mon(:,:,n1) + wgt2 * t_srf_mon(:,:,n2)

    DEALLOCATE(t_srf_mon)

    CALL jsb_netcdf_close_file(input_file)

    t_srf => srf_mem%t_srf%ptr
    ! Note: t_srf now points to srf_mem%t_srf%ptr and can be
    ! used below to fill %t_srf for all other files:
    !    tile%mem%srf%energy_state%t_srf%ptr(:,:) = t_srf(:,:)

    srf_mem%t_srf_old%ptr        (:,:) = t_srf (:,:)
    srf_mem%t_srf_unfilt%ptr     (:,:) = t_srf (:,:)
    srf_mem%t_srf_unfilt_old%ptr (:,:) = t_srf (:,:)
    srf_mem%t_eff_srf%ptr        (:,:) = t_srf (:,:)
!!$    srf_mem%t_eff_srf%ptr        (:,:) = t_srf (:,:)
    srf_mem%s_srf%ptr            (:,:) = 0._wp
    srf_mem%qsat_srf%ptr         (:,:) = -1._wp  ! Negative values will be used to test for start of experiment
!!$    srf_mem%qsat_srf_old%ptr     (:,:) = 0._wp
    srf_mem%tte_corr%ptr         (:,:) = 0._wp
    srf_mem%fact_q_air%ptr       (:,:) = 0.5_wp
    srf_mem%fact_qsat_srf%ptr    (:,:) = 0.5_wp
    srf_mem%fact_qsat_trans_srf%ptr (:,:) = 0.5_wp

    ! Initialize water state

    srf_mem%wfract_srf%ptr  (:,:) = 0._wp
    srf_mem%wsr_srf%ptr     (:,:) = 0._wp

    srf_mem%sfract_srf%ptr  (:,:) = 0._wp
    srf_mem%wsn_srf%ptr     (:,:) = -1._wp    ! -1 used to determine start of experiment
    srf_mem%wsn_glacier%ptr (:,:) = 20._wp    ! Temporary initialization of glacier thickness

    ! Initialize energy fluxes
    !

    srf_mem%sw_srf_net%ptr     (:,:) = 0._wp
    srf_mem%swvis_srf_net%ptr  (:,:) = 0._wp
    srf_mem%swnir_srf_net%ptr  (:,:) = 0._wp
    srf_mem%lw_srf_net%ptr     (:,:) = 0._wp
    srf_mem%rad_srf_net%ptr    (:,:) = 0._wp
    srf_mem%latent_hflx%ptr    (:,:) = 0._wp
    srf_mem%sensible_hflx%ptr  (:,:) = 0._wp

    ! Initialize water fluxes
    !

    srf_mem%evapopot%ptr      (:,:) = 0._wp
    srf_mem%evapotrans%ptr    (:,:) = 0._wp
    srf_mem%trans%ptr         (:,:) = 0._wp

    ! Propagate initial conditions down the tile structure, starting with the tile with which this subroutine was called initially
    ! TBD: sub-tiles should not read the netcdf files again
!!$    IF (tile%no_of_sub_tiles > 0) THEN
!!$      DO i=1,tile%no_of_sub_tiles
!!$        CALL srf_init_ic(tile%sub_tiles(i))
!!$      END DO
!!$    END IF

  END SUBROUTINE srf_init_ic

#endif
END MODULE mo_srf_init
