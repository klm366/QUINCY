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

  dsl4jsb_Use_processes SRF, HD
  dsl4jsb_Use_config(SRF)
  dsl4jsb_Use_config(HD)
  dsl4jsb_Use_memory(SRF)

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
    dsl4jsb_Def_config(SRF)
    dsl4jsb_Def_config(HD)
    dsl4jsb_Def_memory(SRF)

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
    dsl4jsb_Get_memory(SRF)

    ! Get srf config
    dsl4jsb_Get_config(SRF)
    dsl4jsb_Get_config(HD)

    IF (.NOT. dsl4jsb_Config(SRF)%active) RETURN

    ! Initialize parameters
    !
    CALL message(TRIM(routine), &
         'Initializing srf memory for tile '//TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(SRF)%bc_sso_filename))
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SRF)%bc_sso_filename), model%id)

    ! Orography
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='orostd',                   &
      & fill_array = dsl4jsb_var2D_ptr(SRF,oro_stddev))
    dsl4jsb_var2D_onDomain(SRF,oro_stddev) = &
      MERGE(dsl4jsb_var2D_onDomain(SRF,oro_stddev), 0._wp, dsl4jsb_var2D_onDomain(SRF,oro_stddev) >= 0._wp)

    CALL jsb_netcdf_close_file(input_file)

    CALL message(TRIM(routine), 'Initializing srf memory for tile '//TRIM(tile%name)//' from bc_land_frac.nc')
    input_file = jsb_netcdf_open_input('bc_land_frac.nc', model%id)     ! Temporary
!!$    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SRF)%bc_filename), model%id)

    ! Land mask (including glacies and lakes)
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='notsea',                   &        !< ocean 0., land 1.
      & fill_array = dsl4jsb_var2D_ptr(SRF,land))

    ! Glacier mask
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='glac',                     &
      & fill_array = dsl4jsb_var2D_ptr(SRF,glacier))
    ! set glacier mask to zero on ocean grid cells (instead of missing value)
    WHERE (dsl4jsb_var2D_onDomain(SRF,land) == 0._wp)
      dsl4jsb_var2D_onDomain(SRF,glacier) = 0._wp
    END WHERE

    ! Lake mask
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='lake',                     &
      & fill_array = dsl4jsb_var2D_ptr(SRF,lake))

    CALL jsb_netcdf_close_file(input_file)

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SRF)%bc_filename), model%id)

    ! Surface roughness length
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='roughness_length',         &
      & fill_array = dsl4jsb_var2D_ptr(SRF,rough_m_srf))
    dsl4jsb_var2D_onDomain(SRF,rough_m_srf) = &
      MERGE(dsl4jsb_var2D_onDomain(SRF,rough_m_srf), 1._wp, dsl4jsb_var2D_onDomain(SRF,rough_m_srf) >= 0._wp)

    ! Albedo
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='albedo',                   &
      & fill_array = dsl4jsb_var2D_ptr(SRF,albedo_srf))
    dsl4jsb_var2D_onDomain(SRF,albedo_srf) = &
      MERGE(dsl4jsb_var2D_onDomain(SRF,albedo_srf), 0._wp, dsl4jsb_var2D_onDomain(SRF,albedo_srf) >= 0._wp)

    dsl4jsb_var2D_onDomain(SRF,alb_vis_srf) = dsl4jsb_var2D_onDomain(SRF,albedo_srf)
    dsl4jsb_var2D_onDomain(SRF,alb_nir_srf) = dsl4jsb_var2D_onDomain(SRF,albedo_srf)

    CALL jsb_netcdf_close_file(input_file)

    ! Initialize energy state
    !

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SRF)%ic_filename), model%id)

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
    dsl4jsb_var2D_onDomain(SRF,t_srf) = wgt1 * t_srf_mon(:,:,n1) + wgt2 * t_srf_mon(:,:,n2)

    DEALLOCATE(t_srf_mon)

    CALL jsb_netcdf_close_file(input_file)

    t_srf => dsl4jsb_var2D_ptr(SRF,t_srf)
    ! Note: t_srf now points to srf_mem%t_srf%ptr and can be
    ! used below to fill %t_srf for all other files:
    !    tile%mem%srf%energy_state%t_srf%ptr(:,:) = t_srf(:,:)

    dsl4jsb_var2D_onDomain(SRF,t_srf_old)            = t_srf (:,:)
    dsl4jsb_var2D_onDomain(SRF,t_srf_unfilt)         = t_srf (:,:)
    dsl4jsb_var2D_onDomain(SRF,t_srf_unfilt_old)     = t_srf (:,:)
    dsl4jsb_var2D_onDomain(SRF,t_eff_srf)            = t_srf (:,:)
!!$    srf_mem%t_eff_srf%ptr        (:,:) = t_srf (:,:)
    dsl4jsb_var2D_onDomain(SRF,s_srf)                = 0._wp
    dsl4jsb_var2D_onDomain(SRF,qsat_srf)             = -1._wp  ! Negative values will be used to test for start of experiment
!!$    srf_mem%qsat_srf_old%ptr     (:,:) = 0._wp
    dsl4jsb_var2D_onDomain(SRF,q_snocpymlt)          = 0._wp
    dsl4jsb_var2D_onDomain(SRF,fact_q_air)           = 0.5_wp
    dsl4jsb_var2D_onDomain(SRF,fact_qsat_srf)        = 0.5_wp
    dsl4jsb_var2D_onDomain(SRF,fact_qsat_trans_srf)  = 0.5_wp

    ! Initialize water state

    dsl4jsb_var2D_onDomain(SRF,wfract_srf)   = 0._wp
    dsl4jsb_var2D_onDomain(SRF,wsr_srf)      = 0._wp
    dsl4jsb_var2D_onDomain(SRF,sfract_srf)   = 0._wp
    dsl4jsb_var2D_onDomain(SRF,wsn_srf)      = -1._wp    ! -1 used to determine start of experiment
    dsl4jsb_var2D_onDomain(SRF,wsn_glacier)  = 20._wp    ! Temporary initialization of glacier thickness

    ! Initialize energy fluxes
    !

    dsl4jsb_var2D_onDomain(SRF,sw_srf_net)      = 0._wp
    dsl4jsb_var2D_onDomain(SRF,swvis_srf_net)   = 0._wp
    dsl4jsb_var2D_onDomain(SRF,swnir_srf_net)   = 0._wp
    dsl4jsb_var2D_onDomain(SRF,lw_srf_net)      = 0._wp
    dsl4jsb_var2D_onDomain(SRF,rad_srf_net)     = 0._wp
    dsl4jsb_var2D_onDomain(SRF,latent_hflx)     = 0._wp
    dsl4jsb_var2D_onDomain(SRF,sensible_hflx)   = 0._wp

    ! Initialize water fluxes
    !

    dsl4jsb_var2D_onDomain(SRF,evapopot)       = 0._wp
    dsl4jsb_var2D_onDomain(SRF,evapotrans)     = 0._wp
    dsl4jsb_var2D_onDomain(SRF,trans)          = 0._wp

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
