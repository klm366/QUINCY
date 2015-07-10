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

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile
  USE mo_srf_memory_class, ONLY: t_srf_memory, t_srf_phys_state, t_srf_energy_state, t_srf_energy_flux, &
                                  t_srf_water_state, t_srf_water_flux
  USE mo_jsb_base,          ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_close_file, jsb_netcdf_open_input, &
                                  jsb_netcdf_read_2d, jsb_netcdf_read_2d_time
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: srf_init_ic, srf_init_bc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_srf_init'

CONTAINS

  SUBROUTINE srf_init_bc(tile)

    TYPE(t_jsb_tile), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model
    TYPE(t_srf_memory), POINTER :: srf

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':srf_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get srf memory of the tile
    srf => tile%mem%srf

  END SUBROUTINE srf_init_bc

  RECURSIVE SUBROUTINE srf_init_ic(tile)

    USE mo_jsb_time, ONLY: get_time_interpolation_weights

    TYPE(t_jsb_tile), INTENT(inout) :: tile

    TYPE(t_jsb_model),        POINTER :: model
    TYPE(t_srf_memory),       POINTER :: srf
    TYPE(t_srf_phys_state),   POINTER :: phys
    TYPE(t_srf_energy_state), POINTER :: energy
    TYPE(t_srf_energy_flux),  POINTER :: eflux
    TYPE(t_srf_water_state),  POINTER :: water
    TYPE(t_srf_water_flux),   POINTER :: wflux

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
    srf => tile%mem%srf

    IF (.NOT. ASSOCIATED(srf)) RETURN

    ! Initialize physical state
    !
    phys => srf%phys_state

    CALL message(TRIM(routine), &
         'Initializing srf memory for tile '//TRIM(tile%name)//' from '//TRIM(model%config%srf%bc_sso_filename))
    input_file = jsb_netcdf_open_input(TRIM(model%config%srf%bc_sso_filename), model%id)

    ! Orography
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='orostd',                   &
      & fill_array = phys%oro_stddev)
    phys%oro_stddev(:,:) = MERGE(phys%oro_stddev(:,:), 0._wp, phys%oro_stddev(:,:) >= 0._wp)

    CALL jsb_netcdf_close_file(input_file)

    CALL message(TRIM(routine), 'Initializing srf memory for tile '//TRIM(tile%name)//' from bc_land_frac.nc')
    input_file = jsb_netcdf_open_input('bc_land_frac.nc', model%id)     ! Temporary
!!$    input_file = jsb_netcdf_open_input(TRIM(model%config%srf%bc_filename), model%id)

    ! Land mask (including glacies and lakes)
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='notsea',                   &        !< ocean 0., land 1.
      & fill_array = phys%land)

    ! Glacier mask
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='glac',                     &
      & fill_array = phys%glacier)
    ! set glacier mask to zero on ocean grid cells (instead of missing value)
    WHERE (phys%land(:,:) == 0._wp)
      phys%glacier(:,:) = 0._wp
    END WHERE

    ! Lake mask
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='lake',                     &
      & fill_array = phys%lake)

    CALL jsb_netcdf_close_file(input_file)

    input_file = jsb_netcdf_open_input('r2b4_amip.nc', model%id)    !vg Name !!!

    ! coastal cells of the land and ocean
    cell_sea_land_mask => jsb_netcdf_read_2d( &
      & input_file,                                &
      & variable_name='cell_sea_land_mask')
    phys%coast_land(:,:)  = MERGE(1._wp, 0._wp, cell_sea_land_mask(:,:) == 1)
    phys%coast_ocean(:,:) = MERGE(1._wp, 0._wp, cell_sea_land_mask(:,:) == -1)
    DEALLOCATE(cell_sea_land_mask)

    CALL jsb_netcdf_close_file(input_file)

    input_file = jsb_netcdf_open_input(TRIM(model%config%srf%bc_filename), model%id)

    ! Surface roughness length
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='roughness_length',         &
      & fill_array = phys%rough_m_srf)
    phys%rough_m_srf(:,:) = MERGE(phys%rough_m_srf(:,:), 1._wp, phys%rough_m_srf(:,:) >= 0._wp)

    ! Albedo
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='albedo',                   &
      & fill_array = phys%albedo_srf)
    phys%albedo_srf(:,:) = MERGE(phys%albedo_srf(:,:), 0._wp, phys%albedo_srf(:,:) >= 0._wp)

    phys%alb_vis_srf(:,:) = phys%albedo_srf(:,:)
    phys%alb_nir_srf(:,:) = phys%albedo_srf(:,:)

    CALL jsb_netcdf_close_file(input_file)

    ! Initialize energy state
    !
    energy => srf%energy_state

    input_file = jsb_netcdf_open_input(TRIM(model%config%srf%ic_filename), model%id)

    ! Initial surface temperature
!!$    return_pointer => jsb_netcdf_read_2d(         &
!!$      & filename=model%config%srf%ic_filename,        &
!!$      & variable_name='surface_temperature',          &
!!$      fill_array=energy%t_srf,                        &
!!$      & model_id=model%id)
!!$    t_srf => energy%t_srf
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
    energy%t_srf(:,:) = wgt1 * t_srf_mon(:,:,n1) + wgt2 * t_srf_mon(:,:,n2)

    DEALLOCATE(t_srf_mon)

    CALL jsb_netcdf_close_file(input_file)

    t_srf => energy%t_srf
    ! Note: t_srf now points to energy%t_srf and can be
    ! used below to fill %energy_state%t_srf for all other files:
    !    tile%mem%srf%energy_state%t_srf(:,:) = t_srf(:,:)

    energy%t_srf_old        (:,:) = t_srf (:,:)
    energy%t_srf_unfilt     (:,:) = t_srf (:,:)
    energy%t_srf_unfilt_old (:,:) = t_srf (:,:)
    energy%t_eff_srf        (:,:) = t_srf (:,:)
!!$    energy%t_eff_srf        (:,:) = t_srf (:,:)
    energy%s_srf            (:,:) = 0._wp
    energy%qsat_srf         (:,:) = -1._wp  ! Negative values will be used to test for start of experiment
!!$    energy%qsat_srf_old     (:,:) = 0._wp
    energy%tte_corr         (:,:) = 0._wp
    energy%fact_q_air       (:,:) = 0.5_wp
    energy%fact_qsat_srf    (:,:) = 0.5_wp
    energy%fact_qsat_trans_srf (:,:) = 0.5_wp

    ! Initialize water state

    water => srf%water_state

    water%wfract_srf  (:,:) = 0._wp
    water%wsr_srf     (:,:) = 0._wp

    water%sfract_srf  (:,:) = 0._wp
    water%wsn_srf     (:,:) = -1._wp    ! -1 used to determine start of experiment
    water%wsn_glacier (:,:) = 20._wp    ! Temporary initialization of glacier thickness

    ! Initialize energy fluxes
    !
    eflux => srf%energy_flux

    eflux%sw_srf_net     (:,:) = 0._wp
    eflux%swvis_srf_net  (:,:) = 0._wp
    eflux%swnir_srf_net  (:,:) = 0._wp
    eflux%lw_srf_net     (:,:) = 0._wp
    eflux%rad_srf_net    (:,:) = 0._wp
    eflux%latent_hflx    (:,:) = 0._wp
    eflux%sensible_hflx  (:,:) = 0._wp

    ! Initialize water fluxes
    !
    wflux => srf%water_flux

    wflux%evapopot      (:,:) = 0._wp
    wflux%evapotrans    (:,:) = 0._wp
    wflux%trans         (:,:) = 0._wp

    ! Propagate initial conditions down the tile structure, starting with the tile with which this subroutine was called initially
    ! TBD: sub-tiles should not read the netcdf files again
    IF (tile%no_of_sub_tiles > 0) THEN
      DO i=1,tile%no_of_sub_tiles
        CALL srf_init_ic(tile%sub_tiles(i))
      END DO
    END IF

  END SUBROUTINE srf_init_ic

#endif
END MODULE mo_srf_init
