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
  USE mo_jsb_io_netcdf,     ONLY: jsb_netcdf_read_2d, jsb_netcdf_read_2d_time
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: srf_init_ic, srf_init_bc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_srf_init'

CONTAINS

  SUBROUTINE srf_init_bc(tile)

    TYPE(t_jsb_tile), INTENT(in) :: tile

    TYPE(t_jsb_model),       POINTER :: model
    TYPE(t_srf_memory), POINTER :: srf

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':srf_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get srf memory of the tile
    srf => tile%mem%srf

  END SUBROUTINE srf_init_bc

  RECURSIVE SUBROUTINE srf_init_ic(tile)

    USE mo_jsb_time, ONLY: wgt1, wgt2, nmw1,nmw2

    TYPE(t_jsb_tile), INTENT(in) :: tile

    TYPE(t_jsb_model),        POINTER :: model
    TYPE(t_srf_memory),       POINTER :: srf
    TYPE(t_srf_phys_state),   POINTER :: phys
    TYPE(t_srf_energy_state), POINTER :: energy
    TYPE(t_srf_energy_flux),  POINTER :: eflux
    TYPE(t_srf_water_state),  POINTER :: water
    TYPE(t_srf_water_flux),   POINTER :: wflux

    REAL(wp), POINTER :: t_srf(:,:), t_srf_mon(:,:,:)
    REAL(wp), POINTER :: &
      & return_pointer (:,  :)       !< temporary pointer

    INTEGER :: i, n1, n2

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':srf_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get srf memory of the tile
    srf => tile%mem%srf

    IF (.NOT. ASSOCIATED(srf)) RETURN

    CALL message(TRIM(routine), 'Initializing srf memory for tile '//TRIM(tile%name)//' from '//TRIM(model%config%srf%ic_filename))

    ! Initialize physical state
    !
    phys => srf%phys_state

    ! Glacier mask
    return_pointer => jsb_netcdf_read_2d(         &
!!$      & filename=model%config%srf%bc_filename,        &
      & filename='bc_land_frac.nc',                   & ! Temporary
      & variable_name='glac',                         &
      & fill_array = phys%glacier,                    &
      & model_id=model%id)
    phys%glacier(:,:) = MERGE(phys%glacier(:,:), 0._wp, phys%glacier(:,:) >= 0._wp)

    ! Orography
    return_pointer => jsb_netcdf_read_2d(         &
      & filename=model%config%srf%bc_sso_filename,    &
      & variable_name='orostd',                       &
      & fill_array = phys%oro_stddev,                 &
      & model_id=model%id)
    phys%oro_stddev(:,:) = MERGE(phys%oro_stddev(:,:), 0._wp, phys%oro_stddev(:,:) >= 0._wp)

    ! Surface roughness length
    return_pointer => jsb_netcdf_read_2d(         &
      & filename=model%config%srf%bc_filename,        &
      & variable_name='roughness_length',             &
      & fill_array = phys%rough_m_srf,                &
      & model_id=model%id)
    phys%rough_m_srf(:,:) = MERGE(phys%rough_m_srf(:,:), 1._wp, phys%rough_m_srf(:,:) >= 0._wp)

    ! Albedo
    return_pointer => jsb_netcdf_read_2d(         &
      & filename=model%config%srf%bc_filename,        &
      & variable_name='albedo',                       &
      & fill_array = phys%albedo_srf,                 &
      & model_id=model%id)
    phys%albedo_srf(:,:) = MERGE(phys%albedo_srf(:,:), 0._wp, phys%albedo_srf(:,:) >= 0._wp)

    phys%alb_vis_srf(:,:) = phys%albedo_srf(:,:)
    phys%alb_nir_srf(:,:) = phys%albedo_srf(:,:)

    ! Initialize energy state
    !
    energy => srf%energy_state

    ! Initial surface temperature
!!$    return_pointer => jsb_netcdf_read_2d(         &
!!$      & filename=model%config%srf%ic_filename,        &
!!$      & variable_name='surface_temperature',          &
!!$      fill_array=energy%t_srf,                        &
!!$      & model_id=model%id)
!!$    t_srf => energy%t_srf
    t_srf_mon => jsb_netcdf_read_2d_time(         &   ! Temporary
      & filename=model%config%srf%ic_filename,        &
      & variable_name='surf_temp',                    &
      & model_id=model%id)
    t_srf_mon(:,:,:) = MERGE(t_srf_mon(:,:,:), 280._wp, t_srf_mon(:,:,:) > 0._wp)
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
    t_srf => energy%t_srf
    ! Note: t_srf now points to energy%t_srf and can be
    ! used below to fill %energy_state%t_srf for all other files:
    !    tile%mem%srf%energy_state%t_srf(:,:) = t_srf(:,:)

    energy%t_srf_old        (:,:) = t_srf (:,:)
    energy%t_srf_unfilt     (:,:) = t_srf (:,:)
    energy%t_srf_unfilt_old (:,:) = t_srf (:,:)
    energy%t_rad_srf        (:,:) = t_srf (:,:)
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

    water%rel_hum_srf (:,:) = 0._wp

    ! Initialize energy fluxes
    !
    eflux => srf%energy_flux

    eflux%swrad_srf_down (:,:) = 0._wp
    eflux%lwrad_srf_net  (:,:) = 0._wp
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
