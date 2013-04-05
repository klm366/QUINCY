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

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile
  USE mo_srf_memory_class, ONLY: t_srf_memory, t_srf_phys_state, t_srf_energy_state, t_srf_energy_flux, &
                                  t_srf_water_state, t_srf_water_flux
  USE mo_jsb_base,          ONLY: get_model
  USE mo_netcdf_read,       ONLY: netcdf_read_oncells_2d_time, netcdf_read_oncells_2d
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

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':srf_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get srf memory of the tile
    srf => tile%mem%srf

    IF (.NOT. ASSOCIATED(srf)) RETURN

    CALL message(TRIM(routine), 'Initializing srf memory for tile '//TRIM(tile%name)//' from '//TRIM(model%config%srf%ic_filename))

    ! Initialize physical state
    !
    phys => srf%phys_state

    ! Orography
    phys%oro_stddev(:,:) = 140.606_wp   ! Temporary

    ! Albedo
    return_pointer => netcdf_read_oncells_2d(         &
      & filename=model%config%srf%bc_filename,        &
      & variable_name='albedo',                       &
      & fill_array = phys%albedo_srf,                 &
      & patch=model%grid%patch)
    phys%albedo_srf(:,:) = MERGE(phys%albedo_srf(:,:), 0._wp, phys%albedo_srf(:,:) >= 0._wp)

    phys%alb_vis_srf(:,:) = phys%albedo_srf(:,:)
    phys%alb_nir_srf(:,:) = phys%albedo_srf(:,:)

    ! Initialize energy state
    !
    energy => srf%energy_state

    ! Set initial surface temperature to SST
    ! Temporary for comparing with Thomas' jsbach implementation
    t_srf_mon => netcdf_read_oncells_2d_time(         &
      & filename=model%config%srf%ic_filename,        &
      & variable_name='SST',                          &
      !& fill_array=energy%t_srf,                      &
      & patch=model%grid%patch,                       &
      & start_timestep=1, end_timestep=12)
    energy%t_srf(:,:) = t_srf_mon(:,:,1)
    t_srf => energy%t_srf
    DEALLOCATE(t_srf_mon)
    ! Note: t_srf now points to energy%t_srf and can be
    ! used below to fill %energy_state%t_srf for all other files:
    !    tile%mem%srf%energy_state%t_srf(:,:) = t_srf(:,:)

!!$    energy%t_srf_old      (:,:) = t_srf (:,:)
    energy%t_srf_unfilt   (:,:) = t_srf (:,:)
    energy%t_rad_srf      (:,:) = t_srf (:,:)
!!$    energy%t_eff_srf      (:,:) = t_srf (:,:)
    energy%s_srf          (:,:) = 0._wp
    energy%qsat_srf       (:,:) = -1._wp  ! Negative values will be used to test for start of experiment
!!$    energy%qsat_srf_old   (:,:) = 0._wp
    energy%fact_q_air     (:,:) = 0.5_wp
    energy%fact_qsat_srf  (:,:) = 0.5_wp

    ! Initialize water state

    water => srf%water_state

    water%wfract_srf (:,:) = 0._wp
    water%wsr_srf    (:,:) = 0._wp

    water%sfract_srf (:,:) = 0._wp
    water%wsn_srf    (:,:) = 0._wp

    ! Initialize energy fluxes
    !
    eflux => srf%energy_flux

    eflux%swrad_srf_net (:,:) = 0._wp
    eflux%lwrad_srf_net (:,:) = 0._wp
    eflux%rad_srf_net   (:,:) = 0._wp
    eflux%latent_hflx   (:,:) = 0._wp

    ! Initialize water fluxes
    !
    wflux => srf%water_flux

    wflux%evap_srf      (:,:) = 0._wp

    ! Propagate initial conditions down the tile structure, starting with the tile with which this subroutine was called initially
    ! TBD: sub-tiles should not read the netcdf files again
    IF (tile%no_of_sub_tiles > 0) THEN
      DO i=1,tile%no_of_sub_tiles
        CALL srf_init_ic(tile%sub_tiles(i))
      END DO
    END IF

  END SUBROUTINE srf_init_ic

END MODULE mo_srf_init
