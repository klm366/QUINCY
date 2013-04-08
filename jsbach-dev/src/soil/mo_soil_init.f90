!>
!! @brief Initialization of the the soil memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_soil_init

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile
  USE mo_soil_memory_class, ONLY: t_soil_memory, t_soil_phys_state, t_soil_energy_state, t_soil_energy_flux, &
                                  t_soil_water_state, t_soil_water_flux
  USE mo_jsb_base,          ONLY: get_model
  USE mo_netcdf_read,       ONLY: netcdf_read_oncells_2d_time, netcdf_read_oncells_2d
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: soil_init_ic, soil_init_bc

  REAL(wp), PARAMETER :: &
    fao_thermal_diff(0:6) = (/7.4e-7_wp ,8.7e-7_wp ,8.0e-7_wp ,7.4e-7_wp ,7.1e-7_wp ,6.7e-7_wp ,-1._wp/)

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_init'

CONTAINS

  SUBROUTINE soil_init_bc(tile)

    TYPE(t_jsb_tile), INTENT(in) :: tile

    TYPE(t_jsb_model),       POINTER :: model
    TYPE(t_soil_memory), POINTER :: soil

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get soil memory of the tile
    soil => tile%mem%soil

  END SUBROUTINE soil_init_bc

  RECURSIVE SUBROUTINE soil_init_ic(tile)

    USE mo_soil_util, ONLY: soil_depth_to_layers_2d

    TYPE(t_jsb_tile), INTENT(in) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    TYPE(t_soil_memory),       POINTER :: soil
    TYPE(t_soil_phys_state),   POINTER :: phys
    TYPE(t_soil_energy_state), POINTER :: energy
    TYPE(t_soil_energy_flux),  POINTER :: eflux
    TYPE(t_soil_water_state),  POINTER :: water
    TYPE(t_soil_water_flux),   POINTER :: wflux

    INTEGER :: isoil
    REAL(wp), POINTER :: t_srf(:,:), t_srf_mon(:,:,:)
    REAL(wp), POINTER :: &
      & vpwp           (:,  :),    & !< Volumetric soil permanent wilting point [m/m]
      & vpor           (:,  :),    & !< Volumetric soil porosity following Cosby et al. [m/m]
      & vfc            (:,  :),    & !< Volumetric soil field capacity [m/m]
      & fao            (:,  :),    & !< FAO soil types
      & return_pointer (:,  :)       !< temporary pointer

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get soil memory of the tile
    soil => tile%mem%soil

    IF (.NOT. ASSOCIATED(soil)) RETURN

    CALL message(TRIM(routine), 'Initializing soil memory for tile '//TRIM(tile%name))

    ! Initialize physical state
    !
    phys => soil%phys_state

    ! Total soil depth until bedrock

!!$    return_pointer => netcdf_read_oncells_2d(         &
!!$      & filename=model%config%soil%bc_filename,       &
!!$      & variable_name='soil_depth',                   &
!!$      & fill_array = phys%soil_depth, &
!!$      & patch=model%grid%patch)
    phys%soil_depth(:,:) = 2.72603_wp   ! Temporary

    ! Compute actual soil layer depths based on the fixed layer thicknesses from namelist and
    ! the bedrock depth just read in.
    phys%soil_depth_l(:,:,:) = soil_depth_to_layers_2d( &
      & phys%soil_depth(:,:),                         & ! Total soil depth until bedrock (from textures)
      & model%config%soil%dz_water(1:model%config%soil%nsoil_water))    ! Soil layer thicknesses from namelist

    ! Volumetric soil porosity
!!$    return_pointer => netcdf_read_oncells_2d(            &
!!$      & filename=model%config%soil%bc_filename,   &
!!$      & fill_array = phys%vpor,       &
!!$      & variable_name='soil_porosity',            &
!!$      & patch=model%grid%patch)
    phys%vpor(:,:) = 0.441515_wp   ! Temporary

    ! Saturated hydraulic conductivity
!!$    return_pointer => netcdf_read_oncells_2d(         &
!!$      & filename=model%config%soil%bc_filename,       &
!!$      & variable_name='hyd_cond_sat',                 &
!!$      & fill_array = phys%shc,        &
!!$      & patch=model%grid%patch)
    phys%shc(:,:) = 4.60588E-6_wp  ! Temporary

    ! Matrix potential
!!$    return_pointer => netcdf_read_oncells_2d(         &
!!$      & filename=model%config%soil%bc_filename,       &
!!$      & variable_name='moisture_pot',                 &
!!$      & fill_array = phys%mpot,       &
!!$      & patch=model%grid%patch)
    phys%mpot(:,:) = 0.275715_wp  ! Temporary

    ! Exponent B in Clapp and Hornberger
!!$    return_pointer => netcdf_read_oncells_2d(         &
!!$      & filename=model%config%soil%bc_filename,       &
!!$      & variable_name='bclapp',                       &
!!$      & fill_array = phys%bclapp,     &
!!$      & patch=model%grid%patch)
    phys%bclapp(:,:) = 8.58824_wp  ! Temporary

    ! Pore size distribution index
!!$    return_pointer => netcdf_read_oncells_2d(         &
!!$      & filename=model%config%soil%bc_filename,       &
!!$      & variable_name='pore_size_index',              &
!!$      & fill_array = phys%psdi,       &
!!$      & patch=model%grid%patch)
    phys%psdi(:,:) = 0.20747_wp  ! Temporary

    ! Thermal diffusivity of soil
!!$    fao => netcdf_read_oncells_2d(                    &
!!$      & filename=model%config%soil%bc_filename,       &
!!$      & variable_name='fao',                          &
!!$      & patch=model%grid%patch)
!!$    fao(:,:) = MAX(REAL(LBOUND(fao_thermal_diff,1),wp), MIN(REAL(UBOUND(fao_thermal_diff,1),wp), fao(:,:)))
!!$    phys%hdiff(:,:) = fao_thermal_diff(NINT(fao(:,:)))
    phys%hdiff(:,:) = 7.E-7_wp  ! Temporary

    ! Volumetric heat capacity of dry soil
!!$    return_pointer => netcdf_read_oncells_2d(         &
!!$      & filename=model%config%soil%bc_filename,       &
!!$      & variable_name='heat_capacity',                &
!!$      & fill_array = phys%vhcap,      &
!!$      & patch=model%grid%patch)
    phys%vhcap(:,:) = 2.E+6_wp  ! Temporary

    return_pointer => netcdf_read_oncells_2d(         &
      & filename=model%config%soil%bc_filename,       &
      & variable_name='albedo_soil_vis',              &
      & fill_array = phys%alb_vis_soil,               &
      & patch=model%grid%patch)
    phys%alb_vis_soil(:,:) = MERGE(phys%alb_vis_soil(:,:), 0._wp, phys%alb_vis_soil(:,:) >= 0._wp)

    return_pointer => netcdf_read_oncells_2d(         &
      & filename=model%config%soil%bc_filename,       &
      & variable_name='albedo_soil_nir',              &
      & fill_array = phys%alb_nir_soil,               &
      & patch=model%grid%patch)
    phys%alb_nir_soil(:,:) = MERGE(phys%alb_nir_soil(:,:), 0._wp, phys%alb_nir_soil(:,:) >= 0._wp)

    ! Initialize energy state
    !
    energy => soil%energy_state

    DO isoil=1,model%config%soil%nsoil_energy
      energy%t_soil       (:,isoil,:) = tile%mem%srf%energy_state%t_srf(:,:) ! Temporary. Later call to initemp
      energy%t_soil_acoef (:,isoil,:) = 0._wp
      energy%t_soil_bcoef (:,isoil,:) = 1._wp
    END DO

    energy%hcap_grnd (:,:) = 1._wp

    ! Initialize water state

    water => soil%water_state

    water%wfract_soil (:,:) = 0._wp
    water%wsr_soil    (:,:) = 0._wp

    water%sfract_soil (:,:) = 0._wp
    water%wsn_soil    (:,:) = 0._wp

    ! Temporary:
    water%ws_l(:,1,:) = 0.0109355_wp
    water%ws_l(:,2,:) = 0.0427327_wp
    water%ws_l(:,3,:) = 0.153602_wp
    water%ws_l(:,4,:) = 0.251354_wp
    water%ws_l(:,5,:) = 0._wp
    water%ws  (:,  :) = 0.446431_wp

    ! Water content at field capacity
    !
    ! Read volumetric field capacity from boundary condition file
!!$    vfc => netcdf_read_oncells_2d(             &
!!$      & filename=model%config%soil%bc_filename,   &
!!$      & variable_name='soil_field_cap',           &
!!$      & patch=model%grid%patch)

    ! Get field capacity by multiplying with soil layer thicknesses
    DO i=1,model%config%soil%nsoil_water
!!$      water%ws_fc_l(:,i,:) = vfc(:,:) * phys%soil_depth_l(:,i,:)
      water%ws_fc_l(:,i,:) = 0.318517_wp * phys%soil_depth_l(:,i,:) ! Temporary
    END DO
!!$    DEALLOCATE(vfc)

    ! Water content at permanent wilting point
    !
    ! Read volumetric wilting point from boundary condition file
!!$    vpwp => netcdf_read_oncells_2d(            &
!!$      & filename=model%config%soil%bc_filename,   &
!!$      & variable_name='wilting_point',            &
!!$      & patch=model%grid%patch)
    ! Get permanent wilting point by multiplying with soil layer thicknesses
    DO i=1,model%config%soil%nsoil_water
!!$      water&ws_pwp_l(:,i,:) = vpwp(:,:) * phys%soil_depth_l(:,i,:)
      water%ws_pwp_l(:,i,:) = 0.185245_wp * phys%soil_depth_l(:,i,:)     ! Temporary
    END DO
!!$    DEALLOCATE(vpwp)

    ! Water content at saturation (water holding capacity)
    !
    ! Get water holding capacity by multiplying volumetric porosity with soil layer thicknesses
    DO i=1,model%config%soil%nsoil_water
      water%ws_sat_l(:,i,:) = phys%vpor(:,:) * phys%soil_depth_l(:,i,:)
    END DO

    ! Initialize energy fluxes
    !
    eflux => soil%energy_flux

    eflux%grnd_hflx     (:,:) = 0._wp

    ! Initialize water fluxes
    !
    wflux => soil%water_flux

    wflux%runoff        (:,:) = 0._wp
    wflux%drainage      (:,:) = 0._wp

    ! Propagate initial conditions down the tile structure, starting with the tile with which this subroutine was called initially
    ! TBD: sub-tiles should not read the netcdf files again
    IF (tile%no_of_sub_tiles > 0) THEN
      DO i=1,tile%no_of_sub_tiles
        CALL soil_init_ic(tile%sub_tiles(i))
      END DO
    END IF

  END SUBROUTINE soil_init_ic

END MODULE mo_soil_init
