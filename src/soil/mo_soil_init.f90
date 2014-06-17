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
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile
  USE mo_soil_memory_class, ONLY: t_soil_memory, t_soil_phys_state, t_soil_energy_state, t_soil_energy_flux, &
                                  t_soil_water_state, t_soil_water_flux
  USE mo_jsb_base,          ONLY: get_model
  USE mo_netcdf_read,       ONLY: netcdf_read_oncells_2d_time, netcdf_read_oncells_2d, netcdf_read_oncells_2d_extdim
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

    IF (.NOT. ASSOCIATED(soil)) RETURN

    CALL message(TRIM(routine), 'Setting  soil boundary conditions for tile '//TRIM(tile%name))

  END SUBROUTINE soil_init_bc

  RECURSIVE SUBROUTINE soil_init_ic(tile)

    USE mo_soil_process, ONLY: init_soil_temperature
    USE mo_soil_util,    ONLY: soil_depth_to_layers_2d

    TYPE(t_jsb_tile), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    TYPE(t_soil_memory),       POINTER :: soil
    TYPE(t_soil_phys_state),   POINTER :: phys
    TYPE(t_soil_energy_state), POINTER :: energy
    TYPE(t_soil_energy_flux),  POINTER :: eflux
    TYPE(t_soil_water_state),  POINTER :: water
    TYPE(t_soil_water_flux),   POINTER :: wflux

    INTEGER :: isoil
    REAL(wp), POINTER :: &
      & t_srf_mon      (:,:,:),    &
      & vpwp           (:,  :),    & !< Volumetric soil permanent wilting point [m/m]
      & vpor           (:,  :),    & !< Volumetric soil porosity following Cosby et al. [m/m]
      & vfc            (:,  :),    & !< Volumetric soil field capacity [m/m]
      & fao            (:,  :),    & !< FAO soil types
      & return_pointer (:,  :),    & !< temporary pointer
      & return_pointer3d(:,:,:)

    INTEGER :: i, n1, n2

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

    return_pointer => netcdf_read_oncells_2d(         &
      & filename=model%config%soil%bc_filename,       &
      & variable_name='soildepth',                    &
      & fill_array = phys%soil_depth, &
      & patch=model%grid%patch)
    phys%soil_depth(:,:) = MAX(0.1_wp, MERGE(phys%soil_depth(:,:), 0._wp, phys%soil_depth(:,:) >= 0._wp))

    ! Compute actual soil layer depths based on the fixed layer thicknesses from namelist and
    ! the bedrock depth just read in.
    phys%soil_depth_l(:,:,:) = soil_depth_to_layers_2d( &
      & phys%soil_depth(:,:),                         & ! Total soil depth until bedrock (from textures)
      & model%config%soil%dz_water(1:model%config%soil%nsoil_water))    ! Soil layer thicknesses from namelist

    ! Volumetric soil porosity
    return_pointer => netcdf_read_oncells_2d(         &
      & filename=model%config%soil%bc_filename,       &
      & fill_array = phys%vpor,                       &
      & variable_name='porosity',                     &
      & patch=model%grid%patch)
    phys%vpor(:,:) = MERGE(phys%vpor(:,:), 0.2_wp, phys%vpor(:,:) >= 0._wp)

    ! Saturated hydraulic conductivity
    return_pointer => netcdf_read_oncells_2d(         &
      & filename=model%config%soil%bc_filename,       &
      & variable_name='hydcond',                      &
      & fill_array = phys%shc,                        &
      & patch=model%grid%patch)
    phys%shc(:,:) = MERGE(phys%shc(:,:), 4.e-6_wp, phys%shc(:,:) >= 0._wp)

    ! Matrix potential
    return_pointer => netcdf_read_oncells_2d(         &
      & filename=model%config%soil%bc_filename,       &
      & variable_name='moisture_pot',                 &
      & fill_array = phys%mpot,       &
      & patch=model%grid%patch)
    phys%mpot(:,:) = MERGE(phys%mpot(:,:), 0.2_wp, phys%mpot(:,:) >= 0._wp)

    ! Exponent B in Clapp and Hornberger
    return_pointer => netcdf_read_oncells_2d(         &
      & filename=model%config%soil%bc_filename,       &
      & variable_name='bclapp',                       &
      & fill_array = phys%bclapp,                     &
      & patch=model%grid%patch)
    phys%bclapp(:,:) = MERGE(phys%bclapp(:,:), 6._wp, phys%bclapp(:,:) >= 0._wp)

    ! Pore size distribution index
    return_pointer => netcdf_read_oncells_2d(         &
      & filename=model%config%soil%bc_filename,       &
      & variable_name='poresize',                     &
      & fill_array = phys%psdi,                       &
      & patch=model%grid%patch)
    phys%psdi(:,:) = MERGE(phys%psdi(:,:), 0.25_wp, phys%psdi(:,:) >= 0._wp)

    ! Thermal diffusivity of soil
    fao => netcdf_read_oncells_2d(                    &
      & filename=model%config%soil%bc_filename,       &
      & variable_name='fao',                          &
      & patch=model%grid%patch)
    fao(:,:) = MERGE(fao(:,:), 0._wp, fao(:,:) >= 0._wp)
    fao(:,:) = MAX(REAL(LBOUND(fao_thermal_diff,1),wp), MIN(REAL(UBOUND(fao_thermal_diff,1),wp), fao(:,:)))

    DO i=1,SIZE(fao,2)
      phys%hdiff(:,i) = fao_thermal_diff(NINT(fao(:,i)))
    END DO
!!$    phys%hdiff(:,:) = RESHAPE(fao_thermal_diff(PACK(NINT(fao(:,:)), fao(:,:) > -1._wp)), /( SIZE(fao(:,:),1),SIZE(fao(:,:),2) )/ )

    ! Volumetric heat capacity of dry soil
    return_pointer => netcdf_read_oncells_2d(         &
      & filename=model%config%soil%bc_filename,       &
      & variable_name='heatcapacity',                 &
      & fill_array = phys%vhcap,                      &
      & patch=model%grid%patch)
    phys%vhcap(:,:) = MERGE(phys%vhcap(:,:), 2.e6_wp, phys%vhcap(:,:) >= 0._wp)

    ! Albedo
    return_pointer => netcdf_read_oncells_2d(         &
!!$      & filename=model%config%soil%bc_filename,       &
      & filename='bc_land_phys.nc',                   &  ! Temporary
      & variable_name='albedo_soil_vis',              &
      & fill_array = phys%alb_vis_soil,               &
      & patch=model%grid%patch)
    phys%alb_vis_soil(:,:) = MERGE(phys%alb_vis_soil(:,:), 0._wp, phys%alb_vis_soil(:,:) >= 0._wp)

    return_pointer => netcdf_read_oncells_2d(         &
!!$      & filename=model%config%soil%bc_filename,       &
      & filename='bc_land_phys.nc',                   &  ! Temporary
      & variable_name='albedo_soil_nir',              &
      & fill_array = phys%alb_nir_soil,               &
      & patch=model%grid%patch)
    phys%alb_nir_soil(:,:) = MERGE(phys%alb_nir_soil(:,:), 0._wp, phys%alb_nir_soil(:,:) >= 0._wp)

    ! Initialize energy state
    !
    energy => soil%energy_state

    ! Soil temperatures
    !
    ! Read climatological surface temperature
    t_srf_mon => netcdf_read_oncells_2d_time(         &   ! Temporary
      & filename=model%config%srf%ic_filename,        &
      & variable_name='surf_temp',                    &
      & patch=model%grid%patch)
    t_srf_mon(:,:,:) = MERGE(t_srf_mon(:,:,:), 280._wp, t_srf_mon(:,:,:) > 0._wp)

    DO i=1,SIZE(phys%soil_depth(:,:), DIM=2)
      CALL init_soil_temperature( t_srf_mon(:,i,1:12), &
                                & energy%t_soil(:,1:model%config%soil%nsoil_energy,i), &
                                & tile%mem%srf%energy_state%t_srf(:,i))
    END DO
    energy%t_soil_acoef (:,:,:) = 0._wp
    energy%t_soil_bcoef (:,:,:) = 1._wp

    energy%hcap_grnd (:,:) = 0._wp

    ! Initialize water state

    water => soil%water_state

    water%wfract_soil (:,:) = 0._wp
    water%wsr_soil    (:,:) = 0._wp

    water%sfract_soil (:,:) = 0._wp

    ! Initial snow depth
    return_pointer => netcdf_read_oncells_2d(         &
      & filename=model%config%soil%ic_filename,       &
      & variable_name='snow',                         &
      & fill_array = water%wsn_soil,                  &
      & patch=model%grid%patch)
    water%wsn_soil(:,:) = MERGE(water%wsn_soil(:,:), 0._wp, water%wsn_soil(:,:) >= 0._wp)

    ! Initial soil moisture
    return_pointer => netcdf_read_oncells_2d(         &
      & filename=model%config%soil%ic_filename,       &
      & variable_name='init_moist_layered',           &
      & fill_array = water%ws,                        &
      & patch=model%grid%patch)
    water%ws(:,:) = MERGE(water%ws(:,:), 0._wp, water%ws(:,:) >= 0._wp)

    DO isoil=1,model%config%soil%nsoil_water
      return_pointer3d => netcdf_read_oncells_2d_extdim(  &
        & filename=model%config%soil%ic_filename,       &
        & variable_name='layer_moist',                  &
!!$        & fill_array = water%ws_l,                      &
        & patch=model%grid%patch,                       &
        & start_extdim=isoil, end_extdim=isoil, extdim_name='lev')
      water%ws_l(:,isoil,:) = MERGE(return_pointer3d(:,:,1), 0._wp, return_pointer3d(:,:,1) >= 0._wp)
    END DO
!!$      water%ws_l(:,:,:) = MERGE(water%ws_l(:,:,:), 0._wp, water%ws_l(:,:,:) >= 0._wp)

!!$    ! Temporary: Distribute water content of whole column to layers. Should be read in from ini file.
!!$    DO isoil=1,model%config%soil%nsoil_water
!!$      WHERE (phys%soil_depth(:,:) > 0._wp)
!!$        water%ws_l(:,isoil,:) = water%ws(:,:) * phys%soil_depth_l(:,isoil,:) / phys%soil_depth(:,:)
!!$      ELSEWHERE
!!$        water%ws_l(:,isoil,:) = 0._wp
!!$      END WHERE
!!$    END DO

    ! Water content at field capacity
    !
    ! Read volumetric field capacity from boundary condition file
    vfc => netcdf_read_oncells_2d(                         &
      & filename=model%config%soil%bc_filename,            &
      & variable_name='fieldcap',                          &
      & patch=model%grid%patch)
    vfc(:,:) = MERGE(vfc(:,:), 0._wp, vfc(:,:) >= 0._wp)

    ! Get field capacity by multiplying with soil layer thicknesses
    DO i=1,model%config%soil%nsoil_water
      water%ws_fc_l(:,i,:) = vfc(:,:) * phys%soil_depth_l(:,i,:)
    END DO

    DEALLOCATE(vfc)

    ! Water content at permanent wilting point
    !
    ! Read volumetric wilting point from boundary condition file
    vpwp => netcdf_read_oncells_2d(                         &
      & filename=model%config%soil%bc_filename,             &
      & variable_name='wiltpoint',                          &
      & patch=model%grid%patch)
    vpwp(:,:) = MERGE(vpwp(:,:), 0._wp, vpwp(:,:) >= 0._wp)

    ! Get permanent wilting point by multiplying with soil layer thicknesses
    DO i=1,model%config%soil%nsoil_water
      water%ws_pwp_l(:,i,:) = vpwp(:,:) * phys%soil_depth_l(:,i,:)
    END DO

    DEALLOCATE(vpwp)

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

#endif
END MODULE mo_soil_init
