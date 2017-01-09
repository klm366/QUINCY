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

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_close_file, jsb_netcdf_open_input, &
                                  jsb_netcdf_read_2d, jsb_netcdf_read_2d_time, jsb_netcdf_read_2d_extdim

  dsl4jsb_Use_processes SOIL, SRF
  dsl4jsb_Use_config(SOIL)
  dsl4jsb_Use_config(SRF)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SRF)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: soil_init

  REAL(wp), PARAMETER :: &
    fao_thermal_diff(0:6) = (/7.4e-7_wp ,8.7e-7_wp ,8.0e-7_wp ,7.4e-7_wp ,7.1e-7_wp ,6.7e-7_wp ,-1._wp/)

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_init'

CONTAINS

  !
  !> Intialize soil process (after memory has been set up)
  !
  SUBROUTINE soil_init(tile)

    TYPE(t_jsb_tile), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_init'

    CALL soil_init_bc(tile)
    CALL soil_init_ic(tile)

  END SUBROUTINE soil_init

  SUBROUTINE soil_init_bc(tile)

    TYPE(t_jsb_tile), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model
    dsl4jsb_Def_memory(SOIL)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SOIL)

    IF (.NOT. ASSOCIATED(soil__mem)) RETURN

    CALL message(TRIM(routine), 'Setting  soil boundary conditions for tile '//TRIM(tile%name))

  END SUBROUTINE soil_init_bc

  RECURSIVE SUBROUTINE soil_init_ic(tile)

    USE mo_soil_process, ONLY: init_soil_temperature
    USE mo_soil_util,    ONLY: soil_depth_to_layers_2d

    TYPE(t_jsb_tile), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    dsl4jsb_Def_config(SOIL)
    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_config(SRF)
    dsl4jsb_Def_memory(SRF)

    INTEGER :: isoil
    REAL(wp), POINTER :: &
      & t_srf_mon      (:,:,:),    &
      & vpwp           (:,  :),    & !< Volumetric soil permanent wilting point [m/m]
      & vpor           (:,  :),    & !< Volumetric soil porosity following Cosby et al. [m/m]
      & vfc            (:,  :),    & !< Volumetric soil field capacity [m/m]
      & fao            (:,  :),    & !< FAO soil types
      & return_pointer (:,  :),    & !< temporary pointer
      & return_pointer3d(:,:,:)

    TYPE(t_input_file) :: input_file
    INTEGER :: i, n1, n2

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(SRF)

    ! Get soil config
    dsl4jsb_Get_config(SOIL)
    dsl4jsb_Get_config(SRF)

    IF (.NOT. dsl4jsb_Config(SOIL)%active) RETURN

    CALL message(TRIM(routine), 'Initializing soil memory for tile '//TRIM(tile%name))

    ! Initialize physical state
    !

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SOIL)%bc_filename), model%id)

    ! Total soil depth until bedrock

    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='soil_depth',               &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, soil_depth))
    dsl4jsb_var2D_onDomain(SOIL, soil_depth) = &
      MAX(0.1_wp, MERGE(dsl4jsb_var2D_onDomain(SOIL, soil_depth), 0._wp, dsl4jsb_var2D_onDomain(SOIL, soil_depth) >= 0._wp))

    ! Compute actual soil layer depths based on the fixed layer thicknesses from namelist and
    ! the bedrock depth just read in.
    dsl4jsb_var3D_onDomain(SOIL, soil_depth_l) = soil_depth_to_layers_2d( &
      & dsl4jsb_var2D_onDomain(SOIL, soil_depth),                         & ! Total soil depth until bedrock (from textures)
      & dsl4jsb_Config(SOIL)%dz_water(1:dsl4jsb_Config(SOIL)%nsoil_water))    ! Soil layer thicknesses from namelist

    ! Volumetric soil porosity
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, vpor),  &
      & variable_name='soil_porosity')
    dsl4jsb_var2D_onDomain(SOIL, vpor) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL, vpor), 0.2_wp, dsl4jsb_var2D_onDomain(SOIL, vpor) >= 0._wp)

    ! Saturated hydraulic conductivity
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='hyd_cond_sat',             &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, shc))
    dsl4jsb_var2D_onDomain(SOIL, shc) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL, shc), 4.e-6_wp, dsl4jsb_var2D_onDomain(SOIL, shc) >= 0._wp)

    ! Matrix potential
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='moisture_pot',             &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, mpot))
    dsl4jsb_var2D_onDomain(SOIL, mpot) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL, mpot), 0.2_wp, dsl4jsb_var2D_onDomain(SOIL, mpot) >= 0._wp)

    ! Exponent B in Clapp and Hornberger
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='bclapp',                   &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, bclapp))
    dsl4jsb_var2D_onDomain(SOIL, bclapp) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL, bclapp), 6._wp, dsl4jsb_var2D_onDomain(SOIL, bclapp) >= 0._wp)

    ! Pore size distribution index
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='pore_size_index',          &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, psdi))
    dsl4jsb_var2D_onDomain(SOIL, psdi) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL, psdi), 0.25_wp, dsl4jsb_var2D_onDomain(SOIL, psdi) >= 0._wp)

    ! Thermal diffusivity of soil
    fao => jsb_netcdf_read_2d(                    &
      & input_file,                               &
      & variable_name='fao')
    fao(:,:) = MERGE(fao(:,:), 0._wp, fao(:,:) >= 0._wp)
    fao(:,:) = MAX(REAL(LBOUND(fao_thermal_diff,1),wp), MIN(REAL(UBOUND(fao_thermal_diff,1),wp), fao(:,:)))

    DO i=1,SIZE(fao,2)
      dsl4jsb_var2D_ptr(SOIL,hdiff)(:,i) = fao_thermal_diff(NINT(fao(:,i)))
    END DO
!!$    soil_mem%hdiff%ptr(:,:) = RESHAPE(fao_thermal_diff(PACK(NINT(fao(:,:)), fao(:,:) > -1._wp)), /( SIZE(fao(:,:),1),SIZE(fao(:,:),2) )/ )

    ! Volumetric heat capacity of dry soil
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='heat_capacity',            &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, vhcap))
    dsl4jsb_var2D_onDomain(SOIL, vhcap) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL, vhcap), 2.e6_wp, dsl4jsb_var2D_onDomain(SOIL, vhcap) >= 0._wp)

    CALL jsb_netcdf_close_file(input_file)                              ! Temporary
    input_file = jsb_netcdf_open_input('bc_land_phys.nc', model%id)     ! Temporary

    ! Albedo
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='albedo_soil_vis',          &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, alb_vis_soil))
    dsl4jsb_var2D_onDomain(SOIL, alb_vis_soil) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL,alb_vis_soil), 0._wp, dsl4jsb_var2D_onDomain(SOIL,alb_vis_soil) >= 0._wp)

    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='albedo_soil_nir',          &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, alb_nir_soil))
    dsl4jsb_var2D_onDomain(SOIL,alb_nir_soil) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL,alb_nir_soil), 0._wp, dsl4jsb_var2D_onDomain(SOIL,alb_nir_soil) >= 0._wp)

    CALL jsb_netcdf_close_file(input_file)

    ! Initialize energy state
    !

    input_file = jsb_netcdf_open_input(dsl4jsb_Config(SRF)%ic_filename, model%id)

    ! Soil temperatures
    !
    ! Read climatological surface temperature
    t_srf_mon => jsb_netcdf_read_2d_time(         &   ! Temporary
      & input_file,                               &
      & variable_name='surf_temp')
    t_srf_mon(:,:,:) = MERGE(t_srf_mon(:,:,:), 280._wp, t_srf_mon(:,:,:) > 0._wp)

    DO i=1,SIZE(dsl4jsb_var2D_ptr(SOIL,soil_depth), DIM=2)
      CALL init_soil_temperature( t_srf_mon(:,i,1:12), &
        & dsl4jsb_var2D_ptr(SOIL,t_soil)(:,1:dsl4jsb_Config(SOIL)%nsoil_energy,i), &
        & dsl4jsb_var2D_ptr(SRF,t_srf)(:,i))
    END DO
    dsl4jsb_var3D_onDomain(SOIL, t_soil_acoef) = 0._wp
    dsl4jsb_var3D_onDomain(SOIL, t_soil_bcoef) = 1._wp

    dsl4jsb_var2D_onDomain(SOIL, hcap_grnd)     = 0._wp
    dsl4jsb_var2D_onDomain(SOIL, hcap_grnd_old) = 0._wp

    ! Initialize water state

    dsl4jsb_var2D_onDomain(SOIL, wfract_soil) = 0._wp
    dsl4jsb_var2D_onDomain(SOIL, wsr_soil)    = 0._wp

    dsl4jsb_var2D_onDomain(SOIL, sfract_soil) = 0._wp

    CALL jsb_netcdf_close_file(input_file)
    input_file = jsb_netcdf_open_input(dsl4jsb_Config(SOIL)%ic_filename, model%id)

    ! Initial snow depth
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='snow',                     &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, wsn_soil))
    dsl4jsb_var2D_onDomain(SOIL, wsn_soil) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL,wsn_soil), 0._wp, dsl4jsb_var2D_onDomain(SOIL,wsn_soil) >= 0._wp)

    ! Initial soil moisture
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='init_moist',               &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, ws))
    dsl4jsb_var2D_onDomain(SOIL,ws) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL,ws), 0._wp, dsl4jsb_var2D_onDomain(SOIL,ws) >= 0._wp)

    ! Initial overflown water
    dsl4jsb_var2D_onDomain(SOIL,overflow) = 0._wp

    DO isoil=1,dsl4jsb_Config(SOIL)%nsoil_water
      return_pointer3d => jsb_netcdf_read_2d_extdim(  &
        & input_file,                                 &
        & variable_name='layer_moist',                &
!!$        & fill_array = soil_mem%ws_l,                      &
        & start_extdim=isoil, end_extdim=isoil, extdim_name='soillev')
      dsl4jsb_var3D_ptr(SOIL,ws_l) (:,isoil,:) = &
        MERGE(return_pointer3d(:,:,1), 0._wp, return_pointer3d(:,:,1) >= 0._wp)
    END DO
!!$      soil_mem%ws_l%ptr(:,:,:) = MERGE(soil_mem%ws_l%ptr(:,:,:), 0._wp, soil_mem%ws_l%ptr(:,:,:) >= 0._wp)

!!$    ! Temporary: Distribute water content of whole column to layers. Should be read in from ini file.
!!$    DO isoil=1,soil_conf%nsoil_water
!!$      WHERE (soil_mem%soil_depth%ptr(:,:) > 0._wp)
!!$        soil_mem%ws_l%ptr(:,isoil,:) = soil_mem%ws%ptr(:,:) * soil_mem%soil_depth_l%ptr(:,isoil,:) / soil_mem%soil_depth%ptr(:,:)
!!$      ELSEWHERE
!!$        soil_mem%ws_l%ptr(:,isoil,:) = 0._wp
!!$      END WHERE
!!$    END DO

    CALL jsb_netcdf_close_file(input_file)
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SOIL)%bc_filename), model%id)

    ! Water content at field capacity
    !
    ! Read volumetric field capacity from boundary condition file
    vfc => jsb_netcdf_read_2d(                         &
      & input_file,                                    &
      & variable_name='soil_field_cap')
    vfc(:,:) = MERGE(vfc(:,:), 0._wp, vfc(:,:) >= 0._wp)

    ! Get field capacity by multiplying with soil layer thicknesses
    DO i=1,dsl4jsb_Config(SOIL)%nsoil_water
      dsl4jsb_var3D_ptr(SOIL,ws_fc_l) (:,i,:) = vfc(:,:) * dsl4jsb_var3D_ptr(SOIL,soil_depth_l) (:,i,:)
    END DO

    DEALLOCATE(vfc)

    ! Water content at permanent wilting point
    !
    ! Read volumetric wilting point from boundary condition file
    vpwp => jsb_netcdf_read_2d(                         &
      & input_file,                                     &
      & variable_name='wilting_point')
    vpwp(:,:) = MERGE(vpwp(:,:), 0._wp, vpwp(:,:) >= 0._wp)

    ! Get permanent wilting point by multiplying with soil layer thicknesses
    DO i=1,dsl4jsb_Config(SOIL)%nsoil_water
      dsl4jsb_var3D_ptr(SOIL,ws_pwp_l) (:,i,:) = vpwp(:,:) * dsl4jsb_var3D_ptr(SOIL,soil_depth_l) (:,i,:)
    END DO

    DEALLOCATE(vpwp)

    ! Water content at saturation (water holding capacity)
    !
    ! Get water holding capacity by multiplying volumetric porosity with soil layer thicknesses
    DO i=1,dsl4jsb_Config(SOIL)%nsoil_water
      dsl4jsb_var3D_ptr(SOIL,ws_sat_l) (:,i,:) = &
        dsl4jsb_var2D_onDomain(SOIL,vpor) * dsl4jsb_var3D_ptr(SOIL,soil_depth_l) (:,i,:)
    END DO

    CALL jsb_netcdf_close_file(input_file)

    ! Initialize energy fluxes
    !
    dsl4jsb_var2D_onDomain(SOIL,grnd_hflx)     = 0._wp
    dsl4jsb_var2D_onDomain(SOIL,grnd_hflx_old) = 0._wp

    ! Initialize water fluxes
    !
    dsl4jsb_var2D_onDomain(SOIL,runoff)   = 0._wp
    dsl4jsb_var2D_onDomain(SOIL,drainage) = 0._wp

    ! Propagate initial conditions down the tile structure, starting with the tile with which this subroutine was called initially
    ! TBD: sub-tiles should not read the netcdf files again
!!$    IF (tile%no_of_sub_tiles > 0) THEN
!!$      DO i=1,tile%no_of_sub_tiles
!!$        CALL soil_init_ic(tile%sub_tiles(i))
!!$      END DO
!!$    END IF

  END SUBROUTINE soil_init_ic

#endif
END MODULE mo_soil_init
