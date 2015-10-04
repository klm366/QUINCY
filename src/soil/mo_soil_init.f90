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

#include "dsl_jsbach.inc"

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_close_file, jsb_netcdf_open_input, &
                                  jsb_netcdf_read_2d, jsb_netcdf_read_2d_time, jsb_netcdf_read_2d_extdim

  dsl$jsb Use_processes SOIL, SRF
  dsl$jsb Use_config(soil)
  dsl$jsb Use_config(srf)
  dsl$jsb Use_memory(soil)
  dsl$jsb Use_memory(srf)

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
    dsl$jsb Def_memory(soil_mem, soil)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get soil memory of the tile
    dsl$jsb Get_memory(soil_mem, SOIL)

    IF (.NOT. ASSOCIATED(soil_mem)) RETURN

    CALL message(TRIM(routine), 'Setting  soil boundary conditions for tile '//TRIM(tile%name))

  END SUBROUTINE soil_init_bc

  RECURSIVE SUBROUTINE soil_init_ic(tile)

    USE mo_soil_process, ONLY: init_soil_temperature
    USE mo_soil_util,    ONLY: soil_depth_to_layers_2d

    TYPE(t_jsb_tile), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    dsl$jsb Def_config(soil_conf, soil)
    dsl$jsb Def_memory(soil_mem, soil)
    dsl$jsb Def_config(srf_conf, srf)
    dsl$jsb Def_memory(srf_mem, srf)

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
    dsl$jsb Get_memory(soil_mem, SOIL)
    dsl$jsb Get_memory(srf_mem, SRF)

    ! Get soil config
    dsl$jsb Get_config(soil_conf, SOIL)
    dsl$jsb Get_config(srf_conf, SRF)

    IF (.NOT. soil_conf%active) RETURN

    CALL message(TRIM(routine), 'Initializing soil memory for tile '//TRIM(tile%name))

    ! Initialize physical state
    !

    input_file = jsb_netcdf_open_input(TRIM(soil_conf%bc_filename), model%id)

    ! Total soil depth until bedrock

    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='soil_depth',               &
      & fill_array = soil_mem%soil_depth%ptr)
    soil_mem%soil_depth%ptr(:,:) = MAX(0.1_wp, MERGE(soil_mem%soil_depth%ptr(:,:), 0._wp, soil_mem%soil_depth%ptr(:,:) >= 0._wp))

    ! Compute actual soil layer depths based on the fixed layer thicknesses from namelist and
    ! the bedrock depth just read in.
    soil_mem%soil_depth_l%ptr(:,:,:) = soil_depth_to_layers_2d( &
      & soil_mem%soil_depth%ptr(:,:),                         & ! Total soil depth until bedrock (from textures)
      & soil_conf%dz_water(1:soil_conf%nsoil_water))    ! Soil layer thicknesses from namelist

    ! Volumetric soil porosity
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & fill_array = soil_mem%vpor%ptr,                   &
      & variable_name='soil_porosity')
    soil_mem%vpor%ptr(:,:) = MERGE(soil_mem%vpor%ptr(:,:), 0.2_wp, soil_mem%vpor%ptr(:,:) >= 0._wp)

    ! Saturated hydraulic conductivity
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='hyd_cond_sat',             &
      & fill_array = soil_mem%shc%ptr)
    soil_mem%shc%ptr(:,:) = MERGE(soil_mem%shc%ptr(:,:), 4.e-6_wp, soil_mem%shc%ptr(:,:) >= 0._wp)

    ! Matrix potential
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='moisture_pot',             &
      & fill_array = soil_mem%mpot%ptr)
    soil_mem%mpot%ptr(:,:) = MERGE(soil_mem%mpot%ptr(:,:), 0.2_wp, soil_mem%mpot%ptr(:,:) >= 0._wp)

    ! Exponent B in Clapp and Hornberger
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='bclapp',                   &
      & fill_array = soil_mem%bclapp%ptr)
    soil_mem%bclapp%ptr(:,:) = MERGE(soil_mem%bclapp%ptr(:,:), 6._wp, soil_mem%bclapp%ptr(:,:) >= 0._wp)

    ! Pore size distribution index
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='pore_size_index',          &
      & fill_array = soil_mem%psdi%ptr)
    soil_mem%psdi%ptr(:,:) = MERGE(soil_mem%psdi%ptr(:,:), 0.25_wp, soil_mem%psdi%ptr(:,:) >= 0._wp)

    ! Thermal diffusivity of soil
    fao => jsb_netcdf_read_2d(                    &
      & input_file,                               &
      & variable_name='fao')
    fao(:,:) = MERGE(fao(:,:), 0._wp, fao(:,:) >= 0._wp)
    fao(:,:) = MAX(REAL(LBOUND(fao_thermal_diff,1),wp), MIN(REAL(UBOUND(fao_thermal_diff,1),wp), fao(:,:)))

    DO i=1,SIZE(fao,2)
      soil_mem%hdiff%ptr(:,i) = fao_thermal_diff(NINT(fao(:,i)))
    END DO
!!$    soil_mem%hdiff%ptr(:,:) = RESHAPE(fao_thermal_diff(PACK(NINT(fao(:,:)), fao(:,:) > -1._wp)), /( SIZE(fao(:,:),1),SIZE(fao(:,:),2) )/ )

    ! Volumetric heat capacity of dry soil
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='heat_capacity',            &
      & fill_array = soil_mem%vhcap%ptr)
    soil_mem%vhcap%ptr(:,:) = MERGE(soil_mem%vhcap%ptr(:,:), 2.e6_wp, soil_mem%vhcap%ptr(:,:) >= 0._wp)

    CALL jsb_netcdf_close_file(input_file)                              ! Temporary
    input_file = jsb_netcdf_open_input('bc_land_phys.nc', model%id)     ! Temporary

    ! Albedo
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='albedo_soil_vis',          &
      & fill_array = soil_mem%alb_vis_soil%ptr)
    soil_mem%alb_vis_soil%ptr(:,:) = MERGE(soil_mem%alb_vis_soil%ptr(:,:), 0._wp, soil_mem%alb_vis_soil%ptr(:,:) >= 0._wp)

    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='albedo_soil_nir',          &
      & fill_array = soil_mem%alb_nir_soil%ptr)
    soil_mem%alb_nir_soil%ptr(:,:) = MERGE(soil_mem%alb_nir_soil%ptr(:,:), 0._wp, soil_mem%alb_nir_soil%ptr(:,:) >= 0._wp)

    CALL jsb_netcdf_close_file(input_file)

    ! Initialize energy state
    !

    input_file = jsb_netcdf_open_input(srf_conf%ic_filename, model%id)

    ! Soil temperatures
    !
    ! Read climatological surface temperature
    t_srf_mon => jsb_netcdf_read_2d_time(         &   ! Temporary
      & input_file,                               &
      & variable_name='surf_temp')
    t_srf_mon(:,:,:) = MERGE(t_srf_mon(:,:,:), 280._wp, t_srf_mon(:,:,:) > 0._wp)

    DO i=1,SIZE(soil_mem%soil_depth%ptr(:,:), DIM=2)
      CALL init_soil_temperature( t_srf_mon(:,i,1:12), &
        & soil_mem%t_soil%ptr(:,1:soil_conf%nsoil_energy,i), &
        & srf_mem%t_srf%ptr(:,i))
    END DO
    soil_mem%t_soil_acoef%ptr (:,:,:) = 0._wp
    soil_mem%t_soil_bcoef%ptr (:,:,:) = 1._wp

    soil_mem%hcap_grnd%ptr     (:,:) = 0._wp
    soil_mem%hcap_grnd_old%ptr (:,:) = 0._wp

    ! Initialize water state

    soil_mem%wfract_soil%ptr (:,:) = 0._wp
    soil_mem%wsr_soil%ptr    (:,:) = 0._wp

    soil_mem%sfract_soil%ptr (:,:) = 0._wp

    CALL jsb_netcdf_close_file(input_file)
    input_file = jsb_netcdf_open_input(soil_conf%ic_filename, model%id)

    ! Initial snow depth
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='snow',                     &
      & fill_array = soil_mem%wsn_soil%ptr)
    soil_mem%wsn_soil%ptr(:,:) = MERGE(soil_mem%wsn_soil%ptr(:,:), 0._wp, soil_mem%wsn_soil%ptr(:,:) >= 0._wp)

    ! Initial soil moisture
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='init_moist',               &
      & fill_array = soil_mem%ws%ptr)
    soil_mem%ws%ptr(:,:) = MERGE(soil_mem%ws%ptr(:,:), 0._wp, soil_mem%ws%ptr(:,:) >= 0._wp)

    ! Initial overflown water
    soil_mem%overflow%ptr(:,:) = 0._wp

    DO isoil=1,soil_conf%nsoil_water
      return_pointer3d => jsb_netcdf_read_2d_extdim(  &
        & input_file,                                 &
        & variable_name='layer_moist',                &
!!$        & fill_array = soil_mem%ws_l,                      &
        & start_extdim=isoil, end_extdim=isoil, extdim_name='soillev')
      soil_mem%ws_l%ptr(:,isoil,:) = MERGE(return_pointer3d(:,:,1), 0._wp, return_pointer3d(:,:,1) >= 0._wp)
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
    input_file = jsb_netcdf_open_input(TRIM(soil_conf%bc_filename), model%id)

    ! Water content at field capacity
    !
    ! Read volumetric field capacity from boundary condition file
    vfc => jsb_netcdf_read_2d(                         &
      & input_file,                                    &
      & variable_name='soil_field_cap')
    vfc(:,:) = MERGE(vfc(:,:), 0._wp, vfc(:,:) >= 0._wp)

    ! Get field capacity by multiplying with soil layer thicknesses
    DO i=1,soil_conf%nsoil_water
      soil_mem%ws_fc_l%ptr(:,i,:) = vfc(:,:) * soil_mem%soil_depth_l%ptr(:,i,:)
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
    DO i=1,soil_conf%nsoil_water
      soil_mem%ws_pwp_l%ptr(:,i,:) = vpwp(:,:) * soil_mem%soil_depth_l%ptr(:,i,:)
    END DO

    DEALLOCATE(vpwp)

    ! Water content at saturation (water holding capacity)
    !
    ! Get water holding capacity by multiplying volumetric porosity with soil layer thicknesses
    DO i=1,soil_conf%nsoil_water
      soil_mem%ws_sat_l%ptr(:,i,:) = soil_mem%vpor%ptr(:,:) * soil_mem%soil_depth_l%ptr(:,i,:)
    END DO

    CALL jsb_netcdf_close_file(input_file)

    ! Initialize energy fluxes
    !
    soil_mem%grnd_hflx%ptr     (:,:) = 0._wp
    soil_mem%grnd_hflx_old%ptr (:,:) = 0._wp

    ! Initialize water fluxes
    !
    soil_mem%runoff%ptr        (:,:) = 0._wp
    soil_mem%drainage%ptr      (:,:) = 0._wp

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
