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
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input

  dsl4jsb_Use_processes SOIL, SEB, VEG
  dsl4jsb_Use_config(SOIL)
  dsl4jsb_Use_config(SEB)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(VEG)

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

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_init'

    CALL soil_init_bc(tile)
    CALL soil_init_ic(tile)

  END SUBROUTINE soil_init

  SUBROUTINE soil_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SOIL)
    dsl4jsb_Def_memory(SOIL)

    REAL(wp), POINTER :: &
      & return_pointer (:,  :)    !< temporary pointer

    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get soil config
    dsl4jsb_Get_config(SOIL)

    IF (.NOT. dsl4jsb_Config(SOIL)%active) RETURN

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SOIL)

    CALL message(TRIM(routine), 'Setting  soil boundary conditions for tile '//TRIM(tile%name))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SOIL)%bc_sso_filename), model%grid_id)

    ! Orography
    return_pointer => input_file%Read_2d(         &
!!$      & input_file,                               &
      & variable_name='orostd',                   &
      & fill_array = dsl4jsb_var2D_ptr(SOIL,oro_stddev))
    dsl4jsb_var2D_onDomain(SOIL,oro_stddev) = &
      & MERGE(dsl4jsb_var2D_onDomain(SOIL,oro_stddev), 0._wp, dsl4jsb_var2D_onDomain(SOIL,oro_stddev) >= 0._wp)

    CALL input_file%Close()

  END SUBROUTINE soil_init_bc

  SUBROUTINE soil_init_ic(tile)

    USE mo_soil_process, ONLY: init_soil_temperature
    USE mo_soil_util,    ONLY: soil_depth_to_layers_2d

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    dsl4jsb_Def_config(SOIL)
    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_config(SEB)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(VEG)

    INTEGER :: isoil
    REAL(wp), POINTER :: &
      & t_srf_mon      (:,:,:),    &
      & vpor           (:,  :),    & !< Volumetric soil porosity following Cosby et al. [m/m]
      & fao            (:,  :),    & !< FAO soil types
      & return_pointer (:,  :),    & !< temporary pointer
      & return_pointer3d(:,:,:)

    TYPE(t_input_file) :: input_file
    INTEGER :: i, n1, n2

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(SEB)
    IF (tile%Is_process_active(VEG)) THEN
      dsl4jsb_Get_memory(VEG)
    END IF

    ! Get soil config
    dsl4jsb_Get_config(SOIL)
    dsl4jsb_Get_config(SEB)

    IF (.NOT. dsl4jsb_Config(SOIL)%active) RETURN

    CALL message(TRIM(routine), 'Initializing soil memory for tile '//TRIM(tile%name))

    ! Initialize physical state
    !

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SOIL)%bc_filename), model%grid_id)

    ! Total soil depth until bedrock

    return_pointer => input_file%Read_2d(         &
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
    return_pointer => input_file%Read_2d(         &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, vpor),  &
      & variable_name='soil_porosity')
    dsl4jsb_var2D_onDomain(SOIL, vpor) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL, vpor), 0.2_wp, dsl4jsb_var2D_onDomain(SOIL, vpor) >= 0._wp)

    ! Saturated hydraulic conductivity
    return_pointer => input_file%Read_2d(         &
      & variable_name='hyd_cond_sat',             &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, shc))
    dsl4jsb_var2D_onDomain(SOIL, shc) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL, shc), 4.e-6_wp, dsl4jsb_var2D_onDomain(SOIL, shc) >= 0._wp)

    ! Matrix potential
    return_pointer => input_file%Read_2d(         &
      & variable_name='moisture_pot',             &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, mpot))
    dsl4jsb_var2D_onDomain(SOIL, mpot) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL, mpot), 0.2_wp, dsl4jsb_var2D_onDomain(SOIL, mpot) >= 0._wp)

    ! Exponent B in Clapp and Hornberger
    return_pointer => input_file%Read_2d(         &
      & variable_name='bclapp',                   &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, bclapp))
    dsl4jsb_var2D_onDomain(SOIL, bclapp) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL, bclapp), 6._wp, dsl4jsb_var2D_onDomain(SOIL, bclapp) >= 0._wp)

    ! Pore size distribution index
    return_pointer => input_file%Read_2d(         &
      & variable_name='pore_size_index',          &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, psdi))
    dsl4jsb_var2D_onDomain(SOIL, psdi) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL, psdi), 0.25_wp, dsl4jsb_var2D_onDomain(SOIL, psdi) >= 0._wp)

    ! Thermal diffusivity of soil
    fao => input_file%Read_2d(                    &
      & variable_name='fao')
    fao(:,:) = MERGE(fao(:,:), 0._wp, fao(:,:) >= 0._wp)
    fao(:,:) = MAX(REAL(LBOUND(fao_thermal_diff,1),wp), MIN(REAL(UBOUND(fao_thermal_diff,1),wp), fao(:,:)))

    DO i=1,SIZE(fao,2)
      dsl4jsb_var2D_ptr(SOIL,hdiff)(:,i) = fao_thermal_diff(NINT(fao(:,i)))
    END DO
!!$    soil_mem%hdiff%ptr(:,:) = RESHAPE(fao_thermal_diff(PACK(NINT(fao(:,:)), fao(:,:) > -1._wp)), /( SIZE(fao(:,:),1),SIZE(fao(:,:),2) )/ )

    ! Volumetric heat capacity of dry soil
    return_pointer => input_file%Read_2d(         &
      & variable_name='heat_capacity',            &
      & fill_array = dsl4jsb_var2D_ptr(SOIL, vhcap))
    dsl4jsb_var2D_onDomain(SOIL, vhcap) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL, vhcap), 2.e6_wp, dsl4jsb_var2D_onDomain(SOIL, vhcap) >= 0._wp)

    ! Root depth
    IF (tile%Is_process_active(VEG)) THEN
      return_pointer => input_file%Read_2d(         &
        & variable_name='root_depth',               &
        & fill_array = dsl4jsb_var2D_ptr(VEG, root_depth))
      dsl4jsb_var2D_onDomain(VEG, root_depth) = &
        MERGE(dsl4jsb_var2D_onDomain(VEG, root_depth), 0._wp, dsl4jsb_var2D_onDomain(VEG, root_depth) >= 0._wp)
      dsl4jsb_var2D_onDomain(VEG, root_depth) = &           ! Temporary fix (root depth sometimes larger then soil depth)
        MIN(dsl4jsb_var2D_onDomain(VEG, root_depth), dsl4jsb_var2D_onDomain(SOIL, soil_depth))
!!$    dsl4jsb_var2D_onDomain(VEG, root_depth) = soil_mem%soil_depth%ptr(:,:) * 0.75_wp  ! Temporary
      dsl4jsb_var3D_onDomain(VEG, root_depth_l) = soil_depth_to_layers_2d(               &
        & dsl4jsb_var2D_onDomain(VEG, root_depth),                                         & ! Total rooting depth
        & dsl4jsb_Config(SOIL)%dz_water(1:dsl4jsb_Config(SOIL)%nsoil_water))    ! Soil layer thicknesses from namelist
    END IF

    CALL input_file%Close()                              ! Temporary

    ! Initialize energy state
    !

    input_file = jsb_netcdf_open_input(dsl4jsb_Config(SEB)%ic_filename, model%grid_id)

    ! Soil temperatures
    !
    ! Read climatological surface temperature
    t_srf_mon => input_file%Read_2d_time(         &   ! Temporary
      & variable_name='surf_temp')
    t_srf_mon(:,:,:) = MERGE(t_srf_mon(:,:,:), 280._wp, t_srf_mon(:,:,:) > 0._wp)

    DO i=1,SIZE(dsl4jsb_var2D_ptr(SOIL,soil_depth), DIM=2)
      CALL init_soil_temperature( t_srf_mon(:,i,1:12), &
        & dsl4jsb_var2D_ptr(SOIL,t_soil)(:,1:dsl4jsb_Config(SOIL)%nsoil_energy,i), &
        & dsl4jsb_var2D_ptr(SEB,t)(:,i))
    END DO
! R    dsl4jsb_var3D_onDomain(SOIL, t_soil_acoef) = 0._wp
! R    dsl4jsb_var3D_onDomain(SOIL, t_soil_bcoef) = 1._wp

! R    dsl4jsb_var2D_onDomain(SOIL, hcap_grnd)     = 0._wp
! R    dsl4jsb_var2D_onDomain(SOIL, hcap_grnd_old) = 0._wp


    CALL input_file%Close()


    ! Initialize energy fluxes
    !
! R    dsl4jsb_var2D_onDomain(SOIL,grnd_hflx)     = 0._wp
! R    dsl4jsb_var2D_onDomain(SOIL,grnd_hflx_old) = 0._wp


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
