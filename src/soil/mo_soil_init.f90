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
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input

  dsl4jsb_Use_processes SOIL_, SEB_, VEG_
  dsl4jsb_Use_config(SOIL_)
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_memory(SOIL_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(VEG_)

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

    dsl4jsb_Def_config(SOIL_)
    dsl4jsb_Def_memory(SOIL_)

    REAL(wp), POINTER :: &
      & return_pointer (:,  :)    !< temporary pointer

    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get soil config
    dsl4jsb_Get_config(SOIL_)

    IF (.NOT. dsl4jsb_Config(SOIL_)%active) RETURN

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SOIL_)

    CALL message(TRIM(routine), 'Setting  soil boundary conditions for tile '//TRIM(tile%name))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SOIL_)%bc_sso_filename), model%grid_id)

    ! Orography ...
    return_pointer => input_file%Read_2d(         &
      & variable_name='elevation',                &
      & fill_array = dsl4jsb_var2D_ptr(SOIL_,elevation))
    dsl4jsb_var2D_onDomain(SOIL_,elevation) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL_,elevation), 0._wp, dsl4jsb_var2D_onDomain(SOIL_,elevation) >= 0._wp)

    return_pointer => input_file%Read_2d(         &
      & variable_name='orostd',                   &
      & fill_array = dsl4jsb_var2D_ptr(SOIL_,oro_stddev))
    dsl4jsb_var2D_onDomain(SOIL_,oro_stddev) = &
      & MERGE(dsl4jsb_var2D_onDomain(SOIL_,oro_stddev), 0._wp, dsl4jsb_var2D_onDomain(SOIL_,oro_stddev) >= 0._wp)

    CALL input_file%Close()

  END SUBROUTINE soil_init_bc

  SUBROUTINE soil_init_ic(tile)

    USE mo_soil_process, ONLY: init_soil_temperature
    USE mo_soil_util,    ONLY: soil_depth_to_layers_2d

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    dsl4jsb_Def_config(SOIL_)
    dsl4jsb_Def_memory(SOIL_)
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(VEG_)

    INTEGER :: isoil
    REAL(wp), POINTER :: &
      & t_srf_mon      (:,:,:),    &
      & vpor           (:,  :),    & !< Volumetric soil porosity following Cosby et al. [m/m]
      & fao            (:,  :),    & !< FAO soil types
      & return_pointer (:,  :),    & !< temporary pointer
      & return_pointer3d(:,:,:)

    TYPE(t_input_file) :: input_file
    INTEGER :: iblk, n1, n2

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(SEB_)
    IF (tile%Is_process_active(VEG_)) THEN
      dsl4jsb_Get_memory(VEG_)
    END IF

    ! Get soil config
    dsl4jsb_Get_config(SOIL_)
    dsl4jsb_Get_config(SEB_)

    IF (.NOT. dsl4jsb_Config(SOIL_)%active) RETURN

    CALL message(TRIM(routine), 'Initializing soil memory for tile '//TRIM(tile%name))

    ! Initialize physical state
    !

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SOIL_)%bc_filename), model%grid_id)

    ! Total soil depth until bedrock

    return_pointer => input_file%Read_2d(         &
      & variable_name='soil_depth',               &
      & fill_array = dsl4jsb_var2D_ptr(SOIL_, soil_depth))
    dsl4jsb_var2D_onDomain(SOIL_, soil_depth) = &
      MAX(0.1_wp, MERGE(dsl4jsb_var2D_onDomain(SOIL_, soil_depth), 0._wp, dsl4jsb_var2D_onDomain(SOIL_, soil_depth) >= 0._wp))

    ! Compute actual soil layer depths based on the fixed layer thicknesses from namelist and
    ! the bedrock depth just read in.
    dsl4jsb_var3D_onDomain(SOIL_, soil_depth_l) = soil_depth_to_layers_2d( &
      & dsl4jsb_var2D_onDomain(SOIL_, soil_depth),                         & ! Total soil depth until bedrock (from textures)
      & dsl4jsb_Config(SOIL_)%dz_water(1:dsl4jsb_Config(SOIL_)%nsoil_water))    ! Soil layer thicknesses from namelist

    ! Volumetric soil porosity
    return_pointer => input_file%Read_2d(         &
      & fill_array = dsl4jsb_var2D_ptr(SOIL_, vpor),  &
      & variable_name='soil_porosity')
    dsl4jsb_var2D_onDomain(SOIL_, vpor) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL_, vpor), 0.2_wp, dsl4jsb_var2D_onDomain(SOIL_, vpor) >= 0._wp)

    ! Saturated hydraulic conductivity
    return_pointer => input_file%Read_2d(         &
      & variable_name='hyd_cond_sat',             &
      & fill_array = dsl4jsb_var2D_ptr(SOIL_, shc))
    dsl4jsb_var2D_onDomain(SOIL_, shc) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL_, shc), 4.e-6_wp, dsl4jsb_var2D_onDomain(SOIL_, shc) >= 0._wp)

    ! Matrix potential
    return_pointer => input_file%Read_2d(         &
      & variable_name='moisture_pot',             &
      & fill_array = dsl4jsb_var2D_ptr(SOIL_, mpot))
    dsl4jsb_var2D_onDomain(SOIL_, mpot) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL_, mpot), 0.2_wp, dsl4jsb_var2D_onDomain(SOIL_, mpot) >= 0._wp)

    ! Exponent B in Clapp and Hornberger
    return_pointer => input_file%Read_2d(         &
      & variable_name='bclapp',                   &
      & fill_array = dsl4jsb_var2D_ptr(SOIL_, bclapp))
    dsl4jsb_var2D_onDomain(SOIL_, bclapp) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL_, bclapp), 6._wp, dsl4jsb_var2D_onDomain(SOIL_, bclapp) >= 0._wp)

    ! Pore size distribution index
    return_pointer => input_file%Read_2d(         &
      & variable_name='pore_size_index',          &
      & fill_array = dsl4jsb_var2D_ptr(SOIL_, psdi))
    dsl4jsb_var2D_onDomain(SOIL_, psdi) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL_, psdi), 0.25_wp, dsl4jsb_var2D_onDomain(SOIL_, psdi) >= 0._wp)

    ! Thermal diffusivity of soil
    fao => input_file%Read_2d(                    &
      & variable_name='fao')
    fao(:,:) = MERGE(fao(:,:), 0._wp, fao(:,:) >= 0._wp)
    fao(:,:) = MAX(REAL(LBOUND(fao_thermal_diff,1),wp), MIN(REAL(UBOUND(fao_thermal_diff,1),wp), fao(:,:)))

    DO iblk=1,SIZE(fao,2)
      dsl4jsb_var2D_ptr(SOIL_,hdiff)(:,iblk) = fao_thermal_diff(NINT(fao(:,iblk)))
    END DO
!!$    soil_mem%hdiff%ptr(:,:) = RESHAPE(fao_thermal_diff(PACK(NINT(fao(:,:)), fao(:,:) > -1._wp)), /( SIZE(fao(:,:),1),SIZE(fao(:,:),2) )/ )

    ! Volumetric heat capacity of dry soil
    return_pointer => input_file%Read_2d(         &
      & variable_name='heat_capacity',            &
      & fill_array = dsl4jsb_var2D_ptr(SOIL_, vhcap))
    dsl4jsb_var2D_onDomain(SOIL_, vhcap) = &
      MERGE(dsl4jsb_var2D_onDomain(SOIL_, vhcap), 2.e6_wp, dsl4jsb_var2D_onDomain(SOIL_, vhcap) > 0._wp)

    ! Root depth
    IF (tile%Is_process_active(VEG_)) THEN
      return_pointer => input_file%Read_2d(         &
        & variable_name='root_depth',               &
        & fill_array = dsl4jsb_var2D_ptr(VEG_, root_depth))
      dsl4jsb_var2D_onDomain(VEG_, root_depth) = &
        MERGE(dsl4jsb_var2D_onDomain(VEG_, root_depth), 0._wp, dsl4jsb_var2D_onDomain(VEG_, root_depth) >= 0._wp)
      dsl4jsb_var2D_onDomain(VEG_, root_depth) = &           ! Temporary fix (root depth sometimes larger then soil depth)
        MIN(dsl4jsb_var2D_onDomain(VEG_, root_depth), dsl4jsb_var2D_onDomain(SOIL_, soil_depth))
!!$    dsl4jsb_var2D_onDomain(VEG_, root_depth) = soil_mem%soil_depth%ptr(:,:) * 0.75_wp  ! Temporary
      dsl4jsb_var3D_onDomain(VEG_, root_depth_l) = soil_depth_to_layers_2d(               &
        & dsl4jsb_var2D_onDomain(VEG_, root_depth),                                         & ! Total rooting depth
        & dsl4jsb_Config(SOIL_)%dz_water(1:dsl4jsb_Config(SOIL_)%nsoil_water))    ! Soil layer thicknesses from namelist
    END IF

    CALL input_file%Close()                              ! Temporary

    ! Initialize energy state
    !

    input_file = jsb_netcdf_open_input(dsl4jsb_Config(SEB_)%ic_filename, model%grid_id)

    ! Soil temperatures
    !
    ! Read climatological surface temperature
    t_srf_mon => input_file%Read_2d_time(         &   ! Temporary
      & variable_name='surf_temp')
    t_srf_mon(:,:,:) = MERGE(t_srf_mon(:,:,:), 280._wp, t_srf_mon(:,:,:) > 0._wp)

    DO iblk=1,SIZE(dsl4jsb_var2D_ptr(SOIL_,soil_depth), DIM=2)
      CALL init_soil_temperature( t_srf_mon(:,iblk,1:12), &
        & dsl4jsb_var2D_ptr(SOIL_,t_soil)(:,1:dsl4jsb_Config(SOIL_)%nsoil_energy,iblk), &
        & dsl4jsb_var2D_ptr(SEB_,t)(:,iblk))
    END DO
! R    dsl4jsb_var3D_onDomain(SOIL_, t_soil_acoef) = 0._wp
! R    dsl4jsb_var3D_onDomain(SOIL_, t_soil_bcoef) = 1._wp

! R    dsl4jsb_var2D_onDomain(SOIL_, hcap_grnd)     = 0._wp
! R    dsl4jsb_var2D_onDomain(SOIL_, hcap_grnd_old) = 0._wp


    CALL input_file%Close()


    ! Initialize energy fluxes
    !
! R    dsl4jsb_var2D_onDomain(SOIL_,grnd_hflx)     = 0._wp
! R    dsl4jsb_var2D_onDomain(SOIL_,grnd_hflx_old) = 0._wp


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
