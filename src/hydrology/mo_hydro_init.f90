!>
!! @brief Initialization of the the hydrology memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-09)
!!
MODULE mo_hydro_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_lct_class,     ONLY: LAND_TYPE, BARE_TYPE, VEG_TYPE, GLACIER_TYPE, LAKE_TYPE, Contains_lct
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input

  dsl4jsb_Use_processes HYDRO_, SOIL_, VEG_
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_config(SOIL_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(SOIL_)
  dsl4jsb_Use_memory(VEG_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: hydro_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_init'

CONTAINS

  !
  !> Intialize hydrology process (after memory has been set up)
  !
  SUBROUTINE hydro_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':hydro_init'

    CALL hydro_init_bc(tile)
    CALL hydro_init_ic(tile)

  END SUBROUTINE hydro_init

  SUBROUTINE hydro_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model

    dsl4jsb_Def_memory(HYDRO_)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':hydro_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get hydrology memory of the tile
    dsl4jsb_Get_memory(HYDRO_)

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    !CALL message(TRIM(routine), 'Setting  hydrology boundary conditions for tile '//TRIM(tile%name))

    ! R: nothing to do in the moment ...

  END SUBROUTINE hydro_init_bc


  SUBROUTINE hydro_init_ic(tile)

    !-----------------------------------------------------------------------
    !  DECLARATIONS
    USE mo_soil_util,   ONLY: soil_depth_to_layers_2d
    USE mo_veg_util,    ONLY: get_water_in_root_zone

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(SOIL_)

    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SOIL_)
    dsl4jsb_Def_memory(VEG_)

    INTEGER :: ihydrology
    REAL(wp), POINTER ::         &
      & return_pointer  (:,  :), & !< temporary pointer
      & return_pointer3d(:,:,:), &
      & prt_2D          (:,  :)

    TYPE(t_input_file) :: input_file
    INTEGER :: i, isoil

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':hydro_init_ic'

    !-----------------------------------------------------------------------
    ! CONTENT

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    CALL message(TRIM(routine), 'Initializing hydrology memory for tile '//TRIM(tile%name))

    model => Get_model(tile%owner_model_id)

    ! Get hydrology config
    dsl4jsb_Get_config(HYDRO_)

    ! Get hydrology memory of the tile
    dsl4jsb_Get_memory(HYDRO_)

    !----------

    IF (Contains_lct(tile%lcts, VEG_TYPE) .OR. Contains_lct(tile%lcts, BARE_TYPE) &
      & .OR. Contains_lct(tile%lcts, LAND_TYPE)) THEN

      dsl4jsb_Get_config(SOIL_)
      dsl4jsb_Get_memory(SOIL_)

      ! Initialize hydrology variables from bc file

      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HYDRO_)%bc_filename), model%grid_id)

      ! Water content at field capacity
      !
      ! Read volumetric field capacity from boundary condition file
      prt_2D => input_file%Read_2d(                         &
       & variable_name='soil_field_cap')  !< Volumetric soil field capacity [m/m]
                                          ! TBD: R: please change name in infile to w_soil_fcvol or another name that indicates "volumetric"
      prt_2D(:,:) = MERGE(prt_2D(:,:), 0._wp, prt_2D(:,:) >= 0._wp)

      ! Get field capacity by multiplying with soil layer thicknesses
      DO i=1,dsl4jsb_Config(SOIL_)%nsoil_water
        dsl4jsb_var3D_ptr(HYDRO_,w_soil_fc_sl) (:,i,:) = prt_2D(:,:) * dsl4jsb_var3D_ptr(SOIL_,soil_depth_l) (:,i,:)
      END DO
      DEALLOCATE(prt_2D)

      ! Water content at permanent wilting point
      !
      ! Read volumetric wilting point from boundary condition file
      prt_2D => input_file%Read_2d(                         &
        & variable_name='wilting_point') !< Volumetric soil permanent wilting point [m/m]
      prt_2D(:,:) = MERGE(prt_2D(:,:), 0._wp, prt_2D(:,:) >= 0._wp)

      ! Get permanent wilting point by multiplying with soil layer thicknesses
      DO i=1,dsl4jsb_Config(SOIL_)%nsoil_water
        dsl4jsb_var3D_ptr(HYDRO_,w_soil_pwp_sl) (:,i,:) = prt_2D(:,:) * dsl4jsb_var3D_ptr(SOIL_,soil_depth_l) (:,i,:)
      END DO
      DEALLOCATE(prt_2D)

      ! Water content at saturation (water holding capacity)
      !
      ! Get water holding capacity by multiplying volumetric porosity with soil layer thicknesses
      DO i=1,dsl4jsb_Config(SOIL_)%nsoil_water
        dsl4jsb_var3D_ptr(HYDRO_,w_soil_sat_sl) (:,i,:) = &
          dsl4jsb_var2D_onDomain(SOIL_,vpor) * dsl4jsb_var3D_ptr(SOIL_,soil_depth_l) (:,i,:)
      END DO

      CALL input_file%Close()

      ! Initialize hydrology variables from ic file

      input_file = jsb_netcdf_open_input(dsl4jsb_Config(HYDRO_)%ic_filename, model%grid_id)

      ! Initial snow depth
      return_pointer => input_file%Read_2d(         &
        & variable_name='snow',                     &
        & fill_array = dsl4jsb_var2D_ptr(HYDRO_, w_snow_soil))
      dsl4jsb_var2D_onDomain(HYDRO_, w_snow_soil) = &
        & MERGE(dsl4jsb_var2D_onDomain(HYDRO_,w_snow_soil), 0._wp, dsl4jsb_var2D_onDomain(HYDRO_,w_snow_soil) >= 0._wp)

      ! Initial soil moisture
      return_pointer => input_file%Read_2d(         &
        & variable_name='init_moist',               &
        & fill_array = dsl4jsb_var2D_ptr(HYDRO_, w_soil_column))
      dsl4jsb_var2D_onDomain(HYDRO_,w_soil_column) = &
        & MERGE(dsl4jsb_var2D_onDomain(HYDRO_,w_soil_column), 0._wp, dsl4jsb_var2D_onDomain(HYDRO_,w_soil_column) >= 0._wp)

      ! Initial w_soil_overflown water
      dsl4jsb_var2D_onDomain(HYDRO_,w_soil_overflow) = 0._wp
      DO isoil=1,dsl4jsb_Config(SOIL_)%nsoil_water
        return_pointer3d => input_file%Read_2d_extdim(  &
          & variable_name='layer_moist',                &
!!$        & fill_array = soil_mem%w_soil_sl,                      &
          & start_extdim=isoil, end_extdim=isoil, extdim_name='soillev')
        dsl4jsb_var3D_ptr(HYDRO_,w_soil_sl) (:,isoil,:) = &
          & MERGE(return_pointer3d(:,:,1), 0._wp, return_pointer3d(:,:,1) >= 0._wp)
      END DO
!!$      soil_mem%w_soil_sl%ptr(:,:,:) = MERGE(soil_mem%w_soil_sl%ptr(:,:,:), 0._wp, soil_mem%w_soil_sl%ptr(:,:,:) >= 0._wp)

!!$    ! Temporary: Distribute water content of whole column to layers. Should be read in from ini file.
!!$    DO isoil=1,soil_conf%nsoil_water
!!$      WHERE (soil_mem%soil_depth%ptr(:,:) > 0._wp)
!!$        soil_mem%w_soil_sl%ptr(:,isoil,:) = soil_mem%ws%ptr(:,:) * soil_mem%soil_depth_l%ptr(:,isoil,:) / soil_mem%soil_depth%ptr(:,:)
!!$      ELSEWHERE
!!$        soil_mem%w_soil_sl%ptr(:,isoil,:) = 0._wp
!!$      END WHERE
!!$    END DO
      CALL input_file%Close()

    END IF

    IF (Contains_lct(tile%lcts, VEG_TYPE)) THEN
      ! Calculate water content and field capacity in the root zone
      dsl4jsb_Get_memory(VEG_)
      DO i=1,SIZE(dsl4jsb_var2D_ptr(HYDRO_, w_soil_root), 2)
        dsl4jsb_var2D_ptr(HYDRO_, w_soil_root)       (:,  i)    = get_water_in_root_zone( &
          & dsl4jsb_var3D_ptr(HYDRO_, w_soil_sl)     (:,:,i), &
          & dsl4jsb_var3D_ptr(SOIL_, soil_depth_l)   (:,:,i), &
          & dsl4jsb_var3D_ptr(VEG_, root_depth_l)    (:,:,i))
        dsl4jsb_var2D_ptr(HYDRO_, w_soil_root_fc)    (:,  i)    = get_water_in_root_zone( &
          & dsl4jsb_var3D_ptr(HYDRO_, w_soil_fc_sl)  (:,:,i), &
          & dsl4jsb_var3D_ptr(SOIL_, soil_depth_l)   (:,:,i), &
          & dsl4jsb_var3D_ptr(VEG_, root_depth_l)    (:,:,i))
        dsl4jsb_var2D_ptr(HYDRO_, w_soil_root_pwp)   (:,  i)    = get_water_in_root_zone( &
          & dsl4jsb_var3D_ptr(HYDRO_, w_soil_pwp_sl) (:,:,i), &
          & dsl4jsb_var3D_ptr(SOIL_, soil_depth_l)   (:,:,i), &
          & dsl4jsb_var3D_ptr(VEG_, root_depth_l)    (:,:,i))
      END DO
    END IF

  END SUBROUTINE hydro_init_ic

#endif
END MODULE mo_hydro_init
