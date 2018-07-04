!>
!! @brief Initialization of the the radiation memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-09)
!!
MODULE mo_rad_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,  ONLY: t_jsb_config_p
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input

  dsl4jsb_Use_processes RAD_
  dsl4jsb_Use_config(RAD_)
  dsl4jsb_Use_memory(RAD_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: rad_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_init'

CONTAINS

  !
  !> Intialize radiation process (after memory has been set up)
  !
  SUBROUTINE rad_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':rad_init'

    CALL rad_init_bc(tile)
    CALL rad_init_ic(tile)

  END SUBROUTINE rad_init

  SUBROUTINE rad_init_bc(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    TYPE(t_jsb_model),       POINTER :: model
    dsl4jsb_Def_memory(RAD_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':rad_init_bc'

    !-----------------------------------------------------------------------
    ! CONTENT

    model => get_model(tile%owner_model_id)

    ! Get radiation memory of the tile
    dsl4jsb_Get_memory(RAD_)

    IF (.NOT. ASSOCIATED(rad__mem)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting  radiation boundary conditions for tile '//TRIM(tile%name))

    ! R: nothing to do in the moment ...

  END SUBROUTINE rad_init_bc

  RECURSIVE SUBROUTINE rad_init_ic(tile)

!    USE mo_rad_process, ONLY: init_radiation_temperature
!    USE mo_rad_util,    ONLY: radiation_depth_to_layers_2d
    USE mo_rad_constants,  ONLY: AlbedoCanopySnow_age, AlbedoCanopySnow_temp, AlbedoVisInitial, AlbedoNirInitial, &
      &                          AlbedoGlacierVisMax, AlbedoGlacierNirMax, AlbedoLakeWater

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    dsl4jsb_Def_config(RAD_)
    dsl4jsb_Def_memory(RAD_)

    INTEGER :: iradiation
    REAL(wp), POINTER ::         &
      & return_pointer  (:,  :), & !< temporary pointer
      & return_pointer3d(:,:,:)
    REAL(WP) :: &
      albedo_age_weight

    TYPE(t_input_file) :: input_file
    INTEGER :: i, n1, n2
                         
    CHARACTER(len=*), PARAMETER :: routine = modname//':rad_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get radiation memory of the tile
    dsl4jsb_Get_memory(RAD_)

    ! Get radiation config
    dsl4jsb_Get_config(RAD_)

    IF (.NOT. dsl4jsb_Config(RAD_)%active) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing radiation memory for tile '//TRIM(tile%name))

    ! Initialize physical state
    !

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(RAD_)%bc_filename), model%grid_id)

    ! TBD: adapt names of variables in input files to names in memory

    ! background albedo (without snow or canopy)
    return_pointer => input_file%Read_2d( &
      & variable_name='albedo',           &
      & fill_array = dsl4jsb_var_ptr(RAD_,alb_background))
    dsl4jsb_var2D_onDomain(RAD_,alb_background) = &
      MERGE(dsl4jsb_var2D_onDomain(RAD_,alb_background), 0._wp, dsl4jsb_var2D_onDomain(RAD_,alb_background) >= 0._wp)

    IF (tile%is_lake) THEN
      dsl4jsb_var2D_onDomain(RAD_,alb_vis) = AlbedoLakeWater
      dsl4jsb_var2D_onDomain(RAD_,alb_nir) = AlbedoLakeWater
    ELSE IF (tile%is_glacier) THEN
      dsl4jsb_var2D_onDomain(RAD_,alb_vis)     = AlbedoGlacierVisMax
      dsl4jsb_var2D_onDomain(RAD_,alb_nir)     = AlbedoGlacierNirMax
      dsl4jsb_var2D_onDomain(RAD_,alb_vis_lnd) = AlbedoGlacierVisMax
      dsl4jsb_var2D_onDomain(RAD_,alb_nir_lnd) = AlbedoGlacierNirMax
    ELSE
      dsl4jsb_var2D_onDomain(RAD_,alb_vis)     = AlbedoVisInitial
      dsl4jsb_var2D_onDomain(RAD_,alb_nir)     = AlbedoNirInitial
      dsl4jsb_var2D_onDomain(RAD_,alb_vis_lnd) = AlbedoVisInitial
      dsl4jsb_var2D_onDomain(RAD_,alb_nir_lnd) = AlbedoNirInitial
    END IF

    ! soil albedo
    IF (tile%contains_soil) THEN
      return_pointer => input_file%Read_2d( &
        & variable_name='albedo_soil_vis',  &
        & fill_array = dsl4jsb_var_ptr(RAD_, alb_vis_soil))
      dsl4jsb_var2D_onDomain(RAD_, alb_vis_soil) = &
        MERGE(dsl4jsb_var2D_onDomain(RAD_,alb_vis_soil), 0._wp, dsl4jsb_var2D_onDomain(RAD_,alb_vis_soil) >= 0._wp)

      return_pointer => input_file%Read_2d( &
        & variable_name='albedo_soil_nir',  &
        & fill_array = dsl4jsb_var_ptr(RAD_, alb_nir_soil))
      dsl4jsb_var2D_onDomain(RAD_,alb_nir_soil) = &
        MERGE(dsl4jsb_var2D_onDomain(RAD_,alb_nir_soil), 0._wp, dsl4jsb_var2D_onDomain(RAD_,alb_nir_soil) >= 0._wp)
    END IF

    IF (tile%contains_vegetation) THEN
      ! Get canopy albedo
      IF (dsl4jsb_Config(RAD_)%use_alb_canopy .OR. tile%lcts(1)%lib_id == 0) THEN
        return_pointer => input_file%Read_2d( &
          & variable_name='albedo_veg_vis',   &
          & fill_array = dsl4jsb_var_ptr(RAD_, alb_vis_can))
        dsl4jsb_var2D_onDomain(RAD_, alb_vis_can) = &
          & MERGE(dsl4jsb_var2D_onDomain(RAD_, alb_vis_can), 0._wp, dsl4jsb_var2D_onDomain(RAD_, alb_vis_can) >= 0._wp)
      ELSE
        dsl4jsb_var2D_onDomain(RAD_, alb_vis_can) = dsl4jsb_Lctlib_param(AlbedoCanopyVIS)
      END IF
      IF (dsl4jsb_Config(RAD_)%use_alb_canopy .OR. tile%lcts(1)%lib_id == 0) THEN
        return_pointer => input_file%Read_2d( &
          & variable_name='albedo_veg_nir',   &
          & fill_array = dsl4jsb_var_ptr(RAD_, alb_nir_can))
        dsl4jsb_var2D_onDomain(RAD_, alb_nir_can) = &
          & MERGE(dsl4jsb_var2D_onDomain(RAD_, alb_nir_can), 0._wp, dsl4jsb_var2D_onDomain(RAD_, alb_nir_can) >= 0._wp)
      ELSE
        dsl4jsb_var2D_onDomain(RAD_, alb_nir_can) = dsl4jsb_Lctlib_param(AlbedoCanopyNIR)
      END IF

      ! The albedo of snow on canopy:
      ! R: In the model the albedo of snow is weighted with albedo_age_weight between two albedo schemes (age and temperature).
      !    For snow on soil this is done with calculations of snow age and temperature and with differentiation between vis and nir
      !    (see mo_rad_process in subroutine calc_snow_albedo).
      !    However, for snow ON CANOPY the albedo from the age scheme and the albedo from the temperature scheme are simply assumed
      !    to be constant (see the following lines).
      !    Therefore as all inputs to this equation are constants, it was not necessary to calculate it for each time step and block
      !    again (as it was in JSBACH3 in mo_land_surface: update_albedo_snowage_temp). Therefore it put it here and not in
      !    merge_albedos_of_vegtile.
      albedo_age_weight = dsl4jsb_Config(RAD_)%albedo_age_weight

      ! Weight between age and temperature impact on snow albedo for snow on canopy
      dsl4jsb_Config(RAD_)%AlbedoCanopySnow = &
        & albedo_age_weight * AlbedoCanopySnow_age + (1._wp - albedo_age_weight) * AlbedoCanopySnow_temp

    END IF

    CALL input_file%Close()

  END SUBROUTINE rad_init_ic

#endif
END MODULE mo_rad_init
