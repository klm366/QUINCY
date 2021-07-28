!
! This file was automatically generated from JSBACH source file:
!    src/hydrology/mo_hydro_init.f90
! with md5 hash: 72e10c0d48ee564b66837147a18bd362
!
!>
!! @brief Initialization of the the hydrology memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-09)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_hydro_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_grid_class,    ONLY: t_jsb_vgrid
  USE mo_jsb_grid,          ONLY: Get_vgrid
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input

  USE mo_jsb_process_class, ONLY: hydro_, pheno_
  USE mo_hydro_config_class, ONLY: t_hydro_config
  USE mo_hydro_memory_class, ONLY: t_hydro_memory
  USE mo_pheno_memory_class, ONLY: t_pheno_memory

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

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init'

    CALL hydro_init_bc(tile)
    CALL hydro_init_ic(tile)

  END SUBROUTINE hydro_init

  SUBROUTINE hydro_init_bc(tile)

    USE mo_util,    ONLY: soil_depth_to_layers_2d

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_hydro_config), POINTER :: hydro__conf
    TYPE(t_hydro_memory), POINTER :: hydro__mem

    REAL(wp), POINTER :: &
      & ptr_2D(:,  :) !,   &
      !& ptr_3D(:,:,:)

    TYPE(t_input_file) :: input_file
    TYPE(t_jsb_vgrid), POINTER :: soil_w
    INTEGER :: i, nsoil

    TYPE(t_jsb_model),       POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init_bc'

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting  hydrology boundary conditions for tile '//TRIM(tile%name))

    model => Get_model(tile%owner_model_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    ! Get hydrology config
    SELECT TYPE (cc => model%processes(hydro_)%p%config); TYPE IS (t_hydro_config) ; hydro__conf => cc ; END SELECT

    ! Get hydrology memory of the tile
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    !----------

    input_file = jsb_netcdf_open_input(TRIM(hydro__conf%bc_sso_filename), model%grid_id)

    ! Orography ...
    ptr_2D => input_file%Read_2d(  &
      & variable_name='elevation', &
      & fill_array = hydro__mem%elevation%ptr)
    hydro__mem%elevation%ptr(:,:) = &
      MERGE(hydro__mem%elevation%ptr(:,:), 0._wp, hydro__mem%elevation%ptr(:,:) >= 0._wp)

    ptr_2D => input_file%Read_2d(  &
      & variable_name='orostd',    &
      & fill_array = hydro__mem%oro_stddev%ptr)
    hydro__mem%oro_stddev%ptr(:,:) = &
      & MERGE(hydro__mem%oro_stddev%ptr(:,:), 0._wp, hydro__mem%oro_stddev%ptr(:,:) >= 0._wp)


    !$acc update device(hydro__mem%elevation%ptr(:,:), &
    !$acc hydro__mem%oro_stddev%ptr(:,:))
    CALL input_file%Close()

    IF (tile%contains_soil) THEN

      ! Initialize hydrology variables from bc file

      input_file = jsb_netcdf_open_input(TRIM(hydro__conf%bc_filename), model%grid_id)

      ! Total soil depth until bedrock
      ptr_2D => input_file%Read_2d(   &
        & variable_name='soil_depth', &
        & fill_array = hydro__mem%soil_depth%ptr)
      hydro__mem%soil_depth%ptr(:,:) = &
        MAX(0.1_wp, MERGE(hydro__mem%soil_depth%ptr(:,:), 0._wp, hydro__mem%soil_depth%ptr(:,:) >= 0._wp))

      ! Compute actual soil layer depths based on the fixed layer thicknesses from input file and
      ! the bedrock depth just read in.
      hydro__mem%soil_depth_sl%ptr(:,:,:) = soil_depth_to_layers_2d( &
        & hydro__mem%soil_depth%ptr(:,:),                          & ! Total soil depth until bedrock (from textures)
        & soil_w%dz(:))                                                          ! Soil layer thicknesses

      ! Maximum root zone soil moisture
      ptr_2D => input_file%Read_2d(                        &
        & fill_array = hydro__mem%max_moist%ptr, &
        & variable_name='maxmoist')
      ! max_moist is zero over glaciers ... but since we do this only for tiles with soil, we set max_moist to 0.2 over glaciers
      hydro__mem%max_moist%ptr(:,:) = &
        MERGE(hydro__mem%max_moist%ptr(:,:), 0.2_wp, hydro__mem%max_moist%ptr(:,:) > 0._wp)
      IF (hydro__conf%w_soil_limit > 0._wp) THEN
        hydro__mem%max_moist%ptr(:,:) = &
          & MIN(hydro__mem%max_moist%ptr(:,:), hydro__conf%w_soil_limit)
      END IF

      ! Volumetric soil porosity
      ptr_2D => input_file%Read_2d(                            &
        & fill_array = hydro__mem%vol_porosity%ptr,  &
        & variable_name='soil_porosity')
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) >= 0._wp)

      ! Saturated hydraulic conductivity
      ptr_2D => input_file%Read_2d(                 &
        & variable_name='hyd_cond_sat',             &
        & fill_array = hydro__mem%hyd_cond_sat%ptr)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 4.e-6_wp, ptr_2D(:,:) >= 0._wp)

      ! Matrix potential
      ptr_2D => input_file%Read_2d(                 &
        & variable_name='moisture_pot',             &
        & fill_array = hydro__mem%matrix_pot%ptr)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) >= 0._wp)

      ! Exponent B in Clapp and Hornberger
      ptr_2D => input_file%Read_2d(                 &
        & variable_name='bclapp',                   &
        & fill_array = hydro__mem%bclapp%ptr)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 6._wp, ptr_2D(:,:) >= 0._wp)

      ! Pore size distribution index
      ptr_2D => input_file%Read_2d(                 &
        & variable_name='pore_size_index',          &
        & fill_array = hydro__mem%pore_size_index%ptr)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.25_wp, ptr_2D(:,:) >= 0._wp)

      ! Volumetric field capacity [m/m]
      ptr_2D => input_file%Read_2d(       &
        & variable_name='soil_field_cap', & ! TODO change name in infile to w_soil_fcvol or another name that indicates "volumetric"
        & fill_array = hydro__mem%vol_field_cap%ptr)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

      ! Get water content at field capacity by multiplying volumetric field capacity with soil layer thicknesses
      DO i=1,nsoil
        hydro__mem%w_soil_fc_sl%ptr (:,i,:) = ptr_2D(:,:) * hydro__mem%soil_depth_sl%ptr (:,i,:)
      END DO
      NULLIFY(ptr_2D)

      ! Volumetric permanent wilting point
      ! @todo: Do we need this for bare soil tiles?
      ptr_2D => input_file%Read_2d(                                 &
        & variable_name='wilting_point',                            &
        & fill_array=hydro__mem%vol_p_wilt%ptr)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

      ! Get permanent wilting point by multiplying with soil layer thicknesses
      DO i=1,nsoil
        hydro__mem%w_soil_pwp_sl%ptr (:,i,:) = ptr_2D(:,:) * hydro__mem%soil_depth_sl%ptr (:,i,:)
      END DO
      NULLIFY(ptr_2D)

      ! Get water holding capacity (water content at saturation) by multiplying volumetric porosity with soil layer thicknesses
      DO i=1,nsoil
        hydro__mem%w_soil_sat_sl%ptr (:,i,:) = &
          hydro__mem%vol_porosity%ptr(:,:) * hydro__mem%soil_depth_sl%ptr (:,i,:)
      END DO

      !$acc update device(hydro__mem%w_soil_sat_sl%ptr, &
      !$acc hydro__mem%w_soil_pwp_sl%ptr, &
      !$acc hydro__mem%w_soil_fc_sl%ptr, &
      !$acc hydro__mem%vol_field_cap%ptr, &
      !$acc hydro__mem%vol_p_wilt%ptr, &
      !$acc hydro__mem%pore_size_index%ptr, &
      !$acc hydro__mem%bclapp%ptr, &
      !$acc hydro__mem%matrix_pot%ptr, &
      !$acc hydro__mem%hyd_cond_sat%ptr, &
      !$acc hydro__mem%vol_porosity%ptr, &
      !$acc hydro__mem%max_moist%ptr(:,:), &
      !$acc hydro__mem%soil_depth_sl%ptr(:,:,:), &
      !$acc hydro__mem%soil_depth%ptr(:,:) )
      CALL input_file%Close()

    END IF

    IF (tile%contains_vegetation) THEN

      input_file = jsb_netcdf_open_input(TRIM(hydro__conf%bc_filename), model%grid_id)

      ptr_2D => input_file%Read_2d(   &
        & variable_name='root_depth', &
        & fill_array = hydro__mem%root_depth%ptr)
      hydro__mem%root_depth%ptr(:,:) = &
        MERGE(hydro__mem%root_depth%ptr(:,:), 0._wp, hydro__mem%root_depth%ptr(:,:) >= 0._wp)
      hydro__mem%root_depth%ptr(:,:) = &           ! Temporary fix (root depth sometimes larger then soil depth)
        MIN(hydro__mem%root_depth%ptr(:,:), hydro__mem%soil_depth%ptr(:,:))
      hydro__mem%root_depth_sl%ptr(:,:,:) = soil_depth_to_layers_2d( &
        & hydro__mem%root_depth%ptr(:,:),                          & ! Total rooting depth
        & soil_w%dz(:))    ! Soil layer thicknesses

      !$acc update device(hydro__mem%root_depth%ptr(:,:), &
      !$acc hydro__mem%root_depth_sl%ptr(:,:,:))
      CALL input_file%Close()

    END IF

  END SUBROUTINE hydro_init_bc


  SUBROUTINE hydro_init_ic(tile)

    !-----------------------------------------------------------------------
    !  DECLARATIONS
    USE mo_hydro_constants, ONLY: vol_porosity_org, vol_field_cap_org, vol_p_wilt_org
    USE mo_hydro_util,      ONLY: get_water_in_root_zone

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    TYPE(t_hydro_config), POINTER :: hydro__conf

    TYPE(t_hydro_memory), POINTER :: hydro__mem
    TYPE(t_pheno_memory), POINTER :: pheno__mem

    REAL(wp), POINTER :: &
      & ptr_2D(:,  :),   & !< temporary pointers
      & ptr_3D(:,:,:)

    TYPE(t_input_file) :: input_file
    TYPE(t_jsb_vgrid), POINTER :: soil_w
    INTEGER :: i, isoil, nsoil

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init_ic'

    !-----------------------------------------------------------------------
    ! CONTENT

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing hydrology memory for tile '//TRIM(tile%name))

    model => Get_model(tile%owner_model_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    ! Get hydrology config
    SELECT TYPE (cc => model%processes(hydro_)%p%config); TYPE IS (t_hydro_config) ; hydro__conf => cc ; END SELECT

    ! Get memory of the tile
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT
    IF (tile%contains_vegetation) THEN
      SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT
    END IF

    !----------

    IF (tile%contains_soil) THEN

      ! Initialize hydrology variables from ic file

      input_file = jsb_netcdf_open_input(hydro__conf%ic_filename, model%grid_id)

      ! Initial snow depth
      IF (input_file%Has_var('snow')) THEN
        ptr_2D => input_file%Read_2d( &
          & variable_name='snow',     &
          & fill_array = hydro__mem%w_snow_soil%ptr)
        hydro__mem%w_snow_soil%ptr(:,:) = &
          & MERGE(hydro__mem%w_snow_soil%ptr(:,:), 0._wp, hydro__mem%w_snow_soil%ptr(:,:) >= 0._wp)
      ELSE
        hydro__mem%w_snow_soil%ptr(:,:) = 0._wp
        CALL message(TRIM(routine), 'Initializing snow to zero')
      END IF
      hydro__mem%w_snow%ptr(:,:) = hydro__mem%w_snow_soil%ptr(:,:)

      ! Initial soil moisture
      ! ptr_2D => input_file%Read_2d(         &
      !   & variable_name='init_moist',               &
      !   & fill_array = hydro__mem%w_soil_column%ptr)
      ! hydro__mem%w_soil_column%ptr(:,:) = &
      !   & MERGE(hydro__mem%w_soil_column%ptr(:,:), 0._wp, hydro__mem%w_soil_column%ptr(:,:) >= 0._wp)

      ! Initial w_soil_overflown water
      hydro__mem%w_soil_overflow%ptr(:,:) = 0._wp
      DO isoil=1,nsoil
        ptr_3D => input_file%Read_2d_extdim( &
          & variable_name='layer_moist',     &
          & start_extdim=isoil, end_extdim=isoil, extdim_name='soillev')
        hydro__mem%w_soil_sl%ptr (:,isoil,:) = &
          & MERGE(ptr_3D(:,:,1), 0._wp, ptr_3D(:,:,1) >= 0._wp)
      END DO

      CALL input_file%Close()

      !$acc update device(hydro__mem%w_snow%ptr(:,:),       &
      !$acc               hydro__mem%w_snow_soil%ptr(:,:),  &
      !$acc               hydro__mem%w_soil_overflow%ptr(:,:))
      ! avoid compiler warnings about pointers not being dereferenced

      IF (ASSOCIATED(ptr_2D)) CONTINUE
      NULLIFY(ptr_2D)

    END IF

    ! Initialize organic soil layer fractions and modify saturated, field capacity and perm. wilting point capacities accordingly
    ! Note: at the moment, fract_org_sl doesn't change over time, but it could
    IF (hydro__conf%l_organic .AND. tile%contains_soil) THEN
      hydro__mem%fract_org_sl%ptr(:,:,:) = 0._wp
      IF (tile%is_vegetation .AND. tile%lcts(1)%lib_id > 0 .AND. model%lctlib(tile%lcts(1)%lib_id)%ForestFlag) THEN
        ! For PFT tiles
        hydro__mem%fract_org_sl%ptr(:,1,:) = pheno__mem%fract_fpc_max%ptr(:,:)
      ELSE IF (tile%is_vegetation) THEN
        ! For general VEG tile
        hydro__mem%fract_org_sl%ptr(:,1,:) = pheno__mem%fract_forest%ptr(:,:)
      END IF
      ! The following can change later depending on the organic layers and soil ice content.
      hydro__mem%w_soil_sat_sl%ptr(:,:,:) =                                                        &
        & (1 - hydro__mem%fract_org_sl%ptr(:,:,:)) * hydro__mem%w_soil_sat_sl%ptr(:,:,:) &
        & +    hydro__mem%fract_org_sl%ptr(:,:,:)  * hydro__mem%soil_depth_sl%ptr(:,:,:) * vol_porosity_org
      hydro__mem%w_soil_fc_sl%ptr(:,:,:) =                                                         &
        & (1 - hydro__mem%fract_org_sl%ptr(:,:,:)) * hydro__mem%w_soil_fc_sl%ptr(:,:,:)  &
        & +    hydro__mem%fract_org_sl%ptr(:,:,:)  * hydro__mem%soil_depth_sl%ptr(:,:,:) * vol_field_cap_org
      hydro__mem%w_soil_pwp_sl%ptr(:,:,:) =                                                        &
        & (1 - hydro__mem%fract_org_sl%ptr(:,:,:)) * hydro__mem%w_soil_pwp_sl%ptr(:,:,:) &
        & +    hydro__mem%fract_org_sl%ptr(:,:,:)  * hydro__mem%soil_depth_sl%ptr(:,:,:) * vol_p_wilt_org
      
      !$acc update device(hydro__mem%fract_org_sl%ptr(:,:,:), &
      !$acc hydro__mem%w_soil_sat_sl%ptr(:,:,:), &
      !$acc hydro__mem%w_soil_fc_sl%ptr(:,:,:), &
      !$acc hydro__mem%w_soil_pwp_sl%ptr(:,:,:))
    END IF

    IF (tile%contains_soil) THEN
      ! @todo: The interpolation of the initial layer moisture from T63 to the ICON grid leads to positive values where the
      ! water content should actually be zero. Therefore, make sure that water content is not larger than field capacity.
      hydro__mem%w_soil_sl%ptr(:,:,:) = &
        & MIN(hydro__mem%w_soil_sl%ptr(:,:,:), hydro__mem%w_soil_fc_sl%ptr(:,:,:))
      !$acc update device(hydro__mem%w_soil_sl%ptr(:,:,:))

      DO i=1,SIZE(hydro__mem%w_soil_root%ptr, 2)
          CALL get_water_in_root_zone(                       &
          & hydro__mem%w_soil_sl%ptr     (:,:,i),  &
          & hydro__mem%soil_depth_sl%ptr (:,:,i),  &
          & hydro__mem%root_depth_sl%ptr (:,:,i),  &
          & hydro__mem%w_soil_column%ptr (:,  i)   )
      END DO
      !$acc update host(hydro__mem%w_soil_column%ptr)
      WHERE (hydro__mem%w_soil_column%ptr(:,:) > hydro__mem%max_moist%ptr(:,:))
        hydro__mem%w_soil_column%ptr(:,:) = hydro__mem%max_moist%ptr(:,:)
      END WHERE

      WHERE (hydro__mem%max_moist%ptr(:,:) > 0._wp)
        hydro__mem%w_soil_rel%ptr(:,:) = &
          & hydro__mem%w_soil_column%ptr(:,:) / hydro__mem%max_moist%ptr(:,:)
      ELSE WHERE
        hydro__mem%w_soil_rel%ptr(:,:) = 0._wp
      END WHERE

      
      !$acc update device(hydro__mem%w_soil_rel%ptr(:,:),    &
      !$acc               hydro__mem%w_soil_column%ptr(:,:), &
      !$acc               hydro__mem%w_soil_sl%ptr(:,:,:))

      CALL init_soil_properties(tile)
    END IF

    IF (tile%contains_vegetation) THEN
      ! Calculate water content and field capacity in the root zone
      DO i=1,SIZE(hydro__mem%w_soil_root%ptr, 2)
          CALL get_water_in_root_zone( &
          & hydro__mem%w_soil_sl%ptr      (:,:,i), & 
          & hydro__mem%soil_depth_sl%ptr  (:,:,i), &
          & hydro__mem%root_depth_sl%ptr  (:,:,i), &
          & hydro__mem%w_soil_root%ptr    (:,  i)  )
 
          CALL get_water_in_root_zone( & 
          & hydro__mem%w_soil_fc_sl%ptr   (:,:,i), & 
          & hydro__mem%soil_depth_sl%ptr  (:,:,i), &
          & hydro__mem%root_depth_sl%ptr  (:,:,i), &
          & hydro__mem%w_soil_root_fc%ptr (:,  i)  )
          
          CALL get_water_in_root_zone( & 
          & hydro__mem%w_soil_pwp_sl%ptr  (:,:,i), & 
          & hydro__mem%soil_depth_sl%ptr  (:,:,i), &
          & hydro__mem%root_depth_sl%ptr  (:,:,i), &
          & hydro__mem%w_soil_root_pwp%ptr(:,  i)  )
      END DO
      !$acc update host(hydro__mem%w_soil_root%ptr,     &
      !$acc             hydro__mem%w_soil_root_fc%ptr,  &
      !$acc             hydro__mem%w_soil_root_pwp%ptr)
    END IF

  END SUBROUTINE hydro_init_ic

  SUBROUTINE init_soil_properties(tile)

    USE mo_hydro_constants, ONLY: vol_porosity_org, hyd_cond_sat_org, bclapp_org, matrix_pot_org, pore_size_index_org, &
      & vol_field_cap_org, vol_p_wilt_org

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! Local variables
    !
    TYPE(t_hydro_config), POINTER :: hydro__conf
    TYPE(t_hydro_memory), POINTER :: hydro__mem

    ! Pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & hyd_cond_sat,          &
      & vol_porosity,          &
      & bclapp,                &
      & matrix_pot,            &
      & pore_size_index,       &
      & vol_field_cap,         &
      & vol_p_wilt
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & fract_org_sl,          &
      & hyd_cond_sat_sl,       &
      & vol_porosity_sl,       &
      & bclapp_sl,             &
      & matrix_pot_sl,         &
      & pore_size_index_sl,    &
      & vol_field_cap_sl,      &
      & vol_p_wilt_sl

    ! Locally allocated vectors
    !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_w

    INTEGER  :: nsoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_soil_properties'

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (tile%is_lake) RETURN

    model => Get_model(tile%owner_model_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    SELECT TYPE (cc => model%processes(hydro_)%p%config); TYPE IS (t_hydro_config) ; hydro__conf => cc ; END SELECT

    ! Get reference to variables for current block
    !
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    hyd_cond_sat => hydro__mem%hyd_cond_sat%ptr(:,:)
    vol_porosity => hydro__mem%vol_porosity%ptr(:,:)
    bclapp => hydro__mem%bclapp%ptr(:,:)
    matrix_pot => hydro__mem%matrix_pot%ptr(:,:)
    pore_size_index => hydro__mem%pore_size_index%ptr(:,:)
    vol_field_cap => hydro__mem%vol_field_cap%ptr(:,:)
    vol_p_wilt => hydro__mem%vol_p_wilt%ptr(:,:)
    IF (hydro__conf%l_organic) fract_org_sl => hydro__mem%fract_org_sl%ptr(:,:,:)
    hyd_cond_sat_sl => hydro__mem%hyd_cond_sat_sl%ptr(:,:,:)
    vol_porosity_sl => hydro__mem%vol_porosity_sl%ptr(:,:,:)
    bclapp_sl => hydro__mem%bclapp_sl%ptr(:,:,:)
    matrix_pot_sl => hydro__mem%matrix_pot_sl%ptr(:,:,:)
    pore_size_index_sl => hydro__mem%pore_size_index_sl%ptr(:,:,:)
    vol_field_cap_sl => hydro__mem%vol_field_cap_sl%ptr(:,:,:)
    vol_p_wilt_sl => hydro__mem%vol_p_wilt_sl%ptr(:,:,:)

    hyd_cond_sat_sl   (:,:,:) = SPREAD(hyd_cond_sat   (:,:), DIM=2, ncopies=nsoil)
    vol_porosity_sl   (:,:,:) = SPREAD(vol_porosity   (:,:), DIM=2, ncopies=nsoil)
    bclapp_sl         (:,:,:) = SPREAD(bclapp         (:,:), DIM=2, ncopies=nsoil)
    matrix_pot_sl     (:,:,:) = SPREAD(matrix_pot     (:,:), DIM=2, ncopies=nsoil)
    pore_size_index_sl(:,:,:) = SPREAD(pore_size_index(:,:), DIM=2, ncopies=nsoil)
    vol_field_cap_sl  (:,:,:) = SPREAD(vol_field_cap  (:,:), DIM=2, ncopies=nsoil)
    vol_p_wilt_sl     (:,:,:) = SPREAD(vol_p_wilt     (:,:), DIM=2, ncopies=nsoil)
    IF (hydro__conf%l_organic) THEN
      hyd_cond_sat_sl   (:,:,:) = (1 - fract_org_sl(:,:,:)) * hyd_cond_sat_sl   (:,:,:) + fract_org_sl(:,:,:) * hyd_cond_sat_org
      vol_porosity_sl   (:,:,:) = (1 - fract_org_sl(:,:,:)) * vol_porosity_sl   (:,:,:) + fract_org_sl(:,:,:) * vol_porosity_org
      bclapp_sl         (:,:,:) = (1 - fract_org_sl(:,:,:)) * bclapp_sl         (:,:,:) + fract_org_sl(:,:,:) * bclapp_org
      matrix_pot_sl     (:,:,:) = (1 - fract_org_sl(:,:,:)) * matrix_pot_sl     (:,:,:) + fract_org_sl(:,:,:) * matrix_pot_org
      pore_size_index_sl(:,:,:) = (1 - fract_org_sl(:,:,:)) * pore_size_index_sl(:,:,:) + fract_org_sl(:,:,:) * pore_size_index_org
      vol_field_cap_sl  (:,:,:) = (1 - fract_org_sl(:,:,:)) * vol_field_cap_sl  (:,:,:) + fract_org_sl(:,:,:) * vol_field_cap_org
      vol_p_wilt_sl     (:,:,:) = (1 - fract_org_sl(:,:,:)) * vol_p_wilt_sl     (:,:,:) + fract_org_sl(:,:,:) * vol_p_wilt_org
    END IF

    !$acc update device(hyd_cond_sat_sl, vol_porosity_sl, bclapp_sl, &
    !$acc matrix_pot_sl, pore_size_index_sl, vol_field_cap_sl, &
    !$acc vol_p_wilt_sl)

  END SUBROUTINE init_soil_properties

#endif
END MODULE mo_hydro_init
