!
! This file was automatically generated from:
!    /scratch/snx3000/rschnur/icon-rc-testing-plus-mc8.test/externals/jsbach/src/soil_snow_energy/mo_sse_init.f90
!    with md5 hash: b4009f8ada900180c4688731b922c2ea
!
!>
!! @brief Initialization of the the soil and snow memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_sse_init
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

  USE mo_jsb_process_class, ONLY: sse_, seb_, hydro_
  USE mo_sse_config_class, ONLY: t_sse_config
  USE mo_hydro_config_class, ONLY: t_hydro_config
  USE mo_sse_memory_class, ONLY: t_sse_memory
  USE mo_seb_memory_class, ONLY: t_seb_memory
  USE mo_hydro_memory_class, ONLY: t_hydro_memory

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sse_init

  REAL(wp), PARAMETER :: &
    & fao_thermal_diff(0:6) = (/7.4e-7_wp ,8.7e-7_wp ,8.0e-7_wp ,7.4e-7_wp ,7.1e-7_wp ,6.7e-7_wp ,-1._wp/), &
    & fao_vol_hcap    (0:6) = (/2.25e+6_wp,1.93e+6_wp,2.10e+6_wp,2.25e+6_wp,2.36e+6_wp,2.48e+6_wp,-1._wp/)

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_init'

CONTAINS

  !
  !> Intialize soil and snow energy process (after memory has been set up)
  !
  SUBROUTINE sse_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init'

    CALL sse_init_bc(tile)
    CALL sse_init_ic(tile)

  END SUBROUTINE sse_init

  SUBROUTINE sse_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_sse_config), POINTER :: sse__conf
    TYPE(t_sse_memory), POINTER :: sse__mem

    REAL(wp), POINTER :: &
      & fao    (:,  :), & !< FAO soil types
      & ptr_2D (:,  :)    !< temporary pointer

    INTEGER :: iblk

    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get soil config
    SELECT TYPE (cc => model%processes(sse_)%p%config); TYPE IS (t_sse_config) ; sse__conf => cc ; END SELECT

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting  soil energy boundary conditions for tile '//TRIM(tile%name))

    ! Get soil memory of the tile
    SELECT TYPE (mm => tile%mem(sse_)%p); TYPE IS (t_sse_memory); sse__mem => mm ; END SELECT

    input_file = jsb_netcdf_open_input(TRIM(sse__conf%bc_filename), model%grid_id)

    ! Thermal diffusivity of soil
    fao => input_file%Read_2d(                    &
      & variable_name='fao')
    fao(:,:) = MERGE(fao(:,:), 0._wp, fao(:,:) >= 0._wp)
    fao(:,:) = MAX(REAL(LBOUND(fao_thermal_diff,1),wp), MIN(REAL(UBOUND(fao_thermal_diff,1),wp), fao(:,:)))

    DO iblk=1,SIZE(fao,2)
      sse__mem%thermal_diffusivity%ptr(:,iblk) = fao_thermal_diff(NINT(fao(:,iblk)))
    END DO

    ! Volumetric heat capacity of dry soil [J/(m^3*K)]
    IF (sse__conf%l_heat_cap_map) THEN
      ptr_2D => input_file%Read_2d(         &
        & variable_name='heat_capacity',            &
        & fill_array = sse__mem%vol_heat_cap%ptr)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 2.e6_wp, ptr_2D(:,:) > 0._wp)
    ELSE
      DO iblk=1,SIZE(fao,2)
        sse__mem%vol_heat_cap%ptr(:,iblk) = fao_vol_hcap(NINT(fao(:,iblk)))
      END DO
    END IF

    ! Heat conductivity of mineral soil
    IF (sse__conf%l_heat_cond_map) THEN
      ptr_2D => input_file%Read_2d(         &
        & variable_name='heat_conductivity',         &
        & fill_array = sse__mem%heat_cond%ptr)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) > 0._wp)
    ELSE
      sse__mem%heat_cond%ptr(:,:) = &
        & sse__mem%vol_heat_cap%ptr(:,:) * sse__mem%thermal_diffusivity%ptr(:,:)
    END IF

    !$acc data present(sse__mem%thermal_diffusivity%ptr(:,:), &
    !$acc sse__mem%vol_heat_cap%ptr(:,:), &
    !$acc sse__mem%heat_cond%ptr(:,:))
    !$acc update device(sse__mem%thermal_diffusivity%ptr(:,:), &
    !$acc sse__mem%vol_heat_cap%ptr(:,:), &
    !$acc sse__mem%heat_cond%ptr(:,:))
    !$acc end data

    CALL input_file%Close()                              ! Temporary

  END SUBROUTINE sse_init_bc

  SUBROUTINE sse_init_ic(tile)

    USE mo_sse_process, ONLY: init_soil_temperature

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    TYPE(t_sse_config), POINTER :: sse__conf
    TYPE(t_sse_memory), POINTER :: sse__mem
    TYPE(t_seb_memory), POINTER :: seb__mem

    REAL(wp), POINTER :: &
      & t_srf_mon (:,:,:) !,    & !< Climatological monthly values of surface temperature
      !& ptr_2D    (:,  :),    & !< temporary pointer
      !& ptr_3D    (:,:,:)

    TYPE(t_input_file) :: input_file
    TYPE(t_jsb_vgrid), POINTER :: soil_e
    INTEGER :: iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init_ic'

    model => get_model(tile%owner_model_id)

    soil_e => Get_vgrid('soil_depth_energy')

    ! Get soil memory of the tile
    SELECT TYPE (mm => tile%mem(sse_)%p); TYPE IS (t_sse_memory); sse__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT

    ! Get soil config
    SELECT TYPE (cc => model%processes(sse_)%p%config); TYPE IS (t_sse_config) ; sse__conf => cc ; END SELECT

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing soil energy memory for tile '//TRIM(tile%name))

    input_file = jsb_netcdf_open_input(sse__conf%ic_filename, model%grid_id)

    ! Soil and surface temperatures
    !
    ! Read climatological surface temperature
    t_srf_mon => input_file%Read_2d_time(         &   ! Temporary
      & variable_name='surf_temp')
    t_srf_mon(:,:,:) = MERGE(t_srf_mon(:,:,:), 280._wp, t_srf_mon(:,:,:) > 0._wp)

    DO iblk=1,SIZE(t_srf_mon, DIM=2)
      CALL init_soil_temperature(                     &
        & soil_e%levels(:),                           &
        & t_srf_mon(:,iblk,1:12),                     &
        & sse__mem%t_soil%ptr(:,:,iblk),   &
        & seb__mem%t%ptr(:,iblk))
    END DO
    seb__mem%t_unfilt%ptr(:,:) = seb__mem%t%ptr(:,:)

    !$acc data present(seb__mem%t_unfilt%ptr(:,:), & 
    !$acc seb__mem%t%ptr(:,:), &
    !$acc sse__mem%t_soil%ptr(:,:,:))
    !$acc update device(seb__mem%t_unfilt%ptr(:,:), & 
    !$acc seb__mem%t%ptr(:,:), &
    !$acc sse__mem%t_soil%ptr(:,:,:))
    !$acc end data

    CALL input_file%Close()

    CALL init_soil_properties(tile)

  END SUBROUTINE sse_init_ic

  SUBROUTINE init_soil_properties(tile)

    USE mo_jsb_physical_constants, ONLY: ci, dens_snow_min
    USE mo_sse_constants,          ONLY: vol_hcap_ice, vol_hcap_snow, & ! vol_hcap_snow_min,
      &                                  hcond_ice, hcond_snow !, hcond_snow_min, hcond_org
    USE mo_sse_process,            ONLY: calc_vol_heat_capacity, calc_thermal_conductivity

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! Local variables
    !
    TYPE(t_sse_config), POINTER :: sse__conf
    TYPE(t_sse_memory), POINTER :: sse__mem
    TYPE(t_hydro_config), POINTER :: hydro__conf
    TYPE(t_hydro_memory), POINTER :: hydro__mem

    ! Pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & vol_heat_cap, &
      & heat_cond, &
      & vol_heat_cap_snow, &
      & heat_cond_snow
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & soil_depth_sl, &
      & fract_org_sl, &
      & vol_heat_cap_sl, &
      & heat_cond_sl, &
      & w_soil_sl, &
      & w_ice_sl, &
      & vol_porosity_sl, &
      & vol_field_cap_sl

    ! Locally allocated vectors
    ! 
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & fract_water,           &
      & fract_ice,             &
      & soil_e_dz_tmp,         & 
      & vol_heat_cap_tmp,      &
      & heat_cond_tmp

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_e

    INTEGER  :: nsoil, nproma, nblks, i, j, b

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_soil_properties'

    model => Get_model(tile%owner_model_id)

    soil_e => Get_vgrid('soil_depth_energy')
    nsoil = soil_e%n_levels

    SELECT TYPE (cc => model%processes(sse_)%p%config); TYPE IS (t_sse_config) ; sse__conf => cc ; END SELECT
    SELECT TYPE (cc => model%processes(hydro_)%p%config); TYPE IS (t_hydro_config) ; hydro__conf => cc ; END SELECT

    ! Get reference to variables for domain
    !
    SELECT TYPE (mm => tile%mem(sse_)%p); TYPE IS (t_sse_memory); sse__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    vol_heat_cap_sl => sse__mem%vol_heat_cap_sl%ptr(:,:,:)
    heat_cond_sl => sse__mem%heat_cond_sl%ptr(:,:,:)
    vol_heat_cap_snow => sse__mem%vol_heat_cap_snow%ptr(:,:)
    heat_cond_snow => sse__mem%heat_cond_snow%ptr(:,:)

    ! On glacier tile, set variables and return
    IF (tile%is_glacier) THEN
      !$claw expand parallel update
      vol_heat_cap_sl  (:,:,:) = vol_hcap_ice
      !$claw expand parallel update
      heat_cond_sl     (:,:,:) = hcond_ice
      !$cla$w end expand
      !$claw expand parallel update
      vol_heat_cap_snow(:  ,:) = vol_hcap_snow
      !$claw expand parallel update
      heat_cond_snow   (:  ,:) = hcond_snow
      !$cla$w end expand
      RETURN
    END IF

    w_soil_sl => hydro__mem%w_soil_sl%ptr(:,:,:)
    w_ice_sl => hydro__mem%w_ice_sl%ptr(:,:,:)
    vol_porosity_sl => hydro__mem%vol_porosity_sl%ptr(:,:,:)
    vol_field_cap_sl => hydro__mem%vol_field_cap_sl%ptr(:,:,:)
    soil_depth_sl => hydro__mem%soil_depth_sl%ptr(:,:,:)
    vol_heat_cap => sse__mem%vol_heat_cap%ptr(:,:)
    heat_cond => sse__mem%heat_cond%ptr(:,:)
    nproma = SIZE(heat_cond,1)
    nblks  = SIZE(heat_cond,2)
    IF (hydro__conf%l_organic .AND. tile%contains_soil) THEN
      fract_org_sl => hydro__mem%fract_org_sl%ptr(:,:,:)
    ELSE
      ALLOCATE(fract_org_sl(nproma, nsoil, nblks))
      !$claw acc enter data create(fract_org_sl)
      !$claw expand parallel update
      fract_org_sl(:,:,:) = 0._wp
    END IF

    ALLOCATE(soil_e_dz_tmp(nproma,nsoil,nblks),vol_heat_cap_tmp(nproma,nsoil,nblks),heat_cond_tmp(nproma,nsoil,nblks))
    !$claw acc enter data pcreate(soil_e_dz_tmp,vol_heat_cap_tmp,heat_cond_tmp)
    soil_e_dz_tmp(:,:,:) = SPREAD(SPREAD(soil_e%dz(:), DIM=1, ncopies=nproma), DIM=3, ncopies=nblks)
    vol_heat_cap_tmp(:,:,:) = SPREAD( vol_heat_cap(:,:), DIM=2, ncopies=nsoil)
    heat_cond_tmp(:,:,:) = SPREAD( heat_cond(:,:),DIM=2, ncopies=nsoil)
    !$claw acc update device(soil_e_dz_tmp,vol_heat_cap_tmp,heat_cond_tmp)

    IF (sse__conf%l_heat_cap_dyn) THEN
      ALLOCATE(fract_water(nproma, nsoil, nblks), fract_ice(nproma, nsoil, nblks))
      !$claw acc enter data pcreate(fract_water,fract_ice)

      !$claw acc data present(soil_depth_sl,fract_water,fract_ice,w_soil_sl,w_ice_sl)
      !$claw acc update device(soil_depth_sl,w_soil_sl,w_ice_sl)
      !$claw acc parallel loop gang vector collapse(3)
      DO i = 1, nproma
        DO j = 1, nsoil
          DO b = 1, nblks
            IF (soil_depth_sl(i,j,b) > 0._wp) THEN
              fract_water(i,j,b) = w_soil_sl(i,j,b) / soil_depth_sl(i,j,b)
              fract_ice  (i,j,b) = w_ice_sl(i,j,b) / soil_depth_sl(i,j,b)
            ELSE
              fract_water(i,j,b) = 0._wp
              fract_ice  (i,j,b) = 0._wp
            END IF
          END DO
        END DO
      END DO
      !$claw acc update host(fract_water,fract_ice)
      !$claw acc end data
      
      !$claw expand parallel update savepoint=ser_jsb_sse_init_calc_vol_heat_capacity_1 tile_name=tile%name
      vol_heat_cap_sl (:,:,:) = calc_vol_heat_capacity( &
        !&                         SPREAD( vol_heat_cap     (:,:),   DIM=2, ncopies=nsoil)                        , &
        &                                 vol_heat_cap_tmp (:,:,:)                                               , &
        &                                 vol_porosity_sl  (:,:,:)                                               , &
        &                                 fract_org_sl     (:,:,:)                                               , &
        &                                 fract_water      (:,:,:)                                               , &
        &                                 fract_ice        (:,:,:)                                               , &
        &                                 soil_depth_sl    (:,:,:)                                               , &
        &                                 soil_e_dz_tmp    (:,:,:)                                                 &
        !&                   SPREAD(SPREAD(soil_e%dz        (:),     DIM=1, ncopies=nproma), DIM=3, ncopies=nblks)  &
        & )

      !$claw acc exit data delete(fract_water, fract_ice)
      DEALLOCATE(fract_water, fract_ice)
    ELSE

      !$claw acc data present(vol_heat_cap_sl,vol_heat_cap)
      !$claw acc update device(vol_heat_cap)
      !$claw acc parallel loop gang vector collapse(3)
      DO i = 1, nproma
        DO j = 1, nsoil
          DO b=1,nblks
            vol_heat_cap_sl(i,j,b) = vol_heat_cap(i,b)
          END DO
        END DO
      END DO
      !$claw acc update host(vol_heat_cap_sl)
      !$claw acc end data

    END IF

    IF (sse__conf%l_heat_cond_dyn) THEN
      !$claw expand parallel update savepoint=ser_jsb_sse_init_calc_thermal_conductivity_1 tile_name=tile%name
      heat_cond_sl (:,:,:) = calc_thermal_conductivity( &
        !&                              SPREAD( heat_cond(:,:),                    DIM=2, ncopies=nsoil)               , &
        &                                      heat_cond_tmp(:,:,:)                                                   , &
        &                                      w_soil_sl(:,:,:)                                                       , &
        &                                      w_ice_sl(:,:,:)                                                        , &
        &                                      vol_porosity_sl(:,:,:)                                                 , &
        &                                      vol_field_cap_sl(:,:,:)                                                , &
        &                                      fract_org_sl(:,:,:)                                                    , &
        &                                      soil_depth_sl(:,:,:)                                                   , &
        &                                      soil_e_dz_tmp(:,:,:)                                                     &
        !&                   SPREAD(SPREAD(soil_e%dz        (:),     DIM=1, ncopies=nproma), DIM=3, ncopies=nblks)  &
        & )
    ELSE
      !$claw acc data present(heat_cond_sl,heat_cond)
      !$claw acc update device(heat_cond)
      !$claw acc parallel loop gang vector collapse(3)
      DO i = 1, nproma
        DO j = 1, nsoil
          DO b=1,nblks
            heat_cond_sl(i,j,b) = heat_cond(i,b)
          END DO
        END DO
      END DO
      !$claw acc update host(heat_cond_sl)
      !$claw acc end data
    END IF

    IF (.NOT. (hydro__conf%l_organic .AND. tile%contains_soil)) THEN
      !$claw acc exit data delete(fract_org_sl)
      DEALLOCATE(fract_org_sl)
    END IF

    IF (sse__conf%l_snow) THEN
      IF (model%config%l_compat401) THEN
        !$claw expand parallel update
        vol_heat_cap_snow(:,:) = 634500._wp
        !$claw expand parallel update
        heat_cond_snow   (:,:) = 0.31_wp
        !$cla$w end expand
      ELSE
        IF (sse__conf%l_dynsnow) THEN
          ! USe snow density values for fresh snow for initialization
          !$claw expand parallel update
          vol_heat_cap_snow(:,:) = ci * dens_snow_min                               ! Use specific heat of ice (Verseghy 1991)
          !$claw expand parallel update
          heat_cond_snow   (:,:) = 2.5e-6_wp*dens_snow_min**2._wp - 0.123e-3_wp*dens_snow_min + 0.024_wp ! Calonne et al. (2011)
          !$cla$w end expand
          ! vol_heat_cap_snow(:,:) = vol_hcap_snow_min
          ! heat_cond_snow   (:,:) = hcond_snow_min
        ELSE
          !$claw expand parallel update
          vol_heat_cap_snow(:,:) = vol_hcap_snow
          !$claw expand parallel update
          heat_cond_snow   (:,:) = hcond_snow
          !$cla$w end expand
        END IF
      END IF
    ! ELSE
    !   vol_heat_cap_snow and heat_cond_snow are fixed in calc_soiltemp_old
    END IF

    !$acc exit data delete(soil_e_dz_tmp,vol_heat_cap_tmp,heat_cond_tmp)
    DEALLOCATE(soil_e_dz_tmp,vol_heat_cap_tmp,heat_cond_tmp)

  END SUBROUTINE init_soil_properties

#endif
END MODULE mo_sse_init
