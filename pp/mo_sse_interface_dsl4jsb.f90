!
! This file was automatically generated from:
!    /users/rschnur/icon-les-dev-staging.pgi_gpu/externals/jsbach/src/soil_snow_energy/mo_sse_interface.f90
!    with md5 hash: 511c25ee6d452ce35969545512e3fdd8
!
!>
!! @brief Contains the interfaces to the soil and snow energy processes
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

MODULE mo_sse_interface
#ifndef __NO_JSBACH__

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_jsb_time,        ONLY: is_time_experiment_start

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  USE mo_jsb_process_class, ONLY: sse_, seb_, hydro_ 
  USE mo_sse_config_class, ONLY: t_sse_config
  USE mo_hydro_config_class, ONLY: t_hydro_config
  USE mo_sse_memory_class, ONLY: t_sse_memory
  USE mo_seb_memory_class, ONLY: t_seb_memory
  USE mo_hydro_memory_class, ONLY: t_hydro_memory

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_sse_tasks

  !> Type definition for soil_temp_and_hydro task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_and_snow_temp
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_and_snow_temperature
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_and_snow_temperature
  END TYPE tsk_soil_and_snow_temp

  !> Constructor interface for soil_temperature task
  INTERFACE tsk_soil_and_snow_temp
    PROCEDURE Create_task_soil_and_snow_temperature
  END INTERFACE tsk_soil_and_snow_temp

  !> Type definition for soil_properties
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_and_snow_properties
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_and_snow_properties
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_and_snow_properties
  END TYPE tsk_soil_and_snow_properties

  !> Constructor interface for soil_properties
  INTERFACE tsk_soil_and_snow_properties
    PROCEDURE Create_task_soil_and_snow_properties
  END INTERFACE tsk_soil_and_snow_properties

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_interface'

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for soil_temperature
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "soil_temperature"
  !!
  FUNCTION Create_task_soil_and_snow_temperature(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_soil_and_snow_temp::return_ptr)
    CALL return_ptr%Construct(name='soil_and_snow_temperature', process_id=SSE_, owner_model_id=model_id)

  END FUNCTION Create_task_soil_and_snow_temperature

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for soil_properties task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "soil_properties"
  !!
  FUNCTION Create_task_soil_and_snow_properties(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_soil_and_snow_properties::return_ptr)
    CALL return_ptr%Construct(name='soil_and_snow_properties', process_id=SSE_, owner_model_id=model_id)

  END FUNCTION Create_task_soil_and_snow_properties

  SUBROUTINE Register_sse_tasks(process, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,              INTENT(in)    :: model_id

    CALL process%Register_task(tsk_soil_and_snow_temp(model_id))
    CALL process%Register_task(tsk_soil_and_snow_properties(model_id))

  END SUBROUTINE Register_sse_tasks

  SUBROUTINE update_soil_and_snow_temperature(tile, options)

    USE mo_sse_process,            ONLY: calc_soil_temperature, calc_snow_temperature, calc_soiltemp_old
    USE mo_jsb_physical_constants, ONLY: alf, rhoh2o, dens_snow
    USE mo_sse_constants,          ONLY: snow_depth_min

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_sse_config), POINTER :: sse__conf
    TYPE(t_sse_memory), POINTER :: sse__mem
    TYPE(t_seb_memory), POINTER :: seb__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem

    ! Pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: t_unfilt
    REAL(wp), POINTER, DIMENSION(:) :: snowmelt
    REAL(wp), POINTER, DIMENSION(:,:) :: t_soil
    REAL(wp), POINTER, DIMENSION(:,:) :: t_soil_acoef
    REAL(wp), POINTER, DIMENSION(:,:) :: t_soil_bcoef
    REAL(wp), POINTER, DIMENSION(:,:) :: vol_heat_cap_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: heat_cond_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: t_snow
    REAL(wp), POINTER, DIMENSION(:,:) :: t_snow_acoef
    REAL(wp), POINTER, DIMENSION(:,:) :: t_snow_bcoef
    REAL(wp), POINTER, DIMENSION(:) :: vol_heat_cap_snow
    REAL(wp), POINTER, DIMENSION(:) :: heat_cond_snow
    REAL(wp), POINTER, DIMENSION(:) :: hcap_grnd
    REAL(wp), POINTER, DIMENSION(:) :: hcap_grnd_old
    REAL(wp), POINTER, DIMENSION(:) :: grnd_hflx
    REAL(wp), POINTER, DIMENSION(:) :: grnd_hflx_old
    REAL(wp), POINTER, DIMENSION(:,:) :: snow_depth_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: soil_depth_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: w_soil_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: w_ice_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: w_soil_freeze_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: w_ice_melt_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: vol_porosity_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: matrix_pot_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: bclapp_sl
    REAL(wp), POINTER, DIMENSION(:) :: thaw_depth
    REAL(wp), POINTER, DIMENSION(:) :: heat_cond
    REAL(wp), POINTER, DIMENSION(:) :: vol_heat_cap
    REAL(wp), POINTER, DIMENSION(:) :: w_snow_soil
    REAL(wp), POINTER, DIMENSION(:) :: snow_soil_dens
    REAL(wp), POINTER, DIMENSION(:) :: fract_snow_soil


    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc) :: &
      & snow_depth,                    &
      & snow_fract,                    &
      & t_unfilt_corrected,            &
      & t_soil_top,                    &
      & hcap_grnd_snow,                &
      & hcap_grnd_soil,                &
      & grnd_hflx_snow,                &
      & grnd_hflx_soil
    LOGICAL, DIMENSION(options%nc)  :: &
      & ldglac
    REAL(wp), ALLOCATABLE :: ws_max_sl(:,:)

    INTEGER, DIMENSION(options%nc) ::  &
      & itop, itop_old

    INTEGER :: isnow

    LOGICAL  :: lstart
    INTEGER  :: iblk, ics, ice, nc, i
    REAL(wp) :: dtime

    TYPE(t_jsb_vgrid), POINTER :: soil_e, snow_e
    INTEGER                    :: nsoil, nsnow

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_and_snow_temperature'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    model => Get_model(tile%owner_model_id)
    SELECT TYPE (cc => model%processes(sse_)%p%config); TYPE IS (t_sse_config) ; sse__conf => cc ; END SELECT

    soil_e => Get_vgrid('soil_depth_energy')
    nsoil = soil_e%n_levels
    IF (sse__conf%l_snow) THEN
      snow_e => Get_vgrid('snow_depth_energy')
      nsnow = snow_e%n_levels
    END IF

    ! Get reference to variables for current block
    !
    SELECT TYPE (mm => tile%mem(sse_)%p); TYPE IS (t_sse_memory); sse__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    snowmelt => hydro__mem%snowmelt%ptr(ics:ice,iblk)
    w_snow_soil => hydro__mem%w_snow_soil%ptr(ics:ice,iblk)
    snow_soil_dens => hydro__mem%snow_soil_dens%ptr(ics:ice,iblk)
    fract_snow_soil => hydro__mem%fract_snow_soil%ptr(ics:ice,iblk)

    t_soil => sse__mem%t_soil%ptr(ics:ice,:,iblk)            ! inout
    t_soil_acoef => sse__mem%t_soil_acoef%ptr(ics:ice,:,iblk)      ! inout
    t_soil_bcoef => sse__mem%t_soil_bcoef%ptr(ics:ice,:,iblk)      ! inout
    vol_heat_cap_sl => sse__mem%vol_heat_cap_sl%ptr(ics:ice,:,iblk)   ! in
    heat_cond_sl => sse__mem%heat_cond_sl%ptr(ics:ice,:,iblk)      ! in
    vol_heat_cap_snow => sse__mem%vol_heat_cap_snow%ptr(ics:ice,iblk) ! in
    heat_cond_snow => sse__mem%heat_cond_snow%ptr(ics:ice,iblk)    ! in
    hcap_grnd => sse__mem%hcap_grnd%ptr(ics:ice,iblk)         ! out
    hcap_grnd_old => sse__mem%hcap_grnd_old%ptr(ics:ice,iblk)     ! out
    grnd_hflx => sse__mem%grnd_hflx%ptr(ics:ice,iblk)         ! out
    grnd_hflx_old => sse__mem%grnd_hflx_old%ptr(ics:ice,iblk)     ! out
    IF (sse__conf%l_snow) THEN
      snow_depth_sl => sse__mem%snow_depth_sl%ptr(ics:ice,:,iblk)   ! out
      t_snow => sse__mem%t_snow%ptr(ics:ice,:,iblk)          ! inout
      t_snow_acoef => sse__mem%t_snow_acoef%ptr(ics:ice,:,iblk)    ! inout
      t_snow_bcoef => sse__mem%t_snow_bcoef%ptr(ics:ice,:,iblk)    ! inout
    END IF
    IF (.NOT. tile%is_glacier) THEN
      soil_depth_sl => hydro__mem%soil_depth_sl%ptr(ics:ice,:,iblk)    ! in
      vol_porosity_sl => hydro__mem%vol_porosity_sl%ptr(ics:ice,:,iblk)  ! in
      matrix_pot_sl => hydro__mem%matrix_pot_sl%ptr(ics:ice,:,iblk)    ! in
      bclapp_sl => hydro__mem%bclapp_sl%ptr(ics:ice,:,iblk)        ! in
      w_soil_sl => hydro__mem%w_soil_sl%ptr(ics:ice,:,iblk)        ! inout
      w_ice_sl => hydro__mem%w_ice_sl%ptr(ics:ice,:,iblk)         ! inout
      w_soil_freeze_sl => hydro__mem%w_soil_freeze_sl%ptr(ics:ice,:,iblk) ! out
      w_ice_melt_sl => hydro__mem%w_ice_melt_sl%ptr(ics:ice,:,iblk)    ! out
    END IF
    thaw_depth => sse__mem%thaw_depth%ptr(ics:ice,iblk)      ! out

    heat_cond => sse__mem%heat_cond%ptr(ics:ice,iblk)
    vol_heat_cap => sse__mem%vol_heat_cap%ptr(ics:ice,iblk)

    t_unfilt => seb__mem%t_unfilt%ptr(ics:ice,iblk)        ! in

    lstart = is_time_experiment_start()

    !$claw acc enter data create(snow_depth,snow_fract,t_unfilt_corrected,t_soil_top,hcap_grnd_snow,hcap_grnd_soil,grnd_hflx_snow,grnd_hflx_soil,ldglac,itop,itop_old)

    ! Save old ground heat flux and heat capacity before updating
    ! The old flux and capacity (which have been used in the implicit calculation of the new surface temperature) are
    ! given back to the atmosphere so that the surface energy balance can be computed there
    !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_0 tile_name=tile%name
    grnd_hflx_old(:) = grnd_hflx(:)
    hcap_grnd_old(:) = hcap_grnd(:)
    !$claw end expand

    ! Correction of surface temperature for snow and ice melt
    ! t_unfilt is later corrected in task "snowmelt_correction" of process SEB_, but we need it already here on the sub-tiles
    IF (lstart) THEN
      !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_1 tile_name=tile%name
      t_unfilt_corrected(:) = t_unfilt(:)
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_2 tile_name=tile%name
      t_unfilt_corrected(:) = t_unfilt(:) - snowmelt(:) * dtime * alf / hcap_grnd(:)
    END IF

    IF (.NOT. sse__conf%l_snow) THEN
      ldglac(:) = SPREAD(tile%is_glacier, DIM=1, ncopies=nc)
      !$claw acc update device(ldglac)

      CALL calc_soiltemp_old(                          &
        & nc,                                          &
        & nsoil,                                       &
        & soil_e%dz(:),                                &
        & soil_e%levels(:),                            &
        & dtime,                                       &
        & lstart,                                      &
        & t_unfilt_corrected(:)                     ,  &
        & w_snow_soil(:)                            ,  &
        & heat_cond(:)                              ,  &
        & vol_heat_cap(:)                           ,  &
        & t_soil_acoef                          (:,:), &
        & t_soil_bcoef                          (:,:), &
        & t_soil                                (:,:), &
        & hcap_grnd                             (:),   &
        & grnd_hflx                             (:),   &
        & ldglac                                (:)    &
        !& SPREAD(tile%is_glacier, DIM=1, ncopies=nc)   &
        & )


      !$claw acc exit data delete(snow_depth,snow_fract,t_unfilt_corrected,t_soil_top,hcap_grnd_snow,hcap_grnd_soil,grnd_hflx_snow,grnd_hflx_soil,ldglac,itop,itop_old)  
      RETURN

    END IF

    ! l_snow=.TRUE. below here

    IF (.NOT. tile%is_glacier) THEN

      ! Calculate snow depth in m (no longer in m water equivalent)
      IF (sse__conf%l_dynsnow) THEN
        !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_3 tile_name=tile%name
        snow_depth(:) = w_snow_soil(:) * rhoh2o / snow_soil_dens(:)
      ELSE
        !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_4 tile_name=tile%name
        snow_depth(:) = w_snow_soil(:) * rhoh2o / dens_snow
      END IF

      ! Calculate top snow layer for previous time step
      !$claw acc data present(itop_old,snow_depth_sl)
      !$claw acc parallel loop gang vector
      DO i = 1,nc
        itop_old(i) = nsnow + 1 - COUNT(snow_depth_sl(i,:) > 0._wp)
      END DO
      !$claw acc end data

      ! Calculate depths of snow layers
      !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_5 tile_name=tile%name
      snow_depth_sl(:,:) = 0._wp
      !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_6 tile_name=tile%name
      snow_fract(:) = fract_snow_soil(:)


      !$claw acc data present(snow_depth,snow_fract)
      !$claw acc parallel loop gang vector
      DO i = 1, size(snow_depth, dim=1)
        IF (snow_depth(i) > snow_depth_min .AND. snow_fract(i) > EPSILON(1._wp)) THEN
          snow_depth(i) = snow_depth(i) / snow_fract(i)
        ELSE
          snow_depth(i) = 0._wp
          snow_fract(i) = 0._wp
        END IF
      END DO
      !$claw acc end data

      !$claw acc data present(snow_depth,snow_depth_sl)
      !$claw acc parallel
      !$claw acc loop seq
      DO isnow=nsnow,1,-1
        !$claw acc loop vector
        DO i = 1, size(snow_depth, dim=1)
          IF (snow_depth(i) > 0._wp) THEN
            snow_depth_sl(i,isnow) = MIN(snow_depth(i), snow_e%dz(isnow))
            snow_depth(i) = snow_depth(i) - snow_depth_sl(i,isnow)
          END IF
        END DO
      END DO
      !$claw acc end parallel
      !$claw acc end data

      !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_7 tile_name=tile%name
      snow_depth_sl(:,nsnow) = snow_depth_sl(:,nsnow) + snow_depth(:)

      ! Calculate top snow layer for current time step

      !$claw acc data present(itop,snow_depth_sl)
      !$claw acc parallel loop gang vector
      DO i = 1,nc
        itop(i) = nsnow + 1 - COUNT(snow_depth_sl(i,:) > 0._wp)
      END DO
      !$claw acc end data

      CALL calc_snow_temperature( &
        ! in
        & nc,                     &
        & nsnow,                  &
        & soil_e%dz(1:2),         &
        & dtime,                  &
        & lstart,                 &
        & snow_depth_sl    (:,:), &
        & itop_old           (:), &
        & itop               (:), &
        & t_unfilt_corrected (:), &
        & vol_heat_cap_snow  (:), &
        & heat_cond_snow     (:), &
        & vol_heat_cap_sl  (:,1), &
        & heat_cond_sl     (:,1), &
        & t_soil           (:,1), &
        & t_soil_acoef     (:,1), &
        & t_soil_bcoef     (:,1), &
        ! inout
        & t_snow           (:,:), &
        & t_snow_acoef     (:,:), &
        & t_snow_bcoef     (:,:), &
        ! out
        & t_soil_top         (:), &
        & hcap_grnd_snow     (:), &
        & grnd_hflx_snow     (:)  &
        & )

      DO i = 1, size(snow_depth_sl, dim=1)
        IF (snow_depth_sl(i, nsnow) > 0._wp) THEN  ! At least one snow layer above soil
          t_soil_top(i) = snow_fract(i) * t_soil_top(i) + (1._wp - snow_fract(i)) * t_unfilt_corrected(i)
        ELSE                                ! No snow
          t_soil_top(i) = t_unfilt_corrected(i)
        END IF
      END DO
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_8 tile_name=tile%name
      snow_depth_sl(:,:) = 0._wp

      !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_9 tile_name=tile%name
      snow_depth(:) = w_snow_soil(:) * rhoh2o / dens_snow
      snow_fract(:) = 0._wp
      t_soil_top(:) = t_unfilt_corrected(:)
      grnd_hflx_snow(:) = 0._wp
      hcap_grnd_snow(:) = 0._wp
      !$claw end expand
    END IF

    ! @todo: move water/ice calculations to HYDRO_ process
    IF (tile%is_glacier) THEN
      CALL calc_soil_temperature( &
        ! in
        & nc                                       , &
        & nsoil                                    , &
        & soil_e%dz(:)                             , &
        & soil_e%levels(:)                         , &
        & dtime                                    , &
        & lstart                                   , &
        & sse__conf%l_freeze            , &
        & sse__conf%l_supercool         , &
        & t_soil_top                          (:  ), &
        & vol_heat_cap_sl                     (:,:), &
        & heat_cond_sl                        (:,:), &
        ! inout
        & t_soil                              (:,:), &
        & t_soil_acoef                        (:,:), &
        & t_soil_bcoef                        (:,:), &
        ! out
        & hcap_grnd_soil                      (:  ), &
        & grnd_hflx_soil                      (:  ), &
        & thaw_depth                          (:  )  &
        & )
    ELSE
      ALLOCATE(ws_max_sl(nc,nsoil))
      !$claw acc enter data create(ws_max_sl)
      !$claw expand parallel update
      ws_max_sl(:,:) = soil_depth_sl(:,:) * vol_porosity_sl(:,:)

      CALL calc_soil_temperature( &
        ! in
        & nc                                       , &
        & nsoil                                    , &
        & soil_e%dz(:)                             , &
        & soil_e%levels(:)                         , &
        & dtime                                    , &
        & lstart                                   , &
        & sse__conf%l_freeze            , &
        & sse__conf%l_supercool         , &
        & t_soil_top                          (:  ), &
        & vol_heat_cap_sl                     (:,:), &
        & heat_cond_sl                        (:,:), &
        ! inout
        & t_soil                              (:,:), &
        & t_soil_acoef                        (:,:), &
        & t_soil_bcoef                        (:,:), &
        ! out
        & hcap_grnd_soil                      (:  ), &
        & grnd_hflx_soil                      (:  ), &
        & thaw_depth                          (:  ), &
        ! optional, in
        & ws_max_sl                           (:,:), &
        & matrix_pot_sl                       (:,:), &
        & bclapp_sl                           (:,:), &
        ! optional, inout
        & w_soil_sl                           (:,:), &
        & w_ice_sl                            (:,:), &
        ! optional, out
        & w_soil_freeze_sl                    (:,:), &
        & w_ice_melt_sl                       (:,:)  &
        & )
      !$claw acc exit data delete(ws_max_sl)
      DEALLOCATE(ws_max_sl)
    END IF

    !$claw acc data present(grnd_hflx,hcap_grnd,snow_fract,grnd_hflx_snow,grnd_hflx_soil,hcap_grnd_snow,hcap_grnd_soil)
    !$claw acc parallel loop gang vector
    DO i = 1, size(snow_depth_sl, dim=1)
      IF (snow_depth_sl(i,nsnow) > 0._wp) THEN
        grnd_hflx(i) = snow_fract(i) * grnd_hflx_snow(i) + (1._wp - snow_fract(i)) * grnd_hflx_soil(i)
        hcap_grnd(i) = snow_fract(i) * hcap_grnd_snow(i) + (1._wp - snow_fract(i)) * hcap_grnd_soil(i)
      ELSE
        grnd_hflx(i) = grnd_hflx_soil(i)
        hcap_grnd(i) = hcap_grnd_soil(i)
      END IF
    END DO
    !$claw acc end data

    !$claw acc exit data delete(snow_depth,snow_fract,t_unfilt_corrected,t_soil_top,hcap_grnd_snow,hcap_grnd_soil,grnd_hflx_snow,grnd_hflx_soil,ldglac,itop,itop_old)  

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_and_snow_temperature

  SUBROUTINE aggregate_soil_and_snow_temperature(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_sse_config), POINTER :: sse__conf
    TYPE(t_sse_memory), POINTER :: sse__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_and_snow_temperature'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    SELECT TYPE (cc => model%processes(sse_)%p%config); TYPE IS (t_sse_config) ; sse__conf => cc ; END SELECT
    SELECT TYPE (mm => tile%mem(sse_)%p); TYPE IS (t_sse_memory); sse__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, sse__mem%t_soil_acoef,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, sse__mem%t_soil_bcoef,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, sse__mem%t_soil,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, sse__mem%grnd_hflx,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, sse__mem%grnd_hflx_old,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, sse__mem%hcap_grnd,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, sse__mem%hcap_grnd_old,ics,ice,iblk,routine)
    IF (sse__conf%l_snow) THEN
      CALL weighted_by_fract%Aggregate(tile, sse__mem%snow_depth_sl,ics,ice,iblk,routine)
      CALL weighted_by_fract%Aggregate(tile, sse__mem%t_snow_acoef,ics,ice,iblk,routine)
      CALL weighted_by_fract%Aggregate(tile, sse__mem%t_snow_bcoef,ics,ice,iblk,routine)
      CALL weighted_by_fract%Aggregate(tile, sse__mem%t_snow,ics,ice,iblk,routine)
    END IF
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_soil_freeze_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_ice_melt_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, sse__mem%thaw_depth,ics,ice,iblk,routine)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_and_snow_temperature

  SUBROUTINE update_soil_and_snow_properties(tile, options)

    USE mo_jsb_physical_constants, ONLY: ci, dens_snow_min ! , cs
    USE mo_sse_constants,          ONLY: vol_hcap_ice, vol_hcap_snow, hcond_ice, hcond_snow, hcond_snow_min !, vol_hcap_snow_min
    USE mo_sse_process,            ONLY: calc_vol_heat_capacity, calc_thermal_conductivity

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_sse_config), POINTER :: sse__conf
    TYPE(t_sse_memory), POINTER :: sse__mem
    TYPE(t_hydro_config), POINTER :: hydro__conf
    TYPE(t_hydro_memory), POINTER :: hydro__mem

    ! Pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: &
      & snow_soil_dens,       &
      & vol_heat_cap,         &
      & heat_cond,            &
      & vol_heat_cap_snow,    &
      & heat_cond_snow
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & soil_depth_sl,        &
      & fract_org_sl,         &
      & vol_heat_cap_sl,      &
      & heat_cond_sl,         & 
      & w_soil_sl,            &
      & w_ice_sl,             &
      & vol_porosity_sl,      &
      & vol_field_cap_sl

    ! Locally allocated vectors
    !
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & fract_water,          &
      & fract_ice

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_e

    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: & 
      & vol_heat_cap_tmp,                    & 
      & soil_e_dz_tmp,                       &
      & heat_cond_tmp

    INTEGER  :: iblk, ics, ice, nc, i, j
    INTEGER  :: nsoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_and_snow_properties'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    soil_e => Get_vgrid('soil_depth_energy')
    nsoil = soil_e%n_levels

    SELECT TYPE (cc => model%processes(sse_)%p%config); TYPE IS (t_sse_config) ; sse__conf => cc ; END SELECT
    SELECT TYPE (cc => model%processes(hydro_)%p%config); TYPE IS (t_hydro_config) ; hydro__conf => cc ; END SELECT

    ! Get reference to variables for domain
    !
    SELECT TYPE (mm => tile%mem(sse_)%p); TYPE IS (t_sse_memory); sse__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    vol_heat_cap_sl => sse__mem%vol_heat_cap_sl%ptr(ics:ice,:,iblk)   ! out
    heat_cond_sl => sse__mem%heat_cond_sl%ptr(ics:ice,:,iblk)      ! out
    vol_heat_cap_snow => sse__mem%vol_heat_cap_snow%ptr(ics:ice,iblk) ! out
    heat_cond_snow => sse__mem%heat_cond_snow%ptr(ics:ice,iblk) ! out


    ! On glacier tile, set variables and return
    IF (tile%is_glacier) THEN
      IF (model%config%l_compat401) THEN
        !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_properties_expand_0 tile_name=tile%name
        vol_heat_cap_sl(:,:) = 2.09e+6_wp
        heat_cond_sl(:,:) = 2.508_wp
        !$claw end expand
        !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_properties_expand_1 tile_name=tile%name
        vol_heat_cap_snow(:) = 634500._wp
        heat_cond_snow(:) = 0.31_wp
        !$claw end expand
      ELSE
        !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_properties_expand_2 tile_name=tile%name
        vol_heat_cap_sl  (:,:) = vol_hcap_ice
        heat_cond_sl     (:,:) = hcond_ice
        !$claw end expand
        !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_properties_expand_3 tile_name=tile%name
        vol_heat_cap_snow(:  ) = vol_hcap_snow
        heat_cond_snow   (:  ) = hcond_snow
        !$claw end expand
      END IF
      IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
      RETURN
    END IF

    soil_depth_sl => hydro__mem%soil_depth_sl%ptr(ics:ice,:,iblk)
    vol_heat_cap => sse__mem%vol_heat_cap%ptr(ics:ice,iblk)
    heat_cond => sse__mem%heat_cond%ptr(ics:ice,iblk)
    w_soil_sl => hydro__mem%w_soil_sl%ptr(ics:ice,:,iblk)
    w_ice_sl => hydro__mem%w_ice_sl%ptr(ics:ice,:,iblk)
    vol_porosity_sl => hydro__mem%vol_porosity_sl%ptr(ics:ice,:,iblk)
    vol_field_cap_sl => hydro__mem%vol_field_cap_sl%ptr(ics:ice,:,iblk)  

    IF (hydro__conf%l_organic .AND. tile%contains_soil) THEN
      fract_org_sl => hydro__mem%fract_org_sl%ptr(ics:ice,:,iblk)
    ELSE
      ALLOCATE(fract_org_sl(nc, nsoil))
      !$claw acc enter data pcreate(fract_org_sl)
      !$claw expand parallel
      fract_org_sl(:,:) = 0._wp
    END IF

    IF (sse__conf%l_heat_cap_dyn) THEN
      ALLOCATE(fract_water(nc, nsoil), fract_ice(nc, nsoil))
      !$claw acc enter data pcreate(fract_water, fract_ice)

      !$claw acc data present(soil_depth_sl,fract_water,fract_ice,w_soil_sl,w_ice_sl)
      !$claw acc parallel loop gang vector collapse(2)
      DO i = 1, size(soil_depth_sl, dim=1)
        DO j = 1, size(soil_depth_sl, dim=2)
          IF (soil_depth_sl(i,j) > 0._wp) THEN
            fract_water(i,j) = w_soil_sl(i,j) / soil_depth_sl(i,j)
            fract_ice  (i,j) = w_ice_sl(i,j) / soil_depth_sl(i,j)
          ELSE
            fract_water(i,j) = 0._wp
            fract_ice  (i,j) = 0._wp
          END IF
        END DO
      END DO
      !$claw acc end data

      ALLOCATE(vol_heat_cap_tmp(nc,nsoil), soil_e_dz_tmp(nc,nsoil))
      !$claw acc enter data pcreate(vol_heat_cap_tmp,soil_e_dz_tmp)
      !$claw acc update host(vol_heat_cap(:))
      vol_heat_cap_tmp(:,:) = SPREAD( vol_heat_cap(:), DIM=2, ncopies=nsoil)
      soil_e_dz_tmp(:,:) = SPREAD( soil_e%dz(:), DIM=1, ncopies=nc)
      !$claw acc update device(vol_heat_cap_tmp,soil_e_dz_tmp)
 
      !$claw expand parallel update savepoint=ser_jsb_sse_calc_vol_heat_capacity_1 tile_name=tile%name
      vol_heat_cap_sl (:,:) = calc_vol_heat_capacity( &
        !&                        SPREAD( vol_heat_cap(:),                 DIM=2, ncopies=nsoil)             , &
        &                                vol_heat_cap_tmp(:,:)                                              , &
        &                                vol_porosity_sl(:,:)                                               , &
        &                                fract_org_sl(:,:)                                                  , &
        &                                fract_water(:,:), &
        &                                fract_ice  (:,:), &
        &                                soil_depth_sl(:,:)                                                 , &
        !&                        SPREAD( soil_e%dz(:),                    DIM=1, ncopies=nc)                  &
        &                                soil_e_dz_tmp(:,:)                                                   &
        & )
      
      !$claw acc exit data delete(fract_water,fract_ice,vol_heat_cap_tmp,soil_e_dz_tmp)
      DEALLOCATE(fract_water, fract_ice, vol_heat_cap_tmp, soil_e_dz_tmp)
    ELSE
      !$claw acc data present(vol_heat_cap_sl,vol_heat_cap)
      !$claw acc parallel loop gang vector collapse(2)
      DO i = 1, nc
        DO j = 1, nsoil
          vol_heat_cap_sl(i,j) = vol_heat_cap(i)
        END DO
      END DO
      !$claw acc end data
    ENDIF

    IF (sse__conf%l_heat_cond_dyn) THEN
      
      ALLOCATE(heat_cond_tmp(nc,nsoil), soil_e_dz_tmp(nc,nsoil))
      !$claw acc enter data pcreate(heat_cond_tmp,soil_e_dz_tmp)
      !$claw acc update host(heat_cond(:))
      heat_cond_tmp(:,:) = SPREAD( heat_cond(:), DIM=2, ncopies=nsoil)
      soil_e_dz_tmp(:,:) = SPREAD( soil_e%dz(:), DIM=1, ncopies=nc)
      !$claw acc update device(heat_cond_tmp, soil_e_dz_tmp)

      !$claw expand parallel update savepoint=ser_jsb_sse_calc_thermal_conductivity_1 tile_name=tile%name
      heat_cond_sl (:,:) = calc_thermal_conductivity( &
        !&                        SPREAD( heat_cond(:),                    DIM=2, ncopies=nsoil)             , &
        &                                heat_cond_tmp(:,:)                                                 , &
        &                                w_soil_sl(:,:)                                                     , &
        &                                w_ice_sl(:,:)                                                      , &
        &                                vol_porosity_sl(:,:)                                               , &
        &                                vol_field_cap_sl(:,:)                                              , &
        &                                fract_org_sl(:,:)                                                  , &
        &                                soil_depth_sl(:,:)                                                 , &
        &                                soil_e_dz_tmp(:,:)                                                   &
        !&                        SPREAD( soil_e%dz(:),                    DIM=1, ncopies=nc)                  &
        & )

      !$claw acc exit data delete(heat_cond_tmp,soil_e_dz_tmp)
      DEALLOCATE(heat_cond_tmp, soil_e_dz_tmp)
    ELSE
      !$claw acc data present(heat_cond_sl,heat_cond)
      !$claw acc parallel loop gang vector collapse(2)
      DO i = 1, nc
        DO j = 1, nsoil
          heat_cond_sl(i,j) = heat_cond(i)
        END DO
      END DO
      !$claw acc end data
    END IF

    IF (.NOT. (hydro__conf%l_organic .AND. tile%contains_soil)) THEN
      !$claw acc exit data delete(fract_org_sl)
      DEALLOCATE(fract_org_sl)
    END IF

    IF (model%config%l_compat401) THEN
      !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_properties_expand_4 tile_name=tile%name
      vol_heat_cap_snow(:) = 634500._wp
      heat_cond_snow(:) = 0.31_wp
      !$claw end expand
    ELSE
      IF (sse__conf%l_snow) THEN
        IF (sse__conf%l_dynsnow) THEN
          snow_soil_dens => hydro__mem%snow_soil_dens%ptr(ics:ice,iblk)
          !$claw acc data present(snow_soil_dens,vol_heat_cap_snow,heat_cond_snow)
          !$claw acc parallel loop gang vector
          DO i = 1, size(snow_soil_dens)
            IF (snow_soil_dens(i) > dens_snow_min) THEN     ! "old" or mixture of old and fresh snow
              vol_heat_cap_snow(i) = ci * snow_soil_dens(i)                               ! Use specific heat of ice (Verseghy 1991)
              ! heat_cond_snow   (:) = 2.576e-6_wp * snow_soil_dens(:) ** 2._wp + 0.074_wp                      ! Verseghy (1991)
              ! heat_cond_snow   (:) = 2.9e-6_wp*snow_soil_dens(:)**2._wp                                       ! Goodrich (1982)
              ! Calonne et al. (2011):
              heat_cond_snow   (i) = 2.5e-6_wp*snow_soil_dens(i)**2._wp - 0.123e-3_wp*snow_soil_dens(i) + 0.024_wp
            ELSE                                         ! fresh snow
              ! @todo: these minimum values from mo_sse_constants.f90 don't match with the above equations and minimum snow density
              ! vol_heat_cap_snow(:) = vol_hcap_snow_min
              heat_cond_snow   (i) = hcond_snow_min
              vol_heat_cap_snow(i) = ci * dens_snow_min
            END IF
          END DO
          !$claw acc end data
        ELSE
          !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_properties_expand_5  tile_name=tile%name
          vol_heat_cap_snow(:) = vol_hcap_snow                                     ! Use specific heat of snow (cs * dens_snow)
          heat_cond_snow   (:) = hcond_snow
          !$claw end expand
        END IF
      ELSE
        !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_properties_expand_6  tile_name=tile%name
        vol_heat_cap_snow(:) = vol_hcap_snow                                     ! Use specific heat of snow (cs * dens_snow)
        heat_cond_snow   (:) = hcond_snow
        !$claw end expand
      END IF
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_and_snow_properties
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "soil_properties"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_soil_and_snow_properties(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_sse_memory), POINTER :: sse__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_and_snow_properties'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(sse_)%p); TYPE IS (t_sse_memory); sse__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, sse__mem%vol_heat_cap_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, sse__mem%heat_cond_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, sse__mem%vol_heat_cap_snow,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, sse__mem%heat_cond_snow,ics,ice,iblk,routine)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_and_snow_properties

#endif
END MODULE mo_sse_interface
