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
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_jsb_time,        ONLY: is_time_experiment_start

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes SSE_, SEB_, HYDRO_, VEG_, MIX_
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(VEG_)
  dsl4jsb_Use_memory(MIX_)

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

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: SSE_ process
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Register_sse_tasks(process, model_id)

    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,              INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model

    model   => Get_model(model_id)

    CALL process%Register_task(tsk_soil_and_snow_temp(model_id))
    CALL process%Register_task(tsk_soil_and_snow_properties(model_id))

  END SUBROUTINE Register_sse_tasks

  !-----------------------------------------------------------------------------------------------------
  !> update_soil_and_snow_temperature
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_soil_and_snow_temperature(tile, options)

    USE mo_sse_process,            ONLY: calc_soil_temperature, calc_snow_temperature, calc_soiltemp_old
    USE mo_jsb_physical_constants, ONLY: alf, rhoh2o, dens_snow
    USE mo_sse_constants,          ONLY: snow_depth_min
    USE mo_jsb_time,               ONLY: is_newyear


    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(MIX_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: snowmelt
    dsl4jsb_Real3D_onChunk :: t_soil_sl
    dsl4jsb_Real3D_onChunk :: t_soil_acoef
    dsl4jsb_Real3D_onChunk :: t_soil_bcoef
    dsl4jsb_Real3D_onChunk :: vol_heat_cap_sl
    dsl4jsb_Real3D_onChunk :: heat_cond_sl
    dsl4jsb_Real3D_onChunk :: t_snow
    dsl4jsb_Real3D_onChunk :: t_snow_acoef
    dsl4jsb_Real3D_onChunk :: t_snow_bcoef
    dsl4jsb_Real2D_onChunk :: vol_heat_cap_snow
    dsl4jsb_Real2D_onChunk :: heat_cond_snow
    dsl4jsb_Real2D_onChunk :: hcap_grnd
    dsl4jsb_Real2D_onChunk :: hcap_grnd_old
    dsl4jsb_Real2D_onChunk :: grnd_hflx
    dsl4jsb_Real2D_onChunk :: grnd_hflx_old
    dsl4jsb_Real3D_onChunk :: snow_depth_sl
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl
    dsl4jsb_Real3D_onChunk :: ice_soil_sl
    dsl4jsb_Real3D_onChunk :: wtr_freeze_sl
    dsl4jsb_Real3D_onChunk :: ice_melt_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_scool_sl
    dsl4jsb_Real3D_onChunk :: vol_porosity_sl
    dsl4jsb_Real3D_onChunk :: matrix_pot_sl
    dsl4jsb_Real3D_onChunk :: bclapp_sl
    dsl4jsb_Real2D_onChunk :: thaw_depth
    dsl4jsb_Real2D_onChunk :: thaw_depth_max
    dsl4jsb_Real2D_onChunk :: thaw_depth_max_ym1
    dsl4jsb_Real2D_onChunk :: thaw_depth_diag1
    dsl4jsb_Real2D_onChunk :: thaw_depth_max_diag1
    dsl4jsb_Real2D_onChunk :: thaw_depth_max_ym1_diag1
    dsl4jsb_Real2D_onChunk :: thaw_depth_diag2
    dsl4jsb_Real2D_onChunk :: thaw_depth_max_diag2
    dsl4jsb_Real2D_onChunk :: thaw_depth_max_ym1_diag2
    dsl4jsb_Real2D_onChunk :: thaw_depth_diag3
    dsl4jsb_Real2D_onChunk :: thaw_depth_max_diag3
    dsl4jsb_Real2D_onChunk :: thaw_depth_max_ym1_diag3
    dsl4jsb_Real2D_onChunk :: heat_cond
    dsl4jsb_Real2D_onChunk :: vol_heat_cap
    dsl4jsb_Real2D_onChunk :: weq_snow_soil
    dsl4jsb_Real2D_onChunk :: snow_soil_dens
    dsl4jsb_Real2D_onChunk :: fract_snow_soil

    dsl4jsb_Real3D_onChunk :: heat_flx_lat

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
    INTEGER  :: isoil   ! looping
    INTEGER  :: iblk, ics, ice, nc, i
    REAL(wp) :: dtime

    TYPE(t_jsb_vgrid), POINTER :: soil_e, snow_e
    TYPE(t_jsb_vgrid), POINTER :: vgrid_soil_sb                ! Vertical grid
    INTEGER                    :: nsoil, nsnow
    INTEGER                    :: nsoil_sb

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_and_snow_temperature'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (.NOT. tile%Is_process_calculated(SSE_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(SSE_)

    soil_e => Get_vgrid('soil_depth_energy')
    nsoil = soil_e%n_levels
    IF (dsl4jsb_Config(SSE_)%l_snow) THEN
      snow_e => Get_vgrid('snow_depth_energy')
      nsnow = snow_e%n_levels
    END IF
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      = vgrid_soil_sb%n_levels

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(MIX_)

    dsl4jsb_Get_var2D_onChunk(HYDRO_, snowmelt)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_soil)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snow_soil_dens)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_soil)

    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)         ! inout
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_acoef)      ! inout
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_bcoef)      ! inout
    dsl4jsb_Get_var3D_onChunk(SSE_, vol_heat_cap_sl)   ! in
    dsl4jsb_Get_var3D_onChunk(SSE_, heat_cond_sl)      ! in
    dsl4jsb_Get_var3D_onChunk(MIX_, heat_flx_lat)      ! in
    dsl4jsb_Get_var2D_onChunk(SSE_, vol_heat_cap_snow) ! in
    dsl4jsb_Get_var2D_onChunk(SSE_, heat_cond_snow)    ! in
    dsl4jsb_Get_var2D_onChunk(SSE_, hcap_grnd)         ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, hcap_grnd_old)     ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, grnd_hflx)         ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, grnd_hflx_old)     ! out
    IF (dsl4jsb_Config(SSE_)%l_snow) THEN
      dsl4jsb_Get_var3D_onChunk(SSE_, snow_depth_sl)   ! out
      dsl4jsb_Get_var3D_onChunk(SSE_, t_snow)          ! inout
      dsl4jsb_Get_var3D_onChunk(SSE_, t_snow_acoef)    ! inout
      dsl4jsb_Get_var3D_onChunk(SSE_, t_snow_bcoef)    ! inout
    END IF
    ! IF (.NOT. tile%is_glacier) THEN
      dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)    ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl)  ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, matrix_pot_sl)    ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, bclapp_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)      ! inout
      dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)      ! inout
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_freeze_sl)    ! out
      dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_melt_sl)      ! out
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pot_scool_sl) ! out
    ! END IF
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth_max)      ! out
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth_max_ym1)  ! out
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth)      ! out
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth_max_diag1)      ! out
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth_max_ym1_diag1)  ! out
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth_diag1)      ! out
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth_max_diag2)      ! out
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth_max_ym1_diag2)  ! out
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth_diag2)      ! out
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth_max_diag3)      ! out
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth_max_ym1_diag3)  ! out
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth_diag3)      ! out
    dsl4jsb_Get_var2D_onChunk(SSE_,   heat_cond)
    dsl4jsb_Get_var2D_onChunk(SSE_,   vol_heat_cap)
    dsl4jsb_Get_var2D_onChunk(SEB_,   t_unfilt)        ! in

    lstart = is_time_experiment_start(options%current_datetime)

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

    IF (.NOT. dsl4jsb_Config(SSE_)%l_snow) THEN
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
        & weq_snow_soil(:)                          ,  &
        & heat_cond(:)                              ,  &
        & vol_heat_cap(:)                           ,  &
        & t_soil_acoef                          (:,:), &
        & t_soil_bcoef                          (:,:), &
        & t_soil_sl                             (:,:), &
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
      IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
        !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_3 tile_name=tile%name
        snow_depth(:) = weq_snow_soil(:) * rhoh2o / snow_soil_dens(:)
      ELSE
        !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_4 tile_name=tile%name
        snow_depth(:) = weq_snow_soil(:) * rhoh2o / dens_snow
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
        !$claw acc loop gang vector
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
        & t_soil_sl        (:,:), &
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

      !$claw acc parallel loop default(none) present(snow_depth_sl,snow_fract,t_unfilt_corrected,t_soil_top)
      DO i = 1, nc
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
      snow_depth(:) = weq_snow_soil(:) * rhoh2o / dens_snow
      snow_fract(:) = 0._wp
      t_soil_top(:) = t_unfilt_corrected(:)
      grnd_hflx_snow(:) = 0._wp
      hcap_grnd_snow(:) = 0._wp
      !$claw end expand
    END IF

    IF (.NOT. dsl4jsb_Config(SSE_)%l_heat_lat) THEN
      heat_flx_lat(:,:) = 0._wp
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
        & dsl4jsb_Config(SSE_)%l_freeze            , &
        & dsl4jsb_Config(SSE_)%l_supercool         , &
        & dsl4jsb_Config(SSE_)%i_thaw              , &
        & t_soil_top                          (:  ), &
        & soil_depth_sl                       (:,:), &
        & vol_heat_cap_sl                     (:,:), &
        & heat_cond_sl                        (:,:), &
        & heat_flx_lat                        (:,:), &
        ! inout
        & t_soil_sl                           (:,:), &
        & t_soil_acoef                        (:,:), &
        & t_soil_bcoef                        (:,:), &
        ! out
        & hcap_grnd_soil                      (:  ), &
        & grnd_hflx_soil                      (:  ), &
        & thaw_depth                          (:  ),  &
        & thaw_depth_diag1                          (:  ),  &
        & thaw_depth_diag2                          (:  ),  &
        & thaw_depth_diag3                          (:  )   &
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
        & dsl4jsb_Config(SSE_)%l_freeze            , &
        & dsl4jsb_Config(SSE_)%l_supercool         , &
        & dsl4jsb_Config(SSE_)%i_thaw              , &
        & t_soil_top                          (:  ), &
        & soil_depth_sl                       (:,:), &
        & vol_heat_cap_sl                     (:,:), &
        & heat_cond_sl                        (:,:), &
        & heat_flx_lat                        (:,:), &
        ! inout
        & t_soil_sl                           (:,:), &
        & t_soil_acoef                        (:,:), &
        & t_soil_bcoef                        (:,:), &
        ! out
        & hcap_grnd_soil                      (:  ), &
        & grnd_hflx_soil                      (:  ), &
        & thaw_depth                          (:  ),  &
        & thaw_depth_diag1                          (:  ),  &
        & thaw_depth_diag2                          (:  ),  &
        & thaw_depth_diag3                          (:  ),  &
        ! optional, in
        & ws_max_sl                           (:,:), &
        & matrix_pot_sl                       (:,:), &
        & bclapp_sl                           (:,:), &
        ! optional, inout
        & wtr_soil_sl                         (:,:), &
        & ice_soil_sl                         (:,:), &
        ! optional, out
        & wtr_freeze_sl                       (:,:), &
        & ice_melt_sl                         (:,:), &
        & wtr_soil_pot_scool_sl               (:,:)  &
        & )
      !$claw acc exit data delete(ws_max_sl)
      DEALLOCATE(ws_max_sl)
    END IF

    !pdv > calculating maximum thaw_depths
    IF ((.NOT. lstart) .AND. is_newyear(options%current_datetime, options%dtime)) THEN
      thaw_depth_max_ym1(:) = thaw_depth_max(:)
      thaw_depth_max(:) = 0._wp
    ELSE
      thaw_depth_max(:) = MAX(thaw_depth_max(:),thaw_depth(:))
    END IF

    IF ((.NOT. lstart) .AND. is_newyear(options%current_datetime, options%dtime)) THEN
      thaw_depth_max_ym1_diag1(:) = thaw_depth_max_diag1(:)
      thaw_depth_max_diag1(:) = 0._wp
    ELSE
      thaw_depth_max_diag1(:) = MAX(thaw_depth_max_diag1(:),thaw_depth_diag1(:))
    END IF

    IF ((.NOT. lstart) .AND. is_newyear(options%current_datetime, options%dtime)) THEN
      thaw_depth_max_ym1_diag2(:) = thaw_depth_max_diag2(:)
      thaw_depth_max_diag2(:) = 0._wp
    ELSE
      thaw_depth_max_diag2(:) = MAX(thaw_depth_max_diag2(:),thaw_depth_diag2(:))
    END IF

    IF ((.NOT. lstart) .AND. is_newyear(options%current_datetime, options%dtime)) THEN
      thaw_depth_max_ym1_diag3(:) = thaw_depth_max_diag3(:)
      thaw_depth_max_diag3(:) = 0._wp
    ELSE
      thaw_depth_max_diag3(:) = MAX(thaw_depth_max_diag3(:),thaw_depth_diag3(:))
    END IF
    !pdv <

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
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_and_snow_temperature'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(SSE_,   t_soil_acoef,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   t_soil_bcoef,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   t_soil_sl,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   grnd_hflx,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   grnd_hflx_old,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   hcap_grnd,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   hcap_grnd_old,    weighted_by_fract)
    IF (dsl4jsb_Config(SSE_)%l_snow) THEN
      dsl4jsb_Aggregate_onChunk(SSE_,   snow_depth_sl,    weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_,   t_snow_acoef,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_,   t_snow_bcoef,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_,   t_snow,           weighted_by_fract)
    END IF
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_freeze_sl,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, ice_melt_sl,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_pot_scool_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth_max,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth_max_ym1,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth_diag1,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth_max_diag1,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth_max_ym1_diag1,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth_diag2,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth_max_diag2,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth_max_ym1_diag2,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth_diag3,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth_max_diag3,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth_max_ym1_diag3,    weighted_by_fract)
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_and_snow_temperature

  !-----------------------------------------------------------------------------------------------------
  !> update_soil_and_snow_properties
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_soil_and_snow_properties(tile, options)

    USE mo_jsb_physical_constants, ONLY: ci, dens_snow_min ! , cs

    USE mo_sse_constants,          ONLY: vol_hcap_ice, vol_hcap_snow, vol_hcap_org_top, vol_hcap_org_below, &
                                       & hcond_ice, hcond_snow, hcond_snow_min, &
                                       & hcond_org_top, hcond_dry_org_top, hcond_org_below, hcond_dry_org_below !, vol_hcap_snow_min
    USE mo_hydro_constants,        ONLY: vol_porosity_org_top, vol_porosity_org_below
    USE mo_sse_process,            ONLY: calc_vol_heat_capacity, calc_thermal_conductivity

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & snow_soil_dens,       &
      & vol_heat_cap,         &
      & heat_cond,            &
      & vol_heat_cap_snow,    &
      & heat_cond_snow
    dsl4jsb_Real3D_onChunk :: &
      & soil_depth_sl,        &
      & fract_org_sl,         &
      & vol_heat_cap_sl,      &
      & heat_cond_sl,         &
      & wtr_soil_sl,          &
      & ice_soil_sl,          &
      & vol_porosity_sl,      &
      & vol_field_cap_sl

    ! Locally allocated vectors
    !
    dsl4jsb_Real3D_onChunk :: &
      & fract_water,          &
      & fract_ice

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_e

    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: &
      & vol_heat_cap_tmp,                    &
      & vol_heat_cap_org_tmp,                &
      & soil_e_dz_tmp,                       &
      & vol_porosity_org_tmp,                &
      & heat_cond_tmp,                       &
      & hcond_org_tmp,                       &
      & hcond_dry_org_tmp

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

    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(HYDRO_)

    ! Get reference to variables for domain
    !
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var3D_onChunk(SSE_, vol_heat_cap_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(SSE_, heat_cond_sl)      ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, vol_heat_cap_snow) ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, heat_cond_snow) ! out


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

    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_var2D_onChunk(SSE_, vol_heat_cap)
    dsl4jsb_Get_var2D_onChunk(SSE_, heat_cond)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_field_cap_sl)

    IF (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil) THEN
      dsl4jsb_Get_var3D_onChunk(HYDRO_, fract_org_sl)
    ELSE
      ALLOCATE(fract_org_sl(nc, nsoil))
      !$claw acc enter data pcreate(fract_org_sl)
      !$claw expand parallel
      fract_org_sl(:,:) = 0._wp
    END IF

    IF (dsl4jsb_Config(SSE_)%l_heat_cap_dyn) THEN
      ALLOCATE(fract_water(nc, nsoil), fract_ice(nc, nsoil))
      !$claw acc enter data pcreate(fract_water, fract_ice)

      !$claw acc data present(soil_depth_sl,fract_water,fract_ice,wtr_soil_sl,ice_soil_sl)
      !$claw acc parallel loop gang vector collapse(2)
      DO i = 1, size(soil_depth_sl, dim=1)
        DO j = 1, size(soil_depth_sl, dim=2)
          IF (soil_depth_sl(i,j) > 0._wp) THEN
            fract_water(i,j) = wtr_soil_sl(i,j) / soil_depth_sl(i,j)
            fract_ice  (i,j) = ice_soil_sl(i,j) / soil_depth_sl(i,j)
          ELSE
            fract_water(i,j) = 0._wp
            fract_ice  (i,j) = 0._wp
          END IF
        END DO
      END DO
      !$claw acc end data

      ALLOCATE(vol_heat_cap_tmp(nc,nsoil), vol_heat_cap_org_tmp(nc,nsoil), soil_e_dz_tmp(nc,nsoil))
      !$claw acc enter data pcreate(vol_heat_cap_tmp,vol_heat_cap_org_tmp,soil_e_dz_tmp)
      !$claw acc update host(vol_heat_cap(:))
      vol_heat_cap_tmp(:,:) = SPREAD( vol_heat_cap(:), DIM=2, ncopies=nsoil)

      vol_heat_cap_org_tmp(:,1) = vol_hcap_org_top
      DO j = 2, nsoil
        vol_heat_cap_org_tmp(:,j) = vol_hcap_org_below
      END DO

      soil_e_dz_tmp(:,:) = SPREAD( soil_e%dz(:), DIM=1, ncopies=nc)
      !$claw acc update device(vol_heat_cap_tmp,vol_heat_cap_org_tmp,soil_e_dz_tmp)

      !$claw expand parallel update savepoint=ser_jsb_sse_calc_vol_heat_capacity_1 tile_name=tile%name
      vol_heat_cap_sl (:,:) = calc_vol_heat_capacity( &
        !&                        SPREAD( vol_heat_cap(:),                 DIM=2, ncopies=nsoil) , &
        &                                vol_heat_cap_tmp(:,:)     , &
        &                                vol_heat_cap_org_tmp(:,:) , &
        &                                vol_porosity_sl(:,:)      , &
        &                                fract_org_sl(:,:)         , &
        &                                fract_water(:,:)          , &
        &                                fract_ice  (:,:)          , &
        &                                soil_depth_sl(:,:)        , &
        !&                        SPREAD( soil_e%dz(:),                    DIM=1, ncopies=nc) &
        &                                soil_e_dz_tmp(:,:)          &
        & )

      !$claw acc exit data delete(fract_water,fract_ice,vol_heat_cap_tmp, vol_heat_cap_org_tmp,soil_e_dz_tmp)
      DEALLOCATE(fract_water, fract_ice, vol_heat_cap_tmp, vol_heat_cap_org_tmp, soil_e_dz_tmp)
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

    IF (dsl4jsb_Config(SSE_)%l_heat_cond_dyn) THEN

      ALLOCATE(heat_cond_tmp(nc,nsoil), vol_porosity_org_tmp(nc,nsoil), &
             & hcond_org_tmp(nc,nsoil), hcond_dry_org_tmp(nc,nsoil), soil_e_dz_tmp(nc,nsoil))
      !$claw acc enter data pcreate(heat_cond_tmp,vol_porosity_org_tmp,hcond_org_tmp,hcond_dry_org_tmp,soil_e_dz_tmp)
      !$claw acc update host(heat_cond(:))
      heat_cond_tmp(:,:) = SPREAD( heat_cond(:), DIM=2, ncopies=nsoil)

      vol_porosity_org_tmp(:,1)   = vol_porosity_org_top
      hcond_org_tmp(:,1)          = hcond_org_top
      hcond_dry_org_tmp(:,1)      = hcond_dry_org_top
      DO j = 2, nsoil
        vol_porosity_org_tmp(:,j) = vol_porosity_org_below
        hcond_org_tmp(:,j)        = hcond_org_below
        hcond_dry_org_tmp(:,j)    = hcond_dry_org_below
      END DO

      soil_e_dz_tmp(:,:) = SPREAD( soil_e%dz(:), DIM=1, ncopies=nc)
      !$claw acc update device(heat_cond_tmp,vol_porosity_org_tmp,hcond_org_tmp,hcond_dry_org_tmp,soil_e_dz_tmp)

      !$claw expand parallel update savepoint=ser_jsb_sse_calc_thermal_conductivity_1 tile_name=tile%name
      heat_cond_sl (:,:) = calc_thermal_conductivity( &
        !&                        SPREAD( heat_cond(:),                    DIM=2, ncopies=nsoil)             , &
        &                                heat_cond_tmp(:,:)                                                 , &
        &                                wtr_soil_sl(:,:)                                                   , &
        &                                ice_soil_sl(:,:)                                                   , &
        &                                vol_porosity_sl(:,:)                                               , &
        &                                vol_field_cap_sl(:,:)                                              , &
        &                                fract_org_sl(:,:)                                                  , &
        &                                vol_porosity_org_tmp(:,:)                                          , &
        &                                hcond_org_tmp(:,:)                                                 , &
        &                                hcond_dry_org_tmp(:,:)                                             , &
        &                                soil_depth_sl(:,:)                                                 , &
        &                                soil_e_dz_tmp(:,:)                                                   &
        !&                        SPREAD( soil_e%dz(:),                    DIM=1, ncopies=nc)                  &
        & )

      !$claw acc exit data delete(heat_cond_tmp,vol_porosity_org_tmp,hcond_org_tmp,hcond_dry_org_tmp,soil_e_dz_tmp)
      DEALLOCATE(heat_cond_tmp, vol_porosity_org_tmp,hcond_org_tmp, hcond_dry_org_tmp, soil_e_dz_tmp)
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

    IF (.NOT. (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil)) THEN
      !$claw acc exit data delete(fract_org_sl)
      DEALLOCATE(fract_org_sl)
    END IF

    IF (model%config%l_compat401) THEN
      !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_properties_expand_4 tile_name=tile%name
      vol_heat_cap_snow(:) = 634500._wp
      heat_cond_snow(:) = 0.31_wp
      !$claw end expand
    ELSE
      IF (dsl4jsb_Config(SSE_)%l_snow) THEN
        IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
          dsl4jsb_Get_var2D_onChunk(HYDRO_, snow_soil_dens)
          !$claw acc data present(snow_soil_dens,vol_heat_cap_snow,heat_cond_snow)
          !$claw acc parallel loop gang vector
          DO i = 1, nc
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

    dsl4jsb_Def_memory(SSE_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_and_snow_properties'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(SSE_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(SSE_, vol_heat_cap_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_, heat_cond_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_, vol_heat_cap_snow, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_, heat_cond_snow,    weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_and_snow_properties

#endif
END MODULE mo_sse_interface
