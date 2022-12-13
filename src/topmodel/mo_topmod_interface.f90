!>
!! @brief Contains the interfaces to the soil and snow energy processes
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2022-12-06
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_top_interface
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

  dsl4jsb_Use_processes TOP_, SSE_, SEB_, HYDRO_, VEG_ 
  dsl4jsb_Use_config(TOP_)
  dsl4jsb_Use_memory(TOP_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_top_tasks

  !> Type definition for soil_temp_and_hydro task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_top
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_topmodel
    PROCEDURE, NOPASS :: Aggregate => aggregate_topmodel
  END TYPE tsk_top

  !> Constructor interface for topmodel task
  INTERFACE tsk_top
    PROCEDURE Create_task_topmodel
  END INTERFACE tsk_top

  CHARACTER(len=*), PARAMETER :: modname = 'mo_top_interface'

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for topmodel
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "topmodel"
  !!
  FUNCTION Create_task_topmodel(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_top::return_ptr)
    CALL return_ptr%Construct(name='topmodel', process_id=TOP_, owner_model_id=model_id)

  END FUNCTION Create_task_topmodel

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: TOP_ process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Register_top_tasks(process, model_id)

    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,              INTENT(in)    :: model_id

!    TYPE(t_jsb_model), POINTER :: model
    
!    model   => Get_model(model_id)

    CALL process%Register_task(tsk_top(model_id))

  END SUBROUTINE Register_top_tasks

  !-----------------------------------------------------------------------------------------------------
  !> update_topmodel
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_topmodel(tile, options)

!    USE mo_sse_process,            ONLY: calc_soil_temperature, calc_snow_temperature, calc_soiltemp_old
!    USE mo_jsb_physical_constants, ONLY: alf, rhoh2o, dens_snow
!    USE mo_sse_constants,          ONLY: snow_depth_min

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(TOP_)
    dsl4jsb_Def_memory(TOP_)
!    dsl4jsb_Def_memory(VEG_)

    ! Pointers to variables in memory
!    dsl4jsb_Real2D_onChunk :: t_unfilt

    ! Locally allocated vectors
    !
!    REAL(wp), DIMENSION(options%nc) :: &
!      & snow_depth,                    &
!    LOGICAL, DIMENSION(options%nc)  :: &
!      & ldglac
!    REAL(wp), ALLOCATABLE :: ws_max_sl(:,:)

!    INTEGER, DIMENSION(options%nc) ::  &
!      & itop, itop_old

    LOGICAL  :: lstart
    INTEGER  :: iblk, ics, ice, nc, i
    REAL(wp) :: dtime

!    TYPE(t_jsb_vgrid), POINTER :: soil_e, snow_e
!    TYPE(t_jsb_vgrid), POINTER :: vgrid_soil_sb                ! Vertical grid
!    INTEGER                    :: nsoil, nsnow
!    INTEGER                    :: nsoil_sb

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_topmodel'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (.NOT. tile%Is_process_calculated(TOP_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(TOP_)

!   soil_e => Get_vgrid('soil_depth_energy')
!    nsoil = soil_e%n_levels
!    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
!    nsoil_sb      = vgrid_soil_sb%n_levels

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(TOP_)
!    dsl4jsb_Get_memory(SSE_)

!    dsl4jsb_Get_var2D_onChunk(HYDRO_, snowmelt)

!    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)         ! inout

    lstart = is_time_experiment_start(options%current_datetime)

    !$claw acc enter data create(snow_depth,snow_fract,t_unfilt_corrected,t_soil_top,hcap_grnd_snow,hcap_grnd_soil,grnd_hflx_snow,grnd_hflx_soil,ldglac,itop,itop_old)

    ! Save old ground heat flux and heat capacity before updating
    ! The old flux and capacity (which have been used in the implicit calculation of the new surface temperature) are
    ! given back to the atmosphere so that the surface energy balance can be computed there
    !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_0 tile_name=tile%name
!    grnd_hflx_old(:) = grnd_hflx(:)
    !$claw end expand

    ! Correction of surface temperature for snow and ice melt
!    IF (lstart) THEN
!      !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_1 tile_name=tile%name
!      t_unfilt_corrected(:) = t_unfilt(:)
!    ELSE
!      !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_2 tile_name=tile%name
!      t_unfilt_corrected(:) = t_unfilt(:) - snowmelt(:) * dtime * alf / hcap_grnd(:)
!    END IF

!    IF (.NOT. dsl4jsb_Config(SSE_)%l_snow) THEN
!      ldglac(:) = SPREAD(tile%is_glacier, DIM=1, ncopies=nc)
!      !$claw acc update device(ldglac)

!      CALL calc_soiltemp_old(                          &
!        & nc,                                          &
!        & nsoil,                                       &
!        & soil_e%dz(:),                                &
!        & soil_e%levels(:),                            &
!        & dtime,                                       &
!        & lstart,                                      &
!        & t_unfilt_corrected(:)                     ,  &
!        & weq_snow_soil(:)                          ,  &
!        & heat_cond(:)                              ,  &
!        & vol_heat_cap(:)                           ,  &
!        & t_soil_acoef                          (:,:), &
!        & t_soil_bcoef                          (:,:), &
!        & t_soil_sl                             (:,:), &
!        & hcap_grnd                             (:),   &
!        & grnd_hflx                             (:),   &
!        & ldglac                                (:)    &
!        !& SPREAD(tile%is_glacier, DIM=1, ncopies=nc)   &
!        & )


      !$claw acc exit data delete(snow_depth,snow_fract,t_unfilt_corrected,t_soil_top,hcap_grnd_snow,hcap_grnd_soil,grnd_hflx_snow,grnd_hflx_soil,ldglac,itop,itop_old)  
!      RETURN

!    END IF

!      IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
        !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_3 tile_name=tile%name
!        snow_depth(:) = weq_snow_soil(:) * rhoh2o / snow_soil_dens(:)
!      ELSE
        !$claw expand parallel update savepoint=ser_jsb_sse_update_soil_and_snow_temperature_expand_4 tile_name=tile%name
!        snow_depth(:) = weq_snow_soil(:) * rhoh2o / dens_snow
!      END IF

      ! Calculate top snow layer for previous time step
      !$claw acc data present(itop_old,snow_depth_sl)
      !$claw acc parallel loop gang vector
!      DO i = 1,nc
!        itop_old(i) = nsnow + 1 - COUNT(snow_depth_sl(i,:) > 0._wp)
!      END DO
      !$claw acc end data

      !$claw acc data present(snow_depth,snow_depth_sl)
      !$claw acc parallel
      !$claw acc loop seq
!      DO isnow=nsnow,1,-1
        !$claw acc loop gang vector
!        DO i = 1, size(snow_depth, dim=1)
!          IF (snow_depth(i) > 0._wp) THEN
!            snow_depth_sl(i,isnow) = MIN(snow_depth(i), snow_e%dz(isnow))
!          END IF
!        END DO
!      END DO
      !$claw acc end parallel
      !$claw acc end data

      ! Calculate top snow layer for current time step

!      CALL calc_snow_temperature( &
        ! in
!        & nc,                     &
!        & nsnow,                  &
!        & soil_e%dz(1:2),         &
!        & dtime,                  &
!        & lstart,                 &
!        & snow_depth_sl    (:,:), &
!        & )

      !$claw acc parallel loop default(none) present(snow_depth_sl,snow_fract,t_unfilt_corrected,t_soil_top)
!      DO i = 1, nc
!        IF (snow_depth_sl(i, nsnow) > 0._wp) THEN  ! At least one snow layer above soil
!          t_soil_top(i) = snow_fract(i) * t_soil_top(i) + (1._wp - snow_fract(i)) * t_unfilt_corrected(i)
!        ELSE                                ! No snow
!          t_soil_top(i) = t_unfilt_corrected(i)
!        END IF
!      END DO
      !$claw end expand

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_topmodel

  SUBROUTINE aggregate_topmodel(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(TOP_)
    dsl4jsb_Def_memory(TOP_)
!    dsl4jsb_Def_memory(SSE_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_topmod'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(TOP_)
!    dsl4jsb_Get_memory(SSE_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

!    dsl4jsb_Aggregate_onChunk(SSE_,   t_soil_acoef,     weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_topmodel

#endif
END MODULE mo_top_interface
