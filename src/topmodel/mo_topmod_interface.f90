!>
!! @brief Contains the interfaces to the topmodel wetland diagnostics
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

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message, finish, message_text
  USE mo_jsb_control,            ONLY: debug_on
  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_grid_class,         ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,               ONLY: Get_vgrid, Get_grid
  USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,      ONLY: t_jsb_process, TOP_, A2L_, RAD_, SSE_, HYDRO_
  USE mo_jsb_task_class,         ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_math_constants,     ONLY: eps8
  USE mo_jsb_physical_constants, ONLY: tmelt

  dsl4jsb_Use_processes TOP_
  dsl4jsb_Use_config(TOP_)
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_memory(TOP_)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_top_tasks

  !> Type definition for topmodel task
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
    CALL return_ptr%Construct(name='inundated_area', process_id=TOP_, owner_model_id=model_id)

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

    USE mo_jsb_time,    ONLY: is_time_experiment_start, is_time_restart, is_newyear, is_newmonth, is_newday, &
      &                       timesteps_per_day, get_date_components, get_year_length, get_month_length
    USE mo_top_process, ONLY: update_topmodel_runmean, update_peat_runmean, update_water_table, update_inun_frac, &
      &                       update_aet, update_water_balance_peatland, get_climatological_inun

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(TOP_)
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(TOP_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    !
    ! Inputs
    dsl4jsb_Real2D_onChunk :: rain
    dsl4jsb_Real2D_onChunk :: snow
    dsl4jsb_Real2D_onChunk :: t_air
    dsl4jsb_Real2D_onChunk :: rad_srf_net
    dsl4jsb_Real2D_onChunk :: vol_field_cap
    dsl4jsb_Real2D_onChunk :: soil_depth
    dsl4jsb_Real3D_onChunk :: t_soil_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl
    dsl4jsb_Real3D_onChunk :: ice_soil_sl
    dsl4jsb_Real2D_onChunk :: fract_snow
    dsl4jsb_Real2D_onChunk :: snowmelt
    ! Calculated by topmodel
    dsl4jsb_Real2D_onChunk :: water_table
    dsl4jsb_Real2D_onChunk :: frac_liquid_top
    dsl4jsb_Real2D_onChunk :: frac_liquid_top_uncorr
    dsl4jsb_Real2D_onChunk :: precip_sum_day
    dsl4jsb_Real2D_onChunk :: precip_mean_day
    dsl4jsb_Real2D_onChunk :: frac_inun
    dsl4jsb_Real2D_onChunk :: frac_inun_uncorr
    dsl4jsb_Real2D_onChunk :: frac_inun_pre
    dsl4jsb_Real2D_onChunk :: frac_inun_sum_month
    dsl4jsb_Real2D_onChunk :: frac_inun_mean_month
    dsl4jsb_Real2D_onChunk :: frac_inun_min_year
    dsl4jsb_Real2D_onChunk :: frac_inun_min_lt
    dsl4jsb_Real2D_onChunk :: frac_inun_summ_year
    dsl4jsb_Real2D_onChunk :: frac_inun_summ_lt
    dsl4jsb_Real2D_onChunk :: frac_permanent_wetland
    dsl4jsb_Real2D_onChunk :: wt_inun
    dsl4jsb_Real2D_onChunk :: wt_inun_sum_month
    dsl4jsb_Real2D_onChunk :: wt_inun_mean_month
    dsl4jsb_Real2D_onChunk :: wt_inun_min_year
    dsl4jsb_Real2D_onChunk :: wt_inun_min_lt
    dsl4jsb_Real2D_onChunk :: aet
    dsl4jsb_Real2D_onChunk :: aet_sum_day
    dsl4jsb_Real2D_onChunk :: aet_mean_day
    dsl4jsb_Real2D_onChunk :: poaet_mean_day
    dsl4jsb_Real2D_onChunk :: poaet_sum_year
    dsl4jsb_Real2D_onChunk :: poaet_mean_year
    dsl4jsb_Real2D_onChunk :: poaet_lt
    dsl4jsb_Real2D_onChunk :: peat_water_table
    dsl4jsb_Real2D_onChunk :: peat_water_volume
    dsl4jsb_Real2D_onChunk :: peat_water_flux
    dsl4jsb_Real3D_onChunk :: frac_inun_clim

    ! Pointers to non-variables
    !
    REAL(wp), POINTER :: lat(:)
    REAL(wp), POINTER :: cti_mean(:)
    REAL(wp), POINTER :: cti_basin(:)
    REAL(wp), POINTER :: cti_openwater(:)
    REAL(wp), POINTER :: cti_phi(:)
    REAL(wp), POINTER :: cti_chi(:)
    REAL(wp), POINTER :: cti_mu(:)

    CLASS(t_jsb_tile_abstract), POINTER :: ptr_tile
    TYPE(t_jsb_grid), POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e, soil_w
    REAL(wp), ALLOCATABLE :: t_soil(:), wtr_soil(:), ice_soil(:)
    REAL(wp) :: water_table_t, frac_liquid_top_t
    LOGICAL  :: lstart, lrestart, lnew_year, lnew_month, lnew_day, lstop
    INTEGER  :: tspd, iblk, ics, ice, nc, nsoil, i
    INTEGER  :: year, month, summer_state, ylen, mlen
    REAL(wp) :: dtime, sw, si, swmx, simx, wtmx

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_topmodel'

    IF (.NOT. tile%Is_process_calculated(TOP_)) RETURN
    IF (.NOT. is_newday(options%current_datetime,options%dtime)) RETURN ! Topmodel process is updated only once per day

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    lstart     = is_time_experiment_start(options%current_datetime)
    lrestart   = is_time_restart(options%current_datetime)

    dtime = options%dtime
 
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    tspd  = timesteps_per_day(dtime)

    model  => Get_model(tile%owner_model_id)
    grid   => Get_grid(model%grid_id)
    soil_e => Get_vgrid('soil_depth_energy')
    soil_w => Get_vgrid('soil_depth_water')
    nsoil  = soil_e%n_levels
    lStop  = .FALSE.
    IF (nsoil /= soil_w%n_levels) lStop = .TRUE.
    IF (.NOT. lStop) lStop = ANY(ABS(soil_e%dz - soil_w%dz) > eps8)
    IF (lStop) CALL finish(routine,'Sorry, present implementation of topmodel requires water and energy soil levels to coincide')

    dsl4jsb_Get_Config(TOP_)
    dsl4jsb_Get_Config(SSE_)

    ! Allocate local work arrays
    !
    ALLOCATE(t_soil          (nsoil))
    ALLOCATE(wtr_soil        (nsoil))
    ALLOCATE(ice_soil        (nsoil))

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(TOP_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)

    ! A2L_ variables
    dsl4jsb_Get_var2D_onChunk(A2L_, rain)
    dsl4jsb_Get_var2D_onChunk(A2L_, snow)
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)
    ! RAD_ variables
    dsl4jsb_Get_var2D_onChunk(RAD_, rad_srf_net)
    ! SSE_ variables
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)
    ! HYDRO_ variables
    dsl4jsb_Get_var2D_onChunk(HYDRO_, soil_depth) ! StW: May be simplified by using soil_depth_sl instead
    dsl4jsb_Get_var2D_onChunk(HYDRO_, vol_field_cap)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snowmelt)
    ! TOP_ variables
    dsl4jsb_Get_var2D_onChunk(TOP_, frac_inun)
    dsl4jsb_Get_var2D_onChunk(TOP_, frac_inun_pre)
    dsl4jsb_Get_var2D_onChunk(TOP_, wt_inun)
    dsl4jsb_Get_var2D_onChunk(TOP_, aet)
    dsl4jsb_Get_var2D_onChunk(TOP_, frac_liquid_top)

    IF (lstart .AND. (dsl4jsb_Config(TOP_)%dbg_topmodel)) call message(routine,'First instance of update_topmodel')

    IF (dsl4jsb_Config(TOP_)%lpeat) THEN
      ! Variables only available if lpeat
      dsl4jsb_Get_var2D_onChunk(TOP_, precip_sum_day)
      dsl4jsb_Get_var2D_onChunk(TOP_, precip_mean_day)
      dsl4jsb_Get_var2D_onChunk(TOP_, aet_sum_day)
      dsl4jsb_Get_var2D_onChunk(TOP_, aet_mean_day)
      dsl4jsb_Get_var2D_onChunk(TOP_, poaet_mean_day)
      dsl4jsb_Get_var2D_onChunk(TOP_, poaet_sum_year)
      dsl4jsb_Get_var2D_onChunk(TOP_, poaet_mean_year)
      dsl4jsb_Get_var2D_onChunk(TOP_, poaet_lt)
      dsl4jsb_Get_var2D_onChunk(TOP_, peat_water_table)
      dsl4jsb_Get_var2D_onChunk(TOP_, peat_water_volume)
      dsl4jsb_Get_var2D_onChunk(TOP_, peat_water_flux)

      CALL update_aet(tmelt,t_air,rad_srf_net,aet)

! StW: Access to frac_liquid_top is questionable if inun_prescribed - likely seg-fault
      CALL update_water_balance_peatland(tspd, dtime,      & ! Time parameters
        & dsl4jsb_Config(TOP_)%peat_runoff_fac,            & ! Configuration parameters
        & dsl4jsb_Config(TOP_)%peat_runoff_exp,            &
        & rain(:)+snow(:), frac_liquid_top(:),             & ! Input variables from other parts of the model
        & fract_snow(:), snowmelt(:),                      &
        & aet(:),                                          &
        & peat_water_volume(:), peat_water_table(:),       & ! In/out variables
        & peat_water_flux(:)                               &
      )
    ENDIF

    IF (dsl4jsb_Config(TOP_)%inun_prescribed) THEN
      IF (is_newmonth(options%current_datetime,dtime)) THEN
        CALL get_date_components(options%current_datetime,MONTH=month)
        frac_inun_clim => dsl4jsb_memory(TOP_)%frac_inun_clim(ics:ice,iblk,:)
        CALL get_climatological_inun(frac_inun_clim(:,month),frac_inun(:),wt_inun(:))
      ENDIF
      IF (lstart) frac_inun_pre(:) = frac_inun(:) 
    ELSE

      ! Variables only available if not inun_prescribed
      dsl4jsb_Get_var2D_onChunk(TOP_, water_table)
      dsl4jsb_Get_var2D_onChunk(TOP_, frac_liquid_top_uncorr)
      dsl4jsb_Get_var2D_onChunk(TOP_, frac_inun_uncorr)
      dsl4jsb_Get_var2D_onChunk(TOP_, frac_inun_sum_month)
      dsl4jsb_Get_var2D_onChunk(TOP_, frac_inun_mean_month)
      dsl4jsb_Get_var2D_onChunk(TOP_, frac_inun_min_year)
      dsl4jsb_Get_var2D_onChunk(TOP_, frac_inun_min_lt)
      dsl4jsb_Get_var2D_onChunk(TOP_, frac_inun_summ_year)
      dsl4jsb_Get_var2D_onChunk(TOP_, frac_inun_summ_lt)
      dsl4jsb_Get_var2D_onChunk(TOP_, frac_permanent_wetland)
      dsl4jsb_Get_var2D_onChunk(TOP_, wt_inun_sum_month)
      dsl4jsb_Get_var2D_onChunk(TOP_, wt_inun_mean_month)
      dsl4jsb_Get_var2D_onChunk(TOP_, wt_inun_min_year)
      dsl4jsb_Get_var2D_onChunk(TOP_, wt_inun_min_lt)

      ! Referenced constants
      lat           =>                 grid%lat          (ics:ice,iblk)
      cti_mean      => dsl4jsb_memory(TOP_)%cti_mean     (ics:ice,iblk)
      cti_basin     => dsl4jsb_memory(TOP_)%cti_basin    (ics:ice,iblk)
      cti_openwater => dsl4jsb_memory(TOP_)%cti_openwater(ics:ice,iblk)
      cti_phi       => dsl4jsb_memory(TOP_)%cti_phi      (ics:ice,iblk)
      cti_chi       => dsl4jsb_memory(TOP_)%cti_chi      (ics:ice,iblk)
      cti_mu        => dsl4jsb_memory(TOP_)%cti_mu       (ics:ice,iblk)

      ! Save inundated fraction from last time step for the methane process
      IF (.NOT. lrestart) THEN
        frac_inun_pre(:) = frac_inun(:)
      ENDIF

      ! Calculate water table
      !
      DO i = 1, nc
        water_table      (i  ) = 0._wp
        frac_liquid_top  (i  ) = 0._wp
      ENDDO

      DO i = 1, nc
        ! Secure calculations are done on contiguous data - may be faster
        t_soil  (:) = t_soil_sl  (i,:)
        wtr_soil(:) = wtr_soil_sl(i,:)
        ice_soil(:) = ice_soil_sl(i,:)

        ! StW: Very misleading that ubounds/lbounds referes to the indicies and not to the physical direction
        CALL update_water_table(dsl4jsb_Config(SSE_)%l_freeze,tmelt,                            &
          &  soil_depth(i),dsl4jsb_Config(TOP_)%soil_depth_min,                                 &
          &  soil_e%dz(:),soil_e%ubounds(:),dsl4jsb_Config(TOP_)%saturation_threshold,          &
          &  vol_field_cap(i),                                                                  &
          &  t_soil(:),wtr_soil(:),ice_soil(:),                                                 &
          &  water_table_t, frac_liquid_top_t                                                   &
          &  )

        water_table    (i  ) = water_table_t
        frac_liquid_top(i  ) = frac_liquid_top_t
!        frac_liquid_sl (i,:) = frac_liquid_t_sl(:)
      ENDDO
      DO i=1, nc
        frac_liquid_top_uncorr(i) =            frac_liquid_top(i)
        frac_liquid_top       (i) = MAX(0.1_wp,frac_liquid_top(i))
      ENDDO

      ! Calculate inundated area
      CALL update_inun_frac(                         &
        & dsl4jsb_Config(TOP_)%soil_depth_min,       & ! Scalar parameters
        & dsl4jsb_Config(TOP_)%use_cti_basin,        &
        & dsl4jsb_Config(TOP_)%use_cti_flood,        &
        & dsl4jsb_Config(TOP_)%use_cti_max,          &
        & dsl4jsb_Config(TOP_)%use_wt_max,           &
        & dsl4jsb_Config(TOP_)%use_openwater,        &
        & dsl4jsb_Config(TOP_)%cti_min,              &
        & dsl4jsb_Config(TOP_)%cti_max,              &
        & dsl4jsb_Config(TOP_)%cti_min_flooding,     &
        & dsl4jsb_Config(TOP_)%efold,                &
        & dsl4jsb_Config(TOP_)%wt_min,               &
        & dsl4jsb_Config(TOP_)%wt_max,               &
        & soil_depth(:),                             & ! Geographical dependent parameters
        & cti_mean(:),cti_basin(:),cti_openwater(:), &
        & cti_phi(:),cti_mu(:),cti_chi(:),           &
        & water_table(:),frac_liquid_top(:),         & ! Variable inputs
        & frac_inun(:),frac_inun_uncorr(:),          & ! Outputs
        & wt_inun(:) &
      )

      IF (lstart) frac_inun_pre(:) = frac_inun(:) 

    ENDIF

    IF (dsl4jsb_Config(TOP_)%ltime_means) THEN

       ! Get needed time state
      lstart     = is_time_experiment_start(options%current_datetime)
      lnew_year  = is_newyear (options%current_datetime,dtime)
      lnew_month = is_newmonth(options%current_datetime,dtime)
      lnew_day   = is_newday  (options%current_datetime,dtime)
      summer_state = -1
      ylen = 365 ! Only need to be correct if lnew_year
      mlen = 30  ! Only need to be correct if lnew_month
      CALL get_date_components(options%current_datetime,YEAR=year,MONTH=month)
      IF (lnew_month) THEN
        ! Get the previous month and evt. year
        month = MOD(month-1,12)
        IF (month==0) THEN
          month = 12
          year = year - 1
        ENDIF
        IF ((month >= dsl4jsb_Config(TOP_)%month_summer_start) .AND. &
          & (month <  dsl4jsb_Config(TOP_)%month_summer_end  )) THEN
          ! Include month values in summer season sum
          summer_state = 0
        ELSEIF (month == dsl4jsb_Config(TOP_)%month_summer_end) THEN
          ! Include month values in summer season sum and compute the average
          summer_state = dsl4jsb_Config(TOP_)%month_summer_end - dsl4jsb_Config(TOP_)%month_summer_start + 1
        ENDIF
        mlen = get_month_length(year,month)
        IF (lnew_year) ylen = get_year_length(year)
      ENDIF

      IF (.NOT. dsl4jsb_Config(TOP_)%inun_prescribed) THEN
        CALL update_topmodel_runmean(lstart,dsl4jsb_Config(TOP_)%init_runmean,                      &
          &                          lnew_year,lnew_month,lnew_day,                                 &
          &                          mlen,ylen,tspd,summer_state,                                   &
          &                          REAL(dsl4jsb_Config(TOP_)%years_long_term,wp),                 &
          &                          lat(:)>dsl4jsb_Config(TOP_)%lat_const_wetland_type,            &
          &                          frac_inun(:), frac_inun_sum_month(:), frac_inun_mean_month(:), &
          &                          frac_inun_min_year(:), frac_inun_min_lt(:),                    &
          &                          frac_inun_summ_year(:), frac_inun_summ_lt(:),                  &
          &                          frac_permanent_wetland(:),                                     &
          &                          wt_inun(:), wt_inun_sum_month(:), wt_inun_mean_month(:),       &
          &                          wt_inun_min_year(:), wt_inun_min_lt(:)                         &
        )
      ENDIF

      IF (dsl4jsb_Config(TOP_)%lpeat) THEN

        CALL update_peat_runmean(lstart,dsl4jsb_Config(TOP_)%init_runmean,                    &
          &                      lnew_year,lnew_month,lnew_day,                               &
          &                      mlen,ylen,tspd,                                              &
          &                      REAL(dsl4jsb_Config(TOP_)%years_long_term,wp),               &
          &                      rain(:),snow(:),precip_sum_day(:),precip_mean_day(:),        &
          &                      aet(:),                                                      &
          &                      aet_sum_day(:),aet_mean_day(:),                              &
          &                      poaet_mean_day(:),                                           &
          &                      poaet_sum_year(:),poaet_mean_year(:),poaet_lt(:)             &
        )
      ENDIF

    ENDIF

    ! Deallocate local work arrays
    !
    DEALLOCATE(t_soil)
    DEALLOCATE(wtr_soil)
    DEALLOCATE(ice_soil)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_topmodel

  SUBROUTINE aggregate_topmodel(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(TOP_)
    dsl4jsb_Def_memory(TOP_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_topmod'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(TOP_)
    dsl4jsb_Get_memory(TOP_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(TOP_,  precip_mean_day       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TOP_,  water_table           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TOP_,  frac_liquid_top       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TOP_,  frac_liquid_top_uncorr, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun_uncorr      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TOP_,  wt_inun               , weighted_by_fract)
    IF (dsl4jsb_Config(TOP_)%ltime_means) THEN
      dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun_min_year,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun_summ_year,    weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun_mean_month,   weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun_min_lt,       weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun_summ_lt,      weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  wt_inun_min_year,       weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  wt_inun_mean_month,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  wt_inun_min_lt,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  frac_permanent_wetland, weighted_by_fract)
    ENDIF
    IF (dsl4jsb_Config(TOP_)%lpeat) THEN
      dsl4jsb_Aggregate_onChunk(TOP_,  aet,                    weighted_by_fract)
      IF (dsl4jsb_Config(TOP_)%ltime_means) THEN
        dsl4jsb_Aggregate_onChunk(TOP_,  precip_mean_day,        weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  aet_mean_day,           weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  poaet_mean_day,         weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  poaet_mean_year,        weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  poaet_lt,               weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  peat_water_flux,        weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  peat_water_volume,      weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  peat_water_table,       weighted_by_fract)
      ENDIF
    ENDIF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_topmodel

#endif
END MODULE mo_top_interface
