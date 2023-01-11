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
use mo_util, only: real2string, int2string, logical2string

  dsl4jsb_Use_processes TOP_
  dsl4jsb_Use_config(TOP_)
  dsl4jsb_Use_memory(TOP_)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(HYDRO_)

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

    USE mo_jsb_time,    ONLY: is_time_experiment_start, is_newyear, is_newmonth, is_newday, timesteps_per_day, &
      &                       get_date_components, get_year_length, get_month_length
    USE mo_top_process, ONLY: update_topmodel_runmean, update_peat_runmean, update_water_table, update_inun_frac, &
      &                       update_aet, update_water_balance_peatland
    USE mo_jsb_physical_constants, ONLY: stbo

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(TOP_)
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
    dsl4jsb_Real3D_onChunk :: frac_liquid_sl
    dsl4jsb_Real2D_onChunk :: precip_sum_day
    dsl4jsb_Real2D_onChunk :: precip_mean_day
    dsl4jsb_Real2D_onChunk :: frac_inun
    dsl4jsb_Real2D_onChunk :: frac_inun_uncorr
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
    dsl4jsb_Real2D_onChunk :: peat_water_balance

    ! Pointers to non-variables
    !
    REAL(wp), POINTER :: lat(:)
    REAL(wp), POINTER :: cti_mean(:)
    REAL(wp), POINTER :: cti_basin(:)
    REAL(wp), POINTER :: cti_openwater(:)
    REAL(wp), POINTER :: cti_phi(:)
    REAL(wp), POINTER :: cti_chi(:)
    REAL(wp), POINTER :: cti_mu(:)

    ! Locally allocated vectors
    !
!    REAL(wp), DIMENSION(options%nc) :: &
!      & cover_fract, cover_fract_sum

    CLASS(t_jsb_tile_abstract), POINTER :: ptr_tile
    TYPE(t_jsb_grid), POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e, soil_w
    REAL(wp), ALLOCATABLE :: t_soil(:), wtr_soil(:), ice_soil(:), frac_liquid_t_sl(:)
    REAL(wp) :: water_table_t, frac_liquid_top_t
    LOGICAL  :: lstart, lnew_year, lnew_month, lnew_day, lstop
    INTEGER  :: tspd, iblk, ics, ice, nc, nsoil, i
    INTEGER  :: year, month, summer_state, ylen, mlen
    REAL(wp) :: dtime
    INTEGER, SAVE :: n = 0
integer :: loc(1)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_topmodel'

    IF (.NOT. tile%Is_process_calculated(TOP_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime
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

    dsl4jsb_Get_config(TOP_)

    ! Allocate local work arrays
    !
    ALLOCATE(t_soil          (nsoil))
    ALLOCATE(wtr_soil        (nsoil))
    ALLOCATE(ice_soil        (nsoil))
    ALLOCATE(frac_liquid_t_sl(nsoil))

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
      dsl4jsb_Get_var2D_onChunk(TOP_, peat_water_balance)

      CALL update_aet(tmelt,t_air,rad_srf_net-stbo*t_air**4._wp,aet) ! StW: Radiation likely wrong

! StW: Access to frac_liquid_top is questionable if inun_prescribed - likely seg-fault
      CALL update_water_balance_peatland(tspd, dtime,      & ! Time parameters
        & dsl4jsb_Config(TOP_)%peat_runoff_fac,            & ! Configuration parameters
        & dsl4jsb_Config(TOP_)%peat_runoff_exp,            &
        & rain(:)+snow(:), frac_liquid_top(:),             & ! Input variables from other parts of the model
        & fract_snow(:), snowmelt(:),                      &
        & aet(:),                                          &
        & peat_water_volume(:), peat_water_table(:),       & ! In/out variables
        & peat_water_balance(:)                            &
      )
    ENDIF

    IF (dsl4jsb_Config(TOP_)%inun_prescribed) THEN
    ELSE

      ! Variables only available if not inun_prescribed
      dsl4jsb_Get_var2D_onChunk(TOP_, water_table)
      dsl4jsb_Get_var2D_onChunk(TOP_, frac_liquid_top_uncorr)
      dsl4jsb_Get_var3D_onChunk(TOP_, frac_liquid_sl)
      dsl4jsb_Get_var2D_onChunk(TOP_, frac_inun)
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

      ! Calculate water table
      !
      DO i = 1, nc
        water_table      (i  ) = 0._wp
        frac_liquid_top  (i  ) = 0._wp
        frac_liquid_sl   (i,:) = 0._wp
      ENDDO
      DO i = 1, nc
        ! Secure calculations are done on contiguous data - may be faster
        t_soil  (:) = t_soil_sl  (i,:)
        wtr_soil(:) = wtr_soil_sl(i,:)
        ice_soil(:) = ice_soil_sl(i,:)

        CALL update_water_table(.TRUE.,tmelt,soil_depth(i),dsl4jsb_Config(TOP_)%soil_depth_min, &
          &  soil_e%dz(:),soil_e%lbounds(:),dsl4jsb_Config(TOP_)%saturation_threshold,          &
          &  vol_field_cap(i),                                                                  &
          &  t_soil(:),wtr_soil(:),ice_soil(:),                                                 &
          &  water_table_t, frac_liquid_top_t, frac_liquid_t_sl(:)                              &
          &  )

        water_table    (i  ) = water_table_t
        frac_liquid_top(i  ) = frac_liquid_top_t
        frac_liquid_sl (i,:) = frac_liquid_t_sl(:)
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
!call message(routine,'Blk, inun_max: '//int2string(iblk)//' '//real2string(maxval(frac_inun)))

if (any(frac_inun > 1._wp)) then
  loc = maxloc(frac_inun)
  call message(routine,'Large inun, blk: '//int2string(iblk))
  call message('A:',real2string(cti_mean(loc(1))))
  call message('B:',real2string(cti_basin(loc(1))))
  call message('C:',real2string(cti_openwater(loc(1))))
  call message('D:',real2string(cti_phi(loc(1))))
  call message('E:',real2string(cti_chi(loc(1))))
  call message('F:',real2string(cti_mu(loc(1))))
  call message('G:',real2string(water_table(loc(1))))
  call message('H:',real2string(frac_liquid_top(loc(1))))
endif

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

      ! debug
      IF (lnew_month) then
        call message(routine,'ndom: '//int2string(mlen))
      endif
      IF (lnew_year)  then
        call message(routine,'ndoy: '//int2string(ylen))
      endif
    ENDIF


    if (n<24) then
            call message(routine,'blck: '//int2string(iblk)//' '//int2string(ics)//' '//int2string(ice)//' '//int2string(nc))
            write(message_text,*) shape(rain)
            call message(routine,trim(message_text))
            call message(routine,'test: '//int2string(n)//' '//logical2string(lnew_day)//' '// &
            real2string(maxval(rain+snow))//' '//real2string(maxval(precip_sum_day))//' '//real2string(maxval(precip_mean_day)))
    endif
    n = n + 1
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
    dsl4jsb_Aggregate_onChunk(TOP_,  frac_liquid_sl        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun_uncorr      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TOP_,  wt_inun               , weighted_by_fract)
    IF (dsl4jsb_Config(TOP_)%ltime_means) THEN
      dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun_min_year,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun_summ_year,    weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun_sum_month,    weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun_mean_month,   weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun_min_lt,       weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  frac_inun_summ_lt,      weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  wt_inun_min_year,       weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  wt_inun_sum_month,      weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  wt_inun_mean_month,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  wt_inun_min_lt,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(TOP_,  frac_permanent_wetland, weighted_by_fract)
    ENDIF
    IF (dsl4jsb_Config(TOP_)%lpeat) THEN
      dsl4jsb_Aggregate_onChunk(TOP_,  aet,                    weighted_by_fract)
      IF (dsl4jsb_Config(TOP_)%ltime_means) THEN
        dsl4jsb_Aggregate_onChunk(TOP_,  precip_sum_day,         weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  precip_mean_day,        weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  aet_sum_day,            weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  aet_mean_day,           weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  poaet_mean_day,         weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  poaet_sum_year,         weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  poaet_mean_year,        weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  poaet_lt,               weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  peat_water_balance,     weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  peat_water_volume,      weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(TOP_,  peat_water_table,       weighted_by_fract)
      ENDIF
    ENDIF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_topmodel

#endif
END MODULE mo_top_interface
