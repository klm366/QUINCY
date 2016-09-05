!>
!! @brief 
!!    generic routines providing various ways for the averaging of variables over time for all model processes
!!
!! @author
!!    Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!!    First version by Jan Engel (2016-09-01) - the basic routines were previously implemented in SLM mo_veg_memory
!!
MODULE mo_shared_update_time_averages

  USE mo_kind,                    ONLY: wp, eps8
  USE mo_exception,               ONLY: message, finish
  USE mo_jsb_process_class,       ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_time_averages


  CHARACTER(len=*), PARAMETER :: modname = 'mo_shared_update_time_averages'

CONTAINS 


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update time averages across processes
  !! 
  !! 
  !! 1.0 daytime averages \n
  !! 2.0 moving averages \n
  !! 3.0 ...
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_time_averages(tile, time_day, timestep_prop_day)

    USE mo_base_memory,             ONLY: t_slm_tile
    USE mo_atm2land_memory_class,   ONLY: t_atm2land_memory

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)        :: tile      !< one tile
    REAL(wp), INTENT(in)                  :: time_day           !< time in days
    REAL(wp), INTENT(in)                  :: timestep_prop_day  !< length of one time step (proportion of one day)
    ! 0.2 Local
    ! 0.3 Memory allocation
    TYPE(t_atm2land_memory), POINTER :: atm2land__mem   ! dsl4jsb_Def_memory(ATM2LAND_)
    REAL(wp), POINTER                :: swpar_srf_down  ! light at surface level
    REAL(wp), POINTER                :: t_air
    REAL(wp), POINTER                :: q_air
    REAL(wp), POINTER                :: press_srf
    REAL(wp), POINTER                :: co2_mixing_ratio
    REAL(wp), POINTER                :: ga

    SELECT TYPE (mm => tile%mem(ATM2LAND_)%p)            ! dsl4jsb_Get_memory(ATM2LAND_)
      TYPE IS(t_atm2land_memory)
        atm2land__mem => mm
    END SELECT

    ! atm2land
    swpar_srf_down    => atm2land__mem%swpar_srf_down%ptr
    t_air             => atm2land__mem%t_air%ptr
    q_air             => atm2land__mem%q_air%ptr
    press_srf         => atm2land__mem%press_srf%ptr
    co2_mixing_ratio  => atm2land__mem%co2_mixing_ratio%ptr
    ga                => atm2land__mem%ga%ptr



    !> 1.0 daytime averages
    !! 
    CALL calc_daytime_averages(tile%plant, time_day, timestep_prop_day, swpar_srf_down)


    !> 2.0 moving averages
    !! 
    CALL calc_moving_averages(tile%plant, t_air, q_air, press_srf, co2_mixing_ratio, ga)



    !> 3.0 .... averages < open for implementation
    !! 


  END SUBROUTINE update_time_averages


  !-----------------------------------------------------------------------------------------------------
  ! Sub-Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate daytime averages
  !!
  !! daytime: the time of the day over which light is available
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_daytime_averages(veg, time_day, timestep_prop_day, swpar_srf_down)

    USE mo_base_constants,      ONLY: ncanop
    USE mo_veg_memory,          ONLY: t_slm_plant

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_plant), INTENT(inout)      :: veg                !< veg variables per tile
    REAL(wp), INTENT(in)                  :: time_day           !< time in days
    REAL(wp), INTENT(in)                  :: timestep_prop_day  !< length of one time step (proportion of one day)
    REAL(wp), INTENT(in)                  :: swpar_srf_down     !< light at surface level
    ! 0.2 Local
    INTEGER                               :: i
    ! 0.3 Memory allocation


    !> 1.0 if light is available (i.e. daytime)
    !! 
    IF(swpar_srf_down > eps8) THEN
      ! 1.1 count number of values added to var_daytime
      veg%daytime_counter = veg%daytime_counter + 1.0_wp
      ! 1.2 add values and count number of values
      DO i=1,ncanop
        veg%daily_average%curr_val_par_sunlit_cl(i) = veg%daily_average%curr_val_par_sunlit_cl(i)     + veg%par_sunlit_cl(i)
        veg%daily_average%curr_val_par_shaded_cl(i) = veg%daily_average%curr_val_par_shaded_cl(i)     + veg%par_shaded_cl(i)
        veg%daily_average%curr_val_fleaf_sunlit_cl(i) = veg%daily_average%curr_val_fleaf_sunlit_cl(i) + veg%fleaf_sunlit_cl(i)
      END DO
      veg%daily_average%curr_val_t_jmax_opt         = veg%daily_average%curr_val_t_jmax_opt           + veg%t_jmax_opt
    END IF


    !> 2.0 at the end of the day (midnight), calc daytime averages (of the ending day)
    !! 
    !! calculate the average of the current day and pass this value to the according variable for the "previous day" value
    IF( ABS(((time_day + timestep_prop_day) - REAL(NINT(time_day)))) < eps8 ) THEN
      ! check if at least one value had been assigned at the current day, avoiding division by zero 
      ! calculate the average of the current day and pass it to the according avg_prev_* variable; 
      IF(veg%daytime_counter > eps8) THEN
        DO i=1,ncanop
          veg%daily_average%par_sunlit_cl(i)    = veg%daily_average%curr_val_par_sunlit_cl(i)   / veg%daytime_counter
          veg%daily_average%par_shaded_cl(i)    = veg%daily_average%curr_val_par_shaded_cl(i)   / veg%daytime_counter
          veg%daily_average%fleaf_sunlit_cl(i)  = veg%daily_average%curr_val_fleaf_sunlit_cl(i) / veg%daytime_counter
        END DO
        veg%daily_average%t_jmax_opt            = veg%daily_average%curr_val_t_jmax_opt         / veg%daytime_counter
      ELSE  ! if veg%daytime_counter is zero
        veg%daily_average%par_sunlit_cl(:)    = 0.0_wp
        veg%daily_average%par_shaded_cl(:)    = 0.0_wp
        veg%daily_average%fleaf_sunlit_cl(:)  = 0.0_wp
        veg%daily_average%t_jmax_opt          = 0.0_wp
      END IF

!      write(999,*) time_day, veg%daytime_counter, &   ! might be used for testing
!                  veg%daily_average%par_sunlit_cl(1), &
!                  veg%daily_average%par_shaded_cl(1), &
!                  veg%daily_average%fleaf_sunlit_cl(1), &
!                  veg%daily_average%t_jmax_opt, ' |-| '

      ! set all values of the current day to zero
      veg%daily_average%curr_val_par_sunlit_cl(:)      = 0.0_wp
      veg%daily_average%curr_val_par_shaded_cl(:)      = 0.0_wp
      veg%daily_average%curr_val_fleaf_sunlit_cl(:)    = 0.0_wp
      veg%daily_average%curr_val_t_jmax_opt            = 0.0_wp
      veg%daytime_counter                              = 0.0_wp
    END IF

  END SUBROUTINE calc_daytime_averages


  !-----------------------------------------------------------------------------------------------------
  !> calculate and return the moving average of a single variable or an array
  !!
  !! used in Sub-Task calc_moving_averages()
  !! 
  !! per default the avg_period_length has the unit 'day'; 
  !! 
  !! calculates the average across a certain period preceding to the current time step; \n 
  !! actually, the average gets rather estimated than exactly calculated using the weighted arithmetic mean, \n
  !! if the optional parameter 'do_calc' is present: the moving average is calculated if do_calc = TRUE !
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_mavg(current_avg, new_value, avg_period_length, do_calc, avg_period_unit) RESULT(new_avg)

    USE mo_base_constants,      ONLY: dt, one_day, one_year

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in)                      :: current_avg        !< current average value
    REAL(wp), INTENT(in)                      :: new_value          !< new value
    REAL(wp), INTENT(in)                      :: avg_period_length  !< length of the period over which the variable is averaged
    LOGICAL, OPTIONAL, INTENT(in)             :: do_calc            !< optional, do calculate new_avg from the new_value T/F
    CHARACTER(LEN=4), OPTIONAL, INTENT(in)    :: avg_period_unit    !< optional, for different units than 'day' for avg_period length
    REAL(wp)                                  :: new_avg            !< function output: the new average
    ! 0.2 Local
    REAL(wp)                                  :: new_value_weight   ! weight of the new value
    REAL(wp)                                  :: avg_period2seconds
    LOGICAL                                   :: l_do_calc          ! True/ False of the Optional variable do_calc present


    !> 1.0 optional arguments
    !! 

    !> 1.1 set TRUE per default, i.e. can only be FALSE if 'do_calc' is is present
    !!
    l_do_calc = .TRUE.
    ! check whether the optional statement is provided
    IF(PRESENT(do_calc))  l_do_calc = do_calc

    !> 1.2 calc period length in seconds
    !! 
    avg_period2seconds = one_day  ! (one_day == 86400._wp seconds)
    IF(l_do_calc .AND. PRESENT(avg_period_unit)) THEN
      SELECT CASE(TRIM(avg_period_unit))
      ! CASE('day') nothing to do as day is considered default anyway
      CASE('week')
        avg_period2seconds = avg_period2seconds * 7.0_wp
      CASE('year')
        avg_period2seconds = avg_period2seconds * one_year ! (one_year == 365.0_wp days)
!      CASE DEFAULT
!        CALL finish('mo_shared_update_time_averages','Invalid option for avg_period_unit argument to calc_mavg(). Use either day, week, or year.')
      END SELECT  ! select SLM model !
    END IF


    !> 2.0 calc weight of the new_value and the average
    !! calc the new average as weighted mean of the current_avg and the new_value weighted by new_value_weight
    IF(l_do_calc) THEN
      new_value_weight = 1.0_wp / (avg_period_length * avg_period2seconds / dt )  ! dt: length of one time step in seconds
      new_avg = (current_avg * (1.0_wp - new_value_weight)) + (new_value * new_value_weight)
    ELSE
      new_avg = current_avg
    END IF
  END FUNCTION calc_mavg


  !-----------------------------------------------------------------------------------------------------
  ! Sub-Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate moving averages
  !! 
  !! variable-specific averaging period, previous to the current timestep
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_moving_averages(veg, t_air, q_air, press_srf, co2_mixing_ratio, ga)

    USE mo_veg_memory,          ONLY: t_slm_plant, t_moving_average

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_plant), INTENT(inout), TARGET  :: veg    ! @slmtemp - remove TARGET and mvarg pointer later
    REAL(wp), INTENT(in)                :: t_air
    REAL(wp), INTENT(in)                :: q_air
    REAL(wp), INTENT(in)                :: press_srf
    REAL(wp), INTENT(in)                :: co2_mixing_ratio
    REAL(wp), INTENT(in)                :: ga
    ! 0.2 Local
    TYPE(t_moving_average), POINTER           :: mavg  ! @slmtemp
    REAl(wp)                                  :: hlp



    ! bi-directional pointer  @slmtemp
    mavg => veg%moving_average


    ! moving_avg(current average, new value, length of avg_period, do_calc=.FALSE., avg_period_unit='day')
    !           returns: new current average

    mavg%par_sunlit_cl(:)       = calc_mavg(mavg%par_sunlit_cl(:), veg%daily_average%par_sunlit_cl(:), mavg%p_optNfrac)
    mavg%par_shaded_cl(:)       = calc_mavg(mavg%par_shaded_cl(:), veg%daily_average%par_shaded_cl(:), mavg%p_optNfrac)

    mavg%fleaf_sunlit_cl(:)     = calc_mavg(mavg%fleaf_sunlit_cl(:), veg%daily_average%fleaf_sunlit_cl(:), mavg%p_optNfrac)
    mavg%beta_soil_gs           = calc_mavg(mavg%beta_soil_gs, veg%beta_soil_gs, mavg%p_optNfrac)
    mavg%beta_soil_ps           = calc_mavg(mavg%beta_soil_ps, veg%beta_soil_ps, mavg%p_optNfrac)
    mavg%beta_air               = calc_mavg(mavg%beta_air, veg%beta_air, mavg%p_optNfrac)

    mavg%t_jmax_opt             = calc_mavg(mavg%t_jmax_opt, veg%daily_average%t_jmax_opt, mavg%p_weekly)
    mavg%gpp_week               = calc_mavg(mavg%gpp_week, veg%flux%gpp, mavg%p_weekly)
    mavg%maint_respiration_week = calc_mavg(mavg%maint_respiration_week, veg%flux%maint_respiration_pot, mavg%p_weekly)

    mavg%growth_req_n_week      = calc_mavg(mavg%growth_req_n_week, veg%growth_req_n, mavg%p_weekly, &
                                            do_calc=(veg%growth_req_n > eps8))
    mavg%growth_req_p_week      = calc_mavg(mavg%growth_req_p_week, veg%growth_req_p, mavg%p_weekly, &
                                            do_calc=(veg%growth_req_p > eps8))
    mavg%boc_an                 = calc_mavg(mavg%boc_an, veg%boc_an, mavg%p_boc_an, &
                                            do_calc=(veg%pool%leaf%carbon > eps8))
    mavg%unit_npp               = calc_mavg(mavg%unit_npp, veg%flux%unit_npp, mavg%p_opt_alloc, & 
                                            do_calc=veg%growing_season)

    IF (veg%growth_req_n > 0._wp)THEN
       hlp = 1._wp/veg%growth_req_n
    ELSE
       hlp=0._wp
    ENDIF
    mavg%growth_cn_plant  = calc_mavg(mavg%growth_cn_plant, hlp, mavg%p_opt_alloc, do_calc=veg%growing_season)
    
    IF (veg%growth_req_n > 0._wp .AND. veg%growth_req_p > 0._wp)THEN
       hlp = 1._wp/(veg%growth_req_n * veg%growth_req_p)
    ELSE
       hlp=0._wp
    ENDIF
    mavg%growth_cp_plant  = calc_mavg(mavg%growth_cp_plant, hlp, mavg%p_opt_alloc, do_calc=veg%growing_season)

    mavg%allometry_npp                = calc_mavg(mavg%allometry_npp, veg%flux%npp, mavg%p_opt_alloc, &
                                                  do_calc=veg%growing_season)
    mavg%allometry_labile%nitrogen    = calc_mavg(mavg%allometry_labile%nitrogen, veg%pool%labile%nitrogen, mavg%p_opt_alloc, &
                                                  do_calc=veg%growing_season)
    mavg%allometry_labile%phosphorus  = calc_mavg(mavg%allometry_labile%phosphorus, veg%pool%labile%phosphorus, mavg%p_opt_alloc, &
                                                  do_calc=veg%growing_season)
    mavg%allometry_labile%carbon      = calc_mavg(mavg%allometry_labile%carbon, veg%pool%labile%carbon, mavg%p_opt_alloc, &
                                                  do_calc=veg%growing_season)
 
    mavg%demand_uptake_n              = calc_mavg(mavg%demand_uptake_n, veg%f_n_demand, mavg%p_demand_uptake_n, &
                                                  do_calc=(veg%pool%fine_root%carbon > eps8))
    mavg%demand_uptake_p              = calc_mavg(mavg%demand_uptake_p, veg%f_p_demand, mavg%p_demand_uptake_p, &
                                                  do_calc=(veg%pool%fine_root%carbon > eps8))
    mavg%unit_uptake_n                = calc_mavg(mavg%unit_uptake_n, veg%flux%unit_uptake_n, mavg%p_opt_alloc, &
                                                  do_calc=veg%growing_season)
    mavg%unit_uptake_p                = calc_mavg(mavg%unit_uptake_p, veg%flux%unit_uptake_p, mavg%p_opt_alloc, &
                                                  do_calc=veg%growing_season)

    mavg%leaf2root_mass_ratio         = calc_mavg(mavg%leaf2root_mass_ratio, veg%pool%leaf%carbon/veg%pool%fine_root%carbon, &
                                                  mavg%p_yearly, do_calc=(veg%pool%leaf%carbon > eps8))



    mavg%t_air            = calc_mavg(mavg%t_air,             t_air,            mavg%p_optNfrac)
    mavg%ga               = calc_mavg(mavg%ga,                ga,               mavg%p_optNfrac)
    mavg%co2_mixing_ratio = calc_mavg(mavg%co2_mixing_ratio,  co2_mixing_ratio, mavg%p_optNfrac)
    mavg%press_srf        = calc_mavg(mavg%press_srf,         press_srf,        mavg%p_optNfrac)
    mavg%q_air            = calc_mavg(mavg%q_air,             q_air,            mavg%p_optNfrac)

    mavg%t_air_week       = calc_mavg(mavg%t_air_week,        t_air, mavg%p_weekly)
    mavg%t_air_month      = calc_mavg(mavg%t_air_month,       t_air, mavg%p_monthly)
    mavg%t_air_longterm   = calc_mavg(mavg%t_air_longterm,    t_air, mavg%p_longterm)

  END SUBROUTINE calc_moving_averages


END MODULE mo_shared_update_time_averages
