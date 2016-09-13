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
  USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element


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
  SUBROUTINE update_time_averages(tile, lctlib, time_day, timestep_prop_day)

    USE mo_base_memory,             ONLY: t_slm_tile
    USE mo_atm2land_memory_class,   ONLY: t_atm2land_memory

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)        :: tile      !< one tile
    TYPE(t_lctlib_element),INTENT(in)     :: lctlib    !< land-cover-type library - parameters for one pft
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
    CALL calc_moving_averages(lctlib, tile%plant, t_air, q_air, press_srf, co2_mixing_ratio, ga)



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
        veg%par_sunlit_daytime_dacc_cl(i)   = veg%par_sunlit_daytime_dacc_cl(i)     + veg%par_sunlit_cl(i)
        veg%par_shaded_daytime_dacc_cl(i)   = veg%par_shaded_daytime_dacc_cl(i)     + veg%par_shaded_cl(i)
        veg%fleaf_sunlit_daytime_dacc_cl(i) = veg%fleaf_sunlit_daytime_dacc_cl(i)   + veg%fleaf_sunlit_cl(i)
      END DO
      veg%t_jmax_opt_daytime_dacc           = veg%t_jmax_opt_daytime_dacc           + veg%t_jmax_opt
    END IF


    !> 2.0 at the end of the day (midnight), calc daytime averages (of the ending day)
    !! 
    !! calculate the average of the current day and pass this value to the according variable for the "previous day" value
    IF( ABS(((time_day + timestep_prop_day) - REAL(NINT(time_day)))) < eps8 ) THEN
      ! check if at least one value had been assigned at the current day, avoiding division by zero 
      ! calculate the average of the current day and pass it to the according avg_prev_* variable; 
      IF(veg%daytime_counter > eps8) THEN
        DO i=1,ncanop
          veg%par_sunlit_daytime_cl(i)    = veg%par_sunlit_daytime_dacc_cl(i)   / veg%daytime_counter
          veg%par_shaded_daytime_cl(i)    = veg%par_shaded_daytime_dacc_cl(i)   / veg%daytime_counter
          veg%fleaf_sunlit_daytime_cl(i)  = veg%fleaf_sunlit_daytime_dacc_cl(i) / veg%daytime_counter
        END DO
        veg%t_jmax_opt_daytime            = veg%t_jmax_opt_daytime_dacc         / veg%daytime_counter
      ELSE  ! if veg%daytime_counter is zero
        veg%par_sunlit_daytime_cl(:)    = 0.0_wp
        veg%par_shaded_daytime_cl(:)    = 0.0_wp
        veg%fleaf_sunlit_daytime_cl(:)  = 0.0_wp
        veg%t_jmax_opt_daytime          = 0.0_wp
      END IF

!      write(999,*) time_day, veg%daytime_counter, &   ! might be used for testing
!                  veg%par_sunlit_daytime_cl(1), &
!                  veg%par_shaded_daytime_cl(1), &
!                  veg%fleaf_sunlit_daytime_cl(1), &
!                  veg%t_jmax_opt_daytime, ' |-| '

      ! set all values of the current day to zero
      veg%par_sunlit_daytime_dacc_cl(:)      = 0.0_wp
      veg%par_shaded_daytime_dacc_cl(:)      = 0.0_wp
      veg%fleaf_sunlit_daytime_dacc_cl(:)    = 0.0_wp
      veg%t_jmax_opt_daytime_dacc            = 0.0_wp
      veg%daytime_counter                    = 0.0_wp
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
  !!
  !! avg_period_unit can be either day (which is the default anyway), week or year
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
    REAL(wp)                                  :: avg_period2seconds ! for converting the period length (day, week, year) into seconds
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
  SUBROUTINE calc_moving_averages(lctlib, veg, t_air, q_air, press_srf, co2_mixing_ratio, ga)

    USE mo_shared_constants     ! moving average period lengths are defined there as REAL(wp),PARAMETER :: mavg_period_...
    USE mo_veg_memory,          ONLY: t_slm_plant

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),INTENT(in)   :: lctlib             !< land-cover-type library - parameters for one pft
    TYPE(t_slm_plant), INTENT(inout)    :: veg                !< vegetation memory
    REAL(wp), INTENT(in)                :: t_air              !< air temperature ?
    REAL(wp), INTENT(in)                :: q_air              !< air .. ?
    REAL(wp), INTENT(in)                :: press_srf          !< air pressure at surface level ?
    REAL(wp), INTENT(in)                :: co2_mixing_ratio   !< .. ?
    REAL(wp), INTENT(in)                :: ga                 !< aerodynamic resistance
    ! 0.2 Local
    REAl(wp)                            :: hlp


    ! moving_avg(current average, new value, length of avg_period, do_calc=.FALSE., avg_period_unit='day')
    !           returns: new current average

    veg%par_sunlit_mavg_cl(:)       = calc_mavg(veg%par_sunlit_mavg_cl(:), veg%par_sunlit_daytime_cl(:), mavg_period_opt_Nfrac)
    veg%par_shaded_mavg_cl(:)       = calc_mavg(veg%par_shaded_mavg_cl(:), veg%par_shaded_daytime_cl(:), mavg_period_opt_Nfrac)

    veg%fleaf_sunlit_mavg_cl(:)     = calc_mavg(veg%fleaf_sunlit_mavg_cl(:), veg%fleaf_sunlit_daytime_cl(:), mavg_period_opt_Nfrac)
    veg%beta_soil_gs_mavg           = calc_mavg(veg%beta_soil_gs_mavg, veg%beta_soil_gs, mavg_period_opt_Nfrac)
    veg%beta_soil_ps_mavg           = calc_mavg(veg%beta_soil_ps_mavg, veg%beta_soil_ps, mavg_period_opt_Nfrac)
    veg%beta_air_mavg               = calc_mavg(veg%beta_air_mavg, veg%beta_air, mavg_period_opt_Nfrac)

    veg%t_jmax_opt_mavg             = calc_mavg(veg%t_jmax_opt_mavg, veg%t_jmax_opt_daytime, mavg_period_weekly)
    veg%gpp_mavg                    = calc_mavg(veg%gpp_mavg, veg%gpp, mavg_period_weekly)
    veg%an_mavg                     = calc_mavg(veg%an_mavg, veg%an, mavg_period_weekly)
    veg%npp_mavg                    = calc_mavg(veg%npp_mavg, veg%npp, mavg_period_weekly)
    veg%maint_respiration_mavg      = calc_mavg(veg%maint_respiration_mavg, veg%maint_respiration_pot, mavg_period_weekly)
    veg%fmaint_rate_mavg            = calc_mavg(veg%fmaint_rate_mavg, veg%fmaint_rate, mavg_period_weekly)

    veg%growth_req_n_mavg           = calc_mavg(veg%growth_req_n_mavg, veg%growth_req_n, mavg_period_weekly, &
                                            do_calc=(veg%growth_req_n > eps8))
    veg%growth_req_p_mavg           = calc_mavg(veg%growth_req_p_mavg, veg%growth_req_p, mavg_period_weekly, &
                                            do_calc=(veg%growth_req_p > eps8))
    veg%boc_an_mavg                 = calc_mavg(veg%boc_an_mavg, veg%boc_an, mavg_period_boc_an, &
                                            do_calc=(veg%pool%leaf%carbon > eps8))
    veg%unit_npp_mavg               = calc_mavg(veg%unit_npp_mavg, veg%unit_npp, mavg_period_opt_alloc, & 
                                            do_calc=veg%growing_season)

    IF (veg%growth_req_n > 0._wp)THEN
       hlp = 1._wp/veg%growth_req_n
    ELSE
       hlp=0._wp
    ENDIF
    veg%growth_cn_plant_mavg  = calc_mavg(veg%growth_cn_plant_mavg, hlp, mavg_period_opt_alloc, do_calc=veg%growing_season)
    
    IF (veg%growth_req_n > 0._wp .AND. veg%growth_req_p > 0._wp)THEN
       hlp = 1._wp/(veg%growth_req_n * veg%growth_req_p)
    ELSE
       hlp=0._wp
    ENDIF
    veg%growth_cp_plant_mavg  = calc_mavg(veg%growth_cp_plant_mavg, hlp, mavg_period_opt_alloc, do_calc=veg%growing_season)

    veg%allometry_npp_mavg                = calc_mavg(veg%allometry_npp_mavg, veg%npp, mavg_period_opt_alloc, &
                                                  do_calc=veg%growing_season)
    veg%allometry_labile_mavg%nitrogen    = calc_mavg(veg%allometry_labile_mavg%nitrogen, veg%pool%labile%nitrogen, & 
                                                  mavg_period_opt_alloc, do_calc=veg%growing_season)
    veg%allometry_labile_mavg%phosphorus  = calc_mavg(veg%allometry_labile_mavg%phosphorus, veg%pool%labile%phosphorus, & 
                                                  mavg_period_opt_alloc, do_calc=veg%growing_season)
    veg%allometry_labile_mavg%carbon      = calc_mavg(veg%allometry_labile_mavg%carbon, veg%pool%labile%carbon, & 
                                                  mavg_period_opt_alloc, do_calc=veg%growing_season)
 
    veg%demand_uptake_n_mavg              = calc_mavg(veg%demand_uptake_n_mavg, veg%f_n_demand, mavg_period_demand_uptake_n, &
                                                  do_calc=(veg%pool%fine_root%carbon > eps8))
    veg%demand_uptake_p_mavg              = calc_mavg(veg%demand_uptake_p_mavg, veg%f_p_demand, mavg_period_demand_uptake_p, &
                                                  do_calc=(veg%pool%fine_root%carbon > eps8))
    veg%unit_uptake_n_talloc_mavg         = calc_mavg(veg%unit_uptake_n_talloc_mavg, veg%unit_uptake_n, mavg_period_opt_alloc, &
                                                  do_calc=veg%growing_season)
    veg%unit_uptake_p_talloc_mavg         = calc_mavg(veg%unit_uptake_p_talloc_mavg, veg%unit_uptake_p, mavg_period_opt_alloc, &
                                                  do_calc=veg%growing_season)
    veg%unit_uptake_n_troot_mavg          = calc_mavg(veg%unit_uptake_n_troot_mavg, veg%unit_uptake_n, lctlib%tau_fine_root, &
                                                  do_calc=veg%growing_season, avg_period_unit='year')
    veg%unit_uptake_p_troot_mavg          = calc_mavg(veg%unit_uptake_p_troot_mavg, veg%unit_uptake_p, lctlib%tau_fine_root, &
                                                  do_calc=veg%growing_season, avg_period_unit='year')
!@todo should the average period be tau_root?
    veg%cost_n_uptake_root_mavg           = veg%cost_n_uptake_root
    veg%cost_n_fixation_mavg              = veg%cost_n_fixation
!    veg%cost_n_uptake_root_mavg           = calc_mavg(veg%cost_n_uptake_root_mavg, veg%cost_n_uptake_root, mavg_period_yearly, &
!                                                  do_calc=veg%growing_season) 
!    veg%cost_n_fixation_mavg              = calc_mavg(veg%cost_n_fixation_mavg, veg%cost_n_fixation, mavg_period_yearly, &
!                                                  do_calc=veg%growing_season) 

    veg%leaf2root_mass_ratio_mavg         = calc_mavg(veg%leaf2root_mass_ratio_mavg, &
                                                  veg%pool%leaf%carbon/veg%pool%fine_root%carbon, &
                                                  mavg_period_yearly, do_calc=(veg%pool%leaf%carbon > eps8))



    veg%t_air_mavg            = calc_mavg(veg%t_air_mavg,             t_air,            mavg_period_opt_Nfrac)
    veg%ga_mavg               = calc_mavg(veg%ga_mavg,                ga,               mavg_period_opt_Nfrac)
    veg%co2_mixing_ratio_mavg = calc_mavg(veg%co2_mixing_ratio_mavg,  co2_mixing_ratio, mavg_period_opt_Nfrac)
    veg%press_srf_mavg        = calc_mavg(veg%press_srf_mavg,         press_srf,        mavg_period_opt_Nfrac)
    veg%q_air_mavg            = calc_mavg(veg%q_air_mavg,             q_air,            mavg_period_opt_Nfrac)

    veg%t_air_week_mavg       = calc_mavg(veg%t_air_week_mavg,        t_air, mavg_period_weekly)
    veg%t_air_month_mavg      = calc_mavg(veg%t_air_month_mavg,       t_air, mavg_period_monthly)
    veg%t_air_longterm_mavg   = calc_mavg(veg%t_air_longterm_mavg,    t_air, mavg_period_longterm)

  END SUBROUTINE calc_moving_averages


END MODULE mo_shared_update_time_averages
