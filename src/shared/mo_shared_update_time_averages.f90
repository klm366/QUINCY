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

#include "dsl_jsbach.inc"

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

    ! Use of process memories
    dsl4jsb_Use_memory(ATM2LAND_)   !  USE mo_atm2land_memory_class, ONLY: t_atm2land_memory

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)        :: tile      !< one tile
    TYPE(t_lctlib_element),INTENT(in)     :: lctlib    !< land-cover-type library - parameters for one pft
    REAL(wp), INTENT(in)                  :: time_day           !< time in days
    REAL(wp), INTENT(in)                  :: timestep_prop_day  !< length of one time step (proportion of one day)
    ! 0.2 Local
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ATM2LAND_)                   !     TYPE(t_atm2land_memory), POINTER :: atm2land__mem
    REAL(wp), POINTER                :: swpar_srf_down  ! light at surface level
    REAL(wp), POINTER                :: t_air
    REAL(wp), POINTER                :: q_air
    REAL(wp), POINTER                :: press_srf
    REAL(wp), POINTER                :: co2_mixing_ratio
    REAL(wp), POINTER                :: ga
    dsl4jsb_Get_memory(ATM2LAND_)
!    SELECT TYPE (mm => tile%mem(ATM2LAND_)%p) 
!      TYPE IS(t_atm2land_memory)
!        atm2land__mem => mm
!    END SELECT
    swpar_srf_down    => atm2land__mem%swpar_srf_down%ptr
    t_air             => atm2land__mem%t_air%ptr
    q_air             => atm2land__mem%q_air%ptr
    press_srf         => atm2land__mem%press_srf%ptr
    co2_mixing_ratio  => atm2land__mem%co2_mixing_ratio%ptr
    ga                => atm2land__mem%ga%ptr


    !> 1.0 daytime averages
    !! 
    CALL calc_daytime_averages(tile%veg, time_day, timestep_prop_day, swpar_srf_down, & 
                               t_air, press_srf, co2_mixing_ratio, ga)


    !> 2.0 moving averages
    !! 
    CALL calc_moving_averages(lctlib, tile%veg, t_air, press_srf, co2_mixing_ratio, ga)



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
  SUBROUTINE calc_daytime_averages(veg, time_day, timestep_prop_day, swpar_srf_down, t_air, & 
                                   press_srf, co2_mixing_ratio, ga)

    USE mo_base_constants,      ONLY: ncanop
    USE mo_veg_memory,          ONLY: t_slm_plant

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_plant), INTENT(inout)      :: veg                !< veg static memory per tile
    REAL(wp), INTENT(in)                  :: time_day           !< time in days
    REAL(wp), INTENT(in)                  :: timestep_prop_day  !< length of one time step (proportion of one day)
    REAL(wp), INTENT(in)                  :: swpar_srf_down     !< shortwave downward flux in PAR band 
    REAL(wp), INTENT(in)                  :: t_air              !< air temperature (K)
    REAL(wp), INTENT(in)                  :: press_srf          !< surface pressure (Pa)
    REAL(wp), INTENT(in)                  :: co2_mixing_ratio   !< atmospheric CO2 concentration (ppm)
    REAL(wp), INTENT(in)                  :: ga                 !< aerodynamic conductance 
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
      veg%t_air_daytime_dacc            = veg%t_air_daytime_dacc            + t_air
      veg%press_srf_daytime_dacc        = veg%press_srf_daytime_dacc        + press_srf
      veg%co2_mixing_ratio_daytime_dacc = veg%co2_mixing_ratio_daytime_dacc + co2_mixing_ratio
      veg%ga_daytime_dacc               = veg%ga_daytime_dacc               + ga
      veg%beta_air_daytime_dacc         = veg%beta_air_daytime_dacc         + veg%beta_air
      veg%beta_soil_ps_daytime_dacc     = veg%beta_soil_ps_daytime_dacc     + veg%beta_soil_ps
      veg%beta_soil_gs_daytime_dacc     = veg%beta_soil_gs_daytime_dacc     + veg%beta_soil_gs
      veg%t_jmax_opt_daytime_dacc       = veg%t_jmax_opt_daytime_dacc       + veg%t_jmax_opt
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
        veg%t_air_daytime            = veg%t_air_daytime_dacc              / veg%daytime_counter
        veg%press_srf_daytime        = veg%press_srf_daytime_dacc          / veg%daytime_counter
        veg%co2_mixing_ratio_daytime = veg%co2_mixing_ratio_daytime_dacc   / veg%daytime_counter
        veg%ga_daytime               = veg%ga_daytime_dacc                 / veg%daytime_counter
        veg%beta_air_daytime         = veg%beta_air_daytime_dacc           / veg%daytime_counter
        veg%beta_soil_ps_daytime     = veg%beta_soil_ps_daytime_dacc       / veg%daytime_counter
        veg%beta_soil_gs_daytime     = veg%beta_soil_gs_daytime_dacc       / veg%daytime_counter
        veg%t_jmax_opt_daytime       = veg%t_jmax_opt_daytime_dacc         / veg%daytime_counter
      ELSE  ! if veg%daytime_counter is zero
        veg%par_sunlit_daytime_cl(:)    = 0.0_wp
        veg%par_shaded_daytime_cl(:)    = 0.0_wp
        veg%fleaf_sunlit_daytime_cl(:)  = 0.0_wp
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
      veg%t_air_daytime_dacc                 = 0.0_wp
      veg%press_srf_daytime_dacc             = 0.0_wp
      veg%co2_mixing_ratio_daytime_dacc      = 0.0_wp
      veg%ga_daytime_dacc                    = 0.0_wp
      veg%beta_air_daytime_dacc              = 0.0_wp
      veg%beta_soil_ps_daytime_dacc          = 0.0_wp
      veg%beta_soil_gs_daytime_dacc          = 0.0_wp
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
  SUBROUTINE calc_moving_averages(lctlib, veg, t_air, press_srf, co2_mixing_ratio, ga)

    USE mo_shared_constants     ! moving average period lengths are defined there as REAL(wp),PARAMETER :: mavg_period_...
    USE mo_veg_memory,          ONLY: t_slm_plant

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),INTENT(in)   :: lctlib             !< land-cover-type library - parameters for one pft
    TYPE(t_slm_plant), INTENT(inout)    :: veg                !< vegetation memory
    REAL(wp), INTENT(in)                :: t_air              !< air temperature ?
    REAL(wp), INTENT(in)                :: press_srf          !< air pressure at surface level ?
    REAL(wp), INTENT(in)                :: co2_mixing_ratio   !< .. ?
    REAL(wp), INTENT(in)                :: ga                 !< aerodynamic resistance
    ! 0.2 Local
    REAl(wp)                            :: hlp


    !> 1.0
    !! moving_avg(current average, new value, length of avg_period, do_calc=.FALSE., avg_period_unit='day')
    !!           returns: new current average

    !> 1.1 averages at the timescale of the labile pool
    veg%gpp_tlabile_mavg                = calc_mavg(veg%gpp_tlabile_mavg, veg%gpp, mavg_period_tlabile)
    veg%maint_respiration_tlabile_mavg  = calc_mavg(veg%maint_respiration_tlabile_mavg, veg%maint_respiration_pot, mavg_period_tlabile)
    veg%growth_req_n_tlabile_mavg       = calc_mavg(veg%growth_req_n_tlabile_mavg, veg%growth_req_n, mavg_period_tlabile, &
                                                    do_calc=(veg%growth_req_n > eps8))
    veg%growth_req_p_tlabile_mavg       = calc_mavg(veg%growth_req_p_tlabile_mavg, veg%growth_req_p, mavg_period_tlabile, &
                                                    do_calc=(veg%growth_req_p > eps8))
    
    !> 1.2 averages at the time-scale of phenology
    veg%t_air_tphen_mavg                 = calc_mavg(veg%t_air_tphen_mavg, t_air, mavg_period_tphen)
    veg%beta_soil_gs_tphen_mavg          = calc_mavg(veg%beta_soil_gs_tphen_mavg, veg%beta_soil_gs, mavg_period_tphen)

    !> 1.3 averages at the nutrient-uptake demand time-scale
    veg%npp_tuptake_mavg                 = calc_mavg(veg%npp_tuptake_mavg, veg%npp, mavg_period_tuptake, &
                                                     do_calc=(veg%pool%fine_root%carbon > eps8))
    veg%demand_uptake_n_tuptake_mavg     = calc_mavg(veg%demand_uptake_n_tuptake_mavg, veg%f_n_demand, mavg_period_tuptake, &
                                                     do_calc=(veg%pool%fine_root%carbon > eps8))
    veg%demand_uptake_p_tuptake_mavg     = calc_mavg(veg%demand_uptake_p_tuptake_mavg, veg%f_p_demand, mavg_period_tuptake, &
                                                     do_calc=(veg%pool%fine_root%carbon > eps8))
    veg%cost_n_fixation_tuptake_mavg     = calc_mavg(veg%cost_n_fixation_tuptake_mavg, veg%cost_n_fixation, mavg_period_tuptake, & 
                                                     do_calc=veg%growing_season) 
    veg%growth_req_n_tuptake_mavg        = calc_mavg(veg%growth_req_n_tuptake_mavg, veg%growth_req_n, mavg_period_tuptake, &
                                                     do_calc=(veg%growth_req_n > eps8))
    veg%growth_req_p_tuptake_mavg        = calc_mavg(veg%growth_req_p_tuptake_mavg, veg%growth_req_p, mavg_period_tuptake, &
                                                     do_calc=(veg%growth_req_p > eps8))

    !> 1.4 averages at the time-scale of within-leaf N allocation fractions 
    !!
    veg%par_sunlit_tfrac_mavg_cl(:)       = calc_mavg(veg%par_sunlit_tfrac_mavg_cl(:), veg%par_sunlit_daytime_cl(:), mavg_period_tfrac)
    veg%par_shaded_tfrac_mavg_cl(:)       = calc_mavg(veg%par_shaded_tfrac_mavg_cl(:), veg%par_shaded_daytime_cl(:), mavg_period_tfrac)
    veg%fleaf_sunlit_tfrac_mavg_cl(:)     = calc_mavg(veg%fleaf_sunlit_tfrac_mavg_cl(:), veg%fleaf_sunlit_daytime_cl(:), mavg_period_tfrac)
    veg%t_air_tfrac_mavg                  = calc_mavg(veg%t_air_tfrac_mavg, veg%t_air_daytime, mavg_period_tfrac)
    veg%press_srf_tfrac_mavg              = calc_mavg(veg%press_srf_tfrac_mavg, veg%press_srf_daytime, mavg_period_tfrac)
    veg%co2_mixing_ratio_tfrac_mavg       = calc_mavg(veg%co2_mixing_ratio_tfrac_mavg, veg%co2_mixing_ratio_daytime, mavg_period_tfrac)
    veg%ga_tfrac_mavg                     = calc_mavg(veg%ga_tfrac_mavg, veg%ga_daytime, mavg_period_tfrac)
    veg%beta_soil_gs_tfrac_mavg           = calc_mavg(veg%beta_soil_gs_tfrac_mavg, veg%beta_soil_gs_daytime, mavg_period_tfrac)
    veg%beta_soil_ps_tfrac_mavg           = calc_mavg(veg%beta_soil_ps_tfrac_mavg, veg%beta_soil_ps_daytime, mavg_period_tfrac)
    veg%beta_air_tfrac_mavg               = calc_mavg(veg%beta_air_tfrac_mavg, veg%beta_air_daytime, mavg_period_tfrac)

    !> 1.5 averages at the time-scale of leaf N allocation fractions 
    !!
    veg%par_sunlit_tcnl_mavg_cl(:)       = calc_mavg(veg%par_sunlit_tcnl_mavg_cl(:), veg%par_sunlit_daytime_cl(:), mavg_period_tcnl)
    veg%par_shaded_tcnl_mavg_cl(:)       = calc_mavg(veg%par_shaded_tcnl_mavg_cl(:), veg%par_shaded_daytime_cl(:), mavg_period_tcnl)
    veg%fleaf_sunlit_tcnl_mavg_cl(:)     = calc_mavg(veg%fleaf_sunlit_tcnl_mavg_cl(:), veg%fleaf_sunlit_daytime_cl(:), mavg_period_tcnl)
    veg%t_air_tcnl_mavg                  = calc_mavg(veg%t_air_tcnl_mavg, veg%t_air_daytime, mavg_period_tcnl)
    veg%press_srf_tcnl_mavg              = calc_mavg(veg%press_srf_tcnl_mavg, veg%press_srf_daytime, mavg_period_tcnl)
    veg%co2_mixing_ratio_tcnl_mavg       = calc_mavg(veg%co2_mixing_ratio_tcnl_mavg, veg%co2_mixing_ratio_daytime, mavg_period_tcnl)
    veg%ga_tcnl_mavg                     = calc_mavg(veg%ga_tcnl_mavg, veg%ga_daytime, mavg_period_tcnl)
    veg%beta_soil_gs_tcnl_mavg           = calc_mavg(veg%beta_soil_gs_tcnl_mavg, veg%beta_soil_gs_daytime, mavg_period_tcnl)
    veg%beta_soil_ps_tcnl_mavg           = calc_mavg(veg%beta_soil_ps_tcnl_mavg, veg%beta_soil_ps_daytime, mavg_period_tcnl)
    veg%beta_air_tcnl_mavg               = calc_mavg(veg%beta_air_tcnl_mavg, veg%beta_air_daytime, mavg_period_tcnl)

    !> 1.6 time-scale of plant allocation calculation
    !!
    veg%npp_talloc_mavg                  = calc_mavg(veg%npp_talloc_mavg, veg%npp, mavg_period_talloc, &
                                                     do_calc=veg%growing_season)
    veg%unit_npp_talloc_mavg             = calc_mavg(veg%unit_npp_talloc_mavg, veg%unit_npp, mavg_period_talloc, &
                                                     do_calc=veg%growing_season)
    veg%transpiration_talloc_mavg        = calc_mavg(veg%transpiration_talloc_mavg, veg%transpiration, mavg_period_talloc, &
                                                     do_calc=veg%growing_season)
    veg%unit_transpiration_talloc_mavg   = calc_mavg(veg%unit_transpiration_talloc_mavg, veg%unit_transpiration, mavg_period_talloc, &
                                                     do_calc=veg%growing_season)
    veg%dphi_talloc_mavg                 = calc_mavg(veg%dphi_talloc_mavg, veg%dphi, mavg_period_talloc, &
                                                     do_calc=veg%growing_season)
    veg%unit_uptake_n_talloc_mavg        = calc_mavg(veg%unit_uptake_n_talloc_mavg, veg%unit_uptake_n, mavg_period_talloc, &
                                                     do_calc=veg%growing_season)
    veg%unit_uptake_p_talloc_mavg        = calc_mavg(veg%unit_uptake_p_talloc_mavg, veg%unit_uptake_p, mavg_period_talloc, &
                                                     do_calc=veg%growing_season)
    veg%labile_carbon_talloc_mavg        = calc_mavg(veg%labile_carbon_talloc_mavg, veg%pool%labile%carbon, & 
                                                     mavg_period_talloc, do_calc=veg%growing_season)
    veg%labile_nitrogen_talloc_mavg      = calc_mavg(veg%labile_nitrogen_talloc_mavg, veg%pool%labile%nitrogen, & 
                                                     mavg_period_talloc, do_calc=veg%growing_season)
    veg%labile_phosphorus_talloc_mavg    = calc_mavg(veg%labile_phosphorus_talloc_mavg, veg%pool%labile%phosphorus, & 
                                                     mavg_period_talloc, do_calc=veg%growing_season)
    IF (veg%growth_req_n > 0._wp)THEN
       hlp = 1._wp/veg%growth_req_n
    ELSE
       hlp=0._wp
    ENDIF
    veg%growth_cn_talloc_mavg  = calc_mavg(veg%growth_cn_talloc_mavg, hlp, mavg_period_talloc, do_calc=veg%growing_season)
    
    IF (veg%growth_req_n > 0._wp .AND. veg%growth_req_p > 0._wp)THEN
       hlp = 1._wp/(veg%growth_req_n * veg%growth_req_p)
    ELSE
       hlp=0._wp
    ENDIF
    veg%growth_cp_talloc_mavg  = calc_mavg(veg%growth_cp_talloc_mavg, hlp, mavg_period_talloc, do_calc=veg%growing_season)

    !> 1.7 averages at the time-scale of root life-time
    !!
    IF (veg%pool%fine_root%carbon > eps8)THEN
       hlp = veg%pool%leaf%carbon/veg%pool%fine_root%carbon
    ELSE
       hlp=0._wp
    ENDIF
    veg%leaf2root_troot_mavg       = calc_mavg(veg%leaf2root_troot_mavg, hlp, lctlib%tau_fine_root, &
                                               do_calc=(veg%pool%leaf%carbon > eps8) .AND. (veg%pool%fine_root%carbon > eps8), & 
                                               avg_period_unit='year')
    veg%unit_npp_troot_mavg        = calc_mavg(veg%unit_npp_troot_mavg, veg%unit_npp, lctlib%tau_fine_root, &
                                               do_calc=veg%growing_season, avg_period_unit='year') 
    veg%unit_uptake_n_troot_mavg   = calc_mavg(veg%unit_uptake_n_troot_mavg, veg%unit_uptake_n, lctlib%tau_fine_root, &
                                               do_calc=veg%growing_season, avg_period_unit='year')
    veg%unit_uptake_p_troot_mavg   = calc_mavg(veg%unit_uptake_p_troot_mavg, veg%unit_uptake_p, lctlib%tau_fine_root, &
                                               do_calc=veg%growing_season, avg_period_unit='year')
    veg%fmaint_rate_troot_mavg     = calc_mavg(veg%fmaint_rate_troot_mavg, veg%fmaint_rate, lctlib%tau_fine_root, &
                                               do_calc=veg%growing_season, avg_period_unit='year') 

    !> 1.8 averages at the vegetation dynamics time-scale 
    veg%an_boc_tvegdyn_mavg   = calc_mavg(veg%an_boc_tvegdyn_mavg, veg%an_boc, mavg_period_tvegdyn, &
                                          do_calc=(veg%pool%leaf%carbon > eps8))
    !> 1.9 averages with a weekly, monthly timescale 
    !!
    veg%t_air_week_mavg       = calc_mavg(veg%t_air_week_mavg, t_air, mavg_period_weekly)
    veg%t_air_month_mavg      = calc_mavg(veg%t_air_month_mavg, t_air, mavg_period_monthly)
    veg%t_jmax_opt_mavg       = calc_mavg(veg%t_jmax_opt_mavg, veg%t_jmax_opt_daytime, mavg_period_weekly)




  END SUBROUTINE calc_moving_averages


END MODULE mo_shared_update_time_averages
