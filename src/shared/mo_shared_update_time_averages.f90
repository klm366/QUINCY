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
  USE mo_jsb_process_class,       ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
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
  SUBROUTINE update_time_averages(tile, time_day, timestep_prop_day)

    USE mo_jsb_class,               ONLY: Get_model
    USE mo_base_memory,             ONLY: t_slm_tile, t_model

    ! Use of process memories
    dsl4jsb_Use_memory(ATM2LAND_)   !  USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,    ONLY: t_pheno_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory


    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)        :: tile      !< one tile
    REAL(wp), INTENT(in)                  :: time_day           !< time in days
    REAL(wp), INTENT(in)                  :: timestep_prop_day  !< length of one time step (proportion of one day)
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_time_averages'
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ATM2LAND_)                   !     TYPE(t_atm2land_memory), POINTER :: atm2land__mem
    dsl4jsb_Def_memory(ASSIMI_)                     !     TYPE(t_assimi_memory),   POINTER :: assimi__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(PHENO_)                      !     TYPE(t_pheno_memory),    POINTER :: pheno__mem
    dsl4jsb_Def_memory(RAD_)                        !     TYPE(t_rad_memory),      POINTER :: rad__mem
    TYPE(t_assimi_memory),   POINTER :: assimi      ! local pointer to assimi__mem for SLM usability
    TYPE(t_atm2land_memory), POINTER :: atm2land    ! local pointer to atm2land__mem for SLM usability
    TYPE(t_veg_memory),      POINTER :: veg         ! local pointer to veg__mem for SLM usability
    TYPE(t_pheno_memory),    POINTER :: pheno       ! local pointer to pheno__mem for SLM usability
    TYPE(t_rad_memory),      POINTER :: rad         ! local pointer to rad__mem for SLM usability
    dsl4jsb_Get_memory(ATM2LAND_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(RAD_)

    assimi    => assimi__mem 
    atm2land  => atm2land__mem
    veg       => veg__mem 
    pheno     => pheno__mem  
    rad       => rad__mem

    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)


    !> 1.0 daytime averages
    !! 
    CALL calc_daytime_averages(assimi, rad, atm2land, veg, time_day, timestep_prop_day)


    !> 2.0 moving averages
    !! 
    CALL calc_moving_averages(lctlib, assimi, rad, atm2land, veg, pheno)



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
  SUBROUTINE calc_daytime_averages(assimi, rad, atm2land, veg, time_day, timestep_prop_day)

    USE mo_assimi_config_class,   ONLY: assimi_config

    ! Use of process memories
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(ATM2LAND_)   !  USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_assimi_memory),  INTENT(inout) :: assimi             !< assimi memory
    TYPE(t_rad_memory),     INTENT(inout) :: rad                !< rad memory
    TYPE(t_atm2land_memory),INTENT(inout) :: atm2land           !< atm2land memory
    TYPE(t_veg_memory),     INTENT(inout) :: veg                !< veg memory
    REAL(wp), INTENT(in)                  :: time_day           !< time in days
    REAL(wp), INTENT(in)                  :: timestep_prop_day  !< length of one time step (proportion of one day)
    ! 0.2 Local
    INTEGER                               :: icanopy
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_daytime_averages'


    !> 1.0 if light is available (i.e. daytime)
    !! 
    IF(atm2land%swpar_srf_down%ptr > eps8) THEN
      ! 1.1 count number of values added to var_daytime
      veg%daytime_counter%ptr = veg%daytime_counter%ptr + 1.0_wp
      ! 1.2 add values and count number of values
      DO icanopy=1,assimi_config%ncanopy
        rad%ppfd_sunlit_daytime_dacc_cl%ptr(icanopy,1)  = rad%ppfd_sunlit_daytime_dacc_cl%ptr(icanopy,1)   + &
                                                           rad%ppfd_sunlit_cl%ptr(icanopy,1)
        rad%ppfd_shaded_daytime_dacc_cl%ptr(icanopy,1)  = rad%ppfd_shaded_daytime_dacc_cl%ptr(icanopy,1)   + & 
                                                           rad%ppfd_shaded_cl%ptr(icanopy,1)
        veg%fleaf_sunlit_daytime_dacc_cl%ptr(icanopy,1) = veg%fleaf_sunlit_daytime_dacc_cl%ptr(icanopy,1) + & 
                                                           veg%fleaf_sunlit_cl%ptr(icanopy,1)
      END DO
      veg%t_air_daytime_dacc%ptr            = veg%t_air_daytime_dacc%ptr            + atm2land%t_air%ptr
      veg%press_srf_daytime_dacc%ptr        = veg%press_srf_daytime_dacc%ptr        + atm2land%press_srf%ptr
      veg%co2_mixing_ratio_daytime_dacc%ptr = veg%co2_mixing_ratio_daytime_dacc%ptr + atm2land%co2_mixing_ratio%ptr
      veg%ga_daytime_dacc%ptr               = veg%ga_daytime_dacc%ptr               + atm2land%ga%ptr
      assimi%beta_air_daytime_dacc%ptr      = assimi%beta_air_daytime_dacc%ptr      + assimi%beta_air%ptr
      assimi%beta_soil_ps_daytime_dacc%ptr  = assimi%beta_soil_ps_daytime_dacc%ptr  + assimi%beta_soil_ps%ptr
      veg%beta_sinklim_ps_daytime_dacc%ptr  = veg%beta_sinklim_ps_daytime_dacc%ptr  + veg%beta_sinklim_ps%ptr
      assimi%beta_soil_gs_daytime_dacc%ptr  = assimi%beta_soil_gs_daytime_dacc%ptr  + assimi%beta_soil_gs%ptr
      veg%t_jmax_opt_daytime_dacc%ptr       = veg%t_jmax_opt_daytime_dacc%ptr       + assimi%t_jmax_opt%ptr
    END IF


    !> 2.0 at the end of the day (midnight), calc daytime averages (of the ending day)
    !! 
    !! calculate the average of the current day and pass this value to the according variable for the "previous day" value
    IF( ABS(((time_day + timestep_prop_day) - REAL(NINT(time_day)))) < eps8 ) THEN
      ! check if at least one value had been assigned at the current day, avoiding division by zero 
      ! calculate the average of the current day and pass it to the according avg_prev_* variable; 
      IF(veg%daytime_counter%ptr > eps8) THEN
        DO icanopy=1,assimi_config%ncanopy
          rad%ppfd_sunlit_daytime_cl%ptr(icanopy,1)   = rad%ppfd_sunlit_daytime_dacc_cl%ptr(icanopy,1)   / veg%daytime_counter%ptr
          rad%ppfd_shaded_daytime_cl%ptr(icanopy,1)   = rad%ppfd_shaded_daytime_dacc_cl%ptr(icanopy,1)   / veg%daytime_counter%ptr
          veg%fleaf_sunlit_daytime_cl%ptr(icanopy,1)  = veg%fleaf_sunlit_daytime_dacc_cl%ptr(icanopy,1) / veg%daytime_counter%ptr
        END DO
        veg%t_air_daytime%ptr            = veg%t_air_daytime_dacc%ptr              / veg%daytime_counter%ptr
        veg%press_srf_daytime%ptr        = veg%press_srf_daytime_dacc%ptr          / veg%daytime_counter%ptr
        veg%co2_mixing_ratio_daytime%ptr = veg%co2_mixing_ratio_daytime_dacc%ptr   / veg%daytime_counter%ptr
        veg%ga_daytime%ptr               = veg%ga_daytime_dacc%ptr                 / veg%daytime_counter%ptr
        assimi%beta_air_daytime%ptr      = assimi%beta_air_daytime_dacc%ptr        / veg%daytime_counter%ptr
        assimi%beta_soil_ps_daytime%ptr  = assimi%beta_soil_ps_daytime_dacc%ptr    / veg%daytime_counter%ptr
        veg%beta_sinklim_ps_daytime%ptr  = veg%beta_sinklim_ps_daytime_dacc%ptr    / veg%daytime_counter%ptr
        assimi%beta_soil_gs_daytime%ptr  = assimi%beta_soil_gs_daytime_dacc%ptr    / veg%daytime_counter%ptr
        veg%t_jmax_opt_daytime%ptr       = veg%t_jmax_opt_daytime_dacc%ptr         / veg%daytime_counter%ptr
      ELSE  ! if veg%daytime_counter is zero
        rad%ppfd_sunlit_daytime_cl%ptr   = 0.0_wp
        rad%ppfd_shaded_daytime_cl%ptr   = 0.0_wp
        veg%fleaf_sunlit_daytime_cl%ptr  = 0.0_wp
      END IF

!      write(999,*) time_day, veg%daytime_counter%ptr, &   ! might be used for testing
!                  rad%ppfd_sunlit_daytime_cl%ptr(1,1), &
!                  rad%ppfd_shaded_daytime_cl%ptr(1,1), &
!                  veg%fleaf_sunlit_daytime_cl%ptr(1,1), &
!                  veg%t_jmax_opt_daytime%ptr, ' |-| '

      ! set all values of the current day to zero
      rad%ppfd_sunlit_daytime_dacc_cl%ptr        = 0.0_wp
      rad%ppfd_shaded_daytime_dacc_cl%ptr        = 0.0_wp
      veg%fleaf_sunlit_daytime_dacc_cl%ptr       = 0.0_wp
      veg%t_air_daytime_dacc%ptr                 = 0.0_wp
      veg%press_srf_daytime_dacc%ptr             = 0.0_wp
      veg%co2_mixing_ratio_daytime_dacc%ptr      = 0.0_wp
      veg%ga_daytime_dacc%ptr                    = 0.0_wp
      assimi%beta_air_daytime_dacc%ptr           = 0.0_wp
      assimi%beta_soil_ps_daytime_dacc%ptr       = 0.0_wp
      veg%beta_sinklim_ps_daytime_dacc%ptr       = 0.0_wp
      assimi%beta_soil_gs_daytime_dacc%ptr       = 0.0_wp
      veg%t_jmax_opt_daytime_dacc%ptr            = 0.0_wp
      veg%daytime_counter%ptr                    = 0.0_wp
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

    USE mo_base_constants,      ONLY: dtime, one_day, one_year

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
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_mavg'


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
      new_value_weight = 1.0_wp / (avg_period_length * avg_period2seconds / dtime )  ! dtime: length of one time step in seconds
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
  SUBROUTINE calc_moving_averages(lctlib, assimi, rad, atm2land, veg, pheno)

    USE mo_shared_constants     ! moving average period lengths are defined there as REAL(wp),PARAMETER :: mavg_period_...

    ! Use of process memories
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(ATM2LAND_)   !  USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,    ONLY: t_pheno_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory


    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element), INTENT(in)      :: lctlib             !< land-cover-type library - parameters for one pft
    TYPE(t_assimi_memory),  INTENT(inout)   :: assimi             !< assimi memory
    TYPE(t_rad_memory),     INTENT(inout)   :: rad                !< rad memory
    TYPE(t_atm2land_memory),INTENT(inout)   :: atm2land           !< atm2land memory
    TYPE(t_veg_memory),     INTENT(inout)   :: veg                !< veg memory
    TYPE(t_pheno_memory),   INTENT(inout)   :: pheno              !< pheno memory
    ! 0.2 Local
    REAl(wp)                                :: hlp
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_moving_averages'


    !> 1.0
    !! moving_avg(current average, new value, length of avg_period, do_calc=.FALSE., avg_period_unit='day')
    !!           returns: new current average

    !> 1.1 averages at the timescale of the labile pool
    veg%gpp_tlabile_mavg%ptr                = calc_mavg(veg%gpp_tlabile_mavg%ptr, assimi%gross_assimilation%ptr, &
                                                    mavg_period_tlabile)
    veg%maint_respiration_tlabile_mavg%ptr  = calc_mavg(veg%maint_respiration_tlabile_mavg%ptr, veg%maint_respiration_pot%ptr, &
                                                    mavg_period_tlabile)
    veg%growth_req_n_tlabile_mavg%ptr       = calc_mavg(veg%growth_req_n_tlabile_mavg%ptr, veg%growth_req_n%ptr, &
                                                    mavg_period_tlabile, do_calc=(veg%growth_req_n%ptr > eps8))
    veg%growth_req_p_tlabile_mavg%ptr       = calc_mavg(veg%growth_req_p_tlabile_mavg%ptr, veg%growth_req_p%ptr, &
                                                    mavg_period_tlabile, do_calc=(veg%growth_req_p%ptr > eps8))
    
    !> 1.2 averages at the time-scale of phenology
    veg%t_air_tphen_mavg%ptr           = calc_mavg(veg%t_air_tphen_mavg%ptr, atm2land%t_air%ptr, mavg_period_tphen)
    assimi%beta_soil_gs_tphen_mavg%ptr = calc_mavg(assimi%beta_soil_gs_tphen_mavg%ptr, assimi%beta_soil_gs%ptr, mavg_period_tphen)

    !> 1.3 averages at the nutrient-uptake demand time-scale
    veg%npp_tuptake_mavg%ptr                 = calc_mavg(veg%npp_tuptake_mavg%ptr, veg%npp%ptr, mavg_period_tuptake, &
                                                     do_calc=(veg%pool%fine_root%carbon%ptr > eps8))
    veg%demand_uptake_n_tuptake_mavg%ptr     = calc_mavg(veg%demand_uptake_n_tuptake_mavg%ptr, veg%f_n_demand%ptr, &
                                                     mavg_period_tuptake, do_calc=(veg%pool%fine_root%carbon%ptr > eps8))
    veg%demand_uptake_p_tuptake_mavg%ptr     = calc_mavg(veg%demand_uptake_p_tuptake_mavg%ptr, veg%f_p_demand%ptr, &
                                                     mavg_period_tuptake, do_calc=(veg%pool%fine_root%carbon%ptr > eps8))
    veg%cost_n_fixation_tuptake_mavg%ptr     = calc_mavg(veg%cost_n_fixation_tuptake_mavg%ptr, veg%cost_n_fixation%ptr, &
                                                     mavg_period_tuptake, do_calc=pheno%growing_season%ptr) 
    veg%growth_req_n_tuptake_mavg%ptr        = calc_mavg(veg%growth_req_n_tuptake_mavg%ptr, veg%growth_req_n%ptr, &
                                                     mavg_period_tuptake, do_calc=(veg%growth_req_n%ptr > eps8))
    veg%growth_req_p_tuptake_mavg%ptr        = calc_mavg(veg%growth_req_p_tuptake_mavg%ptr, veg%growth_req_p%ptr, &
                                                     mavg_period_tuptake, do_calc=(veg%growth_req_p%ptr > eps8))

    !> 1.4 averages at the time-scale of within-leaf N allocation fractions 
    !!
    rad%ppfd_sunlit_tfrac_mavg_cl%ptr     = calc_mavg(rad%ppfd_sunlit_tfrac_mavg_cl%ptr, rad%ppfd_sunlit_daytime_cl%ptr, &
                                                  mavg_period_tfrac)
    rad%ppfd_shaded_tfrac_mavg_cl%ptr     = calc_mavg(rad%ppfd_shaded_tfrac_mavg_cl%ptr, rad%ppfd_shaded_daytime_cl%ptr, &
                                                  mavg_period_tfrac)
    veg%fleaf_sunlit_tfrac_mavg_cl%ptr    = calc_mavg(veg%fleaf_sunlit_tfrac_mavg_cl%ptr, veg%fleaf_sunlit_daytime_cl%ptr, &
                                                  mavg_period_tfrac)
    veg%t_air_tfrac_mavg%ptr              = calc_mavg(veg%t_air_tfrac_mavg%ptr, veg%t_air_daytime%ptr, &
                                                  mavg_period_tfrac)
    veg%press_srf_tfrac_mavg%ptr          = calc_mavg(veg%press_srf_tfrac_mavg%ptr, veg%press_srf_daytime%ptr, &
                                                  mavg_period_tfrac)
    veg%co2_mixing_ratio_tfrac_mavg%ptr   = calc_mavg(veg%co2_mixing_ratio_tfrac_mavg%ptr, veg%co2_mixing_ratio_daytime%ptr, &
                                                  mavg_period_tfrac)
    veg%ga_tfrac_mavg%ptr                 = calc_mavg(veg%ga_tfrac_mavg%ptr, veg%ga_daytime%ptr, &
                                                  mavg_period_tfrac)
    assimi%beta_soil_ps_tfrac_mavg%ptr    = calc_mavg(assimi%beta_soil_ps_tfrac_mavg%ptr, assimi%beta_soil_ps_daytime%ptr, &
                                                  mavg_period_tfrac)
    veg%beta_sinklim_ps_tfrac_mavg%ptr    = calc_mavg(veg%beta_sinklim_ps_tfrac_mavg%ptr, veg%beta_sinklim_ps_daytime%ptr, &
                                                  mavg_period_tfrac)
    assimi%beta_soil_gs_tfrac_mavg%ptr    = calc_mavg(assimi%beta_soil_gs_tfrac_mavg%ptr, assimi%beta_soil_gs_daytime%ptr, &
                                                  mavg_period_tfrac)
    assimi%beta_air_tfrac_mavg%ptr        = calc_mavg(assimi%beta_air_tfrac_mavg%ptr, assimi%beta_air_daytime%ptr, &
                                                  mavg_period_tfrac)

    !> 1.5 averages at the time-scale of leaf N allocation fractions 
    !!
    rad%ppfd_sunlit_tcnl_mavg_cl%ptr      = calc_mavg(rad%ppfd_sunlit_tcnl_mavg_cl%ptr, rad%ppfd_sunlit_daytime_cl%ptr, &
                                                  mavg_period_tcnl)
    rad%ppfd_shaded_tcnl_mavg_cl%ptr      = calc_mavg(rad%ppfd_shaded_tcnl_mavg_cl%ptr, rad%ppfd_shaded_daytime_cl%ptr, &
                                                  mavg_period_tcnl)
    veg%fleaf_sunlit_tcnl_mavg_cl%ptr     = calc_mavg(veg%fleaf_sunlit_tcnl_mavg_cl%ptr, veg%fleaf_sunlit_daytime_cl%ptr, &
                                                  mavg_period_tcnl)
    veg%t_air_tcnl_mavg%ptr               = calc_mavg(veg%t_air_tcnl_mavg%ptr, veg%t_air_daytime%ptr, &
                                                  mavg_period_tcnl)
    veg%press_srf_tcnl_mavg%ptr           = calc_mavg(veg%press_srf_tcnl_mavg%ptr, veg%press_srf_daytime%ptr, &
                                                  mavg_period_tcnl)
    veg%co2_mixing_ratio_tcnl_mavg%ptr    = calc_mavg(veg%co2_mixing_ratio_tcnl_mavg%ptr, veg%co2_mixing_ratio_daytime%ptr, &
                                                  mavg_period_tcnl)
    veg%ga_tcnl_mavg%ptr                  = calc_mavg(veg%ga_tcnl_mavg%ptr, veg%ga_daytime%ptr, &
                                                  mavg_period_tcnl)
    assimi%beta_soil_ps_tcnl_mavg%ptr     = calc_mavg(assimi%beta_soil_ps_tcnl_mavg%ptr, assimi%beta_soil_ps_daytime%ptr, &
                                                  mavg_period_tcnl)
    veg%beta_sinklim_ps_tcnl_mavg%ptr     = calc_mavg(veg%beta_sinklim_ps_tcnl_mavg%ptr, veg%beta_sinklim_ps_daytime%ptr, &
                                                  mavg_period_tcnl)
    assimi%beta_soil_gs_tcnl_mavg%ptr     = calc_mavg(assimi%beta_soil_gs_tcnl_mavg%ptr, assimi%beta_soil_gs_daytime%ptr, &
                                                  mavg_period_tcnl)
    assimi%beta_air_tcnl_mavg%ptr         = calc_mavg(assimi%beta_air_tcnl_mavg%ptr, assimi%beta_air_daytime%ptr, &
                                                  mavg_period_tcnl)

    !> 1.6 time-scale of plant allocation calculation (veg process)
    !!
    veg%npp_talloc_mavg%ptr                  = calc_mavg(veg%npp_talloc_mavg%ptr, veg%npp%ptr, &
                                                      mavg_period_talloc, do_calc=pheno%growing_season%ptr)
    veg%unit_npp_talloc_mavg%ptr             = calc_mavg(veg%unit_npp_talloc_mavg%ptr, veg%unit_npp%ptr, &
                                                      mavg_period_talloc, do_calc=pheno%growing_season%ptr)
    veg%transpiration_talloc_mavg%ptr        = calc_mavg(veg%transpiration_talloc_mavg%ptr, veg%transpiration%ptr, &
                                                      mavg_period_talloc, do_calc=pheno%growing_season%ptr)
    veg%unit_transpiration_talloc_mavg%ptr   = calc_mavg(veg%unit_transpiration_talloc_mavg%ptr, veg%unit_transpiration%ptr, &
                                                      mavg_period_talloc, do_calc=pheno%growing_season%ptr)
    veg%dphi_talloc_mavg%ptr                 = calc_mavg(veg%dphi_talloc_mavg%ptr, veg%dphi%ptr, &
                                                      mavg_period_talloc, do_calc=pheno%growing_season%ptr)
    veg%unit_uptake_n_act_talloc_mavg%ptr    = calc_mavg(veg%unit_uptake_n_act_talloc_mavg%ptr, veg%unit_uptake_n_act%ptr, &
                                                      mavg_period_talloc, do_calc=pheno%growing_season%ptr)
    veg%unit_uptake_p_act_talloc_mavg%ptr    = calc_mavg(veg%unit_uptake_p_act_talloc_mavg%ptr, veg%unit_uptake_p_act%ptr, &
                                                      mavg_period_talloc, do_calc=pheno%growing_season%ptr)
    veg%labile_carbon_talloc_mavg%ptr        = calc_mavg(veg%labile_carbon_talloc_mavg%ptr, veg%pool%labile%carbon%ptr, & 
                                                      mavg_period_talloc, do_calc=pheno%growing_season%ptr)
    veg%labile_nitrogen_talloc_mavg%ptr      = calc_mavg(veg%labile_nitrogen_talloc_mavg%ptr, veg%pool%labile%nitrogen%ptr, & 
                                                      mavg_period_talloc, do_calc=pheno%growing_season%ptr)
    veg%labile_phosphorus_talloc_mavg%ptr    = calc_mavg(veg%labile_phosphorus_talloc_mavg%ptr, veg%pool%labile%phosphorus%ptr, & 
                                                      mavg_period_talloc, do_calc=pheno%growing_season%ptr)
    veg%reserve_carbon_talloc_mavg%ptr   = calc_mavg(veg%reserve_carbon_talloc_mavg%ptr, veg%pool%reserve%carbon%ptr, & 
                                                      mavg_period_talloc, do_calc=pheno%growing_season%ptr)
    veg%reserve_nitrogen_talloc_mavg%ptr      = calc_mavg(veg%reserve_nitrogen_talloc_mavg%ptr, veg%pool%reserve%nitrogen%ptr, & 
                                                      mavg_period_talloc, do_calc=pheno%growing_season%ptr)
    veg%reserve_phosphorus_talloc_mavg%ptr    = calc_mavg(veg%reserve_phosphorus_talloc_mavg%ptr, veg%pool%reserve%phosphorus%ptr, & 
                                                      mavg_period_talloc, do_calc=pheno%growing_season%ptr)
    IF (veg%growth_req_n%ptr > 0._wp)THEN
       hlp = 1._wp/veg%growth_req_n%ptr
    ELSE
       hlp=0._wp
    ENDIF
    veg%growth_cn_talloc_mavg%ptr  = calc_mavg(veg%growth_cn_talloc_mavg%ptr, hlp, &
                                            mavg_period_tgrowth, do_calc=pheno%growing_season%ptr)
    
    IF (veg%growth_req_n%ptr > 0._wp .AND. veg%growth_req_p%ptr > 0._wp)THEN
       hlp = 1._wp/(veg%growth_req_n%ptr * veg%growth_req_p%ptr)
    ELSE
       hlp=0._wp
    ENDIF
    
    veg%growth_cp_talloc_mavg%ptr  = calc_mavg(veg%growth_cp_talloc_mavg%ptr, hlp, &
                                            mavg_period_tgrowth, do_calc=pheno%growing_season%ptr)

    !> 1.7 averages at the time-scale of root life-time
    !!
    IF (veg%pool%fine_root%carbon%ptr > eps8)THEN
       hlp = veg%pool%leaf%carbon%ptr / veg%pool%fine_root%carbon%ptr
    ELSE
       hlp=0._wp
    ENDIF
    veg%leaf2root_troot_mavg%ptr      = calc_mavg(veg%leaf2root_troot_mavg%ptr, hlp, lctlib%tau_fine_root, &
                                               do_calc=(veg%pool%leaf%carbon%ptr > eps8) .AND. &
                                                       (veg%pool%fine_root%carbon%ptr > eps8), & 
                                               avg_period_unit='year')
    veg%unit_npp_troot_mavg%ptr       = calc_mavg(veg%unit_npp_troot_mavg%ptr, veg%unit_npp%ptr, lctlib%tau_fine_root, &
                                               do_calc=pheno%growing_season%ptr, avg_period_unit='year') 
    veg%unit_uptake_n_pot_troot_mavg%ptr  = calc_mavg(veg%unit_uptake_n_pot_troot_mavg%ptr, veg%unit_uptake_n_pot%ptr, lctlib%tau_fine_root, &
                                               do_calc=pheno%growing_season%ptr, avg_period_unit='year')
    veg%unit_uptake_p_pot_troot_mavg%ptr  = calc_mavg(veg%unit_uptake_p_pot_troot_mavg%ptr, veg%unit_uptake_p_pot%ptr, lctlib%tau_fine_root, &
                                               do_calc=pheno%growing_season%ptr, avg_period_unit='year')
    veg%fmaint_rate_troot_mavg%ptr    = calc_mavg(veg%fmaint_rate_troot_mavg%ptr, veg%fmaint_rate%ptr, lctlib%tau_fine_root, &
                                               do_calc=pheno%growing_season%ptr, avg_period_unit='year') 

    !> 1.8 averages at the vegetation dynamics time-scale 
    veg%an_boc_tvegdyn_mavg%ptr   = calc_mavg(veg%an_boc_tvegdyn_mavg%ptr, assimi%net_assimilation_boc%ptr, mavg_period_tvegdyn, &
                                          do_calc=(veg%pool%leaf%carbon%ptr > eps8))
    !> 1.9 averages at the respiration acclimation time-scale
    veg%t_air_tacclim_mavg%ptr    = calc_mavg(veg%t_air_tacclim_mavg%ptr, atm2land%t_air%ptr, mavg_period_tacclim)
    !> 1.10 averages with a weekly, monthly timescale 
    !!
    veg%t_air_week_mavg%ptr       = calc_mavg(veg%t_air_week_mavg%ptr, atm2land%t_air%ptr, mavg_period_weekly)
    veg%t_air_month_mavg%ptr      = calc_mavg(veg%t_air_month_mavg%ptr, atm2land%t_air%ptr, mavg_period_monthly)
    veg%t_jmax_opt_mavg%ptr       = calc_mavg(veg%t_jmax_opt_mavg%ptr, veg%t_jmax_opt_daytime%ptr, mavg_period_weekly)




  END SUBROUTINE calc_moving_averages


END MODULE mo_shared_update_time_averages
