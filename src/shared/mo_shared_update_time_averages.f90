!>
!! @brief 
!!    generic routines providing various ways for the averaging of variables over time for all model processes
!!
!! @author
!!    Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!!    First version by Jan Engel (2016-09-01) - the basic routines were previously implemented in SLM mo_veg_memory
!!    2nd version by Jan Engel (20170120) - model structure moved from Main-Task into the Sub-Tasks
!!
MODULE mo_shared_update_time_averages

#include "dsl_jsbach.inc"

  USE mo_kind,                    ONLY: wp, eps8
  USE mo_exception,               ONLY: message, finish
  USE mo_jsb_process_class,       ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
  USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
  USE mo_base_memory,             ONLY: t_slm_tile, t_model, t_jsb_task_options


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_time_averages


  CHARACTER(len=*), PARAMETER :: modname = 'mo_shared_update_time_averages'

CONTAINS 


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update time averages across processes at the end of a time step
  !! 
  !! no aggregation of any avg variable to Top Tile (not using Integrate Aggregate Routines)
  !! 
  !! 1.0 daytime averages \n
  !! 2.0 moving averages \n
  !! 3.0 ...
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_time_averages(tile, options, time_day, timestep_prop_day)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout) :: tile             !< one tile
    TYPE(t_jsb_task_options), INTENT(in)  :: options            !< model options
    REAL(wp), INTENT(in)                  :: time_day           !< time in days
    REAL(wp), INTENT(in)                  :: timestep_prop_day  !< length of one time step (proportion of one day)
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_time_averages'



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 daytime averages
    !! 
    CALL calc_daytime_averages(tile, options, time_day, timestep_prop_day)


    !> 2.0 moving averages
    !! 
    CALL calc_moving_averages(tile, options)


    !> 3.0 .... averages < open for implementation, just add another routine to this module and call it here
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
  SUBROUTINE calc_daytime_averages(tile, options, time_day, timestep_prop_day)

    USE mo_jsb_class,               ONLY: Get_model
    USE mo_assimi_config_class,     ONLY: assimi_config

    ! Use of process memories
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(A2L_)   !  USE mo_a2l_memory_class, ONLY: t_a2l_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout) :: tile             !< one tile
    TYPE(t_jsb_task_options), INTENT(in)  :: options            !< model options
    REAL(wp), INTENT(in)                  :: time_day           !< time in days
    REAL(wp), INTENT(in)                  :: timestep_prop_day  !< length of one time step (proportion of one day)
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: icanopy    !< loop over ncanopy
    INTEGER                               :: ichunk     !< loop over chunks ics-ice
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_daytime_averages'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(A2L_)                   !     TYPE(t_a2l_memory), POINTER :: a2l__mem
    dsl4jsb_Def_memory(ASSIMI_)                     !     TYPE(t_assimi_memory),   POINTER :: assimi__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(RAD_)                        !     TYPE(t_rad_memory),      POINTER :: rad__mem
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk      :: t_air                        ! A2L_
    dsl4jsb_Real2D_onChunk      :: swpar_srf_down
    dsl4jsb_Real2D_onChunk      :: press_srf
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio
    dsl4jsb_Real2D_onChunk      :: ga
    dsl4jsb_Real2D_onChunk      :: beta_air                     ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs
    dsl4jsb_Real2D_onChunk      :: beta_air_daytime
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_daytime
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_daytime
    dsl4jsb_Real2D_onChunk      :: beta_air_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_cl               ! RAD_
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_daytime_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_daytime_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_daytime_dacc_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_daytime_dacc_cl
    dsl4jsb_Real2D_onChunk      :: daytime_counter              ! VEG_
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps
    dsl4jsb_Real2D_onChunk      :: t_air_daytime
    dsl4jsb_Real2D_onChunk      :: press_srf_daytime
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_daytime
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_daytime
    dsl4jsb_Real2D_onChunk      :: ga_daytime
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_daytime
    dsl4jsb_Real2D_onChunk      :: press_srf_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: ga_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: t_air_daytime_dacc
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_daytime_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_daytime_dacc_cl
!! *slmtemp
    TYPE(t_veg_memory),      POINTER :: veg         ! local pointer to veg__mem for SLM usability
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(VEG_)) RETURN   ! ask for all processes considered here ?
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(PHENO_)
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(RAD_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)
    dsl4jsb_Get_var2D_onChunk(A2L_, co2_mixing_ratio)
    dsl4jsb_Get_var2D_onChunk(A2L_, ga)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_daytime)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_daytime)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_daytime)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, t_jmax_opt)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_daytime_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_daytime_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_daytime_dacc_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_daytime_dacc_cl)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, daytime_counter)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_daytime_dacc)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_daytime_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_daytime_dacc_cl)


!! *slmtemp
    veg       => veg__mem 



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 if light is available (i.e. daytime)
    !! 
    DO ichunk=1,nc
      IF(swpar_srf_down(ichunk) > eps8) THEN
        ! 1.1 count number of values added to var_daytime
        daytime_counter(ichunk) = daytime_counter(ichunk) + 1.0_wp
        ! 1.2 add values and count number of values
        DO icanopy=1,assimi_config%ncanopy
          ppfd_sunlit_daytime_dacc_cl(ichunk,icanopy)  = ppfd_sunlit_daytime_dacc_cl(ichunk,icanopy)   + &
                                                             ppfd_sunlit_cl(ichunk,icanopy)
          ppfd_shaded_daytime_dacc_cl(ichunk,icanopy)  = ppfd_shaded_daytime_dacc_cl(ichunk,icanopy)   + & 
                                                             ppfd_shaded_cl(ichunk,icanopy)
          fleaf_sunlit_daytime_dacc_cl(ichunk,icanopy) = fleaf_sunlit_daytime_dacc_cl(ichunk,icanopy) + & 
                                                             fleaf_sunlit_cl(ichunk,icanopy)
        END DO
        t_air_daytime_dacc(ichunk)            = t_air_daytime_dacc(ichunk)            + t_air(ichunk)
        press_srf_daytime_dacc(ichunk)        = press_srf_daytime_dacc(ichunk)        + press_srf(ichunk)
        co2_mixing_ratio_daytime_dacc(ichunk) = co2_mixing_ratio_daytime_dacc(ichunk) + co2_mixing_ratio(ichunk)
        ga_daytime_dacc(ichunk)               = ga_daytime_dacc(ichunk)               + ga(ichunk)
        beta_air_daytime_dacc(ichunk)         = beta_air_daytime_dacc(ichunk)         + beta_air(ichunk)
        beta_soil_ps_daytime_dacc(ichunk)     = beta_soil_ps_daytime_dacc(ichunk)     + beta_soil_ps(ichunk)
        beta_sinklim_ps_daytime_dacc(ichunk)  = beta_sinklim_ps_daytime_dacc(ichunk)  + beta_sinklim_ps(ichunk)
        beta_soil_gs_daytime_dacc(ichunk)     = beta_soil_gs_daytime_dacc(ichunk)     + beta_soil_gs(ichunk)
        t_jmax_opt_daytime_dacc(ichunk)       = t_jmax_opt_daytime_dacc(ichunk)       + t_jmax_opt(ichunk)
      END IF
    END DO


    !> 2.0 at the end of the day (midnight), calc daytime averages (of the ending day)
    !! 
    !! calculate the average of the current day and pass this value to the according variable for the "previous day" value
    IF( ABS(((time_day + timestep_prop_day) - REAL(NINT(time_day)))) < eps8 ) THEN
      ! check if at least one value had been assigned at the current day, avoiding division by zero 
      ! calculate the average of the current day and pass it to the according avg_prev_* variable; 
      DO ichunk=1,nc
        IF(daytime_counter(ichunk) > eps8) THEN
          DO icanopy=1,assimi_config%ncanopy
            ppfd_sunlit_daytime_cl(ichunk,icanopy)   = ppfd_sunlit_daytime_dacc_cl(ichunk,icanopy)   / daytime_counter(ichunk)
            ppfd_shaded_daytime_cl(ichunk,icanopy)   = ppfd_shaded_daytime_dacc_cl(ichunk,icanopy)   / daytime_counter(ichunk)
            fleaf_sunlit_daytime_cl(ichunk,icanopy)  = fleaf_sunlit_daytime_dacc_cl(ichunk,icanopy)  / daytime_counter(ichunk)
          END DO
          t_air_daytime(ichunk)            = t_air_daytime_dacc(ichunk)              / daytime_counter(ichunk)
          press_srf_daytime(ichunk)        = press_srf_daytime_dacc(ichunk)          / daytime_counter(ichunk)
          co2_mixing_ratio_daytime(ichunk) = co2_mixing_ratio_daytime_dacc(ichunk)   / daytime_counter(ichunk)
          ga_daytime(ichunk)               = ga_daytime_dacc(ichunk)                 / daytime_counter(ichunk)
          beta_air_daytime(ichunk)         = beta_air_daytime_dacc(ichunk)           / daytime_counter(ichunk)
          beta_soil_ps_daytime(ichunk)     = beta_soil_ps_daytime_dacc(ichunk)       / daytime_counter(ichunk)
          beta_sinklim_ps_daytime(ichunk)  = beta_sinklim_ps_daytime_dacc(ichunk)    / daytime_counter(ichunk)
          beta_soil_gs_daytime(ichunk)     = beta_soil_gs_daytime_dacc(ichunk)       / daytime_counter(ichunk)
          t_jmax_opt_daytime(ichunk)       = t_jmax_opt_daytime_dacc(ichunk)         / daytime_counter(ichunk)
        ELSE  ! if daytime_counter is zero
          ppfd_sunlit_daytime_cl(ichunk,:)  = 0.0_wp
          ppfd_shaded_daytime_cl(ichunk,:)  = 0.0_wp
          fleaf_sunlit_daytime_cl(ichunk,:) = 0.0_wp
        END IF
      END DO  ! loop ics:ice

!      write(999,*) time_day, daytime_counter(1), &   ! might be used for testing
!                  ppfd_sunlit_daytime_cl(1,1), &
!                  ppfd_shaded_daytime_cl(1,1), &
!                  fleaf_sunlit_daytime_cl(1,1), &
!                  t_jmax_opt_daytime(1), ' |-| '

      ! set all values of the current day to zero
      ppfd_sunlit_daytime_dacc_cl(:,:)      = 0.0_wp
      ppfd_shaded_daytime_dacc_cl(:,:)      = 0.0_wp
      fleaf_sunlit_daytime_dacc_cl(:,:)     = 0.0_wp
      t_air_daytime_dacc(:)                 = 0.0_wp
      press_srf_daytime_dacc(:)             = 0.0_wp
      co2_mixing_ratio_daytime_dacc(:)      = 0.0_wp
      ga_daytime_dacc(:)                    = 0.0_wp
      beta_air_daytime_dacc(:)              = 0.0_wp
      beta_soil_ps_daytime_dacc(:)          = 0.0_wp
      beta_sinklim_ps_daytime_dacc(:)       = 0.0_wp
      beta_soil_gs_daytime_dacc(:)          = 0.0_wp
      t_jmax_opt_daytime_dacc(:)            = 0.0_wp
      daytime_counter(:)                    = 0.0_wp
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
  SUBROUTINE calc_moving_averages(tile, options)

    USE mo_jsb_class,               ONLY: Get_model
    USE mo_shared_constants           ! moving average period lengths are defined there as REAL(wp),PARAMETER :: mavg_period_...

    ! Use of process memories
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(A2L_)   !  USE mo_a2l_memory_class, ONLY: t_a2l_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,    ONLY: t_pheno_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout) :: tile             !< one tile
    TYPE(t_jsb_task_options), INTENT(in)  :: options            !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    REAl(wp),DIMENSION(options%nc)        :: hlp
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_moving_averages'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(A2L_)                   !     TYPE(t_a2l_memory), POINTER :: a2l__mem
    dsl4jsb_Def_memory(ASSIMI_)                     !     TYPE(t_assimi_memory),   POINTER :: assimi__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(PHENO_)                      !     TYPE(t_pheno_memory),    POINTER :: pheno__mem
    dsl4jsb_Def_memory(RAD_)                        !     TYPE(t_rad_memory),      POINTER :: rad__mem
    ! Declare pointers to variables in memory
    dsl4jsb_Logical2D_onChunk   :: growing_season             ! PHENO_
    dsl4jsb_Real2D_onChunk      :: t_air                      ! A2L_
    dsl4jsb_Real2D_onChunk      :: gross_assimilation         ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: net_assimilation_boc
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_daytime
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_daytime
    dsl4jsb_Real2D_onChunk      :: beta_air_daytime
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_air_tcnl_mavg
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_daytime_cl     ! RAD_
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_daytime_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tcnl_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tcnl_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tfrac_mavg_cl
    dsl4jsb_Real2D_onChunk      :: n_fixation                 ! VEG_
    dsl4jsb_Real2D_onChunk      :: press_srf_daytime
    dsl4jsb_Real2D_onChunk      :: t_air_daytime
    dsl4jsb_Real2D_onChunk      :: cost_n_fixation
    dsl4jsb_Real2D_onChunk      :: f_p_demand
    dsl4jsb_Real2D_onChunk      :: f_n_demand
    dsl4jsb_Real2D_onChunk      :: npp
    dsl4jsb_Real2D_onChunk      :: maint_respiration_pot
    dsl4jsb_Real2D_onChunk      :: growth_req_n
    dsl4jsb_Real2D_onChunk      :: growth_req_p
    dsl4jsb_Real2D_onChunk      :: unit_npp
    dsl4jsb_Real2D_onChunk      :: transpiration
    dsl4jsb_Real2D_onChunk      :: unit_transpiration
    dsl4jsb_Real2D_onChunk      :: dphi
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_act
    dsl4jsb_Real2D_onChunk      :: unit_uptake_p_act
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_pot
    dsl4jsb_Real2D_onChunk      :: unit_uptake_p_pot
    dsl4jsb_Real2D_onChunk      :: fmaint_rate
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_daytime
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_daytime
    dsl4jsb_Real2D_onChunk      :: ga_daytime
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_daytime
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_month_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_week_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tacclim_mavg
    dsl4jsb_Real2D_onChunk      :: an_boc_tvegdyn_mavg
    dsl4jsb_Real2D_onChunk      :: fmaint_rate_troot_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_p_pot_troot_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_pot_troot_mavg
    dsl4jsb_Real2D_onChunk      :: unit_npp_troot_mavg
    dsl4jsb_Real2D_onChunk      :: leaf2root_troot_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cn_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_npp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: n_fixation_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: npp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: ga_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: press_srf_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: gpp_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: maint_respiration_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_p_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: npp_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: demand_uptake_n_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: demand_uptake_p_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: cost_n_fixation_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_p_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: press_srf_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: ga_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: transpiration_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_transpiration_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: dphi_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_act_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_p_act_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_carbon_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_carbon_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_phosphorus_talloc_mavg
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tcnl_mavg_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_daytime_cl
!! *slmtemp
    TYPE(t_veg_memory),      POINTER :: veg         ! local pointer to veg__mem for SLM usability
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(VEG_)) RETURN   ! ask for all processes considered here ?
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(PHENO_)
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(RAD_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(PHENO_, growing_season)           ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, gross_assimilation)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, net_assimilation_boc)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_daytime)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_daytime)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_daytime)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tphen_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_tcnl_mavg)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_daytime_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_daytime_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tcnl_mavg_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tcnl_mavg_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tfrac_mavg_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tfrac_mavg_cl)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, cost_n_fixation)
    dsl4jsb_Get_var2D_onChunk(VEG_, f_p_demand)
    dsl4jsb_Get_var2D_onChunk(VEG_, f_n_demand)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_pot)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_npp)
    dsl4jsb_Get_var2D_onChunk(VEG_, transpiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_transpiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, dphi)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_act)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_act)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_pot)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_pot)
    dsl4jsb_Get_var2D_onChunk(VEG_, fmaint_rate)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_daytime)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_month_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_week_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tacclim_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, an_boc_tvegdyn_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, fmaint_rate_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_pot_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_pot_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_npp_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, leaf2root_troot_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cn_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_npp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, gpp_tlabile_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_tlabile_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tlabile_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tlabile_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tphen_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, demand_uptake_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, demand_uptake_p_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, cost_n_fixation_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, dphi_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_act_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_act_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_carbon_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_nitrogen_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_phosphorus_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_carbon_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_nitrogen_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_phosphorus_talloc_mavg)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tcnl_mavg_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_daytime_cl)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

!! *slmtemp
    veg       => veg__mem 



    !> 1.0
    !! calc_mavg(current average, new value, length of avg_period, do_calc=.FALSE., avg_period_unit='day')
    !!           returns: new current average

    !> 1.1 averages at the timescale of the labile pool
    gpp_tlabile_mavg(:)                = calc_mavg(gpp_tlabile_mavg(:), gross_assimilation(:), &
                                                    mavg_period_tlabile)
    maint_respiration_tlabile_mavg(:)  = calc_mavg(maint_respiration_tlabile_mavg(:), maint_respiration_pot(:), &
                                                    mavg_period_tlabile)
    growth_req_n_tlabile_mavg(:)       = calc_mavg(growth_req_n_tlabile_mavg(:), growth_req_n(:), &
                                                    mavg_period_tlabile, do_calc=(growth_req_n(:) > eps8))
    growth_req_p_tlabile_mavg(:)       = calc_mavg(growth_req_p_tlabile_mavg(:), growth_req_p(:), &
                                                    mavg_period_tlabile, do_calc=(growth_req_p(:) > eps8))
    
    !> 1.2 averages at the time-scale of phenology
    t_air_tphen_mavg(:)           = calc_mavg(t_air_tphen_mavg(:), t_air(:), mavg_period_tphen)
    beta_soil_gs_tphen_mavg(:) = calc_mavg(beta_soil_gs_tphen_mavg(:), beta_soil_gs(:), mavg_period_tphen)

    !> 1.3 averages at the nutrient-uptake demand time-scale
    npp_tuptake_mavg(:)                 = calc_mavg(npp_tuptake_mavg(:), npp(:), mavg_period_tuptake, &
                                                     do_calc=(veg%pool%fine_root%carbon%ptr > eps8))
    demand_uptake_n_tuptake_mavg(:)     = calc_mavg(demand_uptake_n_tuptake_mavg(:), f_n_demand(:), &
                                                     mavg_period_tuptake, do_calc=(veg%pool%fine_root%carbon%ptr > eps8))
    demand_uptake_p_tuptake_mavg(:)     = calc_mavg(demand_uptake_p_tuptake_mavg(:), f_p_demand(:), &
                                                     mavg_period_tuptake, do_calc=(veg%pool%fine_root%carbon%ptr > eps8))
    cost_n_fixation_tuptake_mavg(:)     = calc_mavg(cost_n_fixation_tuptake_mavg(:), cost_n_fixation(:), &
                                                     mavg_period_tuptake, do_calc=growing_season(:)) 
    growth_req_n_tuptake_mavg(:)        = calc_mavg(growth_req_n_tuptake_mavg(:), growth_req_n(:), &
                                                     mavg_period_tuptake, do_calc=(growth_req_n(:) > eps8))
    growth_req_p_tuptake_mavg(:)        = calc_mavg(growth_req_p_tuptake_mavg(:), growth_req_p(:), &
                                                     mavg_period_tuptake, do_calc=(growth_req_p(:) > eps8))

    !> 1.4 averages at the time-scale of within-leaf N allocation fractions 
    !!
    ppfd_sunlit_tfrac_mavg_cl(:,:)     = calc_mavg(ppfd_sunlit_tfrac_mavg_cl(:,:), ppfd_sunlit_daytime_cl(:,:), &
                                                  mavg_period_tfrac)
    ppfd_shaded_tfrac_mavg_cl(:,:)     = calc_mavg(ppfd_shaded_tfrac_mavg_cl(:,:), ppfd_shaded_daytime_cl(:,:), &
                                                  mavg_period_tfrac)
    fleaf_sunlit_tfrac_mavg_cl(:,:)    = calc_mavg(fleaf_sunlit_tfrac_mavg_cl(:,:), fleaf_sunlit_daytime_cl(:,:), &
                                                  mavg_period_tfrac)
    t_air_tfrac_mavg(:)              = calc_mavg(t_air_tfrac_mavg(:), t_air_daytime(:), &
                                                  mavg_period_tfrac)
    press_srf_tfrac_mavg(:)          = calc_mavg(press_srf_tfrac_mavg(:), press_srf_daytime(:), &
                                                  mavg_period_tfrac)
    co2_mixing_ratio_tfrac_mavg(:)   = calc_mavg(co2_mixing_ratio_tfrac_mavg(:), co2_mixing_ratio_daytime(:), &
                                                  mavg_period_tfrac)
    ga_tfrac_mavg(:)                 = calc_mavg(ga_tfrac_mavg(:), ga_daytime(:), &
                                                  mavg_period_tfrac)
    beta_soil_ps_tfrac_mavg(:)    = calc_mavg(beta_soil_ps_tfrac_mavg(:), beta_soil_ps_daytime(:), &
                                                  mavg_period_tfrac)
    beta_sinklim_ps_tfrac_mavg(:)    = calc_mavg(beta_sinklim_ps_tfrac_mavg(:), beta_sinklim_ps_daytime(:), &
                                                  mavg_period_tfrac)
    beta_soil_gs_tfrac_mavg(:)    = calc_mavg(beta_soil_gs_tfrac_mavg(:), beta_soil_gs_daytime(:), &
                                                  mavg_period_tfrac)
    beta_air_tfrac_mavg(:)        = calc_mavg(beta_air_tfrac_mavg(:), beta_air_daytime(:), &
                                                  mavg_period_tfrac)

    !> 1.5 averages at the time-scale of leaf N allocation fractions 
    !!
    ppfd_sunlit_tcnl_mavg_cl(:,:)      = calc_mavg(ppfd_sunlit_tcnl_mavg_cl(:,:), ppfd_sunlit_daytime_cl(:,:), &
                                                  mavg_period_tcnl)
    ppfd_shaded_tcnl_mavg_cl(:,:)      = calc_mavg(ppfd_shaded_tcnl_mavg_cl(:,:), ppfd_shaded_daytime_cl(:,:), &
                                                  mavg_period_tcnl)
    fleaf_sunlit_tcnl_mavg_cl(:,:)     = calc_mavg(fleaf_sunlit_tcnl_mavg_cl(:,:), fleaf_sunlit_daytime_cl(:,:), &
                                                  mavg_period_tcnl)
    t_air_tcnl_mavg(:)               = calc_mavg(t_air_tcnl_mavg(:), t_air_daytime(:), &
                                                  mavg_period_tcnl)
    press_srf_tcnl_mavg(:)           = calc_mavg(press_srf_tcnl_mavg(:), press_srf_daytime(:), &
                                                  mavg_period_tcnl)
    co2_mixing_ratio_tcnl_mavg(:)    = calc_mavg(co2_mixing_ratio_tcnl_mavg(:), co2_mixing_ratio_daytime(:), &
                                                  mavg_period_tcnl)
    ga_tcnl_mavg(:)                  = calc_mavg(ga_tcnl_mavg(:), ga_daytime(:), &
                                                  mavg_period_tcnl)
    beta_soil_ps_tcnl_mavg(:)     = calc_mavg(beta_soil_ps_tcnl_mavg(:), beta_soil_ps_daytime(:), &
                                                  mavg_period_tcnl)
    beta_sinklim_ps_tcnl_mavg(:)     = calc_mavg(beta_sinklim_ps_tcnl_mavg(:), beta_sinklim_ps_daytime(:), &
                                                  mavg_period_tcnl)
    beta_soil_gs_tcnl_mavg(:)     = calc_mavg(beta_soil_gs_tcnl_mavg(:), beta_soil_gs_daytime(:), &
                                                  mavg_period_tcnl)
    beta_air_tcnl_mavg(:)         = calc_mavg(beta_air_tcnl_mavg(:), beta_air_daytime(:), &
                                                  mavg_period_tcnl)

    !> 1.6 time-scale of plant allocation calculation (veg process)
    !!
    npp_talloc_mavg(:)                  = calc_mavg(npp_talloc_mavg(:), npp(:), &
                                                     mavg_period_talloc, do_calc=growing_season(:))
    n_fixation_talloc_mavg(:)           = calc_mavg(n_fixation_talloc_mavg(:), n_fixation(:), &
                                                      mavg_period_talloc, do_calc=growing_season(:))
    unit_npp_talloc_mavg(:)             = calc_mavg(unit_npp_talloc_mavg(:), unit_npp(:), &
                                                      mavg_period_talloc, do_calc=growing_season(:))
    transpiration_talloc_mavg(:)        = calc_mavg(transpiration_talloc_mavg(:), transpiration(:), &
                                                      mavg_period_talloc, do_calc=growing_season(:))
    unit_transpiration_talloc_mavg(:)   = calc_mavg(unit_transpiration_talloc_mavg(:), unit_transpiration(:), &
                                                      mavg_period_talloc, do_calc=growing_season(:))
    dphi_talloc_mavg(:)                 = calc_mavg(dphi_talloc_mavg(:), dphi(:), &
                                                      mavg_period_talloc, do_calc=growing_season(:))
    unit_uptake_n_act_talloc_mavg(:)    = calc_mavg(unit_uptake_n_act_talloc_mavg(:), unit_uptake_n_act(:), &
                                                      mavg_period_talloc, do_calc=growing_season(:))
    unit_uptake_p_act_talloc_mavg(:)    = calc_mavg(unit_uptake_p_act_talloc_mavg(:), unit_uptake_p_act(:), &
                                                      mavg_period_talloc, do_calc=growing_season(:))
    labile_carbon_talloc_mavg(:)        = calc_mavg(labile_carbon_talloc_mavg(:), veg%pool%labile%carbon%ptr, & 
                                                      mavg_period_talloc, do_calc=growing_season(:))
    labile_nitrogen_talloc_mavg(:)      = calc_mavg(labile_nitrogen_talloc_mavg(:), veg%pool%labile%nitrogen%ptr, & 
                                                      mavg_period_talloc, do_calc=growing_season(:))
    labile_phosphorus_talloc_mavg(:)    = calc_mavg(labile_phosphorus_talloc_mavg(:), veg%pool%labile%phosphorus%ptr, & 
                                                      mavg_period_talloc, do_calc=growing_season(:))
    reserve_carbon_talloc_mavg(:)       = calc_mavg(reserve_carbon_talloc_mavg(:), veg%pool%reserve%carbon%ptr, & 
                                                      mavg_period_talloc, do_calc=growing_season(:))
    reserve_nitrogen_talloc_mavg(:)     = calc_mavg(reserve_nitrogen_talloc_mavg(:), veg%pool%reserve%nitrogen%ptr, & 
                                                      mavg_period_talloc, do_calc=growing_season(:))
    reserve_phosphorus_talloc_mavg(:)   = calc_mavg(reserve_phosphorus_talloc_mavg(:), veg%pool%reserve%phosphorus%ptr, & 
                                                      mavg_period_talloc, do_calc=growing_season(:))
    WHERE(growth_req_n(:) > 0._wp)
       hlp(:) = 1._wp/growth_req_n(:)
    ELSE WHERE
       hlp(:) = 0._wp
    END WHERE
    growth_cn_talloc_mavg(:)  = calc_mavg(growth_cn_talloc_mavg(:), hlp(:), &
                                            mavg_period_tgrowth, do_calc=growing_season(:))
    
    WHERE(growth_req_n(:) > 0._wp .AND. growth_req_p(:) > 0._wp)
       hlp(:) = 1._wp/(growth_req_n(:) * growth_req_p(:))
    ELSE WHERE
       hlp(:) = 0._wp
    END WHERE
    
    growth_cp_talloc_mavg(:)  = calc_mavg(growth_cp_talloc_mavg(:), hlp(:), &
                                            mavg_period_tgrowth, do_calc=growing_season(:))

    !> 1.7 averages at the time-scale of root life-time
    !!
    IF (veg%pool%fine_root%carbon%ptr > eps8) THEN
       hlp(:) = veg%pool%leaf%carbon%ptr / veg%pool%fine_root%carbon%ptr
    ELSE
       hlp(:) = 0._wp
    END IF
    leaf2root_troot_mavg(:)      = calc_mavg(leaf2root_troot_mavg(:), hlp(:), lctlib%tau_fine_root, &
                                               do_calc=(veg%pool%leaf%carbon%ptr > eps8) .AND. &
                                                       (veg%pool%fine_root%carbon%ptr > eps8), & 
                                               avg_period_unit='year')
    unit_npp_troot_mavg(:)       = calc_mavg(unit_npp_troot_mavg(:), unit_npp(:), lctlib%tau_fine_root, &
                                               do_calc=growing_season(:), avg_period_unit='year') 
    unit_uptake_n_pot_troot_mavg(:)  = calc_mavg(unit_uptake_n_pot_troot_mavg(:), unit_uptake_n_pot(:), lctlib%tau_fine_root, &
                                               do_calc=growing_season(:), avg_period_unit='year')
    unit_uptake_p_pot_troot_mavg(:)  = calc_mavg(unit_uptake_p_pot_troot_mavg(:), unit_uptake_p_pot(:), lctlib%tau_fine_root, &
                                               do_calc=growing_season(:), avg_period_unit='year')
    fmaint_rate_troot_mavg(:)    = calc_mavg(fmaint_rate_troot_mavg(:), fmaint_rate(:), lctlib%tau_fine_root, &
                                               do_calc=growing_season(:), avg_period_unit='year') 

    !> 1.8 averages at the vegetation dynamics time-scale 
    an_boc_tvegdyn_mavg(:)   = calc_mavg(an_boc_tvegdyn_mavg(:), net_assimilation_boc(:), mavg_period_tvegdyn, &
                                          do_calc=(veg%pool%leaf%carbon%ptr > eps8))
    !> 1.9 averages at the respiration acclimation time-scale
    t_air_tacclim_mavg(:)    = calc_mavg(t_air_tacclim_mavg(:), t_air(:), mavg_period_tacclim)
    !> 1.10 averages with a weekly, monthly timescale 
    !!
    t_air_week_mavg(:)       = calc_mavg(t_air_week_mavg(:), t_air(:), mavg_period_weekly)
    t_air_month_mavg(:)      = calc_mavg(t_air_month_mavg(:), t_air(:), mavg_period_monthly)
    t_jmax_opt_mavg(:)       = calc_mavg(t_jmax_opt_mavg(:), t_jmax_opt_daytime(:), mavg_period_weekly)




  END SUBROUTINE calc_moving_averages


END MODULE mo_shared_update_time_averages
