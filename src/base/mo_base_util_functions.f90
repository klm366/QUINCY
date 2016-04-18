!>
!! @brief generic functions useful across the whole model
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-02)
!!
MODULE mo_base_util_functions

  USE mo_kind,            ONLY: wp, eps8
  USE mo_base_constants,  ONLY: ncanop, ntile

  IMPLICIT NONE
 
  PRIVATE
 
  PUBLIC :: moving_avrg, daily_average, test_real_for_NaN, test_int_for_NaN

  !-----------------------------------------------------------------------------------------------------
  !> Type for storing and calculating averages per day for multiple variables
  !!
  !! includes the PROCEDURE's necessary to do so, makes the average of the previous day globally 
  !! available;
  !! for each tile, one object is initialized
  !-----------------------------------------------------------------------------------------------------
  TYPE :: t_daily_average
    ! the 1st index of an array can be addressed with 1 (i.e. not 0 !)
    ! averages
    REAL(wp), DIMENSION(ncanop)     :: avrg_prev_par_sunlit_day      !<  average of the previous day
    REAL(wp), DIMENSION(ncanop)     :: avrg_prev_par_shaded_day      !<  average of the previous day
    REAL(wp), DIMENSION(ncanop)     :: avrg_prev_fleaf_sunlit_day    !<  average of the previous day
    ! accumulated value of the current day & number of values accumulated
    REAL(wp), DIMENSION(ncanop)     :: curr_val_par_sunlit_day        !<  (1) value accumulated over current day
    REAL(wp), DIMENSION(ncanop)     :: curr_val_par_shaded_day        !<  (1) value accumulated over current day
    REAL(wp), DIMENSION(ncanop)     :: curr_val_fleaf_sunlit_day      !<  (1) value accumulated over current day
    REAL(wp), DIMENSION(ncanop)     :: curr_count_par_sunlit_day      !<  (2) counter of acc values over current day
    REAL(wp), DIMENSION(ncanop)     :: curr_count_par_shaded_day      !<  (2) counter of acc values over current day
    REAL(wp), DIMENSION(ncanop)     :: curr_count_fleaf_sunlit_day    !<  (2) counter of acc values over current day
    ! condition per variable and timestep
    LOGICAL                         :: l_par_sunlit_day              !< whether sunlight is available at canopy top
    LOGICAL                         :: l_par_shaded_day              !< whether sunlight is available at canopy top
    LOGICAL                         :: l_fleaf_sunlit_day            !< whether sunlight is available at canopy top
  CONTAINS
    PROCEDURE :: init              => initialization_of_variables
    PROCEDURE :: update_conditions => update_conditions_current_timestep
    PROCEDURE :: accumulate        => accumulate_values_over_current_day
    PROCEDURE :: calc_avrg         => calculate_avrg_previous_day  ! at the end of calculations of the last timestep per day
  END TYPE t_daily_average

  TYPE(t_daily_average), SAVE, DIMENSION(ntile) :: daily_average 

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> test if the passed REAL variable is NaN
  !! 
  !! returns TRUE if the variable is NaN; IF (X /= X) THEN X is an NaN beacuse NaN is equal to nothing including itself 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION test_real_for_NaN(test_var) RESULT(l_result)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in)      :: test_var    !< variable to be tested
    LOGICAL                   :: l_result

    l_result = .FALSE.
    IF(test_var /= test_var) l_result = .TRUE.
  END FUNCTION test_real_for_NaN

  !-----------------------------------------------------------------------------------------------------
  !> test if the passed INTEGER variable is NaN
  !!
  !! returns TRUE if the variable is NaN; IF (X /= X) THEN X is an NaN beacuse NaN is equal to nothing including itself 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION test_int_for_NaN(test_var) RESULT(l_result)
    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in)       :: test_var    !< variable to be tested
    LOGICAL                   :: l_result

    l_result = .FALSE.
    IF(test_var /= test_var) l_result = .TRUE.
  END FUNCTION test_int_for_NaN

  !-----------------------------------------------------------------------------------------------------
  !> moving average 
  !!
  !! calculates the average across a certain period preceding to the current value; 
  !! actually, the average gets rather estimated than exactly calculated using the weighted arithmetic mean 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION moving_avrg(cav,nv,w,mask) RESULT(nav)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in)            :: cav     !< current average value
    REAL(wp), INTENT(in)            :: nv      !< new value
    REAL(wp), INTENT(in)            :: w       !< weight of the new value
    LOGICAL, OPTIONAL, INTENT(in)   :: mask    !< bool, whether the new average is calculated from the nv, or av is used as nav
    REAL(wp)                        :: nav     !< function output: the new average
    ! 0.2 Local
    LOGICAL                         :: l_mask  ! True/ False of the Optional variable mask present

    ! check whether the optional statement is provided
    l_mask = .FALSE.
    IF(PRESENT(mask))  l_mask = mask

    ! calc average as weighted mean
    IF(l_mask) THEN
      nav = cav
    ELSE
      nav = cav * (1.0_wp - w) + nv * w
    END IF
  END FUNCTION moving_avrg

  !-----------------------------------------------------------------------------------------------------
  !> intialize all variables (TYPE t_daily_average)
  !! 
  !! just set to zero or FALSE, called in the driver
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE initialization_of_variables (this)
    IMPLICIT NONE 
    CLASS(t_daily_average), INTENT(inout) :: this

    this%avrg_prev_par_sunlit_day    = 0.0_wp    ! this sets the complete array to 0.0_wp
    this%avrg_prev_par_shaded_day    = 0.0_wp
    this%avrg_prev_fleaf_sunlit_day  = 0.0_wp

    this%curr_val_par_sunlit_day      = 0.0_wp
    this%curr_val_par_shaded_day      = 0.0_wp
    this%curr_val_fleaf_sunlit_day    = 0.0_wp
    this%curr_count_par_sunlit_day    = 0.0_wp
    this%curr_count_par_shaded_day    = 0.0_wp
    this%curr_count_fleaf_sunlit_day  = 0.0_wp

    this%l_par_sunlit_day    = .FALSE.
    this%l_par_shaded_day    = .FALSE.
    this%l_fleaf_sunlit_day  = .FALSE.
  END SUBROUTINE initialization_of_variables

  !-----------------------------------------------------------------------------------------------------
  !> definition and check of conditions for each variable (TYPE t_daily_average)
  !! 
  !! for each variable a condition is defined whether the value of the current time step is used for
  !! average calculation or not;
  !! during model-runtime this subroutine is used to check the condition and set the particular
  !! variable TRUE or FALSE, valid for one timestep
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_conditions_current_timestep (this, plant)
    USE mo_base_memory,    ONLY: PlantClass
    USE mo_atmland_memory, ONLY: atm2land

    IMPLICIT NONE 
    TYPE(PlantClass), INTENT(in) :: plant          !< mo_base_memory: TYPE(PlantClass) :: plant; initialized in TYPE TileClass
    CLASS(t_daily_average), INTENT(inout) :: this  !< TYPE t_daily_average

    ! condition type 1: light > zero
    IF(atm2land%swpar_srf_down > eps8) THEN
      this%l_par_sunlit_day    = .TRUE.
      this%l_par_shaded_day    = .TRUE.
      this%l_fleaf_sunlit_day  = .TRUE.
    ELSE
      this%l_par_sunlit_day    = .FALSE.
      this%l_par_shaded_day    = .FALSE.
      this%l_fleaf_sunlit_day  = .FALSE.      
    END IF
    ! condition type 2:
        ! ... currently none
  END SUBROUTINE update_conditions_current_timestep

  !-----------------------------------------------------------------------------------------------------
  !> accumulation of values over one day (TYPE t_daily_average)
  !! 
  !! per variable the values are summed up and the number of summed values is stored
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE accumulate_values_over_current_day (this, plant)
    USE mo_base_memory,   ONLY: PlantClass
    IMPLICIT NONE 
    ! 0.1 InOut
    TYPE(PlantClass), INTENT(in) :: plant          !< mo_base_memory: TYPE(PlantClass) :: plant; initialized in TYPE TileClass
    CLASS(t_daily_average), INTENT(inout) :: this  !< TYPE t_daily_average
    ! 0.2 Local
    INTEGER  :: i                                  ! loop over canopy layers

    CALL this%update_conditions(plant)

    IF(this%l_par_sunlit_day) THEN
      DO i=1,ncanop
        this%curr_val_par_sunlit_day(i) = this%curr_val_par_sunlit_day(i) + plant%canopy_layer(i)%par_sunlit
        this%curr_count_par_sunlit_day(i) = this%curr_count_par_sunlit_day(i) + 1.0_wp
      END DO
    END IF

    IF(this%l_par_shaded_day) THEN
      DO i=1,ncanop
        this%curr_val_par_shaded_day(i) = this%curr_val_par_shaded_day(i) + plant%canopy_layer(i)%par_shaded
        this%curr_count_par_shaded_day(i) = this%curr_count_par_shaded_day(i) + 1.0_wp
      END DO
    END IF

    IF(this%l_fleaf_sunlit_day) THEN
      DO i=1,ncanop
        this%curr_val_fleaf_sunlit_day(i) = this%curr_val_fleaf_sunlit_day(i) + plant%canopy_layer(i)%fleaf_sunlit
        this%curr_count_fleaf_sunlit_day(i) = this%curr_count_fleaf_sunlit_day(i) + 1.0_wp
      END DO
    END IF

  END SUBROUTINE accumulate_values_over_current_day

  !-----------------------------------------------------------------------------------------------------
  !> calc the average of variables for the previous day (TYPE t_daily_average)
  !! 
  !! calculates the average at the end of the day (end of last time step of day), and sets all other
  !! variables to zero then
  !! 
  !! is called in the driver with: CALL daily_average(:)%calc_avrg(tile(:)%plant)
  !! daily_average(:) has the same dimension as tile(:)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calculate_avrg_previous_day (this, plant, time_day, timestep_prop_day)
    USE mo_base_memory,     ONLY: PlantClass
    IMPLICIT NONE 
    ! 0.1 InOut
    TYPE(PlantClass), INTENT(in) :: plant        !< mo_base_memory: TYPE(PlantClass) :: plant; initialized in TYPE TileClass
    REAL(wp), INTENT(in)         :: time_day           !< time in days
    REAL(wp), INTENT(in)         :: timestep_prop_day  !< length of one time step (proportion of one day)
    CLASS(t_daily_average), INTENT(inout) :: this      !< TYPE t_daily_average
    ! 0.2 Local
    INTEGER  :: i                                   ! loop over canopy layers

    ! the accumulation over the current day
    CALL this%accumulate(plant)   ! calls this%update_conditions(plant)

    ! if this is the last timestep of the current day (i.e. next timestep belongs to the next day):
    ! actually it is tested if the next timestep would be midnight: next_step(time_day) - round_to_int(time_day) < eps8
    IF( ABS(((time_day + timestep_prop_day) - REAL(NINT(time_day)))) < eps8 ) THEN
      ! calculate the average of the current day and pass it to the according avrg_prev_* variable; 
      ! check if at least one value had been assigned at the current day, avoiding division by zero 
      DO i=1,ncanop
        ! avrg_prev_par_sunlit_day
        IF( this%curr_count_par_sunlit_day(i) < eps8 ) THEN
          this%avrg_prev_par_sunlit_day(i)    = 0.0_wp
          ELSE
          this%avrg_prev_par_sunlit_day(i)    = this%curr_val_par_sunlit_day(i) / this%curr_count_par_sunlit_day(i)
        END IF
        ! avrg_prev_par_shaded_day
        IF( this%curr_count_par_shaded_day(i) < eps8 ) THEN
          this%avrg_prev_par_shaded_day(i)    = 0.0_wp
          ELSE
          this%avrg_prev_par_shaded_day(i)    = this%curr_val_par_shaded_day(i)   / this%curr_count_par_shaded_day(i)
        END IF
        ! avrg_prev_fleaf_sunlit_day
        IF( this%curr_count_fleaf_sunlit_day(i) < eps8 ) THEN
          this%avrg_prev_fleaf_sunlit_day(i)  = 0.0_wp
          ELSE
          this%avrg_prev_fleaf_sunlit_day(i)  = this%curr_val_fleaf_sunlit_day(i) / this%curr_count_fleaf_sunlit_day(i)
        END IF
      END DO
      
      !write(*,*) time_day, &   ! may enable for testing
      !           NINT(this%curr_count_par_sunlit_day(1)), this%avrg_prev_par_sunlit_day(1), &
      !           NINT(this%curr_count_par_shaded_day(1)), this%avrg_prev_par_shaded_day(1), &
      !           NINT(this%curr_count_fleaf_sunlit_day(1)), this%avrg_prev_fleaf_sunlit_day(1)

      ! set all values of the current day to zero (code can be copied from initialization_of_variables above)
      this%curr_val_par_sunlit_day      = 0.0_wp
      this%curr_val_par_shaded_day      = 0.0_wp
      this%curr_val_fleaf_sunlit_day    = 0.0_wp
      this%curr_count_par_sunlit_day    = 0.0_wp
      this%curr_count_par_shaded_day    = 0.0_wp
      this%curr_count_fleaf_sunlit_day  = 0.0_wp
    END IF

  END SUBROUTINE calculate_avrg_previous_day

END MODULE mo_base_util_functions
