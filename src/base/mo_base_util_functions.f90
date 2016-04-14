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

  !> Type for storing and calculataing averages per day for multiple variables
  TYPE :: t_daily_average
    ! the 1st index of an array can be addressed with 1 (i.e. not 0 !)
    ! averages
    REAL(wp), DIMENSION(ncanop)     :: avrg_prev_par_sunlit_day      !<  average of the previous day
    REAL(wp), DIMENSION(ncanop)     :: avrg_prev_par_shaded_day      !<  average of the previous day
    REAL(wp), DIMENSION(ncanop)     :: avrg_prev_fleaf_sunlit_day    !<  average of the previous day
    ! accumulated value of the current day & number of values accumulated
    REAL(wp), DIMENSION(ncanop, 2)  :: curr_par_sunlit_day           !<  (1) value (2) counter; acc over current day
    REAL(wp), DIMENSION(ncanop, 2)  :: curr_par_shaded_day           !<  (1) value (2) counter; acc over current day
    REAL(wp), DIMENSION(ncanop, 2)  :: curr_fleaf_sunlit_day         !<  (1) value (2) counter; acc over current day
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

  !> test if the passed REAL variable is NaN; 
  !! returns TRUE if the variable is NaN; IF (X /= X) THEN X is an NaN beacuse NaN is equal to nothing including itself 
  FUNCTION test_real_for_NaN(test_var) RESULT(l_result)
    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp), INTENT(in)      :: test_var
    ! 0.2 Output
    LOGICAL                   :: l_result

    l_result = .FALSE.
    IF(test_var /= test_var) l_result = .TRUE.
  END FUNCTION test_real_for_NaN

  !> test if the passed INTEGER variable is NaN; 
  !! returns TRUE if the variable is NaN; IF (X /= X) THEN X is an NaN beacuse NaN is equal to nothing including itself 
  FUNCTION test_int_for_NaN(test_var) RESULT(l_result)
    IMPLICIT NONE
    ! 0.1 Input
    INTEGER, INTENT(in)       :: test_var
    ! 0.2 Output
    LOGICAL                   :: l_result

    l_result = .FALSE.
    IF(test_var /= test_var) l_result = .TRUE.
  END FUNCTION test_int_for_NaN

  !> moving average; 
  !! calculates the average across a certain period preceding to the current value; 
  !! actually, the average gets rather estimated than exactly calculated using the weighted arithmetic mean 
  PURE ELEMENTAL FUNCTION moving_avrg(cav,nv,w,mask) RESULT(nav)
    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp), INTENT(in)            :: cav     !< current average value
    REAL(wp), INTENT(in)            :: nv      !< new value
    REAL(wp), INTENT(in)            :: w       !< weight of the new value
    LOGICAL, OPTIONAL, INTENT(in)   :: mask    !< bool, whether the new average is calculated from the nv, or av is used as nav
    ! 0.2 Output
    REAL(wp)                        :: nav     !< function output: the new average
    ! 0.3 Local
    LOGICAL                         :: l_mask  ! True/ False of the Optional variable mask present

    l_mask = .FALSE.
    IF(PRESENT(mask))  l_mask = mask

    IF(l_mask) THEN
      nav = cav
    ELSE
      nav = cav * (1.0_wp - w) + nv * w  ! weighted mean
    END IF

  END FUNCTION moving_avrg

  SUBROUTINE initialization_of_variables (this)
    IMPLICIT NONE 
    CLASS(t_daily_average), INTENT(inout) :: this

    this%avrg_prev_par_sunlit_day    = 0.0_wp    ! this sets the complete array to 0.0_wp
    this%avrg_prev_par_shaded_day    = 0.0_wp
    this%avrg_prev_fleaf_sunlit_day  = 0.0_wp

    this%curr_par_sunlit_day    = 0.0_wp
    this%curr_par_shaded_day    = 0.0_wp
    this%curr_fleaf_sunlit_day  = 0.0_wp

    this%l_par_sunlit_day    = .FALSE.
    this%l_par_shaded_day    = .FALSE.
    this%l_fleaf_sunlit_day  = .FALSE.

  END SUBROUTINE initialization_of_variables

  SUBROUTINE update_conditions_current_timestep (this, plant)
    USE mo_base_memory,   ONLY: PlantClass
    IMPLICIT NONE 
    TYPE(PlantClass), INTENT(in) :: plant        !< mo_base_memory: TYPE(PlantClass) :: plant; initialized in TYPE TileClass
    CLASS(t_daily_average), INTENT(inout) :: this

    ! condition type 1: light > zero
    IF(plant%canopy_layer(1)%par_sunlit + plant%canopy_layer(1)%par_shaded > eps8) THEN
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

  SUBROUTINE accumulate_values_over_current_day (this, plant)
    USE mo_base_memory,   ONLY: PlantClass
    IMPLICIT NONE 
    TYPE(PlantClass), INTENT(in) :: plant        !< mo_base_memory: TYPE(PlantClass) :: plant; initialized in TYPE TileClass
    CLASS(t_daily_average), INTENT(inout) :: this
    INTEGER  :: i

    CALL this%update_conditions(plant)

    IF(this%l_par_sunlit_day) THEN
      DO i=1,ncanop
        this%curr_par_sunlit_day(i,1) = this%curr_par_sunlit_day(i,1) + plant%canopy_layer(i)%par_sunlit
        this%curr_par_sunlit_day(i,2) = this%curr_par_sunlit_day(i,2) + 1.0_wp
      END DO
    END IF

    IF(this%l_par_shaded_day) THEN
      DO i=1,ncanop
        this%curr_par_shaded_day(i,1) = this%curr_par_shaded_day(i,1) + plant%canopy_layer(i)%par_shaded
        this%curr_par_shaded_day(i,2) = this%curr_par_shaded_day(i,2) + 1.0_wp
      END DO
    END IF

    IF(this%l_fleaf_sunlit_day) THEN
      DO i=1,ncanop
        this%curr_fleaf_sunlit_day(i,1) = this%curr_fleaf_sunlit_day(i,1) + plant%canopy_layer(i)%fleaf_sunlit
        this%curr_fleaf_sunlit_day(i,2) = this%curr_fleaf_sunlit_day(i,2) + 1.0_wp
      END DO
    END IF

  END SUBROUTINE accumulate_values_over_current_day

  !< need to be called with: CALL daily_average(:)%calc_avrg(tile(:)%plant)
  !! daily_average(:) has the same dimension like tile(:)
  SUBROUTINE calculate_avrg_previous_day (this, plant, time_day, timestep_prop_day)
    USE mo_base_memory,     ONLY: PlantClass
    IMPLICIT NONE 
    TYPE(PlantClass), INTENT(in) :: plant        !< mo_base_memory: TYPE(PlantClass) :: plant; initialized in TYPE TileClass
    REAL(wp)                     :: time_day
    REAL(wp)                     :: timestep_prop_day
    CLASS(t_daily_average), INTENT(inout) :: this
    INTEGER  :: i

    ! the accumulation over the current day
    CALL this%accumulate(plant)   ! calls this%update_conditions(plant)

    !< if this is the last timestep of the current day (i.e. next timestep belongs to the next day):
    !! actually it is tested if the next timestep would be midnight: next_step(time_day) - round_to_int(time_day) < eps8
    IF( ABS(((time_day + timestep_prop_day) - REAL(NINT(time_day)))) < eps8 ) THEN
      !< calculate the average of the current day and pass it to the according avrg_prev_* variable; 
      !! check if at least one value had been assigned at the current day, avoiding division by zero 
      DO i=1,ncanop
        ! avrg_prev_par_sunlit_day
        IF( this%curr_par_sunlit_day(i,2) < eps8 ) THEN
          this%avrg_prev_par_sunlit_day(i)    = 0.0_wp
          ELSE
          this%avrg_prev_par_sunlit_day(i)    = this%curr_par_sunlit_day(i,1) / this%curr_par_sunlit_day(i,2)
        END IF
        ! avrg_prev_par_shaded_day
        IF( this%curr_par_shaded_day(i,2) < eps8 ) THEN
          this%avrg_prev_par_shaded_day(i)    = 0.0_wp
          ELSE
          this%avrg_prev_par_shaded_day(i)    = this%curr_par_shaded_day(i,1)   / this%curr_par_shaded_day(i,2)
        END IF
        ! avrg_prev_fleaf_sunlit_day
        IF( this%curr_fleaf_sunlit_day(i,2) < eps8 ) THEN
          this%avrg_prev_fleaf_sunlit_day(i)  = 0.0_wp
          ELSE
          this%avrg_prev_fleaf_sunlit_day(i)  = this%curr_fleaf_sunlit_day(i,1) / this%curr_fleaf_sunlit_day(i,2)
        END IF
      END DO
      !test
      !write(*,*) time_day, &
      !           NINT(this%curr_par_sunlit_day(1,2)), this%avrg_prev_par_sunlit_day(1), &
      !           NINT(this%curr_par_shaded_day(1,2)), this%avrg_prev_par_shaded_day(1), &
      !           NINT(this%curr_fleaf_sunlit_day(1,2)), this%avrg_prev_fleaf_sunlit_day(1)

      !< set all values of the current day to zero (code can be copied from initialization_of_variables above)
      this%curr_par_sunlit_day    = 0.0_wp
      this%curr_par_shaded_day    = 0.0_wp
      this%curr_fleaf_sunlit_day  = 0.0_wp
    END IF

  END SUBROUTINE calculate_avrg_previous_day

END MODULE mo_base_util_functions
