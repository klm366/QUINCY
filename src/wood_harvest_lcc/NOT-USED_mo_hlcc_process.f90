! Meri's to-be-deleted notes
! "dummy" process for wood harvesting for the initial testing
! ultimately change to read data ad deal with that in the init (just like alcc)
! Then I think the process module won't even be needed
! -> ACTUALLY NO REASON TO USE THIS EVEN NOW -> HARVESTED_FRAC DIRECTLY IN INTERFACE (only kept just in case)
!
! NOW: this should just harvests 20% in the beginning of each year
!-----------------------------------------------------

MODULE mo_hlcc_process
    #ifndef __NO_JSBACH__
    
      USE mo_kind,      ONLY: wp
      USE mo_exception, ONLY: finish, message_text
      USE mo_jsb_time,  ONLY: is_newyear !is_newday, get_year, get_year_length

    
      IMPLICIT NONE
      PRIVATE
    
      PUBLIC :: harvested_fract_jsbach
    
      CHARACTER(len=*), PARAMETER :: modname = 'mo_hlcc_process'
    
    CONTAINS
    
      PURE ELEMENTAL SUBROUTINE harvested_fract_jsbach( &
        & dtime,                                        & ! in
        & current_datetime,                             & ! in
        & harvested_fract                               & ! inout
        & )
    
    ! Input Arguments
        REAL(wp), INTENT(in) ::     &
        & dtime,                    &
        & current_datetime


        ! Output Arguments
        REAL(wp), INTENT(inout) :: harvested_fract

        !"Dummy harvesting" harvesting in the beginning of each year.

        IF(is_newyear(current_datetime,dtime)) THEN !MR-TODO : This could hit a snag if current_datetime not as expected (other scritps: options%current_datetime)
          harvested_fract = 0.2
        END IF


    
      END SUBROUTINE harvested_fract_jsbach
    
    #endif
    END MODULE mo_hlcc_process