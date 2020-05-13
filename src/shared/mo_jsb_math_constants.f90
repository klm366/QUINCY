!! @brief
!!  Contains mathematical constants
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!!  First version by Jan Engel copying parameter from mo_base_constants (20190114)
!!
MODULE mo_jsb_math_constants
 
  USE mo_kind,                  ONLY: wp
  USE mo_math_constants

  IMPLICIT NONE
  PUBLIC


  ! time control
  REAL(wp), PARAMETER   :: dtime     = 1800._wp               !< length of one time step in seconds
  REAL(wp), PARAMETER   :: one_year  = 365._wp                !< length of one year in days 
  REAL(wp), PARAMETER   :: one_day   = 86400._wp              !< length of one day in seconds
  ! conversion of degree to radians
  REAL(wp), PARAMETER   :: degcircle = 360._wp                !< degrees per circle
  REAL(wp), PARAMETER   :: hourday   = 24.0_wp                !< hours per day
  ! QUINCY math helpers
  REAL(wp), PARAMETER   :: eps12     = 0.000000000001_wp
  REAL(wp), PARAMETER   :: eps8      = 0.00000001_wp
  REAL(wp), PARAMETER   :: eps4      = 0.0001_wp
  REAL(wp), PARAMETER   :: zero      = 0.0_wp                 !< just zero
  

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_jsb_math_constants'

!CONTAINS


END MODULE mo_jsb_math_constants
