!>
!! @brief fundamental constants of the land memory 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_base_constants

  USE mo_kind, ONLY : wp,eps4

!  INTEGER, PARAMETER    :: nsoil      = 15        !< number of soil layer simulated by the model
  ! convenience 
  REAL(wp), PARAMETER   :: zero       = 0.0_wp    !< just zero
  ! time control
  REAL(wp), SAVE        :: dtime                  !< length of one time step in seconds
  REAL(wp), SAVE        :: one_year               !< length of one year in days 
  REAL(wp), SAVE        :: one_day                !< length of one day in seconds
  ! conversion of degree to radians
  REAL(wp), SAVE        :: pi                     !< the number PI
  REAL(wp), SAVE        :: deg2rad                !< conversion of degree to radians
  REAL(wp), PARAMETER   :: degcircle  = 360._wp   !< degrees per circle
  REAL(wp), PARAMETER   :: hourday    = 24.0_wp   !< hours in a day

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_base_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialise base constants 
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_base_constants

  IMPLICIT NONE

    ! 0.1 InOut
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_base_constants'

    !< 1.0 degree/radians conversion
    !! 
    pi      = 4.0_wp * ATAN(1.0_wp)
    deg2rad = pi / 180.0_wp

    !< 2.0 time parameters
    !! 
    one_day   = 86400._wp
    one_year  = 365._wp

  END SUBROUTINE init_base_constants

END MODULE mo_base_constants
