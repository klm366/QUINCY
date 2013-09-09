!>
!! @brief Contains constants for the vegetation processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_veg_constants
#ifdef __JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER :: &
    !
    ! Parameters for the computation of canopy conductance/resistance using Eq. 3.3.2.12 in ECHAM3 manual
    & conductance_k = 0.9_wp,     & !>
    & conductance_a = 5000._wp,   & !> Jm-3
    & conductance_b = 10._wp,     & !> Wm-2
    & conductance_c = 100._wp,    & !> ms-1

    ! Parameters for albedo
    albedo_snow_can = 0.2_wp,     & !> Albedo of snow on canopy
    sky_view_factor = 1._wp         !> Sky view factor

#endif
END MODULE mo_veg_constants
