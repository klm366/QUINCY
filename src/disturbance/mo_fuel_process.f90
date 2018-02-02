!>
!! @brief Contains the routines for the fuel processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-07-29)
!!
MODULE mo_fuel_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_fuel_jsbach

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fuel_process'

CONTAINS

  ! Calculate fuel (JSBACH3 algorithm)
  PURE ELEMENTAL SUBROUTINE calc_fuel_jsbach( &
    & Ycpool_acid_ag1,       & ! in
    & YCpool_water_ag1,      & ! in
    & YCpool_ethanol_ag1,    & ! in
    & YCpool_nonsoluble_ag1, & ! in
    & YCpool_acid_ag2,       & ! in
    & YCpool_water_ag2,      & ! in
    & YCpool_ethanol_ag2,    & ! in
    & YCpool_nonsoluble_ag2, & ! in
    & fuel                   & ! out
    & )

    ! Input Arguments
    REAL(wp), INTENT(in)  :: Ycpool_acid_ag1
    REAL(wp), INTENT(in)  :: YCpool_water_ag1
    REAL(wp), INTENT(in)  :: YCpool_ethanol_ag1
    REAL(wp), INTENT(in)  :: YCpool_nonsoluble_ag1
    REAL(wp), INTENT(in)  :: Ycpool_acid_ag2
    REAL(wp), INTENT(in)  :: YCpool_water_ag2
    REAL(wp), INTENT(in)  :: YCpool_ethanol_ag2
    REAL(wp), INTENT(in)  :: YCpool_nonsoluble_ag2

    ! Output Arguments
    REAL(wp), INTENT(OUT)   :: fuel

    ! Locals

    ! ---------------------------
    ! Go

    fuel =  Ycpool_acid_ag1    + YCpool_water_ag1      &
        & + YCpool_ethanol_ag1 + YCpool_nonsoluble_ag1 & 
        & + YCpool_acid_ag2    + YCpool_water_ag2      &
        & + YCpool_ethanol_ag2 + YCpool_nonsoluble_ag2

  END SUBROUTINE calc_fuel_jsbach

#endif
END MODULE mo_fuel_process
