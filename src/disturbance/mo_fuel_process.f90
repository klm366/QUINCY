!>
!! @brief Contains the routines for the fuel processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-07-29)
!! 2018-04-25, StW: Changed routine name for consistency with namelist settings
!!
MODULE mo_fuel_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_fuel_jsb3

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fuel_process'

CONTAINS

  ! Calculate fuel (JSBACH3 algorithm)
  PURE ELEMENTAL SUBROUTINE calc_fuel_jsb3( &
    & c_acid_ag1,       & ! in
    & c_water_ag1,      & ! in
    & c_ethanol_ag1,    & ! in
    & c_nonsoluble_ag1, & ! in
    & c_acid_ag2,       & ! in
    & c_water_ag2,      & ! in
    & c_ethanol_ag2,    & ! in
    & c_nonsoluble_ag2, & ! in
    & fuel                   & ! out
    & )

    ! Input Arguments
    REAL(wp), INTENT(in)  :: c_acid_ag1
    REAL(wp), INTENT(in)  :: c_water_ag1
    REAL(wp), INTENT(in)  :: c_ethanol_ag1
    REAL(wp), INTENT(in)  :: c_nonsoluble_ag1
    REAL(wp), INTENT(in)  :: c_acid_ag2
    REAL(wp), INTENT(in)  :: c_water_ag2
    REAL(wp), INTENT(in)  :: c_ethanol_ag2
    REAL(wp), INTENT(in)  :: c_nonsoluble_ag2

    ! Output Arguments
    REAL(wp), INTENT(OUT)   :: fuel

    ! Locals

    ! ---------------------------
    ! Go

    fuel =  c_acid_ag1    + c_water_ag1      &
        & + c_ethanol_ag1 + c_nonsoluble_ag1 & 
        & + c_acid_ag2    + c_water_ag2      &
        & + c_ethanol_ag2 + c_nonsoluble_ag2

  END SUBROUTINE calc_fuel_jsb3

#endif
END MODULE mo_fuel_process
