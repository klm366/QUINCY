!>
!! @brief Contains the routines for the fuel processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-07-29)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
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
    & c_acid_ag1_ta,                          & ! in
    & c_water_ag1_ta,                         & ! in
    & c_ethanol_ag1_ta,                       & ! in
    & c_nonsoluble_ag1_ta,                    & ! in
    & c_acid_ag2_ta,                          & ! in
    & c_water_ag2_ta,                         & ! in
    & c_ethanol_ag2_ta,                       & ! in
    & c_nonsoluble_ag2_ta,                    & ! in
    & fuel                                    & ! out
    & )

    ! Input Arguments
    REAL(wp), INTENT(in)  :: c_acid_ag1_ta
    REAL(wp), INTENT(in)  :: c_water_ag1_ta
    REAL(wp), INTENT(in)  :: c_ethanol_ag1_ta
    REAL(wp), INTENT(in)  :: c_nonsoluble_ag1_ta
    REAL(wp), INTENT(in)  :: c_acid_ag2_ta
    REAL(wp), INTENT(in)  :: c_water_ag2_ta
    REAL(wp), INTENT(in)  :: c_ethanol_ag2_ta
    REAL(wp), INTENT(in)  :: c_nonsoluble_ag2_ta

    ! Output Arguments
    REAL(wp), INTENT(OUT)   :: fuel

    ! Locals

    ! ---------------------------
    ! Go

    fuel =      c_acid_ag1_ta    + c_water_ag1_ta       &
      &       + c_ethanol_ag1_ta + c_nonsoluble_ag1_ta  &
      &       + c_acid_ag2_ta    + c_water_ag2_ta       &
      &       + c_ethanol_ag2_ta + c_nonsoluble_ag2_ta

  END SUBROUTINE calc_fuel_jsbach

#endif
END MODULE mo_fuel_process
