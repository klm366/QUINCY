!>
!! @brief Contains constants for the jedi processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-15)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jedi_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_jsb_physical_constants, ONLY: amco2  !< Molar weight of CO2 [g/mol]

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER :: molarMassC_kg           = 12.01_wp * 1.e-3_wp  !< 12.01e-3_wp   -  Mass of 1 mod C   in kg
  REAL(wp), PARAMETER :: molarMassCO2_kg         = amco2 * 1.e-3_wp     !< 44.0095e-3_wp -  Mass of 1 mol CO2 in kg

  REAL(wp), PARAMETER :: sec_per_day             = 86400._wp            !< seconds per day
  REAL(wp), PARAMETER :: days_per_year           = 365.25_wp
  REAL(wp), PARAMETER :: sec_per_year            = days_per_year * sec_per_day
  
  REAL(wp), PARAMETER :: max_pipe_transpir = 2.33_wp / 86400._wp

#endif
END MODULE mo_jedi_constants
