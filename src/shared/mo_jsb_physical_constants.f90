!>
!! @brief Contains physical constants
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2016-03-18)
!!
MODULE mo_jsb_physical_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  USE mo_physical_constants_iface

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER :: &
    molarMassDryAir = amd   * 1.e-3_wp, & ! 28.970e-3_wp. Mass [kg] per amount of substance [mol] of dry air
    molarMassCO2    = amco2 * 1.e-3_wp    ! 44.0095e-3_wp. Mass [kg] per amount of substance [mol] of CO2

#endif
END MODULE mo_jsb_physical_constants
