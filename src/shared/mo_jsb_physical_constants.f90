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
    & dens_soil          = 2700._wp,   & !< Density of soil particulate (Hillel 1980, Johansen 1977)           [kg/m3]
    & dens_org           = 1300._wp,   & !< Density of organic part of soil                                    [kg/m3]
    & dens_ice           = rhoi,       & !< Density of frozen water in soil (ice)                              [kg/m3]
    & dens_snow          = rhos,       & !< =300.0  Density of snow                                            [kg/m3]
    & dens_snow_min      = 50._wp,     & !< Minimum density of snow (for fresh snow)                           [kg/m3]
    & dens_snow_max      = 300._wp       !< Maximum density of snow                                            [kg/m3]

  REAL(wp), PARAMETER :: &
    molarMassDryAir = amd   * 1.e-3_wp, & ! 28.970e-3_wp. Mass [kg] per amount of substance [mol] of dry air
    molarMassCO2    = amco2 * 1.e-3_wp    ! 44.0095e-3_wp. Mass [kg] per amount of substance [mol] of CO2

#endif
END MODULE mo_jsb_physical_constants
