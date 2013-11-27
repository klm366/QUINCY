!>
!! @brief Contains adapters to physcial schemes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_phy_schemes_iface
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

#ifdef __ICON__
  USE mo_satad, ONLY: specific_humidity => spec_humi, &
                      sat_pres_water
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: specific_humidity !, qsat

  CHARACTER(len=*), PARAMETER :: modname = 'mo_phy_schemes_iface'

!!$CONTAINS
!!$
!!$  REAL(wp) ELEMENTAL FUNCTION qsat()
!!$
!!$  END FUNCTION qsat

#endif
END MODULE mo_phy_schemes_iface
