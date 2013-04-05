!>
!! @brief Contains the interface to io 
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_io_iface

#ifdef __ICON__
  USE mo_cf_convention, ONLY: t_cf_var
  USE mo_grib2,         ONLY: t_grib2_var
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_cf_var, t_grib2_var
  
  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io_iface'

END MODULE mo_jsb_io_iface
