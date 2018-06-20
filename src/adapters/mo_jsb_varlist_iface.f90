!>
!! @brief Contains the interface to varlists/streams 
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_varlist_iface
#ifndef __NO_JSBACH__

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: VARNAME_LEN

  INTEGER, PARAMETER :: VARNAME_LEN = 64

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_varlist_iface'

#endif
END MODULE mo_jsb_varlist_iface
