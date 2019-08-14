!>
!! @par Copyright
!!  This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!!  Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
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
