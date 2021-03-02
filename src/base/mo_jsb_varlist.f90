!>
!! @brief Contains routines for var_lists/streams
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-12-04)
!! Copied from JSBACH4 (git revison e3001993a1da854518008510d6947611a9a942af) by Jan Engel and removed most of the code (2021-02-23)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_varlist

  USE mo_kind,              ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: BASIC, MEDIUM, FULL

  ENUM, BIND(C)
    ENUMERATOR :: BASIC=1, MEDIUM, FULL
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_varlist'

!CONTAINS

END MODULE mo_jsb_varlist
