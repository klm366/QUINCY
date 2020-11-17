!>
!! @par Copyright
!!  This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!!  Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
!! @brief Contains interface to spatial domains from ICON/ECHAM6
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                              by Reiner Schnur (2013-04-04)
!! New subroutines and some adaptations for jsb4 standalone   by Julia Nabel (2016-02-20)
!!
MODULE mo_jsb_domain_iface
#ifndef __NO_JSBACH__

  USE mo_decomposition,      ONLY: pe_decomposed

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_patch, pe_decomposed

#ifdef __ICON__
  ! Dummy type for ECHAM decomposition
  TYPE pe_decomposed
    INTEGER :: nlon = -1
  END TYPE pe_decomposed

  TYPE(pe_decomposed), SAVE :: global_decomposition(1)
#else
  ! Dummy type for ICON patch
  TYPE t_patch
     INTEGER :: id = -1
  END type t_patch

!  TYPE(t_patch), TARGET, SAVE :: p_patch(1)
#endif

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_domain_iface'

#endif
END MODULE mo_jsb_domain_iface
