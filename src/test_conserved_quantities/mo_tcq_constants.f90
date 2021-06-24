!> Summary: constants for the tcq process (test conserved quantities)
!>
!> Copyright (c) 2021 Julia Nabel
!>
!>#### Revision History
!> First version   by Julia Nabel (2021-06-02)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
MODULE mo_tcq_constants
#ifndef __NO_JSBACH__

  USE mo_jsb_varlist,            ONLY: VARNAME_LEN

  IMPLICIT NONE
  PUBLIC

  ! Definitions for lcc

  !> list of potential active variables for the tcq process: For these variables the tcq interface provides active relocation
  CHARACTER(len=VARNAME_LEN) :: tcq_potential_active_vars(1) = [character(len=VARNAME_LEN) :: 'a_veg_c_tcq' ]

  !> list of required passive vars if vars are active (because these are used to relocate the active variables)
  CHARACTER(len=VARNAME_LEN) :: tcq_required_passive_vars(1) = [character(len=VARNAME_LEN) :: 'a_dead_c_tcq' ]

#endif
END MODULE mo_tcq_constants
