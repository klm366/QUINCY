!>
!! @brief Contains the interface to the phenology processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20161214)
!!
MODULE mo_pheno_interface

#include "dsl_jsbach.inc"

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish

  USE mo_update_phenology

  IMPLICIT NONE

!  PRIVATE
!  PUBLIC :: Register_assimi_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_interface'

!CONTAINS

  ! Task
  ! ...





END MODULE mo_pheno_interface
