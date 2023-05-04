!>
!! @brief Contains constants for the phenology processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-09-09)
!! Julia Nabel (2018-01-10): added params for optional functionality: forest regrowth (l_forestRegrowth)
!!                           (derive maxLAI from biomass via allometric relationships)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_pheno_constants
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp

  IMPLICIT NONE
  PUBLIC

  !< photosynthetic pathways and growth forms in the QUINCY model
  !<
  !< TODO re-write as ENUM
  !<
  INTEGER, PARAMETER ::                     &
      ievergreen              = 1         , &
      isummergreen            = 2         , &
      iraingreen              = 3         , &
      iperennial              = 4

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_pheno_constants'

#endif
END MODULE mo_pheno_constants
