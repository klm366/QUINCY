!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Contains constants for the phenology process
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20161214)
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
