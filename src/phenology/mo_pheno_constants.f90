!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., and Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
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

  USE mo_kind,                 ONLY: wp

  IMPLICIT NONE
  PUBLIC


  !< phenology parameters
  !! 
  REAL(wp), SAVE ::   & 
      gdd_t_air_threshold             !< temperature threshold for the accumulation of growing degree days [deg K]



  !< some useful numbers for easy identification of photosynthetic pathways and growth forms
  !! need to be parameters (i.e. constant across model runtime)
  !! 
  INTEGER, PARAMETER ::           &
      ievergreen    = 1         , &
      isummergreen  = 2         , &
      iraingreen    = 3         , &
      iperennial    = 4

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_pheno_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants used in the phenology process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_pheno_constants

    USE mo_jsb_phys_constants,    ONLY: Tzero

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    ! 
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_pheno_constants'

    


    gdd_t_air_threshold       = 5.0_wp + Tzero            !< ... citation


  END SUBROUTINE init_pheno_constants

#endif
END MODULE mo_pheno_constants
