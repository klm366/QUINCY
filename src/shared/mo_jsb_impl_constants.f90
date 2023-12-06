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
!! @brief
!!  Contains parameter for implementation of some specific model code
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!!  First version by Jan Engel (20190114)
!!  updated against ICON-Land by Jan Engel (2023-02-24)
!!
MODULE mo_jsb_impl_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp

  IMPLICIT NONE
  PUBLIC

  ! general parameters
  INTEGER, PARAMETER :: SHORT_NAME_LEN = 10
  INTEGER, PARAMETER :: MEDIUM_NAME_LEN = 20

  ! using the functionality of LOGICAL var with REAL var
  REAL(wp), PARAMETER :: false            = 0.0_wp    !< just zero
  REAL(wp), PARAMETER :: true             = 1.0_wp    !< just one
  REAL(wp), PARAMETER :: test_false_true  = 0.5_wp    !< just 0.5, inbetween false and true, for easy IF (FALSE/TRUE)

  INTEGER,  PARAMETER :: ifs_nsoil         = 4                                  !< Number of soil levels in IFS initial file
  REAL(wp), PARAMETER :: ifs_soil_depth(4) = (/ 0.07_wp,0.28_wp,1._wp,2.89_wp/) !< Soil layer depths in IFS initial file [m]

  REAL(wp), PARAMETER :: def_parameters = -9999.0_wp  !< the default value of model parameters (not constants) prior to their init

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_jsb_impl_constants'

! CONTAINS


#endif
END MODULE mo_jsb_impl_constants
