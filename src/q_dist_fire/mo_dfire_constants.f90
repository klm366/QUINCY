!>
!>  @par Copyright
!>    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!>    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!>    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!>
!>    This file is part of the QUINCY model.
!>    The QUINCY model is free software: you can redistribute it and/or modify
!>    it under the terms of the GNU General Public License version 3,
!>    as published by the Free Software Foundation.
!>
!>    The QUINCY model is distributed in the hope that it will be useful,
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!>    GNU General Public License for more details.
!>
!>    You should have received a copy of the GNU General Public License
!>    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!> 
!> @brief Contains constants for the dfire process
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
!>
MODULE mo_dfire_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                 ONLY: wp

  IMPLICIT NONE

  PUBLIC

  !>
  !> dfire parameters
  !>
  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_dfire_constants'
  !> Add here model constants

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Initialize constants used in the dfire process
  !>
  SUBROUTINE init_dfire_constants
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_dfire_constants'
    ! ----------------------------------------------------------------------------------------------------- !
  END SUBROUTINE init_dfire_constants

#endif
END MODULE mo_dfire_constants
