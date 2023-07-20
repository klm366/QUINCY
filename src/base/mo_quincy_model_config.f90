
!> Summary: Handels the qmodel id (representing quincy slm model set-ups)
!>
!>### Provides qmodel ids and routines to map quincy slm model names and qmodel ids
!>
!> Copyright (c) 2021 Julia Nabel
!>
!>#### Revision History
!>
!> First version   by Julia Nabel   (2021-10-13)
!>
!>#### License
!>
!> QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!> Copyright (C) 2021 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!>
!> This file is part of the QUINCY model.
!> The QUINCY model is free software: you can redistribute it and/or modify
!> it under the terms of the GNU General Public License version 3,
!> as published by the Free Software Foundation.
!>
!> The QUINCY model is distributed in the hope that it will be useful,
!> but WITHOUT ANY WARRANTY; without even the implied warranty of
!> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!> GNU General Public License for more details.
!>
!> You should have received a copy of the GNU General Public License
!> along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!>
MODULE mo_quincy_model_config

  USE mo_util,             ONLY: tolower
  USE mo_exception,        ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Get_quincy_model_config_name, Get_quincy_model_config_id
  PUBLIC :: QLAND, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION

  !> Available quincy model ids
  !>
  ENUM, BIND(C)
    ENUMERATOR :: QLAND=1, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = ' mo_quincy_model_config'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Returns the QUINCY model name for this id
  !
  FUNCTION Get_quincy_model_config_name(id) RESULT(return_value)

    INTEGER, INTENT(in) :: id
    CHARACTER(len=:), ALLOCATABLE :: return_value

    SELECT CASE(id)
    CASE (QLAND)
      return_value = 'land'
    CASE (QPLANT)
      return_value = 'plant'
    CASE (QSOIL)
      return_value = 'soil'
    CASE (QCANOPY)
      return_value = 'canopy'
    CASE (Q_TEST_CANOPY)
      return_value = 'test_canopy'
    CASE (Q_TEST_RADIATION)
      return_value = 'test_radiation'
    CASE DEFAULT
      return_value = ''
    END SELECT

  END FUNCTION Get_quincy_model_config_name

  ! ====================================================================================================== !
  !
  !> Returns the id for this QUINCY model name
  !
  FUNCTION Get_quincy_model_config_id(name) RESULT(return_value)

    CHARACTER(len=*), INTENT(in) :: name
    INTEGER                      :: return_value

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_quincy_model_config_id'

    SELECT CASE(tolower(TRIM(name)))  
    CASE ('land')
      return_value = QLAND
    CASE ('plant')
      return_value = QPLANT
    CASE ('soil')
      return_value = QSOIL
    CASE ('canopy')
      return_value = QCANOPY
    CASE ('test_canopy')
      return_value = Q_TEST_CANOPY
    CASE ('test_radiation')
      return_value = Q_TEST_RADIATION
    CASE DEFAULT
      CALL finish(routine,'Invalid QUINCY model name. Try either land, plant, soil, canopy, test_canopy or test_radiation.')
    END SELECT

  END FUNCTION Get_quincy_model_config_id

END MODULE mo_quincy_model_config
