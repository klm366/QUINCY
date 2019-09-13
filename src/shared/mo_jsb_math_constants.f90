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
!!  Contains mathematical constants
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!!  First version by Jan Engel copying parameter from mo_base_constants (20190114)
!!
MODULE mo_jsb_math_constants
 
  USE mo_kind,                  ONLY: wp

  IMPLICIT NONE
  PUBLIC


  ! time control
  REAL(wp), SAVE        :: dtime                  !< length of one time step in seconds
  REAL(wp), SAVE        :: one_year               !< length of one year in days 
  REAL(wp), SAVE        :: one_day                !< length of one day in seconds
  ! conversion of degree to radians
  REAL(wp), SAVE        :: pi                     !< the number PI
  REAL(wp), SAVE        :: deg2rad                !< conversion of degree to radians
  REAL(wp), PARAMETER   :: degcircle  = 360._wp   !< degrees per circle
  REAL(wp), PARAMETER   :: hourday    = 24.0_wp   !< hours per day
  ! QUINCY math helpers
  REAL(wp), PARAMETER   :: eps12 = 0.000000000001_wp
  REAL(wp), PARAMETER   :: eps8  = 0.00000001_wp
  REAL(wp), PARAMETER   :: eps4  = 0.0001_wp
  REAL(wp), PARAMETER   :: zero  = 0.0_wp           !< just zero
  

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_jsb_math_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize mathematical parameters
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_jsb_math_constants

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    ! 
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_jsb_math_constants'


    !< 1.0 time parameters
    !! 
    !! dtime is read from namelist into model%config%dtime and initialized with quincy_init()
    one_day   = 86400._wp
    one_year  = 365._wp
    
    !< 2.0 degree/radians conversion
    !! 
    pi      = 4.0_wp * ATAN(1.0_wp)
    deg2rad = pi / 180.0_wp


  END SUBROUTINE init_jsb_math_constants

END MODULE mo_jsb_math_constants
