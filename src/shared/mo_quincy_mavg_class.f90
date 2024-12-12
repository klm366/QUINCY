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
!!    Pure moving average
!!
!! @author
!!    Albert Jornet Puig, Leipzig University / MPI-BGC Jena
!!    Ana Bastos, Leipzig University / MPI-BGC Jena
!!
MODULE mo_quincy_mavg_class


  USE mo_kind,                    ONLY: wp
  USE mo_quincy_queue,            ONLY: t_queue

  IMPLICIT NONE
  PUBLIC :: t_mavg

  CHARACTER(len=*), PARAMETER :: modname = 'mo_quincy_mavg_class'


  TYPE :: t_mavg
    REAL(wp) :: dtime_sec       !< length of timestep in seconds [s]
    INTEGER :: ndays            !< how many days to capture [days]
    INTEGER :: ntimesteps       !< how many timesteps per day

    TYPE(t_queue) :: vqueue     !< data
  CONTAINS
    !< constructor (length, dtime)
    PROCEDURE :: add => add_mavg !< add()
    PROCEDURE :: get => get_mavg !< get()
  END TYPE t_mavg

  INTERFACE t_mavg
    MODULE PROCEDURE :: create_t_mavg
  END INTERFACE t_mavg

CONTAINS

  FUNCTION create_t_mavg(dtime, ndays) RESULT (new_mavg)
    REAL(wp), INTENT(in) :: dtime !< length of ts [s]
    INTEGER, INTENT(in) :: ndays !< how many days [days]
    TYPE(t_mavg) :: new_mavg
    new_mavg % dtime_sec = dtime
    new_mavg % ndays = ndays
    new_mavg % ntimesteps = ndays * 3600 * 24 / dtime
    new_mavg % vqueue = t_queue(new_mavg % ntimesteps)
  END FUNCTION create_t_mavg

  ! calculate moving average
  FUNCTION get_mavg(this) RESULT(mavg)
    CLASS(t_mavg), INTENT(in) :: this
    REAL(wp) :: mavg !< output
    REAL(wp), ALLOCATABLE :: values(:)
    INTEGER :: nvalues
    IF (this % vqueue % isEmpty()) THEN
      write(*,*) modname//":: get_mavg:: error: no values to calculate mavg"
      STOP
    ENDIF
    values = this % vqueue % getAll()
    nvalues = size(values)
    mavg = SUM(values) / nvalues
  END FUNCTION get_mavg

  ! add new value into the moving average
  SUBROUTINE add_mavg(this, newval)
    CLASS(t_mavg), INTENT(inout) :: this
    REAL(wp), INTENT(in)         :: newval
    REAL(wp)                     :: oldval
    IF (this % vqueue % isFull()) THEN
      ! drop oldest value
      oldval = this % vqueue % get()
    ENDIF
    call this % vqueue % put(newval)
  END SUBROUTINE add_mavg

END MODULE mo_quincy_mavg_class
