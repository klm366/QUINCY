!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_quincy_mavg_class


  USE mo_kind,                    ONLY: lp
  USE mo_quincy_queue,            ONLY: t_queue

  IMPLICIT NONE
  PUBLIC :: t_mavg

  CHARACTER(len=*), PARAMETER :: modname = 'mo_quincy_mavg_class'


  TYPE :: t_mavg
    REAL(lp) :: dtime_sec       !< length of timestep in seconds [s]
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

  ELEMENTAL FUNCTION create_t_mavg(dtime, ndays) RESULT (new_mavg)
    REAL(lp), INTENT(in) :: dtime !< length of ts [s]
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
    REAL(lp) :: mavg !< output
    REAL(lp), ALLOCATABLE :: values(:)
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
    REAL(lp), INTENT(in)         :: newval
    REAL(lp)                     :: oldval
    IF (this % vqueue % isFull()) THEN
      ! drop oldest value
      oldval = this % vqueue % get()
    ENDIF
    call this % vqueue % put(newval)
  END SUBROUTINE add_mavg

END MODULE mo_quincy_mavg_class
