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
!!    Implementation of a queue structure
!!
!! @author
!!    Albert Jornet Puig, Leipzig University / MPI-BGC Jena
!!    Ana Bastos, Leipzig University / MPI-BGC Jena
!!
MODULE mo_quincy_queue

  USE mo_kind,                    ONLY: wp

  IMPLICIT NONE
  PUBLIC :: t_queue

  CHARACTER(len=*), PARAMETER :: modname = 'mo_quincy_queue'

  ! array queue
  TYPE :: t_queue
    REAL(wp), ALLOCATABLE :: elements(:)
    INTEGER :: ifront, irear
    INTEGER :: capacity
  CONTAINS
    PROCEDURE :: put !put(item) : This will put the item inside the queue.
    PROCEDURE :: get !get(): This will return you an item from the queue.
    PROCEDURE :: getAll !getAll(): return all queue elements as an array
    PROCEDURE :: isEmpty !isEmpty(): It will return true if the queue is empty and false if items are present.
    PROCEDURE :: getQsize !getQsize(): returns the size of the queue.
    PROCEDURE :: isFull !isFull(): returns true if the queue is full, otherwise false.
  END TYPE t_queue


  INTERFACE t_queue
    MODULE PROCEDURE :: create_t_queue
  END INTERFACE t_queue

CONTAINS

  FUNCTION create_t_queue(capacity) RESULT (new_queue)
    INTEGER, INTENT(in) :: capacity
    TYPE(t_queue) :: new_queue !< output

    new_queue%capacity = capacity
    new_queue%ifront = 0
    new_queue%irear = 0
    ALLOCATE(new_queue%elements(capacity))
  END FUNCTION create_t_queue


  SUBROUTINE put(this, item)
    CLASS(t_queue), INTENT(inout) :: this
    REAL(wp), INTENT(in) :: item
    !write(*,*) "mo_quincy_mavg_class:: put:: new item ", item, "..."

    IF (this % isFull()) THEN
      WRITE(*,*) modname//":: put:: error: queue is full"
      STOP
    ENDIF

    this % irear = this % irear + 1

    this % elements(this%irear) = item
    !write(*,*) "mo_quincy_mavg_class:: put:: done"
  END SUBROUTINE put


  FUNCTION get(this) RESULT (item)
    CLASS(t_queue), INTENT(inout) :: this
    REAL(wp) :: item !< output

    !!write(*,*) "mo_quincy_mavg_class:: get:: starting..."

    !IF (this % ifront == this % irear) THEN
    IF (this % isEmpty()) THEN
      WRITE(*,*) modname//":: get:: error: queue is empty"
      STOP
    ENDIF

    item = this % elements(1) !< always pick 1st position

    this % elements = cshift(this % elements, 1)
    this % irear = this % irear - 1
  END FUNCTION get


  INTEGER FUNCTION getQsize(this)
    CLASS(t_queue), INTENT(in) :: this
    getQsize = this % irear
  END FUNCTION getQsize


  LOGICAL FUNCTION isEmpty(this)
    CLASS(t_queue), INTENT(in) :: this
    isEmpty = (0 == this % irear)
  END FUNCTION isEmpty


  LOGICAL FUNCTION isFull(this)
    CLASS(t_queue), INTENT(in) :: this
    isFull = (this % irear == this % capacity)
  END FUNCTION isFull


  FUNCTION getAll(this) RESULT (allelems)
    CLASS(t_queue), INTENT(in) :: this
    REAL(wp), ALLOCATABLE :: allelems(:)
    !!is empty, then error
    IF (this % isEmpty()) THEN
      WRITE(*,*) modname//":: getAll:: error: queue is empty"
      STOP
    ENDIF
    ALLOCATE(allelems, source = this % elements(1:this%irear))
  END FUNCTION getAll

END MODULE mo_quincy_queue
