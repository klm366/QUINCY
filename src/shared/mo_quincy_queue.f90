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

  ELEMENTAL FUNCTION create_t_queue(capacity) RESULT (new_queue)
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
