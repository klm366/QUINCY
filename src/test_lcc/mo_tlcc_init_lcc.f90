!>
!! @brief Initialization of the tlcc lcc structure.
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2021-04-17)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_tlcc_init_lcc
#ifndef __NO_JSBACH__

  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN
  USE mo_jsb_lcc_class,       ONLY: t_jsb_lcc_proc
  USE mo_jsb_lcc,             ONLY: init_lcc
  USE mo_jsb_cqt_class,       ONLY: LIVE_CARBON_CQ_TYPE, DEAD_CARBON_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: tlcc_init_lcc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tlcc_init_lcc_class'
  CHARACTER(len=*), PARAMETER :: procname = 'tlcc'
  
  INTEGER, PARAMETER :: numberOfActiveCQTs = 1 
  INTEGER, PARAMETER :: numberOfPassiveCQTs = 1
  INTEGER, PARAMETER :: numberOfInvolvedChildTiles = 5

CONTAINS

  SUBROUTINE tlcc_init_lcc(tile)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    
    INTEGER :: activeCQTs(numberOfActiveCQTs)
    INTEGER :: passiveCQTs(numberOfPassiveCQTs)
    CHARACTER(len=SHORT_NAME_LEN) :: involved_child_tiles(numberOfInvolvedChildTiles)

    ! JN-TODO: reconsider what to define here and how
    IF (.NOT. tile%name .EQ. 'veg') RETURN

    activeCQTs = (/ LIVE_CARBON_CQ_TYPE /)
    passiveCQTs = (/ DEAD_CARBON_CQ_TYPE /)

    involved_child_tiles = (/ 'pft03', 'pft04' , 'pft05', 'pft06', 'pft11' /)

    CALL init_lcc(procname, tile, activeCQTs, passiveCQTs, involved_child_tiles)

  END SUBROUTINE tlcc_init_lcc

#endif
END MODULE mo_tlcc_init_lcc
