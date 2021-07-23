!> Summary: alcc (anthropogenic lcc) lcc structure initialisation
!>
!>### Initialization of the alcc lcc structure
!>
!>
!> Copyright (c) 2021 Julia Nabel
!>
!>#### Revision History
!> First version   by Julia Nabel (2021-06-23)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
MODULE mo_alcc_init_lcc
#ifndef __NO_JSBACH__

  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN
  USE mo_jsb_lcc,             ONLY: init_lcc
  USE mo_jsb_cqt_class,       ONLY: LIVE_CARBON_CQ_TYPE, DEAD_CARBON_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: alcc_init_lcc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_alcc_init_lcc_class'
  CHARACTER(len=*), PARAMETER :: procname = 'alcc'
  
  INTEGER, PARAMETER :: nr_of_active_cqts = 1 
  INTEGER, PARAMETER :: nr_of_passive_cqts = 2

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialise lcc structure for alcc (anthropogenic lcc) process
  !
  SUBROUTINE alcc_init_lcc(tile)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract),  POINTER :: current_tile
    INTEGER :: active_cqts(nr_of_active_cqts)
    INTEGER :: passive_cqts(nr_of_passive_cqts)
    INTEGER :: nr_of_involved_tiles, i_tile

    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: involved_child_tiles(:)
    ! -------------------------------------------------------------------------------------------------- !
    
    ! JN-TODO: reconsider what to define here and how
    IF (.NOT. tile%name .EQ. 'veg') RETURN

    active_cqts = (/ LIVE_CARBON_CQ_TYPE /)
    passive_cqts = (/ PRODUCT_CARBON_CQ_TYPE, DEAD_CARBON_CQ_TYPE /)

    ! for alcc all child tiles of the veg tile are involved

    nr_of_involved_tiles = tile%Get_no_of_children()
    ALLOCATE(involved_child_tiles(nr_of_involved_tiles))

    current_tile => tile%Get_first_child_tile()
    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile))
      i_tile = i_tile + 1
      involved_child_tiles(i_tile) = current_tile%name
      current_tile => current_tile%Get_next_sibling_tile()
    ENDDO

    CALL init_lcc(procname, tile, active_cqts, passive_cqts, involved_child_tiles)

    DEALLOCATE(involved_child_tiles)

  END SUBROUTINE alcc_init_lcc

#endif
END MODULE mo_alcc_init_lcc
