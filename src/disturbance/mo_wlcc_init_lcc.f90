!> Summary: wlcc (lcc due to wind) initialisation of lcc structure
!>
!>### Initialization of the lcc structure for wlcc (lcc due to wind)
!>
!>
!> Copyright (c) 2022 Julia Nabel
!>
!>#### Revision History
!>
!> First version and ongoing changes   by Julia Nabel (2022-01-06)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
MODULE mo_wlcc_init_lcc
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN
  USE mo_jsb_lcc,             ONLY: init_lcc
  USE mo_jsb_cqt_class,       ONLY: LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, &
    &                               BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: wlcc_init_lcc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_wlcc_init_lcc_class'
  CHARACTER(len=*), PARAMETER :: procname = 'wlcc'
  
  INTEGER, PARAMETER :: nr_of_active_cqts = 1 
  INTEGER, PARAMETER :: nr_of_passive_cqts = 4

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialise lcc structure for wlcc (lcc due to wind) process
  !
  SUBROUTINE wlcc_init_lcc(tile)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':wlcc_init_lcc'

    INTEGER :: active_cqts(nr_of_active_cqts)
    INTEGER :: passive_cqts(nr_of_passive_cqts)
    INTEGER :: nr_of_involved_tiles

    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: involved_child_tiles(:)
    ! -------------------------------------------------------------------------------------------------- !
    
    ! JN-TODO: reconsider what to define here and how
    IF (.NOT. tile%name .EQ. 'veg') RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    active_cqts = (/ LIVE_CARBON_CQ_TYPE /)
    passive_cqts = (/ AG_DEAD_C_CQ_TYPE, BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE /)

    nr_of_involved_tiles = tile%Get_no_of_children()
    ALLOCATE(involved_child_tiles(nr_of_involved_tiles))
    CALL tile%Get_children_names(involved_child_tiles)

    CALL init_lcc(procname, tile, active_cqts, passive_cqts, involved_child_tiles)

    DEALLOCATE(involved_child_tiles)

  END SUBROUTINE wlcc_init_lcc

#endif
END MODULE mo_wlcc_init_lcc
