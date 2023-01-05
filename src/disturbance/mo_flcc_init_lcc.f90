!> Summary: flcc (lcc due to fire) initialisation of lcc structure
!>
!>### Initialization of the lcc structure for flcc (lcc due to fire)
!>
!>
!> Copyright (c) 2021 Julia Nabel
!>
!>#### Revision History
!>
!> - First version and ongoing changes   by Julia Nabel (2021-10-29)
!> - added forest age classes            by Julia Nabel (2022-10-06)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
MODULE mo_flcc_init_lcc
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN
  USE mo_jsb_lcc,             ONLY: init_lcc, count_descendants_for_up_to_2layers, &
    &                               collect_names_of_descendant_leaves_up_to_2layers
  USE mo_jsb_cqt_class,       ONLY: LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, &
    &                               BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE

  dsl4jsb_Use_processes FLCC_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: flcc_init_lcc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_flcc_init_lcc_class'
  CHARACTER(len=*), PARAMETER :: procname = 'flcc'
  
  INTEGER, PARAMETER :: nr_of_active_cqts = 2 
  INTEGER, PARAMETER :: nr_of_passive_cqts = 4

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialise lcc structure for flcc (lcc due to fire) process
  !
  SUBROUTINE flcc_init_lcc(tile)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: involved_tiles(:)
    INTEGER :: nr_of_descendants

    INTEGER :: active_cqts(nr_of_active_cqts)
    INTEGER :: passive_cqts(nr_of_passive_cqts)

    CHARACTER(len=*), PARAMETER :: routine = modname//':flcc_init_lcc'
    ! -------------------------------------------------------------------------------------------------- !
    
    IF (.NOT. tile%Is_process_calculated(FLCC_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    active_cqts = (/ LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE /)
    passive_cqts = (/ AG_DEAD_C_CQ_TYPE, BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE /)

    ! FLCC burns area on the leave tiles
    IF(tile%Has_children()) THEN
      CALL count_descendants_for_up_to_2layers(tile, nr_of_descendants)
      ALLOCATE(involved_tiles(nr_of_descendants))
      CALL collect_names_of_descendant_leaves_up_to_2layers(tile, involved_tiles)
    ELSE
      ALLOCATE(involved_tiles(1))
      involved_tiles(1) = tile%name
    ENDIF
    CALL init_lcc(procname, tile, active_cqts, passive_cqts, involved_tiles)

    DEALLOCATE(involved_tiles)

  END SUBROUTINE flcc_init_lcc

#endif
END MODULE mo_flcc_init_lcc
