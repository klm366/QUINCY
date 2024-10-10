! Meri's to-be-deleted development notes:
! No other than name changes relative to alcc? 
! -> copied alcc init_lcc, swapped alcc to hlcc

! TODO: update a check to follow wlcc process (running on leaf not veg tile): (comparison to WLCC: calls dsl4jsb_Use_prcesses_WLCC_ and checks Is_process_calculated(WLCC_) rather than if VEG tile (IF (.NOT. tile%name .EQ. 'veg')) -> alcc approach (copied here) not applicable!?)
! TODO: filenames (ic, bc, filename_prefix) 
!
!NEEDS: Calls init lcc with wood harvesting option (procname)

!-------------------------------------------
!> hlcc (wood harvesting lcc) lcc structure initialisation
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Initialization of the hlcc lcc structure
!>
MODULE mo_hlcc_init_lcc
    #ifndef __NO_JSBACH__
    
      USE mo_exception,           ONLY: message
      USE mo_jsb_control,         ONLY: debug_on
      USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
      USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN
      USE mo_jsb_lcc,             ONLY: init_lcc, count_descendants_for_up_to_2layers, &
        &                               collect_names_of_descendant_leaves_up_to_2layers
      USE mo_jsb_cqt_class,       ONLY: LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, &
        &                               BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE
    
      IMPLICIT NONE
      PRIVATE
    
      PUBLIC :: hlcc_init_lcc
    
      CHARACTER(len=*), PARAMETER :: modname = 'mo_hlcc_init_lcc_class'
      CHARACTER(len=*), PARAMETER :: procname = 'hlcc'
    
      INTEGER, PARAMETER :: nr_of_active_cqts = 1
      INTEGER, PARAMETER :: nr_of_passive_cqts = 4
    
    CONTAINS
    
      ! ====================================================================================================== !
      !
      !> Initialise lcc structure for hlcc (anthropogenic lcc) process
      !
      SUBROUTINE hlcc_init_lcc(tile)
    
        ! -------------------------------------------------------------------------------------------------- !
        CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
        ! -------------------------------------------------------------------------------------------------- !
        CHARACTER(len=*), PARAMETER :: routine = modname//':hlcc_init_lcc'
    
        INTEGER :: active_cqts(nr_of_active_cqts)
        INTEGER :: passive_cqts(nr_of_passive_cqts)
        INTEGER :: nr_of_descendants
    
        CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: involved_tiles(:)
        ! -------------------------------------------------------------------------------------------------- !
    
        ! JN-TODO: reconsider what to define here and how
        IF (.NOT. tile%name .EQ. 'veg') RETURN
    
        IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    
        active_cqts = (/ LIVE_CARBON_CQ_TYPE /)
        passive_cqts = (/ PRODUCT_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, BG_DEAD_C_CQ_TYPE, FLUX_C_CQ_TYPE /)
    
        ! hlcc runs on the leaf tiles
        IF (tile%Has_children()) THEN
          CALL count_descendants_for_up_to_2layers(tile, nr_of_descendants)
          ALLOCATE(involved_tiles(nr_of_descendants))
          CALL collect_names_of_descendant_leaves_up_to_2layers(tile, involved_tiles)
        ELSE
          ALLOCATE(involved_tiles(1))
          involved_tiles(1) = tile%name
        END IF
    
        CALL init_lcc(procname, tile, active_cqts, passive_cqts, involved_tiles)
    
        DEALLOCATE(involved_tiles)
    
      END SUBROUTINE hlcc_init_lcc
    
    #endif
    END MODULE mo_hlcc_init_lcc
    