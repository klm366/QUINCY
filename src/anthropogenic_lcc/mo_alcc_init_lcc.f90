!> Initialisation of the lcc (land cover change) structure for the ALCC (anthropogenic lcc) process
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
!>#### Initialisation of the lcc (land cover change) structure for the ALCC (anthropogenic lcc) process
!>
!> The ALCC process is designed such, that matter is relocated from one tile to another. To assure conservation
!> we make use of the so called lcc structure, which is initialized in routine [[alcc_init_lcc]].
!> All conserved quantity types (CQTs) needed for active or passive matter relocations within the ALCC process
!>  as well as the affected tiles are specified here. Corresponding lists are provided to the lcc framework to
!>  initialize the lcc structure for the ALCC process.
!>
MODULE mo_alcc_init_lcc
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN
  USE mo_jsb_model_class,     ONLY: t_jsb_model, MODEL_JSBACH, MODEL_QUINCY
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_lcc,             ONLY: init_lcc, count_descendants_for_up_to_2layers, &
    &                               collect_names_of_descendant_leaves_up_to_2layers
  USE mo_jsb_cqt_class,       ONLY: LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, BG_DEAD_C_CQ_TYPE, BGLEV_DEAD_C_CQ_TYPE, &
    &                               PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE, &
    &                               IQ_SL_POOL_CQ_TYPE, IQ_SL_FLUX_CQ_TYPE, IQ_1L2D_POOL_CQ_TYPE, IQ_FLUX_CQ_TYPE

  dsl4jsb_Use_processes ALCC_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: alcc_init_lcc

  CHARACTER(len=*), PARAMETER :: modname  = 'mo_alcc_init_lcc_class'
  CHARACTER(len=*), PARAMETER :: procname = 'alcc'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialise the lcc structure of the ALCC process
  !
  SUBROUTINE alcc_init_lcc(tile)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile            !< Current tile
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model           !< Current instance of the model

    INTEGER, ALLOCATABLE :: active_cqts(:)                          !< IDs of CQTs for active matter relocation with ALCC
    INTEGER, ALLOCATABLE :: passive_cqts(:)                         !< IDs of CQTs for passive matter relocation with ALCC
    INTEGER :: nr_of_descendants                                    !< Number of children of the current tile
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: involved_tiles(:) !< List with the names of the tiles that might be affected
                                                                    !< by cover fraction changes due to ALCC

    CHARACTER(len=*), PARAMETER :: routine = modname//':alcc_init_lcc' !< Name of this routine
    ! -------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)

    ! Note: as it is currently implemented ALCC needs to run on the veg tile.
    IF (.NOT. ((TRIM(tile%name) .EQ. 'veg') .AND. tile%Is_process_calculated(ALCC_))) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! The CQTs IDs are different with JSBACH or QUINCY
    IF (model%config%model_scheme == MODEL_JSBACH) THEN
      ALLOCATE(active_cqts(1))
      ALLOCATE(passive_cqts(4))

      active_cqts = (/ LIVE_CARBON_CQ_TYPE /)
      passive_cqts = (/ PRODUCT_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, BG_DEAD_C_CQ_TYPE, BGLEV_DEAD_C_CQ_TYPE, &
        &               FLUX_C_CQ_TYPE /)
    ELSE IF (model%config%model_scheme == MODEL_QUINCY) THEN
      ! Note: due to the complexity of the BGCMs the lcc framework is not used for active relocation with QUINCY.
      !       Active variables are instead explicitly relocated within the ALCC update task!
      ALLOCATE(active_cqts(0))
      ALLOCATE(passive_cqts(4))
      passive_cqts = (/ IQ_SL_POOL_CQ_TYPE, IQ_SL_FLUX_CQ_TYPE, IQ_1L2D_POOL_CQ_TYPE, IQ_FLUX_CQ_TYPE /)
    END IF

    ! The ALCC lcc task runs on all leave tiles of the veg tile:
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

  END SUBROUTINE alcc_init_lcc

#endif
END MODULE mo_alcc_init_lcc
