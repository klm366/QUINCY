!> QUINCY agriculture variables init
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### initialization of agriculture memory variables using, e.g., ic & bc input files
!>
MODULE mo_q_agr_init
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message
  USE mo_jsb_control,         ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: q_agr_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_agr_init'

CONTAINS

  ! ======================================================================================================= !
  !> Run agriculture init
  !>
  SUBROUTINE q_agr_init(tile)
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_grid,              ONLY: Get_grid
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_grid_class,        ONLY: t_jsb_grid
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_process_class,     ONLY: Q_AGR_
    USE mo_q_agr_constants,       ONLY: def_n_fertiliser, def_crop_type, ix_wheat, ix_tecorn, ix_trcorn, ix_trsoy
    USE mo_q_assimi_constants,    ONLY: ic3phot
    dsl4jsb_Use_memory(Q_AGR_)
    dsl4jsb_Use_config(Q_AGR_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid), POINTER       :: grid
    INTEGER                         :: iblk, ic, is         !< loop dimensions
    INTEGER                         :: nproma, nblks        !< dimensions
    INTEGER                         :: lctlib_ps_pathway
    CHARACTER(len=*), PARAMETER :: routine = modname//':q_agr_init'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_AGR_)
    dsl4jsb_Real2D_onDomain :: n_fertiliser_c3
    dsl4jsb_Real2D_onDomain :: n_fertiliser_c4
    dsl4jsb_Real2D_onDomain :: crop_type_index
    dsl4jsb_Def_config(Q_AGR_)
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    IF (.NOT. tile%Is_process_calculated(Q_AGR_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of agr memory (quincy) for tile '// &
      &                          TRIM(tile%name))
    ! ----------------------------------------------------------------------------------------------------- !
    grid   => Get_grid(model%grid_id)
    nblks  =  grid%nblks
    nproma =  grid%nproma
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_AGR_)
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, n_fertiliser_c3)            ! out
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, n_fertiliser_c4)            ! out
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, crop_type_index)            ! out
    dsl4jsb_Get_config(Q_AGR_)
    ! ----------------------------------------------------------------------------------------------------- !

    !> Determine crop type index
    !crop_type_index(:,:) = def_crop_type
    !> set defaults based on photosynthesis type and latitude
    lctlib_ps_pathway = dsl4jsb_Lctlib_param(ps_pathway)
    IF(lctlib_ps_pathway==ic3phot) THEN ! C3 Crops
      DO ic = 1,nproma
        DO iblk = 1,nblks
          IF(grid%lat(ic,iblk) < 30.0_wp .AND. grid%lat(ic,iblk) > -30.0_wp ) THEN
            crop_type_index(ic,iblk) = ix_trsoy ! tropical C3 crop = soybean
          ELSE
            crop_type_index(ic,iblk) = ix_wheat ! temperate C3 crop = wheat
          END IF
        END DO
      END DO
    ELSE ! C4 crops
      DO ic = 1,nproma
        DO iblk = 1,nblks
          IF(grid%lat(ic,iblk) < 30.0_wp .AND. grid%lat(ic,iblk) > -30.0_wp ) THEN
            crop_type_index(ic,iblk) = ix_trcorn ! tropical C4 crop = corn
          ELSE
            crop_type_index(ic,iblk) = ix_tecorn ! temperate C4 crop = corn
          END IF
        END DO
      END DO
    END IF

    !> default fertiliser (if the maps are read per tile, this could be simplified to just read
    ! one variable which points to a different variable in the fertiliser file depending on
    ! crop type rather than allocating two different variables. Leaving this for now, as I don't
    ! know how the reading will look.
    n_fertiliser_c3(:,:) = def_n_fertiliser
    n_fertiliser_c4(:,:) = def_n_fertiliser
    ! ----------------------------------------------------------------------------------------------------- !

  END SUBROUTINE q_agr_init

#endif
END MODULE mo_q_agr_init
