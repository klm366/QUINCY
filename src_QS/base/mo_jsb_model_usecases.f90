!> Contains methods for initialization of JSBACH usecases.
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
!! This module defines various model usecases that can be enabled individually via namelist.
!! A model usecase defines:
!!    * the tiles that are enabled in the model and their nested structure
!!    * the processes that are active on each tile or its leaves
!!    * the fraction of each tile
!!    * the land-cover-type of each tile

MODULE mo_jsb_model_usecases
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_process_class,   ONLY: A2L_, L2A_, Q_ASSIMI_, VEG_, Q_PHENO_, SPQ_, SB_, DFIRE_, PHYD_, &
    &                               RAD_, TURB_, SEB_, SSE_, HYDRO_, Q_SYL_, Q_AGR_, &
    &                               ON_LEAFS_, ON_TILE_, AGGREGATE_
  USE mo_jsb_tile,            ONLY: t_jsb_tile
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct, &
    &                               LAND_TYPE, VEG_TYPE, GLACIER_TYPE, LAKE_TYPE !, BARE_TYPE   ! define the type of the tile
  USE mo_veg_config_class,    ONLY: t_veg_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_usecase

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_usecases'

CONTAINS

  SUBROUTINE init_usecase(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase'

      CALL init_usecase_lite(model)

  END SUBROUTINE init_usecase

  SUBROUTINE init_usecase_lite(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER   :: box_tile
    INTEGER                      :: pft_id
    CHARACTER(len=*), PARAMETER  :: routine = TRIM(modname)//':init_usecase_lite'

    ! Create just one tile
    box_tile => t_jsb_tile('veg', 'quincy_standalone_pft', &
      & processes       = (/A2L_,     L2A_,     Q_ASSIMI_, VEG_,      Q_PHENO_,  SPQ_,      SB_      , &
      &   DFIRE_   , PHYD_    , RAD_     , TURB_    , SEB_     , SSE_     , HYDRO_ /), &
      & process_actions = (/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, &
      &   ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_ /), &
      & lct=t_jsb_lct(VEG_TYPE), &
      & model_id=model%id, grid_id=model%grid_id, fract_varname='notsea')

      SELECT TYPE (config => model%processes(VEG_)%p%config)
      TYPE IS (t_veg_config)
        pft_id = config%pft_id
      END SELECT

      IF (model%lctlib(pft_id)%ForestFlag) THEN
        IF (model%processes(Q_SYL_)%p%config%active) THEN
          CALL box_tile%Set_process_action(Q_SYL_, ON_TILE_, AGGREGATE_)
        END IF
      END IF

      IF (model%lctlib(pft_id)%CropFlag) THEN
        IF (model%processes(Q_AGR_)%p%config%active) THEN
          CALL box_tile%Set_process_action(Q_AGR_, ON_TILE_, AGGREGATE_)
        END IF
      END IF


    model%top => box_tile

    NULLIFY(box_tile)

  END SUBROUTINE init_usecase_lite

#endif
END MODULE mo_jsb_model_usecases
