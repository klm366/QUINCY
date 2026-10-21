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

  USE mo_exception,           ONLY: message, finish, message_text

  USE mo_jsb_model_class,     ONLY: t_jsb_model, MODEL_JSBACH, MODEL_QUINCY
  USE mo_jsb_process_class,   ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, RAD_, HD_, ASSIMI_, PHENO_, CARBON_, DISTURB_, &
    &                               CLIMBUF_, FUEL_, FAGE_, ALCC_, FLCC_, WLCC_, NLCC_, PPLCC_, VEG_, SB_, &
    &                               Q_PHENO_, Q_SYL_, Q_AGR_, Q_ASSIMI_
  USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, AGGREGATE_ !, ON_SUBTREE_
  USE mo_jsb_tile,            ONLY: t_jsb_tile
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct, &
    &                               LAND_TYPE, VEG_TYPE, GLACIER_TYPE, LAKE_TYPE !, BARE_TYPE   ! define the type of the tile

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_usecase
  PUBLIC :: validate_single_pft_restart

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_usecases'

CONTAINS

  ! ======================================================================================================= !
  !>initialize the usecase that determines the tiles and what processes are running where
  !>
  SUBROUTINE init_usecase(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase'

    IF (model%config%use_lakes) THEN
      CALL message(TRIM(routine), 'Initializing model usecase "'//TRIM(model%config%usecase)//'" with lakes')
    ELSE
      CALL message(TRIM(routine), 'Initializing model usecase "'//TRIM(model%config%usecase)//'"')
    END IF

    ! select usecase based on namelist option
    SELECT CASE (TRIM(model%config%usecase))

    ! jsbach physics
    CASE ('jsbach_lite')
      IF (model%config%use_tmx) THEN
        CALL init_usecase_lite_tmx(model)
      ELSE
        CALL init_usecase_lite(model)
      END IF

    ! jsbach with PFTs
    CASE ('jsbach_pfts')
      IF (model%config%use_tmx) THEN
        CALL init_usecase_pfts_tmx(model)
      ELSE
        CALL init_usecase_pfts(model)
      END IF

    CASE ('jsbach_forest_age_classes')
      ! jsbach 11 PFTs + forest age-classes
      CALL init_usecase_pfts_with_forest_age_classes(model)

    ! QUINCY 13 PFTs
    CASE ('quincy_13_pfts')
      CALL init_usecase_quincy_13_pfts(model)

    CASE ('quincy_single_pft')
      CALL init_usecase_quincy_single_pft(model)

    ! catch wrong usecase setting
    CASE DEFAULT
      CALL finish(TRIM(routine), 'No such usecase defined.')
    END SELECT

  END SUBROUTINE init_usecase

  ! ======================================================================================================= !
  !>init routine for usecase lite with TMX (explicit coupling)
  !>
  SUBROUTINE init_usecase_lite(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_lite'

    ! overview of tile structure:
    !
    !  box
    !   |
    !   |-------
    !   |      |
    !  land  lakes
    !   |
    !   |-------------
    !   |            |
    !  vegetation  glacier

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_,      HYDRO_,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile
    !$ACC ENTER DATA COPYIN(box_tile)

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

      ! Create sub-tiles
      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile, fract_varname='fract_lake')
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
        !$ACC ENTER DATA COPYIN(lake_tile, land_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), &
          & parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
        !$ACC ENTER DATA COPYIN(land_tile)
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', processes=(/SSE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile, &
            & fract_varname='fract_glac')
          !$ACC ENTER DATA COPYIN(glacier_tile)
        END IF

          veg_tile => t_jsb_tile('veg', 'Veg tile', processes=(/SSE_, PHENO_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')
          !$ACC ENTER DATA COPYIN(veg_tile)

    NULLIFY(box_tile, land_tile, veg_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
    IF (model%config%use_lakes) NULLIFY(lake_tile)

  END SUBROUTINE init_usecase_lite

  SUBROUTINE init_usecase_lite_tmx(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_lite_tmx'

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      ! & processes      =(/A2L_,     L2A_,     RAD_,      HYDRO_,    TURB_    /), &
      ! & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & processes      =(/A2L_,     L2A_,     SEB_,      RAD_,      HYDRO_,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile
    !$ACC ENTER DATA COPYIN(box_tile)

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

      ! Create sub-tiles
      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          ! & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile, fract_varname='fract_lake')
        land_tile => t_jsb_tile('land', 'Land tile', &
          ! & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & processes=(/SSE_/), process_actions=(/ON_LEAFS_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
        !$ACC ENTER DATA COPYIN(lake_tile, land_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          ! & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & processes=(/SSE_/), process_actions=(/ON_LEAFS_/), &
          & lct=t_jsb_lct(LAND_TYPE), &
          & parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
        !$ACC ENTER DATA COPYIN(land_tile)
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile, &
            & fract_varname='fract_glac')

          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes=(/PHENO_/), process_actions=(/ON_TILE_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')
          !$ACC ENTER DATA COPYIN(glacier_tile, veg_tile)
        ELSE
          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes=(/PHENO_/), process_actions=(/ON_TILE_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='equal')
          !$ACC ENTER DATA COPYIN(veg_tile)
        END IF

    NULLIFY(box_tile, land_tile, veg_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
    IF (model%config%use_lakes) NULLIFY(lake_tile)

  END SUBROUTINE init_usecase_lite_tmx

  ! ======================================================================================================= !
  !>init routine for usecase pfts
  !>
  SUBROUTINE init_usecase_pfts(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile, pft_tile

    CHARACTER(len=5)  :: pft_shortnames(model%config%npft) !< Short names of the PFTs
    CHARACTER(len=30) :: pft_longnames(model%config%npft)  !< Long names of the PFTs
    INTEGER :: pft_ids(model%config%npft)                  !< IDs of PFTs in lctlib (shifted)
    INTEGER :: i                                           !< Looping index

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_pfts'

    ! Define PFT names and indices depending on the number of PFTs
    CALL set_jsbach_pfts(model%config%npft, pft_shortnames(:), pft_longnames(:), pft_ids(:))

    ! overview of tile structure:
    !
    !  box
    !   |
    !   |-------
    !   |      |
    !  land  lakes
    !   |
    !   |-------------
    !   |            |
    !  vegetation  glacier
    !   |
    !   |-----------------
    !   |      |      |
    !  pft01  pft02  pftXX

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_  ,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile
    !$ACC ENTER DATA COPYIN(box_tile)

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

    IF (model%processes(PPLCC_)%p%config%active) THEN
      CALL box_tile%Set_process_action(PPLCC_, ON_LEAFS_)
    END IF

      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          & processes=(/SEB_, HYDRO_/), process_actions=(/ON_TILE_, ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile)
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
        !$ACC ENTER DATA COPYIN(lake_tile, land_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
        !$ACC ENTER DATA COPYIN(land_tile)
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
            & processes = (/HYDRO_, SSE_/), &
            & process_actions = (/ON_TILE_, ON_TILE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile)
            ! & lct=t_jsb_lct(GLACIER_TYPE, lib_id=1), parent=land_tile)
          !$ACC ENTER DATA COPYIN(glacier_tile)
        END IF

        ! HYDRO_ is active on veg tile only to save runtime
        veg_tile => t_jsb_tile('veg', 'Veg tile', &
          & processes       = (/HYDRO_,   SSE_,     ASSIMI_,   PHENO_/),    &
          & process_actions = (/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
          & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
          & fract_varname='fract_veg')   ! assumes that this usecase always runs with glacier tile (see also init_usecase_lite)
        !$ACC ENTER DATA COPYIN(veg_tile)

        IF (model%processes(ALCC_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(ALCC_, ON_TILE_, AGGREGATE_)
        END IF

        ! Note: FLCC and WLCC are required for the fage usecase. Here, they can be activated for testing,
        !   however, for the pft usecase it currently makes more sense using DISTURB without FLCC and WLCC
        IF (model%processes(FLCC_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(FLCC_, ON_TILE_, AGGREGATE_)
        END IF
        IF (model%processes(WLCC_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(WLCC_, ON_TILE_, AGGREGATE_)
        END IF

        IF (model%processes(CLIMBUF_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(CLIMBUF_, ON_TILE_, AGGREGATE_)
        END IF
        IF (model%processes(NLCC_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(NLCC_, ON_TILE_, AGGREGATE_)
        END IF
        !
        IF (model%processes(CARBON_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(CARBON_, ON_LEAFS_, AGGREGATE_)
        END IF
        IF (model%processes(FUEL_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(FUEL_, ON_TILE_, AGGREGATE_)
        END IF
        IF (model%processes(DISTURB_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(DISTURB_, ON_LEAFS_, AGGREGATE_)
        END IF

          DO i=1,model%config%npft
            pft_tile => t_jsb_tile(pft_shortnames(i), TRIM(pft_longnames(i)), &
              & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=veg_tile)
            !$ACC ENTER DATA COPYIN(pft_tile)
            NULLIFY(pft_tile)
          END DO

    NULLIFY(box_tile, land_tile, veg_tile, pft_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
    IF (model%config%use_lakes) NULLIFY(lake_tile)

  END SUBROUTINE init_usecase_pfts

  ! ======================================================================================================= !
  !>init routine for usecase pfts with TMX (explicit coupling)
  !>
  SUBROUTINE init_usecase_pfts_tmx(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile, pft_tile

    CHARACTER(len=5)  :: pft_shortnames(model%config%npft) !< Short names of the PFTs
    CHARACTER(len=30) :: pft_longnames(model%config%npft)  !< Long names of the PFTs
    INTEGER :: pft_ids(model%config%npft)                  !< IDs of PFTs in lctlib (shifted)
    INTEGER :: i                                           !< Looping index

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_pfts_tmx'

    ! Define PFT names and indices depending on the number of PFTs
    CALL set_jsbach_pfts(model%config%npft, pft_shortnames(:), pft_longnames(:), pft_ids(:))

    ! overview of tile structure:
    !
    !  box
    !   |
    !   |-------
    !   |      |
    !  land  lakes
    !   |
    !   |-------------
    !   |            |
    !  vegetation  glacier
    !   |
    !   |-----------------
    !   |      |      |
    !  pft01  pft02  pftXX

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_     /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile
    !$ACC ENTER DATA COPYIN(box_tile)

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

    IF (model%processes(PPLCC_)%p%config%active) THEN
      CALL box_tile%Set_process_action(PPLCC_, ON_LEAFS_)
    END IF

      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          & processes=(/SEB_, TURB_, HYDRO_/), process_actions=(/ON_TILE_, ON_TILE_, ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile)
        land_tile => t_jsb_tile('land', 'Land tile', &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
        !$ACC ENTER DATA COPYIN(lake_tile, land_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
        !$ACC ENTER DATA COPYIN(land_tile)
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
            & processes       = (/SEB_,     TURB_,    HYDRO_,   SSE_/), &
            & process_actions = (/ON_TILE_, ON_TILE_, ON_TILE_, ON_TILE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile)
          !$ACC ENTER DATA COPYIN(glacier_tile)
        END IF

        veg_tile => t_jsb_tile('veg', 'Veg tile', &
          & processes       = (/SEB_,     TURB_,    HYDRO_,   SSE_,     ASSIMI_,   PHENO_/),    &
          & process_actions = (/ON_TILE_, ON_TILE_, ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
          & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
          & fract_varname='fract_veg')
        !$ACC ENTER DATA COPYIN(veg_tile)

        IF (model%processes(ALCC_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(ALCC_, ON_TILE_, AGGREGATE_)
        END IF

        ! Note: FLCC and WLCC are required for the fage usecase. Here, they can be activated for testing,
        !   however, for the pft usecase it currently makes more sense using DISTURB without FLCC and WLCC
        IF (model%processes(FLCC_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(FLCC_, ON_TILE_, AGGREGATE_)
        END IF
        IF (model%processes(WLCC_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(WLCC_, ON_TILE_, AGGREGATE_)
        END IF

        IF (model%processes(CLIMBUF_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(CLIMBUF_, ON_TILE_, AGGREGATE_)
        END IF
        IF (model%processes(NLCC_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(NLCC_, ON_TILE_, AGGREGATE_)
        END IF
        !
        IF (model%processes(CARBON_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(CARBON_, ON_LEAFS_, AGGREGATE_)
        END IF
        IF (model%processes(FUEL_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(FUEL_, ON_TILE_, AGGREGATE_)
        END IF
        IF (model%processes(DISTURB_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(DISTURB_, ON_LEAFS_, AGGREGATE_)
        END IF

          DO i=1,model%config%npft
            pft_tile => t_jsb_tile(pft_shortnames(i), TRIM(pft_longnames(i)), &
              & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=veg_tile)
            !$ACC ENTER DATA COPYIN(pft_tile)
            NULLIFY(pft_tile)
          END DO

    NULLIFY(box_tile, land_tile, veg_tile, pft_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
    IF (model%config%use_lakes) NULLIFY(lake_tile)

  END SUBROUTINE init_usecase_pfts_tmx

  ! ======================================================================================================= !
  !>init routine for usecase pfts with forest age classes
  !>
  SUBROUTINE init_usecase_pfts_with_forest_age_classes(model)

    dsl4jsb_Use_config(FAGE_)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile, pft_tile, fage_tile

    dsl4jsb_Def_config(FAGE_)

    INTEGER :: i, j
    INTEGER :: nacs
    INTEGER :: pft_ids(model%config%npft)                   !< IDs of PFTs in lctlib (shifted)
    CHARACTER(len=30) :: pft_longnames(model%config%npft)   !< Long names of the PFTs
    CHARACTER(len=30) :: pft_shortdummy(model%config%npft)  !< PFT short names (not used)
    CHARACTER(len=5)  :: pft_shortname                      !< Short name of the PFT
    CHARACTER(len=10) :: age_class_name                     !< Name of the age class
    CHARACTER(len=10) :: init_ac_fracts_scheme              !< Scheme for initializing age class fractions

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_pfts_with_forest_age_classes'

    ! Define PFT names and indices depending on the number of PFTs
    CALL set_jsbach_pfts(model%config%npft, pft_shortdummy(:), pft_longnames(:), pft_ids(:))

    ! overview of tile structure:
    !
    !  box
    !   |
    !   |-------
    !   |      |
    !  land  lakes
    !   |
    !   |-------------
    !   |            |
    !  vegetation  glacier
    !   |
    !   |---------------...------------------------...---
    !   |                           |       |           |
    !  pft01                       pft04   pft05        pft11
    !   |                           |
    !   |---------------...         |---------------...
    !   |           |               |           |
    !  pft01_fc01  pft01_fc02 ...  pft04_fc01  pft04_fc02 ...

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_  ,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile
    !$ACC ENTER DATA COPYIN(box_tile)

    dsl4jsb_Get_config(FAGE_)

    ! FAGE_ is the process that organises ageing in the age classes.
    ! It does not make sense to run age classes without FAGE_ and vice versa!
    IF (.NOT. model%processes(FAGE_)%p%config%active) THEN
      CALL finish(TRIM(routine), "It does not make sense to use age classes without FAGE_ and vice versa!")
    END IF
    IF (.NOT. model%processes(CARBON_)%p%config%active) THEN
      CALL finish(TRIM(routine), "Current forest age class implementation assumes CARBON_ process to be active!")
    END IF
    IF (.NOT. model%processes(PHENO_)%p%config%active) THEN
      CALL finish(TRIM(routine), "Current forest age class implementation assumes PHENO_ process to be active!")
    END IF

    nacs = dsl4jsb_Config(FAGE_)%nacs
    init_ac_fracts_scheme = dsl4jsb_Config(FAGE_)%init_ac_fracts_scheme

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

    IF (model%processes(PPLCC_)%p%config%active) THEN
      CALL box_tile%Set_process_action(PPLCC_, ON_LEAFS_)
    END IF

      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          & processes=(/SEB_, HYDRO_/), process_actions=(/ON_TILE_, ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile)
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
        !$ACC ENTER DATA COPYIN(lake_tile, land_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
        !$ACC ENTER DATA COPYIN(land_tile)
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
            & processes = (/HYDRO_, SSE_/), &
            & process_actions = (/ON_TILE_, ON_TILE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile)
            ! & lct=t_jsb_lct(GLACIER_TYPE, lib_id=1), parent=land_tile)
          !$ACC ENTER DATA COPYIN(glacier_tile)
        END IF

        !>veg_tile
        !!
        ! default jsb4: HYDRO_ & SSE_ are active at veg tile only to save runtime
          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes       = (/HYDRO_,   SSE_,     ASSIMI_,   PHENO_/),    &
            & process_actions = (/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')
          !$ACC ENTER DATA COPYIN(veg_tile)

          IF (model%processes(ALCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(ALCC_, ON_TILE_, AGGREGATE_)
          END IF
          !
          IF (model%processes(CLIMBUF_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(CLIMBUF_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(NLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(NLCC_, ON_TILE_, AGGREGATE_)
          END IF

          IF (model%processes(CARBON_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(CARBON_, ON_LEAFS_, AGGREGATE_)
          END IF
          IF (model%processes(FUEL_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(FUEL_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(DISTURB_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(DISTURB_, ON_LEAFS_, AGGREGATE_)
          END IF

            DO i=1,model%config%npft

              WRITE(pft_shortname, "(A3,I2.2)") "pft", i
              pft_tile => t_jsb_tile(pft_shortname, TRIM(pft_longnames(i)), &
                & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=veg_tile)
              !$ACC ENTER DATA COPYIN(pft_tile)

                IF (model%processes(FLCC_)%p%config%active) THEN
                  CALL pft_tile%Set_process_action(FLCC_, ON_TILE_, AGGREGATE_)
                END IF
                IF (model%processes(WLCC_)%p%config%active) THEN
                  CALL pft_tile%Set_process_action(WLCC_, ON_TILE_, AGGREGATE_)
                END IF

              IF (model%lctlib(pft_ids(i))%ForestFlag) THEN

                CALL pft_tile%Set_process_action(FAGE_, ON_TILE_, AGGREGATE_)

                DO j=1,nacs
                  WRITE(age_class_name, "(A5,A3,I2.2)") pft_shortname, "_fc", j
                  fage_tile => t_jsb_tile(age_class_name, TRIM(age_class_name // ' tile'), &
                    & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=pft_tile, &
                    & fract_varname=TRIM(init_ac_fracts_scheme))
                  !$ACC ENTER DATA COPYIN(fage_tile)
                END DO
              END IF
            END DO

    ! ! Next 3 lines only to avoid compiler warnings about pointers not being dereferenced
    ! IF (ASSOCIATED(glacier_tile)) CONTINUE
    ! IF (ASSOCIATED(lake_tile)) CONTINUE
    ! IF (ASSOCIATED(pft_tile)) CONTINUE

    NULLIFY(box_tile, land_tile, veg_tile, pft_tile, fage_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
    IF (model%config%use_lakes) NULLIFY(lake_tile)

  END SUBROUTINE init_usecase_pfts_with_forest_age_classes

  !>
  !> Set jsbach PFT names and indices
  !>
  SUBROUTINE set_jsbach_pfts(npft, pft_shortnames, pft_longnames, pft_ids)

    INTEGER,           INTENT(in)  :: npft                  !< Number of PFTs
    CHARACTER(len=5),  INTENT(out) :: pft_shortnames(npft)  !< PFT names (e.g. in output)
    CHARACTER(len=30), INTENT(out) :: pft_longnames(npft)   !< PFT description
    INTEGER,           INTENT(out) :: pft_ids(npft)         !< IDs of the PFTs

    CHARACTER(len=*), PARAMETER :: routine = modname//':set_jsbach_pfts'

    ! Define PFT names and indices based on the number of PFTs
    IF (npft == 11) THEN
      pft_shortnames(:) = [character(len=5) ::               &
        & 'pft01', 'pft02', 'pft03', 'pft04', 'pft05', 'pft06', &
        & 'pft07', 'pft08', 'pft09', 'pft10', 'pft11' ]

      pft_longnames(:) = [character(len=30) :: &
        & 'Tropical broadleaf evergreen', &
        & 'Tropical broadleaf deciduous', &
        & 'Extra-tropical evergreen',     &
        & 'Extra-tropical deciduous',     &
        & 'Raingreen shrubs',             &
        & 'Deciduous shrubs',             &
        & 'C3 grass',                     &
        & 'C4 grass',                     &
        & 'C3 pasture',                   &
        & 'C4 pasture',                   &
        & 'C3 crops'                     ]

      ! Index of PFT for lctlib columns
      ! Note: The IDs defined here are each one less than in the lctlib because glacier (first
      !       column in lctlib) is not a PFT and ignored here.
      pft_ids(:) = (/ 1,2,3,4,9,10,11,12,14,15,19 /)

    ELSE IF (npft == 12) THEN
      pft_shortnames(:) = [character(len=5) ::               &
        & 'pft01', 'pft02', 'pft03', 'pft04', 'pft05', 'pft06', &
        & 'pft07', 'pft08', 'pft09', 'pft10', 'pft11', 'pft12' ]

      pft_longnames(:) = [character(len=30) :: &
        & 'Tropical broadleaf evergreen', &
        & 'Tropical broadleaf deciduous', &
        & 'Extra-tropical evergreen',     &
        & 'Extra-tropical deciduous',     &
        & 'Raingreen shrubs',             &
        & 'Deciduous shrubs',             &
        & 'C3 grass',                     &
        & 'C4 grass',                     &
        & 'C3 pasture',                   &
        & 'C4 pasture',                   &
        & 'C3 crops',                     &
        & 'C4 crops'                     ]

      ! Index of PFT for lctlib columns
      ! Note: The IDs defined here are each one less than in the lctlib because glacier (first
      !       column in lctlib) is not a PFT and ignored here.
      pft_ids(:) = (/ 1,2,3,4,9,10,11,12,14,15,19,20 /)

    ELSE
      WRITE (message_text,*) 'Unsupported number of PFTs: ', npft
      CALL finish(routine, message_text)
    END IF

    END SUBROUTINE set_jsbach_pfts

  ! ======================================================================================================= !
  !>init routine for usecase quincy_13_pfts
  !>
  !>  works only in combination with 'model%config%model_scheme = MODEL_QUINCY'
  !>
  !>  may be used for both: standalone and coupled simulation setup
  !>
  !>  the data/lctlib_quincy_nlct14.def is used for this usecase (hardcoded in: mo_jsb_base:jsbach_setup_tiles)
  !>
  SUBROUTINE init_usecase_quincy_13_pfts(model, selected_pft)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model                    !< model instance
    INTEGER, OPTIONAL, INTENT(in) :: selected_pft !< Original PFT index for one-child construction [-]

    CLASS(t_jsb_tile), POINTER    :: box_tile, land_tile, &
      &                              veg_tile, pft_tile, &
      &                              lake_tile, glacier_tile              !< tiles used with this usecase
    INTEGER                       :: i                                    !< for looping
    CHARACTER(len=5) :: fraction_name !< Empty: original PFT field; equal: full parent area
    ! number of pft
    INTEGER, PARAMETER            :: npft = 13
    ! pft short names
    CHARACTER(len=5), PARAMETER   :: pft_shortnames(npft) = [character(len=5) :: &
      &                              'pft01', 'pft02', 'pft03', 'pft04', &
      &                              'pft05', 'pft06', 'pft07', 'pft08', &
      &                              'pft09', 'pft10', 'pft11', 'pft12', 'pft13'  ]
    ! pft long names
    CHARACTER(len=30), PARAMETER  :: pft_longnames(npft) = [character(len=30) ::    &
      &                              'BEM:  moist broadleaved evergreen',           &
      &                              'BED:  dry broadleaved evergreen',             &
      &                              'BDR:  rain green broadleaved deciduous',      &
      &                              'BDS:  summer green broadleaved deciduous',    &
      &                              'NE :  needle-leaved evergreen',               &
      &                              'NS :  summer green needle-leaved',            &
      &                              'TeH:  C3 grass',                              &
      &                              'TrH:  C4 grass',                              &
      &                              'TeP:  C3 pasture',                            &
      &                              'TrP:  C4 pasture',                            &
      &                              'TeC:  C3 crop',                               &
      &                              'TrC:  C4 crop',                               &
      &                              'BSO:  bare soil' ]
    ! index of PFT (columns) in the lctlib file
    INTEGER, PARAMETER            :: pft_ids(npft) = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13 /)
    CHARACTER(len=*), PARAMETER   :: routine = modname//':init_usecase_quincy_13_pfts'
    !--------------------------------------------------------------------------------------------------------
    ! overview of tile structure:
    !
    !  box
    !   |
    !   |-------
    !   |      |
    !  land  lakes
    !   |
    !   |-------------
    !   |            |
    !  vegetation  glacier
    !   |
    !   |-----------------
    !   |      |      |        |
    !  pft01  pft02  pftXX  baresoil
    !--------------------------------------------------------------------------------------------------------
    ! QUINCY model must be enabled for this usecase
    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_QUINCY)
      CALL message(TRIM(routine), 'Using the model_scheme MODEL_QUINCY')
    CASE DEFAULT
      CALL finish(TRIM(routine), 'The quincy_13_pfts model usecase does work only with the model_scheme MODEL_QUINCY.')
    END SELECT

    !------------------------------------------------------------------------------------------------------ !
    ! box tile
    !
    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_     , TURB_     /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_ /), &
      & model_id=model%id, fract_varname='notsea')  ! fraction: all that "is not sea"
    model%top => box_tile
    !$ACC ENTER DATA COPYIN(box_tile)

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

    ! IF PPLCC is active, then PPLCC should run for Quincy
    ! Note: it is not only used for LCC processes but also to provide tree area to the Q_SYL_ process
    IF (model%processes(PPLCC_)%p%config%active) THEN
      CALL box_tile%Set_process_action(PPLCC_, ON_LEAFS_)
    END IF

    !------------------------------------------------------------------------------------------------------ !
    ! land and lake tile
    !
    IF (model%config%use_lakes) THEN
      lake_tile => t_jsb_tile('lake', 'Lake tile', &
        & processes=(/SEB_, HYDRO_/), process_actions=(/ON_TILE_, ON_TILE_/), &
        & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile)
      land_tile => t_jsb_tile('land', 'Land tile', &
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
      !$ACC ENTER DATA COPYIN(lake_tile)
    ELSE
      land_tile => t_jsb_tile('land', 'Land tile', &
        ! add processes here only to make model init and jsb4_forcing working
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        ! none of the tasks of these processes is running via task queue
        & lct=t_jsb_lct(LAND_TYPE), parent=box_tile, &
        & fract_varname='equal')   ! Fraction of land tile = 1
                                   ! "equal" specifically means, no fraction input file is read and fraction is identical with above tile in hierarchy
    END IF
    !$ACC ENTER DATA COPYIN(land_tile)

    !------------------------------------------------------------------------------------------------------ !
    ! glacier tile
    !
    IF (model%config%use_glacier) THEN
      glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
        & processes = (/HYDRO_, SSE_/), &
        & process_actions = (/ON_TILE_, ON_TILE_/), &
        & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile)
      !$ACC ENTER DATA COPYIN(glacier_tile)
    END IF

    !------------------------------------------------------------------------------------------------------ !
    ! vegetation tile with PFT tiles as leafs
    !   fraction is "equal" because it is the only tile at this level (same as land tile)
    !
    ! on veg tile:  SSE_ and HYDRO_
    IF (model%config%use_glacier) THEN
      veg_tile => t_jsb_tile('veg', 'Veg tile', &
        & processes       = (/Q_ASSIMI_, Q_PHENO_,  VEG_,      SB_,       SSE_,      HYDRO_   /), &
        & process_actions = (/ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_TILE_,  ON_TILE_ /), &
        & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
        & fract_varname='fract_veg')  ! fraction of land tile shared with glacier tile
    ELSE
      veg_tile => t_jsb_tile('veg', 'Veg tile', &
        & processes       = (/Q_ASSIMI_, Q_PHENO_,  VEG_,      SB_,       SSE_,      HYDRO_   /), &
        & process_actions = (/ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_TILE_,  ON_TILE_ /), &
        & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
        & fract_varname='equal')    ! fraction = 1.0
    END IF
    !$ACC ENTER DATA COPYIN(veg_tile)

    IF (model%processes(ALCC_)%p%config%active) THEN
      CALL veg_tile%Set_process_action(ALCC_, ON_TILE_, AGGREGATE_)
    END IF

    ! create PFT tiles 1:npft
    DO i = 1,npft
      fraction_name = ''
      IF (PRESENT(selected_pft)) THEN
        IF (i /= selected_pft) CYCLE
        ! The sole child covers vegetation, not the whole grid cell.
        fraction_name = 'equal'
      END IF
      ! via constructor of t_jsb_tile
      pft_tile => t_jsb_tile(pft_shortnames(i), &
        &                    TRIM(pft_longnames(i)), &
        &                    lct = t_jsb_lct(VEG_TYPE, lib_id = pft_ids(i)), &
        &                    parent = veg_tile, fract_varname = fraction_name)

      IF (model%lctlib(pft_ids(i))%ForestFlag) THEN
        IF (model%processes(Q_SYL_)%p%config%active) THEN
          CALL pft_tile%Set_process_action(Q_SYL_, ON_TILE_, AGGREGATE_)
        END IF
      END IF

      IF (model%lctlib(pft_ids(i))%CropFlag) THEN
        IF (model%processes(Q_AGR_)%p%config%active) THEN
          CALL pft_tile%Set_process_action(Q_AGR_, ON_TILE_, AGGREGATE_)
        END IF
      END IF
      !$ACC ENTER DATA COPYIN(pft_tile)
    END DO

    ! catch compiler warnings and cleanup pointers
    !
    ! ! Next line only to avoid compiler warnings about pointers not being dereferenced
    ! IF (ASSOCIATED(pft_tile)) CONTINUE
    NULLIFY(box_tile, land_tile, veg_tile, pft_tile)
    IF (model%config%use_lakes)   NULLIFY(lake_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
  END SUBROUTINE init_usecase_quincy_13_pfts

  !> Construct one dominant vegetated PFT before trait-dependent initialization.
  !! Read_1d broadcasts the tiny global field through host I/O, so selection is
  !! independent of local cell ownership. The input ranks candidates only; the
  !! winner inherits all parent vegetation area via the existing equal rule.
  SUBROUTINE init_usecase_quincy_single_pft(model)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    USE mo_kind, ONLY: lp
    USE mo_jsb_grid, ONLY: Get_grid
    USE mo_jsb_grid_class, ONLY: t_jsb_grid
    USE mo_jsb_io_netcdf, ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_quincy_model_config, ONLY: QSOIL
    USE mo_sb_config_class, ONLY: t_sb_config

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model !< Model with configured grid and QUINCY library
    TYPE(t_jsb_grid), POINTER :: grid !< Global cell count, independent of local ownership
    TYPE(t_input_file) :: input_file !< Existing initial fractions
    REAL(lp), POINTER :: fraction(:) !< One-cell input fraction [-]
    REAL(lp) :: largest_fraction !< Largest positive vegetated share in the input [-]
    INTEGER :: i, winner !< Original PFT indices; zero means no candidate [-]
    CHARACTER(len=11) :: field_name !< Original fraction field, fract_pftNN
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_quincy_single_pft'

    IF (model%config%model_scheme /= MODEL_QUINCY) &
      & CALL finish(routine, 'quincy_single_pft requires the QUINCY model scheme.')
    grid => Get_grid(model%grid_id)
    IF (.NOT. model%config%is_single_site .OR. grid%ntotal_g /= 1) &
      & CALL finish(routine, 'quincy_single_pft requires a single-cell site-level grid.')
    IF (model%Do_fractions_change()) &
      & CALL finish(routine, 'quincy_single_pft requires fixed land cover; disable land-cover change.')
    ! External QSOIL forcing still assumes a complete, contiguous PFT inventory
    ! (including pft01) in its cache. Internal SB replay does not use that cache.
    IF (model%config%qmodel_id == QSOIL) &
      & CALL finish(routine, 'quincy_single_pft does not yet support file-driven QSOIL; use coupled SB replay.')
    IF (model%Is_process_enabled(SB_)) THEN
      SELECT TYPE (config => model%processes(SB_)%p%config)
      TYPE IS (t_sb_config)
        ! Reject rather than silently override the existing execution selector.
        ! All common forcing, scheduling and restart bookkeeping stays in use.
        IF (config%flag_sb_model_spinup_accelerator .AND. &
          & config%sb_model_spinup_accelerator_replay_mode /= 'per_pft') &
          & CALL finish(routine, &
          & 'quincy_single_pft requires sb_model_spinup_accelerator_replay_mode="per_pft" for SB replay.')
      CLASS DEFAULT
        CALL finish(routine, 'Invalid SB configuration for quincy_single_pft.')
      END SELECT
    END IF

    winner = 0
    largest_fraction = 0._lp
    input_file = jsb_netcdf_open_input(TRIM(model%config%fract_filename))
    DO i = 1,13
      WRITE(field_name,'(A,I2.2)') 'fract_pft',i
      IF (.NOT. input_file%Has_var(field_name)) &
        & CALL finish(routine, 'Missing initial fraction: '//field_name)
      fraction => input_file%Read_1d(field_name)
      IF (SIZE(fraction) /= 1) CALL finish(routine, 'Expected one cell in '//field_name)
      IF (.NOT. ieee_is_finite(fraction(1))) CALL finish(routine, 'Nonfinite fraction in '//field_name)
      IF (fraction(1) < 0._lp .OR. fraction(1) > 1._lp) &
        & CALL finish(routine, 'Fraction outside [0,1] in '//field_name)
      ! Original library identity is retained. Strict > resolves exact ties
      ! to the lowest index; bare soil is not a vegetated candidate.
      IF (.NOT. model%lctlib(i)%BareSoilFlag .AND. fraction(1) > largest_fraction) THEN
        winner = i
        largest_fraction = fraction(1)
      END IF
      DEALLOCATE(fraction)
    END DO
    CALL input_file%Close()
    IF (winner == 0) CALL finish(routine, 'No positive vegetated PFT fraction in initial input.')
    ! This expected identity is recomputed from immutable initial input, not
    ! overwritten by restart loading. The independently stored marker must match.
    model%single_pft_index = winner

    WRITE(message_text,'(A,I2.2,A,ES16.8,A)') 'Selected pft',winner, &
      & '; original fraction=',largest_fraction,'; assigned fraction=full parent vegetation area.'
    CALL message(routine, TRIM(message_text))
    CALL message(routine, 'Selection source: '//TRIM(model%config%fract_filename))
    CALL init_usecase_quincy_13_pfts(model, selected_pft=winner)
  END SUBROUTINE init_usecase_quincy_single_pft

  !> Check loaded identity before process-state overrides, BGCM synchronization or replay.
  !! Fixed fractions are not restart fields. A mandatory box-level marker stops
  !! both accidental 13-PFT imports and changed dominant-PFT input on continuation.
  !! Only owned cells are tested; rank padding is not physical state.
  SUBROUTINE validate_single_pft_restart(model)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    USE mo_kind, ONLY: lp
    USE mo_jsb_grid, ONLY: Get_grid
    USE mo_jsb_grid_class, ONLY: t_jsb_grid
    TYPE(t_jsb_model), INTENT(in) :: model !< Model after restart file loading
    TYPE(t_jsb_grid), POINTER :: grid !< Local owned-cell bounds
    INTEGER :: iblk, ic !< Block and cell indices [-]
    REAL(lp) :: stored_id !< Restart marker in original lctlib indexing [-]
    CHARACTER(len=*), PARAMETER :: routine = modname//':validate_single_pft_restart'

    IF (model%config%usecase /= 'quincy_single_pft') RETURN
    IF (.NOT. ASSOCIATED(model%single_pft_restart_id)) &
      & CALL finish(routine, 'Missing single-PFT identity registration.')
    grid => Get_grid(model%grid_id)
    DO iblk = grid%Get_blk_start(), grid%Get_blk_end()
      DO ic = grid%Get_col_start(iblk), grid%Get_col_end(iblk)
        stored_id = model%single_pft_restart_id(ic,iblk)
        IF (.NOT. ieee_is_finite(stored_id)) CALL finish(routine, 'Nonfinite single-PFT restart identity.')
        ! Indices 1:13 are exactly representable; fractional or zero markers
        ! must not be rounded into a valid identity. No tolerance is appropriate.
        IF (stored_id /= REAL(model%single_pft_index,lp)) THEN
          WRITE(message_text,*) 'Single-PFT restart identity mismatch: stored/expected=', &
            & stored_id,model%single_pft_index,'; do not reuse 13-PFT restarts or change selected vegetation.'
          CALL finish(routine, TRIM(message_text))
        END IF
      END DO
    END DO
    CALL message(routine, 'Single-PFT restart identity matches initial selection.')
  END SUBROUTINE validate_single_pft_restart

#endif
END MODULE mo_jsb_model_usecases
