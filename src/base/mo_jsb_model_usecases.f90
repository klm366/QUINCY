!> Contains methods for initialization of JSBACH usecases.
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
!! This module defines various model usecases that can be enabled individually via namelist.
!! A model usecase defines:
!!    * the tiles that are enabled in the model and their nested structure
!!    * the processes that are active on each tile or its leaves
!!    * the fraction of each tile
!!    * the land-cover-type of each tile

MODULE mo_jsb_model_usecases
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message, finish

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_process_class,   ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, RAD_, HD_, ASSIMI_, PHENO_, CARBON_, DISTURB_, &
                                    FUEL_, FAGE_, ALCC_, FLCC_, WLCC_, NLCC_, VEG_, PPLCC_, MIX_
  USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, AGGREGATE_ !, ON_SUBTREE_
  USE mo_jsb_tile,            ONLY: t_jsb_tile
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct, &
    &                               LAND_TYPE, VEG_TYPE, GLACIER_TYPE, LAKE_TYPE, & !, BARE_TYPE
    &                               UPLAND_TYPE, LOWLAND_TYPE, FLOODPLAIN_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_usecase

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_usecases'

CONTAINS

  SUBROUTINE init_usecase(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase'

    IF (model%config%use_lakes) THEN
      CALL message(TRIM(routine), 'Initializing JSBACH usecase "'//TRIM(model%config%usecase)//'" with lakes')
    ELSE
      CALL message(TRIM(routine), 'Initializing JSBACH usecase "'//TRIM(model%config%usecase)//'"')
    END IF


    SELECT CASE (TRIM(model%config%usecase))

    CASE ('jsbach_lite')
      IF (model%config%use_tmx) THEN
        CALL init_usecase_lite_new(model)
      ELSE
        CALL init_usecase_lite(model)
      END IF

    CASE ('jsbach_pfts')
      CALL init_usecase_pfts(model)

    CASE ('jsbach_hydrotiles')
      CALL init_usecase_hydrotiles(model)

    CASE ('jsbach_forest_age_classes')
      ! jsbach4 11 PFTs + forest age-classes
      CALL init_usecase_pfts_with_forest_age_classes(model)

    CASE DEFAULT
      CALL finish(TRIM(routine), 'JSBACH usecase not defined.')

    END SELECT

  END SUBROUTINE init_usecase

  SUBROUTINE init_usecase_lite(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_lite'

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

  SUBROUTINE init_usecase_lite_new(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_lite_new'

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

  END SUBROUTINE init_usecase_lite_new

  SUBROUTINE init_usecase_pfts(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile, pft_tile

    INTEGER :: i

    INTEGER, PARAMETER :: npft = 11

    CHARACTER(len=5), PARAMETER :: pft_shortnames(npft) = [character(len=5) ::  &
      & 'pft01', 'pft02', 'pft03', 'pft04', 'pft05', 'pft06', &
      & 'pft07', 'pft08', 'pft09', 'pft10', 'pft11' ]

    CHARACTER(len=30), PARAMETER :: pft_longnames(npft) = [character(len=30) :: &
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
      & 'C3 crops'                    ]

    ! Index of PFT into lctlib columns
    ! id is one less than in lctlib file because the first column (glacier) has been thrown away
    INTEGER, PARAMETER :: pft_ids(npft) = (/ 1,2,3,4,9,10,11,12,14,15,19 /)

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_pfts'


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

        ! default jsb4: HYDRO_ is active at veg tile only to save runtime
        IF (.NOT. model%config%use_quincy) THEN
          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes       = (/HYDRO_,   SSE_,     ASSIMI_,   PHENO_/),    &
            & process_actions = (/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
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

          IF (model%processes(NLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(NLCC_, ON_TILE_, AGGREGATE_)
          END IF
          !
          ! use_quincy:
        !  HYDRO_ is active at PFT tiles to make all the HYDRO_ variables available to the PFT tiles
        !    note: HYDRO_ results differ only marginally between ON_TILE_ & ON_LEAFS_ (tested by Reiner)
        !  VEG_ is active
        ELSE
          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes       = (/HYDRO_,    SSE_,     ASSIMI_,   PHENO_,    VEG_   /), &
            & process_actions = (/ON_LEAFS_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')
          !$ACC ENTER DATA COPYIN(veg_tile)
        ENDIF

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

            DO i=1,npft
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
  !> This usecase consideres hydrologically different settings for all tiles below the land tile
  !! Here, the tile tree branches into several copys for, e.g. upland, lowland and floodplain settings
  SUBROUTINE init_usecase_hydrotiles(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, &
      & veg_tile, pft_tile

    INTEGER :: i,t

    INTEGER, PARAMETER :: ntopo = 3

    CHARACTER(len=10), PARAMETER :: topo_shortnames(ntopo) = [character(len=10) ::  &
      & 'up', 'low', 'plain']

    CHARACTER(len=50), PARAMETER :: topo_longnames(ntopo) = [character(len=50) :: &
      & 'upland fraction',  &
      & 'lowland fraction', &
      & 'floodplain fraction']

    INTEGER, PARAMETER :: macrotopo_types(ntopo) = [ &
      & UPLAND_TYPE, LOWLAND_TYPE, FLOODPLAIN_TYPE]

    INTEGER, PARAMETER :: npft = 11

    CHARACTER(len=5), PARAMETER :: pft_shortnames(npft) = [character(len=5) ::  &
      & 'pft01', 'pft02', 'pft03', 'pft04', 'pft05', 'pft06', &
      & 'pft07', 'pft08', 'pft09', 'pft10', 'pft11' ]

    CHARACTER(len=30), PARAMETER :: pft_longnames(npft) = [character(len=30) :: &
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
      & 'C3 crops'                    ]

    ! Index of PFT into lctlib columns
    ! id is one less than in lctlib file because the first column (glacier) has been thrown away
    INTEGER, PARAMETER :: pft_ids(npft) = (/ 1,2,3,4,9,10,11,12,14,15,19 /)

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_pfts'

    ! overview of tile structure:
    !
    !  box
    !   |
    !   |-------
    !   |      |
    !  land  lakes
    !   |
    !   |----------------------------------------------------------------------------------
    !   |                          |                            |                         |
    !  veg_upland                  veg_lowland                  veg_floodplain            glacier
    !   |                          |                            |
    !   |-------------             |--------------              |-----------------
    !   |            |             |             |              |                |
    !  pft01_upland  pftXX_upland  pft01_lowland pftXX_lowland  pft01_floodplain pftXX_floodplain
    !
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
          & processes       = (/SEB_,     HYDRO_  /), &
          & process_actions = (/ON_TILE_, ON_TILE_/), &
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
            & processes       = (/HYDRO_,   MIX_,     SSE_/), &
            & process_actions = (/ON_TILE_, ON_TILE_, ON_TILE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile)
          !$ACC ENTER DATA COPYIN(glacier_tile)
        END IF

        ! Create identical tile trees for the different hydrological settings
        DO t=1,ntopo

          ! MIX_, HYDRO_ and SSE_ are active at veg tiles only to save runtime
          veg_tile => t_jsb_tile('veg_'//TRIM(topo_shortnames(t)),                      &
            & 'Vegetation tile ('//TRIM(topo_longnames(t))//')', &
            & processes       = (/MIX_,     HYDRO_,   SSE_,     ASSIMI_,   PHENO_/),    &
            & process_actions = (/ON_TILE_, ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg_'//TRIM(topo_shortnames(t)))
          CALL veg_tile%Add_lct(t_jsb_lct(macrotopo_types(t)))
          !$ACC ENTER DATA COPYIN(veg_tile)

          IF (model%processes(ALCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(ALCC_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(FLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(FLCC_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(WLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(WLCC_, ON_TILE_, AGGREGATE_)
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

            DO i=1,npft
              pft_tile => t_jsb_tile(pft_shortnames(i)//'_'//TRIM(topo_shortnames(t)), &
                & TRIM(pft_longnames(i)//' ('//TRIM(topo_longnames(t))//')'),          &
                & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=veg_tile,         &
                & fract_varname='fract_'//pft_shortnames(i))
              !$ACC ENTER DATA COPYIN(pft_tile)
              NULLIFY(pft_tile)
            END DO

        ENDDO

    ! ! Next 3 lines only to avoid compiler warnings about pointers not being dereferenced
    ! IF (ASSOCIATED(glacier_tile)) CONTINUE
    ! IF (ASSOCIATED(lake_tile)) CONTINUE
    ! IF (ASSOCIATED(pft_tile)) CONTINUE

    NULLIFY(box_tile, land_tile, veg_tile, pft_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
    IF (model%config%use_lakes) NULLIFY(lake_tile)

  END SUBROUTINE init_usecase_hydrotiles

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
    INTEGER, PARAMETER :: npft = 11

    CHARACTER(len=5)  :: pft_shortname
    CHARACTER(len=10) :: age_class_name
    CHARACTER(len=10) :: init_ac_fracts_scheme

    CHARACTER(len=30), PARAMETER :: pft_longnames(npft) = [character(len=30) :: &
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
      & 'C3 crops'                    ]

    ! Index of PFT into lctlib columns
    ! id is one less than in lctlib file because the first column (glacier) has been thrown away
    INTEGER, PARAMETER :: pft_ids(npft) = (/ 1,2,3,4,9,10,11,12,14,15,19 /)

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_pfts_with_forest_age_classes'

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
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
            & processes = (/HYDRO_, SSE_/), &
            & process_actions = (/ON_TILE_, ON_TILE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile)
            ! & lct=t_jsb_lct(GLACIER_TYPE, lib_id=1), parent=land_tile)
        END IF

        !>veg_tile
        !!
        ! default jsb4: HYDRO_ & SSE_ are active at veg tile only to save runtime
          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes       = (/HYDRO_,   SSE_,     ASSIMI_,   PHENO_/),    &
            & process_actions = (/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')

          IF (model%processes(ALCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(ALCC_, ON_TILE_, AGGREGATE_)
          END IF
          !
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

            DO i=1,npft

              WRITE(pft_shortname, "(A3,I2.2)") "pft", i
              pft_tile => t_jsb_tile(pft_shortname, TRIM(pft_longnames(i)), &
                & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=veg_tile)

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

    ! CALL test_pools()

  END SUBROUTINE init_usecase_pfts_with_forest_age_classes

#endif
END MODULE mo_jsb_model_usecases
