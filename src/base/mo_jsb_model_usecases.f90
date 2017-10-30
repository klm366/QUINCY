!>
!! @brief Contains methods for initialization of JSBACH usecases.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_model_usecases
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message, finish
  USE mo_kind,                ONLY: wp

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_process_class,   ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, RAD_, HD_, ASSIMI_, PHENO_
  USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, ON_SUBTREE_
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid !, get_nblks, get_nproma
  USE mo_jsb_grid,            ONLY: Get_Grid
  USE mo_jsb_tile,            ONLY: t_jsb_tile
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct, &
    &                               LAND_TYPE, BARE_TYPE, VEG_TYPE, GLACIER_TYPE, LAKE_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_usecase

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_usecases'

CONTAINS

  SUBROUTINE init_usecase(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase'

    CALL message(TRIM(routine), 'Initializing JSBACH usecase "'//TRIM(model%config%usecase)//'"')

    SELECT CASE (TRIM(model%config%usecase))

    CASE ('jsbach_lite')
      CALL init_usecase_lite(model)

    CASE ('jsbach_lite_nolake')
      CALL init_usecase_lite_nolake(model)

    CASE ('jsbach_lite_with_hd')
      CALL init_usecase_lite_with_hd(model)

    CASE ('jsbach_lite_nolake_with_hd')
      CALL init_usecase_lite_nolake_with_hd(model)

    CASE ('jsbach_lite_tpe_open')
      CALL init_usecase_lite_tpe_open(model)

    CASE ('jsbach_lite_tpe_closed')
      CALL init_usecase_lite_tpe_closed(model)

    CASE ('jsbach_pfts')
      CALL init_usecase_pfts(model)

    CASE ('jsbach_pfts_with_hd')
      CALL init_usecase_pfts_with_hd(model)

    CASE DEFAULT
      CALL finish(TRIM(routine), 'JSBACH usecase not defined.')

    END SELECT

  END SUBROUTINE init_usecase

  SUBROUTINE init_usecase_lite(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile
    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite'

    grid => Get_grid(model%grid_id)

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_,      HYDRO_,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, grid_id=model%grid_id, fract_varname='notsea')
    model%top => box_tile

      ! Create sub-tiles
      lake_tile => t_jsb_tile('lake', 'Lake tile', &
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile, fract_varname='fract_lake')

      land_tile => t_jsb_tile('land', 'Land tile', &
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)

        glacier_tile => t_jsb_tile('glac', 'Glacier tile', processes=(/SSE_/), &
          & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile, &
          & fract_varname='fract_glac')

        veg_tile => t_jsb_tile('veg', 'Veg tile', processes=(/SSE_, PHENO_/), &
          & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
          & fract_varname='fract_veg')

    NULLIFY(box_tile)

  END SUBROUTINE init_usecase_lite

  SUBROUTINE init_usecase_lite_nolake(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, land_tile, glacier_tile, veg_tile
    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite_nolake'

    grid => Get_grid(model%grid_id)

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_,      HYDRO_,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, grid_id=model%grid_id, fract_varname='notsea')
    model%top => box_tile

      ! Create sub-tiles
      land_tile => t_jsb_tile('land', 'Land tile', &
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        & lct=t_jsb_lct(LAND_TYPE), &
        & parent=box_tile, &
        & fract_varname='equal')   ! Fraction = 1.

        glacier_tile => t_jsb_tile('glac', 'Glacier tile', processes=(/SSE_/), &
          & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile, &
          & fract_varname='fract_glac')

        veg_tile => t_jsb_tile('veg', 'Veg tile', processes=(/SSE_, PHENO_/), &
          & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
          & fract_varname='fract_veg')

    NULLIFY(box_tile)

  END SUBROUTINE init_usecase_lite_nolake

  SUBROUTINE init_usecase_lite_with_hd(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile),  POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile
    TYPE(t_jsb_grid),   POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite_with_hd'

    grid => Get_grid(model%grid_id)

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_     , HYDRO_   , TURB_    , HD_     /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_TILE_/), &
      & lct=t_jsb_lct(LAND_TYPE), &
      & model_id=model%id, grid_id=model%grid_id, fract_varname='notsea')
    model%top => box_tile

      ! Create sub-tiles
      lake_tile => t_jsb_tile('lake', 'Lake tile', &
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile, fract_varname='fract_lake')

      land_tile => t_jsb_tile('land', 'Land tile', &
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        & lct=t_jsb_lct(LAND_TYPE), &
        & parent=box_tile, &
        & fract_varname='fract_land')

        glacier_tile => t_jsb_tile('glac', 'Glacier tile', processes=(/SSE_/), &
          & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile, &
          & fract_varname='fract_glac')

        veg_tile => t_jsb_tile('veg', 'Veg tile', processes=(/SSE_, PHENO_/), &
          & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
          & fract_varname='fract_veg')

    NULLIFY(box_tile)

  END SUBROUTINE init_usecase_lite_with_hd

  SUBROUTINE init_usecase_lite_nolake_with_hd(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, land_tile, glacier_tile, veg_tile
    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite_nolake_with_hd'

    grid => Get_grid(model%grid_id)

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_,      HYDRO_,    TURB_    , HD_     /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_TILE_/), &
      & model_id=model%id, grid_id=model%grid_id, fract_varname='notsea')
    model%top => box_tile

      ! Create sub-tiles
      land_tile => t_jsb_tile('land', 'Land tile', &
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        & lct=t_jsb_lct(LAND_TYPE), &
        & parent=box_tile, &
        & fract_varname='equal')   ! Fraction = 1.

        glacier_tile => t_jsb_tile('glac', 'Glacier tile', processes=(/SSE_/), &
          & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile, &
          & fract_varname='fract_glac')

        veg_tile => t_jsb_tile('veg', 'Veg tile', processes=(/SSE_, PHENO_/), &
          & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
          & fract_varname='fract_veg')

    NULLIFY(box_tile)

  END SUBROUTINE init_usecase_lite_nolake_with_hd

  SUBROUTINE init_usecase_lite_tpe_open(model)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite_tpe_open'

#ifdef HAVE_F2003
    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model
#else
    TYPE(t_jsb_model), POINTER                :: model
#endif

    CALL init_usecase_lite(model)

  END SUBROUTINE init_usecase_lite_tpe_open

  SUBROUTINE init_usecase_lite_tpe_closed(model)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite_tpe_closed'

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CALL init_usecase_lite(model)

  END SUBROUTINE init_usecase_lite_tpe_closed

  SUBROUTINE init_usecase_pfts(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile, pft_tile
    TYPE(t_jsb_grid),  POINTER :: grid

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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_pfts'

    grid => Get_grid(model%grid_id)

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_  ,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, grid_id=model%grid_id, fract_varname='notsea')
    model%top => box_tile

      lake_tile => t_jsb_tile('lake', 'Lake tile', &
        & processes=(/SEB_, HYDRO_/), process_actions=(/ON_TILE_, ON_TILE_/), &
        & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile)

      land_tile => t_jsb_tile('land', 'Land tile', &
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)

        glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
          & processes = (/HYDRO_, SSE_/), &
          & process_actions = (/ON_TILE_, ON_TILE_/), &
          & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile)
          ! & lct=t_jsb_lct(GLACIER_TYPE, lib_id=1), parent=land_tile)

        veg_tile => t_jsb_tile('veg', 'Veg tile', &
          & processes       = (/HYDRO_,   SSE_,     ASSIMI_,   PHENO_/), &
          & process_actions = (/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
          & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
          & fract_varname='fract_veg')

          DO i=1,npft
            pft_tile => t_jsb_tile(pft_shortnames(i), TRIM(pft_longnames(i)), &
             & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=veg_tile)
          END DO
          ! DO i=1,npft
          !   pft_tile => t_jsb_tile(pft_shortnames(i), TRIM(pft_longnames(i)), &
          !     & processes=(/ASSIMI_, PHENO_/), &
          !     & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=land_tile)
          ! END DO

    NULLIFY(box_tile)
    NULLIFY(lake_tile)
    NULLIFY(land_tile)
    NULLIFY(glacier_tile)
    NULLIFY(veg_tile)
    NULLIFY(pft_tile)

  END SUBROUTINE init_usecase_pfts

  SUBROUTINE init_usecase_pfts_with_hd(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile, pft_tile
    TYPE(t_jsb_grid),  POINTER :: grid

    INTEGER :: i

    INTEGER, PARAMETER :: npft = 11

    CHARACTER(len=5), PARAMETER :: pft_shortnames(npft) = [character(len=5) :: &
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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_pfts_with_hd'

    grid => Get_grid(model%grid_id)

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_  ,    TURB_    , HD_     /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_TILE_/), &
      & model_id=model%id, grid_id=model%grid_id, fract_varname='notsea')
    model%top => box_tile

      lake_tile => t_jsb_tile('lake', 'Lake tile', &
        & processes=(/SEB_, HYDRO_/), process_actions=(/ON_TILE_, ON_TILE_/), &
        & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile)

      land_tile => t_jsb_tile('land', 'Land tile', &
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)

        glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
          & processes = (/HYDRO_, SSE_/), &
          & process_actions = (/ON_TILE_, ON_TILE_/), &
          & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile)
          ! & lct=t_jsb_lct(GLACIER_TYPE, lib_id=1), parent=land_tile)

        veg_tile => t_jsb_tile('veg', 'Veg tile', &
          & processes       = (/HYDRO_,   SSE_,     ASSIMI_,   PHENO_/), &
          & process_actions = (/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
          & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
          & fract_varname='fract_veg')

          DO i=1,npft
            pft_tile => t_jsb_tile(pft_shortnames(i), TRIM(pft_longnames(i)), &
             & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=veg_tile)
          END DO
          ! DO i=1,npft
          !   pft_tile => t_jsb_tile(pft_shortnames(i), TRIM(pft_longnames(i)), &
          !     & processes=(/ASSIMI_, PHENO_/), &
          !     & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=land_tile)
          ! END DO

    NULLIFY(box_tile)
    NULLIFY(lake_tile)
    NULLIFY(land_tile)
    NULLIFY(glacier_tile)
    NULLIFY(veg_tile)
    NULLIFY(pft_tile)

  END SUBROUTINE init_usecase_pfts_with_hd

#endif
END MODULE mo_jsb_model_usecases
