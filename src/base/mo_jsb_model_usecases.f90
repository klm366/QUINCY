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
  USE mo_jsb_process_class,   ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, RAD_, HD_, ASSIMI_, PHENO_
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

    CASE ('jsbach_lite_tpe_open')
      CALL init_usecase_lite_tpe_open(model)

    CASE ('jsbach_lite_tpe_closed')
      CALL init_usecase_lite_tpe_closed(model)

    CASE ('jsbach_tiles')
      CALL init_usecase_tiles(model)

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
      & processes      =(/ATM2LAND_, LAND2ATM_, RAD_,       HYDRO_,     TURB_     /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, grid_id=model%grid_id, fract_filename='bc_land_frac.nc', fract_varname='notsea')
    model%top => box_tile

      ! Create sub-tiles
      lake_tile => t_jsb_tile('lake', 'Lake tile', &
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile, fract_filename='bc_land_frac.nc', fract_varname='fract_lake')

      land_tile => t_jsb_tile('land', 'Land tile', &
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        & lct=t_jsb_lct(LAND_TYPE), &
        & parent=box_tile, &
        & fract_filename='bc_land_frac.nc', fract_varname='fract_land')

        glacier_tile => t_jsb_tile('glac', 'Glacier tile', processes=(/SOIL_/), &
          & lct=t_jsb_lct(GLACIER_TYPE, lib_id=1), parent=land_tile, &
          & fract_filename='bc_land_frac.nc', fract_varname='fract_glac')

        veg_tile => t_jsb_tile('veg', 'Veg tile', processes=(/SOIL_, VEG_, PHENO_/), &
          & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
          & fract_filename='bc_land_frac.nc', fract_varname='landnoglacier')

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
      & processes      =(/ATM2LAND_, LAND2ATM_, RAD_,      HYDRO_,    TURB_     /), &
      & process_actions=(/ON_TILE_,  ON_TILE_,  ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, grid_id=model%grid_id, fract_filename='bc_land_frac.nc', fract_varname='notsea')
    model%top => box_tile

      ! Create sub-tiles
      land_tile => t_jsb_tile('land', 'Land tile', &
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        & lct=t_jsb_lct(LAND_TYPE), &
        & parent=box_tile, &
        & fract_filename='bc_land_frac.nc', fract_varname='fract_land')

        glacier_tile => t_jsb_tile('glac', 'Glacier tile', processes=(/SOIL_/), &
          & lct=t_jsb_lct(GLACIER_TYPE, lib_id=1), parent=land_tile, &
          & fract_filename='bc_land_frac.nc', fract_varname='fract_glac')

        veg_tile => t_jsb_tile('veg', 'Veg tile', processes=(/SOIL_, VEG_, PHENO_/), &
          & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
          & fract_filename='bc_land_frac.nc', fract_varname='landnoglacier')

    NULLIFY(box_tile)

  END SUBROUTINE init_usecase_lite_nolake

  SUBROUTINE init_usecase_lite_with_hd(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile),  POINTER :: box_tile, lake_tile, glacier_tile, veg_tile
    TYPE(t_jsb_grid),   POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite_with_hd'

    grid => Get_grid(model%grid_id)

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/ATM2LAND_, LAND2ATM_, RAD_,       HYDRO_,     TURB_,      SEB_,      HD_/), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_TILE_, ON_TILE_/), &
      & lct=t_jsb_lct(LAND_TYPE), &
      & model_id=model%id, grid_id=model%grid_id, fract_filename='bc_land_frac.nc', fract_varname='notsea')
    model%top => box_tile
    CALL box_tile%Print()

      lake_tile => t_jsb_tile('lake', 'Lake tile', &
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile, fract_filename='bc_land_frac.nc')
      CALL lake_tile%Print()

      glacier_tile => t_jsb_tile('glac', 'Glacier tile', processes=(/SOIL_/), &
        & lct=t_jsb_lct(GLACIER_TYPE, lib_id=1), parent=box_tile, fract_filename='bc_land_frac.nc')
      CALL glacier_tile%Print()

      veg_tile => t_jsb_tile('veg', 'Veg tile', processes=(/SOIL_, VEG_, PHENO_/), &
        & lct=t_jsb_lct(VEG_TYPE), parent=box_tile,  &
        & fract_filename='bc_land_frac.nc', fract_varname='fract_land')
      CALL veg_tile%Print()


    ! box_tile => t_jsb_tile('box', 'Top tile',            &
    !   & processes=(/ATM2LAND_,LAND2ATM_,RAD_,SEB_,TURB_,SOIL_,HYDRO_,PHENO_,VEG_,HD_/), &
    !   & lct=t_jsb_lct(LAND_TYPE), &
    !   & model_id=model%id, grid_id=model%grid_id, fract_filename='bc_land_frac.nc', fract_varname='notsea')

    ! model%top => box_tile

    ! CALL box_tile%Print()

    NULLIFY(box_tile)

  END SUBROUTINE init_usecase_lite_with_hd

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

  SUBROUTINE init_usecase_tiles(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile, pft_tile
    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_tiles'

    grid => Get_grid(model%grid_id)

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/ATM2LAND_, LAND2ATM_, RAD_,       HYDRO_,     TURB_,      HD_      /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_TILE_/), &
      & model_id=model%id, grid_id=model%grid_id, fract_filename='bc_land_frac.nc', fract_varname='notsea')
    model%top => box_tile

      land_tile => t_jsb_tile('land', 'Land tile', &
        & processes=(/SEB_, SOIL_/), process_actions=(/ON_TILE_, ON_LEAFS_/), &
!!$        & processes=(/SRF,RAD_/), process_actions=(/ON_TILE_,ON_LEAFS_/), &
!!$                            & processes=(/SRF/), process_actions=(/ON_SUBTREE_/), &
        & lct=t_jsb_lct(LAND_TYPE), &
        & parent=box_tile, &
        & fract_filename='bc_land_frac.nc')

        glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
          & lct=t_jsb_lct(GLACIER_TYPE, lib_id=1), parent=land_tile, fract_filename='bc_land_frac.nc')

        veg_tile => t_jsb_tile('veg', 'Veg tile', processes=(/VEG_, ASSIMI_, PHENO_/), &
          & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
          & fract_filename='bc_land_frac.nc', fract_varname='fract_land')

          pft_tile => t_jsb_tile('pft1', 'Tropical evergreen trees', &
            & lct=t_jsb_lct(VEG_TYPE, lib_id=2), parent=veg_tile, fract_filename='')

          pft_tile => t_jsb_tile('pft2', 'Tropical deciduous trees', &
            & lct=t_jsb_lct(VEG_TYPE, lib_id=3), parent=veg_tile, fract_filename='')

      lake_tile => t_jsb_tile('lake', 'Lake tile', &
        & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
        & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile, fract_filename='bc_land_frac.nc')

    NULLIFY(box_tile)

  END SUBROUTINE init_usecase_tiles

#endif
END MODULE mo_jsb_model_usecases
