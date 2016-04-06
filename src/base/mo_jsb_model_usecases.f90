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
  USE mo_jsb_process_factory, ONLY: ATM2LAND, LAND2ATM, SEB, TURB, SOIL, HYDRO, VEG, RAD, HD
  USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, ON_SUBTREE_
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid !, get_nblks, get_nproma
  USE mo_jsb_grid,            ONLY: Get_Grid
  USE mo_jsb_tile,            ONLY: t_jsb_tile
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct, &
    &                               BOX_TYPE, LAND_TYPE, GLACIER_TYPE, PFT_TYPE, LAKE_TYPE

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

    CASE ('jsbach_lite_with_hd')
      CALL init_usecase_lite_with_hd(model)

    CASE ('jsbach_lite_tpe_open')
      CALL init_usecase_lite_tpe_open(model)

    CASE ('jsbach_lite_tpe_closed')
      CALL init_usecase_lite_tpe_closed(model)

    CASE ('jsbach_lite_with_lake')
      CALL init_usecase_lite_with_lake(model)

    CASE ('jsbach_tiles')
      CALL init_usecase_tiles(model)

    CASE DEFAULT
      CALL finish(TRIM(routine), 'JSBACH usecase not defined.')

    END SELECT

  END SUBROUTINE init_usecase

  SUBROUTINE init_usecase_lite(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile
    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite'

    grid => Get_grid(model%grid_id)

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile',            &
      & processes=(/ATM2LAND,LAND2ATM,RAD,SEB,TURB,SOIL,HYDRO,VEG/), &
      & varlist_name=model%shortname, &
      & lct=t_jsb_lct(LAND_TYPE), &
      & model_id=model%id, grid_id=model%grid_id)

    model%top => box_tile

    CALL box_tile%Print()

    NULLIFY(box_tile)

  END SUBROUTINE init_usecase_lite

  SUBROUTINE init_usecase_lite_with_hd(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile),  POINTER :: box_tile
    TYPE(t_jsb_grid),   POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite_with_hd'

    grid => Get_grid(model%grid_id)

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile',            &
      & processes=(/ATM2LAND,LAND2ATM,RAD,SEB,TURB,SOIL,HYDRO,VEG,HD/), &
      & varlist_name=model%shortname, &
      & lct=t_jsb_lct(LAND_TYPE), &
      & model_id=model%id, grid_id=model%grid_id)

    model%top => box_tile

    CALL box_tile%Print()

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

  SUBROUTINE init_usecase_lite_with_lake(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile),  POINTER :: box_tile, lake_tile, land_tile
    TYPE(t_jsb_grid),   POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite_with_lake'

    grid => Get_grid(model%grid_id)

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile',            &
      & processes=(/ATM2LAND,LAND2ATM,RAD,TURB/), &
      & varlist_name=model%shortname, &
      & lct=t_jsb_lct(LAND_TYPE), &
      & model_id=model%id, grid_id=model%grid_id)
    model%top => box_tile

      land_tile => t_jsb_tile('land', 'Land tile', &
        & processes=(/SEB,HYDRO,SOIL,VEG/), process_actions=(/ON_TILE_,ON_TILE_,ON_LEAFS_,ON_LEAFS_/), &
        & varlist_name=model%shortname, lct=t_jsb_lct(LAND_TYPE), parent=box_tile)

      lake_tile => t_jsb_tile('lake', 'Lake tile', &
        & processes=(/SEB/), process_actions=(/ON_TILE_/), &
        & varlist_name=model%shortname, lct=t_jsb_lct(LAKE_TYPE), parent=box_tile)

    CALL box_tile%Print()
      CALL land_tile%Print()
      CALL lake_tile%Print()

    NULLIFY(box_tile)

  END SUBROUTINE init_usecase_lite_with_lake

  SUBROUTINE init_usecase_tiles(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, soilveg_tile
    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_tiles'

    grid => Get_grid(model%grid_id)

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', processes=(/ATM2LAND,LAND2ATM,RAD,TURB/), &
      varlist_name=model%shortname, lct=t_jsb_lct(BOX_TYPE), model_id=model%id, grid_id=model%grid_id)
    model%top => box_tile

      land_tile => t_jsb_tile('land', 'Land tile', &
        & processes=(/SEB,HYDRO/), process_actions=(/ON_TILE_,ON_TILE_/), &
!!$        & processes=(/SRF,RAD/), process_actions=(/ON_TILE_,ON_LEAFS_/), &
!!$                            & processes=(/SRF/), process_actions=(/ON_SUBTREE_/), &
        & varlist_name=model%shortname, lct=t_jsb_lct(LAND_TYPE), parent=box_tile)

        glacier_tile => t_jsb_tile('glac', 'Glacier tile', processes=(/SOIL/), &
          & varlist_name=model%shortname, lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile)

        soilveg_tile => t_jsb_tile('soilveg', 'Soil/veg/glacier tile', processes=(/SOIL,VEG/), &
          & varlist_name=model%shortname, lct=t_jsb_lct(PFT_TYPE), parent=land_tile)

      lake_tile => t_jsb_tile('lake', 'Lake tile', &
        & processes=(/SEB/), process_actions=(/ON_TILE_/), &
        & varlist_name=model%shortname, lct=t_jsb_lct(LAKE_TYPE), parent=box_tile)

    CALL box_tile%Print()
      CALL land_tile%Print()
        CALL glacier_tile%Print()
        CALL soilveg_tile%Print()
      CALL lake_tile%Print()

    NULLIFY(box_tile)

  END SUBROUTINE init_usecase_tiles

#endif
END MODULE mo_jsb_model_usecases
