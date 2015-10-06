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

  USE mo_exception,       ONLY: message, finish
  USE mo_kind,            ONLY: wp

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_process_factory, ONLY: ATM2LAND, LAND2ATM, SRF, SOIL, VEG, GLACIER, LAKE, HD
  USE mo_hsm_class,          ONLY: t_Hsm, t_State
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid !, get_nblks, get_nproma
  USE mo_jsb_grid,           ONLY: Get_Grid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile

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

    CLASS(t_State),   POINTER :: box_tile
    TYPE(t_jsb_grid), POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite'

    grid => Get_grid(model%grid_id)

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', processes=(/ATM2LAND,SRF,SOIL,VEG,HD/), &
      varlist_name=model%shortname, model_id=model%id, grid_id=model%grid_id)
    model%hsm => t_Hsm('hsm', top=box_tile)

    NULLIFY(box_tile)

  END SUBROUTINE init_usecase_lite

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

    CLASS(t_State),   POINTER :: box_tile, glacier_tile, lake_tile
    TYPE(t_jsb_grid), POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_tiles'

    grid => Get_grid(model%grid_id)

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', processes=(/ATM2LAND,SRF,SOIL,VEG,HD/), &
                           varlist_name=model%shortname, model_id=model%id, grid_id=model%grid_id)
    model%hsm => t_Hsm('hsm', box_tile)

!!$    glacier_tile => t_jsb_tile('glacier', 'Glacier tile', processes=(/SRF,GLACIER/), &
!!$                               varlist_name=model%shortname, parent=box_tile)
!!$
!!$    lake_tile => t_jsb_tile('lake', 'Lake tile', processes=(/SRF,GLACIER/), &
!!$                            varlist_name=model%shortname, parent=box_tile)

    NULLIFY(box_tile)

  END SUBROUTINE init_usecase_tiles

#endif
END MODULE mo_jsb_model_usecases
