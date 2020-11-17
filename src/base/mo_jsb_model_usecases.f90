!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
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

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_process_class,   ONLY: A2L_, L2A_, RAD_, ASSIMI_, VEG_, PHENO_, SSE_, SB_, &
                                    ON_LEAFS_, ON_TILE_
  USE mo_jsb_tile,            ONLY: t_jsb_tile

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

    CLASS(t_jsb_tile), POINTER :: box_tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite'

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_,      ASSIMI_,   VEG_,      PHENO_,    SSE_,     SB_/), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, grid_id=model%grid_id, fract_varname='notsea')
    model%top => box_tile

    NULLIFY(box_tile)

  END SUBROUTINE init_usecase_lite

#endif
END MODULE mo_jsb_model_usecases
