!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_nix_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  USE mo_jsb_model_class,   ONLY: t_jsb_model, MODEL_QUINCY, MODEL_JSBACH
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid
  USE mo_jsb_grid,          ONLY: Get_grid
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,            ONLY: missval

  dsl4jsb_Use_processes NIX_
  dsl4jsb_Use_config(NIX_)
!  dsl4jsb_Use_memory(NIX_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: nix_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_nix_init'

CONTAINS

  !
  !> Intialize nix process (after memory has been set up)
  !
  SUBROUTINE nix_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':nix_init'

    model => get_model(tile%owner_model_id)

    ! ---------------------------
    ! Go   


    CALL message (TRIM(modname), 'Hello NIX!')


  END SUBROUTINE nix_init




#endif
END MODULE mo_nix_init








