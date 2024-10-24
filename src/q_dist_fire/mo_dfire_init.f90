!>
!>  @par Copyright
!>    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!>    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!>    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!>
!>    This file is part of the QUINCY model.
!>    The QUINCY model is free software: you can redistribute it and/or modify
!>    it under the terms of the GNU General Public License version 3,
!>    as published by the Free Software Foundation.
!>
!>    The QUINCY model is distributed in the hope that it will be useful,
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!>    GNU General Public License for more details.
!>
!>    You should have received a copy of the GNU General Public License
!>    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!> 
!> @brief Initialization of the the dfire memory.
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
!> @par Revision History
!>
MODULE mo_dfire_init
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message
  USE mo_jsb_control,             ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dfire_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_init'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Intialize disturbance fire process (after memory has been set up)
  !>
  SUBROUTINE dfire_init(tile, options)
    USE mo_jsb_class,               ONLY: get_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: DFIRE_
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_dfire_constants,         ONLY: live_frac_allocation
    dsl4jsb_Use_config(DFIRE_)
    dsl4jsb_Use_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib !< land-cover-type library - parameter across pfts
    INTEGER                                      :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':dfire_init'
    ! ----------------------------------------------------------------------------------------------------- !
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(DFIRE_)
    dsl4jsb_Def_memory(DFIRE_)                          !   TYPE(t_dfire_memory), POINTER      :: dfire__mem
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: dfire_pm_id                                  !< Selected precipitation method
    ! ----------------------------------------------------------------------------------------------------- !
!    dsl4jsb_Real2D_onChunk      :: nesterov_index
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_active(DFIRE_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    model  => get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(DFIRE_)
    dsl4jsb_Get_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
!    dsl4jsb_Get_var2D_onChunk(DFIRE_, nesterov_index)
    ! ----------------------------------------------------------------------------------------------------- !
!    nesterov_index(:) = 0._wp

  END SUBROUTINE dfire_init
#endif
END MODULE mo_dfire_init
