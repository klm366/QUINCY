!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Initialization of the the soil_biogeochemistry memory.
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version for SLM by Jan Engel (20161208) based on the JBACH4 version by Reiner Schnur (2013-04-04)
!!
MODULE mo_sb_init

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message, finish, message_text
  USE mo_jsb_control,             ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sb_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize soil_biogeochemistry process (after memory has been set up)
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE sb_init(tile)

    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_class,             ONLY: get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: SB_
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid


    ! Use of process memories
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct

    TYPE(t_jsb_model),      POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib        !< land-cover-type library - parameter across pft's
    REAL(wp)                                     :: hlp1
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init'

    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)

    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of vegetation memory for tile '//TRIM(tile%name))

    ! define var

    model  => Get_model(tile%owner_model_id)
    IF (.NOT. model%config%use_quincy) RETURN       ! only with QUINCY model
    IF (tile%lcts(1)%lib_id == 0) RETURN            ! run this init only if the present tile is a pft
    lctlib => model%lctlib(tile%lcts(1)%lib_id)

    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)

    ! get var



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    hlp1 = 0.5_wp


  END SUBROUTINE sb_init

END MODULE mo_sb_init
