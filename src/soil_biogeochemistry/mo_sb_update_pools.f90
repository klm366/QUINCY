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
!! @brief Contains the routines for the Task update_sb_pools
!!
!! @author
!!  Silvia Caldararu, Sönke Zaehle 
!!
!! @par Revision History
!! First version by Sönke Zaehle (20170105): copyied from vegetation/mo_veg_soil.f90
!!
MODULE mo_sb_update_pools

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message, message_text, finish
  USE mo_jsb_math_constants,      ONLY: eps8

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_sb_pools

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_update_pools'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil biogeochemical pools 
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sb_pools(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element    
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    ! USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_process_class,     ONLY: SB_
    USE mo_jsb_math_constants,    ONLY: dtime, eps8

    ! Use of process configurations and memory
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)

    IMPLICIT NONE

    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options

    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib               !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_pools'

    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)

    ! Declare pointers to variables in memory
    ! SB_ 2D
    dsl4jsb_Real2D_onChunk :: total_flux_carbon
    dsl4jsb_Real2D_onChunk :: total_flux_nitrogen
    dsl4jsb_Real2D_onChunk :: total_flux_phosphorus
    dsl4jsb_Real2D_onChunk :: total_flux_carbon13
    dsl4jsb_Real2D_onChunk :: total_flux_carbon14
    dsl4jsb_Real2D_onChunk :: total_flux_nitrogen15

    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc

    IF (.NOT. tile%Is_process_active(SB_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)

    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)

    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon)           ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_phosphorus)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon13)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon14)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen15)       ! out

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    total_flux_carbon     = 0.5_wp
    total_flux_nitrogen   = 0.5_wp
    total_flux_phosphorus = 0.5_wp
    total_flux_carbon13   = 0.5_wp
    total_flux_carbon14   = 0.5_wp
    total_flux_nitrogen15 = 0.5_wp

  END SUBROUTINE update_sb_pools

END MODULE mo_sb_update_pools
