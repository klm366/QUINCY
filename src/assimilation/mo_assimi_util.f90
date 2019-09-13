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
!! @brief diverse subroutines and functions for the assimilation process
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2018-04-06)
!!
MODULE mo_assimi_util

#include "dsl_jsbach.inc"

  USE mo_kind,        ONLY: wp
  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: reset_assimi_fluxes    ! called in the land model interface prior to all Tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_util'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> init/reset assimilation fluxes (with/to zero)
  !! 
  !! called in the land model interface prior to all Tasks
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE reset_assimi_fluxes(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_math_constants,    ONLY: zero
    
    ! Use of process memories
    dsl4jsb_Use_memory(ASSIMI_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local 
    TYPE(t_jsb_model),      POINTER       :: model
    !TYPE(t_lctlib_element), POINTER       :: lctlib         !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_assimi_fluxes'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ASSIMI_) 
    ! Declare pointers to variables in memory
    ! ASSIMI_ 2D
    dsl4jsb_Real2D_onChunk      :: gross_assimilation
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C13
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C14
    dsl4jsb_Real2D_onChunk      :: net_assimilation
    dsl4jsb_Real2D_onChunk      :: net_assimilation_boc
    dsl4jsb_Real2D_onChunk      :: maint_respiration_leaf
    ! ASSIMI_ 3D
    dsl4jsb_Real3D_onChunk      :: gross_assimilation_cl
    dsl4jsb_Real3D_onChunk      :: net_assimilation_cl
    dsl4jsb_Real3D_onChunk      :: maint_respiration_leaf_cl
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(ASSIMI_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(PHENO_)
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_) 
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    !lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  gross_assimilation)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  gross_assimilation_C13)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  gross_assimilation_C14)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  net_assimilation)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  net_assimilation_boc)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  maint_respiration_leaf)
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,  gross_assimilation_cl)
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,  net_assimilation_cl)
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,  maint_respiration_leaf_cl)



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! set zero all flux variables of the assimilation process
    gross_assimilation(:)           = zero
    gross_assimilation_C13(:)       = zero
    gross_assimilation_C14(:)       = zero
    net_assimilation(:)             = zero
    net_assimilation_boc(:)         = zero
    maint_respiration_leaf(:)       = zero

    gross_assimilation_cl(:,:)      = zero
    net_assimilation_cl(:,:)        = zero
    maint_respiration_leaf_cl(:,:)  = zero


  END SUBROUTINE reset_assimi_fluxes

END MODULE mo_assimi_util
