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
!! @brief Contains the utilities for the soil_biogeochemistry process
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20170104)
!!
MODULE mo_sb_util

  USE mo_kind,        ONLY: wp
  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message, finish

  ! Use of processes in this module
  dsl4jsb_Use_processes SB_

  ! Use of process memories
  dsl4jsb_Use_memory(SB_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: reset_sb_fluxes

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_util'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> init/reset sb fluxes (with/to zero)
  !! 
  !! called prior to all Tasks
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE reset_sb_fluxes(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_math_constants,    ONLY: zero

    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options

    TYPE(t_jsb_model),      POINTER       :: model
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_sb_fluxes'

    dsl4jsb_Def_memory(SB_)                               ! TYPE(t_sb_memory),       POINTER :: sb__mem

    ! 2D SB_
    dsl4jsb_Real2D_onChunk :: total_flux_carbon
    dsl4jsb_Real2D_onChunk :: total_flux_nitrogen
    dsl4jsb_Real2D_onChunk :: total_flux_phosphorus
    dsl4jsb_Real2D_onChunk :: total_flux_carbon13
    dsl4jsb_Real2D_onChunk :: total_flux_carbon14
    dsl4jsb_Real2D_onChunk :: total_flux_nitrogen15

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(SB_)

    model  => Get_model(tile%owner_model_id)

    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon)           ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_phosphorus)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon13)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon14)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen15)       ! out

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 1.0 initialise/reset fluxes
    !!
    ! 2D fields  
    total_flux_carbon(:)         = zero
    total_flux_nitrogen(:)       = zero
    total_flux_phosphorus(:)     = zero
    total_flux_carbon13(:)       = zero
    total_flux_carbon14(:)       = zero
    total_flux_nitrogen15(:)     = zero

  END SUBROUTINE reset_sb_fluxes

!#endif
END MODULE mo_sb_util
