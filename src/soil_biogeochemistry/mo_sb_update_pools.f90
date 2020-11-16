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
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_jsb_process_class,     ONLY: SB_   ! A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
    USE mo_jsb_pool_class,        ONLY: ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15 
    USE mo_jsb_math_constants,    ONLY: dtime, eps8

    ! Use of pools
    dsl4jsb_Use_pools(SB_) t_sb_main_pool, t_sb_components_pool, t_sb_elements_pool
    
    ! Use of process configurations and memory
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)

    IMPLICIT NONE

    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options

    TYPE(t_jsb_model),      POINTER         :: model
    TYPE(t_lctlib_element), POINTER         :: lctlib               !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER         :: vgrid_soil_sb        !< Vertical grid
    INTEGER                                 :: nsoil_sb             !< number of soil layers as used by the SB_ process (see SSE_ config)
    INTEGER                                 :: iblk, ics, ice, nc
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: hlp1
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: hlp10
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_pools'

    ! Declare process configuration, memory, and pool-structure Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_pool(t_sb_components_pool)  :: sb_pool
    dsl4jsb_Def_pool(t_sb_components_pool)  :: sb_formation
    dsl4jsb_Def_pool(t_sb_components_pool)  :: sb_loss
    dsl4jsb_Def_pool(t_sb_components_pool)  :: sb_transport
    dsl4jsb_Def_pool(t_sb_elements_pool) :: sb_delta_pools
    dsl4jsb_Def_pool(t_sb_elements_pool) :: sb_mycorrhiza_export

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
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels

    ALLOCATE(hlp1(nc, nsoil_sb))
    ALLOCATE(hlp10(nc, nsoil_sb,1))     ! (nc, vgrid, iblk)

    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_pool,               sb_pool)
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_formation,          sb_formation)
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_loss,               sb_loss)
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_transport,          sb_transport)
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_delta_pools,        sb_delta_pools)
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_mycorrhiza_export,  sb_mycorrhiza_export)

    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon)           ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_phosphorus)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon13)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon14)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen15)       ! out

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !!! playground for testing !!!

    !! variables
    ! real2d variables
    total_flux_carbon           = 0.5_wp
    total_flux_nitrogen         = 0.5_wp
    total_flux_phosphorus       = 0.5_wp
    total_flux_carbon13         = 0.5_wp
    total_flux_carbon14         = 0.5_wp
    total_flux_nitrogen15       = 0.5_wp

    !! pools
    hlp1(:,:)                   = 2.0_wp
    hlp10(:,:,:)                = 4.0_wp
    sb_pool%dom                 = hlp1(:,:)
    sb_pool%fungi               = 5.0_wp
    sb_pool%residue             = sb_pool%dom + sb_pool%fungi - 1.0_wp                                ! 6.0
    sb_pool%mycorrhiza          = (sb_pool%residue / hlp1(:,:)) + 17.0_wp - (5.0_wp * sb_pool%fungi)  ! -5.0

    sb_loss                     = 0.2_wp
    sb_loss                     = hlp1(:,:)
    sb_loss                     = hlp10(:,:,1)
    sb_loss                     = sb_loss + 0.9_wp

    sb_pool                     = 4.4_wp
    sb_loss                     = 2.0_wp
    sb_loss                     = (5.0_wp + sb_loss) - (sb_pool / hlp1(:,:)) * sb_pool - 1.0_wp      ! 7 - 2.2 * 4.4 - 1 == 7 - 9.68 - 1 = -3.68

  END SUBROUTINE update_sb_pools

END MODULE mo_sb_update_pools
