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

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_jsb_process_class,     ONLY: SB_   ! A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
    
    ! Use of pools
    dsl4jsb_Use_pools(SB_) t_sb_components_pool, t_sb_elements_pool

    ! Use of process configurations and memory
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout)    :: tile          !< one tile with data structure for one lct

    TYPE(t_jsb_model),      POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib        !< land-cover-type library - parameter across pft's
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init'

    ! Declare process configuration, memory, and pool-structure Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_pool(t_sb_components_pool)  :: sb_pool
    dsl4jsb_Def_pool(t_sb_components_pool)  :: sb_formation
    dsl4jsb_Def_pool(t_sb_components_pool)  :: sb_loss
    dsl4jsb_Def_pool(t_sb_components_pool)  :: sb_transport
    dsl4jsb_Def_pool(t_sb_elements_pool)    :: sb_delta_pools
    dsl4jsb_Def_pool(t_sb_elements_pool)    :: sb_mycorrhiza_export

    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of vegetation memory for tile '//TRIM(tile%name))
    model  => Get_model(tile%owner_model_id)
    IF (.NOT. model%config%use_quincy) RETURN       ! only with QUINCY model
    IF (tile%lcts(1)%lib_id == 0) RETURN            ! run this init only if the present tile is a pft
    lctlib => model%lctlib(tile%lcts(1)%lib_id)

    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_pool,               sb_pool)
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_formation,          sb_formation)
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_loss,               sb_loss)
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_transport,          sb_transport)
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_delta_pools,        sb_delta_pools)
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_mycorrhiza_export,  sb_mycorrhiza_export)
    ! associate domain variables
    CALL sb_pool%Associate_domain_var_ptr()
    CALL sb_formation%Associate_domain_var_ptr()
    CALL sb_loss%Associate_domain_var_ptr()
    CALL sb_transport%Associate_domain_var_ptr()
    CALL sb_delta_pools%Associate_domain_var_ptr()
    CALL sb_mycorrhiza_export%Associate_domain_var_ptr()

    ! get var


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    ! testing pools
    hlp1 = 0.5_wp

    sb_pool       = hlp1
    sb_formation  = 3.5_wp
    sb_loss       = sb_pool + sb_formation
    
    sb_transport      = 1.2_wp
    sb_transport%dom  = 2.0_wp
    
    sb_delta_pools        = 10.0_wp
    sb_mycorrhiza_export  = sb_delta_pools + 5.0_wp
    
    ! print*, ">>>>>>>>>>>>>>>>> sb init <<<<<<<<<<<<<<<<<<<<<"
    ! print*, "0.5 0.5 & 1.0 " , sb_pool%dom%domain_c(1:2,1,1), " + ", sb_pool%dom%domain_c(2,1,1) + sb_pool%dom%domain_c(2,1,1)
    ! print*, "4.0 " , sb_loss%fungi%domain_n(1:2,1,1)
    ! print*, "15.0 " , sb_mycorrhiza_export%domain_p(1:2,1,1)
    ! print*, ">>>>>>>>>>>>>>>>> --      <<<<<<<<<<<<<<<<<<<<<"
    
  END SUBROUTINE sb_init

END MODULE mo_sb_init
