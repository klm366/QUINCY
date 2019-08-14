!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., and Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Initialization of the the assimilation memory.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20161213)
!!
MODULE mo_assimi_init

#include "dsl_jsbach.inc"

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: assimi_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize assimilation process (after memory has been set up)
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE assimi_init(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_assimi_constants,      ONLY: t_jmax_opt_min

    ! Use of process configurations (t_PROC_config)
    !
    ! USE PROC_memory
    dsl4jsb_Use_memory(ASSIMI_)


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib     !< land-cover-type library - parameter across pft's
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':assimi_init'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ASSIMI_)
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onDomain      :: beta_air
    dsl4jsb_Real2D_onDomain      :: beta_soa
    dsl4jsb_Real2D_onDomain      :: soa_tsoa_mavg
    dsl4jsb_Real2D_onDomain      :: beta_soa_tphen_mavg
    dsl4jsb_Real2D_onDomain      :: beta_air_daytime
    dsl4jsb_Real2D_onDomain      :: beta_air_daytime_dacc
    dsl4jsb_Real2D_onDomain      :: beta_air_tfrac_mavg
    dsl4jsb_Real2D_onDomain      :: beta_air_tcnl_mavg
    dsl4jsb_Real2D_onDomain      :: beta_soil_ps
    dsl4jsb_Real2D_onDomain      :: beta_soil_ps_daytime
    dsl4jsb_Real2D_onDomain      :: beta_soil_ps_daytime_dacc
    dsl4jsb_Real2D_onDomain      :: beta_soil_ps_tfrac_mavg
    dsl4jsb_Real2D_onDomain      :: beta_soil_ps_tcnl_mavg
    dsl4jsb_Real2D_onDomain      :: beta_soil_gs
    dsl4jsb_Real2D_onDomain      :: beta_soil_gs_daytime
    dsl4jsb_Real2D_onDomain      :: beta_soil_gs_daytime_dacc
    dsl4jsb_Real2D_onDomain      :: beta_soil_gs_tphen_mavg
    dsl4jsb_Real2D_onDomain      :: beta_soil_gs_tfrac_mavg
    dsl4jsb_Real2D_onDomain      :: beta_soil_gs_tcnl_mavg
    dsl4jsb_Real2D_onDomain      :: t_jmax_opt
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(ASSIMI_)) RETURN
    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! Get process config
    !
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_air)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soa)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, soa_tsoa_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soa_tphen_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_air_daytime)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_air_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_air_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_air_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_ps)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_ps_daytime)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_ps_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_ps_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_ps_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_gs)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_gs_daytime)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_gs_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_gs_tphen_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_gs_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_gs_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, t_jmax_opt)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    beta_air(:,:)                     = 1.0_wp
    beta_soa(:,:)                     = 1.0_wp
    soa_tsoa_mavg(:,:)                = -10.0_wp
    beta_soa_tphen_mavg(:,:)          = 1.0_wp
    beta_soil_ps(:,:)                 = 1.0_wp
    beta_soil_gs(:,:)                 = 1.0_wp
    
    beta_air_daytime(:,:)             = 1.0_wp
    beta_air_daytime_dacc(:,:)        = 1.0_wp
    beta_soil_ps_daytime(:,:)         = 1.0_wp
    beta_soil_ps_daytime_dacc(:,:)    = 1.0_wp
    beta_soil_gs_daytime(:,:)         = 1.0_wp
    beta_soil_gs_daytime_dacc(:,:)    = 1.0_wp
    
    beta_air_tfrac_mavg(:,:)          = 1.0_wp
    beta_air_tcnl_mavg(:,:)           = 1.0_wp
    beta_soil_ps_tfrac_mavg(:,:)      = 1.0_wp
    beta_soil_ps_tcnl_mavg(:,:)       = 1.0_wp
    beta_soil_gs_tphen_mavg(:,:)      = 1.0_wp
    beta_soil_gs_tfrac_mavg(:,:)      = 1.0_wp
    beta_soil_gs_tcnl_mavg(:,:)       = 1.0_wp

    t_jmax_opt(:,:)                   = t_jmax_opt_min 

  END SUBROUTINE assimi_init

END MODULE mo_assimi_init
