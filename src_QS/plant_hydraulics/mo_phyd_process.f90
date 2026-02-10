!> QUINCY calculate plant hydraulics
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### Update the plant hydraulic architectural states 
!>
MODULE mo_phyd_process
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_impl_constants,  ONLY: true, false, test_false_true
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_plant_hydraulics

  CHARACTER(len=*), PARAMETER :: modname = 'mo_phyd_process'

CONTAINS

  ! ======================================================================================================= !
  !> ...
  !>
  !>
  SUBROUTINE update_plant_hydraulics(tile, options)
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: A2L_, PHYD_, HYDRO_
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(PHYD_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER   :: model                  !< the model
    TYPE(t_lctlib_element), POINTER   :: lctlib                 !< land-cover-type library - parameter across pft's
    INTEGER                           :: iblk, ics, ice, nc     !< dimensions
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_plant_hydraulics'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(PHYD_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_ 2D
    dsl4jsb_Real2D_onChunk      :: t_air
    ! HYDRO_ 3D 
    dsl4jsb_Real3D_onChunk      :: wtr_soil_pot_sl               
    ! PHYD_2D 
    dsl4jsb_Real2D_onChunk      :: psi_leaf           
    dsl4jsb_Real2D_onChunk      :: psi_stem           
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(PHYD_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)                 RETURN  ! run task only if the present tile is a pft
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(PHYD_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_ 2D
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)                  ! in
    ! HYDRO_ 3D
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pot_sl)      ! in  [MPa]
    ! PHYD_2D
    dsl4jsb_Get_var2D_onChunk(PHYD_, psi_leaf)            ! inout
    dsl4jsb_Get_var2D_onChunk(PHYD_, psi_stem)            ! inout 
    ! ----------------------------------------------------------------------------------------------------- !

    !> 1.0 Calculate plant hydraulic states
    !>
    CALL calc_plant_hydraulics(              &
       & t_air(:),                           &
       & psi_leaf(:),                        &
       & psi_stem(:)                            )    
  END SUBROUTINE update_plant_hydraulics 


  ! ======================================================================================================= !
  !> calculate plant hydraulics 
  !> 
  ELEMENTAL SUBROUTINE calc_plant_hydraulics(   &
    & t_air,                                    &
    & psi_leaf,                                 &
    & psi_stem                                    )
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_math_constants,      ONLY: one_day
    USE mo_jsb_physical_constants,  ONLY: Tzero
    USE mo_jsb_impl_constants,      ONLY: true, false, test_false_true
    USE mo_phyd_constants,        ONLY: mpa_to_pa, pa_to_mpa, kg_h2o_to_mol
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp),     INTENT(in)    :: t_air                            !< air temp
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp),     INTENT(inout) :: psi_leaf                         !< leaf water potential
    REAL(wp),     INTENT(inout) :: psi_stem                         !< stem water potential
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_plant_hydraulics'
    ! ----------------------------------------------------------------------------------------------------- !

    ! Dummy code just to have some output

    psi_leaf = 1.1_wp
    psi_stem = t_air - 3.0_wp

  END SUBROUTINE calc_plant_hydraulics 

#endif
END MODULE mo_phyd_process
