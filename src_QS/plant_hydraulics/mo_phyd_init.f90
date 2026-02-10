!> QUINCY plant hydraulic variables init
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
!>#### initialization of plant hydraulic memory variables using, e.g., ic & bc input files
!>
MODULE mo_phyd_init
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract

  dsl4jsb_Use_processes PHYD_
  dsl4jsb_Use_config(PHYD_)
  dsl4jsb_Use_memory(PHYD_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: phyd_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_phyd_init'

CONTAINS

  ! ======================================================================================================= !
  !> run quincy plant hydraulics init
  !> 
  !> 
  SUBROUTINE phyd_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':phyd_init'

    model => Get_model(tile%owner_model_id)

    ! init boundary conditions
    CALL phyd_init_bc(tile)
  END SUBROUTINE phyd_init

  ! ======================================================================================================= !
  !> Intialize plant hydraulic variables for boundary conditions (bc)
  !>
  !>
  SUBROUTINE phyd_init_bc(tile)

    ! see also the USE statements in the module header
    USE mo_jsb_impl_constants,    ONLY: false, true
    
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),       POINTER  :: model
    
    CHARACTER(len=*), PARAMETER :: routine = modname//':phyd_init_bc'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(PHYD_)
    dsl4jsb_Def_memory(PHYD_)
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onDomain    :: psi_leaf
    dsl4jsb_Real2D_onDomain    :: psi_stem
    ! ---------------------------
    ! 0.4 Debug Option
    IF (debug_on()) CALL message(routine, 'Setting boundary conditions of phyd memory (quincy) for tile '// &
      &                          TRIM(tile%name))
    ! ---------------------------
    ! 0.5 Get Memory
    model => get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(PHYD_)
    ! Get process memories
    dsl4jsb_Get_memory(PHYD_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onDomain(PHYD_, psi_leaf)
    dsl4jsb_Get_var2D_onDomain(PHYD_, psi_stem)

    ! Initialse Stem and leaf water potentials with zero assuming super well saturated
    ! leaves and stem
    psi_leaf(:,:) = 0.0_wp
    psi_stem(:,:) = 0.0_wp

    IF (dsl4jsb_Config(PHYD_)%use_plant_hydraulics) THEN
      CALL message(TRIM(routine), ' Using plant hydraulics')
    ELSE
      CALL message(TRIM(routine), ' Not using plant hydraulics')
    ENDIF
  END SUBROUTINE phyd_init_bc

#endif
END MODULE mo_phyd_init
