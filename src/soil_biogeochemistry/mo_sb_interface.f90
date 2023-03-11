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
!! @brief Contains the interface to the soil_biogeochemistry processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20170104)
!! 2nd version by Jan Engel (20200618), add to jsbach4
!! code additions for bgcm store solution by Julia Nabel, MPI BGC (20230117)
!!
MODULE mo_sb_interface
#ifndef __NO_JSBACH__

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process

  ! "Integrate" routines
  USE mo_sb_rate_modifier,            ONLY: update_sb_rate_modifier
  USE mo_sb_simple_model,             ONLY: update_sb_simple_model
  ! USE mo_sb_jsm_main,                 ONLY: update_sb_jsm
  USE mo_sb_update_pools,             ONLY: update_sb_pools

  ! Use of processes in this module
  dsl4jsb_Use_processes SB_, VEG_
  ! Use of process configurations
  dsl4jsb_Use_config(SB_)
  ! Use of process memories
  dsl4jsb_Use_memory(SB_)
  dsl4jsb_Use_memory(VEG_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_sb_tasks

  ! ======================================================================================================= !
  !>Type definition: reset_sb_fluxes task
  !> 
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_sb_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_reset_sb_fluxes  !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_reset_sb_fluxes    !< Aggregates computed task variables
  END TYPE tsk_reset_sb_fluxes

  ! ======================================================================================================= !
  !>Type definition: sb_rate_modifier task
  !> 
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_sb_rate_modifier
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sb_rate_modifier       !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_sb_rate_modifier    !< Aggregates computed task variables
  END TYPE tsk_sb_rate_modifier
  
  ! ======================================================================================================= !
  !>Type definition: sb_simple_model task
  !> 
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_sb_simple_model
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sb_simple_model        !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_sb_simple_model     !< Aggregates computed task variables
  END TYPE tsk_sb_simple_model
  
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Type definition: update_sb_jsm task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_sb_jsm
  ! CONTAINS
  !   PROCEDURE, NOPASS :: Integrate => update_sb_jsm      !< Advances task computation for one timestep
  !   PROCEDURE, NOPASS :: Aggregate => aggr_update_sb_jsm !< Aggregates computed task variables
  ! END TYPE tsk_update_sb_jsm

  ! ======================================================================================================= !
  !>Type definition: sb_pools task
  !> 
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_sb_pools
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sb_pools         !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_sb_pools      !< Aggregates computed task variables
  END TYPE tsk_sb_pools  

  ! ======================================================================================================= !
  !>Type definition: time_average_soilbiogeochemistry task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_time_average_soilbiogeochemistry
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_time_average_soilbiogeochemistry    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_time_average_soilbiogeochemistry !< Aggregates computed task variables
  END TYPE tsk_time_average_soilbiogeochemistry
  
    ! ======================================================================================================= !
    !>Type definition: copy_pool_to_bgcm_test task 
    !>
    !> for testing only !!!
    !>
    TYPE, EXTENDS(t_jsb_process_task) :: tsk_copy_pool_to_bgcm_test
    CONTAINS
      PROCEDURE, NOPASS :: Integrate => update_copy_pool_to_bgcm_test    !< Advances task computation for one timestep
      PROCEDURE, NOPASS :: Aggregate => aggregate_copy_pool_to_bgcm_test !< Aggregates computed task variables
    END TYPE tsk_copy_pool_to_bgcm_test
    
    ! ======================================================================================================= !
    !>Type definition: sb_pool_test task 
    !>
    !> for testing only !!!
    !>
    TYPE, EXTENDS(t_jsb_process_task) :: tsk_sb_pool_test
    CONTAINS
      PROCEDURE, NOPASS :: Integrate => update_sb_pool_test    !< Advances task computation for one timestep
      PROCEDURE, NOPASS :: Aggregate => aggregate_sb_pool_test !< Aggregates computed task variables
    END TYPE tsk_sb_pool_test
    
    ! ======================================================================================================= !
    !>Type definition: sb_bgcm_pool_test task 
    !>
    !> for testing only !!!
    !>
    TYPE, EXTENDS(t_jsb_process_task) :: tsk_sb_bgcm_pool_test
    CONTAINS
      PROCEDURE, NOPASS :: Integrate => update_sb_bgcm_pool_test    !< Advances task computation for one timestep
      PROCEDURE, NOPASS :: Aggregate => aggregate_sb_bgcm_pool_test !< Aggregates computed task variables
    END TYPE tsk_sb_bgcm_pool_test

  ! ======================================================================================================= !
  !>Constructor interface: reset_sb_fluxes task
  !> 
  INTERFACE tsk_reset_sb_fluxes
    PROCEDURE Create_task_reset_sb_fluxes         !< Constructor function for task
  END INTERFACE tsk_reset_sb_fluxes

  ! ======================================================================================================= !
  !>Constructor interface: sb_rate_modifier task
  !> 
  INTERFACE tsk_sb_rate_modifier
    PROCEDURE Create_task_sb_rate_modifier         !< Constructor function for task
  END INTERFACE tsk_sb_rate_modifier
  
  ! ======================================================================================================= !
  !>Constructor interface: sb_simple_model task
  !> 
  INTERFACE tsk_sb_simple_model
    PROCEDURE Create_task_sb_simple_model         !< Constructor function for task
  END INTERFACE tsk_sb_simple_model
  
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor interface: update_sb_jsm task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! INTERFACE tsk_update_sb_jsm
  !   PROCEDURE Create_task_update_sb_jsm         !< Constructor function for task
  ! END INTERFACE tsk_update_sb_jsm

  ! ======================================================================================================= !
  !>Constructor interface: sb_pools task
  !> 
  INTERFACE tsk_sb_pools
    PROCEDURE Create_task_sb_pools         !< Constructor function for task
  END INTERFACE tsk_sb_pools

  ! ======================================================================================================= !
  !>Constructor interface: time_average_soilbiogeochemistry task
  !>
  INTERFACE tsk_time_average_soilbiogeochemistry
    PROCEDURE Create_task_update_time_average_soilbiogeochemistry         !< Constructor function for task
  END INTERFACE tsk_time_average_soilbiogeochemistry

    ! ======================================================================================================= !
    !>Constructor interface: copy_pool_to_bgcm_test task
    !> 
    INTERFACE tsk_copy_pool_to_bgcm_test
      PROCEDURE Create_task_copy_pool_to_bgcm_test         !< Constructor function for task
    END INTERFACE tsk_copy_pool_to_bgcm_test
    
    ! ======================================================================================================= !
    !>Constructor interface: sb_pool_test task
    !> 
    INTERFACE tsk_sb_pool_test
      PROCEDURE Create_task_sb_pool_test         !< Constructor function for task
    END INTERFACE tsk_sb_pool_test
    
    ! ======================================================================================================= !
    !>Constructor interface: sb_bgcm_pool_test task
    !> 
    INTERFACE tsk_sb_bgcm_pool_test
      PROCEDURE Create_task_sb_bgcm_pool_test         !< Constructor function for task
    END INTERFACE tsk_sb_bgcm_pool_test

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_interface'

CONTAINS

  ! ======================================================================================================= !
  !>Register tasks: SB_ process
  !> 
  SUBROUTINE Register_sb_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_reset_sb_fluxes(model_id))
    CALL this%Register_task(tsk_sb_rate_modifier(model_id))
    CALL this%Register_task(tsk_sb_simple_model(model_id))
    ! CALL this%Register_task(tsk_update_sb_jsm(model_id))
    CALL this%Register_task(tsk_sb_pools(model_id))
    CALL this%Register_task(tsk_time_average_soilbiogeochemistry(model_id))
    
    ! testing
    CALL this%Register_task(tsk_copy_pool_to_bgcm_test(model_id))
    CALL this%Register_task(tsk_sb_pool_test(model_id))
    CALL this%Register_task(tsk_sb_bgcm_pool_test(model_id))

  END SUBROUTINE Register_sb_tasks

  ! ======================================================================================================= !
  !>Constructor: reset_sb_fluxes task
  !>
  FUNCTION Create_task_reset_sb_fluxes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_reset_sb_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_sb_fluxes', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_sb_fluxes

  ! ======================================================================================================= !
  !>Constructor: sb_rate_modifier task
  !> 
  FUNCTION Create_task_sb_rate_modifier(model_id) RESULT(return_ptr)
  
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr
  
    ALLOCATE(tsk_sb_rate_modifier::return_ptr)
    CALL return_ptr%Construct(name='sb_rate_modifier', process_id=SB_, owner_model_id=model_id)
  
  END FUNCTION Create_task_sb_rate_modifier
  
  ! ======================================================================================================= !
  !>Constructor: sb_simple_model task
  !> 
  FUNCTION Create_task_sb_simple_model(model_id) RESULT(return_ptr)
  
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr
  
    ALLOCATE(tsk_sb_simple_model::return_ptr)
    CALL return_ptr%Construct(name='sb_simple_model', process_id=SB_, owner_model_id=model_id)
  
  END FUNCTION Create_task_sb_simple_model
  
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor: update_sb_jsm task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! FUNCTION Create_task_update_sb_jsm(model_id) RESULT(return_ptr)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   INTEGER,                   INTENT(in) :: model_id
  !   CLASS(t_jsb_process_task), POINTER    :: return_ptr
  ! 
  ! 
  !   ALLOCATE(tsk_update_sb_jsm::return_ptr)
  !   CALL return_ptr%Construct(name='update_sb_jsm', process_id=SB_, owner_model_id=model_id)
  ! 
  ! END FUNCTION Create_task_update_sb_jsm

  ! ======================================================================================================= !
  !>Constructor: sb_pools task
  !> 
  FUNCTION Create_task_sb_pools(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_sb_pools::return_ptr)
    CALL return_ptr%Construct(name='sb_pools', process_id=SB_, owner_model_id=model_id)
    
  END FUNCTION Create_task_sb_pools

  ! ======================================================================================================= !
  !>Constructor: tsk_time_average_soilbiogeochemistry task
  !> 
  FUNCTION Create_task_update_time_average_soilbiogeochemistry(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_time_average_soilbiogeochemistry::return_ptr)
    CALL return_ptr%Construct(name='tavrg_soilbiogeochemistry', process_id=SB_, owner_model_id=model_id)
    
  END FUNCTION Create_task_update_time_average_soilbiogeochemistry
  
    ! ======================================================================================================= !
    !>Constructor: tsk_copy_pool_to_bgcm_test task
    !> 
    FUNCTION Create_task_copy_pool_to_bgcm_test(model_id) RESULT(return_ptr)
      INTEGER,                   INTENT(in) :: model_id
      CLASS(t_jsb_process_task), POINTER    :: return_ptr
  
      ALLOCATE(tsk_copy_pool_to_bgcm_test::return_ptr)
      CALL return_ptr%Construct(name='copy_pool_to_bgcm_test', process_id=SB_, owner_model_id=model_id)
    END FUNCTION Create_task_copy_pool_to_bgcm_test
    
    ! ======================================================================================================= !
    !>Constructor: tsk_sb_pool_test task
    !> 
    FUNCTION Create_task_sb_pool_test(model_id) RESULT(return_ptr)
      INTEGER,                   INTENT(in) :: model_id
      CLASS(t_jsb_process_task), POINTER    :: return_ptr
  
      ALLOCATE(tsk_sb_pool_test::return_ptr)
      CALL return_ptr%Construct(name='sb_pool_test', process_id=SB_, owner_model_id=model_id)
    END FUNCTION Create_task_sb_pool_test
    
    ! ======================================================================================================= !
    !>Constructor: tsk_sb_bgcm_pool_test task
    !> 
    FUNCTION Create_task_sb_bgcm_pool_test(model_id) RESULT(return_ptr)
      INTEGER,                   INTENT(in) :: model_id
      CLASS(t_jsb_process_task), POINTER    :: return_ptr
  
      ALLOCATE(tsk_sb_bgcm_pool_test::return_ptr)
      CALL return_ptr%Construct(name='sb_bgcm_pool_test', process_id=SB_, owner_model_id=model_id)
    END FUNCTION Create_task_sb_bgcm_pool_test

  ! ------------------------------------------------------------------------------------------------------- !
  ! update_* routines
  ! ------------------------------------------------------------------------------------------------------- !

  ! ======================================================================================================= !
  !>Implementation of "integrate": reset_sb_fluxes task
  !>
  SUBROUTINE update_reset_sb_fluxes(tile, options)

    USE mo_sb_util, ONLY: real_reset_sb_fluxes => reset_sb_fluxes
    dsl4jsb_Use_processes SB_

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_reset_sb_fluxes'

    iblk = options%iblk
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL real_reset_sb_fluxes(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_reset_sb_fluxes

  ! ======================================================================================================= !
  !>calculate time moving averages and daytime averages for SB_
  !>
  SUBROUTINE update_time_average_soilbiogeochemistry(tile, options)

    USE mo_jsb_control,         ONLY: debug_on
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,      ONLY: t_jsb_task_options
    USE mo_sb_util,             ONLY: calculate_time_average_soilbiogeochemistry
    dsl4jsb_Use_processes SB_
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_time_average_soilbiogeochemistry'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)          RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. options%iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL calculate_time_average_soilbiogeochemistry(tile, options)
    
    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_time_average_soilbiogeochemistry

  ! ======================================================================================================= !
  !> test that bgcm implementation and application is working
  !> copy all values from sb_pool to sb_bgcm_pool
  !>
  SUBROUTINE update_copy_pool_to_bgcm_test(tile, options)
    USE mo_kind,                  ONLY: wp
    USE mo_jsb_control,           ONLY: debug_on
    USE mo_exception,             ONLY: message, message_text, finish
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid, t_jsb_grid
    USE mo_jsb_grid,              ONLY: Get_vgrid, Get_grid
    USE mo_lnd_bgcm_class,        ONLY: ELEM_C_ID, ELEM_N_ID, ELEM_P_ID, ELEM_C13_ID, ELEM_C14_ID, ELEM_N15_ID, &
      &                                 LAST_ELEM_ID
    USE mo_lnd_bgcm_store_class,  ONLY: SB_BGCM_POOL_ID, VEG_BGCM_POOL_ID
    USE mo_lnd_bgcm_store,        ONLY: t_lnd_bgcm_store
    USE mo_veg_config_class,      ONLY: get_number_of_veg_compartments, &
      &                                 VEG_PART_LEAF_IDX, VEG_PART_FINE_ROOT_IDX, VEG_PART_COARSE_ROOT_IDX, &
      &                                 VEG_PART_SAP_WOOD_IDX, VEG_PART_HEART_WOOD_IDX, VEG_PART_LABILE_IDX, &
      &                                 VEG_PART_RESERVE_IDX, VEG_PART_FRUIT_IDX, VEG_PART_SEED_BED_IDX
    USE mo_sb_config_class,       ONLY: get_number_of_sb_compartments, &
      &                                 SB_PART_DOM_IDX, SB_PART_DOM_ASSOC_IDX, SB_PART_SOLUBLE_LITTER_IDX, &
      &                                 SB_PART_POLYMERIC_LITTER_IDX, SB_PART_WOODY_LITTER_IDX, SB_PART_FUNGI_IDX, &
      &                                 SB_PART_MYCORRHIZA_IDX, SB_PART_MICROBIAL_IDX, SB_PART_RESIDUE_IDX, &
      &                                 SB_PART_RESIDUE_ASSOC_IDX
    dsl4jsb_Use_processes SB_, VEG_
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_pools(SB_)  t_sb_components_pool, t_sb_bgcm_with_components, t_sb_bgcm_with_elements
    dsl4jsb_Use_pools(VEG_) t_veg_tissues_pool
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_jsb_grid),       POINTER       :: hgrid                !< Horizontal grid
    TYPE(t_jsb_vgrid),      POINTER       :: vgrid_soil_sb        !< Vertical grid
    TYPE(t_lnd_bgcm_store), POINTER       :: bgcm_store           !< the bgcm store of this tile
    INTEGER                               :: nsoil_sb             !< number of soil layers as used/defined by the SB_ process
    INTEGER                               :: elem_idx_map(LAST_ELEM_ID - 1)  !< element mapper ID -> IND
    INTEGER                               :: isoil                !< loop over nsoil_sb
    INTEGER                               :: iblk, ics, ice, nc   !< dimensions block/chunk
    INTEGER                               :: nr_elements          !< number of elements
    INTEGER                               :: nr_sb_parts          !< number of sb compartments
    INTEGER                               :: nr_veg_parts         !< number of veg compartments
    REAL(wp), ALLOCATABLE                 :: sb_bgcm_pool_mt(:,:,:,:)
    REAL(wp), ALLOCATABLE                 :: veg_bgcm_pool_mt(:,:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_copy_pool_to_bgcm_test'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(VEG_)
    !------------------------------------------------------------------------------------------------------ !
    dsl4jsb_Def_pool(t_sb_bgcm_with_components)  :: sb_bgcm_pool
    dsl4jsb_Def_pool(t_sb_bgcm_with_components)  :: sb_bgcm_formation
    dsl4jsb_Def_pool(t_sb_bgcm_with_elements)    :: sb_bgcm_delta_pools
    dsl4jsb_Def_pool(t_sb_components_pool)       :: sb_pool
    dsl4jsb_Def_pool(t_veg_tissues_pool)         :: veg_pool
    !------------------------------------------------------------------------------------------------------ !
    ! ! SB_ 2D
    ! dsl4jsb_Real2D_onChunk :: leaching_nh4_solute
    ! dsl4jsb_Real2D_onChunk :: leaching_nh4_n15_solute
    ! dsl4jsb_Real2D_onChunk :: leaching_no3_solute
    ! dsl4jsb_Real2D_onChunk :: leaching_po4_solute
    ! dsl4jsb_Real2D_onChunk :: leaching_no3_n15_solute
    ! ! SB_ 3D
    ! dsl4jsb_Real3D_onChunk :: noy
    ! dsl4jsb_Real3D_onChunk :: n2o
    ! dsl4jsb_Real3D_onChunk :: n2
    !------------------------------------------------------------------------------------------------------ !
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    !------------------------------------------------------------------------------------------------------ !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)          RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. options%iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    !------------------------------------------------------------------------------------------------------ !
    model         => Get_model(tile%owner_model_id)
    hgrid         => Get_grid(model%grid_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    !------------------------------------------------------------------------------------------------------ !
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(VEG_)
    !------------------------------------------------------------------------------------------------------ !
    dsl4jsb_Get_pool(SB_, sb_bgcm_main%sb_bgcm_pool,        sb_bgcm_pool)
    dsl4jsb_Get_pool(SB_, sb_bgcm_main%sb_bgcm_formation,   sb_bgcm_formation)
    dsl4jsb_Get_pool(SB_, sb_bgcm_main%sb_bgcm_delta_pools, sb_bgcm_delta_pools)
    dsl4jsb_Get_pool(SB_,  sb_main_pool%sb_pool,            sb_pool)
    dsl4jsb_Get_pool(VEG_, veg_main_pool%veg_pool,          veg_pool)
    CALL sb_pool%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_pool%Associate_chunk_var_ptr(ics, ice, iblk)
    !------------------------------------------------------------------------------------------------------ !
    elem_idx_map(:) = model%config%elements_index_map(:)
    nr_elements     = model%config%nr_of_used_elements
    nr_sb_parts     = get_number_of_sb_compartments()
    nr_veg_parts    = get_number_of_veg_compartments()
    !------------------------------------------------------------------------------------------------------ !
    ALLOCATE(sb_bgcm_pool_mt(nr_sb_parts, nr_elements, nc, nsoil_sb))
    ALLOCATE(veg_bgcm_pool_mt(nr_veg_parts, nr_elements, nc))

    bgcm_store => tile%bgcm_store
    CALL bgcm_store%Get_matrix_from_bgcm_2l(SB_BGCM_POOL_ID, ics, ice, iblk, tile%name, sb_bgcm_pool_mt)
    CALL bgcm_store%Get_matrix_from_bgcm_2l(VEG_BGCM_POOL_ID, ics, ice, iblk, tile%name, veg_bgcm_pool_mt)
    !------------------------------------------------------------------------------------------------------ !
    ! ! SB_ 2D
    ! dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_solute) 
    ! dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_n15_solute) 
    ! dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_solute) 
    ! dsl4jsb_Get_var2D_onChunk(SB_,       leaching_po4_solute) 
    ! dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_n15_solute) 
    ! ! SB_ 3D
    ! dsl4jsb_Get_var3D_onChunk(SB_,       noy) 
    ! dsl4jsb_Get_var3D_onChunk(SB_,       n2o) 
    ! dsl4jsb_Get_var3D_onChunk(SB_,       n2) 

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------------------------ !
    !> showcase for bgcm
    !> simple test using the ptr3d to show that the bgcm are working
    !>
    !   does not work with p%ptrm1(:,:) because it is not a component of 't_jsb_var'
    !   but only of 'TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real3d'
    !   (would need a 'SELECT TYPE' to work with ptrm1(:,:))
    ! sb_bgcm_pool%pool_list(SB_PART_DOM_IDX)%p%element_list(elem_idx_map(ELEM_C_ID))%p%ptr3d(ics:ice, :, iblk)  = 20.6_wp
    ! sb_bgcm_pool%pool_list(SB_PART_DOM_IDX)%p%element_list(elem_idx_map(ELEM_N_ID))%p%ptr3d(ics:ice, :, iblk)  = &
    !   & sb_bgcm_pool%pool_list(SB_PART_DOM_IDX)%p%element_list(elem_idx_map(ELEM_C_ID))%p%ptr3d(ics:ice, :, iblk) / 2.0_wp   ! = 10.3_wp
    ! sb_bgcm_pool%pool_list(SB_PART_DOM_IDX)%p%element_list(elem_idx_map(ELEM_P_ID))%p%ptr3d(ics:ice, :, iblk)  = n2(:,:) * 2.0_wp
    ! printing for debug
    ! print*, " --------------- "
    ! print*, TRIM(sb_bgcm_pool%pool_list(SB_PART_DOM_IDX)%p%shortname), " ", &
    !   & TRIM(sb_bgcm_pool%pool_list(SB_PART_DOM_IDX)%p%element_list(elem_idx_map(ELEM_C_ID))%p%name), " ", &
    !   & TRIM(sb_bgcm_pool%pool_list(SB_PART_DOM_IDX)%p%element_list(elem_idx_map(ELEM_N_ID))%p%name), " ", &
    !   & TRIM(sb_bgcm_pool%pool_list(SB_PART_DOM_IDX)%p%element_list(elem_idx_map(ELEM_P_ID))%p%name), " n2"
    ! print*, sb_bgcm_pool%pool_list(SB_PART_DOM_IDX)%p%element_list(elem_idx_map(ELEM_C_ID))%p%ptr3d(ics:(ics + 4), 1, iblk)
    ! print*, sb_bgcm_pool%pool_list(SB_PART_DOM_IDX)%p%element_list(elem_idx_map(ELEM_N_ID))%p%ptr3d(ics:(ics + 4), 1, iblk)
    ! print*, sb_bgcm_pool%pool_list(SB_PART_DOM_IDX)%p%element_list(elem_idx_map(ELEM_P_ID))%p%ptr3d(ics:(ics + 4), 1, iblk)
    ! print*, n2(1:4, 1)

    !------------------------------------------------------------------------------------------------------ !
    !> copy data
    !> copy data from sb_pool  -> sb_bgcm_pool
    !> copy data from veg_pool -> veg_bgcm_pool
    !>
    !  from:  TYPE, EXTENDS(t_jsb_pool_elem3d) :: t_sb_components_pool         | sb_pool
    !  to:    TYPE, EXTENDS(t_jsb_pool)        :: t_sb_bgcm_with_components | sb_bgcm_pool
    !  both these pools are used for additional testing in the tasks:
    !   sb_pool_test
    !   sb_bgcm_pool_test
    !  
    !  same for VEG_ veg_pool -> veg_bgcm_pool
    !
    ! SB_
    sb_bgcm_pool_mt(SB_PART_DOM_IDX, elem_idx_map(ELEM_C_ID), :, :)           = sb_pool%dom%carbon(:,:)
    sb_bgcm_pool_mt(SB_PART_DOM_IDX, elem_idx_map(ELEM_N_ID), :, :)           = sb_pool%dom%nitrogen(:,:)
    sb_bgcm_pool_mt(SB_PART_DOM_IDX, elem_idx_map(ELEM_P_ID), :, :)           = sb_pool%dom%phosphorus(:,:)
    sb_bgcm_pool_mt(SB_PART_DOM_IDX, elem_idx_map(ELEM_C13_ID), :, :)         = sb_pool%dom%carbon13(:,:)
    sb_bgcm_pool_mt(SB_PART_DOM_IDX, elem_idx_map(ELEM_C14_ID), :, :)         = sb_pool%dom%carbon14(:,:)
    sb_bgcm_pool_mt(SB_PART_DOM_IDX, elem_idx_map(ELEM_N15_ID), :, :)         = sb_pool%dom%nitrogen15(:,:)
    sb_bgcm_pool_mt(SB_PART_DOM_ASSOC_IDX, elem_idx_map(ELEM_C_ID), :, :)          = sb_pool%dom_assoc%carbon(:,:)
    sb_bgcm_pool_mt(SB_PART_DOM_ASSOC_IDX, elem_idx_map(ELEM_N_ID), :, :)          = sb_pool%dom_assoc%nitrogen(:,:)
    sb_bgcm_pool_mt(SB_PART_DOM_ASSOC_IDX, elem_idx_map(ELEM_P_ID), :, :)          = sb_pool%dom_assoc%phosphorus(:,:)
    sb_bgcm_pool_mt(SB_PART_DOM_ASSOC_IDX, elem_idx_map(ELEM_C13_ID), :, :)        = sb_pool%dom_assoc%carbon13(:,:)
    sb_bgcm_pool_mt(SB_PART_DOM_ASSOC_IDX, elem_idx_map(ELEM_C14_ID), :, :)        = sb_pool%dom_assoc%carbon14(:,:)
    sb_bgcm_pool_mt(SB_PART_DOM_ASSOC_IDX, elem_idx_map(ELEM_N15_ID), :, :)        = sb_pool%dom_assoc%nitrogen15(:,:)
    sb_bgcm_pool_mt(SB_PART_SOLUBLE_LITTER_IDX, elem_idx_map(ELEM_C_ID), :, :)     = sb_pool%soluable_litter%carbon(:,:)
    sb_bgcm_pool_mt(SB_PART_SOLUBLE_LITTER_IDX, elem_idx_map(ELEM_N_ID), :, :)     = sb_pool%soluable_litter%nitrogen(:,:)
    sb_bgcm_pool_mt(SB_PART_SOLUBLE_LITTER_IDX, elem_idx_map(ELEM_P_ID), :, :)     = sb_pool%soluable_litter%phosphorus(:,:)
    sb_bgcm_pool_mt(SB_PART_SOLUBLE_LITTER_IDX, elem_idx_map(ELEM_C13_ID), :, :)   = sb_pool%soluable_litter%carbon13(:,:)
    sb_bgcm_pool_mt(SB_PART_SOLUBLE_LITTER_IDX, elem_idx_map(ELEM_C14_ID), :, :)   = sb_pool%soluable_litter%carbon14(:,:)
    sb_bgcm_pool_mt(SB_PART_SOLUBLE_LITTER_IDX, elem_idx_map(ELEM_N15_ID), :, :)   = sb_pool%soluable_litter%nitrogen15(:,:)
    sb_bgcm_pool_mt(SB_PART_POLYMERIC_LITTER_IDX, elem_idx_map(ELEM_C_ID), :, :)   = sb_pool%polymeric_litter%carbon(:,:)
    sb_bgcm_pool_mt(SB_PART_POLYMERIC_LITTER_IDX, elem_idx_map(ELEM_N_ID), :, :)   = sb_pool%polymeric_litter%nitrogen(:,:)
    sb_bgcm_pool_mt(SB_PART_POLYMERIC_LITTER_IDX, elem_idx_map(ELEM_P_ID), :, :)   = sb_pool%polymeric_litter%phosphorus(:,:)
    sb_bgcm_pool_mt(SB_PART_POLYMERIC_LITTER_IDX, elem_idx_map(ELEM_C13_ID), :, :) = sb_pool%polymeric_litter%carbon13(:,:)
    sb_bgcm_pool_mt(SB_PART_POLYMERIC_LITTER_IDX, elem_idx_map(ELEM_C14_ID), :, :) = sb_pool%polymeric_litter%carbon14(:,:)
    sb_bgcm_pool_mt(SB_PART_POLYMERIC_LITTER_IDX, elem_idx_map(ELEM_N15_ID), :, :) = sb_pool%polymeric_litter%nitrogen15(:,:)
    sb_bgcm_pool_mt(SB_PART_WOODY_LITTER_IDX, elem_idx_map(ELEM_C_ID), :, :)    = sb_pool%woody_litter%carbon(:,:)
    sb_bgcm_pool_mt(SB_PART_WOODY_LITTER_IDX, elem_idx_map(ELEM_N_ID), :, :)    = sb_pool%woody_litter%nitrogen(:,:)
    sb_bgcm_pool_mt(SB_PART_WOODY_LITTER_IDX, elem_idx_map(ELEM_P_ID), :, :)    = sb_pool%woody_litter%phosphorus(:,:)
    sb_bgcm_pool_mt(SB_PART_WOODY_LITTER_IDX, elem_idx_map(ELEM_C13_ID), :, :)  = sb_pool%woody_litter%carbon13(:,:)
    sb_bgcm_pool_mt(SB_PART_WOODY_LITTER_IDX, elem_idx_map(ELEM_C14_ID), :, :)  = sb_pool%woody_litter%carbon14(:,:)
    sb_bgcm_pool_mt(SB_PART_WOODY_LITTER_IDX, elem_idx_map(ELEM_N15_ID), :, :)  = sb_pool%woody_litter%nitrogen15(:,:)
    sb_bgcm_pool_mt(SB_PART_FUNGI_IDX, elem_idx_map(ELEM_C_ID), :, :)           = sb_pool%fungi%carbon(:,:)
    sb_bgcm_pool_mt(SB_PART_FUNGI_IDX, elem_idx_map(ELEM_N_ID), :, :)           = sb_pool%fungi%nitrogen(:,:)
    sb_bgcm_pool_mt(SB_PART_FUNGI_IDX, elem_idx_map(ELEM_P_ID), :, :)           = sb_pool%fungi%phosphorus(:,:)
    sb_bgcm_pool_mt(SB_PART_FUNGI_IDX, elem_idx_map(ELEM_C13_ID), :, :)         = sb_pool%fungi%carbon13(:,:)
    sb_bgcm_pool_mt(SB_PART_FUNGI_IDX, elem_idx_map(ELEM_C14_ID), :, :)         = sb_pool%fungi%carbon14(:,:)
    sb_bgcm_pool_mt(SB_PART_FUNGI_IDX, elem_idx_map(ELEM_N15_ID), :, :)         = sb_pool%fungi%nitrogen15(:,:)
    sb_bgcm_pool_mt(SB_PART_MYCORRHIZA_IDX, elem_idx_map(ELEM_C_ID), :, :)      = sb_pool%mycorrhiza%carbon(:,:)
    sb_bgcm_pool_mt(SB_PART_MYCORRHIZA_IDX, elem_idx_map(ELEM_N_ID), :, :)      = sb_pool%mycorrhiza%nitrogen(:,:)
    sb_bgcm_pool_mt(SB_PART_MYCORRHIZA_IDX, elem_idx_map(ELEM_P_ID), :, :)      = sb_pool%mycorrhiza%phosphorus(:,:)
    sb_bgcm_pool_mt(SB_PART_MYCORRHIZA_IDX, elem_idx_map(ELEM_C13_ID), :, :)    = sb_pool%mycorrhiza%carbon13(:,:)
    sb_bgcm_pool_mt(SB_PART_MYCORRHIZA_IDX, elem_idx_map(ELEM_C14_ID), :, :)    = sb_pool%mycorrhiza%carbon14(:,:)
    sb_bgcm_pool_mt(SB_PART_MYCORRHIZA_IDX, elem_idx_map(ELEM_N15_ID), :, :)    = sb_pool%mycorrhiza%nitrogen15(:,:)
    sb_bgcm_pool_mt(SB_PART_MICROBIAL_IDX, elem_idx_map(ELEM_C_ID), :, :)       = sb_pool%microbial%carbon(:,:)
    sb_bgcm_pool_mt(SB_PART_MICROBIAL_IDX, elem_idx_map(ELEM_N_ID), :, :)       = sb_pool%microbial%nitrogen(:,:)
    sb_bgcm_pool_mt(SB_PART_MICROBIAL_IDX, elem_idx_map(ELEM_P_ID), :, :)       = sb_pool%microbial%phosphorus(:,:)
    sb_bgcm_pool_mt(SB_PART_MICROBIAL_IDX, elem_idx_map(ELEM_C13_ID), :, :)     = sb_pool%microbial%carbon13(:,:)
    sb_bgcm_pool_mt(SB_PART_MICROBIAL_IDX, elem_idx_map(ELEM_C14_ID), :, :)     = sb_pool%microbial%carbon14(:,:)
    sb_bgcm_pool_mt(SB_PART_MICROBIAL_IDX, elem_idx_map(ELEM_N15_ID), :, :)     = sb_pool%microbial%nitrogen15(:,:)
    sb_bgcm_pool_mt(SB_PART_RESIDUE_IDX, elem_idx_map(ELEM_C_ID), :, :)         = sb_pool%residue%carbon(:,:)
    sb_bgcm_pool_mt(SB_PART_RESIDUE_IDX, elem_idx_map(ELEM_N_ID), :, :)         = sb_pool%residue%nitrogen(:,:)
    sb_bgcm_pool_mt(SB_PART_RESIDUE_IDX, elem_idx_map(ELEM_P_ID), :, :)         = sb_pool%residue%phosphorus(:,:)
    sb_bgcm_pool_mt(SB_PART_RESIDUE_IDX, elem_idx_map(ELEM_C13_ID), :, :)       = sb_pool%residue%carbon13(:,:)
    sb_bgcm_pool_mt(SB_PART_RESIDUE_IDX, elem_idx_map(ELEM_C14_ID), :, :)       = sb_pool%residue%carbon14(:,:)
    sb_bgcm_pool_mt(SB_PART_RESIDUE_IDX, elem_idx_map(ELEM_N15_ID), :, :)       = sb_pool%residue%nitrogen15(:,:)
    sb_bgcm_pool_mt(SB_PART_RESIDUE_ASSOC_IDX, elem_idx_map(ELEM_C_ID), :, :)   = sb_pool%residue_assoc%carbon(:,:)
    sb_bgcm_pool_mt(SB_PART_RESIDUE_ASSOC_IDX, elem_idx_map(ELEM_N_ID), :, :)   = sb_pool%residue_assoc%nitrogen(:,:)
    sb_bgcm_pool_mt(SB_PART_RESIDUE_ASSOC_IDX, elem_idx_map(ELEM_P_ID), :, :)   = sb_pool%residue_assoc%phosphorus(:,:)
    sb_bgcm_pool_mt(SB_PART_RESIDUE_ASSOC_IDX, elem_idx_map(ELEM_C13_ID), :, :) = sb_pool%residue_assoc%carbon13(:,:)
    sb_bgcm_pool_mt(SB_PART_RESIDUE_ASSOC_IDX, elem_idx_map(ELEM_C14_ID), :, :) = sb_pool%residue_assoc%carbon14(:,:)
    sb_bgcm_pool_mt(SB_PART_RESIDUE_ASSOC_IDX, elem_idx_map(ELEM_N15_ID), :, :) = sb_pool%residue_assoc%nitrogen15(:,:)
    ! VEG_
    veg_bgcm_pool_mt(VEG_PART_LEAF_IDX, elem_idx_map(ELEM_C_ID), :)          = veg_pool%leaf%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_LEAF_IDX, elem_idx_map(ELEM_N_ID), :)          = veg_pool%leaf%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_LEAF_IDX, elem_idx_map(ELEM_P_ID), :)          = veg_pool%leaf%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_LEAF_IDX, elem_idx_map(ELEM_C13_ID), :)        = veg_pool%leaf%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_LEAF_IDX, elem_idx_map(ELEM_C14_ID), :)        = veg_pool%leaf%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_LEAF_IDX, elem_idx_map(ELEM_N15_ID), :)        = veg_pool%leaf%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_FINE_ROOT_IDX, elem_idx_map(ELEM_C_ID), :)     = veg_pool%fine_root%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_FINE_ROOT_IDX, elem_idx_map(ELEM_N_ID), :)     = veg_pool%fine_root%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_FINE_ROOT_IDX, elem_idx_map(ELEM_P_ID), :)     = veg_pool%fine_root%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_FINE_ROOT_IDX, elem_idx_map(ELEM_C13_ID), :)   = veg_pool%fine_root%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_FINE_ROOT_IDX, elem_idx_map(ELEM_C14_ID), :)   = veg_pool%fine_root%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_FINE_ROOT_IDX, elem_idx_map(ELEM_N15_ID), :)   = veg_pool%fine_root%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_COARSE_ROOT_IDX, elem_idx_map(ELEM_C_ID), :)   = veg_pool%coarse_root%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_COARSE_ROOT_IDX, elem_idx_map(ELEM_N_ID), :)   = veg_pool%coarse_root%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_COARSE_ROOT_IDX, elem_idx_map(ELEM_P_ID), :)   = veg_pool%coarse_root%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_COARSE_ROOT_IDX, elem_idx_map(ELEM_C13_ID), :) = veg_pool%coarse_root%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_COARSE_ROOT_IDX, elem_idx_map(ELEM_C14_ID), :) = veg_pool%coarse_root%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_COARSE_ROOT_IDX, elem_idx_map(ELEM_N15_ID), :) = veg_pool%coarse_root%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_SAP_WOOD_IDX, elem_idx_map(ELEM_C_ID), :)      = veg_pool%sap_wood%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_SAP_WOOD_IDX, elem_idx_map(ELEM_N_ID), :)      = veg_pool%sap_wood%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_SAP_WOOD_IDX, elem_idx_map(ELEM_P_ID), :)      = veg_pool%sap_wood%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_SAP_WOOD_IDX, elem_idx_map(ELEM_C13_ID), :)    = veg_pool%sap_wood%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_SAP_WOOD_IDX, elem_idx_map(ELEM_C14_ID), :)    = veg_pool%sap_wood%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_SAP_WOOD_IDX, elem_idx_map(ELEM_N15_ID), :)    = veg_pool%sap_wood%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_HEART_WOOD_IDX, elem_idx_map(ELEM_C_ID), :)    = veg_pool%heart_wood%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_HEART_WOOD_IDX, elem_idx_map(ELEM_N_ID), :)    = veg_pool%heart_wood%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_HEART_WOOD_IDX, elem_idx_map(ELEM_P_ID), :)    = veg_pool%heart_wood%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_HEART_WOOD_IDX, elem_idx_map(ELEM_C13_ID), :)  = veg_pool%heart_wood%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_HEART_WOOD_IDX, elem_idx_map(ELEM_C14_ID), :)  = veg_pool%heart_wood%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_HEART_WOOD_IDX, elem_idx_map(ELEM_N15_ID), :)  = veg_pool%heart_wood%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_LABILE_IDX, elem_idx_map(ELEM_C_ID), :)        = veg_pool%labile%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_LABILE_IDX, elem_idx_map(ELEM_N_ID), :)        = veg_pool%labile%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_LABILE_IDX, elem_idx_map(ELEM_P_ID), :)        = veg_pool%labile%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_LABILE_IDX, elem_idx_map(ELEM_C13_ID), :)      = veg_pool%labile%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_LABILE_IDX, elem_idx_map(ELEM_C14_ID), :)      = veg_pool%labile%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_LABILE_IDX, elem_idx_map(ELEM_N15_ID), :)      = veg_pool%labile%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_RESERVE_IDX, elem_idx_map(ELEM_C_ID), :)       = veg_pool%reserve%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_RESERVE_IDX, elem_idx_map(ELEM_N_ID), :)       = veg_pool%reserve%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_RESERVE_IDX, elem_idx_map(ELEM_P_ID), :)       = veg_pool%reserve%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_RESERVE_IDX, elem_idx_map(ELEM_C13_ID), :)     = veg_pool%reserve%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_RESERVE_IDX, elem_idx_map(ELEM_C14_ID), :)     = veg_pool%reserve%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_RESERVE_IDX, elem_idx_map(ELEM_N15_ID), :)     = veg_pool%reserve%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_FRUIT_IDX, elem_idx_map(ELEM_C_ID), :)         = veg_pool%fruit%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_FRUIT_IDX, elem_idx_map(ELEM_N_ID), :)         = veg_pool%fruit%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_FRUIT_IDX, elem_idx_map(ELEM_P_ID), :)         = veg_pool%fruit%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_FRUIT_IDX, elem_idx_map(ELEM_C13_ID), :)       = veg_pool%fruit%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_FRUIT_IDX, elem_idx_map(ELEM_C14_ID), :)       = veg_pool%fruit%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_FRUIT_IDX, elem_idx_map(ELEM_N15_ID), :)       = veg_pool%fruit%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_SEED_BED_IDX, elem_idx_map(ELEM_C_ID), :)      = veg_pool%seed_bed%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_SEED_BED_IDX, elem_idx_map(ELEM_N_ID), :)      = veg_pool%seed_bed%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_SEED_BED_IDX, elem_idx_map(ELEM_P_ID), :)      = veg_pool%seed_bed%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_SEED_BED_IDX, elem_idx_map(ELEM_C13_ID), :)    = veg_pool%seed_bed%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_SEED_BED_IDX, elem_idx_map(ELEM_C14_ID), :)    = veg_pool%seed_bed%carbon(:)
    veg_bgcm_pool_mt(VEG_PART_SEED_BED_IDX, elem_idx_map(ELEM_N15_ID), :)    = veg_pool%seed_bed%carbon(:)

    ! write back data from bgcm_matrix to bgc_material, and de-allocate local bgcm_matrix
    CALL bgcm_store%Write_matrix_to_bgcm_2l(SB_BGCM_POOL_ID, ics, ice, iblk, tile%name, sb_bgcm_pool_mt)
    CALL bgcm_store%Write_matrix_to_bgcm_2l(VEG_BGCM_POOL_ID, ics, ice, iblk, tile%name, veg_bgcm_pool_mt)
    DEALLOCATE(sb_bgcm_pool_mt, veg_bgcm_pool_mt)

    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_copy_pool_to_bgcm_test

  ! ======================================================================================================= !
  !> test runtime of bgc_material structures using the sb_bgcm_pool
  !> runtime comparison with update_sb_pool_test
  !>
  SUBROUTINE update_sb_pool_test(tile, options)
    USE mo_kind,                  ONLY: wp
    USE mo_jsb_control,           ONLY: debug_on
    USE mo_exception,             ONLY: message, message_text, finish
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid, t_jsb_grid
    USE mo_jsb_grid,              ONLY: Get_vgrid, Get_grid
    dsl4jsb_Use_processes SB_, VEG_
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_pools(SB_)  t_sb_components_pool
    dsl4jsb_Use_pools(VEG_) t_veg_tissues_pool
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_jsb_grid),       POINTER       :: hgrid                !< Horizontal grid
    TYPE(t_jsb_vgrid),      POINTER       :: vgrid_soil_sb        !< Vertical grid
    INTEGER                               :: nsoil_sb             !< number of soil layers as used/defined by the SB_ process
    INTEGER                               :: isoil                !< loop over nsoil_sb
    INTEGER                               :: iblk, ics, ice, nc   !< dimensions block/chunk
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_sb_pool_test'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(VEG_)
    !------------------------------------------------------------------------------------------------------ !
    dsl4jsb_Def_pool(t_sb_components_pool)  :: sb_pool
    dsl4jsb_Def_pool(t_veg_tissues_pool)    :: veg_pool
    !------------------------------------------------------------------------------------------------------ !
    ! ! SB_ 2D
    ! dsl4jsb_Real2D_onChunk :: leaching_nh4_solute
    ! dsl4jsb_Real2D_onChunk :: leaching_no3_solute
    ! ! SB_ 3D
    ! dsl4jsb_Real3D_onChunk :: noy
    ! dsl4jsb_Real3D_onChunk :: n2
    !------------------------------------------------------------------------------------------------------ !
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    !------------------------------------------------------------------------------------------------------ !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)          RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. options%iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    !------------------------------------------------------------------------------------------------------ !
    model         => Get_model(tile%owner_model_id)
    hgrid         => Get_grid(model%grid_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    !------------------------------------------------------------------------------------------------------ !
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(VEG_)
    !------------------------------------------------------------------------------------------------------ !
    dsl4jsb_Get_pool(SB_,  sb_main_pool%sb_pool,                sb_pool)
    dsl4jsb_Get_pool(VEG_, veg_main_pool%veg_pool,              veg_pool)
    CALL sb_pool%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_pool%Associate_chunk_var_ptr(ics, ice, iblk)
    !------------------------------------------------------------------------------------------------------ !
    ! ! SB_ 2D
    ! dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_solute) 
    ! dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_solute) 
    ! ! SB_ 3D
    ! dsl4jsb_Get_var3D_onChunk(SB_,       noy) 
    ! dsl4jsb_Get_var3D_onChunk(SB_,       n2) 

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------------------------ !
    !> some test calculations for runtime testing
    !> 
    ! sb_pool
    sb_pool = sb_pool * 2.0_wp
    sb_pool = sb_pool * 3.0_wp
    sb_pool = sb_pool + 20.0_wp
    sb_pool = sb_pool + 40.0_wp
    sb_pool = sb_pool - 40.0_wp
    sb_pool = sb_pool - 20.0_wp
    sb_pool = sb_pool / 3.0_wp
    sb_pool = sb_pool / 2.0_wp
    ! veg_pool
    veg_pool = veg_pool * 2.0_wp
    veg_pool = veg_pool * 3.0_wp
    veg_pool = veg_pool + 20.0_wp
    veg_pool = veg_pool + 40.0_wp
    veg_pool = veg_pool - 40.0_wp
    veg_pool = veg_pool - 20.0_wp
    veg_pool = veg_pool / 3.0_wp
    veg_pool = veg_pool / 2.0_wp

    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_sb_pool_test

  ! ======================================================================================================= !
  !> test runtime of bgc_material structures using the sb_bgcm_pool matrix
  !> runtime comparison with update_sb_pool_test
  !>
  SUBROUTINE update_sb_bgcm_pool_test(tile, options)
    USE mo_kind,                  ONLY: wp
    USE mo_jsb_control,           ONLY: debug_on
    USE mo_exception,             ONLY: message, message_text, finish
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid, t_jsb_grid
    USE mo_jsb_grid,              ONLY: Get_vgrid, Get_grid
    USE mo_lnd_bgcm_class,        ONLY: ELEM_C_ID, ELEM_N_ID, ELEM_P_ID, ELEM_C13_ID, ELEM_C14_ID, ELEM_N15_ID, &
      &                                 LAST_ELEM_ID
    USE mo_lnd_bgcm_store_class,  ONLY: SB_BGCM_POOL_ID, VEG_BGCM_POOL_ID
    USE mo_lnd_bgcm_store,        ONLY: t_lnd_bgcm_store
    USE mo_veg_config_class,      ONLY: get_number_of_veg_compartments, &
      &                                 VEG_PART_LEAF_IDX , VEG_PART_FINE_ROOT_IDX, VEG_PART_COARSE_ROOT_IDX, &
      &                                 VEG_PART_SAP_WOOD_IDX, VEG_PART_HEART_WOOD_IDX, VEG_PART_LABILE_IDX, &
      &                                 VEG_PART_RESERVE_IDX, VEG_PART_FRUIT_IDX, VEG_PART_SEED_BED_IDX
    USE mo_sb_config_class,       ONLY: get_number_of_sb_compartments,                                             &
      &                                 SB_PART_DOM_IDX, SB_PART_DOM_ASSOC_IDX, SB_PART_SOLUBLE_LITTER_IDX,        &
      &                                 SB_PART_POLYMERIC_LITTER_IDX, SB_PART_WOODY_LITTER_IDX, SB_PART_FUNGI_IDX, &
      &                                 SB_PART_MYCORRHIZA_IDX, SB_PART_MICROBIAL_IDX, SB_PART_RESIDUE_IDX,        &
      &                                 SB_PART_RESIDUE_ASSOC_IDX
    dsl4jsb_Use_processes SB_, VEG_
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_config(VEG_)
    ! dsl4jsb_Use_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_jsb_grid),       POINTER       :: hgrid                !< Horizontal grid
    TYPE(t_jsb_vgrid),      POINTER       :: vgrid_soil_sb        !< Vertical grid
    TYPE(t_lnd_bgcm_store), POINTER       :: bgcm_store           !< the bgcm store of this tile
    INTEGER                               :: nsoil_sb             !< number of soil layers as used/defined by the SB_ process
    INTEGER                               :: elem_idx_map(LAST_ELEM_ID - 1)  !< element mapper ID -> IND
    INTEGER                               :: isoil                !< loop over nsoil_sb
    INTEGER                               :: iblk, ics, ice, nc   !< dimensions block/chunk
    INTEGER                               :: nr_elements          !< number of elements
    INTEGER                               :: nr_sb_parts          !< number of sb compartments
    INTEGER                               :: nr_veg_parts         !< number of veg compartments
    INTEGER                               :: ixC                  !< index of the Carbon element in the bgcm_matrix
    INTEGER                               :: ixWooLit             !< short for SB_PART_WOODY_LITTER_IDX
    INTEGER                               :: ixCoaRoo             !< short for VEG_PART_COARSE_ROOT_IDX

    REAL(wp), ALLOCATABLE                 :: sb_bgcm_pool_mt(:,:,:,:)
    REAL(wp), ALLOCATABLE                 :: veg_bgcm_pool_mt(:,:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_sb_bgcm_pool_test'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_config(VEG_)
    ! dsl4jsb_Def_memory(SB_)
    ! dsl4jsb_Def_memory(VEG_)
    !------------------------------------------------------------------------------------------------------ !
    ! ! SB_ 2D
    ! dsl4jsb_Real2D_onChunk :: leaching_nh4_solute
    ! dsl4jsb_Real2D_onChunk :: leaching_no3_solute
    ! ! SB_ 3D
    ! dsl4jsb_Real3D_onChunk :: noy
    ! dsl4jsb_Real3D_onChunk :: n2
    !------------------------------------------------------------------------------------------------------ !
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    !------------------------------------------------------------------------------------------------------ !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)          RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. options%iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    !------------------------------------------------------------------------------------------------------ !
    model         => Get_model(tile%owner_model_id)
    hgrid         => Get_grid(model%grid_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    !------------------------------------------------------------------------------------------------------ !
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_config(VEG_)
    !------------------------------------------------------------------------------------------------------ !
    elem_idx_map(:) = model%config%elements_index_map(:)
    nr_elements     = model%config%nr_of_used_elements
    nr_sb_parts     = get_number_of_sb_compartments()
    nr_veg_parts    = get_number_of_veg_compartments()
    !------------------------------------------------------------------------------------------------------ !
    ixC       = elem_idx_map(ELEM_C_ID)
    ixWooLit  = SB_PART_WOODY_LITTER_IDX
    ixCoaRoo  = VEG_PART_COARSE_ROOT_IDX
    !------------------------------------------------------------------------------------------------------ !
    ALLOCATE(sb_bgcm_pool_mt(nr_sb_parts, nr_elements, nc, nsoil_sb))
    ALLOCATE(veg_bgcm_pool_mt(nr_veg_parts, nr_elements, nc))

    bgcm_store => tile%bgcm_store
    CALL bgcm_store%Get_matrix_from_bgcm_2l(SB_BGCM_POOL_ID, ics, ice, iblk, tile%name, sb_bgcm_pool_mt)
    CALL bgcm_store%Get_matrix_from_bgcm_2l(VEG_BGCM_POOL_ID, ics, ice, iblk, tile%name, veg_bgcm_pool_mt)

    !------------------------------------------------------------------------------------------------------ !
    ! ! SB_ 2D
    ! dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_solute) 
    ! dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_solute) 
    ! ! SB_ 3D
    ! dsl4jsb_Get_var3D_onChunk(SB_,       noy) 
    ! dsl4jsb_Get_var3D_onChunk(SB_,       n2) 

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------------------------ !
    !> some test calculations for runtime testing
    !> the same operations as with the "old pools" in update_sb_pool_test()
    !> 
    ! sb_pool
    ! some calculations, no changes
    sb_bgcm_pool_mt(:,:,:,:) = sb_bgcm_pool_mt(:,:,:,:) * 2.0_wp
    sb_bgcm_pool_mt(:,:,:,:) = sb_bgcm_pool_mt(:,:,:,:) * 3.0_wp
    sb_bgcm_pool_mt(:,:,:,:) = sb_bgcm_pool_mt(:,:,:,:) + 20.0_wp
    sb_bgcm_pool_mt(:,:,:,:) = sb_bgcm_pool_mt(:,:,:,:) + 40.0_wp
    sb_bgcm_pool_mt(:,:,:,:) = sb_bgcm_pool_mt(:,:,:,:) - 40.0_wp
    sb_bgcm_pool_mt(:,:,:,:) = sb_bgcm_pool_mt(:,:,:,:) - 20.0_wp
    sb_bgcm_pool_mt(:,:,:,:) = sb_bgcm_pool_mt(:,:,:,:) / 3.0_wp
    sb_bgcm_pool_mt(:,:,:,:) = sb_bgcm_pool_mt(:,:,:,:) / 2.0_wp
    ! calc only with: 
    !   SB_  woody litter, Carbon, sl 1
    !   VEG_ coarse roots, Carbon
    sb_bgcm_pool_mt(ixWooLit, ixC, :, 1) = sb_bgcm_pool_mt(ixWooLit, ixC, :, 1) + &
      &                                            (veg_bgcm_pool_mt(ixCoaRoo, ixC, :) * 2.5_wp)
    sb_bgcm_pool_mt(ixWooLit, ixC, :, 1) = sb_bgcm_pool_mt(ixWooLit, ixC, :, 1) - &
      &                                            (veg_bgcm_pool_mt(ixCoaRoo, ixC, :) * 2.5_wp)
    ! double the values (does not affect the next timestep: see 'update_copy_pool_to_bgcm_test()')
    sb_bgcm_pool_mt(:,:,:,:) = sb_bgcm_pool_mt(:,:,:,:) * 2.0_wp
    ! veg_pool
    veg_bgcm_pool_mt(:,:,:)  = veg_bgcm_pool_mt(:,:,:)  * 2.0_wp
    veg_bgcm_pool_mt(:,:,:)  = veg_bgcm_pool_mt(:,:,:)  * 3.0_wp
    veg_bgcm_pool_mt(:,:,:)  = veg_bgcm_pool_mt(:,:,:)  + 20.0_wp
    veg_bgcm_pool_mt(:,:,:)  = veg_bgcm_pool_mt(:,:,:)  + 40.0_wp
    veg_bgcm_pool_mt(:,:,:)  = veg_bgcm_pool_mt(:,:,:)  - 40.0_wp
    veg_bgcm_pool_mt(:,:,:)  = veg_bgcm_pool_mt(:,:,:)  - 20.0_wp
    veg_bgcm_pool_mt(:,:,:)  = veg_bgcm_pool_mt(:,:,:)  / 3.0_wp
    veg_bgcm_pool_mt(:,:,:)  = veg_bgcm_pool_mt(:,:,:)  / 2.0_wp

    ! write back data from bgcm_matrix to bgc_material, and de-allocate local bgcm_matrix
    CALL bgcm_store%Write_matrix_to_bgcm_2l(SB_BGCM_POOL_ID, ics, ice, iblk, tile%name, sb_bgcm_pool_mt)
    CALL bgcm_store%Write_matrix_to_bgcm_2l(VEG_BGCM_POOL_ID, ics, ice, iblk, tile%name, veg_bgcm_pool_mt)
    DEALLOCATE(sb_bgcm_pool_mt, veg_bgcm_pool_mt)

    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_sb_bgcm_pool_test

  ! ------------------------------------------------------------------------------------------------------- !
  ! aggregate_* routines
  ! ------------------------------------------------------------------------------------------------------- !

  ! ======================================================================================================= !
  !>Implementation of "aggregate": reset_sb_fluxes task
  !>
  SUBROUTINE aggregate_reset_sb_fluxes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_reset_sb_fluxes'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)

    ! nothing to aggregate, yet

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE aggregate_reset_sb_fluxes

  ! ======================================================================================================= !
  !>Implementation of "aggregate": sb_rate_modifier task
  !>
  SUBROUTINE aggregate_sb_rate_modifier(tile, options)
  
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_sb_rate_modifier'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)

    dsl4jsb_Aggregate_onChunk(SB_, rtm_decomposition              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_depolymerisation           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_mic_uptake                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_plant_uptake               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_sorption                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_desorption                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_hsc                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_nitrification              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_denitrification            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_gasdiffusion               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_decomposition              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_depolymerisation           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_mic_uptake                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_plant_uptake               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_sorption                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_desorption                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_hsc                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_nitrification              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_gasdiffusion               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, anaerobic_volume_fraction_sl   , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  
  END SUBROUTINE aggregate_sb_rate_modifier
  
  ! ======================================================================================================= !
  !>Implementation of "aggregate": sb_simple_model task
  !>
  SUBROUTINE aggregate_sb_simple_model(tile, options)
  
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_sb_simple_model'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(VEG_)

    ! VEG 2D
    dsl4jsb_Aggregate_onChunk(VEG_, unit_uptake_n_pot               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_uptake_p_pot               , weighted_by_fract)
! TODO  add missing variables for aggregation
    ! SB_ 2D
    ! dsl4jsb_Aggregate_onChunk(SB_,    , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_,    , weighted_by_fract)
    ! SB_ 3D
    dsl4jsb_Aggregate_onChunk(SB_,  plant_uptake_nh4_sl             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  plant_uptake_no3_sl             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  plant_uptake_po4_sl             , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_,    , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_,    , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  
  END SUBROUTINE aggregate_sb_simple_model

  ! !-----------------------------------------------------------------------------------------------------
  ! !> Implementation of "aggregate": update_sb_jsm task
  ! !!  actually NOT using this routine, i.e., no var is aggregated here
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE aggr_update_sb_jsm(tile, options)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
  !   TYPE(t_jsb_task_options),   INTENT(in)    :: options
  !   ! ---------------------------
  !   ! 0.2 Local
  !   TYPE(t_jsb_model),        POINTER         :: model
  !   CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
  !   INTEGER                                   :: iblk, ics, ice, nc
  !   !X if necessary: REAL(wp) :: dtime, steplen
  !   CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_sb_jsm'
  !   ! ---------------------------
  !   ! 0.3 Declare Memory
  !   ! Declare process configuration and memory Pointers
  !   dsl4jsb_Def_config(SB_)
  !   dsl4jsb_Def_memory(SB_)
  !   ! Get local variables from options argument
  !   iblk    = options%iblk
  !   ics     = options%ics
  !   ice     = options%ice
  !   nc      = options%nc
  !   !X if necessary: dtime   = options%dtime
  !   !X if necessary: steplen = options%steplen
  !   ! ---------------------------
  !   ! 0.4 Process Activity, Debug Option
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
  !   ! ---------------------------
  !   ! 0.5 Get Memory
  !   ! Get process variables (Set pointers to variables in memory)
  !   model => Get_model(tile%owner_model_id)
  !   ! Get process config
  !   dsl4jsb_Get_config(SB_)
  !   ! Get process memories
  !   dsl4jsb_Get_memory(SB_)
  ! 
  ! 
  !   ! ------------------------------------------------------------------------------------------------------------
  !   ! Go aggregate
  !   ! ------------------------------------------------------------------------------------------------------------
  ! 
  ! 
  !   IF (tile%Has_children()) THEN
  !     weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
  !     ! do not aggregate after this routine, only after update_sb_pools()
  !   END IF
  ! 
  !   !X Implementation: Start your process scheme here...
  !   !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.
  ! 
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  ! 
  ! END SUBROUTINE aggr_update_sb_jsm

  ! ======================================================================================================= !
  !>Implementation of "aggregate": sb_pools task
  !> 
  SUBROUTINE aggregate_sb_pools(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_sb_pools'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)

! TODO  add missing variables for aggregation
    ! 2D
    dsl4jsb_Aggregate_onChunk(SB_, total_flux_carbon            , weighted_by_fract)
    ! 3D
    dsl4jsb_Aggregate_onChunk(SB_, nh4_solute                   , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_sb_pools

  ! ======================================================================================================= !
  !>Implementation of "aggregate": time_average_assimilation task
  !>
  SUBROUTINE aggregate_time_average_soilbiogeochemistry(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_time_average_soilbiogeochemistry'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)

    dsl4jsb_Aggregate_onChunk(SB_, microbial_cue_eff_tmic_mavg    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, microbial_nue_eff_tmic_mavg    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, microbial_pue_eff_tmic_mavg    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_c_tmyc_mavg_sl      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_n_tmyc_mavg_sl      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_n_tlabile_mavg_sl   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_p_tlabile_mavg_sl   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, dom_cn_mavg                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, dom_cp_mavg                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, enzyme_frac_poly_c_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, enzyme_frac_poly_n_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, enzyme_frac_poly_p_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_c            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_n            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_p            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_c13          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_c14          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_n15          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_c         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_n         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_p         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_c13       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_c14       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_n15       , weighted_by_fract)


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_time_average_soilbiogeochemistry
  
  ! ======================================================================================================= !
  !> Implementation of "aggregate": copy_pool_to_bgcm_test task
  !>
  !> for testing only
  !>
  SUBROUTINE aggregate_copy_pool_to_bgcm_test(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_copy_pool_to_bgcm_test'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model             => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)

    ! so far nothing TODO

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_copy_pool_to_bgcm_test
  
  ! ======================================================================================================= !
  !> Implementation of "aggregate": sb_pool_test task
  !>
  !> for testing only
  !>
  SUBROUTINE aggregate_sb_pool_test(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_sb_pool_test'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model             => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)

    ! so far nothing TODO

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_sb_pool_test
  
  ! ======================================================================================================= !
  !> Implementation of "aggregate": sb_pool_test task
  !>
  !> for testing only
  !>
  SUBROUTINE aggregate_sb_bgcm_pool_test(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_sb_bgcm_pool_test'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model             => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)

    ! so far nothing TODO

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_sb_bgcm_pool_test

#endif
END MODULE mo_sb_interface
