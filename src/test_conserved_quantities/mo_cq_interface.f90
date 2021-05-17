!> Summary: interface to the cq process (test conserved quantities)
!>
!>### Contains the interfaces to the cq process which is a small test process containing a couple of conserved quantities 
!>
!> Note: 1. This module takes advantage of the handling of the integrate and aggregate routines, which both
!>          call the same routines in this module (integrate = aggregate). 
!>          Thereby each tile is processed by the same routine for the tasks of this process.
!>       2. currently only 2D variables are dealt with -- if also 3D variables need to be conserved,
!>          these will anyway have to specially be dealt with in the concept and thus need to be 
!>          identified with an own cqt e.g. distinguish WATER_2D_CQ_TYPE and WATER_3D_CQ_TYPE?
!>          Also, integration of pool-structure might need some additional work.
!>
!> Copyright (c) 2021 Julia Nabel
!>
!>#### Revision History
!> First version and ongoing changes   by Julia Nabel (2021-03-09)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
MODULE mo_cq_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model    
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes cq_, pheno_

  ! Use of process configurations
  dsl4jsb_Use_config(cq_)

  ! Use of process memories
  dsl4jsb_Use_memory(cq_)
  dsl4jsb_Use_memory(pheno_)

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Register_cq_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_cq_interface'
  
  !> Type definition for cq task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_cq
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_cq     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_cq  !< Aggregates computed task variables
  END TYPE tsk_cq
  
  !> Constructor interface for cq task
  INTERFACE tsk_cq
    PROCEDURE Create_task_cq                        !< Constructor function for task
  END INTERFACE tsk_cq
  
CONTAINS

  ! ====================================================================================================== !
  !
  !> Constructor for cq task
  !
  FUNCTION Create_task_cq(model_id) RESULT(return_ptr)

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id    !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "cq"
    ! -------------------------------------------------------------------------------------------------- !

    ALLOCATE(tsk_cq::return_ptr)
    CALL return_ptr%Construct(name='cq', process_id=cq_, owner_model_id=model_id)

  END FUNCTION Create_task_cq
    
  ! ====================================================================================================== !
  !
  !> Register tasks for cq process
  !
  SUBROUTINE Register_cq_tasks(this, model_id)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_process), INTENT(inout) :: this        !< Instance of cq process class
    INTEGER,               INTENT(in)   :: model_id    !< Model id
    ! -------------------------------------------------------------------------------------------------- !

    CALL this%Register_task(tsk_cq(model_id))
    
  END SUBROUTINE Register_cq_tasks

  ! ====================================================================================================== !
  !
  !> Implementation of "update" for task "cq"
  !>
  !> Calculates a ta variable for each cq variable
  !
  SUBROUTINE update_cq(tile, options)

    USE mo_jsb_parallel,      ONLY: Get_omp_thread
    USE mo_jsb_time,          ONLY: is_newday

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), POINTER :: current_tile
    TYPE(t_jsb_model), POINTER          :: model

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: dtime
    REAL(wp), DIMENSION(options%nc) :: moved_fract, current_fract
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_cq'

    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_memory(cq_)
    dsl4jsb_Def_memory(PHENO_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk ::  a_dead_c_cq
    dsl4jsb_Real2D_onChunk ::  another_dead_c_cq
    dsl4jsb_Real2D_onChunk ::  a_veg_c_cq
    dsl4jsb_Real2D_onChunk ::  an_implicit_scaling_cq

    dsl4jsb_Real2D_onChunk ::  a_dead_c_ta_cq
    dsl4jsb_Real2D_onChunk ::  another_dead_c_ta_cq
    dsl4jsb_Real2D_onChunk ::  a_veg_c_ta_cq
    dsl4jsb_Real2D_onChunk ::  an_implicit_scaling_ta_cq

    dsl4jsb_Real2D_onChunk ::  veg_fract_correction
    dsl4jsb_Real2D_onChunk ::  fract_fpc_max    
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    
    model => Get_model(tile%owner_model_id)

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(cq_)) RETURN

    dsl4jsb_Get_memory(CQ_)
    dsl4jsb_Get_memory(PHENO_)

    dsl4jsb_Get_var2D_onChunk(CQ_,  a_dead_c_cq )
    dsl4jsb_Get_var2D_onChunk(CQ_,  another_dead_c_cq )
    dsl4jsb_Get_var2D_onChunk(CQ_,  a_veg_c_cq )
    dsl4jsb_Get_var2D_onChunk(CQ_,  an_implicit_scaling_cq )

    dsl4jsb_Get_var2D_onChunk(CQ_,  a_dead_c_ta_cq )
    dsl4jsb_Get_var2D_onChunk(CQ_,  another_dead_c_ta_cq )
    dsl4jsb_Get_var2D_onChunk(CQ_,  a_veg_c_ta_cq )
    dsl4jsb_Get_var2D_onChunk(CQ_,  an_implicit_scaling_ta_cq )

    dsl4jsb_Get_var2D_onChunk(PHENO_,  fract_fpc_max )
    dsl4jsb_Get_var2D_onChunk(PHENO_,  veg_fract_correction )

    a_veg_c_ta_cq = a_veg_c_cq * veg_fract_correction * fract_fpc_max
    a_dead_c_ta_cq = a_dead_c_cq * veg_fract_correction * fract_fpc_max
    another_dead_c_ta_cq = another_dead_c_cq * veg_fract_correction * fract_fpc_max
    an_implicit_scaling_ta_cq = an_implicit_scaling_cq * veg_fract_correction * fract_fpc_max
    
  END SUBROUTINE update_cq
  
  ! ====================================================================================================== !
  !
  !> Implementation of "aggregate" for task "cq"
  !>
  !> Aggregates ta variables associated with the cq variables
  !
  SUBROUTINE aggregate_cq(tile, options)

    USE mo_jsb_time,          ONLY: is_newday

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(CQ_)
    dsl4jsb_Def_memory(CQ_)

    TYPE(t_jsb_model),       POINTER :: model
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice
    REAL(wp) :: new_sum(options%nc)
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_cq'
    
    dsl4jsb_Real2D_onChunk ::  old_sum
    dsl4jsb_Real2D_onChunk ::  conservation_test_field
    dsl4jsb_Real2D_onChunk ::  a_dead_c_ta_cq
    dsl4jsb_Real2D_onChunk ::  another_dead_c_ta_cq
    dsl4jsb_Real2D_onChunk ::  a_veg_c_ta_cq
    dsl4jsb_Real2D_onChunk ::  an_implicit_scaling_ta_cq
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    
    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(CQ_)
    dsl4jsb_Get_memory(CQ_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(CQ_, a_dead_c_ta_cq, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CQ_, a_veg_c_ta_cq, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CQ_, another_dead_c_ta_cq, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CQ_, an_implicit_scaling_ta_cq, weighted_by_fract)

    IF(is_newday(options%current_datetime,options%dtime)) THEN
      IF ( tile%name .EQ. 'veg') THEN
        dsl4jsb_Get_var2D_onChunk(CQ_,  old_sum)
        dsl4jsb_Get_var2D_onChunk(CQ_,  conservation_test_field)
        dsl4jsb_Get_var2D_onChunk(CQ_,  a_veg_c_ta_cq)
        dsl4jsb_Get_var2D_onChunk(CQ_,  a_dead_c_ta_cq)
        dsl4jsb_Get_var2D_onChunk(CQ_,  another_dead_c_ta_cq)
        dsl4jsb_Get_var2D_onChunk(CQ_,  an_implicit_scaling_ta_cq)
        new_sum = a_veg_c_ta_cq + a_dead_c_ta_cq + another_dead_c_ta_cq + an_implicit_scaling_ta_cq
        
        conservation_test_field = new_sum - old_sum
        old_sum = new_sum
      ENDIF
    ENDIF


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_cq

#endif
END MODULE mo_cq_interface
