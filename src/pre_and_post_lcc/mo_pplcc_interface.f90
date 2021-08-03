!> Summary: interface to the pplcc process (pre- and postprocessing for all land cover change processes)
!>
!>### Contains the interfaces to the pplcc process which does the pre- and postprocessing for all lcc processes
!>
!> Note: 
!>
!>       1. This module takes advantage of the handling of the integrate and aggregate routines, which both
!>          call the same routines in this module (integrate = aggregate). 
!>          Thereby each tile is processed by the same routine for the tasks of this process.
!>
!>       2. currently only 2D variables are dealt with -- if also 3D variables need to be conserved,
!>          these will anyway have to specially be dealt with in the concept and thus need to be 
!>          identified with an own cqt e.g. distinguish WATER_2D_CQ_TYPE and WATER_3D_CQ_TYPE?
!>          Also, integration of pool-structure might need some additional work.
!>
!> Copyright (c) 2021 Julia Nabel
!>
!>#### Revision History
!> First version and ongoing changes   by Julia Nabel (2021-02-02 - present)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
MODULE mo_pplcc_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model    
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  !USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes PPLCC_, PHENO_

  ! Use of process configurations
  !  dsl4jsb_Use_config(PPLCC_)

  ! Use of process memories
!  dsl4jsb_Use_memory(PPLCC_)
  dsl4jsb_Use_memory(PHENO_)  

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_pplcc_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pplcc_interface'
  
  !> Type definition for pplcc pre task (preprocessing for lcc)
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_pplcc_pre_lcc
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_pplcc_pre_lcc  !< prepares lcc processes for leaves
    PROCEDURE, NOPASS :: Aggregate => update_pplcc_pre_lcc  !< prepares lcc processes for all other nodes
  END TYPE tsk_pplcc_pre_lcc
  
  !> Constructor interface for pplcc pre task (preprocessing for lcc)
  INTERFACE tsk_pplcc_pre_lcc
    PROCEDURE Create_task_pplcc_pre_lcc                      !< Constructor function for task
  END INTERFACE tsk_pplcc_pre_lcc

  !> Type definition for pplcc post task (postprocessing for lcc)
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_pplcc_post_lcc
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_pplcc_post_lcc  !< does postprocessing for lcc processes on leaves
    PROCEDURE, NOPASS :: Aggregate => update_pplcc_post_lcc  !< does postprocessing for lcc processes on other nodes
  END TYPE tsk_pplcc_post_lcc
  
  !> Constructor interface for pplcc post task (postprocessing for lcc)
  INTERFACE tsk_pplcc_post_lcc
    PROCEDURE Create_task_pplcc_post_lcc                        !< Constructor function for task
  END INTERFACE tsk_pplcc_post_lcc
  
CONTAINS

  ! ====================================================================================================== !
  !
  !> Constructor for pplcc pre lcc task (preprocessing for lcc)
  !
  FUNCTION Create_task_pplcc_pre_lcc(model_id) RESULT(return_ptr)

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id    !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "pplcc_pre_lcc"
    ! -------------------------------------------------------------------------------------------------- !

    ALLOCATE(tsk_pplcc_pre_lcc::return_ptr)
    CALL return_ptr%Construct(name='pplcc_pre_lcc', process_id=PPLCC_, owner_model_id=model_id)

  END FUNCTION Create_task_pplcc_pre_lcc

  ! ====================================================================================================== !
  !
  !> Constructor for pplcc post lcc task (postprocessing for lcc)
  !
  FUNCTION Create_task_pplcc_post_lcc(model_id) RESULT(return_ptr)

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id    !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "pplcc_post_lcc"
    ! -------------------------------------------------------------------------------------------------- !

    ALLOCATE(tsk_pplcc_post_lcc::return_ptr)
    CALL return_ptr%Construct(name='pplcc_post_lcc', process_id=PPLCC_, owner_model_id=model_id)

  END FUNCTION Create_task_pplcc_post_lcc

  ! ====================================================================================================== !
  !
  !> Register tasks for pplcc process
  !
  SUBROUTINE Register_pplcc_tasks(this, model_id)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_process), INTENT(inout) :: this      !< Instance of pplcc process class
    INTEGER,              INTENT(in)    :: model_id  !< Model id
    ! -------------------------------------------------------------------------------------------------- !

    CALL this%Register_task(tsk_pplcc_pre_lcc(model_id))
    CALL this%Register_task(tsk_pplcc_post_lcc(model_id))
    
  END SUBROUTINE Register_pplcc_tasks

  ! ====================================================================================================== !
  !
  !> Implementation of "update" for task "pplcc_pre_lcc"
  !>
  !> Task "pplcc_pre_lcc" prepares lcc processes 
  !> -- for now: convert all cq vars to per grid cell.
  !
  SUBROUTINE update_pplcc_pre_lcc(tile, options)

    USE mo_jsb_time,          ONLY: is_newday
    USE mo_jsb_cqt_class,     ONLY: Get_cqt_name, LIVE_CARBON_CQ_TYPE, DEAD_CARBON_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER          :: model

    REAL(wp), DIMENSION(options%nc) :: conversion_factor

    INTEGER  :: iblk, ics, ice, cqt_id, i_cqt, i_cq_var
    REAL(wp) :: dtime
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_pplcc_pre_lcc'

    CHARACTER(len=:), ALLOCATABLE :: thisUnit, expectedUnit
    
    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_memory(PHENO_)
!    dsl4jsb_Def_memory(PPLCC_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: veg_fract_correction
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    dtime   = options%dtime
    
    model => Get_model(tile%owner_model_id)

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(PPLCC_)) RETURN

    ! If not newday, do nothing
    IF( .NOT. is_newday(options%current_datetime,dtime)) RETURN
    
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    IF( tile%Has_conserved_quantities()) THEN
      
      ! Iterate over all cq types and their vars and convert
      DO i_cqt = 1,tile%nr_of_cqts

        ! conversion of m2 tile to m2 box is required for several types of conserved quantities
        CALL tile%Get_fraction(ics, ice, iblk, fract=conversion_factor(:))

        ! conversion depends on the cq type
        cqt_id = tile%conserved_quantities(i_cqt)%p%type_id
        SELECT CASE (cqt_id)
          CASE (LIVE_CARBON_CQ_TYPE, DEAD_CARBON_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE)

            ! in the current implementation of jsb4 carbon variables are only available on the pft tiles
            IF (.NOT. index(tile%name, 'pft') /= 0) THEN
              CALL finish(TRIM(routine), 'Violation of precondition: carbon cq types only expected for pft tiles, ' &
                & //' but found on '// trim(tile%name))
            ENDIF

            ! in jsb4 the carbon variables are specified as mol per canopy
            expectedUnit = 'mol(C) m-2(canopy)'

            ! therefore the conversion factor needs to be adapted
            ! -> Assert that PHENO is active for this tile
            IF (.NOT. tile%Is_process_active(PHENO_)) THEN
              CALL finish(TRIM(routine), 'Violation of precondition: carbon cq types require pheno vars for conversion, ' &
                & //' but pheno is not active for '// trim(tile%name))
            ENDIF

            !- Set pointers to memory on the current tile
            dsl4jsb_Get_memory(PHENO_)       
            ! and get required variables on this tile 
            dsl4jsb_Get_var2D_onChunk(PHENO_, veg_fract_correction)    ! in
            dsl4jsb_Get_var2D_onChunk(PHENO_, fract_fpc_max)           ! in

            conversion_factor(:) = conversion_factor(:) * veg_fract_correction(:) * fract_fpc_max(:)

          ! CASE (WATER_CQ_TYPE)
          CASE DEFAULT
            !>
            !> Assertion: for each CQT a conversion needs to be specified in pplcc_pre_lcc 
            !> 
            CALL finish(TRIM(routine), 'Conversion factors unspecified for '//Get_cqt_name(cqt_id))
        END SELECT
          
        ! Convert all variables of this type
        DO i_cq_var=1,tile%conserved_quantities(i_cqt)%p%no_of_vars
          ! Assertion: variable needs to have the unit expected for this type
          thisUnit = tile%conserved_quantities(i_cqt)%p%cq_vars_2D(i_cq_var)%p%unit
          IF (.NOT. TRIM(thisUnit) .EQ. expectedUnit) THEN
            CALL finish(TRIM(routine), 'Violation of assertion: variable does not have the expected unit ' &
              & // expectedUnit // ' but '// trim(thisUnit))
          ENDIF

          IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), '..... and var '// &
            & TRIM(tile%conserved_quantities(i_cqt)%p%cq_vars_2D(i_cq_var)%p%name))

          !JN-TODO: make a function in a new process module? -> DSL?!
          tile%conserved_quantities(i_cqt)%p%cq_vars_2D(i_cq_var)%p%ptr2d(ics:ice,iblk) &
            & = tile%conserved_quantities(i_cqt)%p%cq_vars_2D(i_cq_var)%p%ptr2d(ics:ice,iblk) &
            &   * conversion_factor(:)

          ! for bookkeeping I would have wanted to change the unit, but seems not to be intended to be changed
        ENDDO
      ENDDO
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_pplcc_pre_lcc

  ! ====================================================================================================== !
  !
  !> Implementation of "update" for task "pplcc_post_lcc"
  !>
  !> Task "pplcc_post_lcc" does a postprocessing after all lcc processes were executed 
  !> -- for now: convert all cq vars back to their initial reference area.
  !
  SUBROUTINE update_pplcc_post_lcc(tile, options)

    USE mo_jsb_time,          ONLY: is_newday
    USE mo_jsb_cqt_class,     ONLY: Get_cqt_name, LIVE_CARBON_CQ_TYPE, DEAD_CARBON_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER          :: model

    REAL(wp), DIMENSION(options%nc) :: conversion_factor

    INTEGER  :: iblk, ics, ice, i_cqt, i_cq_var, cqt_id
    REAL(wp) :: dtime
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_pplcc_post_lcc'

    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_memory(PHENO_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: veg_fract_correction
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    dtime   = options%dtime
    
    model => Get_model(tile%owner_model_id)

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(PPLCC_)) RETURN

    ! If not newday, do nothing
    IF( .NOT. is_newday(options%current_datetime,dtime)) RETURN
    
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    IF( tile%Has_conserved_quantities()) THEN

      ! Iterate over all cq types and their vars and convert and change unit
      DO i_cqt = 1,tile%nr_of_cqts
      ! current tile fraction is required for several types of conserved quantities
        CALL tile%Get_fraction(ics, ice, iblk, fract=conversion_factor(:))

        ! conversion depends on the cq type
        cqt_id = tile%conserved_quantities(i_cqt)%p%type_id

        SELECT CASE (cqt_id)
          CASE (LIVE_CARBON_CQ_TYPE, DEAD_CARBON_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE)

            ! in the current implementation of jsb4 carbon variables are only available on the pft tiles
            IF (.NOT. index(tile%name, 'pft') /= 0) THEN
              CALL finish(TRIM(routine), 'Violation of precondition: carbon cq types only expected for pft tiles, ' &
                & //' but found on '// trim(tile%name))
            ENDIF

            ! in jsb4 the carbon variables are specified as mol per canopy, therefore adapt the conversion factor
            ! -> Assert that PHENO is active for this tile
            IF (.NOT. tile%Is_process_active(PHENO_)) THEN
              CALL finish(TRIM(routine), 'Violation of precondition: carbon cq types require pheno vars for conversion, ' &
                & //' but pheno is not active for '// trim(tile%name))
            ENDIF

            !- Set pointers to memory on the current tile
            dsl4jsb_Get_memory(PHENO_)       
            ! and get required variables on this tile 
            dsl4jsb_Get_var2D_onChunk(PHENO_, veg_fract_correction)    ! in
            dsl4jsb_Get_var2D_onChunk(PHENO_, fract_fpc_max)           ! in
            conversion_factor(:) = conversion_factor(:) * veg_fract_correction(:) * fract_fpc_max(:)

          ! CASE (WATER_CQ_TYPE)
          CASE DEFAULT
            !>
            !> Assertion: for each CQT a conversion needs to be specified in pplcc_post_lcc 
            !> 
            CALL finish(TRIM(routine), 'Conversion factors unspecified for '//Get_cqt_name(cqt_id))
        END SELECT
        
        ! Convert all variables of this type
        DO i_cq_var=1,tile%conserved_quantities(i_cqt)%p%no_of_vars

          IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), '..... and var '// &
            & TRIM(tile%conserved_quantities(i_cqt)%p%cq_vars_2D(i_cq_var)%p%name))

          !JN-TODO: make a function in new process module?
          ! only do the conversion where the conversion factor is not zero
          WHERE(conversion_factor(:) > 0.0)
            tile%conserved_quantities(i_cqt)%p%cq_vars_2D(i_cq_var)%p%ptr2d(ics:ice,iblk) &
              & = tile%conserved_quantities(i_cqt)%p%cq_vars_2D(i_cq_var)%p%ptr2d(ics:ice,iblk) &
              &   / conversion_factor(:)
          ELSEWHERE
            tile%conserved_quantities(i_cqt)%p%cq_vars_2D(i_cq_var)%p%ptr2d(ics:ice,iblk) = 0.0
          ENDWHERE

        ENDDO
      ENDDO
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_pplcc_post_lcc

#endif
END MODULE mo_pplcc_interface
