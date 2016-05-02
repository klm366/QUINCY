!>
!! @brief Contains the interfaces to the assimi process
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-05-30)
!!
MODULE mo_assimi_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module
  
  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message
  
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model  
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module (Get_assimi_memory and Get_assimi_config)
  dsl4jsb_Use_processes ASSIMI
  dsl4jsb_Use_processes VEG

  ! Use of process configurations (t_assimi_config)
  dsl4jsb_Use_config(ASSIMI)

  ! Use of process memories (t_assimi_memory)
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
  dsl4jsb_Use_memory(ASSIMI)
  dsl4jsb_Use_memory(VEG)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables 
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_assimi_tasks !,t_assimi_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_interface'
  
  !> Type definition for assimi process
!   TYPE, EXTENDS(t_jsb_process) :: t_assimi_process
!   CONTAINS
!     PROCEDURE :: Register_tasks => Register_assimi_tasks        !< Registers tasks for assimi process
!   END TYPE t_assimi_process
  
  !> Type definition for assimi task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_assimi
    !LOGICAL      :: waterLimitationFlag
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => assimilation_waterlimited   !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_assimi            !< Aggregates computed task variables
  END TYPE tsk_assimi
  
  !> Constructor interface for assimi task
  INTERFACE tsk_assimi
    PROCEDURE Create_task_assimi                        !< Constructor function for task
  END INTERFACE tsk_assimi

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks
  
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for assimi task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "assimi"
  !!
  FUNCTION Create_task_assimi(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_assimi::return_ptr)

    return_ptr%name = 'assimi'
    return_ptr%process_id = ASSIMI
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_assimi
  
  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for assimi process
  !!
  !! @param[in,out] this      Instance of assimi process class  
  !! @param[in]     model_id  Model id  
  !!
  SUBROUTINE Register_assimi_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_assimi(model_id))

  END SUBROUTINE Register_assimi_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation to calculate the height scaling factors of plant enzyms for assimilation
  !! Task "assimi" calculates the carbon assimilation of plants.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  ! The scaling factors represent the enzymatic distribution of plants in their canopy layers. This correltates
  ! with the amount of nitrogen in these layers (therefore called nitrogen_scaling_factors in JSBACH3).
  ! Eqs. (106), (107) in Knorr (1997)
  ! Note that this routine should only be applied to special vegetation classes    
  SUBROUTINE assimilation_scaling_factors(tile, options)
  
    ! Declarations
    !USE mo_orbit,                 ONLY: inquire_declination
     ! R: Returns declination calculated in last orbit call.
     !    Problem: FATAL ERROR in inquire_declination: Not Initialized
     !    This is because ICON did not call the calculation of orbital parameters in SUBROUTINE orbit_vsop87 or 
     !    alternatively in SUBROUTINE orbit_kepler. Therefore this seems to be a problem solved with the development of ICON...
     !    Note: problem still exists in the 26.10.2015 version of icon-les-dev !
     !    However taken from mo_psrad_orbit works... perhaps precompiler command necessary here...
    USE mo_psrad_orbit,           ONLY: inquire_declination
     
    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER        :: model
    INTEGER                           :: iblk, ics, ice, nc
    TYPE(t_jsb_grid),  POINTER        :: grid
    REAL(wp)                          :: declination ! Solar Declination angle
    INTEGER                           :: i
    
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':assimilation_scaling_factors'
    
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI)
    dsl4jsb_Def_memory(ASSIMI)
    dsl4jsb_Def_memory(VEG)
    
    ! Declare pointers to variables in memory
    INTEGER, POINTER       :: ncanopy          ! number of canopy layers
    dsl4jsb_Real3D_onChunk :: canopy_bound_lai
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: lat    
    
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    
    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(ASSIMI)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    
    model => Get_model(tile%owner_model_id)
    
    ! Get process config
    dsl4jsb_Get_config(ASSIMI) !assimi_conf => Get_assimi_config(model) -->  model%processes(iproc)%p%config
    
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI)
    dsl4jsb_Get_memory(VEG)
    
    ! Get process variables
    ncanopy                  => dsl4jsb_Config(ASSIMI)%config%ncanopy          ! IN
    canopy_bound_lai         => dsl4jsb_Config(ASSIMI)%%canopy_bound_lai ! IN
    grid => Get_grid(model%grid_id)                           ! IN
    lat  => grid%lat(ics:ice, iblk)                           ! IN
    
    dsl4jsb_Get_var2D_onChunk(VEG, lai)                       ! IN  
    dsl4jsb_Get_var3D_onChunk(ASSIMI, assimi_scaling_fact_cl) ! OUT
    
    ! ---------------------------
    ! Go         
    
    CALL inquire_declination(declination) ! R: do not put this in elemental procedure as it would not be pure...
                                          ! e.g. declination = 0.006918_wp 

    DO i=1,ncanopy
       CALL calc_assimilation_scaling_factors(declination,                 & ! input
                                              canopy_bound_lai(i-1),       & ! input
                                              lai(:),                      & ! input
                                              lat(:),                      & ! input  
                                              assimi_scaling_fact_cl(:,i)  & ! output
                                              )
    END DO
  
  RETURN  
      
!     R: to do when the jsbach4 "case" structure is included...
!     DO i=1,ncanopy
!        WHERE(.NOT. mask(:) .OR. .NOT. lctlib%NitrogenScalingFlag(lct(:)))
!           assimi_scaling_fact_cl(:,i) = 1._wp ! to some vegetation classes nitrogen scaling is not applied
!        END WHERE
!     END DO    
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
    
  END SUBROUTINE update_assimi
  
  ! ================================================================================================================================
  !>
  !! @brief Implementation of assimilation for water unlimit conditions
  !! Task "assimi" calculates the carbon assimilation of plants.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE assimilation_waterunlimited(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    TYPE(t_jsb_model), POINTER                :: model
 
    ! Local variables
    INTEGER  :: iblk, ics, ice, nc

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_assimi'
   
    ! Declare process memories
    dsl4jsb_Def_memory(ASSIMI)

    ! Declare pointers to variables in memory
    !X e.g: dsl4jsb_Real2D_onChunk :: t_srf

    
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    
    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(ASSIMI)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    
    model => Get_model(tile%owner_model_id)
    
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI)

    ! Get process variables
    !X e.g: dsl4jsb_Get_var2D_onChunk(ASSIMI, t_srf)     ! IN

    !X Implementation: Start your process scheme here...
    ! ---------------------------
    ! Go
    
    
    
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
    
  END SUBROUTINE update_assimi
    
  
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "assimi"
  !! 
  !! @param[in,out] tile    Tile for which routine is executed.  
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.  
  !!
  SUBROUTINE aggregate_assimi(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    TYPE(t_jsb_model), POINTER                :: model
    
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_assimi'

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen
    
    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    
    !X Implementation: Start your process scheme here...
    
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_assimi

  
  ! ================================================================================================================================
  !>
  !! @brief Another task for "assimi"
  !!
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.  
  !! @param[in,out] tile    Tile for which routine is executed.
  !  ....
  
#endif
END MODULE mo_assimi_interface
