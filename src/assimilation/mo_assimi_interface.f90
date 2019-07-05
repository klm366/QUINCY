!>
!! @brief Contains the interfaces to the assimi process
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-05-30)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_assimi_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  !USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module (Get_assimi_memory and Get_assimi_config)
  dsl4jsb_Use_processes ASSIMI_, HYDRO_, PHENO_, RAD_, A2L_, SEB_

  ! Use of process configurations (t_assimi_config)
  dsl4jsb_Use_config(ASSIMI_)

  ! Use of process memories (t_assimi_memory)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(SEB_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_assimi_tasks !,t_assimi_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_interface'

  !> Type definition for assimilation_scaling_factors task
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_assimilation_scaling_factors
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_assimilation_scaling_factors    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_assimilation_scaling_factors !< Aggregates computed task variables
  END TYPE   tsk_assimilation_scaling_factors

  !> Constructor interface for assimilation_scaling_factors task
  INTERFACE   tsk_assimilation_scaling_factors
    PROCEDURE Create_task_assimilation_scaling_factors         !< Constructor function for task
  END INTERFACE   tsk_assimilation_scaling_factors

  !> Type definition for canopy_cond_unstressed task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_cond_unstressed
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_cond_unstressed    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_cond_unstressed !< Aggregates computed task variables
  END TYPE tsk_canopy_cond_unstressed

  !> Constructor interface for canopy_cond_unstressed task
  INTERFACE tsk_canopy_cond_unstressed
    PROCEDURE Create_task_canopy_cond_unstressed     !< Constructor function for task
  END INTERFACE tsk_canopy_cond_unstressed

  !> Type definition for canopy_cond_stressed task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_cond_stressed
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_cond_stressed    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_cond_stressed !< Aggregates computed task variables
  END TYPE tsk_canopy_cond_stressed

  !> Constructor interface for canopy_cond_stressed task
  INTERFACE tsk_canopy_cond_stressed
    PROCEDURE Create_task_canopy_cond_stressed     !< Constructor function for task
  END INTERFACE tsk_canopy_cond_stressed

  !> Type definition for assimilation task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_assimilation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_assimilation    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_assimilation !< Aggregates computed task variables
  END TYPE tsk_assimilation

  !> Constructor interface for assimilation task
  INTERFACE tsk_assimilation
    PROCEDURE Create_task_assimilation     !< Constructor function for task
  END INTERFACE tsk_assimilation

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for assimi task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "assimi"
  !!
  FUNCTION Create_task_assimilation_scaling_factors(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_assimilation_scaling_factors::return_ptr)
    CALL return_ptr%Construct(name='assimilation_scaling_factors', process_id=ASSIMI_, owner_model_id=model_id)

  END FUNCTION Create_task_assimilation_scaling_factors

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for canopy_cond_unstressed task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "canopy_cond_unstressed"
  !!
  FUNCTION Create_task_canopy_cond_unstressed(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_canopy_cond_unstressed::return_ptr)
    CALL return_ptr%Construct(name='canopy_cond_unstressed', process_id=ASSIMI_, owner_model_id=model_id)

  END FUNCTION Create_task_canopy_cond_unstressed

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for canopy_cond_stressed task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "canopy_cond_stressed"
  !!
  FUNCTION Create_task_canopy_cond_stressed(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_canopy_cond_stressed::return_ptr)
    CALL return_ptr%Construct(name='canopy_cond_stressed', process_id=ASSIMI_, owner_model_id=model_id)

  END FUNCTION Create_task_canopy_cond_stressed

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for assimilation task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "assimilation"
  !!
  FUNCTION Create_task_assimilation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_assimilation::return_ptr)
    CALL return_ptr%Construct(name='assimilation', process_id=ASSIMI_, owner_model_id=model_id)

  END FUNCTION Create_task_assimilation

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for assimi process
  !!
  !! @param[in,out] this      Instance of assimi process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_assimi_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_assimilation_scaling_factors(model_id))
    CALL this%Register_task(tsk_canopy_cond_unstressed(model_id))
    CALL this%Register_task(tsk_canopy_cond_stressed(model_id))
    CALL this%Register_task(tsk_assimilation(model_id))

  END SUBROUTINE Register_assimi_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation to calculate the height scaling factors of plant enzyms for assimilation
  !! Task "assimi" calculates the carbon assimilation of plants.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  ! The scaling factors represent the enzymatic distribution of plants in their canopy layers. This correlates
  ! with the amount of nitrogen in these layers (therefore called nitrogen_scaling_factors in JSBACH3).
  ! Eqs. (106), (107) in Knorr (1997)
  ! Note that this routine should only be applied to special vegetation classes
  SUBROUTINE update_assimilation_scaling_factors(tile, options)

    ! Use declarations
    !USE mo_orbit,                 ONLY: inquire_declination
     ! R: Returns declination calculated in last orbit call.
     !    Problem: FATAL ERROR in inquire_declination: Not Initialized
     !    This is because ICON did not call the calculation of orbital parameters in SUBROUTINE orbit_vsop87 or
     !    alternatively in SUBROUTINE orbit_kepler. Therefore this seems to be a problem solved with the development of ICON...
     !    Note: problem still exists in the 26.10.2015 version of icon-les-dev !
     !    However taken from mo_psrad_orbit works... perhaps precompiler command necessary here...
#ifdef __ICON__
    USE mo_psrad_orbit,       ONLY: inquire_declination
#else
    USE mo_orbit,             ONLY: inquire_declination
#endif
    USE mo_jsb_grid,          ONLY: Get_grid
    USE mo_assimi_process,    ONLY: calc_assimilation_scaling_factors

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER   :: model
    TYPE(t_jsb_grid),  POINTER   :: grid
    INTEGER                      :: iblk, ics, ice !, nc
    REAL(wp)                     :: lat(options%nc)
    REAL(wp)                     :: declination ! Solar Declination angle
    INTEGER                      :: i
    INTEGER                      :: ncanopy          ! number of canopy layers
    REAL(wp),          POINTER   :: canopy_bound_lai(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':assimilation_scaling_factors'

    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(PHENO_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real3D_onChunk :: assimi_scaling_fact_cl
    dsl4jsb_Real2D_onChunk :: lai

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    !nc      = options%nc

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(ASSIMI_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get process config
    dsl4jsb_Get_config(ASSIMI_)

    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(PHENO_)

    ! Get process variables
    ncanopy          =  dsl4jsb_Config(ASSIMI_)%ncanopy
    canopy_bound_lai(0:ncanopy) => dsl4jsb_Config(ASSIMI_)%canopy_bound_lai(0:ncanopy)

    grid => Get_grid(model%grid_id)
    lat(:)  = grid%lat(ics:ice, iblk)

    dsl4jsb_Get_var2D_onChunk(PHENO_, lai)                     ! in
    dsl4jsb_Get_var3D_onChunk(ASSIMI_, assimi_scaling_fact_cl) ! out

    ! ---------------------------
    ! Go

    CALL inquire_declination(declination) ! R: do not put this in elemental procedure as it would not be pure...
                                          ! e.g. declination = 0.006918_wp

    IF (dsl4jsb_Lctlib_param(NitrogenScalingFlag)) THEN
      DO i=1,ncanopy
        CALL calc_assimilation_scaling_factors( &
          & declination,                        & ! input
          & canopy_bound_lai(i-1),              & ! input
          & lai(:),                             & ! input
          & lat(:),                             & ! input
          & assimi_scaling_fact_cl(:,i)         & ! output
          & )
      END DO
    ELSE                                           ! to some vegetation classes nitrogen scaling is not applied
      assimi_scaling_fact_cl(:,:) = 1._wp
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_assimilation_scaling_factors

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "assimilation_scaling_factors"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_assimilation_scaling_factors(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(ASSIMI_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_assimilation_scaling_factors'

    INTEGER  :: iblk, ics, ice !, nc

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    !nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(ASSIMI_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(ASSIMI_, assimi_scaling_fact_cl, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_assimilation_scaling_factors

  ! ================================================================================================================================
  !>
  !! @brief Implementation of task canopy_cond_unstressed
  !! Task "canopy_cond_unstressed" calculates the carbon assimilation of plants and canopy conductance without water stress.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_canopy_cond_unstressed(tile, options)

    USE mo_assimi_process,      ONLY: calc_assimilation_waterunlimited

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER                :: model
    !REAL(wp)                                  :: dtime
    INTEGER                                   :: iblk, ics, ice !, nc
    INTEGER                                   :: icanopy

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_canopy_cond_unstressed'

    ! Declare process memories
    dsl4jsb_Def_memory(ASSIMI_)
    !dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(A2L_)

    dsl4jsb_Def_config(ASSIMI_)

    ! Declare pointers to variables in memory
    INTEGER :: ncanopy                 ! number of canopy layers

    !dsl4jsb_Real2D_onChunk :: C4flag ! Photosynthetic pathway (C3: 0; C4: 1)
                                      ! R: nur falls man C4flag doch als REAL haben wollte

    dsl4jsb_Real2D_onChunk :: t_air                 ! Atmosphere temperature (lowest layer) in Kelvin!
    dsl4jsb_Real2D_onChunk :: press_srf             ! Surface pressure
    dsl4jsb_Real2D_onChunk :: par_down_mol          ! Downward PAR flux in mol (photons)/(m^2 s)
    dsl4jsb_Real3D_onChunk :: apar_per_lai_cl       ! Input. Absorbed PAR of canopy layer.
    dsl4jsb_Real2D_onChunk :: canopy_cond_unlimited ! Canopy conductance no limitation
    !dsl4jsb_Real2D_onChunk :: canopy_cond_limited   ! Canopy conductance water limited

    dsl4jsb_Real2D_onChunk :: CO2_air_mol           ! CO2 particle mixing ratio [molCO2/molDryAir]

    dsl4jsb_Real3D_onChunk :: assimi_scaling_fact_cl
    dsl4jsb_Real3D_onChunk :: canopy_cond_cl
    dsl4jsb_Real3D_onChunk :: lai_cl                  ! Leaf area index [-]


    ! Get local variables from options argument
    !dtime   = options%dtime
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    !nc      = options%nc

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(ASSIMI_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(ASSIMI_)

    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    !dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(A2L_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(A2L_,      t_air)                   ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      press_srf)               ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      CO2_air_mol)             ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,      par_down_mol)            ! in
    dsl4jsb_Get_var3D_onChunk(RAD_,      apar_per_lai_cl)         ! in
    dsl4jsb_Get_var3D_onChunk(RAD_,      lai_cl)                  ! in

    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   assimi_scaling_fact_cl)  ! in
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   canopy_cond_cl)          ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   canopy_cond_unlimited)   ! out

    ncanopy = dsl4jsb_Config(ASSIMI_)%ncanopy

    DO icanopy=1,ncanopy
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Using jedi's Vc,max value to compute photosynthesis
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    IF (TRIM(model%config%usecase) .eq. 'jsbach_jedi' ) THEN
    
    
    
      CALL calc_assimilation_waterunlimited( &
        & dsl4jsb_Lctlib_param(C4Flag),      &
        & dsl4jsb_Lctlib_param(JeDiCarboxRate),  &
        & dsl4jsb_Lctlib_param(JeDiCarboxRate)*1.9_wp,  &  ! ECtransport rate = 1.9 Vc,max
        & t_air(:),                          & ! in
        & press_srf(:),                      & ! in
        & par_down_mol(:),                   & ! in
        & apar_per_lai_cl(:,icanopy),        & ! in
        & CO2_air_mol(:),                    & ! in
        & assimi_scaling_fact_cl(:,icanopy), & ! in
        & canopy_cond_cl(:,icanopy)          & ! out
        & )
    
    
    
    ELSE  ! Oringinal call for assimilation under water unlimited case in JB4
    
    write(*,*) "you are in a wrong routine for jedi assimilation computation!!!!!!!!"
!      CALL calc_assimilation_waterunlimited( &
!        & dsl4jsb_Lctlib_param(C4Flag),      &
!        & dsl4jsb_Lctlib_param(CarboxRate),  &
!        & dsl4jsb_Lctlib_param(ETransport),  &
!        & t_air(:),                          & ! in
!        & press_srf(:),                      & ! in
!        & par_down_mol(:),                   & ! in
!        & apar_per_lai_cl(:,icanopy),        & ! in
!        & CO2_air_mol(:),                    & ! in
!        & assimi_scaling_fact_cl(:,icanopy), & ! in
!        & canopy_cond_cl(:,icanopy)          & ! out
!        & )
        
    ENDIF
    
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    
    END DO

    ! Compute unlimited conductance for the whole canopy (all layers) considering the lai
    canopy_cond_unlimited(:) = MAX(1.e-20_wp,SUM(canopy_cond_cl(:,:) * lai_cl(:,:),DIM=2))

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_canopy_cond_unstressed

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "canopy_cond_unstressed"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_canopy_cond_unstressed(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(ASSIMI_)
    !dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_cond_unstressed'

    INTEGER  :: iblk, ics, ice !, nc

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    !nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(ASSIMI_)
    !dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(ASSIMI_, canopy_cond_cl,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, canopy_cond_unlimited, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_canopy_cond_unstressed

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! @brief Implementation of "update" for task "canopy_cond_stressed"
  !! Task "update_canopy_cond_stressed" calculates canopy conductance under water stress
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_canopy_cond_stressed(tile, options)
    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_hydro_process, ONLY: get_canopy_conductance
    USE mo_phy_schemes,   ONLY: qsat

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(ASSIMI_)

    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(HYDRO_)

    dsl4jsb_Real2D_onChunk :: canopy_cond_unlimited
    dsl4jsb_Real2D_onChunk :: canopy_cond_limited
    dsl4jsb_Real3D_onChunk :: canopy_cond_cl
    dsl4jsb_Real2D_onChunk :: water_stress
    dsl4jsb_Real2D_onChunk :: t
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: press_srf

    INTEGER :: iblk, ics, ice
    INTEGER :: icanopy, ncanopy

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_canopy_cond_stressed'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (.NOT. tile%Is_process_active(ASSIMI_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(ASSIMI_)

    ncanopy = dsl4jsb_Config(ASSIMI_)%ncanopy

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onChunk(A2L_,      q_air)               ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      press_srf)           ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t)                   ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    water_stress)        ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   canopy_cond_unlimited) ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   canopy_cond_limited) ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   canopy_cond_cl)       ! inout

    ! Compute (actual) canopy (stomatal) conductance under water stress.
    canopy_cond_limited(:) = get_canopy_conductance(canopy_cond_unlimited(:),  & ! in, unstressed canopy conductance
                                                    water_stress(:),           & ! in, water stress factor
                                                    q_air(:) >  qsat(t(:),press_srf(:)) & ! in, atmosphere saturated?
                                                    )
    ! Scale unlimited conductance per leaf area with the ratio of limited to unlimited conductance for
    ! whole canopy, ie. compute new limited conductance per leaf area and canopy layer
    DO icanopy=1,ncanopy
      canopy_cond_cl(:,icanopy) = canopy_cond_cl(:,icanopy) * canopy_cond_limited(:) / canopy_cond_unlimited(:)
    END DO

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_canopy_cond_stressed

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "canopy_cond_stressed"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_canopy_cond_stressed(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(ASSIMI_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_cond_stressed'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(ASSIMI_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(ASSIMI_, canopy_cond_limited, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, canopy_cond_cl,      weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_canopy_cond_stressed

  ! ================================================================================================================================
  !>
  !! @brief Implementation of assimilation for water-limited conditions
  !! Task "assimilation" calculates the carbon assimilation of plants.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_assimilation(tile, options)

    USE mo_assimi_process,      ONLY: calc_assimilation_waterlimited, calc_NPP_pot_rate

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER :: model
    INTEGER                    :: iblk, ics, ice !, nc
    REAL(wp)                   :: dtime
    INTEGER                    :: icanopy

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_assimilation'

    ! Declare pointers to process configuration and memory
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_config(ASSIMI_)

    INTEGER :: ncanopy                 ! number of canopy layers

    ! Declare pointers to variables in memory

    dsl4jsb_Real2D_onChunk :: t_air                 ! Atmosphere temperature (lowest layer) in Kelvin!
    dsl4jsb_Real2D_onChunk :: press_srf             ! Surface pressure
    dsl4jsb_Real2D_onChunk :: par_down_mol          ! Downward PAR flux in mol (photons)/(m^2 s)
    dsl4jsb_Real3D_onChunk :: apar_per_lai_cl       ! Input. Absorbed PAR of canopy layer.

    dsl4jsb_Real2D_onChunk :: CO2_air_mol           ! CO2 particle mixing ratio [molCO2/molDryAir]
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: veg_fract_correction  ! factor to convert from tile to canopy area
    dsl4jsb_Real2D_onChunk :: gross_assimilation_ca
    dsl4jsb_Real2D_onChunk :: gross_assimilation
    dsl4jsb_Real2D_onChunk :: dark_respiration_ca
    dsl4jsb_Real3D_onChunk :: dark_respiration_cl
    dsl4jsb_Real2D_onChunk :: dark_respiration
    dsl4jsb_Real2D_onChunk :: NPP_pot_rate_ca
    dsl4jsb_Real3D_onChunk :: assimi_scaling_fact_cl ! Input
    dsl4jsb_Real2D_onChunk :: fract_fpc_max          ! in
    dsl4jsb_Real3D_onChunk :: canopy_cond_cl         ! InOut for water limited case here.
    dsl4jsb_Real3D_onChunk :: gross_assimilation_cl  ! Output.
    dsl4jsb_Real3D_onChunk :: carbox_rate_max_cl     ! Output. Maximum carboxylation rate of canopy layer (=VCmax)
                                                         ! [Mol(CO2)/m^2(canopy)/s]. carbox_rate_max_leaf adapted to canopy layer
                                                         ! and actual temperature.
    dsl4jsb_Real3D_onChunk :: e_transport_rate_max_cl ! Output
    dsl4jsb_Real3D_onChunk :: carbox_rate_cl          ! Output
    dsl4jsb_Real3D_onChunk :: e_transport_rate_cl     ! Output
    dsl4jsb_Real3D_onChunk :: CO2_conc_leaf_cl        ! CO2 concentration inside leaf  [mol(CO2)/mol(Air)]:
                                                      ! Output only for the water limited case here
    dsl4jsb_Real3D_onChunk :: lai_cl                  ! Leaf area index [-]
    dsl4jsb_Real2D_onChunk :: day_NPP_sum

    ! Set local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    !nc      = options%nc
    dtime   = options%dtime

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(ASSIMI_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Set process config and memories
    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(ASSIMI_)

    ncanopy = dsl4jsb_Config(ASSIMI_)%ncanopy

    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(PHENO_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(A2L_,      t_air)                   ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      CO2_air_mol)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      press_srf)               ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,      par_down_mol)            ! in
    dsl4jsb_Get_var3D_onChunk(RAD_,      apar_per_lai_cl)         ! in
    dsl4jsb_Get_var3D_onChunk(RAD_,      lai_cl)                  ! in

    dsl4jsb_Get_var2D_onChunk(PHENO_,    veg_fract_correction)    ! in
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   assimi_scaling_fact_cl)  ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_,    fract_fpc_max )          ! in
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   canopy_cond_cl)          ! in

    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   gross_assimilation_ca)   ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   gross_assimilation)      ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   dark_respiration_ca)     ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   dark_respiration)        ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   NPP_pot_rate_ca)         ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   gross_assimilation_cl)   ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   dark_respiration_cl)     ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   carbox_rate_max_cl)      ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   e_transport_rate_max_cl) ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   e_transport_rate_cl)     ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   carbox_rate_cl)          ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   CO2_conc_leaf_cl)        ! out

    dsl4jsb_Get_var2D_onChunk(PHENO_,    day_NPP_sum)             ! inout

    ! ---------------------------
    ! Go

    ! R: Diese Berechung hier habe ich raus genommen, da CO2_conc_leaf_unlimited_acc weder für assimilation noch wo anderst
    !    verwendet wird. Wird nur ausgegeben. Braucht man wahrscheinlich nicht mehr, auch nicht im Output.
    ! bethy%CO2_conc_leaf_unlimited_acc(kidx0:kidx1,itile) = bethy%CO2_conc_leaf_unlimited_acc(kidx0:kidx1,itile) + &
    !                                                             CO2_conc_leaf_cl * dtime

    DO icanopy=1,ncanopy
    
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Using jedi's Vc,max value to compute photosynthesis
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    IF (TRIM(model%config%usecase) .eq. 'jsbach_jedi' ) THEN
        
    
      CALL calc_assimilation_waterlimited(    &
        & dsl4jsb_Lctlib_param(C4Flag),       &
        & dsl4jsb_Lctlib_param(JeDiCarboxRate),   &
        & dsl4jsb_Lctlib_param(JeDiCarboxRate)*1.9_wp,   &
        & t_air(:),                           & ! in
        & press_srf(:),                       & ! in
        & par_down_mol(:),                    & ! in
        & apar_per_lai_cl(:,icanopy),         & ! in
        & CO2_air_mol(:),                     & ! in
        & assimi_scaling_fact_cl(:,icanopy),  & ! in
        & CO2_conc_leaf_cl(:,icanopy),        & ! out only for the water limited case here
        & canopy_cond_cl(:,icanopy),          & ! in  only for the water limited case here
        & gross_assimilation_cl(:,icanopy),   & ! out
        & dark_respiration_cl(:,icanopy),     & ! out
        & carbox_rate_max_cl(:,icanopy),      & ! out
        & e_transport_rate_max_cl(:,icanopy), & ! out
        & carbox_rate_cl(:,icanopy),          & ! out
        & e_transport_rate_cl(:,icanopy)      & ! out
        & )
   
    ELSE
    write(*,*) "you are in a wrong routine for jedi assimilation computation!!!!!!!!"
!      CALL calc_assimilation_waterlimited(    &
!        & dsl4jsb_Lctlib_param(C4Flag),       &
!        & dsl4jsb_Lctlib_param(CarboxRate),   &
!        & dsl4jsb_Lctlib_param(ETransport),   &
!        & t_air(:),                           & ! in
!        & press_srf(:),                       & ! in
!        & par_down_mol(:),                    & ! in
!        & apar_per_lai_cl(:,icanopy),         & ! in
!        & CO2_air_mol(:),                     & ! in
!        & assimi_scaling_fact_cl(:,icanopy),  & ! in
!        & CO2_conc_leaf_cl(:,icanopy),        & ! out only for the water limited case here
!        & canopy_cond_cl(:,icanopy),          & ! in  only for the water limited case here
!        & gross_assimilation_cl(:,icanopy),   & ! out
!        & dark_respiration_cl(:,icanopy),     & ! out
!        & carbox_rate_max_cl(:,icanopy),      & ! out
!        & e_transport_rate_max_cl(:,icanopy), & ! out
!        & carbox_rate_cl(:,icanopy),          & ! out
!        & e_transport_rate_cl(:,icanopy)      & ! out
!        & )
        
        
    ENDIF
    
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
 
        
        
    END DO

    ! Same for gross_assimilation_ca
    ! R: dies ist eine alternative Schreibweise die verwendet wurde da der PGI Kompiler probleme gemacht hat.
    !gross_assimilation(:) = 0._wp
    !DO icanopy=1,ncanopy
    !   gross_assimilation(:) = gross_assimilation(:) + gross_assimilation_cl(:,icanopy) * lai_cl(:,icanopy)
    !END DO

    ! convert from leaf area of layer to canopy ground area
    gross_assimilation_ca(:) = SUM(gross_assimilation_cl(:,:) * lai_cl(:,:),DIM=2)
    ! convert from canopy area to tile box area
    gross_assimilation(:) = gross_assimilation_ca(:) * veg_fract_correction(:) * fract_fpc_max(:)

    ! Same for dark_respiration_ca
    ! R: dies ist eine alternative Schreibweise die verwendet wurde da der PGI Kompiler probleme gemacht hat.
    !dark_respiration_ca(:) = 0._wp
    !DO icanopy=1,ncanopy
    !   dark_respiration_ca(:) =  dark_respiration_ca(:) + dark_respiration_cl(:,icanopy) * lai_cl(:,icanopy)
    !END DO

    dark_respiration_ca(:) = SUM(dark_respiration_cl(:,:) * lai_cl(:,:),DIM=2)
    ! convert from canopy area to tile box area
    dark_respiration(:) =  dark_respiration_ca(:)  * veg_fract_correction(:) * fract_fpc_max(:)

    CALL calc_NPP_pot_rate(gross_assimilation_ca, dark_respiration_ca, NPP_pot_rate_ca)

    ! Accumulate NPP_pot_rate_ca. day_NPP_sum is divided by the day length and reset to zero in
    ! update_phenology_logrop at the beginning of each day.
    day_NPP_sum(:) = day_NPP_sum(:) + NPP_pot_rate_ca(:) * dtime

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_assimilation

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "assimilation"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_assimilation(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(ASSIMI_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice !, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_assimilation'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    !nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(ASSIMI_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(ASSIMI_, canopy_cond_cl,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, gross_assimilation_ca,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, gross_assimilation,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, dark_respiration_ca,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, dark_respiration,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, NPP_pot_rate_ca,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, CO2_conc_leaf_cl,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, carbox_rate_cl,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, carbox_rate_max_cl,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, e_transport_rate_cl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, e_transport_rate_max_cl,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_assimilation


#endif
END MODULE mo_assimi_interface
