!
! This file was automatically generated from JSBACH source file:
!    src/assimilation/mo_assimi_interface.f90
! with md5 hash: 21d1829089bd4a78675ed9c6ac4f6010
!
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
  USE mo_jsb_process_class, ONLY: assimi_, hydro_, pheno_, rad_, a2l_, seb_

  ! Use of process configurations (t_assimi_config)
  USE mo_assimi_config_class, ONLY: t_assimi_config

  ! Use of process memories (t_assimi_memory)
  USE mo_a2l_memory_class, ONLY: t_a2l_memory
  USE mo_assimi_memory_class, ONLY: t_assimi_memory
  USE mo_hydro_memory_class, ONLY: t_hydro_memory
  USE mo_pheno_memory_class, ONLY: t_pheno_memory
  USE mo_rad_memory_class, ONLY: t_rad_memory
  USE mo_seb_memory_class, ONLY: t_seb_memory

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
    TYPE(t_assimi_config), POINTER :: assimi__conf
    TYPE(t_assimi_memory), POINTER :: assimi__mem
    TYPE(t_pheno_memory), POINTER :: pheno__mem

    ! Declare pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:,:) :: assimi_scaling_fact_cl
    REAL(wp), POINTER, DIMENSION(:) :: lai

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
    SELECT TYPE (cc => model%processes(assimi_)%p%config); TYPE IS (t_assimi_config) ; assimi__conf => cc ; END SELECT

    ! Get process memories
    SELECT TYPE (mm => tile%mem(assimi_)%p); TYPE IS (t_assimi_memory); assimi__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT

    ! Get process variables
    ncanopy          =  assimi__conf%ncanopy
    canopy_bound_lai(0:ncanopy) => assimi__conf%canopy_bound_lai(0:ncanopy)

    grid => Get_grid(model%grid_id)
    lat(:)  = grid%lat(ics:ice, iblk)

    lai => pheno__mem%lai%ptr(ics:ice,iblk)                     ! in
    assimi_scaling_fact_cl => assimi__mem%assimi_scaling_fact_cl%ptr(ics:ice,:,iblk) ! out

    ! ---------------------------
    ! Go

    !$claw acc data create(lat, canopy_bound_lai)
    !$claw acc update device(lat, canopy_bound_lai) ! LOCAL VARIABLE

    CALL inquire_declination(declination) ! R: do not put this in elemental procedure as it would not be pure...
                                          ! e.g. declination = 0.006918_wp

    IF (model%lctlib(tile%lcts(1)%lib_id)%NitrogenScalingFlag) THEN
      DO i=1,ncanopy
        !$claw sca forward update savepoint=ser_jsb_assim_calc_assimilation_scaling_factors_1 tile_name=tile%name canopy_lvl=i
        CALL calc_assimilation_scaling_factors( &
          & declination,                        & ! input
          & canopy_bound_lai(i-1),              & ! input
          & lai(:),                             & ! input
          & lat(:),                             & ! input
          & assimi_scaling_fact_cl(:,i)         & ! output
          & )
      END DO
    ELSE                                           ! to some vegetation classes nitrogen scaling is not applied
      !$claw expand parallel update savepoint=ser_jsb_assim_update_assimilation_scaling_factors_expand_0 tile_name=tile%name 
      assimi_scaling_fact_cl(:,:) = 1._wp
    END IF

    !$claw acc end data

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

    TYPE(t_assimi_memory), POINTER :: assimi__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_assimilation_scaling_factors'

    INTEGER  :: iblk, ics, ice !, nc

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    !nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(assimi_)%p); TYPE IS (t_assimi_memory); assimi__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, assimi__mem%assimi_scaling_fact_cl,ics,ice,iblk,routine)

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
    TYPE(t_assimi_memory), POINTER :: assimi__mem
    !TYPE(t_hydro_memory), POINTER :: hydro__mem
    TYPE(t_rad_memory), POINTER :: rad__mem
    TYPE(t_a2l_memory), POINTER :: a2l__mem

    TYPE(t_assimi_config), POINTER :: assimi__conf

    ! Declare pointers to variables in memory
    INTEGER :: ncanopy                 ! number of canopy layers

    !REAL(wp), POINTER, DIMENSION(:) :: C4flag ! Photosynthetic pathway (C3: 0; C4: 1)
                                      ! R: nur falls man C4flag doch als REAL haben wollte

    REAL(wp), POINTER, DIMENSION(:) :: t_air                 ! Atmosphere temperature (lowest layer) in Kelvin!
    REAL(wp), POINTER, DIMENSION(:) :: press_srf             ! Surface pressure
    REAL(wp), POINTER, DIMENSION(:) :: par_down_mol          ! Downward PAR flux in mol (photons)/(m^2 s)
    REAL(wp), POINTER, DIMENSION(:,:) :: apar_per_lai_cl       ! Input. Absorbed PAR of canopy layer.
    REAL(wp), POINTER, DIMENSION(:) :: canopy_cond_unlimited ! Canopy conductance no limitation
    !REAL(wp), POINTER, DIMENSION(:) :: canopy_cond_limited   ! Canopy conductance water limited

    REAL(wp), POINTER, DIMENSION(:) :: co2_air_mol           ! CO2 particle mixing ratio [molCO2/molDryAir]

    REAL(wp), POINTER, DIMENSION(:,:) :: assimi_scaling_fact_cl
    REAL(wp), POINTER, DIMENSION(:,:) :: canopy_cond_cl
    REAL(wp), POINTER, DIMENSION(:,:) :: lai_cl                  ! Leaf area index [-]


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
    SELECT TYPE (cc => model%processes(assimi_)%p%config); TYPE IS (t_assimi_config) ; assimi__conf => cc ; END SELECT

    ! Get process memories
    SELECT TYPE (mm => tile%mem(assimi_)%p); TYPE IS (t_assimi_memory); assimi__mem => mm ; END SELECT
    !SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(rad_)%p); TYPE IS (t_rad_memory); rad__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT

    ! Get process variables
    t_air => a2l__mem%t_air%ptr(ics:ice,iblk)                   ! in
    press_srf => a2l__mem%press_srf%ptr(ics:ice,iblk)               ! in
    co2_air_mol => a2l__mem%co2_air_mol%ptr(ics:ice,iblk)             ! in
    par_down_mol => rad__mem%par_down_mol%ptr(ics:ice,iblk)            ! in
    apar_per_lai_cl => rad__mem%apar_per_lai_cl%ptr(ics:ice,:,iblk)         ! in
    lai_cl => rad__mem%lai_cl%ptr(ics:ice,:,iblk)                  ! in

    assimi_scaling_fact_cl => assimi__mem%assimi_scaling_fact_cl%ptr(ics:ice,:,iblk)  ! in
    canopy_cond_cl => assimi__mem%canopy_cond_cl%ptr(ics:ice,:,iblk)          ! out
    canopy_cond_unlimited => assimi__mem%canopy_cond_unlimited%ptr(ics:ice,iblk)   ! out

    ncanopy = assimi__conf%ncanopy

    DO icanopy=1,ncanopy
      !$claw sca forward update savepoint=ser_jsb_assim_calc_assimilation_waterunlimited_1 tile_name=tile%name canopy_lvl=icanopy
      CALL calc_assimilation_waterunlimited( &
        & model%lctlib(tile%lcts(1)%lib_id)%C4Flag,      &
        & model%lctlib(tile%lcts(1)%lib_id)%CarboxRate,  &
        & model%lctlib(tile%lcts(1)%lib_id)%ETransport,  &
        & t_air(:),                          & ! in
        & press_srf(:),                      & ! in
        & par_down_mol(:),                   & ! in
        & apar_per_lai_cl(:,icanopy),        & ! in
        & CO2_air_mol(:),                    & ! in
        & assimi_scaling_fact_cl(:,icanopy), & ! in
        & canopy_cond_cl(:,icanopy)          & ! out
        & )
    END DO

    ! Compute unlimited conductance for the whole canopy (all layers) considering the lai
    !$claw expand parallel update savepoint=ser_jsb_assim_update_canopy_cond_unstressed_expand_0 tile_name=tile%name
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

    TYPE(t_assimi_memory), POINTER :: assimi__mem
    !TYPE(t_hydro_memory), POINTER :: hydro__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_cond_unstressed'

    INTEGER  :: iblk, ics, ice !, nc

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    !nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(assimi_)%p); TYPE IS (t_assimi_memory); assimi__mem => mm ; END SELECT
    !SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, assimi__mem%canopy_cond_cl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, assimi__mem%canopy_cond_unlimited,ics,ice,iblk,routine)

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

    TYPE(t_assimi_config), POINTER :: assimi__conf

    TYPE(t_assimi_memory), POINTER :: assimi__mem
    TYPE(t_seb_memory), POINTER :: seb__mem
    TYPE(t_a2l_memory), POINTER :: a2l__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem

    REAL(wp), POINTER, DIMENSION(:) :: canopy_cond_unlimited
    REAL(wp), POINTER, DIMENSION(:) :: canopy_cond_limited
    REAL(wp), POINTER, DIMENSION(:,:) :: canopy_cond_cl
    REAL(wp), POINTER, DIMENSION(:) :: water_stress
    REAL(wp), POINTER, DIMENSION(:) :: t
    REAL(wp), POINTER, DIMENSION(:) :: q_air
    REAL(wp), POINTER, DIMENSION(:) :: press_srf

    INTEGER :: iblk, ics, ice, nc
    INTEGER :: icanopy, ncanopy
    ! LOGICAL, DIMENSION(options%nc) :: q_air_gt_qsat_tmp
    LOGICAL, ALLOCATABLE :: q_air_gt_qsat_tmp(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_canopy_cond_stressed'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_active(ASSIMI_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    SELECT TYPE (cc => model%processes(assimi_)%p%config); TYPE IS (t_assimi_config) ; assimi__conf => cc ; END SELECT

    ncanopy = assimi__conf%ncanopy

    ! Get reference to variables for current block
    !
    SELECT TYPE (mm => tile%mem(assimi_)%p); TYPE IS (t_assimi_memory); assimi__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    q_air => a2l__mem%q_air%ptr(ics:ice,iblk)               ! in
    press_srf => a2l__mem%press_srf%ptr(ics:ice,iblk)           ! in
    t => seb__mem%t%ptr(ics:ice,iblk)                   ! in
    water_stress => hydro__mem%water_stress%ptr(ics:ice,iblk)        ! in
    canopy_cond_unlimited => assimi__mem%canopy_cond_unlimited%ptr(ics:ice,iblk) ! in
    canopy_cond_limited => assimi__mem%canopy_cond_limited%ptr(ics:ice,iblk) ! out
    canopy_cond_cl => assimi__mem%canopy_cond_cl%ptr(ics:ice,:,iblk)       ! inout

    ! Compute (actual) canopy (stomatal) conductance under water stress.
    !$cla$w acc data create(q_air_gt_qsat_tmp)
    ALLOCATE(q_air_gt_qsat_tmp(nc))
    !$claw acc enter data create(q_air_gt_qsat_tmp)
    !$claw expand parallel update savepoint=ser_jsb_assimi_update_canopy_cond_stressed_expand_0 tile_name=tile%name
    q_air_gt_qsat_tmp(:) = q_air(:) > qsat(t(:),press_srf(:))      
    canopy_cond_limited(:) = get_canopy_conductance(canopy_cond_unlimited(:),  & ! in, unstressed canopy conductance
                                                    water_stress(:),           & ! in, water stress factor
                                                    q_air_gt_qsat_tmp(:)       & ! in, atmosphere saturated?
                                                    )
    !$claw end expand
    !$cla$w acc end data
    !$claw acc exit data delete(q_air_gt_qsat_tmp)
    DEALLOCATE(q_air_gt_qsat_tmp)

    ! Scale unlimited conductance per leaf area with the ratio of limited to unlimited conductance for
    ! whole canopy, ie. compute new limited conductance per leaf area and canopy layer
    DO icanopy=1,ncanopy
      !$claw expand parallel update savepoint=ser_jsb_assim_update_canopy_cond_stressed_expand_1 tile_name=tile%name canopy_lvl=icanopy
      canopy_cond_cl(:,icanopy) = canopy_cond_cl(:,icanopy) * canopy_cond_limited(:) / canopy_cond_unlimited(:)
    END DO
    !$claw end expand

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

    TYPE(t_assimi_memory), POINTER :: assimi__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_cond_stressed'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(assimi_)%p); TYPE IS (t_assimi_memory); assimi__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, assimi__mem%canopy_cond_limited,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, assimi__mem%canopy_cond_cl,ics,ice,iblk,routine)

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
    TYPE(t_assimi_memory), POINTER :: assimi__mem
    TYPE(t_rad_memory), POINTER :: rad__mem
    TYPE(t_a2l_memory), POINTER :: a2l__mem
    TYPE(t_pheno_memory), POINTER :: pheno__mem
    TYPE(t_assimi_config), POINTER :: assimi__conf

    INTEGER :: ncanopy                 ! number of canopy layers

    ! Declare pointers to variables in memory

    REAL(wp), POINTER, DIMENSION(:) :: t_air                 ! Atmosphere temperature (lowest layer) in Kelvin!
    REAL(wp), POINTER, DIMENSION(:) :: press_srf             ! Surface pressure
    REAL(wp), POINTER, DIMENSION(:) :: par_down_mol          ! Downward PAR flux in mol (photons)/(m^2 s)
    REAL(wp), POINTER, DIMENSION(:,:) :: apar_per_lai_cl       ! Input. Absorbed PAR of canopy layer.

    REAL(wp), POINTER, DIMENSION(:) :: co2_air_mol           ! CO2 particle mixing ratio [molCO2/molDryAir]

    REAL(wp), POINTER, DIMENSION(:) :: veg_fract_correction  ! factor to convert from tile to canopy area
    REAL(wp), POINTER, DIMENSION(:) :: gross_assimilation_ca
    REAL(wp), POINTER, DIMENSION(:) :: gross_assimilation
    REAL(wp), POINTER, DIMENSION(:) :: dark_respiration_ca
    REAL(wp), POINTER, DIMENSION(:,:) :: dark_respiration_cl
    REAL(wp), POINTER, DIMENSION(:) :: dark_respiration
    REAL(wp), POINTER, DIMENSION(:) :: npp_pot_rate_ca
    REAL(wp), POINTER, DIMENSION(:,:) :: assimi_scaling_fact_cl ! Input
    REAL(wp), POINTER, DIMENSION(:) :: fract_fpc_max          ! in
    REAL(wp), POINTER, DIMENSION(:,:) :: canopy_cond_cl         ! InOut for water limited case here.
    REAL(wp), POINTER, DIMENSION(:,:) :: gross_assimilation_cl  ! Output.
    REAL(wp), POINTER, DIMENSION(:,:) :: carbox_rate_max_cl     ! Output. Maximum carboxylation rate of canopy layer (=VCmax)
                                                         ! [Mol(CO2)/m^2(canopy)/s]. carbox_rate_max_leaf adapted to canopy layer
                                                         ! and actual temperature.
    REAL(wp), POINTER, DIMENSION(:,:) :: e_transport_rate_max_cl ! Output
    REAL(wp), POINTER, DIMENSION(:,:) :: carbox_rate_cl          ! Output
    REAL(wp), POINTER, DIMENSION(:,:) :: e_transport_rate_cl     ! Output
    REAL(wp), POINTER, DIMENSION(:,:) :: co2_conc_leaf_cl        ! CO2 concentration inside leaf  [mol(CO2)/mol(Air)]:
                                                      ! Output only for the water limited case here
    REAL(wp), POINTER, DIMENSION(:,:) :: lai_cl                  ! Leaf area index [-]
    REAL(wp), POINTER, DIMENSION(:) :: day_npp_sum

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

    SELECT TYPE (cc => model%processes(assimi_)%p%config); TYPE IS (t_assimi_config) ; assimi__conf => cc ; END SELECT

    ncanopy = assimi__conf%ncanopy

    ! Get process memories
    SELECT TYPE (mm => tile%mem(assimi_)%p); TYPE IS (t_assimi_memory); assimi__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(rad_)%p); TYPE IS (t_rad_memory); rad__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT

    ! Get process variables
    t_air => a2l__mem%t_air%ptr(ics:ice,iblk)                   ! in
    co2_air_mol => a2l__mem%co2_air_mol%ptr(ics:ice,iblk)             ! in
    press_srf => a2l__mem%press_srf%ptr(ics:ice,iblk)               ! in
    par_down_mol => rad__mem%par_down_mol%ptr(ics:ice,iblk)            ! in
    apar_per_lai_cl => rad__mem%apar_per_lai_cl%ptr(ics:ice,:,iblk)         ! in
    lai_cl => rad__mem%lai_cl%ptr(ics:ice,:,iblk)                  ! in

    veg_fract_correction => pheno__mem%veg_fract_correction%ptr(ics:ice,iblk)    ! in
    assimi_scaling_fact_cl => assimi__mem%assimi_scaling_fact_cl%ptr(ics:ice,:,iblk)  ! in
    fract_fpc_max => pheno__mem%fract_fpc_max%ptr(ics:ice,iblk)          ! in
    canopy_cond_cl => assimi__mem%canopy_cond_cl%ptr(ics:ice,:,iblk)          ! in

    gross_assimilation_ca => assimi__mem%gross_assimilation_ca%ptr(ics:ice,iblk)   ! out
    gross_assimilation => assimi__mem%gross_assimilation%ptr(ics:ice,iblk)      ! out
    dark_respiration_ca => assimi__mem%dark_respiration_ca%ptr(ics:ice,iblk)     ! out
    dark_respiration => assimi__mem%dark_respiration%ptr(ics:ice,iblk)        ! out
    npp_pot_rate_ca => assimi__mem%npp_pot_rate_ca%ptr(ics:ice,iblk)         ! out
    gross_assimilation_cl => assimi__mem%gross_assimilation_cl%ptr(ics:ice,:,iblk)   ! out
    dark_respiration_cl => assimi__mem%dark_respiration_cl%ptr(ics:ice,:,iblk)     ! out
    carbox_rate_max_cl => assimi__mem%carbox_rate_max_cl%ptr(ics:ice,:,iblk)      ! out
    e_transport_rate_max_cl => assimi__mem%e_transport_rate_max_cl%ptr(ics:ice,:,iblk) ! out
    e_transport_rate_cl => assimi__mem%e_transport_rate_cl%ptr(ics:ice,:,iblk)     ! out
    carbox_rate_cl => assimi__mem%carbox_rate_cl%ptr(ics:ice,:,iblk)          ! out
    co2_conc_leaf_cl => assimi__mem%co2_conc_leaf_cl%ptr(ics:ice,:,iblk)        ! out

    day_npp_sum => pheno__mem%day_npp_sum%ptr(ics:ice,iblk)             ! inout

    ! ---------------------------
    ! Go

    ! R: Diese Berechung hier habe ich raus genommen, da CO2_conc_leaf_unlimited_acc weder für assimilation noch wo anderst
    !    verwendet wird. Wird nur ausgegeben. Braucht man wahrscheinlich nicht mehr, auch nicht im Output.
    ! bethy%CO2_conc_leaf_unlimited_acc(kidx0:kidx1,itile) = bethy%CO2_conc_leaf_unlimited_acc(kidx0:kidx1,itile) + &
    !                                                             CO2_conc_leaf_cl * dtime

    DO icanopy=1,ncanopy
      !$claw sca forward update savepoint=ser_jsb_assim_calc_assimilation_waterlimited_1 tile_name=tile%name canopy_lvl=icanopy
      CALL calc_assimilation_waterlimited(    &
        & model%lctlib(tile%lcts(1)%lib_id)%C4Flag,       &
        & model%lctlib(tile%lcts(1)%lib_id)%CarboxRate,   &
        & model%lctlib(tile%lcts(1)%lib_id)%ETransport,   &
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
    END DO

    !$claw expand parallel update savepoint=ser_jsb_assim_update_assimilation_expand_0 tile_name=tile%name

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
    !$claw end expand 

    !$claw sca forward update savepoint=ser_jsb_assim_calc_npp_pot_rate_1 tile_name=tile%name
    CALL calc_NPP_pot_rate(gross_assimilation_ca(:), dark_respiration_ca(:), NPP_pot_rate_ca(:))

    ! Accumulate NPP_pot_rate_ca. day_NPP_sum is divided by the day length and reset to zero in
    ! update_phenology_logrop at the beginning of each day.
    !$claw expand parallel update savepoint=ser_jsb_assim_day_npp_sum_expand_0 tile_name=tile%name
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

    TYPE(t_assimi_memory), POINTER :: assimi__mem

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

    SELECT TYPE (mm => tile%mem(assimi_)%p); TYPE IS (t_assimi_memory); assimi__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, assimi__mem%canopy_cond_cl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, assimi__mem%gross_assimilation_ca,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, assimi__mem%gross_assimilation,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, assimi__mem%dark_respiration_ca,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, assimi__mem%dark_respiration,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, assimi__mem%npp_pot_rate_ca,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, assimi__mem%co2_conc_leaf_cl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, assimi__mem%carbox_rate_cl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, assimi__mem%carbox_rate_max_cl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, assimi__mem%e_transport_rate_cl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, assimi__mem%e_transport_rate_max_cl,ics,ice,iblk,routine)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_assimilation


#endif
END MODULE mo_assimi_interface
