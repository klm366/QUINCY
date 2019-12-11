!
! This file was automatically generated from JSBACH source file:
!    src/hydrology/mo_hydro_interface.f90
! with md5 hash: 23fceb08b0389ac1bf1b05d1da968cf9
!
!>
!! @brief Contains the interfaces to the hydro processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-29)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_hydro_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_jsb_time,           ONLY: is_time_experiment_start
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  USE mo_jsb_process_class, ONLY: seb_, sse_, turb_, pheno_, a2l_, hydro_

  ! Use of process configurations
  USE mo_seb_config_class, ONLY: t_seb_config
  USE mo_sse_config_class, ONLY: t_sse_config
  USE mo_hydro_config_class, ONLY: t_hydro_config

  ! Use of process memories
  USE mo_a2l_memory_class, ONLY: t_a2l_memory
  USE mo_hydro_memory_class, ONLY: t_hydro_memory
  USE mo_seb_memory_class, ONLY: t_seb_memory
  USE mo_sse_memory_class, ONLY: t_sse_memory
  USE mo_turb_memory_class, ONLY: t_turb_memory
  USE mo_pheno_memory_class, ONLY: t_pheno_memory

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_hydro_tasks !,t_hydro_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_interface'

  !> Type definition for surface_hydrology
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_hydrology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_hydrology
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_hydrology
  END TYPE tsk_surface_hydrology

  !> Constructor interface for surface_hydrology
  INTERFACE tsk_surface_hydrology
    PROCEDURE Create_task_surface_hydrology
  END INTERFACE tsk_surface_hydrology

  !> Type definition for soil_properties
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_properties
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_properties
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_properties
  END TYPE tsk_soil_properties

  !> Constructor interface for soil_properties
  INTERFACE tsk_soil_properties
    PROCEDURE Create_task_soil_properties
  END INTERFACE tsk_soil_properties

  !> Type definition for soil_hydrology
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_hydrology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_hydrology
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_hydrology
  END TYPE tsk_soil_hydrology

  !> Constructor interface for soil_hydrology
  INTERFACE tsk_soil_hydrology
    PROCEDURE Create_task_soil_hydrology
  END INTERFACE tsk_soil_hydrology

  !> Type definition for evaporation task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_evaporation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_evaporation
    PROCEDURE, NOPASS :: Aggregate => aggregate_evaporation
  END TYPE tsk_evaporation

  !> Constructor interface for evaporation task
  INTERFACE tsk_evaporation
    PROCEDURE Create_task_evaporation
  END INTERFACE tsk_evaporation

  !> Type definition for canopy_cond_unstressed task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_cond_unstressed
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_cond_unstressed
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_cond_unstressed
  END TYPE tsk_canopy_cond_unstressed

  !> Constructor interface for canopy_cond_unstressed task
  INTERFACE tsk_canopy_cond_unstressed
    PROCEDURE Create_task_canopy_cond_unstressed
  END INTERFACE tsk_canopy_cond_unstressed

  !> Type definition for water_stress task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_water_stress
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_water_stress
    PROCEDURE, NOPASS :: Aggregate => aggregate_water_stress
  END TYPE tsk_water_stress

  !> Constructor interface for water_stress task
  INTERFACE tsk_water_stress
    PROCEDURE Create_task_water_stress
  END INTERFACE tsk_water_stress

  !> Type definition for canopy_cond_stressed task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_cond_stressed
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_cond_stressed
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_cond_stressed
  END TYPE tsk_canopy_cond_stressed

  !> Constructor interface for canopy_cond_stressed task
  INTERFACE tsk_canopy_cond_stressed
    PROCEDURE Create_task_canopy_cond_stressed
  END INTERFACE tsk_canopy_cond_stressed

  !> Type definition for snow_hydrology task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_snow_and_ice_hydrology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_snow_and_ice_hydrology
    PROCEDURE, NOPASS :: Aggregate => aggregate_snow_and_ice_hydrology
  END TYPE tsk_snow_and_ice_hydrology

  !> Constructor interface for snow_hydrology task
  INTERFACE tsk_snow_and_ice_hydrology
    PROCEDURE Create_task_snow_and_ice_hydrology
  END INTERFACE tsk_snow_and_ice_hydrology

  !> Type definition for snow_and_skin_fraction task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_snow_and_skin_fraction
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_snow_and_skin_fraction    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_snow_and_skin_fraction !< Aggregates computed task variables
  END TYPE tsk_snow_and_skin_fraction

  !> Constructor interface for snow_and_skin_fraction task
  INTERFACE tsk_snow_and_skin_fraction
    PROCEDURE Create_task_snow_and_skin_fraction                       !< Constructor function for task
  END INTERFACE tsk_snow_and_skin_fraction

  !> Type definition for water_balance task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_water_balance
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_water_balance
    PROCEDURE, NOPASS :: Aggregate => aggregate_water_balance
  END TYPE tsk_water_balance

  !> Constructor interface for water balance task
  INTERFACE tsk_water_balance
    PROCEDURE Create_task_water_balance
  END INTERFACE tsk_water_balance

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_hydrology task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "surface_hydrology"
  !!
  FUNCTION Create_task_surface_hydrology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_hydrology::return_ptr)
    CALL return_ptr%Construct(name='surface_hydrology', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_hydrology

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for soil_properties task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "soil_properties"
  !!
  FUNCTION Create_task_soil_properties(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_soil_properties::return_ptr)
    CALL return_ptr%Construct(name='soil_properties', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_soil_properties

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for soil_hydrology task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "soil_hydrology"
  !!
  FUNCTION Create_task_soil_hydrology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_soil_hydrology::return_ptr)
    CALL return_ptr%Construct(name='soil_hydrology', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_soil_hydrology

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
    CALL return_ptr%Construct(name='canopy_cond_unstressed', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_canopy_cond_unstressed

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for water_stress task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "water_stress"
  !!
  FUNCTION Create_task_water_stress(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_water_stress::return_ptr)
    CALL return_ptr%Construct(name='water_stress', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_water_stress

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
    CALL return_ptr%Construct(name='canopy_cond_stressed', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_canopy_cond_stressed

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for evaporation task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "evaporation"
  !!
  FUNCTION Create_task_evaporation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_evaporation::return_ptr)
    CALL return_ptr%Construct(name='evaporation', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_evaporation

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for snow hydrology task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "snow_hydrology"
  !!
  FUNCTION Create_task_snow_and_ice_hydrology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_snow_and_ice_hydrology::return_ptr)
    CALL return_ptr%Construct(name='snow_ice_hydrology', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_snow_and_ice_hydrology

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for snow_and_skin_fraction task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "tsk_snow_and_skin_fraction"
  !!
  FUNCTION Create_task_snow_and_skin_fraction(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_snow_and_skin_fraction::return_ptr)
    CALL return_ptr%Construct(name='snow_and_skin_fraction', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_snow_and_skin_fraction

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for water_balance task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "water_balance"
  !!
  FUNCTION Create_task_water_balance(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_water_balance::return_ptr)
    CALL return_ptr%Construct(name='water_balance', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_water_balance

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for hydrology process
  !!
  !! @param[in]     model_id  Model id
  !! @param[in,out] this      Instance of HYDRO_ process class
  !!
  SUBROUTINE Register_hydro_tasks(process, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,              INTENT(in)    :: model_id

    CALL process%Register_task( tsk_surface_hydrology      (model_id))
    CALL process%Register_task( tsk_soil_properties        (model_id))
    CALL process%Register_task( tsk_soil_hydrology         (model_id))
    CALL process%Register_task( tsk_water_stress           (model_id))
    CALL process%Register_task( tsk_canopy_cond_unstressed (model_id))
    CALL process%Register_task( tsk_canopy_cond_stressed   (model_id))
    CALL process%Register_task( tsk_evaporation            (model_id))
    !CALL process%Register_task( tsk_snow_and_ice_hydrology (model_id))
    CALL process%Register_task( tsk_snow_and_skin_fraction (model_id))
    CALL process%Register_task( tsk_water_balance          (model_id))

  END SUBROUTINE Register_hydro_tasks
  !
  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "surface_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_surface_hydrology(tile, options)

    USE mo_hydro_process,          ONLY: calc_surface_hydrology_land, calc_surface_hydrology_glacier
    USE mo_jsb_physical_constants, ONLY: dens_snow

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_hydro_config), POINTER :: hydro__conf
    TYPE(t_sse_config), POINTER :: sse__conf
    TYPE(t_hydro_memory), POINTER :: hydro__mem
    TYPE(t_sse_memory), POINTER :: sse__mem
    TYPE(t_seb_memory), POINTER :: seb__mem
    TYPE(t_pheno_memory), POINTER :: pheno__mem
    TYPE(t_a2l_memory), POINTER :: a2l__mem

    ! Pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: t_unfilt
    REAL(wp), POINTER, DIMENSION(:) :: q_snocpymlt
    REAL(wp), POINTER, DIMENSION(:) :: t_air
    REAL(wp), POINTER, DIMENSION(:) :: wind_10m
    REAL(wp), POINTER, DIMENSION(:) :: rain
    REAL(wp), POINTER, DIMENSION(:) :: snow
    REAL(wp), POINTER, DIMENSION(:) :: w_skin
    REAL(wp), POINTER, DIMENSION(:) :: fract_water
    REAL(wp), POINTER, DIMENSION(:) :: w_snow_soil
    REAL(wp), POINTER, DIMENSION(:) :: snow_soil_dens
    REAL(wp), POINTER, DIMENSION(:) :: w_snow
    REAL(wp), POINTER, DIMENSION(:) :: snow_accum
    REAL(wp), POINTER, DIMENSION(:) :: fract_snow
    REAL(wp), POINTER, DIMENSION(:) :: w_snow_can
    REAL(wp), POINTER, DIMENSION(:) :: evapotrans_soil
    REAL(wp), POINTER, DIMENSION(:) :: evapo_skin
    REAL(wp), POINTER, DIMENSION(:) :: evapo_deficit
    REAL(wp), POINTER, DIMENSION(:) :: lai
    REAL(wp), POINTER, DIMENSION(:) :: fract_fpc_max
    REAL(wp), POINTER, DIMENSION(:) :: hcap_grnd
    REAL(wp), POINTER, DIMENSION(:) :: evapotrans
    REAL(wp), POINTER, DIMENSION(:) :: trans
    REAL(wp), POINTER, DIMENSION(:) :: evapopot
    REAL(wp), POINTER, DIMENSION(:) :: evapo_snow
    REAL(wp), POINTER, DIMENSION(:) :: snowmelt
    REAL(wp), POINTER, DIMENSION(:) :: runoff
    REAL(wp), POINTER, DIMENSION(:) :: drainage
    REAL(wp), POINTER, DIMENSION(:) :: runoff_glac
    REAL(wp), POINTER, DIMENSION(:) :: water_excess
    REAL(wp), POINTER, DIMENSION(:) :: w_glac


    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc) :: &
      & w_skin_max, &
      & w_skin_canopy_max, &
      & trans_tmp

    LOGICAL  :: lstart
    INTEGER  :: iblk, ics, ice
    REAL(wp) :: dtime, w_skin_max_config

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_hydrology'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    dtime = options%dtime

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    SELECT TYPE (cc => model%processes(hydro_)%p%config); TYPE IS (t_hydro_config) ; hydro__conf => cc ; END SELECT
    SELECT TYPE (cc => model%processes(sse_)%p%config); TYPE IS (t_sse_config) ; sse__conf => cc ; END SELECT

    ! Get reference to variables for current block
    !
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT

    IF (tile%is_lake) THEN
      ! compute runoff (P-E) for lakes here and exit
      ! Note: for the purpose of the water balance it is assumed here that all rain and snow fall is going into runoff
      ! @todo: compute actual water balance (snow/ice melt) in lake model and use for runoff. The corresponding energy
      ! fluxes and snow/ice budgets are already computed in the lake model.
      rain => a2l__mem%rain%ptr(ics:ice,iblk)               ! in
      snow => a2l__mem%snow%ptr(ics:ice,iblk)               ! in
      evapotrans => hydro__mem%evapotrans%ptr(ics:ice,iblk)         ! in
      evapopot => hydro__mem%evapopot%ptr(ics:ice,iblk)           ! in
      runoff => hydro__mem%runoff%ptr(ics:ice,iblk)             ! out
      drainage => hydro__mem%drainage%ptr(ics:ice,iblk)           ! out
      w_snow => hydro__mem%w_snow%ptr(ics:ice,iblk)             ! out
      evapo_snow => hydro__mem%evapo_snow%ptr(ics:ice,iblk)         ! out

      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_0 tile_name=tile%name
      runoff  (:) = rain(:) + snow(:) + evapotrans(:)
      drainage(:) = 0._wp
      w_snow(:) = 0._wp
      evapo_snow(:) = evapopot(:)
      !$claw end expand

      RETURN
    END IF

    !$claw acc data pcreate(w_skin_max, w_skin_canopy_max, trans_tmp)

    SELECT TYPE (mm => tile%mem(sse_)%p); TYPE IS (t_sse_memory); sse__mem => mm ; END SELECT
    IF (tile%contains_vegetation) THEN
      SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT
    END IF

    t_air => a2l__mem%t_air%ptr(ics:ice,iblk)          ! in
    wind_10m => a2l__mem%wind_10m%ptr(ics:ice,iblk)       ! in
    rain => a2l__mem%rain%ptr(ics:ice,iblk)           ! in
    snow => a2l__mem%snow%ptr(ics:ice,iblk)           ! in

    hcap_grnd => sse__mem%hcap_grnd%ptr(ics:ice,iblk)      ! out

    w_snow => hydro__mem%w_snow%ptr(ics:ice,iblk)         ! inout
    w_snow_soil => hydro__mem%w_snow_soil%ptr(ics:ice,iblk)    ! out
    snow_soil_dens => hydro__mem%snow_soil_dens%ptr(ics:ice,iblk) ! inout
    fract_snow => hydro__mem%fract_snow%ptr(ics:ice,iblk)     ! in
    evapotrans => hydro__mem%evapotrans%ptr(ics:ice,iblk)     ! in
    evapopot => hydro__mem%evapopot%ptr(ics:ice,iblk)       ! in
    evapo_snow => hydro__mem%evapo_snow%ptr(ics:ice,iblk)     ! out
    q_snocpymlt => hydro__mem%q_snocpymlt%ptr(ics:ice,iblk)    ! out
    snowmelt => hydro__mem%snowmelt%ptr(ics:ice,iblk)       ! out
    runoff => hydro__mem%runoff%ptr(ics:ice,iblk)         ! out
    drainage => hydro__mem%drainage%ptr(ics:ice,iblk)       ! out

    t_unfilt => seb__mem%t_unfilt%ptr(ics:ice,iblk)       ! inout

    IF (tile%is_glacier) THEN
      w_glac => hydro__mem%w_glac%ptr(ics:ice,iblk)        ! inout
      runoff_glac => hydro__mem%runoff_glac%ptr(ics:ice,iblk)   ! inout
    ELSE
      water_excess => hydro__mem%water_excess%ptr(ics:ice,iblk)    ! out
      fract_water => hydro__mem%fract_water%ptr(ics:ice,iblk)     ! in
      w_skin => hydro__mem%w_skin%ptr(ics:ice,iblk)          ! inout
      snow_accum => hydro__mem%snow_accum%ptr(ics:ice,iblk)      ! out
      evapotrans_soil => hydro__mem%evapotrans_soil%ptr(ics:ice,iblk) ! out
      evapo_skin => hydro__mem%evapo_skin%ptr(ics:ice,iblk)      ! out
      evapo_deficit => hydro__mem%evapo_deficit%ptr(ics:ice,iblk)   ! out
      IF (tile%is_vegetation .OR. tile%is_land) THEN
        w_snow_can => hydro__mem%w_snow_can%ptr(ics:ice,iblk)    ! inout
        trans => hydro__mem%trans%ptr(ics:ice,iblk)         ! in
        lai => pheno__mem%lai%ptr(ics:ice,iblk)           ! in
        fract_fpc_max => pheno__mem%fract_fpc_max%ptr(ics:ice,iblk) ! in
      END IF
    END IF

    IF (is_time_experiment_start()) THEN          ! Start of experiment
      lstart = .TRUE.
      w_snow(:) = 0._wp
    ELSE
      lstart = .FALSE.
    END IF

    ! Maximum capacity of skin reservoir (soil + canopy)
    w_skin_max_config = hydro__conf%w_skin_max
    IF (tile%is_vegetation .OR. tile%is_land) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_1 tile_name=tile%name
      w_skin_canopy_max(:) = w_skin_max_config * lai(:) * fract_fpc_max(:)
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_2 tile_name=tile%name
      w_skin_canopy_max(:) = 0._wp
    END IF
    IF (tile%contains_soil) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_3 tile_name=tile%name
      w_skin_max(:) =  w_skin_max_config + w_skin_canopy_max(:)
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_4 tile_name=tile%name
      w_skin_max(:) = 0._wp
    END IF

    ! w_snow_soil_old(:) = w_snow_soil(:) ! Remember old snow on ground for use in snow model

    IF (tile%is_glacier) THEN
      !$claw sca forward update savepoint=ser_jsb_hydro_calc_surface_hydrology_glacier_1 tile_name=tile%name
      CALL calc_surface_hydrology_glacier( &
        & lstart,                          & !< in
        & dtime,                           & !< in
        & t_unfilt(:),                     & !< in
        & fract_snow(:),                   & !< in
        & hcap_grnd(:),                    & !< in
        & evapotrans(:),                   & !< in
        & evapopot(:),                     & !< in
        & rain(:),                         & !< in
        & snow(:),                         & !< in
        & w_glac(:),                       & !< inout
        & q_snocpymlt(:),                  & !< out
        & snowmelt(:),                     & !< out
        & runoff_glac(:),                  & !< out
        & runoff(:)                        & !< out, P-E ... used as runoff for HD model
        & )
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_5 tile_name=tile%name
      drainage(:) = 0._wp
      w_snow_soil(:) = 0._wp
      evapo_snow(:) = evapopot(:)
      !$claw end expand
      IF (model%config%l_compat401) THEN
        !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_6 tile_name=tile%name
        snow_soil_dens(:) = 330._wp
      ELSE
        !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_7 tile_name=tile%name
        snow_soil_dens(:) = dens_snow
      END IF
    ELSE
      IF (tile%is_vegetation) THEN
        !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_8 tile_name=tile%name
        trans_tmp(:) = trans(:)
      ELSE
        !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_9 tile_name=tile%name
        trans_tmp(:) = 0._wp
      END IF
      !$claw sca forward update savepoint=ser_jsb_hydro_calc_surface_hydrology_land_1 tile_name=tile%name
      CALL calc_surface_hydrology_land(   &
        & lstart,                         & !< in
        & dtime,                          & !< in
        & sse__conf%l_dynsnow, & !< in
        & t_unfilt(:),                    & !< in
        & wind_10m(:),                    & !< in
        & t_air(:),                       & !< in
        & w_skin_canopy_max(:),           & !< in
        & w_skin_max(:),                  & !< in
        & fract_snow(:),                  & !< in
        & fract_water(:),                 & !< in
        & hcap_grnd(:),                   & !< in
        & evapotrans(:),                  & !< in
        & evapopot(:),                    & !< in
        & trans_tmp(:),                   & !< in
        & rain(:),                        & !< in
        & snow(:),                        & !< in
        & w_skin(:),                      & !< inout
        & w_snow_soil(:),                 & !< inout
        & w_snow_can(:),                  & !< inout
        & snow_soil_dens(:),              & !< inout
        & q_snocpymlt(:),                 & !< out
        & snow_accum(:),                  & !< out
        & snowmelt(:),                    & !< out
        & evapotrans_soil(:),             & !< out
        & evapo_skin(:),                  & !< out
        & evapo_snow(:),                  & !< out
        & water_excess(:),                & !< out
        & evapo_deficit(:)                & !< out
        & )

      IF (model%config%l_compat401) THEN
        !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_10 tile_name=tile%name
        snow_soil_dens(:) = 330._wp
      END IF

    END IF

    IF (tile%is_glacier) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_11 tile_name=tile%name
      w_snow(:) = 0._wp
    ELSE IF (tile%is_bare) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_12 tile_name=tile%name
      w_snow(:) = w_snow_soil(:)
    ELSE IF (tile%is_vegetation .OR. tile%is_land) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_13 tile_name=tile%name
      w_snow(:) = w_snow_soil(:) + w_snow_can(:)
    END IF

    !$claw acc end data

  END SUBROUTINE update_surface_hydrology
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_surface_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options


    TYPE(t_hydro_memory), POINTER :: hydro__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_surface_hydrology'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_skin,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_snow,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_snow_soil,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_snow_can,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%snowmelt,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%q_snocpymlt,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_glac,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%snow_accum,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%snow_soil_dens,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%evapotrans_soil,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%evapo_skin,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%evapo_snow,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%water_excess,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%evapo_deficit,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%runoff,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%drainage,ics,ice,iblk,routine)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_hydrology
  !
  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "soil_properties"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_soil_properties(tile, options)

    USE mo_hydro_constants, ONLY: vol_porosity_org, hyd_cond_sat_org, bclapp_org, matrix_pot_org, pore_size_index_org, &
      & vol_field_cap_org, vol_p_wilt_org

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_hydro_config), POINTER :: hydro__conf
    TYPE(t_hydro_memory), POINTER :: hydro__mem

    ! Pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: &
      & hyd_cond_sat,         &
      & vol_porosity,         &
      & bclapp,               &
      & matrix_pot,           &
      & pore_size_index,      &
      & vol_field_cap,        &
      & vol_p_wilt
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & fract_org_sl,         &
      & hyd_cond_sat_sl,      &
      & vol_porosity_sl,      &
      & bclapp_sl,            &
      & matrix_pot_sl,        &
      & pore_size_index_sl,   &
      & vol_field_cap_sl,     &
      & vol_p_wilt_sl

    ! Locally allocated vectors
    !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_w

    INTEGER  :: iblk, ics, ice, i
    INTEGER  :: nsoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_properties'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (tile%is_lake .OR. tile%is_glacier) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    SELECT TYPE (cc => model%processes(hydro_)%p%config); TYPE IS (t_hydro_config) ; hydro__conf => cc ; END SELECT

    ! Get reference to variables for current block
    !
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    hyd_cond_sat => hydro__mem%hyd_cond_sat%ptr(ics:ice,iblk)
    vol_porosity => hydro__mem%vol_porosity%ptr(ics:ice,iblk)
    bclapp => hydro__mem%bclapp%ptr(ics:ice,iblk)
    matrix_pot => hydro__mem%matrix_pot%ptr(ics:ice,iblk)
    pore_size_index => hydro__mem%pore_size_index%ptr(ics:ice,iblk)
    vol_field_cap => hydro__mem%vol_field_cap%ptr(ics:ice,iblk)
    vol_p_wilt => hydro__mem%vol_p_wilt%ptr(ics:ice,iblk)
    IF (hydro__conf%l_organic) fract_org_sl => hydro__mem%fract_org_sl%ptr(ics:ice,:,iblk)
    hyd_cond_sat_sl => hydro__mem%hyd_cond_sat_sl%ptr(ics:ice,:,iblk)
    vol_porosity_sl => hydro__mem%vol_porosity_sl%ptr(ics:ice,:,iblk)
    bclapp_sl => hydro__mem%bclapp_sl%ptr(ics:ice,:,iblk)
    matrix_pot_sl => hydro__mem%matrix_pot_sl%ptr(ics:ice,:,iblk)
    pore_size_index_sl => hydro__mem%pore_size_index_sl%ptr(ics:ice,:,iblk)
    vol_field_cap_sl => hydro__mem%vol_field_cap_sl%ptr(ics:ice,:,iblk)
    vol_p_wilt_sl => hydro__mem%vol_p_wilt_sl%ptr(ics:ice,:,iblk)

    DO i=1,nsoil
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_soil_properties_expand_0 tile_name=tile%name soil=i
      hyd_cond_sat_sl   (:,i) = hyd_cond_sat(:)
      vol_porosity_sl   (:,i) = vol_porosity(:)
      bclapp_sl         (:,i) = bclapp(:)
      matrix_pot_sl     (:,i) = matrix_pot(:)
      pore_size_index_sl(:,i) = pore_size_index(:)
      vol_field_cap_sl  (:,i) = vol_field_cap(:)
      vol_p_wilt_sl     (:,i) = vol_p_wilt(:)
      !$claw end expand

    ! hyd_cond_sat_sl   (:,:) = SPREAD(hyd_cond_sat   (:), DIM=2, ncopies=nsoil)
    ! vol_porosity_sl   (:,:) = SPREAD(vol_porosity   (:), DIM=2, ncopies=nsoil)
    ! bclapp_sl         (:,:) = SPREAD(bclapp         (:), DIM=2, ncopies=nsoil)
    ! matrix_pot_sl     (:,:) = SPREAD(matrix_pot     (:), DIM=2, ncopies=nsoil)
    ! pore_size_index_sl(:,:) = SPREAD(pore_size_index(:), DIM=2, ncopies=nsoil)
    ! vol_field_cap_sl  (:,:) = SPREAD(vol_field_cap  (:), DIM=2, ncopies=nsoil)
    ! vol_p_wilt_sl     (:,:) = SPREAD(vol_p_wilt     (:), DIM=2, ncopies=nsoil)

    END DO
    IF (hydro__conf%l_organic) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_soil_properties_expand_1 tile_name=tile%name soil=i
      hyd_cond_sat_sl   (:,:) = (1 - fract_org_sl(:,:)) * hyd_cond_sat_sl   (:,:) + fract_org_sl(:,:) * hyd_cond_sat_org
      vol_porosity_sl   (:,:) = (1 - fract_org_sl(:,:)) * vol_porosity_sl   (:,:) + fract_org_sl(:,:) * vol_porosity_org
      bclapp_sl         (:,:) = (1 - fract_org_sl(:,:)) * bclapp_sl         (:,:) + fract_org_sl(:,:) * bclapp_org
      matrix_pot_sl     (:,:) = (1 - fract_org_sl(:,:)) * matrix_pot_sl     (:,:) + fract_org_sl(:,:) * matrix_pot_org
      pore_size_index_sl(:,:) = (1 - fract_org_sl(:,:)) * pore_size_index_sl(:,:) + fract_org_sl(:,:) * pore_size_index_org
      vol_field_cap_sl  (:,:) = (1 - fract_org_sl(:,:)) * vol_field_cap_sl  (:,:) + fract_org_sl(:,:) * vol_field_cap_org
      vol_p_wilt_sl     (:,:) = (1 - fract_org_sl(:,:)) * vol_p_wilt_sl     (:,:) + fract_org_sl(:,:) * vol_p_wilt_org
      !$claw end expand
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_properties
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "soil_properties"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_soil_properties(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_hydro_memory), POINTER :: hydro__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_properties'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, hydro__mem%hyd_cond_sat_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%vol_porosity_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%bclapp_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%matrix_pot_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%pore_size_index_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%vol_field_cap_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%vol_p_wilt_sl,ics,ice,iblk,routine)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_properties
  !
  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "soil_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_soil_hydrology(tile, options)

    USE mo_hydro_process,          ONLY: calc_soil_hydrology

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_hydro_config), POINTER :: hydro__conf
    TYPE(t_sse_memory), POINTER :: sse__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem

    ! Pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:,:) :: fract_org_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: t_soil

    REAL(wp), POINTER, DIMENSION(:) :: w_soil_rel
    REAL(wp), POINTER, DIMENSION(:) :: w_soil_column
    REAL(wp), POINTER, DIMENSION(:) :: max_moist


    ! Locally allocated vectors

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_w

    LOGICAL  :: ltpe_open
    LOGICAL  :: ltpe_closed

    INTEGER  :: iblk, ics, ice, nc, i
    REAL(wp) :: dtime
    INTEGER  :: nsoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_hydrology'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (tile%is_lake) RETURN

    ! Runoff and drainage for glaciers already computed in update_surface_hydrology
    ! Glaciers have no water in soil
    IF (tile%is_glacier) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    grid => Get_grid(model%grid_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    SELECT TYPE (cc => model%processes(hydro_)%p%config); TYPE IS (t_hydro_config) ; hydro__conf => cc ; END SELECT

    ltpe_open = .FALSE.
    ltpe_closed = .FALSE.
    IF (model%config%tpe_scheme == 'open')   ltpe_open = .TRUE.
    IF (model%config%tpe_scheme == 'closed') ltpe_closed = .TRUE.

    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(sse_)%p); TYPE IS (t_sse_memory); sse__mem => mm ; END SELECT

    w_soil_rel => hydro__mem%w_soil_rel%ptr(ics:ice,iblk)
    w_soil_column => hydro__mem%w_soil_column%ptr(ics:ice,iblk)
    max_moist => hydro__mem%max_moist%ptr(ics:ice,iblk)

    t_soil => sse__mem%t_soil%ptr(ics:ice,:,iblk)     ! in

    IF (hydro__conf%l_organic) THEN
      fract_org_sl => hydro__mem%fract_org_sl%ptr(ics:ice,:,iblk)
    ELSE
      ALLOCATE(fract_org_sl(nc, nsoil))
      !$claw acc enter data create(fract_org_sl)
      !$claw expand parallel update
      fract_org_sl(:,:) = 0._wp
    END IF

    !$claw acc update device(soil_w%dz) ! TODO move if possible
    !$claw sca forward update savepoint=ser_jsb_hydro_calc_soil_hydrology_1 tile_name=tile%name
    CALL calc_soil_hydrology( &
      ! in
      & nc, &
      & nsoil, &
      & soil_w%dz(:), &
      & grid%nlat_g, & !< (Effective) number of latitudes for steepness parameter in infiltration computation<
      & dtime, &
      & ltpe_closed, ltpe_open,                             &
      & hydro__mem%soil_depth_sl%ptr(ics:ice,:,iblk), &
      & hydro__mem%root_depth_sl%ptr(ics:ice,:,iblk), &
      &                                fract_org_sl  (:,:), &
      & hydro__mem%hyd_cond_sat_sl%ptr(ics:ice,:,iblk), &
      & hydro__mem%matrix_pot_sl%ptr(ics:ice,:,iblk), &
      & hydro__mem%bclapp_sl%ptr(ics:ice,:,iblk), &
      & hydro__mem%pore_size_index_sl%ptr(ics:ice,:,iblk), &
      & hydro__mem%vol_porosity_sl%ptr(ics:ice,:,iblk), &
      & hydro__mem%vol_field_cap_sl%ptr(ics:ice,:,iblk), &
      & hydro__mem%vol_p_wilt_sl%ptr(ics:ice,:,iblk), &
      & hydro__mem%max_moist%ptr(ics:ice,iblk), &
      & hydro__mem%oro_stddev%ptr(ics:ice,iblk), &
      &                                t_soil        (:,1), &
      & hydro__mem%water_excess%ptr(ics:ice,iblk), &
      & hydro__mem%evapotrans_soil%ptr(ics:ice,iblk), &
      & hydro__mem%trans%ptr(ics:ice,iblk), &
      ! inout
      & hydro__mem%w_soil_column%ptr(ics:ice,iblk), &
      & hydro__mem%w_soil_sl%ptr(ics:ice,:,iblk), &
      & hydro__mem%w_ice_sl%ptr(ics:ice,:,iblk), &
      & hydro__mem%w_soil_overflow%ptr(ics:ice,iblk), &
      ! out
      & hydro__mem%w_soil_sat_sl%ptr(ics:ice,:,iblk), &
      & hydro__mem%w_soil_fc_sl%ptr(ics:ice,:,iblk), &
      & hydro__mem%w_soil_pwp_sl%ptr(ics:ice,:,iblk), &
      & hydro__mem%runoff%ptr(ics:ice,iblk), &
      & hydro__mem%drainage%ptr(ics:ice,iblk) &
      & )

    !$claw acc data present(max_moist, w_soil_rel, w_soil_column)
    !$claw acc parallel loop gang vector
    DO i = 1, nc
      IF (max_moist(i) > 0._wp) THEN
        w_soil_rel(i) = MIN(w_soil_column(i) / max_moist(i), 1._wp)
          ! @todo - in jsbach3, taking the MIN with 1 is not nessecary. Here, w_soil_column can be larger than
          ! max_moist (in some places) ... this should be investigated, maybe it has to do with the interpolated input fields
      ELSE
        ! max_moist is zero over glaciers
        w_soil_rel(i) = 0._wp
      END IF
    END DO
    !$claw acc end data

    IF (.NOT. hydro__conf%l_organic) THEN
      !$claw acc exit data delete(fract_org_sl)
      DEALLOCATE(fract_org_sl)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_hydrology
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "soil_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_soil_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_hydro_memory), POINTER :: hydro__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_hydrology'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_soil_column,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_soil_rel,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_soil_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_ice_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_soil_sat_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_soil_fc_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_soil_pwp_sl,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%runoff,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%drainage,ics,ice,iblk,routine)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_hydrology
  !
  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "canopy_cond_unstressed"
  !! Task "update_canopy_cond_unstressed" calculates canopy conductance without water stress
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_canopy_cond_unstressed(tile, options)
    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_hydro_process, ONLY: get_canopy_conductance

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_hydro_memory), POINTER :: hydro__mem
    TYPE(t_pheno_memory), POINTER :: pheno__mem
    TYPE(t_a2l_memory), POINTER :: a2l__mem

    REAL(wp), POINTER, DIMENSION(:) :: lai
    REAL(wp), POINTER, DIMENSION(:) :: swpar_srf_down
    REAL(wp), POINTER, DIMENSION(:) :: canopy_cond_unlimited

    INTEGER :: iblk, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_canopy_cond_unstressed'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (.NOT. tile%Is_process_active(HYDRO_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get reference to variables for current block
    !
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT

    swpar_srf_down => a2l__mem%swpar_srf_down%ptr(ics:ice,iblk)        ! IN
    lai => pheno__mem%lai%ptr(ics:ice,iblk)                   ! IN
    canopy_cond_unlimited => hydro__mem%canopy_cond_unlimited%ptr(ics:ice,iblk) ! OUT

    ! Compute (max) canopy (stomatal) conductance with no water stress
    !$claw expand parallel
    canopy_cond_unlimited(:) = get_canopy_conductance( lai(:), swpar_srf_down(:) )

  END SUBROUTINE update_canopy_cond_unstressed
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "canopy_cond_unstressed"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_canopy_cond_unstressed(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_hydro_memory), POINTER :: hydro__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_cond_unstressed'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, hydro__mem%canopy_cond_unlimited,ics,ice,iblk,routine)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_canopy_cond_unstressed

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! @brief Implementation of "update" for task "water_stress"
  !! Task "update_water_stress" calculates water stress
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_water_stress(tile, options)
    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_hydro_process, ONLY: get_water_stress_factor
    USE mo_hydro_util,    ONLY: get_water_in_root_zone

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_hydro_config), POINTER :: hydro__conf

    TYPE(t_hydro_memory), POINTER :: hydro__mem

    REAL(wp), POINTER, DIMENSION(:,:) :: soil_depth_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: root_depth_sl
    REAL(wp), POINTER, DIMENSION(:) :: w_soil_column       !< Water content in soil column
    REAL(wp), POINTER, DIMENSION(:,:) :: w_soil_sl           !< Water content in soil layers
    REAL(wp), POINTER, DIMENSION(:,:) :: w_soil_fc_sl        !< Water content of soil layers at field capacity
    REAL(wp), POINTER, DIMENSION(:,:) :: w_soil_pwp_sl       !< Water content of soil layers at permanent wilting point
    REAL(wp), POINTER, DIMENSION(:) :: w_soil_root         !< Water content in root zone
    REAL(wp), POINTER, DIMENSION(:) :: w_soil_root_fc      !< Water content in root zone at field capacity
    REAL(wp), POINTER, DIMENSION(:) :: w_soil_root_pwp     !< Water content in root zone at permanent wilting point
    REAL(wp), POINTER, DIMENSION(:) :: max_moist
    REAL(wp), POINTER, DIMENSION(:) :: water_stress

    INTEGER :: iblk, ics, ice
    REAL(wp) :: config_w_soil_crit_fract, config_w_soil_wilt_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_water_stress'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (.NOT. tile%Is_process_active(HYDRO_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    model => Get_model(tile%owner_model_id)
    SELECT TYPE (cc => model%processes(hydro_)%p%config); TYPE IS (t_hydro_config) ; hydro__conf => cc ; END SELECT

    ! Get reference to variables for current block
    !
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT
    soil_depth_sl => hydro__mem%soil_depth_sl%ptr(ics:ice,:,iblk)       ! in
    w_soil_sl => hydro__mem%w_soil_sl%ptr(ics:ice,:,iblk)           ! in
    w_soil_fc_sl => hydro__mem%w_soil_fc_sl%ptr(ics:ice,:,iblk)        ! in
    w_soil_pwp_sl => hydro__mem%w_soil_pwp_sl%ptr(ics:ice,:,iblk)       ! in
    max_moist => hydro__mem%max_moist%ptr(ics:ice,iblk)           ! in
    w_soil_column => hydro__mem%w_soil_column%ptr(ics:ice,iblk)       ! in
    root_depth_sl => hydro__mem%root_depth_sl%ptr(ics:ice,:,iblk)       ! in
    w_soil_root => hydro__mem%w_soil_root%ptr(ics:ice,iblk)         ! out
    w_soil_root_fc => hydro__mem%w_soil_root_fc%ptr(ics:ice,iblk)      ! out
    w_soil_root_pwp => hydro__mem%w_soil_root_pwp%ptr(ics:ice,iblk)     ! out
    water_stress => hydro__mem%water_stress%ptr(ics:ice,iblk)        ! out

    ! Get water stress

    ! Calculate water content and field capacity in the root zone
    ! @todo - do we need these?

    config_w_soil_crit_fract = hydro__conf%w_soil_crit_fract
    config_w_soil_wilt_fract = hydro__conf%w_soil_wilt_fract

    !TODO CLAW: This part is computed on CPU because temp arrays are issued by compiler. 
    !$claw acc data present(w_soil_root, w_soil_sl, soil_depth_sl, root_depth_sl, &
    !$claw acc w_soil_root_fc, w_soil_fc_sl, w_soil_root_pwp, w_soil_pwp_sl, w_soil_pwp_sl)
    !$claw acc update host(w_soil_sl, soil_depth_sl, root_depth_sl, w_soil_fc_sl, w_soil_pwp_sl)
    w_soil_root   (:)  = get_water_in_root_zone(w_soil_sl     (:,:), soil_depth_sl(:,:), root_depth_sl(:,:))
    w_soil_root_fc(:)  = get_water_in_root_zone(w_soil_fc_sl  (:,:), soil_depth_sl(:,:), root_depth_sl(:,:))
    w_soil_root_pwp(:) = get_water_in_root_zone(w_soil_pwp_sl (:,:), soil_depth_sl(:,:), root_depth_sl(:,:))
    !$claw acc update device(w_soil_root,w_soil_root_fc,w_soil_root_pwp) ! COMPUTED ON CPU
    !$claw acc end data

    ! Calculate water stress factor
    ! @todo: maybe use the distributed permanent wilting point/field cap from above and change get_water_stress_factor function
    !$claw expand parallel update savepoint=ser_jsb_hydro_update_water_stress_expand_0 tile_name=tile%name
    water_stress(:) = &
      & get_water_stress_factor(w_soil_column(:), max_moist(:), &
          &                     config_w_soil_crit_fract, config_w_soil_wilt_fract)
    !$claw end expand

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_water_stress

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "water_stress"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_water_stress(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_hydro_memory), POINTER :: hydro__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_water_stress'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_soil_root,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_soil_root_fc,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%w_soil_root_pwp,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%water_stress,ics,ice,iblk,routine)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_water_stress

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

    ! TYPE(t_jsb_model), POINTER :: model

    TYPE(t_seb_memory), POINTER :: seb__mem
    TYPE(t_a2l_memory), POINTER :: a2l__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem

    REAL(wp), POINTER, DIMENSION(:) :: canopy_cond_unlimited
    REAL(wp), POINTER, DIMENSION(:) :: canopy_cond_limited
    REAL(wp), POINTER, DIMENSION(:) :: water_stress
    REAL(wp), POINTER, DIMENSION(:) :: t
    REAL(wp), POINTER, DIMENSION(:) :: q_air
    REAL(wp), POINTER, DIMENSION(:) :: press_srf

    INTEGER :: iblk, ics, ice, nc
    ! LOGICAL, DIMENSION(options%nc) :: q_air_gt_qsat_tmp
    LOGICAL, ALLOCATABLE :: q_air_gt_qsat_tmp(:)


    CHARACTER(len=*), PARAMETER :: routine = modname//':update_canopy_cond_stressed'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_active(HYDRO_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    q_air => a2l__mem%q_air%ptr(ics:ice,iblk)               ! in
    press_srf => a2l__mem%press_srf%ptr(ics:ice,iblk)           ! in
    t => seb__mem%t%ptr(ics:ice,iblk)                   ! in
    water_stress => hydro__mem%water_stress%ptr(ics:ice,iblk)        ! in
    canopy_cond_unlimited => hydro__mem%canopy_cond_unlimited%ptr(ics:ice,iblk) ! in
    canopy_cond_limited => hydro__mem%canopy_cond_limited%ptr(ics:ice,iblk) ! out

    !$cla$w acc data create(q_air_gt_qsat_tmp)
    ALLOCATE(q_air_gt_qsat_tmp(nc))
    !$claw acc enter data create(q_air_gt_qsat_tmp)
    !$claw expand parallel update savepoint=ser_jsb_hydro_update_canopy_cond_stressed_expand_0 tile_name=tile%name
    q_air_gt_qsat_tmp(:) = q_air(:) > qsat(t(:),press_srf(:))
    ! Compute (actual) canopy (stomatal) conductance under water stress.
    canopy_cond_limited(:) = get_canopy_conductance(canopy_cond_unlimited(:),  & ! in, unstressed canopy conductance
                                                    water_stress(:),           & ! in, water stress factor
                                                    q_air_gt_qsat_tmp(:)       & ! in, atmosphere saturated?
                                                    )
    !$claw end expand
    !$cla$w acc end data
    !$claw acc exit data delete(q_air_gt_qsat_tmp)
    DEALLOCATE(q_air_gt_qsat_tmp)

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

    TYPE(t_hydro_memory), POINTER :: hydro__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_cond_stressed'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, hydro__mem%canopy_cond_limited,ics,ice,iblk,routine)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_canopy_cond_stressed

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! @brief Implementation of "update" for task "evaporation"
  !! Task "update_evaporation" calculates evapo(transpi)ration
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_evaporation(tile, options)

    USE mo_phy_schemes,            ONLY: q_effective, heat_transfer_coef

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_seb_config), POINTER :: seb__conf
    TYPE(t_a2l_memory), POINTER :: a2l__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem
    TYPE(t_seb_memory), POINTER :: seb__mem
    TYPE(t_turb_memory), POINTER :: turb__mem

    ! Pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: &
      & q_acoef, &
      & q_bcoef, &
      & q_acoef_wtr, &
      & q_bcoef_wtr, &
      & q_acoef_ice, &
      & q_bcoef_ice, &
      & drag_srf, &
      & drag_wtr, &
      & drag_ice, &
      & qsat_star, &
      & qsat_lwtr, &
      & qsat_lice, &
      & fact_qsat_srf, &
      & fact_qsat_trans_srf, &
      & fact_q_air, &
      & fract_lice, &
      & evapotrans_lnd, &
      & evapo_wtr, &
      & evapo_ice, &
      & evapopot,       &
      & evapotrans, &
      & trans

    ! Locally allocated vectors
    !
    ! REAL(wp) ::                  &
    !   & q_air        (options%nc), &  !< Humidity at lowest atmospheric level [kg kg-1]
    !   & q_air_eff    (options%nc), &  !< Effective humidity at lowest atmospheric level
    !   & qsat_srf_eff (options%nc), &  !< Effective surface saturation specific humidity
    !   & heat_tcoef   (options%nc)     !< Heat transfer coefficient (rho*C_h*|v|)
    REAL(wp), ALLOCATABLE ::                  &
      & q_air        (:), &  !< Humidity at lowest atmospheric level [kg kg-1]
      & q_air_eff    (:), &  !< Effective humidity at lowest atmospheric level
      & qsat_srf_eff (:), &  !< Effective surface saturation specific humidity
      & heat_tcoef   (:)     !< Heat transfer coefficient (rho*C_h*|v|)

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: steplen

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_evaporation'

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    !$cla$w acc data create(q_air, q_air_eff, qsat_srf_eff, heat_tcoef)
    ALLOCATE(q_air(nc),q_air_eff(nc),qsat_srf_eff(nc),heat_tcoef(nc))
    !$claw acc enter data create(q_air, q_air_eff, qsat_srf_eff, heat_tcoef)

    ! Get reference to variables for current block
    !
    SELECT TYPE (cc => model%processes(seb_)%p%config); TYPE IS (t_seb_config) ; seb__conf => cc ; END SELECT
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(turb_)%p); TYPE IS (t_turb_memory); turb__mem => mm ; END SELECT

    ! Pointers to variables in memory
    evapopot => hydro__mem%evapopot%ptr(ics:ice,iblk)   ! out
    evapotrans => hydro__mem%evapotrans%ptr(ics:ice,iblk) ! out

    IF (tile%contains_lake) THEN

      q_acoef_wtr => a2l__mem%q_acoef_wtr%ptr(ics:ice,iblk)   ! in
      q_bcoef_wtr => a2l__mem%q_bcoef_wtr%ptr(ics:ice,iblk)   ! in
      drag_wtr => a2l__mem%drag_wtr%ptr(ics:ice,iblk)      ! in
      qsat_lwtr => seb__mem%qsat_lwtr%ptr(ics:ice,iblk)     ! in
      evapo_wtr => hydro__mem%evapo_wtr%ptr(ics:ice,iblk)     ! out

      !$claw expand parallel update savepoint=ser_jsb_hydro_update_evaporation_expand_1 tile_name=tile%name
      q_air(:) = q_acoef_wtr(:) * qsat_lwtr(:) + q_bcoef_wtr(:) ! New moisture at lowest atmospheric level by back-substitution
      heat_tcoef(:) = heat_transfer_coef(drag_wtr(:), steplen)  ! Transfer coefficient
      evapo_wtr (:) = heat_tcoef(:) * (q_air(:) - qsat_lwtr(:)) ! Potential evaporation
      !$claw end expand

      IF (seb__conf%l_ice_on_lakes) THEN
        q_acoef_ice => a2l__mem%q_acoef_ice%ptr(ics:ice,iblk)       ! in
        q_bcoef_ice => a2l__mem%q_bcoef_ice%ptr(ics:ice,iblk)       ! in
        drag_ice => a2l__mem%drag_ice%ptr(ics:ice,iblk)          ! in
        qsat_lice => seb__mem%qsat_lice%ptr(ics:ice,iblk)         ! in
        fract_lice => seb__mem%fract_lice%ptr(ics:ice,iblk)        ! in
        evapo_ice => hydro__mem%evapo_ice%ptr(ics:ice,iblk)         ! out

        !$claw expand parallel update savepoint=ser_jsb_hydro_update_evaporation_expand_2 tile_name=tile%name
        q_air(:) = q_acoef_ice(:) * qsat_lice(:) + q_bcoef_ice(:) ! New moisture at lowest atmospheric level by back-substitution
        heat_tcoef(:) = heat_transfer_coef(drag_ice(:), steplen)  ! Transfer coefficient
        evapo_ice (:) = heat_tcoef(:) * (q_air(:) - qsat_lice(:)) ! Potential evaporation
        evapopot  (:) = (1._wp - fract_lice(:)) * evapo_wtr   (:) + fract_lice(:) * evapo_ice   (:)
        !$claw end expand
      ELSE
        !$claw expand parallel update savepoint=ser_jsb_hydro_update_evaporation_expand_3 tile_name=tile%name
        evapopot  (:) = evapo_wtr(:)
      END IF
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_evaporation_expand_4 tile_name=tile%name
      evapotrans(:) = evapopot(:)

    ELSE IF (tile%contains_soil .OR. tile%contains_glacier) THEN

      q_acoef => a2l__mem%q_acoef%ptr(ics:ice,iblk)             ! in
      q_bcoef => a2l__mem%q_bcoef%ptr(ics:ice,iblk)             ! in
      drag_srf => a2l__mem%drag_srf%ptr(ics:ice,iblk)            ! in
      qsat_star => seb__mem%qsat_star%ptr(ics:ice,iblk)           ! in
      fact_q_air => turb__mem%fact_q_air%ptr(ics:ice,iblk)          ! in
      fact_qsat_srf => turb__mem%fact_qsat_srf%ptr(ics:ice,iblk)       ! in
      fact_qsat_trans_srf => turb__mem%fact_qsat_trans_srf%ptr(ics:ice,iblk) ! in

      evapotrans_lnd => hydro__mem%evapotrans_lnd%ptr(ics:ice,iblk)           ! out
      IF (tile%contains_vegetation) THEN
        trans => hydro__mem%trans%ptr(ics:ice,iblk)               ! out
      END IF

      !$claw expand parallel update savepoint=ser_jsb_hydro_update_evaporation_expand_5 tile_name=tile%name
      q_air(:) = q_acoef(:) * qsat_star(:) + q_bcoef(:) ! New moisture at lowest atmospheric level by back-substitution
      heat_tcoef(:) = heat_transfer_coef(drag_srf(:), steplen)  ! Transfer coefficient

      ! Compute effective air moisture and surface saturation humidity
      q_air_eff(:) = q_effective( SPREAD(0._wp, NCOPIES=nc, DIM=1), q_air(:), &
                              & SPREAD(1._wp, NCOPIES=nc, DIM=1), SPREAD(0._wp, NCOPIES=nc, DIM=1))
      qsat_srf_eff(:) = q_effective(qsat_star(:), q_air(:), fact_qsat_srf(:), fact_q_air(:))

      evapotrans_lnd(:)          = heat_tcoef(:) * (q_air_eff(:) - qsat_srf_eff(:)) ! Evapotranspiration
      evapotrans(:) = evapotrans_lnd(:)
      evapopot(:)           = heat_tcoef(:) * (q_air    (:) - qsat_star   (:)) ! Potential evaporation
      !$claw end expand
      IF (tile%contains_vegetation) THEN
        !$claw expand parallel update savepoint=ser_jsb_hydro_update_evaporation_expand_6 tile_name=tile%name
        trans(:)              = fact_qsat_trans_srf(:) * evapopot(:)             ! Transpiration
      END IF

    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

    !$cla$w acc end data
    !$claw acc exit data delete(q_air, q_air_eff, qsat_srf_eff, heat_tcoef)
    DEALLOCATE(q_air,q_air_eff,qsat_srf_eff,heat_tcoef)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_evaporation
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "evaporation"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_evaporation(tile, options)

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_seb_config), POINTER :: seb__conf
    TYPE(t_hydro_memory), POINTER :: hydro__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_evaporation'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    SELECT TYPE (cc => model%processes(seb_)%p%config); TYPE IS (t_seb_config) ; seb__conf => cc ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, hydro__mem%evapotrans,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%evapopot,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%trans,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%evapotrans_lnd,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%evapo_wtr,ics,ice,iblk,routine)
    IF (seb__conf%l_ice_on_lakes) THEN
      CALL weighted_by_fract%Aggregate(tile, hydro__mem%evapo_ice,ics,ice,iblk,routine)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_evaporation
  !
  !>
  !! @brief Implementation of "update" for task "snow_and_ice_hydrology"
  !! Task "update_snow_hydrology" calculates the new snow depth and fraction
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_snow_and_ice_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_snow_and_ice_hydrology'

    IF (options%nc > 0) CONTINUE ! avoid compiler warnings about dummy arguments not being used

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (tile%contains_lake) THEN
      !! Currently done in SEB process
    ELSE IF (tile%contains_land) THEN
!!$      CALL update_snow_and_ice_hydrology_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_snow_and_ice_hydrology
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "snow_and_ice_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_snow_and_ice_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_snow_and_ice_hydrology'

    INTEGER :: iblk !, ics, ice

    iblk = options%iblk
    !ics  = options%ics
    !ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_snow_and_ice_hydrology
  !
  ! ================================================================================================================================
  !
  !! @brief Implementation of "update" for task "snow_and_skin_fraction"
  !! Task "update_snow_and_skin_fraction" calculates the new snow and wet skin fractions
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_snow_and_skin_fraction( tile, options)

    ! Use declarations
    USE mo_hydro_process,          ONLY: calc_wskin_fractions_lice, calc_wskin_fractions_veg, calc_wskin_fractions_bare
    ! sollte jetzt in dieses File hier kommen: USE mo_hydro_process,  ONLY: calc_hydro_snow_and_skin_fraction

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    !----------------------------
    ! Local variables
    TYPE(t_jsb_model), POINTER :: model
    REAL(wp), DIMENSION(options%nc) :: w_skin_canopy_max, w_skin_max

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_snow_and_skin_fraction'

    INTEGER  :: iblk, ics, ice
    REAL(wp) :: dtime, steplen, config_w_skin_max
    !$claw verbatim INTEGER, SAVE :: task_counter = 1

    ! Declare pointers for process configuration and memory
    TYPE(t_seb_config), POINTER :: seb__conf
    TYPE(t_a2l_memory), POINTER :: a2l__mem
    TYPE(t_hydro_config), POINTER :: hydro__conf
    TYPE(t_seb_memory), POINTER :: seb__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem
    TYPE(t_pheno_memory), POINTER :: pheno__mem

    ! Declare pointers for variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: lai
    REAL(wp), POINTER, DIMENSION(:) :: fract_snow
    REAL(wp), POINTER, DIMENSION(:) :: fract_water
    REAL(wp), POINTER, DIMENSION(:) :: fract_fpc_max
    REAL(wp), POINTER, DIMENSION(:) :: fract_snow_can
    REAL(wp), POINTER, DIMENSION(:) :: fract_snow_soil
    REAL(wp), POINTER, DIMENSION(:) :: fract_snow_lice
    REAL(wp), POINTER, DIMENSION(:) :: w_snow_can
    REAL(wp), POINTER, DIMENSION(:) :: w_snow_soil
    REAL(wp), POINTER, DIMENSION(:) :: w_snow_lice
    REAL(wp), POINTER, DIMENSION(:) :: w_skin
    REAL(wp), POINTER, DIMENSION(:) :: oro_stddev
    REAL(wp), POINTER, DIMENSION(:) :: drag_srf
    REAL(wp), POINTER, DIMENSION(:) :: t
    REAL(wp), POINTER, DIMENSION(:) :: press_srf
    REAL(wp), POINTER, DIMENSION(:) :: q_air

    ! ---------------------------
    ! Go

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    dtime   = options%dtime
    steplen = options%steplen

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Set pointers to process configs and memory
    SELECT TYPE (cc => model%processes(seb_)%p%config); TYPE IS (t_seb_config) ; seb__conf => cc ; END SELECT
    SELECT TYPE (cc => model%processes(hydro_)%p%config); TYPE IS (t_hydro_config) ; hydro__conf => cc ; END SELECT

    ! Use simple scalar for GPU
    config_w_skin_max = hydro__conf%w_skin_max

    !$claw acc enter data pcreate(w_skin_max, w_skin_canopy_max)
    !$claw expand parallel update
    w_skin_max(:) = config_w_skin_max

    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    ! First handle LAKE_TYPE lct
    IF (tile%is_lake) THEN
      IF (seb__conf%l_ice_on_lakes) THEN
        w_snow_lice => hydro__mem%w_snow_lice%ptr(ics:ice,iblk)      ! in
        fract_snow_lice => hydro__mem%fract_snow_lice%ptr(ics:ice,iblk)  ! out
        !$claw sca forward update savepoint=ser_jsb_hydro_calc_wskin_fractions_lice_1 tile_name=tile%name dtime=dtime steplen=steplen task_counter=task_counter
        CALL calc_wskin_fractions_lice( &
          & w_snow_lice(:),             & ! in
          & fract_snow_lice(:)          & ! out
          & )
      END IF
      !$claw acc exit data delete(w_skin_max, w_skin_canopy_max)  
      RETURN
    END IF

    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT
    IF (tile%is_vegetation .OR. tile%is_land) THEN
      SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT
    END IF

    ! Set pointers
    drag_srf => a2l__mem%drag_srf%ptr(ics:ice,iblk)         ! in
    q_air => a2l__mem%q_air%ptr(ics:ice,iblk)            ! in
    press_srf => a2l__mem%press_srf%ptr(ics:ice,iblk)        ! in

    t => seb__mem%t%ptr(ics:ice,iblk)                ! in

    IF (tile%is_bare .OR. tile%is_vegetation .OR. tile%is_land) THEN
      oro_stddev => hydro__mem%oro_stddev%ptr(ics:ice,iblk)       ! in
      fract_water => hydro__mem%fract_water%ptr(ics:ice,iblk)  ! out
    END IF
    IF (tile%is_vegetation .OR. tile%is_land) THEN
      lai => pheno__mem%lai%ptr(ics:ice,iblk)              ! in
      fract_fpc_max => pheno__mem%fract_fpc_max%ptr(ics:ice,iblk)    ! in
      fract_snow_can => hydro__mem%fract_snow_can%ptr(ics:ice,iblk)   ! out
    END IF

    IF (tile%is_vegetation .OR. tile%is_land .OR. tile%is_bare) THEN
      w_skin => hydro__mem%w_skin%ptr(ics:ice,iblk)           ! in
    END IF

    fract_snow => hydro__mem%fract_snow%ptr(ics:ice,iblk)       ! out
    w_snow_soil => hydro__mem%w_snow_soil%ptr(ics:ice,iblk)      ! in
    fract_snow_soil => hydro__mem%fract_snow_soil%ptr(ics:ice,iblk)  ! out

    ! Maximum capacity of skin reservoir (soil + canopy)
    IF (tile%is_vegetation .OR. tile%is_land) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_snow_and_skin_fraction_expand_0 tile_name=tile%name dtime=dtime steplen=steplen task_counter=task_counter
      w_skin_canopy_max(:) = config_w_skin_max * lai(:) * fract_fpc_max(:)
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_snow_and_skin_fraction_expand_1 tile_name=tile%name dtime=dtime steplen=steplen task_counter=task_counter
      w_skin_canopy_max(:) = 0._wp
    END IF
    IF (tile%contains_soil) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_snow_and_skin_fraction_expand_2 tile_name=tile%name dtime=dtime steplen=steplen task_counter=task_counter
      w_skin_max(:) = config_w_skin_max + w_skin_canopy_max(:)
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_snow_and_skin_fraction_expand_3 tile_name=tile%name dtime=dtime steplen=steplen task_counter=task_counter
      w_skin_max(:) = 0._wp
    END IF

    IF (tile%is_glacier) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_snow_and_skin_fraction_expand_4 tile_name=tile%name dtime=dtime steplen=steplen task_counter=task_counter
      fract_snow(:)      = 1._wp
      fract_snow_soil(:) = 1._wp
      !$claw end expand
    END IF

    IF (tile%is_vegetation .OR. tile%is_land) THEN
      w_snow_can => hydro__mem%w_snow_can%ptr(ics:ice,iblk)       ! in
      !$claw sca forward update savepoint=ser_jsb_hydro_calc_wskin_fractions_veg_1 tile_name=tile%name dtime=dtime steplen=steplen task_counter=task_counter
      CALL calc_wskin_fractions_veg( &
        & dtime,                     & ! in
        & steplen,                   & ! in
        & w_skin_max(:),             & ! in
        & oro_stddev(:),             & ! in
        & t(:),                      & ! in, from the previous time step as long as this is called before the asselin filter
        & press_srf(:),              & ! in
        & drag_srf(:),               & ! in
        & q_air(:),                  & ! in
        & w_skin(:),                 & ! in
        & w_snow_soil(:),            & ! in
        & w_snow_can(:),             & ! in
        & fract_snow_can(:),         & ! out
        & fract_water(:),            & ! out
        & fract_snow_soil(:)         & ! out
        & )
      !$claw expand parallel update
      fract_snow(:) = fract_snow_soil(:)
   END IF


    IF (tile%is_bare) THEN
      !$claw sca forward update savepoint=ser_jsb_hydro_calc_wskin_fractions_bare_1 tile_name=tile%name dtime=dtime steplen=steplen  task_counter=task_counter
      CALL calc_wskin_fractions_bare( &
        & dtime,                      & ! in
        & steplen,                    & ! in
        & w_skin_max,                 & ! in
        & oro_stddev(:),              & ! in
        & t(:),                       & ! in
        & press_srf(:),               & ! in
        & drag_srf(:),                & ! in
        & q_air(:),                   & ! in
        & w_skin(:),                  & ! in
        & w_snow_soil(:),             & ! in
        & fract_water(:),             & ! out
        & fract_snow_soil(:)          & ! out
        & )
      !$claw expand parallel update
      fract_snow(:)  = fract_snow_soil(:)
   END IF

   !$claw acc exit data delete(w_skin_max, w_skin_canopy_max)  
   !$claw verbatim  task_counter = task_counter + 1
  END SUBROUTINE update_snow_and_skin_fraction

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "snow_and_skin_fraction"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_snow_and_skin_fraction(tile, options)

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_seb_config), POINTER :: seb__conf
    TYPE(t_hydro_memory), POINTER :: hydro__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_snow_and_skin_fraction'

    INTEGER  :: iblk , ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    SELECT TYPE (cc => model%processes(seb_)%p%config); TYPE IS (t_seb_config) ; seb__conf => cc ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, hydro__mem%fract_water,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%fract_snow,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%fract_snow_soil,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, hydro__mem%fract_snow_can,ics,ice,iblk,routine)
    IF (seb__conf%l_ice_on_lakes) THEN
      CALL weighted_by_fract%Aggregate(tile, hydro__mem%fract_snow_lice,ics,ice,iblk,routine)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_snow_and_skin_fraction

  SUBROUTINE update_water_balance( tile, options)

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_grid),  POINTER    :: grid

    TYPE(t_hydro_memory), POINTER :: hydro__mem
    TYPE(t_a2l_memory), POINTER :: a2l__mem

    ! Pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: &
      & rain, &
      & snow, &
      & runoff, &
      & drainage, &
      & evapotrans, &
      & w_skin, &
      & w_snow, &
      & water_flux, &    ! [m3 s-1]
      & water_content, & ! [m3]
      & water_budget     ! [m3]
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & w_soil_sl, &
      & w_ice_sl

    REAL(wp), POINTER :: &
      & area(:)

    ! Locally allocated vectors
    !
    REAL(wp) :: &
      & tile_fract (options%nc), &
      & water_content_new(options%nc)

    INTEGER  :: iblk, ics, ice
    REAL(wp) :: dtime

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_water_balance'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    dtime = options%dtime

    ! IF (tile%lcts() RETURN   ! Check water balance only on root tile
    ! IF (tile%is_lake) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    area => grid%area (ics:ice, iblk)

    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT

    rain => a2l__mem%rain%ptr(ics:ice,iblk)                ! in
    snow => a2l__mem%snow%ptr(ics:ice,iblk)                ! in
    runoff => hydro__mem%runoff%ptr(ics:ice,iblk)              ! in
    drainage => hydro__mem%drainage%ptr(ics:ice,iblk)            ! in
    IF (tile%contains_soil .AND. .NOT. tile%is_glacier) THEN
      w_soil_sl => hydro__mem%w_soil_sl%ptr(ics:ice,:,iblk)           ! in
      w_ice_sl => hydro__mem%w_ice_sl%ptr(ics:ice,:,iblk)            ! in
      w_skin => hydro__mem%w_skin%ptr(ics:ice,iblk)              ! in
      w_snow => hydro__mem%w_snow%ptr(ics:ice,iblk)              ! in
    END IF
    evapotrans => hydro__mem%evapotrans%ptr(ics:ice,iblk)  ! in
    water_flux => hydro__mem%water_flux%ptr(ics:ice,iblk)          ! out
    water_content => hydro__mem%water_content%ptr(ics:ice,iblk)       ! out
    water_budget => hydro__mem%water_budget%ptr(ics:ice,iblk)        ! out

    !$claw acc data create(tile_fract, water_content_new)
    CALL tile%Get_fraction(ics, ice, iblk, fract=tile_fract(:))
    !$claw acc update device(tile_fract)
    !$claw expand parallel update savepoint=ser_jsb_hydro_update_water_balance_expand_0 tile_name=tile%name
    water_flux(:) = rain(:) + snow(:) + evapotrans(:) - runoff(:) - drainage(:)
    water_flux(:) = water_flux(:) * tile_fract(:) * area(:) / rhoh2o                    ! kg m-2 s-1 -> m3/s
    !$claw end expand

    IF (tile%contains_soil .AND. .NOT. tile%is_glacier) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_water_balance_expand_1 tile_name=tile%name
      water_content_new(:) = (  w_skin(:) + w_snow(:)                                   &
        &                     + SUM(w_soil_sl(:,:), DIM=2) + SUM(w_ice_sl(:,:), DIM=2)  &
        &                    ) * tile_fract(:) * area(:)  ! m -> m^3
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_water_balance_expand_2 tile_name=tile%name
      water_content_new(:) = 0._wp
    END IF

    IF (.NOT. is_time_experiment_start()) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_water_balance_expand_3 tile_name=tile%name
      water_budget(:) = water_content(:) + water_flux(:) * dtime - water_content_new(:)
    END IF

    !$claw expand parallel update savepoint=ser_jsb_hydro_update_water_balance_expand_4 tile_name=tile%name
    water_content(:) = water_content_new(:)

    !$claw acc end data

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_water_balance

  SUBROUTINE aggregate_water_balance(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    !TYPE(t_hydro_memory), POINTER :: hydro__mem

    !CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_water_balance'

    INTEGER :: iblk !, ics, ice

    iblk = options%iblk
    !ics  = options%ics
    !ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    !SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    ! weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    ! Don't aggregate, but explicitely compute water balance on each tile
    CALL update_water_balance(tile, options)
    ! CALL weighted_by_fract%Aggregate(tile, hydro__mem%water_flux,ics,ice,iblk,routine)
    ! CALL weighted_by_fract%Aggregate(tile, hydro__mem%water_content,ics,ice,iblk,routine)
    ! CALL weighted_by_fract%Aggregate(tile, hydro__mem%water_budget,ics,ice,iblk,routine)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_water_balance

#endif
END MODULE mo_hydro_interface
