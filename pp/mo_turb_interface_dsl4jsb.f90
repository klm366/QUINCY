!
! This file was automatically generated from:
!    /scratch/snx3000/rschnur/icon-rc-testing-plus-mc8.test/externals/jsbach/src/turbulence/mo_turb_interface.f90
!    with md5 hash: d4ff71ba345892354090c04e0b711ea3
!
!>
!! @brief Contains the interfaces to the turbulence process
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2016-04-01)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_turb_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  USE mo_jsb_process_class, ONLY: turb_, seb_, hydro_, pheno_, assimi_, sse_, a2l_

  ! Use process configurations
  !USE mo_turb_config_class, ONLY: t_turb_config
  USE mo_hydro_config_class, ONLY: t_hydro_config
  USE mo_seb_config_class, ONLY: t_seb_config

  ! Use process memories
  USE mo_a2l_memory_class, ONLY: t_a2l_memory
  USE mo_turb_memory_class, ONLY: t_turb_memory
  USE mo_seb_memory_class, ONLY: t_seb_memory
  USE mo_hydro_memory_class, ONLY: t_hydro_memory
  USE mo_pheno_memory_class, ONLY: t_pheno_memory
  USE mo_assimi_memory_class, ONLY: t_assimi_memory

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_turb_tasks !,t_turb_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_turb_interface'

  !> Type definition for surface_roughness task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_roughness
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_roughness     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_roughness  !< Aggregates computed task variables
  END TYPE tsk_surface_roughness

  !> Constructor interface for surface_roughness task
  INTERFACE tsk_surface_roughness
    PROCEDURE Create_task_surface_roughness                        !< Constructor function for task
  END INTERFACE tsk_surface_roughness

  !> Type definition for humidity_scaling task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_humidity_scaling
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_humidity_scaling     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_humidity_scaling  !< Aggregates computed task variables
  END TYPE tsk_humidity_scaling

  !> Constructor interface for humidity scaling task
  INTERFACE tsk_humidity_scaling
    PROCEDURE Create_task_humidity_scaling                        !< Constructor function for task
  END INTERFACE tsk_humidity_scaling

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_roughness task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "tsk_surface_roughness"
  !!
  FUNCTION Create_task_surface_roughness(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_roughness::return_ptr)
    CALL return_ptr%Construct(name='surface_roughness', process_id=TURB_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_roughness

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for humidity_scaling task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "tsk_humidity_scaling"
  !!
  FUNCTION Create_task_humidity_scaling(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_humidity_scaling::return_ptr)
    CALL return_ptr%Construct(name='humidity_scaling', process_id=TURB_, owner_model_id=model_id)

  END FUNCTION Create_task_humidity_scaling

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for turbulence process
  !!
  !! @param[in]     model_id  Model id
  !! @param[in,out] this      Instance of TURB process class
  !!
  SUBROUTINE Register_turb_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_surface_roughness(model_id))
    CALL this%Register_task(tsk_humidity_scaling(model_id))

  END SUBROUTINE Register_Turb_Tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "surface_roughness"
  !! Task "<PROCESS_NAME_LOWER_CASE>" <EXPLANATIION_OF_PROCESS>.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_roughness(tile, options)

    USE mo_turb_constants, ONLY: blending_height, roughness_snow

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_seb_config), POINTER :: seb__conf

    ! Declare process memories
    TYPE(t_turb_memory), POINTER :: turb__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem

    ! Declare pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: rough_m
    REAL(wp), POINTER, DIMENSION(:) :: rough_h
    REAL(wp), POINTER, DIMENSION(:) :: fract_snow_soil
    REAL(wp), POINTER, DIMENSION(:) :: fract_snow_lice

    ! Local variables
    TYPE(t_jsb_model), POINTER :: model
    INTEGER  :: iblk, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_roughness'

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(TURB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    SELECT TYPE (cc => model%processes(seb_)%p%config); TYPE IS (t_seb_config) ; seb__conf => cc ; END SELECT

    ! Get process memories
    SELECT TYPE (mm => tile%mem(turb_)%p); TYPE IS (t_turb_memory); turb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    ! Get process variables
    rough_m => turb__mem%rough_m%ptr(ics:ice,iblk)         ! IN
    rough_h => turb__mem%rough_h%ptr(ics:ice,iblk)         ! OUT

    ! RMS TBD: This has to be revised!
    IF (tile%Is_process_active(SSE_)) THEN
      fract_snow_soil => hydro__mem%fract_snow_soil%ptr(ics:ice,iblk) ! IN
      !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_0 tile_name=tile%name
      rough_h(:) =          fract_snow_soil(:)  / (LOG(blending_height /        roughness_snow)  ** 2)  &
        &        + (1._wp - fract_snow_soil(:)) / (LOG(blending_height / MIN(1._wp, rough_m(:))) ** 2)
    ELSE IF (tile%contains_lake) THEN
      IF (seb__conf%l_ice_on_lakes) THEN
        fract_snow_lice => hydro__mem%fract_snow_lice%ptr(ics:ice,iblk) ! IN
        !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_1 tile_name=tile%name
        rough_h(:) =          fract_snow_lice(:)  / (LOG(blending_height /        roughness_snow)  ** 2)  &
          &        + (1._wp - fract_snow_lice(:)) / (LOG(blending_height / MIN(1._wp, rough_m(:))) ** 2)
      ELSE
        !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_2 tile_name=tile%name
        rough_h(:) = LOG(blending_height / MIN(1._wp, rough_m(:))) ** 2
      END IF
    ELSE
      CALL finish(TRIM(routine), 'Don"t know what to do.')
    END IF

    ! Limit roughness length for heat by 1 m (see below Eq.5.6 in ECHAM5 manual)
    !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_3 tile_name=tile%name
    rough_h(:) = MIN(1._wp, blending_height * EXP(-1._wp / SQRT(rough_h(:))))

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_roughness

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_roughness"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_roughness(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_turb_memory), POINTER :: turb__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_roughness'

    ! Local variables
    INTEGER  :: iblk, ics, ice

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(turb_)%p); TYPE IS (t_turb_memory); turb__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, turb__mem%rough_h,ics,ice,iblk,routine)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_roughness

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "humidity_scaling"
  !! Task "<PROCESS_NAME_LOWER_CASE>" <EXPLANATIION_OF_PROCESS>.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_humidity_scaling(tile, options)

    !USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_sse_process,            ONLY: relative_humidity_soil

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_hydro_config), POINTER :: hydro__conf

    ! Declare pointers to process configs and memory
    TYPE(t_turb_memory), POINTER :: turb__mem
    TYPE(t_seb_memory), POINTER :: seb__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem
    TYPE(t_pheno_memory), POINTER :: pheno__mem
    TYPE(t_assimi_memory), POINTER :: assimi__mem
    TYPE(t_a2l_memory), POINTER :: a2l__mem

    ! Declare pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: &
      & fract_fpc,            &
      & fract_snow,           &
      & fract_water,          &
      & max_moist,            &
      & w_soil_column,        &
      & fact_q_air,           &
      & fact_qsat_srf,        &
      & fact_qsat_trans_srf,  &
      & qsat_star,            &
      & canopy_cond_limited,  &
      & pch,                  &
      & wind_air,             &
      & q_air

    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & w_soil_sl,            &
      & w_soil_fc_sl,         &
      & w_soil_pwp_sl

    ! Local variables
    REAL(wp) ::                    &
      fact_qsat_veg(options%nc),   &
      fact_qsat_trans(options%nc), &
      fact_qsat_soil(options%nc),  &
      fact_qair_veg(options%nc),   &
      fact_qair_soil(options%nc),  &
      rel_hum(options%nc),         & ! relative humidity of uppermost soil layer as compared to its field capacity
      zfract_fpc(options%nc),      &
      w_soil_wilt_fract_tmp
    LOGICAL :: &
      mask(options%nc)

    INTEGER :: iblk, ics, ice, i

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_humidity_scaling'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (.NOT. tile%Is_process_active(TURB_) .OR. tile%contains_lake) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    SELECT TYPE (cc => model%processes(hydro_)%p%config); TYPE IS (t_hydro_config) ; hydro__conf => cc ; END SELECT

    !$claw acc data create(fact_qsat_veg,fact_qsat_trans,fact_qsat_soil, & 
    !$claw acc fact_qair_veg,fact_qair_soil,rel_hum,zfract_fpc,mask)

    ! Get pointers to process configs and memory
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(turb_)%p); TYPE IS (t_turb_memory); turb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT
    IF (tile%contains_vegetation) THEN
      SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT
      IF (tile%Is_process_active(ASSIMI_)) THEN
        SELECT TYPE (mm => tile%mem(assimi_)%p); TYPE IS (t_assimi_memory); assimi__mem => mm ; END SELECT
      END IF
    END IF

    ! Set pointers to variables in memory
    pch => a2l__mem%pch%ptr(ics:ice,iblk)                 ! in
    wind_air => a2l__mem%wind_air%ptr(ics:ice,iblk)            ! in
    q_air => a2l__mem%q_air%ptr(ics:ice,iblk)               ! in
    fract_snow => hydro__mem%fract_snow%ptr(ics:ice,iblk)          ! in
    fact_q_air => turb__mem%fact_q_air%ptr(ics:ice,iblk)          ! OUT
    fact_qsat_srf => turb__mem%fact_qsat_srf%ptr(ics:ice,iblk)       ! OUT
    fact_qsat_trans_srf => turb__mem%fact_qsat_trans_srf%ptr(ics:ice,iblk) ! OUT
    IF (tile%contains_soil) THEN
      qsat_star => seb__mem%qsat_star%ptr(ics:ice,iblk)           ! in
      fract_water => hydro__mem%fract_water%ptr(ics:ice,iblk)         ! in
      max_moist => hydro__mem%max_moist%ptr(ics:ice,iblk)           ! in
      w_soil_column => hydro__mem%w_soil_column%ptr(ics:ice,iblk)       ! in
      w_soil_sl => hydro__mem%w_soil_sl%ptr(ics:ice,:,iblk)           ! in
      w_soil_fc_sl => hydro__mem%w_soil_fc_sl%ptr(ics:ice,:,iblk)        ! in
      w_soil_pwp_sl => hydro__mem%w_soil_pwp_sl%ptr(ics:ice,:,iblk)       ! in
    END IF

    ! ---------------------------
    ! Go

    IF (tile%contains_vegetation) THEN
      fract_fpc => pheno__mem%fract_fpc%ptr(ics:ice,iblk)           ! in
      IF (tile%Is_process_active(ASSIMI_)) THEN
        canopy_cond_limited => assimi__mem%canopy_cond_limited%ptr(ics:ice,iblk)         ! in
      ELSE
        canopy_cond_limited => hydro__mem%canopy_cond_limited%ptr(ics:ice,iblk)         ! in
      END IF

      IF (model%config%l_compat401) THEN
        !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_0 tile_name=tile%name
        mask(:) = w_soil_sl(:,1) > w_soil_pwp_sl(:,1)
      ELSE
        w_soil_wilt_fract_tmp = hydro__conf%w_soil_wilt_fract
        !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_1 tile_name=tile%name
        mask(:) = w_soil_column(:) > w_soil_wilt_fract_tmp * max_moist(:)
      END IF

      !$claw acc data present(fract_snow,fract_water,pch,wind_air,canopy_cond_limited,fact_qsat_veg,fact_qsat_trans,mask)
      !$claw acc parallel loop gang vector
      DO i = 1, size(mask, dim=1)
        IF (mask(i)) THEN
          fact_qsat_veg(i) = fract_snow(i) + (1._wp - fract_snow(i)) *   &
            & (fract_water(i) + ( 1._wp - fract_water(i)) /                &
            & (1._wp + pch(i) * MAX(1._wp, wind_air(i)) / MAX(1.E-20_wp, canopy_cond_limited(i))))
          fact_qsat_trans(i) = (1._wp - fract_snow(i)) * &
            & (1._wp - fract_water(i)) /                  &
            & (1._wp + pch(i) * MAX(1._wp, wind_air(i)) / MAX(1.E-20_wp, canopy_cond_limited(i)))
        ELSE 
          fact_qsat_veg(i)   = fract_snow(i) + (1._wp - fract_snow(i)) * fract_water(i)
          fact_qsat_trans(i) = 0._wp
        END IF
      END DO
      !$claw acc end data

      !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_2 tile_name=tile%name
      fact_qair_veg(:) = fact_qsat_veg(:)
      zfract_fpc(:)      = fract_fpc(:)
      !$claw end expand
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_3 tile_name=tile%name
      fact_qsat_veg(:)   = 0._wp
      fact_qsat_trans(:) = 0._wp
      fact_qair_veg(:)   = 0._wp
      zfract_fpc(:)      = 0._wp
      !$claw end expand
    END IF

    IF (tile%contains_soil) THEN
      !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_4 tile_name=tile%name
      rel_hum(:) = relative_humidity_soil(w_soil_sl(:,1), w_soil_fc_sl(:,1))

      !$claw acc data present(fact_qsat_soil,fact_qair_soil,qsat_star,rel_hum,q_air,fract_snow,fract_water)
      !$claw acc parallel loop gang vector
      DO i = 1, size(qsat_star,dim=1)
        IF (qsat_star(i) * rel_hum(i) > q_air(i) .AND. rel_hum(i) > 1.e-10_wp) THEN
          fact_qsat_soil(i) = fract_snow(i) + (1._wp - fract_snow(i)) * &
            & (fract_water(i) + (1._wp - fract_water(i)) * rel_hum(i))
          fact_qair_soil(i) = 1._wp
        ELSE
          fact_qsat_soil(i) = fract_snow(i) + (1._wp - fract_snow(i)) * fract_water(i)
          fact_qair_soil(i) = fact_qsat_soil(i)
        END IF
        IF (q_air(i) > qsat_star(i)) THEN
          fact_qsat_soil(i) = 1._wp
          fact_qair_soil(i) = 1._wp
        END IF
      END DO
      !$claw acc end data

    ELSE IF (tile%contains_glacier) THEN
      !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_5 tile_name=tile%name
      fact_qsat_soil(:) = 1._wp
      fact_qair_soil(:) = 1._wp
      !$claw end expand
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_6 tile_name=tile%name
      fact_qsat_soil(:) = 0._wp
      fact_qair_soil(:) = 0._wp
      !$claw end expand
    END IF

    !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_7 tile_name=tile%name
    fact_qsat_srf(:) = zfract_fpc(:) * fact_qsat_veg(:) + (1._wp - zfract_fpc(:)) * fact_qsat_soil(:)
    fact_q_air(:)    = zfract_fpc(:) * fact_qair_veg(:) + (1._wp - zfract_fpc(:)) * fact_qair_soil(:)
    fact_qsat_trans_srf(:) = zfract_fpc(:) * fact_qsat_trans(:)
    !$claw end expand

    !$claw acc end data

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_humidity_scaling

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "humidity_scaling"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_humidity_scaling(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_turb_memory), POINTER :: turb__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_humidity_scaling'

    ! Local variables
    INTEGER  :: iblk, ics, ice

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(turb_)%p); TYPE IS (t_turb_memory); turb__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, turb__mem%fact_q_air,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, turb__mem%fact_qsat_srf,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, turb__mem%fact_qsat_trans_srf,ics,ice,iblk,routine)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_humidity_scaling

  ! ================================================================================================================================
  !>
  !! @brief Another task for "<PROCESS_NAME_LOWER_CASE>"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !  ....

#endif
END MODULE mo_turb_interface
