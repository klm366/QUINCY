!
! This file was automatically generated from JSBACH source file:
!    src/phenology/mo_pheno_interface.f90
! with md5 hash: df947f7cbff275dedf30229ecd90c535
!
!>
!! @brief Contains the interfaces to the phenlogy processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-09-09)
!! Julia Nabel (2018-01-10): added optional functionality: forest regrowth (l_forestRegrowth)
!!                           (derive maxLAI from biomass via allometric relationships)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_pheno_interface
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

  ! Use of processes in this module (USE mo_jsb_process_class, ONLY: )
  USE mo_jsb_process_class, ONLY: seb_, a2l_, hydro_, assimi_, pheno_, carbon_

  ! Use of process configurations (t_<process>_config)
  USE mo_pheno_config_class, ONLY: t_pheno_config

  ! Use of process memories (t_<process>_memory)
  USE mo_a2l_memory_class, ONLY: t_a2l_memory
  USE mo_seb_memory_class, ONLY: t_seb_memory
  USE mo_hydro_memory_class, ONLY: t_hydro_memory
  USE mo_pheno_memory_class, ONLY: t_pheno_memory
  USE mo_assimi_memory_class, ONLY: t_assimi_memory
  USE mo_carbon_memory_class, ONLY: t_carbon_memory ! JN 10.01.18 only required with l_forestRegrowth = true

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_pheno_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_interface'

  !> Type definition for update_phenology task
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_phenology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_phenology    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_phenology !< Aggregates computed task variables
  END TYPE   tsk_phenology

  !> Constructor interface for update_phenology task
  INTERFACE   tsk_phenology
    PROCEDURE Create_task_phenology         !< Constructor function for task
  END INTERFACE   tsk_phenology

  !> Type definition for foliage_projected_cover task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_foliage_projected_cover
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_fpc
    PROCEDURE, NOPASS :: Aggregate => aggregate_fpc
  END TYPE tsk_foliage_projected_cover

  !> Constructor interface for foliage_projected_cover task
  INTERFACE tsk_foliage_projected_cover
    PROCEDURE Create_task_foliage_projected_cover
  END INTERFACE tsk_foliage_projected_cover

CONTAINS
  ! ================================================================================================================================
  !! Constructors for tasks

  !> Constructor for update_phenology task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "update_phenology"
  !!
  FUNCTION Create_task_phenology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_phenology::return_ptr)
    CALL return_ptr%Construct(name='phenology', process_id=PHENO_, owner_model_id=model_id)

  END FUNCTION Create_task_phenology

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for foliage_projected_cover task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "foliage_projected_cover"
  !!
  FUNCTION Create_task_foliage_projected_cover(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_foliage_projected_cover::return_ptr)
    CALL return_ptr%Construct(name='foliage_projected_cover', process_id=PHENO_, owner_model_id=model_id)

  END FUNCTION Create_task_foliage_projected_cover

  ! =======================================================================================================
  !> Register tasks for pheno process
  !!
  !! @param[in,out] this      Instance of pheno process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_pheno_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_phenology(model_id))
    CALL this%Register_task(tsk_foliage_projected_cover(model_id))

  END SUBROUTINE Register_pheno_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation to update the phenology state
  !! Task "phenology" calculates the phenology state of plants.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE update_phenology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER             :: model
    TYPE(t_pheno_config), POINTER :: pheno__conf

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_phenology'

    IF (.NOT. tile%Is_process_active(PHENO_)) RETURN
    model => Get_model(tile%owner_model_id)
    SELECT TYPE (cc => model%processes(pheno_)%p%config); TYPE IS (t_pheno_config) ; pheno__conf => cc ; END SELECT

    IF (tile%lcts(1)%lib_id == 0 .AND. TRIM(pheno__conf%scheme) /= 'climatology') THEN
      CALL message(TRIM(routine), 'Vegetation tile has no PFTs, using climatology for LAI')
      CALL update_phenology_climatology(tile, options)
    ELSE
      SELECT CASE (TRIM(pheno__conf%scheme))
      CASE ('logrop')
        CALL update_phenology_logrop(tile, options)
      CASE ('climatology')
        CALL update_phenology_climatology(tile, options)
      CASE DEFAULT
        CALL finish('Invalid scheme for phenology')
      END SELECT
    END IF

  END SUBROUTINE update_phenology

  ! ================================================================================================================================
  !>
  !! @brief Implementation of 'LOGROP' phenology scheme
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  !
  SUBROUTINE update_phenology_logrop(tile, options)

    ! Declarations
    USE mo_carbon_interface,     ONLY: rescale_carbon_upon_reference_area_change, &
      &                                calculate_current_c_ta_state_sum, check_carbon_conservation
    USE mo_pheno_process,        ONLY: calc_summergreen_phenology, calc_evergreen_phenology, &
      &                                calc_crop_phenology, calc_raingreen_phenology, calc_grass_phenology,         &
      &                                track_max_green_pool, derive_maxlai_according_to_allometric_relationships
    USE mo_jsb_time_iface,       ONLY: get_year_day
    USE mo_jsb_time,             ONLY: is_newyear, is_newday, timestep_in_days
    USE mo_pheno_constants,      ONLY: PhenoParam
    USE mo_jsb_grid_class,       ONLY: t_jsb_grid
    USE mo_jsb_grid,             ONLY: Get_grid
    USE mo_physical_constants,   ONLY: tmelt

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER :: model
    INTEGER                    :: iblk, ics, ice

    INTEGER :: day_of_year

    REAL(wp), DIMENSION(options%nc) :: &
      & w_soil_filling,                &
      & t_air_in_Celcius,              &
      & LaiMax_corrected,              &
      & veg_fract_correction_old,      &
      & old_c_state_sum, dummy_flux

    TYPE(t_jsb_grid), POINTER :: grid

    ! Declare pointers for process configuration and memory
    TYPE(t_pheno_config), POINTER :: pheno__conf
    TYPE(t_pheno_memory), POINTER :: pheno__mem
    TYPE(t_a2l_memory), POINTER :: a2l__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem
    TYPE(t_assimi_memory), POINTER :: assimi__mem
    TYPE(t_seb_memory), POINTER :: seb__mem
    TYPE(t_carbon_memory), POINTER :: carbon__mem ! JN 10.01.18 only required with l_forestRegrowth = true

    ! Declare pointers for variables in memory
    !INTEGER                :: ncanopy     ! number of canopy layers

    !INTEGER, POINTER                 :: acc_counter ! counter for each delta_time step since the last variable output
    !LOGICAL, POINTER                 :: use_alb_veg_simple
    !LOGICAL, POINTER                 :: first_call_phenology_this_day

    REAL(wp), POINTER, DIMENSION(:) :: t_air
    REAL(wp), POINTER, DIMENSION(:) :: npp_pot_rate_ca
    REAL(wp), POINTER, DIMENSION(:) :: w_soil_rel
    REAL(wp), POINTER, DIMENSION(:) :: lat
    REAL(wp), POINTER, DIMENSION(:) :: heat_sum_SG
    REAL(wp), POINTER, DIMENSION(:) :: heat_sum_EG
    REAL(wp), POINTER, DIMENSION(:) :: heat_sum_CRP
    REAL(wp), POINTER, DIMENSION(:) :: heat_sum_winter
    REAL(wp), POINTER, DIMENSION(:) :: days_since_growth_begin_SG
    REAL(wp), POINTER, DIMENSION(:) :: days_since_growth_begin_EG
    REAL(wp), POINTER, DIMENSION(:) :: growth_phase_SG
    REAL(wp), POINTER, DIMENSION(:) :: growth_phase_EG
    REAL(wp), POINTER, DIMENSION(:) :: growth_phase_CRP
    REAL(wp), POINTER, DIMENSION(:) :: previous_day_temp_mean
    REAL(wp), POINTER, DIMENSION(:) :: pseudo_soil_temp
    REAL(wp), POINTER, DIMENSION(:) :: chill_days_SG
    REAL(wp), POINTER, DIMENSION(:) :: chill_days_EG
    REAL(wp), POINTER, DIMENSION(:) :: lai
    REAL(wp), POINTER, DIMENSION(:) :: day_npp_sum
    REAL(wp), POINTER, DIMENSION(:) :: previous_day_temp_min
    REAL(wp), POINTER, DIMENSION(:) :: previous_day_temp_max
    REAL(wp), POINTER, DIMENSION(:) :: previous_day_npp_pot_rate_ca

    ! JN 10.01.18 only required with l_forestRegrowth = true
    REAL(wp), POINTER, DIMENSION(:) :: c_green
    REAL(wp), POINTER, DIMENSION(:) :: c_reserve
    REAL(wp), POINTER, DIMENSION(:) :: c_woods

    REAL(wp), POINTER, DIMENSION(:) :: current_max_green
    REAL(wp), POINTER, DIMENSION(:) :: veg_carbon_at_max_green

    REAL(wp), POINTER, DIMENSION(:) :: number_of_individuals
    REAL(wp), POINTER, DIMENSION(:) :: biomass_per_individual
    REAL(wp), POINTER, DIMENSION(:) :: maxlai_allom
    REAL(wp), POINTER, DIMENSION(:) :: veg_fract_correction

    REAL(wp), POINTER, DIMENSION(:) :: cconservation_allom

    ! Local variables
    LOGICAL :: l_newday
    REAL(wp) :: timestep  ! Timestep in days
    INTEGER :: i
    REAL(wp) :: MaxLAI
    REAL(wp) :: LAI_negligible_tmp


    CHARACTER(len=*), PARAMETER  :: routine = modname//':update_phenology_logrop'

    ! Set local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(PHENO_)) RETURN

    !$claw acc data create(w_soil_filling, t_air_in_Celcius)  & 
    !$claw create(LaiMax_corrected, veg_fract_correction_old) &
    !$claw create(old_c_state_sum, dummy_flux)
    IF ( .NOT. tile%Is_process_active(ASSIMI_)) THEN
      CALL finish(routine, 'The LOGROP phenology model depends on the assimilation process being activated, too')
    ENDIF

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Set model
    model => Get_model(tile%owner_model_id)

    ! Set pointers to process configs and memory
    SELECT TYPE (cc => model%processes(pheno_)%p%config); TYPE IS (t_pheno_config) ; pheno__conf => cc ; END SELECT
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(assimi_)%p); TYPE IS (t_assimi_memory); assimi__mem => mm ; END SELECT

    !ncanopy =  assimi__conf%ncanopy

    grid    => Get_grid(model%grid_id)
    lat     => grid%lat(ics:ice, iblk)

    !first_call_phenology_this_day  => model%processes(PHENO_)%config%first_call_phenology_this_day(ics:ice,iblk)
    IF (pheno__conf%l_forestRegrowth) THEN

      !ASSERTION: l_forestRegrowth only possible when the CARBON process is active
      IF ( .NOT. tile%Is_process_active(CARBON_)) THEN
        CALL finish(TRIM(routine), 'Forest regrowth can only be switched on when the CARBON process is active, too')
      ENDIF

      SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT
      c_green => carbon__mem%c_green%ptr(ics:ice,iblk)             ! in
      c_reserve => carbon__mem%c_reserve%ptr(ics:ice,iblk)           ! in
      c_woods => carbon__mem%c_woods%ptr(ics:ice,iblk)             ! in

      current_max_green => carbon__mem%current_max_green%ptr(ics:ice,iblk)       ! inout
      veg_carbon_at_max_green => carbon__mem%veg_carbon_at_max_green%ptr(ics:ice,iblk) ! inout

      number_of_individuals => pheno__mem%number_of_individuals%ptr(ics:ice,iblk)   ! inout
      biomass_per_individual => pheno__mem%biomass_per_individual%ptr(ics:ice,iblk)  ! inout
      maxlai_allom => pheno__mem%maxlai_allom%ptr(ics:ice,iblk)            ! inout
      veg_fract_correction => pheno__mem%veg_fract_correction%ptr(ics:ice,iblk)    ! inout

      cconservation_allom => pheno__mem%cconservation_allom%ptr(ics:ice,iblk)     ! inout

    END IF

    ! Set pointers to variables in memory
    t_air => a2l__mem%t_air%ptr(ics:ice,iblk)                  ! in
    npp_pot_rate_ca => assimi__mem%npp_pot_rate_ca%ptr(ics:ice,iblk)        ! in
    w_soil_rel => hydro__mem%w_soil_rel%ptr(ics:ice,iblk)             ! in

    previous_day_temp_mean => seb__mem%previous_day_temp_mean%ptr(ics:ice,iblk) ! in
    previous_day_temp_min => seb__mem%previous_day_temp_min%ptr(ics:ice,iblk)  ! in
    previous_day_temp_max => seb__mem%previous_day_temp_max%ptr(ics:ice,iblk)  ! in
    pseudo_soil_temp => seb__mem%pseudo_soil_temp%ptr(ics:ice,iblk)       ! in

    lai => pheno__mem%lai%ptr(ics:ice,iblk)                    ! inout
    heat_sum_SG => pheno__mem%heat_sum_SG%ptr(ics:ice,iblk)            ! inout
    heat_sum_EG => pheno__mem%heat_sum_EG%ptr(ics:ice,iblk)            ! inout
    heat_sum_CRP => pheno__mem%heat_sum_CRP%ptr(ics:ice,iblk)           ! inout
    heat_sum_winter => pheno__mem%heat_sum_winter%ptr(ics:ice,iblk)        ! inout

    days_since_growth_begin_SG => pheno__mem%days_since_growth_begin_SG%ptr(ics:ice,iblk) ! inout
    days_since_growth_begin_EG => pheno__mem%days_since_growth_begin_EG%ptr(ics:ice,iblk) ! inout


    growth_phase_SG => pheno__mem%growth_phase_SG%ptr(ics:ice,iblk)        ! inout
    growth_phase_EG => pheno__mem%growth_phase_EG%ptr(ics:ice,iblk)        ! inout
    growth_phase_CRP => pheno__mem%growth_phase_CRP%ptr(ics:ice,iblk)       ! inout
    chill_days_SG => pheno__mem%chill_days_SG%ptr(ics:ice,iblk)          ! inout
    chill_days_EG => pheno__mem%chill_days_EG%ptr(ics:ice,iblk)          ! inout

    day_npp_sum => pheno__mem%day_npp_sum%ptr(ics:ice,iblk)            ! inout
    previous_day_npp_pot_rate_ca => pheno__mem%previous_day_npp_pot_rate_ca%ptr(ics:ice,iblk) ! out

    ! ---------------------------
    ! Go

    ! R: NUR ZUR SICHERHEIT:
    !    Diese Abfrage sollte unnötig sein, da der Prozess Phänologie nur für veg tiles aufgerufen werden sollte.
    IF (.NOT. tile%is_vegetation) THEN
       CALL finish(TRIM(routine), 'The programm tries to calculate phenology on a non vegetated tile...')
    END IF

    ! Initialization of local variables
    l_newday = is_newday(options%current_datetime, options%dtime)
    timestep = timestep_in_days(options%dtime)

    !$claw expand parallel update
    t_air_in_Celcius(:) = t_air(:) - tmelt  ! convert Kelvin in Celcius

    day_of_year  = INT(get_year_day(options%current_datetime))

    !$claw expand parallel update
    w_soil_filling(:) = w_soil_rel(:)  ! Initialization of w_soil_filling for crops and raingreen
                                       ! R: As in JS3 all soil layers are set to the same relative_moisture(GP,tile) and only the average
                                       !    over all layers or the uppermost layer is used, we do not need soil layers in the phenology.

    ! R: Note, in JSBACH3 the timestep was 1 day. Except of the former update_growth_phase which is now in the calc_XXX_phenology
    !    procedures the phenology including get_letItDie and get_letItGrow can be called several times a day. Within the
    !    calc_XXX_phenology procedures the is_newday-variable decides when to do the update_growth_phase stuff.

    ! R: While the time issue is independent of the model or block, some calculations should not be done for each tile but only
    !    once for the whole hsm (tile%name == 'box'):
    ! R: Die Berechnungen von update_previous_day_variables in JSBACH3 muß man bis auf NPP nicht für jedes Tile machen. Daher
    !    habe ich die Berechnungen von NPP von den anderen Berechnungen in update_previous_day_variables für JSBACH4 getrennt.
    !    Die Berechnungen die nicht für jedes Tile gemacht werden müssen sind in calc_previous_day_variables und die NPP-
    !    Berechnungen werden im Folgenden vorgenommen. Außerdem ist calc_previous_day_variables nun dem SEB Prozess zugeordnet.
    !    calc_pseudo_soil_temp habe ich auch dem SEB Prozess zugeordnet.

    ! day_NPP_sum contains the accumulated sum of NPP_pot_rate_ca and is computed in update_assimilation each time step
    ! after this routine. At the start of a new day, we here divide the sum by the day length to get
    ! previous_day_NPP_pot_rate_ca (needed in this subroutine) and reset it to zero.
    IF (l_newday) THEN
      !$claw expand parallel update savepoint=ser_jsb_pheno_day_npp_sum_expand_0 tile_name=tile%name
      previous_day_NPP_pot_rate_ca(:) = day_NPP_sum(:)/86400._wp
      !$claw expand parallel update savepoint=ser_jsb_pheno_day_npp_sum_expand_1 tile_name=tile%name
      day_NPP_sum(:) = 0._wp
    END IF

    ! JN 10.01.18
    IF ((pheno__conf%l_forestRegrowth) .AND. (model%lctlib(tile%lcts(1)%lib_id)%ForestFlag)) THEN
      ! to represent forest regrowth LAI needs to be dependant on the available leaf carbon
      ! to keep the overall structure of pheno we mimic this dependency by deriving the maxLAI
      ! from the available biomass via allometric relationships
      ! The maxLai gets updated at the beginning of each year with the biomass at the maximum green pool
      IF(is_newyear(options%current_datetime, options%dtime)) THEN
        !$claw expand parallel update savepoint=ser_jsb_PHENO_update_phenology_logrop_expand_1 tile_name=tile%name
        veg_fract_correction_old(:) = veg_fract_correction(:)

        !$claw sca forward update savepoint=ser_jsb_PHENO_cderive_maxLai_according_to_allometric_relationships_1 tile_name=tile%name
        CALL derive_maxLai_according_to_allometric_relationships( &
          & model%lctlib(tile%lcts(1)%lib_id)%alpha_nr_ind,                   & ! in
          & model%lctlib(tile%lcts(1)%lib_id)%beta_nr_ind,                    & ! in
          & model%lctlib(tile%lcts(1)%lib_id)%alpha_leaf,                     & ! in
          & model%lctlib(tile%lcts(1)%lib_id)%beta_leaf,                      & ! in
          & model%lctlib(tile%lcts(1)%lib_id)%clumpinessFactor,               & ! in
          & model%lctlib(tile%lcts(1)%lib_id)%specificLeafArea_C,             & ! in
          & c_woods(:),                                           & ! in
          & current_max_green(:),                                 & ! InOut
          & veg_carbon_at_max_green(:),                           & ! InOut
          & number_of_individuals(:),                             & ! InOut
          & biomass_per_individual(:),                            & ! InOut
          & veg_fract_correction(:),                              & ! InOut
          & maxLAI_allom(:)                                       & ! InOut
          & )

        ! redistribute
        CALL calculate_current_c_ta_state_sum(tile, options, old_c_state_sum(:))

        CALL rescale_carbon_upon_reference_area_change(tile, options, veg_fract_correction_old, veg_fract_correction)

        !$claw expand parallel update
        dummy_flux(:) = 0._wp
        CALL check_carbon_conservation(tile, options, old_c_state_sum(:), &
          & dummy_flux(:), cconservation_allom(:))
      END IF

      ! carbon dynamics are only calculated on a day by day basis, thus update only required once a day
      IF(is_newday(options%current_datetime, options%dtime)) THEN
        !$claw sca forward update savepoint=ser_jsb_PHENO_track_max_green_pool_1 tile_name=tile%name
        CALL track_max_green_pool(      &
          & c_green(:),                 & ! in
          & c_reserve(:),               & ! in
          & c_woods(:),                 & ! in
          & current_max_green(:),       & ! InOut
          & veg_carbon_at_max_green(:)  & ! InOut
          & )
      ENDIF
    ENDIF


    ! JN 10.01.18: if l_forestRegrowth = true, the maxLai of forest pfts can be variable over time and space, therefore
    !              LaiMax_corrected is now a vector!
    ! R: laiMax_dyn was included to avoid longterm negative NPP. It is reduced if NPP is negative over a year and increased
    !    if NPP is positive (subroutine update_maximum_lai). So if NPP is negative the lai gets very small.
    !    As this scheme wasn't used, it is unnecessary to implement it into JSBACH4.
    !      (LaiMax in JSBACH3 is initialized for vegetation always with a value from lctlib. When calling update_phenology from
    !       theLand%Vegetation%lai_max the laiMax_stat is created. Then NOT used in update_phenology:
    !       IF (ANY(laiMax_dyn(...) < 0._dp)) laiMax_dyn(...) =
    !       laiMax_stat(...) CALL update_maximum_lai(lat(...),laiMax_stat(...),lai(...)) )
    !
    !    Therefore I replaced laiMax_dyn with LaiMax_corrected, which is not dynamic and only used to secure that
    !    the maximum lai is at least PhenoParam%all%LAI_negligible if lctlib%LaiMax < PhenoParam%all%LAI_negligible
    !    However, as lctlib%LaiMax will always be higher than PhenoParam%all%LAI_negligible we could skip this part. But to remind
    !    us that we should think about reimplementing a dynamic lai I keep this part in a way making reimplementation easy...
    !
    ! Set initially the dynamical maximum LAI to the static maximum LAI
    !IF (ANY(laiMax_dyn(1545,:) < 0._dp)) laiMax_dyn(1545,:) = laiMax_stat(1:kidx,:)
    ! Test if lai <= LaiMax_dyn
    ! IF (debug) THEN
    !   DO i=1,kidx
    !      k=kidx0+i-1
    !      DO j=1,ntiles
    !         IF(lai(i,j) > laiMax_dyn(k,j) + EPSILON(1.0_dp)) THEN
    !            CALL message("update_phenology", "lai: "//real2string(lai(i,j))//" laiMax: "//real2string(laiMax_dyn(k,j)))
    !            CALL finish("update_phenology()","ERROR: On input one LAI is larger than maximum LAI.")
    !         END IF
    !      END DO
    !   END DO
    ! ELSE
    !   IF (ANY(lai(1:kidx,1:ntiles) > laiMax_dyn(1545,1:ntiles))) &
    !        CALL finish("update_phenology","ERROR: LAI larger than maximum LAI on input")
    ! END IF

    IF (pheno__conf%l_forestRegrowth) THEN
       !$claw expand parallel update savepoint=ser_jsb_pheno_update_phenology_logrop_expand_2 tile_name=tile%name
       LaiMax_corrected(:) = pheno__mem%maxlai_allom%ptr(ics:ice,iblk)
    ELSE
     ! R: Rename MaxLAI in lctlib to LAI_max in JS3 comes from theLand%Vegetation%LAI_max,
     ! which is for each vegetated tile taken from MaxLAI in the lct library.
       MaxLAI = model%lctlib(tile%lcts(1)%lib_id)%MaxLAI
       !$claw expand parallel update savepoint=ser_jsb_pheno_update_phenology_logrop_expand_3 tile_name=tile%name
       LaiMax_corrected(:) = MaxLAI
    END IF

    ! R: wird jedes mal (d.h. für jedes Tile) aufgerufen, wobei automatisch der Tile spezifische Wert aus dem lctlib genommen wird:
    ! Secure minimum Max-LAI

    LAI_negligible_tmp = PhenoParam%all%LAI_negligible
    !$claw acc data present(LaiMax_corrected, lai)
    !$claw acc parallel loop gang vector
    DO i = 1, SIZE(LaiMax_corrected)
      IF (LaiMax_corrected(i) < LAI_negligible_tmp) THEN
        LaiMax_corrected(i) = LAI_negligible_tmp
        lai(i) = LAI_negligible_tmp - EPSILON(1.0_wp)
      END IF
    END DO
    !$claw acc end data

    ! R: The following PHENO_ process shall be calculated only on the lowest child tiles (=ON_LEAFS).
    !    This is set within the FUNCTION Construct_tile and makes further statements here unnecessary.

    ! none: 0; evergreen: 1; summergreen: 2; raingreen: 3; grasses: 4; crops: 5
    IF (model%lctlib(tile%lcts(1)%lib_id)%PhenologyType ==  2) THEN

      IF (debug_on() .AND. iblk == 1) CALL message('     ','... summergreen ...')
      !$claw sca forward update savepoint=ser_jsb_pheno_calc_summergreen_phenology_1 tile_name=tile%name
      CALL calc_summergreen_phenology(    &
         & day_of_year,                   & ! in
         & l_newday,                      & ! in
         & timestep,                      & ! in
         & lat(:),                        & ! in
         & previous_day_temp_mean(:),     & ! in
         & pseudo_soil_temp(:),           & ! in
         & LaiMax_corrected(:),           & ! in
         & days_since_growth_begin_SG(:), & ! inout
         & chill_days_SG(:),              & ! inout
         & heat_sum_SG(:),                & ! inout
         & growth_phase_SG(:),            & ! inout
         & lai(:)                         & ! inout
         & )

    END IF

    ! none: 0; evergreen: 1; summergreen: 2; raingreen: 3; grasses: 4; crops: 5
    IF (model%lctlib(tile%lcts(1)%lib_id)%PhenologyType ==  1) THEN
      IF (debug_on() .AND. iblk == 1) CALL message('     ','... evergreen ...')
       !$claw sca forward update savepoint=ser_jsb_pheno_calc_evergreen_phenology_1 tile_name=tile%name
       CALL calc_evergreen_phenology(     &
         & day_of_year,                   & ! in
         & l_newday,                      & ! in
         & timestep,                      & ! in
         & lat(:),                        & ! in
         & pseudo_soil_temp(:),           & ! in
         & LaiMax_corrected(:),           & ! in
         & days_since_growth_begin_EG(:), & ! inout
         & chill_days_EG(:),              & ! inout
         & heat_sum_EG(:),                & ! inout
         & growth_phase_EG(:),            & ! inout
         & lai(:)                         & ! inout
         & )
    END IF

    ! none: 0; evergreen: 1; summergreen: 2; raingreen: 3; grasses: 4; crops: 5
    IF (model%lctlib(tile%lcts(1)%lib_id)%PhenologyType ==  5) THEN
       ! R: (Abweichungen zu JSBACH3 nur beim LAI und dort erst nach 9. Stelle nach Komma.)
      IF (debug_on() .AND. iblk == 1) CALL message('     ','... crop ...')
      !$claw sca forward update savepoint=ser_jsb_pheno_calc_crop_phenology_1 tile_name=tile%name
      CALL calc_crop_phenology(            &
        & day_of_year,                     & ! in
        & l_newday,                        & ! in
        & timestep,                        & ! in
        & NPP_pot_rate_ca(:),              & ! in
        & lat(:),                          & ! in
        & pseudo_soil_temp(:),             & ! in
        & LaiMax_corrected(:),             & ! in
        & previous_day_NPP_pot_rate_ca(:), & ! in
        & w_soil_filling(:),               & ! in
        & previous_day_temp_min(:),        & ! in
        & previous_day_temp_max(:),        & ! in
        & model%lctlib(tile%lcts(1)%lib_id)%specificLeafArea_C, & ! in
        & heat_sum_CRP(:),                 & ! inout
        & heat_sum_winter,                 & ! inout
        & growth_phase_CRP(:),             & ! inout
        & lai(:)                           & ! inout
        & )
    END IF

    ! none: 0; evergreen: 1; summergreen: 2; raingreen: 3; grasses: 4; crops: 5
    IF (model%lctlib(tile%lcts(1)%lib_id)%PhenologyType ==  3) THEN

      IF (debug_on() .AND. iblk == 1) CALL message('     ','... raingreen ...')
      !$claw sca forward update savepoint=ser_jsb_pheno_calc_raingreen_phenology_1 tile_name=tile%name
      CALL calc_raingreen_phenology(        &
        & timestep,                         & ! in
        & w_soil_filling(:),                & ! in
        & LaiMax_corrected(:),              & ! in
        & previous_day_NPP_pot_rate_ca(:),  & ! in
        & lai(:)                            & ! inout
        & )
    END IF

    ! none: 0; evergreen: 1; summergreen: 2; raingreen: 3; grasses: 4; crops: 5
    IF (model%lctlib(tile%lcts(1)%lib_id)%PhenologyType ==  4) THEN
      IF (debug_on() .AND. iblk == 1) CALL message('     ','... grasses ...')
       !$claw sca forward update savepoint=ser_jsb_pheno_calc_grass_phenology_1 tile_name=tile%name
       CALL calc_grass_phenology(           &
         & t_air_in_Celcius(:),             & ! in
         & timestep,                        & ! in
         & w_soil_filling(:),               & ! in
         & LaiMax_corrected(:),             & ! in
         & previous_day_NPP_pot_rate_ca(:), & ! in
         & lai(:)                           & ! inout
         & )
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

    !$claw acc end data

   END SUBROUTINE update_phenology_logrop

  ! ================================================================================================================================
  !>
  !! @brief Implementation of 'climatology' phenology scheme
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  !
  SUBROUTINE update_phenology_climatology(tile, options)

    USE mo_jsb_time, ONLY: get_time_interpolation_weights

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice
    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2

    TYPE(t_pheno_memory), POINTER :: pheno__mem

    ! Declare pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: lai
    REAL(wp), POINTER      :: lai_mon(:,:)

    CHARACTER(len=*), PARAMETER            :: routine = modname//':update_phenology_climatology'

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(PHENO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')


    SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT

    lai => pheno__mem%lai%ptr(ics:ice,iblk)                              ! out
    lai_mon(1:,0:) => pheno__mem%lai_mon (ics:ice, iblk, :) ! in

    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
    !$claw expand parallel update savepoint=ser_jsb_pheno_update_phenology_climatology_expand_1 tile_name=tile%name
    lai(:) = REAL(wgt1,wp) * lai_mon(:,nmw1) + REAL(wgt2,wp) * lai_mon(:,nmw2)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_phenology_climatology

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "phenology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_phenology(tile, options)

    USE mo_jsb_time,             ONLY: is_newyear

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER             :: model
    TYPE(t_pheno_config), POINTER :: pheno__conf
    TYPE(t_pheno_memory), POINTER :: pheno__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_phenology'

    ! Local variables
    INTEGER  :: iblk, ics, ice

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    SELECT TYPE (cc => model%processes(pheno_)%p%config); TYPE IS (t_pheno_config) ; pheno__conf => cc ; END SELECT
    SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, pheno__mem%lai,ics,ice,iblk,routine)

    IF (pheno__conf%l_forestRegrowth) THEN
      IF (is_newyear(options%current_datetime, options%dtime)) THEN
        CALL weighted_by_fract%Aggregate(tile, pheno__mem%maxlai_allom,ics,ice,iblk,routine)
        CALL weighted_by_fract%Aggregate(tile, pheno__mem%cconservation_allom,ics,ice,iblk,routine)

        !NOTE: despite changes, no aggregation of cpools necessary:
        ! with l_forestRegrowth PHENO changes Cpools only once a year at the beginning of a timestep,
        ! and other processes (NPP alloc of carbon) aggregate them
      ENDIF
    ENDIF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob u. wie die einzeln. pheno Tasks aggregiert werden sollen.

    IF (TRIM(pheno__conf%scheme) == 'logrop') THEN

    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_phenology

  SUBROUTINE update_fpc(tile, options)

    USE mo_pheno_process, ONLY: get_foliage_projected_cover
    USE mo_jsb_time,      ONLY: get_time_interpolation_weights

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    !----------------------------
    ! Local variables

    TYPE(t_pheno_config), POINTER :: pheno__conf
    TYPE(t_pheno_memory), POINTER :: pheno__mem

    ! Declare pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: lai
    REAL(wp), POINTER, DIMENSION(:) :: fract_fpc
    REAL(wp), POINTER, DIMENSION(:) :: fract_fpc_max
    REAL(wp), POINTER      :: fract_fpc_mon(:,:)


    INTEGER  :: iblk, ics, ice
    REAL(wp) :: wgt1, wgt2, ClumpinessFactor
    INTEGER  :: nmw1,nmw2

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_fpc'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (.NOT. tile%Is_process_active(PHENO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    SELECT TYPE (cc => model%processes(pheno_)%p%config); TYPE IS (t_pheno_config) ; pheno__conf => cc ; END SELECT
    SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT

    lai => pheno__mem%lai%ptr(ics:ice,iblk)            ! in
    fract_fpc_max => pheno__mem%fract_fpc_max%ptr(ics:ice,iblk)  ! in
    fract_fpc => pheno__mem%fract_fpc%ptr(ics:ice,iblk)      ! out

    ! Update foliage projected cover
    !
    IF (        TRIM(pheno__conf%scheme) == 'climatology' &
         & .OR. tile%lcts(1)%lib_id == 0) THEN         ! General vegetation tile, no specific PFT
      IF (model%config%l_compat401) THEN
        !$claw expand parallel update savepoint=ser_jsb_pheno_update_fpc_expand_1 tile_name=tile%name
        fract_fpc(:) = get_foliage_projected_cover(fract_fpc_max(:), lai(:), 1._wp)
      ELSE
        fract_fpc_mon(1:,0:) => pheno__mem%fract_fpc_mon (ics:ice, iblk, :) ! in
        CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
        !$claw expand parallel update savepoint=ser_jsb_pheno_update_fpc_expand_2 tile_name=tile%name
        fract_fpc(:) = REAL(wgt1,wp) * fract_fpc_mon(:,nmw1) + REAL(wgt2,wp) * fract_fpc_mon(:,nmw2)
      END IF
    ELSE        
      ! We have PFTs and %scheme /= climatology
      ClumpinessFactor = model%lctlib(tile%lcts(1)%lib_id)%ClumpinessFactor ! Do a copy to avoid deep copy on accelerator
      !$claw expand parallel update savepoint=ser_jsb_pheno_update_fpc_expand_3 tile_name=tile%name
      fract_fpc(:) = get_foliage_projected_cover(fract_fpc_max(:), lai(:), ClumpinessFactor)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_fpc

  SUBROUTINE aggregate_fpc(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_pheno_memory), POINTER :: pheno__mem

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_fpc'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, pheno__mem%fract_fpc,ics,ice,iblk,routine)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_fpc

#endif
END MODULE mo_pheno_interface
