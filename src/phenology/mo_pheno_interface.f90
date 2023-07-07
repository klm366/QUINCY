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

  USE mo_jsb_model_class,    ONLY: t_jsb_model, MODEL_JSBACH, MODEL_QUINCY
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module (USE mo_jsb_process_class, ONLY: )
  dsl4jsb_Use_processes SEB_, A2L_, HYDRO_, ASSIMI_, PHENO_, CARBON_, VEG_

  ! Use of process configurations (t_<process>_config)
  dsl4jsb_Use_config(PHENO_)

  ! Use of process memories (t_<process>_memory)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(VEG_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_pheno_tasks, global_phenology_diagnostics, merge_pheno_variables_upon_ageing

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_interface'

  !> Type definition for the phenology task
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

  ! ======================================================================================================= !
  !> Register tasks for pheno process
  !>
  SUBROUTINE Register_pheno_tasks(this, model_id)
    USE mo_jsb_model_class,   ONLY: t_jsb_model, MODEL_JSBACH
    USE mo_jsb_class,         ONLY: Get_model
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in) :: model_id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    ! ----------------------------------------------------------------------------------------------------- !
    model   => Get_model(model_id)

    CALL this%Register_task(tsk_phenology(model_id))
    ! register "foliage_projected_cover" only for the JSBACH model
    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_JSBACH)
      CALL this%Register_task(tsk_foliage_projected_cover(model_id))
    END SELECT
  END SUBROUTINE Register_pheno_tasks

  ! ======================================================================================================= !
  !>Implementation of "update" of the "phenology" task
  !> 
  !>  calls a model specific routine
  !> 
  SUBROUTINE update_phenology(tile, options)
    USE mo_jsb_model_class,         ONLY: t_jsb_model, MODEL_JSBACH, MODEL_QUINCY
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_pheno_update_phenology,  ONLY: update_phenology_quincy
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER             :: model
    dsl4jsb_Def_config(PHENO_)
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_phenology'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(PHENO_)) RETURN
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(PHENO_)

    ! call diffent routines depending on the model
    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_JSBACH)
      CALL update_phenology_jsbach(tile, options)
    CASE (MODEL_QUINCY)
      CALL update_phenology_quincy(tile, options)
    END SELECT
  END SUBROUTINE update_phenology

  ! ======================================================================================================= !
  !>
  !! @brief Implementation to update the phenology state
  !! Task "phenology" calculates the phenology state of plants.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE update_phenology_jsbach(tile, options)
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_class,               ONLY: Get_model
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER             :: model
    dsl4jsb_Def_config(PHENO_)
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_phenology_jsbach'
    ! ----------------------------------------------------------------------------------------------------- !
    model        => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(PHENO_)

    IF (tile%lcts(1)%lib_id == 0 .AND. TRIM(dsl4jsb_Config(PHENO_)%scheme_jsbach) /= 'climatology') THEN
      CALL message(TRIM(routine), 'Vegetation tile has no PFTs, using climatology for LAI')
      CALL update_phenology_climatology(tile, options)
    ELSE
      SELECT CASE (TRIM(dsl4jsb_Config(PHENO_)%scheme_jsbach))
      CASE ('logrop')
        CALL update_phenology_logrop(tile, options)
      CASE ('climatology')
        CALL update_phenology_climatology(tile, options)
      CASE DEFAULT
        CALL finish('Invalid scheme for phenology')
      END SELECT
    END IF
  END SUBROUTINE update_phenology_jsbach

  ! ======================================================================================================= !
  !>
  !! @brief Implementation of 'LOGROP' phenology scheme
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  !
  SUBROUTINE update_phenology_logrop(tile, options)

    ! Declarations
    USE mo_carbon_interface,        ONLY: rescale_carbon_upon_reference_area_change, &
      &                                   calculate_current_c_ta_state_sum, check_carbon_conservation, &
      &                                   recalc_carbon_per_tile_vars
    USE mo_pheno_process,           ONLY: calc_summergreen_phenology, calc_evergreen_phenology, &
      &                                   calc_crop_phenology, calc_raingreen_phenology, calc_grass_phenology,         &
      &                                   track_max_green_pool, derive_maxlai_according_to_allometric_relationships
    USE mo_jsb_time_iface,          ONLY: get_year_day
    USE mo_jsb_time,                ONLY: is_newyear, is_newday, timestep_in_days
    USE mo_pheno_parameters,        ONLY: pheno_param_jsbach
    USE mo_jsb_grid_class,          ONLY: t_jsb_grid
    USE mo_jsb_grid,                ONLY: Get_grid
    USE mo_jsb_physical_constants,  ONLY: tmelt

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER :: model
    INTEGER                    :: iblk, ics, ice

    INTEGER :: day_of_year

    REAL(wp), DIMENSION(options%nc) :: &
      & t_air_in_Celcius,              &
      & LaiMax_corrected,              &
      & veg_fract_correction_old,      &
      & old_c_state_sum_ta, dummy_flux

    TYPE(t_jsb_grid), POINTER :: grid

    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(CARBON_) ! JN 10.01.18 only required with l_forestRegrowth = true

    ! Declare pointers for variables in memory
    !INTEGER                :: ncanopy     ! number of canopy layers

    !INTEGER, POINTER                 :: acc_counter ! counter for each delta_time step since the last variable output
    !LOGICAL, POINTER                 :: use_alb_veg_simple
    !LOGICAL, POINTER                 :: first_call_phenology_this_day

    dsl4jsb_Real2D_onChunk :: t_air
    dsl4jsb_Real2D_onChunk :: NPP_pot_rate_ca
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_rel
    dsl4jsb_Real2D_onChunk :: lat
    dsl4jsb_Real2D_onChunk :: heat_sum_SG
    dsl4jsb_Real2D_onChunk :: heat_sum_EG
    dsl4jsb_Real2D_onChunk :: heat_sum_CRP
    dsl4jsb_Real2D_onChunk :: heat_sum_winter
    dsl4jsb_Real2D_onChunk :: days_since_growth_begin_SG
    dsl4jsb_Real2D_onChunk :: days_since_growth_begin_EG
    dsl4jsb_Real2D_onChunk :: growth_phase_SG
    dsl4jsb_Real2D_onChunk :: growth_phase_EG
    dsl4jsb_Real2D_onChunk :: growth_phase_CRP
    dsl4jsb_Real2D_onChunk :: previous_day_temp_mean
    dsl4jsb_Real2D_onChunk :: pseudo_soil_temp
    dsl4jsb_Real2D_onChunk :: chill_days_SG
    dsl4jsb_Real2D_onChunk :: chill_days_EG
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: lai_ta
    dsl4jsb_Real2D_onChunk :: day_NPP_sum
    dsl4jsb_Real2D_onChunk :: previous_day_temp_min
    dsl4jsb_Real2D_onChunk :: previous_day_temp_max
    dsl4jsb_Real2D_onChunk :: previous_day_NPP_pot_rate_ca
    dsl4jsb_Real2D_onChunk :: veg_fract_correction
    dsl4jsb_Real2D_onChunk :: fract_fpc_max

    ! JN 10.01.18 only required with l_forestRegrowth = true
    dsl4jsb_Real2D_onChunk :: c_green
    dsl4jsb_Real2D_onChunk :: c_reserve
    dsl4jsb_Real2D_onChunk :: c_woods

    dsl4jsb_Real2D_onChunk :: current_max_green
    dsl4jsb_Real2D_onChunk :: veg_carbon_at_max_green

    dsl4jsb_Real2D_onChunk :: number_of_individuals
    dsl4jsb_Real2D_onChunk :: biomass_per_individual
    dsl4jsb_Real2D_onChunk :: maxLAI_allom

    dsl4jsb_Real2D_onChunk :: cconservation_allom

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

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(PHENO_)) RETURN

    !$claw acc data create(t_air_in_Celcius)  &
    !$claw create(LaiMax_corrected, veg_fract_correction_old) &
    !$claw create(old_c_state_sum_ta, dummy_flux)
    IF ( .NOT. tile%Has_process_memory(ASSIMI_)) THEN
      CALL finish(routine, 'The LOGROP phenology model depends on the assimilation process being activated, too')
    ENDIF

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Set model
    model => Get_model(tile%owner_model_id)

    ! Set pointers to process configs and memory
    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(ASSIMI_)

    !ncanopy =  dsl4jsb_Config(ASSIMI_)%ncanopy

    grid    => Get_grid(model%grid_id)
    lat     => grid%lat(ics:ice, iblk)

    !first_call_phenology_this_day  => model%processes(PHENO_)%config%first_call_phenology_this_day(ics:ice,iblk)
    IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN

      !ASSERTION: l_forestRegrowth only possible when the CARBON process is active
      IF ( .NOT. tile%Has_process_memory(CARBON_)) THEN
        CALL finish(TRIM(routine), 'Forest regrowth can only be switched on when the CARBON process is active, too')
      ENDIF

      dsl4jsb_Get_memory(CARBON_)
      dsl4jsb_Get_var2D_onChunk(CARBON_,   c_green)             ! in
      dsl4jsb_Get_var2D_onChunk(CARBON_,   c_reserve)           ! in
      dsl4jsb_Get_var2D_onChunk(CARBON_,   c_woods)             ! in

      dsl4jsb_Get_var2D_onChunk(CARBON_,    current_max_green)       ! inout
      dsl4jsb_Get_var2D_onChunk(CARBON_,    veg_carbon_at_max_green) ! inout

      dsl4jsb_Get_var2D_onChunk(PHENO_,    number_of_individuals)   ! inout
      dsl4jsb_Get_var2D_onChunk(PHENO_,    biomass_per_individual)  ! inout
      dsl4jsb_Get_var2D_onChunk(PHENO_,    maxLAI_allom)            ! inout

      dsl4jsb_Get_var2D_onChunk(PHENO_,    cconservation_allom)     ! inout

    END IF

    ! Set pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(A2L_,      t_air)                  ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   NPP_pot_rate_ca)        ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_rel)       ! in

    dsl4jsb_Get_var2D_onChunk(SEB_,      previous_day_temp_mean) ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      previous_day_temp_min)  ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      previous_day_temp_max)  ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      pseudo_soil_temp)       ! in

    dsl4jsb_Get_var2D_onChunk(PHENO_,    fract_fpc_max)          ! in

    dsl4jsb_Get_var2D_onChunk(PHENO_,    lai)                    ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_,    veg_fract_correction)   ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_,    heat_sum_SG)            ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_,    heat_sum_EG)            ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_,    heat_sum_CRP)           ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_,    heat_sum_winter)        ! inout

    dsl4jsb_Get_var2D_onChunk(PHENO_,    days_since_growth_begin_SG) ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_,    days_since_growth_begin_EG) ! inout


    dsl4jsb_Get_var2D_onChunk(PHENO_,    growth_phase_SG)        ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_,    growth_phase_EG)        ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_,    growth_phase_CRP)       ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_,    chill_days_SG)          ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_,    chill_days_EG)          ! inout

    dsl4jsb_Get_var2D_onChunk(PHENO_,    day_NPP_sum)            ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_,    previous_day_NPP_pot_rate_ca) ! out
    dsl4jsb_Get_var2D_onChunk(PHENO_,    lai_ta)                 ! out

    ! ---------------------------
    ! Go

    ! R: NUR ZUR SICHERHEIT:
    !    Diese Abfrage sollte unnoetig sein, da der Prozess Phaenologie nur fuer veg tiles aufgerufen werden sollte.
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
    ! R: As in JS3 all soil layers are set to the same relative_moisture(GP,tile) and only the average
    !    over all layers or the uppermost layer is used, we do not need soil layers in the phenology.

    ! R: Note, in JSBACH3 the timestep was 1 day. Except of the former update_growth_phase which is now in the calc_XXX_phenology
    !    procedures the phenology including get_letItDie and get_letItGrow can be called several times a day. Within the
    !    calc_XXX_phenology procedures the is_newday-variable decides when to do the update_growth_phase stuff.

    ! R: While the time issue is independent of the model or block, some calculations should not be done for each tile but only
    !    once for the whole hsm (tile%name == 'box'):
    ! R: Die Berechnungen von update_previous_day_variables in JSBACH3 muss man bis auf NPP nicht fuer jedes Tile machen. Daher
    !    habe ich die Berechnungen von NPP von den anderen Berechnungen in update_previous_day_variables fuer JSBACH4 getrennt.
    !    Die Berechnungen die nicht fuer jedes Tile gemacht werden muessen sind in calc_previous_day_variables und die NPP-
    !    Berechnungen werden im Folgenden vorgenommen. Ausserdem ist calc_previous_day_variables nun dem SEB Prozess zugeordnet.
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
    IF ((dsl4jsb_Config(PHENO_)%l_forestRegrowth) .AND. (dsl4jsb_Lctlib_param(ForestFlag))) THEN
      ! to represent forest regrowth LAI needs to be dependant on the available leaf carbon
      ! to keep the overall structure of pheno we mimic this dependency by deriving the maxLAI
      ! from the available biomass via allometric relationships
      ! The maxLai gets updated at the beginning of each year with the biomass at the maximum green pool
      IF(is_newyear(options%current_datetime, options%dtime)) THEN
        !$claw expand parallel update savepoint=ser_jsb_PHENO_update_phenology_logrop_expand_1 tile_name=tile%name
        veg_fract_correction_old(:) = veg_fract_correction(:)

        !$claw sca forward update savepoint=ser_jsb_PHENO_cderive_maxLai_according_to_allometric_relationships_1 tile_name=tile%name
        CALL derive_maxLai_according_to_allometric_relationships( &
          & dsl4jsb_Lctlib_param(alpha_nr_ind),                   & ! in
          & dsl4jsb_Lctlib_param(beta_nr_ind),                    & ! in
          & dsl4jsb_Lctlib_param(alpha_leaf),                     & ! in
          & dsl4jsb_Lctlib_param(beta_leaf),                      & ! in
          & dsl4jsb_Lctlib_param(clumpinessFactor),               & ! in
          & dsl4jsb_Lctlib_param(specificLeafArea_C),             & ! in
          & c_woods(:),                                           & ! in
          & current_max_green(:),                                 & ! InOut
          & veg_carbon_at_max_green(:),                           & ! InOut
          & number_of_individuals(:),                             & ! InOut
          & biomass_per_individual(:),                            & ! InOut
          & veg_fract_correction(:),                              & ! InOut
          & maxLAI_allom(:)                                       & ! InOut
          & )

        CALL recalc_carbon_per_tile_vars(tile, options)
        CALL calculate_current_c_ta_state_sum(tile, options, old_c_state_sum_ta(:))

        ! redistribute
        CALL rescale_carbon_upon_reference_area_change(tile, options, veg_fract_correction_old, veg_fract_correction)

        !$claw expand parallel update
        dummy_flux(:) = 0._wp
        ! todo: adapt for GPU
        CALL check_carbon_conservation(tile, options, old_c_state_sum_ta(:), &
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
    !    the maximum lai is at least pheno_param_jsbach%all%LAI_negligible if lctlib%LaiMax < pheno_param_jsbach%all%LAI_negligible
    !    However, as lctlib%LaiMax will always be higher than pheno_param_jsbach%all%LAI_negligible we could skip this part. But to remind
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

    IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
       !$claw expand parallel update savepoint=ser_jsb_pheno_update_phenology_logrop_expand_2 tile_name=tile%name
       LaiMax_corrected(:) = dsl4jsb_var2D_onChunk(PHENO_, maxLAI_allom )
    ELSE
     ! R: Rename MaxLAI in lctlib to LAI_max in JS3 comes from theLand%Vegetation%LAI_max,
     ! which is for each vegetated tile taken from MaxLAI in the lct library.
       MaxLAI = dsl4jsb_Lctlib_param(MaxLAI)
       !$claw expand parallel update savepoint=ser_jsb_pheno_update_phenology_logrop_expand_3 tile_name=tile%name
       LaiMax_corrected(:) = MaxLAI
    END IF

    ! R: wird jedes mal (d.h. fuer jedes Tile) aufgerufen, wobei automatisch der Tile spezifische Wert aus dem lctlib genommen wird:
    ! Secure minimum Max-LAI

    LAI_negligible_tmp = pheno_param_jsbach%all%LAI_negligible
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
    IF (dsl4jsb_Lctlib_param(PhenologyType) ==  2) THEN

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
    IF (dsl4jsb_Lctlib_param(PhenologyType) ==  1) THEN
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
    IF (dsl4jsb_Lctlib_param(PhenologyType) ==  5) THEN
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
        & wtr_rootzone_rel(:),             & ! in
        & previous_day_temp_min(:),        & ! in
        & previous_day_temp_max(:),        & ! in
        & dsl4jsb_Lctlib_param(specificLeafArea_C), & ! in
        & heat_sum_CRP(:),                 & ! inout
        & heat_sum_winter,                 & ! inout
        & growth_phase_CRP(:),             & ! inout
        & lai(:)                           & ! inout
        & )
    END IF

    ! none: 0; evergreen: 1; summergreen: 2; raingreen: 3; grasses: 4; crops: 5
    IF (dsl4jsb_Lctlib_param(PhenologyType) ==  3) THEN

      IF (debug_on() .AND. iblk == 1) CALL message('     ','... raingreen ...')
      !$claw sca forward update savepoint=ser_jsb_pheno_calc_raingreen_phenology_1 tile_name=tile%name
      CALL calc_raingreen_phenology(        &
        & timestep,                         & ! in
        & wtr_rootzone_rel(:),              & ! in
        & LaiMax_corrected(:),              & ! in
        & previous_day_NPP_pot_rate_ca(:),  & ! in
        & lai(:)                            & ! inout
        & )
    END IF

    ! none: 0; evergreen: 1; summergreen: 2; raingreen: 3; grasses: 4; crops: 5
    IF (dsl4jsb_Lctlib_param(PhenologyType) ==  4) THEN
      IF (debug_on() .AND. iblk == 1) CALL message('     ','... grasses ...')
       !$claw sca forward update savepoint=ser_jsb_pheno_calc_grass_phenology_1 tile_name=tile%name
       CALL calc_grass_phenology(           &
         & t_air_in_Celcius(:),             & ! in
         & timestep,                        & ! in
         & wtr_rootzone_rel(:),             & ! in
         & LaiMax_corrected(:),             & ! in
         & previous_day_NPP_pot_rate_ca(:), & ! in
         & lai(:)                           & ! inout
         & )
    END IF

    ! Calculated LAI relative to the tile area
    !$claw expand parallel update savepoint=ser_jsb_pheno_update_phenology_logrop_expand_4 tile_name=tile%name
    lai_ta(:) = lai(:) * veg_fract_correction(:) * fract_fpc_max(:)


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

    dsl4jsb_Def_memory(PHENO_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: lai
    REAL(wp), POINTER      :: lai_mon(:,:)

    CHARACTER(len=*), PARAMETER            :: routine = modname//':update_phenology_climatology'

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(PHENO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')


    dsl4jsb_Get_memory(PHENO_)

    dsl4jsb_Get_var2D_onChunk(PHENO_, lai)                              ! out
    lai_mon(1:,0:) => dsl4jsb_memory(PHENO_)%lai_mon (ics:ice, iblk, :) ! in

    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
    !$claw expand parallel update savepoint=ser_jsb_pheno_update_phenology_climatology_expand_1 tile_name=tile%name
    lai(:) = REAL(wgt1,wp) * lai_mon(:,nmw1) + REAL(wgt2,wp) * lai_mon(:,nmw2)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_phenology_climatology

  ! ======================================================================================================= !
  !>Implementation of "aggregate" of the "phenology" task
  !> 
  SUBROUTINE aggregate_phenology(tile, options)
    USE mo_jsb_model_class,   ONLY: t_jsb_model, MODEL_JSBACH, MODEL_QUINCY
    USE mo_jsb_class,         ONLY: Get_model
    USE mo_jsb_time,          ONLY: is_newyear
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_phenology'
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
    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)
    ! dsl4jsb_Get_memory(VEG_) <- get it in the QUINCY-specific code-part (as VEG_ is not available in JSBACH)

    ! call aggregate for different variables depending on the model
    SELECT CASE (model%config%model_scheme)
    ! ----------------------------------------------------------------------------------------------------- !
    !>JSBACH
    !>
    CASE (MODEL_JSBACH)
      dsl4jsb_Aggregate_onChunk(PHENO_, lai,           weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(PHENO_, lai_ta,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(PHENO_, fract_fpc_max, weighted_by_fract)

      IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
        IF (is_newyear(options%current_datetime, options%dtime)) THEN
          dsl4jsb_Aggregate_onChunk(PHENO_, maxLAI_allom,    weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(PHENO_, cconservation_allom,    weighted_by_fract)

          !NOTE: despite changes, no aggregation of cpools necessary:
          ! with l_forestRegrowth PHENO changes Cpools only once a year at the beginning of a timestep,
          ! and other processes (NPP alloc of carbon) aggregate them
        ENDIF
      ENDIF

      !X Implementation: Start your process scheme here...
      !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob u. wie die einzeln. pheno Tasks aggregiert werden sollen.

    ! ----------------------------------------------------------------------------------------------------- !
    !>QUINCY
    !>
    CASE (MODEL_QUINCY)
      ! get VEG_ memory
      dsl4jsb_Get_memory(VEG_)
      ! PHENO_
      dsl4jsb_Aggregate_onChunk(PHENO_, growing_season              ,    weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(PHENO_, gdd                         ,    weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(PHENO_, nd_dormance                 ,    weighted_by_fract)
      ! VEG_
      dsl4jsb_Aggregate_onChunk(VEG_,   growth_req_n_tlabile_mavg   ,    weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_,   growth_req_p_tlabile_mavg   ,    weighted_by_fract)
    END SELECT

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_phenology

  ! ======================================================================================================= !
  !>update foliage projected cover
  !>
  SUBROUTINE update_fpc(tile, options)

    USE mo_pheno_process, ONLY: get_foliage_projected_cover
    USE mo_jsb_time,      ONLY: get_time_interpolation_weights

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    !----------------------------
    ! Local variables

    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: fract_fpc
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    REAL(wp), POINTER      :: fract_fpc_mon(:,:)


    INTEGER  :: iblk, ics, ice
    REAL(wp) :: wgt1, wgt2, ClumpinessFactor
    INTEGER  :: nmw1,nmw2

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_fpc'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (.NOT. tile%Is_process_calculated(PHENO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)

    dsl4jsb_Get_var2D_onChunk(PHENO_, lai)            ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_, fract_fpc_max)  ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_, fract_fpc)      ! out

    ! Update foliage projected cover
    !
    IF (TRIM(dsl4jsb_Config(PHENO_)%scheme_jsbach) == 'climatology' &
         & .OR. tile%lcts(1)%lib_id == 0) THEN         ! General vegetation tile, no specific PFT
      IF (model%config%l_compat401) THEN
        !$claw expand parallel update savepoint=ser_jsb_pheno_update_fpc_expand_1 tile_name=tile%name
        fract_fpc(:) = get_foliage_projected_cover(fract_fpc_max(:), lai(:), 1._wp)
      ELSE
        fract_fpc_mon(1:,0:) => dsl4jsb_memory(PHENO_)%fract_fpc_mon (ics:ice, iblk, :) ! in
        CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
        !$claw expand parallel update savepoint=ser_jsb_pheno_update_fpc_expand_2 tile_name=tile%name
        fract_fpc(:) = REAL(wgt1,wp) * fract_fpc_mon(:,nmw1) + REAL(wgt2,wp) * fract_fpc_mon(:,nmw2)
      END IF
    ELSE        
      ! We have PFTs and %scheme /= climatology
      ClumpinessFactor = dsl4jsb_Lctlib_param(ClumpinessFactor) ! Do a copy to avoid deep copy on accelerator
      !$claw expand parallel update savepoint=ser_jsb_pheno_update_fpc_expand_3 tile_name=tile%name
      fract_fpc(:) = get_foliage_projected_cover(fract_fpc_max(:), lai(:), ClumpinessFactor)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_fpc

  SUBROUTINE aggregate_fpc(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(PHENO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_fpc'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(PHENO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(PHENO_, fract_fpc, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_fpc


  ! ====================================================================================================== !
  !
  !> Handles merging for those pheno variables which need to be merged upon ageing (FAGE process)
  !
  SUBROUTINE merge_pheno_variables_upon_ageing(target, source, options, moved_area)

    USE mo_fage_process,     ONLY : weighted_avg_var_upon_area_movement, &
      &                             weighted_avg_per_canopy_var_upon_area_movement
    USE mo_carbon_interface, ONLY : rescale_carbon_vars_upon_ageing_induced_ref_area_change
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract),      INTENT(inout) :: target     !< target age class
    CLASS(t_jsb_tile_abstract),      INTENT(inout) :: source     !< source age class
    TYPE(t_jsb_task_options),        INTENT(in)    :: options    !< Additional run-time parameters
    REAL(wp), DIMENSION(options%nc), INTENT(in)    :: moved_area !< area moved from source to target ac
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER  :: iblk, ics, ice

    dsl4jsb_Def_memory_tile(PHENO_, target)
    dsl4jsb_Def_memory_tile(PHENO_, source)

    dsl4jsb_Real2D_onChunk :: veg_fract_correction_source
    dsl4jsb_Real2D_onChunk :: veg_fract_correction_target
    dsl4jsb_Real2D_onChunk :: lai_source
    dsl4jsb_Real2D_onChunk :: lai_target
    dsl4jsb_Real2D_onChunk :: day_npp_sum_source
    dsl4jsb_Real2D_onChunk :: day_npp_sum_target

    REAL(wp), DIMENSION(options%nc) :: old_veg_fract_correction_target, target_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':merge_pheno_variables_upon_ageing'
    ! -------------------------------------------------------------------------------------------------- !
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For target '//TRIM(target%name)//' ...')

    CALL target%Get_fraction(ics, ice, iblk, fract=target_fract(:))

    dsl4jsb_Get_memory_tile(PHENO_, source)
    dsl4jsb_Get_memory_tile(PHENO_, target)

    dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, veg_fract_correction, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, veg_fract_correction, target)

    dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, lai, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, lai, target)
    dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, day_npp_sum, source)
    dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, day_npp_sum, target)

    ! There are pheno variables that need to be merged upon ageing 
    ! because they are used for bookkeeping by various processes but are not conserved
    ! (i.e. they are not automatically dealed with by pplcc and the lcc framework)
    ! For the pheno process these are lai and the day_npp_sum
    ! [According to the specifications in Reick et al. (2021) the LAI is to be interpreted as leaf area 
    ! "per canopy area", where canopy area is defined as area covered by a hypothetical tile without canopy gaps.]
    CALL weighted_avg_per_canopy_var_upon_area_movement(target_fract(:), moved_area(:), &
      &  veg_fract_correction_source(:), veg_fract_correction_target(:), lai_source(:), lai_target(:))
    CALL weighted_avg_per_canopy_var_upon_area_movement(target_fract(:), moved_area(:), &
      &  veg_fract_correction_source(:), veg_fract_correction_target(:), day_npp_sum_source(:), day_npp_sum_target(:))

    ! Upon ageing also the veg_fract_correction which accounts for gaps in the canopy in JSBACH needs to be merged
    old_veg_fract_correction_target(:) = veg_fract_correction_target(:)
    CALL weighted_avg_var_upon_area_movement( target_fract(:), moved_area(:), &
      &                                       veg_fract_correction_source(:), veg_fract_correction_target(:))

    ! This change in reference area, however, induces changes in those per canopy variables that need
    ! to be merged upon ageing because they are used for bookkeeping by various processes
    ! but are not conserved (i.e. not automatically dealed with by pplcc and the lcc framework)
    ! As stated above, for the pheno process these are lai and the day_npp_sum
    lai_target(:)          = lai_target(:)          &
      &                      * (old_veg_fract_correction_target(:) / veg_fract_correction_target(:))
    day_npp_sum_target(:)  = day_npp_sum_target(:)  &
      &                      * (old_veg_fract_correction_target(:) / veg_fract_correction_target(:))

    ! And for carbon these are also several
    CALL rescale_carbon_vars_upon_ageing_induced_ref_area_change(target, options, &
      &                                                        old_veg_fract_correction_target(:))

  END SUBROUTINE merge_pheno_variables_upon_ageing

  !-----------------------------------------------------------------------------------------------------
  !> Calculations of diagnostic global land mean phenology output
  !!
  !! The routine is called from jsbach_finish_timestep, after the loop over the nproma blocks.
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE global_phenology_diagnostics(tile)

#ifndef __ICON__
    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_phenology_diagnostics'
    IF (debug_on()) CALL message(TRIM(routine), 'Global diagnostics only available with ICON')
#else

    USE mo_sync,                  ONLY: global_sum_array
    USE mo_jsb_grid,              ONLY: Get_grid
    USE mo_jsb_grid_class,        ONLY: t_jsb_grid

    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile

    ! Local variables
    !
    dsl4jsb_Def_memory(PHENO_)

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_phenology_diagnostics'

    ! Pointers to variables in memory

    dsl4jsb_Real2D_onDomain :: lai_ta
    dsl4jsb_Real2D_onDomain :: fract_fpc

    REAL(wp), POINTER       :: lai_ta_gmean(:)
    REAL(wp), POINTER       :: fract_fpc_gmean(:)

    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid

    REAL(wp), POINTER      :: area(:,:)
    REAL(wp), POINTER      :: notsea(:,:)
    LOGICAL,  POINTER      :: is_in_domain(:,:) ! T: cell in domain (not halo)
    REAL(wp), ALLOCATABLE  :: in_domain (:,:)   ! 1: cell in domain, 0: halo cell
    REAL(wp), ALLOCATABLE  :: scaling (:,:)
    REAL(wp)               :: global_land_area


    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_var2D_onDomain(PHENO_,  lai_ta)             ! in
    dsl4jsb_Get_var2D_onDomain(PHENO_,  fract_fpc)          ! in

    lai_ta_gmean        => PHENO__mem%lai_ta_gmean%ptr(:)   ! out
    fract_fpc_gmean  => PHENO__mem%fract_fpc_gmean%ptr(:)   ! out

    model => Get_model(tile%owner_model_id)
    grid  => Get_grid(model%grid_id)
    area         => grid%area(:,:)
    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)
    notsea       => tile%fract(:,:)   ! fraction of the box tile: notsea


    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')


    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    ! Domain Mask - to mask all halo cells for global sums (otherwise these
    ! cells are counted twice)
    ALLOCATE (in_domain(grid%nproma,grid%nblks))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSEWHERE
      in_domain = 0._wp
    END WHERE

    ALLOCATE (scaling(grid%nproma,grid%nblks))

    ! Calculate global land means phenology variables, if requested for output
    global_land_area = global_sum_array(area(:,:) * notsea(:,:) * in_domain(:,:))
    scaling(:,:) = notsea(:,:) * area(:,:) * in_domain(:,:)
    IF (PHENO__mem%lai_ta_gmean%is_in_output)     &
      &  lai_ta_gmean     = global_sum_array(lai_ta(:,:)    * scaling(:,:)) / global_land_area
    IF (PHENO__mem%fract_fpc_gmean%is_in_output)  &
      &  fract_fpc_gmean  = global_sum_array(fract_fpc(:,:) * scaling(:,:)) / global_land_area

    DEALLOCATE (scaling, in_domain)

#endif
  END SUBROUTINE global_phenology_diagnostics

#endif
END MODULE mo_pheno_interface
