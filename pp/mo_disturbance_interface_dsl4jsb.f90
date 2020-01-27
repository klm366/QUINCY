!
! This file was automatically generated from JSBACH source file:
!    src/disturbance/mo_disturbance_interface.f90
! with md5 hash: 8b15bd12e750c80c49bab44f18c70572
!
!>
!! @brief Contains the interfaces to the disturb process
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-05-04)
!! Julia Nabel (2018):
!!      ++ c conservation test
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_disturb_interface
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

  ! Use of processes in this module (Get_disturb_memory and Get_disturb_config)
  USE mo_jsb_process_class, ONLY: disturb_, fuel_, assimi_, a2l_, carbon_, seb_

  ! Use of process configurations (t_disturb_config)
  USE mo_disturb_config_class, ONLY: t_disturb_config
  USE mo_carbon_config_class, ONLY: t_carbon_config

  ! Use of process memories (t_disturb_memory)
  USE mo_a2l_memory_class, ONLY: t_a2l_memory
  USE mo_assimi_memory_class, ONLY: t_assimi_memory
  USE mo_disturb_memory_class, ONLY: t_disturb_memory
  USE mo_fuel_memory_class, ONLY: t_fuel_memory
  USE mo_carbon_memory_class, ONLY: t_carbon_memory
  USE mo_seb_memory_class, ONLY: t_seb_memory

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  ! R: Disturbance sollte nicht zu Dynveg, da sie auch stattfinden soll, wenn dynveg aus!
  !    Disturbance sollte nicht zu CARBON, da er nur zum Teil mit C zu tun hat.
  !    So verändert er über die disturbed fractions wenn dynveg an ist evtl. die
  !    cover fractions der pfts. Disturbance ist eigentlich ein Vegetations-Prozess neben Dynveg.
  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_disturb_tasks !,t_disturb_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_interface'

  !> Type definition for natural_disturbances task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_natural_disturbances
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_natural_disturbances    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_natural_disturbances !< Aggregates computed task variables
  END TYPE tsk_natural_disturbances

  !> Constructor interface for natural_disturbances task
  INTERFACE tsk_natural_disturbances
    PROCEDURE Create_task_natural_disturbances     !< Constructor function for task
  END INTERFACE tsk_natural_disturbances

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for C_natural_disturbance task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "natural_disturbances"
  !!
  FUNCTION Create_task_natural_disturbances(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_natural_disturbances::return_ptr)
    CALL return_ptr%Construct(name='natural_disturbances', process_id=DISTURB_, owner_model_id=model_id)

  END FUNCTION Create_task_natural_disturbances

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for disturb process
  !!
  !! @param[in,out] this      Instance of disturb process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_disturb_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_natural_disturbances(model_id))

  END SUBROUTINE Register_disturb_tasks


  ! ================================================================================================================================
  !>
  !! @brief Implementation of task natural_disturbances
  !! Task "natural_disturbances" calculates the disturbance of plants
  !! due to fire and windbreak.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_natural_disturbances(tile, options)

    USE mo_disturb_process,       ONLY: burned_fract_jsbach, broken_woody_fract_jsbach, get_relative_humidity_air
    USE mo_carbon_process,        ONLY: relocate_carbon_fire, relocate_carbon_damage
    USE mo_carbon_interface,      ONLY: recalc_per_tile_vars, calculate_current_c_ag_1_and_bg_sums, &
                                      & calculate_current_c_ta_state_sum, check_carbon_conservation
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_time,              ONLY: is_time_experiment_start, timesteps_per_day, is_newday
    USE mo_disturbance_constants, ONLY: persist_rel_hum, persist_wind_10m

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER                :: model
    REAL(wp), DIMENSION(options%nc)           :: old_c_state_sum, dummy_flux

    INTEGER :: iblk, nc, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_natural_disturbances'

    ! Declare process memories
    TYPE(t_a2l_memory), POINTER :: a2l__mem
    TYPE(t_assimi_memory), POINTER :: assimi__mem
    TYPE(t_carbon_memory), POINTER :: carbon__mem
    TYPE(t_disturb_memory), POINTER :: disturb__mem
    TYPE(t_fuel_memory), POINTER :: fuel__mem
    TYPE(t_seb_memory), POINTER :: seb__mem
    TYPE(t_disturb_config), POINTER :: disturb__conf
    TYPE(t_carbon_config), POINTER :: carbon__conf

    ! Declare pointers to variables in memory
    !INTEGER :: ncanopy                 ! number of canopy layers

    !REAL(wp), POINTER, DIMENSION(:) :: C4flag ! Photosynthetic pathway (C3: 0; C4: 1)
                                      ! R: nur falls man C4flag doch als REAL haben wollte


!     REAL(wp), POINTER, DIMENSION(:) :: co2_air_mol           ! CO2 particle mixing ratio [molCO2/molDryAir]

    REAL(wp), POINTER, DIMENSION(:) :: fuel                   ! Amount of fuel on which the fire was based.
    REAL(wp), POINTER, DIMENSION(:) :: q_air
    REAL(wp), POINTER, DIMENSION(:) :: press_srf
    REAL(wp), POINTER, DIMENSION(:) :: t_unfilt
    REAL(wp), POINTER, DIMENSION(:) :: wind_10m               ! 10m wind speed [m/s]
    REAL(wp), POINTER, DIMENSION(:) :: max_wind_10m_act
    REAL(wp), POINTER, DIMENSION(:) :: prev_day_max_wind_10m
    REAL(wp), POINTER, DIMENSION(:) :: max_wind_10m
    REAL(wp), POINTER, DIMENSION(:) :: q_rel_air_climbuf
    REAL(wp), POINTER, DIMENSION(:) :: q_rel_air_climbuf_yDay

    REAL(wp), POINTER, DIMENSION(:) :: cconservation_fire
    REAL(wp), POINTER, DIMENSION(:) :: cconservation_wind

    REAL(wp), POINTER, DIMENSION(:) :: burned_fract
    REAL(wp), POINTER, DIMENSION(:) :: damaged_fract
    !
    REAL(wp), POINTER, DIMENSION(:) :: c_green
    REAL(wp), POINTER, DIMENSION(:) :: c_woods
    REAL(wp), POINTER, DIMENSION(:) :: c_reserve
    REAL(wp), POINTER, DIMENSION(:) :: cflux_dist_greenreserve_2_soil
    REAL(wp), POINTER, DIMENSION(:) :: cflux_dist_woods_2_soil
    !
    REAL(wp), POINTER, DIMENSION(:) :: c_acid_ag1
    REAL(wp), POINTER, DIMENSION(:) :: c_water_ag1
    REAL(wp), POINTER, DIMENSION(:) :: c_ethanol_ag1
    REAL(wp), POINTER, DIMENSION(:) :: c_nonsoluble_ag1
    REAL(wp), POINTER, DIMENSION(:) :: c_acid_ag2
    REAL(wp), POINTER, DIMENSION(:) :: c_water_ag2
    REAL(wp), POINTER, DIMENSION(:) :: c_ethanol_ag2
    REAL(wp), POINTER, DIMENSION(:) :: c_nonsoluble_ag2
    ! below ground
    REAL(wp), POINTER, DIMENSION(:) :: c_acid_bg1
    REAL(wp), POINTER, DIMENSION(:) :: c_water_bg1
    REAL(wp), POINTER, DIMENSION(:) :: c_ethanol_bg1
    REAL(wp), POINTER, DIMENSION(:) :: c_nonsoluble_bg1
    REAL(wp), POINTER, DIMENSION(:) :: c_acid_bg2
    REAL(wp), POINTER, DIMENSION(:) :: c_water_bg2
    REAL(wp), POINTER, DIMENSION(:) :: c_ethanol_bg2
    REAL(wp), POINTER, DIMENSION(:) :: c_nonsoluble_bg2
    !
    REAL(wp), POINTER, DIMENSION(:) :: c_humus_1
    REAL(wp), POINTER, DIMENSION(:) :: c_humus_2
    REAL(wp), POINTER, DIMENSION(:) :: cflux_fire_all_2_atm
    !
    REAL(wp), POINTER, DIMENSION(:) :: cover_fract_pot
    !
    REAL(wp), POINTER, DIMENSION(:) :: cflux_fire_all_2_atm_ta

    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc) :: &
      & q_rel_air

    !LOGICAL  :: pft_mask ! True for pfts to which the disturbance should be applied.
    REAL(wp) :: persist ! Factor (as fraction) determining the persistance of q_rel_air_climbuf
                        ! relative to the time-step-actual q_rel_air
    REAL(wp) :: fire_minimum, fire_tau

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(DISTURB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get process configurations
    SELECT TYPE (cc => model%processes(disturb_)%p%config); TYPE IS (t_disturb_config) ; disturb__conf => cc ; END SELECT
    SELECT TYPE (cc => model%processes(carbon_)%p%config); TYPE IS (t_carbon_config) ; carbon__conf => cc ; END SELECT

    ! Get process memories
    SELECT TYPE (mm => tile%mem(disturb_)%p); TYPE IS (t_disturb_memory); disturb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(fuel_)%p); TYPE IS (t_fuel_memory); fuel__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT
    !SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(assimi_)%p); TYPE IS (t_assimi_memory); assimi__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT

    ! Get process variables
    fuel => fuel__mem%fuel%ptr(ics:ice,iblk)                     ! in
    q_air => a2l__mem%q_air%ptr(ics:ice,iblk)                    ! in
    press_srf => a2l__mem%press_srf%ptr(ics:ice,iblk)                ! in
    wind_10m => a2l__mem%wind_10m%ptr(ics:ice,iblk)                 ! in
    t_unfilt => seb__mem%t_unfilt%ptr(ics:ice,iblk)                 ! in
    max_wind_10m_act => disturb__mem%max_wind_10m_act%ptr(ics:ice,iblk)         ! inout
    prev_day_max_wind_10m => disturb__mem%prev_day_max_wind_10m%ptr(ics:ice,iblk)    ! inout
    max_wind_10m => disturb__mem%max_wind_10m%ptr(ics:ice,iblk)             ! inout
    q_rel_air_climbuf => disturb__mem%q_rel_air_climbuf%ptr(ics:ice,iblk)        ! inout
    q_rel_air_climbuf_yDay => disturb__mem%q_rel_air_climbuf_yDay%ptr(ics:ice,iblk)   ! out

    burned_fract => disturb__mem%burned_fract%ptr(ics:ice,iblk)             ! out
    damaged_fract => disturb__mem%damaged_fract%ptr(ics:ice,iblk)            ! out
    !
    cconservation_wind => disturb__mem%cconservation_wind%ptr(ics:ice,iblk)     ! out
    cconservation_fire => disturb__mem%cconservation_fire%ptr(ics:ice,iblk)     ! out
    !
    c_green => carbon__mem%c_green%ptr(ics:ice,iblk)                  ! inout
    c_woods => carbon__mem%c_woods%ptr(ics:ice,iblk)                  ! inout
    c_reserve => carbon__mem%c_reserve%ptr(ics:ice,iblk)                ! inout
    cflux_dist_greenreserve_2_soil => carbon__mem%cflux_dist_greenreserve_2_soil%ptr(ics:ice,iblk)! inout
    cflux_dist_woods_2_soil => carbon__mem%cflux_dist_woods_2_soil%ptr(ics:ice,iblk)  ! inout
    !
    c_acid_ag1 => carbon__mem%c_acid_ag1%ptr(ics:ice,iblk)              ! inout
    c_water_ag1 => carbon__mem%c_water_ag1%ptr(ics:ice,iblk)             ! inout
    c_ethanol_ag1 => carbon__mem%c_ethanol_ag1%ptr(ics:ice,iblk)           ! inout
    c_nonsoluble_ag1 => carbon__mem%c_nonsoluble_ag1%ptr(ics:ice,iblk)        ! inout
    c_acid_ag2 => carbon__mem%c_acid_ag2%ptr(ics:ice,iblk)              ! inout
    c_water_ag2 => carbon__mem%c_water_ag2%ptr(ics:ice,iblk)             ! inout
    c_ethanol_ag2 => carbon__mem%c_ethanol_ag2%ptr(ics:ice,iblk)           ! inout
    c_nonsoluble_ag2 => carbon__mem%c_nonsoluble_ag2%ptr(ics:ice,iblk)        ! inout
    ! below ground
    c_acid_bg1 => carbon__mem%c_acid_bg1%ptr(ics:ice,iblk)              ! inout
    c_water_bg1 => carbon__mem%c_water_bg1%ptr(ics:ice,iblk)             ! inout
    c_ethanol_bg1 => carbon__mem%c_ethanol_bg1%ptr(ics:ice,iblk)           ! inout
    c_nonsoluble_bg1 => carbon__mem%c_nonsoluble_bg1%ptr(ics:ice,iblk)        ! inout
    c_acid_bg2 => carbon__mem%c_acid_bg2%ptr(ics:ice,iblk)              ! inout
    c_water_bg2 => carbon__mem%c_water_bg2%ptr(ics:ice,iblk)             ! inout
    c_ethanol_bg2 => carbon__mem%c_ethanol_bg2%ptr(ics:ice,iblk)           ! inout
    c_nonsoluble_bg2 => carbon__mem%c_nonsoluble_bg2%ptr(ics:ice,iblk)        ! inout
    !
    c_humus_1 => carbon__mem%c_humus_1%ptr(ics:ice,iblk)               ! inout
    c_humus_2 => carbon__mem%c_humus_2%ptr(ics:ice,iblk)               ! inout
    cflux_fire_all_2_atm => carbon__mem%cflux_fire_all_2_atm%ptr(ics:ice,iblk)    ! inout
    cover_fract_pot => assimi__mem%cover_fract_pot%ptr(ics:ice,iblk)         ! in

    !ta variable required for the carbon conservation test
    cflux_fire_all_2_atm_ta => carbon__mem%cflux_fire_all_2_atm_ta%ptr(ics:ice,iblk) ! out

    ! ---------------------------
    ! Go

    q_rel_air(:) = get_relative_humidity_air(nc, q_air(:), t_unfilt(:), press_srf(:))

    ! For each time step update the climate buffer

    ! R: This may be later in a climbuf process for dynveg and disturbance.
    ! Maximum daily wind speed and relative air humidity smoothed in time
    persist = persist_rel_hum ** (1._wp/REAL(timesteps_per_day(options%dtime),wp))
    q_rel_air_climbuf(:) = q_rel_air_climbuf(:) * persist  +  MIN(q_rel_air(:),100._wp) * (1._wp - persist)

    IF (.NOT. is_newday(options%datetime, options%dtime) .OR. is_time_experiment_start(options%datetime)) THEN
      WHERE (wind_10m(:) > max_wind_10m_act(:))
         max_wind_10m_act(:) = wind_10m(:)
      END WHERE
    ELSE
      prev_day_max_wind_10m(:) = max_wind_10m_act(:)
      max_wind_10m(:) = max_wind_10m(:) * persist_wind_10m + max_wind_10m_act(:) * (1._wp - persist_wind_10m)
      max_wind_10m_act(:) = wind_10m(:)
      q_rel_air_climbuf_yDay(:) = q_rel_air_climbuf(:)
    END IF

    IF (is_newday(options%datetime, options%dtime)) THEN ! only once per day
      cflux_dist_greenreserve_2_soil(:) = 0._wp
      cflux_dist_woods_2_soil(:) = 0._wp

      burned_fract(:) = 0._wp
      IF (      model%lctlib(tile%lcts(1)%lib_id)%dynamic_PFT &
          .OR. (disturb__conf%lburn_pasture .AND. model%lctlib(tile%lcts(1)%lib_id)%pasture_PFT) &
         ) THEN

        SELECT CASE (disturb__conf%fire_algorithm)
          CASE (0) !! No fire algorithm
          CASE (1) !! jsbach algorithm

            ! Calculate burned area. Woody and grass types call the subroutine with different parameters
            IF (model%lctlib(tile%lcts(1)%lib_id)%woody_PFT) THEN  ! Woody Type
              fire_minimum = disturb__conf%fire_minimum_woody
              fire_tau     = disturb__conf%fire_tau_woody
            ELSE                             ! all other types
              fire_minimum = disturb__conf%fire_minimum_grass
              fire_tau     = disturb__conf%fire_tau_grass
            ENDIF

            CALL burned_fract_jsbach(                            &
              & disturb__conf%fire_rel_hum_threshold, &
              & disturb__conf%fire_litter_threshold,  &
              & fire_minimum,                                    & ! in
              & fire_tau,                                        & ! in
              & q_rel_air_climbuf(:),                            & ! in
              & fuel(:),                                         & ! in
              & burned_fract(:)                                  & ! out
              & )

          CASE (2) !! Arora & Boer algorithm
            CALL finish('disturbed_fract','Arora & Boer algorithm not implemented yet.')
          CASE (3)
            CALL finish('disturbed_fract','Thornicke spit fire not implemented yet.')
          CASE (4) !! Read burned_fract
          CASE DEFAULT
            CALL finish('disturbed_fract','Unknown fire algorithm')
        END SELECT

      END IF

      damaged_fract(:) = 0._wp
      IF (model%lctlib(tile%lcts(1)%lib_id)%woody_PFT .AND. model%lctlib(tile%lcts(1)%lib_id)%dynamic_PFT) THEN

        SELECT CASE (disturb__conf%windbreak_algorithm)
          CASE (0) !! No windbreak algorithm
          CASE (1) !! jsbach algorithm
            CALL broken_woody_fract_jsbach(                 & ! in
              & disturb__conf%wnd_threshold,     & ! in
              & disturb__conf%wnd_damage_scale,  & ! in
              & cover_fract_pot(:),                         & ! in
              & prev_day_max_wind_10m(:),                   & ! in
              & max_wind_10m(:),                            & ! in
              & damaged_fract(:)                            & ! out
              & )
          CASE DEFAULT
            CALL finish('disturbed_frac','Unknown windbreak algorithm')
        END SELECT

      END IF

      IF ( tile%Is_process_active(CARBON_) ) THEN

        ! Calculate carbon relocation connected with damages and fires
        ! R: hier sollen nur die Variablen übergeben werden, die von dieser Prozedur hier im Folgenden benötigt werden.
        !    Andere Variablen werden erst später bereit gestellt, denn in diesem Fall rufen wir ein weiteres
        !    Interface beim CARBON Prozess auf (SUBROUTINE relocate_disturbed_carbon).
        !    Erst dort sollen dann alle Variablen an die eigentlichen Prozess-Prozeduren übergeben werden.

        SELECT CASE (disturb__conf%windbreak_algorithm)
          CASE (0) !! No windbreak algorithm
          CASE (1) !! jsbach algorithm

            !Calculate current sum before operation for c conservation test
            CALL calculate_current_c_ta_state_sum(tile, options, old_c_state_sum(:))

            CALL relocate_carbon_damage(                 &
              & damaged_fract(:),                        & ! in
              & model%lctlib(tile%lcts(1)%lib_id)%LeafLit_coef(1:5), & ! in
              & model%lctlib(tile%lcts(1)%lib_id)%WoodLit_coef(1:5), & ! in
              & c_green(:),                              & ! inout
              & c_woods(:),                              & ! inout
              & c_reserve(:),                            & ! inout
              & cflux_dist_greenreserve_2_soil(:),       & ! inout
              & cflux_dist_woods_2_soil(:),              & ! inout
              & c_acid_ag1(:),                           & ! inout
              & c_water_ag1(:),                          & ! inout
              & c_ethanol_ag1(:),                        & ! inout
              & c_nonsoluble_ag1(:),                     & ! inout
              & c_acid_ag2(:),                           & ! inout
              & c_water_ag2(:),                          & ! inout
              & c_ethanol_ag2(:),                        & ! inout
              & c_nonsoluble_ag2(:),                     & ! inout
              & c_acid_bg1(:),                           & ! inout
              & c_water_bg1(:),                          & ! inout
              & c_ethanol_bg1(:),                        & ! inout
              & c_nonsoluble_bg1(:),                     & ! inout
              & c_acid_bg2(:),                           & ! inout
              & c_water_bg2(:),                          & ! inout
              & c_ethanol_bg2(:),                        & ! inout
              & c_nonsoluble_bg2(:),                     & ! inout
              & c_humus_1(:),                            & ! inout
              & c_humus_2(:)                             & ! inout
              & )

            CALL calculate_current_c_ag_1_and_bg_sums(tile, options)
            CALL recalc_per_tile_vars(tile, options,                                      &
              & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:),      &
              & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:),                 &
              & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
              & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:),                 &
              & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
              & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:),                 &
              & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
              & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:),                 &
              & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
              & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:),                       &
              & cflux_dist_greenreserve_2_soil = cflux_dist_greenreserve_2_soil(:),       &
              & cflux_dist_woods_2_soil = cflux_dist_woods_2_soil(:))

            ! wind has no extra flux to atmos, only redistribution among cpools - therefore cflux = 0._wp
            dummy_flux = 0._wp
            CALL check_carbon_conservation(tile, options, old_c_state_sum(:), &
              & dummy_flux(:), cconservation_wind(:))

          CASE DEFAULT
            CALL finish('disturbed_fract','Unknown windbreak algorithm')
        END SELECT

        cflux_fire_all_2_atm = 0._wp
        SELECT CASE (disturb__conf%fire_algorithm)
          CASE (0) !! No fire algorithm
          CASE (1) !! jsbach algorithm

           !Calculate current sum before operation for c conservation test
            CALL calculate_current_c_ta_state_sum(tile, options, old_c_state_sum(:))

            CALL relocate_carbon_fire(                    &
              & burned_fract(:),                          & ! in, only for windbreak
              & carbon__conf%fire_fract_wood_2_atmos , & ! in
              & model%lctlib(tile%lcts(1)%lib_id)%LeafLit_coef(1:5),  & ! in
              & model%lctlib(tile%lcts(1)%lib_id)%WoodLit_coef(1:5),  & ! in
              & c_green(:),                               & ! inout
              & c_woods(:),                               & ! inout
              & c_reserve(:),                             & ! inout
              & cflux_dist_greenreserve_2_soil(:),        & ! inout
              & cflux_dist_woods_2_soil(:),               & ! inout
              & c_acid_ag1(:),                            & ! inout
              & c_water_ag1(:),                           & ! inout
              & c_ethanol_ag1(:),                         & ! inout
              & c_nonsoluble_ag1(:),                      & ! inout
              & c_acid_ag2(:),                            & ! inout
              & c_water_ag2(:),                           & ! inout
              & c_ethanol_ag2(:),                         & ! inout
              & c_nonsoluble_ag2(:),                      & ! inout
              & c_acid_bg1(:),                            & ! inout
              & c_water_bg1(:),                           & ! inout
              & c_ethanol_bg1(:),                         & ! inout
              & c_nonsoluble_bg1(:),                      & ! inout
              & c_acid_bg2(:),                            & ! inout
              & c_water_bg2(:),                           & ! inout
              & c_ethanol_bg2(:),                         & ! inout
              & c_nonsoluble_bg2(:),                      & ! inout
              & c_humus_1(:),                             & ! inout
              & c_humus_2(:),                             & ! inout
              & cflux_fire_all_2_atm                      & ! out
              & )

            CALL calculate_current_c_ag_1_and_bg_sums(tile, options)
            CALL recalc_per_tile_vars(tile, options,                                      &
              & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:),      &
              & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:),                 &
              & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
              & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:),                 &
              & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
              & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:),                 &
              & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
              & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:),                 &
              & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
              & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:),                       &
              & cflux_dist_greenreserve_2_soil = cflux_dist_greenreserve_2_soil(:),       &
              & cflux_dist_woods_2_soil = cflux_dist_woods_2_soil(:),                     &
              & cflux_fire_all_2_atm = cflux_fire_all_2_atm(:))

            !@todo: shouldn't cflux_fire_all_2_atm be negative since its a flux away from land?
            !For conservation test: negate currently positve cflux_fire_all_2_atm (flux away from land!)
            CALL check_carbon_conservation(tile, options, old_c_state_sum(:), &
              & -cflux_fire_all_2_atm_ta(:), cconservation_fire(:))

          CASE (2) !! Arora & Boer algorithm
            CALL finish('disturbed_fract','Arora & Boer algorithm not implemented yet.')
          CASE (3)
            CALL finish('disturbed_fract','Thornicke spit fire not implemented yet.')
          CASE (4) !! Read burned_fract
          CASE DEFAULT
            CALL finish('disturbed_fract','Unknown fire algorithm')
        END SELECT

        ! R: This corresponds to JS3, but I did not understand it...
        !    disturbance%disturb_2_atmos corresponds to cflux_fire_all_2_atm => carbon__mem%cflux_fire_all_2_atm%ptr(ics:ice,iblk)
        !    However, I guess it is better not to take this from JS3.
        !IF ( tile%Is_process_active(CARBON_) ) THEN
        !  ! Calculate CO2 and other fluxes to the atmosphere (updated once a day)
        !  CO2_emission(:) = disturbance%disturb_2_atmos(kidx0:kidx1)   * molarMassCO2_kg/86400._wp
        !ENDIF

      ENDIF

    ENDIF ! is_newday

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_natural_disturbances

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "natural_disturbances"
  !!
  !! @param[in,out] tile    Tile for which aggregation of child tiles is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_natural_disturbances(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_disturb_memory), POINTER :: disturb__mem
    TYPE(t_carbon_memory), POINTER :: carbon__mem

    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_natural_disturbance'

    INTEGER  :: iblk, ics, ice

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    SELECT TYPE (mm => tile%mem(disturb_)%p); TYPE IS (t_disturb_memory); disturb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%Aggregate(tile, disturb__mem%burned_fract,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, disturb__mem%damaged_fract,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, disturb__mem%q_rel_air_climbuf,ics,ice,iblk,routine)
    CALL weighted_by_fract%Aggregate(tile, disturb__mem%q_rel_air_climbuf_yDay,ics,ice,iblk,routine)

    ! R: I do not think these variables should be aggregated:
    !CALL weighted_by_fract%Aggregate(tile, disturb__mem%max_wind_10m_act,ics,ice,iblk,routine)
    !CALL weighted_by_fract%Aggregate(tile, disturb__mem%prev_day_max_wind_10m,ics,ice,iblk,routine)
    !CALL weighted_by_fract%Aggregate(tile, disturb__mem%max_wind_10m,ics,ice,iblk,routine)

    IF ( tile%Is_process_active(CARBON_) ) THEN
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_green_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_woods_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_reserve_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%cflux_dist_greenreserve_2_soil_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%cflux_dist_woods_2_soil_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_acid_ag1_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_water_ag1_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_ethanol_ag1_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_nonsoluble_ag1_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_acid_ag2_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_water_ag2_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_ethanol_ag2_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_nonsoluble_ag2_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_humus_1_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%c_humus_2_ta,ics,ice,iblk,routine)
       CALL weighted_by_fract%Aggregate(tile, carbon__mem%cflux_fire_all_2_atm_ta,ics,ice,iblk,routine)
    ENDIF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_natural_disturbances

#endif
END MODULE mo_disturb_interface
