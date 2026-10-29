!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_q_dfire_sf_init
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: lp
  USE mo_exception,               ONLY: message, finish
  USE mo_jsb_control,             ONLY: debug_on

  USE mo_jsb_model_class,         ONLY: t_jsb_model
  USE mo_jsb_grid_class,          ONLY: t_jsb_grid
  USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract

  USE mo_q_dfire_sf_constants,    ONLY: ratio_as_dead_fuel
#ifndef __QUINCY_STANDALONE__
  USE mo_jsb_time,                ONLY: timestep_in_days
#endif
  USE mo_jsb_process_class,       ONLY: Q_DFIRE_SF_, SB_, Q_DFIRE_FUEL_SF_

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: q_dfire_sf_init

#ifndef __QUINCY_STANDALONE__
  PUBLIC ::  provide_lightning_events
#endif

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_init'

  TYPE :: t_q_dfire_init_vars
    REAL(lp), POINTER ::      &
        pop_dens(:,:) => NULL(), &
        a_nd(:,:) => NULL(), &
        lightning_events(:,:,:) => NULL()
  END TYPE t_q_dfire_init_vars

  TYPE(t_q_dfire_init_vars) :: dfire_init_vars

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Intialize disturbance fire process (after memory has been set up)
  !>
  SUBROUTINE q_dfire_sf_init(tile)
    USE mo_jsb_class,               ONLY: get_model
    USE mo_jsb_grid,                ONLY: Get_grid

    USE mo_q_dfire_sf_constants,    ONLY: Q_DFIRE_PRECIP_METHOD_DAYTIME, & !< Method to compute precipitation
                                          Q_DFIRE_PRECIP_METHOD_APROX_MAVG, &
                                          Q_DFIRE_PRECIP_METHOD_MAVG

    dsl4jsb_Use_config(Q_DFIRE_SF_)
    dsl4jsb_Use_memory(Q_DFIRE_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER              :: model
    TYPE(t_jsb_grid),       POINTER              :: hgrid              !< Horizontal grid
    INTEGER                                      :: iblk, ics, ice, nc, npft
    INTEGER :: nblks, nproma
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':q_dfire_sf_init'
    ! ----------------------------------------------------------------------------------------------------- !
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(Q_DFIRE_SF_)
    dsl4jsb_Def_memory(Q_DFIRE_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: dfire_pm_id                                  !< Selected precipitation method
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real2D_onDomain      :: mask_is_fire_intense_enough
    dsl4jsb_Real2D_onDomain      :: precip_daytime
    dsl4jsb_Real2D_onDomain      :: precip_daytime_old
    dsl4jsb_Real2D_onDomain      :: precip_qc_mavg
    dsl4jsb_Real2D_onDomain      :: precip_daily_mmd
    dsl4jsb_Real2D_onDomain      :: precip_mavg

    dsl4jsb_Real2D_onDomain      :: t_air_min
    dsl4jsb_Real2D_onDomain      :: t_air_max
    dsl4jsb_Real2D_onDomain      :: t_air_min_old
    dsl4jsb_Real2D_onDomain      :: t_air_max_old

    dsl4jsb_Real2D_onDomain      :: nesterov_index
    dsl4jsb_Real2D_onDomain      :: fdi
    dsl4jsb_Real2D_onDomain      :: ros_f
    dsl4jsb_Real2D_onDomain      :: ros_b
    dsl4jsb_Real2D_onDomain      :: mean_fire_area
    dsl4jsb_Real2D_onDomain      :: fire_duration

    dsl4jsb_Real2D_onDomain      :: area_burned
    dsl4jsb_Real2D_onDomain      :: area_burned_fract
    dsl4jsb_Real2D_onDomain      :: area_burned_output
    dsl4jsb_Real2D_onDomain      :: area_burned_fract_output
    dsl4jsb_Real2D_onDomain      :: i_surface

    dsl4jsb_Real3D_onDomain      :: fuel_moisture

    dsl4jsb_Real2D_onDomain      :: e_n_ignitions
    dsl4jsb_Real2D_onDomain      :: e_n_humans
    dsl4jsb_Real2D_onDomain      :: e_n_lightning
    dsl4jsb_Real2D_onDomain      :: ab_moisture

    dsl4jsb_Real2D_onDomain      :: reaction_intensity
    dsl4jsb_Real2D_onDomain      :: prop_flux_ratio
    dsl4jsb_Real2D_onDomain      :: wind_factor
    dsl4jsb_Real2D_onDomain      :: eff_heating_number
    dsl4jsb_Real2D_onDomain      :: preig_heat

    dsl4jsb_Real2D_onDomain      :: opt_react_velo
    dsl4jsb_Real2D_onDomain      :: net_fuel_load
    dsl4jsb_Real2D_onDomain      :: moist_damp_coeff
    dsl4jsb_Real2D_onDomain      :: mineral_damp_coeff

    dsl4jsb_Real3D_onDomain      :: dead_fuel_consumption_1hr
    dsl4jsb_Real3D_onDomain      :: dead_fuel_consumption_10hr
    dsl4jsb_Real3D_onDomain      :: dead_fuel_consumption_100hr
    dsl4jsb_Real3D_onDomain      :: dead_fuel_consumption_1000hr

    dsl4jsb_Real2D_onDomain      :: distance_travelled
    dsl4jsb_Real2D_onDomain      :: l2b_ratio

    dsl4jsb_Real3D_onDomain      :: frac_consump_fcl
    dsl4jsb_Real2D_onDomain      :: ic_pop_dens
    dsl4jsb_Real3D_onDomain      :: ic_lightning_events_monthly
    dsl4jsb_Real2D_onDomain      :: a_nd
    ! ----------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_active(Q_DFIRE_SF_)) RETURN
    !IF (tile%lcts(1)%lib_id == 0) RETURN                !< run this init only if the present tile is a pft
    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model  => get_model(tile%owner_model_id)
    hgrid  => Get_grid(model%grid_id)
    nblks  =  hgrid%nblks
    nproma =  hgrid%nproma
    ! ----------------------------------------------------------------------------------------------------- !
    ! CALL model%processes(Q_DFIRE_FUEL_SF_)%p%Test_init_status_asked_by(Q_DFIRE_SF_)    !< Test if process dependency is already initialized
    ! CALL model%processes(SB_)%p%Test_init_status_asked_by(Q_DFIRE_SF_)    !< Test if process dependency is already initialized
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_DFIRE_SF_)
    dsl4jsb_Get_memory(Q_DFIRE_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, nesterov_index)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, fdi)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, ros_f)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, ros_b)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, mean_fire_area)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, fire_duration)

    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, e_n_ignitions)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, e_n_humans)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, e_n_lightning)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, ab_moisture)

    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, reaction_intensity)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, prop_flux_ratio)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, wind_factor)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, eff_heating_number)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, preig_heat)

    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, opt_react_velo)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, net_fuel_load)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, moist_damp_coeff)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, mineral_damp_coeff)

    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, mask_is_fire_intense_enough)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, precip_qc_mavg)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, precip_daily_mmd)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, precip_mavg)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, precip_daytime)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, precip_daytime_old)

    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, t_air_min)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, t_air_max)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, t_air_min_old)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, t_air_max_old)

    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, area_burned)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, area_burned_fract)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, area_burned_output)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, area_burned_fract_output)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, i_surface)

    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_SF_, fuel_moisture)

    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_SF_, dead_fuel_consumption_1hr)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_SF_, dead_fuel_consumption_10hr)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_SF_, dead_fuel_consumption_100hr)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_SF_, dead_fuel_consumption_1000hr)

    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, distance_travelled)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, l2b_ratio)

    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_SF_, frac_consump_fcl)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, ic_pop_dens)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_SF_, ic_lightning_events_monthly)
    dsl4jsb_Get_var2D_onDomain(Q_DFIRE_SF_, a_nd)
    ! ----------------------------------------------------------------------------------------------------- !
    nesterov_index = 0._lp
    fdi = 0._lp
    ros_f = 0._lp
    ros_b = 0._lp
    mean_fire_area = 0._lp
    fire_duration = 0._lp

    e_n_ignitions = 0._lp
    e_n_humans = 0._lp
    e_n_lightning = 0._lp
    ab_moisture = 0._lp

    mask_is_fire_intense_enough = 0._lp !< TRUE=1._lp, FALSE=0._lp
    precip_daytime = 0._lp
    precip_daytime_old = 0._lp
    precip_qc_mavg = 0._lp
    precip_daily_mmd = 0._lp
    precip_mavg = 0._lp

    t_air_min = 0._lp
    t_air_max = 0._lp
    t_air_min_old = 0._lp
    t_air_max_old = 0._lp

    area_burned = 0._lp
    area_burned_fract = 0._lp
    area_burned_output = 0._lp
    area_burned_fract_output = 0._lp
    i_surface = 0._lp
    dead_fuel_consumption_1hr = 0._lp
    dead_fuel_consumption_10hr = 0._lp
    dead_fuel_consumption_100hr = 0._lp
    dead_fuel_consumption_1000hr = 0._lp

    fuel_moisture = 0._lp

    reaction_intensity = 0._lp
    prop_flux_ratio = 0._lp
    wind_factor = 0._lp
    eff_heating_number = 0._lp
    preig_heat = 0._lp

    distance_travelled = 0._lp
    l2b_ratio = 0._lp
    frac_consump_fcl = 0._lp

    ic_pop_dens = 0._lp
    a_nd = 0._lp

    !< Chose daily precipitation. Those are calculated in 3 different ways.
    dfire_pm_id = dsl4jsb_Config(Q_DFIRE_SF_)%dfire_precip_method_id
    IF (dfire_pm_id == Q_DFIRE_PRECIP_METHOD_DAYTIME) THEN
      !< daytime
      dsl4jsb_var_ptr(Q_DFIRE_SF_,precip_daily) => dsl4jsb_var_ptr(Q_DFIRE_SF_,precip_daytime)
    ELSE IF (dfire_pm_id == Q_DFIRE_PRECIP_METHOD_MAVG) THEN
      !< moving average
      dsl4jsb_var_ptr(Q_DFIRE_SF_,precip_daily) => dsl4jsb_var_ptr(Q_DFIRE_SF_,precip_mavg)
    ELSE IF (dfire_pm_id == Q_DFIRE_PRECIP_METHOD_APROX_MAVG) THEN
      !< aproximation to moving average (quincy default method)
      dsl4jsb_var_ptr(Q_DFIRE_SF_,precip_daily) => dsl4jsb_var_ptr(Q_DFIRE_SF_,precip_qc_mavg)
    ENDIF

    !< Initialize ratio of sb pools to dead fuel
    ratio_as_dead_fuel % polymeric      = dsl4jsb_Config(Q_DFIRE_SF_)%ratio_polymeric_litter_to_dead_fuel
    ratio_as_dead_fuel % dom_assoc      = dsl4jsb_Config(Q_DFIRE_SF_)%ratio_dom_assoc_to_dead_fuel
    ratio_as_dead_fuel % dom            = dsl4jsb_Config(Q_DFIRE_SF_)%ratio_dom_to_dead_fuel
    ratio_as_dead_fuel % soluable       = dsl4jsb_Config(Q_DFIRE_SF_)%ratio_soluable_litter_to_dead_fuel
    ratio_as_dead_fuel % residue        = dsl4jsb_Config(Q_DFIRE_SF_)%ratio_residue_to_dead_fuel
    ratio_as_dead_fuel % residue_assoc  = dsl4jsb_Config(Q_DFIRE_SF_)%ratio_residue_assoc_to_dead_fuel

#ifdef __QUINCY_STANDALONE__
    !> Quincy standalone only
    !> Population density and lightning events inputs are both loaded in q_dfire_read_init_vars
    CALL q_dfire_read_init_vars(tile)

    ic_pop_dens = dfire_init_vars%pop_dens
    a_nd = dfire_init_vars%a_nd
#else
    !> IQ only
    !> -population density is loaded in q_dfire_read_init_vars
    !> -lightning events are loaded in provide_lightning_events
    ! call only at root tile, i.e., the only tile without associated parent tile (to avoid unnessary i/o)
    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL q_dfire_read_init_vars(tile)
    END IF
    ic_pop_dens = dfire_init_vars%pop_dens(1,1)
    a_nd = dfire_init_vars%a_nd(1,1)
#endif
    ! Mark process as successfully initialized, MR 428
    !CALL model%processes(Q_DFIRE_SF_)%p%Set_init_status_true
  END SUBROUTINE q_dfire_sf_init


#ifdef __QUINCY_STANDALONE__
  SUBROUTINE q_dfire_read_init_vars(tile)
    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: get_model
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid
    USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
    USE mo_jsb_io,             ONLY: missval

    dsl4jsb_Use_config(Q_DFIRE_SF_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(Q_DFIRE_SF_)

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid
    INTEGER                     :: nproma
    INTEGER                     :: nblks
    CHARACTER(len=*), PARAMETER :: routine = modname//':q_dfire_read_init_vars'

    IF (debug_on()) CALL message(routine, 'Reading/setting DFIRE_ init vars')

    model  => Get_model(tile%owner_model_id)
    hgrid  => Get_grid(model%grid_id)
    nproma = hgrid%Get_nproma()
    nblks  = hgrid%Get_nblks()

    ALLOCATE( &
      & dfire_init_vars%pop_dens               (nproma, nblks),       &
      & dfire_init_vars%a_nd                   (nproma, nblks),       &
      & dfire_init_vars%lightning_events       (nproma, nblks, 12),   & ! 12=months in a year
      & )

    dfire_init_vars%pop_dens              (:,:) = missval
    dfire_init_vars%a_nd                  (:,:) = missval
    dfire_init_vars%lightning_events    (:,:,:) = missval

    dsl4jsb_Get_config(Q_DFIRE_SF_)
    ! ----------------------------------------------------------------------------------------------------- !

    dfire_init_vars%pop_dens(:,:) = dsl4jsb_Config(Q_DFIRE_SF_)%density_human_population
    dfire_init_vars%lightning_events(:,:,:) = dsl4jsb_Config(Q_DFIRE_SF_)%lightning_frequency
    dfire_init_vars%a_nd(:,:) = dsl4jsb_Config(Q_DFIRE_SF_)%a_nd
  END SUBROUTINE q_dfire_read_init_vars
#else
  ! ======================================================================================================= !
  !> read DFIRE_ init vars from input file
  !>
  SUBROUTINE q_dfire_read_init_vars(tile)
    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: get_model
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid
    USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_jsb_io,             ONLY: missval

    dsl4jsb_Use_config(Q_DFIRE_SF_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(Q_DFIRE_SF_)

    REAL(lp), POINTER           :: ptr_2D(:, :)
    REAL(lp), POINTER           :: ptr_2D_time(:, :, :)
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid
    TYPE(t_input_file)          :: input_file
    INTEGER                     :: nproma
    INTEGER                     :: nblks
    CHARACTER(len=*), PARAMETER :: routine = modname//':q_dfire_read_init_vars'

    IF (debug_on()) CALL message(routine, 'Reading/setting DFIRE_ init vars')

    model  => Get_model(tile%owner_model_id)
    hgrid  => Get_grid(model%grid_id)
    nproma = hgrid%Get_nproma()
    nblks  = hgrid%Get_nblks()

    ALLOCATE( &
        & dfire_init_vars%pop_dens               (nproma, nblks),       &
        & dfire_init_vars%a_nd                   (nproma, nblks))

    dfire_init_vars%pop_dens              (:,:) = missval
    dfire_init_vars%a_nd                  (:,:) = missval

    dsl4jsb_Get_config(Q_DFIRE_SF_)
    ! ----------------------------------------------------------------------------------------------------- !

    ! Read data from populations density file
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(Q_DFIRE_SF_)%ic_popdens_filename), model%grid_id)

    !> pop dens
    ptr_2D => input_file%Read_2d(  &
        & variable_name='Band1', &
        & fill_array = dfire_init_vars%pop_dens)
    ! Clean: set negative values to 0
    dfire_init_vars%pop_dens = MERGE(ptr_2D, 0._lp, ptr_2D >= 0._lp)

    CALL input_file%Close()

    ! Read data from human ignitions file
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(Q_DFIRE_SF_)%ic_human_igns_filename), model%grid_id)

    !> a_nd
    ptr_2D => input_file%Read_2d(  &
        & variable_name='a_nd', &
        & fill_array = dfire_init_vars%a_nd)
    ! Clean: set negative values to 0
    dfire_init_vars%a_nd = MERGE(ptr_2D, 0._lp, ptr_2D >= 0._lp)

    CALL input_file%Close()

  END SUBROUTINE q_dfire_read_init_vars


  ! ====================================================================================================== !
  !
  !> Provide lightning events depending on the month
  !> Call at the start of the month to update lightning events data (monthly forcing data)
  !
  SUBROUTINE provide_lightning_events(model_id, current_datetime, dtime)

    USE mo_jsb_time_iface,         ONLY: t_datetime
    USE mo_jsb_time,               ONLY: is_newyear, get_year, get_month, is_time_experiment_start
    USE mo_io_units,               ONLY: filename_max
    USE mo_jsb_io_netcdf,          ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_jsb_grid,               ONLY: Get_grid
    USE mo_jsb_class,              ONLY: get_model
    USE mo_jsb_parallel,           ONLY: Get_omp_thread

    dsl4jsb_Use_memory(Q_DFIRE_SF_)
    dsl4jsb_Use_config(Q_DFIRE_SF_)
    dsl4jsb_Use_config(SB_)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id
    TYPE(t_datetime), POINTER, INTENT(in) :: current_datetime
    REAL(lp),                  INTENT(in) :: dtime
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),           POINTER :: model
    TYPE(t_jsb_grid),            POINTER :: hgrid
    CLASS(t_jsb_tile_abstract),  POINTER :: tile
    CLASS(t_jsb_tile_abstract),  POINTER :: box_tile
    LOGICAL, SAVE                        :: is_read = .FALSE. !< TRUE if forcing file has beed read
    ! -------------------------------------------------------------------------------------------------- !! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_DFIRE_SF_)
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(Q_DFIRE_SF_)
    dsl4jsb_Def_memory_tile(Q_DFIRE_SF_, box_tile)
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real3D_onDomain :: ic_lightning_events_monthly ! out
    dsl4jsb_Real2D_onDomain :: ic_lightning_events         ! out
    ! -------------------------------------------------------------------------------------------------- !
    REAL(lp),       POINTER :: ptr_3D(:,:,:)  ! tmp pointer
    LOGICAL :: is_experiment_start, read_deposition_data

    INTEGER :: current_month, nproma, nblks, i_month, no_omp_thread
    TYPE(t_input_file)          :: input_file
    CHARACTER(len=*), PARAMETER :: routine = modname//':provide_lightning_events'
    ! -------------------------------------------------------------------------------------------------- !
    model => Get_model(model_id)
    CALL model%Get_top_tile(tile)
    no_omp_thread = Get_omp_thread()
    hgrid => Get_grid(model%grid_id)

    nproma = hgrid%nproma
    nblks  = hgrid%nblks
    CALL model%Get_top_tile(box_tile)

    is_experiment_start = is_time_experiment_start(current_datetime)
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_DFIRE_SF_)
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(Q_DFIRE_SF_)
    dsl4jsb_Get_var3D_onDomain(Q_DFIRE_SF_, ic_lightning_events_monthly) ! yearly lightning events
    dsl4jsb_Get_memory_tile(Q_DFIRE_SF_, box_tile)
    dsl4jsb_Get_var2D_onDomain_tile(Q_DFIRE_SF_, ic_lightning_events, box_tile) ! yearly lightning events
    ! -------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')

      ! Check if data needs to be read
    IF (is_experiment_start &
          & .OR. (is_newyear(current_datetime, dtime))) THEN
        read_deposition_data = .TRUE.
    ELSE
        read_deposition_data = .FALSE.
    END IF

    IF (read_deposition_data) THEN
        ! In this case we search a file ending on the current year
        input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(Q_DFIRE_SF_)%ic_lightning_events_filename), model%grid_id)

        !>
        ptr_3D => input_file%Read_2D_time(variable_name='HRMC_COM_FR')
        DO i_month = 1,12
          ! Data is loaded as #events km-2 day-1
          ! Units conversion from km-2 to ha-1
          ic_lightning_events_monthly(:,i_month,:) = ptr_3D(:,:,i_month) / 100.
        END DO
        CALL input_file%Close()

    END IF ! (read_deposition_data) THEN

    ! Assign the deposition values of this month to the A2L variables
    current_month = get_month(current_datetime)

    ic_lightning_events(:,:) = ic_lightning_events_monthly(:,current_month,:)

    IF (debug_on()) CALL message(TRIM(routine), 'Finishing routine')
  END SUBROUTINE provide_lightning_events
#endif

#endif
END MODULE mo_q_dfire_sf_init
