!> Contains the interfaces for the surface energy balance process.
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
!>#### Interfaces of the surface energy balance process
!>
!> The surface energy balance (SEB_) process computes the surface temperature from the balance of
!> the energy fluxes at the surface,
!> ```
!> C dT/dt = R_net + H + LE + G + phase change
!> ```
!> with net radiation R_net, sensible and latent heat fluxes H and LE, the heat flux G from the
!> ground below (`forc_hflx`) and the heat capacity C of the surface layer. All fluxes are positive
!> downward. The balance is solved implicitly, coupled to the vertical diffusion in the atmosphere
!> through the Richtmyer-Morton coefficients (`t_acoef`, `t_bcoef`, `q_acoef`, `q_bcoef`) from A2L_
!> (Schulz et al., 2001, doi:10.1175/1520-0450(2001)040<0642:OTLSAC>2.0.CO;2). The science is in
!> [[mo_seb_land]] (land, including glaciers) and [[mo_seb_lake]] (lakes).
!>
!> The process has four tasks; in a time step they run in this order (see
!> [[mo_jsb_model_init:jsbach_init]]):
!>
!> 1. `surface_energy` ([[update_surface_energy]]): new surface temperature, after the HYDRO_ tasks
!>    `snow_and_wet_fraction` and `soil_properties`,
!> 2. `surface_fluxes` ([[update_surface_fluxes]]): sensible and latent heat fluxes from the new
!>    surface state, after the HYDRO_ task `evaporation`,
!> 3. `asselin_filter` ([[update_asselin]]): time filter of the surface temperature, after all
!>    other tasks of the physics,
!> 4. `energy_balance` ([[update_energy_balance]]): diagnostic residual of the energy fluxes,
!>    directly after `asselin_filter` since it uses the filtered temperature `t`.
!>
!> In addition, [[seb_check_temperature_range]] and [[global_seb_diagnostics]] are called from
!> `jsbach_finish_timestep`.
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_seb_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_jsb_control,     ONLY: debug_on, acc_stream
  USE mo_jsb_time,        ONLY: is_time_experiment_start
  USE mo_kind,            ONLY: lp
  USE mo_exception,       ONLY: message, finish, message_text

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes SEB_, HYDRO_, SSE_, A2L_, RAD_
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(SSE_)

  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(RAD_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_seb_tasks, seb_check_temperature_range
#ifndef __QUINCY_STANDALONE__
  PUBLIC :: global_seb_diagnostics
#endif
#ifdef __QUINCY_STANDALONE__
  ! USEd by the mo_qs_model_interface
  PUBLIC :: update_surface_energy, update_asselin, update_surface_fluxes, update_energy_balance
#endif

  !> Type definition for surface_energy task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_energy
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_energy
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_energy
  END TYPE tsk_surface_energy

  !> Constructor interface for surface_energy task
  INTERFACE tsk_surface_energy
    PROCEDURE Create_task_surface_energy
  END INTERFACE tsk_surface_energy

  !> Type definition for asselin_filter task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_asselin_filter
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_asselin
    PROCEDURE, NOPASS :: Aggregate => aggregate_asselin
  END TYPE tsk_asselin_filter

  !> Constructor interface for asselin_filter task
  INTERFACE tsk_asselin_filter
    PROCEDURE Create_task_asselin_filter
  END INTERFACE tsk_asselin_filter

  !> Type definition for surface_fluxes task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_fluxes
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_fluxes
  END TYPE tsk_surface_fluxes

  !> Constructor interface for surface_fluxes task
  INTERFACE tsk_surface_fluxes
    PROCEDURE Create_task_surface_fluxes
  END INTERFACE tsk_surface_fluxes

  !> Type definition for energy_balance task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_energy_balance
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_energy_balance
    PROCEDURE, NOPASS :: Aggregate => aggregate_energy_balance
  END TYPE tsk_energy_balance

  !> Constructor interface for energy_balance task
  INTERFACE tsk_energy_balance
    PROCEDURE Create_task_energy_balance
  END INTERFACE tsk_energy_balance

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_interface'

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_energy task
  !>
  FUNCTION Create_task_surface_energy(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id   !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr !< Instance of process task "surface_energy"

    ALLOCATE(tsk_surface_energy::return_ptr)
    CALL return_ptr%Construct(name='surface_energy', process_id=SEB_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_energy

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for asselin_filter task
  !>
  FUNCTION Create_task_asselin_filter(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id   !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr !< Instance of process task "asselin_filter"

    ALLOCATE(tsk_asselin_filter::return_ptr)
    CALL return_ptr%Construct(name='asselin_filter', process_id=SEB_, owner_model_id=model_id)

  END FUNCTION Create_task_asselin_filter

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_fluxes task
  !>
  FUNCTION Create_task_surface_fluxes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id   !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr !< Instance of process task "surface_fluxes"

    ALLOCATE(tsk_surface_fluxes::return_ptr)
    CALL return_ptr%Construct(name='surface_fluxes', process_id=SEB_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_fluxes

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for energy_balance task (energy balance residual)
  !>
  FUNCTION Create_task_energy_balance(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id   !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr !< Instance of process task "energy_balance"

    ALLOCATE(tsk_energy_balance::return_ptr)
    CALL return_ptr%Construct(name='energy_balance', process_id=SEB_, owner_model_id=model_id)

  END FUNCTION Create_task_energy_balance

  ! =======================================================================================================
  !> Register tasks for surface energy process
  !>
  SUBROUTINE Register_seb_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this     !< Instance of surface energy process class
    INTEGER,              INTENT(in)    :: model_id !< Model id

    CALL this%Register_task(tsk_surface_energy(model_id))
    CALL this%Register_task(tsk_asselin_filter(model_id))
    CALL this%Register_task(tsk_surface_fluxes(model_id))
    CALL this%Register_task(tsk_energy_balance(model_id))

  END SUBROUTINE Register_seb_tasks

  ! ================================================================================================================================
  !>
  !>#### Implementation of "update" for task "surface_energy"
  !>
  !> Task "surface_energy" calculates the new surface temperature from the surface energy balance
  !> (net radiation, sensible and latent heat, ground heat, phase changes). Depending on the tile,
  !> [[mo_seb_lake:update_surface_energy_lake]] or [[mo_seb_land:update_surface_energy_land]] is
  !> called; with `use_tmx`, the latter is used on all tiles but the lake.
  !>
  SUBROUTINE update_surface_energy(tile, options)

    USE mo_seb_land, ONLY: update_surface_energy_land
    USE mo_seb_lake, ONLY: update_surface_energy_lake

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile    !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< Additional run-time parameters

    TYPE(t_jsb_model), POINTER :: model !< This instance of ICON-Land

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_energy'

    IF (.NOT. tile%Is_process_calculated(SEB_)) RETURN

    model => Get_model(tile%owner_model_id)

    IF (tile%is_lake) THEN
      CALL update_surface_energy_lake(tile, options)
    ELSE IF (tile%is_land .AND. .NOT. model%config%use_tmx) THEN
      CALL update_surface_energy_land(tile, options)
    ELSE IF (model%config%use_tmx) THEN
      CALL update_surface_energy_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type '//TRIM(tile%lcts(1)%name)//' on tile '//TRIM(tile%name))
    END IF

  END SUBROUTINE update_surface_energy

  ! -------------------------------------------------------------------------------------------------------
  !>
  !> Implementation of "aggregate" for task "surface_energy"
  !>
  !> Aggregates the surface temperatures, the surface state passed to the atmosphere and the heat
  !> capacity and heat fluxes of the surface layer, weighted by tile fractions.
  !>
  SUBROUTINE aggregate_surface_energy(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile    !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< Additional run-time parameters

    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_config(SEB_)

    TYPE(t_jsb_model), POINTER :: model                   !< This instance of ICON-Land
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract !< Aggregation method: area weighted means

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_surface_energy'

    INTEGER :: &
      & iblk, & !< Current block index
      & ics,  & !< Index of first cell in chunk
      & ice     !< Index of last cell in chunk

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(SEB_, t,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, t_old,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, t_unfilt,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, t_srf,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, t_eff4,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, qsat_star,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, s_star,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, forc_hflx,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, heat_cap,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, le_phase_change, weighted_by_fract)

    IF (model%config%use_lakes) THEN
      dsl4jsb_Aggregate_onChunk(SEB_, t_lwtr,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SEB_, s_lwtr,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SEB_, qsat_lwtr,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SEB_, fract_lice, weighted_by_fract)
      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Aggregate_onChunk(SEB_, t_lice,     weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(SEB_, s_lice,     weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(SEB_, qsat_lice,  weighted_by_fract)
      END IF
    END IF

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_energy

  ! ================================================================================================================================
  !>
  !>#### Implementation of "update" for task "asselin_filter"
  !>
  !> Task "asselin_filter" applies the Asselin time filter to the new surface temperature, if
  !> applicable, see [[mo_seb_land:update_asselin_land]]. There is no filter on lake tiles.
  !>
  SUBROUTINE update_asselin(tile, options)

    USE mo_seb_land, ONLY: update_asselin_land

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile    !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< Additional run-time parameters

    TYPE(t_jsb_model), POINTER :: model !< This instance of ICON-Land

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_asselin'

    IF (.NOT. tile%Is_process_calculated(SEB_)) RETURN

    model => Get_model(tile%owner_model_id)

    IF (tile%contains_lake) THEN
      ! No Asselin filter for lake tile
    ELSE IF (tile%is_land .AND. .NOT. model%config%use_tmx) THEN
      CALL update_asselin_land(tile, options)
    ELSE IF (model%config%use_tmx) THEN
      CALL update_asselin_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_asselin

  ! -------------------------------------------------------------------------------------------------------
  !>
  !> Implementation of "aggregate" for task "asselin_filter"
  !>
  SUBROUTINE aggregate_asselin(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile    !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< Additional run-time parameters

    dsl4jsb_Def_memory(SEB_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract !< Aggregation method: area weighted means

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_asselin'

    INTEGER :: &
      & iblk, &                         !< Current block index
      & ics,  &                         !< Index of first cell in chunk
      & ice                             !< Index of last cell in chunk
    TYPE(t_jsb_model), POINTER :: model !< This instance of ICON-Land

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(SEB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(SEB_, t,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, t_filt, weighted_by_fract)
    IF (model%config%use_tmx) THEN
      dsl4jsb_Aggregate_onChunk(SEB_, t_rad4, weighted_by_fract)
    END IF

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_asselin

  ! ================================================================================================================================
  !>
  !>#### Implementation of "update" for task "surface_fluxes"
  !>
  !> Task "surface_fluxes" computes the sensible and latent heat fluxes from the new surface state
  !> and the evaporation computed by HYDRO_, see [[mo_seb_land:update_surface_fluxes_land]] and
  !> [[mo_seb_lake:update_surface_fluxes_lake]].
  !>
  SUBROUTINE update_surface_fluxes(tile, options)

    USE mo_seb_land, ONLY: update_surface_fluxes_land
    USE mo_seb_lake, ONLY: update_surface_fluxes_lake

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile    !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< Additional run-time parameters

    TYPE(t_jsb_model), POINTER :: model !< This instance of ICON-Land

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_fluxes'

    IF (.NOT. tile%Is_process_calculated(SEB_)) RETURN

    model => Get_model(tile%owner_model_id)

    IF (tile%contains_lake) THEN
      CALL update_surface_fluxes_lake(tile, options)
    ELSE IF (tile%is_land .AND. .NOT. model%config%use_tmx) THEN
      CALL update_surface_fluxes_land(tile, options)
    ELSE IF (model%config%use_tmx) THEN
      CALL update_surface_fluxes_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_surface_fluxes

  ! -------------------------------------------------------------------------------------------------------
  !>
  !> Implementation of "aggregate" for task "surface_fluxes"
  !>
  SUBROUTINE aggregate_surface_fluxes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile    !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< Additional run-time parameters

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)

    TYPE(t_jsb_model),       POINTER :: model             !< This instance of ICON-Land
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract !< Aggregation method: area weighted means

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_surface_fluxes'

    INTEGER :: &
      & iblk, & !< Current block index
      & ics,  & !< Index of first cell in chunk
      & ice     !< Index of last cell in chunk

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(SEB_,   latent_hflx,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_,   sensible_hflx,      weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(SEB_,   latent_hflx_lnd,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_,   sensible_hflx_lnd,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_,   latent_hflx_wtr,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_,   sensible_hflx_wtr,  weighted_by_fract)
    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      dsl4jsb_Aggregate_onChunk(SEB_,   latent_hflx_ice,    weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SEB_,   sensible_hflx_ice,  weighted_by_fract)
    END IF

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_fluxes

  ! -------------------------------------------------------------------------------------------------------
  !>
  !>#### Implementation of "update" for task "energy_balance"
  !>
  !> This diagnostic computes `net_srf_rad_balance`, the residual of all energy fluxes into the
  !> surface: net radiation, sensible and latent heat fluxes, `forc_hflx` and the latent energy of
  !> phase changes. The change of heat storage is not taken into account, so a large residual does
  !> not necessarily indicate a problem. If `eb_threshold` > 0, all fluxes are printed for cells
  !> with a residual above the threshold or a surface temperature outside 50 to 400 K (CPU only).
  !> Nothing is computed at the first time step of the experiment.
  !>
  SUBROUTINE update_energy_balance(tile, options)

    USE mo_jsb_grid,               ONLY: Get_grid
    USE mo_jsb_grid_class,         ONLY: t_jsb_grid

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile    !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< Additional run-time parameters

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model !< This instance of ICON-Land
    TYPE(t_jsb_grid),  POINTER    :: grid  !< Horizontal grid

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(RAD_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & rad_srf_net,        &               !< Surface net radiation                               [W m-2]
      & lw_srf_down,        &               !< Longwave downwards radiation                        [W m-2]
      & swvis_srf_down,     &               !< Shortwave visible downwards radiation               [W m-2]
      & swnir_srf_down,     &               !< Shortwave near infrared downwards radiation         [W m-2]
      & swpar_srf_down,     &               !< Shortwave downwards photosynthetic active radiation [W m-2]
      & lw_srf_net,         &               !< Longwave net radiation on land surface              [W m-2]
      & sw_srf_net,         &               !< Shortwave net radiation on land surface             [W m-2]
      & sensible_hflx,      &               !< Sensible heat flux                                  [W m-2]
      & latent_hflx,        &               !< Latent heat flux                                    [W m-2]
      & forc_hflx,          &               !< Ground heat and residual phase change flux          [W m-2]
      & le_phase_change,    &               !< Latent energy required/released due to phase change [J m-2]
      & t,                  &               !< Current time steps temperature                          [K]
      & t_old,              &               !< Previous time steps temperature                         [K]
      & net_srf_rad_balance                 !< Net surface radiation balance                       [W m-2]

    REAL(lp), POINTER :: &
      & tile_fract(:),    & !< Tile fraction []
      & lat(:),           & !< Latitude of the cells [deg]
      & lon(:)              !< Longitude of the cells [deg]

    ! Locally allocated vectors
    !

    REAL(lp), DIMENSION(options%nc) :: &
      & lw_srf_up,       & !< Upward longwave radiation (< 0) [W m-2]
      & sw_srf_up,       & !< Upward (reflected) shortwave radiation (< 0) [W m-2]
      & sw_srf_down,     & !< Downward shortwave radiation [W m-2]
      & delta_t_srf        !< Change of the surface temperature in the time step [K]

    INTEGER :: &
      & iblk, & !< Current block index
      & ics,  & !< Index of first cell in chunk
      & ice,  & !< Index of last cell in chunk
      & nc,   & !< Number of cells in chunk
      & ic      !< Cell index

    INTEGER :: &
      & n_cells, &                   !< Number of cells exceeding eb_threshold or the temperature range
      & ie                           !< Index in the list of these cells
    INTEGER  :: cell_idx(options%nc) !< Indices of these cells in the chunk

    REAL(lp) :: dtime               !< Time step length [s]
    LOGICAL  :: is_experiment_start !< First time step of the experiment

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_energy_balance'
    CHARACTER(len=5000)         :: message_text_long !< Text of the diagnostic message

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime
    is_experiment_start = is_time_experiment_start(options%current_datetime)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    tile_fract => tile%fract(ics:ice, iblk)
    lat        => grid%lat(ics:ice, iblk)
    lon        => grid%lon(ics:ice, iblk)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(RAD_)

    dsl4jsb_Get_var2D_onChunk(A2L_,      lw_srf_down)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      swvis_srf_down)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      swnir_srf_down)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      swpar_srf_down)             ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,      rad_srf_net)                ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,      lw_srf_net)                 ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,      sw_srf_net)                 ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      latent_hflx)                ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      sensible_hflx)              ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      forc_hflx)                  ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      le_phase_change)            ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t)                          ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_old)                      ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      net_srf_rad_balance)        ! out

    n_cells = 0   ! initialize counter for cells exceeding energy balance threshold

    !$ACC DATA CREATE(lw_srf_up, sw_srf_up, sw_srf_down, delta_t_srf)  ASYNC(acc_stream)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      ! Compute fluxes used in log message
      lw_srf_up(ic)   = lw_srf_net(ic) - lw_srf_down(ic)
      sw_srf_down(ic) = swvis_srf_down(ic) + swnir_srf_down(ic) + swpar_srf_down(ic)
      sw_srf_up(ic)   = sw_srf_net(ic) - sw_srf_down(ic)

      ! Compute surface/skin temperature change
      delta_t_srf(ic) = t(ic) - t_old(ic)

      IF (.NOT. is_experiment_start .AND. tile_fract(ic) > EPSILON(1.0_lp)) THEN
        ! Compute energy flux residual (without surface/skin temperature change)
        net_srf_rad_balance(ic) = rad_srf_net(ic) + sensible_hflx(ic) + latent_hflx(ic) &
          &                     + forc_hflx(ic) + le_phase_change(ic) / dtime
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    ! Note that contrary to the water balance, the net surface radiation balance just
    ! checks the residual of the fluxes while ignoring the heat storages. Therefore,
    ! a large residual can still be acceptable if the heat storage is large enough to
    ! compensate for it, and large values do not necessarily indicate a problem.
#ifndef _OPENACC
    DO ic=1,nc
      IF (dsl4jsb_Config(SEB_)%eb_threshold > 0._lp) THEN
        IF (ABS(net_srf_rad_balance(ic)) > dsl4jsb_Config(SEB_)%eb_threshold .OR. &
          & ((t(ic) > 400._lp .OR. t(ic) < 50._lp).AND. tile_fract(ic) > EPSILON(1.0_lp))) THEN
          n_cells = n_cells + 1       ! Count number of cells with energy balance issues of this time step
          cell_idx(n_cells) = ic      ! Index of the respective grid cell (on chunk)
        END IF
      END IF
    END DO

    DO ie = 1, n_cells
      ic = cell_idx(ie)

      WRITE (message_text_long,*) 'Net surface radiation balance [W m-2]', NEW_LINE('a'), &
        & 'on ',TRIM(tile%name),' tile at', lat(ic),'N and ',lon(ic),'E',  NEW_LINE('a'), &
        & '(ic: ',ic,' iblk: ',iblk, ' tile_fract: ',tile_fract(ic),'):',  NEW_LINE('a'), &
        & 'Net surface radiation balance: ', net_srf_rad_balance(ic),      NEW_LINE('a'), &
        & 'Longwave downwards:            ', lw_srf_down(ic),              NEW_LINE('a'), &
        & 'Longwave upwards:              ', lw_srf_up(ic),                NEW_LINE('a'), &
        & 'Shortwave downwards:           ', sw_srf_down(ic),              NEW_LINE('a'), &
        & 'Shortwave upwards:             ', sw_srf_up(ic),                NEW_LINE('a'), &
        & 'Net surface radiation:         ', rad_srf_net(ic),              NEW_LINE('a'), &
        & 'Sensible heat flux:            ', sensible_hflx(ic),            NEW_LINE('a'), &
        & 'Latent heat flux:              ', latent_hflx(ic),              NEW_LINE('a'), &
        & 'Ground heat flux:              ', forc_hflx(ic),                NEW_LINE('a'), &
        & 'Phase change flux:             ', le_phase_change(ic) / dtime,  NEW_LINE('a'), &
        & 'Temperature change [K]:        ', delta_t_srf(ic)
      CALL message (TRIM(routine), message_text_long, all_print=.TRUE.)
    END DO
#endif

    !$ACC END DATA

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_energy_balance

  ! -------------------------------------------------------------------------------------------------------
  !>
  !> Implementation of "aggregate" for task "energy_balance"
  !>
  !> The residual is not aggregated but computed with [[update_energy_balance]] on each tile.
  !>
  SUBROUTINE aggregate_energy_balance(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile    !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< Additional run-time parameters

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_energy_balance'

    INTEGER :: iblk !< Current block index

    iblk = options%iblk

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Don't aggregate, but explicitely compute energy balance on each tile
    CALL update_energy_balance(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_energy_balance

  !-----------------------------------------------------------------------------------------------------
  !> Make sure temperature is within a realistic range and finish the simulation with a meaningful
  !> message if not.
  !>
  !> The routine is called at the end of each time step from `jsbach_finish_timestep`, if the SEB_
  !> process is active. It stops the model if the surface temperature `t` on any tile with
  !> non-zero fraction is below 50 K or above 400 K.
  !>
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE seb_check_temperature_range(model_id)

    USE mo_jsb_grid,            ONLY: Get_grid
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid

    ! Argument
    INTEGER, INTENT(in) :: model_id !< Model id

    ! Local variables
    TYPE(t_jsb_model), POINTER           :: model !< This instance of ICON-Land
    TYPE(t_jsb_grid),  POINTER           :: grid  !< Horizontal grid
    CLASS(t_jsb_tile_abstract), POINTER  :: tile  !< Current tile

    dsl4jsb_Def_memory(SEB_)

    INTEGER :: &
      & iblk,  &           !< Current block index
      & ic,    &           !< Cell index
      & itile              !< Tile index
    REAL(lp), POINTER :: &
      & lat(:,:),        & !< Latitude of the cells [deg]
      & lon(:,:),        & !< Longitude of the cells [deg]
      & tile_fract(:,:)    !< Tile fraction []
    CHARACTER(len=*),  PARAMETER  :: routine = modname//':seb_check_temperature_range'

    dsl4jsb_Real2D_onDomain :: t !< Surface temperature [K]

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    model => Get_model(model_id)
    grid  => Get_grid(model%grid_id)
    lat   => grid%lat(:,:)
    lon   => grid%lon(:,:)

#ifndef __QUINCY_STANDALONE__
    DO itile=1,SIZE(model%tiles_postorder)
      tile => model%tiles_postorder(itile)%p
      IF (.NOT. tile%Is_process_active(SEB_)) CYCLE
#endif

      dsl4jsb_Get_memory(SEB_)
      dsl4jsb_Get_var2D_onDomain(SEB_, t)

      tile_fract => tile%fract(:,:)
      IF (ANY(t(:,:) < 50._lp .OR. t(:,:) > 400._lp)) THEN
        DO iblk = 1, grid%nblks
          DO ic = 1, grid%nproma
            IF ((t(ic,iblk) < 50._lp .OR. t(ic,iblk) > 400._lp) .AND. tile_fract(ic,iblk) > 0._lp) THEN
              WRITE (message_text,*) 'Temperature out of bound: ', t(ic,iblk), 'K',             NEW_LINE('a'), &
                & 'on ',TRIM(tile%name),' tile at', lat(ic,iblk), 'N and ', lon(ic,iblk), 'E',  NEW_LINE('a'), &
                & '(ic: ',ic,' iblk: ',iblk, ' tile_fract:', tile_fract(ic,iblk),'):',          NEW_LINE('a'), &
                & 'One thing to check: consistency of land-sea mask and forcing.'
              CALL finish(TRIM(routine), TRIM(message_text))
            END IF
          END DO
        END DO
      END IF
#ifndef __QUINCY_STANDALONE__
    ENDDO
#endif
  END SUBROUTINE seb_check_temperature_range

  !-----------------------------------------------------------------------------------------------------
  !> Calculations of diagnostic global land mean output
  !>
  !> Computes the global land mean surface temperature `t_gmean` on the root tile, if it is in
  !> the output. The routine is called from `jsbach_finish_timestep`, after the loop over the
  !> nproma blocks.
  !>
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE global_seb_diagnostics(tile, global_land_area, land_area_scaling)

    USE mo_util,                    ONLY: scaled_global_sum_2D
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile                   !< Current tile
    REAL(lp),                   INTENT(in) :: global_land_area       !< Global land area
    REAL(lp),                   INTENT(in) :: land_area_scaling(:,:) !< Land area scaling factor (accounting for halo cells)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SEB_)

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_seb_diagnostics'

    dsl4jsb_Real2D_onDomain :: t          !< Surface temperature on the root tile [K]
    REAL(lp), POINTER       :: t_gmean(:) !< Global land mean surface temperature [K]
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_var2D_onDomain(SEB_,  t) ! in
    t_gmean => SEB__mem%t_gmean%ptr(:)   ! out

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')
    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    IF (SEB__mem%t_gmean%is_in_output) THEN
      CALL scaled_global_sum_2D(t_gmean, t, land_area_scaling, 1.0_lp / global_land_area, acc_stream)
      !$ACC UPDATE ASYNC(acc_stream) DEVICE(SEB__mem%t_gmean%ptr(:))
    END IF

  END SUBROUTINE global_seb_diagnostics

  !
  ! ================================================================================================================================
  !
#endif
END MODULE mo_seb_interface
