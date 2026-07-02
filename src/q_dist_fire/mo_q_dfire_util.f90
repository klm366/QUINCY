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
MODULE mo_q_dfire_util
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_jsb_math_constants,    ONLY: eps8, eps1

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calculate_time_average_q_dist_fire, calculate_daily_weather, calculate_pure_moving_averages
  PUBLIC :: utils_linspace
  PUBLIC :: dfire_update_pools

  PUBLIC :: logical2real, real2logical

  !< Rate of spread related methods
  PUBLIC :: calc_ros_opt_react_velocity, calc_ros_prop_flux_ratio, calc_ros_wind_factor, calc_ros_fbd

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_util'

CONTAINS

  ! ======================================================================================================= !
  !> calculate time moving averages and daytime averages for Q_DFIRE_
  !>
  SUBROUTINE calculate_time_average_q_dist_fire(tile, options)

    USE mo_jsb_impl_constants,      ONLY: test_false_true
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_lnd_time_averages,       ONLY: mavg_period_dfire, calc_time_mavg
    !------------------------------------------------------------------------------------------------------ !
    dsl4jsb_Use_processes A2L_, Q_DFIRE_
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(Q_DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER         :: model                !< instance of the model
    TYPE(t_jsb_vgrid),      POINTER         :: vgrid_canopy         !< Vertical grid
    INTEGER                                 :: ncanopy              !< number of canopy layers, from vgrid
    INTEGER                                 :: iblk, ics, ice, nc   !< dimensions
    REAL(wp)                                :: dtime                !< timestep
    CHARACTER(len=*), PARAMETER             :: routine = modname//':calculate_time_average_q_dist_fire'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(Q_DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: rain
    dsl4jsb_Real2D_onChunk      :: snow
    ! Q_DFIRE_
    dsl4jsb_Real2D_onChunk      :: precip_qc_mavg
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_DFIRE_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    vgrid_canopy  => Get_vgrid('q_canopy_layer')
    ncanopy       =  vgrid_canopy%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(Q_DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Get_var2D_onChunk(A2L_, rain)
    dsl4jsb_Get_var2D_onChunk(A2L_, snow)
    ! Q_DFIRE_
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, precip_qc_mavg)
    ! ----------------------------------------------------------------------------------------------------- !

    precip_qc_mavg(:) = calc_time_mavg(dtime, precip_qc_mavg(:), (rain(:)+snow(:)), mavg_period_dfire)

  END SUBROUTINE calculate_time_average_q_dist_fire



  ! ======================================================================================================= !
  !>
  !> Subtask to estimate daily weather variables: precipitation, t_air_max and t_air_min
  !>
  !> A day is defined from the lowest point of the sun of the previous day to the current one.
  !> The lowest point of the sun of the current day is estimated from the zenith.
  !>
  SUBROUTINE calculate_daily_weather(tile, options)

    USE mo_kind,                    ONLY: wp
    USE mo_jsb_math_constants,      ONLY: one_day
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model

    ! Use of process memories
    dsl4jsb_Use_processes A2L_, Q_DFIRE_
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(Q_DFIRE_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile               !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)    :: options            !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER         :: model              !< instance of the model
    INTEGER                                 :: iblk, ics, ice, nc, ic
    REAL(wp)                                :: dtime              !< timestep
    LOGICAL                                 :: is_midnight(options%nc) !< is peak of the day
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calculate_daily_weather'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(Q_DFIRE_)
    ! Declare pointers to variables in memory
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: t_air
    dsl4jsb_Real2D_onChunk      :: rain
    dsl4jsb_Real2D_onChunk      :: snow
    dsl4jsb_Real2D_onChunk      :: local_time_day_seconds
    ! Q_DFIRE_ 2D
    dsl4jsb_Real2D_onChunk      :: precip_daytime         !<
    dsl4jsb_Real2D_onChunk      :: t_air_max              !<
    dsl4jsb_Real2D_onChunk      :: t_air_min              !<
    dsl4jsb_Real2D_onChunk      :: precip_daytime_old
    dsl4jsb_Real2D_onChunk      :: t_air_max_old
    dsl4jsb_Real2D_onChunk      :: t_air_min_old
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(Q_DFIRE_)) RETURN   ! ask for all processes considered here ?
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(tile%owner_model_id)
    ! Get process config
    !..
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(Q_DFIRE_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, rain)
    dsl4jsb_Get_var2D_onChunk(A2L_, snow)
    dsl4jsb_Get_var2D_onChunk(A2L_, local_time_day_seconds)         ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, precip_daytime)
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, t_air_max)
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, t_air_min)
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, precip_daytime_old)
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, t_air_max_old)
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, t_air_min_old)

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 1.0 init
    !!
    is_midnight(:) = .FALSE.

    ! 1.1 update min/max values
    t_air_max_old(:) = max(t_air_max_old(:), t_air(:))
    t_air_min_old(:) = min(t_air_min_old(:), t_air(:))
    ! 1.2 update accumulation values
    precip_daytime_old(:) = precip_daytime_old(:) + ((rain(:) + snow(:)) / (one_day / dtime))

    !> 2.0 at the timestep after local midnight (1st timestep of the new day) calc daytime averages
    !> 2.1 1D variables
    is_midnight = (ABS(local_time_day_seconds - dtime) < eps1)

    DO ic=1, nc
      IF (is_midnight(ic)) THEN
        ! check if at least one value had been assigned at the current day, avoiding division by zero
        ! calculate the average over the previous day and pass it to the according *_daytime variable;
        t_air_max(ic)            = t_air_max_old(ic)
        t_air_min(ic)            = t_air_min_old(ic)
        precip_daytime(ic)       = precip_daytime_old(ic)
        ! init variables after min/max daily value
        t_air_max_old(ic)        = t_air(ic)  !< current temperature
        t_air_min_old(ic)        = t_air(ic)
        precip_daytime_old(ic)   = 0._wp      !< no precipitation
      END IF
    END DO

    !    IF (debug_on() .AND. iblk == 1) THEN
    IF (is_midnight(1)) THEN
      ! check if at least one value had been assigned at the current day, avoiding division by zero
      ! calculate the average over the previous day and pass it to the according *_daytime variable;
      IF (t_air_max(1) < t_air_min(1))  THEN
        write(*,*) "error: calc_daily_weather:: t_air_min, t_air_max =", t_air_min, t_air_max
        CALL finish("calc_daily_weather", "t_air_min cant be higher than t_air_max")
      ENDIF
    ENDIF

  END SUBROUTINE calculate_daily_weather


  ! ======================================================================================================= !
  !>
  !> Subtask to calculate pure moving average
  !>
  !>
  SUBROUTINE calculate_pure_moving_averages(tile, options)

    USE mo_kind,                    ONLY: wp
    USE mo_jsb_math_constants,      ONLY: eps8
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_quincy_mavg_class,       ONLY: t_mavg

    ! Use of process memories
    dsl4jsb_Use_processes A2L_, Q_DFIRE_
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(Q_DFIRE_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile               !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)    :: options            !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER         :: model                !< instance of the model
    INTEGER                                 :: iblk, ics, ice, nc, ic
    REAL(wp)                                :: dtime                !< timestep
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calculate_pure_moving_averages'
    TYPE(t_mavg), SAVE, ALLOCATABLE :: precip_mavg_type(:)
    ! ---------------------------
    !$omp threadprivate( precip_mavg_type)
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(Q_DFIRE_)
    ! Declare pointers to variables in memory
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: rain
    dsl4jsb_Real2D_onChunk      :: snow
    ! Q_DFIRE_ 2D
    dsl4jsb_Real2D_onChunk      :: precip_mavg
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(Q_DFIRE_)) RETURN   ! ask for all processes considered here ?
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(tile%owner_model_id)
    ! Get process config
    !..
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(Q_DFIRE_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(A2L_, rain)
    dsl4jsb_Get_var2D_onChunk(A2L_, snow)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(Q_DFIRE_, precip_mavg)
    IF (.NOT. ALLOCATED(precip_mavg_type)) THEN
      ALLOCATE(precip_mavg_type(nc))
      precip_mavg_type = t_mavg(dtime, 1) !< 1=to daily
    ENDIF

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    DO ic = 1, nc
      CALL precip_mavg_type(ic)%add((snow(ic) + rain(ic)))
      precip_mavg(ic) = precip_mavg_type(ic)%get()
    END DO

  END SUBROUTINE calculate_pure_moving_averages


  ! ======================================================================================================= !
  !>
  !> Same method as numpy linspace.
  !> https://gist.github.com/ivan-pi/f4b4741d7ed54ceff787c85d6ba22a5a
  !>
  FUNCTION utils_linspace(start,end,num,endpoint,step) RESULT (samples)

    ! PARAMETERS
    REAL(wp), INTENT(in) :: start
    !! The starting value of the sequence.
    REAL(wp), INTENT(in) :: end
    !! The end value of the sequence, unless `endpoint` is set to `.false.`.
    !! In that case, the sequence consists of all but the last of `num + 1`
    !! evenly spaced samples, so that `end` is excluded. Note that the
    !! step size changes when `endpoint` is `.false.`.
    INTEGER, INTENT(in), OPTIONAL :: num
    !! Number of samples to generate. Default value is 50.
    LOGICAL, INTENT(in), OPTIONAL :: endpoint
    !! If `.true.`, `end` is the last sample. Otherwise, it is not included. Default is `.true.`.
    REAL(wp), INTENT(out), OPTIONAL :: step
    !! If present, `step` is the size of spacing between samples.

    ! RETURNS
    REAL(wp), ALLOCATABLE :: samples(:)
    !! There are `num` equally spaced samples in the closed interval `[start, stop]` or
    !! the half-open interval `[start, stop)` (depending on whether `endpoint` is `.true.` or `.false.`).

    INTEGER :: num_, i
    LOGICAL :: endpoint_
    REAL(wp) :: step_

    num_ = 50
    IF (present(num)) num_ = num

    endpoint_ = .TRUE.
    IF (present(endpoint)) endpoint_ = endpoint

    ! find step size
    IF (endpoint_) then
      step_ = (end - start)/REAL(num_-1,wp)
    ELSE
      step_ = (end - start)/REAL(num_,wp)
    END IF

    IF (present(step)) step = step_

    ALLOCATE(samples(num_))
    DO i = 1, num_
      samples(i) = start + (i-1)*step_
    END DO
  END FUNCTION utils_linspace


  ! ======================================================================================================= !
  !>
  !> Convert logical into real(wp)
  !> TRUE=1._wp, FALSE=0._wp
  !>
  ELEMENTAL PURE REAL(wp) FUNCTION logical2real(val)
    LOGICAL, INTENT(in) :: val

    IF (val) THEN
      logical2real = 1._wp
    ELSE
      logical2real = 0._wp
    ENDIF

  END FUNCTION logical2real


  ! ======================================================================================================= !
  !>
  !> Convert real to logical:
  !> TRUE=1._wp, FALSE=0._wp
  !>
  PURE ELEMENTAL LOGICAL FUNCTION real2logical(val)
    REAL(wp), INTENT(in) :: val

    IF (val == 1._wp) THEN
      real2logical = .TRUE.
    ELSE ! (val == 0._wp)
      real2logical = .FALSE.
    ENDIF
  END FUNCTION real2logical


  ! ======================================================================================================= !
  !>
  !> Apply dfire_update_pools
  !> 1. Remove fire emissions from sb_pool
  !> 2. Transfer litter pools (residue) into polymeric litter
  !>
  SUBROUTINE dfire_update_pools(sb_fire_emissions, sb_fire_litter_ash, sb_pool)
    USE mo_lnd_bgcm_idx, ONLY: ix_woody_litter, ix_polymeric_litter, ix_soluable_litter, &
                               ix_residue, ix_residue_assoc, ix_dom, ix_dom_assoc

    REAL(wp), INTENT(in) :: sb_fire_emissions(:,:,:)
    REAL(wp), INTENT(in) :: sb_fire_litter_ash(:,:,:)
    REAL(wp), INTENT(inout) :: sb_pool(:,:,:)

    !> 1. Remove fire emissions from sb_pool_mt
    sb_pool = sb_pool - sb_fire_emissions

    !> Does it make sense? not really, to keep as purely informative.
    !CALL pool_transfer(sb_fire_litter_ash_mt(ix_polymeric_litter,:,:,df_sf_spq_layer_top_sl), & ! amount (IN)
    !                   sb_pool_mt(ix_polymeric_litter,:,:,df_sf_spq_layer_top_sl), & ! pool from (INOUT)
    !                   sb_pool_mt(ix_polymeric_litter,:,:,df_sf_spq_layer_top_sl)) ! pool to (INOUT)

    !> 2. Woody litter residue is transformed into polymeric
    CALL pool_transfer(sb_fire_litter_ash(ix_woody_litter,:,:), & ! amount (IN)
                        sb_pool(ix_woody_litter,:,:), & ! pool from (INOUT)
                        sb_pool(ix_polymeric_litter,:,:)) ! pool to (INOUT)
    CALL pool_transfer(sb_fire_litter_ash(ix_soluable_litter,:,:), & ! amount (IN)
                        sb_pool(ix_soluable_litter,:,:), & ! pool from (INOUT)
                        sb_pool(ix_polymeric_litter,:,:)) ! pool to (INOUT)
    CALL pool_transfer(sb_fire_litter_ash(ix_residue,:,:), & ! amount (IN)
                        sb_pool(ix_residue,:,:), & ! pool from (INOUT)
                        sb_pool(ix_polymeric_litter,:,:)) ! pool to (INOUT)
    CALL pool_transfer(sb_fire_litter_ash(ix_residue_assoc,:,:), & ! amount (IN)
                        sb_pool(ix_residue_assoc,:,:), & ! pool from (INOUT)
                        sb_pool(ix_polymeric_litter,:,:)) ! pool to (INOUT)
    CALL pool_transfer(sb_fire_litter_ash(ix_dom,:,:), & ! amount (IN)
                        sb_pool(ix_dom,:,:), & ! pool from (INOUT)
                        sb_pool(ix_polymeric_litter,:,:)) ! pool to (INOUT)
    CALL pool_transfer(sb_fire_litter_ash(ix_dom_assoc,:,:), & ! amount (IN)
                        sb_pool(ix_dom_assoc,:,:), & ! pool from (INOUT)
                        sb_pool(ix_polymeric_litter,:,:)) ! pool to (INOUT)
  END SUBROUTINE dfire_update_pools


  ! ======================================================================================================= !
  !>
  !> Small helper routine to transfer matter between pools
  !>
  SUBROUTINE pool_transfer(amount, from, to)
    REAL(wp), INTENT(in) :: amount(:,:)
    REAL(wp), INTENT(inout) :: from(:,:)
    REAL(wp), INTENT(inout) :: to(:,:)

    from = from - amount
    to = to + amount
  END SUBROUTINE pool_transfer


  ! ======================================================================================================= !
  !>
  !> Calculate wind factor (Eq. A5 from Thonicke et al. (2010)), used to calculate rate of spread
  !>
  ELEMENTAL PURE FUNCTION calc_ros_wind_factor(u_f_wind, sigma, beta, beta_op) RESULT (wind_factor)
    REAL(wp), INTENT(in)        :: u_f_wind     !< Forward wind [m min-1]
    REAL(wp), INTENT(in)        :: sigma        !< Surface area to volume ratio of the fuel [cm-1]
    REAL(wp), INTENT(in)        :: beta         !< Packing ratio [-]
    REAL(wp), INTENT(in)        :: beta_op      !< Optimum packing ratio [-]
    REAL(wp)                    :: wind_factor  !< Output [-]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                    :: C_var, B_var, E_var
    REAL(wp)                    :: beta_ratio
    ! ----------------------------------------------------------------------------------------------------- !
    C_var = 7.47_wp * exp(-0.8711_wp * (sigma**0.55_wp)) !< Eq. A8
    B_var = 0.15988_wp * (sigma ** 0.54_wp) !< Eq. A7
    E_var = 0.7515_wp * exp(-0.01094_wp * sigma) !< Eq. A9
    beta_ratio = beta / beta_op
    wind_factor = C_var * ((3.281_wp * u_f_wind) ** B_var) * (beta_ratio ** (-E_var))
  END FUNCTION calc_ros_wind_factor


  ! ======================================================================================================= !
  !>
  !> Calculate propagation flux ratio (Eq. A2 from Thonicke et al. (2010))
  !>
  !> Used to calculate rate of spread
  !>
  ELEMENTAL PURE FUNCTION calc_ros_prop_flux_ratio(sigma, beta) RESULT (prop_flux_ratio)
    REAL(wp), INTENT(in) :: sigma           !< surface-area to volume ratio of the fuel [cm-1]
    REAL(wp), INTENT(in) :: beta            !< packing ratio [-]
    REAL(wp)             :: prop_flux_ratio !< output [-]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: e_calc
    ! ----------------------------------------------------------------------------------------------------- !
    e_calc = (0.792_wp + 3.7597_wp * sqrt(sigma)) * (beta + 0.1_wp)
    prop_flux_ratio = exp(e_calc) / (192._wp + 7.9095_wp * sigma)
  END FUNCTION calc_ros_prop_flux_ratio


  ! ======================================================================================================= !
  !>
  !> Calculate The optimum reaction velocity (Table A1, 1st from Thonicke et al. (2010))
  !>
  !> Used to calculate rate of spread
  !>
  ELEMENTAL PURE FUNCTION calc_ros_opt_react_velocity(max_react_vel, beta, beta_op, A_var) RESULT (opt_react_vel)
    REAL(wp), INTENT(in)        :: max_react_vel   !< Maximum reaction velocity [min-1]
    REAL(wp), INTENT(in)        :: beta            !< packing ratio [kg m-3]
    REAL(wp), INTENT(in)        :: beta_op         !< Eq. A6 [-]
    REAL(wp), INTENT(in)        :: A_var           !< [-]
    REAL(wp)                    :: opt_react_vel   !< output [min-1]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: beta_ratio
    ! ----------------------------------------------------------------------------------------------------- !
    beta_ratio = beta / beta_op
    opt_react_vel = max_react_vel * (beta_ratio ** A_var) * exp(A_var * (1 - beta_ratio))
  END FUNCTION calc_ros_opt_react_velocity


  ! ======================================================================================================= !
  !>
  !> Calculate the weighted average of the Fuel bulk density for the current PFT.
  !> Using the Pfeiffer et al. (2013) paper, equation A17. Remaining equations should be implemented
  !> when accounting for the other pfts.
  !> Thonicke's et al. (2010) is not used as the equation does not seems correct when accounting for
  !> 10 and 100 hours fuel bulk densities (slightly more).
  !>
  !> Used to calculate rate of spread
  !>
  FUNCTION calc_ros_fbd(nfclasses, fuel_bulk_density, dead_fuel) RESULT (fbd)
    USE mo_q_dfire_constants, ONLY: df_sf_fuel_class_factors

    INTEGER, INTENT(in)         :: nfclasses !< The number of used fuel classes
    REAL(wp), INTENT(in)        :: fuel_bulk_density !< [kg m-3]
    REAL(wp), INTENT(in)        :: dead_fuel(nfclasses) !< 1, 10, 100 hours [gDM m-2]
    REAL(wp)                    :: fbd !< output [kg m-3]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                    :: total_dead_fuel
    REAL(wp)                    :: fbd_accounting_fc
    ! ----------------------------------------------------------------------------------------------------- !
    total_dead_fuel = sum(dead_fuel)
    fbd_accounting_fc = fuel_bulk_density * SUM(dead_fuel * df_sf_fuel_class_factors)
    fbd = fbd_accounting_fc / total_dead_fuel
  END FUNCTION calc_ros_fbd
#endif
END MODULE mo_q_dfire_util
