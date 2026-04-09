!> Contains the interfaces to the nlcc process
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
MODULE mo_nlcc_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on, acc_stream
  USE mo_jsb_time,        ONLY: is_time_experiment_start
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes A2L_, PHENO_, ASSIMI_, CARBON_, NLCC_, DISTURB_, FAGE_

  ! Use of process configurations
  dsl4jsb_Use_config(NLCC_)
  dsl4jsb_Use_config(FAGE_)

  ! Use of process memories
!! TR$  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(NLCC_)
  dsl4jsb_Use_memory(DISTURB_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_nlcc_tasks !, t_nlcc_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_nlcc_interface'

  !> Type definition for nlcc task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_nlcc
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_nlcc     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_nlcc  !< Aggregates computed task variables
  END TYPE tsk_nlcc

  !> Constructor interface for nlcc task
  INTERFACE tsk_nlcc
    PROCEDURE Create_task_nlcc                        !< Constructor function for task
  END INTERFACE tsk_nlcc

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for nlcc task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "nlcc"
  !!
  FUNCTION Create_task_nlcc(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_nlcc::return_ptr)
    CALL return_ptr%Construct(name='nlcc', process_id=NLCC_, owner_model_id=model_id)

  END FUNCTION Create_task_nlcc

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for nlcc process
  !!
  !! @param[in,out] this      Instance of nlcc process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_nlcc_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_nlcc(model_id))

  END SUBROUTINE Register_nlcc_tasks

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "nlcc"
  !! Task "nlcc" <EXPLANATIION_OF_PROCESS>.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_nlcc(tile, options)

    USE mo_nlcc_process,   ONLY: calc_climbuf, bioclim_limits, potential_tree_fpc, fpc_to_cover_fract_pot, &
                                 cover_fract_pot_to_cover_fract, desert_fraction, fpc_daily, scale_fpc, fract_small
    USE mo_jsb_time,       ONLY: is_newyear, is_newmonth, is_newday, get_day_length, get_previous_month_length
    USE mo_jsb_grid,       ONLY: Get_grid
    USE mo_jsb_grid_class, ONLY: t_jsb_grid
    USE mo_jsb_tile_class, ONLY: t_jsb_tile_abstract

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER  ::  model
    TYPE(t_jsb_grid),  POINTER  ::  grid
    CLASS(t_jsb_tile_abstract),  POINTER  ::  ptr_tile
    LOGICAL  :: lstart, init_running_means, init_mmtemp20
    INTEGER  :: iblk, ics, ice, nc, no_children
    INTEGER  :: ic, itile
    REAL(wp) :: seconds_day   ! number of seconds per day
    REAL(wp) :: seconds_month ! number of seconds in previous month
    REAL(wp) :: dtime
    REAL(wp) :: accelerate_nlcc
    REAL(wp), ALLOCATABLE :: fuel(:) !!$ TR this should be the sum of the YASSO AWEN pools
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_nlcc'
    REAL(wp), POINTER  ::  land_cover_class(:,:)
    LOGICAL, ALLOCATABLE  :: is_woody(:)
    LOGICAL, ALLOCATABLE  :: is_grass(:)
    LOGICAL, ALLOCATABLE  :: is_pasture(:)
    LOGICAL, ALLOCATABLE  :: is_crop(:)
    LOGICAL, ALLOCATABLE  :: is_dynamic(:)
    REAL(wp), POINTER  ::  bclimit_min_cold_mmtemp(:,:)
    REAL(wp), POINTER  ::  bclimit_max_cold_mmtemp(:,:)
    REAL(wp), POINTER  ::  bclimit_max_warm_mmtemp(:,:)
    REAL(wp), POINTER  ::  bclimit_min_temprange(:,:)
    REAL(wp), POINTER  ::  bclimit_min_gdd(:,:)
    REAL(wp), POINTER  ::  tau_c_woods(:,:)
    REAL(wp), POINTER  ::  NPP_mean_5year(:,:)
    REAL(wp), POINTER  ::  burned_fract(:,:)
    REAL(wp), POINTER  ::  damaged_fract(:,:)
    REAL(wp), POINTER  ::  max_green_bio(:,:)
    REAL(wp), POINTER  ::  sla(:,:)
    REAL(wp), ALLOCATABLE  ::  cover_fract_pot_previous_year(:,:)
    REAL(wp), ALLOCATABLE  ::  cover_fract_sum(:)
    REAL(wp), ALLOCATABLE  :: fract_fpc_max(:,:)

    ! Declare process memories
    dsl4jsb_Def_config(NLCC_)
    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(NLCC_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(DISTURB_)

    ! Declare pointers to variables in nlcc memory
    dsl4jsb_Real2D_onChunk :: temp_sum_day
    dsl4jsb_Real2D_onChunk :: temp_sum_month
    dsl4jsb_Real2D_onChunk :: min_mmtemp_of_yr
    dsl4jsb_Real2D_onChunk :: max_mmtemp_of_yr
    dsl4jsb_Real2D_onChunk :: min_mmtemp20
    dsl4jsb_Real2D_onChunk :: max_mmtemp20
    dsl4jsb_Real3D_onChunk :: act_fpc
    dsl4jsb_Real2D_onChunk :: bare_fpc
    dsl4jsb_Real3D_onChunk :: pot_fpc
    dsl4jsb_Real3D_onChunk :: cover_fract_pot
    dsl4jsb_Real3D_onChunk :: cover_fract
    dsl4jsb_Real3D_onChunk :: bio_exist
    dsl4jsb_Real2D_onChunk :: gdd_sum_year
    dsl4jsb_Real2D_onChunk :: gdd_prev_year
    dsl4jsb_Real2D_onChunk :: sum_green_bio_memory
    dsl4jsb_Real2D_onChunk :: desert_fpc

    ! Declare pointers to variables in foreign memory
    dsl4jsb_Real2D_onChunk :: t_air

    ! Local variables
    LOGICAL :: new_year, new_month, new_day

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(NLCC_)
    dsl4jsb_Get_config(FAGE_)

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(NLCC_)) RETURN

    ! Consistency check
    IF (model%Is_process_enabled(FAGE_)) THEN
      CALL finish(routine, 'Setup ERROR: NLCC is not possible with forest age (FAGE) enabled.')
    END IF

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get process memories
    dsl4jsb_Get_memory(NLCC_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(DISTURB_)

    ! import lct information from PFT tiles
    !TODO: collect_var_real2d has to be replaced by appropriate dsl4jsb invocations!

    ! Define masks for different vegetation types
    no_children = tile%Get_no_of_children()
    land_cover_class => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%land_cover_class)
    !TODO: only 1d vector for is_woody etc
    ALLOCATE(is_woody(no_children))
    ALLOCATE(is_grass(no_children))
    ALLOCATE(is_pasture(no_children))
    ALLOCATE(is_crop(no_children))
    ALLOCATE(is_dynamic(no_children))
    ALLOCATE(cover_fract_sum(nc))
    ALLOCATE(cover_fract_pot_previous_year(nc,no_children))
    ALLOCATE(fuel(nc))
    ALLOCATE(fract_fpc_max(nc,no_children))
    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(is_woody, is_grass, is_pasture, is_crop, is_dynamic) &
    !$ACC   CREATE(cover_fract_sum, cover_fract_pot_previous_year, fuel, fract_fpc_max)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO itile = 1, no_children
      ! Point to the correct PFT
      IF (itile == 1) THEN
        ptr_tile => tile%Get_first_child_tile()
      ELSE
        ptr_tile => ptr_tile%Get_next_sibling_tile()
      END IF
      is_woody(itile)   = .FALSE.
      is_grass(itile)   = .FALSE.
      is_pasture(itile) = .FALSE.
      is_crop(itile)    = .FALSE.
      is_dynamic(itile) = .FALSE.

      IF (  dsl4jsb_Lctlib_param_tile(LandcoverClass, ptr_tile) == 3 .OR. &
          & dsl4jsb_Lctlib_param_tile(LandcoverClass, ptr_tile) == 5)       is_woody(itile)   = .TRUE.
      IF (  dsl4jsb_Lctlib_param_tile(LandcoverClass, ptr_tile) == 4)       is_grass(itile)   = .TRUE.
      IF (  dsl4jsb_Lctlib_param_tile(LandcoverClass, ptr_tile) == 7)       is_pasture(itile) = .TRUE.
      IF (  dsl4jsb_Lctlib_param_tile(LandcoverClass, ptr_tile) == 6)       is_crop(itile)    = .TRUE.
      IF (  dsl4jsb_Lctlib_param_tile(LandcoverClass, ptr_tile) == 3 .OR. &
          & dsl4jsb_Lctlib_param_tile(LandcoverClass, ptr_tile) == 4 .OR. &
          & dsl4jsb_Lctlib_param_tile(LandcoverClass, ptr_tile) == 5)       is_dynamic(itile) = .TRUE.
    END DO
    !$ACC END PARALLEL LOOP

    ! import bioclimatic limits from lct library
    bclimit_min_cold_mmtemp => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%bclimit_min_cold_mmtemp)
    bclimit_max_cold_mmtemp => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%bclimit_max_cold_mmtemp)
    bclimit_max_warm_mmtemp => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%bclimit_max_warm_mmtemp)
    bclimit_min_temprange => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%bclimit_min_temprange)
    bclimit_min_gdd => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%bclimit_min_gdd)
    tau_c_woods => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%tau_c_woods)

    ! import specific leave area from lct library
    sla => carbon__mem%collect_var_real2d(ics, ice, iblk, no_children, carbon__mem%sla)

    ! import annual maximum carbon content of the green pool
    max_green_bio => carbon__mem%collect_var_real2d(ics, ice, iblk, no_children, carbon__mem%max_green_bio)

    ! import 5 year running average of potential NPP
    NPP_mean_5year => assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, assimi__mem%NPP_mean_5year)

    ! import burned fractions and damaged fractions
    burned_fract => disturb__mem%collect_var_real2d(ics, ice, iblk, no_children, disturb__mem%burned_fract)
    damaged_fract => disturb__mem%collect_var_real2d(ics, ice, iblk, no_children, disturb__mem%damaged_fract)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(NLCC_,temp_sum_day)     ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,temp_sum_month)   ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,min_mmtemp_of_yr) ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,max_mmtemp_of_yr) ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,min_mmtemp20)     ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,max_mmtemp20)     ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,gdd_sum_year)     ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,gdd_prev_year)    ! inout
    dsl4jsb_Get_var3D_onChunk(NLCC_,act_fpc)          ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,bare_fpc)         ! inout
    dsl4jsb_Get_var3D_onChunk(NLCC_,pot_fpc)          ! inout
    dsl4jsb_Get_var3D_onChunk(NLCC_,cover_fract_pot)  ! inout
    dsl4jsb_Get_var3D_onChunk(NLCC_,cover_fract)      ! inout
    dsl4jsb_Get_var3D_onChunk(NLCC_,bio_exist)        ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,sum_green_bio_memory) ! inout
    dsl4jsb_Get_var2D_onChunk(NLCC_,desert_fpc)       ! inout
    dsl4jsb_Get_var2D_onChunk(A2L_,   t_air)          ! in

    IF (is_time_experiment_start(options%current_datetime) .OR. dsl4jsb_Config(NLCC_)%init_nlcc) THEN
      lstart = .TRUE.
    ELSE
      lstart = .FALSE.
    END IF
    new_year = is_newyear(options%current_datetime, options%dtime)
    new_month = is_newmonth(options%current_datetime, options%dtime)
    new_day = is_newday(options%current_datetime, options%dtime)

    seconds_day   = get_day_length()
    seconds_month = get_previous_month_length(options%current_datetime)

    ! Scale cover fractions relative to the veg tile fraction
    ! The sum of all child PFT fractions gives the patrent tile fraction, i.e. the veg tile fraction.
    ! TODO: Directly use the veg tile fraction

    ! Import cover fractions from grid
    grid => Get_grid(model%grid_id)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO itile = 1, no_children
      ! Point to the correct PFT
      IF (itile == 1) THEN
        ptr_tile => tile%Get_first_child_tile()
      ELSE
        ptr_tile => ptr_tile%Get_next_sibling_tile()
      END IF
      ! Get cover fractions
      CALL ptr_tile%Get_fraction(ics, ice, iblk, fract=cover_fract(:,itile))
      !$ACC LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic = 1, nc
        cover_fract_sum(ic) = cover_fract_sum(ic) + cover_fract(ic,itile)
      END DO
      !$ACC END LOOP
    END DO ! itile
    !$ACC END PARALLEL LOOP

    ! Do the scaling
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO itile = 1, no_children
      DO ic = 1, nc
        cover_fract(ic,itile) = cover_fract(ic,itile) / MAX(cover_fract_sum(ic), fract_small)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    ! Initialize act_fpc, bare fpc, and pot_fpc
    IF (lstart) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO itile = 1, no_children
        DO ic = 1, nc
         act_fpc(ic,itile) = cover_fract(ic,itile)
         bare_fpc(ic) = 0._wp
         ! TODO: replace no_children by npft, no need to pass as argument
         CALL scale_fpc(routine, nc, iblk, no_children, is_dynamic, act_fpc, bare_fpc)
         pot_fpc(ic,itile) = act_fpc(ic,itile) !!$ TR to be improved (scale, non-woodies = 0, etc.)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    ENDIF

    accelerate_nlcc    = dsl4jsb_Config(NLCC_)%accelerate_nlcc
    init_running_means = dsl4jsb_Config(NLCC_)%init_running_means

    CALL calc_climbuf(lstart, init_running_means, new_day, new_month, new_year, dtime, t_air, &
                      seconds_day, seconds_month, temp_sum_day, temp_sum_month, &
                      min_mmtemp_of_yr, max_mmtemp_of_yr, min_mmtemp20, max_mmtemp20, &
                      gdd_sum_year, gdd_prev_year)

    ! Annual calculations
    IF (new_year) THEN
      DO itile = 1,no_children
         CALL bioclim_limits (is_dynamic(itile), bclimit_min_cold_mmtemp(:,itile), bclimit_max_cold_mmtemp(:,itile),    &
                              bclimit_max_warm_mmtemp(:,itile), bclimit_min_temprange(:,itile), bclimit_min_gdd(:,itile), &
                              min_mmtemp20(:), max_mmtemp20(:), gdd_prev_year(:), bio_exist(:,itile))
      END DO
      IF (.NOT. lstart) THEN
         CALL potential_tree_fpc (nc, no_children, is_dynamic(:), is_woody(:), bio_exist(:,:), &
                                  NPP_mean_5year(:,:), act_fpc(:,:), pot_fpc(:,:))
      ENDIF
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
      DO itile = 1, no_children
        DO ic = 1, nc
          cover_fract_pot_previous_year(ic,itile) = cover_fract_pot(ic,itile)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      CALL fpc_to_cover_fract_pot(nc, no_children, is_dynamic(:), is_woody(:), &
                                   act_fpc(:,:), bare_fpc(:), cover_fract_pot(:,:))
      IF (.NOT. lstart) THEN
        ! new cover fractions
        CALL cover_fract_pot_to_cover_fract(nc, no_children, is_dynamic(:), is_woody(:), is_pasture(:), is_crop(:), &
                                            cover_fract_pot(:,:), cover_fract_pot_previous_year(:,:), cover_fract(:,:))
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO itile = 1, no_children
          ! Point to the correct pft
          IF (itile == 1) THEN
             ptr_tile => tile%Get_first_child_tile()
          ELSE
             ptr_tile => ptr_tile%Get_next_sibling_tile()
          END IF
          ! Set cover fractions - and rescale from veg tile to the box tile fraction
          CALL ptr_tile%Set_fraction(ics, ice, iblk, fract=cover_fract(:,itile) * cover_fract_sum(:))
        END DO ! itile
        !$ACC END PARALLEL LOOP

        ! Calculate desert fraction
        CALL desert_fraction(nc, no_children, init_running_means, accelerate_nlcc,   &
                             is_woody(:), is_dynamic(:),                             &
                             act_fpc(:,:), bare_fpc(:),                              &
                             sla(:,:), max_green_bio(:,:),                           &
                             sum_green_bio_memory(:), desert_fpc(:))

        ! Calculate fract_fpc_max and hand it back to child tiles (pft tiles)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
        DO itile = 1, no_children
          DO ic = 1, nc
            fract_fpc_max(ic,itile) = MAX(fract_small, 1._wp-desert_fpc(ic))
          END DO
        END DO
        !$ACC END PARALLEL LOOP
        CALL pheno__mem%HandDown_var_real2d(ics, ice, iblk, no_children, pheno__mem%fract_fpc_max, fract_fpc_max)
      END IF ! not lstart
    END IF ! new_year

    ! Daily calculations
    IF (new_day) THEN
       CALL fpc_daily(nc, no_children, accelerate_nlcc, is_woody(:), is_dynamic(:), tau_c_woods(:,:),           &
                      NPP_mean_5year(:,:), bio_exist(:,:), burned_fract(:,:), damaged_fract(:,:), pot_fpc(:,:), &
                      act_fpc(:,:), bare_fpc(:), fuel(:))
    ENDIF

    !$ACC END DATA

    DEALLOCATE(is_woody, is_grass, is_pasture, is_crop, is_dynamic)
    DEALLOCATE(bclimit_min_cold_mmtemp, bclimit_max_cold_mmtemp, &
               bclimit_max_warm_mmtemp, bclimit_min_temprange, bclimit_min_gdd, tau_c_woods)
    DEALLOCATE(fract_fpc_max)
    DEALLOCATE(NPP_mean_5year)
    DEALLOCATE(burned_fract, damaged_fract)
    DEALLOCATE(max_green_bio, sla)
    DEALLOCATE(fuel)
    DEALLOCATE(cover_fract_pot_previous_year)
    DEALLOCATE(cover_fract_sum)
    DEALLOCATE(land_cover_class)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_nlcc

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "nlcc"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_nlcc(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(NLCC_)

    ! Local variables
    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract
    INTEGER  :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_nlcc'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(NLCC_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(NLCC_, min_mmtemp20,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, max_mmtemp20,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, gdd_prev_year,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, bare_fpc,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, desert_fpc,            weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(NLCC_, act_fpc,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, pot_fpc,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, cover_fract_pot,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, cover_fract,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(NLCC_, bio_exist,             weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_nlcc

#endif
END MODULE mo_nlcc_interface
