!> Contains the interfaces to the qlcc process
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_qlcc_interface
! #ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
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
!   dsl4jsb_Use_processes A2L_, PHENO_, ASSIMI_, CARBON_, QLCC_, DISTURB_
  dsl4jsb_Use_processes A2L_, QLCC_, Q_ASSIMI_, VEG_, DISTURB_

  ! Use of process configurations
  dsl4jsb_Use_config(QLCC_)

  ! Use of process memories
!! TR$  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
  dsl4jsb_Use_memory(A2L_)
!   dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(Q_ASSIMI_)
!   dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(QLCC_)
  dsl4jsb_Use_memory(DISTURB_)
  dsl4jsb_Use_memory(VEG_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_qlcc_tasks !, t_qlcc_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_qlcc_interface'

  !> Type definition for qlcc task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_qlcc
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_qlcc     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_qlcc  !< Aggregates computed task variables
  END TYPE tsk_qlcc

  !> Constructor interface for qlcc task
  INTERFACE tsk_qlcc
    PROCEDURE Create_task_qlcc                        !< Constructor function for task
  END INTERFACE tsk_qlcc

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for qlcc task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "qlcc"
  !!
  FUNCTION Create_task_qlcc(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_qlcc::return_ptr)
    CALL return_ptr%Construct(name='qlcc', process_id=QLCC_, owner_model_id=model_id)

  END FUNCTION Create_task_qlcc

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for qlcc process
  !!
  !! @param[in,out] this      Instance of qlcc process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_qlcc_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_qlcc(model_id))

  END SUBROUTINE Register_qlcc_tasks

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "qlcc"
  !! Task "qlcc" <EXPLANATIION_OF_PROCESS>.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_qlcc(tile, options)

    USE mo_qlcc_process,   ONLY: calc_q_climbuf, bioclim_q_limits, potential_tree_fpc, fpc_to_cover_fract_pot, &
                                 cover_fract_pot_to_cover_fract, desert_fraction, scale_fpc, fpc_daily
    USE mo_jsb_time,       ONLY: is_newyear, is_newmonth, is_newday
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
    LOGICAL  :: lstart, init_running_means !, init_mmtemp20
    INTEGER  :: iblk, ics, ice, nc, no_children
    INTEGER  :: i, itile
    REAL(wp) :: dtime
    REAL(wp) :: accelerate_qlcc
    REAL(wp), ALLOCATABLE :: fuel(:) !!$ TR this should be the sum of the YASSO AWEN pools
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_qlcc'
    REAL(wp), POINTER  ::  var_collect(:,:)
    LOGICAL, ALLOCATABLE  :: is_woody(:,:)
    LOGICAL, ALLOCATABLE  :: is_grass(:,:)
    LOGICAL, ALLOCATABLE  :: is_pasture(:,:)
    LOGICAL, ALLOCATABLE  :: is_crop(:,:)
    LOGICAL, ALLOCATABLE  :: is_dynamic(:,:)
    LOGICAL, ALLOCATABLE  :: is_bare(:,:)
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
    dsl4jsb_Def_config(QLCC_)
    dsl4jsb_Def_memory(QLCC_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Def_memory(DISTURB_)
    dsl4jsb_Def_memory(VEG_)

    ! Declare pointers to variables in qlcc memory
    dsl4jsb_Real2D_onChunk :: seconds_day
    dsl4jsb_Real2D_onChunk :: seconds_month
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
    dsl4jsb_Real3D_onChunk :: npp_5yr
    dsl4jsb_Real3D_onChunk :: max_bio
    dsl4jsb_Real3D_onChunk :: bio_exist
    dsl4jsb_Real2D_onChunk :: gdd_sum_year
    dsl4jsb_Real2D_onChunk :: gdd_prev_year
    dsl4jsb_Real2D_onChunk :: sum_green_bio_memory
    dsl4jsb_Real2D_onChunk :: desert_fpc
    dsl4jsb_Real2D_onChunk :: t_air

    ! Local variables
    LOGICAL :: new_year, new_month, new_day
    INTEGER :: bare_tile_no                  ! number of bare tile

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    !X if necessary: steplen = options%steplen

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    no_children = tile%Get_no_of_children()
    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(QLCC_)

    ! Get process memories
    dsl4jsb_Get_memory(QLCC_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(DISTURB_)

    ! import lct information from PFT tiles
    ! build masks for different vegetation types (is_dynamic, is_woody, etc)
    var_collect => q_assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, q_assimi__mem%land_cover_class)
    ALLOCATE(is_woody(nc,no_children))
    ALLOCATE(is_grass(nc,no_children))
    ALLOCATE(is_pasture(nc,no_children))
    ALLOCATE(is_crop(nc,no_children))
    ALLOCATE(is_dynamic(nc,no_children))
    ALLOCATE(is_bare(nc,no_children))
    is_woody(:,:) = .FALSE.
    is_grass(:,:) = .FALSE.
    is_pasture(:,:) = .FALSE.
    is_crop(:,:) = .FALSE.
    is_dynamic(:,:) = .FALSE.
    is_bare(:,:) = .FALSE.
    DO i = 1,nc
      DO itile = 1,no_children
        IF (NINT(var_collect(i,itile)) == 3 .OR. NINT(var_collect(i,itile)) == 5) is_woody(i,itile) = .TRUE.
        IF (NINT(var_collect(i,itile)) == 4) is_grass(i,itile) = .TRUE.
        IF (NINT(var_collect(i,itile)) == 7) is_pasture(i,itile) = .TRUE.
        IF (NINT(var_collect(i,itile)) == 6) is_crop(i,itile) = .TRUE.
        IF (NINT(var_collect(i,itile)) == 0) is_bare(i,itile) = .TRUE.
        IF (NINT(var_collect(i,itile)) == 3 .OR. NINT(var_collect(i,itile)) == 4 .OR. &
            NINT(var_collect(i,itile)) == 5) is_dynamic(i,itile) = .TRUE.
      END DO
    END DO
    DO itile = 1,no_children
      IF (ANY(is_bare(:,itile))) bare_tile_no = itile
    END DO

    ! import bioclimatic limits from lct library
    bclimit_min_cold_mmtemp => q_assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, q_assimi__mem%bclimit_min_cold_mmtemp)
    bclimit_max_cold_mmtemp => q_assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, q_assimi__mem%bclimit_max_cold_mmtemp)
    bclimit_max_warm_mmtemp => q_assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, q_assimi__mem%bclimit_max_warm_mmtemp)
    bclimit_min_temprange => q_assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, q_assimi__mem%bclimit_min_temprange)
    bclimit_min_gdd => q_assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, q_assimi__mem%bclimit_min_gdd)
    tau_c_woods => q_assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, q_assimi__mem%tau_c_woods)

    ! import specific leave area from lct library
    sla => q_assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, q_assimi__mem%sla)

    ! import annual maximum carbon content of the green pool
    max_green_bio => veg__mem%collect_var_real2d(ics, ice, iblk, no_children, veg__mem%max_green_bio)

    ! import 5 year running average of potential NPP
    NPP_mean_5year => q_assimi__mem%collect_var_real2d(ics, ice, iblk, no_children, q_assimi__mem%NPP_mean_5year)

    ! import burned fractions and damaged fractions
    burned_fract => disturb__mem%collect_var_real2d(ics, ice, iblk, no_children, disturb__mem%burned_fract)
    damaged_fract => disturb__mem%collect_var_real2d(ics, ice, iblk, no_children, disturb__mem%damaged_fract)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(QLCC_, seconds_day)      ! inout
    dsl4jsb_Get_var2D_onChunk(QLCC_, seconds_month)    ! inout
    dsl4jsb_Get_var2D_onChunk(QLCC_, temp_sum_day)     ! inout
    dsl4jsb_Get_var2D_onChunk(QLCC_, temp_sum_month)   ! inout
    dsl4jsb_Get_var2D_onChunk(QLCC_, min_mmtemp_of_yr) ! inout
    dsl4jsb_Get_var2D_onChunk(QLCC_, max_mmtemp_of_yr) ! inout
    dsl4jsb_Get_var2D_onChunk(QLCC_, min_mmtemp20)     ! inout
    dsl4jsb_Get_var2D_onChunk(QLCC_, max_mmtemp20)     ! inout
    dsl4jsb_Get_var2D_onChunk(QLCC_, gdd_sum_year)     ! inout
    dsl4jsb_Get_var2D_onChunk(QLCC_, gdd_prev_year)    ! inout
    dsl4jsb_Get_var3D_onChunk(QLCC_, act_fpc)          ! inout
    dsl4jsb_Get_var2D_onChunk(QLCC_, bare_fpc)         ! inout
    dsl4jsb_Get_var3D_onChunk(QLCC_, pot_fpc)          ! inout
    dsl4jsb_Get_var3D_onChunk(QLCC_, cover_fract_pot)  ! inout
    dsl4jsb_Get_var3D_onChunk(QLCC_, cover_fract)      ! inout
    dsl4jsb_Get_var3D_onChunk(QLCC_, npp_5yr)          ! inout
    dsl4jsb_Get_var3D_onChunk(QLCC_, max_bio)          ! inout
    dsl4jsb_Get_var3D_onChunk(QLCC_, bio_exist)        ! inout
    dsl4jsb_Get_var2D_onChunk(QLCC_, sum_green_bio_memory) ! inout
    dsl4jsb_Get_var2D_onChunk(QLCC_, desert_fpc)       ! inout
    dsl4jsb_Get_var2D_onChunk(A2L_,  t_air)            ! in

    IF (is_time_experiment_start(options%current_datetime) .OR. dsl4jsb_Config(QLCC_)%init_qlcc) THEN
      lstart = .TRUE.
    ELSE
      lstart = .FALSE.
    END IF
    new_year = is_newyear(options%current_datetime, options%dtime)
    new_month = is_newmonth(options%current_datetime, options%dtime)
    new_day = is_newday(options%current_datetime, options%dtime)

    ! qlcc options
    accelerate_qlcc    = dsl4jsb_Config(QLCC_)%accelerate_qlcc
    init_running_means = dsl4jsb_Config(QLCC_)%init_running_means


    ! import cover fractions from grid
!     ALLOCATE (cover_fract_sum(nc))
!     cover_fract_sum(:) = 0._wp
    DO itile = 1,no_children
      ! Point to the correct pft
      IF (itile == 1) THEN
        ptr_tile => tile%Get_first_child_tile()
      ELSE
        ptr_tile => ptr_tile%Get_next_sibling_tile()
      END IF
      ! get cover fractions
      CALL ptr_tile%Get_fraction(ics, ice, iblk, fract=cover_fract(:,itile))
!       cover_fract_sum(:) = cover_fract_sum(:) + cover_fract(:,itile)
    END DO ! itile

    ! set diagnostic 5 year NPP
    DO itile = 1,no_children
      npp_5yr(:,itile) = NPP_mean_5year(:,itile)
      max_bio(:,itile) = max_green_bio(:,itile)
    END DO ! itile

    ! scale cover_fract to fractions relative to the veg_tile
!XXX modify to account for bare soil tile    ! XXX
!     cover_fract(:,:) = cover_fract(:,:) / MAX(SPREAD(cover_fract_sum(:), NCOPIES=no_children, DIM=2),1.e-10_wp)
!!$    WRITE (*,*) 'sum of cover_fract: ',SUM(cover_fract(:,:),DIM=2)

    ! allocate more local variables
    ALLOCATE (cover_fract_pot_previous_year(nc,no_children))
    ALLOCATE (fuel(nc))

    ! initialize act_fpc, bare fpc, and pot_fpc
    IF (lstart) THEN
       act_fpc(:,:) = cover_fract(:,:)
       bare_fpc(:) = 0._wp
       CALL scale_fpc(routine, nc, iblk, no_children, is_dynamic, act_fpc, bare_fpc)
       pot_fpc(:,:) = act_fpc(:,:) !!$ TR to be improved (scale, non-woodies = 0, etc.)
    ENDIF

    ! calculate climate buffer
    CALL calc_q_climbuf(lstart, init_running_means, new_day, new_month, new_year, dtime, t_air, &
                      seconds_day, seconds_month, temp_sum_day, temp_sum_month,                 &
                      min_mmtemp_of_yr, max_mmtemp_of_yr, min_mmtemp20, max_mmtemp20,           &
                      gdd_sum_year, gdd_prev_year)

    ! annual calculations
    IF (new_year) THEN
       DO itile = 1,no_children
          CALL bioclim_q_limits (is_dynamic(:,itile), bclimit_min_cold_mmtemp(:,itile), bclimit_max_cold_mmtemp(:,itile),  &
                               bclimit_max_warm_mmtemp(:,itile), bclimit_min_temprange(:,itile), bclimit_min_gdd(:,itile), &
                               min_mmtemp20(:), max_mmtemp20(:), gdd_prev_year(:), bio_exist(:,itile))
       ENDDO
       IF (.NOT. lstart) THEN
          CALL potential_tree_fpc (nc, no_children, is_dynamic(:,:), is_woody(:,:), bio_exist(:,:), &
                                   NPP_mean_5year(:,:), act_fpc(:,:), pot_fpc(:,:))
       ENDIF

       cover_fract_pot_previous_year(:,:) = cover_fract_pot(:,:)
       CALL fpc_to_cover_fract_pot(nc, no_children, is_dynamic(:,:), is_woody(:,:), &
                                   act_fpc(:,:), bare_fpc(:), cover_fract_pot(:,:))

       IF (.NOT. lstart) THEN
          ! new desert fraction
          CALL desert_fraction(nc, no_children, init_running_means, accelerate_qlcc,   &
                               is_woody(:,:), is_dynamic(:,:), is_bare(:,:),           &
                               act_fpc(:,:), bare_fpc(:),                              &
                               sla(:,:), max_green_bio(:,:),                           &
                               sum_green_bio_memory(:), desert_fpc(:))

          act_fpc(:, bare_tile_no) = desert_fpc(:)

          ! new cover fractions
          CALL cover_fract_pot_to_cover_fract(nc, no_children, is_dynamic(:,:), is_woody(:,:), is_pasture(:,:), is_crop(:,:), &
                                              is_bare(:,:), cover_fract_pot(:,:), cover_fract_pot_previous_year(:,:), cover_fract(:,:))

          DO itile = 1,no_children
             ! Point to the correct pft
             IF (itile == 1) THEN
                ptr_tile => tile%Get_first_child_tile()
             ELSE
                ptr_tile => ptr_tile%Get_next_sibling_tile()
             END IF
             ! set cover fractions
             CALL ptr_tile%Set_fraction(ics, ice, iblk, fract=cover_fract(:,itile))
          END DO ! itile

       ENDIF ! not lstart

    ENDIF ! new_year

    ! PPLCC fraction diagnostics are calculated on the box tile, although the PPLCC process is calculated
    ! on leafs (compare usecase definition). We only enter this routine on the box tile, because
    ! update_pplcc_post_lcc is not only defined as Integrate but also as Aggregate procedure and thus is
    ! called from the box tile on the way up the tile tree.
    IF (.NOT. ASSOCIATED(tile%parent)) CALL qlcc_fraction_diagnostics(tile, options)

    ! daily calculations
    IF (new_day) THEN
       CALL fpc_daily(nc, no_children, accelerate_qlcc, is_woody(:,:), is_dynamic(:,:), is_bare(:,:), tau_c_woods(:,:),       &
                      NPP_mean_5year(:,:), bio_exist(:,:), burned_fract(:,:), damaged_fract(:,:), pot_fpc(:,:), &
                      act_fpc(:,:), bare_fpc(:), fuel(:))
    ENDIF

    DEALLOCATE(is_woody, is_grass, is_pasture, is_crop, is_dynamic, is_bare)
    DEALLOCATE(bclimit_min_cold_mmtemp, bclimit_max_cold_mmtemp, &
               bclimit_max_warm_mmtemp, bclimit_min_temprange, bclimit_min_gdd, tau_c_woods)
!     IF (new_year .AND. .NOT. lstart) DEALLOCATE(fract_fpc_max)
    DEALLOCATE(NPP_mean_5year)
    DEALLOCATE(burned_fract, damaged_fract)
    DEALLOCATE(max_green_bio, sla)
    DEALLOCATE(fuel)
    DEALLOCATE(cover_fract_pot_previous_year)
!     DEALLOCATE(cover_fract_sum)
    DEALLOCATE(var_collect)


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_qlcc

  !-----------------------------------------------------------------------------------------------------
  !> Post LCC diagnostics of land cover fractions, i.e. aggregation of different land cover types
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE qlcc_fraction_diagnostics(tile, options)

    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in) :: options  !< Additional run-time parameters

    ! Local variables
    TYPE(t_jsb_model), POINTER           :: model
    CLASS(t_jsb_tile_abstract),  POINTER :: ptr_tile, ptr_pft
    REAL(wp), DIMENSION(options%nc)      :: pft_fract, veg_fract, bare_fract
    INTEGER                              :: iblk, ics, ice
    INTEGER                              :: ilct, nlct
    INTEGER                              :: ipft, npft
    REAL(wp)                             :: fract_fpc_max

    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(QLCC_)
!     dsl4jsb_Def_memory(PHENO_)

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':qlcc_fraction_diagnostics'

    ! Pointers to variables in memory
    !
!     dsl4jsb_Real2D_onChunk :: fract_fpc_max

    dsl4jsb_Real2D_onChunk :: tree_fract
    dsl4jsb_Real2D_onChunk :: shrub_fract
    dsl4jsb_Real2D_onChunk :: grass_fract
    dsl4jsb_Real2D_onChunk :: crop_fract
    dsl4jsb_Real2D_onChunk :: pasture_fract
    dsl4jsb_Real2D_onChunk :: baresoil_fract
    dsl4jsb_Real2D_onChunk :: C3pft_fract
    dsl4jsb_Real2D_onChunk :: C4pft_fract


    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    ! TK: Hack to minimise code changes: fract_fpc_max = 1 everywhere.
    ! Should really be removed.
    fract_fpc_max = 1._wp

    dsl4jsb_Get_memory(QLCC_)
!     dsl4jsb_Get_memory(PHENO_)

    model => Get_model(tile%owner_model_id)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
!     CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! The fraction variables are defined as fractions relative to the grid cell fraction. Thus
    ! only calculations on the box tile make sense.
    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    ! Get required variables on this tile
!     dsl4jsb_Get_var2D_onChunk(PHENO_, fract_fpc_max)        ! in

    dsl4jsb_Get_var2D_onChunk(QLCC_, tree_fract)           ! out
    dsl4jsb_Get_var2D_onChunk(QLCC_, shrub_fract)          ! out
    dsl4jsb_Get_var2D_onChunk(QLCC_, grass_fract)          ! out
    dsl4jsb_Get_var2D_onChunk(QLCC_, crop_fract)           ! out
    dsl4jsb_Get_var2D_onChunk(QLCC_, pasture_fract)        ! out
    dsl4jsb_Get_var2D_onChunk(QLCC_, baresoil_fract)       ! out
    dsl4jsb_Get_var2D_onChunk(QLCC_, C3pft_fract)          ! out
    dsl4jsb_Get_var2D_onChunk(QLCC_, C4pft_fract)          ! out

    ! Initialization
    tree_fract     = 0._wp
    shrub_fract    = 0._wp
    grass_fract    = 0._wp
    crop_fract     = 0._wp
    pasture_fract  = 0._wp
    baresoil_fract = 0._wp
    C3pft_fract    = 0._wp
    C4pft_fract    = 0._wp


    ! To get the PFT fractions we first need a pointer to the land tile,
    ! and from the land tile to the veg tile, which has PFT children

    ! Get pointer to the land tile
    nlct = tile%Get_no_of_children()
    DO ilct=1,nlct
      IF (ilct == 1) THEN
        ptr_tile => tile%Get_first_child_tile()
      ELSE
        ptr_tile => ptr_tile%Get_next_sibling_tile()
      END IF
      ! Only the land tile needs to be considered
      IF (TRIM(ptr_tile%name) == 'land') EXIT
    END DO

    ! Get pointer to the veg tile
    nlct = ptr_tile%Get_no_of_children()
    DO ilct=1,nlct
      IF (ilct == 1) THEN
        ptr_tile => ptr_tile%Get_first_child_tile()
      ELSE
        ptr_tile => ptr_tile%Get_next_sibling_tile()
      END IF
      ! Only the veg tile needs to be considered
      IF (TRIM(ptr_tile%name) == 'veg') EXIT
    END DO

    npft = ptr_tile%Get_no_of_children()
    DO ipft=1,npft

      ! Point to the correct pft
      IF (ipft == 1) THEN
        ptr_pft => ptr_tile%Get_first_child_tile()
      ELSE
        ptr_pft => ptr_pft%Get_next_sibling_tile()
      END IF

      ! Get fraction for this PFT
      CALL ptr_pft%Get_fraction(ics, ice, iblk, fract=pft_fract(:))

      IF (model%lctlib(ptr_pft%lcts(1)%lib_id)%NaturalVegFlag) THEN
        IF (model%lctlib(ptr_pft%lcts(1)%lib_id)%ForestFlag) THEN
          tree_fract(:) = tree_fract(:) + pft_fract(:) * fract_fpc_max
        ELSE IF (model%lctlib(ptr_pft%lcts(1)%lib_id)%GrassFlag) THEN
          grass_fract(:) = grass_fract(:) + pft_fract(:) * fract_fpc_max
        ELSE   ! Current lctlib file does not include a shrub flag
          shrub_fract(:) = shrub_fract(:) + pft_fract(:) * fract_fpc_max
        END IF
      ELSE IF (model%lctlib(ptr_pft%lcts(1)%lib_id)%CropFlag) THEN
        crop_fract(:) = crop_fract(:) + pft_fract(:) * fract_fpc_max
      ELSE IF (model%lctlib(ptr_pft%lcts(1)%lib_id)%PastureFlag) THEN
        pasture_fract(:) = pasture_fract(:) + pft_fract(:) * fract_fpc_max
      END IF
      IF (model%lctlib(ptr_pft%lcts(1)%lib_id)%BareSoilflag) THEN
        baresoil_fract(:) = baresoil_fract(:) + pft_fract(:) * fract_fpc_max
      ELSE
        IF (model%lctlib(ptr_pft%lcts(1)%lib_id)%ps_pathway == 2) THEN
          C4pft_fract(:) = C4pft_fract(:) + pft_fract(:) * fract_fpc_max
        ELSE
          C3pft_fract(:) = C3pft_fract(:) + pft_fract(:) * fract_fpc_max
        END IF
      END IF

    END DO ! ipft

!     ! Fraction of the veg tile
!     CALL ptr_tile%Get_fraction(ics, ice, iblk, fract=veg_fract(:))

!     ! Bare land tile
!     ! Note: In the current usecase the land tile does not have a bare soil child, thus
!     !       bare_fract is set to zero.
!     bare_fract(:)=0._wp
!
!     ! Baresoil fraction
!     baresoil_fract(:) = bare_fract(:) + veg_fract(:) * (1._wp-fract_fpc_max(:))


  END SUBROUTINE qlcc_fraction_diagnostics

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "qlcc"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_qlcc(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(QLCC_)

    ! Local variables
    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract
    INTEGER  :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_qlcc'
    REAL(wp) :: fract_fpc_max

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
!     CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! PPLCC fraction diagnostics are calculated on the box tile, although the PPLCC process is calculated
    ! on leafs (compare usecase definition). We only enter this routine on the box tile, because
    ! update_pplcc_post_lcc is not only defined as Integrate but also as Aggregate procedure and thus is
    ! called from the box tile on the way up the tile tree.
    IF (.NOT. ASSOCIATED(tile%parent)) CALL qlcc_fraction_diagnostics(tile, options)


    ! -------------------------------------------------------------------------------
    ! originally in mo_qlcc_process...
    ! copied from pplcc
    ! SUBROUTINE qlcc_fraction_diagnostics

    dsl4jsb_Get_memory(QLCC_)

    fract_fpc_max = 1._wp

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(QLCC_, min_mmtemp20,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(QLCC_, max_mmtemp20,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(QLCC_, gdd_prev_year,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(QLCC_, bare_fpc,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(QLCC_, desert_fpc,            weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(QLCC_, act_fpc,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(QLCC_, pot_fpc,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(QLCC_, cover_fract_pot,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(QLCC_, cover_fract,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(QLCC_, bio_exist,             weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_qlcc

! #endif
END MODULE mo_qlcc_interface
