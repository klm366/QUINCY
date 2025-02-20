!> QUINCY simple harvest process interface
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### definition and init of tasks for the quincy simple harvest process
!>
MODULE mo_q_sh_interface
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_math_constants,  ONLY: eps8
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_jsb_time,            ONLY: is_newyear, is_newday, get_year, get_year_length
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,   ONLY: t_jsb_process
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes Q_SH_, VEG_
  dsl4jsb_Use_config(Q_SH_)
  dsl4jsb_Use_config(VEG_)
  dsl4jsb_Use_memory(Q_SH_)

#ifndef __QUINCY_STANDALONE__
  dsl4jsb_Use_processes PPLCC_
  dsl4jsb_Use_memory(PPLCC_)
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_q_sh_tasks
  PUBLIC ::  update_simple_harvest

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_sh_interface'

  ! ======================================================================================================= !
  !> Type definition: simple_harvest task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_simple_harvest
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_simple_harvest
    PROCEDURE, NOPASS :: Aggregate => aggregate_simple_harvest
  END TYPE tsk_simple_harvest

  !> Constructor interface: update_simple_harvest task
  !>
  INTERFACE tsk_simple_harvest
    PROCEDURE Create_task_simple_harvest
  END INTERFACE tsk_simple_harvest

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: Q_SH_
  !>
  SUBROUTINE Register_q_sh_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_simple_harvest(model_id))

  END SUBROUTINE Register_q_sh_tasks

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_simple_harvest task
  !>
  FUNCTION Create_task_simple_harvest(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(tsk_simple_harvest::return_ptr)
    CALL return_ptr%Construct(name='simple_harvest', process_id=Q_SH_, owner_model_id=model_id)
  END FUNCTION Create_task_simple_harvest

  ! ======================================================================================================= !
  !> Implementation of "update": simple_harvest task
  !>
  SUBROUTINE update_simple_harvest(tile, options)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), POINTER :: box_tile
    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(Q_SH_)

#ifndef __QUINCY_STANDALONE__
    dsl4jsb_Def_memory_tile(PPLCC_, box_tile)
    dsl4jsb_Real2D_onChunk :: tree_fract_box_tile
#endif

    dsl4jsb_Def_memory(Q_SH_)
    dsl4jsb_Real2D_onChunk :: fract_harvest   ! fraction harvested in this timestep (TODO: unit / rel to?)
    dsl4jsb_Real2D_onChunk :: fract_wood_to_slash
    dsl4jsb_Real2D_onChunk :: fract_wood_to_pp_fuel
    dsl4jsb_Real2D_onChunk :: fract_wood_to_pp_paper
    dsl4jsb_Real2D_onChunk :: fract_wood_to_pp_fiberboard
    dsl4jsb_Real2D_onChunk :: fract_wood_to_pp_oirw
    dsl4jsb_Real2D_onChunk :: fract_wood_to_pp_pv
    dsl4jsb_Real2D_onChunk :: fract_wood_to_pp_sawnwood
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                         :: ic, iblk, ics, ice, nc !< loop counter / dimensions
    REAL(wp)                        :: dtime                  !< timestep length
    INTEGER                         :: number_of_days         !< number of days in this year
    INTEGER                         :: current_year           !< the current year
    REAL(wp)                        :: fract_sum              !< sum of allocation fractions to different pools
    REAL(wp), DIMENSION(options%nc) :: cover_fraction         !< cover fraction of given tile
    REAL(wp), DIMENSION(options%nc) :: fact_harvest_share     !< factor to determine the harvest share of this tile
                                                              !< depends on cover fraction and total area with PFTs with trees
    REAL(wp), DIMENSION(options%nc) :: box_fract_harvest      !< total area relative to grid cell harvest from all forested PFTS
    TYPE(t_jsb_model), POINTER      :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_simple_harvest'

    ! ----------------------------------------------------------------------------------------------------- !
    ! Get pointers to process configs and memory
    model => Get_model(tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(Q_SH_)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_SH_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    nc    = options%nc
    ics   = options%ics
    ice   = options%ice
    iblk  = options%iblk
    dtime = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_SH_)
    dsl4jsb_Get_var2D_onChunk(Q_SH_, fract_harvest)               ! out
    IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
      dsl4jsb_Get_var2D_onChunk(Q_SH_, fract_wood_to_slash)         ! in
      dsl4jsb_Get_var2D_onChunk(Q_SH_, fract_wood_to_pp_fuel)       ! in
      dsl4jsb_Get_var2D_onChunk(Q_SH_, fract_wood_to_pp_paper)      ! in
      dsl4jsb_Get_var2D_onChunk(Q_SH_, fract_wood_to_pp_fiberboard) ! in
      dsl4jsb_Get_var2D_onChunk(Q_SH_, fract_wood_to_pp_oirw)       ! in
      dsl4jsb_Get_var2D_onChunk(Q_SH_, fract_wood_to_pp_pv)         ! in
      dsl4jsb_Get_var2D_onChunk(Q_SH_, fract_wood_to_pp_sawnwood)   ! in
    END IF
    ! ----------------------------------------------------------------------------------------------------- !
    ! Since harvest only is happening daily or annually (or in QS even just once at the stand replacing year)
    ! the default would be that nothing is harvested in this timestep
    fract_harvest(:) = 0._wp

    ! TODO: change this later to testing if running quincy standalone (and of course change the message accordingly...)
    IF (dsl4jsb_Config(Q_SH_)%l_read_harvest_map) THEN
      CALL finish( routine, 'Reading of harvest maps is not yet implemented. ')
    END IF

    ! Check if some harvest is to be done or if we can RETURN here
    IF (dsl4jsb_Config(Q_SH_)%flag_stand_harvest) THEN
      IF (dsl4jsb_Config(Q_SH_)%flag_stand_harvest_event) THEN
        fract_harvest(:) =  1._wp
      ELSE
        ! If we run QS and use stand replacing harvest, but do not have a harvest event, then we can return
      END IF
    ELSE IF (dsl4jsb_Config(Q_SH_)%l_daily_harvest) THEN
      ! In case of daily harvest we have to check if this is a new day
      IF (.NOT. is_newday(options%current_datetime, dtime)) THEN
        ! If not then we do not have a harvest event and can return
        RETURN
      END IF
    ELSE
      ! In case of annual harvest we have to check if this is a new year
      IF (.NOT. is_newyear(options%current_datetime, dtime)) THEN
        ! If not then we do not have a harvest event and can return
        RETURN
      END IF
    END IF

    ! If this is not a stand replacing harvest event in which the whole tile is to be harvested,
    ! the to be harvested area still needs to be determined
    IF (.NOT. dsl4jsb_Config(Q_SH_)%flag_stand_harvest) THEN

      IF (dsl4jsb_Config(Q_SH_)%l_read_harvest_map) THEN
        ! TODO: get read maps when implemented
      ELSE
        ! If not maps have been read, the to be harvested area is derived from the constant value from the config
        box_fract_harvest(:) = dsl4jsb_Config(Q_SH_)%harvest_fraction
      END IF

      ! Calculate the factor with which to derive the share to be harvested from this tile
      CALL tile%Get_fraction(ics, ice, iblk, fract=cover_fraction)
      fact_harvest_share(:) = cover_fraction(:)

#ifndef __QUINCY_STANDALONE__
      CALL model%Get_top_tile(box_tile)
      dsl4jsb_Get_memory_tile(PPLCC_, box_tile)
      dsl4jsb_Get_var2d_onChunk_tile_name(PPLCC_, tree_fract, box_tile) ! in

      DO ic = 1,nc
        IF (box_fract_harvest(ic) > tree_fract_box_tile(ic)) THEN
          IF (debug_on()) THEN
            WRITE(message_text,*) 'The harvested fraction (', box_fract_harvest(ic), &
              &                    ') exceeds the tree area, harvesting available area only...'
            CALL message(TRIM(routine),message_text)
          END IF
          box_fract_harvest(ic) = tree_fract_box_tile(ic) - eps8
        END IF

        IF (tree_fract_box_tile(ic) > eps8) THEN
          fact_harvest_share(ic) = cover_fraction(ic) / tree_fract_box_tile(ic)
        ELSE
          fact_harvest_share(ic) = 0.0_wp ! TODO: think about this
        END IF
      END DO
#endif

      IF (dsl4jsb_Config(Q_SH_)%l_daily_harvest) THEN
        ! In case of daily harvest the harvest needs to be distributed over the number of days in this year
        current_year  = get_year(options%current_datetime)
        number_of_days = get_year_length(current_year)
        fact_harvest_share(:) = fact_harvest_share(:) / REAL(number_of_days, wp)
      END IF

      ! If not maps have been read, the to be harvested area is derived from the constant value from the config
      fract_harvest(:) = box_fract_harvest(:) * fact_harvest_share(:)

      ! Assertions
      DO ic = 1,nc
        ! Assert: harvest fraction should not be negative
        IF (fract_harvest(ic) < 0._wp) THEN
          WRITE(message_text,*) 'The harvested fraction (', fract_harvest(ic), &
            &                    ') should not be negative. Please check!'
          CALL finish(TRIM(routine),message_text)
        END IF

        IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
          ! Assert: sum of allocation fractions to different product pools should sum up to 1._
          !         !NOTE: i.e. not including the slash fraction!
          fract_sum = fract_wood_to_pp_fuel(ic) + fract_wood_to_pp_paper(ic)        &
            &         + fract_wood_to_pp_fiberboard(ic) + fract_wood_to_pp_oirw(ic) + fract_wood_to_pp_pv(ic) &
            &         + fract_wood_to_pp_sawnwood(ic)
          IF (ABS(fract_sum - 1._wp) >= eps8) THEN
            WRITE(message_text,*) 'A sum of allocation fractions to different pools does not sum up to one, but was: ', &
              &                    fract_sum,'! Please check!'
            CALL finish(TRIM(routine),message_text)
          END IF
        END IF
      END DO

      ! Assert: harvest fraction is not allowed to exceed the available cover fraction
      DO ic = 1,nc
        IF (fract_harvest(ic) > cover_fraction(ic)) THEN
          IF (fract_harvest(ic) - eps8 < cover_fraction(ic)) THEN
            fract_harvest(ic) = cover_fraction(ic)
          ELSE
            WRITE(message_text,*) 'The harvested fraction (', fract_harvest(ic), &
              & ') is not allowed to exceed the actually available fraction in this tile (', &
              & cover_fraction(ic),' in ', TRIM(tile%name),'). Please check!'
            CALL finish(TRIM(routine),message_text)
          END IF
        END IF
      END DO
    END IF ! .NOT. dsl4jsb_Config(Q_SH_)%flag_stand_harvest

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_simple_harvest

  ! ======================================================================================================= !
  !> Implementation of "aggregate": simple_harvest task
  !>
  SUBROUTINE aggregate_simple_harvest(tile, options)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_SH_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_simple_harvest'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_SH_)
    dsl4jsb_Aggregate_onChunk(Q_SH_, fract_harvest, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE aggregate_simple_harvest

#endif
END MODULE mo_q_sh_interface
