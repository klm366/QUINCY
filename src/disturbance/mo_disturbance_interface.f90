!>
!! @brief Contains the interfaces to the disturb process
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-05-04)
!!
MODULE mo_disturb_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process, DISTURB_
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module (Get_disturb_memory and Get_disturb_config)
  dsl4jsb_Use_processes DISTURB_, ASSIMI_, HYDRO_, PHENO_, RAD_, A2L_, CARBON_, SEB_

  ! Use of process configurations (t_disturb_config)
  dsl4jsb_Use_config(ASSIMI_)
  dsl4jsb_Use_config(DISTURB_)

  ! Use of process memories (t_disturb_memory)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(DISTURB_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(SEB_)

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

    return_ptr%name = 'natural_disturbances'
    return_ptr%process_id = DISTURB_
    return_ptr%owner_model_id = model_id

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

    USE mo_disturb_process,       ONLY: disturbed_frac, get_relative_humidity_air
    USE mo_carbon_interface,      ONLY: relocate_disturbed_carbon
    USE mo_disturbance_constants, ONLY: DIST_FIRE, DIST_WINDBREAK
    USE mo_jsb_time,              ONLY: is_newday
    USE mo_jsb_tile_class,        ONLY: t_jsb_aggregator_weighted_by_fract

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER                :: model

    LOGICAL, allocatable                      :: pft_mask(:), &            ! Dimension: (no_children).
                                                                           ! True for pfts to which the disturbance should be applied.
                                                 is_woody(:)               ! Dimension: no_children. True for woody pfts.

    INTEGER, allocatable                      :: no_children, &            ! Number of existing sub tiles
                                                 i, j                      ! Index numbers over children of the current tile

    ! R: The DISTURB_ process is calculated on the veg tile, because fire needs "fuel", which is a sum over the child tiles!
    !    Furthermore, the climate-buffers can be handled somewhat easier this way, as they have to be called only once...
    !    However, windbreak as a process could have been calculated also on the child tiles.
    !
    ! R: These variables are needed for our calculations on the parent (=veg) tile, but we do not want them in the
    !    memory of the parent as a collection of its child tiles. We want them on each child tile! Therefore will give
    !    these parent results back to the corresponding childs at the end of this procedure.
    ! The following variables all have the dimensions (ics:ice , no_children). ics:ice values do not change.
    ! They all collect the corresponding variables of the child pfts
    REAL(wp), allocatable                     :: veg_frac_correction_st(:,:),   &
                                                 fract_st(:,:),                 &
                                                 cover_frac_pot_st(:,:),        &
                                                 burned_frac_st(:,:),           &
                                                 damaged_frac_st(:,:),          &
                                                 !
                                                 Cpool_green_st(:,:),           &
                                                 Cpool_woods_st(:,:),           &
                                                 Cpool_reserve_st(:,:),         &
                                                 carbon_2_GreenLitterPools_st(:,:), &
                                                 carbon_2_WoodLitterPools_st(:,:),  &
                                                 !
                                                 YCpool_acid_ag1_st(:,:),       &
                                                 YCpool_water_ag1_st(:,:),      &
                                                 YCpool_ethanol_ag1_st(:,:),    &
                                                 YCpool_nonsoluble_ag1_st(:,:), &
                                                 !
                                                 YCpool_acid_ag2_st(:,:),       &
                                                 YCpool_water_ag2_st(:,:),      &
                                                 YCpool_ethanol_ag2_st(:,:),    &
                                                 YCpool_nonsoluble_ag2_st(:,:), &
                                                 !
                                                 YCpool_humus_1_st(:,:),        &
                                                 YCpool_humus_2_st(:,:),        &
                                                 carbon_2_atmos_st_ta(:,:),     &
                                                 !
                                                 YCpool_acid_bg1_st(:,:),       &
                                                 YCpool_water_bg1_st(:,:),      &
                                                 YCpool_ethanol_bg1_st(:,:),    &
                                                 YCpool_nonsoluble_bg1_st(:,:), &
                                                 YCpool_acid_bg2_st(:,:),       &
                                                 YCpool_water_bg2_st(:,:),      &
                                                 YCpool_ethanol_bg2_st(:,:),    &
                                                 YCpool_nonsoluble_bg2_st(:,:)


    INTEGER, allocatable                      :: iblk, ics, ice, nc

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_natural_disturbances'

    ! Declare process memories
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(DISTURB_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_config(DISTURB_)

    ! Declare pointers to variables in memory
    !INTEGER :: ncanopy                 ! number of canopy layers

    !dsl4jsb_Real2D_onChunk :: C4flag ! Photosynthetic pathway (C3: 0; C4: 1)
                                      ! R: nur falls man C4flag doch als REAL haben wollte


!     dsl4jsb_Real2D_onChunk :: CO2_air_mol           ! CO2 particle mixing ratio [molCO2/molDryAir]

    dsl4jsb_Real2D_onChunk :: fuel                   ! Amount of fuel on which the fire was based.
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: press_srf
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: q_rel_air
    dsl4jsb_Real2D_onChunk :: wind_10m               ! 10m wind speed [m/s]
    dsl4jsb_Real2D_onChunk :: max_wind_10m_act
    dsl4jsb_Real2D_onChunk :: prev_day_max_wind_10m
    dsl4jsb_Real2D_onChunk :: max_wind_10m
    dsl4jsb_Real2D_onChunk :: q_rel_air_climbuf
    dsl4jsb_Real2D_onChunk :: q_rel_air_climbuf_yDay

    dsl4jsb_Real2D_onChunk :: burned_frac
    dsl4jsb_Real2D_onChunk :: damaged_frac
    !
    dsl4jsb_Real2D_onChunk :: Cpool_green
    dsl4jsb_Real2D_onChunk :: Cpool_woods
    dsl4jsb_Real2D_onChunk :: Cpool_reserve
    dsl4jsb_Real2D_onChunk :: carbon_2_GreenLitterPools
    dsl4jsb_Real2D_onChunk :: carbon_2_WoodLitterPools
    !
    dsl4jsb_Real2D_onChunk :: YCpool_acid_ag1
    dsl4jsb_Real2D_onChunk :: YCpool_water_ag1
    dsl4jsb_Real2D_onChunk :: YCpool_ethanol_ag1
    dsl4jsb_Real2D_onChunk :: YCpool_nonsoluble_ag1
    dsl4jsb_Real2D_onChunk :: YCpool_acid_ag2
    dsl4jsb_Real2D_onChunk :: YCpool_water_ag2
    dsl4jsb_Real2D_onChunk :: YCpool_ethanol_ag2
    dsl4jsb_Real2D_onChunk :: YCpool_nonsoluble_ag2
    ! below ground
    dsl4jsb_Real2D_onChunk :: YCpool_acid_bg1
    dsl4jsb_Real2D_onChunk :: YCpool_water_bg1
    dsl4jsb_Real2D_onChunk :: YCpool_ethanol_bg1
    dsl4jsb_Real2D_onChunk :: YCpool_nonsoluble_bg1
    dsl4jsb_Real2D_onChunk :: YCpool_acid_bg2
    dsl4jsb_Real2D_onChunk :: YCpool_water_bg2
    dsl4jsb_Real2D_onChunk :: YCpool_ethanol_bg2
    dsl4jsb_Real2D_onChunk :: YCpool_nonsoluble_bg2
    !
    dsl4jsb_Real2D_onChunk :: YCpool_humus_1
    dsl4jsb_Real2D_onChunk :: YCpool_humus_2
    dsl4jsb_Real2D_onChunk :: carbon_2_atmos
    !
    dsl4jsb_Real2D_onChunk :: veg_frac_correction
    !
    dsl4jsb_Real2D_onChunk :: Cpool_green_ta
    dsl4jsb_Real2D_onChunk :: Cpool_woods_ta
    dsl4jsb_Real2D_onChunk :: Cpool_reserve_ta
    dsl4jsb_Real2D_onChunk :: carbon_2_GreenLitterPools_ta
    dsl4jsb_Real2D_onChunk :: carbon_2_WoodLitterPools_ta
    !
    dsl4jsb_Real2D_onChunk :: YCpool_acid_ag1_ta
    dsl4jsb_Real2D_onChunk :: YCpool_water_ag1_ta
    dsl4jsb_Real2D_onChunk :: YCpool_ethanol_ag1_ta
    dsl4jsb_Real2D_onChunk :: YCpool_nonsoluble_ag1_ta
    dsl4jsb_Real2D_onChunk :: YCpool_acid_ag2_ta
    dsl4jsb_Real2D_onChunk :: YCpool_water_ag2_ta
    dsl4jsb_Real2D_onChunk :: YCpool_ethanol_ag2_ta
    dsl4jsb_Real2D_onChunk :: YCpool_nonsoluble_ag2_ta
    ! below ground
    dsl4jsb_Real2D_onChunk :: YCpool_acid_bg1_ta
    dsl4jsb_Real2D_onChunk :: YCpool_water_bg1_ta
    dsl4jsb_Real2D_onChunk :: YCpool_ethanol_bg1_ta
    dsl4jsb_Real2D_onChunk :: YCpool_nonsoluble_bg1_ta
    dsl4jsb_Real2D_onChunk :: YCpool_acid_bg2_ta
    dsl4jsb_Real2D_onChunk :: YCpool_water_bg2_ta
    dsl4jsb_Real2D_onChunk :: YCpool_ethanol_bg2_ta
    dsl4jsb_Real2D_onChunk :: YCpool_nonsoluble_bg2_ta
    !
    dsl4jsb_Real2D_onChunk :: YCpool_humus_1_ta
    dsl4jsb_Real2D_onChunk :: YCpool_humus_2_ta
    dsl4jsb_Real2D_onChunk :: carbon_2_atmos_ta


    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(ASSIMI_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get model
    model   => Get_model(tile%owner_model_id)

    ! Get process configurations
    dsl4jsb_Get_config(DISTURB_)

    ! Get process memories
    dsl4jsb_Get_memory(DISTURB_)
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SEB_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   fuel)                     ! in, note: only on parent tile
    dsl4jsb_Get_var2D_onChunk(A2L_,       q_air)                    ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,       press_srf)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,       wind_10m)                 ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,       t_unfilt)                 ! in
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   max_wind_10m_act)         ! inout
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   prev_day_max_wind_10m)    ! inout
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   max_wind_10m)             ! inout
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   q_rel_air_climbuf)        ! inout
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   q_rel_air_climbuf_yDay)   ! out

    dsl4jsb_Get_var2D_onChunk(DISTURB_,   burned_frac)              ! out
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   damaged_frac)             ! out
    !
    dsl4jsb_Get_var2D_onChunk(CARBON_,    Cpool_green)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    Cpool_woods)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    Cpool_reserve)            ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    carbon_2_GreenLitterPools)! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    carbon_2_WoodLitterPools) ! inout
    !
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_acid_ag1)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_water_ag1)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_ethanol_ag1)       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_nonsoluble_ag1)    ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_acid_ag2)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_water_ag2)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_ethanol_ag2)       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_nonsoluble_ag2)    ! inout
    ! below ground
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_acid_bg1)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_water_bg1)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_ethanol_bg1)       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_nonsoluble_bg1)    ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_acid_bg2)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_water_bg2)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_ethanol_bg2)       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_nonsoluble_bg2)    ! inout
    !
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_humus_1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_humus_2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    carbon_2_atmos)           ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_,     veg_frac_correction )     ! in
    !
    dsl4jsb_Get_var2D_onChunk(CARBON_,    Cpool_green_ta)              ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    Cpool_woods_ta)              ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    Cpool_reserve_ta)            ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    carbon_2_GreenLitterPools_ta)! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    carbon_2_WoodLitterPools_ta) ! out
    !
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_acid_ag1_ta)          ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_water_ag1_ta)         ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_ethanol_ag1_ta)       ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_nonsoluble_ag1_ta)    ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_acid_ag2_ta)          ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_water_ag2_ta)         ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_ethanol_ag2_ta)       ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_nonsoluble_ag2_ta)    ! out
    ! below ground
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_acid_bg1_ta)          ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_water_bg1_ta)         ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_ethanol_bg1_ta)       ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_nonsoluble_bg1_ta)    ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_acid_bg2_ta)          ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_water_bg2_ta)         ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_ethanol_bg2_ta)       ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_nonsoluble_bg2_ta)    ! out
    !
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_humus_1_ta)           ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    YCpool_humus_2_ta)           ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    carbon_2_atmos_ta)           ! out

    ! ---------------------------
    ! Go

    IF (is_newday()) THEN ! only once per day
      burned_frac  = 0._wp
      damaged_frac = 0._wp

     ! JS3:  CALL average_tiles(soil%surface_temperature_unfiltered(kidx0:kidx1,:), &
     !                               surface%is_present(kidx0:kidx1,:) .AND. .NOT. surface%is_lake(kidx0:kidx1,:), &
     !                               surface%cover_fract(kidx0:kidx1,:),        &
     !                               surface_temperature_unfilt_avg)
     ! t_unfilt_avgtiles= Average_over_children_real2d(....) ! R: Not necessary because the DISTURB_ process runs on the veg tile
                                                             !    and should inherit its t_unfilt from the SEB process running
                                                             !    on land! => t_unfilt_avgtiles =t_unfilt
      q_rel_air(:) = get_relative_humidity_air(nc, q_air(:), t_unfilt(:), press_srf(:))

      ! R: The lctlib parameters have to be collected for all child pfts. This function is not present in dsl4jsb in the moment
      !    (dsl4jsb_Lctlib_param(dynamic_pft) is related to the current tile). However, this is one possibility:
      no_children = tile%Get_no_of_children()
      DO i=1,no_children
         j=i+1  ! lcts(1) is the current parent tile, lcts(2) is the first child tile

         ! Gather woody pfts
         IF (model%lctlib(tile%lcts(j)%lib_id)%WOODY_PFT) THEN
           is_woody(i) = .TRUE.
         ELSE
           is_woody(i) = .FALSE.
         END IF

         ! Gather those pfts where disturbance should be applied
         ! R: Same pft_mask is used for fire and windbreak in JS3. We could change this in the future.
         IF ( model%lctlib(tile%lcts(j)%lib_id)%dynamic_pft .OR. &  ! 0 in lctlib means .FALSE.
             (dsl4jsb_Config(DISTURB_)%lburn_pasture .AND. model%lctlib(tile%lcts(j)%lib_id)%pasture_pft) ) THEN
           pft_mask(i) = .TRUE.
         ELSE
           pft_mask(i) = .FALSE.
         END IF
      END DO

      ! Get fractions of the childs
      SELECT TYPE (agg=>tile%aggregators(1)%p) ! First aggregator is in the moment aggregator_weighted_by_fract
          CLASS IS (t_jsb_aggregator_weighted_by_fract)
              fract_st(:,:) = agg%fractions(ics:ice,iblk,:) ! or also fract_st(:,:) = agg(ics:ice,iblk,:)
      END SELECT

      ! Collect variables from our child tiles
      SELECT TYPE (parentmem=>tile%mem(ASSIMI_)%p)
         CLASS IS (t_assimi_memory)
            cover_frac_pot_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children, parentmem%cover_frac_pot ) ! in
      END SELECT
      SELECT TYPE (parentmem=>tile%mem(PHENO_)%p)
         CLASS IS (t_pheno_memory)
            veg_frac_correction_st(:,:) = parentmem%Collect_var_real2d(ics, ice, iblk, no_children, parentmem%veg_frac_correction) ! in
      END SELECT
      SELECT TYPE (parentmem=>tile%mem(DISTURB_)%p)
         CLASS IS (t_disturb_memory)
            burned_frac_st(:,:) = parentmem%Collect_var_real2d(ics, ice, iblk, no_children,  parentmem%burned_frac ) ! inout
            damaged_frac_st(:,:) = parentmem%Collect_var_real2d(ics, ice, iblk, no_children,  parentmem%damaged_frac ) ! inout
      END SELECT

      IF ( tile%Is_process_active(CARBON_) ) THEN
        SELECT TYPE (parentmem=>tile%mem(CARBON_)%p)
         CLASS IS (t_carbon_memory)
            Cpool_green_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children,  parentmem%Cpool_green )     ! inout
            Cpool_woods_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children,  parentmem%Cpool_woods )     ! inout
            Cpool_reserve_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children,  parentmem%Cpool_reserve ) ! inout
            carbon_2_GreenLitterPools_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children,      &
                                                                              parentmem%carbon_2_GreenLitterPools )       ! inout
            carbon_2_WoodLitterPools_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children,       &
                                                                             parentmem%carbon_2_WoodLitterPools )         ! inout
            !
            YCpool_acid_ag1_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children, parentmem%YCpool_acid_ag1 )   ! inout
            YCpool_water_ag1_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children, parentmem%YCpool_water_ag1 ) ! inout
            YCpool_ethanol_ag1_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children,            &
                                                                       parentmem%YCpool_ethanol_ag1 )                     ! inout
            YCpool_nonsoluble_ag1_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children,         &
                                                                          parentmem%YCpool_nonsoluble_ag1 )               ! inout
            YCpool_acid_ag2_st(:,:) = parentmem%Collect_var_real2d( ics, ice,  iblk, no_children, parentmem%YCpool_acid_ag2 )  ! inout
            YCpool_water_ag2_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children, parentmem%YCpool_water_ag2 ) ! inout
            YCpool_ethanol_ag2_st(:,:) = parentmem%Collect_var_real2d( ics, ice,  iblk, no_children,            &
                                                                       parentmem%YCpool_ethanol_ag2 )                     ! inout
            YCpool_nonsoluble_ag2_st(:,:) = parentmem%Collect_var_real2d( ics, ice,  iblk, no_children,         &
                                                                          parentmem%YCpool_nonsoluble_ag2 )               ! inout
            YCpool_humus_1_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children, parentmem%YCpool_humus_1 ) ! inout
            YCpool_humus_2_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children, parentmem%YCpool_humus_2 ) ! inout
            !carbon_2_atmos_st_ta(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children, parentmem%carbon_2_atmos ) ! inout
            ! below ground
            YCpool_acid_bg1_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children, parentmem%YCpool_acid_bg1 )  ! inout
            YCpool_water_bg1_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children, parentmem%YCpool_water_bg1) ! inout
            YCpool_ethanol_bg1_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children,             &
                                                                       parentmem%YCpool_ethanol_bg1 )                          ! inout
            YCpool_nonsoluble_bg1_st(:,:) = parentmem%Collect_var_real2d( ics, ice,  iblk, no_children,         &
                                                                          parentmem%YCpool_nonsoluble_bg1 )                    ! inout
            YCpool_acid_bg2_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children,  parentmem%YCpool_acid_bg2)   ! inout
            YCpool_water_bg2_st(:,:) = parentmem%Collect_var_real2d( ics, ice, iblk, no_children,  parentmem%YCpool_water_bg2) ! inout
            YCpool_ethanol_bg2_st(:,:) = parentmem%Collect_var_real2d( ics, ice,  iblk, no_children,            &
                                                                       parentmem%YCpool_ethanol_bg2 )                          ! inout
            YCpool_nonsoluble_bg2_st(:,:) = parentmem%Collect_var_real2d( ics, ice,  iblk, no_children,         &
                                                                          parentmem%YCpool_nonsoluble_bg2 ) ! inout
        END SELECT
      ENDIF

      ! R: X1: Note: in the following we hand "tile" to the procedures disturbed_frac and to relocate_disturbed_carbon
      !    because we need the model configs. For disturbed_frac we also need the model options. Furthermore for
      !    relocate_disturbed_carbon we need it for the lctlib data of the subtiles. This should normally be done by
      !    handing down all these variables.
      !    However, in this case this would make it more difficult to find errors and later (post CMIP6) to change 
      !    the code. (Later there will be much more variables to hand down and the calling structure meight be 
      !    changed also.) Therefore for now it is better this way.

      CALL disturbed_frac(tile,                                 & ! in
                          !ics, ice, iblk,                      & ! in
                          DIST_FIRE,                            & ! in
                          pft_mask(:),                          & ! in
                          is_woody(:),                          & ! in
                          no_children,                          & ! in
                          fract_st(:,:),                        & ! in
                          cover_frac_pot_st(:,:),               & ! in
                          veg_frac_correction_st(:,:),          & ! in
                          q_rel_air(:),                         & ! in
                          wind_10m(:),                          & ! in
                          !tile%mem(CARBON_)%p ,                 & ! in R: only to find FUNCTION Average_over_children
                          max_wind_10m_act(:),                  & ! inout
                          prev_day_max_wind_10m(:),             & ! inout
                          max_wind_10m(:),                      & ! inout
                          q_rel_air_climbuf(:),                 & ! inout
                          q_rel_air_climbuf_yDay(:),            & ! out
                          burned_frac_st(:,:),                  & ! inout
                          !
                          YCpool_acid_ag1_st(:,:),              & ! in, optional , necessary for fire
                          YCpool_water_ag1_st(:,:),             & ! in, optional , necessary for fire
                          YCpool_ethanol_ag1_st(:,:),           & ! in, optional , necessary for fire
                          YCpool_nonsoluble_ag1_st(:,:),        & ! in, optional , necessary for fire
                          YCpool_acid_ag2_st(:,:),              & ! in, optional , necessary for fire
                          YCpool_water_ag2_st(:,:),             & ! in, optional , necessary for fire
                          YCpool_ethanol_ag2_st(:,:),           & ! in, optional , necessary for fire
                          YCpool_nonsoluble_ag2_st(:,:),        & ! in, optional , necessary for fire
                          fuel(:)                               & ! inout, optional , necessary for fire
                          )

      CALL disturbed_frac(tile,                                 & ! in
                         !ics, ice, iblk,                       & ! in
                          DIST_WINDBREAK,                       & ! in
                          pft_mask(:),                          & ! in
                          is_woody(:),                          & ! in
                          no_children,                          & ! in
                          fract_st(:,:),                        & ! in
                          cover_frac_pot_st(:,:),               & ! in
                          veg_frac_correction_st(:,:),          & ! in
                          q_rel_air(:),                         & ! in
                          wind_10m(:),                          & ! in
                          !tile%mem(CARBON_)%p ,                 & ! in R: only to find FUNCTION Average_over_children
                          max_wind_10m_act(:),                  & ! inout
                          prev_day_max_wind_10m(:),             & ! inout
                          max_wind_10m(:),                      & ! inout
                          q_rel_air_climbuf(:),                 & ! inout
                          q_rel_air_climbuf_yDay(:),            & ! out
                          damaged_frac_st(:,:)                  & ! inout
                          )


      ! Cumulate disturbed areas for diagnostics
!       disturbance%box_burned_frac_avg     (kidx0:kidx1,:) =                                                    &
!       disturbance%box_burned_frac_avg     (kidx0:kidx1,:) + disturbance%burned_frac       (kidx0:kidx1,:)      &
!                                                             * fract_st(:,:)      &
!                                                             * SPREAD(surface%veg_ratio_max(kidx0:kidx1),DIM=2, &
!                                                                      ncopies=no_children)                   &
!                                                             * 86400._wp
!       IF (.NOT. ASSOCIATED(disturbance%box_burned_frac_diag_avg,disturbance%box_burned_frac_avg))                &
!         disturbance%box_burned_frac_diag_avg(kidx0:kidx1,:) =                                                    &
!         disturbance%box_burned_frac_diag_avg(kidx0:kidx1,:) + disturbance%burned_frac_diag  (kidx0:kidx1,:)      &
!                                                               * fract_st(:,:)     &
!                                                               * SPREAD(surface%veg_ratio_max(kidx0:kidx1),DIM=2, &
!                                                                        ncopies=no_children)                   &
!                                                               * 86400._wp
!       disturbance%box_damaged_frac_avg    (kidx0:kidx1,:) =                                                    &
!       disturbance%box_damaged_frac_avg    (kidx0:kidx1,:) + disturbance%damaged_frac      (kidx0:kidx1,:)      &
!                                                             * fract_st(:,:)     &
!                                                             * SPREAD(surface%veg_ratio_max(kidx0:kidx1),DIM=2, &
!                                                                      ncopies=no_children)                   &
!                                                             * 86400._wp


      IF ( tile%Is_process_active(CARBON_) ) THEN

        ! Calculate carbon relocation connected with damages and fires
        ! R: hier sollen nur die Variablen übergeben werden, die von dieser Prozedur hier im Folgenden benötigt werden.
        !    Andere Variablen werden erst später bereit gestellt, denn in diesem Fall rufen wir ein weiteres
        !    Interface beim CARBON Prozess auf (SUBROUTINE relocate_disturbed_carbon).
        !    Erst dort sollen dann alle Variablen an die eigentlichen Prozess-Prozeduren übergeben werden.
        CALL relocate_disturbed_carbon(tile, options,                 & ! in, we need tile to find
                                       DIST_WINDBREAK,                & ! in
                                       no_children,                   & ! in
                                       damaged_frac_st(:,:),          & ! in, only for windbreak
                                       fract_st(:,:),                 & ! in
                                       veg_frac_correction_st(:,:),   & ! in
                                       Cpool_green_st(:,:),           & ! inout
                                       Cpool_woods_st(:,:),           & ! inout
                                       Cpool_reserve_st(:,:),         & ! inout
                                       carbon_2_GreenLitterPools_st(:,:), & ! inout
                                       carbon_2_WoodLitterPools_st(:,:),  & ! inout
                                       YCpool_acid_ag1_st(:,:),       & ! inout
                                       YCpool_water_ag1_st(:,:),      & ! inout
                                       YCpool_ethanol_ag1_st(:,:),    & ! inout
                                       YCpool_nonsoluble_ag1_st(:,:), & ! inout
                                       YCpool_acid_ag2_st(:,:),       & ! inout
                                       YCpool_water_ag2_st(:,:),      & ! inout
                                       YCpool_ethanol_ag2_st(:,:),    & ! inout
                                       YCpool_nonsoluble_ag2_st(:,:), & ! inout
                                       !
                                       YCpool_acid_bg1_st,       & ! inout
                                       YCpool_water_bg1_st,      & ! inout
                                       YCpool_ethanol_bg1_st,    & ! inout
                                       YCpool_nonsoluble_bg1_st, & ! inout
                                       YCpool_acid_bg2_st,       & ! inout
                                       YCpool_water_bg2_st,      & ! inout
                                       YCpool_ethanol_bg2_st,    & ! inout
                                       YCpool_nonsoluble_bg2_st, & ! inout
                                       !
                                       YCpool_humus_1_st(:,:),   & ! inout
                                       YCpool_humus_2_st(:,:),   & ! inout
                                       carbon_2_atmos_st_ta      & ! out, is used only with fire
                                       )

        CALL relocate_disturbed_carbon(tile, options,                 & ! in
                                       DIST_FIRE,                     & ! in
                                       no_children,                   & ! in
                                       burned_frac_st(:,:),           & ! in, only for fire
                                       fract_st(:,:),                 & ! in
                                       veg_frac_correction_st(:,:),   & ! in
                                       Cpool_green_st(:,:),           & ! inout
                                       Cpool_woods_st(:,:),           & ! inout
                                       Cpool_reserve_st(:,:),         & ! inout
                                       carbon_2_GreenLitterPools_st(:,:), & ! inout
                                       carbon_2_WoodLitterPools_st(:,:),  & ! inout
                                       YCpool_acid_ag1_st(:,:),       & ! inout
                                       YCpool_water_ag1_st(:,:),      & ! inout
                                       YCpool_ethanol_ag1_st(:,:),    & ! inout
                                       YCpool_nonsoluble_ag1_st(:,:), & ! inout
                                       YCpool_acid_ag2_st(:,:),       & ! inout
                                       YCpool_water_ag2_st(:,:),      & ! inout
                                       YCpool_ethanol_ag2_st(:,:),    & ! inout
                                       YCpool_nonsoluble_ag2_st(:,:), & ! inout
                                         !
                                       YCpool_acid_bg1_st,       & ! inout
                                       YCpool_water_bg1_st,      & ! inout
                                       YCpool_ethanol_bg1_st,    & ! inout
                                       YCpool_nonsoluble_bg1_st, & ! inout
                                       YCpool_acid_bg2_st,       & ! inout
                                       YCpool_water_bg2_st,      & ! inout
                                       YCpool_ethanol_bg2_st,    & ! inout
                                       YCpool_nonsoluble_bg2_st, & ! inout
                                       !
                                       YCpool_humus_1_st(:,:),    & ! inout
                                       YCpool_humus_2_st(:,:),    & ! inout
                                       carbon_2_atmos_st_ta(:,:)  & ! out, is used only with fire
                                       )

      ENDIF ! is_process_active(CARBON_)


        ! R: This corresponds to JS3, but I did not understand it...
        !    disturbance%disturb_2_atmos corresponds to dsl4jsb_Get_var2D_onChunk(CARBON_,  carbon_2_atmos )
        !    However, I guess it is better not to take this from JS3.
        !IF ( tile%Is_process_active(CARBON_) ) THEN
        !  ! Calculate CO2 and other fluxes to the atmosphere (updated once a day)
        !  CO2_emission(:) = disturbance%disturb_2_atmos(kidx0:kidx1)   * molarMassCO2_kg/86400._wp
        !ENDIF

      ! R: Bei ähnlichen Prozessen, bei denen auf Parent Dinge berechnet werden, die dann auf die Childes über ein festes Schema (nicht 
      !    individuell parametrisiert) zurück gegeben werden, sollte sowas langfristig eigentlich von einer INHERIT Funktion übernommen 
      !    werden, die der Aggregate Prozedur entspricht.
      !    Das Problem wäre, daß das Message System das Aggregate alternativ zum Integrate ausführt und nicht die gesamte HSM nach dem
      !    Integrage für das Agregate durchläuft. Daher kann man ein INHERT für einen Task nicht genau so hinzufügen wie ein AGGREGATE!
      !    Langfristig wäre es aber strukturell vollständiger die HSM jeweils für INTEGRATE, AGGREGATE und -je nach Task evtl. -
      !    INHERIT  zu durchlaufen. Dazu müsste man in FUNCTION Process_task_on_tile unter case(INTEGRATE) nur nach dem letzten Sibling
      !    nicht auf das Parent springen und "AGGREGATE" senden (und entsprechend dort nicht auf INTEGRATE zurück schalten).
      !    Stattdessen müsste man in SUBROUTINE Run_tasks die HSM nochmal für AGGREGATE und INHERIT starten!
      !    Wie auch immer: bei Berechnungen, die die Childs parametrisieren sollte man es wie im Folgenden machen, da sich eine Funktion 
      !    wohl nicht lohnt:
      !
      ! Now we put our results calculated on the current "parent" tile for its corresponding childs 
      ! * into the "parent" variable
      ! * down to the corresponding child tiles

      burned_frac = SUM( burned_frac_st(:,:) * fract_st(:,:) , DIM=2) ! Note, this works because the fract_st sum over DIM=2 
      SELECT TYPE (parentmem=>tile%mem(DISTURB_)%p)
        CLASS IS (t_disturb_memory)! is always = 1 !
          CALL parentmem%HandDown_var_real3d(ics, ice, iblk, no_children, parentmem%burned_frac,   &
                                            burned_frac_st(:,:)   ) ! inout
          damaged_frac = SUM( damaged_frac_st(:,:) * fract_st(:,:) , DIM=2)
          CALL parentmem%HandDown_var_real3d(ics, ice, iblk, no_children, parentmem%damaged_frac,  &
                                            damaged_frac_st(:,:)  ) ! inout
      END SELECT

      IF ( tile%Is_process_active(CARBON_) ) THEN
        SELECT TYPE (parentmem=>tile%mem(CARBON_)%p)
          CLASS IS (t_carbon_memory)
            Cpool_green    = SUM( Cpool_green_st(:,:) * fract_st(:,:) , DIM=2)
            Cpool_green_ta = Cpool_green              * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%Cpool_green,   &
                                            Cpool_green_st(:,:)       )

            Cpool_woods    = SUM( Cpool_woods_st(:,:) * fract_st(:,:) , DIM=2)
            Cpool_woods_ta = Cpool_woods              * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%Cpool_woods,   &
                                                Cpool_woods_st(:,:)      )

            Cpool_reserve    = SUM( Cpool_reserve_st(:,:) * fract_st(:,:) , DIM=2)
            Cpool_reserve_ta = Cpool_reserve              * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%Cpool_reserve, &
                                                Cpool_reserve_st(:,:)    )

            carbon_2_GreenLitterPools = SUM( carbon_2_GreenLitterPools_st(:,:) * fract_st(:,:) , DIM=2)
            carbon_2_GreenLitterPools_ta = carbon_2_GreenLitterPools           * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%carbon_2_GreenLitterPools, &
                                                carbon_2_GreenLitterPools_st(:,:)    )

            carbon_2_WoodLitterPools = SUM( carbon_2_WoodLitterPools_st(:,:) * fract_st(:,:) , DIM=2)
            carbon_2_WoodLitterPools_ta = carbon_2_WoodLitterPools           * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%carbon_2_WoodLitterPools, &
                                                carbon_2_WoodLitterPools_st(:,:)    )

            YCpool_acid_ag1 = SUM( YCpool_acid_ag1_st(:,:) * fract_st(:,:) , DIM=2)
            YCpool_acid_ag1_ta = YCpool_acid_ag1           * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%YCpool_acid_ag1,      &
                                                YCpool_acid_ag1_st(:,:)       )

            YCpool_water_ag1 = SUM( YCpool_water_ag1_st(:,:) * fract_st(:,:) , DIM=2)
            YCpool_water_ag1_ta = YCpool_water_ag1           * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%YCpool_water_ag1,     &
                                                YCpool_water_ag1_st(:,:)      )

            YCpool_ethanol_ag1 = SUM( YCpool_ethanol_ag1_st(:,:) * fract_st(:,:) , DIM=2)
            YCpool_ethanol_ag1_ta = YCpool_ethanol_ag1           * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%YCpool_ethanol_ag1,   &
                                                YCpool_ethanol_ag1_st(:,:)    )

            YCpool_nonsoluble_ag1 = SUM( YCpool_nonsoluble_ag1_st(:,:) * fract_st(:,:) , DIM=2)
            YCpool_nonsoluble_ag1_ta = YCpool_nonsoluble_ag1           * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%YCpool_nonsoluble_ag1, &
                                                YCpool_nonsoluble_ag1_st(:,:) )

            YCpool_acid_ag2 = SUM( YCpool_acid_ag2_st(:,:) * fract_st(:,:) , DIM=2)
            YCpool_acid_ag2_ta = YCpool_acid_ag2           * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%YCpool_acid_ag2,       &
                                                YCpool_acid_ag2_st(:,:)       )

            YCpool_water_ag2 = SUM( YCpool_water_ag2_st(:,:) * fract_st(:,:) , DIM=2)
            YCpool_water_ag2_ta = YCpool_water_ag2           * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%YCpool_water_ag2,      &
                                                YCpool_water_ag2_st(:,:)      )

            YCpool_ethanol_ag2 = SUM( YCpool_ethanol_ag2_st(:,:) * fract_st(:,:) , DIM=2)
            YCpool_ethanol_ag2_ta = YCpool_ethanol_ag2           * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%YCpool_ethanol_ag2,     &
                                                YCpool_ethanol_ag2_st(:,:)    )

            YCpool_nonsoluble_ag2 = SUM( YCpool_nonsoluble_ag2_st(:,:) * fract_st(:,:) , DIM=2)
            YCpool_nonsoluble_ag2_ta = YCpool_nonsoluble_ag2           * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%YCpool_nonsoluble_ag2, &
                                                YCpool_nonsoluble_ag2_st(:,:) )

            YCpool_humus_1 = SUM( YCpool_humus_1_st(:,:) * fract_st(:,:) , DIM=2)
            YCpool_humus_1_ta = YCpool_humus_1           * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%YCpool_humus_1,        &
                                                YCpool_humus_1_st(:,:) )

            YCpool_humus_2 = SUM( YCpool_humus_2_st(:,:) * fract_st(:,:) , DIM=2)
            YCpool_humus_2_ta = YCpool_humus_2           * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%YCpool_humus_2,        &
                                                YCpool_humus_2_st(:,:) )

            carbon_2_atmos = SUM( carbon_2_atmos_st_ta(:,:) * fract_st(:,:) , DIM=2)
            carbon_2_atmos_ta = carbon_2_atmos           * veg_frac_correction
            CALL parentmem%HandDown_var_real3d( ics, ice, iblk,  no_children,   parentmem%carbon_2_atmos,        &
                                                carbon_2_atmos_st_ta(:,:)    )
        END SELECT
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

    dsl4jsb_Def_memory(DISTURB_)
    dsl4jsb_Def_memory(CARBON_)
    !dsl4jsb_Def_config(DISTURB_)

    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract

    TYPE(t_jsb_model), POINTER                :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_natural_disturbance'

    INTEGER  :: iblk, ics, ice, nc

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(DISTURB_)
    dsl4jsb_Get_memory(CARBON_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(DISTURB_, burned_frac,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(DISTURB_, damaged_frac,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(DISTURB_, fuel,                   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(DISTURB_, q_rel_air_climbuf,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(DISTURB_, q_rel_air_climbuf_yDay, weighted_by_fract)

    ! R: I do not think these variables should be aggregated:
    !dsl4jsb_Aggregate_onChunk(DISTURB_, max_wind_10m_act,       weighted_by_fract)
    !dsl4jsb_Aggregate_onChunk(DISTURB_, prev_day_max_wind_10m,  weighted_by_fract)
    !dsl4jsb_Aggregate_onChunk(DISTURB_, max_wind_10m,           weighted_by_fract)
    !dsl4jsb_Aggregate_onChunk(DISTURB_, q_rel_air_climbuf_yDay, weighted_by_fract)

    ! R: Die Effekte von Disturbance auf Carbon-Variablen wurden durch 
    !    CALL relocate_disturbed_carbon ausgerechnet. Entsprechend wurden die C-Variablen
    !    für das veg und dessen child tiles auch geändert. Nun müssen diese Änderungen aber
    !    auch nach oben akkumuliert werden:
    IF ( tile%Is_process_active(CARBON_) ) THEN
       dsl4jsb_Aggregate_onChunk(CARBON_,  Cpool_green,               weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  Cpool_woods,               weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  Cpool_reserve,             weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  carbon_2_GreenLitterPools, weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  carbon_2_WoodLitterPools,  weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_acid_ag1,           weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_water_ag1,          weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_ethanol_ag1,        weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_nonsoluble_ag1,     weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_acid_ag2,           weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_water_ag2,          weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_ethanol_ag2,        weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_nonsoluble_ag2,     weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_humus_1,            weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_humus_2,            weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  carbon_2_atmos,            weighted_by_fract)

       dsl4jsb_Aggregate_onChunk(CARBON_,  Cpool_green_ta,               weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  Cpool_woods_ta,               weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  Cpool_reserve_ta,             weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  carbon_2_GreenLitterPools_ta, weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  carbon_2_WoodLitterPools_ta,  weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_acid_ag1_ta,           weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_water_ag1_ta,          weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_ethanol_ag1_ta,        weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_nonsoluble_ag1_ta,     weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_acid_ag2_ta,           weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_water_ag2_ta,          weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_ethanol_ag2_ta,        weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_nonsoluble_ag2_ta,     weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_humus_1_ta,            weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  YCpool_humus_2_ta,            weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  carbon_2_atmos_ta,            weighted_by_fract)
    ENDIF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_natural_disturbances

#endif
END MODULE mo_disturb_interface
