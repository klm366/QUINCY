!> Summary: Contains the interfaces to the HD model
!>
!> @author
!>  Veronika Gayler, MPI-M Hamburg
!>
!> @par Revision History
!> First version by Veronika Gayler (2013-04-20)
!>
!> @par Copyright
!> This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!> Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!>

MODULE mo_hd_interface
!#ifndef __NO_JSBACH__
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: finish, message, message_text
  USE mo_sync,             ONLY: SYNC_C, sync_patch_array, global_sum_array
  USE mo_sync,             ONLY: global_max

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_time,           ONLY: is_time_experiment_start, is_newday
  USE mo_jsb_math_constants, ONLY: deg2rad, pi
  USE mo_jsb_parallel,       ONLY: Get_omp_thread

  USE mo_model_domain,       ONLY: p_patch

  dsl4jsb_Use_processes HD_, HYDRO_
  dsl4jsb_Use_config(HD_)
  dsl4jsb_Use_memory(HD_)
  dsl4jsb_Use_memory(HYDRO_)

  USE mo_hd_reservoir_cascade, ONLY: reservoir_cascade, evapo_wetskin
  USE mo_jsb_physical_constants, ONLY: rhoh2o

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_hd, hd_lateral_flow, lake_lateral_flow, hd_check_water_budget

  PUBLIC :: Register_hd_tasks

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_hd_local_outflow
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_hd
    PROCEDURE, NOPASS :: Aggregate => aggregate_hd
  END TYPE tsk_hd_local_outflow

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_dynlake
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_dynlake
    PROCEDURE, NOPASS :: Aggregate => aggregate_dynlake
  END TYPE tsk_dynlake

  INTERFACE tsk_hd_local_outflow
    PROCEDURE Create_task_hd_local_outflow
  END INTERFACE tsk_hd_local_outflow

  INTERFACE tsk_dynlake
    PROCEDURE Create_task_dynlake
  END INTERFACE tsk_dynlake

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_interface'

CONTAINS

  SUBROUTINE Register_hd_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_hd_local_outflow(model_id))
    CALL this%Register_task(tsk_dynlake(model_id))

  END SUBROUTINE Register_hd_tasks

  FUNCTION Create_task_hd_local_outflow(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_hd_local_outflow::return_ptr)
    CALL return_ptr%Construct(name='hd_local_outflow', process_id=HD_, owner_model_id=model_id)

  END FUNCTION Create_task_hd_local_outflow


  FUNCTION Create_task_dynlake(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_dynlake::return_ptr)
    CALL return_ptr%Construct(name='dynlake', process_id=HD_,owner_model_id=model_id)

  END FUNCTION Create_task_dynlake

  !------------------------------------------------------------------------------------------------
  SUBROUTINE update_hd(tile, options)
  !------------------------------------------------------------------------------------------------
  ! the routine is called within the block-loop from jsbach_interface. There is no information on
  ! neighbouring grid cells (on other patches) availabe.
  ! We calculate here the local outflow of each grid cell.
  !------------------------------------------------------------------------------------------------

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_tile_abstract), POINTER       :: box_tile, lake_tile, land_tile, veg_tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: day_1y

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HD_)
    dsl4jsb_Def_memory(HYDRO_)

    TYPE(t_hydro_memory), POINTER :: lake_hydro__mem, land_hydro__mem

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: runoff
    dsl4jsb_Real2D_onChunk :: drainage
    dsl4jsb_Real2D_onChunk :: discharge
    dsl4jsb_Real2D_onChunk :: w_river
    dsl4jsb_Real2D_onChunk :: w_overl
    dsl4jsb_Real2D_onChunk :: evapo_river
    dsl4jsb_Real2D_onChunk :: evapo_overl
    dsl4jsb_Real2D_onChunk :: oro_stddev
    dsl4jsb_Real2D_onChunk :: outflow_outlet
    dsl4jsb_Real2D_onChunk :: outflow_runoff
    dsl4jsb_Real2D_onChunk :: outflow_drainage
    dsl4jsb_Real2D_onChunk :: outflow_rivers
    dsl4jsb_Real2D_onChunk :: ret_overlflow
    dsl4jsb_Real2D_onChunk :: nres_overlflow
    dsl4jsb_Real2D_onChunk :: ret_baseflow
    dsl4jsb_Real2D_onChunk :: nres_baseflow
    dsl4jsb_Real2D_onChunk :: ret_riverflow
    dsl4jsb_Real2D_onChunk :: nres_riverflow
    dsl4jsb_Real2D_onChunk :: basin_mask
    dsl4jsb_Real2D_onChunk :: runoff_lake
    dsl4jsb_Real2D_onChunk :: w_river_res
    dsl4jsb_Real2D_onChunk :: w_overl_res
    dsl4jsb_Real2D_onChunk :: w_river_res_old
    dsl4jsb_Real2D_onChunk :: w_overl_res_old
    dsl4jsb_Real2D_onChunk :: fract_river
    dsl4jsb_Real2D_onChunk :: fract_overl
    dsl4jsb_Real2D_onChunk :: volume_river_seas_mean
    dsl4jsb_Real2D_onChunk :: volume_overl_seas_mean
    dsl4jsb_Real3D_onChunk :: overlflow_res
    dsl4jsb_Real3D_onChunk :: baseflow_res
    dsl4jsb_Real3D_onChunk :: riverflow_res
    dsl4jsb_Real3D_onChunk :: volume_river_seas
    dsl4jsb_Real3D_onChunk :: volume_overl_seas
    dsl4jsb_Real2D_onChunk :: local_fluxes
    dsl4jsb_Real2D_onChunk :: local_budget
    dsl4jsb_Real2D_onChunk :: local_wbal_error

    REAL(wp), POINTER :: area(:)
    REAL(wp), POINTER :: lake_runoff(:)
    REAL(wp), POINTER :: land_runoff(:)
    REAL(wp), POINTER :: lake_fract(:)
    REAL(wp), POINTER :: land_fract(:)
    REAL(wp), POINTER :: veg_fract(:)

    REAL(wp), DIMENSION(options%nc) :: &
      & tile_fract, &
      & inflow,     &
      & local_budget_old

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_hd'

    CHARACTER(len=32) :: routing_scheme

    INTEGER  :: iblk, ics, ice, nc, ndays
    REAL(wp) :: steplen

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen

    IF (ASSOCIATED(tile%parent)) &
      & CALL finish(TRIM(routine), 'HD model works on root tile only!')

    model => Get_model(tile%owner_model_id)

    IF (.NOT. tile%Is_process_calculated(HD_)) RETURN

    dsl4jsb_Get_config(HD_)

    IF (.NOT. model%Is_process_enabled(HD_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    grid => Get_grid(model%grid_id)

    ! Configuration parameters from the namelist
    routing_scheme    = dsl4jsb_Config(HD_)%routing_scheme

    CALL model%Get_top_tile(box_tile)
    lake_tile => box_tile%first_child_tile
    land_tile => lake_tile%next_sibling_tile
    veg_tile => land_tile%first_child_tile%next_sibling_tile

    day_1y => Get_vgrid('days_in_year')
    ndays = day_1y%n_levels    

    ! Get reference to variables in memory
    !
    dsl4jsb_Get_memory(HD_)
    !dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    SELECT TYPE (mm => lake_tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); lake_hydro__mem => mm; END SELECT
    SELECT TYPE (mm => land_tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); land_hydro__mem => mm; END SELECT

    lake_runoff => lake_hydro__mem%runoff%ptr(ics:ice,iblk)
    land_runoff => land_hydro__mem%runoff%ptr(ics:ice,iblk)
    lake_fract => lake_tile%fract(ics:ice,iblk)
    land_fract => land_tile%fract(ics:ice,iblk)
    veg_fract => veg_tile%fract(ics:ice,iblk)

    area => grid %area (ics:ice, iblk)

    CALL tile%Get_fraction(ics, ice, iblk, fract=tile_fract(:))

    dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff)      ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, drainage)    ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, discharge)   ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_river) ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_overl) ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, w_river)     ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, w_overl)     ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, oro_stddev)  ! in

    dsl4jsb_Get_var2D_onChunk(HD_, outflow_outlet)    ! in
    dsl4jsb_Get_var2D_onChunk(HD_, ret_overlflow)    ! in
    dsl4jsb_Get_var2D_onChunk(HD_, nres_overlflow)   ! in
    dsl4jsb_Get_var2D_onChunk(HD_, ret_baseflow)     ! in
    dsl4jsb_Get_var2D_onChunk(HD_, nres_baseflow)    ! in
    dsl4jsb_Get_var2D_onChunk(HD_, ret_riverflow)    ! in
    dsl4jsb_Get_var2D_onChunk(HD_, nres_riverflow)   ! in
    dsl4jsb_Get_var2D_onChunk(HD_, basin_mask)       ! in
    dsl4jsb_Get_var2D_onChunk(HD_, runoff_lake)      ! in
    dsl4jsb_Get_var2D_onChunk(HD_, volume_river_seas_mean)! inout
    dsl4jsb_Get_var2D_onChunk(HD_, volume_overl_seas_mean)! inout
    dsl4jsb_Get_var3D_onChunk(HD_, overlflow_res)    ! inout
    dsl4jsb_Get_var3D_onChunk(HD_, baseflow_res)     ! inout
    dsl4jsb_Get_var3D_onChunk(HD_, riverflow_res)    ! inout
    dsl4jsb_Get_var3D_onChunk(HD_, volume_river_seas)! inout
    dsl4jsb_Get_var3D_onChunk(HD_, volume_overl_seas)! inout
    dsl4jsb_Get_var2D_onChunk(HD_, outflow_runoff)   ! out
    dsl4jsb_Get_var2D_onChunk(HD_, outflow_drainage) ! out
    dsl4jsb_Get_var2D_onChunk(HD_, outflow_rivers)   ! out
    dsl4jsb_Get_var2D_onChunk(HD_, w_river_res)      ! out
    dsl4jsb_Get_var2D_onChunk(HD_, w_overl_res)      ! out
    dsl4jsb_Get_var2D_onChunk(HD_, w_river_res_old)  ! out
    dsl4jsb_Get_var2D_onChunk(HD_, w_overl_res_old)  ! out
    dsl4jsb_Get_var2D_onChunk(HD_, fract_river)      ! out
    dsl4jsb_Get_var2D_onChunk(HD_, fract_overl)      ! out
    dsl4jsb_Get_var2D_onChunk(HD_, local_fluxes)     ! out
    dsl4jsb_Get_var2D_onChunk(HD_, local_budget)     ! out
    dsl4jsb_Get_var2D_onChunk(HD_, local_wbal_error) ! out

    SELECT CASE (TRIM(routing_scheme))
    CASE ('full')

      ! local reservoir content for water balance check
      local_budget_old(:) = local_budget(:)

      IF (model%config%use_wetskins) THEN
        ! reservoir wetskin evaporation [m3]
        CALL evapo_wetskin (w_river(:)*area, veg_fract(:), nc, riverflow_res(:,:))
        CALL evapo_wetskin (w_overl(:)*area, veg_fract(:), nc, overlflow_res(:,:))
      ENDIF


      ! unit conversions: [kg m-2 s-1] -> [m3/s]
      ! use_dynlake = .FALSE.:  basin_mask is set to zero
      WHERE (basin_mask(:) > 0._wp)
        runoff_lake(:) = lake_runoff(:) * lake_fract(:) * area / rhoh2o
        inflow(:)      = land_runoff(:) * land_fract(:) * area / rhoh2o
      ELSEWHERE
        inflow(:)      = runoff(:) * tile_fract(:) * area / rhoh2o
      END WHERE
      inflow(:)  = inflow(:) + outflow_outlet(:)
      CALL reservoir_cascade (inflow(:), nc, steplen, ret_overlflow(:), nres_overlflow(:), overlflow_res(:,:), outflow_runoff(:))

      inflow(:) = drainage(:) * tile_fract(:) * area / rhoh2o
      CALL reservoir_cascade (inflow(:), nc, steplen, ret_baseflow(:),  nres_baseflow(:),  baseflow_res(:,:),  outflow_drainage(:))

      inflow(:) = discharge(:)
      CALL reservoir_cascade (inflow(:), nc, steplen, ret_riverflow(:), nres_riverflow(:), riverflow_res(:,:), outflow_rivers(:))

      ! river reservoir water content: [m3] 
      WHERE (veg_fract(:) > 0._wp)
         w_river_res(:) = SUM(riverflow_res(:,:), DIM=2, MASK=riverflow_res(:,:)>0._wp)/area
         w_overl_res(:) = SUM(overlflow_res(:,:), DIM=2, MASK=overlflow_res(:,:)>0._wp)/area
      ELSE WHERE
         w_river_res(:) = 0._wp
         w_overl_res(:) = 0._wp
      END WHERE

      IF (is_newday(options%current_datetime, options%dtime)) THEN
        volume_river_seas(:,1:ndays-1) = volume_river_seas(:,2:ndays)
        volume_overl_seas(:,1:ndays-1) = volume_overl_seas(:,2:ndays)
        volume_river_seas(:,ndays)     = 0._wp
        volume_overl_seas(:,ndays)     = 0._wp
      END IF

      volume_river_seas(:,ndays) = volume_river_seas(:,ndays) + (w_river_res(:) - w_river_res_old(:)) 
      volume_overl_seas(:,ndays) = volume_overl_seas(:,ndays) + (w_overl_res(:) - w_overl_res_old(:)) 


      WHERE(w_river_res(:)>EPSILON(1._wp))
        volume_river_seas_mean(:) = MAX(0._wp, SUM(volume_river_seas(:,:),DIM=2))
      ELSE WHERE
        volume_river_seas_mean(:) = 0._wp
      END WHERE
      WHERE(w_overl_res(:)>EPSILON(1._wp))
        volume_overl_seas_mean(:) = MAX(0._wp, SUM(volume_overl_seas(:,:),DIM=2))
      ELSE WHERE
        volume_overl_seas_mean(:) = 0._wp
      END WHERE

      IF (model%config%use_wetskins) THEN
        ! update reservoir wetskin fraction
        CALL update_wetskin (fract_river(:),volume_river_seas_mean(:),oro_stddev(:),veg_fract(:),area,nc)
        CALL update_wetskin (fract_overl(:),volume_overl_seas_mean(:),oro_stddev(:),veg_fract(:),area,nc)
      ENDIF

      w_river_res_old(:) = w_river_res(:)
      w_overl_res_old(:) = w_overl_res(:)


      ! 1. Check water balance locally

      ! Reservoir content  [m3]
      local_budget(:) =  SUM(overlflow_res(:,:), DIM=2) & !< water in overland flow reservoirs
        &              + SUM(baseflow_res (:,:), DIM=2) & !< water in baseflow reservoir
        &              + SUM(riverflow_res(:,:), DIM=2)   !< water in riverflow reservoir

      IF (.NOT. is_time_experiment_start(options%current_datetime)) THEN
        ! Water fluxes [m3/s]
        local_fluxes(:) =  runoff(:)   * tile_fract(:) * area(:) / rhoh2o         &
          &              + drainage(:) * tile_fract(:) * area(:) / rhoh2o         &
          &              + discharge(:)                                           &
          &              + evapo_river(:) * area(:) / steplen                     & !< river evaporation [m] -> [m3/s] 
          &              + evapo_overl(:) * area(:) / steplen                     & !< river evaporation [m] -> [m3/s] 
          &              - lake_runoff(:) * lake_fract(:) * area / rhoh2o         &
          &              - outflow_runoff(:) - outflow_drainage(:) - outflow_rivers(:)

        ! Local water balance error [m3]
        local_wbal_error(:) = local_budget_old(:) + local_fluxes(:) * steplen - local_budget(:)
      END IF

    CASE ('weighted_to_coast')
      !
      ! dummy version of HD model
      !
      ! unit conversion: [kg m-2 s-1] -> [m3/s]
      discharge(:) = tile_fract(:) * (runoff(:) + drainage(:)) * area(:) / rhoh2o

    CASE ('zero')
      !
      ! even more dummy: set discharge to zero
      !
      discharge(:) = 0._wp

    CASE DEFAULT
      WRITE (message_text,*)  'no valid parameter for routing_scheme: ', TRIM(routing_scheme)
      CALL finish(TRIM(routine), message_text)

    END SELECT

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_hd

  SUBROUTINE aggregate_hd(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    IF (options%nc > 0) CONTINUE ! only to avoid compiler warnings
    IF (tile%level > 0) CONTINUE ! only to avoid compiler warnings

    ! Nothing to do because the HD model only runs on the top tile

  END SUBROUTINE aggregate_hd


  SUBROUTINE update_wetskin(fract,height,oro_stddev,fract_veg,cell_area,nc)

    ! Inputdata
    REAL(wp), INTENT(in) :: height(:)               !< Average flooding volume        [m]
    REAL(wp), INTENT(in) :: oro_stddev(:)           !< Orographic standard deviation  [m]
    REAL(wp), INTENT(in) :: fract_veg(:)            !< Vegetation tile fraction       [/]
    REAL(wp), INTENT(in) :: cell_area(:)            !< Grid cell area                 [m2]
    INTEGER, Intent(in)  :: nc                      !< Number of grid cells           [/]

    ! Result
    REAL(wp), INTENT(inout) :: fract(:)             !< New flooding fraction          [/]

    ! Local variables
    REAL(wp)             :: volume_max(nc)          !< Orographic scaling parameter   [m3]
    REAL(wp)             :: volume(nc)              !< Orographic scaling parameter   [m3]

    ! Local parameters
    REAL(wp), PARAMETER  :: eps = EPSILON(1.0_wp)   !< Really small value
    REAL(wp), PARAMETER  :: shape_par = 2._wp       !< shape parameter (1: cone, 2: parabolic)

    ! Orographic scaling parameter
    volume_max  = MAX(1._wp, 2._wp * oro_stddev * cell_area) / (1._wp + 2._wp / shape_par)
    volume      = height * cell_area

    WHERE (fract_veg < 0._wp)                      ! valid grid cell
      fract = 0._wp
    ELSE WHERE

      WHERE (volume < eps)  
        ! 1. Wetland dessicated
        fract = 0.0_wp

      ELSE WHERE

        fract = MIN(fract_veg, MAX(0.0_wp, (volume / volume_max) ** (2._wp / (2._wp + shape_par))))

      END WHERE

    END WHERE

  END SUBROUTINE update_wetskin


  SUBROUTINE update_dynlake(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    CLASS(t_jsb_tile_abstract), POINTER       :: box_tile, lake_tile, land_tile, fract_tile

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: dl_l

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HD_)

    ! Pointers to variables in memory

    dsl4jsb_Real2D_onChunk ::  &
      & height_oro_min,     &
      & basin_mask,         & 
      & height_lake,        &
      & area_lake,          &
      & fract_lake,          &
      & volume_lake,        &
      & height_dyn

    dsl4jsb_Real3D_onChunk ::  &
      & heights_lake,       &
      & areas_lake,         &
      & volumes_lake

    REAL(wp), POINTER :: fract_lake_old(:)
    REAL(wp), POINTER :: fract_land_old(:)
    REAL(wp), POINTER :: area(:)

    LOGICAL, ALLOCATABLE  ::  &
      & volume_mask_low(:,:), &        !> dummy field 
      & volume_mask_high(:,:),&        !> dummy field 
      & volume_mask(:)                 !> dummy field 

    REAL, ALLOCATABLE  ::  &
      & height_low(:),   &        !> dummy field 
      & height_high(:),  &        !> dummy field 
      & area_low(:),     &        !> dummy field 
      & area_high(:),    &        !> dummy field 
      & volume_low(:),   &        !> dummy field 
      & volume_high(:)               !> dummy field 

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_dynlake'

    CHARACTER(len=32) :: routing_scheme        ! river routing scheme

    INTEGER  :: iblk, ics, ice, nc, nlev
    INTEGER  :: no_omp_thread

    REAL(wp) :: steplen

    REAL(wp), DIMENSION(options%nc) :: fract_lake_changes

    no_omp_thread = Get_omp_thread()

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen

    IF (ASSOCIATED(tile%parent)) &
      & CALL finish(TRIM(routine), 'HD model works on root tile only!')

    model => Get_model(tile%owner_model_id)

    IF (.NOT. tile%Is_process_active(HD_)) RETURN

    dsl4jsb_Get_config(HD_)

    IF (.NOT. dsl4jsb_Config(HD_)%active) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    grid => Get_grid(model%grid_id)
    area => grid%area(ics:ice,iblk)

    dl_l => Get_vgrid('dynamic_lake_level')
    nlev    = dl_l%n_levels

    ! Configuration parameters from the namelist
    routing_scheme    = dsl4jsb_Config(HD_)%routing_scheme

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HD_)

    dsl4jsb_Get_var3D_onChunk(HD_, heights_lake)     ! IN
    dsl4jsb_Get_var3D_onChunk(HD_, areas_lake)       ! IN
    dsl4jsb_Get_var3D_onChunk(HD_, volumes_lake)     ! IN
    dsl4jsb_Get_var2D_onChunk(HD_, height_oro_min)   ! IN
    dsl4jsb_Get_var2D_onChunk(HD_, basin_mask)       ! IN
    dsl4jsb_Get_var2D_onChunk(HD_, height_lake)      ! OUT
    dsl4jsb_Get_var2D_onChunk(HD_, area_lake)        ! OUT
    dsl4jsb_Get_var2D_onChunk(HD_, fract_lake)        ! OUT
    dsl4jsb_Get_var2D_onChunk(HD_, volume_lake)      ! OUT
    dsl4jsb_Get_var2D_onChunk(HD_, height_dyn)       ! OUT

    CALL model%Get_top_tile(box_tile)
    lake_tile => box_tile%first_child_tile
    land_tile => lake_tile%next_sibling_tile

    fract_lake_old => lake_tile%fract(ics:ice,iblk)
    fract_land_old => land_tile%fract(ics:ice,iblk)

    ! Configuration parameters from the namelist
    routing_scheme = dsl4jsb_Config(HD_)%routing_scheme

    SELECT CASE (TRIM(routing_scheme))
    CASE ('full')

      ALLOCATE (volume_mask_low(nc,nlev))
      ALLOCATE (volume_mask_high(nc,nlev))
      ALLOCATE (volume_mask(nc))
      ALLOCATE (height_low(nc))
      ALLOCATE (height_high(nc))
      ALLOCATE (area_low(nc))
      ALLOCATE (area_high(nc))
      ALLOCATE (volume_low(nc))
      ALLOCATE (volume_high(nc))

      ! update grid cell lake height [m] and area [m3]
      volume_mask_low (:,:)=SPREAD(volume_lake(:), DIM=2, NCOPIES=100)>=volumes_lake(:,:)
      volume_mask_high(:,:)=SPREAD(volume_lake(:), DIM=2, NCOPIES=100)<=volumes_lake(:,:)
      volume_mask(:) = basin_mask(:) > 0._wp .AND. volume_lake(:)>=MINVAL(volumes_lake(:,:), DIM=2)

      WHERE ( basin_mask(:) > 0._wp  .AND. volume_lake(:)<MINVAL(volumes_lake(:,:), DIM=2) .AND. volume_lake(:)>0._wp )
         height_lake(:)= volume_lake(:) * heights_lake(:,1) / volumes_lake(:,1)
         fract_lake  (:)= volume_lake(:) * areas_lake  (:,1) / volumes_lake(:,1)
      ELSE WHERE (volume_mask(:) .AND. volume_lake(:)<MAXVAL(volumes_lake(:,:), DIM=2))
        height_low (:)=MAXVAL(heights_lake(:,:), DIM=2, MASK=volume_mask_low (:,:))
        height_high(:)=MINVAL(heights_lake(:,:), DIM=2, MASK=volume_mask_high(:,:))
        area_low   (:)=MAXVAL(areas_lake  (:,:), DIM=2, MASK=volume_mask_low (:,:))
        area_high  (:)=MINVAL(areas_lake  (:,:), DIM=2, MASK=volume_mask_high(:,:))
        volume_low (:)=MAXVAL(volumes_lake(:,:), DIM=2, MASK=volume_mask_low (:,:))
        volume_high(:)=MINVAL(volumes_lake(:,:), DIM=2, MASK=volume_mask_high(:,:))
        height_lake(:)= height_low(:) + (volume_lake(:)-volume_low(:)) * & 
                       (height_high(:)-height_low(:)) / &
                       (volume_high(:)-volume_low(:)) 
        fract_lake (:)= area_low  (:) + (volume_lake(:)-volume_low(:)) * & 
                       (area_high  (:)-area_low  (:)) / &
                       (volume_high(:)-volume_low(:))
      ELSE WHERE ( volume_mask(:) .AND. volume_lake(:)>=MAXVAL(volumes_lake(:,:), DIM=2) )
        height_lake(:) = MAXVAL(heights_lake(:,:), DIM=2) + & 
                         (volume_lake(:)-MAXVAL(volumes_lake(:,:), DIM=2)) / & 
                         MAXVAL(areas_lake(:,:), DIM=2)
        fract_lake(:)  = MAXVAL(areas_lake(:,:), DIM=2)
      ELSE WHERE
        height_lake(:) = 0._wp
        fract_lake(:)   = 0._wp
      END WHERE
      height_dyn(:) = height_oro_min(:) + height_lake(:)
      area_lake (:) = fract_lake(:) * area

      ! update lake fraction [%]
      WHERE (basin_mask(:) > 0._wp)
        ! lake and land fraction changes [%]
        fract_lake_changes(:) = MIN( MAX(fract_lake(:), 0._wp), 1._wp) - fract_lake_old(:)
      ELSE WHERE
        fract_lake_changes(:) = 0._wp
      END WHERE

      ! update fractions of lake tile, land tile and their children
      CALL model%Get_top_tile(fract_tile)
      IF (ASSOCIATED(fract_tile%parent)) CALL finish(routine, 'Should only run on box tile!')
      IF (model%config%use_sizelakes) THEN
      CALL do_stuff_on_tree(fract_tile, options, iblk, ics, ice, fract_lake_changes, fract_land_old(:))
      END IF

    CASE ('weighted_to_coast')
      WRITE (message_text,*)  'no valid parameter for dynlake: ', TRIM(routing_scheme)
      CALL finish(TRIM(routine), message_text)

    CASE ('zero')
      WRITE (message_text,*)  'no valid parameter for dynlake: ', TRIM(routing_scheme)
      CALL finish(TRIM(routine), message_text)

    CASE DEFAULT
      WRITE (message_text,*)  'no valid parameter for dynlake: ', TRIM(routing_scheme)
      CALL finish(TRIM(routine), message_text)

    END SELECT

  END SUBROUTINE update_dynlake

  SUBROUTINE aggregate_dynlake(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_dynlake' 
 
    CALL message(routine, 'aggregate dynlake task')

  END SUBROUTINE aggregate_dynlake


  RECURSIVE SUBROUTINE do_stuff_on_tree(cur, options, iblk, ics, ice, fract_lake_changes, fract_land) 

    CLASS(t_jsb_tile_abstract), POINTER :: cur      !< current tile
    CLASS(t_jsb_tile_abstract), POINTER :: child    !< child of cur

    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER, INTENT(in) :: iblk, ics, ice
    INTEGER :: i

    REAL(wp),INTENT(in) :: fract_lake_changes(:)    !< lake fraction changes
    REAL(wp),INTENT(in) :: fract_land(:)            !< lake fraction
    REAL(wp)            :: fract(options%nc)        !< fraction of current tile

    ! skipt box tile
    IF (TRIM(cur%name) /= 'box') THEN

      IF (iblk==1) CALL message('Dynlake at:', TRIM(cur%name))

      ! adjust fraction for lake, land and all children
      CALL cur%Get_fraction(ics, ice, iblk, fract=fract(:))
      IF(cur%is_lake) THEN
        fract(:) = fract(:) + fract_lake_changes(:)
      ELSE IF(cur%is_land) THEN
        fract(:) = fract(:) - fract_lake_changes(:)
      ELSE
        fract(:) = fract(:) - fract_lake_changes(:) * fract(:) / MAX(EPSILON(1._wp),fract_land(:))
      ENDIF
      CALL cur%Set_fraction(ics, ice, iblk, fract=fract(:))

    END IF

    IF (cur%Has_children()) THEN
      child => cur%first_child_tile
      DO i=1,cur%Get_no_of_children()
        CALL do_stuff_on_tree(child, options, iblk, ics, ice, fract_lake_changes, fract_land)
        child => child%next_sibling_tile
      END DO
    END IF

  END SUBROUTINE do_stuff_on_tree

  !------------------------------------------------------------------------------------------------
  SUBROUTINE hd_lateral_flow(tile, options)
  !------------------------------------------------------------------------------------------------
  ! The routine is called from jsbach_finish_timestep, after the loop over the nproma blocks.
  ! It is thus possible to synchronize the halos. And to use data from surounding grid cells.
  !------------------------------------------------------------------------------------------------

    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile
    TYPE(t_jsb_task_options),   INTENT(in) :: options

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER    :: grid

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HD_)
    dsl4jsb_Def_memory(HYDRO_)

    INTEGER                       :: jc, blk, n
    INTEGER                       :: jc_us, blk_us

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':hd_lateral_flow'

    ! Pointers to variables in memory

    dsl4jsb_Real2D_onDomain :: &
      & hd_mask,               &
      & basin_mask,            &
      & coast_ocean,           &
      & outflow_runoff,        &
      & outflow_drainage,      &
      & outflow_rivers,        &
      & outflow_resid,         &
      & outflow_count,         &
      & discharge,             &
      & discharge_ocean,       &
      & internal_drain,        &
      & discharge_lake,        &
      & riverflow_lake,        &
      & overlflow_lake,        &
      & baseflow_lake,         &
      & runoff_lake,           &
      & volume_lake,           &
      & water_budget_lake,     &
      & nsplit

    dsl4jsb_Real3D_onDomain :: &
      & nidx_upstream,         &
      & bidx_upstream

    REAL(wp) ::                 &
      & global_discharge,       &  ! global sum of discharge to the ocean
      & global_outflow,         &
      & global_internal_drain,  &
      & global_discharge_lake,      &
      & global_discharge_ocean, &
      & sum_weights,            &
      & steplen

    INTEGER  ::         &
      & nproma,         &
      & nblks,          &
      & blks, blke, jcs, jce

    REAL(wp), ALLOCATABLE :: cos_lat(:,:) ! Cosines of latitudes of cell centers

    LOGICAL, POINTER  ::  &
      & is_in_domain    (:,:) ! T: cell in domain (not halo)

    CHARACTER(len=32) :: &
      & routing_scheme        ! river routing scheme

    REAL(wp), ALLOCATABLE  :: &
      & in_domain       (:,:) ! 1: cell in domain, 0: halo cell

    REAL(wp) :: global_resid

    LOGICAL  :: debug_hd, use_bifurcated_rivers

    IF (options%nc > 0) CONTINUE ! only to avoid compiler warnings

    model => Get_model(tile%owner_model_id)

    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(HD_)

    IF (.NOT. model%Is_process_enabled(HD_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    IF (ASSOCIATED(tile%parent)) &
      & CALL finish(TRIM(routine), 'HD model works on root tile only!')

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HD_)
    !dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)

    nproma =  grid%nproma
    nblks  =  grid%nblks
    steplen = options%steplen

    dsl4jsb_Get_var2D_onDomain(HD_, hd_mask)          ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, basin_mask)       ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, outflow_runoff)   ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, outflow_drainage) ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, discharge_lake)   ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, riverflow_lake)   ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, overlflow_lake)   ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, baseflow_lake)    ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, runoff_lake)      ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, volume_lake)      ! INOUT
    dsl4jsb_Get_var2D_onDomain(HD_, water_budget_lake)! OUT
    dsl4jsb_Get_var2D_onDomain(HD_, outflow_rivers)   ! IN
    dsl4jsb_Get_var3D_onDomain(HD_, nidx_upstream)    ! IN
    dsl4jsb_Get_var3D_onDomain(HD_, bidx_upstream)    ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, nsplit)           ! IN


    dsl4jsb_Get_var2D_onDomain(HD_, coast_ocean)      ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, outflow_resid)    ! OUT
    ! dsl4jsb_Get_var2D_onDomain(HD_, outflow_count)    ! OUT

    dsl4jsb_Get_var2D_onDomain(HYDRO_, discharge)       ! OUT
    dsl4jsb_Get_var2D_onDomain(HYDRO_, discharge_ocean) ! OUT
    dsl4jsb_Get_var2D_onDomain(HYDRO_, internal_drain)  ! OUT

    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)

    ! Configuration parameters from the namelist
    routing_scheme        = dsl4jsb_Config(HD_)%routing_scheme
    debug_hd              = dsl4jsb_Config(HD_)%debug_hd
    use_bifurcated_rivers = dsl4jsb_Config(HD_)%use_bifurcated_rivers


    ! Domain Mask - to mask all halo cells for global sums (otherwise these cells are counted twice)
    ALLOCATE (in_domain(nproma,nblks))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSEWHERE
      in_domain = 0._wp
    END WHERE

    SELECT CASE (TRIM(routing_scheme))
    CASE ('full')

      ! Synchronize halos
      CALL sync_patch_array(sync_c, grid%patch, outflow_runoff(:,:))
      CALL sync_patch_array(sync_c, grid%patch, outflow_drainage(:,:))
      CALL sync_patch_array(sync_c, grid%patch, outflow_rivers(:,:))

      discharge(:,:) = 0._wp
      discharge_ocean(:,:) = 0._wp
      internal_drain(:,:) = 0._wp
      discharge_lake(:,:) = 0._wp
      riverflow_lake(:,:) = 0._wp

      ! Enable the next line for debugging ... but this only works on one (!) processor, i.e. no halos
      ! outflow_count(:,:) = nsplit(:,:)

      global_outflow  = global_sum_array((outflow_runoff(:,:)    & 
                                        + outflow_drainage(:,:)  & 
                                        + outflow_rivers(:,:)    & 
                                        + runoff_lake(:,:)) * in_domain(:,:))

      WHERE (hd_mask(:,:) == 2._wp .AND.  basin_mask(:,:) > 0._wp)
        ! lake water drained into internal lake
        riverflow_lake(:,:) = outflow_rivers(:,:)
        outflow_rivers(:,:) = 0._wp
      END WHERE

      IF (debug_hd) THEN
        ! Consistency check: There should not be outflow in pure ocean cells
        WHERE (hd_mask(:,:) == -1._wp)
          outflow_resid(:,:) = (outflow_runoff(:,:)+outflow_drainage(:,:)+outflow_rivers(:,:))*in_domain(:,:)
        END WHERE
        global_resid = global_sum_array(outflow_resid(:,:)*in_domain(:,:))
        IF (global_resid /= 0._wp) THEN
          WRITE (message_text,*) 'WARNING: Outflow to non-inflow ocean cells: ', global_resid, ' m3/s'
          CALL message (routine, message_text)
        END IF
      END IF

      blks = grid%get_blk_start()
      blke = grid%get_blk_end()
      DO blk = blks, blke
        jcs = grid%get_col_start(blk)
        jce = grid%get_col_end(blk)
        DO jc = jcs, jce
          ! inflow from grid cells upstream: outflow from runnoff, drainage and riverflow
          DO n = 1, hd__mem%nneigh
            ! local indices of grid cells upstream
            jc_us  =  nidx_upstream(jc,n,blk)
            blk_us =  bidx_upstream(jc,n,blk)
            IF (jc_us > 0 .and. blk_us > 0) THEN
            ! Workaround for NVIDIA compiler 21.3 internal compiler error
            ! Note that nsplit is set to 1._wp for use_bifurcated_rivers=false
#if defined(__NVCOMPILER)
              discharge(jc,blk) =  discharge(jc,blk)                      &
                                 + ( ( outflow_runoff(jc_us,blk_us)       &
                                     + outflow_drainage(jc_us,blk_us)     &
                                     + outflow_rivers(jc_us,blk_us)   )   &
                                   / nsplit(jc_us,blk_us)               )
#else
              IF (use_bifurcated_rivers) THEN
                discharge(jc,blk) =  discharge(jc,blk)                      &
                                   + ( ( outflow_runoff(jc_us,blk_us)       &
                                       + outflow_drainage(jc_us,blk_us)     &
                                       + outflow_rivers(jc_us,blk_us)   )   &
                                     / nsplit(jc_us,blk_us)               )
              ELSE
                discharge(jc,blk) =  discharge(jc,blk)               &
                                   + outflow_runoff(jc_us,blk_us)    &
                                   + outflow_drainage(jc_us,blk_us)  &
                                   + outflow_rivers(jc_us,blk_us)
              END IF
#endif
              ! Enable the next six lines for debugging ... but this only works on one (!) processor, i.e. no halos
              !outflow_count(jc_us,blk_us) = outflow_count(jc_us,blk_us) - 1._wp
              !IF (outflow_count(jc_us,blk_us) == 0) THEN
              !  outflow_runoff  (jc_us,blk_us) = 0._wp
              !  outflow_drainage(jc_us,blk_us) = 0._wp
              !  outflow_rivers  (jc_us,blk_us) = 0._wp
              !END IF
            END IF
          END DO
        END DO
      END DO

      ! Enable this block for debugging ... but this only works on one (!) processor, i.e. no halos
      !dbg>>
      ! print*, 'Outflow residual from lateral flow:'
      ! WHERE (hd_mask(:,:) == 2._wp)
      !   outflow_resid(:,:) = outflow_runoff(:,:) + outflow_drainage(:,:) + outflow_rivers(:,:)
      ! END WHERE
      ! print*, '  o  internal drainage cells    : ', global_sum_array(outflow_resid(:,:)*in_domain(:,:))
      ! outflow_resid(:,:) = 0._wp
      ! WHERE (hd_mask(:,:) == 1._wp)
      !   outflow_resid(:,:) = outflow_runoff(:,:) + outflow_drainage(:,:) + outflow_rivers(:,:)
      ! END WHERE
      ! print*, '  o  hd cells (mask=1)          : ', global_sum_array(outflow_resid(:,:)*in_domain(:,:))
      ! outflow_resid(:,:) = 0._wp
      ! WHERE (hd_mask(:,:) == 0._wp)
      !   outflow_resid(:,:) = outflow_runoff(:,:) + outflow_drainage(:,:) + outflow_rivers(:,:)
      ! END WHERE
      ! print*, '  o  ocean inflow cells (mask=0): ', global_sum_array(outflow_resid(:,:)*in_domain(:,:))
      ! outflow_resid(:,:) = 0._wp
      ! WHERE (hd_mask(:,:) == -1._wp)
      !   outflow_resid(:,:) = outflow_runoff(:,:) + outflow_drainage(:,:) + outflow_rivers(:,:)
      ! END WHERE
      ! print*, '  o  pure ocean cells (mask=-1) : ', global_sum_array(outflow_resid(:,:)*in_domain(:,:))
      ! WHERE (hd_mask(:,:) == 1._wp)
      !   outflow_resid(:,:) = outflow_runoff(:,:) + outflow_drainage(:,:) + outflow_rivers(:,:)
      ! END WHERE
      ! print*, '   => Problematic: mask = 1 : ', global_sum_array(outflow_resid(:,:)*in_domain(:,:))
      ! DO blk = blks, blke
      !   jcs = grid%get_col_start(blk)
      !   jce = grid%get_col_end(blk)
      !   DO jc = jcs, jce
      !     IF (outflow_resid(jc,blk) /= 0._wp) print*, 'AAA outflow_resid: ', jc,blk, &
      !       &      ' lon: ', grid%lon(jc,blk), 'lat: ', grid%lat(jc,blk), outflow_resid(jc,blk)
      !   END DO
      ! END DO
      !dbg<<

      WHERE (hd_mask(:,:) <= 0._wp)
        ! ocean inflow cells
        discharge_ocean(:,:) = outflow_runoff(:,:) + outflow_drainage(:,:) + discharge(:,:)
        outflow_runoff(:,:) = 0._wp
        outflow_drainage(:,:) = 0._wp
        discharge(:,:) = 0._wp
      ELSE WHERE (hd_mask(:,:) == 2._wp .AND. basin_mask(:,:) == 0._wp)
        ! internal drainage cells
        internal_drain(:,:) = outflow_runoff(:,:) + outflow_drainage(:,:) + discharge(:,:)
        outflow_runoff(:,:) = 0._wp
        outflow_drainage(:,:) = 0._wp
        discharge(:,:) = 0._wp
      ELSE WHERE (hd_mask(:,:) == 2._wp .AND.  basin_mask(:,:) > 0._wp)
        ! lake water drained into internal lake
        discharge(:,:) = discharge(:,:)
        baseflow_lake(:,:) = outflow_drainage(:,:)
        overlflow_lake(:,:) = outflow_runoff(:,:)
        outflow_runoff(:,:) = 0._wp
        outflow_drainage(:,:) = 0._wp
      END WHERE
      ! inflow at lake sink points
      discharge_lake(:,:) =  riverflow_lake(:,:) + overlflow_lake(:,:) + baseflow_lake(:,:) + runoff_lake(:,:)
      ! lake water storage [m3]
      water_budget_lake(:,:)  = riverflow_lake(:,:) + overlflow_lake(:,:) + baseflow_lake(:,:) + runoff_lake(:,:)
      volume_lake(:,:) = volume_lake(:,:) + water_budget_lake(:,:)

      ! global integrals
      global_discharge       = global_sum_array(discharge(:,:)*in_domain(:,:))
      global_internal_drain  = global_sum_array(internal_drain(:,:)*in_domain(:,:))
      global_discharge_lake  = global_sum_array(discharge_lake(:,:)*in_domain(:,:))
      global_discharge_ocean = global_sum_array(discharge_ocean(:,:)*in_domain(:,:))

      IF (debug_hd) THEN
        !
        ! 2. Check water conservation of lateral flow
        !
        WRITE (message_text,*) 'outflow from reservoirs: ', global_outflow, ' m3/s'
        CALL message (routine, message_text)
        WRITE (message_text,*) 'global discharge (land): ', global_discharge, ' m3/s'
        CALL message (routine, message_text)
        WRITE (message_text,*) 'discharge to the ocean : ', global_discharge_ocean, ' m3/s'
        CALL message (routine, message_text)
        WRITE (message_text,*) 'internal drainage      : ', global_internal_drain, ' m3/s'
        CALL message (routine, message_text)
        WRITE (message_text,*) 'Water budget error     : ', &
          & global_outflow - global_discharge - global_discharge_ocean - global_internal_drain - &
          & global_discharge_lake, ' m3/s'
        CALL message (routine, message_text)
      END IF

      ! distribute internal drainage to ocean inflow cells (weighted with ocean inflow)
      IF (global_discharge_ocean > 0._wp) THEN
        discharge_ocean(:,:) = discharge_ocean(:,:) + &
          & (discharge_ocean(:,:) * global_internal_drain / global_discharge_ocean)
      ELSE IF (global_internal_drain > 0._wp) THEN
        WRITE (message_text,*) 'no discharge to ocean! Internal drainage of ', &
          & global_internal_drain, ' m3/s will be lost.'
        CALL finish (routine, message_text)
      END IF

      IF (debug_hd) THEN
        ! water conservation test
        CALL message (routine, 'check re-distribution of water from internal drainage')
        global_discharge_ocean = global_sum_array(discharge_ocean(:,:)*in_domain(:,:))
        WRITE (message_text,*) 'discharge to ocean now : ', global_discharge_ocean, ' m3/s'
        CALL message (routine, message_text)
      END IF

    CASE ('weighted_to_coast')

      ! Dummy discharge: the discharge is distributed to all coastal ocean cells, weighted by cosine of latitude

      global_discharge = global_sum_array(discharge(:,:) * in_domain(:,:))

      ALLOCATE(cos_lat(SIZE(is_in_domain,1),SIZE(is_in_domain,2)))
      cos_lat(:,:) = MERGE(COS(deg2rad * grid%lat(:,:)), 0._wp, coast_ocean(:,:) > 0.5_wp)  ! coast_ocean is either 0. or 1.
      sum_weights = global_sum_array(coast_ocean(:,:) * cos_lat(:,:) * in_domain(:,:))
      WHERE (coast_ocean(:,:) > 0.5_wp)
        discharge_ocean(:,:) = global_discharge * cos_lat(:,:) / sum_weights
      ELSEWHERE
        discharge_ocean(:,:) = 0._wp
      END WHERE
      ! discharge(:,:) = 0._wp

    CASE ('zero')

      discharge_ocean(:,:) = 0._wp

    CASE DEFAULT
      WRITE (message_text,*)  'no valid parameter for routing_scheme: ', TRIM(routing_scheme)
      CALL finish(TRIM(routine), message_text)

    END SELECT

    DEALLOCATE (in_domain)

  END SUBROUTINE hd_lateral_flow


  !------------------------------------------------------------------------------------------------
  SUBROUTINE lake_lateral_flow(tile, options)
  !------------------------------------------------------------------------------------------------
  ! The routine is called from jsbach_finish_timestep, after the loop over the
  ! nproma blocks.
  ! It is thus possible to synchronize the halos. And to use data from
  ! surounding grid cells.
  !------------------------------------------------------------------------------------------------

    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile
    TYPE(t_jsb_task_options),   INTENT(in) :: options

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HD_)

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':lake_lateral_flow'

    ! Pointers to variables in memory

    dsl4jsb_Real2D_onDomain :: &
      & volume_lake_max,       &
      & basin_mask,            &
      & outlets,               &
      & volume_lake,           &
      & water_budget_lake,     &
      & outflow_outlet,        &
      & height_lake,           &
      & height_dyn,            &
      & inflow_lake,           &
      & outflow_lake,          &
      & overflow_lake,         &
      & inflow_outlet

    dsl4jsb_Real3D_onDomain :: &
      & outflow_min

    REAL(wp) ::                &
      & height_lake_diff,      &
      & u_lake,                &
      & k_lake,                &
      & global_lake_volume,    &  ! water budget change within the time step
      & global_lake_volume_error,  &  ! water budget change within the time step
      & global_lateral_flow_error, &  ! water budget change within the time step
      & global_lateral_flow_abs, &  ! water budget change within the time step
      & outlets_maxs,          &
      & neigh_dist,            &
      & steplen

    INTEGER  :: jc, jcs, jce       ! grid cells
    INTEGER  :: blk, blks, blke    ! blocks
    INTEGER  :: njc(3), nblk(3)    ! beighbor grid cell block and cell index
    INTEGER  :: n, m, l, nm, nl, nnm, ncount ! count neighbors
    INTEGER  :: nnjc(3), nnblk(3)
    INTEGER  :: nproma, nblks

    CHARACTER(len=32) :: &
      & routing_scheme        ! river routing scheme

    LOGICAL,  POINTER ::    &
      & is_in_domain (:,:)      ! T: cell in domain (not halo)

    REAL(wp),  POINTER ::    &
      & lat (:,:),           & ! latitude in degree
      & lon (:,:)              ! longitude in degree

    LOGICAL, ALLOCATABLE  ::  &
      & neigh_mask(:,:),       & ! mask for neighbor cells 
      & nneigh_mask(:,:)         ! mask for neighbor's neighbor cells 

    REAL(wp), ALLOCATABLE  ::  &
      & in_domain       (:,:), &  ! 1: cell in domain, 0: halo cell
      & volume_from_lake(:,:), &
      & outflow_mins(:,:)

    REAL(wp) :: outlets_max(p_patch(1)%nblks_c)

    LOGICAL  :: debug_hd

    IF (options%nc > 0) CONTINUE ! only to avoid compiler warnings

    model => Get_model(tile%owner_model_id)
    grid  => Get_grid(model%grid_id)
    is_in_domain => grid%patch%cells%decomp_info%owner_mask  (:,:)
    lat  => grid%lat(:,:)
    lon  => grid%lon(:,:)

    nproma  = grid%nproma
    nblks   = grid%nblks
    steplen = options%steplen

    dsl4jsb_Get_config(HD_)

    IF (.NOT. dsl4jsb_Config(HD_)%active) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    IF (ASSOCIATED(tile%parent)) &
      & CALL finish(TRIM(routine), 'HD model works on root tile only!')

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HD_)

    dsl4jsb_Get_var2D_onDomain(HD_, volume_lake_max)  ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, basin_mask)       ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, outlets)          ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, volume_lake)      ! INOUT
    dsl4jsb_Get_var2D_onDomain(HD_, water_budget_lake)! OUT
    dsl4jsb_Get_var2D_onDomain(HD_, outflow_outlet)   ! INOUT
    dsl4jsb_Get_var2D_onDomain(HD_, height_lake)      ! OUT
    dsl4jsb_Get_var2D_onDomain(HD_, height_dyn)       ! OUT
    dsl4jsb_Get_var2D_onDomain(HD_, inflow_lake)      ! OUT
    dsl4jsb_Get_var2D_onDomain(HD_, outflow_lake)     ! OUT
    dsl4jsb_Get_var2D_onDomain(HD_, inflow_outlet)    ! OUT
    dsl4jsb_Get_var2D_onDomain(HD_, overflow_lake)    ! OUT
    dsl4jsb_Get_var3D_onDomain(HD_, outflow_min)      ! IN

    ! Configuration parameters from the namelist
    routing_scheme = dsl4jsb_Config(HD_)%routing_scheme
    debug_hd       = dsl4jsb_Config(HD_)%debug_hd

    ! Mask out all halo cells for the global sums (otherwise these cells are
    ! counted twice)
    ALLOCATE (in_domain(nproma, nblks))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSE WHERE
      in_domain = 0._wp
    END WHERE

    SELECT CASE (TRIM(routing_scheme))
    CASE ('full')

      ! Synchronize halos
      CALL sync_patch_array(sync_c, grid%patch, basin_mask(:,:))
      CALL sync_patch_array(sync_c, grid%patch, height_dyn(:,:))
      CALL sync_patch_array(sync_c, grid%patch, height_lake(:,:))
      CALL sync_patch_array(sync_c, grid%patch, volume_lake(:,:))
      CALL sync_patch_array(sync_c, grid%patch, outflow_min(:,:,:))



      ALLOCATE (neigh_mask(nproma, nblks))
      ALLOCATE (nneigh_mask(nproma, nblks))
      ALLOCATE (outflow_mins(nproma, nblks))

      outflow_lake(:,:) = 0._wp
      inflow_lake(:,:) = 0._wp

      ! update dynamic lake water flow between lake grid cells

      ! block index (blk)
      blks = grid%get_blk_start()
      blke = grid%get_blk_end()
      DO blk = blks, blke

        ! grid cell index (jc)
        jcs = grid%get_col_start(blk)
        jce = grid%get_col_end(blk)
        DO jc = jcs, jce

          ! 3 neighbors cells: block and cell index 
          njc =  grid%patch%cells%neighbor_idx(jc,blk,:)
          nblk = grid%patch%cells%neighbor_blk(jc,blk,:)
          neigh_mask(:,:) = .FALSE.
          ! mask neighbors grid cells of the current basin with .TRUE.
          DO nm = 1, SIZE(njc)
            IF (basin_mask(jc,blk)==basin_mask(njc(nm),nblk(nm))) THEN
              neigh_mask(njc(nm),nblk(nm)) = .TRUE.
            END IF
          END DO          

          DO n = 1, SIZE(njc)

            IF(basin_mask(jc,blk) > 0._wp                       .AND. &    
               basin_mask(jc,blk) == basin_mask(njc(n),nblk(n))) THEN
          
              ! neighbor's 3 neighbor cells: block and cell index 
              nnjc =  grid%patch%cells%neighbor_idx(njc(n),nblk(n),:)
              nnblk = grid%patch%cells%neighbor_blk(njc(n),nblk(n),:)
              nneigh_mask(:,:) = .FALSE.
              ! mask neighbor's neighbor grid cells of the current basin with .TRUE.
              DO nnm = 1, SIZE(nnjc)
                IF (basin_mask(jc,blk)==basin_mask(nnjc(nnm),nnblk(nnm))) THEN
                  nneigh_mask(nnjc(nnm),nnblk(nnm)) = .TRUE.
                END IF
              END DO
              ! neighbor cell minimum outflow height
              outflow_mins(:,:) = 0._wp
              l = FINDLOC(grid%patch%cells%neighbor_idx(njc(n),nblk(n),:), jc, DIM=1)
              IF (outflow_min(njc(n),l,nblk(n)) /= 0._wp) THEN
                 outflow_mins(njc(n),nblk(n))= outflow_min(njc(n),l,nblk(n))
                 ! neighbor's neighbor cells minimum outflow height
                 DO nnm = 1, SIZE(nnjc)
                   nl =  FINDLOC(grid%patch%cells%neighbor_idx(nnjc(nnm),nnblk(nnm),:), njc(n), DIM=1)
                   IF (outflow_min(nnjc(nnm),nl,nnblk(nnm)) /= 0._wp) THEN
                     outflow_mins(nnjc(nnm),nnblk(nnm)) = outflow_min(nnjc(nnm),nl,nnblk(nnm))
                   END IF
                 END DO
              END IF

              ! inflow 
              IF(height_dyn(jc,blk) == MINVAL(PACK(height_dyn(:,:),nneigh_mask(:,:))) .AND. &
                 height_dyn(jc,blk) < height_dyn(njc(n),nblk(n))                      .AND. &
                 height_dyn(njc(n),nblk(n)) > outflow_mins(njc(n),nblk(n))            .AND. &
                 volume_lake(njc(n),nblk(n)) > 0._wp)                                  THEN
                height_lake_diff = MIN(height_lake(njc(n),nblk(n)), height_dyn(njc(n),nblk(n)) - height_dyn(jc,blk))
                CALL haversine_distance(lat(jc,blk),lon(jc,blk),lat(njc(n),nblk(n)),lon(njc(n),nblk(n)), neigh_dist)
                u_lake           = 100._wp * (height_lake_diff/neigh_dist)**0.5_wp * (height_lake_diff)**(2._wp/3._wp)
                k_lake           = neigh_dist/u_lake
                ncount = COUNT(height_dyn(jc,blk)   == PACK(height_dyn(:,:)  ,nneigh_mask(:,:)) .AND. &
                               outflow_mins(jc,blk) >= PACK(outflow_mins(:,:),nneigh_mask(:,:)))
                inflow_lake(jc,blk) = inflow_lake(jc,blk) + MIN( volume_lake(njc(n),nblk(n))/k_lake, &
                                                                 volume_lake(njc(n),nblk(n))/steplen) / ncount
              END IF
                
              ! outflow
              IF ( height_dyn(njc(n),nblk(n)) == MINVAL(PACK(height_dyn(:,:),neigh_mask(:,:))) .AND. &
                   height_dyn(njc(n),nblk(n)) < height_dyn(jc,blk)                             .AND. & 
                   height_dyn(jc,blk) > outflow_mins(jc,blk)                                   .AND. &
                   volume_lake(jc,blk)> 0._wp)                                                  THEN               
                height_lake_diff = MIN(height_lake(jc,blk), height_dyn(jc,blk) - height_dyn(njc(n),nblk(n)))
                CALL haversine_distance(lat(jc,blk),lon(jc,blk),lat(njc(n),nblk(n)),lon(njc(n),nblk(n)), neigh_dist)
                u_lake           = 100._wp * (height_lake_diff/neigh_dist)**0.5_wp * (height_lake_diff)**(2._wp/3._wp)
                k_lake           = neigh_dist/u_lake
                ncount = COUNT(height_dyn(njc(n),nblk(n))   == PACK(height_dyn(:,:)  ,neigh_mask(:,:)) .AND. &
                               outflow_mins(njc(n),nblk(n)) >= PACK(outflow_mins(:,:),neigh_mask(:,:)))
                outflow_lake(jc,blk) = outflow_lake(jc,blk) - MIN(volume_lake(jc,blk)/k_lake, &
                                                                  volume_lake(jc,blk)/steplen) / ncount
              END IF

            END IF
          END DO
        END DO
      END DO

      ! checkout water budget conservation: lateral flow only
      water_budget_lake(:,:) = water_budget_lake(:,:) + outflow_lake(:,:) + inflow_lake(:,:)
      volume_lake(:,:) = volume_lake(:,:) + ( outflow_lake(:,:) + inflow_lake(:,:) ) * steplen
      global_lateral_flow_error = global_sum_array((outflow_lake(:,:)+inflow_lake(:,:))*in_domain(:,:))
      global_lateral_flow_abs   = global_sum_array((ABS(outflow_lake(:,:))+ABS(inflow_lake(:,:)))*in_domain(:,:))
      IF (ABS(global_lateral_flow_error)/MAX(EPSILON(1._wp),global_lateral_flow_abs) > 10._wp*EPSILON(1._wp)) THEN
        WRITE (message_text,*) 'Lake water conservation problem: lateral flow imbalance: ', &
             global_lateral_flow_error,' m3'
        CALL finish (TRIM(routine), message_text)
      END IF


      ! outlet overflow to riverflow
      ALLOCATE (volume_from_lake(nproma, nblks))

      inflow_outlet(:,:)=0._wp
      inflow_outlet(:,:)=outflow_outlet
      overflow_lake(:,:)=0._wp

      DO blk = 1, p_patch(1)%nblks_c ! block centers
        outlets_max(blk) = MAXVAL(outlets(:,blk))
      END DO
      ! number of outlets
      global_lake_volume = global_sum_array(volume_lake(:,:)*in_domain(:,:))
      outlets_maxs = global_max(MAXVAL(outlets_max))  ! Get max over all PEs
      DO n = 1, outlets_maxs
        volume_from_lake(:,:) = 0._wp
        WHERE (n==basin_mask(:,:).AND.volume_lake(:,:)>volume_lake_max(:,:))
          volume_from_lake(:,:) = volume_lake(:,:) - volume_lake_max(:,:)
          volume_lake(:,:)      = volume_lake(:,:) - volume_from_lake(:,:)
          overflow_lake(:,:)    = (overflow_lake(:,:) - volume_from_lake(:,:)) / steplen
        END WHERE
        WHERE(n==outlets)
          outflow_outlet(:,:)=global_sum_array(volume_from_lake(:,:)/steplen*in_domain(:,:))
        END WHERE
      END DO

      water_budget_lake(:,:)= water_budget_lake(:,:) + overflow_lake(:,:)

      ! checkout water budget conservation: outlet overflow
      global_lake_volume_error=global_lake_volume-global_sum_array((volume_lake(:,:)+outflow_outlet(:,:)*steplen)*in_domain(:,:))
      IF (ABS(global_lake_volume_error)/MAX(global_lake_volume, EPSILON(1._wp)) > 10._wp*EPSILON(1._wp)) THEN
        WRITE (message_text,*) 'Lake water conservation problem: outlet overflow imbalance: ', &
             global_lake_volume_error,' m3'
        CALL finish (TRIM(routine), message_text)
      END IF
 
    CASE ('weighted_to_coast')
      WRITE (message_text,*)  'no valid parameter for dynlake: ', TRIM(routing_scheme)
      CALL finish(TRIM(routine), message_text)

    CASE ('zero')
      WRITE (message_text,*)  'no valid parameter for dynlake: ', TRIM(routing_scheme)
      CALL finish(TRIM(routine), message_text)

    CASE DEFAULT
      WRITE (message_text,*)  'no valid parameter for dynlake: ', TRIM(routing_scheme)
      CALL finish(TRIM(routine), message_text)

    END SELECT

  END SUBROUTINE lake_lateral_flow


  SUBROUTINE haversine_distance(clat, clon, nlat, nlon, dist)

    REAL(wp), INTENT(in) :: &
      & clat,               &
      & clon,               &
      & nlat,               &
      & nlon

    REAL(wp), INTENT(out) :: &
      & dist

    REAL(wp) :: r, a, delta_phi, delta_lambda

    r = 6371220.0_wp ! Reference Earth's Radius (m)

    delta_phi = (nlat - clat) * deg2rad
    delta_lambda = (nlon - clon) * deg2rad
    a = SIN(delta_phi / 2)**2 + COS(clat*deg2rad) * COS(nlat*deg2rad) * SIN(delta_lambda / 2)**2
    dist = r * (2 * ATAN2(SQRT(a), SQRT(1 - a)))

  END SUBROUTINE haversine_distance

  !------------------------------------------------------------------------------------------------
  SUBROUTINE hd_check_water_budget(tile, options)
  !------------------------------------------------------------------------------------------------
  ! calculation of the global water budget at the end of a time step, and compare it to the
  ! budget from the previous time step.
  !------------------------------------------------------------------------------------------------

    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile
    TYPE(t_jsb_task_options),   INTENT(in) :: options

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER    :: grid

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HD_)
    dsl4jsb_Def_memory(HYDRO_)

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':hd_check_water_budget'

    INTEGER :: nres_o_max, nres_b_max, nres_r_max
    INTEGER :: nproma, nblks

    REAL(wp), POINTER ::       &
      & area            (:,:)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: &
      & runoff          , &
      & drainage        , &
      & discharge       , &
      & inflow_outlet,    &
      & discharge_ocean , &
      & discharge_lake,   &
      & evapo_river,      &
      & evapo_overl,      &
      & water_budget    , &
      & water_budget_change

    dsl4jsb_Real3D_onDomain :: &
      & overlflow_res ,        &
      & baseflow_res  ,        &
      & riverflow_res

    LOGICAL,  POINTER ::    &
      & is_in_domain (:,:)      ! T: cell in domain (not halo)

    REAL(wp), ALLOCATABLE  ::     &
      & tile_fract       (:,:),   &
      & in_domain        (:,:),   &  ! 1: cell in domain, 0: halo cell
      & water_budget_old (:,:),   &  ! water budget of previous time step
      & water_flux       (:,:)       ! water fluxes within the time step

    REAL(wp) ::                   &
      & steplen,                  &
      & global_water_budget,      &  ! water budget of current step
      & global_water_budget_old,  &  ! water budget of previous step
      & global_water_error,       &  ! water budget change within the time step
      & global_water_flux            ! water fluxes changing the budget

    REAL(wp), POINTER :: water_error_gsum(:)

    LOGICAL :: debug_hd, diag_water_budget

    model => Get_model(tile%owner_model_id)

    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(HD_)

    IF (.NOT. model%Is_process_enabled(HD_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    IF (ASSOCIATED(tile%parent)) &
      & CALL finish(TRIM(routine), 'HD model works on root tile only!')

    steplen = options%steplen

    dsl4jsb_Get_memory(HD_)
    dsl4jsb_Get_memory(HYDRO_)

    ! Configuration parameters from the namelist
    debug_hd          = dsl4jsb_Config(HD_)%debug_hd
    diag_water_budget = dsl4jsb_Config(HD_)%diag_water_budget

    ! Get reference to variables in memory

    nres_o_max = hd__mem%nres_o_max
    nres_b_max = hd__mem%nres_b_max
    nres_r_max = hd__mem%nres_r_max

    nproma          = grid%nproma
    nblks           = grid%nblks

    area            => grid%area                                (:,:)
    is_in_domain    => grid%patch%cells%decomp_info%owner_mask  (:,:)

    dsl4jsb_Get_var3D_onDomain(HD_, overlflow_res) ! IN
    dsl4jsb_Get_var3D_onDomain(HD_, baseflow_res) ! IN
    dsl4jsb_Get_var3D_onDomain(HD_, riverflow_res) ! IN
    dsl4jsb_Get_var2D_onDomain(HD_, inflow_outlet)  ! IN

    dsl4jsb_Get_var2D_onDomain(HYDRO_, runoff)          ! IN
    dsl4jsb_Get_var2D_onDomain(HYDRO_, drainage)        ! IN
    dsl4jsb_Get_var2D_onDomain(HYDRO_, discharge)       ! IN
    dsl4jsb_Get_var2D_onDomain(HYDRO_, discharge_ocean) ! IN
    dsl4jsb_Get_var2D_onDomain(HYDRO_, evapo_river)     ! IN
    dsl4jsb_Get_var2D_onDomain(HYDRO_, evapo_overl)     ! IN

    dsl4jsb_Get_var2D_onDomain(HD_, discharge_lake)     ! in
    dsl4jsb_Get_var2D_onDomain(HD_, water_budget)       ! out
    dsl4jsb_Get_var2D_onDomain(HD_, water_budget_change)! out

    water_error_gsum => HD__mem%water_error_gsum%ptr(:) ! out

    ! Mask out all halo cells for the global sums (otherwise these cells are counted twice)
    ALLOCATE (in_domain(nproma, nblks))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSE WHERE
      in_domain = 0._wp
    END WHERE

    IF (TRIM(dsl4jsb_Config(HD_)%routing_scheme) /= 'full') THEN
      IF (diag_water_budget) THEN
        WRITE (message_text,*) '  global discharge to the ocean: ', global_sum_array(discharge_ocean*in_domain), ' m3/s'
        CALL message (TRIM(routine),message_text)
      END IF
      RETURN
    END IF

    ALLOCATE (water_budget_old(nproma, nblks))
    water_budget_old(:,:) = water_budget(:,:) * in_domain(:,:)
    global_water_budget_old = global_sum_array(water_budget_old(:,:))

    ! Calculate global land water budget
    !
    water_budget(:,:) = &  ! [m3]
      &   SUM(overlflow_res(:,:,:), DIM=2) * in_domain(:,:) & !< water in overland flow reservoirs
      & + SUM(baseflow_res (:,:,:), DIM=2) * in_domain(:,:) & !< water in baseflow reservoir
      & + SUM(riverflow_res(:,:,:), DIM=2) * in_domain(:,:) & !< water in riverflow reservoirs
      & + discharge(:,:) * in_domain(:,:) * steplen           !< river discharge [m3/s -> m3]
    global_water_budget = global_sum_array(water_budget(:,:))

    IF (debug_hd) THEN
      CALL message(routine, 'global water budget:')
      WRITE(message_text,*) '-   overlandflow reservoir: ', &
           global_sum_array(overlflow_res*SPREAD(in_domain,DIM=2,NCOPIES=nres_o_max)), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   baseflow reservoir    : ', &
           global_sum_array(baseflow_res*SPREAD(in_domain,DIM=2,NCOPIES=nres_b_max)), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   riverflow reservoir   : ', &
           global_sum_array(riverflow_res*SPREAD(in_domain,DIM=2,NCOPIES=nres_r_max)), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   discharge             : ', &
           global_sum_array(discharge*in_domain)*steplen, ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-->          total budget : ', global_water_budget, ' m3'
      CALL message(routine, message_text)
    END IF

    ALLOCATE(tile_fract(nproma, nblks))
    CALL tile%Get_fraction(fract=tile_fract(:,:))
    ALLOCATE(water_flux(nproma, nblks))
    water_flux(:,:) = &
      &   runoff         (:,:) * area(:,:)/rhoh2o*tile_fract(:,:)*in_domain(:,:)  & !< runoff   [kg m-2 s-1] -> [m3/s]
      & + drainage       (:,:) * area(:,:)/rhoh2o*tile_fract(:,:)*in_domain(:,:)  & !< drainage [kg m-2 s-1] -> [m3/s]
      & + evapo_river    (:,:) * area(:,:)/steplen*in_domain(:,:)                 & !< river evaporation [m] -> [m3/s] 
      & + evapo_overl    (:,:) * area(:,:)/steplen*in_domain(:,:)                 & !< river evaporation [m] -> [m3/s] 
      & + inflow_outlet  (:,:) * in_domain(:,:)                                   & !< drainge from the lakes (outlet overflow)
      & - discharge_ocean(:,:) * in_domain(:,:)                                   & !< discharge to the ocean
      & - discharge_lake (:,:) * in_domain(:,:)                                     !< discharge to the lakes

    global_water_flux = global_sum_array(water_flux(:,:))

    IF (debug_hd) THEN
      CALL message(routine, 'global water fluxes during time interval:')
      WRITE(message_text,*) '-   runoff (incl. P-E on glacier/lakes) : ', &
           global_sum_array(runoff/rhoh2o*area*tile_fract*in_domain*steplen), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   drainage              : ', &
           global_sum_array(drainage/rhoh2o*area*tile_fract*in_domain*steplen), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   river evaporation     : ', &
           global_sum_array(evapo_river*area*in_domain), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   overfland evaporation : ', &
           global_sum_array(evapo_overl*area*in_domain), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   outlet                : ', &
           global_sum_array(inflow_outlet*in_domain*steplen), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   discharge_ocean (neg.): ', &
           global_sum_array(discharge_ocean*in_domain*steplen), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   discharge_lake (neg.) : ', &
           global_sum_array(discharge_lake*in_domain*steplen), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-->      sum of the fluxes: ', global_water_flux*steplen, ' m3'
      CALL message(routine, message_text)
    END IF

    IF (diag_water_budget) THEN
      WRITE (message_text,*) '  global discharge to the ocean: ', global_sum_array(discharge_ocean*in_domain), ' m3/s'
      CALL message (TRIM(routine),message_text)
    END IF

    IF (is_time_experiment_start(options%current_datetime)) RETURN

    !
    ! 3. Global water balance check
    !
    water_budget_change(:,:) = water_budget_old(:,:) + water_flux(:,:) * steplen - water_budget(:,:)
    global_water_error = global_sum_array(water_budget_change(:,:))

    IF (debug_hd) THEN
      WRITE (message_text,*) 'Water budget imbalance during time step: ', global_water_error,' m3'
      CALL message (routine,message_text)
    END IF

    IF (global_water_budget_old /= 0._wp .AND. global_water_budget /= 0._wp) THEN        !vg: zero at the very first time step
      IF (ABS(global_water_error/global_water_budget) > 10._wp*EPSILON(1._wp)) THEN
        WRITE (message_text,*) 'Water conservation problem: budget imbalance: ', &
             global_water_error,' m3'
        IF (dsl4jsb_Config(HD_)%enforce_water_budget) THEN
          CALL finish (TRIM(routine), message_text)
        ELSE
          CALL message (TRIM(routine), message_text)
        END IF
      END IF
    END IF

    ! No unit transformation [m3]: 1
    IF (HD__mem%water_error_gsum%is_in_output) &
      & water_error_gsum = global_water_error

    DEALLOCATE (in_domain, tile_fract, water_budget_old, water_flux)

  END SUBROUTINE hd_check_water_budget

#endif
END MODULE mo_hd_interface
