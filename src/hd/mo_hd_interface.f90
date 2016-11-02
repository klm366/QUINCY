!>
!! @brief Contains the interfaces to the HD model
!!
!! @author
!!  Veronika Gayler, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Veronika Gayler (2013-04-20)
!!
#include "dsl_jsbach.inc"

MODULE mo_hd_interface
!#ifndef __NO_JSBACH__
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: finish, message, message_text
  USE mo_sync,             ONLY: SYNC_C, sync_patch_array, global_sum_array

  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_jsb_class,        ONLY: get_model
  USE mo_jsb_grid_class,   ONLY: t_jsb_grid
  USE mo_jsb_grid,         ONLY: Get_grid
  USE mo_jsb_domain_iface, ONLY: get_blk_start, get_blk_end, get_col_start, get_col_end
  USE mo_jsb_tile_class,   ONLY: t_jsb_tile
  USE mo_jsb_control,      ONLY: debug_on

  dsl4jsb_Use_processes HD, SRF, SOIL
  dsl4jsb_Use_config(HD)
  dsl4jsb_Use_memory(HD)
  dsl4jsb_Use_memory(SRF)
  dsl4jsb_Use_memory(SOIL)

  USE mo_hd_reservoir_cascade, ONLY: reservoir_cascade
  USE mo_srf_constants,   ONLY: rhoh2o

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_hd, hd_lateral_flow, hd_check_water_budget
  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_interface'

CONTAINS

  !------------------------------------------------------------------------------------------------
  SUBROUTINE update_hd(tile, iblk, ics, ice, nc, steplen)
  !------------------------------------------------------------------------------------------------
  ! the routine is called within the block-loop from jsbach_interface. There is no information on
  ! neighbouring grid cells (on other patches) availabe.
  ! We calculate here the local outflow of each grid cell.
  !------------------------------------------------------------------------------------------------

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc
    REAL(wp),          INTENT(in) :: steplen

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    dsl4jsb_Def_config(HD)
    dsl4jsb_Def_memory(HD)
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: land
    dsl4jsb_Real2D_onChunk :: lake
    dsl4jsb_Real2D_onChunk :: glacier
    dsl4jsb_Real2D_onChunk :: runoff
    dsl4jsb_Real2D_onChunk :: drainage
    dsl4jsb_Real2D_onChunk :: pme_glacier
    dsl4jsb_Real2D_onChunk :: discharge
    dsl4jsb_Real2D_onChunk :: outflow_runoff
    dsl4jsb_Real2D_onChunk :: outflow_drainage
    dsl4jsb_Real2D_onChunk :: outflow_rivers
    dsl4jsb_Real2D_onChunk :: ret_overlflow
    dsl4jsb_Real2D_onChunk :: nres_overlflow
    dsl4jsb_Real2D_onChunk :: ret_baseflow
    dsl4jsb_Real2D_onChunk :: nres_baseflow
    dsl4jsb_Real2D_onChunk :: ret_riverflow
    dsl4jsb_Real2D_onChunk :: nres_riverflow
    dsl4jsb_Real3D_onChunk :: overlflow_res
    dsl4jsb_Real3D_onChunk :: baseflow_res
    dsl4jsb_Real3D_onChunk :: riverflow_res

    REAL(wp), POINTER :: area(:)

    REAL(wp), ALLOCATABLE :: inflow(:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_hd'

    CHARACTER(len=32) :: routing_scheme

    IF (.NOT. tile%Is_process_active(HD)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    model => get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(HD)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HD)
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)

    area => grid %area (ics:ice, iblk)

    dsl4jsb_Get_var2D_onChunk(SRF,  land,        IN)
    dsl4jsb_Get_var2D_onChunk(SRF,  lake,        IN)
    dsl4jsb_Get_var2D_onChunk(SRF,  glacier,     IN)
    dsl4jsb_Get_var2D_onChunk(SOIL, runoff,      IN)
    dsl4jsb_Get_var2D_onChunk(SOIL, drainage,    IN)
    dsl4jsb_Get_var2D_onChunk(SOIL, pme_glacier, IN)
    dsl4jsb_Get_var2D_onChunk(SOIL, discharge,   INOUT)

    dsl4jsb_Get_var2D_onChunk(HD, ret_overlflow,    IN)
    dsl4jsb_Get_var2D_onChunk(HD, nres_overlflow,   IN)
    dsl4jsb_Get_var2D_onChunk(HD, ret_baseflow,     IN)
    dsl4jsb_Get_var2D_onChunk(HD, nres_baseflow,    IN)
    dsl4jsb_Get_var2D_onChunk(HD, ret_riverflow,    IN)
    dsl4jsb_Get_var2D_onChunk(HD, nres_riverflow,   IN)
    dsl4jsb_Get_var3D_onChunk(HD, overlflow_res,    INOUT)
    dsl4jsb_Get_var3D_onChunk(HD, baseflow_res,     INOUT)
    dsl4jsb_Get_var3D_onChunk(HD, riverflow_res,    INOUT)
    dsl4jsb_Get_var2D_onChunk(HD, outflow_runoff,   OUT)
    dsl4jsb_Get_var2D_onChunk(HD, outflow_drainage, OUT)
    dsl4jsb_Get_var2D_onChunk(HD, outflow_rivers,   OUT)


    ! Configuration parameters from the namelist
    routing_scheme    = dsl4jsb_Config(HD)%routing_scheme


    SELECT CASE (TRIM(routing_scheme))
    CASE ('full')

      ALLOCATE(inflow(nc))

      ! unit conversions: [kg m-2 s-1] -> [m3/s]
      inflow(:) = (runoff(:) + pme_glacier(:) * glacier) * land * area / rhoh2o
      CALL reservoir_cascade (inflow(:), nc, steplen, ret_overlflow(:), nres_overlflow(:), overlflow_res(:,:), outflow_runoff(:))

      inflow(:) = drainage(:) * land * area / rhoh2o
      CALL reservoir_cascade (inflow(:), nc, steplen, ret_baseflow(:),  nres_baseflow(:),  baseflow_res(:,:),  outflow_drainage(:))

      inflow(:) = discharge(:)
      CALL reservoir_cascade (inflow(:), nc, steplen, ret_riverflow(:), nres_riverflow(:), riverflow_res(:,:), outflow_rivers(:))

      DEALLOCATE (inflow)

    CASE ('weighted_to_coast', 'weighted_to_ocean', 'equally_to_coast', 'equally_to_ocean')
      !
      ! dummy version of HD model
      !
      discharge(:) = land(:) * (runoff(:) + drainage(:) + pme_glacier(:))

    CASE DEFAULT
      WRITE (message_text,*)  'no valid parameter for routing_scheme: ', TRIM(routing_scheme)
      CALL finish(TRIM(routine), message_text)

    END SELECT

  END SUBROUTINE update_hd

  !------------------------------------------------------------------------------------------------
  SUBROUTINE hd_lateral_flow(model_id)
  !------------------------------------------------------------------------------------------------
  ! The routine is called from jsbach_finish_timestep, after the loop over the nproma blocks.
  ! It is thus possible to synchronize the halos. And to use data from surounding grid cells.
  !------------------------------------------------------------------------------------------------

    INTEGER, INTENT(in) :: model_id

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER    :: grid
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_tile),  POINTER    :: tile

    dsl4jsb_Def_config(HD)
    dsl4jsb_Def_memory(HD)
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)

    INTEGER                       :: jc, blk, n
    INTEGER                       :: jc_us, blk_us

    CHARACTER(len=*),  PARAMETER  :: routine = TRIM(modname)//':hd_lateral_flow'

    ! Pointers to variables in memory
    REAL(wp), POINTER :: area(:,:)

    dsl4jsb_Real2D_onDomain :: &
      & land,                  &
      & hd_mask,               &
      & coast_ocean,           &
      & outflow_runoff,        &
      & outflow_drainage,      &
      & outflow_rivers,        &
      & outflow_resid,         &
      & discharge,             &
      & discharge_ocean,       &
      & internal_drain

    dsl4jsb_Real3D_onDomain :: &
      & nidx_upstream,         &
      & bidx_upstream

    REAL(wp) ::                 &
      & global_discharge,       &  ! global sum of discharge to the ocean
      & global_outflow,         &
      & global_outflow_resid,   &
      & global_internal_drain,  &
      & global_discharge_ocean, &
      & sum_weights

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

    LOGICAL  :: debug_hd

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    model => get_model(model_id)
    tile => model%Get_top_tile()        ! HD works on the root tile
    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(HD)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HD)
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)

    area   => grid%area(:,:)
    nproma =  grid%nproma
    nblks  =  grid%nblks

    dsl4jsb_Get_var2D_onDomain(HD, hd_mask,          IN)
    dsl4jsb_Get_var2D_onDomain(HD, outflow_runoff,   IN)
    dsl4jsb_Get_var2D_onDomain(HD, outflow_drainage, IN)
    dsl4jsb_Get_var2D_onDomain(HD, outflow_rivers,   IN)
    dsl4jsb_Get_var2D_onDomain(HD, outflow_resid,    OUT)
    dsl4jsb_Get_var3D_onDomain(HD, nidx_upstream,    IN)
    dsl4jsb_Get_var3D_onDomain(HD, bidx_upstream,    IN)

    dsl4jsb_Get_var2D_onDomain(SRF, land,        IN)
    dsl4jsb_Get_var2D_onDomain(SRF, coast_ocean, IN)

    dsl4jsb_Get_var2D_onDomain(SOIL, discharge,       OUT)
    dsl4jsb_Get_var2D_onDomain(SOIL, discharge_ocean, OUT)
    dsl4jsb_Get_var2D_onDomain(SOIL, internal_drain,  OUT)

    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)

    ! Configuration parameters from the namelist
    routing_scheme = dsl4jsb_Config(HD)%routing_scheme
    debug_hd       = dsl4jsb_Config(HD)%debug_hd


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

      global_outflow        = global_sum_array((outflow_runoff(:,:)+outflow_drainage(:,:)+outflow_rivers(:,:))*in_domain(:,:))

      blks = get_blk_start(grid%patch)
      blke = get_blk_end(grid%patch)
      DO blk = blks, blke
        jcs = get_col_start(blk, grid%patch)
        jce = get_col_end(blk, grid%patch)
        DO jc = jcs, jce
          ! inflow from grid cells upstream: outflow from runnoff, drainage and riverflow
          DO n = 1, hd__mem%nneigh
            ! local indices of grid cells upstream
            jc_us  =  nidx_upstream(jc,n,blk)
            blk_us =  bidx_upstream(jc,n,blk)
            IF (jc_us > 0 .and. blk_us > 0) THEN
              discharge(jc,blk) =  discharge(jc,blk)               &
                                 + outflow_runoff(jc_us,blk_us)    &
                                 + outflow_drainage(jc_us,blk_us)  &
                                 + outflow_rivers(jc_us,blk_us)
              ! Enable the next three lines for debugging ... but this only works on one (!) processor, i.e. no halos
              ! outflow_runoff  (jc_us,blk_us) = 0._wp
              ! outflow_drainage(jc_us,blk_us) = 0._wp
              ! outflow_rivers  (jc_us,blk_us) = 0._wp
            END IF
          END DO
        END DO
      END DO

      ! Enable this block for debugging ... but this only works on one (!) processor, i.e. no halos
      ! outflow_resid(:,:) = outflow_runoff(:,:) + outflow_drainage(:,:) + outflow_rivers(:,:)
      ! DO blk = blks, blke
      !   jcs = get_col_start(blk, grid%patch)
      !   jce = get_col_end(blk, grid%patch)
      !   DO jc = jcs, jce
      !     IF (outflow_resid(jc,blk) /= 0._wp) print*,'AAA ', jc,blk, outflow_resid(jc,blk)
      !   END DO
      ! END DO

      WHERE (hd_mask(:,:) == 0._wp)
        ! ocean inflow cells
        discharge_ocean(:,:) = discharge(:,:)
        discharge(:,:) = 0._wp
      ELSE WHERE (hd_mask(:,:) == 2._wp)
        ! internal drainage cells
        internal_drain(:,:) = discharge(:,:)
        discharge(:,:) = 0._wp
      END WHERE

      ! global integrals
      ! Switch the next two lines for debugging ... but this only works on one (!) processor, i.e. no halos
      ! global_outflow_resid  = global_sum_array(outflow_resid(:,:)*in_domain(:,:))
      global_outflow_resid   = 0._wp
      global_discharge       = global_sum_array(discharge(:,:)*in_domain(:,:))
      global_internal_drain  = global_sum_array(internal_drain(:,:)*in_domain(:,:))
      global_discharge_ocean = global_sum_array(discharge_ocean(:,:)*in_domain(:,:))

      IF (debug_hd) THEN
        ! water conservation test
        WRITE (message_text,*) 'outflow from reservoirs: ', global_outflow, ' m3/s'
        CALL message (routine, message_text)
        ! Enable the next two lines for debugging ... but this only works on one (!) processor, i.e. no halos
        ! WRITE (message_text,*) 'outflow residual: ', global_outflow_resid, ' m3/s'
        ! CALL message (routine, message_text)
        WRITE (message_text,*) 'global discharge (land): ', global_discharge, ' m3/s'
        CALL message (routine, message_text)
        WRITE (message_text,*) 'discharge to the ocean : ', global_discharge_ocean, ' m3/s'
        CALL message (routine, message_text)
        WRITE (message_text,*) 'internal drainage      : ', global_internal_drain, ' m3/s'
        CALL message (routine, message_text)
        WRITE (message_text,*) 'Water budget error     : ', &
          & global_outflow - global_outflow_resid - global_discharge - global_discharge_ocean - global_internal_drain, ' m3/s'
        CALL message (routine, message_text)
      END IF

      ! distribute internal drainage to ocean inflow cells (weighted with ocean inflow)
      IF (global_discharge_ocean > 0._wp) THEN
        discharge_ocean(:,:) = discharge_ocean(:,:) + &
          & (discharge_ocean(:,:) * (global_internal_drain+global_outflow_resid) / global_discharge_ocean)
      ELSE IF (global_internal_drain+global_outflow_resid > 0._wp) THEN
        WRITE (message_text,*) 'no discharge to ocean! Internal drainage of ', &
          & global_internal_drain+global_outflow_resid, ' m3/s will be lost.'
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

      ! Dummy discharge: the discharge is equally distributed to all coastal ocean cells

      global_discharge = global_sum_array(discharge(:,:)*land(:,:)*area(:,:)*in_domain(:,:))

      ALLOCATE(cos_lat(SIZE(is_in_domain,1),SIZE(is_in_domain,2)))
      cos_lat(:,:) = MERGE(COS(grid%lat(:,:)), 0._wp, coast_ocean(:,:) == 1._wp)
      sum_weights = global_sum_array(coast_ocean(:,:)*cos_lat(:,:)*area(:,:)*in_domain(:,:))
      WHERE (coast_ocean(:,:) == 1._wp)
        discharge_ocean(:,:) = global_discharge * cos_lat(:,:) / sum_weights
      ELSEWHERE
        discharge_ocean(:,:) = 0._wp
      END WHERE
      discharge(:,:) = 0._wp

    CASE DEFAULT
      WRITE (message_text,*)  'no valid parameter for routing_scheme: ', TRIM(routing_scheme)
      CALL finish(TRIM(routine), message_text)

    END SELECT

    DEALLOCATE (in_domain)

  END SUBROUTINE hd_lateral_flow

  !------------------------------------------------------------------------------------------------
  SUBROUTINE hd_check_water_budget(model_id, steplen)
  !------------------------------------------------------------------------------------------------
  ! calculation of the global water budget at the end of a time step, and compare it to the
  ! budget from the previous time step.
  !------------------------------------------------------------------------------------------------

    INTEGER,  INTENT(in) :: model_id
    REAL(wp), INTENT(in) :: steplen

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER    :: grid
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_tile),  POINTER    :: tile

    dsl4jsb_Def_config(HD)
    dsl4jsb_Def_memory(HD)
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)

    CHARACTER(len=*),  PARAMETER  :: routine = TRIM(modname)//':hd_check_water_budget'

    INTEGER :: nres_o_max, nres_b_max, nres_r_max
    INTEGER :: nproma, nblks

    REAL(wp), POINTER ::       &
      & area            (:,:)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: &
      & land            , &
      & glacier         , &
      & runoff          , &
      & drainage        , &
      & pme_glacier     , &
      & outflow_resid   , &
      & discharge       , &
      & discharge_ocean

    dsl4jsb_Real3D_onDomain :: &
      & overlflow_res ,        &
      & baseflow_res  ,        &
      & riverflow_res 

    LOGICAL,  POINTER ::    &
      & is_in_domain (:,:)      ! T: cell in domain (not halo)

    REAL(wp), ALLOCATABLE  ::   &
      & in_domain       (:,:) ! 1: cell in domain, 0: halo cell

    REAL(wp) ::              &
      & water_budget,        &  ! water budget of current step
      & water_budget_old,    &  ! water budget of previous step
      & water_budget_change, &  ! water budget change within the time step
      & water_fluxes            ! water fluxes changeing the budget

    LOGICAL :: debug_hd, diag_water_budget

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    model => get_model(model_id)
    tile => model%Get_top_tile()        ! HD works on the root tile
    grid => Get_grid(model%grid_id)


    dsl4jsb_Get_config(HD)

    dsl4jsb_Get_memory(HD)
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)

    ! Get reference to variables in memory

    nres_o_max = hd__mem%nres_o_max
    nres_b_max = hd__mem%nres_b_max
    nres_r_max = hd__mem%nres_r_max

    nproma          = grid%nproma
    nblks           = grid%nblks

    area            => grid%area                                (:,:)
    is_in_domain    => grid%patch%cells%decomp_info%owner_mask  (:,:)

    dsl4jsb_Get_var2D_onDomain(SRF, land,    IN)
    dsl4jsb_Get_var2D_onDomain(SRF, glacier, IN)

    dsl4jsb_Get_var3D_onDomain(HD, overlflow_res, IN)
    dsl4jsb_Get_var3D_onDomain(HD, baseflow_res,  IN)
    dsl4jsb_Get_var3D_onDomain(HD, riverflow_res, IN)

    dsl4jsb_Get_var2D_onDomain(SOIL, runoff,          IN)
    dsl4jsb_Get_var2D_onDomain(SOIL, drainage,        IN)
    dsl4jsb_Get_var2D_onDomain(SOIL, pme_glacier,     IN)
    dsl4jsb_Get_var2D_onDomain(HD,   outflow_resid,   IN)
    dsl4jsb_Get_var2D_onDomain(SOIL, discharge,       IN)
    dsl4jsb_Get_var2D_onDomain(SOIL, discharge_ocean, IN)

    water_budget   =  dsl4jsb_var2D_ptr(HD,water_budget)  (1,1)

    ! Configuration parameters from the namelist
    debug_hd          = dsl4jsb_Config(HD)%debug_hd
    diag_water_budget = dsl4jsb_Config(HD)%diag_water_budget

    ! Mask out all halo cells for the global sums (otherwise these cells are counted twice)
    ALLOCATE (in_domain(nproma, nblks))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSEWHERE
      in_domain = 0._wp
    END WHERE

    water_budget_old = water_budget

    ! re-calculate global land water budget
    !
    water_budget = &
           global_sum_array(overlflow_res*SPREAD(in_domain,DIM=2,NCOPIES=nres_o_max)) & !< water in overland flow reservoirs [m3]
         + global_sum_array(baseflow_res *SPREAD(in_domain,DIM=2,NCOPIES=nres_b_max)) & !< water in baseflow reservoir       [m3]
         + global_sum_array(riverflow_res*SPREAD(in_domain,DIM=2,NCOPIES=nres_r_max)) & !< water in riverflow reservoirs     [m3]
         + global_sum_array(discharge*in_domain) * steplen                              !< river discharge [m3/s -> m3]

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
      WRITE(message_text,*) '-->          total budget : ', water_budget, ' m3'
      CALL message(routine, message_text)
    END IF

    water_fluxes = &
           global_sum_array(runoff             *area/rhoh2o*land*in_domain)  & !< runoff                 [kg m-2 s-1] -> [m3/s]
         + global_sum_array(pme_glacier*glacier*area/rhoh2o*land*in_domain)  & !< P-E over glaciers      [kg m-2 s-1] -> [m3/s]
         + global_sum_array(drainage           *area/rhoh2o*land*in_domain)  & !< drainage               [kg m-2 s-1] -> [m3/s]
         - global_sum_array(discharge_ocean    *in_domain)                     !< discharge to the ocean

    IF (debug_hd) THEN
      CALL message(routine, 'global water fluxes in time interval:')
      WRITE(message_text,*) '-   runoff                : ', &
           global_sum_array(runoff/rhoh2o*area*land*in_domain*steplen), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   pme_glacier           : ', &
           global_sum_array(pme_glacier/rhoh2o*glacier*area*land*in_domain*steplen), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   drainage              : ', &
           global_sum_array(drainage/rhoh2o*area*land*in_domain*steplen), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-   discharge_ocean (neg.): ', &
           global_sum_array(discharge_ocean*in_domain*steplen), ' m3'
      CALL message(routine, message_text)
      WRITE(message_text,*) '-->      sum of the fluxes: ', water_fluxes*steplen, ' m3'
      CALL message(routine, message_text)
    END IF

    ! water conservation test
    !
    water_budget_change =  water_budget_old - water_budget + water_fluxes*steplen

    IF (debug_hd) THEN
      WRITE (message_text,*) 'Water budget change in time step: ', water_budget_change,' m3'
      CALL message (routine,message_text)
    END IF

    IF (water_budget_old /= 0._wp .AND. water_budget /= 0._wp) THEN        !vg: zero at the very first time step
      IF (ABS(water_budget_change/water_budget) > 10._wp*EPSILON(1._wp)) THEN
        WRITE (message_text,*) 'Water conservation problem: budget change: ', &
             water_budget_change,' m3'
        CALL finish (routine, message_text)
      END IF
    END IF
    IF (diag_water_budget) THEN
      WRITE (message_text,*) '  global discharge to the ocean: ', global_sum_array(discharge_ocean*in_domain), ' m3/s'
      CALL message (TRIM(routine),message_text)
    END IF

    dsl4jsb_var2D_ptr(HD,water_budget) = water_budget

    DEALLOCATE (in_domain)

  END SUBROUTINE hd_check_water_budget

#endif
END MODULE mo_hd_interface
