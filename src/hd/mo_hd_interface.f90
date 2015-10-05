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

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish, message, message_text
  USE mo_sync,            ONLY: SYNC_C, sync_patch_array, global_sum_array

  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_class,       ONLY: get_model
  USE mo_jsb_grid_class,  ONLY: t_jsb_grid
  USE mo_jsb_grid,        ONLY: Get_grid
  USE mo_jsb_tile_class,  ONLY: t_jsb_tile
  USE mo_jsb_control,     ONLY: debug_on

  dsl4jsb_Use_processes HD, SRF, SOIL
  dsl4jsb_Use_config(HD)
  dsl4jsb_Use_memory(HD)
  dsl4jsb_Use_memory(SRF)
  dsl4jsb_Use_memory(SOIL)

  USE mo_hd_reservoir_cascade, ONLY: reservoir_cascade

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_hd, hd_riverflow, hd_init_water_budget, hd_check_water_budget

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_interface'

  REAL(wp), SAVE  :: water_budget

CONTAINS

  !------------------------------------------------------------------------------------------------
  SUBROUTINE update_hd(tile, iblk, ics, ice, nc, steplen)
  !------------------------------------------------------------------------------------------------
  ! the routine is called within the block-loop from jsbach_interface. There is no information on 
  ! neighbouring grid cells (on other paches) availabe.
  ! We calculate here the local discharge of each grid cell.
  !------------------------------------------------------------------------------------------------

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc
    REAL(wp),          INTENT(in) :: steplen

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(HD)
    dsl4jsb_Def_memory(HD)
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: area
    dsl4jsb_Real2D_onChunk :: land
    dsl4jsb_Real2D_onChunk :: lake
    dsl4jsb_Real2D_onChunk :: glacier
    dsl4jsb_Real2D_onChunk :: runoff
    dsl4jsb_Real2D_onChunk :: drainage
    dsl4jsb_Real2D_onChunk :: pme_glacier
    dsl4jsb_Real2D_onChunk :: discharge

    dsl4jsb_Real2D_onChunk :: ret_overlflow
    dsl4jsb_Real2D_onChunk :: nres_overlflow
    dsl4jsb_Real2D_onChunk :: ret_baseflow
    dsl4jsb_Real2D_onChunk :: nres_baseflow
    dsl4jsb_Real2D_onChunk :: ret_riverflow
    dsl4jsb_Real2D_onChunk :: nres_riverflow

    dsl4jsb_Real3D_onChunk :: overlflow_res
    dsl4jsb_Real3D_onChunk :: baseflow_res
    dsl4jsb_Real3D_onChunk :: riverflow_res

    REAL(wp), ALLOCATABLE :: inflow(:)
    REAL(wp), ALLOCATABLE :: outflow_runoff(:)
    REAL(wp), ALLOCATABLE :: outflow_drainage(:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_hd'

    LOGICAL           :: debug_hd
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

    area           => grid               %area                (ics:ice, iblk)

    dsl4jsb_Get_var2D_onChunk(SRF,  land,        OUT)
    dsl4jsb_Get_var2D_onChunk(SRF,  lake,        OUT)
    dsl4jsb_Get_var2D_onChunk(SRF,  glacier,     OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, runoff,      OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, drainage,    OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, pme_glacier, OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, discharge,   OUT)

    dsl4jsb_Get_var3D_onChunk(HD, overlflow_res,  OUT)
    dsl4jsb_Get_var3D_onChunk(HD, baseflow_res,   OUT)
    dsl4jsb_Get_var3D_onChunk(HD, riverflow_res,  OUT)
    dsl4jsb_Get_var2D_onChunk(HD, ret_overlflow,  OUT)
    dsl4jsb_Get_var2D_onChunk(HD, nres_overlflow, OUT)
    dsl4jsb_Get_var2D_onChunk(HD, ret_baseflow,   OUT)
    dsl4jsb_Get_var2D_onChunk(HD, nres_baseflow,  OUT)
    dsl4jsb_Get_var2D_onChunk(HD, ret_riverflow,  OUT)
    dsl4jsb_Get_var2D_onChunk(HD, nres_riverflow, OUT)


    ! Configuration parameters from the namelist
    debug_hd          = dsl4jsb_Config(HD)%debug_hd
    routing_scheme    = dsl4jsb_Config(HD)%routing_scheme


    SELECT CASE (TRIM(routing_scheme))
    CASE ('full')
      CALL finish(TRIM(routine), 'full HD model not yet implemented')

      ALLOCATE(inflow(nc))
      ALLOCATE(outflow_runoff(nc))
      ALLOCATE(outflow_drainage(nc))

      inflow(:) = runoff(:) + pme_glacier(:)  !vg eigentlich [m3/s]
      CALL reservoir_cascade (inflow(:), nc, steplen, ret_overlflow(:), nres_overlflow(:), overlflow_res(:,:), outflow_runoff(:))

      inflow(:) = drainage(:)  !vg eigentlich [m3/s]
      CALL reservoir_cascade (inflow(:), nc, steplen, ret_baseflow(:), nres_baseflow(:), baseflow_res(:,:), outflow_drainage(:))

      inflow(:) = outflow_runoff(:) + outflow_drainage(:)
      CALL reservoir_cascade (inflow(:), nc, steplen, ret_riverflow(:), nres_riverflow(:), riverflow_res(:,:), discharge(:))

      DEALLOCATE (inflow, outflow_runoff, outflow_drainage)

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
  SUBROUTINE hd_riverflow(model_id)
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
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)

    CHARACTER(len=*),  PARAMETER  :: routine = TRIM(modname)//':hd_riverflow'

    ! Pointers to variables in memory
    REAL(wp), POINTER :: area(:,:)

    dsl4jsb_Real2D_onDomain :: &
      & land,                  &
      & coast_ocean,           &
      & ocean,                 &
      & discharge,             &
      & discharge_ocean

    REAL(wp) ::            &
      & global_discharge,  &  ! global sum of discharge to the ocean
      & area_coast_ocean,  &  ! total area of all coastal ocean grid cells
      & area_ocean,        &  ! total area of all ocean grid cells
      & sum_weights

    REAL(wp), ALLOCATABLE :: cos_lat(:,:) ! Cosines of latitudes of cell centers

    LOGICAL, POINTER  ::   &
      & is_in_domain    (:,:) ! T: cell in domain (not halo)

    CHARACTER(len=32) :: &
      & routing_scheme        ! river routing scheme

    REAL(wp), ALLOCATABLE  ::   &
      & in_domain       (:,:) ! 1: cell in domain, 0: halo cell 

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    model => get_model(model_id)
    tile => model%Get_top_tile()        ! HD works on the root tile
    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(HD)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)

    area => grid%area(:,:)

    dsl4jsb_Get_var2D_onDomain(SRF, land,        IN)
    dsl4jsb_Get_var2D_onDomain(SRF, coast_ocean, IN)
    dsl4jsb_Get_var2D_onDomain(SRF, ocean,       IN)

    dsl4jsb_Get_var2D_onDomain(SOIL, discharge,       OUT)
    dsl4jsb_Get_var2D_onDomain(SOIL, discharge_ocean, OUT)

    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)

    ! Configuration parameters from the namelist
    routing_scheme = dsl4jsb_Config(HD)%routing_scheme


    ! Dummy discharge: the discharge is equally distributed to either
    !     - all coastal ocean cells (routing_scheme = 'equally_to_coast')  or
    !     - all ocean cells (routing_scheme = 'equally_to_ocean')
    ! (routing_scheme = 'full' for real river routing not yet implemented)

    ! Mask out all halo cells for the global sums (otherwise these cells are counted twice)
    ALLOCATE (in_domain(SIZE(is_in_domain,1),SIZE(is_in_domain,2)))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSEWHERE
      in_domain = 0._wp
    END WHERE
    global_discharge = global_sum_array(discharge(:,:)*land(:,:)*area(:,:)*in_domain(:,:))
    area_coast_ocean = global_sum_array(coast_ocean(:,:)*area(:,:)*in_domain(:,:))
    area_ocean = global_sum_array(ocean(:,:)*area(:,:)*in_domain(:,:))

    SELECT CASE (TRIM(routing_scheme))
    CASE ('full')
      CALL finish(TRIM(routine), 'full HD model not yet implemented')

!vg   CALL river_routing(discharge) -> discharge_ocean

    CASE ('weighted_to_coast')
      ALLOCATE(cos_lat(SIZE(is_in_domain,1),SIZE(is_in_domain,2)))
      cos_lat(:,:) = MERGE(COS(grid%lat(:,:)), 0._wp, coast_ocean(:,:) == 1._wp)
      sum_weights = global_sum_array(coast_ocean(:,:)*cos_lat(:,:)*area(:,:)*in_domain(:,:))
      WHERE (coast_ocean(:,:) == 1._wp)
        discharge_ocean(:,:) = global_discharge * cos_lat(:,:) / sum_weights
      ELSEWHERE
        discharge_ocean(:,:) = 0._wp
      END WHERE
    CASE ('weighted_to_ocean')
      ALLOCATE(cos_lat(SIZE(is_in_domain,1),SIZE(is_in_domain,2)))
      cos_lat(:,:) = MERGE(COS(grid%lat(:,:)), 0._wp, ocean(:,:) == 1._wp)
      sum_weights = global_sum_array(ocean(:,:)*cos_lat(:,:)*area(:,:)*in_domain(:,:))
      WHERE (ocean(:,:) == 1._wp)
        discharge_ocean(:,:) = global_discharge * cos_lat(:,:) / sum_weights
      ELSEWHERE
        discharge_ocean(:,:) = 0._wp
      END WHERE
    CASE ('equally_to_coast')
      WHERE (coast_ocean(:,:) == 1._wp)
        discharge_ocean(:,:) = global_discharge/area_coast_ocean
      ELSEWHERE
        discharge_ocean(:,:) = 0._wp
      END WHERE

    CASE ('equally_to_ocean')
      WHERE (ocean(:,:) == 1._wp)
        discharge_ocean(:,:) = global_discharge/area_ocean
      ELSEWHERE
        discharge_ocean(:,:) = 0._wp
      END WHERE

    CASE DEFAULT
      WRITE (message_text,*)  'no valid parameter for routing_scheme: ', TRIM(routing_scheme)
      CALL finish(TRIM(routine), message_text)

    END SELECT

    DEALLOCATE (in_domain)

!vg     ! Synchronize halos
!vg     IF (PRESENT(subset_range)) THEN
!vg       IF (.NOT. subset_range%is_in_domain) CALL sync_patch_array(sync_c, patch, discharge_ocean)
!vg     END IF

  END SUBROUTINE hd_riverflow

  !------------------------------------------------------------------------------------------------
  SUBROUTINE hd_init_water_budget(model_id)
  !------------------------------------------------------------------------------------------------
  ! calculation of the water budget before riverflow calculations
  !------------------------------------------------------------------------------------------------

    INTEGER, INTENT(in) :: model_id

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER    :: grid
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_tile),  POINTER    :: tile

    dsl4jsb_Def_config(HD)
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)

    CHARACTER(len=*),  PARAMETER  :: routine = TRIM(modname)//':hd_init_water_budget'

    ! Pointers to variables in memory
    REAL(wp), POINTER ::    &
      & area         (:,:)

    dsl4jsb_Real2D_onDomain :: &
      & land,                  &
      & glacier,               &
      & runoff,                &
      & drainage,              &
      & pme_glacier,           &
      & discharge

    LOGICAL, POINTER  ::   &
      & is_in_domain    (:,:) ! T: cell in domain (not halo) 

    REAL(wp), ALLOCATABLE  ::   &
      & in_domain       (:,:) ! 1: cell in domain, 0: halo cell

    LOGICAL :: debug_hd

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    model => get_model(model_id)
    tile => model%Get_top_tile()        ! HD works on the root tile
    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(HD)

    ! Get reference to variables in memory
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)

    area           => grid                        %area        (:,:)
    is_in_domain   => grid%patch%cells%decomp_info%owner_mask  (:,:)

    dsl4jsb_Get_var2D_onDomain(SRF,land,    IN)
    dsl4jsb_Get_var2D_onDomain(SRF,glacier, IN)

    dsl4jsb_Get_var2D_onDomain(SOIL,runoff,      IN)
    dsl4jsb_Get_var2D_onDomain(SOIL,drainage,    IN)
    dsl4jsb_Get_var2D_onDomain(SOIL,pme_glacier, IN)
    dsl4jsb_Get_var2D_onDomain(SOIL,discharge,   IN)

    ! Configuration parameters from the namelist
    debug_hd = dsl4jsb_Config(HD)%debug_hd

    ! Mask out all halo cells for the global sums (otherwise these cells are counted twice)
    ALLOCATE (in_domain(SIZE(is_in_domain(:,1)),SIZE(is_in_domain(1,:))))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSEWHERE
      in_domain = 0._wp
    END WHERE
    !
    ! calculate global land water budget
    !
    water_budget = &
          global_sum_array(runoff*area*land*in_domain)                    & !< runoff
        + global_sum_array(drainage*area*land*in_domain)                  & !< drainage
        + global_sum_array(pme_glacier*glacier*area*land*in_domain)   !vg & !< P-E over glaciers
!vg         + global_sum_array(flfmem)                      & !< water in the overland flow reservoirs
!vg         + global_sum_array(fgmem)                       & !< water in the baseflow reservoir
!vg         + global_sum_array(frfmem)                      & !< water in the riverflow reservoirs
!vg         + global_sum_array(finfl)                         !< inflow data from the hd restart file

    DEALLOCATE (in_domain)

  END SUBROUTINE hd_init_water_budget

  !------------------------------------------------------------------------------------------------
  SUBROUTINE hd_check_water_budget(model_id)
  !------------------------------------------------------------------------------------------------
  ! calculation of the global water budget at the end of a time step, and compare it to the
  ! budget from the beginning of the time step.
  !------------------------------------------------------------------------------------------------

    INTEGER, INTENT(in) :: model_id

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER    :: grid
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_tile),  POINTER    :: tile

    dsl4jsb_Def_config(HD)
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)

    CHARACTER(len=*),  PARAMETER  :: routine = TRIM(modname)//':hd_check_water_budget'

    ! Pointers to variables in memory
    REAL(wp), POINTER ::       &
      & area            (:,:)

    dsl4jsb_Real2D_onDomain :: &
      & land            (:,:), &
      & discharge_ocean (:,:)

    LOGICAL,  POINTER ::    &
      & is_in_domain (:,:)      ! T: cell in domain (not halo)

    REAL(wp), ALLOCATABLE  ::   &
      & in_domain       (:,:) ! 1: cell in domain, 0: halo cell

    REAL(wp) ::             &
      & water_budget_old,   &   ! water budget at the beginning of the time step
      & water_budget_change     ! water budget change within the time step

    LOGICAL :: debug_hd, diag_water_budget

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    model => get_model(model_id)
    tile => model%Get_top_tile()        ! HD works on the root tile
    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(HD)

    ! Get reference to variables in memory
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)

    area            => grid         %area            (:,:)
    is_in_domain    => grid%patch%cells%decomp_info%owner_mask (:,:)

    dsl4jsb_Get_var2D_onDomain(SRF,land, IN)

    dsl4jsb_Get_var2D_onDomain(SOIL, discharge_ocean, IN)

    ! Configuration parameters from the namelist
    debug_hd          = dsl4jsb_Config(HD)%debug_hd
    diag_water_budget = dsl4jsb_Config(HD)%diag_water_budget

    ! Mask out all halo cells for the global sums (otherwise these cells are counted twice)
    ALLOCATE (in_domain(SIZE(is_in_domain(:,1)),SIZE(is_in_domain(1,:))))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSEWHERE
      in_domain = 0._wp
    END WHERE

    water_budget_old = water_budget

    ! re-calculate global land water budget
    !
    water_budget = &
           global_sum_array(discharge_ocean*area*in_domain)  !vg  & !< river discharge
!vg          + global_sum_array(flfmem)                    & !< water in the overland flow reservoirs
!vg          + global_sum_array(fgmem)                     & !< water in the baseflow reservoir
!vg          + global_sum_array(frfmem)                    & !< water in the riverflow reservoirs
!vg          + global_sum_array(finfl)                       !< inflow data

    ! water conservation test
    !
    water_budget_change = water_budget_old - water_budget
    IF (diag_water_budget) THEN
      WRITE (message_text,*) 'Water budget change: ', water_budget_change,' kg/s'
    END IF
    CALL message (TRIM(routine),message_text)
    IF (water_budget /= 0._wp) THEN
      IF (ABS(water_budget_change/water_budget) > 1.E-7_wp) THEN
        WRITE (message_text,*) 'Water conservation problem: budget change: ', &
             water_budget_change,' kg/s'
        CALL finish (TRIM(routine),message_text)
      END IF
    END IF
    IF (diag_water_budget) THEN
      WRITE (message_text,*) '  global discharge to the ocean: ', global_sum_array(discharge_ocean*area*in_domain), ' kg/s'
      CALL message (TRIM(routine),message_text)
    END IF

    DEALLOCATE (in_domain)

  END SUBROUTINE hd_check_water_budget


#endif
END MODULE mo_hd_interface
