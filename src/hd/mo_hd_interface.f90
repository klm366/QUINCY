!>
!! @brief Contains the interfaces to the HD model
!!
!! @author
!!  Veronika Gayler, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Veronika Gayler (2013-04-20)
!!
MODULE mo_hd_interface
!#ifndef __NO_JSBACH__
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

#include "dsl_jsbach.inc"

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish, message, message_text
  USE mo_sync,            ONLY: SYNC_C, sync_patch_array, global_sum_array

  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_class,       ONLY: get_model
  USE mo_jsb_grid_class,  ONLY: t_jsb_grid
  USE mo_jsb_grid,        ONLY: Get_grid
  USE mo_jsb_tile_class,  ONLY: t_jsb_tile
  USE mo_jsb_control,     ONLY: debug_on

  dsl$jsb Use_processes HD, SRF, SOIL
  dsl$jsb Use_config(hd)
  dsl$jsb Use_memory(hd)
  dsl$jsb Use_memory(srf)
  dsl$jsb Use_memory(soil)

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

    dsl$jsb Def_config(hd_conf, hd)
    dsl$jsb Def_memory(hd_mem, hd)
    dsl$jsb Def_memory(srf_mem, srf)
    dsl$jsb Def_memory(soil_mem, soil)

    ! Pointers to variables in memory
    dsl$jsb onChunk_Real2D :: area
    dsl$jsb onChunk_Real2D :: land
    dsl$jsb onChunk_Real2D :: lake
    dsl$jsb onChunk_Real2D :: glacier
    dsl$jsb onChunk_Real2D :: runoff
    dsl$jsb onChunk_Real2D :: drainage
    dsl$jsb onChunk_Real2D :: pme_glacier
    dsl$jsb onChunk_Real2D :: discharge

    dsl$jsb onChunk_Real2D :: ret_overlflow
    dsl$jsb onChunk_Real2D :: nres_overlflow
    dsl$jsb onChunk_Real2D :: ret_baseflow
    dsl$jsb onChunk_Real2D :: nres_baseflow
    dsl$jsb onChunk_Real2D :: ret_riverflow
    dsl$jsb onChunk_Real2D :: nres_riverflow

    dsl$jsb onChunk_Real3D :: overlflow_res
    dsl$jsb onChunk_Real3D :: baseflow_res
    dsl$jsb onChunk_Real3D :: riverflow_res

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

    dsl$jsb Get_config(hd_conf, HD)

    ! Get reference to variables for current block
    !
    dsl$jsb Get_memory(hd_mem, HD)
    dsl$jsb Get_memory(srf_mem, SRF)
    dsl$jsb Get_memory(soil_mem, SOIL)

    area           => grid               %area                (ics:ice, iblk)
    dsl$jsb Get_var2D_onChunk( srf_mem,  land,           OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,  lake,           OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,  glacier,        OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem, runoff,         OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem, drainage,       OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem, pme_glacier,    OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem, discharge,      OUT)

    dsl$jsb Get_var3D_onChunk( hd_mem,   overlflow_res,  OUT)
    dsl$jsb Get_var3D_onChunk( hd_mem,   baseflow_res,   OUT)
    dsl$jsb Get_var3D_onChunk( hd_mem,   riverflow_res,  OUT)

    dsl$jsb Get_var2D_onChunk( hd_mem,   ret_overlflow,  OUT)
    dsl$jsb Get_var2D_onChunk( hd_mem,   nres_overlflow, OUT)
    dsl$jsb Get_var2D_onChunk( hd_mem,   ret_baseflow,   OUT)
    dsl$jsb Get_var2D_onChunk( hd_mem,   nres_baseflow,  OUT)
    dsl$jsb Get_var2D_onChunk( hd_mem,   ret_riverflow,  OUT)
    dsl$jsb Get_var2D_onChunk( hd_mem,   nres_riverflow, OUT)


    ! Configuration parameters from the namelist
    debug_hd          = hd_conf%debug_hd
    routing_scheme    = hd_conf%routing_scheme


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

    dsl$jsb Def_config(hd_conf, hd)
    dsl$jsb Def_memory(srf_mem, srf)
    dsl$jsb Def_memory(soil_mem, soil)

    CHARACTER(len=*),  PARAMETER  :: routine = TRIM(modname)//':hd_riverflow'

    ! Pointers to variables in memory
    REAL(wp), POINTER ::   &
      & area            (:,:), &
      & land            (:,:), &
      & coast_ocean     (:,:), &
      & ocean           (:,:), &
      & discharge       (:,:), &
      & discharge_ocean (:,:)

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

    dsl$jsb Get_config(hd_conf, HD)

    ! Get reference to variables for current block
    !
    dsl$jsb Get_memory(srf_mem, SRF)
    dsl$jsb Get_memory(soil_mem, SOIL)

    area            => grid     %area                (:,:)
    land            => srf_mem  %land%ptr            (:,:)
    coast_ocean     => srf_mem  %coast_ocean%ptr     (:,:)
    ocean           => srf_mem  %ocean%ptr           (:,:)
    discharge       => soil_mem %discharge%ptr       (:,:)
    discharge_ocean => soil_mem %discharge_ocean%ptr (:,:)

    is_in_domain    => grid%patch%cells%decomp_info%owner_mask  (:,:)

    ! Configuration parameters from the namelist
    routing_scheme = hd_conf%routing_scheme


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

    dsl$jsb Def_config(hd_conf, hd)
    dsl$jsb Def_memory(srf_mem, srf)
    dsl$jsb Def_memory(soil_mem, soil)

    CHARACTER(len=*),  PARAMETER  :: routine = TRIM(modname)//':hd_init_water_budget'

    ! Pointers to variables in memory
    REAL(wp), POINTER ::    &
      & area         (:,:), &
      & land         (:,:), &
      & glacier      (:,:), &
      & runoff       (:,:), &
      & drainage     (:,:), &
      & pme_glacier  (:,:), &
      & discharge    (:,:)

    LOGICAL, POINTER  ::   &
      & is_in_domain    (:,:) ! T: cell in domain (not halo) 

    REAL(wp), ALLOCATABLE  ::   &
      & in_domain       (:,:) ! 1: cell in domain, 0: halo cell

    LOGICAL :: debug_hd

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    model => get_model(model_id)
    tile => model%Get_top_tile()        ! HD works on the root tile
    grid => Get_grid(model%grid_id)

    dsl$jsb Get_config(hd_conf, HD)

    ! Get reference to variables in memory
    dsl$jsb Get_memory(srf_mem, SRF)
    dsl$jsb Get_memory(soil_mem, SOIL)

    area           => grid                        %area        (:,:)
    land           => srf_mem   %land%ptr        (:,:)
    glacier        => srf_mem   %glacier%ptr     (:,:)
    runoff         => soil_mem  %runoff%ptr      (:,:)
    drainage       => soil_mem  %drainage%ptr    (:,:)
    pme_glacier    => soil_mem  %pme_glacier%ptr (:,:)
    discharge      => soil_mem  %discharge%ptr   (:,:)
    is_in_domain   => grid%patch%cells%decomp_info%owner_mask  (:,:)

    ! Configuration parameters from the namelist
    debug_hd          = hd_conf%debug_hd

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

    Def_config(hd_conf, hd)
    Def_memory(srf_mem, srf)
    Def_memory(soil_mem, soil)

    CHARACTER(len=*),  PARAMETER  :: routine = TRIM(modname)//':hd_check_water_budget'

    ! Pointers to variables in memory
    REAL(wp), POINTER ::       &
      & area            (:,:), &
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

    dsl$jsb Get_config(hd_conf, HD)

    ! Get reference to variables in memory
    dsl$jsb Get_memory(srf_mem, SRF)
    dsl$jsb Get_memory(soil_mem, SOIL)

    area            => grid         %area            (:,:)
    land            => srf_mem  %land%ptr            (:,:)
    discharge_ocean => soil_mem %discharge_ocean%ptr (:,:)
    is_in_domain    => grid%patch%cells%decomp_info%owner_mask (:,:)

    ! Configuration parameters from the namelist
    debug_hd          = hd_conf%debug_hd
    diag_water_budget = hd_conf%diag_water_budget

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
