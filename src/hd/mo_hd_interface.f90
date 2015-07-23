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
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_tile_class,  ONLY: t_jsb_tile
  USE mo_jsb_base,        ONLY: get_model
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish, message, message_text
  USE mo_sync,            ONLY: SYNC_C, sync_patch_array, global_sum_array
  USE mo_impl_constants,  ONLY: vname_len
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
    TYPE(t_jsb_model), POINTER    :: model

    ! Pointers to variables in memory
    REAL(wp), POINTER :: &
      & area        (:), &
      & land        (:), &
      & lake        (:), &
      & glacier     (:), &
      & runoff      (:), &
      & drainage    (:), &
      & pme_glacier (:), &
      & discharge   (:)

    REAL(wp), POINTER ::    &
      & ret_overlflow  (:), &
      & nres_overlflow (:), &
      & ret_baseflow   (:), &
      & nres_baseflow  (:), &
      & ret_riverflow  (:), &
      & nres_riverflow (:)

    REAL(wp), POINTER ::     &
      & overlflow_res (:,:), &
      & baseflow_res  (:,:), &
      & riverflow_res (:,:)

    REAL(wp), ALLOCATABLE :: inflow(:)
    REAL(wp), ALLOCATABLE :: outflow_runoff(:)
    REAL(wp), ALLOCATABLE :: outflow_drainage(:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_hd'

    LOGICAL                  :: debug_hd
    CHARACTER(len=vname_len) :: routing_scheme

    IF (.NOT. ASSOCIATED(tile%mem%hd)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    area           => model%grid                    %area                (ics:ice, iblk)
    land           => tile%mem%srf  %phys_state     %land                (ics:ice, iblk)
    lake           => tile%mem%srf  %phys_state     %lake                (ics:ice, iblk)
    glacier        => tile%mem%srf  %phys_state     %glacier             (ics:ice, iblk)
    runoff         => tile%mem%soil %water_flux     %runoff              (ics:ice, iblk)
    drainage       => tile%mem%soil %water_flux     %drainage            (ics:ice, iblk)
    pme_glacier    => tile%mem%soil %water_flux     %pme_glacier         (ics:ice, iblk)
    discharge      => tile%mem%soil %water_flux     %discharge           (ics:ice, iblk)

    overlflow_res  => tile%mem%hd   %hd_state       %overlflow_res       (ics:ice, :, iblk)
    baseflow_res   => tile%mem%hd   %hd_state       %baseflow_res        (ics:ice, :, iblk)
    riverflow_res  => tile%mem%hd   %hd_state       %riverflow_res       (ics:ice, :, iblk)

    ret_overlflow  => tile%mem%hd   %hd_param       %ret_overlflow       (ics:ice, iblk)
    nres_overlflow => tile%mem%hd   %hd_param       %nres_overlflow      (ics:ice, iblk)
    ret_baseflow   => tile%mem%hd   %hd_param       %ret_baseflow        (ics:ice, iblk)
    nres_baseflow  => tile%mem%hd   %hd_param       %nres_baseflow       (ics:ice, iblk)
    ret_riverflow  => tile%mem%hd   %hd_param       %ret_riverflow       (ics:ice, iblk)
    nres_riverflow => tile%mem%hd   %hd_param       %nres_riverflow      (ics:ice, iblk)
 

    ! Configuration parameters from the namelist
    debug_hd          = model%config%hd%debug_hd
    routing_scheme    = model%config%hd%routing_scheme


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

    CASE ('equally_to_coast', 'equally_to_ocean')
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
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_tile),  POINTER    :: tile
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
      & area_ocean            ! total area of all ocean grid cells

    LOGICAL, POINTER  ::   &
      & is_in_domain    (:,:) ! T: cell in domain (not halo)

    CHARACTER(len=vname_len) :: &
      & routing_scheme        ! river routing scheme

    REAL(wp), ALLOCATABLE  ::   &
      & in_domain       (:,:) ! 1: cell in domain, 0: halo cell 

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    model => get_model(model_id)
    tile => model%land        ! HD works on the root tile

    ! Get reference to variables for current block
    !
    area            => model%grid                    %area            (:,:)
    land            => tile%mem%srf  %phys_state     %land            (:,:)
    coast_ocean     => tile%mem%srf  %phys_state     %coast_ocean     (:,:)
    ocean           => tile%mem%srf  %phys_state     %ocean           (:,:)
    discharge       => tile%mem%soil %water_flux     %discharge       (:,:)
    discharge_ocean => tile%mem%soil %water_flux     %discharge_ocean (:,:)

    is_in_domain    => model%grid%patch%cells%decomp_info%owner_mask  (:,:)

    ! Configuration parameters from the namelist
    routing_scheme = model%config%hd%routing_scheme


    ! Dummy discharge: the discharge is equally distributed to either
    !     - all coastal ocean cells (routing_scheme = 'equally_to_coast')  or
    !     - all ocean cells (routing_scheme = 'equally_to_ocean')
    ! (routing_scheme = 'full' for real river routing not yet implemented)

    ! Mask out all halo cells for the global sums (otherwise these cells are counted twice)
    ALLOCATE (in_domain(SIZE(is_in_domain(:,1)),SIZE(is_in_domain(1,:))))
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
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_tile),  POINTER    :: tile
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
    tile => model%land        ! HD works on the root tile

    ! Get reference to variables in memory
    area           => model%grid                       %area        (:,:)
    land           => tile%mem%srf   %phys_state       %land        (:,:)
    glacier        => tile%mem%srf   %phys_state       %glacier     (:,:)
    runoff         => tile%mem%soil  %water_flux       %runoff      (:,:)
    drainage       => tile%mem%soil  %water_flux       %drainage    (:,:)
    pme_glacier    => tile%mem%soil  %water_flux       %pme_glacier (:,:)
    discharge      => tile%mem%soil  %water_flux       %discharge   (:,:)
    is_in_domain   => model%grid%patch%cells%decomp_info%owner_mask (:,:)

    ! Configuration parameters from the namelist
    debug_hd          = model%config%hd%debug_hd

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
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_tile),  POINTER    :: tile
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
    tile => model%land        ! HD works on the root tile

    ! Get reference to variables in memory
    area            => model%grid                    %area            (:,:)
    land            => tile%mem%srf  %phys_state     %land            (:,:)
    discharge_ocean => tile%mem%soil %water_flux     %discharge_ocean (:,:)
    is_in_domain    => model%grid%patch%cells%decomp_info%owner_mask  (:,:)

    ! Configuration parameters from the namelist
    debug_hd          = model%config%hd%debug_hd
    diag_water_budget = model%config%hd%diag_water_budget

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
