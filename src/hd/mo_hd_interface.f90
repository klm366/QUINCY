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
  USE mo_hd_constants,    ONLY: nl, nb
  USE mo_sync,            ONLY: SYNC_C, sync_patch_array, global_sum_array

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_hd, hd_riverflow, hd_init_water_budget, hd_check_water_budget
  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_interface'

CONTAINS

  !------------------------------------------------------------------------------------------------
  SUBROUTINE update_hd(tile, iblk, ics, ice, nc)
  !------------------------------------------------------------------------------------------------
  ! the routine is called within the block-loop from jsbach_interface. There is no information on 
  ! neighbouring grid cells (on other paches) availabe.
  ! We calculate here the local discharge of each grid cell.
  !------------------------------------------------------------------------------------------------

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model

    ! Pointers to variables in memory
    REAL(wp), POINTER :: &
      & area        (:), &
      & land        (:), &
      & lake        (:), &
      & coast_ocean (:), &
      & glacier     (:), &
      & runoff      (:), &
      & drainage    (:), &
      & pme_glacier (:), &
      & discharge   (:)

    REAL(wp), ALLOCATABLE :: flfmem(:,:), fgmem(:,:), frfmem(:,:), finfl(:,:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_hd'

    LOGICAL :: debug_hd

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    area           => model%grid                    %area                (ics:ice, iblk)
    land           => tile%mem%srf  %phys_state     %land                (ics:ice, iblk)
    lake           => tile%mem%srf  %phys_state     %lake                (ics:ice, iblk)
    glacier        => tile%mem%srf  %phys_state     %glacier             (ics:ice, iblk)
    coast_ocean    => tile%mem%srf  %phys_state     %coast_ocean         (ics:ice, iblk)
    runoff         => tile%mem%soil %water_flux     %runoff              (ics:ice, iblk)
    drainage       => tile%mem%soil %water_flux     %drainage            (ics:ice, iblk)
    pme_glacier    => tile%mem%soil %water_flux     %pme_glacier         (ics:ice, iblk)
    discharge      => tile%mem%soil %water_flux     %discharge           (ics:ice, iblk)

    ! Configuration parameters from the namelist
    debug_hd          = model%config%hd%debug_hd

    !>>>
    ALLOCATE (flfmem(nl,nb))
    flfmem = 0._wp
    ALLOCATE (fgmem(nl,nb))
    fgmem = 0._wp
    ALLOCATE (frfmem(nl,nb))
    frfmem = 0._wp
    ALLOCATE (finfl(nl,nb))
    finfl = 0._wp
    !<<<

    ! dummy version of HD model
    !
    discharge(:) = land(:) * (runoff(:) + drainage(:) + pme_glacier(:))

    !>>>
    DEALLOCATE (flfmem, fgmem, frfmem, finfl)
    !<<<

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
      & discharge       (:,:), &
      & discharge_ocean (:,:)

    REAL(wp) ::            &
      & global_discharge,  &  ! global sum of discharge to the ocean
      & area_coast_ocean      ! total area of coastal ocean grid cells

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    model => get_model(model_id)
    tile => model%land        ! HD works on the root tile

    ! Get reference to variables for current block
    !
    area            => model%grid                    %area            (:,:)
    land            => tile%mem%srf  %phys_state     %land            (:,:)
    coast_ocean     => tile%mem%srf  %phys_state     %coast_ocean     (:,:)
    discharge       => tile%mem%soil %water_flux     %discharge       (:,:)
    discharge_ocean => tile%mem%soil %water_flux     %discharge_ocean (:,:)

!vg     ! Synchronize halos
!vg     CALL sync_patch_array(SYNC_C, patch, discharge(:,:))

    ! Dummy discharge: the discharge is equally distributed to the coastal ocean cells
    !
    global_discharge = global_sum_array(discharge(:,:)*land(:,:)*area(:,:))
    area_coast_ocean = global_sum_array(coast_ocean(:,:)*area(:,:))
    WHERE (coast_ocean(:,:) == 1._wp)
      discharge_ocean(:,:) = area(:,:) * global_discharge/area_coast_ocean
    ELSEWHERE
      discharge_ocean(:,:) = 0._wp
    END WHERE

  END SUBROUTINE hd_riverflow

  !------------------------------------------------------------------------------------------------
  SUBROUTINE hd_init_water_budget(model_id)
  !------------------------------------------------------------------------------------------------
  ! calculation of the global water budget at the beginning of the timestep 
  !    For the moment we only check the water treated by the HD model
  !------------------------------------------------------------------------------------------------

    INTEGER, INTENT(in) :: model_id

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_tile),  POINTER    :: tile
    CHARACTER(len=*),  PARAMETER  :: routine = TRIM(modname)//':hd_init_water_budget'

    ! Pointers to variables in memory
    REAL(wp), POINTER :: &
      & area        (:,:), &
      & land        (:,:), &
      & coast_ocean (:,:), &
      & glacier     (:,:), &
      & runoff      (:,:), &
      & drainage    (:,:), &
      & pme_glacier (:,:), &
      & discharge   (:,:)

    REAL(wp) ::            &
      & water_budget          ! global sum of water playing a role in the HD mode

    LOGICAL :: debug_hd, diag_water_budget

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')

    model => get_model(model_id)
    tile => model%land        ! HD works on the root tile

    ! Get reference to variables in memory
    area           => model%grid                    %area        (:,:)
    land           => tile%mem%srf  %phys_state     %land        (:,:)
    coast_ocean    => tile%mem%srf  %phys_state     %coast_ocean (:,:)
    glacier        => tile%mem%srf  %phys_state     %glacier     (:,:)
    runoff         => tile%mem%soil %water_flux     %runoff      (:,:)
    drainage       => tile%mem%soil %water_flux     %drainage    (:,:)
    pme_glacier    => tile%mem%soil %water_flux     %pme_glacier (:,:)
    discharge      => tile%mem%soil %water_flux     %discharge   (:,:)

    ! Configuration parameters from the namelist
    debug_hd          = model%config%hd%debug_hd
    diag_water_budget = model%config%hd%diag_water_budget

    !
    ! calculate global land water budget
    !
    water_budget = &
          global_sum_array(runoff*area*land)                    & !< runoff
        + global_sum_array(drainage*area*land)                  & !< drainage
        + global_sum_array(pme_glacier*glacier*area*land)   !vg & !< P-E over glaciers
!vg         + global_sum_array(flfmem)                      & !< water in the overland flow reservoirs
!vg         + global_sum_array(fgmem)                       & !< water in the baseflow reservoir
!vg         + global_sum_array(frfmem)                      & !< water in the riverflow reservoirs
!vg         + global_sum_array(finfl)                         !< inflow data from the hd restart file

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
    REAL(wp), POINTER :: &
      & area        (:,:), &
      & land        (:,:), &
      & coast_ocean (:,:), &
      & glacier     (:,:), &
      & runoff      (:,:), &
      & drainage    (:,:), &
      & pme_glacier (:,:), &
      & discharge_ocean (:,:)

    REAL(wp), SAVE  ::  water_budget       ! global sum of water playing a role in the HD mode
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
    coast_ocean     => tile%mem%srf  %phys_state     %coast_ocean     (:,:)
    glacier         => tile%mem%srf  %phys_state     %glacier         (:,:)
    runoff          => tile%mem%soil %water_flux     %runoff          (:,:)
    drainage        => tile%mem%soil %water_flux     %drainage        (:,:)
    pme_glacier     => tile%mem%soil %water_flux     %pme_glacier     (:,:)
    discharge_ocean => tile%mem%soil %water_flux     %discharge_ocean (:,:)

    ! Configuration parameters from the namelist
    debug_hd          = model%config%hd%debug_hd
    diag_water_budget = model%config%hd%diag_water_budget

    water_budget_old = water_budget

    ! re-calculate global land water budget
    !
    water_budget = &
           global_sum_array(discharge_ocean*area) !vg  & !< river discharge
!vg          + global_sum_array(flfmem)                    & !< water in the overland flow reservoirs
!vg          + global_sum_array(fgmem)                     & !< water in the baseflow reservoir
!vg          + global_sum_array(frfmem)                    & !< water in the riverflow reservoirs
!vg          + global_sum_array(finfl)                       !< inflow data

    ! water conservation test
    !
    water_budget_change = water_budget_old - water_budget
    IF (diag_water_budget) THEN
      WRITE (message_text,*) 'Water budget change: ', water_budget_change,' m3/s'
    END IF
    CALL message (TRIM(routine),message_text)
    IF (water_budget /= 0._wp) THEN
      IF (ABS(water_budget_change/water_budget) > 1.E-7_wp) THEN
        WRITE (message_text,*) 'Water conservation problem: budget change: ', &
             water_budget_change,' m3/s'
!vg       CALL finish (TRIM(routine),message_text)
      END IF
    END IF
    IF (diag_water_budget) THEN
      WRITE (message_text,*) '  global discharge to the ocean: ', SUM(discharge_ocean*area), ' m3/s'
      CALL message (TRIM(routine),message_text)
    END IF   

  END SUBROUTINE hd_check_water_budget


#endif
END MODULE mo_hd_interface
