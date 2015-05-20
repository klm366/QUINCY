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

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_hd
  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_interface'

CONTAINS

  SUBROUTINE update_hd(tile, iblk, ics, ice, nc)

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model

    ! Pointers to variables in memory
    REAL(wp), POINTER :: &
      & land        (:), &
      & lake        (:), &
      & glacier     (:), &
      & runoff      (:), &
      & drainage    (:), &
      & evapotrans  (:), &
      & pme_glacier (:), &
      & area        (:)

    REAL(wp) :: water_budget   ! global sum of water playing a role in this routine. Needed for water conservation test.

    !>>>
    REAL(wp), ALLOCATABLE :: flfmem(:,:), fgmem(:,:), frfmem(:,:), finfl(:,:)
    REAL(wp), ALLOCATABLE :: discharge(:)
    REAL(wp)              :: conservation_test
    !<<<

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_hd'

    LOGICAL :: debug_hd, diag_water_budget

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of hydrological discharge model')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    land           => tile%mem%srf  %phys_state     %land                (ics:ice, iblk)
    lake           => tile%mem%srf  %phys_state     %lake                (ics:ice, iblk)
    glacier        => tile%mem%srf  %phys_state     %glacier             (ics:ice, iblk)
    runoff         => tile%mem%soil %water_flux     %runoff              (ics:ice, iblk)
    drainage       => tile%mem%soil %water_flux     %drainage            (ics:ice, iblk)
    evapotrans     => tile%mem%srf  %water_flux     %evapotrans          (ics:ice, iblk)

    area           => model%grid%area                                    (ics:ice, iblk)

    ! Configuration parameters from the namelist
    debug_hd          = model%config%hd%debug_hd
    diag_water_budget = model%config%hd%diag_water_budget

    !>>>
    ALLOCATE (flfmem(nl,nb))
    flfmem = 0._wp
    ALLOCATE (fgmem(nl,nb))
    fgmem = 0._wp
    ALLOCATE (frfmem(nl,nb))
    frfmem = 0._wp
    ALLOCATE (finfl(nl,nb))
    finfl = 0._wp
    ALLOCATE (discharge(nc))
    discharge = 0._wp
    !<<<

print*, 'vg: runoff  : ',  MINVAL(runoff),   MAXVAL(runoff),   SUM(runoff)
print*, 'vg: drainage: ',  MINVAL(drainage), MAXVAL(drainage), SUM(drainage)

    ! calculate global land water budget
    !
    water_budget = &
           SUM(runoff*area*land)            & !< runoff
         + SUM(drainage*area*land)          & !< drainage
         + SUM(flfmem)                      & !< water in the overland flow reservoirs
         + SUM(fgmem)                       & !< water in the baseflow reservoir
         + SUM(frfmem)                      & !< water in the riverflow reservoirs
         + SUM(finfl)                         !< inflow data from the hd restart file
    conservation_test=water_budget

    ! dummy version of HD model
    !
    discharge(:) = land(:) * (runoff(:) + drainage(:))


    ! re-calculate global land water budget
    !
    water_budget = &
           SUM(discharge*area*land)         & !< river discharge
         + SUM(flfmem)                      & !< water in the overland flow reservoirs
         + SUM(fgmem)                       & !< water in the baseflow reservoir
         + SUM(frfmem)                      & !< water in the riverflow reservoirs
         + SUM(finfl)                         !< inflow data from the hd restart file

    ! water conservation test
    !
    conservation_test = conservation_test - water_budget
    IF (diag_water_budget) THEN
      WRITE (message_text,*) 'Water budget change: ', conservation_test,' m3/s'
    END IF
    CALL message (TRIM(routine),message_text)
    IF (water_budget /= 0._wp) THEN
    IF (ABS(conservation_test/water_budget) > 1.E-7_wp) THEN
      WRITE (message_text,*) 'Water conservation problem: budget change: ', &
           conservation_test,' m3/s'
      CALL finish (TRIM(routine),message_text)
    END IF
    END IF
    IF (diag_water_budget) THEN
      WRITE (message_text,*) '  global discharge: ', SUM(discharge*area*land), ' m3/s'
      CALL message (TRIM(routine),message_text)
    END IF

    !>>>
    DEALLOCATE (flfmem, fgmem, frfmem, finfl, discharge)
    !<<<

  END SUBROUTINE update_hd


#endif
END MODULE mo_hd_interface
