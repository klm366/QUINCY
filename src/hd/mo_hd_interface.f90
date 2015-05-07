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
  USE mo_exception,       ONLY: message
  USE mo_util,            ONLY: one_of
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
      & pme_glacier (:)

    REAL(wp) :: water_budget   ! global sum of water playing a role in this routine. Needed for water conservation test.

    !>>>
    REAL(wp), ALLOCATABLE :: flfmem(:,:), fgmem(:,:), frfmem(:,:), finfl(:,:)
    REAL(wp), ALLOCATABLE :: discharge(:), gridarea(:)
    REAL(wp)              :: conservation_test
    !<<<

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_hd'

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

    !>>>
    ALLOCATE (flfmem(nl,nb))
    ALLOCATE (fgmem(nl,nb))
    ALLOCATE (frfmem(nl,nb))
    ALLOCATE (finfl(nl,nb))
    ALLOCATE (discharge(nc))
    ALLOCATE (gridarea(nc))
    !<<<

    ! calculate global land water budget
    !
    water_budget = &
           SUM(runoff*gridarea*land)        & !< runoff
         + SUM(drainage*gridarea*land)      & !< drainage
         + SUM(flfmem)                      & !< water in the overland flow reservoirs
         + SUM(fgmem)                       & !< water in the baseflow reservoir
         + SUM(frfmem)                      & !< water in the riverflow reservoirs
         + SUM(finfl)                         !< inflow data from the hd restart file
    conservation_test=water_budget


    water_budget = &
           SUM(discharge*gridarea*land)     & !< river discharge
         + SUM(flfmem)                      & !< water in the overland flow reservoirs
         + SUM(fgmem)                       & !< water in the baseflow reservoir
         + SUM(frfmem)                      & !< water in the riverflow reservoirs
         + SUM(finfl)                         !< inflow data from the hd restart file


  END SUBROUTINE update_hd


#endif
END MODULE mo_hd_interface
