!>
!! @brief Contains the routines for the vegetation processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2019-09-11)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_carbon_veg
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp

  USE mo_jsb_lctlib_class, ONLY: t_lctlib_element

  dsl4jsb_Use_pools(CARBON_) t_plant_pool

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_allocation

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_veg'

CONTAINS

  ! -------------------------------------------------------------------------------------------------------------------------------
  !>
  !! Potential plant carbon allocation (exception: actual allocation for NPP < 0)
  !!
  SUBROUTINE calc_allocation( &
    & update_interval,      & ! in
    & with_nitrogen,        & ! in
    & lctlib,               & ! in
    & plant,                & ! in
    & LAI,                  & ! in
    & maxLAI,               & ! in
    & NPP_pot,              & ! in
    & alloc,                & ! inout
    & root_exudates,        & ! out
    & NPP_act,              & ! out
    & NPP_flux_correction   & ! out
    & )

    REAL(wp),               INTENT(in)    :: update_interval  !< Frequency of updates [s]
    LOGICAL,                INTENT(in)    :: with_nitrogen
    TYPE(t_lctlib_element), INTENT(in)    :: lctlib
    TYPE(t_plant_pool),     INTENT(in)    :: plant            !< Plant state pools [mol/m^2(canopy)]
    REAL(wp),               INTENT(in)    :: LAI(:)           !< Leaf area index []
    REAL(wp),               INTENT(in)    :: maxLAI(:)        !< Maximum LAI []
    REAL(wp),               INTENT(in)    :: NPP_pot(:)       !< Potential NPP rate (from assimilation) [mol(C)/m^2/s]
    TYPE(t_plant_pool),     INTENT(inout) :: alloc            !< Flux pool for allocation [mol/m^2/s]
    REAL(wp),               INTENT(out)   :: &
      & root_exudates(:), &                                   !< Part of NPP_pot that goes into root exudates [mol(C)/m^2/s]
      & NPP_act(:), &                                         !< Actual NPP rate [mol(C)/m^2/s] after N-limitation
                                                              !  and excess carbon drop
      & NPP_flux_correction(:)                                !< NPP flux correction [mol(C)/m^2/s]

    REAL(wp), DIMENSION(SIZE(LAI)) :: &
      & c_reserve_optimal, &
      & exess_carbon

    INTEGER :: ic

    root_exudates(:) = 0._wp
    NPP_act(:) = 0._wp
    NPP_flux_correction(:) = 0._wp

    c_reserve_optimal(:) = lctlib%reserveC2leafC * maxLAI(:) / lctlib%specificLeafArea_C

    DO ic=1,SIZE(LAI)

      IF (NPP_pot(ic) <= 0._wp) THEN
        ! De-allocation of reserve carbon.
        ! No allocation to green and wood pools, or root exudates.
        ! Therefore, NPP_act = alloc%reserve%carbon

        $alloc%green%carbon(ic)$   = 0._wp
        $alloc%wood%carbon(ic)$    = 0._wp
        root_exudates(ic)          = 0._wp
        $alloc%reserve%carbon(ic)$ = MAX(NPP_pot(ic), - $plant%reserve%carbon(ic)$ / update_interval)  ! Eq. 6.16

      END IF

    END DO

  END SUBROUTINE calc_allocation

#endif
END MODULE mo_carbon_veg
