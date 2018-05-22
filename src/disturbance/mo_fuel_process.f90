!>
!! @brief Contains the routines for the fuel processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-07-29)
!! 2018-04-25, StW: Changed routine name for consistency with namelist settings
!!
MODULE mo_fuel_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_fuel_jsb3
  PUBLIC :: calc_c2g_lightning
  PUBLIC :: calc_nesterov
  PUBLIC :: calc_human_ignition

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fuel_process'

CONTAINS

  ! Calculate fuel (JSBACH3 algorithm)
  PURE ELEMENTAL SUBROUTINE calc_fuel_jsb3( &
    & c_acid_ag1,       & ! in
    & c_water_ag1,      & ! in
    & c_ethanol_ag1,    & ! in
    & c_nonsoluble_ag1, & ! in
    & c_acid_ag2,       & ! in
    & c_water_ag2,      & ! in
    & c_ethanol_ag2,    & ! in
    & c_nonsoluble_ag2, & ! in
    & fuel                   & ! out
    & )

    ! Input Arguments
    REAL(wp), INTENT(in)  :: c_acid_ag1
    REAL(wp), INTENT(in)  :: c_water_ag1
    REAL(wp), INTENT(in)  :: c_ethanol_ag1
    REAL(wp), INTENT(in)  :: c_nonsoluble_ag1
    REAL(wp), INTENT(in)  :: c_acid_ag2
    REAL(wp), INTENT(in)  :: c_water_ag2
    REAL(wp), INTENT(in)  :: c_ethanol_ag2
    REAL(wp), INTENT(in)  :: c_nonsoluble_ag2

    ! Output Arguments
    REAL(wp), INTENT(OUT)   :: fuel

    ! Locals

    ! ---------------------------
    ! Go

    fuel =  c_acid_ag1    + c_water_ag1      &
        & + c_ethanol_ag1 + c_nonsoluble_ag1 & 
        & + c_acid_ag2    + c_water_ag2      &
        & + c_ethanol_ag2 + c_nonsoluble_ag2

  END SUBROUTINE calc_fuel_jsb3

  SUBROUTINE calc_c2g_lightning(var_light,var_c2g,lc2g,frq,c2g,lat)
#ifdef HAVE_MOFFLEDAC
  USE mo_input,          ONLY: input_var_list, InputVarIsUpdated, InputVarIsPresent

  TYPE (input_var_list), POINTER :: var_light, var_c2g
  LOGICAL,         INTENT(in)    :: lc2g
  REAL(wp),        INTENT(inout) :: frq(:)
  REAL(wp),        INTENT(in)    :: c2g(:), lat(:)

    IF (InputVarIsUpdated(var_light) .AND. .NOT. lc2g) THEN
      IF (InputVarIsPresent(var_c2g)) THEN
        frq(:) = frq(:) / c2g(:)
      ELSE
        frq(:) = frq(:)*0.1_wp*(1._wp+lat(:)*lat(:)/900._wp)
      ENDIF
   ENDIF
#endif

  END SUBROUTINE calc_c2g_lightning

  ! To be called once a day
  PURE ELEMENTAL SUBROUTINE calc_nesterov(NI_acc,precip,temp_min,temp_max)
  REAL(wp), INTENT(inout) :: NI_acc
  REAL(wp), INTENT(in)    :: precip, & ! Expected in [mm (/day)]
                             temp_min, temp_max ! Expected in [degC]

  REAL(wp) :: dNI

    IF ((precip < 3._wp) .AND. (temp_min > 4._wp)) THEN
      dNI = temp_max * (temp_max - (temp_min - 4._wp))
      NI_acc = NI_acc + dNI
    ELSE
      NI_acc = 0._wp
    ENDIF

  END SUBROUTINE calc_nesterov

  ! Calculate human ignitions based on population density and kultural propensity of igniting fires
  PURE ELEMENTAL SUBROUTINE calc_human_ignition(popden,a_nd,ignitions)
  REAL(wp), INTENT(in)  :: popden
  REAL(wp), INTENT(in)  :: a_nd
  REAL(wp), INTENT(out) :: ignitions

    ! using equation from c code: parameter=6.8, now in namelist default= 5
    ignitions = EXP(-0.5_wp * (popden**0.5_wp)) * a_nd / 100._wp * popden

  END SUBROUTINE calc_human_ignition

#endif
END MODULE mo_fuel_process
