!>
!! @brief Contains the routines for the hydro processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-29)
!!
MODULE mo_hydro_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_snow_and_skin_fraction

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_process'

CONTAINS

  ELEMENTAL SUBROUTINE calc_snow_and_skin_fraction( &
    & dtime,          & ! in
    & steplen,        & ! in
    & wsr_capacity,   & ! in
    & lai,            & ! in
    & oro_stddev,     & ! in
    & fract_fpc_max,  & ! in
    & t_srf,          & ! in
    & press_srf,      & ! in
    & drag_srf,       & ! in
    & q_air,          & ! in
    & w_skin_soil,    & ! in
    & w_skin_can,     & ! in
    & w_snow_soil,    & ! in
    & w_snow_can,     & ! in
    & w_snow_ice,     & ! in
    & fract_snow_can, & ! inout
    & fract_water_can,    & ! out
    & fract_water_soil,   & ! out
    & fract_water_srf,    & ! out
    & fract_snow_soil,    & ! out
    & fract_snow_srf,     & ! out
    & fract_snow_ice      & ! out
    & )

    USE mo_soil_process,           ONLY: snow_fract_noforest
    USE mo_phy_schemes,            ONLY: heat_transfer_coef, qsat
    USE mo_jsb_physical_constants, ONLY: rhoh2o

      REAL(wp),  INTENT(in)    ::                &
                                  dtime,         &
                                  steplen,       &
                                  lai,           &
                                  oro_stddev,    &
                                  fract_fpc_max, &
                                  t_srf,         &
                                  press_srf,     &
                                  drag_srf,      &
                                  q_air,         &
                                  wsr_capacity,  &
                                  w_skin_soil,   &
                                  w_skin_can,    &
                                  w_snow_soil,   &
                                  w_snow_can,    &
                                  w_snow_ice

      REAL(wp),  INTENT(inout) ::                    &
                                 fract_snow_can
      REAL(wp),  INTENT(out)   ::                    &
                                 fract_water_can,    &
                                 fract_water_soil,   &
                                 fract_water_srf,    &
                                 fract_snow_soil,    &
                                 fract_snow_srf,     &
                                 fract_snow_ice

     REAL(wp) ::               &
                 t_srf_old,    &
                 qsat_srf_old, &
                 evap_pot,     &
                 wsr_max

    ! Maximum water holding capacity of skin reservoir (soil + canopy)
    wsr_max = wsr_capacity * (1._wp + lai * fract_fpc_max)

    ! Wet skin reservoir area fraction on soil
    fract_water_soil = MERGE(MIN(1._wp, (w_skin_soil + w_skin_can) / wsr_max), 0._wp, wsr_max > EPSILON(1._wp))

    ! Wet skin fraction on canopy
!!$    fract_water_can = MERGE(MIN(1._wp, w_skin_can / wsr_max), 0._wp, wsr_max > EPSILON(1._wp))
    ! TBD
    fract_water_can = 0._wp

    ! Snow fraction on soil
    fract_snow_soil = snow_fract_noforest(w_snow_soil, oro_stddev)
    fract_snow_soil = MERGE(fract_snow_can, fract_snow_soil, &
                               fract_snow_soil < EPSILON(1._wp) .AND. fract_snow_can > EPSILON(1._wp))

    ! Snow fraction on canopy
    fract_snow_can = MERGE(MIN(1._wp, w_snow_can / wsr_max), 0._wp, wsr_max > EPSILON(1._wp))

    ! Snow fraction on (lake) ice
    fract_snow_ice = TANH(w_snow_ice * 100._wp)

    ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than
    ! equivalent snow water content from soil and canopy; same for skin reservoir
    !
    ! Potential evaporation using old values of air and surface humidity
    !
    t_srf_old = t_srf
    qsat_srf_old = qsat(t_srf_old, press_srf)
    evap_pot =  heat_transfer_coef(drag_srf, steplen) * (qsat_srf_old - q_air) ! Positive upwards

    IF (fract_snow_soil > 0._wp .AND. w_snow_soil+w_snow_can > EPSILON(1._wp)) THEN
      fract_snow_soil = fract_snow_soil / MAX(1._wp, fract_snow_soil * evap_pot * dtime / (rhoh2o * &
        & MAX(EPSILON(1._wp),(w_snow_soil + w_snow_can))))
    END IF
    IF (fract_water_soil > 0._wp .AND. w_skin_soil+w_skin_can > EPSILON(1._wp)) THEN
      fract_water_soil = fract_water_soil / MAX(1._wp, (1._wp - fract_snow_soil) * evap_pot * dtime / (rhoh2o * &
        & MAX(EPSILON(1._wp), (w_skin_soil + w_skin_can))))
    END IF

    fract_water_srf = fract_water_soil + fract_water_can ! <---- dadurch lookup table overflow
    fract_snow_srf = fract_snow_soil ! + fract_snow_can    TBD


   END SUBROUTINE calc_snow_and_skin_fraction

#endif
END MODULE mo_hydro_process
