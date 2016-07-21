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

  PUBLIC :: calc_wskin_fractions_lice, calc_wskin_fractions_veg, calc_wskin_fractions_bare

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_process'

CONTAINS

  ELEMENTAL SUBROUTINE calc_wskin_fractions_lice( &
    & w_snow_lice,                                & ! in
    & fract_snow_lice                             & ! out
    & )

    REAL(wp),  INTENT(in)    :: w_snow_lice
    REAL(wp),  INTENT(out)   :: fract_snow_lice ! Snow fraction on (lake) ice

    fract_snow_lice = TANH(w_snow_lice * 100._wp)

  END SUBROUTINE calc_wskin_fractions_lice
   
   
  ELEMENTAL SUBROUTINE calc_wskin_fractions_veg( &
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
    & fract_snow_can, & ! inout
    & fract_water_can,    & ! out
    & fract_water_soil,   & ! out
    & fract_water_veg,    & ! out
    & fract_snow_soil,    & ! out
    & fract_snow_veg      & ! out
    & )

    USE mo_soil_process,           ONLY: Get_snow_fract_noforest
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
                                  w_snow_can

      REAL(wp),  INTENT(inout) ::                  &
                                 fract_snow_can
      REAL(wp),  INTENT(out)   ::                  &
                                 fract_water_can,  &
                                 fract_water_soil, &
                                 fract_water_veg,  &
                                 fract_snow_soil,  &
                                 fract_snow_veg

     REAL(wp) ::               &
                 t_srf_old,    &
                 qsat_srf_old, &
                 evap_pot,     &
                 wsr_max
    
    ! Maximum water holding capacity of skin reservoir (soil + canopy)
    wsr_max = wsr_capacity * (1._wp + lai * fract_fpc_max) 

    ! Wet skin reservoir area fraction on soil (between 0 and 1)
    fract_water_soil = MERGE(MIN(1._wp, (w_skin_soil + w_skin_can) / wsr_max), 0._wp, wsr_max > EPSILON(1._wp))

    ! Wet skin fraction on canopy (between 0 and 1)
!!$    fract_water_can = MERGE(MIN(1._wp, w_skin_can / wsr_max), 0._wp, wsr_max > EPSILON(1._wp))
    ! TBD
    fract_water_can = 0._wp

    
    ! Snow fraction on soil
    fract_snow_soil = Get_snow_fract_noforest(w_snow_soil, oro_stddev) ! snow on soil below forest
    fract_snow_soil = MERGE(fract_snow_can, fract_snow_soil, &
                               fract_snow_soil < EPSILON(1._wp) .AND. fract_snow_can > EPSILON(1._wp))

    ! Snow fraction on canopy
    fract_snow_can = MERGE(MIN(1._wp, w_snow_can / wsr_max), 0._wp, wsr_max > EPSILON(1._wp))

    
    ! Potential evaporation using old values of air and surface humidity
    t_srf_old = t_srf
    qsat_srf_old = qsat(t_srf_old, press_srf)
    evap_pot =  heat_transfer_coef(drag_srf, steplen) * (qsat_srf_old - q_air) ! Positive upwards
    
    ! Reduce snow cover on soil if snow loss during the time step due to potential evaporation is larger than
    ! snow water content from soil and canopy together

    IF (fract_snow_soil > 0._wp .AND. w_snow_soil+w_snow_can > EPSILON(1._wp)) THEN
      fract_snow_soil = fract_snow_soil / MAX(1._wp, fract_snow_soil * evap_pot * dtime / (rhoh2o * &
        & MAX(EPSILON(1._wp),(w_snow_soil + w_snow_can))))
    END IF
    ! Same for water cover on soil
    IF (fract_water_soil > 0._wp .AND. w_skin_soil+w_skin_can > EPSILON(1._wp)) THEN
      fract_water_soil = fract_water_soil / MAX(1._wp, (1._wp - fract_snow_soil) * evap_pot * dtime / (rhoh2o * &
        & MAX(EPSILON(1._wp), (w_skin_soil + w_skin_can))))
    END IF

    ! Merge water on and below canopy
    fract_water_veg = fract_water_soil + fract_water_can ! <---- dadurch lookup table overflow
    fract_snow_veg  = fract_snow_soil ! + fract_snow_can    TBD

   END SUBROUTINE calc_wskin_fractions_veg



  ELEMENTAL SUBROUTINE calc_wskin_fractions_bare( &
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
    & w_snow_soil,    & ! in
    & fract_water_soil,  & ! out
    & fract_snow_soil    & ! out
    & )

    USE mo_soil_process,           ONLY: Get_snow_fract_noforest
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
                                  w_snow_soil

      REAL(wp),  INTENT(out)   ::                  &
                                 fract_water_soil, &
                                 fract_snow_soil

     REAL(wp) ::               &
                 t_srf_old,    &
                 qsat_srf_old, &
                 evap_pot,     &
                 wsr_max


    ! Maximum water holding capacity of skin reservoir
    ! R: TBD: Note, wsr_max in JSBACH3 is thought to be the holding capacity for the whole (!) skin reservoir: soil + canopy
    !         In JSBACH3 there is no bare soil tile but only bare soil as a part of a plant tile. In JSBACH3 there is no
    !         influence on wsr_max when the canopy of the tile is nearly zero. Thus, the following line of JSBACH4 is consistent
    !         with JSBACH3. BUT: in JSBACH3 it was very unlikely that the canopy of the tile gets to zero and in the following line
    !         this unrealistic case is taken as standard for all bare soil surface types.
    !         fract_fpc_max = Maximum fraction of tile that can be covered by leaves (assuming infinite lai)   
    wsr_max = wsr_capacity * (1._wp + lai * fract_fpc_max)

    ! Wet skin reservoir area fraction on soil (between 0 and 1)
    ! R: Deleted the contribution of w_skin_can for bare soil. Consistent to JSBACH3 (see comment above) and in this line it seems
    !    reasonable.
    fract_water_soil = MERGE(MIN(1._wp, w_skin_soil / wsr_max), 0._wp, wsr_max > EPSILON(1._wp))


    ! Snow fraction on soil
    fract_snow_soil = Get_snow_fract_noforest(w_snow_soil, oro_stddev) ! snow on soil below forest
    ! R: The following lines of JSBACH3 was to guarantee that when the soil is without snow, at least the snow from the canopy
    !    that falls on the soil is represented...
    !fract_snow_soil = MERGE(fract_snow_can, fract_snow_soil, 
    !                           fract_snow_soil < EPSILON(1._wp) .AND. fract_snow_can > EPSILON(1._wp))


    ! Potential evaporation using old values of air and surface humidity
    t_srf_old = t_srf
    qsat_srf_old = qsat(t_srf_old, press_srf)
    evap_pot =  heat_transfer_coef(drag_srf, steplen) * (qsat_srf_old - q_air) ! Positive upwards

    ! Reduce snow cover on soil if snow loss during the time step due to potential evaporation is larger than
    ! snow water content from soil and canopy together
    ! R: In the following lines I just kicked out the contribution of w_skin_can as compared to JSBACH3. For the intension of
    !    the parametrization this can be done so easily for the bare type and it is consistent with JSBACH3.
    IF (fract_snow_soil > 0._wp .AND. w_snow_soil > EPSILON(1._wp)) THEN
      fract_snow_soil = fract_snow_soil / MAX(1._wp, fract_snow_soil * evap_pot * dtime / (rhoh2o * &
        & MAX(EPSILON(1._wp),(w_snow_soil))))
    END IF
    ! Same for water cover on soil
    IF (fract_water_soil > 0._wp .AND. w_skin_soil > EPSILON(1._wp)) THEN
      fract_water_soil = fract_water_soil / MAX(1._wp, (1._wp - fract_snow_soil) * evap_pot * dtime / (rhoh2o * &
        & MAX(EPSILON(1._wp), (w_skin_soil))))
    END IF

   END SUBROUTINE calc_wskin_fractions_bare
#endif
END MODULE mo_hydro_process
