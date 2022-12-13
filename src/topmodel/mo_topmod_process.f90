!>
!! @brief Contains the routines for the soil and snow energy processes
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2022-12-06
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_top_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: relative_humidity_soil, init_soil_temperature, &
    & calc_vol_heat_capacity, calc_thermal_conductivity, calc_soil_temperature, calc_snow_temperature, &
    & calc_soiltemp_old, Get_liquid_max

  CHARACTER(len=*), PARAMETER :: modname = 'mo_top_process'

CONTAINS

  !!
  ! @brief      Calculates the volumetric heat capacity of soil (de Vries, 1963)
  !
  ! @param [in]  vol_heat_cap_soil_dry_mineral  Heat capacity of dry mineral soil solid
  
  ! @return     The volumetric heat capacity of soil.
  !!
!  ELEMENTAL PURE FUNCTION calc_vol_heat_capacity( &
!    & vol_heat_cap_soil_dry_mineral,              &
!    & ) RESULT(vol_heat_cap)

!  END FUNCTION calc_thermal_conductivity

!  SUBROUTINE init_soil_temperature(cmid, tclim, tsoil)

!  END SUBROUTINE init_soil_temperature

  !
  ! @brief      Calculates the soil temperature.
  !
  ! @param      nidx          The nidx
  ! @return     { description_of_the_return_value } !
  !
!  SUBROUTINE calc_soil_temperature(       &
!    & )
    !$claw expand parallel update savepoint=ser_jsb_sse_calc_soil_temperature_expand_4 counter=ser_counter
    ! z1(:,nsoil-1) = zdz2(:,nsoil) + zdz1(:,nsoil-1)
!    t_soil_acoef(:,nsoil-1) = zdz2(:,nsoil) * t_soil_sl(:,nsoil) / z1(:,nsoil-1)
!    t_soil_bcoef(:,nsoil-1) = zdz1(:,nsoil-1) / z1(:,nsoil-1)
    !$claw end expand

    ! soil layers above
    !$claw acc data present(z1,zdz1, zdz2,t_soil_acoef,t_soil_bcoef,t_soil_sl)
    !$claw acc parallel
    !$claw acc loop gang vector
!    DO jl=1,nidx
      !$claw acc loop seq
!      DO jk = nsoil-1, 2, -1
!        z1(jl,jk-1) = 1._wp / (zdz2(jl,jk) + zdz1(jl,jk-1) + zdz1(jl,jk) * (1._wp - t_soil_bcoef(jl,jk)))
!        t_soil_acoef(jl,jk-1) = (t_soil_sl(jl,jk) * zdz2(jl,jk) + zdz1(jl,jk) * t_soil_acoef(jl,jk)) * z1(jl,jk-1)
!        t_soil_bcoef(jl,jk-1) = zdz1(jl,jk-1) * z1(jl,jk-1)
!      END DO
!    END DO  
    !$claw acc end parallel
    !$claw acc end data

!#if defined(_CLAW) && defined(SERIALIZE)
    !$claw acc data create(z1_sum,t_soil_acoef_sum,t_soil_bcoef_sum)
    !$claw expand parallel update savepoint=ser_jsb_sse_calc_soil_temperature_expand_5 counter=ser_counter
!    z1_sum(:)           = SUM(z1(:,1:nsoil-1), DIM=2)
!    t_soil_acoef_sum(:) = SUM(t_soil_acoef(:,1:nsoil-1), DIM=2)
!    t_soil_bcoef_sum(:) = SUM(t_soil_bcoef(:,1:nsoil-1), DIM=2)
    !$claw end expand
    !$claw acc end data
!#endif

    !------------------------------------------------------------------------------------------------
    ! Computation of surface diffusive heat flux from the ground and heat capacity of the ground
    !------------------------------------------------------------------------------------------------
    !$claw expand parallel update savepoint=ser_jsb_sse_calc_soil_temperature_expand_6 counter=ser_counter
!   grnd_hflx(:) = zdz1(:,1) * (t_soil_acoef(:,1) + (t_soil_bcoef(:,1) - 1._wp) * t_soil_sl(:,1))
!   hcap_grnd(:) = (zdz2(:,1) * delta_time + delta_time * (1._wp - t_soil_bcoef(:,1)) * zdz1(:,1))
    !$claw end expand

    !$claw acc end data
    !$claw verbatim ser_counter = ser_counter + 1
 !END SUBROUTINE calc_soil_temperature

#endif
END MODULE mo_top_process
