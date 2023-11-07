!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief functions needed to calculate variables affecting the atm-land interface 
!!
!! @author
!!  Soenke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Soenke Zaehle (2015-12-19)
!!
MODULE mo_atmland_util
#ifndef __NO_JSBACH__

  USE mo_kind

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calc_spec_humidity_sat

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_util'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> calculation of atmospheric water content at saturation
  !!
  !! constants for the calculation of atmospheric water content at saturation
  !! from Campell and Norman's Introduction to Environmental Biophysics
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_spec_humidity_sat(t_air, press_srf) RESULT(q_sat)

    USE mo_atmland_constants,       ONLY: eps_vpd
    USE mo_jsb_physical_constants,  ONLY: Tzero

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: t_air           !< air temperature (K)
    REAL(wp), INTENT(in) :: press_srf       !< atmospheric pressure (Pa)
    REAL(wp)             :: q_sat           !< saturating specific humidity (g/g)
    ! ---------------------------
    ! 0.2 Local
    REAL(wp) :: esat            ! saturating vapour pressure (Pa)
    REAL(wp) :: a               ! Pa
    REAL(wp) :: b               ! no unit / einheitenlos
    REAL(wp) :: c               ! degC
    REAL(wp) :: tmin            ! K
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_spec_humidity_sat'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    a    = 611.0_wp
    b    = 17.502_wp
    c    = 240.97_wp
    tmin = 150._wp   ! for temperatures lower 50 K esat gets unreasonable numbers from the below calculations (very large / very small)

    ! calculating saturating vapour pressure [Pa] (Campbell & Norman)
    IF(t_air >= tmin) THEN
       esat = a * EXP(b * (t_air - Tzero) / (t_air - Tzero + c))
    ELSE
       esat = a * EXP(b * (tmin - Tzero) / (tmin - Tzero + c))
    ENDIF
    ! conversion to saturating specific humidity [kg kg-1]
    q_sat = eps_vpd * esat / (press_srf - (1.0_wp - eps_vpd) * esat)

  END FUNCTION calc_spec_humidity_sat

#endif
END MODULE mo_atmland_util
