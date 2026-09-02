!> Calculate insect damage for each timestep
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
MODULE mo_q_dinsect_damage
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message, finish, warning
  USE mo_jsb_math_constants,      ONLY: one_day, one_year, eps8, eps4
  USE mo_q_dinsect_util,            ONLY: Write2DMatrixToFile
  USE mo_util,                    ONLY: int2string

! #ifdef __QUINCY_WITH_NETCDF__
!   USE netcdf
! #endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: derive_flux_frac_with_sigmoid_from_nml


  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dinsect_damage'

CONTAINS

  ! ======================================================================================================= !
  !> Calculate the flux_frac for each time step using sigmoid function
  !> It can use the given yearly damage fraction to calcute the distribution of damage over a year,
  !> then, return the value for specific time step
  ! ----------------------------------------------------------------------------------------------------- !
  SUBROUTINE derive_flux_frac_with_sigmoid_from_nml(start_doy, &
                                                    duration, &
                                                    current_year, &
                                                    current_doy, &
                                                    current_hour, &
                                                    yearly_damage_frac, &
                                                    current_flux_frac)
    USE mo_exception,   ONLY: message

    INTEGER,          INTENT(in)  :: start_doy          !< Start day of disturbance
    INTEGER,          INTENT(in)  :: duration           !< Duration in days
    INTEGER,          INTENT(in)  :: current_year       !< Current year
    INTEGER,          INTENT(in)  :: current_doy        !< Current day of year
    REAL(wp),         INTENT(in)  :: current_hour       !< Current hour
    REAL(wp),         INTENT(in)  :: yearly_damage_frac          !< Prescribed damage fraction

    REAL(wp),         INTENT(out) :: current_flux_frac     !< Flux fraction

    !> Local variables
    REAL(wp) :: a, b, c                           !< Sigmoid parameters
    REAL(wp) :: sigma, sigma_prev                 !< Current and previous sigmoid values
    REAL(wp) :: time
    REAL(wp) :: flux_frac
    INTEGER  :: year, doy, hour_step              !< Time counters
    INTEGER  :: timestep, timesteps               !< Matrix indices
    INTEGER  :: index                             !< Index of current time step
    REAL(wp) :: hour                             !< Hour of day
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: damageMatrix
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':derive_flux_frac_with_sigmoid'

    !> 1.0 Initialize variables and allocate matrix
    !>
    timesteps = 365 * 48                         !< Total timesteps (half-hourly for one year)
    timestep = 1
    ALLOCATE(damageMatrix(4, timesteps))

    !> 2.0 Calculate sigmoid parameters
    !>
    a = yearly_damage_frac * (1 + yearly_damage_frac**2)**(1 + yearly_damage_frac**2)     !< Maximum amplitude
    b = 1.0_wp / (duration / 20.0_wp)           !< Slope of sigmoid
    c = duration / 2.0_wp                       !< Midpoint of disturbance
    sigma_prev = 0.0_wp

    !> 3.0 Calculate damage for each timestep
    !>
    day_loop: DO doy = 1, 365
      hour_loop: DO hour_step = 0, 47
        hour = REAL(hour_step) * 0.5_wp
        time = REAL(doy - start_doy) + hour / 24.0_wp  ! Time relative to the start of the disturbance

        !> 3.1 Calculate damage during disturbance period
        !>
        IF (doy >= start_doy .AND. doy <= start_doy + duration) THEN
          sigma = a / (1.0_wp + EXP(-b * (time - c)))

          !> Avoid initial spike
          IF (timestep == 1 .OR. (doy == start_doy .AND. hour_step == 0)) THEN
            flux_frac = 0.0_wp
          ELSE
            flux_frac = sigma - sigma_prev
          END IF
          sigma_prev = sigma
        ELSE
          flux_frac = 0.0_wp
        END IF

        !> 3.2 Store results in matrix
        !>
        damageMatrix(1, timestep) = current_year
        damageMatrix(2, timestep) = doy
        damageMatrix(3, timestep) = hour
        damageMatrix(4, timestep) = flux_frac

        timestep = timestep + 1
        IF (doy > 365) EXIT hour_loop
      END DO hour_loop
    END DO day_loop

  !> 4. get the flux_frac for the current time step
    index = (current_doy - 1) * 48 + current_hour * 2 + 1
    IF(damageMatrix(2,index)==current_doy) THEN
      IF(damageMatrix(3,index)==current_hour) THEN
        current_flux_frac = damageMatrix(4, index)
      ELSE
      CALL warning("derive_flux_frac_with_sigmoid_from_nml", &
        "Error: need check sigmoid function hour calculation!")
      END IF
    ELSE
    CALL warning("derive_flux_frac_with_sigmoid_from_nml", &
      "Error: need check sigmoid function doy calculation!")
    END IF

    !> 5.0 Print diagnostic damage matrix for debug
    IF (current_doy == 1 .AND. current_hour == 0) THEN
      !< I want to print the damageMatrix
      CALL Write2DMatrixToFile("damageMatrix.txt", damageMatrix)   ! Here is used for printing out the damageMatrix for debug
    END IF

  END SUBROUTINE derive_flux_frac_with_sigmoid_from_nml


#endif
END MODULE mo_q_dinsect_damage