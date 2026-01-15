!>
!! @brief Contains constants for the assimilation processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-05-11)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_assimi_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_physical_constants,  ONLY: Dwv, &            !< diffusion coefficient for water vapour in air
                                        Dco2              !< diffusion coefficient for CO2 in air

  IMPLICIT NONE
  PUBLIC

  ! ======================================================================================================= !
  !> MODEL_JSBACH
  !>

  ! ======================================================================================================= !
  !> MODEL_QUINCY
  !>
  REAL(wp), SAVE :: &
    ! ratio of the molecular diffusion constants of water and CO2 in air (Nobel, Plant Phys. 2005, p 380)
    & Dwv2co2_air  = Dwv / Dco2, &                      !< ratio of diffusion coefficient for H2O and CO2 in air
    ! ratio of the diffusion constant of water and CO2 in air, accounting for turbulent transfer (Nobel, Plant Phys. 2005, p 380)
    & Dwv2co2_turb = (Dwv / Dco2) ** (2.0_wp / 3.0_wp)  !< ratio of diffusion coefficient for H2O and CO2 in turbulent air

  ! some useful numbers for easy identification of photosynthetic pathways and growth forms
  ! need to be parameters (i.e. constant across model runtime)
  ! TODO re-write as ENUM
  INTEGER, PARAMETER ::    &
    & ic3phot   = 1, &
    & ic4phot   = 2

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_assimi_constants'

#endif
END MODULE mo_assimi_constants
